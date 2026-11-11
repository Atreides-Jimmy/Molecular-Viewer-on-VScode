import { FchkBasis, CavitySphere } from '../parsers/fchkParser';

/**
 * Molecular-orbital and cavity isosurfaces, computed from the file's own data.
 *
 * Two facts shape this module:
 *
 *  - **A molecular orbital is a sum of contracted Gaussians**, and the
 *    coefficients a checkpoint stores belong to *normalized* contracted
 *    functions. Rebuilding the field therefore means expanding each shell into
 *    its Cartesian components, normalizing each primitive analytically, and
 *    normalizing the contraction once by its own self-overlap. That last step
 *    is what makes the result independent of the writer's convention: if the
 *    file's coefficients already describe normalized contractions the factor
 *    is 1 to within rounding, and if they do not, the factor corrects it.
 *    Within a shell the components are not normalized alike — a normalized
 *    `xy` primitive is sqrt(3) times a normalized `xx` one — so each component
 *    carries its own factor, and the pure (spherical) d/f functions are the
 *    real solid harmonics renormalized in that metric. The check is the
 *    orbital's own norm: the d-richest occupied orbital of the 147-atom test
 *    checkpoint integrates to 1.00002 after these factors and read 0.982
 *    without them, and the frontier orbitals read 1.0001.
 *
 *  - **The field is only needed near the surface.** A shell contributes
 *    nothing at a point farther from its centre than its own decay radius, so
 *    the sum is taken shell by shell over the grid points inside that radius
 *    instead of point by point over every shell. For a 432-atom basis that is
 *    a hundredfold saving, and it is what keeps a 4680-function orbital inside
 *    a second on one core.
 *
 * The surface itself is marching *tetrahedra* — the cube split into six
 * tetrahedra, sixteen cases, no lookup tables — because a table-driven
 * marching cubes needs a 256-entry triangle table this project would then have
 * to own and audit. Every vertex carries the field gradient as its normal, so
 * the surface shades smoothly without welding the mesh.
 *
 * Everything is computed in Bohr (the unit of the coefficient block and of
 * every length in it) and converted to Angstrom on the way out, which keeps
 * the amplitude scale of `0.02` meaning what it means in GaussView.
 */

/** Bohr to Angstrom, CODATA 2018. */
const BOHR = 0.529177210903;

/** A triangle mesh with per-vertex normals, ready for `BufferGeometry`. */
export interface SurfaceMesh {
    positions: Float32Array;
    normals: Float32Array;
    indices: Uint32Array;
}

/** What a task has done so far, for the panel's progress line. */
export interface SurfaceProgress {
    /** 0..1 over the whole task. */
    fraction: number;
    /** What it is doing now, in the panel's own words. */
    phase: string;
}

/** A running or finished surface computation. */
export interface SurfaceTask {
    /** Work for at most `budgetMs`; true when there is nothing left to do. */
    step(budgetMs: number): boolean;
    progress(): SurfaceProgress;
    /** The meshes, once `step` has returned true. */
    surfaces(): { sign: number; mesh: SurfaceMesh }[];
    /** How the grid was chosen, for the panel to state. */
    gridNote(): string;
}

/** Grid resolution presets, in Bohr. */
export const GRID_STEPS = { coarse: 0.62, normal: 0.42, fine: 0.30 } as const;
/** Point ceilings, so a huge molecule cannot ask for gigabytes. */
const GRID_CAP = { coarse: 2500000, normal: 5000000, fine: 9000000 } as const;

export type GridQuality = keyof typeof GRID_STEPS;

/** The sampled field: a regular grid in Bohr. */
export interface FieldGrid {
    ox: number;
    oy: number;
    oz: number;
    h: number;
    nx: number;
    ny: number;
    nz: number;
    values: Float32Array;
}

/** A box in Bohr. */
export interface FieldBox {
    min: [number, number, number];
    max: [number, number, number];
}

/** A grid over `box` at (about) `step`, coarsened until it fits `cap`. */
export function buildGrid(box: FieldBox, step: number, cap: number): FieldGrid {
    let h = step;
    for (let guard = 0; guard < 64; guard++) {
        const nx = Math.max(2, Math.ceil((box.max[0] - box.min[0]) / h) + 1);
        const ny = Math.max(2, Math.ceil((box.max[1] - box.min[1]) / h) + 1);
        const nz = Math.max(2, Math.ceil((box.max[2] - box.min[2]) / h) + 1);
        if (nx * ny * nz <= cap) {
            return { ox: box.min[0], oy: box.min[1], oz: box.min[2], h, nx, ny, nz, values: new Float32Array(nx * ny * nz) };
        }
        h *= 1.12;
    }
    throw new Error('the structure is too large to sample on a grid');
}

/** One primitive: exponent, with its normalization folded into `d`. */
interface EvalPrim {
    alpha: number;
    d: number;
}

/** One Cartesian component of a shell: powers and its MO coefficient. */
interface EvalComp {
    lx: number;
    ly: number;
    lz: number;
    c: number;
}

/** One contracted shell, ready to be evaluated. */
interface EvalShell {
    x: number;
    y: number;
    z: number;
    /** Distance beyond which the shell cannot change the field by more than `eps`. */
    extent: number;
    prims: EvalPrim[];
    comps: EvalComp[];
}

/** Double factorial, for the primitive normalization. */
function dfact(n: number): number {
    let r = 1;
    for (let k = n; k > 1; k -= 2) r *= k;
    return r;
}

/** Normalization constant of a Cartesian Gaussian primitive. */
function primitiveNorm(alpha: number, lx: number, ly: number, lz: number): number {
    const l = lx + ly + lz;
    const denom = dfact(2 * lx - 1) * dfact(2 * ly - 1) * dfact(2 * lz - 1);
    return Math.pow(2 * alpha / Math.PI, 0.75) * Math.sqrt(Math.pow(4 * alpha, l) / denom);
}

/** `int t^n exp(-p t^2) dt` over the whole line. */
function axisMoment(n: number, p: number): number {
    if (n % 2 === 1) return 0;
    return dfact(n - 1) / Math.pow(2 * p, n / 2) * Math.sqrt(Math.PI / p);
}

/**
 * Self-overlap of a contraction.
 *
 * For one axis the moment of a Gaussian about its own centre is
 * `(n-1)!! / (2p)^(n/2) * sqrt(pi/p)` for even `n` and zero for odd, so the
 * whole integral is a product of three such moments times the primitive norms
 * and the contraction coefficients. Only same-centre pairs appear (this is the
 * contraction with itself), so there is no centre-separation factor to carry.
 */
export function contractionSelfOverlap(prims: { alpha: number; d: number }[], lx: number, ly: number, lz: number): number {
    let total = 0;
    for (let i = 0; i < prims.length; i++) {
        const pi = prims[i];
        const ni = primitiveNorm(pi.alpha, lx, ly, lz);
        for (let j = 0; j < prims.length; j++) {
            const pj = prims[j];
            const nj = primitiveNorm(pj.alpha, lx, ly, lz);
            const p = pi.alpha + pj.alpha;
            // The three one-dimensional moments already carry the (pi/p)^(3/2)
            // of the Gaussian integral: multiplying by it again would shrink
            // every contraction by an arbitrary factor.
            total += pi.d * pj.d * ni * nj *
                axisMoment(2 * lx, p) * axisMoment(2 * ly, p) * axisMoment(2 * lz, p);
        }
    }
    return total;
}

/** Angular powers of one Cartesian component, in the order a shell lists them. */
function cartesianPowers(l: number): number[][] {
    switch (l) {
        case 0: return [[0, 0, 0]];
        case 1: return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        // Gaussian's 6D order: XX YY ZZ XY XZ YZ.
        case 2: return [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]];
        // Gaussian's 10F order: XXX YYY ZZZ XYY XXY XXZ XZZ YZZ YYZ XYZ.
        case 3: return [
            [3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 2, 0], [2, 1, 0],
            [2, 0, 1], [1, 0, 2], [0, 1, 2], [0, 2, 1], [1, 1, 1]
        ];
        default: return [];
    }
}

/**
 * Coefficients of the five real d functions in the six Cartesian d functions.
 *
 * The pure functions are the real spherical harmonics times r^2, expressed in
 * the *normalized* Cartesian functions: `d(z2) = (2 zz - xx - yy)/sqrt(6)`,
 * `d(x2-y2) = (xx - yy)/sqrt(2)`, and `xz`, `yz`, `xy` themselves. Order
 * matches Gaussian: m = 0, +1, -1, +2, -2.
 */
const PURE_D_TO_CARTESIAN: number[][] = [
    // xx,                 yy,                 zz,                 xy, xz, yz
    [-1 / Math.sqrt(6), -1 / Math.sqrt(6), 2 / Math.sqrt(6), 0, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
    [1 / Math.sqrt(2), -1 / Math.sqrt(2), 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0]
];

/**
 * The same for the seven real f functions in the ten Cartesian ones, in
 * Gaussian's order (m = 0, +1, -1, +2, -2, +3, -3 over XXX YYY ZZZ XYY XXY XXZ
 * XZZ YZZ YYZ XYZ). Each row is one spherical function's mixture.
 */
const PURE_F_TO_CARTESIAN: number[][] = [
    // xxx,               yyy,               zzz,               xyy,               xxy,               xxz,               xzz,               yzz,                yyz,               xyz
    [0, 0, 2 / Math.sqrt(22), 0, 0, -3 / Math.sqrt(22), 0, 0, -3 / Math.sqrt(22), 0],
    [-1 / Math.sqrt(18), 0, 0, -1 / Math.sqrt(18), 0, 0, 4 / Math.sqrt(18), 0, 0, 0],
    [0, -1 / Math.sqrt(18), 0, 0, -1 / Math.sqrt(18), 0, 0, 4 / Math.sqrt(18), 0, 0],
    [0, 0, 0, 0, 0, 1 / Math.sqrt(2), 0, 0, -1 / Math.sqrt(2), 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [1 / Math.sqrt(10), 0, 0, -3 / Math.sqrt(10), 0, 0, 0, 0, 0, 0],
    [0, -1 / Math.sqrt(10), 0, 0, 3 / Math.sqrt(10), 0, 0, 0, 0, 0]
];

/**
 * Overlap of two normalized Cartesian components of one shell.
 *
 * Both are normalized primitives at the same centre with the same exponent, so
 * the integral factors over the three axes and the exponent cancels:
 * `(l1+l2-1)!! / sqrt((2l1-1)!!(2l2-1)!!)` per axis, and zero when an axis's
 * total power is odd. The result is the metric that the rows above are measured
 * in: the rows are the real solid harmonics, which are orthonormal *as
 * functions*, not with respect to the plain coefficient vector.
 */
function cartesianOverlap(powers: number[][]): number[][] {
    const n = powers.length;
    const m: number[][] = [];
    for (let i = 0; i < n; i++) {
        m.push(new Array<number>(n).fill(0));
        for (let j = 0; j < n; j++) {
            let v = 1;
            for (let a = 0; a < 3; a++) {
                const l1 = powers[i][a];
                const l2 = powers[j][a];
                if ((l1 + l2) % 2 === 1) { v = 0; break; }
                v *= dfact(l1 + l2 - 1) / Math.sqrt(dfact(2 * l1 - 1) * dfact(2 * l2 - 1));
            }
            m[i][j] = v;
        }
    }
    return m;
}

/**
 * The pure-shell mixtures, renormalized in the Cartesian metric.
 *
 * A row such as `(2 zz - xx - yy) / sqrt(6)` is the real solid harmonic, but
 * the six Cartesian components are not orthogonal to each other (xx·yy is 1/3
 * of xx·xx), so that row's norm is sqrt(2/3), not 1. Dividing a row by its norm
 * in this metric is what makes each pure d/f basis function a normalized
 * function, which is what the file's coefficients multiply.
 */
const PURE_MIX_CACHE: { [l: number]: number[][] } = {};
function pureMixture(l: number): number[][] {
    const cached = PURE_MIX_CACHE[l];
    if (cached) return cached;
    const raw = l === 2 ? PURE_D_TO_CARTESIAN : PURE_F_TO_CARTESIAN;
    const overlap = cartesianOverlap(cartesianPowers(l));
    const out: number[][] = [];
    const n = overlap.length;
    for (const row of raw) {
        let q = 0;
        for (let i = 0; i < n; i++) {
            if (row[i] === 0) continue;
            for (let j = 0; j < n; j++) q += row[i] * row[j] * overlap[i][j];
        }
        const s = Math.sqrt(q > 0 ? q : 1);
        out.push(row.slice(0, n).map(v => v / s));
    }
    PURE_MIX_CACHE[l] = out;
    return out;
}

/** Cartesian functions one shell type contributes. */
export function shellFunctionCount(type: number): number {
    switch (type) {
        case 0: return 1;
        case 1: return 3;
        case -1: return 4;
        case 2: return 6;
        case -2: return 5;
        case 3: return 10;
        case -3: return 7;
        case 4: return 15;
        case -4: return 9;
        default: return 0;
    }
}

/** The angular momentum a shell type stands for. */
function shellAngular(type: number): number {
    return type === 0 ? 0 : Math.min(4, Math.abs(type));
}

/** Largest angular momentum this engine evaluates. */
export const MAX_SUPPORTED_L = 3;

/**
 * Expand one Cartesian shell into its components with their MO coefficients.
 *
 * Returns null for an angular momentum this engine cannot represent (g and
 * beyond), so the caller refuses the picture instead of drawing a wrong one.
 */
function expandCartesianShell(l: number, pure: boolean, coefficients: ArrayLike<number>, offset: number): EvalComp[] | null {
    if (l > MAX_SUPPORTED_L) return null;
    const powers = cartesianPowers(l);
    if (!pure || l < 2) {
        const comps: EvalComp[] = [];
        for (let k = 0; k < powers.length; k++) {
            comps.push({ lx: powers[k][0], ly: powers[k][1], lz: powers[k][2], c: coefficients[offset + k] || 0 });
        }
        return comps;
    }
    const mix = pureMixture(l);
    const cart = new Array<number>(powers.length).fill(0);
    for (let m = 0; m < mix.length; m++) {
        const c = coefficients[offset + m] || 0;
        if (c === 0) continue;
        const row = mix[m];
        for (let k = 0; k < cart.length; k++) cart[k] += c * row[k];
    }
    const comps: EvalComp[] = [];
    for (let k = 0; k < powers.length; k++) {
        comps.push({ lx: powers[k][0], ly: powers[k][1], lz: powers[k][2], c: cart[k] });
    }
    return comps;
}

/**
 * Radius beyond which a primitive cannot change the field by more than `eps`.
 *
 * `N d exp(-a r^2) r^l < eps` is solved for r by a few fixed-point steps; the
 * polynomial factor `r^l` is what makes a diffuse d function's cutoff a little
 * larger than a tight s function's, which matters at the 1e-5 level this is
 * asked for.
 */
function primitiveExtent(alpha: number, l: number, scale: number, eps: number): number {
    const e = Math.log(Math.max(scale, 1e-300) / eps);
    if (e <= 0) return 0;
    let r = Math.sqrt(e / alpha);
    for (let k = 0; k < 4; k++) {
        const v = (e + l * Math.log(Math.max(r, 1e-12))) / alpha;
        if (!(v > 0)) break;
        r = Math.sqrt(v);
    }
    return r;
}

/**
 * Turn a basis plus one orbital's coefficients into evaluable shells.
 *
 * `eps` is the absolute amplitude a shell is allowed to be dropped at; 1e-5 is
 * a twentieth of a percent of the default isovalue, well below what any contour
 * can show.
 */
export function buildEvalShells(basis: FchkBasis, coefficients: ArrayLike<number>, eps = 1e-5): EvalShell[] {
    const shells: EvalShell[] = [];
    let offset = 0;
    let prim = 0;
    for (let s = 0; s < basis.shellType.length; s++) {
        const type = basis.shellType[s];
        const nprim = basis.shellPrim[s] || 0;
        const nfun = shellFunctionCount(type);
        if (nfun === 0) throw new Error('this checkpoint uses a shell type this viewer does not know (' + type + ')');
        const x = basis.centers[s * 3] || 0;
        const y = basis.centers[s * 3 + 1] || 0;
        const z = basis.centers[s * 3 + 2] || 0;
        const l = shellAngular(type);

        const prims: { alpha: number; d: number }[] = [];
        for (let k = 0; k < nprim; k++) {
            const alpha = basis.exponents[prim + k];
            if (!(alpha > 0)) continue;
            prims.push({ alpha, d: basis.coefficients[prim + k] || 0 });
        }
        if (type === -1) {
            // An sp shell is one s function and three p functions sharing the
            // exponents but not the contraction coefficients.
            const sComps: EvalComp[] = [{ lx: 0, ly: 0, lz: 0, c: coefficients[offset] || 0 }];
            const pComps: EvalComp[] = [];
            const powers = cartesianPowers(1);
            for (let k = 0; k < 3; k++) {
                pComps.push({ lx: powers[k][0], ly: powers[k][1], lz: powers[k][2], c: coefficients[offset + 1 + k] || 0 });
            }
            const pPrims: { alpha: number; d: number }[] = [];
            for (let k = 0; k < nprim; k++) {
                const alpha = basis.exponents[prim + k];
                if (!(alpha > 0)) continue;
                pPrims.push({ alpha, d: basis.spCoefficients[prim + k] || 0 });
            }
            shells.push(makeShell(x, y, z, 0, prims, sComps, eps));
            shells.push(makeShell(x, y, z, 1, pPrims, pComps, eps));
        } else {
            // The sign of the shell type is the file's own statement of the
            // shell's angular form (-2 = five pure d, -3 = seven pure f), and it
            // is also what fixes how many coefficients the shell consumes. The
            // job-level `Pure/Cartesian` scalars are not used here: a file whose
            // per-shell sign and global scalar disagree would otherwise read the
            // wrong number of coefficients and shift every shell after it.
            const pure = type < 0 && l >= 2;
            const comps = expandCartesianShell(l, pure, coefficients, offset);
            if (!comps) {
                throw new Error('this basis uses g functions (angular momentum ' + l + '), which this viewer does not evaluate');
            }
            shells.push(makeShell(x, y, z, l, prims, comps, eps));
        }
        offset += nfun;
        prim += nprim;
    }
    return shells;
}

/** Assemble one evaluable shell, normalized end to end. */
function makeShell(x: number, y: number, z: number, l: number, prims: { alpha: number; d: number }[], comps: EvalComp[], eps: number): EvalShell {
    // Every component of a shell shares one radial part and, because the ratio
    // of their primitive normalizations does not depend on the exponent, one
    // contraction normalization. What differs is the (2l-1)!! factor of each
    // Cartesian component: a normalized xy primitive is sqrt(3) times a
    // normalized xx one, and a Gaussian's coefficients multiply normalized
    // functions. The ratio is folded into the component's coefficient once
    // here; without it every xy/xz/yz lobe of a d shell is drawn 1/sqrt(3) too
    // small and a pure d mixture is distorted (measured: a d-rich orbital then
    // integrates to 0.982 instead of 1).
    const ref = comps.length > 0 ? comps[0] : { lx: l, ly: 0, lz: 0, c: 1 };
    const self = contractionSelfOverlap(prims, ref.lx, ref.ly, ref.lz);
    const norm = self > 0 ? 1 / Math.sqrt(self) : 1;
    const refDenom = dfact(2 * ref.lx - 1) * dfact(2 * ref.ly - 1) * dfact(2 * ref.lz - 1);
    const scaled: EvalComp[] = [];
    for (const c of comps) {
        const denom = dfact(2 * c.lx - 1) * dfact(2 * c.ly - 1) * dfact(2 * c.lz - 1);
        scaled.push({ lx: c.lx, ly: c.ly, lz: c.lz, c: c.c * Math.sqrt(refDenom / denom) });
    }

    const out: EvalPrim[] = [];
    let extent = 0;
    for (const p of prims) {
        if (p.d === 0) continue;
        const d = p.d * primitiveNorm(p.alpha, ref.lx, ref.ly, ref.lz) * norm;
        out.push({ alpha: p.alpha, d });
        let scale = Math.abs(d);
        for (const c of scaled) {
            const s2 = Math.abs(d * c.c);
            if (s2 > scale) scale = s2;
        }
        const r = primitiveExtent(p.alpha, l, scale, eps);
        if (r > extent) extent = r;
    }
    return { x, y, z, extent, prims: out, comps: scaled };
}

/**
 * Add the field of `shells[from..to)` to a grid.
 *
 * The innermost loop is the whole feature's cost centre: component powers and
 * coefficients are hoisted into typed arrays, and the exponent loop is the
 * only thing between a grid point and its amplitude.
 */
export function sampleOrbitalRange(shells: EvalShell[], from: number, to: number, grid: FieldGrid): void {
    const { ox, oy, oz, h, nx, ny, nz, values } = grid;
    for (let s = from; s < to; s++) {
        const shell = shells[s];
        if (shell.extent <= 0 || shell.prims.length === 0) continue;
        const i0 = Math.max(0, Math.ceil((shell.x - shell.extent - ox) / h));
        const i1 = Math.min(nx - 1, Math.floor((shell.x + shell.extent - ox) / h));
        const j0 = Math.max(0, Math.ceil((shell.y - shell.extent - oy) / h));
        const j1 = Math.min(ny - 1, Math.floor((shell.y + shell.extent - oy) / h));
        const k0 = Math.max(0, Math.ceil((shell.z - shell.extent - oz) / h));
        const k1 = Math.min(nz - 1, Math.floor((shell.z + shell.extent - oz) / h));
        if (i0 > i1 || j0 > j1 || k0 > k1) continue;
        const ext2 = shell.extent * shell.extent;
        const comps = shell.comps;
        const nPrim = shell.prims.length;
        const nComp = comps.length;
        const px = new Int32Array(nComp);
        const py = new Int32Array(nComp);
        const pz = new Int32Array(nComp);
        const pc = new Float64Array(nComp);
        for (let q = 0; q < nComp; q++) {
            px[q] = comps[q].lx; py[q] = comps[q].ly; pz[q] = comps[q].lz; pc[q] = comps[q].c;
        }
        const alpha = new Float64Array(nPrim);
        const dcoef = new Float64Array(nPrim);
        for (let q = 0; q < nPrim; q++) { alpha[q] = shell.prims[q].alpha; dcoef[q] = shell.prims[q].d; }

        for (let k = k0; k <= k1; k++) {
            const dz = oz + k * h - shell.z;
            const dz2 = dz * dz;
            for (let j = j0; j <= j1; j++) {
                const dy = oy + j * h - shell.y;
                const dy2 = dz2 + dy * dy;
                if (dy2 > ext2) continue;
                let idx = i0 + nx * (j + ny * k);
                for (let i = i0; i <= i1; i++, idx++) {
                    const dx = ox + i * h - shell.x;
                    const r2 = dy2 + dx * dx;
                    if (r2 > ext2) continue;
                    let radial = 0;
                    for (let q = 0; q < nPrim; q++) radial += dcoef[q] * Math.exp(-alpha[q] * r2);
                    if (radial === 0) continue;
                    let angular = 0;
                    for (let q = 0; q < nComp; q++) {
                        let term = pc[q];
                        if (term === 0) continue;
                        if (px[q] === 1) term *= dx;
                        else if (px[q] === 2) term *= dx * dx;
                        else if (px[q] === 3) term *= dx * dx * dx;
                        if (py[q] === 1) term *= dy;
                        else if (py[q] === 2) term *= dy * dy;
                        else if (py[q] === 3) term *= dy * dy * dy;
                        if (pz[q] === 1) term *= dz;
                        else if (pz[q] === 2) term *= dz * dz;
                        else if (pz[q] === 3) term *= dz * dz * dz;
                        angular += term;
                    }
                    if (angular !== 0) values[idx] += radial * angular;
                }
            }
        }
    }
}

/* ===== Marching tetrahedra ===== */

/** Growable mesh buffers. */
class MeshBuilder {
    private pos = new Float32Array(3 * 2048);
    private nrm = new Float32Array(3 * 2048);
    private idx = new Uint32Array(2048);
    private vcount = 0;
    private icount = 0;

    /** Push a vertex at grid position (gx, gy, gz), normal from the field. */
    vertexAt(grid: FieldGrid, gx: number, gy: number, gz: number): number {
        if ((this.vcount + 1) * 3 > this.pos.length) {
            const p = new Float32Array(this.pos.length * 2);
            p.set(this.pos);
            this.pos = p;
            const n = new Float32Array(this.nrm.length * 2);
            n.set(this.nrm);
            this.nrm = n;
        }
        const nx = (sampleGrid(grid, gx + 1, gy, gz) - sampleGrid(grid, gx - 1, gy, gz)) * 0.5;
        const ny = (sampleGrid(grid, gx, gy + 1, gz) - sampleGrid(grid, gx, gy - 1, gz)) * 0.5;
        const nz = (sampleGrid(grid, gx, gy, gz + 1) - sampleGrid(grid, gx, gy, gz - 1)) * 0.5;
        const len = Math.sqrt(nx * nx + ny * ny + nz * nz);
        const at = this.vcount * 3;
        this.pos[at] = (gx * grid.h + grid.ox) * BOHR;
        this.pos[at + 1] = (gy * grid.h + grid.oy) * BOHR;
        this.pos[at + 2] = (gz * grid.h + grid.oz) * BOHR;
        if (len > 1e-12) {
            this.nrm[at] = nx / len; this.nrm[at + 1] = ny / len; this.nrm[at + 2] = nz / len;
        } else {
            this.nrm[at] = 0; this.nrm[at + 1] = 0; this.nrm[at + 2] = 0;
        }
        return this.vcount++;
    }

    triangle(a: number, b: number, c: number): void {
        if (this.icount + 3 > this.idx.length) {
            const i = new Uint32Array(this.idx.length * 2);
            i.set(this.idx);
            this.idx = i;
        }
        this.idx[this.icount++] = a;
        this.idx[this.icount++] = b;
        this.idx[this.icount++] = c;
    }

    finish(): SurfaceMesh {
        return {
            positions: this.pos.subarray(0, this.vcount * 3),
            normals: this.nrm.subarray(0, this.vcount * 3),
            indices: this.idx.subarray(0, this.icount)
        };
    }
}

/** Trilinear sample of the grid at fractional grid coordinates; clamped at the border. */
function sampleGrid(grid: FieldGrid, x: number, y: number, z: number): number {
    const { nx, ny, nz, values } = grid;
    let i = Math.floor(x), j = Math.floor(y), k = Math.floor(z);
    if (i < 0) i = 0; else if (i > nx - 2) i = nx - 2;
    if (j < 0) j = 0; else if (j > ny - 2) j = ny - 2;
    if (k < 0) k = 0; else if (k > nz - 2) k = nz - 2;
    const fx = Math.min(1, Math.max(0, x - i));
    const fy = Math.min(1, Math.max(0, y - j));
    const fz = Math.min(1, Math.max(0, z - k));
    const v000 = values[i + nx * (j + ny * k)];
    const v100 = values[i + 1 + nx * (j + ny * k)];
    const v010 = values[i + nx * (j + 1 + ny * k)];
    const v110 = values[i + 1 + nx * (j + 1 + ny * k)];
    const v001 = values[i + nx * (j + ny * (k + 1))];
    const v101 = values[i + 1 + nx * (j + ny * (k + 1))];
    const v011 = values[i + nx * (j + 1 + ny * (k + 1))];
    const v111 = values[i + 1 + nx * (j + 1 + ny * (k + 1))];
    const a = v000 + (v100 - v000) * fx;
    const b = v010 + (v110 - v010) * fx;
    const c = v001 + (v101 - v001) * fx;
    const d = v011 + (v111 - v011) * fx;
    const e = a + (b - a) * fy;
    const f = c + (d - c) * fy;
    return e + (f - e) * fz;
}

/** The six tetrahedra a cube splits into, all sharing the 0-7 diagonal. */
const TETS: number[][] = [
    [0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7],
    [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]
];

/** Corner offsets of the cube, in the order TETS indexes them. */
const CORNERS: number[][] = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]
];

/** The six edges of a tetrahedron, as corner pairs. */
const TET_EDGES: number[][] = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];

/**
 * Extract one isosurface from the grid.
 *
 * `sign` decides which side of zero is contoured (`iso * sign`), because a
 * molecular orbital has two lobes of opposite phase and they are drawn as two
 * separate surfaces. The winding of a triangle is not fixed: every vertex
 * carries an explicit normal and the material draws both faces, so the shading
 * is right either way and this code has one class of bug fewer.
 */
export function extractSurface(grid: FieldGrid, iso: number, sign: number): SurfaceMesh {
    const level = iso * sign;
    const builder = new MeshBuilder();
    const { ox, oy, oz, h, nx, ny, nz, values } = grid;
    const v = new Float64Array(8);
    const p = new Float64Array(12);
    const cross = new Float64Array(12);

    for (let k = 0; k < nz - 1; k++) {
        for (let j = 0; j < ny - 1; j++) {
            let base = nx * (j + ny * k);
            for (let i = 0; i < nx - 1; i++, base++) {
                // The eight corners, and a quick rejection: the overwhelming
                // majority of cells lie wholly on one side of the level.
                v[0] = values[base];
                v[1] = values[base + 1];
                v[2] = values[base + nx];
                v[3] = values[base + nx + 1];
                v[4] = values[base + nx * ny];
                v[5] = values[base + nx * ny + 1];
                v[6] = values[base + nx * ny + nx];
                v[7] = values[base + nx * ny + nx + 1];
                let below = 0;
                for (let c = 0; c < 8; c++) if (v[c] < level) below++;
                if (below === 0 || below === 8) continue;

                for (let t = 0; t < 6; t++) {
                    const tet = TETS[t];
                    let mask = 0;
                    for (let c = 0; c < 4; c++) {
                        const idx = tet[c];
                        p[c * 3] = ox + (i + CORNERS[idx][0]) * h;
                        p[c * 3 + 1] = oy + (j + CORNERS[idx][1]) * h;
                        p[c * 3 + 2] = oz + (k + CORNERS[idx][2]) * h;
                        if (v[idx] < level) mask |= 1 << c;
                    }
                    if (mask === 0 || mask === 15) continue;

                    // Crossing points on the six edges of this tetrahedron.
                    let count = 0;
                    for (let e = 0; e < 6; e++) {
                        const a = TET_EDGES[e][0];
                        const b = TET_EDGES[e][1];
                        if (((mask >> a) & 1) === ((mask >> b) & 1)) continue;
                        const va = v[tet[a]], vb = v[tet[b]];
                        const denom = vb - va;
                        const s = Math.abs(denom) > 1e-30 ? (level - va) / denom : 0.5;
                        const ax = p[a * 3], ay = p[a * 3 + 1], az = p[a * 3 + 2];
                        cross[count * 3] = ax + (p[b * 3] - ax) * s;
                        cross[count * 3 + 1] = ay + (p[b * 3 + 1] - ay) * s;
                        cross[count * 3 + 2] = az + (p[b * 3 + 2] - az) * s;
                        count++;
                    }
                    if (count === 3) {
                        builder.triangle(
                            gridVertex(builder, grid, cross[0], cross[1], cross[2]),
                            gridVertex(builder, grid, cross[3], cross[4], cross[5]),
                            gridVertex(builder, grid, cross[6], cross[7], cross[8]));
                    } else if (count === 4) {
                        // A quad: split along its shorter diagonal, so the two
                        // triangles stay as close to the true surface as the
                        // four sampled points allow.
                        const d02 = Math.hypot(cross[0] - cross[6], cross[1] - cross[7], cross[2] - cross[8]);
                        const d13 = Math.hypot(cross[3] - cross[9], cross[4] - cross[10], cross[5] - cross[11]);
                        if (d02 <= d13) {
                            builder.triangle(
                                gridVertex(builder, grid, cross[0], cross[1], cross[2]),
                                gridVertex(builder, grid, cross[3], cross[4], cross[5]),
                                gridVertex(builder, grid, cross[6], cross[7], cross[8]));
                            builder.triangle(
                                gridVertex(builder, grid, cross[0], cross[1], cross[2]),
                                gridVertex(builder, grid, cross[6], cross[7], cross[8]),
                                gridVertex(builder, grid, cross[9], cross[10], cross[11]));
                        } else {
                            builder.triangle(
                                gridVertex(builder, grid, cross[0], cross[1], cross[2]),
                                gridVertex(builder, grid, cross[3], cross[4], cross[5]),
                                gridVertex(builder, grid, cross[9], cross[10], cross[11]));
                            builder.triangle(
                                gridVertex(builder, grid, cross[3], cross[4], cross[5]),
                                gridVertex(builder, grid, cross[6], cross[7], cross[8]),
                                gridVertex(builder, grid, cross[9], cross[10], cross[11]));
                        }
                    }
                }
            }
        }
    }
    return builder.finish();
}

/** Push one vertex given its position in Bohr-space grid coordinates. */
function gridVertex(builder: MeshBuilder, grid: FieldGrid, x: number, y: number, z: number): number {
    return builder.vertexAt(grid, (x - grid.ox) / grid.h, (y - grid.oy) / grid.h, (z - grid.oz) / grid.h);
}

/* ===== The two tasks ===== */

/**
 * Shared phase machine: sample the field shell batch by shell batch, then cut
 * one surface per requested level. Every `step` is bounded by wall time, so the
 * host's message loop keeps running and the panel's progress bar moves.
 */
abstract class GridTask implements SurfaceTask {
    protected grid: FieldGrid;
    protected meshes: { sign: number; mesh: SurfaceMesh }[] = [];
    protected done = false;
    protected sampled = 0;
    protected level = 0;
    private note = '';

    protected constructor(box: FieldBox, quality: GridQuality) {
        this.grid = buildGrid(box, GRID_STEPS[quality], GRID_CAP[quality]);
        this.note = (this.grid.h * BOHR).toFixed(2) + ' A grid, ' +
            (this.grid.nx * this.grid.ny * this.grid.nz / 1e6).toFixed(2) + ' M points';
    }

    /** Total units of sampling work, and one unit's worth of it. */
    protected abstract sampleUnits(): number;
    protected abstract sampleUnit(unit: number): void;
    protected abstract isolevels(): number[];

    step(budgetMs: number): boolean {
        if (this.done) return true;
        const end = Date.now() + budgetMs;
        const units = Math.max(1, this.sampleUnits());
        while (Date.now() < end) {
            if (this.sampled < units) {
                // A batch, not a single unit: the clock read costs more than a
                // tiny shell does, so the loop amortises it over a few.
                const batch = Math.min(units, this.sampled + 8);
                for (let u = this.sampled; u < batch; u++) this.sampleUnit(u);
                this.sampled = batch;
                continue;
            }
            const levels = this.isolevels();
            if (this.level < levels.length) {
                const at = levels[this.level];
                const mesh = extractSurface(this.grid, Math.abs(at), at < 0 ? -1 : 1);
                if (mesh.indices.length > 0) this.meshes.push({ sign: at < 0 ? -1 : 1, mesh });
                this.level++;
                continue;
            }
            this.done = true;
            return true;
        }
        return false;
    }

    progress(): SurfaceProgress {
        if (this.done) return { fraction: 1, phase: 'Done' };
        const levels = this.isolevels().length;
        const units = Math.max(1, this.sampleUnits());
        const total = 1 + levels;
        const sampling = Math.min(1, this.sampled / units);
        return {
            fraction: Math.min(0.999, (sampling + this.level) / total),
            phase: sampling < 1 ? 'Sampling the field' : 'Cutting the surface'
        };
    }

    surfaces(): { sign: number; mesh: SurfaceMesh }[] { return this.meshes; }
    gridNote(): string { return this.note; }
}

/** One molecular orbital, both lobes. */
export class OrbitalSurfaceTask extends GridTask {
    private shells: EvalShell[];

    constructor(basis: FchkBasis, coefficients: ArrayLike<number>, atoms: { x: number; y: number; z: number }[], private readonly iso: number, quality: GridQuality) {
        super(boxOf(atoms, 6), quality);
        this.shells = buildEvalShells(basis, coefficients);
        if (this.shells.length === 0) throw new Error('this orbital has no basis functions to evaluate');
    }

    protected sampleUnits(): number { return this.shells.length; }

    protected sampleUnit(unit: number): void {
        sampleOrbitalRange(this.shells, unit, unit + 1, this.grid);
    }

    protected isolevels(): number[] { return [this.iso, -this.iso]; }
}

/** The solvation cavity: the union of the spheres the file lists. */
export class CavitySurfaceTask extends GridTask {
    private spheres: CavitySphere[];

    constructor(spheres: CavitySphere[], quality: GridQuality) {
        super(boxOfSpheres(spheres), quality);
        this.spheres = spheres;
    }

    protected sampleUnits(): number { return this.spheres.length; }

    /**
     * Paint one sphere into the field as `distance - radius` minima.
     *
     * Every point starts far outside, so the zero level of the resulting
     * distance field is the outer envelope of the union — which is the cavity
     * GePol builds, before its smoothing.
     */
    protected sampleUnit(unit: number): void {
        const grid = this.grid;
        const { ox, oy, oz, h, nx, ny, nz, values } = grid;
        if (this.sampled === 0 && unit === 0) values.fill(1e9);
        const sphere = this.spheres[unit];
        const x = sphere.x / BOHR, y = sphere.y / BOHR, z = sphere.z / BOHR, r = sphere.r / BOHR;
        const reach = r + 2 * h;
        const i0 = Math.max(0, Math.floor((x - reach - ox) / h));
        const i1 = Math.min(nx - 1, Math.ceil((x + reach - ox) / h));
        const j0 = Math.max(0, Math.floor((y - reach - oy) / h));
        const j1 = Math.min(ny - 1, Math.ceil((y + reach - oy) / h));
        const k0 = Math.max(0, Math.floor((z - reach - oz) / h));
        const k1 = Math.min(nz - 1, Math.ceil((z + reach - oz) / h));
        for (let k = k0; k <= k1; k++) {
            const dz = oz + k * h - z;
            for (let j = j0; j <= j1; j++) {
                const dy = oy + j * h - y;
                let idx = i0 + nx * (j + ny * k);
                for (let i = i0; i <= i1; i++, idx++) {
                    const dx = ox + i * h - x;
                    const d = Math.sqrt(dx * dx + dy * dy + dz * dz) - r;
                    if (d < values[idx]) values[idx] = d;
                }
            }
        }
    }

    protected isolevels(): number[] { return [0]; }
}

/** The box a set of atoms occupies, plus `pad` Angstrom. */
function boxOf(atoms: { x: number; y: number; z: number }[], pad: number): FieldBox {
    let minX = Infinity, minY = Infinity, minZ = Infinity;
    let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
    for (const a of atoms) {
        const x = a.x / BOHR, y = a.y / BOHR, z = a.z / BOHR;
        if (x < minX) minX = x; if (x > maxX) maxX = x;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
        if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
    }
    const p = pad / BOHR;
    return { min: [minX - p, minY - p, minZ - p], max: [maxX + p, maxY + p, maxZ + p] };
}

/** The box a set of spheres occupies, plus two grid steps of room. */
function boxOfSpheres(spheres: CavitySphere[]): FieldBox {
    let minX = Infinity, minY = Infinity, minZ = Infinity;
    let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
    for (const s of spheres) {
        const x = s.x / BOHR, y = s.y / BOHR, z = s.z / BOHR, r = s.r / BOHR;
        minX = Math.min(minX, x - r); maxX = Math.max(maxX, x + r);
        minY = Math.min(minY, y - r); maxY = Math.max(maxY, y + r);
        minZ = Math.min(minZ, z - r); maxZ = Math.max(maxZ, z + r);
    }
    const pad = 1.5;
    return { min: [minX - pad, minY - pad, minZ - pad], max: [maxX + pad, maxY + pad, maxZ + pad] };
}
