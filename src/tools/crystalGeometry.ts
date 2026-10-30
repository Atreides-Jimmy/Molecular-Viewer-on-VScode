/**
 * Geometry shared by the periodic-structure readers.
 *
 * A CIF and a VESTA file describe the same crystal in two encodings, so the
 * arithmetic that turns a cell and a symmetry operation into atoms has to be
 * *one* implementation. When the two readers kept their own copies they drifted:
 * the CIF reader learned that 0.5 A is far too coarse a coincidence tolerance
 * (it folds genuinely distinct sites together), while the VESTA reader still
 * used the old value, and the same structure read through the two extensions
 * came back with different atom counts.
 *
 * Everything here is pure arithmetic — no file-format knowledge — so both
 * readers can share it without one depending on the other.
 */

// ===== Cell =====

export function cellToLatticeVectors(a: number, b: number, c: number, alpha: number, beta: number, gamma: number): number[][] {
    const rad = Math.PI / 180;
    const ar = alpha * rad;
    const br = beta * rad;
    const gr = gamma * rad;

    const cosA = Math.cos(ar);
    const cosB = Math.cos(br);
    const cosG = Math.cos(gr);
    const sinG = Math.sin(gr);

    const aVec = [a, 0, 0];
    const bVec = [b * cosG, b * sinG, 0];
    const cX = c * cosB;
    const cY = c * (cosA - cosB * cosG) / sinG;
    const cZ2 = 1 - cosA * cosA - cosB * cosB - cosG * cosG + 2 * cosA * cosB * cosG;
    const cZ = c * Math.sqrt(Math.max(0, cZ2)) / sinG;
    const cVec = [cX, cY, cZ];

    return [aVec, bVec, cVec];
}

export function fracToCart(frac: number[], lattice: number[][]): [number, number, number] {
    const x = frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0];
    const y = frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1];
    const z = frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2];
    return [x, y, z];
}

export function cartToFrac(cart: number[], lattice: number[][]): number[] {
    const [aVec, bVec, cVec] = lattice;
    const det = aVec[0] * (bVec[1] * cVec[2] - bVec[2] * cVec[1]) -
                aVec[1] * (bVec[0] * cVec[2] - bVec[2] * cVec[0]) +
                aVec[2] * (bVec[0] * cVec[1] - bVec[1] * cVec[0]);
    if (Math.abs(det) < 1e-10) return [0, 0, 0];
    const invDet = 1 / det;
    const inv = [
        [
            (bVec[1] * cVec[2] - bVec[2] * cVec[1]) * invDet,
            (bVec[2] * cVec[0] - bVec[0] * cVec[2]) * invDet,
            (bVec[0] * cVec[1] - bVec[1] * cVec[0]) * invDet
        ],
        [
            (aVec[2] * cVec[1] - aVec[1] * cVec[2]) * invDet,
            (aVec[0] * cVec[2] - aVec[2] * cVec[0]) * invDet,
            (aVec[1] * cVec[0] - aVec[0] * cVec[1]) * invDet
        ],
        [
            (aVec[1] * bVec[2] - aVec[2] * bVec[1]) * invDet,
            (aVec[2] * bVec[0] - aVec[0] * bVec[2]) * invDet,
            (aVec[0] * bVec[1] - aVec[1] * bVec[0]) * invDet
        ]
    ];
    return [
        inv[0][0] * cart[0] + inv[0][1] * cart[1] + inv[0][2] * cart[2],
        inv[1][0] * cart[0] + inv[1][1] * cart[1] + inv[1][2] * cart[2],
        inv[2][0] * cart[0] + inv[2][1] * cart[1] + inv[2][2] * cart[2]
    ];
}

/**
 * A fractional coordinate folded into [0, 1).
 *
 * A value a rounding error below 1 — the usual result of `-x + 1` for `x` close
 * to 0 — would otherwise name the same atom in the *next* cell, which every
 * comparison downstream would then treat as a different point. Folding it back
 * to 0 is what makes "this image is the atom at (0,0,0)" true again.
 */
export function wrapFrac(v: number): number {
    let r = v - Math.floor(v);
    if (r >= 1 - 1e-6) r -= 1;
    if (r < 0) r += 1;
    return r;
}

// ===== Symmetry =====

export interface SymOp {
    R: number[][];
    t: number[];
}

export function parseSymOp(s: string): SymOp {
    const parts = s.split(',').map(p => p.trim());
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const t = [0, 0, 0];

    for (let idx = 0; idx < 3; idx++) {
        if (idx >= parts.length) continue;
        const expr = parts[idx].replace(/\s+/g, '').toLowerCase();
        let sign = 1;
        let k = 0;
        while (k < expr.length) {
            const ch = expr[k];
            if (ch === '+') { sign = 1; k++; continue; }
            if (ch === '-') { sign = -1; k++; continue; }
            if (ch === 'x' || ch === 'y' || ch === 'z') {
                const col = ch === 'x' ? 0 : ch === 'y' ? 1 : 2;
                R[idx][col] += sign;
                sign = 1;
                k++;
                continue;
            }
            const numMatch = expr.substring(k).match(/^(\d+\/\d+|\d+\.?\d*)/);
            if (numMatch) {
                let v: number;
                if (numMatch[1].includes('/')) {
                    const [p, q] = numMatch[1].split('/').map(Number);
                    v = p / q;
                } else {
                    v = parseFloat(numMatch[1]);
                }
                t[idx] += sign * v;
                sign = 1;
                k += numMatch[1].length;
                continue;
            }
            k++;
        }
    }

    return { R, t };
}

export function applySymOp(op: SymOp, frac: number[]): number[] {
    const result = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        result[i] = op.R[i][0] * frac[0] + op.R[i][1] * frac[1] + op.R[i][2] * frac[2] + op.t[i];
    }
    return result;
}

/**
 * `1/2`, not `0.5`, where the value is a symmetry fraction.
 *
 * Crystallographic translations are almost always halves, thirds, quarters or
 * sixths, and those are what a file writes back (`-x+1/2, y, -z+1/2`). Keeping
 * the rational spelling makes an operation read out of a numeric table look like
 * the one the same structure's CIF states — which is what lets the VESTA reader
 * report a traceable operation and lets the two readers agree on equality.
 */
function formatSymFraction(v: number): string {
    if (Math.abs(v) < 1e-8) return '';
    if (Math.abs(v - Math.round(v)) < 1e-8) return String(Math.round(v));
    for (let den = 2; den <= 24; den++) {
        const num = Math.round(v * den);
        if (Math.abs(v - num / den) < 1e-6) {
            let p = Math.abs(num);
            let q = den;
            const gcd = (x: number, y: number): number => y === 0 ? x : gcd(y, x % y);
            const g = gcd(p, q);
            p /= g;
            q /= g;
            return (num < 0 ? '-' : '') + (q === 1 ? String(p) : p + '/' + q);
        }
    }
    return v.toFixed(4).replace(/0+$/, '').replace(/\.$/, '');
}

/** One row of an operation: `-x+1/2`, `y`, `2z-1/3`, `0`. */
function symOpRow(cx: number, cy: number, cz: number, t: number): string {
    let out = '';
    const vars = [cx, cy, cz];
    const names = ['x', 'y', 'z'];
    for (let i = 0; i < 3; i++) {
        const c = vars[i];
        if (Math.abs(c) < 1e-8) continue;
        const a = Math.abs(c);
        const body = (a === 1 ? '' : String(a)) + names[i];
        if (out === '') out = c < 0 ? '-' + body : body;
        else out += c < 0 ? '-' + body : '+' + body;
    }
    const ts = formatSymFraction(t);
    if (ts !== '') {
        if (out === '') out = ts;
        else out += ts.charAt(0) === '-' ? ts : '+' + ts;
    }
    return out === '' ? '0' : out;
}

/** The CIF spelling of an operation (`-x+1/2, y, -z+1/2`). */
export function symOpToString(op: SymOp): string {
    const rows: string[] = [];
    for (let i = 0; i < 3; i++) {
        rows.push(symOpRow(op.R[i][0], op.R[i][1], op.R[i][2], op.t[i]));
    }
    return rows.join(', ');
}

/**
 * A key two spellings of one operation share, for de-duplication. Operations
 * are compared on their numbers, not on their text, so `x+1/2` and `1/2+x`
 * cannot sneak into a table twice.
 */
export function symOpKey(op: SymOp): string {
    return op.R.map((row, i) => row.join(',') + '|' + op.t[i].toFixed(6)).join(';');
}

// ===== Coincident sites =====

/**
 * How close two images must be before they are the same atom, in Angstrom.
 *
 * This has to be far below a bond length: a structure can legitimately list two
 * sites 0.3 A apart (a split position, a partially occupied pair) and folding
 * them into one loses an atom that the file says is there. A tenth of an
 * angstrom only catches what symmetry expansion produces on purpose — the
 * images of a site that sits on a special position — plus genuine duplicates.
 */
export const DEDUPE_TOL = 0.1;
export const DEDUPE_TOL2 = DEDUPE_TOL * DEDUPE_TOL;

/**
 * Two images closer than this are the same point, not two points that happen to
 * be near each other: the difference is round-off in the fractional -> cartesian
 * round trip (a few 1e-9 A at most). Only merges *above* it say anything about
 * the structure rather than about the arithmetic, so only those are reported.
 */
export const EXACT_MERGE = 1e-4;

export interface NearIndex {
    add(x: number, y: number, z: number): void;
    /** Nearest atom within DEDUPE_TOL of (x,y,z), with its distance. */
    find(x: number, y: number, z: number): { index: number; dist: number } | null;
}

/**
 * A grid over the atoms placed so far, one cell per tolerance.
 *
 * Symmetry expansion compares every image against every atom already placed, so
 * a framework (thousands of sites times hundreds of operations) turns the naive
 * scan into billions of distance checks. A grid makes each lookup local; the
 * neighbourhood is 3x3x3 cells because two points closer than one cell size can
 * only ever land in adjacent cells. The result is the same as the exhaustive
 * scan — the distance test still decides — only it no longer depends on how
 * many atoms the cell already holds.
 */
export function createNearIndex(): NearIndex {
    const buckets = new Map<number, number[]>();
    const xs: number[] = [];
    const ys: number[] = [];
    const zs: number[] = [];
    const cellOf = (v: number) => Math.floor(v / DEDUPE_TOL);
    // A number hash of the three cell indices. Collisions only cost a distance
    // check, never a wrong answer.
    const bucketKey = (ix: number, iy: number, iz: number) =>
        (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791);

    return {
        add(x: number, y: number, z: number): void {
            const k = bucketKey(cellOf(x), cellOf(y), cellOf(z));
            const bucket = buckets.get(k);
            if (bucket) bucket.push(xs.length);
            else buckets.set(k, [xs.length]);
            xs.push(x);
            ys.push(y);
            zs.push(z);
        },
        find(x: number, y: number, z: number): { index: number; dist: number } | null {
            const ix = cellOf(x), iy = cellOf(y), iz = cellOf(z);
            let best = -1;
            let bestD2 = DEDUPE_TOL2;
            for (let dx = -1; dx <= 1; dx++) {
                for (let dy = -1; dy <= 1; dy++) {
                    for (let dz = -1; dz <= 1; dz++) {
                        const bucket = buckets.get(bucketKey(ix + dx, iy + dy, iz + dz));
                        if (!bucket) continue;
                        for (let n = 0; n < bucket.length; n++) {
                            const i = bucket[n];
                            const ex = xs[i] - x, ey = ys[i] - y, ez = zs[i] - z;
                            const d2 = ex * ex + ey * ey + ez * ez;
                            if (d2 < bestD2) {
                                bestD2 = d2;
                                best = i;
                            }
                        }
                    }
                }
            }
            return best < 0 ? null : { index: best, dist: Math.sqrt(bestD2) };
        }
    };
}
