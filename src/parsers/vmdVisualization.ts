import { AtomGroup, VmdMolecule, VmdRep, VmdScene } from '../types';

/**
 * Everything about the visualization files this toolchain produces: the VMD
 * TCL scripts that describe how a structure should be drawn, and the bare
 * coordinate files whose temperature-factor column carries the per-atom scalar
 * a visualization is made of.
 *
 * Three layers live here, in this order:
 *   1. the reference renderer's palettes and radii (colour ids, element
 *      colours, van der Waals radii, the RGB scalar gradient) — measured
 *      constants, nothing here needs the reference renderer installed;
 *   2. the resolution of a parsed script into drawing properties — which
 *      atoms and bonds a representation covers, in which colour, at which
 *      size, and which surface patches it asks for;
 *   3. the parsing of the script itself, and the recognition of a bare
 *      coordinate file that carries a scalar channel.
 */

/**
 * VMD's built-in colour table as used by `mol color ColorID <i>` (the index is
 * taken modulo 33). Entries 0-13, 16 and 17 are the measured RGB values of the
 * reference palette; the remaining indices keep the previously used colours,
 * which are close approximations rather than measurements.
 */
export const VMD_COLOR_IDS: string[] = [
    '#0000FF', // 0  blue
    '#FF0000', // 1  red
    '#595959', // 2  gray
    '#FF8000', // 3  orange
    '#FFFF00', // 4  yellow
    '#808033', // 5  tan
    '#999999', // 6  silver
    '#00FF00', // 7  green
    '#FFFFFF', // 8  white
    '#FF9999', // 9  pink
    '#40BFBF', // 10 cyan
    '#A600A6', // 11 purple
    '#80E666', // 12 lime
    '#E666B3', // 13 mauve
    '#CC7722', // 14
    '#5F9EA0', // 15
    '#000000', // 16 black
    '#E0F705', // 17 yellow2
    '#999900', // 18
    '#00CC00', // 19
    '#009900', // 20
    '#00CCCC', // 21
    '#009999', // 22
    '#0000CC', // 23
    '#000099', // 24
    '#CC00CC', // 25
    '#990099', // 26
    '#FF00FF', // 27
    '#CC0099', // 28
    '#CC0000', // 29
    '#990000', // 30
    '#CC6600', // 31
    '#996600'  // 32
];

/** Colour used for an index or element the palette does not know. */
export const VMD_FALLBACK_COLOR = '#FF1493';

/**
 * `mol color Name` / `mol color Element`: the colour the reference renderer
 * assigns to an element name. C and H are measured; the rest follow the
 * reference renderer's element colour names (alkali metals purple, alkaline
 * earths green, other metals silver, …). Elements missing here keep the
 * viewer's own element colour.
 */
export const VMD_NAME_COLORS: { [element: string]: string } = {
    H: '#FFFFFF', He: '#999999',
    Li: '#A600A6', Be: '#00FF00', B: '#FF9999', C: '#40BFBF', N: '#0000FF',
    O: '#FF0000', F: '#00FF00', Ne: '#999999',
    Na: '#A600A6', Mg: '#00FF00', Al: '#999999', Si: '#808033', P: '#808033',
    S: '#FFFF00', Cl: '#00FF00', Ar: '#999999',
    K: '#A600A6', Ca: '#00FF00', Sc: '#999999', Ti: '#999999', V: '#999999',
    Cr: '#999999', Mn: '#999999', Fe: '#FF8000', Co: '#FF8000', Ni: '#FF8000',
    Cu: '#FF8000', Zn: '#999999', Ga: '#999999', Ge: '#999999', As: '#00FF00',
    Se: '#FFFF00', Br: '#A62929', Kr: '#999999',
    Rb: '#A600A6', Sr: '#00FF00', Y: '#999999', Zr: '#999999', Nb: '#999999',
    Mo: '#999999', Ru: '#999999', Rh: '#999999', Pd: '#999999', Ag: '#999999',
    Cd: '#999999', In: '#999999', Sn: '#999999', Sb: '#999999', Te: '#FFFF00',
    I: '#A600A6', Xe: '#999999',
    Cs: '#A600A6', Ba: '#00FF00', La: '#999999', Ce: '#999999',
    Hf: '#999999', Ta: '#999999', W: '#999999', Re: '#999999', Os: '#999999',
    Ir: '#999999', Pt: '#999999', Au: '#FF8000', Hg: '#999999',
    Tl: '#999999', Pb: '#999999', Bi: '#999999', Th: '#999999', U: '#999999'
};

/**
 * van der Waals radii in Å, as reported by the reference renderer's per-atom
 * radius property. Used with the sphere scale of a CPK representation.
 */
export const VDW_RADII: { [element: string]: number } = {
    H: 1.00, He: 1.40, Li: 1.82, Be: 1.53, B: 1.92, C: 1.70, N: 1.55, O: 1.52,
    F: 1.47, Ne: 1.54, Na: 2.27, Mg: 1.73, Al: 1.84, Si: 2.10, P: 1.80,
    S: 1.80, Cl: 1.75, Ar: 1.88, K: 2.75, Ca: 2.31, Sc: 2.30, Ti: 2.15,
    V: 2.05, Cr: 2.05, Mn: 2.05, Fe: 2.05, Co: 2.00, Ni: 2.00, Cu: 2.00,
    Zn: 2.10, Ga: 1.87, Ge: 2.11, As: 1.85, Se: 1.90, Br: 1.85, Kr: 2.02,
    Rb: 3.03, Sr: 2.49, Y: 2.40, Zr: 2.30, Nb: 2.15, Mo: 2.10, Ru: 2.05,
    Rh: 2.00, Pd: 2.05, Ag: 2.10, Cd: 2.20, In: 2.20, Sn: 2.25, Sb: 2.20,
    Te: 2.06, I: 1.98, Xe: 2.16, Cs: 3.43, Ba: 2.68, La: 2.50, Ce: 2.48,
    Hf: 2.25, Ta: 2.20, W: 2.10, Re: 2.05, Os: 2.00, Ir: 2.00, Pt: 2.05,
    Au: 2.10, Hg: 2.05, Tl: 2.20, Pb: 2.30, Bi: 2.30, Th: 2.40, U: 2.30
};

export const DEFAULT_VDW_RADIUS = 1.50;

/**
 * `mol color Beta` gradient as the reference renderer actually renders it.
 * Not a plain red → green → blue ramp: the two passive channels sit on a
 * 15/255 floor and the active channel ramps at twice the nominal rate, so the
 * values here are what a rendered sphere really shows in the reference
 * renderer, not the nominal endpoints written into its colour scale.
 */
const GRADIENT_BASE = 15 / 255;
const GRADIENT_RATE = 2.0;

function toHex(channel: number): string {
    let v = Math.round(channel * 255);
    if (!isFinite(v)) v = 0;
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    const s = v.toString(16).toUpperCase();
    return s.length < 2 ? '0' + s : s;
}

function clamp01(v: number): number {
    if (!isFinite(v)) return 0;
    return v < 0 ? 0 : v > 1 ? 1 : v;
}

/**
 * Colour of a scalar on the reference renderer's RGB gradient, clamped to the
 * given range. A range that collapses to a single value paints everything with
 * the low end, which is what the reference renderer does.
 */
export function scalarColor(value: number, min: number, max: number): string {
    let t: number;
    if (!isFinite(min) || !isFinite(max) || max === min) {
        t = 0;
    } else {
        t = (value - min) / (max - min);
        if (!isFinite(t)) t = 0;
        t = clamp01(t);
    }
    let r: number, g: number, b: number;
    if (t <= 0.5) {
        r = clamp01(GRADIENT_BASE + GRADIENT_RATE * (0.5 - t));
        g = clamp01(GRADIENT_BASE + GRADIENT_RATE * t);
        b = GRADIENT_BASE;
    } else {
        r = GRADIENT_BASE;
        g = clamp01(GRADIENT_BASE + GRADIENT_RATE * (1 - t));
        b = clamp01(GRADIENT_BASE + GRADIENT_RATE * (t - 0.5));
    }
    return '#' + toHex(r) + toHex(g) + toHex(b);
}

/** Colour table entry for `mol color ColorID <id>`, honouring script overrides. */
export function vmdColor(colorId: number, scene?: VmdScene): string {
    let id = Math.floor(colorId);
    if (!isFinite(id)) id = 0;
    id = ((id % VMD_COLOR_IDS.length) + VMD_COLOR_IDS.length) % VMD_COLOR_IDS.length;
    if (scene && scene.colorOverrides && scene.colorOverrides[id]) {
        return scene.colorOverrides[id];
    }
    return VMD_COLOR_IDS[id] || VMD_FALLBACK_COLOR;
}

function normalizeElement(element: string): string {
    const e = (element || '').trim();
    if (!e) return 'C';
    return e.charAt(0).toUpperCase() + e.slice(1).toLowerCase();
}

/** Colour of an element under `mol color Name`, or null to keep the viewer's. */
export function vmdNameColor(element: string): string | null {
    return VMD_NAME_COLORS[normalizeElement(element)] || null;
}

export function vdwRadius(element: string): number {
    const r = VDW_RADII[normalizeElement(element)];
    return typeof r === 'number' ? r : DEFAULT_VDW_RADIUS;
}

/** How a single representation draws the selected atoms. */
export interface VmdRepKind {
    /** Draws atoms and/or bonds at all (false for pure surface styles). */
    atomic: boolean;
    /**
     * Sphere size as a multiplier on the viewer's own sphere radius (0 = the
     * representation draws no spheres). The script's scale is relative to the
     * reference renderer's default sphere scale, not an absolute radius: this
     * viewer draws ball-and-stick, so the script's *relative* sizes are kept
     * (`CPK 0.7` is smaller than `CPK 1.0`) without turning the structure into
     * a space-filling model.
     */
    sphereScale: number | null;
    /** Bond cylinder radius in Å, rebased the same way (null = no bonds). */
    bondRadius: number | null;
    /** Surface probe radius in Å — absolute, because a probe is a physical size. */
    probe: number | null;
}

const SURFACE_STYLES: { [style: string]: boolean } = {
    surf: true, surface: true, quicksurf: true, msms: true, isosurface: true
};

/**
 * Bond radius of the reference renderer's own default representation, and the
 * stick thickness this viewer uses by default: a script's bond radius is
 * rescaled by the ratio, so `CPK … 0.3` (the reference default) keeps the
 * viewer's usual sticks and `CPK … 0.1` draws thinner ones.
 */
const REFERENCE_BOND_RADIUS = 0.3;
const VIEWER_BOND_RADIUS = 0.12;

function arg(args: number[], index: number, fallback: number): number {
    const v = args[index];
    return typeof v === 'number' && isFinite(v) ? v : fallback;
}

function rebaseBondRadius(value: number): number {
    return Math.max(0, VIEWER_BOND_RADIUS * (value / REFERENCE_BOND_RADIUS));
}

/**
 * Translate `mol representation <style> <args…>` into drawing parameters.
 * Sizes are taken relative to the reference renderer's defaults (CPK: sphere
 * scale 1.0, bond radius 0.3) so an argument-less representation renders with
 * the viewer's normal ball-and-stick geometry.
 */
export function resolveRep(rep: VmdRep): VmdRepKind {
    const style = (rep.style || '').trim().toLowerCase();
    const args = rep.styleArgs || [];

    if (SURFACE_STYLES[style]) {
        return { atomic: false, sphereScale: null, bondRadius: null, probe: Math.max(0.01, arg(args, 0, 1.4)) };
    }
    if (style === 'vdw') {
        // Spheres only, at the viewer's scale — no sticks.
        return { atomic: true, sphereScale: Math.max(0, arg(args, 0, 1.0)), bondRadius: null, probe: null };
    }
    if (style === 'licorice') {
        return { atomic: true, sphereScale: 0, bondRadius: rebaseBondRadius(arg(args, 0, 0.3)), probe: null };
    }
    if (style === 'lines' || style === 'bonds') {
        // The argument of Lines/Bonds is a pixel width, not a radius.
        return { atomic: true, sphereScale: 0, bondRadius: 0.03, probe: null };
    }
    // CPK, and every style whose first two numbers are a sphere scale and a bond
    // radius (the reference renderer's uniform representation arguments).
    return {
        atomic: true,
        sphereScale: Math.max(0, arg(args, 0, 1.0)),
        bondRadius: rebaseBondRadius(arg(args, 1, REFERENCE_BOND_RADIUS)),
        probe: null
    };
}

export interface SelectionResult {
    /** Atom indices, or null for "the whole molecule". */
    indices: number[] | null;
    /** Set when the expression could not be resolved and was widened to "all". */
    warning: string | null;
}

function parseIndexList(list: string, atomCount: number): number[] | null {
    const tokens = list.replace(/,/g, ' ').split(/\s+/).filter(function (t) { return t !== ''; });
    const picked: number[] = [];
    let i = 0;
    while (i < tokens.length) {
        const start = parseInt(tokens[i], 10);
        if (isNaN(start)) return null;
        if (i + 1 < tokens.length && /^to$/i.test(tokens[i + 1])) {
            const end = i + 2 < tokens.length ? parseInt(tokens[i + 2], 10) : NaN;
            if (isNaN(end)) return null;
            // Expand only inside the structure: a typo like `index 0 to 1e9`
            // must not build a billion-entry list before the range check below.
            const lo = Math.max(0, Math.min(start, end));
            const hi = Math.min(atomCount - 1, Math.max(start, end));
            for (let k = lo; k <= hi; k++) picked.push(k);
            i += 3;
        } else {
            picked.push(start);
            i += 1;
        }
    }
    const out: number[] = [];
    for (let k = 0; k < picked.length; k++) {
        const n = picked[k];
        if (n >= 0 && n < atomCount && out.indexOf(n) < 0) out.push(n);
    }
    return out;
}

function parseElementList(list: string, atoms: { element: string }[]): number[] | null {
    const tokens = list.replace(/,/g, ' ').split(/\s+/).filter(function (t) { return t !== ''; });
    if (tokens.length === 0) return null;
    const wanted: { [name: string]: boolean } = {};
    for (let i = 0; i < tokens.length; i++) {
        const t = tokens[i];
        if (!/^[A-Za-z]{1,2}$/.test(t)) return null;
        wanted[t.toUpperCase()] = true;
    }
    const out: number[] = [];
    for (let i = 0; i < atoms.length; i++) {
        if (wanted[(atoms[i].element || '').toUpperCase()]) out.push(i);
    }
    return out;
}

/**
 * Resolve a `mol selection` expression to atom indices. `null` means the whole
 * molecule. Selections that need residue or serial metadata the viewer does not
 * keep are widened to the whole molecule and reported as a warning instead of
 * silently drawing the wrong atoms.
 */
export function resolveSelection(selection: string, atoms: { element: string }[]): SelectionResult {
    const raw = (selection || '').trim();
    if (raw === '' || /^all$/i.test(raw)) return { indices: null, warning: null };
    if (/^none$/i.test(raw)) return { indices: [], warning: null };

    let body = raw;
    let negate = false;
    const notMatch = body.match(/^not\s+([\s\S]+)$/i);
    if (notMatch) {
        negate = true;
        body = notMatch[1].trim();
    }

    // Only a single clause is resolved; a compound expression keeps its first
    // clause (the widest one) and says so, rather than dropping the rep.
    let warning: string | null = null;
    const compound = body.match(/^([\s\S]*?)\s+(and|or)\s+[\s\S]+$/i);
    if (compound) {
        body = compound[1].trim();
        warning = 'compound selection "' + raw + '" simplified to "' + body + '"';
    }

    const clause = body.match(/^([A-Za-z_]+)\s+([\s\S]+)$/);
    if (!clause) {
        return { indices: null, warning: 'selection "' + raw + '" is not understood — drawn as the whole molecule' };
    }
    const keyword = clause[1].toLowerCase();
    const list = clause[2].trim();

    let picked: number[] | null;
    if (keyword === 'index') {
        picked = parseIndexList(list, atoms.length);
    } else if (keyword === 'element' || keyword === 'name' || keyword === 'type') {
        picked = parseElementList(list, atoms);
    } else {
        picked = null;
    }

    if (picked === null) {
        return { indices: null, warning: 'selection "' + raw + '" needs residue metadata the viewer does not keep — drawn as the whole molecule' };
    }
    if (!negate) return { indices: picked, warning: warning };

    const excluded: { [index: number]: boolean } = {};
    for (let i = 0; i < picked.length; i++) excluded[picked[i]] = true;
    const rest: number[] = [];
    for (let i = 0; i < atoms.length; i++) if (!excluded[i]) rest.push(i);
    return { indices: rest, warning: warning };
}

/** A translucent isosurface to build around a set of atoms. */
export interface VmdSurfaceSegment {
    color: string;
    indices: number[];
    /** Probe radius in Å: the surface sits this far outside the van der Waals surface. */
    probe: number;
    label: string;
}

export interface VmdResolvedScene {
    /** Per-atom colour, null = keep the element colour. */
    atomColors: (string | null)[];
    /** Per-atom sphere size as a multiplier on the viewer's own radius, null = viewer default, 0 = not drawn. */
    atomScales: (number | null)[];
    /** Per-bond cylinder radius in Å, null = the viewer's own radius, 0 = not drawn. */
    bondRadii: (number | null)[];
    segments: VmdSurfaceSegment[];
    warnings: string[];
}

interface AtomLike {
    element: string;
    beta?: number;
    occupancy?: number;
}

interface BondLike {
    atom1: number;
    atom2: number;
}

type ColorResolver = (atom: AtomLike) => string | null;

/**
 * Build the colour function of one representation.
 *
 * `hasExplicitRange` is set when the script wrote both ends of the colour scale
 * (`mol scaleminmax`). That is an unconditional statement about the column: the
 * representation is meant to be read as a scalar channel, so the colouring is
 * applied even when every value in the column is the same — the reference
 * renderer then paints the whole molecule with the middle of the ramp. The
 * column's own values only decide *what* the layout looks like when the script
 * leaves the range open, in which case the range comes from the data the
 * representation selects.
 */
function makeColorResolver(rep: VmdRep, indices: number[], atoms: AtomLike[], scene: VmdScene, hasExplicitRange: boolean): ColorResolver | null {
    const mode = (rep.colorMode || '').trim().toLowerCase();

    if (mode === 'colorid') {
        const color = vmdColor(typeof rep.colorArg === 'number' ? rep.colorArg : 0, scene);
        return function () { return color; };
    }
    if (mode === 'name' || mode === 'element') {
        return function (atom: AtomLike) { return vmdNameColor(atom.element); };
    }
    if (mode === 'beta' || mode === 'occupancy') {
        const field: 'beta' | 'occupancy' = mode === 'beta' ? 'beta' : 'occupancy';
        const writtenMin = typeof rep.scaleMin === 'number' ? rep.scaleMin : null;
        const writtenMax = typeof rep.scaleMax === 'number' ? rep.scaleMax : null;
        let min: number;
        let max: number;
        if (hasExplicitRange && writtenMin !== null && writtenMax !== null) {
            // The script's own scale, used as written — including for a column
            // whose values are all identical.
            min = writtenMin;
            max = writtenMax;
        } else {
            // No written range: the layout comes from what the selection holds.
            let lo = Infinity, hi = -Infinity;
            for (let k = 0; k < indices.length; k++) {
                const v = atoms[indices[k]][field];
                if (typeof v === 'number' && isFinite(v)) {
                    if (v < lo) lo = v;
                    if (v > hi) hi = v;
                }
            }
            if (lo === Infinity) return null;
            min = writtenMin === null ? lo : writtenMin;
            max = writtenMax === null ? hi : writtenMax;
        }
        return function (atom: AtomLike) {
            const v = atom[field];
            if (typeof v !== 'number' || !isFinite(v)) return null;
            // Always (value - min) / (max - min) against the range in use: a
            // column that happens to be constant still sits at its own place on
            // a written scale (a written min === max is the only collapsed case,
            // and scalarColor then paints the low end).
            return scalarColor(v, min, max);
        };
    }
    return null;
}

/**
 * Turn one molecule's representation list into concrete drawing properties for
 * the atoms and bonds of that molecule.
 *
 * Representations are applied in file order. For an atom, the representation
 * drawing the largest sphere decides both its size and its colour — that is the
 * one visible in the composed image when two representations nest. For a bond,
 * the last representation covering both of its endpoints decides.
 */
export function resolveVmdScene(
    atoms: AtomLike[],
    bonds: BondLike[],
    molecule: VmdMolecule,
    scene: VmdScene
): VmdResolvedScene {
    const n = atoms.length;
    const atomColors: (string | null)[] = new Array(n).fill(null);
    const atomScales: (number | null)[] = new Array(n).fill(null);
    const bondRadii: (number | null)[] = new Array(bonds.length).fill(null);
    const segments: VmdSurfaceSegment[] = [];
    const warnings: string[] = [];

    const dominantScale: number[] = new Array(n).fill(-1);

    const wholeMolecule: number[] = [];
    for (let i = 0; i < n; i++) wholeMolecule.push(i);

    for (let ri = 0; ri < molecule.reps.length; ri++) {
        const rep = molecule.reps[ri];
        const kind = resolveRep(rep);
        const selection = resolveSelection(rep.selection, atoms);
        if (selection.warning) warnings.push('rep ' + ri + ': ' + selection.warning);
        const indices = selection.indices === null ? wholeMolecule : selection.indices;
        // A script that writes an explicit colour-scale range is stating that
        // the file's temperature-factor column is a scalar channel, even when
        // the values happen to be constant — the reference renderer then paints
        // the whole molecule with the middle of the ramp. The data-derived
        // range is only a fallback for a representation that leaves the range
        // open.
        const hasExplicitRange =
            typeof rep.scaleMin === 'number' && typeof rep.scaleMax === 'number';
        // Built for every representation, surfaces included: a surface patch
        // takes its colour from the representation that asks for it.
        const resolver = makeColorResolver(rep, indices, atoms, scene, hasExplicitRange);

        if (kind.atomic && kind.sphereScale !== null && kind.sphereScale > 0) {
            for (let k = 0; k < indices.length; k++) {
                const idx = indices[k];
                if (idx < 0 || idx >= n) continue;
                const scale = kind.sphereScale;
                if (scale >= dominantScale[idx]) {
                    dominantScale[idx] = scale;
                    atomScales[idx] = scale;
                    if (resolver) atomColors[idx] = resolver(atoms[idx]);
                }
            }
        } else if (kind.atomic && kind.sphereScale !== null) {
            // Bond-only styles (Licorice, Lines): no spheres, but they still
            // contribute a colour where no sphere is drawn at all.
            for (let k = 0; k < indices.length; k++) {
                const idx = indices[k];
                if (idx < 0 || idx >= n) continue;
                if (dominantScale[idx] <= 0) {
                    dominantScale[idx] = 0;
                    atomScales[idx] = 0;
                    if (resolver) atomColors[idx] = resolver(atoms[idx]);
                }
            }
        }

        if (kind.atomic && kind.bondRadius !== null && kind.bondRadius > 0) {
            const inSelection: { [index: number]: boolean } = {};
            for (let k = 0; k < indices.length; k++) inSelection[indices[k]] = true;
            for (let bi = 0; bi < bonds.length; bi++) {
                const b = bonds[bi];
                if (inSelection[b.atom1] && inSelection[b.atom2]) bondRadii[bi] = kind.bondRadius;
            }
        }

        if (kind.probe !== null) {
            let color = VMD_FALLBACK_COLOR;
            if (resolver && indices.length > 0) {
                // One colour has to stand for the whole patch. For a
                // scalar-driven colour mode that is the colour of the patch's
                // mean scalar; for every other mode the colour is constant.
                const clone: AtomLike = { element: atoms[indices[0]].element };
                let betaSum = 0, betaCount = 0, occSum = 0, occCount = 0;
                for (let k = 0; k < indices.length; k++) {
                    const a = atoms[indices[k]];
                    if (typeof a.beta === 'number' && isFinite(a.beta)) { betaSum += a.beta; betaCount++; }
                    if (typeof a.occupancy === 'number' && isFinite(a.occupancy)) { occSum += a.occupancy; occCount++; }
                }
                if (betaCount > 0) clone.beta = betaSum / betaCount;
                if (occCount > 0) clone.occupancy = occSum / occCount;
                color = resolver(clone) || resolver(atoms[indices[0]]) || color;
            } else if (resolver) {
                color = resolver(atoms[0]) || color;
            }
            segments.push({ color: color, indices: indices, probe: kind.probe, label: rep.style + ' ' + (ri + 1) });
        }
    }

    // What the script does not draw is not drawn: atoms no sphere-drawing
    // representation covers get scale 0 (hidden), bonds no bond-drawing
    // representation covers get radius 0 as well. A molecule without any
    // representation keeps the viewer's own defaults, so a script that only
    // names a structure still shows it.
    if (molecule.reps.length > 0) {
        for (let i = 0; i < n; i++) {
            if (dominantScale[i] < 0) {
                dominantScale[i] = 0;
                atomScales[i] = 0;
                atomColors[i] = null;
            }
        }
        for (let i = 0; i < bonds.length; i++) {
            if (bondRadii[i] === null) bondRadii[i] = 0;
        }
    }

    return {
        atomColors: atomColors,
        atomScales: atomScales,
        bondRadii: bondRadii,
        segments: segments,
        warnings: warnings
    };
}

/**
 * Parser for VMD visualization scripts (`mol new` + `mol representation` /
 * `mol color` / `mol selection` / `mol addrep` …).
 *
 * The script is a display recipe for a structure stored elsewhere, so parsing
 * produces two things: the referenced structure file(s) to load, and the full
 * scene description — every molecule, every representation with its style,
 * colour mode, selection, material and explicit colour-scale range, plus the
 * display background and any palette overrides.
 */

const DISPLAY_BACKGROUND: { [name: string]: string } = {
    white: '#FFFFFF',
    black: '#000000',
    gray: '#808080',
    grey: '#808080',
    blue: '#0000FF',
    red: '#FF0000',
    green: '#00FF00'
};

export interface TclParseResult {
    /** Structure file referenced by the first `mol new`. */
    sourceFile: string;
    /** `type` argument of the first `mol new`. */
    sourceType: string;
    /**
     * Legacy per-group colour overrides: one entry per representation that both
     * picks a ColorID colour and selects atoms by explicit index. Kept because
     * the structure-diff view colours its atoms from this list.
     */
    groups: AtomGroup[];
    /** Complete display recipe. */
    scene: VmdScene;
}

/** Compatibility accessor for the ColorID palette (modulo 33). */
export function getVmdColor(colorId: number): string {
    return vmdColor(colorId);
}

/**
 * Split a script line into words, honouring `{…}` and `"…"` grouping. Inside a
 * group a backslash escapes the next character; outside one it only escapes
 * whitespace, so Windows-style paths survive intact.
 */
function tokenize(line: string): string[] {
    const out: string[] = [];
    let i = 0;
    while (i < line.length) {
        const ch = line[i];
        if (ch === ' ' || ch === '\t') { i++; continue; }
        if (ch === '{' || ch === '"') {
            // Braces nest; a double-quoted word ends at the next unescaped quote.
            const closer = ch === '{' ? '}' : '"';
            let depth = 1;
            let j = i + 1;
            let buf = '';
            while (j < line.length) {
                const d = line[j];
                if (d === '\\' && j + 1 < line.length) {
                    const n = line[j + 1];
                    buf += n === 'n' ? '\n' : n === 't' ? '\t' : n === 'r' ? '\r' : n;
                    j += 2;
                    continue;
                }
                if (ch === '{' && d === '{') depth++;
                else if (d === closer) {
                    depth--;
                    if (depth === 0) break;
                }
                buf += d;
                j++;
            }
            out.push(buf);
            i = j + 1;
            continue;
        }
        let j = i;
        let buf = '';
        while (j < line.length && !/\s/.test(line[j])) {
            if (line[j] === '\\' && j + 1 < line.length && /\s/.test(line[j + 1])) {
                buf += line[j + 1];
                j += 2;
                continue;
            }
            buf += line[j];
            j++;
        }
        out.push(buf);
        i = j;
    }
    return out;
}

function toHexByte(v: number): string {
    let n = Math.round(v * 255);
    if (!isFinite(n)) n = 0;
    if (n < 0) n = 0;
    if (n > 255) n = 255;
    const s = n.toString(16).toUpperCase();
    return s.length < 2 ? '0' + s : s;
}

function newRep(): VmdRep {
    return { style: 'CPK', styleArgs: [], colorMode: 'Name', selection: 'all', material: 'Opaque' };
}

function parseRep(rep: VmdRep, tokens: string[], from: number): void {
    rep.style = tokens[from];
    const args: number[] = [];
    for (let i = from + 1; i < tokens.length; i++) {
        const v = parseFloat(tokens[i]);
        if (!isNaN(v)) args.push(v);
    }
    rep.styleArgs = args;
}

export function parseTcl(content: string): TclParseResult {
    // A byte-order mark rides in front of the first command in files saved by
    // some editors and would turn "mol new …" into an unrecognisable word.
    const lines = content.replace(/^\uFEFF/, '').split(/\r?\n/);
    const scene: VmdScene = { molecules: [], colorScaleMethod: 'RGB' };
    const groups: AtomGroup[] = [];

    let sourceFile = '';
    let sourceType = '';
    let current: VmdMolecule | null = null;
    let pending: VmdRep | null = null;
    // `mol new` also creates one default representation, which scripts almost
    // always delete again with `mol delrep 0 top`. Keeping track of it is what
    // makes representation indices address the same representation as in the
    // reference renderer: with the default still present the first `addrep`
    // is number 1, without it it is number 0.
    let defaultRepExists = false;
    const repIndex = function (vmdIndex: number): number {
        return defaultRepExists ? vmdIndex - 1 : vmdIndex;
    };

    for (let li = 0; li < lines.length; li++) {
        const trimmed = lines[li].replace(/^\s+/, '');
        if (trimmed === '' || trimmed.charAt(0) === '#') continue;
        const tokens = tokenize(trimmed);
        if (tokens.length === 0) continue;
        const verb = tokens[0].toLowerCase();

        if (verb === 'color') {
            const sub = (tokens[1] || '').toLowerCase();
            if (sub === 'display' && (tokens[2] || '').toLowerCase() === 'background') {
                const name = (tokens[3] || '').toLowerCase();
                if (DISPLAY_BACKGROUND[name]) {
                    scene.background = DISPLAY_BACKGROUND[name];
                } else if (/^#[0-9a-f]{6}$/i.test(tokens[3] || '')) {
                    scene.background = tokens[3].toUpperCase();
                } else if (tokens.length >= 6) {
                    const r = parseFloat(tokens[3]), g = parseFloat(tokens[4]), b = parseFloat(tokens[5]);
                    if (!isNaN(r) && !isNaN(g) && !isNaN(b)) {
                        scene.background = '#' + toHexByte(r) + toHexByte(g) + toHexByte(b);
                    }
                }
            } else if (sub === 'scale' && (tokens[2] || '').toLowerCase() === 'method') {
                scene.colorScaleMethod = tokens[3] || 'RGB';
            } else if (sub === 'change' && (tokens[2] || '').toLowerCase() === 'rgb') {
                const id = parseInt(tokens[3], 10);
                const r = parseFloat(tokens[4]), g = parseFloat(tokens[5]), b = parseFloat(tokens[6]);
                if (!isNaN(id) && !isNaN(r) && !isNaN(g) && !isNaN(b)) {
                    if (!scene.colorOverrides) scene.colorOverrides = {};
                    const size = VMD_COLOR_IDS.length;
                    scene.colorOverrides[((id % size) + size) % size] = '#' + toHexByte(r) + toHexByte(g) + toHexByte(b);
                }
            }
            continue;
        }

        if (verb !== 'mol') continue;
        const op = (tokens[1] || '').toLowerCase();

        if (op === 'new') {
            const file = tokens[2] || '';
            let type = '';
            for (let i = 3; i < tokens.length - 1; i++) {
                if (tokens[i].toLowerCase() === 'type') { type = tokens[i + 1]; break; }
            }
            current = { sourceFile: file, sourceType: type, label: file, reps: [] };
            scene.molecules.push(current);
            defaultRepExists = true;
            if (!sourceFile) {
                sourceFile = file;
                sourceType = type;
            }
            pending = null;
            continue;
        }

        if (op === 'addrep') {
            if (current) current.reps.push(pending || newRep());
            pending = null;
            continue;
        }

        if (op === 'delrep') {
            const idx = parseInt(tokens[2], 10);
            if (current && !isNaN(idx) && idx >= 0) {
                if (defaultRepExists) {
                    if (idx === 0) defaultRepExists = false;
                    else if (idx - 1 < current.reps.length) current.reps.splice(idx - 1, 1);
                } else if (idx < current.reps.length) {
                    current.reps.splice(idx, 1);
                }
            }
            continue;
        }

        if (op === 'representation') {
            if (tokens.length < 3) continue;
            if (!pending) pending = newRep();
            parseRep(pending, tokens, 2);
            continue;
        }

        if (op === 'color') {
            if (!pending) pending = newRep();
            pending.colorMode = tokens[2] || 'Name';
            if (tokens.length >= 4) {
                const v = parseInt(tokens[3], 10);
                if (!isNaN(v)) pending.colorArg = v;
            }
            continue;
        }

        if (op === 'selection') {
            if (!pending) pending = newRep();
            pending.selection = tokens.slice(2).join(' ').trim() || 'all';
            continue;
        }

        if (op === 'material') {
            if (!pending) pending = newRep();
            pending.material = tokens[2] || 'Opaque';
            continue;
        }

        if (op === 'scaleminmax') {
            const tail = tokens.slice(3);
            if (current && tail.length >= 2) {
                const min = parseFloat(tail[tail.length - 2]);
                const max = parseFloat(tail[tail.length - 1]);
                let target = current.reps.length - 1;
                if (tail.length >= 3) {
                    const written = parseInt(tail[tail.length - 3], 10);
                    if (!isNaN(written)) target = repIndex(written);
                }
                if (!isNaN(min) && !isNaN(max) && target >= 0 && target < current.reps.length) {
                    current.reps[target].scaleMin = min;
                    current.reps[target].scaleMax = max;
                }
            }
            continue;
        }

        // Rep modifications address an existing representation: `mol modstyle
        // <rep> <mol> <style> args…`, `mol modcolor <rep> <mol> <mode> [arg]`,
        // `mol modselect <rep> <mol> <selection…>`, `mol modmaterial …`.
        if (op === 'modstyle' || op === 'modcolor' || op === 'modselect' || op === 'modmaterial') {
            if (!current) continue;
            const written = parseInt(tokens[2], 10);
            if (isNaN(written)) continue;
            const idx = repIndex(written);
            if (idx < 0 || idx >= current.reps.length) continue;
            const rep = current.reps[idx];
            if (op === 'modstyle' && tokens.length >= 5) {
                parseRep(rep, tokens, 4);
            } else if (op === 'modcolor') {
                rep.colorMode = tokens[4] || rep.colorMode;
                if (tokens.length >= 6) {
                    const v = parseInt(tokens[5], 10);
                    if (!isNaN(v)) rep.colorArg = v;
                }
            } else if (op === 'modselect') {
                rep.selection = tokens.slice(4).join(' ').trim() || rep.selection;
            } else if (op === 'modmaterial' && tokens.length >= 5) {
                rep.material = tokens[4];
            }
            continue;
        }
    }

    // Legacy group list: only representations that name a ColorID colour and
    // select atoms by explicit index, in file order. Selections that need
    // residue metadata are skipped here (they are handled by the scene path).
    for (let mi = 0; mi < scene.molecules.length; mi++) {
        const molecule = scene.molecules[mi];
        for (let ri = 0; ri < molecule.reps.length; ri++) {
            const rep = molecule.reps[ri];
            if ((rep.colorMode || '').toLowerCase() !== 'colorid') continue;
            const indices = parseRawIndexList(rep.selection);
            if (indices === null || indices.length === 0) continue;
            const colorId = typeof rep.colorArg === 'number' ? rep.colorArg : 0;
            groups.push({ colorId: colorId, color: vmdColor(colorId, scene), indices: indices });
        }
    }

    return { sourceFile: sourceFile, sourceType: sourceType, groups: groups, scene: scene };
}

/**
 * Atom indices written literally in an `index …` selection, without an atom
 * count to validate against (the structure is not loaded yet when the script is
 * parsed). Returns null when the expression is not an index list.
 */
function parseRawIndexList(selection: string): number[] | null {
    const body = (selection || '').trim().replace(/^not\s+/i, '');
    const clause = body.match(/^index\s+([\s\S]+)$/i);
    if (!clause) return null;
    const tokens = clause[1].replace(/,/g, ' ').split(/\s+/).filter(function (t) { return t !== ''; });
    const picked: number[] = [];
    let i = 0;
    while (i < tokens.length) {
        const start = parseInt(tokens[i], 10);
        if (isNaN(start)) return null;
        if (i + 1 < tokens.length && /^to$/i.test(tokens[i + 1])) {
            const end = i + 2 < tokens.length ? parseInt(tokens[i + 2], 10) : NaN;
            if (isNaN(end)) return null;
            // No atom count to clamp against here, so an absurd range is
            // rejected outright instead of being expanded.
            const lo = Math.min(start, end), hi = Math.max(start, end);
            if (hi - lo > 1000000) return null;
            for (let k = lo; k <= hi; k++) picked.push(k);
            i += 3;
        } else {
            picked.push(start);
            i += 1;
        }
    }
    return picked;
}

/**
 * Record types that mark a real crystallographic PDB (deposited structures,
 * files exported with an experimental header). A file containing any record
 * outside the bare coordinate set below is treated as such.
 */
const COORDINATE_ONLY_RECORDS: { [record: string]: boolean } = {
    ATOM: true, HETATM: true, CONECT: true, TER: true, END: true,
    MODEL: true, ENDMDL: true, MASTER: true
};

/**
 * True when a PDB is a bare coordinate file whose temperature-factor column
 * carries a per-atom scalar rather than a crystallographic B-factor.
 *
 * The visualization output of this toolchain is exactly that: a PDB with no
 * header records at all, holding a unit index, a per-atom force error or a
 * round-to-round difference in columns 61-66 (see `test/samples/SPEC.md`), so
 * that opening the file — with or without the script that wrote it — shows the
 * structure coloured by that scalar. A file that does carry header records, or
 * whose column is constant, is left alone.
 */
export function hasScalarChannel(content: string, atoms: { beta?: number }[]): boolean {
    const lines = content.split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
        const record = lines[i].substring(0, 6).trim().toUpperCase();
        if (record === '') continue;
        if (!COORDINATE_ONLY_RECORDS[record]) return false;
    }
    let min = Infinity;
    let max = -Infinity;
    let count = 0;
    for (let i = 0; i < atoms.length; i++) {
        const beta = atoms[i].beta;
        if (typeof beta !== 'number' || !isFinite(beta)) continue;
        count++;
        if (beta < min) min = beta;
        if (beta > max) max = beta;
    }
    return count >= 2 && max > min;
}
