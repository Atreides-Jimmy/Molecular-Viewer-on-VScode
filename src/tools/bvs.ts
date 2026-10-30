/**
 * Bond-valence sums (BVS) for a periodic structure.
 *
 * The bond valence of one contact is Brown and Altermatt's exponential
 *   s = exp((R0 - d) / B)
 * with R0 and B tabulated per cation-anion pair. Summing every contact of a
 * site gives its bond-valence sum, which an ionic crystal predicts to equal
 * the site's formal valence: this is the read-out that says whether a site is
 * under- or over-coordinated without guessing bonds from distances.
 *
 * The module is read-only: it touches nothing but the coordinates and the cell
 * it is handed, and reports numbers. The viewer's bond list is not consulted
 * and not modified.
 *
 * ---------------------------------------------------------------------------
 * Parameter data
 * ---------------------------------------------------------------------------
 * The table below is derived from I. David Brown's accumulated bond-valence
 * parameter file bvparm2013.cif (2013-09-11), available from the IUCr:
 *   http://www.iucr.org/resources/data/datasets/bond-valence-parameters
 *
 * COPYRIGHT NOTICE, reproduced as the data licence requires:
 *   This table may be used and distributed without fee for non-profit
 *   purposes providing (1) that this copyright notice is included and (2) no
 *   fee is charged for the table and (3) details of any changes made in this
 *   list by anyone other than the copyright owner are suitably noted. Please
 *   consult the copyright owner regarding any other uses. The copyright is
 *   owned by I. David Brown, Brockhouse Institute for Materials Research,
 *   McMaster University, Hamilton, Ontario, Canada.
 *
 * Changes made here: the table keeps only the cation states and anion
 * elements listed below, drops rows flagged 'unchecked' and rows fitted to
 * gas-phase or small-sample data, and keeps one fit per pair. The preference
 * is Brown and Altermatt (1985, reference a), then Brese and O'Keeffe (1991,
 * reference b), then other published fits, then the fit closest to B = 0.37.
 * Two pairs use Brese and O'Keeffe's Table 2 values directly because the
 * accumulated file lists no such fit (Al-O) or only a superseded one (Si-O):
 * Si-O 1.624 and Al-O 1.651, both with B = 0.37.
 *
 * Scope: main-group cations and the closed-shell cations Ag(I), Zn(II) and
 * Cd(II). Transition metals with several accessible oxidation states and
 * variable spin states are deliberately absent: a BVS built on one chosen
 * state of such a metal is not a check, it is a guess. The anion centres are
 * the usual closed-shell anions; their sums use the cation parameters of
 * their neighbours, with each neighbour's oxidation state chosen as the one
 * whose own sum best matches it.
 *
 * Both the extension host and the webview run the same code: BVS_SOURCE
 * serialises this module into one expression the page evaluates, exactly as
 * the bond-perception module does, so the report cannot drift between them.
 */

/**
 * A minimal atom for the analysis: what the sum needs, nothing more. The
 * viewer's Atom is assignable to it.
 */
export interface BvsAtom {
    element: string;
    x: number;
    y: number;
    z: number;
    ghost?: boolean;
}

export interface BvsAtomResult {
    kind: 'cation' | 'anion';
    element: string;
    /** Cation oxidation state chosen for this site; null on an anion site. */
    state: number | null;
    /** The valence the sum is checked against: state for a cation, charge for an anion. */
    valence: number | null;
    /** Bond-valence sum over every periodic image inside the pair cutoffs. */
    sum: number | null;
    /** |sum - valence|. */
    delta: number | null;
    /** Number of contacts that contributed. */
    contacts: number;
    status: 'ok' | 'no-contacts' | 'no-params';
}

export interface BvsCounts {
    atoms: number;
    ok: number;
    noParams: number;
    noContacts: number;
    skipped: number;
    cations: number;
    anions: number;
}

export interface BvsReport {
    ok: boolean;
    reason: string;
    parameters: string;
    maxCut: number;
    /** One entry per input atom, in input order; null for a non-centre or a ghost. */
    perAtom: (BvsAtomResult | null)[];
    warnings: string[];
    counts: BvsCounts;
}

// ---------------------------------------------------------------------------
// 0. The parameter table
// ---------------------------------------------------------------------------

/**
 * R0 and B per cation state and anion element. Every pair uses the same B in
 * the curated set (0.37 A) except two Ga/Bi fits where the published value is
 * 0.35; the pair function reads both numbers rather than assuming 0.37.
 */
const BVS_CATIONS: { [element: string]: { [state: string]: { [anion: string]: [number, number] } } } = {
    Li: {
        '1': { N: [1.61, 0.37], O: [1.466, 0.37], F: [1.36, 0.37], S: [1.94, 0.37], Cl: [1.91, 0.37], Se: [2.09, 0.37], Br: [2.02, 0.37], Te: [2.3, 0.37], I: [2.22, 0.37] },
    },
    Be: {
        '2': { N: [1.5, 0.37], O: [1.381, 0.37], F: [1.281, 0.37], S: [1.83, 0.37], Cl: [1.76, 0.37], Se: [1.97, 0.37], Br: [1.9, 0.37], Te: [2.21, 0.37], I: [2.1, 0.37] },
    },
    B: {
        '3': { N: [1.47, 0.37], O: [1.371, 0.37], F: [1.281, 0.37], S: [1.82, 0.37], Cl: [1.74, 0.37], Se: [1.95, 0.37], Br: [1.88, 0.37], Te: [2.2, 0.37], I: [2.1, 0.37] },
    },
    Na: {
        '1': { N: [1.93, 0.37], O: [1.803, 0.37], F: [1.677, 0.37], S: [2.3, 0.37], Cl: [2.15, 0.37], Se: [2.41, 0.37], Br: [2.33, 0.37], Te: [2.64, 0.37], I: [2.56, 0.37] },
    },
    Mg: {
        '2': { N: [1.85, 0.37], O: [1.693, 0.37], F: [1.578, 0.37], S: [2.18, 0.37], Cl: [2.08, 0.37], Se: [2.32, 0.37], Br: [2.28, 0.37], Te: [2.53, 0.37], I: [2.46, 0.37] },
    },
    Al: {
        '3': { N: [1.79, 0.37], O: [1.651, 0.37], F: [1.545, 0.37], S: [2.13, 0.37], Cl: [2.032, 0.37], Se: [2.27, 0.37], Br: [2.2, 0.37], Te: [2.48, 0.37], I: [2.41, 0.37] },
    },
    Si: {
        '4': { N: [1.724, 0.37], O: [1.624, 0.37], F: [1.58, 0.37], S: [2.126, 0.37], Cl: [2.03, 0.37], Se: [2.26, 0.37], Br: [2.2, 0.37], Te: [2.49, 0.37], I: [2.41, 0.37] },
    },
    P: {
        '5': { N: [1.704, 0.37], O: [1.617, 0.37], S: [2.145, 0.37], Cl: [2.02, 0.37] },
    },
    K: {
        '1': { N: [2.26, 0.37], O: [2.132, 0.37], F: [1.992, 0.37], S: [2.59, 0.37], Cl: [2.519, 0.37], Se: [2.72, 0.37], Br: [2.66, 0.37], Te: [2.93, 0.37], I: [2.88, 0.37] },
    },
    Ca: {
        '2': { N: [2.14, 0.37], O: [1.967, 0.37], F: [1.842, 0.37], S: [2.45, 0.37], Cl: [2.37, 0.37], Se: [2.56, 0.37], Br: [2.49, 0.37], Te: [2.76, 0.37], I: [2.72, 0.37] },
    },
    Ga: {
        '3': { O: [1.73, 0.37], F: [1.62, 0.37], S: [2.163, 0.37], Cl: [2.07, 0.37], Br: [2.2, 0.35] },
    },
    Ge: {
        '4': { O: [1.748, 0.37], F: [1.66, 0.37], S: [2.217, 0.37], Cl: [2.14, 0.37] },
    },
    As: {
        '3': { O: [1.789, 0.37], F: [1.7, 0.37], S: [2.272, 0.37], Cl: [2.16, 0.37] },
        '5': { O: [1.767, 0.37], F: [1.62, 0.37] },
    },
    Rb: {
        '1': { N: [2.37, 0.37], O: [2.263, 0.37], F: [2.16, 0.37], S: [2.7, 0.37], Cl: [2.652, 0.37], Se: [2.81, 0.37], Br: [2.78, 0.37], Te: [3, 0.37], I: [3.01, 0.37] },
    },
    Sr: {
        '2': { N: [2.23, 0.37], O: [2.118, 0.37], F: [2.019, 0.37], S: [2.59, 0.37], Cl: [2.51, 0.37], Se: [2.72, 0.37], Br: [2.68, 0.37], Te: [2.87, 0.37], I: [2.88, 0.37] },
    },
    In: {
        '3': { O: [1.902, 0.37], F: [1.792, 0.37], S: [2.37, 0.37], Cl: [2.28, 0.37] },
    },
    Sn: {
        '2': { O: [1.984, 0.37], F: [1.925, 0.37], Cl: [2.36, 0.37] },
        '4': { O: [1.905, 0.37], F: [1.843, 0.37], S: [2.399, 0.37], Cl: [2.276, 0.37] },
    },
    Sb: {
        '3': { N: [2.108, 0.37], O: [1.955, 0.37], F: [1.883, 0.37], S: [2.474, 0.37], Cl: [2.35, 0.37] },
        '5': { O: [1.942, 0.37], F: [1.797, 0.37], Cl: [2.3, 0.37] },
    },
    Cs: {
        '1': { N: [2.53, 0.37], O: [2.417, 0.37], F: [2.33, 0.37], S: [2.89, 0.37], Cl: [2.791, 0.37], Se: [2.98, 0.37], Br: [2.95, 0.37], Te: [3.16, 0.37], I: [3.18, 0.37] },
    },
    Ba: {
        '2': { N: [2.47, 0.37], O: [2.285, 0.37], F: [2.188, 0.37], S: [2.769, 0.37], Cl: [2.69, 0.37], Se: [2.88, 0.37], Br: [2.88, 0.37], Te: [3.08, 0.37], I: [3.13, 0.37] },
    },
    Tl: {
        '1': { O: [2.124, 0.37], F: [2.15, 0.37], S: [2.545, 0.37], Cl: [2.56, 0.37], I: [2.822, 0.37] },
        '3': { O: [2.003, 0.37], F: [1.88, 0.37], Cl: [2.32, 0.37] },
    },
    Pb: {
        '2': { O: [2.112, 0.37], F: [2.03, 0.37], S: [2.541, 0.37], Cl: [2.53, 0.37] },
        '4': { O: [2.042, 0.37], F: [1.94, 0.37], Cl: [2.43, 0.37] },
    },
    Bi: {
        '3': { O: [2.094, 0.37], F: [1.99, 0.37], S: [2.57, 0.37], Cl: [2.48, 0.37], Br: [2.597, 0.37] },
        '5': { O: [2.06, 0.37], F: [1.97, 0.37], Cl: [2.44, 0.37] },
    },
    Ag: {
        '1': { O: [1.842, 0.37], F: [1.8, 0.37], S: [2.119, 0.37], Cl: [2.09, 0.37] },
    },
    Zn: {
        '2': { N: [1.77, 0.37], O: [1.704, 0.37], F: [1.62, 0.37], S: [2.09, 0.37], Cl: [2.01, 0.37], Se: [2.22, 0.37], Br: [2.15, 0.37], Te: [2.45, 0.37], I: [2.36, 0.37] },
    },
    Cd: {
        '2': { N: [1.96, 0.37], O: [1.904, 0.37], F: [1.811, 0.37], S: [2.304, 0.37], Cl: [2.212, 0.37], Se: [2.4, 0.37], Br: [2.35, 0.37], Te: [2.59, 0.37], I: [2.57, 0.37] },
    }
};

/** Formal charge magnitudes of the anion centres the sum knows. */
const BVS_ANIONS: { [element: string]: number } = {
    N: 3, O: 2, F: 1, S: 2, Cl: 1, Se: 2, Br: 1, Te: 2, I: 1
};

/**
 * The longest distance any pair can contribute from. A contact is inside its
 * own pair cutoff when s = exp((R0 - d) / B) is still worth summing; the
 * per-pair cutoff below is R0 + 3.3 B, where one contact contributes 0.037.
 * This ceiling only bounds the lattice search.
 */
export const BVS_MAX_CUT = 4.0;

/** How far past R0 a contact is still counted, in units of B. */
const BVS_TAIL = 3.3;

/**
 * Above this the all-pairs lattice search is no longer an on-demand read: a
 * framework cell of a few thousand sites makes the report a background job,
 * and there is no progress to show for it. The report says so instead.
 */
export const BVS_MAX_ATOMS = 1500;

/** How many lattice repeats one axis may need; a guard against a tiny cell. */
const BVS_MAX_REPEAT = 12;

export const BVS_PARAMETERS = 'bvparm2013';

// ---------------------------------------------------------------------------
// 1. Lookup
// ---------------------------------------------------------------------------

export function bvsStatesOf(element: string): number[] | null {
    const t = BVS_CATIONS[element];
    if (!t) return null;
    const keys = Object.keys(t);
    const out: number[] = [];
    for (let i = 0; i < keys.length; i++) out.push(parseInt(keys[i], 10));
    return out;
}

export function bvsAnionValence(element: string): number | null {
    return Object.prototype.hasOwnProperty.call(BVS_ANIONS, element) ? BVS_ANIONS[element] : null;
}

/**
 * What a site of this element can be checked as. A cation is known by the
 * oxidation states the table has parameters for; an anion by its charge. An
 * element that is both in this table's scope lists as a cation (none does).
 */
export function bvsClassify(element: string): { kind: string; states: number[] | null; valence: number | null } {
    const states = bvsStatesOf(element);
    if (states) return { kind: 'cation', states: states, valence: null };
    const v = bvsAnionValence(element);
    if (v !== null) return { kind: 'anion', states: null, valence: v };
    return { kind: 'none', states: null, valence: null };
}

export function bvsPair(element: string, state: number, anion: string): [number, number] | null {
    const t = BVS_CATIONS[element];
    if (!t) return null;
    const s = t[String(state)];
    if (!s) return null;
    const p = s[anion];
    return p || null;
}

export function bvsPairCutoff(pair: [number, number]): number {
    const c = pair[0] + BVS_TAIL * pair[1];
    return c < BVS_MAX_CUT ? c : BVS_MAX_CUT;
}

// ---------------------------------------------------------------------------
// 2. Lattice geometry
// ---------------------------------------------------------------------------

/**
 * Fractional-from-cartesian transform for the viewer's lattice convention:
 * the lattice vectors are rows, so a Cartesian position is f times the matrix
 * and the transform back is the transpose of the matrix inverse. The client
 * uses this same arrangement, which is what keeps the module's fractional
 * coordinates identical to the ones the page computes elsewhere.
 */
function bvsInvert(lv: number[][]): number[][] | null {
    const det = lv[0][0] * (lv[1][1] * lv[2][2] - lv[1][2] * lv[2][1]) -
        lv[0][1] * (lv[1][0] * lv[2][2] - lv[1][2] * lv[2][0]) +
        lv[0][2] * (lv[1][0] * lv[2][1] - lv[1][1] * lv[2][0]);
    if (Math.abs(det) < 1e-12) return null;
    const s = 1 / det;
    return [
        [(lv[1][1] * lv[2][2] - lv[1][2] * lv[2][1]) * s, (lv[1][2] * lv[2][0] - lv[1][0] * lv[2][2]) * s, (lv[1][0] * lv[2][1] - lv[1][1] * lv[2][0]) * s],
        [(lv[0][2] * lv[2][1] - lv[0][1] * lv[2][2]) * s, (lv[0][0] * lv[2][2] - lv[0][2] * lv[2][0]) * s, (lv[0][1] * lv[2][0] - lv[0][0] * lv[2][1]) * s],
        [(lv[0][1] * lv[1][2] - lv[0][2] * lv[1][1]) * s, (lv[0][2] * lv[1][0] - lv[0][0] * lv[1][2]) * s, (lv[0][0] * lv[1][1] - lv[0][1] * lv[1][0]) * s]
    ];
}

/** Distance between adjacent lattice planes of each family: how far a repeat moves. */
function bvsPlaneWidths(lv: number[][]): number[] {
    const a = lv[0], b = lv[1], c = lv[2];
    const bcx = b[1] * c[2] - b[2] * c[1], bcy = b[2] * c[0] - b[0] * c[2], bcz = b[0] * c[1] - b[1] * c[0];
    const cax = c[1] * a[2] - c[2] * a[1], cay = c[2] * a[0] - c[0] * a[2], caz = c[0] * a[1] - c[1] * a[0];
    const abx = a[1] * b[2] - a[2] * b[1], aby = a[2] * b[0] - a[0] * b[2], abz = a[0] * b[1] - a[1] * b[0];
    const det = Math.abs(a[0] * bcx + a[1] * bcy + a[2] * bcz);
    return [det / Math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz),
        det / Math.sqrt(cax * cax + cay * cay + caz * caz),
        det / Math.sqrt(abx * abx + aby * aby + abz * abz)];
}

/**
 * A lower bound on the Cartesian distance between an atom and any image of
 * another: the displacement along each reciprocal normal cannot shrink below
 * its wrapped fractional component times that family's plane spacing, so the
 * largest of the three is a bound. Rejecting a pair whose bound already
 * exceeds the cutoff avoids walking the repeats for every distant pair, which
 * in a framework is nearly all of them.
 */
function bvsLowerBound(fx: number, fy: number, fz: number, widths: number[]): number {
    let m = Math.abs(fx - Math.round(fx)) * widths[0];
    const y = Math.abs(fy - Math.round(fy)) * widths[1];
    if (y > m) m = y;
    const z = Math.abs(fz - Math.round(fz)) * widths[2];
    if (z > m) m = z;
    return m;
}

/**
 * Sum s over every image of the displacement within the cutoff. The search
 * window is the same bound the shortest-image search uses: the repeats a
 * displacement of this length can reach across each plane family. A zero
 * distance is the atom's own image of itself, not a contact.
 */
function bvsSumImages(fx: number, fy: number, fz: number, cut: number, lv: number[][],
    widths: number[], r0: number, b: number): { s: number; c: number } {
    const n: number[] = [1, 1, 1];
    for (let k = 0; k < 3; k++) {
        const q = widths[k] > 1e-9 ? cut / widths[k] : 0;
        let v = q > 1 ? Math.ceil(q + 1e-9) : 1;
        if (!(v >= 1)) v = 1;
        if (v > BVS_MAX_REPEAT) v = BVS_MAX_REPEAT;
        n[k] = v;
    }
    const loX = Math.ceil(fx - n[0] - 1e-9), hiX = Math.floor(fx + n[0] + 1e-9);
    const loY = Math.ceil(fy - n[1] - 1e-9), hiY = Math.floor(fy + n[1] + 1e-9);
    const loZ = Math.ceil(fz - n[2] - 1e-9), hiZ = Math.floor(fz + n[2] + 1e-9);
    const cut2 = cut * cut;
    let s = 0, c = 0;
    for (let sx = loX; sx <= hiX; sx++) {
        const gx = fx - sx;
        for (let sy = loY; sy <= hiY; sy++) {
            const gy = fy - sy;
            for (let sz = loZ; sz <= hiZ; sz++) {
                const gz = fz - sz;
                const ax = gx * lv[0][0] + gy * lv[1][0] + gz * lv[2][0];
                const ay = gx * lv[0][1] + gy * lv[1][1] + gz * lv[2][1];
                const az = gx * lv[0][2] + gy * lv[1][2] + gz * lv[2][2];
                const d2 = ax * ax + ay * ay + az * az;
                if (d2 <= cut2 && d2 > 1e-8) {
                    s += Math.exp((r0 - Math.sqrt(d2)) / b);
                    c++;
                }
            }
        }
    }
    return { s: s, c: c };
}

// ---------------------------------------------------------------------------
// 3. The analysis
// ---------------------------------------------------------------------------

/**
 * Bond-valence sums for every site of one cell.
 *
 * Cations are summed first, because an anion's sum depends on which oxidation
 * state each neighbouring cation was read as. A site with several tabulated
 * states keeps the state whose sum is closest to that state: Pb(IV) fits a
 * PbO2 site and Pb(II) fits a PbO site, and nothing else had to be told which
 * one it was looking at. An anion's contacts then use each cation neighbour's
 * chosen state, and are compared with the anion's own charge.
 *
 * The lattice search counts every periodic image inside each pair's cutoff,
 * so the sum is periodic whether or not the cell happens to contain both ends
 * of a contact: one NaCl formula unit still gives sodium six chloride
 * neighbours.
 */
export function bvsAnalyse(atoms: BvsAtom[], lv: number[][]): BvsReport {
    const empty: BvsReport = {
        ok: false, reason: '', parameters: BVS_PARAMETERS, maxCut: BVS_MAX_CUT,
        perAtom: [], warnings: [],
        counts: { atoms: 0, ok: 0, noParams: 0, noContacts: 0, skipped: 0, cations: 0, anions: 0 }
    };
    const n = atoms ? atoms.length : 0;
    empty.counts.atoms = n;
    if (!n) { empty.reason = 'the structure has no atoms'; return empty; }
    if (!lv || lv.length < 3) { empty.reason = 'the structure has no unit cell'; return empty; }
    if (n > BVS_MAX_ATOMS) { empty.reason = 'too many atoms for an on-demand report (' + n + '; the limit is ' + BVS_MAX_ATOMS + ')'; return empty; }
    const inv = bvsInvert(lv);
    if (!inv) { empty.reason = 'the unit cell is degenerate'; return empty; }
    const widths = bvsPlaneWidths(lv);

    const frac: number[][] = new Array(n);
    for (let i = 0; i < n; i++) {
        const a = atoms[i];
        frac[i] = [inv[0][0] * a.x + inv[0][1] * a.y + inv[0][2] * a.z,
            inv[1][0] * a.x + inv[1][1] * a.y + inv[1][2] * a.z,
            inv[2][0] * a.x + inv[2][1] * a.y + inv[2][2] * a.z];
    }

    const kinds: string[] = new Array(n);
    const statesOf: (number[] | null)[] = new Array(n);
    const valenceOf: (number | null)[] = new Array(n);
    const cations: number[] = [];
    const anions: number[] = [];
    for (let i = 0; i < n; i++) {
        const a = atoms[i];
        if (a.ghost) { kinds[i] = 'none'; statesOf[i] = null; valenceOf[i] = null; continue; }
        const cls = bvsClassify(a.element);
        kinds[i] = cls.kind;
        statesOf[i] = cls.states;
        valenceOf[i] = cls.valence;
        if (cls.kind === 'cation') cations.push(i);
        else if (cls.kind === 'anion') anions.push(i);
    }

    const perAtom: (BvsAtomResult | null)[] = new Array(n);
    for (let i = 0; i < n; i++) perAtom[i] = null;

    // Pass 1: every cation site, with its best-fitting oxidation state.
    for (let ci = 0; ci < cations.length; ci++) {
        const i = cations[ci];
        const el = atoms[i].element;
        const states = statesOf[i] as number[];
        const fi = frac[i];
        let hasPair = false;
        let best: { state: number; sum: number; count: number; delta: number } | null = null;
        for (let si = 0; si < states.length; si++) {
            const state = states[si];
            let sum = 0, count = 0;
            for (let aj = 0; aj < anions.length; aj++) {
                const j = anions[aj];
                const pair = bvsPair(el, state, atoms[j].element);
                if (!pair) continue;
                hasPair = true;
                const cut = bvsPairCutoff(pair);
                const fj = frac[j];
                const fx = fj[0] - fi[0], fy = fj[1] - fi[1], fz = fj[2] - fi[2];
                if (bvsLowerBound(fx, fy, fz, widths) > cut) continue;
                const r = bvsSumImages(fx, fy, fz, cut, lv, widths, pair[0], pair[1]);
                sum += r.s;
                count += r.c;
            }
            if (count > 0) {
                const delta = Math.abs(sum - state);
                if (!best || delta < best.delta) best = { state: state, sum: sum, count: count, delta: delta };
            }
        }
        if (best) {
            perAtom[i] = { kind: 'cation', element: el, state: best.state, valence: best.state, sum: best.sum, delta: best.delta, contacts: best.count, status: 'ok' };
        } else {
            perAtom[i] = { kind: 'cation', element: el, state: null, valence: null, sum: null, delta: null, contacts: 0, status: (hasPair || !anions.length) ? 'no-contacts' : 'no-params' };
        }
    }

    // Pass 2: every anion site, against the cations' chosen states.
    for (let ai = 0; ai < anions.length; ai++) {
        const j = anions[ai];
        const el = atoms[j].element;
        const valence = valenceOf[j] as number;
        const fj = frac[j];
        let sum = 0, count = 0, hasPair = false;
        for (let ci = 0; ci < cations.length; ci++) {
            const i = cations[ci];
            const st = perAtom[i];
            if (!st || st.status !== 'ok') continue;
            const pair = bvsPair(atoms[i].element, st.state as number, el);
            if (!pair) continue;
            hasPair = true;
            const cut = bvsPairCutoff(pair);
            const fi = frac[i];
            const fx = fi[0] - fj[0], fy = fi[1] - fj[1], fz = fi[2] - fj[2];
            if (bvsLowerBound(fx, fy, fz, widths) > cut) continue;
            const r = bvsSumImages(fx, fy, fz, cut, lv, widths, pair[0], pair[1]);
            sum += r.s;
            count += r.c;
        }
        if (count === 0 && !hasPair && cations.length) {
            // A pair can exist even when no cation reached it: every cation was
            // itself without contacts, so none has a chosen state. The anion's
            // answer is still "nothing in range", not "no parameters".
            for (let ci = 0; ci < cations.length && !hasPair; ci++) {
                const sts = statesOf[cations[ci]] as number[];
                for (let si = 0; si < sts.length; si++) {
                    if (bvsPair(atoms[cations[ci]].element, sts[si], el)) { hasPair = true; break; }
                }
            }
        }
        if (count > 0) {
            perAtom[j] = { kind: 'anion', element: el, state: null, valence: valence, sum: sum, delta: Math.abs(sum - valence), contacts: count, status: 'ok' };
        } else {
            perAtom[j] = { kind: 'anion', element: el, state: null, valence: valence, sum: null, delta: null, contacts: 0, status: (hasPair || !cations.length) ? 'no-contacts' : 'no-params' };
        }
    }

    // Counts and warnings, aggregated by element so a framework's forty
    // unparameterised oxygens say one thing rather than forty.
    const counts: BvsCounts = { atoms: n, ok: 0, noParams: 0, noContacts: 0, skipped: 0, cations: cations.length, anions: anions.length };
    const skipByEl: { [el: string]: number } = {};
    const noContactByEl: { [el: string]: number } = {};
    const noParamByEl: { [el: string]: number } = {};
    for (let i = 0; i < n; i++) {
        const r = perAtom[i];
        if (!r) {
            if (!atoms[i].ghost) {
                counts.skipped++;
                skipByEl[atoms[i].element] = (skipByEl[atoms[i].element] || 0) + 1;
            }
            continue;
        }
        if (r.status === 'ok') counts.ok++;
        else if (r.status === 'no-contacts') { counts.noContacts++; noContactByEl[r.element] = (noContactByEl[r.element] || 0) + 1; }
        else { counts.noParams++; noParamByEl[r.element] = (noParamByEl[r.element] || 0) + 1; }
    }

    const describe = (m: { [el: string]: number }, total: number, tail: string): string => {
        const els = Object.keys(m).sort((a, b) => m[b] - m[a]);
        let shown = '';
        for (let i = 0; i < els.length && i < 6; i++) {
            if (shown !== '') shown += ', ';
            shown += els[i] + ' (' + m[els[i]] + ')';
        }
        if (els.length > 6) shown += ', and ' + (els.length - 6) + ' more elements';
        return total + ' site(s) ' + tail + ': ' + shown + '.';
    };
    const warnings: string[] = [];
    if (counts.skipped) warnings.push(describe(skipByEl, counts.skipped, 'skipped for want of bond-valence parameters'));
    if (counts.noParams) warnings.push(describe(noParamByEl, counts.noParams, 'have no usable parameter pair'));
    if (counts.noContacts) warnings.push(describe(noContactByEl, counts.noContacts, 'have parameters but no contacts inside the cutoffs'));

    const report: BvsReport = {
        ok: true, reason: '', parameters: BVS_PARAMETERS, maxCut: BVS_MAX_CUT,
        perAtom: perAtom, warnings: warnings, counts: counts
    };
    return report;
}

// ---------------------------------------------------------------------------
// 4. Serialisation for the webview
// ---------------------------------------------------------------------------

const BVS_MODULE_FUNCTIONS: Function[] = [
    bvsStatesOf, bvsAnionValence, bvsClassify, bvsPair, bvsPairCutoff,
    bvsInvert, bvsPlaneWidths, bvsLowerBound, bvsSumImages, bvsAnalyse
];

/**
 * The module as one self-contained expression, for the webview to evaluate.
 * The table and the constants are emitted as JSON and the functions as their
 * own source text inside one scope, so the names they reference resolve to
 * the declarations above them. Both sides run this text.
 */
export const BVS_SOURCE: string = (function buildSource(): string {
    const parts: string[] = [];
    parts.push('(function(){');
    parts.push('var BVS_CATIONS=' + JSON.stringify(BVS_CATIONS) + ';');
    parts.push('var BVS_ANIONS=' + JSON.stringify(BVS_ANIONS) + ';');
    parts.push('var BVS_MAX_CUT=' + BVS_MAX_CUT + ';');
    parts.push('var BVS_TAIL=' + BVS_TAIL + ';');
    parts.push('var BVS_MAX_ATOMS=' + BVS_MAX_ATOMS + ';');
    parts.push('var BVS_MAX_REPEAT=' + BVS_MAX_REPEAT + ';');
    parts.push('var BVS_PARAMETERS=' + JSON.stringify(BVS_PARAMETERS) + ';');
    for (let i = 0; i < BVS_MODULE_FUNCTIONS.length; i++) {
        // The compiled module writes its exported constants as `exports.X`.
        // The page has no `exports`; inside this scope every such name is one
        // of the local variables emitted above, so qualifying it is removed
        // and the function body binds to that local exactly as the host's
        // body binds to its module variable. A reference to any export not
        // emitted above would become a free variable in the page, so it is
        // refused here, at module load, rather than when a panel opens.
        const body = BVS_MODULE_FUNCTIONS[i].toString();
        const refs = body.match(/exports\.\w+/g) || [];
        for (let r = 0; r < refs.length; r++) {
            if (refs[r] !== 'exports.BVS_MAX_CUT' && refs[r] !== 'exports.BVS_MAX_ATOMS' && refs[r] !== 'exports.BVS_PARAMETERS') {
                throw new Error('bond-valence source refers to an export the page does not define: ' + refs[r]);
            }
        }
        parts.push(body.split('exports.').join(''));
    }
    parts.push('return {analyse:bvsAnalyse, classify:bvsClassify, pair:bvsPair, cutoff:bvsPairCutoff, states:bvsStatesOf, anionValence:bvsAnionValence, parameters:BVS_PARAMETERS, maxCut:BVS_MAX_CUT};');
    parts.push('})()');
    const source = parts.join('\n');
    if (source.indexOf('`') >= 0 || source.indexOf('${') >= 0) {
        throw new Error('bond-valence source contains a template-literal delimiter');
    }
    return source;
})();
