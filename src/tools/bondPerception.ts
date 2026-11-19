/**
 * Bond perception from bare 3D coordinates: connectivity, bond orders, and
 * formal charges.
 *
 * One pass decides the connectivity; two checks then decide the orders, and
 * both are needed because each covers what the other cannot.
 *
 * The geometric check reads each bond length against a per-element-pair table
 * of (order, length, tolerance) entries and takes the entry the measured
 * distance actually supports. It is fast, and it is the only one of the two
 * that can express an aromatic bond (order 1.5). Its blind spot is that it
 * judges one bond at a time: a centre whose neighbours add up to more than its
 * valence ceiling (an oxygen that also accepts a hydrogen-bond contact, a
 * nitrogen in a nitro group) cannot be repaired by re-reading lengths, because
 * every candidate it is allowed to consider is already the one the geometry
 * asks for.
 *
 * The valence check works on the graph instead. Every atom is given a target
 * valence drawn from the element's allowed set; unsaturated atoms are paired
 * up over the connectivity graph by a maximum-cardinality matching, each
 * matched pair gaining one unit of bond order, until no unsaturated pair
 * remains; a candidate assignment is accepted only when every atom's order sum
 * equals its target valence exactly, and — when the total charge of the
 * molecule is known — when the formal charges derived from those valences sum
 * to it. It never looks at a bond length, so it can settle a structure the
 * length table cannot, and it is the only path that yields formal charges.
 *
 * Their order is fixed, and that order is what keeps the result monotonic. The
 * geometric check runs first and its result is kept unless some atom in the
 * structure ends up over its ceiling, in which case the valence check
 * re-derives that connected component and takes over only if it finds an
 * assignment that passes every check — and only when the answer leaves the
 * structure no worse on the measure the check exists to improve. A component
 * the valence check cannot settle keeps the geometric result, which the local
 * relaxation has already brought as close to sane as lengths allow. So a
 * structure whose geometry is self-consistent is bit-identical to what the
 * viewer drew before, and a structure whose geometry is not gets a chemically
 * consistent answer instead of a silently wrong one.
 *
 * The whole module is written so that a single copy of it serves both the
 * extension host and the webview: every function is self-contained (it
 * references only its own locals, the other functions here, and the tables
 * declared at the top of this file), and `BOND_PERCEPTION_SOURCE` serialises
 * the tables and the functions into one expression the webview can evaluate.
 * The client therefore runs the same code the host does, not a hand-kept copy
 * of it.
 */

import { Atom, Bond } from '../types';

// ---------------------------------------------------------------------------
// 0. Element data
// ---------------------------------------------------------------------------

/** Covalent radii in Å, used only where no per-pair cutoff is tabulated. */
const COV_R: { [key: string]: number } = {
    H: 0.31, He: 0.28, Li: 1.28, Be: 0.96, B: 0.85, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Ne: 0.58,
    Na: 1.66, Mg: 1.41, Al: 1.21, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02, Ar: 1.06,
    K: 2.03, Ca: 1.76, Sc: 1.70, Ti: 1.60, V: 1.53, Cr: 1.39, Mn: 1.39, Fe: 1.32,
    Co: 1.26, Ni: 1.24, Cu: 1.32, Zn: 1.22, Ga: 1.22, Ge: 1.20, As: 1.19, Se: 1.20,
    Br: 1.20, Kr: 1.16, I: 1.39,
    Rb: 2.20, Sr: 1.95, Y: 1.90, Zr: 1.75, Nb: 1.64, Mo: 1.54, Tc: 1.47, Ru: 1.46,
    Rh: 1.42, Pd: 1.39, Ag: 1.45, Cd: 1.44, In: 1.42, Sn: 1.39, Sb: 1.39, Te: 1.38,
    Xe: 1.40,
    Cs: 2.44, Ba: 2.15, La: 2.07, Ce: 2.04, Pr: 2.03, Nd: 2.01, Pm: 1.99, Sm: 1.98,
    Eu: 1.98, Gd: 1.96, Tb: 1.94, Dy: 1.92, Ho: 1.92, Er: 1.89, Tm: 1.90, Yb: 1.87,
    Lu: 1.87, Hf: 1.75, Ta: 1.70, W: 1.62, Re: 1.51, Os: 1.44, Ir: 1.41, Pt: 1.36,
    Au: 1.36, Hg: 1.32, Tl: 1.45, Pb: 1.46, Bi: 1.48, Po: 1.40, At: 1.50, Rn: 1.50
};

/** Distance beyond which an element pair can never bond, when tabulated. */
const PAIR_CUT: { [key: string]: number } = {
    // H-H is the shortest real bond there is — 0.74 Å in H2, 0.9 Å in a
    // dihydrogen complex — and this entry used to be 0.0, which cut the pair off
    // before any distance could qualify: a hydrogen molecule, a bridging
    // hydride and a dihydrogen ligand all came out as two separate atoms, with
    // no bond, no valence and no force-field term. 1.1 Å covers every bonded
    // H···H and stops well short of the shortest hydrogen-bond contact (1.4 Å
    // in a very strong one, 1.6–2.0 Å in the ordinary case).
    HH: 1.1, CH: 1.3, HO: 1.2, HN: 1.3,
    CC: 1.9, CO: 1.7, CN: 1.7, NN: 1.7,
    NO: 1.8, CF: 1.6, CS: 2.0
};

interface OrderSpec { o: number; l: number; t: number; }

/**
 * Measured length each bond order implies for an element pair. This is the
 * only place bond lengths turn into chemistry, and the only table that can
 * produce the aromatic order 1.5.
 *
 * The aromatic entry is listed before the double bond of the same pair, and
 * that order is not what decides an alkene from an aromatic bond: both bands
 * bracket a length of 1.38–1.40, so the length alone cannot tell them apart.
 * The ring is what tells them apart, and `bpAromaticRing` is where that is
 * decided — a 1.5 that does not lie on a cycle of 1.5 bonds is downgraded to
 * the best order that is not aromatic. Listing the aromatic entry first
 * therefore only means "an aromatic ring is the reading a length in the overlap
 * supports", not "every bond of that length is aromatic".
 */
const SPECS: { [key: string]: OrderSpec[] } = {
    'C+C': [
        { o: 3.0, l: 1.20, t: 0.05 },
        { o: 1.5, l: 1.39, t: 0.05 },
        { o: 2.0, l: 1.38, t: 0.05 },
        { o: 1.0, l: 1.51, t: 0.10 }
    ],
    'C+N': [
        { o: 2.0, l: 1.26, t: 0.05 },
        // The aromatic window runs to 1.42 A on purpose: a pyrrole's or an
        // indole's C–N is 1.37 A, an N-aryl bond is 1.40 A, and a C–N of 1.42 A
        // is still part of a ring's delocalised system — while a 1.43 A one is
        // an amine's single bond. A narrower window dropped such a ring onto the
        // single-bond band and turned the whole ring into a Kekulé structure.
        { o: 1.5, l: 1.36, t: 0.06 },
        { o: 1.0, l: 1.43, t: 0.10 },
        { o: 3.0, l: 1.16, t: 0.06 }
    ],
    'C+O': [
        { o: 2.0, l: 1.24, t: 0.05 },
        { o: 1.0, l: 1.39, t: 0.05 }
    ],
    'N+N': [
        { o: 1.0, l: 1.41, t: 0.10 },
        { o: 2.0, l: 1.25, t: 0.06 },
        { o: 3.0, l: 1.10, t: 0.06 }
    ],
    'N+O': [
        { o: 2.0, l: 1.20, t: 0.06 },
        { o: 1.5, l: 1.30, t: 0.06 },
        { o: 1.0, l: 1.40, t: 0.15 }
    ],
    'O+O': [
        { o: 2.0, l: 1.21, t: 0.06 },
        { o: 1.0, l: 1.48, t: 0.15 }
    ],
    'C+S': [
        { o: 1.5, l: 1.73, t: 0.06 },
        { o: 2.0, l: 1.60, t: 0.10 },
        { o: 1.0, l: 1.82, t: 0.15 }
    ],
    // Sulphur–oxygen is the pair the second-row octet cannot describe, and the
    // length is the only evidence the file gives: a sulphoxide's S=O (1.48 Å),
    // a sulphone's two (1.44 Å) and a sulphate's four (1.47 Å) all read as
    // double bonds. What separates a sulphate from a sulphone is not a length
    // but the count — four double bonds put sulphur at eight, over its ceiling
    // of six, and the relaxation then brings it back to the two-double/two-single
    // pair of valences six; two double bonds already sum to six and are left
    // alone. The single-bond band is deliberately wide (1.42–1.82 Å) so the
    // reduction has somewhere to go: a sulphate's 1.47 Å S–O is 0.17 Å from the
    // single ideal and only a tolerance of 0.20 brackets it. Without that, the
    // ladder of a sulphate's S–O bond held one entry (the double), a centre two
    // units over its ceiling had nothing to reach for, and sulphur came out at
    // eight with a +6 formal charge.
    'O+S': [
        { o: 2.0, l: 1.46, t: 0.08 },
        { o: 1.0, l: 1.62, t: 0.20 }
    ],
    // Phosphorus–oxygen: phosphine oxide P=O is 1.48 Å, a phosphate's four P–O
    // are 1.51–1.55 Å, and the double band stops at 1.51 so the anion keeps its
    // single bonds (and the dative charge that balances its four oxygens).
    'O+P': [
        { o: 2.0, l: 1.47, t: 0.04 },
        { o: 1.0, l: 1.59, t: 0.10 }
    ],
    // A silicon–oxygen bond is a single bond at every length it is found at
    // (1.58–1.68 Å, siloxane and silicate alike). The entry exists so the pair
    // has a length at all: without one, `bpSpecDev` reported a perfect fit for
    // any order the relaxation proposed and the pass was free to invent a
    // one-and-a-half bond out of a distance that never argued for one.
    'O+Si': [{ o: 1.0, l: 1.63, t: 0.10 }],
    // Chlorine–oxygen is drawn single on purpose: perchlorate and chlorate are
    // the dative reading, in which Cl carries +3 and each oxygen −1, and that
    // is the only reading whose charges add up to the anion's own charge.
    'Cl+O': [{ o: 1.0, l: 1.70, t: 0.14 }],
    'C+F': [{ o: 1.0, l: 1.33, t: 0.10 }],
    'C+H': [{ o: 1.0, l: 0.97, t: 0.15 }],
    'N+H': [{ o: 1.0, l: 0.88, t: 0.15 }],
    'O+H': [{ o: 1.0, l: 0.85, t: 0.15 }]
};

/** Hard valence ceiling per element. Never exceeded by a sane assignment. */
const CEIL: { [key: string]: number } = {
    // The halogens are not all monovalent: perchlorate, chlorate, the chlorine
    // fluorides and their bromine and iodine analogues are ordinary chemistry,
    // and a ceiling of 1 made every one of them permanently over-coordinated.
    // The relaxation then had a centre it could never satisfy — it spent its
    // budget on it in every pass, and the diagnostic listed it for ever. A
    // ceiling of 7 is the most a halogen is ever seen to carry; the bonds it
    // actually has are still decided by the lengths.
    //
    // Nitrogen is four, not three: an ammonium, a pyridinium and a nitro group
    // all carry four bonds, and the answer for them is a positive formal charge
    // (which `bpAtomicCharge` states), not a repair. At a ceiling of three every
    // one of them counted as over-coordinated, and the relaxation then spent its
    // budget lowering the bonds of a perfectly good aromatic ring to bring the
    // nitrogen down — a pyrrole's or a pyridinium's two aromatic bonds came out
    // single and the ring stopped being aromatic, from a pass whose whole job is
    // to remove excess that was never there. A nitrogen measured with five bonds
    // is still over, and still repaired.
    H: 1, C: 4, N: 4, O: 2, F: 1, S: 6, P: 5, Cl: 7, Br: 7, I: 7, B: 3
};

/**
 * Elements that carry exactly one bond, so they have no bond order left to
 * reassign: their single sigma bond is their chemistry, and lowering it would
 * draw an empty valence shell instead of a plain single bond.
 */
const TERMINAL: { [key: string]: number } = {
    H: 1, F: 1, Cl: 1, Br: 1, I: 1, He: 1, Ne: 1, Ar: 1, Kr: 1, Xe: 1, Rn: 1
};

/**
 * Atomic number per element symbol, for every element the viewer can draw.
 *
 * The table used to stop at the main group, so every transition metal read as
 * −1, `bpPossibleValences` answered nothing for it, and the valence solver
 * refused the *whole component* an iron atom belonged to — a ferrocene or a
 * metalloprotein kept whatever the geometric pass produced even where the
 * organic part needed repairing. The numbers are here now; what a metal has no
 * entry for is a valence (see `bpValenceModelled`), which is a different
 * question and is answered by leaving such an atom frozen rather than by
 * disabling its neighbours.
 */
const Z_OF: { [key: string]: number } = {
    H: 1, He: 2, Li: 3, Be: 4, B: 5, C: 6, N: 7, O: 8, F: 9, Ne: 10,
    Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, Ar: 18,
    K: 19, Ca: 20, Sc: 21, Ti: 22, V: 23, Cr: 24, Mn: 25, Fe: 26, Co: 27,
    Ni: 28, Cu: 29, Zn: 30, Ga: 31, Ge: 32, As: 33, Se: 34, Br: 35, Kr: 36,
    Rb: 37, Sr: 38, Y: 39, Zr: 40, Nb: 41, Mo: 42, Tc: 43, Ru: 44, Rh: 45,
    Pd: 46, Ag: 47, Cd: 48, In: 49, Sn: 50, Sb: 51, Te: 52, I: 53, Xe: 54,
    Cs: 55, Ba: 56, La: 57, Ce: 58, Pr: 59, Nd: 60, Pm: 61, Sm: 62, Eu: 63,
    Gd: 64, Tb: 65, Dy: 66, Ho: 67, Er: 68, Tm: 69, Yb: 70, Lu: 71, Hf: 72,
    Ta: 73, W: 74, Re: 75, Os: 76, Ir: 77, Pt: 78, Au: 79, Hg: 80,
    Tl: 81, Pb: 82, Bi: 83, Po: 84, At: 85, Rn: 86
};

/**
 * Elements whose bonds have no length table and no valence model: the d and f
 * block.
 *
 * Their bond *lengths* are still read from the covalent radii (see
 * `bpMetalSpec`) — a metal–oxygen bond at 2.05 Å is a bond, and its order is
 * one — but their *valence* is not a number a file can be asked about: an iron
 * is six-coordinate in one complex and four in the next, and which of its
 * ligands are "bonds" is the chemistry of the compound, not of the distance.
 * The solver therefore treats a metal as frozen (see `bpDetermineOrders`) and
 * the formal charge as unknown, which is what the module promises for an
 * element it has no entry for.
 */
const METALS: { [key: string]: number } = {
    Sc: 1, Ti: 1, V: 1, Cr: 1, Mn: 1, Fe: 1, Co: 1, Ni: 1, Cu: 1, Zn: 1,
    Y: 1, Zr: 1, Nb: 1, Mo: 1, Tc: 1, Ru: 1, Rh: 1, Pd: 1, Ag: 1, Cd: 1,
    La: 1, Ce: 1, Pr: 1, Nd: 1, Pm: 1, Sm: 1, Eu: 1, Gd: 1, Tb: 1, Dy: 1,
    Ho: 1, Er: 1, Tm: 1, Yb: 1, Lu: 1, Hf: 1, Ta: 1, W: 1, Re: 1, Os: 1,
    Ir: 1, Pt: 1, Au: 1, Hg: 1, Ac: 1, Th: 1, Pa: 1, U: 1, Np: 1, Pu: 1
};

/** Target valences an atom is allowed to take, keyed by atomic number. */
const VALENCE: { [key: number]: number[] } = {
    1: [1],
    5: [3, 4],
    6: [4],
    7: [3, 4],
    8: [2, 1, 3],
    9: [1],
    14: [4],
    15: [5, 3],
    // 4 is the sulphoxide: one S=O and two S–C add up to four, which is neither
    // the sulphone's six nor a thioether's two, and leaving it out of the set
    // made every sulphoxide read as off its valence in the diagnostic.
    16: [6, 4, 3, 2, 1],
    // The halogens the oxoacids and the interhalogens are built from: 1 in a
    // halide, 3 in ClF3, 4 in the dative perchlorate, 5 in ClF5, 7 in IF7. The
    // set is what the solver may aim at and what the diagnostic measures
    // against; without these values every hypervalent halogen was reported off
    // its valence by an amount the pass could never repair.
    17: [1, 3, 4, 5, 7],
    32: [4],
    35: [1, 3, 4, 5, 7],
    53: [1, 3, 4, 5, 7]
};

/** Target valences for elements outside VALENCE, by periodic-table group. */
const FALLBACK_VALENCE: { [key: number]: number[] } = {
    3: [1], 4: [2], 11: [1], 12: [2], 13: [3], 19: [1], 20: [2],
    31: [3], 33: [3, 5], 34: [2, 6], 37: [1], 38: [2], 49: [3], 50: [4],
    51: [3, 5], 52: [2, 6], 55: [1], 56: [2], 81: [3], 82: [4], 83: [3, 5]
};

/** Outermost electron count, the input of the generic formal-charge formula. */
const OUTER_E: { [key: number]: number } = {
    1: 1, 2: 2, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8,
    11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8, 19: 1, 20: 2,
    31: 3, 32: 4, 33: 5, 34: 6, 35: 7, 36: 8, 37: 1, 38: 2, 49: 3, 50: 4,
    51: 5, 52: 6, 53: 7, 54: 8, 55: 1, 56: 2, 81: 3, 82: 4, 83: 5, 84: 6,
    85: 7, 86: 8
};

/** Structures above this many atoms skip the local relaxation sweep. */
const MAX_RELAX_ATOMS = 20000;

/** Components above this many atoms keep the geometric result outright. */
const MAX_SOLVER_ATOMS = 3000;

/** Valence combinations a component may enumerate before it gives up. */
const MAX_COMBOS = 200000;

/**
 * Element operations one component may spend on those combinations.
 *
 * A bound on the count alone is not a bound on the work: the sweep copies a
 * `k`×`k` order matrix and runs a blossom matching inside every combination it
 * looks at, so 200 000 combinations of a 3 000-atom component is 200 000 × 10⁷
 * operations — minutes to hours on one core, from a single structure that
 * merely happens to contain a component the solver was offered. The budget is
 * charged in the units the work is actually made of (`k` for a combination's
 * per-atom checks, `k`×`k` for the matrix, `k`×`k` per matching round), so a
 * small molecule still walks the whole combination space — a 30-atom component
 * costs about 900 per combination and 200 000 of them fit — while a component
 * large enough for the copy to be expensive is cut short and keeps the
 * geometric answer, which is exactly what `declined` reports.
 */
const MAX_COMBO_WORK = 40000000;

/**
 * Neighbours an over-ceiling centre may lower at once.
 *
 * Three covers every polyvalent centre the tables describe (a sulphate's
 * sulphur needs two, a chlorate's one) and bounds the subset search at a few
 * hundred plans even for a ten-coordinate metal.
 */
const MAX_SHED_BONDS = 3;

/**
 * Distance below which two atoms can never be bonded: shorter than the
 * shortest real bond (H-H, 0.74 Å), so a pair this close is a coordinate
 * defect, not chemistry.
 */
const MIN_BOND_DISTANCE = 0.4;

// ---------------------------------------------------------------------------
// 1. Element helpers
// ---------------------------------------------------------------------------

function bpNorm(el: string): string {
    // A missing or empty element symbol is not something the parsers produce,
    // but every table lookup downstream is keyed by this string, so one bad
    // record would otherwise throw from the middle of a whole-structure pass.
    // An empty symbol instead degrades the way an unknown element does: no
    // table entry, so the generic fallbacks apply.
    if (typeof el !== 'string' || el.length === 0) return '';
    return el.charAt(0).toUpperCase() + el.slice(1).toLowerCase();
}

function bpPairKey(e1: string, e2: string): string {
    const a = e1.toUpperCase();
    const b = e2.toUpperCase();
    return a < b ? a + b : b + a;
}

function bpSpecKey(e1: string, e2: string): string {
    const a = bpNorm(e1);
    const b = bpNorm(e2);
    return a < b ? a + '+' + b : b + '+' + a;
}

function bpDist(a: Atom, b: Atom): number {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    const dz = a.z - b.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function bpIsTerminal(el: string): boolean {
    return TERMINAL[bpNorm(el)] === 1;
}

function bpCeiling(el: string): number {
    const c = CEIL[bpNorm(el)];
    return c === undefined ? 100 : c;
}

function bpZOf(el: string): number {
    const z = Z_OF[bpNorm(el)];
    return z === undefined ? -1 : z;
}

/** Whether the element is a transition or inner-transition metal. */
function bpIsMetal(el: string): boolean {
    return METALS[bpNorm(el)] === 1;
}

/**
 * Whether a valence model exists for an element at all.
 *
 * This is the question the solver asks, and it is not the same as "does the
 * element have a Z". A metal has one and no valence — which is why it is
 * frozen inside the solve rather than used as a reason to refuse the whole
 * component. An element with neither (a dummy, a superheavy placeholder) is
 * frozen the same way.
 */
function bpValenceModelled(z: number): boolean {
    return bpPossibleValences(z) !== null;
}

/**
 * The length table for an element pair, with the metals filled in.
 *
 * A metal has no tabulated lengths — its bonds are 2.0–2.5 Å depending on the
 * ligand, the oxidation state and the coordination number — so the pair would
 * otherwise read as "no evidence", which the relaxation took for a perfect fit
 * and the refinement took for freedom to move the order by half a unit. A
 * single bond at the covalent-radius sum, with the widest tolerance that still
 * is not a second bond, states the one thing that is true of every metal
 * contact: it is a single bond, and a length far from that sum is worth
 * nothing. Nothing above a single bond is ever proposed, which is what keeps a
 * ferrocene's five Fe–C contacts single rather than double.
 */
function bpMetalSpec(el1: string, el2: string): OrderSpec[] | null {
    if (!bpIsMetal(el1) && !bpIsMetal(el2)) return null;
    const r1 = COV_R[bpNorm(el1)];
    const r2 = COV_R[bpNorm(el2)];
    if (r1 === undefined || r2 === undefined) return null;
    return [{ o: 1.0, l: r1 + r2, t: 0.45 }];
}

/** The specs of a pair: the table's own, else the metal fallback, else none. */
function bpSpecList(el1: string, el2: string): OrderSpec[] | null {
    const specs = SPECS[bpSpecKey(el1, el2)];
    if (specs) return specs;
    return bpMetalSpec(el1, el2);
}

/**
 * Whether an atom is a ghost (a position with no nucleus, electrons or bonds).
 *
 * Bond perception must answer this before anything else: a ghost that was
 * allowed into the connectivity pass would be a centre of unknown radius
 * (the 1.5 Å fallback), so it would bond to whatever sits within ~3.5 Å of it,
 * take a share of the valence solver's attention and — in a counterpoise
 * structure, where the ghosts stand exactly on the other fragment's atoms —
 * draw a full set of sticks into points that have no atom at all.
 *
 * Both spellings mean the same thing: the `ghost` flag, and the dummy element
 * `X` every parser normalises a dummy to. Kept identical to `isGhostAtom` in
 * `types.ts`; this copy exists because the module is serialised into the page
 * as source, where nothing can be imported.
 */
function bpIsGhost(a: Atom | null | undefined): boolean {
    if (!a) return true;
    return a.ghost === true || a.element === 'X';
}

// ---------------------------------------------------------------------------
// 2. Connectivity and geometric bond orders
// ---------------------------------------------------------------------------

/**
 * The best order a pair's length supports that is *not* the aromatic one.
 *
 * Used only to downgrade a bond the ring test has ruled out. A tolerance band
 * that brackets the distance wins over the nearest tabulated length, because
 * the two disagree exactly where it matters: a C–C of 1.44 Å is inside the
 * single-bond band (1.41–1.61) and outside the double-bond one (1.33–1.43), yet
 * it is nearer to 1.38 than to 1.51 — so a nearest-fit rule turned the one
 * single bond of a furan's Kekulé form into a second double bond and left the
 * carbon between them five-valent. Where no band brackets the distance the
 * nearest length is the only evidence there is, and it is used.
 */
function bpNonAromaticOrder(el1: string, el2: string, d: number): number {
    const specs = bpSpecList(el1, el2);
    if (!specs) return 1.0;
    let bandOrder = 0;
    let bandDev = Infinity;
    let nearest = 1.0;
    let nearestDiff = Infinity;
    for (let i = 0; i < specs.length; i++) {
        if (Math.abs(specs[i].o - 1.5) < 1e-9) continue;
        const diff = Math.abs(d - specs[i].l);
        if (diff <= specs[i].t && diff < bandDev) {
            bandDev = diff;
            bandOrder = specs[i].o;
        }
        if (diff < nearestDiff) {
            nearestDiff = diff;
            nearest = specs[i].o;
        }
    }
    return bandOrder > 0 ? bandOrder : nearest;
}

/**
 * Take the aromatic order 1.5 off every bond that is not part of an aromatic
 * ring.
 *
 * The C–C length table brackets 1.34–1.44 Å with *both* the aromatic and the
 * double order, so the length alone cannot tell an aromatic bond from an
 * alkene: a four-substituted C=C at 1.40 Å, a ring alkene, a ring imine and a
 * conjugated diene all land in the aromatic band. What the distance cannot
 * state is the ring, and the reading a length in the overlap supports is only
 * "could be aromatic" — the graph has to agree.
 *
 * **The rule is a cycle of aromatic bonds, not merely a cycle.** An earlier
 * version of this pass asked only whether a bond lay on *some* cycle (a bridge
 * search over the whole connectivity), which accepted exactly what it was meant
 * to reject: one alkene bond of a cyclohexene ring is on a cycle, so it stayed
 * 1.5, went into the MM parameters as a half order and was exported to MOL2 as
 * "ar". A ring alkene, a ring imine and a bond of a partly saturated fused
 * system are all on a cycle; none of them is aromatic.
 *
 * Running the same bridge search over the *aromatic edges alone* answers the
 * right question: a 1.5 lies on a cycle all of whose edges are 1.5 if and only
 * if it is not a bridge of that subgraph. Every bond of benzene, pyridine,
 * pyrrole, thiophene, naphthalene and a porphyrin's rings survives, because
 * each of them closes a loop of aromatic bonds. A cyclohexene's single 1.5 is a
 * bridge of a one-bond subgraph and is downgraded to the best order its length
 * supports that is not aromatic; so is a furan's pair of C=C, because its C–O
 * bonds cannot carry 1.5 without putting the oxygen at three (the Kekulé form
 * of furan is two C=C and two C–O, which is what comes back); so is a
 * cyclooctatetraene's, whose alternating 1.33/1.46 ring is a polyene.
 *
 * It also removes the flicker a trajectory had. Under the whole-graph test, any
 * edge crossing the cutoff anywhere in the structure — a hydrogen-bond contact
 * appearing between two molecules, a ligand drifting past 3 Å — could close a
 * ring and re-label bonds that had nothing to do with it, so orders jumped
 * between frames near the cutoff. An edge of order 1 can never join a cycle of
 * order-1.5 edges, so a contact appearing or disappearing now changes nothing
 * about which bonds are aromatic; only the ring the bond belongs to breaking
 * does, which is a real change in the structure rather than an artefact of the
 * threshold.
 *
 * Iterative Tarjan with an explicit stack: a 20 000-atom structure must not be
 * able to overflow the call stack, and the pass is O(V+E) over the aromatic
 * edges, which are a small fraction of the bonds.
 */
function bpAromaticRing(bondMap: Map<number, Map<number, number>>, atoms: Atom[]): number {
    const n = atoms.length;

    // The aromatic subgraph, edges numbered in (i, j) insertion order.
    const from: number[] = [];
    const to: number[] = [];
    const adj: number[][] = new Array(n);
    for (let i = 0; i < n; i++) adj[i] = [];
    bondMap.forEach(function (neighbors, i) {
        neighbors.forEach(function (order, j) {
            if (j <= i || Math.abs(order - 1.5) > 1e-9) return;
            const id = from.length;
            from.push(i);
            to.push(j);
            adj[i].push(id);
            adj[j].push(id);
        });
    });
    const m = from.length;
    if (m === 0) return 0;

    const isBridge = new Uint8Array(m);
    const disc = new Int32Array(n).fill(-1);
    const low = new Int32Array(n);
    const stackNode: number[] = [];
    const stackEdge: number[] = [];
    const stackNext: number[] = [];
    let timer = 0;
    for (let root = 0; root < n; root++) {
        if (disc[root] !== -1) continue;
        disc[root] = timer;
        low[root] = timer;
        timer++;
        stackNode.push(root);
        stackEdge.push(-1);
        stackNext.push(0);
        while (stackNode.length > 0) {
            const top = stackNode.length - 1;
            const u = stackNode[top];
            if (stackNext[top] < adj[u].length) {
                const e = adj[u][stackNext[top]];
                stackNext[top]++;
                if (e === stackEdge[top]) continue;
                const v = from[e] === u ? to[e] : from[e];
                if (disc[v] === -1) {
                    disc[v] = timer;
                    low[v] = timer;
                    timer++;
                    stackNode.push(v);
                    stackEdge.push(e);
                    stackNext.push(0);
                } else if (disc[v] < low[u]) {
                    low[u] = disc[v];
                }
            } else {
                stackNode.pop();
                stackNext.pop();
                const parentEdge = stackEdge.pop() as number;
                if (top > 0) {
                    const p = stackNode[top - 1];
                    if (low[u] < low[p]) low[p] = low[u];
                    // u cannot reach above its parent without using the edge it
                    // came in on, so that edge closes no loop of aromatic bonds.
                    if (low[u] > disc[p]) isBridge[parentEdge] = 1;
                }
            }
        }
    }

    let changed = 0;
    for (let e = 0; e < m; e++) {
        if (!isBridge[e]) continue;
        const i = from[e];
        const j = to[e];
        const ni = bondMap.get(i);
        if (!ni || ni.get(j) !== 1.5) continue;
        const order = bpNonAromaticOrder(atoms[i].element, atoms[j].element, bpDist(atoms[i], atoms[j]));
        ni.set(j, order);
        const nj = bondMap.get(j);
        if (nj) nj.set(i, order);
        changed++;
    }
    return changed;
}

/**
 * Distance above which a (normalised) element pair can never bond: the
 * per-pair cutoff when tabulated, otherwise the covalent radii sum plus half
 * an angstrom.
 */
function bpPairCutoff(e1: string, e2: string): number {
    const c = PAIR_CUT[bpPairKey(e1, e2)];
    if (c !== undefined) return c;
    const r1 = COV_R[bpNorm(e1)];
    const r2 = COV_R[bpNorm(e2)];
    return (r1 === undefined ? 1.5 : r1) + (r2 === undefined ? 1.5 : r2) + 0.5;
}

/**
 * Bond order the measured length implies for an element pair. The first spec
 * whose tolerance brackets the distance wins; failing that, the pair's closest
 * fit; failing that (no table entry at all), a covalent-radius ratio. Returns 0
 * when the distance is beyond the pair's cutoff.
 */
function bpGeometricOrder(el1: string, el2: string, d: number): number {
    // A floor as well as a ceiling. The shortest real bond is H-H at 0.74 A,
    // so anything under this is not a bond — it is two atoms in the same place,
    // which is a data defect rather than a structure. Without the floor such a
    // defect is catastrophic rather than merely wrong: every atom bonds to
    // every other one, so a 3000-atom pile-up produces a nine-million-entry
    // bond map and then an all-pairs search per atom in the relaxation. With
    // it, a pile-up simply has no bonds at all and the pass returns at once.
    if (!(d >= MIN_BOND_DISTANCE)) return 0;
    if (d > bpPairCutoff(bpNorm(el1), bpNorm(el2))) return 0;

    const specs = bpSpecList(el1, el2);
    if (specs) {
        for (let i = 0; i < specs.length; i++) {
            if (Math.abs(d - specs[i].l) <= specs[i].t) return specs[i].o;
        }
        let order = 1.0;
        let minDiff = Infinity;
        for (let i = 0; i < specs.length; i++) {
            const diff = Math.abs(d - specs[i].l);
            if (diff < minDiff) {
                minDiff = diff;
                order = specs[i].o;
            }
        }
        return order;
    }

    // No tabulated lengths for this pair means no evidence for anything above
    // a single bond: a covalent-radius ratio cannot tell an ionic Na-Cl pair
    // from a triple bond (it classified Na-Cl as double and Li-F as triple).
    return 1.0;
}

/**
 * The order whose tolerance brackets this length with the smallest deviation,
 * or 0 when the length falls inside no tolerance band of the pair.
 *
 * Unlike bpGeometricOrder this never falls back to a nearest fit, so the answer
 * says whether the length is *inside* a band rather than merely closest to one.
 * That distinction is what lets a caller compare a bond's situation before and
 * after a move and decide whether re-deriving it could change anything.
 */
function bpSupportedOrder(el1: string, el2: string, d: number): number {
    const specs = bpSpecList(el1, el2);
    if (!specs) return 0;
    let order = 0;
    let best = Infinity;
    for (let i = 0; i < specs.length; i++) {
        const dev = Math.abs(d - specs[i].l);
        if (dev <= specs[i].t && dev < best) {
            best = dev;
            order = specs[i].o;
        }
    }
    return order;
}

/**
 * How far an assigned order sits from the length the table expects for it.
 *
 * **−1 means the pair has no length evidence at all**, and that is a different
 * answer from 0. It used to return 0 for an untabulated pair, which says
 * "perfect fit" — so a silicon–oxygen or a sulphur–oxygen bond, neither of
 * which had a table entry, read as fitting its length perfectly at *every*
 * order, and the refinement was free to move such a bond by half a unit
 * (`bpSpecDev(next) <= bpSpecDev(current)` was 0 ≤ 0 for any `next`). That is
 * where the half-integer Si–O and S–O orders came from. The caller now has to
 * say what no evidence means for it, and every one of them reads it as "this
 * bond cannot be moved on length grounds".
 *
 * `Infinity` is the second special answer: the pair is tabulated but the order
 * in question is not one of its values (a 1.5 on a C–H bond, say), so no length
 * could justify it either.
 */
function bpSpecDev(el1: string, el2: string, order: number, d: number): number {
    const specs = bpSpecList(el1, el2);
    if (!specs) return -1;
    let best = Infinity;
    for (let i = 0; i < specs.length; i++) {
        if (Math.abs(specs[i].o - order) > 1e-9) continue;
        const dev = Math.abs(d - specs[i].l);
        if (dev < best) best = dev;
    }
    return best;
}

/**
 * Every bond the coordinates imply, as a symmetric index map of bond orders.
 *
 * A spatial hash grid (3 Å cells) replaces the all-pairs loop: a bond needs
 * d <= cutoff(element pair), so per element the largest cutoff against every
 * element present gives a search radius that safely covers every pair that
 * could bond. Candidate pairs come out in (i asc, j asc) order, so the
 * insertion sequence — and therefore everything downstream — is the same the
 * all-pairs loop would produce.
 */
function bpConnectivity(atoms: Atom[]): Map<number, Map<number, number>> {
    const n = atoms.length;
    const bondMap = new Map<number, Map<number, number>>();
    if (n === 0) return bondMap;

    const CELL = 3;
    const normEls: string[] = new Array(n);
    const distinct: string[] = [];
    const seenEl: { [key: string]: number } = {};
    for (let i = 0; i < n; i++) {
        const el = bpNorm(atoms[i].element);
        normEls[i] = el;
        if (!seenEl[el]) {
            seenEl[el] = 1;
            distinct.push(el);
        }
    }

    const rangeFor: { [key: string]: number } = {};
    for (let i = 0; i < distinct.length; i++) {
        let maxCut = 0;
        for (let j = 0; j < distinct.length; j++) {
            const c = bpPairCutoff(distinct[i], distinct[j]);
            if (c > maxCut) maxCut = c;
        }
        rangeFor[distinct[i]] = Math.floor(maxCut / CELL) + 1;
    }

    const grid = new Map<number, number[]>();
    const cellX = new Int32Array(n);
    const cellY = new Int32Array(n);
    const cellZ = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        if (bpIsGhost(atoms[i])) continue;
        const cx = Math.floor(atoms[i].x / CELL);
        const cy = Math.floor(atoms[i].y / CELL);
        const cz = Math.floor(atoms[i].z / CELL);
        cellX[i] = cx;
        cellY[i] = cy;
        cellZ[i] = cz;
        // XOR hash collisions only merge buckets (extra candidates are filtered
        // by the distance test), never lose one.
        const key = (cx * 73856093) ^ (cy * 19349663) ^ (cz * 83492791);
        let bucket = grid.get(key);
        if (!bucket) {
            bucket = [];
            grid.set(key, bucket);
        }
        bucket.push(i);
    }

    const candidates: number[] = [];
    for (let i = 0; i < n; i++) {
        if (bpIsGhost(atoms[i])) continue;
        const range = rangeFor[normEls[i]];
        const cx = cellX[i];
        const cy = cellY[i];
        const cz = cellZ[i];
        candidates.length = 0;
        for (let dx = -range; dx <= range; dx++) {
            for (let dy = -range; dy <= range; dy++) {
                for (let dz = -range; dz <= range; dz++) {
                    const key = ((cx + dx) * 73856093) ^ ((cy + dy) * 19349663) ^ ((cz + dz) * 83492791);
                    const bucket = grid.get(key);
                    if (!bucket) continue;
                    for (let k = 0; k < bucket.length; k++) {
                        const j = bucket[k];
                        if (j > i) candidates.push(j);
                    }
                }
            }
        }
        if (candidates.length === 0) continue;
        candidates.sort((a, b) => a - b);
        let prev = -1;
        for (let c = 0; c < candidates.length; c++) {
            const j = candidates[c];
            if (j === prev) continue;
            prev = j;
            const order = bpGeometricOrder(atoms[i].element, atoms[j].element, bpDist(atoms[i], atoms[j]));
            if (order > 0) {
                if (!bondMap.has(i)) bondMap.set(i, new Map<number, number>());
                if (!bondMap.has(j)) bondMap.set(j, new Map<number, number>());
                (bondMap.get(i) as Map<number, number>).set(j, order);
                (bondMap.get(j) as Map<number, number>).set(i, order);
            }
        }
    }
    // Aromaticity is a property of a ring, and the length table cannot tell a
    // ring bond from an alkene of the same length: the order the geometry
    // suggests stands only where the graph closes a loop of aromatic bonds
    // around it.
    if (bondMap.size > 0) bpAromaticRing(bondMap, atoms);
    return bondMap;
}

// ---------------------------------------------------------------------------
// 3. Snapping orders onto the values the table allows
// ---------------------------------------------------------------------------

/** Every order the element pair may carry, read out of the length table. */
function bpAllowedOrders(el1: string, el2: string): number[] {
    const specs = bpSpecList(el1, el2);
    if (!specs) return [1.0];
    const out: number[] = [];
    for (let i = 0; i < specs.length; i++) {
        if (out.indexOf(specs[i].o) < 0) out.push(specs[i].o);
    }
    return out;
}

function bpNearestAllowed(allowed: number[], order: number): number {
    let best = allowed[0];
    let bestDiff = Infinity;
    for (let i = 0; i < allowed.length; i++) {
        const diff = Math.abs(order - allowed[i]);
        if (diff < bestDiff) {
            bestDiff = diff;
            best = allowed[i];
        }
    }
    return best;
}

/**
 * Collapse every order the element pair does not list onto the closest one it
 * does. An order the length table has no entry for is meaningless for that
 * pair, so it cannot survive this pass.
 */
function bpFixOrders(bondMap: Map<number, Map<number, number>>, atoms: Atom[]): void {
    bondMap.forEach(function (neighbors, i) {
        neighbors.forEach(function (order, j) {
            if (i > j) return;
            const allowed = bpAllowedOrders(atoms[i].element, atoms[j].element);
            if (allowed.indexOf(order) >= 0) return;
            const snapped = bpNearestAllowed(allowed, order);
            neighbors.set(j, snapped);
            const nj = bondMap.get(j);
            if (nj) nj.set(i, snapped);
        });
    });
}

// ---------------------------------------------------------------------------
// 4. Maximum-cardinality matching (general graph)
// ---------------------------------------------------------------------------

/**
 * Maximum-cardinality matching of a general graph: mate[v] is v's partner, or
 * -1 when v is unmatched. General graphs may contain odd cycles (five-membered
 * rings, fused systems), which is why this is the blossom algorithm and not a
 * bipartite one.
 */
function bpMaximumMatching(n: number, edges: Array<[number, number]>): number[] {
    const adj: number[][] = [];
    for (let i = 0; i < n; i++) adj.push([]);
    for (let e = 0; e < edges.length; e++) {
        adj[edges[e][0]].push(edges[e][1]);
        adj[edges[e][1]].push(edges[e][0]);
    }

    const match: number[] = new Array(n);
    const p: number[] = new Array(n);
    const base: number[] = new Array(n);
    const used: boolean[] = new Array(n);
    const blossom: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) {
        match[i] = -1;
        p[i] = -1;
        base[i] = i;
        used[i] = false;
        blossom[i] = false;
    }
    let queue: number[] = [];

    function lca(a: number, b: number): number {
        const seen: boolean[] = new Array(n);
        for (let i = 0; i < n; i++) seen[i] = false;
        let x = a;
        for (;;) {
            x = base[x];
            seen[x] = true;
            if (match[x] === -1) break;
            x = p[match[x]];
        }
        let y = b;
        for (;;) {
            y = base[y];
            if (seen[y]) return y;
            y = p[match[y]];
        }
    }

    function markPath(v: number, b: number, child: number): void {
        let x = v;
        let c = child;
        while (base[x] !== b) {
            blossom[base[x]] = true;
            blossom[base[match[x]]] = true;
            p[x] = c;
            c = match[x];
            x = p[match[x]];
        }
    }

    function findPath(root: number): number {
        for (let i = 0; i < n; i++) {
            used[i] = false;
            p[i] = -1;
            base[i] = i;
        }
        queue = [root];
        used[root] = true;
        let head = 0;
        while (head < queue.length) {
            const v = queue[head];
            head++;
            const nbr = adj[v];
            for (let a = 0; a < nbr.length; a++) {
                const to = nbr[a];
                if (base[v] === base[to] || match[v] === to) continue;
                if (to === root || (match[to] !== -1 && p[match[to]] !== -1)) {
                    const curbase = lca(v, to);
                    for (let i = 0; i < n; i++) blossom[i] = false;
                    markPath(v, curbase, to);
                    markPath(to, curbase, v);
                    for (let i = 0; i < n; i++) {
                        if (blossom[base[i]]) {
                            base[i] = curbase;
                            if (!used[i]) {
                                used[i] = true;
                                queue.push(i);
                            }
                        }
                    }
                } else if (p[to] === -1) {
                    p[to] = v;
                    if (match[to] === -1) return to;
                    const next = match[to];
                    used[next] = true;
                    queue.push(next);
                }
            }
        }
        return -1;
    }

    for (let i = 0; i < n; i++) {
        if (match[i] !== -1) continue;
        let v = findPath(i);
        while (v !== -1) {
            const pv = p[v];
            const ppv = match[pv];
            match[v] = pv;
            match[pv] = v;
            v = ppv;
        }
    }
    return match;
}

// ---------------------------------------------------------------------------
// 5. Valence-conservation solver
// ---------------------------------------------------------------------------

function bpPossibleValences(z: number): number[] | null {
    if (VALENCE[z]) return VALENCE[z].slice();
    if (FALLBACK_VALENCE[z]) return FALLBACK_VALENCE[z].slice();
    return null;
}

/**
 * Whether a centre's total bond valence fits under its ceiling.
 *
 * The ceiling is deliberately the plain per-element maximum and not the
 * element's allowed valence set. Judging against the set sounds more correct —
 * a four-bonded nitrogen is a nitrogen with a positive charge, not an error —
 * but it makes the relaxation repair a centre by *raising* it (3.5 to 4 is a
 * legal value), and the half unit it adds lands on a neighbour, pushing that
 * one off its own value in turn. Measured on the bundled corpus that cascaded:
 * the count of centres off their value went the wrong way (120 to 156 on the
 * largest structure, 3 to 7 on another) and eleven bonds moved up where none
 * had before. Repairs that only ever remove excess do not have that failure
 * mode, and a genuine four-bonded nitrogen is still reachable — the valence
 * solver assigns it, with the charge that goes with it.
 */
function bpValenceOk(v: number, ceil: number): boolean {
    return v <= ceil + 0.1;
}

/**
 * How far a centre's valence sits from the nearest value its element may take,
 * counted in both directions.
 *
 * This is the measure the relaxation above deliberately does *not* use — it
 * repairs by removing excess, and judging it against the allowed set made it
 * repair by adding instead, which pushed neighbours off their own values. Here
 * it is the right question, because it is used only to compare two finished
 * assignments: it notices a half-unit that a ceiling comparison cannot (a
 * carbon at 3.5 is under 4 and yet is not a valence carbon takes), which is
 * exactly the kind of damage a re-assignment can do while looking legal.
 */
function bpValenceGap(z: number, v: number, ceil: number): number {
    const set = bpPossibleValences(z);
    if (!set) return Math.max(0, v - ceil);
    let best = Infinity;
    for (let i = 0; i < set.length; i++) {
        const gap = Math.abs(set[i] - v);
        if (gap < best) best = gap;
    }
    return best;
}

/**
 * Tighten each atom's candidate valences before the combination sweep.
 *
 * An atom of degree d can reach at most d plus the spare capacity of its
 * neighbours, so any candidate outside that band is impossible and is dropped.
 * Each tightening shortens the round after it, so the pass runs to a fixed
 * point instead of once.
 *
 * A frozen atom — one with no valence model, a metal — has exactly one
 * candidate: the degree it already has. That is not a target valence, it is the
 * statement that its bonds are not up for negotiation, and it is what lets the
 * rest of its component be solved without the metal being consulted.
 */
function bpPruneValences(zs: number[], deg: number[], nbr: number[][], frozen: boolean[]): number[][] | null {
    const k = zs.length;
    let cur: number[][] = new Array(k);
    for (let i = 0; i < k; i++) {
        if (frozen[i]) {
            cur[i] = [deg[i]];
            continue;
        }
        const v = bpPossibleValences(zs[i]);
        if (!v) return null;
        cur[i] = v;
    }

    for (let round = 0; round < 100; round++) {
        const next: number[][] = new Array(k);
        let same = true;
        for (let i = 0; i < k; i++) {
            if (cur[i].length === 1) {
                next[i] = cur[i].slice();
                if (next[i].length !== cur[i].length) same = false;
                continue;
            }
            let avail = 0;
            for (let q = 0; q < nbr[i].length; q++) {
                const j = nbr[i][q];
                let mx = 0;
                for (let t = 0; t < cur[j].length; t++) {
                    if (cur[j][t] > mx) mx = cur[j][t];
                }
                avail += mx - deg[j];
            }
            const keep: number[] = [];
            for (let t = 0; t < cur[i].length; t++) {
                const val = cur[i][t];
                if (val >= deg[i] && val <= deg[i] + avail) keep.push(val);
            }
            if (keep.length === 0) return null;
            next[i] = keep;
            if (keep.length !== cur[i].length) {
                same = false;
            } else {
                for (let t2 = 0; t2 < keep.length; t2++) {
                    if (keep[t2] !== cur[i][t2]) {
                        same = false;
                        break;
                    }
                }
            }
        }
        cur = next;
        if (same) break;
    }
    return cur;
}

/**
 * Formal charge implied by an element, its total bond valence and how many
 * neighbours it has. The generic branch counts the octet the atom has filled
 * against its outermost electrons; the exceptions are the centres the octet
 * cannot describe.
 *
 * The valence alone is not enough for three of them, and the neighbour count is
 * what makes the difference:
 *
 *   - **Nitrogen.** Four units of bond order make an ammonium when there are
 *     four neighbours and no charge at all when there are three: a pyrrole, an
 *     indole or an amide nitrogen contributes its lone pair to the ring, so the
 *     aromatic pair in its two ring bonds is electrons the atom already had, not
 *     a fourth bond. The octet formula charged every one of them +1.
 *   - **Sulphur.** Six of bond order is neutral whatever the split (a sulphone's
 *     two S=O, a sulphate's two S=O and two S–O); three units over two
 *     neighbours is a thiophene (neutral), three over three is the dative
 *     sulphoxide's S(+1)–O(−1), and four over three is a sulphoxide drawn with
 *     a double bond, which is neutral. The formula charged the dative cases +1
 *     and +2, which is a charge the file never stated.
 *   - **Phosphorus.** Five units is neutral (phosphine oxide, PF5); four is the
 *     phosphate dative count, which the formula already gets right.
 */
function bpAtomicCharge(z: number, valence: number, degree: number): number | null {
    if (z === 1) return 1 - valence;
    if (z === 5) return 3 - valence;
    if (z === 7 && valence === 4) return degree === 4 ? 1 : 0;
    if (z === 15 && valence === 5) return 0;
    if (z === 16) {
        if (valence === 6) return 0;
        if (valence === 3 && degree === 2) return 0;
        if (valence === 3 && degree === 3) return 1;
        if (valence === 4 && degree === 3) return 0;
    }
    const outer = OUTER_E[z];
    if (outer === undefined) return null;
    return outer - 8 + valence;
}

function bpCheckValency(order: number[], valency: number[]): boolean {
    for (let i = 0; i < valency.length; i++) {
        if (valency[i] > order[i]) return false;
    }
    return true;
}

function bpCheckSaturation(order: number[], valency: number[]): boolean {
    for (let i = 0; i < order.length; i++) {
        if (order[i] > valency[i]) return false;
    }
    return true;
}

/**
 * Whether the formal charges implied by a valence assignment sum to the
 * molecule's total charge.
 *
 * Two carbon cases are folded in — a two-coordinate carbon carrying a valence of
 * two is a carbene and contributes one extra unit, and a three-coordinate carbon
 * leaves one unit of its charge to the molecule — and both are decided from the
 * whole assignment rather than from a running total. The three-coordinate case
 * used to be tested as `total + 1 < charge` *while* the total was still being
 * accumulated, so whether a carbon contributed its extra units depended on how
 * many atoms happened to be counted before it: renumbering the same structure
 * could turn a component the solver accepted into one it rejected. The base
 * total is summed first and the carbons are added afterwards, in a loop whose
 * length is a count and not a position.
 *
 * An atom with no valence model (a metal) contributes nothing and is not
 * consulted: it has no formal charge the tables can state, and one iron must
 * not decide whether its component's organic part balances.
 */
function bpCheckCharge(zs: number[], degrees: number[], valency: number[], charge: number): boolean {
    let base = 0;
    let carbenes = 0;
    let threeCoord = 0;
    for (let i = 0; i < zs.length; i++) {
        const z = zs[i];
        if (!bpValenceModelled(z)) continue;
        const c = bpAtomicCharge(z, valency[i], degrees[i]);
        if (c === null) return false;
        base += c;
        if (z === 6) {
            if (degrees[i] === 2 && valency[i] === 2) carbenes++;
            else if (degrees[i] === 3) threeCoord++;
        }
    }
    let total = base + carbenes;
    for (let t = 0; t < threeCoord; t++) {
        if (total + 1 >= charge) break;
        total += 2;
    }
    return total === charge;
}

interface SolverResult {
    orders: number[][];
    valency: number[];
    exact: boolean;
}

/** Outcome of one solve attempt, including why it gave up. */
interface SolverAttempt {
    result: SolverResult | null;
    reason: string;
    combos: number;
}

/**
 * Search the valence combinations of one component for an assignment that
 * saturates every atom and matches the total charge.
 *
 * From the all-single-bond graph, atoms short of their target valence are
 * paired over the connectivity graph by a maximum-cardinality matching and
 * each matched pair gains one unit of order, repeated until no unsaturated
 * pair remains. An assignment is accepted when no atom exceeds its target and
 * all of them reach it; one that fits the ceiling but stays short is kept as
 * the best approximation, the largest total order winning. Nothing is returned
 * when no legal combination survives pruning, or when the sweep runs out of the
 * work budget (see `MAX_COMBO_WORK`) — the caller then keeps the geometric
 * result, and the reason says which of the two happened.
 *
 * **An atom with no valence model is frozen, not fatal.** A metal has no target
 * valence the tables can state, and the first version of this function read
 * that as "this component cannot be solved", which handed every structure
 * containing an iron back to the geometric pass — a ferrocene's cyclopentadienyl
 * rings, a magnesium complex's ligands, a zinc enzyme's whole organic part. The
 * atom now keeps the bonds the geometry gave it: its target is its own degree,
 * which makes it saturated by construction, so the matching never touches it and
 * no bond of its is moved. Only a component in which *every* atom is unmodelled
 * still declines, which is the right answer for a bare metal cluster.
 */
function bpDetermineOrders(
    zs: number[],
    edges: Array<[number, number]>,
    charge: number,
    k: number
): SolverAttempt {
    const frozen: boolean[] = new Array(k);
    let modelled = 0;
    for (let i = 0; i < k; i++) {
        frozen[i] = !bpValenceModelled(zs[i]);
        if (!frozen[i]) modelled++;
    }
    if (modelled === 0) return { result: null, reason: 'unmodelled', combos: 0 };

    const con: number[][] = [];
    const val0: number[] = new Array(k);
    const nbr: number[][] = [];
    for (let i = 0; i < k; i++) {
        const row: number[] = new Array(k);
        for (let t = 0; t < k; t++) row[t] = 0;
        con.push(row);
        val0[i] = 0;
        nbr.push([]);
    }
    for (let e = 0; e < edges.length; e++) {
        const a = edges[e][0];
        const b = edges[e][1];
        con[a][b]++;
        con[b][a]++;
        val0[a]++;
        val0[b]++;
        nbr[a].push(b);
        nbr[b].push(a);
    }
    const deg = val0.slice();

    const combos = bpPruneValences(zs, deg, nbr, frozen);
    if (!combos) return { result: null, reason: 'no-legal-valence', combos: 0 };

    let total = 1;
    for (let i = 0; i < k; i++) {
        if (combos[i].length === 0) return { result: null, reason: 'no-legal-valence', combos: 0 };
        total *= combos[i].length;
        if (total > MAX_COMBOS) return { result: null, reason: 'combo-overflow', combos: total };
    }

    let best: SolverResult | null = null;
    let bestSum = -1;
    let work = 0;
    const idx: number[] = new Array(k);
    for (let i = 0; i < k; i++) idx[i] = 0;

    for (let iter = 0; iter < total; iter++) {
        // The budget is charged before the combination is looked at, so a
        // component whose combinations are individually expensive stops at the
        // first one that would overrun instead of after a hundred of them.
        work += k;
        if (work > MAX_COMBO_WORK) {
            return {
                result: best,
                reason: best ? 'approximate' : 'work-budget',
                combos: total
            };
        }

        const order: number[] = new Array(k);
        for (let i = 0; i < k; i++) order[i] = combos[i][idx[i]];

        let unsat = false;
        for (let i = 0; i < k; i++) {
            if (order[i] > val0[i]) {
                unsat = true;
                break;
            }
        }

        if (!unsat) {
            if (bpCheckValency(order, val0) && bpCheckCharge(zs, deg, val0, charge) && bpCheckSaturation(order, val0)) {
                return { result: { orders: con, valency: val0, exact: true }, reason: 'ok', combos: total };
            }
        } else {
            const ordMat: number[][] = new Array(k);
            for (let i = 0; i < k; i++) ordMat[i] = con[i].slice();
            const val = val0.slice();
            work += k * k;

            // The matching round count is bounded by the structure as well as
            // by the constant: every round that changes something adds at least
            // one unit of bond order, and the total is bounded by the sum of the
            // ceilings.
            const guardLimit = Math.min(200, 4 + 2 * k);
            for (let guard = 0; guard < guardLimit; guard++) {
                let changed = false;
                const ua: number[] = [];
                for (let i = 0; i < k; i++) {
                    if (order[i] > val[i]) ua.push(i);
                }
                const pairs: Array<[number, number]> = [];
                for (let a = 0; a < ua.length; a++) {
                    for (let b = a + 1; b < ua.length; b++) {
                        if (con[ua[a]][ua[b]] > 0) pairs.push([ua[a], ua[b]]);
                    }
                }
                if (pairs.length) {
                    work += k + ua.length * ua.length;
                    if (work > MAX_COMBO_WORK) break;
                    const mate = bpMaximumMatching(k, pairs);
                    for (let v = 0; v < k; v++) {
                        if (mate[v] !== -1 && v < mate[v]) {
                            changed = true;
                            ordMat[v][mate[v]]++;
                            val[v]++;
                            val[mate[v]]++;
                        }
                    }
                }
                if (!changed) break;
            }

            if (bpCheckValency(order, val) && bpCheckCharge(zs, deg, val, charge)) {
                if (bpCheckSaturation(order, val)) {
                    return { result: { orders: ordMat, valency: val, exact: true }, reason: 'ok', combos: total };
                }
                let s = 0;
                for (let i = 0; i < k; i++) s += val[i];
                if (s > bestSum) {
                    bestSum = s;
                    best = { orders: ordMat, valency: val, exact: false };
                }
            }
        }

        for (let i = 0; i < k; i++) {
            idx[i]++;
            if (idx[i] < combos[i].length) break;
            idx[i] = 0;
        }
    }
    return {
        result: best,
        reason: best ? 'approximate' : 'no-assignment',
        combos: total
    };
}

// ---------------------------------------------------------------------------
// 6. Local relaxation (the geometric fallback)
// ---------------------------------------------------------------------------

/** Candidate orders for one bond that the measured distance actually supports. */
function bpRelaxLadder(el1: string, el2: string, d: number): number[] {
    const specs = bpSpecList(el1, el2);
    const out: number[] = [];
    if (!specs) return out;
    let best = -1;
    let bestDiff = Infinity;
    for (let i = 0; i < specs.length; i++) {
        const diff = Math.abs(d - specs[i].l);
        if (diff <= specs[i].t) out.push(i);
        if (diff < bestDiff) {
            bestDiff = diff;
            best = i;
        }
    }
    if (best >= 0 && out.indexOf(best) < 0) out.push(best);
    return out;
}

/**
 * Re-plan one over-coordinated centre from scratch, reversibly.
 *
 * The question is not which neighbour to downgrade but which combination of
 * allowed orders puts the neighbours closest to the lengths they were measured
 * at. Candidate orders are limited to the values the distance supports, so a
 * long non-bonded contact can never be promoted into a multiple bond; a plan
 * must keep the centre within its ceiling, must not push any neighbour over
 * its own (each neighbour contributes its headroom), and is accepted only when
 * it strictly lowers the total deviation. Both directions of movement are
 * available, and a balanced pair move is searched alongside single transfers.
 * Terminal neighbours are pinned to a single bond and excluded from the
 * budget. Returns whether anything changed.
 *
 * An over-ceiling centre is searched differently, because it is over by an
 * amount that is not always one unit: a sulphate's sulphur with four S=O bonds
 * measured at 1.47 Å is two units above its ceiling of six, and a plan that
 * moves one neighbour at a time can never reach it — the first version of this
 * pass could only ever shed a single unit, so the largest repair it could make
 * was from 7 to 6 and a centre at 8 stayed at 8. The over-ceiling branch
 * therefore searches subsets of up to three neighbours, one step down each,
 * and keeps the fitting plan whose bonds sit closest to their lengths. Three
 * is enough for every polyvalent centre the tables contain (a sulphate needs
 * two, a perphosphonate two, a sulphone none) and bounds the search at a few
 * hundred plans for the widest coordination.
 */
function bpRelaxCenter(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[],
    i: number,
    entries: Array<[number, number]>,
    curDev: number,
    curSum: number
): boolean {
    const el = bpNorm(atoms[i].element);
    const maxV = bpCeiling(el);

    const free: number[] = [];
    let terminalSum = 0;
    for (let k = 0; k < entries.length; k++) {
        const j = entries[k][0];
        if (bpIsTerminal(atoms[j].element)) terminalSum += 1.0;
        else free.push(k);
    }
    if (free.length === 0) return false;
    const nf = free.length;

    // Candidate orders and their deviations, built together so the two arrays
    // can never get out of step. They used to be filled by two separate pushes,
    // keyed only by the assumption that every candidate is new: for a pair whose
    // table lists the same order twice (or a current order already on the
    // ladder) the deviation of one bond was read for another.
    const cand: number[][] = [];
    const devs: number[][] = [];
    const curIdx: number[] = [];
    for (let fi = 0; fi < nf; fi++) {
        const j = entries[free[fi]][0];
        const cur = entries[free[fi]][1];
        const d = bpDist(atoms[i], atoms[j]);
        const specs = bpSpecList(el, atoms[j].element);
        const ladder: number[] = [];
        const dlist: number[] = [];
        const add = function (o: number, dev: number): void {
            if (!(o > 0)) return;
            const at = ladder.indexOf(o);
            if (at >= 0) {
                if (dev < dlist[at]) dlist[at] = dev;
                return;
            }
            ladder.push(o);
            dlist.push(dev);
        };
        if (specs) {
            const idxs = bpRelaxLadder(el, atoms[j].element, d);
            for (let q = 0; q < idxs.length; q++) {
                add(specs[idxs[q]].o, Math.abs(d - specs[idxs[q]].l));
            }
            if (ladder.indexOf(cur) < 0) {
                // The order already assigned, when the length no longer argues
                // for it: it stays reachable so a plan can keep it, priced at
                // the deviation of the nearest tabulated order.
                let cd = Infinity;
                for (let q = 0; q < specs.length; q++) {
                    const df = Math.abs(d - specs[q].l);
                    if (df < cd) cd = df;
                }
                add(cur, cd);
            }
        } else {
            add(cur, 0);
        }
        cand.push(ladder);
        devs.push(dlist);
        const at = ladder.indexOf(cur);
        curIdx.push(at >= 0 ? at : 0);
    }

    const headroom: number[] = [];
    for (let fi = 0; fi < nf; fi++) {
        const j = entries[free[fi]][0];
        const cap = bpCeiling(atoms[j].element);
        const own = bondMap.get(j);
        let used = 0;
        if (own) {
            own.forEach(function (o, k2) {
                if (k2 !== i) used += o;
            });
        }
        headroom.push(cap - used);
    }

    let chosen = curIdx.slice();
    let bestSigned = 0;
    let bestAbs = 0;
    let found = false;

    // A centre whose valence is not one its element may take is judged on
    // valence first: any plan landing on an allowed value is acceptable, and
    // among those the one leaving its neighbours closest to their measured
    // lengths wins. That is the only way such a centre can be repaired at all —
    // the nitrogen measured at 3.5 because one of its bonds fell inside the
    // pair's aromatic tolerance has to give that half unit up even though the
    // distance is happier with it. A centre already sitting on an allowed value
    // keeps the stricter rule: a plan must lower the total deviation, which is
    // what makes the sweep idempotent and stops it inventing multiple bonds
    // out of long contacts.
    const overCenter = !bpValenceOk(curSum, maxV);

    const evaluate = function (plan: number[]): void {
        let sum = terminalSum;
        for (let k = 0; k < nf; k++) sum += cand[k][plan[k]];
        if (!bpValenceOk(sum, maxV)) return;

        let signed = 0;
        let abs = 0;
        for (let k = 0; k < nf; k++) {
            const v = cand[k][plan[k]];
            if (v - cand[k][curIdx[k]] > headroom[k] + 0.001) return;
            signed += devs[k][plan[k]] - devs[k][curIdx[k]];
            abs += devs[k][plan[k]];
        }

        if (overCenter) {
            if (!found || abs < bestAbs - 1e-9) {
                bestSigned = signed;
                bestAbs = abs;
                chosen = plan.slice();
                found = true;
            }
            return;
        }

        if (sum < curSum - 0.001) return;
        if (!found || signed < bestSigned - 1e-9 || (Math.abs(signed - bestSigned) <= 1e-9 && abs < bestAbs - 1e-9)) {
            bestSigned = signed;
            bestAbs = abs;
            chosen = plan.slice();
            found = true;
        }
    };

    const scratch = curIdx.slice();

    // Phase 0 — shed units, changing the centre's total. The balanced phases
    // below always leave the total where it was, so without this a centre above
    // its ceiling has nothing to reach for: every plan they offer still adds up
    // to the same over-large sum. Only an over-ceiling centre runs this phase,
    // because moving a neighbour with nothing balancing it is exactly what the
    // strict-improvement rule forbids everywhere else.
    if (overCenter) {
        const plan = curIdx.slice();
        const walk = function (from: number, depth: number): void {
            if (depth > 0) evaluate(plan);
            if (depth >= MAX_SHED_BONDS) return;
            for (let fi = from; fi < nf; fi++) {
                for (let q = 0; q < cand[fi].length; q++) {
                    if (cand[fi][q] >= cand[fi][curIdx[fi]]) continue;
                    const keep = plan[fi];
                    plan[fi] = q;
                    walk(fi + 1, depth + 1);
                    plan[fi] = keep;
                }
            }
        };
        walk(0, 0);
    } else {
        // Phase 1 — one neighbour, any direction.
        for (let fi = 0; fi < nf; fi++) {
            for (let i1 = 0; i1 < cand[fi].length; i1++) {
                if (i1 === curIdx[fi]) continue;
                scratch[fi] = i1;
                evaluate(scratch);
                scratch[fi] = curIdx[fi];
            }
        }

        // Phase 2 and 3 — a balanced pair, so the centre's total is unchanged
        // while its neighbours' lengths are served better. Both orderings of the
        // pair are tried, because the first loop may move the *earlier* bond up
        // while the later one comes down.
        for (let fi = 0; fi < nf; fi++) {
            for (let i1 = 0; i1 < cand[fi].length; i1++) {
                const delta = cand[fi][i1] - cand[fi][curIdx[fi]];
                if (delta === 0) continue;
                for (let fj = 0; fj < nf; fj++) {
                    if (fj === fi) continue;
                    for (let i2 = 0; i2 < cand[fj].length; i2++) {
                        if (Math.abs(cand[fj][i2] - cand[fj][curIdx[fj]] + delta) > 0.001) continue;
                        scratch[fi] = i1;
                        scratch[fj] = i2;
                        evaluate(scratch);
                        scratch[fi] = curIdx[fi];
                        scratch[fj] = curIdx[fj];
                    }
                }
            }
        }
    }

    // Outside the over-ceiling case a plan that is not strictly better than
    // leaving the atom alone is rejected, which is what makes the whole sweep
    // idempotent. An over-ceiling centre is accepted as soon as one of its
    // plans fits, because no plan that does not fit was ever allowed above.
    if (!found) return false;
    if (!overCenter && bestSigned > -1e-9) return false;

    let improved = false;
    for (let fi = 0; fi < nf; fi++) {
        const k = free[fi];
        const next = cand[fi][chosen[fi]];
        if (next === entries[k][1]) continue;
        const j = entries[k][0];
        const ni = bondMap.get(i);
        const nj = bondMap.get(j);
        if (ni) ni.set(j, next);
        if (nj) nj.set(i, next);
        improved = true;
    }
    return improved;
}

/**
 * Bring every centre whose valence is not one its element may take back to a
 * value it may take, to a fixed point.
 *
 * A worklist, not a full sweep per round. Only a centre that is off its
 * valence — and a neighbour an accepted plan has just moved — is ever
 * examined, so the cost follows the atoms that actually need attention instead
 * of the size of the structure. The distinction is not academic: a sweep per
 * round is O(rounds x atoms), and the previous shape of this pass got away
 * with it only because it could almost never repair anything, so it stopped
 * after a single round. Once an over-coordinated centre can actually be
 * brought down, a pathological input — thousands of atoms clustered inside
 * bonding range — keeps finding work and re-scans everything each round: that
 * measured a 4.7x regression on a 2000-atom cluster (0.87 s to 4.06 s) while
 * the real corpus moved by 1.12x.
 *
 * The number of examinations is capped as well, so no input can make the loop
 * unbounded; a structure with more violations than the cap allows simply keeps
 * whatever the earlier visits settled.
 */
function bpRelaxValence(bondMap: Map<number, Map<number, number>>, atoms: Atom[]): void {
    const n = atoms.length;
    const queue: number[] = [];
    const queued: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) queued[i] = false;

    const off = function (i: number): boolean {
        const nb = bondMap.get(i);
        if (!nb || nb.size === 0) return false;
        let sum = 0;
        nb.forEach(function (o) {
            sum += o;
        });
        return !bpValenceOk(sum, bpCeiling(atoms[i].element));
    };

    for (let i = 0; i < n; i++) {
        if (bondMap.has(i) && off(i)) {
            queue.push(i);
            queued[i] = true;
        }
    }
    if (queue.length === 0) return;

    // A budget on examinations, not on rounds: bounded so that a degenerate
    // input cannot spin, and generous enough that a real structure never
    // reaches it.
    const budget = Math.min(200000, 8 * n + 2000);
    let head = 0;
    let spent = 0;

    while (head < queue.length && spent < budget) {
        const i = queue[head];
        head++;
        queued[i] = false;
        spent++;
        if (!off(i)) continue;

        const nb = bondMap.get(i);
        if (!nb || nb.size === 0) continue;
        const entries: Array<[number, number]> = [];
        nb.forEach(function (o, j) {
            entries.push([j, o]);
        });
        entries.sort(function (a, b) {
            return a[0] - b[0];
        });

        let curSum = 0;
        let curDev = 0;
        for (let e = 0; e < entries.length; e++) {
            curSum += entries[e][1];
            const d = bpDist(atoms[i], atoms[entries[e][0]]);
            const specs = bpSpecList(atoms[i].element, atoms[entries[e][0]].element);
            let ideal = d;
            if (specs) {
                let minD = Infinity;
                for (let s = 0; s < specs.length; s++) {
                    if (specs[s].o === entries[e][1] && Math.abs(d - specs[s].l) < minD) {
                        minD = Math.abs(d - specs[s].l);
                        ideal = specs[s].l;
                    }
                }
            }
            curDev += Math.abs(d - ideal);
        }

        if (!bpRelaxCenter(bondMap, atoms, i, entries, curDev, curSum)) continue;

        // The accepted plan moved bonds, so every neighbour whose order changed
        // is a candidate for needing its own repair; and the centre itself may
        // still be off its valence.
        for (let e = 0; e < entries.length; e++) {
            const j = entries[e][0];
            const after = nb.get(j);
            if (after === undefined || after === entries[e][1] || queued[j]) continue;
            queue.push(j);
            queued[j] = true;
        }
        if (!queued[i] && off(i)) {
            queue.push(i);
            queued[i] = true;
        }
    }
}

// 7. Putting the geometric verdict back
// ---------------------------------------------------------------------------

/** Total bond order around one atom. */
function bpValenceOf(bondMap: Map<number, Map<number, number>>, i: number): number {
    const nb = bondMap.get(i);
    if (!nb) return 0;
    let sum = 0;
    nb.forEach(function (o) {
        sum += o;
    });
    return sum;
}

/** Total distance of a component's valences from the values they may take. */
function bpComponentGap(bondMap: Map<number, Map<number, number>>, atoms: Atom[], comp: number[]): number {
    let total = 0;
    for (let c = 0; c < comp.length; c++) {
        const i = comp[c];
        total += bpValenceGap(bpZOf(atoms[i].element), bpValenceOf(bondMap, i), bpCeiling(atoms[i].element));
    }
    return total;
}

/**
 * Restore the aromatic bonds a solver pass could only express as an
 * alternating single/double pattern.
 *
 * The geometric pass recognises an aromatic bond from its length, and it is the
 * only thing in this module that can produce order 1.5. The valence solver
 * works on the graph and can only produce whole orders, so a ring it re-derived
 * comes back as alternating single and double bonds — the same chemistry, drawn
 * as a Kekulé structure. Putting the 1.5 back is a **group** move by
 * construction: raising one ring bond on its own would change a ring atom's
 * valence, so the bonds are grouped by shared atom and a group is applied only
 * when the changes cancel *exactly* at every atom it touches — which is
 * precisely the case for a ring whose bonds alternate.
 *
 * The cancellation is the whole point, and it is what makes this pass safe by
 * construction rather than merely tidy: because every atom ends up exactly
 * where it was, no atom can be moved off the valence already established, so a
 * component either gains its aromatic rings back or keeps the answer it had. A
 * restore without that condition — one that put each bond back on the order
 * its own length argued for, independently — would not have the property:
 * undoing a bond that was lowered to pay for an over-coordinated centre puts
 * that centre straight back over its ceiling, which costs more than the ring
 * notation is worth, and it leaves rings half-restored besides. Hence the
 * group condition, and hence a group that does not balance is left precisely
 * as it was found.
 */
function bpRestoreAromatic(bondMap: Map<number, Map<number, number>>, atoms: Atom[], comp: number[]): number {
    const inComp: { [key: number]: number } = {};
    for (let c = 0; c < comp.length; c++) inComp[comp[c]] = 1;

    // Bonds the geometry calls aromatic but an earlier pass lowered.
    const cand: Array<[number, number, number]> = [];
    for (let c = 0; c < comp.length; c++) {
        const i = comp[c];
        const nb = bondMap.get(i);
        if (!nb) continue;
        const partners: number[] = [];
        nb.forEach(function (o, j) {
            partners.push(j);
        });
        partners.sort(function (a, b) {
            return a - b;
        });
        for (let p = 0; p < partners.length; p++) {
            const j = partners[p];
            if (j <= i || !inComp[j]) continue;
            const cur = nb.get(j);
            if (cur === undefined || Math.abs(cur - 1.5) < 1e-9) continue;
            if (bpGeometricOrder(atoms[i].element, atoms[j].element, bpDist(atoms[i], atoms[j])) !== 1.5) continue;
            cand.push([i, j, cur]);
        }
    }
    if (cand.length === 0) return 0;

    // Group the candidates by shared atom.
    const group: { [key: number]: number } = {};
    let groups = 0;
    for (let e = 0; e < cand.length; e++) {
        const a = cand[e][0];
        const b = cand[e][1];
        const ga = group[a];
        const gb = group[b];
        if (ga === undefined && gb === undefined) {
            group[a] = groups;
            group[b] = groups;
            groups++;
        } else if (ga === undefined) {
            group[a] = gb as number;
        } else if (gb === undefined) {
            group[b] = ga;
        } else if (ga !== gb) {
            for (var k in group) {
                if (group[k] === gb) group[k] = ga;
            }
        }
    }

    // Net change each group would make at each atom it touches.
    const net: Array<{ [key: number]: number }> = [];
    for (let g = 0; g < groups; g++) net.push({});
    for (let e = 0; e < cand.length; e++) {
        const g = group[cand[e][0]];
        const shift = 1.5 - cand[e][2];
        const ng = net[g];
        ng[cand[e][0]] = (ng[cand[e][0]] || 0) + shift;
        ng[cand[e][1]] = (ng[cand[e][1]] || 0) + shift;
    }
    const balanced: boolean[] = [];
    for (let g = 0; g < groups; g++) {
        let ok = true;
        const ng = net[g];
        for (var ka in ng) {
            if (Math.abs(ng[+ka]) > 1e-6) {
                ok = false;
                break;
            }
        }
        balanced.push(ok);
    }

    let restored = 0;
    for (let e = 0; e < cand.length; e++) {
        if (!balanced[group[cand[e][0]]]) continue;
        const na = bondMap.get(cand[e][0]);
        const nb2 = bondMap.get(cand[e][1]);
        if (na) na.set(cand[e][1], 1.5);
        if (nb2) nb2.set(cand[e][0], 1.5);
        restored++;
    }
    return restored;
}

/**
 * A last look at valence, once the orders are settled.
 *
 * For every atom whose valence is not a value its element may take, try
 * adjusting the orders of the bonds around it, and keep the adjustment only
 * when the whole neighbourhood comes out strictly better.
 *
 * The measure is the one the diagnostic reports and the safety valve compares:
 * the total distance of an atom's valence from the values its element may take.
 * What differs is the target. The local relaxation removes excess and stops
 * there, because raising an order to fill an under-filled centre would push the
 * half unit onto a neighbour and leave that one wrong in turn. Here the
 * neighbourhood is scored as a whole, and a move is kept only when the score
 * improves — so a repair that costs a neighbour more than it gains is rejected
 * on the spot rather than trusted. That is what makes it safe to run by
 * default: it cannot make a structure worse, because "worse" is the thing it
 * refuses.
 *
 * Two rules keep it small:
 *
 *   - **Hydrogen is pinned.** A hydrogen has one bond and its order is 1, so
 *     every bond to a hydrogen is held fixed. That removes the largest source
 *     of pointless churn — a large share of the bonds in a typical structure —
 *     and it means the pass never reconsiders a valence that was never in
 *     question.
 *   - **One neighbourhood at a time.** A centre and the atoms it bonds to are
 *     scored together and never beyond. The corrections this looks for are
 *     local, and a wider search would cost more than it returns.
 *
 * A bond may only take a value its own length supports (or keep the value it
 * has), so an order is never invented out of a distance that argues against it.
 *
 * Returns how many bonds it moved.
 */
function bpRefineValences(bondMap: Map<number, Map<number, number>>, atoms: Atom[]): number {
    const n = atoms.length;
    let moved = 0;
    // Examinations are capped, so a structure where nothing can be improved
    // costs a bounded amount and then stops.
    const budget = Math.min(4000, 8 * n + 200);
    let spent = 0;
    // Which atoms are hydrogen, decided once instead of per bond examined.
    const isH: boolean[] = new Array(n);
    for (let i0 = 0; i0 < n; i0++) isH[i0] = bpNorm(atoms[i0].element) === 'H';

    const gapOf = function (i: number): number {
        return bpValenceGap(bpZOf(atoms[i].element), bpValenceOf(bondMap, i), bpCeiling(atoms[i].element));
    };
    const scoreOf = function (i: number, around: number[]): number {
        let s = gapOf(i);
        for (let k = 0; k < around.length; k++) s += gapOf(around[k]);
        return s;
    };
    // How far the centre's own bonds sit from the lengths the table expects.
    // A move that fixes a valence but drags a bond well off its length is not
    // a gain, so the neighbourhood is only allowed to trade one for the other
    // when the length fit does not get worse. A pair with no table entry has no
    // fit to compare and contributes nothing — which is why the *admission* of
    // a candidate below cannot be decided here.
    const devOf = function (i: number, around: number[]): number {
        const nb = bondMap.get(i);
        if (!nb) return 0;
        let s = 0;
        const a = atoms[i];
        for (let k = 0; k < around.length; k++) {
            const j = around[k];
            const o = nb.get(j);
            if (o === undefined) continue;
            const dev = bpSpecDev(a.element, atoms[j].element, o, bpDist(a, atoms[j]));
            if (dev >= 0 && dev < Infinity) s += dev;
        }
        return s;
    };
    /**
     * Whether the pair's own table says a move to a new order fits this length
     * at least as well as the current one does.
     *
     * A pair with no table at all answers *no*, and that is the whole point.
     * bpSpecDev used to answer 0 for such a pair — "perfect fit" — so a
     * silicon–oxygen or sulphur–oxygen bond, neither of which had an entry, was
     * freely movable by half a unit and came out at 1.5. No evidence is not
     * evidence for every order; it is evidence for none, and the bond keeps the
     * order the geometry gave it.
     */
    const moveAllowed = function (el1: string, el2: string, d: number, cur: number, next: number): boolean {
        const curDev = bpSpecDev(el1, el2, cur, d);
        if (!(curDev >= 0)) return false;
        const dev = bpSpecDev(el1, el2, next, d);
        return dev >= 0 && dev <= curDev + 1e-9;
    };
    const setOrder = function (i: number, j: number, v: number): void {
        const ni = bondMap.get(i);
        if (ni) ni.set(j, v);
        const nj = bondMap.get(j);
        if (nj) nj.set(i, v);
    };

    for (let round = 0; round < 4; round++) {
        let any = false;
        for (let i = 0; i < n && spent < budget; i++) {
            if (isH[i]) continue;
            const nb = bondMap.get(i);
            if (!nb || nb.size < 2) continue;
            if (gapOf(i) <= 1e-6) continue;
            spent++;

            // Cheap probe before anything is allocated: count how many ways
            // the movable bonds could be arranged. When every one of them has
            // exactly one candidate — its current order, and nothing else its
            // length is at least as happy with — there is nothing to try, and
            // that is the common case for a centre that is off its valence for
            // a reason no re-arrangement can fix. Skipping here is what stops
            // such a structure from paying for the enumeration below.
            let choices = 1;
            nb.forEach(function (o, j) {
                if (choices > 4096) return;
                if (isH[j]) return;
                const d = bpDist(atoms[i], atoms[j]);
                const s = bpSupportedOrder(atoms[i].element, atoms[j].element, d);
                let c = 1;
                if (s > 0 && Math.abs(s - o) > 1e-9) c++;
                for (let t = 0; t < 2; t++) {
                    const v = o + (t === 0 ? -0.5 : 0.5);
                    if (v < 1 || v > 3) continue;
                    if (Math.abs(v - s) < 1e-9) continue;
                    if (moveAllowed(atoms[i].element, atoms[j].element, d, o, v)) c++;
                }
                choices *= c;
            });
            if (choices <= 1) continue;

            const around: number[] = [];
            const movable: number[] = [];
            nb.forEach(function (o, j) {
                around.push(j);
                if (!isH[j]) movable.push(j);
            });
            if (movable.length === 0) continue;

            const before = scoreOf(i, around);

            // Candidate orders per movable bond: what its own length supports,
            // plus whatever it currently is, so a repair never depends on the
            // length table having an entry for the pair.
            const cands: number[][] = [];
            const orig: number[] = [];
            let combos = 1;
            for (let m = 0; m < movable.length; m++) {
                const j = movable[m];
                const cur = nb.get(j) as number;
                orig.push(cur);
                const d = bpDist(atoms[i], atoms[j]);
                const s = bpSupportedOrder(atoms[i].element, atoms[j].element, d);
                const list: number[] = [];
                if (s > 0) list.push(s);
                if (list.indexOf(cur) < 0) list.push(cur);
                // The neighbouring steps as well. A centre can be off its
                // valence because a bond sits a whole unit away from what the
                // geometry alone would pick, and the only way to see whether
                // paying that unit is worth it is to try it. The score decides
                // — but only among values that do not fit the length *worse*
                // than the current one, which is the same condition the sweep
                // applies at the end. Filtering here rather than there is what
                // keeps the pass cheap: on a structure where every centre is
                // off its valence for a reason no re-arrangement can fix, every
                // bond drops to a single candidate, the combination count
                // collapses to one, and the whole neighbourhood is skipped
                // instead of enumerated.
                const extra = [cur - 0.5, cur + 0.5];
                for (let t = 0; t < extra.length; t++) {
                    const v = extra[t];
                    if (v < 1 || v > 3 || list.indexOf(v) >= 0) continue;
                    if (moveAllowed(atoms[i].element, atoms[j].element, d, cur, v)) list.push(v);
                }
                cands.push(list);
                combos *= list.length;
            }
            if (combos <= 1 || combos > 4096) continue;

            // Sweep the combinations in place, so every score is read from the
            // real state, then put the neighbourhood back the way it was.
            const beforeDev = devOf(i, around);
            let bestScore = before;
            const pick: number[] = [];
            const cur: number[] = [];
            for (let m = 0; m < movable.length; m++) cur.push(0);
            const walk = function (m: number): void {
                if (m === movable.length) {
                    for (let k = 0; k < movable.length; k++) setOrder(i, movable[k], cands[k][cur[k]]);
                    const sc = scoreOf(i, around);
                    // A gain in valence that costs the bond lengths is not a
                    // gain — the lengths are the measurement, and a repair is
                    // only worth taking when it is free in that currency.
                    if (sc < bestScore - 1e-9 && devOf(i, around) <= beforeDev + 1e-9) {
                        bestScore = sc;
                        pick.length = 0;
                        for (let k = 0; k < movable.length; k++) pick.push(cur[k]);
                    }
                    return;
                }
                for (let c = 0; c < cands[m].length; c++) {
                    cur[m] = c;
                    walk(m + 1);
                }
            };
            walk(0);

            for (let m = 0; m < movable.length; m++) setOrder(i, movable[m], orig[m]);
            if (pick.length !== movable.length) continue;

            for (let m = 0; m < movable.length; m++) {
                const v = cands[m][pick[m]];
                if (Math.abs(v - orig[m]) > 1e-9) {
                    setOrder(i, movable[m], v);
                    moved++;
                }
            }
            any = true;
        }
        if (!any) break;
    }
    return moved;
}

// ---------------------------------------------------------------------------
// 8. Components and diagnostics
// ---------------------------------------------------------------------------

/** Connected components of the bond graph, as ascending index lists. */
function bpComponents(bondMap: Map<number, Map<number, number>>, n: number): number[][] {
    const seen: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) seen[i] = false;
    const comps: number[][] = [];
    for (let start = 0; start < n; start++) {
        if (seen[start] || !bondMap.has(start)) continue;
        const comp: number[] = [];
        const stack: number[] = [start];
        seen[start] = true;
        while (stack.length) {
            const i = stack.pop() as number;
            comp.push(i);
            const nb = bondMap.get(i);
            if (!nb) continue;
            nb.forEach(function (o, j) {
                if (!seen[j]) {
                    seen[j] = true;
                    stack.push(j);
                }
            });
        }
        comp.sort(function (a, b) {
            return a - b;
        });
        comps.push(comp);
    }
    return comps;
}

/** Whether any atom of the component sits above its valence ceiling. */
function bpHasViolation(bondMap: Map<number, Map<number, number>>, atoms: Atom[], comp: number[]): boolean {
    for (let c = 0; c < comp.length; c++) {
        const i = comp[c];
        const nb = bondMap.get(i);
        if (!nb) continue;
        if (!bpValenceOk(bpValenceOf(bondMap, i), bpCeiling(atoms[i].element))) return true;
    }
    return false;
}

interface ComponentSolution {
    orders: Array<[number, number, number]>;
    charges: number[];
    valency: number[];
    exact: boolean;
}

interface ComponentAttempt {
    solution: ComponentSolution | null;
    reason: string;
    combos: number;
}

/**
 * Re-derive one component with the valence solver, in the component's own
 * index space so the matching never sees an index it does not own. Returns
 * null when the component has no legal assignment at all.
 *
 * The index map is the caller's (`localOf`), not a fresh array: it used to be
 * allocated here, `new Array(atoms.length)` for *every* component, so a
 * structure of 3 000 fragments paid 3 000 × 3 000 slots for a map that holds one
 * component at a time. It is allocated once per pass and the entries this
 * component wrote are cleared before returning, which leaves it exactly as it
 * was found.
 */
function bpSolveComponent(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[],
    comp: number[],
    charge: number,
    localOf: number[]
): ComponentAttempt {
    const k = comp.length;
    for (let i = 0; i < k; i++) localOf[comp[i]] = i;

    const zs: number[] = new Array(k);
    for (let i = 0; i < k; i++) {
        zs[i] = bpZOf(atoms[comp[i]].element);
    }

    const edges: Array<[number, number]> = [];
    for (let i = 0; i < k; i++) {
        const nb = bondMap.get(comp[i]);
        if (!nb) continue;
        const partners: number[] = [];
        nb.forEach(function (o, j) {
            partners.push(j);
        });
        for (let p = 0; p < partners.length; p++) {
            const lj = localOf[partners[p]];
            if (lj < 0 || lj <= i) continue;
            edges.push([i, lj]);
        }
    }

    const attempt = bpDetermineOrders(zs, edges, charge, k);

    // Hand the map back in the state it arrived in, whatever the outcome: a
    // stale entry would make a later component see an index it does not own.
    for (let i = 0; i < k; i++) localOf[comp[i]] = -1;
    if (!attempt.result) {
        return { solution: null, reason: attempt.reason, combos: attempt.combos };
    }

    const res = attempt.result;
    const orders: Array<[number, number, number]> = [];
    for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
            if (res.orders[i][j] > 0) orders.push([comp[i], comp[j], res.orders[i][j]]);
        }
    }

    // Neighbours as well as the valence: the formal charge of a nitrogen or a
    // sulphur is not a function of its bond order alone (pyrrole N and
    // pyridinium N have the same sum over bonds and differ only in how many
    // neighbours it takes), so the degree has to be the component's own.
    const degree = new Array<number>(k).fill(0);
    for (let e = 0; e < edges.length; e++) {
        degree[edges[e][0]]++;
        degree[edges[e][1]]++;
    }
    const charges: number[] = new Array(k);
    for (let i = 0; i < k; i++) {
        const c = bpAtomicCharge(zs[i], res.valency[i], degree[i]);
        charges[i] = c === null ? 0 : c;
    }

    return {
        solution: { orders: orders, charges: charges, valency: res.valency, exact: res.exact },
        reason: attempt.reason,
        combos: attempt.combos
    };
}

// ---------------------------------------------------------------------------
// 9. Entry point
// ---------------------------------------------------------------------------

export interface PerceptionOptions {
    /** Total molecular charge, when the file states one. */
    charge?: number;
    /**
     * Hand components the geometric pass left off their allowed valence to the
     * valence solver.
     *
     * Off by default, and that default is the point: the geometric pass is the
     * viewer's own behaviour, and the solver answers a *component*-wide
     * question — it re-derives every bond in the component, so it can move
     * bonds that were never in doubt (a ring comes back as an alternating
     * single/double pattern before the aromatic restoration puts it right, and
     * a bond next to a repaired centre is re-assigned to keep the books
     * balanced). That is more correct chemistry, but it is a visible change
     * the default has no business making on its own. Diagnostics always report
     * how many components the solver *could* settle, so a caller that wants
     * them can ask.
     */
    solve?: boolean;
    /** Components above this atom count keep the geometric result. */
    maxSolverAtoms?: number;
    /**
     * Also produce the diagnostics: the off-valence list and the components the
     * solver could settle.
     *
     * Off by default, because the common caller wants the bond list and nothing
     * else — a re-derivation after an edit, the webview's client-side entry
     * point, a structure being loaded. Those extra results are O(atoms) work
     * over the finished result (a component walk, a valence pass per component),
     * which is nothing on a small molecule and very visible on a big one: on a
     * dense few-thousand-atom system they were most of the cost of a
     * re-derivation that only ever read `.bonds`. Ask for them when you are
     * going to look at them.
     *
     * The formal charges are *not* gated on this: one charge per atom is the
     * valence that already had to be summed for the bond list, and `charges`
     * promises them unconditionally.
     */
    extras?: boolean;
}

export interface PerceptionDiagnostics {
    /**
     * Atoms whose valence is not one their element may take, counted in both
     * directions — a carbon at 3.5 is listed, and so is a nitrogen at 5. This
     * is a stricter test than a ceiling comparison, which would pass both a
     * half-unit and a nitrogen that is really a cation.
     */
    offValence: number[];
    /** Components held off their allowed valence that the solver was not run on. */
    repairable: number;
    /** Components the valence solver settled exactly. */
    solved: number;
    /** Components it settled only approximately. */
    approximate: number;
    /** Components it was offered but could not settle at all. */
    declined: number;
    /** Of those, how many were refused for exceeding the combination budget. */
    declinedByCombos: number;
    /** Of those, how many ran out of the work budget before finding an answer. */
    declinedByWork: number;
    /** Answers that were discarded because they were worse than what stood. */
    rolledBack: number;
    /** Bonds moved back onto the order their measured length argues for. */
    restoredBonds: number;
    /** Bonds moved by the final valence pass. */
    refinedBonds: number;
    /**
     * Whether the local relaxation was skipped because the structure is above
     * `MAX_RELAX_ATOMS`. The pass used to return silently, so a 30 000-atom
     * trajectory was drawn with over-coordinated centres and nothing said why.
     */
    relaxSkipped: boolean;
    /** Components the solver was not offered because of `maxSolverAtoms`. */
    oversizeComponents: number;
    /** Whether the total charge was supplied rather than assumed neutral. */
    chargeKnown: boolean;
}

export interface PerceptionResult {
    bonds: Bond[];
    /**
     * Formal charge per atom, implied by the final bond orders. Null for an
     * element with no entry in the tables. This is derived whatever the
     * options say — it costs nothing, changes no bond, and is the one piece of
     * the valence machinery that is pure addition to what the viewer already
     * displayed.
     */
    charges: (number | null)[];
    diagnostics: PerceptionDiagnostics;
}

/**
 * The atom list the passes below see: ghosts are read as the dummy element.
 *
 * A ghost has no valence to satisfy, so leaving it with the element its basis
 * belongs to would make it a permanent valence violation — a ghost carbon with
 * no bonds is a carbon four bonds short — and the solver would spend its budget
 * on a component that cannot be fixed. As `X` it has no target valence at all,
 * which is exactly true. The copy is made only when the structure has a ghost.
 */
function bpGhostView(rawAtoms: Atom[]): Atom[] {
    let has = false;
    for (let i = 0; i < rawAtoms.length; i++) {
        const a = rawAtoms[i];
        if (a && a.ghost === true && a.element !== 'X') { has = true; break; }
    }
    if (!has) return rawAtoms;
    const out: Atom[] = new Array(rawAtoms.length);
    for (let i = 0; i < rawAtoms.length; i++) {
        const a = rawAtoms[i];
        out[i] = (a && a.ghost === true) ? { ...a, element: 'X' } : a;
    }
    return out;
}

/**
 * Perceive bonds, orders and formal charges for a bare coordinate set.
 *
 * The geometric check runs first — the element-pair length table, the order
 * snap, and the local relaxation — and is what runs by default; its result
 * stands untouched unless the valence check is asked for. That check is
 * offered only the connected components that still hold an atom off its
 * allowed valence, and only a solution that passes every check replaces what
 * the geometry produced.
 *
 * The formal charges and the diagnostics are produced either way: they read
 * the finished orders and write none, and they are where the result can be
 * inspected — how many atoms the geometric check left off their allowed
 * valence, and how many components the valence check could settle if asked.
 *
 * A ghost atom is not part of any of it: it is left out of the connectivity
 * pass, so it ends up with no bond and every other atom's bonds are exactly
 * what they would be without it (see `bpGhostView`).
 */
export function perceiveBonds(rawAtoms: Atom[], opts?: PerceptionOptions): PerceptionResult {
    const atoms = bpGhostView(rawAtoms);
    const options = opts || {};
    const n = atoms.length;
    const diagnostics: PerceptionDiagnostics = {
        offValence: [],
        repairable: 0,
        solved: 0,
        approximate: 0,
        declined: 0,
        declinedByCombos: 0,
        declinedByWork: 0,
        rolledBack: 0,
        restoredBonds: 0,
        refinedBonds: 0,
        relaxSkipped: false,
        oversizeComponents: 0,
        chargeKnown: typeof options.charge === 'number'
    };
    if (n === 0) return { bonds: [], charges: [], diagnostics: diagnostics };

    const maxSolver = options.maxSolverAtoms === undefined ? MAX_SOLVER_ATOMS : options.maxSolverAtoms;
    const totalCharge = typeof options.charge === 'number' ? options.charge : 0;

    // 1) Connectivity and the geometric order for every bond.
    const bondMap = bpConnectivity(atoms);

    // 2) Drop any order the element pair has no entry for.
    bpFixOrders(bondMap, atoms);

    // 3) Bring over-coordinated centres as close to their lengths as the
    //    lengths allow. This is the pass the viewer has always run. It is
    //    skipped above MAX_RELAX_ATOMS, and the skip is reported rather than
    //    silent: a caller that sees over-coordinated centres in a big structure
    //    can now tell whether the pass declined or found nothing to do.
    if (n <= MAX_RELAX_ATOMS) bpRelaxValence(bondMap, atoms);
    else diagnostics.relaxSkipped = true;

    // 4) Hand each component still off its allowed valence to the solver, when
    //    it has been asked for. Without the solver the geometric result stands
    //    exactly as the viewer has always produced it, and the components the
    //    solver *could* have settled are counted instead — so a caller can see
    //    what is available before opting in.
    const useSolver = options.solve === true;
    const wantExtras = options.extras === true;
    // The component split serves the valence check and the diagnostics. A
    // caller who asks for neither is not made to walk the whole graph for it.
    const comps = (useSolver || wantExtras) ? bpComponents(bondMap, n) : [];

    if (!useSolver) {
        for (let c = 0; c < comps.length; c++) {
            if (comps[c].length > maxSolver) {
                diagnostics.oversizeComponents++;
                continue;
            }
            if (bpHasViolation(bondMap, atoms, comps[c])) diagnostics.repairable++;
        }
    }

    // One index map for every component, cleared by the solve as it finishes
    // with it (see bpSolveComponent).
    const localOf: number[] = useSolver ? new Array(n).fill(-1) : [];

    for (let c = 0; c < comps.length && useSolver; c++) {
        const comp = comps[c];
        if (comp.length > maxSolver) {
            diagnostics.oversizeComponents++;
            continue;
        }
        if (!bpHasViolation(bondMap, atoms, comp)) continue;

        // The total charge describes the whole structure, not each fragment.
        // Giving +1 to every component of a two-fragment +1 system asked each
        // one to be a cation; only a component that covers every atom can
        // inherit the stated total, and a fragment of a larger structure is
        // solved as neutral (its own charge is not stated anywhere).
        const compCharge = comp.length === n ? totalCharge : 0;
        const attempt = bpSolveComponent(bondMap, atoms, comp, compCharge, localOf);
        if (!attempt.solution) {
            diagnostics.declined++;
            if (attempt.reason === 'combo-overflow') diagnostics.declinedByCombos++;
            else if (attempt.reason === 'work-budget') diagnostics.declinedByWork++;
            continue;
        }
        const sol = attempt.solution;

        // Safety valve. The solver's answer is adopted only when it leaves the
        // component no worse on the measure this pass exists to improve — the
        // total distance of its atoms' valences from the values they may take.
        // An exact answer lands every atom on its target, so it can only come
        // out at zero; an approximate one (fits every ceiling without
        // saturating anything) has no such guarantee, and this is what keeps
        // such an answer from being taken when the geometric pass was already
        // better. The property this buys is the one the whole pass is for:
        // no re-derivation can come out worse than the pass it replaces.
        const gapBefore = bpComponentGap(bondMap, atoms, comp);
        const snapshot: Array<[number, number, number]> = [];
        for (let c = 0; c < comp.length; c++) {
            const nb = bondMap.get(comp[c]);
            if (!nb) continue;
            nb.forEach(function (o, j) {
                if (j > comp[c]) snapshot.push([comp[c], j, o]);
            });
        }

        for (let e = 0; e < sol.orders.length; e++) {
            const a = sol.orders[e][0];
            const b = sol.orders[e][1];
            const na = bondMap.get(a);
            const nb = bondMap.get(b);
            if (na) na.set(b, sol.orders[e][2]);
            if (nb) nb.set(a, sol.orders[e][2]);
        }

        if (bpComponentGap(bondMap, atoms, comp) > gapBefore + 1e-9) {
            for (let e = 0; e < snapshot.length; e++) {
                const na = bondMap.get(snapshot[e][0]);
                const nb = bondMap.get(snapshot[e][1]);
                if (na) na.set(snapshot[e][1], snapshot[e][2]);
                if (nb) nb.set(snapshot[e][0], snapshot[e][2]);
            }
            diagnostics.rolledBack++;
            continue;
        }

        if (sol.exact) diagnostics.solved++;
        else diagnostics.approximate++;
    }

    // 5) Give a ring its aromatic bonds back where the solver could only
    //    express it as an alternating single/double pattern. Only the solver
    //    costs a component its 1.5 bonds, so only the solver's run needs this.
    if (useSolver) {
        for (let c = 0; c < comps.length; c++) {
            diagnostics.restoredBonds += bpRestoreAromatic(bondMap, atoms, comps[c]);
        }
    }

    // 5b) One last look at valence. A centre can still be off a value its
    //     element may take — a half unit left where a table entry brackets two
    //     orders, or a component the valence check declined — and a bond around
    //     such a centre can sometimes be moved to fix it without the
    //     neighbourhood paying for it elsewhere. Only strictly-improving moves
    //     are kept, so this cannot leave a structure worse than the passes
    //     above did; and hydrogen bonds are pinned, so it cannot churn the
    //     largest group of bonds in the structure for nothing.
    diagnostics.refinedBonds = bpRefineValences(bondMap, atoms);

    // 5c) The ring rule, once more, as the guarantee: a pass above may have
    //     re-planned a centre onto 1.5 (the local relaxation does that when a
    //     centre is over-coordinated, and the aromatic restore of the solver
    //     path looks at the length alone), and an order that states a ring must
    //     not survive on a bond that is not in one. In the ordinary path this
    //     changes nothing — the connectivity pass already applied it — so it
    //     costs a second O(V+E) walk and no behaviour.
    bpAromaticRing(bondMap, atoms);

    // 6) Collect the result. The formal charges and the off-valence list are
    //    derived from the finished orders here — they read, they never write,
    //    so they carry no risk to the bonds the passes above settled.
    const bonds: Bond[] = [];
    const seen: { [key: string]: number } = {};
    const off = wantExtras ? new Set<number>() : null;
    // One entry per atom, always. The result type says the charges are derived
    // whatever the options say, and a caller that read .charges got an empty
    // array whenever it had not also asked for the diagnostics — the one field
    // of the two that costs nothing (the valence is already summed for the bond
    // list) was the one held back.
    const charges: (number | null)[] = new Array(n).fill(null);

    bondMap.forEach(function (neighbors, i) {
        let valence = 0;
        neighbors.forEach(function (order, j) {
            valence += order;
            const key = i < j ? i + '-' + j : j + '-' + i;
            if (seen[key]) return;
            seen[key] = 1;
            bonds.push({ atom1: i < j ? i : j, atom2: i < j ? j : i, order: order });
        });
        const z = bpZOf(atoms[i].element);
        if (off && bpValenceGap(z, valence, bpCeiling(atoms[i].element)) > 1e-6) {
            (off as Set<number>).add(i);
        }
        charges[i] = bpAtomicCharge(z, valence, neighbors.size);
    });

    bonds.sort(function (a, b) {
        return a.atom1 - b.atom1 || a.atom2 - b.atom2;
    });
    if (off) {
        (off as Set<number>).forEach(function (i) {
            diagnostics.offValence.push(i);
        });
        diagnostics.offValence.sort(function (a, b) {
            return a - b;
        });
    }

    return {
        bonds: bonds,
        charges: charges,
        diagnostics: diagnostics
    };
}

// ---------------------------------------------------------------------------
// 10. One source for both sides
// ---------------------------------------------------------------------------

/**
 * The functions that make up the module, in dependency order. The webview gets
 * them by name, so this list is the contract: a function used anywhere above
 * has to appear here, or the client copy would reference something undefined.
 */
const MODULE_FUNCTIONS: Array<(...args: never[]) => unknown> = [
    bpNorm, bpPairKey, bpSpecKey, bpDist, bpIsTerminal, bpCeiling, bpZOf,
    bpIsMetal, bpMetalSpec, bpSpecList, bpValenceModelled, bpIsGhost,
    bpNonAromaticOrder, bpAromaticRing, bpPairCutoff, bpGeometricOrder, bpConnectivity,
    bpAllowedOrders, bpNearestAllowed, bpFixOrders,
    bpMaximumMatching,
    bpPossibleValences, bpValenceOk, bpValenceGap, bpPruneValences, bpAtomicCharge,
    bpCheckValency, bpCheckSaturation, bpCheckCharge, bpDetermineOrders,
    bpRelaxLadder, bpRelaxCenter, bpRelaxValence,
    bpValenceOf, bpComponentGap, bpRestoreAromatic,
    bpSupportedOrder, bpSpecDev, bpRefineValences,
    bpComponents, bpHasViolation, bpSolveComponent,
    bpGhostView,
    perceiveBonds
];

/** Every module-level constant the emitted source declares. */
const MODULE_CONSTANTS: string[] = [
    'COV_R', 'PAIR_CUT', 'SPECS', 'CEIL', 'TERMINAL', 'Z_OF', 'METALS', 'VALENCE',
    'FALLBACK_VALENCE', 'OUTER_E', 'MAX_RELAX_ATOMS', 'MAX_SOLVER_ATOMS',
    'MAX_COMBOS', 'MAX_COMBO_WORK', 'MAX_SHED_BONDS', 'MIN_BOND_DISTANCE'
];

/**
 * A function body with its comments removed. The guard below reads names, and
 * prose is not a reference: the word `XOR` in a comment once looked like a
 * module constant to a naive scan.
 */
function bpWithoutComments(body: string): string {
    return body
        .replace(/\/\*[\s\S]*?\*\//g, '')
        .replace(/(^|[^:])\/\/[^\n]*/g, '$1');
}

/** The names a body declares for itself: locals and parameters. */
function bpLocalNames(code: string, body: string): { [name: string]: boolean } {
    const local: { [name: string]: boolean } = {};
    const decl = /(?:const|let|var)\s+([A-Za-z_$][\w$]*)/g;
    let m: RegExpExecArray | null = decl.exec(code);
    while (m !== null) {
        local[m[1]] = true;
        m = decl.exec(code);
    }
    const sig = /\(([^)]*)\)/.exec(body);
    if (sig !== null) {
        const params = sig[1].split(',');
        for (let i = 0; i < params.length; i++) {
            const name = params[i].trim().replace(/=.*$/, '').trim();
            if (name !== '') local[name] = true;
        }
    }
    return local;
}

/**
 * The module as one self-contained expression, for the webview to evaluate.
 *
 * The tables are emitted as JSON and the functions as their own source text,
 * inside a single scope, so the names the functions reference resolve to the
 * tables declared alongside them. Nothing here is hand-written a second time:
 * both sides run this text, which is what keeps the host detector and the
 * viewport's re-detection from drifting apart.
 */
export const BOND_PERCEPTION_SOURCE: string = (function buildSource(): string {
    const parts: string[] = [];
    parts.push('(function(){');
    parts.push('var COV_R=' + JSON.stringify(COV_R) + ';');
    parts.push('var PAIR_CUT=' + JSON.stringify(PAIR_CUT) + ';');
    parts.push('var SPECS=' + JSON.stringify(SPECS) + ';');
    parts.push('var CEIL=' + JSON.stringify(CEIL) + ';');
    parts.push('var TERMINAL=' + JSON.stringify(TERMINAL) + ';');
    parts.push('var Z_OF=' + JSON.stringify(Z_OF) + ';');
    parts.push('var METALS=' + JSON.stringify(METALS) + ';');
    parts.push('var VALENCE=' + JSON.stringify(VALENCE) + ';');
    parts.push('var FALLBACK_VALENCE=' + JSON.stringify(FALLBACK_VALENCE) + ';');
    parts.push('var OUTER_E=' + JSON.stringify(OUTER_E) + ';');
    parts.push('var MAX_RELAX_ATOMS=' + MAX_RELAX_ATOMS + ';');
    parts.push('var MAX_SOLVER_ATOMS=' + MAX_SOLVER_ATOMS + ';');
    parts.push('var MAX_COMBOS=' + MAX_COMBOS + ';');
    parts.push('var MAX_COMBO_WORK=' + MAX_COMBO_WORK + ';');
    parts.push('var MAX_SHED_BONDS=' + MAX_SHED_BONDS + ';');
    parts.push('var MIN_BOND_DISTANCE=' + MIN_BOND_DISTANCE + ';');
    // Every name the emitted functions reference has to be declared in the
    // same scope: the constants emitted above or one of the functions in
    // MODULE_FUNCTIONS. A missing entry would become a free variable in the
    // page, where it surfaces only when a panel opens or a bond is re-detected,
    // so the contract is checked here instead, at module load. Comments are
    // stripped first, and names a body declares itself are not references.
    const declared: { [name: string]: boolean } = {};
    for (let i = 0; i < MODULE_CONSTANTS.length; i++) declared[MODULE_CONSTANTS[i]] = true;
    for (let i = 0; i < MODULE_FUNCTIONS.length; i++) {
        const name = MODULE_FUNCTIONS[i].name;
        if (name === '') throw new Error('bond perception source contains an unnamed function');
        declared[name] = true;
    }
    for (let i = 0; i < MODULE_FUNCTIONS.length; i++) {
        const body = MODULE_FUNCTIONS[i].toString();
        const code = bpWithoutComments(body);
        const local = bpLocalNames(code, body);
        const caps = code.match(/\b[A-Z][A-Z0-9_]{2,}\b/g) || [];
        for (let r = 0; r < caps.length; r++) {
            if (!declared[caps[r]] && !local[caps[r]]) {
                throw new Error('bond perception source refers to a name the page does not define: ' + caps[r]);
            }
        }
        const calls = code.match(/\bbp[A-Z]\w*\b/g) || [];
        for (let r = 0; r < calls.length; r++) {
            if (!declared[calls[r]] && !local[calls[r]]) {
                throw new Error('bond perception source refers to a function outside MODULE_FUNCTIONS: ' + calls[r]);
            }
        }
        parts.push(body);
    }
    parts.push('return {perceiveBonds:perceiveBonds, pairCutoff:bpPairCutoff, supportedOrder:bpSupportedOrder};');
    parts.push('})()');
    const source = parts.join('\n');
    // The webview receives this inside its own template literal, so a backtick
    // or a dollar-brace anywhere in the module — including in a comment — would
    // silently truncate the script it is embedded in. Failing loudly here is
    // the only way that mistake cannot reach a user.
    if (source.indexOf('`') >= 0 || source.indexOf('${') >= 0) {
        throw new Error('bond perception source contains a template-literal delimiter');
    }
    // The page puts this text inside a <script> element. The HTML parser ends
    // that element at the first `</script`, even inside a JavaScript string or
    // comment, so the rest of the module would be handed to the markup parser.
    if (/<\/script/i.test(source)) {
        throw new Error('bond perception source contains a </script> sequence');
    }
    return source;
})();

/**
 * Internals exposed for the test harness only.
 *
 * These are the pieces whose termination and correctness cannot be observed
 * through `perceiveBonds` alone — the matching in particular, which is only
 * reachable when a component both needs the solver and fits its budget. They
 * are deliberately NOT part of `BOND_PERCEPTION_SOURCE`, so nothing here
 * reaches the webview or costs it a byte.
 */
export const __internals = {
    matching: bpMaximumMatching,
    pruneValences: bpPruneValences,
    determineOrders: bpDetermineOrders,
    connectivity: bpConnectivity,
    components: bpComponents,
    restAromatic: bpRestoreAromatic,
    relaxValence: bpRelaxValence,
    valenceGap: bpValenceGap,
    aromaticRing: bpAromaticRing,
    atomicCharge: bpAtomicCharge,
    specDev: bpSpecDev,
    specList: bpSpecList
};
