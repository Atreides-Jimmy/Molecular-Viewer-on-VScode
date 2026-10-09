/**
 * Molecular mechanics geometry optimisation ("MM Optimize")
 * ========================================================
 * An edit places atoms where the user asked them to be, not where the chemistry
 * wants them; this module settles the difference. It is a small, self-contained
 * force field — bond stretching, angle bending, torsion, out-of-plane
 * (improper) and a Lennard-Jones non-bonded term — whose reference lengths and
 * angles are the plugin's own experience tables (`empiricalBondLength` /
 * `empiricalBondAngle` / `empiricalVdwRadius`), the same windows the hydrogen
 * completion reads and the same numbers the SMILES -> 3D engine projects its
 * generated structures onto, so a relaxed bond is the bond generation would have
 * produced rather than a second, drifting set of numbers.
 *
 * What is shared with that engine is the numbers, not the method, and the two
 * are deliberately kept apart: the engine's relaxation is a *constraint
 * projection* (every bond length and 1-3 distance pushed onto its table value,
 * rings and sp2 centres flattened by projection, non-bonded contacts pushed
 * apart) that exists to reproduce the reference generator's coordinates
 * exactly; this module is an *energy minimiser* that moves a named selection
 * against frozen neighbours and can only ever lower the energy, which is what
 * lets a hand-edited structure be settled without disturbing the rest. Neither
 * algorithm can stand in for the other — a projection has no notion of frozen
 * atoms or of an energy that must not rise, an energy minimiser cannot promise
 * the reference's coordinates — so each caller reaches for the one it needs,
 * and both reach for the same tables.
 *
 * Only the atoms the caller names as mobile may move; every other atom is a
 * fixed source of force. That is what makes the operation usable as a *partial*
 * optimisation — a fragment can be relaxed against the rest of a structure
 * without disturbing it — and it is also why the energy reported back is the
 * energy of the mobile subsystem's own terms: every term that involves at least
 * one mobile atom is summed, a term between two fixed atoms is omitted because
 * it cannot change and would only add a constant.
 *
 * The optimiser is a preconditioned conjugate-gradient descent. Each mobile
 * atom's gradient is scaled by the inverse of its own stiffness (the bond and
 * angle force constants around it), the step is capped, and a step is accepted
 * only when it lowers the energy — so a run is monotone by construction and the
 * geometry it returns is never worse than the one it was given, converged or
 * not. The direction is the Polak-Ribiere conjugate one on top of that
 * preconditioned gradient (falling back to the plain gradient whenever the
 * conjugate direction is not a descent one), the prefactor follows the line
 * search — it grows while steps keep being accepted and is halved until one is
 * — and a run stops when the largest force falls under the tolerance, when the
 * energy stops improving, when the step budget or the wall-clock budget is
 * spent, or when no step can be found.
 *
 * Deliberate restraints, all of them there to keep a bad or unusual input from
 * being "improved" into a worse one:
 *
 *   - an angle whose current value is more than 30 degrees off its theoretical
 *     target carries no force at all, so a strained ring (cyclopropane's 60
 *     degree internal angle) is left exactly as drawn instead of being pulled
 *     open against the bonds that hold it;
 *   - a torsion is only set up where four non-hydrogen atoms define it, so a
 *     methyl group's rotation stays free (the angle and non-bonded terms
 *     already settle it) and only a chain's backbone and a planar group's
 *     conjugation are actually restrained;
 *   - the non-bonded term is switched off smoothly between 6 and 8 Å (a bare
 *     cutoff would leave a step of its energy at the boundary, and a pair
 *     crossing it would stall a line search), excludes bonded and 1-3 pairs,
 *     scales a 1-4 pair to half, and is continued by its tangent below 0.2 σ
 *     so that a pile-up of atoms stays finite while still being pushed apart;
 *     a 1-3 pair that comes inside 0.45 of its van der Waals sum is re-admitted
 *     as pure repulsion, because an angle squeezed that far sits in the angle
 *     term's capped window and two neighbours would otherwise be free to
 *     collapse onto each other;
 *   - the whole call is bounded by atom count, step count and wall-clock time.
 *
 * Contract: a structured result, never throws, never mutates its input, no
 * randomness — one structure always yields one answer.
 */

import {
    empiricalBondAngle,
    empiricalBondLength,
    empiricalVdwRadius,
} from './smilesTo3d';

/** Name the panel shows next to the parameter list. */
const FF_NAME = 'MVFF (built-in)';

/** Above this many mobile atoms a call is refused instead of run. */
const MAX_MOBILE_ATOMS = 3000;
const DEFAULT_MAX_ITERATIONS = 500;
const DEFAULT_TIME_BUDGET_MS = 5000;
/** Convergence: largest per-atom force component, kcal/(mol·Å). */
const DEFAULT_GRADIENT_TOL = 0.05;
/** Consecutive accepted steps gaining less than this are a plateau (kcal/mol). */
const ENERGY_PLATEAU = 1e-7;
const PLATEAU_STEPS = 6;
/** Largest distance one atom may travel in one accepted step (Å). */
const MAX_STEP = 0.15;
/** Step prefactor: starts here, grows on accepted steps, halves on rejected. */
const LAMBDA_START = 1.0;
const LAMBDA_MAX = 4.0;
const LAMBDA_MIN = 1e-10;
const MAX_BACKTRACK = 40;
/** Angle term: a *compressed* angle beyond this carries no force, so a strained
 *  ring (cyclopropane's 60° internal angle) is left as drawn. */
const ANGLE_CAP_COMPRESSED_RAD = (30.0 * Math.PI) / 180.0;
/** A *widened* angle beyond this carries no force either; the window is wider
 *  on this side because an angle that is too open is not a clash — it is a
 *  drawing that should be pulled back — while an over-compressed one only ever
 *  appears in a strained ring that must be left alone. */
const ANGLE_CAP_OPEN_RAD = (60.0 * Math.PI) / 180.0;
/** Non-bonded term. */
const VDW_CUTOFF = 8.0;
/** The taper starts here and the pair is gone at VDW_CUTOFF (see evaluate). */
const VDW_TAPER_START = 6.0;
/** The pair list is built to here, so it stays valid between rebuilds. */
const VDW_LIST_CUTOFF = 10.0;
/** Pair list is rebuilt once any atom has drifted this far (Å) since the last
 *  build. Shorter than the 2 Å buffer the list is built with, so a pair that
 *  comes into range is always seen before it matters. */
const VDW_REBUILD_DRIFT = 0.5;
/** A 1-3 pair (two atoms sharing a neighbour) counts at full strength below
 *  this fraction of its van der Waals sum — see buildNonbondedList. */
const CLASH_1_3_FRACTION = 0.45;
/** A 1-4 pair (three bonds apart) counts at this fraction of full strength. */
const FOURTEEN_SCALE = 0.5;
/** σ/r ceiling: keeps a pile-up finite instead of Infinity/NaN. */
const MIN_DISTANCE_RATIO = 0.2;
/** Grid cell of the non-bonded search (one cell covers the whole list cutoff). */
const GRID_CELL = VDW_LIST_CUTOFF;

/** Preconditioner floor: an atom with no terms at all is not stepped blindly. */
const STIFFNESS_FLOOR = 30.0;

/** Bond stretch: k = CHARGE_K · Z_i · Z_j / r0^3, clamped to this band. */
const CHARGE_K = 664.12;
const BOND_K_MIN = 50.0;
const BOND_K_MAX = 2000.0;

/** UFF-style effective charges (the Z in the bond force constant). */
const EFFECTIVE_CHARGE: { [el: string]: number } = {
    H: 0.712, Li: 1.026, Be: 1.565, B: 1.755, C: 1.912, N: 2.544, O: 2.300,
    F: 1.735, Na: 1.081, Mg: 1.787, Al: 1.792, Si: 2.323, P: 2.863, S: 2.703,
    Cl: 2.348, K: 1.165, Ca: 1.986, Se: 2.516, Te: 2.470, As: 2.517,
    Br: 2.519, I: 2.650,
};
const CHARGE_DEFAULT = 2.0;

/** Angle bend force constants by the central element (kcal/mol/rad²). */
const ANGLE_FORCE: { [el: string]: number } = {
    H: 30.0, B: 80.0, C: 100.0, N: 100.0, O: 120.0, F: 90.0, Si: 70.0,
    P: 95.0, S: 90.0, Cl: 95.0, Se: 90.0, Te: 90.0, As: 95.0, Br: 95.0, I: 95.0,
};
const ANGLE_FORCE_DEFAULT = 90.0;
/** Angle stiffness rises with the centre's hybridisation (sp3, sp2, sp). */
const HYBRIDISATION_ANGLE_SCALE = [1.0, 1.35, 1.6];

/** Three-fold torsion barriers by element (kcal/mol); pair value is the mean. */
const TORSION_BARRIER: { [el: string]: number } = {
    C: 1.0, N: 0.5, O: 0.4, S: 0.8, P: 0.8, Si: 0.8,
};
const TORSION_DEFAULT = 0.5;
/** A single bond touching an sp2 centre prefers planarity (conjugation). */
const CONJUGATED_TORSION_BARRIER = 1.0;
/** A double bond's planarity restraint; cis and trans are both minima. */
const DOUBLE_TORSION_BARRIER = 8.0;
/** An aromatic bond's ring-planarity restraint. */
const AROMATIC_TORSION_BARRIER = 1.5;

/** Out-of-plane restraint at an sp2 centre (kcal/mol/Å²-equivalent). */
const IMPROPER_FORCE: { [el: string]: number } = {
    B: 120.0, C: 200.0, N: 150.0, O: 120.0,
};
const IMPROPER_DEFAULT = 100.0;
/** Fixed reference area the out-of-plane volume measure is divided by (Å²). */
const IMPROPER_AREA_REFERENCE = 2.0;

/** Lennard-Jones well depths (kcal/mol) paired with the van der Waals radii. */
const LJ_EPSILON: { [el: string]: number } = {
    H: 0.044, B: 0.083, C: 0.095, N: 0.077, O: 0.060, F: 0.053, Si: 0.113,
    P: 0.160, S: 0.174, Cl: 0.148, Se: 0.170, Te: 0.200, As: 0.160,
    Br: 0.168, I: 0.210,
};
const LJ_DEFAULT = 0.10;

/** One mobile atom's answer, in the order the caller named the atoms. */
export interface MmMovedAtom {
    index: number;
    x: number;
    y: number;
    z: number;
}

/** How many of each term involved at least one mobile atom. */
export interface MmTermCounts {
    bonds: number;
    angles: number;
    torsions: number;
    impropers: number;
    nonbonded: number;
}

export interface MmOptimizeResult {
    ok: boolean;
    /** One-line summary for the status bar. */
    message: string;
    /** The structure the answer was computed from; the caller checks it. */
    atomCount: number;
    /** The atoms that were allowed to move, in the caller's order. */
    indices: number[];
    /** Final coordinates of `indices`, in the same order. */
    coords: number[][];
    /** Energy of the mobile subsystem's terms, before and after (kcal/mol). */
    energyBefore: number;
    energyAfter: number;
    /** Accepted geometry steps. */
    iterations: number;
    /** Largest / RMS per-atom force at the final geometry, kcal/(mol·Å). */
    maxGrad: number;
    rmsGrad: number;
    /** Largest distance any atom travelled (Å). */
    maxMove: number;
    converged: boolean;
    /** Why the run stopped: 'gradient' | 'energy' | 'iterations' | 'time' | 'step'. */
    stopReason: string;
    counts: MmTermCounts;
    /** Selected entries that were not usable (bad index, no coordinates). */
    invalid: number;
}

/**
 * Settings a caller may change for one run. Every field is optional and every
 * value is validated and clamped here, so a panel can offer choices without
 * being able to put the optimiser into an unreasonable state.
 */
export interface MmOptimizeOptions {
    /** Accepted geometry steps before the run stops (1–100000). */
    maxIterations?: number;
    /** Wall-clock budget in ms (10–60000). */
    timeBudgetMs?: number;
    /** Convergence: largest per-atom force component, kcal/(mol·Å). */
    gradientTol?: number;
    /** Largest distance one atom may travel in one accepted step (Å). */
    maxStep?: number;
    /** Include the Lennard-Jones non-bonded term (default true). */
    nonbonded?: boolean;
    /** Include torsion terms (default true). */
    torsions?: boolean;
    /** Include the out-of-plane restraint at sp2 centres (default true). */
    impropers?: boolean;
}

/** The three term switches, resolved to their values for one run. */
export interface MmTermFlags {
    nonbonded: boolean;
    torsions: boolean;
    impropers: boolean;
}

/** One term of the force field, with everything the evaluation needs. */
interface BondTerm {
    i: number;
    j: number;
    k: number;
    r0: number;
}
interface AngleTerm {
    c: number;
    i: number;
    j: number;
    k: number;
    theta0: number;
}
interface TorsionTerm {
    i: number;
    j: number;
    k: number;
    l: number;
    n: number;
    phase: number;
    v: number;
}
interface ImproperTerm {
    c: number;
    i: number;
    j: number;
    k: number;
    kf: number;
}
interface NonbondedTerm {
    i: number;
    j: number;
    epsilon: number;
    sigma: number;
    scale: number;
    /** A 1-3 pair re-admitted at clash range: r^-12 only, never attractive. */
    repulsionOnly?: boolean;
}

/** Everything the force field is made of, built once per call. */
interface Field {
    bonds: BondTerm[];
    angles: AngleTerm[];
    torsions: TorsionTerm[];
    impropers: ImproperTerm[];
    nonbonded: NonbondedTerm[];
    counts: MmTermCounts;
    /** Per-atom stiffness (preconditioner input). */
    stiffness: Float64Array;
}

/** Normalised element symbol: first letter upper, rest lower. */
function sym(raw: unknown): string {
    if (typeof raw !== 'string') {
        return '';
    }
    const s = raw.trim();
    if (s.length === 0) {
        return '';
    }
    return s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
}

function isFinite3(x: unknown, y: unknown, z: unknown): boolean {
    return typeof x === 'number' && typeof y === 'number' && typeof z === 'number' &&
        isFinite(x) && isFinite(y) && isFinite(z);
}

/** Bond order snapped to the values the force field parameterises. */
function classifyOrder(order: unknown): number {
    const o = typeof order === 'number' && isFinite(order) ? order : 1;
    if (o < 1.25) {
        return 1;
    }
    if (o < 1.75) {
        return 1.5;
    }
    if (o < 2.5) {
        return 2;
    }
    return 3;
}

/** The bond-length table's order argument: 0 means aromatic. */
function tableOrder(order: number): number {
    return order === 1.5 ? 0 : order;
}

function effectiveCharge(el: string): number {
    const v = EFFECTIVE_CHARGE[el];
    return v === undefined ? CHARGE_DEFAULT : v;
}

function angleForce(el: string): number {
    const v = ANGLE_FORCE[el];
    return v === undefined ? ANGLE_FORCE_DEFAULT : v;
}

function torsionBarrier(el: string): number {
    const v = TORSION_BARRIER[el];
    return v === undefined ? TORSION_DEFAULT : v;
}

function improperForce(el: string): number {
    const v = IMPROPER_FORCE[el];
    return v === undefined ? IMPROPER_DEFAULT : v;
}

function ljEpsilon(el: string): number {
    const v = LJ_EPSILON[el];
    return v === undefined ? LJ_DEFAULT : v;
}

/**
 * Hybridisation code the engine uses: 2 = sp, 1 = sp2, 0 = sp3. Read off the
 * bond orders an atom already has, exactly the way the engine's own passes read
 * a bond list (a triple bond is sp, a double or an aromatic bond is sp2).
 */
function hybridisationOf(orders: number[]): number {
    let max = 0;
    for (let i = 0; i < orders.length; i++) {
        if (orders[i] > max) {
            max = orders[i];
        }
    }
    if (max >= 2.5) {
        return 2;
    }
    if (max >= 1.25) {
        return 1;
    }
    return 0;
}

/** Index of the first neighbour of `atom` that is not `skip` (or -1). */
function firstOtherNeighbour(nb: number[][], atom: number, skip: number): number {
    const list = nb[atom];
    if (!list) {
        return -1;
    }
    for (let i = 0; i < list.length; i++) {
        if (list[i] !== skip) {
            return list[i];
        }
    }
    return -1;
}

/**
 * Build every term that involves at least one mobile atom. Terms between two
 * fixed atoms are skipped: they cannot change, and leaving them out keeps the
 * reported energy a difference-preserving quantity. A fixed atom is otherwise a
 * full participant — its bonds, angles, torsions, impropers and non-bonded
 * contacts pull on every mobile atom it shares a term with, which is what makes
 * this a *restrained* optimisation of the selection against its surroundings
 * rather than an optimisation of the selection on its own.
 */
function buildField(
    elements: string[],
    x: Float64Array,
    y: Float64Array,
    z: Float64Array,
    bonds: { i: number; j: number; order: number }[],
    mobile: boolean[],
    n: number,
    flags: MmTermFlags
): Field {
    const bondsOut: BondTerm[] = [];
    const anglesOut: AngleTerm[] = [];
    const torsionsOut: TorsionTerm[] = [];
    const impropersOut: ImproperTerm[] = [];
    const stiffness = new Float64Array(n);

    const nb: number[][] = [];
    const nbOrders: number[][] = [];
    for (let i = 0; i < n; i++) {
        nb.push([]);
        nbOrders.push([]);
    }
    for (let b = 0; b < bonds.length; b++) {
        const bond = bonds[b];
        nb[bond.i].push(bond.j);
        nb[bond.j].push(bond.i);
        nbOrders[bond.i].push(bond.order);
        nbOrders[bond.j].push(bond.order);
    }

    // ---- bond stretch ----
    for (let b = 0; b < bonds.length; b++) {
        const bond = bonds[b];
        if (!mobile[bond.i] && !mobile[bond.j]) {
            continue;
        }
        const r0 = empiricalBondLength(elements[bond.i], elements[bond.j], tableOrder(bond.order));
        const k = CHARGE_K * effectiveCharge(elements[bond.i]) * effectiveCharge(elements[bond.j]) /
            (r0 * r0 * r0);
        const kClamped = Math.min(BOND_K_MAX, Math.max(BOND_K_MIN, k));
        bondsOut.push({ i: bond.i, j: bond.j, k: kClamped, r0: r0 });
        stiffness[bond.i] += kClamped;
        stiffness[bond.j] += kClamped;
    }

    // ---- angle bend ----
    for (let c = 0; c < n; c++) {
        const list = nb[c];
        if (list.length < 2) {
            continue;
        }
        const hyb = hybridisationOf(nbOrders[c]);
        const theta0 = (empiricalBondAngle(elements[c], hyb) * Math.PI) / 180.0;
        const k = angleForce(elements[c]) * HYBRIDISATION_ANGLE_SCALE[hyb];
        for (let a = 0; a < list.length - 1; a++) {
            for (let b2 = a + 1; b2 < list.length; b2++) {
                const i = list[a];
                const j = list[b2];
                if (!mobile[c] && !mobile[i] && !mobile[j]) {
                    continue;
                }
                anglesOut.push({ c: c, i: i, j: j, k: k, theta0: theta0 });
                stiffness[c] += 0.5 * k;
                stiffness[i] += 0.25 * k;
                stiffness[j] += 0.25 * k;
            }
        }
    }

    // ---- torsion ----
    // Only four non-hydrogen atoms define one, and a triple bond has none.
    const isH: boolean[] = [];
    for (let i = 0; i < n; i++) {
        isH.push(elements[i] === 'H');
    }
    for (let b = 0; flags.torsions && b < bonds.length; b++) {
        const bond = bonds[b];
        const j = bond.i;
        const k = bond.j;
        if (bond.order >= 2.5 || isH[j] || isH[k]) {
            continue;
        }
        const na = firstOtherNeighbour(nb, j, k);
        const nc = firstOtherNeighbour(nb, k, j);
        if (na < 0 || nc < 0 || isH[na] || isH[nc]) {
            continue;
        }
        // A dihedral's energy depends on all four of its atoms, so a term is
        // live as soon as ANY of them is mobile — including a mobile flanking
        // atom hanging off a fixed central pair, which is exactly the restrained
        // case (a terminal group being relaxed on a frozen backbone).
        if (!mobile[na] && !mobile[j] && !mobile[k] && !mobile[nc]) {
            continue;
        }
        const hybJ = hybridisationOf(nbOrders[j]);
        const hybK = hybridisationOf(nbOrders[k]);
        let terms: TorsionTerm | null = null;
        if (bond.order >= 1.75) {
            // A double bond: two-fold, so both cis and trans are minima — the
            // term can only straighten the group, never flip it inside out.
            terms = { i: na, j: j, k: k, l: nc, n: 2, phase: 0, v: DOUBLE_TORSION_BARRIER };
        } else if (bond.order >= 1.25) {
            // An aromatic bond: ring planarity, at a barrier below a real double
            // bond's because the ring's own angles already hold most of it.
            terms = { i: na, j: j, k: k, l: nc, n: 2, phase: 0, v: AROMATIC_TORSION_BARRIER };
        } else if (hybJ === 2 || hybK === 2) {
            // A single bond to a linear centre: effectively free rotation.
            continue;
        } else if (hybJ === 1 || hybK === 1) {
            // A single bond touching an sp2 centre: conjugated, prefers planar
            // (this is what keeps a vinyl or a biphenyl group coplanar).
            terms = { i: na, j: j, k: k, l: nc, n: 2, phase: 0, v: CONJUGATED_TORSION_BARRIER };
        } else {
            const v = 0.5 * (torsionBarrier(elements[j]) + torsionBarrier(elements[k]));
            terms = { i: na, j: j, k: k, l: nc, n: 3, phase: Math.PI / 3.0, v: v };
        }
        torsionsOut.push(terms);
        const add = terms.n * terms.n * terms.v * 0.25;
        stiffness[terms.j] += add;
        stiffness[terms.k] += add;
    }

    // ---- out-of-plane (improper) at an sp2 centre with three neighbours ----
    for (let c = 0; flags.impropers && c < n; c++) {
        const list = nb[c];
        if (list.length !== 3) {
            continue;
        }
        if (hybridisationOf(nbOrders[c]) !== 1) {
            continue;
        }
        const i = list[0];
        const j = list[1];
        const k = list[2];
        if (!mobile[c] && !mobile[i] && !mobile[j] && !mobile[k]) {
            continue;
        }
        const kf = improperForce(elements[c]);
        impropersOut.push({ c: c, i: i, j: j, k: k, kf: kf });
        stiffness[c] += 0.5 * kf;
    }

    // ---- non-bonded pair list (mobile against anything within the cutoff) ----
    const nonbondedOut = flags.nonbonded
        ? buildNonbondedList(x, y, z, elements, nb, mobile, n)
        : [];

    return {
        bonds: bondsOut,
        angles: anglesOut,
        torsions: torsionsOut,
        impropers: impropersOut,
        nonbonded: nonbondedOut,
        counts: {
            bonds: bondsOut.length,
            angles: anglesOut.length,
            torsions: torsionsOut.length,
            impropers: impropersOut.length,
            nonbonded: nonbondedOut.length,
        },
        stiffness: stiffness,
    };
}

/** Topological relation of a pair: 0 = bonded, 1 = 1-3, 2 = 1-4, 3 = further. */
function pairRelation(nb: number[][], bonded: Set<number>[], i: number, j: number): number {
    if (bonded[i].has(j)) {
        return 0;
    }
    const ni = nb[i];
    for (let a = 0; a < ni.length; a++) {
        const m = ni[a];
        if (bonded[m].has(j)) {
            return 1;
        }
        const nm = nb[m];
        for (let b = 0; b < nm.length; b++) {
            if (bonded[nm[b]].has(j)) {
                return 2;
            }
        }
    }
    return 3;
}

/**
 * Every non-bonded pair within VDW_LIST_CUTOFF that involves a mobile atom,
 * with the 1-4 pair scale already folded into its strength. Bonded and 1-3
 * pairs are skipped (their chemistry is carried by the bond and angle terms).
 */
function buildNonbondedList(
    x: Float64Array,
    y: Float64Array,
    z: Float64Array,
    elements: string[],
    nb: number[][],
    mobile: boolean[],
    n: number
): NonbondedTerm[] {
    const bonded: Set<number>[] = [];
    for (let i = 0; i < n; i++) {
        const s = new Set<number>();
        const list = nb[i];
        for (let k = 0; k < list.length; k++) {
            s.add(list[k]);
        }
        bonded.push(s);
    }
    // Cell keying: the three integer cell coordinates are joined into a string
    // key. A packed or hashed key is not usable here — two different cells
    // mapping to one key would make the same bucket be scanned twice (each pair
    // counted twice) and would change the answer the moment an atom crossed a
    // cell boundary, so the pair list must not depend on a hash at all.
    const cells = new Map<string, number[]>();
    const cellOf = new Int32Array(n * 3);
    for (let i = 0; i < n; i++) {
        const cx = Math.floor(x[i] / GRID_CELL);
        const cy = Math.floor(y[i] / GRID_CELL);
        const cz = Math.floor(z[i] / GRID_CELL);
        cellOf[i * 3] = cx;
        cellOf[i * 3 + 1] = cy;
        cellOf[i * 3 + 2] = cz;
        const key = cx + ',' + cy + ',' + cz;
        const bucket = cells.get(key);
        if (bucket) {
            bucket.push(i);
        } else {
            cells.set(key, [i]);
        }
    }
    const out: NonbondedTerm[] = [];
    const cut2 = VDW_LIST_CUTOFF * VDW_LIST_CUTOFF;
    for (let i = 0; i < n; i++) {
        if (!mobile[i]) {
            continue;
        }
        const cx = cellOf[i * 3];
        const cy = cellOf[i * 3 + 1];
        const cz = cellOf[i * 3 + 2];
        for (let dx = -1; dx <= 1; dx++) {
            for (let dy = -1; dy <= 1; dy++) {
                for (let dz = -1; dz <= 1; dz++) {
                    const key = (cx + dx) + ',' + (cy + dy) + ',' + (cz + dz);
                    const bucket = cells.get(key);
                    if (!bucket) {
                        continue;
                    }
                    for (let q = 0; q < bucket.length; q++) {
                        const j = bucket[q];
                        if (j === i) {
                            continue;
                        }
                        // Each pair is added by exactly one of its two mobile
                        // ends (the lower index), because the list is scanned
                        // from every mobile atom and a mobile-mobile pair would
                        // otherwise be counted twice.
                        if (mobile[j] && j < i) {
                            continue;
                        }
                        const rx = x[i] - x[j];
                        const ry = y[i] - y[j];
                        const rz = z[i] - z[j];
                        if (rx * rx + ry * ry + rz * rz > cut2) {
                            continue;
                        }
                        const rel = pairRelation(nb, bonded, i, j);
                        if (rel === 0) {
                            // Bonded: the bond term owns this pair.
                            continue;
                        }
                        const epsilon = Math.sqrt(ljEpsilon(elements[i]) * ljEpsilon(elements[j]));
                        const sigma = empiricalVdwRadius(elements[i]) + empiricalVdwRadius(elements[j]);
                        if (rel === 1) {
                            // A 1-3 pair is carried by its angle term and normally
                            // excluded — but not at clash range: an angle squeezed
                            // past the point where no exclusion can be physical
                            // (a hydrogen pushed onto its neighbour's hydrogen by a
                            // drag, say) would otherwise sit in the angle term's
                            // capped, force-free window and could collapse onto the
                            // atom. Below this fraction of the van der Waals sum the
                            // pair comes back as pure repulsion, which can only ever
                            // push the two apart.
                            const r = Math.sqrt(rx * rx + ry * ry + rz * rz);
                            if (!(r < CLASH_1_3_FRACTION * sigma)) {
                                continue;
                            }
                            out.push({ i: i, j: j, epsilon: epsilon, sigma: sigma, scale: 1.0, repulsionOnly: true });
                            continue;
                        }
                        out.push({
                            i: i,
                            j: j,
                            epsilon: epsilon,
                            sigma: sigma,
                            scale: rel === 2 ? FOURTEEN_SCALE : 1.0,
                        });
                    }
                }
            }
        }
    }
    return out;
}

/** Scratch buffers of one evaluation, allocated once per call. */
interface Gradient {
    gx: Float64Array;
    gy: Float64Array;
    gz: Float64Array;
}

/**
 * Energy of every term that involves a mobile atom, and the gradient of that
 * energy over the mobile atoms (other entries are left at zero). The two are
 * always evaluated together: the optimiser needs both, and evaluating them in
 * one pass guarantees they describe the same geometry.
 */
function evaluate(
    field: Field,
    x: Float64Array,
    y: Float64Array,
    z: Float64Array,
    mobile: boolean[],
    grad: Gradient | null
): number {
    const gx = grad ? grad.gx : null;
    const gy = grad ? grad.gy : null;
    const gz = grad ? grad.gz : null;
    // The buffers are reused across steps, so each evaluation starts from zero
    // (the accumulators below are the only writers).
    if (gx && gy && gz) {
        gx.fill(0);
        gy.fill(0);
        gz.fill(0);
    }
    const addGrad = (atom: number, vx: number, vy: number, vz: number): void => {
        if (gx === null || gy === null || gz === null || !mobile[atom]) {
            return;
        }
        gx[atom] += vx;
        gy[atom] += vy;
        gz[atom] += vz;
    };

    let energy = 0.0;

    // ---- bond stretch ----
    for (let t = 0; t < field.bonds.length; t++) {
        const term = field.bonds[t];
        const dx = x[term.i] - x[term.j];
        const dy = y[term.i] - y[term.j];
        const dz = z[term.i] - z[term.j];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        // Two atoms on top of each other leave the force direction undefined;
        // a fixed sub-nanometre nudge along +x makes the answer deterministic
        // and still pushes the pair apart.
        if (r < 1e-6) {
            const dev0 = 1e-6 - term.r0;
            energy += 0.5 * term.k * dev0 * dev0;
            addGrad(term.i, term.k * dev0, 0, 0);
            addGrad(term.j, -term.k * dev0, 0, 0);
            continue;
        }
        const dev = r - term.r0;
        energy += 0.5 * term.k * dev * dev;
        const f = (term.k * dev) / r;
        addGrad(term.i, dx * f, dy * f, dz * f);
        addGrad(term.j, -dx * f, -dy * f, -dz * f);
    }

    // ---- angle bend ----
    for (let t = 0; t < field.angles.length; t++) {
        const term = field.angles[t];
        let ux = x[term.i] - x[term.c];
        let uy = y[term.i] - y[term.c];
        let uz = z[term.i] - z[term.c];
        let vx = x[term.j] - x[term.c];
        let vy = y[term.j] - y[term.c];
        let vz = z[term.j] - z[term.c];
        const lu = Math.sqrt(ux * ux + uy * uy + uz * uz);
        const lv = Math.sqrt(vx * vx + vy * vy + vz * vz);
        if (lu < 1e-9 || lv < 1e-9) {
            continue;
        }
        ux /= lu; uy /= lu; uz /= lu;
        vx /= lv; vy /= lv; vz /= lv;
        let cosTheta = ux * vx + uy * vy + uz * vz;
        if (cosTheta > 1.0) {
            cosTheta = 1.0;
        } else if (cosTheta < -1.0) {
            cosTheta = -1.0;
        }
        const theta = Math.acos(cosTheta);
        const dev = theta - term.theta0;
        if (dev > ANGLE_CAP_OPEN_RAD) {
            // Off by more than the cap on the open side: a drawing that is too
            // wide to be trusted, so the term keeps its capped energy but exerts
            // no force and the optimiser never fights whatever holds it there.
            energy += 0.5 * term.k * ANGLE_CAP_OPEN_RAD * ANGLE_CAP_OPEN_RAD;
            continue;
        }
        if (dev < -ANGLE_CAP_COMPRESSED_RAD) {
            // Squeezed past the compressed cap: a strained ring, left as drawn.
            energy += 0.5 * term.k * ANGLE_CAP_COMPRESSED_RAD * ANGLE_CAP_COMPRESSED_RAD;
            continue;
        }
        energy += 0.5 * term.k * dev * dev;
        const sinTheta = Math.sqrt(1.0 - cosTheta * cosTheta);
        if (sinTheta < 1e-6) {
            continue;
        }
        // dθ/dx_i = −(v − cosθ·u) / (|r_i|·sinθ): moving i changes only u, and
        // d(u·v)/du = v.
        const f = (term.k * dev) / sinTheta;
        const giX = (-(vx - cosTheta * ux) * f) / lu;
        const giY = (-(vy - cosTheta * uy) * f) / lu;
        const giZ = (-(vz - cosTheta * uz) * f) / lu;
        const gjX = (-(ux - cosTheta * vx) * f) / lv;
        const gjY = (-(uy - cosTheta * vy) * f) / lv;
        const gjZ = (-(uz - cosTheta * vz) * f) / lv;
        addGrad(term.i, giX, giY, giZ);
        addGrad(term.j, gjX, gjY, gjZ);
        addGrad(term.c, -(giX + gjX), -(giY + gjY), -(giZ + gjZ));
    }

    // ---- torsion ----
    for (let t = 0; t < field.torsions.length; t++) {
        const term = field.torsions[t];
        const b1x = x[term.j] - x[term.i];
        const b1y = y[term.j] - y[term.i];
        const b1z = z[term.j] - z[term.i];
        const b2x = x[term.k] - x[term.j];
        const b2y = y[term.k] - y[term.j];
        const b2z = z[term.k] - z[term.j];
        const b3x = x[term.l] - x[term.k];
        const b3y = y[term.l] - y[term.k];
        const b3z = z[term.l] - z[term.k];
        const n1x = b1y * b2z - b1z * b2y;
        const n1y = b1z * b2x - b1x * b2z;
        const n1z = b1x * b2y - b1y * b2x;
        const n2x = b2y * b3z - b2z * b3y;
        const n2y = b2z * b3x - b2x * b3z;
        const n2z = b2x * b3y - b2y * b3x;
        const l2 = Math.sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
        const ln1 = Math.sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
        const ln2 = Math.sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
        if (l2 < 1e-9 || ln1 < 1e-9 || ln2 < 1e-9) {
            continue;
        }
        // The signed angle between the two plane normals about b2: its tangent
        // component is (n1×n2)·b̂2 and its cosine n1·n2. Both normals are
        // perpendicular to b2, so this atan2 is the dihedral itself, and the
        // gradient below is the exact derivative of *this* expression.
        const crossX = n1y * n2z - n1z * n2y;
        const crossY = n1z * n2x - n1x * n2z;
        const crossZ = n1x * n2y - n1y * n2x;
        const xDot = n1x * n2x + n1y * n2y + n1z * n2z;
        const yDot = (crossX * b2x + crossY * b2y + crossZ * b2z) / l2;
        const phi = Math.atan2(yDot, xDot);
        const arg = term.n * (phi - term.phase);
        energy += 0.5 * term.v * (1.0 - Math.cos(arg));
        if (!gx) {
            continue;
        }
        // d(phi)/d(atom), derived from the two plane normals: phi is the signed
        // angle about b2 between n1 and n2, so its derivative w.r.t. a normal is
        // a rotation vector in the plane (u1 below points the direction phi
        // decreases in), and the atoms' derivatives follow by the chain rule
        // through b1 = r_j − r_i, b2 = r_k − r_j, b3 = r_l − r_k. Translation
        // invariance is exact: the four derivatives sum to zero.
        const invLen = 1.0 / l2;
        const u1x = -(b2y * n1z - b2z * n1y) * invLen / (ln1 * ln1);
        const u1y = -(b2z * n1x - b2x * n1z) * invLen / (ln1 * ln1);
        const u1z = -(b2x * n1y - b2y * n1x) * invLen / (ln1 * ln1);
        const u2x = (b2y * n2z - b2z * n2y) * invLen / (ln2 * ln2);
        const u2y = (b2z * n2x - b2x * n2z) * invLen / (ln2 * ln2);
        const u2z = (b2x * n2y - b2y * n2x) * invLen / (ln2 * ln2);
        const gb1x = b2y * u1z - b2z * u1y;
        const gb1y = b2z * u1x - b2x * u1z;
        const gb1z = b2x * u1y - b2y * u1x;
        const gb3x = u2y * b2z - u2z * b2y;
        const gb3y = u2z * b2x - u2x * b2z;
        const gb3z = u2x * b2y - u2y * b2x;
        const gb2x = (u1y * b1z - u1z * b1y) + (b3y * u2z - b3z * u2y);
        const gb2y = (u1z * b1x - u1x * b1z) + (b3z * u2x - b3x * u2z);
        const gb2z = (u1x * b1y - u1y * b1x) + (b3x * u2y - b3y * u2x);
        const dEdPhi = 0.5 * term.v * term.n * Math.sin(arg);
        addGrad(term.i, -dEdPhi * gb1x, -dEdPhi * gb1y, -dEdPhi * gb1z);
        addGrad(term.j, dEdPhi * (gb1x - gb2x), dEdPhi * (gb1y - gb2y), dEdPhi * (gb1z - gb2z));
        addGrad(term.k, dEdPhi * (gb2x - gb3x), dEdPhi * (gb2y - gb3y), dEdPhi * (gb2z - gb3z));
        addGrad(term.l, dEdPhi * gb3x, dEdPhi * gb3y, dEdPhi * gb3z);
    }

    // ---- out-of-plane (improper) ----
    for (let t = 0; t < field.impropers.length; t++) {
        const term = field.impropers[t];
        const ux = x[term.i] - x[term.c];
        const uy = y[term.i] - y[term.c];
        const uz = z[term.i] - z[term.c];
        const vx = x[term.j] - x[term.c];
        const vy = y[term.j] - y[term.c];
        const vz = z[term.j] - z[term.c];
        const wx = x[term.k] - x[term.c];
        const wy = y[term.k] - y[term.c];
        const wz = z[term.k] - z[term.c];
        // Signed volume of the parallelepiped: zero exactly when planar.
        const v6 = ux * (vy * wz - vz * wy) +
            uy * (vz * wx - vx * wz) +
            uz * (vx * wy - vy * wx);
        // Scaled by a fixed reference area, so the energy of a given
        // out-of-plane height does not depend on the centre's own bond lengths
        // (which also keeps the gradient an exact derivative of this energy).
        const scaled = v6 / IMPROPER_AREA_REFERENCE;
        energy += 0.5 * term.kf * scaled * scaled;
        if (!gx) {
            continue;
        }
        const f = (term.kf * scaled) / IMPROPER_AREA_REFERENCE;
        // d(v6)/d(atom), exact for all four atoms: v6 is the scalar triple
        // product u·(v×w), so each derivative is the cross product of the other
        // two vectors and the centre takes the negative sum of the three.
        const giX = vy * wz - vz * wy;
        const giY = vz * wx - vx * wz;
        const giZ = vx * wy - vy * wx;
        const gjX = wy * uz - wz * uy;
        const gjY = wz * ux - wx * uz;
        const gjZ = wx * uy - wy * ux;
        const gkX = uy * vz - uz * vy;
        const gkY = uz * vx - ux * vz;
        const gkZ = ux * vy - uy * vx;
        addGrad(term.i, f * giX, f * giY, f * giZ);
        addGrad(term.j, f * gjX, f * gjY, f * gjZ);
        addGrad(term.k, f * gkX, f * gkY, f * gkZ);
        addGrad(term.c, -f * (giX + gjX + gkX), -f * (giY + gjY + gkY), -f * (giZ + gjZ + gkZ));
    }

    // ---- non-bonded ----
    for (let t = 0; t < field.nonbonded.length; t++) {
        const term = field.nonbonded[t];
        let dx = x[term.i] - x[term.j];
        let dy = y[term.i] - y[term.j];
        let dz = z[term.i] - z[term.j];
        let r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r >= VDW_CUTOFF) {
            continue;
        }
        // Coincident non-bonded atoms: same deterministic nudge as the bond
        // term, so a pile-up separates instead of sitting in a 0/0 direction.
        if (r < 1e-6) {
            dx = 1e-6;
            dy = 0;
            dz = 0;
            r = 1e-6;
        }
        let ratio = term.sigma / r;
        if (term.repulsionOnly) {
            // A 1-3 clash: pure r^-12 repulsion, kept finite (and still
            // separating) by the same tangent continuation the main path uses.
            const c = 1.0 / MIN_DISTANCE_RATIO;
            if (ratio > c) {
                const s6c = Math.pow(c, 6);
                const s12c = s6c * s6c;
                const rc = term.sigma / c;
                const slope = -(12.0 * s12c) / rc;
                energy += term.scale * term.epsilon * (s12c + slope * (r - rc));
                if (gx) {
                    const dedr = term.scale * term.epsilon * slope;
                    addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
                    addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
                }
                continue;
            }
            const s6r = Math.pow(ratio, 6);
            const s12r = s6r * s6r;
            energy += term.scale * term.epsilon * s12r;
            if (gx) {
                const dedr = -(12.0 * term.scale * term.epsilon * s12r) / r;
                addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
                addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
            }
            continue;
        }
        if (ratio > 1.0 / MIN_DISTANCE_RATIO) {
            // Closer than MIN_DISTANCE_RATIO·σ the potential is continued by the
            // tangent at that radius: the energy stays bounded (a pile-up cannot
            // reach Infinity or NaN) while its slope still pushes the pair
            // apart, which a hard clamp would not — a clamped LJ is flat, and a
            // flat region with a large force leaves a line search with nothing
            // to descend.
            const rc = term.sigma * MIN_DISTANCE_RATIO;
            const ratioC = 1.0 / MIN_DISTANCE_RATIO;
            const s6c = Math.pow(ratioC, 6);
            const ljc = s6c * s6c - 2.0 * s6c;
            const slope = (12.0 * (s6c - s6c * s6c)) / rc;
            const value = ljc + slope * (r - rc);
            energy += term.scale * term.epsilon * value;
            if (!gx) {
                continue;
            }
            const dedr = term.scale * term.epsilon * slope;
            addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
            addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
            continue;
        }
        const s6 = Math.pow(ratio, 6);
        const s12 = s6 * s6;
        // A switching function tapers the pair to zero between TAPER_START and
        // the cutoff instead of chopping it off: a bare cutoff leaves a step of
        // the LJ energy at 8 Å (≈6e-4 kcal/mol) and a pair crossing it makes the
        // energy discontinuous, which no line search can step through.
        let switchValue = 1.0;
        let dSwitch = 0.0;
        if (r > VDW_TAPER_START) {
            const rOn2 = VDW_TAPER_START * VDW_TAPER_START;
            const rOff2 = VDW_CUTOFF * VDW_CUTOFF;
            const denom = rOff2 - rOn2;
            const u = rOff2 - r * r;
            switchValue = (u * u * (rOff2 + 2.0 * r * r - 3.0 * rOn2)) /
                (denom * denom * denom);
            dSwitch = (12.0 * r * u * (rOn2 - r * r)) / (denom * denom * denom);
        }
        const lj = s12 - 2.0 * s6;
        energy += term.scale * term.epsilon * lj * switchValue;
        if (!gx) {
            continue;
        }
        // dE/dr of (LJ·S): the LJ part as before, plus the taper's own slope.
        const dedr = term.scale * term.epsilon *
            ((12.0 * (s6 - s12) * switchValue) / r + lj * dSwitch);
        addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
        addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
    }

    return energy;
}

/** Largest and RMS per-atom force over the mobile atoms. */
function forceStats(grad: Gradient, mobile: boolean[], n: number): { max: number; rms: number } {
    let max = 0.0;
    let sum = 0.0;
    let count = 0;
    for (let i = 0; i < n; i++) {
        if (!mobile[i]) {
            continue;
        }
        const fx = grad.gx[i];
        const fy = grad.gy[i];
        const fz = grad.gz[i];
        const f = Math.sqrt(fx * fx + fy * fy + fz * fz);
        if (f > max) {
            max = f;
        }
        sum += f * f;
        count++;
    }
    return { max: max, rms: count > 0 ? Math.sqrt(sum / count) : 0.0 };
}

/** Normalise the raw payloads into arrays the field builder can trust. */
function readInput(rawAtoms: unknown, rawBonds: unknown, rawSelected: unknown): {
    n: number;
    elements: string[];
    x: Float64Array;
    y: Float64Array;
    z: Float64Array;
    bonds: { i: number; j: number; order: number }[];
    selected: number[];
    invalid: number;
} {
    const atomList: unknown[] = Array.isArray(rawAtoms) ? rawAtoms : [];
    const n = atomList.length;
    const elements: string[] = new Array(n);
    const x = new Float64Array(n);
    const y = new Float64Array(n);
    const z = new Float64Array(n);
    const usable: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) {
        const a = atomList[i] as { element?: unknown; x?: unknown; y?: unknown; z?: unknown } | null;
        const el = a ? sym(a.element) : '';
        elements[i] = el;
        const finite = !!a && isFinite3(a.x, a.y, a.z);
        usable[i] = finite;
        x[i] = finite ? (a as { x: number }).x : 0.0;
        y[i] = finite ? (a as { y: number }).y : 0.0;
        z[i] = finite ? (a as { z: number }).z : 0.0;
    }
    const bonds: { i: number; j: number; order: number }[] = [];
    const rawBondList: unknown[] = Array.isArray(rawBonds) ? rawBonds : [];
    for (let b = 0; b < rawBondList.length; b++) {
        const bond = rawBondList[b] as { atom1?: unknown; atom2?: unknown; order?: unknown } | null;
        if (!bond) {
            continue;
        }
        const i = bond.atom1;
        const j = bond.atom2;
        if (typeof i !== 'number' || typeof j !== 'number') {
            continue;
        }
        if (!isFinite(i) || !isFinite(j)) {
            continue;
        }
        const ii = Math.trunc(i);
        const jj = Math.trunc(j);
        if (ii === jj || ii < 0 || jj < 0 || ii >= n || jj >= n) {
            continue;
        }
        if (!usable[ii] || !usable[jj]) {
            continue;
        }
        bonds.push({ i: ii, j: jj, order: classifyOrder(bond.order) });
    }
    const selected: number[] = [];
    let invalid = 0;
    const seen = new Set<number>();
    const rawList: unknown[] = Array.isArray(rawSelected) ? rawSelected : [];
    for (let s = 0; s < rawList.length; s++) {
        const v = rawList[s];
        if (typeof v !== 'number' || !isFinite(v)) {
            invalid++;
            continue;
        }
        const idx = Math.trunc(v);
        if (idx < 0 || idx >= n || !usable[idx] || seen.has(idx)) {
            invalid++;
            continue;
        }
        seen.add(idx);
        selected.push(idx);
    }
    return { n: n, elements: elements, x: x, y: y, z: z, bonds: bonds, selected: selected, invalid: invalid };
}

/** Result of a single-point evaluation (the force field's public access). */
export interface MmEvaluation {
    energy: number;
    gradient: number[][];
    maxGrad: number;
    rmsGrad: number;
}

/** Every option, resolved against its default and clamped to its sane range. */
interface MMResolvedOptions {
    maxIterations: number;
    timeBudgetMs: number;
    gradientTol: number;
    maxStep: number;
    flags: MmTermFlags;
}

/**
 * Validate a caller's settings. Anything missing, of the wrong type, or outside
 * its range falls back to the default, so an arbitrary message payload cannot
 * make a run unreasonable (and a panel's choices are the only thing that can
 * change the behaviour).
 */
function readOptions(rawOptions: unknown): MMResolvedOptions {
    const o = (rawOptions && typeof rawOptions === 'object' ? rawOptions : {}) as MmOptimizeOptions;
    const num = (v: unknown, fallback: number, lo: number, hi: number): number =>
        typeof v === 'number' && isFinite(v) ? Math.max(lo, Math.min(hi, v)) : fallback;
    const flag = (v: unknown, fallback: boolean): boolean => typeof v === 'boolean' ? v : fallback;
    return {
        maxIterations: Math.trunc(num(o.maxIterations, DEFAULT_MAX_ITERATIONS, 1, 100000)),
        timeBudgetMs: num(o.timeBudgetMs, DEFAULT_TIME_BUDGET_MS, 10, 60000),
        gradientTol: num(o.gradientTol, DEFAULT_GRADIENT_TOL, 1e-6, 1000),
        maxStep: num(o.maxStep, MAX_STEP, 0.005, 2.0),
        flags: {
            nonbonded: flag(o.nonbonded, true),
            torsions: flag(o.torsions, true),
            impropers: flag(o.impropers, true),
        },
    };
}

/**
 * Energy and per-atom forces at the given geometry, over every term that
 * involves at least one mobile atom. `rawSelected` names the mobile atoms; an
 * empty or missing selection means "everything is mobile", which is the whole
 * structure's energy. `rawOptions` may switch the term blocks off, exactly as
 * the optimiser does.
 */
export function evaluateForceField(
    rawAtoms: unknown,
    rawBonds: unknown,
    rawSelected?: unknown,
    rawOptions?: unknown
): MmEvaluation {
    const settings = readOptions(rawOptions);
    const input = readInput(rawAtoms, rawBonds, rawSelected);
    const n = input.n;
    const mobile: boolean[] = new Array(n);
    if (!Array.isArray(rawSelected)) {
        for (let i = 0; i < n; i++) {
            mobile[i] = true;
        }
    } else {
        for (let i = 0; i < n; i++) {
            mobile[i] = false;
        }
        for (let s = 0; s < input.selected.length; s++) {
            mobile[input.selected[s]] = true;
        }
    }
    const field = buildField(input.elements, input.x, input.y, input.z, input.bonds, mobile, n, settings.flags);
    const grad: Gradient = {
        gx: new Float64Array(n),
        gy: new Float64Array(n),
        gz: new Float64Array(n),
    };
    const energy = evaluate(field, input.x, input.y, input.z, mobile, grad);
    const stats = forceStats(grad, mobile, n);
    const gradient: number[][] = [];
    for (let i = 0; i < n; i++) {
        gradient.push([grad.gx[i], grad.gy[i], grad.gz[i]]);
    }
    return { energy: energy, gradient: gradient, maxGrad: stats.max, rmsGrad: stats.rms };
}

/**
 * Relax the selected atoms of a structure. Returns the new coordinates of the
 * selected atoms only — every unselected atom is left exactly where it was —
 * together with the energy before and after and the reason the run stopped.
 */
export function optimizeSelection(
    rawAtoms: unknown,
    rawBonds: unknown,
    rawSelected: unknown,
    rawOptions?: unknown
): MmOptimizeResult {
    const settings = readOptions(rawOptions);
    const maxIterations = settings.maxIterations;
    const timeBudgetMs = settings.timeBudgetMs;
    const gradientTol = settings.gradientTol;
    const maxStep = settings.maxStep;

    const input = readInput(rawAtoms, rawBonds, rawSelected);
    const n = input.n;
    const empty: MmOptimizeResult = {
        ok: false,
        message: '',
        atomCount: n,
        indices: [],
        coords: [],
        energyBefore: 0,
        energyAfter: 0,
        iterations: 0,
        maxGrad: 0,
        rmsGrad: 0,
        maxMove: 0,
        converged: false,
        stopReason: 'none',
        counts: { bonds: 0, angles: 0, torsions: 0, impropers: 0, nonbonded: 0 },
        invalid: input.invalid,
    };
    if (n === 0) {
        empty.message = 'Molecular Mechanics: the structure has no atoms';
        return empty;
    }
    if (input.selected.length === 0) {
        empty.message = 'Molecular Mechanics: select at least one atom to optimise';
        return empty;
    }
    if (input.selected.length > MAX_MOBILE_ATOMS) {
        empty.message = 'Molecular Mechanics: ' + input.selected.length +
            ' atoms selected — the limit is ' + MAX_MOBILE_ATOMS + ' atoms per run';
        return empty;
    }

    const mobile: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) {
        mobile[i] = false;
    }
    for (let s = 0; s < input.selected.length; s++) {
        mobile[input.selected[s]] = true;
    }

    const field = buildField(input.elements, input.x, input.y, input.z, input.bonds, mobile, n, settings.flags);
    if (field.bonds.length === 0 && field.angles.length === 0 &&
        field.torsions.length === 0 && field.impropers.length === 0 &&
        field.nonbonded.length === 0) {
        empty.ok = true;
        empty.converged = true;
        empty.stopReason = 'gradient';
        empty.message = 'Molecular Mechanics: nothing to optimise — the selection has no bonds or contacts';
        empty.indices = input.selected.slice();
        empty.coords = input.selected.map((i) => [input.x[i], input.y[i], input.z[i]]);
        empty.counts = field.counts;
        return empty;
    }

    // Working copies: only the mobile atoms are ever written back, but the
    // buffers are whole-structure so every term reads the same arrays.
    const x = new Float64Array(input.x);
    const y = new Float64Array(input.y);
    const z = new Float64Array(input.z);
    const trialX = new Float64Array(n);
    const trialY = new Float64Array(n);
    const trialZ = new Float64Array(n);
    const grad: Gradient = {
        gx: new Float64Array(n),
        gy: new Float64Array(n),
        gz: new Float64Array(n),
    };
    const precondition = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        precondition[i] = 1.0 / Math.max(field.stiffness[i], STIFFNESS_FLOOR);
    }
    // Conjugate-gradient bookkeeping: the search direction of the last accepted
    // step, the gradient it was built from, and the current direction. A plain
    // steepest descent converges far too slowly on a flat non-bonded landscape
    // (hundreds of steps to polish off the last tenth of a kcal), while the
    // monotone acceptance below still guarantees that every step lowers the
    // energy — the direction only decides how fast, never whether.
    const dirX = new Float64Array(n);
    const dirY = new Float64Array(n);
    const dirZ = new Float64Array(n);
    const prevDirX = new Float64Array(n);
    const prevDirY = new Float64Array(n);
    const prevDirZ = new Float64Array(n);
    const prevGradX = new Float64Array(n);
    const prevGradY = new Float64Array(n);
    const prevGradZ = new Float64Array(n);
    let cgValid = false;

    const startTime = Date.now();
    let lambda = LAMBDA_START;
    let energy = evaluate(field, x, y, z, mobile, grad);
    const energyBefore = energy;
    let lastEnergy = energy;
    let stats = forceStats(grad, mobile, n);
    let iterations = 0;
    let plateau = 0;
    let converged = false;
    let stopReason = 'gradient';
    let resetUsed = false;
    // Per-atom drift since the pair list was built: how far each atom has moved
    // in total, so the list is rebuilt before a pair can slip into range unseen.
    const drift = new Float64Array(n);

    while (true) {
        if (stats.max <= gradientTol) {
            converged = true;
            stopReason = 'gradient';
            break;
        }
        if (iterations >= maxIterations) {
            stopReason = 'iterations';
            break;
        }
        if (Date.now() - startTime > timeBudgetMs) {
            stopReason = 'time';
            break;
        }

        // Search direction: Polak-Ribiere conjugate gradient when a previous
        // direction exists, the plain gradient otherwise (and whenever the
        // conjugate direction is not a descent direction).
        let beta = 0.0;
        if (cgValid) {
            let num = 0.0;
            let den = 0.0;
            for (let i = 0; i < n; i++) {
                if (!mobile[i]) {
                    continue;
                }
                num += grad.gx[i] * (grad.gx[i] - prevGradX[i]) +
                    grad.gy[i] * (grad.gy[i] - prevGradY[i]) +
                    grad.gz[i] * (grad.gz[i] - prevGradZ[i]);
                den += prevGradX[i] * prevGradX[i] + prevGradY[i] * prevGradY[i] +
                    prevGradZ[i] * prevGradZ[i];
            }
            beta = den > 1e-12 ? Math.max(0.0, num / den) : 0.0;
            if (!(beta > 0.0)) {
                cgValid = false;
            }
        }
        let slope = 0.0;
        for (let i = 0; i < n; i++) {
            if (!mobile[i]) {
                continue;
            }
            dirX[i] = grad.gx[i] + (cgValid ? beta * prevDirX[i] : 0.0);
            dirY[i] = grad.gy[i] + (cgValid ? beta * prevDirY[i] : 0.0);
            dirZ[i] = grad.gz[i] + (cgValid ? beta * prevDirZ[i] : 0.0);
            slope += dirX[i] * grad.gx[i] + dirY[i] * grad.gy[i] + dirZ[i] * grad.gz[i];
        }
        if (!(slope > 0.0)) {
            // Not a descent direction (or the structure has no mobile terms at
            // all): fall back to the gradient itself.
            cgValid = false;
            for (let i = 0; i < n; i++) {
                dirX[i] = grad.gx[i];
                dirY[i] = grad.gy[i];
                dirZ[i] = grad.gz[i];
            }
        }

        // Trust radius: the *whole* preconditioned direction is scaled down so its
// largest atom step is the cap, rather than clipping each atom on its own.
// A per-atom clip would decouple the trial from the line-search factor — every
// attempt would move the clipped atom by the same distance, so halving the
// factor could not shrink the step, and a step whose curvature cost exceeds its
// gradient gain would keep failing until the factor fell below the clip
// threshold. Scaling the vector keeps the trial length proportional to the
// factor, which is what makes the search able to walk down at all.
        let scaledMax = 0.0;
        for (let i = 0; i < n; i++) {
            if (!mobile[i]) {
                continue;
            }
            const sx = precondition[i] * dirX[i];
            const sy = precondition[i] * dirY[i];
            const sz = precondition[i] * dirZ[i];
            const len = Math.sqrt(sx * sx + sy * sy + sz * sz);
            if (len > scaledMax) {
                scaledMax = len;
            }
        }
        const trustScale = scaledMax > maxStep ? maxStep / scaledMax : 1.0;

        // A trial step: the preconditioned descent direction inside that trust
        // radius, line-searched by halving until the energy actually drops.
        let accepted = false;
        let acceptedEnergy = energy;
        for (let attempt = 0; attempt < MAX_BACKTRACK; attempt++) {
            for (let i = 0; i < n; i++) {
                trialX[i] = x[i];
                trialY[i] = y[i];
                trialZ[i] = z[i];
            }
            let moved = 0.0;
            for (let i = 0; i < n; i++) {
                if (!mobile[i]) {
                    continue;
                }
                const dx = -lambda * trustScale * precondition[i] * dirX[i];
                const dy = -lambda * trustScale * precondition[i] * dirY[i];
                const dz = -lambda * trustScale * precondition[i] * dirZ[i];
                trialX[i] = x[i] + dx;
                trialY[i] = y[i] + dy;
                trialZ[i] = z[i] + dz;
                const step = Math.sqrt(dx * dx + dy * dy + dz * dz);
                if (step > moved) {
                    moved = step;
                }
            }
            if (moved < 1e-9) {
                break;
            }
            // The non-bonded pair list is NOT refreshed here: it describes the
            // geometry the search is descending from, and rebuilding it from a
            // trial geometry would change the energy function between attempts
            // (the clash guard's membership depends on distance), which leaves a
            // line search comparing incomparable numbers. It is refreshed after
            // an accepted step instead, once an atom has drifted far enough.
            const trialEnergy = evaluate(field, trialX, trialY, trialZ, mobile, null);
            if (trialEnergy < energy) {
                energy = trialEnergy;
                accepted = true;
                acceptedEnergy = trialEnergy;
                lambda = Math.min(LAMBDA_MAX, lambda * 1.6);
                break;
            }
            lambda *= 0.5;
            if (lambda < LAMBDA_MIN) {
                break;
            }
        }
        if (!accepted) {
            if (cgValid && !resetUsed) {
                // The conjugate direction led nowhere; drop it, restore a step
                // size that a plain gradient step can use, and try once more.
                resetUsed = true;
                cgValid = false;
                lambda = LAMBDA_START;
                continue;
            }
            stopReason = 'step';
            break;
        }

        // Commit the accepted step: keep the direction and the gradient it was
        // built from for the next step's conjugation, then re-evaluate.
        let driftMax = 0.0;
        for (let i = 0; i < n; i++) {
            if (!mobile[i]) {
                continue;
            }
            prevDirX[i] = dirX[i];
            prevDirY[i] = dirY[i];
            prevDirZ[i] = dirZ[i];
            prevGradX[i] = grad.gx[i];
            prevGradY[i] = grad.gy[i];
            prevGradZ[i] = grad.gz[i];
            const mx = trialX[i] - x[i];
            const my = trialY[i] - y[i];
            const mz = trialZ[i] - z[i];
            x[i] = trialX[i];
            y[i] = trialY[i];
            z[i] = trialZ[i];
            drift[i] += Math.sqrt(mx * mx + my * my + mz * mz);
            if (drift[i] > driftMax) {
                driftMax = drift[i];
            }
        }
        cgValid = true;
        resetUsed = false;
        // Rebuild the non-bonded pair list once an atom has drifted far enough
        // for a pair to have entered or left the range the list has to cover.
        // Built from the committed geometry, so every energy the next line
        // search compares comes from the same list.
        if (settings.flags.nonbonded && driftMax >= VDW_REBUILD_DRIFT) {
            field.nonbonded = buildNonbondedList(
                x, y, z, input.elements, buildNeighbourLists(input.bonds, n), mobile, n);
            field.counts.nonbonded = field.nonbonded.length;
            for (let i = 0; i < n; i++) {
                drift[i] = 0.0;
            }
        }
        // The gradient is evaluated on the committed geometry (it is a full pass,
        // so its energy is the accepted one) and that evaluation is what the
        // next line search compares against — never the pre-rebuild value. A
        // freshly built pair list can hold a clash the previous one could not
        // see (the guard's membership depends on distance), and comparing the
        // next trial against a stale basis would look like a step-independent
        // jump of the pair's energy that no line search can walk down.
        energy = evaluate(field, x, y, z, mobile, grad);
        stats = forceStats(grad, mobile, n);
        iterations++;

        // A converging run keeps gaining less and less; stop when the gain has
        // been flat for a while rather than burning the whole step budget.
        if (lastEnergy - energy < ENERGY_PLATEAU) {
            plateau++;
        } else {
            plateau = 0;
        }
        lastEnergy = energy;
        if (plateau >= PLATEAU_STEPS) {
            stopReason = 'energy';
            break;
        }
    }
    if (energy > energyBefore) {
        // Cannot happen — steps are only committed on a strict decrease — but
        // the returned geometry must never be a worse one, so this is enforced.
        energy = energyBefore;
        for (let i = 0; i < n; i++) {
            x[i] = input.x[i];
            y[i] = input.y[i];
            z[i] = input.z[i];
        }
    }

    const indices: number[] = [];
    const coords: number[][] = [];
    let maxMove = 0.0;
    for (let s = 0; s < input.selected.length; s++) {
        const i = input.selected[s];
        indices.push(i);
        coords.push([x[i], y[i], z[i]]);
        const mx = x[i] - input.x[i];
        const my = y[i] - input.y[i];
        const mz = z[i] - input.z[i];
        const d = Math.sqrt(mx * mx + my * my + mz * mz);
        if (d > maxMove) {
            maxMove = d;
        }
    }

    const stopText: { [key: string]: string } = {
        gradient: 'converged',
        energy: 'energy converged',
        iterations: 'step budget reached',
        time: 'time budget reached',
        step: 'no further gain',
    };
    const delta = energy - energyBefore;
    const message = 'Molecular Mechanics: ' + indices.length + ' atom(s) in ' + iterations +
        ' steps — E ' + energyBefore.toFixed(3) + ' → ' + energy.toFixed(3) +
        ' kcal/mol (Δ ' + delta.toFixed(3) + ', ' + (stopText[stopReason] || stopReason) + ')';

    return {
        ok: true,
        message: message,
        atomCount: n,
        indices: indices,
        coords: coords,
        energyBefore: energyBefore,
        energyAfter: energy,
        iterations: iterations,
        maxGrad: stats.max,
        rmsGrad: stats.rms,
        maxMove: maxMove,
        converged: converged,
        stopReason: stopReason,
        counts: field.counts,
        invalid: input.invalid,
    };
}

/** Neighbour lists of a bond list (the pair-list rebuild needs them). */
function buildNeighbourLists(bonds: { i: number; j: number; order: number }[], n: number): number[][] {
    const nb: number[][] = [];
    for (let i = 0; i < n; i++) {
        nb.push([]);
    }
    for (let b = 0; b < bonds.length; b++) {
        nb[bonds[b].i].push(bonds[b].j);
        nb[bonds[b].j].push(bonds[b].i);
    }
    return nb;
}

/** One choice of a settable parameter, with the value that travels with a run. */
export interface MmControlOption {
    value: number | boolean;
    label: string;
}

/** One settable parameter, rendered by the panel as a selection control. */
export interface MmControl {
    /** Key name in the run request (`MmOptimizeOptions`). */
    key: string;
    label: string;
    options: MmControlOption[];
    /** The value selected by default — the module's own default. */
    defaultValue: number | boolean;
}

/** What the panel shows before anything has been run. */
export interface MmFieldDescription {
    name: string;
    /** Parameter description lines, rendered in the panel exactly as they come. */
    lines: string[];
    /** The settings a caller may change, rendered as selection controls. */
    controls: MmControl[];
}

/**
 * The force field's parameter set in prose plus the settings a caller may
 * change, both generated from the constants above so the panel can never
 * describe a field the optimiser does not use and can never offer a default
 * that differs from the one the module would apply on its own.
 */
export function describeForceField(): MmFieldDescription {
    const bondK = (elA: number, elB: number, r0: number): number =>
        Math.min(BOND_K_MAX, Math.max(BOND_K_MIN, (CHARGE_K * elA * elB) / (r0 * r0 * r0)));
    const kLow = bondK(0.712, 1.912, 1.54);
    const kHigh = bondK(2.863, 2.863, 1.60);
    const kC = bondK(1.912, 1.912, 1.20);
    const angleLow = Math.min(ANGLE_FORCE.H, ANGLE_FORCE.Si) * HYBRIDISATION_ANGLE_SCALE[0];
    const angleHigh = Math.max(ANGLE_FORCE.O, ANGLE_FORCE.N) * HYBRIDISATION_ANGLE_SCALE[2];
    const dflt = ' (default)';
    return {
        name: FF_NAME,
        lines: [
            'Bond stretch k = ' + CHARGE_K.toFixed(2) + '·Z₁Z₂/r₀³ ≈ ' + Math.round(kLow) + '–' +
            Math.round(kHigh) + ' kcal/(mol·Å²) (' + Math.round(kC) + ' for C=C); angle bend towards the ' +
            'element\'s sp/sp²/sp³ target (k = ' + Math.round(angleLow) + '–' + Math.round(angleHigh) +
            ' kcal/(mol·rad²), no force beyond 30° off)',
            'Torsion 3-fold staggered / 2-fold planar (V = 0.4–8), improper at three-coordinate sp² ' +
            'centres, Lennard-Jones σ = the van der Waals radii, ε = 0.044–0.21 kcal/mol, switched off ' +
            'between ' + VDW_TAPER_START.toFixed(1) + ' and ' + VDW_CUTOFF.toFixed(1) + ' Å, no charges',
            'Restrained optimisation: only the selected atoms move, but every fixed atom still pulls on ' +
            'them through all of its terms (bonds, angles, torsions, impropers, non-bonded contacts).',
        ],
        controls: [
            {
                key: 'nonbonded', label: 'Non-bonded', defaultValue: true,
                options: [
                    { value: true, label: 'on · 1-4 x0.5 · 6-8 Å' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'torsions', label: 'Torsion terms', defaultValue: true,
                options: [
                    { value: true, label: 'on' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'impropers', label: 'Improper terms', defaultValue: true,
                options: [
                    { value: true, label: 'on' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'gradientTol', label: 'Convergence max|F|', defaultValue: DEFAULT_GRADIENT_TOL,
                options: [
                    { value: DEFAULT_GRADIENT_TOL, label: DEFAULT_GRADIENT_TOL.toFixed(2) + ' kcal/mol/Å' + dflt },
                    { value: 0.01, label: '0.01 kcal/mol/Å (tight)' },
                    { value: 0.2, label: '0.2 kcal/mol/Å (loose)' },
                ],
            },
            {
                key: 'maxIterations', label: 'Max steps', defaultValue: DEFAULT_MAX_ITERATIONS,
                options: [
                    { value: DEFAULT_MAX_ITERATIONS, label: String(DEFAULT_MAX_ITERATIONS) + dflt },
                    { value: 100, label: '100' },
                    { value: 1000, label: '1000' },
                ],
            },
            {
                key: 'timeBudgetMs', label: 'Time budget', defaultValue: DEFAULT_TIME_BUDGET_MS,
                options: [
                    { value: DEFAULT_TIME_BUDGET_MS, label: (DEFAULT_TIME_BUDGET_MS / 1000) + ' s' + dflt },
                    { value: 1000, label: '1 s' },
                    { value: 2000, label: '2 s' },
                    { value: 10000, label: '10 s' },
                ],
            },
            {
                key: 'maxStep', label: 'Max step', defaultValue: MAX_STEP,
                options: [
                    { value: MAX_STEP, label: MAX_STEP.toFixed(2) + ' Å' + dflt },
                    { value: 0.05, label: '0.05 Å (gentle)' },
                    { value: 0.3, label: '0.30 Å (fast)' },
                ],
            },
        ],
    };
}