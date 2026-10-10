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
 *   - the non-bonded term is switched off smoothly over the 2 Å inside the
 *     cutoff (6 → 8 Å at the default, settable 5–12 Å) (a bare
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
/** Non-bonded term: the pair is switched off here. Settable per run (the panel
 *  offers 5-12 Å); the default is the value every earlier release used. */
const DEFAULT_CUTOFF = 8.0;
/** The taper starts this far inside the cutoff, so a pair reaches zero exactly
 *  at the cutoff. Fixed width, so every cutoff keeps the same smooth shoulder. */
const DEFAULT_TAPER_WIDTH = 2.0;
/** The pair list is built this far past the cutoff, so it stays valid between
 *  rebuilds (with the drift threshold below as the guarantee). */
const LIST_BUFFER = 2.0;
/** Pair list is rebuilt once any atom has drifted this far (Å) since the last
 *  build. Shorter than the buffer the list is built with, so a pair that
 *  comes into range is always seen before it matters. */
const VDW_REBUILD_DRIFT = 0.5;
/** A 1-3 pair (two atoms sharing a neighbour) counts at full strength below
 *  this fraction of its van der Waals sum — see buildNonbondedList. */
const CLASH_1_3_FRACTION = 0.45;
/** A 1-4 pair (three bonds apart) counts at this fraction of full strength. */
const FOURTEEN_SCALE = 0.5;
/** σ/r ceiling: keeps a pile-up finite instead of Infinity/NaN. */
const MIN_DISTANCE_RATIO = 0.2;
/** One history point is kept per accepted step up to this many; beyond it the
 *  array is halved so a long run's curve stays bounded (see MmSession.pushPoint). */
const HISTORY_CAP = 4000;

/** Preconditioner floor: an atom with no terms at all is not stepped blindly. */
const STIFFNESS_FLOOR = 30.0;

/** Bond stretch: k = CHARGE_K · Z_i · Z_j / r0^3, clamped to this band. */
const CHARGE_K = 664.12;
const BOND_K_MIN = 50.0;
const BOND_K_MAX = 2000.0;

/**
 * UFF-style effective charges (the Z in the bond force constant). The elements
 * below the first rows — the main-group elements of periods 4 and 5 that a
 * structure may name (Ga Ge Rb Sr In Sn Sb) — are added by the module's own
 * rule: a new element takes the value of its *nearest tabulated member in the
 * same group* (the periodic row above), so Ga/In follow Al, Ge/Sn follow Si, Sb
 * follows As and Rb/Sr follow K/Ca. These are group-analogue estimates, with no
 * literature value claimed for them; they exist so the heavier main-group
 * elements get the same treatment as the first rows (smilesTo3d already carries
 * their covalent and van der Waals radii) instead of falling through to the
 * table default. No already-tabulated value is ever changed.
 */
const EFFECTIVE_CHARGE: { [el: string]: number } = {
    H: 0.712, Li: 1.026, Be: 1.565, B: 1.755, C: 1.912, N: 2.544, O: 2.300,
    F: 1.735, Na: 1.081, Mg: 1.787, Al: 1.792, Si: 2.323, P: 2.863, S: 2.703,
    Cl: 2.348, K: 1.165, Ca: 1.986, Se: 2.516, Te: 2.470, As: 2.517,
    Br: 2.519, I: 2.650,
    Ga: 1.792, Ge: 2.323, Rb: 1.165, Sr: 1.986, In: 1.792, Sn: 2.323, Sb: 2.517,
};
const CHARGE_DEFAULT = 2.0;

/**
 * Angle bend force constants by the central element (kcal/mol/rad²). The
 * period-4/5 main-group centres follow their group analogue above, the same rule
 * the effective charges use: Al/Ga/In follow B, Ge/Sn follow Si, Sb follows As.
 * Group-analogue estimates, no literature value claimed; existing values kept.
 */
const ANGLE_FORCE: { [el: string]: number } = {
    H: 30.0, B: 80.0, C: 100.0, N: 100.0, O: 120.0, F: 90.0, Si: 70.0,
    P: 95.0, S: 90.0, Cl: 95.0, Se: 90.0, Te: 90.0, As: 95.0, Br: 95.0, I: 95.0,
    Al: 80.0, Ga: 80.0, In: 80.0, Ge: 70.0, Sn: 70.0, Sb: 95.0,
};
const ANGLE_FORCE_DEFAULT = 90.0;
/** Angle stiffness rises with the centre's hybridisation (sp3, sp2, sp). */
const HYBRIDISATION_ANGLE_SCALE = [1.0, 1.35, 1.6];

/**
 * Three-fold torsion barriers by element (kcal/mol); pair value is the mean.
 * Period-4/5 additions follow their group analogue (Ge/Sn follow Si, As/Sb
 * follow P, Se/Te follow S) — group-analogue estimates, no literature value.
 */
const TORSION_BARRIER: { [el: string]: number } = {
    C: 1.0, N: 0.5, O: 0.4, S: 0.8, P: 0.8, Si: 0.8,
    Ge: 0.8, Sn: 0.8, As: 0.8, Sb: 0.8, Se: 0.8, Te: 0.8,
};
const TORSION_DEFAULT = 0.5;
/** A single bond touching an sp2 centre prefers planarity (conjugation). */
const CONJUGATED_TORSION_BARRIER = 1.0;
/** A double bond's planarity restraint; cis and trans are both minima. */
const DOUBLE_TORSION_BARRIER = 8.0;
/** An aromatic bond's ring-planarity restraint. */
const AROMATIC_TORSION_BARRIER = 1.5;

/**
 * Out-of-plane restraint at an sp2 centre (kcal/mol/Å²-equivalent). Period-4/5
 * additions follow their group analogue (Al/Ga/In follow B, Ge/Sn follow C,
 * As/Sb follow N, Se/Te follow O) — group-analogue estimates, no literature.
 */
const IMPROPER_FORCE: { [el: string]: number } = {
    B: 120.0, C: 200.0, N: 150.0, O: 120.0,
    Al: 120.0, Ga: 120.0, In: 120.0, Ge: 200.0, Sn: 200.0,
    As: 150.0, Sb: 150.0, Se: 120.0, Te: 120.0,
};
const IMPROPER_DEFAULT = 100.0;
/** Fixed reference area the out-of-plane volume measure is divided by (Å²). */
const IMPROPER_AREA_REFERENCE = 2.0;

/**
 * Lennard-Jones well depths (kcal/mol) paired with the van der Waals radii (the
 * σ comes from smilesTo3d, which already tabulates every element named here).
 * Period-4/5 additions follow their group analogue (Al/Ga/In follow B, Ge/Sn
 * follow Si, Sb follows As) — group-analogue estimates, no literature value.
 */
const LJ_EPSILON: { [el: string]: number } = {
    H: 0.044, B: 0.083, C: 0.095, N: 0.077, O: 0.060, F: 0.053, Si: 0.113,
    P: 0.160, S: 0.174, Cl: 0.148, Se: 0.170, Te: 0.200, As: 0.160,
    Br: 0.168, I: 0.210,
    Al: 0.083, Ga: 0.083, In: 0.083, Ge: 0.113, Sn: 0.113, Sb: 0.160,
};
const LJ_DEFAULT = 0.10;

/**
 * ---- Electrostatics (a light Coulomb term with a distance-dependent
 * dielectric) ----
 * A structure's partial charges are not supplied by the caller, so they are
 * derived here from the elements and bond orders alone: every bond shifts a
 * little charge from the less electronegative atom to the more electronegative
 * one, weighted by the Pauling electronegativity difference and, lightly, by the
 * bond order, and the whole set is then shifted so the structure is neutral. No
 * formal charges are needed — a bondless ion pair simply stays neutral. The
 * interaction is a Coulomb term with a distance-dependent dielectric ε_r = r, so
 * it is short-ranged (∝ 1/r²) and needs no cut-off of its own: it rides the same
 * non-bonded pair list, the same taper, excludes bonded and 1-3 pairs and halves
 * a 1-4 pair, exactly as the Lennard-Jones term does. A tangent continuation
 * below MIN_DISTANCE_RATIO·σ keeps a pile-up finite.
 */
const COULOMB_K = 332.0637;
const CHARGE_TRANSFER = 0.16;
const ELECTRONEGATIVITY: { [el: string]: number } = {
    H: 2.20, Li: 0.98, Be: 1.57, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98,
    Na: 0.93, Mg: 1.31, Al: 1.61, Si: 1.90, P: 2.19, S: 2.58, Cl: 3.16,
    K: 0.82, Ca: 1.00, Ga: 1.81, Ge: 2.01, As: 2.18, Se: 2.55, Br: 2.96,
    Rb: 0.82, Sr: 0.95, In: 1.78, Sn: 1.96, Sb: 2.05, Te: 2.10, I: 2.66,
    Sc: 1.36, Ti: 1.54, V: 1.63, Cr: 1.66, Mn: 1.55, Fe: 1.83, Co: 1.88,
    Ni: 1.91, Cu: 1.90, Zn: 1.65, Zr: 1.33, Nb: 1.60, Mo: 2.16, Ru: 2.20,
    Rh: 2.28, Pd: 2.20, Ag: 1.93, Cd: 1.69, Pt: 2.28, Au: 2.54, Hg: 2.00,
};
const ELECTRONEGATIVITY_DEFAULT = 2.20;
/** Bond-order weighting of a bond's charge transfer (by classified order). */
function chargeTransferScale(order: number): number {
    if (order >= 3) {
        return 1.6;
    }
    if (order >= 2) {
        return 1.5;
    }
    if (order >= 1.5) {
        return 1.35;
    }
    return 1.0;
}

/**
 * ---- Polarisation (damped induced-dipole attraction) ----
 * Two atoms with charges q_i, q_j and static dipole polarizabilities α_i, α_j
 * attract through the dipoles their charges induce by
 * U = −k·(α_i·q_j² + α_j·q_i²)/r⁴, damped so the term falls smoothly to zero
 * below ≈0.6 σ (a pile-up cannot reach −∞) and to zero at the non-bonded cutoff.
 * k is calibrated so α = 1.45 Å³, q = 0.4 e at r = 2 Å gives ≈ −0.5 kcal/mol.
 */
const POLARIZABILITY: { [el: string]: number } = {
    H: 0.667, Li: 24.3, Be: 5.60, B: 3.03, C: 1.76, N: 1.10, O: 0.802, F: 0.557,
    Na: 24.1, Mg: 10.6, Al: 6.80, Si: 5.38, P: 3.63, S: 2.90, Cl: 2.18,
    K: 43.4, Ca: 22.8, Ga: 8.12, Ge: 6.07, As: 4.31, Se: 3.77, Br: 3.05,
    Rb: 47.3, Sr: 27.6, In: 10.2, Sn: 7.70, Sb: 6.60, Te: 5.50, I: 5.35,
    Sc: 17.8, Ti: 14.6, V: 12.4, Cr: 11.6, Mn: 9.4, Fe: 8.4, Co: 7.5,
    Ni: 6.8, Cu: 6.2, Zn: 5.75, Zr: 12.5, Nb: 10.8, Mo: 9.6, Ru: 8.6,
    Rh: 7.5, Pd: 6.6, Ag: 6.5, Cd: 6.4, Pt: 6.2, Au: 5.8, Hg: 5.7,
};
const POLARIZABILITY_DEFAULT = 1.0;
const POL_K = 17.24;
/** Polarisation damping begins below this fraction of the pair's σ. */
const POL_DAMP_FRACTION = 0.6;

/**
 * ---- Hydrogen bond ----
 * A donor–H···acceptor contact (D and A each one of N, O, F) is described by a
 * 10-12 well blended over the hydrogen–acceptor Lennard-Jones pair: the well
 * replaces the LJ inside roughly 1.6–2.8 Å and within a D–H···A angle window
 * (≈150°–120°), and the blend is smooth in both r and the angle, so the total is
 * continuous and, with the hydrogen-bond flag off, is exactly the plain LJ.
 */
const HBOND_EPSILON = 2.0;
const HBOND_R0 = 1.90;
const HBOND_DIST_ON = 1.60;
const HBOND_DIST_OFF = 2.80;
const HBOND_ANGLE_ON = (150.0 * Math.PI) / 180.0;
const HBOND_ANGLE_OFF = (120.0 * Math.PI) / 180.0;

/**
 * ---- Metal coordination ----
 * The generic angle table knows nothing about metals, so a metal centre would be
 * bent toward the table default (≈109°). With the coordination term on, a metal's
 * ligand–metal–ligand angles are instead restrained toward the ideal geometry of
 * its coordination number — linear (2), trigonal (3), tetrahedral (4, or square
 * planar for the d⁸ metals) — and, for five or six ligands, toward a term with
 * minima at both cis (90°) and trans (180°) so an octahedron is not compressed.
 * Five-coordinate centres are placed on the same cis/trans term (a light
 * approximation of the trigonal-bipyramidal ideal). The metal–ligand distances
 * are left to the ordinary bond term, which already uses the empirical target
 * this shell is measured with.
 */
const METALS: { [el: string]: boolean } = {
    Sc: true, Ti: true, V: true, Cr: true, Mn: true, Fe: true, Co: true,
    Ni: true, Cu: true, Zn: true, Zr: true, Nb: true, Mo: true, Ru: true,
    Rh: true, Pd: true, Ag: true, Cd: true, Pt: true, Au: true, Hg: true,
};
/** d⁸ metals that prefer a square-planar four-coordinate geometry. */
const SQUARE_PLANAR_METALS: { [el: string]: boolean } = {
    Ni: true, Pd: true, Pt: true, Au: true,
};
const COORD_ANGLE_K = 30.0;
const COORD_OCTA_K = 15.0;

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
    /** Include the light Coulomb term (default true). */
    electrostatic?: boolean;
    /** Include the hydrogen-bond term (default true). */
    hbond?: boolean;
    /** Include the metal coordination geometry term (default false). */
    coordination?: boolean;
    /** Include the damped induced-dipole polarisation term (default false). */
    polarization?: boolean;
    /** Non-bonded cutoff in Å (4–14, default 8): the pair is switched off here,
     *  tapering from 2 Å inside it. The pair list is built 2 Å past it. */
    cutoff?: number;
    /** Name of a preset to take as the *base* for every key this object does not
     *  set itself ('low' | 'medium' | 'high'); an explicit key always wins, and
     *  an unknown name is ignored. */
    preset?: string;
}

/** One precision preset: a named, ready-made set of run settings. */
export interface MmPreset {
    key: string;
    label: string;
    options: MmOptimizeOptions;
}

/**
 * The three precision presets the panel offers, in the module (not the client),
 * so a preset cannot drift from what the optimiser actually does. Medium is
 * exactly the module's own defaults and the two others are the same knobs moved
 * together, the heavier terms switched in as the level rises:
 *
 *   low     bonds and angles only (every distance- and dihedral-dependent term
 *           off), a wide convergence target and a short budget — a rough
 *           pre-optimisation before a real calculation;
 *   medium  the full geometric field — Lennard-Jones, the light Coulomb term and
 *           the hydrogen-bond term — at the module's own defaults, i.e. the
 *           balanced setting the panel starts on;
 *   high    everything medium has plus the metal coordination and polarisation
 *           terms, a five-times tighter convergence target, a shorter step and a
 *           wider non-bonded cutoff (a geometry to hand to a QM program).
 *
 * Every value here is one of the values the corresponding control offers, so a
 * preset is always exactly what the panel then shows.
 */
export const MM_PRESETS: MmPreset[] = [
    {
        key: 'low',
        label: 'Low — quick pre-optimisation',
        options: {
            nonbonded: false, torsions: false, impropers: true,
            electrostatic: false, hbond: false, coordination: false, polarization: false,
            cutoff: 6,
            gradientTol: 0.4, maxIterations: 120, timeBudgetMs: 1500, maxStep: 0.30,
        },
    },
    {
        key: 'medium',
        label: 'Medium — balanced (default)',
        options: {
            nonbonded: true, torsions: true, impropers: true,
            electrostatic: true, hbond: true, coordination: false, polarization: false,
            cutoff: DEFAULT_CUTOFF,
            gradientTol: DEFAULT_GRADIENT_TOL, maxIterations: DEFAULT_MAX_ITERATIONS,
            timeBudgetMs: DEFAULT_TIME_BUDGET_MS, maxStep: MAX_STEP,
        },
    },
    {
        key: 'high',
        label: 'High — careful',
        options: {
            nonbonded: true, torsions: true, impropers: true,
            electrostatic: true, hbond: true, coordination: true, polarization: true,
            cutoff: 10,
            gradientTol: 0.01, maxIterations: 2000, timeBudgetMs: 20000, maxStep: 0.05,
        },
    },
];

/** The preset a request names, or an empty set when it names none. */
function presetOf(raw: unknown): MmOptimizeOptions {
    if (typeof raw !== 'string') {
        return {};
    }
    for (let i = 0; i < MM_PRESETS.length; i++) {
        if (MM_PRESETS[i].key === raw) {
            return MM_PRESETS[i].options;
        }
    }
    return {};
}

/** The term switches, resolved to their values for one run. */
export interface MmTermFlags {
    nonbonded: boolean;
    torsions: boolean;
    impropers: boolean;
    /** Light Coulomb term with a distance-dependent dielectric (default on). */
    electrostatic: boolean;
    /** Donor–H···acceptor hydrogen-bond term (default on). */
    hbond: boolean;
    /** Metal coordination geometry term (default off). */
    coordination: boolean;
    /** Damped induced-dipole polarisation term (default off). */
    polarization: boolean;
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

/**
 * One ligand–metal–ligand angle restraint: `kind` 0 is a plain harmonic toward
 * `theta0`, kind 1 is the cis/trans term with minima at 90° and 180° (used for
 * the higher coordination numbers, where the ideal polyhedron has both).
 */
interface CoordTerm {
    m: number;
    i: number;
    j: number;
    k: number;
    kind: number;
    theta0: number;
}

/** Everything the force field is made of, built once per call. */
interface Field {
    bonds: BondTerm[];
    angles: AngleTerm[];
    torsions: TorsionTerm[];
    impropers: ImproperTerm[];
    nonbonded: NonbondedTerm[];
    coordination: CoordTerm[];
    counts: MmTermCounts;
    /** Per-atom stiffness (preconditioner input). */
    stiffness: Float64Array;
    /** Non-bonded settings this list was built with: the pair is gone at
     *  `cutoff` and the list itself reaches `cutoff + LIST_BUFFER`. */
    cutoff: number;
    taperStart: number;
    /** The resolved term switches, so the evaluation can read them. */
    modes: MmTermFlags;
    /** Per-atom partial charge (e), from the electronegativity model. */
    charges: Float64Array;
    /** Per-atom static dipole polarizability (Å³). */
    polarizability: Float64Array;
    /** For an H atom, the index of the N/O/F it is bonded to, else −1. */
    donorOf: Int32Array;
    /** 1 when the atom is a hydrogen-bond acceptor (N, O or F). */
    isAcceptor: Uint8Array;
    /** 1 when the atom is a hydrogen. */
    isHydrogen: Uint8Array;
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

function electronegativity(el: string): number {
    const v = ELECTRONEGATIVITY[el];
    return v === undefined ? ELECTRONEGATIVITY_DEFAULT : v;
}

function polarizabilityOf(el: string): number {
    const v = POLARIZABILITY[el];
    return v === undefined ? POLARIZABILITY_DEFAULT : v;
}

/** A hydrogen-bond donor/acceptor element (N, O or F). */
function isDonorAcceptorElement(el: string): boolean {
    return el === 'N' || el === 'O' || el === 'F';
}

function isMetal(el: string): boolean {
    return METALS[el] === true;
}

/**
 * Partial charges from the elements and bond orders alone: each bond shifts a
 * little charge toward its more electronegative end, and the whole set is then
 * shifted so the structure is neutral (an ion pair with no bonds stays neutral).
 */
function computeCharges(
    elements: string[],
    bonds: { i: number; j: number; order: number }[],
    n: number
): Float64Array {
    const q = new Float64Array(n);
    for (let b = 0; b < bonds.length; b++) {
        const bond = bonds[b];
        // d > 0 when the first atom is the more electronegative one, and the
        // charge moves toward it: that atom gains electron density (goes negative).
        const d = CHARGE_TRANSFER *
            (electronegativity(elements[bond.i]) - electronegativity(elements[bond.j])) *
            chargeTransferScale(bond.order);
        q[bond.i] -= d;
        q[bond.j] += d;
    }
    let sum = 0.0;
    for (let i = 0; i < n; i++) {
        sum += q[i];
    }
    const mean = n > 0 ? sum / n : 0.0;
    for (let i = 0; i < n; i++) {
        q[i] -= mean;
    }
    return q;
}

/**
 * The ideal ligand–metal–ligand angle for a coordination number: kind 0 is a
 * plain harmonic toward `theta0`, kind 1 is the cis/trans term (minima at 90°
 * and 180°) used where the ideal polyhedron has both.
 */
function coordinationIdeal(el: string, cn: number): { kind: number; theta0: number } {
    if (cn === 2) {
        return { kind: 0, theta0: Math.PI };
    }
    if (cn === 3) {
        return { kind: 0, theta0: (2.0 * Math.PI) / 3.0 };
    }
    if (cn === 4 && !SQUARE_PLANAR_METALS[el]) {
        return { kind: 0, theta0: (109.47122063449069 * Math.PI) / 180.0 };
    }
    return { kind: 1, theta0: 0.0 };
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
    settings: MMResolvedOptions
): Field {
    const flags = settings.flags;
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

    // ---- per-atom topology for the distance-dependent terms ----
    // Charges, polarizabilities and the hydrogen-bond donor/acceptor flags come
    // from the elements and the bond list alone, so they never change during a
    // run and are computed once here.
    const charges = computeCharges(elements, bonds, n);
    const polarizability = new Float64Array(n);
    const isAcceptor = new Uint8Array(n);
    const isHydrogen = new Uint8Array(n);
    const donorOf = new Int32Array(n);
    donorOf.fill(-1);
    for (let i = 0; i < n; i++) {
        const el = elements[i];
        polarizability[i] = polarizabilityOf(el);
        isAcceptor[i] = isDonorAcceptorElement(el) ? 1 : 0;
        isHydrogen[i] = el === 'H' ? 1 : 0;
    }
    // A hydrogen bonded to exactly one N/O/F is a donor hydrogen; the donor is
    // that neighbour, which the hydrogen-bond angle window is measured through.
    for (let i = 0; i < n; i++) {
        if (isHydrogen[i] === 0) {
            continue;
        }
        const list = nb[i];
        if (list.length === 1 && isDonorAcceptorElement(elements[list[0]])) {
            donorOf[i] = list[0];
        }
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
        // A metal centre's angles belong to the coordination term when it is on:
        // the generic table has no metal entries and would bend the ligands
        // toward its default, which is not the metal's geometry.
        if (flags.coordination && isMetal(elements[c])) {
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

    // ---- metal coordination (ligand–metal–ligand geometry) ----
    const coordinationOut: CoordTerm[] = [];
    if (flags.coordination) {
        for (let m = 0; m < n; m++) {
            if (!isMetal(elements[m])) {
                continue;
            }
            const list = nb[m];
            const cn = list.length;
            if (cn < 2) {
                continue;
            }
            const ideal = coordinationIdeal(elements[m], cn);
            const kc = ideal.kind === 0 ? COORD_ANGLE_K : COORD_OCTA_K;
            for (let a = 0; a < cn - 1; a++) {
                for (let b2 = a + 1; b2 < cn; b2++) {
                    const i = list[a];
                    const j = list[b2];
                    if (!mobile[m] && !mobile[i] && !mobile[j]) {
                        continue;
                    }
                    coordinationOut.push({ m: m, i: i, j: j, k: kc, kind: ideal.kind, theta0: ideal.theta0 });
                    stiffness[m] += 0.5 * kc;
                    stiffness[i] += 0.25 * kc;
                    stiffness[j] += 0.25 * kc;
                }
            }
        }
    }

    // ---- non-bonded pair list ----
    // Built whenever any distance-dependent term is on (Lennard-Jones, Coulomb,
    // hydrogen bond or polarisation), since they all ride the same pair list and
    // the same exclusions.
    const anyDistance = flags.nonbonded || flags.electrostatic || flags.hbond || flags.polarization;
    const nonbondedOut = anyDistance
        ? buildNonbondedList(x, y, z, elements, nb, mobile, n, settings.listCutoff)
        : [];

    return {
        bonds: bondsOut,
        angles: anglesOut,
        torsions: torsionsOut,
        impropers: impropersOut,
        nonbonded: nonbondedOut,
        coordination: coordinationOut,
        counts: {
            bonds: bondsOut.length,
            angles: anglesOut.length,
            torsions: torsionsOut.length,
            impropers: impropersOut.length,
            nonbonded: nonbondedOut.length,
        },
        stiffness: stiffness,
        cutoff: settings.cutoff,
        taperStart: settings.taperStart,
        modes: flags,
        charges: charges,
        polarizability: polarizability,
        donorOf: donorOf,
        isAcceptor: isAcceptor,
        isHydrogen: isHydrogen,
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
 * Every non-bonded pair within `listCutoff` (the run's cutoff plus the buffer
 * that keeps the list valid between rebuilds) that involves a mobile atom,
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
    n: number,
    listCutoff: number
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
    // One cell covers the whole list cutoff, so the 27-cell sweep below sees
    // every candidate pair whatever the cutoff is.
    const cells = new Map<string, number[]>();
    const cellOf = new Int32Array(n * 3);
    for (let i = 0; i < n; i++) {
        const cx = Math.floor(x[i] / listCutoff);
        const cy = Math.floor(y[i] / listCutoff);
        const cz = Math.floor(z[i] / listCutoff);
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
    const cut2 = listCutoff * listCutoff;
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

    // ---- metal coordination ----
    // Same angle geometry as the angle-bend term, but kind 1 uses a term with
    // minima at both 90° and 180° so a cis/trans polyhedron is not compressed.
    for (let t = 0; t < field.coordination.length; t++) {
        const term = field.coordination[t];
        let ux = x[term.i] - x[term.m];
        let uy = y[term.i] - y[term.m];
        let uz = z[term.i] - z[term.m];
        let vx = x[term.j] - x[term.m];
        let vy = y[term.j] - y[term.m];
        let vz = z[term.j] - z[term.m];
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
        let termEnergy: number;
        let f = 0.0;
        if (term.kind === 1) {
            const c = cosTheta;
            termEnergy = term.k * c * c * (c + 1.0) * (c + 1.0);
            // For a cosine-basis energy the 1/sinθ of the angle gradient cancels,
            // so f = (dE/dθ)/sinθ = −dE/dc here (no clamps needed).
            f = -2.0 * term.k * c * (c + 1.0) * (2.0 * c + 1.0);
        } else {
            const sinT = Math.sqrt(1.0 - cosTheta * cosTheta);
            const dev = Math.acos(cosTheta) - term.theta0;
            termEnergy = 0.5 * term.k * dev * dev;
            f = sinT < 1e-6 ? 0.0 : (term.k * dev) / sinT;
        }
        energy += termEnergy;
        if (!gx || f === 0.0) {
            continue;
        }
        const giX = (-(vx - cosTheta * ux) * f) / lu;
        const giY = (-(vy - cosTheta * uy) * f) / lu;
        const giZ = (-(vz - cosTheta * uz) * f) / lu;
        const gjX = (-(ux - cosTheta * vx) * f) / lv;
        const gjY = (-(uy - cosTheta * vy) * f) / lv;
        const gjZ = (-(uz - cosTheta * vz) * f) / lv;
        addGrad(term.i, giX, giY, giZ);
        addGrad(term.j, gjX, gjY, gjZ);
        addGrad(term.m, -(giX + gjX), -(giY + gjY), -(giZ + gjZ));
    }

    // ---- non-bonded ----
    for (let t = 0; t < field.nonbonded.length; t++) {
        const term = field.nonbonded[t];
        let dx = x[term.i] - x[term.j];
        let dy = y[term.i] - y[term.j];
        let dz = z[term.i] - z[term.j];
        let r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r >= field.cutoff) {
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
        const useLJ = field.modes.nonbonded;
        let ratio = term.sigma / r;
        if (term.repulsionOnly) {
            // A 1-3 clash belongs to the Lennard-Jones family: with the LJ term
            // off there is nothing to add for it, and it carries no attractive
            // term of its own either.
            if (!useLJ) {
                continue;
            }
            // Pure r^-12 repulsion, kept finite (and still separating) by the
            // same tangent continuation the main path uses.
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
        // A switching function tapers every distance-dependent term to zero
        // between the taper start and the cutoff instead of chopping it off: a
        // bare cutoff leaves a step of the energy at the boundary (≈6e-4 kcal/mol
        // for the LJ at 8 Å) and a pair crossing it makes the energy
        // discontinuous, which no line search can step through.
        let switchValue = 1.0;
        let dSwitch = 0.0;
        if (r > field.taperStart) {
            const rOn2 = field.taperStart * field.taperStart;
            const rOff2 = field.cutoff * field.cutoff;
            const denom = rOff2 - rOn2;
            const u = rOff2 - r * r;
            switchValue = (u * u * (rOff2 + 2.0 * r * r - 3.0 * rOn2)) /
                (denom * denom * denom);
            dSwitch = (12.0 * r * u * (rOn2 - r * r)) / (denom * denom * denom);
        }
        // ---- Lennard-Jones (the pair's own contribution) ----
        let ljContribution = 0.0;
        let dLjContribution = 0.0;
        if (useLJ && ratio > 1.0 / MIN_DISTANCE_RATIO) {
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
            ljContribution = term.scale * term.epsilon * (ljc + slope * (r - rc));
            dLjContribution = term.scale * term.epsilon * slope;
        } else if (useLJ) {
            const s6 = Math.pow(ratio, 6);
            const s12 = s6 * s6;
            const lj = s12 - 2.0 * s6;
            ljContribution = term.scale * term.epsilon * lj * switchValue;
            // dE/dr of (LJ·S): the LJ part, plus the taper's own slope.
            dLjContribution = term.scale * term.epsilon *
                ((12.0 * (s6 - s12) * switchValue) / r + lj * dSwitch);
        }
        energy += ljContribution;
        if (gx) {
            addGrad(term.i, (dLjContribution * dx) / r, (dLjContribution * dy) / r,
                (dLjContribution * dz) / r);
            addGrad(term.j, (-dLjContribution * dx) / r, (-dLjContribution * dy) / r,
                (-dLjContribution * dz) / r);
        }

        // ---- light Coulomb term (distance-dependent dielectric ε_r = r) ----
        if (field.modes.electrostatic) {
            const qi = field.charges[term.i];
            const qj = field.charges[term.j];
            // Below MIN_DISTANCE_RATIO·σ the value is clamped at that radius so a
            // pile-up stays finite; the gradient is zero there (the LJ's own
            // tangent is what pushes such a pair apart).
            const floor = term.sigma * MIN_DISTANCE_RATIO;
            const rE = r > floor ? r : floor;
            const coulomb = (COULOMB_K * qi * qj) / (rE * rE);
            energy += term.scale * switchValue * coulomb;
            if (gx && r > floor) {
                const dCoulDr = (-2.0 * COULOMB_K * qi * qj) / (rE * rE * rE);
                const dedr = term.scale * (dCoulDr * switchValue + coulomb * dSwitch);
                addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
                addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
            }
        }

        // ---- damped induced-dipole polarisation ----
        if (field.modes.polarization) {
            const qi = field.charges[term.i];
            const qj = field.charges[term.j];
            const coef = field.polarizability[term.i] * qj * qj +
                field.polarizability[term.j] * qi * qi;
            if (coef > 0.0) {
                const rD = POL_DAMP_FRACTION * term.sigma;
                let damp = 1.0;
                let dDampDr = 0.0;
                if (r < rD) {
                    // Smoothstep in t = r/r_D, squared so the 1/r⁴ stays finite
                    // as r → 0 (U tends to a constant, not −∞).
                    const t = r / rD;
                    const s = t * t * (3.0 - 2.0 * t);
                    damp = s * s;
                    dDampDr = (2.0 * s * (6.0 * t * (1.0 - t))) / rD;
                }
                const inv = 1.0 / (r * r * r * r);
                const base = -POL_K * coef * inv;
                energy += term.scale * switchValue * base * damp;
                if (gx) {
                    const dBaseDr = (4.0 * POL_K * coef * inv) / r;
                    const dedr = term.scale *
                        ((dBaseDr * damp + base * dDampDr) * switchValue + base * damp * dSwitch);
                    addGrad(term.i, (dedr * dx) / r, (dedr * dy) / r, (dedr * dz) / r);
                    addGrad(term.j, (-dedr * dx) / r, (-dedr * dy) / r, (-dedr * dz) / r);
                }
            }
        }

        // ---- hydrogen bond (10-12 well blended over the LJ pair) ----
        if (field.modes.hbond) {
            let h = -1;
            let a = -1;
            let d = -1;
            if (field.isHydrogen[term.i] === 1 && field.isAcceptor[term.j] === 1 &&
                field.donorOf[term.i] >= 0) {
                h = term.i; a = term.j; d = field.donorOf[term.i];
            } else if (field.isHydrogen[term.j] === 1 && field.isAcceptor[term.i] === 1 &&
                field.donorOf[term.j] >= 0) {
                h = term.j; a = term.i; d = field.donorOf[term.j];
            }
            if (h >= 0) {
                // The D–H···A angle is measured at the hydrogen.
                const ux = x[d] - x[h];
                const uy = y[d] - y[h];
                const uz = z[d] - z[h];
                const vx = x[a] - x[h];
                const vy = y[a] - y[h];
                const vz = z[a] - z[h];
                const lu = Math.sqrt(ux * ux + uy * uy + uz * uz);
                const lv = Math.sqrt(vx * vx + vy * vy + vz * vz);
                if (lu > 1e-9 && lv > 1e-9) {
                    let cosT = (ux * vx + uy * vy + uz * vz) / (lu * lv);
                    if (cosT > 1.0) {
                        cosT = 1.0;
                    } else if (cosT < -1.0) {
                        cosT = -1.0;
                    }
                    const theta = Math.acos(cosT);
                    let wR = 1.0;
                    let dwR = 0.0;
                    if (r >= HBOND_DIST_OFF) {
                        wR = 0.0;
                    } else if (r > HBOND_DIST_ON) {
                        const t = (HBOND_DIST_OFF - r) / (HBOND_DIST_OFF - HBOND_DIST_ON);
                        wR = t * t * (3.0 - 2.0 * t);
                        dwR = -(6.0 * t * (1.0 - t)) / (HBOND_DIST_OFF - HBOND_DIST_ON);
                    }
                    let wA = 1.0;
                    let dwA = 0.0;
                    if (theta <= HBOND_ANGLE_OFF) {
                        wA = 0.0;
                    } else if (theta < HBOND_ANGLE_ON) {
                        const t = (theta - HBOND_ANGLE_OFF) / (HBOND_ANGLE_ON - HBOND_ANGLE_OFF);
                        wA = t * t * (3.0 - 2.0 * t);
                        dwA = (6.0 * t * (1.0 - t)) / (HBOND_ANGLE_ON - HBOND_ANGLE_OFF);
                    }
                    const w = wR * wA;
                    if (w > 0.0) {
                        // Blend the 10-12 well over the pair's LJ: since the LJ
                        // contribution is already added, only the difference is
                        // added here, so the term is continuous and reduces to
                        // the plain LJ wherever the window is closed.
                        const rho = HBOND_R0 / r;
                        const rho10 = Math.pow(rho, 10);
                        const wellShape = 5.0 * rho10 * rho * rho - 6.0 * rho10;
                        const well = term.scale * HBOND_EPSILON * wellShape * switchValue;
                        const dWell = term.scale * HBOND_EPSILON *
                            ((60.0 * HBOND_R0 * (rho10 / rho - rho10 * rho)) / (r * r) * switchValue +
                                wellShape * dSwitch);
                        const c = well - ljContribution;
                        const dCdr = dWell - dLjContribution;
                        energy += w * c;
                        if (gx) {
                            const dwdr = wA * dwR;
                            const dwdtheta = wR * dwA;
                            const radial = w * dCdr + c * dwdr;
                            const ahx = x[h] - x[a];
                            const ahy = y[h] - y[a];
                            const ahz = z[h] - z[a];
                            addGrad(h, (radial * ahx) / r, (radial * ahy) / r, (radial * ahz) / r);
                            addGrad(a, (-radial * ahx) / r, (-radial * ahy) / r, (-radial * ahz) / r);
                            const angular = c * dwdtheta;
                            if (angular !== 0.0) {
                                const sinT = Math.sqrt(Math.max(0.0, 1.0 - cosT * cosT));
                                if (sinT > 1e-6) {
                                    const g = angular / sinT;
                                    const ddX = (-(vx / lv - cosT * ux / lu) * g) / lu;
                                    const ddY = (-(vy / lv - cosT * uy / lu) * g) / lu;
                                    const ddZ = (-(vz / lv - cosT * uz / lu) * g) / lu;
                                    const daX = (-(ux / lu - cosT * vx / lv) * g) / lv;
                                    const daY = (-(uy / lu - cosT * vy / lv) * g) / lv;
                                    const daZ = (-(uz / lu - cosT * vz / lv) * g) / lv;
                                    addGrad(d, ddX, ddY, ddZ);
                                    addGrad(a, daX, daY, daZ);
                                    addGrad(h, -(ddX + daX), -(ddY + daY), -(ddZ + daZ));
                                }
                            }
                        }
                    }
                }
            }
        }
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
function readInput(rawAtoms: unknown, rawBonds: unknown, rawSelected: unknown): MmInput {
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
    /** Non-bonded cutoff (Å): the pair is gone here. */
    cutoff: number;
    /** The taper starts here (cutoff − taper width). */
    taperStart: number;
    /** The pair list is built to here (cutoff + buffer). */
    listCutoff: number;
    flags: MmTermFlags;
}

/**
 * Validate a caller's settings. Anything missing, of the wrong type, or outside
 * its range falls back to `preset`'s value when one was named, and to the
 * module's own default otherwise, so an arbitrary message payload cannot make a
 * run unreasonable (and a panel's choices are the only thing that can change
 * the behaviour). An explicitly given key always wins over the preset.
 */
function readOptions(rawOptions: unknown): MMResolvedOptions {
    const o = (rawOptions && typeof rawOptions === 'object' ? rawOptions : {}) as MmOptimizeOptions;
    const base = presetOf(o.preset);
    const num = (v: unknown, fallback: number, lo: number, hi: number): number =>
        typeof v === 'number' && isFinite(v) ? Math.max(lo, Math.min(hi, v)) : fallback;
    const flag = (v: unknown, fallback: boolean): boolean => typeof v === 'boolean' ? v : fallback;
    const pick = (key: 'maxIterations' | 'timeBudgetMs' | 'gradientTol' | 'maxStep' | 'cutoff' |
        'nonbonded' | 'torsions' | 'impropers' | 'electrostatic' | 'hbond' |
        'coordination' | 'polarization'): unknown =>
        o[key] !== undefined ? o[key] : base[key];
    const cutoff = num(pick('cutoff'), DEFAULT_CUTOFF, 4, 14);
    return {
        maxIterations: Math.trunc(num(pick('maxIterations'), DEFAULT_MAX_ITERATIONS, 1, 100000)),
        timeBudgetMs: num(pick('timeBudgetMs'), DEFAULT_TIME_BUDGET_MS, 10, 60000),
        gradientTol: num(pick('gradientTol'), DEFAULT_GRADIENT_TOL, 1e-6, 1000),
        maxStep: num(pick('maxStep'), MAX_STEP, 0.005, 2.0),
        cutoff: cutoff,
        taperStart: cutoff - DEFAULT_TAPER_WIDTH,
        listCutoff: cutoff + LIST_BUFFER,
        flags: {
            nonbonded: flag(pick('nonbonded'), true),
            torsions: flag(pick('torsions'), true),
            impropers: flag(pick('impropers'), true),
            electrostatic: flag(pick('electrostatic'), true),
            hbond: flag(pick('hbond'), true),
            coordination: flag(pick('coordination'), false),
            polarization: flag(pick('polarization'), false),
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
    const field = buildField(input.elements, input.x, input.y, input.z, input.bonds, mobile, n, settings);
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
 * Relax the selected atoms of a structure — as a *stateful run*.
 * ============================================================
 * The optimiser used to be one call that ran to completion; a caller that has
 * to watch it work (a progress curve), interrupt it (pause), continue it, or
 * abandon it (the view moved on) needs the phases apart. The search itself is
 * unchanged — one implementation, driven two ways:
 *
 *   beginOptimization(...)   validate, build the force field, take the first
 *                            energy and force reading; it either refuses, with
 *                            the complete answer that says why, or hands back
 *                            a session to steer;
 *   session.step(budget)     run at most `budget` passes of the search and
 *                            report where it stands (energy, forces, steps and
 *                            the current coordinates of the mobile atoms);
 *   session.finish()         the complete result — valid at any time: after
 *                            convergence, after a budget ran out, or halfway
 *                            through, which is exactly what a pause reports.
 *
 * session.pause()/resume() stop the wall-clock budget counting while nobody is
 * stepping, and session.abandon() drops a run outright. `optimizeSelection`
 * below is the old single call expressed through those phases
 * (begin → step(∞) → finish), so a caller that only wants an answer still gets
 * one — and the two paths cannot disagree, because they are the same code.
 */

/** One point of a run's convergence history (the live curve's source). */
export interface MmHistoryPoint {
    /** Accepted steps taken; 0 is the geometry the run was given. */
    step: number;
    energy: number;
    maxGrad: number;
    rmsGrad: number;
}

/** What beginning a run produced: a session to step, or the answer itself. */
export interface MmBeginResult {
    /** The run to step — null when `result` already holds the whole answer
     *  (a refused call, or a selection with nothing to optimise). */
    session: MmSession | null;
    result: MmOptimizeResult | null;
}

/** The state of a run after at most `budget` passes of its search. */
export interface MmStepResult {
    /** Accepted geometry steps taken so far. */
    step: number;
    energy: number;
    /** Largest / RMS per-atom force at the current geometry, kcal/(mol·Å). */
    maxGrad: number;
    rmsGrad: number;
    /** The run has stopped; stepping a finished session does nothing. */
    done: boolean;
    converged: boolean;
    stopReason: string;
    /** Current coordinates of the run's atoms, in the caller's order. */
    coords: number[][];
}

/** Everything `readInput` produced, as the session consumes it. */
interface MmInput {
    n: number;
    elements: string[];
    x: Float64Array;
    y: Float64Array;
    z: Float64Array;
    bonds: { i: number; j: number; order: number }[];
    selected: number[];
    invalid: number;
}

/** The working geometry and the buffers one search pass writes. */
interface Work {
    x: Float64Array;
    y: Float64Array;
    z: Float64Array;
    trialX: Float64Array;
    trialY: Float64Array;
    trialZ: Float64Array;
    grad: Gradient;
    /** Preconditioner input: inverse stiffness per atom. */
    precondition: Float64Array;
    /** Per-atom drift since the pair list was built: how far each atom has
     *  moved in total, so the list is rebuilt before a pair can slip into
     *  range unseen. */
    drift: Float64Array;
    /** Step prefactor: starts at LAMBDA_START, grows on accepted steps, halves
     *  on rejected ones. It is part of the run's state, not of one pass. */
    lambda: number;
}

/**
 * Conjugate-gradient bookkeeping: the search direction of the last accepted
 * step, the gradient it was built from, and whether that direction may be
 * conjugated at all. A plain steepest descent converges far too slowly on a
 * flat non-bonded landscape (hundreds of steps to polish off the last tenth of
 * a kcal), while the monotone acceptance still guarantees that every step
 * lowers the energy — the direction only decides how fast, never whether.
 */
interface CgState {
    dirX: Float64Array;
    dirY: Float64Array;
    dirZ: Float64Array;
    prevDirX: Float64Array;
    prevDirY: Float64Array;
    prevDirZ: Float64Array;
    prevGradX: Float64Array;
    prevGradY: Float64Array;
    prevGradZ: Float64Array;
    valid: boolean;
    /** A conjugate direction that led nowhere is retried once as a plain one. */
    resetUsed: boolean;
}

/**
 * One run in progress: the input it was built from, the force field, the
 * working geometry and the bookkeeping of the search. A caller reads only the
 * public fields it displays (`indices`, `atomCount`, `history`, `iterations`,
 * `energy`, `energyBefore`, `done`, `converged`, `stopReason`) and calls
 * `step` / `finish` / `pause` / `resume` / `abandon`; every buffer stays
 * private, so a half-finished run cannot be poked into an inconsistent state
 * from outside.
 */
export class MmSession {
    /** Atoms that may move, in the caller's order. */
    readonly indices: number[];
    /** Whole-structure atom count the run was built from. */
    readonly atomCount: number;
    /** One point per accepted step; step 0 is the input geometry. */
    readonly history: MmHistoryPoint[] = [];
    /** Accepted geometry steps so far. */
    iterations = 0;
    /** True once the run has stopped, whatever stopped it. */
    done = false;
    /** True when it stopped on the force tolerance. */
    converged = false;
    /** Why it stopped: 'gradient' | 'energy' | 'iterations' | 'time' | 'step'. */
    stopReason = 'gradient';
    /** Energy of the mobile subsystem's terms at the current geometry. */
    energy: number;
    /** Energy at the geometry the run was given. */
    readonly energyBefore: number;

    private readonly settings: MMResolvedOptions;
    private readonly elements: string[];
    private readonly bonds: { i: number; j: number; order: number }[];
    private readonly mobile: boolean[];
    private readonly inputX: Float64Array;
    private readonly inputY: Float64Array;
    private readonly inputZ: Float64Array;
    private readonly field: Field;
    private readonly work: Work;
    private readonly cg: CgState;
    /** Selected entries that were not usable (bad index, no coordinates). */
    private readonly invalidCount: number;
    private startTime: number;
    /** When the run was paused (0 = not paused), so resume can shift the clock. */
    private pauseStarted = 0;
    private lastEnergy: number;
    private maxForce: number;
    private rmsForce: number;
    private plateau = 0;

    constructor(input: MmInput, settings: MMResolvedOptions) {
        const n = input.n;
        this.settings = settings;
        this.atomCount = n;
        this.indices = input.selected.slice();
        this.invalidCount = input.invalid;
        this.elements = input.elements;
        this.bonds = input.bonds;
        const mobile: boolean[] = new Array(n);
        for (let i = 0; i < n; i++) {
            mobile[i] = false;
        }
        for (let s = 0; s < input.selected.length; s++) {
            mobile[input.selected[s]] = true;
        }
        this.mobile = mobile;
        this.inputX = new Float64Array(input.x);
        this.inputY = new Float64Array(input.y);
        this.inputZ = new Float64Array(input.z);
        const field = buildField(input.elements, input.x, input.y, input.z, input.bonds, mobile, n, settings);
        this.field = field;
        const precondition = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            precondition[i] = 1.0 / Math.max(field.stiffness[i], STIFFNESS_FLOOR);
        }
        // Working copies of the geometry: only the mobile atoms are ever written
        // back, but the buffers are whole-structure so every term reads the same
        // arrays the input was read into.
        this.work = {
            x: new Float64Array(input.x),
            y: new Float64Array(input.y),
            z: new Float64Array(input.z),
            trialX: new Float64Array(n),
            trialY: new Float64Array(n),
            trialZ: new Float64Array(n),
            grad: {
                gx: new Float64Array(n),
                gy: new Float64Array(n),
                gz: new Float64Array(n),
            },
            precondition: precondition,
            drift: new Float64Array(n),
            lambda: LAMBDA_START,
        };
        this.cg = {
            dirX: new Float64Array(n), dirY: new Float64Array(n), dirZ: new Float64Array(n),
            prevDirX: new Float64Array(n), prevDirY: new Float64Array(n), prevDirZ: new Float64Array(n),
            prevGradX: new Float64Array(n), prevGradY: new Float64Array(n), prevGradZ: new Float64Array(n),
            valid: false,
            resetUsed: false,
        };
        this.startTime = Date.now();
        this.lastEnergy = 0.0;
        this.maxForce = 0.0;
        this.rmsForce = 0.0;
        if (this.isEmptyField()) {
            // Nothing to optimise: no pass is ever run and the history stays
            // empty; beginOptimization turns this into its own answer.
            this.energy = 0.0;
            this.energyBefore = 0.0;
            return;
        }
        this.energy = evaluate(field, this.work.x, this.work.y, this.work.z, mobile, this.work.grad);
        this.energyBefore = this.energy;
        this.lastEnergy = this.energy;
        this.readStats();
        this.pushPoint();
    }

    /** True when the built field has no live term at all. */
    private isEmptyField(): boolean {
        const c = this.field.counts;
        return c.bonds === 0 && c.angles === 0 && c.torsions === 0 &&
            c.impropers === 0 && c.nonbonded === 0;
    }

    /** Read the current force statistics out of the gradient buffer. */
    private readStats(): void {
        const stats = forceStats(this.work.grad, this.mobile, this.atomCount);
        this.maxForce = stats.max;
        this.rmsForce = stats.rms;
    }

    /** Record one point of the convergence history (step 0 is the input). */
    private pushPoint(): void {
        if (this.history.length >= HISTORY_CAP) {
            // Halve the array — every other point, with step 0 and the newest
            // kept — so a very long run's curve stays bounded. A reader that
            // tracks the last step it saw only ever sees the survivors once.
            const shrunk: MmHistoryPoint[] = [];
            for (let i = 0; i < this.history.length; i += 2) {
                shrunk.push(this.history[i]);
            }
            const last = this.history[this.history.length - 1];
            if (shrunk[shrunk.length - 1] !== last) {
                shrunk.push(last);
            }
            this.history.length = 0;
            for (let i = 0; i < shrunk.length; i++) {
                this.history.push(shrunk[i]);
            }
        }
        this.history.push({
            step: this.iterations,
            energy: this.energy,
            maxGrad: this.maxForce,
            rmsGrad: this.rmsForce,
        });
    }

    /** Current coordinates of the run's atoms, in the caller's order. */
    coords(): number[][] {
        const out: number[][] = [];
        for (let s = 0; s < this.indices.length; s++) {
            const i = this.indices[s];
            out.push([this.work.x[i], this.work.y[i], this.work.z[i]]);
        }
        return out;
    }

    /** The answer of a selection with no live terms at all — null otherwise. */
    trivialResult(): MmOptimizeResult | null {
        if (!this.isEmptyField()) {
            return null;
        }
        return {
            ok: true,
            message: 'Molecular Mechanics: nothing to optimise — the selection has no bonds or contacts',
            atomCount: this.atomCount,
            indices: this.indices.slice(),
            coords: this.coords(),
            energyBefore: 0,
            energyAfter: 0,
            iterations: 0,
            maxGrad: 0,
            rmsGrad: 0,
            maxMove: 0,
            converged: true,
            stopReason: 'gradient',
            counts: this.field.counts,
            invalid: this.invalidCount,
        };
    }

    /** Stop the wall-clock budget counting while nobody is stepping the run. */
    pause(): void {
        if (this.pauseStarted === 0) {
            this.pauseStarted = Date.now();
        }
    }

    /** Resume after a pause: the time spent paused is added back to the budget. */
    resume(): void {
        if (this.pauseStarted !== 0) {
            this.startTime += Date.now() - this.pauseStarted;
            this.pauseStarted = 0;
        }
    }

    /** Drop the run: it cannot be stepped or finished into anything new. */
    abandon(): void {
        this.done = true;
        this.stopReason = 'cancelled';
    }

    /**
     * Run up to `budget` passes of the search (Infinity runs to completion) and
     * report the state it is in. A stopped run answers at once, so a caller may
     * keep asking without doing further work.
     */
    step(budget: number): MmStepResult {
        let used = 0;
        while (used < budget && !this.done) {
            this.iterate();
            used++;
        }
        return {
            step: this.iterations,
            energy: this.energy,
            maxGrad: this.maxForce,
            rmsGrad: this.rmsForce,
            done: this.done,
            converged: this.converged,
            stopReason: this.stopReason,
            coords: this.coords(),
        };
    }

    /**
     * One pass of the search: the stop checks, and — when none of them fires —
     * one full trial step. Returns false when the run has just stopped.
     */
    private iterate(): boolean {
        // Buffer aliases keep the pass below identical to the single-call form; the
        // arrays are the session's own, so writing through an alias is writing
        // into the run.
        const n = this.atomCount;
        const mobile = this.mobile;
        const field = this.field;
        const work = this.work;
        const cg = this.cg;
        const x = work.x;
        const y = work.y;
        const z = work.z;
        const trialX = work.trialX;
        const trialY = work.trialY;
        const trialZ = work.trialZ;
        const grad = work.grad;
        const precondition = work.precondition;
        const drift = work.drift;
        const dirX = cg.dirX;
        const dirY = cg.dirY;
        const dirZ = cg.dirZ;
        const prevDirX = cg.prevDirX;
        const prevDirY = cg.prevDirY;
        const prevDirZ = cg.prevDirZ;
        const prevGradX = cg.prevGradX;
        const prevGradY = cg.prevGradY;
        const prevGradZ = cg.prevGradZ;
        if (this.maxForce <= this.settings.gradientTol) {
            this.converged = true;
            this.stopReason = 'gradient';
            this.done = true;
            return false;
        }
        if (this.iterations >= this.settings.maxIterations) {
            this.stopReason = 'iterations';
            this.done = true;
            return false;
        }
        if (Date.now() - this.startTime > this.settings.timeBudgetMs) {
            this.stopReason = 'time';
            this.done = true;
            return false;
        }

        // Search direction: Polak-Ribiere conjugate gradient when a previous
        // direction exists, the plain gradient otherwise (and whenever the
        // conjugate direction is not a descent direction).
        let beta = 0.0;
        if (cg.valid) {
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
                cg.valid = false;
            }
        }
        let slope = 0.0;
        for (let i = 0; i < n; i++) {
            if (!mobile[i]) {
                continue;
            }
            dirX[i] = grad.gx[i] + (cg.valid ? beta * prevDirX[i] : 0.0);
            dirY[i] = grad.gy[i] + (cg.valid ? beta * prevDirY[i] : 0.0);
            dirZ[i] = grad.gz[i] + (cg.valid ? beta * prevDirZ[i] : 0.0);
            slope += dirX[i] * grad.gx[i] + dirY[i] * grad.gy[i] + dirZ[i] * grad.gz[i];
        }
        if (!(slope > 0.0)) {
            // Not a descent direction (or the structure has no mobile terms at
            // all): fall back to the gradient itself.
            cg.valid = false;
            for (let i = 0; i < n; i++) {
                dirX[i] = grad.gx[i];
                dirY[i] = grad.gy[i];
                dirZ[i] = grad.gz[i];
            }
        }

        // Trust radius: the *whole* preconditioned direction is scaled down so its
        // largest atom step is the cap, rather than clipping each atom on its own.
        // A per-atom clip would decouple the trial from the line-search factor —
        // every attempt would move the clipped atom by the same distance, so
        // halving the factor could not shrink the step, and a step whose curvature
        // cost exceeds its gradient gain would keep failing until the factor fell
        // below the clip threshold. Scaling the vector keeps the trial length
        // proportional to the factor, which is what makes the search able to walk
        // down at all.
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
        const trustScale = scaledMax > this.settings.maxStep ? this.settings.maxStep / scaledMax : 1.0;

        // A trial step: the preconditioned descent direction inside that trust
        // radius, line-searched by halving until the energy actually drops.
        let accepted = false;
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
                const dx = -work.lambda * trustScale * precondition[i] * dirX[i];
                const dy = -work.lambda * trustScale * precondition[i] * dirY[i];
                const dz = -work.lambda * trustScale * precondition[i] * dirZ[i];
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
            if (trialEnergy < this.energy) {
                this.energy = trialEnergy;
                accepted = true;
                work.lambda = Math.min(LAMBDA_MAX, work.lambda * 1.6);
                break;
            }
            work.lambda *= 0.5;
            if (work.lambda < LAMBDA_MIN) {
                break;
            }
        }
        if (!accepted) {
            if (cg.valid && !cg.resetUsed) {
                // The conjugate direction led nowhere; drop it, restore a step
                // size that a plain gradient step can use, and try once more.
                cg.resetUsed = true;
                cg.valid = false;
                work.lambda = LAMBDA_START;
                return true;
            }
            this.stopReason = 'step';
            this.done = true;
            return false;
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
        cg.valid = true;
        cg.resetUsed = false;
        // Rebuild the non-bonded pair list once an atom has drifted far enough
        // for a pair to have entered or left the range the list has to cover.
        // Built from the committed geometry, so every energy the next line
        // search compares comes from the same list.
        // The pair list is shared by every distance-dependent term, so it is
        // rebuilt whenever any of them is on.
        const anyDistance = this.settings.flags.nonbonded || this.settings.flags.electrostatic ||
            this.settings.flags.hbond || this.settings.flags.polarization;
        if (anyDistance && driftMax >= VDW_REBUILD_DRIFT) {
            field.nonbonded = buildNonbondedList(
                x, y, z, this.elements, buildNeighbourLists(this.bonds, n), mobile, n,
                this.settings.listCutoff);
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
        this.energy = evaluate(field, x, y, z, mobile, grad);
        this.readStats();
        this.iterations++;
        this.pushPoint();

        // A converging run keeps gaining less and less; stop when the gain has
        // been flat for a while rather than burning the whole step budget.
        if (this.lastEnergy - this.energy < ENERGY_PLATEAU) {
            this.plateau++;
        } else {
            this.plateau = 0;
        }
        this.lastEnergy = this.energy;
        if (this.plateau >= PLATEAU_STEPS) {
            this.stopReason = 'energy';
            this.done = true;
        }
        return true;
    }

    /**
     * The complete answer, valid at any time — after convergence, after a budget
     * ran out, or halfway through (which is what a pause reports). The geometry
     * returned is never a worse one than the input: every accepted step lowered
     * the energy, and the guard here enforces it even if that ever stopped
     * being true.
     */
    finish(): MmOptimizeResult {
        const useInitial = this.energy > this.energyBefore;
        const indices: number[] = [];
        const coords: number[][] = [];
        let maxMove = 0.0;
        for (let s = 0; s < this.indices.length; s++) {
            const i = this.indices[s];
            const px = useInitial ? this.inputX[i] : this.work.x[i];
            const py = useInitial ? this.inputY[i] : this.work.y[i];
            const pz = useInitial ? this.inputZ[i] : this.work.z[i];
            indices.push(i);
            coords.push([px, py, pz]);
            const mx = px - this.inputX[i];
            const my = py - this.inputY[i];
            const mz = pz - this.inputZ[i];
            const d = Math.sqrt(mx * mx + my * my + mz * mz);
            if (d > maxMove) {
                maxMove = d;
            }
        }
        const energy = useInitial ? this.energyBefore : this.energy;
        const stopText: { [key: string]: string } = {
            gradient: 'converged',
            energy: 'energy converged',
            iterations: 'step budget reached',
            time: 'time budget reached',
            step: 'no further gain',
        };
        const delta = energy - this.energyBefore;
        const message = 'Molecular Mechanics: ' + indices.length + ' atom(s) in ' + this.iterations +
            ' steps — E ' + this.energyBefore.toFixed(3) + ' → ' + energy.toFixed(3) +
            ' kcal/mol (Δ ' + delta.toFixed(3) + ', ' + (stopText[this.stopReason] || this.stopReason) + ')';
        return {
            ok: true,
            message: message,
            atomCount: this.atomCount,
            indices: indices,
            coords: coords,
            energyBefore: this.energyBefore,
            energyAfter: energy,
            iterations: this.iterations,
            maxGrad: this.maxForce,
            rmsGrad: this.rmsForce,
            maxMove: maxMove,
            converged: this.converged,
            stopReason: this.stopReason,
            counts: this.field.counts,
            invalid: this.invalidCount,
        };
    }
}

/**
 * Begin a run: validate the request, build the force field and take the first
 * energy and force reading. Either it refuses — with a complete result that
 * says why, so a caller has one shape to handle — or it returns the session to
 * step. Nothing is run here beyond that first evaluation, so the call is safe
 * to make from a UI thread.
 */
export function beginOptimization(
    rawAtoms: unknown,
    rawBonds: unknown,
    rawSelected: unknown,
    rawOptions?: unknown
): MmBeginResult {
    const settings = readOptions(rawOptions);
    const input = readInput(rawAtoms, rawBonds, rawSelected);
    const n = input.n;
    const refused: MmOptimizeResult = {
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
        refused.message = 'Molecular Mechanics: the structure has no atoms';
        return { session: null, result: refused };
    }
    if (input.selected.length === 0) {
        refused.message = 'Molecular Mechanics: select at least one atom to optimise';
        return { session: null, result: refused };
    }
    if (input.selected.length > MAX_MOBILE_ATOMS) {
        refused.message = 'Molecular Mechanics: ' + input.selected.length +
            ' atoms selected — the limit is ' + MAX_MOBILE_ATOMS + ' atoms per run';
        return { session: null, result: refused };
    }
    const session = new MmSession(input, settings);
    const trivial = session.trivialResult();
    if (trivial) {
        return { session: null, result: trivial };
    }
    return { session: session, result: null };
}

/**
 * Relax the selected atoms of a structure. Returns the new coordinates of the
 * selected atoms only — every unselected atom is left exactly where it was —
 * together with the energy before and after and the reason the run stopped.
 * This is the single-call form of the session API (begin → step → finish) and
 * the two can never disagree: they are the same implementation.
 */
export function optimizeSelection(
    rawAtoms: unknown,
    rawBonds: unknown,
    rawSelected: unknown,
    rawOptions?: unknown
): MmOptimizeResult {
    const began = beginOptimization(rawAtoms, rawBonds, rawSelected, rawOptions);
    if (!began.session) {
        return began.result as MmOptimizeResult;
    }
    began.session.step(Number.POSITIVE_INFINITY);
    return began.session.finish();
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
    value: number | boolean | string;
    label: string;
}

/** One settable parameter, rendered by the panel as a selection control. */
export interface MmControl {
    /** Key name in the run request (`MmOptimizeOptions`). */
    key: string;
    label: string;
    options: MmControlOption[];
    /** The value selected by default — the module's own default. */
    defaultValue: number | boolean | string;
}

/** What the panel shows before anything has been run. */
export interface MmFieldDescription {
    name: string;
    /** Parameter description lines, rendered in the panel exactly as they come. */
    lines: string[];
    /** The settings a caller may change, rendered as selection controls. */
    controls: MmControl[];
    /** The precision presets, each with the settings it stands for, so a panel
     *  can apply one without holding a second copy of the numbers. */
    presets: MmPreset[];
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
            'between ' + (DEFAULT_CUTOFF - DEFAULT_TAPER_WIDTH).toFixed(1) + ' and ' +
            DEFAULT_CUTOFF.toFixed(1) + ' Å at the default cutoff (settable, 5–12 Å)',
            'Light electrostatics: partial charges from the Pauling electronegativity difference ' +
            '(weighted by bond order, neutralised per structure), Coulomb with a distance-dependent ' +
            'dielectric εᵣ = r; hydrogen bond: a 10-12 well over a donor–H···acceptor contact ' +
            '(N/O/F) inside ≈1.6–2.8 Å and a D–H···A window of ≈120°–150°',
            'Metal coordination (high): ligand–metal–ligand angles toward the ideal geometry of the ' +
            'metal\'s coordination number; polarisation (high): a damped induced-dipole attraction ∝ 1/r⁴',
            'Restrained optimisation: only the selected atoms move, but every fixed atom still pulls on ' +
            'them through all of its terms (bonds, angles, torsions, impropers, non-bonded contacts).',
        ],
        controls: [
            {
                key: 'preset', label: 'Precision preset', defaultValue: 'medium',
                options: MM_PRESETS.map((p) => ({ value: p.key, label: p.label })).concat([
                    { value: 'custom', label: 'Custom — manual settings' },
                ]),
            },
            {
                key: 'nonbonded', label: 'Non-bonded', defaultValue: true,
                options: [
                    { value: true, label: 'on · 1-4 x0.5 · tapered cutoff' },
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
                key: 'electrostatic', label: 'Electrostatics', defaultValue: true,
                options: [
                    { value: true, label: 'on · εᵣ = r · 1-4 ×0.5' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'hbond', label: 'Hydrogen bond', defaultValue: true,
                options: [
                    { value: true, label: 'on · N/O/F donor–acceptor' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'coordination', label: 'Metal coordination', defaultValue: false,
                options: [
                    { value: true, label: 'on · geometry by coordination number' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'polarization', label: 'Polarisation', defaultValue: false,
                options: [
                    { value: true, label: 'on · damped induced dipole' },
                    { value: false, label: 'off' },
                ],
            },
            {
                key: 'gradientTol', label: 'Convergence max|F|', defaultValue: DEFAULT_GRADIENT_TOL,
                options: [
                    { value: DEFAULT_GRADIENT_TOL, label: DEFAULT_GRADIENT_TOL.toFixed(2) + ' kcal/mol/Å' + dflt },
                    { value: 0.01, label: '0.01 kcal/mol/Å (tight)' },
                    { value: 0.2, label: '0.2 kcal/mol/Å (loose)' },
                    { value: 0.4, label: '0.4 kcal/mol/Å (rough)' },
                ],
            },
            {
                key: 'maxIterations', label: 'Max steps', defaultValue: DEFAULT_MAX_ITERATIONS,
                options: [
                    { value: DEFAULT_MAX_ITERATIONS, label: String(DEFAULT_MAX_ITERATIONS) + dflt },
                    { value: 120, label: '120 (quick)' },
                    { value: 1000, label: '1000' },
                    { value: 2000, label: '2000 (careful)' },
                ],
            },
            {
                key: 'timeBudgetMs', label: 'Time budget', defaultValue: DEFAULT_TIME_BUDGET_MS,
                options: [
                    { value: DEFAULT_TIME_BUDGET_MS, label: (DEFAULT_TIME_BUDGET_MS / 1000) + ' s' + dflt },
                    { value: 1000, label: '1 s' },
                    { value: 1500, label: '1.5 s' },
                    { value: 2000, label: '2 s' },
                    { value: 10000, label: '10 s' },
                    { value: 20000, label: '20 s' },
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
            {
                key: 'cutoff', label: 'Non-bonded cutoff', defaultValue: DEFAULT_CUTOFF,
                options: [
                    { value: DEFAULT_CUTOFF, label: DEFAULT_CUTOFF.toFixed(0) + ' Å' + dflt },
                    { value: 5, label: '5 Å (tight)' },
                    { value: 6, label: '6 Å' },
                    { value: 10, label: '10 Å' },
                    { value: 12, label: '12 Å (wide)' },
                ],
            },
        ],
        presets: MM_PRESETS,
    };
}