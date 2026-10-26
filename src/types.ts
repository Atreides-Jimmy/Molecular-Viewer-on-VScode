export interface Atom {
    element: string;
    x: number;
    y: number;
    z: number;
    index: number;
    occupancy?: number;
    // PDB columns 61-66. Not a B-factor: the visualization toolchain writes an
    // arbitrary per-atom scalar (unit index, force error, …) into that column,
    // so it is kept as a plain number the viewer can colour by.
    beta?: number;
    baseIdx?: number;
    cellI?: number;
    cellJ?: number;
    cellK?: number;
}

export interface Bond {
    atom1: number;
    atom2: number;
    order: number;
    crossCell?: boolean;
    shift?: [number, number, number];
}

export interface CrystalData {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
    latticeVectors: number[][];
    spaceGroup?: string;
    symmetryOps: string[];
    baseAtoms: Atom[];
    baseBonds: Bond[];
}

export interface GjfMeta {
    link0Lines: string[];
    routeLine: string;
    titleLines: string[];
    chargeMultLine: string;
    afterConnectContent: string;
}

export interface AtomGroup {
    colorId: number;
    color: string;
    indices: number[];
}

/** One `mol addrep` block of a VMD visualization script. */
export interface VmdRep {
    /** `mol representation <style>` — CPK, Surf, Licorice, VDW, Lines, … */
    style: string;
    /** Numeric arguments of the representation (sphere scale, bond radius, …). */
    styleArgs: number[];
    /** `mol color <mode>` — Beta, ColorID, Name, Element, … */
    colorMode: string;
    /** Argument of the color mode (the ColorID index when colorMode is ColorID). */
    colorArg?: number;
    /** Raw `mol selection` expression. */
    selection: string;
    /** `mol material` — Opaque, Transparent, … */
    material: string;
    /** `mol scaleminmax <mol> <rep> <min> <max>` colour-scale range, if written. */
    scaleMin?: number;
    scaleMax?: number;
}

/** One `mol new` molecule of a VMD visualization script. */
export interface VmdMolecule {
    /** File the script loads (`mol new <file>`). */
    sourceFile: string;
    /** `mol new … type <t>`, empty when the script relies on the extension. */
    sourceType: string;
    /** Display label (source file name, plus the molecule ordinal). */
    label: string;
    reps: VmdRep[];
}

/** Everything a VMD visualization script says about how to draw the structure. */
export interface VmdScene {
    molecules: VmdMolecule[];
    /** `color Display Background <name>` — a colour name or #rrggbb. */
    background?: string;
    /** `color scale method <name>` — RGB is the only one implemented. */
    colorScaleMethod: string;
    /** `color change rgb <id> <r> <g> <b>` overrides, keyed by colour id. */
    colorOverrides?: { [id: number]: string };
}

export interface OptStep {
    step: number;
    energy?: number;
    maxForce?: number;
    rmsForce?: number;
    maxDisplacement?: number;
    rmsDisplacement?: number;
}

/**
 * The convergence targets the program itself states in its output.
 *
 * The convergence panel draws them as dashed target lines. They are read from
 * the file instead of being hard-coded per program, because the same program
 * can be told to converge tighter (ORCA's `TolMaxG`/`TolRMSD`, a Gaussian
 * `opt=(tight)` route) and a chart that drew the default tolerances would then
 * show a curve crossing the wrong line. A program that states none gets no
 * lines at all rather than a borrowed scale.
 */
export interface OptThresholds {
    /** Hartree/Bohr — largest single force component allowed at convergence. */
    maxForce?: number;
    /** Hartree/Bohr — RMS force allowed at convergence. */
    rmsForce?: number;
    /** Bohr — largest single displacement of the last step. */
    maxDisplacement?: number;
    /** Bohr — RMS displacement of the last step. */
    rmsDisplacement?: number;
}

export interface NormalMode {
    index: number;
    frequency: number;
    symmetry?: string;
    reducedMass?: number;
    forceConstant?: number;
    irIntensity?: number;
    displacements: number[][];
}

export interface MolecularData {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    filePath?: string;
    gjfMeta?: GjfMeta;
    charge?: number;
    multiplicity?: number;
    atomGroups?: AtomGroup[];
    /** Present when the structure was loaded through a VMD visualization script. */
    vmd?: VmdScene;
    crystal?: CrystalData;
}
