export interface Atom {
    element: string;
    x: number;
    y: number;
    z: number;
    index: number;
    occupancy?: number;
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

export interface OptStep {
    step: number;
    energy?: number;
    maxForce?: number;
    rmsForce?: number;
    maxDisplacement?: number;
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
    crystal?: CrystalData;
}
