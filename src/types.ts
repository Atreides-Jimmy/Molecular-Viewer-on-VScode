export interface Atom {
    element: string;
    x: number;
    y: number;
    z: number;
    index: number;
}

export interface Bond {
    atom1: number;
    atom2: number;
    order: number;
}

export interface GjfMeta {
    link0Lines: string[];
    routeLine: string;
    titleLines: string[];
    chargeMultLine: string;
    afterConnectContent: string;
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
}
