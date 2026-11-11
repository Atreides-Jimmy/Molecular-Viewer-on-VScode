// Payload shapes exchanged between the host and the webview page, and the
// element colours both sides use.
import { MolecularData, NormalMode, OptStep, OptThresholds } from '../types';
import { LogFrame, OrcaFrame, RouteSection, LogChargeSet, LogSummary } from '../parsers/index';

export const ATOM_COLORS: { [key: string]: string } = {
    H: '#FFFFFF', He: '#D9FFFF', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5',
    C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050', Ne: '#B3E3F5',
    Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6', Si: '#F0C8A0', P: '#FF8000',
    S: '#FFFF30', Cl: '#1FF01F', Ar: '#80D1E3', K: '#8F40D4', Ca: '#3DFF00',
    Sc: '#E6E6E6', Ti: '#BFC2C7', V: '#A6A6AB', Cr: '#8A99C7', Mn: '#9C7AC7',
    Fe: '#E06633', Co: '#F090A0', Ni: '#50D050', Cu: '#C88033', Zn: '#7D80B0',
    Ga: '#C28F8F', Ge: '#668F8F', As: '#BD80E3', Se: '#FFA100', Br: '#A62929',
    Kr: '#5CB8D1', Rb: '#702EB0', Sr: '#00FF00', Y: '#94FFFF', Zr: '#94E0E0',
    Nb: '#73C2C9', Mo: '#54B5B5', Tc: '#3B9E9E', Ru: '#248F8F', Rh: '#0A7D8C',
    Pd: '#006985', Ag: '#C0C0C0', Cd: '#FFD98F', In: '#A67573', Sn: '#668080',
    Sb: '#9E63B5', Te: '#D47A00', I: '#940094', Xe: '#429EB0', Cs: '#57178F',
    Ba: '#00C900', La: '#70D4FF', Ce: '#FFFFC7', Pr: '#D9FFC7', Nd: '#C7FFC7',
    Pm: '#A3FFC7', Sm: '#8FFFC7', Eu: '#61FFC7', Gd: '#45FFC7', Tb: '#30FFC7',
    Dy: '#1FFFC7', Ho: '#00FF9C', Er: '#00E675', Tm: '#00D452', Yb: '#00BF38',
    Lu: '#00AB24', Hf: '#4DC2FF', Ta: '#4DA6FF', W: '#2194D6', Re: '#267DAB',
    Os: '#266696', Ir: '#175487', Pt: '#D0D0E0', Au: '#FFD123', Hg: '#B8B8D0',
    Tl: '#A6544D', Pb: '#575961', Bi: '#9E4FB5', Po: '#AB5C00', At: '#754F45',
    Rn: '#428296', Fr: '#420066', Ra: '#007D00', Ac: '#70ABFA', Th: '#00BAFF',
    Pa: '#00A1FF', U: '#008FFF', Np: '#0080FF', Pu: '#006BFF', Am: '#545CF2',
    Cm: '#785CE3', Bk: '#8A4FE3', Cf: '#A136D4', Es: '#B31FD4', Fm: '#B31FBA',
    Md: '#B30DA6', No: '#BD0D87', Lr: '#C70066', Rf: '#CC0059', Db: '#D9004F',
    Sg: '#E00045', Bh: '#E6002E', Hs: '#EB0026',
    // Jmol's table ends at Hs and most copies of it stop at Mt. The remaining
    // ten follow the same transactinide ramp (Avogadro's continuation), so a
    // superheavy atom keeps a real element colour instead of the magenta
    // "unknown element" fallback - both when a file is opened and when one is
    // imported, which read the same table.
    Mt: '#E80026', Ds: '#EA0023', Rg: '#ED0021', Cn: '#EF001E', Nh: '#F2001C',
    Fl: '#F40019', Mc: '#F70016', Lv: '#F90014', Ts: '#FC0011', Og: '#FF000F'
};

/** 导入消息里的原子载荷 (客户端 importStructure 逐项写入 MD.atoms)。 */
export interface ImportAtomPayload {
    element: string;
    x: number;
    y: number;
    z: number;
    color: string;
    /** 该原子是 ghost/哑点 (见 Atom.ghost); 多份拷贝时也必须保留。 */
    ghost?: boolean;
}

/** 导入消息里的键载荷 (索引相对本次导入的原子数组)。 */
export interface ImportBondPayload {
    atom1: number;
    atom2: number;
    order: number;
}

/**
 * Payload of one molecule of a VMD visualization script: the structure the
 * script loads plus every drawing property the script's representations imply.
 * Everything a representation decides is resolved host-side (colours, sphere
 * radii, bond radii, surface patches) so the webview only has to draw it.
 */
export interface VmdAtomPayload {
    element: string;
    x: number;
    y: number;
    z: number;
    color: string;
    beta?: number;
    occupancy?: number;
    /** A ghost/dummy atom of the loaded structure (see `Atom.ghost`). */
    ghost?: boolean;
    /** Sphere size as a multiplier on the viewer's own radius; 0 = not drawn. */
    vmdScale?: number;
    /** van der Waals radius in Å — the reference length a `Surf` probe adds to. */
    vdw: number;
}

export interface VmdBondPayload {
    atom1: number;
    atom2: number;
    order: number;
    /** Cylinder radius in Å from the script; undefined = the viewer's own. */
    vmdR?: number;
}

/** One translucent surface patch (a `Surf` representation). */
export interface VmdSegmentPayload {
    color: string;
    indices: number[];
    probe: number;
    /** The egg this patch draws, when the structure sits next to an egg/unit file. */
    egg?: { name: string; units: number[]; level?: string };
}

export interface VmdRepPayload {
    style: string;
    colorMode: string;
    colorArg?: number;
    scaleMin?: number;
    scaleMax?: number;
    material: string;
    selection: string;
    /** Number of atoms this representation selects. */
    atoms: number;
}

export interface VmdMoleculePayload {
    label: string;
    /** Resolved absolute path of the structure file the script loads. */
    sourceFile: string;
    atoms: VmdAtomPayload[];
    bonds: VmdBondPayload[];
    segments: VmdSegmentPayload[];
    reps: VmdRepPayload[];
    warnings: string[];
    /** The molecule's own electronic state, as its file states it. */
    charge?: number;
    multiplicity?: number;
    /** Whether the molecule's bonds were authored rather than perceived. */
    hasExplicitBonds: boolean;
    /** The route and title card of a Gaussian input, for a later Save As. */
    gjfMeta?: MolecularData['gjfMeta'];
    /** Parser warnings of the molecule's file (distinct from `warnings`). */
    parseWarnings?: string[];
}

export interface VmdPayload {
    molecules: VmdMoleculePayload[];
    colorScaleMethod: string;
    warnings: string[];
}

/** Everything one structure file contributes, whichever way it was reached. */
export interface LoadedStructure {
    data: MolecularData;
    frames: (LogFrame | OrcaFrame)[];
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    routes?: RouteSection[];
    optThresholds?: OptThresholds;
    /** The population analyses and the job results a log states (see
     *  `LogChargeSet` / `LogSummary`); absent for every other format. */
    charges?: LogChargeSet[];
    summary?: LogSummary[];
    logSource?: string;
    error?: string;
}

/**
 * One frame of a binary trajectory, compact enough to inline in the page.
 *
 * The frames of a log or an XYZ carry a whole atom object list because they can
 * change composition. A trajectory cannot: every frame has the atoms of the
 * topology, so only the coordinates travel and the page pairs them with the
 * element template it already has. That keeps a 500-frame run from being
 * serialised as half a million element names.
 */
export interface TrajectoryFramePayload {
    stepLabel: string;
    /** `[x0, y0, z0, x1, …]`, Angstrom. */
    xyz: number[];
    cell?: number[][];
    time?: number;
    step?: number;
}

/** A trajectory and the template structure its frames are drawn with. */
export interface TrajectoryDocument {
    data: MolecularData;
    trajectoryFrames: TrajectoryFramePayload[];
}
