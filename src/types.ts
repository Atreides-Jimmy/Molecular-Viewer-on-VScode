export interface Atom {
    element: string;
    x: number;
    y: number;
    z: number;
    index: number;
    /**
     * A ghost atom: the coordinates of an atom without a nucleus, electrons or
     * bonds. `element` still names the atom the ghost stands for, because that
     * is what a counterpoise calculation needs (`C-Bq` in Gaussian, `C:` in
     * ORCA); a bare dummy point has `element === 'X'`.
     *
     * A ghost never carries a bond, never counts in a formula or an electron
     * count and is skipped by bond perception — every one of those rules is
     * keyed off `isGhostAtom`, so a structure cannot end up with a ghost that
     * some pass treats as a real atom and another does not.
     */
    ghost?: boolean;
    occupancy?: number;
    // PDB columns 61-66. Not a B-factor: the visualization toolchain writes an
    // arbitrary per-atom scalar (unit index, force error, …) into that column,
    // so it is kept as a plain number the viewer can colour by.
    beta?: number;
    baseIdx?: number;
    cellI?: number;
    cellJ?: number;
    cellK?: number;
    /**
     * The atom-site label the file gave the asymmetric-unit site this atom was
     * expanded from (`C12`, `O1W`). Only the CIF parser writes it — the other
     * formats have no labels — and it travels on the atom, so the passes that
     * rebuild `baseAtoms` (delete, disorder removal) keep site and label
     * together instead of relying on a parallel array staying aligned.
     */
    siteLabel?: string;
    /**
     * The asymmetric-unit site this cell atom is an image of (`asymIdx`, an
     * index into `CrystalData.asymAtoms`) and the operation that placed it
     * (`symOp`, an index into `CrystalData.symmetryOps`). The pair is what
     * answers "where did this atom come from" for a tooltip or a site table, and
     * both parsers write it while expanding, so it is exact rather than
     * reconstructed by re-running the expansion. Absent on a structure that
     * never had an asymmetric unit to begin with.
     */
    asymIdx?: number;
    symOp?: number;
    /**
     * POSCAR "Selective dynamics" flags for this atom, one letter per axis
     * (`T T F`). Kept verbatim so a read-save-read round trip states exactly
     * what the file stated; the viewer itself has no use for the constraint.
     */
    constraint?: string;
    /**
     * The file's disorder group (`_atom_site_disorder_group`), when it states
     * one. Empty/undefined means "not part of any group". The page uses it to
     * show one component of a disordered site at a time; it is never inferred
     * from occupancy, because two groups can each be fully occupied in
     * different parts of a structure and occupancy alone cannot tell them apart.
     */
    disorderGroup?: string;
}

export interface Bond {
    atom1: number;
    atom2: number;
    order: number;
    crossCell?: boolean;
    shift?: [number, number, number];
    /**
     * `_geom_bond_valence` as the file states it, when it states one. A bond
     * valence is not a bond order — it is the share of an atom's valence this
     * bond carries, so a single bond is 1 and a delocalised one can be 1.5 —
     * and it is kept verbatim so a read-save-read round trip says what the
     * file said instead of a number the viewer guessed.
     */
    valence?: number;
}

/**
 * One atom site exactly as the file lists it — the asymmetric unit — together
 * with the symmetry operations that actually placed an image of it.
 *
 * `baseAtoms` is the *expanded* cell: every site multiplied by every operation
 * and deduplicated. That is what the viewer draws, but it cannot answer "what
 * does the file itself contain", which is the question a save-back to CIF has
 * to answer (writing the expanded cell under the full operation list would
 * expand it again on the next read). Keeping the original sites alongside makes
 * that answerable without re-reading the file.
 *
 * Nothing but the parser writes this. Every other pass may read it; a pass that
 * edits coordinates must keep editing `baseAtoms`, because the asymmetric unit
 * is a record of the file, not of the session.
 */
export interface AsymAtom {
    element: string;
    /** Site label as written (`C12`, `O1W`); empty when the file has none. */
    label: string;
    /** Cartesian position of the site (the identity image), in Angstrom. */
    x: number;
    y: number;
    z: number;
    /** Fractional position as listed, wrapped into [0, 1). */
    frac: [number, number, number];
    /** Position in `asymAtoms` — the index an `asymBonds` entry refers to. */
    index: number;
    occupancy?: number;
    /**
     * Indices into `CrystalData.symmetryOps` whose image of this site survived
     * deduplication, in increasing order. A site on a special position lists
     * fewer operations than the file declares; a site whose images all landed
     * on atoms already placed (a duplicate site) lists none. Identity is index 0
     * only when the file happens to list it first.
     */
    symmetryOps: number[];
    /** `_atom_site_disorder_group` for this site, when the file states one. */
    disorderGroup?: string;
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
    /**
     * Where the operations came from: the file's own symmetry table, or the
     * built-in space-group table consulted because the file stated only a
     * space-group name. Absent on formats that carry no operations at all.
     */
    symmetrySource?: 'file' | 'builtin';
    baseAtoms: Atom[];
    baseBonds: Bond[];
    /**
     * The asymmetric unit the file lists (see `AsymAtom`). Optional because
     * only the formats that state symmetry fill it: a POSCAR or a VESTA file is
     * already one cell with no operations, and a cube file has no cell at all.
     */
    asymAtoms?: AsymAtom[];
    /**
     * Bonds among `asymAtoms`, by their `index`. Present only when a format
     * states bonds for the asymmetric unit itself; a CIF's `_geom_bond` loop is
     * expanded into the base cell instead (see `fileBonds`), so this stays
     * undefined rather than being invented.
     */
    asymBonds?: Bond[];
    /**
     * The bonds the file itself states (`_geom_bond`), already expanded into
     * `baseAtoms` indices with the periodic image each bond points at. Only the
     * CIF reader fills it, because only CIF carries a bond table among the
     * formats read here.
     *
     * This is evidence, not the session's bond list: `baseBonds` stays the
     * geometry-perceived list, and the page merges this table over it once, on
     * the first rebuild. The table states which of the detected bonds the file
     * agrees with and with what order; a bond the table lists that perception
     * did not find is reported, not drawn, because a single bad row must not
     * replace the picture the structure's geometry supports.
     *
     * An entry's `order` is 0 when the table lists the bond but states no order
     * (`_geom_bond_order` absent or missing for that row): the file has said
     * nothing about the order, and the perceived one stands.
     */
    fileBonds?: Bond[];
    /**
     * How many whole molecules (finite connected components) the unit cell
     * contains — the crystallographic Z for a molecular crystal. Written only
     * by the page's `rebuildCrystal` as a session observation on its own copy;
     * it is never part of the payload and the host never reads it. A cell whose
     * only component is an infinite network (a framework, diamond) stays
     * without a meaningful Z, which the page shows by omitting the line.
     */
    moleculeCount?: number;
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
    /**
     * Things the parser had to work around, in the file's own terms ("no
     * symmetry operations in the file: read as P1"). The viewer shows them as
     * one line in the structure panel; they are never a substitute for the
     * structure itself, only an explanation of what it is.
     *
     * A parser reports a problem here instead of throwing or silently
     * substituting a default: a CIF whose cell is missing comes back as an
     * empty structure with a warning, not as a structure with a 1 A cube.
     */
    warnings?: string[];
    /** Present when the structure was loaded through a VMD visualization script. */
    vmd?: VmdScene;
    crystal?: CrystalData;
}

/**
 * Whether an atom is a ghost (a position with no nucleus, electrons or bonds).
 *
 * Two spellings mean the same thing, and every pass in the extension and in the
 * page asks this one question rather than testing either: the `ghost` flag, and
 * the dummy element `X` that every file format without a ghost concept writes
 * (`X`, `Bq`, `Du`, `DA`, `Xx`, `q` — the parsers all normalise to `X`). The
 * page carries the same predicate in its own script (the viewer's client code
 * is inlined, so it cannot import this one); the two are kept identical.
 */
export function isGhostAtom(a: { element?: unknown; ghost?: unknown } | null | undefined): boolean {
    if (!a) return false;
    if (a.ghost === true) return true;
    return a.element === 'X';
}
