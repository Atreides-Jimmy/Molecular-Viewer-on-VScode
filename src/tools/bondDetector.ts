import { Atom, Bond, MolecularData, isGhostAtom } from '../types';
import { perceiveBonds } from './bondPerception';
import { stripGhostBonds } from './ghostAtoms';

/**
 * Bond detection for structures that arrive without a bond table.
 *
 * The work itself lives in `bondPerception`, which is the single copy of the
 * logic the extension host and the webview both run — the webview evaluates
 * the same source rather than keeping a hand-maintained duplicate. This module
 * is the thin entry point the parsers and the extension call: it decides
 * whether a structure needs detecting at all, and forwards the file's total
 * charge when the format states one.
 *
 * A charge matters because the valence solver uses it: with the total charge
 * known it can accept or reject an assignment on the formal charges it
 * implies, which is what separates a neutral amine from an ammonium and a
 * nitro group from a nitrite. Without one the solver assumes neutrality, which
 * is right for most structures and is reported in the diagnostics when it is
 * not.
 */

/** Detect the bonds of a bare coordinate set. */
export function detectBonds(atoms: Atom[], charge?: number): Bond[] {
    return perceiveBonds(atoms, charge === undefined ? undefined : { charge: charge }).bonds;
}

/**
 * The bonds a file states *beside* its bond table, merged over the list it
 * belongs to: a stated bond the list already holds changes nothing (the
 * geometry's order stands), a stated bond it lacks is added, and a bond that
 * touches a ghost — or an atom that is not there — is dropped.
 *
 * The returned array is the input itself when nothing was added, so a structure
 * whose file states no extra bond is byte-identical to what it was.
 */
function mergeStatedBonds(atoms: Atom[], bonds: Bond[], stated: Bond[]): Bond[] {
    const known: { [key: string]: true } = {};
    for (const b of bonds) known[b.atom1 + '_' + b.atom2] = true;
    let added = 0;
    let merged: Bond[] | null = null;
    for (const s of stated) {
        const i1 = s.atom1;
        const i2 = s.atom2;
        if (!isFinite(i1) || !isFinite(i2) || i1 === i2) continue;
        if (i1 < 0 || i2 < 0 || i1 >= atoms.length || i2 >= atoms.length) continue;
        if (isGhostAtom(atoms[i1]) || isGhostAtom(atoms[i2])) continue;
        const lo = Math.min(i1, i2);
        const hi = Math.max(i1, i2);
        const key = lo + '_' + hi;
        if (known[key]) continue;
        known[key] = true;
        if (merged === null) merged = bonds.slice();
        merged.push({ atom1: lo, atom2: hi, order: s.order >= 1 ? s.order : 1 });
        added++;
    }
    return added > 0 && merged !== null ? merged : bonds;
}

/**
 * Give a structure its bonds unless its file already stated them.
 *
 * Authored bonds — a Gaussian `connect` section, a MOL2 bond block, PDB
 * `CONECT` records, a SMILES import — are never re-derived: the file said what
 * the bonds are, so the detector has nothing to add, and re-deriving them
 * could only lose information the parser had.
 *
 * A file may also state a bond without stating the structure's whole bond list
 * — a PDB `SSBOND` record is one — and `statedBonds` is that case: the bonds
 * are merged over whichever list this structure ends up with, so a protein
 * keeps its perceived bonds and gains the disulfide the file declares. The
 * field is cleared here, because after this pass the merged list is the truth.
 *
 * A bond that touches a ghost atom is dropped either way: the ghost has no
 * valence, so a file's bond table naming one (a `CONECT` record, a MOL2 bond
 * block) describes something the viewer cannot draw and the force field cannot
 * use. This is the one place every load path passes through, so no parser has
 * to remember it.
 */
export function ensureBonds(data: MolecularData): MolecularData {
    const stated = data.statedBonds && data.statedBonds.length > 0 ? data.statedBonds : undefined;
    const authored = stripGhostBonds(data.atoms, data.bonds);
    // `hasExplicitBonds` is the file's statement that its bond list is the
    // whole list, and an empty list is a statement too: a structure whose file
    // says it has no bonds (a MOL2 `@<TRIPOS>BOND` block with no records, an
    // empty `connect` section) must stay unbonded. Testing the list's length
    // here let that case fall through to perception, which invented bonds the
    // file had explicitly denied.
    if (data.hasExplicitBonds) {
        if (authored === data.bonds && stated === undefined) return data;
        const bonds = stated === undefined ? authored : mergeStatedBonds(data.atoms, authored, stated);
        return { ...data, bonds, statedBonds: undefined };
    }

    const perceived = detectBonds(data.atoms, data.charge);
    return {
        ...data,
        bonds: stated === undefined ? perceived : mergeStatedBonds(data.atoms, perceived, stated),
        hasExplicitBonds: false,
        statedBonds: undefined
    };
}
