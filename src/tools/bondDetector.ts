import { Atom, Bond, MolecularData } from '../types';
import { perceiveBonds } from './bondPerception';
import { stripGhostBonds } from '../parsers/ghostAtoms';

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
 * Give a structure its bonds unless its file already stated them.
 *
 * Authored bonds — a Gaussian `connect` section, a MOL2 bond block, PDB
 * `CONECT` records, a SMILES import — are never re-derived: the file said what
 * the bonds are, so the detector has nothing to add, and re-deriving them
 * could only lose information the parser had.
 *
 * A bond that touches a ghost atom is dropped either way: the ghost has no
 * valence, so a file's bond table naming one (a `CONECT` record, a MOL2 bond
 * block) describes something the viewer cannot draw and the force field cannot
 * use. This is the one place every load path passes through, so no parser has
 * to remember it.
 */
export function ensureBonds(data: MolecularData): MolecularData {
    const authored = stripGhostBonds(data.atoms, data.bonds);
    if (data.hasExplicitBonds && data.bonds.length > 0) {
        if (authored === data.bonds) return data;
        return { ...data, bonds: authored };
    }

    const bonds = detectBonds(data.atoms, data.charge);
    return {
        ...data,
        bonds,
        hasExplicitBonds: false
    };
}
