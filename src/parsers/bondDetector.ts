import { Atom, Bond, MolecularData } from '../types';
import { perceiveBonds } from './bondPerception';

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
 */
export function ensureBonds(data: MolecularData): MolecularData {
    if (data.hasExplicitBonds && data.bonds.length > 0) {
        return data;
    }

    const bonds = detectBonds(data.atoms, data.charge);
    return {
        ...data,
        bonds,
        hasExplicitBonds: false
    };
}
