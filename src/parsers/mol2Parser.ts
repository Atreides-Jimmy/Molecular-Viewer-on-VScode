import { Atom, Bond, MolecularData } from '../types';
import { ghostFromLabel } from '../tools/ghostAtoms';

export function parseMol2(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let title = '';
    let hasExplicitBonds = false;

    let i = 0;
    // The counter runs across the whole file, not per `@<TRIPOS>ATOM` block: a
    // MOL2 may hold several MOLECULE blocks, and restarting the numbering at
    // each block made two atoms share `index` and let a bond from the second
    // block address the first block's atoms.
    let atomIndex = 0;
    // MOL2 bonds name atoms by `atom_id` (the first column), which need not be
    // contiguous — a coordinate row this reader has to skip would otherwise
    // shift every later bond onto the wrong atom. The map is built as the atom
    // rows are accepted, so a bond resolves to the atom that actually exists.
    // It is block-local: a second MOLECULE block usually restarts its ids at 1,
    // and a file-wide map would resolve those bonds to the first block's atoms.
    // `atomIndex` alone stays continuous across blocks.
    let idToIndex: { [id: number]: number } = {};
    let unresolvedBonds = 0;

    while (i < lines.length) {
        const line = lines[i].trim();

        if (line.startsWith('@<TRIPOS>MOLECULE')) {
            i++;
            if (i < lines.length) {
                title = lines[i].trim();
            }
        } else if (line.startsWith('@<TRIPOS>ATOM')) {
            i++;
            idToIndex = {};
            while (i < lines.length && lines[i].trim() !== '' && !lines[i].trim().startsWith('@<TRIPOS>')) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 6) {
                    // Tripos writes a dummy atom with the type `Du`; the type
                    // field is also the only place a ghost can be named here,
                    // since a MOL2 atom line has no separate ghost marker.
                    const elementRaw = parts[5].split('.')[0];
                    const ghost = ghostFromLabel(elementRaw);
                    const element = ghost
                        ? ghost.element
                        : elementRaw.charAt(0).toUpperCase() + elementRaw.slice(1).toLowerCase();
                    const x = parseFloat(parts[2]);
                    const y = parseFloat(parts[3]);
                    const z = parseFloat(parts[4]);

                    if (element && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        const atom: Atom = { element, x, y, z, index: atomIndex };
                        if (ghost) atom.ghost = true;
                        atoms.push(atom);
                        const id = parseInt(parts[0], 10);
                        if (!isNaN(id) && idToIndex[id] === undefined) idToIndex[id] = atomIndex;
                        atomIndex++;
                    }
                }
                i++;
            }
            continue;
        } else if (line.startsWith('@<TRIPOS>BOND')) {
            hasExplicitBonds = true;
            i++;
            while (i < lines.length && lines[i].trim() !== '' && !lines[i].trim().startsWith('@<TRIPOS>')) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    const id1 = parseInt(parts[1], 10);
                    const id2 = parseInt(parts[2], 10);
                    const bondTypeStr = parts[3].toLowerCase();
                    let order = 1;
                    if (bondTypeStr === '2') {
                        order = 2;
                    } else if (bondTypeStr === 'ar') {
                        order = 1.5;
                    } else if (bondTypeStr === 'am') {
                        order = 1;
                    } else if (bondTypeStr === '3') {
                        order = 3;
                    } else if (bondTypeStr === 'du') {
                        order = 1;
                    } else {
                        const parsed = parseFloat(bondTypeStr);
                        if (!isNaN(parsed)) order = parsed;
                    }

                    const atom1 = id1 >= 0 ? idToIndex[id1] : undefined;
                    const atom2 = id2 >= 0 ? idToIndex[id2] : undefined;
                    if (atom1 === undefined || atom2 === undefined) {
                        // A bond whose id names no atom this file listed cannot
                        // be drawn; falling back to `id - 1` would silently join
                        // whatever atom now sits at that position. Count it so
                        // the file's own defect can be reported.
                        unresolvedBonds++;
                    } else if (atom1 !== atom2) {
                        const exists = bonds.some(b =>
                            (b.atom1 === atom1 && b.atom2 === atom2) ||
                            (b.atom1 === atom2 && b.atom2 === atom1)
                        );
                        if (!exists) {
                            bonds.push({ atom1, atom2, order });
                        }
                    }
                }
                i++;
            }
            continue;
        }

        i++;
    }

    const warnings: string[] = [];
    if (unresolvedBonds > 0) {
        warnings.push(unresolvedBonds + ' bond' + (unresolvedBonds === 1 ? '' : 's') +
            ' name an atom id the file does not list: those bonds were left out');
    }
    return warnings.length > 0
        ? { atoms, bonds, title, hasExplicitBonds, warnings }
        : { atoms, bonds, title, hasExplicitBonds };
}
