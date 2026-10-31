import { Atom, Bond, MolecularData } from '../types';
import { ghostFromLabel } from '../tools/ghostAtoms';

/**
 * One atom label of an ORCA coordinate block.
 *
 * ORCA spells a dummy `DA` / `X` / `Xx` and a ghost `<element>:` (a colon right
 * after the symbol), so the ghost reading comes first; everything else keeps
 * the old rule — a one- or two-letter element symbol.
 */
function orcaAtomLabel(raw: string): { element: string; ghost: boolean } | null {
    const ghost = ghostFromLabel(raw);
    if (ghost) return { element: ghost.element, ghost: true };
    if (!/^[A-Za-z]{1,2}$/.test(raw)) return null;
    return { element: raw.charAt(0).toUpperCase() + raw.slice(1).toLowerCase(), ghost: false };
}

export function parseOrcainp(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge = 0;
    let multiplicity = 1;

    let inXyz = false;
    let inCoordsBlock = false;
    let inCoordsBohr = false;
    let atomIndex = 0;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        if (line.startsWith('#')) continue;

        if (line.startsWith('*')) {
            const lowerLine = line.toLowerCase();
            if (lowerLine.match(/\*\s*xyz/) || lowerLine.match(/\*\s*xyzfile/)) {
                inXyz = true;
                const parts = line.split(/\s+/);
                for (let p = 0; p < parts.length; p++) {
                    if (parts[p].toLowerCase() === 'xyz' || parts[p].toLowerCase() === 'xyzfile') {
                        if (p + 2 < parts.length) {
                            charge = parseInt(parts[p + 1], 10) || 0;
                            multiplicity = parseInt(parts[p + 2], 10) || 1;
                        }
                    }
                }
                continue;
            } else if (inXyz) {
                inXyz = false;
                continue;
            }
        }

        if (inXyz && line !== '') {
            const parts = line.split(/\s+/);
            if (parts.length >= 4) {
                const label = orcaAtomLabel(parts[0]);
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (label && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                    const atom: Atom = { element: label.element, x, y, z, index: atomIndex };
                    if (label.ghost) atom.ghost = true;
                    atoms.push(atom);
                    atomIndex++;
                }
            }
        }

        if (line.toLowerCase().startsWith('%coords')) {
            inCoordsBlock = true;
            inCoordsBohr = false;
            continue;
        }

        if (inCoordsBlock) {
            if (line.toLowerCase() === 'end') {
                inCoordsBlock = false;
                inCoordsBohr = false;
                continue;
            }
            const lowerLine = line.toLowerCase();
            if (lowerLine.startsWith('xyz')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 3) {
                    charge = parseInt(parts[1], 10) || 0;
                    multiplicity = parseInt(parts[2], 10) || 1;
                }
                continue;
            }
            if (lowerLine.startsWith('units')) {
                if (lowerLine.includes('bohr') || lowerLine.includes('angs')) {
                    inCoordsBohr = lowerLine.includes('bohr');
                }
                continue;
            }
            if (line !== '' && !line.startsWith('%')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 4) {
                    const label = orcaAtomLabel(parts[0]);
                    let x = parseFloat(parts[1]);
                    let y = parseFloat(parts[2]);
                    let z = parseFloat(parts[3]);
                    if (label && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        if (inCoordsBohr) {
                            const BOHR = 0.529177249;
                            x *= BOHR; y *= BOHR; z *= BOHR;
                        }
                        const atom: Atom = { element: label.element, x, y, z, index: atomIndex };
                        if (label.ghost) atom.ghost = true;
                        atoms.push(atom);
                        atomIndex++;
                    }
                }
            }
        }
    }

    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
