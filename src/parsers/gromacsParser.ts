import { MolecularData } from '../types';
import { elementFromAtomName, normalizeElementSymbol } from '../tools/periodicTable';
import { buildMolecularStructure, buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * GROMACS `.gro` reader.
 *
 * The format is fixed-width, positions are in nanometres, and a triclinic box
 * line lists its nine numbers in an order that is *not* row-major
 * (`v1x v2y v3z v1y v1z v2x v2z v3x v3y`). Reading that line as three rows
 * gives a cell with the tilt on the wrong axis, which is why the nine values
 * are reassembled explicitly here. A file whose box line is absent or all zero
 * is a molecule, not a 0 Angstrom crystal.
 */

const NM_TO_ANGSTROM = 10;

function cellFromBoxValues(box: number[]): number[][] | null {
    if (box.length === 0) return null;
    if (box.length >= 9) {
        // v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z)
        const a = [box[0], box[3], box[4]];
        const b = [box[5], box[1], box[6]];
        const c = [box[7], box[8], box[2]];
        return [a, b, c].map(v => v.map(x => x * NM_TO_ANGSTROM));
    }
    if (box.length >= 3) {
        return [
            [box[0] * NM_TO_ANGSTROM, 0, 0],
            [0, box[1] * NM_TO_ANGSTROM, 0],
            [0, 0, box[2] * NM_TO_ANGSTROM]
        ];
    }
    return null;
}

/** Float from a fixed-width slice, or undefined when the slice is not one. */
function sliceNumber(line: string, start: number, end: number): number | undefined {
    if (line.length < end) return undefined;
    const text = line.slice(start, end).trim();
    const v = parseFloat(text);
    return isFinite(v) ? v : undefined;
}

interface GroAtom {
    label: string;
    residueName: string;
    x: number;
    y: number;
    z: number;
}

function readAtomLines(lines: string[], count: number): { atoms: GroAtom[]; next: number } {
    const atoms: GroAtom[] = [];
    let i = 2;
    for (; i < lines.length && atoms.length < count; i++) {
        const line = lines[i];
        if (line.trim() === '') continue;
        let label = '';
        let residueName = '';
        let x: number | undefined;
        let y: number | undefined;
        let z: number | undefined;
        if (line.length >= 44) {
            label = line.slice(10, 15).trim();
            residueName = line.slice(5, 10).trim();
            x = sliceNumber(line, 20, 28);
            y = sliceNumber(line, 28, 36);
            z = sliceNumber(line, 36, 44);
        }
        if (label === '' || x === undefined || y === undefined || z === undefined) {
            // A writer that collapsed the padding still produces a readable
            // token line: a residue (number and name, possibly merged), an atom
            // name, an integer serial, then x y z. Finding the serial is what
            // keeps `1 SOL OW 1 x y z` from reading the serial as x.
            const parts = line.trim().split(/\s+/);
            label = '';
            x = undefined;
            y = undefined;
            z = undefined;
            for (let k = 2; k + 3 < parts.length; k++) {
                if (!/^\d+$/.test(parts[k])) continue;
                const px = parseFloat(parts[k + 1]);
                const py = parseFloat(parts[k + 2]);
                const pz = parseFloat(parts[k + 3]);
                if (!isFinite(px) || !isFinite(py) || !isFinite(pz)) continue;
                label = parts[k - 1];
                residueName = parts[k - 2];
                x = px;
                y = py;
                z = pz;
                break;
            }
            if (label === '' || x === undefined || y === undefined || z === undefined) {
                if (parts.length < 6) continue;
                label = parts[1];
                x = parseFloat(parts[3]);
                y = parseFloat(parts[4]);
                z = parseFloat(parts[5]);
                if (!isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
            }
        }
        atoms.push({
            label,
            residueName,
            x: x * NM_TO_ANGSTROM,
            y: y * NM_TO_ANGSTROM,
            z: z * NM_TO_ANGSTROM
        });
    }
    return { atoms, next: i };
}

function atomElement(atom: GroAtom): string {
    // The atom name is the best hint (`OW`, `HW1`, `NA`, `CL`, `Mg`); the
    // residue name is only a fallback because solvent residues are named `SOL`.
    return elementFromAtomName(atom.label) ||
        normalizeElementSymbol(atom.residueName) ||
        'X';
}

export function parseGro(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    if (lines.length < 2) {
        return {
            atoms: [], bonds: [], title: 'GROMACS gro', hasExplicitBonds: false,
            warnings: ['too few lines for a .gro file: no structure to show']
        };
    }
    const title = oneLineTitle(lines[0], 'GROMACS gro');
    const count = parseInt(lines[1].trim(), 10);
    if (!isFinite(count) || count < 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the atom count line of the .gro file is unreadable: no structure to show']
        };
    }
    const atomBlock = readAtomLines(lines, count);
    const atoms = atomBlock.atoms;
    if (atoms.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the .gro file lists no atoms: no structure to show']
        };
    }
    if (atoms.length < count) {
        warnings.push('the .gro file declares ' + count + ' atoms but only ' + atoms.length + ' are readable');
    }

    const sites: Site[] = atoms.map(a => {
        const element = atomElement(a);
        return { element, cart: [a.x, a.y, a.z], label: a.label };
    });
    const unknown = sites.filter(s => s.element === 'X').length;
    if (unknown > 0) warnings.push(unknown + ' atom' + (unknown === 1 ? '' : 's') + ' have no recognisable element name: read as X');

    // The box follows the atom block, not line `2 + count`: a file with a blank
// line among its atoms still reads the same rows, and its box would otherwise
// be one of the atoms.
    let boxIndex = atomBlock.next;
    while (boxIndex < lines.length && lines[boxIndex].trim() === '') boxIndex++;
    const boxLine = boxIndex < lines.length ? lines[boxIndex] : '';
    const boxValues = boxLine.trim().split(/\s+/).map(Number).filter(v => isFinite(v));
    const lattice = cellFromBoxValues(boxValues);
    if (lattice && isUsableLattice(lattice)) {
        return buildPeriodicStructure(lattice, sites, title, warnings);
    }
    if (lattice && boxValues.some(v => v !== 0)) {
        warnings.push('the .gro box is degenerate: the structure is read as a molecule');
    }
    return buildMolecularStructure(sites, title, warnings);
}

export function isGro(content: string): boolean {
    const lines = content.split(/\r?\n/);
    if (lines.length < 3) return false;
    const count = lines[1].trim();
    if (!/^\d+$/.test(count)) return false;
    // The fixed-width format runs the residue number and name together
    // (`    1SOL`), but a collapsed writer puts a space between them
    // (`1 SOL OW 1 ...`); both are the same residue/atom/serial/x-y-z row.
    return /^\s*\d+\S*(?:\s+\S+)?\s+\S+\s+\d+\s+[-+0-9.]+\s+[-+0-9.]+\s+[-+0-9.]+/.test(lines[2]);
}
