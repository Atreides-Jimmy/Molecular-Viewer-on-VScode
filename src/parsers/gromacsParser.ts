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

/** The first three numbers of a text, whatever the field widths around them. */
function firstThreeNumbers(text: string): [number, number, number] | undefined {
    const re = /[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/g;
    const values: number[] = [];
    let m: RegExpExecArray | null;
    while ((m = re.exec(text)) !== null) {
        const v = parseFloat(m[0]);
        if (!isFinite(v)) continue;
        values.push(v);
        if (values.length === 3) return [values[0], values[1], values[2]];
    }
    return undefined;
}

/**
 * The three coordinate columns of a `.gro` atom line.
 *
 * The format is fixed-width, but the width of a coordinate field is not:
 * `-ndec 6` widens every field, and a value that overflows its own field
 * (`-1000.000`) swallows the space before the next one. Reading the fixed
 * 8-column slices therefore landed inside the wrong field for both, while
 * scanning *numbers* is independent of the widths and of the token count.
 *
 * The scan stays anchored on the serial, because a row may carry three
 * velocity columns after the positions (a `.gro` written by `mdrun`) and those
 * are not coordinates: the three numbers that follow the serial are the
 * position, whoever wrote the file.
 */
function coordinateTriple(line: string, fixedPrefix: boolean): [number, number, number] | undefined {
    if (fixedPrefix) {
        // The serial really is in its column, so the coordinates start right
        // after it — even when the first of them ran into the serial's field.
        const tail = firstThreeNumbers(line.slice(20));
        if (tail) return tail;
    }
    // A writer that collapsed the padding (or a prefix that overflowed its
    // columns) has no column 20 to cut at: the serial token is the anchor.
    const parts = line.trim().split(/\s+/);
    for (let k = 2; k + 3 < parts.length; k++) {
        if (!/^\d+$/.test(parts[k])) continue;
        const triple = firstThreeNumbers(parts.slice(k + 1).join(' '));
        if (triple) return triple;
    }
    return undefined;
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
        // The label block is fixed-width too, but it is only trusted where the
        // serial really is in its column: a residue number of six digits or an
        // atom name longer than five characters moves it, and a wrong label
        // would name the wrong element.
        const fixedPrefix = line.length >= 20 && /^\s*\d+\s*$/.test(line.slice(15, 20));
        if (fixedPrefix) {
            residueName = line.slice(5, 10).trim();
            label = line.slice(10, 15).trim();
        }
        if (label === '' || residueName === '') {
            // A writer that collapsed the padding still produces a readable
            // token line: a residue (number and name, possibly merged), an atom
            // name, an integer serial, then x y z. Finding the serial is what
            // keeps `1 SOL OW 1 x y z` from reading the serial as x.
            const parts = line.trim().split(/\s+/);
            label = '';
            residueName = '';
            for (let k = 2; k + 3 < parts.length; k++) {
                if (!/^\d+$/.test(parts[k])) continue;
                const px = parseFloat(parts[k + 1]);
                const py = parseFloat(parts[k + 2]);
                const pz = parseFloat(parts[k + 3]);
                if (!isFinite(px) || !isFinite(py) || !isFinite(pz)) continue;
                label = parts[k - 1];
                residueName = parts[k - 2];
                break;
            }
            if (label === '' || residueName === '') {
                if (parts.length < 6) continue;
                label = parts[1];
                residueName = parts[0];
            }
        }
        const coords = coordinateTriple(line, fixedPrefix);
        if (!coords) continue;
        atoms.push({
            label,
            residueName,
            x: coords[0] * NM_TO_ANGSTROM,
            y: coords[1] * NM_TO_ANGSTROM,
            z: coords[2] * NM_TO_ANGSTROM
        });
    }
    return { atoms, next: i };
}

/**
 * The monoatomic-ion residues whose names are not element spellings.
 *
 * CHARMM (and the GROMACS ports of its force fields) writes its ions as `SOD`,
 * `CLA`, `CAL`, `POT`, `CES`, `BAR` and `IOD`, and the topology's atom name is
 * the element the ion carries (`SOD  NA`, `CLA  CL`). Neither name is an
 * element symbol, so neither of the rules below can read them: a three-letter
 * residue fails the "residue name is the element" test, and `NA` then matches
 * the protein atom-name pattern and becomes nitrogen. A sodium drawn as a
 * nitrogen takes nitrogen's radius, colour and mass, and bonds to the wrong
 * neighbours. The table states the element the residue is, which is the only
 * place that fact exists.
 */
const ION_RESIDUE_ELEMENTS: { [residue: string]: string } = {
    SOD: 'Na', CLA: 'Cl', CAL: 'Ca', POT: 'K', CES: 'Cs', BAR: 'Ba', IOD: 'I'
};

function atomElement(atom: GroAtom): string {
    const label = atom.label.trim();
    const residue = atom.residueName.trim();

    // An ion residue is checked first: its atom name is a label for the ion,
    // which may well read like a protein atom (`NA`, `CL`). The table is only
    // trusted when the atom name agrees with it — it repeats the residue or
    // spells the ion's element — so a residue that merely happens to be called
    // `CAL` and carries a carbon is left to the ordinary rules.
    const ion = ION_RESIDUE_ELEMENTS[residue.toUpperCase()];
    if (ion && (label.toUpperCase() === residue.toUpperCase() || normalizeElementSymbol(label) === ion)) {
        return ion;
    }

    // A monoatomic ion's residue name is the element itself and its atom name
    // repeats it (`CA` in `CA`, `CL` in `CL`, `NA` in `NA`). The residue name
    // is the only place that can say so: `CA` is also a protein's C-alpha.
    if (/^[A-Za-z]{1,2}$/.test(residue) && label.toUpperCase() === residue.toUpperCase()) {
        const residueElement = normalizeElementSymbol(residue);
        if (residueElement) return residueElement;
    }

    // Topology atom names are labels, not element spellings: CA is C-alpha,
    // CD is C-delta, HG/HG1/HD11 is H-gamma. A name whose first letter is an
    // organic element takes that element; CL, BR, FE and ZN do not match this
    // shape and still resolve as their real symbols.
    if (/^[HCNOSP](?:[ABDEGZ][0-9]*|[0-9]+|')?$/.test(label)) {
        return label.charAt(0).toUpperCase();
    }

    // Everything else (OW, HW1, SIL, …) goes through the shared atom-name rules.
    return elementFromAtomName(label) ||
        normalizeElementSymbol(residue) ||
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

/**
 * A GROMACS run log. It shares `.log` with Gaussian/ORCA/xTB and with ABACUS,
 * but contains no coordinates: recognising it here keeps the log readers from
 * claiming it and lets the viewer say which file to open instead.
 */
export function isGromacsLog(content: string): boolean {
    return /^:-\)\s*GROMACS\s+-/m.test(content) ||
        /^GROMACS version:/im.test(content) ||
        /^\s*GROMACS\s*-\s*gmx\b/im.test(content);
}
