import { Atom, Bond, MolecularData } from '../types';
import { elementFromZ, normalizeElementSymbol } from '../tools/periodicTable';

const ATOMIC_NUMBER_TO_ELEMENT: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn'
};

function parseMopacElement(token: string): string {
    const num = parseInt(token, 10);
    if (!isNaN(num) && num > 0 && ATOMIC_NUMBER_TO_ELEMENT[num]) {
        return ATOMIC_NUMBER_TO_ELEMENT[num];
    }
    const elem = token.charAt(0).toUpperCase() + token.slice(1).toLowerCase();
    if (/^[A-Za-z]{1,2}$/.test(elem)) return elem;
    return 'C';
}

interface MopacOutputRow {
    element: string;
    ghost: boolean;
    x: number;
    y: number;
    z: number;
}

/** One element label of a MOPAC output table (`C`, `Cl`, `X`, or an atomic number). */
function mopacOutputElement(token: string): { element: string; ghost: boolean } | undefined {
    const upper = token.toUpperCase();
    if (upper === 'X' || upper === 'XX') return { element: 'X', ghost: true };
    if (/^\d+$/.test(token)) {
        const element = elementFromZ(parseInt(token, 10));
        return element ? { element, ghost: false } : undefined;
    }
    const element = normalizeElementSymbol(token);
    return element ? { element, ghost: false } : undefined;
}

/**
 * One row of a MOPAC output coordinate table.
 *
 * MOPAC writes `1 C -4.1918 2.5617 1.3518` — the serial number comes before
 * the element label. Hand-written tables drop the serial, and an old table may
 * use the atomic number instead of the symbol; all three spellings are read.
 */
function mopacOutputRow(line: string): MopacOutputRow | undefined {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) return undefined;
    let at = 0;
    if (/^\d+$/.test(parts[0])) {
        // With a leading serial the element must follow it: treating the number
        // as the element would turn `1 2 3 4` into a hydrogen at (2, 3, 4).
        if (parts.length < 5 || !mopacOutputElement(parts[1])) return undefined;
        at = 1;
    }
    const element = mopacOutputElement(parts[at]);
    if (!element) return undefined;
    const x = parseFloat(parts[at + 1]);
    const y = parseFloat(parts[at + 2]);
    const z = parseFloat(parts[at + 3]);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return undefined;
    return { element: element.element, ghost: element.ghost, x, y, z };
}

export function parseMopac(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let hasExplicitBonds = false;
    let title = '';
    let charge: number | undefined;
    let multiplicity: number | undefined;
    let atomIndex = 0;

    let inCoords = false;
    let coordLines: string[] = [];
    let foundChargeMult = false;
    // Fallback strategy state (previously a separate full-file scan that ran
    // only when the primary strategy found no atoms). The first fallback
    // header's block is now captured inline during the same traversal; both
    // collectors see every line exactly as the two passes did.
    let fallbackPastHeader = false;
    const fallbackAtoms: Atom[] = [];

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        // --- Fallback block rows: captured without consuming the line, so
        // the primary strategy still processes it as before ---
        if (fallbackPastHeader && line !== '' && !line.startsWith('-') && !line.startsWith('=')) {
            const row = mopacOutputRow(line);
            if (row) {
                const atom: Atom = {
                    element: row.element, x: row.x, y: row.y, z: row.z,
                    index: fallbackAtoms.length
                };
                if (row.ghost) atom.ghost = true;
                fallbackAtoms.push(atom);
            }
        }

        if (!fallbackPastHeader) {
            if (line.match(/CARTESIAN COORDINATES/i) || line.match(/ATOM\s+X\s+Y\s+Z/i) || line.match(/FINAL GEOMETRY/i)) {
                fallbackPastHeader = true;
            }
        }

        if (i === 0) {
            // The first line is the keyword line (`PM7 CHARGE=1 MS=2`), not
            // the title; the title is the line after it. Reading the keywords
            // as the title lost the input's charge and multiplicity.
            const cmMatch = line.match(/\bCHARGE\s*=\s*(-?\d+)/i);
            if (cmMatch) charge = parseInt(cmMatch[1], 10);
            const msMatch = line.match(/\b(?:MULT|MS|MULTIPlicity)\s*=\s*(\d+)/i);
            if (msMatch) multiplicity = parseInt(msMatch[1], 10);
            continue;
        }

        if (i === 1) {
            title = line || 'MOPAC Calculation';
            continue;
        }

        // `foundChargeMult` only stops the charge scan; it must not stop the
        // coordinate-block detection, or a charge on the keyword line (read
        // above) would leave every geometry row unparsed.
        if (!foundChargeMult) {
            const cmMatch = line.match(/\bCHARGE\s*=\s*(-?\d+)/i);
            if (cmMatch) charge = parseInt(cmMatch[1], 10);

            const msMatch = line.match(/\b(?:MULT|MS|MULTIPlicity)\s*=\s*(\d+)/i);
            if (msMatch) multiplicity = parseInt(msMatch[1], 10);

            if (cmMatch || msMatch) {
                foundChargeMult = true;
            }
        }

        if (!inCoords &&
            (line.match(/^\s*[A-Za-z]{1,2}\s+/) || line.match(/^\s*\d+\s+/))) {
            inCoords = true;
        }

        if (inCoords) {
            if (line === '' || line.toLowerCase().startsWith('old') || line.toLowerCase().startsWith('===')) {
                if (coordLines.length > 0) {
                    // Primary block complete. The previous implementation broke
                    // out of the loop here; freezing these flags and continuing
                    // keeps the primary result identical while the fallback
                    // header search still runs inline.
                    inCoords = false;
                    foundChargeMult = true;
                }
                continue;
            }
            coordLines.push(line);
        }
    }

    for (const cl of coordLines) {
        const parts = cl.split(/\s+/);
        if (parts.length < 4) continue;

        const elem = parseMopacElement(parts[0]);
        const xVal = parseFloat(parts[1]);
        const yVal = parseFloat(parts[3]);
        const zVal = parseFloat(parts[5]);

        if (isNaN(xVal) || isNaN(yVal) || isNaN(zVal)) continue;
        let x: number, y: number, z: number;
        const xLabel = parts.length > 2 ? parts[2] : '';
        const yLabel = parts.length > 4 ? parts[4] : '';
        const zLabel = parts.length > 6 ? parts[6] : '';

        const isAngstrom = xLabel.toUpperCase() === '1' || xLabel.toUpperCase() === 'A' || xLabel === 'Å';
        const isBohr = xLabel.toUpperCase() === '0' || xLabel.toUpperCase() === 'B';

        if (isAngstrom) {
            x = xVal; y = yVal; z = zVal;
        } else if (isBohr) {
            x = xVal * 0.529177249;
            y = yVal * 0.529177249;
            z = zVal * 0.529177249;
        } else {
            x = xVal; y = yVal; z = zVal;
        }

        // MOPAC's dummy atom is `XX` (and `X` on newer versions); both resolve
        // to the element `X` here, which is the viewer's dummy spelling.
        const atom: Atom = { element: elem, x, y, z, index: atomIndex };
        if (elem === 'X' || elem === 'Xx' || elem === 'XX') {
            atom.element = 'X';
            atom.ghost = true;
        }
        atoms.push(atom);
        atomIndex++;
    }

    if (atoms.length === 0 && fallbackAtoms.length > 0) {
        for (const fa of fallbackAtoms) atoms.push(fa);
    }

    const data: MolecularData = { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
    if (atoms.length === 0 && coordLines.length > 0) {
        // MOPAC also accepts Z-matrix (internal-coordinate) input, whose rows
        // have distances and angles rather than x/y/z triples. Saying so is
        // better than an unexplained empty canvas.
        data.warnings = ['the MOPAC input states internal coordinates: only a Cartesian geometry can be shown'];
    }
    return data;
}

/**
 * Whether the text is a MOPAC calculation output rather than a MOPAC input.
 *
 * MOPAC inputs share `.dat`/`.mop` with data files and never print a version
 * banner; outputs always announce the program and always print at least one
 * `CARTESIAN COORDINATES` table.
 */
export function isMopacOutput(content: string): boolean {
    return /\bMOPAC\s*:?\s*(?:Version\s*)?\(?\s*(?:20\d\d|[6-9]\d?(?:\.\d+)?)/i.test(content) &&
        /^\s*CARTESIAN COORDINATES\b/im.test(content);
}

/**
 * Read the final geometry from a MOPAC output.
 *
 * Every geometry change is followed by a fresh `CARTESIAN COORDINATES` table,
 * so the last readable block is the structure the calculation ended with.
 * Rows are Angstrom unless the heading says Bohr. A MOPAC output that stops
 * before any table produces an explanatory warning instead of a blank canvas.
 */
export function parseMopacOutput(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    let atoms: Atom[] = [];
    for (let i = 0; i < lines.length; i++) {
        const header = lines[i];
        if (!/^\s*CARTESIAN COORDINATES\b/i.test(header)) continue;
        const scale = /BOHR/i.test(header) ? 0.529177249 : 1;
        const block: Atom[] = [];
        let skipped = 0;
        for (let j = i + 1; j < lines.length; j++) {
            const raw = lines[j];
            const line = raw.trim();
            if (line === '') continue;
            // The next geometry block starts here; the outer loop reads it.
            if (/^\s*CARTESIAN COORDINATES\b/i.test(raw)) break;
            const row = mopacOutputRow(line);
            if (row) {
                const atom: Atom = {
                    element: row.element,
                    x: row.x * scale,
                    y: row.y * scale,
                    z: row.z * scale,
                    index: block.length
                };
                if (row.ghost) atom.ghost = true;
                block.push(atom);
                continue;
            }
            // The column heading and a few banner lines sit between the
            // heading and the first row; after the table has started, the
            // first non-row line ends it.
            if (block.length > 0 || ++skipped > 8) break;
        }
        if (block.length > 0) atoms = block;
    }

    if (atoms.length === 0) {
        return {
            atoms: [], bonds: [], title: 'MOPAC output', hasExplicitBonds: false,
            warnings: ['the MOPAC output holds no readable Cartesian coordinate table']
        };
    }
    return { atoms, bonds: [], title: 'MOPAC output', hasExplicitBonds: false };
}
