import { Atom, Bond, MolecularData, GjfMeta } from '../types';
import { ghostFromLabel } from '../tools/ghostAtoms';

const ATOMIC_NUMBER_MAP: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
};

function resolveElement(raw: string): string {
    let cleaned = raw.replace(/\(.*?\)/g, '');
    const asNum = parseInt(cleaned, 10);
    if (!isNaN(asNum) && asNum > 0 && ATOMIC_NUMBER_MAP[asNum]) {
        return ATOMIC_NUMBER_MAP[asNum];
    }
    cleaned = cleaned.replace(/[0-9]/g, '');
    const capitalized = cleaned.charAt(0).toUpperCase() + cleaned.slice(1).toLowerCase();
    if (capitalized.length >= 1 && capitalized.length <= 2 && /[A-Z][a-z]?/.test(capitalized)) {
        return capitalized;
    }
    return cleaned;
}

function parseAtomLine(parts: string[]): { element: string; x: number; y: number; z: number; ghost?: boolean } | null {
    if (parts.length < 4) return null;
    // Gaussian writes a ghost as `<element>-Bq` (and a bare dummy as `Bq`).
    // The marker has to be read before the element is resolved: the ordinary
    // path strips digits and would hand back the meaningless symbol `C-bq`.
    const ghost = ghostFromLabel(parts[0].replace(/\(.*?\)/g, ''));
    const element = ghost ? ghost.element : resolveElement(parts[0]);

    if (parts.length === 4) {
        const x = parseFloat(parts[1]);
        const y = parseFloat(parts[2]);
        const z = parseFloat(parts[3]);
        if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
            const out: { element: string; x: number; y: number; z: number; ghost?: boolean } = { element, x, y, z };
            if (ghost) out.ghost = true;
            return out;
        }
        return null;
    }

    // 5+ columns: extra columns are fixed-atom markers (0/-1/1).
    // Markers can appear after the element (e.g., "C -1 x y z") or after coordinates (e.g., "C x y z -1").
    // In Gaussian GJF, coordinates always have decimal points; markers are plain integers.
    // Strategy: filter out marker values (plain integers with value 0/-1/1), take first 3 remaining.
    const isMarker = (s: string): boolean => {
        if (s.includes('.') || /[eE]/.test(s)) return false; // has decimal point or exponent → coordinate
        const v = parseFloat(s);
        return !isNaN(v) && (v === 0 || v === -1 || v === 1);
    };

    const coords: number[] = [];
    for (let pi = 1; pi < parts.length; pi++) {
        if (isMarker(parts[pi])) continue;
        const v = parseFloat(parts[pi]);
        if (!isNaN(v)) coords.push(v);
    }

    if (coords.length >= 3) {
        const out: { element: string; x: number; y: number; z: number; ghost?: boolean } =
            { element, x: coords[0], y: coords[1], z: coords[2] };
        if (ghost) out.ghost = true;
        return out;
    }

    // Fallback: all values were plain integers — use positional logic.
    // If first value after element is 0/-1/1, treat as marker and skip.
    const allNums: number[] = [];
    for (let pi = 1; pi < parts.length; pi++) {
        const v = parseFloat(parts[pi]);
        if (!isNaN(v)) allNums.push(v);
    }
    if (allNums.length >= 4 && (allNums[0] === 0 || allNums[0] === -1 || allNums[0] === 1)) {
        const out: { element: string; x: number; y: number; z: number; ghost?: boolean } =
            { element, x: allNums[1], y: allNums[2], z: allNums[3] };
        if (ghost) out.ghost = true;
        return out;
    }
    if (allNums.length >= 3) {
        const out: { element: string; x: number; y: number; z: number; ghost?: boolean } =
            { element, x: allNums[0], y: allNums[1], z: allNums[2] };
        if (ghost) out.ghost = true;
        return out;
    }

    return null;
}

/**
 * Whether a line is the charge/multiplicity line rather than the first atom row.
 *
 * Gaussian writes the pair (or one pair per fragment) as plain integers, and
 * the atom block follows it. A file that omits the pair therefore used to lose
 * its first atom: the row was read as the pair, failed to parse as two
 * integers, and was consumed anyway. A line only counts here when every field
 * is a plain integer, the count is even, and every multiplicity field (the
 * odd-indexed one of each pair) is plausible. That rejects an atomic-number
 * atom row with integral coordinates such as `6 0 0 0` (multiplicity 0 and 0)
 * while still accepting a multi-fragment list such as `1 1 0 1`.
 */
function looksLikeChargeMultLine(line: string): boolean {
    const tokens = line.trim().split(/\s+/);
    if (tokens.length < 2 || tokens.length % 2 !== 0) return false;
    if (!tokens.every(t => /^[+-]?\d+$/.test(t))) return false;
    const values = tokens.map(t => parseInt(t, 10));
    if (values.some(v => !isFinite(v))) return false;
    for (let i = 0; i < values.length; i += 2) {
        if (Math.abs(values[i]) > 20) return false;
        if (values[i + 1] < 1 || values[i + 1] > 10) return false;
    }
    return true;
}

export function parseGjf(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge: number | undefined;
    let multiplicity: number | undefined;
    const warnings: string[] = [];

    const link0Lines: string[] = [];
    const routeLines: string[] = [];
    const titleLines: string[] = [];
    let chargeMultLine = '';
    let afterConnectContent = '';

    let i = 0;

    while (i < lines.length && lines[i].trim() !== '') {
        if (lines[i].trim().startsWith('%')) {
            link0Lines.push(lines[i]);
        } else {
            // A Gaussian route section may span several lines. Every one of
            // them is kept (not only the last), so saving the file back does
            // not drop the keywords an earlier line stated.
            routeLines.push(lines[i].trim());
        }
        i++;
    }
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }

    while (i < lines.length && lines[i].trim() !== '') {
        titleLines.push(lines[i]);
        i++;
    }
    title = titleLines.map(l => l.trim()).join(' ');
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }

    if (i < lines.length) {
        if (looksLikeChargeMultLine(lines[i])) {
            chargeMultLine = lines[i];
            const cmParts = lines[i].trim().split(/\s+/);
            const chrg = parseInt(cmParts[0], 10);
            const mult = parseInt(cmParts[1], 10);
            if (!isNaN(chrg)) charge = chrg;
            if (!isNaN(mult)) multiplicity = mult;
            i++;
        } else {
            // The atom block starts here; keeping the line for the loop below
            // is what preserves the first atom.
            warnings.push('the charge/multiplicity line is missing: the atom block was read from line ' + (i + 1));
        }
    }

    let atomIndex = 0;
    while (i < lines.length) {
        const line = lines[i].trim();
        if (line === '' || line.startsWith('--')) {
            break;
        }

        const parts = line.split(/\s+/);
        const result = parseAtomLine(parts);
        if (result) {
            const atom: Atom = { element: result.element, x: result.x, y: result.y, z: result.z, index: atomIndex };
            if (result.ghost) atom.ghost = true;
            atoms.push(atom);
            atomIndex++;
        }
        i++;
    }

    let connectStartLine = -1;
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }

    if (i < lines.length) {
        connectStartLine = i;
        let allNumeric = true;
        let lineCount = 0;

        for (let li = i; li < lines.length; li++) {
            const tl = lines[li].trim();
            if (tl === '' || tl.startsWith('--')) break;
            lineCount++;
            const tp = tl.split(/\s+/);
            const firstNum = parseInt(tp[0], 10);
            if (isNaN(firstNum)) { allNumeric = false; break; }
            for (let ti = 1; ti < tp.length; ti++) {
                if (isNaN(parseFloat(tp[ti]))) { allNumeric = false; break; }
            }
            if (!allNumeric) break;
        }

        if (allNumeric && lineCount > 0) {
            let connectEndLine = i;
            for (let li = i; li < lines.length; li++) {
                const tl = lines[li].trim();
                if (tl === '' || tl.startsWith('--')) { connectEndLine = li; break; }
                const cparts = tl.split(/\s+/);
                const atom1Num = parseInt(cparts[0], 10);
                if (atom1Num < 1 || atom1Num > atoms.length) {
                    connectEndLine = li + 1;
                    continue;
                }
                parseConnectLine(tl, atoms.length, bonds);
                connectEndLine = li + 1;
            }
            // `hasExplicitBonds` must state what the atom block actually got:
            // a line whose every bond index is out of range adds nothing, and
            // claiming an explicit table then makes the viewer skip its bond
            // perception and show a structure with no bonds at all.
            if (bonds.length > 0) hasExplicitBonds = true;
            let afterIdx = connectEndLine;
            while (afterIdx < lines.length && lines[afterIdx].trim() === '') {
                afterIdx++;
            }
            if (afterIdx < lines.length) {
                afterConnectContent = lines.slice(afterIdx).join('\n');
            }
        } else {
            afterConnectContent = lines.slice(i).join('\n');
        }
    }

    const gjfMeta: GjfMeta = {
        link0Lines,
        routeLine: routeLines.join('\n'),
        titleLines,
        chargeMultLine,
        afterConnectContent
    };

    // A GJF may hold several jobs joined by `--Link1--`. The atom block ends
    // at the first separator, so only the first structure is drawn; saying so
    // keeps the rest of the file from looking like a parsing failure.
    let linkedJobs = 0;
    for (let li = 0; li < lines.length; li++) {
        if (/^\s*--Link1--\s*$/i.test(lines[li])) linkedJobs++;
    }
    if (linkedJobs > 0) {
        warnings.push('the file holds ' + (linkedJobs + 1) +
            ' Gaussian jobs joined by --Link1--: only the first structure is shown (split the file to open another)');
    }
    const result: MolecularData = { atoms, bonds, title, hasExplicitBonds, gjfMeta, charge, multiplicity };
    if (warnings.length > 0) result.warnings = warnings;
    return result;
}

function parseConnectLine(line: string, totalAtoms: number, bonds: Bond[]): void {
    const cparts = line.split(/\s+/);
    if (cparts.length < 2) return;

    const atom1 = parseInt(cparts[0], 10) - 1;
    if (isNaN(atom1) || atom1 < 0 || atom1 >= totalAtoms) return;

    let j = 1;
    while (j + 1 < cparts.length) {
        const atom2 = parseInt(cparts[j], 10) - 1;
        const bondOrder = parseFloat(cparts[j + 1]) || 1;

        if (!isNaN(atom2) && atom2 >= 0 && atom2 < totalAtoms && atom1 !== atom2) {
            const exists = bonds.some(b =>
                (b.atom1 === atom1 && b.atom2 === atom2) ||
                (b.atom1 === atom2 && b.atom2 === atom1)
            );
            if (!exists) {
                const lo = Math.min(atom1, atom2), hi = Math.max(atom1, atom2);
                bonds.push({ atom1: lo, atom2: hi, order: bondOrder });
            }
        }
        j += 2;
    }

    if (cparts.length === 2) {
        const atom2 = parseInt(cparts[1], 10) - 1;
        if (!isNaN(atom2) && atom2 >= 0 && atom2 < totalAtoms && atom1 !== atom2) {
            const exists = bonds.some(b =>
                (b.atom1 === atom1 && b.atom2 === atom2) ||
                (b.atom1 === atom2 && b.atom2 === atom1)
            );
            if (!exists) {
                const lo2 = Math.min(atom1, atom2), hi2 = Math.max(atom1, atom2);
                bonds.push({ atom1: lo2, atom2: hi2, order: 1 });
            }
        }
    }
}
