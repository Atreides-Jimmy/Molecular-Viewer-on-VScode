import { Atom, Bond, MolecularData } from '../types';
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

export function resolveElement(raw: string): string {
    const trimmed = raw.trim();
    const asNum = parseInt(trimmed, 10);
    if (!isNaN(asNum) && asNum > 0 && ATOMIC_NUMBER_MAP[asNum]) {
        return ATOMIC_NUMBER_MAP[asNum];
    }
    const cleaned = trimmed.replace(/[0-9]/g, '');
    const capitalized = cleaned.charAt(0).toUpperCase() + cleaned.slice(1).toLowerCase();
    if (capitalized.length >= 1 && /[A-Z][a-z]?/.test(capitalized)) {
        return capitalized;
    }
    return trimmed;
}

export function parseXyz(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];

    let atomCount = 0;
    if (lines.length > 0) {
        atomCount = parseInt(lines[0].trim(), 10) || 0;
    }

    const title = lines.length > 1 ? lines[1].trim() : '';

    const startLine = 2;
    let atomIndex = 0;
    for (let i = startLine; i < Math.min(startLine + atomCount, lines.length); i++) {
        const line = lines[i].trim();
        if (line === '') continue;

        const parts = line.split(/\s+/);
        if (parts.length >= 4) {
            // A ghost marker in the element field (`X`, `Bq`, `C-Bq`, `Du`, …)
            // is read before the ordinary element resolution: `resolveElement`
            // would turn `C-Bq` into the element "C-bq" and the atom would then
            // be a real, unknown, freely bonding centre.
            const ghost = ghostFromLabel(parts[0]);
            const element = ghost ? ghost.element : resolveElement(parts[0]);
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);

            if (element && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                const atom: Atom = { element, x, y, z, index: atomIndex };
                if (ghost) atom.ghost = true;
                atoms.push(atom);
                atomIndex++;
            }
        }
    }

    return { atoms, bonds, title, hasExplicitBonds: false };
}

/** Frames a multi-frame XYZ may contribute before the rest is ignored.
 *
 * A file is read whole into memory, and a structure file that a viewer can open
 * in one click can be a molecular-dynamics dump with hundreds of thousands of
 * frames — reading all of them would lock the window up for a trajectory nobody
 * steps through by hand. The limit is far above any quantum-chemistry
 * trajectory (an optimization, an IRC, a surface scan) and matches the
 * exposure a log file of the same length already has. */
const MAX_TRAJECTORY_FRAMES = 2000;

/** One structure of a multi-frame XYZ file. */
export interface XyzFrame {
    atoms: Atom[];
    title: string;
    /** Energy read from the comment line, in Hartree, when the comment states
     *  one in a recognizable form. */
    energy?: number;
}

export interface XyzTrajectory {
    frames: XyzFrame[];
    /** Comment of the first frame (`MD.title`). */
    title: string;
}

/**
 * Energy in an XYZ comment line.
 *
 * The trajectory files the usual programs write put the run's energy in the
 * comment — ORCA's `Coordinates from ORCA-job h2o E -76.321097057517`, xTB's
 * `energy: -17.6 gnorm: 0.001`, and the many tools that copy either style. The
 * marker must be a standalone `E`/`energy`/`Etot` (with an optional `=` or `:`)
 * followed by a signed decimal, so a comment that merely contains a number does
 * not become a curve.
 */
const COMMENT_ENERGY_RE = /(?:^|[\s,(])(?:E|Etot|Energy|TOTAL ENERGY)\s*[=:]?\s*(-?\d+\.\d+)/i;

function energyFromComment(comment: string): number | undefined {
    const m = COMMENT_ENERGY_RE.exec(comment);
    if (!m) return undefined;
    const v = parseFloat(m[1]);
    return isFinite(v) ? v : undefined;
}

/**
 * Parse a multi-frame XYZ file — the trajectory format ORCA writes as
 * `<job>_trj.xyz`, and the one xTB and many converters write as well.
 *
 * Frames are read back to back (`count`, comment, `count` atom rows, repeat).
 * Parsing stops at the first frame that does not fit that shape, and the frames
 * read so far are returned: a trajectory whose tail was cut off by a killed job
 * still shows every geometry that was written in full, instead of failing
 * outright. A single-frame file is returned as one frame — the caller decides
 * whether that is a trajectory or just a structure.
 */
export function parseXyzTrajectory(content: string): XyzTrajectory {
    const lines = content.split(/\r?\n/);
    const frames: XyzFrame[] = [];
    let i = 0;

    while (i < lines.length && frames.length < MAX_TRAJECTORY_FRAMES) {
        // Skip blank lines between frames (some writers separate them).
        while (i < lines.length && lines[i].trim() === '') i++;
        if (i >= lines.length) break;

        const count = parseInt(lines[i].trim(), 10);
        if (!isFinite(count) || count <= 0 || count > 1000000) break;
        const comment = (lines[i + 1] || '').trim();
        const start = i + 2;
        if (start + count > lines.length) break;

        const atoms: Atom[] = [];
        for (let k = 0; k < count; k++) {
            const row = lines[start + k].trim();
            const parts = row.split(/\s+/);
            if (parts.length < 4) { atoms.length = 0; break; }
            const ghost = ghostFromLabel(parts[0]);
            const element = ghost ? ghost.element : resolveElement(parts[0]);
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);
            if (!element || isNaN(x) || isNaN(y) || isNaN(z)) { atoms.length = 0; break; }
            const atom: Atom = { element, x, y, z, index: atoms.length };
            if (ghost) atom.ghost = true;
            atoms.push(atom);
        }
        // A frame that did not read completely ends the file; the partial frame
        // is dropped rather than shown as a truncated molecule.
        if (atoms.length !== count) break;

        frames.push({ atoms, title: comment, energy: energyFromComment(comment) });
        i = start + count;
    }

    return { frames, title: frames.length > 0 ? frames[0].title : '' };
}
