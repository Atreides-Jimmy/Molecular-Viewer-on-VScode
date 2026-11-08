import { Atom, Bond, MolecularData } from '../types';
import { ghostFromLabel } from '../tools/ghostAtoms';
import { normalizeElementSymbol } from '../tools/periodicTable';

/**
 * MDL CTfile (`.mol`) — V2000 and V3000.
 *
 * Save As writes MDL Mol, so a saved file must be openable by the same program
 * that wrote it. V2000 was never read here (a `.mol` fell through to the XYZ
 * reader and came back as an empty structure), and the writer's own output is
 * not the only V2000 a user will meet: RDKit, OpenBabel, ChemDraw and Marvin
 * all write one with three header lines, fixed-width counts and a three-column
 * atom/bond table.
 *
 * The counts line is located by its version marker rather than assumed to be
 * line four. Files written before this reader included only a title before the
 * counts line, and a reader that hard-codes the line number would lose the
 * whole file; the marker (`... V2000` / `V3000`) is what every writer states.
 *
 * V3000 states atom IDs that need not be contiguous, so bonds are resolved
 * through an explicit id -> index map rather than by subtracting one.
 */

interface MolElement {
    element: string;
    ghost?: boolean;
}

/** `*`, `LP` and the other spellings a molfile uses for a dummy point. */
function molElement(symbol: string): MolElement | undefined {
    const text = String(symbol == null ? '' : symbol).trim();
    if (text === '') return undefined;
    // V2000's wildcard atom (and V3000's) is `*`: a point with no element.
    if (text === '*') return { element: 'X', ghost: true };
    // A lone pair is drawn as a point, not an atom.
    if (text.toUpperCase() === 'LP') return { element: 'X', ghost: true };
    // Deuterium and tritium are hydrogen in every structure file that names
    // them by isotope; the viewer has no separate isotope state.
    if (/^(D|T)$/i.test(text)) return { element: 'H' };
    const ghost = ghostFromLabel(text);
    if (ghost) return { element: ghost.element, ghost: true };
    const element = normalizeElementSymbol(text);
    if (element) return { element: element };
    // A query atom (`A`, `Q`, `L`, `R#`, …) names no element. Writing it as the
    // dummy point is honest; calling it carbon would put a non-existent atom
    // into the formula and into bond perception.
    return { element: 'X', ghost: true };
}

function findCountsLine(lines: string[]): number {
    // Three header lines are the format's own, but the writer shipped here used
    // one for a while, and deposited files exist without the program/comment
    // lines. The version marker is the reliable anchor, and a few extra lines
    // of preamble cost nothing.
    const limit = Math.min(lines.length, 24);
    for (let i = 0; i < limit; i++) {
        const line = lines[i];
        if (!/\bV2000\b/i.test(line) && !/\bV3000\b/i.test(line)) continue;
        if (!/^\s*\d/.test(line)) continue;
        return i;
    }
    return -1;
}

/** Whether the text is an MDL molfile this reader understands. */
export function isMdlMol(content: string): boolean {
    const lines = String(content || '').replace(/^\uFEFF/, '').split(/\r?\n/);
    if (findCountsLine(lines) < 0) return false;
    // A CTfile ends with `M  END`; the marker keeps a stray "V2000" in an
    // unrelated output from claiming the file.
    return /\bM\s*END\b/i.test(content);
}

function v2000Counts(line: string): { atoms: number; bonds: number } {
    // Padded (` 12 13  0…`) and merely spaced (`12 13 0…`) both give the two
    // counts as the leading integers; the fixed-width slice would read a
    // leading-blank file as `12` and ` 1`.
    const m = /(\d+)\s+(\d+)/.exec(line);
    if (!m) return { atoms: 0, bonds: 0 };
    return { atoms: parseInt(m[1], 10), bonds: parseInt(m[2], 10) };
}

interface MolAtomRow {
    element: MolElement;
    x: number;
    y: number;
    z: number;
}

/** One V2000 atom line: fixed columns first, whitespace tokens as a fallback. */
function v2000AtomRow(line: string): MolAtomRow | undefined {
    const x = parseFloat(line.substring(0, 10));
    const y = parseFloat(line.substring(10, 20));
    const z = parseFloat(line.substring(20, 30));
    const fixedSymbol = line.substring(31, 34).trim();
    // The symbol field must look like one. A space-separated line can leave a
    // digit in that column, and accepting it would turn the atom into a ghost
    // at coordinates read from the wrong fields instead of taking the tokens.
    if (isFinite(x) && isFinite(y) && isFinite(z) && fixedSymbol !== '' && !/[0-9]/.test(fixedSymbol)) {
        const element = molElement(fixedSymbol);
        if (element) return { element: element, x, y, z };
    }
    // A hand-written or reformatted file may drop the fixed columns. The
    // standard order is x y z symbol; the symbol-first form exists in old
    // listings and is accepted when the token positions say so.
    const parts = line.trim().split(/\s+/);
    if (parts.length >= 4) {
        const coords = [parseFloat(parts[0]), parseFloat(parts[1]), parseFloat(parts[2])];
        if (coords.every(v => isFinite(v))) {
            const element = molElement(parts[3]);
            if (element) return { element: element, x: coords[0], y: coords[1], z: coords[2] };
        }
        const symbol = molElement(parts[0]);
        const rest = [parseFloat(parts[1]), parseFloat(parts[2]), parseFloat(parts[3])];
        if (symbol && rest.every(v => isFinite(v))) {
            return { element: symbol, x: rest[0], y: rest[1], z: rest[2] };
        }
    }
    return undefined;
}

/** V2000 bond type -> the viewer's bond order (4 is aromatic). */
function v2000BondOrder(type: number): number {
    if (type === 2) return 2;
    if (type === 3) return 3;
    if (type === 4) return 1.5;
    if (type === 7) return 2;
    // 5, 6 and 8 are query/delocalised single bonds; a single stick is what
    // their geometry supports.
    return 1;
}

interface MolBondRow {
    atom1: number;
    atom2: number;
    order: number;
}

/** One V2000 bond line, again with the fixed columns preferred. */
function v2000BondRow(line: string): MolBondRow | undefined {
    const f1 = line.substring(0, 3);
    const f2 = line.substring(3, 6);
    const f3 = line.substring(6, 9);
    // Each fixed field is a right-justified integer. Checking the field shape
    // keeps a space-separated line (`1 2 1 0`, whose first nine characters would
    // also parse as numbers) from being read as a 1-1 bond.
    const integerField = (f: string) => /^\s*\d+\s*$/.test(f);
    if (integerField(f1) && integerField(f2) && integerField(f3)) {
        return { atom1: parseInt(f1, 10), atom2: parseInt(f2, 10), order: v2000BondOrder(parseInt(f3, 10)) };
    }
    const parts = line.trim().split(/\s+/);
    if (parts.length < 3) return undefined;
    const p1 = parseInt(parts[0], 10);
    const p2 = parseInt(parts[1], 10);
    const p3 = parseInt(parts[2], 10);
    if (!isFinite(p1) || !isFinite(p2) || !isFinite(p3)) return undefined;
    return { atom1: p1, atom2: p2, order: v2000BondOrder(p3) };
}

function pushBond(bonds: Bond[], seen: Set<string>, atom1: number, atom2: number, order: number): boolean {
    if (!isFinite(atom1) || !isFinite(atom2)) return false;
    if (atom1 < 0 || atom2 < 0 || atom1 === atom2) return false;
    const lo = atom1 < atom2 ? atom1 : atom2;
    const hi = atom1 < atom2 ? atom2 : atom1;
    const key = lo + '-' + hi;
    if (seen.has(key)) return false;
    seen.add(key);
    bonds.push({ atom1, atom2, order });
    return true;
}

function parseV2000(lines: string[], countsAt: number, title: string): MolecularData {
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    const warnings: string[] = [];
    const counts = v2000Counts(lines[countsAt]);
    // A row that could not be read must not shift the numbering of the rows
    // after it: bonds address the row position in the file, so every declared
    // atom keeps its own slot and an unreadable one is marked as absent.
    const slotOf: number[] = [];
    let cursor = countsAt + 1;
    for (let i = 0; i < counts.atoms && cursor < lines.length; i++, cursor++) {
        const row = v2000AtomRow(lines[cursor]);
        if (!row) {
            slotOf.push(-1);
            continue;
        }
        const atom: Atom = { element: row.element.element, x: row.x, y: row.y, z: row.z, index: atoms.length };
        if (row.element.ghost) atom.ghost = true;
        atoms.push(atom);
        slotOf.push(atoms.length - 1);
    }
    const unreadableAtoms = slotOf.filter(s => s < 0).length;
    if (unreadableAtoms > 0) {
        warnings.push(unreadableAtoms + ' atom row(s) of ' + counts.atoms +
            ' could not be read: their positions are left empty');
    }
    let skippedBonds = 0;
    const seen = new Set<string>();
    for (let i = 0; i < counts.bonds && cursor < lines.length; i++, cursor++) {
        const row = v2000BondRow(lines[cursor]);
        if (!row) { skippedBonds++; continue; }
        const a1 = row.atom1 - 1;
        const a2 = row.atom2 - 1;
        if (a1 < 0 || a1 >= slotOf.length || a2 < 0 || a2 >= slotOf.length) {
            skippedBonds++;
            continue;
        }
        const m1 = slotOf[a1];
        const m2 = slotOf[a2];
        if (m1 < 0 || m2 < 0) { skippedBonds++; continue; }
        if (!pushBond(bonds, seen, m1, m2, row.order)) skippedBonds++;
    }
    if (skippedBonds > 0) {
        warnings.push(skippedBonds + ' bond row(s) named a missing atom or repeated an existing bond: dropped');
    }
    return {
        atoms, bonds, title,
        hasExplicitBonds: true,
        warnings: warnings.length > 0 ? warnings : undefined
    };
}

function parseV3000(
    lines: string[], countsAt: number, title: string
): MolecularData {
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    const warnings: string[] = [];
    const idToIndex: { [id: number]: number } = {};
    const seen = new Set<string>();
    let section = '';
    let declaredAtoms = 0;
    let declaredBonds = 0;
    let skippedBonds = 0;
    for (let i = countsAt; i < lines.length; i++) {
        const raw = lines[i];
        const m = /^(?:\s*)M\s+V30\s+(.*)$/i.exec(raw);
        if (!m) continue;
        // V3000 continues a long record with a trailing `-`; the next `M  V30`
        // line is its remainder.
        let text = m[1];
        while (text.endsWith('-') && i + 1 < lines.length) {
            const next = /^(?:\s*)M\s+V30\s+(.*)$/i.exec(lines[i + 1]);
            if (!next) break;
            text = text.slice(0, -1) + next[1];
            i++;
        }
        const parts = text.trim().split(/\s+/);
        if (parts.length === 0) continue;
        const head = parts[0].toUpperCase();
        if (head === 'COUNTS') {
            declaredAtoms = parseInt(parts[1], 10) || 0;
            declaredBonds = parseInt(parts[2], 10) || 0;
            continue;
        }
        if (head === 'BEGIN') {
            section = (parts[1] || '').toUpperCase();
            continue;
        }
        if (head === 'END') {
            section = '';
            continue;
        }
        if (section === 'ATOM') {
            const id = parseInt(parts[0], 10);
            const element = molElement(parts[1]);
            const x = parseFloat(parts[2]);
            const y = parseFloat(parts[3]);
            const z = parseFloat(parts[4]);
            if (!isFinite(id) || !element || !isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
            const atom: Atom = { element: element.element, x, y, z, index: atoms.length };
            if (element.ghost) atom.ghost = true;
            idToIndex[id] = atoms.length;
            atoms.push(atom);
            continue;
        }
        if (section === 'BOND') {
            const type = parseInt(parts[1], 10);
            const id1 = parseInt(parts[2], 10);
            const id2 = parseInt(parts[3], 10);
            const i1 = idToIndex[id1];
            const i2 = idToIndex[id2];
            if (!isFinite(type) || i1 === undefined || i2 === undefined) { skippedBonds++; continue; }
            if (!pushBond(bonds, seen, i1, i2, v2000BondOrder(type))) skippedBonds++;
        }
    }
    if (declaredAtoms > 0 && atoms.length < declaredAtoms) {
        warnings.push('the file declares ' + declaredAtoms + ' atoms but only ' +
            atoms.length + ' atom record(s) could be read');
    }
    if (skippedBonds > 0) {
        warnings.push(skippedBonds + ' bond record(s) named a missing atom or repeated an existing bond: dropped');
    }
    return {
        atoms, bonds, title,
        hasExplicitBonds: true,
        warnings: warnings.length > 0 ? warnings : undefined
    };
}

export function parseMol(content: string): MolecularData {
    const lines = String(content || '').replace(/^\uFEFF/, '').split(/\r?\n/);
    const countsAt = findCountsLine(lines);
    const title = (lines[0] || '').trim() || 'MDL Mol';
    if (countsAt < 0) {
        return {
            atoms: [], bonds: [], title,
            hasExplicitBonds: false,
            warnings: ['not an MDL molfile: no V2000 or V3000 counts line was found']
        };
    }
    if (/\bV3000\b/i.test(lines[countsAt])) {
        return parseV3000(lines, countsAt, title);
    }
    return parseV2000(lines, countsAt, title);
}
