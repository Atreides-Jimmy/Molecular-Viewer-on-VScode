import { Atom, Bond, MolecularData } from '../types';
import { elementFromZ, normalizeElementSymbol } from '../tools/periodicTable';

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

/**
 * One row of a MOPAC *input* atom list (`MOPAC` geometry section of a
 * `.mop`/`.dat` input).
 *
 * The classic spelling puts an optimisation flag after each coordinate —
 * `C -4.1918 1 2.5617 0 1.3518 1`, where a `0` means "optimise this
 * coordinate" and a `1` means "hold it fixed". Those flags are NOT a unit
 * marker: reading a `0` as "Bohr" scaled that coordinate by 0.529 (a molecule
 * whose atoms were all frozen came back at 0.529 of its size) and a mixed flag
 * set mixed two units inside one structure, which no shape survived.
 * A MOPAC geometry block is stated in Å, so no unit is derived from it.
 *
 * A bare `C x y z` row (no flags) is legal as well and used to be dropped,
 * because the reader looked for the coordinates in the three odd columns of
 * the flagged form. Both shapes are read positionally now.
 */
function mopacInputRow(line: string): MopacOutputRow | undefined {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) return undefined;
    const upper = parts[0].toUpperCase();
    let element: { element: string; ghost: boolean } | undefined;
    if (upper === 'X' || upper === 'XX') {
        element = { element: 'X', ghost: true };
    } else if (!/^\d+$/.test(parts[0])) {
        // An input row names the element; a leading number would be an output
        // table's serial, which this reader does not take.
        const sym = normalizeElementSymbol(parts[0]);
        if (sym) element = { element: sym, ghost: false };
    }
    if (!element) return undefined;

    const num = (t: string | undefined): number => (t === undefined ? NaN : parseFloat(t));
    // `x f y f z f` — the flags sit between the coordinates, one per axis.
    const flagged = parts.length >= 7 &&
        isMopacFlag(parts[2]) && isMopacFlag(parts[4]) && isMopacFlag(parts[6]);
    // `x y z f f f` — the same three flags written as a group after the
    // coordinates, which some writers produce. A row that matches neither is
    // an internal-coordinate (Z-matrix) row and is left to the caller.
    const grouped = !flagged && parts.length >= 7 &&
        isMopacNumber(parts[1]) && isMopacNumber(parts[2]) && isMopacNumber(parts[3]) &&
        isMopacFlag(parts[4]) && isMopacFlag(parts[5]) && isMopacFlag(parts[6]);
    const x = num(parts[1]);
    const y = flagged ? num(parts[3]) : num(parts[2]);
    const z = flagged ? num(parts[5]) : num(parts[3]);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return undefined;
    if (parts.length >= 7 && !flagged && !grouped) return undefined;
    if (parts.length === 5 || parts.length === 6) return undefined;
    return { element: element.element, ghost: element.ghost, x, y, z };
}

/** `0` / `1` — MOPAC's per-coordinate optimisation flag. */
function isMopacFlag(token: string): boolean {
    return token === '0' || token === '1';
}

/** A token that parses as a coordinate (used only to tell the two row shapes apart). */
function isMopacNumber(token: string): boolean {
    return /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?$/.test(token);
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
        const row = mopacInputRow(cl);
        if (!row) continue;

        // MOPAC's dummy atom is `XX` (and `X` on newer versions); both resolve
        // to the element `X` here, which is the viewer's dummy spelling.
        const atom: Atom = { element: row.element, x: row.x, y: row.y, z: row.z, index: atomIndex };
        if (row.ghost) atom.ghost = true;
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
