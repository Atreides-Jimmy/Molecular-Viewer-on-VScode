import { Atom, AsymAtom, Bond, CrystalData, MolecularData } from '../types';

/**
 * CIF reader.
 *
 * A CIF can state its symmetry under three different dictionary generations,
 * and the same structure is written under all three in the wild:
 *
 *   _symmetry_equiv_pos_as_xyz         CIF 1.1, still what most files write
 *   _space_group_symop_operation_xyz   CIF 1.2 / the core dictionary
 *   _space_group_symop.operation_xyz   CIF 2, category and object split by a dot
 *
 * The same three-generational split applies to the space-group name, and CIF 2
 * writes `_atom_site.fract_x` for what CIF 1.1 calls `_atom_site_fract_x`. A
 * reader that knows one spelling does not fail on the others — it silently
 * reads the file as P1 and shows one quarter of the cell. So every lookup here
 * goes through a normalised tag (see `normTag`) and a priority list of
 * spellings, and a file whose symmetry cannot be read says so through
 * `MolecularData.warnings` rather than passing for a P1 structure.
 *
 * Two deliberate policy decisions:
 *
 *  - A missing value (`?`, `.`, an absent key) is *missing*, never a value.
 *    `?` is not a space group called "?", and a cell of 1 A is not what a file
 *    with no `_cell_length_a` means. `parseNumber` returns null, and the cell
 *    checks turn that into an empty structure plus a warning.
 *  - `warnings` describes what the parser had to decide. It is not a log: a
 *    perfectly ordinary P1 file with its single `x, y, z` operation gets none.
 */

interface CifLoop {
    keys: string[];
    rows: string[][];
}

interface CifData {
    items: { [key: string]: string };
    loops: CifLoop[];
}

// ===== Values =====

/**
 * CIF's two ways of writing "there is no value here": `?` (unknown) and `.`
 * (not applicable). Both mean the field is absent, and nothing downstream ever
 * sees them — an empty string is missing too, which is what a bare `_key` line
 * tokenises to.
 */
function isMissingValue(token: string | undefined | null): boolean {
    if (token === undefined || token === null) return true;
    const t = token.trim();
    return t === '' || t === '?' || t === '.';
}

/**
 * One CIF numeric field, or null when the field says nothing.
 *
 * The parenthesised part of `10.7307(7)` is the standard uncertainty and never
 * part of the value, so the number ends at the first `(` — cutting there also
 * survives an esd whose closing bracket was lost to a clipped column, where
 * deleting the brackets instead would have glued the digits on (`107.307`).
 *
 * A non-string argument is a missing field, not a programming error: the caller
 * is usually `parseNumber(items[key])` and the key may simply not be there.
 */
function parseNumber(token: unknown): number | null {
    if (typeof token !== 'string' || isMissingValue(token)) return null;
    let s = token.trim();
    const paren = s.indexOf('(');
    if (paren >= 0) s = s.substring(0, paren);
    s = s.trim();
    if (s === '') return null;
    const v = parseFloat(s);
    return Number.isFinite(v) ? v : null;
}

/**
 * An element symbol out of whatever `_atom_site_type_symbol` holds.
 *
 * The column is meant to hold a bare symbol, but files carry the oxidation
 * state with it (`Fe3+`, `O2-`, `Cu(II)`, `N+`). Left alone, that string becomes
 * the element name and every table lookup — radius, colour, valence, mass —
 * misses, so the atom is drawn as an unknown grey ball and never bonds. The
 * symbol is the leading run of letters, which is also what makes a label used
 * as a fallback (`C12`, `O1W`) come out as `C` and `O`.
 */
function normaliseElement(raw: string): string {
    const s = String(raw === undefined || raw === null ? '' : raw).trim();
    if (s === '' || s === '?' || s === '.') return '';
    // 'Cu(II)' -> 'Cu'; the state may also be written '(+2)' or '(2+)'.
    const bare = s.replace(/\([^)]*\)/g, ' ');
    const m = bare.match(/[A-Za-z]{1,2}/);
    if (!m) return '';
    return m[0].charAt(0).toUpperCase() + m[0].slice(1).toLowerCase();
}

// ===== Tags =====

/**
 * A CIF data name in one spelling.
 *
 * CIF 2 splits a data name into category and object with a dot and CIF 1.1 uses
 * an underscore for the same boundary, so `_space_group_symop.operation_xyz`
 * and `_space_group_symop_operation_xyz` are one tag written two ways. Lookups
 * compare this form, which also means they stop guessing with substrings:
 * `_geom_bond_atom_site_label_1` must not answer a question about
 * `_atom_site_label`, and `_atom_site_aniso_label` must not answer one about
 * `label`.
 */
function normTag(tag: string): string {
    const s = tag.trim().toLowerCase();
    const dot = s.indexOf('.');
    return dot < 0 ? s : s.substring(0, dot) + '_' + s.substring(dot + 1);
}

/** Column of the first of `tags` that this loop carries, or -1. */
function columnOf(loop: CifLoop, tags: string[]): number {
    for (const tag of tags) {
        const want = normTag(tag);
        for (let i = 0; i < loop.keys.length; i++) {
            if (normTag(loop.keys[i]) === want) return i;
        }
    }
    return -1;
}

/**
 * The first loop carrying one of `tags`, searched tag by tag so the caller
 * states a priority instead of the file's order deciding it. The column comes
 * back with the loop because the tag may sit anywhere in the key list.
 */
function findLoopAny(data: CifData, tags: string[]): { loop: CifLoop; col: number } | null {
    for (const tag of tags) {
        for (const loop of data.loops) {
            const col = columnOf(loop, [tag]);
            if (col >= 0) return { loop, col };
        }
    }
    return null;
}

/**
 * The first of `tags` that has a value in the header items (not in a loop).
 * Missing values are skipped, so `_symmetry_space_group_name_H-M ?` falls
 * through to the next spelling instead of becoming the string '?'.
 */
function findItemAny(data: CifData, tags: string[]): string | undefined {
    for (const tag of tags) {
        const want = normTag(tag);
        for (const key of Object.keys(data.items)) {
            if (normTag(key) !== want) continue;
            const v = data.items[key];
            if (!isMissingValue(v)) return v;
        }
    }
    return undefined;
}

/** Single-line form of a header value: semicolon blocks arrive multi-line. */
function oneLine(value: string): string {
    return value.split(/\r?\n/).map(s => s.trim()).filter(s => s !== '').join(' ');
}

function tokenize(content: string): CifData {
    const data: CifData = { items: {}, loops: [] };
    const lines = content.split(/\r?\n/);
    let i = 0;
    let inSemicolon = false;
    let semicolonBuffer: string[] = [];
    let semicolonKey: string | null = null;
    // A bare `_key` line: its value is the semicolon block that may follow.
    let pendingKey: string | null = null;

    function flushSemicolon() {
        if (semicolonKey !== null) {
            data.items[semicolonKey] = semicolonBuffer.join('\n').trim();
        }
        inSemicolon = false;
        semicolonBuffer = [];
        semicolonKey = null;
    }

    function readValue(rest: string): { value: string; consumed: number } {
        const trimmed = rest.trimStart();
        if (trimmed.startsWith("'")) {
            const end = trimmed.indexOf("'", 1);
            if (end >= 0) return { value: trimmed.substring(1, end), consumed: (rest.length - trimmed.length) + end + 1 };
        }
        if (trimmed.startsWith('"')) {
            const end = trimmed.indexOf('"', 1);
            if (end >= 0) return { value: trimmed.substring(1, end), consumed: (rest.length - trimmed.length) + end + 1 };
        }
        const parts = trimmed.split(/\s+/);
        return { value: parts[0], consumed: rest.length - trimmed.length + parts[0].length };
    }

    while (i < lines.length) {
        const rawLine = lines[i];
        const line = rawLine.trim();

        if (inSemicolon) {
            if (line.startsWith(';')) {
                flushSemicolon();
                i++;
                continue;
            }
            semicolonBuffer.push(rawLine);
            i++;
            continue;
        }

        if (line === '' || line.startsWith('#')) {
            i++;
            continue;
        }

        if (line.startsWith(';')) {
            inSemicolon = true;
            semicolonBuffer = [];
            // A text block belongs to the data name written above it, which is
            // how a multi-line value is spelled in CIF:
            //     _chemical_name_systematic
            //     ;
            //     ...
            //     ;
            // Without this the block would be parsed and then dropped, and the
            // item above it would read as empty.
            semicolonKey = pendingKey;
            pendingKey = null;
            const after = line.substring(1).trim();
            if (after !== '') semicolonBuffer.push(after);
            i++;
            continue;
        }

        if (line.startsWith('data_')) {
            pendingKey = null;
            i++;
            continue;
        }

        if (line.startsWith('save_')) {
            pendingKey = null;
            i++;
            continue;
        }

        if (line.startsWith('loop_')) {
            pendingKey = null;
            const loop: CifLoop = { keys: [], rows: [] };
            i++;
            let startedRows = false;
            while (i < lines.length) {
                const l = lines[i].trim();
                if (l === '' || l.startsWith('#')) { i++; continue; }
                if (l.startsWith('_')) {
                    if (startedRows) break;
                    const parts = l.split(/\s+/);
                    if (parts.length > 1) {
                        // Key with inline value = scalar item, not part of loop
                        break;
                    }
                    loop.keys.push(parts[0]);
                    i++;
                } else if (l.startsWith('loop_') || l.startsWith('data_') || l.startsWith('save_')) {
                    break;
                } else {
                    startedRows = true;
                    const row: string[] = [];
                    let lineRest = lines[i];
                    for (let k = 0; k < loop.keys.length; k++) {
                        if (lineRest.trim() === '') break;
                        if (lineRest.trim().startsWith(';')) {
                            const buf: string[] = [];
                            const firstPart = lineRest.trim().substring(1);
                            if (firstPart !== '') buf.push(firstPart);
                            i++;
                            while (i < lines.length && !lines[i].trim().startsWith(';')) {
                                buf.push(lines[i]);
                                i++;
                            }
                            row.push(buf.join('\n').trim());
                            lineRest = '';
                            if (i < lines.length) {
                                lineRest = lines[i].replace(/^\s*;\s*/, '');
                            }
                        } else {
                            const { value, consumed } = readValue(lineRest);
                            row.push(value);
                            lineRest = lineRest.substring(consumed);
                        }
                    }
                    if (row.length > 0) loop.rows.push(row);
                    i++;
                }
            }
            if (loop.keys.length > 0) {
                if (loop.rows.length === 0 || loop.rows[0].length === loop.keys.length) {
                    data.loops.push(loop);
                }
            }
            continue;
        }

        if (line.startsWith('_')) {
            const match = line.match(/^(_\S+)\s+(.*)$/);
            if (match) {
                let value = match[2].trim();
                if (value.startsWith(';')) {
                    inSemicolon = true;
                    semicolonKey = match[1];
                    semicolonBuffer = [];
                    const after = value.substring(1);
                    if (after.trim() !== '') semicolonBuffer.push(after);
                    i++;
                    continue;
                }
                if (value.startsWith("'") && value.endsWith("'") && value.length >= 2) {
                    value = value.substring(1, value.length - 1);
                } else if (value.startsWith('"') && value.endsWith('"') && value.length >= 2) {
                    value = value.substring(1, value.length - 1);
                } else {
                    value = value.split(/\s+/)[0];
                }
                data.items[match[1]] = value;
                pendingKey = null;
            } else {
                const key = line.split(/\s+/)[0];
                data.items[key] = '';
                pendingKey = key;
            }
            i++;
            continue;
        }

        pendingKey = null;
        i++;
    }

    flushSemicolon();
    return data;
}

// ===== Cell =====

function cellToLatticeVectors(a: number, b: number, c: number, alpha: number, beta: number, gamma: number): number[][] {
    const rad = Math.PI / 180;
    const ar = alpha * rad;
    const br = beta * rad;
    const gr = gamma * rad;

    const cosA = Math.cos(ar);
    const cosB = Math.cos(br);
    const cosG = Math.cos(gr);
    const sinG = Math.sin(gr);

    const aVec = [a, 0, 0];
    const bVec = [b * cosG, b * sinG, 0];
    const cX = c * cosB;
    const cY = c * (cosA - cosB * cosG) / sinG;
    const cZ2 = 1 - cosA * cosA - cosB * cosB - cosG * cosG + 2 * cosA * cosB * cosG;
    const cZ = c * Math.sqrt(Math.max(0, cZ2)) / sinG;
    const cVec = [cX, cY, cZ];

    return [aVec, bVec, cVec];
}

function fracToCart(frac: number[], lattice: number[][]): [number, number, number] {
    const x = frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0];
    const y = frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1];
    const z = frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2];
    return [x, y, z];
}

// ===== Symmetry =====

interface SymOp {
    R: number[][];
    t: number[];
}

function parseSymOp(s: string): SymOp {
    const parts = s.split(',').map(p => p.trim());
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const t = [0, 0, 0];

    for (let idx = 0; idx < 3; idx++) {
        if (idx >= parts.length) continue;
        const expr = parts[idx].replace(/\s+/g, '').toLowerCase();
        let sign = 1;
        let k = 0;
        while (k < expr.length) {
            const ch = expr[k];
            if (ch === '+') { sign = 1; k++; continue; }
            if (ch === '-') { sign = -1; k++; continue; }
            if (ch === 'x' || ch === 'y' || ch === 'z') {
                const col = ch === 'x' ? 0 : ch === 'y' ? 1 : 2;
                R[idx][col] += sign;
                sign = 1;
                k++;
                continue;
            }
            const numMatch = expr.substring(k).match(/^(\d+\/\d+|\d+\.?\d*)/);
            if (numMatch) {
                let v: number;
                if (numMatch[1].includes('/')) {
                    const [p, q] = numMatch[1].split('/').map(Number);
                    v = p / q;
                } else {
                    v = parseFloat(numMatch[1]);
                }
                t[idx] += sign * v;
                sign = 1;
                k += numMatch[1].length;
                continue;
            }
            k++;
        }
    }

    return { R, t };
}

function applySymOp(op: SymOp, frac: number[]): number[] {
    const result = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        result[i] = op.R[i][0] * frac[0] + op.R[i][1] * frac[1] + op.R[i][2] * frac[2] + op.t[i];
    }
    return result;
}

function wrapFrac(v: number): number {
    let r = v - Math.floor(v);
    if (r >= 1 - 1e-6) r -= 1;
    if (r < 0) r += 1;
    return r;
}

/**
 * Whether a string is an operation this parser can apply: `x`, `y` and `z` each
 * appear, and there are at least three comma-separated parts. Filenames, unit
 * cells and quoted notes sometimes end up in the same loop when a file is
 * malformed, and counting those as operations would place atoms at nonsense.
 */
function looksLikeSymOp(s: string): boolean {
    if (s.indexOf(',') < 0) return false;
    if (!/[xyzXYZ]/.test(s)) return false;
    if (!/x/i.test(s) || !/y/i.test(s) || !/z/i.test(s)) return false;
    return /^[xyzXYZ0-9+\-/.,\s]+$/.test(s);
}

/**
 * The operations the file states, and which spelling carried them.
 *
 * A file states its operations once, so the first spelling that yields a
 * readable operation wins; reading the other spellings too would double every
 * site's images in the files that carry both. `sawTag` distinguishes "the file
 * says nothing about symmetry" from "the file says something this parser could
 * not read" — both fall back to P1, but only one of them is the file's fault.
 */
const SYMOP_TAGS = [
    '_symmetry_equiv_pos_as_xyz',       // CIF 1.1 — what most files still write
    '_space_group_symop_operation_xyz', // CIF 1.2 / core dictionary
    '_space_group_symop.operation_xyz'  // CIF 2 spelling
];

const SPACE_GROUP_TAGS = [
    '_space_group_name_H-M_alt',        // core dictionary (and CIF 2, dotted)
    '_space_group.name_H-M_alt',
    '_symmetry_space_group_name_H-M',   // CIF 1.1
    '_space_group_name_H-M_ref'
];

function collectSymOps(data: CifData): { ops: string[]; tag: string | null; sawTag: string | null } {
    let sawTag: string | null = null;
    for (const tag of SYMOP_TAGS) {
        const raw: string[] = [];
        const loopHit = findLoopAny(data, [tag]);
        if (loopHit) {
            if (sawTag === null) sawTag = tag;
            for (const row of loopHit.loop.rows) {
                const v = row[loopHit.col];
                if (typeof v === 'string') raw.push(v);
            }
        }
        // A structure with one operation is sometimes written as a plain item
        // rather than a one-row loop.
        const scalar = findItemAny(data, [tag]);
        if (scalar !== undefined) {
            if (sawTag === null) sawTag = tag;
            raw.push(scalar);
        }

        const ops: string[] = [];
        const seen: { [key: string]: true } = {};
        for (const chunk of raw) {
            // A semicolon block holds one operation per line.
            for (const line of chunk.split(/\r?\n/)) {
                const s = line.trim().replace(/^;+/, '').replace(/;+$/, '').trim();
                if (!looksLikeSymOp(s)) continue;
                const key = s.toLowerCase().replace(/\s+/g, '');
                if (seen[key]) continue;
                seen[key] = true;
                ops.push(s);
            }
        }
        if (ops.length > 0) return { ops, tag, sawTag };
    }
    return { ops: [], tag: null, sawTag };
}

// ===== Coincident sites =====

/**
 * How close two images must be before they are the same atom, in Angstrom.
 *
 * This has to be far below a bond length: a structure can legitimately list two
 * sites 0.3 A apart (a split position, a partially occupied pair) and folding
 * them into one loses an atom that the file says is there. A tenth of an
 * angstrom only catches what symmetry expansion produces on purpose — the
 * images of a site that sits on a special position — plus genuine duplicates.
 */
const DEDUPE_TOL = 0.1;
const DEDUPE_TOL2 = DEDUPE_TOL * DEDUPE_TOL;

/**
 * Two images closer than this are the same point, not two points that happen to
 * be near each other: the difference is round-off in the fractional -> cartesian
 * round trip (a few 1e-9 A at most). Only merges *above* it say anything about
 * the structure rather than about the arithmetic, so only those are reported.
 */
const EXACT_MERGE = 1e-4;

interface NearIndex {
    add(x: number, y: number, z: number): void;
    /** Nearest atom within DEDUPE_TOL of (x,y,z), with its distance. */
    find(x: number, y: number, z: number): { index: number; dist: number } | null;
}

/**
 * A grid over the atoms placed so far, one cell per tolerance.
 *
 * Symmetry expansion compares every image against every atom already placed, so
 * a framework (thousands of sites times hundreds of operations) turns the naive
 * scan into billions of distance checks. A grid makes each lookup local; the
 * neighbourhood is 3x3x3 cells because two points closer than one cell size can
 * only ever land in adjacent cells. The result is the same as the exhaustive
 * scan — the distance test still decides — only it no longer depends on how
 * many atoms the cell already holds.
 */
function createNearIndex(): NearIndex {
    const buckets = new Map<number, number[]>();
    const xs: number[] = [];
    const ys: number[] = [];
    const zs: number[] = [];
    const cellOf = (v: number) => Math.floor(v / DEDUPE_TOL);
    // A number hash of the three cell indices. Collisions only cost a distance
    // check, never a wrong answer.
    const bucketKey = (ix: number, iy: number, iz: number) =>
        (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791);

    return {
        add(x: number, y: number, z: number): void {
            const k = bucketKey(cellOf(x), cellOf(y), cellOf(z));
            const bucket = buckets.get(k);
            if (bucket) bucket.push(xs.length);
            else buckets.set(k, [xs.length]);
            xs.push(x);
            ys.push(y);
            zs.push(z);
        },
        find(x: number, y: number, z: number): { index: number; dist: number } | null {
            const ix = cellOf(x), iy = cellOf(y), iz = cellOf(z);
            let best = -1;
            let bestD2 = DEDUPE_TOL2;
            for (let dx = -1; dx <= 1; dx++) {
                for (let dy = -1; dy <= 1; dy++) {
                    for (let dz = -1; dz <= 1; dz++) {
                        const bucket = buckets.get(bucketKey(ix + dx, iy + dy, iz + dz));
                        if (!bucket) continue;
                        for (let n = 0; n < bucket.length; n++) {
                            const i = bucket[n];
                            const ex = xs[i] - x, ey = ys[i] - y, ez = zs[i] - z;
                            const d2 = ex * ex + ey * ey + ez * ez;
                            if (d2 < bestD2) {
                                bestD2 = d2;
                                best = i;
                            }
                        }
                    }
                }
            }
            return best < 0 ? null : { index: best, dist: Math.sqrt(bestD2) };
        }
    };
}

// ===== Entry point =====

/** Every site merged away by an inexact match is worth one line, but no more. */
const MAX_MERGE_WARNINGS = 6;

/**
 * The atom-site loop, in fractional coordinates when the file has them.
 *
 * A loop only counts when it carries all three axes: `_atom_site_aniso_label`
 * and friends live in loops of their own, and picking one of those would read
 * the file's anisotropy table as its atoms.
 */
function atomSiteLoop(data: CifData): { loop: CifLoop; useFract: boolean } | null {
    for (const useFract of [true, false]) {
        const axis = useFract ? '_atom_site_fract_' : '_atom_site_Cartn_';
        const hit = findLoopAny(data, [axis + 'x']);
        if (!hit) continue;
        if (columnOf(hit.loop, [axis + 'y']) < 0) continue;
        if (columnOf(hit.loop, [axis + 'z']) < 0) continue;
        return { loop: hit.loop, useFract };
    }
    return null;
}

/**
 * A CIF with a cell but no sites. The cell that was read still comes back — the
 * wireframe and the cell panel are what the file says, and are worth showing on
 * their own — but there is nothing to draw inside it and the warning says why.
 *
 * A file missing its cell lengths does *not* come through here: there is no
 * cell to show, and inventing a 1 A cube would put a picture on screen that the
 * file never described (see `parseCif`).
 */
function emptyCellCif(
    title: string,
    warnings: string[],
    cell: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number },
    spaceGroup: string | undefined,
    symmetryOps: string[]
): MolecularData {
    const lattice = cellToLatticeVectors(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
    return {
        atoms: [],
        bonds: [],
        title,
        hasExplicitBonds: false,
        warnings,
        crystal: {
            a: cell.a, b: cell.b, c: cell.c,
            alpha: cell.alpha, beta: cell.beta, gamma: cell.gamma,
            latticeVectors: lattice,
            spaceGroup,
            symmetryOps,
            baseAtoms: [],
            baseBonds: [],
            asymAtoms: []
        }
    };
}

export function parseCif(content: string): MolecularData {
    const data = tokenize(content);
    const warnings: string[] = [];

    // ----- cell -----
    // Lengths are required to build a cell at all; angles default to 90 with a
    // warning, because an orthogonal cell is what a file that omits them means
    // and refusing the structure would be worse than saying so.
    const aVal = parseNumber(findItemAny(data, ['_cell_length_a']));
    const bVal = parseNumber(findItemAny(data, ['_cell_length_b']));
    const cVal = parseNumber(findItemAny(data, ['_cell_length_c']));
    const alphaVal = parseNumber(findItemAny(data, ['_cell_angle_alpha']));
    const betaVal = parseNumber(findItemAny(data, ['_cell_angle_beta']));
    const gammaVal = parseNumber(findItemAny(data, ['_cell_angle_gamma']));

    const missingLengths: string[] = [];
    if (!(aVal !== null && aVal > 0)) missingLengths.push('a');
    if (!(bVal !== null && bVal > 0)) missingLengths.push('b');
    if (!(cVal !== null && cVal > 0)) missingLengths.push('c');

    const a = aVal !== null && aVal > 0 ? aVal : 1;
    const b = bVal !== null && bVal > 0 ? bVal : 1;
    const c = cVal !== null && cVal > 0 ? cVal : 1;
    const alpha = alphaVal !== null && alphaVal > 0 ? alphaVal : 90;
    const beta = betaVal !== null && betaVal > 0 ? betaVal : 90;
    const gamma = gammaVal !== null && gammaVal > 0 ? gammaVal : 90;
    const cell = { a, b, c, alpha, beta, gamma };

    const sgRaw = findItemAny(data, SPACE_GROUP_TAGS);
    const spaceGroup = sgRaw === undefined ? undefined : oneLine(sgRaw);
    const spaceGroupOut = spaceGroup === undefined || spaceGroup === '' ? undefined : spaceGroup;

    const scan = collectSymOps(data);
    const symmetryOps = scan.ops.length > 0 ? scan.ops.slice() : ['x, y, z'];

    // The cell comes first: without one there is nothing to expand, and a
    // warning about symmetry would only bury the reason the view is empty.
    // Nothing is substituted for the missing lengths — a structure drawn in a
    // 1 A cube would be a picture of something the file never said.
    if (missingLengths.length > 0) {
        warnings.push('unit cell incomplete (' + missingLengths.join(', ') + ' missing): no structure to show');
        return {
            atoms: [],
            bonds: [],
            title: 'No unit cell in CIF',
            hasExplicitBonds: false,
            warnings
        };
    }
    if (alphaVal === null || betaVal === null || gammaVal === null) {
        warnings.push('cell angles missing: the missing ones read as 90 degrees');
    }

    if (scan.ops.length === 0) {
        if (scan.sawTag) {
            warnings.push('no readable symmetry operation under ' + scan.sawTag + ': read as P1');
        } else {
            warnings.push('no symmetry operations in the file: read as P1' +
                (spaceGroupOut ? ' (declared space group ' + spaceGroupOut + ')' : ''));
        }
    }

    const lattice = cellToLatticeVectors(a, b, c, alpha, beta, gamma);

    // ----- sites -----
    const site = atomSiteLoop(data);
    if (!site) {
        warnings.push('no atom sites in the file');
        return emptyCellCif('No atom sites found in CIF', warnings, cell, spaceGroupOut, symmetryOps);
    }

    const atomLoop = site.loop;
    const useFract = site.useFract;
    const axis = useFract ? '_atom_site_fract_' : '_atom_site_Cartn_';
    const labelCol = columnOf(atomLoop, ['_atom_site_label']);
    const typeCol = columnOf(atomLoop, ['_atom_site_type_symbol']);
    const xCol = columnOf(atomLoop, [axis + 'x']);
    const yCol = columnOf(atomLoop, [axis + 'y']);
    const zCol = columnOf(atomLoop, [axis + 'z']);
    const occCol = columnOf(atomLoop, ['_atom_site_occupancy']);

    interface RawAtom {
        element: string;
        label: string;
        frac: number[];
        occupancy: number;
    }
    const rawAtoms: RawAtom[] = [];

    const cols = [labelCol, typeCol, xCol, yCol, zCol].filter(col => col >= 0);
    const minCols = cols.length > 0 ? Math.max(...cols) + 1 : 0;
    let skipped = 0;

    for (const row of atomLoop.rows) {
        if (row.length < minCols) { skipped++; continue; }
        const label = labelCol >= 0 ? (row[labelCol] || '') : '';
        const typeSym = typeCol >= 0 ? row[typeCol] : '';
        // The label is the fallback the dictionary itself suggests (`C12` is a
        // carbon), and it is also what a file with no type column has.
        const element = normaliseElement(isMissingValue(typeSym) ? label : typeSym) ||
            normaliseElement(label) || 'C';
        const xv = parseNumber(row[xCol]);
        const yv = parseNumber(row[yCol]);
        const zv = parseNumber(row[zCol]);
        if (xv === null || yv === null || zv === null) { skipped++; continue; }
        const occupancy = occCol >= 0 ? (parseNumber(row[occCol]) ?? 1) : 1;

        const frac = useFract ? [xv, yv, zv] : cartToFrac([xv, yv, zv], lattice);
        rawAtoms.push({ element, label, frac, occupancy });
    }
    if (skipped > 0) {
        warnings.push(skipped + ' atom site' + (skipped === 1 ? '' : 's') +
            ' skipped: the row is short or its coordinates are not numbers');
    }

    const parsedSymOps = symmetryOps.map(parseSymOp);

    // ----- expand -----
    const baseAtoms: Atom[] = [];
    const asymAtoms: AsymAtom[] = [];
    const nearMisses: { label: string; dist: number }[] = [];
    const index = createNearIndex();

    for (let siteIdx = 0; siteIdx < rawAtoms.length; siteIdx++) {
        const raw = rawAtoms[siteIdx];
        const hitOps: number[] = [];
        for (let opIdx = 0; opIdx < parsedSymOps.length; opIdx++) {
            const transformed = applySymOp(parsedSymOps[opIdx], raw.frac);
            const wrapped = transformed.map(wrapFrac);
            const [cx, cy, cz] = fracToCart(wrapped, lattice);

            const near = index.find(cx, cy, cz);
            if (near) {
                // An image that lands on an atom already there is what symmetry
                // does to a site on a special position — normal, and not worth a
                // line. An image that lands *near* one means two sites the file
                // lists separately are being folded together, which is the case
                // the tolerance exists for and the one the user must hear about.
                if (near.dist > EXACT_MERGE) {
                    nearMisses.push({ label: raw.label || raw.element + (siteIdx + 1), dist: near.dist });
                }
                continue;
            }

            index.add(cx, cy, cz);
            baseAtoms.push({
                element: raw.element,
                x: cx, y: cy, z: cz,
                index: baseAtoms.length,
                occupancy: raw.occupancy,
                baseIdx: baseAtoms.length
            });
            hitOps.push(opIdx);
        }

        // The asymmetric unit keeps the site as listed, in the same place as
        // its identity image, so a save-back to CIF can write it verbatim.
        const siteFrac = [wrapFrac(raw.frac[0]), wrapFrac(raw.frac[1]), wrapFrac(raw.frac[2])];
        const [sx, sy, sz] = fracToCart(siteFrac, lattice);
        asymAtoms.push({
            element: raw.element,
            label: raw.label,
            x: sx, y: sy, z: sz,
            frac: [siteFrac[0], siteFrac[1], siteFrac[2]],
            index: asymAtoms.length,
            occupancy: raw.occupancy,
            symmetryOps: hitOps
        });
    }

    for (let i = 0; i < nearMisses.length && i < MAX_MERGE_WARNINGS; i++) {
        const m = nearMisses[i];
        warnings.push('site ' + m.label + ': an image ' + m.dist.toFixed(3) +
            ' A from an existing atom was merged into it');
    }
    if (nearMisses.length > MAX_MERGE_WARNINGS) {
        warnings.push('and ' + (nearMisses.length - MAX_MERGE_WARNINGS) +
            ' more images merged into an atom less than ' + DEDUPE_TOL + ' A away');
    }

    const atoms: Atom[] = baseAtoms.map(at => ({
        element: at.element,
        x: at.x, y: at.y, z: at.z,
        index: at.index,
        occupancy: at.occupancy,
        baseIdx: at.index,
        cellI: 0, cellJ: 0, cellK: 0
    }));

    const crystal: CrystalData = {
        a, b, c, alpha, beta, gamma,
        latticeVectors: lattice,
        spaceGroup: spaceGroupOut,
        symmetryOps,
        baseAtoms: baseAtoms.map(at => ({ ...at, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: [] as Bond[],
        asymAtoms
    };

    // The title chain skips missing values, so a file whose `_chemical_name_common`
    // is `?` is titled by its formula instead of by a question mark.
    const title = findItemAny(data, ['_chemical_name_common', '_chemical_formula_sum']) || 'CIF Structure';

    return {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        warnings,
        crystal
    };
}

function cartToFrac(cart: number[], lattice: number[][]): number[] {
    const [aVec, bVec, cVec] = lattice;
    const det = aVec[0] * (bVec[1] * cVec[2] - bVec[2] * cVec[1]) -
                aVec[1] * (bVec[0] * cVec[2] - bVec[2] * cVec[0]) +
                aVec[2] * (bVec[0] * cVec[1] - bVec[1] * cVec[0]);
    if (Math.abs(det) < 1e-10) return [0, 0, 0];
    const invDet = 1 / det;
    const inv = [
        [
            (bVec[1] * cVec[2] - bVec[2] * cVec[1]) * invDet,
            (bVec[2] * cVec[0] - bVec[0] * cVec[2]) * invDet,
            (bVec[0] * cVec[1] - bVec[1] * cVec[0]) * invDet
        ],
        [
            (aVec[2] * cVec[1] - aVec[1] * cVec[2]) * invDet,
            (aVec[0] * cVec[2] - aVec[2] * cVec[0]) * invDet,
            (aVec[1] * cVec[0] - aVec[0] * cVec[1]) * invDet
        ],
        [
            (aVec[1] * bVec[2] - aVec[2] * bVec[1]) * invDet,
            (aVec[2] * bVec[0] - aVec[0] * bVec[2]) * invDet,
            (aVec[0] * bVec[1] - aVec[1] * bVec[0]) * invDet
        ]
    ];
    return [
        inv[0][0] * cart[0] + inv[0][1] * cart[1] + inv[0][2] * cart[2],
        inv[1][0] * cart[0] + inv[1][1] * cart[1] + inv[1][2] * cart[2],
        inv[2][0] * cart[0] + inv[2][1] * cart[1] + inv[2][2] * cart[2]
    ];
}

export function cartesianToFractional(cart: number[], lattice: number[][]): number[] {
    return cartToFrac(cart, lattice);
}

export function fractionalToCartesian(frac: number[], lattice: number[][]): [number, number, number] {
    return fracToCart(frac, lattice);
}
