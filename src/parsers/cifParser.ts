import { Atom, AsymAtom, Bond, CrystalData, CrystalSymmetryAlternative, CrystalSymmetryCheck, MolecularData } from '../types';
import {
    applySymOp,
    cartToFrac,
    cellToLatticeVectors,
    createNearIndex,
    DEDUPE_TOL,
    EXACT_MERGE,
    fracToCart,
    parseSymOp,
    SymOp,
    wrapFrac
} from '../tools/crystalGeometry';
import { lookupSpaceGroupOps, SpaceGroupOps } from '../tools/spaceGroupTable';
import { AnisoTensor, rotateAnisoTensor, symOpCartesianRotation } from '../tools/adp';

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

/**
 * One `data_` block. A CIF may hold several structures (a two-phase powder
 * refinement, two compounds in one supporting-information file), and each block
 * starts a fresh item namespace — merging them would let the second block's
 * `_cell_length_a` overwrite the first and append its sites to the first
 * block's loop, producing a structure that is neither.
 */
interface CifBlock extends CifData {
    name: string;
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

function tokenize(content: string): CifBlock[] {
    const blocks: CifBlock[] = [];
    let data: CifBlock = { name: '', items: {}, loops: [] };
    const lines = content.split(/\r?\n/);
    let i = 0;
    let inSemicolon = false;
    let semicolonBuffer: string[] = [];
    let semicolonKey: string | null = null;
    // A bare `_key` line: its value is the semicolon block that may follow.
    let pendingKey: string | null = null;

    function startBlock(name: string) {
        // Tokens written before the first `data_` name belong to the first
        // block, exactly as they did when the whole file was one namespace.
        if (blocks.length === 0 && data.name === '') {
            data.name = name;
            return;
        }
        blocks.push(data);
        data = { name, items: {}, loops: [] };
    }

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
            startBlock(line.substring(5).trim());
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
    blocks.push(data);
    return blocks;
}

// ===== Symmetry =====

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

/**
 * The warning a file gets when its space-group name, not its own table, supplied
 * the operations. It names the group and setting the table used, because the
 * coordinates in the file were written for *some* setting: a reader that picks
 * a different origin or cell choice would draw a different, equally plausible
 * crystal. The origin caveat is only added when the choice really was a guess.
 */
function builtInSymmetryWarning(hit: SpaceGroupOps): string {
    let s = 'no symmetry operations in the file: applied ' + hit.ops.length +
        ' operation' + (hit.ops.length === 1 ? '' : 's') + ' of ' + hit.symbol +
        ' (No. ' + hit.number + (hit.settingLabel === '' ? '' : ', ' + hit.settingLabel) +
        ') from the built-in space-group table';
    if (hit.assumedOrigin) {
        s += '; the file does not state which origin its coordinates use, and the other origin would expand them differently';
    }
    return s;
}

// ===== Coincident sites =====
// The tolerance, its grid index and the "only round-off is quiet" rule are
// shared with the VESTA reader — see crystalGeometry.ts. Two readers expanding
// the same structure must use one ruler, or the same crystal gets two atom
// counts.

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
    symmetryOps: string[],
    symmetrySource?: 'file' | 'builtin'
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
            symmetrySource,
            baseAtoms: [],
            baseBonds: [],
            asymAtoms: []
        }
    };
}

// ===== Expanding the asymmetric unit =====

/**
 * One atom site as the file lists it, before symmetry places its images.
 * Module level because a file whose name and own operations disagree is
 * expanded twice (see the space-group check below).
 */
interface CifRawAtom {
    element: string;
    label: string;
    frac: number[];
    occupancy: number;
    disorderGroup?: string;
    uAniso?: AnisoTensor;
}

interface CifExpansion {
    baseAtoms: Atom[];
    asymAtoms: AsymAtom[];
    nearMisses: { label: string; dist: number }[];
}

/**
 * Multiply every listed site by every operation and deduplicate the images,
 * keeping the sites themselves as the asymmetric unit.
 *
 * The file's operations and the declared name's operations both run through
 * this function, so the two atom counts a mismatch report states are computed
 * one way and are comparable.
 */
function expandCifSites(rawAtoms: CifRawAtom[], parsedSymOps: SymOp[], lattice: number[][]): CifExpansion {
    // One Cartesian rotation per operation, shared by every site: an ADP is a
    // tensor, so the same operation must rotate all of them the same way.
    const symCartRot = parsedSymOps.map(op => symOpCartesianRotation(op, lattice));

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
                baseIdx: baseAtoms.length,
                // The file's own name for the site this image came from; empty
                // labels stay undefined so a file without labels adds nothing.
                siteLabel: raw.label || undefined,
                // Where this image came from and which operation placed it —
                // the two numbers a tooltip needs to say "site O1, op #7".
                asymIdx: asymAtoms.length,
                symOp: opIdx,
                disorderGroup: raw.disorderGroup,
                uAniso: raw.uAniso ? rotateAnisoTensor(raw.uAniso, symCartRot[opIdx]) : undefined
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
            symmetryOps: hitOps,
            disorderGroup: raw.disorderGroup,
            uAniso: raw.uAniso
        });
    }

    return { baseAtoms, asymAtoms, nearMisses };
}

// ===== Space-group self-consistency check =====
//
// A file may state its operations and, separately, a space-group name. The
// operations are the authority — the coordinates were refined under them — so
// an inconsistent name changes nothing in the picture. But the pair disagreeing
// means the file contradicts itself, and silently drawing from one half of it
// is how a structure ends up on screen in a setting nobody can reproduce. The
// check compares the two operation sets and, when they differ, expands the
// sites under the name too, so the warning can state both atom counts and the
// page can offer a one-time switch.
//
// The comparison is exact in 24ths of a cell (the denominator every space group
// uses) and ignores a difference that is a lattice centring translation: a
// centered group's operations come in centring copies, and a file that lists
// them in a different order, or with one copy shifted by the centring vector,
// still describes the same group.

interface ExactCifSymOp {
    R: number[][];
    t: number[];
}

const CIF_SYM_UNITS = 24;

function cifExactSymOp(op: SymOp): ExactCifSymOp | null {
    const R: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const t = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            const v = Math.round(op.R[i][j]);
            // A rotation coefficient that is not an integer is not a space-group
            // operation; it cannot equal any operation of the table.
            if (Math.abs(op.R[i][j] - v) > 1e-6) return null;
            R[i][j] = v;
        }
        const q = Math.round(op.t[i] * CIF_SYM_UNITS);
        t[i] = ((q % CIF_SYM_UNITS) + CIF_SYM_UNITS) % CIF_SYM_UNITS;
    }
    return { R, t };
}

/** Distinct exact operations, keyed in 24ths so a repeated row collapses. */
function cifExactOps(ops: SymOp[]): ExactCifSymOp[] | null {
    const out: ExactCifSymOp[] = [];
    const seen: { [key: string]: true } = {};
    for (const op of ops) {
        const exact = cifExactSymOp(op);
        if (exact === null) return null;
        const key = exact.R[0].join(',') + ';' + exact.R[1].join(',') + ';' + exact.R[2].join(',') +
            '|' + exact.t.join(',');
        if (seen[key]) continue;
        seen[key] = true;
        out.push(exact);
    }
    return out;
}

function cifRotationsEqual(a: number[][], b: number[][]): boolean {
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            if (a[i][j] !== b[i][j]) return false;
        }
    }
    return true;
}

function cifIsIdentityRotation(R: number[][]): boolean {
    return R[0][0] === 1 && R[0][1] === 0 && R[0][2] === 0 &&
        R[1][0] === 0 && R[1][1] === 1 && R[1][2] === 0 &&
        R[2][0] === 0 && R[2][1] === 0 && R[2][2] === 1;
}

/**
 * The translations that leave the cell's lattice invariant, read from the two
 * sets themselves: a pure translation among a group's operations is a centring
 * vector, so no group letter has to be trusted.
 */
function cifCenteringVectors(groups: ExactCifSymOp[][]): number[][] {
    const out: number[][] = [[0, 0, 0]];
    const seen: { [key: string]: true } = { '0,0,0': true };
    for (const group of groups) {
        for (const op of group) {
            if (!cifIsIdentityRotation(op.R)) continue;
            const key = op.t.join(',');
            if (seen[key]) continue;
            seen[key] = true;
            out.push(op.t.slice());
        }
    }
    return out;
}

function cifOpsMatch(from: ExactCifSymOp[], to: ExactCifSymOp[], centering: number[][]): boolean {
    for (const f of from) {
        let hit = false;
        for (const t of to) {
            if (!cifRotationsEqual(f.R, t.R)) continue;
            const d = [
                ((f.t[0] - t.t[0]) % CIF_SYM_UNITS + CIF_SYM_UNITS) % CIF_SYM_UNITS,
                ((f.t[1] - t.t[1]) % CIF_SYM_UNITS + CIF_SYM_UNITS) % CIF_SYM_UNITS,
                ((f.t[2] - t.t[2]) % CIF_SYM_UNITS + CIF_SYM_UNITS) % CIF_SYM_UNITS
            ];
            for (const c of centering) {
                if (d[0] === c[0] && d[1] === c[1] && d[2] === c[2]) { hit = true; break; }
            }
            if (hit) break;
        }
        if (!hit) return false;
    }
    return true;
}

/**
 * Whether two operation sets describe the same group up to lattice centring.
 * Both directions are checked, so neither a missing nor a spurious operation
 * can pass; the identity/centring translations of either set are the tolerance.
 */
function symOpSetsEquivalent(a: SymOp[], b: SymOp[]): boolean {
    const ea = cifExactOps(a);
    const eb = cifExactOps(b);
    if (ea === null || eb === null) return false;
    const centering = cifCenteringVectors([ea, eb]);
    return cifOpsMatch(ea, eb, centering) && cifOpsMatch(eb, ea, centering);
}

function symmetryCheckWarning(check: CrystalSymmetryCheck): string {
    const label = check.symbol + ' (No. ' + check.number +
        (check.settingLabel === '' ? '' : ', ' + check.settingLabel) + ')';
    return 'space-group check: the file names ' + label + ', whose ' + check.tableOps.length +
        ' operation' + (check.tableOps.length === 1 ? '' : 's') + ' expand the sites to ' +
        check.tableAtoms + ' atom' + (check.tableAtoms === 1 ? '' : 's') + ', but its own ' +
        check.fileOps.length + ' operation' + (check.fileOps.length === 1 ? '' : 's') +
        ' expand them to ' + check.fileAtoms + ' atom' + (check.fileAtoms === 1 ? '' : 's') +
        '; the file\'s operations are used, and the crystal panel offers a switch to the name\'s';
}

// ===== The file's own bond table =====
//
// `_geom_bond` lists the bonds the depositor refined, in the site labels of the
// asymmetric unit, with a symmetry code on each end that names the operation
// and the cell the bonded image sits in (`2_766`: operation 2, translation
// (2, 1, 1)). The table is the file's own statement, and it is worth reading:
// it can disagree with distance-based perception on bond order, and it is the
// only source for bonds whose length falls outside a covalent-radius cutoff.
//
// What comes out is expressed in `baseAtoms` — the expanded cell the viewer
// draws and edits — because that, not the asymmetric unit, is what every page
// pass addresses. A row whose two ends are the same site image is skipped: a
// bond from an atom to its own image has no drawable length or direction, and
// the viewer's bond graph holds one edge per pair anyway.

interface SiteSymmetry {
    op: SymOp;
    shift: number[];
}

/** `.` and a missing code both mean the listed site itself, not operation 1. */
function parseSiteSymmetry(code: string, ops: SymOp[]): SiteSymmetry | null {
    const s = String(code === undefined || code === null ? '' : code).trim();
    if (isMissingValue(s)) {
        return { op: parseSymOp('x, y, z'), shift: [0, 0, 0] };
    }
    const m = /^(-?\d+)_([0-9])([0-9])([0-9])$/.exec(s);
    if (!m) return null;
    const n = parseInt(m[1], 10);
    if (!(n >= 1 && n <= ops.length)) return null;
    return { op: ops[n - 1], shift: [parseInt(m[2], 10) - 5, parseInt(m[3], 10) - 5, parseInt(m[4], 10) - 5] };
}

function bondTableLoop(data: CifData): CifLoop | null {
    for (const loop of data.loops) {
        if (columnOf(loop, ['_geom_bond_atom_site_label_1']) >= 0 &&
            columnOf(loop, ['_geom_bond_atom_site_label_2']) >= 0) return loop;
    }
    return null;
}

type BondEndpoint = { baseIdx: number; offset: number[] } | 'label' | 'code' | 'image';

function readCifBondTable(
    data: CifData,
    parsedSymOps: SymOp[],
    asymAtoms: AsymAtom[],
    baseAtoms: Atom[],
    lattice: number[][]
): { bonds: Bond[]; warnings: string[] } {
    const warnings: string[] = [];
    const loop = bondTableLoop(data);
    if (!loop) {
        // A loop carrying `_geom_bond_` keys but no label pair is a table the
        // reader cannot use. Saying so beats reading the structure as if its
        // file had never claimed any bonds.
        const stray = data.loops.some(l => l.keys.some(k => normTag(k).indexOf('_geom_bond_') === 0));
        if (stray) warnings.push('a _geom_bond loop has no atom-site label columns: its bond table was ignored');
        return { bonds: [], warnings };
    }

    const col1 = columnOf(loop, ['_geom_bond_atom_site_label_1']);
    const col2 = columnOf(loop, ['_geom_bond_atom_site_label_2']);
    const code1Col = columnOf(loop, ['_geom_bond_site_symmetry_1']);
    const code2Col = columnOf(loop, ['_geom_bond_site_symmetry_2']);
    const orderCol = columnOf(loop, ['_geom_bond_order']);
    const valenceCol = columnOf(loop, ['_geom_bond_valence']);
    const minCols = Math.max(col1, col2) + 1;

    const siteOf: { [label: string]: AsymAtom } = {};
    let duplicateLabels = 0;
    let duplicateExample = '';
    for (const site of asymAtoms) {
        const label = (site.label || '').trim();
        if (label === '') continue;
        if (siteOf[label]) {
            duplicateLabels++;
            if (duplicateExample === '') duplicateExample = label;
            continue;
        }
        siteOf[label] = site;
    }

    // Every image of a site as the cell holds it, with the fractional position
    // the image search compares against. A site on a special position has
    // fewer images than the operation list, and the merged-away operations all
    // land on the images that did survive — which is why a bond row may name
    // any operation the file lists.
    const imagesBySite: { [idx: number]: { baseIdx: number; frac: number[] }[] } = {};
    for (const at of baseAtoms) {
        const idx = at.asymIdx;
        if (idx === undefined) continue;
        const baseIdx = at.baseIdx !== undefined ? at.baseIdx : at.index;
        const frac = cartToFrac([at.x, at.y, at.z], lattice);
        (imagesBySite[idx] || (imagesBySite[idx] = [])).push({ baseIdx, frac });
    }

    function resolve(label: string, code: string | undefined): BondEndpoint {
        const site = siteOf[label];
        if (!site) return 'label';
        const sym = parseSiteSymmetry(code === undefined ? '' : code, parsedSymOps);
        if (!sym) return 'code';
        const moved = applySymOp(sym.op, site.frac);
        const target = [moved[0] + sym.shift[0], moved[1] + sym.shift[1], moved[2] + sym.shift[2]];
        const images = imagesBySite[site.index] || [];
        for (const img of images) {
            const offset = [
                Math.round(target[0] - img.frac[0]),
                Math.round(target[1] - img.frac[1]),
                Math.round(target[2] - img.frac[2])
            ];
            let residual = 0;
            for (let axis = 0; axis < 3; axis++) {
                const r = Math.abs(target[axis] - img.frac[axis] - offset[axis]);
                if (r > residual) residual = r;
            }
            // The target is a position this reader computed from the file's own
            // operations, so a real match is round-off; a tenth of a percent of
            // a cell is far below the separation of any two distinct images.
            if (residual < 1e-3) return { baseIdx: img.baseIdx, offset };
        }
        return 'image';
    }

    const bonds: Bond[] = [];
    const seen: { [key: string]: true } = {};
    let badLabels = 0;
    let labelExample = '';
    let badCodes = 0;
    let codeExample = '';
    let badImages = 0;
    let selfBonds = 0;

    for (const row of loop.rows) {
        if (row.length < minCols) continue;
        const label1 = (row[col1] || '').trim();
        const label2 = (row[col2] || '').trim();
        if (isMissingValue(label1) || isMissingValue(label2)) {
            badLabels++;
            if (labelExample === '') labelExample = label1 || label2;
            continue;
        }
        const code1 = code1Col >= 0 ? row[code1Col] : '';
        const code2 = code2Col >= 0 ? row[code2Col] : '';
        const end1 = resolve(label1, code1);
        const end2 = resolve(label2, code2);
        if (end1 === 'label' || end2 === 'label') {
            badLabels++;
            if (labelExample === '') labelExample = end1 === 'label' ? label1 : label2;
            continue;
        }
        if (end1 === 'code' || end2 === 'code') {
            badCodes++;
            if (codeExample === '') codeExample = (end1 === 'code' ? code1 : code2).trim();
            continue;
        }
        if (end1 === 'image' || end2 === 'image') {
            badImages++;
            continue;
        }
        if (end1.baseIdx === end2.baseIdx) {
            selfBonds++;
            continue;
        }
        let a1 = end1.baseIdx;
        let a2 = end2.baseIdx;
        // The cell stores one edge per pair, atom1 < atom2, so the shift is
        // stated for that direction: the image of atom2 in cell c+shift bonds
        // with the image of atom1 in cell c.
        let offset = [
            end2.offset[0] - end1.offset[0],
            end2.offset[1] - end1.offset[1],
            end2.offset[2] - end1.offset[2]
        ];
        if (a1 > a2) {
            const t = a1; a1 = a2; a2 = t;
            offset = [-offset[0], -offset[1], -offset[2]];
        }
        const key = a1 + '_' + a2 + '_' + offset[0] + '_' + offset[1] + '_' + offset[2];
        if (seen[key]) continue;
        seen[key] = true;

        const orderRaw = orderCol >= 0 ? parseNumber(row[orderCol]) : null;
        // A file that lists no order column has not said the bond is single; it
        // has said nothing about the order. Zero means exactly that, so the
        // merge only overrides perception when the table really states one.
        const order = orderRaw !== null && orderRaw >= 1 && orderRaw <= 3 ? orderRaw : 0;
        const valenceRaw = valenceCol >= 0 ? parseNumber(row[valenceCol]) : null;
        const bond: Bond = { atom1: a1, atom2: a2, order, shift: [offset[0], offset[1], offset[2]] };
        if (valenceRaw !== null && valenceRaw > 0) bond.valence = valenceRaw;
        bonds.push(bond);
    }

    if (duplicateLabels > 0) {
        warnings.push('bond table: site label ' + duplicateExample + ' names more than one site: its bonds use the first');
    }
    if (badLabels > 0) {
        warnings.push('bond table: ' + badLabels + ' row' + (badLabels === 1 ? '' : 's') +
            ' name a site the atom list does not have' + (labelExample ? ' (e.g. ' + labelExample + ')' : ''));
    }
    if (badCodes > 0) {
        warnings.push('bond table: ' + badCodes + ' row' + (badCodes === 1 ? '' : 's') +
            ' carry a symmetry code that names no operation in the file' +
            (codeExample ? ' (e.g. ' + codeExample + ')' : ''));
    }
    if (badImages > 0) {
        warnings.push('bond table: ' + badImages + ' row' + (badImages === 1 ? '' : 's') +
            ' name an image of a site that is not in the expanded cell');
    }
    if (selfBonds > 0) {
        warnings.push('bond table: ' + selfBonds + ' bond' + (selfBonds === 1 ? '' : 's') +
            ' join a site to its own symmetry image and were skipped');
    }
    return { bonds, warnings };
}

// ===== Anisotropic displacement parameters =====
//
// `_atom_site_aniso_U_ij` is a table of its own, keyed by site label, and the
// six numbers are a Cartesian tensor in Angstrom squared. The loop is read
// separately from the atom sites because the two need not list the same rows
// (a file may refine ADPs for the heavy atoms only), and it is joined back by
// the label the dictionary itself defines as the join key.
//
// A U_iso-only loop is not this table: it states one isotropic value per site,
// which is not an ellipsoid, and drawing it as one would invent an orientation.
// A record missing a diagonal element is skipped rather than completed from
// the other two — the tensor would no longer be the one the file refined.
interface AnisoTable {
    byLabel: { [label: string]: AnisoTensor };
    warnings: string[];
}

function readAnisoTable(data: CifData): AnisoTable {
    const byLabel: { [label: string]: AnisoTensor } = {};
    const warnings: string[] = [];
    let skipped = 0;
    let skippedExample = '';
    let zeroedOffDiagonal = 0;
    let negative = 0;
    let negativeExample = '';

    for (const loop of data.loops) {
        const labelCol = columnOf(loop, ['_atom_site_aniso_label']);
        if (labelCol < 0) continue;
        const diagCols = ['_atom_site_aniso_U_11', '_atom_site_aniso_U_22', '_atom_site_aniso_U_33']
            .map(tag => columnOf(loop, [tag]));
        if (diagCols.some(col => col < 0)) continue;
        const offCols = ['_atom_site_aniso_U_12', '_atom_site_aniso_U_13', '_atom_site_aniso_U_23']
            .map(tag => columnOf(loop, [tag]));
        const rowWidth = Math.max(labelCol, ...diagCols, ...offCols) + 1;

        for (const row of loop.rows) {
            if (row.length < rowWidth) {
                skipped++;
                if (skippedExample === '') skippedExample = 'a short row';
                continue;
            }
            const label = (row[labelCol] || '').trim();
            if (label === '' || isMissingValue(label)) { skipped++; continue; }
            const d0 = parseNumber(row[diagCols[0]]);
            const d1 = parseNumber(row[diagCols[1]]);
            const d2 = parseNumber(row[diagCols[2]]);
            if (d0 === null || d1 === null || d2 === null) {
                skipped++;
                if (skippedExample === '') skippedExample = 'site ' + label;
                continue;
            }
            // A symmetry-constrained tensor often writes the off-diagonals as
            // `?`/`.`; they are zero by symmetry, and reading them as zero is
            // what the file means, not a guess.
            const off = offCols.map(col => col >= 0 ? parseNumber(row[col]) : null);
            if (off.some(v => v === null)) zeroedOffDiagonal++;
            const u: AnisoTensor = [d0, d1, d2, off[0] ?? 0, off[1] ?? 0, off[2] ?? 0];
            if (d0 < 0 || d1 < 0 || d2 < 0) {
                negative++;
                if (negativeExample === '') {
                    const pairs: [string, number][] = [['U_11', d0], ['U_22', d1], ['U_33', d2]];
                    const bad = pairs.filter(p => p[1] < 0)[0];
                    negativeExample = 'site ' + label + ': ' + bad[0] + '=' + bad[1].toFixed(4);
                }
            }
            // The first record for a label wins: a duplicated row is a file
            // error, and taking the later one would make two readers disagree
            // depending on which they happened to look at.
            if (!byLabel[label]) byLabel[label] = u;
        }
    }

    if (skipped > 0) {
        warnings.push('anisotropic displacement: ' + skipped + ' record' + (skipped === 1 ? '' : 's') +
            ' skipped: ' + (skippedExample === 'a short row' ? 'the row is short' : 'a diagonal value is missing') +
            (skippedExample && skippedExample !== 'a short row' ? ' (' + skippedExample + ')' : ''));
    }
    if (zeroedOffDiagonal > 0) {
        warnings.push('anisotropic displacement: ' + zeroedOffDiagonal + ' record' +
            (zeroedOffDiagonal === 1 ? '' : 's') + ' lack off-diagonal U values, read as 0');
    }
    if (negative > 0) {
        warnings.push('anisotropic displacement: ' + negative + ' record' + (negative === 1 ? ' has' : 's have') +
            ' a negative diagonal element (e.g. ' + negativeExample +
            '): the ellipsoid is drawn from the data that is there');
    }
    return { byLabel, warnings };
}

export function parseCif(content: string): MolecularData {
    const blocks = tokenize(content);
    const warnings: string[] = [];

    // Several data blocks mean several structures. Merging them lets block 2's
    // cell lengths overwrite block 1's and appends its sites to block 1's loop,
    // so one file would be read as a structure that no block describes. Read
    // the first block that carries atom sites and say the others exist.
    let data = blocks[0];
    if (blocks.length > 1) {
        const withSites = blocks.findIndex(block => {
            const site = atomSiteLoop(block);
            return site !== null && site.loop.rows.length > 0;
        });
        const chosen = withSites >= 0 ? withSites : 0;
        data = blocks[chosen];
        const label = (b: CifBlock, k: number) => b.name === '' ? 'block ' + (k + 1) : "'" + b.name + "'";
        const others = blocks.length - 1;
        warnings.push('file holds ' + blocks.length + ' data blocks: reading ' + label(data, chosen) +
            ' (' + others + ' other' + (others === 1 ? '' : 's') + ' ignored)');
    }

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
    let symmetryOps = scan.ops.length > 0 ? scan.ops.slice() : ['x, y, z'];
    // Provenance travels with the operations: the page can then say that the
    // expansion came from a name rather than from the file's own table.
    let symmetrySource: 'file' | 'builtin' | undefined = scan.ops.length > 0 ? 'file' : undefined;

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
            // The file states no operations at all. If it names a group, the
            // built-in table can supply them; a file with its own (even
            // unreadable) operation tag never reaches this branch.
            const hit = spaceGroupOut === undefined ? null : lookupSpaceGroupOps(spaceGroupOut, cell);
            if (hit !== null && hit.number > 1) {
                symmetryOps = hit.ops.slice();
                symmetrySource = 'builtin';
                warnings.push(builtInSymmetryWarning(hit));
            } else if (spaceGroupOut !== undefined && hit === null) {
                warnings.push('no symmetry operations in the file: read as P1 (declared space group ' +
                    spaceGroupOut + ', which is not in the built-in space-group table)');
            } else {
                warnings.push('no symmetry operations in the file: read as P1' +
                    (spaceGroupOut ? ' (declared space group ' + spaceGroupOut + ')' : ''));
            }
        }
    }

    const lattice = cellToLatticeVectors(a, b, c, alpha, beta, gamma);

    // ----- sites -----
    const site = atomSiteLoop(data);
    if (!site) {
        warnings.push('no atom sites in the file');
        return emptyCellCif('No atom sites found in CIF', warnings, cell, spaceGroupOut, symmetryOps, symmetrySource);
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
    const disorderCol = columnOf(atomLoop, ['_atom_site_disorder_group']);

    const aniso = readAnisoTable(data);
    warnings.push(...aniso.warnings);

    const rawAtoms: CifRawAtom[] = [];

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
        const groupRaw = disorderCol >= 0 ? (row[disorderCol] || '').trim() : '';
        // `?` and `.` both mean "no group": a file that marks none must not
        // grow a group called "?" that the panel would offer to filter.
        const disorderGroup = groupRaw && !isMissingValue(groupRaw) ? groupRaw : undefined;

        const frac = useFract ? [xv, yv, zv] : cartToFrac([xv, yv, zv], lattice);
        // The join is by site label, the dictionary's own key. A site whose
        // label is empty cannot be named by the aniso table and keeps none.
        const uAniso = label.trim() !== '' ? aniso.byLabel[label.trim()] : undefined;
        rawAtoms.push({ element, label, frac, occupancy, disorderGroup, uAniso });
    }
    if (skipped > 0) {
        warnings.push(skipped + ' atom site' + (skipped === 1 ? '' : 's') +
            ' skipped: the row is short or its coordinates are not numbers');
    }
    const anisoLabels = Object.keys(aniso.byLabel);
    if (anisoLabels.length > 0) {
        const siteLabels: { [label: string]: true } = {};
        for (const raw of rawAtoms) {
            const l = raw.label.trim();
            if (l !== '') siteLabels[l] = true;
        }
        const unmatched = anisoLabels.filter(l => !siteLabels[l]);
        if (unmatched.length > 0) {
            warnings.push('anisotropic displacement: ' + unmatched.length + ' label' +
                (unmatched.length === 1 ? '' : 's') +
                ' name a site the atom list does not have (e.g. ' + unmatched[0] + ')');
        }
    }

    const parsedSymOps = symmetryOps.map(parseSymOp);

    // ----- expand -----
    const expansion = expandCifSites(rawAtoms, parsedSymOps, lattice);
    const baseAtoms = expansion.baseAtoms;
    const asymAtoms = expansion.asymAtoms;
    const nearMisses = expansion.nearMisses;

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
        cellI: 0, cellJ: 0, cellK: 0,
        uAniso: at.uAniso
    }));

    // The file's bond table is read last: it needs the expanded cell (to name
    // the image each end sits in) and the asymmetric unit (to know what its
    // site labels mean). Its warnings join the file's own, so the panel's one
    // summary line accounts for them.
    const bondTable = readCifBondTable(data, parsedSymOps, asymAtoms, baseAtoms, lattice);
    warnings.push(...bondTable.warnings);

    // The file's operations are what the picture is built from — that does not
    // change. When the file also names a group the built-in table knows, the
    // two statements are compared once; a disagreement is reported with the
    // atom count of each expansion, and the name's expansion is kept so the
    // page can offer a one-time switch.
    let symmetryCheck: CrystalSymmetryCheck | undefined;
    let symmetryAlternative: CrystalSymmetryAlternative | undefined;
    if (scan.ops.length > 0 && spaceGroupOut !== undefined) {
        const hit = lookupSpaceGroupOps(spaceGroupOut, cell);
        if (hit !== null) {
            const nameSymOps = hit.ops.map(parseSymOp);
            if (!symOpSetsEquivalent(parsedSymOps, nameSymOps)) {
                const nameExpansion = expandCifSites(rawAtoms, nameSymOps, lattice);
                symmetryCheck = {
                    declared: spaceGroupOut,
                    symbol: hit.symbol,
                    settingLabel: hit.settingLabel,
                    number: hit.number,
                    fileOps: symmetryOps.slice(),
                    tableOps: hit.ops.slice(),
                    fileAtoms: baseAtoms.length,
                    tableAtoms: nameExpansion.baseAtoms.length
                };
                symmetryAlternative = {
                    ops: hit.ops.slice(),
                    symbol: hit.symbol,
                    settingLabel: hit.settingLabel,
                    baseAtoms: nameExpansion.baseAtoms,
                    asymAtoms: nameExpansion.asymAtoms
                };
                // First, because the panel's warning line shows only the head
                // of the list and a self-contradicting file is the one thing
                // here that changes what the structure is.
                warnings.unshift(symmetryCheckWarning(symmetryCheck));
            }
        }
    }

    const crystal: CrystalData = {
        a, b, c, alpha, beta, gamma,
        latticeVectors: lattice,
        spaceGroup: spaceGroupOut,
        symmetryOps,
        symmetrySource,
        symmetryCheck,
        symmetryAlternative,
        baseAtoms: baseAtoms.map(at => ({ ...at, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: [] as Bond[],
        asymAtoms,
        fileBonds: bondTable.bonds.length > 0 ? bondTable.bonds : undefined
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

// The two conversions are part of this module's published surface (callers use
// them to move between the two coordinate systems a CIF can list), but the
// arithmetic is shared with the other crystal readers.
export { cartToFrac as cartesianToFractional, fracToCart as fractionalToCartesian };
