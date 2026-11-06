import { MolecularData } from '../types';
import { elementFromZ, normalizeElementSymbol } from '../tools/periodicTable';
import { BOHR_TO_ANGSTROM, buildMolecularStructure, buildPeriodicStructure, isUsableLattice, Site } from './periodicUtils';
import { parseXyzTrajectory } from './xyzParser';

/**
 * ABINIT reader, for both the input (`.abi`) and the output (`.abo`).
 *
 * ABINIT writes its geometry as three ingredients that never appear in the
 * same place twice: `acell` (three lengths), `rprim` (dimensionless primitive
 * vectors), and positions as `xred` (fractional), `xcart` (Bohr) or `xangst`
 * (Angstrom). The output repeats all of them under `-outvars` and, after a
 * relaxation, prints the *final* cell and coordinates in a different block
 * again. The reader therefore has one tokenizer and two entry points:
 *
 *  - the input uses the tokenizer on the file;
 *  - the output uses it on the `-outvars` echo as a fallback, then prefers the
 *    last `R(1)=...` lattice and the `cartesian coordinates (angstrom) at end`
 *    table, which are what the run actually converged to.
 *
 * A cell is built from `rprim` rows scaled by the matching `acell` entry, not
 * by one common length: ABINIT allows three different `acell` values and using
 * only the first would silently sheare an orthorhombic cell.
 */

interface TokenizedInput {
    /** Keyword -> every value token after it, in file order. */
    values: { [key: string]: string[] };
    /** Keyword -> position in the file, for the xred/xcart/xangst choice. */
    order: { [key: string]: number };
    /** The file's own leading comments, joined. */
    comments: string[];
}

const POSITION_KEYWORDS = ['xangst', 'xcart', 'xred'];

function isKeywordToken(token: string): boolean {
    return /^[A-Za-z][A-Za-z0-9_]*$/.test(token) && !/^(bohr|angstrom|ang|nm|hartree|ev|ry|rydberg)$/i.test(token);
}

/**
 * Tokenize the simple `keyword value...` grammar ABINIT uses.
 *
 * A continuation line holds only values, so a line whose first token is not a
 * number is a new keyword — except unit words (`Bohr`, `Angstrom`), which are
 * values even though they are alphabetic. Comments (`#`, `!`) are dropped, and
 * quoted strings stay one token so a pseudo path cannot look like a keyword.
 */
function tokenizeAbinit(content: string): TokenizedInput {
    const values: { [key: string]: string[] } = {};
    const order: { [key: string]: number } = {};
    const comments: string[] = [];
    let current = '';
    let counter = 0;

    const lines = content.split(/\r?\n/);
    for (const rawLine of lines) {
        const trimmedRaw = rawLine.trim();
        if (trimmedRaw.startsWith('#')) {
            // Collect the comment before its line is stripped away: a full-line
            // comment is file metadata, and the first one titles the viewer.
            const text = trimmedRaw.replace(/^#+\s*/, '');
            if (text !== '') comments.push(text);
            continue;
        }
        let line = rawLine;
        const hash = line.indexOf('#');
        if (hash >= 0) line = line.slice(0, hash);
        const bang = line.indexOf('!');
        if (bang >= 0) line = line.slice(0, bang);
        line = line.trim();
        if (line === '') continue;

        const tokens = line.match(/"[^"]*"|'[^']*'|\S+/g) || [];
        if (tokens.length === 0) continue;
        const first = tokens[0] || '';
        if (isKeywordToken(first)) {
            // A keyword stated again replaces its earlier values, as ABINIT's
            // own parser keeps the last definition: appending would mix a
            // stale geometry into the live one.
            current = first.toLowerCase();
            values[current] = [];
            order[current] = counter++;
            for (let i = 1; i < tokens.length; i++) values[current].push(tokens[i]);
        } else if (current !== '') {
            for (const t of tokens) values[current].push(t);
        }
    }
    return { values, order, comments };
}

/**
 * Expand ABINIT's repeat notation (`3*10.0` = three tens) into its values.
 *
 * The writer uses it for a cell that repeats a number — `acell 3*10.0` is an
 * ordinary way to state a cubic cell — and a reader that parses the token as
 * the number 3 would draw a 3 Bohr cell instead.
 */
function expandRepeats(tokens: string[] | undefined): string[] {
    if (!tokens) return [];
    const out: string[] = [];
    for (const t of tokens) {
        const m = /^(\d+)\*(\S+)$/.exec(t);
        if (m) {
            const n = parseInt(m[1], 10);
            if (isFinite(n) && n > 0 && n <= 100000) {
                for (let i = 0; i < n; i++) out.push(m[2]);
                continue;
            }
        }
        out.push(t);
    }
    return out;
}

/** A number, accepting the Fortran `1.0d-3` exponent ABINIT files may carry. */
function parseAbinitNumber(token: string): number {
    const m = /^([-+]?(?:\d+\.?\d*|\.\d+))[dD]([-+]?\d+)$/.exec(token);
    if (m) return parseFloat(m[1] + 'e' + m[2]);
    return parseFloat(token);
}

function numberList(tokens: string[] | undefined): number[] {
    const out: number[] = [];
    for (const t of expandRepeats(tokens)) {
        const v = parseAbinitNumber(t);
        if (isFinite(v)) out.push(v);
    }
    return out;
}

/**
 * The element for every `znucl` entry, keeping the type indices aligned.
 *
 * A charge with no element (a saved ghost's 0) and an unreadable symbol become
 * the dummy X rather than being dropped: one dropped entry would shift every
 * later `typat` index onto the wrong element, which is a silent relabelling of
 * the whole cell.
 */
function znuclElements(tokens: string[] | undefined): string[] {
    const out: string[] = [];
    for (const t of expandRepeats(tokens)) {
        const z = parseAbinitNumber(t);
        if (isFinite(z)) out.push(elementFromZ(z) || 'X');
        else out.push(normalizeElementSymbol(t) || 'X');
    }
    return out;
}

function unitOf(tokens: string[] | undefined): 'bohr' | 'angstrom' {
    if (!tokens) return 'bohr';
    for (const t of tokens) {
        const s = t.toLowerCase();
        if (s.indexOf('angstrom') >= 0 || s === 'ang') return 'angstrom';
    }
    return 'bohr';
}

function datasetKey(input: TokenizedInput, key: string, dataset: number): string[] | undefined {
    const suffix = key + String(dataset);
    if (dataset > 1 && input.values[suffix]) return input.values[suffix];
    return input.values[key];
}

function datasetOrder(input: TokenizedInput, key: string, dataset: number): number {
    const suffix = key + String(dataset);
    if (dataset > 1 && input.order[suffix] !== undefined) return input.order[suffix];
    return input.order[key] !== undefined ? input.order[key] : -1;
}

interface AbinitGeometry {
    lattice: number[][];
    sites: Site[];
    warnings: string[];
    natom: number;
}

/** The dataset whose geometry the input describes (ABINIT's last one wins). */
function abinitDataset(input: TokenizedInput): number {
    const ndtsetRaw = numberList(input.values['ndtset']);
    return ndtsetRaw.length > 0 && ndtsetRaw[0] > 1 ? Math.floor(ndtsetRaw[0]) : 1;
}

interface AbinitCell {
    lattice: number[][];
    warnings: string[];
    /** Whether the file states a cell rather than relying on ABINIT's 1-Bohr default. */
    stated: boolean;
}

function cellFromInput(input: TokenizedInput, dataset: number): AbinitCell {
    const warnings: string[] = [];
    const acellTokens = datasetKey(input, 'acell', dataset);
    const acellValues = numberList(acellTokens);
    if (acellValues.length === 0) {
        warnings.push('no acell in the ABINIT input: the cell is read in Bohr with acell = 1');
    }
    const unit = unitOf(acellTokens);
    const factor = unit === 'angstrom' ? 1 : BOHR_TO_ANGSTROM;
    let acell = [1, 1, 1];
    if (acellValues.length === 1) acell = [acellValues[0], acellValues[0], acellValues[0]];
    else if (acellValues.length >= 3) acell = [acellValues[0], acellValues[1], acellValues[2]];
    acell = acell.map(v => v * factor);

    const rprimValues = numberList(datasetKey(input, 'rprim', dataset));
    let rprim = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (rprimValues.length >= 9) {
        rprim = [
            [rprimValues[0], rprimValues[1], rprimValues[2]],
            [rprimValues[3], rprimValues[4], rprimValues[5]],
            [rprimValues[6], rprimValues[7], rprimValues[8]]
        ];
    } else if (rprimValues.length > 0) {
        warnings.push('rprim needs nine numbers: the identity is used instead');
    }
    const lattice = rprim.map((row, i) => row.map(v => v * acell[i]));
    return { lattice, warnings, stated: acellValues.length > 0 || rprimValues.length > 0 };
}

function geometryFromInput(input: TokenizedInput): AbinitGeometry | null {
    const dataset = abinitDataset(input);
    const cell = cellFromInput(input, dataset);
    const warnings = cell.warnings;
    const lattice = cell.lattice;

    const natomRaw = numberList(datasetKey(input, 'natom', dataset));
    let natom = natomRaw.length > 0 ? Math.floor(natomRaw[0]) : 0;

    let typat = numberList(datasetKey(input, 'typat', dataset)).map(v => Math.floor(v));
    const znucl = znuclElements(datasetKey(input, 'znucl', dataset));

    // Whichever of xred/xcart/xangst the file states last wins, as ABINIT's
    // own parser does; preferring a fixed member of the family would read a
    // stale geometry when the input overrides one with another.
    let chosen = '';
    let chosenOrder = -1;
    for (const key of POSITION_KEYWORDS) {
        const o = datasetOrder(input, key, dataset);
        if (o >= 0 && o > chosenOrder && numberList(datasetKey(input, key, dataset)).length >= 3) {
            chosen = key;
            chosenOrder = o;
        }
    }
    if (chosen === '') {
        warnings.push('the ABINIT input states no xred, xcart or xangst: no sites to show');
        return { lattice, sites: [], warnings, natom };
    }

    if (natom <= 0) {
        // With no natom, the chosen block alone decides the count; taking it
        // from a stale member of the coordinate family would clip the live one.
        natom = Math.floor(numberList(datasetKey(input, chosen, dataset)).length / 3);
    }
    if (typat.length === 0 && natom > 0) {
        // One type for the whole cell is the only reading that needs no typat.
        typat = new Array(natom).fill(1);
    }

    const positions = numberList(datasetKey(input, chosen, dataset));
    const sites: Site[] = [];
    const count = natom > 0 ? Math.min(natom, Math.floor(positions.length / 3)) : Math.floor(positions.length / 3);
    for (let i = 0; i < count; i++) {
        const type = typat[i] !== undefined ? typat[i] : 1;
        const element = znucl[type - 1] || znucl[0] || 'X';
        const x = positions[3 * i];
        const y = positions[3 * i + 1];
        const z = positions[3 * i + 2];
        if (chosen === 'xred') sites.push({ element, frac: [x, y, z] });
        else if (chosen === 'xcart') sites.push({ element, cart: [x * BOHR_TO_ANGSTROM, y * BOHR_TO_ANGSTROM, z * BOHR_TO_ANGSTROM] });
        else sites.push({ element, cart: [x, y, z] });
    }
    if (znucl.length === 0) warnings.push('no znucl in the ABINIT input: every site is read as the dummy element X');
    if (sites.length === 0) warnings.push('no readable sites in the ABINIT input');
    return { lattice, sites, warnings, natom };
}

export interface AbinitExternalStructure {
    content: string;
    fileName: string;
}

export interface AbinitXyzfileScan {
    /** The first `xyzfile` the input states: the geometry the run starts from. */
    xyzfile?: string;
    /** Every distinct `xyzfile` name, in file order. */
    files: string[];
}

/** The XYZ structure files an ABINIT input names; ABINIT resolves them from its run folder. */
export function scanAbinitXyzfile(content: string): AbinitXyzfileScan {
    const input = tokenizeAbinit(content);
    const entries: { name: string; order: number }[] = [];
    for (const key of Object.keys(input.values)) {
        if (!/^xyzfile\d*$/.test(key)) continue;
        const tokens = input.values[key];
        const raw = tokens.length > 0 ? tokens[0] : '';
        const name = raw.replace(/^['"]|['"]$/g, '').trim();
        if (name === '') continue;
        entries.push({ name, order: input.order[key] !== undefined ? input.order[key] : 0 });
    }
    entries.sort((a, b) => a.order - b.order);
    const files: string[] = [];
    for (const entry of entries) {
        if (files.indexOf(entry.name) < 0) files.push(entry.name);
    }
    return { files, xyzfile: files.length > 0 ? files[0] : undefined };
}

export function parseAbinitInput(content: string, external?: AbinitExternalStructure): MolecularData {
    const input = tokenizeAbinit(content);
    const title = input.comments.length > 0 ? input.comments[0] : 'ABINIT input';
    const xyz = scanAbinitXyzfile(content);

    if (xyz.xyzfile) {
        if (!external) {
            return {
                atoms: [], bonds: [], title, hasExplicitBonds: false,
                warnings: ['the ABINIT input reads its structure from "' + xyz.xyzfile +
                    '", but that file could not be read']
            };
        }
        const trajectory = parseXyzTrajectory(external.content);
        const frame = trajectory.frames.length > 0 ? trajectory.frames[0] : undefined;
        if (!frame || frame.atoms.length === 0) {
            return {
                atoms: [], bonds: [], title, hasExplicitBonds: false,
                warnings: ['xyzfile "' + external.fileName + '": the coordinate file has no atoms']
            };
        }
        const sites: Site[] = frame.atoms.map(atom => {
            const site: Site = { element: atom.element, cart: [atom.x, atom.y, atom.z] };
            if (atom.ghost) site.ghost = true;
            return site;
        });
        const displayTitle = input.comments.length === 0 && frame.title !== '' ? frame.title : title;
        const cell = cellFromInput(input, abinitDataset(input));
        if (!cell.stated) {
            // ABINIT's 1-Bohr default cell would squeeze a molecular xyzfile
            // into a half-Angstrom cube; the molecule is better shown free.
            return buildMolecularStructure(sites, displayTitle, [
                'the ABINIT input states no cell: the xyzfile structure is shown without a cell'
            ]);
        }
        if (!isUsableLattice(cell.lattice)) {
            return {
                atoms: [], bonds: [], title: displayTitle, hasExplicitBonds: false,
                warnings: cell.warnings.concat('the ABINIT input cell is degenerate: no structure to show')
            };
        }
        return buildPeriodicStructure(cell.lattice, sites, displayTitle, cell.warnings);
    }

    const geometry = geometryFromInput(input);
    if (!geometry || geometry.sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry ? geometry.warnings : ['unreadable ABINIT input: no structure to show']
        };
    }
    if (!isUsableLattice(geometry.lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry.warnings.concat('the ABINIT input cell is degenerate: no structure to show')
        };
    }
    return buildPeriodicStructure(geometry.lattice, geometry.sites, title, geometry.warnings);
}

interface OutvarGeometry {
    lattice: number[][] | null;
    natom: number;
    typat: number[];
    znucl: string[];
    sites: Site[];
    warnings: string[];
}

function parseOutvarBlock(content: string): OutvarGeometry {
    const warnings: string[] = [];
    const lines = content.split(/\r?\n/);
    let start = -1;
    for (let i = 0; i < lines.length; i++) {
        // A multi-dataset run echoes the resolved variables once per dataset;
        // the last echo describes the geometry the run finished with.
        if (/-outvars: echo values of preprocessed input variables/i.test(lines[i])) start = i;
    }
    let block = '';
    if (start >= 0) {
        const chunk: string[] = [];
        for (let i = start; i < lines.length; i++) {
            const line = lines[i];
            if (i > start && /^\s*={5,}/.test(line)) break;
            if (i > start && /^-{2,}\s*$/.test(line.trim())) break;
            chunk.push(line);
        }
        block = chunk.join('\n');
    }
    const input = tokenizeAbinit(block);

    const acellValues = numberList(input.values['acell']);
    const acellUnit = unitOf(input.values['acell']);
    const factor = acellUnit === 'angstrom' ? 1 : BOHR_TO_ANGSTROM;
    let acell = [1, 1, 1];
    if (acellValues.length === 1) acell = [acellValues[0], acellValues[0], acellValues[0]];
    else if (acellValues.length >= 3) acell = [acellValues[0], acellValues[1], acellValues[2]];
    acell = acell.map(v => v * factor);

    const rprimValues = numberList(input.values['rprim']);
    let rprim = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (rprimValues.length >= 9) {
        rprim = [
            [rprimValues[0], rprimValues[1], rprimValues[2]],
            [rprimValues[3], rprimValues[4], rprimValues[5]],
            [rprimValues[6], rprimValues[7], rprimValues[8]]
        ];
    }
    let lattice: number[][] | null = rprimValues.length >= 9 || acellValues.length > 0
        ? rprim.map((row, i) => row.map(v => v * acell[i]))
        : null;

    const natomRaw = numberList(input.values['natom']);
    let natom = natomRaw.length > 0 ? Math.floor(natomRaw[0]) : 0;

    const typat = numberList(input.values['typat']).map(v => Math.floor(v));
    const znucl = znuclElements(input.values['znucl']);

    // The echoed xred is the last geometry when no Cartesian table follows.
    let sites: Site[] = [];
    const xred = numberList(input.values['xred']);
    const xcart = numberList(input.values['xcart']);
    const xangst = numberList(input.values['xangst']);
    const positions = xangst.length >= 3 ? { list: xangst, kind: 'angstrom' }
        : xcart.length >= 3 ? { list: xcart, kind: 'bohr' }
            : xred.length >= 3 ? { list: xred, kind: 'frac' } : null;
    if (positions) {
        const n = natom > 0 ? Math.min(natom, Math.floor(positions.list.length / 3)) : Math.floor(positions.list.length / 3);
        if (natom <= 0) natom = n;
        for (let i = 0; i < n; i++) {
            const type = typat[i] !== undefined ? typat[i] : 1;
            const element = znucl[type - 1] || znucl[0] || 'X';
            const x = positions.list[3 * i];
            const y = positions.list[3 * i + 1];
            const z = positions.list[3 * i + 2];
            if (positions.kind === 'frac') sites.push({ element, frac: [x, y, z] });
            else if (positions.kind === 'bohr') sites.push({ element, cart: [x * BOHR_TO_ANGSTROM, y * BOHR_TO_ANGSTROM, z * BOHR_TO_ANGSTROM] });
            else sites.push({ element, cart: [x, y, z] });
        }
    }

    // Final lattice: the last R(1)/R(2)/R(3) group is in Bohr and is the cell
    // the run finished with, which matters for a vc-relax.
    const rVectors: number[][][] = [];
    for (let i = 0; i < lines.length; i++) {
        if (!/space primitive vectors/i.test(lines[i])) continue;
        const rows: number[][] = [];
        for (let j = i + 1; j < lines.length && rows.length < 3; j++) {
            const m = /^\s*R\(\d\)\s*=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)/.exec(lines[j]);
            if (!m) break;
            rows.push([parseFloat(m[1]) * BOHR_TO_ANGSTROM, parseFloat(m[2]) * BOHR_TO_ANGSTROM, parseFloat(m[3]) * BOHR_TO_ANGSTROM]);
        }
        if (rows.length === 3) rVectors.push(rows);
    }
    if (rVectors.length > 0) lattice = rVectors[rVectors.length - 1];

    const finalAcell = /Final cell parameters:\s*acell\s*=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*(\w*)/i.exec(content);
    // The final acell wins over the echoed input cell — but only when no R()
    // table stated the complete final cell, which is the better source.
    if (rVectors.length === 0 && finalAcell) {
        const u = /angstrom/i.test(finalAcell[4] || '') ? 1 : BOHR_TO_ANGSTROM;
        const a = [parseFloat(finalAcell[1]) * u, parseFloat(finalAcell[2]) * u, parseFloat(finalAcell[3]) * u];
        lattice = rprim.map((row, i) => row.map(v => v * a[i]));
    }

    // Final positions: the numbered block after a relaxation is Cartesian
    // Angstrom and supersedes the echoed xred.
    let finalCart: number[][] = [];
    for (let i = lines.length - 1; i >= 0; i--) {
        if (!/cartesian coordinates\s*\((angstrom|bohr)\)\s*at end/i.test(lines[i])) continue;
        const bohr = /bohr/i.test(lines[i]);
        const rows: number[][] = [];
        for (let j = i + 1; j < lines.length; j++) {
            const parts = lines[j].trim().split(/\s+/);
            if (parts.length < 4) break;
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
            const f = bohr ? BOHR_TO_ANGSTROM : 1;
            rows.push([x * f, y * f, z * f]);
        }
        if (rows.length > 0) { finalCart = rows; break; }
    }
    if (finalCart.length > 0) {
        if (natom <= 0) natom = finalCart.length;
        sites = finalCart.slice(0, natom).map((row, i) => {
            const type = typat[i] !== undefined ? typat[i] : 1;
            return {
                element: znucl[type - 1] || znucl[0] || 'X',
                cart: row
            } as Site;
        });
    }

    if (!lattice) warnings.push('the ABINIT output states no cell: no structure to show');
    if (sites.length === 0) warnings.push('the ABINIT output states no coordinates: no structure to show');
    if (znucl.length === 0 && sites.length > 0) warnings.push('the ABINIT output names no elements: every site is read as the dummy element X');
    return { lattice, natom, typat, znucl, sites, warnings };
}

export function parseAbinitOutput(content: string): MolecularData {
    const geometry = parseOutvarBlock(content);
    const titleMatch = /^\s*\.ABINIT\s+([^\s]+)/m.exec(content);
    const title = titleMatch ? 'ABINIT ' + titleMatch[1] : 'ABINIT output';
    if (!geometry.lattice || geometry.sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry.warnings
        };
    }
    if (!isUsableLattice(geometry.lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry.warnings.concat('the ABINIT output cell is degenerate: no structure to show')
        };
    }
    const warnings = geometry.warnings.slice();
    return buildPeriodicStructure(geometry.lattice, geometry.sites, title, warnings);
}

export function isAbinitOutput(content: string): boolean {
    return /-outvars: echo values of preprocessed input variables/i.test(content) ||
        /^\s*\.ABINIT\b/m.test(content);
}

export function isAbinitInput(content: string): boolean {
    // An output echoes its input variables (`natom`, `ecut`, ...) in the
    // `-outvars` block, so the output banner has to win before the key scan.
    if (isAbinitOutput(content)) return false;
    return /^\s*(acell|rprim|xred|xcart|xangst|natom|typat|znucl|xyzfile\d*|ecut|ndtset|pseudos|kptopt|ngkpt)\b/m.test(content) ||
        /^\s*(natom|ntypat)\s+\d+/m.test(content);
}
