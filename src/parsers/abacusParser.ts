import { MolecularData } from '../types';
import { elementFromAtomName, normalizeElementSymbol } from '../tools/periodicTable';
import { BOHR_TO_ANGSTROM, buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * ABACUS reader: the input structure (`STRU`), the SCF log and the MD dump.
 *
 * The three files describe the same cell in three layouts:
 *
 *  - `STRU` is section oriented: `ATOMIC_SPECIES`, `LATTICE_CONSTANT`,
 *    `LATTICE_VECTORS`, then `ATOMIC_POSITIONS` with a mode line and one
 *    element block per species (element, magnetization, count, sites).
 *  - `running_scf.log` echoes the cell it read: lattice constant, lattice
 *    vectors in units of a0, an atom-type table and a `DIRECT COORDINATES`
 *    table.
 *  - `MD_dump` repeats the whole cell once per `MDSTEP`, with positions in
 *    Angstrom and force/velocity columns after them.
 *
 * Units are the trap in all three: `Cartesian` means *multiples of
 * LATTICE_CONSTANT* in a STRU, `Cartesian_angstrom` means Angstrom, and the log
 * states its lattice in a0 while the dump states it in Angstrom. Every one is
 * converted here; a STRU read with the wrong unit would draw a cell thousands
 * of Angstroms across and nothing downstream could tell it was wrong.
 */

const STRU_SECTIONS = /^(ATOMIC_SPECIES|NUMERICAL_ORBITAL|NUMERICAL_DESCRIPTOR|LATTICE_CONSTANT|LATTICE_VECTORS|ATOMIC_POSITIONS|K_POINTS|PAW|ORBITAL|DFT_ORBITAL|EXTERNAL_ORBITAL|EXTERNAL_CHARGE|INIT_CHARGE|ATOMIC_MAGMOM|ATOMIC_VELOCITY|ATOMIC_FORCE|CONSTRAINT|MD_INPUT)\b/;

/** A value line whose leading token is a number and whose unit words are known. */
function parseLengthUnit(tokens: string[], defaultUnit: 'bohr' | 'angstrom'): { value: number; unit: 'bohr' | 'angstrom' } | null {
    const value = parseFloat(tokens[0]);
    if (!isFinite(value)) return null;
    let unit = defaultUnit;
    for (let i = 1; i < tokens.length; i++) {
        const t = tokens[i].toLowerCase();
        if (t.indexOf('angstrom') >= 0 || t === 'a') unit = 'angstrom';
        else if (t.indexOf('bohr') >= 0 || t === 'au' || t === 'a.u.') unit = 'bohr';
    }
    return { value, unit };
}

function toAngstrom(value: number, unit: 'bohr' | 'angstrom'): number {
    return unit === 'bohr' ? value * BOHR_TO_ANGSTROM : value;
}

/** Remove ABACUS's two comment styles without touching a `#` inside a path. */
function stripComment(line: string): string {
    let out = line;
    const slash = out.indexOf('//');
    if (slash >= 0) out = out.slice(0, slash);
    const hash = out.indexOf('#');
    if (hash >= 0) out = out.slice(0, hash);
    const bang = out.indexOf('!');
    if (bang >= 0) out = out.slice(0, bang);
    return out.trim();
}

function readNumbers(lines: string[], start: number, count: number): number[] {
    const out: number[] = [];
    for (let i = start; i < lines.length && out.length < count; i++) {
        const tokens = stripComment(lines[i]).split(/\s+/).filter(t => t !== '');
        if (tokens.length === 0) continue;
        for (const t of tokens) {
            const v = parseFloat(t);
            if (!isFinite(v)) return out;
            out.push(v);
            if (out.length >= count) break;
        }
    }
    return out;
}

interface StruPositionMode {
    kind: 'direct' | 'cartesian';
    unit: 'angstrom' | 'bohr' | 'lattice';
    /** Axes a Cartesian mode centres on, e.g. `xy`; empty for a corner origin. */
    center: string;
}

function parsePositionMode(line: string): StruPositionMode | null {
    // ABACUS writes modes as `Direct`, `Cartesian_bohr` or the centered
    // `Cartesian_angstrom_center_xy`, so the spelling is matched by prefix
    // rather than by equality: an unknown suffix would otherwise turn every
    // Cartesian coordinate into a fractional one.
    const t = line.trim().toLowerCase().replace(/[\s_-]/g, '');
    if (t === 'direct' || t === 'fractional') return { kind: 'direct', unit: 'angstrom', center: '' };
    if (t !== 'angstrom' && t !== 'bohr' && !t.startsWith('cartesian')) return null;
    const unit = t.indexOf('bohr') >= 0 ? 'bohr'
        : t === 'angstrom' || t.indexOf('angstrom') >= 0 || t.indexOf('cartesianang') === 0 ? 'angstrom'
        : 'lattice';
    const centered = /center([xyz]+)$/.exec(t);
    return { kind: 'cartesian', unit, center: centered ? centered[1] : '' };
}

export function parseAbacusStru(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const title = 'ABACUS STRU';

    let latticeConstant = 1.0;
    let lattice: number[][] | null = null;
    let mode: StruPositionMode | null = null;
    const sites: Site[] = [];

    let i = 0;
    while (i < lines.length) {
        const raw = stripComment(lines[i]);
        if (raw === '') { i++; continue; }
        const head = raw.split(/\s+/)[0].toUpperCase();

        if (head === 'LATTICE_CONSTANT') {
            i++;
            while (i < lines.length && stripComment(lines[i]) === '') i++;
            if (i < lines.length) {
                const tokens = stripComment(lines[i]).split(/\s+/);
                const parsed = parseLengthUnit(tokens, 'bohr');
                if (parsed) {
                    latticeConstant = toAngstrom(parsed.value, parsed.unit);
                } else {
                    warnings.push('LATTICE_CONSTANT is not a number: read as 1.0');
                }
            }
            i++;
            continue;
        }

        if (head === 'LATTICE_VECTORS') {
            i++;
            const values = readNumbers(lines, i, 9);
            if (values.length === 9) {
                lattice = [
                    [values[0], values[1], values[2]],
                    [values[3], values[4], values[5]],
                    [values[6], values[7], values[8]]
                ];
                // The STRU states the vectors dimensionless and scales them by
                // LATTICE_CONSTANT; the constant was already converted to
                // Angstrom above, so the product is Angstrom.
                for (let r = 0; r < 3; r++) {
                    for (let c = 0; c < 3; c++) lattice[r][c] *= latticeConstant;
                }
                // One line may hold all nine numbers; skip only as many lines as
                // those numbers actually occupied.
                let need = 9;
                while (i < lines.length && need > 0) {
                    const tokens = stripComment(lines[i]).split(/\s+/).filter(t => t !== '');
                    if (tokens.length === 0) { i++; continue; }
                    need -= tokens.length;
                    i++;
                }
            } else {
                warnings.push('LATTICE_VECTORS is incomplete: no unit cell to show');
            }
            continue;
        }

        if (head === 'ATOMIC_SPECIES') {
            i++;
            while (i < lines.length) {
                const line = stripComment(lines[i]);
                if (line === '') { i++; continue; }
                if (STRU_SECTIONS.test(line.split(/\s+/)[0].toUpperCase() + ' ')) break;
                const tokens = line.split(/\s+/);
                if (tokens.length >= 1 && /^[A-Za-z]/.test(tokens[0])) {
                    // The species table is not used to place atoms — the
                    // ATOMIC_POSITIONS blocks name their element themselves —
                    // so only the symbol's spelling is checked here.
                    normalizeElementSymbol(tokens[0]);
                }
                i++;
            }
            continue;
        }

        if (head === 'ATOMIC_POSITIONS') {
            // The mode may share the keyword's line (`ATOMIC_POSITIONS Direct`).
            const rest = raw.split(/\s+/).slice(1).join(' ');
            if (rest !== '') mode = parsePositionMode(rest);
            i++;
            while (i < lines.length && stripComment(lines[i]) === '') i++;
            if (mode === null && i < lines.length) {
                const possible = parsePositionMode(stripComment(lines[i]));
                if (possible) { mode = possible; i++; }
            }
            if (mode === null) {
                mode = { kind: 'direct', unit: 'angstrom', center: '' };
                warnings.push('ATOMIC_POSITIONS states no mode: coordinates read as Direct');
            }

            while (i < lines.length) {
                const line = stripComment(lines[i]);
                if (line === '') { i++; continue; }
                const tokens = line.split(/\s+/);
                if (STRU_SECTIONS.test(tokens[0].toUpperCase() + ' ')) break;
                if (!/^[A-Za-z]/.test(tokens[0]) || tokens.length > 1) break;

                const element = normalizeElementSymbol(tokens[0]) || tokens[0];
                i++;
                while (i < lines.length && stripComment(lines[i]) === '') i++;

                // Standard STRU: element / magnetization / count / coordinates.
                // A file that omits the magnetization line is still readable:
                // the count line is a lone integer, while the line after a
                // magnetization is either the count or — when the magnetization
                // was omitted after all — the first coordinate row.
                const singleNumber = (text: string): number | null => {
                    const parts = stripComment(text).split(/\s+/).filter(t => t !== '');
                    if (parts.length !== 1) return null;
                    const v = parseFloat(parts[0]);
                    return isFinite(v) ? v : null;
                };
                let count = NaN;
                if (i < lines.length) {
                    const first = singleNumber(lines[i]);
                    if (first !== null) {
                        const firstText = stripComment(lines[i]);
                        const decimal = /[.eE]/.test(firstText);
                        const nextIdx = (() => {
                            let j = i + 1;
                            while (j < lines.length && stripComment(lines[j]) === '') j++;
                            return j;
                        })();
                        const second = nextIdx < lines.length ? singleNumber(lines[nextIdx]) : null;
                        if (second !== null && Number.isInteger(second) && second >= 0) {
                            // Two lone numbers: magnetization + count, or (both
                            // integers and the first no larger) count + next
                            // element's magnetization. The first reads as the
                            // magnetization in either case; the count is the
                            // integer after it.
                            count = second;
                            i = nextIdx + 1;
                        } else if (!decimal && Number.isInteger(first) && first >= 0) {
                            // A lone integer followed by a coordinate row is the
                            // count of a block that states no magnetization.
                            count = first;
                            i++;
                        } else {
                            count = first;
                            i++;
                        }
                    } else {
                        const parts = stripComment(lines[i]).split(/\s+/);
                        // Magnetization written as `mag 2.0` on its own line.
                        if (parts[0].toLowerCase() === 'mag') {
                            i++;
                            while (i < lines.length && stripComment(lines[i]) === '') i++;
                            if (i < lines.length) {
                                count = parseFloat(stripComment(lines[i]).split(/\s+/)[0]);
                                i++;
                            }
                        }
                    }
                }
                if (!isFinite(count) || count < 0) {
                    warnings.push('element ' + element + ': atom count is unreadable; the block is skipped');
                    break;
                }
                count = Math.floor(count);

                for (let k = 0; k < count && i < lines.length; k++) {
                    while (i < lines.length && stripComment(lines[i]) === '') i++;
                    if (i >= lines.length) break;
                    const row = stripComment(lines[i]).split(/\s+/);
                    const x = parseFloat(row[0]);
                    const y = parseFloat(row[1]);
                    const z = parseFloat(row[2]);
                    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) {
                        warnings.push('element ' + element + ': site ' + (k + 1) + ' is unreadable');
                        break;
                    }
                    let cart: [number, number, number];
                    let frac: [number, number, number] | undefined;
                    if (mode.kind === 'direct') {
                        frac = [x, y, z];
                        cart = [0, 0, 0];
                    } else {
                        const factor = mode.unit === 'angstrom' ? 1 : mode.unit === 'bohr' ? BOHR_TO_ANGSTROM : latticeConstant;
                        cart = [x * factor, y * factor, z * factor];
                        // A centered mode states its origin as a Direct point
                        // (0.5, 0.5, 0.0) for `_center_xy`; moving the origin
                        // back to the cell corner is one half-cell shift along
                        // the named axes.
                        if (mode.center && lattice) {
                            const half = [0, 0, 0];
                            for (let axis = 0; axis < 3; axis++) {
                                if (mode.center.indexOf('xyz'.charAt(axis)) >= 0) half[axis] = 0.5;
                            }
                            for (let axis = 0; axis < 3; axis++) {
                                cart[axis] += half[0] * lattice[0][axis] + half[1] * lattice[1][axis] + half[2] * lattice[2][axis];
                            }
                        }
                    }
                    const site: Site = { element, cart: mode.kind === 'cartesian' ? cart : undefined, frac };
                    // `m 1 0 1` after the coordinates is ABACUS's move flag
                    // (1 = free, 0 = fixed); keeping it lets a save write the
                    // fixed axes back instead of silently freeing them.
                    const moveAt = row.indexOf('m');
                    if (moveAt >= 0) {
                        const flags = row.slice(moveAt + 1, moveAt + 4);
                        if (flags.length === 3 && flags.every(f => /^[01]$/.test(f))) {
                            site.constraint = flags.map(f => f === '1' ? 'T' : 'F').join(' ');
                        }
                    }
                    sites.push(site);
                    i++;
                }
            }
            continue;
        }

        i++;
    }

    if (!lattice || !isUsableLattice(lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['no usable LATTICE_VECTORS in the STRU: no structure to show']
        };
    }
    if (sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['no ATOMIC_POSITIONS sites in the STRU: no structure to show']
        };
    }
    const data = buildPeriodicStructure(lattice, sites, title, warnings);
    data.title = oneLineTitle(title, 'ABACUS STRU');
    return data;
}

interface LogAtomType {
    element: string;
    count: number;
}

function collectAbacusTypes(lines: string[]): LogAtomType[] {
    const types: LogAtomType[] = [];
    for (let i = 0; i < lines.length; i++) {
        const m = /READING ATOM TYPE\s+(\d+)/i.exec(lines[i]);
        if (!m) continue;
        let element = '';
        let count = 0;
        for (let j = i + 1; j < Math.min(lines.length, i + 20); j++) {
            const label = /atom label\s*=\s*(\S+)/i.exec(lines[j]);
            if (label) element = normalizeElementSymbol(label[1]) || label[1];
            const number = /number of atom for this type\s*=\s*(\d+)/i.exec(lines[j]);
            if (number) count = parseInt(number[1], 10);
        }
        if (element || count > 0) types.push({ element, count });
    }
    return types;
}

function readLogLattice(lines: string[]): number[][] | null {
    // The last constant pairs with the last vectors: a variable-cell run prints
    // both again for every ionic step, and taking the first constant would
    // scale the final cell by the starting lattice parameter.
    let a0 = NaN;
    for (const line of lines) {
        const ang = /lattice constant\s*\(Angstrom\)\s*=\s*([-+0-9.eE]+)/i.exec(line);
        if (ang) a0 = parseFloat(ang[1]);
    }
    if (!isFinite(a0) || a0 <= 0) {
        for (const line of lines) {
            const bohr = /lattice constant\s*\(Bohr\)\s*=\s*([-+0-9.eE]+)/i.exec(line);
            if (bohr) a0 = parseFloat(bohr[1]) * BOHR_TO_ANGSTROM;
        }
    }
    if (!isFinite(a0) || a0 <= 0) return null;

    let header = -1;
    for (let i = lines.length - 1; i >= 0; i--) {
        if (/lattice vectors:/i.test(lines[i])) { header = i; break; }
    }
    if (header < 0) return null;
    const rows: number[][] = [];
    for (let i = header + 1; i < lines.length && rows.length < 3; i++) {
        const parts = lines[i].trim().split(/\s+/).map(Number);
        if (parts.length === 3 && parts.every(v => isFinite(v))) rows.push(parts);
    }
    if (rows.length !== 3) return null;
    return rows.map(r => r.map(v => v * a0));
}

function abacusLabelElement(label: string): string | undefined {
    const cleaned = label
        .replace(/^taud?_/i, '')
        .replace(/^tauc?_/i, '')
        .replace(/^\d+/, '')
        .replace(/\d+$/, '');
    return normalizeElementSymbol(cleaned) || elementFromAtomName(cleaned);
}

function abacusTypeElement(types: LogAtomType[], index: number): string {
    // The per-type list is a lookup, not an array the size of the log's atom
    // count: a corrupted `number of atom` must not allocate a huge table.
    let remaining = index;
    for (const t of types) {
        const n = Math.max(t.count, 1);
        // An ATOMIC_SPECIES row with no readable label leaves the type
        // unnamed; the site is then a dummy rather than an empty element.
        if (remaining < n) return t.element || 'X';
        remaining -= n;
    }
    return 'X';
}

export function parseAbacusLog(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const title = 'ABACUS SCF';
    const lattice = readLogLattice(lines);
    if (!lattice || !isUsableLattice(lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['no lattice vectors in the ABACUS log: no structure to show']
        };
    }

    const types = collectAbacusTypes(lines);

    // The last complete DIRECT COORDINATES table is the final geometry: a log
    // that ran ionic steps prints the starting cell before each step.
    const sites: Site[] = [];
    const tableStarts: number[] = [];
    for (let i = 0; i < lines.length; i++) {
        if (/DIRECT COORDINATES/i.test(lines[i])) tableStarts.push(i);
    }
    for (let t = tableStarts.length - 1; t >= 0 && sites.length === 0; t--) {
        let j = tableStarts[t] + 1;
        // Skip the column header.
        if (j < lines.length && /atom\s+x\s+y\s+z/i.test(lines[j])) j++;
        for (; j < lines.length; j++) {
            const row = lines[j].trim();
            if (row === '') {
                if (sites.length > 0) break;
                continue;
            }
            const parts = row.split(/\s+/);
            if (parts.length < 4) break;
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);
            // The label may be capitalised (`Si1`), so the table ends when the
            // coordinate columns stop being numbers, not on an uppercase letter.
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
            let element = abacusLabelElement(parts[0]);
            if (!element) element = abacusTypeElement(types, sites.length);
            sites.push({ element, frac: [x, y, z], label: parts[0] });
        }
    }

    // A log that prints only Cartesian coordinates still carries the cell.
    if (sites.length === 0) {
        let start = -1;
        for (let i = lines.length - 1; i >= 0; i--) {
            if (/CARTESIAN COORDINATES/i.test(lines[i])) { start = i; break; }
        }
        if (start >= 0) {
            // ABACUS prints `CARTESIAN COORDINATES ( UNIT = <lat0> Bohr )`:
            // the values are in units of lat0, so multiplying by the Bohr
            // constant alone would shrink the structure by lat0.
            const unit = /UNIT\s*=\s*([-+0-9.eE]+)\s*(Bohr|Angstrom|Ang|A\.U\.?)/i.exec(lines[start]);
            let factor = BOHR_TO_ANGSTROM;
            if (unit) {
                const value = parseFloat(unit[1]);
                if (isFinite(value) && value > 0) {
                    factor = /bohr|a\.u/i.test(unit[2]) ? value * BOHR_TO_ANGSTROM : value;
                }
            }
            let j = start + 1;
            if (j < lines.length && /atom\s+x\s+y\s+z/i.test(lines[j])) j++;
            for (; j < lines.length; j++) {
                const parts = lines[j].trim().split(/\s+/);
                if (parts.length < 4) break;
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
                const element = abacusLabelElement(parts[0]) || abacusTypeElement(types, sites.length);
                sites.push({ element, cart: [x * factor, y * factor, z * factor], label: parts[0] });
            }
        }
    }

    if (sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the ABACUS log states a cell but no coordinate table: no structure to show']
        };
    }
    if (tableStarts.length > 1) {
        warnings.push(tableStarts.length + ' coordinate tables in the log: showing the last');
    }
    return buildPeriodicStructure(lattice, sites, title, warnings);
}

interface MdDumpFrame {
    step: string;
    /** Angstrom, already multiplied by LATTICE_CONSTANT. */
    lattice: number[][];
    sites: Site[];
    constant?: number;
}

export function parseAbacusMdDump(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    // The viewer draws one structure, so only the first complete MDSTEP is
    // kept; a long dump is counted rather than accumulated frame by frame.
    let firstFrame: MdDumpFrame | null = null;
    let frameCount = 0;
    const isComplete = (candidate: MdDumpFrame): boolean => {
        if (candidate.sites.length === 0 || candidate.lattice.length !== 3) return false;
        frameCount++;
        return true;
    };
    const warnings: string[] = [];

    let frame: MdDumpFrame | null = null;
    let i = 0;
    while (i < lines.length) {
        const line = lines[i].trim();
        const stepMatch = /^MDSTEP:\s*(\S*)/i.exec(line);
        if (stepMatch) {
            if (frame && isComplete(frame) && firstFrame === null) firstFrame = frame;
            frame = { step: stepMatch[1] || String(frameCount), lattice: [], sites: [] };
            i++;
            continue;
        }
        if (!frame) { i++; continue; }

        if (/^LATTICE_CONSTANT:/i.test(line)) {
            const m = /^LATTICE_CONSTANT:\s*([-+0-9.eE]+)\s*(\S*)/i.exec(line);
            if (m) {
                const value = parseFloat(m[1]);
                const unit = /angstrom/i.test(m[2] || '') ? 'angstrom' : 'bohr';
                frame.constant = toAngstrom(value, unit);
            }
            i++;
            continue;
        }
        if (/^LATTICE_VECTORS/i.test(line)) {
            i++;
            const constant = frame.constant !== undefined ? frame.constant : 1;
            for (let r = 0; r < 3 && i < lines.length; r++, i++) {
                const parts = lines[i].trim().split(/\s+/).map(Number);
                if (parts.length >= 3 && parts.slice(0, 3).every(v => isFinite(v))) {
                    frame.lattice.push([parts[0] * constant, parts[1] * constant, parts[2] * constant]);
                } else {
                    break;
                }
            }
            continue;
        }
        if (/^INDEX\s+LABEL/i.test(line)) {
            // The table header states the coordinate unit; ABACUS writes
            // Angstrom today, but a Bohr dump would otherwise be off by 1.89.
            const positionUnit = /POSITION\s*\(([^)]*)\)/i.exec(line);
            const positionFactor = positionUnit && /bohr|a\.u\./i.test(positionUnit[1])
                ? BOHR_TO_ANGSTROM
                : 1;
            i++;
            while (i < lines.length) {
                const row = lines[i].trim();
                if (row === '' || /^MDSTEP:/i.test(row) || /^[A-Z_]+\s+\(/i.test(row) || /^VIRIAL/i.test(row)) break;
                const parts = row.split(/\s+/);
                if (parts.length < 5) break;
                const x = parseFloat(parts[2]);
                const y = parseFloat(parts[3]);
                const z = parseFloat(parts[4]);
                if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
                const element = normalizeElementSymbol(parts[1]) || elementFromAtomName(parts[1]) || 'X';
                frame.sites.push({ element, cart: [x * positionFactor, y * positionFactor, z * positionFactor], label: parts[1] });
                i++;
            }
            continue;
        }
        i++;
    }
    if (frame && isComplete(frame) && firstFrame === null) firstFrame = frame;

    if (firstFrame === null) {
        return {
            atoms: [], bonds: [], title: 'ABACUS MD dump', hasExplicitBonds: false,
            warnings: ['no complete MDSTEP block in the ABACUS dump: no structure to show']
        };
    }
    if (frameCount > 1) {
        warnings.push(frameCount + ' MD steps in the dump: showing the first (step ' + firstFrame.step + ')');
    }
    const first = firstFrame;
    if (!isUsableLattice(first.lattice)) {
        return {
            atoms: [], bonds: [], title: 'ABACUS MD dump', hasExplicitBonds: false,
            warnings: ['the MDSTEP lattice is unusable: no structure to show']
        };
    }
    return buildPeriodicStructure(first.lattice, first.sites, 'ABACUS MD dump (step ' + first.step + ')', warnings);
}

/** Whether a text is any ABACUS output this reader understands. */
export function isAbacusLog(content: string): boolean {
    return /READING UNITCELL INFORMATION/i.test(content) ||
        /Atomic-orbital Based Ab-initio Computation at UStc/i.test(content) ||
        /ABACUS v\d/i.test(content);
}

export function isAbacusStru(content: string): boolean {
    return /^ATOMIC_SPECIES\s*$/m.test(content) &&
        /^LATTICE_VECTORS\s*$/m.test(content) &&
        /^ATOMIC_POSITIONS\b/m.test(content);
}

export function isAbacusMdDump(content: string): boolean {
    return /^MDSTEP:/m.test(content) && /^LATTICE_VECTORS\s*$/m.test(content);
}
