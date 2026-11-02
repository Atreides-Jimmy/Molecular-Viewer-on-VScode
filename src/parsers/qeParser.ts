import { MolecularData } from '../types';
import { elementFromAtomName, normalizeElementSymbol } from '../tools/periodicTable';
import {
    BOHR_TO_ANGSTROM,
    buildPeriodicStructure,
    isUsableLattice,
    oneLineTitle,
    Site
} from './periodicUtils';

/**
 * Quantum ESPRESSO reader, for pw.x / cp.x inputs and outputs.
 *
 * The input states the cell in one of two ways: a `CELL_PARAMETERS` block, or
 * a Bravais-lattice index (`ibrav` 1..14) with `celldm`/`A,B,C`. Positions are
 * `ATOMIC_POSITIONS` in `crystal`, `alat`, `bohr` or `angstrom` units. The
 * output echoes both in a third spelling again — `crystal axes` scaled by
 * `alat` and a `tau` table in alat units — and after a relax or vc-relax it
 * prints the converged cell and coordinates in a `Begin final coordinates`
 * block. A reader that understands only one spelling silently draws a cell
 * with the wrong shape or in the wrong units, so every spelling is handled
 * here and the final block always wins over the echoed one.
 */

function parseFloatSafe(text: string): number | undefined {
    const v = parseFloat(text);
    return isFinite(v) ? v : undefined;
}

/** Split a namelist line into `key = value` assignments. */
function parseNamelist(text: string): { [key: string]: string } {
    const out: { [key: string]: string } = {};
    // A QE namelist is Fortran: assignments end at a comma or a line end, and
    // `!` starts a comment. Commas inside quoted strings (a pseudo path with a
    // comma) are rare but a quote-aware scan keeps them from splitting a value.
    let cleaned = '';
    let quote = '';
    for (let i = 0; i < text.length; i++) {
        const ch = text[i];
        if (quote !== '') {
            cleaned += ch;
            if (ch === quote) quote = '';
            continue;
        }
        if (ch === "'" || ch === '"') { quote = ch; cleaned += ch; continue; }
        if (ch === '!') { while (i < text.length && text[i] !== '\n') i++; cleaned += '\n'; continue; }
        cleaned += ch;
    }
    const parts = cleaned.split(/[,\n]/);
    for (const part of parts) {
        const eq = part.indexOf('=');
        if (eq < 0) continue;
        let key = part.slice(0, eq).trim().toLowerCase();
        if (key === '' || key.startsWith('&')) continue;
        key = key.replace(/\s+/g, '');
        const value = part.slice(eq + 1).trim().replace(/^['"]|['"]$/g, '').trim();
        if (value !== '') out[key] = value;
    }
    return out;
}

function namelistBlock(lines: string[], name: string): string {
    const start = lines.findIndex(l => new RegExp('^\\s*&' + name + '\\b', 'i').test(l));
    if (start < 0) return '';
    const chunk: string[] = [];
    for (let i = start + 1; i < lines.length; i++) {
        // Any new &NAME line ends the block, including one on the very next
        // line: a namelist with no body is empty, not a continuation of this one.
        if (/^\s*&\w+/.test(lines[i])) break;
        if (/^\s*\/\s*$/.test(lines[i])) break;
        chunk.push(lines[i]);
    }
    return chunk.join('\n');
}

function numValue(map: { [key: string]: string }, key: string): number | undefined {
    const v = map[key];
    if (v === undefined) return undefined;
    return parseFloatSafe(v);
}

function celldmValue(map: { [key: string]: string }, index: number): number | undefined {
    return numValue(map, 'celldm(' + index + ')');
}

interface QeSystem {
    ibrav: number;
    /** Lengths in Angstrom. */
    a: number;
    b: number;
    c: number;
    cosAB: number;
    cosAC: number;
    cosBC: number;
    /** Cartesian positions in `alat` units are scaled by this (Angstrom). */
    alat: number;
}

function readQeSystem(lines: string[]): QeSystem {
    const map = parseNamelist(namelistBlock(lines, 'SYSTEM'));
    const ibrav = (() => {
        const v = parseInt(map['ibrav'] || '', 10);
        return isFinite(v) ? v : 0;
    })();
    const celldm1 = celldmValue(map, 1);
    const celldm2 = celldmValue(map, 2);
    const celldm3 = celldmValue(map, 3);
    const celldm4 = celldmValue(map, 4);
    const celldm5 = celldmValue(map, 5);
    const celldm6 = celldmValue(map, 6);

    let a = numValue(map, 'a');
    let b = numValue(map, 'b');
    let c = numValue(map, 'c');
    // A/B/C are already Angstrom; celldm is in Bohr, and celldm(2)/(3) are
    // ratios to celldm(1), not lengths in their own right.
    const aFromCelldm = celldm1 !== undefined && celldm1 > 0 ? celldm1 * BOHR_TO_ANGSTROM : undefined;
    if (a === undefined) a = aFromCelldm;
    if (a === undefined) a = 1;
    if (b === undefined) b = a * (celldm2 !== undefined && celldm2 > 0 ? celldm2 : 1);
    if (c === undefined) c = a * (celldm3 !== undefined && celldm3 > 0 ? celldm3 : 1);

    let cosAB = numValue(map, 'cosab');
    let cosAC = numValue(map, 'cosac');
    let cosBC = numValue(map, 'cosbc');
    // QE's angle parameters are per-Bravais, not one fixed meaning: celldm(4)
    // is cos(ab) for the monoclinic/trigonal lattices but cos(bc) for the
    // triclinic one, while celldm(5) is cos(ac) and celldm(6) cos(ab). The
    // mapping below is the one latgen.f90 uses (celldm2abc).
    if (cosAB === undefined) {
        if (ibrav === 14) cosAB = celldm6;
        else if (ibrav === -12 || ibrav === -13) cosAB = 0;
        else cosAB = celldm4;
    }
    if (cosAC === undefined) {
        cosAC = (ibrav === 14 || ibrav === -12 || ibrav === -13) ? celldm5 : 0;
    }
    if (cosBC === undefined) cosBC = ibrav === 14 ? celldm4 : 0;

    return {
        ibrav,
        a, b, c,
        cosAB: cosAB !== undefined ? cosAB : 0,
        cosAC: cosAC !== undefined ? cosAC : 0,
        cosBC: cosBC !== undefined ? cosBC : 0,
        alat: aFromCelldm !== undefined ? aFromCelldm : a
    };
}

/**
 * The primitive vectors of an `ibrav` cell, in Angstrom.
 *
 * The fourteen Bravais lattices are QE's own numbering, not the crystallographic
 * one, so each is written out from the QE manual rather than derived from
 * `cellToLatticeVectors`: a rhombohedral ibrav 5 has its first vector along x,
 * while the generic cell builder would put b in the xy plane and c somewhere
 * else, and a save-back to `CELL_PARAMETERS` must reproduce the cell the file
 * meant.
 */
function latticeFromIbrav(sys: QeSystem): number[][] | null {
    const { a, b, c } = sys;
    const cosAB = sys.cosAB;
    const cosAC = sys.cosAC;
    const cosBC = sys.cosBC;
    // The formulas are the ones pw.x builds in latgen.f90, including the
    // negative indices that select the alternative axis of a monoclinic or
    // rhombohedral cell; a different (even isometric) primitive cell would make
    // `crystal` positions mean another site.
    switch (sys.ibrav) {
        case 1:
            return [[a, 0, 0], [0, a, 0], [0, 0, a]];
        case 2: {
            const t = a / 2;
            return [[-t, 0, t], [0, t, t], [-t, t, 0]];
        }
        case 3: {
            const t = a / 2;
            return [[t, t, t], [-t, t, t], [-t, -t, t]];
        }
        case -3: {
            const t = a / 2;
            return [[-t, t, t], [t, -t, t], [t, t, -t]];
        }
        case 4:
            return [[a, 0, 0], [-a / 2, a * Math.sqrt(3) / 2, 0], [0, 0, c]];
        case 5:
        case -5: {
            const term1 = Math.sqrt(Math.max(0, 1 + 2 * cosAB));
            const term2 = Math.sqrt(Math.max(0, 1 - cosAB));
            if (sys.ibrav === 5) {
                // threefold axis along c
                const x = a * term2 / Math.SQRT2;
                const y = -x / Math.sqrt(3);
                const az = a * term1 / Math.sqrt(3);
                const ay = Math.SQRT2 * a * term2 / Math.sqrt(3);
                return [[x, y, az], [0, ay, az], [-x, y, az]];
            }
            // threefold axis along (111)
            const x = a * (term1 - 2 * term2) / 3;
            const y = a * (term1 + term2) / 3;
            return [[x, y, y], [y, x, y], [y, y, x]];
        }
        case 6:
            return [[a, 0, 0], [0, a, 0], [0, 0, c]];
        case 7: {
            const t = a / 2;
            return [[t, -t, c / 2], [t, t, c / 2], [-t, -t, c / 2]];
        }
        case 8:
            return [[a, 0, 0], [0, b, 0], [0, 0, c]];
        case 9:
            return [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]];
        case -9:
            return [[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]];
        case 91:
            return [[a, 0, 0], [0, b / 2, -c / 2], [0, b / 2, c / 2]];
        case 10:
            return [[a / 2, 0, c / 2], [a / 2, b / 2, 0], [0, b / 2, c / 2]];
        case 11:
            return [[a / 2, b / 2, c / 2], [-a / 2, b / 2, c / 2], [-a / 2, -b / 2, c / 2]];
        case 12: {
            const sb = Math.sqrt(Math.max(0, 1 - cosAB * cosAB));
            return [[a, 0, 0], [b * cosAB, b * sb, 0], [0, 0, c]];
        }
        case -12: {
            const sb = Math.sqrt(Math.max(0, 1 - cosAC * cosAC));
            return [[a, 0, 0], [0, b, 0], [c * cosAC, 0, c * sb]];
        }
        case 13: {
            const sb = Math.sqrt(Math.max(0, 1 - cosAB * cosAB));
            return [[a / 2, 0, -c / 2], [b * cosAB, b * sb, 0], [a / 2, 0, c / 2]];
        }
        case -13: {
            const sb = Math.sqrt(Math.max(0, 1 - cosAC * cosAC));
            return [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [c * cosAC, 0, c * sb]];
        }
        case 14: {
            const sb = Math.sqrt(Math.max(0, 1 - cosAB * cosAB)) || 1;
            const cy = c * (cosBC - cosAB * cosAC) / sb;
            const cz2 = 1 - cosAB * cosAB - cosAC * cosAC - cosBC * cosBC + 2 * cosAB * cosAC * cosBC;
            const cz = c * Math.sqrt(Math.max(0, cz2)) / sb;
            return [[a, 0, 0], [b * cosAB, b * sb, 0], [c * cosAC, cy, cz]];
        }
        default:
            return null;
    }
}

function alignmentUnit(header: string): string {
    const m = /(?:CELL_PARAMETERS|ATOMIC_POSITIONS)\s*(?:\(|\{)?\s*(alat|bohr|angstrom|ang|crystal|au)/i.exec(header);
    if (!m) return '';
    const u = m[1].toLowerCase();
    return u === 'ang' ? 'angstrom' : u === 'au' ? 'bohr' : u;
}

/**
 * The cards that can follow a coordinate block. A species label (`Fe1`, `Si`)
 * is not one of them, which is why the test is a fixed word list rather than
 * `starts with a capital letter`.
 */
function isQeCardLine(line: string): boolean {
    return /^\s*(?:&[A-Za-z]|CELL_PARAMETERS|ATOMIC_POSITIONS|ATOMIC_SPECIES|K_POINTS|OCCUPATIONS|CONSTRAINTS|ATOMIC_VELOCITIES|ATOMIC_FORCES|HUBBARD|ADDITIONAL_K_POINTS|SOLVENTS|END)\b/i.test(line);
}

function readThreeByThree(lines: string[], start: number): { rows: number[][]; next: number } {
    const rows: number[][] = [];
    let i = start;
    while (i < lines.length && rows.length < 3) {
        const line = lines[i].trim();
        if (line === '') { i++; continue; }
        const parts = line.split(/\s+/);
        if (parts.length < 3) break;
        const v = [parseFloat(parts[0]), parseFloat(parts[1]), parseFloat(parts[2])];
        if (v.some(x => !isFinite(x))) break;
        rows.push(v);
        i++;
    }
    return { rows, next: i };
}

/** A CELL_PARAMETERS block, converted to Angstrom. */
function readCellParameters(lines: string[], start: number, alat: number): { lattice: number[][] | null; next: number } {
    const header = lines[start];
    const unit = alignmentUnit(header);
    const read = readThreeByThree(lines, start + 1);
    if (read.rows.length !== 3) return { lattice: null, next: start + 1 };
    // QE's default for CELL_PARAMETERS is alat, the same as for positions;
    // reading an unqualified block as Angstrom would rescale a cell written in
    // units of the lattice parameter by the lattice parameter itself.
    const factor = unit === 'bohr' ? BOHR_TO_ANGSTROM : unit === 'angstrom' ? 1 : alat;
    return {
        lattice: read.rows.map(r => r.map(v => v * factor)),
        next: read.next
    };
}

function readAtomicSpecies(lines: string[]): { [label: string]: string } {
    const map: { [label: string]: string } = {};
    for (let i = 0; i < lines.length; i++) {
        if (!/^\s*ATOMIC_SPECIES\b/i.test(lines[i])) continue;
        for (let j = i + 1; j < lines.length; j++) {
            const line = lines[j].trim();
            if (line === '') continue;
            if (isQeCardLine(line)) break;
            const parts = line.split(/\s+/);
            if (parts.length < 2) break;
            const label = parts[0];
            const element = normalizeElementSymbol(label) || elementFromAtomName(label);
            if (element) map[label] = element;
        }
        break;
    }
    return map;
}

function readAtomicPositions(
    lines: string[],
    start: number,
    alat: number,
    species: { [label: string]: string },
    nat: number
): { rows: Site[]; next: number } {
    const header = lines[start];
    const explicit = alignmentUnit(header);
    // QE's default for ATOMIC_POSITIONS is alat; crystal is the only unit that
    // means fractional coordinates, and `alat` is Cartesian in units of a.
    const unit = explicit || 'alat';
    const rows: Site[] = [];
    let i = start + 1;
    while (i < lines.length && (nat <= 0 || rows.length < nat)) {
        const line = lines[i].trim();
        if (line === '') { i++; continue; }
        if (isQeCardLine(line)) break;
        const parts = line.split(/\s+/);
        if (parts.length < 4) break;
        const x = parseFloat(parts[1]);
        const y = parseFloat(parts[2]);
        const z = parseFloat(parts[3]);
        if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
        const label = parts[0];
        const element = species[label] || normalizeElementSymbol(label) || elementFromAtomName(label) || label;
        const constraint = (() => {
            if (parts.length < 7) return undefined;
            const flags = parts.slice(4, 7);
            if (!flags.every(f => /^[01]$/.test(f))) return undefined;
            return flags.map(f => f === '1' ? 'T' : 'F').join(' ');
        })();
        let site: Site;
        if (unit === 'crystal') {
            site = { element, frac: [x, y, z], label };
        } else {
            const factor = unit === 'bohr' ? BOHR_TO_ANGSTROM : unit === 'alat' ? alat : 1;
            site = { element, cart: [x * factor, y * factor, z * factor], label };
        }
        if (constraint) site.constraint = constraint;
        rows.push(site);
        i++;
    }
    return { rows, next: i };
}

function natFromSystem(lines: string[]): number {
    const map = parseNamelist(namelistBlock(lines, 'SYSTEM'));
    const raw = parseInt(map['nat'] || '', 10);
    return isFinite(raw) && raw > 0 ? raw : 0;
}

export function parseQeInput(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const sys = readQeSystem(lines);
    const species = readAtomicSpecies(lines);
    const nat = natFromSystem(lines);
    const title = (() => {
        const control = parseNamelist(namelistBlock(lines, 'CONTROL'));
        const prefix = control['prefix'];
        return oneLineTitle(prefix ? 'Quantum ESPRESSO input: ' + prefix : 'Quantum ESPRESSO input', 'Quantum ESPRESSO input');
    })();

    let lattice: number[][] | null = sys.ibrav === 0 ? null : latticeFromIbrav(sys);

    // A CELL_PARAMETERS block overrides the Bravais cell (ibrav 0 requires it).
    for (let i = 0; i < lines.length; i++) {
        if (!/^\s*CELL_PARAMETERS\b/i.test(lines[i])) continue;
        const read = readCellParameters(lines, i, sys.alat);
        if (read.lattice) lattice = read.lattice;
        else warnings.push('CELL_PARAMETERS is incomplete: the cell could not be read');
    }
    if (!lattice) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the Quantum ESPRESSO input states no usable cell: no structure to show']
        };
    }

    let sites: Site[] = [];
    let sawPositions = false;
    for (let i = 0; i < lines.length; i++) {
        if (!/^\s*ATOMIC_POSITIONS\b/i.test(lines[i])) continue;
        sawPositions = true;
        // `crystal_sg` lists the symmetry-unique sites; the rest are generated
        // from `space_group`. The reader does not expand them, so it says so
        // rather than drawing a cell with only the listed atoms.
        if (/crystal_sg/i.test(lines[i])) {
            warnings.push('ATOMIC_POSITIONS crystal_sg lists symmetry-unique sites: only the listed atoms are shown');
        }
        const read = readAtomicPositions(lines, i, sys.alat, species, nat);
        if (read.rows.length > 0) sites = read.rows;
    }
    if (!sawPositions || sites.length === 0) {
        warnings.push('the Quantum ESPRESSO input states no ATOMIC_POSITIONS: no atoms to show');
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings
        };
    }
    if (!isUsableLattice(lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the Quantum ESPRESSO input cell is degenerate: no structure to show']
        };
    }
    return buildPeriodicStructure(lattice, sites, title, warnings);
}

interface OutputGeometry {
    lattice: number[][] | null;
    sites: Site[];
    warnings: string[];
}

/** The crystal axes echoed by pw.x, in Angstrom. */
function readOutputAxes(lines: string[]): { lattice: number[][] | null; alat: number } {
    let alat = NaN;
    for (const line of lines) {
        const m = /lattice parameter\s*\(alat\)\s*=\s*([-+0-9.eE]+)\s*a\.u\./i.exec(line);
        if (m) alat = parseFloat(m[1]) * BOHR_TO_ANGSTROM;
    }
    if (!isFinite(alat) || alat <= 0) return { lattice: null, alat: 1 };

    let header = -1;
    for (let i = lines.length - 1; i >= 0; i--) {
        if (/crystal axes:/i.test(lines[i])) { header = i; break; }
    }
    if (header < 0) return { lattice: null, alat };
    const rows: number[][] = [];
    for (let i = header + 1; i < lines.length && rows.length < 3; i++) {
        const m = /a\(\d\)\s*=\s*\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)/i.exec(lines[i]);
        if (m) rows.push([parseFloat(m[1]) * alat, parseFloat(m[2]) * alat, parseFloat(m[3]) * alat]);
    }
    return { lattice: rows.length === 3 ? rows : null, alat };
}

/** The `site n. atom tau(i) = (...)` table, in alat units. */
function readTauTable(lines: string[], alat: number): Site[] {
    const sites: Site[] = [];
    for (let i = lines.length - 1; i >= 0 && sites.length === 0; i--) {
        if (!/positions\s*\(alat units\)/i.test(lines[i])) continue;
        for (let j = i + 1; j < lines.length; j++) {
            const m = /^\s*(\d+)\s+(\S+)\s+tau\(\s*\d+\s*\)\s*=\s*\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)/i.exec(lines[j]);
            if (!m) {
                if (sites.length > 0 && lines[j].trim() === '') break;
                if (sites.length > 0) break;
                continue;
            }
            const element = normalizeElementSymbol(m[2]) || elementFromAtomName(m[2]) || m[2];
            sites.push({
                element,
                cart: [parseFloat(m[3]) * alat, parseFloat(m[4]) * alat, parseFloat(m[5]) * alat],
                label: m[2]
            });
        }
    }
    return sites;
}

function parseOutputGeometry(content: string): OutputGeometry {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const axes = readOutputAxes(lines);
    let lattice = axes.lattice;

    // The converged cell and coordinates, when the run relaxed. Both may be
    // present, or only the positions (a fixed-cell relax).
    let beginFinal = -1;
    let endFinal = lines.length;
    for (let i = 0; i < lines.length; i++) {
        if (/Begin final coordinates/i.test(lines[i])) beginFinal = i;
        if (beginFinal >= 0 && /End final coordinates/i.test(lines[i])) { endFinal = i; break; }
    }
    let finalSites: Site[] = [];
    if (beginFinal >= 0) {
        const species: { [label: string]: string } = {};
        const nat = 0;
        for (let i = beginFinal; i < endFinal; i++) {
            if (/^\s*CELL_PARAMETERS\b/i.test(lines[i])) {
                const read = readCellParameters(lines, i, axes.alat);
                if (read.lattice) lattice = read.lattice;
            } else if (/^\s*ATOMIC_POSITIONS\b/i.test(lines[i])) {
                const read = readAtomicPositions(lines, i, axes.alat, species, nat);
                if (read.rows.length > 0) finalSites = read.rows;
            }
        }
    }

    let sites = finalSites;
    if (sites.length === 0) {
        // A single-point run (and any run without a final block) draws the
        // geometry pw.x echoed in its header.
        sites = readTauTable(lines, axes.alat > 0 ? axes.alat : 1);
        if (sites.length === 0) {
            // cp.x does not print a tau table; it does print ATOMIC_POSITIONS
            // blocks in its own output, so the last one is the best fallback.
            for (let i = lines.length - 1; i >= 0 && sites.length === 0; i--) {
                if (!/^\s*ATOMIC_POSITIONS\b/i.test(lines[i])) continue;
                const read = readAtomicPositions(lines, i, axes.alat, {}, 0);
                sites = read.rows;
            }
        }
        if (sites.length > 0 && finalSites.length === 0 && beginFinal >= 0) {
            warnings.push('the final coordinate block cannot be read: showing the geometry the output echoed');
        }
    }

    if (!lattice) warnings.push('the Quantum ESPRESSO output states no crystal axes: no cell to show');
    if (sites.length === 0) warnings.push('the Quantum ESPRESSO output states no coordinate table: no atoms to show');
    return { lattice, sites, warnings };
}

export function parseQeOutput(content: string): MolecularData {
    const geometry = parseOutputGeometry(content);
    const titleMatch = /^\s*Program\s+(PWSCF|CP)\s+v\.?\s*([0-9.]+)/im.exec(content);
    const title = titleMatch ? 'Quantum ESPRESSO ' + titleMatch[1].toUpperCase() + ' ' + titleMatch[2] : 'Quantum ESPRESSO output';
    if (!geometry.lattice || geometry.sites.length === 0 || !isUsableLattice(geometry.lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry.warnings.length > 0 ? geometry.warnings : ['no structure to show']
        };
    }
    return buildPeriodicStructure(geometry.lattice, geometry.sites, title, geometry.warnings);
}

export function isQeInput(content: string): boolean {
    return /^\s*&SYSTEM\b/im.test(content) &&
        (/^\s*ATOMIC_POSITIONS\b/im.test(content) || /^\s*K_POINTS\b/im.test(content) || /^\s*CELL_PARAMETERS\b/im.test(content));
}

export function isQeOutput(content: string): boolean {
    if (/^\s*Program\s+(PWSCF|CP)\b/im.test(content)) return true;
    return /bravais-lattice index\s*=/i.test(content) &&
        /lattice parameter\s*\(alat\)/i.test(content) &&
        /crystal axes:/i.test(content);
}
