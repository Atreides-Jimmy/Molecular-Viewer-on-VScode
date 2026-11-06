import { MolecularData } from '../types';
import { elementFromZ, normalizeElementSymbol } from '../tools/periodicTable';
import { BOHR_TO_ANGSTROM, buildMolecularStructure, buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * PySCF output and Molden readers.
 *
 * PySCF prints a geometry three ways, and only one of them is the geometry the
 * run finished with:
 *
 *  - the `[INPUT]` table, once, at the start (both Angstrom and Bohr columns);
 *  - `Cartesian coordinates (Angstrom)` inside every optimization cycle, with
 *    a dX dY dZ column that is easy to mistake for a position;
 *  - `New geometry`, after each cycle, again with both unit columns.
 *
 * The last `New geometry` block wins, the dX/dY/dZ table never does, and the
 * `[INPUT]` table is the fallback for a single-point run. A periodic Cell also
 * logs `lattice vectors a1/a2/a3 [...]`; PySCF prints those in **Bohr** (its
 * internal unit) even when the Cell was built in Angstrom, so they are
 * converted here. The text carries no unit marker, which is why the conversion
 * is not conditional.
 *
 * Molden files state their atom unit in the `[Atoms]` header (`AU` or `Angs`);
 * the rows after it are `element index atomicNumber x y z`.
 */

interface PyscfGeometry {
    lattice: number[][] | null;
    sites: Site[];
    warnings: string[];
}

const GEOMETRY_ROW = /^\s*(?:\[INPUT\]\s*)?(\d+)\s+([A-Za-z]{1,3})\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+(AA|Angstrom|Ang|AU|Bohr)\b/i;

function unitFactor(unit: string): number {
    return /^(AA|Angstrom|Ang)$/i.test(unit) ? 1 : BOHR_TO_ANGSTROM;
}

function geometryRow(line: string): Site | null {
    const m = GEOMETRY_ROW.exec(line);
    if (!m) return null;
    const x = parseFloat(m[3]);
    const y = parseFloat(m[4]);
    const z = parseFloat(m[5]);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return null;
    const f = unitFactor(m[6]);
    return {
        element: normalizeElementSymbol(m[2]) || m[2],
        cart: [x * f, y * f, z * f],
        label: m[2]
    };
}

function readNewGeometryBlock(lines: string[], start: number): Site[] {
    const rows: Site[] = [];
    for (let i = start + 1; i < lines.length; i++) {
        if (lines[i].trim() === '') {
            if (rows.length > 0) break;
            continue;
        }
        const site = geometryRow(lines[i]);
        if (!site) break;
        rows.push(site);
    }
    return rows;
}

function readLattice(lines: string[]): number[][] | null {
    const vectors: number[][] = [];
    for (let i = 0; i < lines.length; i++) {
        const m = /lattice vectors\s+a1\s*\[([^\]]+)\]/i.exec(lines[i]);
        if (!m) continue;
        const first = m[1].split(',').map(Number);
        const second = /a2\s*\[([^\]]+)\]/i.exec(lines[i + 1] || '');
        const third = /a3\s*\[([^\]]+)\]/i.exec(lines[i + 2] || '');
        if (first.length < 3 || !second || !third) continue;
        const b = second[1].split(',').map(Number);
        const c = third[1].split(',').map(Number);
        if (b.length < 3 || c.length < 3) continue;
        // The log prints the internal (Bohr) cell even for an Angstrom Cell.
        vectors.push([
            first[0] * BOHR_TO_ANGSTROM, first[1] * BOHR_TO_ANGSTROM, first[2] * BOHR_TO_ANGSTROM
        ], [
            b[0] * BOHR_TO_ANGSTROM, b[1] * BOHR_TO_ANGSTROM, b[2] * BOHR_TO_ANGSTROM
        ], [
            c[0] * BOHR_TO_ANGSTROM, c[1] * BOHR_TO_ANGSTROM, c[2] * BOHR_TO_ANGSTROM
        ]);
    }
    if (vectors.length >= 3) {
        const lattice = vectors.slice(vectors.length - 3);
        return isUsableLattice(lattice) ? lattice : null;
    }
    return null;
}

export function isPyscfOutput(content: string): boolean {
    return /PySCF version/i.test(content) ||
        /^\*\*\*\*\*\*\*\* <class 'pyscf\./m.test(content) ||
        /^\s*\[INPUT\]\s+num\. atoms\s*=/m.test(content);
}

export function isMolden(content: string): boolean {
    return /^\s*\[Molden Format\]/im.test(content) && /^\s*\[Atoms\]/im.test(content);
}

export function parseMolden(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    let title = 'Molden file';
    if (lines.length > 1 && lines[0].trim() !== '' && !/^\s*\[/.test(lines[1])) {
        title = oneLineTitle(lines[1], title);
    }

    let start = -1;
    let bohr = true;
    for (let i = 0; i < lines.length; i++) {
        // The unit follows the section name either in parentheses (`(AU)`,
        // `(Angs)`) or bare (`[Atoms] Angs`); both spellings are common, and a
        // header without a unit means Bohr, which is Molden's default.
        const m = /^\s*\[Atoms\]\s*\(?\s*([A-Za-z.]*)\s*\)?/i.exec(lines[i]);
        if (!m) continue;
        start = i + 1;
        bohr = !/ang/i.test(m[1] || 'AU');
        break;
    }
    if (start < 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the Molden file has no [Atoms] section: no structure to show']
        };
    }

    const sites: Site[] = [];
    for (let i = start; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line === '') continue;
        if (line.startsWith('[')) break;
        const parts = line.split(/\s+/);
        if (parts.length < 6) break;
        const z = parseInt(parts[2], 10);
        const x = parseFloat(parts[3]);
        const y = parseFloat(parts[4]);
        const zc = parseFloat(parts[5]);
        if (!isFinite(x) || !isFinite(y) || !isFinite(zc)) break;
        const element = normalizeElementSymbol(parts[0]) || elementFromZ(z) || parts[0];
        const f = bohr ? BOHR_TO_ANGSTROM : 1;
        sites.push({ element, cart: [x * f, y * f, zc * f], label: parts[0] });
    }
    if (sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the Molden [Atoms] section has no readable rows: no structure to show']
        };
    }
    return buildMolecularStructure(sites, title, warnings);
}

function readPyscfGeometry(content: string): PyscfGeometry {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const lattice = readLattice(lines);

    // The run states its atom count once, in the opening [INPUT] block. A
    // geometry block with fewer rows than that is a table cut short (the job
    // was killed, the log was truncated); letting it replace the last complete
    // geometry would show a fragment as the finished structure.
    let expected = 0;
    for (const line of lines) {
        const m = /^\s*\[INPUT\]\s*num\.\s*atoms\s*=\s*(\d+)/i.exec(line);
        if (m) { expected = parseInt(m[1], 10); break; }
    }

    // The last `New geometry` is the one the optimizer finished with.
    let sites: Site[] = [];
    let partial = false;
    for (let i = 0; i < lines.length; i++) {
        if (lines[i].trim() !== 'New geometry') continue;
        const rows = readNewGeometryBlock(lines, i);
        if (rows.length === 0) continue;
        if (expected > 0 && rows.length < expected) { partial = true; continue; }
        sites = rows;
    }
    if (sites.length === 0) {
        // A single-point run: the starting table is the geometry.
        for (const line of lines) {
            if (line.indexOf('[INPUT]') < 0) continue;
            const site = geometryRow(line);
            if (site) sites.push(site);
            if (expected > 0 && sites.length >= expected) break;
        }
    }
    if (partial) {
        warnings.push('the last New geometry table is shorter than the stated atom count: ' +
            'the previous complete geometry is shown');
    }
    return { lattice, sites, warnings };
}

export function parsePyscfOutput(content: string): MolecularData {
    const geometry = readPyscfGeometry(content);
    const title = 'PySCF output';
    if (geometry.sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the PySCF output states no atom table: no structure to show']
        };
    }
    if (geometry.lattice) {
        return buildPeriodicStructure(geometry.lattice, geometry.sites, title, geometry.warnings);
    }
    return buildMolecularStructure(geometry.sites, title, geometry.warnings);
}
