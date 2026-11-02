import { MolecularData } from '../types';
import { elementFromMass, normalizeElementSymbol, typeElementLabel } from '../tools/periodicTable';
import { fracToCart } from '../tools/crystalGeometry';
import { buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * LAMMPS reader: the data file and the dump trajectory.
 *
 * LAMMPS names atoms by *type*, not by element, so the element a type stands for
 * lives in another file: a `Masses` section (the mass nearly always names the
 * element), a `mass` command, or a `dump_modify ... element` list in the input
 * script. The loader reads those companions and passes them in; this reader
 * keeps the honest `T3` label when nothing names the type instead of guessing.
 *
 * The two file families also state their box differently: a data file lists the
 * true triclinic bounds, while a dump lists the *bounding* box of a tilted cell
 * (`xlo_bound`, `xhi_bound`, `xy`). Applying the data-file reading to a dump
 * shifts a triclinic cell sideways; the correction below is LAMMPS's own
 * definition of the two.
 */

export interface LammpsParseOptions {
    /** Type index (1-based) -> element, from a sibling data file or input script. */
    typeElements?: { [type: number]: string };
    /** The run's `units` value, when a companion input script states it. */
    units?: string;
    /** The `atom_style` a companion input script states for this data file. */
    atomStyle?: string;
}

const UNIT_TO_ANGSTROM: { [unit: string]: number } = {
    metal: 1,
    real: 1,
    nano: 10,
    micro: 10000,
    si: 1e10,
    cgs: 1e8,
    electron: 0.529177210903,
    lj: 1
};

function lengthFactor(units: string | undefined): number {
    if (!units) return 1;
    const f = UNIT_TO_ANGSTROM[units.toLowerCase()];
    return f !== undefined ? f : 1;
}

function isNumeric(token: string | undefined): boolean {
    return token !== undefined && token !== '' && isFinite(parseFloat(token));
}

const SECTION_START = /^(Masses|Atoms|Velocities|Bonds|Angles|Dihedrals|Impropers|Pair Coeffs|PairIJ Coeffs|Bond Coeffs|Angle Coeffs|Dihedral Coeffs|Improper Coeffs|BondBond Coeffs|BondAngle Coeffs|MiddleBondTorsion Coeffs|EndBondTorsion Coeffs|AngleTorsion Coeffs|AngleAngleTorsion Coeffs|BondBond13 Coeffs|AngleAngle Coeffs|Ellipsoids|Lines|Triangles|Bodies|Molecules)\b/i;

interface DataBox {
    lattice: number[][] | null;
    /** True box bounds in file units, before conversion. */
    xlo: number;
    xhi: number;
    ylo: number;
    yhi: number;
    zlo: number;
    zhi: number;
    tilt: [number, number, number];
}

function readDataBox(lines: string[]): DataBox {
    const bound: { [key: string]: number } = {};
    let tilt: [number, number, number] = [0, 0, 0];
    for (const line of lines) {
        const m = /^\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+(xlo|ylo|zlo)\s+(xhi|yhi|zhi)\b/i.exec(line);
        if (m) {
            const lo = m[3].toLowerCase();
            const hi = m[4].toLowerCase();
            bound[lo] = parseFloat(m[1]);
            bound[hi] = parseFloat(m[2]);
            continue;
        }
        const t = /^\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+xy\s+xz\s+yz\b/i.exec(line);
        if (t) tilt = [parseFloat(t[1]), parseFloat(t[2]), parseFloat(t[3])];
    }
    const xlo = bound['xlo'] !== undefined ? bound['xlo'] : 0;
    const xhi = bound['xhi'] !== undefined ? bound['xhi'] : 0;
    const ylo = bound['ylo'] !== undefined ? bound['ylo'] : 0;
    const yhi = bound['yhi'] !== undefined ? bound['yhi'] : 0;
    const zlo = bound['zlo'] !== undefined ? bound['zlo'] : 0;
    const zhi = bound['zhi'] !== undefined ? bound['zhi'] : 0;
    const lx = xhi - xlo;
    const ly = yhi - ylo;
    const lz = zhi - zlo;
    let lattice: number[][] | null = null;
    if (lx > 0 && ly > 0 && lz > 0) {
        lattice = [[lx, 0, 0], [tilt[0], ly, 0], [tilt[1], tilt[2], lz]];
    }
    return { lattice, xlo, xhi, ylo, yhi, zlo, zhi, tilt };
}

function parseMasses(lines: string[]): { [type: number]: string } {
    const out: { [type: number]: string } = {};
    for (let i = 0; i < lines.length; i++) {
        if (!/^Masses\b/i.test(lines[i].trim())) continue;
        for (let j = i + 1; j < lines.length; j++) {
            const line = lines[j].trim();
            if (line === '') {
                // A blank line ends the section unless more rows follow it in
                // the same section; the next section header is the real end.
                let k = j + 1;
                while (k < lines.length && lines[k].trim() === '') k++;
                if (k >= lines.length || SECTION_START.test(lines[k].trim())) break;
                j = k - 1;
                continue;
            }
            if (SECTION_START.test(line)) break;
            const parts = line.split(/\s+/);
            if (parts.length < 2) continue;
            const type = parseInt(parts[0], 10);
            const mass = parseFloat(parts[1]);
            if (!isFinite(type) || type <= 0) continue;
            const comment = /#\s*([A-Za-z][A-Za-z0-9_+-]*)/.exec(line);
            const tag = comment ? comment[1] : '';
            // A zero-mass type is a dummy point; the spellings the other
            // formats use for a ghost name it explicitly in the comment.
            const dummy = /^(x|xx|du|dummy|bq)$/i.test(tag);
            const fromComment = dummy ? 'X' : (tag ? normalizeElementSymbol(tag) : undefined);
            const element = fromComment ||
                (isFinite(mass) && mass === 0 ? 'X' : (isFinite(mass) ? elementFromMass(mass) : undefined));
            if (element) out[type] = element;
        }
        break;
    }
    return out;
}

export function lammpsTypeElements(content: string): { [type: number]: string } {
    return parseMasses(content.split(/\r?\n/));
}

interface AtomStyle {
    kind: 'atomic' | 'charge' | 'full' | 'molecular';
    /** Index of the x column in a split row. */
    x: number;
    /** Index of the type column (always 1 in the styles below). */
    type: number;
}

/** The layouts this reader knows by name; any other style is read by width. */
function atomStyleFromName(name: string): AtomStyle | null {
    if (name === 'atomic') return { kind: 'atomic', x: 2, type: 1 };
    if (name === 'charge') return { kind: 'charge', x: 3, type: 1 };
    if (name === 'full') return { kind: 'full', x: 4, type: 2 };
    // `id mol type x y z`: unlike `full`, molecular has no charge column, so
    // x is the fourth field, not the fifth.
    if (name === 'molecular' || name === 'bond' || name === 'angle') return { kind: 'molecular', x: 3, type: 2 };
    return null;
}

function isIntegerToken(token: string | undefined): boolean {
    return token !== undefined && /^[+-]?\d+$/.test(token);
}

/**
 * Widths 6 and 9 hold either `id type q x y z [flags]` (charge) or
 * `id mol type x y z [flags]` (molecular). The coordinates start at column 3
 * in both, so only the type column differs. A fractional or non-positive
 * third column cannot be an atom type, and a value above the declared type
 * count says which column is the molecule index. When neither test decides,
 * charge wins: a molecular file states its style in the Atoms header or the
 * input script in practice, while integral ionic charges are common.
 */
function chargeOrMolecular(rows: string[][], atomTypes: number | undefined): AtomStyle {
    let molecular = 0;
    let charge = 0;
    for (const row of rows) {
        const third = row[2];
        if (!isIntegerToken(third)) { charge += 2; continue; }
        const thirdValue = parseInt(third || '', 10);
        if (thirdValue <= 0) { charge += 2; continue; }
        if (atomTypes !== undefined && atomTypes > 0) {
            if (parseInt(row[1] || '', 10) > atomTypes) molecular++;
            if (thirdValue > atomTypes) charge++;
        }
    }
    return molecular > charge
        ? { kind: 'molecular', x: 3, type: 2 }
        : { kind: 'charge', x: 3, type: 1 };
}

function atomStyleOf(header: string, rows: string[][], hint?: string, atomTypes?: number): AtomStyle {
    const comment = /#\s*(\S+)/.exec(header);
    // The file's own comment describes the file being read, so it wins over a
    // hint from the input script.
    const named = atomStyleFromName(comment ? comment[1].toLowerCase() : '');
    if (named) return named;
    const hinted = hint ? atomStyleFromName(hint.toLowerCase()) : null;
    if (hinted) return hinted;
    // No style comment: the count of numeric columns decides. Image flags add
    // three integer columns after the style's own, so a width one step above a
    // common style can still be that style with flags; reading an atomic row
    // with flags as `full` would shift x onto y.
    let width = 0;
    for (const row of rows) width = Math.max(width, row.length);
    if (width >= 10) return { kind: 'full', x: 4, type: 2 };
    if (width === 9 || width === 6) return chargeOrMolecular(rows, atomTypes);
    if (width === 8) return { kind: 'atomic', x: 2, type: 1 };
    if (width === 7) return { kind: 'full', x: 4, type: 2 };
    return { kind: 'atomic', x: 2, type: 1 };
}

function declaredAtomTypes(lines: string[]): number | undefined {
    for (const line of lines) {
        const m = /^\s*(\d+)\s+atom\s+types\s*$/i.exec(line);
        if (m) return parseInt(m[1], 10);
    }
    return undefined;
}

interface DataAtoms {
    rows: { id: number; type: number; x: number; y: number; z: number }[];
    style: AtomStyle;
    next: number;
}

function readDataAtoms(lines: string[], start: number, styleHint?: string, atomTypes?: number): DataAtoms {
    const header = lines[start].trim();
    const rawRows: string[][] = [];
    let i = start + 1;
    for (; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line === '') {
            let k = i + 1;
            while (k < lines.length && lines[k].trim() === '') k++;
            if (k >= lines.length || SECTION_START.test(lines[k].trim())) { i = k; break; }
            continue;
        }
        // A comment row is not an atom and must not seed the style guess.
        if (line.startsWith('#')) continue;
        if (SECTION_START.test(line)) break;
        const parts = line.split(/\s+/);
        if (parts.length < 5) continue;
        rawRows.push(parts);
    }
    const style: AtomStyle = rawRows.length > 0
        ? atomStyleOf(header, rawRows, styleHint, atomTypes)
        : { kind: 'atomic', x: 2, type: 1 };
    const rows: { id: number; type: number; x: number; y: number; z: number }[] = [];
    for (const parts of rawRows) {
        const id = parseInt(parts[0], 10);
        const type = parseInt(parts[style.type], 10);
        const x = parseFloat(parts[style.x]);
        const y = parseFloat(parts[style.x + 1]);
        const z = parseFloat(parts[style.x + 2]);
        if (!isFinite(id) || !isFinite(type) || !isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
        rows.push({ id, type, x, y, z });
    }
    return { rows, style, next: i };
}

export function parseLammpsData(content: string, options?: LammpsParseOptions): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const title = oneLineTitle(lines[0] || '', 'LAMMPS data');
    const unitFactor = lengthFactor(options && options.units);
    const masses = parseMasses(lines);
    const typeElements: { [type: number]: string } = { ...(options && options.typeElements ? options.typeElements : {}) };
    for (const key of Object.keys(masses)) {
        const type = Number(key);
        if (!typeElements[type]) typeElements[type] = masses[type];
    }

    const box = readDataBox(lines);
    if (!box.lattice) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the LAMMPS data file states no usable box: no structure to show']
        };
    }

    let atomsSection = -1;
    for (let i = 0; i < lines.length; i++) {
        if (/^Atoms\b/i.test(lines[i].trim())) { atomsSection = i; break; }
    }
    if (atomsSection < 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the LAMMPS data file has no Atoms section: no structure to show']
        };
    }
    const data = readDataAtoms(lines, atomsSection, options && options.atomStyle, declaredAtomTypes(lines));
    if (data.rows.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the Atoms section of the LAMMPS data file has no readable rows: no structure to show']
        };
    }

    // The box bounds are where the cell's origin sits; an atom coordinate is
    // absolute, so subtracting the origin expresses it against the cell the
    // lattice below is built at.
    const sites: Site[] = data.rows.map(row => {
        const element = typeElements[row.type] || typeElementLabel(row.type);
        return {
            element,
            cart: [
                (row.x - box.xlo) * unitFactor,
                (row.y - box.ylo) * unitFactor,
                (row.z - box.zlo) * unitFactor
            ],
            label: element
        };
    });
    const missing = data.rows.filter(row => !typeElements[row.type]);
    if (missing.length > 0) {
        const types = Array.from(new Set(missing.map(row => row.type))).sort((a, b) => a - b);
        warnings.push('atom type' + (types.length === 1 ? ' ' : 's ') + types.join(', ') +
            ' has no mass or element hint: read as ' + types.map(t => typeElementLabel(t)).join(', '));
    }
    const lattice = box.lattice.map(row => row.map(v => v * unitFactor));
    if (!isUsableLattice(lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the LAMMPS box is degenerate: no structure to show']
        };
    }
    return buildPeriodicStructure(lattice, sites, title, warnings);
}

interface DumpBox {
    lattice: number[][] | null;
    /** The three raw header lines, already converted to Angstrom. */
    xlo: number; xhi: number;
    ylo: number; yhi: number;
    zlo: number; zhi: number;
    tilt: [number, number, number];
}

function readDumpBox(lines: string[], start: number, factor: number): { box: DumpBox; next: number } {
    const values: number[][] = [];
    let i = start + 1;
    while (i < lines.length && values.length < 3) {
        const parts = lines[i].trim().split(/\s+/).map(Number);
        if (parts.length < 2 || parts.some(v => !isFinite(v))) break;
        values.push(parts);
        i++;
    }
    let xlo = 0, xhi = 0, ylo = 0, yhi = 0, zlo = 0, zhi = 0;
    let tilt: [number, number, number] = [0, 0, 0];
    if (values.length >= 3) {
        if (values[0].length >= 3) {
            // Triclinic dump: xlo_bound xhi_bound xy / ... / ...
            const xy = values[0][2] * factor;
            const xz = values[1][2] * factor;
            const yz = values[2][2] * factor;
            const xloB = values[0][0] * factor;
            const xhiB = values[0][1] * factor;
            const yloB = values[1][0] * factor;
            const yhiB = values[1][1] * factor;
            const zloB = values[2][0] * factor;
            const zhiB = values[2][1] * factor;
            // LAMMPS: xlo = xlo_bound - min(0, xy, xz, xy+xz), and the upper
            // bound subtracts the maximum instead.
            xlo = xloB - Math.min(0, xy, xz, xy + xz);
            xhi = xhiB - Math.max(0, xy, xz, xy + xz);
            ylo = yloB - Math.min(0, yz);
            yhi = yhiB - Math.max(0, yz);
            zlo = zloB;
            zhi = zhiB;
            tilt = [xy, xz, yz];
        } else {
            xlo = values[0][0] * factor;
            xhi = values[0][1] * factor;
            ylo = values[1][0] * factor;
            yhi = values[1][1] * factor;
            zlo = values[2][0] * factor;
            zhi = values[2][1] * factor;
        }
    }
    const lx = xhi - xlo;
    const ly = yhi - ylo;
    const lz = zhi - zlo;
    const lattice = lx > 0 && ly > 0 && lz > 0
        ? [[lx, 0, 0], [tilt[0], ly, 0], [tilt[1], tilt[2], lz]]
        : null;
    return { box: { lattice, xlo, xhi, ylo, yhi, zlo, zhi, tilt }, next: i };
}

function findColumn(columns: string[], names: string[]): number {
    for (const name of names) {
        const idx = columns.indexOf(name);
        if (idx >= 0) return idx;
    }
    return -1;
}

interface DumpFrame {
    timestep: string;
    lattice: number[][] | null;
    sites: Site[];
    /** The box in Angstrom, for scaled coordinates. */
    xlo: number; xhi: number; ylo: number; yhi: number; zlo: number; zhi: number;
    tilt: [number, number, number];
}

function parseDumpFrames(content: string, options?: LammpsParseOptions): { first: DumpFrame | null; count: number } {
    const lines = content.split(/\r?\n/);
    const factor = lengthFactor(options && options.units);
    const typeElements = options && options.typeElements ? options.typeElements : {};
    // The viewer draws one frame; a long trajectory is counted, not held in
    // memory frame by frame.
    let first: DumpFrame | null = null;
    let frameCount = 0;
    let i = 0;
    while (i < lines.length) {
        if (!/^ITEM:\s*TIMESTEP/i.test(lines[i].trim())) { i++; continue; }
        const timestep = (lines[i + 1] || '').trim();
        const numberHeader = i + 2;
        if (numberHeader >= lines.length || !/^ITEM:\s*NUMBER OF ATOMS/i.test(lines[numberHeader].trim())) {
            i = numberHeader;
            continue;
        }
        const count = parseInt((lines[numberHeader + 1] || '').trim(), 10);
        const boxHeader = numberHeader + 2;
        if (boxHeader >= lines.length || !/^ITEM:\s*BOX BOUNDS/i.test(lines[boxHeader].trim())) {
            i = boxHeader;
            continue;
        }
        const boxRead = readDumpBox(lines, boxHeader, factor);
        i = boxRead.next;
        if (i >= lines.length || !/^ITEM:\s*ATOMS/i.test(lines[i].trim())) continue;
        const columns = lines[i].trim().replace(/^ITEM:\s*ATOMS\s*/i, '').trim().split(/\s+/);
        i++;
        const col = (names: string[]) => findColumn(columns, names);
        const cId = col(['id']);
        const cType = col(['type']);
        const cElement = col(['element']);
        // A dump may carry both wrapped (`x`) and scaled (`xs`) columns; one
        // consistent triple decides both the index and whether the numbers are
        // fractions. Picking each axis separately could read `x` as scaled just
        // because an unused `xs` column exists elsewhere in the header.
        const groups: { names: string[]; scaled: boolean }[] = [
            { names: ['x', 'y', 'z'], scaled: false },
            { names: ['xu', 'yu', 'zu'], scaled: false },
            { names: ['xs', 'ys', 'zs'], scaled: true },
            { names: ['xsu', 'ysu', 'zsu'], scaled: true }
        ];
        const group = groups.find(g => g.names.every(n => columns.indexOf(n) >= 0));
        const cX = group ? columns.indexOf(group.names[0]) : -1;
        const cY = group ? columns.indexOf(group.names[1]) : -1;
        const cZ = group ? columns.indexOf(group.names[2]) : -1;
        const xsScaled = group ? group.scaled : false;
        if (cX < 0 || cY < 0 || cZ < 0) {
            const skip = isFinite(count) ? Math.max(0, Math.floor(count)) : 0;
            i += skip;
            continue;
        }
        const sites: Site[] = [];
        const limit = isFinite(count) && count > 0 ? count : Infinity;
        for (let row = 0; row < limit && i < lines.length; row++, i++) {
            const parts = lines[i].trim().split(/\s+/);
            if (parts.length <= cZ) break;
            if (cId >= 0 && !isNumeric(parts[cId])) break;
            const px = parseFloat(parts[cX]);
            const py = parseFloat(parts[cY]);
            const pz = parseFloat(parts[cZ]);
            if (!isFinite(px) || !isFinite(py) || !isFinite(pz)) break;
            const type = cType >= 0 ? parseInt(parts[cType], 10) : 1;
            const label = cElement >= 0 ? parts[cElement] : '';
            const element = (label ? normalizeElementSymbol(label) : undefined) ||
                typeElements[type] || typeElementLabel(isFinite(type) ? type : 1);
            let x: number, y: number, z: number;
            if (xsScaled) {
                // Scaled coordinates are fractions of the *restricted* box, so
                // the cell origin has to be added back; omitting it would move
                // every atom of a box whose xlo/ylo/zlo is not zero.
                const cart = fracToCart([px, py, pz], boxRead.box.lattice || [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
                x = cart[0] + boxRead.box.xlo;
                y = cart[1] + boxRead.box.ylo;
                z = cart[2] + boxRead.box.zlo;
            } else {
                // Like a data file, the dump's coordinates are absolute and the
                // lattice below is built at the box origin.
                x = px * factor - boxRead.box.xlo;
                y = py * factor - boxRead.box.ylo;
                z = pz * factor - boxRead.box.zlo;
            }
            sites.push({ element, cart: [x, y, z], label: label || element });
        }
        if (sites.length > 0) {
            frameCount++;
            if (first === null) {
                first = {
                    timestep,
                    lattice: boxRead.box.lattice,
                    sites,
                    xlo: boxRead.box.xlo, xhi: boxRead.box.xhi,
                    ylo: boxRead.box.ylo, yhi: boxRead.box.yhi,
                    zlo: boxRead.box.zlo, zhi: boxRead.box.zhi,
                    tilt: boxRead.box.tilt
                };
            }
        }
    }
    return { first, count: frameCount };
}

export function parseLammpsDump(content: string, options?: LammpsParseOptions): MolecularData {
    const parsed = parseDumpFrames(content, options);
    const warnings: string[] = [];
    const title = 'LAMMPS dump';
    const first = parsed.first;
    if (first === null) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['no complete frame in the LAMMPS dump: no structure to show']
        };
    }
    if (parsed.count > 1) {
        warnings.push(parsed.count + ' frames in the dump: showing the first (timestep ' + first.timestep + ')');
    }
    if (!first.lattice || !isUsableLattice(first.lattice)) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the dump states no usable box: no structure to show']
        };
    }
    const missing = first.sites.filter(s => /^T\d+$/.test(s.element));
    if (missing.length > 0) {
        const labels = Array.from(new Set(missing.map(s => s.element))).sort();
        warnings.push('atom type' + (labels.length === 1 ? ' ' : 's ') + labels.join(', ') +
            ' has no element column or companion data: read as typed labels');
    }
    return buildPeriodicStructure(first.lattice, first.sites, title + ' (timestep ' + first.timestep + ')', warnings);
}

export interface LammpsHints {
    units?: string;
    typeMasses: { [type: number]: number };
    typeElements: { [type: number]: string };
    readData?: string;
    atomStyle?: string;
}

/**
 * The type names and units a LAMMPS input script states.
 *
 * `dump_modify ... element Si Mg` is the one place a dump's types are named
 * explicitly; `mass 1 28.0855` names them by mass; `units metal` decides
 * whether the numbers in the dump are Angstrom, nanometres or metres.
 */
export function parseLammpsInputHints(content: string): LammpsHints {
    const hints: LammpsHints = { typeMasses: {}, typeElements: {} };
    const lines = content.split(/\r?\n/);
    for (const raw of lines) {
        const line = raw.replace(/#.*$/, '').trim();
        if (line === '') continue;
        const unit = /^units\s+(\S+)/i.exec(line);
        if (unit) { hints.units = unit[1].toLowerCase(); continue; }
        const mass = /^mass\s+(\d+)\s+([-+0-9.eE]+)/i.exec(line);
        if (mass) {
            const type = parseInt(mass[1], 10);
            const value = parseFloat(mass[2]);
            if (isFinite(type) && isFinite(value)) {
                hints.typeMasses[type] = value;
                const element = elementFromMass(value);
                if (element) hints.typeElements[type] = element;
            }
            continue;
        }
        const style = /^atom_style\s+(\S+)/i.exec(line);
        if (style) { hints.atomStyle = style[1].toLowerCase(); continue; }
        const readData = /^read_data\s+(\S+)/i.exec(line);
        if (readData) { hints.readData = readData[1]; continue; }
        const elements = /^dump_modify\b[^\n]*?\belement\s+(.+)$/i.exec(line);
        if (elements) {
            const list = elements[1].trim().split(/\s+/);
            for (let i = 0; i < list.length; i++) {
                const element = normalizeElementSymbol(list[i]);
                // The element list is positional and ends at the next
                // dump_modify keyword; a trailing `format %0.8f` would
                // otherwise name a further type fluorine.
                if (!element) break;
                hints.typeElements[i + 1] = element;
            }
        }
    }
    return hints;
}

export function isLammpsData(content: string): boolean {
    return /^\s*\d+\s+atoms\s*$/im.test(content) &&
        /^\s*\d+\s+atom types\s*$/im.test(content) &&
        /^\s*Atoms\b/im.test(content) &&
        /^\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+xlo\s+xhi\s*$/im.test(content);
}

export function isLammpsDump(content: string): boolean {
    return /^ITEM:\s*TIMESTEP\s*$/im.test(content) &&
        /^ITEM:\s*ATOMS\b/im.test(content) &&
        /^ITEM:\s*BOX BOUNDS\b/im.test(content);
}
