import { Atom, Bond, MolecularData } from '../types';
import { elementFromMass, normalizeElementSymbol, typeElementLabel } from '../tools/periodicTable';
import { cartToFrac, fracToCart } from '../tools/crystalGeometry';
import { buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';
import { chooseTrajectoryFrameIndexes } from './trajectoryParser';

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
    /**
     * Order the atoms by their LAMMPS ID before building the structure.
     *
     * A binary `dump dcd`/`dump xtc` always writes atoms sorted by ID, so a
     * data file whose Atoms section is not in ID order would otherwise pair
     * another atom's element and bonds with each trajectory position.
     */
    orderById?: boolean;
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

/**
 * Angstrom per LAMMPS length unit for a `units` value.
 *
 * A text dump and a data file state coordinates in the run's units, so the
 * loader scales them with this. A binary `dump dcd` also writes raw units;
 * `dump xtc` always writes nanometres (LAMMPS converts even in `units lj`),
 * so its reader keeps the XTC default of ten.
 */
export function lammpsLengthFactor(units: string | undefined): number {
    if (!units) return 1;
    const f = UNIT_TO_ANGSTROM[units.toLowerCase()];
    return f !== undefined ? f : 1;
}

function isNumeric(token: string | undefined): boolean {
    return token !== undefined && token !== '' && isFinite(parseFloat(token));
}

/** A LAMMPS line with its `#` comment removed; a `#` inside quotes is text. */
function stripLammpsComment(line: string): string {
    let quote = '';
    for (let i = 0; i < line.length; i++) {
        const ch = line[i];
        if (quote) {
            if (ch === quote) quote = '';
        } else if (ch === '"' || ch === "'") {
            quote = ch;
        } else if (ch === '#') {
            return line.slice(0, i);
        }
    }
    return line;
}

/** Words of a LAMMPS command, with quoted arguments kept as one word. */
function lammpsWords(line: string): string[] {
    const words: string[] = [];
    const re = /"([^"]*)"|'([^']*)'|(\S+)/g;
    let m: RegExpExecArray | null;
    while ((m = re.exec(line)) !== null) {
        words.push(m[1] !== undefined ? m[1] : m[2] !== undefined ? m[2] : m[3]);
    }
    return words;
}

/** `${name}` / `$name` substitution from the script's own `variable` commands. */
function expandLammpsVariables(token: string, variables: { [name: string]: string }): string {
    let value = token;
    for (let pass = 0; pass < 10; pass++) {
        const next = value
            .replace(/\$\{([A-Za-z_][A-Za-z0-9_]*)(?::[^}]*)?\}/g, (whole, name) =>
                Object.prototype.hasOwnProperty.call(variables, name) ? variables[name] : whole)
            .replace(/\$([A-Za-z_][A-Za-z0-9_]*)/g, (whole, name) =>
                Object.prototype.hasOwnProperty.call(variables, name) ? variables[name] : whole);
        if (next === value) break;
        value = next;
    }
    return value;
}

function defineLammpsVariable(words: string[], variables: { [name: string]: string }): boolean {
    if (words.length < 3 || words[0].toLowerCase() !== 'variable') return false;
    const style = words[2].toLowerCase();
    if (style !== 'string' && style !== 'index') return false;
    // An `index` variable holds the first value; a `string` keeps the rest.
    const value = style === 'index' ? (words[3] || '') : words.slice(3).join(' ');
    const expanded = expandLammpsVariables(value, variables);
    if (expanded !== '') variables[words[1]] = expanded;
    return true;
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

function lammpsDataLine(line: string): string {
    return stripLammpsComment(line).trim();
}

function readDataBox(lines: string[]): DataBox {
    const bound: { [key: string]: number } = {};
    let tilt: [number, number, number] = [0, 0, 0];
    for (const raw of lines) {
        const line = lammpsDataLine(raw);
        const m = /^([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+(xlo|ylo|zlo)\s+(xhi|yhi|zhi)\b/i.exec(line);
        if (m) {
            const lo = m[3].toLowerCase();
            const hi = m[4].toLowerCase();
            bound[lo] = parseFloat(m[1]);
            bound[hi] = parseFloat(m[2]);
            continue;
        }
        const t = /^([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+xy\s+xz\s+yz\b/i.exec(line);
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
    for (const raw of lines) {
        const line = lammpsDataLine(raw);
        const m = /^(\d+)\s+atom\s+types$/i.exec(line);
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
        // A trailing comment is not a column: keeping it would inflate the
        // row's width and let the style guess read x from the wrong field.
        const line = lammpsDataLine(lines[i]);
        if (line === '') continue;
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

/** The line index of the `Bonds` section, or -1. */
function findBondSection(lines: string[]): number {
    // The section header is the word alone (optionally followed by a `#`
    // comment). A title line that merely starts with "Bonds" is not a section.
    for (let i = 1; i < lines.length; i++) {
        if (/^Bonds(?:\s*#.*)?\s*$/i.test(lines[i].trim())) return i;
    }
    return -1;
}

/** The bond count the file states (`N bonds`), or undefined when it states none. */
function declaredBonds(lines: string[]): number | undefined {
    for (const raw of lines) {
        const line = lammpsDataLine(raw);
        const m = /^(\d+)\s+bonds$/i.exec(line);
        if (m) return parseInt(m[1], 10);
    }
    return undefined;
}

/**
 * A data file's `Bonds` rows: `id type atom1 atom2` (the four leading columns
 * are the same in the atomic/bond/molecular/full styles and in hybrid styles).
 *
 * LAMMPS states connectivity but not a chemical bond order: bond *types* are
 * force-field parameter indices. A file that says nothing about order therefore
 * carries order 0, which the crystal page reads as "the perceived order stands".
 * The writer here labels its own table `Bonds # order`, and only that label
 * makes the type an order (1 single, 2 double, 3 triple, 4 aromatic) — reading
 * it on a force-field file would relabel every bond.
 */
function readDataBonds(
    lines: string[], start: number, idToIndex: { [id: number]: number }
): { atom1: number; atom2: number; order: number }[] {
    const byOrder = /#\s*order\b/i.test(lines[start] || '');
    const rows: { atom1: number; atom2: number; order: number }[] = [];
    const seen = new Set<string>();
    for (let i = start + 1; i < lines.length; i++) {
        const line = lammpsDataLine(lines[i]);
        if (line === '') continue;
        if (SECTION_START.test(line)) break;
        const parts = line.split(/\s+/);
        if (parts.length < 4) continue;
        const id1 = parseInt(parts[2], 10);
        const id2 = parseInt(parts[3], 10);
        const a1 = idToIndex[id1];
        const a2 = idToIndex[id2];
        if (a1 === undefined || a2 === undefined || a1 === a2) continue;
        let order = 0;
        if (byOrder) {
            const type = parseInt(parts[1], 10);
            order = type === 2 ? 2 : type === 3 ? 3 : type === 4 ? 1.5 : 1;
        }
        const key1 = a1 < a2 ? a1 : a2;
        const key2 = a1 < a2 ? a2 : a1;
        const key = key1 + '-' + key2;
        if (seen.has(key)) continue;
        seen.add(key);
        rows.push({ atom1: key1, atom2: key2, order });
    }
    return rows;
}

/**
 * The lattice repeat that puts atom b nearest atom a — the shift a periodic
 * bond stores (`atom1` in cell c bonds `atom2` in cell c+s`).
 *
 * A data file has no shift column: LAMMPS applies the minimum image itself, so
 * the same convention is recovered here by enumerating the rounded offset and
 * its neighbours (a skewed cell's nearest image can be one repeat off the
 * rounding, which is why the 3×3×3 neighbourhood is scanned).
 */
function minimumImageShift(a: Atom, b: Atom, lattice: number[][]): [number, number, number] {
    const fa = cartToFrac([a.x, a.y, a.z], lattice);
    const fb = cartToFrac([b.x, b.y, b.z], lattice);
    const df = [fb[0] - fa[0], fb[1] - fa[1], fb[2] - fa[2]];
    const base = [-Math.round(df[0]), -Math.round(df[1]), -Math.round(df[2])];
    let best: [number, number, number] = [0, 0, 0];
    let bestD2 = Infinity;
    for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
            for (let dz = -1; dz <= 1; dz++) {
                const s: [number, number, number] = [base[0] + dx, base[1] + dy, base[2] + dz];
                const p = fracToCart([df[0] + s[0], df[1] + s[1], df[2] + s[2]], lattice);
                const d2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
                if (d2 < bestD2 - 1e-12) {
                    bestD2 = d2;
                    best = s;
                }
            }
        }
    }
    return best;
}

/** A data file's bonds as the cell's file bond table (`crystal.fileBonds`). */
function dataBondsToFileBonds(
    rows: { atom1: number; atom2: number; order: number }[],
    atoms: Atom[], lattice: number[][]
): Bond[] {
    const bonds: Bond[] = [];
    for (const row of rows) {
        const a = atoms[row.atom1];
        const b = atoms[row.atom2];
        if (!a || !b) continue;
        const shift = minimumImageShift(a, b, lattice);
        bonds.push({ atom1: row.atom1, atom2: row.atom2, order: row.order, shift });
    }
    return bonds;
}

export function parseLammpsData(content: string, options?: LammpsParseOptions): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    const title = oneLineTitle(lines[0] || '', 'LAMMPS data');
    const unitFactor = lammpsLengthFactor(options && options.units);
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
    if (options && options.orderById) {
        data.rows.sort((a, b) => a.id - b.id);
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
    const structure = buildPeriodicStructure(lattice, sites, title, warnings);
    // A data file states connectivity in its `Bonds` section; the format the
    // cell is built in carries it as the file bond table, which the crystal page
    // merges over its perceived list (and keeps the file's bond orders where the
    // table states them). Without this a saved LAMMPS file re-opened with no
    // bonds at all.
    const bondSection = findBondSection(lines);
    if (bondSection >= 0 && structure.crystal) {
        const idToIndex: { [id: number]: number } = {};
        data.rows.forEach((row, i) => { idToIndex[row.id] = i; });
        const rows = readDataBonds(lines, bondSection, idToIndex);
        const declared = declaredBonds(lines);
        if (declared !== undefined && declared > rows.length) {
            warnings.push('the file declares ' + declared + ' bonds but only ' + rows.length +
                ' could be read: the rest name atoms the Atoms section does not contain');
        }
        if (rows.length > 0) {
            structure.crystal.fileBonds = dataBondsToFileBonds(rows, structure.atoms, lattice);
        }
    }
    // The reader can add to the warnings after the structure was built, so the
    // list is attached here as well (it is the same array in either place).
    if (warnings.length > 0) structure.warnings = warnings;
    return structure;
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
    let lattice: number[][] | null = null;
    const general = /\babc\s+origin\b/i.test(lines[start] || '');
    if (values.length >= 3) {
        if (general || (values[0].length >= 4 && values[1].length >= 4 && values[2].length >= 4)) {
            // General triclinic (`dump_modify triclinic/general`): each line is
            // one edge vector followed by that axis's origin component,
            //   ax ay az originx / bx by bz originy / cx cy cz originz
            // and the vectors are already in Cartesian form, not lower
            // triangular, so they are taken as the lattice exactly as written.
            const a = values[0].slice(0, 3).map(v => v * factor);
            const b = values[1].slice(0, 3).map(v => v * factor);
            const c = values[2].slice(0, 3).map(v => v * factor);
            lattice = isUsableLattice([a, b, c]) ? [a, b, c] : null;
            xlo = (values[0][3] || 0) * factor;
            ylo = (values[1][3] || 0) * factor;
            zlo = (values[2][3] || 0) * factor;
            xhi = xlo + Math.hypot(a[0], a[1], a[2]);
            yhi = ylo + Math.hypot(b[0], b[1], b[2]);
            zhi = zlo + Math.hypot(c[0], c[1], c[2]);
            tilt = [b[0], c[0], c[1]];
        } else if (values[0].length >= 3) {
            // Restricted triclinic: xlo_bound xhi_bound xy / ... / ...
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
    if (!lattice) {
        const lx = xhi - xlo;
        const ly = yhi - ylo;
        const lz = zhi - zlo;
        lattice = lx > 0 && ly > 0 && lz > 0
            ? [[lx, 0, 0], [tilt[0], ly, 0], [tilt[1], tilt[2], lz]]
            : null;
    }
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
    /** Whether the block held every atom its header declared. */
    complete: boolean;
}

interface DumpHeader {
    timestep: string;
    count: number;
    box: DumpBox;
    /** First atom row (the line after `ITEM: ATOMS`). */
    body: number;
    columns: string[];
    /** Angstrom per coordinate unit, from the run's `units`. */
    factor: number;
}

/**
 * Every complete `ITEM: TIMESTEP … BOX … ATOMS` header in the file.
 *
 * Only the headers are read; the atom rows are left for the frames the caller
 * actually wants, which is what makes a multi-frame dump affordable to open.
 * The scan always looks for the next `ITEM: TIMESTEP` instead of stepping by
 * the declared atom count, so one short or overlong block does not hide the
 * frames after it.
 */
function scanDumpHeaders(lines: string[], factor: number): DumpHeader[] {
    const headers: DumpHeader[] = [];
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
        if (!isFinite(count) || count <= 0) { i = numberHeader + 2; continue; }
        const boxHeader = numberHeader + 2;
        if (boxHeader >= lines.length || !/^ITEM:\s*BOX BOUNDS/i.test(lines[boxHeader].trim())) {
            i = boxHeader;
            continue;
        }
        const boxRead = readDumpBox(lines, boxHeader, factor);
        const atomsHeader = boxRead.next;
        if (atomsHeader >= lines.length || !/^ITEM:\s*ATOMS/i.test(lines[atomsHeader].trim())) {
            i = atomsHeader;
            continue;
        }
        const columns = lines[atomsHeader].trim().replace(/^ITEM:\s*ATOMS\s*/i, '').trim().split(/\s+/);
        headers.push({ timestep, count, box: boxRead.box, body: atomsHeader + 1, columns, factor });
        i = atomsHeader + 1;
    }
    return headers;
}

function buildDumpFrame(lines: string[], header: DumpHeader, typeElements: { [type: number]: string }): DumpFrame {
    const columns = header.columns;
    const box = header.box;
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
    const rows: { id: number; site: Site }[] = [];
    if (!group) {
        return {
            timestep: header.timestep, lattice: box.lattice, sites: [],
            xlo: box.xlo, xhi: box.xhi, ylo: box.ylo, yhi: box.yhi,
            zlo: box.zlo, zhi: box.zhi, tilt: box.tilt, complete: false
        };
    }
    const cX = columns.indexOf(group.names[0]);
    const cY = columns.indexOf(group.names[1]);
    const cZ = columns.indexOf(group.names[2]);
    const xsScaled = group.scaled;
    const factor = header.factor;
    let i = header.body;
    for (let row = 0; row < header.count && i < lines.length; row++, i++) {
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
            // A scaled coordinate is a fraction of the restricted box, and the
            // lattice is built at that box's origin, so the fraction is
            // already relative to the cell: adding xlo back would translate
            // every atom of a box that does not start at zero.
            const cart = fracToCart([px, py, pz], box.lattice || [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
            x = cart[0];
            y = cart[1];
            z = cart[2];
        } else {
            // Like a data file, the dump's coordinates are absolute and the
            // lattice below is built at the box origin.
            x = px * factor - box.xlo;
            y = py * factor - box.ylo;
            z = pz * factor - box.zlo;
        }
        rows.push({ id: cId >= 0 ? parseInt(parts[cId], 10) : row, site: { element, cart: [x, y, z], label: label || element } });
    }
    // An ID column gives every frame the same atom order even when LAMMPS
    // wrote the rows out of order; without one the file order is the only
    // correspondence the dump offers. A data-file topology read with
    // `orderById` lines up with the sorted order.
    if (cId >= 0) rows.sort((a, b) => a.id - b.id);
    const sites = rows.map(r => r.site);
    return {
        timestep: header.timestep, lattice: box.lattice, sites,
        xlo: box.xlo, xhi: box.xhi, ylo: box.ylo, yhi: box.yhi,
        zlo: box.zlo, zhi: box.zhi, tilt: box.tilt,
        complete: sites.length === header.count
    };
}

/** The first frame with readable atoms, and how many frame headers exist. */
function parseDumpFrames(content: string, options?: LammpsParseOptions): { first: DumpFrame | null; count: number } {
    const lines = content.split(/\r?\n/);
    const factor = lammpsLengthFactor(options && options.units);
    const typeElements = options && options.typeElements ? options.typeElements : {};
    const headers = scanDumpHeaders(lines, factor);
    for (const header of headers) {
        const frame = buildDumpFrame(lines, header, typeElements);
        if (frame.sites.length > 0) return { first: frame, count: headers.length };
    }
    return { first: null, count: headers.length };
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
    if (!first.complete) {
        warnings.push('the frame at timestep ' + first.timestep + ' is incomplete: showing the atoms it held');
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

export interface LammpsDumpFrame {
    /** 0-based index of the frame among the blocks in the file. */
    index: number;
    timestep: string;
    /** The timestep as a number, when it is an integer. */
    step?: number;
    sites: Site[];
    lattice: number[][];
}

export interface LammpsDumpTrajectory {
    frames: LammpsDumpFrame[];
    /** Atom count of the first frame; every kept frame has this count. */
    natoms: number;
    /** Frames the dump contains, before sampling. */
    frameCount: number;
    sampled: boolean;
    truncated: boolean;
    warnings: string[];
}

export interface LammpsDumpTrajectoryOptions extends LammpsParseOptions {
    maxFrames?: number;
    maxTotalAtoms?: number;
}

/**
 * A text LAMMPS dump as a trajectory.
 *
 * A dump is ASCII, so there is no random access to seek through; the scan
 * still reads only the three `ITEM:` header lines of every block and decodes
 * the atom rows of the sampled frames. Frames that are incomplete, state no
 * usable box, or change the atom count are skipped with a warning rather than
 * producing a frame whose positions would be paired with the wrong atoms.
 */
export function parseLammpsDumpTrajectory(
    content: string,
    options?: LammpsDumpTrajectoryOptions
): LammpsDumpTrajectory {
    const lines = content.split(/\r?\n/);
    const factor = lammpsLengthFactor(options && options.units);
    const typeElements = options && options.typeElements ? options.typeElements : {};
    const headers = scanDumpHeaders(lines, factor);
    const result: LammpsDumpTrajectory = {
        frames: [],
        natoms: headers.length > 0 ? headers[0].count : 0,
        frameCount: headers.length,
        sampled: false,
        truncated: false,
        warnings: []
    };
    if (headers.length === 0) {
        result.warnings.push('no complete frame in the LAMMPS dump: no structure to show');
        return result;
    }
    const last = headers[headers.length - 1];
    result.truncated = last.body + last.count > lines.length;
    if (result.truncated) {
        result.warnings.push('the dump ends in an incomplete frame: showing the frames that could be read');
    }
    const chosen = chooseTrajectoryFrameIndexes(
        headers.length, result.natoms,
        options && options.maxFrames, options && options.maxTotalAtoms
    );
    result.sampled = chosen.length < headers.length;

    let incomplete = 0;
    let varying = 0;
    let noBox = 0;
    for (const index of chosen) {
        const header = headers[index];
        if (header.count !== result.natoms) { varying++; continue; }
        if (!header.box.lattice || !isUsableLattice(header.box.lattice)) { noBox++; continue; }
        const frame = buildDumpFrame(lines, header, typeElements);
        if (frame.sites.length !== header.count) { incomplete++; continue; }
        const parsed: LammpsDumpFrame = {
            index,
            timestep: frame.timestep,
            sites: frame.sites,
            lattice: frame.lattice as number[][]
        };
        if (/^[+-]?\d+$/.test(frame.timestep)) parsed.step = parseInt(frame.timestep, 10);
        result.frames.push(parsed);
    }
    if (incomplete > 0) {
        result.warnings.push(incomplete + ' incomplete frame' + (incomplete === 1 ? '' : 's') +
            ' in the dump: skipped');
    }
    if (varying > 0) {
        result.warnings.push(varying + ' frame' + (varying === 1 ? ' has' : 's have') +
            ' a different atom count: skipped');
    }
    if (noBox > 0) {
        result.warnings.push(noBox + ' frame' + (noBox === 1 ? '' : 's') +
            ' states no usable box: skipped');
    }
    const missing = result.frames.length > 0 ? result.frames[0].sites.filter(s => /^T\d+$/.test(s.element)) : [];
    if (missing.length > 0) {
        const labels = Array.from(new Set(missing.map(s => s.element))).sort();
        result.warnings.push('atom type' + (labels.length === 1 ? ' ' : 's ') + labels.join(', ') +
            ' has no element column or companion data: read as typed labels');
    }
    return result;
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
 * whether the numbers in the dump are Angstrom, nanometres or metres. Variables
 * defined in the script are expanded in the same words the command has, so
 * `variable data string si.data` followed by `read_data ${data}` names the file
 * the run actually read.
 */
export function parseLammpsInputHints(content: string, inherited?: { [name: string]: string }): LammpsHints {
    const hints: LammpsHints = { typeMasses: {}, typeElements: {} };
    const variables: { [name: string]: string } = { ...(inherited || {}) };
    for (const raw of content.split(/\r?\n/)) {
        const words = lammpsWords(stripLammpsComment(raw));
        if (words.length === 0) continue;
        if (defineLammpsVariable(words, variables)) continue;
        const command = words[0].toLowerCase();
        if (command === 'units' && words.length >= 2) {
            hints.units = expandLammpsVariables(words[1], variables).toLowerCase();
            continue;
        }
        if (command === 'mass' && words.length >= 3) {
            const type = parseInt(expandLammpsVariables(words[1], variables), 10);
            const value = parseFloat(expandLammpsVariables(words[2], variables));
            if (isFinite(type) && isFinite(value)) {
                hints.typeMasses[type] = value;
                const element = elementFromMass(value);
                if (element) hints.typeElements[type] = element;
            }
            continue;
        }
        if (command === 'atom_style' && words.length >= 2) {
            hints.atomStyle = expandLammpsVariables(words[1], variables).toLowerCase();
            continue;
        }
        if (command === 'read_data' && words.length >= 2) {
            hints.readData = expandLammpsVariables(words[1], variables);
            continue;
        }
        if (command === 'dump_modify') {
            let at = -1;
            for (let i = 1; i < words.length; i++) {
                if (words[i].toLowerCase() === 'element') { at = i; break; }
            }
            if (at < 0) continue;
            // The element list is positional and ends at the next dump_modify
            // keyword; a trailing `format %0.8f` would otherwise name a further
            // type fluorine.
            for (let i = at + 1; i < words.length; i++) {
                const element = normalizeElementSymbol(words[i]);
                if (!element) break;
                hints.typeElements[i - at] = element;
            }
        }
    }
    return hints;
}

export interface LammpsInputScan {
    /** Variables defined (or inherited) when the file finished. */
    variables: { [name: string]: string };
    /** `read_data` targets in command order, variables expanded. */
    readData: string[];
    /** `include` targets in command order, variables expanded. */
    includes: string[];
    /** The last `read_restart` target, variables expanded. */
    readRestart?: string;
}

/**
 * What a LAMMPS input script says about files outside it.
 *
 * The reader never touches the file system: the host follows the returned
 * names, reads those files, and hands the structure back through
 * `externalStructureFile`. Only the two variable styles whose value is a
 * literal (`string` and `index`) are expanded; an arithmetic `equal` or an
 * immediate `$(...)` expression is left as written and reported as an
 * unresolved reference instead of being guessed.
 */
export function scanLammpsInput(content: string, inherited?: { [name: string]: string }): LammpsInputScan {
    const variables: { [name: string]: string } = { ...(inherited || {}) };
    const readData: string[] = [];
    const includes: string[] = [];
    let readRestart: string | undefined;
    for (const raw of content.split(/\r?\n/)) {
        const words = lammpsWords(stripLammpsComment(raw));
        if (words.length === 0) continue;
        if (defineLammpsVariable(words, variables)) continue;
        const command = words[0].toLowerCase();
        if (command === 'read_data' && words.length >= 2) {
            readData.push(expandLammpsVariables(words[1], variables));
        } else if (command === 'include' && words.length >= 2) {
            includes.push(expandLammpsVariables(words[1], variables));
        } else if (command === 'read_restart' && words.length >= 2) {
            readRestart = expandLammpsVariables(words[1], variables);
        }
    }
    return { variables, readData, includes, readRestart };
}

// Commands that only a LAMMPS input uses at the start of a line. The words
// deliberately left out (`lattice`, `region`, `group`, `velocity`, `run`,
// `set`, `thermo`, `neighbor`, ...) also open output lines of ABACUS, QE,
// PySCF, GROMACS and GPUMD, so they cannot be evidence on their own.
const LAMMPS_STRONG_COMMANDS = new Set([
    'units', 'atom_style', 'atom_modify', 'read_data', 'read_restart', 'read_dump',
    'create_box', 'create_atoms', 'pair_style', 'pair_coeff', 'pair_modify',
    'bond_style', 'bond_coeff', 'angle_style', 'angle_coeff', 'dihedral_style',
    'dihedral_coeff', 'improper_style', 'improper_coeff', 'kspace_style',
    'kspace_modify', 'special_bonds', 'neigh_modify', 'comm_modify', 'timestep',
    'fix', 'unfix', 'compute', 'dump', 'undump', 'dump_modify', 'thermo_style',
    'thermo_modify', 'min_style', 'min_modify', 'minimize', 'reset_timestep',
    'write_data', 'write_restart', 'write_dump', 'delete_atoms', 'displace_atoms',
    'change_box', 'replicate', 'delete_bonds', 'create_bonds'
]);

/**
 * Whether a file is a LAMMPS input script rather than a structure.
 *
 * A data file starts with `N atoms` / `Atoms`, a dump with `ITEM:`, and neither
 * uses LAMMPS commands at line starts. One structural command (`read_data`,
 * `create_box`, ...) or two strong commands distinguish a script; one weak word
 * or an echoed `key = value` line does not, because ABACUS, QE, PySCF and
 * GROMACS outputs open lines with the same words.
 */
export function isLammpsLog(content: string): boolean {
    // Every LAMMPS log opens with the version banner; the same words in an
    // input script are a comment, so the test is anchored at the line start.
    return /^LAMMPS\s*\(/im.test(content);
}

export function isLammpsInput(content: string): boolean {
    // A LAMMPS log echoes the input commands it ran, so command-looking lines
    // are not evidence that this file is the input itself.
    if (isLammpsLog(content)) return false;
    let strong = 0;
    for (const raw of content.split(/\r?\n/)) {
        const line = stripLammpsComment(raw).trim();
        if (line === '') continue;
        // `lattice constant (Bohr) = 10.26` is printed by ABACUS/QE, not a
        // command; an input line rarely carries `=` in its first words.
        if (/^\S+(?:\s+\S+){0,3}\s*=\s*\S/.test(line)) continue;
        // LAMMPS commands are lowercase and case-sensitive. Accepting any case
        // made Gaussian's "Compute ..." / "Fix ..." section lines look like a
        // script, and this test runs before the output readers.
        const command = line.split(/\s+/)[0];
        if (!LAMMPS_STRONG_COMMANDS.has(command)) continue;
        if (/^(read_data|read_restart|create_box|create_atoms)$/.test(command)) return true;
        if (++strong >= 2) return true;
    }
    return false;
}

/**
 * A LAMMPS input as a structure: the data file its `read_data` names.
 *
 * The host is the only part that can resolve the path (relative, included,
 * variable-expanded), so it hands the referenced file in. When it could not,
 * the warning names what the script asked for instead of drawing a fragment of
 * the script itself.
 */
export function parseLammpsInputStructure(
    content: string,
    external: { content: string; fileName: string } | undefined,
    options?: LammpsParseOptions
): MolecularData {
    const title = oneLineTitle((content.split(/\r?\n/)[0] || '').replace(/^#\s*/, ''), 'LAMMPS input');
    if (external) {
        const hints = parseLammpsInputHints(content);
        const lammps: LammpsParseOptions = { ...(options || {}) };
        if (!lammps.units && hints.units) lammps.units = hints.units;
        if (!lammps.atomStyle && hints.atomStyle) lammps.atomStyle = hints.atomStyle;
        const typeElements = { ...hints.typeElements, ...(lammps.typeElements || {}) };
        if (Object.keys(typeElements).length > 0) lammps.typeElements = typeElements;
        return parseLammpsData(external.content, lammps);
    }

    const scan = scanLammpsInput(content);
    let warning: string;
    if (scan.readData.length > 0) {
        warning = 'the LAMMPS input reads its structure from "' + scan.readData[0] +
            '", but that file could not be read';
    } else if (scan.readRestart) {
        warning = 'the LAMMPS input continues from the restart file "' + scan.readRestart +
            '": a restart file holds no readable atom list';
    } else {
        warning = 'the LAMMPS input builds its structure at run time (no read_data): no structure file to show';
    }
    return { atoms: [], bonds: [], title, hasExplicitBonds: false, warnings: [warning] };
}

export function isLammpsData(content: string): boolean {
    let atoms = false;
    let types = false;
    let section = false;
    let box = false;
    for (const raw of content.split(/\r?\n/)) {
        const line = lammpsDataLine(raw);
        if (line === '') continue;
        if (!atoms && /^\d+\s+atoms$/i.test(line)) atoms = true;
        else if (!types && /^\d+\s+atom types$/i.test(line)) types = true;
        else if (!section && /^Atoms\b/i.test(line)) section = true;
        else if (!box && /^[-+0-9.eE]+\s+[-+0-9.eE]+\s+xlo\s+xhi/i.test(line)) box = true;
        if (atoms && types && section && box) return true;
    }
    return false;
}

export function isLammpsDump(content: string): boolean {
    return /^ITEM:\s*TIMESTEP\s*$/im.test(content) &&
        /^ITEM:\s*ATOMS\b/im.test(content) &&
        /^ITEM:\s*BOX BOUNDS\b/im.test(content);
}
