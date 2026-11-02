import { MolecularData } from '../types';
import { elementFromAtomName, normalizeElementSymbol } from '../tools/periodicTable';
import { cellToLatticeVectors, fracToCart } from '../tools/crystalGeometry';
import { parseXyz } from './xyzParser';
import { parsePdb } from './pdbParser';
import {
    BOHR_TO_ANGSTROM,
    buildMolecularStructure,
    buildPeriodicStructure,
    isUsableLattice,
    latticeToCellParams,
    oneLineTitle,
    Site
} from './periodicUtils';

/**
 * CP2K reader: the input (`.inp`) and the output (`.out`).
 *
 * CP2K states a structure in nested `&SECTION` blocks, so the reader keeps the
 * block structure instead of flattening the file: `&CELL` and `&COORD` are
 * what describe the geometry, and their keys (`ABC`, `A/B/C`, `SCALED`, `UNIT`)
 * decide what the numbers mean. An input that gets its atoms from
 * `&TOPOLOGY COORD_FILE_NAME` has no inline coordinates at all; the companion
 * file is read by the loader and handed in as `externalCoordinateFile`.
 *
 * A cell whose `PERIODIC` is `NONE` is a molecule in a box, not a crystal —
 * the water example is exactly that — so it is returned without a `CrystalData`
 * rather than drawn with a 20 Angstrom periodic frame the file says is absent.
 */

export interface Cp2kParseOptions {
    externalCoordinateFile?: { content: string; fileName: string };
}

interface Cp2kSection {
    name: string;
    arg: string;
    parent: Cp2kSection | null;
    keys: { [key: string]: string[] };
    raw: string[];
    children: Cp2kSection[];
}

function newSection(name: string, arg: string, parent: Cp2kSection | null): Cp2kSection {
    return { name: name.toUpperCase(), arg, parent, keys: {}, raw: [], children: [] };
}

function parseCp2kSections(content: string): Cp2kSection[] {
    const roots: Cp2kSection[] = [];
    const stack: Cp2kSection[] = [];
    const lines = content.split(/\r?\n/);
    for (const rawLine of lines) {
        let line = rawLine;
        const hash = line.indexOf('#');
        if (hash >= 0) line = line.slice(0, hash);
        // CP2K writes comments with `!` as well; a trailing comment after a
        // value is not part of the value.
        const bang = line.indexOf('!');
        if (bang >= 0) line = line.slice(0, bang);
        line = line.trim();
        if (line === '') continue;
        if (/^&END\b/i.test(line)) {
            if (stack.length > 0) stack.pop();
            continue;
        }
        if (/^&[A-Za-z]/.test(line)) {
            const m = /^&([A-Za-z0-9_]+)\s*(.*)$/.exec(line);
            if (!m) continue;
            const section = newSection(m[1], m[2] || '', stack.length > 0 ? stack[stack.length - 1] : null);
            if (stack.length > 0) stack[stack.length - 1].children.push(section);
            else roots.push(section);
            stack.push(section);
            continue;
        }
        const current = stack.length > 0 ? stack[stack.length - 1] : null;
        if (!current) continue;
        const tokens = line.split(/\s+/);
        // Inside &COORD an atom row is `element x y z ...`; a key row is
        // `KEYWORD value...`, and the three coordinate columns tell them apart.
        // Everywhere else a line like `A 5.0 0.0 0.0` is a key with a vector
        // value, so only &COORD may take the atom-row reading.
        if (current.name === 'COORD' && tokens.length >= 4 &&
            isFinite(parseFloat(tokens[1])) &&
            isFinite(parseFloat(tokens[2])) &&
            isFinite(parseFloat(tokens[3]))) {
            current.raw.push(line);
            continue;
        }
        if (/^[A-Za-z][A-Za-z0-9_]*$/.test(tokens[0])) {
            const key = tokens[0].toUpperCase();
            current.keys[key] = tokens.slice(1);
        } else {
            current.raw.push(line);
        }
    }
    return roots;
}

function flattenSections(roots: Cp2kSection[]): Cp2kSection[] {
    const out: Cp2kSection[] = [];
    const walk = (s: Cp2kSection) => {
        out.push(s);
        for (const c of s.children) walk(c);
    };
    for (const r of roots) walk(r);
    return out;
}

function lastSection(sections: Cp2kSection[], name: string): Cp2kSection | undefined {
    for (let i = sections.length - 1; i >= 0; i--) {
        if (sections[i].name === name) return sections[i];
    }
    return undefined;
}

/** The keywords CP2K accepts as parameters on a section's own line. */
const SECTION_PARAMS = ['SCALED', 'UNIT', 'PERIODIC', 'A', 'B', 'C', 'ABC', 'ALPHA_BETA_GAMMA'];

/**
 * A section parameter written on the `&NAME` line (`&COORD SCALED`,
 * `&CELL UNIT bohr`), read as if it were a key in the body: CP2K treats the
 * two spellings alike, and the body-only reader would miss the flag entirely.
 */
function sectionArgKeys(section: Cp2kSection): { [key: string]: string[] } {
    const out: { [key: string]: string[] } = {};
    const tokens = section.arg.trim() === '' ? [] : section.arg.trim().split(/\s+/);
    let i = 0;
    while (i < tokens.length) {
        const key = tokens[i].toUpperCase();
        i++;
        if (SECTION_PARAMS.indexOf(key) < 0) continue;
        const values: string[] = [];
        while (i < tokens.length && SECTION_PARAMS.indexOf(tokens[i].toUpperCase()) < 0) {
            values.push(tokens[i]);
            i++;
        }
        out[key] = values;
    }
    return out;
}

/** A section key, whether it was written in the body or on the section line. */
function sectionKey(section: Cp2kSection, key: string): string[] | undefined {
    if (section.keys[key] !== undefined) return section.keys[key];
    return sectionArgKeys(section)[key];
}

function sectionNumbers(key: string[] | undefined): number[] {
    if (!key) return [];
    const out: number[] = [];
    for (const t of key) {
        const v = parseFloat(t);
        if (isFinite(v)) out.push(v);
    }
    return out;
}

function cellLattice(sections: Cp2kSection[]): number[][] | null {
    const cell = lastSection(sections, 'CELL');
    if (!cell) return null;
    // &CELL UNIT applies to the cell vectors and defaults to Angstrom; a
    // SCALED &COORD is fractional regardless of it. Reading a Bohr cell as
    // Angstrom would draw the whole crystal at 1.89x its true size.
    const unitTokens = sectionKey(cell, 'UNIT');
    const unitText = unitTokens ? (unitTokens[0] || '').toLowerCase() : 'angstrom';
    const factor = unitText.indexOf('bohr') >= 0 || unitText.indexOf('a.u') >= 0 || unitText === 'au'
        ? BOHR_TO_ANGSTROM
        : unitText.indexOf('nm') >= 0 ? 10
        : unitText.indexOf('pm') >= 0 ? 0.01
        : 1;
    const scale = (lattice: number[][]) => lattice.map(row => row.map(v => v * factor));
    // `&CELL A 5 0 0` may be written on the section line, like any other
    // section parameter; reading only the body would miss that spelling.
    const aVec = sectionNumbers(sectionKey(cell, 'A'));
    const bVec = sectionNumbers(sectionKey(cell, 'B'));
    const cVec = sectionNumbers(sectionKey(cell, 'C'));
    if (aVec.length >= 3 && bVec.length >= 3 && cVec.length >= 3) {
        return scale([aVec.slice(0, 3), bVec.slice(0, 3), cVec.slice(0, 3)]);
    }
    const abc = sectionNumbers(sectionKey(cell, 'ABC'));
    if (abc.length >= 3) {
        const angles = sectionNumbers(sectionKey(cell, 'ALPHA_BETA_GAMMA'));
        const alpha = angles.length >= 1 ? angles[0] : 90;
        const beta = angles.length >= 2 ? angles[1] : 90;
        const gamma = angles.length >= 3 ? angles[2] : 90;
        return scale(cellToLatticeVectors(abc[0], abc[1], abc[2], alpha, beta, gamma));
    }
    return null;
}

function isPeriodic(sections: Cp2kSection[]): boolean {
    const cell = lastSection(sections, 'CELL');
    if (!cell) return false;
    const periodic = sectionKey(cell, 'PERIODIC');
    if (!periodic || periodic.length === 0) return true;
    const text = periodic.join(' ').toUpperCase();
    if (/NONE/.test(text)) return false;
    // CP2K also spells the periodic directions as a bitmask (1 = X, 2 = Y,
    // 4 = Z). Any non-zero value is a periodic cell.
    if (/^\d+$/.test(text.trim())) return parseInt(text.trim(), 10) !== 0;
    return /X|Y|Z/.test(text);
}

function coordSites(section: Cp2kSection): Site[] {
    const sites: Site[] = [];
    const scaledTokens = sectionKey(section, 'SCALED');
    // A bare `SCALED` is CP2K's shorthand for `.TRUE.`.
    const scaled = scaledTokens !== undefined &&
        (scaledTokens.length === 0 || /^\.?T/i.test(scaledTokens[0] || ''));
    const unitTokens = sectionKey(section, 'UNIT');
    const unitText = unitTokens ? (unitTokens[0] || '').toLowerCase() : 'angstrom';
    const factor = unitText.indexOf('bohr') >= 0 || unitText.indexOf('a.u') >= 0 || unitText === 'au'
        ? BOHR_TO_ANGSTROM
        : unitText.indexOf('nm') >= 0 ? 10
        : unitText.indexOf('pm') >= 0 ? 0.01
        : 1;

    for (const line of section.raw) {
        const parts = line.split(/\s+/);
        if (parts.length < 4) continue;
        const x = parseFloat(parts[1]);
        const y = parseFloat(parts[2]);
        const z = parseFloat(parts[3]);
        if (!isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
        const element = normalizeElementSymbol(parts[0]) || elementFromAtomName(parts[0]) || parts[0];
        if (scaled) sites.push({ element, frac: [x, y, z], label: parts[0] });
        else sites.push({ element, cart: [x * factor, y * factor, z * factor], label: parts[0] });
    }
    return sites;
}

/** The coordinate-file format the input states, else the file's extension. */
function cp2kCoordFormat(sections: Cp2kSection[], fileName: string): string {
    const topology = lastSection(sections, 'TOPOLOGY');
    const tokens = topology ? sectionKey(topology, 'COORD_FILE_FORMAT') : undefined;
    const declared = tokens ? tokens.join(' ').trim().toLowerCase() : '';
    // CP2K's default is OFF, which makes it detect the format from the name.
    if (declared !== '' && declared !== 'off') return declared;
    return fileName.toLowerCase().split('.').pop() || '';
}

/** Read a CP2K companion coordinate file in the formats an input can name. */
function cp2kCoordinateSites(content: string, format: string): Site[] {
    const kind = format.toLowerCase().replace(/^\./, '');
    if (kind === 'xyz') {
        const xyz = parseXyz(content);
        return xyz.atoms.map(atom => {
            const site: Site = { element: atom.element, cart: [atom.x, atom.y, atom.z] };
            if (atom.ghost) site.ghost = true;
            return site;
        });
    }
    if (kind === 'pdb' || kind === 'ent') {
        const pdb = parsePdb(content);
        return pdb.atoms.map(atom => ({ element: atom.element, cart: [atom.x, atom.y, atom.z] }));
    }
    if (kind === 'cp2k' || kind === 'crd' || kind === 'coord') {
        const sites: Site[] = [];
        for (const raw of content.split(/\r?\n/)) {
            const line = raw.replace(/#.*$/, '').trim();
            if (line === '') continue;
            const parts = line.split(/\s+/);
            // A leading atom count (the CP2K/CRD layout) is not a site row.
            if (parts.length < 4 || !/^[A-Za-z]/.test(parts[0])) continue;
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
            const element = normalizeElementSymbol(parts[0]) || elementFromAtomName(parts[0]) || 'X';
            sites.push({ element, cart: [x, y, z], label: parts[0] });
        }
        return sites;
    }
    return [];
}

export function parseCp2kInput(content: string, options?: Cp2kParseOptions): MolecularData {
    const warnings: string[] = [];
    const sections = flattenSections(parseCp2kSections(content));
    const title = (() => {
        const global = lastSection(sections, 'GLOBAL');
        const project = global && global.keys['PROJECT'];
        return oneLineTitle(project ? project.join(' ') : '', 'CP2K input');
    })();

    const coord = lastSection(sections, 'COORD');
    let sites = coord ? coordSites(coord) : [];
    const topology = lastSection(sections, 'TOPOLOGY');
    const coordName = topology && topology.keys['COORD_FILE_NAME']
        ? topology.keys['COORD_FILE_NAME'].join(' ').trim().replace(/^['"]|['"]$/g, '').trim()
        : '';
    const external = options && options.externalCoordinateFile ? options.externalCoordinateFile : undefined;
    if (coordName !== '' || external) {
        if (external) {
            const format = cp2kCoordFormat(sections, external.fileName);
            const externalSites = cp2kCoordinateSites(external.content, format);
            if (externalSites.length > 0) {
                // A stated COORD_FILE_NAME is where CP2K reads the geometry;
                // an inline &COORD block is then not what the run uses.
                sites = externalSites;
            } else if (['xyz', 'pdb', 'ent', 'cp2k', 'crd', 'coord'].indexOf(format) < 0) {
                warnings.push('COORD_FILE_NAME ' + (coordName || external.fileName) +
                    ': the coordinate format "' + (format || 'unknown') + '" is not read');
            } else {
                warnings.push('COORD_FILE_NAME ' + (coordName || external.fileName) +
                    ': the coordinate file has no readable atoms');
            }
        } else if (coordName !== '' && sites.length === 0) {
            warnings.push('the CP2K input reads its structure from "' + coordName +
                '", but that file could not be read');
        } else if (coordName !== '') {
            warnings.push('the CP2K input reads its structure from "' + coordName +
                '", but that file could not be read: showing the inline &COORD block');
        }
    }
    if (sites.length === 0 && coordName === '' && !external) {
        warnings.push('the CP2K input has no &COORD block: no atoms to show');
    }

    const lattice = cellLattice(sections);
    const periodic = lattice !== null && isPeriodic(sections);
    if (sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: warnings.length > 0 ? warnings : ['the CP2K input has no atoms to show']
        };
    }
    if (periodic && lattice && isUsableLattice(lattice)) {
        return buildPeriodicStructure(lattice, sites, title, warnings);
    }
    if (lattice && !periodic) {
        // PERIODIC NONE: the ABC box is a bounding box, not a cell. SCALED
        // coordinates are fractions of that box, so they still need converting
        // before the molecule is drawn without a lattice.
        const molecular = sites.map(s => s.frac
            ? { element: s.element, label: s.label, ghost: s.ghost, cart: fracToCart(s.frac, lattice) }
            : s);
        return buildMolecularStructure(molecular, title, warnings);
    }
    if (!lattice) warnings.push('the CP2K input states no &CELL: the structure is read as a molecule');
    else warnings.push('the CP2K &CELL is degenerate: the structure is read as a molecule');
    if (sites.some(s => s.frac)) {
        // Fractional coordinates are meaningless without the cell that scales
        // them; drawing the raw 0..1 numbers as Angstrom would stack every
        // atom at the origin instead of admitting the input is incomplete.
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: warnings.concat('&COORD is SCALED but the input states no usable &CELL: the fractional coordinates cannot be placed')
        };
    }
    return buildMolecularStructure(sites, title, warnings);
}

interface Cp2kOutGeometry {
    lattice: number[][] | null;
    sites: Site[];
    periodic: boolean;
    warnings: string[];
}

function parseCp2kOutputGeometry(content: string): Cp2kOutGeometry {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    let lattice: number[][] | null = null;
    let periodic = true;

    // A run over several steps repeats the cell; the last complete a/b/c
    // group is the final geometry. Requiring the group to arrive in order
    // keeps a truncated final print from mixing vectors of two cells.
    let current: number[][] = [];
    for (const line of lines) {
        const m = /CELL\|\s*Vector\s+([abc])\s*\[(angstrom|bohr)\]:\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)/i.exec(line);
        if (m) {
            const unit = m[2].toLowerCase();
            const f = unit === 'bohr' ? BOHR_TO_ANGSTROM : 1;
            const vector = [parseFloat(m[3]) * f, parseFloat(m[4]) * f, parseFloat(m[5]) * f];
            const letter = m[1].toLowerCase();
            if (letter === 'a') current = [vector];
            else if (letter === 'b' && current.length === 1) current.push(vector);
            else if (letter === 'c' && current.length === 2) {
                current.push(vector);
                lattice = current;
                current = [];
            } else current = [];
            continue;
        }
        const per = /CELL\|\s*Periodicity\s+(\S+)/i.exec(line);
        if (per && /NONE/i.test(per[1])) periodic = false;
    }

    // The last coordinate table is the geometry the run finished with.
    let sites: Site[] = [];
    for (let i = lines.length - 1; i >= 0 && sites.length === 0; i--) {
        const header = /ATOMIC COORDINATES IN\s+(ANGSTROM|A\.U\.|BOHR)/i.exec(lines[i]);
        if (!header) continue;
        const factor = /angstrom/i.test(header[1]) ? 1 : BOHR_TO_ANGSTROM;
        const rows: Site[] = [];
        for (let j = i + 1; j < lines.length; j++) {
            const line = lines[j].trim();
            if (line === '') {
                if (rows.length > 0) break;
                continue;
            }
            const parts = line.split(/\s+/);
            // The table CP2K prints is `Atom Kind Element Z X Y Z Z(eff) Mass`
            // with an unlabeled atomic-number column after the element; a
            // table without that column (some versions, older print sections)
            // is `Atom Kind Element X Y Z Z(eff) Mass`. Reading x from the
            // atomic number would put every atom at Z=8/14 of its real place.
            let xCol = 3;
            if (parts.length >= 9 && /^\d+$/.test(parts[3] || '')) xCol = 4;
            if (parts.length < xCol + 3) {
                if (rows.length > 0) break;
                continue;
            }
            const index = parseInt(parts[0], 10);
            const x = parseFloat(parts[xCol]);
            const y = parseFloat(parts[xCol + 1]);
            const z = parseFloat(parts[xCol + 2]);
            if (!isFinite(index) || !isFinite(x) || !isFinite(y) || !isFinite(z)) {
                if (rows.length > 0) break;
                continue;
            }
            const element = normalizeElementSymbol(parts[2]) || elementFromAtomName(parts[2]) || parts[2];
            rows.push({ element, cart: [x * factor, y * factor, z * factor], label: parts[2] });
        }
        if (rows.length > 0) sites = rows;
    }

    if (!lattice) warnings.push('the CP2K output states no cell vectors: the structure is read as a molecule');
    if (sites.length === 0) warnings.push('the CP2K output states no coordinate table: no structure to show');
    if (!periodic && lattice) {
        // A `PERIODIC NONE` run's cell is a box around the molecule; keeping it
        // would turn a water molecule into a 20 A crystal, so it is dropped.
        lattice = null;
    } else if (lattice && !isUsableLattice(lattice)) {
        warnings.push('the CP2K output cell vectors are degenerate: the structure is read as a molecule');
    }
    return { lattice, sites, periodic, warnings };
}

export function parseCp2kOutput(content: string): MolecularData {
    const geometry = parseCp2kOutputGeometry(content);
    const titleMatch = /GLOBAL\|\s*Project name\s+(\S+)/i.exec(content);
    const title = titleMatch ? 'CP2K ' + titleMatch[1] : 'CP2K output';
    if (geometry.sites.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: geometry.warnings
        };
    }
    if (geometry.lattice && geometry.periodic && isUsableLattice(geometry.lattice)) {
        return buildPeriodicStructure(geometry.lattice, geometry.sites, title, geometry.warnings);
    }
    return buildMolecularStructure(geometry.sites, title, geometry.warnings);
}

export function isCp2kOutput(content: string): boolean {
    return /CP2K\|?\s*version string/i.test(content) ||
        (/MODULE QUICKSTEP:?\s+ATOMIC COORDINATES/i.test(content) && /CP2K/i.test(content));
}

export function isCp2kInput(content: string): boolean {
    return /^\s*&(GLOBAL|FORCE_EVAL|SUBSYS|COORD|CELL)\b/m.test(content) &&
        /^\s*&END\b/m.test(content);
}

export { latticeToCellParams };
