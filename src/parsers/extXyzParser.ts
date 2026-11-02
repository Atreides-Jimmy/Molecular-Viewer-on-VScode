import { MolecularData } from '../types';
import { normalizeElementSymbol, elementFromAtomName, elementFromMass } from '../tools/periodicTable';
import { buildMolecularStructure, buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * Extended XYZ / GPUMD `model.xyz` reader.
 *
 * The comment line carries the cell (`Lattice="..."`, row-major) and the
 * column layout (`Properties=species:S:1:pos:R:3:...`). A plain XYZ parser
 * reads the first four columns and throws the cell away, which turns every
 * GPUMD or ASE/OVITO file into a molecule; this reader keeps the cell when the
 * line states one and the `pbc` flags allow it. The property layout is honoured
 * by name, so a file whose `pos` is not the first triple after `species`
 * (charge before position, force before position) still reads correctly.
 *
 * Only the first frame is shown: the viewer draws one structure, and a
 * periodic trajectory has no cell per frame in the frame list it navigates.
 * The count of frames beyond the first is reported instead of being lost.
 */

interface XyzProperty {
    name: string;
    type: string;
    count: number;
}

interface XyzComment {
    lattice: number[][] | null;
    properties: XyzProperty[] | null;
    pbc: boolean[] | null;
    frameTitle: string;
}

function readQuotedValue(comment: string, key: string): string | null {
    const re = new RegExp('(?:^|\\s)' + key + '\\s*=\\s*', 'i');
    const m = re.exec(comment);
    if (!m) return null;
    const start = m.index + m[0].length;
    const quote = comment.charAt(start);
    if (quote === '"' || quote === "'") {
        const end = comment.indexOf(quote, start + 1);
        return comment.slice(start + 1, end >= 0 ? end : comment.length);
    }
    // An unquoted vector (`Lattice=5.431 0 0 ...`) runs to the next key; taking
    // only the first token would drop the rest of the cell.
    const rest = comment.slice(start);
    const next = /\s+[A-Za-z_][A-Za-z0-9_]*\s*=/.exec(rest);
    return (next ? rest.slice(0, next.index) : rest).trim();
}

function parseProperties(comment: string): XyzProperty[] | null {
    const text = readQuotedValue(comment, 'Properties');
    if (!text) return null;
    const tokens = text.trim().split(':');
    if (tokens.length < 3 || tokens.length % 3 !== 0) return null;
    const out: XyzProperty[] = [];
    for (let i = 0; i + 2 < tokens.length; i += 3) {
        const count = parseInt(tokens[i + 2], 10);
        out.push({ name: tokens[i].toLowerCase(), type: tokens[i + 1], count: isFinite(count) && count > 0 ? count : 1 });
    }
    return out;
}

function parseLattice(comment: string): number[][] | null {
    const text = readQuotedValue(comment, 'Lattice');
    if (!text) return null;
    const values = text.trim().split(/\s+/).map(Number).filter(v => isFinite(v));
    if (values.length < 9) return null;
    const lattice = [
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]]
    ];
    return isUsableLattice(lattice) ? lattice : null;
}

function parsePbc(comment: string): boolean[] | null {
    const text = readQuotedValue(comment, 'pbc');
    if (!text) return null;
    const tokens = text.trim().split(/\s+/);
    if (tokens.length < 3) return null;
    return tokens.slice(0, 3).map(t => /^(t|true|1)$/i.test(t));
}

function parseComment(comment: string): XyzComment {
    return {
        lattice: parseLattice(comment),
        properties: parseProperties(comment),
        pbc: parsePbc(comment),
        frameTitle: comment.trim()
    };
}

export function isExtXyz(content: string): boolean {
    const lines = content.split(/\r?\n/);
    if (lines.length < 2) return false;
    if (!/^\s*\d+\s*$/.test(lines[0])) return false;
    const comment = lines[1];
    return /(?:^|\s)(?:Lattice|Properties|pbc)\s*=/i.test(comment);
}

function countFrames(lines: string[]): number {
    let frames = 0;
    let i = 0;
    while (i < lines.length) {
        const line = lines[i].trim();
        if (line === '') { i++; continue; }
        const count = parseInt(line, 10);
        if (!isFinite(count) || count < 0) break;
        frames++;
        i += 2 + count;
    }
    return frames;
}

export function parseExtXyz(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    if (lines.length < 2) {
        return {
            atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false,
            warnings: ['too few lines for an extended XYZ file: no structure to show']
        };
    }
    const count = parseInt(lines[0].trim(), 10);
    if (!isFinite(count) || count <= 0) {
        return {
            atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false,
            warnings: ['the first line of the extended XYZ file is not an atom count: no structure to show']
        };
    }
    const comment = parseComment(lines[1]);
    const properties = comment.properties;
    // Without a Properties declaration the file is a plain XYZ with extra
    // columns; species and pos are the first four columns, as usual.
    const layout: XyzProperty[] = properties || [
        { name: 'species', type: 'S', count: 1 },
        { name: 'pos', type: 'R', count: 3 }
    ];
    // Find species, position and mass columns in one pass over the layout.
    let posIndex = -1;
    let speciesCol = -1;
    let massCol = -1;
    let cursor = 0;
    for (const prop of layout) {
        if (prop.name === 'species' || prop.name === 'element' || prop.name === 'symbol') speciesCol = cursor;
        if (prop.name === 'pos' || prop.name === 'position' || prop.name === 'positions') posIndex = cursor;
        if (prop.name === 'mass') massCol = cursor;
        cursor += prop.count;
    }
    if (speciesCol < 0) {
        return {
            atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false,
            warnings: ['the Properties line names no species column: no structure to show']
        };
    }
    if (posIndex < 0) {
        return {
            atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false,
            warnings: ['the Properties line names no position column: no structure to show']
        };
    }

    const sites: Site[] = [];
    // The frame owns exactly `count` non-empty rows after the comment line.
    // Walking past that boundary when a row is unreadable would take the next
    // frame's atoms (its count and comment lines are just skipped) and show
    // them as frame one.
    let seen = 0;
    for (let i = 2; i < lines.length && seen < count; i++) {
        const line = lines[i].trim();
        if (line === '') continue;
        seen++;
        const parts = line.split(/\s+/);
        if (parts.length <= Math.max(speciesCol, posIndex + 2)) continue;
        const label = parts[speciesCol];
        let element = normalizeElementSymbol(label) || elementFromAtomName(label) || 'X';
        // A mass column lets a generic `X` label become the element the mass
        // names (GPUMD sometimes writes `X` with the mass beside it).
        if (element === 'X' && massCol >= 0) {
            const mass = parseFloat(parts[massCol]);
            element = elementFromMass(mass) || 'X';
        }
        const x = parseFloat(parts[posIndex]);
        const y = parseFloat(parts[posIndex + 1]);
        const z = parseFloat(parts[posIndex + 2]);
        if (!isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
        sites.push({ element, cart: [x, y, z], label });
    }
    if (sites.length === 0) {
        return {
            atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false,
            warnings: ['the extended XYZ file has no readable atom rows: no structure to show']
        };
    }

    const title = oneLineTitle(comment.frameTitle, 'Extended XYZ');
    const frames = countFrames(lines);
    if (frames > 1) warnings.push(frames + ' frames in the file: showing the first');

    const periodic = comment.lattice !== null &&
        (comment.pbc === null || comment.pbc.some(v => v));
    if (periodic && comment.lattice) {
        return buildPeriodicStructure(comment.lattice, sites, title, warnings);
    }
    if (comment.lattice && comment.pbc && !comment.pbc.some(v => v)) {
        warnings.push('pbc is F F F: the Lattice is a box around a molecule and is not drawn');
    }
    return buildMolecularStructure(sites, title, warnings);
}
