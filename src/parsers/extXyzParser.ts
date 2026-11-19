import { MolecularData, StructureStep } from '../types';
import { normalizeElementSymbol, elementFromAtomName, elementFromMass } from '../tools/periodicTable';
import {
    buildMolecularStructure, buildPeriodicStructure, isUsableLattice, oneLineTitle, structureSteps,
    StepGeometry, Site
} from './periodicUtils';

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
    /** `energy=` — the frame's total energy, kept so it is not dropped silently. */
    energy?: number;
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
    const tokens = text.trim().split(/\s+/).filter(v => v !== '');
    // The nine numbers are the three cell vectors in order. Filtering the
    // unreadable ones out instead of refusing them shifted every number after
    // the bad slot, so a cell with one typo became a different, plausible
    // lattice. A tenth value (some writers append the cell's volume) is
    // ignored; anything else that is not a number rejects the cell.
    if (tokens.length < 9) return null;
    const values: number[] = [];
    for (let i = 0; i < 9; i++) {
        const v = Number(tokens[i]);
        if (!isFinite(v)) return null;
        values.push(v);
    }
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
    const parsed: XyzComment = {
        lattice: parseLattice(comment),
        properties: parseProperties(comment),
        pbc: parsePbc(comment),
        frameTitle: comment.trim()
    };
    const energy = readQuotedValue(comment, 'energy');
    if (energy !== null) {
        const value = Number(energy);
        if (isFinite(value)) parsed.energy = value;
    }
    return parsed;
}

/** The comment-line keys that make a plain XYZ an *extended* one. */
const EXTXYZ_KEYS = /(?:^|\s)(?:Lattice|Properties|pbc)\s*=/i;

export function isExtXyz(content: string): boolean {
    // A byte-order mark and leading blank lines are not part of the format but
    // are written by enough editors that reading only the first two lines
    // turned a valid extended XYZ into a plain one (whose cell is then thrown
    // away, because a plain XYZ has no comment keys at all).
    const lines = content.replace(/^\uFEFF/, '').split(/\r?\n/);
    let i = 0;
    while (i < lines.length && lines[i].trim() === '') i++;
    if (i + 1 >= lines.length) return false;
    if (!/^\d+$/.test(lines[i].trim())) return false;
    return EXTXYZ_KEYS.test(lines[i + 1]);
}

/** One frame of the file: its comment, its atoms, and what did not read. */
interface ExtXyzFrame {
    comment: XyzComment;
    sites: Site[];
    /** Atom rows the frame's count line declares. */
    declared: number;
    /** Rows that were present but unreadable. */
    badRows: number;
    /** 1-based line number of the first unreadable row (0 when there is none). */
    firstBadLine: number;
    /** The file ended before the frame's declared rows did. */
    short: boolean;
}

/**
 * Split the file into frames and read each one's atoms.
 *
 * A frame is its count line, exactly one comment line, and that many non-empty
 * atom rows; nothing else delimits it. A row that cannot be read still consumes
 * its slot — walking past it towards the next line that parses is what let one
 * frame's atoms be shown as another's, and dropped the row without a word — so
 * the count of unreadable rows is reported instead.
 */
function readExtXyzFrames(lines: string[], start: number, speciesCol: number, posIndex: number, massCol: number): ExtXyzFrame[] {
    const frames: ExtXyzFrame[] = [];
    let i = start;
    while (i < lines.length) {
        while (i < lines.length && lines[i].trim() === '') i++;
        if (i >= lines.length) break;
        const countText = lines[i].trim();
        if (!/^\d+$/.test(countText)) break;
        const declared = parseInt(countText, 10);
        if (!isFinite(declared) || declared <= 0) break;
        i++;

        const comment = parseComment(i < lines.length ? lines[i] : '');
        i++;

        const sites: Site[] = [];
        let badRows = 0;
        let firstBadLine = 0;
        while (sites.length + badRows < declared) {
            if (i >= lines.length) break;
            const line = lines[i].trim();
            const lineNo = i + 1;
            i++;
            if (line === '') continue;
            const parts = line.split(/\s+/);
            if (parts.length <= Math.max(speciesCol, posIndex + 2)) {
                badRows++;
                if (!firstBadLine) firstBadLine = lineNo;
                continue;
            }
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
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) {
                badRows++;
                if (!firstBadLine) firstBadLine = lineNo;
                continue;
            }
            sites.push({ element, cart: [x, y, z], label });
        }
        const short = sites.length + badRows < declared;
        frames.push({ comment, sites, declared, badRows, firstBadLine, short });
        // A frame the file cut short ends the file: the rows after it cannot be
        // attributed without its missing atoms.
        if (short) break;
    }
    return frames;
}

/**
 * A GPUMD `run.in` control script. GPUMD is configured through this file,
 * which names the potential and the files to dump but holds no coordinates;
 * recognising it lets the viewer point at the `model.xyz` or `dump.xyz` it
 * writes instead of showing an unexplained empty canvas.
 */
export function isGpumdInput(content: string): boolean {
    // CP2K inputs use the same words (`POTENTIAL`, `ENSEMBLE`) inside their
    // sections; the section marker identifies them first.
    if (/^\s*&[A-Za-z]/m.test(content)) return false;
    let keywords = 0;
    for (const raw of content.split(/\r?\n/)) {
        const line = raw.replace(/#.*$/, '').trim();
        if (line === '' || /^\S+(?:\s+\S+){0,3}\s*=/.test(line)) continue;
        // GPUMD's commands are lower-case; keeping the case keeps CP2K's
        // upper-case keywords from matching here. `velocity` and `run` are
        // shared with LAMMPS, so they never count.
        const command = line.split(/\s+/)[0];
        if (/^(potential|ensemble|time_step|dump_thermo|dump_xyz|dump_exyz|dump_restart|dump_combined)$/.test(command)) {
            if (++keywords >= 2) return true;
        }
    }
    return false;
}

export function parseExtXyz(content: string): MolecularData {
    const lines = content.replace(/^\uFEFF/, '').split(/\r?\n/);
    const warnings: string[] = [];
    const nothing = (message: string): MolecularData => ({
        atoms: [], bonds: [], title: 'Extended XYZ', hasExplicitBonds: false, warnings: [message]
    });
    if (lines.length < 2) return nothing('too few lines for an extended XYZ file: no structure to show');

    // The count line may sit behind blank lines (and behind the byte-order mark
    // removed above); the comment line is the one right after it.
    let first = 0;
    while (first < lines.length && lines[first].trim() === '') first++;
    if (first >= lines.length) return nothing('the extended XYZ file is empty: no structure to show');
    const count = parseInt(lines[first].trim(), 10);
    if (!isFinite(count) || count <= 0) {
        return nothing('the first line of the extended XYZ file is not an atom count: no structure to show');
    }

    const comment = parseComment(first + 1 < lines.length ? lines[first + 1] : '');
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
        return nothing('the Properties line names no species column: no structure to show');
    }
    if (posIndex < 0) {
        return nothing('the Properties line names no position column: no structure to show');
    }

    const frames = readExtXyzFrames(lines, first, speciesCol, posIndex, massCol);
    const firstFrame = frames.length > 0 ? frames[0] : null;
    const sites = firstFrame ? firstFrame.sites : [];
    if (!firstFrame || sites.length === 0) {
        return nothing('the extended XYZ file has no readable atom rows: no structure to show');
    }
    // A row that did not read is a statement about the file, not something to
    // drop in silence: the declared count is what the frame promised, and the
    // structure that opens now holds fewer atoms than that.
    if (firstFrame.badRows > 0) {
        warnings.push(firstFrame.badRows + ' atom row' + (firstFrame.badRows === 1 ? '' : 's') +
            ' of the first frame could not be read (first at line ' + firstFrame.firstBadLine +
            '): ' + (firstFrame.badRows === 1 ? 'it is' : 'they are') + ' skipped');
    }
    if (firstFrame.short) {
        warnings.push('the file ends after ' + (sites.length + firstFrame.badRows) + ' of the ' +
            firstFrame.declared + ' atom rows the first frame declares');
    }
    let laterBad = 0;
    for (let k = 1; k < frames.length; k++) laterBad += frames[k].badRows;
    if (laterBad > 0) {
        warnings.push(laterBad + ' atom row' + (laterBad === 1 ? '' : 's') +
            ' in the frames after the first could not be read: ' + (laterBad === 1 ? 'it is' : 'they are') + ' skipped');
    }

    const title = oneLineTitle(comment.frameTitle, 'Extended XYZ');

    // A cell is drawn only when every direction is periodic, because the
    // viewer's crystal model is three-dimensional: it wraps atoms into the
    // cell and searches periodic images in all three vectors. Treating a
    // partially periodic cell (`pbc="T F F"`, a slab or a polymer) as a full
    // one invented bonds and a density in the directions the file said were
    // not periodic. Such a file is shown as the molecule it is, and the
    // warning names the directions that were dropped.
    const pbc = comment.pbc;
    const periodic = comment.lattice !== null && (pbc === null || pbc.every(v => v));

    // Every geometry after the first is a step of the same system: an ASE or
    // OVITO file that states one is a trajectory, and its frames used to be
    // counted and thrown away. `structureSteps` drops a geometry that repeats
    // the one before it and samples a long run, exactly as the MD readers do.
    let steps: StructureStep[] = [];
    let stepsInFile = 0;
    if (frames.length > 1) {
        const natoms = sites.length;
        const geometries: StepGeometry[] = [];
        let mismatched = 0;
        for (let k = 0; k < frames.length; k++) {
            const frame = frames[k];
            if (frame.sites.length !== natoms) {
                mismatched++;
                continue;
            }
            let cell: number[][] | undefined;
            if (periodic) {
                const own = frame.comment.lattice;
                const framePbc = frame.comment.pbc || pbc;
                const framePeriodic = own ? (framePbc === null || framePbc.every(v => v)) : true;
                cell = framePeriodic ? (own || comment.lattice || undefined) : undefined;
            }
            geometries.push({ label: 'Frame ' + (k + 1), sites: frame.sites, lattice: cell });
        }
        const collected = structureSteps(geometries, natoms);
        steps = collected.steps;
        stepsInFile = collected.total;
        if (mismatched > 0) {
            // A frame list has one element per frame: a geometry with another
            // atom count cannot be a frame of this structure, and saying so is
            // better than letting it disappear from the count.
            warnings.push(mismatched + ' frame' + (mismatched === 1 ? '' : 's') + ' state' +
                (mismatched === 1 ? 's' : '') + ' an atom count other than ' + natoms +
                ': only the matching frames are navigable');
        }
    }
    if (frames.length > 1 && steps.length <= 1) {
        warnings.push(frames.length + ' frames in the file: showing the first');
    }
    // The columns and keys a frame states beside its geometry are not drawn, so
    // the panel says they are there instead of letting them vanish.
    const energies = frames.filter(frame => frame.comment.energy !== undefined).length;
    const statesForces = layout.some(prop => /^(force|forces|f)$/i.test(prop.name));
    if (energies > 0 || statesForces) {
        const stated = energies > 0 && statesForces
            ? 'a total energy and per-atom forces'
            : energies > 0 ? 'a total energy' : 'per-atom forces';
        warnings.push('the file states ' + stated + ' beside the coordinates: the viewer draws the geometry, ' +
            'not ' + (energies > 0 && statesForces ? 'those columns' : 'that column'));
    }

    let data: MolecularData;
    if (periodic && comment.lattice) {
        data = buildPeriodicStructure(comment.lattice, sites, title, warnings);
    } else {
        if (comment.lattice && pbc && !pbc.some(v => v)) {
            warnings.push('pbc is F F F: the Lattice is a box around a molecule and is not drawn');
        } else if (comment.lattice && pbc) {
            const axes = ['x', 'y', 'z'].filter((_, i) => pbc[i]).join(', ');
            warnings.push('pbc is ' + pbc.map(v => (v ? 'T' : 'F')).join(' ') +
                ': only the ' + axes + ' direction' + (axes.indexOf(',') >= 0 ? 's are' : ' is') +
                ' periodic, and this viewer draws a full 3D cell. The cell is not applied and the ' +
                'atoms are shown as written, without periodic bonds');
        }
        data = buildMolecularStructure(sites, title, warnings);
    }
    if (steps.length > 1) {
        data.steps = steps;
        data.stepsInFile = stepsInFile;
    }
    return data;
}

