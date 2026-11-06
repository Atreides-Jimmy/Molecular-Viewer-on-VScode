import { MolecularData } from '../types';
import { elementFromZ, normalizeElementSymbol, typeElementLabel } from '../tools/periodicTable';
import { buildMolecularStructure, buildPeriodicStructure, isUsableLattice, oneLineTitle, Site } from './periodicUtils';

/**
 * DeepMD / DP-GEN dataset reader, for `coord.raw` and `coord.npy`.
 *
 * A dataset is a directory, not one file: `coord` holds the frames, `box` the
 * cell of each frame, `type.raw` the type of every atom (a single list shared
 * by all frames), `type_map.raw` the element each type stands for, and an
 * empty `nopbc` file marks a non-periodic system. The loader therefore reads
 * the companions and hands them in; this reader only turns them into a
 * structure and refuses to guess when a companion is missing:
 *
 *  - a coordinate row of 3N numbers is N atoms of three Cartesian components;
 *  - a box row of nine numbers is the row-major cell, and an all-zero box plus
 *    `nopbc` means the frame is a molecule in vacuum;
 *  - types are 0-based indexes into `type_map.raw`, so the first type is the
 *    first element named there, not `T0`.
 *
 * The viewer shows one structure, so the first frame is drawn and the number
 * of remaining frames is reported as a warning rather than silently dropped.
 */

export interface DeepmdParseOptions {
    /** Element per type index, from `type_map.raw`. */
    typeMap?: string[];
    /** Type index per atom, from `type.raw`. */
    types?: number[];
    /** The first frame's nine box numbers, from `box.raw` or `box.npy`. */
    box?: number[];
    /** The dataset has a `nopbc` marker: a molecule, even with a box. */
    nopbc?: boolean;
    /** Raw bytes when the primary file is a `.npy`. */
    binary?: Uint8Array;
}

export interface NpyArray {
    shape: number[];
    fortranOrder: boolean;
    dtype: string;
    /** Values in C (row-major) order. */
    values: number[];
}

function npyHeader(data: Uint8Array): { offset: number; shape: number[]; fortranOrder: boolean; dtype: string } | null {
    if (data.length < 10) return null;
    if (!(data[0] === 0x93 && data[1] === 0x4e && data[2] === 0x55 && data[3] === 0x4d && data[4] === 0x50 && data[5] === 0x59)) return null;
    const major = data[6];
    let headerStart = 10;
    let headerLength = data[8] | (data[9] << 8);
    if (major >= 2) {
        if (data.length < 12) return null;
        headerStart = 12;
        headerLength = (data[8] | (data[9] << 8) | (data[10] << 16) | (data[11] << 24)) >>> 0;
    }
    const end = headerStart + headerLength;
    if (end > data.length) return null;
    let text = '';
    for (let i = headerStart; i < end; i++) text += String.fromCharCode(data[i]);

    const descr = /['"]descr['"]\s*:\s*['"]([^'"]+)['"]/.exec(text);
    if (!descr) return null;
    const order = /['"]fortran_order['"]\s*:\s*(True|False)/.exec(text);
    const shapeMatch = /['"]shape['"]\s*:\s*\(([^)]*)\)/.exec(text);
    const shape: number[] = [];
    if (shapeMatch) {
        for (const part of shapeMatch[1].split(',')) {
            const v = parseInt(part.trim(), 10);
            if (isFinite(v)) shape.push(v);
        }
    }
    return {
        offset: end,
        shape,
        fortranOrder: !!order && order[1] === 'True',
        dtype: descr[1]
    };
}

/** Decode a `.npy` array; returns null when the bytes are not a numeric array. */
export function readNpy(data: Uint8Array, maxValues = Infinity): NpyArray | null {
    const header = npyHeader(data);
    if (!header) return null;
    const m = /^([<>=|])([ifu])(\d+)$/.exec(header.dtype);
    if (!m) return null;
    const littleEndian = m[1] !== '>';
    const kind = m[2];
    const bytes = parseInt(m[3], 10);
    const width: { [key: string]: number } = { i: 1, u: 1, f: 1 };
    if (!width[kind] || ![1, 2, 4, 8].includes(bytes)) return null;
    // Float16 has no DataView accessor here; returning values filled with NaN
    // would turn a readable dtype into an all-NaN structure.
    if (kind === 'f' && bytes !== 4 && bytes !== 8) return null;
    if (header.offset + bytes > data.length) return null;

    let total = 1;
    for (const n of header.shape) total *= n;
    if (total <= 0 || !isFinite(total)) {
        total = Math.floor((data.length - header.offset) / bytes);
    }
    const available = Math.floor((data.length - header.offset) / bytes);
    const count = Math.min(total, available, maxValues);
    const view = new DataView(data.buffer, data.byteOffset + header.offset, count * bytes);
    const values: number[] = [];
    for (let i = 0; i < count; i++) {
        const at = i * bytes;
        if (kind === 'f') {
            if (bytes === 4) values.push(view.getFloat32(at, littleEndian));
            else values.push(view.getFloat64(at, littleEndian));
        } else if (kind === 'u') {
            if (bytes === 1) values.push(view.getUint8(at));
            else if (bytes === 2) values.push(view.getUint16(at, littleEndian));
            else if (bytes === 4) values.push(view.getUint32(at, littleEndian));
            else {
                const lowAt = littleEndian ? at : at + 4;
                const highAt = littleEndian ? at + 4 : at;
                const low = view.getUint32(lowAt, littleEndian);
                const high = view.getUint32(highAt, littleEndian);
                values.push(high * 4294967296 + low);
            }
        } else {
            if (bytes === 1) values.push(view.getInt8(at));
            else if (bytes === 2) values.push(view.getInt16(at, littleEndian));
            else if (bytes === 4) values.push(view.getInt32(at, littleEndian));
            else if (bytes === 8) {
                // The word order follows the array's endianness, not the host's:
                // a big-endian i8 stores the high word first.
                const lowAt = littleEndian ? at : at + 4;
                const highAt = littleEndian ? at + 4 : at;
                const low = view.getUint32(lowAt, littleEndian);
                const high = view.getInt32(highAt, littleEndian);
                values.push(high * 4294967296 + low);
            } else values.push(NaN);
        }
    }

    // DeepMD writes C order, but an array saved by another tool may not.
    if (header.fortranOrder && header.shape.length === 2) {
        const rows = header.shape[0];
        const cols = header.shape[1];
        if (rows * cols <= values.length) {
            const reordered: number[] = new Array(rows * cols);
            for (let r = 0; r < rows; r++) {
                for (let c = 0; c < cols; c++) reordered[r * cols + c] = values[c * rows + r];
            }
            return { shape: header.shape, fortranOrder: true, dtype: header.dtype, values: reordered };
        }
    }
    return { shape: header.shape, fortranOrder: header.fortranOrder, dtype: header.dtype, values };
}

/** The first coordinate frame of a DeepMD file, and how many frames it holds. */
export interface DeepmdFrames {
    first: number[] | null;
    count: number;
}

/**
 * The first frame row of a DeepMD coordinate file, raw text or `.npy`.
 *
 * `natomsHint` is the length of `type.raw` (the atom count shared by every
 * frame). A two-dimensional `(natoms, 3)` array is one frame written per atom,
 * not `natoms` frames of a single atom; without the type list the two readings
 * are indistinguishable.
 */
export function deepmdCoordinateFrames(content: string, binary?: Uint8Array, natomsHint?: number): DeepmdFrames {
    if (binary && binary.length > 0) {
        const header = npyHeader(binary);
        if (header && header.shape.length > 0) {
            // The first axis is normally the frame index and everything after
            // it is one frame's 3N coordinates: `(frames, 3N)` and
            // `(frames, natoms, 3)` are both single frames per row. Only that
            // row is decoded for a C-order array; a Fortran-order file has to
            // be read in full to recover its first row.
            let frameWidth: number;
            let frameCount: number;
            if (header.shape.length >= 2) {
                const rowWidth = header.shape.slice(1).reduce((a, b) => a * b, 1);
                const oneFramePerAtom = header.shape.length === 2 && rowWidth === 3 &&
                    natomsHint !== undefined && natomsHint > 0 && header.shape[0] === natomsHint;
                frameWidth = oneFramePerAtom ? header.shape[0] * rowWidth : rowWidth;
                frameCount = oneFramePerAtom ? 1 : header.shape[0];
            } else {
                frameWidth = header.shape[0];
                frameCount = 1;
            }
            if (!isFinite(frameWidth) || frameWidth <= 0) return { first: null, count: 0 };
            const array = readNpy(binary, header.fortranOrder ? Infinity : frameWidth);
            if (array && array.values.length >= frameWidth) {
                return {
                    first: array.values.slice(0, frameWidth),
                    count: isFinite(frameCount) ? Math.max(1, Math.floor(frameCount)) : 1
                };
            }
        }
        return { first: null, count: 0 };
    }
    let first: number[] | null = null;
    let count = 0;
    for (const line of content.split(/\r?\n/)) {
        const text = line.trim();
        if (text === '') continue;
        const row = text.split(/\s+/).map(Number);
        if (row.length === 0 || row.some(v => !isFinite(v))) continue;
        count++;
        if (first === null) first = row;
    }
    return { first, count };
}

/** Parse the integer type list of `type.raw`. */
export function parseDeepmdTypes(content: string): number[] {
    const out: number[] = [];
    for (const line of content.split(/\r?\n/)) {
        for (const token of line.trim().split(/\s+/)) {
            if (token === '') continue;
            const v = parseInt(token, 10);
            if (isFinite(v)) out.push(v);
        }
    }
    return out;
}

/** Parse `type_map.raw` into element symbols, keeping the file's order. */
export function parseDeepmdTypeMap(content: string): string[] {
    const out: string[] = [];
    for (const line of content.split(/\r?\n/)) {
        const text = line.trim();
        if (text === '') continue;
        // A map entry may be a symbol or an atomic number.
        const z = parseInt(text, 10);
        const byZ = isFinite(z) ? elementFromZ(z) : undefined;
        out.push(normalizeElementSymbol(text) || byZ || text);
    }
    return out;
}

/** The first frame's box row of a `box.raw`/`box.npy` companion. */
export function deepmdFirstBox(content: string, binary?: Uint8Array): number[] | undefined {
    if (binary && binary.length > 0) {
        const array = readNpy(binary);
        if (array && array.values.length >= 9) return array.values.slice(0, 9);
    }
    for (const line of content.split(/\r?\n/)) {
        const text = line.trim();
        if (text === '') continue;
        const row = text.split(/\s+/).map(Number);
        if (row.length >= 9 && row.slice(0, 9).every(v => isFinite(v))) return row.slice(0, 9);
    }
    return undefined;
}

export function isDeepmdDatasetFile(fileName: string): boolean {
    return /^coord\.(raw|npy)$/i.test(fileName);
}

/**
 * DP-GEN's `model_devi.out`: the deviation table `dpgen run` writes during
 * model deviation. It shares `.out` with the electronic-structure programs but
 * holds no coordinates, so the viewer identifies it and points at the frames.
 */
export function isDpgenModelDevi(content: string): boolean {
    return /^\s*#?\s*step\s+max_devi_v\s+min_devi_v\s+avg_devi_v\s+max_devi_f\s+min_devi_f\s+avg_devi_f\b/mi.test(content);
}

function boxLattice(values: number[] | undefined): number[][] | null {
    if (!values || values.length < 9) return null;
    const lattice = [
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]]
    ];
    return isUsableLattice(lattice) ? lattice : null;
}

export function parseDeepmdDataset(content: string, fileName: string, options?: DeepmdParseOptions): MolecularData {
    const warnings: string[] = [];
    const title = oneLineTitle('DeepMD ' + fileName, 'DeepMD dataset');
    const dataset = deepmdCoordinateFrames(
        content,
        options && options.binary,
        options && options.types && options.types.length > 0 ? options.types.length : undefined
    );
    if (!dataset.first || dataset.first.length === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['no coordinate rows in the DeepMD dataset: no structure to show']
        };
    }
    const frame = dataset.first;
    if (!frame.every(v => isFinite(v))) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the first DeepMD coordinate frame holds non-finite numbers: no structure to show']
        };
    }
    const width = frame.length;
    if (width === 0 || width % 3 !== 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['a DeepMD coordinate row holds ' + width + ' numbers, which is not three per atom: no structure to show']
        };
    }

    const types = options && options.types ? options.types : [];
    const typeMap = options && options.typeMap ? options.typeMap : [];
    let natom = width / 3;
    if (types.length > 0 && types.length !== natom) {
        warnings.push('type.raw lists ' + types.length + ' atoms but the first frame holds ' + natom + ': the shorter list decides the elements');
        natom = Math.min(natom, types.length);
    }

    const sites: Site[] = [];
    for (let i = 0; i < natom; i++) {
        const type = i < types.length ? types[i] : 0;
        const element = typeMap[type] || typeElementLabel(type + 1);
        sites.push({ element, cart: [frame[3 * i], frame[3 * i + 1], frame[3 * i + 2]] });
    }
    if (typeMap.length === 0) {
        warnings.push('no type_map.raw next to the dataset: atoms are named by their type index');
    }

    if (dataset.count > 1) warnings.push(dataset.count + ' frames in the dataset: showing the first');

    const lattice = boxLattice(options && options.box);
    const nonPeriodic = !!(options && options.nopbc);
    if (nonPeriodic) {
        if (lattice) warnings.push('the dataset has a nopbc marker: the box is not drawn');
        return buildMolecularStructure(sites, title, warnings);
    }
    if (!lattice) {
        warnings.push('no usable box for the dataset: the frame is read as a molecule');
        return buildMolecularStructure(sites, title, warnings);
    }
    return buildPeriodicStructure(lattice, sites, title, warnings);
}
