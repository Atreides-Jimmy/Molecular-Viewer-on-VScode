import { cellToLatticeVectors } from '../tools/crystalGeometry';

/**
 * Binary MD trajectories, read without a native dependency.
 *
 * Three formats cover the files a molecular-dynamics run leaves beside its
 * structure: GROMACS `.xtc` (compressed positions) and `.trr` (full precision
 * positions, velocities and forces), and the CHARMM/NAMD/LAMMPS `.dcd`
 * (Fortran-record positions with an optional unit cell). They are all raw
 * coordinate streams: a trajectory never says which element an atom is and
 * never states bonds, so the caller pairs it with a topology file it reads
 * through the ordinary structure readers.
 *
 * The reader is deliberately two-pass. The first pass walks the record
 * headers and records where each frame starts and how big it is; only then is
 * a frame selected for decoding. A run with a hundred thousand saved frames
 * therefore costs a header scan, not a full decode, and a truncated last
 * record costs that record rather than the whole file.
 */

export type TrajectoryFormat = 'xtc' | 'trr' | 'dcd';

export interface TrajectoryFrame {
    /** `[x0, y0, z0, x1, …]`, in Angstrom. */
    positions: Float64Array;
    /** The frame's cell vectors (rows a, b, c), in Angstrom, when it states one. */
    cell?: number[][];
    /** The time the file states, in the format's own unit (ps for XTC/TRR). */
    time?: number;
    /** The integration step the file states. */
    step?: number;
}

export interface TrajectoryReadOptions {
    /** Keep at most this many frames (the rest are sampled evenly). */
    maxFrames?: number;
    /** Keep at most this many atom coordinates across all kept frames. */
    maxTotalAtoms?: number;
    /**
     * Angstrom per stored length unit. XTC and TRR store nanometres, so the
     * default is 10; DCD stores the program's own unit, which is Angstrom for
     * CHARMM, NAMD and LAMMPS `real`/`metal`, so its default is 1. LAMMPS
     * `units si`/`electron`/… pass their own factor here.
     */
    lengthScale?: number;
}

export interface BinaryTrajectory {
    format: TrajectoryFormat;
    /** Atom count the file states. */
    natoms: number;
    frames: TrajectoryFrame[];
    /** Frames the file contains, before sampling. */
    frameCount: number;
    /** Whether some frames were skipped to stay inside the limits. */
    sampled: boolean;
    /** Whether the file ended in the middle of a record. */
    truncated: boolean;
    /** Angstrom per stored length unit that was applied. */
    lengthScale: number;
    warnings: string[];
}

const DEFAULT_MAX_FRAMES = 500;
const DEFAULT_MAX_TOTAL_ATOMS = 1_500_000;
/**
 * Above this many atoms in one frame the trajectory is reported as too large:
 * the page has to hold every coordinate of a frame in its own heap and the
 * browser is the limiting factor, not this scan.
 */
const MAX_SINGLE_FRAME_ATOMS = 5_000_000;

class RecordError extends Error {}

function ensure(dv: DataView, offset: number, size: number): void {
    if (!isFinite(offset) || offset < 0 || size < 0 || offset + size > dv.byteLength) {
        throw new RecordError('record extends past the end of the file');
    }
}

function finitePositive(v: number): boolean {
    return isFinite(v) && v > 0;
}

function usableCell(lattice: number[][] | undefined): number[][] | undefined {
    if (!lattice || lattice.length !== 3) return undefined;
    for (const row of lattice) {
        if (row.length !== 3 || row.some(v => !isFinite(v))) return undefined;
    }
    const [a, b, c] = lattice;
    const det =
        a[0] * (b[1] * c[2] - b[2] * c[1]) -
        a[1] * (b[0] * c[2] - b[2] * c[0]) +
        a[2] * (b[0] * c[1] - b[1] * c[0]);
    if (!isFinite(det) || Math.abs(det) < 1e-9) return undefined;
    return lattice;
}

function scaleCell(lattice: number[][], scale: number): number[][] {
    return lattice.map(row => row.map(v => v * scale));
}

function allFinite(values: Float64Array): boolean {
    for (let i = 0; i < values.length; i++) {
        if (!isFinite(values[i])) return false;
    }
    return true;
}

/* ------------------------------------------------------------------------- *
 * XTC
 *
 * Magic 1995, big-endian XDR. A frame is a fixed header (step, time, 3x3 box),
 * a coordinate count, and either that many raw float triples or a
 * self-contained compressed block (precision, integer bounds, small-integer
 * table index, byte count, opaque bytes). The compression is GROMACS's integer
 * scheme; the table and the decoder below are its documented constants.
 * ------------------------------------------------------------------------- */

const XTC_MAGIC = 1995;
const TRR_MAGIC = 1993;

const MAGIC_INTS = new Uint32Array([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216
]);
const XTC_FIRST_IDX = 9;

function sizeOfInt(size: number): number {
    let num = 1;
    let bits = 0;
    while (size >= num && bits < 32) {
        bits++;
        num *= 2;
    }
    return bits;
}

function sizeOfInts(sizes: number[]): number {
    const bytes = new Uint8Array(32);
    let numOfBytes = 1;
    bytes[0] = 1;
    for (let i = 0; i < sizes.length; i++) {
        let tmp = 0;
        let bytecnt = 0;
        for (; bytecnt < numOfBytes; bytecnt++) {
            tmp += bytes[bytecnt] * sizes[i];
            bytes[bytecnt] = tmp & 0xff;
            tmp = Math.floor(tmp / 256);
        }
        while (tmp !== 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp = Math.floor(tmp / 256);
        }
        numOfBytes = bytecnt;
    }
    let num = 1;
    numOfBytes--;
    let numOfBits = 0;
    while (bytes[numOfBytes] >= num) {
        numOfBits++;
        num *= 2;
    }
    return numOfBits + numOfBytes * 8;
}

function decodeBits(buf: Int32Array, cbuf: Uint8Array, numOfBits: number, buf2: Uint32Array): number {
    const mask = numOfBits >= 32 ? 0xffffffff : (1 << numOfBits) - 1;
    let lastBB0 = buf2[1];
    let lastBB1 = buf2[2];
    let cnt = buf[0];
    let num = 0;
    let bits = numOfBits;

    while (bits >= 8) {
        // A corrupt byte count must fail the frame, not read `undefined`
        // (which bitwise arithmetic turns into zero) and yield coordinates
        // from a stream that ended early.
        if (cnt >= cbuf.length) throw new RecordError('compressed coordinate block is truncated');
        lastBB1 = (lastBB1 << 8) | cbuf[cnt++];
        num |= (lastBB1 >> lastBB0) << (bits - 8);
        bits -= 8;
    }
    if (bits > 0) {
        if (lastBB0 < bits) {
            if (cnt >= cbuf.length) throw new RecordError('compressed coordinate block is truncated');
            lastBB0 += 8;
            lastBB1 = (lastBB1 << 8) | cbuf[cnt++];
        }
        lastBB0 -= bits;
        num |= (lastBB1 >> lastBB0) & ((1 << bits) - 1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastBB0;
    buf[2] = lastBB1;
    return num;
}

function decodeInts(
    buf: Int32Array,
    cbuf: Uint8Array,
    numOfInts: number,
    numOfBits: number,
    sizes: Uint32Array,
    nums: Int32Array,
    buf2: Uint32Array
): void {
    // Three small integers multiply out to at most 16777216^3, nine bytes.
    const bytes = new Int32Array(32);
    let numOfBytes = 0;
    let bits = numOfBits;
    while (bits > 8) {
        bytes[numOfBytes++] = decodeBits(buf, cbuf, 8, buf2);
        bits -= 8;
    }
    if (bits > 0) bytes[numOfBytes++] = decodeBits(buf, cbuf, bits, buf2);

    for (let i = numOfInts - 1; i > 0; i--) {
        let num = 0;
        for (let j = numOfBytes - 1; j >= 0; j--) {
            num = num * 256 + bytes[j];
            const p = Math.floor(num / sizes[i]);
            bytes[j] = p;
            num -= p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] + bytes[1] * 256 + bytes[2] * 65536 + bytes[3] * 16777216;
}

interface XtcFrameInfo {
    offset: number;
    natoms: number;
    step: number;
    time: number;
    box: number[];
    lsize: number;
    raw: boolean;
    dataAt: number;
    precision?: number;
    minInt?: number[];
    maxInt?: number[];
    smallIdx?: number;
    compressedBytes?: number;
}

function scanXtc(dv: DataView): { infos: XtcFrameInfo[]; truncated: boolean } {
    const infos: XtcFrameInfo[] = [];
    let offset = 0;
    let truncated = false;

    while (offset < dv.byteLength) {
        try {
            ensure(dv, offset, 12);
            const magic = dv.getInt32(offset, false);
            if (magic !== XTC_MAGIC) throw new RecordError('not an XTC frame');
            const natoms = dv.getInt32(offset + 4, false);
            const step = dv.getInt32(offset + 8, false);
            if (natoms <= 0 || natoms > MAX_SINGLE_FRAME_ATOMS) throw new RecordError('bad atom count');
            ensure(dv, offset + 12, 4 + 36 + 4);
            const time = dv.getFloat32(offset + 12, false);
            const box: number[] = [];
            for (let i = 0; i < 9; i++) box.push(dv.getFloat32(offset + 16 + i * 4, false));
            let p = offset + 16 + 36;
            const lsize = dv.getInt32(p, false);
            p += 4;
            if (lsize <= 0 || lsize > MAX_SINGLE_FRAME_ATOMS) throw new RecordError('bad coordinate count');

            if (lsize <= 9) {
                ensure(dv, p, lsize * 12);
                infos.push({
                    offset, natoms, step, time, box, lsize,
                    raw: true, dataAt: p
                });
                p += lsize * 12;
            } else {
                ensure(dv, p, 4 + 24 + 4 + 4);
                const precision = dv.getFloat32(p, false);
                const minInt = [
                    dv.getInt32(p + 4, false),
                    dv.getInt32(p + 8, false),
                    dv.getInt32(p + 12, false)
                ];
                const maxInt = [
                    dv.getInt32(p + 16, false),
                    dv.getInt32(p + 20, false),
                    dv.getInt32(p + 24, false)
                ];
                const smallIdx = dv.getInt32(p + 28, false);
                const compressedBytes = dv.getInt32(p + 32, false);
                if (!(precision > 0) || compressedBytes < 0) throw new RecordError('bad compressed block');
                const dataAt = p + 36;
                const padded = (compressedBytes + 3) & ~3;
                ensure(dv, dataAt, padded);
                infos.push({
                    offset, natoms, step, time, box, lsize,
                    raw: false, dataAt, precision, minInt, maxInt, smallIdx, compressedBytes
                });
                p = dataAt + padded;
            }
            offset = p;
        } catch {
            truncated = offset < dv.byteLength;
            break;
        }
    }
    return { infos, truncated };
}

function decodeXtcFrame(bytes: Uint8Array, dv: DataView, info: XtcFrameInfo): Float64Array {
    const natoms = info.natoms;
    if (info.raw) {
        const available = Math.min(info.lsize, natoms);
        const out = new Float64Array(natoms * 3);
        for (let i = 0; i < available * 3; i++) {
            out[i] = dv.getFloat32(info.dataAt + i * 4, false);
        }
        if (available < natoms) throw new RecordError('short raw coordinate block');
        return out;
    }

    const lsize = info.lsize;
    if (lsize < natoms) throw new RecordError('short compressed coordinate block');

    const sizeInt = new Uint32Array(3);
    sizeInt[0] = (info.maxInt![0] - info.minInt![0] + 1) >>> 0;
    sizeInt[1] = (info.maxInt![1] - info.minInt![1] + 1) >>> 0;
    sizeInt[2] = (info.maxInt![2] - info.minInt![2] + 1) >>> 0;
    if (sizeInt[0] === 0 || sizeInt[1] === 0 || sizeInt[2] === 0) {
        throw new RecordError('degenerate coordinate range');
    }

    const bitsizeInt = new Uint32Array(3);
    let bitsize = 0;
    if (sizeInt[0] > 0xffffff || sizeInt[1] > 0xffffff || sizeInt[2] > 0xffffff) {
        bitsizeInt[0] = sizeOfInt(sizeInt[0]);
        bitsizeInt[1] = sizeOfInt(sizeInt[1]);
        bitsizeInt[2] = sizeOfInt(sizeInt[2]);
    } else {
        bitsize = sizeOfInts([sizeInt[0], sizeInt[1], sizeInt[2]]);
    }

    let smallIdx = info.smallIdx!;
    if (smallIdx < 0 || smallIdx >= MAGIC_INTS.length) throw new RecordError('bad small-integer index');
    let tmpIdx = smallIdx - 1;
    tmpIdx = XTC_FIRST_IDX > tmpIdx ? XTC_FIRST_IDX : tmpIdx;
    let smaller = Math.floor(MAGIC_INTS[tmpIdx] / 2);
    let smallnum = Math.floor(MAGIC_INTS[smallIdx] / 2);
    const sizeSmall = new Uint32Array(3);
    sizeSmall[0] = sizeSmall[1] = sizeSmall[2] = MAGIC_INTS[smallIdx];

    const cbuf = new Uint8Array(bytes.buffer, bytes.byteOffset + info.dataAt, info.compressedBytes!);
    const buf = new Int32Array(3);
    const buf2 = new Uint32Array(buf.buffer);
    const thiscoord = new Int32Array(3);
    const prevcoord = new Int32Array(3);
    const out = new Float64Array(natoms * 3);
    const invPrecision = 1 / info.precision!;
    let lfp = 0;
    let run = 0;
    let i = 0;

    while (i < lsize) {
        if (bitsize === 0) {
            thiscoord[0] = decodeBits(buf, cbuf, bitsizeInt[0], buf2);
            thiscoord[1] = decodeBits(buf, cbuf, bitsizeInt[1], buf2);
            thiscoord[2] = decodeBits(buf, cbuf, bitsizeInt[2], buf2);
        } else {
            decodeInts(buf, cbuf, 3, bitsize, sizeInt, thiscoord, buf2);
        }
        i++;
        thiscoord[0] += info.minInt![0];
        thiscoord[1] += info.minInt![1];
        thiscoord[2] += info.minInt![2];
        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];

        const flag = decodeBits(buf, cbuf, 1, buf2);
        let isSmaller = 0;
        if (flag === 1) {
            run = decodeBits(buf, cbuf, 5, buf2);
            isSmaller = run % 3;
            run -= isSmaller;
            isSmaller--;
        }

        if (run > 0) {
            thiscoord[0] = thiscoord[1] = thiscoord[2] = 0;
            for (let k = 0; k < run; k += 3) {
                decodeInts(buf, cbuf, 3, smallIdx, sizeSmall, thiscoord, buf2);
                i++;
                thiscoord[0] += prevcoord[0] - smallnum;
                thiscoord[1] += prevcoord[1] - smallnum;
                thiscoord[2] += prevcoord[2] - smallnum;
                if (k === 0) {
                    // GROMACS swaps the first two atoms of a run: water's
                    // hydrogens then differ by a smaller integer.
                    let tmp = thiscoord[0]; thiscoord[0] = prevcoord[0]; prevcoord[0] = tmp;
                    tmp = thiscoord[1]; thiscoord[1] = prevcoord[1]; prevcoord[1] = tmp;
                    tmp = thiscoord[2]; thiscoord[2] = prevcoord[2]; prevcoord[2] = tmp;
                    if (lfp + 3 > out.length) throw new RecordError('coordinate overrun');
                    out[lfp++] = prevcoord[0] * invPrecision;
                    out[lfp++] = prevcoord[1] * invPrecision;
                    out[lfp++] = prevcoord[2] * invPrecision;
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                if (lfp + 3 > out.length) throw new RecordError('coordinate overrun');
                out[lfp++] = thiscoord[0] * invPrecision;
                out[lfp++] = thiscoord[1] * invPrecision;
                out[lfp++] = thiscoord[2] * invPrecision;
            }
        } else {
            if (lfp + 3 > out.length) throw new RecordError('coordinate overrun');
            out[lfp++] = thiscoord[0] * invPrecision;
            out[lfp++] = thiscoord[1] * invPrecision;
            out[lfp++] = thiscoord[2] * invPrecision;
        }

        smallIdx += isSmaller;
        if (isSmaller < 0) {
            smallnum = smaller;
            smaller = smallIdx > XTC_FIRST_IDX ? Math.floor(MAGIC_INTS[smallIdx - 1] / 2) : 0;
        } else if (isSmaller > 0) {
            smaller = smallnum;
            smallnum = Math.floor(MAGIC_INTS[smallIdx] / 2);
        }
        if (smallIdx < 0 || smallIdx >= MAGIC_INTS.length) throw new RecordError('small-integer index out of range');
        sizeSmall[0] = sizeSmall[1] = sizeSmall[2] = MAGIC_INTS[smallIdx];
    }
    if (lfp !== out.length) throw new RecordError('coordinate count mismatch');
    return out;
}

/* ------------------------------------------------------------------------- *
 * TRR
 *
 * Magic 1993, big-endian XDR. The version string is an XDR string (a length
 * followed by that many bytes, padded to four); thirteen size integers follow.
 * Each frame states its own byte sizes, so a frame can be skipped without
 * understanding the blocks it does not carry.
 * ------------------------------------------------------------------------- */

interface TrrFrameInfo {
    offset: number;
    natoms: number;
    step: number;
    time: number;
    box: number[];
    floatSize: number;
    xAt: number;
    xBytes: number;
    xSize: number;
}

function scanTrr(dv: DataView): { infos: TrrFrameInfo[]; truncated: boolean } {
    const infos: TrrFrameInfo[] = [];
    let offset = 0;
    let truncated = false;

    while (offset < dv.byteLength) {
        try {
            ensure(dv, offset, 12);
            const magic = dv.getInt32(offset, false);
            if (magic !== TRR_MAGIC) throw new RecordError('not a TRR frame');
            const stringLength = dv.getInt32(offset + 8, false);
            if (stringLength < 0 || stringLength > 4096) throw new RecordError('bad version string');
            const paddedString = (stringLength + 3) & ~3;
            let p = offset + 12 + paddedString;
            ensure(dv, p, 52);
            const boxSize = dv.getInt32(p + 8, false);
            const virSize = dv.getInt32(p + 12, false);
            const presSize = dv.getInt32(p + 16, false);
            const xSize = dv.getInt32(p + 28, false);
            const vSize = dv.getInt32(p + 32, false);
            const fSize = dv.getInt32(p + 36, false);
            const natoms = dv.getInt32(p + 40, false);
            const step = dv.getInt32(p + 44, false);
            const candidateSizes = [boxSize > 0 ? boxSize / 9 : 0, xSize > 0 ? xSize / (natoms * 3) : 0,
                vSize > 0 ? vSize / (natoms * 3) : 0, fSize > 0 ? fSize / (natoms * 3) : 0];
            const floatSize = candidateSizes.find(v => v === 4 || v === 8) || 0;
            if (!floatSize || natoms <= 0 || natoms > MAX_SINGLE_FRAME_ATOMS) {
                throw new RecordError('bad TRR frame header');
            }
            if (boxSize < 0 || virSize < 0 || presSize < 0 || xSize < 0 || vSize < 0 || fSize < 0) {
                throw new RecordError('negative TRR block size');
            }
            for (const size of [boxSize, virSize, presSize, xSize, vSize, fSize]) {
                if (size % floatSize !== 0) throw new RecordError('unaligned TRR block');
            }
            p += 52;
            ensure(dv, p, 2 * floatSize);
            const time = floatSize === 8 ? dv.getFloat64(p, false) : dv.getFloat32(p, false);
            p += 2 * floatSize;
            ensure(dv, p, boxSize);
            const box: number[] = [];
            if (boxSize >= 9 * floatSize) {
                for (let i = 0; i < 9; i++) {
                    box.push(floatSize === 8 ? dv.getFloat64(p + i * 8, false) : dv.getFloat32(p + i * 4, false));
                }
            }
            p += boxSize + virSize + presSize;
            ensure(dv, p, xSize + vSize + fSize);
            if (xSize > 0) {
                infos.push({ offset, natoms, step, time, box, floatSize, xAt: p, xBytes: xSize, xSize });
            }
            p += xSize + vSize + fSize;
            offset = p;
        } catch {
            truncated = offset < dv.byteLength;
            break;
        }
    }
    return { infos, truncated };
}

function decodeTrrFrame(dv: DataView, info: TrrFrameInfo): Float64Array {
    const count = info.natoms * 3;
    if (info.xBytes !== count * info.floatSize) throw new RecordError('short TRR coordinate block');
    const out = new Float64Array(count);
    for (let i = 0; i < count; i++) {
        out[i] = info.floatSize === 8
            ? dv.getFloat64(info.xAt + i * 8, false)
            : dv.getFloat32(info.xAt + i * 4, false);
    }
    return out;
}

/* ------------------------------------------------------------------------- *
 * DCD
 *
 * Fortran unformatted records with 32- or 64-bit length markers, in either
 * byte order. The header states the atom count, the fixed-atom count and the
 * CHARMM extensions; a CHARMM unit-cell block (six doubles per frame) and a
 * CHARMM 4th-dimension block are skipped or read according to those flags.
 * ------------------------------------------------------------------------- */

interface DcdReader {
    littleEndian: boolean;
    record64: boolean;
    markerSize: number;
}

interface DcdHeader {
    natoms: number;
    nsets: number;
    istart: number;
    nsavc: number;
    delta: number;
    nfixed: number;
    charmm: boolean;
    extraBlock: boolean;
    has4dims: boolean;
    freeIndexes: number[];
    firstOffset: number;
}

interface DcdFrameInfo {
    offset: number;
    step: number;
    time: number;
    cell?: number[][];
    xAt: number;
    yAt: number;
    zAt: number;
    /** Floats per axis in this frame (N, or N - fixed for later frames). */
    axisCount: number;
}

function asciiAt(bytes: Uint8Array, offset: number, text: string): boolean {
    if (offset + text.length > bytes.length) return false;
    for (let i = 0; i < text.length; i++) {
        if (bytes[offset + i] !== text.charCodeAt(i)) return false;
    }
    return true;
}

function detectDcd(bytes: Uint8Array, dv: DataView): DcdReader | null {
    if (bytes.length < 12) return null;
    if (asciiAt(bytes, 4, 'CORD')) {
        const le = dv.getInt32(0, true);
        const be = dv.getInt32(0, false);
        if (le === 84) return { littleEndian: true, record64: false, markerSize: 4 };
        if (be === 84) return { littleEndian: false, record64: false, markerSize: 4 };
        return null;
    }
    if (asciiAt(bytes, 8, 'CORD')) {
        let le = false;
        let be = false;
        try {
            le = Number(dv.getBigInt64(0, true)) === 84;
            be = Number(dv.getBigInt64(0, false)) === 84;
        } catch {
            return null;
        }
        if (le) return { littleEndian: true, record64: true, markerSize: 8 };
        if (be) return { littleEndian: false, record64: true, markerSize: 8 };
    }
    return null;
}

function dcdMarker(dv: DataView, reader: DcdReader, offset: number): number {
    ensure(dv, offset, reader.markerSize);
    const value = reader.record64
        ? Number(dv.getBigInt64(offset, reader.littleEndian))
        : dv.getInt32(offset, reader.littleEndian);
    if (!isFinite(value) || value < 0) throw new RecordError('bad DCD record marker');
    return value;
}

function readDcdHeader(bytes: Uint8Array, dv: DataView, reader: DcdReader): DcdHeader {
    const content = reader.markerSize;
    if (!asciiAt(bytes, content, 'CORD')) throw new RecordError('DCD magic not found');
    const header = content + 4;
    ensure(dv, header, 80);
    const int32 = (at: number) => dv.getInt32(header + at, reader.littleEndian);
    const version = int32(76);
    const charmm = version !== 0;
    const extraBlock = charmm && int32(40) !== 0;
    const has4dims = charmm && int32(44) === 1;
    const natoms = 0;
    const nsets = int32(0);
    const istart = int32(4);
    const nsavc = int32(8);
    const nfixed = int32(32);
    const delta = charmm
        ? dv.getFloat32(header + 36, reader.littleEndian)
        : dv.getFloat64(header + 36, reader.littleEndian);
    let offset = content + 84 + reader.markerSize;

    // Title record.
    const titleSize = dcdMarker(dv, reader, offset);
    offset += reader.markerSize;
    ensure(dv, offset, titleSize);
    offset += titleSize + reader.markerSize;

    // Atom-count record.
    const countSize = dcdMarker(dv, reader, offset);
    if (countSize !== 4) throw new RecordError('bad DCD atom-count record');
    offset += reader.markerSize;
    ensure(dv, offset, 4);
    const natomsValue = dv.getInt32(offset, reader.littleEndian);
    offset += 4 + reader.markerSize;

    // Fixed-atom free-index record.
    const freeIndexes: number[] = [];
    if (nfixed > 0 && nfixed < natomsValue) {
        const indexBytes = dcdMarker(dv, reader, offset);
        if (indexBytes !== (natomsValue - nfixed) * 4) throw new RecordError('bad DCD fixed-atom record');
        offset += reader.markerSize;
        ensure(dv, offset, indexBytes);
        for (let i = 0; i < natomsValue - nfixed; i++) {
            freeIndexes.push(dv.getInt32(offset + i * 4, reader.littleEndian));
        }
        offset += indexBytes + reader.markerSize;
    }
    return {
        natoms: natomsValue,
        nsets, istart, nsavc, delta,
        nfixed: nfixed > 0 && nfixed < natomsValue ? nfixed : 0,
        charmm, extraBlock, has4dims, freeIndexes,
        firstOffset: offset
    };
}

function dcdCellLattice(values: number[]): number[][] | undefined {
    const a = values[0], b = values[2], c = values[5];
    if (!finitePositive(a) || !finitePositive(b) || !finitePositive(c)) return undefined;
    let alpha = values[4];
    let beta = values[3];
    let gamma = values[1];
    if (alpha >= -1 && alpha <= 1 && beta >= -1 && beta <= 1 && gamma >= -1 && gamma <= 1) {
        // CHARMM and NAMD >= 2.5 write the cosines of the angles; the
        // asin formulation leaves orthogonal cells at exactly 90 degrees.
        alpha = 90 - Math.asin(alpha) * 90 / (Math.PI / 2);
        beta = 90 - Math.asin(beta) * 90 / (Math.PI / 2);
        gamma = 90 - Math.asin(gamma) * 90 / (Math.PI / 2);
    }
    if (![alpha, beta, gamma].every(isFinite)) return undefined;
    return usableCell(cellToLatticeVectors(a, b, c, alpha, beta, gamma));
}

function skipDcdRecord(dv: DataView, reader: DcdReader, offset: number): number {
    const size = dcdMarker(dv, reader, offset);
    ensure(dv, offset + reader.markerSize, size + reader.markerSize);
    return offset + reader.markerSize + size + reader.markerSize;
}

function scanDcd(dv: DataView, reader: DcdReader, header: DcdHeader): { infos: DcdFrameInfo[]; truncated: boolean } {
    const infos: DcdFrameInfo[] = [];
    let offset = header.firstOffset;
    let truncated = false;

    const limit = header.nsets > 0 ? header.nsets : Infinity;
    while (offset < dv.byteLength && infos.length < limit) {
        try {
            let cell: number[][] | undefined;
            if (header.extraBlock) {
                const size = dcdMarker(dv, reader, offset);
                offset += reader.markerSize;
                if (size === 48) {
                    ensure(dv, offset, 48);
                    const values: number[] = [];
                    for (let i = 0; i < 6; i++) values.push(dv.getFloat64(offset + i * 8, reader.littleEndian));
                    cell = dcdCellLattice(values);
                    offset += 48 + reader.markerSize;
                } else {
                    ensure(dv, offset, size + reader.markerSize);
                    offset += size + reader.markerSize;
                }
            }
            const axisCount = header.nfixed > 0 && infos.length > 0
                ? header.natoms - header.nfixed
                : header.natoms;
            // Each axis is its own Fortran record: leading marker, floats,
            // trailing marker. The next axis starts after both markers.
            const expected = axisCount * 4;
            const xMarker = offset;
            if (dcdMarker(dv, reader, xMarker) !== expected) throw new RecordError('bad DCD X record');
            const xAt = xMarker + reader.markerSize;
            ensure(dv, xAt, expected + reader.markerSize);
            const yMarker = xAt + expected + reader.markerSize;
            if (dcdMarker(dv, reader, yMarker) !== expected) throw new RecordError('bad DCD Y record');
            const yAt = yMarker + reader.markerSize;
            ensure(dv, yAt, expected + reader.markerSize);
            const zMarker = yAt + expected + reader.markerSize;
            if (dcdMarker(dv, reader, zMarker) !== expected) throw new RecordError('bad DCD Z record');
            const zAt = zMarker + reader.markerSize;
            ensure(dv, zAt, expected + reader.markerSize);
            offset = zAt + expected + reader.markerSize;
            if (header.has4dims) offset = skipDcdRecord(dv, reader, offset);

            // ISTART/NSAVC place the frame on the run's step axis; DELTA is the time
            // per step, so the frame's time is that step times DELTA. Using the
            // frame index alone put every DCD frame at its own index in time.
            const step = header.istart + infos.length * header.nsavc;
            const time = step * header.delta;
            infos.push({ offset: xMarker, step, time, cell, xAt, yAt, zAt, axisCount });
        } catch {
            truncated = offset < dv.byteLength;
            break;
        }
    }
    return { infos, truncated };
}

function readDcdAxis(dv: DataView, reader: DcdReader, at: number, count: number, out: Float64Array, axis: number): void {
    for (let i = 0; i < count; i++) {
        out[i * 3 + axis] = dv.getFloat32(at + i * 4, reader.littleEndian);
    }
}

function decodeDcdFrame(
    dv: DataView,
    reader: DcdReader,
    header: DcdHeader,
    info: DcdFrameInfo,
    firstFrame?: Float64Array
): Float64Array {
    const natoms = header.natoms;
    const out = new Float64Array(natoms * 3);
    if (header.nfixed === 0 || info.axisCount === natoms) {
        readDcdAxis(dv, reader, info.xAt, natoms, out, 0);
        readDcdAxis(dv, reader, info.yAt, natoms, out, 1);
        readDcdAxis(dv, reader, info.zAt, natoms, out, 2);
        return out;
    }
    if (!firstFrame) throw new RecordError('fixed-atom coordinates need the first frame');
    out.set(firstFrame.subarray(0, natoms * 3));
    const free = info.axisCount;
    const temp = new Float64Array(free * 3);
    readDcdAxis(dv, reader, info.xAt, free, temp, 0);
    readDcdAxis(dv, reader, info.yAt, free, temp, 1);
    readDcdAxis(dv, reader, info.zAt, free, temp, 2);
    for (let i = 0; i < free; i++) {
        const target = header.freeIndexes[i] - 1;
        if (target >= 0 && target < natoms) {
            out[target * 3] = temp[i * 3];
            out[target * 3 + 1] = temp[i * 3 + 1];
            out[target * 3 + 2] = temp[i * 3 + 2];
        }
    }
    return out;
}

/* ------------------------------------------------------------------------- *
 * Public entry points
 * ------------------------------------------------------------------------- */

/** The binary trajectory format the bytes are, or null for anything else. */
export function detectBinaryTrajectory(bytes: Uint8Array): TrajectoryFormat | null {
    if (bytes.length < 12) return null;
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const magic = view.getInt32(0, false);
    if (magic === XTC_MAGIC) return 'xtc';
    if (magic === TRR_MAGIC) return 'trr';
    if (detectDcd(bytes, view)) return 'dcd';
    return null;
}

/**
 * The frame indexes to keep when a trajectory is larger than the limits.
 *
 * Frames are taken at a uniform stride and the last frame is always kept, so
 * the sampled view still ends where the run ended. A text dump uses the same
 * policy as the binary readers, which is why it is exported.
 */
export function chooseTrajectoryFrameIndexes(
    count: number,
    natoms: number,
    maxFrames?: number,
    maxTotalAtoms?: number
): number[] {
    if (count <= 0) return [];
    const frameLimit = maxFrames !== undefined && maxFrames > 0 ? maxFrames : DEFAULT_MAX_FRAMES;
    const atomLimit = maxTotalAtoms !== undefined && maxTotalAtoms > 0 ? maxTotalAtoms : DEFAULT_MAX_TOTAL_ATOMS;
    let stride = 1;
    while (Math.ceil(count / stride) > frameLimit ||
        Math.ceil(count / stride) * natoms > atomLimit) {
        stride++;
        if (stride > count) return [0];
    }
    const chosen: number[] = [];
    for (let i = 0; i < count; i += stride) chosen.push(i);
    if (chosen[chosen.length - 1] !== count - 1) chosen[chosen.length - 1] = count - 1;
    return chosen;
}

function readXtc(bytes: Uint8Array, dv: DataView, scale: number, options: TrajectoryReadOptions, result: BinaryTrajectory): void {
    const scan = scanXtc(dv);
    result.truncated = scan.truncated;
    result.frameCount = scan.infos.length;
    if (scan.infos.length === 0) {
        result.warnings.push('the XTC file holds no complete frame');
        return;
    }
    result.natoms = scan.infos[0].natoms;
    const maxFrames = options.maxFrames || DEFAULT_MAX_FRAMES;
    const maxAtoms = options.maxTotalAtoms || DEFAULT_MAX_TOTAL_ATOMS;
    const chosen = chooseTrajectoryFrameIndexes(scan.infos.length, result.natoms, maxFrames, maxAtoms);
    result.sampled = chosen.length < scan.infos.length;

    for (const index of chosen) {
        const info = scan.infos[index];
        if (info.natoms !== result.natoms) {
            result.warnings.push('XTC frame ' + (index + 1) + ' has a different atom count: skipped');
            continue;
        }
        try {
            const raw = decodeXtcFrame(bytes, dv, info);
            const positions = new Float64Array(raw.length);
            for (let i = 0; i < raw.length; i++) positions[i] = raw[i] * scale;
            if (!allFinite(positions)) {
                result.warnings.push('XTC frame ' + (index + 1) + ' holds a non-finite coordinate: skipped');
                continue;
            }
            const frame: TrajectoryFrame = {
                positions, step: info.step,
                time: isFinite(info.time) ? info.time : undefined
            };
            const cell = usableCell(matrixFromNine(info.box));
            if (cell) frame.cell = scaleCell(cell, scale);
            result.frames.push(frame);
        } catch {
            result.warnings.push('XTC frame ' + (index + 1) + ' could not be decoded: skipped');
        }
    }
}

function readTrr(dv: DataView, scale: number, options: TrajectoryReadOptions, result: BinaryTrajectory): void {
    const scan = scanTrr(dv);
    result.truncated = scan.truncated;
    result.frameCount = scan.infos.length;
    if (scan.infos.length === 0) {
        result.warnings.push('the TRR file holds no frame with coordinates');
        return;
    }
    result.natoms = scan.infos[0].natoms;
    const maxFrames = options.maxFrames || DEFAULT_MAX_FRAMES;
    const maxAtoms = options.maxTotalAtoms || DEFAULT_MAX_TOTAL_ATOMS;
    const chosen = chooseTrajectoryFrameIndexes(scan.infos.length, result.natoms, maxFrames, maxAtoms);
    result.sampled = chosen.length < scan.infos.length;

    for (const index of chosen) {
        const info = scan.infos[index];
        if (info.natoms !== result.natoms) {
            result.warnings.push('TRR frame ' + (index + 1) + ' has a different atom count: skipped');
            continue;
        }
        try {
            const raw = decodeTrrFrame(dv, info);
            const positions = new Float64Array(raw.length);
            for (let i = 0; i < raw.length; i++) positions[i] = raw[i] * scale;
            if (!allFinite(positions)) {
                result.warnings.push('TRR frame ' + (index + 1) + ' holds a non-finite coordinate: skipped');
                continue;
            }
            const frame: TrajectoryFrame = {
                positions, step: info.step,
                time: isFinite(info.time) ? info.time : undefined
            };
            const cell = usableCell(matrixFromNine(info.box));
            if (cell) frame.cell = scaleCell(cell, scale);
            result.frames.push(frame);
        } catch {
            result.warnings.push('TRR frame ' + (index + 1) + ' could not be decoded: skipped');
        }
    }
}

function readDcd(bytes: Uint8Array, dv: DataView, reader: DcdReader, scale: number, options: TrajectoryReadOptions, result: BinaryTrajectory): void {
    let header: DcdHeader;
    try {
        header = readDcdHeader(bytes, dv, reader);
    } catch (err) {
        // A damaged header is a statement about this file, not a reason to
        // fail the document: the warning path still opens the panel.
        const msg = err instanceof Error ? err.message : String(err);
        result.warnings.push('the DCD header could not be read (' + msg + '): no trajectory frames');
        return;
    }
    if (header.natoms <= 0 || header.natoms > MAX_SINGLE_FRAME_ATOMS) {
        result.warnings.push('the DCD header states ' + header.natoms +
            ' atoms: the file is too large or corrupt to show');
        return;
    }
    result.natoms = header.natoms;
    const scan = scanDcd(dv, reader, header);
    result.truncated = scan.truncated;
    result.frameCount = scan.infos.length;
    if (scan.infos.length === 0) {
        result.warnings.push('the DCD file holds no complete frame');
        return;
    }
    const maxFrames = options.maxFrames || DEFAULT_MAX_FRAMES;
    const maxAtoms = options.maxTotalAtoms || DEFAULT_MAX_TOTAL_ATOMS;
    const chosen = chooseTrajectoryFrameIndexes(scan.infos.length, header.natoms, maxFrames, maxAtoms);
    result.sampled = chosen.length < scan.infos.length;

    // Fixed atoms are stored only in the first frame; every later frame must
    // borrow them, even when the first frame itself is not on display.
    let firstRaw: Float64Array | undefined;
    if (header.nfixed > 0) {
        try {
            firstRaw = decodeDcdFrame(dv, reader, header, scan.infos[0]);
        } catch {
            result.warnings.push('the first DCD frame could not be decoded: fixed atoms are missing');
        }
    }

    for (const index of chosen) {
        const info = scan.infos[index];
        try {
            const raw = decodeDcdFrame(dv, reader, header, info, firstRaw);
            const positions = new Float64Array(raw.length);
            for (let i = 0; i < raw.length; i++) positions[i] = raw[i] * scale;
            if (!allFinite(positions)) {
                result.warnings.push('DCD frame ' + (index + 1) + ' holds a non-finite coordinate: skipped');
                continue;
            }
            const frame: TrajectoryFrame = {
                positions, step: info.step,
                time: isFinite(info.time) ? info.time : undefined
            };
            if (info.cell) frame.cell = scaleCell(info.cell, scale);
            result.frames.push(frame);
        } catch {
            result.warnings.push('DCD frame ' + (index + 1) + ' could not be decoded: skipped');
        }
    }
}

/**
 * XTC/TRR store the frame's box as the nine numbers of the 3x3 matrix, rows
 * first (a, b, c), which is exactly the lattice the viewer draws. A `.gro`
 * file packs the same box differently (v1x v2y v3z v1y v1z v2x v2z v3x v3y),
 * so this value must not be read the way a `.gro` line is.
 */
function matrixFromNine(values: number[]): number[][] | undefined {
    if (values.length !== 9) return undefined;
    return [[values[0], values[1], values[2]],
            [values[3], values[4], values[5]],
            [values[6], values[7], values[8]]];
}

/**
 * Read a binary trajectory. Frames are sampled when the file would exceed the
 * limits; decoding failures cost those frames and leave a warning, so a file
 * with one bad record still opens on the rest.
 */
export function readBinaryTrajectory(bytes: Uint8Array, options: TrajectoryReadOptions = {}): BinaryTrajectory {
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const format = detectBinaryTrajectory(bytes);
    if (!format) {
        throw new Error('not an XTC, TRR or DCD trajectory');
    }
    const defaultScale = format === 'dcd' ? 1 : 10;
    const scale = options.lengthScale !== undefined && isFinite(options.lengthScale) && options.lengthScale > 0
        ? options.lengthScale
        : defaultScale;
    const result: BinaryTrajectory = {
        format,
        natoms: 0,
        frames: [],
        frameCount: 0,
        sampled: false,
        truncated: false,
        lengthScale: scale,
        warnings: []
    };
    if (format === 'xtc') {
        readXtc(bytes, view, scale, options, result);
    } else if (format === 'trr') {
        readTrr(view, scale, options, result);
    } else {
        const reader = detectDcd(bytes, view)!;
        readDcd(bytes, view, reader, scale, options, result);
    }
    if (result.truncated) {
        result.warnings.push('the trajectory ends in an incomplete record: showing the frames that could be read');
    }
    return result;
}

/** A file name whose extension names a binary trajectory. */
export function isTrajectoryFileName(fileName: string): boolean {
    return /\.(xtc|trr|dcd)$/i.test(fileName);
}
