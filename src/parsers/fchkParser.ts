import { Atom, CavitySphere, SolventCavity } from '../types';
import { elementFromZ } from '../tools/periodicTable';
import { LogChargeSet, LogEnergy, scrfModelOf } from './logParser';

/**
 * A Gaussian formatted checkpoint (`.fchk`) reader.
 *
 * A checkpoint carries what a log does not: the wavefunction, the basis set
 * each molecular orbital is expanded in, and the cavity the solvation model
 * built. It is also *large* — a 432-atom job's MO coefficient block alone is
 * 4680 x 4680 numbers and 333 MB of text — so this reader never turns the file
 * into a string, never splits it into lines and never holds the coefficients:
 * it streams the bytes once, keeps the small arrays it needs, and records the
 * byte offset of the big ones. One molecular orbital is then read back as a
 * single contiguous run of that file (the coefficients of one MO are adjacent,
 * because the array is written in Fortran order), which is what makes an
 * isosurface possible without a 350 MB heap.
 *
 * Two things are deliberately NOT derived here:
 *
 *  - **Occupations.** A formatted checkpoint states electron counts, not a
 *    per-orbital occupation vector. It is rebuilt from the counts and the
 *    job's spin kind (see `occupationsOf`), which is exact for a closed-shell
 *    or a UHF job and the standard ROHF pattern for an `IROHF` one.
 *  - **Anything the file does not print.** A field that is absent stays
 *    absent; nothing is defaulted, so a panel can never show a number the file
 *    did not state.
 */

/** Bohr to Angstrom, CODATA 2018 — the unit of every coordinate in an fchk. */
const BOHR = 0.529177210903;

/** Where one numeric array of the file lives, so a column can be read later. */
export interface FchkArrayLayout {
    /** Byte offset of the first value. */
    start: number;
    /** Values one line carries (5 for R arrays, 6 for I arrays). */
    perLine: number;
    /** Bytes one line occupies, line terminator included. */
    lineBytes: number;
    /** Bytes one value occupies. */
    fieldWidth: number;
    /** Bytes of line terminator (1 for LF, 2 for CRLF). */
    terminator: number;
    /** Values the array states. */
    count: number;
    /**
     * Whether every line this reader walked had exactly `lineBytes` and the
     * field width divides it. A file whose layout is not uniform cannot be
     * addressed by arithmetic, so the column reader refuses it instead of
     * returning numbers from the wrong place.
     */
    uniform: boolean;
}

/** One spin channel's orbitals: what the file states, plus where to read more. */
export interface FchkOrbitalSpin {
    /** `Alpha` or `Beta`, as the file names it. */
    label: string;
    /** Orbital energies in Hartree, in the file's own order. */
    energies: number[];
    /** Per-orbital occupation, rebuilt (see `occupationsOf`). */
    occupations: number[];
    /** How many of them the electron count says are occupied. */
    occupied: number;
    /** Where the coefficient columns live, when the file carries them. */
    layout?: FchkArrayLayout;
}

/** The contracted shells a molecular orbital is expanded in. */
export interface FchkBasis {
    /** Per shell: 0 s, 1 p, -1 sp, 2 cartesian d, -2 pure d, 3/10f, -3/7f. */
    shellType: number[];
    /** Per shell: primitives it is contracted from. */
    shellPrim: number[];
    /**
     * Per shell: the atom it is centred on, **0-based**.
     *
     * The file writes this one 1-based ("Shell to atom map" numbers the atoms
     * from one), and it is converted here so a reader can index the atom list
     * with it directly. Nothing else in the viewer reads the raw array, and the
     * orbital mesh needs it whenever a checkpoint states no shell coordinates
     * and the shells have to be placed from their atom.
     */
    shellAtom: number[];
    /** Per primitive, in shell order. */
    exponents: number[];
    /** Per primitive: the s part (or the whole contraction for s/p/d/f). */
    coefficients: number[];
    /** Per primitive: the p part of an `sp` shell; empty when the file has none. */
    spCoefficients: number[];
    /** Per shell, three coordinates in Bohr. */
    centers: number[];
    /** Largest angular momentum the file states. */
    maxL: number;
    /** Whether d shells are 6 Cartesian functions (true) or 5 pure ones. */
    cartesianD: boolean;
    /** Same for f shells. */
    cartesianF: boolean;
}

/** Everything needed to draw an orbital of this file. */
export interface FchkOrbitals {
    nBasis: number;
    spins: FchkOrbitalSpin[];
    basis: FchkBasis;
    /** File the coefficient columns are read from. */
    path: string;
    /** `Method/Basis` as the checkpoint's force line words it. */
    model: string;
}

/** One sphere of the solvation cavity, in Angstrom. */
export { CavitySphere, SolventCavity } from '../types';

export interface FchkResult {
    title: string;
    atoms: Atom[];
    charge?: number;
    multiplicity?: number;
    nAtoms: number;
    nBasis: number;
    nAlpha: number;
    nBeta: number;
    /** `Method/Basis`, from the checkpoint's force line. */
    model: string;
    /** Every energy the file states, under its own label. */
    energies: LogEnergy[];
    rmsGradient?: number;
    /** The population analyses the file carries, in the shape a log's are. */
    charges: LogChargeSet[];
    orbitals?: FchkOrbitals;
    cavity?: SolventCavity;
    /** What the reader had to work around, in the file's own terms. */
    warnings: string[];
}

export interface FchkReadOptions {
    /** Bytes the file holds, so progress can be a fraction. */
    size?: number;
    /** Called as the stream advances, at most a few dozen times. */
    onProgress?: (bytesRead: number) => void;
    /** Path recorded on `orbitals.path`, so a column can be read back. */
    path?: string;
}

/** A formatted checkpoint section header: `Name  <type>  N= <count>`. */
const ARRAY_HEADER_RE = /^(\S.*?[^ \t])\s+([IRC])\s+N=\s*(\d+)\s*$/;
/** A formatted checkpoint scalar: `Name  <type>  <value>`. */
const SCALAR_HEADER_RE = /^(\S.*?[^ \t])\s+([IR])\s+(\S+)\s*$/;

/** Values one line of an array carries, by the array's type letter. */
const PER_LINE: { [type: string]: number } = { R: 5, I: 6, C: 5 };

/**
 * Arrays whose name ends in ` Charges` are per-atom tables, but these are not
 * populations of the wavefunction and must not appear as a charge scheme.
 */
const CHARGE_BLOCKLIST = [
    'nuclear charges', 'mm charges', 'atom modified mm charges',
    'oniom charges', 'point charges', 'link charges'
];

/**
 * Parse a formatted checkpoint.
 *
 * The input is any async iterable of byte chunks (a `fs.createReadStream`), so
 * a 900 MB checkpoint costs a 64 kB read buffer plus whatever the small arrays
 * need: the coefficient block is skipped by counting its lines, and its offset
 * is recorded for `readFchkColumn`.
 */
export async function parseFchk(input: AsyncIterable<Uint8Array>, options: FchkReadOptions = {}): Promise<FchkResult> {
    const scan = new FchkScan(options);
    for await (const chunk of input) {
        scan.feed(chunk instanceof Buffer ? chunk : Buffer.from(chunk));
    }
    return scan.finish();
}

/**
 * Parse a checkpoint that is already in memory.
 *
 * The streaming reader is the one every caller should use — it is what keeps a
 * 900 MB checkpoint from becoming a string — but `parseFile` only ever receives
 * text, so a short checkpoint opened under a generic extension still gets the
 * same reader rather than the XYZ fallback. The bytes are taken as latin1,
 * which is what the line-oriented scan expects.
 */
export function parseFchkText(content: string): FchkResult {
    const scan = new FchkScan({});
    scan.feed(Buffer.from(content, 'latin1'));
    return scan.finish();
}

/** One wanted array being collected while its value lines stream past. */
interface Collector {
    name: string;
    type: string;
    /** Value lines the array occupies. */
    linesLeft: number;
    /** Values it states. */
    count: number;
    numbers: number[];
    text: string[];
}

/**
 * The byte-level state machine.
 *
 * Three modes, one per thing a line can be: a section header (`header`), a
 * value line of an array this reader keeps (`collect`), or a line of an array
 * it does not (`skip`). Skipping never builds a string and never parses a
 * number — it counts line terminators, which is what keeps a 333 MB
 * coefficient block free.
 */
class FchkScan {
    /** Bytes of the last chunk that did not end on a line boundary. */
    private pending: Buffer = Buffer.alloc(0);
    /** Absolute byte offset of `pending[0]`. */
    private base = 0;
    /** Bytes handed to `feed` so far, for progress. */
    private read = 0;
    /** When the last progress callback fired, so the stream is not throttled. */
    private lastReport = 0;

    private mode: 'header' | 'collect' | 'skip' = 'header';

    // --- what has been read so far ---
    private nAtoms = 0;
    private charge?: number;
    private multiplicity?: number;
    private nBasis = 0;
    private nAlpha = 0;
    private nBeta = 0;
    private irohf = 0;
    private atomicNumbers: number[] = [];
    private coordinates: number[] = [];
    private title = '';
    private model = '';
    private scalars: { name: string; value: number }[] = [];
    private shellType: number[] = [];
    private shellPrim: number[] = [];
    private shellAtom: number[] = [];
    private exponents: number[] = [];
    private coefficients: number[] = [];
    private spCoefficients: number[] = [];
    private shellCenters: number[] = [];
    private maxL = 0;
    private cartesianD = true;
    private cartesianF = false;
    private alphaEnergies: number[] = [];
    private betaEnergies: number[] = [];
    private alphaLayout?: FchkArrayLayout;
    private betaLayout?: FchkArrayLayout;
    private singleLayout?: FchkArrayLayout;
    private pcmRadii: number[] = [];
    private pcmCenters: number[] = [];
    private chargeTables: { label: string; values: number[] }[] = [];
    private warnings: string[] = [];

    // --- the collector in flight ---
    private collector: Collector | null = null;

    // --- the skip in flight ---
    private skipName = '';
    private skipLeft = 0;
    private skipCount = 0;
    private skipLineBytes = 0;
    private skipTerm = 1;
    private skipUniform = true;
    private skipStart = 0;
    /** Where the layout of the array being skipped goes when it ends. */
    private skipTarget: 'alpha' | 'beta' | 'single' | null = null;

    constructor(private readonly options: FchkReadOptions) {}

    feed(chunk: Buffer): void {
        this.read += chunk.length;
        const now = Date.now();
        if (this.options.onProgress && now - this.lastReport > 150) {
            this.lastReport = now;
            this.options.onProgress(this.read);
        }

        const buf = this.pending.length > 0 ? Buffer.concat([this.pending, chunk]) : chunk;
        const length = buf.length;
        let pos = 0;
        scan: for (;;) {
            if (this.mode === 'skip') {
                while (this.skipLeft > 0) {
                    const nl = buf.indexOf(10, pos);
                    // No complete line left: everything from `pos` is a partial
                    // line and stays in `pending` for the next chunk.
                    if (nl < 0) break scan;
                    const lineBytes = nl - pos + 1;
                    if (this.skipLineBytes === 0) {
                        this.skipLineBytes = lineBytes;
                        this.skipTerm = nl > pos && buf[nl - 1] === 13 ? 2 : 1;
                    } else if (lineBytes !== this.skipLineBytes && this.skipLeft > 1) {
                        // The last line of an array is short unless the count is
                        // a multiple of the stride, so it alone may differ.
                        this.skipUniform = false;
                    }
                    this.skipLeft--;
                    pos = nl + 1;
                }
                this.finishSkip();
                continue;
            }
            if (this.mode === 'collect') {
                const collector = this.collector;
                if (!collector) {
                    this.mode = 'header';
                    continue;
                }
                while (collector.linesLeft > 0) {
                    const nl = buf.indexOf(10, pos);
                    if (nl < 0) break scan;
                    let end = nl;
                    if (end > pos && buf[end - 1] === 13) end--;
                    this.takeValueLine(buf.toString('latin1', pos, end), collector);
                    collector.linesLeft--;
                    pos = nl + 1;
                }
                this.commit(collector);
                this.collector = null;
                this.mode = 'header';
                continue;
            }
            const nl = buf.indexOf(10, pos);
            if (nl < 0) break;
            let end = nl;
            if (end > pos && buf[end - 1] === 13) end--;
            this.onHeaderLine(buf.toString('latin1', pos, end), this.base + nl + 1);
            pos = nl + 1;
        }

        this.base += pos;
        this.pending = pos >= length ? Buffer.alloc(0) : Buffer.from(buf.subarray(pos));
    }

    /** Append one value line to the collector in flight. */
    private takeValueLine(line: string, collector: Collector): void {
        if (collector.type === 'C') {
            // Character arrays hold fixed 12-byte fields; joining the raw text
            // and trimming at the end is what the title and the route need.
            collector.text.push(line);
            return;
        }
        const parts = line.trim();
        if (parts === '') return;
        const fields = parts.split(/\s+/);
        for (let k = 0; k < fields.length; k++) {
            if (collector.numbers.length >= collector.count) break;
            const v = parseFloat(fields[k].replace(/[Dd]/, 'E'));
            if (Number.isFinite(v)) collector.numbers.push(v);
        }
    }

    /** Read one header line: a scalar, an array header, or nothing. */
    private onHeaderLine(text: string, dataStart: number): void {
        const line = text.replace(/\s+$/, '');
        if (line === '') return;

        const array = ARRAY_HEADER_RE.exec(line);
        if (array) {
            this.onArrayHeader(array[1].trim(), array[2], parseInt(array[3], 10), dataStart);
            return;
        }
        // A value line that the state machine should never see here; it must
        // not be mistaken for a section header.
        if (/^[ \t]/.test(text) || /^[-+0-9.]/.test(text)) return;
        // The checkpoint's second line is the force/method line: no type
        // column, two runs of text, and the only place a method is named.
        if (this.model === '' && /^Force\s/.test(line)) {
            this.model = line.slice(5).trim().replace(/\s{2,}/g, ' ');
            return;
        }

        const scalar = SCALAR_HEADER_RE.exec(line);
        if (!scalar) {
            if (this.title === '') this.title = line.trim();
            return;
        }
        const value = parseFloat(scalar[3].replace(/[Dd]/, 'E'));
        if (!Number.isFinite(value)) return;
        this.noteScalar(scalar[1].trim(), value);
    }

    private noteScalar(name: string, value: number): void {
        switch (name) {
            case 'Number of atoms': this.nAtoms = Math.round(value); return;
            case 'Charge': this.charge = Math.round(value); return;
            case 'Multiplicity': this.multiplicity = Math.round(value); return;
            case 'Number of basis functions': this.nBasis = Math.round(value); return;
            case 'Number of alpha electrons': this.nAlpha = Math.round(value); return;
            case 'Number of beta electrons': this.nBeta = Math.round(value); return;
            case 'IROHF': this.irohf = Math.round(value); return;
            case 'Highest angular momentum': this.maxL = Math.round(value); return;
            case 'Pure/Cartesian d shells': this.cartesianD = Math.round(value) === 1; return;
            case 'Pure/Cartesian f shells': this.cartesianF = Math.round(value) === 1; return;
            default: break;
        }
        // Every other scalar is kept when it reads like a result: the summary
        // shows energies, so they are collected under their own labels and
        // nothing is invented for a name this reader does not know.
        if (/Energy$/.test(name) || name === 'RMS Force' || name === 'Virial Ratio') {
            this.scalars.push({ name, value });
        }
    }

    private onArrayHeader(name: string, type: string, count: number, dataStart: number): void {
        const perLine = PER_LINE[type] || 5;
        if (type === 'C') {
            if (name === 'Full Title' || name === 'Molecule Name' || name === 'Route') {
                this.collector = { name, type, linesLeft: Math.ceil(count / perLine), count, numbers: [], text: [] };
                this.mode = 'collect';
                return;
            }
            this.beginSkip(name, count, perLine, dataStart, null);
            return;
        }
        if (!this.wantsArray(name, count)) {
            this.beginSkip(name, count, perLine, dataStart, this.layoutTargetFor(name));
            return;
        }
        this.collector = { name, type, linesLeft: Math.ceil(count / perLine), count, numbers: [], text: [] };
        this.mode = 'collect';
    }

    /** Whether a numeric array is one this reader keeps. */
    private wantsArray(name: string, count: number): boolean {
        switch (name) {
            case 'Atomic numbers':
            case 'Current cartesian coordinates':
            case 'Shell types':
            case 'Number of primitives per shell':
            case 'Shell to atom map':
            case 'Primitive exponents':
            case 'Contraction coefficients':
            case 'P(S=P) Contraction coefficients':
            case 'Coordinates of each shell':
            case 'Alpha Orbital Energies':
            case 'Beta Orbital Energies':
            case 'Orbital Energies':
            case 'PCM-SphereRadii':
            case 'PCM-SphereLoc':
                return true;
            default: break;
        }
        if (/Charges$/.test(name)) {
            if (CHARGE_BLOCKLIST.indexOf(name.toLowerCase()) >= 0) return false;
            // A per-atom population is normally exactly as long as the
            // structure. A different length is still collected rather than
            // skipped: `chargesOf` fills the missing slots with null and says
            // how many the file stated, so a truncated table (or a QM-region
            // population of an ONIOM job) stays visible as partial data
            // instead of disappearing from the panel without a word.
            return this.nAtoms > 0 && count > 0;
        }
        return false;
    }

    /** The coefficient layout a skipped array describes, if it is one. */
    private layoutTargetFor(name: string): 'alpha' | 'beta' | 'single' | null {
        if (name === 'Alpha MO coefficients') return 'alpha';
        if (name === 'Beta MO coefficients') return 'beta';
        if (name === 'MO coefficients') return 'single';
        return null;
    }

    private beginSkip(name: string, count: number, perLine: number, dataStart: number, target: 'alpha' | 'beta' | 'single' | null): void {
        this.skipName = name;
        this.skipCount = count;
        this.skipLeft = Math.ceil(count / perLine);
        this.skipLineBytes = 0;
        this.skipTerm = 1;
        this.skipUniform = true;
        this.skipStart = dataStart;
        this.skipTarget = target;
        this.mode = 'skip';
    }

    /** Turn a finished skip into a layout, when the array was a coefficient block. */
    private finishSkip(): void {
        if (this.skipTarget) {
            const perLine = 5;
            const lineBytes = this.skipLineBytes;
            const fieldWidth = lineBytes > 0 ? Math.round((lineBytes - this.skipTerm) / perLine) : 0;
            const layout: FchkArrayLayout = {
                start: this.skipStart,
                perLine,
                lineBytes,
                fieldWidth,
                terminator: this.skipTerm,
                count: this.skipCount,
                uniform: this.skipUniform && fieldWidth > 0 && fieldWidth * perLine + this.skipTerm === lineBytes
            };
            if (this.skipTarget === 'alpha') this.alphaLayout = layout;
            else if (this.skipTarget === 'beta') this.betaLayout = layout;
            else this.singleLayout = layout;
        }
        this.skipTarget = null;
        this.skipName = '';
        this.mode = 'header';
    }

    private commit(collector: Collector): void {
        const name = collector.name;
        if (collector.type === 'C') {
            const text = collector.text.join('').replace(/\0/g, '').replace(/\s+$/, '');
            if (name === 'Route') {
                if (text.trim() !== '') this.model = text.trim();
            } else if (text.trim() !== '') {
                this.title = text.trim();
            }
            return;
        }
        const v = collector.numbers;
        switch (name) {
            case 'Atomic numbers': this.atomicNumbers = v; return;
            case 'Current cartesian coordinates': this.coordinates = v; return;
            case 'Shell types': this.shellType = v.map(x => Math.round(x)); return;
            case 'Number of primitives per shell': this.shellPrim = v.map(x => Math.round(x)); return;
            case 'Shell to atom map': this.shellAtom = v.map(x => Math.max(0, Math.round(x) - 1)); return;
            case 'Primitive exponents': this.exponents = v; return;
            case 'Contraction coefficients': this.coefficients = v; return;
            case 'P(S=P) Contraction coefficients': this.spCoefficients = v; return;
            case 'Coordinates of each shell': this.shellCenters = v; return;
            case 'Alpha Orbital Energies':
            case 'Orbital Energies': this.alphaEnergies = v; return;
            case 'Beta Orbital Energies': this.betaEnergies = v; return;
            case 'PCM-SphereRadii': this.pcmRadii = v; return;
            case 'PCM-SphereLoc': this.pcmCenters = v; return;
            default:
                if (/Charges$/.test(name)) {
                    // The file's own label, without the plural the panel's rows
                    // do not need: `Mulliken Charges` is the scheme's name.
                    this.chargeTables.push({ label: name.replace(/\s+Charges$/, ''), values: v });
                }
                return;
        }
    }

    finish(): FchkResult {
        // A file need not end with a newline. The last line has then never been
        // handed to the state machine, and dropping it would silently truncate
        // the last array (a charge table one value short is refused whole).
        if (this.pending.length > 0) {
            const raw = this.pending.toString('latin1');
            const line = raw.endsWith('\r') ? raw.slice(0, -1) : raw;
            this.pending = Buffer.alloc(0);
            if (this.mode === 'skip') {
                if (this.skipLeft <= 1) {
                    this.skipLeft = 0;
                    this.finishSkip();
                } else {
                    this.skipLeft--;
                }
            } else if (this.mode === 'collect' && this.collector) {
                this.takeValueLine(line, this.collector);
                this.collector.linesLeft--;
                if (this.collector.linesLeft <= 0) {
                    this.commit(this.collector);
                    this.collector = null;
                    this.mode = 'header';
                }
            } else {
                this.onHeaderLine(line, this.base);
            }
        }
        if (this.collector) {
            // A file that ends inside an array still has whatever it printed; a
            // truncated job is a real case and must open rather than throw.
            this.commit(this.collector);
            this.warnings.push('the file ended inside the "' + this.collector.name + '" section: what it stated was read');
            this.collector = null;
        }
        if (this.nAtoms <= 0 || this.atomicNumbers.length === 0 || this.coordinates.length < this.nAtoms * 3) {
            throw new Error('no structure in this formatted checkpoint: it states ' +
                this.nAtoms + ' atoms, ' + this.atomicNumbers.length + ' atomic numbers and ' +
                Math.floor(this.coordinates.length / 3) + ' positions');
        }

        const atoms: Atom[] = [];
        for (let i = 0; i < this.nAtoms; i++) {
            atoms.push({
                element: elementFromZ(this.atomicNumbers[i]) || 'X',
                x: this.coordinates[i * 3] * BOHR,
                y: this.coordinates[i * 3 + 1] * BOHR,
                z: this.coordinates[i * 3 + 2] * BOHR,
                index: i
            });
        }

        const result: FchkResult = {
            title: this.title || 'Formatted checkpoint',
            atoms,
            charge: this.charge,
            multiplicity: this.multiplicity,
            nAtoms: this.nAtoms,
            nBasis: this.nBasis || this.alphaEnergies.length,
            nAlpha: this.nAlpha,
            nBeta: this.nBeta,
            model: this.model,
            energies: this.energiesOf(),
            charges: this.chargesOf(atoms),
            warnings: this.warnings
        };
        for (const s of this.scalars) {
            if (s.name === 'RMS Force') result.rmsGradient = s.value;
        }

        const orbitals = this.orbitalsOf();
        if (orbitals) result.orbitals = orbitals;
        const cavity = this.cavityOf(atoms);
        if (cavity) result.cavity = cavity;
        return result;
    }

    /** The energies the file states, de-duplicated by value. */
    private energiesOf(): LogEnergy[] {
        const out: LogEnergy[] = [];
        for (const s of this.scalars) {
            if (!/Energy$/.test(s.name)) continue;
            // `Total Energy` repeats the method's own line in a job that also
            // printed `SCF Energy`; one number under two names is noise.
            if (out.some(e => e.value === s.value)) continue;
            out.push({ label: s.name, value: s.value });
        }
        return out;
    }

    /** The population analyses, in the shape the charge panel already reads. */
    private chargesOf(atoms: Atom[]): LogChargeSet[] {
        const out: LogChargeSet[] = [];
        for (const table of this.chargeTables) {
            const values: (number | null)[] = [];
            for (let i = 0; i < atoms.length; i++) {
                const q = i < table.values.length ? table.values[i] : NaN;
                values.push(Number.isFinite(q) ? q : null);
            }
            // A table whose value line was truncated or held a field that is
            // not a number is still a table the file printed; its unreadable
            // slots stay null so the panel can show "(n/total)" and let the
            // user see which atoms have no value. Skipping the table silently
            // made the analysis disappear instead.
            if (table.values.length !== atoms.length) {
                this.warnings.push('the "' + table.label + '" charge table states ' + table.values.length +
                    ' values for ' + atoms.length + ' atoms: ' +
                    (table.values.length < atoms.length
                        ? 'the missing values are left blank'
                        : 'the extra values are ignored'));
            }
            const covered = values.filter(q => q !== null).length;
            out.push({
                label: table.label,
                values,
                elements: atoms.map(a => a.element),
                covered,
                // The file can state a table that is shorter than the atom list
                // or whose unparsable fields were left as null. That is a partial
                // table, and the panel's `(n/total)` suffix is what tells the
                // user so; claiming completeness would hide the missing values.
                complete: covered === atoms.length,
                frame: 0
            });
        }
        return out;
    }

    private orbitalsOf(): FchkOrbitals | undefined {
        const alphaLayout = this.alphaLayout || this.singleLayout;
        if (this.alphaEnergies.length === 0 && !alphaLayout) return undefined;
        const nBasis = this.nBasis || this.alphaEnergies.length;
        const basis: FchkBasis = {
            shellType: this.shellType,
            shellPrim: this.shellPrim,
            shellAtom: this.shellAtom,
            exponents: this.exponents,
            coefficients: this.coefficients,
            spCoefficients: this.spCoefficients,
            centers: this.shellCenters,
            maxL: Math.max(this.maxL, highestAngularMomentum(this.shellType)),
            cartesianD: this.cartesianD,
            cartesianF: this.cartesianF
        };
        // A closed shell is a file with ONE set of orbitals. Equal alpha and
        // beta electron counts do not prove it: a broken-symmetry UHF/UKS
        // singlet also has nAlpha === nBeta and states multiplicity 1, and
        // reading it as closed shell threw its beta orbitals away and forced
        // two electrons into every alpha level. The checkpoint states the spin
        // kind by storing a second set of MO coefficients (and its own orbital
        // energies), so the presence of that channel is what decides.
        const hasBetaChannel = this.betaLayout !== undefined || this.betaEnergies.length > 0;
        const closed = !hasBetaChannel && this.nBeta === this.nAlpha &&
            (this.multiplicity === undefined || this.multiplicity === 1) && this.irohf === 0;
        const spins: FchkOrbitalSpin[] = [];
        if (this.alphaEnergies.length > 0) {
            spins.push({
                label: 'Alpha',
                energies: this.alphaEnergies,
                occupations: occupationsOf(this.alphaEnergies.length, this.nAlpha, this.nBeta, closed, false, this.irohf !== 0),
                occupied: this.nAlpha,
                layout: alphaLayout
            });
        }
        if (this.betaLayout && !closed) {
            // The beta coefficient block can exist while the file states no
            // beta orbital energies. Reusing the alpha list there put alpha
            // numbers under the Beta heading, and an isosurface chosen by that
            // list could be a different orbital than the one the user read.
            // A channel whose levels the file does not state is not shown;
            // the missing block is named instead of invented.
            if (this.betaEnergies.length === 0) {
                this.warnings.push('the checkpoint stores beta MO coefficients but no beta orbital energies: ' +
                    'the beta channel is not shown, because its levels would be guessed');
            } else {
                spins.push({
                    label: 'Beta',
                    energies: this.betaEnergies,
                    occupations: occupationsOf(this.betaEnergies.length, this.nBeta, this.nBeta, false, true, this.irohf !== 0),
                    occupied: this.nBeta,
                    layout: this.betaLayout
                });
            }
        }
        if (spins.length === 0) return undefined;
        return { nBasis, spins, basis, path: this.options.path || '', model: this.model };
    }

    private cavityOf(atoms: Atom[]): SolventCavity | undefined {
        if (this.pcmRadii.length === 0 || this.pcmCenters.length < this.pcmRadii.length * 3) return undefined;
        const spheres: CavitySphere[] = [];
        for (let i = 0; i < this.pcmRadii.length; i++) {
            spheres.push({
                x: this.pcmCenters[i * 3] * BOHR,
                y: this.pcmCenters[i * 3 + 1] * BOHR,
                z: this.pcmCenters[i * 3 + 2] * BOHR,
                r: this.pcmRadii[i] * BOHR
            });
        }
        if (spheres.length !== atoms.length) {
            this.warnings.push('the cavity lists ' + spheres.length + ' spheres for ' + atoms.length +
                ' atoms: the spheres are drawn where the file put them');
        }
        // The route names the solver (`SCRF=(SMD, …)`); when it does not, the
        // file still stores its spheres in the `PCM-*` arrays, which is the
        // model family those arrays belong to.
        const model = scrfModelOf(this.model) || 'PCM';
        return { spheres, sphereCount: spheres.length, model, source: 'fchk' };
    }
}

/** The largest angular momentum a shell-type list implies. */
function highestAngularMomentum(types: number[]): number {
    let max = 0;
    for (const t of types) {
        const l = t === 0 ? 0 : Math.abs(t);
        if (l > max) max = l;
    }
    return max;
}

/**
 * Per-orbital occupations this reader can defend.
 *
 * A formatted checkpoint states `Number of alpha electrons` and `Number of
 * beta electrons`, never an occupation vector, so the vector is rebuilt from
 * the counts:
 *
 *  - closed shell (`nAlpha === nBeta`, singlet, not ROHF, no beta channel):
 *    the first `nAlpha` orbitals hold two electrons, the rest none;
 *  - unrestricted: each spin's first `nAlpha`/`nBeta` orbitals hold one;
 *  - ROHF: the first `nBeta` alpha orbitals are the doubly occupied pair, the
 *    remaining `nAlpha - nBeta` are singly occupied, and the beta channel lists
 *    the same pairs holding ONE beta electron each — writing 2 there counted
 *    every paired beta electron twice, since the alpha row already states the
 *    pair.
 */
function occupationsOf(count: number, nAlpha: number, nBeta: number, closed: boolean, beta: boolean, rohf: boolean): number[] {
    const out = new Array<number>(count).fill(0);
    if (closed) {
        for (let i = 0; i < Math.min(count, nAlpha); i++) out[i] = 2;
        return out;
    }
    if (rohf) {
        if (beta) {
            // One electron per occupied beta spin orbital.
            for (let i = 0; i < Math.min(count, nBeta); i++) out[i] = 1;
            return out;
        }
        for (let i = 0; i < Math.min(count, nBeta); i++) out[i] = 2;
        for (let i = nBeta; i < Math.min(count, nAlpha); i++) out[i] = 1;
        return out;
    }
    const occupied = beta ? nBeta : nAlpha;
    for (let i = 0; i < Math.min(count, occupied); i++) out[i] = 1;
    return out;
}

/**
 * Read one molecular orbital's coefficient column, in the file's AO order.
 *
 * The coefficient array is written in Fortran order (all AOs of MO 1, then all
 * of MO 2, …), so one orbital is one contiguous run of `perLine`-value lines.
 * The run is located by arithmetic on the stride the checkpoint's own lines
 * state, which is why the layout is verified while the section streams past: a
 * file whose lines are not uniform is refused rather than read at a guessed
 * offset.
 *
 * `read` is the caller's ranged reader (a buffered `fs.read`), kept out of this
 * module so the parser owns no file handle.
 */
export async function readFchkColumn(
    layout: FchkArrayLayout,
    nBasis: number,
    orbital: number,
    read: (start: number, length: number) => Promise<Buffer>
): Promise<Float64Array> {
    if (!layout.uniform) {
        throw new Error('the coefficient block of this checkpoint is not written in a uniform line layout, ' +
            'so one orbital cannot be addressed by offset');
    }
    if (nBasis <= 0 || orbital < 0 || (orbital + 1) * nBasis > layout.count) {
        throw new Error('orbital ' + (orbital + 1) + ' is outside the ' +
            (nBasis > 0 ? Math.floor(layout.count / nBasis) : 0) + ' orbitals this checkpoint stores');
    }
    const first = orbital * nBasis;
    const last = first + nBasis - 1;
    const firstLine = Math.floor(first / layout.perLine);
    const lastLine = Math.floor(last / layout.perLine);
    const lineCount = Math.ceil(layout.count / layout.perLine);
    const lastLineValues = layout.count - (lineCount - 1) * layout.perLine;
    let length = (lastLine - firstLine + 1) * layout.lineBytes;
    if (lastLine === lineCount - 1) {
        // The array's last line holds only the values its count leaves over, so
        // its bytes stop before the next section begins. Reading a full line
        // there would parse the *next* section's characters as coefficients.
        length = (lastLine - firstLine) * layout.lineBytes + lastLineValues * layout.fieldWidth + layout.terminator;
    }
    const bytes = await read(layout.start + firstLine * layout.lineBytes, length);
    // The final line of the array may be the last line of the file, where a
    // writer can omit the line terminator; the values themselves must still be
    // present. Anything shorter is a short read (or a truncated file), and it
    // used to be filled with zeros here — the orbital was then drawn from a
    // column that was mostly invented. Refusing it is the only honest answer.
    const needed = length - (lastLine === lineCount - 1 ? layout.terminator : 0);
    if (bytes.length < needed) {
        throw new Error('the coefficient block of this checkpoint ended after ' + bytes.length +
            ' of the ' + needed + ' bytes orbital ' + (orbital + 1) + ' needs; the file is truncated');
    }
    const out = new Float64Array(nBasis);
    for (let k = 0; k < nBasis; k++) {
        const global = first + k;
        const line = Math.floor(global / layout.perLine) - firstLine;
        const field = global % layout.perLine;
        const at = line * layout.lineBytes + field * layout.fieldWidth;
        const text = bytes.toString('latin1', at, at + layout.fieldWidth).replace(/[Dd]/, 'E');
        const v = parseFloat(text);
        if (!Number.isFinite(v)) {
            // A field that is not a number (damage, a shifted layout) is not a
            // zero coefficient: refusing it keeps a wrong orbital off the screen.
            throw new Error('coefficient ' + (k + 1) + ' of orbital ' + (orbital + 1) +
                ' is not a number ("' + text.trim() + '"); the checkpoint may be damaged');
        }
        out[k] = v;
    }
    return out;
}
