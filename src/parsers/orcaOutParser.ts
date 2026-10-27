import { Atom, Bond, NormalMode, OptStep, OptThresholds } from '../types';
import { RouteKeyword, RouteSection } from './logParser';

/**
 * ORCA output (`.out`) parser.
 *
 * ORCA writes the same information a Gaussian log does, in a different shape,
 * and the viewer's optimization/frequency features were written against the
 * Gaussian shape. This module produces that shape from an ORCA file, so the
 * convergence panel, the normal-mode playback and the route-keyword panel work
 * on ORCA runs without the webview knowing which program wrote the file:
 *
 *   Gaussian                          ORCA
 *   ------------------------------    ------------------------------------------
 *   `Standard orientation:`           `CARTESIAN COORDINATES (ANGSTROEM)` + `(A.U.)`
 *   `SCF Done: E(...) = ...`          `FINAL SINGLE POINT ENERGY      ...`
 *   Item/Value/Converged? table       `|Geometry convergence|` table
 *   route card `# ...`                input echo `|  1> ! ...`
 *   `Harmonic frequencies ... cm**-1` `VIBRATIONAL FREQUENCIES` + `NORMAL MODES`
 *
 * Everything is found in ONE pass over the file. Each block is recognized by
 * its own header line, and every block parser consumes exactly its own lines
 * and hands back the index of the first line it did not use, so no block is
 * scanned twice and no header is seen from the middle of another block.
 *
 * The parser is deliberately tolerant about what it does not recognize: an
 * unknown line is skipped, a block whose rows do not parse is dropped rather
 * than half-read, and a file with no recognizable block yields zero frames
 * (which the caller reports as "no structure found") instead of throwing.
 */

/** ORCA writes A.U. coordinates; Bohr → Å uses the same constant the viewer has
 *  always used for the A.U. blocks, so switching to that block as the fallback
 *  cannot shift an existing structure. */
const BOHR_TO_ANGSTROM = 0.529177249;

/** Every element ORCA can print in a coordinate block. Used to tell a real atom
 *  row from a column-header row in the A.U. block, whose header starts with
 *  `NO` / `LB` / `FRAG` — tokens that would otherwise parse as symbols (`NO` is
 *  nobelium, `Y` is yttrium). */
const ELEMENT_SYMBOLS = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ' +
    'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce ' +
    'Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ' +
    'Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl ' +
    'Mc Lv Ts Og'
).split(' ');
const ELEMENT_SET: { [sym: string]: true } = {};
for (const s of ELEMENT_SYMBOLS) ELEMENT_SET[s] = true;

/** `O`, `cl`, `CL`, `H` → `O`, `Cl`, `Cl`, `H`; anything that is not an element
 *  symbol → null. Case is folded because ORCA's own output is capitalized but
 *  hand-written and third-party writers are not always. */
function normalizeElement(raw: string): string | null {
    if (!/^[A-Za-z]{1,3}$/.test(raw)) return null;
    const sym = raw.charAt(0).toUpperCase() + raw.slice(1).toLowerCase();
    return ELEMENT_SET[sym] ? sym : null;
}

/** A pure dashed separator (`-----`, ORCA's 26-to-80 character underlines). */
function isDashLine(t: string): boolean {
    return t.length >= 3 && /^-+$/.test(t);
}

export interface OrcaFrame {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    stepLabel: string;
    charge?: number;
    multiplicity?: number;
    /** Single-point energy (Hartree) of the geometry in this frame.
     *
     *  ORCA prints the geometry of a cycle BEFORE the energy computed at it
     *  (`CARTESIAN COORDINATES` … then `FINAL SINGLE POINT ENERGY`), so the
     *  energy is bound backwards to the frames still waiting for one — see
     *  `bindEnergy`. The convergence panel uses this to put its highlight on the
     *  step whose structure is on screen. */
    energy?: number;
}

export interface OrcaLogResult {
    frames: OrcaFrame[];
    /** The job's simple-input line (`! B3LYP def2-SVP Opt freq`) — ORCA's
     *  counterpart of a Gaussian job title. */
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    routes?: RouteSection[];
    optThresholds?: OptThresholds;
    charge?: number;
    multiplicity?: number;
}

/**
 * Whether a file is ORCA output.
 *
 * Each test below is a string only ORCA writes, so a Gaussian log can never be
 * taken for an ORCA one (and vice versa) — which matters because both programs
 * are routinely told to write `.out`. The tests run in the order the strings
 * appear in a run, so the common case decides on the first one or two.
 */
export function isOrcaOutput(content: string): boolean {
    if (content.indexOf('FINAL SINGLE POINT ENERGY') >= 0) return true;
    if (/^\s*\*\s*O\s+R\s+C\s+A\s*\*\s*$/m.test(content)) return true;   // the banner
    if (/CARTESIAN COORDINATES \((?:ANGSTROEM|ANGSTROM|A\.U\.)\)/.test(content)) return true;
    if (content.indexOf('ORCA TERMINATED NORMALLY') >= 0) return true;
    return false;
}

/**
 * Parse one `CARTESIAN COORDINATES (ANGSTROEM)` block.
 *
 * Rows are `ELEM x y z` (with an optional leading index in some hand-written
 * and third-party variants). The block ends at the first line that is not such
 * a row — a blank line, the next underline, or the next block's header — and the
 * index of that line is returned so the caller resumes exactly there.
 *
 * `complete` reports whether the block stopped on such a terminator. It is false
 * only when the file itself ended inside the block, which is what a killed job
 * leaves behind — the caller uses it to drop that one truncated geometry
 * instead of drawing a molecule with its last atoms missing.
 */
function parseAngstromBlock(lines: string[], start: number): { atoms: Atom[]; end: number } {
    const atoms: Atom[] = [];
    let i = start + 1;
    let skipped = 0;

    while (i < lines.length) {
        const t = lines[i].trim();
        if (t === '' || isDashLine(t)) {
            // Before the first atom these are the header's own underline and
            // the blank line around it; after it they end the block.
            if (atoms.length > 0) break;
            if (++skipped > 3) break;
            i++;
            continue;
        }
        const parts = t.split(/\s+/);
        // Optional leading integer index: `1 O 0.000000 0.000000 0.000000`.
        const off = parts.length >= 5 && /^\d+$/.test(parts[0]) ? 1 : 0;
        const element = normalizeElement(parts[off] || '');
        if (!element) {
            if (atoms.length === 0 && ++skipped <= 3) { i++; continue; }
            break;
        }
        const x = parseFloat(parts[off + 1]);
        const y = parseFloat(parts[off + 2]);
        const z = parseFloat(parts[off + 3]);
        if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
        atoms.push({ element, x, y, z, index: atoms.length });
        i++;
    }
    return { atoms, end: i };
}

/**
 * Parse one `CARTESIAN COORDINATES (A.U.)` block.
 *
 * Rows are `index ELEM ZA FRAG MASS x y z` in ORCA 4/5/6 (older versions omit
 * `FRAG`), so the element is the first token that is a real symbol and the
 * coordinates are the last three numbers on the line — a shape that survives
 * every column variant ORCA has used. The column header (`NO LB ZA FRAG MASS
 * X Y Z`) has no numbers after its `NO`, so it falls through the same path and
 * is skipped as a header rather than read as a one-atom row.
 */
function parseAuBlock(lines: string[], start: number): { atoms: Atom[]; end: number } {
    const atoms: Atom[] = [];
    let i = start + 1;
    let skipped = 0;

    while (i < lines.length) {
        const t = lines[i].trim();
        if (t === '' || isDashLine(t)) {
            if (atoms.length > 0) break;
            if (++skipped > 3) break;
            i++;
            continue;
        }
        const parts = t.split(/\s+/);
        let element: string | null = null;
        let at = -1;
        for (let k = 0; k < Math.min(parts.length, 4); k++) {
            element = normalizeElement(parts[k]);
            if (element) { at = k; break; }
        }
        if (!element || at < 0) {
            if (atoms.length === 0 && ++skipped <= 3) { i++; continue; }
            break;
        }
        const nums: number[] = [];
        for (let k = at + 1; k < parts.length; k++) {
            const v = parseFloat(parts[k]);
            if (isFinite(v)) nums.push(v);
        }
        if (nums.length < 3) {
            if (atoms.length === 0 && ++skipped <= 3) { i++; continue; }
            break;
        }
        atoms.push({
            element,
            x: nums[nums.length - 3] * BOHR_TO_ANGSTROM,
            y: nums[nums.length - 2] * BOHR_TO_ANGSTROM,
            z: nums[nums.length - 1] * BOHR_TO_ANGSTROM,
            index: atoms.length
        });
        i++;
    }
    return { atoms, end: i };
}

/**
 * Skip a coordinate block without parsing it, returning the index of the first
 * line the caller should look at again. Used for the A.U. half of a cycle whose
 * Ångström half was already read: it describes the same geometry in another
 * unit, and parsing it would double the cost of every optimization cycle for a
 * result that is always discarded.
 */
function skipCoordinateBlock(lines: string[], start: number): number {
    let i = start + 1;
    while (i < lines.length && (lines[i].trim() === '' || isDashLine(lines[i].trim()))) i++;
    while (i < lines.length) {
        const t = lines[i].trim();
        if (t === '' || isDashLine(t)) break;
        i++;
    }
    return i;
}

/** How far back one energy may travel. ORCA prints one geometry per energy, but
 *  a file that prints a second geometry block per cycle must not have its older
 *  frames re-labelled by a later energy. */
const MAX_FRAMES_PER_ENERGY = 8;

/**
 * Keep only the frames that describe the molecule the run is about.
 *
 * Every geometry of one ORCA output has the same atom count, so a frame with a
 * different one is not a geometry of this run: it is a block a killed job left
 * half written (the file simply ends in the middle of it, and the rows that did
 * get written parse into a shorter "molecule"), or a stray block from a source
 * concatenated onto the file. Both are worse than showing nothing — a truncated
 * structure silently loses bonds and atoms.
 *
 * The count that wins is the most frequent one, so a run whose FIRST block was
 * the truncated one is treated exactly like a run whose last block was, and the
 * tie-break takes the larger count (the complete molecule) rather than the
 * smaller. A file where every block has the same count — the ordinary case —
 * passes through untouched.
 */
function dropOddSizedFrames<T extends { atoms: Atom[] }>(frames: T[]): T[] {
    if (frames.length <= 1) return frames;
    const seen = new Map<number, number>();
    for (const f of frames) seen.set(f.atoms.length, (seen.get(f.atoms.length) || 0) + 1);
    let best = -1;
    let bestCount = 0;
    seen.forEach((count, size) => {
        if (count > bestCount || (count === bestCount && size > best)) {
            best = size;
            bestCount = count;
        }
    });
    if (best < 0) return frames;
    return frames.filter(f => f.atoms.length === best);
}

/** Bind a freshly read single-point energy to the geometry it was computed at.
 *  ORCA's order is geometry → energy, so the frame waiting for it is the last
 *  one that still has none; frames that already carry an energy are stepped
 *  over, never overwritten. */
function bindEnergy(frames: OrcaFrame[], energy: number): void {
    for (let k = frames.length - 1; k >= 0 && k >= frames.length - MAX_FRAMES_PER_ENERGY; k--) {
        if (frames[k].energy == null) frames[k].energy = energy;
    }
}

/** `     6:    1637.64 cm**-1` — one line of the frequency list. Imaginary modes
 *  print as negative numbers, so the sign is preserved. */
const FREQ_LINE_RE = /^(\d+):\s+(-?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+cm\*\*-1$/;

/** One row of the `|Geometry convergence|` table. Only the value column is used;
 *  the tolerance column is read from the run's own `Convergence Tolerances`
 *  block (see `optThresholds`), which is what the panel draws. */
const CONV_ROW_RE = /^(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+(\S+)(?:\s+(\S+))?/;

/** `TolMAXG  ....  3.0000e-04` — one line of the tolerances block. */
const TOL_LINE_RE = /Tol(MAXG|RMSG|MAXD|RMSD)\s+\.+\s*(\S+)/;

/** `|  1> ! B3LYP def2-SVP Opt freq` — one echoed input line. */
const INPUT_ECHO_RE = /^\|\s*\d+>\s?(.*)$/;

/**
 * Read the echoed input file, which gives the three things the rest of the
 * output only states indirectly: the simple-input keywords (the route panel),
 * and the charge/multiplicity of the run (the info panel and bond perception).
 *
 * ORCA 4/5/6 all echo the input as `|  n> <text>`; the echo repeats for every
 * job in a compound run, so the values are simply kept up to date as they are
 * met. Returns whether the line belonged to the echo.
 */
interface EchoState {
    keywords: string[];
    raw: string;
    charge?: number;
    multiplicity?: number;
}

function readEchoLine(line: string, state: EchoState): boolean {
    const m = INPUT_ECHO_RE.exec(line);
    if (!m) return false;
    const text = m[1].trim();

    if (text.charAt(0) === '!') {
        // The simple input line; wrapped `!` lines are joined in order.
        state.raw = state.raw ? state.raw + ' ' + text : text;
        for (const tok of text.slice(1).trim().split(/\s+/)) {
            if (tok) state.keywords.push(tok);
        }
        return true;
    }

    // `* xyz <charge> <multiplicity>` and `* xyzfile <charge> <mult> <file>`.
    if (/^\*\s*(?:xyz|xyzfile)\b/i.test(text)) {
        const tokens = text.replace(/^\*/, '').trim().split(/\s+/);
        const c = parseInt(tokens[1], 10);
        const s = parseInt(tokens[2], 10);
        if (!isNaN(c)) state.charge = c;
        if (!isNaN(s)) state.multiplicity = s;
    }
    return true;
}

/** The keyword list as the route panel expects it. ORCA's simple input has no
 *  `key=value` option syntax (options live in `%block` directives), so every
 *  token is a keyword with no options. */
function buildOrcaRoute(state: EchoState): RouteSection | undefined {
    if (state.keywords.length === 0) return undefined;
    const keywords: RouteKeyword[] = state.keywords.map(k => ({ name: k, options: [] }));
    const hasOpt = state.keywords.some(k => {
        const n = k.toLowerCase();
        // `opt`, `OptTS`, `Opt Freq`, and `GOAT` (ORCA's global optimizer).
        return n.indexOf('opt') === 0 || n === 'goat';
    });
    return { raw: state.raw, keywords, hasOpt };
}

/**
 * Parse one `VIBRATIONAL FREQUENCIES` block into a map keyed by ORCA's own
 * (0-based) mode number.
 *
 * ORCA lists all 3N modes, the first 3N−6 of which are the translations and
 * rotations it prints as exactly `0.00`. They are kept here and dropped later by
 * `buildModes`, which can also look at their displacement vectors.
 */
function parseFrequencies(lines: string[], start: number): { freqs: Map<number, number>; end: number } {
    const freqs = new Map<number, number>();
    let i = start + 1;
    while (i < lines.length) {
        const t = lines[i].trim();
        const m = FREQ_LINE_RE.exec(t);
        if (m) {
            freqs.set(parseInt(m[1], 10), parseFloat(m[2]));
            i++;
            continue;
        }
        // Before the first entry: the underline, the blank line and the
        // `Scaling factor for frequencies = ... (already applied!)` banner. The
        // walk is bounded, so a section that was suppressed cannot make this
        // consume the rest of the file.
        if (freqs.size === 0) {
            if (i - start > 40) return { freqs, end: start + 1 };
            i++;
            continue;
        }
        break;
    }
    return { freqs, end: i };
}

/** How many lines the `NORMAL MODES` preamble may occupy before the parser
 *  decides the matrix is missing and hands the block back untouched. */
const NORMAL_MODES_PREAMBLE_MAX = 40;

/**
 * Parse one `NORMAL MODES` matrix into one vector of degrees of freedom per
 * mode (`dofs[mode][dof]`), in the order ORCA prints them.
 *
 * The block is a set of column panels, each introduced by a header line that
 * carries the mode numbers of that panel:
 *
 * ```
 *                   0          1          2
 *       0       0.000000   0.000000   0.000000
 *       1       0.000000   0.000000   0.069601
 *                   …          (next panel: modes 3, 4, 5)
 * ```
 *
 * Row numbers are degrees of freedom (0 … 3N−1), so atom a's vector is rows
 * 3a, 3a+1, 3a+2. Grouping them into atoms happens once the run's atom count is
 * settled (`vectorsForAtoms`), not here: a block that was half written must not
 * decide how many atoms the trajectory has.
 *
 * A panel header is integers only, while a data row always has a decimal point
 * — that is what tells the two apart. When no panel starts within the preamble
 * budget, the matrix is taken to be absent and `end` is left at `start + 1`, so
 * the outer scan continues normally and the sections that follow (IR SPECTRUM,
 * thermochemistry) are still read.
 */
function parseNormalModes(lines: string[], start: number): { dofs: number[][]; end: number } {
    interface Panel { modes: number[]; rows: { dof: number; values: number[] }[] }
    const panels: Panel[] = [];
    let current: Panel | null = null;
    let i = start + 1;
    const preambleEnd = start + 1 + NORMAL_MODES_PREAMBLE_MAX;

    const closePanel = () => {
        if (current && current.rows.length > 0) panels.push(current);
        current = null;
    };

    while (i < lines.length) {
        // The preamble budget bounds the search for the first panel header, so a
        // matrix that is not there cannot make this walk into the rest of the
        // file.
        if (current === null && panels.length === 0 && i > preambleEnd) {
            return { dofs: [], end: start + 1 };
        }
        const t = lines[i].trim();
        if (t === '' || isDashLine(t)) {
            closePanel();
            i++;
            continue;
        }
        const parts = t.split(/\s+/);
        if (/^\d+$/.test(parts[0])) {
            let allInts = true;
            for (const p of parts) { if (!/^\d+$/.test(p)) { allInts = false; break; } }
            if (allInts) {
                closePanel();
                current = { modes: parts.map(p => parseInt(p, 10)), rows: [] };
                i++;
                continue;
            }
            if (current) {
                const values: number[] = [];
                for (let k = 1; k < parts.length; k++) {
                    const v = parseFloat(parts[k]);
                    values.push(isFinite(v) ? v : 0);
                }
                current.rows.push({ dof: parseInt(parts[0], 10), values });
                i++;
                continue;
            }
        }
        // A line that is neither a header, a row nor an underline is the
        // explanatory preamble while no panel has started yet.
        closePanel();
        if (panels.length > 0) break;
        i++;
    }
    closePanel();
    if (panels.length === 0) return { dofs: [], end: start + 1 };

    // Scatter the panels into one vector per mode. Values are addressed through
    // each panel's own mode numbers, so a file whose panels are labelled out of
    // order still lands every column on the mode it belongs to.
    const dofByMode: number[][] = [];
    for (const panel of panels) {
        for (const row of panel.rows) {
            for (let c = 0; c < panel.modes.length; c++) {
                const mode = panel.modes[c];
                let dofs = dofByMode[mode];
                if (!dofs) { dofs = []; dofByMode[mode] = dofs; }
                dofs[row.dof] = row.values[c] !== undefined ? row.values[c] : 0;
            }
        }
    }

    return { dofs: dofByMode, end: i };
}

/**
 * Group one mode's degrees of freedom into per-atom vectors.
 *
 * A mode whose DOF count does not cover the structure is left without vectors
 * instead of being cut short, so the panel reports "no displacements for this
 * structure" rather than animating a truncated one.
 */
function vectorsForAtoms(dofs: number[] | undefined, atomCount: number): number[][] {
    if (!dofs || atomCount <= 0 || dofs.length < atomCount * 3) return [];
    const perAtom: number[][] = [];
    for (let a = 0; a < atomCount; a++) {
        perAtom.push([dofs[a * 3] || 0, dofs[a * 3 + 1] || 0, dofs[a * 3 + 2] || 0]);
    }
    return perAtom;
}

/**
 * Turn a frequency list and its displacement matrix into the panel's modes.
 *
 * Two things ORCA does differently from Gaussian are settled here:
 *
 *  * **The rigid-body modes are dropped.** ORCA prints all 3N modes, so every
 *    water frequency list starts with six `0.00 cm**-1` rows that animate
 *    nothing. Gaussian prints only the 3N−6 (3N−5 for a linear molecule)
 *    vibrational ones, and a list that matches Gaussian is what the panel was
 *    built for. A mode is dropped only when its frequency is zero **and** its
 *    displacement vector is zero, so a genuinely soft mode that rounds to `0.00`
 *    survives.
 *  * **Imaginary modes keep their sign.** ORCA prints them as negative
 *    frequencies; the panel already marks a negative frequency with `*`, exactly
 *    as it marks a Gaussian imaginary mode.
 */
function buildModes(
    freqs: Map<number, number>,
    dofsByMode: number[][],
    atomCount: number,
    intensities: Map<number, number>
): NormalMode[] | undefined {
    const modeNumbers: number[] = [];
    freqs.forEach((_v, k) => modeNumbers.push(k));
    if (modeNumbers.length === 0) return undefined;
    modeNumbers.sort((a, b) => a - b);

    const modes: NormalMode[] = [];
    for (const k of modeNumbers) {
        const freq = freqs.get(k) as number;
        const disp = vectorsForAtoms(dofsByMode[k], atomCount);
        let maxMag = 0;
        for (const d of disp) {
            const m = Math.abs(d[0]) + Math.abs(d[1]) + Math.abs(d[2]);
            if (m > maxMag) maxMag = m;
        }
        if (Math.abs(freq) < 1e-6 && maxMag < 1e-9) continue;   // translation/rotation
        modes.push({
            index: modes.length + 1,
            frequency: freq,
            irIntensity: intensities.get(k),
            displacements: disp
        });
    }
    return modes.length > 0 ? modes : undefined;
}

export function parseOrcaOut(content: string): OrcaLogResult {
    const lines = content.split(/\r?\n/);

    const frames: OrcaFrame[] = [];
    /** The Ångström blocks are the primary geometry source; the A.U. blocks are
     *  parsed only when no Ångström block exists at all (a run that suppressed
     *  them with `%output`), because every cycle prints both and they describe
     *  the same geometry. */
    const auFrames: OrcaFrame[] = [];

    const echo: EchoState = { keywords: [], raw: '' };
    let charge: number | undefined;
    let multiplicity: number | undefined;

    /** Label the next geometry block gets, set by the banner that precedes it. */
    let pendingLabel = '';
    let lastEnergy: number | undefined;
    const optSteps: OptStep[] = [];
    let stepCounter = 0;
    const thresholds: OptThresholds = {};

    /** The frequency set being assembled: `freqs` comes from `VIBRATIONAL
     *  FREQUENCIES`, `dofs` from the `NORMAL MODES` matrix after it and
     *  `intensities` from the `IR SPECTRUM` after that. A run with two frequency
     *  jobs replaces the whole set with the last one — that is the section which
     *  belongs to the structure the file ends on, and the panel has one list. */
    let pendingFreqs: Map<number, number> | null = null;
    let pendingDofs: number[][] | null = null;
    const pendingIntensities = new Map<number, number>();

    /** Energies are bound to whichever geometry list is in use. Ångström blocks
     *  come before their energy, so `frames` is already populated by the time
     *  the first energy is printed unless the run suppressed them. */
    const liveFrames = () => (frames.length > 0 ? frames : auFrames);

    let i = 0;
    while (i < lines.length) {
        const line = lines[i];

        // --- Echoed input file (keywords, charge, multiplicity) --------------
        // The char check keeps the regex off the overwhelming majority of lines.
        if (line.charCodeAt(0) === 124 /* '|' */ && readEchoLine(line, echo)) {
            if (echo.charge !== undefined) charge = echo.charge;
            if (echo.multiplicity !== undefined) multiplicity = echo.multiplicity;
            i++;
            continue;
        }

        // --- Banners that name the geometry following them -------------------
        // `*  GEOMETRY OPTIMIZATION CYCLE   1  *`
        if (line.indexOf('GEOMETRY OPTIMIZATION CYCLE') >= 0) {
            const m = line.match(/CYCLE\s+(\d+)/);
            pendingLabel = m ? 'Cycle ' + m[1] : 'Cycle';
            i++;
            continue;
        }
        // `*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***` — the
        // geometry the run actually converged to, printed once after HURRAY.
        if (line.indexOf('FINAL ENERGY EVALUATION AT THE STATIONARY POINT') >= 0) {
            pendingLabel = 'Final';
            i++;
            continue;
        }

        // --- Geometry blocks --------------------------------------------------
        if (line.indexOf('CARTESIAN COORDINATES') >= 0) {
            if (line.indexOf('ANGSTROEM') >= 0 || line.indexOf('ANGSTROM') >= 0) {
                const block = parseAngstromBlock(lines, i);
                if (block.atoms.length > 0) {
                    const label = pendingLabel || 'Geometry ' + (frames.length + 1);
                    frames.push({
                        atoms: block.atoms,
                        bonds: [],
                        title: label,
                        hasExplicitBonds: false,
                        stepLabel: label,
                        charge,
                        multiplicity
                    });
                }
                i = block.end;
                continue;
            }
            if (line.indexOf('A.U.') >= 0) {
                if (frames.length > 0) {
                    // The Ångström half of this same geometry was already read.
                    i = skipCoordinateBlock(lines, i);
                    continue;
                }
                const block = parseAuBlock(lines, i);
                if (block.atoms.length > 0) {
                    const label = pendingLabel || 'Geometry ' + (auFrames.length + 1);
                    auFrames.push({
                        atoms: block.atoms,
                        bonds: [],
                        title: label,
                        hasExplicitBonds: false,
                        stepLabel: label,
                        charge,
                        multiplicity
                    });
                }
                i = block.end;
                continue;
            }
        }

        // --- Single-point energies -------------------------------------------
        // `FINAL SINGLE POINT ENERGY       -76.321097057517`
        if (line.indexOf('FINAL SINGLE POINT ENERGY') >= 0) {
            const m = line.match(/(-?\d+\.\d+(?:[EeDd][-+]?\d+)?)/);
            if (m) {
                lastEnergy = parseFloat(m[1].replace(/[Dd]/, 'E'));
                bindEnergy(liveFrames(), lastEnergy);
            }
            i++;
            continue;
        }
        // Correlated methods restate the same number as `Total Energy : ... Eh`;
        // it is only ever used when no `FINAL SINGLE POINT ENERGY` was printed.
        if (line.indexOf('Total Energy') >= 0) {
            const m = line.match(/Total Energy\s*:\s*(-?\d+\.\d+)\s*Eh/);
            if (m) {
                lastEnergy = parseFloat(m[1]);
                bindEnergy(liveFrames(), lastEnergy);
            }
            i++;
            continue;
        }

        // --- Optimization convergence table ----------------------------------
        // `----------------------|Geometry convergence|-------------------------`
        // The same five rows are printed by ORCA 4, 5 and 6. `Energy change` is a
        // *difference*, not an energy, so it is not stored as one; the step's
        // energy is the `FINAL SINGLE POINT ENERGY` printed for its geometry.
        if (line.indexOf('Geometry convergence') >= 0) {
            const step: OptStep = { step: stepCounter + 1, energy: lastEnergy };
            let j = i + 1;
            let matched = 0;
            let guard = 0;
            while (j < lines.length && matched < 5 && guard++ < 24) {
                const t = lines[j].trim();
                const m = CONV_ROW_RE.exec(t);
                if (m) {
                    const value = parseFloat(m[2].replace(/[Dd]/, 'E'));
                    if (isFinite(value)) {
                        if (m[1] === 'MAX gradient') step.maxForce = value;
                        else if (m[1] === 'RMS gradient') step.rmsForce = value;
                        else if (m[1] === 'MAX step') step.maxDisplacement = value;
                        else if (m[1] === 'RMS step') step.rmsDisplacement = value;
                        // The table's own Tolerance column is the same number the
                        // `Convergence Tolerances` block states. It is used only
                        // when that block was not printed (older versions, or a
                        // run that suppressed it), so a chart always draws the
                        // targets the run was actually held to.
                        const tol = parseFloat((m[3] || '').replace(/[Dd]/, 'E'));
                        if (isFinite(tol)) {
                            if (m[1] === 'MAX gradient' && thresholds.maxForce === undefined) thresholds.maxForce = tol;
                            else if (m[1] === 'RMS gradient' && thresholds.rmsForce === undefined) thresholds.rmsForce = tol;
                            else if (m[1] === 'MAX step' && thresholds.maxDisplacement === undefined) thresholds.maxDisplacement = tol;
                            else if (m[1] === 'RMS step' && thresholds.rmsDisplacement === undefined) thresholds.rmsDisplacement = tol;
                        }
                        matched++;
                    }
                } else if (matched > 0 && (t === '' || isDashLine(t))) {
                    break;
                }
                j++;
            }
            if (matched > 0) {
                optSteps.push(step);
                stepCounter++;
                i = j;
            } else {
                // Not a real table (a mention of one, or a print-suppressed run):
                // step over the header only, never over the lines that follow it.
                i++;
            }
            continue;
        }

        // --- The run's own convergence targets --------------------------------
        // `Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr`
        if (line.indexOf('TolMAXG') >= 0 || line.indexOf('TolRMSG') >= 0 ||
            line.indexOf('TolMAXD') >= 0 || line.indexOf('TolRMSD') >= 0) {
            const m = TOL_LINE_RE.exec(line);
            if (m) {
                const v = parseFloat(m[2].replace(/[Dd]/, 'E'));
                if (isFinite(v)) {
                    if (m[1] === 'MAXG') thresholds.maxForce = v;
                    else if (m[1] === 'RMSG') thresholds.rmsForce = v;
                    else if (m[1] === 'MAXD') thresholds.maxDisplacement = v;
                    else if (m[1] === 'RMSD') thresholds.rmsDisplacement = v;
                }
            }
            i++;
            continue;
        }

        // --- Vibrational frequencies -----------------------------------------
        if (line.indexOf('VIBRATIONAL FREQUENCIES') >= 0) {
            // A second frequency job replaces the first: the structure the file
            // ends on is the one its last frequencies belong to.
            const res = parseFrequencies(lines, i);
            pendingFreqs = res.freqs.size > 0 ? res.freqs : null;
            pendingDofs = null;
            pendingIntensities.clear();
            i = res.end;
            continue;
        }

        // --- Normal mode displacement matrix ---------------------------------
        if (line.indexOf('NORMAL MODES') >= 0) {
            const res = parseNormalModes(lines, i);
            pendingDofs = res.dofs;
            i = res.end;
            continue;
        }

        // --- IR intensities ---------------------------------------------------
        // `  6:   1637.64   0.010942   55.30  0.002085  (...)`
        // `Int` is in km/mol — the same unit as Gaussian's `IR Inten`, so the
        // field means the same thing whichever program wrote the file.
        if (line.indexOf('IR SPECTRUM') >= 0) {
            // The rows sit under a two-line table header and a rule:
            //     Mode   freq       eps      Int      T**2   TX   TY   TZ
            //            cm**-1   L/(mol*cm) km/mol    a.u.
            //     -----------------------------------------------------------------
            //       6:   1637.64   0.010942   55.30  0.002085  (...)
            // so the first row is located within a small window, and only the
            // consecutive rows after it are read. Finding nothing in that window
            // leaves the scan where it was, so a mention of "IR SPECTRUM" with no
            // table under it cannot swallow the lines that follow.
            let j = i + 1;
            let found = -1;
            while (j < lines.length && j <= i + 12 && found < 0) {
                if (/^\s*\d+:/.test(lines[j])) found = j;
                j++;
            }
            if (found < 0) { i++; continue; }
            j = found;
            let guard = 0;
            while (j < lines.length && guard++ < 400) {
                const t = lines[j].trim();
                const m = /^(\d+):\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/.exec(t);
                if (!m) break;
                pendingIntensities.set(parseInt(m[1], 10), parseFloat(m[4]));
                j++;
            }
            i = j;
            continue;
        }

        i++;
    }

    // --- Assemble the result -------------------------------------------------
    // The A.U. fallback is used only when no Ångström geometry was found; a
    // half-written block at either end of the file is dropped from whichever
    // list is used.
    const source = dropOddSizedFrames(frames.length > 0 ? frames : auFrames);

    // A coordinate block comes before the `SCF SETTINGS` section that states the
    // charge, so a frame created before that point is filled in here. A frame
    // that already carries a value keeps it: a compound run can change the
    // charge between jobs and each job's frames must keep their own.
    if (charge !== undefined || multiplicity !== undefined) {
        for (const f of source) {
            if (f.charge === undefined && charge !== undefined) f.charge = charge;
            if (f.multiplicity === undefined && multiplicity !== undefined) f.multiplicity = multiplicity;
        }
    }

    // Labels are composed once, when every fact they mention is known — the same
    // `(Charge=0 Mult=1)` suffix a Gaussian log's `Step k` frames carry, so a
    // trajectory reads the same whichever program wrote it.
    const suffix = (charge !== undefined && multiplicity !== undefined)
        ? ' (Charge=' + charge + ' Mult=' + multiplicity + ')'
        : '';
    for (const f of source) {
        f.stepLabel = f.stepLabel + suffix;
        f.title = f.stepLabel;
    }

    // The frequency set was assembled while scanning; it is built once here, when
    // the IR intensities that follow it are known as well.
    const modeAtomCount = source.length > 0 ? source[source.length - 1].atoms.length : 0;
    const modes = pendingFreqs
        ? buildModes(pendingFreqs, pendingDofs || [], modeAtomCount, pendingIntensities)
        : undefined;

    const route = buildOrcaRoute(echo);
    const hasThresholds = thresholds.maxForce !== undefined || thresholds.rmsForce !== undefined ||
        thresholds.maxDisplacement !== undefined || thresholds.rmsDisplacement !== undefined;

    return {
        frames: source,
        title: echo.raw,
        optSteps: optSteps.length > 0 ? optSteps : undefined,
        normalModes: modes,
        routes: route ? [route] : undefined,
        optThresholds: hasThresholds ? thresholds : undefined,
        charge,
        multiplicity
    };
}
