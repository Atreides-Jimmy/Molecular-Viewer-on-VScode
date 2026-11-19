import { Atom, Bond, MolecularData, OptStep, NormalMode, SolventCavity, OrbitalLevels, OrbitalSpinLevels, ScanPoint, ScanResult, IrcResult, IrcPoint, IrcRun, IrcDirection, LogPathways } from '../types';

const ATOMIC_NUMBER_MAP: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
};

export interface LogFrame {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    stepLabel: string;
    charge?: number;
    multiplicity?: number;
    /** The cell vectors (rows a, b, c, in Angstrom) a trajectory frame states. */
    cell?: number[][];
    /** The time the trajectory states, in the format's own unit. */
    time?: number;
    /** The integration step the trajectory states. */
    step?: number;
    /** Single-point energy (Hartree) of the optimization step this frame is the
     *  geometry of, filled by the parser for opt logs so the webview can both
     *  label the frame and highlight its point in the convergence charts.
     *  [2026-09-12 | Atreides-Jimmy] Added for the convergence-chart highlight
     *  (feature 2) and value tooltip (feature 3). */
    energy?: number;
    /** The optimizer step (1-based, in the file's own series) whose geometry
     *  this frame is, when the job that printed it runs an optimization.
     *
     *  The frame INDEX is not that number: a log that prints the input
     *  orientation as a block of its own repeats the first step's geometry as a
     *  frame of its own (the input pose, then the standard pose), so every
     *  later frame sits one slot past its step. The convergence charts follow
     *  this link instead of the index. Absent on frames no optimizer produced —
     *  the converged structure an `opt freq` run hands to its `Freq` job, an
     *  IRC path point, a trajectory frame.
     *  [2026-11-12 | Atreides-Jimmy] Added for the step guide of feature 2. */
    optStep?: number;
    /** The scan point this geometry belongs to (1-based), for a relaxed scan.
     *  Every step of a point carries it, so the scan panel can mark the point
     *  whose relaxation is on screen and not only its converged end.
     *  [2026-11-12 | Atreides-Jimmy] Added with the Scan panel. */
    scanPoint?: number;
    /** The optimizer step within that scan point (1-based). */
    scanStep?: number;
}

export interface GaussianLogResult {
    frames: LogFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    routes?: RouteSection[];
    /** The population analyses of the final job, one per scheme the file
     *  printed (Mulliken, ESP, APT, …). Absent when the job states none. */
    charges?: LogChargeSet[];
    /** One block of headline results per job in the file. */
    summary?: LogSummary[];
    /** The last eigenvalue dump of the file: what the level diagram draws.
     *  A log states energies and no coefficients, so this can never become a
     *  surface — only a formatted checkpoint can. */
    orbitals?: OrbitalLevels;
    /** The solvation cavity the file describes (type, solvent, area, volume).
     *  A log names no per-sphere radii; those are in the checkpoint. */
    cavity?: SolventCavity;
    /** A relaxed scan and/or an IRC path the file holds. Both are sequences of
     *  structures along one coordinate rather than an optimizer series, so the
     *  page shows them in place of the convergence charts.
     *  [2026-11-12 | Atreides-Jimmy] Added with the Scan and IRC panels. */
    pathways?: LogPathways;
}

/**
 * One per-atom charge table a log prints: a population-analysis scheme read
 * from the block's own rows.
 *
 * A block that lists fewer atoms than the structure — Gaussian's "Mulliken
 * charges with hydrogens summed into heavy atoms" is the common one — is kept
 * as it stands: `values` holds `null` at the atoms it says nothing about, so a
 * view can paint exactly what the file stated and never fill a gap with a
 * guessed zero.
 */
export interface LogChargeSet {
    /** The block's own name, as printed and without the trailing colon. */
    label: string;
    /** Charge per atom in electron units, by 0-based atom index. */
    values: (number | null)[];
    /** The element of each atom the values were read for, in the same order.
     *  A page only paints a value onto an atom that is still the atom it was
     *  read for: an edit, or a frame of another job, changes this list. */
    elements: string[];
    /** How many atoms the block states a value for. */
    covered: number;
    /** True when every atom of that frame has a value. */
    complete: boolean;
    /** 0-based index of the frame whose geometry the analysis belongs to. */
    frame: number;
}

/** One energy a job printed, under the file's own label (`E(RPM7)`). */
export interface LogEnergy {
    label: string;
    /** Hartree, as printed. */
    value: number;
    /** Which SCF of the job printed it (1-based). An optimisation prints one
     *  per step; the last is the converged geometry. */
    scf?: number;
}

/** The dipole moment a job states, in Debye. */
export interface LogDipole {
    x: number;
    y: number;
    z: number;
    total: number;
}

/**
 * The headline results of one job in a log: the block a reader checks before
 * looking at the structure — what energy, how well converged, what dipole and
 * point group, what it cost.
 *
 * Every field is what the FILE states. Nothing here is derived, defaulted or
 * converted (an energy stays in Hartree, a dipole in Debye, a time in the
 * program's own words), so a summary can never show a number the file does not
 * contain; a field the file omits stays absent and is not shown at all.
 */
export interface LogSummary {
    /** 1-based job number in the file (a `--Link1--` file holds several). */
    job: number;
    nAtoms?: number;
    nBasis?: number;
    charge?: number;
    multiplicity?: number;
    /** Every SCF energy the job printed, last value per label, in printed order. */
    energies: LogEnergy[];
    /** Cycles the last SCF took. */
    scfCycles?: number;
    /** The last 'RMS Gradient Norm' the job printed. */
    rmsGradient?: number;
    /** The label that line carried, when the file words it differently. */
    rmsGradientLabel?: string;
    dipole?: LogDipole;
    pointGroup?: string;
    /** 'Job cpu time' / 'Elapsed time', as the file words them. */
    cpuTime?: string;
    elapsedTime?: string;
    termination?: 'normal' | 'error';
    /**
     * The SCF (1-based, within this job) each geometry-level value was last
     * printed at. They are all the same for a job that finished — the converged
     * geometry — and differ when the file only printed a value for an earlier
     * step (a run that was cut short), which the panel marks instead of
     * presenting two geometries as one.
     */
    scfOf?: {
        energy?: number;
        rmsGradient?: number;
        dipole?: number;
    };
}

/**
 * How far a table header may be scanned for its dashed rules. Gaussian's
 * orientation header is five lines; a limit keeps a damaged file whose dashes
 * are missing from swallowing every later frame, convergence table and
 * frequency section to EOF.
 */
const MAX_DASHED_HEADER_LINES = 24;

function skipDashedLines(lines: string[], startIdx: number, count: number): number {
    let i = startIdx;
    let skipped = 0;
    const limit = Math.min(lines.length, startIdx + MAX_DASHED_HEADER_LINES);
    while (i < limit && skipped < count) {
        const text = lines[i].trim();
        if (text.includes('---')) {
            skipped++;
            i++;
            continue;
        }
        // The header is over as soon as an atom row appears, even when one or
        // both of its dashed rules is missing: counting dashes only would scan
        // past the geometry and lose the frame (and everything after it).
        if (/^\d/.test(text)) break;
        i++;
    }
    return i;
}

/** The shape of a route card's first line: '#' plus the keywords, with or
 *  without a print level. Gaussian prints `#p B3LYP/6-31G(d) opt freq` and
 *  `#B3LYP/6-31G(d) opt freq` for the same request, so the '#' may be followed
 *  by whitespace, by nothing, by a single-letter print level (#p/#P/#t/#n,
 *  case-insensitive), or directly by the first keyword.
 *
 *  This is a SHAPE test only — '# OF POINTS ALONG THE PATH = 1' (an IRC
 *  banner) has it too. What separates a card from such text is the dash box
 *  it is printed in, and the scan below therefore accepts a boxed line of this
 *  shape, and an unboxed one only when it states a method/basis pair (a '/',
 *  which every printed card carries) or is closed by a rule. */
const ROUTE_START_RE = /^#(?:\s|[A-Za-z]\b|$|[A-Za-z0-9][^\s]*\/)/;
/** Pure dashed separator line (Gaussian prints 50-70 dashes). Route cards are
 *  always enclosed between two of these, single-line or wrapped. */
const DASH_RE = /^\s*-{5,}\s*$/;

export interface RouteKeyword {
    name: string;       // keyword as written (case preserved), e.g. 'opt', 'SCRF'
    options: string[];  // option list as written, e.g. ['SMD','Solvent=Acetone']
}

export interface RouteSection {
    raw: string;             // route text joined into one line
    keywords: RouteKeyword[];
    hasOpt: boolean;         // the section requests an optimization
}



/**
 * Join route lines into a single string.  Handles Gaussian column-80 word
 * wrapping: if the previous line ends with an alphanumeric character and the
 * next line starts with one, they are joined without a space (mid-word wrap).
 */
function joinRouteLines(lines: string[]): string {
    if (lines.length === 0) return '';
    const trimmed = lines.map(l => l.trimStart());
    let result = trimmed[0];
    for (let i = 1; i < trimmed.length; i++) {
        const cur = trimmed[i];
        if (/\w$/.test(result) && /^\w/.test(cur)) {
            result += cur;
        } else {
            result += ' ' + cur;
        }
    }
    return result.trim();
}

/**
 * Structured parse of one route section, used both for the optimization gate
 * and for the webview's Route panel. Tokenization is paren- and quote-aware
 * at the top level, so scrf=(SMD, Solvent=Acetone) and external="python3
 * /usr/bin/xtb --alpb water" each stay a single token, and a bare '=(calcall)'
 * continuation token (route wrapping can split 'opt' from its option list) is
 * folded into the preceding keyword's options. Full-width characters are
 * normalized for parsing; the raw text is kept for display.
 */
function parseRouteSection(routeText: string): RouteSection {
    const r = routeText
        .replace(/（/g, '(').replace(/）/g, ')')
        .replace(/＝/g, '=').replace(/，/g, ',');
    // Top-level tokens: whitespace splits only outside parens/quotes
    const tokens: string[] = [];
    let cur = '';
    let depth = 0;
    let inQuote = false;
    for (let c = 0; c < r.length; c++) {
        const ch = r[c];
        if (inQuote) { cur += ch; if (ch === '"') inQuote = false; continue; }
        if (ch === '"') { inQuote = true; cur += ch; continue; }
        if (ch === '(') depth++;
        else if (ch === ')') { if (depth > 0) depth--; }
        else if (depth === 0 && /\s/.test(ch)) { if (cur !== '') { tokens.push(cur); cur = ''; } continue; }
        cur += ch;
    }
    if (cur !== '') tokens.push(cur);
    // Split tokens into keyword entries. The '=' separating a keyword from
    // its value is the first '=' OUTSIDE parentheses/quotes — iop(5/17=17)
    // is a bare-paren keyword, not 'iop(5/17' = '17)'.
    const eqOutside = (s: string): number => {
        let d = 0;
        let q = false;
        for (let c = 0; c < s.length; c++) {
            const ch = s[c];
            if (q) { if (ch === '"') q = false; continue; }
            if (ch === '"') q = true;
            else if (ch === '(') d++;
            else if (ch === ')') { if (d > 0) d--; }
            else if (ch === '=' && d === 0) return c;
        }
        return -1;
    };
    const keywords: RouteKeyword[] = [];
    const pushValue = (kw: RouteKeyword | null, value: string) => {
        let v = value.trim();
        if (v.startsWith('(') && v.endsWith(')')) v = v.slice(1, -1);
        const opts = v.split(',').map(s => s.trim()).filter(s => s !== '');
        if (kw && opts.length > 0) kw.options.push(...opts);
    };
    for (const t of tokens) {
        if (t[0] === '=') {
            // Wrapped continuation: '=(calcall)' belongs to the previous keyword
            pushValue(keywords.length > 0 ? keywords[keywords.length - 1] : null, t.slice(1));
            continue;
        }
        const eq = eqOutside(t);
        if (eq > 0) {
            const kw: RouteKeyword = { name: t.slice(0, eq), options: [] };
            pushValue(kw, t.slice(eq + 1));
            keywords.push(kw);
            continue;
        }
        const pi = t.indexOf('(');
        if (pi > 0 && t.endsWith(')')) {
            // Only split when the part before '(' is a simple alphabetic keyword
            // (e.g. iop, scrf, TD).  Basis-set names like 6-31G(d,p) or
            // method/basis combos like B3LYP/6-31G(d,p) contain digits, '-',
            // '+', or '/' and must be kept as a single token.
            const namePart = t.slice(0, pi);
            if (/^[A-Za-z]+$/.test(namePart)) {
                keywords.push({
                    name: namePart,
                    options: t.slice(pi + 1, -1).split(',').map(s => s.trim()).filter(s => s !== '')
                });
                continue;
            }
        }
        keywords.push({ name: t, options: [] });
    }
    // Optimization request — the same semantics as the panel-visibility gate:
    // whole-token opt / optimization (any spelling, with or without options),
    // or calcall (standalone non-standard form, or as a keyword option such as
    // freq=calcall which implies an optimization).
    const hasOpt = keywords.some(k => {
        const n = k.name.toLowerCase();
        return n === 'opt' || n === 'optimization' || n === 'calcall' ||
            k.options.some(o => o.toLowerCase() === 'calcall');
    });
    return { raw: routeText.replace(/\s+/g, ' ').trim(), keywords, hasOpt };
}

/**
 * Parse ONE harmonic-frequencies section starting AT the header line.
 * Returns the modes found and the line index where the section ends so the
 * caller's single-pass scan can resume there. Stops early at another section
 * header, so each section is parsed independently and the last one wins
 * (matching the previous "scan for the last header, then parse" behavior).
 */
function parseModesFromHeader(lines: string[], headerIdx: number): { modes: NormalMode[]; endIdx: number } {
    const modes: NormalMode[] = [];
    let i = headerIdx + 1;

    // Skip blank lines and column header text until we reach the first frequency block
    while (i < lines.length) {
        const raw = lines[i];
        // A new section header ends this one (caller re-checks and parses there)
        if (raw.includes('Harmonic frequencies') && raw.includes('cm**-1')) break;

        const line = raw.trim();

        // Detect end of normal modes section
        if (line === '' && modes.length > 0) {
            // Check if next non-blank line is another block or end
            let j = i;
            while (j < lines.length && lines[j].trim() === '') j++;
            if (j >= lines.length) break;
            const nextLine = lines[j].trim();
            // If next line looks like another block header (starts with numbers), continue
            if (!/^\d/.test(nextLine) && !nextLine.includes('Frequencies')) break;
        }

        // Block header line: "1                      2                      3" (1-3 columns)
        const colMatch = line.match(/^(\d+)(?:\s+(\d+))?(?:\s+(\d+))?\s*$/);
        if (!colMatch) {
            i++;
            continue;
        }

        // Next line: symmetry labels "A                      A                      A"
        const symLine = lines[i + 1] || '';
        const symParts = symLine.trim().split(/\s+/).filter(function(s){return s.length>0});

        // Frequencies line: "Frequencies --   1583.7822              3674.5814              3795.6092"
        const freqLine = lines[i + 2] || '';
        if (!freqLine.includes('Frequencies')) {
            i++;
            continue;
        }
        const freqParts = freqLine.split(/Frequencies\s+--/)[1] || '';
        const freqs = freqParts.trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));

        if (freqs.length === 0) { i++; continue; }

        // Reduced masses, force constants, IR intensities
        let redMasses: number[] = [];
        let frcConsts: number[] = [];
        let irIntens: number[] = [];
        let atomStartIdx = i + 3;
        for (let k = i + 3; k < Math.min(lines.length, i + 12); k++) {
            const kl = lines[k];
            if (kl.includes('Red. masses')) {
                const parts = kl.split(/Red\. masses\s+--/);
                if (parts[1]) redMasses = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('Frc consts')) {
                const parts = kl.split(/Frc consts\s+--/);
                if (parts[1]) frcConsts = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('IR Inten')) {
                const parts = kl.split(/IR Inten\s+--/);
                if (parts[1]) irIntens = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('Atom') && kl.includes('AN')) {
                atomStartIdx = k + 1;
                break;
            }
        }

        // Parse displacement vectors until blank line or non-atom line
        const nCols = freqs.length;
        const displacementsByMode: number[][][] = Array.from({ length: nCols }, () => []);

        let rowIdx = atomStartIdx;
        while (rowIdx < lines.length) {
            const rl = lines[rowIdx].trim();
            if (rl === '' || rl.includes('---') || rl.includes('Frequencies') || /^(\d+)(?:\s+(\d+))?(?:\s+(\d+))?\s*$/.test(rl)) break;

            const parts = rl.split(/\s+/).map(parseFloat);
            // Format: atomIdx atomicNum X1 Y1 Z1 [X2 Y2 Z2] [X3 Y3 Z3]
            if (parts.length >= 2 + nCols * 3 && !parts.some(p => isNaN(p))) {
                for (let c = 0; c < nCols; c++) {
                    const base = 2 + c * 3;
                    displacementsByMode[c].push([parts[base], parts[base + 1], parts[base + 2]]);
                }
            }
            rowIdx++;
        }

        for (let c = 0; c < freqs.length; c++) {
            const globalIdx = modes.length + 1;
            modes.push({
                index: globalIdx,
                frequency: freqs[c],
                symmetry: symParts[c] || undefined,
                reducedMass: redMasses[c],
                forceConstant: frcConsts[c],
                irIntensity: irIntens[c],
                displacements: displacementsByMode[c] || []
            });
        }

        i = rowIdx;
    }

    return { modes, endIdx: i };
}

/**
 * Attach a freshly parsed single-point energy to the geometry frame(s) that
 * still wait for one, i.e. the frames emitted after the previous energy was
 * consumed by a convergence table.
 *
 * [2026-09-12 | Atreides-Jimmy] Added for the convergence-chart highlight and
 * value tooltip: the webview can only point at "the step this structure
 * belongs to" if the parser records which energy a frame's geometry was
 * computed at.
 *
 * Gaussian prints the geometry of a cycle BEFORE the energy computed at it
 * (`Standard orientation` … `SCF Done` … convergence table), and an
 * `opt=external` run prints `Energy= … NIter=` in the same order, so the energy
 * is bound BACKWARDS to the frames printed since the previous one. The previous
 * form read the energy into the next frame at creation time and therefore
 * produced 139 frames with no energy at all on a real optimisation (the value
 * was consumed by the convergence table before the next geometry arrived).
 *
 * The backward scan is bounded by `MAX_FRAMES_PER_ENERGY` and never reaches
 * past `floor` (the first frame of the current job), so a pathological log
 * cannot turn this into a quadratic scan and a `--Link1--` job cannot inherit
 * its predecessor's energy. Frames that already carry an energy are stepped
 * over, never overwritten.
 */
const MAX_FRAMES_PER_ENERGY = 8;
function assignPendingEnergyToFrames(frames: LogFrame[], energy: number, floor: number): void {
    const stop = Math.max(floor, frames.length - MAX_FRAMES_PER_ENERGY);
    for (let k = frames.length - 1; k >= stop; k--) {
        if (frames[k].energy == null) frames[k].energy = energy;
    }
}

/** A number as these programs print one: `-0.1234`, `.10134D-02`, `1.5e-3`.
 *  Every value the summary reads carries this shape, so the patterns below
 *  state it once. */
const LOG_NUMBER = '([-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[DdEe][-+]?\\d+)?)';

/** Parse one number in a program's own notation (Fortran `D` exponents). */
function parseLogNumber(text: string): number {
    const v = parseFloat(text.replace(/[Dd]/, 'E'));
    return isFinite(v) ? v : NaN;
}

// The patterns the single pass asks every line about, built once: the scan runs
// over every line of a log, and a regex compiled per line is the most expensive
// thing in it. `LOG_NUMBER` is spliced in so a value in `D` notation reads the
// same as one in plain decimals.
const SCF_DONE_RE = new RegExp('SCF Done:\\s+E\\(([^)]+)\\)\\s*=\\s*' + LOG_NUMBER);
const SCF_CYCLES_RE = /after\s+(\d+)\s+cycles/;
const CHARGE_MULT_RE = /Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)/;
const RMS_GRADIENT_RE = new RegExp('RMS Gradient Norm\\s*=\\s*' + LOG_NUMBER);
const POINT_GROUP_RE = /Point Group\s*=\s*(\S+)/;
const FULL_POINT_GROUP_RE = /Full point group\s+(\S+)/;
const DIPOLE_RE = new RegExp('X=\\s*' + LOG_NUMBER + '\\s+Y=\\s*' + LOG_NUMBER +
    '\\s+Z=\\s*' + LOG_NUMBER + '\\s+Tot=\\s*' + LOG_NUMBER);
const NATOMS_RE = /NAtoms=\s*(\d+)/;
const NBUSE_RE = /NBsUse=\s*(\d+)/;
/** 'Job cpu time:  0 days  4 hours 25 minutes 40.3 seconds.' — the tail is
 *  kept as written; only the label and the closing period are dropped. */
const JOB_CPU_RE = /Job cpu time:\s*(.+?)\s*\.?\s*$/;
const ELAPSED_RE = /Elapsed time:\s*(.+?)\s*\.?\s*$/;

// --- A relaxed surface scan and an IRC path -------------------------------
// ` Step number   1 out of a maximum of   57 on scan point     1 out of    32`
// The `on scan point` tail is what makes the table a scan point's step instead
// of a step of one optimisation; a plain Berny run prints the line without it
// and its tables keep the old meaning.
const SCAN_STEP_RE = /Step number\s+(\d+)\s+out of a maximum of\s+(\d+)\s+on scan point\s+(\d+)\s+out of\s+(\d+)/;
/** ` Point Number:   0          Path Number:   1` — an IRC point's own header. */
const IRC_POINT_RE = /Point Number:\s*(\d+)\s+Path Number:\s*(\d+)/;
/** One row of the Berny parameter table: name, definition, value, derivative. */
const SCAN_PARAM_ROW_RE = /^!\s+([A-Za-z]+\d+)\s+([A-Za-z]\([\d,]+\))\s+(-?\d+(?:\.\d+)?)\s+(.*)$/;
/** The derivative column of a scanned parameter. `Scan` alone, not `Scanning`. */
const SCAN_MARK_RE = /^Scan\b/;
/** `Summary of Optimized Potential Surface Scan (add -269.0 to energies):` */
const SCAN_TABLE_BASE_RE = /\(add\s+([-+]?\d+(?:\.\d+)?)\s+to energies\)/i;
/** A group's column header — the point numbers, and nothing else. */
const SCAN_TABLE_HEAD_RE = /^\d+(\s+\d+)+$/;
/** `Eigenvalues --    -0.65139  -0.65102 …` */
const SCAN_TABLE_ENERGY_RE = /^Eigenvalues\s*-+\s+(.*)$/;
/** `R2           3.09780   3.04780 …` — a parameter row of the summary. */
const SCAN_TABLE_ROW_RE = /^([A-Za-z]+\d+)\s+(.*)$/;
/** The summary packs negative values without a blank: `-168.89124-167.99169`. */
const SCAN_TABLE_NUMBER_RE = /-?\d+\.\d+/g;

/** The most atoms a charge table may address. A row that names an index past
 *  this is not a charge table (or not one this reader understands), and the
 *  table is refused rather than used to build an array of that length. */
const MAX_CHARGE_ATOMS = 200000;

/**
 * Read one per-atom charge table starting at its header line.
 *
 * Gaussian prints a population analysis as a short header naming the scheme
 * followed by one row per atom:
 *
 *      Mulliken charges:
 *                    1
 *          1  N   -0.654835
 *          2  H    0.274990
 *
 * The column header between the two is a bare atom number in some versions and
 * a dashed rule in others, and a block that sums hydrogens into heavy atoms
 * lists fewer rows than the structure has atoms. Both are read as they are: the
 * row count is not assumed, and the value of an atom the block omits stays
 * `null`.
 *
 * A row must be exactly (atom number, label, one number) — the shape every
 * simple charge table has — so a table with extra columns (Hirshfeld's spin
 * densities, a fit's per-atom errors) is refused instead of having one of its
 * other columns read as a charge.
 *
 * Returns null when the line is not such a table, so a mention of charges in
 * a settings block ("Polarization charges : Total charges.") costs one failed
 * attempt and nothing else.
 */
function readChargeBlock(
    lines: string[],
    headerIdx: number
): { label: string; values: (number | null)[]; endIdx: number } | null {
    const header = lines[headerIdx].trim();
    if (header.length < 6 || header.length > 120) return null;
    const label = header.replace(/:\s*$/, '').replace(/\s+/g, ' ').trim();
    if (label.length < 5) return null;
    let i = headerIdx + 1;
    // The column header sits between the title and the first row: a bare atom
    // number, a dashed rule, or nothing at all. A block whose rows do not
    // begin within that window is not a table.
    for (let k = 0; k < 4 && i < lines.length; k++, i++) {
        const text = lines[i].trim();
        if (text === '' || /^-{4,}/.test(text) || /^\d+$/.test(text)) continue;
        break;
    }
    const values: (number | null)[] = [];
    let rows = 0;
    let maxIdx = 0;
    while (i < lines.length) {
        const text = lines[i].trim();
        if (text === '' || text.length > 80) break;
        const m = text.match(/^(\d+)\s+([A-Za-z]{1,3}\d*)\s+([-+]?(?:\d+\.?\d*|\.\d+)(?:[DdEe][-+]?\d+)?)$/);
        if (!m) break;
        const idx = parseInt(m[1], 10);
        const value = parseLogNumber(m[3]);
        if (!isFinite(value) || idx < 1 || idx > MAX_CHARGE_ATOMS) break;
        if (idx > maxIdx) {
            maxIdx = idx;
            while (values.length < maxIdx) values.push(null);
        }
        values[idx - 1] = value;
        rows++;
        i++;
    }
    if (rows === 0) return null;
    return { label, values, endIdx: i };
}

/** Whether a job summary carries anything worth showing. A job whose only
 *  record is its charge/multiplicity line is not a result block. */
function summaryHasData(s: LogSummary): boolean {
    return s.energies.length > 0 || s.rmsGradient != null || !!s.dipole ||
        !!s.pointGroup || !!s.cpuTime || !!s.elapsedTime || !!s.termination ||
        s.nAtoms != null || s.nBasis != null;
}

/** A charge table before the frame it belongs to is known. */
interface ChargeDraft {
    label: string;
    values: (number | null)[];
    /** Frames on the list when the table ended — the geometry it describes. */
    anchor: number;
    /** Which job printed it (1-based). Two jobs of a `--Link1--` file can print
     *  the same scheme for different structures, and they are both kept. */
    job: number;
}

/** The charge tables of the file, in printed order. */
interface ChargeState {
    drafts: ChargeDraft[];
    /** The job being read, 1-based: a route card opens a new one. */
    job: number;
}

/** The headline results of every job in the file, and the one being filled. */
interface SummaryState {
    list: LogSummary[];
    cur: LogSummary | null;
}

/** The job being read, created (and numbered) on its first recorded value. */
function summaryOf(st: SummaryState): LogSummary {
    if (!st.cur) {
        st.cur = { job: st.list.length + 1, energies: [] };
        st.list.push(st.cur);
    }
    return st.cur;
}

/** Record the SCF (1-based) a geometry-level value was last printed at.
 *
 *  A summary keeps the last value of each field, and for a finished
 *  optimisation that is one geometry: the converged one. It is not always so —
 *  a run that was cut short prints no analysis for its final step, so its
 *  dipole is the first step's while its energy is the last one's. Stamping the
 *  SCF lets the panel say which values do not belong to the final geometry
 *  instead of presenting two structures as one. */
function noteScf(st: SummaryState, field: 'energy' | 'rmsGradient' | 'dipole', scf: number): void {
    const s = summaryOf(st);
    if (!s.scfOf) s.scfOf = {};
    s.scfOf[field] = scf;
}

/** Record one SCF energy under the file's own label, keeping the LAST value
 *  printed for it: an optimisation prints one per step, and the converged one
 *  is what a summary is for. `scf` is the number of that SCF in the job. */
function noteEnergy(st: SummaryState, label: string, value: number, scf: number): void {
    const s = summaryOf(st);
    noteScf(st, 'energy', scf);
    for (let k = 0; k < s.energies.length; k++) {
        if (s.energies[k].label === label) {
            s.energies[k].value = value;
            s.energies[k].scf = scf;
            return;
        }
    }
    s.energies.push({ label, value, scf });
}

/** Record one charge table, replacing an earlier table of the same scheme in
 *  the same job: a population analysis printed at every step keeps the last
 *  one, which is the geometry the run ended on. A table of ANOTHER job is a
 *  different structure's analysis and is kept, not replaced — a `--Link1--`
 *  file that analysed both jobs used to lose every table but the last. */
function noteChargeTable(st: ChargeState, label: string, values: (number | null)[], anchor: number): void {
    for (let k = 0; k < st.drafts.length; k++) {
        if (st.drafts[k].label === label && st.drafts[k].job === st.job) {
            st.drafts[k] = { label: label, values: values, anchor: anchor, job: st.job };
            return;
        }
    }
    st.drafts.push({ label: label, values: values, anchor: anchor, job: st.job });
}

/**
 * One eigenvalue dump: the two spin channels of one SCF.
 *
 * Gaussian prints every occupied alpha eigenvalue, then every virtual alpha
 * one, then the same for beta; each block is a run of lines that all repeat the
 * `Alpha  occ. eigenvalues --` banner. The occupied count is therefore known
 * the moment the virtual block starts, which is what separates HOMO from LUMO
 * without a second pass over the numbers.
 */
interface OrbitalDump {
    alpha: number[];
    /** Where the alpha virtual block begins; -1 while no virtual was printed. */
    alphaOcc: number;
    beta: number[];
    betaOcc: number;
    /** The geometry whose SCF printed this dump (pre-filter frame index). */
    frame: number;
}

interface OrbitalState {
    /** The last dump whose virtual block was reached — the one that counts. */
    last: OrbitalDump | null;
    cur: OrbitalDump | null;
}

/** ` Alpha  occ. eigenvalues --  -77.28310 ...` — one line of an eigenvalue dump. */
const ORBITAL_EIGEN_RE = /^\s*(Alpha|Beta)\s+(occ|virt)\.\s*eigenvalues\s*--(.*)$/;

/**
 * Record one eigenvalue line.
 *
 * An optimisation prints one dump per SCF, and the panel must describe the
 * geometry the run ended on, so the last complete dump wins. A dump is complete
 * once its alpha *virtual* block starts — that is the only boundary the file
 * states, and it is why the previous dump is retired there rather than at the
 * next occupied line: an occupied block runs over dozens of lines, and treating
 * each as a new dump would leave only the last five eigenvalues.
 */
function noteOrbitalLine(st: OrbitalState, spin: string, kind: string, text: string, frame: number): void {
    const tokens = text.trim().split(/\s+/);
    const values: number[] = [];
    for (let k = 0; k < tokens.length; k++) {
        const v = parseLogNumber(tokens[k]);
        if (isFinite(v)) values.push(v);
    }
    if (values.length === 0) return;
    // A dump starts at the first line printed after a geometry's SCF closed, so
    // the frame on the list when it begins is the one it describes.
    const fresh = (): OrbitalDump => ({ alpha: [], alphaOcc: -1, beta: [], betaOcc: -1, frame });
    if (!st.cur) st.cur = fresh();
    if (spin === 'Alpha' && kind === 'occ' && st.cur.alphaOcc >= 0) {
        st.last = st.cur;
        st.cur = fresh();
    }
    if (spin === 'Alpha') {
        if (kind === 'virt' && st.cur.alphaOcc < 0) st.cur.alphaOcc = st.cur.alpha.length;
        for (let k = 0; k < values.length; k++) st.cur.alpha.push(values[k]);
    } else {
        if (kind === 'virt' && st.cur.betaOcc < 0) st.cur.betaOcc = st.cur.beta.length;
        for (let k = 0; k < values.length; k++) st.cur.beta.push(values[k]);
    }
}

/**
 * The last dump as the level diagram needs it.
 *
 * Occupations are the file's own occupied/virtual split: two electrons per
 * orbital when the file printed no beta block (a closed-shell job), one per
 * spin orbital otherwise. A log prints the split, not a population, so nothing
 * finer is claimed.
 */
function orbitalLevelsOf(st: OrbitalState): OrbitalLevels | undefined {
    // A dump whose virtual block was reached is the converged one; a file that
    // was cut inside its occupied block still has those numbers, so it is used
    // only when nothing better exists.
    const complete = st.cur && st.cur.alphaOcc >= 0 ? st.cur : null;
    const dump = complete || st.last || st.cur;
    if (!dump || dump.alpha.length === 0) return undefined;
    const closed = dump.beta.length === 0;
    const alphaOcc = dump.alphaOcc >= 0 ? dump.alphaOcc : dump.alpha.length;
    const betaOcc = dump.betaOcc >= 0 ? dump.betaOcc : dump.beta.length;
    const spins: OrbitalSpinLevels[] = [];
    spins.push({
        label: 'Alpha',
        energies: dump.alpha,
        occupations: occupationsOfLevels(dump.alpha.length, alphaOcc, closed ? 2 : 1),
        occupied: alphaOcc
    });
    if (!closed) {
        spins.push({
            label: 'Beta',
            energies: dump.beta,
            occupations: occupationsOfLevels(dump.beta.length, betaOcc, 1),
            occupied: betaOcc
        });
    }
    const levels: OrbitalLevels = { spins, closedShell: closed };
    if (dump.frame >= 0) levels.frame = dump.frame;
    return levels;
}

function occupationsOfLevels(count: number, occupied: number, perOrbital: number): number[] {
    const out = new Array<number>(count).fill(0);
    for (let i = 0; i < Math.min(count, occupied); i++) out[i] = perOrbital;
    return out;
}

/** The cavity block's own lines, read in one window (see the scan below). */
function noteCavityLines(cavity: SolventCavity, lines: string[], from: number, to: number): void {
    for (let k = from; k < to; k++) {
        const line = lines[k];
        // The continuum banner names the family the solver belongs to:
        // `Polarizable Continuum Model (PCM)`, `(CPCM)`, `(IEFPCM)`, …
        if (cavity.model === undefined && line.indexOf('Polarizable Continuum Model') >= 0) {
            const m = /Polarizable Continuum Model\s*\(([^)]+)\)/.exec(line);
            if (m) cavity.model = m[1].trim();
        }
        if (cavity.type === undefined && line.indexOf('Cavity type') >= 0) {
            const m = line.match(/Cavity type\s*:\s*(.+?)\.?\s*$/);
            if (m) {
                cavity.type = m[1].trim();
                const alpha = /Alpha=\s*([0-9.]+)/.exec(m[1]);
                if (alpha) cavity.alpha = parseLogNumber(alpha[1]);
            }
        }
        if (cavity.solvent === undefined && line.indexOf('Solvent') >= 0 && line.indexOf('Eps') >= 0) {
            const m = line.match(/Solvent\s*:\s*([^,]+),/);
            if (m) cavity.solvent = m[1].trim();
        }
        if (cavity.sphereCount === undefined && line.indexOf('NSphG=') >= 0) {
            const m = line.match(/NSphG=\s*(\d+)/);
            if (m) cavity.sphereCount = parseInt(m[1], 10);
        }
        if (cavity.area === undefined && line.indexOf('Cavity surface area') >= 0) {
            const m = line.match(/=\s*([0-9.]+)/);
            if (m) cavity.area = parseLogNumber(m[1]);
        }
        if (cavity.volume === undefined && line.indexOf('Cavity volume') >= 0) {
            const m = line.match(/=\s*([0-9.]+)/);
            if (m) cavity.volume = parseLogNumber(m[1]);
        }
    }
}

/**
 * Forget the cavity of the previous geometry.
 *
 * An optimisation builds the cavity again at every step and prints one PCM
 * block per geometry, so the panel — which describes the structure on screen —
 * must keep the LAST block. Each banner starts a new geometry's cavity, and
 * clearing here is what stops the first step's area and volume from winning a
 * `=== undefined` race for the rest of the file.
 */
function clearCavity(cavity: SolventCavity): void {
    cavity.model = undefined;
    cavity.type = undefined;
    cavity.solvent = undefined;
    cavity.alpha = undefined;
    cavity.area = undefined;
    cavity.volume = undefined;
    cavity.sphereCount = undefined;
}

/** Route tokens that are parameters of `SCRF`, not the model itself. */
const SCRF_PARAMETER_TOKENS = [
    'SOLVENT', 'READ', 'EXTERNALITERATION', 'NONE', 'EPS', 'EPSINF', 'RSOLV',
    'CAVITY', 'SURFACETYPE', 'SURFACERADIUS', 'ADDITIVE', 'TENSER'
];

/**
 * The solvation model named by the route, or undefined when it names none.
 *
 * The continuum banner says `PCM` for every IEFPCM-family job — SMD included —
 * so the route is the only place that distinguishes `SMD` from `PCM`; the LAST
 * route of a `--Link1--` file is the one that describes the cavity kept below.
 */
export function scrfModelOf(content: string): string | undefined {
    const re = /SCRF\s*=\s*\(?\s*([A-Za-z][A-Za-z0-9-]*)/gi;
    let match: RegExpExecArray | null;
    let last: string | undefined;
    while ((match = re.exec(content)) !== null) {
        const token = match[1];
        if (SCRF_PARAMETER_TOKENS.indexOf(token.toUpperCase()) >= 0) continue;
        last = token;
    }
    return last;
}

/** How far into a log the title card is looked for. Gaussian's own banner, the
 *  route card and the IOp block put it well inside the first 120 lines; the
 *  window stops a 60 MB log from paying a second full pass for a one-line field. */
const TITLE_SCAN_LINES = 400;

/**
 * The job title a Gaussian log echoes: the line between two dashed rules that
 * immediately precede `Symbolic Z-matrix:`.
 *
 * The title card sits AFTER the route section, not before it — the line above
 * a route card is the rule that opens the card, so the previous anchor (the
 * last non-empty, non-`%` line before the first `#`) captured `-----…` and put
 * 51 hyphens into every GJF/XYZ title and every Diff label the file produced.
 * Anchoring on `Symbolic Z-matrix:` also survives a `--Link1--` job whose
 * title card is preceded by a different banner.
 */
function gaussianTitleOf(lines: string[]): string {
    const limit = Math.min(lines.length, TITLE_SCAN_LINES);
    let symbolic = -1;
    for (let i = 0; i < limit; i++) {
        if (lines[i].indexOf('Symbolic Z-matrix') >= 0) { symbolic = i; break; }
    }
    if (symbolic >= 0) {
        // Walk up over the blank line and the closing rule, then read the title
        // back to the rule that opens it. A title card may be empty, which is a
        // statement (no title), not a reason to grab the line above it.
        let k = symbolic - 1;
        while (k >= 0 && lines[k].trim() === '') k--;
        if (k < 0 || !DASH_RE.test(lines[k].trim())) return '';
        k--;
        const end = k;
        // A title card is one line (Gaussian joins a multi-line title with
        // blanks); the walk is bounded so a box with no opening rule cannot
        // swallow the banner above it.
        for (let guard = 0; guard < 4 && k >= 0; guard++) {
            const t = lines[k].trim();
            if (t === '' || DASH_RE.test(t)) break;
            k--;
        }
        const parts: string[] = [];
        for (let j = k + 1; j <= end; j++) parts.push(lines[j].trim());
        return parts.join(' ').trim();
    }

    // No Z-matrix section (a checkpoint-restart continuation, a shape written
    // by another tool): keep the older anchor, but a rule is a rule, never a
    // title.
    let routeIdx = -1;
    for (let i = 0; i < limit; i++) {
        if (ROUTE_START_RE.test(lines[i].trim())) { routeIdx = i; break; }
    }
    if (routeIdx < 0) return '';
    for (let i = routeIdx - 1; i >= 0; i--) {
        const t = lines[i].trim();
        if (t === '' || t.startsWith('%') || DASH_RE.test(t)) continue;
        return t;
    }
    return '';
}

/**
 * One scan point while the single pass reads it.
 *
 * A relaxed scan is a series of *optimisations*, not a series of steps: each
 * point relaxes everything but the driven coordinate, prints its own tables and
 * ends on its own stationary point. The draft keeps the point's own numbers —
 * its last table, the geometry that table describes and how many steps it took
 * — so the finished series is one entry per point and never a 106-step
 * "convergence" made of eight separate relaxations.
 */
interface ScanDraft {
    /** 1-based scan point, as the file numbers it. */
    point: number;
    /** The step number the point's latest table carried. */
    step: number;
    /** How many tables the point printed. */
    steps: number;
    /** How many of the point's four convergence rows said `YES`. */
    yes: number;
    /** True when the file stated `Optimization completed.` for this point. */
    converged: boolean;
    /** PRE-FILTER index of the geometry of the point's last step, or -1. */
    frame: number;
    /** The point's last convergence table. */
    last?: OptStep;
}

/** The internal coordinate a scan drives, as the file's parameter table names it. */
interface ScanCoordInfo {
    /** `R2`, `D13` — the name the summary table repeats. */
    name: string;
    /** `R(1,6)`, `D(4,3,5,6)` — the definition a reader recognises. */
    label: string;
    /** Angstrom for a bond length, degrees for an angle or a dihedral. */
    unit: 'A' | 'deg';
}

/** The run's closing summary: the energy and the driven coordinate per point. */
interface ScanTable {
    /** Point number → energy in Hartree (the table's own base added back). */
    energy: { [point: number]: number };
    /** Parameter name → point number → value, in the file's own units. */
    values: { [name: string]: { [point: number]: number } };
}

interface ScanState {
    /** The point being read. */
    cur: ScanDraft | null;
    /** Finished points, in the order the file took them. */
    list: ScanDraft[];
    /** How many points the run set out to take, when its step line says so. */
    total: number;
    /** The coordinate(s) whose derivative column reads `Scan`. */
    params: ScanCoordInfo[];
    /** The closing summary table, when the run printed one. */
    table: ScanTable | null;
}

/** One IRC path point, recorded when its geometry frame is built. */
interface IrcDraft {
    path: number;
    point: number;
    /** PRE-FILTER index of the frame holding that point's geometry. */
    frame: number;
    /** Which IRC calculation of the file this point belongs to (1-based):
     *  a `--Link1--` file can run the analysis more than once, and each run is
     *  one reaction path of its own. */
    run: number;
    /** Which half of the path, when the file says (a banner or its route card). */
    dir?: IrcDirection;
}

/**
 * Role a geometry block plays in the frame list, kept beside `frames` because
 * the filter that drops the redundant blocks runs AFTER the pass and the labels
 * are no longer a reliable statement of it (a scan point's geometry is labelled
 * by its point, not by `Step`).
 */
interface FrameKinds {
    /** True for an `Input orientation` block. */
    input: boolean[];
    /** True for a frame a scan point or an IRC path point claims. */
    bound: boolean[];
}

/**
 * Read the `! Initial Parameters !` block a Berny scan opens with and return
 * the rows whose derivative column says `Scan` — the coordinate(s) the run
 * drives.
 *
 *     ! R2    R(1,6)                  3.0978         Scan                            !
 *
 * The block is bounded by its own dashed rules, so the read stops at the first
 * line that is neither a row nor a rule: a damaged file cannot make it swallow
 * the geometry below.
 */
function readScanParameters(lines: string[], from: number): { params: ScanCoordInfo[]; endIdx: number } {
    const params: ScanCoordInfo[] = [];
    const limit = Math.min(lines.length, from + 4000);
    let rows = 0;
    let i = from + 1;
    for (; i < limit; i++) {
        const t = lines[i].trim();
        if (t === '') break;
        if (t.charAt(0) === '!') {
            const m = SCAN_PARAM_ROW_RE.exec(t);
            if (m) {
                rows++;
                if (SCAN_MARK_RE.test(m[4])) {
                    const name = m[1];
                    params.push({
                        name,
                        label: m[2],
                        unit: name.charAt(0).toUpperCase() === 'R' ? 'A' : 'deg'
                    });
                }
            }
            continue;
        }
        if (t.charAt(0) !== '-') break;
        // The rule above the first row opens the table; a rule after it closes.
        if (rows > 0) { i++; break; }
    }
    return { params, endIdx: i };
}

/**
 * Read the `Summary of Optimized Potential Surface Scan` table.
 *
 * The table is printed in groups of five points, each group repeating the
 * parameter rows under its own column header:
 *
 *     Summary of Optimized Potential Surface Scan (add -269.0 to energies):
 *                            1         2         3         4         5
 *          Eigenvalues --    -0.65139  -0.65102  -0.65047  -0.64973  -0.64879
 *                R2           3.09780   3.04780   2.99780   2.94780   2.89780
 *
 * Values are packed without a separating blank when two are negative
 * (`-168.89124-167.99169`), so they are read by their decimal shape rather than
 * by whitespace. The header's `add -269.0 to energies` base is put back, so an
 * energy read here is a total energy like every other one in the viewer.
 */
function readScanTable(lines: string[], from: number): ScanTable {
    const table: ScanTable = { energy: {}, values: {} };
    const baseMatch = lines[from].match(SCAN_TABLE_BASE_RE);
    const base = baseMatch ? parseFloat(baseMatch[1]) : 0;
    let cols: number[] = [];
    let rows = 0;
    const limit = Math.min(lines.length, from + 40000);
    for (let i = from + 1; i < limit; i++) {
        const t = lines[i].trim();
        if (t === '') break;
        // A group's column header: the point numbers the rows below belong to.
        if (SCAN_TABLE_HEAD_RE.test(t)) {
            cols = t.split(/\s+/).map(Number);
            continue;
        }
        if (cols.length === 0) break;
        const eg = SCAN_TABLE_ENERGY_RE.exec(t);
        if (eg) {
            putScanValues(table.energy, cols, eg[1], base);
            rows++;
            continue;
        }
        const row = SCAN_TABLE_ROW_RE.exec(t);
        if (!row) break;
        const target = table.values[row[1]] || (table.values[row[1]] = {});
        putScanValues(target, cols, row[2], 0);
        rows++;
    }
    return table;
}

/** Lay the numbers of one summary row onto the points its column header named. */
function putScanValues(target: { [point: number]: number }, cols: number[], text: string, offset: number): void {
    const nums = text.match(SCAN_TABLE_NUMBER_RE);
    if (!nums) return;
    const n = Math.min(nums.length, cols.length);
    for (let k = 0; k < n; k++) target[cols[k]] = parseFloat(nums[k]) + offset;
}

/**
 * Where a PRE-FILTER frame index ended up after the input/standard filter.
 * Undefined when the frame was dropped, which a caller treats as "no structure
 * to point at" rather than as frame 0.
 */
function remapFrameIndex(remap: number[] | null, idx: number, preCount: number): number | undefined {
    if (idx < 0 || idx >= preCount) return undefined;
    if (!remap) return idx;
    let k = idx;
    while (k >= 0 && remap[k] < 0) k--;
    return k < 0 ? undefined : remap[k];
}

/** The label of a scan point's geometry: `Point 3 · step 2 (Charge=1 Mult=1)`. */
function scanFrameLabel(point: number, step: number, chargeMult: string): string {
    return 'Point ' + point + ' · step ' + step + (chargeMult ? ' (' + chargeMult + ')' : '');
}

/**
 * Name every frame an IRC path claims by the row the card lists it under.
 *
 * The parser can only label a geometry when it meets it, and the transition
 * state is printed above the path banner with no point header at all — it used
 * to keep the generic `Input` label while the card called the same structure
 * `TS`. Once the paths and the run that owns them are known, the labels are
 * rewritten from those links: `TS`, `R12`, `F3`, with the run number when a
 * `--Link1--` file holds more than one IRC calculation. The electronic state
 * the file printed stays in the suffix.
 */
function labelIrcFrames(frames: LogFrame[], irc: IrcResult | undefined): void {
    if (!irc) return;
    const dirOf = new Map<number, IrcDirection>();
    for (const p of irc.paths) if (p.dir) dirOf.set(p.path, p.dir);
    const runOfPath = new Map<number, number>();
    for (let r = 0; r < irc.runs.length; r++) {
        for (const pn of irc.runs[r].paths) runOfPath.set(pn, r + 1);
    }
    const manyRuns = irc.runs.length > 1;
    const set = (frame: number | undefined, name: string, run: number | undefined): void => {
        if (typeof frame !== 'number' || frame < 0 || frame >= frames.length) return;
        const f = frames[frame];
        const cm = (f.charge != null && f.multiplicity != null)
            ? ' (Charge=' + f.charge + ' Mult=' + f.multiplicity + ')' : '';
        const label = (manyRuns && run ? 'Run ' + run + ' · ' : '') + name + cm;
        f.stepLabel = label;
        f.title = label;
    };
    for (const p of irc.points) {
        const dir = dirOf.get(p.path);
        const name = (dir === 'reverse' ? 'R' : dir === 'forward' ? 'F' : '') + p.point;
        set(p.frame, name, runOfPath.get(p.path));
    }
    for (let r = 0; r < irc.runs.length; r++) {
        set(irc.runs[r].tsFrame, 'TS', manyRuns ? r + 1 : undefined);
    }
}

/** Two frames describe the same geometry, element for element. */
function sameGeometry(a: Atom[], b: Atom[]): boolean {
    if (a.length !== b.length) return false;
    for (let k = 0; k < a.length; k++) {
        if (a[k].element !== b[k].element) return false;
        if (Math.abs(a[k].x - b[k].x) > 1e-6) return false;
        if (Math.abs(a[k].y - b[k].y) > 1e-6) return false;
        if (Math.abs(a[k].z - b[k].z) > 1e-6) return false;
    }
    return true;
}

/**
 * Turn the scan drafts into the series the page draws, in frame-list indices.
 *
 * The two sources are combined rather than chosen: the tables give every point's
 * energy, forces and geometry, and the closing summary gives the driven
 * coordinate and a value for a point whose table the run never reached (a job
 * that stopped early). A file whose tables are missing entirely — a scan read
 * only through its summary, or a log whose head is gone — still yields its
 * points, with the energies the summary states.
 */
function buildScanResult(st: ScanState, frames: LogFrame[], remap: number[] | null): ScanResult | undefined {
    const preCount = remap ? remap.length : frames.length;
    const info = st.params.length > 0 ? st.params[0] : null;
    const points: ScanPoint[] = [];
    for (const d of st.list) {
        const p: ScanPoint = { point: d.point };
        if (d.last) {
            p.energy = d.last.energy;
            p.maxForce = d.last.maxForce;
            p.rmsForce = d.last.rmsForce;
            p.maxDisplacement = d.last.maxDisplacement;
            p.rmsDisplacement = d.last.rmsDisplacement;
        }
        if (p.energy == null && st.table && st.table.energy[d.point] != null) {
            p.energy = st.table.energy[d.point];
        }
        if (d.steps > 0) p.steps = d.steps;
        const fr = remapFrameIndex(remap, d.frame, preCount);
        if (fr !== undefined) p.frame = fr;
        // The point is converged when the file said so, and equally when all
        // four of its last table's rows said YES: a run cut off before its
        // "Optimization completed." line still states where it got to.
        p.converged = d.converged || (d.last !== undefined && d.yes === 4);
        if (info && st.table) {
            const col = st.table.values[info.name];
            if (col && col[d.point] != null) p.coord = col[d.point];
        }
        points.push(p);
    }
    // No table of its own: the run's summary still states its points.
    if (points.length === 0 && st.table) {
        const nums = Object.keys(st.table.energy).map(Number).sort((a, b) => a - b);
        for (const n of nums) {
            const p: ScanPoint = { point: n, energy: st.table.energy[n] };
            if (info) {
                const col = st.table.values[info.name];
                if (col && col[n] != null) p.coord = col[n];
            }
            points.push(p);
        }
    }
    if (points.length === 0) return undefined;
    points.sort((a, b) => a.point - b.point);
    const result: ScanResult = { points };
    // A log whose head is gone has no step line to state the run's size; the
    // summary's highest point number is then the size it actually reached.
    let total = st.total;
    if (total <= 0 && st.table) {
        const nums = Object.keys(st.table.energy).map(Number);
        if (nums.length > 0) total = Math.max.apply(null, nums);
    }
    if (total > 0) result.total = total;
    if (info) result.coord = { label: info.label, unit: info.unit };
    return result;
}

/** The path directions an `irc=(…)` route card asks for, in its own order. */
function ircRouteDirections(sec: RouteSection): IrcDirection[] {
    const dirs: IrcDirection[] = [];
    for (const k of sec.keywords) {
        if (k.name.toLowerCase() !== 'irc') continue;
        for (const o of k.options) {
            const v = o.toLowerCase();
            if (v === 'forward' && dirs.indexOf('forward') < 0) dirs.push('forward');
            else if (v === 'reverse' && dirs.indexOf('reverse') < 0) dirs.push('reverse');
        }
    }
    return dirs;
}

/** Turn the IRC path points into the series the page draws.
 *
 * A path number is only unique inside its own IRC job: a `--Link1--` file that
 * ran the analysis twice prints `Path Number: 1` again for the second run.
 * Each run of consecutive drafts that keeps the same file-stated path number is
 * therefore given its own series number — 1, 2, 3, … in file order. A single
 * IRC job keeps its own numbers exactly (1, 2), and two jobs can never merge
 * their points into one chart.
 *
 * The drafts are grouped by the run they belong to as well: the page joins one
 * run's forward and reverse halves into ONE reaction path, and a file with two
 * IRC jobs still gets two charts. The transition state the run started from is
 * carried with it, so the joined curve has a single barrier between the two
 * halves instead of two near-degenerate points where they meet. */
function buildIrcResult(drafts: IrcDraft[], frames: LogFrame[], remap: number[] | null, runTs: { run: number; frame: number }[]): IrcResult | undefined {
    if (drafts.length === 0) return undefined;
    const preCount = remap ? remap.length : frames.length;
    const points: IrcPoint[] = [];
    const paths: { path: number; count: number; dir?: IrcDirection }[] = [];
    const pointSeries: number[] = [];
    let prevPath = -1;
    let series = 0;
    for (const d of drafts) {
        if (d.path !== prevPath) {
            prevPath = d.path;
            series = paths.length + 1;
            paths.push({ path: series, count: 1 });
        } else {
            paths[paths.length - 1].count++;
        }
        pointSeries.push(series);
        if (d.dir && !paths[paths.length - 1].dir) paths[paths.length - 1].dir = d.dir;
        const p: IrcPoint = { path: series, point: d.point };
        const fr = remapFrameIndex(remap, d.frame, preCount);
        if (fr !== undefined) {
            p.frame = fr;
            const e = frames[fr].energy;
            if (e != null) p.energy = e;
        }
        points.push(p);
    }

    // One run per consecutive group of drafts. A run that computed both halves
    // names each half; a file that stated neither direction still printed the
    // halves in forward/reverse order, which is how they can be joined.
    const runs: IrcRun[] = [];
    let start = 0;
    while (start < drafts.length) {
        const runNo = drafts[start].run;
        let end = start;
        while (end < drafts.length && drafts[end].run === runNo) end++;
        const runPathNumbers: number[] = [];
        for (let k = start; k < end; k++) {
            if (runPathNumbers.indexOf(pointSeries[k]) < 0) runPathNumbers.push(pointSeries[k]);
        }
        const runPaths: { path: number; count: number; dir?: IrcDirection }[] = [];
        for (const pn of runPathNumbers) {
            const meta = paths.filter(p => p.path === pn)[0];
            if (meta) runPaths.push(meta);
        }
        if (runPaths.length === 2) {
            const a = runPaths[0];
            const b = runPaths[1];
            if (!a.dir && !b.dir) { a.dir = 'forward'; b.dir = 'reverse'; }
            else if (!a.dir && b.dir === 'reverse') a.dir = 'forward';
            else if (!b.dir && a.dir === 'forward') b.dir = 'reverse';
            else if (!a.dir && b.dir === 'forward') a.dir = 'reverse';
            else if (!b.dir && a.dir === 'reverse') b.dir = 'forward';
        }
        const run: IrcRun = { paths: runPathNumbers };
        const ts = runTs.filter(t => t.run === runNo)[0];
        if (ts && (remap === null || remap[ts.frame] >= 0)) {
            const fr = remapFrameIndex(remap, ts.frame, preCount);
            const claimed = fr !== undefined && drafts.some(d => remapFrameIndex(remap, d.frame, preCount) === fr);
            if (fr !== undefined && !claimed) {
                run.tsFrame = fr;
                const e = frames[fr].energy;
                if (e != null) run.tsEnergy = e;
            }
        }
        runs.push(run);
        start = end;
    }
    return { points, paths, runs };
}

/**
 * Order an IRC calculation's frames along its reaction path.
 *
 * The file prints the transition state first, then one half of the path and
 * then the other, so its frame list begins in the middle of the reaction
 * coordinate. The viewer plays frames in list order: it would start at the
 * barrier and jump from one far end to the other when it crossed between the
 * halves. A reaction path has a natural order — the reverse half from its far
 * end to the barrier, the barrier itself, then the forward half away from it —
 * and that is the order the IRC card lists. Putting the frames in that order
 * makes frame 1 the leftmost structure, so opening the file shows the start of
 * the path and Play runs it from one end to the other.
 *
 * A run that computed one direction keeps the barrier on the side the half
 * starts from: a forward-only run begins at the barrier, a reverse-only one
 * ends there. A file with several IRC jobs keeps each run's frames together
 * where the run began. Every frame reference — the path points, each run's
 * transition state, scan points and charge tables — is carried across the
 * move, so a click on a path row, a charge map and the chart highlight still
 * name the geometry they describe.
 */
function orderIrcFrames(
    frames: LogFrame[],
    irc: IrcResult | undefined,
    scan: ScanResult | undefined,
    charges: LogChargeSet[],
    levels?: OrbitalLevels
): void {
    if (!irc || irc.runs.length === 0 || frames.length < 2) return;

    const members = new Set<number>();
    const sequences: number[][] = [];
    for (const run of irc.runs) {
        let fwd: number | null = null;
        let rev: number | null = null;
        for (const pn of run.paths) {
            const meta = irc.paths.filter(p => p.path === pn)[0];
            if (!meta || !meta.dir) continue;
            if (meta.dir === 'forward' && fwd === null) fwd = pn;
            else if (meta.dir === 'reverse' && rev === null) rev = pn;
        }
        const seq: number[] = [];
        const seen = new Set<number>();
        const push = (fr: number | undefined): void => {
            if (fr === undefined || fr < 0 || fr >= frames.length || seen.has(fr)) return;
            seen.add(fr);
            seq.push(fr);
        };
        const pointsOf = (pn: number): IrcPoint[] =>
            irc.points.filter(p => p.path === pn);
        if (rev !== null) {
            // Point 0 is the end nearest the barrier; the far end is where the
            // joined path starts.
            const pts = pointsOf(rev).sort((a, b) => b.point - a.point);
            for (const p of pts) push(p.frame);
        }
        push(run.tsFrame);
        if (fwd !== null) {
            const pts = pointsOf(fwd).sort((a, b) => a.point - b.point);
            for (const p of pts) push(p.frame);
        }
        // A series the file neither named nor paired (an unusual payload with
        // more than two paths) stays with its run, after the named halves,
        // rather than being left behind.
        const rest = irc.points
            .filter(p => run.paths.indexOf(p.path) >= 0)
            .sort((a, b) => a.path - b.path || a.point - b.point);
        for (const p of rest) push(p.frame);
        if (seq.length === 0) continue;
        for (const fr of seq) members.add(fr);
        sequences.push(seq);
    }
    if (sequences.length === 0 || members.size < 2) return;

    // Each run is emitted where its first structure was; every frame outside
    // the runs keeps its place. A run's frames are contiguous in the file, so
    // this reorders the path without moving another job's frames around it.
    const byAnchor = new Map<number, number[]>();
    for (const seq of sequences) {
        const anchor = Math.min(...seq);
        if (byAnchor.has(anchor)) return;   // overlapping runs: leave the file alone
        byAnchor.set(anchor, seq);
    }
    const next: LogFrame[] = [];
    const newIndex: number[] = new Array(frames.length).fill(-1);
    const placed = new Set<number>();
    for (let k = 0; k < frames.length; k++) {
        if (placed.has(k)) continue;
        const seq = byAnchor.get(k);
        if (seq) {
            for (const fr of seq) {
                if (placed.has(fr)) return;
                placed.add(fr);
                newIndex[fr] = next.length;
                next.push(frames[fr]);
            }
        } else if (!members.has(k)) {
            placed.add(k);
            newIndex[k] = next.length;
            next.push(frames[k]);
        }
    }
    // Every frame must be accounted for; otherwise the list is left untouched
    // rather than silently dropping or duplicating a structure.
    if (next.length !== frames.length || newIndex.some(v => v < 0)) return;

    for (const p of irc.points) {
        if (typeof p.frame === 'number') p.frame = newIndex[p.frame];
    }
    for (const run of irc.runs) {
        if (typeof run.tsFrame === 'number') run.tsFrame = newIndex[run.tsFrame];
    }
    if (scan) {
        for (const p of scan.points) {
            if (typeof p.frame === 'number') p.frame = newIndex[p.frame];
        }
    }
    for (const c of charges) c.frame = newIndex[c.frame];
    if (levels && typeof levels.frame === 'number') levels.frame = newIndex[levels.frame];
    frames.length = 0;
    for (const f of next) frames.push(f);
}

export function parseGaussianLog(content: string, onProgress?: (fraction: number) => void): GaussianLogResult {
    const lines = content.split(/\r?\n/);
    const lineTotal = lines.length || 1;
    let lastProgress = 0;
    const frames: LogFrame[] = [];
    const title = gaussianTitleOf(lines);

    // The charge and multiplicity printed with the geometry (`Charge = 0
    // Multiplicity = 1`). A `--Link1--` file can run several jobs with
    // different values, so the scan updates them as it passes each job's line
    // and every frame keeps its own job's values; the old single match at the
    // top of the file stamped the first job onto all of them.
    let chargeMultLine = '';
    let logCharge: number | undefined;
    let logMultiplicity: number | undefined;

    // Single pass over the file: geometry frames, convergence-table steps and
    // harmonic-frequency sections are all detected in ONE traversal (these
    // used to be three separate full scans). Consumed branches skip only
    // atom-coordinate rows / frequency tables, which never match the other
    // detectors, so the merged scan yields exactly the previous results.
    let stdFrames = 0;          // standard-orientation frame counter
    const optSteps: OptStep[] = [];
    let pendingEnergy: number | undefined;
    let stepCounter = 0;        // convergence-table step counter
    let modes: NormalMode[] | undefined;
    let fallbackSeen = false;
    let fallbackFrame: LogFrame | undefined;
    // Optimization steps are the optimizer's per-step Item/Value/Converged?
    // tables — but SP jobs that compute gradients (Force keyword) and freq
    // jobs print the SAME gradient banner and the same 4-metric table without
    // running an optimizer. The reliable trigger is the route card of the job
    // that printed the table: only a job that asks for an optimization (opt /
    // optimization / freq=calcall) runs Berny's step loop, so a table is taken
    // only while such a job is the one being read. The distinction matters for
    // the `opt freq` composite every Gaussian run writes as two jobs: the
    // `Freq` job repeats the table once at the geometry it was handed, under
    // '! Non-Optimized Parameters !', and appending it would add a pseudo-step
    // whose force is ~1e-6 — a drop the optimization never made. A job that
    // stated no card at all (a log whose head is gone) cannot be classified, so
    // its tables are taken. All route sections are collected (parsed into
    // RouteSection[]) for the Route panel.
    let jobHasOpt = false;      // the current job's card requests an optimization
    let jobCardSeen = false;    // the current job stated a card at all
    const routes: RouteSection[] = [];
    // The per-atom charge tables of the job currently being read, and the
    // headline results of every job in the file. Both are plain boxes handed to
    // the helpers below rather than variables closed over by them: this loop
    // visits every line of the file, and a function that captures its state
    // makes the engine reach for a context object on each access inside it.
    const charges: ChargeState = { drafts: [], job: 0 };
    const jobs: SummaryState = { list: [], cur: null };
    /** The eigenvalue dumps of the job being read (the last one wins). */
    const orbitals: OrbitalState = { last: null, cur: null };
    /** The cavity the job's solvation model states, when it states one. */
    const cavity: SolventCavity = { source: 'log' };
    /** How many SCFs of the job being read have closed (1-based for the next
     *  one). A geometry-level value is stamped with it, so a summary can tell
     *  a value of the final geometry from one the file only printed earlier. */
    let scfCount = 0;
    /** Where the current job's frames begin, so an energy can never be bound to
     *  a frame the previous job of a `--Link1--` file printed. */
    let jobFrameStart = 0;
    /** The relaxed scan being read: its points so far, the coordinate its
     *  parameter table marks `Scan`, and the summary table it may close with. */
    const scanState: ScanState = { cur: null, list: [], total: 0, params: [], table: null };
    /** The IRC point whose geometry comes next, once its own header was read. */
    let ircPending: { path: number; point: number; dir?: IrcDirection } | null = null;
    /** The path number the last `Point Number` header stated. */
    let ircPathNumber = 0;
    /** The path whose starting geometry comes next, when the file announces a
     *  half of the path by a banner instead of a point header: the reverse path
     *  begins right after `Beginning calculation of the REVERSE path.` with a
     *  geometry that carries no `Point Number` of its own. It is point 0. */
    let ircStartPath: number | null = null;
    /** The IRC calculation being read: a file may hold several. */
    let ircRun = 0;
    /** True once the current run has a start banner or a bound point, so a
     *  geometry that follows a job boundary without a banner opens a new one. */
    let ircRunOpen = false;
    /** Which half the next path point belongs to, when the file says so. */
    let ircDirection: IrcDirection | undefined = undefined;
    /** The directions the current job's own route card states, if any. */
    let jobIrcDirs: IrcDirection[] | null = null;
    /** The transition-state frame each run starts from (PRE-FILTER indices). */
    const ircRunTs: { run: number; frame: number }[] = [];
    const ircDrafts: IrcDraft[] = [];
    /** The role each frame block played, kept beside `frames` for the filter
     *  that runs after the pass (labels no longer state it). */
    const kinds: FrameKinds = { input: [], bound: [] };

    let i = 0;
    while (i < lines.length) {
        // A whole log is one synchronous pass, so the only honest progress is
        // how far into its lines the scan is. Reporting every 2048 lines keeps
        // the callback (and its IPC) off the per-line path.
        if (onProgress && i - lastProgress >= 2048) {
            onProgress(i / lineTotal);
            lastProgress = i;
        }
        const line = lines[i];

        // --- Charge/multiplicity of the job this frame belongs to ---
        // Gaussian prints the line once per job, before that job's geometries,
        // so updating the state here gives each frame the values of its own job.
        if (line.indexOf('Charge') >= 0 && line.indexOf('Multiplicity') >= 0) {
            const cm = line.match(CHARGE_MULT_RE);
            if (cm) {
                logCharge = parseInt(cm[1], 10);
                logMultiplicity = parseInt(cm[2], 10);
                chargeMultLine = `Charge=${cm[1]} Mult=${cm[2]}`;
                const s = summaryOf(jobs);
                s.charge = logCharge;
                s.multiplicity = logMultiplicity;
                i++;
                continue;
            }
        }

        // --- Route card (keywords line) at the start of each job section ---
        // Route cards are enclosed between two dashed separator lines (the
        // section may wrap onto continuation lines, e.g. '... Fre' / 'q' or
        // '... opt' / '=(calcall) freq'), so the dash+'#...' pairing anchors
        // the section and the closing dash bounds it; everything between the
        // dashes is route content. A '#...' line without its top dash (only
        // possible in non-standard output) still triggers collection bounded
        // by blank/dash lines, preserving the previous behavior.
        // Pre-filter: only '-'/'#'-leading lines (after blanks/tabs) can
        // participate, so long logs skip the regex work for most lines.
        {
            let p0 = 0;
            const L0 = line.length;
            while (p0 < L0 && (line.charCodeAt(p0) === 32 || line.charCodeAt(p0) === 9)) p0++;
            const c0 = p0 < L0 ? line[p0] : '';
            if (c0 === '-' || c0 === '#') {
                const tr = line.trim();
                let routeStart = -1;
                if (DASH_RE.test(tr)) {
                    if (i + 1 < lines.length && ROUTE_START_RE.test(lines[i + 1].trim())) routeStart = i + 1;
                } else if (ROUTE_START_RE.test(tr) &&
                           (tr.indexOf('/') >= 0 ||
                            (i + 1 < lines.length && DASH_RE.test(lines[i + 1].trim())))) {
                    // A card whose OPENING rule is missing (hand-edited or
                    // non-standard output). It must still look like a card: it
                    // states a method/basis pair, or its closing rule sits right
                    // below it. The bare '#…' shape alone is not enough — an IRC
                    // banner's `# OF POINTS ALONG THE PATH = 1` has that shape
                    // and used to be collected as one route section per point.
                    routeStart = i;
                }
                if (routeStart >= 0) {
                    let j = routeStart + 1;
                    while (j < lines.length) {
                        const t2 = lines[j].trim();
                        if (t2 === '' || DASH_RE.test(t2)) break;
                        if (ROUTE_START_RE.test(t2)) {
                            // A '#' line after a blank/dash line is a new section;
                            // a '#' line directly after a collected route line is a
                            // continuation (Gaussian wraps long route cards and may
                            // start the continuation with '#p ...' at the next line).
                            if (j > routeStart + 1) {
                                const prev = lines[j - 1].trim();
                                if (prev === '' || DASH_RE.test(prev)) break;
                            }
                        }
                        j++;
                    }
                    const sec = parseRouteSection(joinRouteLines(lines.slice(routeStart, j).map(l => l.trim())));
                    routes.push(sec);
                    // Which kind of job the following tables belong to.
                    jobHasOpt = sec.hasOpt;
                    jobCardSeen = true;
                    // The directions this job's own card states (an IRC run that
                    // computed only one half leaves no direction banner, so its
                    // card is the only place the file says which half it is).
                    jobIrcDirs = ircRouteDirections(sec);
                    // A route card is the start of a job: the results and the
                    // population analysis that follow belong to it, not to the
                    // job whose card preceded them. This is the only boundary a
                    // `--Link1--` file states, and the first card of a file
                    // finds nothing to close.
                    jobs.cur = null;
                    // The tables already read belong to the job that printed
                    // them and stay: a `--Link1--` file that analysed two jobs
                    // gave the panel only the last one's schemes, although the
                    // earlier tables still name a frame that is on the list.
                    // The job counter is what keeps a later table of the same
                    // scheme from overwriting an earlier job's.
                    charges.job++;
                    orbitals.last = null;
                    orbitals.cur = null;
                    clearCavity(cavity);
                    scfCount = 0;
                    // A job that ended without a convergence table (a single
                    // point) leaves its SCF energy pending; carrying it into the
                    // next job would stamp it onto that job's first geometry.
                    pendingEnergy = undefined;
                    // A scan belongs to the job that ran it: the point still
                    // open ends here rather than absorbing the next job's steps.
                    if (scanState.cur) { scanState.list.push(scanState.cur); scanState.cur = null; }
                    // An IRC path belongs to its job just as strictly: a point
                    // header waiting for a geometry, an announced path start or a
                    // path number from the previous job must not reach the next.
                    ircPending = null;
                    ircStartPath = null;
                    ircPathNumber = 0;
                    ircDirection = undefined;
                    ircRunOpen = false;
                    jobFrameStart = frames.length;
                    i = j;
                    continue;
                }
            }
        }

        // --- Harmonic frequencies section (parsed inline; the LAST section wins) ---
        if (line.includes('Harmonic frequencies') && line.includes('cm**-1')) {
            const res = parseModesFromHeader(lines, i);
            // A mode's displacement list is one 3-vector per atom of the
            // geometry the section was printed for — the structure above it, not
            // whatever the file ends on. A matrix cut short (or a section that
            // belongs to another job) gives lists of the wrong length, and those
            // modes cannot be animated on any structure: they are dropped, and a
            // section that yields nothing usable does not replace the set an
            // earlier job printed.
            const atomCount = frames.length > 0 ? frames[frames.length - 1].atoms.length : 0;
            const usable = atomCount > 0
                ? res.modes.filter(m => m.displacements.length === 0 || m.displacements.length === atomCount)
                : res.modes;
            if (usable.length > 0) modes = usable;
            i = res.endIdx;
            continue;
        }

        // --- Geometry frames (opt/freq logs) ---
        if (line.includes('Standard orientation:') || line.includes('Input orientation:')) {
            const isStandard = line.includes('Standard orientation:');

            i = skipDashedLines(lines, i + 1, 2);

            const atoms: Atom[] = [];
            while (i < lines.length) {
                const coordLine = lines[i].trim();
                if (coordLine === '' || coordLine.includes('---')) break;

                const parts = coordLine.split(/\s+/);
                if (parts.length >= 6) {
                    const atomicNum = parseInt(parts[1], 10);
                    if (isNaN(atomicNum)) { i++; continue; }
                    const element = ATOMIC_NUMBER_MAP[atomicNum] || 'X';
                    const x = parseFloat(parts[3]);
                    const y = parseFloat(parts[4]);
                    const z = parseFloat(parts[5]);

                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ element, x, y, z, index: atoms.length });
                    }
                }
                i++;
            }

            if (atoms.length > 0) {
                const cmText = chargeMultLine ? ' (' + chargeMultLine + ')' : '';
                // An IRC point states its own number BEFORE its geometry, so the
                // block that follows is that point's structure — label, frame
                // link and all. A path's starting geometry has no such header
                // when a banner announced the path (`ircStartPath`): that block
                // is point 0 of the announced path. Everything else is a step or
                // the input pose.
                const irc = ircPending ||
                    (ircStartPath !== null ? { path: ircStartPath, point: 0, dir: ircDirection } : null);
                const label = irc
                    ? 'Path ' + irc.path + ' · point ' + irc.point + cmText
                    : isStandard ? `Step ${stdFrames + 1}${cmText}` : 'Input' + cmText;
                // Which optimizer step this geometry is: its own table comes
                // AFTER the geometry (`Standard orientation` … `SCF Done` … the
                // four rows), so the number it will take is the next one. A job
                // whose tables are not a series (a Freq job re-reading the
                // converged geometry, an IRC) stamps nothing, and the page then
                // has no step to mark on those frames.
                //
                // A scan's tables are not one series either: a point's geometry
                // is stamped with that point's own step number when its table
                // arrives (see the convergence-table branch), so from the first
                // point on no frame claims a position in `optSteps`.
                const stepOfFrame = (!irc && scanState.cur === null && (jobHasOpt || !jobCardSeen))
                    ? optSteps.length + 1 : undefined;
                frames.push({
                    atoms,
                    bonds: [],
                    title: label,
                    hasExplicitBonds: false,
                    stepLabel: label,
                    charge: logCharge,
                    multiplicity: logMultiplicity,
                    optStep: stepOfFrame,
                    // No energy yet: this program prints a cycle's geometry
                    // BEFORE the energy computed at it (`Standard orientation` …
                    // `SCF Done` … convergence table), so the value is bound back
                    // to this frame by `assignPendingEnergyToFrames` when the
                    // energy line arrives. Reading `pendingEnergy` here instead
                    // was the bug: it handed the *previous* cycle's energy to any
                    // frame printed before the next one, and the ordinary opt log
                    // — where the convergence table clears the value in between —
                    // left every frame without one. [2026-11-12 | Atreides-Jimmy]
                    energy: undefined
                });
                kinds.input.push(!isStandard);
                // A frame an IRC point claims is never the redundant copy of a
                // step: it is the point's only structure, so the filter below
                // must keep it even if it is an `Input orientation` block.
                kinds.bound.push(irc !== null);
                if (irc) {
                    if (!ircRunOpen) { ircRun++; ircRunOpen = true; }
                    ircDrafts.push({
                        path: irc.path, point: irc.point, frame: frames.length - 1,
                        run: ircRun, dir: irc.dir
                    });
                    ircPending = null;
                    ircStartPath = null;
                }
                if (isStandard) stdFrames++;
            }
            i++; // matches the previous outer-loop increment
            continue;
        }

        // --- Fallback frame for old-style logs: only the first occurrence is
        // considered, and it is used only when no orientation frames exist ---
        if (!fallbackSeen && line.includes('Coordinates (Angstroms)')) {
            fallbackSeen = true;
            i = skipDashedLines(lines, i + 1, 2);

            const atoms: Atom[] = [];
            while (i < lines.length) {
                const coordLine = lines[i].trim();
                if (coordLine === '' || coordLine.includes('---')) break;

                const parts = coordLine.split(/\s+/);
                if (parts.length >= 4) {
                    const element = parts[1].replace(/[0-9]/g, '');
                    const el = element.charAt(0).toUpperCase() + element.slice(1).toLowerCase();
                    const x = parseFloat(parts[2]);
                    const y = parseFloat(parts[3]);
                    const z = parseFloat(parts[4]);

                    if (el && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ element: el, x, y, z, index: atoms.length });
                    }
                }
                i++;
            }

            if (atoms.length > 0) {
                fallbackFrame = {
                    atoms,
                    bonds: [],
                    title: 'Coordinates',
                    hasExplicitBonds: false,
                    stepLabel: 'Coordinates',
                    charge: logCharge,
                    multiplicity: logMultiplicity
                };
            }
            continue;
        }

        // --- Convergence tables, charge tables and the job summary ---
        //
        // Every value read below is printed as a banner: indented, and starting
        // with a capital letter. That shape is the gate — two character reads —
        // so the literal scans run on about a sixth of the lines instead of on
        // every coordinate row, frequency table and eigenvalue dump a log is
        // mostly made of. Nothing is lost by it: these programs indent their
        // banners, and a line they did not would be one this reader never
        // matched either.
        let gateP = 0;
        while (gateP < line.length && line.charCodeAt(gateP) === 32) gateP++;
        const gateC = gateP < line.length ? line.charCodeAt(gateP) : 0;
        // A Berny parameter table opens with `!    Initial Parameters    !`; its
        // rows are the only place the file says WHICH coordinate a scan drives.
        // The '!' is checked before the capital-letter gate, and the rows are
        // read in one bounded block, so a coordinate row costs a char code test.
        if (gateC === 33 && line.indexOf('Initial Parameters') >= 0) {
            const block = readScanParameters(lines, i);
            if (block.params.length > 0) scanState.params = block.params;
            i = block.endIdx;
            continue;
        }
        // The banner that opens an IRC calculation starts its path numbering
        // afresh. It is printed with leading asterisks, so it is checked before
        // the capital-letter gate; the forward path states its own
        // `Point Number` header right below, which takes precedence.
        if (line.indexOf('Start new reaction path') >= 0) {
            ircPathNumber = 0;
            // A new IRC calculation has its own path numbering AND its own
            // transition state: the geometry printed above the banner is the
            // barrier both halves start from. A run whose card names one
            // direction only would otherwise hide which half it computed.
            ircDirection = (jobIrcDirs && jobIrcDirs.length === 1) ? jobIrcDirs[0] : undefined;
            ircRun++;
            ircRunOpen = true;
            const tsFrame = frames.length - 1;
            if (tsFrame >= jobFrameStart && !ircDrafts.some(d => d.frame === tsFrame)) {
                ircRunTs.push({ run: ircRun, frame: tsFrame });
                // The path's start geometry is one of the structures the card
                // shows, so the input-orientation filter must keep it even when
                // an earlier job of a --Link1-- file printed an input pose.
                kinds.bound[tsFrame] = true;
            }
            if (ircPending === null) ircStartPath = 1;
            i++;
            continue;
        }
        if (gateC >= 65 && gateC <= 90) {

            // --- A relaxed scan: which point the tables below belong to ---
            // `Step number   1 out of a maximum of   57 on scan point  1 out of 32`
            // The line precedes a point's first table, so the point is known by
            // the time that table is read.
            if (gateC === 83 && line.indexOf('Step number') >= 0) {
                const sp = SCAN_STEP_RE.exec(line);
                if (sp) {
                    const point = parseInt(sp[3], 10);
                    const stepNo = parseInt(sp[1], 10);
                    const firstPoint = scanState.cur === null && scanState.list.length === 0;
                    if (!scanState.cur || scanState.cur.point !== point) {
                        if (scanState.cur) scanState.list.push(scanState.cur);
                        scanState.cur = {
                            point, step: stepNo, steps: 0, yes: 0,
                            converged: false, frame: -1
                        };
                    } else {
                        scanState.cur.step = stepNo;
                    }
                    // The geometry a step's table describes is printed ABOVE the
                    // step line (Gaussian prints the structure, its SCF, then
                    // announces the next step). Binding it here — not only when
                    // the table arrives — keeps the structure of a point whose
                    // run was cut off between the two, which would otherwise sit
                    // on screen with no point naming it.
                    const fr = frames.length - 1;
                    if (fr >= 0 && fr >= jobFrameStart) {
                        scanState.cur.frame = fr;
                        const f = frames[fr];
                        f.scanPoint = point;
                        f.scanStep = stepNo;
                        f.optStep = undefined;
                        f.stepLabel = scanFrameLabel(point, stepNo, chargeMultLine);
                        kinds.bound[fr] = true;
                    }
                    // The pose the scan starts from is printed before point 1's
                    // own geometry, so no table ever claims it. It is where the
                    // first point's relaxation began, so it is listed with that
                    // point: opening the file then marks a row and the chart
                    // instead of showing a frame no row describes.
                    if (firstPoint) {
                        for (let k = jobFrameStart; k < frames.length - 1; k++) {
                            if (kinds.input[k] && frames[k].scanPoint === undefined) {
                                frames[k].scanPoint = point;
                                break;
                            }
                        }
                    }
                    const total = parseInt(sp[4], 10);
                    if (total > 0) scanState.total = total;
                    i++;
                    continue;
                }
            }

            // --- The second half of an IRC path: its start geometry comes next ---
            // The file prints no `Point Number` header for the reverse path's
            // starting structure, so the banner itself names the point: it is
            // point 0 of the next path. An explicit header still wins, since the
            // binding below only uses this when nothing is pending.
            if (line.indexOf('Beginning calculation of the') >= 0) {
                const bd = /Beginning calculation of the (FORWARD|REVERSE) path/i.exec(line);
                if (bd) ircDirection = bd[1].toLowerCase() === 'reverse' ? 'reverse' : 'forward';
                // A file that skipped the start banner still opens its run here;
                // the geometry above this line is that run's transition state.
                if (!ircRunOpen) {
                    ircRun++;
                    ircRunOpen = true;
                    const tsFrame = frames.length - 1;
                    if (tsFrame >= jobFrameStart && !ircDrafts.some(d => d.frame === tsFrame)) {
                        ircRunTs.push({ run: ircRun, frame: tsFrame });
                        kinds.bound[tsFrame] = true;
                    }
                }
                if (ircPending === null) ircStartPath = ircPathNumber + 1;
                i++;
                continue;
            }

            // --- End of a path, or of the whole IRC ---
            // A point the run announced and then hit the step limit on never
            // printed a geometry: its pending header must not claim the next
            // path's start structure. The same holds when the calculation ends
            // (normally or not), so all three markers drop the pending point.
            // The announced path start is kept: a file that printed the reverse
            // banner before its own closing line has already named where the
            // next geometry belongs.
            if (line.indexOf('Maximum number of steps reached') >= 0 ||
                (line.indexOf('Calculation of') >= 0 && line.indexOf('path complete') >= 0) ||
                line.indexOf('Reaction path calculation complete') >= 0) {
                // `Calculation of FORWARD/REVERSE path complete.` names the half
                // just finished. A run that printed no direction banner (a
                // forward-only job normally prints none) is labelled from here.
                const cd = /Calculation of (FORWARD|REVERSE) path complete/i.exec(line);
                if (cd) {
                    const dir: IrcDirection = cd[1].toLowerCase() === 'reverse' ? 'reverse' : 'forward';
                    for (let d = ircDrafts.length - 1;
                         d >= 0 && ircDrafts[d].path === ircPathNumber && !ircDrafts[d].dir; d--) {
                        ircDrafts[d].dir = dir;
                    }
                }
                ircPending = null;
                i++;
                continue;
            }

            // --- An IRC path point: the geometry below is that point's ---
            if (line.indexOf('Point Number:') >= 0) {
                const ip = IRC_POINT_RE.exec(line);
                if (ip) {
                    ircPathNumber = parseInt(ip[2], 10);
                    ircStartPath = null;
                    ircPending = { path: ircPathNumber, point: parseInt(ip[1], 10), dir: ircDirection };
                    i++;
                    continue;
                }
            }

            // --- The scan point's own optimisation finished ---
            if (gateC === 79 && line.indexOf('Optimization completed') >= 0) {
                if (scanState.cur) scanState.cur.converged = true;
                i++;
                continue;
            }

            // --- The run's closing scan summary: coordinate and energy per point ---
            if (gateC === 83 && line.indexOf('Summary of') >= 0 && line.indexOf('Surface Scan') >= 0) {
                const table = readScanTable(lines, i);
                if (Object.keys(table.energy).length > 0) scanState.table = table;
                i += 1;
                continue;
            }

            // SCF Done:  E(RHF) =  -76.0107469158     A.U. after   10 cycles
            if (line.indexOf('SCF Done') >= 0) {
                const scfMatch = line.match(SCF_DONE_RE);
                if (scfMatch) {
                    // The geometry of this step was printed just above, so the
                    // energy is bound back to the frame(s) still waiting for one
                    // — and only to this job's frames. The copy in
                    // `pendingEnergy` is what the convergence table below (and
                    // an `opt=external` run's `Energy=` line) consumes.
                    // [2026-11-12 | Atreides-Jimmy]
                    pendingEnergy = parseLogNumber(scfMatch[2]);
                    scfCount++;
                    // E(RPM7), E(RB3LYP), … — the label is the file's own.
                    noteEnergy(jobs, 'E(' + scfMatch[1] + ')', pendingEnergy, scfCount);
                    assignPendingEnergyToFrames(frames, pendingEnergy, jobFrameStart);
                    const cycles = line.match(SCF_CYCLES_RE);
                    if (cycles) summaryOf(jobs).scfCycles = parseInt(cycles[1], 10);
                    i++;
                    continue;
                }
            }

            // --- A per-atom charge table (population analysis) ---
            // ` Mulliken charges:`, ` ESP charges:`, ` APT charges:`,
            // ` Mulliken charges with hydrogens summed into heavy atoms:` —
            // whatever the file printed, by the file's own name. A line that
            // merely mentions charges yields no rows and costs one attempt.
            if (line.indexOf('harges') >= 0) {
                const block = readChargeBlock(lines, i);
                if (block) {
                    noteChargeTable(charges, block.label, block.values, frames.length);
                    i = block.endIdx;
                    continue;
                }
            }

            // --- Orbital eigenvalues: the level diagram's only source ---
            // ` Alpha  occ. eigenvalues --  -77.28310 …` and the `virt.` block
            // that follows it, one pair per spin channel. A log prints these
            // for every SCF, so the reader keeps them all and the last one —
            // the converged geometry's — is what the panel draws.
            if (line.indexOf('eigenvalues --') >= 0) {
                const eig = ORBITAL_EIGEN_RE.exec(line);
                if (eig) {
                    noteOrbitalLine(orbitals, eig[1], eig[2], eig[3], frames.length - 1);
                    i++;
                    continue;
                }
            }

            // --- The solvation cavity the job built ---
            // A block per geometry: an optimisation prints the cavity it builds
            // at every step, and the panel describes the structure on screen, so
            // each banner starts the cavity afresh and the last one wins. The
            // numbers that describe it (type, solvent, sphere count, area,
            // volume) follow within a few dozen lines, so they are read in one
            // window rather than by testing every line for five literals.
            if (line.indexOf('Polarizable Continuum Model') >= 0) {
                clearCavity(cavity);
                noteCavityLines(cavity, lines, i, Math.min(lines.length, i + 90));
            }

            // --- RMS gradient of the final SCF ---
            // Only the full label is taken: `RMS Force = .10134D-02` is the
            // SCF's DIIS convergence test, not the gradient of the structure.
            if (line.indexOf('RMS Gradient Norm') >= 0) {
                const g = line.match(RMS_GRADIENT_RE);
                if (g) {
                    const s = summaryOf(jobs);
                    s.rmsGradient = parseLogNumber(g[1]);
                    s.rmsGradientLabel = 'RMS Gradient Norm';
                    noteScf(jobs, 'rmsGradient', scfCount);
                    i++;
                    continue;
                }
            }

            // --- Dipole moment: the label line, then X / Y / Z / Tot ---
            if (line.indexOf('ipole moment') >= 0 || line.indexOf('ipole Moment') >= 0) {
                const d = (lines[i + 1] || '').match(DIPOLE_RE);
                if (d) {
                    summaryOf(jobs).dipole = {
                        x: parseLogNumber(d[1]), y: parseLogNumber(d[2]),
                        z: parseLogNumber(d[3]), total: parseLogNumber(d[4])
                    };
                    noteScf(jobs, 'dipole', scfCount);
                    i += 2;
                    continue;
                }
            }

            // --- Atoms, basis functions ---
            if (line.indexOf('NAtoms=') >= 0) {
                const n = line.match(NATOMS_RE);
                if (n) { summaryOf(jobs).nAtoms = parseInt(n[1], 10); i++; continue; }
            }
            if (line.indexOf('NBsUse=') >= 0) {
                const n = line.match(NBUSE_RE);
                if (n) { summaryOf(jobs).nBasis = parseInt(n[1], 10); i++; continue; }
            }

            // --- Job cpu time / elapsed time, kept in the file's own words ---
            if (line.indexOf('time:') >= 0) {
                if (line.indexOf('Job cpu time') >= 0) {
                    const t = line.match(JOB_CPU_RE);
                    if (t) { summaryOf(jobs).cpuTime = t[1]; i++; continue; }
                }
                if (line.indexOf('Elapsed time') >= 0) {
                    const t = line.match(ELAPSED_RE);
                    if (t) { summaryOf(jobs).elapsedTime = t[1]; i++; continue; }
                }
            }

            // --- Point group ('Point Group = C1', 'Full point group  C1  NOp 1') ---
            if (line.indexOf('oint group') >= 0) {
                const pg = line.match(POINT_GROUP_RE) || line.match(FULL_POINT_GROUP_RE);
                if (pg) {
                    summaryOf(jobs).pointGroup = pg[1];
                    i++;
                    continue;
                }
            }

            // --- Termination: the one line that says the job finished ---
            if (line.indexOf('ermination') >= 0) {
                if (line.indexOf('Normal termination') >= 0) {
                    summaryOf(jobs).termination = 'normal';
                    ircPending = null;
                    ircStartPath = null;
                    i++;
                    continue;
                }
                if (line.indexOf('Error termination') >= 0) {
                    summaryOf(jobs).termination = 'error';
                    // A run cut off after a `Point Number` header leaves that
                    // header waiting; nothing may bind to it.
                    ircPending = null;
                    ircStartPath = null;
                    i++;
                    continue;
                }
            }
        }

        // External energy (opt=external, e.g. Gaussian + xTB/XO combined runs
        // have no SCF Done; the optimizer prints the external energy as):
        //  Energy=    -9905.57592     NIter=   0.
        // Anchored so lines like "Predicted change in Energy=-5.282751D-01"
        // do not match. The literal test in front of the pattern is what that
        // anchor already requires, so it can only save work.
        if (line.indexOf('Energy=') >= 0) {
            const extEnergyMatch = line.match(/^\s*Energy=\s*(-?\d+\.\d+)\s+NIter=/);
            if (extEnergyMatch) {
                pendingEnergy = parseFloat(extEnergyMatch[1]);
                // [2026-09-12 | Atreides-Jimmy] opt=external logs print the cycle's
                // geometry BEFORE its energy — bind it here (same purpose as above),
                // to this job's frames only.
                assignPendingEnergyToFrames(frames, pendingEnergy, jobFrameStart);
                i++;
                continue;
            }
        }

        // Convergence criteria block starts with "Maximum Force" data line
        if (line.includes('Maximum Force') && /-?\d+\.\d/.test(line)) {
            // The series belongs to the job that printed the table: only a job
            // whose route asks for an optimization runs Berny's step loop (see
            // `jobHasOpt`). A job that stated no card at all — a log whose head
            // is gone — cannot be classified, and its tables are taken, because
            // a file that plainly holds a Berny run is an optimization whatever
            // its first lines say.
            const take = jobHasOpt || !jobCardSeen;
            // A scan point's table is that point's step, not a step of one
            // series: eight points of four steps each are eight relaxations of
            // four steps, never a 32-step optimisation. The point takes its own
            // table (and the geometry that table describes) and the file's
            // series stays empty, which is what tells the page to draw the scan
            // instead of a convergence curve.
            const inScan = take && scanState.cur !== null;
            const step: OptStep | null = take
                ? { step: inScan && scanState.cur ? scanState.cur.step : stepCounter + 1, energy: pendingEnergy }
                : null;
            if (step && !inScan) stepCounter++;
            pendingEnergy = undefined;
            let yes = 0;

            // Parse this block: 4 lines (Maximum Force, RMS Force, Maximum Displacement, RMS Displacement)
            for (let k = 0; k < 4 && i < lines.length; k++, i++) {
                if (!step) continue;    // a non-optimizing job's table: consume it whole
                const cur = lines[i];
                const valMatch = cur.match(/(-?\d+\.\d+(?:[EDed][-+]?\d+)?)/);
                if (!valMatch) continue;
                const val = parseFloat(valMatch[1].replace(/[EDed]/i, 'e'));
                if (cur.includes('Maximum Force') && !cur.includes('RMS')) {
                    step.maxForce = val;
                } else if (cur.includes('RMS') && cur.includes('Force')) {
                    step.rmsForce = val;
                } else if (cur.includes('Maximum Displacement')) {
                    step.maxDisplacement = val;
                } else if (cur.includes('RMS') && cur.includes('Displacement')) {
                    step.rmsDisplacement = val;
                }
                // The file's own verdict on that row, used for a scan point
                // whose run was cut off before it could print "completed".
                if (cur.indexOf('YES') >= 0) yes++;
            }
            if (step) {
                const cur = scanState.cur;
                if (inScan && cur) {
                    cur.steps++;
                    cur.yes = yes;
                    cur.last = step;
                    cur.step = step.step;
                    const fr = frames.length - 1;
                    if (fr >= 0) {
                        cur.frame = fr;
                        const f = frames[fr];
                        f.scanPoint = cur.point;
                        f.scanStep = step.step;
                        f.optStep = undefined;
                        f.stepLabel = scanFrameLabel(cur.point, step.step, chargeMultLine);
                        kinds.bound[fr] = true;
                    }
                } else {
                    optSteps.push(step);
                }
            }
            continue;
        }

        i++;
    }

    // If the log contains both "Input orientation" and "Standard orientation" frames,
    // keep only the first "Input" frame (initial geometry) and skip all subsequent ones.
    // "Input orientation" of step N+1 is the same molecular structure as "Standard orientation"
    // of step N, just in a different coordinate system (original input vs principal axes).
    // This creates visually redundant frames in the 3D viewer (which normalizes orientation).
    //
    // The role is read from `kinds`, not from the label: a scan point's geometry
    // is labelled by its point and an IRC point's by its path, and neither says
    // which block it came from. A frame a scan point or an IRC point claims is
    // kept whatever its block is — it is that point's only structure.
    const hasInput = kinds.input.indexOf(true) >= 0;
    const hasStandard = kinds.input.indexOf(false) >= 0;
    // Where each frame ended up after the filter, so a charge table recorded
    // against a pre-filter position still names the frame it describes.
    let frameRemap: number[] | null = null;
    if (hasInput && hasStandard) {
        let inputSeen = false;
        const kept: LogFrame[] = [];
        const remap: number[] = [];
        for (let k = 0; k < frames.length; k++) {
            const f = frames[k];
            if (kinds.input[k] && !kinds.bound[k]) {
                if (inputSeen) { remap.push(-1); continue; }
                inputSeen = true;
            }
            remap.push(kept.length);
            kept.push(f);
        }
        frames.length = 0;
        frames.push(...kept);
        frameRemap = remap;
    }

    if (frames.length === 0 && fallbackFrame) {
        frames.push(fallbackFrame);
    }

    // The last scan point has no following point to close it.
    if (scanState.cur) {
        scanState.list.push(scanState.cur);
        scanState.cur = null;
    }
    // A finished scan or IRC prints its converged geometry once more after its
    // summary, with no table of its own: the block repeats the last point's
    // structure and would end the frame list on a copy of the frame above it.
    // It is dropped, so the list is exactly the steps the run took. Only an
    // unbound trailing frame can go — every frame a point claims stays, and no
    // link can name a frame that is not there.
    //
    // The count before the drop is kept: a charge table or an orbital dump
    // recorded against one of these frames must still be addressable after it.
    const preDropFrameCount = frames.length;
    if ((scanState.list.length > 0 || ircDrafts.length > 0) && frames.length > 1) {
        const lastIdx = frames.length - 1;
        const last = frames[lastIdx];
        let claimed = last.scanPoint != null;
        const preCount = frameRemap ? frameRemap.length : frames.length;
        for (let d = 0; !claimed && d < ircDrafts.length; d++) {
            if (remapFrameIndex(frameRemap, ircDrafts[d].frame, preCount) === lastIdx) claimed = true;
        }
        if (!claimed && sameGeometry(last.atoms, frames[lastIdx - 1].atoms)) {
            frames.pop();
        }
    }
    // The two coordinate-driven analyses, in the frame list's own indices. A
    // file that holds neither leaves both undefined, which is what keeps the
    // convergence panel on an ordinary optimisation.
    const scan = buildScanResult(scanState, frames, frameRemap);
    const irc = buildIrcResult(ircDrafts, frames, frameRemap, ircRunTs);
    const pathways: LogPathways | undefined = (scan || irc) ? { scan, irc } : undefined;

    // Population analyses: a table describes the geometry printed before it
    // (Gaussian closes a step's SCF with its charges, then prints the next
    // geometry), so it takes the last frame that existed when the table ended —
    // and the frame at the end of a single-point job, where no geometry follows
    // the analysis at all. Values are laid out by atom index, padded to the
    // frame's atom list so a view can read one by index without a second table.
    //
    // The anchor is a PRE-FILTER index: the table was recorded while every
    // `Input orientation` block was still on the list. It must therefore be
    // clamped to the pre-filter count and carried through the filter BEFORE it
    // is clamped to the final list — clamping to the shorter post-filter list
    // first landed on an index that had been dropped, backed the lookup up to a
    // mid-file `Standard orientation`, and labelled the table with a geometry
    // from the middle of the run.
    const chargeSets: LogChargeSet[] = [];
    if (frames.length > 0) {
        const chargePreCount = frameRemap ? frameRemap.length : preDropFrameCount;
        /** The job each table came from, parallel to `chargeSets`. */
        const chargeJobs: number[] = [];
        for (const draft of charges.drafts) {
            let fr = draft.anchor > 0 ? draft.anchor - 1 : 0;
            if (fr > chargePreCount - 1) fr = chargePreCount - 1;
            const mapped = remapFrameIndex(frameRemap, fr, chargePreCount);
            fr = mapped === undefined ? 0 : mapped;
            // The trailing copy a scan/IRC prints after its summary is gone by
            // now; a table anchored on it belongs to the structure before it.
            if (fr > frames.length - 1) fr = frames.length - 1;
            const frameAtoms = frames[fr].atoms;
            const values: (number | null)[] = [];
            let covered = 0;
            for (let k = 0; k < frameAtoms.length; k++) {
                const v = k < draft.values.length ? draft.values[k] : null;
                values.push(v == null ? null : v);
                if (v != null) covered++;
            }
            chargeSets.push({
                label: draft.label,
                values,
                elements: frameAtoms.map(a => a.element),
                covered,
                complete: frameAtoms.length > 0 && covered === frameAtoms.length,
                frame: fr
            });
            chargeJobs.push(draft.job);
        }
        // The panel describes the structure the file ends on, so the LAST job's
        // tables are the ones its dropdown opens on; the earlier jobs' tables
        // follow, each still naming its own frame. Within one job the printed
        // order is kept. A `--Link1--` file that analysed two jobs therefore
        // offers both schemes, and a scheme that appears in two jobs says which
        // one it belongs to instead of showing twice under one name.
        if (chargeSets.length > 1) {
            const order = chargeSets.map((_s, k) => k);
            order.sort((a, b) => (chargeJobs[b] - chargeJobs[a]) || (a - b));
            const sortedSets = order.map(k => chargeSets[k]);
            const sortedJobs = order.map(k => chargeJobs[k]);
            chargeSets.length = 0;
            for (const s of sortedSets) chargeSets.push(s);
            chargeJobs.length = 0;
            for (const j of sortedJobs) chargeJobs.push(j);
        }
        // Two jobs can print the same scheme for different structures; the
        // label says which job it belongs to rather than showing twice under
        // one name. A single-job file keeps the file's own wording.
        const distinctJobs = new Set(chargeJobs);
        if (distinctJobs.size > 1) {
            const counts = new Map<string, number>();
            for (const s of chargeSets) counts.set(s.label, (counts.get(s.label) || 0) + 1);
            for (let k = 0; k < chargeSets.length; k++) {
                if ((counts.get(chargeSets[k].label) || 0) > 1) {
                    chargeSets[k].label = chargeSets[k].label + ' (job ' + Math.max(1, chargeJobs[k]) + ')';
                }
            }
        }
    }

    // The orbitals the level diagram shows are the last SCF's, so they belong
    // to the frame that SCF was computed at. Its index takes the same mapping
    // as a charge table; a checkpoint's orbitals carry no frame.
    const levels = orbitalLevelsOf(orbitals);
    if (levels && typeof levels.frame === 'number') {
        const preCount = frameRemap ? frameRemap.length : preDropFrameCount;
        const mapped = remapFrameIndex(frameRemap, levels.frame, preCount);
        levels.frame = mapped === undefined ? undefined : Math.min(mapped, frames.length - 1);
    }

    // An IRC's frame list is its reaction path: reverse far end, transition
    // state, forward far end. The move happens after every frame-bound record
    // is built so the path points, the run's TS, the charge tables and the
    // orbitals can be carried across with it.
    orderIrcFrames(frames, irc, scan, chargeSets, levels);

    // Frames a path claims are named by the row that lists them (`TS`, `R24`,
    // `F0`), not by the block they were printed in. This runs after the move,
    // because the indices are the final ones.
    labelIrcFrames(frames, irc);

    const summary = jobs.list.filter(summaryHasData);
    // A cavity object with nothing but its source is not a statement about the
    // file: no solvation model ran, so the panel must not claim one.
    const hasCavity = cavity.type !== undefined || cavity.solvent !== undefined ||
        cavity.sphereCount !== undefined || cavity.area !== undefined || cavity.volume !== undefined;
    // The banner calls every continuum job PCM, SMD included; when the route
    // names the model, its word is the more specific one and wins.
    if (hasCavity) {
        const routeModel = scrfModelOf(content);
        if (routeModel) cavity.model = routeModel;
    }

    onProgress?.(1);
    return {
        frames,
        title,
        // A table is only taken from a job that asked to optimize (or from a
        // file that states no route at all), so a non-empty list is already the
        // answer — there is no file-level opt flag left to consult.
        optSteps: optSteps.length > 0 ? optSteps : undefined,
        normalModes: modes,
        routes: routes.length > 0 ? routes : undefined,
        charges: chargeSets.length > 0 ? chargeSets : undefined,
        summary: summary.length > 0 ? summary : undefined,
        orbitals: levels,
        cavity: hasCavity ? cavity : undefined,
        pathways
    };
}
