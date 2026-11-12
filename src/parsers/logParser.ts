import { Atom, Bond, MolecularData, OptStep, NormalMode, SolventCavity, OrbitalLevels, OrbitalSpinLevels } from '../types';

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
}

/** The charge tables of the job being read. */
interface ChargeState {
    drafts: ChargeDraft[];
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
 *  one, which is the geometry the run ended on. */
function noteChargeTable(st: ChargeState, label: string, values: (number | null)[], anchor: number): void {
    for (let k = 0; k < st.drafts.length; k++) {
        if (st.drafts[k].label === label) {
            st.drafts[k] = { label: label, values: values, anchor: anchor };
            return;
        }
    }
    st.drafts.push({ label: label, values: values, anchor: anchor });
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
function noteOrbitalLine(st: OrbitalState, spin: string, kind: string, text: string): void {
    const tokens = text.trim().split(/\s+/);
    const values: number[] = [];
    for (let k = 0; k < tokens.length; k++) {
        const v = parseLogNumber(tokens[k]);
        if (isFinite(v)) values.push(v);
    }
    if (values.length === 0) return;
    if (!st.cur) st.cur = { alpha: [], alphaOcc: -1, beta: [], betaOcc: -1 };
    if (spin === 'Alpha' && kind === 'occ' && st.cur.alphaOcc >= 0) {
        st.last = st.cur;
        st.cur = { alpha: [], alphaOcc: -1, beta: [], betaOcc: -1 };
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
    return { spins, closedShell: closed };
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
    const charges: ChargeState = { drafts: [] };
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
                    // A route card is the start of a job: the results and the
                    // population analysis that follow belong to it, not to the
                    // job whose card preceded them. This is the only boundary a
                    // `--Link1--` file states, and the first card of a file
                    // finds nothing to close.
                    jobs.cur = null;
                    charges.drafts.length = 0;
                    orbitals.last = null;
                    orbitals.cur = null;
                    clearCavity(cavity);
                    scfCount = 0;
                    // A job that ended without a convergence table (a single
                    // point) leaves its SCF energy pending; carrying it into the
                    // next job would stamp it onto that job's first geometry.
                    pendingEnergy = undefined;
                    jobFrameStart = frames.length;
                    i = j;
                    continue;
                }
            }
        }

        // --- Harmonic frequencies section (parsed inline; the LAST section wins) ---
        if (line.includes('Harmonic frequencies') && line.includes('cm**-1')) {
            const res = parseModesFromHeader(lines, i);
            modes = res.modes.length > 0 ? res.modes : undefined;
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
                const label = isStandard
                    ? `Step ${stdFrames + 1}${chargeMultLine ? ' (' + chargeMultLine + ')' : ''}`
                    : `Input${chargeMultLine ? ' (' + chargeMultLine + ')' : ''}`;
                frames.push({
                    atoms,
                    bonds: [],
                    title: label,
                    hasExplicitBonds: false,
                    stepLabel: label,
                    charge: logCharge,
                    multiplicity: logMultiplicity,
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
        if (gateC >= 65 && gateC <= 90) {

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
                    noteOrbitalLine(orbitals, eig[1], eig[2], eig[3]);
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
                    i++;
                    continue;
                }
                if (line.indexOf('Error termination') >= 0) {
                    summaryOf(jobs).termination = 'error';
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
            const step: OptStep | null = take ? { step: stepCounter + 1, energy: pendingEnergy } : null;
            if (step) stepCounter++;
            pendingEnergy = undefined;

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
            }
            if (step) optSteps.push(step);
            continue;
        }

        i++;
    }

    // If the log contains both "Input orientation" and "Standard orientation" frames,
    // keep only the first "Input" frame (initial geometry) and skip all subsequent ones.
    // "Input orientation" of step N+1 is the same molecular structure as "Standard orientation"
    // of step N, just in a different coordinate system (original input vs principal axes).
    // This creates visually redundant frames in the 3D viewer (which normalizes orientation).
    const hasInput = frames.some(f => f.stepLabel.startsWith('Input'));
    const hasStandard = frames.some(f => f.stepLabel.startsWith('Step'));
    // Where each frame ended up after the filter, so a charge table recorded
    // against a pre-filter position still names the frame it describes.
    let frameRemap: number[] | null = null;
    if (hasInput && hasStandard) {
        let inputSeen = false;
        const kept: LogFrame[] = [];
        const remap: number[] = [];
        for (let k = 0; k < frames.length; k++) {
            const f = frames[k];
            if (f.stepLabel.startsWith('Input')) {
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

    // Population analyses: a table describes the geometry printed before it
    // (Gaussian closes a step's SCF with its charges, then prints the next
    // geometry), so it takes the last frame that existed when the table ended —
    // and the frame at the end of a single-point job, where no geometry follows
    // the analysis at all. Values are laid out by atom index, padded to the
    // frame's atom list so a view can read one by index without a second table.
    const chargeSets: LogChargeSet[] = [];
    if (frames.length > 0) {
        for (const draft of charges.drafts) {
            let fr = draft.anchor > 0 ? draft.anchor - 1 : 0;
            if (fr > frames.length - 1) fr = frames.length - 1;
            if (frameRemap) {
                while (fr > 0 && frameRemap[fr] < 0) fr--;
                fr = frameRemap[fr] < 0 ? 0 : frameRemap[fr];
            }
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
        }
    }

    const summary = jobs.list.filter(summaryHasData);
    const levels = orbitalLevelsOf(orbitals);
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
        cavity: hasCavity ? cavity : undefined
    };
}
