import { Atom, Bond, MolecularData, OptStep, NormalMode } from '../types';

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
    101: 'Md', 102: 'No', 103: 'Lr'
};

export interface LogFrame {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    stepLabel: string;
    charge?: number;
    multiplicity?: number;
}

export interface GaussianLogResult {
    frames: LogFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    routes?: RouteSection[];
}

function skipDashedLines(lines: string[], startIdx: number, count: number): number {
    let i = startIdx;
    let skipped = 0;
    while (i < lines.length && skipped < count) {
        if (lines[i].includes('---')) {
            skipped++;
        }
        i++;
    }
    return i;
}

/** A route card starts with '#' optionally followed by a single-letter print
 *  level (#p/#P/#t/#n, case-insensitive) or a bare '#'; the keyword part must
 *  follow. '#pop'-style strings are not route cards. */
const ROUTE_START_RE = /^#(\s|[A-Za-z]\b|$)/;
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
            // Keyword with a bare parenthesized option list, e.g. iop(5/17=17)
            keywords.push({
                name: t.slice(0, pi),
                options: t.slice(pi + 1, -1).split(',').map(s => s.trim()).filter(s => s !== '')
            });
            continue;
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

export function parseGaussianLog(content: string): GaussianLogResult {
    const lines = content.split(/\r?\n/);
    const frames: LogFrame[] = [];
    let title = '';

    const titleMatch = content.match(/#\s*[A-Za-z]/);
    if (titleMatch) {
        const routeIdx = content.indexOf(titleMatch[0]);
        const afterRoute = content.substring(0, routeIdx);
        const titleLines = afterRoute.split(/\r?\n/).filter(l => l.trim() !== '' && !l.trim().startsWith('%'));
        if (titleLines.length > 0) title = titleLines[titleLines.length - 1].trim();
    }

    const chargeMultMatch = content.match(/Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)/);
    let chargeMultLine = '';
    let logCharge: number | undefined;
    let logMultiplicity: number | undefined;
    if (chargeMultMatch) {
        chargeMultLine = `Charge=${chargeMultMatch[1]} Mult=${chargeMultMatch[2]}`;
        logCharge = parseInt(chargeMultMatch[1], 10);
        logMultiplicity = parseInt(chargeMultMatch[2], 10);
    }

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
    // running an optimizer. The reliable trigger is the route card: emit
    // optSteps only when a '#...' route section in the file actually requests
    // an optimization (opt / optimization / freq=calcall). All route sections
    // are collected (parsed into RouteSection[]) for the Route panel.
    let routeHasOpt = false;
    const routes: RouteSection[] = [];

    let i = 0;
    while (i < lines.length) {
        const line = lines[i];

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
                } else if (ROUTE_START_RE.test(tr)) {
                    routeStart = i;
                }
                if (routeStart >= 0) {
                    let j = routeStart + 1;
                    while (j < lines.length) {
                        const t2 = lines[j].trim();
                        if (t2 === '' || DASH_RE.test(t2) || ROUTE_START_RE.test(t2)) break;
                        j++;
                    }
                    const sec = parseRouteSection(lines.slice(routeStart, j).map(l => l.trim()).join(' '));
                    routes.push(sec);
                    if (sec.hasOpt) routeHasOpt = true;
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
                    multiplicity: logMultiplicity
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

        // --- Convergence table steps ---
        // SCF Done:  E(RHF) =  -76.0107469158     A.U. after   10 cycles
        const scfMatch = line.match(/SCF Done:\s+E\([^)]+\)\s*=\s*(-?\d+\.\d+)/);
        if (scfMatch) {
            pendingEnergy = parseFloat(scfMatch[1]);
            i++;
            continue;
        }

        // External energy (opt=external, e.g. Gaussian + xTB/XO combined runs
        // have no SCF Done; the optimizer prints the external energy as):
        //  Energy=    -9905.57592     NIter=   0.
        // Anchored so lines like "Predicted change in Energy=-5.282751D-01"
        // do not match.
        const extEnergyMatch = line.match(/^\s*Energy=\s*(-?\d+\.\d+)\s+NIter=/);
        if (extEnergyMatch) {
            pendingEnergy = parseFloat(extEnergyMatch[1]);
            i++;
            continue;
        }

        // Convergence criteria block starts with "Maximum Force" data line
        if (line.includes('Maximum Force') && /-?\d+\.\d/.test(line)) {
            const step: OptStep = { step: stepCounter + 1, energy: pendingEnergy };
            stepCounter++;
            pendingEnergy = undefined;

            // Parse this block: 4 lines (Maximum Force, RMS Force, Maximum Displacement, RMS Displacement)
            for (let k = 0; k < 4 && i < lines.length; k++, i++) {
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
            optSteps.push(step);
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
    if (hasInput && hasStandard) {
        let inputSeen = false;
        const filtered = frames.filter(f => {
            if (f.stepLabel.startsWith('Input')) {
                if (inputSeen) return false;
                inputSeen = true;
                return true;
            }
            return true;
        });
        frames.length = 0;
        frames.push(...filtered);
    }

    if (frames.length === 0 && fallbackFrame) {
        frames.push(fallbackFrame);
    }

    return {
        frames,
        title,
        optSteps: routeHasOpt && optSteps.length > 0 ? optSteps : undefined,
        normalModes: modes,
        routes: routes.length > 0 ? routes : undefined
    };
}
