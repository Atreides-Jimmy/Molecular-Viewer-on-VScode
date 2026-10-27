import { MolecularData, AtomGroup, OptStep, NormalMode, OptThresholds } from '../types';
import { parseGjf } from './gjfParser';
import { parseXyz, parseXyzTrajectory, XyzFrame, XyzTrajectory } from './xyzParser';
import { parseMol2 } from './mol2Parser';
import { parseGaussianLog, LogFrame, GaussianLogResult, RouteSection, RouteKeyword } from './logParser';
import { parseCoord } from './coordParser';
import { parseOrcainp } from './orcaInpParser';
import { parseOrcaOut, isOrcaOutput, OrcaFrame, OrcaLogResult } from './orcaOutParser';
import { parsePdb } from './pdbParser';
import { parseMopac } from './mopacParser';
import { parseCif } from './cifParser';
import { parseVasp } from './vaspParser';
import { parseCube } from './cubeParser';
import { parseVesta } from './vestaParser';
import { parseXtbLog, isXtbLog, XtbLogResult } from './xtbParser';
import { parseTcl, TclParseResult, resolveVmdScene, VmdResolvedScene, hasScalarChannel, resolveSelection, vdwRadius, parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg } from './vmdVisualization';

export { parseGaussianLog, LogFrame, GaussianLogResult };
export { RouteSection, RouteKeyword };
export { parseOrcaOut, isOrcaOutput, OrcaFrame, OrcaLogResult };
export { parseXyzTrajectory, XyzFrame, XyzTrajectory };
export { parseCif };
export { parseVasp };
export { parseCube };
export { parseVesta };
export { parseXtbLog, isXtbLog, XtbLogResult };
// Visualization files: the VMD scripts and the bare scalar-carrying coordinate
// files they describe (see src/parsers/vmdVisualization.ts).
export { parseTcl, TclParseResult };
export { resolveVmdScene, VmdResolvedScene };
export { hasScalarChannel };
export { resolveSelection, vdwRadius };
export { parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg };

export interface LogFileResult {
    frames: LogFrame[] | OrcaFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    // Route-card sections (Gaussian's keywords lines, ORCA's simple-input line),
    // shown in the webview's Route panel. Other formats leave it undefined.
    routes?: RouteSection[];
    // The convergence targets the file itself states, when it states them.
    // The panel draws them as dashed target lines; a program that states none
    // (xtb's gnorm criterion, a trajectory with no table) gets no lines rather
    // than a borrowed scale.
    optThresholds?: OptThresholds;
    // Which log format produced this result. The webview uses it to label the
    // panel; the numeric scales come from `optThresholds`.
    source?: 'gaussian' | 'xtb' | 'orca';
}

export function parseFile(content: string, fileName: string): MolecularData {
    const ext = fileName.toLowerCase().split('.').pop() || '';

    switch (ext) {
        case 'gjf':
        case 'gjf03':
        case 'gjf09':
        case 'gjf16':
        case 'com':
            return parseGjf(content);
        case 'xyz':
            return parseXyz(content);
        case 'mol2':
            return parseMol2(content);
        case 'log':
        case 'out':
            return parseLogAsSingleFrame(content);
        case 'coord':
            return parseCoord(content);
        case 'inp': {
            // `.inp` is ORCA's input extension, but Gaussian users name their
            // inputs that way too. ORCA's own syntax is unambiguous (`* xyz` /
            // `* xyzfile` / `%coords`); anything else that starts with a route
            // card is read as a Gaussian input instead of coming back empty.
            const orca = parseOrcainp(content);
            if (orca.atoms.length > 0) return orca;
            if (/^\s*#/m.test(content)) return parseGjf(content);
            return orca;
        }
        case 'pdb':
        case 'ent':
            return parsePdb(content);
        case 'mop':
        case 'mopac':
        case 'dat':
            return parseMopac(content);
        case 'cif':
            return parseCif(content);
        case 'vasp':
        case 'poscar':
        case 'contcar':
            return parseVasp(content);
        case 'cube':
            return parseCube(content);
        case 'vesta':
            return parseVesta(content);
        default:
            return tryAutoParse(content, fileName);
    }
}

/**
 * Parse a log/trajectory file into the frames and the analysis data the viewer
 * needs.
 *
 * The format is decided by the CONTENT, not by the extension: Gaussian and ORCA
 * are both routinely told to write `.out`, and xTB shares `.log` with Gaussian.
 * Every detector below matches a string only one program writes, so the sniffing
 * cannot pick the wrong parser; the extension is not consulted at all, which
 * means a mis-named file (an ORCA run saved as `.log`, a Gaussian log saved as
 * `.out`) still opens correctly.
 */
export function parseLogFile(content: string, fileName?: string): LogFileResult {
    if (isOrcaOutput(content)) {
        const result = parseOrcaOut(content);
        return {
            frames: result.frames,
            title: result.title,
            optSteps: result.optSteps,
            normalModes: result.normalModes,
            routes: result.routes,
            optThresholds: result.optThresholds,
            source: 'orca'
        };
    }
    if (isXtbLog(content)) {
        const xResult = parseXtbLog(content);
        return {
            frames: xResult.frames,
            title: xResult.title,
            optSteps: xResult.optSteps,
            normalModes: xResult.normalModes,
            source: 'xtb'
        };
    }
    const gResult = parseGaussianLog(content);
    return {
        frames: gResult.frames,
        title: gResult.title,
        optSteps: gResult.optSteps,
        normalModes: gResult.normalModes,
        routes: gResult.routes,
        source: 'gaussian'
    };
}

function parseLogAsSingleFrame(content: string): MolecularData {
    const result = parseLogFile(content);
    if (result.frames.length > 0) {
        const f = result.frames[0];
        return {
            atoms: f.atoms,
            bonds: f.bonds,
            title: f.title,
            hasExplicitBonds: f.hasExplicitBonds,
            charge: f.charge,
            multiplicity: f.multiplicity
        };
    }
    return { atoms: [], bonds: [], title: 'Empty', hasExplicitBonds: false };
}

export function tryAutoParse(content: string, fileName: string = ''): MolecularData {
    const lines = content.split(/\r?\n/).filter(l => l.trim() !== '');

    if (content.includes('$coord')) {
        return parseCoord(content);
    }

    if (content.match(/\*\s*xyz/i) || content.match(/\*\s*xyzfile/i)) {
        return parseOrcainp(content);
    }

    if (content.match(/^(ATOM|HETATM)/m)) {
        return parsePdb(content);
    }

    if (content.match(/CARTESIAN COORDINATES/i) && content.match(/MOPAC/i)) {
        return parseMopac(content);
    }

    if (content.match(/^data_/m) && content.includes('_cell_length_a')) {
        return parseCif(content);
    }

    if (fileName.match(/\.vasp$/i) || fileName.match(/poscar/i) || fileName.match(/contcar/i)) {
        return parseVasp(content);
    }

    if (fileName.match(/\.cube$/i)) {
        return parseCube(content);
    }

    if (fileName.match(/\.vesta$/i) || content.includes('#VESTA_FORMAT_VERSION')) {
        return parseVesta(content);
    }

    if (isOrcaOutput(content)) {
        const result = parseOrcaOut(content);
        if (result.frames.length > 0) {
            return {
                atoms: result.frames[0].atoms,
                bonds: result.frames[0].bonds,
                title: result.frames[0].title,
                hasExplicitBonds: result.frames[0].hasExplicitBonds,
                charge: result.charge,
                multiplicity: result.multiplicity
            };
        }
    }

    if (lines.length > 0) {
        const firstLine = lines[0].trim();
        const possibleCount = parseInt(firstLine, 10);
        if (!isNaN(possibleCount) && possibleCount > 0 && possibleCount < 100000) {
            return parseXyz(content);
        }
    }

    if (content.includes('@<TRIPOS>')) {
        return parseMol2(content);
    }

    if (content.includes('Standard orientation:') || content.includes('Input orientation:')) {
        return parseLogAsSingleFrame(content);
    }

    if (content.includes('--Link1--') || content.match(/#\s*[A-Za-z]/)) {
        return parseGjf(content);
    }

    return parseXyz(content);
}
