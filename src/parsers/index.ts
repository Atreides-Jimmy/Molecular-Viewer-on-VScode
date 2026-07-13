import { MolecularData, AtomGroup, OptStep, NormalMode } from '../types';
import { parseGjf } from './gjfParser';
import { parseXyz } from './xyzParser';
import { parseMol2 } from './mol2Parser';
import { parseGaussianLog, LogFrame, GaussianLogResult } from './logParser';
import { parseCoord } from './coordParser';
import { parseOrcainp } from './orcaInpParser';
import { parseOrcaOut, OrcaFrame } from './orcaOutParser';
import { parsePdb } from './pdbParser';
import { parseMopac } from './mopacParser';
import { parseTcl, TclParseResult } from './tclParser';
import { parseCif } from './cifParser';
import { parseVasp } from './vaspParser';
import { parseCube } from './cubeParser';
import { parseVesta } from './vestaParser';

export { parseGaussianLog, LogFrame, GaussianLogResult };
export { parseOrcaOut, OrcaFrame };
export { parseTcl, TclParseResult };
export { parseCif };
export { parseVasp };
export { parseCube };
export { parseVesta };

export interface LogFileResult {
    frames: LogFrame[] | OrcaFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
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
        case 'inp':
            return parseOrcainp(content);
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

export function parseLogFile(content: string, fileName?: string): LogFileResult {
    const ext = (fileName || '').toLowerCase().split('.').pop() || '';
    if (ext === 'out') {
        const result = parseOrcaOut(content);
        return { frames: result.frames, title: result.title };
    }
    const gResult = parseGaussianLog(content);
    return {
        frames: gResult.frames,
        title: gResult.title,
        optSteps: gResult.optSteps,
        normalModes: gResult.normalModes
    };
}

function parseLogAsSingleFrame(content: string): MolecularData {
    const result = parseGaussianLog(content);
    if (result.frames.length > 0) {
        return {
            atoms: result.frames[0].atoms,
            bonds: result.frames[0].bonds,
            title: result.frames[0].title,
            hasExplicitBonds: result.frames[0].hasExplicitBonds
        };
    }
    return { atoms: [], bonds: [], title: 'Empty', hasExplicitBonds: false };
}

function tryAutoParse(content: string, fileName: string = ''): MolecularData {
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

    if (content.includes('CARTESIAN COORDINATES (ANGSTROEM)')) {
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
