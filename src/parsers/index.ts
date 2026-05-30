import { MolecularData, AtomGroup } from '../types';
import { parseGjf } from './gjfParser';
import { parseXyz } from './xyzParser';
import { parseMol2 } from './mol2Parser';
import { parseGaussianLog, LogFrame } from './logParser';
import { parseCoord } from './coordParser';
import { parseOrcainp } from './orcaInpParser';
import { parseOrcaOut, OrcaFrame } from './orcaOutParser';
import { parsePdb } from './pdbParser';
import { parseMopac } from './mopacParser';
import { parseTcl, TclParseResult } from './tclParser';

export { parseGaussianLog, LogFrame };
export { parseOrcaOut, OrcaFrame };
export { parseTcl, TclParseResult };

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
        default:
            return tryAutoParse(content);
    }
}

export function parseLogFile(content: string, fileName?: string): { frames: LogFrame[] | OrcaFrame[], title: string } {
    const ext = (fileName || '').toLowerCase().split('.').pop() || '';
    if (ext === 'out') {
        const result = parseOrcaOut(content);
        return { frames: result.frames, title: result.title };
    }
    return parseGaussianLog(content);
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

function tryAutoParse(content: string): MolecularData {
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
