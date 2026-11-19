// Load and parse one structure file the way the VMD scene and the Diff/Import
// paths consume it.
import * as vscode from 'vscode';
import * as fs from 'fs';
import { MolecularData } from '../types';
import {
    detectBinaryTrajectory, parseFile, parseLogFile, parseStructuredProgramOutput,
    parseFchk, isTrajectoryFileName, looksLikeBinaryTrajectory, FchkResult
} from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { LoadedStructure } from './payloads';
import { buildParseOptions, applySetupWarnings } from './parseSetup';
import { loadBinaryTrajectoryDocument } from './trajectoryDocuments';
import { decodeText } from './loading';

/** A checkpoint read through the workspace API is buffered whole; refuse the ones that cannot be. */
const FCHK_BUFFER_LIMIT = 256 * 1024 * 1024;

/**
 * Read a formatted checkpoint without turning it into a string.
 *
 * `loadStructureFile` follows the VMD/Diff/Import path and only needs the
 * structure, but a checkpoint can be most of a gigabyte: it is streamed from
 * disk exactly like the direct-open reader, and only the small arrays are kept.
 * A non-local URI has no path `fs` can open, so it is read through the
 * workspace API and refused above the buffer limit rather than allocating the
 * whole file.
 */
export async function loadFchkStructure(
    uri: vscode.Uri,
    fileName: string,
    report?: (fraction: number, message?: string) => void
): Promise<MolecularData> {
    let size = 0;
    try {
        size = (await vscode.workspace.fs.stat(uri)).size || 0;
    } catch {
        size = 0;
    }
    let result: FchkResult;
    if (uri.scheme === 'file') {
        let lastPercent = -1;
        result = await parseFchk(fs.createReadStream(uri.fsPath, { highWaterMark: 1 << 16 }), {
            size,
            path: uri.fsPath,
            onProgress: (read) => {
                if (!report || size <= 0) return;
                const percent = Math.min(99, Math.floor((read / size) * 100));
                if (percent === lastPercent) return;
                lastPercent = percent;
                report(percent / 100, 'Reading the checkpoint · ' + percent + '%');
            }
        });
    } else {
        if (size > FCHK_BUFFER_LIMIT) {
            throw new Error('this formatted checkpoint is ' + Math.round(size / (1024 * 1024)) +
                ' MB and is not a local file: it can only be read through the workspace file API, ' +
                'which would hold all of it in memory at once');
        }
        const bytes = await vscode.workspace.fs.readFile(uri);
        result = await parseFchk((async function* () { yield bytes; })(), {
            size: bytes.length,
            path: ''
        });
    }
    report?.(1, 'Preparing the checkpoint');
    const data = ensureBonds({
        atoms: result.atoms,
        bonds: [],
        title: /^molecule name$/i.test(result.title.trim()) ? fileName : result.title,
        hasExplicitBonds: false,
        charge: result.charge,
        multiplicity: result.multiplicity,
        warnings: result.warnings.length > 0 ? result.warnings.slice() : undefined
    });
    data.filePath = uri.fsPath;
    return data;
}

/** Load and parse one structure file (any supported format). */
export async function loadStructureFile(uri: vscode.Uri): Promise<LoadedStructure> {
    const fileName = uri.path.split('/').pop() || 'unknown.xyz';
    const ext = fileName.toLowerCase().split('.').pop() || '';
    // A formatted checkpoint is streamed by its own reader, not read as text:
    // it is what keeps a VMD script naming a large `.fchk` from allocating it
    // twice before the parser sees it.
    if (ext === 'fchk' || ext === 'fch') {
        return { data: await loadFchkStructure(uri, fileName), frames: [] };
    }
    const content = await vscode.workspace.fs.readFile(uri);
    const binaryFormat = detectBinaryTrajectory(content);
    const binaryUnsupported = !binaryFormat &&
        (isTrajectoryFileName(fileName) || ext === 'trj') &&
        looksLikeBinaryTrajectory(content);
    if (binaryUnsupported) {
        return {
            data: {
                atoms: [], bonds: [], title: fileName, hasExplicitBonds: false,
                warnings: ['this trajectory file is binary but is not an XTC, TRR or DCD file: ' +
                    'it cannot be read as coordinates (export it as DCD/XTC/TRR or as a text trajectory)']
            },
            frames: [],
            error: 'unsupported binary trajectory'
        };
    }
    if (binaryFormat) {
        const loaded = await loadBinaryTrajectoryDocument(uri, fileName, content, binaryFormat);
        return { data: loaded.data, frames: [] };
    }
    const text = decodeText(content);
    const setup = await buildParseOptions(uri, fileName, content, text);
    const options = setup && setup.options;

    // A `.log`/`.out` is a Gaussian/ORCA/xTB trajectory only when a program
    // this reader knows has not written it: QE, CP2K, PySCF, ABINIT, VASP
    // (xml/OUTCAR), LAMMPS and ABACUS all use the same extensions. The
    // structured result is kept and used below instead of being probed here
    // and parsed again by `parseFile`.
    let structuredData: MolecularData | undefined;
    if (ext === 'log' || ext === 'out') {
        structuredData = parseStructuredProgramOutput(text, options) || undefined;
    }
    if ((ext === 'log' || ext === 'out') && !structuredData) {
        const logResult = parseLogFile(text, fileName);
        if (logResult.frames.length === 0) {
            return { data: { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false,
                             warnings: ['no Gaussian, ORCA or xTB structure frame was found in this log: ' +
                                        'open the output file that holds the geometry'] },
                     frames: [], error: 'no structure found' };
        }
        const f = logResult.frames[0];
        return {
            data: ensureBonds({
                atoms: f.atoms, bonds: f.bonds, title: f.title,
                hasExplicitBonds: f.hasExplicitBonds, charge: f.charge, multiplicity: f.multiplicity
            }),
            frames: logResult.frames,
            optSteps: logResult.optSteps,
            normalModes: logResult.normalModes,
            routes: logResult.routes,
            optThresholds: logResult.optThresholds,
            charges: logResult.charges,
            summary: logResult.summary,
            // A relaxed scan or an IRC path is a sequence of structures the
            // file states, not a series of optimizer steps: the page draws
            // them in the convergence card's place, so they travel with the
            // rest of what a log states.
            pathways: logResult.pathways,
            logSource: logResult.source
        };
    }

    const parsedText = setup && setup.text !== undefined ? setup.text : text;
    const parsed = structuredData || parseFile(parsedText, fileName, options);
    let data: MolecularData;
    if (parsed.crystal) {
        parsed.bonds = [];
        parsed.crystal.baseBonds = [];
        data = parsed;
    } else {
        data = ensureBonds(parsed);
    }
    applySetupWarnings(data, setup);
    data.filePath = setup && setup.reference ? setup.reference.uri.fsPath : uri.fsPath;
    return { data: data, frames: [] };
}
