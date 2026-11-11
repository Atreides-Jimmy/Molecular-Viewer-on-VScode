// Load and parse one structure file the way the VMD scene and the Diff/Import
// paths consume it.
import * as vscode from 'vscode';
import { MolecularData } from '../types';
import { detectBinaryTrajectory, parseFile, parseLogFile, isStructuredProgramOutput } from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { LoadedStructure } from './payloads';
import { buildParseOptions, applySetupWarnings } from './parseSetup';
import { loadBinaryTrajectoryDocument } from './trajectoryDocuments';
import { decodeText } from './loading';

/** Load and parse one structure file (any supported format). */
export async function loadStructureFile(uri: vscode.Uri): Promise<LoadedStructure> {
    const content = await vscode.workspace.fs.readFile(uri);
    const fileName = uri.path.split('/').pop() || 'unknown.xyz';
    const binaryFormat = detectBinaryTrajectory(content);
    if (binaryFormat) {
        const loaded = await loadBinaryTrajectoryDocument(uri, fileName, content, binaryFormat);
        return { data: loaded.data, frames: [] };
    }
    const text = decodeText(content);
    const ext = fileName.toLowerCase().split('.').pop() || '';
    const setup = await buildParseOptions(uri, fileName, content, text);
    const options = setup && setup.options;

    // A `.log`/`.out` is a Gaussian/ORCA/xTB trajectory only when a program
    // this reader knows has not written it: QE, CP2K, PySCF, ABINIT, VASP
    // (xml/OUTCAR), LAMMPS and ABACUS all use the same extensions.
    if ((ext === 'log' || ext === 'out') && !isStructuredProgramOutput(text)) {
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
            logSource: logResult.source
        };
    }

    const parsedText = setup && setup.text !== undefined ? setup.text : text;
    const parsed = parseFile(parsedText, fileName, options);
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
