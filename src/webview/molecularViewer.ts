// The custom editor provider: open a document, load the structures a VMD scene
// names, and route the page's messages. The page itself is built in page.ts.
import * as vscode from 'vscode';
import * as path from 'path';
import * as fs from 'fs';
import {
    parseFile, parseLogFile, parseXyzTrajectory, LogFrame, OrcaFrame, RouteSection,
    LogChargeSet, LogSummary, parseTcl, hasScalarChannel,
    isFrameLogOutput, isExtXyz, isLammpsDump, parseLammpsDumpTrajectory, parseFchk,
    readFchkColumn, FchkOrbitals, detectBinaryTrajectory, isTrajectoryFileName,
    looksLikeBinaryTrajectory, parseStructuredProgramOutput
} from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { addHydrogens } from '../tools/addHydrogens';
import { beginOptimization, MmHistoryPoint, MmOptimizeResult, MmSession } from '../tools/mmOptimizer';
import { smilesTo3dSafe, validateSmiles } from '../tools/smilesTo3d';
import { CavitySurfaceTask, GridQuality, OrbitalSurfaceTask, SurfaceTask } from '../tools/orbitalMesh';
import { MolecularData, AtomGroup, OptStep, NormalMode, OptThresholds, VmdScene, OrbitalLevels, SolventCavity, LogPathways } from '../types';
import {
    LoadProgress, decodeText, readFileWithProgress, LOAD_PROGRESS_MIN_BYTES,
    VIEW_PROGRESS_MIN_ATOMS, fileSize, wantsLoadProgress, yieldToHost, progressRange,
    withLoadProgress, withToolProgress
} from './loading';
import {
    ATOM_COLORS, LoadedStructure, TrajectoryFramePayload, VmdPayload
} from './payloads';
import { buildParseOptions, applySetupWarnings, readToolFile } from './parseSetup';
import {
    loadBinaryTrajectoryDocument, loadDumpTrajectoryDocument, loadStructureStepsDocument,
    parseTextTrajectoryFrames, pickTextTrajectoryFrame, latticeCellParameters, PickedTextFrame
} from './trajectoryDocuments';
import { loadStructureFile, loadFchkStructure } from './structureLoader';
import { buildVmdMoleculePayload, resolveVmdSource, readEggMetadata, EggMeta } from './vmdPayload';
import { askImportChoice, postImportCopies } from './importHelpers';
import { renderViewerPage } from './page';

/**
 * One Molecular Mechanics run in flight for a webview panel.
 *
 * The optimiser is a session that can be stepped, so the extension host runs it
 * in small chunks, yielding to the message loop between them: every progress
 * chunk is posted to the panel (which draws it live) and — because the host's
 * message loop is only blocked for the length of one chunk — Pause / Resume /
 * Cancel and a new run are serviced between chunks instead of after the whole
 * optimisation. `lastStep` is the step of the newest history point the panel
 * has been sent, so a chunk only carries what the panel has not seen (the
 * module downsamples a very long history, so points are matched by step, not by
 * index); `lastPost` throttles the stream to a readable frame rate.
 */
interface MmRunState {
    id: string;
    session: MmSession;
    paused: boolean;
    stopped: boolean;
    /**
     * Whether a pump wake-up is already queued. A run is stepped by exactly one
     * timer chain: a Resume that arrives while this run is already running (or
     * while a queued wake-up has not fired yet) must not add a second chain, or
     * the same session is stepped twice per interval.
     */
    scheduled: boolean;
    lastStep: number;
    lastPost: number;
}

/** Search passes per chunk: short enough that a pause is felt immediately. */
const MM_CHUNK_STEPS = 10;
/**
 * The chunk is not a fixed number of passes but a fixed *amount of work*: the
 * whole point of stepping the session is that the message loop gets a turn
 * between chunks, so a Pause is felt within one chunk — and one pass costs more
 * the more atoms the selection has. Measured on this machine, ten passes cost
 * ~3 ms on a 72-atom structure but ~42 ms (worst chunk 66 ms) at 356 atoms, so
 * at the 3000-atom cap a fixed ten would hold the loop for a large fraction of a
 * second and a Pause would land that much late. The count therefore starts at
 * ten and is adapted from how long the last chunk actually took; it is only ever
 * raised back towards ten, so a cheap run is not sliced more finely than it
 * needs to be. The answer does not depend on it — a run stepped in chunks of 1,
 * 2, 3, 7, 10 or 1000 passes was verified to be byte-identical.
 */
const MM_CHUNK_TARGET_MS = 60;
/** A chunk this much over target asks for fewer passes; this much under, more. */
const MM_CHUNK_OVER = 1.5;
const MM_CHUNK_UNDER = 0.5;
/** Shortest gap between two progress messages (ms), unless the run just ended. */
const MM_PROGRESS_MS = 75;

/**
 * Other suffixes a Save As filter offers for the same content, so a name that
 * is legitimately `.poscar` or `.dat` is not questioned for a `.vasp`/`.xy`
 * payload. The first entry is the canonical extension the page sends.
 */
const SAVE_EXTENSION_ALIASES: { [extension: string]: string[] } = {
    '.gjf': ['.gjf', '.com', '.gjf03', '.gjf09', '.gjf16'],
    '.xyz': ['.xyz', '.extxyz'],
    '.vasp': ['.vasp', '.poscar', '.contcar'],
    '.lmp': ['.lmp', '.data'],
    '.pwi': ['.pwi', '.in'],
    '.xy': ['.xy', '.dat']
};

/**
 * Conventional extension-less file names a format is read from; a VASP cell
 * named `POSCAR` must stay `POSCAR`, not become `POSCAR.vasp`.
 */
const SAVE_EXTENSIONLESS_NAMES: { [extension: string]: string[] } = {
    '.vasp': ['poscar', 'contcar']
};

/**
 * The bridge between a document and the panel showing it.
 *
 * A custom editor's save, revert and backup methods are handed a document and
 * nothing else, while the content being edited lives in the page: these are the
 * few calls the provider needs from the panel that owns one.
 */
interface EditorPanel {
    document: MolecularDocument;
    post(message: { [key: string]: unknown }): void;
    /**
     * Save what the page holds. `dialog` opens the page's own format dialog
     * (or re-writes the format and path it last exported, see the page's
     * `hostSave`); `path` builds content for the extension the editor's own
     * dialog chose and writes it there. Resolves the written path, or
     * undefined when the user backed out (which keeps the document dirty).
     */
    save(request: { mode: 'dialog' } | { mode: 'path'; target: vscode.Uri; format: string }): Promise<string | undefined>;
    /** The page's hot-exit state (see `backupCustomDocument`). */
    session(): Promise<unknown>;
    /** Rebuild the page from the document's current content. */
    reload(): Promise<void>;
    dispose(): void;
}

/** The extensions a system Save As can carry, mapped to the content builder
 *  that writes them. An extension this table does not know falls back to the
 *  viewer's own default (a plain `.xyz`), which every reader can open. */
const SAVE_AS_FORMATS: { [extension: string]: string } = {
    '.xyz': 'xyz',
    '.extxyz': 'extxyz',
    '.gjf': 'gjf',
    '.com': 'gjf',
    '.gjf03': 'gjf',
    '.gjf09': 'gjf',
    '.gjf16': 'gjf',
    '.cif': 'cif',
    '.vasp': 'vasp',
    '.poscar': 'vasp',
    '.contcar': 'vasp',
    '.cube': 'cube',
    '.vesta': 'vesta',
    '.coord': 'coord',
    '.mol2': 'mol2',
    '.mol': 'mol',
    '.pdb': 'pdb',
    '.ent': 'pdb',
    '.mop': 'mop',
    '.abi': 'abi',
    '.gro': 'gro',
    '.lmp': 'lmp',
    '.data': 'lmp',
    '.pwi': 'qe',
    '.in': 'qe',
    '.inp': 'orca'
};

export class MolecularViewerProvider implements vscode.CustomEditorProvider<MolecularDocument> {
    /** The panel showing each open document (one, see the registration). */
    private readonly panels = new Map<MolecularDocument, EditorPanel>();
    private readonly changeEmitter = new vscode.EventEmitter<vscode.CustomDocumentEditEvent<MolecularDocument>>();
    /**
     * Fired once per edit the page makes, so the editor knows the document is
     * modified, offers the close prompt, and routes its own Undo/Redo commands
     * back into the page's history. The page keeps the undo stack; the two
     * callbacks are how the editor's commands reach it.
     */
    readonly onDidChangeCustomDocument = this.changeEmitter.event;

    constructor(private readonly context: vscode.ExtensionContext) {}

    dispose(): void {
        this.changeEmitter.dispose();
        for (const panel of this.panels.values()) panel.dispose();
        this.panels.clear();
    }

    async openCustomDocument(
        uri: vscode.Uri,
        openContext: vscode.CustomDocumentOpenContext,
        token: vscode.CancellationToken
    ): Promise<MolecularDocument> {
        const document = await this.openDocument(uri, token);
        // A window that was closed with unsaved edits hands back the backup the
        // page wrote (see `backupCustomDocument`). It is applied by the page on
        // top of the freshly parsed file, so a backup that cannot be read only
        // costs the edits — never the document.
        if (openContext.backupId) {
            const session = await this.readBackup(openContext.backupId);
            if (session) {
                document.restoreSession = session;
            } else {
                document.data.warnings = (document.data.warnings || []).concat(
                    'the unsaved edits of this file could not be restored from the editor backup: the file itself is shown');
            }
        }
        return document;
    }

    /** Read the JSON session a previous window wrote, or undefined. */
    private async readBackup(backupId: string): Promise<unknown> {
        try {
            const text = decodeText(await vscode.workspace.fs.readFile(vscode.Uri.parse(backupId)));
            const parsed = JSON.parse(text);
            // `v` names the shape; a session of another version is not guessed
            // at, because half-applied structure state is worse than none.
            return parsed && parsed.v === 1 && Array.isArray(parsed.atoms) ? parsed : undefined;
        } catch {
            return undefined;
        }
    }

    /**
     * Read and parse one document, dispatching on its name.
     *
     * Everything below this branch works on a file already in memory (or on a
     * checkpoint streamed by its own reader), so the dispatch has to happen
     * before either.
     */
    private async openDocument(uri: vscode.Uri, token: vscode.CancellationToken): Promise<MolecularDocument> {
        if (token.isCancellationRequested) throw new vscode.CancellationError();
        const fileName = uri.path.split('/').pop() || 'unknown.xyz';
        const ext0 = fileName.toLowerCase().split('.').pop() || '';
        // A formatted checkpoint is the one file this reader streams instead of
        // reading: a 432-atom job's is 850 MB, and turning that into a string
        // (then into lines) would cost two gigabytes of heap before any parser
        // ran. The branch is here, above `readFile`, because everything below
        // assumes the file is already in memory.
        if (ext0 === 'fchk' || ext0 === 'fch') {
            return this.openFchkDocument(uri, fileName, token);
        }
        const size = await fileSize(uri);
        if (!wantsLoadProgress(size, ext0)) {
            return this.readDocument(uri, fileName, size, token);
        }
        return withLoadProgress('Molecular Viewer: opening ' + fileName, token, (progress, isCancelled) =>
            this.readDocument(uri, fileName, size, token, progress, isCancelled));
    }

    /**
     * The editor's Save (Ctrl+S, Save All, the close prompt's Save) and Save As.
     *
     * Both write through the page, because the page is where the structure is:
     * Save re-uses the last format the user exported (or asks for one), and
     * Save As honours the extension of the file the editor's own dialog chose.
     * A cancelled dialog rejects with `CancellationError`, which is how the
     * editor is told the save did not happen and the document stays dirty.
     */
    async saveCustomDocument(document: MolecularDocument, _cancellation: vscode.CancellationToken): Promise<void> {
        const panel = this.panels.get(document);
        if (!panel) throw new Error('the viewer is not showing this document, so there is nothing to save from it');
        const written = await panel.save({ mode: 'dialog' });
        if (!written) throw new vscode.CancellationError();
    }

    async saveCustomDocumentAs(document: MolecularDocument, destination: vscode.Uri, _cancellation: vscode.CancellationToken): Promise<void> {
        const panel = this.panels.get(document);
        if (!panel) throw new Error('the viewer is not showing this document, so there is nothing to save from it');
        const name = destination.path.split('/').pop() || '';
        const dot = name.lastIndexOf('.');
        const ext = dot > 0 ? name.substring(dot).toLowerCase() : '';
        const format = SAVE_AS_FORMATS[ext] || 'xyz';
        const written = await panel.save({ mode: 'path', target: destination, format: format });
        if (!written) throw new vscode.CancellationError();
    }

    async revertCustomDocument(document: MolecularDocument, token: vscode.CancellationToken): Promise<void> {
        const fresh = await this.openDocument(document.uri, token);
        document.setContent(fresh);
        const panel = this.panels.get(document);
        if (panel) await panel.reload();
    }

    /**
     * Hot exit: keep what the page holds when the window is closed with unsaved
     * edits.
     *
     * The session is the page's own state — the structure on screen with its
     * electronic state, and the frames that were edited — not a file in any
     * format, so nothing is lost to a round trip through XYZ. The page writes
     * it, the host only stores it beside VS Code's own backups; a window that
     * reopens the file reads it back (see `openCustomDocument`).
     *
     * A document with no panel (a tab whose webview is gone) has no state to
     * back up, and says so rather than pretending the edits are safe.
     */
    async backupCustomDocument(
        document: MolecularDocument,
        context: vscode.CustomDocumentBackupContext,
        _cancellation: vscode.CancellationToken
    ): Promise<vscode.CustomDocumentBackup> {
        const panel = this.panels.get(document);
        const session = panel ? await panel.session() : undefined;
        if (!session) throw new Error('the viewer is not showing this document, so its edits cannot be backed up');
        const target = context.destination;
        if (target.scheme === 'file') {
            // The parent of the suggested backup path is not guaranteed to
            // exist (the API says so); creating it costs one call.
            await vscode.workspace.fs.createDirectory(vscode.Uri.file(path.dirname(target.fsPath)));
        }
        await vscode.workspace.fs.writeFile(target, new TextEncoder().encode(JSON.stringify(session)));
        return {
            id: target.toString(),
            delete: () => { void vscode.workspace.fs.delete(target, { useTrash: false }).then(undefined, () => undefined); }
        };
    }

    /**
     * Read and parse one document, optionally behind a progress notification.
     *
     * Everything below the checkpoint branch lives here so the plain path and
     * the reported path run the same code: a small file takes the first with no
     * notification, and a large one gets the reading percentage and the phase
     * it is in without a second copy of the parse.
     */
    private async readDocument(
        uri: vscode.Uri,
        fileName: string,
        size: number,
        token: vscode.CancellationToken,
        progress?: LoadProgress,
        isCancelled?: () => boolean
    ): Promise<MolecularDocument> {
        const wasCancelled = (): boolean => token.isCancellationRequested || (!!isCancelled && isCancelled());
        // Reading owns the first 30% of the bar (in whole points, so a bar
        // really fills), decoding/preparing the next 5, and the parser the rest.
        const readRange = progress ? progressRange(progress, 0, 30) : undefined;
        const content = readRange
            ? await readFileWithProgress(uri, size, readRange, wasCancelled)
            : await vscode.workspace.fs.readFile(uri);
        // The tab may have been closed while the file was being read; the
        // parse below can be expensive and its result has nowhere to go.
        if (wasCancelled()) throw new vscode.CancellationError();

        const ext = fileName.toLowerCase().split('.').pop() || '';
        // A binary trajectory is detected by its magic before anything decodes
        // it as text; the companion scripts and topology beside it are read by
        // the trajectory loader itself. Decoding first would also copy a
        // multi-megabyte binary file through a string only to discard it.
        const binaryFormat = detectBinaryTrajectory(content);
        // `.trj` is text for xTB and binary for CHARMM/NAMD. A binary file of a
        // format this reader does not know (no XTC/TRR/DCD magic) must not be
        // decoded into mojibake and reported as a parse failure: it gets a
        // warning that names the real reason instead.
        const binaryUnsupported = !binaryFormat &&
            (isTrajectoryFileName(fileName) || ext === 'trj') &&
            looksLikeBinaryTrajectory(content);
        const textContent = binaryFormat || binaryUnsupported ? '' : decodeText(content);
        if (progress) {
            progress.report({
                increment: 5,
                message: binaryFormat ? 'Reading the trajectory'
                    : binaryUnsupported ? 'Checking the trajectory'
                    : 'Preparing the parse'
            });
        }
        const setup = binaryFormat || binaryUnsupported ? undefined : await buildParseOptions(uri, fileName, content, textContent);
        const parseOptions = setup && setup.options;
        const parseText = setup && setup.text !== undefined ? setup.text : textContent;
        if (wasCancelled()) throw new vscode.CancellationError();
        const parseRange = progress ? progressRange(progress, 35, 100) : undefined;
        if (parseRange) {
            // The parse is a synchronous pass, so the phase is named and the
            // message loop gets one turn to paint it before the host is busy.
            parseRange(0, binaryFormat ? 'Reading the trajectory'
                : ext === 'tcl' ? 'Loading the script'
                : (ext === 'log' || ext === 'out') ? 'Reading the output'
                : 'Parsing the structure');
            await yieldToHost();
            // Cancel may have been pressed during that turn; do not start a
            // synchronous parse whose answer nobody wants.
            if (wasCancelled()) throw new vscode.CancellationError();
        }
        let data: MolecularData;
        let frames: (LogFrame | OrcaFrame)[] = [];
        let trajectoryFrames: TrajectoryFramePayload[] | undefined;
        let trajectoryDocument = false;
        let atomGroups: AtomGroup[] | undefined;
        let optSteps: OptStep[] | undefined;
        let normalModes: NormalMode[] | undefined;
        let routes: RouteSection[] | undefined;
        let optThresholds: OptThresholds | undefined;
        let logSource: string | undefined;
        let charges: LogChargeSet[] | undefined;
        let summary: LogSummary[] | undefined;
        // A log states eigenvalues and the cavity model, never the coefficients:
        // the level diagram can be drawn from these two, a surface cannot.
        let orbitalLevels: OrbitalLevels | undefined;
        let cavity: SolventCavity | undefined;
        // A relaxed scan and/or an IRC path: structures along one coordinate,
        // drawn in the convergence card's place.
        let pathways: LogPathways | undefined;
        let vmdPayload: VmdPayload | undefined;
        let autoBeta = false;

        try {
        // Structured-output detection is a full pass per detector, so it owns
        // the first half of the parser's slice; a frame log is recognised by
        // cheap markers and leaves the whole slice to its parser. The detector
        // returns the parsed structure, not only a boolean: `parseFile` would
        // otherwise run the same full passes again to find the reader.
        let logStructured = false;
        let logDetected = false;
        let structuredData: MolecularData | undefined;
        if (ext === 'log' || ext === 'out') {
            if (isFrameLogOutput(textContent)) {
                logStructured = false;
            } else {
                logDetected = true;
                structuredData = parseStructuredProgramOutput(textContent, parseOptions,
                    parseRange ? (fraction) => parseRange(fraction * 0.5) : undefined) || undefined;
                logStructured = !!structuredData;
            }
        }
        if (binaryUnsupported) {
            data = {
                atoms: [], bonds: [], title: fileName, hasExplicitBonds: false,
                warnings: ['this trajectory file is binary but is not an XTC, TRR or DCD file: ' +
                    'it cannot be read as coordinates. Export it as DCD/XTC/TRR, or as a text ' +
                    'trajectory (e.g. .xyz) that this viewer can open']
            };
        } else if (binaryFormat) {
            const loaded = await loadBinaryTrajectoryDocument(uri, fileName, content, binaryFormat);
            data = loaded.data;
            trajectoryFrames = loaded.trajectoryFrames;
            trajectoryDocument = true;
        } else if (ext === 'tcl') {
            // A script is a display recipe: every `mol new` it performs becomes
            // one selectable molecule in the viewer, and the first one that
            // loads is also the document's structure so Save As, Diff and the
            // editing tools keep working on a real coordinate file.
            const tclResult = parseTcl(textContent);
            // The script is a few kilobytes; the structures it names are the
            // work. One large source (or several that add up) gets the same
            // notification an opened structure does, while a small scene still
            // opens without a notification flashing up.
            let sceneBytes = 0;
            for (const molecule of tclResult.scene.molecules) {
                const sourceUri = resolveVmdSource(uri, molecule.sourceFile);
                if (sourceUri) sceneBytes += await fileSize(sourceUri);
            }
            const loaded = parseRange
                ? await this.loadTclScene(tclResult.scene, uri, parseRange)
                : sceneBytes >= LOAD_PROGRESS_MIN_BYTES
                    ? await withLoadProgress('Molecular Viewer: opening ' + fileName, token, (sceneProgress) =>
                        this.loadTclScene(tclResult.scene, uri, progressRange(sceneProgress, 0, 100)))
                    : await this.loadTclScene(tclResult.scene, uri);
            vmdPayload = loaded.payload && loaded.payload.molecules.length > 0 ? loaded.payload : undefined;
            data = loaded.data;
            frames = loaded.frames;
            optSteps = loaded.optSteps;
            normalModes = loaded.normalModes;
            routes = loaded.routes;
            logSource = loaded.logSource;
            // The first molecule of a scene can be a program output, and a
            // scene is otherwise the only way that file is read. Everything
            // the file states about it travels, not only its geometry: the
            // Summary and Charges cards, the convergence thresholds and the
            // scan/IRC card all read the fields below, and a script naming a
            // Gaussian log used to open with none of them.
            optThresholds = loaded.optThresholds;
            charges = loaded.charges;
            summary = loaded.summary;
            pathways = loaded.pathways;
            if (!vmdPayload) {
                // Legacy path only: with a resolved scene the payload carries the
                // colours, and it is the one that decides them.
                atomGroups = tclResult.groups.length > 0 ? tclResult.groups : undefined;
                data.atomGroups = atomGroups;
            }
        } else if ((ext === 'log' || ext === 'out') && !logStructured) {
            const logBase = logDetected ? 0.5 : 0;
            const logResult = parseLogFile(textContent, fileName,
                parseRange ? (fraction) => parseRange(logBase + fraction * (1 - logBase),
                    'Reading the output · ' + Math.round(fraction * 100) + '%') : undefined);
            frames = logResult.frames;
            optSteps = logResult.optSteps;
            normalModes = logResult.normalModes;
            routes = logResult.routes;
            optThresholds = logResult.optThresholds;
            logSource = logResult.source;
            charges = logResult.charges;
            summary = logResult.summary;
            orbitalLevels = logResult.orbitals;
            cavity = logResult.cavity;
            pathways = logResult.pathways;
            if (frames.length > 0) {
                data = ensureBonds({
                    atoms: frames[0].atoms,
                    bonds: frames[0].bonds,
                    title: frames[0].title,
                    hasExplicitBonds: frames[0].hasExplicitBonds,
                    charge: frames[0].charge,
                    multiplicity: frames[0].multiplicity
                });
            } else {
                data = { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false,
                         warnings: ['no Gaussian, ORCA or xTB structure frame was found in this log: ' +
                                    'open the output file that holds the geometry'] };
            }
        } else if (isLammpsDump(textContent)) {
            // A dump with more than one frame is a trajectory, not a crystal:
            // its frames can change cell and the molecules cross the box, so it
            // is drawn outside the symmetry machinery. A single-frame dump keeps
            // the periodic-structure path it has always had. The decision uses
            // the file's frame count, not the sampled list: a long dump sampled
            // down to one frame is still a trajectory (and says so in its
            // warnings), not the first frame of a crystal.
            const dumpOptions = parseOptions && parseOptions.lammps ? parseOptions.lammps : undefined;
            const dumpTrajectory = parseLammpsDumpTrajectory(textContent, dumpOptions);
            if (dumpTrajectory.frameCount > 1) {
                const loaded = await loadDumpTrajectoryDocument(uri, fileName, dumpTrajectory, dumpOptions);
                data = loaded.data;
                trajectoryFrames = loaded.trajectoryFrames;
                trajectoryDocument = true;
            } else {
                data = parseFile(parseText, fileName, parseOptions);
                if (data.crystal) {
                    data.bonds = [];
                    data.crystal.baseBonds = [];
                } else {
                    data = ensureBonds(data);
                }
            }
        } else if ((ext === 'xyz' || ext === 'trj') && !isExtXyz(textContent)) {
            // A multi-frame XYZ is a trajectory — ORCA writes one for every
            // optimization as `<job>_trj.xyz`, with each frame's energy in its
            // comment line — so it gets the same frame navigation (and, when the
            // comment states an energy per frame, the same convergence chart) a
            // log does. xTB's `xtb.trj` uses the same layout. A single-frame
            // file takes the ordinary structure path. GPUMD/ASE extended XYZ
            // keeps its cell on the comment line, which the trajectory reader
            // cannot carry, so it takes the structure path.
            const traj = parseXyzTrajectory(textContent);
            if (traj.frames.length > 1) {
                frames = traj.frames.map((f, k) => ({
                    atoms: f.atoms,
                    bonds: [],
                    title: f.title,
                    hasExplicitBonds: false,
                    stepLabel: 'Frame ' + (k + 1),
                    energy: f.energy
                }));
                // The energy curve needs a value for every frame: a partial
                // series would put step k at frame k of a different numbering.
                const steps: OptStep[] = [];
                for (let k = 0; k < traj.frames.length; k++) {
                    const e = traj.frames[k].energy;
                    if (e === undefined) { steps.length = 0; break; }
                    steps.push({ step: k + 1, energy: e });
                }
                optSteps = steps.length > 1 ? steps : undefined;
                data = ensureBonds({
                    atoms: frames[0].atoms,
                    bonds: [],
                    title: frames[0].title,
                    hasExplicitBonds: false
                });
            } else {
                data = ensureBonds(parseFile(parseText, fileName, parseOptions));
            }
        } else if (structuredData) {
            // The structured reader's own result, already parsed during the
            // single detection pass above. Re-parsing through `parseFile`
            // would scan the text for the same detectors a second time.
            data = structuredData;
            if (data.crystal) {
                data.bonds = [];
                data.crystal.baseBonds = [];
            } else {
                data = ensureBonds(data);
            }
        } else {
            data = parseFile(parseText, fileName, parseOptions);
            if (data.crystal) {
                data.bonds = [];
                data.crystal.baseBonds = [];
            } else {
                data = ensureBonds(data);
            }
            // Some visualizations ship as a bare PDB instead of a script: the
            // per-atom scalar itself *is* the visualization, so it is painted
            // straight away (the toolbar toggle turns it off again).
            if (ext === 'pdb' || ext === 'ent') {
                autoBeta = hasScalarChannel(textContent, data.atoms);
            }
        }
        // A reader that states several geometries (a CP2K molecular-dynamics
        // output, an ABACUS `MD_dump`, a Quantum ESPRESSO relax or MD) describes
        // one system over time: the document gets the frame navigation a text
        // dump already has, and each frame brings the cell it was computed in. A
        // file whose reader states one geometry keeps the crystal path it had.
        const statedSteps = data.steps;
        if (!trajectoryDocument && statedSteps && statedSteps.length > 1 && data.atoms.length > 0) {
            const document = loadStructureStepsDocument(
                fileName, data, statedSteps, data.stepsInFile || statedSteps.length
            );
            data = document.data;
            trajectoryFrames = document.trajectoryFrames;
            trajectoryDocument = true;
        }
        } catch (err) {
            // A Cancel on a scene load is not a parse error: it must abort the
            // open rather than leave a tab showing "the file could not be
            // parsed".
            if (err instanceof vscode.CancellationError) throw err;
            const msg = err instanceof Error ? err.message : String(err);
            data = { atoms: [], bonds: [], title: 'Parse error: ' + msg, hasExplicitBonds: false,
                     warnings: ['the file could not be parsed: ' + msg] };
            frames = [];
            optSteps = undefined;
            normalModes = undefined;
            optThresholds = undefined;
            charges = undefined;
            summary = undefined;
            orbitalLevels = undefined;
            cavity = undefined;
            pathways = undefined;
        }

        // The parse is synchronous, so a Cancel pressed while it ran could only
        // be seen here: the whole result is dropped rather than written into a
        // document whose tab is gone.
        if (wasCancelled()) throw new vscode.CancellationError();
        applySetupWarnings(data, setup);
        if (trajectoryDocument) {
            // A trajectory is its own file: a read_data an accompanying script
            // names is a topology, and Save As must write beside the trajectory.
            data.filePath = uri.fsPath;
        } else if (ext !== 'tcl') {
            data.filePath = setup && setup.reference ? setup.reference.uri.fsPath : uri.fsPath;
        }
        // Formats whose parser has no line callback still finish the bar.
        parseRange?.(1, 'Preparing the document');

        return new MolecularDocument(uri, data, frames, optSteps, normalModes, logSource, routes, vmdPayload, autoBeta, optThresholds, trajectoryFrames, charges, summary, undefined, cavity, orbitalLevels, pathways);
    }

    /**
     * Read a Gaussian formatted checkpoint.
     *
     * The file is streamed, so the cost is one pass over its bytes and whatever
     * the small arrays need; the MO coefficient block is skipped by line count
     * and only its byte offset is kept, which is what lets an orbital be read
     * back one column at a time (see `readFchkColumn`).
     *
     * A checkpoint is local to the machine that wrote it: `fs.createReadStream`
     * needs a real path, and a `vscode.workspace.fs` stream does not exist for
     * a virtual scheme. Anything that is not a `file:` URI is read through the
     * workspace API instead (a remote workspace buffers it, which is slower but
     * correct) only when it is small enough to be safe.
     */
    private async openFchkDocument(uri: vscode.Uri, fileName: string, token: vscode.CancellationToken): Promise<MolecularDocument> {
        let size = 0;
        try {
            const stat = await vscode.workspace.fs.stat(uri);
            size = stat.size || 0;
        } catch {
            size = 0;
        }
        // A checkpoint that is not a local file can only be read through the
        // workspace API, which returns the whole file as one buffer. That is how
        // a 900 MB checkpoint would become a 900 MB allocation, so the fallback
        // is offered only up to a size a buffer can take; beyond it the tab says
        // why instead of failing with an allocation error.
        const FCHK_BUFFER_LIMIT = 256 * 1024 * 1024;
        if (uri.scheme !== 'file' && size > FCHK_BUFFER_LIMIT) {
            return new MolecularDocument(uri, {
                atoms: [], bonds: [], title: fileName, hasExplicitBonds: false,
                filePath: uri.fsPath,
                warnings: ['this formatted checkpoint is ' + Math.round(size / (1024 * 1024)) +
                    ' MB and is not a local file: it can only be read through the workspace file API, ' +
                    'which would hold all of it in memory at once. Open it from a local folder.']
            });
        }

        let result: Awaited<ReturnType<typeof parseFchk>>;
        try {
            result = await vscode.window.withProgress({
                location: vscode.ProgressLocation.Notification,
                title: 'Molecular Viewer: reading ' + fileName,
                cancellable: true
            }, async (progress, progressToken) => {
                progress.report({ message: 'reading the checkpoint' });
                let reported = 0;
                const onProgress = (read: number) => {
                    if (size > 0) {
                        const pct = Math.min(99, Math.floor((read / size) * 100));
                        if (pct > reported) {
                            progress.report({ increment: pct - reported, message: pct + '%' });
                            reported = pct;
                        }
                    }
                };
                if (uri.scheme === 'file') {
                    const stream = fs.createReadStream(uri.fsPath, { highWaterMark: 1 << 16 });
                    const cancelled = new Promise<never>((_, reject) => {
                        progressToken.onCancellationRequested(() => {
                            stream.destroy();
                            reject(new vscode.CancellationError());
                        });
                    });
                    return await Promise.race([
                        parseFchk(stream, { size, path: uri.fsPath, onProgress }),
                        cancelled
                    ]);
                }
                const bytes = await vscode.workspace.fs.readFile(uri);
                return await parseFchk((async function* () { yield bytes; })(), { size: bytes.length, path: '' });
            });
        } catch (e: any) {
            // A truncated or non-checkpoint file must open the viewer with the
            // reader's own words on the canvas, the way every other parse
            // failure does, instead of failing the tab with "Unknown error".
            if (e instanceof vscode.CancellationError) throw e;
            return new MolecularDocument(uri, {
                atoms: [], bonds: [], title: fileName, hasExplicitBonds: false,
                filePath: uri.fsPath,
                warnings: ['this formatted checkpoint could not be read: ' + (e && e.message ? String(e.message) : String(e))]
            });
        }

        if (token.isCancellationRequested) throw new vscode.CancellationError();

        // The coefficients are read back one column at a time by byte offset, so
        // only a checkpoint that is a real local file can become a surface. A
        // document read through the workspace API (a remote or virtual scheme)
        // was buffered once for this pass and cannot be addressed again; the
        // levels are still its own, so they stay and only the surface is refused.
        if (result.orbitals && !result.orbitals.path) {
            for (const spin of result.orbitals.spins) spin.layout = undefined;
            result.warnings.push('this checkpoint was read without a local file path, so its orbital coefficients cannot be read back: the orbital levels are shown, isosurfaces need the file on disk');
        }

        const title = /^molecule name$/i.test(result.title.trim()) ? fileName : result.title;
        const data: MolecularData = ensureBonds({
            atoms: result.atoms,
            bonds: [],
            title,
            hasExplicitBonds: false,
            charge: result.charge,
            multiplicity: result.multiplicity,
            warnings: result.warnings.length > 0 ? result.warnings.slice() : undefined
        });
        data.filePath = uri.fsPath;

        // The headline results, in the shape the summary panel already reads. A
        // checkpoint states no point group, no cpu time and no termination, so
        // those lines simply do not appear.
        const summary: LogSummary[] = [];
        if (result.energies.length > 0 || result.rmsGradient !== undefined) {
            summary.push({
                job: 1,
                nAtoms: result.nAtoms,
                nBasis: result.nBasis,
                charge: result.charge,
                multiplicity: result.multiplicity,
                energies: result.energies,
                rmsGradient: result.rmsGradient,
                rmsGradientLabel: result.rmsGradient !== undefined ? 'RMS Force' : undefined
            });
        }

        return new MolecularDocument(uri, data, [], undefined, undefined, undefined, undefined, undefined, false, undefined, undefined,
            result.charges.length > 0 ? result.charges : undefined,
            summary.length > 0 ? summary : undefined,
            result.orbitals, result.cavity, result.orbitals ? orbitalLevelsOf(result.orbitals) : undefined);
    }

    /**
     * Load every structure a VMD script refers to and resolve each one's
     * representations into drawable properties. Molecules that cannot be loaded
     * are skipped with a warning instead of failing the whole script, and the
     * first one that does load becomes the document's structure.
     */
    private async loadTclScene(scene: VmdScene, scriptUri: vscode.Uri, reportFraction?: (fraction: number, message?: string) => void): Promise<{
        payload: VmdPayload;
        data: MolecularData;
        frames: (LogFrame | OrcaFrame)[];
        optSteps?: OptStep[];
        normalModes?: NormalMode[];
        routes?: RouteSection[];
        optThresholds?: OptThresholds;
        charges?: LogChargeSet[];
        summary?: LogSummary[];
        pathways?: LogPathways;
        logSource?: string;
    }> {
        const payload: VmdPayload = {
            molecules: [],
            colorScaleMethod: scene.colorScaleMethod || 'RGB',
            warnings: []
        };

        let first: LoadedStructure | undefined;

        // eggFile candidates are scored and picked once per (directory, structure)
        // pair, so a multi-molecule script sharing one structure asks a single time.
        const eggMetaCache = new Map<string, Promise<EggMeta | undefined>>();

        for (let i = 0; i < scene.molecules.length; i++) {
            const molecule = scene.molecules[i];
            if (reportFraction) {
                reportFraction(i / scene.molecules.length,
                    'Loading molecule ' + (i + 1) + ' of ' + scene.molecules.length +
                    (molecule.sourceFile ? ': ' + molecule.sourceFile : ''));
            }
            const sourceUri = resolveVmdSource(scriptUri, molecule.sourceFile);
            if (!sourceUri) {
                payload.warnings.push('mol new without a structure file was skipped');
                continue;
            }
            let loaded: LoadedStructure;
            try {
                loaded = await loadStructureFile(sourceUri);
            } catch (e) {
                const msg = e instanceof Error ? e.message : String(e);
                payload.warnings.push('could not read ' + molecule.sourceFile + ': ' + msg);
                continue;
            }
            if (loaded.error || loaded.data.atoms.length === 0) {
                payload.warnings.push('no structure found in ' + molecule.sourceFile);
                continue;
            }
            if (!first) first = loaded;
            // 只有画 egg 的脚本才需要去找 eggFile/unitfile。
            const wantsEggs = molecule.reps.some(function (r) { return /^surf/i.test(r.style); });
            const eggMeta = wantsEggs ? await readEggMetadata(sourceUri, eggMetaCache, {
                atoms: loaded.data.atoms,
                bonds: loaded.data.bonds,
                molecule: molecule,
                scene: scene
            }) : undefined;
            payload.molecules.push(buildVmdMoleculePayload(molecule, loaded.data, scene, function (element) {
                return ATOM_COLORS[element] || '#FF1493';
            }, eggMeta));
        }

        if (!first) {
            reportFraction?.(1);
            const reference = scene.molecules.length > 0 ? scene.molecules[0].sourceFile : '';
            const title = reference
                ? 'Failed to load source file: ' + reference
                : 'No source file specified in TCL';
            return { payload: payload, data: { atoms: [], bonds: [], title: title, hasExplicitBonds: false }, frames: [] };
        }

        reportFraction?.(1);
        return {
            payload: payload,
            data: first.data,
            frames: first.frames,
            optSteps: first.optSteps,
            normalModes: first.normalModes,
            routes: first.routes,
            optThresholds: first.optThresholds,
            charges: first.charges,
            summary: first.summary,
            pathways: first.pathways,
            logSource: first.logSource
        };
    }

    async resolveCustomEditor(
        document: MolecularDocument,
        webviewPanel: vscode.WebviewPanel,
        token: vscode.CancellationToken
    ): Promise<void> {
        webviewPanel.webview.options = {
            enableScripts: true,
        };

        if (token.isCancellationRequested) return;
        // The panel outlives this call, and a closed tab must not be posted to
        // (an assignment or a message to a disposed webview throws).
        let disposed = false;
        const post = (message: { [key: string]: unknown }): void => {
            if (disposed) return;
            try {
                webviewPanel.webview.postMessage(message);
            } catch {
                disposed = true;
            }
        };
        const buildHtml = (report?: (fraction: number, message?: string) => void) => renderViewerPage(
            this.context.extensionUri,
            {
                data: document.data,
                frames: document.frames,
                optSteps: document.optSteps,
                normalModes: document.normalModes,
                logSource: document.logSource,
                routes: document.routes,
                vmdPayload: document.vmdPayload,
                autoBeta: document.autoBeta,
                optThresholds: document.optThresholds,
                trajectoryFrames: document.trajectoryFrames,
                charges: document.charges,
                summary: document.summary,
                orbitalLevels: document.orbitalLevels,
                orbitals: document.orbitals,
                cavity: document.cavity,
                pathways: document.pathways,
                restoreSession: document.restoreSession,
                onProgress: report
            }
        );
        // A large trajectory is serialised into the page's payload atom by atom;
        // that can take a moment after the document itself is ready, and the
        // page's own loading line only appears once the HTML is assigned. The
        // same wrapper serves the first build and a reload, so a re-read of a
        // big file reports the same way the first read did.
        const renderPage = async (): Promise<string> => {
            let viewAtoms = document.data.atoms.length;
            for (const frame of document.frames) viewAtoms += frame.atoms.length;
            if (document.trajectoryFrames) {
                for (const frame of document.trajectoryFrames) viewAtoms += Math.floor(frame.xyz.length / 3);
            }
            return viewAtoms >= VIEW_PROGRESS_MIN_ATOMS
                ? await vscode.window.withProgress({
                    location: vscode.ProgressLocation.Notification,
                    title: 'Molecular Viewer: preparing the 3D view',
                    cancellable: false
                }, async (progress) => {
                    const report = progressRange(progress, 0, 100);
                    report(0, 'Building the scene (' + viewAtoms + ' atoms)');
                    await yieldToHost();
                    const page = await buildHtml(report);
                    report(1, 'Ready');
                    return page;
                })
                : await buildHtml();
        };
        const html = await renderPage();
        // Reading the renderer is asynchronous; if the editor was closed in
        // the meantime, assigning html throws "Webview is disposed" and the
        // extension log fills with a failure the user never caused.
        if (token.isCancellationRequested) return;
        try {
            webviewPanel.webview.html = html;
        } catch {
            return;
        }

        // --- The editor's own channels: save, revert, backup, undo/redo ---
        //
        // What is being edited lives in the page, so each of these is a message
        // into it and an answer back. A save the editor asked for waits for that
        // answer: the page's dialog decides whether anything was written, and a
        // save nobody answers must not hang the close that triggered it.
        let saveWaiter: { resolve: (written: string | undefined) => void } | null = null;
        let saveTimer: NodeJS.Timeout | undefined;
        /** Where a host-requested save has to go (the editor's Save As target). */
        let forcedSaveTarget: vscode.Uri | undefined;
        const finishSave = (written: string | undefined): void => {
            if (saveTimer) {
                clearTimeout(saveTimer);
                saveTimer = undefined;
            }
            forcedSaveTarget = undefined;
            const waiter = saveWaiter;
            saveWaiter = null;
            if (waiter) waiter.resolve(written);
        };
        let sessionWaiter: ((session: unknown) => void) | null = null;
        let sessionTimer: NodeJS.Timeout | undefined;
        const finishSession = (session: unknown): void => {
            if (sessionTimer) {
                clearTimeout(sessionTimer);
                sessionTimer = undefined;
            }
            const waiter = sessionWaiter;
            sessionWaiter = null;
            if (waiter) waiter(session);
        };
        /** Paths this panel wrote, and when: a save through the viewer is not
         *  "the file changed on disk" to the watcher below. */
        const selfWrites = new Map<string, number>();
        const reloadPage = async (): Promise<void> => {
            const page = await renderPage();
            if (disposed || token.isCancellationRequested) return;
            try {
                webviewPanel.webview.html = page;
            } catch {
                disposed = true;
            }
        };
        const bridge: EditorPanel = {
            document,
            post,
            save: (request) => new Promise<string | undefined>((resolve) => {
                if (saveWaiter) finishSave(undefined);
                forcedSaveTarget = request.mode === 'path' ? request.target : undefined;
                // The waiter is in place before the page is asked: an answer
                // that arrives at once (a page that has nothing to save) must
                // find someone to resolve.
                saveWaiter = { resolve };
                post(request.mode === 'path'
                    ? { command: 'hostSaveTo', format: request.format }
                    : { command: 'hostSave' });
                // A page that never answers (a webview that failed to start)
                // cancels the save instead of leaving the editor waiting.
                saveTimer = setTimeout(() => finishSave(undefined), 5 * 60 * 1000);
            }),
            session: () => new Promise<unknown>((resolve) => {
                if (sessionWaiter) finishSession(undefined);
                sessionWaiter = resolve;
                post({ command: 'hostSession' });
                sessionTimer = setTimeout(() => finishSession(undefined), 10 * 1000);
            }),
            reload: reloadPage,
            dispose: () => {
                finishSave(undefined);
                finishSession(undefined);
            }
        };
        this.panels.set(document, bridge);

        // --- The file on disk ---
        //
        // The view holds a parsed copy of the file, so another program -- a job
        // still appending to its own log, an editor, a script -- can change the
        // file while the page shows the version it opened. The watcher never
        // reloads on its own: an optimisation log changes every few seconds and
        // a prompt per change would be noise. It tells the page, and the page
        // shows one notice with a Reload button.
        let fileWatcher: vscode.FileSystemWatcher | undefined;
        let noticeTimer: NodeJS.Timeout | undefined;
        if (document.uri.scheme === 'file' && document.uri.fsPath) {
            const watched = document.uri.fsPath;
            const announce = (deleted: boolean): void => {
                if (disposed) return;
                const written = selfWrites.get(watched);
                if (written !== undefined && Date.now() - written < 2000) return;
                if (noticeTimer) clearTimeout(noticeTimer);
                noticeTimer = setTimeout(() => {
                    noticeTimer = undefined;
                    post({ command: 'fileChanged', deleted: deleted });
                }, 300);
            };
            try {
                fileWatcher = vscode.workspace.createFileSystemWatcher(
                    new vscode.RelativePattern(vscode.Uri.file(path.dirname(watched)), path.basename(watched))
                );
                fileWatcher.onDidChange(() => announce(false));
                fileWatcher.onDidCreate(() => announce(false));
                fileWatcher.onDidDelete(() => announce(true));
            } catch {
                // A scheme or a workspace that cannot be watched keeps the
                // behaviour it had: the view shows the copy it opened.
                fileWatcher = undefined;
            }
        }

        // The panel's Molecular Mechanics run (at most one at a time). The pump
        // steps it in chunks and yields between them, so the panel keeps
        // receiving progress and its Pause / Resume / Cancel -- or a brand new
        // run -- are handled while the optimisation is still working.
        let mmRun: MmRunState | null = null;

        // --- Orbitals and the cavity: the surface work runs here ---
        //
        // Sampling a 4680-function orbital on a million points costs about a
        // second of pure arithmetic, and cutting its two lobes out of that grid
        // another fraction of one. That is far too long to do inside a message
        // handler — the panel would freeze with no way to show progress — so the
        // work is a task stepped in slices with a yield between them: the same
        // shape as the Molecular Mechanics pump, and the reason the panel can
        // draw a progress bar at all. One run at a time per panel; a new request
        // cancels the one in flight.
        //
        // The meshes travel as base64 of their typed arrays. A JSON array of
        // 300 000 floats is 2 MB of digits; the same bytes as base64 are 1.6 MB
        // and decode with one pass, which matters because a message is copied
        // through the webview's own serialisation.
        interface SurfaceJob {
            key: string;
            label: string;
            /** Colour pair index, so several orbitals stay tellable apart. */
            hue: number;
            task: SurfaceTask;
        }
        let surfaceRun: { id: number; reqId: number; cancelled: boolean; jobs: SurfaceJob[] } | null = null;
        let surfaceRunId = 0;

        const encodeFloats = (a: Float32Array): string =>
            Buffer.from(a.buffer, a.byteOffset, a.byteLength).toString('base64');
        const encodeIndices = (a: Uint32Array): string =>
            Buffer.from(a.buffer, a.byteOffset, a.byteLength).toString('base64');

        /** One orbital's coefficient column, straight out of the checkpoint. */
        const readColumn = async (spinIndex: number, orbital: number): Promise<Float64Array> => {
            if (!document.orbitals) throw new Error('this document carries no orbital coefficients');
            const spin = document.orbitals.spins[spinIndex];
            if (!spin || !spin.layout) throw new Error('this spin channel has no coefficients in the file');
            const handle = await fs.promises.open(document.orbitals.path, 'r');
            try {
                return await readFchkColumn(spin.layout, document.orbitals.nBasis, orbital, async (start, length) => {
                    const buffer = Buffer.alloc(length);
                    let filled = 0;
                    // A read may legally return fewer bytes than asked for, so a
                    // short read is retried at the offset it stopped at rather
                    // than handed on to be parsed (it used to become zeros). A
                    // file that ends inside the column returns what it had, and
                    // readFchkColumn reports the truncation.
                    while (filled < length) {
                        const { bytesRead } = await handle.read(buffer, filled, length - filled, start + filled);
                        if (bytesRead <= 0) break;
                        filled += bytesRead;
                    }
                    return buffer.subarray(0, filled);
                });
            } finally {
                await handle.close();
            }
        };

        /** Build every task the request asks for; the first failure aborts all. */
        const buildJobs = async (message: { [key: string]: any }, run: { id: number }): Promise<SurfaceJob[] | string> => {
            const quality: GridQuality = message.quality === 'coarse' || message.quality === 'fine' ? message.quality : 'normal';
            const iso = typeof message.iso === 'number' && isFinite(message.iso) && message.iso > 0 ? message.iso : 0.02;
            const jobs: SurfaceJob[] = [];
            if (message.kind === 'cavity') {
                if (!document.cavity || !document.cavity.spheres || document.cavity.spheres.length === 0) {
                    return 'this file states no cavity spheres to draw';
                }
                jobs.push({
                    key: 'cavity', label: 'Solvent cavity', hue: 0,
                    task: new CavitySurfaceTask(document.cavity.spheres, quality)
                });
                return jobs;
            }
            if (!Array.isArray(message.orbitals) || message.orbitals.length === 0) return 'no orbital was selected';
            if (message.orbitals.length > 10) return 'at most ten orbitals can be drawn at once';
            for (const spec of message.orbitals) {
                const spinIndex = spec && typeof spec.spin === 'number' ? spec.spin : 0;
                const index = spec && typeof spec.index === 'number' ? spec.index : -1;
                const label = spec && typeof spec.label === 'string' ? spec.label : '';
                if (!document.orbitals || !document.orbitals.spins[spinIndex]) return 'no such spin channel';
                const energy = document.orbitals.spins[spinIndex].energies[index];
                if (index < 0 || energy === undefined) return 'no such orbital';
                const hue = typeof spec.hue === 'number' && spec.hue >= 0 ? spec.hue % 10 : jobs.length % 10;
                const coefficients = await readColumn(spinIndex, index);
                let task: SurfaceTask;
                try {
                    task = new OrbitalSurfaceTask(document.orbitals.basis, coefficients, document.data.atoms, iso, quality);
                } catch (e: any) {
                    return e && e.message ? String(e.message) : String(e);
                }
                jobs.push({
                    key: 'mo-' + spinIndex + '-' + index,
                    label: label || ('MO ' + (index + 1)),
                    hue, task
                });
            }
            if (run.id !== surfaceRunId) return '';   // superseded while reading
            return jobs;
        };

        /** Step the jobs one at a time, posting progress and finally the meshes. */
        const pumpSurfaces = async (run: { id: number; reqId: number; cancelled: boolean; jobs: SurfaceJob[] }): Promise<void> => {
            const post = (payload: { [key: string]: unknown }) => {
                if (!run.cancelled) {
                    try { webviewPanel.webview.postMessage(payload); } catch { run.cancelled = true; }
                }
            };
            try {
                for (let j = 0; j < run.jobs.length; j++) {
                    const job = run.jobs[j];
                    if (run.cancelled) return;
                    let lastPost = 0;
                    for (;;) {
                        const done = job.task.step(28);
                        if (run.cancelled) return;
                        const now = Date.now();
                        if (done || now - lastPost > 110) {
                            lastPost = now;
                            const progress = job.task.progress();
                            post({
                                command: 'surfaceProgress',
                                reqId: run.reqId,
                                id: run.id,
                                job: j,
                                jobs: run.jobs.length,
                                label: job.label,
                                fraction: progress.fraction,
                                phase: progress.phase
                            });
                        }
                        if (done) break;
                        // A turn for the message loop: this is where Cancel (or a
                        // brand new request) is seen while the grid is being filled.
                        await new Promise<void>(resolve => setTimeout(resolve, 0));
                    }
                    // One orbital at a time, not the whole batch at the end: ten
                    // orbitals at Fine are tens of megabytes of base64, and
                    // sending them one by one lets each lobe appear as soon as
                    // its grid is cut (and keeps the message below the size a
                    // webview can copy comfortably).
                    post({
                        command: 'surfaceResult', reqId: run.reqId, id: run.id, ok: true,
                        // Every answer of a multi-orbital run is appended; the
                        // `surfaceDone` that follows retires the progress line.
                        partial: run.jobs.length > 1,
                        surfaces: [{
                            key: job.key,
                            label: job.label,
                            hue: job.hue,
                            grid: job.task.gridNote(),
                            lobes: job.task.surfaces().map(lobe => ({
                                sign: lobe.sign,
                                pos: encodeFloats(lobe.mesh.positions),
                                nrm: encodeFloats(lobe.mesh.normals),
                                idx: encodeIndices(lobe.mesh.indices),
                                vertices: lobe.mesh.positions.length / 3,
                                triangles: lobe.mesh.indices.length / 3
                            }))
                        }]
                    });
                }
                post({ command: 'surfaceDone', reqId: run.reqId, id: run.id });
            } catch (e: any) {
                // A task that throws must still answer the page, or its progress
                // bar would sit at whatever fraction it had reached forever.
                post({
                    command: 'surfaceResult', reqId: run.reqId, id: run.id, ok: false,
                    message: 'The surface could not be drawn: ' + (e && e.message ? e.message : e)
                });
            } finally {
                // The tasks hold the sampled grid (tens of megabytes); releasing
                // the run here is what lets it go as soon as it is done.
                if (surfaceRun === run) surfaceRun = null;
            }
        };

        /** Everything the panel needs to draw a run's current state. */
        const mmStateMessage = (run: MmRunState, command: string): { [key: string]: unknown } => {
            const snap = run.session.finish();
            return {
                command: command,
                runId: run.id,
                step: snap.iterations,
                energy: snap.energyAfter,
                energyBefore: snap.energyBefore,
                maxGrad: snap.maxGrad,
                rmsGrad: snap.rmsGrad,
                coords: snap.coords,
                indices: snap.indices,
                points: mmNewPoints(run),
                done: run.session.done,
                converged: snap.converged,
                stopReason: snap.stopReason,
                message: snap.message
            };
        };

        /** History points the panel has not been sent yet, oldest first. */
        const mmNewPoints = (run: MmRunState): MmHistoryPoint[] => {
            const history = run.session.history;
            const points: MmHistoryPoint[] = [];
            for (let i = 0; i < history.length; i++) {
                if (history[i].step > run.lastStep) {
                    points.push(history[i]);
                }
            }
            if (points.length > 0) {
                run.lastStep = points[points.length - 1].step;
            }
            return points;
        };

        /** The finished result, in the shape the panel's last-run block reads. */
        const mmResultMessage = (runId: string, result: MmOptimizeResult): { [key: string]: unknown } => ({
            command: 'mmResult',
            runId: runId,
            ok: result.ok,
            message: result.message,
            atomCount: result.atomCount,
            indices: result.indices,
            coords: result.coords,
            energyBefore: result.energyBefore,
            energyAfter: result.energyAfter,
            iterations: result.iterations,
            maxGrad: result.maxGrad,
            rmsGrad: result.rmsGrad,
            maxMove: result.maxMove,
            converged: result.converged,
            stopReason: result.stopReason,
            counts: result.counts,
            invalid: result.invalid
        });

        /** Drop the run without finishing it (cancel, replace, dispose). */
        const mmDiscardRun = (reason: 'cancelled' | 'replaced'): MmRunState | null => {
            const run = mmRun;
            if (!run) {
                return null;
            }
            run.stopped = true;
            run.session.abandon();
            mmRun = null;
            if (reason === 'cancelled') {
                webviewPanel.webview.postMessage({ command: 'mmCancelled', runId: run.id });
            }
            return run;
        };

        // Each chunk is scheduled on a timer rather than chained through
        // `setImmediate`: the message loop must get a turn between two chunks,
        // because that is where the panel's Pause / Resume / Cancel land. A
        // timer also lets the host's own I/O be serviced first, so a Pause is
        // seen within one chunk instead of after the whole optimisation.
        //
        // The pump is handed the run it was scheduled for, and the run is what
        // the wake-up checks: a chunk chain left over from a replaced run wakes
        // up, sees it no longer owns the panel's run, and stops. Reading the
        // current run instead would let that stale chain step whatever run is
        // there — so after a cancel-and-restart the new run would be driven by
        // two chains at once, at twice the chunk rate it was asked for.
        const mmSchedulePump = (run: MmRunState): void => {
            // One chain per run. The flag is cleared at the moment the wake-up
            // actually runs, so a Resume (or any other request) arriving while
            // a wake-up is queued joins that chain instead of starting another.
            if (run.scheduled) {
                return;
            }
            run.scheduled = true;
            setTimeout(() => {
                run.scheduled = false;
                mmPump(run);
            }, 0);
        };

        // Passes the next chunk will ask for, adapted from the last one's cost
        // (see MM_CHUNK_TARGET_MS). Per panel, so one big run cannot slow the
        // next run's first chunk down.
        let mmChunkSteps = MM_CHUNK_STEPS;
        /** Size one chunk by how much work the previous one turned out to be. */
        const mmAdaptChunk = (costMs: number): void => {
            if (costMs > MM_CHUNK_TARGET_MS * MM_CHUNK_OVER) {
                mmChunkSteps = Math.max(1, mmChunkSteps > 2 ? Math.floor(mmChunkSteps / 2) : 1);
            } else if (costMs < MM_CHUNK_TARGET_MS * MM_CHUNK_UNDER && mmChunkSteps < MM_CHUNK_STEPS) {
                mmChunkSteps++;
            }
        };

        const mmPump = (run: MmRunState): void => {
            if (mmRun !== run || run.stopped || run.paused) {
                return;
            }
            try {
                const chunkStart = Date.now();
                const state = run.session.step(mmChunkSteps);
                mmAdaptChunk(Date.now() - chunkStart);
                if (mmRun !== run) {
                    // Replaced while this chunk ran (cannot happen between two
                    // synchronous statements, but the guard keeps the pump from
                    // driving a run nobody owns).
                    return;
                }
                const now = Date.now();
                if (!state.done && now - run.lastPost < MM_PROGRESS_MS) {
                    mmSchedulePump(run);
                    return;
                }
                run.lastPost = now;
                webviewPanel.webview.postMessage({
                    command: 'mmProgress',
                    runId: run.id,
                    step: state.step,
                    energy: state.energy,
                    maxGrad: state.maxGrad,
                    rmsGrad: state.rmsGrad,
                    coords: state.coords,
                    points: mmNewPoints(run),
                    done: state.done
                });
                if (state.done) {
                    run.stopped = true;
                    mmRun = null;
                    webviewPanel.webview.postMessage(mmResultMessage(run.id, run.session.finish()));
                    return;
                }
                mmSchedulePump(run);
            } catch (e: any) {
                const run2 = mmRun;
                const id = run2 ? run2.id : run.id;
                if (run2) {
                    run2.stopped = true;
                    run2.session.abandon();
                    mmRun = null;
                }
                webviewPanel.webview.postMessage({
                    command: 'mmResult',
                    runId: id,
                    ok: false,
                    message: 'Molecular Mechanics failed: ' + (e && e.message ? e.message : e)
                });
            }
        };

        webviewPanel.onDidDispose(() => {
            disposed = true;
            // Whatever the editor was waiting for can no longer be answered.
            bridge.dispose();
            if (this.panels.get(document) === bridge) this.panels.delete(document);
            if (fileWatcher) {
                fileWatcher.dispose();
                fileWatcher = undefined;
            }
            if (noticeTimer) {
                clearTimeout(noticeTimer);
                noticeTimer = undefined;
            }
            mmDiscardRun('replaced');
            // A surface run in flight has nowhere to post its meshes; stopping
            // it here is what keeps a closed tab from finishing a second of
            // arithmetic for nobody.
            if (surfaceRun) {
                surfaceRun.cancelled = true;
                surfaceRun = null;
            }
        });

        webviewPanel.webview.onDidReceiveMessage(async (message) => {
            // The page is the only sender, but a message it did not shape
            // (null, a string) must not throw before its own case runs.
            if (!message || typeof message !== 'object' || typeof message.command !== 'string') return;
            switch (message.command) {
                case 'saveFile': {
                    // One write route for every save: the toolbar's dialog, a
                    // silent re-write of the path the page last exported, and
                    // the editor's Save As target all end in writePageContent.
                    const written = await this.writePageContent(message, document.uri, forcedSaveTarget);
                    if (written) {
                        selfWrites.set(written.fsPath, Date.now());
                        // The page keeps the path: the next Save As starts
                        // beside the file just written, and a later Save
                        // re-writes it without a dialog.
                        post({ command: 'saveResult', ok: true, path: written.fsPath });
                    } else {
                        post({ command: 'saveResult', ok: false });
                    }
                    // A save the editor asked for is answered either way; a
                    // dismissed dialog leaves the document dirty.
                    if (saveWaiter) finishSave(written ? written.fsPath : undefined);
                    break;
                }
                // A data file the page exported (the powder pattern's
                // two-column text): written through the same route, but it is
                // not the structure, so it does not become the document's save
                // target and no save is reported back.
                case 'saveDataFile':
                    await this.writePageContent(message, document.uri);
                    break;
                // The page could not produce content for a save the editor
                // asked for (a dismissed dialog, an unavailable format): the
                // editor is told so, and the document stays modified.
                case 'saveCancelled':
                    if (saveWaiter) finishSave(undefined);
                    break;
                // One message per edit that survives: it is what marks the
                // document modified, shows the close prompt and gives the
                // editor's own Undo/Redo commands something to call.
                case 'edit':
                    this.changeEmitter.fire({
                        document,
                        label: typeof message.label === 'string' && message.label ? message.label : 'Edit structure',
                        undo: () => post({ command: 'hostUndo' }),
                        redo: () => post({ command: 'hostRedo' })
                    });
                    break;
                // The page's state for a hot-exit backup (see
                // backupCustomDocument).
                case 'session':
                    if (sessionWaiter) finishSession(message.session);
                    break;
                // The page's own Reload button. Reading the file again is
                // exactly the editor's Revert, and going through it is what
                // also clears the modified state this page's edits created —
                // reloading by hand would leave a clean page behind a tab that
                // still claims to be edited.
                case 'reloadRequest':
                    try {
                        await vscode.commands.executeCommand('workbench.action.files.revert');
                    } catch (e: any) {
                        // A host without the command still gets the reload; the
                        // tab keeps its dot, which the user can clear by saving.
                        try {
                            const fresh = await this.openDocument(document.uri, new vscode.CancellationTokenSource().token);
                            document.setContent(fresh);
                            await reloadPage();
                        } catch (e2: any) {
                            vscode.window.showErrorMessage('Reload failed: ' + (e2 && e2.message ? e2.message : e2));
                        }
                    }
                    break;
                // A PNG snapshot travels as base64: a webview message is the
                // page's only channel to the host, and bytes are what the file
                // has to be written as. The decode, the PNG signature and the
                // save dialog all live in savePng(), beside the Save As helper
                // that decides where the dialog opens.
                case 'savePng':
                    await this.savePng(message, document.uri);
                    break;
                case 'diffFile':
                    try {
                        const result = await vscode.window.showOpenDialog({
                            canSelectMany: false,
                            openLabel: 'Select File to Compare',
                            filters: {
                                'Molecular Files': ['gjf', 'xyz', 'extxyz', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol', 'mol2', 'log', 'out', 'fchk', 'fch', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'tcl', 'cif', 'vasp', 'poscar', 'contcar', 'cube', 'vesta', 'abi', 'abo', 'stru', 'gro', 'xtc', 'trr', 'dcd', 'trj', 'molden', 'pwi', 'pwo', 'lammpstrj', 'dump', 'lmp', 'lammps', 'data', 'md_dump', 'xdatcar', 'outcar', 'xml'],
                                'All Files': ['*']
                            }
                        });
                        if (!result || result.length === 0) {
                            webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                            break;
                        }
                        const diffUri = result[0];
                        const diffFileName = diffUri.path.split('/').pop() || 'unknown.xyz';
                        const diffExt = diffFileName.toLowerCase().split('.').pop() || '';
                        // A formatted checkpoint is compared by its structure,
                        // but it is never read as text: a checkpoint can be a
                        // gigabyte, and the streaming reader already produces
                        // the atoms the comparison needs.
                        const diffFchk = diffExt === 'fchk' || diffExt === 'fch';
                        const diffSize = await fileSize(diffUri);
                        const diffHeavy = wantsLoadProgress(diffSize, diffExt);
                        const diffFchkData = diffFchk
                            ? await withToolProgress('Molecular Viewer: reading ' + diffFileName, diffHeavy,
                                async (progress, isCancelled) => {
                                    const range = progressRange(progress, 0, 100);
                                    range(0, 'Reading the checkpoint');
                                    await yieldToHost();
                                    if (isCancelled()) throw new vscode.CancellationError();
                                    const fchk = await loadFchkStructure(diffUri, diffFileName, range);
                                    range(1);
                                    return fchk;
                                })
                            : undefined;
                        const diffLoad = diffFchk ? undefined
                            : diffHeavy
                                ? await withLoadProgress('Molecular Viewer: reading ' + diffFileName, undefined,
                                    (progress, cancelled) => readToolFile(
                                        diffUri, diffFileName, diffSize, progressRange(progress, 0, 100), cancelled))
                                : await readToolFile(diffUri, diffFileName, diffSize);
                        if (diffLoad && diffLoad.binaryUnsupported) {
                            vscode.window.showWarningMessage('Diff: this trajectory file is binary but is not an ' +
                                'XTC, TRR or DCD file, so it cannot be read as coordinates.');
                            webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                            break;
                        }
                        const diffContent = diffLoad ? diffLoad.content : new Uint8Array(0);
                        // A binary trajectory is compared by a frame the user
                        // picks: decoding it as text would only yield mojibake.
                        const diffBinary = diffLoad ? diffLoad.binary : null;
                        const diffText = diffLoad ? diffLoad.text : '';
                        const diffSetup = diffLoad ? diffLoad.setup : undefined;
                        const diffOptions = diffSetup && diffSetup.options;
                        const diffParseText = diffSetup && diffSetup.text !== undefined ? diffSetup.text : diffText;
                        // Detecting structured output is a full pass per detector;
                        // a large log would otherwise sit silent between the read
                        // and the parse. The parsed result is kept: probing and
                        // then letting `parseFile` detect again would scan the
                        // same text twice.
                        let diffStructuredData: MolecularData | undefined;
                        if (diffLoad && (diffExt === 'log' || diffExt === 'out') && diffText !== '') {
                            if (!isFrameLogOutput(diffText)) {
                                diffStructuredData = diffHeavy
                                    ? await withToolProgress('Molecular Viewer: reading ' + diffFileName, true, async (progress, isCancelled) => {
                                        const range = progressRange(progress, 0, 100);
                                        range(0, 'Reading the output');
                                        await yieldToHost();
                                        if (isCancelled()) throw new vscode.CancellationError();
                                        const structured = parseStructuredProgramOutput(diffText, diffOptions,
                                            (fraction) => range(fraction));
                                        range(1);
                                        return structured || undefined;
                                    })
                                    : parseStructuredProgramOutput(diffText, diffOptions) || undefined;
                            }
                        }

                        let diffData: MolecularData;
                        let diffFrames: (LogFrame | OrcaFrame)[] = [];
                        let diffVmdColors: string[] | null = null;
                        // A crystal comparison needs the *cell* of the right-hand
                        // file, which the display-only atom list cannot state:
                        // the base sites, their lattice and (on the page) the
                        // periodic bonds. Captured before the payload is reduced
                        // to a display list.
                        let diffCrystalPayload: {
                            a: number; b: number; c: number;
                            alpha: number; beta: number; gamma: number;
                            latticeVectors: number[][];
                            spaceGroup: string;
                            symmetryOps: string[];
                            symmetrySource?: 'file' | 'builtin';
                            baseAtoms: { element: string; x: number; y: number; z: number; occupancy?: number }[];
                            // The compared file's own bond table, by base-site
                            // index: the right half merges it over perception
                            // exactly as the left half merges the open file's.
                            fileBonds?: { atom1: number; atom2: number; order: number; shift?: number[]; valence?: number }[];
                            /** One trajectory frame's cell, not a second cell to compare. */
                            trajectory?: boolean;
                        } | null = null;

                        if (diffFchkData) {
                            // Parsed by the streaming checkpoint reader above;
                            // the comparison only needs its atoms and bonds.
                            diffData = diffFchkData;
                        } else if (diffBinary) {
                            const loaded = await withToolProgress('Molecular Viewer: reading ' + diffFileName, diffHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Decoding the trajectory');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                const decoded = await loadBinaryTrajectoryDocument(diffUri, diffFileName, diffContent, diffBinary);
                                range(1);
                                return decoded;
                            });
                            if (loaded.trajectoryFrames.length === 0 || loaded.data.atoms.length === 0) {
                                vscode.window.showErrorMessage('Diff: no structure found in ' + diffFileName);
                                webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                                break;
                            }
                            let frameIdx = 0;
                            if (loaded.trajectoryFrames.length > 1) {
                                const picked = await vscode.window.showQuickPick(
                                    loaded.trajectoryFrames.map((f, i) => ({
                                        label: `Frame ${i + 1}: ${f.stepLabel || ''}`,
                                        index: i
                                    })),
                                    { placeHolder: `Select frame to compare (${loaded.trajectoryFrames.length} frames found)`, canPickMany: false }
                                );
                                if (!picked) {
                                    webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                                    break;
                                }
                                frameIdx = picked.index;
                            }
                            const xyz = loaded.trajectoryFrames[frameIdx].xyz;
                            diffData = {
                                atoms: loaded.data.atoms.map((a, i) => ({
                                    ...a, x: xyz[i * 3], y: xyz[i * 3 + 1], z: xyz[i * 3 + 2]
                                })),
                                bonds: loaded.data.bonds,
                                title: diffFileName + ' — ' + loaded.trajectoryFrames[frameIdx].stepLabel,
                                hasExplicitBonds: loaded.data.hasExplicitBonds,
                                warnings: loaded.data.warnings
                            };
                        } else if (diffExt === 'tcl') {
                            const tclResult = parseTcl(diffText);
                            // The script is tiny; the structures it names decide
                            // whether loading them is worth a notification.
                            let tclBytes = 0;
                            for (const molecule of tclResult.scene.molecules) {
                                const sourceUri = resolveVmdSource(diffUri, molecule.sourceFile);
                                if (sourceUri) tclBytes += await fileSize(sourceUri);
                            }
                            const tclHeavy = diffHeavy || tclBytes >= LOAD_PROGRESS_MIN_BYTES;
                            const loaded = await withToolProgress('Molecular Viewer: loading ' + diffFileName, tclHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Loading the script');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                return this.loadTclScene(tclResult.scene, diffUri, range);
                            });
                            diffData = loaded.data;
                            diffFrames = loaded.frames;
                            if (loaded.payload.molecules.length > 0) {
                                // The comparison side shows the molecules the
                                // script draws, coloured the way the script
                                // colours them.
                                const molecule = loaded.payload.molecules[0];
                                diffData = {
                                    atoms: molecule.atoms.map(function (a, i) {
                                        return { element: a.element, x: a.x, y: a.y, z: a.z, index: i, beta: a.beta, occupancy: a.occupancy };
                                    }),
                                    bonds: molecule.bonds.map(function (b) {
                                        return { atom1: b.atom1, atom2: b.atom2, order: b.order };
                                    }),
                                    title: loaded.data.title,
                                    hasExplicitBonds: loaded.data.hasExplicitBonds,
                                    filePath: loaded.data.filePath
                                };
                                diffVmdColors = molecule.atoms.map(function (a) { return a.color; });
                            } else {
                                diffData = { atoms: [], bonds: [], title: 'No source', hasExplicitBonds: false };
                            }
                        } else if ((diffExt === 'log' || diffExt === 'out') && !diffStructuredData) {
                            const logResult = await withToolProgress('Molecular Viewer: reading ' + diffFileName, diffHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Reading the optimization steps');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                return parseLogFile(diffText, diffFileName,
                                    (fraction) => range(fraction,
                                        'Reading the optimization steps · ' + Math.round(fraction * 100) + '%'));
                            });
                            diffFrames = logResult.frames;
                            if (diffFrames.length > 0) {
                                diffData = ensureBonds({
                                    atoms: diffFrames[0].atoms, bonds: diffFrames[0].bonds,
                                    title: diffFrames[0].title, hasExplicitBonds: diffFrames[0].hasExplicitBonds,
                                    charge: diffFrames[0].charge, multiplicity: diffFrames[0].multiplicity
                                });
                            } else {
                                diffData = { atoms: [], bonds: [], title: 'No structures', hasExplicitBonds: false };
                            }
                        } else {
                            // A multi-frame XYZ/xTB trajectory or LAMMPS dump is
                            // a trajectory when it is opened; comparing its first
                            // frame behind the user's back would compare a
                            // geometry they never picked.
                            const frameChoices = await withToolProgress('Molecular Viewer: reading ' + diffFileName, diffHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Reading the trajectory frames');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                const frames = parseTextTrajectoryFrames(diffParseText, diffFileName, diffExt, diffOptions);
                                range(1);
                                return frames;
                            });
                            let pickedFrame: PickedTextFrame | null = null;
                            if (frameChoices.length > 1) {
                                const picked = await pickTextTrajectoryFrame(frameChoices, 'compare');
                                if (picked === 'cancelled') {
                                    webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                                    break;
                                }
                                pickedFrame = picked;
                            }
                            if (pickedFrame) {
                                diffData = ensureBonds({
                                    atoms: pickedFrame.atoms,
                                    bonds: [],
                                    title: diffFileName + ' — ' + pickedFrame.label,
                                    hasExplicitBonds: false
                                });
                                if (pickedFrame.lattice) {
                                    const cell = latticeCellParameters(pickedFrame.lattice);
                                    diffCrystalPayload = {
                                        a: cell.a, b: cell.b, c: cell.c,
                                        alpha: cell.alpha, beta: cell.beta, gamma: cell.gamma,
                                        latticeVectors: pickedFrame.lattice,
                                        spaceGroup: '',
                                        symmetryOps: [],
                                        baseAtoms: pickedFrame.atoms.map(a => ({
                                            element: a.element, x: a.x, y: a.y, z: a.z
                                        })),
                                        trajectory: true
                                    };
                                }
                            } else {
                                diffData = await withToolProgress('Molecular Viewer: reading ' + diffFileName, diffHeavy, async (progress, isCancelled) => {
                                    const range = progressRange(progress, 0, 100);
                                    range(0, 'Parsing the structure');
                                    await yieldToHost();
                                    if (isCancelled()) throw new vscode.CancellationError();
                                    const parsed = diffStructuredData || parseFile(diffParseText, diffFileName, diffOptions);
                                    // A crystal is drawn through its own cell and
                                    // states no perceived bonds here; a molecule's
                                    // bonds are part of what the comparison reads.
                                    const result = parsed.crystal ? parsed : ensureBonds(parsed);
                                    range(1);
                                    return result;
                                });
                                if (diffData.crystal) {
                                    const c = diffData.crystal;
                                    diffCrystalPayload = {
                                        a: c.a, b: c.b, c: c.c,
                                        alpha: c.alpha, beta: c.beta, gamma: c.gamma,
                                        latticeVectors: c.latticeVectors,
                                        spaceGroup: c.spaceGroup || '',
                                        symmetryOps: c.symmetryOps || [],
                                        symmetrySource: c.symmetrySource,
                                        baseAtoms: c.baseAtoms.map(a => ({
                                            element: a.element, x: a.x, y: a.y, z: a.z,
                                            occupancy: a.occupancy
                                        })),
                                        fileBonds: (c.fileBonds || []).map(b => ({
                                            atom1: b.atom1, atom2: b.atom2, order: b.order,
                                            shift: b.shift, valence: b.valence
                                        }))
                                    };
                                    diffData.bonds = [];
                                    diffData.crystal.baseBonds = [];
                                }
                            }
                            diffData.filePath = diffUri.fsPath;
                            applySetupWarnings(diffData, diffSetup);
                            if (diffSetup && diffSetup.reference) {
                                diffData.filePath = diffSetup.reference.uri.fsPath;
                            }
                        }

                        let selectedFrameIdx = 0;
                        if (diffFrames.length > 1) {
                            const items = diffFrames.map((f, i) => ({
                                label: `Frame ${i + 1}: ${f.stepLabel || ''}`,
                                description: `${f.atoms.length} atoms`,
                                index: i
                            }));
                            const picked = await vscode.window.showQuickPick(items, {
                                placeHolder: `Select frame to compare (${diffFrames.length} frames found)`,
                                canPickMany: false
                            });
                            if (!picked) {
                                webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                                break;
                            }
                            selectedFrameIdx = picked.index;
                            const sf = diffFrames[selectedFrameIdx];
                            diffData = ensureBonds({
                                atoms: sf.atoms, bonds: sf.bonds, title: sf.title,
                                hasExplicitBonds: sf.hasExplicitBonds, charge: sf.charge, multiplicity: sf.multiplicity
                            });
                            // Frames come from the loaded file, so the script's
                            // own colour recipe no longer describes them.
                            diffVmdColors = null;
                        }

                        const diffAtomColors: { [key: string]: string } = {
                            H: '#FFFFFF', C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050',
                            S: '#FFFF30', Cl: '#1FF01F', Br: '#A62929', I: '#940094', P: '#FF8000',
                            Na: '#AB5CF2', Mg: '#8AFF00', K: '#8F40D4', Ca: '#3DFF00', Fe: '#E06633',
                            Cu: '#C88033', Zn: '#7D80B0', Ag: '#C0C0C0', Au: '#FFD123', Si: '#F0C8A0',
                            B: '#FFB5B5', Li: '#CC80FF', Be: '#C2FF00', Al: '#BFA6A6'
                        };
                        const diffAtoms = diffData.atoms.map((a, i) => ({
                            element: a.element, x: a.x, y: a.y, z: a.z,
                            color: (diffVmdColors && diffVmdColors[i]) || diffAtomColors[a.element] || '#FF1493',
                            ghost: a.ghost === true ? true : undefined
                        }));
                        if (diffData.atomGroups && diffData.atomGroups.length > 0) {
                            const override: { [key: number]: string } = {};
                            for (const g of diffData.atomGroups) for (const idx of g.indices) override[idx] = g.color;
                            for (let i = 0; i < diffAtoms.length; i++) if (override[i]) diffAtoms[i].color = override[i];
                        }
                        const diffBonds = diffData.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order }));
                        const diffFramesData = diffFrames.map(f => ({
                            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: diffAtomColors[a.element] || '#FF1493', ghost: a.ghost === true ? true : undefined })),
                            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
                            stepLabel: f.stepLabel
                        }));

                        webviewPanel.webview.postMessage({
                            command: 'diffResult',
                            cancelled: false,
                            fileName: diffFileName,
                            atoms: diffAtoms,
                            bonds: diffBonds,
                            frames: diffFramesData,
                            selectedFrame: selectedFrameIdx,
                            title: diffData.title,
                            // The comparison side is a different file; its total charge
                            // is what its own bonds should be checked against, not the
                            // document on screen.
                            charge: diffData.charge,
                            // A cell, when the file stated one: the page compares it
                            // against the cell on screen before it looks at a graph.
                            crystal: diffCrystalPayload
                        });
                    } catch (e: any) {
                        // A Cancel pressed on the read/parse notification is not
                        // a failed comparison: it posts the same dismissal the
                        // page already understands.
                        if (e instanceof vscode.CancellationError) {
                            webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                            break;
                        }
                        vscode.window.showErrorMessage('Diff failed: ' + (e.message || e));
                        webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                    }
                    break;
                case 'importFile':
                    try {
                        // Step 0 — one dialog collects both the number of copies and the
                        // source: its input box (which a quick pick renders anyway) takes
                        // the copy count, and the two entries below it pick the source, so
                        // the count is settled before the structure. Dismissing the dialog
                        // cancels the import, restoring the webview status text like every
                        // other cancel path. [2026-09-16]
                        const choice = await askImportChoice();
                        if (!choice) {
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }

                        if (choice.source === 'smiles') {
                            // SMILES import: an input box validates the string on every
                            // keystroke (parse + implicit hydrogens + valence — the same
                            // lightweight check the engine runs, without the geometry
                            // stage, so it stays instant); Enter builds the 3D structure
                            // with the full SMILES -> 3D engine and feeds it into the
                            // regular importResult path (same message shape as a file
                            // import, so placement / undo / bond merging are unchanged).
                            const smiles = await vscode.window.showInputBox({
                                prompt: 'SMILES string (e.g. CCO, c1ccccc1, CC(=O)Oc1ccccc1C(=O)O)',
                                placeHolder: 'CC(=O)Oc1ccccc1C(=O)O',
                                ignoreFocusOut: true,
                                validateInput: (value: string): string | undefined => {
                                    const trimmed = value.trim();
                                    if (!trimmed) {
                                        return 'Enter a SMILES string';
                                    }
                                    const check = validateSmiles(trimmed);
                                    if (!check.ok) {
                                        return check.errorText;
                                    }
                                    return undefined;
                                }
                            });
                            if (smiles === undefined || smiles.trim() === '') {
                                webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                break;
                            }
                            const trimmed = smiles.trim();
                            const gen = smilesTo3dSafe(trimmed);
                            if (!gen.ok || !gen.info) {
                                // Cannot normally happen (validateInput already rejected
                                // bad input); kept as a guard so the view never hangs.
                                vscode.window.showErrorMessage('SMILES import failed: ' + gen.errorText);
                                webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                break;
                            }
                            const info = gen.info;
                            const impAtoms = info.atoms.map((pair) => ({
                                element: pair[1],
                                x: info.coordinates[pair[0]][0],
                                y: info.coordinates[pair[0]][1],
                                z: info.coordinates[pair[0]][2],
                                color: ATOM_COLORS[pair[1]] || '#FF1493'
                            }));
                            // Aromatic bonds carry order 1.5 — the same convention the
                            // MOL2 parser uses for 'ar' bonds — so they render as
                            // 1 solid + 1 dashed line and export as aromatic in MOL2/mol.
                            const impBonds = info.bondDetails.map((b) => ({
                                atom1: b.a,
                                atom2: b.b,
                                order: b.aromatic ? 1.5 : b.order
                            }));
                            const shown = trimmed.length > 48 ? trimmed.substring(0, 45) + '...' : trimmed;
                            postImportCopies(webviewPanel.webview, 'SMILES: ' + shown,
                                             impAtoms, impBonds, true, choice.copies, info.totalCharge);
                            break;
                        }

                        const result = await vscode.window.showOpenDialog({
                            canSelectMany: false,
                            openLabel: 'Select Structure to Import',
                            filters: {
                                'Molecular Files': ['gjf', 'xyz', 'extxyz', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol', 'mol2', 'log', 'out', 'fchk', 'fch', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'tcl', 'cif', 'vasp', 'poscar', 'contcar', 'cube', 'vesta', 'abi', 'abo', 'stru', 'gro', 'xtc', 'trr', 'dcd', 'trj', 'molden', 'pwi', 'pwo', 'lammpstrj', 'dump', 'lmp', 'lammps', 'data', 'md_dump', 'xdatcar', 'outcar', 'xml'],
                                'All Files': ['*']
                            }
                        });
                        if (!result || result.length === 0) {
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        const impUri = result[0];
                        const impFileName = impUri.path.split('/').pop() || 'unknown.xyz';
                        const impExt = impFileName.toLowerCase().split('.').pop() || '';
                        if (impExt === 'tcl') {
                            vscode.window.showErrorMessage('Import: TCL scripts are not supported (use Diff for TCL).');
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        // A formatted checkpoint is imported through the same
                        // streaming reader the direct-open path uses. Reading
                        // it as text first would allocate the whole file and
                        // then fail to parse it.
                        const impFchk = impExt === 'fchk' || impExt === 'fch';
                        const impSize = await fileSize(impUri);
                        const impHeavy = wantsLoadProgress(impSize, impExt);
                        const impFchkData = impFchk
                            ? await withToolProgress('Molecular Viewer: reading ' + impFileName, impHeavy,
                                async (progress, isCancelled) => {
                                    const range = progressRange(progress, 0, 100);
                                    range(0, 'Reading the checkpoint');
                                    await yieldToHost();
                                    if (isCancelled()) throw new vscode.CancellationError();
                                    const fchk = await loadFchkStructure(impUri, impFileName, range);
                                    range(1);
                                    return fchk;
                                })
                            : undefined;
                        const impLoad = impFchk ? undefined
                            : impHeavy
                                ? await withLoadProgress('Molecular Viewer: reading ' + impFileName, undefined,
                                    (progress, cancelled) => readToolFile(
                                        impUri, impFileName, impSize, progressRange(progress, 0, 100), cancelled))
                                : await readToolFile(impUri, impFileName, impSize);
                        if (impLoad && impLoad.binaryUnsupported) {
                            vscode.window.showWarningMessage('Import: this trajectory file is binary but is not an ' +
                                'XTC, TRR or DCD file, so it cannot be read as coordinates.');
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        const impContent = impLoad ? impLoad.content : new Uint8Array(0);
                        const impBinary = impLoad ? impLoad.binary : null;
                        const impText = impLoad ? impLoad.text : '';
                        const impSetup = impLoad ? impLoad.setup : undefined;
                        const impOptions = impSetup && impSetup.options;
                        const impParseText = impSetup && impSetup.text !== undefined ? impSetup.text : impText;
                        // Same structured-output pass as Diff: it must not run
                        // between two notifications without one of its own. The
                        // parsed result is reused below so the same text is not
                        // scanned a second time by `parseFile`.
                        let impStructuredData: MolecularData | undefined;
                        if (impLoad && (impExt === 'log' || impExt === 'out') && impText !== '') {
                            if (!isFrameLogOutput(impText)) {
                                impStructuredData = impHeavy
                                    ? await withToolProgress('Molecular Viewer: reading ' + impFileName, true, async (progress, isCancelled) => {
                                        const range = progressRange(progress, 0, 100);
                                        range(0, 'Reading the output');
                                        await yieldToHost();
                                        if (isCancelled()) throw new vscode.CancellationError();
                                        const structured = parseStructuredProgramOutput(impText, impOptions,
                                            (fraction) => range(fraction));
                                        range(1);
                                        return structured || undefined;
                                    })
                                    : parseStructuredProgramOutput(impText, impOptions) || undefined;
                            }
                        }

                        let impData: MolecularData | undefined = impFchkData;
                        if (impFchkData) {
                            // Already parsed above; nothing left to read.
                        } else if (impBinary) {
                            const loaded = await withToolProgress('Molecular Viewer: reading ' + impFileName, impHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Decoding the trajectory');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                const decoded = await loadBinaryTrajectoryDocument(impUri, impFileName, impContent, impBinary);
                                range(1);
                                return decoded;
                            });
                            if (loaded.trajectoryFrames.length === 0 || loaded.data.atoms.length === 0) {
                                vscode.window.showErrorMessage('Import: no structure found in ' + impFileName);
                                webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                break;
                            }
                            let frameIdx = 0;
                            if (loaded.trajectoryFrames.length > 1) {
                                const picked = await vscode.window.showQuickPick(
                                    loaded.trajectoryFrames.map((f, i) => ({
                                        label: `Frame ${i + 1}: ${f.stepLabel || ''}`,
                                        index: i
                                    })),
                                    { placeHolder: `Select frame to import (${loaded.trajectoryFrames.length} frames found)`, canPickMany: false }
                                );
                                if (!picked) {
                                    webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                    break;
                                }
                                frameIdx = picked.index;
                            }
                            const xyz = loaded.trajectoryFrames[frameIdx].xyz;
                            impData = {
                                atoms: loaded.data.atoms.map((a, i) => ({
                                    ...a, x: xyz[i * 3], y: xyz[i * 3 + 1], z: xyz[i * 3 + 2]
                                })),
                                bonds: loaded.data.bonds,
                                title: impFileName + ' — ' + loaded.trajectoryFrames[frameIdx].stepLabel,
                                hasExplicitBonds: loaded.data.hasExplicitBonds,
                                warnings: loaded.data.warnings
                            };
                        } else if ((impExt === 'log' || impExt === 'out') && !impStructuredData) {
                            const logResult = await withToolProgress('Molecular Viewer: reading ' + impFileName, impHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Reading the optimization steps');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                return parseLogFile(impText, impFileName,
                                    (fraction) => range(fraction,
                                        'Reading the optimization steps · ' + Math.round(fraction * 100) + '%'));
                            });
                            if (logResult.frames.length === 0) {
                                vscode.window.showErrorMessage('Import: no structure found in ' + impFileName);
                                webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                break;
                            }
                            let frameIdx = 0;
                            if (logResult.frames.length > 1) {
                                const items = logResult.frames.map((f, i) => ({
                                    label: `Frame ${i + 1}: ${f.stepLabel || ''}`,
                                    description: `${f.atoms.length} atoms`,
                                    index: i
                                }));
                                const picked = await vscode.window.showQuickPick(items, {
                                    placeHolder: `Select frame to import (${logResult.frames.length} frames found)`,
                                    canPickMany: false
                                });
                                if (!picked) {
                                    webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                    break;
                                }
                                frameIdx = picked.index;
                            }
                            const sf = logResult.frames[frameIdx];
                            impData = ensureBonds({
                                atoms: sf.atoms, bonds: sf.bonds, title: sf.title,
                                hasExplicitBonds: sf.hasExplicitBonds, charge: sf.charge, multiplicity: sf.multiplicity
                            });
                        } else {
                            // A multi-frame XYZ/xTB trajectory or LAMMPS dump is
                            // opened as a trajectory, so importing it asks which
                            // frame to bring in rather than taking the first.
                            const frameChoices = await withToolProgress('Molecular Viewer: reading ' + impFileName, impHeavy, async (progress, isCancelled) => {
                                const range = progressRange(progress, 0, 100);
                                range(0, 'Reading the trajectory frames');
                                await yieldToHost();
                                if (isCancelled()) throw new vscode.CancellationError();
                                const frames = parseTextTrajectoryFrames(impParseText, impFileName, impExt, impOptions);
                                range(1);
                                return frames;
                            });
                            let pickedFrame: PickedTextFrame | null = null;
                            if (frameChoices.length > 1) {
                                const picked = await pickTextTrajectoryFrame(frameChoices, 'import');
                                if (picked === 'cancelled') {
                                    webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                                    break;
                                }
                                pickedFrame = picked;
                            }
                            if (pickedFrame) {
                                impData = ensureBonds({
                                    atoms: pickedFrame.atoms,
                                    bonds: [],
                                    title: impFileName + ' — ' + pickedFrame.label,
                                    hasExplicitBonds: false
                                });
                            } else {
                                impData = await withToolProgress('Molecular Viewer: reading ' + impFileName, impHeavy, async (progress, isCancelled) => {
                                    const range = progressRange(progress, 0, 100);
                                    range(0, 'Parsing the structure');
                                    await yieldToHost();
                                    if (isCancelled()) throw new vscode.CancellationError();
                                    const parsed = impStructuredData || parseFile(impParseText, impFileName, impOptions);
                                    // Crystal files (CIF/VASP/VESTA/Cube): import the
                                    // unit-cell content as a molecule; periodicity is
                                    // dropped (the current structure stays molecular).
                                    let result: MolecularData;
                                    if (parsed.crystal) {
                                        const cellAtoms = (parsed.atoms.length > 0 ? parsed.atoms : parsed.crystal.baseAtoms)
                                            .map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, index: 0 }));
                                        result = ensureBonds({ atoms: cellAtoms, bonds: [], title: parsed.title, hasExplicitBonds: false });
                                    } else {
                                        result = ensureBonds(parsed);
                                    }
                                    range(1);
                                    return result;
                                });
                            }
                        }
                        if (!impData || impData.atoms.length === 0) {
                            vscode.window.showErrorMessage('Import: no atoms found in ' + impFileName);
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        applySetupWarnings(impData, impSetup);

                        const impAtoms = impData.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: ATOM_COLORS[a.element] || '#FF1493', ghost: a.ghost === true ? true : undefined }));
                        const impBonds = impData.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order }));
                        postImportCopies(webviewPanel.webview, impFileName,
                                         impAtoms, impBonds, !!impData.hasExplicitBonds, choice.copies,
                                         impData.charge, impData.multiplicity);
                    } catch (e: any) {
                        // A Cancel pressed on the read/parse notification is not
                        // a failed import: it posts the same dismissal the page
                        // already understands.
                        if (e instanceof vscode.CancellationError) {
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        vscode.window.showErrorMessage('Import failed: ' + (e.message || e));
                        webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                    }
                    break;
                case 'addHydrogens':
                    // Completion runs here, where the SMILES -> 3D engine's
                    // experience tables live, so an added X-H bond has the same
                    // length, and an added H-X-Y angle the same value, that a
                    // generated structure would carry. The structure is read
                    // only: the answer is the new hydrogens, never a moved atom.
                    try {
                        const hResult = addHydrogens(message.atoms, message.bonds, message.selected);
                        webviewPanel.webview.postMessage({
                            command: 'hydrogenResult',
                            ok: true,
                            // The request's identity and the structure version it
                            // was computed from travel back, so the page can drop
                            // an answer whose frame/molecule/edit is gone.
                            requestId: typeof message.requestId === 'number' ? message.requestId : undefined,
                            revision: typeof message.revision === 'string' ? message.revision : undefined,
                            atomCount: Array.isArray(message.atoms) ? message.atoms.length : 0,
                            atoms: hResult.atoms,
                            bonds: hResult.bonds,
                            centres: hResult.centres,
                            message: hResult.message
                        });
                    } catch (e: any) {
                        webviewPanel.webview.postMessage({
                            command: 'hydrogenResult',
                            ok: false,
                            requestId: typeof message.requestId === 'number' ? message.requestId : undefined,
                            revision: typeof message.revision === 'string' ? message.revision : undefined,
                            message: 'Add H failed: ' + (e && e.message ? e.message : e)
                        });
                    }
                    break;
                case 'info':
                    vscode.window.showInformationMessage(message.text);
                    break;
                case 'mmOptimize':
                    // Molecular mechanics runs here, in the extension host: the
                    // work is bounded but heavy for a large selection, and
                    // running it off the webview keeps the view responsive while
                    // it works. It is stepped in chunks (see mmPump), so the
                    // panel watches the energy and the coordinates instead of
                    // waiting for one final answer, and it can stop the run —
                    // the answer carries new coordinates for the selected atoms
                    // only; every other atom is a fixed source of force and is
                    // not returned at all.
                    try {
                        // A new run replaces one already in flight for this
                        // panel: the old session is dropped (its answer is no
                        // longer wanted) and told to the panel as cancelled.
                        mmDiscardRun('replaced');
                        const runId = typeof message.runId === 'string' && message.runId.length > 0
                            ? message.runId
                            : 'mm-' + Date.now();
                        // The panel's settings travel with the request; the module
                        // validates every one of them, so an unexpected payload
                        // simply falls back to its own defaults.
                        const began = beginOptimization(message.atoms, message.bonds, message.selected, message.options);
                        if (began.result !== null) {
                            // Refused, or a selection with nothing to optimise:
                            // one complete answer, exactly as before.
                            webviewPanel.webview.postMessage(mmResultMessage(runId, began.result));
                            break;
                        }
                        const session = began.session as MmSession;
                        const run: MmRunState = {
                            id: runId,
                            session: session,
                            paused: false,
                            stopped: false,
                            scheduled: false,
                            lastStep: -1,
                            lastPost: Date.now()
                        };
                        mmRun = run;
                        // A new run starts optimistic and shrinks its chunk only
                        // if its own first chunks turn out to be expensive.
                        mmChunkSteps = MM_CHUNK_STEPS;
                        // The first message carries the input geometry (step 0
                        // of the history), so the curve starts where the run did.
                        webviewPanel.webview.postMessage(mmStateMessage(run, 'mmProgress'));
                        mmSchedulePump(run);
                    } catch (e: any) {
                        webviewPanel.webview.postMessage({
                            command: 'mmResult',
                            runId: typeof message.runId === 'string' ? message.runId : '',
                            ok: false,
                            message: 'Molecular Mechanics failed: ' + (e && e.message ? e.message : e)
                        });
                    }
                    break;
                case 'mmPause': {
                    // Stop pumping and keep the session: the panel is told where
                    // the run stands, and `mmResume` continues from exactly this
                    // state (the search's step size, conjugate direction and
                    // pair list are all part of it, not of an iteration).
                    // Idempotent, and it always answers: a pause whose reply was
                    // lost (or that arrived twice) can never leave the panel
                    // waiting for an answer that will not come.
                    const run = mmRun;
                    if (run && !run.stopped && run.id === message.runId) {
                        if (!run.paused) {
                            run.paused = true;
                            run.session.pause();
                        }
                        webviewPanel.webview.postMessage(mmStateMessage(run, 'mmPaused'));
                    } else if (typeof message.runId === 'string' && message.runId.length > 0) {
                        // No run of that id is in flight here — it already
                        // finished, or a newer run replaced it between the panel's
                        // click and this message. The panel flipped to Pause
                        // optimistically, so it is waiting to be told something;
                        // answering that the run is over is what lets it return to
                        // Start instead of sitting on a button that does nothing.
                        webviewPanel.webview.postMessage({ command: 'mmCancelled', runId: message.runId });
                    }
                    break;
                }
                case 'mmResume': {
                    const run = mmRun;
                    if (run && !run.stopped && run.id === message.runId) {
                        if (run.paused) {
                            run.paused = false;
                            run.lastPost = 0;
                            // The time spent paused is added back to the run's
                            // wall-clock budget, so a pause cannot eat it.
                            run.session.resume();
                        }
                        // Answer the request either way, so the panel hears back
                        // even when it asked to resume a run that was already
                        // going (its optimistic state flip is then confirmed).
                        webviewPanel.webview.postMessage(mmStateMessage(run, 'mmProgress'));
                        // The scheduler refuses a second chain when one is
                        // already queued (run.scheduled), so resuming a run
                        // that never stopped cannot step the session twice.
                        mmSchedulePump(run);
                    } else if (typeof message.runId === 'string' && message.runId.length > 0) {
                        // Same as the pause above: a resume for a run that is no
                        // longer here is answered, so the panel's optimistic
                        // "running" cannot outlive the run it belongs to.
                        webviewPanel.webview.postMessage({ command: 'mmCancelled', runId: message.runId });
                    }
                    break;
                }
                case 'mmCancel': {
                    // The run is dropped; whatever geometry it produced is
                    // already on screen, and the panel is told so it can stop
                    // waiting for progress.
                    const run = mmRun;
                    if (run && run.id === message.runId) {
                        mmDiscardRun('cancelled');
                    }
                    break;
                }
                case 'surfaceDraw': {
                    // A new request supersedes the one in flight: its answer is
                    // no longer wanted and the panel has already cleared the
                    // picture, so it must not arrive later and draw something
                    // the user did not ask for.
                    if (surfaceRun) surfaceRun.cancelled = true;
                    const id = ++surfaceRunId;
                    const reqId = typeof message.reqId === 'number' ? message.reqId : 0;
                    const run: { id: number; reqId: number; cancelled: boolean; jobs: SurfaceJob[] } = { id, reqId, cancelled: false, jobs: [] };
                    surfaceRun = run;
                    try {
                        const built = await buildJobs(message, run);
                        if (typeof built === 'string') {
                            // '' means "superseded while reading", which is not
                            // a failure: the newer request will answer.
                            if (built !== '') {
                                webviewPanel.webview.postMessage({ command: 'surfaceResult', reqId, id, ok: false, message: built });
                            }
                            break;
                        }
                        run.jobs = built;
                        void pumpSurfaces(run);
                    } catch (e: any) {
                        webviewPanel.webview.postMessage({
                            command: 'surfaceResult', reqId, id,
                            ok: false,
                            message: 'The surface could not be drawn: ' + (e && e.message ? e.message : e)
                        });
                    }
                    break;
                }
                case 'surfaceCancel':
                    if (surfaceRun) {
                        surfaceRun.cancelled = true;
                        surfaceRun = null;
                    }
                    break;
                case 'error':
                    vscode.window.showErrorMessage(message.text);
                    break;
            }
        });
    }

    /** Where a Save dialog starts: the folder the structure came from, with the suggested name. */
    private saveDialogUri(filePath: unknown, suggestedName: unknown, fallbackName: string, documentUri?: vscode.Uri): vscode.Uri {
        const srcPath = typeof filePath === 'string' ? filePath : '';
        const srcDir = srcPath ? srcPath.substring(0, srcPath.replace(/\\/g, '/').lastIndexOf('/')) : '';
        const name = typeof suggestedName === 'string' && suggestedName ? suggestedName : fallbackName;
        if (srcDir) return vscode.Uri.file(srcDir + '/' + name);
        // A structure built in the page (a merged file, a pasted fragment) has
        // no source path; a bare file name would put the dialog at the drive
        // root (or the process's current directory). Start beside the document
        // the panel was opened from, then in the workspace folder.
        if (documentUri && documentUri.scheme === 'file' && documentUri.fsPath) {
            return vscode.Uri.joinPath(vscode.Uri.file(path.dirname(documentUri.fsPath)), name);
        }
        const folder = vscode.workspace.workspaceFolders && vscode.workspace.workspaceFolders[0];
        if (folder) return vscode.Uri.joinPath(folder.uri, name);
        return vscode.Uri.file(name);
    }

    /**
     * The file a Save As actually writes.
     *
     * The page states the extension of the content it generated, because the
     * save dialog's filter is only a suggestion: a user can rename a Gaussian
     * input to `water.xyz` or save it through All Files, and the dialog does
     * not tell the host which filter was current. A name typed without an
     * extension gets the content's own extension, and a name whose extension
     * contradicts the content is confirmed rather than written as a file whose
     * suffix lies about what it holds. A format with no extension of its own
     * (STRU) is written exactly as named, and a page that sends no extension
     * keeps the old behaviour.
     */
    private async resolveSaveTarget(uri: vscode.Uri, extension: unknown): Promise<vscode.Uri | undefined> {
        const raw = typeof extension === 'string' ? extension.trim() : '';
        const expected = /^\.[A-Za-z0-9]{1,16}$/.test(raw) ? raw.toLowerCase() : '';
        if (!expected) return uri;
        const pathText = uri.path;
        const slash = pathText.lastIndexOf('/');
        const name = pathText.substring(slash + 1);
        const dot = name.lastIndexOf('.');
        const actual = dot >= 0 ? name.substring(dot).toLowerCase() : '';
        // A filter can offer more than one suffix for the same format (a VASP
        // cell may be `.vasp`, `.poscar` or `.contcar`), and a name without any
        // extension can still be the format's conventional one (`POSCAR`).
        const accepted = SAVE_EXTENSION_ALIASES[expected] || [expected];
        if (accepted.indexOf(actual) >= 0) return uri;
        const dir = pathText.substring(0, slash + 1);
        if (actual === '') {
            const conventional = SAVE_EXTENSIONLESS_NAMES[expected] || [];
            if (conventional.indexOf(name.toLowerCase()) >= 0) return uri;
            return uri.with({ path: dir + name + expected });
        }
        const corrected = dir + name.substring(0, name.length - actual.length) + expected;
        const keepLabel = 'Keep ' + actual;
        const choice = await vscode.window.showWarningMessage(
            'The file name ends in "' + actual + '", but the content is ' +
            expected.substring(1).toUpperCase() + ' format.',
            'Save as ' + expected, keepLabel);
        if (choice === 'Save as ' + expected) return uri.with({ path: corrected });
        if (choice === keepLabel) return uri;
        // A dismissed warning writes no file, exactly like Cancel.
        return undefined;
    }

    /**
     * Write what the page produced.
     *
     * Three callers, one write: the toolbar's Save As (a dialog with the
     * page's own format filters), a save the editor asked for where the page
     * chose to re-write the file it last exported (a target, no dialog), and
     * the editor's Save As (a target the editor's own dialog chose). The
     * written path goes back to the page, which is what makes the next Save As
     * start beside the file just written instead of beside the file opened.
     */
    private async writePageContent(
        message: { content?: unknown; extension?: unknown; suggestedName?: unknown; filePath?: unknown; targetPath?: unknown },
        documentUri: vscode.Uri,
        forcedTarget?: vscode.Uri
    ): Promise<vscode.Uri | undefined> {
        try {
            // The page is the only sender, but content that is not text would
            // be encoded as "[object Object]" and leave a file that no reader
            // can open.
            if (typeof message.content !== 'string') {
                vscode.window.showErrorMessage('Save failed: the page sent no file content.');
                return undefined;
            }
            let target = forcedTarget;
            if (!target) {
                // A path the page already holds from an earlier save is
                // written as it stands: the format was chosen then, and asking
                // again would only repeat the dialog.
                const known = typeof message.targetPath === 'string' ? message.targetPath : '';
                if (known) {
                    target = vscode.Uri.file(known);
                } else {
                    const uri = await vscode.window.showSaveDialog({
                        defaultUri: this.saveDialogUri(message.filePath, message.suggestedName, 'molecule.xyz', documentUri),
                        filters: {
                            'XYZ Files': ['xyz'],
                            'Gaussian Input': ['gjf'],
                            'CIF Files': ['cif'],
                            'VASP POSCAR': ['vasp', 'poscar', 'contcar'],
                            'Gaussian Cube': ['cube'],
                            'VESTA': ['vesta'],
                            'Turbomole Coord': ['coord'],
                            'ORCA Input': ['inp'],
                            'MOL2': ['mol2'],
                            'MDL Mol': ['mol'],
                            'PDB': ['pdb'],
                            'MOPAC': ['mop'],
                            'ABACUS STRU': ['STRU'],
                            'ABINIT Input': ['abi'],
                            'CP2K Input': ['inp'],
                            'Quantum ESPRESSO Input': ['in', 'pwi'],
                            'LAMMPS Data': ['lmp', 'data'],
                            'GROMACS': ['gro'],
                            'Extended XYZ': ['xyz'],
                            // The powder XRD panel exports a two-column
                            // 2theta/intensity text file.
                            'XY Data': ['xy', 'dat'],
                            'All Files': ['*']
                        }
                    });
                    if (!uri) return undefined;
                    target = await this.resolveSaveTarget(uri, message.extension);
                    if (!target) return undefined;
                }
            }
            await vscode.workspace.fs.writeFile(target, new TextEncoder().encode(message.content));
            vscode.window.showInformationMessage('Saved: ' + target.fsPath);
            return target;
        } catch (e: any) {
            vscode.window.showErrorMessage('Save failed: ' + (e && e.message ? e.message : e));
            return undefined;
        }
    }

    /**
     * Write the PNG the page captured.
     *
     * The page sends base64 without a data-URL prefix, because that is what a
     * canvas hands out and the message is its only transport. The bytes are
     * written exactly as decoded — the page is the only thing that knows what it
     * drew — but they have to be a PNG: a truncated message, an empty canvas
     * read or a corrupted string must fail loudly instead of leaving a file the
     * user will find unopenable later.
     */
    private async savePng(message: { data?: unknown; suggestedName?: unknown; filePath?: unknown }, documentUri?: vscode.Uri): Promise<void> {
        const base64 = typeof message.data === 'string' ? message.data : '';
        if (!base64) {
            vscode.window.showErrorMessage('Snapshot failed: the page sent no image data.');
            return;
        }
        // 96 MB of base64 is ~72 MB of PNG — already far past a 20 MP capture,
        // so a larger message is a corrupted one rather than a big picture.
        if (base64.length > 96 * 1024 * 1024) {
            vscode.window.showErrorMessage('Snapshot failed: the image is too large to write.');
            return;
        }
        try {
            const bytes = Buffer.from(base64, 'base64');
            // 0x89 'P' 'N' 'G': the eight-byte signature every PNG file starts with.
            if (bytes.length < 8 || bytes[0] !== 0x89 || bytes[1] !== 0x50 || bytes[2] !== 0x4e || bytes[3] !== 0x47) {
                throw new Error('the captured data is not a PNG image');
            }
            const uri = await vscode.window.showSaveDialog({
                defaultUri: this.saveDialogUri(message.filePath, message.suggestedName, 'structure-snapshot.png', documentUri),
                filters: {
                    'PNG Image': ['png'],
                    'All Files': ['*']
                }
            });
            if (!uri) {
                return;
            }
            // A name typed without an extension would otherwise become a PNG with
            // no .png on it, which no image viewer would offer to open.
            const target = /[^/\\]+\.[^/\\]+$/.test(uri.fsPath) ? uri : vscode.Uri.file(uri.fsPath + '.png');
            await vscode.workspace.fs.writeFile(target, bytes);
            vscode.window.showInformationMessage('Saved: ' + target.fsPath);
        } catch (e: any) {
            vscode.window.showErrorMessage('Snapshot failed: ' + (e && e.message ? e.message : e));
        }
    }
}

class MolecularDocument implements vscode.CustomDocument {
    constructor(
        public readonly uri: vscode.Uri,
        public data: MolecularData,
        public frames: (LogFrame | OrcaFrame)[] = [],
        public optSteps?: OptStep[],
        public normalModes?: NormalMode[],
        public logSource?: string,
        public routes?: RouteSection[],
        public vmdPayload?: VmdPayload,
        public autoBeta: boolean = false,
        public optThresholds?: OptThresholds,
        public trajectoryFrames?: TrajectoryFramePayload[],
        public charges?: LogChargeSet[],
        public summary?: LogSummary[],
        /**
         * The orbitals a checkpoint carries, coefficients included. Held on the
         * document rather than on the payload: the coefficient block is a
         * scalar product of the basis size, and the page never needs it — a
         * surface arrives as a mesh (see `pumpSurfaces`).
         */
        public orbitals?: FchkOrbitals,
        /** The cavity the file states, spheres included when it lists them. */
        public cavity?: SolventCavity,
        /** Energies and occupations only, for the level diagram. */
        public orbitalLevels?: OrbitalLevels,
        /**
         * A relaxed scan and/or an IRC path. One slot rather than two, because
         * the page draws them in the same card as the convergence charts: what
         * a file holds decides which of the three that card shows.
         */
        public pathways?: LogPathways,
        /**
         * What the page held when the window was closed with unsaved edits —
         * the structure on screen, the electronic state and the frames that
         * were edited. Absent unless VS Code handed a backup id back to
         * `openCustomDocument`, in which case the page applies it instead of
         * showing the file's first frame.
         */
        public restoreSession?: unknown
    ) {}

    /**
     * Take on another read of the same file, in place.
     *
     * A panel holds this object, not a snapshot of its fields: a revert and the
     * external-change reload both re-read the file and then rebuild the page
     * from the same document, so the new content has to land on the instance
     * the panel already has.
     */
    setContent(other: MolecularDocument): void {
        this.data = other.data;
        this.frames = other.frames;
        this.optSteps = other.optSteps;
        this.normalModes = other.normalModes;
        this.logSource = other.logSource;
        this.routes = other.routes;
        this.vmdPayload = other.vmdPayload;
        this.autoBeta = other.autoBeta;
        this.optThresholds = other.optThresholds;
        this.trajectoryFrames = other.trajectoryFrames;
        this.charges = other.charges;
        this.summary = other.summary;
        this.orbitals = other.orbitals;
        this.cavity = other.cavity;
        this.orbitalLevels = other.orbitalLevels;
        this.pathways = other.pathways;
        // A re-read is the file as it is now: the backup that was applied when
        // the tab opened belongs to the previous reading and must not be put
        // back on top of a revert or a reload.
        this.restoreSession = undefined;
    }

    dispose(): void {}
}

/** The level diagram's view of a checkpoint's orbitals. */
function orbitalLevelsOf(orbitals: FchkOrbitals): OrbitalLevels {
    return {
        spins: orbitals.spins.map(s => ({
            label: s.label,
            energies: s.energies,
            occupations: s.occupations,
            occupied: s.occupied
        })),
        closedShell: orbitals.spins.length < 2
    };
}
