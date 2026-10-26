import * as vscode from 'vscode';
import * as path from 'path';
import { parseFile, parseLogFile, parseXyzTrajectory, LogFrame, OrcaFrame, RouteSection, parseTcl, resolveVmdScene, resolveSelection, vdwRadius, hasScalarChannel, parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg } from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { BOND_PERCEPTION_SOURCE } from '../tools/bondPerception';
import { addHydrogens } from '../tools/addHydrogens';
import { beginOptimization, describeForceField, MmHistoryPoint, MmOptimizeResult, MmSession } from '../tools/mmOptimizer';
import { smilesTo3dSafe, validateSmiles } from '../tools/smilesTo3d';
import { MolecularData, AtomGroup, OptStep, NormalMode, OptThresholds, VmdScene, VmdMolecule } from '../types';

const ATOM_COLORS: { [key: string]: string } = {
    H: '#FFFFFF', He: '#D9FFFF', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5',
    C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050', Ne: '#B3E3F5',
    Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6', Si: '#F0C8A0', P: '#FF8000',
    S: '#FFFF30', Cl: '#1FF01F', Ar: '#80D1E3', K: '#8F40D4', Ca: '#3DFF00',
    Sc: '#E6E6E6', Ti: '#BFC2C7', V: '#A6A6AB', Cr: '#8A99C7', Mn: '#9C7AC7',
    Fe: '#E06633', Co: '#F090A0', Ni: '#50D050', Cu: '#C88033', Zn: '#7D80B0',
    Ga: '#C28F8F', Ge: '#668F8F', As: '#BD80E3', Se: '#FFA100', Br: '#A62929',
    Kr: '#5CB8D1', Rb: '#702EB0', Sr: '#00FF00', Y: '#94FFFF', Zr: '#94E0E0',
    Nb: '#73C2C9', Mo: '#54B5B5', Tc: '#3B9E9E', Ru: '#248F8F', Rh: '#0A7D8C',
    Pd: '#006985', Ag: '#C0C0C0', Cd: '#FFD98F', In: '#A67573', Sn: '#668080',
    Sb: '#9E63B5', Te: '#D47A00', I: '#940094', Xe: '#429EB0', Cs: '#57178F',
    Ba: '#00C900', La: '#70D4FF', Ce: '#FFFFC7', Pr: '#D9FFC7', Nd: '#C7FFC7',
    Pm: '#A3FFC7', Sm: '#8FFFC7', Eu: '#61FFC7', Gd: '#45FFC7', Tb: '#30FFC7',
    Dy: '#1FFFC7', Ho: '#00FF9C', Er: '#00E675', Tm: '#00D452', Yb: '#00BF38',
    Lu: '#00AB24', Hf: '#4DC2FF', Ta: '#4DA6FF', W: '#2194D6', Re: '#267DAB',
    Os: '#266696', Ir: '#175487', Pt: '#D0D0E0', Au: '#FFD123', Hg: '#B8B8D0',
    Tl: '#A6544D', Pb: '#575961', Bi: '#9E4FB5', Po: '#AB5C00', At: '#754F45',
    Rn: '#428296', Fr: '#420066', Ra: '#007D00', Ac: '#70ABFA', Th: '#00BAFF',
    Pa: '#00A1FF', U: '#008FFF', Np: '#0080FF', Pu: '#006BFF', Am: '#545CF2',
    Cm: '#785CE3', Bk: '#8A4FE3', Cf: '#A136D4', Es: '#B31FD4', Fm: '#B31FBA',
    Md: '#B30DA6', No: '#BD0D87', Lr: '#C70066', Rf: '#CC0059', Db: '#D9004F',
    Sg: '#E00045', Bh: '#E6002E', Hs: '#EB0026'
};

/** 导入消息里的原子载荷 (客户端 importStructure 逐项写入 MD.atoms)。 */
interface ImportAtomPayload {
    element: string;
    x: number;
    y: number;
    z: number;
    color: string;
}

/** 导入消息里的键载荷 (索引相对本次导入的原子数组)。 */
interface ImportBondPayload {
    atom1: number;
    atom2: number;
    order: number;
}

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
 * Payload of one molecule of a VMD visualization script: the structure the
 * script loads plus every drawing property the script's representations imply.
 * Everything a representation decides is resolved host-side (colours, sphere
 * radii, bond radii, surface patches) so the webview only has to draw it.
 */
interface VmdAtomPayload {
    element: string;
    x: number;
    y: number;
    z: number;
    color: string;
    beta?: number;
    occupancy?: number;
    /** Sphere size as a multiplier on the viewer's own radius; 0 = not drawn. */
    vmdScale?: number;
    /** van der Waals radius in Å — the reference length a `Surf` probe adds to. */
    vdw: number;
}

interface VmdBondPayload {
    atom1: number;
    atom2: number;
    order: number;
    /** Cylinder radius in Å from the script; undefined = the viewer's own. */
    vmdR?: number;
}

/** One translucent surface patch (a `Surf` representation). */
interface VmdSegmentPayload {
    color: string;
    indices: number[];
    probe: number;
    /** The egg this patch draws, when the structure sits next to an egg/unit file. */
    egg?: { name: string; units: number[]; level?: string };
}

interface VmdRepPayload {
    style: string;
    colorMode: string;
    colorArg?: number;
    scaleMin?: number;
    scaleMax?: number;
    material: string;
    selection: string;
    /** Number of atoms this representation selects. */
    atoms: number;
}

interface VmdMoleculePayload {
    label: string;
    /** Resolved absolute path of the structure file the script loads. */
    sourceFile: string;
    atoms: VmdAtomPayload[];
    bonds: VmdBondPayload[];
    segments: VmdSegmentPayload[];
    reps: VmdRepPayload[];
    warnings: string[];
}

interface VmdPayload {
    molecules: VmdMoleculePayload[];
    colorScaleMethod: string;
    warnings: string[];
}

/** Everything one structure file contributes, whichever way it was reached. */
interface LoadedStructure {
    data: MolecularData;
    frames: (LogFrame | OrcaFrame)[];
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    routes?: RouteSection[];
    optThresholds?: OptThresholds;
    logSource?: string;
    error?: string;
}

const EMPTY_STRUCTURE: MolecularData = { atoms: [], bonds: [], title: 'Empty', hasExplicitBonds: false };

/** Load and parse one structure file (any supported format). */
async function loadStructureFile(uri: vscode.Uri): Promise<LoadedStructure> {
    const content = await vscode.workspace.fs.readFile(uri);
    const text = new TextDecoder().decode(content);
    const fileName = uri.path.split('/').pop() || 'unknown.xyz';
    const ext = fileName.toLowerCase().split('.').pop() || '';

    if (ext === 'log' || ext === 'out') {
        const logResult = parseLogFile(text, fileName);
        if (logResult.frames.length === 0) {
            return { data: { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false },
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
            logSource: logResult.source
        };
    }

    const parsed = parseFile(text, fileName);
    let data: MolecularData;
    if (parsed.crystal) {
        parsed.bonds = [];
        parsed.crystal.baseBonds = [];
        data = parsed;
    } else {
        data = ensureBonds(parsed);
    }
    data.filePath = uri.fsPath;
    return { data: data, frames: [] };
}

/** Egg metadata read from the directory next to a structure file. */
interface EggMeta { eggFile?: string; unitFile?: string }

/** What the eggFile chooser needs to score each candidate against the structure. */
interface EggMatchContext {
    atoms: MolecularData['atoms'];
    bonds: MolecularData['bonds'];
    molecule: VmdMolecule;
    scene: VmdScene;
}

/**
 * Read the toolchain's egg metadata from the directory of a structure file: a
 * `unitfile` (every unit's atoms) and an `eggFile*` (the eggs and the unit
 * numbers they wrap). Either file is optional; anything unreadable or absent
 * yields undefined, which just leaves the surface patches unlabelled.
 * When several eggFile* candidates sit next to the structure, usually only one
 * of them describes it, so the user is asked which one to use — each candidate
 * is scored against the structure's surface segments first and the match count
 * is shown as the hint. Cancelling just leaves the patches unlabelled. The
 * (directory, structure) pair is cached by the caller so a multi-molecule
 * script asks once instead of once per molecule.
 */
async function readEggMetadata(
    sourceUri: vscode.Uri | undefined,
    cache?: Map<string, Promise<EggMeta | undefined>>,
    matchContext?: EggMatchContext
): Promise<EggMeta | undefined> {
    if (!sourceUri) return undefined;
    const dirPath = path.dirname(sourceUri.fsPath);
    if (cache) {
        const key = dirPath + '|' + sourceUri.fsPath;
        let pending = cache.get(key);
        if (!pending) {
            pending = readEggMetadataFromDir(dirPath, matchContext);
            cache.set(key, pending);
        }
        return pending;
    }
    return readEggMetadataFromDir(dirPath, matchContext);
}

async function readEggMetadataFromDir(dirPath: string, matchContext?: EggMatchContext): Promise<EggMeta | undefined> {
    try {
        const dir = vscode.Uri.file(dirPath);
        const entries = await vscode.workspace.fs.readDirectory(dir);
        let unitName: string | undefined;
        const eggNames: string[] = [];
        for (const [name, type] of entries) {
            if (type !== vscode.FileType.File) continue;
            if (/^unitfile$/i.test(name)) unitName = name;
            else if (/^eggfile/i.test(name)) eggNames.push(name);
        }
        if (!unitName && eggNames.length === 0) return undefined;
        const readText = async (name: string): Promise<string | undefined> => {
            try {
                return new TextDecoder().decode(await vscode.workspace.fs.readFile(vscode.Uri.joinPath(dir, name)));
            } catch {
                return undefined;
            }
        };
        let eggName: string | undefined;
        if (eggNames.length === 1) {
            eggName = eggNames[0];
        } else if (eggNames.length > 1) {
            eggNames.sort();
            const unitText = unitName ? await readText(unitName) : undefined;
            const units = unitText !== undefined ? parseUnitFile(unitText) : null;
            const segments = matchContext
                ? resolveVmdScene(matchContext.atoms, matchContext.bonds, matchContext.molecule, matchContext.scene).segments
                : null;
            const items: vscode.QuickPickItem[] = [];
            for (const name of eggNames) {
                let description = '';
                if (segments) {
                    const text = await readText(name);
                    const eggs = text !== undefined ? parseEggFile(text) : null;
                    const matched = eggs && units
                        ? matchEggsToSegments(segments, eggs, units).filter(function (e) { return e; }).length
                        : 0;
                    description = (units ? '' : 'no unitfile, ')
                        + matched + '/' + segments.length + ' patches match';
                }
                items.push({ label: name, description: description });
            }
            const picked = await vscode.window.showQuickPick(items, {
                placeHolder: 'Several eggFile* files sit next to ' + path.basename(dirPath)
                    + ' - pick the one describing this structure',
                canPickMany: false
            });
            eggName = picked ? picked.label : undefined;
        }
        const meta: EggMeta = {};
        if (eggName) meta.eggFile = await readText(eggName);
        if (unitName) meta.unitFile = await readText(unitName);
        return meta;
    } catch {
        return undefined;
    }
}

/**
 * Resolve one molecule of a VMD script into a drawable payload. The source file
 * path is relative to the script unless it is absolute, exactly as in VMD.
 */
function resolveVmdSource(scriptUri: vscode.Uri, sourceFile: string): vscode.Uri | undefined {
    if (!sourceFile) return undefined;
    if (path.isAbsolute(sourceFile)) return vscode.Uri.file(sourceFile);
    return vscode.Uri.file(path.resolve(path.dirname(scriptUri.fsPath), sourceFile));
}

/** Build the webview payload for one already-loaded VMD molecule. */function buildVmdMoleculePayload(
    molecule: VmdMolecule,
    data: MolecularData,
    scene: VmdScene,
    colorOf: (element: string) => string,
    eggMeta?: { eggFile?: string; unitFile?: string }
): VmdMoleculePayload {
    const resolved = resolveVmdScene(data.atoms, data.bonds, molecule, scene);

    // A script that draws eggs usually sits next to the toolchain's egg/unit
    // files; matching atom sets attaches each egg's name and unit numbers to
    // the patch that draws it. Files that are missing or malformed simply
    // leave the patches unlabelled.
    let eggs: (VmdEgg | null)[] | null = null;
    if (eggMeta && resolved.segments.length > 0) {
        eggs = matchEggsToSegments(resolved.segments,
            eggMeta.eggFile ? parseEggFile(eggMeta.eggFile) : null,
            eggMeta.unitFile ? parseUnitFile(eggMeta.unitFile) : null);
    }

    const atoms: VmdAtomPayload[] = data.atoms.map(function (a, i) {
        const item: VmdAtomPayload = {
            element: a.element, x: a.x, y: a.y, z: a.z,
            color: resolved.atomColors[i] || colorOf(a.element),
            vdw: vdwRadius(a.element)
        };
        if (typeof a.beta === 'number') item.beta = a.beta;
        if (typeof a.occupancy === 'number') item.occupancy = a.occupancy;
        if (typeof resolved.atomScales[i] === 'number') item.vmdScale = resolved.atomScales[i] as number;
        return item;
    });

    const bonds: VmdBondPayload[] = data.bonds.map(function (b, i) {
        const item: VmdBondPayload = { atom1: b.atom1, atom2: b.atom2, order: b.order };
        if (typeof resolved.bondRadii[i] === 'number') item.vmdR = resolved.bondRadii[i] as number;
        return item;
    });

    const segments: VmdSegmentPayload[] = resolved.segments.map(function (s, i) {
        const item: VmdSegmentPayload = { color: s.color, indices: s.indices.slice(), probe: s.probe };
        const egg = eggs ? eggs[i] : null;
        if (egg) {
            item.egg = { name: egg.name, units: egg.units.slice() };
            if (egg.level !== undefined) item.egg.level = egg.level;
        }
        return item;
    });

    const total = data.atoms.length;
    const reps: VmdRepPayload[] = molecule.reps.map(function (rep) {
        const selection = resolveSelection(rep.selection, data.atoms);
        const item: VmdRepPayload = {
            style: rep.style,
            colorMode: rep.colorMode,
            material: rep.material,
            selection: rep.selection,
            atoms: selection.indices === null ? total : selection.indices.length
        };
        if (typeof rep.colorArg === 'number') item.colorArg = rep.colorArg;
        if (typeof rep.scaleMin === 'number') item.scaleMin = rep.scaleMin;
        if (typeof rep.scaleMax === 'number') item.scaleMax = rep.scaleMax;
        return item;
    });

    return {
        label: molecule.label || molecule.sourceFile,
        sourceFile: data.filePath || molecule.sourceFile,
        atoms: atoms,
        bonds: bonds,
        segments: segments,
        reps: reps,
        warnings: resolved.warnings
    };
}

export class MolecularViewerProvider implements vscode.CustomReadonlyEditorProvider<MolecularDocument> {
    constructor(private readonly context: vscode.ExtensionContext) {}

    async openCustomDocument(
        uri: vscode.Uri,
        _openContext: vscode.CustomDocumentOpenContext,
        _token: vscode.CancellationToken
    ): Promise<MolecularDocument> {
        const content = await vscode.workspace.fs.readFile(uri);
        const textContent = new TextDecoder().decode(content);
        const fileName = uri.path.split('/').pop() || 'unknown.xyz';

        const ext = fileName.toLowerCase().split('.').pop() || '';
        let data: MolecularData;
        let frames: (LogFrame | OrcaFrame)[] = [];
        let atomGroups: AtomGroup[] | undefined;
        let optSteps: OptStep[] | undefined;
        let normalModes: NormalMode[] | undefined;
        let routes: RouteSection[] | undefined;
        let optThresholds: OptThresholds | undefined;
        let logSource: string | undefined;
        let vmdPayload: VmdPayload | undefined;
        let autoBeta = false;

        try {
        if (ext === 'tcl') {
            // A script is a display recipe: every `mol new` it performs becomes
            // one selectable molecule in the viewer, and the first one that
            // loads is also the document's structure so Save As, Diff and the
            // editing tools keep working on a real coordinate file.
            const tclResult = parseTcl(textContent);
            const loaded = await this.loadTclScene(tclResult.scene, uri);
            vmdPayload = loaded.payload && loaded.payload.molecules.length > 0 ? loaded.payload : undefined;
            data = loaded.data;
            frames = loaded.frames;
            optSteps = loaded.optSteps;
            normalModes = loaded.normalModes;
            routes = loaded.routes;
            logSource = loaded.logSource;
            if (!vmdPayload) {
                // Legacy path only: with a resolved scene the payload carries the
                // colours, and it is the one that decides them.
                atomGroups = tclResult.groups.length > 0 ? tclResult.groups : undefined;
                data.atomGroups = atomGroups;
            }
        } else if (ext === 'log' || ext === 'out') {
            const logResult = parseLogFile(textContent, fileName);
            frames = logResult.frames;
            optSteps = logResult.optSteps;
            normalModes = logResult.normalModes;
            routes = logResult.routes;
            optThresholds = logResult.optThresholds;
            logSource = logResult.source;
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
                data = { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false };
            }
        } else if (ext === 'xyz') {
            // A multi-frame XYZ is a trajectory — ORCA writes one for every
            // optimization as `<job>_trj.xyz`, with each frame's energy in its
            // comment line — so it gets the same frame navigation (and, when the
            // comment states an energy per frame, the same convergence chart) a
            // log does. A single-frame file takes the ordinary structure path.
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
                data = ensureBonds(parseFile(textContent, fileName));
            }
        } else {
            data = parseFile(textContent, fileName);
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
        } catch (err) {
            const msg = err instanceof Error ? err.message : String(err);
            data = { atoms: [], bonds: [], title: 'Parse error: ' + msg, hasExplicitBonds: false };
            frames = [];
            optSteps = undefined;
            normalModes = undefined;
            optThresholds = undefined;
        }

        if (ext !== 'tcl') {
            data.filePath = uri.fsPath;
        }

        return new MolecularDocument(uri, data, frames, optSteps, normalModes, logSource, routes, vmdPayload, autoBeta, optThresholds);
    }

    /**
     * Load every structure a VMD script refers to and resolve each one's
     * representations into drawable properties. Molecules that cannot be loaded
     * are skipped with a warning instead of failing the whole script, and the
     * first one that does load becomes the document's structure.
     */
    private async loadTclScene(scene: VmdScene, scriptUri: vscode.Uri): Promise<{
        payload: VmdPayload;
        data: MolecularData;
        frames: (LogFrame | OrcaFrame)[];
        optSteps?: OptStep[];
        normalModes?: NormalMode[];
        routes?: RouteSection[];
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
            const reference = scene.molecules.length > 0 ? scene.molecules[0].sourceFile : '';
            const title = reference
                ? 'Failed to load source file: ' + reference
                : 'No source file specified in TCL';
            return { payload: payload, data: { atoms: [], bonds: [], title: title, hasExplicitBonds: false }, frames: [] };
        }

        return {
            payload: payload,
            data: first.data,
            frames: first.frames,
            optSteps: first.optSteps,
            normalModes: first.normalModes,
            routes: first.routes,
            logSource: first.logSource
        };
    }

    async resolveCustomEditor(
        document: MolecularDocument,
        webviewPanel: vscode.WebviewPanel,
        _token: vscode.CancellationToken
    ): Promise<void> {
        webviewPanel.webview.options = {
            enableScripts: true,
        };

        webviewPanel.webview.html = await this.getHtmlForWebview(webviewPanel.webview, document.data, document.frames, document.optSteps, document.normalModes, document.logSource, document.routes, document.vmdPayload, document.autoBeta, document.optThresholds);

        // The panel's Molecular Mechanics run (at most one at a time). The pump
        // steps it in chunks and yields between them, so the panel keeps
        // receiving progress and its Pause / Resume / Cancel -- or a brand new
        // run -- are handled while the optimisation is still working.
        let mmRun: MmRunState | null = null;

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
            setTimeout(() => mmPump(run), 0);
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
            mmDiscardRun('replaced');
        });

        webviewPanel.webview.onDidReceiveMessage(async (message) => {
            switch (message.command) {
                case 'saveFile':
                    try {
                        const srcPath = message.filePath || '';
                        const srcDir = srcPath ? srcPath.substring(0, srcPath.replace(/\\/g, '/').lastIndexOf('/')) : '';
                        const defaultName = message.suggestedName || 'molecule.xyz';
                        const defaultUri = srcDir ? vscode.Uri.file(srcDir + '/' + defaultName) : vscode.Uri.file(defaultName);
                        const uri = await vscode.window.showSaveDialog({
                            defaultUri: defaultUri,
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
                                'All Files': ['*']
                            }
                        });
                        if (uri) {
                            await vscode.workspace.fs.writeFile(uri, new TextEncoder().encode(message.content));
                            vscode.window.showInformationMessage('Saved: ' + uri.fsPath);
                        }
                    } catch (e: any) {
                        vscode.window.showErrorMessage('Save failed: ' + (e.message || e));
                    }
                    break;
                case 'diffFile':
                    try {
                        const result = await vscode.window.showOpenDialog({
                            canSelectMany: false,
                            openLabel: 'Select File to Compare',
                            filters: {
                                'Molecular Files': ['gjf', 'xyz', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol2', 'log', 'out', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'tcl', 'cif', 'vasp', 'poscar', 'contcar', 'cube', 'vesta'],
                                'All Files': ['*']
                            }
                        });
                        if (!result || result.length === 0) {
                            webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                            break;
                        }
                        const diffUri = result[0];
                        const diffContent = await vscode.workspace.fs.readFile(diffUri);
                        const diffText = new TextDecoder().decode(diffContent);
                        const diffFileName = diffUri.path.split('/').pop() || 'unknown.xyz';
                        const diffExt = diffFileName.toLowerCase().split('.').pop() || '';

                        let diffData: MolecularData;
                        let diffFrames: (LogFrame | OrcaFrame)[] = [];
                        let diffVmdColors: string[] | null = null;

                        if (diffExt === 'tcl') {
                            const tclResult = parseTcl(diffText);
                            const loaded = await this.loadTclScene(tclResult.scene, diffUri);
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
                        } else if (diffExt === 'log' || diffExt === 'out') {
                            const logResult = parseLogFile(diffText, diffFileName);
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
                            diffData = parseFile(diffText, diffFileName);
                            if (diffData.crystal) {
                                diffData.bonds = [];
                                diffData.crystal.baseBonds = [];
                            } else {
                                diffData = ensureBonds(diffData);
                            }
                            diffData.filePath = diffUri.fsPath;
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
                            color: (diffVmdColors && diffVmdColors[i]) || diffAtomColors[a.element] || '#FF1493'
                        }));
                        if (diffData.atomGroups && diffData.atomGroups.length > 0) {
                            const override: { [key: number]: string } = {};
                            for (const g of diffData.atomGroups) for (const idx of g.indices) override[idx] = g.color;
                            for (let i = 0; i < diffAtoms.length; i++) if (override[i]) diffAtoms[i].color = override[i];
                        }
                        const diffBonds = diffData.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order }));
                        const diffFramesData = diffFrames.map(f => ({
                            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: diffAtomColors[a.element] || '#FF1493' })),
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
                            charge: diffData.charge
                        });
                    } catch (e: any) {
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
                                             impAtoms, impBonds, true, choice.copies);
                            break;
                        }

                        const result = await vscode.window.showOpenDialog({
                            canSelectMany: false,
                            openLabel: 'Select Structure to Import',
                            filters: {
                                'Molecular Files': ['gjf', 'xyz', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol2', 'log', 'out', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'cif', 'vasp', 'poscar', 'contcar', 'cube', 'vesta'],
                                'All Files': ['*']
                            }
                        });
                        if (!result || result.length === 0) {
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }
                        const impUri = result[0];
                        const impContent = await vscode.workspace.fs.readFile(impUri);
                        const impText = new TextDecoder().decode(impContent);
                        const impFileName = impUri.path.split('/').pop() || 'unknown.xyz';
                        const impExt = impFileName.toLowerCase().split('.').pop() || '';

                        if (impExt === 'tcl') {
                            vscode.window.showErrorMessage('Import: TCL scripts are not supported (use Diff for TCL).');
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }

                        let impData: MolecularData | undefined;
                        if (impExt === 'log' || impExt === 'out') {
                            const logResult = parseLogFile(impText, impFileName);
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
                            impData = parseFile(impText, impFileName);
                            // Crystal files (CIF/VASP/VESTA/Cube): import the
                            // unit-cell content as a molecule; periodicity is
                            // dropped (the current structure stays molecular).
                            if (impData.crystal) {
                                const cellAtoms = (impData.atoms.length > 0 ? impData.atoms : impData.crystal.baseAtoms)
                                    .map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, index: 0 }));
                                impData = ensureBonds({ atoms: cellAtoms, bonds: [], title: impData.title, hasExplicitBonds: false });
                            } else {
                                impData = ensureBonds(impData);
                            }
                        }
                        if (!impData || impData.atoms.length === 0) {
                            vscode.window.showErrorMessage('Import: no atoms found in ' + impFileName);
                            webviewPanel.webview.postMessage({ command: 'importResult', cancelled: true });
                            break;
                        }

                        const impAtoms = impData.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: ATOM_COLORS[a.element] || '#FF1493' }));
                        const impBonds = impData.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order }));
                        postImportCopies(webviewPanel.webview, impFileName,
                                         impAtoms, impBonds, !!impData.hasExplicitBonds, choice.copies);
                    } catch (e: any) {
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
                case 'error':
                    vscode.window.showErrorMessage(message.text);
                    break;
            }
        });
    }

    private async getHtmlForWebview(webview: vscode.Webview, data: MolecularData, frames: (LogFrame | OrcaFrame)[] = [], optSteps?: OptStep[], normalModes?: NormalMode[], logSource?: string, routes?: RouteSection[], vmdPayload?: VmdPayload, autoBeta: boolean = false, optThresholds?: OptThresholds): Promise<string> {
        const nonce = getNonce();

        const threeJsBytes = await vscode.workspace.fs.readFile(
            vscode.Uri.joinPath(this.context.extensionUri, 'media', 'three.min.js')
        );
        const threeJsContent = new TextDecoder().decode(threeJsBytes);

        const atomColors = ATOM_COLORS;

        const atomData = data.atoms.map(a => {
            let color = atomColors[a.element] || '#FF1493';
            return {
                element: a.element, x: a.x, y: a.y, z: a.z,
                color: color,
                beta: a.beta,
                occupancy: a.occupancy,
                baseIdx: a.baseIdx,
                cellI: a.cellI, cellJ: a.cellJ, cellK: a.cellK
            };
        });

        if (data.atomGroups && data.atomGroups.length > 0) {
            const atomColorOverride: { [key: number]: string } = {};
            for (const group of data.atomGroups) {
                for (const idx of group.indices) {
                    atomColorOverride[idx] = group.color;
                }
            }
            for (let i = 0; i < atomData.length; i++) {
                if (atomColorOverride[i]) {
                    atomData[i].color = atomColorOverride[i];
                }
            }
        }

        const bondData = data.bonds.map(b => ({
            atom1: b.atom1, atom2: b.atom2, order: b.order
        }));

        const framesData = frames.map(f => ({
            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: atomColors[a.element] || '#FF1493' })),
            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
            stepLabel: f.stepLabel,
            // Single-point energy of this frame's structure (opt logs only), used
            // by the convergence charts to highlight the selected step and to
            // report its exact value. [2026-09-12 | Atreides-Jimmy]
            energy: (f as { energy?: number }).energy
        }));

        const atomGroupsData = data.atomGroups ? data.atomGroups.map(g => ({
            colorId: g.colorId, color: g.color, indices: g.indices
        })) : [];

        const crystalData = data.crystal ? {
            a: data.crystal.a, b: data.crystal.b, c: data.crystal.c,
            alpha: data.crystal.alpha, beta: data.crystal.beta, gamma: data.crystal.gamma,
            latticeVectors: data.crystal.latticeVectors,
            spaceGroup: data.crystal.spaceGroup || '',
            symmetryOps: data.crystal.symmetryOps,
            baseAtoms: data.crystal.baseAtoms.map(a => ({
                element: a.element, x: a.x, y: a.y, z: a.z,
                index: a.index, occupancy: a.occupancy, baseIdx: a.baseIdx
            })),
            baseBonds: data.crystal.baseBonds.map(b => ({
                atom1: b.atom1, atom2: b.atom2, order: b.order
            }))
        } : null;

        const optStepsData = optSteps ? optSteps.map(s => ({
            step: s.step, energy: s.energy, maxForce: s.maxForce, rmsForce: s.rmsForce,
            maxDisplacement: s.maxDisplacement, rmsDisplacement: s.rmsDisplacement
        })) : null;

        const normalModesData = normalModes ? normalModes.map(m => ({
            index: m.index, frequency: m.frequency, symmetry: m.symmetry || '',
            reducedMass: m.reducedMass, forceConstant: m.forceConstant, irIntensity: m.irIntensity,
            displacements: m.displacements
        })) : null;

        const routesData = routes ? routes.map(r => ({
            raw: r.raw, hasOpt: r.hasOpt,
            keywords: r.keywords.map(k => ({ name: k.name, options: k.options }))
        })) : null;

        const vmdData = vmdPayload && vmdPayload.molecules.length > 0 ? vmdPayload : null;

        const optThresholdsData = optThresholds ? {
            maxForce: optThresholds.maxForce, rmsForce: optThresholds.rmsForce,
            maxDisplacement: optThresholds.maxDisplacement, rmsDisplacement: optThresholds.rmsDisplacement
        } : null;

        const jsonData = JSON.stringify({ atoms: atomData, bonds: bondData, title: data.title, atomColors: atomColors, filePath: data.filePath || '', frames: framesData, gjfMeta: data.gjfMeta || null, charge: data.charge, multiplicity: data.multiplicity, atomGroups: atomGroupsData, crystal: crystalData, optSteps: optStepsData, normalModes: normalModesData, optSource: logSource || null, optThresholds: optThresholdsData, routes: routesData, vmd: vmdData, autoBeta: !!(autoBeta && !vmdData), hasExplicitBonds: !!data.hasExplicitBonds, mmParams: describeForceField() });

        // The bond perception module, as source. Both sides run this text —
        // the host calls the compiled function directly, the webview evaluates
        // the serialisation below — so the viewport's re-detection can never
        // drift away from what the parser produced.
        const bondPerceptionSource = BOND_PERCEPTION_SOURCE;

        return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'nonce-${nonce}'; style-src 'nonce-${nonce}';">
<title>Molecular Viewer</title>
<style nonce="${nonce}">
*{margin:0;padding:0;box-sizing:border-box}
html{width:100%;height:100%;overflow:hidden}
body{width:100%;height:100%;overflow:hidden;display:flex;flex-direction:column;background:var(--vscode-editor-background,#1e1e1e);font-family:var(--vscode-font-family,sans-serif);color:var(--vscode-editor-foreground,#ccc)}
/* The viewer is not a document: nothing on the page is text to be selected, so
   Ctrl+A can only mean Select All of the *structure*, and a drag that starts on
   the canvas and ends over a panel never paints the toolbar blue on the way.
   This is a structural guarantee — it holds even where the key never reaches the
   page's own handler. user-select inherits, so one rule covers every panel;
   form fields get their text behaviour back, and the load-error block is
   deliberately selectable so a failure can be copied into a bug report. */
body{-webkit-user-select:none;user-select:none}
input,textarea{-webkit-user-select:text;user-select:text}
#error-msg{-webkit-user-select:text;user-select:text}
#toolbar{flex-shrink:0;background:var(--vscode-editor-background,#1e1e1e);border-bottom:1px solid var(--vscode-panel-border,#444);display:flex;flex-direction:column;padding:2px 0;z-index:20}
#toolbar .trow{display:flex;align-items:center;gap:2px;min-height:34px;padding:0 8px;overflow-x:auto;overflow-y:hidden;scrollbar-width:thin}
#toolbar .trow::-webkit-scrollbar{height:3px}
#toolbar .trow::-webkit-scrollbar-track{background:transparent}
#toolbar .trow::-webkit-scrollbar-thumb{background:rgba(255,255,255,0.15);border-radius:2px}
.tbtn{background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff);border:1px solid var(--vscode-panel-border,#444);padding:3px 10px;border-radius:3px;cursor:pointer;font-size:11px;white-space:nowrap;flex-shrink:0}
.tbtn:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
.tbtn.active{background:var(--vscode-button-background,#0e639c);border-color:var(--vscode-button-background,#0e639c)}
.tsep{width:1px;height:20px;background:var(--vscode-panel-border,#444);margin:0 4px;flex-shrink:0}
#status-bar{height:24px;flex-shrink:0;background:var(--vscode-statusBar-background,#007acc);color:var(--vscode-statusBar-foreground,#fff);display:flex;align-items:center;padding:0 10px;font-size:11px;z-index:20;gap:12px;overflow:hidden}
#selection-info{white-space:nowrap;overflow:hidden;text-overflow:ellipsis;flex:1;min-width:0}
#mode-info{white-space:nowrap;flex-shrink:0}
#container{flex:1;position:relative;overflow:hidden;min-height:0;outline:none}
#mol-info{position:absolute;top:8px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.55);padding:6px 10px;border-radius:4px;z-index:25;pointer-events:none;line-height:1.6;max-width:calc(100% - 16px);overflow:hidden}
#mm-panel{display:none;position:absolute;top:8px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.62);border:1px solid var(--vscode-panel-border,#444);border-radius:4px;padding:8px 10px;z-index:26;line-height:1.5;width:340px;max-width:calc(100% - 16px);max-height:calc(100% - 60px);overflow-y:auto;pointer-events:auto}
#mm-panel.show{display:block}
#mm-panel .mm-head{display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:5px;font-weight:600}
#mm-panel .mm-close{cursor:pointer;opacity:0.7;padding:0 2px}
#mm-panel .mm-close:hover{opacity:1}
#mm-panel .mm-param{color:var(--vscode-descriptionForeground,#a9a9a9);font-size:10px;line-height:1.45;margin:1px 0}
#mm-panel .mm-sec{margin:7px 0 3px;padding-top:6px;border-top:1px solid rgba(255,255,255,0.15);font-size:10px;color:var(--vscode-descriptionForeground,#999);text-transform:uppercase;letter-spacing:0.4px}
#mm-panel .mm-ctl{display:flex;align-items:center;gap:6px;margin:3px 0}
#mm-panel .mm-ctl-label{flex:0 0 96px;color:var(--vscode-descriptionForeground,#a9a9a9)}
#mm-panel .mm-ctl-sel{flex:1 1 auto;min-width:0;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;font-size:11px;padding:1px 2px}
#mm-panel .mm-sel{margin-top:6px;padding-top:5px;border-top:1px solid rgba(255,255,255,0.15);word-break:break-word}
#mm-panel .mm-run{margin-top:6px;width:100%;background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border:none;border-radius:3px;padding:4px 6px;font-size:11px;cursor:pointer}
#mm-panel .mm-run:hover:not(:disabled){background:var(--vscode-button-hoverBackground,#1177bb)}
#mm-panel .mm-run:disabled{opacity:0.5;cursor:default}
#mm-panel .mm-res{margin-top:6px;padding-top:5px;border-top:1px solid rgba(255,255,255,0.15);font-size:10px;line-height:1.45}
#mm-panel .mm-res-head{color:var(--vscode-textLink-foreground,#3794ff);font-weight:600}
#mm-panel .mm-err{color:var(--vscode-errorForeground,#f66);margin-top:5px;font-size:10px;line-height:1.4}
#mm-panel .mm-btns{display:flex;gap:6px;margin-top:6px}
#mm-panel .mm-btns button{flex:1 1 0;min-width:0;background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#ccc);border:1px solid var(--vscode-button-border,#555);border-radius:3px;padding:4px 6px;font-size:11px;cursor:pointer}
#mm-panel .mm-btns button:hover:not(:disabled){background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#mm-panel .mm-btns button:disabled{opacity:0.5;cursor:default}
#mm-panel .mm-prog{margin-top:6px;padding-top:5px;border-top:1px solid rgba(255,255,255,0.15);font-size:10px;line-height:1.45;white-space:pre-wrap}
#mm-panel .mm-curve-title{color:var(--vscode-descriptionForeground,#999);font-size:10px;margin:5px 0 1px}
#mm-panel .mm-curve{width:100%;height:78px;display:block}
#mm-reopen{display:none;position:absolute;top:36px;left:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#mm-reopen.show{display:block}
#mm-reopen:hover{background:var(--vscode-button-background,#0e639c)}
#mm-reopen.busy{border-color:var(--vscode-textLink-foreground,#3794ff);color:var(--vscode-textLink-foreground,#3794ff)}
#axes-indicator{position:absolute;bottom:12px;left:12px;width:90px;height:90px;z-index:25;pointer-events:none}
#axes-indicator svg{width:100%;height:100%}
#crystal-panel{position:absolute;top:8px;right:8px;background:rgba(0,0,0,0.6);border:1px solid var(--vscode-panel-border,#444);border-radius:4px;padding:8px 10px;font-size:11px;color:var(--vscode-editor-foreground,#ccc);z-index:25;display:none;min-width:180px;max-width:calc(100% - 16px)}
#crystal-panel h4{margin:0 0 6px 0;font-size:11px;font-weight:600}
#crystal-panel .bnd-row{display:flex;align-items:center;gap:4px;margin:3px 0}
#crystal-panel .bnd-row label{width:42px;text-align:right}
#crystal-panel .bnd-row input{width:48px;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 3px;font-size:11px;text-align:center}
#crystal-panel .bnd-btn{margin-top:6px;width:100%;background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border:none;border-radius:3px;padding:4px 6px;font-size:11px;cursor:pointer}
#crystal-panel .bnd-btn:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#egg-panel{position:absolute;top:8px;right:8px;background:rgba(0,0,0,0.6);border:1px solid var(--vscode-panel-border,#444);border-radius:4px;padding:8px 10px;font-size:11px;color:var(--vscode-editor-foreground,#ccc);z-index:25;display:none;min-width:200px;max-width:calc(100% - 16px);max-height:calc(100% - 70px);overflow-y:auto}
#egg-panel.show{display:block}
#egg-panel .ep-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:5px;font-weight:600}
#egg-panel .ep-head .ep-actions{display:flex;gap:4px;align-items:center}
#egg-panel .ep-head button{background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff);border:1px solid var(--vscode-panel-border,#444);border-radius:3px;padding:1px 6px;font-size:10px;cursor:pointer}
#egg-panel .ep-head button:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#egg-panel .ep-close{cursor:pointer;opacity:0.7;padding:0 2px}
#egg-panel .ep-close:hover{opacity:1}
#egg-panel .ep-row{display:flex;align-items:center;gap:6px;margin:3px 0;padding:3px 6px 3px 5px;cursor:pointer;white-space:nowrap;border-left:3px solid var(--vscode-descriptionForeground,#888);border-radius:2px}
#egg-panel .ep-row:hover{background:rgba(255,255,255,0.06)}
#egg-panel .ep-row input{margin:0;cursor:pointer}
#egg-panel .ep-row.ep-hidden .ep-name,#egg-panel .ep-row.ep-hidden .ep-units{opacity:0.45;text-decoration:line-through}
#egg-panel .ep-name{font-weight:600}
#egg-panel .ep-units{color:var(--vscode-descriptionForeground,#999)}
#select-panel{display:none;position:absolute;bottom:40px;left:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:280px;max-width:calc(100% - 16px);box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#select-panel .sp-title{font-size:12px;font-weight:bold;margin-bottom:4px}
#select-panel input[type=text]{width:100%;padding:6px 8px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:12px;box-sizing:border-box}
#select-panel .sp-btns{display:flex;gap:4px;margin-top:4px}
#select-panel .sp-btns button{flex:1}
#rotate-panel{display:none;position:absolute;bottom:40px;left:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:320px;max-width:calc(100% - 16px);box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#rotate-panel .rp-title{font-size:12px;font-weight:bold;margin-bottom:6px}
#rotate-panel .rp-field{margin:4px 0}
#rotate-panel .rp-field label{display:block;font-size:11px;color:var(--vscode-descriptionForeground,#999);margin-bottom:2px}
#rotate-panel input[type=text]{width:100%;padding:6px 8px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:12px;box-sizing:border-box}
#rotate-panel .rp-btns{display:flex;gap:4px;margin-top:6px}
#rotate-panel .rp-btns button{flex:1}
#rotate-panel .rp-slider-row{display:flex;align-items:center;gap:6px;margin:4px 0}
#rotate-panel .rp-slider-row label{font-size:11px;white-space:nowrap}
#rotate-panel .rp-slider-row input[type=range]{flex:1}
#rotate-panel .rp-angle-val{font-size:12px;font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff);min-width:46px;text-align:right}
#rotate-panel .rp-slider-row input[type=number]{width:58px;flex-shrink:0;padding:3px 4px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:11px;box-sizing:border-box}
#rotate-panel .rp-error{color:var(--vscode-errorForeground,#f66);font-size:11px;margin-top:4px;line-height:1.4}
#rotate-panel .rp-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:4px;line-height:1.4}
/* The View-mode controls are two cards, not one: Auto Rotate and Edit
   Coordinates each carry their own collision box in layoutPanels(), so a tall
   neighbour (the Diff results panel above all) can push them apart instead of
   one of them ending up underneath something else. */
#view-panel,#bake-panel{display:block;position:absolute;bottom:40px;right:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:210px;max-width:calc(100% - 16px);box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#view-panel .vp-title,#bake-panel .vp-title{font-size:12px;font-weight:bold;margin-bottom:6px}
#view-panel .tbtn,#bake-panel .tbtn{width:100%;text-align:center}
#view-panel .vp-hint,#bake-panel .vp-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:6px;line-height:1.4}
canvas{display:block}
#atom-tooltip{position:absolute;display:none;color:var(--vscode-editor-foreground,#ccc);font-size:12px;background:var(--vscode-editor-background,#1e1e1e);padding:4px 8px;border-radius:3px;border:1px solid var(--vscode-panel-border,#444);pointer-events:none;z-index:30}
#box-select-rect{position:fixed;display:none;border:1px dashed #ffeb3b;background:rgba(255,235,59,0.12);pointer-events:none;z-index:40}
#modal-overlay{position:fixed;top:0;left:0;right:0;bottom:0;background:rgba(0,0,0,0.5);z-index:100;display:none;align-items:center;justify-content:center}
#modal-overlay.show{display:flex}
#modal{background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-panel-border,#444);border-radius:6px;padding:16px 20px;min-width:min(320px,calc(100vw - 32px));max-width:min(420px,calc(100vw - 32px));box-shadow:0 8px 30px rgba(0,0,0,0.5)}
#modal h3{font-size:14px;margin-bottom:10px;color:var(--vscode-editor-foreground,#ccc)}
#modal label{font-size:12px;display:block;margin:6px 0 2px}
#modal input[type=number],#modal select{width:100%;padding:4px 8px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:12px}
#modal input[type=range]{width:100%;margin:4px 0}
#modal .modal-row{display:flex;gap:8px;align-items:center;margin:4px 0}
#modal .modal-row label{margin:0;white-space:nowrap;min-width:60px}
#modal .modal-row input,#modal .modal-row select{flex:1}
#modal .modal-btns{display:flex;justify-content:flex-end;gap:8px;margin-top:14px}
#modal .mbtn{padding:5px 16px;border-radius:3px;cursor:pointer;font-size:12px;border:1px solid var(--vscode-panel-border,#444)}
#modal .mbtn-ok{background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border-color:var(--vscode-button-background,#0e639c)}
#modal .mbtn-ok:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#modal .mbtn-cancel{background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff)}
#modal .mbtn-cancel:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#modal .mbtn-danger{background:#c33;border-color:#c33;color:#fff}
#modal .current-val{font-size:13px;color:var(--vscode-descriptionForeground,#999);margin-bottom:6px}
#loading{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--vscode-editor-foreground,#ccc);font-size:14px}
#error-msg{display:none;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:#f66;padding:20px;font-size:13px;text-align:center;max-width:80%;z-index:50}
.hidden{display:none!important}
#frame-nav{display:none;align-items:center;gap:2px;flex-shrink:0}
#frame-nav.show{display:flex}
#frame-num{width:40px;padding:2px 4px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:11px;text-align:center}
#frame-info{color:var(--vscode-statusBar-foreground,#fff);font-size:11px;padding:0 4px;white-space:nowrap;max-width:160px;overflow:hidden;text-overflow:ellipsis}
#auto-play.playing{background:var(--vscode-button-background,#0e639c);border-color:var(--vscode-button-background,#0e639c)}
#mol-nav{display:none;align-items:center;gap:2px;flex-shrink:0}
#mol-nav.show{display:flex}
#mol-label{color:var(--vscode-statusBar-foreground,#fff);font-size:11px;padding:0 4px;white-space:nowrap;max-width:180px;overflow:hidden;text-overflow:ellipsis}
#beta-btn.active{background:var(--vscode-button-background,#0e639c);border-color:var(--vscode-button-background,#0e639c)}
#diff-panel{display:none;position:absolute;bottom:40px;right:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.65);padding:8px 12px;border-radius:4px;z-index:25;max-width:40%;max-height:60%;overflow-y:auto;pointer-events:auto;line-height:1.5}
#diff-panel.show{display:block}
#diff-panel h4{font-size:12px;margin-bottom:4px;color:var(--vscode-textLink-foreground,#3794ff)}
#diff-panel .diff-row{margin:2px 0}
#diff-panel .diff-close{position:absolute;top:4px;right:8px;cursor:pointer;color:var(--vscode-descriptionForeground,#999);font-size:14px}
#diff-panel .diff-close:hover{color:var(--vscode-errorForeground,#f66)}
#diff-reopen{display:none;position:absolute;bottom:40px;right:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#diff-reopen.show{display:block}
#diff-reopen:hover{background:var(--vscode-button-background,#0e639c)}
#diff-label,#diff-label-right{display:none;position:absolute;top:8px;color:var(--vscode-editor-foreground,#ccc);font-size:13px;font-weight:bold;background:rgba(0,0,0,0.55);padding:4px 10px;border-radius:4px;z-index:25;pointer-events:none;max-width:40%;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
#diff-label.show,#diff-label-right.show{display:block}
#diff-label{left:25%;transform:translateX(-50%)}
#diff-label-right{left:75%;transform:translateX(-50%)}
#diff-panel .diff-stats{margin-top:8px;padding-top:6px;border-top:1px solid rgba(255,255,255,0.18)}
#diff-panel .diff-stats-head{cursor:pointer;font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff);user-select:none;-webkit-user-select:none}
#diff-panel .diff-stats-head:hover{text-decoration:underline}
#diff-panel .diff-stats-arrow{display:inline-block;width:11px}
#diff-panel .diff-stats-body{margin-top:5px}
#diff-panel .diff-stats.collapsed .diff-stats-body{display:none}
#diff-panel .diff-stat-note{font-size:10px;line-height:1.35;color:var(--vscode-descriptionForeground,#999);margin-bottom:5px}
#diff-panel .diff-stat-title{font-size:11px;color:#f0c674;margin:7px 0 2px}
#diff-panel .diff-stat-scroll{overflow-x:auto;overflow-y:hidden}
#diff-panel .diff-stat-table{border-collapse:collapse;font-size:10px;white-space:nowrap;min-width:100%}
#diff-panel .diff-stat-table th,#diff-panel .diff-stat-table td{padding:1px 7px;text-align:right;font-weight:normal}
#diff-panel .diff-stat-table th:first-child,#diff-panel .diff-stat-table td:first-child{text-align:left}
#diff-panel .diff-stat-table th{color:var(--vscode-descriptionForeground,#999);border-bottom:1px solid rgba(255,255,255,0.18)}
#diff-panel .diff-stat-table tr.diff-stat-total td{border-top:1px solid rgba(255,255,255,0.18);font-weight:bold}
#opt-panel{display:none;position:absolute;top:36px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.78);padding:8px 10px;border-radius:4px;z-index:25;width:340px;max-width:calc(100% - 16px);max-height:70%;overflow:hidden;pointer-events:auto;line-height:1.4}
#opt-panel.show{display:flex;flex-direction:column}
#opt-panel .opt-head{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;padding-bottom:4px;border-bottom:1px solid rgba(255,255,255,0.15)}
#opt-panel .opt-title{font-size:12px;font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff)}
#opt-panel .opt-close{cursor:pointer;color:var(--vscode-descriptionForeground,#999);font-size:14px;padding:0 4px}
#opt-panel .opt-close:hover{color:var(--vscode-errorForeground,#f66)}
#opt-panel .opt-resize{height:6px;cursor:ns-resize;background:rgba(255,255,255,0.1);border-radius:3px;margin:4px 0}
#opt-panel .opt-resize:hover{background:rgba(255,255,255,0.25)}
#opt-panel .opt-body{overflow-y:scroll;overflow-x:hidden;flex:1 1 auto;min-height:0;width:100%;box-sizing:border-box}
#opt-panel canvas{display:block;width:100%;max-width:100%;margin:2px 0 6px}
#opt-panel .opt-chart-title{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:4px}
#opt-reopen{display:none;position:absolute;top:36px;left:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#opt-reopen.show{display:block}
#opt-reopen:hover{background:var(--vscode-button-background,#0e639c)}
#freq-panel{display:none;position:absolute;top:36px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.78);padding:8px 10px;border-radius:4px;z-index:25;width:300px;max-width:calc(100% - 16px);max-height:70%;overflow:hidden;pointer-events:auto;line-height:1.4}
#freq-panel.show{display:flex;flex-direction:column}
#freq-panel .freq-head{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;padding-bottom:4px;border-bottom:1px solid rgba(255,255,255,0.15)}
#freq-panel .freq-title{font-size:12px;font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff)}
#freq-panel .freq-close{cursor:pointer;color:var(--vscode-descriptionForeground,#999);font-size:14px;padding:0 4px}
#freq-panel .freq-close:hover{color:var(--vscode-errorForeground,#f66)}
#freq-panel .freq-body{overflow-y:scroll;overflow-x:hidden;flex:1 1 auto;min-height:0;width:100%;box-sizing:border-box}
#freq-panel .freq-row{display:flex;justify-content:space-between;align-items:center;padding:3px 4px;cursor:pointer;border-radius:3px}
#freq-panel .freq-row:hover{background:rgba(255,255,255,0.1)}
#freq-panel .freq-row.playing{background:rgba(255,200,0,0.25);color:#ffd700}
#freq-panel .freq-idx{width:28px;color:var(--vscode-descriptionForeground,#999)}
#freq-panel .freq-val{flex:1;font-family:monospace}
#freq-panel .freq-sym{width:30px;text-align:right;color:var(--vscode-descriptionForeground,#999)}
/* A mode has two independent affordances: the circle selects it for the
   manual-displacement slider, the ▶ plays it. They are separate controls on
   purpose — watching a mode must never arm the slider, and picking one must
   never start a playback. */
#freq-panel .freq-pick{width:11px;height:11px;flex:0 0 11px;margin:0 6px 0 2px;border-radius:50%;border:1px solid #767676;cursor:pointer;box-sizing:border-box;align-self:center}
#freq-panel .freq-pick:hover{border-color:#ffb86b;box-shadow:0 0 0 2px rgba(255,159,67,0.18)}
#freq-panel .freq-pick.on{border-color:#ff9f43;background:#ff9f43;box-shadow:inset 0 0 0 2px rgba(18,18,18,0.85)}
#freq-panel .freq-play{width:18px;text-align:center}
#freq-panel .freq-stop{margin-top:6px;width:100%;padding:4px;border-radius:3px;border:1px solid var(--vscode-button-border,#555);background:var(--vscode-button-background,#0e639c);color:#fff;cursor:pointer;font-size:11px}
#freq-panel .freq-stop:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
/* Manual displacement: one mode written into the coordinates. The block is a
   fixed-height footer of the panel (the mode list scrolls above it), so moving
   the slider cannot resize the panel and make it fight the layout. */
#freq-panel .freq-dis{margin-top:6px;padding-top:6px;border-top:1px solid rgba(255,255,255,0.15);flex:0 0 auto}
#freq-panel .vd-head{display:flex;justify-content:space-between;align-items:baseline;gap:6px;font-size:10px;letter-spacing:0.4px;text-transform:uppercase;color:var(--vscode-descriptionForeground,#999)}
#freq-panel .vd-mode{font-family:monospace;text-transform:none;letter-spacing:0;color:#ffb86b;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:168px}
/* The traceability line: every mode that is in the distortion, and the largest
   displacement of the whole combination. */
#freq-panel .vd-comps{margin-top:2px;font-family:monospace;font-size:10px;color:#c99a5b;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
#freq-panel .vd-slider{width:100%;margin:5px 0 3px;accent-color:#ff9f43}
#freq-panel .vd-slider:disabled{opacity:0.35}
#freq-panel .vd-row{display:flex;align-items:baseline;gap:6px}
#freq-panel .vd-val{flex:0 0 auto;font-family:monospace;font-size:11px;color:#ffb86b}
#freq-panel .vd-max{flex:1 1 auto;font-size:10px;color:var(--vscode-descriptionForeground,#999);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
#freq-panel .vd-reset{flex:0 0 auto;padding:2px 9px;border-radius:3px;border:1px solid var(--vscode-button-border,#555);background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#ccc);font-size:10px;cursor:pointer}
#freq-panel .vd-reset:hover:not(:disabled){background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#freq-panel .vd-reset:disabled{opacity:0.4;cursor:default}
#freq-panel .vd-hint{margin-top:3px;font-size:10px;line-height:1.35;color:var(--vscode-descriptionForeground,#999);min-height:27px}
/* A mode row carries a second, independent state: playing (the animation that
   borrows the geometry) and displaced (the session that writes it). */
#freq-panel .freq-row.displaced{box-shadow:inset 2px 0 0 #ff9f43}
#freq-panel .freq-row.displaced .freq-idx{color:#ffb86b}
#freq-reopen{display:none;position:absolute;top:36px;left:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#freq-reopen.show{display:block}
#freq-reopen:hover{background:var(--vscode-button-background,#0e639c)}
#route-panel{display:none;position:absolute;top:36px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.78);padding:8px 10px;border-radius:4px;z-index:25;width:360px;max-width:calc(100% - 16px);max-height:70%;overflow:hidden;pointer-events:auto;line-height:1.4}
#route-panel.show{display:flex;flex-direction:column}
#route-panel .rt-head{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;padding-bottom:4px;border-bottom:1px solid rgba(255,255,255,0.15)}
#route-panel .rt-title{font-size:12px;font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff)}
#route-panel .rt-close{cursor:pointer;color:var(--vscode-descriptionForeground,#999);font-size:14px;padding:0 4px}
#route-panel .rt-close:hover{color:var(--vscode-errorForeground,#f66)}
#route-panel .rt-body{overflow-y:scroll;overflow-x:hidden;flex:1 1 auto;min-height:0;width:100%;box-sizing:border-box}
#route-panel .rt-job{font-weight:bold;color:var(--vscode-textLink-foreground,#3794ff);margin:6px 0 3px}
#route-panel .rt-badge{display:inline-block;background:rgba(80,200,120,0.25);color:#69db7c;border-radius:3px;padding:0 5px;font-size:10px;font-weight:normal;margin-left:5px}
#route-panel .rt-raw{font-family:monospace;font-size:10px;background:rgba(255,255,255,0.06);border-radius:3px;padding:4px 6px;margin:3px 0 6px;word-break:break-all;color:var(--vscode-descriptionForeground,#aaa)}
#route-panel .rt-kw{display:flex;align-items:baseline;padding:1px 2px}
#route-panel .rt-kname{font-family:monospace;font-weight:bold;color:#ffd700;min-width:60px;flex:none}
#route-panel .rt-kopts{font-family:monospace;color:var(--vscode-descriptionForeground,#aaa);word-break:break-all}
#route-reopen{display:none;position:absolute;top:36px;left:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#route-reopen.show{display:block}
#route-reopen:hover{background:var(--vscode-button-background,#0e639c)}
@media (max-width:600px){.tbtn{padding:3px 6px;font-size:10px}#frame-info{font-size:10px}#toolbar .trow{gap:1px}}
</style>
</head>
<body>
<div id="toolbar">
<div class="trow" id="toolbar-modes">
<button class="tbtn active" data-mode="view">View</button>
<button class="tbtn" data-mode="bondLength">Bond Length</button>
<button class="tbtn" data-mode="bondAngle">Bond Angle</button>
<button class="tbtn" data-mode="dihedral">Dihedral</button>
<button class="tbtn" data-mode="rotateGroup">Rotate Group</button>
<button class="tbtn" data-mode="moveAtoms">Move Atoms</button>
<button class="tbtn" data-mode="mmOpt">MM Optimize</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="addAtom">Add Atom</button>
<button class="tbtn" data-mode="deleteAtom">Delete Atom</button>
<button class="tbtn" data-mode="replaceAtom">Replace Atom</button>
<button class="tbtn" data-mode="addH">Add H</button>
<button class="tbtn" id="import-btn">Import Structure</button>
</div>
<div class="trow" id="toolbar-actions">
<button class="tbtn" data-mode="selectAtoms">Select Atoms</button>
<button class="tbtn" data-mode="boxSelect">Box Select</button>
<div class="tsep"></div>
<button class="tbtn" id="beta-btn" style="display:none">Color by Beta</button>
<button class="tbtn" id="egg-btn" style="display:none">Eggs</button>
<button class="tbtn" id="diff-btn">Diff</button>
<button class="tbtn" id="undo-btn" disabled>Undo</button>
<button class="tbtn" id="save-btn">Save As</button>
<button class="tbtn" id="reset-btn">Reset View</button>
<div class="tsep" id="frame-sep"></div>
<div id="frame-nav">
<button class="tbtn" id="prev-frame">◀</button>
<input type="number" id="frame-num" min="1" value="1">
<span id="frame-info">1/1</span>
<button class="tbtn" id="next-frame">▶</button>
<button class="tbtn" id="auto-play">⏵ Play</button>
</div>
<div class="tsep hidden" id="mol-sep"></div>
<div id="mol-nav">
<button class="tbtn" id="prev-mol">◀</button>
<span id="mol-label">1/1</span>
<button class="tbtn" id="next-mol">▶</button>
</div>
</div>
</div>
<div id="status-bar"><span id="mode-info">View Mode</span><span id="selection-info"></span></div>
<div id="container"><div id="loading">Loading 3D Viewer...</div><div id="mol-info"></div><div id="mm-panel"></div><div id="mm-reopen">MM Optimize</div><div id="egg-panel"></div><div id="diff-label"></div><div id="diff-label-right"></div><div id="diff-panel"></div><div id="diff-reopen">📊 Show Results</div><div id="axes-indicator"></div><div id="crystal-panel"><h4>Supercell Bounds</h4><div class="bnd-row"><label>a min</label><input type="number" id="bnd-a-min" value="0" step="0.1"></div><div class="bnd-row"><label>a max</label><input type="number" id="bnd-a-max" value="1" step="0.1"></div><div class="bnd-row"><label>b min</label><input type="number" id="bnd-b-min" value="0" step="0.1"></div><div class="bnd-row"><label>b max</label><input type="number" id="bnd-b-max" value="1" step="0.1"></div><div class="bnd-row"><label>c min</label><input type="number" id="bnd-c-min" value="0" step="0.1"></div><div class="bnd-row"><label>c max</label><input type="number" id="bnd-c-max" value="1" step="0.1"></div><button id="bnd-remove-disorder" class="bnd-btn">Remove Disorder &lt;0.5</button></div><div id="view-panel"><div class="vp-title">Auto Rotate</div><button class="tbtn" id="spin-btn">Off</button><div class="vp-hint">Drag and release to spin, or hold an arrow key &mdash; the longer the hold, the faster. The next left press stops it.</div></div><div id="bake-panel"><div class="vp-title">Edit Coordinates</div><div class="vp-hint" id="vp-bake-hint">Hold <b>Alt</b> and drag: left rotates the structure, right/middle moves it. The atom coordinates change, so Save writes the new geometry.</div></div><div id="select-panel"><div class="sp-title">Select Atoms</div><input type="text" id="sel-panel-input" placeholder="e.g. 1 3-5 C H 8"><div class="sp-btns"><button class="mbtn mbtn-ok" id="sel-panel-ok">Select</button><button class="mbtn mbtn-cancel" id="sel-panel-clear">Clear</button></div></div><div id="rotate-panel"><div class="rp-title">Rotate Group Around Axis</div><div class="rp-field"><label>Group atoms (click atoms or type, e.g. H581 O380 C44 ...):</label><input type="text" id="rp-group-input" placeholder="e.g. H581 O380 C44 C43 C42 ..."></div><div class="rp-field"><label>Axis atoms (2+ atoms, single-bonded chain):</label><input type="text" id="rp-axis-input" placeholder="e.g. C32 C33 C43 C53  or  C33 C43"></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-apply">Apply &amp; Rotate</button><button class="mbtn mbtn-cancel" id="rp-clear">Clear</button></div><div class="rp-error" id="rp-error"></div><div id="rp-slider-wrap" style="display:none"><div class="rp-slider-row"><label>Angle:</label><input type="range" id="rp-slider" min="0" max="360" step="0.5" value="0"><input type="number" id="rp-angle-input" min="0" max="360" step="0.5" value="0" title="Type an angle (0-360) and press Enter"><span class="rp-angle-val" id="rp-angle-val">0°</span></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-done">Done</button><button class="mbtn mbtn-danger" id="rp-reset">Reset</button></div></div><div class="rp-hint">Axis atoms are kept fixed; the rest of the selected group rotates around the line through the axis atoms (0–360°).</div></div><div id="opt-panel"></div><div id="opt-reopen">📈 Convergence</div><div id="freq-panel"></div><div id="freq-reopen">🎵 Modes</div><div id="route-panel"></div><div id="route-reopen">⚙️ Route</div></div>
<div id="error-msg"></div>
<div id="atom-tooltip"></div>
<div id="box-select-rect"></div>
<div id="modal-overlay"><div id="modal"></div></div>
<script nonce="${nonce}">
${threeJsContent}
</script>
<script nonce="${nonce}">
(function(){
try{
var MD=${jsonData};
var BP=${bondPerceptionSource};
var AN={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Ne:10,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Sc:21,Ti:22,V:23,Cr:24,Mn:25,Fe:26,Co:27,Ni:28,Cu:29,Zn:30,Ga:31,Ge:32,As:33,Se:34,Br:35,Kr:36,Rb:37,Sr:38,Y:39,Zr:40,Nb:41,Mo:42,Tc:43,Ru:44,Rh:45,Pd:46,Ag:47,Cd:48,In:49,Sn:50,Sb:51,Te:52,I:53,Xe:54,Cs:55,Ba:56,La:57,Ce:58,Pr:59,Nd:60,Pm:61,Sm:62,Eu:63,Gd:64,Tb:65,Dy:66,Ho:67,Er:68,Tm:69,Yb:70,Lu:71,Hf:72,Ta:73,W:74,Re:75,Os:76,Ir:77,Pt:78,Au:79,Hg:80,Tl:81,Pb:82,Bi:83,Po:84,At:85,Rn:86};
var CRY=MD.crystal||null;
// ===== VMD visualization script =====
// Present when the document is a .tcl script: every molecule the script
// creates, with the drawing properties its representations resolve to. The
// script's first molecule becomes the structure on screen (so editing, Save As
// and Diff operate on a real coordinate file); the toolbar's molecule selector
// switches between them.
var VMD=MD.vmd||null;
var vmdMolIdx=0;
var vmdSurfaceGroup=null;
var vmdSurfaceMeshes=[];
// Which Surf patches the egg panel has hidden, keyed by the patch's segment
// index (the patch order of the script). Declared here, above the first
// rebuildScene: a var further down would still be undefined there.
var vmdHiddenSegments={};
if(VMD){
    VMD.activeSegments=(VMD.molecules[0]&&VMD.molecules[0].segments)||[];
    MD.atoms=VMD.molecules[0].atoms;
    MD.bonds=VMD.molecules[0].bonds;
    MD.atoms.forEach(function(a,i){a.index=i});
    MD.bonds.forEach(function(b,i){b.index=i});
}
var boundary={aMin:0,aMax:1,bMin:0,bMax:1,cMin:0,cMax:1};
var occTextureCache={};
var CRY_INV=null;
function computeCrystalInv(){
    if(!CRY){CRY_INV=null;return}
    var lv=CRY.latticeVectors;
    var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
    if(Math.abs(det)<1e-12){CRY_INV=null;return}
    var invDet=1/det;
    CRY_INV=[
        [(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])*invDet,(lv[1][2]*lv[2][0]-lv[1][0]*lv[2][2])*invDet,(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0])*invDet],
        [(lv[0][2]*lv[2][1]-lv[0][1]*lv[2][2])*invDet,(lv[0][0]*lv[2][2]-lv[0][2]*lv[2][0])*invDet,(lv[0][1]*lv[2][0]-lv[0][0]*lv[2][1])*invDet],
        [(lv[0][1]*lv[1][2]-lv[0][2]*lv[1][1])*invDet,(lv[0][2]*lv[1][0]-lv[0][0]*lv[1][2])*invDet,(lv[0][0]*lv[1][1]-lv[0][1]*lv[1][0])*invDet]
    ];
}
computeCrystalInv();
function formulaFromAtoms(arr){
    var counts={};
    arr.forEach(function(a){counts[a.element]=(counts[a.element]||0)+1});
    var keys=Object.keys(counts);
    var hasC=counts.C!==undefined;
    keys.sort(function(a,b){
        if(hasC){
            if(a==='C'&&b!=='C')return-1;
            if(b==='C'&&a!=='C')return 1;
            if(a==='H'&&b!=='H')return-1;
            if(b==='H'&&a!=='H')return 1;
        }
        return a<b?-1:1;
    });
    var s='';
    keys.forEach(function(k){s+=k+(counts[k]>1?counts[k]:'')});
    return s;
}
function updateMolInfo(){
    var infoEl=document.getElementById('mol-info');
    if(!infoEl)return;
    var nAtoms=MD.atoms.length;
    var chrg=MD.charge!=null?MD.charge:'-';
    var mult=MD.multiplicity!=null?MD.multiplicity:'-';
    var nElectrons=0;
    MD.atoms.forEach(function(a){nElectrons+=(AN[a.element]||0)});
    if(typeof chrg==='number')nElectrons-=chrg;
    var html='Formula: '+formulaFromAtoms(MD.atoms)+'<br>Atoms: '+nAtoms+'<br>Charge: '+chrg+'<br>Electrons: '+nElectrons+'<br>Multiplicity: '+mult;
    if(CRY){
        var sa=boundary.aMax-boundary.aMin,sb=boundary.bMax-boundary.bMin,sc=boundary.cMax-boundary.cMin;
        html+='<br>Unit Cell: '+formulaFromAtoms(CRY.baseAtoms);
        html+='<br>Cell: '+CRY.a.toFixed(3)+' × '+CRY.b.toFixed(3)+' × '+CRY.c.toFixed(3);
        html+='<br>Angles: '+CRY.alpha.toFixed(1)+'°, '+CRY.beta.toFixed(1)+'°, '+CRY.gamma.toFixed(1)+'°';
        if(CRY.spaceGroup)html+='<br>SG: '+CRY.spaceGroup;
        html+='<br>Supercell: '+sa.toFixed(2)+'×'+sb.toFixed(2)+'×'+sc.toFixed(2);
    }
    if(VMD&&VMD.molecules[vmdMolIdx]){
        var vm=VMD.molecules[vmdMolIdx];
        var styles={};
        vm.reps.forEach(function(r){styles[r.style]=(styles[r.style]||0)+1});
        var styleText=Object.keys(styles).map(function(k){return k+'×'+styles[k]}).join(', ');
        html+='<br>VMD script: '+vm.reps.length+' rep'+(vm.reps.length===1?'':'s')+(styleText?' ('+styleText+')':'');
        html+='<br>Colour: '+VMD.colorScaleMethod;
        if(vm.segments.length)html+='<br>Surfaces: '+vm.segments.length;
        var warns=(VMD.warnings||[]).concat(vm.warnings||[]);
        if(warns.length)html+='<br><span style="color:#e8a33d">⚠ '+warns.length+' script warning'+(warns.length===1?'':'s')+'</span>';
    }
    infoEl.innerHTML=html;
}
updateMolInfo();
var container=document.getElementById('container');
var loadingEl=document.getElementById('loading');
var errorEl=document.getElementById('error-msg');
var tooltipEl=document.getElementById('atom-tooltip');
var boxRectEl=document.getElementById('box-select-rect');
var modeInfoEl=document.getElementById('mode-info');
var selInfoEl=document.getElementById('selection-info');
var modalOverlay=document.getElementById('modal-overlay');
var modalEl=document.getElementById('modal');
var vscodeApi=acquireVsCodeApi();

function showError(msg){if(loadingEl)loadingEl.style.display='none';errorEl.style.display='block';errorEl.textContent=msg}

if(typeof THREE==='undefined'){showError('Three.js library failed to load. Please reinstall the extension.');return}

var cw=container.clientWidth||window.innerWidth;
var ch=container.clientHeight||(window.innerHeight-60);
if(ch<1)ch=window.innerHeight-60;
if(cw<1)cw=window.innerWidth;

var scene=new THREE.Scene();
// The viewer keeps its own background: a script's "color Display Background"
// is a rendering setting of the reference renderer, and a white background here
// would clash with the editor theme and wash out white atoms.
scene.background=new THREE.Color(0x1e1e1e);
var camera=new THREE.PerspectiveCamera(60,cw/ch,0.1,1000);
var renderer=null;
try{
    renderer=new THREE.WebGLRenderer({antialias:true});
}catch(webglErr){
    // three.js throws here only when neither webgl2, webgl nor
    // experimental-webgl could be created at all — which in an editor webview
    // means the *window's* shared GPU process is unavailable (a driver reset,
    // or too many live WebGL contexts), not that this structure is too large.
    // No page can bring that process back, so the message says what can.
    throw new Error((webglErr&&webglErr.message?webglErr.message:'WebGL is unavailable')+
        ' - the GPU process of the editor window is unavailable. Run "Developer: Reload Window" (or restart the editor) and open the file again.');
}
renderer.setSize(cw,ch);
renderer.setPixelRatio(window.devicePixelRatio);
container.appendChild(renderer.domElement);
// A lost context is the other half of the failure: the driver reset takes the
// context away from every webview in the window at once, and a page that does
// not handle it simply goes blank. Keep the canvas alive, say what happened,
// and rebuild the scene if (and when) the context comes back.
renderer.domElement.addEventListener('webglcontextlost',function(ev){
    ev.preventDefault();
    contextLost=true;
    needsRender=false;
    if(loadingEl)loadingEl.style.display='none';
    if(errorEl){
        errorEl.style.display='block';
        errorEl.textContent='WebGL context lost - the graphics driver reset. Waiting for it to come back; if it does not, run "Developer: Reload Window".';
    }
},false);
renderer.domElement.addEventListener('webglcontextrestored',function(){
    contextLost=false;
    if(errorEl)errorEl.style.display='none';
    // three.js re-initialises its own GPU state, and the scene is rebuilt so
    // every buffer is uploaded fresh instead of from a cache that survived the
    // reset. The shared pools survive too: their data is still valid.
    rebuildScene();
    needsRender=true;
},false);
if(loadingEl)loadingEl.style.display='none';

scene.add(new THREE.AmbientLight(0x404040,1.5));
var dl1=new THREE.DirectionalLight(0xffffff,0.8);dl1.position.set(5,10,7);scene.add(dl1);
var dl2=new THREE.DirectionalLight(0xffffff,0.4);dl2.position.set(-5,-3,-5);scene.add(dl2);

var pivotGroup=new THREE.Group();scene.add(pivotGroup);
var moleculeGroup=new THREE.Group();pivotGroup.add(moleculeGroup);
// Surface patches of a VMD script live *inside* the molecule group, so they
// inherit every view transform the atoms get (rotation, panning) instead of
// standing still while the structure turns. rebuildScene keeps this group
// across its clear and rebuilds the meshes inside it.
vmdSurfaceGroup=new THREE.Group();
moleculeGroup.add(vmdSurfaceGroup);

var CX=0,CY=0,CZ=0;
MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length;

var CR={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.84,C:0.76,N:0.71,O:0.66,F:0.57,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,K:2.03,Ca:1.76,Fe:1.32,Cu:1.32,Zn:1.22,Br:1.20,I:1.39};
function getR(el){return(CR[el]||1.50)*0.5}
// Sphere radius of one atom. A VMD script scales the viewer's own radius (so a
// structure always keeps its ball-and-stick look, and 0 means "not drawn");
// without a script the plain default applies.
function atomRadius(a){
    var r=getR(a.element);
    return typeof a.vmdScale==='number'?r*a.vmdScale:r;
}

function getOccupancyTexture(color,occupancy){
    var key=color+'_'+occupancy.toFixed(3);
    if(occTextureCache[key])return occTextureCache[key];
    var canvas=document.createElement('canvas');
    canvas.width=4;canvas.height=256;
    var ctx=canvas.getContext('2d');
    var boundY=Math.round((1-occupancy)*255);
    ctx.fillStyle='#FFFFFF';
    ctx.fillRect(0,0,4,boundY);
    ctx.fillStyle=color;
    ctx.fillRect(0,boundY,4,256-boundY);
    var tex=new THREE.CanvasTexture(canvas);
    tex.needsUpdate=true;
    occTextureCache[key]=tex;
    return tex;
}

function detectCrystalBaseBonds(){
    if(!CRY||!CRY_INV)return;
    var baseAtoms=CRY.baseAtoms;
    var lv=CRY.latticeVectors;
    var inv=CRY_INV;
    var n=baseAtoms.length;
    var CR2={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.85,C:0.76,N:0.71,O:0.66,F:0.57,Ne:0.58,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,Ar:1.06,K:2.03,Ca:1.76,Sc:1.70,Ti:1.60,V:1.53,Cr:1.39,Mn:1.39,Fe:1.32,Co:1.26,Ni:1.24,Cu:1.32,Zn:1.22,Ga:1.22,Ge:1.20,As:1.19,Se:1.20,Br:1.20,Kr:1.16,Rb:2.20,Sr:1.95,Y:1.90,Zr:1.75,Nb:1.64,Mo:1.54,Tc:1.47,Ru:1.46,Rh:1.42,Pd:1.39,Ag:1.45,Cd:1.44,In:1.42,Sn:1.39,Sb:1.39,Te:1.38,Xe:1.40,Cs:2.44,Ba:2.15,La:2.07,Ce:2.04,Pr:2.03,Nd:2.01,Pm:1.99,Sm:1.98,Eu:1.98,Gd:1.96,Tb:1.94,Dy:1.92,Ho:1.92,Er:1.89,Tm:1.90,Yb:1.87,Lu:1.87,Hf:1.75,Ta:1.70,W:1.62,Re:1.51,Os:1.44,Ir:1.41,Pt:1.36,Au:1.36,Hg:1.32,Tl:1.45,Pb:1.46,Bi:1.48,Po:1.40,At:1.50,Rn:1.50,I:1.39};
    var BS={'C+C':[{o:3,l:1.20,t:0.05},{o:1.5,l:1.39,t:0.05},{o:2,l:1.38,t:0.05},{o:1,l:1.51,t:0.10}],'C+N':[{o:2,l:1.26,t:0.05},{o:1.5,l:1.36,t:0.05},{o:1,l:1.43,t:0.10},{o:3,l:1.16,t:0.06}],'C+O':[{o:2,l:1.24,t:0.05},{o:1,l:1.39,t:0.05}],'N+N':[{o:1,l:1.41,t:0.10},{o:2,l:1.25,t:0.06},{o:3,l:1.10,t:0.06}],'N+O':[{o:2,l:1.20,t:0.06},{o:1.5,l:1.30,t:0.06},{o:1,l:1.40,t:0.15}],'O+O':[{o:2,l:1.21,t:0.06},{o:1,l:1.48,t:0.15}],'C+S':[{o:1.5,l:1.73,t:0.06},{o:2,l:1.60,t:0.10},{o:1,l:1.82,t:0.15}],'C+F':[{o:1,l:1.33,t:0.10}],'C+H':[{o:1,l:0.97,t:0.15}],'N+H':[{o:1,l:0.88,t:0.15}],'O+H':[{o:1,l:0.85,t:0.15}]};
    var BC={HH:0,CH:1.3,HO:1.2,HN:1.3,CC:1.9,CO:1.7,CN:1.7,NN:1.7,NO:1.8,CF:1.6,CS:2.0};
    function pk(e1,e2){e1=e1.toUpperCase();e2=e2.toUpperCase();return e1<e2?e1+e2:e2+e1}
    function sk(e1,e2){e1=e1.charAt(0).toUpperCase()+e1.slice(1).toLowerCase();e2=e2.charAt(0).toUpperCase()+e2.slice(1).toLowerCase();return e1<e2?e1+'+'+e2:e2+'+'+e1}
    function gbo(el1,el2,d){
        var p=pk(el1,el2);var co=BC[p];
        if(co!==undefined){if(d>co)return 0}else{var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;if(d>(r1+r2)+0.5)return 0}
        var s=sk(el1,el2);var sp=BS[s];
        if(sp){for(var k=0;k<sp.length;k++){if(Math.abs(d-sp[k].l)<=sp[k].t)return sp[k].o}var bo=1,md=Infinity;for(var k=0;k<sp.length;k++){var df=Math.abs(d-sp[k].l);if(df<md){md=df;bo=sp[k].o}}return bo}
        var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;var rs=r1+r2;var ratio=rs?d/rs:1;
        if(ratio<0.85)return 3;if(ratio<0.90)return 2;return 1;
    }
    function cartToFrac(a){
        return[inv[0][0]*a.x+inv[0][1]*a.y+inv[0][2]*a.z,
               inv[1][0]*a.x+inv[1][1]*a.y+inv[1][2]*a.z,
               inv[2][0]*a.x+inv[2][1]*a.y+inv[2][2]*a.z];
    }
    var bonds=[];
    for(var i=0;i<n;i++){
        var f1=cartToFrac(baseAtoms[i]);
        f1[0]-=Math.floor(f1[0]);f1[1]-=Math.floor(f1[1]);f1[2]-=Math.floor(f1[2]);
        for(var j=i+1;j<n;j++){
            var f2=cartToFrac(baseAtoms[j]);
            f2[0]-=Math.floor(f2[0]);f2[1]-=Math.floor(f2[1]);f2[2]-=Math.floor(f2[2]);
            var fdx=f2[0]-f1[0],fdy=f2[1]-f1[1],fdz=f2[2]-f1[2];
            var sx=Math.round(fdx),sy=Math.round(fdy),sz=Math.round(fdz);
            var wfx=fdx-sx,wfy=fdy-sy,wfz=fdz-sz;
            var dx=wfx*lv[0][0]+wfy*lv[1][0]+wfz*lv[2][0];
            var dy=wfx*lv[0][1]+wfy*lv[1][1]+wfz*lv[2][1];
            var dz=wfx*lv[0][2]+wfy*lv[1][2]+wfz*lv[2][2];
            var dist=Math.sqrt(dx*dx+dy*dy+dz*dz);
            var bo=gbo(baseAtoms[i].element,baseAtoms[j].element,dist);
            if(bo>0){
                bonds.push({atom1:i,atom2:j,order:bo,shift:[-sx,-sy,-sz]});
            }
        }
    }
    CRY.baseBonds=bonds;
}

function rebuildCrystal(){
    if(!CRY)return;
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    var baseAtoms=CRY.baseAtoms;
    var baseBonds=CRY.baseBonds;

    if(!CRY_INV){MD.atoms=[];MD.bonds=[];updateMolInfo();return}
    var inv=CRY_INV;

    if(baseBonds.length===0||!baseBonds[0].shift){
        detectCrystalBaseBonds();
        baseBonds=CRY.baseBonds;
    }

    var adj={};var groupOf={};var groups=[];
    baseBonds.forEach(function(bb){
        if(!adj[bb.atom1])adj[bb.atom1]=[];
        if(!adj[bb.atom2])adj[bb.atom2]=[];
        adj[bb.atom1].push({nb:bb.atom2,shift:bb.shift});
        adj[bb.atom2].push({nb:bb.atom1,shift:bb.shift});
    });
    baseAtoms.forEach(function(ba,bi){
        if(groupOf[bi]!==undefined)return;
        if(!adj[bi]){groupOf[bi]=-1;return}
        var gid=groups.length;var queue=[bi];groupOf[bi]=gid;
        var hasShift=false;
        while(queue.length){
            var cur=queue.shift();
            (adj[cur]||[]).forEach(function(e){
                if(e.shift&&(e.shift[0]||e.shift[1]||e.shift[2]))hasShift=true;
                if(groupOf[e.nb]===undefined){groupOf[e.nb]=gid;queue.push(e.nb)}
            });
        }
        groups.push({hasShift:hasShift});
    });

    var iMin=Math.floor(aMin),iMax=Math.ceil(aMax);
    var jMin=Math.floor(bMin),jMax=Math.ceil(bMax);
    var kMin=Math.floor(cMin),kMax=Math.ceil(cMax);
    var EPS=1e-9;

    var atoms=[];
    var bonds=[];
    var idx=0;
    var cellMap={};

    baseAtoms.forEach(function(ba,bi){
        var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
        var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
        var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
        fx-=Math.floor(fx);fy-=Math.floor(fy);fz-=Math.floor(fz);
        for(var ci=iMin;ci<=iMax;ci++){
            for(var cj=jMin;cj<=jMax;cj++){
                for(var ck=kMin;ck<=kMax;ck++){
                    var ax=fx+ci,ay=fy+cj,az=fz+ck;
                    if(ax<aMin-EPS||ax>aMax+EPS)continue;
                    if(ay<bMin-EPS||ay>bMax+EPS)continue;
                    if(az<cMin-EPS||az>cMax+EPS)continue;
                    var ox=ax*lv[0][0]+ay*lv[1][0]+az*lv[2][0];
                    var oy=ax*lv[0][1]+ay*lv[1][1]+az*lv[2][1];
                    var oz=ax*lv[0][2]+ay*lv[1][2]+az*lv[2][2];
                    var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
                    atoms.push({
                        element:ba.element,
                        x:ox,y:oy,z:oz,
                        index:idx,
                        occupancy:ba.occupancy,
                        color:color,
                        baseIdx:ba.baseIdx!=null?ba.baseIdx:bi,
                        cellI:ci,cellJ:cj,cellK:ck
                    });
                    cellMap[bi+'_'+ci+'_'+cj+'_'+ck]=idx;
                    idx++;
                }
            }
        }
    });

    var activeGroups={};
    atoms.forEach(function(a){
        var g=groupOf[a.baseIdx];
        if(g!==undefined&&g>=0&&!groups[g].hasShift){
            activeGroups[g+'_'+a.cellI+'_'+a.cellJ+'_'+a.cellK]=true;
        }
    });
    baseAtoms.forEach(function(ba,bi){
        var g=groupOf[bi];
        if(g===undefined||g<0||groups[g].hasShift)return;
        var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
        var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
        var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
        fx-=Math.floor(fx);fy-=Math.floor(fy);fz-=Math.floor(fz);
        for(var ci=iMin;ci<=iMax;ci++){
            for(var cj=jMin;cj<=jMax;cj++){
                for(var ck=kMin;ck<=kMax;ck++){
                    if(!activeGroups[g+'_'+ci+'_'+cj+'_'+ck])continue;
                    if(cellMap[bi+'_'+ci+'_'+cj+'_'+ck]!==undefined)continue;
                    var ax=fx+ci,ay=fy+cj,az=fz+ck;
                    var ox=ax*lv[0][0]+ay*lv[1][0]+az*lv[2][0];
                    var oy=ax*lv[0][1]+ay*lv[1][1]+az*lv[2][1];
                    var oz=ax*lv[0][2]+ay*lv[1][2]+az*lv[2][2];
                    var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
                    atoms.push({
                        element:ba.element,
                        x:ox,y:oy,z:oz,
                        index:idx,
                        occupancy:ba.occupancy,
                        color:color,
                        baseIdx:ba.baseIdx!=null?ba.baseIdx:bi,
                        cellI:ci,cellJ:cj,cellK:ck
                    });
                    cellMap[bi+'_'+ci+'_'+cj+'_'+ck]=idx;
                    idx++;
                }
            }
        }
    });

    var seenPairs={};
    var seenSplits={};
    baseBonds.forEach(function(bb){
        var b1=bb.atom1,b2=bb.atom2;
        var s=bb.shift||[0,0,0];
        var hasShift=s[0]||s[1]||s[2];
        for(var ci1=iMin;ci1<=iMax;ci1++){
            for(var cj1=jMin;cj1<=jMax;cj1++){
                for(var ck1=kMin;ck1<=kMax;ck1++){
                    var idx1=cellMap[b1+'_'+ci1+'_'+cj1+'_'+ck1];
                    if(idx1==null)continue;
                    var ci2=ci1+s[0],cj2=cj1+s[1],ck2=ck1+s[2];
                    var idx2=cellMap[b2+'_'+ci2+'_'+cj2+'_'+ck2];
                    if(idx2!=null){
                        var key=idx1<idx2?idx1+'_'+idx2:idx2+'_'+idx1;
                        if(!seenPairs[key]){
                            seenPairs[key]=true;
                            bonds.push({atom1:idx1,atom2:idx2,order:bb.order});
                        }
                    }else{
                        var idx2s=cellMap[b2+'_'+ci1+'_'+cj1+'_'+ck1];
                        if(idx2s!=null){
                            var skey='F_'+idx1+'_'+idx2s+'_'+s[0]+'_'+s[1]+'_'+s[2];
                            if(!seenSplits[skey]){
                                seenSplits[skey]=true;
                                bonds.push({atom1:idx1,atom2:idx2s,order:bb.order,shift:s});
                            }
                        }
                    }
                }
            }
        }
        if(hasShift){
            var rs=[-s[0],-s[1],-s[2]];
            for(var ci2=iMin;ci2<=iMax;ci2++){
                for(var cj2=jMin;cj2<=jMax;cj2++){
                    for(var ck2=kMin;ck2<=kMax;ck2++){
                        var idx2r=cellMap[b2+'_'+ci2+'_'+cj2+'_'+ck2];
                        if(idx2r==null)continue;
                        var ri1=ci2+rs[0],rj1=cj2+rs[1],rk1=ck2+rs[2];
                        var ridx1=cellMap[b1+'_'+ri1+'_'+rj1+'_'+rk1];
                        if(ridx1==null){
                            var idx1r=cellMap[b1+'_'+ci2+'_'+cj2+'_'+ck2];
                            if(idx1r!=null){
                                var skey='R_'+idx2r+'_'+idx1r+'_'+rs[0]+'_'+rs[1]+'_'+rs[2];
                                if(!seenSplits[skey]){
                                    seenSplits[skey]=true;
                                    bonds.push({atom1:idx2r,atom2:idx1r,order:bb.order,shift:rs});
                                }
                            }
                        }
                    }
                }
            }
        }
    });

    MD.atoms=atoms;
    MD.bonds=bonds;
    updateMolInfo();
    needsRender=true;
}

function removeDisorder(){
    if(!CRY)return;
    pushUndo();
    var baseCountBefore=CRY.baseAtoms.length;
    var newBaseAtoms=[];
    var oldToNew={};
    CRY.baseAtoms.forEach(function(ba,i){
        var occ=ba.occupancy!=null?ba.occupancy:1;
        if(occ<0.5)return;
        var newIdx=newBaseAtoms.length;
        oldToNew[i]=newIdx;
        var newOcc=occ>0.5?1:occ;
        var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
        newBaseAtoms.push({
            element:ba.element,x:ba.x,y:ba.y,z:ba.z,
            index:newIdx,occupancy:newOcc,baseIdx:newIdx,color:color
        });
    });
    CRY.baseBonds=CRY.baseBonds.filter(function(b){
        return oldToNew[b.atom1]!=null&&oldToNew[b.atom2]!=null;
    }).map(function(b){
        return{atom1:oldToNew[b.atom1],atom2:oldToNew[b.atom2],order:b.order,shift:b.shift};
    });
    CRY.baseAtoms=newBaseAtoms;
    if(CRY.baseAtoms.length<baseCountBefore)provMark('disorder',1);
    rebuildCrystal();
    rebuildScene();
}

var cellWireframe=null;
function buildCellWireframe(){
    if(cellWireframe){moleculeGroup.remove(cellWireframe);cellWireframe=null}
    if(!CRY)return;
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    var corners=[];
    for(var i=0;i<2;i++)for(var j=0;j<2;j++)for(var k=0;k<2;k++){
        var ai=i?aMax:aMin,aj=j?bMax:bMin,ak=k?cMax:cMin;
        var x=ai*lv[0][0]+aj*lv[1][0]+ak*lv[2][0];
        var y=ai*lv[0][1]+aj*lv[1][1]+ak*lv[2][1];
        var z=ai*lv[0][2]+aj*lv[1][2]+ak*lv[2][2];
        corners.push(new THREE.Vector3(x-CX,y-CY,z-CZ));
    }
    var edges=[[0,4],[1,5],[2,6],[3,7],[0,2],[1,3],[4,6],[5,7],[0,1],[2,3],[4,5],[6,7]];
    var positions=[];
    edges.forEach(function(e){
        positions.push(corners[e[0]].x,corners[e[0]].y,corners[e[0]].z);
        positions.push(corners[e[1]].x,corners[e[1]].y,corners[e[1]].z);
    });
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
    var mat=new THREE.LineBasicMaterial({color:0x88aaff});
    cellWireframe=new THREE.LineSegments(geo,mat);
    moleculeGroup.add(cellWireframe);
}

function updateAxesIndicator(){
    var el=document.getElementById('axes-indicator');
    if(!el)return;
    var vectors;
    if(CRY){
        var lv=CRY.latticeVectors;
        vectors=[
            {dir:new THREE.Vector3(lv[0][0],lv[0][1],lv[0][2]).normalize(),color:'#FF4444',label:'a'},
            {dir:new THREE.Vector3(lv[1][0],lv[1][1],lv[1][2]).normalize(),color:'#44FF44',label:'b'},
            {dir:new THREE.Vector3(lv[2][0],lv[2][1],lv[2][2]).normalize(),color:'#4488FF',label:'c'}
        ];
    }else{
        // Molecular structure: show the Cartesian X/Y/Z axes of the molecule's
        // initial coordinate frame. rotQuat starts as identity, so at load time
        // the indicator matches the world axes; applying rotQuat afterwards
        // shows the current orientation no matter how the molecule is rotated
        // or edited (orientation is purely rotQuat-driven).
        vectors=[
            {dir:new THREE.Vector3(1,0,0),color:'#FF4444',label:'X'},
            {dir:new THREE.Vector3(0,1,0),color:'#44FF44',label:'Y'},
            {dir:new THREE.Vector3(0,0,1),color:'#4488FF',label:'Z'}
        ];
    }
    vectors.forEach(function(v){v.dir.applyQuaternion(rotQuat)});
    vectors.sort(function(a,b){return a.dir.z-b.dir.z});
    var cx=45,cy=45,scale=30;
    var svg='<svg viewBox="0 0 90 90">';
    vectors.forEach(function(v){
        var x2=cx+v.dir.x*scale;
        var y2=cy-v.dir.y*scale;
        var depth=(v.dir.z+1)/2;
        var opacity=0.4+depth*0.6;
        var width=1.5+depth*1.5;
        svg+='<line x1="'+cx+'" y1="'+cy+'" x2="'+x2.toFixed(1)+'" y2="'+y2.toFixed(1)+'" stroke="'+v.color+'" stroke-width="'+width.toFixed(1)+'" stroke-linecap="round" opacity="'+opacity.toFixed(2)+'"/>';
        svg+='<circle cx="'+x2.toFixed(1)+'" cy="'+y2.toFixed(1)+'" r="3" fill="'+v.color+'" opacity="'+opacity.toFixed(2)+'"/>';
        svg+='<text x="'+(x2+4).toFixed(1)+'" y="'+(y2+4).toFixed(1)+'" fill="'+v.color+'" font-size="10" font-weight="bold" opacity="'+opacity.toFixed(2)+'">'+v.label+'</text>';
    });
    svg+='<circle cx="'+cx+'" cy="'+cy+'" r="2" fill="#888"/>';
    svg+='</svg>';
    el.innerHTML=svg;
}

var atomMeshes=[];
var currentBondIdx=-1;
var needsRender=true;
var contextLost=false;      // the WebGL context is gone until it is restored

// ===== Shared geometries and materials =====
// Every drawn thing used to own a geometry and a material. A 2175-atom file
// built 2175 sphere geometries plus 13906 bond cylinders, each of the 16081
// objects with its own material, and the first frame had to create roughly
// 32000 GL buffers. A per-object material is the expensive half: three.js
// re-uploads a material's whole uniform set before every draw call whose
// material id changed, so ~16k distinct materials meant roughly a million
// uniform uploads per frame — enough driver work in one frame to reset the GPU
// process on Windows (driver timeout), and once that shared GPU process is
// gone *no* webview in the window can create a WebGL context again, so every
// later file fails with three.js's "Error creating WebGL context.".
//
// Two pools remove the resource half of that: a sphere geometry is shared by
// every atom with the same radius, a material by every atom with the same
// colour, occupancy and highlight state, and every bond cylinder is drawn by
// one instanced mesh. A pooled resource outlives the meshes that draw with it,
// so it is marked, and disposeMesh leaves it alone.
var geoPool={};
var matPool={};
function pooledGeo(key,make){
    var g=geoPool[key];
    if(!g){g=make();g.userData={shared:true};geoPool[key]=g}
    return g;
}
function pooledMat(key,make){
    var m=matPool[key];
    if(!m){m=make();m.userData={shared:true};matPool[key]=m}
    return m;
}
function sphereGeoFor(r,segs){
    return pooledGeo('s|'+r.toFixed(4)+'|'+segs[0]+'|'+segs[1],function(){return new THREE.SphereGeometry(r,segs[0],segs[1])});
}
// An atom's colour follows the scalar-colour toggle (a scene rebuild), and its
// highlight is a *state*: plain, selected, or diff-highlighted. Both used to be
// written into the atom's own material, which is exactly why every atom needed
// one. Keying the pool by (colour, occupancy, state) keeps the same three
// appearances while a whole structure shares a handful of materials.
function atomMaterial(col,occ,state){
    var o=(occ!=null&&occ<0.999)?occ.toFixed(3):'1';
    return pooledMat('a|'+col+'|'+o+'|'+state,function(){
        var m;
        if(o!=='1')m=new THREE.MeshPhongMaterial({color:0xffffff,map:getOccupancyTexture(col,occ),shininess:80,specular:0x444444});
        else m=new THREE.MeshPhongMaterial({color:new THREE.Color(col),shininess:80,specular:0x444444});
        if(state===1){m.emissive=new THREE.Color(0xffff00);m.emissiveIntensity=0.6}
        else if(state===2){m.emissive=new THREE.Color(0xff6600);m.emissiveIntensity=0.4}
        return m;
    });
}
// Ask for one of those three states. Every writer of the highlight goes
// through here, so no call site can leave a shared material in a changed state
// for the atoms it is shared with.
function setAtomHighlight(mesh,state){
    if(!mesh||!mesh.userData)return;
    mesh.material=atomMaterial(mesh.userData.baseColor!=null?mesh.userData.baseColor:'#FF1493',mesh.userData.occ,state);
}
// The comparison's atoms follow the same idea, with the diff highlight being a
// colour of its own (0xff6600, the way the comparison has always marked a
// difference) rather than a glow over the structure's colour.
function diffAtomMaterial(col,state){
    return pooledMat('d|'+col+'|'+state,function(){
        var m=new THREE.MeshPhongMaterial({color:new THREE.Color(state===2?0xff6600:col),shininess:80,specular:0x444444});
        if(state===1){m.emissive=new THREE.Color(0xffff00);m.emissiveIntensity=0.6}
        else if(state===2){m.emissive=new THREE.Color(0xff6600);m.emissiveIntensity=0.4}
        return m;
    });
}
function setDiffAtomHighlight(mesh,state){
    if(!mesh||!mesh.userData)return;
    mesh.material=diffAtomMaterial(mesh.userData.baseColor!=null?mesh.userData.baseColor:'#FF1493',state);
}

function disposeMesh(ch){
    // A pooled geometry/material is still drawing every other atom or bond
    // with the same shape or colour, so it must not be released with the one
    // mesh that happened to be disposed first.
    if(ch.geometry&&!(ch.geometry.userData&&ch.geometry.userData.shared))ch.geometry.dispose();
    if(ch.material&&!(ch.material.userData&&ch.material.userData.shared)){
        if(ch.material.map)ch.material.map.dispose();
        ch.material.dispose();
    }
}
function sphereSegs(n){
    if(n==null)n=MD.atoms.length;
    if(n<=200)return [32,24];
    if(n<=1000)return [20,16];
    return [12,10];
}

// ===== Bonds: one instanced mesh for every cylinder =====
// A bond used to be two to eight separate meshes — a cylinder per half, six
// more for an aromatic bond's dashed line — so 3309 bonds became 13906 meshes,
// which is where most of the object count came from. A bond piece is now a
// *record*: one shared unit cylinder draws all of them, with the piece's own
// radius and length in the instance matrix and its colour in the instance
// colour buffer. Placement lives in one routine that reads a piece's transform
// off the two atoms it spans — the very arithmetic normal-mode playback
// already used — so building a piece and animating it can no longer drift.
var bondPieceIdx=[];        // bond index -> [piece slot, ...]
var _bpDir=new THREE.Vector3(),_bpQuat=new THREE.Quaternion();
var _bpM=new THREE.Matrix4(),_bpPos=new THREE.Vector3(),_bpScale=new THREE.Vector3();
var _bpCol=new THREE.Color(),_bpY=new THREE.Vector3(0,1,0);

// One instanced cylinder pool: a unit cylinder (radius 1, height 1, centred),
// its own material, and the piece list that fills it. The main structure draws
// every bond from the single pool below; the comparison draws its bonds from
// two, because a highlighted bond carries an orange glow in its material.
function newBondPool(shininess,emiCol,emiInt,parent){
    return {
        pieces:[],
        inst:null,
        cap:0,
        parent:parent,
        mat:pooledMat('bond|'+shininess+'|'+emiCol+'|'+emiInt,function(){
            var m=new THREE.MeshPhongMaterial({color:0xffffff,shininess:shininess,specular:0x222222});
            if(emiCol){m.emissive=new THREE.Color(emiCol);m.emissiveIntensity=emiInt}
            return m;
        })
    };
}
function unitCylinderGeo(){
    return pooledGeo('cyl8',function(){return new THREE.CylinderGeometry(1,1,1,8,1)});
}
// Grow (never shrink) the pool's instanced mesh. The mesh is only recreated
// when the piece count outgrows the capacity, so an edit that only moves atoms
// reuses the buffers it already uploaded.
function poolEnsure(pool,n){
    if(pool.inst&&pool.cap>=n){
        // A scene rebuild clears the group the pool draws into, so the mesh
        // has to be put back before it is drawn again — its buffers are still
        // uploaded and its size is unchanged, which is the whole point of
        // keeping the capacity.
        if(pool.parent&&pool.inst.parent!==pool.parent)pool.parent.add(pool.inst);
        return;
    }
    if(pool.inst){
        if(pool.inst.parent)pool.inst.parent.remove(pool.inst);
        pool.inst.dispose();
        pool.inst=null;
    }
    var cap=1;
    while(cap<n)cap*=2;
    var inst=new THREE.InstancedMesh(unitCylinderGeo(),pool.mat,cap);
    // The instances span the whole structure while the shared unit geometry is
    // one unit wide, so the mesh's own bounding sphere says nothing about where
    // its bonds are and three.js would cull the lot as soon as the origin left
    // the view. One draw call is cheap enough to always draw.
    inst.frustumCulled=false;
    _bpCol.setRGB(1,1,1);
    for(var i=0;i<cap;i++)inst.setColorAt(i,_bpCol);
    inst.count=0;
    pool.inst=inst;pool.cap=cap;
    if(pool.parent)pool.parent.add(inst);
}
// Write one piece's transform and colour. The cylinder's midpoint is the
// instance position, its length the y scale and its radius the x/z scale, so a
// single unit cylinder covers every radius and length in the structure. The
// colour goes straight into the instance colour buffer — no temporary colour
// object per piece per frame.
function poolPut(pool,slot,px,py,pz,dx,dy,dz,l,r,col){
    _bpDir.set(dx,dy,dz).normalize();
    _bpQuat.setFromUnitVectors(_bpY,_bpDir);
    _bpM.compose(_bpPos.set(px,py,pz),_bpQuat,_bpScale.set(r,l,r));
    pool.inst.setMatrixAt(slot,_bpM);
    var arr=pool.inst.instanceColor.array;
    arr[slot*3]=col.r;arr[slot*3+1]=col.g;arr[slot*3+2]=col.b;
}
// A piece that cannot be placed (its bond or an atom is gone, or the two atoms
// sit on top of each other) is collapsed to zero size rather than left wearing
// the transform it had a frame ago.
function poolHide(pool,slot){
    _bpM.makeScale(0,0,0);
    pool.inst.setMatrixAt(slot,_bpM);
}
function poolUpload(pool){
    if(!pool.inst)return;
    // The count field is how many instances the mesh *draws*, and it is the one
    // field three.js does not take from the buffers: a new InstancedMesh starts
    // at the capacity, but a pool that is refilled after an edit or a rebuild
    // has to state its piece count here or it draws whatever it was left at (a
    // fresh mesh is created with count 0, so forgetting this draws no bonds at
    // all).
    pool.inst.count=pool.pieces.length;
    pool.inst.instanceMatrix.needsUpdate=true;
    if(pool.inst.instanceColor)pool.inst.instanceColor.needsUpdate=true;
}
function poolClear(pool){
    if(pool.inst){
        if(pool.inst.parent)pool.inst.parent.remove(pool.inst);
        pool.inst.dispose();
        pool.inst=null;pool.cap=0;
    }
    pool.pieces.length=0;
}

// The main structure's bonds: one pool, one draw call, whatever the size.
var bondPool=newBondPool(40,0,0,moleculeGroup);

// ===== VMD Surf patches =====
// A "mol representation Surf <probe> <isovalue>" patch is the isosurface of
// f(p) = min_i(|p - c_i| - (R_vdw,i + probe)) over the atoms that
// representation selects, i.e. the solvent-accessible surface of that group of
// atoms. It is extracted here as a real mesh by marching tetrahedra (each grid
// cube is split into six tetrahedra, which needs no ambiguous-case table), and
// drawn translucent: a script stacks many patches, and opaque ones would hide
// everything but the outermost.
var VMD_SURF_CELL=0.45;
var VMD_SURF_MAX_DIV=64;
var VMD_SURF_OPACITY=0.38;
// The six tetrahedra of a cube, all sharing the 0-7 diagonal.
var VMD_TETRA=[[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]];
var VMD_CORNER_X=[0,1,0,1,0,1,0,1];
var VMD_CORNER_Y=[0,0,1,1,0,0,1,1];
var VMD_CORNER_Z=[0,0,0,0,1,1,1,1];
var vmdSurfaceSignature='';

var VMD_TET_IDX=new Int32Array(4);
var VMD_TET_IN=new Int32Array(4);
var VMD_TET_OUT=new Int32Array(4);
// Scratch state, so that extracting a patch allocates nothing per tetrahedron.
function vmdEmitTet(pos,px,py,pz,v,a,b,c,d){
    var t,i,ins=0,out=0;
    VMD_TET_IDX[0]=a;VMD_TET_IDX[1]=b;VMD_TET_IDX[2]=c;VMD_TET_IDX[3]=d;
    for(i=0;i<4;i++){t=VMD_TET_IDX[i];if(v[t]<0)VMD_TET_IN[ins++]=t;else VMD_TET_OUT[out++]=t}
    if(ins===0||ins===4)return;
    var p,q,r;
    if(ins===1){
        var i0=VMD_TET_IN[0];
        p=VMD_TET_OUT[0];q=VMD_TET_OUT[1];r=VMD_TET_OUT[2];
        var t0=v[i0]/(v[i0]-v[p]),t1=v[i0]/(v[i0]-v[q]),t2=v[i0]/(v[i0]-v[r]);
        pos.push(px[i0]+(px[p]-px[i0])*t0,py[i0]+(py[p]-py[i0])*t0,pz[i0]+(pz[p]-pz[i0])*t0,
                 px[i0]+(px[q]-px[i0])*t1,py[i0]+(py[q]-py[i0])*t1,pz[i0]+(pz[q]-pz[i0])*t1,
                 px[i0]+(px[r]-px[i0])*t2,py[i0]+(py[r]-py[i0])*t2,pz[i0]+(pz[r]-pz[i0])*t2);
    }else if(ins===3){
        p=VMD_TET_IN[0];q=VMD_TET_IN[1];r=VMD_TET_IN[2];
        var o0=VMD_TET_OUT[0];
        var u0=v[p]/(v[p]-v[o0]),u1=v[q]/(v[q]-v[o0]),u2=v[r]/(v[r]-v[o0]);
        pos.push(px[p]+(px[o0]-px[p])*u0,py[p]+(py[o0]-py[p])*u0,pz[p]+(pz[o0]-pz[p])*u0,
                 px[q]+(px[o0]-px[q])*u1,py[q]+(py[o0]-py[q])*u1,pz[q]+(pz[o0]-pz[q])*u1,
                 px[r]+(px[o0]-px[r])*u2,py[r]+(py[o0]-py[r])*u2,pz[r]+(pz[o0]-pz[r])*u2);
    }else{
        var i1=VMD_TET_IN[0],i2=VMD_TET_IN[1],o1=VMD_TET_OUT[0],o2=VMD_TET_OUT[1];
        var c0=v[i1]/(v[i1]-v[o1]),c1=v[i1]/(v[i1]-v[o2]);
        var c2=v[i2]/(v[i2]-v[o2]),c3=v[i2]/(v[i2]-v[o1]);
        var ax=px[i1]+(px[o1]-px[i1])*c0,ay=py[i1]+(py[o1]-py[i1])*c0,az=pz[i1]+(pz[o1]-pz[i1])*c0;
        var bx=px[i1]+(px[o2]-px[i1])*c1,by=py[i1]+(py[o2]-py[i1])*c1,bz=pz[i1]+(pz[o2]-pz[i1])*c1;
        var gx=px[i2]+(px[o2]-px[i2])*c2,gy=py[i2]+(py[o2]-py[i2])*c2,gz=pz[i2]+(pz[o2]-pz[i2])*c2;
        var hx=px[i2]+(px[o1]-px[i2])*c3,hy=py[i2]+(py[o1]-py[i2])*c3,hz=pz[i2]+(pz[o1]-pz[i2])*c3;
        pos.push(ax,ay,az,bx,by,bz,gx,gy,gz);
        pos.push(ax,ay,az,gx,gy,gz,hx,hy,hz);
    }
}

function buildVmdSurfaceGeometry(indices,probe){
    var i,k,a;
    var cx=[],cy=[],cz=[],cr=[];
    var minx=Infinity,miny=Infinity,minz=Infinity,maxx=-Infinity,maxy=-Infinity,maxz=-Infinity;
    for(i=0;i<indices.length;i++){
        a=MD.atoms[indices[i]];
        if(!a)continue;
        var r=(typeof a.vdw==='number'?a.vdw:1.5)+probe;
        if(!(r>0))continue;
        cx.push(a.x);cy.push(a.y);cz.push(a.z);cr.push(r);
        if(a.x-r<minx)minx=a.x-r;
        if(a.x+r>maxx)maxx=a.x+r;
        if(a.y-r<miny)miny=a.y-r;
        if(a.y+r>maxy)maxy=a.y+r;
        if(a.z-r<minz)minz=a.z-r;
        if(a.z+r>maxz)maxz=a.z+r;
    }
    var count=cx.length;
    if(count===0)return null;
    var spanX=maxx-minx,spanY=maxy-miny,spanZ=maxz-minz;
    var span=Math.max(spanX,spanY,spanZ);
    if(!(span>0))return null;
    var cell=VMD_SURF_CELL;
    var div=Math.ceil(span/cell);
    if(div>VMD_SURF_MAX_DIV){div=VMD_SURF_MAX_DIV;cell=span/div}
    var nx=Math.max(2,Math.ceil(spanX/cell)+1);
    var ny=Math.max(2,Math.ceil(spanY/cell)+1);
    var nz=Math.max(2,Math.ceil(spanZ/cell)+1);
    var field=new Float32Array(nx*ny*nz);
    for(var kz=0;kz<nz;kz++){
        var z=minz+kz*cell;
        for(var ky=0;ky<ny;ky++){
            var y=miny+ky*cell;
            var row=(kz*ny+ky)*nx;
            for(var kx=0;kx<nx;kx++){
                var x=minx+kx*cell;
                var best=Infinity;
                for(k=0;k<count;k++){
                    var dx=x-cx[k],dy=y-cy[k],dz=z-cz[k];
                    var d=Math.sqrt(dx*dx+dy*dy+dz*dz)-cr[k];
                    if(d<best)best=d;
                }
                // A sample of exactly zero would be classified as "outside", and
                // a crossing towards it would land exactly on the grid corner:
                // every crossing of that tetrahedron would collapse onto the
                // same point and produce zero-area triangles. Nudging it keeps
                // every crossing strictly inside its edge — the surface moves by
                // far less than the grid spacing, and the mesh stays closed.
                if(best===0)best=1e-9;
                field[row+kx]=best;
            }
        }
    }
    var pos=[];
    var v=new Array(8),px=new Array(8),py=new Array(8),pz=new Array(8);
    // The field changes by at most the cell diagonal across a cell, so a cell
    // whose smallest corner value exceeds that cannot hold any part of the
    // surface and is skipped.
    var cellDiag=cell*1.7320508;
    for(kz=0;kz<nz-1;kz++){
        for(ky=0;ky<ny-1;ky++){
            for(kx=0;kx<nx-1;kx++){
                var base=(kz*ny+ky)*nx+kx;
                var vx0=minx+kx*cell,vy0=miny+ky*cell,vz0=minz+kz*cell;
                var nearest=Infinity;
                for(k=0;k<8;k++){
                    px[k]=vx0+VMD_CORNER_X[k]*cell;
                    py[k]=vy0+VMD_CORNER_Y[k]*cell;
                    pz[k]=vz0+VMD_CORNER_Z[k]*cell;
                    var gx=kx+VMD_CORNER_X[k],gy=ky+VMD_CORNER_Y[k],gz=kz+VMD_CORNER_Z[k];
                    v[k]=field[(gz*ny+gy)*nx+gx];
                    if(v[k]<nearest)nearest=v[k];
                }
                if(nearest>cellDiag)continue;
                for(k=0;k<6;k++){
                    var q=VMD_TETRA[k];
                    vmdEmitTet(pos,px,py,pz,v,q[0],q[1],q[2],q[3]);
                }
            }
        }
    }
    if(pos.length===0)return null;
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
    geo.computeVertexNormals();
    return geo;
}

// Rebuild the surface patches so they match the current geometry. Cheap when
// nothing moved (the atom coordinates the patches are built from are compared
// with the ones they were built from).
function clearVmdSurfaceMeshes(){
    for(var i=vmdSurfaceMeshes.length-1;i>=0;i--){
        disposeMesh(vmdSurfaceMeshes[i]);
        if(vmdSurfaceGroup)vmdSurfaceGroup.remove(vmdSurfaceMeshes[i]);
    }
    vmdSurfaceMeshes.length=0;
}
function refreshVmdSurfaces(){
    if(!vmdSurfaceGroup)return;
    var segs=(VMD&&VMD.activeSegments)||[];
    var sig='';
    for(var s=0;s<segs.length;s++){
        sig+=segs[s].probe+';';
        var idx=segs[s].indices;
        for(var k=0;k<idx.length;k++){
            var a=MD.atoms[idx[k]];
            if(!a)continue;
            sig+=a.x.toFixed(3)+','+a.y.toFixed(3)+','+a.z.toFixed(3)+';';
        }
    }
    sig+='|'+CX.toFixed(3)+','+CY.toFixed(3)+','+CZ.toFixed(3);
    if(sig===vmdSurfaceSignature&&vmdSurfaceMeshes.length>0){
        vmdSurfaceGroup.visible=true;
        return;
    }
    vmdSurfaceSignature=sig;
    clearVmdSurfaceMeshes();
    for(s=0;s<segs.length;s++){
        var geo=null;
        try{geo=buildVmdSurfaceGeometry(segs[s].indices,segs[s].probe)}catch(e){geo=null}
        if(!geo)continue;
        geo.translate(-CX,-CY,-CZ);
        var mat=new THREE.MeshPhongMaterial({color:new THREE.Color(segs[s].color),transparent:true,opacity:VMD_SURF_OPACITY,depthWrite:false,side:THREE.DoubleSide,shininess:30,specular:0x222222});
        var mesh=new THREE.Mesh(geo,mat);
        mesh.renderOrder=4;
        mesh.userData.vmdSurface=true;
        mesh.userData.segmentIndex=s;
        vmdSurfaceGroup.add(mesh);
        vmdSurfaceMeshes.push(mesh);
    }
    vmdSurfaceGroup.visible=vmdSurfaceMeshes.length>0;
    // The egg panel mirrors the patches: rebuild its rows for the new segment
    // list, then put the panel's shown/hidden choices on the fresh meshes.
    rebuildEggPanel();
    applyVmdSegmentVisibility();
    needsRender=true;
}

// Patches are hidden while a preview moves the atoms (a drag or a slider),
// because rebuilding them on every preview frame would not stay interactive;
// they are brought back by the next refreshVmdSurfaces() at the commit.
function hideVmdSurfaces(){
    if(vmdSurfaceMeshes.length>0&&vmdSurfaceGroup)vmdSurfaceGroup.visible=false;
}

// ===== Egg list =====
// The toolchain writes an eggFile next to the structure ("U12F1 : 1,10-15" =
// the egg named U12F1 wraps units 1 and 10..15) plus a unitfile with every
// unit's atoms; the extension matches each Surf patch to its egg by atom set
// and sends the names along with the payload. The panel lists the eggs and
// toggles the patches; all of them are shown initially.
function escapeHtml(s){
    return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
// "1, 10-15" style: compress consecutive unit numbers into ranges.
function formatUnitList(units){
    var u=units.slice().sort(function(a,b){return a-b});
    var parts=[],i=0;
    while(i<u.length){
        var j=i;
        while(j+1<u.length&&u[j+1]===u[j]+1)j++;
        if(j-i>=2)parts.push(u[i]+'\u2013'+u[j]);        else{for(var k=i;k<=j;k++)parts.push(String(u[k]))}
        i=j+1;
    }
    return parts.join(', ');
}
function eggSegments(){
    var list=[];
    var segs=(VMD&&VMD.activeSegments)||[];
    for(var i=0;i<segs.length;i++){if(segs[i]&&segs[i].egg)list.push({idx:i,seg:segs[i]})}
    return list;
}
function eggShownCount(list){
    var n=0;
    for(var i=0;i<list.length;i++){if(!vmdHiddenSegments[list[i].idx])n++}
    return n;
}
function eggAllHidden(list){
    return list.length>0&&eggShownCount(list)===0;
}
// Apply the panel's choices to the built patches. Called after a rebuild (the
// meshes are new) and on every toggle.
function applyVmdSegmentVisibility(){
    for(var i=0;i<vmdSurfaceMeshes.length;i++){
        var m=vmdSurfaceMeshes[i];
        if(m&&m.userData&&typeof m.userData.segmentIndex==='number')m.visible=!vmdHiddenSegments[m.userData.segmentIndex];
    }
    needsRender=true;
}
function updateEggButton(){
    var btn=document.getElementById('egg-btn');
    if(!btn)return;
    wireEggPanelButton();
    btn.style.display=eggSegments().length>0?'':'none';
}
function updateEggPanelState(list){
    var panel=document.getElementById('egg-panel');
    if(!panel)return;
    var rows=panel.querySelectorAll('.ep-row');
    for(var i=0;i<rows.length;i++){
        var idx=rows[i].getAttribute('data-seg');
        var hidden=!!vmdHiddenSegments[idx];
        rows[i].classList.toggle('ep-hidden',hidden);
        var cb=rows[i].querySelector('input[type=checkbox]');
        if(cb)cb.checked=!hidden;
    }
    var head=document.getElementById('egg-count');
    if(head)head.textContent=eggShownCount(list)+'/'+list.length+' shown';
    var toggleAll=document.getElementById('egg-toggle-all');
    if(toggleAll)toggleAll.textContent=eggAllHidden(list)?'Show all':'Hide all';
}
function rebuildEggPanel(){
    var panel=document.getElementById('egg-panel');
    var btn=document.getElementById('egg-btn');
    if(!panel||!btn)return;
    var list=eggSegments();
    updateEggButton();
    if(list.length===0){
        // The current molecule draws no eggs: no panel, choices forgotten.
        panel.classList.remove('show');
        vmdHiddenSegments={};
        return;
    }
    var html='<div class="ep-head"><span>Eggs <span id="egg-count">'+eggShownCount(list)+'/'+list.length+' shown</span></span>'+
             '<span class="ep-actions"><button id="egg-toggle-all">'+(eggAllHidden(list)?'Show all':'Hide all')+'</button><span id="egg-close" class="ep-close" title="Hide list">\u00d7</span></span></div>';
    for(var i=0;i<list.length;i++){
        var e=list[i].seg.egg;
        var hidden=!!vmdHiddenSegments[list[i].idx];
        html+='<label class="ep-row'+(hidden?' ep-hidden':'')+'" data-seg="'+list[i].idx+'" title="'+escapeHtml(e.name)+'" style="border-left-color:'+list[i].seg.color+'">'+
              '<input type="checkbox" data-seg="'+list[i].idx+'"'+(hidden?'':' checked')+'>'+
              '<span class="ep-name">'+escapeHtml(e.name)+'</span>'+
              '<span class="ep-units">units '+formatUnitList(e.units)+'</span></label>';
    }
    panel.innerHTML=html;
    panel.classList.add('show');
    var boxes=panel.querySelectorAll('input[type=checkbox]');
    for(var b=0;b<boxes.length;b++){
        boxes[b].addEventListener('change',function(){
            var idx=this.getAttribute('data-seg');
            vmdHiddenSegments[idx]=!this.checked;
            applyVmdSegmentVisibility();
            updateEggPanelState(list);
        });
    }
    var toggleAll=document.getElementById('egg-toggle-all');
    if(toggleAll){
        // 标签在点击处理器里也直接刷新一次, 不只依赖后面的对账函数。
        toggleAll.addEventListener('click',function(){
            // 一键全部显示/全部隐藏: 有 egg 显示着就全部隐藏, 全隐藏时全部显示。
            var hide=!eggAllHidden(list);
            for(var i=0;i<list.length;i++)vmdHiddenSegments[list[i].idx]=hide;
            applyVmdSegmentVisibility();
            this.textContent=hide?'Show all':'Hide all';
            updateEggPanelState(list);
        });
    }
    var close=document.getElementById('egg-close');
    // 关闭只隐藏面板; 工具栏的 Eggs 按钮随时可以再打开。
    if(close)close.addEventListener('click',function(){panel.classList.remove('show')});
}
var eggBtnWired=false;
function wireEggPanelButton(){
    if(eggBtnWired)return;
    var btn=document.getElementById('egg-btn');
    if(!btn)return;
    eggBtnWired=true;
    // 工具栏按钮永远"打开"面板 (幂等): 关闭只由面板内的 × 负责, 这样任何
    // 时候都有地方可以再打开。重复点击只是刷新面板内容。
    btn.addEventListener('click',function(){rebuildEggPanel()});
}

function rebuildScene(){
    // The surface group rides inside the molecule group (so that it keeps
    // following the view transform); lift it out before the clear, drop its
    // meshes here, and put it back — refreshVmdSurfaces rebuilds them below.
    if(vmdSurfaceGroup&&vmdSurfaceGroup.parent===moleculeGroup)moleculeGroup.remove(vmdSurfaceGroup);
    clearVmdSurfaceMeshes();
    var oldChildren=moleculeGroup.children.slice();
    oldChildren.forEach(function(ch){disposeMesh(ch);moleculeGroup.remove(ch)});
    if(vmdSurfaceGroup)moleculeGroup.add(vmdSurfaceGroup);
    atomMeshes.length=0;
    cellWireframe=null;
    rotAxisLineMesh=null;
    CX=0;CY=0;CZ=0;
    MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
    if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    var segs=sphereSegs();
    // The pools are caches: an entry outlives the meshes that drew with it, so
    // a rebuild after an edit or a colour toggle does not re-create every
    // material. Scalar colouring is the one case where the number of entries
    // grows with the data — every distinct scalar is its own colour, so a
    // 2175-atom structure carrying a B-factor per atom reaches a few hundred
    // (measured: 524) — so the material pool is dropped once it passes a bound
    // no real structure reaches. Only the cache is emptied: the meshes on
    // screen hold their own references, so nothing about the picture changes.
    if(Object.keys(matPool).length>8192)matPool={};
    MD.atoms.forEach(function(a,i){
        a.index=i;
        var r=atomRadius(a);
        var drawn=r>0.0005;
        // 标量配色是视图状态: 在这里按原子现算, a.color 保持元素颜色不变。
        var col=betaColorOf(a)||a.color;
        var occ=a.occupancy!=null?a.occupancy:1;
        // Geometry and material come from the pools: atoms of the same radius
        // (and colour/occupancy) draw with the very same objects, so a
        // structure of any size costs a handful of GPU buffers and materials
        // instead of one of each per atom.
        var mesh=new THREE.Mesh(sphereGeoFor(drawn?r:0.001,segs),atomMaterial(col,occ,0));
        mesh.position.set(a.x-CX,a.y-CY,a.z-CZ);
        mesh.userData={element:a.element,index:i,baseColor:col,occ:occ};
        mesh.visible=drawn;
        moleculeGroup.add(mesh);
        atomMeshes.push(mesh);
    });
    rebuildBondPieces();
    if(CRY)buildCellWireframe();
    highlightSelected();
    needsRender=true;
    if(diffMode)refreshDiff();
    // 分子组成可能变了 (增删/替换原子、导入、换帧、去无序…), 左上角的
    // 分子式/原子数/电子数面板在这里统一刷新 —— 所有改组成的操作收尾都走
    // rebuildScene, 因此不会漏。
    updateMolInfo();
    // 脚本的表面片 (Surf) 也要跟着新的几何重建; 中心刚算好, 放在最后。
    refreshVmdSurfaces();
}
function updateScenePositions(keepCenter){
    if(!keepCenter){
        CX=0;CY=0;CZ=0;
        MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
        if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    }
    atomMeshes.forEach(function(m,i){var a=MD.atoms[i];if(a)m.position.set(a.x-CX,a.y-CY,a.z-CZ)});
    rebuildBondPieces();
    if(CRY)buildCellWireframe();
    highlightSelected();
    hideVmdSurfaces();
    needsRender=true;
}

function updateAtomMeshPositions(){
    // Light per-frame update used during Move Atoms drags: refresh atom
    // sphere positions + selection highlights only (no bond mesh rebuild).
    atomMeshes.forEach(function(m,i){var a=MD.atoms[i];if(a)m.position.set(a.x-CX,a.y-CY,a.z-CZ)});
    highlightSelected();
    hideVmdSurfaces();
}

// 网格位置都是"相对旋转中心"存放的 (a - C), 所以中心一变, 所有网格只需整体
// 平移 -Δ 就能保持与原子坐标一致 —— 不用重建任何几何。
function shiftAllMeshesForCenter(dx,dy,dz){
    if(!dx&&!dy&&!dz)return;
    var i;
    for(i=0;i<atomMeshes.length;i++){
        atomMeshes[i].position.x-=dx;atomMeshes[i].position.y-=dy;atomMeshes[i].position.z-=dz;
    }
    // A bond piece's transform is computed from the atoms *minus the centre*,
    // so re-placing them all against the new centre is the same shift — and it
    // is one pass over the instances, with no geometry involved.
    placeBondPieces();
    // Surface meshes need no shift: their geometry is written relative to the
    // centre, and refreshVmdSurfaces rebuilds them for the new one.
}
// 重新确定旋转中心 (全部原子的质心)。增删原子、移动/旋转原子、导入结构等凡改变
// 几何的操作在收尾时都要调用它 —— 否则旋转会绕着一个已经不是结构中心的点转,
// 分子看起来就会"绕着圈外的某点在打转"。空结构时中心回到原点。
// 已有的网格按中心位移整体平移 (O(原子数+键数) 次赋值), 不重建几何; 拖拽/预览
// 期间中心保持不变, 只在收尾时重算。
function updateRotationCenter(){
    var oldX=CX,oldY=CY,oldZ=CZ;
    var n=MD.atoms.length;
    CX=0;CY=0;CZ=0;
    if(n>0){
        for(var i=0;i<n;i++){CX+=MD.atoms[i].x;CY+=MD.atoms[i].y;CZ+=MD.atoms[i].z}
        CX/=n;CY/=n;CZ/=n;
    }
    shiftAllMeshesForCenter(CX-oldX,CY-oldY,CZ-oldZ);
    if(rotAxisLineMesh)showRotAxisLine();   // 轴指示线也是相对中心画的
    // 中心变了, 表面片也要按新中心重摆 (签名里含中心, 没变时是空操作)。
    refreshVmdSurfaces();
    needsRender=true;
}

// Rebuild every piece from the current bond list. A piece is pure
// description, so building the list is cheap; the instanced mesh that draws
// them is only ever recreated when the count outgrows the capacity it has.
function rebuildBondPieces(){
    bondPool.pieces.length=0;
    for(var bi=0;bi<MD.bonds.length;bi++){currentBondIdx=bi;createBond(MD.bonds[bi])}
    reindexBondPieces();
    placeBondPieces();
}
// bond index -> the piece slots that draw it, so a drag can re-place exactly
// the bonds that moved without walking the whole list.
function reindexBondPieces(){
    bondPieceIdx=[];
    for(var i=0;i<bondPool.pieces.length;i++){
        var bi=bondPool.pieces[i].bondIdx;
        if(!bondPieceIdx[bi])bondPieceIdx[bi]=[];
        bondPieceIdx[bi].push(i);
    }
}
// Place every piece and upload the instance buffers once. This is the whole
// structure's bond refresh (a scene rebuild, a frame switch, a centre change)
// and also one frame of a normal-mode playback — the arithmetic the two used to
// keep apart is now the same routine.
function placeBondPieces(){
    if(bondPool.pieces.length===0){poolClear(bondPool);needsRender=true;return}
    poolEnsure(bondPool,bondPool.pieces.length);
    for(var i=0;i<bondPool.pieces.length;i++)placePiece(bondPool,i);
    poolUpload(bondPool);
    needsRender=true;
}

// Incremental bond refresh for Move Atoms drags: only bonds with at least one
// endpoint in the moving selection are re-placed; every other piece keeps the
// transform it already has (its endpoints did not move, so it is still exact).
// A drag changes no bond's order, so each affected bond's pieces keep their
// slots and only their instance matrices are rewritten — no geometry at all.
// Runs synchronously on every drag frame so bonds never lag behind atoms.
function rebuildBondsTouching(selSet){
    var sel={};
    for(var i=0;i<selSet.length;i++)sel[selSet[i]]=true;
    for(var bi=0;bi<MD.bonds.length;bi++){
        var b=MD.bonds[bi];
        if(!sel[b.atom1]&&!sel[b.atom2])continue;
        var slots=bondPieceIdx[bi];
        if(!slots)continue;
        for(var k=0;k<slots.length;k++)placePiece(bondPool,slots[k]);
    }
    poolUpload(bondPool);
    needsRender=true;
}

// Release-time localized bond rebuild for Move Atoms: recomputeMovedBonds just
// replaced MD.bonds with a re-derived list, so the pieces of bonds in the
// affected region are stale (bonds may have appeared, vanished or changed
// order) while every other piece is still exact (its endpoints never moved).
// Drops the pieces whose OLD bond touched the region (by old index — the list
// has already been re-derived), renumbers the surviving pieces into the new
// bond array (indices shifted because the removed bonds sat in the middle of
// the old list — placePiece() dereferences bondIdx into MD.bonds, so a stale
// index would place the wrong bond), then adds only the region's new bonds.
function rebuildMovedBondMeshes(oldBonds,selSet){
    var sel={};
    for(var i=0;i<selSet.length;i++)sel[selSet[i]]=true;
    var keep=[];
    for(i=0;i<bondPool.pieces.length;i++){
        var p=bondPool.pieces[i];
        var ob=oldBonds[p.bondIdx];
        if(ob&&(sel[ob.atom1]||sel[ob.atom2]))continue;
        keep.push(p);
    }
    if(oldBonds!==MD.bonds){
        var key2idx={};
        for(var bi=0;bi<MD.bonds.length;bi++){var b2=MD.bonds[bi];key2idx[Math.min(b2.atom1,b2.atom2)+'-'+Math.max(b2.atom1,b2.atom2)]=bi}
        for(i=0;i<keep.length;i++){
            var ob2=oldBonds[keep[i].bondIdx];
            if(!ob2)continue;
            keep[i].bondIdx=key2idx[Math.min(ob2.atom1,ob2.atom2)+'-'+Math.max(ob2.atom1,ob2.atom2)];
        }
    }
    bondPool.pieces.length=0;
    for(i=0;i<keep.length;i++)bondPool.pieces.push(keep[i]);
    for(bi=0;bi<MD.bonds.length;bi++){
        var b=MD.bonds[bi];
        if(!sel[b.atom1]&&!sel[b.atom2])continue;
        currentBondIdx=bi;createBond(b);
    }
    reindexBondPieces();
    placeBondPieces();
}

// Place one piece from the current atom coordinates. The piece's position,
// orientation and size all follow the two atoms it spans, so this is the single
// source of bond geometry: building the scene, dragging an atom, streaming an
// optimisation and playing a normal mode all end up here.
function placePiece(pool,slot){
    var p=pool.pieces[slot];
    var b=MD.bonds[p.bondIdx];
    if(!b||p.bondIdx<0||p.bondIdx>=MD.bonds.length){poolHide(pool,slot);return}
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2){poolHide(pool,slot);return}
    var a2x=a2.x,a2y=a2.y,a2z=a2.z;
    if(b.shift&&CRY){
        var lv=CRY.latticeVectors;
        a2x+=b.shift[0]*lv[0][0]+b.shift[1]*lv[1][0]+b.shift[2]*lv[2][0];
        a2y+=b.shift[0]*lv[0][1]+b.shift[1]*lv[1][1]+b.shift[2]*lv[2][1];
        a2z+=b.shift[0]*lv[0][2]+b.shift[1]*lv[1][2]+b.shift[2]*lv[2][2];
    }
    var sx=a1.x-CX,sy=a1.y-CY,sz=a1.z-CZ;
    var ex=a2x-CX,ey=a2y-CY,ez=a2z-CZ;
    var dx=ex-sx,dy=ey-sy,dz=ez-sz;
    var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(l<0.001){poolHide(pool,slot);return}
    var mpx=(sx+ex)*0.5,mpy=(sy+ey)*0.5,mpz=(sz+ez)*0.5;
    var px=0,py=0,pz=0;
    if(p.perpMag){
        // The extra line of a double or triple bond sits in the plane
        // perpendicular to it: cross the bond with y — or with x for a bond
        // that is itself nearly parallel to y, where that cross product would
        // vanish and the direction would be lost.
        var upx=0,upy=1,upz=0;
        if(Math.abs(dy)>=0.99){upx=1;upy=0;upz=0}
        var qx=dy*upz-dz*upy,qy=dz*upx-dx*upz,qz=dx*upy-dy*upx;
        var ql=Math.sqrt(qx*qx+qy*qy+qz*qz);
        if(ql>0.0001){qx/=ql;qy/=ql;qz/=ql}
        var off=p.perpMag*p.perpSign;
        px=qx*off;py=qy*off;pz=qz*off;
    }
    if(p.isDashed){
        // The dash window is stored as a fraction of the bond, so the dashes
        // follow the bond when it stretches — an aromatic bond that a drag
        // lengthens keeps its pattern instead of leaving a gap.
        var t0=p.dashT0*l,t1=p.dashT1*l;
        var dl=t1-t0;
        if(dl<0.001){poolHide(pool,slot);return}
        var dirx=dx/l,diry=dy/l,dirz=dz/l;
        var dsx=sx+dirx*t0+px,dsy=sy+diry*t0+py,dsz=sz+dirz*t0+pz;
        var dex=sx+dirx*t1+px,dey=sy+diry*t1+py,dez=sz+dirz*t1+pz;
        poolPut(pool,slot,(dsx+dex)*0.5,(dsy+dey)*0.5,(dsz+dez)*0.5,dx,dy,dz,dl,p.r,p.col);
        return;
    }
    var msx,msy,msz,mex,mey,mez;
    if(p.half===0){msx=sx+px;msy=sy+py;msz=sz+pz;mex=mpx+px;mey=mpy+py;mez=mpz+pz}
    else{msx=mpx+px;msy=mpy+py;msz=mpz+pz;mex=ex+px;mey=ey+py;mez=ez+pz}
    var hdx=mex-msx,hdy=mey-msy,hdz=mez-msz;
    var hl=Math.sqrt(hdx*hdx+hdy*hdy+hdz*hdz);
    if(hl<0.001){poolHide(pool,slot);return}
    poolPut(pool,slot,(msx+mex)*0.5,(msy+mey)*0.5,(msz+mez)*0.5,hdx,hdy,hdz,hl,p.r,p.col);
}

// One drawn cylinder of one bond. Nothing is placed here: a piece records
// *what* to draw (radius, colour, which half, the aromatic dash window, the
// multiple-bond offset) and placePiece() reads its transform off the atoms
// every time, so an edit or an animation can never leave a piece behind.
function addBondPiece(r,col,isDashed,dashT0,dashT1,perpSign,perpMag,half){
    bondPool.pieces.push({bondIdx:currentBondIdx,r:r,col:col,isDashed:!!isDashed,
        dashT0:dashT0||0,dashT1:dashT1||0,perpSign:perpSign||0,perpMag:perpMag||0,half:half||0});
}
// The dashed line of an aromatic bond: six short cylinders along the bond,
// shifted onto one side by the same kind of offset a double bond uses.
function addDashPieces(r,col,off){
    var dashes=6,seg=1/dashes,gap=seg*0.35,dashLen=seg-gap;
    for(var k=0;k<dashes;k++){
        var t0=k*seg+gap*0.5;
        var t1=t0+dashLen;
        if(t1>1)t1=1;
        addBondPiece(r,col,true,t0,t1,1,off,0);
    }
}

// Declare the pieces one bond is drawn with. The order thresholds are the same
// the viewer has always used: one line below 1.25, an aromatic line with its
// dashed partner below 1.75, two lines below 2.5, three below 3.5.
function createBond(b){
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2)return;
    // A VMD script sets the bond radius per bond (the representation that
    // covers both endpoints); 0 means the script does not draw it at all.
    var br=typeof b.vmdR==='number'?b.vmdR:0.12;
    if(br<=0.0005)return;
    var ord=b.order||1;
    // Multiple-bond line offsets keep their proportions to the bond radius.
    var offK=br/0.12;
    // 键的两半跟随它连的原子: 标量 (Beta) 配色下原子球是按标量画视图色, 这里
    // 让半键读同一个颜色, 一根键两半就与两端原子一致。betaColorOf 在配色关闭
    // (或原子没有标量) 时返回 null, 退化成元素颜色, 非 Beta 文件的行为不变。
    var c1=new THREE.Color(betaColorOf(a1)||a1.color),c2=new THREE.Color(betaColorOf(a2)||a2.color);
    if(ord<1.25){
        addBondPiece(br,c1,false,0,0,0,0,0);addBondPiece(br,c2,false,0,0,0,0,1);
    }else if(ord<1.75){
        var off=0.10*offK;
        addBondPiece(br,c1,false,0,0,0,0,0);addBondPiece(br,c2,false,0,0,0,0,1);
        addDashPieces(br*0.7,c1,off);
    }else if(ord<2.5){
        var off2=0.12*offK;
        addBondPiece(br*0.6,c1,false,0,0,1,off2,0);addBondPiece(br*0.6,c2,false,0,0,1,off2,1);
        addBondPiece(br*0.6,c1,false,0,0,-1,off2,0);addBondPiece(br*0.6,c2,false,0,0,-1,off2,1);
    }else if(ord<3.5){
        var off3=0.15*offK;
        addBondPiece(br*0.45,c1,false,0,0,0,0,0);addBondPiece(br*0.45,c2,false,0,0,0,0,1);
        addBondPiece(br*0.45,c1,false,0,0,1,off3,0);addBondPiece(br*0.45,c2,false,0,0,1,off3,1);
        addBondPiece(br*0.45,c1,false,0,0,-1,off3,0);addBondPiece(br*0.45,c2,false,0,0,-1,off3,1);
    }else{
        addBondPiece(br,c1,false,0,0,0,0,0);addBondPiece(br,c2,false,0,0,0,0,1);
    }
}

var isRot=false,isPan=false,prevM={x:0,y:0},panX=0,panY=0,camDist=10;
var keyState={ArrowUp:false,ArrowDown:false,ArrowLeft:false,ArrowRight:false};
var keyRotSpeed=0.03;
var rotQuat=new THREE.Quaternion();
var currentMode='view';
var selectedAtoms=[];
var originalCoords=null;
var modalCallback=null;

// ===== Move Atoms mode interactions =====
// In move-atoms mode: plain left-drag / arrow keys rotate the whole view,
// plain right-drag pans the whole view (cursor-locked); Ctrl+left-drag /
// Ctrl+arrow keys rotate the selected fragment, Ctrl+right-drag translates
// the fragment with the cursor locked onto the grabbed point.
var moveDragActive=false;      // a Ctrl+right-drag fragment translate in progress
var moveDragOrig=null;         // {idx:{x,y,z}} positions at drag start
var moveDragMoved=false;       // whether the drag produced any displacement
var moveDragPlaneZ=0;          // world-Z of the movement plane (screen-parallel)
var moveDragStartWorld=null;   // reference point on the movement plane at drag start
var fragRotActive=false;       // a Ctrl+left-drag fragment rotation in progress
var fragRotMoved=false;        // whether the rotation drag produced any change
var fragRotCx=0,fragRotCy=0,fragRotCz=0; // rotation pivot (absolute file coords)
var moveKeyActive=false;       // a Ctrl+arrow fragment rotation session in progress
var moveKeyLatch=false;        // block fragment-key restart until arrows released
var keyRotLatch=false;         // block view rotation until arrows released (after a fragment session ends mid-hold)
var ctrlDown=false;            // live Ctrl modifier state
var panGrabActive=false;       // view pan uses cursor locking
var panGrabOffset={x:0,y:0};   // grabbed point world xy minus pan offset
var panRefZ=0;                 // world z of the pan reference plane
// Diff mode pans each viewport independently, so the right side needs its own
// grab bookkeeping. panGrabSide records which viewport the active grab belongs
// to (matters because the two sides use different cameras and pan offsets).
var panGrabSide='left';
var diffPanGrabOffset={x:0,y:0};
var diffPanRefZ=0;

function recomputeMovedBonds(selOverride){
    // Partial bond re-derivation after a fragment move/rotate session ends.
    //
    // Only the selected atoms changed position, so connectivity — and with it
    // every bond order — can only differ inside the connected components that
    // contain them. Those components are re-derived whole, with the same
    // perception module the host used to load the structure, and everything
    // outside them is kept untouched.
    //
    // Whole components, not a neighbourhood: the module reasons about a
    // component as a unit (its valence pass and its aromatic restoration both
    // work on the component), so handing it a fragment would ask it a
    // different question than the load did and could answer differently. Taking
    // the component whole is also what makes a move reproducible — the bonds
    // that come out are the ones a fresh load of the moved coordinates gives,
    // because the same code sees the same graph.
    //
    // Returns the atoms whose bonds were re-derived (the caller rebuilds only
    // their meshes), an empty array when nothing changed, or null when the
    // caller should fall back to a full scene refresh.
    //
    // selOverride: the atoms that moved, when they are not the current
    // selection (a molecular-mechanics run moves the atoms its own request
    // named, which the user may have re-selected in the meantime).
    var selList=selOverride||selectedAtoms;
    if(selList.length===0)return [];
    var sel={};
    for(var si=0;si<selList.length;si++)sel[selList[si]]=true;
    var atoms=MD.atoms,n=atoms.length;
    if(MD.hasExplicitBonds){
        // Authored bonds keep their orders and are never re-derived — but a
        // bond whose endpoints were pulled beyond the bonding cutoff can no
        // longer be a bond. Only bonds touching the selection can have changed:
        // every other pair kept its distance.
        var keptB=[];
        for(var bx=0;bx<MD.bonds.length;bx++){
            var xb=MD.bonds[bx];
            if(sel[xb.atom1]||sel[xb.atom2]){
                var xu=atoms[xb.atom1],xv=atoms[xb.atom2];
                if(xu&&xv){
                    var xdx=xu.x-xv.x,xdy=xu.y-xv.y,xdz=xu.z-xv.z;
                    if(Math.sqrt(xdx*xdx+xdy*xdy+xdz*xdz)>BP.pairCutoff(xu.element,xv.element))continue;
                }
            }
            keptB.push(xb);
        }
        if(keptB.length===MD.bonds.length)return [];
        MD.bonds=keptB;
        var selArr=[];for(var sk in sel)selArr.push(+sk);
        return selArr;
    }

    // Cheap check first. A move only changes the distances from the moved atoms
    // outward, so only the bonds touching them can have changed geometry — and
    // a bond whose current order is still the one its own length supports
    // cannot come out different. If every such bond is still supported and no
    // new contact has appeared, there is nothing to re-derive and the component
    // walk is skipped altogether; that is what keeps a small drag cheap on a
    // large structure, while a drag that really changes the chemistry still
    // falls through to the full re-derivation below.
    var adjOf={};
    for(var ab=0;ab<MD.bonds.length;ab++){
        var abo=MD.bonds[ab];
        (adjOf[abo.atom1]||(adjOf[abo.atom1]={}))[abo.atom2]=abo.order;
        (adjOf[abo.atom2]||(adjOf[abo.atom2]={}))[abo.atom1]=abo.order;
    }
    var stable=true;
    for(si=0;si<selList.length&&stable;si++){
        var s0=selList[si];
        var nb0=adjOf[s0];
        if(!nb0)continue;
        var a0=atoms[s0];
        if(!a0)continue;
        for(var j0 in nb0){
            var jj=+j0;
            var aj=atoms[jj];
            if(!aj)continue;
            var dx0=a0.x-aj.x,dy0=a0.y-aj.y,dz0=a0.z-aj.z;
            var d0=Math.sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
            if(BP.supportedOrder(a0.element,aj.element,d0)!==nb0[j0]){stable=false;break}
        }
    }

    // Seeds: the moved atoms plus every atom that has just come within bonding
    // range of one. The sweep is the same spatial hash and the same per-element
    // search radii the detector uses, so a contact dragged onto a previously
    // unrelated atom is still found.
    var newContact=false;
    var inA={};
    for(si=0;si<selList.length;si++)inA[selList[si]]=true;
    var CELL=3;
    var normEls=new Array(n);var distinct=[];var seenEl={};
    for(var i=0;i<n;i++){var rawEl=atoms[i].element||'';var ne=rawEl.charAt(0).toUpperCase()+rawEl.slice(1).toLowerCase();normEls[i]=ne;
        if(!seenEl[ne]){seenEl[ne]=1;distinct.push(ne)}}
    var rangeFor={};
    for(var di=0;di<distinct.length;di++){var e1=distinct[di];var mc=0;
        for(var dj=0;dj<distinct.length;dj++){var c2=BP.pairCutoff(e1,distinct[dj]);if(c2>mc)mc=c2}
        rangeFor[e1]=Math.floor(mc/CELL)+1}
    var grid=new Map();
    for(i=0;i<n;i++){var cx=Math.floor(atoms[i].x/CELL),cy=Math.floor(atoms[i].y/CELL),cz=Math.floor(atoms[i].z/CELL);
        var key=(cx*73856093)^(cy*19349663)^(cz*83492791);
        var bkt=grid.get(key);if(!bkt){bkt=[];grid.set(key,bkt)}bkt.push(i)}
    for(si=0;si<selList.length;si++){
        var s=selList[si];if(!atoms[s])continue;
        var range=rangeFor[normEls[s]];
        var csx=Math.floor(atoms[s].x/CELL),csy=Math.floor(atoms[s].y/CELL),csz=Math.floor(atoms[s].z/CELL);
        for(var dx=-range;dx<=range;dx++)for(var dy=-range;dy<=range;dy++)for(var dz=-range;dz<=range;dz++){
            var gk=((csx+dx)*73856093)^((csy+dy)*19349663)^((csz+dz)*83492791);
            var gb=grid.get(gk);if(!gb)continue;
            for(var q=0;q<gb.length;q++){var j=gb[q];if(j===s)continue;
                var ddx=atoms[s].x-atoms[j].x,ddy=atoms[s].y-atoms[j].y,ddz=atoms[s].z-atoms[j].z;
                if(Math.sqrt(ddx*ddx+ddy*ddy+ddz*ddz)<=BP.pairCutoff(atoms[s].element,atoms[j].element)){
                    inA[j]=1;
                    // a pair that had no bond and is now inside the cutoff is a
                    // bond that has to be considered
                    if(!(adjOf[s]&&adjOf[s][j]))newContact=true}}
        }
    }

    // Nothing moved a bond off the order its own length supports, and nothing
    // new came into range: the re-derivation below would reproduce exactly what
    // is already there. This is the common outcome of a small adjustment.
    if(stable&&!newContact)return [];

    // Grow the seeds to whole components along the current bond graph.
    var adj={};
    for(var bi=0;bi<MD.bonds.length;bi++){var b=MD.bonds[bi];
        (adj[b.atom1]||(adj[b.atom1]=[])).push(b.atom2);
        (adj[b.atom2]||(adj[b.atom2]=[])).push(b.atom1)}
    var queue=[];for(var ik in inA)queue.push(+ik);
    while(queue.length){var ii=queue.pop();var nbA=adj[ii];if(!nbA)continue;
        for(var a2=0;a2<nbA.length;a2++){var jj=nbA[a2];if(!inA[jj]){inA[jj]=1;queue.push(jj)}}}
    var comp=[];for(var ck in inA)comp.push(+ck);
    comp.sort(function(a,b){return a-b});
    if(comp.length===0)return [];

    // Re-derive the component in its own index space. The molecule's total
    // charge describes the whole structure, so it is only meaningful when the
    // component is the whole structure; a component of a multi-fragment system
    // gets the neutral assumption, exactly as a bare coordinate file would.
    var subAtoms=new Array(comp.length);
    for(var ci=0;ci<comp.length;ci++)subAtoms[ci]=atoms[comp[ci]];
    var opt=(comp.length===n&&typeof MD.charge==='number')?{charge:MD.charge}:undefined;
    var res=BP.perceiveBonds(subAtoms,opt);

    // Splice: a bond with both ends inside the component is replaced by the
    // freshly derived one; everything else is kept as it is. What the derived
    // bond carries is the chemistry (its endpoints and its order) — anything the
    // *file's own recipe* put on a bond for drawing it, a VMD script's per-bond
    // radius above all, is not the detector's business and is carried across for
    // a pair that already had it, so re-deriving an order cannot silently
    // re-thin a bond the script asked to be drawn thicker (or hide one the
    // script said not to draw at all).
    var drawProp={};
    for(bi=0;bi<MD.bonds.length;bi++){
        b=MD.bonds[bi];
        if(b.vmdR!==undefined)drawProp[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.vmdR;
    }
    var compIdx={};for(ci=0;ci<comp.length;ci++)compIdx[comp[ci]]=1;
    var out=[];
    for(bi=0;bi<MD.bonds.length;bi++){b=MD.bonds[bi];
        if(!compIdx[b.atom1]&&!compIdx[b.atom2])out.push(b)}
    for(var ri=0;ri<res.bonds.length;ri++){var rb=res.bonds[ri];
        var na=comp[rb.atom1],nx=comp[rb.atom2];
        var nbb={atom1:na,atom2:nx,order:rb.order};
        var kp=drawProp[Math.min(na,nx)+'-'+Math.max(na,nx)];
        if(kp!==undefined)nbb.vmdR=kp;
        out.push(nbb)}

    // An unchanged result means there is nothing to rebuild.
    var before={};
    for(bi=0;bi<MD.bonds.length;bi++){b=MD.bonds[bi];before[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.order}
    var after={};
    for(bi=0;bi<out.length;bi++){b=out[bi];after[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.order}
    var keysA=Object.keys(after);
    var same=keysA.length===Object.keys(before).length;
    if(same)for(var ka=0;ka<keysA.length;ka++){
        if(before[keysA[ka]]!==after[keysA[ka]]){same=false;break}}
    if(same)return [];

    MD.bonds=out;
    return comp;
}
// Shared release path for every Move Atoms session end (mouse drag, fragment
// rotation drag, arrow-key session, forced termination): re-derive the bonds
// for the moved region, then rebuild only the affected region's bond meshes.
// A full scene rebuild runs only when the region-scoped re-detection fell
// back to the full detector (null region).
function refreshMovedBondsAndMeshes(selOverride){
    var oldBonds=MD.bonds;
    var region=recomputeMovedBonds(selOverride);
    if(region){if(region.length)rebuildMovedBondMeshes(oldBonds,region)}
    else updateScenePositions(true);
    // 原子被移动过, 质心随之改变 —— 放在区域刷新之后: 上面的网格都还是按旧中心
    // 摆放的, 这里再整体平移一次, 随后旋转就绕新的结构中心进行。
    updateRotationCenter();
}
function endMoveDrag(){
    if(!moveDragActive)return;
    moveDragActive=false;
    moveDragOrig=null;
    // A Ctrl+click without any displacement pushes a redundant snapshot — drop it.
    if(!moveDragMoved&&undoStack.length>0){undoStack.pop();updateUndoBtn()}
    // Re-detect connectivity once on release (not per frame during the drag).
    if(moveDragMoved){provMark('translate',1);refreshMovedBondsAndMeshes();}
    else refreshVmdSurfaces();
    document.body.style.cursor='';
    if(currentMode==='moveAtoms')modeInfoEl.textContent=MODE_INFO.moveAtoms;
}
function endFragRot(){
    if(!fragRotActive)return;
    fragRotActive=false;
    // A Ctrl+click without any rotation pushes a redundant snapshot — drop it.
    if(!fragRotMoved&&undoStack.length>0){undoStack.pop();updateUndoBtn()}
    // Re-detect connectivity once on release (not per frame during the drag).
    if(fragRotMoved){provMark('rotate',1);refreshMovedBondsAndMeshes();}
    else refreshVmdSurfaces();
    document.body.style.cursor='';
    if(currentMode==='moveAtoms')modeInfoEl.textContent=MODE_INFO.moveAtoms;
}
function stopMoveSession(){
    // Force-terminate any in-progress move/rotate (mode switch, frame switch,
    // undo, delete...). An Alt-drag writes the same coordinates a move session
    // writes, and every entry that replaces the structure already passes
    // through here, so it is ended at this one door too. The latch prevents
    // arrow keys (still held down) from instantly starting a new session; it
    // clears once all keys are released.
    if(bakeActive)endBakeSession();
    if(moveDragActive||moveKeyActive||fragRotActive){
        // The geometry these sessions produced stays; the session is only
        // ended, so what it did is still part of the structure that gets saved.
        if(moveDragActive&&moveDragMoved)provMark('translate',1);
        if((fragRotActive||moveKeyActive)&&fragRotMoved)provMark('rotate',1);
        moveDragActive=false;moveKeyActive=false;fragRotActive=false;
        moveKeyLatch=true;
        moveDragOrig=null;
        // Atoms may have moved — re-detect connectivity once.
        refreshMovedBondsAndMeshes();
        document.body.style.cursor='';
    }
}
function applyMoveWorldDelta(dw){
    // dw: world-space displacement. Convert to molecule-local space (the view
    // rotation is applied to moleculeGroup, not the camera) so the group moves
    // exactly along the drag direction on screen regardless of orientation.
    var dl=dw.clone().applyQuaternion(rotQuat.clone().conjugate());
    selectedAtoms.forEach(function(i){
        var a=MD.atoms[i];if(!a)return;
        if(moveDragOrig&&moveDragOrig[i]){
            var o=moveDragOrig[i];
            a.x=o.x+dl.x;a.y=o.y+dl.y;a.z=o.z+dl.z;
        }else{
            a.x+=dl.x;a.y+=dl.y;a.z+=dl.z;
        }
    });
    updateAtomMeshPositions();
    rebuildBondsTouching(selectedAtoms);
}
function applyFragRotDelta(ax,ay){
    // Incremental rotation of the selected fragment around its (fixed) center,
    // in WORLD space: ax = angle around world X, ay = angle around world Y.
    // Same axis/angle convention as the whole-molecule rotation, so the
    // fragment turns in the same direction the view would. The world-space
    // quaternion is converted to molecule-local space via conjugation.
    if(!ax&&!ay)return;
    var Q=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),ay);
    if(ax)Q.premultiply(new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),ax));
    var qComb=rotQuat.clone().conjugate().multiply(Q).multiply(rotQuat);
    var v=new THREE.Vector3();
    selectedAtoms.forEach(function(i){
        var a=MD.atoms[i];if(!a)return;
        v.set(a.x-fragRotCx,a.y-fragRotCy,a.z-fragRotCz).applyQuaternion(qComb);
        a.x=fragRotCx+v.x;a.y=fragRotCy+v.y;a.z=fragRotCz+v.z;
    });
    updateAtomMeshPositions();
    rebuildBondsTouching(selectedAtoms);
}
function screenPointOnPlane(e,planeZ){
    // Intersect the mouse ray with the z=planeZ plane (camera sits on the
    // Z axis looking at the origin, so this plane is parallel to screen).
    // planeZ must stay in front of the camera, otherwise the hit is a
    // mirrored point behind the camera (inverted drag direction).
    var rect=canvas.getBoundingClientRect();
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var dz=raycaster.ray.direction.z;
    if(Math.abs(dz)<1e-6)return null;
    var t=(planeZ-raycaster.ray.origin.z)/dz;
    return raycaster.ray.origin.clone().add(raycaster.ray.direction.clone().multiplyScalar(t));
}
function screenToPlanePoint(e){
    return screenPointOnPlane(e,moveDragPlaneZ);
}
function atomWorldZ(a){
    // World Z of an atom: local (a - center) rotated by the view quaternion.
    // pan offsets only shift x/y, so z is unaffected by panning.
    return new THREE.Vector3(a.x-CX,a.y-CY,a.z-CZ).applyQuaternion(rotQuat).z;
}
function diffAtomWorldZ(a){
    // Same as atomWorldZ, but for the right (comparison) structure, which has
    // its own center and its own rotation quaternion.
    return new THREE.Vector3(a.x-diffCX,a.y-diffCY,a.z-diffCZ).applyQuaternion(diffRotQuat).z;
}
function diffViewportSide(e){
    // Which half of the canvas the pointer is over: 'left' | 'right'.
    var rect=canvas.getBoundingClientRect();
    return (e.clientX-rect.left)<Math.floor(rect.width/2)?'left':'right';
}
function configureDiffCamera(side,dist){
    // Point the shared camera at one diff viewport with the exact aspect and
    // distance used by the render pass for that side, so raycasting from the
    // cursor matches what the user sees on that half of the canvas.
    var rect=canvas.getBoundingClientRect();
    var w=rect.width,h=rect.height,halfW=Math.floor(w/2);
    camera.aspect=(side==='right'?(w-halfW):halfW)/h;
    camera.updateProjectionMatrix();
    camera.position.set(0,0,dist);
    camera.updateMatrixWorld();
}
function sideScreenToPlanePoint(e,side,planeZ){
    // Like screenPointOnPlane, but maps the cursor into the NDC of the given
    // diff viewport instead of the whole canvas, and configures the camera to
    // match. planeZ must stay in front of the camera.
    var rect=canvas.getBoundingClientRect();
    var w=rect.width,h=rect.height,halfW=Math.floor(w/2);
    var localX=e.clientX-rect.left;
    mouse.x=(side==='right'?((localX-halfW)/(w-halfW)):(localX/halfW))*2-1;
    mouse.y=-((e.clientY-rect.top)/h)*2+1;
    configureDiffCamera(side,side==='right'?diffCamDist:camDist);
    raycaster.setFromCamera(mouse,camera);
    var dz=raycaster.ray.direction.z;
    if(Math.abs(dz)<1e-6)return null;
    var t=(planeZ-raycaster.ray.origin.z)/dz;
    return raycaster.ray.origin.clone().add(raycaster.ray.direction.clone().multiplyScalar(t));
}
function initPanGrab(e){
    // Cursor-locked view panning: grab the atom under the cursor (or the
    // molecule-center plane when clicking empty space) so the grabbed point
    // sticks to the cursor while panning — 1:1 cursor speed.
    // In diff mode each viewport pans independently, so the grab is scoped to
    // the side the pointer is over and applied to that side's pan offset.
    panGrabActive=false;
    if(diffMode){
        var side=diffViewportSide(e);
        panGrabSide=side;
        var idxD=getClickedAtom(e); // also sets diffActiveSide + camera for that side
        if(side==='right'){
            var refZR=(idxD>=0&&diffData&&diffData.atoms&&diffData.atoms[idxD])?diffAtomWorldZ(diffData.atoms[idxD]):0;
            var hitR=sideScreenToPlanePoint(e,'right',refZR);
            if(!hitR)return;
            diffPanGrabOffset={x:hitR.x-diffPanX,y:hitR.y-diffPanY};
            diffPanRefZ=refZR;
        }else{
            var refZL=(idxD>=0&&MD.atoms[idxD])?atomWorldZ(MD.atoms[idxD]):0;
            var hitL=sideScreenToPlanePoint(e,'left',refZL);
            if(!hitL)return;
            panGrabOffset={x:hitL.x-panX,y:hitL.y-panY};
            panRefZ=refZL;
        }
        panGrabActive=true;
        return;
    }
    var idx=getClickedAtom(e);
    var refZ=0;
    if(idx>=0&&MD.atoms[idx])refZ=atomWorldZ(MD.atoms[idx]);
    var hit=screenPointOnPlane(e,refZ);
    if(!hit)return;
    panGrabOffset={x:hit.x-panX,y:hit.y-panY};
    panRefZ=refZ;
    panGrabSide='left';
    panGrabActive=true;
}

if(CRY){rebuildCrystal();}
rebuildScene();

var spInit=document.getElementById('select-panel');
if(spInit)spInit.style.display='none';
container.tabIndex=0;
container.focus();

if(CRY){
    var crystalPanel=document.getElementById('crystal-panel');
    if(crystalPanel)crystalPanel.style.display='block';
    layoutPanels();
    ['bnd-a-min','bnd-a-max','bnd-b-min','bnd-b-max','bnd-c-min','bnd-c-max'].forEach(function(id){
        var inp=document.getElementById(id);
        if(inp){
            inp.addEventListener('change',function(){
                var v=parseFloat(inp.value);
                if(isNaN(v))v=0;
                v=Math.max(-5,Math.min(5,v));
                inp.value=v;
                if(id==='bnd-a-min')boundary.aMin=v;
                else if(id==='bnd-a-max')boundary.aMax=v;
                else if(id==='bnd-b-min')boundary.bMin=v;
                else if(id==='bnd-b-max')boundary.bMax=v;
                else if(id==='bnd-c-min')boundary.cMin=v;
                else if(id==='bnd-c-max')boundary.cMax=v;
                if(boundary.aMax<boundary.aMin)boundary.aMax=boundary.aMin;
                if(boundary.bMax<boundary.bMin)boundary.bMax=boundary.bMin;
                if(boundary.cMax<boundary.cMin)boundary.cMax=boundary.cMin;
                rebuildCrystal();
                rebuildScene();
            });
        }
    });
    var rmDisorderBtn=document.getElementById('bnd-remove-disorder');
    if(rmDisorderBtn){
        rmDisorderBtn.addEventListener('click',function(){
            if(!CRY)return;
            var removedCount=0;
            CRY.baseAtoms.forEach(function(ba){
                var occ=ba.occupancy!=null?ba.occupancy:1;
                if(occ<0.5)removedCount++;
            });
            if(removedCount===0){
                showModal('<h3>No Disorder</h3><div class="current-val">No atoms with occupancy &lt; 0.5 found.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
                document.getElementById('m-cancel').addEventListener('click',hideModal);
                return;
            }
            showModal('<h3>Remove Disorder</h3>'+
                '<div class="current-val">Remove '+removedCount+' atom(s) with occupancy &lt; 0.5?<br>Atoms with occupancy &gt; 0.5 will be set to 1.0.<br>Atoms with occupancy = 0.5 will be kept unchanged.</div>'+
                '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok mbtn-danger" id="m-ok">Remove</button></div>',null);
            document.getElementById('m-ok').addEventListener('click',function(){hideModal();removeDisorder()});
            document.getElementById('m-cancel').addEventListener('click',hideModal);
        });
    }
}

// Axes indicator is shown for ALL structures: crystal files display the
// lattice a/b/c axes, molecular files display the Cartesian X/Y/Z axes of
// the initial coordinate frame.
var axesInitEl=document.getElementById('axes-indicator');
if(axesInitEl)axesInitEl.style.display='block';
updateAxesIndicator();

// Camera distance that frames the current structure: the farthest atom from the
// rotation centre, with margin. Recomputed when the structure is swapped.
function computeInitialCamDist(){
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    return maxD*2.5+5;
}
var initCam=computeInitialCamDist();
camera.position.set(0,0,initCam);camera.lookAt(0,0,0);
camDist=initCam;

// ===== Molecular mechanics optimisation =====
// The panel carries the force field's own parameter list (built host-side from
// the module that runs the optimisation, so the two can never disagree), the
// precision presets, the atoms a run would move and the buttons that steer it.
// Selection is the viewer's usual one — a click here, or a selection carried in
// from Select Atoms / Box Select / Replace Atom / Add H / Rotate Group / Move
// Atoms — and only those atoms move.
//
// A run is a *session* the host steps in chunks, so the panel watches it work:
// the host posts mmProgress / mmPaused / mmCancelled / mmResult, and every
// progress chunk's coordinates are written straight into the structure (spheres
// and only the bonds that touch them — the cheap Move Atoms path) while the
// energy and force curves grow from the history that travels with it. **Pause**
// stops the host's pump and keeps the session, so **Resume** continues from
// exactly that state and **Restore initial** writes this run's own starting
// geometry back; leaving the mode — another tool, a frame, a molecule, an undo,
// an import, a comparison — **cancels** the run, because an optimisation must
// never keep writing coordinates into a structure the user has moved on from.
// The run id ties every reply to the run that asked for it, so a late answer
// from a replaced run is ignored instead of applied.
var mmLastResult=null;       // last finished run's numbers, shown under the buttons
var mmErrorMsg='';           // last failure, shown in the panel
var mmOptions={};            // the settings chosen in the panel (key -> value)
var mmRun=null;              // the run being watched (null = none), see mmStartRun
var mmRunSeq=0;              // run id source; the host echoes the id back
var MM_HISTORY_CAP=4000;     // curve points kept; past this, every other one

function mmPanelEl(){return document.getElementById('mm-panel')}
function mmPanelVisible(){var el=mmPanelEl();return !!(el&&el.classList.contains('show'))}
// The controls the panel offers come from the host's field description, so the
// choices and the defaults are the optimiser's own and cannot drift from it.
function mmControls(){
    var p=MD.mmParams;
    return (p&&p.controls)?p.controls:[];
}
// The value a control currently holds: what the user picked, or its default.
function mmOptionValue(ctl){
    var v=mmOptions[ctl.key];
    return (v===undefined||v===null)?ctl.defaultValue:v;
}
function mmOptionSnapshot(){
    var out={};
    var ctls=mmControls();
    for(var i=0;i<ctls.length;i++)out[ctls[i].key]=mmOptionValue(ctls[i]);
    return out;
}

// ---- precision presets: the choice lives here, the values in the module ----
function mmPresets(){
    var p=MD.mmParams;
    return (p&&p.presets)?p.presets:[];
}
function mmFindPreset(key){
    var list=mmPresets();
    for(var i=0;i<list.length;i++)if(list[i].key===key)return list[i];
    return null;
}
// Choosing a preset writes the module's values into the controls, so what runs
// is exactly what the selects show.
function mmApplyPreset(key){
    mmOptions.preset=key;
    var p=mmFindPreset(key);
    if(!p)return;
    var ks=Object.keys(p.options||{});
    for(var i=0;i<ks.length;i++)mmOptions[ks[i]]=p.options[ks[i]];
}
// A manual change keeps the preset label only while every other control still
// matches it; otherwise the choice is a custom one.
function mmSyncPresetState(){
    var cur=mmOptions.preset;
    if(cur==='custom')return;
    var p=mmFindPreset(cur);
    if(!p){mmOptions.preset='custom';return}
    var ctls=mmControls();
    for(var i=0;i<ctls.length;i++){
        var k=ctls[i].key;
        if(k==='preset')continue;
        if(p.options[k]===undefined)continue;
        if(String(mmOptionValue(ctls[i]))!==String(p.options[k])){mmOptions.preset='custom';return}
    }
}
// The preset select is the one control a manual change can move *indirectly*,
// so it is updated in place rather than by re-rendering the panel (a rebuild
// would drop the user's focus and redraw every canvas mid-click).
function mmSyncPresetSelect(){
    var el=mmPanelEl();
    if(!el||typeof el.querySelector!=='function')return;
    var sel=el.querySelector('.mm-ctl-sel[data-mmkey="preset"]');
    if(sel&&String(sel.value)!==String(mmOptions.preset))sel.value=String(mmOptions.preset);
}

// ---- the run being watched ----
function mmRunActive(){return !!(mmRun&&(mmRun.state==='running'||mmRun.state==='paused'))}
function mmRunMatches(msg){return !!(mmRun&&msg&&msg.runId===mmRun.id)}
function mmPoints(){return mmRun?mmRun.history:[]}
// The host sends only the points this panel has not seen; a point already on the
// curve is dropped, so one re-sent after a history thinning cannot double up.
function mmAppendPoints(pts){
    if(!mmRun||!Array.isArray(pts))return;
    var i;
    for(i=0;i<pts.length;i++){
        var p=pts[i];
        if(!p||typeof p.step!=='number'||p.step<=mmRun.lastStep)continue;
        mmRun.history.push({step:p.step,energy:p.energy,maxGrad:p.maxGrad,rmsGrad:p.rmsGrad});
        mmRun.lastStep=p.step;
    }
    if(mmRun.history.length>MM_HISTORY_CAP){
        var thinned=[];
        for(i=0;i<mmRun.history.length;i+=2)thinned.push(mmRun.history[i]);
        var last=mmRun.history[mmRun.history.length-1];
        if(thinned[thinned.length-1]!==last)thinned.push(last);
        mmRun.history=thinned;
    }
}
// Write one chunk's coordinates into the structure. The first write takes the
// undo snapshot, so a single Undo returns the structure to before the run; per
// chunk only the moved atoms' spheres and the bond meshes that touch them are
// rebuilt (the Move Atoms drag path), and the full release pass runs once, when
// the run stops.
function mmApplyLiveCoords(coords){
    if(!mmRun||!Array.isArray(coords))return false;
    var idx=mmRun.indices;
    var n=Math.min(idx.length,coords.length);
    // Which of this chunk's coordinates are actually new, decided *before* any
    // of them is written — so the undo snapshot below is taken once, in front of
    // the first real change, and a chunk that moves nothing (a converged run's
    // first report, or a repeated state) leaves no undo step behind.
    var plan=[];
    for(var k=0;k<n;k++){
        var c=coords[k];
        if(!c||c.length<3||!isFinite(c[0])||!isFinite(c[1])||!isFinite(c[2]))continue;
        var i=idx[k];
        if(typeof i!=='number'||i<0||i>=MD.atoms.length)continue;
        var a=MD.atoms[i];
        if(a.x===c[0]&&a.y===c[1]&&a.z===c[2])continue;
        plan.push([i,c]);
    }
    if(plan.length===0)return false;
    if(!mmRun.undoPushed){pushUndo();mmRun.undoPushed=true}
    var moved=[];
    for(k=0;k<plan.length;k++){
        i=plan[k][0];c=plan[k][1];
        a=MD.atoms[i];
        a.x=c[0];a.y=c[1];a.z=c[2];
        moved.push(i);
        if(mmRun.moved.indexOf(i)<0)mmRun.moved.push(i);
        if(!mmRun.movedAll[i]){mmRun.movedAll[i]=true;mmRun.movedCount++}
    }
    updateAtomMeshPositions();
    rebuildBondsTouching(moved);
    return true;
}
// Re-derive the connectivity of the region the run moved — the Move Atoms
// release path — exactly once per stop. This is the one heavyweight step on the
// client side of a run, and it is guarded: it must never be able to take its
// caller (and with it the panel's state and its buttons) down with it. A
// failure here leaves the spheres exactly where the last chunk put them and
// says so, rather than stranding the panel half-updated.
function mmEndLiveRefresh(){
    if(!mmRun||mmRun.moved.length===0)return;
    var moved=mmRun.moved;
    mmRun.moved=[];
    try{
        refreshMovedBondsAndMeshes(moved);
    }catch(err){
        mmErrorMsg='Molecular Mechanics: the moved region could not be re-perceived - '+
            (err&&err.message?err.message:err);
        updateAtomMeshPositions();
    }
}
function mmUpdateProgressText(){
    var el=document.getElementById('mm-prog');
    if(!el)return;
    if(!mmRun||!mmRun.last){el.textContent='';return}
    var L=mmRun.last;
    var what=mmRun.state==='paused'?'paused':(mmRun.state==='running'?'running':(mmRun.restored?'restored':'finished'));
    var txt='Step '+L.step+' ('+what+') - E '+L.energy.toFixed(3)+' kcal/mol';
    if(typeof mmRun.energyBefore==='number')txt+=' (delta '+(L.energy-mmRun.energyBefore).toFixed(3)+')';
    txt+=' \u00b7 max|F| '+L.maxGrad.toFixed(3)+' \u00b7 RMS '+L.rmsGrad.toFixed(3);
    el.textContent=txt;
}

// ---- the live curves (plain 2-D canvas, drawn from the run's history) ----
function mmFmtValue(v){
    var a=Math.abs(v);
    if(a!==0&&(a<0.01||a>=100000))return v.toExponential(1);
    if(a>=1000)return v.toFixed(0);
    if(a>=10)return v.toFixed(1);
    return v.toFixed(2);
}
function mmDrawCurve(canvas,series,useLog){
    var pts=mmPoints();
    if(!canvas||pts.length<1)return;
    var ctx=canvas.getContext('2d');
    if(!ctx)return;
    var i,s;
    var dpr=window.devicePixelRatio||1;
    // The canvas fills the panel (width:100%), so its own box — not the panel's
    // padding box — is the width to draw into. A pixel width pinned by an earlier
    // pass is cleared first, so the measurement is of the current layout and not
    // of the last draw, and the display box is left to CSS: drawing into the
    // padding box instead used to make the curve wider than the element that
    // shows it, clipping its right edge (the end-of-axis step label with it).
    canvas.style.width='';canvas.style.height='';
    var w=canvas.clientWidth||0;
    if(w<1)w=300;
    var h=canvas.clientHeight||78;
    canvas.width=w*dpr;canvas.height=h*dpr;
    ctx.setTransform(dpr,0,0,dpr,0,0);
    ctx.clearRect(0,0,w,h);
    var padL=46,padR=6,padT=15,padB=14;
    var pw=w-padL-padR,ph=h-padT-padB;
    var vals=[];
    for(s=0;s<series.length;s++){
        for(i=0;i<pts.length;i++){
            var v=pts[i][series[s].field];
            if(typeof v==='number'&&isFinite(v))vals.push(v);
        }
    }
    if(vals.length===0||pw<=4||ph<=4){
        ctx.fillStyle='#888';ctx.font='9px sans-serif';ctx.textAlign='left';ctx.textBaseline='top';
        ctx.fillText('no data',padL,padT+10);
        return;
    }
    var vmin=Math.min.apply(null,vals),vmax=Math.max.apply(null,vals);
    var log=!!useLog&&vmin>0;
    var lo,hi;
    if(log){
        lo=Math.log(Math.max(vmin,1e-12));hi=Math.log(Math.max(vmax,1e-12));
    }else{lo=vmin;hi=vmax}
    if(hi-lo<1e-12)hi=lo+1;
    function yOf(v){
        var t=log?Math.log(Math.max(v,1e-12)):v;
        return padT+ph*(hi-t)/(hi-lo);
    }
    ctx.strokeStyle='rgba(255,255,255,0.12)';ctx.lineWidth=1;
    ctx.fillStyle='#888';ctx.font='9px sans-serif';ctx.textAlign='right';ctx.textBaseline='middle';
    for(var g=0;g<=2;g++){
        var y=padT+ph*g/2;
        ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(padL+pw,y);ctx.stroke();
        var tv=log?Math.exp(hi-(hi-lo)*g/2):hi-(hi-lo)*g/2;
        ctx.fillText(mmFmtValue(tv),padL-3,y);
    }
    ctx.textAlign='center';ctx.textBaseline='top';
    ctx.fillText('0',padL,padT+ph+2);
    ctx.fillText(String(pts[pts.length-1].step),padL+pw,padT+ph+2);
    for(s=0;s<series.length;s++){
        var ser=series[s];
        ctx.strokeStyle=ser.color;ctx.lineWidth=1.5;ctx.beginPath();
        var drew=false;
        for(i=0;i<pts.length;i++){
            var vv=pts[i][ser.field];
            if(typeof vv!=='number'||!isFinite(vv))continue;
            if(log&&vv<=0)continue;
            var x=padL+(pts.length<=1?0:pw*i/(pts.length-1));
            var yy=yOf(vv);
            if(!drew){ctx.moveTo(x,yy);drew=true}else ctx.lineTo(x,yy);
        }
        ctx.stroke();
        // A run that has taken no step yet has exactly one point, and one moveTo
        // draws no line at all — the curve would read as "no data" for the whole
        // first chunk even though the run has started and reported its input
        // energy. The single point is marked with a dot instead.
        if(pts.length===1){
            var v1=pts[0][ser.field];
            if(typeof v1==='number'&&isFinite(v1)&&!(log&&v1<=0)){
                ctx.fillStyle=ser.color;
                ctx.fillRect(padL-1,yOf(v1)-1,3,3);
            }
        }
    }
    ctx.textAlign='left';ctx.textBaseline='top';ctx.font='9px sans-serif';
    for(s=0;s<series.length;s++){
        ctx.fillStyle=series[s].color;
        ctx.fillRect(padL+s*92,padT-11,8,2);
        ctx.fillText(series[s].label,padL+s*92+11,padT-13);
    }
}
function mmDrawCurves(){
    if(currentMode!=='mmOpt')return;
    if(!mmRun)return;
    mmDrawCurve(document.getElementById('mm-cv-energy'),
        [{field:'energy',color:'#3794ff',label:'Energy'}],false);
    mmDrawCurve(document.getElementById('mm-cv-force'),
        [{field:'maxGrad',color:'#ff6b6b',label:'max |F|'},{field:'rmsGrad',color:'#74c0fc',label:'rms |F|'}],true);
}

// ---- the small button that reopens the panel after its x ----
// Same idea as the log file's Convergence / Route buttons: closing the panel
// hides it, this brings it back, and its label doubles as the run's state while
// one is in flight (so a closed panel still shows that the run is alive).
function mmUpdateReopen(){
    var el=document.getElementById('mm-reopen');
    if(!el)return;
    var show=(!diffMode&&!CRY&&currentMode==='mmOpt'&&!mmPanelVisible());
    el.classList.toggle('show',show);
    el.classList.toggle('busy',mmRunActive());
    if(mmRunActive()&&mmRun.last){
        el.textContent='\u2697\uFE0F MM '+(mmRun.state==='paused'?'paused':'running')+' - step '+mmRun.last.step;
    }else{
        el.textContent='\u2697\uFE0F MM Optimize';
    }
}
(function(){
    var el=document.getElementById('mm-reopen');
    if(el)el.addEventListener('click',function(){showMmPanel()});
})();

function showMmPanel(){
    var el=mmPanelEl();
    if(!el)return;
    el.classList.add('show');
    mmUpdateReopen();
    refreshMmPanel();
}
function hideMmPanel(){
    var el=mmPanelEl();
    if(!el)return;
    el.classList.remove('show');
    mmUpdateReopen();
    layoutPanels();
}
function mmStopText(reason){
    if(reason==='gradient')return 'converged';
    if(reason==='energy')return 'energy converged';
    if(reason==='time')return 'time budget reached';
    if(reason==='iterations')return 'step budget reached';
    return 'no further gain';
}
function mmResultHTML(r){
    if(!r)return '';
    var delta=r.energyAfter-r.energyBefore;
    var c=r.counts||{};
    var html='<div class="mm-res"><div class="mm-res-head">Last run</div>';
    html+='Energy: '+r.energyBefore.toFixed(3)+' &rarr; '+r.energyAfter.toFixed(3)+' kcal/mol (&Delta; '+delta.toFixed(3)+')<br>';
    html+='Steps: '+r.iterations+' ('+mmStopText(r.stopReason)+') &middot; max|F| '+r.maxGrad.toFixed(3)+' &middot; RMS '+r.rmsGrad.toFixed(3)+'<br>';
    html+='Moved '+r.moved+' atom(s), largest displacement '+r.maxMove.toFixed(3)+' A<br>';
    html+='Terms: '+(c.bonds||0)+' bond, '+(c.angles||0)+' angle, '+(c.torsions||0)+' torsion, '+(c.impropers||0)+' improper, '+(c.nonbonded||0)+' non-bonded';
    html+='</div>';
    return html;
}
// One selection control per settable parameter, with the current choice
// preselected — the default choice is the one the panel opens with.
function mmControlHTML(ctl){
    var cur=mmOptionValue(ctl);
    var html='<div class="mm-ctl"><label class="mm-ctl-label">'+ctl.label+'</label>';
    html+='<select class="mm-ctl-sel" data-mmkey="'+ctl.key+'">';
    for(var i=0;i<ctl.options.length;i++){
        var opt=ctl.options[i];
        var sel=(String(opt.value)===String(cur))?' selected':'';
        html+='<option value="'+String(opt.value)+'"'+sel+'>'+opt.label+'</option>';
    }
    html+='</select></div>';
    return html;
}
// The buttons of the current state. While a run is going there is **one**
// toggle button that reads Pause and becomes Resume when the run is paused —
// the same single-button pattern the vibration panel uses, so the two states
// sit in the same place and a pause is always one click away from continuing.
// Restore initial writes *this run's* starting geometry back, so it is the way
// out of a run that went further than wanted (and the way back to Start).
function mmButtonsHTML(){
    if(mmRun&&(mmRun.state==='running'||mmRun.state==='paused')){
        var paused=mmRun.state==='paused';
        return '<div class="mm-btns">'+
            '<button id="mm-toggle">'+(paused?'Resume':'Pause')+'</button>'+
            '<button id="mm-restore">Restore initial</button></div>';
    }
    var n=selectedAtoms.length;
    var html='<button class="mm-run" id="mm-run"'+(n===0?' disabled':'')+'>Start Optimization</button>';
    if(mmRun&&!mmRun.restored){
        html+='<div class="mm-btns"><button id="mm-restore">Restore initial</button></div>';
    }
    return html;
}
function refreshMmPanel(){
    var el=mmPanelEl();
    if(!el||currentMode!=='mmOpt')return;
    var params=MD.mmParams||{};
    var lines=params.lines||[];
    var html='<div class="mm-head"><span>Molecular Mechanics &mdash; '+(params.name||'optimisation')+'</span><span class="mm-close" id="mm-close">×</span></div>';
    for(var i=0;i<lines.length;i++)html+='<div class="mm-param">'+lines[i]+'</div>';
    var ctls=mmControls();
    if(ctls.length){
        html+='<div class="mm-sec">Parameters</div>';
        for(i=0;i<ctls.length;i++)html+=mmControlHTML(ctls[i]);
    }
    html+='<div class="mm-prog" id="mm-prog"></div>';
    // The canvases exist as soon as there is a run, so the curves can grow with
    // it — the panel is not necessarily re-rendered while a run is streaming.
    if(mmRun){
        html+='<div class="mm-curve-title">Energy (kcal/mol)</div><canvas class="mm-curve" id="mm-cv-energy"></canvas>';
        html+='<div class="mm-curve-title">Force (kcal/mol/\u00c5, log)</div><canvas class="mm-curve" id="mm-cv-force"></canvas>';
    }
    var n=selectedAtoms.length;
    if(mmRunActive()){
        // While a run is going, the panel's buttons act on the atoms the *run*
        // named — not on whatever the user has clicked since. Changing the
        // selection mid-run is allowed (it is the next run's), so both sets are
        // named instead of showing only the newer one and implying the run
        // follows it.
        html+='<div class="mm-sel">'+(mmRun.state==='paused'?'Paused on':'Optimizing')+
            ' ('+mmRun.indices.length+'): '+formatAtomList(mmRun.indices)+
            '<br>Selected in the view ('+n+'): '+(n?formatAtomList(selectedAtoms):'none')+'</div>';
    }else{
        html+='<div class="mm-sel">Selected atoms ('+n+'): '+(n?formatAtomList(selectedAtoms):'none &mdash; click atoms in the view, or use Select Atoms / Box Select')+'</div>';
    }
    html+=mmButtonsHTML();
    if(mmErrorMsg)html+='<div class="mm-err">'+mmErrorMsg+'</div>';
    html+=mmResultHTML(mmLastResult);
    el.innerHTML=html;
    wireMmPanel();
    mmUpdateProgressText();
    layoutPanels();
    mmDrawCurves();
}
// The panel's controls are re-rendered on every refresh, so their handlers are
// bound once to the panel itself (the pattern the egg panel uses): a click is
// routed by the element's id, a change by the control's key. Nothing is bound
// to the replaced innerHTML, so no listener can survive its element.
var mmPanelWired=false;
function wireMmPanel(){
    var el=mmPanelEl();
    if(!el||mmPanelWired)return;
    mmPanelWired=true;
    el.addEventListener('click',function(e){
        var t=e&&e.target;
        if(!t||!t.id)return;
        if(t.id==='mm-close')hideMmPanel();
        else if(t.id==='mm-run')startMmOptimization();
        else if(t.id==='mm-toggle')mmToggleRun();
        else if(t.id==='mm-restore')mmRestoreInitial();
    });
    el.addEventListener('change',function(e){
        var t=e&&e.target;
        if(!t||typeof t.getAttribute!=='function')return;
        var key=t.getAttribute('data-mmkey');
        if(!key)return;
        if(key==='preset'){
            // A preset rewrites every other control, so the panel is rebuilt to
            // show the values that will run; "custom" only relabels the choice.
            var v=String(t.value);
            if(v==='custom'){mmOptions.preset='custom';return}
            mmApplyPreset(v);
            refreshMmPanel();
            return;
        }
        var ctls=mmControls();
        for(var k=0;k<ctls.length;k++){
            if(ctls[k].key!==key)continue;
            var opts=ctls[k].options;
            for(var o=0;o<opts.length;o++){
                if(String(opts[o].value)===String(t.value)){mmOptions[key]=opts[o].value;break}
            }
            break;
        }
        // Hand-edited values are a custom choice, and the select says so.
        mmSyncPresetState();
        mmSyncPresetSelect();
    });
}
// Every path that changes the selection calls this, so the panel's list and its
// button state follow the view without the panel having to own the selection.
function mmSelectionChanged(){
    if(currentMode==='mmOpt'&&mmPanelVisible())refreshMmPanel();
}
// Start a run: remember the atoms it may move and the geometry they start from
// (Restore's copy), then ask the host for a session. The answer arrives as a
// stream of mmProgress chunks; nothing is written into the structure here.
function startMmOptimization(){
    if(mmRunActive())return;
    // An Add H request still in flight will append atoms when it lands. Starting
    // a run now would mean the two answers describe different structures, so the
    // start waits for it rather than racing it (the round trip is a millisecond
    // or two — this is a guard against the click, not a long wait).
    if(addHPending){modeInfoEl.textContent='Molecular Mechanics - wait for the Add H request to finish first';return}
    if(diffMode){modeInfoEl.textContent='Molecular Mechanics is not available in diff mode';return}
    if(CRY){modeInfoEl.textContent='Molecular Mechanics is not available for crystal structures';return}
    if(selectedAtoms.length===0){
        modeInfoEl.textContent='Molecular Mechanics - select atoms first: click atoms here, or use Select Atoms / Box Select';
        return;
    }
    // Anything else that owns the coordinates has to let go before the run
    // starts: a playing normal mode and a frame auto-play both write the
    // structure every frame, and a Move Atoms drag is mid-session on it. Each is
    // stopped *before* the starting geometry is captured, so what the run calls
    // "where it began" is the geometry on screen at the click, not a half-drawn
    // vibration amplitude. A manual displacement is closed the same way — the
    // run writes those very coordinates, so its slider must not stay armed over
    // them (the run's own undo snapshot is taken in front of its first write).
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    stopFrameAutoPlay();
    vdCommit();
    var i,indices=selectedAtoms.slice(),origin=new Array(indices.length);
    for(i=0;i<indices.length;i++){
        var a=MD.atoms[indices[i]];
        origin[i]=a?[a.x,a.y,a.z]:null;
    }
    mmRunSeq++;
    mmRun={
        id:'mm-'+mmRunSeq+'-'+Date.now(),
        state:'running',
        atomCount:MD.atoms.length,
        indices:indices,
        origin:{indices:indices.slice(),coords:origin},
        history:[],
        lastStep:-1,
        last:null,
        energyBefore:null,
        undoPushed:false,
        restored:false,
        moved:[],
        movedAll:{},
        movedCount:0
    };
    mmErrorMsg='';
    mmLastResult=null;   // the block described a run this one replaces
    modeInfoEl.textContent='Molecular Mechanics: optimizing '+indices.length+' atom(s)...';
    mmUpdateReopen();
    refreshMmPanel();
    vscodeApi.postMessage({
        command:'mmOptimize',
        runId:mmRun.id,
        atoms:MD.atoms.map(function(a){return{element:a.element,x:a.x,y:a.y,z:a.z}}),
        bonds:MD.bonds.map(function(b){return{atom1:b.atom1,atom2:b.atom2,order:b.order}}),
        selected:indices.slice(),
        options:mmOptionSnapshot()
    });
}
// One progress chunk: the structure follows it live and the curves grow.
function mmHandleProgress(msg){
    if(!mmRunMatches(msg))return;
    // A run that has already stopped on this side (its result applied, its
    // geometry restored, or it was cancelled) takes nothing more from the wire.
    if(mmRun.state==='done')return;
    // While the panel holds the view, a chunk still in flight is *recorded* but
    // not drawn onto the structure. The two halves are separated on purpose: the
    // curve is a record of the search and the host advances its "already sent"
    // mark the moment it posts a point, so a point dropped here is gone for
    // good; the coordinates are not the record, they are the thing the pause is
    // there to hold still. The host's own mmPaused reply then carries the
    // geometry it actually stopped at, so what the user paused on is what they
    // keep looking at — while the curve stays complete.
    var held=(mmRun.state!=='running');
    if(typeof msg.energyBefore==='number')mmRun.energyBefore=msg.energyBefore;
    mmAppendPoints(msg.points);
    mmRun.last={
        step:typeof msg.step==='number'?msg.step:0,
        energy:typeof msg.energy==='number'?msg.energy:0,
        maxGrad:typeof msg.maxGrad==='number'?msg.maxGrad:0,
        rmsGrad:typeof msg.rmsGrad==='number'?msg.rmsGrad:0
    };
    if(!held){
        mmApplyLiveCoords(msg.coords);
        modeInfoEl.textContent='Molecular Mechanics: optimizing '+mmRun.indices.length+
            ' atom(s) - step '+mmRun.last.step+' (E '+mmRun.last.energy.toFixed(2)+' kcal/mol)';
    }
    mmUpdateProgressText();
    mmUpdateReopen();
    mmDrawCurves();
}
// The host stopped pumping at a pause: the session is kept, so Resume continues
// from exactly this geometry and these forces.
function mmHandlePaused(msg){
    if(!mmRunMatches(msg))return;
    if(mmRun.state!=='running'&&mmRun.state!=='paused')return;
    mmAppendPoints(msg.points);
    if(typeof msg.energyBefore==='number')mmRun.energyBefore=msg.energyBefore;
    mmApplyLiveCoords(msg.coords);
    mmRun.state='paused';
    mmRun.last={
        step:typeof msg.step==='number'?msg.step:0,
        energy:typeof msg.energy==='number'?msg.energy:0,
        maxGrad:typeof msg.maxGrad==='number'?msg.maxGrad:0,
        rmsGrad:typeof msg.rmsGrad==='number'?msg.rmsGrad:0
    };
    // The panel is put into its paused state first, then the heavy region
    // rebuild runs behind the guard — so the button the user needs is on screen
    // before any work that could fail.
    refreshMmPanel();
    mmEndLiveRefresh();
    mmDrawCurves();
    modeInfoEl.textContent='Molecular Mechanics: paused at step '+mmRun.last.step+
        ' - Resume continues it, Restore initial puts the starting geometry back';
}
// The host dropped the run (this panel asked for it, or it was replaced): the
// geometry already applied stays, and the panel stops waiting for progress.
function mmHandleCancelled(msg){
    if(!mmRunMatches(msg))return;
    if(!mmRunActive())return;
    mmRun.state='done';
    mmEndLiveRefresh();
    modeInfoEl.textContent='Molecular Mechanics: run cancelled';
    mmUpdateProgressText();
    mmUpdateReopen();
    refreshMmPanel();
}
function applyMmResult(msg){
    if(mmRun&&msg&&msg.runId&&msg.runId!==mmRun.id)return;   // a replaced run's answer
    if(!mmRun){
        // No run on this side any more (cancelled, or the view moved on): the
        // answer addresses a structure state that is no longer the one on
        // screen, so nothing is written — the message is only reported.
        if(msg&&msg.message)modeInfoEl.textContent=msg.message;
        return;
    }
    if(mmRunActive()){mmRun.state='done';mmEndLiveRefresh()}
    if(!msg||!msg.ok){
        mmErrorMsg=(msg&&msg.message)||'Molecular Mechanics failed';
        modeInfoEl.textContent=mmErrorMsg;
        mmUpdateProgressText();
        mmUpdateReopen();
        refreshMmPanel();
        return;
    }
    if(diffMode||CRY){
        // The view moved on (a comparison opened, a crystal was loaded) while
        // the request was in flight; the answer addresses the structure it was
        // computed from, which is no longer the one on screen.
        mmErrorMsg='Molecular Mechanics: the view changed while computing - nothing was applied';
        modeInfoEl.textContent=mmErrorMsg;
        refreshMmPanel();
        return;
    }
    if(typeof msg.atomCount==='number'&&msg.atomCount!==MD.atoms.length){
        mmErrorMsg='Molecular Mechanics: the structure changed while computing - nothing was applied';
        modeInfoEl.textContent=mmErrorMsg;
        refreshMmPanel();
        return;
    }
    if(typeof msg.energyBefore==='number')mmRun.energyBefore=msg.energyBefore;
    mmAppendPoints(msg.points);
    // The final coordinates are authoritative (the live stream has already
    // written them, and this re-writes the same values — the cheap path).
    mmApplyLiveCoords(msg.coords);
    mmRun.state='done';
    mmEndLiveRefresh();
    if(mmRun.movedCount>0)provMark('mmOpt',1);
    mmRun.last={
        step:typeof msg.iterations==='number'?msg.iterations:(mmRun.last?mmRun.last.step:0),
        energy:typeof msg.energyAfter==='number'?msg.energyAfter:(mmRun.last?mmRun.last.energy:0),
        maxGrad:typeof msg.maxGrad==='number'?msg.maxGrad:(mmRun.last?mmRun.last.maxGrad:0),
        rmsGrad:typeof msg.rmsGrad==='number'?msg.rmsGrad:(mmRun.last?mmRun.last.rmsGrad:0)
    };
    mmLastResult={
        energyBefore:typeof msg.energyBefore==='number'?msg.energyBefore:0,
        energyAfter:typeof msg.energyAfter==='number'?msg.energyAfter:0,
        iterations:msg.iterations||0,
        maxGrad:typeof msg.maxGrad==='number'?msg.maxGrad:0,
        rmsGrad:typeof msg.rmsGrad==='number'?msg.rmsGrad:0,
        maxMove:typeof msg.maxMove==='number'?msg.maxMove:0,
        stopReason:msg.stopReason||'',
        counts:msg.counts||null,
        moved:mmRun.movedCount
    };
    mmUpdateProgressText();
    mmUpdateReopen();
    refreshMmPanel();
    modeInfoEl.textContent=msg.message||('Molecular Mechanics: optimized '+mmRun.movedCount+' atom(s)');
}
// One button, two states. The flip is applied here, not waited for: the host's
// answer only confirms it. That matters because the panel must never be left
// showing a state it cannot leave — a lost or slow reply would otherwise strand
// the user on a disabled button, with the run neither paused in the UI nor
// continuable. The host is idempotent for the same reason (a pause that arrives
// twice, or after the run already stopped, is harmless).
function mmToggleRun(){
    if(!mmRun)return;
    if(mmRun.state==='running'){
        mmRun.state='paused';
        vscodeApi.postMessage({command:'mmPause',runId:mmRun.id});
    }else if(mmRun.state==='paused'){
        mmRun.state='running';
        vscodeApi.postMessage({command:'mmResume',runId:mmRun.id});
    }else{
        return;
    }
    mmEndLiveRefresh();
    mmUpdateProgressText();
    refreshMmPanel();
    modeInfoEl.textContent=mmRun.state==='paused'
        ?('Molecular Mechanics: paused at step '+(mmRun.last?mmRun.last.step:0)+
            ' - Resume continues it, Restore initial puts the starting geometry back')
        :('Molecular Mechanics: optimizing '+mmRun.indices.length+' atom(s)...');
}
// Put this run's own starting geometry back. A live run is cancelled first (its
// remaining steps describe a geometry nobody wants any more); the write itself
// is one undo snapshot plus the moved-region refresh.
function mmRestoreInitial(){
    if(!mmRun||mmRun.restored)return;
    var origin=mmRun.origin;
    if(mmRunActive())vscodeApi.postMessage({command:'mmCancel',runId:mmRun.id});
    pushUndo();
    // The run's geometry is being taken back, so its clause comes off the title
    // as well (an Undo of this restore puts the run back, clause included).
    provUnmark('mmOpt');
    var moved=[];
    for(var k=0;k<origin.indices.length;k++){
        var c=origin.coords[k];
        var i=origin.indices[k];
        if(!c||typeof i!=='number'||i<0||i>=MD.atoms.length)continue;
        var a=MD.atoms[i];
        if(a.x!==c[0]||a.y!==c[1]||a.z!==c[2])moved.push(i);
        a.x=c[0];a.y=c[1];a.z=c[2];
    }
    mmRun.state='done';
    mmRun.restored=true;
    mmRun.moved=[];
    updateAtomMeshPositions();
    if(moved.length>0)refreshMovedBondsAndMeshes(moved);
    mmUpdateProgressText();
    mmUpdateReopen();
    refreshMmPanel();
    modeInfoEl.textContent='Molecular Mechanics: restored the geometry this run started from';
}
// Leaving the mode — or any edit that replaces the structure — cancels the run,
// so a late chunk can never write into a structure the user has moved on from.
function mmCancelRun(){
    if(!mmRun)return;
    if(mmRunActive())vscodeApi.postMessage({command:'mmCancel',runId:mmRun.id});
    mmEndLiveRefresh();
    mmRun=null;
    mmUpdateProgressText();
    mmUpdateReopen();
}

var MODE_INFO={view:'View Mode',bondLength:'Bond Length - Click 2 atoms',bondAngle:'Bond Angle - Click 3 atoms (central 2nd)',dihedral:'Dihedral - Click 4 atoms',addAtom:'Add Atom - Click anchor atom',deleteAtom:'Delete Atom - Click atom to delete',replaceAtom:'Replace Atom - Click atom(s) to select, click button again to confirm',addH:'Add H - Click atom(s) to select, click the button again to complete their hydrogens',selectAtoms:'Select Atoms - Input indices or element symbols',boxSelect:'Box Select - Drag a rectangle to select atoms',rotateGroup:'Rotate Group - Select group + axis atoms, then adjust angle',mmOpt:'Molecular Mechanics - Click atoms to select them; the panel optimizes only those atoms',moveAtoms:'Move Atoms - Select atoms; drag/arrows rotate view, right-drag pans; Ctrl+drag or Ctrl+arrows rotates group, Ctrl+right-drag moves group'};

function setMode(m){
    // An Alt-drag that is still running belongs to the mode it started in: a
    // mode switch commits the coordinates it has written so far.
    if(bakeActive)endBakeSession();
    if(m==='moveAtoms'&&CRY){
        // Moving atoms in a periodic supercell breaks image consistency;
        // the mode is only offered for molecular structures.
        modeInfoEl.textContent='Move Atoms is not available for crystal structures';
        return;
    }
    if(m==='mmOpt'&&(CRY||diffMode)){
        // The force field relaxes a molecular structure against its own
        // neighbours; a periodic cell would need image atoms it does not have,
        // and a two-structure comparison has two structures to choose from.
        modeInfoEl.textContent='Molecular Mechanics is not available '+(CRY?'for crystal structures':'in diff mode');
        return;
    }
    if(m==='addH'&&(CRY||diffMode)){
        // Completing hydrogens writes new atoms into the structure on screen,
        // which a periodic cell and a two-structure comparison do not have a
        // place for.
        modeInfoEl.textContent='Add H is not available '+(CRY?'for crystal structures':'in diff mode');
        return;
    }
    if(currentMode===m){
        if(m==='selectAtoms'){
            var sp=document.getElementById('select-panel');
            if(sp)sp.style.display=(sp.style.display==='none'?'block':'none');
        }
        if(m==='rotateGroup'){
            var rp=document.getElementById('rotate-panel');
            if(rp)rp.style.display=(rp.style.display==='none'?'block':'none');
        }
        if(m==='replaceAtom'){
            if(selectedAtoms.length>0){showReplaceAtomModal();return}
        }
        if(m==='addH'){
            if(selectedAtoms.length>0){startAddHydrogens();return}
        }
        if(m==='mmOpt'){
            // The toolbar button always opens the panel again (the × only
            // closes it), so there is never a state without a way back.
            var mp0=mmPanelEl();
            if(mp0&&!mp0.classList.contains('show')){showMmPanel()}else{refreshMmPanel()}
            layoutPanels();
            return;
        }
        layoutPanels();
        return;
    }
    var oldMode=currentMode;
    currentMode=m;
    bakeHintShown=false;
    updateBakeHint();
    var preserveSel=((oldMode==='selectAtoms'||oldMode==='boxSelect'||oldMode==='replaceAtom'||oldMode==='addH'||oldMode==='rotateGroup'||oldMode==='moveAtoms'||oldMode==='mmOpt')&&(m==='selectAtoms'||m==='boxSelect'||m==='replaceAtom'||m==='addH'||m==='rotateGroup'||m==='moveAtoms'||m==='mmOpt'));
    if(oldMode==='moveAtoms'&&m!=='moveAtoms')stopMoveSession();
    if(oldMode==='mmOpt'&&m!=='mmOpt')mmCancelRun();
    // Auto Rotate is a View-mode control: leaving the mode takes its panel away,
    // so the spin it started ends with it rather than running on with no way
    // back to the switch.
    if(oldMode==='view'&&m!=='view')stopSpin();
    if(!preserveSel){
        selectedAtoms=[];diffSelectedAtoms=[];originalCoords=null;
    }
    hideModal();highlightSelected();
    modeInfoEl.textContent=MODE_INFO[m]||m;
    if(!preserveSel)selInfoEl.textContent='';
    document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.classList.toggle('active',b.dataset.mode===m)});
    var sp2=document.getElementById('select-panel');
    if(sp2)sp2.style.display=(m==='selectAtoms')?'block':'none';
    var rp2=document.getElementById('rotate-panel');
    if(rp2)rp2.style.display=(m==='rotateGroup')?'block':'none';
    var vp2=document.getElementById('view-panel');
    if(vp2)vp2.style.display=(m==='view')?'block':'none';
    // The Edit Coordinates card is the second half of the View-mode controls and
    // shares its visibility, but it is its own box in layoutPanels().
    var bp2=document.getElementById('bake-panel');
    if(bp2)bp2.style.display=(m==='view')?'block':'none';
    var mp2=mmPanelEl();
    if(mp2){
        if(m==='mmOpt'){showMmPanel()}else if(mp2.classList.contains('show')){mp2.classList.remove('show')}
    }
    if(m==='rotateGroup'){
        resetRotAxisState();
        syncRotatePanelGroupField();
    }
    if(oldMode==='rotateGroup'&&m!=='rotateGroup'){
        removeRotAxisLine();
        // 分组旋转可能已经改变了质心: 离开该模式时重定旋转中心
        updateRotationCenter();
    }
    layoutPanels();
}

function resetSelection(){
    selectedAtoms=[];diffSelectedAtoms=[];originalCoords=null;highlightSelected();
    selInfoEl.textContent='';
    mmSelectionChanged();
}
/**
 * Select every atom on screen — the keyboard's counterpart to dragging a box
 * over the whole view, which in Diff mode means both halves.
 *
 * Nothing here changes the structure, so it is allowed in every mode and takes
 * no undo step; the only thing it waits for is a closed dialog, whose preview
 * was captured for the selection it opened with. A move/rotate session owns the
 * selection it started with, so one in flight is ended first (exactly as
 * releasing the mouse would) rather than being handed a different atom set
 * halfway through.
 */
function selectAllAtoms(){
    if(modalOpen()){
        modeInfoEl.textContent='Select All: close the open dialog first';
        return false;
    }
    if(!MD.atoms||MD.atoms.length===0){
        modeInfoEl.textContent='Select All: the structure is empty';
        return false;
    }
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    selectedAtoms=[];
    for(var i=0;i<MD.atoms.length;i++)selectedAtoms.push(i);
    if(diffMode){
        diffSelectedAtoms=[];
        var dn=(diffData&&diffData.atoms)?diffData.atoms.length:0;
        for(var j=0;j<dn;j++)diffSelectedAtoms.push(j);
        updateDiffSelInfo();
    }else if(currentMode==='selectAtoms'||currentMode==='replaceAtom'||currentMode==='addH'||currentMode==='mmOpt'){
        // These modes word a finished selection themselves, and their label
        // path opens nothing (a two-atom structure must not pop a dialog here).
        checkSelectionComplete();
    }else{
        selInfoEl.textContent='Selected: '+formatAtomList(selectedAtoms)+' ('+selectedAtoms.length+' atoms)';
    }
    highlightSelected();
    mmSelectionChanged();
    modeInfoEl.textContent='Selected all '+selectedAtoms.length+' atom'+(selectedAtoms.length===1?'':'s');
    return true;
}

document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.addEventListener('click',function(){
    var mode=this.dataset.mode;
    if(mode==='deleteAtom'&&selectedAtoms.length>0&&(currentMode==='selectAtoms'||currentMode==='boxSelect')){
        showBatchDeleteModal();
    }else{
        setMode(mode);
    }
})});
(function(){
    var inp=document.getElementById('sel-panel-input');
    var okBtn=document.getElementById('sel-panel-ok');
    var clrBtn=document.getElementById('sel-panel-clear');
    if(!inp||!okBtn||!clrBtn)return;
    function doSelect(){
        var input=inp.value.trim();
        if(!input)return;
        var indices=parseSelectInput(input);
        indices.forEach(function(idx){
            if(selectedAtoms.indexOf(idx)<0)selectedAtoms.push(idx);
        });
        highlightSelected();
        var names=formatAtomList(selectedAtoms);
        selInfoEl.textContent='Selected: '+names+' ('+selectedAtoms.length+' atoms)';
        inp.value='';
        mmSelectionChanged();
    }
    okBtn.addEventListener('click',doSelect);
    inp.addEventListener('keydown',function(e){if(e.key==='Enter'){e.preventDefault();doSelect()}});
    clrBtn.addEventListener('click',function(){
        selectedAtoms=[];diffSelectedAtoms=[];highlightSelected();selInfoEl.textContent='';
        mmSelectionChanged();
    });
})();

// ===== Rotate Group Around Axis =====
var rotActive=false;
var rotateGroupAtoms=[];
var rotateAxisAtoms=[];
var rotAxisParams=null;
var rotAxisLineMesh=null;
var rotCommittedAngle=0;

function resetRotAxisState(){
    rotActive=false;
    rotateGroupAtoms=[];
    rotateAxisAtoms=[];
    rotAxisParams=null;
    rotCommittedAngle=0;
    removeRotAxisLine();
    var ew=document.getElementById('rp-error');if(ew)ew.textContent='';
    var sw=document.getElementById('rp-slider-wrap');if(sw)sw.style.display='none';
    var slider=document.getElementById('rp-slider');if(slider)slider.value=0;
    var ain=document.getElementById('rp-angle-input');if(ain)ain.value='0';
    var av=document.getElementById('rp-angle-val');if(av)av.textContent='0°';
    // 分组旋转会改变质心, 而这里是每个旋转会话的收尾点 (Done / Clear / 换模式 /
    // 撤销 / 换帧 …), 所以统一在这里重定旋转中心。
    updateRotationCenter();
}
function syncRotatePanelGroupField(){
    if(currentMode!=='rotateGroup')return;
    var inp=document.getElementById('rp-group-input');
    if(!inp)return;
    inp.value=selectedAtoms.map(function(i){return i+1}).join(' ');
}
function setRotatePanelError(msg){
    var ew=document.getElementById('rp-error');
    if(ew)ew.textContent=msg||'';
}
function showRotAxisLine(){
    removeRotAxisLine();
    if(!rotAxisParams)return;
    var p=rotAxisParams;
    var len=maxD*2+5;
    var s=new THREE.Vector3(p.ox-CX-p.dx*len,p.oy-CY-p.dy*len,p.oz-CZ-p.dz*len);
    var e=new THREE.Vector3(p.ox-CX+p.dx*len,p.oy-CY+p.dy*len,p.oz-CZ+p.dz*len);
    var dir=new THREE.Vector3().subVectors(e,s);
    var hl=dir.length();
    if(hl<0.001)return;
    var g=new THREE.CylinderGeometry(0.04,0.04,hl,6,1);
    var m=new THREE.MeshPhongMaterial({color:0xff8800,transparent:true,opacity:0.8});
    var mesh=new THREE.Mesh(g,m);
    var mid=new THREE.Vector3().addVectors(s,e).multiplyScalar(0.5);
    mesh.position.copy(mid);
    var axis=new THREE.Vector3(0,1,0);
    mesh.quaternion.copy(new THREE.Quaternion().setFromUnitVectors(axis,dir.normalize()));
    moleculeGroup.add(mesh);
    rotAxisLineMesh=mesh;
    needsRender=true;
}
function removeRotAxisLine(){
    if(rotAxisLineMesh){
        moleculeGroup.remove(rotAxisLineMesh);
        disposeMesh(rotAxisLineMesh);
        rotAxisLineMesh=null;
        needsRender=true;
    }
}
function validateRotationAxis(axisIdx){
    if(!axisIdx||axisIdx.length<2){
        return {valid:false,error:'Rotation axis needs at least 2 atoms.'};
    }
    var seen={};
    for(var k=0;k<axisIdx.length;k++){
        if(seen[axisIdx[k]]){
            return {valid:false,error:'Duplicate atom in axis: '+(axisIdx[k]+1)+'.'};
        }
        seen[axisIdx[k]]=true;
    }
    var axisSet={};
    axisIdx.forEach(function(i){axisSet[i]=true});
    for(var bi=0;bi<MD.bonds.length;bi++){
        var b=MD.bonds[bi];
        if(axisSet[b.atom1]&&axisSet[b.atom2]){
            if(b.order>=1.75){
                var ordStr=b.order>=2.5?'triple':'double';
                return {valid:false,error:'Atoms '+(b.atom1+1)+' and '+(b.atom2+1)+' on the axis are connected by a '+ordStr+' bond. A rotation axis must use single bonds.'};
            }
        }
    }
    return {valid:true};
}
function computeRotAxisParams(axisIdx){
    if(axisIdx.length<2)return null;
    var p0=MD.atoms[axisIdx[0]],p1=MD.atoms[axisIdx[axisIdx.length-1]];
    var dx=p1.x-p0.x,dy=p1.y-p0.y,dz=p1.z-p0.z;
    var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(l<1e-6)return null;
    dx/=l;dy/=l;dz/=l;
    return {ox:p0.x,oy:p0.y,oz:p0.z,dx:dx,dy:dy,dz:dz};
}
function applyGroupRotation(angleDeg){
    if(!rotAxisParams||!originalCoords)return;
    var p=rotAxisParams;
    var rad=angleDeg*Math.PI/180;
    var axisSet={};
    rotateAxisAtoms.forEach(function(i){axisSet[i]=true});
    for(var gi=0;gi<rotateGroupAtoms.length;gi++){
        var idx=rotateGroupAtoms[gi];
        if(axisSet[idx])continue;
        var oc=originalCoords[idx];
        if(!oc)continue;
        var r=rotAroundAxis(oc.x,oc.y,oc.z,p.ox,p.oy,p.oz,p.dx,p.dy,p.dz,rad);
        MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
    }
    updateScenePositions(true);
}
function commitRotAngle(newAngle){
    // Commit one angle adjustment: push an undo snapshot of the state at the
    // previously committed angle (tagged with rotAngle so doUndo can sync the
    // slider display), then apply the new angle.
    if(!rotActive||!rotAxisParams||!originalCoords)return;
    if(isNaN(newAngle))newAngle=rotCommittedAngle;
    newAngle=Math.max(0,Math.min(360,newAngle));
    if(Math.abs(newAngle-rotCommittedAngle)>1e-9){
        restoreOriginal();
        applyGroupRotation(rotCommittedAngle);
        pushUndo();
        provMark('rotate',1);
        undoStack[undoStack.length-1].rotAngle=rotCommittedAngle;
        restoreOriginal();
        applyGroupRotation(newAngle);
        rotCommittedAngle=newAngle;
        // 一次角度调整落定后质心可能变了, 立刻重定旋转中心 (与滑块拖动中的
        // 实时预览不同: 预览期间中心保持不变, 免得整屏跟着跳)。
        updateRotationCenter();
    }
    var slider=document.getElementById('rp-slider');
    if(slider)slider.value=newAngle;
    var av=document.getElementById('rp-angle-val');
    if(av)av.textContent=newAngle.toFixed(1)+'°';
    var ain=document.getElementById('rp-angle-input');
    if(ain)ain.value=newAngle;
    // 表面片在预览期间被隐藏, 提交后立刻恢复 (签名没变时是空操作)。
    refreshVmdSurfaces();
}
(function(){
    var gi=document.getElementById('rp-group-input');
    var ai=document.getElementById('rp-axis-input');
    var applyBtn=document.getElementById('rp-apply');
    var clearBtn=document.getElementById('rp-clear');
    var slider=document.getElementById('rp-slider');
    var doneBtn=document.getElementById('rp-done');
    var resetBtn=document.getElementById('rp-reset');
    if(!gi||!ai||!applyBtn||!clearBtn||!slider||!doneBtn||!resetBtn)return;

    function syncGroupFromField(){
        var parsed=parseSelectInput(gi.value);
        selectedAtoms=parsed;
        highlightSelected();
        var names=formatAtomList(selectedAtoms);
        selInfoEl.textContent=selectedAtoms.length?('Selected: '+names):'';
    }
    gi.addEventListener('keydown',function(e){
        if(e.key==='Enter'){e.preventDefault();syncGroupFromField();}
    });
    ai.addEventListener('keydown',function(e){
        if(e.key==='Enter'){e.preventDefault();applyBtn.click();}
    });
    applyBtn.addEventListener('click',function(){
        var groupField=gi.value.trim();
        var axisField=ai.value.trim();
        var groupIdx=groupField?parseSelectInput(groupField):selectedAtoms.slice();
        var axisIdx=parseSelectInput(axisField);
        var verr=validateRotationAxis(axisIdx);
        if(!verr.valid){setRotatePanelError(verr.error);return;}
        if(groupIdx.length===0){setRotatePanelError('Please specify at least one group atom.');return;}
        var params=computeRotAxisParams(axisIdx);
        if(!params){setRotatePanelError('Axis atoms are too close together to define a direction.');return;}
        setRotatePanelError('');
        rotateGroupAtoms=groupIdx.slice();
        rotateAxisAtoms=axisIdx.slice();
        rotAxisParams=params;
        // Capture the session baseline WITHOUT pushing an undo snapshot here —
        // each committed angle adjustment pushes its own snapshot, so undo
        // steps back one angle adjustment at a time.
        originalCoords=MD.atoms.map(function(a){return{x:a.x,y:a.y,z:a.z}});
        rotCommittedAngle=0;
        var sw=document.getElementById('rp-slider-wrap');if(sw)sw.style.display='block';
        slider.value=0;
        var ain=document.getElementById('rp-angle-input');if(ain)ain.value='0';
        var av=document.getElementById('rp-angle-val');if(av)av.textContent='0°';
        showRotAxisLine();
        rotActive=true;
        applyGroupRotation(0);
        modeInfoEl.textContent='Rotate Group - adjust angle (0-360°)';
    });
    slider.addEventListener('input',function(){
        if(!rotActive||!originalCoords)return;
        var v=parseFloat(slider.value);
        restoreOriginal();
        applyGroupRotation(v);
        var av=document.getElementById('rp-angle-val');
        if(av)av.textContent=v.toFixed(1)+'°';
        var ain=document.getElementById('rp-angle-input');
        if(ain)ain.value=v;
    });
    slider.addEventListener('change',function(){
        if(!rotActive)return;
        commitRotAngle(parseFloat(slider.value));
    });
    var angleInput=document.getElementById('rp-angle-input');
    if(angleInput){
        angleInput.addEventListener('keydown',function(e){
            if(e.key!=='Enter')return;
            e.preventDefault();
            commitRotAngle(parseFloat(angleInput.value));
        });
        angleInput.addEventListener('change',function(){
            if(!rotActive)return;
            commitRotAngle(parseFloat(angleInput.value));
        });
    }
    doneBtn.addEventListener('click',function(){
        restoreOriginal();
        applyGroupRotation(parseFloat(slider.value));
        resetRotAxisState();
        selectedAtoms=[];originalCoords=null;highlightSelected();
        selInfoEl.textContent='';
        gi.value='';ai.value='';
        modeInfoEl.textContent=MODE_INFO.rotateGroup;
    });
    resetBtn.addEventListener('click',function(){
        if(originalCoords)restoreOriginal();
        rebuildScene();
        // Remove all undo snapshots pushed by this rotation session (tagged
        // with rotAngle) so Reset reverts the whole session at once.
        while(undoStack.length>0&&undoStack[undoStack.length-1].rotAngle!=null){
            undoStack.pop();
        }
        updateUndoBtn();
        resetRotAxisState();
    });
    clearBtn.addEventListener('click',function(){
        resetRotAxisState();
        selectedAtoms=[];originalCoords=null;highlightSelected();
        selInfoEl.textContent='';
        gi.value='';ai.value='';
    });
})();
document.getElementById('reset-btn').addEventListener('click',function(){
    stopSpin();
    rotQuat.identity();panX=0;panY=0;
    if(diffMode){
        camDist=camDiffInit();diffCamDist=camDist;
        diffRotQuat.identity();diffPanX=0;diffPanY=0;
    }else{
        camDist=initCam;camera.position.set(0,0,camDist);
    }
    updateTransform();
});
var undoStack=[];var MAX_UNDO=50;
// What a saved file's title line says about where the structure came from and
// what was done to it (the functions that build and record it are at "Save
// provenance", next to doSave). Declared here, with the other view state, so it
// is never undefined when an edit path reads it.
var PROV_PREFIX='MV: ';   // our own marker, so re-saving a file never stacks lines
var PROV_MAX=190;         // characters — a title line, not a paragraph
var provSource=provBaseName(MD.filePath);
var provMerged=[];        // files merged into this structure, in order
var provEdits=[];         // [{k,n}] — what changed it, in the order it happened
var provModes={};         // mode index -> {freq,scale}: displacement kept in the geometry

function pushUndo(){
    // Every coordinate edit opens with a snapshot. An open manual displacement
    // is a *pending* rewrite of the whole structure — its next slider move would
    // put every atom back onto base + amplitude — so it is closed here, once,
    // for every edit path at the one door they all pass through. The geometry it
    // produced stays exactly as it is, and the snapshot below therefore captures
    // the displaced structure the edit is about to change.
    if(vdActive&&!vdSuppress)vdCommit();
    // hasExplicitBonds is part of the restore state: an Import can flip it
    // false→true, and an undo of that import must restore it, otherwise Move
    // Atoms would keep using the explicit-bonds (deletion-only) path on what
    // is an auto-detected structure again.
    var snap={atoms:JSON.parse(JSON.stringify(MD.atoms)),bonds:JSON.parse(JSON.stringify(MD.bonds)),hasExplicitBonds:!!MD.hasExplicitBonds,prov:provClone()};
    if(CRY){snap.baseAtoms=JSON.parse(JSON.stringify(CRY.baseAtoms));snap.baseBonds=JSON.parse(JSON.stringify(CRY.baseBonds))}
    undoStack.push(snap);
    if(undoStack.length>MAX_UNDO)undoStack.shift();
    updateUndoBtn();
}
function updateUndoBtn(){var b=document.getElementById('undo-btn');if(b)b.disabled=undoStack.length===0}
function doUndo(){
    if(undoStack.length===0)return;
    if(mmRun)mmCancelRun();
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||moveKeyActive||fragRotActive)stopMoveSession();
    if(bakeActive)endBakeSession();
    // The snapshot replaces every coordinate, so a displacement session has
    // nothing left to be measured from — it is dropped, not committed.
    vdEnd();
    var snap=undoStack.pop();
    MD.atoms=snap.atoms;MD.bonds=snap.bonds;
    MD.hasExplicitBonds=!!snap.hasExplicitBonds;
    // The provenance belongs to the geometry: an undone edit is not in the
    // structure any more, so its clause comes off the title with it.
    provRestore(snap.prov);
    if(CRY&&snap.baseAtoms){CRY.baseAtoms=snap.baseAtoms;CRY.baseBonds=snap.baseBonds}
    if(CRY){rebuildCrystal()}rebuildScene();updateUndoBtn();
    if(diffMode)resetSelection();
    if(currentMode==='mmOpt'){
        // The panel's last-run block describes a geometry that has just been
        // undone, so it goes with it.
        mmLastResult=null;mmErrorMsg='';
        refreshMmPanel();
    }
    if(currentMode==='rotateGroup'&&rotActive){
        if(snap.rotAngle!=null){
            // Undo one angle adjustment: sync slider and displays to the
            // restored angle so the panel stays consistent with the molecule.
            rotCommittedAngle=snap.rotAngle;
            var slider=document.getElementById('rp-slider');
            if(slider)slider.value=snap.rotAngle;
            var av=document.getElementById('rp-angle-val');
            if(av)av.textContent=snap.rotAngle.toFixed(1)+'°';
            var ain=document.getElementById('rp-angle-input');
            if(ain)ain.value=snap.rotAngle;
            showRotAxisLine();
        }else{
            // Undo stepped past this rotation session's own history (e.g. an
            // earlier unrelated edit) — close the session cleanly.
            resetRotAxisState();
            originalCoords=null;
        }
    }
}
document.getElementById('undo-btn').addEventListener('click',doUndo);
document.getElementById('save-btn').addEventListener('click',doSave);
document.getElementById('import-btn').addEventListener('click',function(){
    if(diffMode){modeInfoEl.textContent='Import is not available in diff mode';return}
    if(CRY){modeInfoEl.textContent='Import is not available for crystal structures';return}
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    modeInfoEl.textContent='Import: choose source...';
    vscodeApi.postMessage({command:'importFile'});
});
document.getElementById('diff-btn').addEventListener('click',function(){
    if(diffMode){
        exitDiff();
    }else{
        modeInfoEl.textContent='Diff: selecting file...';
        vscodeApi.postMessage({command:'diffFile'});
    }
});

var diffMode=false;
var diffData=null;
var diffPivot=null;
var diffMolGroup=null;
var diffAtomMeshes=[];
// The comparison's bonds use the same instanced cylinders as the main
// structure, in two pools: the plain one, and one whose material carries the
// orange glow a highlighted bond is drawn with. A highlighted bond simply
// belongs to the other pool, so the highlight costs a refill of the two lists.
var diffBondPieces=[];
var diffBondPoolPlain=newBondPool(40,0,0,null);
var diffBondPoolHot=newBondPool(40,0xff6600,0.3,null);
var _dbHotCol=new THREE.Color();
function addDiffBondPiece(px,py,pz,dx,dy,dz,l,r,col,key,hi){
    diffBondPieces.push({x:px,y:py,z:pz,dx:dx,dy:dy,dz:dz,l:l,r:r,col:col,key:key,hi:!!hi});
}
// Refill both pools from the piece list: each piece goes to the pool its
// highlight state selects, and its colour follows. Called when the comparison
// is built and whenever the highlighted set changes.
function diffBuildBondPools(){
    _dbHotCol.setRGB(1,0.4,0);
    diffBondPoolPlain.parent=diffMolGroup;
    diffBondPoolHot.parent=diffMolGroup;
    diffBondPoolPlain.pieces.length=0;
    diffBondPoolHot.pieces.length=0;
    for(var i=0;i<diffBondPieces.length;i++){
        var p=diffBondPieces[i];
        (p.hi?diffBondPoolHot:diffBondPoolPlain).pieces.push(p);
    }
    diffWriteBondPool(diffBondPoolPlain);
    diffWriteBondPool(diffBondPoolHot);
    needsRender=true;
}
function diffWriteBondPool(pool){
    if(pool.pieces.length===0){poolClear(pool);return}
    poolEnsure(pool,pool.pieces.length);
    for(var i=0;i<pool.pieces.length;i++){
        var p=pool.pieces[i];
        poolPut(pool,i,p.x,p.y,p.z,p.dx,p.dy,p.dz,p.l,p.r,p.hi?_dbHotCol:p.col);
    }
    poolUpload(pool);
}
var diffCX=0,diffCY=0,diffCZ=0;
var diffRotQuat=new THREE.Quaternion();
var diffThresholdPct=2;
var diffDetBonds1=null;
var diffDetBonds2=null;
var diffPanX=0,diffPanY=0;
var diffCamDist=10;
var diffTransformSide='left';
var diffActiveSide='left';
var diffSelectedAtoms=[];
var diffMapping=null;
var diffReverseMapping=null;
var diffPanelEl=document.getElementById('diff-panel');
var diffReopenEl=document.getElementById('diff-reopen');
var diffPanelHTML='';
var diffResultsHTML='';
var diffStatsExpanded=true;
var diffLabelLeft=document.getElementById('diff-label');
var diffLabelRight=document.getElementById('diff-label-right');

window.addEventListener('message',function(event){
    var msg=event.data;
    if(msg.command==='diffResult'){
        if(msg.cancelled){
            modeInfoEl.textContent='View Mode';
            return;
        }
        startDiff(msg);
    }else if(msg.command==='importResult'){
        if(msg.cancelled){
            modeInfoEl.textContent=MODE_INFO[currentMode]||currentMode;
            return;
        }
        importStructure(msg);
    }else if(msg.command==='hydrogenResult'){
        applyHydrogenResult(msg);
    }else if(msg.command==='mmProgress'){
        mmHandleProgress(msg);
    }else if(msg.command==='mmPaused'){
        mmHandlePaused(msg);
    }else if(msg.command==='mmCancelled'){
        mmHandleCancelled(msg);
    }else if(msg.command==='mmResult'){
        applyMmResult(msg);
    }
});

/**
 * Where a block of incoming atoms goes so that it cannot overlap anything on
 * screen: its centroid is put one bounding-sphere radius plus a buffer away from
 * the existing structure's, along +X.
 *
 * A bounding sphere contains every atom, so this single O(n) test guarantees no
 * incoming atom overlaps ANY existing atom — including blocks added earlier,
 * which are part of "existing" by now. Rotation-invariant, no overlap iterations
 * needed. The two user-visible placements (Import, and paste of a copied
 * fragment) share it, so a pasted copy lands where an import would.
 */
function placeBlockBeside(newAtoms){
    var n=newAtoms.length,i,a;
    var c2x=0,c2y=0,c2z=0;
    for(i=0;i<n;i++){a=newAtoms[i];c2x+=a.x;c2y+=a.y;c2z+=a.z}
    if(n>0){c2x/=n;c2y/=n;c2z/=n}
    var N0=MD.atoms.length;
    if(N0===0)return{x:-c2x,y:-c2y,z:-c2z};
    var c1x=0,c1y=0,c1z=0,R1=0,R2=0;
    for(i=0;i<N0;i++){a=MD.atoms[i];c1x+=a.x;c1y+=a.y;c1z+=a.z}
    c1x/=N0;c1y/=N0;c1z/=N0;
    for(i=0;i<N0;i++){
        a=MD.atoms[i];
        var d1=Math.sqrt((a.x-c1x)*(a.x-c1x)+(a.y-c1y)*(a.y-c1y)+(a.z-c1z)*(a.z-c1z));
        if(d1>R1)R1=d1;
    }
    for(i=0;i<n;i++){
        a=newAtoms[i];
        var d2=Math.sqrt((a.x-c2x)*(a.x-c2x)+(a.y-c2y)*(a.y-c2y)+(a.z-c2z)*(a.z-c2z));
        if(d2>R2)R2=d2;
    }
    // buf=2.5 Å per side covers any atom radius, so the closest two atoms of the
    // two blocks end up 2*2.5+2.0 = 7.0 Å apart (the import has always spaced
    // them exactly this way) — visibly separate, and never bonded by accident.
    var BUF=2.5,GAP=2.0;
    var D=(R1+BUF)+(R2+BUF)+GAP;
    return{x:c1x+D-c2x,y:c1y-c2y,z:c1z-c2z};
}

function importStructure(msg){
    // Append an imported structure (parsed by the extension host) beside the
    // current structure at a guaranteed non-overlapping position. One undo
    // snapshot per import.
    if(diffMode||CRY)return;
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    mmCancelRun();
    var newAtoms=msg.atoms||[];
    var newBonds=msg.bonds||[];
    if(newAtoms.length===0)return;
    pushUndo();
    var N0=MD.atoms.length;
    // --- Non-overlap placement ---
    // The geometry lives in placeBlockBeside, which the fragment paste uses too:
    // one rule, so a pasted copy lands exactly where an imported structure would.
    var t=placeBlockBeside(newAtoms);
    for(var i=0;i<newAtoms.length;i++){
        var a=newAtoms[i];
        MD.atoms.push({element:a.element,x:a.x+t.x,y:a.y+t.y,z:a.z+t.z,color:a.color||((MD.atomColors&&MD.atomColors[a.element])||'#FF1493')});
    }
    // --- Bonds ---
    if(!MD.hasExplicitBonds&&!msg.hasExplicitBonds){
        // Fully automatic bond detection for the combined system (regular
        // behavior for auto-detected files).
        MD.bonds=detectBondsFromAtoms(MD.atoms,MD.charge);
    }else{
        // At least one side has authored bonds: keep them and connect the
        // imported part explicitly instead of re-detecting everything
        // (re-detection could alter authored bond orders).
        MD.hasExplicitBonds=true;
        if(msg.hasExplicitBonds){
            newBonds.forEach(function(b){
                if(b.atom1<newAtoms.length&&b.atom2<newAtoms.length)
                    MD.bonds.push({atom1:N0+b.atom1,atom2:N0+b.atom2,order:b.order==null?1:b.order});
            });
        }else{
            // Imported file had no authored bonds: auto-detect, keeping only
            // bonds that touch a new atom (existing bonds stay as authored).
            var det=detectBondsFromAtoms(MD.atoms,MD.charge);
            det.forEach(function(b){if(b.atom1>=N0||b.atom2>=N0)MD.bonds.push(b)});
        }
    }
    rebuildScene();
    // Zoom out if the combined structure no longer fits the current view.
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    var need=maxD*2.5+5;
    if(camDist<need){camDist=need;camera.position.set(0,0,camDist);camera.lookAt(0,0,0)}
    needsRender=true;
    provAddSource(msg.fileName);
    modeInfoEl.textContent='Imported '+newAtoms.length+' atoms from '+(msg.fileName||'file')+' - placed beside the existing structure';
}

// ===== Fragment clipboard (Ctrl+C / Ctrl+X / Ctrl+V) =====
// The clipboard holds one fragment, in coordinates relative to its own centroid,
// so it carries no position: pasting is a translation to a free spot, and the
// same fragment can be pasted again and again. The bonds *inside* the fragment
// travel with it verbatim — a copy must be a copy, bond orders included — and no
// bond is invented between the copy and whatever it lands beside, which is why a
// paste does not re-detect the combined system the way an import does.
//
// It is the viewer's own clipboard: the system clipboard is left alone, because
// reading it would need a permission the webview may refuse and writing it would
// replace text the user copied elsewhere. A copy is also mirrored into
// localStorage, best effort and re-validated on the way back, so a fragment can
// be pasted into another panel; if storage is unavailable the fragment simply
// stays in memory for this panel.
var CLIP_KEY='molecularViewer3d.fragment';
var PASTE_ATOM_LIMIT=20000;   // the ceiling a single imported block is held to
var clipFrag=null;

function clipStore(){
    if(!clipFrag)return;
    try{window.localStorage.setItem(CLIP_KEY,JSON.stringify(clipFrag))}catch(e){}
}
function clipLoad(){
    clipFrag=null;
    try{
        var raw=window.localStorage.getItem(CLIP_KEY);
        if(raw)clipFrag=clipSanitize(JSON.parse(raw));
    }catch(e){clipFrag=null}
    return clipFrag;
}
/** Whatever comes back from storage is data, not a fragment: it may have been
 *  written by an older version, by another tab, or be something else entirely.
 *  Anything that is not a complete, finite, in-range fragment is discarded
 *  rather than half-used. */
function clipSanitize(c){
    if(!c||typeof c!=='object')return null;
    var src=c.atoms;
    if(!Array.isArray(src)||src.length===0||src.length>PASTE_ATOM_LIMIT)return null;
    var atoms=[],i,a;
    for(i=0;i<src.length;i++){
        a=src[i];
        if(!a||typeof a.element!=='string'||!/^[A-Z][a-z]?$/.test(a.element))return null;
        if(typeof a.x!=='number'||typeof a.y!=='number'||typeof a.z!=='number')return null;
        if(!isFinite(a.x)||!isFinite(a.y)||!isFinite(a.z))return null;
        atoms.push({element:a.element,x:a.x,y:a.y,z:a.z});
    }
    var bonds=[],srcB=Array.isArray(c.bonds)?c.bonds:[];
    for(i=0;i<srcB.length;i++){
        var b=srcB[i];
        if(!b||typeof b.a!=='number'||typeof b.b!=='number')continue;
        if(b.a<0||b.b<0||b.a>=atoms.length||b.b>=atoms.length||b.a===b.b)continue;
        bonds.push({a:b.a,b:b.b,order:(typeof b.order==='number'&&isFinite(b.order))?b.order:1});
    }
    return {v:1,atoms:atoms,bonds:bonds};
}
/** Copy the selected atoms and the bonds among them. False means there was
 *  nothing to copy, so the caller can leave the key to the browser. */
function copySelection(){
    if(selectedAtoms.length===0){
        modeInfoEl.textContent='Copy: select atoms first - click atoms, or use Select Atoms / Box Select';
        return false;
    }
    var idx=selectedAtoms.slice().sort(function(a,b){return a-b});
    var map={},atoms=[],sx=0,sy=0,sz=0,i,a;
    for(i=0;i<idx.length;i++){
        a=MD.atoms[idx[i]];
        if(!a)continue;
        // A selection should not repeat an index, but a repeated one would make
        // the fragment hold the same atom twice and bond the copy to the wrong
        // end, so it is dropped here rather than trusted.
        if(map[idx[i]]!==undefined)continue;
        map[idx[i]]=atoms.length;
        atoms.push({element:a.element,x:a.x,y:a.y,z:a.z});
        sx+=a.x;sy+=a.y;sz+=a.z;
    }
    if(atoms.length===0)return false;
    var cx=sx/atoms.length,cy=sy/atoms.length,cz=sz/atoms.length;
    atoms.forEach(function(q){q.x-=cx;q.y-=cy;q.z-=cz});
    var bonds=[];
    MD.bonds.forEach(function(b){
        var p=map[b.atom1],q=map[b.atom2];
        if(p===undefined||q===undefined)return;
        bonds.push({a:p,b:q,order:b.order});
    });
    clipFrag={v:1,atoms:atoms,bonds:bonds};
    clipStore();
    modeInfoEl.textContent='Copied '+atoms.length+' atom(s)'+(bonds.length?' and '+bonds.length+' bond(s)':'');
    return true;
}
/** Cut = copy, then the same deletion the toolbar's batch delete uses (so a
 *  crystal loses its base atoms and its base bonds exactly as a plain Delete
 *  does), with one undo step for the pair. */
function cutSelection(){
    if(diffMode){modeInfoEl.textContent='Cut is not available in diff mode';return false}
    var n=selectedAtoms.length;
    if(n===0){
        modeInfoEl.textContent='Cut: select atoms first - click atoms, or use Select Atoms / Box Select';
        return false;
    }
    if(!copySelection())return false;
    pushUndo();
    provMark('cutAtoms',n);
    batchDeleteSelected();
    modeInfoEl.textContent='Cut '+n+' atom(s) - Ctrl+V pastes them back';
    return true;
}
/** Paste the clipboard beside the current structure, as one undo step, and
 *  select what arrived. Refused where appending atoms cannot be kept: a
 *  periodic cell (its contents are regenerated from the base atoms, so the copy
 *  would vanish on the next rebuild) and a two-structure comparison (which has
 *  no single structure to grow). */
function pasteClipboard(){
    if(diffMode){modeInfoEl.textContent='Paste is not available in diff mode';return false}
    if(CRY){modeInfoEl.textContent='Paste is not available for crystal structures';return false}
    if(!clipFrag)clipLoad();
    if(!clipFrag){modeInfoEl.textContent='Paste: nothing copied yet';return false}
    var frag=clipFrag;
    if(MD.atoms.length+frag.atoms.length>PASTE_ATOM_LIMIT){
        modeInfoEl.textContent='Paste: the structure would exceed '+PASTE_ATOM_LIMIT+' atoms - nothing was pasted';
        return false;
    }
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    stopFrameAutoPlay();
    pushUndo();
    provMark('pasteAtoms',frag.atoms.length);
    var t=placeBlockBeside(frag.atoms);
    var N0=MD.atoms.length;
    frag.atoms.forEach(function(a){
        MD.atoms.push({element:a.element,x:a.x+t.x,y:a.y+t.y,z:a.z+t.z,
            color:(MD.atomColors&&MD.atomColors[a.element])||'#FF1493'});
    });
    frag.bonds.forEach(function(b){
        MD.bonds.push({atom1:N0+b.a,atom2:N0+b.b,order:b.order});
    });
    rebuildScene();
    // Select what just arrived, so it can be moved, cut or copied straight away.
    selectedAtoms=[];
    for(var i=N0;i<MD.atoms.length;i++)selectedAtoms.push(i);
    highlightSelected();
    selInfoEl.textContent='Selected: '+formatAtomList(selectedAtoms)+' ('+selectedAtoms.length+' atoms)';
    // Zoom out when the copy lands outside the current view (the import's rule).
    var maxD=0;
    MD.atoms.forEach(function(a){
        var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);
        if(dd>maxD)maxD=dd;
    });
    var need=maxD*2.5+5;
    if(camDist<need){camDist=need;camera.position.set(0,0,camDist);camera.lookAt(0,0,0)}
    needsRender=true;
    modeInfoEl.textContent='Pasted '+frag.atoms.length+' atom(s) beside the existing structure';
    return true;
}

// ===== Add H (hydrogen completion) =====
// The geometry is computed in the extension host, by the SMILES -> 3D engine's
// own experience tables (X-H lengths, theoretical angles) with every atom that
// was already in the structure frozen — a hydrogen is placed, an existing bond
// is never moved. The request carries the structure as it stands in the view
// (edits included) and the answer carries only the new hydrogens.
var addHPending=false;
function startAddHydrogens(){
    if(addHPending)return;
    if(diffMode){modeInfoEl.textContent='Add H is not available in diff mode';return}
    if(CRY){modeInfoEl.textContent='Add H is not available for crystal structures';return}
    if(selectedAtoms.length===0){
        modeInfoEl.textContent='Add H - select atoms first: click atoms here, or use Select Atoms / Box Select';
        return;
    }
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    addHPending=true;
    modeInfoEl.textContent='Add H: completing '+selectedAtoms.length+' selected atom(s)...';
    vscodeApi.postMessage({
        command:'addHydrogens',
        atoms:MD.atoms.map(function(a){return{element:a.element,x:a.x,y:a.y,z:a.z}}),
        bonds:MD.bonds.map(function(b){return{atom1:b.atom1,atom2:b.atom2,order:b.order}}),
        selected:selectedAtoms.slice()
    });
}
function applyHydrogenResult(msg){
    addHPending=false;
    if(!msg||!msg.ok){modeInfoEl.textContent=(msg&&msg.message)||'Add H failed';return}
    var newAtoms=Array.isArray(msg.atoms)?msg.atoms:[];
    var newBonds=Array.isArray(msg.bonds)?msg.bonds:[];
    if(newAtoms.length===0){modeInfoEl.textContent=msg.message||'Add H: nothing to add';return}
    if(typeof msg.atomCount==='number'&&msg.atomCount!==MD.atoms.length){
        // The answer names atom indices of the structure it was computed from;
        // if the structure changed while it was in flight, the bonds address
        // atoms that are no longer there.
        modeInfoEl.textContent='Add H: the structure changed while computing - nothing was added';
        return;
    }
    // Appending atoms is the one structure change that does not pass through a
    // mode switch, and the answer lands whenever the host gets to it — by which
    // time the user may have started something else entirely. So every writer is
    // ended here, explicitly, before the structure grows: an optimisation must
    // never keep writing coordinates into a structure it did not start from, a
    // normal-mode playback must not rewrite atoms around the new hydrogens (its
    // own stop would write the pre-add geometry back over them), and a Move
    // Atoms session must not re-derive bonds on a structure that has just gained
    // atoms.
    if(vibActive||vibPaused)stopVibration();
    stopMoveSession();
    mmCancelRun();
    pushUndo();
    provMark('addH',newAtoms.length);
    var N0=MD.atoms.length;
    newAtoms.forEach(function(a){
        MD.atoms.push({element:'H',x:a.x,y:a.y,z:a.z,color:(MD.atomColors&&MD.atomColors.H)||'#FFFFFF'});
    });
    // The bonds arrive with their indices already resolved against the extended
    // atom list; an entry pointing anywhere but at a new hydrogen is dropped
    // rather than stored.
    newBonds.forEach(function(b){
        if(typeof b.atom1!=='number'||typeof b.atom2!=='number')return;
        if(b.atom1<0||b.atom1>=N0)return;
        if(b.atom2<N0||b.atom2>=MD.atoms.length)return;
        MD.bonds.push({atom1:b.atom1,atom2:b.atom2,order:b.order==null?1:b.order});
    });
    rebuildScene();
    modeInfoEl.textContent=msg.message||('Add H: added '+newAtoms.length+' hydrogens');
}

function buildAdjacency(atoms,bonds){
    var adj=[];
    for(var i=0;i<atoms.length;i++)adj[i]=[];
    bonds.forEach(function(b){
        if(adj[b.atom1])adj[b.atom1].push({to:b.atom2,order:b.order});
        if(adj[b.atom2])adj[b.atom2].push({to:b.atom1,order:b.order});
    });
    return adj;
}

function atomSignature(i,atoms,adj){
    var el=atoms[i].element;
    var deg=adj[i].length;
    var nbEls=adj[i].map(function(e){return atoms[e.to].element}).sort().join(',');
    return el+'|'+deg+'|'+nbEls;
}

function findAtomMapping(atoms1,bonds1,atoms2,bonds2){
    var n1=atoms1.length,n2=atoms2.length;
    if(n1!==n2)return null;

    var adj1=buildAdjacency(atoms1,bonds1);
    var adj2=buildAdjacency(atoms2,bonds2);

    var elCount1={},elCount2={};
    atoms1.forEach(function(a){elCount1[a.element]=(elCount1[a.element]||0)+1});
    atoms2.forEach(function(a){elCount2[a.element]=(elCount2[a.element]||0)+1});
    var els1=Object.keys(elCount1).sort(),els2=Object.keys(elCount2).sort();
    if(els1.length!==els2.length)return null;
    for(var k=0;k<els1.length;k++){
        if(els1[k]!==els2[k]||elCount1[els1[k]]!==elCount2[els2[k]])return null;
    }

    var deg1={},deg2={};
    for(var i=0;i<n1;i++){var d=adj1[i].length;deg1[d]=(deg1[d]||0)+1}
    for(var i=0;i<n2;i++){var d=adj2[i].length;deg2[d]=(deg2[d]||0)+1}
    var dk1=Object.keys(deg1).sort(),dk2=Object.keys(deg2).sort();
    if(dk1.length!==dk2.length)return null;
    for(var k=0;k<dk1.length;k++){if(dk1[k]!==dk2[k]||deg1[dk1[k]]!==deg2[dk2[k]])return null}

    // Precompute signatures and candidate lists
    var sig1=[],sig2=[];
    for(var i=0;i<n1;i++)sig1[i]=atomSignature(i,atoms1,adj1);
    for(var j=0;j<n2;j++)sig2[j]=atomSignature(j,atoms2,adj2);
    var candidates=[];
    for(var i=0;i<n1;i++){
        candidates[i]=[];
        for(var j=0;j<n2;j++){
            if(sig1[i]===sig2[j])candidates[i].push(j);
        }
        if(candidates[i].length===0)return null;
    }

    var map=new Array(n1).fill(-1);
    var rmap=new Array(n2).fill(-1);

    function consistent(a1,a2){
        if(atoms1[a1].element!==atoms2[a2].element)return false;
        var nb2set={};
        adj2[a2].forEach(function(e){nb2set[e.to]=true});
        for(var k=0;k<adj1[a1].length;k++){
            var nb1=adj1[a1][k].to;
            if(map[nb1]!==-1){
                if(nb2set[map[nb1]]===undefined)return false;
            }
        }
        var nb1set={};
        adj1[a1].forEach(function(e){nb1set[e.to]=true});
        for(var k=0;k<adj2[a2].length;k++){
            var nb2=adj2[a2][k].to;
            if(rmap[nb2]!==-1){
                if(nb1set[rmap[nb2]]===undefined)return false;
            }
        }
        return true;
    }

    function nextUnassigned(){
        var best=-1,bestKey=null;
        for(var i=0;i<n1;i++){
            if(map[i]!==-1)continue;
            var mappedNbrs=0;
            for(var k=0;k<adj1[i].length;k++){
                if(map[adj1[i][k].to]!==-1)mappedNbrs++;
            }
            var avail=0;
            for(var c=0;c<candidates[i].length;c++){
                if(rmap[candidates[i][c]]===-1)avail++;
            }
            var key=(mappedNbrs>0?1000000:0)-avail*1000+mappedNbrs;
            if(bestKey===null||key>bestKey){bestKey=key;best=i}
        }
        return best;
    }

    var iterCount=0;
    var maxIter=500000;

    function backtrack(depth){
        if(iterCount++>maxIter)return false;
        if(depth===n1)return true;
        var a1=nextUnassigned();
        if(a1===-1)return false;
        for(var c=0;c<candidates[a1].length;c++){
            var a2=candidates[a1][c];
            if(rmap[a2]!==-1)continue;
            if(!consistent(a1,a2))continue;
            map[a1]=a2;rmap[a2]=a1;
            if(backtrack(depth+1))return true;
            map[a1]=-1;rmap[a2]=-1;
        }
        return false;
    }

    if(backtrack(0))return map;
    return null;
}

function dist3d(a,b){var dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;return Math.sqrt(dx*dx+dy*dy+dz*dz)}
function angle3d(a,b,c){
    var v1={x:a.x-b.x,y:a.y-b.y,z:a.z-b.z},v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var d1=Math.sqrt(v1.x*v1.x+v1.y*v1.y+v1.z*v1.z),d2=Math.sqrt(v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(d1*d2);
    return Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI;
}
function dihedral3d(a,b,c,d){
    var v1={x:b.x-a.x,y:b.y-a.y,z:b.z-a.z};
    var v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var v3={x:d.x-c.x,y:d.y-c.y,z:d.z-c.z};
    var n1={x:v1.y*v2.z-v1.z*v2.y,y:v1.z*v2.x-v1.x*v2.z,z:v1.x*v2.y-v1.y*v2.x};
    var n2={x:v2.y*v3.z-v2.z*v3.y,y:v2.z*v3.x-v2.x*v3.z,z:v2.x*v3.y-v2.y*v3.x};
    var d1=Math.sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z),d2=Math.sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(n1.x*n2.x+n1.y*n2.y+n1.z*n2.z)/(d1*d2);
    var sign=v1.x*n2.x+v1.y*n2.y+v1.z*n2.z;
    var ang=Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI;
    return sign<0?-ang:ang;
}

// --- Group statistics for internal-coordinate differences -----------------
// Every compared bond / angle / dihedral - not just those above the display
// threshold - feeds an accumulator, so the reported rms / maxd / mean / std
// describe the whole compared set. Values are signed differences
// (left - right); a torsion's difference is folded into [-180, 180] so a
// chain and its reverse never disagree, while bond and angle differences
// cannot leave that range in the first place. rms is the root mean square of
// the raw differences, maxd the largest absolute difference, mean the signed
// average and std the population standard deviation. Entries are grouped by
// element composition key (e.g. C-C, C-C-H, C-C-C-H) - bonds by their two
// elements, angles by the two outer elements around the centre, torsions by
// the four elements along the chain, each canonicalised so that one internal
// coordinate always falls in the same group - and reported per group plus one
// combined "All" row.

function diffWrap180(v){
    while(v>180)v-=360;
    while(v<-180)v+=360;
    return v;
}

function newDiffStatAcc(){return{n:0,sum:0,sumSq:0,maxAbs:0}}

function collectDiffStat(map,key,delta){
    if(!isFinite(delta))return;
    var acc=map[key];
    if(!acc){acc=newDiffStatAcc();map[key]=acc}
    acc.n++;
    acc.sum+=delta;
    acc.sumSq+=delta*delta;
    var abs=delta<0?-delta:delta;
    if(abs>acc.maxAbs)acc.maxAbs=abs;
}

function finishDiffStat(acc){
    if(!acc||acc.n<=0)return{n:0,rms:0,maxd:0,mean:0,std:0};
    var n=acc.n;
    var mean=acc.sum/n;
    var msq=acc.sumSq/n;
    var variance=msq-mean*mean;
    return{n:n,rms:Math.sqrt(msq>0?msq:0),maxd:acc.maxAbs,mean:mean,std:Math.sqrt(variance>0?variance:0)};
}

function summarizeDiffStats(groupMap){
    var keys=Object.keys(groupMap).sort();
    var groups=[];
    var total=newDiffStatAcc();
    for(var i=0;i<keys.length;i++){
        var acc=groupMap[keys[i]];
        var g=finishDiffStat(acc);
        g.key=keys[i];
        groups.push(g);
        total.n+=acc.n;total.sum+=acc.sum;total.sumSq+=acc.sumSq;
        if(acc.maxAbs>total.maxAbs)total.maxAbs=acc.maxAbs;
    }
    return{groups:groups,total:finishDiffStat(total)};
}

function diffPairKey(e1,e2){return e1<=e2?e1+'-'+e2:e2+'-'+e1}

function diffTripletKey(e1,eCenter,e2){
    return e1<=e2?(e1+'-'+eCenter+'-'+e2):(e2+'-'+eCenter+'-'+e1);
}

// Canonical group key for a torsion. The two central elements keep the order
// they have along the chain and the whole quadruplet is reversed when the
// tail element sorts before the head (a torsion and its reverse are the same
// angle), so one torsion always lands in one group - e.g. C-C-C-H - no matter
// from which end the enumeration happened to reach it.
function diffQuadKey(e1,eCenter1,eCenter2,e4){
    return e1<=e4?(e1+'-'+eCenter1+'-'+eCenter2+'-'+e4):(e4+'-'+eCenter2+'-'+eCenter1+'-'+e1);
}

function computeConformationDiffs(atoms1,bonds1,atoms2,bonds2,mapping,threshPct){
    var diffs={bonds:[],angles:[],dihedrals:[]};
    var tp=threshPct||2;
    var statAcc={bonds:{},angles:{},dihedrals:{}};

    var bondSet2={};
    bonds2.forEach(function(b){bondSet2[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.order});
    bonds1.forEach(function(b){
        var a1m=Math.min(mapping[b.atom1],mapping[b.atom2]);
        var a2m=Math.max(mapping[b.atom1],mapping[b.atom2]);
        var key=a1m+'-'+a2m;
        if(bondSet2[key]!==undefined){
            var d1=dist3d(atoms1[b.atom1],atoms1[b.atom2]);
            var d2=dist3d(atoms2[mapping[b.atom1]],atoms2[mapping[b.atom2]]);
            collectDiffStat(statAcc.bonds,diffPairKey(atoms1[b.atom1].element,atoms1[b.atom2].element),d1-d2);
            var dd=Math.abs(d1-d2);
            var avg=(d1+d2)/2;
            var pct=avg>0.01?dd/avg*100:0;
            if(pct>tp){
                diffs.bonds.push({i1:b.atom1,i2:b.atom2,d1:d1,d2:d2,diff:dd,pct:pct});
            }
        }
    });

    var adj1=buildAdjacency(atoms1,bonds1);
    var angleSeen={};
    for(var j=0;j<atoms1.length;j++){
        var nbs=adj1[j].map(function(e){return e.to});
        for(var k=0;k<nbs.length;k++){
            for(var m=k+1;m<nbs.length;m++){
                var i1=nbs[k],i3=nbs[m];
                // Dedup key over the three atom indices, order-insensitive:
                // the pair is already visited once per centre, this only stops
                // a repeated neighbour entry from counting the same angle twice.
                // (A min/max key over the two outer atoms collides for two
                // different angles sharing the higher atom, silently dropping
                // angles whose centre has a smaller index than its neighbours.)
                var key=[i1,j,i3].sort().join('-');
                if(angleSeen[key])continue;
                angleSeen[key]=true;
                var a1=angle3d(atoms1[i1],atoms1[j],atoms1[i3]);
                var a2=angle3d(atoms2[mapping[i1]],atoms2[mapping[j]],atoms2[mapping[i3]]);
                collectDiffStat(statAcc.angles,diffTripletKey(atoms1[i1].element,atoms1[j].element,atoms1[i3].element),a1-a2);
                var da=Math.abs(a1-a2);
                var avgA=(a1+a2)/2;
                var pctA=avgA>1.0?da/avgA*100:da*100;
                if(pctA>tp){
                    diffs.angles.push({i1:i1,i2:j,i3:i3,a1:a1,a2:a2,diff:da,pct:pctA});
                }
            }
        }
    }

    var dihedralSeen={};
    bonds1.forEach(function(b){
        var a1=b.atom1,a2=b.atom2;
        [[a1,a2],[a2,a1]].forEach(function(pair){
            var start=pair[0],mid=pair[1];
            adj1[mid].forEach(function(e2){
                var end=e2.to;
                if(end===start)return;
                adj1[start].forEach(function(e1){
                    var s0=e1.to;
                    if(s0===mid)return;
                    var key=[s0,start,mid,end].sort().join('-');
                    if(dihedralSeen[key])return;
                    dihedralSeen[key]=true;
                    var d1=dihedral3d(atoms1[s0],atoms1[start],atoms1[mid],atoms1[end]);
                    var d2=dihedral3d(atoms2[mapping[s0]],atoms2[mapping[start]],atoms2[mapping[mid]],atoms2[mapping[end]]);
                    collectDiffStat(statAcc.dihedrals,
                        diffQuadKey(atoms1[s0].element,atoms1[start].element,atoms1[mid].element,atoms1[end].element),
                        diffWrap180(d1-d2));
                    var dd=Math.abs(d1-d2);
                    if(dd>180)dd=360-dd;
                    var avgD=(Math.abs(d1)+Math.abs(d2))/2;
                    var pctD=avgD>1.0?dd/avgD*100:dd*10;
                    if(pctD>tp){
                        diffs.dihedrals.push({i1:s0,i2:start,i3:mid,i4:end,d1:d1,d2:d2,diff:dd,pct:pctD});
                    }
                });
            });
        });
    });

    diffs.stats={
        bonds:summarizeDiffStats(statAcc.bonds),
        angles:summarizeDiffStats(statAcc.angles),
        dihedrals:summarizeDiffStats(statAcc.dihedrals)
    };

    return diffs;
}

function startDiff(msg){
    // The comparison has two structures to choose from, so a molecular-mechanics
    // panel has nothing to act on while it is open; it comes back on exit, and a
    // run that was still going is dropped rather than left writing coordinates
    // into the structure being compared.
    if(currentMode==='mmOpt'){
        var mmPanelInDiff=mmPanelEl();
        if(mmPanelInDiff)mmPanelInDiff.classList.remove('show');
    }
    mmCancelRun();
    diffData={atoms:msg.atoms,bonds:msg.bonds,title:msg.title||msg.fileName,fileName:msg.fileName,charge:msg.charge};

    var detBonds1=detectBondsFromAtoms(MD.atoms,MD.charge);
    var detBonds2=detectBondsFromAtoms(diffData.atoms,diffData.charge);
    diffDetBonds1=detBonds1;
    diffDetBonds2=detBonds2;

    var mapping=findAtomMapping(MD.atoms,detBonds1,diffData.atoms,detBonds2);

    if(!mapping){
        diffMapping=null;diffReverseMapping=null;
        showDiffPanel('<span class="diff-close">×</span>'+
            '<h4>Skeletons Differ</h4>'+
            '<div class="diff-row">Molecular skeletons are different — cannot compare.</div>'+
            '<div class="diff-row">Left: '+MD.atoms.length+' atoms, '+detBonds1.length+' bonds (detected)</div>'+
            '<div class="diff-row">Right: '+diffData.atoms.length+' atoms, '+detBonds2.length+' bonds (detected)</div>'+
            '<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Click Diff again to exit.</div>',true);
        diffPanelEl.classList.add('show');
        diffMode=true;
        diffLabelLeft.textContent='Original: '+(MD.title||'molecule');
        diffLabelLeft.classList.add('show');
        diffLabelRight.textContent='Diff: '+diffData.fileName+' (skeleton differs)';
        diffLabelRight.classList.add('show');
        modeInfoEl.textContent='Diff Mode (skeletons differ)';
        updateBakeHint();
        enterDiffRender();
        // The displacement slider acts on one structure; a comparison shows
        // two, so it is switched off (the session itself survives — the left
        // side is still the geometry the user displaced).
        vdSyncPanel();
        layoutPanels();
        return;
    }

    var diffs=computeConformationDiffs(MD.atoms,detBonds1,diffData.atoms,detBonds2,mapping,diffThresholdPct);
    var totalDiff=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;

    diffMapping=mapping;
    diffReverseMapping=new Array(diffData.atoms.length).fill(-1);
    for(var mi=0;mi<mapping.length;mi++){diffReverseMapping[mapping[mi]]=mi}

    var html=buildDiffResultsHTML(diffs,totalDiff);
    showDiffPanel(html,true);
    wireDiffSlider();
    diffMode=true;
    diffLabelLeft.textContent='Original: '+(MD.title||'molecule');
    diffLabelLeft.classList.add('show');
    diffLabelRight.textContent='Diff: '+diffData.fileName;
    diffLabelRight.classList.add('show');
    modeInfoEl.textContent='Diff Mode ('+totalDiff+' differences)';

    enterDiffRender(diffs,mapping);
    vdSyncPanel();
    layoutPanels();
}

function buildDiffResultsHTML(diffs,totalDiff){
    var html='<span class="diff-close">×</span>';
    html+='<h4>Diff Results</h4>';
    html+='<div class="diff-row" style="margin-bottom:4px">';
    html+='<label style="font-size:11px">Threshold: <b><span id="thresh-val">'+diffThresholdPct.toFixed(1)+'</span>%</b></label>';
    html+='<input type="range" id="thresh-slider" min="0" max="20" step="0.5" value="'+diffThresholdPct+'" style="width:100%;accent-color:var(--vscode-textLink-foreground,#3794ff)">';
    html+='<div style="font-size:10px;color:var(--vscode-descriptionForeground,#999)">Relative difference: Δ/avg × 100%</div>';
    html+='</div>';
    html+='<div id="diff-results-list">';
    html+=buildDiffListHTML(diffs,totalDiff);
    html+='</div>';
    return html;
}

function buildDiffListHTML(diffs,totalDiff){
    var html='';
    if(totalDiff===0){
        html+='<div class="diff-row" style="color:#4ec9b0">Structures are identical (same conformation within thresholds).</div>';
    }else{
        html+='<div class="diff-row">Mapping: '+MD.atoms.length+' atoms matched. '+totalDiff+' differences found.</div>';
    }
    html+=buildDiffStatsHTML(diffs.stats);
    if(totalDiff>0){
        if(diffs.bonds.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Length Differences ('+diffs.bonds.length+')</b></div>';
            diffs.bonds.sort(function(a,b){return b.pct-a.pct});
            diffs.bonds.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+
                    ': '+d.d1.toFixed(3)+' vs '+d.d2.toFixed(3)+' Å (Δ='+d.diff.toFixed(3)+', '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.angles.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Angle Differences ('+diffs.angles.length+')</b></div>';
            diffs.angles.sort(function(a,b){return b.pct-a.pct});
            diffs.angles.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+'-'+MD.atoms[d.i3].element+(d.i3+1)+
                    ': '+d.a1.toFixed(1)+'° vs '+d.a2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+'°, '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.dihedrals.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Dihedral Differences ('+diffs.dihedrals.length+')</b></div>';
            diffs.dihedrals.sort(function(a,b){return b.pct-a.pct});
            diffs.dihedrals.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+'-'+MD.atoms[d.i3].element+(d.i3+1)+'-'+MD.atoms[d.i4].element+(d.i4+1)+
                    ': '+d.d1.toFixed(1)+'° vs '+d.d2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+'°, '+d.pct.toFixed(2)+'%)</div>';
            });
        }
    }
    return html;
}

function diffStatNum(v,decimals){
    if(!isFinite(v))return '-';
    return v.toFixed(decimals);
}

function buildDiffStatsHTML(stats){
    if(!stats)return '';
    var html='<div class="diff-stats">';
    html+='<div class="diff-stats-head"><span class="diff-stats-arrow">'+(diffStatsExpanded?'▼':'▶')+'</span> Internal Coordinate Statistics</div>';
    html+='<div class="diff-stats-body">';
    html+='<div class="diff-stat-note">All compared internal coordinates (signed Δ = left − right), independent of the threshold. rms = √mean(Δ²), maxd = max|Δ|, mean = mean(Δ), std = population standard deviation.</div>';
    html+=buildDiffStatTableHTML('Bond Lengths',stats.bonds,'Å',4);
    html+=buildDiffStatTableHTML('Bond Angles',stats.angles,'°',3);
    html+=buildDiffStatTableHTML('Dihedrals',stats.dihedrals,'°',3);
    html+='</div></div>';
    return html;
}

function buildDiffStatTableHTML(title,report,unit,decimals){
    if(!report||!report.groups||!report.total)return '';
    if(report.total.n===0){
        return '<div class="diff-stat-title">'+title+': no comparable coordinates</div>';
    }
    var html='<div class="diff-stat-title">'+title+' ('+unit+') — n='+report.total.n+'</div>';
    html+='<div class="diff-stat-scroll"><table class="diff-stat-table">';
    html+='<tr><th>Group</th><th>n</th><th>rms</th><th>maxd</th><th>mean</th><th>std</th></tr>';
    for(var i=0;i<report.groups.length;i++){
        var g=report.groups[i];
        html+='<tr><td>'+g.key+'</td><td>'+g.n+'</td><td>'+diffStatNum(g.rms,decimals)+'</td><td>'+diffStatNum(g.maxd,decimals)+
            '</td><td>'+diffStatNum(g.mean,decimals)+'</td><td>'+diffStatNum(g.std,decimals)+'</td></tr>';
    }
    var t=report.total;
    html+='<tr class="diff-stat-total"><td>All</td><td>'+t.n+'</td><td>'+diffStatNum(t.rms,decimals)+'</td><td>'+diffStatNum(t.maxd,decimals)+
        '</td><td>'+diffStatNum(t.mean,decimals)+'</td><td>'+diffStatNum(t.std,decimals)+'</td></tr>';
    html+='</table></div>';
    return html;
}

// The statistics block is collapsible: the panel can hold both the detailed
// difference list and the grouped statistics without one hiding the other.
// The expanded/collapsed choice is kept in diffStatsExpanded, so it survives
// the list re-render that every threshold-slider move triggers.
function wireDiffStatsToggle(rootEl){
    if(!rootEl)return;
    var wrap=rootEl.querySelector('.diff-stats');
    if(!wrap)return;
    var head=wrap.querySelector('.diff-stats-head');
    var arrow=wrap.querySelector('.diff-stats-arrow');
    if(diffStatsExpanded)wrap.classList.remove('collapsed');
    else wrap.classList.add('collapsed');
    if(!head)return;
    head.addEventListener('click',function(){
        diffStatsExpanded=!diffStatsExpanded;
        if(diffStatsExpanded)wrap.classList.remove('collapsed');
        else wrap.classList.add('collapsed');
        if(arrow)arrow.textContent=diffStatsExpanded?'▼':'▶';
        layoutPanels();
    });
}

function wireDiffSlider(){
    var slider=diffPanelEl.querySelector('#thresh-slider');
    if(slider){
        slider.addEventListener('input',function(){
            diffThresholdPct=parseFloat(slider.value);
            var valEl=diffPanelEl.querySelector('#thresh-val');
            if(valEl)valEl.textContent=diffThresholdPct.toFixed(1);
            recomputeDiff();
        });
    }
}

function recomputeDiff(){
    if(!diffMapping||!diffDetBonds1||!diffDetBonds2)return;
    var diffs=computeConformationDiffs(MD.atoms,diffDetBonds1,diffData.atoms,diffDetBonds2,diffMapping,diffThresholdPct);
    var totalDiff=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;
    var listEl=diffPanelEl.querySelector('#diff-results-list');
    if(listEl){
        listEl.innerHTML=buildDiffListHTML(diffs,totalDiff);
        wireDiffStatsToggle(listEl);
    }
    var fullHTML=buildDiffResultsHTML(diffs,totalDiff);
    diffPanelHTML=fullHTML;
    diffResultsHTML=fullHTML;
    modeInfoEl.textContent='Diff Mode ('+totalDiff+' differences)';
    updateDiffHighlights(diffs);
}
function refreshDiff(){
    if(!diffMode)return;
    diffDetBonds1=detectBondsFromAtoms(MD.atoms,MD.charge);
    diffMapping=findAtomMapping(MD.atoms,diffDetBonds1,diffData.atoms,diffDetBonds2);
    if(diffMapping){
        diffReverseMapping=new Array(diffData.atoms.length).fill(-1);
        for(var mi=0;mi<diffMapping.length;mi++){diffReverseMapping[diffMapping[mi]]=mi}
        recomputeDiff();
    }else{
        diffReverseMapping=null;
        diffPanelEl.innerHTML='<span class="diff-close">×</span>'+
            '<h4>Skeletons Differ</h4>'+
            '<div class="diff-row">After edit, molecular skeletons differ — cannot compare.</div>'+
            '<div class="diff-row">Left: '+MD.atoms.length+' atoms, '+diffDetBonds1.length+' bonds (detected)</div>'+
            '<div class="diff-row">Right: '+diffData.atoms.length+' atoms, '+diffDetBonds2.length+' bonds (detected)</div>';
        diffLabelRight.textContent='Diff: '+diffData.fileName+' (skeleton differs)';
        atomMeshes.forEach(function(m){m.userData.diffHi=false;setAtomHighlight(m,0)});
        diffAtomMeshes.forEach(function(m){m.userData.diffHi=false;setDiffAtomHighlight(m,0)});
        diffBondPieces.forEach(function(p){p.hi=false});
        diffBuildBondPools();
        needsRender=true;
    }
}

function updateDiffHighlights(diffs){
    var highlightAtoms={};
    var highlightBonds={};
    diffs.bonds.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;
        highlightBonds[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true;
    });
    diffs.angles.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;
    });
    diffs.dihedrals.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;highlightAtoms[d.i4]=true;
    });

    atomMeshes.forEach(function(m,i){
        m.userData.diffHi=highlightAtoms[i]?true:false;
    });

    diffAtomMeshes.forEach(function(m,i){
        var origI=diffReverseMapping?diffReverseMapping[i]:-1;
        m.userData.diffHi=(origI>=0&&highlightAtoms[origI])?true:false;
    });

    // A highlighted bond moves to the pool whose material carries the orange
    // glow; the two pools are refilled from the piece list, which is a pass
    // over the instances, not a material edit per mesh.
    diffBondPieces.forEach(function(p){
        if(p.key)p.hi=highlightBonds[p.key]?true:false;
    });
    diffBuildBondPools();

    highlightSelected();
}

function showDiffPanel(html,isResults){
    diffPanelHTML=html;
    if(isResults)diffResultsHTML=html;
    diffPanelEl.innerHTML=html;
    diffPanelEl.classList.add('show');
    diffReopenEl.classList.remove('show');
    var closeBtn=diffPanelEl.querySelector('.diff-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            diffPanelEl.classList.remove('show');
            diffReopenEl.classList.add('show');
            layoutPanels();
        });
    }
    wireDiffSlider();
    wireDiffStatsToggle(diffPanelEl);
    layoutPanels();
}

diffReopenEl.addEventListener('click',function(){
    var html=diffResultsHTML||diffPanelHTML;
    diffPanelEl.innerHTML=html;
    diffPanelEl.classList.add('show');
    diffReopenEl.classList.remove('show');
    var closeBtn=diffPanelEl.querySelector('.diff-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            diffPanelEl.classList.remove('show');
            diffReopenEl.classList.add('show');
            layoutPanels();
        });
    }
    wireDiffSlider();
    wireDiffStatsToggle(diffPanelEl);
    layoutPanels();
});

function enterDiffRender(diffs,mapping){
    if(vibActive||vibPaused)stopVibration();
    stopMoveSession();
    // A comparison has two structures and two view states, so the panel stops
    // advertising Alt coordinate editing for as long as it is up.
    updateBakeHint();
    diffPivot=new THREE.Group();scene.add(diffPivot);
    diffMolGroup=new THREE.Group();diffPivot.add(diffMolGroup);

    diffCX=0;diffCY=0;diffCZ=0;
    diffData.atoms.forEach(function(a){diffCX+=a.x;diffCY+=a.y;diffCZ+=a.z});
    if(diffData.atoms.length>0){diffCX/=diffData.atoms.length;diffCY/=diffData.atoms.length;diffCZ/=diffData.atoms.length}

    var highlightAtoms={};
    var highlightBonds={};
    if(diffs&&mapping){
        diffs.bonds.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;
            highlightBonds[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true;
        });
        diffs.angles.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;
        });
        diffs.dihedrals.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;highlightAtoms[d.i4]=true;
        });
    }

    var dSegs=sphereSegs(diffData.atoms.length);
    diffData.atoms.forEach(function(a,i){
        var r=getR(a.element);
        var isHi=highlightAtoms[mapping?mapping.indexOf(i):-1];
        // Same pooled geometry and material scheme as the left-hand structure:
        // a comparison of two large structures draws its two sides from a
        // handful of objects instead of one per atom.
        var mesh=new THREE.Mesh(sphereGeoFor(r,dSegs),diffAtomMaterial(a.color,isHi?2:0));
        mesh.position.set(a.x-diffCX,a.y-diffCY,a.z-diffCZ);
        mesh.userData={element:a.element,index:i,diffHi:isHi?true:false,baseColor:a.color};
        diffMolGroup.add(mesh);
        diffAtomMeshes.push(mesh);
    });

    // The right-hand bonds: two pieces per bond (one per half), each carrying
    // the transform it is drawn with. The comparison's structure never moves,
    // so a piece stores its own numbers and the highlight can be re-applied by
    // refilling the two pools rather than editing a material per mesh.
    diffBondPieces.length=0;
    diffData.bonds.forEach(function(b){
        var a1=diffData.atoms[b.atom1],a2=diffData.atoms[b.atom2];
        if(!a1||!a2)return;
        var sx=a1.x-diffCX,sy=a1.y-diffCY,sz=a1.z-diffCZ;
        var ex=a2.x-diffCX,ey=a2.y-diffCY,ez=a2.z-diffCZ;
        var dx=ex-sx,dy=ey-sy,dz=ez-sz;
        var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
        if(l<0.001)return;
        var mpx=(sx+ex)*0.5,mpy=(sy+ey)*0.5,mpz=(sz+ez)*0.5;
        var origI1=mapping?mapping.indexOf(b.atom1):-1;
        var origI2=mapping?mapping.indexOf(b.atom2):-1;
        var bondKey=Math.min(origI1,origI2)+'-'+Math.max(origI1,origI2);
        var isHi=highlightBonds[bondKey];
        addDiffBondPiece(sx,sy,sz,mpx-sx,mpy-sy,mpz-sz,l/2,0.12,new THREE.Color(a1.color),bondKey,isHi);
        addDiffBondPiece(mpx,mpy,mpz,ex-mpx,ey-mpy,ez-mpz,l/2,0.12,new THREE.Color(a2.color),bondKey,isHi);
    });
    diffBuildBondPools();

    if(diffs&&mapping){
        var leftHi={};
        diffs.bonds.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true});
        diffs.angles.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true;leftHi[d.i3]=true});
        diffs.dihedrals.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true;leftHi[d.i3]=true;leftHi[d.i4]=true});
        atomMeshes.forEach(function(m,i){
            m.userData.diffHi=leftHi[i]?true:false;
        });
        // The flags above are only half of it — the appearance has to be put on
        // the atoms now, because nothing else runs on this path: the comparison
        // is entered without a recompute, so the highlight would otherwise wait
        // for the next click to appear.
        highlightSelected();
    }

    var allAtoms=MD.atoms.concat(diffData.atoms);
    var maxD2=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD2)maxD2=dd});
    diffData.atoms.forEach(function(a){var dx=a.x-diffCX,dy=a.y-diffCY,dz=a.z-diffCZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD2)maxD2=dd});
    var newCam=maxD2*2.5+5;
    camDist=newCam;diffCamDist=newCam;
    camera.position.set(0,0,camDist);
    panX=0;panY=0;
    diffPanX=0;diffPanY=0;
    diffRotQuat.copy(rotQuat);
    diffTransformSide='left';
    updateTransform();
}

function exitDiff(){
    stopSpin();
    diffMode=false;
    updateBakeHint();
    // One structure again: the displacement block takes its mode back.
    vdSyncPanel();
    // A molecular-mechanics session survives the comparison: its panel comes
    // back with the same selection once the second structure is gone.
    if(currentMode==='mmOpt')showMmPanel();
    if(diffPivot){
        scene.remove(diffPivot);
        diffPivot.traverse(function(ch){
            // The comparison draws with the same pooled geometries/materials
            // the structure uses, so releasing them here would break the side
            // that is about to be shown again — only what this mesh owns goes.
            disposeMesh(ch);
        });
        diffPivot=null;diffMolGroup=null;
    }
    poolClear(diffBondPoolPlain);
    poolClear(diffBondPoolHot);
    diffBondPieces.length=0;
    diffAtomMeshes=[];
    diffSelectedAtoms=[];
    diffMapping=null;diffReverseMapping=null;
    diffDetBonds1=null;diffDetBonds2=null;
    diffActiveSide='left';
    diffTransformSide='left';
    diffRotQuat.identity();
    diffPanX=0;diffPanY=0;diffCamDist=10;
    diffPanelEl.classList.remove('show');
    diffPanelEl.innerHTML='';
    diffReopenEl.classList.remove('show');
    diffPanelHTML='';
    diffResultsHTML='';
    diffLabelLeft.classList.remove('show');
    diffLabelRight.classList.remove('show');
    layoutPanels();
    modeInfoEl.textContent='View Mode';
    atomMeshes.forEach(function(m){
        m.userData.diffHi=false;
        setAtomHighlight(m,0);
    });
    highlightSelected();
    pivotGroup.visible=true;
    var w=container.clientWidth||window.innerWidth;
    var h=container.clientHeight||(window.innerHeight-60);
    if(w<1)w=window.innerWidth;
    if(h<1)h=window.innerHeight-60;
    renderer.setScissorTest(false);
    renderer.setViewport(0,0,w,h);
    renderer.setScissor(0,0,w,h);
    renderer.autoClear=true;
    camera.aspect=w/h;
    camera.updateProjectionMatrix();
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    camDist=maxD*2.5+5;
    camera.position.set(0,0,camDist);
    updateTransform();
}

function camDiffInit(){
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    if(diffData){
        diffData.atoms.forEach(function(a){var dx=a.x-diffCX,dy=a.y-diffCY,dz=a.z-diffCZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    }
    return maxD*2.5+5;
}

// ===== Scalar (Beta) colouring =====
// PDB column 61-66 is not a B-factor in this toolchain: it carries an arbitrary
// per-atom scalar (unit index, per-atom force error, …). Any structure that has
// one can be painted by it, which is how the force-error and unit-index files
// are meant to be looked at. Must stay identical to vmdVisualization.scalarColor.
// The gradient is not a plain red-green-blue ramp: the two passive channels sit
// on a 15/255 floor and the active channel ramps at twice the nominal rate,
// which is what a rendered sphere really shows in the reference renderer.
// These two constants must stay above the startup code that colours a bare
// scalar file at load — declared further down they would still be undefined at
// that point, and every atom would come out black.
var GRADIENT_BASE=15/255;
var GRADIENT_RATE=2.0;
function clamp01(v){
    if(!isFinite(v))return 0;
    return v<0?0:v>1?1:v;
}
function scalarColor(v,vmin,vmax){
    var t;
    if(!isFinite(vmin)||!isFinite(vmax)||vmax===vmin){t=0}
    else{
        t=(v-vmin)/(vmax-vmin);
        if(!isFinite(t))t=0;
        t=clamp01(t);
    }
    var r,g,b;
    if(t<=0.5){
        r=clamp01(GRADIENT_BASE+GRADIENT_RATE*(0.5-t));
        g=clamp01(GRADIENT_BASE+GRADIENT_RATE*t);
        b=GRADIENT_BASE;
    }else{
        r=GRADIENT_BASE;
        g=clamp01(GRADIENT_BASE+GRADIENT_RATE*(1-t));
        b=clamp01(GRADIENT_BASE+GRADIENT_RATE*(t-0.5));
    }
    var hex=function(x){var n=Math.round(x*255);if(!isFinite(n))n=0;if(n<0)n=0;if(n>255)n=255;var s=n.toString(16).toUpperCase();return s.length<2?'0'+s:s};
    return '#'+hex(r)+hex(g)+hex(b);
}
var betaColorOn=false;
// 配色使用的范围: 按钮开启时取数据自身的范围。重建场景 (编辑、撤销、导入、
// 换帧、换分子) 时按原子现算颜色, 所以它是纯视图状态 —— 没有需要备份/还原
// 的东西, 组成怎么变都不会错位。
var betaRangeUsed=null;
function betaRange(){
    var lo=Infinity,hi=-Infinity,any=false;
    MD.atoms.forEach(function(a){
        var b=a.beta;
        if(typeof b==='number'&&isFinite(b)){any=true;if(b<lo)lo=b;if(b>hi)hi=b}
    });
    return any&&hi>lo?{lo:lo,hi:hi}:null;
}
function updateBetaButton(){
    var btn=document.getElementById('beta-btn');
    if(!btn)return;
    var ok=!VMD&&!CRY&&betaRange()!==null;
    btn.style.display=ok?'':'none';
}
// 开启配色时某个原子的标量颜色; 不带标量的原子返回 null (保留元素颜色)。
// 在重建场景时读取, 而不是把颜色写进 a.color —— a.color 始终是元素颜色。
function betaColorOf(a){
    if(!betaColorOn||!betaRangeUsed)return null;
    var b=a.beta;
    return (typeof b==='number'&&isFinite(b))?scalarColor(b,betaRangeUsed.lo,betaRangeUsed.hi):null;
}
function applyBetaColoring(on){
    if(on){
        var range=betaRange();
        if(!range)return;
        betaRangeUsed=range;
    }
    betaColorOn=on&&!!betaRangeUsed;
    var btn=document.getElementById('beta-btn');
    if(btn)btn.classList.toggle('active',betaColorOn);
    rebuildScene();
    modeInfoEl.textContent=betaColorOn
        ?('Beta colouring: '+betaRangeUsed.lo.toFixed(3)+' – '+betaRangeUsed.hi.toFixed(3))
        :(MODE_INFO[currentMode]||currentMode);
}
(function(){
    var btn=document.getElementById('beta-btn');
    if(btn)btn.addEventListener('click',function(){applyBetaColoring(!betaColorOn)});
    var prev=document.getElementById('prev-mol');
    var next=document.getElementById('next-mol');
    if(prev)prev.addEventListener('click',function(){switchMolecule(vmdMolIdx-1)});
    if(next)next.addEventListener('click',function(){switchMolecule(vmdMolIdx+1)});
})();

var currentFrame=0;
var totalFrames=MD.frames?MD.frames.length:0;
var autoPlayTimer=null;
var isAutoPlaying=false;
// How long each structure of a trajectory is held while it plays. A frame is a
// still, not a film frame: every step is a whole scene rebuild, and it has to
// stay on screen long enough to be read as one structure — but short enough
// that a run of a few hundred steps is watchable rather than a slideshow.
var FRAME_PLAY_MS=350;
var frameNavEl=document.getElementById('frame-nav');
var frameSepEl=document.getElementById('frame-sep');
if(totalFrames>1){
    frameNavEl.classList.add('show');
    frameSepEl.classList.remove('hidden');
    updateFrameInfo();
}else{
    frameNavEl.classList.remove('show');
    frameSepEl.classList.add('hidden');
}
// Molecule selector of a VMD script (one entry per "mol new"), and the Beta
// colouring toggle for structures that carry a per-atom scalar.
if(VMD){
    setVmdNavEnabled(VMD.molecules.length>1&&!CRY);
    updateMolLabel();
    var vmdWarn=(VMD.warnings||[]).concat((VMD.molecules[vmdMolIdx]&&VMD.molecules[vmdMolIdx].warnings)||[]);
    if(vmdWarn.length>0){
        modeInfoEl.textContent='TCL: '+vmdWarn[0]+(vmdWarn.length>1?' (+'+(vmdWarn.length-1)+' more)':'');
    }
}
updateBetaButton();
// A visualization that ships as a bare PDB rather than a script: paint its
// per-atom scalar straight away, so opening the file shows the picture it was
// written for. The toolbar button stays active and turns it back off.
if(MD.autoBeta)applyBetaColoring(true);
function updateFrameInfo(){
    var el=document.getElementById('frame-info');
    var numEl=document.getElementById('frame-num');
    if(el)el.textContent='/'+totalFrames+(MD.frames[currentFrame]?' - '+MD.frames[currentFrame].stepLabel:'');
    if(numEl)numEl.value=currentFrame+1;
}
function switchFrame(idx){
    if(idx<0||idx>=totalFrames)return;
    // A frame is a structure of its own, so anything still writing the one on
    // screen has to let go before it is replaced: a normal-mode playback
    // rewrites every atom on each animation frame, and its own stop puts the
    // geometry it started from back — done after the swap, that would write the
    // previous frame's coordinates over the frame just loaded.
    if(vibActive||vibPaused)stopVibration();
    stopMoveSession();
    mmCancelRun();
    // A frame is a structure of its own: the equilibrium a displacement session
    // was measured from belongs to the frame being left.
    vdEnd();
    currentFrame=idx;
    var f=MD.frames[idx];
    MD.atoms=f.atoms.map(function(a,i){a.index=i;return a});
    MD.bonds=f.bonds||[];
    if(MD.bonds.length===0){
        MD.bonds=detectBondsFromAtoms(MD.atoms,MD.charge);
    }
    if(VMD){
        // A frame is a structure of its own; the script's recipe was resolved
        // against the coordinates it was loaded from, so its surfaces and its
        // molecule list no longer describe what is on screen.
        VMD.activeSegments=[];
        setVmdNavEnabled(false);
    }
    undoStack.length=0;updateUndoBtn();
    if(diffMode)resetSelection();
    resetRotAxisState();
    if(currentMode==='mmOpt'){mmLastResult=null;refreshMmPanel()}
    rebuildScene();
    updateFrameInfo();
    // [2026-09-12 | Atreides-Jimmy] Follow the structure change in the
    // convergence charts: the point of the step now on screen gets highlighted
    // (renderOptCharts is a no-op while the panel is closed).
    if(typeof renderOptCharts==='function')renderOptCharts();
    // A frame is a structure of its own: whatever was done to the one before
    // is not part of this one, so the saved title starts from the frame again.
    provReset();
}

// ===== VMD molecules =====
// A script may create several molecules ("mol new" repeated). They are shown one
// at a time rather than superimposed, because that is what makes them
// comparable — the coordinates are of the same system at different stages.
function setVmdNavEnabled(on){
    var sep=document.getElementById('mol-sep');
    var nav=document.getElementById('mol-nav');
    if(sep)sep.classList.toggle('hidden',!on);
    if(nav)nav.classList.toggle('show',!!on);
}
function updateMolLabel(){
    var el=document.getElementById('mol-label');
    if(!el)return;
    var m=VMD.molecules[vmdMolIdx];
    el.textContent=(vmdMolIdx+1)+'/'+VMD.molecules.length+' '+(m?(m.label||''):'');
}
function switchMolecule(idx){
    if(!VMD||idx<0||idx>=VMD.molecules.length||idx===vmdMolIdx)return;
    if(diffMode){modeInfoEl.textContent='Molecule switching is not available in diff mode';return}
    if(CRY){modeInfoEl.textContent='Molecule switching is not available for crystal structures';return}
    // Same rule as a frame change: every writer of the current structure ends
    // before it is replaced (see switchFrame), a displacement session included.
    if(vibActive||vibPaused)stopVibration();
    stopMoveSession();
    mmCancelRun();
    vdEnd();
    vmdMolIdx=idx;
    var m=VMD.molecules[idx];
    MD.atoms=m.atoms;
    MD.bonds=m.bonds;
    MD.atoms.forEach(function(a,i){a.index=i});
    MD.bonds.forEach(function(b,i){b.index=i});
    VMD.activeSegments=m.segments||[];
    // 新分子的新补丁列表: 面板选择全部复位 (初始时所有 egg 都显示)。
    vmdHiddenSegments={};
    vmdSurfaceSignature='';
    undoStack.length=0;updateUndoBtn();
    resetRotAxisState();
    if(currentMode==='mmOpt'){mmLastResult=null;refreshMmPanel()}
    rebuildScene();
    updateMolLabel();
    // Frame the new structure: it may sit at a very different position or size.
    var d=computeInitialCamDist();
    initCam=d;camDist=d;
    camera.position.set(0,0,d);
    rotQuat.identity();panX=0;panY=0;
    updateTransform();
    modeInfoEl.textContent=MODE_INFO[currentMode]||currentMode;
    provReset();
}

// ===== Bond perception (client side) =====
// The extension host serialises its bond-perception module into this page as
// BP, so the viewport re-derives bonds with exactly the code the parser used:
// one implementation of the tables, the order assignment and the valence pass,
// in one place. That is what makes a re-derivation reproducible — a move lands
// on the bonds a fresh load of the moved coordinates gives, because nothing
// here can disagree with what produced them.
// The total charge is passed in rather than read off MD here: a re-detection
// does not always concern the structure on screen (the diff's comparison side
// is a different file with its own charge, and it has none, so it gets the
// neutral assumption a bare coordinate set would).
function detectBondsFromAtoms(atoms,charge){
    return BP.perceiveBonds(atoms,(typeof charge==='number')?{charge:charge}:undefined).bonds;
}
// Stopping auto-play is its own step because more than the button needs it: a
// frame auto-play replaces the structure on every tick, and any operation that
// needs the structure to hold still — a Molecular Mechanics run, an edit, a
// normal-mode playback — has to end it first rather than be cancelled by it a
// moment later.
function stopFrameAutoPlay(){
    if(!isAutoPlaying)return;
    if(autoPlayTimer){clearInterval(autoPlayTimer);autoPlayTimer=null}
    isAutoPlaying=false;
    var btn=document.getElementById('auto-play');
    if(btn){btn.textContent='⏵ Play';btn.classList.remove('playing')}
}
document.getElementById('prev-frame').addEventListener('click',function(){if(currentFrame>0)switchFrame(currentFrame-1)});
document.getElementById('next-frame').addEventListener('click',function(){if(currentFrame<totalFrames-1)switchFrame(currentFrame+1)});
document.getElementById('frame-num').addEventListener('change',function(){var n=parseInt(this.value);if(!isNaN(n)&&n>=1&&n<=totalFrames)switchFrame(n-1)});
document.getElementById('auto-play').addEventListener('click',function(){
    if(isAutoPlaying){
        stopFrameAutoPlay();
    }else{
        // Playback replaces the structure on every tick, so everything else
        // writing it ends now rather than on the first tick — the same rule the
        // other starters follow, applied at the moment the user asks for it.
        if(vibActive||vibPaused)stopVibration();
        stopMoveSession();
        mmCancelRun();
        // Playback stops on the last structure instead of wrapping around, so
        // asking for it while that structure is on screen means "run the
        // trajectory again": restart from the first one. Without this the
        // button would do nothing at all there.
        if(totalFrames>1&&currentFrame>=totalFrames-1)switchFrame(0);
        isAutoPlaying=true;this.textContent='⏸ Stop';this.classList.add('playing');
        autoPlayTimer=setInterval(function(){
            var next=currentFrame+1;
            // The last structure is where a trajectory ends, and that is
            // usually the one worth stopping to look at — the geometry the run
            // converged to. Wrapping back to the first one would hide it, so
            // playback stops here and the button returns to Play.
            if(next>=totalFrames){stopFrameAutoPlay();return}
            switchFrame(next);
        },FRAME_PLAY_MS);
    }
});

// ===== Opt Convergence & Frequency Vibration Panels =====
var optPanelEl=document.getElementById('opt-panel');
var optReopenEl=document.getElementById('opt-reopen');
var freqPanelEl=document.getElementById('freq-panel');
var freqReopenEl=document.getElementById('freq-reopen');
var routePanelEl=document.getElementById('route-panel');
var routeReopenEl=document.getElementById('route-reopen');
var OPT_STEPS=MD.optSteps||null;
var OPT_SOURCE=MD.optSource||null;
// The convergence targets to draw as dashed guides. The file's own values win:
// ORCA states its criteria (TolMAXG / TolRMSG / TolMAXD / TolRMSD) in every
// optimization output, and a run told to converge tighter than the default must
// show the line it was actually held to. Gaussian only prints its criteria as
// the Tolerance column of the same table, so its documented defaults stand in
// (force 0.000450/0.000300 Hartree/Bohr, displacement 0.001800/0.001200 Bohr),
// and they apply ONLY to a Gaussian log — xtb's gnorm criterion and a bare
// trajectory state no tolerance at all, and a borrowed scale would draw a target
// the run was never asked to meet.
var OPT_THRESHOLDS=(function(){
    var t=MD.optThresholds;
    if(t&&(t.maxForce!=null||t.rmsForce!=null||t.maxDisplacement!=null||t.rmsDisplacement!=null))return t;
    if(OPT_SOURCE==='gaussian')return {maxForce:0.000450,rmsForce:0.000300,maxDisplacement:0.001800,rmsDisplacement:0.001200};
    return null;
})();
var NORMAL_MODES=MD.normalModes||null;

// ---------------------------------------------------------------------------
// Convergence chart value formatting, axis labels, point highlight and hover
// tooltip.
// [2026-09-12 | Atreides-Jimmy] Added:
//   * feature 1 — energy axis in normal decimal notation: the tick range is
//     derived from the highest/lowest plotted value and the labels keep only the
//     changing decimals, so a run from -47077.6380 to -47077.6713 reads
//     .64 .65 .65 .66 .67 (base -47077 in the title) instead of 4.71e+04.
//   * feature 2 — the point of the optimization step currently on screen is
//     marked with a ring + vertical guide.
//   * feature 3 — hovering a point reports its exact value in a tooltip.
// ---------------------------------------------------------------------------

/** Decimal places that resolve a value on its own: two significant digits of
 *  its magnitude plus one guard digit, capped at 12 (12.16 -> 2, 0.0000451 -> 7,
 *  a value around 6912 -> 0). Used as a floor for the hover tooltip, which must
 *  stay readable even for a single-point series with no span to measure. */
function convValueDecimals(v){
    var mag=Math.abs(v);
    var dec=mag>0?Math.max(0,2-Math.floor(Math.log10(mag))):6;
    if(dec>12)dec=12;
    return dec;
}

/** Format v with an explicit number of decimals, trimming trailing zeros so the
 *  tooltip reads '-6912.8' rather than '-6912.80', and never using scientific
 *  notation. A null dec falls back to the value's own magnitude. */
function convFormat(v,dec){
    if(v==null||isNaN(v))return String(v);
    if(dec==null)dec=convValueDecimals(v);
    var str=v.toFixed(dec);
    if(str.indexOf('.')>=0)str=str.replace(/0+$/,'').replace(/\.$/,'');
    if(str==='-0')str='0';
    return str;
}

/** A value in plain decimal notation with just enough decimals to resolve its
 *  own magnitude (kept under its historical name for the linear-axis fallback). */
function convFullValue(v){return convFormat(v,convValueDecimals(v))}

/** Decimal places needed so that two numbers one tick-step apart never print
 *  the same label: with dec decimals the finest difference a label can show
 *  is 10^-dec, so we keep two significant digits of the STEP. The step is
 *  span/4 (five gridlines) and is taken from the plotted range, NOT from whether
 *  the run's two ends share an integer part — a wide sub-integer run that
 *  crosses an integer boundary (e.g. -6911.977 … -6912.806) used to be forced
 *  to 0 decimals and collapsed the axis to '0 0 0 1 1'.
 *  [2026-09-13 | Atreides-Jimmy] Replaced the "first differing decimal of the
 *  two ends" rule with this step-based one. */
function convStepDecimals(span){
    span=Math.abs(span);
    if(!(span>0)||!isFinite(span))return 6;
    var step=span/4;
    var d=1-Math.floor(Math.log10(step));
    if(d<0)d=0;
    if(d>12)d=12;
    return d;
}

/** Number of decimals for an axis spanning [vmin,vmax] (see convStepDecimals). */
function convTickDecimals(vmin,vmax){return convStepDecimals(vmax-vmin)}

/** Axis tick label: only the part that actually varies over the axis is kept —
 *  the shared integer part (base) is printed once next to the chart title, so
 *  ticks read '.67' / '.66' / … for a run around -47077. A null base falls back
 *  to the full value — always as plain decimals, never scientific notation.
 *  [2026-09-13 | Atreides-Jimmy] A tick landing exactly on the base now prints
 *  '0' (the old code stripped '0.000' down to a bare '000'). */
function convTickLabel(v,base,dec){
    if(base==null)return v.toFixed(dec);
    var frac=Math.abs(v-base).toFixed(dec);
    if(/^0(\.0+)?$/.test(frac))return '0';
    // Drop the leading zero of a pure fraction so the ticks read '.9 .8 .7'
    // (the user-visible convention of this panel) rather than '0.9 0.8 0.7'.
    if(frac.indexOf('0.')===0)frac=frac.substring(1);
    return frac;
}

/** Base (the integer whose fraction the axis prints) chosen so that every tick
 *  lies on ONE side of it, which keeps the printed fraction monotonic from top
 *  to bottom:
 *    - span ≥ 1        → null (the ticks carry their full value),
 *    - all values ≤ 0  → ceil(vmax)  (base at/above the data),
 *    - all values ≥ 0  → floor(vmin) (base at/below the data),
 *    - straddling 0    → null (signed values are already short enough).
 *  For a run that does NOT cross an integer this is exactly the integer part of
 *  the run (ceil(vmin) == ceil(vmax) for negatives, floor likewise for
 *  positives), so the familiar '.67 … .64' tick style is preserved.
 *  [2026-09-13 | Atreides-Jimmy] Changed from "integer part of vmin" to this
 *  side-picking base. For a run that did cross an integer (e.g. -6911.977 …
 *  -6912.806) the old base (-6912) put ticks on both sides and the absolute
 *  fraction folded into a non-monotonic '0 0 0 1 1'. */
function convAxisBase(vmin,vmax){
    if(Math.abs(vmax-vmin)>=1)return null;
    if(vmin<0&&vmax>0)return null;
    return vmax<=0?Math.ceil(vmax):Math.floor(vmin);
}

/** Widen the left gutter so the widest of the five tick labels fits. The labels
 *  are no longer assumed to be short — the decimal count now follows the data,
 *  so a very fine span can want more digits than the old fixed 40px gutter held.
 *  Returns padL between minPad and 45% of the canvas width.
 *  [2026-09-13 | Atreides-Jimmy] Replaces the old pass/fail gutter check, which
 *  reacted to an overflow by dropping the base and printing even wider values. */
function convFitGutter(ctx,labelFor,vmin,vmax,minPad,w){
    ctx.font='9px sans-serif';
    var widest=0;
    for(var g=0;g<=4;g++){
        var t=ctx.measureText(labelFor(vmax-(vmax-vmin)*g/4)).width;
        if(t>widest)widest=t;
    }
    var pad=Math.ceil(widest)+6;
    if(pad<minPad)pad=minPad;
    var cap=Math.floor(w*0.45);
    if(pad>cap)pad=cap;
    return pad;
}

/** [min,max] of one numeric field over the plotted steps (null when empty). */
function convFieldRange(steps,field){
    var lo=Infinity,hi=-Infinity,any=false;
    for(var i=0;i<steps.length;i++){
        var v=steps[i][field];
        if(v==null||isNaN(v)||!isFinite(v))continue;
        if(v<lo)lo=v;
        if(v>hi)hi=v;
        any=true;
    }
    return any?[lo,hi]:null;
}

/** Decimal places that resolve one plotted series over the steps (its own span,
 *  via convStepDecimals), used by the hover tooltip. A total energy around
 *  -6912 with a sub-unit spread used to fall back to 0 decimals and print as a
 *  bare integer; this keeps the decimals that actually separate the steps.
 *  [2026-09-13 | Atreides-Jimmy] */
function convSeriesDecimals(steps,field){
    var r=convFieldRange(steps,field);
    return r?convStepDecimals(r[1]-r[0]):6;
}

/** Tooltip number: enough decimals to (a) resolve the value's own magnitude,
 *  (b) resolve the series' span so successive steps stay distinguishable, and
 *  (c) never drop below millihartree-level detail (3 decimals). Capped at 10 so
 *  a near-degenerate span cannot produce an absurdly wide tooltip.
 *  [2026-09-13 | Atreides-Jimmy] */
function convTooltipValue(v,spanDec){
    var d=convValueDecimals(v);
    if(spanDec!=null&&spanDec>d)d=spanDec;
    if(d<3)d=3;
    if(d>10)d=10;
    return convFormat(v,d);
}

/** Index of the optimization step whose structure is currently displayed, or -1.
 *  Frames of an opt log are the per-step geometries in step order (the parser
 *  numbers them "Step k"), so the frame index IS the step index; the frame's
 *  own single-point energy is only a guard for formats/labels that do not
 *  follow that numbering. */
function convSelectedStep(steps){
    if(typeof currentFrame!=='number'||currentFrame<0||currentFrame>=steps.length)return -1;
    var f=MD.frames&&MD.frames[currentFrame];
    var e=f?f.energy:undefined;
    if(e!=null&&!isNaN(e)&&steps[currentFrame].energy!=null&&Math.abs(steps[currentFrame].energy-e)>1e-9)return -1;
    return currentFrame;
}

/** Draw the value tooltip of the data point nearest to the cursor (feature 3)
 *  and report whether a point was hit.
 *
 *  The tooltip is a pure function of the chart (canvas.__convGeom) plus the last
 *  pointer position (canvas.__convHover), so a repaint either draws one tooltip
 *  or none — it can never leave an older tooltip behind on the canvas. Points of
 *  a different step are ignored so the tooltip never reports a value the user is
 *  not pointing at.
 *  [2026-09-12 | Atreides-Jimmy] Changed: returns whether the tooltip was drawn,
 *  so the pointer handler can tell "still hovering a point" from "just left it"
 *  and repaint on that transition (it used to be painted from mousemove only
 *  while a point was in range, which left the last tooltip on screen until the
 *  pointer left the canvas). */
function convDrawTooltip(ctx,canvas){
    var geom=canvas.__convGeom;
    if(!geom||!canvas.__convHover||!geom.series.length)return false;
    var mx=canvas.__convHover.x,my=canvas.__convHover.y;
    var best=null,bestD=289; // 17² px picking radius
    for(var s=0;s<geom.series.length;s++){
        var pts=geom.series[s].pts;
        for(var p=0;p<pts.length;p++){
            var dx=pts[p].x-mx,dy=pts[p].y-my,d=dx*dx+dy*dy;
            if(d<=bestD){bestD=d;best={i:pts[p].i,x:pts[p].x,y:pts[p].y,v:pts[p].v,color:geom.series[s].color,name:geom.series[s].label,dec:geom.series[s].dec}}
        }
    }
    if(!best)return false;
    var step=steps_get(geom.steps,best.i);
    var lines=[(geom.xLabel||'Step')+' '+(best.i+1)];
    lines.push(best.name+' = '+convTooltipValue(best.v,best.dec));
    if(step&&step.energy!=null&&geom.field!=='energy')lines.push('Energy = '+convTooltipValue(step.energy,geom.energyDec));
    ctx.font='9px sans-serif';
    var wid=0;
    for(var li=0;li<lines.length;li++)wid=Math.max(wid,ctx.measureText(lines[li]).width);
    var bw=wid+10,bh=lines.length*11+7;
    var bx=best.x+10,by=best.y-bh-7;
    if(bx+bw>geom.padL+geom.pw)bx=best.x-bw-10;
    if(bx<1)bx=1;
    if(by<1)by=best.y+10;
    ctx.fillStyle='rgba(12,14,18,0.94)';
    ctx.strokeStyle=best.color;ctx.lineWidth=1;
    ctx.beginPath();
    if(ctx.roundRect)ctx.roundRect(bx,by,bw,bh,3);else ctx.rect(bx,by,bw,bh);
    ctx.fill();ctx.stroke();
    ctx.textAlign='left';ctx.textBaseline='top';
    var ty=by+3;
    for(li=0;li<lines.length;li++){
        ctx.fillStyle=li===lines.length-1&&lines.length>1?'#ccc':best.color;
        ctx.fillText(lines[li],bx+5,ty);
        ty+=11;
    }
    ctx.strokeStyle=best.color;ctx.lineWidth=1;
    ctx.beginPath();ctx.arc(best.x,best.y,3,0,Math.PI*2);ctx.stroke();
    return true;
}

/** Step object at idx, tolerating a missing steps array. */
function steps_get(steps,i){return steps&&steps[i]?steps[i]:null}

/** Repaint a chart with the current pointer state: the selected-step highlight,
 *  plus a tooltip when the pointer is close enough to a point. Every repaint
 *  clears the canvas first (the draw functions do), so the tooltip disappears in
 *  the same frame the pointer leaves its point.
 *  [2026-09-12 | Atreides-Jimmy] Added: single repaint entry point for hover. */
function convApplyHover(canvas,pt){
    canvas.__convHover=pt||null;
    if(typeof canvas.__convRedraw==='function')canvas.__convRedraw();
}

/** Wire the hover interaction once per canvas (feature 3). The repaint is
 *  throttled to real sub-pixel motion; as soon as the pointer moves off the
 *  point it was on, that repaint drops the tooltip (previously the tooltip was
 *  only redrawn when another point was reached, so fast sweeps left several
 *  tooltip boxes stacked on the canvas until the pointer left the chart). */
function convAttachHover(canvas){
    if(canvas.__convHoverBound)return;
    canvas.__convHoverBound=true;
    canvas.addEventListener('mousemove',function(e){
        var r=canvas.getBoundingClientRect();
        var pt={x:e.clientX-r.left,y:e.clientY-r.top};
        var prev=canvas.__convHover;
        if(prev&&Math.abs(prev.x-pt.x)<1&&Math.abs(prev.y-pt.y)<1)return;
        convApplyHover(canvas,pt);
    });
    canvas.addEventListener('mouseleave',function(){
        if(canvas.__convHover)convApplyHover(canvas,null);
    });
}

/** Highlight the point of the displayed step (feature 2), then draw the value
 *  tooltip of whatever point the pointer is currently on (feature 3). Called at
 *  the end of every chart repaint. */
function convDrawState(ctx,canvas){
    var geom=canvas.__convGeom;
    if(!geom)return;
    var sel=geom.selected;
    if(sel!=null&&sel>=0){
        var sx=geom.padL+(geom.n<=1?0:geom.pw*sel/(geom.n-1));
        ctx.save();
        ctx.strokeStyle='rgba(255,214,102,0.55)';ctx.lineWidth=1;ctx.setLineDash([2,2]);
        ctx.beginPath();ctx.moveTo(sx,geom.padT);ctx.lineTo(sx,geom.padT+geom.ph);ctx.stroke();
        ctx.setLineDash([]);
        for(var s=0;s<geom.series.length;s++){
            var pts=geom.series[s].pts;
            for(var p=0;p<pts.length;p++){
                if(pts[p].i!==sel)continue;
                ctx.beginPath();ctx.arc(pts[p].x,pts[p].y,5.5,0,Math.PI*2);
                ctx.fillStyle='rgba(255,214,102,0.22)';ctx.fill();
                ctx.beginPath();ctx.arc(pts[p].x,pts[p].y,3.4,0,Math.PI*2);
                ctx.fillStyle=pts[p].color;ctx.fill();
                ctx.strokeStyle='#ffd666';ctx.lineWidth=1.4;ctx.stroke();
            }
        }
        ctx.restore();
    }
    convDrawTooltip(ctx,canvas);
}

function drawConvergenceChart(canvas,steps,field,label,color,threshold){
    if(!canvas.__convRedraw){
        canvas.__convRedraw=function(){drawConvergenceChart(canvas,steps,field,label,color,threshold)};
        convAttachHover(canvas);
    }
    var ctx=canvas.getContext('2d');
    var dpr=window.devicePixelRatio||1;
    var w=(canvas.parentElement&&canvas.parentElement.clientWidth)||canvas.clientWidth||300;
    var h=canvas.clientHeight||100;
    canvas.width=w*dpr;canvas.height=h*dpr;
    canvas.style.width=w+'px';canvas.style.height=h+'px';canvas.style.maxWidth='100%';canvas.style.display='block';
    ctx.setTransform(dpr,0,0,dpr,0,0);
    ctx.clearRect(0,0,w,h);
    var padL=40,padR=8,padT=14,padB=18;
    var pw=w-padL-padR,ph=h-padT-padB;
    var range=convFieldRange(steps,field);
    if(!range){ctx.fillStyle='#888';ctx.font='10px sans-serif';ctx.fillText('No data',padL,padT+12);return}
    var dmin=range[0],dmax=range[1];
    var vmin=dmin,vmax=dmax;
    if(threshold!=null){vmin=Math.min(vmin,threshold);vmax=Math.max(vmax,threshold)}
    if(vmax-vmin<1e-12){vmax=vmin+Math.abs(vmin)*1e-6+1e-9}
    var n=steps.length;
    var sel=steps.length?convSelectedStep(steps):-1;

    // --- Feature 1: decimal axis labels, integer part shown once -----------
    // Decimals come from the gridline step (span / 4) and the base is picked on
    // ONE side of the plotted range, so every gridline gets a distinct,
    // monotonic label even when the run crosses an integer boundary.
    var dec=convTickDecimals(vmin,vmax);
    var vbase=convAxisBase(vmin,vmax);
    function labelFor(v){return convTickLabel(v,vbase,dec)}
    // Widen the gutter to hold the labels; the plot width then follows.
    padL=convFitGutter(ctx,labelFor,vmin,vmax,40,w);
    pw=w-padL-padR;
    // Chart title: series label + the shared integer part of the run (vbase is
    // exactly that integer, so the two can never disagree).
    ctx.fillStyle='#ccc';ctx.font='10px sans-serif';ctx.textAlign='left';ctx.textBaseline='top';
    ctx.fillText(vbase==null?label:label+'  '+vbase+' +',padL,0);

    // Grid + axis labels
    ctx.strokeStyle='rgba(255,255,255,0.12)';ctx.lineWidth=1;
    ctx.fillStyle='#888';ctx.font='9px sans-serif';ctx.textAlign='right';ctx.textBaseline='middle';
    for(var g=0;g<=4;g++){
        var y=padT+ph*g/4;
        ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(padL+pw,y);ctx.stroke();
        ctx.fillText(labelFor(vmax-(vmax-vmin)*g/4),padL-3,y);
    }
    ctx.textAlign='center';ctx.textBaseline='top';
    for(g=0;g<n;g+=Math.max(1,Math.floor(n/6))){
        var x=padL+(n<=1?0:pw*g/(n-1));
        ctx.fillText(String(g+1),x,padT+ph+2);
    }
    // Threshold line
    if(threshold!=null){
        var ty=padT+ph*(vmax-threshold)/(vmax-vmin);
        if(padL-3-ctx.measureText(labelFor(threshold)).width>=0){
            ctx.fillStyle='#888';ctx.textAlign='right';ctx.textBaseline='middle';
            ctx.fillText(labelFor(threshold),padL-3,ty);
        }
        ctx.strokeStyle='rgba(255,180,0,0.6)';ctx.setLineDash([3,3]);
        ctx.beginPath();ctx.moveTo(padL,ty);ctx.lineTo(padL+pw,ty);ctx.stroke();
        ctx.setLineDash([]);
    }
    // Data line + points
    var pts=[];
    ctx.strokeStyle=color;ctx.lineWidth=1.5;ctx.beginPath();
    var drew=false;
    for(var i=0;i<n;i++){
        var v=steps[i][field];
        if(v==null||isNaN(v))continue;
        x=padL+(n<=1?0:pw*i/(n-1));
        y=padT+ph*(vmax-v)/(vmax-vmin);
        pts.push({i:i,x:x,y:y,v:v,color:color});
        if(!drew){ctx.moveTo(x,y);drew=true}else ctx.lineTo(x,y);
    }
    ctx.stroke();
    ctx.fillStyle=color;
    for(i=0;i<pts.length;i++){
        ctx.beginPath();ctx.arc(pts[i].x,pts[i].y,2,0,Math.PI*2);ctx.fill();
    }
    canvas.__convGeom={padL:padL,padT:padT,pw:pw,ph:ph,n:n,xLabel:'Step',yLabel:label,
        spanVal:vmax-vmin,field:field,steps:steps,
        series:[{label:label,color:color,pts:pts,dec:convSeriesDecimals(steps,field)}],selected:sel};
    convDrawState(ctx,canvas);
}

/** Redraw every convergence chart of the opt panel. Called after the panel is
 *  built and whenever the displayed structure changes, so the selected-step
 *  highlight (feature 2) always follows the 3D view.
 *  [2026-09-12 | Atreides-Jimmy] */
function renderOptCharts(){
    if(!OPT_STEPS||OPT_STEPS.length===0)return;
    // Nothing to redraw while the panel is closed (its canvases have no layout
    // then); buildOptPanel draws them again once the panel is shown.
    if(!optPanelEl.classList.contains('show'))return;
    var canvases=optPanelEl.querySelectorAll('canvas');
    canvases.forEach(function(cv){
        var field=cv.dataset.field;
        if(field==='energy'){
            cv.style.height='90px';
            drawConvergenceChart(cv,OPT_STEPS,'energy','Energy (Hartree)','#3794ff',null);
        }else if(field==='force'){
            cv.style.height='110px';
            drawDualChart(cv,OPT_STEPS,['maxForce','rmsForce'],['Max Force','RMS Force'],['#ff6b6b','#ffa94d'],OPT_THRESHOLDS?OPT_THRESHOLDS.maxForce:null,OPT_THRESHOLDS?OPT_THRESHOLDS.rmsForce:null);
        }else if(field==='disp'){
            cv.style.height='110px';
            drawDualChart(cv,OPT_STEPS,['maxDisplacement','rmsDisplacement'],['Max Disp','RMS Disp'],['#51cf66','#74c0fc'],OPT_THRESHOLDS?OPT_THRESHOLDS.maxDisplacement:null,OPT_THRESHOLDS?OPT_THRESHOLDS.rmsDisplacement:null);
        }
    });
}

function buildOptPanel(){
    if(!OPT_STEPS||OPT_STEPS.length===0)return;
    var hasEnergy=OPT_STEPS.some(function(s){return s.energy!=null});
    var hasForce=OPT_STEPS.some(function(s){return s.maxForce!=null||s.rmsForce!=null});
    var hasDisp=OPT_STEPS.some(function(s){return s.maxDisplacement!=null||s.rmsDisplacement!=null});
    var chartsHTML='';
    if(hasEnergy)chartsHTML+='<div class="opt-chart-title">Energy (Hartree)</div><canvas data-field="energy" data-color="#3794ff" data-label="Energy"></canvas>';
    if(hasForce)chartsHTML+='<div class="opt-chart-title">Force (max + rms)</div><canvas data-field="force" data-color="#ff6b6b" data-label="Force"></canvas>';
    if(hasDisp)chartsHTML+='<div class="opt-chart-title">Displacement (max + rms)</div><canvas data-field="disp" data-color="#51cf66" data-label="Displacement"></canvas>';
    optPanelEl.innerHTML=
        '<div class="opt-head"><span class="opt-title">Optimization Convergence ('+OPT_STEPS.length+' steps)</span>'+
        '<span class="opt-close">×</span></div>'+
        '<div class="opt-resize" id="opt-resize"></div>'+
        '<div class="opt-body" id="opt-body">'+chartsHTML+'</div>';
    optPanelEl.classList.add('show');
    optReopenEl.classList.remove('show');

    // Draw charts (deferred to ensure layout is computed)
    requestAnimationFrame(renderOptCharts);

    // Close handler
    var closeBtn=optPanelEl.querySelector('.opt-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            optPanelEl.classList.remove('show');
            optReopenEl.classList.add('show');
            layoutPanels();
        });
    }
    // Resize handler (vertical drag)
    var resizeH=optPanelEl.querySelector('#opt-resize');
    if(resizeH){
        var dragging=false,startY=0,startH=0;
        resizeH.addEventListener('mousedown',function(e){dragging=true;startY=e.clientY;startH=optPanelEl.offsetHeight;e.preventDefault();document.body.style.cursor='ns-resize'});
        document.addEventListener('mousemove',function(e){
            if(!dragging)return;
            var nh=Math.max(120,Math.min(600,startH-(e.clientY-startY)));
            optPanelEl.style.maxHeight=nh+'px';
            // Redraw charts
            renderOptCharts();
        });
        document.addEventListener('mouseup',function(){if(dragging){dragging=false;document.body.style.cursor='';layoutPanels()}});
    }
    layoutPanels();
}

function drawDualChart(canvas,steps,fields,labels,colors,thresh1,thresh2){
    if(!canvas.__convRedraw){
        canvas.__convRedraw=function(){drawDualChart(canvas,steps,fields,labels,colors,thresh1,thresh2)};
        convAttachHover(canvas);
    }
    var ctx=canvas.getContext('2d');
    var dpr=window.devicePixelRatio||1;
    var w=(canvas.parentElement&&canvas.parentElement.clientWidth)||canvas.clientWidth||300;
    var h=canvas.clientHeight||110;
    canvas.width=w*dpr;canvas.height=h*dpr;
    canvas.style.width=w+'px';canvas.style.height=h+'px';canvas.style.maxWidth='100%';canvas.style.display='block';
    ctx.setTransform(dpr,0,0,dpr,0,0);
    ctx.clearRect(0,0,w,h);
    var padL=40,padR=8,padT=14,padB=18;
    var pw=w-padL-padR,ph=h-padT-padB;
    var allVals=[];
    steps.forEach(function(s){fields.forEach(function(f){if(s[f]!=null&&!isNaN(s[f]))allVals.push(s[f])})});
    if(allVals.length===0){ctx.fillStyle='#888';ctx.font='10px sans-serif';ctx.fillText('No data',padL,padT+12);return}
    // Logarithmic Y axis: force/displacement values span decades across an
    // optimization, so a linear axis squeezes the convergence-threshold lines
    // against the bottom once early steps dominate the range. Non-positive
    // values (e.g. a converged step printed as 0.000000) are clamped to a
    // floor just below the smallest positive value; falls back to linear
    // only when nothing positive remains or a threshold is non-positive.
    var posVals=allVals.filter(function(v){return v>0});
    var useLog=posVals.length>0;
    [thresh1,thresh2].forEach(function(th){if(th!=null&&th<=0)useLog=false});
    var vmin,vmax,vFloor=0;
    if(useLog){
        vFloor=Math.min.apply(null,posVals)*1e-2;
        var lv=allVals.map(function(v){return Math.log10(Math.max(v,vFloor))});
        if(thresh1!=null)lv.push(Math.log10(thresh1));
        if(thresh2!=null)lv.push(Math.log10(thresh2));
        vmin=Math.min.apply(null,lv);vmax=Math.max.apply(null,lv);
        if(vmax-vmin<1e-12)vmax=vmin+1;
    }else{
        vmin=Math.min.apply(null,allVals),vmax=Math.max.apply(null,allVals);
        if(thresh1!=null){vmin=Math.min(vmin,thresh1);vmax=Math.max(vmax,thresh1)}
        if(thresh2!=null){vmin=Math.min(vmin,thresh2);vmax=Math.max(vmax,thresh2)}
        if(vmax-vmin<1e-12)vmax=vmin+1;
    }
    function yOf(v){
        if(useLog)return padT+ph*(vmax-Math.log10(Math.max(v,vFloor)))/(vmax-vmin);
        return padT+ph*(vmax-v)/(vmax-vmin);
    }
    function valAtFrac(frac){
        // Axis value at vertical fraction frac (0=top): inverse of yOf.
        if(useLog)return Math.pow(10,vmax-(vmax-vmin)*frac);
        return vmax-(vmax-vmin)*frac;
    }
    var n=steps.length;
    var sel=convSelectedStep(steps);
    // [2026-09-12 | Atreides-Jimmy] The linear fallback now labels its ticks in
    // the same condensed decimal notation as the energy chart (feature 1); the
    // log axis keeps its log10 tick values by design.
    // [2026-09-13 | Atreides-Jimmy] The linear labels use the true linear axis
    // extremes (which may include a threshold), so base/decimals match the ticks
    // that are actually drawn; the gutter is widened to fit them.
    var dmin=Math.min.apply(null,allVals),dmax=Math.max.apply(null,allVals);
    var decLin=useLog?0:convTickDecimals(vmin,vmax);
    var vbase=useLog?null:convAxisBase(vmin,vmax);
    function linLabel(v){
        return vbase!=null?convTickLabel(v,vbase,decLin):v.toFixed(decLin);
    }
    if(!useLog){
        padL=convFitGutter(ctx,linLabel,vmin,vmax,40,w);
        pw=w-padL-padR;
    }
    // Grid + axis label. On the log axis the tick values ARE the log10
    // values (evenly spaced, e.g. -3.35), with a rotated "log₁₀" axis label
    // on the left.
    ctx.strokeStyle='rgba(255,255,255,0.12)';ctx.lineWidth=1;
    ctx.fillStyle='#888';ctx.font='9px sans-serif';ctx.textAlign='right';ctx.textBaseline='middle';
    for(var g=0;g<=4;g++){
        var y=padT+ph*g/4;
        ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(padL+pw,y);ctx.stroke();
        if(useLog)ctx.fillText((vmax-(vmax-vmin)*g/4).toFixed(2),padL-3,y);
        else ctx.fillText(linLabel(valAtFrac(g/4)),padL-3,y);
    }
    if(useLog){
        ctx.save();
        ctx.translate(9,padT+ph/2);
        ctx.rotate(-Math.PI/2);
        ctx.textAlign='center';ctx.textBaseline='middle';
        ctx.fillText('log\u2081\u2080',0,0);
        ctx.restore();
        ctx.fillStyle='#888';ctx.font='9px sans-serif';ctx.textAlign='right';ctx.textBaseline='middle';
    }
    ctx.textAlign='center';ctx.textBaseline='top';
    for(g=0;g<n;g+=Math.max(1,Math.floor(n/6))){
        var x=padL+(n<=1?0:pw*g/(n-1));
        ctx.fillText(String(g+1),x,padT+ph+2);
    }
    // Threshold lines
    [thresh1,thresh2].forEach(function(th,idx){
        if(th==null)return;
        var ty=yOf(th);
        ctx.strokeStyle='rgba(255,180,0,0.5)';ctx.setLineDash([3,3]);
        ctx.beginPath();ctx.moveTo(padL,ty);ctx.lineTo(padL+pw,ty);ctx.stroke();
        ctx.setLineDash([]);
    });
    // Legend
    ctx.font='9px sans-serif';ctx.textAlign='left';ctx.textBaseline='top';
    for(var fi=0;fi<fields.length;fi++){
        ctx.fillStyle=colors[fi];
        ctx.fillRect(padL+fi*80,padT-12,8,2);
        ctx.fillText(labels[fi],padL+fi*80+10,padT-14);
    }
    // Data lines + points (kept for the hover tooltip / selected-step mark)
    var series=[];
    for(fi=0;fi<fields.length;fi++){
        ctx.strokeStyle=colors[fi];ctx.lineWidth=1.5;ctx.beginPath();
        var drew=false;
        var pts=[];
        for(var i=0;i<n;i++){
            var v=steps[i][fields[fi]];
            if(v==null||isNaN(v))continue;
            x=padL+(n<=1?0:pw*i/(n-1));
            y=yOf(v);
            pts.push({i:i,x:x,y:y,v:v,color:colors[fi]});
            if(!drew){ctx.moveTo(x,y);drew=true}else ctx.lineTo(x,y);
        }
        ctx.stroke();
        ctx.fillStyle=colors[fi];
        for(i=0;i<pts.length;i++){
            ctx.beginPath();ctx.arc(pts[i].x,pts[i].y,2,0,Math.PI*2);ctx.fill();
        }
        series.push({label:labels[fi],color:colors[fi],pts:pts,dec:convSeriesDecimals(steps,fields[fi])});
    }
    // spanVal states how fine the data (not the log axis) is: it drives the
    // tooltip's decimal resolution, so it must come from the plotted values.
    canvas.__convGeom={padL:padL,padT:padT,pw:pw,ph:ph,n:n,xLabel:'Step',yLabel:labels[0],
        spanVal:Math.max.apply(null,allVals)-Math.min.apply(null,allVals),field:fields[0],steps:steps,
        series:series,selected:sel,energyDec:convSeriesDecimals(steps,'energy')};
    convDrawState(ctx,canvas);
}

// Vibration animation
var vibActive=false,vibModeIdx=-1,vibBasePositions=null,vibTimer=null;
var vibPaused=false,vibPausedElapsed=0,vibStartTime=0;
var vibPeriodMs=3000,vibAmplitude=0.35,vibMaxMag=1;
// A playback frame moves every atom and then re-places every bond piece from
// the new coordinates — the same routine the scene build and a drag use, so a
// bond can never disagree with the atoms it spans. Only the instance buffers
// are uploaded; no geometry is rebuilt and no mesh is created per frame.
function updateVibBondMeshes(){
    placeBondPieces();
}
function updateVibrationPositions(){
    for(var i=0;i<MD.atoms.length&&i<atomMeshes.length;i++){
        atomMeshes[i].position.set(MD.atoms[i].x-CX,MD.atoms[i].y-CY,MD.atoms[i].z-CZ);
    }
    updateVibBondMeshes();
}
function startVibration(modeIdx){
    stopVibration();
    if(!NORMAL_MODES||modeIdx<0||modeIdx>=NORMAL_MODES.length)return;
    // Playing a mode is the animation half of a mode's selection and nothing
    // else: the manual-displacement target (vdTarget) is set by its own circle
    // in the list, so watching a mode can never re-aim the slider.
    // A normal mode is a playback, not the user's geometry: it rewrites every
    // atom's position on each frame and puts the geometry it started from back
    // when it stops. A Molecular Mechanics run writing coordinates into the same
    // array would fight it frame by frame — and stopping the mode would then
    // rewind the run's progress. So the run is ended first, exactly as any other
    // change of what the structure is: the geometry it reached so far stays on
    // screen and the mode plays over that.
    var hadRun=mmRunActive();
    mmCancelRun();
    if(hadRun)modeInfoEl.textContent='Molecular Mechanics: the run was ended so the normal mode can play';
    // The two other loops that write the structure end here as well. A frame
    // auto-play replaces the structure on every tick — and since a frame change
    // stops a playback, leaving it on would kill this one a moment after it
    // started — and a Move Atoms session is mid-edit on the selection the
    // playback is about to rewrite.
    stopFrameAutoPlay();
    stopMoveSession();
    if(totalFrames>0&&currentFrame!==totalFrames-1){switchFrame(totalFrames-1)}
    var mode=NORMAL_MODES[modeIdx];
    if(!mode.displacements||mode.displacements.length!==MD.atoms.length)return;
    // A manual displacement is the user's structure, not a pose: the playback
    // always swings from the *undisplaced* equilibrium, so a mode looks the same
    // whatever number the slider was left at, and stopping puts the displaced
    // geometry back (see stopVibration). Without a session the geometry on
    // screen already is the equilibrium.
    vibBasePositions=(vdActive&&vdBase&&vdBase.length===MD.atoms.length)
        ?vdBase.map(function(p){return{x:p.x,y:p.y,z:p.z}})
        :MD.atoms.map(function(a){return{x:a.x,y:a.y,z:a.z}});
    var maxMag=0;
    mode.displacements.forEach(function(d){var m=Math.sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);if(m>maxMag)maxMag=m});
    if(maxMag<1e-8)maxMag=1;
    vibAmplitude=0.35;
    vibMaxMag=maxMag;
    vibActive=true;vibPaused=false;vibModeIdx=modeIdx;
    vibPausedElapsed=0;
    vibStartTime=performance.now();
    vibTimer=requestAnimationFrame(vibTick);
    updateFreqButton('pause');
}
function vibTick(){
    if(!vibActive)return;
    var mode=NORMAL_MODES[vibModeIdx];
    if(!mode||!vibBasePositions)return;
    var elapsed=vibPausedElapsed+(performance.now()-vibStartTime);
    var s=Math.sin(2*Math.PI*(elapsed%vibPeriodMs)/vibPeriodMs)*vibAmplitude/vibMaxMag;
    for(var i=0;i<MD.atoms.length&&i<mode.displacements.length;i++){
        var d=mode.displacements[i];
        MD.atoms[i].x=vibBasePositions[i].x+d[0]*s;
        MD.atoms[i].y=vibBasePositions[i].y+d[1]*s;
        MD.atoms[i].z=vibBasePositions[i].z+d[2]*s;
    }
    updateVibrationPositions();
    vibTimer=requestAnimationFrame(vibTick);
}
function pauseVibration(){
    if(!vibActive)return;
    vibPausedElapsed+=performance.now()-vibStartTime;
    if(vibTimer){cancelAnimationFrame(vibTimer);vibTimer=null}
    vibActive=false;vibPaused=true;
    updateFreqButton('resume');
}
function resumeVibration(){
    if(!vibPaused)return;
    vibActive=true;vibPaused=false;
    vibStartTime=performance.now();
    vibTimer=requestAnimationFrame(vibTick);
    updateFreqButton('pause');
}
function updateFreqButton(state){
    var btn=freqPanelEl.querySelector('#freq-stop');
    if(!btn)return;
    btn.textContent=state==='resume'?'▶ Resume':'⏸ Pause';
}
function stopVibration(){
    if(vibTimer){cancelAnimationFrame(vibTimer);vibTimer=null}
    if(vibBasePositions){
        for(var i=0;i<MD.atoms.length&&i<vibBasePositions.length;i++){
            MD.atoms[i].x=vibBasePositions[i].x;
            MD.atoms[i].y=vibBasePositions[i].y;
            MD.atoms[i].z=vibBasePositions[i].z;
        }
        vibBasePositions=null;
        // A manual displacement is the geometry the user left the structure in;
        // the playback only borrowed it, so it is written back before the
        // animation disappears. Without a session this is a no-op.
        vdWriteAtoms();
        updateVibrationPositions();
    }
    vibActive=false;vibPaused=false;vibModeIdx=-1;
    vibPausedElapsed=0;
    // Update freq panel highlight
    var rows=freqPanelEl.querySelectorAll('.freq-row');
    rows.forEach(function(r){r.classList.remove('playing')});
    updateFreqButton('pause');
}

// ===== Manual displacement along normal modes =====
// The playback above is a *view* of a mode: it borrows the geometry for one
// sine cycle and puts it back. Manual displacement is the other half — modes
// written into the coordinates, so the structure can be pulled along them and
// then measured, rotated, edited or saved like any other.
//
// The distortion is a *combination* of modes, not one at a time. The normal
// modes of one frequency job are an orthogonal basis of the 3N-6 displacements
// (orthogonal under the mass-weighted inner product, which the printed vectors
// already carry), so
//     x = x0 + sum_k (scale_k * dir_k)
// is exactly a point in normal-coordinate space: displacing along a second mode
// is the same arithmetic as a scan along a line in that space, and every
// component keeps the per-atom weights and phases the calculation produced. A
// component is added by dragging the slider of a picked mode and removed by
// dragging it back to 0, so nothing accumulates invisibly — the block lists
// every component that is in the structure, and every contributing row is
// marked.
//
// State:
//   * vdBase — the equilibrium geometry the whole combination is measured
//     from, captured once when the first component is added. Every write is
//     rebuilt from it (base + sum), never from the previous frame, so dragging
//     back and forth cannot accumulate drift.
//   * vdComp — the components, as {idx, scale, dir}, where dir is that mode's
//     unit displacement (largest per-atom magnitude 1 Å).
//   * vdTarget — the mode the slider points at (the one whose circle is
//     filled). The slider edits *that* component; it plays nothing.
//
// The slider is a *fraction* of one mode, not a distance: -1.00 to +1.00, where
// 1.00 moves that mode's most-moving atom by 1 Å. The block prints what each
// component comes to and the largest displacement of the whole combination
// (which is not the sum of the parts — vectors can cancel), so a number is
// never the only thing on offer.
//
// A session is a *pending* rewrite of every atom. Anything that edits the
// structure closes it first (vdCommit, called from pushUndo, which every
// coordinate edit opens with), and anything that replaces the structure drops
// it (vdEnd). The playback is the one deliberate exception: it swings around
// the *pure* vdBase — the mode it shows is the mode the calculation produced,
// not that mode on top of an arbitrary offset — and writes the combination back
// when it stops.
var VD_MAX_SCALE=1.0;    // slider range: -1 … +1 of one mode's displacement
var vdActive=false;      // a displacement session holds the geometry
var vdTarget=-1;         // mode the slider points at (its circle is filled)
var vdComp=[];           // components: {idx, scale, dir}
var vdBase=null;         // equilibrium coordinates [{x,y,z}, ...]
var vdMaxDisp=0;         // largest atomic displacement of the combination, Å
var vdUndoDepth=-1;      // undoStack length right after the session's snapshot
var vdSuppress=false;    // the session's own pushUndo must not close it
var vdMarkedSig='';      // component list the row markers were last built for
var vdPickedRow=-1;      // row whose pick circle is currently filled

function vdAvailable(){
    return !CRY&&!diffMode&&!!NORMAL_MODES&&NORMAL_MODES.length>0&&MD.atoms.length>0;
}
/** Whether one mode can be written into the current structure at all. */
function vdUsable(idx){
    if(!NORMAL_MODES||idx<0||idx>=NORMAL_MODES.length)return false;
    var m=NORMAL_MODES[idx];
    return !!(m&&m.displacements&&m.displacements.length===MD.atoms.length);
}
/** Unit displacement vectors of one mode, scaled so the largest per-atom
 *  magnitude is 1 Å — that is what makes the slider's 1.00 mean "the
 *  most-moving atom travels 1 Å". Null when the mode carries no usable vector
 *  for this structure. */
function vdDirFor(idx){
    if(!vdUsable(idx))return null;
    var d=NORMAL_MODES[idx].displacements,n=MD.atoms.length;
    var maxSq=0,i,v,q;
    for(i=0;i<n;i++){
        v=d[i];
        q=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
        if(q>maxSq)maxSq=q;
    }
    if(!(maxSq>1e-24))return null;
    var inv=1/Math.sqrt(maxSq);
    var dir=new Array(n);
    for(i=0;i<n;i++){
        v=d[i];
        dir[i]=[v[0]*inv,v[1]*inv,v[2]*inv];
    }
    return dir;
}
/** The component of one mode (null when that mode is not in the distortion). */
function vdFindComp(idx){
    for(var i=0;i<vdComp.length;i++)if(vdComp[i].idx===idx)return vdComp[i];
    return null;
}
/** The share of the picked mode: what the slider shows and edits. */
function vdCurScale(){
    var c=vdFindComp(vdTarget);
    return c?c.scale:0;
}
/** Rebuild every atom from the equilibrium plus every component. Returns false
 *  (and writes nothing) when the session and the structure no longer describe
 *  the same atom list. */
function vdWriteAtoms(){
    var n=MD.atoms.length,i,c;
    vdMaxDisp=0;
    if(!vdActive||!vdBase||vdBase.length!==n)return false;
    // Every component must describe this atom list before anything is written:
    // a half-applied combination would be a structure nobody asked for.
    for(c=0;c<vdComp.length;c++){
        var cd=vdComp[c].dir;
        if(!cd||cd.length!==n)return false;
    }
    for(i=0;i<n;i++){
        var a=MD.atoms[i],b=vdBase[i];
        a.x=b.x;a.y=b.y;a.z=b.z;
    }
    for(c=0;c<vdComp.length;c++){
        var comp=vdComp[c],dir=comp.dir,sc=comp.scale;
        for(i=0;i<n;i++){
            var a2=MD.atoms[i],v=dir[i];
            a2.x+=v[0]*sc;a2.y+=v[1]*sc;a2.z+=v[2]*sc;
        }
    }
    // The largest displacement of the *combination* is not the sum of the
    // components' own (vectors can cancel), so it is measured on the result.
    var maxSq=0;
    for(i=0;i<n;i++){
        var a3=MD.atoms[i],b3=vdBase[i];
        var dx=a3.x-b3.x,dy=a3.y-b3.y,dz=a3.z-b3.z;
        var q=dx*dx+dy*dy+dz*dz;
        if(q>maxSq)maxSq=q;
    }
    vdMaxDisp=Math.sqrt(maxSq);
    return true;
}
/** Add a mode to the distortion or change its share; false when the mode
 *  carries no usable displacement for this structure. */
function vdSetComp(idx,scale){
    var c=vdFindComp(idx);
    if(c){c.scale=scale;return true}
    var dir=vdDirFor(idx);
    if(!dir)return false;
    vdComp.push({idx:idx,scale:scale,dir:dir});
    return true;
}
/** Take one mode out of the distortion. */
function vdDropComp(idx){
    for(var i=0;i<vdComp.length;i++){
        if(vdComp[i].idx===idx){vdComp.splice(i,1);return true}
    }
    return false;
}
/** The combination as the block states it, so the structure is always
 *  traceable to the modes it was built from. */
function vdCompsText(){
    if(!vdComp.length)return 'no displacement yet';
    var parts=[];
    for(var i=0;i<vdComp.length;i++){
        var c=vdComp[i];
        parts.push('#'+(c.idx+1)+' '+(c.scale>=0?'+':'')+c.scale.toFixed(2));
    }
    return parts.join(' · ')+' · total '+vdMaxDisp.toFixed(2)+' Å';
}
/** Signature of the component list, for the row markers. */
function vdCompSig(){
    var s='';
    for(var i=0;i<vdComp.length;i++)s+=vdComp[i].idx+',';
    return s;
}
/** End the session bookkeeping; the atoms are left exactly as they are. */
function vdDropSession(popUndo){
    // The snapshot belongs to the session: when the geometry is back at the
    // equilibrium it is not an edit any more, unless something else pushed on
    // top of it in the meantime.
    if(popUndo&&vdUndoDepth===undoStack.length&&undoStack.length>0){undoStack.pop();updateUndoBtn()}
    vdActive=false;vdBase=null;vdComp.length=0;vdMaxDisp=0;vdUndoDepth=-1;
}
/** Open a session: the equilibrium snapshot and the one undo step it owns. */
function vdBegin(){
    vdBase=MD.atoms.map(function(a){return{x:a.x,y:a.y,z:a.z}});
    vdComp.length=0;vdActive=true;
    // One undo step per session, taken in front of the first coordinate it
    // writes — the same rule an Alt-drag follows.
    vdSuppress=true;pushUndo();vdSuppress=false;
    vdUndoDepth=undoStack.length;
}
/** Back to the equilibrium: every component is removed and the atoms return to
 *  the snapshot, so the session's own undo step goes with it. */
function vdClearAll(msg){
    if(!vdActive){vdSyncPanel();return}
    vdComp.length=0;
    vdWriteAtoms();            // nothing left to add: the atoms become vdBase
    updateVibrationPositions();
    vdDropSession(true);
    updateRotationCenter();
    vdSyncPanel();
    modeInfoEl.textContent=msg||'Manual displacement cleared - back to the equilibrium geometry';
}
/** Drop the session without touching the geometry: the structure it was
 *  measured from is being replaced (a frame, a molecule, an Undo). */
function vdEnd(){
    vdDropSession(false);
    vdSyncPanel();
}
/** Close the session and keep the geometry it produced. Called in front of
 *  every coordinate edit, so the edit starts from a structure nothing else is
 *  holding a rewrite for, and the rotation centre describes what is on screen. */
function vdCommit(){
    if(!vdActive)return;
    if(vibActive||vibPaused)stopVibration();
    // The components are about to be cleared, but the displacement stays in the
    // coordinates — so it is kept for the title line before the session closes.
    provCommitModes();
    vdDropSession(false);
    // The centre follows the geometry; meshes are stored relative to it, so
    // this re-centres without moving the picture a pixel.
    updateRotationCenter();
    vdSyncPanel();
}
function vdReset(){
    if(vibActive||vibPaused)stopVibration();
    vdClearAll();
}
/** The manual half of a mode's selection: the circle in the list points the
 *  displacement slider at that mode. It plays nothing and writes nothing — the
 *  first slider move is what applies it — so picking a mode to adjust can never
 *  start an animation, and watching one can never re-aim the slider. */
function vdPick(idx){
    if(!NORMAL_MODES||idx<0||idx>=NORMAL_MODES.length)return;
    vdTarget=idx;
    vdSyncPanel();
    var m=NORMAL_MODES[idx];
    if(!vdUsable(idx)){
        modeInfoEl.textContent='Manual displacement - this mode has no displacements for the current structure';
        return;
    }
    modeInfoEl.textContent='Manual displacement: mode '+(idx+1)+
        (m&&m.frequency!=null?(' ('+m.frequency.toFixed(2)+' cm⁻¹)'):'')+
        ' - drag the slider to move the structure along it';
}
/** The one entry point of the slider: one mode's fraction, -1 … +1. */
function vdSetScale(scale){
    if(!vdAvailable()){vdSyncPanel();return}
    // Two writers may already own the coordinates. A Molecular Mechanics run
    // streams them chunk by chunk and a Move Atoms session is mid-edit; a slider
    // move under either would be overwritten (or would overwrite the run), so it
    // is refused and the thumb is put back where the geometry actually is. The
    // run keeps going — a slider is not a reason to end it.
    if(mmRunActive()){
        modeInfoEl.textContent='Manual displacement - a Molecular Mechanics run is using these coordinates';
        vdSyncPanel();return;
    }
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    if(!isFinite(scale))scale=0;
    if(scale>VD_MAX_SCALE)scale=VD_MAX_SCALE;
    if(scale<-VD_MAX_SCALE)scale=-VD_MAX_SCALE;
    // A playback owns the coordinates while it runs; ending it first puts the
    // combination back, which is the geometry this edit is measured from.
    if(vibActive||vibPaused)stopVibration();
    var idx=vdTarget;
    if(idx<0){vdSyncPanel();return}
    // A repeated event with the same value (input + change, a keyboard repeat)
    // is not a new edit and writes nothing.
    if(Math.abs(scale-vdCurScale())<1e-12){vdSyncPanel();return}
    var zero=Math.abs(scale)<1e-9;
    if(zero&&!vdActive){vdSyncPanel();return}
    if(!vdActive){
        if(!vdUsable(idx)){vdSyncPanel();return}
        vdBegin();
    }
    if(zero){
        // Dragging a component back to 0 takes that mode out of the
        // combination; when it was the last one the structure is back at the
        // equilibrium and the session is over.
        vdDropComp(idx);
        if(!vdComp.length){
            vdClearAll('Manual displacement cleared - back to the equilibrium geometry');
            return;
        }
    }else if(!vdSetComp(idx,scale)){
        if(!vdComp.length)vdDropSession(true);
        vdSyncPanel();return;
    }
    vdWriteAtoms();
    updateVibrationPositions();
    modeInfoEl.textContent='Manual displacement: #'+(idx+1)+' '+(scale>=0?'+':'')+scale.toFixed(2)+
        ' - '+vdComp.length+' mode'+(vdComp.length===1?'':'s')+', total max '+vdMaxDisp.toFixed(2)+' Å';
    vdSyncPanel();
}
/** One place that keeps the panel's displacement block honest: the mode the
 *  slider points at, the component list, whether the slider can be moved at
 *  all, the value, the reset button and the row markers. Safe to call when the
 *  panel is not open, and cheap enough to call on every slider event: nothing
 *  is written to the DOM unless it actually changed. */
function vdSyncPanel(){
    var panel=freqPanelEl;
    if(!panel)return;
    var slider=panel.querySelector('#vd-slider');
    if(!slider)return;
    var modeEl=panel.querySelector('#vd-mode');
    var compsEl=panel.querySelector('#vd-comps');
    var valEl=panel.querySelector('#vd-val');
    var maxEl=panel.querySelector('#vd-max');
    var resetBtn=panel.querySelector('#vd-reset');
    var hintEl=panel.querySelector('#vd-hint');
    var m=(NORMAL_MODES&&vdTarget>=0&&vdTarget<NORMAL_MODES.length)?NORMAL_MODES[vdTarget]:null;
    var usable=vdUsable(vdTarget);
    var enabled=usable&&vdAvailable();
    var off=!enabled;
    if(slider.disabled!==off)slider.disabled=off;
    var cur=vdCurScale();
    // Only write the value when it actually differs: this runs on every input
    // event of a drag, and re-assigning the thumb's own value is pointless work.
    if(parseFloat(slider.value)!==cur)slider.value=String(cur);
    var noReset=!vdActive;
    if(resetBtn&&resetBtn.disabled!==noReset)resetBtn.disabled=noReset;
    if(valEl){
        var valTxt=(cur<0?'-':'+')+Math.abs(cur).toFixed(2);
        if(valEl.textContent!==valTxt)valEl.textContent=valTxt;
    }
    // The fraction alone does not say how big the distortion is, so the same
    // row spells out what this mode's share comes to for the atom that moves
    // most within it (the unit vector makes that exactly |share| Å).
    if(maxEl){
        var maxTxt='max '+Math.abs(cur).toFixed(2)+' Å';
        if(maxEl.textContent!==maxTxt)maxEl.textContent=maxTxt;
    }
    if(compsEl){
        var compsTxt=vdCompsText();
        if(compsEl.textContent!==compsTxt)compsEl.textContent=compsTxt;
    }
    if(modeEl){
        var modeTxt=m?('mode '+(vdTarget+1)+' · '+m.frequency.toFixed(2)+' cm⁻¹'):'no mode selected';
        if(modeEl.textContent!==modeTxt)modeEl.textContent=modeTxt;
    }
    if(hintEl){
        var hintTxt;
        if(!NORMAL_MODES||!NORMAL_MODES.length)hintTxt='No normal modes in this document.';
        else if(!vdAvailable())hintTxt=diffMode?'Not available in diff mode.':'Not available for crystal structures.';
        else if(!m)hintTxt='Pick a mode with the circle, then drag to move the structure along it.';
        else if(!usable)hintTxt='This mode has no displacements for the current structure.';
        else if(vdComp.length&&Math.abs(cur)>1e-9)hintTxt='In the distortion: drag to change its share, 0 takes it out.';
        else if(vdComp.length)hintTxt='Drag to add this mode to the distortion.';
        else hintTxt='Drag along the mode: 1.00 moves the most-moving atom by 1 Å.';
        if(hintEl.textContent!==hintTxt)hintEl.textContent=hintTxt;
    }
    // The markers change only when the set of components or the picked mode
    // changes, so the row list is walked on that transition and not on every
    // slider event.
    var sig=vdCompSig();
    var picked=vdTarget;
    if(sig!==vdMarkedSig||picked!==vdPickedRow){
        vdMarkedSig=sig;vdPickedRow=picked;
        panel.querySelectorAll('.freq-row').forEach(function(r){
            var idx=parseInt(r.dataset.idx,10);
            r.classList.toggle('displaced',!!vdFindComp(idx));
            var pick=r.querySelector('.freq-pick');
            if(pick)pick.classList.toggle('on',idx===picked);
        });
    }
}

function buildFreqPanel(){
    if(!NORMAL_MODES||NORMAL_MODES.length===0)return;
    var rowsHTML=NORMAL_MODES.map(function(m,idx){
        var freqStr=m.frequency.toFixed(2);
        var symStr=m.symmetry||'';
        var sign=m.frequency<0?'*':'';
        return '<div class="freq-row" data-idx="'+idx+'">'+
            '<span class="freq-idx">'+(idx+1)+'</span>'+
            '<span class="freq-val">'+sign+freqStr+' cm⁻¹</span>'+
            '<span class="freq-sym">'+symStr+'</span>'+
            '<span class="freq-pick" data-idx="'+idx+'" data-pick="1" title="Select this mode for manual displacement"></span>'+
            '<span class="freq-play" title="Play this vibration">▶</span></div>';
    }).join('');
    freqPanelEl.innerHTML=
        '<div class="freq-head"><span class="freq-title">Normal Modes ('+NORMAL_MODES.length+')</span>'+
        '<span class="freq-close">×</span></div>'+
        '<div class="freq-body">'+rowsHTML+'</div>'+
        '<button class="freq-stop" id="freq-stop">⏸ Pause</button>'+
        '<div class="freq-dis">'+
        '<div class="vd-head"><span>Manual Displacement</span><span class="vd-mode" id="vd-mode">no mode selected</span></div>'+
        '<div class="vd-comps" id="vd-comps">no displacement yet</div>'+
        '<input type="range" class="vd-slider" id="vd-slider" min="-'+VD_MAX_SCALE+'" max="'+VD_MAX_SCALE+'" step="0.01" value="0">'+
        '<div class="vd-row"><span class="vd-val" id="vd-val">+0.00</span>'+
        '<span class="vd-max" id="vd-max">max 0.00 Å</span>'+
        '<button class="vd-reset" id="vd-reset" disabled>Reset all</button></div>'+
        '<div class="vd-hint" id="vd-hint">Pick a mode with the circle, then drag to move the structure along it.</div>'+
        '</div>';
    freqPanelEl.classList.add('show');
    freqReopenEl.classList.remove('show');

    var rows=freqPanelEl.querySelectorAll('.freq-row');
    rows.forEach(function(r){
        r.addEventListener('click',function(e){
            // The circle inside the row is the manual-selection control; its
            // own handler stops the event, and this guard keeps the two apart
            // even if it does not.
            if(e&&e.target&&e.target.dataset&&e.target.dataset.pick)return;
            var idx=parseInt(r.dataset.idx);
            if(isNaN(idx))return;
            if((vibActive||vibPaused)&&vibModeIdx===idx){
                stopVibration();
            }else{
                stopVibration();
                startVibration(idx);
                if(vibActive){
                    rows.forEach(function(rr){rr.classList.remove('playing')});
                    r.classList.add('playing');
                }
            }
        });
    });
    // The other half of the selection: pick a mode for the displacement slider.
    // It plays nothing and writes nothing; the first slider move does.
    freqPanelEl.querySelectorAll('.freq-pick').forEach(function(p){
        p.addEventListener('click',function(e){
            if(e&&typeof e.stopPropagation==='function')e.stopPropagation();
            var idx=parseInt(this.dataset.idx,10);
            if(isNaN(idx))return;
            vdPick(idx);
        });
    });
    var closeBtn=freqPanelEl.querySelector('.freq-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            stopVibration();
            freqPanelEl.classList.remove('show');
            freqReopenEl.classList.add('show');
            layoutPanels();
        });
    }
    var stopBtn=freqPanelEl.querySelector('#freq-stop');
    if(stopBtn){
        stopBtn.addEventListener('click',function(){
            if(vibActive){pauseVibration()}
            else if(vibPaused){resumeVibration()}
        });
    }
    // Manual displacement: the slider writes the mode currently named in the
    // block. 'input' covers a drag, a click on the track and an arrow key on a
    // focused slider; the handler itself drops repeated identical values.
    var vdSlider=freqPanelEl.querySelector('#vd-slider');
    if(vdSlider){
        vdSlider.addEventListener('input',function(){vdSetScale(parseFloat(this.value))});
        vdSlider.addEventListener('change',function(){vdSetScale(parseFloat(this.value))});
    }
    var vdResetBtn=freqPanelEl.querySelector('#vd-reset');
    if(vdResetBtn)vdResetBtn.addEventListener('click',function(){vdReset()});
    // The rows were just recreated, so the marker caches start empty.
    vdMarkedSig='';
    vdPickedRow=-1;
    vdSyncPanel();
    layoutPanels();
}

optReopenEl.addEventListener('click',function(){buildOptPanel()});
freqReopenEl.addEventListener('click',function(){buildFreqPanel()});

// ===== Route Keywords Panel (Gaussian route cards, ORCA simple input) =====
var ROUTES=MD.routes||null;
function rtEsc(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;')}
function buildRoutePanel(){
    if(!ROUTES||ROUTES.length===0)return;
    var multi=ROUTES.length>1;
    var html='<div class="rt-head"><span class="rt-title">Route Keywords</span><span class="rt-close">×</span></div><div class="rt-body">';
    ROUTES.forEach(function(sec,si){
        if(multi)html+='<div class="rt-job">Job '+(si+1)+(sec.hasOpt?'<span class="rt-badge">opt</span>':'')+'</div>';
        html+='<div class="rt-raw">'+rtEsc(sec.raw)+'</div>';
        html+=sec.keywords.map(function(k){
            var opts=(k.options&&k.options.length>0)?'<span class="rt-kopts">= '+rtEsc(k.options.join(', '))+'</span>':'';
            return '<div class="rt-kw"><span class="rt-kname">'+rtEsc(k.name)+'</span>'+opts+'</div>';
        }).join('');
    });
    html+='</div>';
    routePanelEl.innerHTML=html;
    routePanelEl.classList.add('show');
    routeReopenEl.classList.remove('show');
    var closeBtn=routePanelEl.querySelector('.rt-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            routePanelEl.classList.remove('show');
            routeReopenEl.classList.add('show');
            layoutPanels();
        });
    }
    layoutPanels();
}
routeReopenEl.addEventListener('click',function(){buildRoutePanel()});

// Initial panel display — show reopen buttons only, don't auto-open panels
if(OPT_STEPS&&OPT_STEPS.length>0){
    optReopenEl.classList.add('show');
}
if(NORMAL_MODES&&NORMAL_MODES.length>0){
    freqReopenEl.classList.add('show');
}
if(ROUTES&&ROUTES.length>0){
    routeReopenEl.classList.add('show');
}
layoutPanels();

function formatAtomList(indices,atoms){
    atoms=atoms||MD.atoms;
    if(indices.length===0)return '';
    var list=indices.slice(0,15).map(function(i){var a=atoms[i];return(a?a.element:'?')+(i+1)});
    var str=list.join(', ');
    if(indices.length>15)str+=' ... ('+indices.length+' atoms)';
    return str;
}
function getMoleculeScreenBox(){
    var container=document.getElementById('container');
    if(!container)return null;
    var w=container.clientWidth,h=container.clientHeight;
    if(w<2||h<2)return null;
    if(!atomMeshes||atomMeshes.length===0||!moleculeGroup)return null;
    scene.updateMatrixWorld(true);
    var minX=Infinity,minY=Infinity,maxX=-Infinity,maxY=-Infinity;
    function projectGroup(group,camW,offX){
        if(!group||!group.children.length)return;
        var box=new THREE.Box3().setFromObject(group);
        if(box.isEmpty())return;
        var cs=[
            [box.min.x,box.min.y,box.min.z],[box.min.x,box.min.y,box.max.z],
            [box.min.x,box.max.y,box.min.z],[box.min.x,box.max.y,box.max.z],
            [box.max.x,box.min.y,box.min.z],[box.max.x,box.min.y,box.max.z],
            [box.max.x,box.max.y,box.min.z],[box.max.x,box.max.y,box.max.z]
        ];
        var v=new THREE.Vector3();
        for(var k=0;k<cs.length;k++){
            v.set(cs[k][0],cs[k][1],cs[k][2]);
            v.project(camera);
            if(v.z>-1&&v.z<1){
                var sx=offX+(v.x*0.5+0.5)*camW;
                var sy=(-v.y*0.5+0.5)*h;
                if(sx<minX)minX=sx; if(sx>maxX)maxX=sx;
                if(sy<minY)minY=sy; if(sy>maxY)maxY=sy;
            }
        }
    }
    var sa=camera.aspect,sp=camera.position.clone();
    try{
        if(diffMode&&diffPivot&&diffAtomMeshes.length>0){
            var halfW=Math.floor(w/2);
            camera.aspect=halfW/h;camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);camera.updateMatrixWorld(true);
            projectGroup(moleculeGroup,halfW,0);
            camera.aspect=(w-halfW)/h;camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);camera.updateMatrixWorld(true);
            projectGroup(diffMolGroup,(w-halfW),halfW);
        }else{
            camera.aspect=w/h;camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);camera.updateMatrixWorld(true);
            projectGroup(moleculeGroup,w,0);
        }
    }finally{
        camera.aspect=sa;camera.updateProjectionMatrix();
        camera.position.copy(sp);camera.updateMatrixWorld(true);
    }
    if(minX===Infinity)return null;
    var pad=12;
    return {left:Math.max(0,minX-pad),top:Math.max(0,minY-pad),
            right:Math.min(w,maxX+pad),bottom:Math.min(h,maxY+pad)};
}
function layoutPanels(){
    // The Molecular Mechanics reopen button follows the mode (and the run) the
    // way the other reopen buttons follow their panels, so it is refreshed
    // before the layout measures anything.
    mmUpdateReopen();
    var container=document.getElementById('container');
    if(!container){mmDrawCurves();return;}
    var cRect=container.getBoundingClientRect();
    var cTop=cRect.top,cLeft=cRect.left,cW=cRect.width,cH=cRect.height;
    var molBox=getMoleculeScreenBox();
    var molVR=null;
    if(molBox){
        molVR={left:cLeft+molBox.left,top:cTop+molBox.top,
               right:cLeft+molBox.right,bottom:cTop+molBox.bottom};
    }
    var defs=[
        {id:'mol-info',pr:0,anc:'tl'},
        {id:'mm-panel',pr:1,anc:'tl'},
        {id:'crystal-panel',pr:1,anc:'tr'},
        // The egg list shares the top-right corner with the supercell panel and
        // the Diff reopen button shares the bottom-right with the View cards;
        // both used to be absent from this list, so a panel laid out over them
        // covered them instead of giving way.
        {id:'egg-panel',pr:1,anc:'tr'},
        {id:'diff-panel',pr:2,anc:'br'},
        {id:'diff-reopen',pr:2,anc:'br'},
        // Auto Rotate and Edit Coordinates are two cards, each with its own
        // collision box, stacked bottom-up in that reading order.
        {id:'bake-panel',pr:3,anc:'br'},
        {id:'view-panel',pr:4,anc:'br'},
        {id:'select-panel',pr:3,anc:'bl'},
        {id:'rotate-panel',pr:3,anc:'bl'},
        {id:'diff-label',pr:4,anc:'tc',pct:0.25},
        {id:'diff-label-right',pr:5,anc:'tc',pct:0.75},
        {id:'opt-panel',pr:7,anc:'tl'},
        {id:'freq-panel',pr:8,anc:'tl'},
        {id:'route-panel',pr:8,anc:'tl'},
        {id:'mm-reopen',pr:9,anc:'tl'},
        {id:'opt-reopen',pr:9,anc:'tl'},
        {id:'freq-reopen',pr:10,anc:'tl'},
        {id:'route-reopen',pr:10,anc:'tl'},
        {id:'axes-indicator',pr:6,anc:'bl'}
    ];
    var vis=[];
    defs.forEach(function(d){
        var el=document.getElementById(d.id);
        if(!el||window.getComputedStyle(el).display==='none')return;
        el.style.top='';el.style.bottom='';el.style.left='';el.style.right='';
        el.style.transform='';el.style.maxWidth='';el.style.fontSize='';
        vis.push({d:d,el:el});
    });
    vis.sort(function(a,b){return a.d.pr-b.d.pr;});
    function rect(el){return el.getBoundingClientRect();}
    function ov(a,b){return !(a.right<=b.left+1||a.left>=b.right-1||a.bottom<=b.top+1||a.top>=b.bottom-1);}
    for(var i=0;i<vis.length;i++){
        var cur=vis[i];
        for(var pass=0;pass<8;pass++){
            var cr=rect(cur.el);
            var conf=null,cfR=null,isMol=false;
            if(molVR&&ov(cr,molVR)){conf='mol';cfR=molVR;isMol=true}
            if(!conf){
                for(var j=0;j<i;j++){
                    var rr=rect(vis[j].el);
                    if(ov(cr,rr)){conf=vis[j];cfR=rr;break;}
                }
            }
            if(!conf)break;
            var a=cur.d.anc;
            var resolved=false;
            if(isMol&&a==='tc'&&!resolved){
                var cc=(cr.left+cr.right)/2;
                var contC=cLeft+cW/2;
                if(cc<contC){
                    var tl=molBox.left-cr.width-4;
                    if(tl>=2){cur.el.style.left=tl+'px';cur.el.style.transform='none';cur.el.style.right='auto';resolved=true;}
                }else{
                    var tr=cW-molBox.right-cr.width-4;
                    if(tr>=2){cur.el.style.right=tr+'px';cur.el.style.left='auto';cur.el.style.transform='none';resolved=true;}
                }
            }
            if(!resolved&&!isMol&&(a==='tc'||a==='tl'||a==='tr')){
                var nt=cfR.bottom-cTop+4;
                if(nt+cr.height<cH-8){cur.el.style.top=nt+'px';resolved=true;}
            }
            if(!resolved&&!isMol&&(a==='bl'||a==='br')){
                var nb=cH-(cfR.top-cTop)+4;
                if(nb+cr.height<cH-8){cur.el.style.bottom=nb+'px';resolved=true;}
            }
            if(!resolved){
                var mw=parseInt(window.getComputedStyle(cur.el).maxWidth)||1000;
                if(mw>80){cur.el.style.maxWidth=Math.max(80,mw-30)+'px';resolved=true;}
            }
            if(!resolved){
                var fs=parseInt(window.getComputedStyle(cur.el).fontSize)||11;
                if(fs>8){cur.el.style.fontSize=(fs-1)+'px';resolved=true;}
            }
            if(!resolved&&!isMol&&conf.el){
                var cmw=parseInt(window.getComputedStyle(conf.el).maxWidth)||1000;
                if(cmw>80){conf.el.style.maxWidth=Math.max(80,cmw-30)+'px';resolved=true;}
            }
            if(!resolved&&!isMol&&conf.el){
                var cfs=parseInt(window.getComputedStyle(conf.el).fontSize)||11;
                if(cfs>8){conf.el.style.fontSize=(cfs-1)+'px';resolved=true;}
            }
            if(!resolved)break;
        }
    }
    // A reflow (a window resize, the panel being moved or restyled) changes the
    // box the curves are drawn into, and they are sized from that box — so they
    // are redrawn here rather than left at the width of the previous layout.
    // A no-op unless the Molecular Mechanics panel is up with a run behind it.
    mmDrawCurves();
}
// The highlight is a *state* picked from the shared material pool (plain,
// selected, diff-highlighted), never a field written into a material — every
// atom with the same colour and state draws with the same material, and a
// per-atom material would be one material per atom again.
function highlightSelected(){
    // One lookup for the whole pass instead of one scan per atom: Select All
    // makes the selection the entire structure, and an indexOf per mesh would
    // pay for that with a quadratic walk (135 ms at 20 000 atoms measured).
    // Same answer as before — an index the structure does not have never
    // matches, exactly as indexOf returned -1.
    var sel={},i;
    for(i=0;i<selectedAtoms.length;i++)sel[selectedAtoms[i]]=true;
    atomMeshes.forEach(function(m,i){
        setAtomHighlight(m,sel[i]===true?1:(m.userData.diffHi?2:0));
    });
    if(diffMode){
        var dsel={};
        for(i=0;i<diffSelectedAtoms.length;i++)dsel[diffSelectedAtoms[i]]=true;
        diffAtomMeshes.forEach(function(m,i){
            setDiffAtomHighlight(m,dsel[i]===true?1:(m.userData.diffHi?2:0));
        });
    }
    needsRender=true;
}

function selectAtom(idx){
    if(diffMode&&diffActiveSide==='right'){
        var need=diffSelectedAtoms.length>=requiredCount();
        if(need){diffSelectedAtoms=[idx]}else if(diffSelectedAtoms.indexOf(idx)>=0){return}else{diffSelectedAtoms.push(idx)}
        if(diffMapping){
            var needL=selectedAtoms.length>=requiredCount();
            if(needL){selectedAtoms=[]}
            diffSelectedAtoms.forEach(function(i){
                var li=diffReverseMapping[i];
                if(li>=0&&selectedAtoms.indexOf(li)<0)selectedAtoms.push(li);
            });
        }
        highlightSelected();
        updateDiffSelInfo();
        checkSelectionComplete();
        return;
    }
    if(diffMode){
        var needL=selectedAtoms.length>=requiredCount();
        if(needL){selectedAtoms=[idx]}else if(selectedAtoms.indexOf(idx)>=0){return}else{selectedAtoms.push(idx)}
        if(diffMapping){
            var needR=diffSelectedAtoms.length>=requiredCount();
            if(needR){diffSelectedAtoms=[]}
            selectedAtoms.forEach(function(i){
                var ri=diffMapping[i];
                if(ri>=0&&diffSelectedAtoms.indexOf(ri)<0)diffSelectedAtoms.push(ri);
            });
        }
        highlightSelected();
        updateDiffSelInfo();
        checkSelectionComplete();
        return;
    }
    if(currentMode==='selectAtoms'||currentMode==='replaceAtom'||currentMode==='addH'||currentMode==='rotateGroup'||currentMode==='moveAtoms'||currentMode==='mmOpt'){
        var pos=selectedAtoms.indexOf(idx);
        if(pos>=0){
            selectedAtoms.splice(pos,1);
            highlightSelected();
            var dn=formatAtomList(selectedAtoms);
            selInfoEl.textContent=selectedAtoms.length?('Selected: '+dn):'';
            if(currentMode==='rotateGroup')syncRotatePanelGroupField();
            checkSelectionComplete();
            mmSelectionChanged();
            return;
        }
    }
    if(selectedAtoms.indexOf(idx)>=0)return;
    selectedAtoms.push(idx);
    highlightSelected();
    var names=formatAtomList(selectedAtoms);
    selInfoEl.textContent='Selected: '+names;
    if(currentMode==='rotateGroup')syncRotatePanelGroupField();
    checkSelectionComplete();
    mmSelectionChanged();
}

function updateDiffSelInfo(){
    var lNames=formatAtomList(selectedAtoms);
    var rNames=formatAtomList(diffSelectedAtoms,diffData.atoms);
    selInfoEl.textContent='[Left] '+(lNames||'-')+'   |   [Right] '+(rNames||'-');
}

function requiredCount(){
    if(currentMode==='bondLength')return 2;
    if(currentMode==='bondAngle')return 3;
    if(currentMode==='dihedral')return 4;
    if(currentMode==='addAtom'||currentMode==='deleteAtom')return 1;
    return 999;
}

function checkSelectionComplete(){
    if(diffMode){
        if(currentMode==='bondLength'){
            if(selectedAtoms.length===2||diffSelectedAtoms.length===2)showDiffMeasurement('bondLength');
        }else if(currentMode==='bondAngle'){
            if(selectedAtoms.length===3||diffSelectedAtoms.length===3)showDiffMeasurement('bondAngle');
        }else if(currentMode==='dihedral'){
            if(selectedAtoms.length===4||diffSelectedAtoms.length===4)showDiffMeasurement('dihedral');
        }else if(currentMode==='selectAtoms'||currentMode==='boxSelect'||currentMode==='replaceAtom'){
            updateDiffSelInfo();
        }else if(currentMode==='deleteAtom'&&selectedAtoms.length===1){
            showDeleteAtomModal();
        }else if(currentMode==='addAtom'&&selectedAtoms.length===1){
            showAddAtomModal();
        }
        return;
    }
    if(currentMode==='bondLength'&&selectedAtoms.length===2)showBondLengthModal();
    else if(currentMode==='bondAngle'&&selectedAtoms.length===3)showBondAngleModal();
    else if(currentMode==='dihedral'&&selectedAtoms.length===4)showDihedralModal();
    else if(currentMode==='addAtom'&&selectedAtoms.length===1)showAddAtomModal();
    else if(currentMode==='deleteAtom'&&selectedAtoms.length===1)showDeleteAtomModal();
    else if(currentMode==='selectAtoms'){
        var names=formatAtomList(selectedAtoms);
        selInfoEl.textContent='Selected: '+names;
    }
    else if(currentMode==='replaceAtom'){
        var rnames=formatAtomList(selectedAtoms);
        selInfoEl.textContent='Selected to replace: '+rnames+' (click Replace Atom button to confirm)';
    }
    else if(currentMode==='addH'){
        var hnames=formatAtomList(selectedAtoms);
        selInfoEl.textContent='Selected for Add H: '+hnames+' (click Add H button to complete)';
    }
    else if(currentMode==='mmOpt'){
        var mnames=formatAtomList(selectedAtoms);
        selInfoEl.textContent='Selected for MM optimization: '+mnames+' ('+selectedAtoms.length+' atoms, click Start Optimization)';
    }
}

function fmtMeasurement(kind,sel,atomsArr){
    if(!sel||sel.length===0)return null;
    if(kind==='bondLength'&&sel.length===2){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]];
        var d=dist(a1,a2);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1),value:d.toFixed(4)+' A',num:d};
    }
    if(kind==='bondAngle'&&sel.length===3){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]];
        var ang=angle(a1,a2,a3);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1)+' - '+a3.element+(sel[2]+1),value:ang.toFixed(2)+' deg',num:ang};
    }
    if(kind==='dihedral'&&sel.length===4){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]],a4=atomsArr[sel[3]];
        var dih=dihedral(a1,a2,a3,a4);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1)+' - '+a3.element+(sel[2]+1)+' - '+a4.element+(sel[3]+1),value:dih.toFixed(2)+' deg',num:dih};
    }
    return null;
}

function showDiffMeasurement(kind){
    var kindLabel=kind==='bondLength'?'Bond Length':kind==='bondAngle'?'Bond Angle':'Dihedral Angle';
    var left=fmtMeasurement(kind,selectedAtoms,MD.atoms);
    var right=fmtMeasurement(kind,diffSelectedAtoms,diffData.atoms);
    var html='<span class="diff-close">×</span><h4>'+kindLabel+' (Diff Mode)</h4>';
    html+='<div class="diff-row" style="color:#4ec9b0"><b>Left: '+(MD.title||'original')+'</b></div>';
    if(left){
        html+='<div class="diff-row">  '+left.label+': '+left.value+'</div>';
    }else{
        html+='<div class="diff-row" style="color:var(--vscode-descriptionForeground,#999)">  Select '+(kind==='bondLength'?2:kind==='bondAngle'?3:4)+' atoms on left...</div>';
    }
    html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Right: '+diffData.fileName+'</b></div>';
    if(right){
        html+='<div class="diff-row">  '+right.label+': '+right.value+'</div>';
    }else{
        html+='<div class="diff-row" style="color:var(--vscode-descriptionForeground,#999)">  Select '+(kind==='bondLength'?2:kind==='bondAngle'?3:4)+' atoms on right...</div>';
    }
    if(left&&right){
        var diff=Math.abs(left.num-right.num);
        if(kind==='dihedral'&&diff>180)diff=360-diff;
        html+='<div class="diff-row" style="margin-top:6px;color:#ff6600"><b>Δ = '+diff.toFixed(4)+(kind==='bondLength'?' A':' deg')+'</b></div>';
    }
    html+='<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Selecting atoms on one side auto-selects the corresponding atoms on the other.</div>';
    showDiffPanel(html);
}

function dist(a,b){var dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;return Math.sqrt(dx*dx+dy*dy+dz*dz)}
function angle(a,b,c){var v1={x:a.x-b.x,y:a.y-b.y,z:a.z-b.z},v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var d1=Math.sqrt(v1.x*v1.x+v1.y*v1.y+v1.z*v1.z),d2=Math.sqrt(v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(d1*d2);
    return Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI}

function dihedral(a,b,c,d){
    var b1={x:b.x-a.x,y:b.y-a.y,z:b.z-a.z},b2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z},b3={x:d.x-c.x,y:d.y-c.y,z:d.z-c.z};
    function cross(u,v){return{x:u.y*v.z-u.z*v.y,y:u.z*v.x-u.x*v.z,z:u.x*v.y-u.y*v.x}}
    function dot(u,v){return u.x*v.x+u.y*v.y+u.z*v.z}
    function norm(v){var l=Math.sqrt(v.x*v.x+v.y*v.y+v.z*v.z);return l<1e-10?{x:0,y:0,z:0}:{x:v.x/l,y:v.y/l,z:v.z/l}}
    var n1=cross(b1,b2),n2=cross(b2,b3);
    var m=cross(norm(n1),norm(b2));
    var x=dot(n1,n2),y=dot(m,n2);
    return-Math.atan2(y,x)*180/Math.PI;
}

function getMovable(fixedSet,startIdx){
    var adj=[];for(var i=0;i<MD.atoms.length;i++)adj[i]=[];
    MD.bonds.forEach(function(b){adj[b.atom1].push(b.atom2);adj[b.atom2].push(b.atom1)});
    var visited=new Set(fixedSet);visited.add(startIdx);
    var queue=[startIdx],result=[startIdx];
    while(queue.length>0){var cur=queue.shift();adj[cur].forEach(function(nb){
        if(!visited.has(nb)){visited.add(nb);queue.push(nb);result.push(nb)}})}
    return result;
}

function rotAroundAxis(px,py,pz,ox,oy,oz,dx,dy,dz,angle){
    var c=Math.cos(angle),s=Math.sin(angle);
    var x=px-ox,y=py-oy,z=pz-oz;
    var kx=dx,ky=dy,kz=dz;
    var l=Math.sqrt(kx*kx+ky*ky+kz*kz);if(l<1e-10)return{x:px,y:py,z:pz};
    kx/=l;ky/=l;kz/=l;
    var dot=x*kx+y*ky+z*kz;
    var rx=x*c+(ky*z-kz*y)*s+kx*dot*(1-c);
    var ry=y*c+(kz*x-kx*z)*s+ky*dot*(1-c);
    var rz=z*c+(kx*y-ky*x)*s+kz*dot*(1-c);
    return{x:rx+ox,y:ry+oy,z:rz+oz};
}

function saveOriginal(){originalCoords=MD.atoms.map(function(a){return{x:a.x,y:a.y,z:a.z}});pushUndo()}
function restoreOriginal(){if(!originalCoords)return;originalCoords.forEach(function(c,i){var a=MD.atoms[i];if(!a)return;a.x=c.x;a.y=c.y;a.z=c.z})}

function applyBondLength(targetLen,fixFirst,moveMode){
    var i1=selectedAtoms[0],i2=selectedAtoms[1];
    var o1=originalCoords[i1],o2=originalCoords[i2];
    var dx=o2.x-o1.x,dy=o2.y-o1.y,dz=o2.z-o1.z;
    var curLen=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(curLen<1e-10)return;
    if(fixFirst){
        var newI2x=o1.x+dx*targetLen/curLen;
        var newI2y=o1.y+dy*targetLen/curLen;
        var newI2z=o1.z+dz*targetLen/curLen;
        var tx=newI2x-o2.x,ty=newI2y-o2.y,tz=newI2z-o2.z;
        var fixedSet=new Set([i1]);
        var movable=(moveMode==='atom')?[i2]:getMovable(fixedSet,i2);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            MD.atoms[idx].x=oc.x+tx;
            MD.atoms[idx].y=oc.y+ty;
            MD.atoms[idx].z=oc.z+tz;
        });
    }else{
        var newI1x=o2.x-dx*targetLen/curLen;
        var newI1y=o2.y-dy*targetLen/curLen;
        var newI1z=o2.z-dz*targetLen/curLen;
        var tx=newI1x-o1.x,ty=newI1y-o1.y,tz=newI1z-o1.z;
        var fixedSet=new Set([i2]);
        var movable=(moveMode==='atom')?[i1]:getMovable(fixedSet,i1);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            MD.atoms[idx].x=oc.x+tx;
            MD.atoms[idx].y=oc.y+ty;
            MD.atoms[idx].z=oc.z+tz;
        });
    }
    updateScenePositions(true);
}

function applyBondAngle(targetDeg,fixFirstTwo,moveMode){
    var i1=selectedAtoms[0],i2=selectedAtoms[1],i3=selectedAtoms[2];
    var curDeg=angle(originalCoords[i1],originalCoords[i2],originalCoords[i3]);
    var delta=(targetDeg-curDeg)*Math.PI/180;
    if(Math.abs(delta)<1e-10)return;
    var v1={x:originalCoords[i1].x-originalCoords[i2].x,y:originalCoords[i1].y-originalCoords[i2].y,z:originalCoords[i1].z-originalCoords[i2].z};
    var v2={x:originalCoords[i3].x-originalCoords[i2].x,y:originalCoords[i3].y-originalCoords[i2].y,z:originalCoords[i3].z-originalCoords[i2].z};
    var cx=v1.y*v2.z-v1.z*v2.y,cy=v1.z*v2.x-v1.x*v2.z,cz=v1.x*v2.y-v1.y*v2.x;
    var cl=Math.sqrt(cx*cx+cy*cy+cz*cz);if(cl<1e-10)return;
    cx/=cl;cy/=cl;cz/=cl;
    var ox=originalCoords[i2].x,oy=originalCoords[i2].y,oz=originalCoords[i2].z;
    if(fixFirstTwo){
        var fixedSet=new Set([i1,i2]);
        var movable=(moveMode==='atom')?[i3]:getMovable(fixedSet,i3);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,cx,cy,cz,delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }else{
        var fixedSet=new Set([i2,i3]);
        var movable=(moveMode==='atom')?[i1]:getMovable(fixedSet,i1);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,cx,cy,cz,-delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }
    updateScenePositions(true);
}

function applyDihedral(targetDeg,fixFirstThree,moveMode){
    var i1=selectedAtoms[0],i2=selectedAtoms[1],i3=selectedAtoms[2],i4=selectedAtoms[3];
    var curDeg=dihedral(originalCoords[i1],originalCoords[i2],originalCoords[i3],originalCoords[i4]);
    var delta=(targetDeg-curDeg)*Math.PI/180;
    if(Math.abs(delta)<1e-10)return;
    var ax=originalCoords[i3].x-originalCoords[i2].x,ay=originalCoords[i3].y-originalCoords[i2].y,az=originalCoords[i3].z-originalCoords[i2].z;
    var al=Math.sqrt(ax*ax+ay*ay+az*az);if(al<1e-10)return;
    ax/=al;ay/=al;az/=al;
    var ox=originalCoords[i2].x,oy=originalCoords[i2].y,oz=originalCoords[i2].z;
    if(fixFirstThree){
        var fixedSet=new Set([i1,i2,i3]);
        var movable=(moveMode==='atom')?[i4]:getMovable(fixedSet,i4);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,ax,ay,az,delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }else{
        var fixedSet=new Set([i2,i3,i4]);
        var movable=(moveMode==='atom')?[i1]:getMovable(fixedSet,i1);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,ax,ay,az,-delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }
    updateScenePositions(true);
}

function showModal(html,cb){if(vibActive||vibPaused)stopVibration();stopMoveSession();modalEl.innerHTML=html;modalOverlay.classList.add('show');modalCallback=cb}
function hideModal(){modalOverlay.classList.remove('show');modalCallback=null}

function showBondLengthModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]];
    var cur=dist(a1,a2);
    var existingBond=MD.bonds.find(function(b){return(b.atom1===selectedAtoms[0]&&b.atom2===selectedAtoms[1])||(b.atom1===selectedAtoms[1]&&b.atom2===selectedAtoms[0])});
    var curOrder=existingBond?existingBond.order:1;
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1);
    showModal('<h3>Adjust Bond Length</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(4)+' A, Bond order: '+curOrder+'</div>'+
        '<label>Fix atom:</label><select id="m-fix"><option value="1">Fix '+n1+' (move '+n2+')</option><option value="2">Fix '+n2+' (move '+n1+')</option></select>'+
        '<label>Move:</label><select id="m-move"><option value="group">Translate Group (atom + connected)</option><option value="atom">Translate Atom (only)</option></select>'+
        '<label>Bond order:</label><select id="m-order"><option value="0"'+(curOrder===0?' selected':'')+'>None (0) - Remove bond</option><option value="1"'+(curOrder===1?' selected':'')+'>Single (1.0)</option><option value="1.5"'+(curOrder===1.5?' selected':'')+'>Aromatic (1.5)</option><option value="2"'+(curOrder===2?' selected':'')+'>Double (2.0)</option><option value="3"'+(curOrder===3?' selected':'')+'>Triple (3.0)</option></select>'+
        '<label>Target length (A):</label><input type="number" id="m-val" value="'+cur.toFixed(4)+'" step="0.01" min="0.3" max="6">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(4)+'" min="0.3" max="6" step="0.01">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix'),orderEl=document.getElementById('m-order'),moveEl=document.getElementById('m-move');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyBondLength(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyBondLength(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    fixEl.addEventListener('change',function(){applyBondLength(parseFloat(valEl.value),this.value==='1',moveEl.value)});
    moveEl.addEventListener('change',function(){applyBondLength(parseFloat(valEl.value),fixEl.value==='1',this.value)});
    orderEl.addEventListener('change',function(){
        var newOrder=parseFloat(this.value);
        if(newOrder===0){
            if(existingBond){MD.bonds=MD.bonds.filter(function(b){return b!==existingBond});existingBond=null;rebuildScene()}
        }else{
            if(!existingBond){
                existingBond={atom1:selectedAtoms[0],atom2:selectedAtoms[1],order:newOrder};
                MD.bonds.push(existingBond);
            }else{existingBond.order=newOrder}
            rebuildScene();
        }
    });
    document.getElementById('m-ok').addEventListener('click',function(){provMarkIfMoved('bondLen');hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){
        if(existingBond&&originalCoords){existingBond.order=curOrder}
        restoreOriginal();rebuildScene();undoStack.pop();updateUndoBtn();hideModal();originalCoords=null;resetSelection()
    });
}

function showBondAngleModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]];
    var cur=angle(a1,a2,a3);
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1),n3=a3.element+(selectedAtoms[2]+1);
    showModal('<h3>Adjust Bond Angle</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(2)+' deg</div>'+
        '<label>Fix side:</label><select id="m-fix"><option value="1">Fix '+n1+'-'+n2+' (move '+n3+')</option><option value="2">Fix '+n2+'-'+n3+' (move '+n1+')</option></select>'+
        '<label>Move:</label><select id="m-move"><option value="group">Translate Group (atom + connected)</option><option value="atom">Translate Atom (only)</option></select>'+
        '<label>Target angle (deg):</label><input type="number" id="m-val" value="'+cur.toFixed(2)+'" step="0.5" min="5" max="175">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(2)+'" min="5" max="175" step="0.5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix'),moveEl=document.getElementById('m-move');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyBondAngle(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyBondAngle(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    fixEl.addEventListener('change',function(){applyBondAngle(parseFloat(valEl.value),this.value==='1',moveEl.value)});
    moveEl.addEventListener('change',function(){applyBondAngle(parseFloat(valEl.value),fixEl.value==='1',this.value)});
    document.getElementById('m-ok').addEventListener('click',function(){provMarkIfMoved('bondAngle');hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){restoreOriginal();rebuildScene();undoStack.pop();updateUndoBtn();hideModal();originalCoords=null;resetSelection()});
}

function showDihedralModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]],a4=MD.atoms[selectedAtoms[3]];
    var cur=dihedral(a1,a2,a3,a4);
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1),n3=a3.element+(selectedAtoms[2]+1),n4=a4.element+(selectedAtoms[3]+1);
    showModal('<h3>Adjust Dihedral Angle</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(2)+' deg</div>'+
        '<label>Fix side:</label><select id="m-fix"><option value="1">Fix '+n1+'-'+n2+'-'+n3+' (move '+n4+')</option><option value="2">Fix '+n2+'-'+n3+'-'+n4+' (move '+n1+')</option></select>'+
        '<label>Move:</label><select id="m-move"><option value="group">Translate Group (atom + connected)</option><option value="atom">Translate Atom (only)</option></select>'+
        '<label>Target dihedral (deg):</label><input type="number" id="m-val" value="'+cur.toFixed(2)+'" step="1" min="-180" max="180">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(2)+'" min="-180" max="180" step="0.5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix'),moveEl=document.getElementById('m-move');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyDihedral(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyDihedral(parseFloat(this.value),fixEl.value==='1',moveEl.value)});
    fixEl.addEventListener('change',function(){applyDihedral(parseFloat(valEl.value),this.value==='1',moveEl.value)});
    moveEl.addEventListener('change',function(){applyDihedral(parseFloat(valEl.value),fixEl.value==='1',this.value)});
    document.getElementById('m-ok').addEventListener('click',function(){provMarkIfMoved('dihedral');hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){restoreOriginal();rebuildScene();undoStack.pop();updateUndoBtn();hideModal();originalCoords=null;resetSelection()});
}

var UNIT_CANDIDATES=(function(){
    var arr=[];var N=50;var phi=Math.PI*(3-Math.sqrt(5));
    for(var i=0;i<N;i++){
        var y=1-(i/(N-1))*2;
        var r=Math.sqrt(1-y*y);
        var th=phi*i;
        arr.push({x:Math.cos(th)*r,y:y,z:Math.sin(th)*r});
    }
    return arr;
})();
function pickFreeDirection(neighbors){
    if(!neighbors||neighbors.length===0)return {x:0,y:0,z:1};
    var best=UNIT_CANDIDATES[0];var bestMin=-1;
    for(var ci=0;ci<UNIT_CANDIDATES.length;ci++){
        var c=UNIT_CANDIDATES[ci];
        var minAng=Infinity;
        for(var ni=0;ni<neighbors.length;ni++){
            var nv=neighbors[ni];
            var dot=c.x*nv.x+c.y*nv.y+c.z*nv.z;
            var ang=Math.acos(Math.max(-1,Math.min(1,dot)));
            if(ang<minAng)minAng=ang;
        }
        if(minAng>bestMin){bestMin=minAng;best=c;}
    }
    return best;
}
function showAddAtomModal(){
    var anchorIdx=selectedAtoms[0];
    var anchor=MD.atoms[anchorIdx];
    if(!anchor)return;
    showModal('<h3>Add Atom</h3>'+
        '<label>Element:</label><select id="m-elem">'+
        '<option>H</option><option>He</option>'+
        '<option>Li</option><option>Be</option><option>B</option><option>C</option><option>N</option><option>O</option><option>F</option><option>Ne</option>'+
        '<option>Na</option><option>Mg</option><option>Al</option><option>Si</option><option>P</option><option>S</option><option>Cl</option><option>Ar</option>'+
        '<option>K</option><option>Ca</option><option>Sc</option><option>Ti</option><option>V</option><option>Cr</option><option>Mn</option><option>Fe</option>'+
        '<option>Co</option><option>Ni</option><option>Cu</option><option>Zn</option><option>Ga</option><option>Ge</option><option>As</option><option>Se</option>'+
        '<option>Br</option><option>Kr</option><option>Rb</option><option>Sr</option><option>Y</option><option>Zr</option><option>Nb</option><option>Mo</option>'+
        '<option>Ru</option><option>Rh</option><option>Pd</option><option>Ag</option><option>Cd</option><option>In</option><option>Sn</option><option>Sb</option>'+
        '<option>Te</option><option>I</option><option>Xe</option><option>Cs</option><option>Ba</option><option>La</option><option>Ce</option><option>Pr</option>'+
        '<option>Nd</option><option>Sm</option><option>Eu</option><option>Gd</option><option>Tb</option><option>Dy</option><option>Ho</option><option>Er</option>'+
        '<option>Tm</option><option>Yb</option><option>Lu</option><option>Hf</option><option>Ta</option><option>W</option><option>Re</option><option>Os</option>'+
        '<option>Ir</option><option>Pt</option><option>Au</option><option>Hg</option><option>Tl</option><option>Pb</option><option>Bi</option>'+
        '</select>'+
        '<label>Bond order:</label><select id="m-bond-order"><option value="1">Single (1.0)</option><option value="1.5">Aromatic (1.5)</option><option value="2">Double (2.0)</option><option value="3">Triple (3.0)</option></select>'+
        '<label>Bond length (A):</label><input type="number" id="m-val" value="1.09" step="0.01" min="0.3" max="5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    document.getElementById('m-elem').addEventListener('change',function(){
        var defaults={H:1.09,He:1.30,Li:2.12,Be:1.67,B:1.58,C:1.54,N:1.47,O:1.43,F:1.36,Ne:1.35,
            Na:1.88,Mg:1.63,Al:1.84,Si:1.87,P:1.80,S:1.82,Cl:1.77,Ar:1.74,
            K:2.34,Ca:2.00,Sc:1.88,Ti:1.87,V:1.79,Cr:1.79,Mn:1.83,Fe:1.80,
            Co:1.78,Ni:1.73,Cu:1.84,Zn:1.88,Ga:1.87,Ge:1.88,As:1.87,Se:1.90,
            Br:1.94,Kr:1.90,Rb:2.53,Sr:2.15,Y:2.12,Zr:2.06,Nb:2.04,Mo:2.08,
            Ru:2.07,Rh:2.09,Pd:2.05,Ag:2.10,Cd:2.07,In:2.10,Sn:2.17,Sb:2.12,
            Te:2.14,I:2.14,Xe:2.16,Cs:2.65,Ba:2.22,La:2.32,Ce:2.30,Pr:2.31,
            Nd:2.30,Sm:2.29,Eu:2.29,Gd:2.29,Tb:2.28,Dy:2.28,Ho:2.27,Er:2.27,
            Tm:2.26,Yb:2.26,Lu:2.26,Hf:2.23,Ta:2.22,W:2.18,Re:2.20,Os:2.19,
            Ir:2.16,Pt:2.13,Au:2.14,Hg:2.14,Tl:2.20,Pb:2.22,Bi:2.23};
        document.getElementById('m-val').value=defaults[this.value]||1.5;
    });
    document.getElementById('m-ok').addEventListener('click',function(){
        pushUndo();
        provMark('addAtoms',1);
        var el=document.getElementById('m-elem').value;
        var bl=parseFloat(document.getElementById('m-val').value)||1.5;
        var bondOrder=parseFloat(document.getElementById('m-bond-order').value)||1;
        var dir={x:0,y:0,z:1};
        if(CRY){
            var baseAnchorIdx=anchor.baseIdx!=null?anchor.baseIdx:anchorIdx;
            var baseAnchor=CRY.baseAtoms[baseAnchorIdx];
            var bondedBase=[];
            CRY.baseBonds.forEach(function(b){
                if(b.atom1===baseAnchorIdx)bondedBase.push(b.atom2);
                if(b.atom2===baseAnchorIdx)bondedBase.push(b.atom1);
            });
            if(bondedBase.length>0){
                var nbrs=bondedBase.map(function(bi){
                    var ba=CRY.baseAtoms[bi];
                    var dx=ba.x-baseAnchor.x,dy=ba.y-baseAnchor.y,dz=ba.z-baseAnchor.z;
                    var dl=Math.sqrt(dx*dx+dy*dy+dz*dz);
                    if(dl<1e-10)return{x:0,y:0,z:1};
                    return{x:dx/dl,y:dy/dl,z:dz/dl};
                });
                dir=pickFreeDirection(nbrs);
            }
            var newBaseIdx=CRY.baseAtoms.length;
            CRY.baseAtoms.push({element:el,x:baseAnchor.x+dir.x*bl,y:baseAnchor.y+dir.y*bl,z:baseAnchor.z+dir.z*bl,index:newBaseIdx,baseIdx:newBaseIdx,occupancy:1});
            CRY.baseBonds.push({atom1:baseAnchorIdx,atom2:newBaseIdx,order:bondOrder});
            rebuildCrystal();rebuildScene();hideModal();originalCoords=null;resetSelection();
        }else{
            var bonded=[];
            MD.bonds.forEach(function(b){
                if(b.atom1===anchorIdx)bonded.push(b.atom2);
                if(b.atom2===anchorIdx)bonded.push(b.atom1);
            });
            if(bonded.length>0){
                var nbrs=bonded.map(function(bi){
                    var ba=MD.atoms[bi];
                    var dx=ba.x-anchor.x,dy=ba.y-anchor.y,dz=ba.z-anchor.z;
                    var dl=Math.sqrt(dx*dx+dy*dy+dz*dz);
                    if(dl<1e-10)return{x:0,y:0,z:1};
                    return{x:dx/dl,y:dy/dl,z:dz/dl};
                });
                dir=pickFreeDirection(nbrs);
            }
            var newIdx=MD.atoms.length;
            MD.atoms.push({element:el,x:anchor.x+dir.x*bl,y:anchor.y+dir.y*bl,z:anchor.z+dir.z*bl,color:MD.atomColors[el]||'#FF1493',index:newIdx});
            MD.bonds.push({atom1:anchorIdx,atom2:newIdx,order:bondOrder});
            rebuildScene();hideModal();originalCoords=null;resetSelection();
        }
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
}

function showReplaceAtomModal(){
    if(selectedAtoms.length===0)return;
    var names=formatAtomList(selectedAtoms);
    showModal('<h3>Replace Atom</h3>'+
        '<div class="current-val">Replace '+selectedAtoms.length+' atom(s): '+names+'</div>'+
        '<label>New element:</label><select id="m-elem">'+
        '<option>H</option><option>He</option>'+
        '<option>Li</option><option>Be</option><option>B</option><option>C</option><option>N</option><option>O</option><option>F</option><option>Ne</option>'+
        '<option>Na</option><option>Mg</option><option>Al</option><option>Si</option><option>P</option><option>S</option><option>Cl</option><option>Ar</option>'+
        '<option>K</option><option>Ca</option><option>Sc</option><option>Ti</option><option>V</option><option>Cr</option><option>Mn</option><option>Fe</option>'+
        '<option>Co</option><option>Ni</option><option>Cu</option><option>Zn</option><option>Ga</option><option>Ge</option><option>As</option><option>Se</option>'+
        '<option>Br</option><option>Kr</option><option>Rb</option><option>Sr</option><option>Y</option><option>Zr</option><option>Nb</option><option>Mo</option>'+
        '<option>Ru</option><option>Rh</option><option>Pd</option><option>Ag</option><option>Cd</option><option>In</option><option>Sn</option><option>Sb</option>'+
        '<option>Te</option><option>I</option><option>Xe</option><option>Cs</option><option>Ba</option><option>La</option><option>Ce</option><option>Pr</option>'+
        '<option>Nd</option><option>Sm</option><option>Eu</option><option>Gd</option><option>Tb</option><option>Dy</option><option>Ho</option><option>Er</option>'+
        '<option>Tm</option><option>Yb</option><option>Lu</option><option>Hf</option><option>Ta</option><option>W</option><option>Re</option><option>Os</option>'+
        '<option>Ir</option><option>Pt</option><option>Au</option><option>Hg</option><option>Tl</option><option>Pb</option><option>Bi</option>'+
        '</select>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Replace</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        pushUndo();
        provMark('replaceAtoms',selectedAtoms.length);
        var el=document.getElementById('m-elem').value;
        selectedAtoms.forEach(function(idx){
            var a=MD.atoms[idx];
            a.element=el;
            a.color=MD.atomColors[el]||'#FF1493';
            if(CRY){
                var bi=a.baseIdx!=null?a.baseIdx:idx;
                if(CRY.baseAtoms[bi]){
                    CRY.baseAtoms[bi].element=el;
                    CRY.baseAtoms[bi].color=a.color;
                }
            }
        });
        if(CRY){rebuildCrystal();}
        rebuildScene();hideModal();originalCoords=null;resetSelection();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
}

function showDeleteAtomModal(){
    var idx=selectedAtoms[0];
    var a=MD.atoms[idx];
    if(!a)return;
    var name=a.element+(idx+1);
    showModal('<h3>Delete Atom</h3>'+
        '<div class="current-val">Delete '+name+'?</div>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok mbtn-danger" id="m-ok">Delete</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        pushUndo();
        provMark('delAtoms',1);
        if(CRY){
            var baseIdx=a.baseIdx!=null?a.baseIdx:idx;
            CRY.baseAtoms.splice(baseIdx,1);
            CRY.baseAtoms.forEach(function(ba,i){ba.index=i;ba.baseIdx=i});
            CRY.baseBonds=CRY.baseBonds.filter(function(b){return b.atom1!==baseIdx&&b.atom2!==baseIdx}).map(function(b){
                return{atom1:b.atom1>baseIdx?b.atom1-1:b.atom1,atom2:b.atom2>baseIdx?b.atom2-1:b.atom2,order:b.order,shift:b.shift};
            });
            rebuildCrystal();rebuildScene();hideModal();resetSelection();
        }else{
            MD.atoms.splice(idx,1);
            MD.atoms.forEach(function(a,i){a.index=i});
            MD.bonds=MD.bonds.filter(function(b){return b.atom1!==idx&&b.atom2!==idx}).map(function(b){
                return{atom1:b.atom1>idx?b.atom1-1:b.atom1,atom2:b.atom2>idx?b.atom2-1:b.atom2,order:b.order};
            });
            rebuildScene();hideModal();resetSelection();
        }
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();resetSelection()});
}

function showBatchDeleteModal(){
    if(selectedAtoms.length===0)return;
    stopMoveSession();
    var names=formatAtomList(selectedAtoms);
    showModal('<h3>Delete Atoms</h3>'+
        '<div class="current-val">Delete '+selectedAtoms.length+' atoms: '+names+'?</div>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok mbtn-danger" id="m-ok">Delete</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){pushUndo();provMark('delAtoms',selectedAtoms.length);batchDeleteSelected()});
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();resetSelection()});
}

function batchDeleteSelected(){
    // Atom indices become invalid after deletion — close any active rotate
    // group session or move session so stale state can't be reused.
    stopMoveSession();
    if(rotActive){resetRotAxisState();originalCoords=null}
    var delSet={};
    selectedAtoms.forEach(function(idx){delSet[idx]=true});
    var sortedAsc=selectedAtoms.slice().sort(function(a,b){return a-b});
    if(CRY){
        var baseToDelete={};
        selectedAtoms.forEach(function(idx){
            var a=MD.atoms[idx];
            var bi=a.baseIdx!=null?a.baseIdx:idx;
            baseToDelete[bi]=true;
        });
        var newBaseAtoms=[];
        var oldToNew={};
        CRY.baseAtoms.forEach(function(ba,i){
            if(baseToDelete[i])return;
            var newIdx=newBaseAtoms.length;
            oldToNew[i]=newIdx;
            var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
            newBaseAtoms.push({element:ba.element,x:ba.x,y:ba.y,z:ba.z,
                index:newIdx,occupancy:ba.occupancy!=null?ba.occupancy:1,baseIdx:newIdx,color:color});
        });
        CRY.baseBonds=CRY.baseBonds.filter(function(b){
            return !baseToDelete[b.atom1]&&!baseToDelete[b.atom2];
        }).map(function(b){
            return{atom1:oldToNew[b.atom1],atom2:oldToNew[b.atom2],order:b.order,shift:b.shift};
        });
        CRY.baseAtoms=newBaseAtoms;
        rebuildCrystal();rebuildScene();hideModal();resetSelection();
    }else{
        for(var k=sortedAsc.length-1;k>=0;k--){
            MD.atoms.splice(sortedAsc[k],1);
        }
        MD.atoms.forEach(function(a,i){a.index=i});
        MD.bonds=MD.bonds.filter(function(b){
            return !delSet[b.atom1]&&!delSet[b.atom2];
        }).map(function(b){
            var newA=b.atom1,newB=b.atom2;
            sortedAsc.forEach(function(d){
                if(d<b.atom1)newA--;
                if(d<b.atom2)newB--;
            });
            return{atom1:newA,atom2:newB,order:b.order};
        });
        rebuildScene();hideModal();resetSelection();
    }
}

// ===== Save provenance =====
// A saved file's title line states where the structure came from and what has
// been done to it here: the file it was opened from, the frame of a trajectory
// it was taken out of, and the edits that changed it — atoms added or deleted,
// the structure translated or turned, an optimisation, a displacement along a
// normal mode — with the numbers a reader would otherwise have to reconstruct.
//
// The record is written where the change happens, not inferred at save time,
// because "was this structure edited?" cannot be answered from the coordinates
// alone, and it travels in the undo snapshot, so undoing an edit also takes its
// clause off the title. Only formats that own a title or comment line get it;
// every other format is written exactly as before.
// The state itself is declared with the rest of the view state, above
// pushUndo, because the undo snapshot is the first thing that reads it; the
// functions that fill it in are the ones below.

// Basename without a regex: the separator escape is the one thing in this page
// that is easy to get wrong when the code itself lives in a template literal.
function provBaseName(p){
    var s=String(p||'');
    var cut=-1;
    for(var i=0;i<s.length;i++){
        var c=s.charAt(i);
        if(c==='/'||c===String.fromCharCode(92))cut=i;
    }
    return s.substring(cut+1);
}
function provClone(){
    var edits=[],modes={},i,k;
    for(i=0;i<provEdits.length;i++)edits.push({k:provEdits[i].k,n:provEdits[i].n});
    for(k in provModes){if(provModes.hasOwnProperty(k))modes[k]={freq:provModes[k].freq,scale:provModes[k].scale}}
    return {edits:edits,merged:provMerged.slice(),modes:modes};
}
function provRestore(p){
    var i,k;
    provEdits=[];provMerged=[];provModes={};
    if(!p)return;
    if(p.edits)for(i=0;i<p.edits.length;i++)provEdits.push({k:p.edits[i].k,n:p.edits[i].n});
    if(p.merged)provMerged=p.merged.slice();
    if(p.modes)for(k in p.modes){if(p.modes.hasOwnProperty(k))provModes[k]={freq:p.modes[k].freq,scale:p.modes[k].scale}}
}
/** The structure itself was replaced (another frame, another molecule): the
 *  edits and the merged files belong to the structure that is gone. */
function provReset(){provEdits=[];provMerged=[];provModes={}}
/** Record one kind of change. Counted kinds add up, so ten single deletions
 *  read as "10 atoms deleted" rather than ten clauses. */
function provMark(kind,n){
    for(var i=0;i<provEdits.length;i++){
        if(provEdits[i].k===kind){provEdits[i].n=(provEdits[i].n||0)+(n||1);return}
    }
    provEdits.push({k:kind,n:n||1});
}
/** A change that was taken back without an Undo (MM "Restore initial", a
 *  displacement dragged back to zero) is not an edit any more. */
function provUnmark(kind){
    for(var i=0;i<provEdits.length;i++){
        if(provEdits[i].k===kind){provEdits.splice(i,1);return}
    }
}
function provAddSource(name){
    if(!name)return;
    if(provMerged.indexOf(name)<0)provMerged.push(name);
}
function provFreqOf(idx){
    var m=(typeof NORMAL_MODES!=='undefined'&&NORMAL_MODES)?NORMAL_MODES[idx]:null;
    return (m&&typeof m.frequency==='number')?m.frequency:null;
}
/** Keep a displacement the slider committed: once the session closes it is no
 *  longer in vdComp, but it is still in the coordinates and must still be said. */
function provCommitModes(){
    if(typeof vdComp==='undefined'||!vdComp)return;
    for(var i=0;i<vdComp.length;i++){
        var c=vdComp[i];
        if(Math.abs(c.scale)<1e-9){delete provModes[c.idx];continue}
        provModes[c.idx]={freq:provFreqOf(c.idx),scale:c.scale};
    }
}
function provDropModes(){provModes={}}
function provModesText(){
    var seen={},k,i;
    for(k in provModes){if(provModes.hasOwnProperty(k))seen[k]={freq:provModes[k].freq,scale:provModes[k].scale}}
    // The live session is not committed yet, but its displacement is on screen.
    if(typeof vdActive!=='undefined'&&vdActive&&typeof vdComp!=='undefined'&&vdComp){
        for(i=0;i<vdComp.length;i++){
            var c=vdComp[i];
            if(Math.abs(c.scale)<1e-9)delete seen[c.idx];
            else seen[c.idx]={freq:provFreqOf(c.idx),scale:c.scale};
        }
    }
    var keys=[];
    for(k in seen){if(seen.hasOwnProperty(k))keys.push(k)}
    keys.sort(function(a,b){return a-b});
    var parts=[];
    for(i=0;i<keys.length;i++){
        var e=seen[keys[i]];
        var txt='mode '+(Number(keys[i])+1);
        if(e.freq!=null){
            txt+=' ('+Math.abs(e.freq).toFixed(1)+(e.freq<0?'i':'')+' cm-1)';
        }
        txt+=' '+(e.scale>=0?'+':'-')+Math.abs(e.scale).toFixed(2)+' A';
        parts.push(txt);
    }
    return parts.join(' + ');
}
/** The three geometry dialogs preview live and keep the result on OK; nothing
 *  moved means nothing to report, however many times the slider was touched. */
function provMarkIfMoved(kind){
    if(!originalCoords)return;
    for(var i=0;i<originalCoords.length;i++){
        var oc=originalCoords[i],a=MD.atoms[i];
        if(!oc||!a)continue;
        if(Math.abs(a.x-oc.x)>1e-9||Math.abs(a.y-oc.y)>1e-9||Math.abs(a.z-oc.z)>1e-9){provMark(kind,1);return}
    }
}
/** The supercell is a property of the crystal on screen, not an edit that can be
 *  undone, so it is read from the bounds at save time instead of being marked. */
function provSupercellText(){
    if(typeof CRY==='undefined'||!CRY||typeof boundary==='undefined'||!boundary)return '';
    var ax=[];
    if(boundary.aMin!==0||boundary.aMax!==1)ax.push('a['+boundary.aMin+','+boundary.aMax+']');
    if(boundary.bMin!==0||boundary.bMax!==1)ax.push('b['+boundary.bMin+','+boundary.bMax+']');
    if(boundary.cMin!==0||boundary.cMax!==1)ax.push('c['+boundary.cMin+','+boundary.cMax+']');
    return ax.length?('supercell '+ax.join(' ')):'';
}
/** The title as PDB TITLE records: 60 characters per record, with PDB's
 *  continuation number in columns 9-10 rather than one over-long line. */
function provPdbTitle(text){
    var s=String(text||''),out='',serial=1,i=0;
    while(i<s.length){
        var head=(serial<10)?('TITLE   '+serial+' '):('TITLE  '+serial+' ');
        out+=head+s.substring(i,i+60)+'\\n';
        i+=60;serial++;
    }
    return out;
}
function provEditText(e){
    var n=e.n||1;
    switch(e.k){
        case 'addAtoms':return n===1?'1 atom added':n+' atoms added';
        case 'delAtoms':return n===1?'1 atom deleted':n+' atoms deleted';
        case 'cutAtoms':return n===1?'1 atom cut':n+' atoms cut';
        case 'pasteAtoms':return n===1?'1 atom pasted':n+' atoms pasted';
        case 'replaceAtoms':return n===1?'1 atom replaced':n+' atoms replaced';
        case 'addH':return n===1?'1 H added':n+' H added';
        case 'translate':return 'translated';
        case 'rotate':return 'rotated';
        case 'mmOpt':return 'MM-optimized';
        case 'bondLen':return 'bond length adjusted';
        case 'bondAngle':return 'bond angle adjusted';
        case 'dihedral':return 'dihedral adjusted';
        case 'disorder':return 'disorder removed';
        case 'supercell':return 'supercell rebuilt';
    }
    return '';
}
/** Join as many whole clauses as fit, never cutting one in half. */
function provJoinText(head,parts,max){
    if(!parts.length)return head;
    var whole=(head?head+'; ':'')+parts.join(', ');
    if(whole.length<=max)return whole;
    var used=[];
    for(var i=0;i<parts.length;i++){
        var cand=used.concat([parts[i]]);
        var txt=(head?head+'; ':'')+cand.join(', ');
        if(txt.length>max-6){
            if(used.length)return (head?head+'; ':'')+used.join(', ')+', ...';
            return txt.substring(0,max-4)+' ...';
        }
        used=cand;
    }
    return (head?head+'; ':'')+used.join(', ');
}
/** The file's own title, when it says something the provenance does not. */
function provOriginalTitle(){
    if(typeof MD==='undefined'||!MD)return '';
    // A trajectory's title belongs to its *first* structure — a log's title is
    // the comment of the geometry the file opens on, and switching frames does
    // not rewrite it. The frame clause already names the structure on screen, so
    // appending this would state the opposite of it ("frame 79 of 79 ... | Step
    // 1"): the tail is dropped rather than contradicting the line.
    if(typeof totalFrames!=='undefined'&&totalFrames>1)return '';
    // Same for a script's molecule list: the title describes the first molecule,
    // and the molecule clause names the one on screen.
    if(typeof VMD!=='undefined'&&VMD&&VMD.molecules&&VMD.molecules.length>1&&vmdMolIdx>0)return '';
    var t=String(MD.title||'');
    t=t.split('\\n').join(' ').split('\\r').join(' ').split('\\t').join(' ');
    t=t.trim();
    if(t.indexOf(PROV_PREFIX)===0)return '';
    if(t===''||t==='Modified structure'||t==='Empty'||t==='No structures found')return '';
    return t;
}
function provTitle(fallback,keepOriginal){
    var merged=provMerged.slice();
    var src=provSource;
    if(!src&&merged.length)src=merged.shift();
    var head=[];
    if(src)head.push('from '+src);
    if(totalFrames>1){
        var fr='frame '+(currentFrame+1)+' of '+totalFrames;
        var fl=(MD.frames&&MD.frames[currentFrame]&&MD.frames[currentFrame].stepLabel)?String(MD.frames[currentFrame].stepLabel):'';
        if(fl)fr+=', '+fl;
        head.push(fr);
    }else if(VMD&&VMD.molecules&&VMD.molecules.length>1){
        head.push('molecule '+(vmdMolIdx+1)+' of '+VMD.molecules.length);
    }
    // Merged files come after the frame: the frame belongs to the source file.
    if(merged.length)head.push('merged with '+merged.join(', '));
    var body=[];
    for(var i=0;i<provEdits.length;i++){
        var t=provEditText(provEdits[i]);
        if(t)body.push(t);
    }
    var modeTxt=provModesText();
    if(modeTxt)body.push('displaced along '+modeTxt);
    var scTxt=provSupercellText();
    if(scTxt)body.push(scTxt);
    var line=provJoinText(head.join(', '),body,PROV_MAX);
    if(keepOriginal!==false){
        var orig=provOriginalTitle();
        if(orig){
            if(!line){
                // Nothing to add provenance to: the file's own title stands in.
                line=orig.length>PROV_MAX?(orig.substring(0,PROV_MAX-3)+'...'):orig;
            }else{
                var room=PROV_MAX-line.length-3;
                if(room>=14)line+=' | '+(orig.length>room?orig.substring(0,room-3)+'...':orig);
            }
        }
    }
    if(!line)line=fallback||'Modified structure';
    return PROV_PREFIX+line;
}

function doSave(){
    if(!MD.atoms||MD.atoms.length===0){alert('No atoms to save');return}
    for(var ci=0;ci<MD.atoms.length;ci++){var ca=MD.atoms[ci];if(isNaN(ca.x)||isNaN(ca.y)||isNaN(ca.z)){alert('Atom '+(ci+1)+' has invalid coordinates (NaN). Cannot save.');return}}
    if(CRY&&!CRY_INV){alert('Crystal lattice is degenerate (zero determinant). Cannot save crystal formats.');return}
    try{
    // The line every format below is titled with: where this structure came
    // from and what has been done to it since (see "Save provenance" above).
    // A GJF title section can hold several lines, so it takes the provenance
    // line and then the file's original title lines verbatim; every other
    // format gets one line, with the original title kept after a '|'.
    var saveTitle=provTitle(MD.title||'Modified structure',true);
    var gjfTitle=provTitle(MD.title||'Modified structure',false);
    var xyz=MD.atoms.length+'\\n'+saveTitle+'\\n';
    MD.atoms.forEach(function(a){xyz+=a.element+'  '+a.x.toFixed(6)+'  '+a.y.toFixed(6)+'  '+a.z.toFixed(6)+'\\n'});
    var gjf='';
    var meta=MD.gjfMeta;
    if(meta){
        meta.link0Lines.forEach(function(l){gjf+=l+'\\n'});
        gjf+=meta.routeLine+'\\n\\n';
        gjf+=gjfTitle+'\\n';
        // The file's own title card is kept, minus any line an earlier save of
        // ours wrote there: the provenance is rebuilt, never accumulated.
        meta.titleLines.forEach(function(l){
            var t=String(l);
            if(t.indexOf(PROV_PREFIX)===0)return;
            if(t.trim()==='')return;
            gjf+=t+'\\n';
        });
        gjf+='\\n'+meta.chargeMultLine+'\\n';
    }else{
        gjf='%chk=molecule.chk\\n%mem=4GB\\n%nproc=4\\n# B3LYP/6-31G(d)\\n\\n'+saveTitle+'\\n\\n0 1\\n';
    }
    MD.atoms.forEach(function(a){
        var el=a.element; while(el.length<2) el+=' ';
        var sx=a.x.toFixed(8); while(sx.length<17) sx=' '+sx;
        var sy=a.y.toFixed(8); while(sy.length<17) sy=' '+sy;
        var sz=a.z.toFixed(8); while(sz.length<17) sz=' '+sz;
        gjf+=' '+el+' '+sx+' '+sy+' '+sz+'\\n';
    });
    gjf+='\\n';
    MD.atoms.forEach(function(a,i){
        var parts=[i+1];
        MD.bonds.forEach(function(b){
            if(b.atom1===i)parts.push(b.atom2+1,b.order.toFixed(1));
        });
        gjf+=parts.join(' ')+'\\n';
    });
    gjf+='\\n';
    if(meta&&meta.afterConnectContent){
        var accLines=meta.afterConnectContent.split('\\n').filter(function(ln){
            var t=ln.trim();
            if(t===''||t.startsWith('--'))return true;
            var p=t.split(/\\s+/);
            var firstNum=parseInt(p[0],10);
            if(!isNaN(firstNum)&&firstNum>=1&&firstNum<=MD.atoms.length&&p.every(function(pp){return !isNaN(parseFloat(pp))})){
                return false;
            }
            return true;
        });
        if(accLines.length>0)gjf+=accLines.join('\\n')+'\\n';
    }

    var chrg=MD.charge||0;
    var mult=MD.multiplicity||1;
    var AN2={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Ne:10,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Sc:21,Ti:22,V:23,Cr:24,Mn:25,Fe:26,Co:27,Ni:28,Cu:29,Zn:30,Ga:31,Ge:32,As:33,Se:34,Br:35,Kr:36,Rb:37,Sr:38,Y:39,Zr:40,Nb:41,Mo:42,Tc:43,Ru:44,Rh:45,Pd:46,Ag:47,Cd:48,In:49,Sn:50,Sb:51,Te:52,I:53,Xe:54,Cs:55,Ba:56,La:57,Ce:58,Pr:59,Nd:60,Pm:61,Sm:62,Eu:63,Gd:64,Tb:65,Dy:66,Ho:67,Er:68,Tm:69,Yb:70,Lu:71,Hf:72,Ta:73,W:74,Re:75,Os:76,Ir:77,Pt:78,Au:79,Hg:80,Tl:81,Pb:82,Bi:83,Po:84,At:85,Rn:86};

    var coord='$coord\\n';
    var ANG_TO_BOHR=1.8897259886;
    MD.atoms.forEach(function(a){coord+='  '+(a.x*ANG_TO_BOHR).toFixed(8)+'  '+(a.y*ANG_TO_BOHR).toFixed(8)+'  '+(a.z*ANG_TO_BOHR).toFixed(8)+' '+a.element.toLowerCase()+'\\n'});
    coord+='$end\\n';
    if(chrg!==0)coord+='$chrg '+chrg+'\\n';
    if(mult!==1)coord+='$spin '+(mult-1)/2+'\\n';

    // ORCA's first line is a comment line, so the provenance is one there and
    // the keywords keep their own line below it.
    var orcaInp='# '+saveTitle+'\\n! B3LYP def2-SVP\\n\\n* xyz '+chrg+' '+mult+'\\n';
    MD.atoms.forEach(function(a){orcaInp+=a.element+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+'\\n'});
    orcaInp+='*\\n';

    var mol2='@<TRIPOS>MOLECULE\\n'+saveTitle+'\\n'+MD.atoms.length+' '+MD.bonds.length+' 0 0 0\\nSMALL\\nNO_CHARGES\\n\\n';
    mol2+='@<TRIPOS>ATOM\\n';
    MD.atoms.forEach(function(a,i){mol2+=(i+1)+' '+a.element+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+' '+a.element+' 1 UNK 0.000\\n'});
    mol2+='@<TRIPOS>BOND\\n';
    var bondIdx=1;
    MD.bonds.forEach(function(b){
        var bt='1';
        if(b.order>=2.5)bt='3';
        else if(b.order>=1.75)bt='2';
        else if(b.order>=1.25)bt='ar';
        mol2+=bondIdx+' '+(b.atom1+1)+' '+(b.atom2+1)+' '+bt+'\\n';
        bondIdx++;
    });

    var mol=saveTitle+'\\n '+MD.atoms.length+'  '+MD.bonds.length+'  0  0  0  0  0  0  0  0999 V2000\\n';
    MD.atoms.forEach(function(a){
        var z=AN2[a.element]||0;
        var sx=(a.x*10).toFixed(4);
        var sy=(a.y*10).toFixed(4);
        var sz=(a.z*10).toFixed(4);
        while(sx.length<10)sx=' '+sx;
        while(sy.length<10)sy=' '+sy;
        while(sz.length<10)sz=' '+sz;
        mol+=sx+sy+sz+' '+a.element+' 0  0  0  0  0  0  0  0  0  0  0  0\\n';
    });
    MD.bonds.forEach(function(b){
        var bt=1;
        if(b.order>=2.5)bt=3;
        else if(b.order>=1.75)bt=2;
        else if(b.order>=1.25)bt=4;
        var a1=b.atom1+1,a2=b.atom2+1;
        var s1=''+a1,s2=''+a2,s3=''+bt;
        while(s1.length<3)s1=' '+s1;
        while(s2.length<3)s2=' '+s2;
        while(s3.length<3)s3=' '+s3;
        mol+=s1+s2+s3+'  0  0  0  0\\n';
    });
    mol+='M  END\\n';

    // PDB owns a TITLE record, so the provenance goes in there (wrapped at 60
    // characters per record with a continuation number) instead of being lost.
    var pdb=provPdbTitle(saveTitle);
    MD.atoms.forEach(function(a,i){
        var serial=(i+1).toString();while(serial.length<5)serial=' '+serial;
        var name;
        if(a.element.length===1){name=' '+a.element+'  '}
        else{name=a.element+'  '}
        var sx=a.x.toFixed(3);while(sx.length<8)sx=' '+sx;
        var sy=a.y.toFixed(3);while(sy.length<8)sy=' '+sy;
        var sz=a.z.toFixed(3);while(sz.length<8)sz=' '+sz;
        var el=a.element.charAt(0).toUpperCase();
        if(a.element.length>1)el+=a.element.charAt(1).toLowerCase();
        while(el.length<2)el=' '+el;
        pdb+='ATOM  '+serial+name+'UNK A   1    '+sx+sy+sz+'  1.00  0.00          '+el+'\\n';
    });
    MD.bonds.forEach(function(b){
        var a1=(b.atom1+1).toString();while(a1.length<5)a1=' '+a1;
        var a2=(b.atom2+1).toString();while(a2.length<5)a2=' '+a2;
        pdb+='CONECT'+a1+a2+'\\n';
    });
    pdb+='END\\n';

    var mopac=saveTitle+'\\n';
    var chrgStr=chrg!==0?' CHARGE='+chrg:'';
    var multStr=mult!==1?' MS='+mult:'';
    mopac+='PM7'+chrgStr+multStr+'\\n\\n';
    MD.atoms.forEach(function(a){
        var sx=a.x.toFixed(5);var sy=a.y.toFixed(5);var sz=a.z.toFixed(5);
        mopac+=a.element+' '+sx+' 1 '+sy+' 1 '+sz+' 1\\n';
    });
    mopac+='\\n';

    var cifContent='';
    if(CRY){
        var inv=CRY_INV;
        cifContent='data_modified\\n';
        cifContent+='_cell_length_a '+CRY.a.toFixed(6)+'\\n';
        cifContent+='_cell_length_b '+CRY.b.toFixed(6)+'\\n';
        cifContent+='_cell_length_c '+CRY.c.toFixed(6)+'\\n';
        cifContent+='_cell_angle_alpha '+CRY.alpha.toFixed(6)+'\\n';
        cifContent+='_cell_angle_beta '+CRY.beta.toFixed(6)+'\\n';
        cifContent+='_cell_angle_gamma '+CRY.gamma.toFixed(6)+'\\n';
        if(CRY.spaceGroup)cifContent+='_symmetry_space_group_name_H-M \\''+CRY.spaceGroup+'\\'\\n';
        cifContent+='loop_\\n_symmetry_equiv_pos_as_xyz\\n';
        (CRY.symmetryOps||['x, y, z']).forEach(function(op){cifContent+='\\''+op+'\\'\\n'});
        cifContent+='loop_\\n_atom_site_label\\n_atom_site_type_symbol\\n_atom_site_fract_x\\n_atom_site_fract_y\\n_atom_site_fract_z\\n_atom_site_occupancy\\n';
        CRY.baseAtoms.forEach(function(ba,i){
            var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
            var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
            var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
            var label=ba.element+(i+1);
            var occ=(ba.occupancy!=null?ba.occupancy:1).toFixed(3);
            cifContent+=label+' '+ba.element+' '+fx.toFixed(6)+' '+fy.toFixed(6)+' '+fz.toFixed(6)+' '+occ+'\\n';
        });
    }

    var vaspContent='';
    if(CRY){
        var lv=CRY.latticeVectors;
        var vInv=CRY_INV;
        var elemOrder=[];
        var elemCounts={};
        MD.atoms.forEach(function(a){
            if(elemCounts[a.element]===undefined){elemOrder.push(a.element);elemCounts[a.element]=0}
            elemCounts[a.element]++;
        });
        vaspContent=saveTitle+'\\n';
        vaspContent+='1.0\\n';
        for(var ii=0;ii<3;ii++){
            vaspContent+='  '+lv[ii][0].toFixed(8)+'  '+lv[ii][1].toFixed(8)+'  '+lv[ii][2].toFixed(8)+'\\n';
        }
        vaspContent+=elemOrder.join('  ')+'\\n';
        vaspContent+=elemOrder.map(function(e){return elemCounts[e]}).join('  ')+'\\n';
        vaspContent+='Direct\\n';
        MD.atoms.forEach(function(a){
            var fx=vInv[0][0]*a.x+vInv[0][1]*a.y+vInv[0][2]*a.z;
            var fy=vInv[1][0]*a.x+vInv[1][1]*a.y+vInv[1][2]*a.z;
            var fz=vInv[2][0]*a.x+vInv[2][1]*a.y+vInv[2][2]*a.z;
            vaspContent+='  '+fx.toFixed(8)+'  '+fy.toFixed(8)+'  '+fz.toFixed(8)+'\\n';
        });
    }

    var cubeContent='';
    if(CRY){
        var clv=CRY.latticeVectors;
        var nGrid=2;
        cubeContent='Molecular Viewer export\\n';
        // The cube format has exactly two comment lines before the atom count:
        // the provenance takes the second one, it does not add a third.
        cubeContent+=saveTitle+'\\n';
        cubeContent+='-'+MD.atoms.length+'  0.000000  0.000000  0.000000\\n';
        cubeContent+=nGrid+'  '+(clv[0][0]/nGrid).toFixed(8)+'  '+(clv[0][1]/nGrid).toFixed(8)+'  '+(clv[0][2]/nGrid).toFixed(8)+'\\n';
        cubeContent+=nGrid+'  '+(clv[1][0]/nGrid).toFixed(8)+'  '+(clv[1][1]/nGrid).toFixed(8)+'  '+(clv[1][2]/nGrid).toFixed(8)+'\\n';
        cubeContent+=nGrid+'  '+(clv[2][0]/nGrid).toFixed(8)+'  '+(clv[2][1]/nGrid).toFixed(8)+'  '+(clv[2][2]/nGrid).toFixed(8)+'\\n';
        MD.atoms.forEach(function(a){
            var z=AN2[a.element]||0;
            cubeContent+=z+'  0.000000  '+a.x.toFixed(8)+'  '+a.y.toFixed(8)+'  '+a.z.toFixed(8)+'\\n';
        });
        var nVals=nGrid*nGrid*nGrid;
        for(var vi=0;vi<nVals;vi++){
            cubeContent+='0.0000e+00';
            if((vi+1)%6===0||vi===nVals-1){cubeContent+='\\n'}else{cubeContent+='  '}
        }
    }

    var vestaContent='';
    if(CRY){
        var vInv=CRY_INV;
        vestaContent='#VESTA_FORMAT_VERSION 3.5.4\\n\\n\\nCRYSTAL\\n\\nTITLE\\n'+saveTitle+'\\n\\nGROUP\\n    1   1 P 1\\n\\nSYMOP\\n 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1\\n-1.0 -1.0 -1.0  1 0 0  0 0 0  0 0 0\\n\\nTRANM 0\\n 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1\\n\\nLTRANSL\\n -1\\n 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\\n\\nLORIENT\\n -1   0   0   0   0\\n 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000\\n 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000\\n\\nLMATRIX\\n 1.000000  0.000000  0.000000  0.000000\\n 0.000000  1.000000  0.000000  0.000000\\n 0.000000  0.000000  1.000000  0.000000\\n 0.000000  0.000000  0.000000  1.000000\\n 0.000000  0.000000  0.000000\\n\\nCELLP\\n';
        vestaContent+=' '+CRY.a.toFixed(6)+'  '+CRY.b.toFixed(6)+'  '+CRY.c.toFixed(6)+'  '+CRY.alpha.toFixed(6)+'  '+CRY.beta.toFixed(6)+'  '+CRY.gamma.toFixed(6)+'\\n';
        vestaContent+=' 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\\n\\nSTRUC\\n';
        CRY.baseAtoms.forEach(function(ba,i){
            var fx=vInv[0][0]*ba.x+vInv[0][1]*ba.y+vInv[0][2]*ba.z;
            var fy=vInv[1][0]*ba.x+vInv[1][1]*ba.y+vInv[1][2]*ba.z;
            var fz=vInv[2][0]*ba.x+vInv[2][1]*ba.y+vInv[2][2]*ba.z;
            var label=ba.element+(i+1);
            var occ=(ba.occupancy!=null?ba.occupancy:1).toFixed(4);
            var idxStr=''+(i+1);while(idxStr.length<3)idxStr=' '+idxStr;
            var elStr=ba.element;while(elStr.length<3)elStr=elStr+' ';
            var labelStr=label;while(labelStr.length<4)labelStr=labelStr+' ';
            vestaContent+=idxStr+' '+elStr+' '+labelStr+' '+occ+'   '+fx.toFixed(6)+'   '+fy.toFixed(6)+'   '+fz.toFixed(6)+'    2i       1\\n';
            vestaContent+='                            0.000000   0.000000   0.000000  0.00\\n';
        });
    }
    }catch(err){
        alert('Failed to generate save content: '+(err&&err.message?err.message:String(err)));
        return;
    }

    // Saving a displaced structure is exactly what the feature is for, so the
    // dialog says which modes the geometry is built from instead of leaving the
    // composition to be discovered in the file.
    var vdNote=vdActive
        ?'<div class="current-val">Coordinates include a manual displacement: '+
            vdCompsText()+'.</div>'
        :'';
    showModal('<h3>Save File</h3>'+vdNote+
        '<label>Format:</label><select id="m-fmt">'+
        '<option value="xyz">XYZ (.xyz)</option>'+
        '<option value="gjf">Gaussian Input (.gjf)</option>'+
        (CRY?'<option value="cif">CIF (.cif)</option>':'')+
        (CRY?'<option value="vasp">VASP POSCAR (.vasp)</option>':'')+
        (CRY?'<option value="cube">Gaussian Cube (.cube)</option>':'')+
        (CRY?'<option value="vesta">VESTA (.vesta)</option>':'')+
        '<option value="coord">Turbomole Coord (.coord)</option>'+
        '<option value="inp">ORCA Input (.inp)</option>'+
        '<option value="mol2">MOL2 (.mol2)</option>'+
        '<option value="mol">MDL Mol (.mol)</option>'+
        '<option value="pdb">PDB (.pdb)</option>'+
        '<option value="mop">MOPAC Input (.mop)</option>'+
        '</select>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Save</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        var fmt=document.getElementById('m-fmt').value;
        var content,ext;
        switch(fmt){
            case 'gjf':content=gjf;ext='.gjf';break;
            case 'cif':content=cifContent;ext='.cif';break;
            case 'vasp':content=vaspContent;ext='.vasp';break;
            case 'cube':content=cubeContent;ext='.cube';break;
            case 'vesta':content=vestaContent;ext='.vesta';break;
            case 'coord':content=coord;ext='.coord';break;
            case 'inp':content=orcaInp;ext='.inp';break;
            case 'mol2':content=mol2;ext='.mol2';break;
            case 'mol':content=mol;ext='.mol';break;
            case 'pdb':content=pdb;ext='.pdb';break;
            case 'mop':content=mopac;ext='.mop';break;
            default:content=xyz;ext='.xyz';
        }
        vscodeApi.postMessage({command:'saveFile',content:content,suggestedName:'molecule_modified'+ext,filePath:MD.filePath||''});
        hideModal();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal()});
}

function parseSelectInput(input){
    var result=[];
    var tokens=input.split(/[\\s,;]+/);
    tokens.forEach(function(tok){
        if(!tok)return;
        var rangeMatch=tok.match(/^(\\d+)-(\\d+)$/);
        if(rangeMatch){
            var start=parseInt(rangeMatch[1],10);
            var end=parseInt(rangeMatch[2],10);
            if(!isNaN(start)&&!isNaN(end)){
                for(var k=start;k<=end;k++){
                    var idx=k-1;
                    if(idx>=0&&idx<MD.atoms.length&&result.indexOf(idx)<0)result.push(idx);
                }
            }
            return;
        }
        var elNumMatch=tok.match(/^([A-Za-z]+)(\d+)$/);
        if(elNumMatch){
            var elName=elNumMatch[1].charAt(0).toUpperCase()+elNumMatch[1].slice(1).toLowerCase();
            var en=parseInt(elNumMatch[2],10);
            var ei=en-1;
            if(ei>=0&&ei<MD.atoms.length&&MD.atoms[ei].element===elName&&result.indexOf(ei)<0)result.push(ei);
            return;
        }
        var num=parseInt(tok,10);
        if(!isNaN(num)&&num>0){
            var idx2=num-1;
            if(idx2<MD.atoms.length&&result.indexOf(idx2)<0)result.push(idx2);
            return;
        }
        var el=tok.charAt(0).toUpperCase()+tok.slice(1).toLowerCase();
        MD.atoms.forEach(function(a,i){
            if(a.element===el&&result.indexOf(i)<0)result.push(i);
        });
    });
    return result;
}

function updateTransform(){
    if(diffMode){
        moleculeGroup.quaternion.copy(rotQuat);
        pivotGroup.position.set(panX,panY,0);
        if(diffMolGroup){
            diffMolGroup.quaternion.copy(diffRotQuat);
            diffPivot.position.set(diffPanX,diffPanY,0);
        }
    }else{
        moleculeGroup.quaternion.copy(rotQuat);
        pivotGroup.position.set(panX,panY,0);
    }
    needsRender=true;
}

var canvas=renderer.domElement;
var raycaster=new THREE.Raycaster();
var mouse=new THREE.Vector2();
var isBoxSelecting=false;var boxStart=null;

// ===== Auto-rotate: flick the molecule and it keeps turning =====
// The toggle lives in the View-mode panel (#view-panel, bottom-right, shown
// only while the View button's mode is active) rather than in the toolbar.
// Two inputs feed one angular velocity:
//   - a left-drag that rotates the view (or a one-finger touch drag) arms a
//     spin: on release the molecule is handed the angular velocity the drag
//     was released with (spinLaunch);
//   - an arrow key, held, spins the view up — the longer the hold, the faster
//     (keySpinDrive) — and releasing it hands over the speed it reached.
// Either way the molecule then keeps turning about a fixed world axis until
// something stops it. A right/middle press pans and leaves the spin running.
// This is a view-only rotation — it writes rotQuat (or the grabbed side's
// diffRotQuat in Diff mode) and never touches the coordinates, so it cannot
// fight any edit.
//
// Arbitration: exactly one writer holds a rotation quaternion at a time. A
// left press, Reset View, leaving Diff mode, leaving View mode and switching
// the toggle off each call stopSpin(); an arrow key first stops any spin and
// then drives its own, so the two can never compose with each other into a
// drift. Ctrl+arrow in Move Atoms rotates the selection, not the view, and
// never reaches any of this.
var spinOn=false;                          // toolbar toggle (off by default)
var spinAxis=new THREE.Vector3(0,1,0);     // fixed world axis of the current spin
var spinSpeed=0;                           // rad/s, 0 = not spinning
var spinQuatRef=null;                      // rotQuat or diffRotQuat, whichever spins
var spinSamples=[];                        // recent {t,x,y} of the arming drag
var spinPrevT=0;                           // last animate() timestamp, for the frame step
var SPIN_ANGLE_PER_PX=0.008;               // must equal the drag rotation scale (see mousemove)
var SPIN_WINDOW_MS=120;                    // only the last 120 ms before release set the speed
var SPIN_GAIN=1.0;                         // released speed = drag speed * SPIN_GAIN
var SPIN_MAX=8.0;                          // rad/s cap, so a wild flick cannot blur the view
var SPIN_MIN=0.05;                         // below this a release reads as a stop, not a spin
var KEY_SPIN_BASE=0.5;                     // rad/s an arrow key starts at (Auto Rotate on)
var KEY_SPIN_RATE=2.5;                     // rad/s added per second the key is held, up to SPIN_MAX
var keySpinActive=false;                   // an arrow-key spin-up session is running
var keySpinHold=0;                         // seconds the current arrow-key hold has lasted

function stopSpin(){spinSpeed=0;spinQuatRef=null;spinSamples.length=0;keySpinActive=false;keySpinHold=0}

// One point per rotation mousemove/touchmove. The buffer keeps only the last
// SPIN_WINDOW_MS, so a drag that is held still before release releases at a
// standstill (no spin) rather than throwing whatever speed it reached earlier.
function spinRecordPoint(x,y){
    if(!spinOn)return;
    var t=performance.now();
    spinSamples.push({t:t,x:x,y:y});
    while(spinSamples.length>0&&t-spinSamples[0].t>SPIN_WINDOW_MS)spinSamples.shift();
    if(spinSamples.length>64)spinSamples.shift();
}

// Release: the accumulated drag becomes an angular velocity. Distance over
// time is the release speed, so dragging farther or flicking faster spins the
// molecule faster; the direction is the drag direction, mapped onto the world
// X/Y axes exactly the way the drag itself rotated the view. Every path that
// decides not to spin clears the spin itself, so a stale speed from an earlier
// flick can never leak through a release that should have been a standstill.
function spinLaunch(){
    if(!spinOn){stopSpin();return}
    var now=performance.now();
    while(spinSamples.length>0&&now-spinSamples[0].t>SPIN_WINDOW_MS)spinSamples.shift();
    if(spinSamples.length<2){stopSpin();return}
    var a=spinSamples[0],b=spinSamples[spinSamples.length-1];
    spinSamples.length=0;
    var dt=(b.t-a.t)/1000;
    if(!(dt>0.001)){stopSpin();return}
    var wx=(b.y-a.y)*SPIN_ANGLE_PER_PX*SPIN_GAIN/dt;
    var wy=(b.x-a.x)*SPIN_ANGLE_PER_PX*SPIN_GAIN/dt;
    var mag=Math.sqrt(wx*wx+wy*wy);
    if(!(mag>SPIN_MIN)){stopSpin();return}
    if(mag>SPIN_MAX){var s=SPIN_MAX/mag;wx*=s;wy*=s;mag=SPIN_MAX}
    spinAxis.set(wx/mag,wy/mag,0);
    spinSpeed=mag;
    spinQuatRef=diffMode?(diffTransformSide==='right'?diffRotQuat:rotQuat):rotQuat;
}

// A held arrow key, with Auto Rotate on, drives the very same angular velocity
// a released drag does: the axis is the one the pressed keys name, and the
// speed ramps up the longer they are held, so holding a key longer spins the
// view faster and releasing it hands that speed straight to the free rotation
// (nothing is applied here — the one application site in animate() turns
// whatever spinQuatRef/spinAxis/spinSpeed hold). The key mapping is the mouse
// mapping: Right/Left are +Y/-Y and Down/Up are +X/-X, so a flick and a key
// hold spin the same way. Ctrl+arrows in Move Atoms rotate the selection, not
// the view, and never reach here.
function keySpinDrive(dt){
    if(!keySpinActive){keySpinActive=true;keySpinHold=0;spinSpeed=0}
    if(dt>0)keySpinHold+=dt;
    var kx=0,ky=0;
    if(keyState.ArrowDown)kx+=1;
    if(keyState.ArrowUp)kx-=1;
    if(keyState.ArrowRight)ky+=1;
    if(keyState.ArrowLeft)ky-=1;
    if(kx===0&&ky===0){stopSpin();return}
    var m=Math.sqrt(kx*kx+ky*ky);
    spinAxis.set(kx/m,ky/m,0);
    spinSpeed=Math.min(KEY_SPIN_BASE+KEY_SPIN_RATE*keySpinHold,SPIN_MAX);
    spinQuatRef=diffMode?(diffTransformSide==='right'?diffRotQuat:rotQuat):rotQuat;
}

function updateSpinButton(){
    var b=document.getElementById('spin-btn');
    if(!b)return;
    b.classList.toggle('active',spinOn);
    b.textContent=spinOn?'On':'Off';
}
(function(){
    var b=document.getElementById('spin-btn');
    if(!b)return;
    b.addEventListener('click',function(){
        spinOn=!spinOn;
        if(!spinOn)stopSpin();
        updateSpinButton();
        modeInfoEl.textContent=spinOn?'Auto Rotate on: drag and release to spin':(MODE_INFO[currentMode]||currentMode);
    });
})();

// ===== Alt-drag: write the view transform into the coordinates =====
// A View-mode drag normally turns the view (rotQuat) or slides it (panX/panY)
// and leaves the file alone. Holding Alt during the press turns the very same
// drag into a coordinate edit: the structure keeps the orientation and position
// it was dragged to, and Save writes those numbers. The view state is what
// stays untouched, and that is exactly what keeps the picture identical to a
// plain drag — with p the on-screen position of an atom, R the view quaternion,
// C the pivot (the atoms' centroid) and pan the view offset:
//     plain drag:  p = pan' + R' * (a  - C)
//     Alt drag:    p = pan  + R  * (a' - C),   a' = C + (R^-1 R' R)(a - C)
// Rotation is taken about C, which a rotation about C leaves where it was. A
// translation freezes pan and moves the atoms by R^-1 * dpan; on release the
// pan the drag would have reached is committed too, and because the centre
// follows the atoms by exactly that translation the two cancel on screen — the
// picture stays where the drag put it while the file carries the movement.
// One undo snapshot per drag (dropped again for an Alt-click that moves
// nothing), and no bond is re-derived: a rigid transform changes no distance.
// Only the press decides: releasing Alt mid-drag keeps the drag an edit, and a
// drag that was not started with Alt never becomes one halfway through. Diff
// mode (two structures, two view states) and crystal structures (coordinates
// tied to the cell) keep the view-only behaviour.
var bakeActive=false;        // an Alt-drag that writes coordinates is running
var bakeKind='';             // 'rot' | 'pan'
var bakeButton=-1;           // the mouse button that started that drag
var bakeMoved=false;         // the drag actually changed a coordinate
var bakeOrigX=null,bakeOrigY=null,bakeOrigZ=null;  // coordinates at the press
var bakeCx=0,bakeCy=0,bakeCz=0;                    // pivot, frozen for the drag
var bakeQTotal=new THREE.Quaternion();             // rotation written so far
var bakeLocalDelta=new THREE.Vector3();            // translation written so far
var bakePan0={x:0,y:0};      // pan offset kept frozen for a baked translation
var bakePanTgt={x:0,y:0};    // the pan a plain drag would have reached
var bakeHintShown=false;     // the "not available" notice is shown once per mode
var _bakeV=new THREE.Vector3();

function bakeUnavailableReason(){
    if(diffMode)return 'Alt coordinate editing is not available in diff mode';
    if(CRY)return 'Alt coordinate editing is not available for crystal structures';
    return 'Alt coordinate editing needs a loaded structure';
}
// Everything that would make a baked edit wrong: a periodic cell, a second
// structure with its own view state, or nothing loaded at all.
function bakeAllowed(){
    return currentMode==='view'&&!diffMode&&!CRY&&MD.atoms.length>0;
}
function startBakeSession(kind,button){
    if(!bakeAllowed())return false;
    // Vibration owns the coordinates while it plays; a baked edit replaces them
    // wholesale, so the animation ends first (its own base positions are kept
    // for a rewind, which is what stopVibration does everywhere else).
    if(vibActive||vibPaused)stopVibration();
    // A frame auto-play would replace the structure under the snapshot a
    // moment later, so the drag takes the structure before it starts.
    stopFrameAutoPlay();
    var n=MD.atoms.length;
    bakeOrigX=new Float64Array(n);
    bakeOrigY=new Float64Array(n);
    bakeOrigZ=new Float64Array(n);
    for(var i=0;i<n;i++){
        var a=MD.atoms[i];
        if(!a)continue;
        bakeOrigX[i]=a.x;bakeOrigY[i]=a.y;bakeOrigZ[i]=a.z;
    }
    bakeActive=true;bakeKind=kind;bakeButton=button;bakeMoved=false;
    bakeQTotal.identity();bakeLocalDelta.set(0,0,0);
    bakePan0.x=panX;bakePan0.y=panY;
    bakePanTgt.x=panX;bakePanTgt.y=panY;
    pushUndo();
    // The snapshot may close a manual displacement, which re-centres the
    // structure, so the pivot is read after it rather than before.
    bakeCx=CX;bakeCy=CY;bakeCz=CZ;
    modeInfoEl.textContent=(kind==='rot'?'Alt: rotating ':'Alt: moving ')+n+
        ' atom(s) in the file - release to commit, Undo reverts';
    return true;
}
// Recompute every atom from the snapshot with the transform accumulated so far.
// Working from the snapshot (rather than from the previous frame) keeps the
// drag exact: one mousemove is one absolute result, nothing accumulates.
function bakeWriteAtoms(){
    // The snapshot describes one structure. If the array was replaced under the
    // drag (a frame or molecule switch that did not pass through the session
    // door), writing it would mix two structures — do nothing instead.
    if(!bakeOrigX||MD.atoms.length!==bakeOrigX.length)return;
    var n=MD.atoms.length;
    var q=bakeQTotal;
    var ox=bakeLocalDelta.x,oy=bakeLocalDelta.y,oz=bakeLocalDelta.z;
    var cx=bakeCx,cy=bakeCy,cz=bakeCz,v=_bakeV;
    for(var i=0;i<n;i++){
        var a=MD.atoms[i];
        if(!a){continue}
        v.set(bakeOrigX[i]-cx,bakeOrigY[i]-cy,bakeOrigZ[i]-cz).applyQuaternion(q);
        a.x=cx+v.x+ox;a.y=cy+v.y+oy;a.z=cz+v.z+oz;
    }
    bakeMoved=true;
    // Spheres and bond pieces only: the geometry itself never changes, and CX is
    // frozen for the drag, so both are a re-placement of existing objects.
    updateAtomMeshPositions();
    placeBondPieces();
    needsRender=true;
}
function bakeRotateBy(qx,qy){
    if(!bakeActive||bakeKind!=='rot')return;
    // The view premultiplies qx then qy, i.e. Q = qy*qx in the world frame.
    // Conjugate it into the molecule's own frame to write it into the
    // coordinates while the view quaternion stays exactly as it is.
    var qWorld=qy.clone().multiply(qx);
    var qLocal=rotQuat.clone().conjugate().multiply(qWorld).multiply(rotQuat);
    bakeQTotal.premultiply(qLocal);
    bakeQTotal.normalize();
    bakeWriteAtoms();
}
function bakeTranslateTo(dx,dy){
    if(!bakeActive||bakeKind!=='pan')return;
    // pan lives in the world frame; the atoms have to move along the same
    // world direction, which the view conjugation turns into molecule-local.
    bakeLocalDelta.set(dx,dy,0).applyQuaternion(rotQuat.clone().conjugate());
    bakeWriteAtoms();
}
function endBakeSession(){
    if(!bakeActive)return;
    bakeActive=false;
    if(!bakeMoved){
        // An Alt-click that changed nothing pushed a snapshot it does not need.
        if(undoStack.length>0){undoStack.pop();updateUndoBtn()}
    }else{
        // The centre follows the atoms, so a baked translation changes it by
        // exactly that translation. Committing the pan the drag would have
        // reached makes the two cancel on screen: the picture stays where the
        // drag put it while the coordinates carry the movement. Then the new
        // centroid becomes the rotation pivot again.
        provMark(bakeKind==='pan'?'translate':'rotate',1);
        if(bakeKind==='pan'){panX=bakePanTgt.x;panY=bakePanTgt.y}
        updateRotationCenter();
        if(bakeKind==='pan')updateTransform();
    }
    bakeOrigX=bakeOrigY=bakeOrigZ=null;
    bakeKind='';bakeButton=-1;bakeMoved=false;
    bakeQTotal.identity();bakeLocalDelta.set(0,0,0);
    needsRender=true;
    if(currentMode==='view')modeInfoEl.textContent=MODE_INFO.view;
}
// The panel states the feature where the mode it belongs to is switched on.
function updateBakeHint(){
    var hint=document.getElementById('vp-bake-hint');
    if(!hint)return;
    // CRY never changes within a session; only a comparison can come and go, and
    // every path that enters or leaves one calls this.
    if(diffMode||CRY){
        hint.textContent=bakeUnavailableReason()+'.';
    }else{
        hint.innerHTML='Hold <b>Alt</b> and drag: left rotates the structure, right/middle moves it. '+
            'The atom coordinates change, so Save writes the new geometry.';
    }
}
updateBakeHint();

function getClickedAtom(e){
    var rect=canvas.getBoundingClientRect();
    if(diffMode){
        var w=rect.width;
        var h=rect.height;
        var halfW=Math.floor(w/2);
        var localX=e.clientX-rect.left;
        if(localX<halfW){
            diffActiveSide='left';
            mouse.x=(localX/halfW)*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=halfW/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            var hitsL=raycaster.intersectObjects(atomMeshes);
            if(hitsL.length>0)return hitsL[0].object.userData.index;
        }else{
            diffActiveSide='right';
            mouse.x=((localX-halfW)/(w-halfW))*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=(w-halfW)/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            var hitsR=raycaster.intersectObjects(diffAtomMeshes);
            if(hitsR.length>0)return hitsR[0].object.userData.index;
        }
        return-1;
    }
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var hits=raycaster.intersectObjects(atomMeshes);
    if(hits.length>0)return hits[0].object.userData.index;
    return-1;
}

canvas.addEventListener('mousedown',function(e){
    // Any left press ends a free spin; a right/middle press pans and leaves it
    // running. This is the only mouse path that stops the spin.
    if(e.button===0)stopSpin();
    if(currentMode==='boxSelect'&&e.button===0){
        boxStart={x:e.clientX,y:e.clientY};
        isBoxSelecting=true;
        boxRectEl.style.display='block';
        boxRectEl.style.left=boxStart.x+'px';
        boxRectEl.style.top=boxStart.y+'px';
        boxRectEl.style.width='0px';boxRectEl.style.height='0px';
        e.preventDefault();return;
    }
    if(currentMode==='moveAtoms'&&!diffMode&&selectedAtoms.length>0&&e.ctrlKey&&(e.button===0||e.button===2)){
        // Ctrl+left-drag: rotate the selected fragment around its center.
        // Ctrl+right-drag: translate the fragment, cursor-locked onto the
        // grabbed atom (or the fragment-center plane when clicking space).
        // One undo snapshot per drag; a click without change drops it again.
        if(vibActive||vibPaused)stopVibration();
        moveDragOrig={};
        var mcx=0,mcy=0,mcz=0,mn=0;
        selectedAtoms.forEach(function(i){
            var a=MD.atoms[i];if(!a)return;
            moveDragOrig[i]={x:a.x,y:a.y,z:a.z};
            mcx+=a.x;mcy+=a.y;mcz+=a.z;mn++;
        });
        if(mn===0)return;
        mcx/=mn;mcy/=mn;mcz/=mn;
        pushUndo();
        if(e.button===0){
            // Fragment rotation: pivot is the fragment center, kept fixed in
            // absolute file coordinates for the whole drag.
            fragRotCx=mcx;fragRotCy=mcy;fragRotCz=mcz;
            fragRotActive=true;fragRotMoved=false;
            document.body.style.cursor='move';
            modeInfoEl.textContent='Rotating '+mn+' atom(s) - release Ctrl or mouse button to finish';
        }else{
            // Fragment translation. The movement plane goes through the world
            // position of the atom under the cursor (if it belongs to the
            // selection) so that atom stays exactly under the cursor; through
            // the fragment center otherwise. planeZ is derived from LOCAL
            // coordinates (atom - molecule center): using absolute file
            // coordinates can place the plane behind the camera, which mirrors
            // the ray intersection and inverts the drag direction.
            var hitIdx=getClickedAtom(e);
            var planeZ;
            if(hitIdx>=0&&selectedAtoms.indexOf(hitIdx)>=0&&MD.atoms[hitIdx]){
                planeZ=atomWorldZ(MD.atoms[hitIdx]);
            }else{
                planeZ=new THREE.Vector3(mcx-CX,mcy-CY,mcz-CZ).applyQuaternion(rotQuat).z;
            }
            moveDragPlaneZ=planeZ;
            moveDragStartWorld=screenToPlanePoint(e);
            if(!moveDragStartWorld)moveDragStartWorld=new THREE.Vector3(mcx-CX,mcy-CY,mcz-CZ).applyQuaternion(rotQuat);
            moveDragActive=true;moveDragMoved=false;
            document.body.style.cursor='move';
            modeInfoEl.textContent='Moving '+mn+' atom(s) - release Ctrl or mouse button to finish';
        }
        prevM={x:e.clientX,y:e.clientY};
        e.preventDefault();return;
    }
    // Alt in View mode: this press writes the drag into the coordinates instead
    // of only turning or sliding the view. Left rotates, middle/right translates
    // — the very same mapping the view uses. Decided here, once, so a drag can
    // never change its meaning halfway through; anywhere it is not available
    // (a crystal, a comparison) the drag stays view-only and says why once.
    if(e.altKey&&currentMode==='view'&&(e.button===0||e.button===1||e.button===2)){
        if(!startBakeSession(e.button===0?'rot':'pan',e.button)&&!bakeHintShown){
            bakeHintShown=true;
            modeInfoEl.textContent=bakeUnavailableReason();
        }
    }
    if(currentMode!=='view'&&e.button===0){
        var idx=getClickedAtom(e);
        if(idx>=0){selectAtom(idx);e.preventDefault();return}
    }
    if(e.button===0)isRot=true;
    else if(e.button===1||e.button===2){
        isPan=true;
        // Cursor-locked panning: the point under the cursor sticks to it.
        initPanGrab(e);
    }
    if(diffMode){
        var rect=canvas.getBoundingClientRect();
        var halfW=Math.floor(rect.width/2);
        diffTransformSide=(e.clientX-rect.left)<halfW?'left':'right';
    }
    prevM={x:e.clientX,y:e.clientY};
    e.preventDefault();
});

canvas.addEventListener('mousemove',function(e){
    if(isBoxSelecting){
        var bx=Math.min(boxStart.x,e.clientX),by=Math.min(boxStart.y,e.clientY);
        var bw=Math.abs(e.clientX-boxStart.x),bh=Math.abs(e.clientY-boxStart.y);
        boxRectEl.style.left=bx+'px';boxRectEl.style.top=by+'px';
        boxRectEl.style.width=bw+'px';boxRectEl.style.height=bh+'px';
        return;
    }
    if(moveDragActive){
        if(!e.ctrlKey){endMoveDrag();return}
        var cur=screenToPlanePoint(e);
        if(cur){
            var dw=cur.sub(moveDragStartWorld);
            if(dw.lengthSq()>1e-10)moveDragMoved=true;
            applyMoveWorldDelta(dw);
        }
        prevM={x:e.clientX,y:e.clientY};
        return;
    }
    if(fragRotActive){
        if(!e.ctrlKey){endFragRot();return}
        var rdm={x:e.clientX-prevM.x,y:e.clientY-prevM.y};
        if(rdm.x!==0||rdm.y!==0)fragRotMoved=true;
        applyFragRotDelta(rdm.y*0.008,rdm.x*0.008);
        prevM={x:e.clientX,y:e.clientY};
        return;
    }
    var dm={x:e.clientX-prevM.x,y:e.clientY-prevM.y};
    if(diffMode){
        var side=diffTransformSide;
        var rq=side==='right'?diffRotQuat:rotQuat;
        if(isRot){
            spinRecordPoint(e.clientX,e.clientY);
            var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
            var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
            rq.premultiply(qx);rq.premultiply(qy);rq.normalize();
            updateTransform()
        }
        if(isPan){
            // Pan the viewport that was actually grabbed. In diff mode the two
            // sides keep separate pan offsets, so a right-drag on the right
            // half must move the comparison structure, not the original.
            if(panGrabActive){
                // Cursor-locked: recompute this side's pan so the grabbed point
                // follows the cursor 1:1 (same world-z reference plane as at grab).
                if(panGrabSide==='right'){
                    var pcr=sideScreenToPlanePoint(e,'right',diffPanRefZ);
                    if(pcr){diffPanX=pcr.x-diffPanGrabOffset.x;diffPanY=pcr.y-diffPanGrabOffset.y}
                }else{
                    var pcl=sideScreenToPlanePoint(e,'left',panRefZ);
                    if(pcl){panX=pcl.x-panGrabOffset.x;panY=pcl.y-panGrabOffset.y}
                }
            }else{
                var pk=0.01*((panGrabSide==='right'?diffCamDist:camDist)/20);
                if(panGrabSide==='right'){diffPanX+=dm.x*pk;diffPanY-=dm.y*pk}
                else{panX+=dm.x*pk;panY-=dm.y*pk}
            }
            updateTransform()
        }
    }else{
        if(isRot){
            spinRecordPoint(e.clientX,e.clientY);
            var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
            var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
            if(bakeActive&&bakeKind==='rot'){
                // An Alt-drag writes the rotation into the coordinates and leaves
                // the view quaternion alone, so the two can never turn it twice.
                if(dm.x||dm.y)bakeRotateBy(qx,qy);
            }else{
                rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();
                updateTransform()
            }
        }
        if(isPan){
            var bakePan=bakeActive&&bakeKind==='pan';
            if(panGrabActive){
                // Cursor-locked: recompute pan so the grabbed point follows
                // the cursor 1:1 (same world-z reference plane as at grab).
                var pc2=screenPointOnPlane(e,panRefZ);
                if(pc2){
                    if(bakePan){bakePanTgt.x=pc2.x-panGrabOffset.x;bakePanTgt.y=pc2.y-panGrabOffset.y}
                    else{panX=pc2.x-panGrabOffset.x;panY=pc2.y-panGrabOffset.y}
                }
            }else if(bakePan){
                // No grab reference: accumulate the plain-drag pan a would-be
                // view pan would have reached, keeping panX/panY frozen instead.
                bakePanTgt.x+=dm.x*0.01*(camDist/20);
                bakePanTgt.y-=dm.y*0.01*(camDist/20);
            }else{
                panX+=dm.x*0.01*(camDist/20);panY-=dm.y*0.01*(camDist/20)
            }
            if(bakePan){
                // The whole translation since the press, applied in one go: the
                // pan offset stays where it was, the atoms take the movement.
                var bdx=bakePanTgt.x-bakePan0.x,bdy=bakePanTgt.y-bakePan0.y;
                if(bdx||bdy)bakeTranslateTo(bdx,bdy);
            }else updateTransform()
        }
    }
    prevM={x:e.clientX,y:e.clientY};
    var rect=canvas.getBoundingClientRect();
    if(diffMode){
        var w=rect.width;
        var h=rect.height;
        var halfW=Math.floor(w/2);
        var localX=e.clientX-rect.left;
        var hitsD;
        if(localX<halfW){
            mouse.x=(localX/halfW)*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=halfW/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            hitsD=raycaster.intersectObjects(atomMeshes);
        }else{
            mouse.x=((localX-halfW)/(w-halfW))*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=(w-halfW)/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            hitsD=raycaster.intersectObjects(diffAtomMeshes);
        }
        if(hitsD&&hitsD.length>0){
            var o=hitsD[0].object,i=o.userData.index;
            var a=diffMode&&o.userData.diffHi!==undefined&&localX>=halfW?diffData.atoms[i]:MD.atoms[i];
            tooltipEl.textContent=a.element+(i+1)+' ('+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4)+')';
            tooltipEl.style.display='block';
            tooltipEl.style.left=(e.clientX-container.getBoundingClientRect().left+15)+'px';
            tooltipEl.style.top=(e.clientY-container.getBoundingClientRect().top-10)+'px';
        }else{tooltipEl.style.display='none'}
        return;
    }
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var hits=raycaster.intersectObjects(atomMeshes);
    if(hits.length>0){var o=hits[0].object,i=o.userData.index,a=MD.atoms[i];
        tooltipEl.textContent=a.element+(i+1)+' ('+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4)+')';
        tooltipEl.style.display='block';
        tooltipEl.style.left=(e.clientX-container.getBoundingClientRect().left+15)+'px';
        tooltipEl.style.top=(e.clientY-container.getBoundingClientRect().top-10)+'px';
    }else{tooltipEl.style.display='none'}
});

canvas.addEventListener('mouseup',function(e){
    if(moveDragActive){endMoveDrag()}
    if(fragRotActive){endFragRot()}
    // An Alt-drag commits on the release of the button that started it. A bake
    // drag does not throw the view afterwards: the drag was an edit, and a spin
    // launched from it would keep turning what was just committed.
    var bakedThisRelease=bakeActive&&e.button===bakeButton;
    if(bakedThisRelease)endBakeSession();
    panGrabActive=false;
    if(isBoxSelecting){
        isBoxSelecting=false;boxRectEl.style.display='none';
        var rect=canvas.getBoundingClientRect();
        var bx=Math.min(boxStart.x,e.clientX),by=Math.min(boxStart.y,e.clientY);
        var bw=Math.abs(e.clientX-boxStart.x),bh=Math.abs(e.clientY-boxStart.y);
        if(bw<3&&bh<3){boxStart=null;layoutPanels();return;}
        var v=new THREE.Vector3();
        if(diffMode){
            var dw=rect.width,dh=rect.height,dhalfW=Math.floor(dw/2);
            diffActiveSide='left';
            camera.aspect=dhalfW/dh;camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);camera.updateMatrixWorld(true);
            atomMeshes.forEach(function(m){
                m.getWorldPosition(v);v.project(camera);
                var sx=(v.x*0.5+0.5)*dhalfW+rect.left;
                var sy=(-v.y*0.5+0.5)*dh+rect.top;
                if(sx>=bx&&sx<=bx+bw&&sy>=by&&sy<=by+bh){selectAtom(m.userData.index);}
            });
            if(diffAtomMeshes&&diffAtomMeshes.length>0){
                diffActiveSide='right';
                camera.aspect=(dw-dhalfW)/dh;camera.updateProjectionMatrix();
                camera.position.set(0,0,diffCamDist);camera.updateMatrixWorld(true);
                diffAtomMeshes.forEach(function(m){
                    m.getWorldPosition(v);v.project(camera);
                    var sx=(v.x*0.5+0.5)*(dw-dhalfW)+rect.left+dhalfW;
                    var sy=(-v.y*0.5+0.5)*dh+rect.top;
                    if(sx>=bx&&sx<=bx+bw&&sy>=by&&sy<=by+bh){selectAtom(m.userData.index);}
                });
            }
        }else{
            atomMeshes.forEach(function(m){
                m.getWorldPosition(v);v.project(camera);
                var sx=(v.x*0.5+0.5)*rect.width+rect.left;
                var sy=(-v.y*0.5+0.5)*rect.height+rect.top;
                if(sx>=bx&&sx<=bx+bw&&sy>=by&&sy<=by+bh){selectAtom(m.userData.index);}
            });
        }
        boxStart=null;needsRender=true;layoutPanels();return;
    }
    if(e.button===0&&isRot&&!bakedThisRelease)spinLaunch();
    isRot=false;isPan=false;layoutPanels();
});
canvas.addEventListener('mouseleave',function(){if(bakeActive)endBakeSession();isRot=false;isPan=false;panGrabActive=false;spinSamples.length=0;tooltipEl.style.display='none';if(moveDragActive)endMoveDrag();if(fragRotActive)endFragRot();if(isBoxSelecting){isBoxSelecting=false;boxRectEl.style.display='none';boxStart=null}});
canvas.addEventListener('wheel',function(e){e.preventDefault();
    if(diffMode){
        var rect=canvas.getBoundingClientRect();
        var halfW=Math.floor(rect.width/2);
        var side=(e.clientX-rect.left)<halfW?'left':'right';
        if(side==='right'){
            diffCamDist*=e.deltaY>0?1.1:0.9;
            diffCamDist=Math.max(1,Math.min(500,diffCamDist));
        }else{
            camDist*=e.deltaY>0?1.1:0.9;
            camDist=Math.max(1,Math.min(500,camDist));
        }
    }else{
        camDist*=e.deltaY>0?1.1:0.9;
        camDist=Math.max(1,Math.min(500,camDist));
        camera.position.z=camDist;
    }
    needsRender=true;layoutPanels();
},{passive:false});
canvas.addEventListener('contextmenu',function(e){e.preventDefault()});

var touchSD=0;
canvas.addEventListener('touchstart',function(e){e.preventDefault();
    if(e.touches.length===1){
        stopSpin();
        isRot=true;prevM={x:e.touches[0].clientX,y:e.touches[0].clientY};
        if(diffMode){var rect=canvas.getBoundingClientRect();var halfW=Math.floor(rect.width/2);diffTransformSide=(e.touches[0].clientX-rect.left)<halfW?'left':'right'}
    }
    else if(e.touches.length===2){isRot=false;var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;touchSD=Math.sqrt(dx*dx+dy*dy)}
},{passive:false});
canvas.addEventListener('touchmove',function(e){e.preventDefault();
    if(e.touches.length===1&&isRot){
        spinRecordPoint(e.touches[0].clientX,e.touches[0].clientY);
        var dm={x:e.touches[0].clientX-prevM.x,y:e.touches[0].clientY-prevM.y};
        var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
        var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
        if(diffMode){
            var rq=diffTransformSide==='right'?diffRotQuat:rotQuat;
            rq.premultiply(qx);rq.premultiply(qy);rq.normalize();
        }else{
            rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();
        }
        updateTransform();
        prevM={x:e.touches[0].clientX,y:e.touches[0].clientY}
    }
    else if(e.touches.length===2){var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY,d=Math.sqrt(dx*dx+dy*dy);
        if(touchSD>0){
            if(diffMode){
                var sc=touchSD/d;
                if(diffTransformSide==='right'){diffCamDist*=sc;diffCamDist=Math.max(1,Math.min(500,diffCamDist))}
                else{camDist*=sc;camDist=Math.max(1,Math.min(500,camDist))}
            }else{
                camDist*=touchSD/d;camDist=Math.max(1,Math.min(500,camDist));camera.position.z=camDist;
            }
        }
        touchSD=d}
    needsRender=true;
},{passive:false});
canvas.addEventListener('touchend',function(e){if(e.touches.length===0&&isRot)spinLaunch();isRot=false;if(e.touches.length<2)touchSD=0});

// ===== Keyboard shortcuts =====
// The page keeps one keyboard door. Esc, Ctrl/Cmd+Z, Ctrl/Cmd+S and the
// selection/clipboard keys Ctrl/Cmd+A and Ctrl/Cmd+X/C/V are the keys a hand
// reaches for before the mouse, so they are settled here, in front of the
// viewer's own key handling (arrows, Delete). A keystroke that belongs to a text
// field stays the field's own — Ctrl+Z there is the browser's text undo, Ctrl+A
// its select-all, Ctrl+V its paste — with Esc as the one exception, because Esc
// is how a user gets out of a dialog that a field sits in.
var IS_MAC=/^(Mac|iPhone|iPad|iPod)/.test(navigator.platform||'');
var KEY_MOD=IS_MAC?'Cmd':'Ctrl';
function keyTargetEditable(t){if(!t)return false;var n=t.tagName;return n==='INPUT'||n==='TEXTAREA'||n==='SELECT'||t.isContentEditable===true}
function modalOpen(){return !!(modalOverlay&&modalOverlay.classList.contains('show'))}
/**
 * Esc: back out of exactly one state, in the order the user built them up.
 * Every step takes back something the user started and nothing else — a
 * dismissed dialog, an abandoned rectangle, a tool left for View, a selection
 * let go — so a second press can only undo a smaller thing, never the structure.
 */
function escapeOneLevel(){
    if(modalOpen()){
        // The dialog's own Cancel button is the path that knows what
        // cancelling means (some of them put back the geometry they were
        // previewing); only a dialog without one falls back to plain hiding.
        var cancel=modalEl.querySelector('.mbtn-cancel');
        if(cancel){cancel.click();return true}
        hideModal();return true;
    }
    if(isBoxSelecting){
        // A rectangle still being dragged is thrown away: it has selected
        // nothing yet, so there is nothing to take back.
        isBoxSelecting=false;boxRectEl.style.display='none';boxStart=null;
        needsRender=true;layoutPanels();return true;
    }
    if(currentMode!=='view'){setMode('view');return true}
    if(selectedAtoms.length>0||diffSelectedAtoms.length>0){resetSelection();return true}
    return false;
}
(function(){
    // The two mouse-free actions say their own keys, on the buttons that do them.
    var ub=document.getElementById('undo-btn');if(ub)ub.title='Undo ('+KEY_MOD+'+Z)';
    var sb=document.getElementById('save-btn');if(sb)sb.title='Save As ('+KEY_MOD+'+S)';
})();
document.addEventListener('keydown',function(e){
    if(e.key==='Escape'&&!e.ctrlKey&&!e.metaKey&&!e.altKey){
        if(escapeOneLevel()){e.preventDefault();return}
    }
    var editable=keyTargetEditable(e.target);
    var mod=e.ctrlKey||e.metaKey;
    if(mod&&!e.altKey&&!e.shiftKey){
        // e.key follows the layout (on AZERTY Ctrl+Z arrives as 'z'), e.code
        // does not (on a Cyrillic layout it arrives as 'я'), so a shortcut
        // answers to either: one physical press, one action, any keyboard.
        var lk=(e.key||'').toLowerCase();
        if(lk==='z'||e.code==='KeyZ'){
            // Not while a dialog is open: it previews geometry on top of an
            // undo step of its own, and an undo behind its back would split
            // the two apart — its Cancel would then pop someone else's step.
            if(!editable&&!modalOpen()){e.preventDefault();doUndo();return}
        }else if(lk==='s'||e.code==='KeyS'){
            if(modalOpen()){
                // The Save dialog is the one dialog that can finish its own
                // job from the keyboard: a second Ctrl+S presses its Save
                // button, keeping the format already chosen.
                var fmtEl=document.getElementById('m-fmt'),okEl=document.getElementById('m-ok');
                if(fmtEl&&okEl){e.preventDefault();okEl.click();return}
            }else if(!editable){
                e.preventDefault();doSave();return;
            }
        }else if(lk==='a'||e.code==='KeyA'){
            // Select every atom. Outside a text field the page has no use for
            // select-all, so the key is consumed either way: a refusal (a dialog
            // is open) must not quietly select the page's text instead. Auto-
            // repeat is ignored — each repeat would rebuild the whole selection
            // and re-highlight every atom.
            if(!editable&&!e.repeat){e.preventDefault();selectAllAtoms();return}
        }else if(lk==='c'||e.code==='KeyC'){
            // Copy is read-only, so it is allowed wherever there is a selection
            // — a dialog or a comparison on screen does not make it unsafe.
            // Auto-repeat is ignored: a held key would re-copy the same atoms
            // dozens of times. The page has no selectable text (see the
            // user-select rules in the stylesheet), so there is nothing for the
            // browser's own copy to take here; a field's text is excluded above.
            if(!editable&&!e.repeat&&copySelection()){e.preventDefault();return}
        }else if(lk==='x'||e.code==='KeyX'){
            // Cut deletes, so it waits until no dialog is previewing the
            // structure it is about to change.
            if(!editable&&!e.repeat&&!modalOpen()&&cutSelection()){e.preventDefault();return}
        }else if(lk==='v'||e.code==='KeyV'){
            // Paste grows the structure: never behind a dialog, and never on a
            // key repeat, or holding the key would paste again and again.
            if(!editable&&!e.repeat&&!modalOpen()&&pasteClipboard()){e.preventDefault();return}
        }
    }
    if(editable)return;
    ctrlDown=e.ctrlKey;
    // The first arrow key of a hold-session stops whatever spin is running, so
    // the key takes over; keys pressed while another is already held (and the
    // auto-repeat of a held key) must not, or every repeat would restart the
    // spin-up ramp from its base speed.
    if(e.key in keyState){
        if(!keyState.ArrowLeft&&!keyState.ArrowRight&&!keyState.ArrowUp&&!keyState.ArrowDown)stopSpin();
        keyState[e.key]=true;e.preventDefault();
    }
    if(e.key==='Delete'){
        e.preventDefault();
        if(selectedAtoms.length>0){showBatchDeleteModal()}
        else{setMode('deleteAtom')}
    }
});
document.addEventListener('keyup',function(e){
    ctrlDown=e.ctrlKey;
    if(e.key in keyState){
        keyState[e.key]=false;e.preventDefault();
        // Clear the move-session latches once every arrow key is released, so
        // fresh sessions may start on the next key press. An arrow-key spin
        // ends its ramp here too: the speed it reached stays and becomes the
        // free rotation, while the next hold starts its ramp over again.
        if(!keyState.ArrowLeft&&!keyState.ArrowRight&&!keyState.ArrowUp&&!keyState.ArrowDown){moveKeyLatch=false;keyRotLatch=false;keySpinActive=false;keySpinHold=0}
    }
});
window.addEventListener('blur',function(){
    ctrlDown=false;
    keyState.ArrowLeft=keyState.ArrowRight=keyState.ArrowUp=keyState.ArrowDown=false;
    keySpinActive=false;keySpinHold=0;
    // The mouse-up may never arrive after focus is lost: commit what the drag
    // already wrote rather than leaving the session (and its snapshot) open.
    if(bakeActive)endBakeSession();
});

window.addEventListener('resize',function(){var rw=container.clientWidth||window.innerWidth;var rh=container.clientHeight||(window.innerHeight-60);if(rw<1)rw=window.innerWidth;if(rh<1)rh=window.innerHeight-60;camera.aspect=rw/rh;camera.updateProjectionMatrix();renderer.setSize(rw,rh);needsRender=true;layoutPanels()});

// Hand the context back when the page is really going away. Editors keep their
// webviews (and their contexts) alive while hidden, and the window has a hard
// cap on how many WebGL contexts may exist at once — a viewer that is finished
// with its canvas should not be one of them.
window.addEventListener('pagehide',function(){
    try{
        if(renderer){
            renderer.dispose();
            if(typeof renderer.forceContextLoss==='function')renderer.forceContextLoss();
        }
    }catch(e){}
},false);

function animate(){
    requestAnimationFrame(animate);
    var spinNow=performance.now();
    var spinDt=spinPrevT?(spinNow-spinPrevT)/1000:0;
    spinPrevT=spinNow;
    // One clamp for every rotation source: a tab that sat in the background for
    // a while must not hand the molecule a several-second jump on its first frame.
    if(spinDt>0.1)spinDt=0.1;
    var keyActive=false;
    // A running Alt-drag owns the orientation: an arrow key must not turn the
    // view under it (that would apply the drag's rotation twice) nor start a
    // spin. The keys are simply ignored until the drag commits.
    if(!bakeActive&&(keyState.ArrowLeft||keyState.ArrowRight||keyState.ArrowUp||keyState.ArrowDown)){
        keyActive=true;
        if(moveKeyActive&&!ctrlDown){
            // Ctrl released mid-hold: end the fragment rotation session and
            // latch view rotation until all arrow keys are released.
            moveKeyActive=false;
            keyRotLatch=true;
            refreshMovedBondsAndMeshes();
            if(currentMode==='moveAtoms')modeInfoEl.textContent=MODE_INFO.moveAtoms;
        }
        if(currentMode==='moveAtoms'&&selectedAtoms.length>0&&!diffMode&&!moveKeyLatch&&ctrlDown){
            // Ctrl+arrow keys rotate the selected fragment instead of the
            // view. One undo snapshot per key-hold session; releasing all
            // arrow keys (or Ctrl) ends the session.
            var fAx=0,fAy=0;
            if(keyState.ArrowRight)fAy+=keyRotSpeed;
            if(keyState.ArrowLeft)fAy-=keyRotSpeed;
            if(keyState.ArrowDown)fAx+=keyRotSpeed;
            if(keyState.ArrowUp)fAx-=keyRotSpeed;
            if(fAx!==0||fAy!==0){
                if(!moveKeyActive){
                    moveKeyActive=true;
                    if(vibActive||vibPaused)stopVibration();
                    pushUndo();
                    // Rotation pivot = fragment center, fixed for the session.
                    var fcx=0,fcy=0,fcz=0,fn=0;
                    selectedAtoms.forEach(function(i){var a=MD.atoms[i];if(!a)return;fcx+=a.x;fcy+=a.y;fcz+=a.z;fn++});
                    if(fn>0){fcx/=fn;fcy/=fn;fcz/=fn}
                    fragRotCx=fcx;fragRotCy=fcy;fragRotCz=fcz;
                    modeInfoEl.textContent='Rotating '+fn+' atom(s) - release Ctrl or arrow keys to finish';
                }
                applyFragRotDelta(fAx,fAy);
            }else if(moveKeyActive){
                moveKeyActive=false;
                refreshMovedBondsAndMeshes();
                if(currentMode==='moveAtoms')modeInfoEl.textContent=MODE_INFO.moveAtoms;
            }
        }else if(!diffMode&&!keyRotLatch){
            if(spinOn){
                // Auto Rotate: a held arrow key spins the view up rather than
                // turning it at a fixed rate. Nothing is applied here — the one
                // spin site below turns whatever the drive just set.
                keySpinDrive(spinDt);
            }else{
                if(keyState.ArrowLeft){var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),-keyRotSpeed);rotQuat.premultiply(qx)}
                if(keyState.ArrowRight){var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),keyRotSpeed);rotQuat.premultiply(qx)}
                if(keyState.ArrowUp){var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),-keyRotSpeed);rotQuat.premultiply(qy)}
                if(keyState.ArrowDown){var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),keyRotSpeed);rotQuat.premultiply(qy)}
                rotQuat.normalize();
                updateTransform();
            }
        }else if(diffMode){
            if(spinOn){
                keySpinDrive(spinDt);
            }else{
                var rq=(diffTransformSide==='right')?diffRotQuat:rotQuat;
                if(keyState.ArrowLeft){var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),-keyRotSpeed);rq.premultiply(qx)}
                if(keyState.ArrowRight){var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),keyRotSpeed);rq.premultiply(qx)}
                if(keyState.ArrowUp){var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),-keyRotSpeed);rq.premultiply(qy)}
                if(keyState.ArrowDown){var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),keyRotSpeed);rq.premultiply(qy)}
                rq.normalize();
                updateTransform();
            }
        }
    }
    if(spinQuatRef&&spinSpeed>0){
        // The one place a spin is applied, whatever set it — a released drag
        // (spinLaunch) or a held arrow key (keySpinDrive). Free rotation, with
        // no damping, so the view keeps turning until something stops it; the
        // frame step was clamped above.
        if(spinDt>0){
            var spinQ=new THREE.Quaternion().setFromAxisAngle(spinAxis,spinSpeed*spinDt);
            spinQuatRef.premultiply(spinQ);spinQuatRef.normalize();
            updateTransform();
        }
    }
    if(!needsRender&&!keyActive)return;
    needsRender=false;
    updateAxesIndicator();
    var w=container.clientWidth||window.innerWidth;
    var h=container.clientHeight||(window.innerHeight-60);
    if(w<1)w=window.innerWidth;
    if(h<1)h=window.innerHeight-60;
    if(diffMode){
        var halfW=Math.floor(w/2);
        renderer.autoClear=false;
        renderer.setScissorTest(true);
        renderer.setViewport(0,0,halfW,h);
        renderer.setScissor(0,0,halfW,h);
        camera.aspect=halfW/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,camDist);
        renderer.clear();
        if(diffPivot)diffPivot.visible=false;
        pivotGroup.visible=true;
        renderer.render(scene,camera);
        renderer.setViewport(halfW,0,w-halfW,h);
        renderer.setScissor(halfW,0,w-halfW,h);
        camera.aspect=(w-halfW)/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,diffCamDist);
        renderer.clear();
        if(diffPivot)diffPivot.visible=true;
        pivotGroup.visible=false;
        renderer.render(scene,camera);
        renderer.setScissorTest(false);
        renderer.autoClear=true;
        pivotGroup.visible=true;
        if(diffPivot)diffPivot.visible=true;
    }else{
        camera.aspect=w/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,camDist);
        renderer.render(scene,camera);
    }
}
animate();
}catch(e){var el=document.getElementById('error-msg');var ll=document.getElementById('loading');if(ll)ll.style.display='none';if(el){el.style.display='block';el.textContent='Error: '+e.message}}
})();
</script>
</body>
</html>`;
    }
}

class MolecularDocument implements vscode.CustomDocument {
    constructor(
        public readonly uri: vscode.Uri,
        public readonly data: MolecularData,
        public readonly frames: (LogFrame | OrcaFrame)[] = [],
        public readonly optSteps?: OptStep[],
        public readonly normalModes?: NormalMode[],
        public readonly logSource?: string,
        public readonly routes?: RouteSection[],
        public readonly vmdPayload?: VmdPayload,
        public readonly autoBeta: boolean = false,
        public readonly optThresholds?: OptThresholds
    ) {}

    dispose(): void {}
}

function getNonce(): string {
    let text = '';
    const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    for (let i = 0; i < 32; i++) {
        text += possible.charAt(Math.floor(Math.random() * possible.length));
    }
    return text;
}

// ---------------------------------------------------------------------------
// Import helpers
// ---------------------------------------------------------------------------

/** 一次导入的份数上限, 与新增原子总数上限同时生效。 */
const IMPORT_COPY_LIMIT = 50;
const IMPORT_ATOM_LIMIT = 20000;

/** 对话框里的来源选项: 用自己的字段而不是标签文本做判定。 */
interface ImportSourceItem extends vscode.QuickPickItem {
    source: 'file' | 'smiles';
}

/** 用户选定的来源与份数。 */
interface ImportChoice {
    source: 'file' | 'smiles';
    copies: number;
}

/** 解析对话框输入框里的份数; 非法返回 null (空、非数字、0、负数、超上限)。 */
function parseCopyCount(text: string): number | null {
    const trimmed = text.trim();
    if (!/^[0-9]+$/.test(trimmed)) {
        return null;
    }
    const count = parseInt(trimmed, 10);
    if (count < 1 || count > IMPORT_COPY_LIMIT) {
        return null;
    }
    return count;
}

/**
 * 弹出导入对话框: **输入框填份数, 列表里点来源** —— 先定份数, 再定结构。
 *
 * 输入框是这个对话框本来就有的那个 (快速选择总会在列表上方渲染一个输入框);
 * 两个来源选项都带 alwaysShow, 因此往里打数字不会把选项过滤掉, 它就从"没有用
 * 的筛选框"变成了"份数输入框"。
 *
 * 份数不合法时带着原输入重新弹出, 原因写在标题里 —— 不弹模态提示, 免得打断,
 * 用户也可以直接把数字改掉。返回 undefined 表示取消。
 */
async function askImportChoice(): Promise<ImportChoice | undefined> {
    let preset = '1';
    let complaint = '';
    for (;;) {
        const items: ImportSourceItem[] = [
            { label: 'From File', description: 'pick a structure file (.xyz, .gjf, .mol2, .log, ...)', source: 'file', alwaysShow: true },
            { label: 'From SMILES', description: 'type a SMILES string; coordinates are generated locally', source: 'smiles', alwaysShow: true }
        ];
        const quickPick = vscode.window.createQuickPick<ImportSourceItem>();
        quickPick.title = complaint === '' ? 'Import structure' : 'Import structure — ' + complaint;
        quickPick.placeholder = 'Copies to import (1-' + IMPORT_COPY_LIMIT + '), then pick a source';
        quickPick.value = preset;
        quickPick.ignoreFocusOut = true;
        quickPick.matchOnDescription = false;
        quickPick.items = items;
        quickPick.activeItems = [items[0]];
        const picked = await new Promise<{ item: ImportSourceItem; text: string } | undefined>((resolve) => {
            let accepted = false;
            quickPick.onDidAccept(() => {
                accepted = true;
                const chosen = quickPick.selectedItems.length > 0
                    ? quickPick.selectedItems[0]
                    : (quickPick.activeItems.length > 0 ? quickPick.activeItems[0] : items[0]);
                resolve({ item: chosen, text: quickPick.value });
                quickPick.hide();
            });
            quickPick.onDidHide(() => {
                quickPick.dispose();
                if (!accepted) {
                    resolve(undefined);          // Esc / 被其它界面挤掉 = 取消
                }
            });
            quickPick.show();
        });
        if (picked === undefined) {
            return undefined;
        }
        const copies = parseCopyCount(picked.text);
        if (copies !== null) {
            return { source: picked.item.source, copies: copies };
        }
        preset = picked.text;                    // 保留原输入, 直接改数字即可
        complaint = 'copies must be a whole number from 1 to ' + IMPORT_COPY_LIMIT;
    }
}

/** 相邻两份之间的间隔参数: 与客户端摆放使用的缓冲/间隙保持一致。 */
const IMPORT_PLACEMENT_BUFFER = 2.5;
const IMPORT_PLACEMENT_GAP = 2.0;

/**
 * 把 N 份拷贝拼成一个整体。份与份之间用**同一个**间距 (
 * 两份各自的包围球半径 + 缓冲 + 间隙), 以共同质心为中心均匀排开 —— 于是每两份
 * 之间的最小距离都与"单独导入两次"所保证的一样, 不会重叠。
 *
 * 单份时原样返回 (不复制数组), 保证"只导入一份"的行为与以前逐位相同。
 */
function buildImportBlock(atoms: ImportAtomPayload[], bonds: ImportBondPayload[],
                          copies: number): { atoms: ImportAtomPayload[]; bonds: ImportBondPayload[] } {
    if (copies <= 1 || atoms.length === 0) {
        return { atoms: atoms, bonds: bonds };
    }
    // 单份的包围球 (与客户端摆放判定用的是同一个量)
    let cx = 0.0;
    let cy = 0.0;
    let cz = 0.0;
    for (let i = 0; i < atoms.length; i++) {
        cx += atoms[i].x;
        cy += atoms[i].y;
        cz += atoms[i].z;
    }
    cx /= atoms.length;
    cy /= atoms.length;
    cz /= atoms.length;
    let radius = 0.0;
    for (let i = 0; i < atoms.length; i++) {
        const dx = atoms[i].x - cx;
        const dy = atoms[i].y - cy;
        const dz = atoms[i].z - cz;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > radius) {
            radius = d;
        }
    }
    const stride = 2.0 * (radius + IMPORT_PLACEMENT_BUFFER) + IMPORT_PLACEMENT_GAP;
    const outAtoms: ImportAtomPayload[] = [];
    const outBonds: ImportBondPayload[] = [];
    for (let c = 0; c < copies; c++) {
        const shift = (c - (copies - 1) / 2.0) * stride;
        const base = c * atoms.length;
        for (let i = 0; i < atoms.length; i++) {
            const a = atoms[i];
            outAtoms.push({ element: a.element, x: a.x + shift, y: a.y, z: a.z, color: a.color });
        }
        for (let k = 0; k < bonds.length; k++) {
            const b = bonds[k];
            outBonds.push({ atom1: base + b.atom1, atom2: base + b.atom2, order: b.order });
        }
    }
    return { atoms: outAtoms, bonds: outBonds };
}

/**
 * 按份数拼好后**一次性**发出 (只发一条导入消息)。
 *
 * 位置计算因此只发生一次: 客户端对"整个块"做一次包围球避让摆放 (第 1 份的位置
 * 由这次计算确定), 份与份之间则是上面那个等间距, 客户端还会一次成键、一次场景
 * 重建, 撤销也只要一步 —— 逐份重算位置、逐份重建场景的开销全部省掉。
 *
 * 份数已在对话框里校验; 这里再按"新增原子总数"上限截断 (此时才知道结构有多大)
 * 并告知用户, 免得一次导入把视图撑到无法交互。
 */
function postImportCopies(webview: vscode.Webview, fileName: string,
                          atoms: ImportAtomPayload[], bonds: ImportBondPayload[],
                          hasExplicitBonds: boolean, copies: number): void {
    const affordable = Math.max(1, Math.floor(IMPORT_ATOM_LIMIT / Math.max(1, atoms.length)));
    let count = copies;
    if (count > affordable) {
        count = affordable;
        vscode.window.showWarningMessage('Import: ' + copies + ' copies of a ' + atoms.length +
            '-atom structure would add ' + (copies * atoms.length) + ' atoms, so this import is limited to ' +
            count + ' copies.');
    }
    const block = buildImportBlock(atoms, bonds, count);
    webview.postMessage({
        command: 'importResult',
        cancelled: false,
        fileName: count > 1 ? fileName + ' x ' + count : fileName,
        atoms: block.atoms,
        bonds: block.bonds,
        hasExplicitBonds: hasExplicitBonds
    });
}
