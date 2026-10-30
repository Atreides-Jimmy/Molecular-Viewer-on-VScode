import * as vscode from 'vscode';
import * as path from 'path';
import { parseFile, parseLogFile, parseXyzTrajectory, LogFrame, OrcaFrame, RouteSection, parseTcl, resolveVmdScene, resolveSelection, vdwRadius, hasScalarChannel, parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg } from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { BOND_PERCEPTION_SOURCE } from '../tools/bondPerception';
import { BVS_SOURCE } from '../tools/bvs';
import { addHydrogens } from '../tools/addHydrogens';
import { beginOptimization, describeForceField, MmHistoryPoint, MmOptimizeResult, MmSession } from '../tools/mmOptimizer';
import { smilesTo3dSafe, validateSmiles } from '../tools/smilesTo3d';
import { MolecularData, Atom, AsymAtom, AtomGroup, OptStep, NormalMode, OptThresholds, VmdScene, VmdMolecule } from '../types';

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
    /** A ghost/dummy atom of the loaded structure (see `Atom.ghost`). */
    ghost?: boolean;
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
        if (a.ghost === true) item.ghost = true;
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

/**
 * The payload shape of one base-cell atom. The page keeps its own copy of the
 * fields a tooltip and a save need, so the mapping lives here once and is used
 * by both the file's expansion and a space-group alternative's.
 */
function crystalBasePayload(a: Atom) {
    return {
        element: a.element, x: a.x, y: a.y, z: a.z,
        index: a.index, occupancy: a.occupancy, baseIdx: a.baseIdx,
        siteLabel: a.siteLabel || undefined,
        asymIdx: a.asymIdx, symOp: a.symOp,
        constraint: a.constraint, disorderGroup: a.disorderGroup,
        uAniso: a.uAniso
    };
}

/** The payload shape of one asymmetric-unit site (label, frac, ops). */
function crystalAsymPayload(a: AsymAtom) {
    return {
        element: a.element, label: a.label, x: a.x, y: a.y, z: a.z,
        frac: a.frac, index: a.index, occupancy: a.occupancy,
        symmetryOps: a.symmetryOps, disorderGroup: a.disorderGroup,
        uAniso: a.uAniso
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
                                // The powder XRD panel exports a two-column
                                // 2theta/intensity text file.
                                'XY Data': ['xy', 'dat'],
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
                        } | null = null;

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
                                    }))
                                };
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

                        const impAtoms = impData.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: ATOM_COLORS[a.element] || '#FF1493', ghost: a.ghost === true ? true : undefined }));
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
                ghost: a.ghost === true ? true : undefined,
                baseIdx: a.baseIdx,
                cellI: a.cellI, cellJ: a.cellJ, cellK: a.cellK,
                uAniso: a.uAniso
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
            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: atomColors[a.element] || '#FF1493', ghost: a.ghost === true ? true : undefined })),
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
            symmetrySource: data.crystal.symmetrySource,
            symmetryCheck: data.crystal.symmetryCheck ? {
                declared: data.crystal.symmetryCheck.declared,
                symbol: data.crystal.symmetryCheck.symbol,
                settingLabel: data.crystal.symmetryCheck.settingLabel,
                number: data.crystal.symmetryCheck.number,
                fileOps: data.crystal.symmetryCheck.fileOps,
                tableOps: data.crystal.symmetryCheck.tableOps,
                fileAtoms: data.crystal.symmetryCheck.fileAtoms,
                tableAtoms: data.crystal.symmetryCheck.tableAtoms
            } : undefined,
            symmetryAlternative: data.crystal.symmetryAlternative ? {
                ops: data.crystal.symmetryAlternative.ops,
                symbol: data.crystal.symmetryAlternative.symbol,
                settingLabel: data.crystal.symmetryAlternative.settingLabel,
                baseAtoms: data.crystal.symmetryAlternative.baseAtoms.map(crystalBasePayload),
                asymAtoms: data.crystal.symmetryAlternative.asymAtoms.map(crystalAsymPayload)
            } : undefined,
            baseAtoms: data.crystal.baseAtoms.map(crystalBasePayload),
            baseBonds: data.crystal.baseBonds.map(b => ({
                atom1: b.atom1, atom2: b.atom2, order: b.order, shift: b.shift, valence: b.valence
            })),
            fileBonds: (data.crystal.fileBonds || []).map(b => ({
                atom1: b.atom1, atom2: b.atom2, order: b.order, shift: b.shift, valence: b.valence
            })),
            asymAtoms: (data.crystal.asymAtoms || []).map(crystalAsymPayload)
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

        const jsonData = JSON.stringify({ atoms: atomData, bonds: bondData, title: data.title, atomColors: atomColors, filePath: data.filePath || '', frames: framesData, gjfMeta: data.gjfMeta || null, charge: data.charge, multiplicity: data.multiplicity, atomGroups: atomGroupsData, crystal: crystalData, optSteps: optStepsData, normalModes: normalModesData, optSource: logSource || null, optThresholds: optThresholdsData, routes: routesData, vmd: vmdData, autoBeta: !!(autoBeta && !vmdData), hasExplicitBonds: !!data.hasExplicitBonds, warnings: data.warnings || [], mmParams: describeForceField() });

        // The bond perception module, as source. Both sides run this text —
        // the host calls the compiled function directly, the webview evaluates
        // the serialisation below — so the viewport's re-detection can never
        // drift away from what the parser produced.
        const bondPerceptionSource = BOND_PERCEPTION_SOURCE;

        // The bond-valence module, as source, for the same reason: the report
        // the page shows computes its sums with the module this host ships,
        // not a second copy kept in the template. It carries a reading aid
        // only — nothing it returns is written back to the structure.
        const bvsSource = BVS_SOURCE;

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
#mol-info{position:absolute;top:8px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.55);padding:6px 10px;border-radius:4px;z-index:25;pointer-events:none;line-height:1.6;max-width:min(300px,calc(100% - 16px));overflow:hidden}
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
/* The crystal card is a fixed 210 px column, the same width as the View and
   Edit Coordinates cards below it. Its width used to be shrink-to-fit, and its
   six action buttons are inline-level with width:100%, so intrinsic sizing
   measured all six side by side and stretched the card to their sum (about
   600 px) even though it renders them stacked. A fixed column keeps the card
   over its own controls, not over the structure. */
#crystal-panel{position:absolute;top:8px;right:8px;background:rgba(0,0,0,0.6);border:1px solid var(--vscode-panel-border,#444);border-radius:4px;padding:8px 10px 5px;font-size:11px;color:var(--vscode-editor-foreground,#ccc);z-index:25;display:none;width:210px;max-width:calc(100% - 16px);max-height:calc(100% - 16px);overflow-y:auto}
#crystal-panel h4{margin:0 0 6px 0;font-size:11px;font-weight:600}
#crystal-panel .bnd-row{display:flex;align-items:center;gap:4px;margin:2px 0}
#crystal-panel .bnd-row label{width:52px;text-align:right}
#crystal-panel .bnd-row input{width:48px;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 3px;font-size:11px;text-align:center}
/* a/b/c each carry their lower and upper bound on one line: two inputs fill
   the row around a dash, so the three axes cost three rows instead of six. */
#crystal-panel .bnd-axis label{width:14px;text-align:center;font-weight:600}
#crystal-panel .bnd-axis input{flex:1 1 0;width:auto;min-width:0}
#crystal-panel .bnd-dash{flex:0 0 auto;font-size:10px;color:var(--vscode-descriptionForeground,#999)}
#crystal-panel .bnd-row select{flex:1;min-width:0;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 2px;font-size:11px}
#crystal-panel .bnd-row select:disabled{opacity:.5;cursor:default}
#crystal-panel .bnd-occ{display:flex;align-items:center;gap:4px;margin:3px 0 3px 2px;font-size:10px;line-height:1.3;color:var(--vscode-descriptionForeground,#999);cursor:pointer;white-space:normal}
#crystal-panel .bnd-occ input{width:auto;margin:0;flex:none}
#crystal-panel .bnd-occ.off{opacity:.5;cursor:default}
#crystal-panel .bnd-btn{display:block;margin-top:4px;width:100%;background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border:none;border-radius:3px;padding:3px 6px;font-size:11px;cursor:pointer}
#crystal-panel .bnd-btn:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#crystal-panel .bnd-presets{display:flex;gap:4px;margin-top:4px}
#crystal-panel .bnd-mini{flex:1;background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff);border:1px solid var(--vscode-panel-border,#444);border-radius:3px;padding:3px 4px;font-size:11px;cursor:pointer;white-space:nowrap;min-width:0}
#crystal-panel .bnd-mini:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#crystal-panel .bnd-cell{margin-top:4px;font-size:10px;line-height:1.4;color:var(--vscode-descriptionForeground,#999);white-space:pre-line}
#crystal-panel .bnd-sym{margin-top:6px;font-size:10px;line-height:1.4;color:#e8a33d;white-space:normal}
#crystal-panel .bnd-sym .bnd-btn{margin-top:4px;background:var(--vscode-button-secondaryBackground,#3a3d41)}
#crystal-panel .bnd-sym .bnd-btn:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#crystal-panel .hk-title{margin:6px 0 2px;font-weight:600;font-size:11px}
#crystal-panel .hk-row label{width:34px}
#crystal-panel .hk-row input{flex:1 1 0;width:auto;min-width:0;padding:1px 2px}
#crystal-panel .hk-val{margin:1px 0 0 34px;font-size:10px;line-height:1.35;color:var(--vscode-descriptionForeground,#999);white-space:normal}
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
#neigh-panel,#int-panel{display:none;position:absolute;bottom:40px;left:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:330px;max-width:calc(100% - 16px);box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#neigh-panel .ng-title,#int-panel .ng-title{font-size:12px;font-weight:bold;margin-bottom:4px}
#neigh-panel .ng-line,#int-panel .ng-line{display:flex;align-items:center;gap:6px;margin:4px 0;font-size:11px;min-height:20px}
#neigh-panel .ng-line label,#int-panel .ng-line label{white-space:nowrap;color:var(--vscode-descriptionForeground,#999)}
#neigh-panel .ng-line .ng-chk,#int-panel .ng-line .ng-chk{display:flex;align-items:center;gap:3px;color:var(--vscode-editor-foreground,#ccc)}
#neigh-panel input[type=number],#int-panel input[type=number]{width:64px;padding:3px 5px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:11px;box-sizing:border-box}
#neigh-panel .ng-anchor{font-weight:600;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
#neigh-panel .ng-summary,#int-panel .ng-summary{font-size:11px;color:var(--vscode-descriptionForeground,#999);margin:4px 0;min-height:14px}
#neigh-panel .ng-list,#int-panel .ng-list{max-height:220px;overflow-y:auto;border-top:1px solid rgba(128,128,128,0.25);margin-top:4px}
#neigh-panel .ng-item,#int-panel .ng-item{display:flex;align-items:baseline;gap:8px;padding:2px 3px;font-size:11px;cursor:pointer;white-space:nowrap}
#neigh-panel .ng-item:hover,#int-panel .ng-item:hover{background:rgba(255,255,255,0.08)}
#neigh-panel .ng-item.ng-off{opacity:0.6;cursor:default}
#neigh-panel .ng-name,#int-panel .ng-name{flex:1;overflow:hidden;text-overflow:ellipsis}
#neigh-panel .ng-d,#int-panel .ng-d{color:var(--vscode-descriptionForeground,#999);font-variant-numeric:tabular-nums}
#neigh-panel .ng-cell{color:#e8a33d;font-size:10px}
#neigh-panel .ng-empty,#int-panel .ng-empty{font-size:11px;color:var(--vscode-descriptionForeground,#999);padding:4px 3px}
#neigh-panel .ng-btns{display:flex;gap:4px;margin-top:6px}
#neigh-panel .ng-btns button{flex:1;padding:4px 6px;border-radius:3px;cursor:pointer;font-size:11px;border:1px solid var(--vscode-button-background,#0e639c);background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff)}
#neigh-panel .ng-btns button:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#neigh-panel .ng-hint,#int-panel .ng-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:4px;line-height:1.4}
#int-panel .ng-types{flex-wrap:wrap;gap:4px 8px}
#int-panel .ng-type{flex:none;font-size:9px;font-weight:600;padding:0 3px;border-radius:2px;color:#111}
#int-panel .ng-type-hb{background:#4fc3f7}
#int-panel .ng-type-xb{background:#ffb74d}
#int-panel .ng-type-cb{background:#ffd54f}
#int-panel .ng-type-pn{background:#ce93d8}
#int-panel .ng-type-pp{background:#81c784}
#int-panel .ng-type-chpi{background:#f06292}
#int-panel .ng-a{color:var(--vscode-descriptionForeground,#999);font-variant-numeric:tabular-nums}
#xrd-panel{display:none;position:absolute;bottom:40px;left:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:456px;max-width:calc(100% - 16px);max-height:calc(100% - 16px);overflow-y:auto;box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#xrd-panel .xr-title{display:flex;align-items:center;justify-content:space-between;font-size:12px;font-weight:bold;margin-bottom:4px}
#xrd-panel .xr-close{background:none;border:none;color:var(--vscode-editor-foreground,#ccc);font-size:15px;line-height:1;cursor:pointer;padding:0 2px}
#xrd-panel .xr-close:hover{color:#fff}
#xrd-panel .xr-line{display:flex;align-items:center;gap:6px;margin:4px 0;font-size:11px;min-height:20px}
#xrd-panel .xr-line label{white-space:nowrap;color:var(--vscode-descriptionForeground,#999)}
#xrd-panel select{flex:1;min-width:0;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 2px;font-size:11px}
#xrd-panel select.xr-narrow{flex:0 0 62px}
#xrd-panel .xr-summary{font-size:11px;color:var(--vscode-descriptionForeground,#999);margin:2px 0;min-height:14px}
#xrd-panel canvas{display:block;width:100%;height:150px;border:1px solid rgba(128,128,128,0.25);border-radius:3px;background:#0d1117;cursor:crosshair}
#xrd-panel .xr-note{font-size:10px;color:var(--vscode-descriptionForeground,#999);line-height:1.4;margin:3px 0}
#xrd-panel .xr-list{max-height:132px;overflow-y:auto;border-top:1px solid rgba(128,128,128,0.25);margin-top:3px;font-variant-numeric:tabular-nums}
#xrd-panel .xr-head,#xrd-panel .xr-row{display:flex;gap:6px;font-size:10px;padding:1px 2px;white-space:nowrap}
#xrd-panel .xr-head{color:var(--vscode-descriptionForeground,#999);border-bottom:1px solid rgba(128,128,128,0.18);margin-bottom:1px}
#xrd-panel .xr-c1{width:62px;text-align:right}
#xrd-panel .xr-c2{width:58px;text-align:right}
#xrd-panel .xr-c3{flex:1;text-align:right}
#xrd-panel .xr-c4{width:40px;text-align:right}
#xrd-panel .xr-c5{width:24px;text-align:right;color:#e8a33d}
#xrd-panel .xr-row{color:var(--vscode-editor-foreground,#ccc);cursor:pointer}
#xrd-panel .xr-row:hover{background:rgba(128,128,128,0.12)}
#xrd-panel .xr-row.sel{background:rgba(78,201,176,0.16)}
#xrd-panel .xr-empty{font-size:11px;color:var(--vscode-descriptionForeground,#999);padding:4px 2px}
#xrd-panel .xr-btns{display:flex;gap:4px;margin-top:6px}
#xrd-panel .xr-btns button{flex:1;padding:4px 6px;border-radius:3px;cursor:pointer;font-size:11px;border:1px solid var(--vscode-button-background,#0e639c);background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff)}
#xrd-panel .xr-btns button:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#xrd-panel .xr-btns button.xr-secondary{background:var(--vscode-button-secondaryBackground,#3a3d41);border-color:var(--vscode-panel-border,#444);color:var(--vscode-button-secondaryForeground,#fff)}
#xrd-panel .xr-btns button.xr-secondary:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#xrd-panel .xr-btns button:disabled{opacity:.5;cursor:default}
#xrd-panel .xr-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:4px;line-height:1.4}
/* The cavity report sits in the same floating slot as the powder panel and
   follows its list idiom: a scrollable body with clickable rows. */
#cav-panel{display:none;position:absolute;bottom:40px;left:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:392px;max-width:calc(100% - 16px);max-height:calc(100% - 16px);overflow-y:auto;box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#cav-panel .cv-title{display:flex;align-items:center;justify-content:space-between;font-size:12px;font-weight:bold;margin-bottom:4px}
#cav-panel .cv-close{background:none;border:none;color:var(--vscode-editor-foreground,#ccc);font-size:15px;line-height:1;cursor:pointer;padding:0 2px}
#cav-panel .cv-close:hover{color:#fff}
#cav-panel .cv-line{display:flex;align-items:center;gap:6px;margin:4px 0;font-size:11px}
#cav-panel .cv-line label{white-space:nowrap;color:var(--vscode-descriptionForeground,#999)}
#cav-panel .cv-line input{width:58px;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 3px;font-size:11px}
#cav-panel .cv-line button{padding:3px 8px;border-radius:3px;cursor:pointer;font-size:11px;border:1px solid var(--vscode-button-background,#0e639c);background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff)}
#cav-panel .cv-line button:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#cav-panel .cv-line button.cv-secondary{background:var(--vscode-button-secondaryBackground,#3a3d41);border-color:var(--vscode-panel-border,#444);color:var(--vscode-button-secondaryForeground,#fff);margin-left:auto}
#cav-panel .cv-line button.cv-secondary:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#cav-panel .cv-summary{font-size:11px;color:var(--vscode-editor-foreground,#ccc);line-height:1.5;margin:3px 0}
#cav-panel .cv-sub{font-size:10px;color:var(--vscode-descriptionForeground,#999)}
#cav-panel .cv-note{font-size:10px;color:#e8a33d;margin-top:2px;line-height:1.4;min-height:0}
#cav-panel .cv-list{max-height:176px;overflow-y:auto;border-top:1px solid rgba(128,128,128,0.25);margin-top:4px;font-variant-numeric:tabular-nums}
#cav-panel .cv-head,#cav-panel .cv-row{display:flex;gap:6px;font-size:10px;padding:1px 2px;white-space:nowrap}
#cav-panel .cv-head{color:var(--vscode-descriptionForeground,#999);border-bottom:1px solid rgba(128,128,128,0.18);margin-bottom:1px}
#cav-panel .cv-name{width:64px}
#cav-panel .cv-vol{width:64px;text-align:right}
#cav-panel .cv-c{flex:1;text-align:right}
#cav-panel .cv-r{width:62px;text-align:right}
#cav-panel .cv-row{color:var(--vscode-editor-foreground,#ccc);cursor:pointer}
#cav-panel .cv-row:hover{background:rgba(128,128,128,0.12)}
#cav-panel .cv-row.sel{background:rgba(177,140,255,0.18)}
#cav-panel .cv-empty{font-size:11px;color:var(--vscode-descriptionForeground,#999);padding:4px 2px}
#cav-panel .cv-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:4px;line-height:1.4}
/* The View-mode controls are two cards, not one: Auto Rotate and Edit
   Coordinates each carry their own collision box in layoutPanels(), so a tall
   neighbour (the Diff results panel above all) can push them apart instead of
   one of them ending up underneath something else. */
#view-panel,#bake-panel{display:block;position:absolute;bottom:40px;right:8px;z-index:50;background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:4px;padding:8px;width:210px;max-width:calc(100% - 16px);box-shadow:0 2px 8px rgba(0,0,0,0.3)}
#view-panel .vp-title,#bake-panel .vp-title{font-size:12px;font-weight:bold;margin-bottom:6px}
#view-panel .tbtn,#bake-panel .tbtn{width:100%;text-align:center}
#view-panel .vp-hint,#bake-panel .vp-hint{font-size:10px;color:var(--vscode-descriptionForeground,#999);margin-top:6px;line-height:1.4}
/* The drawing style is the second section of the View card, under Auto Rotate;
   its own title needs the top margin the shared .vp-title does not carry. */
#view-panel .vp-sub{margin-top:10px}
#view-panel select{width:100%;padding:3px 4px;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#3c3c3c);border-radius:3px;font-size:11px;box-sizing:border-box;cursor:pointer}
#view-panel .vp-check{display:flex;align-items:center;gap:4px;margin-top:6px;font-size:11px;color:var(--vscode-editor-foreground,#ccc);cursor:pointer}
#view-panel .vp-check input{margin:0}
#view-panel .vp-check.off{opacity:.45;cursor:default}
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
#modal.modal-wide{max-width:min(760px,calc(100vw - 24px))}
#modal .bvs-scroll{max-height:min(46vh,420px);overflow:auto;margin-top:4px}
#modal .bvs-table{width:100%;border-collapse:collapse;font-size:11px}
#modal .bvs-table th{text-align:left;font-weight:600;opacity:.75;padding:3px 6px;border-bottom:1px solid var(--vscode-panel-border,#444);position:sticky;top:0;background:var(--vscode-editor-background,#1e1e1e);white-space:nowrap}
#modal .bvs-table td{padding:2px 6px;border-bottom:1px solid rgba(255,255,255,0.06);white-space:nowrap}
#modal .bvs-table td.num{text-align:right;font-variant-numeric:tabular-nums}
#modal .bvs-note{font-size:11px;color:var(--vscode-descriptionForeground,#999);margin-top:8px;line-height:1.45}
#modal .bvs-warn{font-size:11px;color:#e8a33d;margin-top:4px;line-height:1.45}
#modal .bvs-meh{color:#e8a33d}
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
/* One header row, one close control and one chip for the panels that used to
   have neither: the x folds the panel into its chip, the chip puts it back,
   and layoutPanels() gives both a collision box, so a panel can never cover
   another panel or the chip that belongs to it. */
.pn-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:6px}
#crystal-panel .pn-head h4{margin:0;font-size:11px;font-weight:600}
#view-panel .pn-head .vp-title,#bake-panel .pn-head .vp-title,
#select-panel .pn-head .sp-title,#neigh-panel .pn-head .ng-title,
#int-panel .pn-head .ng-title,
#rotate-panel .pn-head .rp-title{margin-bottom:0}
.pn-close{background:none;border:none;color:inherit;font:inherit;font-size:15px;line-height:1;cursor:pointer;opacity:0.7;padding:0 2px}
.pn-close:hover{opacity:1;color:var(--vscode-errorForeground,#f66)}
.pn-chip{display:none;position:absolute;z-index:52;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto;white-space:nowrap;max-width:calc(100% - 16px);overflow:hidden;text-overflow:ellipsis}
.pn-chip.show{display:block}
.pn-chip:hover{background:var(--vscode-button-background,#0e639c)}
#crystal-reopen,#egg-reopen{top:8px;right:8px}
#view-reopen,#bake-reopen{right:8px;bottom:40px}
#select-reopen,#neigh-reopen,#int-reopen,#rotate-reopen,#xrd-reopen,#cav-reopen{left:8px;bottom:40px}
/* A panel whose height the layout had to cap scrolls inside itself instead of
   spilling over its neighbours. */
#view-panel,#bake-panel,#select-panel,#neigh-panel,#int-panel,#rotate-panel{overflow:auto}
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
<button class="tbtn" data-mode="neighbours">Neighbours</button>
<button class="tbtn" data-mode="interactions">Interactions</button>
<button class="tbtn" id="diff-btn">Diff</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="selectAtoms">Select Atoms</button>
<button class="tbtn" data-mode="boxSelect">Box Select</button>
<div class="tsep"></div>
<button class="tbtn" id="beta-btn" style="display:none">Color by Beta</button>
<button class="tbtn" id="egg-btn" style="display:none">Eggs</button>
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
<div id="container"><div id="loading">Loading 3D Viewer...</div><div id="mol-info"></div><div id="mm-panel"></div><div id="mm-reopen">MM Optimize</div><div id="egg-panel"></div><div id="diff-label"></div><div id="diff-label-right"></div><div id="diff-panel"></div><div id="diff-reopen">📊 Show Results</div><div id="axes-indicator"></div><div id="crystal-panel"><div class="pn-head"><h4>Supercell Bounds</h4><button class="pn-close" id="crystal-close" type="button" title="Minimize the crystal panel">&times;</button></div><div class="bnd-row bnd-axis"><label for="bnd-a-min">a</label><input type="number" id="bnd-a-min" value="0" step="0.1" title="lower bound along a" aria-label="lower bound along a"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-a-max" value="1" step="0.1" title="upper bound along a" aria-label="upper bound along a"></div><div class="bnd-row bnd-axis"><label for="bnd-b-min">b</label><input type="number" id="bnd-b-min" value="0" step="0.1" title="lower bound along b" aria-label="lower bound along b"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-b-max" value="1" step="0.1" title="upper bound along b" aria-label="upper bound along b"></div><div class="bnd-row bnd-axis"><label for="bnd-c-min">c</label><input type="number" id="bnd-c-min" value="0" step="0.1" title="lower bound along c" aria-label="lower bound along c"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-c-max" value="1" step="0.1" title="upper bound along c" aria-label="upper bound along c"></div><div class="bnd-presets"><button class="bnd-mini" id="bnd-p2" title="Eight unit cells, 0..2 along every axis">2×2×2</button><button class="bnd-mini" id="bnd-p3" title="Twenty-seven unit cells, 0..3 along every axis">3×3×3</button><button class="bnd-mini" id="bnd-p4" title="Sixty-four unit cells, 0..4 along every axis">4×4×4</button></div><div class="bnd-presets"><button class="bnd-mini" id="bnd-reset" title="Back to the single unit cell, 0..1 along every axis">Reset</button><button class="bnd-mini" id="bnd-origin" title="Keep the spans, move the low corner onto the lattice origin">Align Origin</button></div><div class="bnd-cell" id="bnd-cell-line"></div><div class="bnd-sym" id="bnd-sym-note" style="display:none"></div><div class="bnd-row"><label>Disorder</label><select id="bnd-disorder" title="Show the ordered part plus one disorder group; All shows every alternative the cell lists"><option value="">All</option></select></div><label class="bnd-occ off" id="bnd-occ-label" for="bnd-occ-normalize" title="When one disorder group is selected, count its sites at occupancy 1 in the formula, density and saved files, so the cell describes the ordered variant on screen. All keeps every occupancy the file states."><input type="checkbox" id="bnd-occ-normalize" checked disabled>Selected group fully occupied</label><div class="hk-title">Lattice Planes (h k l)</div><div class="bnd-row hk-row"><label>(h k l)</label><input type="number" id="hk-h" value="0" step="1" title="h"><input type="number" id="hk-k" value="0" step="1" title="k"><input type="number" id="hk-l" value="1" step="1" title="l"></div><div class="hk-val" id="hk-d1"></div><div class="bnd-row hk-row"><label>2nd</label><input type="number" id="hk2-h" value="1" step="1" title="h of the second plane"><input type="number" id="hk2-k" value="0" step="1" title="k of the second plane"><input type="number" id="hk2-l" value="0" step="1" title="l of the second plane"></div><div class="hk-val" id="hk-d2"></div><button id="hk-show" class="bnd-btn" title="Draw the (hkl) plane through the centre of the box">Show Plane</button><button id="xrd-open" class="bnd-btn" title="Calculate the powder X-ray diffraction pattern of the unit cell">Powder XRD Pattern...</button><button id="bvs-open" class="bnd-btn" title="Bond-valence sums per site: how well the coordination around each cation and anion matches its formal valence (read-only)">Bond Valence (BVS)...</button><button id="cav-open" class="bnd-btn" title="Probe the unit cell for closed cavities and open channels and report their volume (read-only)">Cavities...</button><button id="bnd-edit-cell" class="bnd-btn" title="Change the cell parameters; every site keeps its fractional coordinates">Edit Cell...</button><button id="bnd-remove-disorder" class="bnd-btn">Remove Disorder &lt;0.5</button></div><div id="view-panel"><div class="pn-head"><span class="vp-title">Auto Rotate</span><button class="pn-close" id="view-close" type="button" title="Minimize the View controls">&times;</button></div><button class="tbtn" id="spin-btn">Off</button><div class="vp-hint">Drag and release to spin, or hold an arrow key &mdash; the longer the hold, the faster. The next left press stops it.</div><div class="vp-title vp-sub">Display</div><select id="style-select" title="How the structure is drawn: ball and stick (the default), van der Waals spheres, sticks without the atom spheres, or coordination polyhedra for a crystal"><option value="ball">Ball &amp; Stick</option><option value="spacefill">Space-filling</option><option value="stick">Stick</option>${(function(){try{return JSON.parse(jsonData).crystal?'<option value="polyhedra">Polyhedra</option>':'';}catch(e){return '';}}())}</select><label class="vp-check off" id="ortep-line" title="Draw atoms that carry anisotropic displacement parameters as 50% probability ellipsoids"><input type="checkbox" id="ortep-chk" disabled> Thermal ellipsoids (50%)</label></div><div id="bake-panel"><div class="pn-head"><span class="vp-title">Edit Coordinates</span><button class="pn-close" id="bake-close" type="button" title="Minimize the Edit Coordinates card">&times;</button></div><div class="vp-hint" id="vp-bake-hint">Hold <b>Alt</b> and drag: left rotates the structure, right/middle moves it. The atom coordinates change, so Save writes the new geometry.</div></div><div id="select-panel"><div class="pn-head"><span class="sp-title">Select Atoms</span><button class="pn-close" id="select-close" type="button" title="Minimize the selection panel">&times;</button></div><input type="text" id="sel-panel-input" placeholder="e.g. 1 3-5 C H 8"><div class="sp-btns"><button class="mbtn mbtn-ok" id="sel-panel-ok">Select</button><button class="mbtn mbtn-cancel" id="sel-panel-clear">Clear</button></div></div><div id="neigh-panel"><div class="pn-head"><span class="ng-title">Neighbours</span><button class="pn-close" id="neigh-close" type="button" title="Minimize the neighbour list">&times;</button></div><div class="ng-line"><label for="ng-radius">Radius (A)</label><input type="number" id="ng-radius" value="3.5" min="0.1" max="12" step="0.1"><label class="ng-chk" for="ng-skip"><input type="checkbox" id="ng-skip" checked>skip bonded</label></div><div class="ng-line"><label>Anchor</label><span class="ng-anchor" id="ng-anchor">click an atom</span></div><div class="ng-summary" id="ng-summary"></div><div class="ng-list" id="ng-list"></div><div class="ng-btns"><button id="ng-select">Select listed</button></div><div class="ng-hint" id="ng-hint">Click a row to select that atom. Crystal images outside the drawn cell are listed with their cell offset.</div></div><div id="int-panel"><div class="pn-head"><span class="ng-title">Non-covalent Interactions</span><button class="pn-close" id="int-close" type="button" title="Minimize the interaction list">&times;</button></div><div class="ng-line ng-types"><label class="ng-chk" for="int-hb"><input type="checkbox" id="int-hb" checked>H-bond</label><label class="ng-chk" for="int-xb"><input type="checkbox" id="int-xb" checked>halogen</label><label class="ng-chk" for="int-cb"><input type="checkbox" id="int-cb" checked>chalcogen</label><label class="ng-chk" for="int-pn"><input type="checkbox" id="int-pn" checked>pnictogen</label><label class="ng-chk" for="int-pp"><input type="checkbox" id="int-pp" checked>&pi;-&pi;</label><label class="ng-chk" for="int-chpi"><input type="checkbox" id="int-chpi" checked>CH-&pi;</label></div><div class="ng-line"><label class="ng-chk" for="int-show"><input type="checkbox" id="int-show" checked>draw dashed lines</label></div><div class="ng-summary" id="int-summary"></div><div class="ng-list" id="int-list"></div><div class="ng-hint" id="int-hint">H-bond: d(H&middot;&middot;&middot;A) &lt; 2.5 &Aring; and D-H&middot;&middot;&middot;A &gt; 120&deg;. Halogen, chalcogen and pnictogen bonds: d(X&middot;&middot;&middot;A) &lt; vdW sum &minus; 0.1 &Aring; and D-X&middot;&middot;&middot;A &gt; 130&deg;; acceptors N, O, S, F, Cl, Br, I. &pi;-&pi;: aromatic 5&ndash;7 rings with a plane angle &le; 30&deg;, centroid distance &le; 4.5 &Aring; and a 2.8&ndash;3.8 &Aring; face gap (parallel/offset), or a &ge; 60&deg; T-shape with one ring edge 3.0&ndash;5.5 &Aring; over the other face. CH-&pi;: H&middot;&middot;&middot;centroid &lt; 3.0 &Aring;, C-H&middot;&middot;&middot;centroid &gt; 120&deg; and the H at least 1.0 &Aring; above the ring face.</div></div><div id="xrd-panel"><div class="xr-title"><span>Powder XRD</span><button class="xr-close" id="xrd-close" title="Close (Esc)">&times;</button></div><div class="xr-line"><label for="xrd-rad">Radiation</label><select id="xrd-rad" title="The X-ray wavelength the pattern is calculated for"></select><label for="xrd-2max">2&theta; max</label><select id="xrd-2max" class="xr-narrow" title="The highest diffraction angle in the pattern"></select></div><div class="xr-summary" id="xrd-summary"></div><canvas id="xrd-canvas" width="440" height="150"></canvas><div class="xr-note" id="xrd-peak-info"></div><div class="xr-note" id="xrd-exp-info"></div><div class="xr-note" id="xrd-note"></div><div class="xr-list" id="xrd-list"></div><div class="xr-btns"><button id="xrd-refresh" class="xr-secondary" title="Recalculate after the structure changed">Recalculate</button><button id="xrd-export" title="Save the drawn curve as a two-column text file">Export .xy</button></div><div class="xr-btns"><button id="xrd-import" class="xr-secondary" title="Overlay a measured two-column 2theta/intensity pattern">Import experimental .xy...</button><button id="xrd-clear-exp" class="xr-secondary" disabled title="Remove the experimental curve">Clear exp.</button></div><input type="file" id="xrd-file" accept=".xy,.txt,.csv,text/plain" style="display:none"><div class="xr-hint">Pattern of the unit cell on screen (the supercell bounds do not change it) with Cromer-Mann form factors, Lorenz-polarisation and a Caglioti Gaussian profile. Click a calculated peak to read its position against the experimental curve.</div></div><div id="cav-panel"><div class="cv-title"><span>Voids and Cavities</span><button class="cv-close" id="cav-close" title="Close (Esc)">&times;</button></div><div class="cv-line"><label for="cav-probe">Probe radius</label><input type="number" id="cav-probe" value="1.2" min="0" max="20" step="0.1" title="Radius of the spherical probe rolled through the cell"><span style="font-size:10px;opacity:.7">&Aring;</span><button id="cav-run" title="Scan the unit cell with this probe">Analyse</button><button id="cav-clear" class="cv-secondary" title="Remove the drawn surface">Clear</button></div><div class="cv-summary" id="cav-summary">Not analysed yet.</div><div class="cv-note" id="cav-note"></div><div class="cv-list" id="cav-list"></div><div class="cv-hint">A point is free when a probe of this radius fits there without touching an atom (van der Waals radii, periodic images included). The free space is split into components: one whose walk winds through a cell face is an open channel, otherwise it is a closed cavity. Click a row to draw that void as a translucent surface. Volumes are per unit cell, so a supercell shows the same numbers; pockets under 0.5 &Aring;&sup3; are ignored. r_free is the largest sphere that fits at the void's deepest point.</div></div><div id="rotate-panel"><div class="pn-head"><span class="rp-title">Rotate Group Around Axis</span><button class="pn-close" id="rotate-close" type="button" title="Minimize the rotate panel">&times;</button></div><div class="rp-field"><label>Group atoms (click atoms or type, e.g. H581 O380 C44 ...):</label><input type="text" id="rp-group-input" placeholder="e.g. H581 O380 C44 C43 C42 ..."></div><div class="rp-field"><label>Axis atoms (2+ atoms, single-bonded chain):</label><input type="text" id="rp-axis-input" placeholder="e.g. C32 C33 C43 C53  or  C33 C43"></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-apply">Apply &amp; Rotate</button><button class="mbtn mbtn-cancel" id="rp-clear">Clear</button></div><div class="rp-error" id="rp-error"></div><div id="rp-slider-wrap" style="display:none"><div class="rp-slider-row"><label>Angle:</label><input type="range" id="rp-slider" min="0" max="360" step="0.5" value="0"><input type="number" id="rp-angle-input" min="0" max="360" step="0.5" value="0" title="Type an angle (0-360) and press Enter"><span class="rp-angle-val" id="rp-angle-val">0°</span></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-done">Done</button><button class="mbtn mbtn-danger" id="rp-reset">Reset</button></div></div><div class="rp-hint">Axis atoms are kept fixed; the rest of the selected group rotates around the line through the axis atoms (0–360°).</div></div><div id="opt-panel"></div><div id="opt-reopen">📈 Convergence</div><div id="freq-panel"></div><div id="freq-reopen">🎵 Modes</div><div id="route-panel"></div><div id="route-reopen">⚙️ Route</div><div id="crystal-reopen" class="pn-chip" title="Open the crystal panel">Crystal</div><div id="egg-reopen" class="pn-chip" title="Open the egg list">Eggs</div><div id="view-reopen" class="pn-chip" title="Open the View controls">View</div><div id="bake-reopen" class="pn-chip" title="Open the Edit Coordinates card">Coordinates</div><div id="select-reopen" class="pn-chip" title="Open the Select Atoms panel">Select Atoms</div><div id="neigh-reopen" class="pn-chip" title="Open the Neighbours panel">Neighbours</div><div id="int-reopen" class="pn-chip" title="Open the Non-covalent Interactions panel">Interactions</div><div id="rotate-reopen" class="pn-chip" title="Open the Rotate Group panel">Rotate Group</div><div id="xrd-reopen" class="pn-chip" title="Open the powder XRD panel">Powder XRD</div><div id="cav-reopen" class="pn-chip" title="Open the cavity and void analysis">Cavities</div></div>
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
var BVS=${bvsSource};
var AN={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Ne:10,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Sc:21,Ti:22,V:23,Cr:24,Mn:25,Fe:26,Co:27,Ni:28,Cu:29,Zn:30,Ga:31,Ge:32,As:33,Se:34,Br:35,Kr:36,Rb:37,Sr:38,Y:39,Zr:40,Nb:41,Mo:42,Tc:43,Ru:44,Rh:45,Pd:46,Ag:47,Cd:48,In:49,Sn:50,Sb:51,Te:52,I:53,Xe:54,Cs:55,Ba:56,La:57,Ce:58,Pr:59,Nd:60,Pm:61,Sm:62,Eu:63,Gd:64,Tb:65,Dy:66,Ho:67,Er:68,Tm:69,Yb:70,Lu:71,Hf:72,Ta:73,W:74,Re:75,Os:76,Ir:77,Pt:78,Au:79,Hg:80,Tl:81,Pb:82,Bi:83,Po:84,At:85,Rn:86};
// Standard atomic weights (IUPAC 2021, abridged), the one table the density
// line of the structure panel needs: a cell's occupancy-weighted mass over its
// volume. It sits beside the atomic numbers it mirrors, above every reader; an
// element the table does not know leaves the density line out rather than
// guessing a mass for it.
var MASS={H:1.008,He:4.0026,Li:6.94,Be:9.0122,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Ne:20.180,
    Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,Ar:39.948,K:39.098,Ca:40.078,
    Sc:44.956,Ti:47.867,V:50.942,Cr:51.996,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,
    Ga:69.723,Ge:72.630,As:74.922,Se:78.971,Br:79.904,Kr:83.798,Rb:85.468,Sr:87.62,Y:88.906,Zr:91.224,
    Nb:92.906,Mo:95.95,Tc:98,Ru:101.07,Rh:102.91,Pd:106.42,Ag:107.87,Cd:112.41,In:114.82,Sn:118.71,
    Sb:121.76,Te:127.60,I:126.90,Xe:131.29,Cs:132.91,Ba:137.33,La:138.91,Ce:140.12,Pr:140.91,Nd:144.24,
    Pm:145,Sm:150.36,Eu:151.96,Gd:157.25,Tb:158.93,Dy:162.50,Ho:164.93,Er:167.26,Tm:168.93,Yb:173.05,
    Lu:174.97,Hf:178.49,Ta:180.95,W:183.84,Re:186.21,Os:190.23,Ir:192.22,Pt:195.08,Au:196.97,Hg:200.59,
    Tl:204.38,Pb:207.2,Bi:208.98,Po:209,At:210,Rn:222};
// ===== Ghost atoms =====
// A ghost is a position without a nucleus, electrons or bonds: the marker a
// counterpoise calculation puts where the other fragment's atoms are. It keeps
// the element it stands for (the 'element' field), because that is what decides
// its basis functions in the two programs that have a ghost concept at all; a
// bare dummy point is element 'X'.
//
// One predicate decides everything below - drawing, formula, electron count,
// label, bond guards - and it is the same expression as isGhostAtom in
// src/types.ts (the page's script is inlined and cannot import it, so the two
// are kept identical by hand). Nothing in this page may test element 'X' or the
// ghost flag on its own: a structure loaded from a file carries either.
var GHOST_COLOR='#8FA8BF';   // grey-blue, deliberately unlike any element
var GHOST_RADIUS=0.28;       // sphere radius in A (a real atom is 0.3-1.0)
function isGhostAtom(a){return !!a&&(a.ghost===true||a.element==='X')}
/** One atom's label everywhere the page names an atom: C3, or C*3 / X5. */
function atomLabel(a,i){
    if(!a)return '?';
    var el=isGhostAtom(a)?(a.element==='X'?'X':a.element+'*'):a.element;
    return el+((typeof i==='number')?(i+1):'');
}
// ===== Display style =====
// A structure can be drawn four ways, and the choice belongs to the view, not
// to a mode: it is one variable here, above every block that reads it, and the
// toolbar switching tools never touches it - which is what makes a style chosen
// in View still be there in Bond Length, Move Atoms, MM Optimize or Diff, and
// what makes a frame, molecule or undo switch keep it too.
//
//   'ball'      the drawing this viewer has always produced: a sphere of half
//               the covalent radius per atom, cylinders between them.
//   'spacefill' each atom as its van der Waals sphere, and no sticks at all -
//               the model space-filling/CPK means.
//   'stick'     the sticks stay, and every sphere shrinks to the thickness of
//               the sticks that meet at it, so the model reads as tubes with no
//               beads (and no notch at a bend, because the joint fills it).
//   'polyhedra' ball and stick plus the coordination hull of every metal or
//               metalloid centre (see the polyhedra section): translucent
//               faces in the centre's colour with their edges, merged into one
//               geometry. A centre with fewer than four neighbours, a planar
//               shell, or a structure whose centres bond only to other centres
//               keeps the ball-and-stick picture, so the style never draws a
//               shape that is not there.
//
// Only two places in the page read it - the atom radius of rebuildScene (and of
// the comparison's right half) and the piece builder createBond - so one part
// of the picture cannot be drawn in one style and another part in another.
var displayStyle='ball';
var DISPLAY_STYLES={ball:1,spacefill:1,stick:1,polyhedra:1};
// van der Waals radii: the same table the host resolves a script's Surf probe
// and the VDW representation against (src/parsers/vmdVisualization.ts). It is
// copied rather than injected because the page's HTML is a template literal
// whose substitution list the test harnesses pin; verify_display_style.js
// compares the two copies element by element, so a drift fails a test rather
// than silently drawing atoms at the wrong size.
var VDW_R={H:1.00,He:1.40,Li:1.82,Be:1.53,B:1.92,C:1.70,N:1.55,O:1.52,
    F:1.47,Ne:1.54,Na:2.27,Mg:1.73,Al:1.84,Si:2.10,P:1.80,S:1.80,Cl:1.75,
    Ar:1.88,K:2.75,Ca:2.31,Sc:2.30,Ti:2.15,V:2.05,Cr:2.05,Mn:2.05,Fe:2.05,
    Co:2.00,Ni:2.00,Cu:2.00,Zn:2.10,Ga:1.87,Ge:2.11,As:1.85,Se:1.90,Br:1.85,
    Kr:2.02,Rb:3.03,Sr:2.49,Y:2.40,Zr:2.30,Nb:2.15,Mo:2.10,Ru:2.05,Rh:2.00,
    Pd:2.05,Ag:2.10,Cd:2.20,In:2.20,Sn:2.25,Sb:2.20,Te:2.06,I:1.98,Xe:2.16,
    Cs:3.43,Ba:2.68,La:2.50,Ce:2.48,Hf:2.25,Ta:2.20,W:2.10,Re:2.05,Os:2.00,
    Ir:2.00,Pt:2.05,Au:2.10,Hg:2.05,Tl:2.20,Pb:2.30,Bi:2.30,Th:2.40,U:2.30};
var DEFAULT_VDW_R=1.50;
// The viewer's own bond radius: the radius of a cylinder whose bond carries
// none, and the reference a script's bond radius is rebased against on the host
// (VIEWER_BOND_RADIUS there, the same 0.12). It is also the joint a stick model
// puts at an atom that has nothing thicker to match.
var BOND_R=0.12;
var STICK_JOINT_R=BOND_R;
// Radius of one bond's cylinders: what a script asked for, or the default. The
// stick model's joints read the same number, so a joint is exactly as thick as
// the sticks it caps and a script that draws thin lines does not grow beads.
function bondPieceRadius(b){
    return (b&&typeof b.vmdR==='number')?b.vmdR:BOND_R;
}
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
// The most atoms a supercell may build. The bounds inputs alone can ask for
// 11x11x11 cells of a large structure, which would be minutes of mesh building
// for a picture nobody asked for — the same ceiling the clipboard refuses
// fragments over, applied before the first sphere is made.
var MAX_CRYSTAL_ATOMS=30000;
// The most cell repeats the periodic neighbour search will enumerate along one
// axis. The count it actually needs comes from the cell's own shape (a repeat
// is only worth trying while it is narrower than the cutoff being searched), so
// a real crystal asks for one or two; the cap exists so a degenerate cell —
// two lattice vectors nearly parallel, say — cannot turn a bond search into a
// few hundred thousand distance tests per atom pair.
var CRYSTAL_MAX_IMAGE_REPEAT=4;
// How far outside the bounds a whole molecule's centre may sit and still be
// drawn (Angstrom, converted per axis with the lattice vector's length). Zero
// is the strict reading: only images centred inside the boundary are drawn.
var CRYSTAL_COPY_MARGIN=0;
// A bond the page itself added (Add Atom) carries no periodic image shift;
// this stands in wherever one is read.
var ZERO_SHIFT=[0,0,0];
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
// The lattice matrix the page works with, from the six cell parameters, in the
// orientation every parser builds it (a along x, b in the xy plane). Edit Cell
// is the one place in the page that changes the cell, so it uses this formula
// and the same six inputs the parsers read.
function crystalLatticeFromParams(a,b,c,alpha,beta,gamma){
    var rad=Math.PI/180;
    var ca=Math.cos(alpha*rad),cb=Math.cos(beta*rad),cg=Math.cos(gamma*rad),sg=Math.sin(gamma*rad);
    var av=[a,0,0];
    var bv=[b*cg,b*sg,0];
    var cx=c*cb;
    var cy=c*(ca-cb*cg)/sg;
    var cz2=1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg;
    var cz=c*Math.sqrt(Math.max(0,cz2))/sg;
    return [av,bv,[cx,cy,cz]];
}
function crystalFracToCart(f,lv){
    return [f[0]*lv[0][0]+f[1]*lv[1][0]+f[2]*lv[2][0],
            f[0]*lv[0][1]+f[1]*lv[1][1]+f[2]*lv[2][1],
            f[0]*lv[0][2]+f[1]*lv[1][2]+f[2]*lv[2][2]];
}
// The six bounds inputs (three a/b/c rows of lower-upper pairs) show the
// boundary in force. It lives outside the
// crystal panel's own block because an undo restores the boundary together with
// the base cell, and has to put the inputs back with it.
function syncBndInputs(){
    [['bnd-a-min','aMin'],['bnd-a-max','aMax'],['bnd-b-min','bMin'],
     ['bnd-b-max','bMax'],['bnd-c-min','cMin'],['bnd-c-max','cMax']].forEach(function(pr){
        var el=document.getElementById(pr[0]);
        if(el)el.value=boundary[pr[1]];
    });
}
// The cell parameters in force, in the supercell panel. Only a display: the
// dialog reads them back from CRY, so a stale line can never be edited from.
function syncCellPanel(){
    var el=document.getElementById('bnd-cell-line');
    if(!el||!CRY)return;
    el.textContent='a '+CRY.a.toFixed(3)+'  b '+CRY.b.toFixed(3)+'  c '+CRY.c.toFixed(3)+
        '\\nα '+CRY.alpha.toFixed(2)+'°  β '+CRY.beta.toFixed(2)+'°  γ '+CRY.gamma.toFixed(2)+'°';
}
// The parser's own check of the file's space-group name against the operations
// the file lists. One line states both expansions; the switch below is offered
// only while the parser kept the name's expansion, and it is one-time (undo
// restores the file's). It follows every rebuild, so an undo that puts the
// file's symmetry back also puts this note back.
function syncSymmetryNote(){
    var el=document.getElementById('bnd-sym-note');
    if(!el)return;
    var check=CRY&&CRY.symmetryCheck;
    if(!check){el.style.display='none';el.innerHTML='';return}
    var sym=escapeHtml(check.symbol)+(check.settingLabel?' ('+escapeHtml(check.settingLabel)+')':'');
    var html='Symmetry: '+sym+' gives '+check.tableOps.length+' op'+(check.tableOps.length===1?'':'s')+
        ' &rarr; '+check.tableAtoms+' atoms; the file lists '+check.fileOps.length+' op'+(check.fileOps.length===1?'':'s')+
        ' &rarr; '+check.fileAtoms+' atoms. The file\\'s operations are in use.';
    if(CRY.symmetryAlternative){
        html+='<button class="bnd-btn" id="bnd-sym-switch" title="Re-expand the cell under the name\\'s operations; Undo restores the file\\'s expansion">Use the name\\'s operations</button>';
    }
    el.innerHTML=html;
    el.style.display='block';
    var btn=document.getElementById('bnd-sym-switch');
    if(btn)btn.addEventListener('click',crystalUseDeclaredSymmetry);
}
// The one-time switch to the declared name's operations. The parser expanded
// the same asymmetric unit under that name already, so the new cell costs no
// re-reading; its bonds, however, refer to atom indices that no longer exist,
// so perception runs again on the new cell. The step is undoable: pushUndo
// carries the symmetry state (operations, sites, file bonds, this very check),
// which no ordinary edit changes and an ordinary snapshot would not restore.
function crystalUseDeclaredSymmetry(){
    if(!CRY||!CRY.symmetryAlternative)return;
    var alt=CRY.symmetryAlternative;
    if(mmRun)mmCancelRun();
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||moveKeyActive||fragRotActive)stopMoveSession();
    if(bakeActive)endBakeSession();
    pushUndo({
        ops:CRY.symmetryOps.slice(),
        source:CRY.symmetrySource,
        asymAtoms:JSON.parse(JSON.stringify(CRY.asymAtoms||[])),
        fileBonds:CRY.fileBonds?JSON.parse(JSON.stringify(CRY.fileBonds)):undefined,
        bondsDetected:!!CRY.bondsDetected,
        check:CRY.symmetryCheck,
        alternative:CRY.symmetryAlternative
    });
    CRY.baseAtoms=alt.baseAtoms;
    CRY.asymAtoms=alt.asymAtoms;
    CRY.symmetryOps=alt.ops;
    CRY.symmetrySource='builtin';
    CRY.symmetryCheck=null;
    CRY.symmetryAlternative=null;
    CRY.fileBonds=undefined;
    CRY.baseBonds=[];
    CRY.bondsDetected=false;
    CRY.bondSource=null;
    // The old display indices name atoms that are gone, so the selection is
    // dropped rather than remapped to whatever now sits at those numbers.
    resetSelection();
    rebuildCrystal();
    syncCellPanel();
    rebuildScene(true);
    layoutPanels();
}
function crystalCellEstimate(nb){
    var na=Math.ceil(nb.aMax)-Math.floor(nb.aMin);
    var nbv=Math.ceil(nb.bMax)-Math.floor(nb.bMin);
    var nc=Math.ceil(nb.cMax)-Math.floor(nb.cMin);
    if(na<1)na=1;
    if(nbv<1)nbv=1;
    if(nc<1)nc=1;
    return na*nbv*nc;
}
// The disorder groups a cell states, in the file's own spelling. An atom with
// no group belongs to the ordered part and is never a group of its own.
function crystalDisorderGroups(){
    var all=(CRY&&CRY.baseAtoms)||[];
    var seen={},out=[];
    for(var i=0;i<all.length;i++){
        var g=all[i].disorderGroup;
        if(g&&!seen[g]){seen[g]=true;out.push(g)}
    }
    out.sort(function(x,y){return x<y?-1:x>y?1:0});
    return out;
}
// The base sites the view shows, saves and counts: every ungrouped site plus
// the alternative the panel names. The filter is a view state only — CRY keeps
// every site at its own index, so edits still commit to the site they name and
// switching back to All (or undoing) loses nothing.
function crystalVisibleBaseAtoms(){
    var all=(CRY&&CRY.baseAtoms)||[];
    var g=CRY&&CRY.disorderFilter;
    if(!g)return all;
    var out=[];
    for(var i=0;i<all.length;i++){
        var a=all[i];
        if(!a.disorderGroup||a.disorderGroup===g)out.push(a);
    }
    return out;
}
// A disorder selection shows one alternative of a partially occupied crystal.
// That alternative is a valid ordered variant, so with this option on its sites
// enter the statistics and the saved files at full occupancy: the cell's
// formula and density then describe the variant on screen instead of the
// half-filled average of every alternative. An ungrouped partial site keeps
// the occupancy the file states, and the All view never changes.
var disorderOccNormalize=true;
function crystalOccOf(ba){
    if(CRY&&CRY.disorderFilter&&disorderOccNormalize&&ba&&ba.disorderGroup===CRY.disorderFilter)return 1;
    return ba&&ba.occupancy!=null?ba.occupancy:1;
}
// The panel's All / group select, rebuilt from the cell in force. It follows
// every rebuild — an edit, an undo, a removed atom — so a group that no longer
// exists leaves the list and the view falls back to All.
function syncDisorderPanel(){
    var sel=document.getElementById('bnd-disorder');
    if(!sel||!CRY)return;
    var groups=crystalDisorderGroups();
    if(CRY.disorderFilter&&groups.indexOf(CRY.disorderFilter)<0){
        CRY.disorderFilter=null;
        // The group the title named is gone with the atoms that carried it;
        // All is the state on screen now, so neither clause is true any more.
        provUnmark('disorderGroup');
        provUnmark('occNormalize');
    }
    var want=CRY.disorderFilter||'';
    var opts='<option value="">All</option>';
    for(var i=0;i<groups.length;i++){
        var g=groups[i];
        opts+='<option value="'+escapeHtml(g)+'"'+(g===want?' selected':'')+'>Group '+escapeHtml(g)+'</option>';
    }
    sel.innerHTML=opts;
    sel.value=want;
    sel.disabled=groups.length===0;
    // The occupancy switch only means something while one group is shown; it
    // keeps its own state so going back to All and into a group again does not
    // silently flip how the variant is counted.
    var occChk=document.getElementById('bnd-occ-normalize');
    if(occChk){
        occChk.checked=disorderOccNormalize;
        occChk.disabled=!CRY.disorderFilter;
    }
    var occLabel=document.getElementById('bnd-occ-label');
    if(occLabel)occLabel.className='bnd-occ'+(CRY.disorderFilter?'':' off');
}
// Edit Cell: the six parameters are the one description of the box a reader
// can check against the file, so they are shown back as typed values. A site's
// fractional coordinates are the invariant across the edit — that is what
// makes a bond along a stretch with a, and only with a. The base cell, its
// bonds and the bounds all stay in their own terms; the picture, the info
// panel and every export read the new cell from CRY.
function applyCellEdit(a,b,c,al,be,ga){
    if(!CRY)return false;
    var vals=[a,b,c,al,be,ga];
    var bad=vals.some(function(v){return !isFinite(v)});
    if(!bad&&!(a>0&&b>0&&c>0))bad=true;
    if(!bad&&(al<=0||al>=180||be<=0||be>=180||ga<=0||ga>=180))bad=true;
    var lv=null;
    if(!bad){
        lv=crystalLatticeFromParams(a,b,c,al,be,ga);
        var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-
                lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+
                lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
        if(!(Math.abs(det)>1e-8))bad=true;
    }
    if(bad){
        showModal('<h3>Invalid Cell</h3><div class="current-val">Lengths must be positive, angles between 0 and 180 degrees, and the three vectors must enclose a nonzero volume.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
        document.getElementById('m-cancel').addEventListener('click',hideModal);
        return false;
    }
    // The picture is rebuilt from the base cell, so the same ceiling that
    // guards the bounds guards a changed cell too.
    var est=CRY.baseAtoms.length*crystalCellEstimate(boundary);
    if(est>MAX_CRYSTAL_ATOMS){
        showModal('<h3>Supercell Too Large</h3><div class="current-val">This cell would build about '+est+' atoms; the limit is '+MAX_CRYSTAL_ATOMS+'.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
        document.getElementById('m-cancel').addEventListener('click',hideModal);
        return false;
    }
    var inv=CRY_INV;
    if(!inv){
        showModal('<h3>Invalid Cell</h3><div class="current-val">The cell in force has no inverse, so its sites have no fractional coordinates to keep.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
        document.getElementById('m-cancel').addEventListener('click',hideModal);
        return false;
    }
    if(moveDragActive||moveKeyActive||fragRotActive)stopMoveSession();
    // Applying the cell that is already in force is not an edit: the sites are
    // put back where they are, and the title has nothing to add.
    var cellChanged=Math.abs(CRY.a-a)>1e-9||Math.abs(CRY.b-b)>1e-9||Math.abs(CRY.c-c)>1e-9||
        Math.abs(CRY.alpha-al)>1e-9||Math.abs(CRY.beta-be)>1e-9||Math.abs(CRY.gamma-ga)>1e-9;
    pushUndo();
    var base=CRY.baseAtoms;
    var fracs=[];
    var i;
    for(i=0;i<base.length;i++){
        var p=base[i];
        fracs.push([
            inv[0][0]*p.x+inv[0][1]*p.y+inv[0][2]*p.z,
            inv[1][0]*p.x+inv[1][1]*p.y+inv[1][2]*p.z,
            inv[2][0]*p.x+inv[2][1]*p.y+inv[2][2]*p.z
        ]);
    }
    CRY.a=a;CRY.b=b;CRY.c=c;CRY.alpha=al;CRY.beta=be;CRY.gamma=ga;
    CRY.latticeVectors=lv;
    computeCrystalInv();
    for(i=0;i<base.length;i++){
        var cart=crystalFracToCart(fracs[i],lv);
        base[i].x=cart[0];base[i].y=cart[1];base[i].z=cart[2];
    }
    // baseBonds keep their fractional shifts: the same two sites in the same
    // neighbouring cells are still joined by the same bond.
    rebuildCrystal();
    syncCellPanel();
    if(cellChanged){
        provSet('cellEdit',a.toFixed(3)+' x '+b.toFixed(3)+' x '+c.toFixed(3)+' A, '+
            al.toFixed(1)+'/'+be.toFixed(1)+'/'+ga.toFixed(1)+' deg');
    }
    rebuildScene(true);
    mmSelectionChanged();
    layoutPanels();
    return true;
}
function formulaFromAtoms(arr){
    var counts={};
    // A ghost contributes no atom to the formula: it has no nucleus, so the
    // composition of a counterpoise structure is the composition of its real
    // atoms (the ghosts are the *other* fragment, already counted there).
    arr.forEach(function(a){if(isGhostAtom(a))return;counts[a.element]=(counts[a.element]||0)+1});
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
    var nGhosts=0;
    MD.atoms.forEach(function(a){
        // A ghost has no nucleus and no electrons, so it adds neither to the
        // electron count nor to the formula below.
        if(isGhostAtom(a)){nGhosts++;return}
        nElectrons+=(AN[a.element]||0);
    });
    if(typeof chrg==='number')nElectrons-=chrg;
    var html='Formula: '+formulaFromAtoms(MD.atoms)+'<br>Atoms: '+nAtoms+
        (nGhosts>0?' ('+nGhosts+' ghost'+(nGhosts===1?'':'s')+')':'')+
        '<br>Charge: '+chrg+'<br>Electrons: '+nElectrons+'<br>Multiplicity: '+mult;
    if(nGhosts>0){
        html+='<br>Ghosts: '+nGhosts+' (no bonds, no basis, not in the formula)';
    }
    if(CRY){
        var sa=boundary.aMax-boundary.aMin,sb=boundary.bMax-boundary.bMin,sc=boundary.cMax-boundary.cMin;
        // A group view shows one ordered variant, so the cell's composition is
        // the composition of that variant; with the occupancy switch on its
        // sites are counted at full occupancy instead of the file's partial
        // value. All hands back the base array itself, so its line is unchanged
        // to the digit.
        var visBase=crystalVisibleBaseAtoms();
        var occNorm=!!(CRY.disorderFilter&&disorderOccNormalize);
        html+='<br>Unit Cell: '+formulaFromAtoms(visBase)+
            (occNorm?' (group '+escapeHtml(CRY.disorderFilter)+' fully occupied)':'');
        html+='<br>Cell: '+CRY.a.toFixed(3)+' × '+CRY.b.toFixed(3)+' × '+CRY.c.toFixed(3);
        html+='<br>Angles: '+CRY.alpha.toFixed(1)+'°, '+CRY.beta.toFixed(1)+'°, '+CRY.gamma.toFixed(1)+'°';
        if(CRY.spaceGroup)html+='<br>SG: '+escapeHtml(CRY.spaceGroup)+(CRY.symmetrySource==='builtin'?' (ops from built-in table)':'');
        html+='<br>Supercell: '+sa.toFixed(2)+'×'+sb.toFixed(2)+'×'+sc.toFixed(2);
        // The file's own bond table, when it states one: how many bonds it
        // lists against how many distance perception found. A row the table
        // states and perception does not is a disagreement worth showing; the
        // file's extra rows are not drawn (see mergeCrystalFileBonds).
        if(CRY.bondSource&&CRY.bondSource.file>0){
            html+='<br>Bonds: file '+CRY.bondSource.file+' / detected '+CRY.bondSource.detected;
            if(CRY.bondSource.fileOnly>0){
                html+=' <span style="color:#e8a33d">('+CRY.bondSource.fileOnly+
                    ' only in the file)</span>';
            }
        }
        // Crystallographic constants of the unit cell itself: its volume, the
        // number of whole molecules it holds (Z — the count rebuildCrystal's
        // unwrap pass produces; a cell whose only component is an infinite
        // network has no meaningful Z, and the line stays out rather than
        // showing a number that would be one), and the density the cell's
        // occupancy-weighted mass implies. Everything reads only what the page
        // already carries, so a harness that runs this function on a stub
        // crystal without a lattice simply sees none of these lines.
        if(CRY.latticeVectors){
            var crV=CRY.latticeVectors;
            var crDet=crV[0][0]*(crV[1][1]*crV[2][2]-crV[1][2]*crV[2][1])-crV[0][1]*(crV[1][0]*crV[2][2]-crV[1][2]*crV[2][0])+crV[0][2]*(crV[1][0]*crV[2][1]-crV[1][1]*crV[2][0]);
            var crVol=Math.abs(crDet);
            if(crVol>0){
                html+='<br>Volume: '+crVol.toFixed(2)+' Å³';
                if(CRY.moleculeCount>0)html+='<br>Z: '+CRY.moleculeCount+' (molecules)';
                if(typeof MASS!=='undefined'&&MASS){
                    var crMass=0,crMassOk=true;
                    for(var crI=0;crI<visBase.length;crI++){
                        var crM=MASS[visBase[crI].element];
                        if(crM==null){crMassOk=false;break}
                        crMass+=crM*crystalOccOf(visBase[crI]);
                    }
                    if(crMassOk&&crMass>0)html+='<br>Density: '+(crMass/(0.602214*crVol)).toFixed(3)+' g/cm³';
                }
            }
        }
        // G9: the asymmetric unit the file lists, with the operations that
        // placed its images. Collapsed by default, because a molecular crystal
        // with hundreds of sites would otherwise push the composition off the
        // panel. This is the same data the tooltip quotes for one atom, shown
        // as the whole site list.
        if(CRY.asymAtoms&&CRY.asymAtoms.length){
            var asu=CRY.asymAtoms;
            var asuRows='';
            for(var ai=0;ai<asu.length;ai++){
                var asuSite=asu[ai];
                var hit=asuSite.symmetryOps||[];
                var ops='';
                var shown=0;
                for(var hj=0;hj<hit.length&&shown<3;hj++){
                    if(ops!=='')ops+='; ';
                    ops+='#'+(hit[hj]+1)+' '+escapeHtml((CRY.symmetryOps&&CRY.symmetryOps[hit[hj]])||'');
                    shown++;
                }
                if(hit.length>shown)ops+=' +'+(hit.length-shown);
                if(ops==='')ops='-';
                asuRows+='<tr><td>'+escapeHtml(asuSite.label||('#'+(ai+1)))+'</td><td>'+escapeHtml(asuSite.element)+
                    '</td><td>'+(asuSite.occupancy!=null?asuSite.occupancy:1)+'</td><td>'+hit.length+'</td><td>'+ops+'</td></tr>';
            }
            html+='<details style="margin-top:5px"><summary style="cursor:pointer">Asymmetric unit: '+
                asu.length+' site'+(asu.length===1?'':'s')+'</summary>'+
                '<table style="font-size:11px;border-collapse:collapse;margin-top:3px">'+
                '<tr style="opacity:.75"><th align="left">site</th><th align="left">el</th><th align="left">occ</th><th align="left">mult</th><th align="left">operations</th></tr>'+
                asuRows+'</table></details>';
        }
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
    // What the parser had to work around ("read as P1", "unit cell incomplete"),
    // one line: the first warning verbatim plus how many others there are. It
    // sits last so it never pushes the composition off the panel, and the text
    // is escaped because a file's own words (a site label, a space-group name)
    // can contain anything. Nothing is added when the parser was happy, so an
    // ordinary structure reads exactly as it always did.
    var parseWarn=MD.warnings||[];
    if(parseWarn.length){
        html+='<br><span style="color:#e8a33d">⚠ '+escapeHtml(parseWarn[0])+
            (parseWarn.length>1?' (+'+(parseWarn.length-1)+' more)':'')+'</span>';
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
// without a script the plain default applies. A ghost has no element to take a
// radius from — it is drawn as the same small marker whatever it stands for, so
// the picture says "point, not atom" before any label is read.
function atomRadius(a){
    if(isGhostAtom(a))return typeof a.vmdScale==='number'?(GHOST_RADIUS*a.vmdScale):GHOST_RADIUS;
    var r=getR(a.element);
    return typeof a.vmdScale==='number'?r*a.vmdScale:r;
}
// Radius the *current display style* draws one atom with. 'ball' is exactly
// atomRadius, so the default drawing cannot drift from what it always was, and
// the value is returned unchanged - not multiplied by 1. A ghost is the same
// small marker in every style: it is a point without a nucleus, so neither a
// van der Waals sphere nor a stick joint is a meaningful size for it.
//
// A space-filling sphere is the element's van der Waals radius, which is the one
// length in the picture that is a property of the atom rather than of the
// drawing. A script's per-atom sphere scale is deliberately *not* applied here:
// that scale is a ratio inside its ball-and-stick recipe (and a zero means "this
// representation draws no sphere for it"), so applying it would keep a
// space-filling model relative to a drawing it is replacing, and would let a
// subset selection empty the picture - the spheres and the sticks are the two
// things the style decides, and the script's colours and surfaces still apply.
//
// A stick joint is jointR, the thickest stick meeting at that atom (see
// stickJointRadii); an atom with no sticks at all is drawn at the default stick
// thickness instead, so an isolated ion does not vanish.
function drawnRadius(a,jointR){
    // The polyhedra style is the ball-and-stick picture with hulls added, so
    // its atoms are drawn exactly as the default style draws them.
    if(displayStyle==='ball'||displayStyle==='polyhedra')return atomRadius(a);
    if(isGhostAtom(a))return atomRadius(a);
    if(displayStyle==='spacefill')return VDW_R[a.element]||DEFAULT_VDW_R;
    return jointR>0?jointR:STICK_JOINT_R;
}
// Thickness of the sticks meeting at each atom, in one pass over the bond list:
// the stick style's joint is exactly as thick as the bond it caps, so a script
// that draws thin lines (or a structure whose bonds are thin) does not end up
// with beads at its atoms. Bonds a script draws at radius 0 do not thicken
// anything - nothing is drawn there to fill.
function stickJointRadii(){
    var out=new Float32Array(MD.atoms.length);
    for(var i=0;i<MD.bonds.length;i++){
        var b=MD.bonds[i];
        var r=bondPieceRadius(b);
        if(!(r>0))continue;
        if(b.atom1<out.length&&r>out[b.atom1])out[b.atom1]=r;
        if(b.atom2<out.length&&r>out[b.atom2])out[b.atom2]=r;
    }
    return out;
}

// ===== Coordination polyhedra (the polyhedra display style) =====
// The style draws the coordination hull of every metal or metalloid centre:
// translucent faces in the centre's colour with their edges, all centres in
// one merged geometry, over the ball-and-stick picture the other styles draw.
// A centre whose shell is entirely made of other centres (a pure metal, an
// intermetallic) is left as ball and stick: its hull would put a shape on every
// atom of the structure and say nothing a bond does not. So is a centre with
// fewer than four neighbours, and one with a planar shell (square planar,
// trigonal, or a collinear detector accident) - there is no solid to draw, and
// a flat sheet would only hide the atoms inside it.
//
// The hull is the convex hull of the neighbour positions. Every triple of
// points proposes a plane; a triple is a hull face when every other point is
// behind its plane (within a tolerance scaled to the shell). Coplanar hull
// vertices are collected per plane, ordered around it and fanned, so a square
// face of an octahedron is two triangles rather than four overlapped ones.
// A shell wider than POLY_NEIGHBOURS_MAX is trimmed to its nearest points (a
// real coordination shell is 4 to 16 neighbours; anything wider is a bond
// detector blur that would cost the hull its quadratic worst case for no
// picture). A structure with more centres than POLY_MAX_CENTERS falls back to
// ball and stick entirely: the style must never stall a rebuild.
//
// The page carries no convex-hull library and must not grow one for this: the
// arithmetic below is a dozen lines of geometry, not a dependency.
var POLY_CENTERS=(function(){
    var s='Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr '+
        'Nb Mo Tc Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho '+
        'Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np '+
        'Pu Am Cm Bk Cf B Si Ge As Sb Te';
    var out={},list=s.split(' ');
    for(var i=0;i<list.length;i++)out[list[i]]=1;
    return out;
})();
var POLY_MAX_CENTERS=8000;
var POLY_NEIGHBOURS_MAX=24;
function isPolyCenter(el){return !!POLY_CENTERS[el]}

// A colour string from an atom record as three 0..1 channels. Anything the
// page does not recognise reads as neutral grey rather than nothing.
function polyRgb(c){
    if(typeof c==='string'&&c.charAt(0)==='#'&&c.length>=7){
        var r=parseInt(c.substr(1,2),16),g=parseInt(c.substr(3,2),16),b=parseInt(c.substr(5,2),16);
        if(isFinite(r)&&isFinite(g)&&isFinite(b))return [r/255,g/255,b/255];
    }
    return [0.72,0.72,0.72];
}
// Andrew's monotone chain over the projected points: the convex polygon of a
// planar face, counter-clockwise in the plane's own basis (which the caller
// orients so that counter-clockwise is outward).
function polyChain(q){
    var m=q.length,i,p;
    if(m<3)return [];
    var idx=[];
    for(i=0;i<m;i++)idx.push(i);
    idx.sort(function(a,b){
        if(q[a][0]!==q[b][0])return q[a][0]-q[b][0];
        return q[a][1]-q[b][1];
    });
    function cross(o,a,b){
        return (q[a][0]-q[o][0])*(q[b][1]-q[o][1])-(q[a][1]-q[o][1])*(q[b][0]-q[o][0]);
    }
    var lower=[];
    for(i=0;i<m;i++){
        p=idx[i];
        while(lower.length>=2&&cross(lower[lower.length-2],lower[lower.length-1],p)<=1e-12)lower.pop();
        lower.push(p);
    }
    var upper=[];
    for(i=m-1;i>=0;i--){
        p=idx[i];
        while(upper.length>=2&&cross(upper[upper.length-2],upper[upper.length-1],p)<=1e-12)upper.pop();
        upper.push(p);
    }
    lower.pop();upper.pop();
    var out=lower.concat(upper);
    return out.length>=3?out:[];
}
// The hull of one shell: faces as index triples and edges as index pairs, both
// into pts, or null when the shell is not a solid. Points are relative to the
// centre they surround, so the tolerances scale with the shell's own size.
function polyHull(pts){
    var n=pts.length,i,j,k,m;
    if(n<4)return null;
    var scale=0;
    for(i=0;i<n;i++){
        var d=Math.sqrt(pts[i][0]*pts[i][0]+pts[i][1]*pts[i][1]+pts[i][2]*pts[i][2]);
        if(d>scale)scale=d;
    }
    if(!(scale>1e-6))return null;
    var eps=1e-6*scale;
    var planes=[];
    for(i=0;i<n-2;i++){
        for(j=i+1;j<n-1;j++){
            for(k=j+1;k<n;k++){
                var ax=pts[i][0],ay=pts[i][1],az=pts[i][2];
                var ux=pts[j][0]-ax,uy=pts[j][1]-ay,uz=pts[j][2]-az;
                var vx=pts[k][0]-ax,vy=pts[k][1]-ay,vz=pts[k][2]-az;
                var nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
                var len=Math.sqrt(nx*nx+ny*ny+nz*nz);
                if(!(len>1e-8*scale*scale))continue;
                nx/=len;ny/=len;nz/=len;
                var pd=nx*ax+ny*ay+nz*az;
                var pos=0,neg=0;
                for(m=0;m<n;m++){
                    if(m===i||m===j||m===k)continue;
                    var sd=nx*pts[m][0]+ny*pts[m][1]+nz*pts[m][2]-pd;
                    if(sd>eps)pos++;
                    else if(sd<-eps)neg++;
                }
                if(pos&&neg)continue;
                if(!pos&&!neg)continue;
                if(pos){nx=-nx;ny=-ny;nz=-nz;pd=-pd}
                var pl=null;
                for(m=0;m<planes.length;m++){
                    var p0=planes[m];
                    if(p0.nx*nx+p0.ny*ny+p0.nz*nz>1-1e-6&&Math.abs(p0.d-pd)<eps){pl=p0;break}
                }
                if(!pl){pl={nx:nx,ny:ny,nz:nz,d:pd};planes.push(pl)}
            }
        }
    }
    if(planes.length<4)return null;
    var faces=[],edges=[],seen={};
    for(m=0;m<planes.length;m++){
        var pl2=planes[m];
        // The plane's own basis: u in the plane, v = n x u, so u x v = n.
        var bux,buy,buz;
        if(Math.abs(pl2.nx)<=Math.abs(pl2.ny)&&Math.abs(pl2.nx)<=Math.abs(pl2.nz)){bux=0;buy=-pl2.nz;buz=pl2.ny}
        else if(Math.abs(pl2.ny)<=Math.abs(pl2.nz)){bux=pl2.nz;buy=0;buz=-pl2.nx}
        else{bux=-pl2.ny;buy=pl2.nx;buz=0}
        var bl=Math.sqrt(bux*bux+buy*buy+buz*buz);
        bux/=bl;buy/=bl;buz/=bl;
        var bvx=pl2.ny*buz-pl2.nz*buy,bvy=pl2.nz*bux-pl2.nx*buz,bvz=pl2.nx*buy-pl2.ny*bux;
        var on=[],q=[];
        for(i=0;i<n;i++){
            var off=pl2.nx*pts[i][0]+pl2.ny*pts[i][1]+pl2.nz*pts[i][2]-pl2.d;
            if(Math.abs(off)>1e-5*scale)continue;
            on.push(i);
            q.push([pts[i][0]*bux+pts[i][1]*buy+pts[i][2]*buz,
                    pts[i][0]*bvx+pts[i][1]*bvy+pts[i][2]*bvz,i]);
        }
        var ring=polyChain(q);
        if(ring.length<3)continue;
        var poly=[];
        for(i=0;i<ring.length;i++)poly.push(on[ring[i]]);
        for(i=1;i+1<poly.length;i++)faces.push([poly[0],poly[i],poly[i+1]]);
        for(i=0;i<poly.length;i++){
            var e1=poly[i],e2=poly[(i+1)%poly.length];
            var key=e1<e2?e1+'_'+e2:e2+'_'+e1;
            if(seen[key])continue;
            seen[key]=1;
            edges.push([e1,e2]);
        }
    }
    if(faces.length<4)return null;
    // A nearly flat shell passes the plane tests as a paper-thin solid; its
    // volume is what tells it apart from a polyhedron worth drawing.
    var vol=0;
    for(i=0;i<faces.length;i++){
        var A=pts[faces[i][0]],B=pts[faces[i][1]],C=pts[faces[i][2]];
        vol+=A[0]*(B[1]*C[2]-B[2]*C[1])-A[1]*(B[0]*C[2]-B[2]*C[0])+A[2]*(B[0]*C[1]-B[1]*C[0]);
    }
    vol=Math.abs(vol)/6;
    if(vol<1e-3*scale*scale*scale)return null;
    return {faces:faces,edges:edges};
}
// Every centre's hull as one merged, indexed body, in group-local coordinates
// (the centre is at cx,cy,cz). A shell point is one vertex, shared by every
// face and every edge that meets there, and the triangles and segments are
// index lists into those vertices: a twelve-vertex polyhedron costs twelve
// vertices, not three per face corner. Each vertex also remembers the atom and
// the drawn lattice image it stands on, which is what lets a drag move the
// picture by rewriting the position buffer instead of rebuilding the hull.
// Pure: the caller draws it. shiftOf(bond) gives the lattice offset of the
// bond's second endpoint, or null when the bond is not periodic; the two
// images are then symmetric, which is what makes the hull of a centre at a
// cell face close across it.
function polyhedraGeometry(atoms,bonds,cx,cy,cz,shiftOf){
    if(!atoms||atoms.length===0||!bonds||bonds.length===0)return null;
    var n=atoms.length,nb=new Array(n),i,j,t;
    function add(ci,nel,ox,oy,oz,px,py,pz,ai,sx,sy,sz){
        var l=nb[ci];
        if(!l){l=[];nb[ci]=l}
        var dx=px-ox,dy=py-oy,dz=pz-oz;
        // The same image twice (two bond records, or a bond and its reverse)
        // is one neighbour; two different images are two.
        for(t=0;t<l.length;t++){
            if(Math.abs(l[t].x-dx)<1e-6&&Math.abs(l[t].y-dy)<1e-6&&Math.abs(l[t].z-dz)<1e-6)return;
        }
        l.push({el:nel,x:dx,y:dy,z:dz,atom:ai,sx:sx,sy:sy,sz:sz});
    }
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        if(!b)continue;
        var a1=atoms[b.atom1],a2=atoms[b.atom2];
        if(!a1||!a2)continue;
        if(isGhostAtom(a1)||isGhostAtom(a2))continue;
        var sh=shiftOf?shiftOf(b):null;
        var sx=sh?sh[0]:0,sy=sh?sh[1]:0,sz=sh?sh[2]:0;
        add(b.atom1,a2.element,a1.x,a1.y,a1.z,a2.x+sx,a2.y+sy,a2.z+sz,b.atom2,sx,sy,sz);
        add(b.atom2,a1.element,a2.x,a2.y,a2.z,a1.x-sx,a1.y-sy,a1.z-sz,b.atom1,-sx,-sy,-sz);
    }
    var eligible=0;
    for(i=0;i<n;i++){
        var a=atoms[i];
        if(!a||isGhostAtom(a)||!isPolyCenter(a.element))continue;
        var list=nb[i];
        if(!list||list.length<4)continue;
        for(j=0;j<list.length;j++){if(!isPolyCenter(list[j].el)){eligible++;break}}
    }
    if(eligible===0||eligible>POLY_MAX_CENTERS)return null;
    var positions=[],faceCols=[],edgeCols=[],faceIdx=[],edgeIdx=[],bindings=[];
    var centers=0,triangles=0,segments=0;
    for(i=0;i<n;i++){
        var c=atoms[i];
        if(!c||isGhostAtom(c)||!isPolyCenter(c.element))continue;
        var shell=nb[i];
        if(!shell||shell.length<4)continue;
        var other=false;
        for(j=0;j<shell.length;j++){if(!isPolyCenter(shell[j].el)){other=true;break}}
        if(!other)continue;
        if(shell.length>POLY_NEIGHBOURS_MAX){
            shell.sort(function(p,q){
                return (p.x*p.x+p.y*p.y+p.z*p.z)-(q.x*q.x+q.y*q.y+q.z*q.z);
            });
            shell.length=POLY_NEIGHBOURS_MAX;
        }
        var pts=[];
        for(j=0;j<shell.length;j++)pts.push([shell[j].x,shell[j].y,shell[j].z]);
        var hull=polyHull(pts);
        if(!hull)continue;
        centers++;
        var rgb=polyRgb(c.color);
        var er=rgb[0]*0.68,eg=rgb[1]*0.68,eb=rgb[2]*0.68;
        // One vertex per shell point, then index lists. A point the hull did
        // not use (a neighbour inside the convex hull) stays as an unreferenced
        // vertex: harmless, and it keeps the mapping trivial.
        var base=bindings.length;
        for(j=0;j<shell.length;j++){
            var sh2=shell[j];
            positions.push(c.x+sh2.x-cx,c.y+sh2.y-cy,c.z+sh2.z-cz);
            faceCols.push(rgb[0],rgb[1],rgb[2]);
            edgeCols.push(er,eg,eb);
            bindings.push([sh2.atom,sh2.sx,sh2.sy,sh2.sz]);
        }
        for(j=0;j<hull.faces.length;j++){
            var f=hull.faces[j];
            faceIdx.push(base+f[0],base+f[1],base+f[2]);
            triangles++;
        }
        for(j=0;j<hull.edges.length;j++){
            var e=hull.edges[j];
            edgeIdx.push(base+e[0],base+e[1]);
            segments++;
        }
    }
    if(centers===0)return null;
    return {positions:positions,faceColors:faceCols,edgeColors:edgeCols,
        faceIndex:faceIdx,edgeIndex:edgeIdx,bindings:bindings,
        centers:centers,triangles:triangles,segments:segments,vertices:bindings.length};
}
// The lattice offset a bond's second endpoint is drawn at (the same shift
// placeBondPieces applies), or null for a bond that stays inside the cell.
function polyShiftVector(b){
    if(!CRY||!b||!b.shift||!CRY.latticeVectors)return null;
    var s=b.shift;
    if(!s[0]&&!s[1]&&!s[2])return null;
    var lv=CRY.latticeVectors;
    return [s[0]*lv[0][0]+s[1]*lv[1][0]+s[2]*lv[2][0],
            s[0]*lv[0][1]+s[1]*lv[1][1]+s[2]*lv[2][1],
            s[0]*lv[0][2]+s[1]*lv[1][2]+s[2]*lv[2][2]];
}
// The coordination shells of a periodic structure, completed for the hulls.
//
// rebuildCrystal() draws one bond per base pair, so a centre on a face or a
// corner of the displayed block sees only the neighbours inside it: its hull
// would be left open, or with three neighbours dropped entirely. That drawing
// is right for the sticks — the image really is outside the block — but the
// solid this style exists to show should still be the centre's coordination
// polyhedron. So for every bonded base pair the code collects every lattice
// image of the partner within the pair cutoff, and binds each image to the
// nearest drawn copy of that site with the repeat difference as the shift.
// Vertices therefore land exactly where the infinite lattice puts them, and a
// drag carries them with a real atom, the way a split bond's stub follows its
// drawn end.
//
// Shells are rebuilt from CRY.baseBonds on each call, so a deleted bond
// disappears at once; a cell with no polyhedron centre skips the O(base^2)
// pair scan outright. Only a periodic component is completed: a finite
// molecule is drawn whole, so it already has every bond its hull needs.
function crystalPolyBonds(atoms){
    if(!CRY||!CRY_INV||!CRY.baseAtoms||!CRY.baseBonds||!CRY.latticeVectors)return null;
    var base=CRY.baseAtoms,baseBonds=CRY.baseBonds,lv=CRY.latticeVectors,inv=CRY_INV;
    var nb=base.length,i,j;
    if(!nb||!atoms||!atoms.length)return null;
    var anyCenter=false;
    for(i=0;i<nb;i++){
        if(base[i]&&isPolyCenter(base[i].element)){anyCenter=true;break}
    }
    if(!anyCenter)return null;

    // Every drawn copy of every base site, by site: which display atom stands
    // there, at which repeat of the site. A site's copies can be a whole cell
    // box (a periodic component) or a scatter (each copy of a molecule), so
    // the exact repeat is looked up in a map first and the nearest of a short
    // list is only a fallback.
    var images=[],imageMap=[],spanLo=[],spanHi=[];
    for(i=0;i<atoms.length;i++){
        var a=atoms[i];
        if(!a||isGhostAtom(a)||a.baseIdx==null||!a.imgShift)continue;
        var bi=a.baseIdx;
        if(bi<0||bi>=nb)continue;
        if(!images[bi]){
            images[bi]=[];imageMap[bi]={};
            spanLo[bi]=[1e9,1e9,1e9];spanHi[bi]=[-1e9,-1e9,-1e9];
        }
        var s0=a.imgShift,key0=s0[0]+'_'+s0[1]+'_'+s0[2];
        if(imageMap[bi][key0]!==undefined)continue;
        imageMap[bi][key0]=i;images[bi].push(i);
        for(var k=0;k<3;k++){
            if(s0[k]<spanLo[bi][k])spanLo[bi][k]=s0[k];
            if(s0[k]>spanHi[bi][k])spanHi[bi][k]=s0[k];
        }
    }
    // The drawn copy of a site closest to a wanted repeat, as
    // [repeat x, repeat y, repeat z, display index].
    function pickImage(site,want){
        var map=imageMap[site],list=images[site];
        if(!map||!list||!list.length)return null;
        var lo=spanLo[site],hi=spanHi[site],t;
        var q0=want[0]<lo[0]?lo[0]:(want[0]>hi[0]?hi[0]:want[0]);
        var q1=want[1]<lo[1]?lo[1]:(want[1]>hi[1]?hi[1]:want[1]);
        var q2=want[2]<lo[2]?lo[2]:(want[2]>hi[2]?hi[2]:want[2]);
        var at=map[q0+'_'+q1+'_'+q2];
        if(at!==undefined)return [q0,q1,q2,at];
        var best=-1,bd=Infinity;
        var lim=list.length<64?list.length:64;
        for(t=0;t<lim;t++){
            var s=atoms[list[t]].imgShift;
            var dx=(want[0]-s[0])*lv[0][0]+(want[1]-s[1])*lv[1][0]+(want[2]-s[2])*lv[2][0];
            var dy=(want[0]-s[0])*lv[0][1]+(want[1]-s[1])*lv[1][1]+(want[2]-s[2])*lv[2][1];
            var dz=(want[0]-s[0])*lv[0][2]+(want[1]-s[1])*lv[1][2]+(want[2]-s[2])*lv[2][2];
            var d2=dx*dx+dy*dy+dz*dz;
            if(d2<bd){bd=d2;best=list[t]}
        }
        if(best<0)return null;
        var sb=atoms[best].imgShift;
        return [sb[0],sb[1],sb[2],best];
    }

    var widths=crystalPlaneWidths(lv);
    if(!(widths[0]>1e-9)||!(widths[1]>1e-9)||!(widths[2]>1e-9))return null;
    var w0=widths[0]*widths[0],w1=widths[1]*widths[1],w2=widths[2]*widths[2];
    var comp=crystalComponents(baseBonds,nb);
    var groupOf=comp.groupOf,groups=comp.groups;
    var shells=new Array(nb);
    function addShell(site,partner,s){
        if(!shells[site])shells[site]=[];
        shells[site].push({partner:partner,s:s});
    }
    var seenPair={};
    for(i=0;i<baseBonds.length;i++){
        var bb=baseBonds[i];
        if(!bb)continue;
        var b1=bb.atom1,b2=bb.atom2;
        if(b1==null||b2==null||b1<0||b2<0||b1>=nb||b2>=nb||b1===b2)continue;
        var pk=b1<b2?b1+'_'+b2:b2+'_'+b1;
        if(seenPair[pk])continue;
        seenPair[pk]=true;
        var gid=groupOf[b1];
        if(gid==null||gid<0||!groups[gid].periodic)continue;
        if(!images[b1]||!images[b2])continue;
        var el1=base[b1].element,el2=base[b2].element;
        var c1=isPolyCenter(el1),c2=isPolyCenter(el2);
        if(!c1&&!c2)continue;
        var cut=BP.pairCutoff(el1,el2);
        if(!(cut>0))continue;
        var ax=base[b2].x-base[b1].x,ay=base[b2].y-base[b1].y,az=base[b2].z-base[b1].z;
        var fx=inv[0][0]*ax+inv[0][1]*ay+inv[0][2]*az;
        var fy=inv[1][0]*ax+inv[1][1]*ay+inv[1][2]*az;
        var fz=inv[2][0]*ax+inv[2][1]*ay+inv[2][2]*az;
        // The same nearest-image lower bound the bond detector uses: the
        // displacement measured across the three plane normals is a distance
        // no image of the pair can beat.
        var qx=fx-Math.round(fx),qy=fy-Math.round(fy),qz=fz-Math.round(fz);
        if(qx*qx*w0+qy*qy*w1+qz*qz*w2>cut*cut)continue;
        var rx=Math.round(fx),ry=Math.round(fy),rz=Math.round(fz);
        var nx=Math.ceil(cut/widths[0]+1e-9),ny=Math.ceil(cut/widths[1]+1e-9),nz=Math.ceil(cut/widths[2]+1e-9);
        if(nx>CRYSTAL_MAX_IMAGE_REPEAT)nx=CRYSTAL_MAX_IMAGE_REPEAT;
        if(ny>CRYSTAL_MAX_IMAGE_REPEAT)ny=CRYSTAL_MAX_IMAGE_REPEAT;
        if(nz>CRYSTAL_MAX_IMAGE_REPEAT)nz=CRYSTAL_MAX_IMAGE_REPEAT;
        for(var sx=rx-nx;sx<=rx+nx;sx++){
            for(var sy=ry-ny;sy<=ry+ny;sy++){
                for(var sz=rz-nz;sz<=rz+nz;sz++){
                    var ox=ax+sx*lv[0][0]+sy*lv[1][0]+sz*lv[2][0];
                    var oy=ay+sx*lv[0][1]+sy*lv[1][1]+sz*lv[2][1];
                    var oz=az+sx*lv[0][2]+sy*lv[1][2]+sz*lv[2][2];
                    if(ox*ox+oy*oy+oz*oz>cut*cut)continue;
                    // The bond points base2 at repeat s relative to base1;
                    // read the other way, base1 sits at -s from base2.
                    if(c1)addShell(b1,b2,[sx,sy,sz]);
                    if(c2)addShell(b2,b1,[-sx,-sy,-sz]);
                }
            }
        }
    }

    // Every image of a bonded partner is a shell point of each drawn centre
    // that is a lattice translation of its site: the centre's own repeat
    // carries the shell with it.
    var extra=[];
    for(i=0;i<atoms.length;i++){
        var c=atoms[i];
        if(!c||isGhostAtom(c)||c.baseIdx==null||!c.imgShift)continue;
        if(!isPolyCenter(c.element))continue;
        var list=shells[c.baseIdx];
        if(!list)continue;
        var ci=c.imgShift;
        for(j=0;j<list.length;j++){
            var e=list[j];
            var want=[ci[0]+e.s[0],ci[1]+e.s[1],ci[2]+e.s[2]];
            var pick=pickImage(e.partner,want);
            if(!pick)continue;
            extra.push({atom1:i,atom2:pick[3],order:1,
                shift:[want[0]-pick[0],want[1]-pick[1],want[2]-pick[2]]});
        }
    }
    return extra.length?extra:null;
}
// Turn a merged indexed body into the two objects the style draws: one mesh
// for the faces and one line set for the edges. They share the vertex buffer
// (the edges add only their own darker colour list and their own index list),
// so the deduplicated positions exist once in memory and are uploaded once a
// drag. Sharing a vertex between faces is what makes the geometry indexed, and
// flat shading is what keeps each facet a plane rather than a smooth blot: the
// material derives its face normals in the shader, so a deformed hull is cor-
// rectly lit while a drag moves its vertices, with no normal pass per frame.
// A body with no segments (impossible in practice, but the caller should not
// have to know) yields a null second mesh.
function polyMeshPair(body,parent){
    var fg=new THREE.BufferGeometry();
    var posAttr=new THREE.Float32BufferAttribute(body.positions,3);
    fg.setAttribute('position',posAttr);
    fg.setAttribute('color',new THREE.Float32BufferAttribute(body.faceColors,3));
    fg.setIndex(body.faceIndex);
    var fm=new THREE.MeshPhongMaterial({vertexColors:true,flatShading:true,transparent:true,opacity:0.36,
        side:THREE.DoubleSide,depthWrite:false,shininess:30,specular:0x222222});
    var faces=new THREE.Mesh(fg,fm);
    // The hulls sit over the atoms and sticks on purpose; drawing them after
    // the opaque geometry blends the faces over the picture, not into it.
    faces.renderOrder=2;
    // The hull is never rebuilt during a drag, so the bounding sphere computed
    // when it was built would be stale by the end of one; two objects are
    // cheaper to draw than to cull wrongly.
    faces.frustumCulled=false;
    parent.add(faces);
    var edges=null;
    if(body.segments>0){
        var eg=new THREE.BufferGeometry();
        eg.setAttribute('position',posAttr);
        eg.setAttribute('color',new THREE.Float32BufferAttribute(body.edgeColors,3));
        eg.setIndex(body.edgeIndex);
        var em=new THREE.LineBasicMaterial({vertexColors:true,transparent:true,opacity:0.9});
        edges=new THREE.LineSegments(eg,em);
        edges.renderOrder=3;
        edges.frustumCulled=false;
        parent.add(edges);
    }
    return [faces,edges];
}
function disposePolyMesh(m){
    if(!m)return;
    if(m.parent)m.parent.remove(m);
    disposeMesh(m);
}
var polyFacesMesh=null;
var polyEdgesMesh=null;
// The built body, kept for its vertex->atom bindings. The meshes on screen hold
// the same positions in their shared attribute; syncing reads the atoms and
// writes that attribute, and the body is what says which atom each vertex
// belongs to.
var polyBody=null;
// A drag rewrites coordinates without rebuilding the scene. Every per-frame
// path that moves atoms sets this: the hulls could follow their atoms exactly
// (polySyncPositions), but the hull *shapes* are now stale, and the commit that
// closes the drag rebuilds them. A centre-only shift does not set it — the
// relative geometry did not change, so moving the vertices is exact.
var polyShapeDirty=false;
function buildPolyhedra(){
    disposePolyMesh(polyFacesMesh);
    disposePolyMesh(polyEdgesMesh);
    polyFacesMesh=null;polyEdgesMesh=null;polyBody=null;
    polyShapeDirty=false;
    if(displayStyle!=='polyhedra')return;
    // A periodic structure's drawn bonds stop at the block faces. The hulls
    // complete them across a face (crystalPolyBonds), so a centre on the edge
    // of the block still shows its polyhedron; a molecule needs no help.
    var bonds=MD.bonds;
    if(CRY){
        var extra=crystalPolyBonds(MD.atoms);
        if(extra)bonds=(MD.bonds||[]).concat(extra);
    }
    var body=polyhedraGeometry(MD.atoms,bonds,CX,CY,CZ,polyShiftVector);
    if(!body)return;
    var pair=polyMeshPair(body,moleculeGroup);
    polyFacesMesh=pair[0];polyEdgesMesh=pair[1];polyBody=body;
}
// Move every hull vertex to the atom (and lattice image) it was built from.
// The topology is untouched: each vertex stays on the shell point it was
// computed for, so the hull deforms with the atoms the way the sticks do, and
// the closing rebuild is the only step that recomputes which points are a hull.
function polySyncPositions(){
    if(!polyBody||!polyFacesMesh)return;
    var attr=polyFacesMesh.geometry.attributes.position;
    var arr=attr.array;
    var bind=polyBody.bindings;
    for(var i=0;i<bind.length;i++){
        var a=MD.atoms[bind[i][0]];
        if(!a)continue;
        var o=i*3;
        arr[o]=a.x+bind[i][1]-CX;
        arr[o+1]=a.y+bind[i][2]-CY;
        arr[o+2]=a.z+bind[i][3]-CZ;
    }
    // The edge geometry shares this very attribute object, so one upload
    // carries both.
    attr.needsUpdate=true;
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

// ===== Periodic bonds: one neighbour search, then the shared detector =====
//
// A crystal file states one cell, but the chemistry happens in the lattice that
// cell is repeated in: a carbonate's carbon is bonded to an oxygen the file has
// wrapped into the next cell, so the bonds cannot be read off the centre cell
// alone. And they have to be the bonds the molecular path would find, or the
// same structure would be drawn with different orders depending on whether it
// arrived as a CIF or as an XYZ.
//
// So no chemistry is decided here. The cell's atoms are laid out the way the
// renderer draws them — a finite molecule whole, at the unwrap offsets
// rebuildCrystal() places it by, and a periodic network as the centre cell plus
// the neighbour images its split bonds point at — and the page's own detector
// (detectBondsFromAtoms — the perception module the host and every molecular
// path already run) perceives that set. Each bond it returns is written back as
// the base pair plus the integer image it points at, which is the form
// rebuildCrystal() draws and every editor writes back into the cell. Orders,
// the aromatic 1.5, the valence relaxation: all of it is whatever the shared
// module says, not a second table kept here.
//
// The neighbour search is enumerated, never rounded. Rounding the fractional
// difference gives the shortest image only for an orthogonal cell; in a skewed
// one it can leave a bond a whole cell long (drawn straight through the
// structure) or choose an image that is not the nearest at all, which loses the
// bond. The enumeration is bounded rather than guessed: a repeat s_i moves a
// displacement only across the planes the other two lattice vectors span, at a
// fixed spacing w_i, so an image is at least |s_i - f_i| * w_i away — every
// image a pair can reach lies inside s_i in [ceil(f_i - n_i), floor(f_i + n_i)]
// with n_i = cut / w_i. For an ordinary cell n_i is 1 and that is the classic
// block of 27 images; a flat or skewed cell widens it to 5x5x5, which is
// exactly the case rounding got wrong.
// The same search, over any cell — the loaded one, or the right-hand cell of a
// crystal comparison, which is not the structure on screen and must not be
// written into CRY. Returns the base bonds (pair + image shift) and touches
// nothing else.
function detectCrystalBaseBondsFor(cry,inv){
    if(!cry||!cry.latticeVectors||!cry.baseAtoms||!inv)return[];
    var baseAtoms=cry.baseAtoms;
    var lv=cry.latticeVectors;
    var n=baseAtoms.length;
    var widths=crystalPlaneWidths(lv);
    var round=Math.round;
    var i,j;

    // The widest cutoff each element can have against anything in this cell.
    // The pair loop runs n^2/2 times, so it reads one number per atom rather
    // than asking the table — which normalises two element strings — per pair;
    // the exact pair cutoff is asked only for the few pairs that survive the
    // cheap rejection below.
    var elements=[],seenEl={};
    for(i=0;i<n;i++){
        var elI=baseAtoms[i].element;
        if(!seenEl[elI]){seenEl[elI]=true;elements.push(elI)}
    }
    var widestOf={};
    for(i=0;i<elements.length;i++){
        var widest=0;
        for(j=0;j<elements.length;j++){
            var wc=BP.pairCutoff(elements[i],elements[j]);
            if(wc>widest)widest=wc;
        }
        widestOf[elements[i]]=widest;
    }
    var capOf=new Float64Array(n);
    for(i=0;i<n;i++)capOf[i]=widestOf[baseAtoms[i].element]||0;

    // The bonds the search found, then the atoms they join — laid out for the
    // detector exactly as rebuildCrystal() draws them. Plain records holding
    // only the fields the detector reads: the perception module runs the same
    // code the host does, and nothing else about an atom is its business.
    var probe=[];
    var probeAt=[];   // per probe atom: the base site and the cell offset of this copy
    var pairs=[];     // the bonded base pairs, in (atom1 asc, atom2 asc) order
    var pairAt={};

    var w0=widths[0]*widths[0],w1=widths[1]*widths[1],w2=widths[2]*widths[2];
    for(i=0;i<n;i++){
        var ai=baseAtoms[i];
        var ax=ai.x,ay=ai.y,az=ai.z;
        var capI=capOf[i];
        for(j=i+1;j<n;j++){
            var bj=baseAtoms[j];
            var dx=bj.x-ax,dy=bj.y-ay,dz=bj.z-az;
            var fx=inv[0][0]*dx+inv[0][1]*dy+inv[0][2]*dz;
            var fy=inv[1][0]*dx+inv[1][1]*dy+inv[1][2]*dz;
            var fz=inv[2][0]*dx+inv[2][1]*dy+inv[2][2]*dz;
            // The pairs that can bond are a small share of all of them, and the
            // cheapest test is arithmetic: no image of either atom can be closer
            // than the displacement measured across the three plane normals.
            var qx=fx-round(fx),qy=fy-round(fy),qz=fz-round(fz);
            var low=qx*qx*w0+qy*qy*w1+qz*qz*w2;
            var cap=capI<capOf[j]?capI:capOf[j];
            if(low>cap*cap)continue;
            var cut=BP.pairCutoff(ai.element,bj.element);
            if(!(cut>0)||low>cut*cut)continue;
            var near=crystalShortestFractional(fx,fy,fz,cut,lv,widths);
            if(near.dist>cut)continue;
            // near.shift is the repeat subtracted from bj's fractional
            // position, so the bond points at the opposite cell.
            var shift=[-near.shift[0],-near.shift[1],-near.shift[2]];
            pairAt[i+'_'+j+'_'+shift[0]+'_'+shift[1]+'_'+shift[2]]=pairs.length;
            pairs.push({i:i,j:j,shift:shift,order:0});
        }
    }

    // The periodic graph's components, walked once with an image offset per
    // atom — the same unwrap rebuildCrystal() uses, so the detector and the
    // renderer cannot disagree about which atoms form a molecule. A finite
    // component is a molecule, and it is laid out once per site at its
    // unwrapped position: that is exactly the molecule the molecular path
    // perceives, so its orders are the molecular path's. A periodic component
    // has no single consistent layout (a cycle comes back a cell over), so it
    // is laid out as the centre cell plus the images the search found — the
    // same cut the renderer draws.
    var comp=crystalComponents(pairs.map(function(p){return{atom1:p.i,atom2:p.j,shift:p.shift}}),n);
    var groupOf=comp.groupOf,groups=comp.groups;
    function addProbe(baseIdx,off){
        var b=baseAtoms[baseIdx];
        probe.push({
            element:b.element,
            x:b.x+off[0]*lv[0][0]+off[1]*lv[1][0]+off[2]*lv[2][0],
            y:b.y+off[0]*lv[0][1]+off[1]*lv[1][1]+off[2]*lv[2][1],
            z:b.z+off[0]*lv[0][2]+off[1]*lv[1][2]+off[2]*lv[2][2]
        });
        probeAt.push({base:baseIdx,shift:off});
    }
    groups.forEach(function(g){
        var k;
        if(g.periodic){
            for(k=0;k<g.members.length;k++)addProbe(g.members[k],ZERO_SHIFT);
        }else{
            for(k=0;k<g.members.length;k++)addProbe(g.members[k],g.offset[g.members[k]]);
        }
    });
    var imageAt={};
    pairs.forEach(function(p){
        var g=groupOf[p.i];
        if(g===undefined||!groups[g].periodic)return;
        if(!(p.shift[0]||p.shift[1]||p.shift[2]))return;
        var key=p.j+'_'+p.shift[0]+'_'+p.shift[1]+'_'+p.shift[2];
        if(imageAt[key])return;
        imageAt[key]=true;
        addProbe(p.j,p.shift);
    });

    // The detector runs on that set exactly as it runs on a molecule. A cell
    // states no total charge, so it gets the neutral assumption a bare
    // coordinate set gets.
    var found=detectBondsFromAtoms(probe);

    // Read every perceived bond back in the cell's own terms: which two base
    // sites it joins, and how many cells apart the copies are. The image search
    // already decided that number for each pair, so a bond whose two copies are
    // further apart than the pair's own shortest image is not this cell's bond
    // and is left out — the pair list holds one bond per base pair, as the
    // renderer and every editor expect.
    for(i=0;i<found.length;i++){
        var b=found[i];
        var at1=probeAt[b.atom1],at2=probeAt[b.atom2];
        var base1=at1.base,base2=at2.base;
        if(base1===base2)continue;
        var sh=[at2.shift[0]-at1.shift[0],at2.shift[1]-at1.shift[1],at2.shift[2]-at1.shift[2]];
        if(base1>base2){
            var t=base1;base1=base2;base2=t;
            sh=[-sh[0],-sh[1],-sh[2]];
        }
        var at=pairAt[base1+'_'+base2+'_'+sh[0]+'_'+sh[1]+'_'+sh[2]];
        // Two images exactly equidistant is a tie the search broke towards the
        // larger repeat, and the walk may have taken the same tie from the
        // other end. The pair with the opposite shift is the bond this pair
        // holds.
        if(at===undefined)at=pairAt[base1+'_'+base2+'_'+(-sh[0])+'_'+(-sh[1])+'_'+(-sh[2])];
        if(at===undefined)continue;
        pairs[at].order=b.order;
    }

    var bonds=[];
    for(i=0;i<pairs.length;i++){
        // A pair the search reached but the detector would not accept — two
        // sites closer than any real bond, a length no table supports — stays
        // out: the order is what the perception says, and here it says nothing.
        if(!(pairs[i].order>0))continue;
        bonds.push({atom1:pairs[i].i,atom2:pairs[i].j,order:pairs[i].order,shift:pairs[i].shift});
    }
    return bonds;
}
function detectCrystalBaseBonds(){
    if(!CRY||!CRY_INV)return;
    CRY.baseBonds=detectCrystalBaseBondsFor(CRY,CRY_INV);
    mergeCrystalFileBonds();
    CRY.bondsDetected=true;
}

// The file's own bond table (_geom_bond), merged over the detected list once.
//
// The table is the depositor's statement, so it wins where it speaks: a bond it
// lists keeps the author's order and valence instead of the perceived one, and
// when perception finds nothing at all (an unusual element, a long bond outside
// every radius cutoff) the table is the only connectivity there is and is used
// as is.
//
// Where it is silent, detection stands. A table that lists only some of the
// bonds — CCDC 142619 lists its heavy-atom bonds and no C-H at all — must not
// orphan every hydrogen in the cell, so the detected bonds the file did not
// mention stay. A row the table lists but perception did not find is counted
// and reported, not drawn: one bad row must not put a bond through a molecule
// that the geometry says is not there, and the panel line is where the file's
// disagreement is visible.
function mergeCrystalFileBonds(){
    if(!CRY)return;
    var file=CRY.fileBonds||[];
    CRY.bondSource=null;
    if(!file.length)return;
    var bonds=CRY.baseBonds||[];
    var i,f;
    if(bonds.length===0){
        var only=[];
        for(f=0;f<file.length;f++){
            var fbOnly=file[f];
            if(!fbOnly||fbOnly.atom1===fbOnly.atom2)continue;
            only.push({atom1:fbOnly.atom1,atom2:fbOnly.atom2,order:fbOnly.order||1,
                shift:fbOnly.shift||ZERO_SHIFT,valence:fbOnly.valence});
        }
        CRY.baseBonds=only;
        CRY.bondSource={file:file.length,detected:0,matched:only.length,fileOnly:0};
        return;
    }
    // One entry per undirected pair and shift, in the (atom1 < atom2) direction
    // the cell stores, so a table row written from either end finds the edge.
    var byPair={};
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        var lo=b.atom1<b.atom2?b.atom1:b.atom2;
        var hi=b.atom1<b.atom2?b.atom2:b.atom1;
        var s=b.shift||ZERO_SHIFT;
        if(b.atom1>b.atom2)s=[-s[0],-s[1],-s[2]];
        byPair[lo+'_'+hi+'_'+s[0]+'_'+s[1]+'_'+s[2]]=b;
    }
    var matched=0;
    for(f=0;f<file.length;f++){
        var fb=file[f];
        if(!fb||fb.atom1===fb.atom2)continue;
        var lo2=fb.atom1<fb.atom2?fb.atom1:fb.atom2;
        var hi2=fb.atom1<fb.atom2?fb.atom2:fb.atom1;
        var fs=fb.shift||ZERO_SHIFT;
        if(fb.atom1>fb.atom2)fs=[-fs[0],-fs[1],-fs[2]];
        var hit=byPair[lo2+'_'+hi2+'_'+fs[0]+'_'+fs[1]+'_'+fs[2]];
        // A tie between two equidistant images is stored with one shift; the
        // file may state the other. It is the same edge either way.
        if(!hit)hit=byPair[lo2+'_'+hi2+'_'+(-fs[0])+'_'+(-fs[1])+'_'+(-fs[2])];
        if(!hit)continue;
        matched++;
        if(fb.order>0)hit.order=fb.order;
        if(fb.valence>0)hit.valence=fb.valence;
    }
    CRY.bondSource={file:file.length,detected:bonds.length,matched:matched,fileOnly:file.length-matched};
}

// The connected components of a base bond list, each walked once with an
// integer image offset per atom — the unwrap. Bond (a1,a2,shift s) means the
// image of a2 in cell c+s bonds with the image of a1 in cell c, so walking the
// edge a1->a2 adds s to the image offset and a2->a1 subtracts it. A cycle that
// comes back shifted marks the whole component periodic: that is what an
// infinite network is, and what a molecule never is. Offsets are exact
// integers, so the consistency check cannot drift.
//
// Two readers share this walk. rebuildCrystal() places a finite component by
// these offsets (a molecule is drawn whole), and detectCrystalBaseBonds() lays
// a finite component out at them so the detector perceives the molecule the
// molecular path would.
function crystalComponents(bondList,n){
    var adj={};var groupOf={};var groups=[];
    bondList.forEach(function(bb){
        if(!adj[bb.atom1])adj[bb.atom1]=[];
        if(!adj[bb.atom2])adj[bb.atom2]=[];
        adj[bb.atom1].push({nb:bb.atom2,s:bb.shift,sign:1});
        adj[bb.atom2].push({nb:bb.atom1,s:bb.shift,sign:-1});
    });
    for(var bi=0;bi<n;bi++){
        if(groupOf[bi]!==undefined)continue;
        if(!adj[bi]){groupOf[bi]=-1;continue}
        var gid=groups.length;var queue=[bi];groupOf[bi]=gid;
        var offset={};offset[bi]=[0,0,0];
        var periodic=false;
        for(var qi=0;qi<queue.length;qi++){
            var cur=queue[qi];
            var oc=offset[cur];
            var nbrs=adj[cur]||[];
            for(var ni=0;ni<nbrs.length;ni++){
                var e=nbrs[ni];
                var s=e.s||ZERO_SHIFT;
                var no=[oc[0]+e.sign*s[0],oc[1]+e.sign*s[1],oc[2]+e.sign*s[2]];
                if(offset[e.nb]!==undefined){
                    var po=offset[e.nb];
                    if(no[0]!==po[0]||no[1]!==po[1]||no[2]!==po[2])periodic=true;
                    continue;
                }
                offset[e.nb]=no;
                groupOf[e.nb]=gid;
                queue.push(e.nb);
            }
        }
        groups.push({members:queue.slice(),offset:offset,periodic:periodic});
    }
    return {groups:groups,groupOf:groupOf};
}

// The spacing between the lattice planes a repeat along each axis moves a
// displacement across: the cell volume over the area the other two vectors
// span. This is the lattice's own measure of "one cell along this axis", and
// what bounds how many repeats a search can possibly need.
function crystalPlaneWidths(lv){
    var a=lv[0],b=lv[1],c=lv[2];
    var bcx=b[1]*c[2]-b[2]*c[1],bcy=b[2]*c[0]-b[0]*c[2],bcz=b[0]*c[1]-b[1]*c[0];
    var cax=c[1]*a[2]-c[2]*a[1],cay=c[2]*a[0]-c[0]*a[2],caz=c[0]*a[1]-c[1]*a[0];
    var abx=a[1]*b[2]-a[2]*b[1],aby=a[2]*b[0]-a[0]*b[2],abz=a[0]*b[1]-a[1]*b[0];
    var det=Math.abs(a[0]*bcx+a[1]*bcy+a[2]*bcz);
    return [det/Math.sqrt(bcx*bcx+bcy*bcy+bcz*bcz),
            det/Math.sqrt(cax*cax+cay*cay+caz*caz),
            det/Math.sqrt(abx*abx+aby*aby+abz*abz)];
}

// The shortest periodic image of a fractional displacement: which repeats to
// subtract, what is left inside one cell, and how long that is.
//
// cut is what bounds the search — the distance the caller cares about — so
// every image that could be within it is enumerated and none beyond it is.
// A displacement exactly midway between two cells is a tie, and the walk runs
// from the largest repeat down so a tie keeps the answer rounding gave.
function crystalShortestFractional(fx,fy,fz,cut,lv,widths){
    var n=[1,1,1],k,q;
    for(k=0;k<3;k++){
        q=widths[k]>1e-9?cut/widths[k]:0;
        n[k]=q>1?Math.ceil(q+1e-9):1;
        if(!(n[k]>=1)||n[k]>CRYSTAL_MAX_IMAGE_REPEAT)n[k]=CRYSTAL_MAX_IMAGE_REPEAT;
    }
    var loX=Math.ceil(fx-n[0]-1e-9),hiX=Math.floor(fx+n[0]+1e-9);
    var loY=Math.ceil(fy-n[1]-1e-9),hiY=Math.floor(fy+n[1]+1e-9);
    var loZ=Math.ceil(fz-n[2]-1e-9),hiZ=Math.floor(fz+n[2]+1e-9);
    if(hiX<loX){loX=hiX=Math.round(fx)}
    if(hiY<loY){loY=hiY=Math.round(fy)}
    if(hiZ<loZ){loZ=hiZ=Math.round(fz)}
    var best=null,bestD2=Infinity;
    for(var sx=hiX;sx>=loX;sx--){
        var gx=fx-sx;
        for(var sy=hiY;sy>=loY;sy--){
            var gy=fy-sy;
            for(var sz=hiZ;sz>=loZ;sz--){
                var gz=fz-sz;
                var ax=gx*lv[0][0]+gy*lv[1][0]+gz*lv[2][0];
                var ay=gx*lv[0][1]+gy*lv[1][1]+gz*lv[2][1];
                var az=gx*lv[0][2]+gy*lv[1][2]+gz*lv[2][2];
                var d2=ax*ax+ay*ay+az*az;
                if(d2<bestD2-1e-12){
                    bestD2=d2;
                    best=[sx,sy,sz,gx,gy,gz];
                }
            }
        }
    }
    return {shift:[best[0],best[1],best[2]],frac:[best[3],best[4],best[5]],dist:Math.sqrt(bestD2)};
}

// ---------------------------------------------------------------------------
// Bond-valence sums (H7)
// ---------------------------------------------------------------------------
// A period-correct reading of how well each site's coordination matches its
// formal valence. The arithmetic lives in the BVS module — the same text the
// extension host runs — and this side gathers the cell as the picture shows
// it, keeps the result until the cell changes, and presents it. Nothing here
// changes the structure: contacts are counted, not turned into bonds, and no
// coordinate is written.
var bvsCache=null;
// Above this many atoms the report is no longer cheap enough to compute on a
// mouse-move; the tooltip stays silent and the button still offers the full
// table, which is where the reader asked for it.
var BVS_HOVER_LIMIT=800;
var BVS_SETTLED=0.25;

function bvsInvalidate(){bvsCache=null}

// The first image in the base cell of every asymmetric-unit site. A general
// position appears once per symmetry operation; every image has the same sum,
// so the table lists the site once and quotes the first image.
function bvsAsymFirstImage(){
    var map={};
    var atoms=(CRY&&CRY.baseAtoms)||[];
    for(var i=0;i<atoms.length;i++){
        var ai=atoms[i].asymIdx;
        if(ai!=null&&map[ai]==null)map[ai]=i;
    }
    return map;
}

// The report is computed over exactly the cell atoms the picture shows: with a
// disorder group selected, the hidden alternative is not a contact and takes
// no part. Each result is addressed back through reportOfBase, so an atom the
// filter removed simply has no entry.
function bvsEnsureReport(){
    if(!CRY)return null;
    if(bvsCache)return bvsCache;
    var src=CRY.baseAtoms||[],atoms=[],baseOf=[];
    for(var i=0;i<src.length;i++){
        if(neighBaseVisible(src[i])){atoms.push(src[i]);baseOf.push(i)}
    }
    var report=BVS.analyse(atoms,CRY.latticeVectors);
    var reportOfBase={};
    for(var j=0;j<baseOf.length;j++)reportOfBase[baseOf[j]]=j;
    bvsCache={report:report,baseOf:baseOf,reportOfBase:reportOfBase};
    return bvsCache;
}

function bvsChargeText(r){
    if(r.kind==='cation')return r.element+(r.state===1?'+':r.state+'+');
    return r.element+(r.valence===1?'-':r.valence+'-');
}

// One site's computed sum, or null when the site is not a bond-valence centre,
// has no usable parameter, or has nothing inside the cutoffs.
function bvsResultOfBase(baseIdx){
    var c=bvsEnsureReport();
    if(!c||!c.report.ok||baseIdx==null)return null;
    var vi=c.reportOfBase[baseIdx];
    if(vi==null)return null;
    var r=c.report.perAtom[vi];
    return r&&r.status==='ok'?r:null;
}

// The hover line adds itself only when a sum exists: an organic structure has
// no centres in the parameter table, and its tooltips stay exactly as they
// were. The silence above BVS_HOVER_LIMIT is the same promise — a hover must
// never sit and compute a framework's worth of contacts.
function bvsTooltipSuffix(a){
    if(!CRY||diffMode||!a||a.baseIdx==null)return '';
    if(!bvsCache&&(CRY.baseAtoms||[]).length>BVS_HOVER_LIMIT)return '';
    var r=bvsResultOfBase(a.baseIdx);
    if(!r)return '';
    return ' · BVS Σs='+r.sum.toFixed(2)+' |Δ|='+r.delta.toFixed(2)+' ('+bvsChargeText(r)+')';
}

function showBvsReport(){
    if(!CRY)return;
    // The sums describe one cell; a comparison shows two structures, and the
    // left one alone would be a report about half the picture.
    if(diffMode){
        showModal('<h3>Bond Valence (BVS)</h3>'+
            '<div class="current-val">Bond-valence sums read one unit cell, and Diff shows two structures. Close the comparison to inspect a cell.</div>'+
            '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-close">Close</button></div>',null,true);
        document.getElementById('m-close').addEventListener('click',hideModal);
        return;
    }
    var wrapped=bvsEnsureReport();
    var rep=wrapped&&wrapped.report;
    if(!rep||!rep.ok){
        showModal('<h3>Bond Valence (BVS)</h3>'+
            '<div class="current-val">'+escapeHtml(rep.reason)+'</div>'+
            '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-close">Close</button></div>',null,true);
        document.getElementById('m-close').addEventListener('click',hideModal);
        return;
    }
    var first=bvsAsymFirstImage();
    var entries=[];
    if(CRY.asymAtoms&&CRY.asymAtoms.length){
        for(var s=0;s<CRY.asymAtoms.length;s++){
            var sa=CRY.asymAtoms[s];
            var b=first[sa.index];
            var vi=b!=null?wrapped.reportOfBase[b]:null;
            entries.push({label:sa.label||('#'+(sa.index+1)),element:sa.element,
                occ:sa.occupancy,r:vi!=null?rep.perAtom[vi]:null});
        }
    }else{
        var src=CRY.baseAtoms||[];
        for(var i=0;i<src.length;i++){
            if(!neighBaseVisible(src[i]))continue;
            var vi2=wrapped.reportOfBase[i];
            entries.push({label:'#'+(i+1),element:src[i].element,
                occ:src[i].occupancy,r:vi2!=null?rep.perAtom[vi2]:null});
        }
    }
    var rows='';
    for(var e=0;e<entries.length;e++){
        var en=entries[e],r=en.r;
        var ion='—',sumT='—',valT='—',deltaT='—',cn='—',note='';
        if(r){
            ion=bvsChargeText(r);
            if(r.status==='ok'){
                sumT=r.sum.toFixed(3);
                valT=String(r.valence);
                deltaT=r.delta.toFixed(3);
                cn=String(r.contacts);
                if(r.delta>BVS_SETTLED)deltaT='<span class="bvs-meh">'+deltaT+'</span>';
            }else note=r.status==='no-contacts'?'no contacts':'no parameters';
        }else{
            var cls=BVS.classify(en.element);
            note=cls.kind==='none'?'no parameters':'not evaluated';
        }
        var occTxt=en.occ!=null&&en.occ<0.999?' <span style="opacity:.6">×'+en.occ.toFixed(2)+'</span>':'';
        rows+='<tr><td>'+escapeHtml(en.label)+'</td><td>'+escapeHtml(en.element)+occTxt+
            '</td><td>'+ion+'</td><td class="num">'+sumT+'</td><td class="num">'+valT+
            '</td><td class="num">'+deltaT+'</td><td class="num">'+cn+'</td><td>'+note+'</td></tr>';
    }
    var warn='';
    for(var w=0;w<rep.warnings.length;w++)warn+='<div class="bvs-warn">⚠ '+escapeHtml(rep.warnings[w])+'</div>';
    showModal('<h3>Bond Valence (BVS)</h3>'+
        '<div class="current-val">'+rep.counts.ok+' of '+rep.counts.atoms+' sites evaluated · parameters '+
        escapeHtml(rep.parameters)+' · pair cutoffs ≤ '+rep.maxCut.toFixed(1)+' A</div>'+
        warn+
        '<div class="bvs-scroll"><table class="bvs-table">'+
        '<tr><th>site</th><th>atom</th><th>ion</th><th>Σs</th><th>V</th><th>|Δ|</th><th>CN</th><th></th></tr>'+
        rows+'</table></div>'+
        '<div class="bvs-note">Σs sums exp((R0 − d)/B) over every periodic image inside the cutoff of each cation-anion pair. V is the formal valence the site is checked against, and CN the number of contacts counted. A |Δ| above '+
        BVS_SETTLED.toFixed(2)+' marks a site whose coordination does not support the valence it was read as. The report is read-only: bonds and coordinates are unchanged.</div>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-close">Close</button></div>',null,true);
    document.getElementById('m-close').addEventListener('click',hideModal);
}

function rebuildCrystal(){
    if(!CRY)return;
    // The base cell, its disorder filter or the cell parameters may be about
    // to change; any bond-valence report held from the previous state is no
    // longer about this structure.
    bvsInvalidate();
    // The void map is a reading of the same cells, so it is dropped for the
    // same reason: what it named may not exist in the geometry that follows.
    cavityInvalidate();
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    var baseAtoms=CRY.baseAtoms;
    var baseBonds=CRY.baseBonds;

    if(!CRY_INV){MD.atoms=[];MD.bonds=[];updateMolInfo();return}
    var inv=CRY_INV;

    // Detect the base connectivity exactly once per session. An empty list is
    // a result too (a structure the detector finds no bond in, or a cell whose
    // bonds the user has all deleted), and re-detecting then would resurrect a
    // bond an edit removed and quietly disagree with the file being saved.
    if(!CRY.bondsDetected){
        detectCrystalBaseBonds();
        baseBonds=CRY.baseBonds;
    }

    // A disorder selection shows the ordered sites plus one alternative. The
    // cell's own bond list is detected once on every site (so switching never
    // re-perceives, and a bond an edit removed stays removed) and then cut down
    // to the bonds whose two endpoints are on screen. The subset is remapped to
    // local indices for the walk below; the display keeps the original baseIdx,
    // because every commit and tooltip addresses the cell the file lists. A
    // filter that names a group no site carries any more (its last atom was
    // deleted) falls back to All here, before anything is placed.
    if(CRY.disorderFilter){
        var vgOk=false;
        for(var vg=0;vg<baseAtoms.length;vg++){
            if(baseAtoms[vg].disorderGroup===CRY.disorderFilter){vgOk=true;break}
        }
        if(!vgOk)CRY.disorderFilter=null;
    }
    if(CRY.disorderFilter){
        var visMap=new Array(baseAtoms.length);
        var visAtoms=[];
        for(var va=0;va<baseAtoms.length;va++){
            var baV=baseAtoms[va];
            if(!baV.disorderGroup||baV.disorderGroup===CRY.disorderFilter){
                visMap[va]=visAtoms.length;visAtoms.push(baV);
            }else visMap[va]=-1;
        }
        var visBonds=[];
        for(var vb=0;vb<baseBonds.length;vb++){
            var bbV=baseBonds[vb];
            var m1=visMap[bbV.atom1],m2=visMap[bbV.atom2];
            if(m1>=0&&m2>=0)visBonds.push({atom1:m1,atom2:m2,order:bbV.order,shift:bbV.shift});
        }
        baseAtoms=visAtoms;baseBonds=visBonds;
    }

    // ===== Connected components, unwrapped =====
    //
    // The old grouping test was "does any bond of this component cross a cell
    // face" (hasShift). In a molecular crystal the file wraps every atom into
    // the cell first, so a molecule that straddles a face owns such a bond
    // whatever its chemistry — CCDC 142619 put all ten of its molecules on the
    // wrong side of that test and every one of them was drawn cut. The test
    // that actually separates a molecule from an infinite chain is whether the
    // component's bond cycles are consistent: walking any cycle of real bonds
    // must come back to the same periodic image. So each component is walked
    // once with an image offset per atom (the unwrap), and a cycle that comes
    // back shifted marks the whole component periodic.
    //
    //   finite component  — a molecule: copied whole, see pass 2 below;
    //   periodic component (diamond, a polymer chain) — the old behaviour:
    //                       per-atom boundary test and bonds cut at the face,
    //                       because drawing it "whole" has no meaning.
    //
    // Bond (a1,a2,shift s) means the image of a2 in cell c+s bonds with the
    // image of a1 in cell c, so walking the edge a1->a2 adds s to the image
    // offset and a2->a1 subtracts it. Offsets are exact integers, so the
    // consistency check cannot drift. The walk itself lives in
    // crystalComponents(), because the bond detector lays a finite molecule out
    // by the same offsets and the two must not disagree.
    var comp=crystalComponents(baseBonds,baseAtoms.length);
    var groupOf=comp.groupOf;var groups=comp.groups;

    var iMin=Math.floor(aMin),iMax=Math.ceil(aMax);
    var jMin=Math.floor(bMin),jMax=Math.ceil(bMax);
    var kMin=Math.floor(cMin),kMax=Math.ceil(cMax);
    var EPS=1e-9;

    // The whole-molecule copy test expands the boundary by this margin
    // (Angstrom, converted per axis with that lattice vector's length), so a
    // molecule whose centre sits just outside can still be drawn complete.
    var marginA=CRYSTAL_COPY_MARGIN/Math.sqrt(lv[0][0]*lv[0][0]+lv[0][1]*lv[0][1]+lv[0][2]*lv[0][2]);
    var marginB=CRYSTAL_COPY_MARGIN/Math.sqrt(lv[1][0]*lv[1][0]+lv[1][1]*lv[1][1]+lv[1][2]*lv[1][2]);
    var marginC=CRYSTAL_COPY_MARGIN/Math.sqrt(lv[2][0]*lv[2][0]+lv[2][1]*lv[2][1]+lv[2][2]*lv[2][2]);
    // Anchor cells are scanned a cell past the boundary on every side: the
    // centroid test does the actual filtering, and the extra ring only makes
    // the scan independent of where unwrapping put each molecule's coordinates.
    var extA=1+Math.ceil(marginA),extB=1+Math.ceil(marginB),extC=1+Math.ceil(marginC);

    var atoms=[];
    var bonds=[];
    var idx=0;
    var cellMap={};
    // base atom index -> its placed image cells, in placement order. The bond
    // pass walks these instead of the boundary's cell range, because a whole
    // molecule's atoms may sit outside it.
    var placedCells={};
    function rememberCell(bi,ci,cj,ck){
        if(!placedCells[bi])placedCells[bi]=[];
        placedCells[bi].push([ci,cj,ck]);
    }

    // Pass 1 — periodic components and unbonded atoms: the per-image boundary
    // test, exactly as it always was. Copying a network "whole" has no meaning,
    // and a lone atom has nothing to complete.
    baseAtoms.forEach(function(ba,bi){
        var g=groupOf[bi];
        if(g>=0&&!groups[g].periodic)return;
        var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
        var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
        var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
        // The integer cell offset this image sits at, measured from the base
        // site's own (unwrapped) fractional position. It is the one number
        // commitCrystalCoords() needs to turn a move of this image back into a
        // move of the site: display - base - imgShift is the edit alone.
        var floorX=Math.floor(fx),floorY=Math.floor(fy),floorZ=Math.floor(fz);
        fx-=floorX;fy-=floorY;fz-=floorZ;
        for(var ci=iMin;ci<=iMax;ci++){
            for(var cj=jMin;cj<=jMax;cj++){
                for(var ck=kMin;ck<=kMax;ck++){
                    var ax=fx+ci,ay=fy+cj,az=fz+ck;
                    // Half-open on the max side, the same rule the whole-molecule
                    // pass uses: the image of an atom that sits exactly on the
                    // upper face belongs to the cell below it, not to both. The
                    // closed test (<= aMax) drew that one site twice — at 0 and at
                    // 1 — for a cell that states it once.
                    if(ax<aMin-EPS||ax>=aMax-1e-12)continue;
                    if(ay<bMin-EPS||ay>=bMax-1e-12)continue;
                    if(az<cMin-EPS||az>=cMax-1e-12)continue;
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
                        cellI:ci,cellJ:cj,cellK:ck,
                        imgShift:[ci-floorX,cj-floorY,ck-floorZ],
                        uAniso:ba.uAniso
                    });
                    cellMap[bi+'_'+ci+'_'+cj+'_'+ck]=idx;
                    rememberCell(bi,ci,cj,ck);
                    idx++;
                }
            }
        }
    });

    // Pass 2 — finite molecules, whole: every image whose unwrapped centroid
    // lands inside the boundary (plus the margin) is placed complete, atoms
    // outside included, so a molecule straddling a face is drawn as the
    // molecule it is instead of as two fragments. The BFS offsets already put
    // the component's atoms into one image; the centroid is that image's
    // centre, and the whole component is shifted by the same integer vector so
    // the centroid sits in [0,1)^3 and the anchor scan can stay within the
    // boundary's own cell range. Each atom's cell label records which periodic
    // image of its base site is on screen, so different images of one site stay
    // tellable apart and a selection can find the same atom again.
    groups.forEach(function(g){
        if(g.periodic)return;
        var members=g.members,n=members.length;
        var u={};
        var ucx=0,ucy=0,ucz=0;
        for(var t=0;t<n;t++){
            var bi=members[t];
            var ba=baseAtoms[bi];
            var o=g.offset[bi];
            var ux=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z+o[0];
            var uy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z+o[1];
            var uz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z+o[2];
            u[bi]=[ux,uy,uz];
            ucx+=ux;ucy+=uy;ucz+=uz;
        }
        ucx/=n;ucy/=n;ucz/=n;
        var tx=Math.floor(ucx),ty=Math.floor(ucy),tz=Math.floor(ucz);
        ucx-=tx;ucy-=ty;ucz-=tz;
        for(var t2=0;t2<n;t2++){
            var uu=u[members[t2]];
            uu[0]-=tx;uu[1]-=ty;uu[2]-=tz;
        }
        for(var ci=iMin-extA;ci<=iMax+extA;ci++){
            for(var cj=jMin-extB;cj<=jMax+extB;cj++){
                for(var ck=kMin-extC;ck<=kMax+extC;ck++){
                    var cx=ucx+ci,cy=ucy+cj,cz=ucz+ck;
                    // Half-open on the max side: a molecule whose centre sits
                    // exactly on the upper face (a centrosymmetric molecule on
                    // 0, for instance) is one image, not two — the closed test
                    // would draw the face molecule at both ends of the range.
                    if(cx<aMin-marginA-EPS||cx>=aMax+marginA-EPS)continue;
                    if(cy<bMin-marginB-EPS||cy>=bMax+marginB-EPS)continue;
                    if(cz<cMin-marginC-EPS||cz>=cMax+marginC-EPS)continue;
                    for(var t3=0;t3<n;t3++){
                        var bi3=members[t3];
                        var ba3=baseAtoms[bi3];
                        var uu3=u[bi3];
                        var ax=uu3[0]+ci,ay=uu3[1]+cj,az=uu3[2]+ck;
                        var ox=ax*lv[0][0]+ay*lv[1][0]+az*lv[2][0];
                        var oy=ax*lv[0][1]+ay*lv[1][1]+az*lv[2][1];
                        var oz=ax*lv[0][2]+ay*lv[1][2]+az*lv[2][2];
                        var o3=g.offset[bi3];
                        var ci3=o3[0]-tx+ci,cj3=o3[1]-ty+cj,ck3=o3[2]-tz+ck;
                        var color=ba3.color||((MD.atomColors&&MD.atomColors[ba3.element])||'#FF1493');
                        atoms.push({
                            element:ba3.element,
                            x:ox,y:oy,z:oz,
                            index:idx,
                            occupancy:ba3.occupancy,
                            color:color,
                            baseIdx:ba3.baseIdx!=null?ba3.baseIdx:bi3,
                            cellI:ci3,cellJ:cj3,cellK:ck3,
                            // The whole-copy placement puts the atom at
                            // base + (offset - t + cell), so that integer is its
                            // image offset (see commitCrystalCoords()).
                            imgShift:[ci3,cj3,ck3],
                            uAniso:ba3.uAniso
                        });
                        cellMap[bi3+'_'+ci3+'_'+cj3+'_'+ck3]=idx;
                        rememberCell(bi3,ci3,cj3,ck3);
                        idx++;
                    }
                }
            }
        }
    });

    // Bonds: for every bond, walk the images actually placed (placedCells)
    // rather than the boundary's cell range — a whole molecule's atoms can sit
    // outside it. A finite component's offsets are consistent by construction
    // (offset of the bonded image = offset of this one + shift), so every one
    // of its bonds finds both endpoints and none is cut; only a periodic
    // component's bonds lose an endpoint at the boundary and become split
    // bonds (shift tells the renderer to draw the stub towards the image).
    var seenPairs={};
    var seenSplits={};
    baseBonds.forEach(function(bb){
        var b1=bb.atom1,b2=bb.atom2;
        var s=bb.shift||ZERO_SHIFT;
        var cells1=placedCells[b1]||[];
        var cells2=placedCells[b2]||[];
        for(var p=0;p<cells1.length;p++){
            var c1=cells1[p];
            var idx1=cellMap[b1+'_'+c1[0]+'_'+c1[1]+'_'+c1[2]];
            if(idx1==null)continue;
            var idx2=cellMap[b2+'_'+(c1[0]+s[0])+'_'+(c1[1]+s[1])+'_'+(c1[2]+s[2])];
            if(idx2!=null){
                var key=idx1<idx2?idx1+'_'+idx2:idx2+'_'+idx1;
                if(!seenPairs[key]){
                    seenPairs[key]=true;
                    bonds.push({atom1:idx1,atom2:idx2,order:bb.order});
                }
            }else{
                var idx2s=cellMap[b2+'_'+c1[0]+'_'+c1[1]+'_'+c1[2]];
                if(idx2s!=null){
                    var skey='F_'+idx1+'_'+idx2s+'_'+s[0]+'_'+s[1]+'_'+s[2];
                    if(!seenSplits[skey]){
                        seenSplits[skey]=true;
                        bonds.push({atom1:idx1,atom2:idx2s,order:bb.order,shift:s});
                    }
                }
            }
        }
        if(s[0]||s[1]||s[2]){
            var rs=[-s[0],-s[1],-s[2]];
            for(var q=0;q<cells2.length;q++){
                var c2=cells2[q];
                var idx2r=cellMap[b2+'_'+c2[0]+'_'+c2[1]+'_'+c2[2]];
                if(idx2r==null)continue;
                var ridx1=cellMap[b1+'_'+(c2[0]+rs[0])+'_'+(c2[1]+rs[1])+'_'+(c2[2]+rs[2])];
                if(ridx1==null){
                    var idx1r=cellMap[b1+'_'+c2[0]+'_'+c2[1]+'_'+c2[2]];
                    if(idx1r!=null){
                        var skey2='R_'+idx2r+'_'+idx1r+'_'+rs[0]+'_'+rs[1]+'_'+rs[2];
                        if(!seenSplits[skey2]){
                            seenSplits[skey2]=true;
                            bonds.push({atom1:idx2r,atom2:idx1r,order:bb.order,shift:rs});
                        }
                    }
                }
            }
        }
    });

    // Z for a molecular crystal: the whole molecules the cell contains — the
    // finite components plus the unbonded single atoms. A cell whose only
    // component is periodic leaves the count at zero and the panel leaves the
    // line out.
    var nMol=0;
    for(var gi=0;gi<groups.length;gi++){if(!groups[gi].periodic)nMol++}
    for(var bl=0;bl<baseAtoms.length;bl++){if(groupOf[bl]===-1)nMol++}
    CRY.moleculeCount=nMol;

    MD.atoms=atoms;
    MD.bonds=bonds;
    syncDisorderPanel();
    syncSymmetryNote();
    updateMolInfo();
    needsRender=true;
}

// ===== Writing an edit of the picture back into the base cell =====
//
// A crystal's file content is its base cell, while the picture is a supercell:
// one base site may be on screen several times, at an integer number of cells
// from where the file lists it. Saving CIF/VESTA writes baseAtoms, and any
// rebuildCrystal() (a bounds change, an undo, a removed atom) redraws from it —
// so a coordinate edit that stays in MD.atoms alone is lost by both, and the
// same session would save two different structures depending on the format.
//
// rebuildCrystal() records each image's integer cell offset as imgShift, so
// the displacement the edit made to one image is exactly
//
//   delta_frac = frac(image) − frac(base) − imgShift
//
// and adding it to the base site carries the edit into the cell. A site that is
// on screen more than once (a 2x2x2 supercell shows each site eight times) gets
// the delta of its **most-moved** image: an edit moves one image's copy, or
// translates a group of them by the same vector, while the untouched copies all
// report zero — summing them would multiply the edit by the number of images.
// displayIdxList names the atoms the edit actually moved. A session that
// commits one step at a time (Rotate Group) keeps the other images at the
// geometry of an earlier commit, and for those the apparent delta is the
// previous step, not this one, so a step must only write back its own atoms.
function commitCrystalCoords(displayIdxList){
    if(!CRY||!CRY_INV)return;
    // Coordinates are about to move; a sum computed before the move describes
    // the old geometry.
    bvsInvalidate();
    cavityInvalidate();
    var inv=CRY_INV;
    var list=displayIdxList||null;
    var n=list?list.length:MD.atoms.length;
    var best={};
    for(var k=0;k<n;k++){
        var i=list?list[k]:k;
        var a=MD.atoms[i];
        if(!a||!a.imgShift||a.baseIdx==null)continue;
        var ba=CRY.baseAtoms[a.baseIdx];
        if(!ba)continue;
        var ax=a.x-ba.x,ay=a.y-ba.y,az=a.z-ba.z;
        var dx=inv[0][0]*ax+inv[0][1]*ay+inv[0][2]*az-a.imgShift[0];
        var dy=inv[1][0]*ax+inv[1][1]*ay+inv[1][2]*az-a.imgShift[1];
        var dz=inv[2][0]*ax+inv[2][1]*ay+inv[2][2]*az-a.imgShift[2];
        var d2=dx*dx+dy*dy+dz*dz;
        // A site that sits exactly where it was built moves nothing; skipping it
        // keeps a run of commits from creeping the base by round-off.
        if(d2<1e-18)continue;
        var cur=best[a.baseIdx];
        if(!cur||d2>cur.d2)best[a.baseIdx]={dx:dx,dy:dy,dz:dz,d2:d2};
    }
    // The base is written back through the same inverse the delta was read
    // with, so cart → frac → cart closes on the value the display already had
    // (to the last bits of the double).
    var lv=CRY.latticeVectors;
    var keys=Object.keys(best);
    for(var q=0;q<keys.length;q++){
        var key=keys[q];
        var b=CRY.baseAtoms[key];
        if(!b)continue;
        var d=best[key];
        var fx=inv[0][0]*b.x+inv[0][1]*b.y+inv[0][2]*b.z+d.dx;
        var fy=inv[1][0]*b.x+inv[1][1]*b.y+inv[1][2]*b.z+d.dy;
        var fz=inv[2][0]*b.x+inv[2][1]*b.y+inv[2][2]*b.z+d.dz;
        b.x=fx*lv[0][0]+fy*lv[1][0]+fz*lv[2][0];
        b.y=fx*lv[0][1]+fy*lv[1][1]+fz*lv[2][1];
        b.z=fx*lv[0][2]+fy*lv[1][2]+fz*lv[2][2];
    }
}

// The bond a crystal state actually holds, for a bond drawn between two images:
// its endpoints are the two sites, and its shift is the number of cells between
// the images the two atoms sit in (a split bond adds the shift of the stub it
// was drawn as). A bond object created by hand (no image labels yet) is read
// from the two display atoms it is about to join.
function crystalBaseBondFromDisplay(b){
    if(!CRY||!b)return null;
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2)return null;
    if(a1.baseIdx==null||a2.baseIdx==null)return null;
    var s=b.shift||ZERO_SHIFT;
    var c1x=a1.cellI||0,c1y=a1.cellJ||0,c1z=a1.cellK||0;
    var c2x=a2.cellI||0,c2y=a2.cellJ||0,c2z=a2.cellK||0;
    return {atom1:a1.baseIdx,atom2:a2.baseIdx,
        shift:[c2x+s[0]-c1x,c2y+s[1]-c1y,c2z+s[2]-c1z]};
}
// The base bond equivalent to the one a display bond describes, or null.
// Direction is not part of a bond's identity, so the reversed spelling with the
// opposite shift is the same bond.
function crystalFindBaseBond(desc){
    if(!CRY||!desc)return null;
    var s=desc.shift;
    for(var i=0;i<CRY.baseBonds.length;i++){
        var x=CRY.baseBonds[i];
        var xs=x.shift||ZERO_SHIFT;
        if(x.atom1===desc.atom1&&x.atom2===desc.atom2&&xs[0]===s[0]&&xs[1]===s[1]&&xs[2]===s[2])return x;
        if(x.atom1===desc.atom2&&x.atom2===desc.atom1&&xs[0]===-s[0]&&xs[1]===-s[1]&&xs[2]===-s[2])return x;
    }
    return null;
}
// A bond drawn as a stub across a cell face: the two atoms on screen are a whole
// cell apart on purpose, so the distance between them is not the bond's length
// and editing it would write back a displacement the bond does not have.
function crystalCrossFaceBondBetween(i1,i2){
    if(!CRY)return null;
    for(var i=0;i<MD.bonds.length;i++){
        var b=MD.bonds[i];
        if(!((b.atom1===i1&&b.atom2===i2)||(b.atom1===i2&&b.atom2===i1)))continue;
        var s=b.shift||ZERO_SHIFT;
        if(s[0]||s[1]||s[2])return b;
    }
    return null;
}
// Direction from "from" to the periodic image of "to" that a bond shift names,
// as a unit vector. A molecule that straddles a cell face holds its neighbours a
// whole cell away in the file, and the raw difference would point a new atom
// straight through the wall instead of away from its neighbours.
function crystalImageDirection(from,to,sign,shift){
    var dx=to.x-from.x,dy=to.y-from.y,dz=to.z-from.z;
    if(CRY_INV){
        var inv=CRY_INV,lv=CRY.latticeVectors;
        var fx=inv[0][0]*(to.x-from.x)+inv[0][1]*(to.y-from.y)+inv[0][2]*(to.z-from.z)+sign*shift[0];
        var fy=inv[1][0]*(to.x-from.x)+inv[1][1]*(to.y-from.y)+inv[1][2]*(to.z-from.z)+sign*shift[1];
        var fz=inv[2][0]*(to.x-from.x)+inv[2][1]*(to.y-from.y)+inv[2][2]*(to.z-from.z)+sign*shift[2];
        // The nearest image of that neighbour is the one the anchor is bonded
        // to as far as the free direction is concerned — found by the same
        // enumerated search the bonds use, because rounding the fractional
        // difference picks the nearest image only in an orthogonal cell.
        var cx=fx*lv[0][0]+fy*lv[1][0]+fz*lv[2][0];
        var cy=fx*lv[0][1]+fy*lv[1][1]+fz*lv[2][1];
        var cz=fx*lv[0][2]+fy*lv[1][2]+fz*lv[2][2];
        var near=crystalShortestFractional(fx,fy,fz,Math.sqrt(cx*cx+cy*cy+cz*cz),lv,crystalPlaneWidths(lv));
        fx=near.frac[0];fy=near.frac[1];fz=near.frac[2];
        dx=fx*lv[0][0]+fy*lv[1][0]+fz*lv[2][0];
        dy=fx*lv[0][1]+fy*lv[1][1]+fz*lv[2][1];
        dz=fx*lv[0][2]+fy*lv[1][2]+fz*lv[2][2];
    }
    var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(l<1e-10)return null;
    return {x:dx/l,y:dy/l,z:dz/l};
}

// Put a crystal's base cell back the way the snapshot on top of the undo stack
// recorded it. A dialog's Cancel calls this just before it pops that snapshot,
// so the base state it hands back is not one the stack keeps a reference to.
function restoreCrystalBaseFromUndoTop(){
    if(!CRY||undoStack.length===0)return;
    var snap=undoStack[undoStack.length-1];
    if(!snap.baseAtoms||!snap.baseBonds)return;
    CRY.baseAtoms=snap.baseAtoms;
    CRY.baseBonds=snap.baseBonds;
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
            index:newIdx,occupancy:newOcc,baseIdx:newIdx,color:color,
            siteLabel:ba.siteLabel,asymIdx:ba.asymIdx,symOp:ba.symOp,
            constraint:ba.constraint,disorderGroup:ba.disorderGroup,
            uAniso:ba.uAniso
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

var cellWireframe=null;      // the bounds' outer box
var cellGridLines=null;      // the interior unit-cell division lines
var cellEdgeLabels=[];       // the a/b/c sprites on the low corner's three edges
// The cell panel's (hkl) calculator and the one plane it can draw. The state
// lives up here because rebuildScene refreshes the plane, and that runs before
// the section that defines the calculator (see "Lattice planes" below).
var hklPlane={h:0,k:0,l:1};
var hklPlane2={h:1,k:0,l:0};
var hklShown=false;
var hklPlaneMesh=null;
var hklPlaneBorder=null;
// One text sprite per call: the letter is drawn onto a small canvas (no font
// file, no image load — the page's CSP allows nothing else) and the sprite
// always faces the camera.
function makeEdgeLabelSprite(text,color){
    var canvas=document.createElement('canvas');
    canvas.width=64;canvas.height=64;
    var ctx=canvas.getContext('2d');
    ctx.font='bold 44px sans-serif';
    ctx.textAlign='center';
    ctx.textBaseline='middle';
    ctx.fillStyle=color;
    ctx.fillText(text,32,34);
    var tex=new THREE.CanvasTexture(canvas);
    var mat=new THREE.SpriteMaterial({map:tex,transparent:true,depthTest:false,depthWrite:false});
    return new THREE.Sprite(mat);
}
function buildCellWireframe(){
    if(cellWireframe){moleculeGroup.remove(cellWireframe);disposeMesh(cellWireframe);cellWireframe=null}
    if(cellGridLines){moleculeGroup.remove(cellGridLines);disposeMesh(cellGridLines);cellGridLines=null}
    cellEdgeLabels.forEach(function(sp){
        moleculeGroup.remove(sp);
        // The sprite's own plane geometry is three.js's shared one — only the
        // material and its canvas texture are this sprite's to release.
        if(sp.material){if(sp.material.map)sp.material.map.dispose();sp.material.dispose()}
    });
    cellEdgeLabels=[];
    if(!CRY)return;
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    // One corner of the subdivision grid, in world space relative to the centre.
    function cellPt(ai,aj,ak){
        var x=ai*lv[0][0]+aj*lv[1][0]+ak*lv[2][0];
        var y=ai*lv[0][1]+aj*lv[1][1]+ak*lv[2][1];
        var z=ai*lv[0][2]+aj*lv[1][2]+ak*lv[2][2];
        return new THREE.Vector3(x-CX,y-CY,z-CZ);
    }
    // The planes the supercell is divided at: its two faces per axis plus every
    // integer plane between them — the unit-cell boundaries themselves. 0..1
    // therefore yields exactly the twelve box edges this always drew, and a
    // 3x3x3 supercell shows its 27 cells instead of one outer shell.
    function gridVals(mn,mx){
        var vals=[mn,mx];
        for(var v=Math.ceil(mn+1e-9);v<=mx-1e-9;v++)vals.push(v);
        return vals;
    }
    function onFace(v,mn,mx){return Math.abs(v-mn)<1e-9||Math.abs(v-mx)<1e-9}
    var aG=gridVals(aMin,aMax),bG=gridVals(bMin,bMax),cG=gridVals(cMin,cMax);
    var outer=[],inner=[];
    function addLine(p0,p1,bucket){
        bucket.push(p0.x,p0.y,p0.z,p1.x,p1.y,p1.z);
    }
    // Lines run along one axis at every grid plane of the other two; only the
    // lines whose both cross planes are boundary faces belong to the box.
    bG.forEach(function(bj){
        cG.forEach(function(ck){
            addLine(cellPt(aMin,bj,ck),cellPt(aMax,bj,ck),(onFace(bj,bMin,bMax)&&onFace(ck,cMin,cMax))?outer:inner);
        });
    });
    aG.forEach(function(ai){
        cG.forEach(function(ck){
            addLine(cellPt(ai,bMin,ck),cellPt(ai,bMax,ck),(onFace(ai,aMin,aMax)&&onFace(ck,cMin,cMax))?outer:inner);
        });
    });
    aG.forEach(function(ai){
        bG.forEach(function(bj){
            addLine(cellPt(ai,bj,cMin),cellPt(ai,bj,cMax),(onFace(ai,aMin,aMax)&&onFace(bj,bMin,bMax))?outer:inner);
        });
    });
    function lineSegmentsFor(arr){
        if(arr.length===0)return null;
        var geo=new THREE.BufferGeometry();
        geo.setAttribute('position',new THREE.Float32BufferAttribute(arr,3));
        return geo;
    }
    var outerGeo=lineSegmentsFor(outer);
    if(outerGeo){
        cellWireframe=new THREE.LineSegments(outerGeo,new THREE.LineBasicMaterial({color:0x88aaff}));
        moleculeGroup.add(cellWireframe);
    }
    var innerGeo=lineSegmentsFor(inner);
    if(innerGeo){
        // The interior division lines stay one step quieter than the outer box,
        // so a 3x3x3 grid reads as structure rather than as a cage over atoms.
        cellGridLines=new THREE.LineSegments(innerGeo,new THREE.LineBasicMaterial({color:0x88aaff,transparent:true,opacity:0.3}));
        moleculeGroup.add(cellGridLines);
    }
    // a/b/c on the three edges that meet at the bounds' low corner, nudged
    // outside the box along the other two axes. The colours are the axes
    // indicator's (a red, b green, c blue), so both readings of "which
    // direction is a" agree.
    var la=Math.sqrt(lv[0][0]*lv[0][0]+lv[0][1]*lv[0][1]+lv[0][2]*lv[0][2]);
    var lb=Math.sqrt(lv[1][0]*lv[1][0]+lv[1][1]*lv[1][1]+lv[1][2]*lv[1][2]);
    var lc=Math.sqrt(lv[2][0]*lv[2][0]+lv[2][1]*lv[2][1]+lv[2][2]*lv[2][2]);
    var minEdge=Math.min(la,lb,lc);
    if(minEdge>1e-9){
        var ua=[lv[0][0]/la,lv[0][1]/la,lv[0][2]/la];
        var ub=[lv[1][0]/lb,lv[1][1]/lb,lv[1][2]/lb];
        var uc=[lv[2][0]/lc,lv[2][1]/lc,lv[2][2]/lc];
        var dOut=minEdge*0.12;
        var labS=minEdge*0.16;
        function edgeLabel(text,color,mid,off){
            var sp=makeEdgeLabelSprite(text,color);
            sp.position.set(mid.x+off[0],mid.y+off[1],mid.z+off[2]);
            sp.scale.set(labS,labS,1);
            moleculeGroup.add(sp);
            cellEdgeLabels.push(sp);
        }
        edgeLabel('a','#FF4444',cellPt((aMin+aMax)/2,bMin,cMin),[-dOut*(ub[0]+uc[0]),-dOut*(ub[1]+uc[1]),-dOut*(ub[2]+uc[2])]);
        edgeLabel('b','#44FF44',cellPt(aMin,(bMin+bMax)/2,cMin),[-dOut*(ua[0]+uc[0]),-dOut*(ua[1]+uc[1]),-dOut*(ua[2]+uc[2])]);
        edgeLabel('c','#4488FF',cellPt(aMin,bMin,(cMin+cMax)/2),[-dOut*(ua[0]+ub[0]),-dOut*(ua[1]+ub[1]),-dOut*(ua[2]+ub[2])]);
    }
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
// An atom's appearance is what its *material* decides, and a material cannot
// vary inside one instanced draw: the element colour rides the instance colour
// buffer, an occupancy rides a texture, and a ghost is translucent. The
// appearance key below is that classification. The highlight is a third axis —
// plain, selected and difference-marked are separate instanced meshes (see the
// atom sphere layer), so selecting one atom moves one instance between two
// buffers and can never change how another atom is drawn.
//
// A ghost is a fourth appearance: faint, and drawn without writing depth, so a
// ghost standing exactly on a real atom (the usual counterpoise layout) reads as
// a shadow *around* that atom instead of hiding it.
function atomInstanceMaterial(appearance,state){
    return pooledMat('ia|'+appearance+'|'+state,function(){
        var m;
        if(appearance==='g'){
            m=new THREE.MeshPhongMaterial({color:new THREE.Color(GHOST_COLOR),shininess:20,specular:0x222222,transparent:true,opacity:0.45,depthWrite:false});
        }else if(appearance.charAt(0)==='o'){
            var bits=appearance.split('|');
            m=new THREE.MeshPhongMaterial({color:0xffffff,map:getOccupancyTexture(bits[1],parseFloat(bits[2])),shininess:80,specular:0x444444});
        }else{
            m=new THREE.MeshPhongMaterial({color:0xffffff,shininess:80,specular:0x444444});
        }
        if(state===1){m.emissive=new THREE.Color(0xffff00);m.emissiveIntensity=0.6}
        else if(state===2){m.emissive=new THREE.Color(0xff6600);m.emissiveIntensity=0.4}
        return m;
    });
}
// Ask for one of those three states. Every writer of the highlight goes through
// here; the handle records the state and moves its own instance between the two
// state meshes of its appearance.
function setAtomHighlight(mesh,state){
    if(!mesh||!mesh.userData)return;
    var ud=mesh.userData;
    var old=ud.hlState||0;
    if(old===state)return;
    ud.hlState=state;
    if(ud.sphereLayer){
        sphereLayerWriteState(ud.sphereLayer,mesh,old,false);
        sphereLayerWriteState(ud.sphereLayer,mesh,state,true);
    }
}
// The comparison's atoms follow the same rule. The difference marker is the
// colour the difference half writes into the instance colour buffer (the
// comparison's orange), not a second material.
function setDiffAtomHighlight(mesh,state){
    setAtomHighlight(mesh,state);
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

// ===== Atoms: one instanced sphere per appearance and state =====
// The bonds have been instanced for a long time; the atoms were not, so a
// crystal's 2x2x2 view was 2352 draw calls for its spheres alone. A sphere is
// now a shared unit geometry scaled per instance, and every atom is a *record*
// — the atomMeshes array holds a lightweight handle (position + userData, never
// added to the scene) instead of a Mesh. Every other pass keeps the contract it
// always had: "position" is the atom minus the rotation centre, "userData"
// carries the element, index and highlight state, "visible" says whether a
// sphere is drawn at all.
//
// The material cannot vary inside one instanced draw, so atoms are grouped by
// the appearance its material decides (ordinary / ghost / occupancy) and, within
// that, by highlight state (plain / selected / difference-marked) into three
// instanced meshes. The three state meshes hold one slot per atom of the
// appearance, so a highlight change is two setMatrixAt writes — it never
// re-indexes the buffers, never touches another atom, and never changes the
// number of draw calls.
function atomHandle(){
    return {
        position:{x:0,y:0,z:0,
            set:function(x,y,z){this.x=x;this.y=y;this.z=z},
            copy:function(v){this.x=v.x;this.y=v.y;this.z=v.z}},
        userData:{},
        visible:true
    };
}
var mainSphereLayer={kind:'main',parent:null,classes:{},count:0};
var diffSphereLayer={kind:'diff',parent:null,classes:{},count:0};
var atomColorCache={};
var ATOM_WHITE=new THREE.Color(0xffffff);
var ATOM_ORANGE=new THREE.Color(0xff6600);
var _atomM=new THREE.Matrix4();
var _atomZero=new THREE.Matrix4().makeScale(0,0,0);
function atomColorOf(col){
    var c=atomColorCache[col];
    if(!c){c=new THREE.Color(col);atomColorCache[col]=c}
    return c;
}
function unitSphereGeo(segs){
    return pooledGeo('us|'+segs[0]+'|'+segs[1],function(){return new THREE.SphereGeometry(1,segs[0],segs[1])});
}
// ===== Thermal ellipsoids (ORTEP) =====
//
// A site with an anisotropic displacement tensor is drawn as the ellipsoid
// that contains 50% of its displacement probability: the surface
// x^T U^-1 x = chi2(3, 0.5), whose semi-axes are
// sqrt(chi2) * sqrt(eigenvalue of U) along U's eigenvectors. The viewer's unit
// sphere is mapped by a matrix whose columns are those three axes, so one
// instanced draw keeps drawing every atom; only the per-instance matrix
// differs from a plain sphere.
//
// The tensor arrives in Cartesian Angstrom squared, already rotated into the
// frame of the image it belongs to (see the CIF reader), so nothing here
// depends on the crystal's symmetry. Eigenvectors are found by Jacobi
// rotations: for a symmetric 3x3 they converge in a handful of sweeps, and no
// closed-form cubic root is needed.
var ORTEP_K=Math.sqrt(2.365973884); // chi-square, 3 degrees of freedom, 50%
var ortepEnabled=false;
var ortepMatrixCache={};
function anisoEigen(u){
    var a=[[u[0],u[3],u[4]],[u[3],u[1],u[5]],[u[4],u[5],u[2]]];
    var v=[[1,0,0],[0,1,0],[0,0,1]];
    for(var sweep=0;sweep<12;sweep++){
        var off=Math.abs(a[0][1])+Math.abs(a[0][2])+Math.abs(a[1][2]);
        if(off<1e-15)break;
        for(var p=0;p<2;p++){
            for(var q=p+1;q<3;q++){
                var apq=a[p][q];
                if(Math.abs(apq)<1e-18)continue;
                var theta=(a[q][q]-a[p][p])/(2*apq);
                var t=(theta>=0?1:-1)/(Math.abs(theta)+Math.sqrt(theta*theta+1));
                var c=1/Math.sqrt(t*t+1),s=t*c;
                for(var k=0;k<3;k++){
                    var akp=a[k][p],akq=a[k][q];
                    a[k][p]=c*akp-s*akq;
                    a[k][q]=s*akp+c*akq;
                }
                for(k=0;k<3;k++){
                    var apk=a[p][k],aqk=a[q][k];
                    a[p][k]=c*apk-s*aqk;
                    a[q][k]=s*apk+c*aqk;
                }
                for(k=0;k<3;k++){
                    var vkp=v[k][p],vkq=v[k][q];
                    v[k][p]=c*vkp-s*vkq;
                    v[k][q]=s*vkp+c*vkq;
                }
            }
        }
    }
    var order=[0,1,2].sort(function(i,j){return a[j][j]-a[i][i]});
    return {
        val:[a[order[0]][order[0]],a[order[1]][order[1]],a[order[2]][order[2]]],
        vec:[[v[0][order[0]],v[0][order[1]],v[0][order[2]]],
             [v[1][order[0]],v[1][order[1]],v[1][order[2]]],
             [v[2][order[0]],v[2][order[1]],v[2][order[2]]]]
    };
}
// Whether a tensor can draw an ellipsoid at all. The trace of U is the sum of
// its eigenvalues, so a positive trace means at least one axis has a positive
// length; a tensor with a non-positive trace would collapse to a point.
function isAnisoUsable(u){
    if(!u||u.length<6)return false;
    for(var i=0;i<6;i++)if(typeof u[i]!=='number'||!isFinite(u[i]))return false;
    return u[0]+u[1]+u[2]>0;
}
function anisoEllipsoid(u){
    if(!isAnisoUsable(u))return null;
    var key=u.join(',');
    var hit=ortepMatrixCache[key];
    if(hit!==undefined)return hit;
    var e=anisoEigen(u);
    var axes=[],maxR=0;
    for(var i=0;i<3;i++){
        // A non-positive eigenvalue is not an ellipsoid along that axis; it is
        // also not a reason to hide the atom. The axis collapses and the atom
        // still shows the two dimensions the file did state.
        var lam=e.val[i]>0?e.val[i]:0;
        var s=ORTEP_K*Math.sqrt(lam);
        if(s>maxR)maxR=s;
        axes.push(new THREE.Vector3(e.vec[0][i]*s,e.vec[1][i]*s,e.vec[2][i]*s));
    }
    var res=maxR>1e-6?{m:new THREE.Matrix4().makeBasis(axes[0],axes[1],axes[2]),r:maxR}:null;
    if(Object.keys(ortepMatrixCache).length>4096)ortepMatrixCache={};
    ortepMatrixCache[key]=res;
    return res;
}
// What the atom's material must know: a ghost is translucent whatever its
// element, an occupancy below 1 carries its own texture, everything else is an
// ordinary element sphere whose colour rides the instance buffer. An atom drawn
// as a thermal ellipsoid gets a suffix of its own, so the class it lands in is
// rebuilt when the ellipsoid is switched on or off and no stale matrix from the
// other mode can survive in a buffer.
function atomAppearanceOf(ud){
    var key;
    if(ud.ghost)key='g';
    else if(ud.occ!=null&&ud.occ<0.999)key='o|'+ud.baseColor+'|'+ud.occ.toFixed(3);
    else key='n';
    // A ghost is a marker, never a surface: even a tensor on it must not move
    // it into the ellipsoid class, whose matrix would then ignore the marker.
    if(ud.anisoM&&!ud.ghost)key+='|a';
    return key;
}
function writeAtomMatrix(inst,slot,h,r){
    if(!inst)return;
    var ud=h.userData;
    if(r>0&&ud&&ud.anisoM){
        _atomM.copy(ud.anisoM);
        _atomM.setPosition(h.position.x,h.position.y,h.position.z);
        inst.setMatrixAt(slot,_atomM);
    }else if(r>0){
        _atomM.makeScale(r,r,r);
        _atomM.setPosition(h.position.x,h.position.y,h.position.z);
        inst.setMatrixAt(slot,_atomM);
    }else inst.setMatrixAt(slot,_atomZero);
}
function atomClassEnsure(layer,key,size,geo,geoKey,parent){
    var c=layer.classes[key];
    if(!c){
        c={key:key,insts:[null,null,null],cap:0,geoKey:'',slots:[],shown:[0,0,0]};
        layer.classes[key]=c;
    }
    // Grow, never shrink: an edit that only rewrites the same structure reuses
    // the buffers it already uploaded. A change of sphere tessellation (the atom
    // count crossed a band) is the one case that has to rebuild the geometry
    // reference.
    if(c.cap<size||c.geoKey!==geoKey){
        for(var s=0;s<3;s++){
            var old=c.insts[s];
            if(old){if(old.parent)old.parent.remove(old);old.dispose()}
            var inst=new THREE.InstancedMesh(geo,atomInstanceMaterial(key,s),size);
            // The instances span the whole structure while the shared unit
            // sphere is one unit wide, so the mesh's own bounds say nothing
            // about where the atoms are; one draw call is cheap enough to always
            // draw, exactly as the bond pools do.
            inst.frustumCulled=false;
            inst.count=0;
            c.insts[s]=inst;
        }
        c.cap=size;c.geoKey=geoKey;
    }
    for(var k=0;k<3;k++){
        var now=c.insts[k];
        if(now.parent!==parent)parent.add(now);
    }
    return c;
}
// Build (or refresh) a layer from the handles and the atom list they describe.
// isDiff is the comparison's right half, where a marked atom is drawn in the
// comparison's orange instead of a glow over its element colour.
function sphereLayerBuild(layer,handles,parent,segs,isDiff){
    layer.parent=parent;
    layer.count=handles.length;
    var used={},order=[],i;
    for(i=0;i<handles.length;i++){
        var h=handles[i],ud=h.userData;
        var app=atomAppearanceOf(ud);
        ud.sphereAppearance=app;
        ud.sphereLayer=layer;
        var list=used[app];
        if(!list){list=[];used[app]=list;order.push(app)}
        ud.sphereSlot=list.length;
        list.push(h);
    }
    // An appearance no atom carries any more (a ghost converted back, a
    // disordered site removed) takes its three meshes with it.
    Object.keys(layer.classes).forEach(function(k){
        if(used[k])return;
        var c=layer.classes[k];
        for(var s=0;s<3;s++){
            var inst=c.insts[s];
            if(inst){if(inst.parent)inst.parent.remove(inst);inst.dispose()}
        }
        delete layer.classes[k];
    });
    if(handles.length===0)return;
    var geo=unitSphereGeo(segs);
    var geoKey=segs[0]+'x'+segs[1];
    order.forEach(function(k){
        var members=used[k];
        var c=atomClassEnsure(layer,k,members.length,geo,geoKey,parent);
        c.slots=members;
        var occupancy=k.charAt(0)==='o';
        for(var s=0;s<3;s++){
            c.insts[s].count=members.length;
            c.shown[s]=0;
        }
        for(var m=0;m<members.length;m++){
            var mm=members[m],mu=mm.userData;
            var state=mu.hlState||0;
            // A ghost and an occupancy material already carry their colour; a
            // plain atom takes the element colour from the instance buffer. The
            // difference half paints every marked atom orange.
            var col=(occupancy||mu.ghost)?ATOM_WHITE:atomColorOf(mu.baseColor||'#FF1493');
            for(s=0;s<3;s++)c.insts[s].setColorAt(m,(isDiff&&s===2)?ATOM_ORANGE:col);
            for(s=0;s<3;s++){
                var show=(s===state)&&mu.drawn!==false;
                writeAtomMatrix(c.insts[s],m,mm,show?mu.radius:0);
                if(show)c.shown[s]++;
            }
        }
        for(s=0;s<3;s++){
            var ins=c.insts[s];
            ins.instanceMatrix.needsUpdate=true;
            if(ins.instanceColor)ins.instanceColor.needsUpdate=true;
            ins.visible=c.shown[s]>0;
        }
    });
}
// One atom changes state: hide its slot in the old mesh, show it in the new one.
// Both writes are constant time, and the object count on screen is unchanged.
function sphereLayerWriteState(layer,h,state,show){
    var ud=h.userData;
    var c=layer&&layer.classes[ud.sphereAppearance];
    if(!c)return;
    var inst=c.insts[state];
    if(!inst)return;
    var vis=!!show&&ud.drawn!==false;
    writeAtomMatrix(inst,ud.sphereSlot,h,vis?ud.radius:0);
    c.shown[state]+=vis?1:-1;
    if(c.shown[state]<0)c.shown[state]=0;
    inst.instanceMatrix.needsUpdate=true;
    inst.visible=c.shown[state]>0;
    needsRender=true;
}
// After a drag, a vibration frame or a re-centring, every handle has moved. Only
// the state a slot is currently shown in carries its position, so one matrix
// write per atom is enough.
function sphereLayerSyncPositions(layer){
    if(!layer||layer.count===0)return;
    Object.keys(layer.classes).forEach(function(k){
        var c=layer.classes[k];
        for(var m=0;m<c.slots.length;m++){
            var h=c.slots[m],ud=h.userData;
            writeAtomMatrix(c.insts[ud.hlState||0],m,h,ud.drawn!==false?ud.radius:0);
        }
        // Only a state that has atoms on screen was written; the other two
        // buffers are still all-zero and do not need re-uploading.
        for(var s=0;s<3;s++)if(c.shown[s]>0)c.insts[s].instanceMatrix.needsUpdate=true;
    });
}
function sphereLayerClear(layer){
    if(!layer)return;
    Object.keys(layer.classes).forEach(function(k){
        var c=layer.classes[k];
        for(var s=0;s<3;s++){
            var inst=c.insts[s];
            if(inst){if(inst.parent)inst.parent.remove(inst);inst.dispose()}
        }
    });
    layer.classes={};layer.count=0;layer.parent=null;
}

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
        // The current molecule draws no eggs: no panel, no chip, choices
        // forgotten.
        panel.classList.remove('show');
        var chip0=document.getElementById('egg-reopen');
        if(chip0)chip0.classList.remove('show');
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
    // A minimized list stays folded while the structure is rebuilt underneath
    // it; only an explicit open (the toolbar button or the chip) unfolds it.
    panel.classList.toggle('show',!panelMinState.egg);
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
    // The x folds the list into its chip; the chip (or the toolbar's Eggs
    // button) unfolds it again.
    if(close)close.addEventListener('click',function(){
        panelMinState.egg=true;
        panel.classList.remove('show');
        syncPanelChips();
        layoutPanels();
    });
    syncPanelChips();
    layoutPanels();
}
function openEggPanel(){
    panelMinState.egg=false;
    rebuildEggPanel();
    layoutPanels();
}
var eggBtnWired=false;
function wireEggPanelButton(){
    if(eggBtnWired)return;
    var btn=document.getElementById('egg-btn');
    if(!btn)return;
    eggBtnWired=true;
    // 工具栏按钮永远"打开"面板 (幂等): 关闭只由面板内的 × 负责, 这样任何
    // 时候都有地方可以再打开。重复点击只是刷新面板内容。
    btn.addEventListener('click',function(){openEggPanel()});
}

function rebuildScene(keepCenter){
    // Every structural change passes through here; the interaction cache is a
    // reading of that structure, so it is retired with it. View-only updates
    // that do not rebuild the scene (a checkbox toggle) keep the cache.
    intStructRev++;
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
    cellGridLines=null;
    cellEdgeLabels=[];
    rotAxisLineMesh=null;
    neighLines=null;
    // The plane is a child of the group the clear above emptied; hklRefresh at
    // the end of the build puts the current one back.
    hklPlaneMesh=null;
    hklPlaneBorder=null;
    // The polyhedra are children of that same group; buildPolyhedra below
    // recreates them from the geometry this build produces.
    polyFacesMesh=null;
    polyEdgesMesh=null;
    // Recomputing the centroid is the default (the structure changed, the
    // centre follows it). A supercell bounds change passes keepCenter: the
    // atoms the old bounds already showed keep their exact screen positions
    // and the new material appears around them, instead of the whole picture
    // jumping to the enlarged structure's new centroid.
    if(!keepCenter){
        CX=0;CY=0;CZ=0;
        MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
        if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    }
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
    // The stick style's joints are the thickness of the sticks meeting at each
    // atom, read off the bond list once for the whole build. The other two
    // styles do not consult it, so the pass is skipped rather than paid for.
    var stickR=(displayStyle==='stick')?stickJointRadii():null;
    MD.atoms.forEach(function(a,i){
        a.index=i;
        var gh=isGhostAtom(a);
        var r=drawnRadius(a,stickR?stickR[i]:0);
        var drawn=r>0.0005;
        // 标量配色是视图状态: 在这里按原子现算, a.color 保持元素颜色不变。
        // A ghost has no element colour and no scalar to show: it keeps the one
        // marker colour, so a counterpoise ghost never looks like the atom it
        // stands for.
        var col=gh?GHOST_COLOR:(betaColorOf(a)||a.color);
        var occ=a.occupancy!=null?a.occupancy:1;
        // One record per atom, not one Mesh per atom: radius, colour and the
        // appearance key are what the instanced layer needs, and every other
        // pass still reads "position", "userData" and "visible" exactly as
        // before.
        var mesh=atomHandle();
        // Polyhedra draws the centre atom as a marker under its hull, not as
        // the ball an ellipsoid stands in for: the ellipsoid is held back
        // there, and the switch's chosen state returns with the ball style.
        var aniso=(ortepEnabled&&displayStyle!=='polyhedra'&&!gh&&a.uAniso)?anisoEllipsoid(a.uAniso):null;
        mesh.position.set(a.x-CX,a.y-CY,a.z-CZ);
        mesh.userData={element:a.element,index:i,baseColor:col,occ:occ,ghost:gh,
                       radius:drawn?(aniso?aniso.r:r):0,drawn:drawn,hlState:0,
                       anisoM:aniso?aniso.m:null};
        mesh.visible=drawn;
        atomMeshes.push(mesh);
    });
    sphereLayerBuild(mainSphereLayer,atomMeshes,moleculeGroup,segs,false);
    rebuildBondPieces();
    buildPolyhedra();
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
    // The contact lines are children of the molecule group the clear above
    // removed: put them back against the geometry that was just built, so the
    // panel and the picture never disagree after an edit, an undo or a switch.
    if(neighActive)neighRefresh();
    // The interaction list is a reading of the same geometry: an edit, an undo
    // or a switch re-detects it against the atoms that are on screen now.
    if(intActive)intRefresh();
    // The lattice-plane cut follows the cell it was drawn in: an edit, an undo
    // or a supercell change rebuilds it against the new box.
    hklRefresh();
    // The powder pattern is a reading of the cell's contents, so an edit that
    // moved an atom or changed a site's element must reach it too. The panel
    // stays closed (and costs nothing) until it is opened.
    xrdRefresh();
    // The cavity surface belongs to the molecule group the clear at the top of
    // this function emptied; put the selected void back against the geometry
    // that was just rebuilt.
    cavityRefresh();
    // The ellipsoid switch only exists for a structure that states ADPs; a
    // molecule read from a coordinate file (or a crystal whose file has no
    // aniso table) leaves it disabled rather than offering a view that would
    // draw the same spheres as before.
    syncOrtepControl();
}
function updateScenePositions(keepCenter){
    if(!keepCenter){
        CX=0;CY=0;CZ=0;
        MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
        if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    }
    atomMeshes.forEach(function(m,i){var a=MD.atoms[i];if(a)m.position.set(a.x-CX,a.y-CY,a.z-CZ)});
    sphereLayerSyncPositions(mainSphereLayer);
    rebuildBondPieces();
    buildPolyhedra();
    if(CRY)buildCellWireframe();
    hklRefresh();
    highlightSelected();
    hideVmdSurfaces();
    // A drag moves the atoms under a surface that describes the geometry
    // before it; the cavity mesh is dropped for the duration and rebuilt by
    // the next rebuildScene (a commit invalidates the report outright).
    cavityHideSurface();
    needsRender=true;
}

function updateAtomMeshPositions(){
    // Light per-frame update used during Move Atoms drags: refresh atom
    // sphere positions + selection highlights only (no bond mesh rebuild).
    atomMeshes.forEach(function(m,i){var a=MD.atoms[i];if(a)m.position.set(a.x-CX,a.y-CY,a.z-CZ)});
    sphereLayerSyncPositions(mainSphereLayer);
    // The hulls follow their atoms for the frame; their shapes are recomputed
    // by the commit that closes the drag (updateRotationCenter sees the flag).
    polySyncPositions();
    if(polyBody)polyShapeDirty=true;
    highlightSelected();
    hideVmdSurfaces();
    cavityHideSurface();
}

// 网格位置都是"相对旋转中心"存放的 (a - C), 所以中心一变, 所有网格只需整体
// 平移 -Δ 就能保持与原子坐标一致 —— 不用重建任何几何。
function shiftAllMeshesForCenter(dx,dy,dz){
    if(!dx&&!dy&&!dz)return;
    var i;
    for(i=0;i<atomMeshes.length;i++){
        atomMeshes[i].position.x-=dx;atomMeshes[i].position.y-=dy;atomMeshes[i].position.z-=dz;
    }
    sphereLayerSyncPositions(mainSphereLayer);
    // A bond piece's transform is computed from the atoms *minus the centre*,
    // so re-placing them all against the new centre is the same shift — and it
    // is one pass over the instances, with no geometry involved.
    placeBondPieces();
    // The hull vertices are written relative to the same centre; a new centre
    // moves them by the same -Δ, again without rebuilding a hull.
    polySyncPositions();
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
    // 接触线同样相对中心存放, 中心一变就按新的几何重画。
    if(neighActive)neighRefresh();
    // The interaction list is a reading of the same geometry: an edit, an undo
    // or a switch re-detects it against the atoms that are on screen now.
    if(intActive)intRefresh();
    // 中心变了, 表面片也要按新中心重摆 (签名里含中心, 没变时是空操作)。
    refreshVmdSurfaces();
    // A drag moved the atoms relative to one another: the hulls followed them
    // for the duration, and this commit recomputes which points form each hull
    // from the coordinates the atoms now have. A centre-only shift was already
    // handled by shiftAllMeshesForCenter's position sync, so it does not pay
    // for a rebuild.
    if(polyShapeDirty)buildPolyhedra();
    needsRender=true;
}

// Rebuild every piece from the current bond list. A piece is pure
// description, so building the list is cheap; the instanced mesh that draws
// them is only ever recreated when the count outgrows the capacity it has.
function rebuildBondPieces(){
    // The aromatic dash sides belong to the structure and are planned before
    // any piece is declared; every other piece is described by its bond alone.
    planAromaticDashes();
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
    // The moved region was re-derived, so the aromatic systems it touches have
    // to be re-planned before the region's pieces are built (a piece reads its
    // side and plane from the plan when it is placed). The kept pieces need no
    // touching: they read the same, now current, plan.
    planAromaticDashes();
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
        //
        // An aromatic dash instead offsets inside the plane of its aromatic
        // system, which the plan knows (see the section below): crossing the
        // bond with that plane's normal keeps the offset perpendicular to the
        // bond *and* inside the plane, which is what makes a ring show all its
        // dashes on one side. The plane and the side are read from the plan
        // here rather than stored in the piece, so a re-planned system — a
        // fragment rotation turns every ring plane with it, an MM run relaxes
        // it — reaches every piece that is placed from then on.
        var nx=0,ny=0,nz=0,side=p.perpSign;
        if(p.isDashed&&aromPlan.sign&&aromPlan.sign[p.bondIdx]){
            nx=aromPlan.nx[p.bondIdx];ny=aromPlan.ny[p.bondIdx];nz=aromPlan.nz[p.bondIdx];
            side=aromPlan.sign[p.bondIdx];
        }
        var qx=0,qy=0,qz=0,ql=0;
        if(nx||ny||nz){
            qx=ny*dz-nz*dy;qy=nz*dx-nx*dz;qz=nx*dy-ny*dx;
            ql=Math.sqrt(qx*qx+qy*qy+qz*qz);
            if(ql>0.0001){qx/=ql;qy/=ql;qz/=ql}
        }
        if(ql<=0.0001){
            // No plane (or a bond that runs along the normal itself, where the
            // projection would be zero): the world-up cross, as before.
            var upx=0,upy=1,upz=0;
            if(Math.abs(dy)>=0.99){upx=1;upy=0;upz=0}
            qx=dy*upz-dz*upy;qy=dz*upx-dx*upz;qz=dx*upy-dy*upx;
            ql=Math.sqrt(qx*qx+qy*qy+qz*qz);
            if(ql>0.0001){qx/=ql;qy/=ql;qz/=ql}
        }
        var off=p.perpMag*side;
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
// shifted onto one side by the same kind of offset a double bond uses. Which
// side, and in which plane, comes from the plan below — read when the piece is
// placed, so a plan that is recomputed (a fragment rotation turns every ring
// plane with it) reaches every piece without touching one.
function addDashPieces(r,col,off){
    var dashes=6,seg=1/dashes,gap=seg*0.35,dashLen=seg-gap;
    for(var k=0;k<dashes;k++){
        var t0=k*seg+gap*0.5;
        var t1=t0+dashLen;
        if(t1>1)t1=1;
        addBondPiece(r,col,true,t0,t1,1,off,0);
    }
}

// ===== Aromatic dashes: which side they sit on =====
//
// An aromatic bond is one solid line plus a six-dash line offset to one side.
// That side used to fall out of the bond alone — the offset was the bond
// direction crossed with the world up axis — so two aromatic bonds meeting at
// an atom could put their dashes on opposite sides. Around a ring that reads
// badly: seen from the face of the ring, half the bonds carry their dashes
// inside and half outside, so only part of the ring looks aromatic, which is
// exactly the cue a reader looks for.
//
// The side is now planned once per bond list, from the aromatic bonds
// themselves, and every bond of one connected aromatic system offsets in that
// system's own plane. The plan is a walk over the system, and a bond's sign
// says which of its two sides the dash goes on:
//
//   * at a two-bond atom the band simply passes through — the second bond
//     continues on the same side of the path (in sign space: one bond arrives,
//     one leaves), which is what makes a ring or a chain show one continuous
//     band;
//   * a three-bond atom is a branch, where the band must pick a side for each
//     arm: each arm takes the side that points most nearly the way the
//     arriving dash points, so a fused ring keeps the arm that stays in its own
//     ring, and the shared bond carries the band across;
//   * crossing a bond that no ring runs through (a biaryl link) starts a new
//     band on the far side — mirrored — because the arms there are a ring of
//     their own and must not inherit the side the band arrived on.
//
// A fused system can still be left with one atom whose two dashes face
// opposite ways (a shared bond belongs to two rings and the walk only knows the
// one it came through); a bounded local repair settles those, flipping only
// bonds whose flip leaves fewer contradicted atoms behind.
//
// Which of the two sides a system settled on is arbitrary — flipping every
// sign of a system keeps it self-consistent and just moves all its dashes to
// the other side — so the plan is flipped when the majority of its bonds would
// then point away from the system's own centre. A ring therefore keeps its
// dashes **inside**, the way it is drawn on paper.
//
// Nothing here reads the atom order a file happens to use: a bond stored the
// other way round has its direction — and so its "left" — flipped, and sigma
// (which endpoint is atom1) accounts for that when the sign is written.
var aromPlan={sign:null,nx:null,ny:null,nz:null};
function isAromaticOrder(o){return o>=1.25&&o<1.75}
// The far endpoint of a bond as it is drawn, with the periodic shift applied
// exactly as placePiece() applies it, so a crystal's wrapped bond is planned
// along the line it is actually drawn on.
function bondEnd(b,out){
    var a2=MD.atoms[b.atom2];
    if(!a2)return false;
    var ex=a2.x,ey=a2.y,ez=a2.z;
    if(b.shift&&CRY){
        var lv=CRY.latticeVectors;
        ex+=b.shift[0]*lv[0][0]+b.shift[1]*lv[1][0]+b.shift[2]*lv[2][0];
        ey+=b.shift[0]*lv[0][1]+b.shift[1]*lv[1][1]+b.shift[2]*lv[2][1];
        ez+=b.shift[0]*lv[0][2]+b.shift[1]*lv[1][2]+b.shift[2]*lv[2][2];
    }
    out[0]=ex;out[1]=ey;out[2]=ez;
    return true;
}
// Unit direction of one bond.
var _bondScratch=[0,0,0];
function bondUnitDir(b,out){
    var a1=MD.atoms[b.atom1];
    if(!a1||!bondEnd(b,_bondScratch))return false;
    var dx=_bondScratch[0]-a1.x,dy=_bondScratch[1]-a1.y,dz=_bondScratch[2]-a1.z;
    var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(l<0.0001)return false;
    out[0]=dx/l;out[1]=dy/l;out[2]=dz/l;
    return true;
}
// One connected aromatic system: find its plane, then walk it and write a side
// (±1) and that plane into every bond it contains.
function planAromaticComponent(comp,adj,dx,dy,dz,sign,nx,ny,nz,mark,stamp,isBridge){
    var bonds=MD.bonds,i,k,p,q;
    // The bonds of the system span its plane, so the cross product of two bonds
    // meeting at an atom is the plane normal (up to sign). Take the pair that
    // is least collinear as the reference — one distorted bond should not tilt
    // a plane it is drawn in — and average every pair aligned to it.
    var atoms=[];
    for(i=0;i<comp.length;i++){
        var cb=bonds[comp[i]];
        if(mark[cb.atom1]!==stamp){mark[cb.atom1]=stamp;atoms.push(cb.atom1)}
        if(mark[cb.atom2]!==stamp){mark[cb.atom2]=stamp;atoms.push(cb.atom2)}
    }
    var rx=0,ry=0,rz=0,refLen=0;
    for(i=0;i<atoms.length;i++){
        var al=adj[atoms[i]];
        if(!al||al.length<2)continue;
        for(p=0;p<al.length-1;p++)for(q=p+1;q<al.length;q++){
            var A=al[p],B=al[q];
            var cx=dy[A]*dz[B]-dz[A]*dy[B],
                cy=dz[A]*dx[B]-dx[A]*dz[B],
                cz=dx[A]*dy[B]-dy[A]*dx[B];
            var cl=cx*cx+cy*cy+cz*cz;
            if(cl>refLen){refLen=cl;rx=cx;ry=cy;rz=cz}
        }
    }
    if(refLen>1e-12){
        var rl=Math.sqrt(refLen);rx/=rl;ry/=rl;rz/=rl;
        var sx=0,sy=0,sz=0;
        for(i=0;i<atoms.length;i++){
            var al2=adj[atoms[i]];
            if(!al2||al2.length<2)continue;
            for(p=0;p<al2.length-1;p++)for(q=p+1;q<al2.length;q++){
                var A2=al2[p],B2=al2[q];
                var cx2=dy[A2]*dz[B2]-dz[A2]*dy[B2],
                    cy2=dz[A2]*dx[B2]-dx[A2]*dz[B2],
                    cz2=dx[A2]*dy[B2]-dy[A2]*dx[B2];
                if(cx2*rx+cy2*ry+cz2*rz<0){cx2=-cx2;cy2=-cy2;cz2=-cz2}
                sx+=cx2;sy+=cy2;sz+=cz2;
            }
        }
        var sl=Math.sqrt(sx*sx+sy*sy+sz*sz);
        if(sl>1e-9){rx=sx/sl;ry=sy/sl;rz=sz/sl}
    }else{
        // Every bond of the system is collinear (a straight chain), so the
        // bonds name no plane: any perpendicular to the chain is as good, and
        // one fixed axis keeps the dashes on one side of it.
        var ax=dx[comp[0]],ay=dy[comp[0]],az=dz[comp[0]];
        var ux=0,uy=1,uz=0;
        if(Math.abs(ay)>0.9){ux=1;uy=0;uz=0}
        rx=ay*uz-az*uy;ry=az*ux-ax*uz;rz=ax*uy-ay*ux;
        var vl=Math.sqrt(rx*rx+ry*ry+rz*rz);
        if(vl<1e-9){rx=1;ry=0;rz=0}else{rx/=vl;ry/=vl;rz/=vl}
    }
    // Walk: sign = +1 puts the dash on the left of atom1 -> atom2, -1 on the
    // other side. At an atom, a bond that arrives there (its dash already on
    // the left of the travel direction) sends every other bond of that atom out
    // on the left of its own direction; at a two-bond atom, a bond that leaves
    // takes the arriving one's place, so the band stays on one side of the path
    // through the whole walk.
    //
    // A three-bond atom is a branch, and there the band cannot simply continue
    // (a bond has two sides and the arm must pick one), so each arm takes the
    // side that keeps the band's own direction: the dash that points most
    // nearly the way the arriving dash points. That is what keeps a fused ring
    // whole — the perimeter arm stays on the side the band was already on,
    // while the shared bond takes whichever of its two sides that is, and the
    // ring on the far side then starts from there.
    //
    // One case is the opposite: when the band arrives over a bond that no ring
    // runs through (a biaryl link) and the arms are ring bonds, the arms are a
    // ring of their own and must not inherit the band's side — that would draw
    // the second ring's dashes on the outside. Those arms are mirrored instead,
    // which puts them on the inside of their own ring. A bridge whose arms are
    // bridges too is a branch of a chain, where there is no inside to prefer,
    // so it keeps the band's direction like any other branch.
    var seed=comp[0],stack=[seed];
    sign[seed]=1;
    while(stack.length){
        var cur=stack.pop();
        var b2=bonds[cur];
        for(var e=0;e<2;e++){
            var X=e===0?b2.atom1:b2.atom2;
            var sg=e===0?1:-1;                    // +1 when X is the bond's atom1
            var incoming=sign[cur]===-sg;
            var alx=adj[X];
            if(!alx)continue;
            var deg=alx.length;
            // The arriving dash's direction, for the branch rule below.
            var oex=0,oey=0,oez=0;
            if(deg>2&&incoming){
                oex=sign[cur]*(ry*dz[cur]-rz*dy[cur]);
                oey=sign[cur]*(rz*dx[cur]-rx*dz[cur]);
                oez=sign[cur]*(rx*dy[cur]-ry*dx[cur]);
            }
            var acrossBridge=deg>2&&incoming&&isBridge&&isBridge[cur];
            for(k=0;k<alx.length;k++){
                var f=alx[k];
                if(f===cur||sign[f])continue;
                var sgF=bonds[f].atom1===X?1:-1;
                if(deg>2&&incoming){
                    var cx3=ry*dz[f]-rz*dy[f],cy3=rz*dx[f]-rx*dz[f],cz3=rx*dy[f]-ry*dx[f];
                    var ahead=(oex*cx3+oey*cy3+oez*cz3)>=0;
                    if(acrossBridge&&!isBridge[f])ahead=!ahead;
                    sign[f]=ahead?sgF:-sgF;
                }else{
                    sign[f]=(incoming||deg>2)?sgF:-sgF;
                }
                stack.push(f);
            }
        }
    }
    // A bond shared by two rings belongs to both, and the walk only knows the
    // ring it came through, so a fused system can be left with one atom whose
    // two dashes face opposite ways. In sign space the rule at a two-bond atom
    // is exactly "one bond in, one bond out" (a bond that leaves offsets to the
    // left of its own direction, a bond that arrives offsets to the left of the
    // direction it came from, which is the same side of the path), and flipping
    // a bond flips only its own two ends — so a violated atom is settled
    // locally, by taking whichever of its two bonds leaves fewer violated atoms
    // behind. A flip is only taken when it strictly reduces that count, so a
    // ring that is already consistent is never touched.
    var stateOf=function(bi,X){return sign[bi]===(bonds[bi].atom1===X?1:-1)?1:-1};
    var badAt=function(X){
        var a=adj[X];
        if(!a||a.length!==2)return 0;
        var be=bonds[a[0]];
        if(be.atom1===be.atom2)return 0;
        return stateOf(a[0],X)===stateOf(a[1],X)?1:0;
    };
    for(var pass=0;pass<3;pass++){
        var flipped=false;
        for(i=0;i<atoms.length;i++){
            var vx=atoms[i];
            if(!badAt(vx))continue;
            var pair=adj[vx];
            var pick=-1,pickGain=0;
            for(k=0;k<2;k++){
                var cand=pair[k],cb2=bonds[cand];
                var was=badAt(cb2.atom1)+badAt(cb2.atom2);
                sign[cand]=-sign[cand];
                var now=badAt(cb2.atom1)+badAt(cb2.atom2);
                sign[cand]=-sign[cand];
                if(was-now>pickGain){pickGain=was-now;pick=cand}
            }
            if(pick>=0){sign[pick]=-sign[pick];flipped=true}
        }
        if(!flipped)break;
    }
    // Which side the walk settled on is arbitrary until it is compared with
    // where the system sits: a ring's bonds outvote a chain's, whose own vote
    // means nothing (a chain's centre lies on the chain itself).
    var centreX=0,centreY=0,centreZ=0;
    for(i=0;i<atoms.length;i++){
        var ca=MD.atoms[atoms[i]];
        centreX+=ca.x;centreY+=ca.y;centreZ+=ca.z;
    }
    centreX/=atoms.length;centreY/=atoms.length;centreZ/=atoms.length;
    var vote=0,e3=[0,0,0];
    for(i=0;i<comp.length;i++){
        var bi=comp[i],bb=bonds[bi];
        var ba1=MD.atoms[bb.atom1];
        if(!sign[bi])sign[bi]=1;
        if(!ba1||!bondEnd(bb,e3))continue;
        vote+=sign[bi]*((ry*dz[bi]-rz*dy[bi])*((ba1.x+e3[0])*0.5-centreX)+
                        (rz*dx[bi]-rx*dz[bi])*((ba1.y+e3[1])*0.5-centreY)+
                        (rx*dy[bi]-ry*dx[bi])*((ba1.z+e3[2])*0.5-centreZ));
    }
    if(vote>0){for(i=0;i<comp.length;i++)sign[comp[i]]=-sign[comp[i]]}
    for(i=0;i<comp.length;i++){
        var si=comp[i];
        nx[si]=rx;ny[si]=ry;nz[si]=rz;
    }
}
// Which aromatic bonds no ring runs through (the bridges of the aromatic
// subgraph). Only asked for when some atom has three aromatic bonds, which is
// the only place the answer is used, so a structure of separate rings and
// chains pays nothing for it. One iterative DFS, disc/low as usual.
function markAromaticBridges(adj,disc,low,isBridge,na){
    var time=0,a,v,bi,b,w,i;
    for(a=0;a<na;a++){
        if(!adj[a]||disc[a])continue;
        var st=[a],stEdge=[-1],stIdx=[0];
        disc[a]=low[a]=++time;
        while(st.length){
            v=st[st.length-1];
            i=stIdx[stIdx.length-1];
            var list=adj[v];
            if(i<list.length){
                stIdx[stIdx.length-1]=i+1;
                bi=list[i];b=MD.bonds[bi];
                w=b.atom1===v?b.atom2:b.atom1;
                if(w===v)continue;
                if(!disc[w]){
                    disc[w]=low[w]=++time;
                    st.push(w);stEdge.push(bi);stIdx.push(0);
                }else if(bi!==stEdge[stEdge.length-1]&&disc[w]<low[v]){
                    low[v]=disc[w];
                }
            }else{
                st.pop();stIdx.pop();
                var pe=stEdge.pop();
                if(st.length){
                    var parent=st[st.length-1];
                    if(low[v]<low[parent])low[parent]=low[v];
                    if(low[v]>disc[parent]&&pe>=0)isBridge[pe]=1;
                }
            }
        }
    }
}
// Plan every aromatic bond of the current bond list. Called wherever the bond
// list is (re)built — the scene build, a frame or molecule switch, an edit, an
// Undo, an Import, an MM commit and a Move Atoms release — and never on a path
// that only re-places existing pieces, because a bond's side is a property of
// the structure, not of the instant it is drawn. Cost is one pass over the
// bonds, one walk per aromatic system, and a bridge search that is only taken
// when some atom has three aromatic bonds (the one place its answer is used).
function planAromaticDashes(){
    var bonds=MD.bonds,nb=bonds.length,i,na=MD.atoms.length;
    aromPlan={sign:null,nx:null,ny:null,nz:null};
    if(nb===0)return;
    // A structure without a single aromatic bond pays one scan and nothing
    // else — no arrays, no adjacency, no walk — which is what the overwhelming
    // majority of files (any protein, any alkane) get from this feature.
    var any=false;
    for(i=0;i<nb;i++){
        if(isAromaticOrder(bonds[i].order)){any=true;break}
    }
    if(!any)return;
    var sign=new Int8Array(nb),nx=new Float32Array(nb),ny=new Float32Array(nb),nz=new Float32Array(nb);
    aromPlan={sign:sign,nx:nx,ny:ny,nz:nz};
    var dx=new Float32Array(nb),dy=new Float32Array(nb),dz=new Float32Array(nb);
    // Atom -> its aromatic bonds, as a plain indexed array: this is walked once
    // per atom per component and once more inside the repair, and an array of
    // arrays stays on the fast path where an object keyed by atom index does
    // not (a structure of tens of thousands of atoms turns that object into a
    // dictionary). Both this and the mark array below are only created once an
    // aromatic bond has actually been seen, so a structure without any pays
    // nothing for the plan.
    var adj=null,mark=null,list=[],d3=[0,0,0],maxDeg=0;
    for(i=0;i<nb;i++){
        var b=bonds[i];
        if(!isAromaticOrder(b.order))continue;
        if(!bondUnitDir(b,d3))continue;
        if(!adj)adj=new Array(na);
        dx[i]=d3[0];dy[i]=d3[1];dz[i]=d3[2];
        var a1=b.atom1,a2=b.atom2;
        if(a1!==a2){
            var l1=adj[a1];if(!l1)l1=adj[a1]=[];
            l1.push(i);
            if(l1.length>maxDeg)maxDeg=l1.length;
        }
        var l2=adj[a2];if(!l2)l2=adj[a2]=[];
        l2.push(i);
        if(l2.length>maxDeg)maxDeg=l2.length;
        list.push(i);
    }
    if(list.length===0)return;
    var isBridge=null;
    if(maxDeg>2){
        isBridge=new Uint8Array(nb);
        markAromaticBridges(adj,new Int32Array(na),new Int32Array(na),isBridge,na);
    }
    var seen=new Uint8Array(nb),stack=[];
    // One mark array for the whole pass, stamped per component, so a structure
    // made of many small aromatic fragments does not allocate one set per ring.
    mark=new Int32Array(na);
    var stamp=0;
    for(var s=0;s<list.length;s++){
        var seed=list[s];
        if(seen[seed])continue;
        var comp=[];
        stack.length=0;stack.push(seed);seen[seed]=1;
        while(stack.length){
            var cur=stack.pop();
            comp.push(cur);
            var cb=bonds[cur];
            for(var e=0;e<2;e++){
                var al=adj[e===0?cb.atom1:cb.atom2];
                if(!al)continue;
                for(var k=0;k<al.length;k++){
                    var nb2=al[k];
                    if(!seen[nb2]){seen[nb2]=1;stack.push(nb2)}
                }
            }
        }
        stamp++;
        planAromaticComponent(comp,adj,dx,dy,dz,sign,nx,ny,nz,mark,stamp,isBridge);
    }
}

// Declare the pieces one bond is drawn with. The order thresholds are the same
// the viewer has always used: one line below 1.25, an aromatic line with its
// dashed partner below 1.75, two lines below 2.5, three below 3.5.
function createBond(b){
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2)return;
    // A ghost has no valence and no bonds: every path that makes or keeps a
    // bond drops the ones that touch one (bond perception, loading, editing,
    // pasting). This is the last door - a bond list restored from an older
    // snapshot, or one built before an atom became a ghost, must not draw a
    // stick into a point either.
    if(isGhostAtom(a1)||isGhostAtom(a2))return;
    // A space-filling model has no sticks: its spheres are van der Waals
    // spheres, which already touch, so a cylinder between two atom centres would
    // only show where it pokes out of them. This is the one door a piece is
    // declared through - the scene build, a moved region's rebuild and a
    // comparison's halves all come here - so no path can leave a stick behind.
    if(displayStyle==='spacefill')return;
    // A VMD script sets the bond radius per bond (the representation that
    // covers both endpoints); 0 means the script does not draw it at all.
    var br=bondPieceRadius(b);
    if(br<=0.0005)return;
    var ord=b.order||1;
    // Multiple-bond line offsets keep their proportions to the bond radius.
    var offK=br/BOND_R;
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
// Which panels are folded into their chip. The x sets a flag, the chip clears
// it, and syncPanelChips() is the one place that turns the flags (plus the
// mode and the loaded data) into display, so a panel and its chip can never
// both be visible and neither can be lost.
var panelMinState={crystal:false,egg:false,view:false,bake:false,select:false,neigh:false,rotate:false,xrd:false,cav:false};
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
    // A rotation turns each ring's plane with the fragment, and an aromatic
    // dash is offset inside that plane: re-plan before re-placing, or the
    // dashes would tilt out of their ring for the length of the drag. The plan
    // is linear in the aromatic bonds and every piece reads it when it is
    // placed, so this costs one pass and touches no piece.
    planAromaticDashes();
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
    // One door for every way the bounds change: the six inputs, the N×N×N
    // presets, Reset and Align Origin all go through applyBoundary, so the size
    // guard, the steady centre and the selection remap behave the same no
    // matter which one moved.
    function boundsEqual(a,b){
        return a.aMin===b.aMin&&a.aMax===b.aMax&&a.bMin===b.bMin&&a.bMax===b.bMax&&a.cMin===b.cMin&&a.cMax===b.cMax;
    }
    function applyBoundary(nb){
        if(!CRY)return false;
        if(boundsEqual(nb,boundary)){syncBndInputs();return true}
        // Scale guard: refuse, say so, and leave the structure (and the six
        // inputs, which keep showing the bounds still in force) as they were.
        var est=CRY.baseAtoms.length*crystalCellEstimate(nb);
        if(est>MAX_CRYSTAL_ATOMS){
            syncBndInputs();
            showModal('<h3>Supercell Too Large</h3><div class="current-val">These bounds would build about '+est+' atoms; the limit is '+MAX_CRYSTAL_ATOMS+'.<br>Narrow the boundary range.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
            document.getElementById('m-cancel').addEventListener('click',hideModal);
            return false;
        }
        // The atoms are about to be re-indexed: record each selected atom by
        // its lattice identity (base site + image cell), then re-find that same
        // atom after the rebuild. An image the new bounds no longer contain
        // drops out of the selection — its base site stays selected wherever it
        // is still on screen.
        var selKeys=selectedAtoms.map(function(i){
            var a=MD.atoms[i];
            return (a&&a.baseIdx!=null&&a.cellI!=null)?(a.baseIdx+'_'+a.cellI+'_'+a.cellJ+'_'+a.cellK):null;
        });
        boundary.aMin=nb.aMin;boundary.aMax=nb.aMax;
        boundary.bMin=nb.bMin;boundary.bMax=nb.bMax;
        boundary.cMin=nb.cMin;boundary.cMax=nb.cMax;
        syncBndInputs();
        rebuildCrystal();
        var newIdx={};
        MD.atoms.forEach(function(a,i){
            if(a.baseIdx!=null&&a.cellI!=null){
                var k=a.baseIdx+'_'+a.cellI+'_'+a.cellJ+'_'+a.cellK;
                if(newIdx[k]===undefined)newIdx[k]=i;
            }
        });
        selectedAtoms=selKeys.filter(function(k){
            return k!=null&&newIdx[k]!==undefined;
        }).map(function(k){return newIdx[k]});
        // The old centre stays the rotation centre (rebuildScene(true)): the
        // picture does not jump when the bounds move.
        rebuildScene(true);
        refreshSelectionText();
        mmSelectionChanged();
        layoutPanels();
        return true;
    }
    function refreshSelectionText(){
        if(selectedAtoms.length===0){selInfoEl.textContent='';return}
        var names=formatAtomList(selectedAtoms);
        if(currentMode==='replaceAtom')selInfoEl.textContent='Selected to replace: '+names+' (click Replace Atom button to confirm)';
        else if(currentMode==='addH')selInfoEl.textContent='Selected for Add H: '+names+' (click Add H button to complete)';
        else if(currentMode==='mmOpt')selInfoEl.textContent='Selected for MM optimization: '+names+' ('+selectedAtoms.length+' atoms, click Start Optimization)';
        else selInfoEl.textContent='Selected: '+names;
    }
    ['bnd-a-min','bnd-a-max','bnd-b-min','bnd-b-max','bnd-c-min','bnd-c-max'].forEach(function(id){
        var inp=document.getElementById(id);
        if(inp){
            inp.addEventListener('change',function(){
                var v=parseFloat(inp.value);
                if(isNaN(v))v=0;
                v=Math.max(-5,Math.min(5,v));
                inp.value=v;
                var nb={aMin:boundary.aMin,aMax:boundary.aMax,bMin:boundary.bMin,bMax:boundary.bMax,cMin:boundary.cMin,cMax:boundary.cMax};
                if(id==='bnd-a-min')nb.aMin=v;
                else if(id==='bnd-a-max')nb.aMax=v;
                else if(id==='bnd-b-min')nb.bMin=v;
                else if(id==='bnd-b-max')nb.bMax=v;
                else if(id==='bnd-c-min')nb.cMin=v;
                else if(id==='bnd-c-max')nb.cMax=v;
                if(nb.aMax<nb.aMin)nb.aMax=nb.aMin;
                if(nb.bMax<nb.bMin)nb.bMax=nb.bMin;
                if(nb.cMax<nb.cMin)nb.cMax=nb.cMin;
                applyBoundary(nb);
            });
        }
    });
    // N×N×N presets: the standard origin-aligned supercell, 0..N along every
    // axis (a 2×2×2 preset is eight unit cells, not nine).
    [['bnd-p2',2],['bnd-p3',3],['bnd-p4',4]].forEach(function(pr){
        var btn=document.getElementById(pr[0]);
        if(btn)btn.addEventListener('click',function(){
            applyBoundary({aMin:0,aMax:pr[1],bMin:0,bMax:pr[1],cMin:0,cMax:pr[1]});
        });
    });
    var bndResetBtn=document.getElementById('bnd-reset');
    if(bndResetBtn)bndResetBtn.addEventListener('click',function(){
        applyBoundary({aMin:0,aMax:1,bMin:0,bMax:1,cMin:0,cMax:1});
    });
    var bndOriginBtn=document.getElementById('bnd-origin');
    if(bndOriginBtn)bndOriginBtn.addEventListener('click',function(){
        // Same spans, low corner moved onto the lattice origin: the supercell
        // then starts at the crystallographic 0,0,0 instead of wherever the
        // slab happened to be floating.
        applyBoundary({aMin:0,aMax:boundary.aMax-boundary.aMin,bMin:0,bMax:boundary.bMax-boundary.bMin,cMin:0,cMax:boundary.cMax-boundary.cMin});
    });
    // Edit Cell: the six parameters are the one description of the box a
    // reader can check against the file, so they are shown back as typed
    // values. A site's fractional coordinates are the invariant across the
    // edit — that is what makes the stretch along a bond along a, and only
    // along a.
    var editCellBtn=document.getElementById('bnd-edit-cell');
    if(editCellBtn)editCellBtn.addEventListener('click',function(){
        if(!CRY)return;
        var rows=[['cell-a','a (Å)',CRY.a],['cell-b','b (Å)',CRY.b],['cell-c','c (Å)',CRY.c],
                  ['cell-alpha','α (°)',CRY.alpha],['cell-beta','β (°)',CRY.beta],['cell-gamma','γ (°)',CRY.gamma]];
        var body='';
        rows.forEach(function(r){
            body+='<div class="modal-row"><label for="'+r[0]+'">'+r[1]+'</label>'+
                '<input type="number" id="'+r[0]+'" step="0.001" value="'+r[2].toFixed(6)+'"></div>';
        });
        showModal('<h3>Edit Cell</h3>'+
            '<div class="current-val">Every site keeps its fractional coordinates: distances along an axis that changes scale with it, perpendicular distances do not.</div>'+
            body+
            '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button>'+
            '<button class="mbtn mbtn-ok" id="m-ok">Apply</button></div>',null);
        document.getElementById('m-cancel').addEventListener('click',hideModal);
        document.getElementById('m-ok').addEventListener('click',function(){
            var val={};
            rows.forEach(function(r){val[r[0]]=parseFloat(document.getElementById(r[0]).value)});
            if(applyCellEdit(val['cell-a'],val['cell-b'],val['cell-c'],val['cell-alpha'],val['cell-beta'],val['cell-gamma']))hideModal();
        });
    });
    // Bond Valence (H7): a read-only reading aid. It counts contacts from the
    // cell's own geometry, including the periodic images the drawn cell does
    // not contain, so one formula unit of NaCl still reports six neighbours.
    var bvsBtn=document.getElementById('bvs-open');
    if(bvsBtn)bvsBtn.addEventListener('click',showBvsReport);
    syncCellPanel();
    // Disorder groups: All shows every alternative the file lists, a group
    // shows the ordered part plus that alternative. The switch redraws the same
    // cell through rebuildCrystal() like any edit, and it takes an undo step so
    // Undo returns to the previous view (the snapshot carries the group state).
    var disorderSel=document.getElementById('bnd-disorder');
    if(disorderSel){
        disorderSel.addEventListener('change',function(){
            if(!CRY)return;
            var want=disorderSel.value||null;
            if((CRY.disorderFilter||null)===want){syncDisorderPanel();return}
            if(moveDragActive||moveKeyActive||fragRotActive)stopMoveSession();
            pushUndo();
            var selKeys=selectedAtoms.map(function(i){
                var a=MD.atoms[i];
                return (a&&a.baseIdx!=null&&a.cellI!=null)?(a.baseIdx+'_'+a.cellI+'_'+a.cellJ+'_'+a.cellK):null;
            });
            CRY.disorderFilter=want;
            rebuildCrystal();
            // The variant on screen is what a saved file will contain, so the
            // title says which alternative it shows; All is the file's own
            // state and drops the clause again.
            if(CRY.disorderFilter){
                provSet('disorderGroup',CRY.disorderFilter);
                if(disorderOccNormalize)provSet('occNormalize');
                else provUnmark('occNormalize');
            }else{
                provUnmark('disorderGroup');
                provUnmark('occNormalize');
            }
            var newIdx={};
            MD.atoms.forEach(function(a,i){
                if(a.baseIdx!=null&&a.cellI!=null){
                    var k=a.baseIdx+'_'+a.cellI+'_'+a.cellJ+'_'+a.cellK;
                    if(newIdx[k]===undefined)newIdx[k]=i;
                }
            });
            selectedAtoms=selKeys.filter(function(k){
                return k!=null&&newIdx[k]!==undefined;
            }).map(function(k){return newIdx[k]});
            rebuildScene(true);
            refreshSelectionText();
            mmSelectionChanged();
            layoutPanels();
        });
    }
    // The occupancy switch counts a group view as the ordered variant it
    // depicts. Nothing in the scene changes, so no rebuild and no undo step —
    // the info panel is re-read now and every export reads the choice when it
    // is written.
    var occNormChk=document.getElementById('bnd-occ-normalize');
    if(occNormChk){
        occNormChk.addEventListener('change',function(){
            disorderOccNormalize=!!occNormChk.checked;
            // The choice changes the cell's formula, density and occupancies
            // in every export, so the saved title states it; with no group on
            // screen it has no effect and no clause.
            if(disorderOccNormalize&&CRY&&CRY.disorderFilter)provSet('occNormalize');
            else provUnmark('occNormalize');
            updateMolInfo();
        });
    }
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
    // A run relaxes the geometry, so a ring's plane turns with it: the plan is
    // refreshed once per chunk (never per frame) and every piece reads it when
    // it is placed, so the dashes stay in their ring while the run streams.
    planAromaticDashes();
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
    // A ghost is a position, not a particle: it has no bonds and no electrons,
    // so there is nothing in it to relax. It is left out of the mobile set (the
    // reference geometry keeps its place) and handed to the force field as the
    // dummy element, whose well depth is zero - so it neither moves nor pulls
    // the atoms around it towards a point that is not there.
    var movable=selectedAtoms.filter(function(idx){return !isGhostAtom(MD.atoms[idx])});
    if(movable.length===0){
        modeInfoEl.textContent='Molecular Mechanics - the selection holds no real atom (ghosts are not particles)';
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
    var i,indices=movable,origin=new Array(indices.length);
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
        atoms:MD.atoms.map(function(a){return{element:isGhostAtom(a)?'X':a.element,x:a.x,y:a.y,z:a.z}}),
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

var MODE_INFO={view:'View Mode',bondLength:'Bond Length - Click 2 atoms',bondAngle:'Bond Angle - Click 3 atoms (central 2nd)',dihedral:'Dihedral - Click 4 atoms',neighbours:'Neighbours - Click an atom to list the atoms within the radius',interactions:'Interactions - Non-covalent contacts detected from geometry: H-bonds, sigma-hole contacts, pi-pi stacking and CH-pi; click a row to select its sites',addAtom:'Add Atom - Click anchor atom',deleteAtom:'Delete Atom - Click atom to delete',replaceAtom:'Replace Atom - Click atom(s) to select, click button again to confirm',addH:'Add H - Click atom(s) to select, click the button again to complete their hydrogens',selectAtoms:'Select Atoms - Input indices or element symbols',boxSelect:'Box Select - Drag a rectangle to select atoms',rotateGroup:'Rotate Group - Select group + axis atoms, then adjust angle',mmOpt:'Molecular Mechanics - Click atoms to select them; the panel optimizes only those atoms',moveAtoms:'Move Atoms - Select atoms; drag/arrows rotate view, right-drag pans; Ctrl+drag or Ctrl+arrows rotates group, Ctrl+right-drag moves group'};

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
    if(m==='addH'&&diffMode){
        // Two structures are on screen and the completion writes into one of
        // them; a comparison has no single structure to complete.
        modeInfoEl.textContent='Add H is not available in diff mode';
        return;
    }
    if(currentMode===m){
        // Pressing the mode's own toolbar button toggles its panel. A hidden
        // panel is a minimized one now, so hiding it leaves the chip that
        // brings it back rather than nothing at all.
        if(m==='selectAtoms'){
            var sp=document.getElementById('select-panel');
            if(sp)panelMinState.select=window.getComputedStyle(sp).display!=='none';
        }
        if(m==='rotateGroup'){
            var rp=document.getElementById('rotate-panel');
            if(rp)panelMinState.rotate=window.getComputedStyle(rp).display!=='none';
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
    // Entering a mode opens its panels again: a chip folded away in an earlier
    // visit belongs to that visit, not to this one.
    if(m==='view'){panelMinState.view=false;panelMinState.bake=false}
    if(m==='selectAtoms')panelMinState.select=false;
    if(m==='neighbours')panelMinState.neigh=false;
    if(m==='interactions')panelMinState.int=false;
    if(m==='rotateGroup')panelMinState.rotate=false;
    bakeHintShown=false;
    updateBakeHint();
    var preserveSel=((oldMode==='selectAtoms'||oldMode==='boxSelect'||oldMode==='replaceAtom'||oldMode==='addH'||oldMode==='rotateGroup'||oldMode==='moveAtoms'||oldMode==='mmOpt'||oldMode==='neighbours'||oldMode==='interactions')&&(m==='selectAtoms'||m==='boxSelect'||m==='replaceAtom'||m==='addH'||m==='rotateGroup'||m==='moveAtoms'||m==='mmOpt'||m==='neighbours'||m==='interactions'));
    if(oldMode==='moveAtoms'&&m!=='moveAtoms')stopMoveSession();
    if(oldMode==='mmOpt'&&m!=='mmOpt')mmCancelRun();
    // Auto Rotate is a View-mode control: leaving the mode takes its panel away,
    // so the spin it started ends with it rather than running on with no way
    // back to the switch.
    if(oldMode==='view'&&m!=='view')stopSpin();
    // The contact lines belong to the Neighbours panel: leaving the mode takes
    // both away, and the panel starts empty the next time it is opened.
    if(oldMode==='neighbours'&&m!=='neighbours')neighEndSession();
    if(oldMode==='interactions'&&m!=='interactions')intEndSession();
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
    var np2=document.getElementById('neigh-panel');
    if(np2)np2.style.display=(m==='neighbours')?'block':'none';
    if(m==='neighbours'&&oldMode!=='neighbours')neighStartSession();
    var ip2=document.getElementById('int-panel');
    if(ip2)ip2.style.display=(m==='interactions')?'block':'none';
    if(m==='interactions'&&oldMode!=='interactions')intStartSession();
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
        // A previewed angle is the geometry the user leaves the mode on, so it
        // is carried into the base cell before the mode's state is dropped.
        if(rotActive&&CRY)commitCrystalCoords(rotMovedDisplayAtoms());
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
// Where in the undo stack this rotation session began. Reset takes the session
// back as a whole and needs its own first snapshot, not an earlier session's.
var rotUndoMark=0;

// The atoms a rotation actually turned (the group minus its axis): what a
// crystal write-back has to carry into the base cell for this step.
function rotMovedDisplayAtoms(){
    var axisSet={};
    rotateAxisAtoms.forEach(function(i){axisSet[i]=true});
    var out=[];
    for(var i=0;i<rotateGroupAtoms.length;i++){
        var idx=rotateGroupAtoms[i];
        if(axisSet[idx])continue;
        out.push(idx);
    }
    return out;
}
// The base cell this rotation session started from: the earliest snapshot the
// session pushed — taken at its first committed angle, before that angle's
// geometry reached the base. Reset drops the session's snapshots, so this has
// to be read before they are popped.
function restoreRotSessionBase(){
    if(!CRY)return false;
    var from=rotUndoMark>0?rotUndoMark:0;
    for(var i=from;i<undoStack.length;i++){
        var s=undoStack[i];
        if(s.rotAngle!=null&&s.baseAtoms&&s.baseBonds){
            CRY.baseAtoms=s.baseAtoms;CRY.baseBonds=s.baseBonds;
            return true;
        }
    }
    return false;
}
function resetRotAxisState(){
    // A previewed angle still on screen when the session ends (Done, Clear, a
    // mode switch) is its last word: carry it into the base cell before the set
    // of moved atoms is cleared.
    if(rotActive&&CRY)commitCrystalCoords(rotMovedDisplayAtoms());
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
        // The new angle is the geometry the file has to keep, and a preview
        // alone never reaches the base cell: write this step back now.
        if(CRY)commitCrystalCoords(rotMovedDisplayAtoms());
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
        // steps back one angle adjustment at a time. The mark is where Reset
        // starts looking for the session's own first snapshot.
        rotUndoMark=undoStack.length;
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
        // A crystal's base cell carries the session's committed angles, so Reset
        // restores it from the session's own first snapshot as well — the display
        // is not the whole state any more.
        if(restoreRotSessionBase())rebuildCrystal();
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

function pushUndo(symmetryState){
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
    var snap={hasExplicitBonds:!!MD.hasExplicitBonds,prov:provClone()};
    if(CRY){
        // A crystal's picture is a pure function of three things — the base
        // cell, its bonds and the bounds — so those three are the snapshot, and
        // the displayed supercell is not stored at all: a 2x2x2 view then costs
        // exactly what a 1x1x1 view does, and doUndo() rebuilds the picture from
        // the cell it restores.
        snap.baseAtoms=JSON.parse(JSON.stringify(CRY.baseAtoms));
        snap.baseBonds=JSON.parse(JSON.stringify(CRY.baseBonds));
        snap.boundary={aMin:boundary.aMin,aMax:boundary.aMax,bMin:boundary.bMin,bMax:boundary.bMax,cMin:boundary.cMin,cMax:boundary.cMax};
        // The cell is part of the state an edit can change (Edit Cell), and a
        // snapshot that restored the atoms but not the box would put them back
        // into a different cell than the one they were cut from.
        snap.cell={a:CRY.a,b:CRY.b,c:CRY.c,alpha:CRY.alpha,beta:CRY.beta,gamma:CRY.gamma,
            latticeVectors:JSON.parse(JSON.stringify(CRY.latticeVectors))};
        // The disorder selection is view state the base cell cannot rebuild
        // from itself, so an undo of a group switch must put it back too. The
        // occupancy switch changes what a group view exports, so it is part of
        // the same state: a snapshot that restored the group but not the way
        // its sites are counted would save a different composition than the
        // one on screen when the snapshot was taken.
        snap.disorderFilter=CRY.disorderFilter||null;
        snap.disorderOccNormalize=!!disorderOccNormalize;
        // Only a symmetry switch passes this: the operations and sites are not
        // part of an ordinary edit, and snapshotting them every time would
        // double the undo memory for a state that never changes.
        if(symmetryState)snap.symmetry=symmetryState;
    }else{
        snap.atoms=JSON.parse(JSON.stringify(MD.atoms));
        snap.bonds=JSON.parse(JSON.stringify(MD.bonds));
    }
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
    MD.hasExplicitBonds=!!snap.hasExplicitBonds;
    // The provenance belongs to the geometry: an undone edit is not in the
    // structure any more, so its clause comes off the title with it.
    provRestore(snap.prov);
    if(CRY&&snap.baseAtoms&&snap.baseBonds){
        // A crystal snapshot holds the base cell and the bounds, and the picture
        // is rebuilt from them — this is the step that keeps a large supercell's
        // undo as cheap as a unit cell's.
        CRY.baseAtoms=snap.baseAtoms;CRY.baseBonds=snap.baseBonds;
        if(snap.boundary){
            boundary.aMin=snap.boundary.aMin;boundary.aMax=snap.boundary.aMax;
            boundary.bMin=snap.boundary.bMin;boundary.bMax=snap.boundary.bMax;
            boundary.cMin=snap.boundary.cMin;boundary.cMax=snap.boundary.cMax;
            syncBndInputs();
        }
        if(snap.cell){
            CRY.a=snap.cell.a;CRY.b=snap.cell.b;CRY.c=snap.cell.c;
            CRY.alpha=snap.cell.alpha;CRY.beta=snap.cell.beta;CRY.gamma=snap.cell.gamma;
            CRY.latticeVectors=snap.cell.latticeVectors;
            computeCrystalInv();
        }
        // A snapshot taken before the disorder filter existed restores All,
        // which is what it was showing.
        CRY.disorderFilter=snap.disorderFilter||null;
        // Same for the occupancy switch: it was on unless the snapshot says
        // otherwise, and rebuildCrystal() below re-syncs the checkbox to it.
        disorderOccNormalize=snap.disorderOccNormalize!==false;
        // A symmetry switch changes the operation set and the sites together,
        // so undoing one must restore them together — the atoms in the snapshot
        // were expanded under the operations it also carries.
        if(snap.symmetry){
            CRY.symmetryOps=snap.symmetry.ops;
            CRY.symmetrySource=snap.symmetry.source;
            CRY.asymAtoms=snap.symmetry.asymAtoms;
            CRY.fileBonds=snap.symmetry.fileBonds;
            CRY.bondsDetected=snap.symmetry.bondsDetected;
            CRY.symmetryCheck=snap.symmetry.check;
            CRY.symmetryAlternative=snap.symmetry.alternative;
        }
        syncCellPanel();
        rebuildCrystal();
    }else{
        MD.atoms=snap.atoms;MD.bonds=snap.bonds;
    }
    rebuildScene();updateUndoBtn();
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
// The comparison's right-hand hulls, in the same style the left half draws.
var diffPolyFacesMesh=null;
var diffPolyEdgesMesh=null;
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
// A crystal comparison runs on the base cells, so the difference list's indices
// are base indices, the right-half display index is the matching base index, and
// the left half has one displayed image per base atom (more in a supercell).
// These four carry that translation; the molecular comparison leaves them at
// their defaults.
var diffBaseIndices=false;
var diffHasCrystal=false;
var diffLabelAtoms=null;
var diffDisplayByBase=null;
var diffReverseBaseMapping=null;
var diffRightCrystal=null;
var diffRightInv=null;
var diffRightBonds=null;
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
        var item={element:a.element,x:a.x+t.x,y:a.y+t.y,z:a.z+t.z,color:a.color||((MD.atomColors&&MD.atomColors[a.element])||'#FF1493')};
        // An imported ghost stays a ghost - its element is what says which basis
        // it stands for, and the flag is what keeps it unbonded.
        if(a.ghost===true||a.element==='X')item.ghost=true;
        MD.atoms.push(item);
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
                // A ghost can carry no bond: an authored bond naming one in the
                // imported file is dropped here as it is on load.
                if(newAtoms[b.atom1]&&isGhostAtom(newAtoms[b.atom1]))return;
                if(newAtoms[b.atom2]&&isGhostAtom(newAtoms[b.atom2]))return;
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
var PASTE_ATOM_LIMIT=30000;   // the same ceiling as MAX_CRYSTAL_ATOMS: a crystal that can be cut can be pasted
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
        // The ghost flag is carried through the clipboard, but only as a real
        // boolean: anything else in that field is data from an older version or
        // another tab and is dropped rather than trusted.
        var item={element:a.element,x:a.x,y:a.y,z:a.z};
        if(a.ghost===true)item.ghost=true;
        atoms.push(item);
    }
    var bonds=[],srcB=Array.isArray(c.bonds)?c.bonds:[];
    for(i=0;i<srcB.length;i++){
        var b=srcB[i];
        if(!b||typeof b.a!=='number'||typeof b.b!=='number')continue;
        if(b.a<0||b.b<0||b.a>=atoms.length||b.b>=atoms.length||b.a===b.b)continue;
        // A ghost cannot carry a bond, so a stored fragment that has one is
        // stripped of it here rather than pasted as something the viewer would
        // refuse to draw.
        if(isGhostAtom(atoms[b.a])||isGhostAtom(atoms[b.b]))continue;
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
        var item={element:a.element,x:a.x,y:a.y,z:a.z};
        if(isGhostAtom(a))item.ghost=true;
        atoms.push(item);
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
        var item={element:a.element,x:a.x+t.x,y:a.y+t.y,z:a.z+t.z,
            color:(MD.atomColors&&MD.atomColors[a.element])||'#FF1493'};
        if(a.ghost===true||a.element==='X')item.ghost=true;
        MD.atoms.push(item);
    });
    frag.bonds.forEach(function(b){
        // A fragment read back from storage may still name a ghost in a bond;
        // a ghost has no valence, so such a bond is not pasted at all.
        if(isGhostAtom(frag.atoms[b.a])||isGhostAtom(frag.atoms[b.b]))return;
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
    if(CRY){startCrystalAddHydrogens();return}
    if(selectedAtoms.length===0){
        modeInfoEl.textContent='Add H - select atoms first: click atoms here, or use Select Atoms / Box Select';
        return;
    }
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    // A ghost has no valence, so it can be neither a centre (nothing to
    // complete) nor a neighbour (no bond to count). It goes to the engine as the
    // dummy element, and out of the selection.
    var hSel=selectedAtoms.filter(function(i){return !isGhostAtom(MD.atoms[i])});
    if(hSel.length===0){
        modeInfoEl.textContent='Add H - the selection holds no real atom (ghosts have no valence)';
        return;
    }
    addHPending=true;
    modeInfoEl.textContent='Add H: completing '+hSel.length+' selected atom(s)...';
    vscodeApi.postMessage({
        command:'addHydrogens',
        atoms:MD.atoms.map(function(a){return{element:isGhostAtom(a)?'X':a.element,x:a.x,y:a.y,z:a.z}}),
        bonds:MD.bonds.map(function(b){return{atom1:b.atom1,atom2:b.atom2,order:b.order}}),
        selected:hSel
    });
}
function applyHydrogenResult(msg){
    addHPending=false;
    if(addHCrystalCtx){
        var cryCtx=addHCrystalCtx;
        addHCrystalCtx=null;
        applyCrystalHydrogenResult(msg,cryCtx);
        return;
    }
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

// ---------------------------------------------------------------------------
// Crystal Add H - the periodic neighbourhood
// ---------------------------------------------------------------------------
// A crystal's completion cannot run on the base cell alone: the hydroxyl
// oxygen of a molecule that straddles a cell face stores its carbon a whole
// cell away, and an engine handed that pair alone would see no bond at all.
// The request therefore carries a local cluster of the selected image - the
// centre, its bonded neighbours out to two bonds (the reach the engine's
// hybridization pass reads), and every atom image within the contact radius of
// the centre (the reach its clash checks read) - laid out in the cell frame of
// the image the user picked. The host's molecular engine then places the
// hydrogens exactly as it does for a molecule, and the answer comes back in
// that same frame.
//
// Each new hydrogen is wrapped into the base cell and stored as a new base
// site; the bond that joins it to its centre carries the cell shift between the
// two images - bond {a1,a2,s} means "image c of a1 bonds image c+s of a2", the
// convention every other crystal editor here reads. The base count is checked
// when the answer lands, so an edit made while the request was in flight cannot
// write the hydrogens into a different structure.
var CRYSTAL_H_BOND_DEPTH=2;
var CRYSTAL_H_CONTACT_RADIUS=6.0;
var addHCrystalCtx=null;

// The base cell as the picture shows it: with a disorder group selected the
// hidden alternatives take no part (a hydrogen must not complete a neighbour
// that is not on screen), and view indices are mapped back to the base indices
// every stored bond and new site addresses.
function crystalHBaseView(){
    var baseAtoms=CRY.baseAtoms,baseBonds=CRY.baseBonds;
    if(!CRY.disorderFilter)return{atoms:baseAtoms,bonds:baseBonds,toBase:null};
    var visMap=new Array(baseAtoms.length),atoms=[],toBase=[];
    for(var i=0;i<baseAtoms.length;i++){
        var ba=baseAtoms[i];
        if(!ba.disorderGroup||ba.disorderGroup===CRY.disorderFilter){
            visMap[i]=atoms.length;toBase.push(i);atoms.push(ba);
        }else visMap[i]=-1;
    }
    var bonds=[];
    for(var b=0;b<baseBonds.length;b++){
        var bb=baseBonds[b];
        if(visMap[bb.atom1]>=0&&visMap[bb.atom2]>=0){
            bonds.push({atom1:visMap[bb.atom1],atom2:visMap[bb.atom2],order:bb.order,shift:bb.shift});
        }
    }
    return{atoms:atoms,bonds:bonds,toBase:toBase};
}

// The cluster one request runs on: atoms carry Cartesian coordinates in the
// selected image's cell frame and bonds the indices of the atoms here.
function crystalHCluster(view,anchors){
    var lv=CRY.latticeVectors,inv=CRY_INV;
    var widths=crystalPlaneWidths(lv);
    var baseAtoms=view.atoms,baseBonds=view.bonds;
    var nodes=[],nodeOf={};
    function keyOf(b,c){return b+'_'+c[0]+'_'+c[1]+'_'+c[2]}
    function addNode(b,c){
        var key=keyOf(b,c);
        var at=nodeOf[key];
        if(at!==undefined)return at;
        var ba=baseAtoms[b];
        if(!ba)return -1;
        at=nodes.length;
        nodes.push({base:b,cell:[c[0],c[1],c[2]],element:ba.element,
            x:ba.x+c[0]*lv[0][0]+c[1]*lv[1][0]+c[2]*lv[2][0],
            y:ba.y+c[0]*lv[0][1]+c[1]*lv[1][1]+c[2]*lv[2][1],
            z:ba.z+c[0]*lv[0][2]+c[1]*lv[1][2]+c[2]*lv[2][2]});
        nodeOf[key]=at;
        return at;
    }
    var adj=[];
    for(var i=0;i<baseAtoms.length;i++)adj[i]=[];
    baseBonds.forEach(function(b){
        if(adj[b.atom1])adj[b.atom1].push({other:b.atom2,sign:1,bond:b});
        if(adj[b.atom2])adj[b.atom2].push({other:b.atom1,sign:-1,bond:b});
    });
    var outAnchors=[];
    var R2=CRYSTAL_H_CONTACT_RADIUS*CRYSTAL_H_CONTACT_RADIUS;
    anchors.forEach(function(an){
        var root=addNode(an.base,an.cell);
        if(root<0)return;
        outAnchors.push({base:an.base,cell:[an.cell[0],an.cell[1],an.cell[2]],node:root});
        // Bonded neighbours, followed through the shifts the cell records
        // rather than the raw distances: a wrapped molecule's bond points a
        // cell away, and the walk has to reach the image it names.
        var depthAt={};depthAt[root]=0;
        var queue=[root];
        for(var qi=0;qi<queue.length;qi++){
            var cur=queue[qi];
            if(depthAt[cur]>=CRYSTAL_H_BOND_DEPTH)continue;
            var cn=nodes[cur];
            var edges=adj[cn.base]||[];
            for(var ei=0;ei<edges.length;ei++){
                var ed=edges[ei];
                var s=ed.bond.shift||ZERO_SHIFT;
                var nc=[cn.cell[0]+ed.sign*s[0],cn.cell[1]+ed.sign*s[1],cn.cell[2]+ed.sign*s[2]];
                var nb=addNode(ed.other,nc);
                if(nb<0||depthAt[nb]!==undefined)continue;
                depthAt[nb]=depthAt[cur]+1;
                queue.push(nb);
            }
        }
        // Non-bonded images within the clash radius. The nearest image comes
        // from the same enumerated search the bond detector uses (rounding is
        // only the nearest image in an orthogonal cell), and the images one
        // more cell out are checked too so a small cell's second contact is
        // not missed.
        var rn=nodes[root];
        var ab=baseAtoms[an.base];
        var afx=inv[0][0]*ab.x+inv[0][1]*ab.y+inv[0][2]*ab.z;
        var afy=inv[1][0]*ab.x+inv[1][1]*ab.y+inv[1][2]*ab.z;
        var afz=inv[2][0]*ab.x+inv[2][1]*ab.y+inv[2][2]*ab.z;
        for(var b2=0;b2<baseAtoms.length;b2++){
            var ba2=baseAtoms[b2];
            var dfx=inv[0][0]*ba2.x+inv[0][1]*ba2.y+inv[0][2]*ba2.z-afx;
            var dfy=inv[1][0]*ba2.x+inv[1][1]*ba2.y+inv[1][2]*ba2.z-afy;
            var dfz=inv[2][0]*ba2.x+inv[2][1]*ba2.y+inv[2][2]*ba2.z-afz;
            var near=crystalShortestFractional(dfx,dfy,dfz,CRYSTAL_H_CONTACT_RADIUS,lv,widths);
            if(near.dist>CRYSTAL_H_CONTACT_RADIUS)continue;
            var c0=an.cell[0]-near.shift[0],c1=an.cell[1]-near.shift[1],c2=an.cell[2]-near.shift[2];
            for(var ox=-2;ox<=2;ox++){
                for(var oy=-2;oy<=2;oy++){
                    for(var oz=-2;oz<=2;oz++){
                        var cell=[c0+ox,c1+oy,c2+oz];
                        if(nodeOf[keyOf(b2,cell)]!==undefined)continue;
                        var px=ba2.x+cell[0]*lv[0][0]+cell[1]*lv[1][0]+cell[2]*lv[2][0];
                        var py=ba2.y+cell[0]*lv[0][1]+cell[1]*lv[1][1]+cell[2]*lv[2][1];
                        var pz=ba2.z+cell[0]*lv[0][2]+cell[1]*lv[1][2]+cell[2]*lv[2][2];
                        var ddx=px-rn.x,ddy=py-rn.y,ddz=pz-rn.z;
                        if(ddx*ddx+ddy*ddy+ddz*ddz<=R2)addNode(b2,cell);
                    }
                }
            }
        }
    });
    if(!outAnchors.length)return null;
    // Every bond of the base list whose two images are both in the cluster.
    // The centre's valence and hybridisation are read off these, so a bond the
    // walk followed above must appear here too.
    var bonds=[],seen={};
    for(var ni=0;ni<nodes.length;ni++){
        var nn=nodes[ni];
        var edges2=adj[nn.base]||[];
        for(var e2=0;e2<edges2.length;e2++){
            var ed2=edges2[e2];
            var s2=ed2.bond.shift||ZERO_SHIFT;
            var j=nodeOf[keyOf(ed2.other,[nn.cell[0]+ed2.sign*s2[0],nn.cell[1]+ed2.sign*s2[1],nn.cell[2]+ed2.sign*s2[2]])];
            if(j===undefined||j===ni)continue;
            var pk=ni<j?ni+'_'+j:j+'_'+ni;
            if(seen[pk])continue;
            seen[pk]=true;
            bonds.push({atom1:ni,atom2:j,order:ed2.bond.order==null?1:ed2.bond.order});
        }
    }
    return{atoms:nodes,bonds:bonds,anchors:outAnchors};
}

function startCrystalAddHydrogens(){
    if(!CRY_INV){modeInfoEl.textContent='Add H is not available: the cell has no usable lattice';return}
    if(selectedAtoms.length===0){
        modeInfoEl.textContent='Add H - select atoms first: click atoms here, or use Select Atoms / Box Select';
        return;
    }
    var view=crystalHBaseView();
    var visOf={};
    for(var vi=0;vi<view.atoms.length;vi++){
        visOf[view.toBase?view.toBase[vi]:vi]=vi;
    }
    // One anchor per base site: a supercell draws several images of one site
    // and the file holds the site once, so completing each image would append
    // two hydrogens to the same base position.
    var anchors=[],seen={};
    for(var s=0;s<selectedAtoms.length;s++){
        var a=MD.atoms[selectedAtoms[s]];
        if(!a||a.baseIdx==null||isGhostAtom(a))continue;
        if(seen[a.baseIdx])continue;
        var viewIdx=visOf[a.baseIdx];
        if(viewIdx===undefined)continue;
        seen[a.baseIdx]=true;
        var g=a.imgShift||ZERO_SHIFT;
        anchors.push({base:viewIdx,cell:[g[0],g[1],g[2]]});
    }
    if(anchors.length===0){
        modeInfoEl.textContent='Add H - the selection holds no crystal site to complete (ghosts have no valence)';
        return;
    }
    if(vibActive||vibPaused)stopVibration();
    if(moveDragActive||fragRotActive||moveKeyActive)stopMoveSession();
    var cluster=crystalHCluster(view,anchors);
    if(!cluster){
        modeInfoEl.textContent='Add H: the selected site has no periodic neighbourhood to complete';
        return;
    }
    addHPending=true;
    addHCrystalCtx={cry:CRY,baseCount:CRY.baseAtoms.length,toBase:view.toBase,
        atomCount:cluster.atoms.length,anchors:cluster.anchors};
    modeInfoEl.textContent='Add H: completing '+anchors.length+' selected site'+(anchors.length===1?'':'s')+'...';
    vscodeApi.postMessage({
        command:'addHydrogens',
        atoms:cluster.atoms.map(function(a){return{element:a.element,x:a.x,y:a.y,z:a.z}}),
        bonds:cluster.bonds.map(function(b){return{atom1:b.atom1,atom2:b.atom2,order:b.order}}),
        selected:cluster.anchors.map(function(a){return a.node})
    });
}

function applyCrystalHydrogenResult(msg,ctx){
    if(!msg||!msg.ok){modeInfoEl.textContent=(msg&&msg.message)||'Add H failed';return}
    if(!CRY||CRY!==ctx.cry)return;
    if(CRY.baseAtoms.length!==ctx.baseCount||
       (typeof msg.atomCount==='number'&&msg.atomCount!==ctx.atomCount)){
        modeInfoEl.textContent='Add H: the structure changed while computing - nothing was added';
        return;
    }
    var newAtoms=Array.isArray(msg.atoms)?msg.atoms:[];
    var newBonds=Array.isArray(msg.bonds)?msg.bonds:[];
    if(newAtoms.length===0){modeInfoEl.textContent=msg.message||'Add H: nothing to add';return}
    var anchorOf={};
    ctx.anchors.forEach(function(a){anchorOf[a.node]=a});
    var inv=CRY_INV,lv=CRY.latticeVectors;
    if(!inv){modeInfoEl.textContent='Add H: the cell has no usable lattice';return}
    var adds=[];
    for(var i=0;i<newBonds.length;i++){
        var b=newBonds[i];
        if(!b)continue;
        var an=anchorOf[b.atom1];
        if(!an)continue;
        var h=newAtoms[b.atom2-ctx.atomCount];
        if(!h)continue;
        var fx=inv[0][0]*h.x+inv[0][1]*h.y+inv[0][2]*h.z;
        var fy=inv[1][0]*h.x+inv[1][1]*h.y+inv[1][2]*h.z;
        var fz=inv[2][0]*h.x+inv[2][1]*h.y+inv[2][2]*h.z;
        if(!isFinite(fx)||!isFinite(fy)||!isFinite(fz))continue;
        var mx=Math.floor(fx),my=Math.floor(fy),mz=Math.floor(fz);
        var wx=fx-mx,wy=fy-my,wz=fz-mz;
        var x=wx*lv[0][0]+wy*lv[1][0]+wz*lv[2][0];
        var y=wx*lv[0][1]+wy*lv[1][1]+wz*lv[2][1];
        var z=wx*lv[0][2]+wy*lv[1][2]+wz*lv[2][2];
        var baseIdx=ctx.toBase?ctx.toBase[an.base]:an.base;
        if(!CRY.baseAtoms[baseIdx])continue;
        // A site that already carries a hydrogen at this base position is not
        // given a second one: the H its bond would name is already there.
        var dup=false;
        for(var k=0;k<CRY.baseAtoms.length;k++){
            var ba=CRY.baseAtoms[k];
            if(ba.element!=='H')continue;
            var ddx=ba.x-x,ddy=ba.y-y,ddz=ba.z-z;
            if(ddx*ddx+ddy*ddy+ddz*ddz<0.25){dup=true;break}
        }
        for(var a2=0;a2<adds.length&&!dup;a2++){
            var ad=adds[a2];
            var ex=ad.x-x,ey=ad.y-y,ez=ad.z-z;
            if(ex*ex+ey*ey+ez*ez<0.25)dup=true;
        }
        if(dup)continue;
        // The image of the new site that the bond points at: the one at the
        // spot the host placed, measured from the selected image's cell.
        adds.push({base:baseIdx,x:x,y:y,z:z,
            shift:[mx-an.cell[0],my-an.cell[1],mz-an.cell[2]]});
    }
    if(adds.length===0){modeInfoEl.textContent=msg.message||'Add H: nothing to add';return}
    // The answer is applied in one step: a single snapshot, one rebuild.
    pushUndo();
    provMark('addH',adds.length);
    var color=(MD.atomColors&&MD.atomColors.H)||'#FFFFFF';
    adds.forEach(function(ad){
        var ni=CRY.baseAtoms.length;
        CRY.baseAtoms.push({element:'H',x:ad.x,y:ad.y,z:ad.z,index:ni,baseIdx:ni,occupancy:1,color:color});
        var desc={atom1:ad.base,atom2:ni,shift:ad.shift};
        if(!crystalFindBaseBond(desc))CRY.baseBonds.push({atom1:ad.base,atom2:ni,order:1,shift:ad.shift});
    });
    rebuildCrystal();
    rebuildScene();
    var text='Add H: added '+adds.length+' H to '+ctx.anchors.length+' crystal site'+
        (ctx.anchors.length===1?'':'s');
    if(adds.length<newAtoms.length)text+=' ('+(newAtoms.length-adds.length)+' already present)';
    modeInfoEl.textContent=text;
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

    // When the two lists are one graph in the same order -- the same file
    // compared with itself, or a file and its unedited copy -- the identity is
    // the correspondence a user means. Answering it directly is also what keeps
    // a self-comparison exact: the backtracking search below may legitimately
    // pick a symmetry permutation instead, and a swapped pair of equivalent
    // atoms would then show a difference the user never made. A reordered file
    // fails the neighbour test and takes the general search.
    var sameOrder=true,i;
    for(i=0;i<n1;i++){
        if(atoms1[i].element!==atoms2[i].element){sameOrder=false;break}
    }
    if(sameOrder){
        var neighbourKey=function(adj,idx){
            var list=adj[idx].map(function(e){return e.to});
            list.sort(function(a,b){return a-b});
            return list.join(',');
        };
        for(i=0;i<n1;i++){
            if(neighbourKey(adj1,i)!==neighbourKey(adj2,i)){sameOrder=false;break}
        }
        if(sameOrder)return atoms1.map(function(_,idx){return idx});
    }

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

// ===== Crystal comparison: the cell's periodic graph, not the displayed copies =====
//
// A crystal's screen is a supercell whose molecules are drawn whole, while the
// file it is compared against is one wrapped cell: at 2x2x2 the two lists are
// 2352 and 294 atoms, and even at 0..1 the same atom sits at a different
// cartesian position. Comparing those lists could only ever answer "Skeletons
// Differ" -- which is exactly what it did. A comparison now happens on what both
// files actually state: their **cell parameters** and their **base cell's
// periodic graph**. Every distance is measured through the bond's own image, so
// a bond across a cell face is its real length and not a cell.
var CRYSTAL_DIFF_LEN_TOL=1e-2;   // Angstrom: two cells this close are the same cell
var CRYSTAL_DIFF_ANG_TOL=0.05;   // degree

// The six cell parameters derived from the lattice vectors themselves (a file's
// a/b/c/angles are the same numbers, but the vectors are the ones every other
// crystal pass uses).
function crystalCellParams(cry){
    var lv=cry&&cry.latticeVectors;
    if(!lv||lv.length<3)return null;
    function len(v){return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])}
    function ang(u,v){
        var lu=len(u),lvv=len(v);
        if(!(lu>1e-12)||!(lvv>1e-12))return 90;
        var d=(u[0]*v[0]+u[1]*v[1]+u[2]*v[2])/(lu*lvv);
        return Math.acos(Math.max(-1,Math.min(1,d)))*180/Math.PI;
    }
    return {a:len(lv[0]),b:len(lv[1]),c:len(lv[2]),
            alpha:ang(lv[1],lv[2]),beta:ang(lv[0],lv[2]),gamma:ang(lv[0],lv[1])};
}
function crystalCellParamsEqual(p,q){
    if(!p||!q)return false;
    var lens=['a','b','c'],angs=['alpha','beta','gamma'],i;
    for(i=0;i<3;i++)if(Math.abs(p[lens[i]]-q[lens[i]])>CRYSTAL_DIFF_LEN_TOL)return false;
    for(i=0;i<3;i++)if(Math.abs(p[angs[i]]-q[angs[i]])>CRYSTAL_DIFF_ANG_TOL)return false;
    return true;
}
function crystalCellDiffHTML(p,q){
    var rows=[['a','a'],['b','b'],['c','c'],['alpha','α'],['beta','β'],['gamma','γ']];
    var html='';
    for(var i=0;i<rows.length;i++){
        var k=rows[i][0],name=rows[i][1];
        var unit=(k==='a'||k==='b'||k==='c')?' Å':'°';
        var dec=(k==='a'||k==='b'||k==='c')?3:2;
        html+='<div class="diff-row">'+name+': '+(p?p[k].toFixed(dec):'?')+unit+
              ' vs '+(q?q[k].toFixed(dec):'?')+unit+'</div>';
    }
    return html;
}
// The inverse lattice of any cell. computeCrystalInv() reads the loaded CRY;
// the right-hand cell of a comparison has no business writing into it.
function crystalInverseOf(lv){
    if(!lv||lv.length<3)return null;
    var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-
            lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+
            lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
    if(Math.abs(det)<1e-12)return null;
    var invDet=1/det;
    return [
        [(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])*invDet,(lv[1][2]*lv[2][0]-lv[1][0]*lv[2][2])*invDet,(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0])*invDet],
        [(lv[0][2]*lv[2][1]-lv[0][1]*lv[2][2])*invDet,(lv[0][0]*lv[2][2]-lv[0][2]*lv[2][0])*invDet,(lv[0][1]*lv[2][0]-lv[0][0]*lv[2][1])*invDet],
        [(lv[0][1]*lv[1][2]-lv[0][2]*lv[1][1])*invDet,(lv[0][2]*lv[1][0]-lv[0][0]*lv[1][2])*invDet,(lv[0][0]*lv[1][1]-lv[0][1]*lv[1][0])*invDet]
    ];
}
function crystalFracsOf(cry,inv){
    return (cry.baseAtoms||[]).map(function(a){
        return [inv[0][0]*a.x+inv[0][1]*a.y+inv[0][2]*a.z,
                inv[1][0]*a.x+inv[1][1]*a.y+inv[1][2]*a.z,
                inv[2][0]*a.x+inv[2][1]*a.y+inv[2][2]*a.z];
    });
}
// One bond per base pair, keyed low-to-high, its shift moving the higher atom's
// image relative to the lower one.
function crystalBondMap(bonds){
    var map={};
    (bonds||[]).forEach(function(b){
        var i=b.atom1,j=b.atom2,s=b.shift||ZERO_SHIFT;
        if(i>j){var t=i;i=j;j=t;s=[-s[0],-s[1],-s[2]]}
        map[i+','+j]=s;
    });
    return map;
}
function crystalBondShift(map,from,to){
    if(from===to)return ZERO_SHIFT;
    var a=from,b=to,neg=false;
    if(a>b){var t=a;a=b;b=t;neg=true}
    var s=map[a+','+b];
    if(!s)return ZERO_SHIFT;
    return neg?[-s[0],-s[1],-s[2]]:s;
}
// The cartesian vector from an atom to the bonded image of another. This is what
// a bond length, an angle and a dihedral are measured through; a bond whose
// shift is (0,0,0) is the plain difference.
function crystalPairVector(cry,frs,map,from,to){
    var lv=cry.latticeVectors;
    var s=crystalBondShift(map,from,to);
    var fa=frs[from],fb=frs[to];
    var dx=fb[0]+s[0]-fa[0],dy=fb[1]+s[1]-fa[1],dz=fb[2]+s[2]-fa[2];
    return [dx*lv[0][0]+dy*lv[1][0]+dz*lv[2][0],
            dx*lv[0][1]+dy*lv[1][1]+dz*lv[2][1],
            dx*lv[0][2]+dy*lv[1][2]+dz*lv[2][2]];
}
function crystalVectorAngle(u,v){
    var du=Math.sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    var dv=Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if(du<1e-10||dv<1e-10)return 0;
    var d=(u[0]*v[0]+u[1]*v[1]+u[2]*v[2])/(du*dv);
    return Math.acos(Math.max(-1,Math.min(1,d)))*180/Math.PI;
}
// A dihedral along a chain of bonded images: the three bond vectors, laid end to
// end. Taking each bond's own image is what keeps a torsion across a face the
// torsion of the molecule, not of the cell.
function crystalVectorDihedral(v1,v2,v3){
    var p1={x:v1[0],y:v1[1],z:v1[2]};
    var p2={x:v1[0]+v2[0],y:v1[1]+v2[1],z:v1[2]+v2[2]};
    var p3={x:p2.x+v3[0],y:p2.y+v3[1],z:p2.z+v3[2]};
    return dihedral3d({x:0,y:0,z:0},p1,p2,p3);
}
// The six measurement closures computeConformationDiffs can run on instead of
// the plain cartesian ones.
function crystalGeometryOf(cry,inv,bonds){
    var frs=crystalFracsOf(cry,inv);
    var map=crystalBondMap(bonds);
    return {
        dist:function(i,j){
            var v=crystalPairVector(cry,frs,map,i,j);
            return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        },
        angle:function(i,j,k){
            var u=crystalPairVector(cry,frs,map,j,i);
            var v=crystalPairVector(cry,frs,map,j,k);
            return crystalVectorAngle(u,v);
        },
        dihedral:function(a,b,c,d){
            var w0=crystalPairVector(cry,frs,map,a,b);
            var w1=crystalPairVector(cry,frs,map,b,c);
            var w2=crystalPairVector(cry,frs,map,c,d);
            return crystalVectorDihedral(w0,w1,w2);
        }
    };
}
function computeCrystalConformationDiffs(leftCry,leftInv,leftBonds,rightCry,rightInv,rightBonds,mapping,threshPct){
    var g1=crystalGeometryOf(leftCry,leftInv,leftBonds);
    var g2=crystalGeometryOf(rightCry,rightInv,rightBonds);
    return computeConformationDiffs(leftCry.baseAtoms,leftBonds,rightCry.baseAtoms,rightBonds,
        mapping,threshPct,
        {dist1:g1.dist,dist2:g2.dist,angle1:g1.angle,angle2:g2.angle,
         dihedral1:g1.dihedral,dihedral2:g2.dihedral});
}
// The whole decision, in one place so the panel and the verification harness
// read the same answer: cells first, then the base graphs, then the geometry.
function crystalStructureDiff(leftCry,leftInv,leftBonds,rightCry,rightInv,rightBonds,threshPct){
    var leftCell=crystalCellParams(leftCry);
    var rightCell=crystalCellParams(rightCry);
    if(!leftCell||!rightCell)return {status:'no-cell'};
    if(!crystalCellParamsEqual(leftCell,rightCell)){
        return {status:'cell-differs',leftCell:leftCell,rightCell:rightCell};
    }
    var mapping=findAtomMapping(leftCry.baseAtoms,leftBonds,rightCry.baseAtoms,rightBonds);
    if(!mapping){
        return {status:'skeleton-differs',leftCell:leftCell,rightCell:rightCell,
                leftAtoms:(leftCry.baseAtoms||[]).length,leftBonds:(leftBonds||[]).length,
                rightAtoms:(rightCry.baseAtoms||[]).length,rightBonds:(rightBonds||[]).length};
    }
    return {status:'ok',mapping:mapping,leftCell:leftCell,rightCell:rightCell,
            diffs:computeCrystalConformationDiffs(leftCry,leftInv,leftBonds,
                rightCry,rightInv,rightBonds,mapping,threshPct)};
}

// geo is optional: six closures that measure one side's internal coordinates.
// The molecular path leaves it out and gets the plain cartesian expressions it
// always used; the crystal comparison supplies a set that measures through each
// bond's periodic image, so the same loops, the same statistics and the same
// difference list serve both — and a crystal can never quietly use a different
// notion of "bond length" from a molecule.
function computeConformationDiffs(atoms1,bonds1,atoms2,bonds2,mapping,threshPct,geo){
    var diffs={bonds:[],angles:[],dihedrals:[]};
    var tp=threshPct||2;
    var statAcc={bonds:{},angles:{},dihedrals:{}};
    var dist1=geo?geo.dist1:function(i,j){return dist3d(atoms1[i],atoms1[j])};
    var dist2=geo?geo.dist2:function(i,j){return dist3d(atoms2[i],atoms2[j])};
    var angle1=geo?geo.angle1:function(i,j,k){return angle3d(atoms1[i],atoms1[j],atoms1[k])};
    var angle2=geo?geo.angle2:function(i,j,k){return angle3d(atoms2[i],atoms2[j],atoms2[k])};
    var dihedral1=geo?geo.dihedral1:function(a,b,c,d){return dihedral3d(atoms1[a],atoms1[b],atoms1[c],atoms1[d])};
    var dihedral2=geo?geo.dihedral2:function(a,b,c,d){return dihedral3d(atoms2[a],atoms2[b],atoms2[c],atoms2[d])};

    var bondSet2={};
    bonds2.forEach(function(b){bondSet2[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.order});
    bonds1.forEach(function(b){
        var a1m=Math.min(mapping[b.atom1],mapping[b.atom2]);
        var a2m=Math.max(mapping[b.atom1],mapping[b.atom2]);
        var key=a1m+'-'+a2m;
        if(bondSet2[key]!==undefined){
            var d1=dist1(b.atom1,b.atom2);
            var d2=dist2(mapping[b.atom1],mapping[b.atom2]);
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
                var a1=angle1(i1,j,i3);
                var a2=angle2(mapping[i1],mapping[j],mapping[i3]);
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
                    var d1=dihedral1(s0,start,mid,end);
                    var d2=dihedral2(mapping[s0],mapping[start],mapping[mid],mapping[end]);
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

// ===== Crystal Diff =====
//
// Two crystals are compared as cells (crystalStructureDiff): the screen's
// supercell and the file's wrapped cell are two drawings of the same object, so
// only the base cell is comparable. The right half is still drawn from
// diffData.atoms (the file's one cell) and the left from the displayed
// supercell; every difference index is a base index, and these helpers translate
// it to each half's display indices.
function crystalDiffDisplayMaps(baseMapping){
    var left=MD.atoms,right=diffData.atoms;
    var disp=new Array(left.length);
    var byBase=[];
    var rev=new Array(right.length).fill(-1);
    for(var i=0;i<left.length;i++){
        disp[i]=-1;
        var bi=left[i]&&left[i].baseIdx;
        if(bi==null||bi<0||bi>=baseMapping.length)continue;
        var rj=baseMapping[bi];
        if(rj==null||rj<0)continue;
        disp[i]=rj;
        if(!byBase[bi])byBase[bi]=[];
        byBase[bi].push(i);
        if(rev[rj]<0)rev[rj]=i;
    }
    var revBase=new Array(right.length).fill(-1);
    for(var b=0;b<baseMapping.length;b++){
        var rb=baseMapping[b];
        if(rb>=0&&rb<revBase.length)revBase[rb]=b;
    }
    diffMapping=disp;
    diffReverseMapping=rev;
    diffDisplayByBase=byBase;
    diffReverseBaseMapping=revBase;
}
function crystalDiffStatusHTML(res){
    var html='<span class="diff-close">×</span>';
    if(res.status==='cell-differs'){
        html+='<h4>Unit Cells Differ</h4>'+
            '<div class="diff-row">The two files state different cells, so their periodic graphs cannot be compared.</div>'+
            crystalCellDiffHTML(res.leftCell,res.rightCell)+
            '<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">The two structures are still shown side by side. Click Diff again to exit.</div>';
    }else if(res.status==='skeleton-differs'){
        html+='<h4>Skeletons Differ</h4>'+
            '<div class="diff-row">The unit cells match, but the base-cell graphs do not.</div>'+
            '<div class="diff-row">Left: '+res.leftAtoms+' atoms, '+res.leftBonds+' bonds (base cell)</div>'+
            '<div class="diff-row">Right: '+res.rightAtoms+' atoms, '+res.rightBonds+' bonds (base cell)</div>'+
            '<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Click Diff again to exit.</div>';
    }else{
        html+='<h4>Not Comparable</h4>'+
            '<div class="diff-row">A crystal is compared with another crystal: one side of this comparison states no unit cell.</div>'+
            '<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Click Diff again to exit.</div>';
    }
    return html;
}
function clearDiffMarks(){
    atomMeshes.forEach(function(m){m.userData.diffHi=false;setAtomHighlight(m,0)});
    diffAtomMeshes.forEach(function(m){m.userData.diffHi=false;setDiffAtomHighlight(m,0)});
    diffBondPieces.forEach(function(p){p.hi=false});
    diffBuildBondPools();
    needsRender=true;
}
function showCrystalDiffStatus(res){
    diffMapping=null;diffReverseMapping=null;diffDisplayByBase=null;diffReverseBaseMapping=null;
    showDiffPanel(crystalDiffStatusHTML(res),true);
    diffPanelEl.classList.add('show');
    var tag=res.status==='cell-differs'?'cell differs':
            (res.status==='skeleton-differs'?'skeleton differs':'not comparable');
    diffLabelRight.textContent='Diff: '+diffData.fileName+' ('+tag+')';
    modeInfoEl.textContent='Diff Mode ('+tag+')';
    clearDiffMarks();
}
function startCrystalDiff(msg){
    diffMode=true;
    diffLabelLeft.textContent='Original: '+(MD.title||'crystal');
    diffLabelLeft.classList.add('show');
    diffLabelRight.textContent='Diff: '+diffData.fileName;
    diffLabelRight.classList.add('show');
    var right=msg.crystal||null;
    var rightInv=right?crystalInverseOf(right.latticeVectors):null;
    if(!CRY||!right||!rightInv){
        diffBaseIndices=false;diffLabelAtoms=null;
        showCrystalDiffStatus({status:'incompatible'});
        updateBakeHint();
        enterDiffRender();
        vdSyncPanel();
        layoutPanels();
        return;
    }
    if(!CRY.bondsDetected)detectCrystalBaseBonds();
    var leftBonds=CRY.baseBonds;
    var rightBonds=detectCrystalBaseBondsFor(right,rightInv);
    diffRightCrystal=right;diffRightInv=rightInv;diffRightBonds=rightBonds;
    diffDetBonds1=leftBonds;diffDetBonds2=rightBonds;
    var res=crystalStructureDiff(CRY,CRY_INV,leftBonds,right,rightInv,rightBonds,diffThresholdPct);
    diffBaseIndices=true;
    diffLabelAtoms=CRY.baseAtoms;
    if(res.status==='ok'){
        crystalDiffDisplayMaps(res.mapping);
        var diffs=res.diffs;
        var total=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;
        showDiffPanel(buildDiffResultsHTML(diffs,total),true);
        wireDiffSlider();
        diffLabelRight.textContent='Diff: '+diffData.fileName;
        modeInfoEl.textContent='Diff Mode ('+total+' differences)';
        enterDiffRender(diffs,res.mapping);
        vdSyncPanel();
        layoutPanels();
        return;
    }
    showCrystalDiffStatus(res);
    updateBakeHint();
    enterDiffRender();
    vdSyncPanel();
    layoutPanels();
}
// A rebuild of the left structure (an edit, a style change, an undo) re-runs
// the cell comparison and updates the open panel in place; the right half is a
// model of the file's cell and does not move.
function crystalDiffRefresh(){
    if(!CRY||!diffRightCrystal||!diffRightInv)return;
    var res=crystalStructureDiff(CRY,CRY_INV,CRY.baseBonds,diffRightCrystal,diffRightInv,
        diffRightBonds,diffThresholdPct);
    if(res.status!=='ok'){showCrystalDiffStatus(res);return}
    diffBaseIndices=true;diffLabelAtoms=CRY.baseAtoms;
    crystalDiffDisplayMaps(res.mapping);
    var diffs=res.diffs;
    var total=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;
    var listEl=diffPanelEl.querySelector('#diff-results-list');
    if(listEl){
        listEl.innerHTML=buildDiffListHTML(diffs,total);
        wireDiffStatsToggle(listEl);
    }
    var fullHTML=buildDiffResultsHTML(diffs,total);
    diffPanelHTML=fullHTML;
    diffResultsHTML=fullHTML;
    modeInfoEl.textContent='Diff Mode ('+total+' differences)';
    diffLabelRight.textContent='Diff: '+diffData.fileName;
    updateDiffHighlights(diffs);
}
function startDiff(msg){
    // A powder pattern or a cavity report describes one cell, and the
    // comparison has two: the reading aids go away before the second structure
    // is drawn, so neither can cover its half of the view.
    closeAnalysisPanelsForDiff();
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
    diffBaseIndices=false;diffHasCrystal=false;diffLabelAtoms=null;diffDisplayByBase=null;
    diffReverseBaseMapping=null;diffRightCrystal=null;diffRightInv=null;diffRightBonds=null;

    // A crystal on either side is compared as a cell, never as the two atom
    // lists the screen and the file happen to hold.
    if(CRY||msg.crystal){
        diffHasCrystal=true;
        startCrystalDiff(msg);
        return;
    }

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
    // A crystal comparison's indices address the base cell, so its labels come
    // from that list; the molecular one keeps MD.atoms.
    var labelAtoms=diffLabelAtoms||MD.atoms;
    if(totalDiff===0){
        html+='<div class="diff-row" style="color:#4ec9b0">Structures are identical (same conformation within thresholds).</div>';
    }else{
        html+='<div class="diff-row">Mapping: '+labelAtoms.length+' atoms matched. '+totalDiff+' differences found.</div>';
    }
    html+=buildDiffStatsHTML(diffs.stats);
    if(totalDiff>0){
        if(diffs.bonds.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Length Differences ('+diffs.bonds.length+')</b></div>';
            diffs.bonds.sort(function(a,b){return b.pct-a.pct});
            diffs.bonds.forEach(function(d){
                html+='<div class="diff-row">  '+atomLabel(labelAtoms[d.i1],d.i1)+'-'+atomLabel(labelAtoms[d.i2],d.i2)+
                    ': '+d.d1.toFixed(3)+' vs '+d.d2.toFixed(3)+' Å (Δ='+d.diff.toFixed(3)+', '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.angles.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Angle Differences ('+diffs.angles.length+')</b></div>';
            diffs.angles.sort(function(a,b){return b.pct-a.pct});
            diffs.angles.forEach(function(d){
                html+='<div class="diff-row">  '+atomLabel(labelAtoms[d.i1],d.i1)+'-'+atomLabel(labelAtoms[d.i2],d.i2)+'-'+atomLabel(labelAtoms[d.i3],d.i3)+
                    ': '+d.a1.toFixed(1)+'° vs '+d.a2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+'°, '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.dihedrals.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Dihedral Differences ('+diffs.dihedrals.length+')</b></div>';
            diffs.dihedrals.sort(function(a,b){return b.pct-a.pct});
            diffs.dihedrals.forEach(function(d){
                html+='<div class="diff-row">  '+atomLabel(labelAtoms[d.i1],d.i1)+'-'+atomLabel(labelAtoms[d.i2],d.i2)+'-'+atomLabel(labelAtoms[d.i3],d.i3)+'-'+atomLabel(labelAtoms[d.i4],d.i4)+
                    ': '+d.d1.toFixed(1)+'° vs '+d.d2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+', '+d.pct.toFixed(2)+'%)</div>';
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
    // A crystal comparison is re-run on the cells, not on the two display lists
    // (whose atom counts only coincide at 0..1).
    if(diffHasCrystal){
        if(CRY&&diffRightCrystal)crystalDiffRefresh();
        return;
    }
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
    // A crystal comparison survives a rebuild of the left structure by
    // re-running the cell comparison; comparing display lists again would
    // reintroduce the supercell-against-cell mismatch this path removes.
    if(diffHasCrystal){
        if(CRY&&diffRightCrystal)crystalDiffRefresh();
        return;
    }
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
        if(diffBaseIndices&&diffDisplayByBase){
            var bi=MD.atoms[i]&&MD.atoms[i].baseIdx;
            m.userData.diffHi=(bi!=null&&highlightAtoms[bi]===true);
        }else m.userData.diffHi=highlightAtoms[i]?true:false;
    });

    diffAtomMeshes.forEach(function(m,i){
        var origI=(diffBaseIndices&&diffReverseBaseMapping)?
            diffReverseBaseMapping[i]:(diffReverseMapping?diffReverseMapping[i]:-1);
        if(origI==null)origI=-1;
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

// Build the comparison's right-hand model: one sphere per atom and two cylinder
// pieces per bond, drawn in the current display style exactly as the left-hand
// half is (that half is the ordinary scene, which rebuildScene rebuilds).
//
// hiAtom is keyed by right-hand atom index and hiBond by the pair of
// left-hand indices a bond maps to - the key this half has always used for its
// bond marks - so a caller can hand over either the sets the difference
// computation just produced or the ones already on screen.
//
// The model is built in place rather than by re-entering the comparison, so the
// two view states, the camera and the pan stay where the user put them.
function buildDiffModel(hiAtom,hiBond){
    if(!diffMolGroup||!diffData)return;
    // The atoms are records now, not children of the group: the instanced layer
    // below is what draws them, so there is nothing to remove or dispose here.
    diffAtomMeshes=[];
    var dSegs=sphereSegs(diffData.atoms.length);
    diffData.atoms.forEach(function(a,i){
        var gh=isGhostAtom(a);
        // The same radius rule as the left half — the default balls, van der
        // Waals spheres, or a joint as thick as the sticks. A comparison's atoms
        // carry no script scale, so there is none to apply.
        var r=drawnRadius({element:a.element,ghost:gh},0);
        var isHi=hiAtom[i]===true;
        // A ghost is the small faint marker it is on the left, so the two halves
        // read the same way. The mark is a state on the handle; the instanced
        // layer draws a marked atom in the comparison's orange.
        var diffCol=gh?GHOST_COLOR:a.color;
        var mesh=atomHandle();
        mesh.position.set(a.x-diffCX,a.y-diffCY,a.z-diffCZ);
        mesh.userData={element:a.element,index:i,diffHi:isHi,baseColor:diffCol,ghost:gh,
                       radius:r,drawn:r>0.0005,hlState:isHi?2:0};
        mesh.visible=r>0.0005;
        diffAtomMeshes.push(mesh);
    });
    sphereLayerBuild(diffSphereLayer,diffAtomMeshes,diffMolGroup,dSegs,true);

    // The right-hand bonds: two pieces per bond (one per half), each carrying
    // the transform it is drawn with. The comparison's structure never moves,
    // so a piece stores its own numbers and the highlight can be re-applied by
    // refilling the two pools rather than editing a material per mesh. A
    // space-filling model draws no sticks at all; this half has its own piece
    // list, so it takes the same decision as createBond rather than inheriting
    // it.
    diffBondPieces.length=0;
    if(displayStyle!=='spacefill')diffData.bonds.forEach(function(b){
        var a1=diffData.atoms[b.atom1],a2=diffData.atoms[b.atom2];
        if(!a1||!a2)return;
        var sx=a1.x-diffCX,sy=a1.y-diffCY,sz=a1.z-diffCZ;
        var ex=a2.x-diffCX,ey=a2.y-diffCY,ez=a2.z-diffCZ;
        var dx=ex-sx,dy=ey-sy,dz=ez-sz;
        var l=Math.sqrt(dx*dx+dy*dy+dz*dz);
        if(l<0.001)return;
        var mpx=(sx+ex)*0.5,mpy=(sy+ey)*0.5,mpz=(sz+ez)*0.5;
        var origI1=diffReverseMapping?diffReverseMapping[b.atom1]:-1;
        var origI2=diffReverseMapping?diffReverseMapping[b.atom2]:-1;
        var bondKey=Math.min(origI1,origI2)+'-'+Math.max(origI1,origI2);
        var isHi=hiBond[bondKey]===true;
        addDiffBondPiece(sx,sy,sz,mpx-sx,mpy-sy,mpz-sz,l/2,BOND_R,new THREE.Color(a1.color),bondKey,isHi);
        addDiffBondPiece(mpx,mpy,mpz,ex-mpx,ey-mpy,ez-mpz,l/2,BOND_R,new THREE.Color(a2.color),bondKey,isHi);
    });
    diffBuildBondPools();
    buildDiffPolyhedra();
}

// The right half wears the same display style as the left: the coordination
// hulls of its metal/metalloid centres over its ball-and-stick picture. The
// comparison's bonds carry no periodic image of their own (that half is built
// against transformed coordinates), so there is no lattice offset to resolve.
function buildDiffPolyhedra(){
    disposePolyMesh(diffPolyFacesMesh);
    disposePolyMesh(diffPolyEdgesMesh);
    diffPolyFacesMesh=null;diffPolyEdgesMesh=null;
    if(displayStyle!=='polyhedra'||!diffData||!diffMolGroup)return;
    var body=polyhedraGeometry(diffData.atoms,diffData.bonds,diffCX,diffCY,diffCZ,null);
    if(!body)return;
    var pair=polyMeshPair(body,diffMolGroup);
    diffPolyFacesMesh=pair[0];diffPolyEdgesMesh=pair[1];
}

// Rebuild that model from the state on screen, for a change of drawing style:
// the marks the difference computation put on the two halves are read back off
// the meshes and pieces before they are replaced, so restyling a comparison
// keeps its differences marked and keeps the right-hand selection glowing.
function rebuildDiffModel(){
    if(!diffMode||!diffMolGroup)return;
    var hiAtom={},hiBond={};
    diffAtomMeshes.forEach(function(m,i){if(m.userData.diffHi)hiAtom[i]=true});
    diffBondPieces.forEach(function(p){if(p.hi)hiBond[p.key]=true});
    buildDiffModel(hiAtom,hiBond);
    highlightSelected();
    needsRender=true;
}

function enterDiffRender(diffs,mapping){
    if(vibActive||vibPaused)stopVibration();
    stopMoveSession();
    // The contact search belongs to one structure; while two are on screen the
    // panel says so and its lines leave the picture.
    if(neighActive)neighRefresh();
    // The interaction list is a reading of the same geometry: an edit, an undo
    // or a switch re-detects it against the atoms that are on screen now.
    if(intActive)intRefresh();
    // The lattice plane is one cell's reading aid and has no side in a
    // comparison, so it leaves the picture with the second structure.
    hklRefresh();
    // A comparison has two structures and two view states, so the panel stops
    // advertising Alt coordinate editing for as long as it is up.
    updateBakeHint();
    diffPivot=new THREE.Group();scene.add(diffPivot);
    diffMolGroup=new THREE.Group();diffPivot.add(diffMolGroup);

    diffCX=0;diffCY=0;diffCZ=0;
    diffData.atoms.forEach(function(a){diffCX+=a.x;diffCY+=a.y;diffCZ+=a.z});
    if(diffData.atoms.length>0){diffCX/=diffData.atoms.length;diffCY/=diffData.atoms.length;diffCZ/=diffData.atoms.length}

    // The right-hand model is built from the marks the difference computation
    // produced, keyed the way that half is addressed (see buildDiffModel).
    var hiAtom={},hiBond={};
    if(diffs&&mapping){
        var markHi=function(idx){var ri=mapping[idx];if(ri>=0)hiAtom[ri]=true};
        diffs.bonds.forEach(function(d){
            markHi(d.i1);markHi(d.i2);
            hiBond[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true;
        });
        diffs.angles.forEach(function(d){markHi(d.i1);markHi(d.i2);markHi(d.i3)});
        diffs.dihedrals.forEach(function(d){markHi(d.i1);markHi(d.i2);markHi(d.i3);markHi(d.i4)});
    }
    buildDiffModel(hiAtom,hiBond);

    if(diffs&&mapping){
        var leftHi={};
        // A crystal difference names a base site; every displayed image of that
        // site wears the mark, because the user sees a supercell.
        var markLeft=function(idx){
            if(diffBaseIndices&&diffDisplayByBase){
                var list=diffDisplayByBase[idx];
                if(list)for(var k=0;k<list.length;k++)leftHi[list[k]]=true;
            }else leftHi[idx]=true;
        };
        diffs.bonds.forEach(function(d){markLeft(d.i1);markLeft(d.i2)});
        diffs.angles.forEach(function(d){markLeft(d.i1);markLeft(d.i2);markLeft(d.i3)});
        diffs.dihedrals.forEach(function(d){markLeft(d.i1);markLeft(d.i2);markLeft(d.i3);markLeft(d.i4)});
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
    sphereLayerClear(diffSphereLayer);
    diffAtomMeshes=[];
    // The hull meshes left with the pivot above; the references must not
    // outlive them, or a later style change would dispose them twice.
    diffPolyFacesMesh=null;diffPolyEdgesMesh=null;
    diffSelectedAtoms=[];
    diffMapping=null;diffReverseMapping=null;
    diffBaseIndices=false;diffHasCrystal=false;diffLabelAtoms=null;
    diffDisplayByBase=null;diffReverseBaseMapping=null;
    diffRightCrystal=null;diffRightInv=null;diffRightBonds=null;
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
    // One structure again: the contact list and its lines come back with it.
    if(neighActive)neighRefresh();
    // The interaction list is a reading of the same geometry: an edit, an undo
    // or a switch re-detects it against the atoms that are on screen now.
    if(intActive)intRefresh();
    // ...and so does the lattice plane, against the one cell left on screen.
    hklRefresh();
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
    sphereLayerSyncPositions(mainSphereLayer);
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
    // A displaced geometry is a different geometry: the aromatic sides are
    // re-planned from it (a displaced ring is still a ring, and its plane may
    // have tilted), which the pieces then read when they are placed. Playback
    // of a mode deliberately keeps the equilibrium plan: it borrows the
    // geometry for one cycle and puts it back, and a plan that changed 60 times
    // a second could flip a side mid-animation.
    planAromaticDashes();
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
    // atomLabel marks a ghost with a star (C*3, X5), so a selection that is
    // about to be ghosted - or one that already is - reads as points, not atoms.
    var list=indices.slice(0,15).map(function(i){return atomLabel(atoms[i],i)});
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
    function projectGroup(group,meshes,camW,offX){
        if(!group||!group.children.length)return;
        var box=new THREE.Box3().setFromObject(group);
        if(box.isEmpty())return;
        // The spheres are instanced from a unit geometry, so their own object
        // bounds say "one unit around the group origin". Each atom's drawn
        // radius is added here instead, which is exactly the footprint the old
        // per-atom Mesh gave Box3.setFromObject.
        if(meshes){
            var pv=new THREE.Vector3();
            for(var mi=0;mi<meshes.length;mi++){
                var mh=meshes[mi];
                if(!mh)continue;
                var mr=mh.userData?mh.userData.radius:0;
                pv.set(mh.position.x,mh.position.y,mh.position.z).applyMatrix4(group.matrixWorld);
                if(pv.x-mr<box.min.x)box.min.x=pv.x-mr;
                if(pv.y-mr<box.min.y)box.min.y=pv.y-mr;
                if(pv.z-mr<box.min.z)box.min.z=pv.z-mr;
                if(pv.x+mr>box.max.x)box.max.x=pv.x+mr;
                if(pv.y+mr>box.max.y)box.max.y=pv.y+mr;
                if(pv.z+mr>box.max.z)box.max.z=pv.z+mr;
            }
        }
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
            projectGroup(moleculeGroup,atomMeshes,halfW,0);
            camera.aspect=(w-halfW)/h;camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);camera.updateMatrixWorld(true);
            projectGroup(diffMolGroup,diffAtomMeshes,(w-halfW),halfW);
        }else{
            camera.aspect=w/h;camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);camera.updateMatrixWorld(true);
            projectGroup(moleculeGroup,atomMeshes,w,0);
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
// One place where a panel and its chip cannot both be visible. The × of a
// panel sets its flag, the chip clears it, and the mode (or the data) decides
// whether the pair applies at all. Rebuilt panels (Molecular Mechanics, eggs,
// Diff, convergence, modes, route) refresh their own chip in the code that
// builds them; the static panels call this helper after any state change.
function panelChip(id,on){
    var el=document.getElementById(id);
    if(el)el.classList.toggle('show',!!on);
}
function panelMinimize(name){
    panelMinState[name]=true;
    syncPanelChips();
    layoutPanels();
}
function panelRestore(name){
    panelMinState[name]=false;
    syncPanelChips();
    layoutPanels();
}
// The powder pattern, the cavity report and the BVS table each read one unit
// cell; a comparison shows two structures, so starting Diff folds them away
// rather than leaving a reading aid over the right-hand half.
function closeAnalysisPanelsForDiff(){
    if(xrdOpen)xrdClosePanel();
    if(cavOpen)cavityClosePanel();
}
function syncPanelChips(){
    var st=panelMinState;
    var cp=document.getElementById('crystal-panel');
    if(cp)cp.style.display=(CRY&&!st.crystal)?'block':'none';
    panelChip('crystal-reopen',!!CRY&&st.crystal);
    // The egg list only exists while the molecule draws patches: without rows
    // there is neither a panel nor a chip to fold it into.
    var eggs=eggSegments().length;
    var ep=document.getElementById('egg-panel');
    if(ep){
        if(!eggs)ep.classList.remove('show');
        else ep.classList.toggle('show',!st.egg);
    }
    panelChip('egg-reopen',eggs>0&&st.egg);
    var inView=currentMode==='view';
    var vp=document.getElementById('view-panel');
    if(vp)vp.style.display=(inView&&!st.view)?'block':'none';
    panelChip('view-reopen',inView&&st.view);
    var bp=document.getElementById('bake-panel');
    if(bp)bp.style.display=(inView&&!st.bake)?'block':'none';
    panelChip('bake-reopen',inView&&st.bake);
    var sp=document.getElementById('select-panel');
    if(sp)sp.style.display=(currentMode==='selectAtoms'&&!st.select)?'block':'none';
    panelChip('select-reopen',currentMode==='selectAtoms'&&st.select);
    var np=document.getElementById('neigh-panel');
    if(np)np.style.display=(currentMode==='neighbours'&&!st.neigh)?'block':'none';
    panelChip('neigh-reopen',currentMode==='neighbours'&&st.neigh);
    var ip=document.getElementById('int-panel');
    if(ip)ip.style.display=(currentMode==='interactions'&&!st.int)?'block':'none';
    panelChip('int-reopen',currentMode==='interactions'&&st.int);
    var rp=document.getElementById('rotate-panel');
    if(rp)rp.style.display=(currentMode==='rotateGroup'&&!st.rotate)?'block':'none';
    panelChip('rotate-reopen',currentMode==='rotateGroup'&&st.rotate);
    var xp=document.getElementById('xrd-panel');
    if(xp)xp.style.display=(xrdOpen&&!st.xrd)?'block':'none';
    panelChip('xrd-reopen',!!(xrdOpen&&CRY&&st.xrd));
    // The cavity report is opened from the crystal panel like the powder
    // pattern, and folds into its own chip the same way.
    var cvp=document.getElementById('cav-panel');
    if(cvp)cvp.style.display=(cavOpen&&!st.cav)?'block':'none';
    panelChip('cav-reopen',!!(cavOpen&&CRY&&st.cav));
}
// Turning a floating panel into a chip and back is one click each way for the
// static panels; the rebuilt panels bind their own in the code that builds
// them. Bound once, by id, so a rebuilt innerHTML cannot leave a dead handler.
(function(){
    function on(id,fn){var el=document.getElementById(id);if(el)el.addEventListener('click',fn)}
    on('crystal-close',function(){panelMinimize('crystal')});
    on('crystal-reopen',function(){panelRestore('crystal')});
    on('view-close',function(){panelMinimize('view')});
    on('view-reopen',function(){panelRestore('view')});
    on('bake-close',function(){panelMinimize('bake')});
    on('bake-reopen',function(){panelRestore('bake')});
    on('select-close',function(){panelMinimize('select')});
    on('select-reopen',function(){panelRestore('select')});
    on('neigh-close',function(){panelMinimize('neigh')});
    on('neigh-reopen',function(){panelRestore('neigh')});
    on('int-close',function(){panelMinimize('int')});
    on('int-reopen',function(){panelRestore('int')});
    on('rotate-close',function(){panelMinimize('rotate')});
    on('rotate-reopen',function(){panelRestore('rotate')});
    on('xrd-reopen',function(){xrdOpenPanel()});
    on('cav-reopen',function(){cavityOpenPanel()});
    on('egg-reopen',function(){openEggPanel()});
})();
/**
 * Place every visible panel without letting two of them (or a panel and a
 * chip) overlap.
 *
 * The panels are visited in priority order and each takes the first free slot
 * on its own edge, nearest that edge first, scanning in one-pixel steps. Before
 * anything is placed, each side of the container is given a width budget (a
 * wider panel is narrowed and wraps or scrolls its content) and each stack a
 * height budget (scrollable panels shrink proportionally, down to a floor), so
 * a tall panel cannot use up the strip its neighbours need. A panel whose
 * natural height still does not fit is capped to the largest free strip at
 * that edge and scrolls inside itself, which is what turns "they cannot all be
 * full size" into "they are all readable", not into "they cover each other".
 * The molecule's screen box is honoured by the two Diff labels, which are
 * moved along the top edge until they clear it; the corner panels keep their
 * corner, as they always have.
 *
 * Every measurement is container-relative: a slot is converted to inline
 * left/top, so a panel cannot drift with a CSS anchor after a resize.
 */
function layoutPanels(){
    // The Molecular Mechanics reopen button follows the mode (and the run) the
    // way the other reopen buttons follow their panels, so it is refreshed
    // before the layout measures anything.
    mmUpdateReopen();
    syncPanelChips();
    var container=document.getElementById('container');
    if(!container){mmDrawCurves();return;}
    var cRect=container.getBoundingClientRect();
    var cW=cRect.width,cH=cRect.height;
    if(cW<2||cH<2){mmDrawCurves();return;}
    var margin=8,gap=2,minH=16;
    var molBox=getMoleculeScreenBox();
    var molRect=molBox?{l:molBox.left,t:molBox.top,r:molBox.right,b:molBox.bottom}:null;
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
        {id:'diff-reopen',pr:2,anc:'br',chip:true},
        {id:'crystal-reopen',pr:2,anc:'tr',chip:true},
        {id:'egg-reopen',pr:2,anc:'tr',chip:true},
        // Auto Rotate and Edit Coordinates are two cards, each with its own
        // collision box, stacked bottom-up in that reading order.
        {id:'bake-panel',pr:3,anc:'br'},
        {id:'view-panel',pr:4,anc:'br'},
        {id:'select-panel',pr:3,anc:'bl'},
        {id:'neigh-panel',pr:3,anc:'bl'},
        {id:'int-panel',pr:3,anc:'bl'},
        {id:'rotate-panel',pr:3,anc:'bl'},
        {id:'select-reopen',pr:3,anc:'bl',chip:true},
        {id:'neigh-reopen',pr:3,anc:'bl',chip:true},
        {id:'int-reopen',pr:3,anc:'bl',chip:true},
        {id:'rotate-reopen',pr:3,anc:'bl',chip:true},
        // The powder panel is a reading aid opened from the crystal panel: it
        // is stacked above the bottom-left tools, so it never covers them.
        {id:'xrd-panel',pr:4,anc:'bl'},
        {id:'xrd-reopen',pr:4,anc:'bl',chip:true},
        // The cavity report is a second reading aid opened from the crystal
        // panel; it takes the same bottom-left stack as the powder one.
        {id:'cav-panel',pr:4,anc:'bl'},
        {id:'cav-reopen',pr:4,anc:'bl',chip:true},
        {id:'diff-label',pr:4,anc:'tc',pct:0.25},
        {id:'diff-label-right',pr:5,anc:'tc',pct:0.75},
        {id:'bake-reopen',pr:3,anc:'br',chip:true},
        {id:'view-reopen',pr:4,anc:'br',chip:true},
        {id:'opt-panel',pr:7,anc:'tl'},
        {id:'freq-panel',pr:8,anc:'tl'},
        {id:'route-panel',pr:8,anc:'tl'},
        {id:'mm-reopen',pr:9,anc:'tl',chip:true},
        {id:'opt-reopen',pr:9,anc:'tl',chip:true},
        {id:'freq-reopen',pr:10,anc:'tl',chip:true},
        {id:'route-reopen',pr:10,anc:'tl',chip:true}
    ];
    var vis=[];
    defs.forEach(function(d){
        var el=document.getElementById(d.id);
        if(!el)return;
        // A cap from the previous layout belongs to the layout, not to the
        // panel (the Molecular Mechanics panel keeps a height the user
        // dragged): the original value goes back before measuring.
        if(el.__lpMaxH!==undefined){el.style.maxHeight=el.__lpMaxH;el.__lpMaxH=undefined;}
        el.style.top='';el.style.bottom='';el.style.left='';el.style.right='';
        el.style.transform='';el.style.maxWidth='';el.style.fontSize='';
        if(window.getComputedStyle(el).display==='none')return;
        vis.push({d:d,el:el,sz:measure(el)});
    });
    vis.sort(function(a,b){return a.d.pr-b.d.pr;});
    var placed=[];
    function measure(el){var r=el.getBoundingClientRect();return {w:r.width,h:r.height};}
    // A gap of 2 px between boxes, so two panels that touch are still two
    // panels and a chip can be clicked without catching its neighbour.
    function hit(r,q){return r.l<q.r+gap&&r.r>q.l-gap&&r.t<q.b+gap&&r.b>q.t-gap;}
    function hitsPanel(r){for(var i=0;i<placed.length;i++){if(hit(r,placed[i]))return true}return false;}
    function hitsMol(r){return !!(molRect&&r.l<molRect.r&&r.r>molRect.l&&r.t<molRect.b&&r.b>molRect.t);}
    function capHeight(el,h){
        if(el.__lpMaxH===undefined)el.__lpMaxH=el.style.maxHeight;
        el.style.maxHeight=Math.max(20,Math.floor(h))+'px';
    }
    function restoreCap(el){
        if(el.__lpMaxH!==undefined){el.style.maxHeight=el.__lpMaxH;el.__lpMaxH=undefined;}
    }
    function place(el,x,y){
        el.style.left=Math.round(x)+'px';
        el.style.top=Math.round(y)+'px';
        el.style.right='auto';
        el.style.bottom='auto';
        el.style.transform='none';
    }
    function xFor(anc,pct,w){
        if(pct!=null)return Math.max(margin,Math.min(cW-margin-w,cW*pct-w/2));
        if(anc==='tr'||anc==='br')return Math.max(margin,cW-margin-w);
        if(anc==='tc')return Math.max(margin,cW/2-w/2);
        return margin;
    }
    // Boxes that hug the same edge share that edge's strip, so the strip is
    // budgeted before anything is placed: the widest boxes on a side are
    // capped to their share of the container width, and the tallest
    // scrollable one gives up height until the stack fits. Without this a
    // tall panel in a narrow window leaves no strip for the cards below it,
    // and those cards are pushed to the far side, on top of the other column.
    var sides={l:[],r:[]};
    vis.forEach(function(v){
        if(v.d.anc==='tl'||v.d.anc==='bl')sides.l.push(v);
        else if(v.d.anc==='tr'||v.d.anc==='br')sides.r.push(v);
    });
    function sideW(list){
        var w=0;
        list.forEach(function(v){if(v.sz.w>w)w=v.sz.w;});
        return w;
    }
    var availW=cW-2*margin-gap,fairW=availW/2;
    var lw=sideW(sides.l),rw=sideW(sides.r);
    var capL=lw,capR=rw;
    if(lw+rw>availW){
        // The narrower side keeps its width; the wider one takes the rest,
        // but never less than an even split when both are wide.
        if(lw>=rw){capL=Math.max(availW-rw,fairW);capR=Math.min(rw,availW-capL);}
        else{capR=Math.max(availW-lw,fairW);capL=Math.min(lw,availW-capR);}
    }
    function capWidth(list,w){
        list.forEach(function(v){
            if(v.sz.w<=w+0.5)return;
            v.el.style.maxWidth=Math.floor(w)+'px';
            v.sz=measure(v.el);
        });
    }
    capWidth(sides.l,capL);
    capWidth(sides.r,capR);
    // The information card and the Diff labels are one-shot text: clipping
    // them would hide words, so they are never height-capped. The other
    // panels scroll their own content.
    function canShrink(v){
        if(v.d.chip)return false;
        if(v.d.id==='mol-info'||v.d.id==='diff-label'||v.d.id==='diff-label-right')return false;
        var oy=window.getComputedStyle(v.el).overflowY;
        return oy==='auto'||oy==='scroll'||oy==='hidden';
    }
    function budgetH(list){
        var availH=cH-2*margin,pass=0;
        function sum(){
            var t=0;
            for(var k=0;k<list.length;k++)t+=list[k].sz.h+(k?gap:0);
            return t;
        }
        while(sum()>availH+0.5&&pass++<8){
            var fixedH=(list.length-1)*gap,shrinkH=0;
            for(var k=0;k<list.length;k++){
                if(canShrink(list[k])&&list[k].sz.h>40)shrinkH+=list[k].sz.h;
                else fixedH+=list[k].sz.h;
            }
            var room=availH-fixedH;
            if(shrinkH<=0||room<=0)break;
            var f=room/shrinkH,changed=false;
            for(var k2=0;k2<list.length;k2++){
                var v=list[k2];
                if(!canShrink(v)||v.sz.h<=40)continue;
                var target=Math.max(40,v.sz.h*f);
                if(target>=v.sz.h-0.5||v.__lpStuck)continue;
                var before=v.sz.h;
                capHeight(v.el,target);
                v.sz=measure(v.el);
                if(v.sz.h<before-0.5)changed=true;
                else v.__lpStuck=true;
            }
            if(!changed)break;
        }
    }
    budgetH(sides.l);
    budgetH(sides.r);
    // Own edge first, then the opposite edge on the same side, then the other
    // side: a panel only leaves its corner when that corner has no free strip
    // left at all.
    function anchorOrder(anc){
        var order=[anc];
        var opp={tl:'bl',bl:'tl',tr:'br',br:'tr',tc:'tl'}[anc];
        if(opp&&order.indexOf(opp)<0)order.push(opp);
        var mir={tl:'tr',tr:'tl',bl:'br',br:'bl',tc:'tr'}[anc];
        if(mir&&order.indexOf(mir)<0)order.push(mir);
        ['tl','tr','bl','br'].forEach(function(a){if(order.indexOf(a)<0)order.push(a)});
        return order;
    }
    // The first free slot along one edge, nearest that edge first. A placement
    // label (a Diff label) also looks for a row clear of the molecule, which is
    // what keeps the labels off the structure when there is room.
    function findSlot(d,anc,w,h,avoidMol){
        var x=xFor(anc,anc===d.anc?d.pct:null,w);
        var last=cH-margin-h;
        if(last<margin)return null;
        var bottom=(anc==='bl'||anc==='br');
        var first=null;
        // One pixel per step: two stacked cards can be left with only a
        // couple of pixels of clearance, and a coarse grid stepped straight
        // past the one row that fits.
        var step=1;
        for(var y=bottom?Math.floor(last):margin;
                bottom?y>=margin-0.01:y<=last+0.01;
                y+=bottom?-step:step){
            var r={l:x,t:y,r:x+w,b:y+h};
            if(r.l<margin-0.5||r.r>cW-margin+0.5)continue;
            if(hitsPanel(r))continue;
            if(!first)first=r;
            if(!avoidMol||!hitsMol(r))return r;
        }
        return first;
    }
    // Free vertical strips in one x column, excluding the panels already
    // placed there (or across it). Used when the natural height does not fit.
    function freeGaps(anc,w){
        var x=xFor(anc,null,w);
        var bars=[];
        for(var i=0;i<placed.length;i++){
            var q=placed[i];
            if(q.r+gap<=x||q.l-gap>=x+w)continue;
            bars.push(q);
        }
        bars.sort(function(p,q){return p.t-q.t});
        var gaps=[],top=margin;
        for(var b=0;b<bars.length;b++){
            if(bars[b].t-gap>top)gaps.push({t:top,b:bars[b].t-gap});
            if(bars[b].b+gap>top)top=bars[b].b+gap;
        }
        if(cH-margin>top)gaps.push({t:top,b:cH-margin});
        return gaps;
    }
    // Cap a panel to a free strip rather than let it overlap: its own CSS (or
    // the overflow rule added for these cards) makes the rest scrollable. The
    // strip nearest the panel's edge wins.
    function findShrunkSlot(el,d,anc,sz){
        if(d.chip)return null;
        var w=Math.min(sz.w,cW-2*margin);
        var gaps=freeGaps(anc,w);
        var bottom=(anc==='bl'||anc==='br');
        gaps.sort(function(p,q){
            var dp=Math.abs(bottom?(cH-margin-p.b):(p.t-margin));
            var dq=Math.abs(bottom?(cH-margin-q.b):(q.t-margin));
            return dp-dq;
        });
        for(var i=0;i<gaps.length;i++){
            var g=gaps[i],room=g.b-g.t;
            if(room<minH)continue;
            capHeight(el,Math.min(sz.h,room));
            var sz2=measure(el);
            var x=xFor(anc,null,sz2.w);
            var y=bottom?g.b-sz2.h:g.t;
            var r={l:x,t:y,r:x+sz2.w,b:y+sz2.h};
            // Capping can add a scrollbar and widen the box; the fresh box is
            // what has to fit, not the one the gaps were measured from.
            if(r.l<margin-0.5||r.r>cW-margin+0.5)continue;
            if(hitsPanel(r))continue;
            return r;
        }
        // The cap belonged to this attempt only; the caller may still find a
        // full-height slot under another anchor or width.
        restoreCap(el);
        return null;
    }
    for(var i=0;i<vis.length;i++){
        var el=vis[i].el,d=vis[i].d;
        var sz=measure(el);
        if(sz.h>cH-2*margin){capHeight(el,cH-2*margin);sz=measure(el);}
        if(sz.w>cW-2*margin){el.style.maxWidth=Math.floor(cW-2*margin)+'px';sz=measure(el);}
        var slot=null;
        var order=anchorOrder(d.anc);
        for(var a=0;a<order.length&&!slot;a++){
            slot=findSlot(d,order[a],sz.w,sz.h,order[a]==='tc');
        }
        if(!slot){
            for(var a2=0;a2<order.length&&!slot;a2++){
                slot=findShrunkSlot(el,d,order[a2],sz);
            }
        }
        if(!slot&&!d.chip){
            // Neither corner nor a shortened strip fit the natural width: a
            // narrower box can slip beside a tall neighbour. Widths come
            // from the edges of what is already placed, widest first, so the
            // least content is hidden.
            var ws=[];
            for(var pi=0;pi<placed.length;pi++){
                var pw=[placed[pi].l-margin-gap,cW-margin-placed[pi].r-gap];
                for(var pj=0;pj<2;pj++){
                    var w0=Math.floor(pw[pj]);
                    if(w0>=120&&w0<sz.w-0.5&&ws.indexOf(w0)<0)ws.push(w0);
                }
            }
            ws.sort(function(a,b){return b-a});
            for(var wi=0;wi<ws.length&&!slot;wi++){
                restoreCap(el);
                el.style.maxWidth=ws[wi]+'px';
                var sw=measure(el);
                for(var wa=0;wa<order.length&&!slot;wa++){
                    slot=findSlot(d,order[wa],sw.w,sw.h,order[wa]==='tc');
                }
                for(var wb=0;wb<order.length&&!slot;wb++){
                    slot=findShrunkSlot(el,d,order[wb],sw);
                }
                if(!slot)el.style.maxWidth='';
            }
        }
        if(!slot){
            // A viewport too small for every box at once: the panel keeps its
            // own corner (undoing any cap) rather than disappearing. Each other
            // panel has already been given a free slot around it, so only a
            // window smaller than the panel itself can leave a residue.
            restoreCap(el);
            el.style.maxWidth='';
            sz=measure(el);
            if(sz.h>cH-2*margin){capHeight(el,cH-2*margin);sz=measure(el);}
            if(sz.w>cW-2*margin){el.style.maxWidth=Math.floor(cW-2*margin)+'px';sz=measure(el);}
            var bw=sz.w,bh=sz.h;
            var bx=xFor(d.anc,d.pct,bw);
            var by=(d.anc==='bl'||d.anc==='br')?cH-margin-bh:margin;
            slot={l:bx,t:by,r:bx+bw,b:by+bh};
        }
        place(el,slot.l,slot.t);
        var fr=el.getBoundingClientRect();
        placed.push({l:fr.left-cRect.left,t:fr.top-cRect.top,
                     r:fr.right-cRect.left,b:fr.bottom-cRect.top});
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
    if(currentMode==='neighbours'){
        // One atom is the centre of the search, not a member of a growing
        // selection: a click moves the anchor and the panel follows it.
        selectedAtoms=[idx];diffSelectedAtoms=[];originalCoords=null;
        highlightSelected();
        selInfoEl.textContent='';
        neighSetAnchor(idx);
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

// ===== Neighbours / contact search =====
// The panel answers "what is around this atom?" with one list, nearest first.
// A crystal search runs in fractional space: every base site is placed at each
// integer image that could reach inside the radius, so a contact across a cell
// face is found even when that image is not drawn. A molecule is just the atoms
// on screen. Bonds are hidden by default — the point of the search is the
// non-covalent shell around an atom — and the checkbox brings them back.
var neighActive=false;
var neighAnchor=null;
var neighRadius=3.5;
var neighSkipBonded=true;
var neighRows=[];
var neighLines=null;
var NEIGH_MAX_ROWS=200;
var NEIGH_MAX_RADIUS=12;

function neighListEl2(){return document.getElementById('ng-list')}
function neighSummaryEl(){return document.getElementById('ng-summary')}

// A site the drawing hides (the other alternative of a disorder group) has no
// image on screen, so it is not a contact the reader can see either.
function neighBaseVisible(ba){
    return !!ba&&(!CRY||!CRY.disorderFilter||!ba.disorderGroup||ba.disorderGroup===CRY.disorderFilter);
}

function neighFractional(p){
    var inv=CRY_INV;
    return [inv[0][0]*p.x+inv[0][1]*p.y+inv[0][2]*p.z,
            inv[1][0]*p.x+inv[1][1]*p.y+inv[1][2]*p.z,
            inv[2][0]*p.x+inv[2][1]*p.y+inv[2][2]*p.z];
}

// Both directions of every base bond: a bond (a, b, shift s) is the contact
// between image c of a and image c+s of b, and asking from b it is s on the
// other side. Exact integer keys, so no tolerance is involved.
function neighCrystalBondSet(){
    var set={};
    var list=(CRY&&CRY.baseBonds)||[];
    for(var i=0;i<list.length;i++){
        var b=list[i],s=b.shift||ZERO_SHIFT;
        set[b.atom1+'_'+b.atom2+'_'+s[0]+','+s[1]+','+s[2]]=true;
        set[b.atom2+'_'+b.atom1+'_'+(-s[0])+','+(-s[1])+','+(-s[2])]=true;
    }
    return set;
}

function neighDisplayMap(){
    var map={};
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a.baseIdx!=null&&a.imgShift)map[a.baseIdx+'_'+a.imgShift.join(',')]=i;
    }
    return map;
}

function neighAnchorDisplayIndex(){
    if(!neighAnchor)return -1;
    if(neighAnchor.bi==null)return MD.atoms[neighAnchor.idx]?neighAnchor.idx:-1;
    var s=neighAnchor.shift;
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a.baseIdx===neighAnchor.bi&&a.imgShift&&
           a.imgShift[0]===s[0]&&a.imgShift[1]===s[1]&&a.imgShift[2]===s[2])return i;
    }
    return -1;
}

// The anchor's own Cartesian position, recomputed from the base site so it is
// exact even if the image carrying it has just left the screen.
function neighAnchorPoint(){
    if(!neighAnchor)return null;
    if(neighAnchor.bi!=null&&CRY&&CRY.baseAtoms&&CRY.latticeVectors){
        var ba=CRY.baseAtoms[neighAnchor.bi],lv=CRY.latticeVectors,s=neighAnchor.shift;
        if(!ba)return null;
        return {x:ba.x+s[0]*lv[0][0]+s[1]*lv[1][0]+s[2]*lv[2][0],
                y:ba.y+s[0]*lv[0][1]+s[1]*lv[1][1]+s[2]*lv[2][1],
                z:ba.z+s[0]*lv[0][2]+s[1]*lv[1][2]+s[2]*lv[2][2]};
    }
    var a=MD.atoms[neighAnchor.idx];
    return a?{x:a.x,y:a.y,z:a.z}:null;
}

function neighByDistance(a,b){
    if(Math.abs(a.dist-b.dist)>1e-12)return a.dist-b.dist;
    if(a.bi!=null&&b.bi!=null){
        if(a.bi!==b.bi)return a.bi-b.bi;
        return (a.s[0]-b.s[0])||(a.s[1]-b.s[1])||(a.s[2]-b.s[2]);
    }
    return (a.idx==null?0:a.idx)-(b.idx==null?0:b.idx);
}

function neighComputeMolecule(){
    var out={rows:[],total:0,bonded:0};
    var ai=neighAnchor.idx,a=MD.atoms[ai];
    if(!a)return out;
    var r2=neighRadius*neighRadius,bonded={},i;
    for(i=0;i<MD.bonds.length;i++){
        var b=MD.bonds[i];
        if(b.atom1===ai)bonded[b.atom2]=true;
        else if(b.atom2===ai)bonded[b.atom1]=true;
    }
    for(i=0;i<MD.atoms.length;i++){
        if(i===ai)continue;
        var q=MD.atoms[i];
        var dx=q.x-a.x,dy=q.y-a.y,dz=q.z-a.z;
        var d2=dx*dx+dy*dy+dz*dz;
        if(d2>r2+1e-9)continue;
        out.total++;
        var isBond=!!bonded[i];
        if(isBond)out.bonded++;
        if(neighSkipBonded&&isBond)continue;
        out.rows.push({idx:i,dist:Math.sqrt(d2),bonded:isBond,px:q.x,py:q.y,pz:q.z});
    }
    out.rows.sort(neighByDistance);
    return out;
}

function neighComputeCrystal(){
    var out={rows:[],total:0,bonded:0};
    var abi=neighAnchor.bi,abase=CRY&&CRY.baseAtoms?CRY.baseAtoms[abi]:null;
    if(!abase||!neighBaseVisible(abase))return out;
    var as=neighAnchor.shift,lv=CRY.latticeVectors;
    var fa=neighFractional(abase);
    var widths=crystalPlaneWidths(lv);
    // A displacement inside the radius cannot stray further than r/width along
    // any axis, so this box holds every image that could qualify and nothing
    // else. The repeat cap is the same one the rest of the page uses.
    var lim=[1,1,1],k;
    for(k=0;k<3;k++){
        var q=widths[k]>1e-9?neighRadius/widths[k]:CRYSTAL_MAX_IMAGE_REPEAT;
        lim[k]=Math.min(CRYSTAL_MAX_IMAGE_REPEAT,Math.max(1,Math.ceil(q+1e-9)));
    }
    var bondSet=neighCrystalBondSet();
    var disp=neighDisplayMap();
    var r2=neighRadius*neighRadius;
    var ax0=abase.x+as[0]*lv[0][0]+as[1]*lv[1][0]+as[2]*lv[2][0];
    var ay0=abase.y+as[0]*lv[0][1]+as[1]*lv[1][1]+as[2]*lv[2][1];
    var az0=abase.z+as[0]*lv[0][2]+as[1]*lv[1][2]+as[2]*lv[2][2];
    for(var bi=0;bi<CRY.baseAtoms.length;bi++){
        var bj=CRY.baseAtoms[bi];
        if(!neighBaseVisible(bj))continue;
        var fj=neighFractional(bj);
        var gx=fj[0]-fa[0],gy=fj[1]-fa[1],gz=fj[2]-fa[2];
        var loX=Math.ceil(-gx-lim[0]-1e-9),hiX=Math.floor(-gx+lim[0]+1e-9);
        var loY=Math.ceil(-gy-lim[1]-1e-9),hiY=Math.floor(-gy+lim[1]+1e-9);
        var loZ=Math.ceil(-gz-lim[2]-1e-9),hiZ=Math.floor(-gz+lim[2]+1e-9);
        for(var nx=loX;nx<=hiX;nx++){
            var dx=gx+nx;
            for(var ny=loY;ny<=hiY;ny++){
                var dy=gy+ny;
                for(var nz=loZ;nz<=hiZ;nz++){
                    var dz=gz+nz;
                    var vx=dx*lv[0][0]+dy*lv[1][0]+dz*lv[2][0];
                    var vy=dx*lv[0][1]+dy*lv[1][1]+dz*lv[2][1];
                    var vz=dx*lv[0][2]+dy*lv[1][2]+dz*lv[2][2];
                    var d2=vx*vx+vy*vy+vz*vz;
                    if(d2>r2+1e-9)continue;
                    if(bi===abi&&nx===0&&ny===0&&nz===0)continue;
                    out.total++;
                    var isBond=!!bondSet[abi+'_'+bi+'_'+nx+','+ny+','+nz];
                    if(isBond)out.bonded++;
                    if(neighSkipBonded&&isBond)continue;
                    var di=disp[bi+'_'+(as[0]+nx)+','+(as[1]+ny)+','+(as[2]+nz)];
                    out.rows.push({bi:bi,s:[nx,ny,nz],dist:Math.sqrt(d2),bonded:isBond,
                        idx:(di===undefined?-1:di),px:ax0+vx,py:ay0+vy,pz:az0+vz});
                }
            }
        }
    }
    out.rows.sort(neighByDistance);
    return out;
}

function neighAnchorText(){
    if(!neighAnchor)return 'click an atom';
    if(neighAnchor.bi!=null&&CRY&&CRY.baseAtoms){
        var ba=CRY.baseAtoms[neighAnchor.bi];
        if(ba)return atomLabel(ba,neighAnchor.bi);
    }else{
        var a=MD.atoms[neighAnchor.idx];
        if(a)return atomLabel(a,neighAnchor.idx);
    }
    return 'click an atom';
}

function neighRowLabel(row){
    if(row.bi!=null&&CRY&&CRY.baseAtoms){
        var ba=CRY.baseAtoms[row.bi];
        var nm=atomLabel(ba,row.bi);
        var s=row.s;
        if(s[0]||s[1]||s[2])nm+=' ['+s[0]+','+s[1]+','+s[2]+']';
        return nm;
    }
    return atomLabel(MD.atoms[row.idx],row.idx);
}

function neighRowTitle(row){
    if(row.bi!=null&&CRY&&CRY.baseAtoms){
        var ba=CRY.baseAtoms[row.bi];
        var t=(ba&&ba.siteLabel)?('site '+ba.siteLabel):('site #'+(row.bi+1));
        return t+(row.idx>=0?' - click to select this image':' - this image is outside the drawn supercell');
    }
    return 'click to select this atom';
}

function neighPanelHint(){
    var el=document.getElementById('ng-hint');
    if(!el)return;
    if(CRY&&!CRY_INV)el.textContent='The cell has no inverse, so only the atoms on screen are searched.';
    else if(CRY)el.textContent='Click a row to select that atom. Crystal images outside the drawn cell are listed with their cell offset.';
    else el.textContent='Click a row to select that atom.';
}

function neighClearLines(){
    if(!neighLines)return;
    if(neighLines.parent)moleculeGroup.remove(neighLines);
    if(neighLines.geometry)neighLines.geometry.dispose();
    if(neighLines.material)neighLines.material.dispose();
    neighLines=null;
    needsRender=true;
}

// One dashed segment per listed contact, from the anchor to the image itself.
// An image outside the drawn cell still gets its line, so the picture shows
// where the contact is even when the supercell does not reach it.
function neighDrawLines(shown){
    neighClearLines();
    if(!neighActive||diffMode)return;
    var p=neighAnchorPoint();
    if(!p)return;
    var pts=[],i;
    for(i=0;i<shown&&i<neighRows.length;i++){
        var r=neighRows[i];
        if(r.px==null)continue;
        pts.push(p.x-CX,p.y-CY,p.z-CZ,r.px-CX,r.py-CY,r.pz-CZ);
    }
    if(pts.length===0)return;
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(pts,3));
    var mat=new THREE.LineDashedMaterial({color:0xffb300,dashSize:0.16,gapSize:0.12,transparent:true,opacity:0.9});
    var lines=new THREE.LineSegments(geo,mat);
    lines.computeLineDistances();
    moleculeGroup.add(lines);
    neighLines=lines;
    needsRender=true;
}

function neighRefresh(){
    if(!neighActive)return;
    var anchorEl=document.getElementById('ng-anchor'),sumEl=neighSummaryEl(),listEl=neighListEl2();
    if(diffMode){
        neighRows=[];
        if(anchorEl)anchorEl.textContent='-';
        if(sumEl)sumEl.textContent='Neighbours are not available in Diff.';
        if(listEl)listEl.innerHTML='<div class="ng-empty">Close the comparison to search contacts.</div>';
        neighClearLines();
        return;
    }
    if(neighAnchor&&neighAnchor.bi!=null){
        var ba=CRY&&CRY.baseAtoms?CRY.baseAtoms[neighAnchor.bi]:null;
        if(!ba||!neighBaseVisible(ba))neighAnchor=null;
    }else if(neighAnchor&&!MD.atoms[neighAnchor.idx]){
        neighAnchor=null;
    }
    if(anchorEl)anchorEl.textContent=neighAnchorText();
    neighPanelHint();
    if(!neighAnchor){
        neighRows=[];
        if(sumEl)sumEl.textContent='Click an atom to choose the search centre.';
        if(listEl)listEl.innerHTML='<div class="ng-empty">No anchor atom.</div>';
        neighClearLines();
        return;
    }
    var res=(CRY&&CRY_INV&&neighAnchor.bi!=null)?neighComputeCrystal():neighComputeMolecule();
    neighRows=res.rows;
    var shown=Math.min(neighRows.length,NEIGH_MAX_ROWS);
    var txt;
    if(neighRows.length===0)txt='No contacts within '+neighRadius.toFixed(2)+' A';
    else txt=neighRows.length+' contact'+(neighRows.length===1?'':'s')+' within '+neighRadius.toFixed(2)+' A';
    if(neighSkipBonded&&res.bonded>0)txt+=' ('+res.bonded+' bonded hidden)';
    if(neighRows.length>shown)txt+=' - nearest '+shown+' shown';
    if(sumEl)sumEl.textContent=txt;
    var html='';
    for(var i=0;i<shown;i++){
        var row=neighRows[i];
        var off=row.idx<0;
        html+='<div class="ng-item'+(off?' ng-off':'')+'" data-row="'+i+'" title="'+escapeHtml(neighRowTitle(row))+'">'+
            '<span class="ng-name">'+escapeHtml(neighRowLabel(row))+'</span>'+
            '<span class="ng-d">'+row.dist.toFixed(3)+' A</span>'+
            (off?'<span class="ng-cell">image</span>':'')+
            '</div>';
    }
    if(listEl)listEl.innerHTML=html||'<div class="ng-empty">No contacts.</div>';
    neighDrawLines(shown);
}

function neighSetAnchor(idx){
    var a=MD.atoms[idx];
    if(!a||!neighActive)return;
    if(CRY&&CRY_INV&&a.baseIdx!=null&&a.imgShift){
        neighAnchor={bi:a.baseIdx,shift:[a.imgShift[0],a.imgShift[1],a.imgShift[2]]};
    }else{
        neighAnchor={idx:idx};
    }
    neighRefresh();
}

function neighStartSession(){
    neighActive=true;
    neighAnchor=null;
    neighRows=[];
    var rEl=document.getElementById('ng-radius');
    if(rEl)rEl.value=neighRadius.toFixed(1);
    var sEl=document.getElementById('ng-skip');
    if(sEl)sEl.checked=neighSkipBonded;
    // An atom already selected when the mode opens is the obvious centre.
    if(selectedAtoms.length===1&&MD.atoms[selectedAtoms[0]])neighSetAnchor(selectedAtoms[0]);
    else neighRefresh();
}

function neighEndSession(){
    neighActive=false;
    neighClearLines();
    neighRows=[];
    neighAnchor=null;
}

function neighRowClick(rowIndex){
    var row=neighRows[rowIndex];
    if(!row)return;
    if(row.idx<0){
        modeInfoEl.textContent='Neighbour '+neighRowLabel(row)+' is outside the drawn supercell - grow the bounds to select it';
        return;
    }
    var ai=neighAnchorDisplayIndex();
    selectedAtoms=(ai>=0&&ai!==row.idx)?[ai,row.idx]:[row.idx];
    diffSelectedAtoms=[];originalCoords=null;
    highlightSelected();
    selInfoEl.textContent='Selected: '+formatAtomList(selectedAtoms)+' ('+row.dist.toFixed(3)+' A from the anchor)';
    modeInfoEl.textContent='Neighbour '+neighRowLabel(row)+' at '+row.dist.toFixed(3)+' A';
}

function neighSelectListed(){
    var picks=[],seen={},i;
    for(i=0;i<neighRows.length;i++){
        var row=neighRows[i];
        if(row.idx<0||seen[row.idx])continue;
        seen[row.idx]=true;
        picks.push(row.idx);
    }
    if(picks.length===0){
        modeInfoEl.textContent='Neighbours: none of the listed images is in the drawn supercell';
        return;
    }
    selectedAtoms=picks;diffSelectedAtoms=[];originalCoords=null;
    highlightSelected();
    selInfoEl.textContent='Selected: '+formatAtomList(selectedAtoms)+' ('+picks.length+' neighbour'+(picks.length===1?'':'s')+')';
    modeInfoEl.textContent='Selected '+picks.length+' listed neighbour'+(picks.length===1?'':'s');
}

(function(){
    var panel=document.getElementById('neigh-panel'),list=neighListEl2();
    var skip=document.getElementById('ng-skip'),radius=document.getElementById('ng-radius');
    var pick=document.getElementById('ng-select');
    if(!panel||!list||!skip||!radius||!pick)return;
    radius.addEventListener('input',function(){
        var v=parseFloat(this.value);
        // Let an empty or half-typed value be: the next keystroke may finish it.
        if(!(v>0)||v>NEIGH_MAX_RADIUS)return;
        neighRadius=v;
        neighRefresh();
    });
    radius.addEventListener('change',function(){
        var v=parseFloat(this.value);
        if(!(v>0))v=3.5;
        if(v>NEIGH_MAX_RADIUS)v=NEIGH_MAX_RADIUS;
        neighRadius=v;this.value=v.toFixed(1);
        neighRefresh();
    });
    skip.addEventListener('change',function(){neighSkipBonded=this.checked;neighRefresh()});
    pick.addEventListener('click',neighSelectListed);
    list.addEventListener('click',function(e){
        var node=e.target;
        while(node&&node!==list&&!(node.classList&&node.classList.contains('ng-item')))node=node.parentNode;
        if(!node||node===list)return;
        var idx=parseInt(node.getAttribute('data-row'),10);
        if(idx>=0)neighRowClick(idx);
    });
})();

// ===== Non-covalent interactions =====
// A geometric reading of the contacts the drawing cannot show: hydrogen bonds,
// and the sigma-hole contacts of halogens, chalcogens and pnictogens. The
// criteria are stated here and nowhere else; the panel lists every contact
// they match, each drawn as a dashed segment coloured by its type. A crystal is
// read through the base cell's bond table and integer images, so a contact
// across a cell face is the contact it really is; a molecule is read from the
// atoms on screen, because that is all there is.
var intActive=false;
var intShowTypes={hb:true,xb:true,cb:true,pn:true,pp:true,chpi:true};
var intShowLines=true;
var intRows=[];
var intLines=null;
var INT_MAX_ROWS=400;
// The pi-pi pair search is the one quadratic loop in the panel: rings are
// paired through a grid now, but a pathological cell could still offer a huge
// number of image pairs, so the number of image checks is capped and the
// summary says when the cap was reached.
var INT_PP_WORK=2000000;
// Per-type rows, valid until the scene is rebuilt. A checkbox toggle re-reads
// the types already computed instead of repeating the search; a structural
// change bumps intStructRev in rebuildScene() and the cache falls with it.
var intStructRev=0;
var intCache=null;
var INT_SEARCH_RADIUS=4.5;      // longest sigma-hole contact, with margin
var INT_HB_DISTANCE=2.5;        // d(H...A)
var INT_HB_ANGLE=120;           // D-H...A
var INT_SIGMA_ANGLE=130;        // D-X...A
var INT_VDW_SLACK=0.1;          // d(X...A) < r_vdw(X) + r_vdw(A) - 0.1
var INT_ACCEPTORS={N:1,O:1,S:1,F:1,Cl:1,Br:1,I:1};
// The elements that carry a sigma-hole, and the interaction each one names.
var INT_SIGMA={Cl:'xb',Br:'xb',I:'xb',S:'cb',Se:'cb',Te:'cb',P:'pn',As:'pn',Sb:'pn',Bi:'pn'};
var INT_SIGMA_TYPE={xb:1,cb:1,pn:1};
var INT_TYPE_ORDER=['hb','xb','cb','pn','pp','chpi'];
var INT_TYPE_NAME={hb:'H-bond',xb:'halogen bond',cb:'chalcogen bond',pn:'pnictogen bond',
    pp:'\u03c0-\u03c0 contact',chpi:'CH-\u03c0 contact'};
var INT_TYPE_PLURAL={hb:'H-bonds',xb:'halogen bonds',cb:'chalcogen bonds',pn:'pnictogen bonds',
    pp:'\u03c0-\u03c0 contacts',chpi:'CH-\u03c0 contacts'};
var INT_TYPE_TAG={hb:'HB',xb:'XB',cb:'CB',pn:'PN',pp:'\u03c0\u03c0',chpi:'CH\u03c0'};
var INT_TYPE_COLOR={hb:0x4fc3f7,xb:0xffb74d,cb:0xffd54f,pn:0xce93d8,pp:0x81c784,chpi:0xf06292};

// Aromatic rings and the pi contacts they make. A ring is a planar chordless
// cycle of 5-7 atoms drawn from the element set below, with every edge at
// aromatic order (1.5) - or, for a file that states only single bonds, with
// bond lengths only a delocalised ring gives (intRingAromatic). The contact
// criteria are stated here and nowhere else, as the sigma-hole ones are:
// parallel rings need a 2.8-3.8 A plane gap and at most 2.5 A of slip; a
// T-shaped pair needs the centroid of the edge-on ring 3+ A above the other
// ring's face, within 2.5 A of its centre; CH-pi is a C-H hydrogen over a ring
// face, pointing at the ring centroid. Rings are read from the drawn atoms for
// a molecule and from the base cell's bonds plus integer images for a crystal.
var INT_RING_MIN=5;
var INT_RING_MAX=7;
var INT_RING_PLANAR=0.30;       // max atom distance from the ring plane
var INT_MAX_RINGS=3000;
var INT_RING_WORK=3000000;      // cycle-search budget, so a huge lattice cannot stall a redraw
var INT_RING_ELEM={B:1,C:1,N:1,O:1,P:1,S:1,Se:1};
var INT_RING_MAX_SINGLE={
    'B-B':1.60,'B-C':1.55,'B-N':1.46,'B-O':1.39,'B-P':1.90,'B-S':1.86,'B-Se':2.02,
    'C-C':1.47,'C-N':1.42,'C-O':1.41,'C-P':1.84,'C-S':1.76,'C-Se':1.92,
    'N-N':1.42,'N-O':1.40,'N-P':1.76,'N-S':1.73,'N-Se':1.89,
    'O-O':1.44,'O-P':1.70,'O-S':1.70,'O-Se':1.86,
    'P-P':2.24,'P-S':2.12,'P-Se':2.26,
    'S-S':2.02,'S-Se':2.14,'Se-Se':2.30};
var INT_RING_MAX_DOUBLE={
    'B-C':1.38,'B-N':1.36,'C-C':1.42,'C-N':1.38,'C-O':1.34,'C-P':1.70,'C-S':1.66,'C-Se':1.82,
    'N-N':1.32,'N-O':1.30,'N-P':1.62,'N-S':1.62,'N-Se':1.78,
    'O-O':1.30,'O-P':1.56,'O-S':1.58,'O-Se':1.74,
    'P-P':2.06,'P-S':1.98,'P-Se':2.12,
    'S-S':1.90,'S-Se':2.02,'Se-Se':2.18};
var INT_PP_ANGLE=30;            // angle between the ring planes, deg
var INT_PP_DISTANCE=4.5;        // centroid-centroid
var INT_PP_MIN_GAP=2.8;         // interplanar distance
var INT_PP_MAX_GAP=3.8;
var INT_PP_LATERAL=2.5;         // in-plane slip
var INT_PP_PD_OFFSET=1.2;       // slip at which a parallel pair is 'displaced'
var INT_PP_ATOM_CUT=3.9;        // closest atom-atom contact
var INT_T_ANGLE=60;             // angle between the ring planes, deg
var INT_T_DISTANCE=5.5;         // centroid-centroid
var INT_T_HEIGHT=3.0;           // edge centroid above the face plane
var INT_T_LATERAL=2.5;          // edge centroid inside the face
var INT_T_ATOM_CUT=4.0;         // closest atom-atom contact
var INT_PI_MIN_ATOM=2.6;        // bonded or fused rings are not a stacking contact
var INT_CHPI_HC=3.0;            // d(H...centroid)
var INT_CHPI_ANGLE=120;         // C-H...centroid
var INT_CHPI_MIN_HEIGHT=1.0;    // H above the ring plane
var INT_CHPI_OVERHANG=0.8;      // H projection within the ring radius plus this

// Angle at the common origin of two vectors, in degrees.
function intAngleVectors(ax,ay,az,bx,by,bz){
    var la=Math.sqrt(ax*ax+ay*ay+az*az),lb=Math.sqrt(bx*bx+by*by+bz*bz);
    if(la<1e-9||lb<1e-9)return NaN;
    var c=(ax*bx+ay*by+az*bz)/(la*lb);
    if(c>1)c=1;else if(c<-1)c=-1;
    return Math.acos(c)*180/Math.PI;
}

function intPairKey2(a,b){return a<b?(a+'-'+b):(b+'-'+a)}

function intShiftText(s){
    if(!s[0]&&!s[1]&&!s[2])return '';
    return ' ['+s[0]+','+s[1]+','+s[2]+']';
}

// An image pair (ring, ring + t) and (ring + t, ring) are the same physical
// contact; of t and -t the lexicographically positive one names it.
function intShiftCanonical(s){
    if(s[0]!==0)return s[0]>0;
    if(s[1]!==0)return s[1]>0;
    return s[2]>0;
}

// A cycle's identity independent of the direction it was walked: atom order
// sorted, with the integer image each atom sits in. The minimum-index atom is
// the root and has shift zero, so the key is the same in both directions.
function intCycleKey(idx,shifts){
    var pairs=[],i;
    for(i=0;i<idx.length;i++)pairs.push({i:idx[i],s:shifts[i]});
    pairs.sort(function(a,b){return a.i-b.i});
    var out=[];
    for(i=0;i<pairs.length;i++)out.push(pairs[i].i+':'+pairs[i].s[0]+','+pairs[i].s[1]+','+pairs[i].s[2]);
    return out.join(';');
}

// Every chordless cycle of INT_RING_MIN..INT_RING_MAX atoms, once per root:
// the smallest atom in the cycle is the root, so a cycle is found only from
// it. A path never takes an atom adjacent to an earlier path atom (that would
// be a chord); the one exception is an atom adjacent to the root, which can
// only be the closing atom. A periodic cycle whose accumulated image shifts
// do not cancel is a helix, not a ring, and is dropped. The callback receives
// (atoms, edge orders, per-atom shifts relative to the root); returning false
// stops the search.
function intRingCycles(adj,fn){
    var n=adj.length,mark=new Int32Array(n),token=0,stop=false,work=0;
    var root=0,pathAtoms=[],pathEdges=[],pathShifts=[];
    function edgeTo(j,i){
        var list=adj[j],k;
        if(!list)return null;
        for(k=0;k<list.length;k++)if(list[k].j===i)return list[k];
        return null;
    }
    function dfs(cur){
        if(stop)return;
        var L=pathAtoms.length,list=adj[cur]||[],e;
        for(e=0;e<list.length;e++){
            if(++work>INT_RING_WORK){stop=true;return;}
            var en=list[e],j=en.j;
            if(j<root||j===root)continue;
            if(mark[j]===token)continue;
            var chord=false,q;
            for(q=1;q<L-1;q++){
                if(edgeTo(j,pathAtoms[q])){chord=true;break;}
            }
            if(chord)continue;
            var closeEdge=(cur===root)?null:edgeTo(j,root);
            var ns=[pathShifts[L-1][0]+en.sh[0],
                    pathShifts[L-1][1]+en.sh[1],
                    pathShifts[L-1][2]+en.sh[2]];
            if(closeEdge){
                if(L+1>=INT_RING_MIN&&L+1<=INT_RING_MAX&&
                   ns[0]+closeEdge.sh[0]===0&&ns[1]+closeEdge.sh[1]===0&&ns[2]+closeEdge.sh[2]===0){
                    if(fn(pathAtoms.concat([j]),pathEdges.concat([en.order,closeEdge.order]),pathShifts.concat([ns]))===false){
                        stop=true;return;
                    }
                }
                continue;
            }
            if(L>=INT_RING_MAX)continue;
            mark[j]=token;
            pathAtoms.push(j);pathEdges.push(en.order);pathShifts.push(ns);
            dfs(j);
            pathAtoms.pop();pathEdges.pop();pathShifts.pop();
            mark[j]=0;
            if(stop)return;
        }
    }
    for(root=0;root<n;root++){
        var rl=adj[root];
        if(!rl||rl.length<2)continue;
        token++;
        mark[root]=token;
        pathAtoms=[root];pathEdges=[];pathShifts=[[0,0,0]];
        dfs(root);
        if(stop)break;
    }
}

function intMakeRingLabel(atoms,idx,shifts){
    var s='ring(',i;
    for(i=0;i<idx.length;i++){
        if(i)s+=',';
        s+=atomLabel(atoms[idx[i]],idx[i]);
        if(shifts&&(shifts[i][0]||shifts[i][1]||shifts[i][2]))s+=intShiftText(shifts[i]);
    }
    return s+')';
}

// The list row names a ring by its first atoms and its size; the full atom
// list stays in the row's title, where it can be read without truncation.
function intMakeRingShortLabel(atoms,idx){
    var s='ring(',lim=Math.min(3,idx.length),i;
    for(i=0;i<lim;i++){
        if(i)s+=',';
        s+=atomLabel(atoms[idx[i]],idx[i]);
    }
    if(idx.length>lim)s+='+'+(idx.length-lim);
    return s+')';
}

// A ring is aromatic when every edge has the perceived aromatic order, or -
// for a file that states only single bonds (a PDB CONECT table, a CIF with no
// bond orders) - when a planar 5/6-ring of sp2-capable atoms has either all
// short single edges or a Kekule alternation, each edge inside the length band
// its element pair allows for that order. Saturated rings are longer than the
// single-bond ceilings (cyclopentane C-C 1.54, THF C-O 1.43), so they cannot
// pass by accident.
function intRingAromatic(elem,orders,lens,deg){
    var n=elem.length,i,allArom=true,ord;
    for(i=0;i<n;i++)if(!isAromaticOrder(orders[i])){allArom=false;break;}
    if(allArom)return true;
    if(n!==5&&n!==6)return false;
    for(i=0;i<n;i++){
        var z=deg[i];
        if(typeof z==='number'&&z>3)return false;
    }
    var allSingle=true,dbl=0,pkey,mx;
    for(i=0;i<n;i++){
        ord=orders[i];
        pkey=intPairKey2(elem[i],elem[(i+1)%n]);
        if(ord>=2.5)return false;
        if(ord>=1.75){
            allSingle=false;dbl++;
            mx=INT_RING_MAX_DOUBLE[pkey];
            if(!(mx&&lens[i]<=mx))return false;
        }else if(ord<1.25){
            mx=INT_RING_MAX_SINGLE[pkey];
            if(!(mx&&lens[i]<=mx))return false;
        }else return false;
    }
    if(allSingle)return true;
    if(dbl!==Math.floor(n/2))return false;
    var shared=new Array(n);
    for(i=0;i<n;i++)shared[i]=0;
    for(i=0;i<n;i++)if(orders[i]>=1.75){shared[i]++;shared[(i+1)%n]++;}
    for(i=0;i<n;i++)if(shared[i]>1)return false;
    return true;
}

// Turn a cycle into a ring record: positions (an image shift applied for a
// crystal), centroid, Newell normal, planarity, radius and the atom offsets
// from the centroid. An off-plane or non-aromatic cycle is not a ring.
function intMakeRing(idx,orders,shifts,cfg){
    var n=idx.length,i,atoms=cfg.atoms;
    if(n<INT_RING_MIN||n>INT_RING_MAX)return null;
    var px=[],py=[],pz=[],elem=[],deg=[];
    for(i=0;i<n;i++){
        var ai=idx[i],a=atoms[ai];
        if(!a||isGhostAtom(a)||!INT_RING_ELEM[a.element])return null;
        var sh=shifts[i],x=a.x,y=a.y,z=a.z;
        if(cfg.lv){
            x+=sh[0]*cfg.lv[0][0]+sh[1]*cfg.lv[1][0]+sh[2]*cfg.lv[2][0];
            y+=sh[0]*cfg.lv[0][1]+sh[1]*cfg.lv[1][1]+sh[2]*cfg.lv[2][1];
            z+=sh[0]*cfg.lv[0][2]+sh[1]*cfg.lv[1][2]+sh[2]*cfg.lv[2][2];
        }
        px.push(x);py.push(y);pz.push(z);
        elem.push(a.element);deg.push(cfg.deg[ai]);
    }
    var cx=0,cy=0,cz=0;
    for(i=0;i<n;i++){cx+=px[i];cy+=py[i];cz+=pz[i];}
    cx/=n;cy/=n;cz/=n;
    var nx=0,ny=0,nz=0,rx=[],ry=[],rz=[];
    for(i=0;i<n;i++){
        var k=(i+1)%n;
        var xi=px[i]-cx,yi=py[i]-cy,zi=pz[i]-cz;
        var xj=px[k]-cx,yj=py[k]-cy,zj=pz[k]-cz;
        rx.push(xi);ry.push(yi);rz.push(zi);
        nx+=yi*zj-zi*yj;ny+=zi*xj-xi*zj;nz+=xi*yj-yi*xj;
    }
    var nl=Math.sqrt(nx*nx+ny*ny+nz*nz);
    if(!(nl>1e-9))return null;
    nx/=nl;ny/=nl;nz/=nl;
    var dev=0,radius=0;
    for(i=0;i<n;i++){
        var d=Math.abs(rx[i]*nx+ry[i]*ny+rz[i]*nz);
        if(d>dev)dev=d;
        var rr=Math.sqrt(rx[i]*rx[i]+ry[i]*ry[i]+rz[i]*rz[i]);
        if(rr>radius)radius=rr;
    }
    if(!(dev<=INT_RING_PLANAR))return null;
    var lens=[];
    for(i=0;i<n;i++){
        var k2=(i+1)%n;
        var dx=px[k2]-px[i],dy=py[k2]-py[i],dz=pz[k2]-pz[i];
        lens.push(Math.sqrt(dx*dx+dy*dy+dz*dz));
    }
    if(!intRingAromatic(elem,orders,lens,deg))return null;
    var ring={
        idx:idx.slice(0),orders:orders.slice(0),shifts:shifts.slice(0),n:n,
        cx:cx,cy:cy,cz:cz,nx:nx,ny:ny,nz:nz,radius:radius,
        rel:[],label:intMakeRingLabel(atoms,idx,shifts),
        shortLabel:intMakeRingShortLabel(atoms,idx)
    };
    for(i=0;i<n;i++)ring.rel.push([rx[i],ry[i],rz[i]]);
    if(cfg.lv){
        var fx=0,fy=0,fz=0;
        for(i=0;i<n;i++){
            var f=cfg.fracs[idx[i]];
            fx+=f[0]+shifts[i][0];fy+=f[1]+shifts[i][1];fz+=f[2]+shifts[i][2];
        }
        ring.fc=[fx/n,fy/n,fz/n];
    }
    return ring;
}

function intCollectRings(adj,cfg){
    var rings=[],seen={};
    intRingCycles(adj,function(idx,orders,shifts){
        if(rings.length>=INT_MAX_RINGS)return false;
        var key=intCycleKey(idx,shifts);
        if(seen[key])return true;
        seen[key]=1;
        var r=intMakeRing(idx,orders,shifts,cfg);
        if(r)rings.push(r);
        return rings.length<INT_MAX_RINGS;
    });
    return rings;
}

// The ring graph keeps only edges whose two ends could belong to an aromatic
// ring, and the full degree counts every bond (hydrogens included) so the
// length fallback can tell an sp2 centre from an sp3 one.
function intMoleculeRingGraph(atoms,bonds){
    var adj=new Array(atoms.length),deg=new Array(atoms.length),i;
    for(i=0;i<atoms.length;i++){adj[i]=[];deg[i]=0;}
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        if(!b||b.atom1==null||b.atom2==null)continue;
        var i1=b.atom1,i2=b.atom2;
        if(i1===i2)continue;
        var a1=atoms[i1],a2=atoms[i2];
        if(!a1||!a2||isGhostAtom(a1)||isGhostAtom(a2))continue;
        deg[i1]++;deg[i2]++;
        if(!INT_RING_ELEM[a1.element]||!INT_RING_ELEM[a2.element])continue;
        var order=(typeof b.order==='number')?b.order:1;
        if(order>=2.75)continue;
        adj[i1].push({j:i2,order:order,sh:ZERO_SHIFT});
        adj[i2].push({j:i1,order:order,sh:ZERO_SHIFT});
    }
    return {adj:adj,deg:deg};
}

function intBuildMoleculeRings(){
    var atoms=MD.atoms||[],bonds=MD.bonds||[];
    if(atoms.length<INT_RING_MIN)return [];
    var g=intMoleculeRingGraph(atoms,bonds);
    return intCollectRings(g.adj,{atoms:atoms,deg:g.deg});
}

function intCrystalRingGraph(base,bonds){
    var n=base.length,adj=new Array(n),deg=new Array(n),i;
    for(i=0;i<n;i++){adj[i]=[];deg[i]=0;}
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        if(!b||b.atom1==null||b.atom2==null)continue;
        var i1=b.atom1,i2=b.atom2;
        if(i1===i2)continue;
        var a1=base[i1],a2=base[i2];
        if(!a1||!a2||isGhostAtom(a1)||isGhostAtom(a2))continue;
        if(!neighBaseVisible(a1)||!neighBaseVisible(a2))continue;
        deg[i1]++;deg[i2]++;
        if(!INT_RING_ELEM[a1.element]||!INT_RING_ELEM[a2.element])continue;
        var order=(typeof b.order==='number')?b.order:1;
        if(order>=2.75)continue;
        var s=b.shift||ZERO_SHIFT;
        adj[i1].push({j:i2,order:order,sh:[s[0],s[1],s[2]]});
        adj[i2].push({j:i1,order:order,sh:[-s[0],-s[1],-s[2]]});
    }
    return {adj:adj,deg:deg};
}

function intBuildCrystalRings(fracs){
    var base=CRY.baseAtoms||[],bonds=CRY.baseBonds||[],lv=CRY.latticeVectors;
    if(!lv||base.length<INT_RING_MIN)return [];
    var g=intCrystalRingGraph(base,bonds);
    return intCollectRings(g.adj,{atoms:base,deg:g.deg,lv:lv,fracs:fracs||intBaseFracs()});
}

// A lower bound on the shortest image of a fractional displacement, from the
// spacing of the lattice planes each axis moves a point across. A pair whose
// bound is already past the cutoff owns no image inside it.
function intFracLowerBound(fx,fy,fz,widths){
    var w0=widths[0]||0,w1=widths[1]||0,w2=widths[2]||0;
    var ax=Math.abs(fx-Math.round(fx)),
        ay=Math.abs(fy-Math.round(fy)),
        az=Math.abs(fz-Math.round(fz));
    return Math.max(ax*w0,ay*w1,az*w2);
}

// A wrapped fractional bin grid over candidate sites, one cell wide in each
// plane spacing. Moving a fractional coordinate by step_k along axis k changes
// the distance to the k-th plane family by step_k*width_k = radius, so two
// sites within the radius of each other are always in the same bin or in a
// neighbouring one - the grid can only skip pairs that are out of reach.
// Bins are counted from floor(f/step) and wrapped at the cell, the convention
// the CH-pi grid already uses.
function intFracGridBuild(bins,widths,radius){
    var step=[],cells=[],offs=[],grid={},i,q;
    for(q=0;q<3;q++){
        var w=widths[q]>1e-9?widths[q]:1;
        step[q]=w>radius?radius/w:1;
        cells[q]=Math.max(1,Math.ceil(1/step[q]-1e-9));
        // A wrap of -1 and +1 is the same bin when the grid has fewer than
        // three cells; the offsets are chosen so every bin is visited once.
        offs[q]=cells[q]>=3?[-1,0,1]:(cells[q]===2?[0,1]:[0]);
    }
    for(i=0;i<bins.length;i++){
        var f=bins[i],s='';
        for(q=0;q<3;q++){
            var v=Math.floor(f[q]/step[q])%cells[q];
            if(v<0)v+=cells[q];
            s+=(q?',':'')+v;
        }
        (grid[s]||(grid[s]=[])).push(i);
    }
    return {step:step,cells:cells,offs:offs,grid:grid,count:bins.length};
}
function intFracGridVisit(g,f,fn){
    if(!g||!g.count)return;
    var st=g.step,ce=g.cells,of=g.offs;
    var bx=Math.floor(f[0]/st[0]),by=Math.floor(f[1]/st[1]),bz=Math.floor(f[2]/st[2]);
    for(var a=0;a<of[0].length;a++){
        var cx=((bx+of[0][a])%ce[0]+ce[0])%ce[0];
        for(var b=0;b<of[1].length;b++){
            var cy=((by+of[1][b])%ce[1]+ce[1])%ce[1];
            for(var c=0;c<of[2].length;c++){
                var cz=((bz+of[2][c])%ce[2]+ce[2])%ce[2];
                var list=g.grid[cx+','+cy+','+cz];
                if(!list)continue;
                for(var i=0;i<list.length;i++)fn(list[i]);
            }
        }
    }
}
// The molecule counterpart: a plain Cartesian grid with no wrapping, since a
// molecule has no periodic images. Buckets hold atom indices themselves.
function intCartGridBuild(atoms,list,cell){
    var grid={},i;
    for(i=0;i<list.length;i++){
        var a=atoms[list[i]];
        var key=Math.floor(a.x/cell)+','+Math.floor(a.y/cell)+','+Math.floor(a.z/cell);
        (grid[key]||(grid[key]=[])).push(list[i]);
    }
    return {cell:cell,grid:grid,count:list.length};
}
function intCartGridVisit(g,x,y,z,fn){
    if(!g||!g.count)return;
    var cell=g.cell,gx=Math.floor(x/cell),gy=Math.floor(y/cell),gz=Math.floor(z/cell);
    for(var ox=-1;ox<=1;ox++)for(var oy=-1;oy<=1;oy++)for(var oz=-1;oz<=1;oz++){
        var list=g.grid[(gx+ox)+','+(gy+oy)+','+(gz+oz)];
        if(!list)continue;
        for(var i=0;i<list.length;i++)fn(list[i]);
    }
}

// Every image of one fractional point inside a sphere around another, with
// its integer shift and Cartesian displacement. The optional budget carries
// the caller's remaining image tests: one is spent per image tried, and the
// walk stops when the budget is gone, so a caller pairing thousands of sites
// cannot stall the page.
function intForEachPointImage(center,site,lim,radius,fn,budget){
    var lv=CRY.latticeVectors,r2=radius*radius;
    var gx=site[0]-center[0],gy=site[1]-center[1],gz=site[2]-center[2];
    var loX=Math.ceil(-gx-lim[0]-1e-9),hiX=Math.floor(-gx+lim[0]+1e-9);
    var loY=Math.ceil(-gy-lim[1]-1e-9),hiY=Math.floor(-gy+lim[1]+1e-9);
    var loZ=Math.ceil(-gz-lim[2]-1e-9),hiZ=Math.floor(-gz+lim[2]+1e-9);
    for(var nx=loX;nx<=hiX;nx++){
        var dx=gx+nx;
        for(var ny=loY;ny<=hiY;ny++){
            var dy=gy+ny;
            for(var nz=loZ;nz<=hiZ;nz++){
                if(budget&&--budget.left<0){budget.hit=true;return}
                var dz=gz+nz;
                var vx=dx*lv[0][0]+dy*lv[1][0]+dz*lv[2][0];
                var vy=dx*lv[0][1]+dy*lv[1][1]+dz*lv[2][1];
                var vz=dx*lv[0][2]+dy*lv[1][2]+dz*lv[2][2];
                var d2=vx*vx+vy*vy+vz*vz;
                if(d2>r2+1e-9)continue;
                fn([nx,ny,nz],vx,vy,vz,Math.sqrt(d2));
            }
        }
    }
}

// Closest atom-atom contact between two rings when the second ring's centroid
// is at dx,dy,dz from the first's. Ring atom offsets are relative to their own
// centroids, so no per-image array has to be built.
function intRingMinContact(a,b,dx,dy,dz){
    var min=Infinity,ra=a.rel,rb=b.rel,i,j;
    for(i=0;i<ra.length;i++){
        var ax=ra[i][0]-dx,ay=ra[i][1]-dy,az=ra[i][2]-dz;
        for(j=0;j<rb.length;j++){
            var vx=ax-rb[j][0],vy=ay-rb[j][1],vz=az-rb[j][2];
            var d2=vx*vx+vy*vy+vz*vz;
            if(d2<min)min=d2;
        }
    }
    return Math.sqrt(min);
}

// The pi-pi criteria: one parallel family (plane angle at most 30 degrees,
// centroid distance at most 4.5 A, plane gap 2.8-3.8 A, slip at most 2.5 A)
// and one T family (planes at least 60 degrees apart, an edge centroid 3+ A
// above the other ring's face and within 2.5 A of its centre, centroid
// distance at most 5.5 A). Both demand a closest atom contact between 2.6 and
// 3.9/4.0 A, which keeps bonded and fused rings out. A parallel pair is called
// displaced once its slip reaches 1.2 A.
function intPiRow(a,b,dx,dy,dz){
    var R=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(!(R<=INT_T_DISTANCE))return null;
    var c=Math.abs(a.nx*b.nx+a.ny*b.ny+a.nz*b.nz);
    if(c>1)c=1;
    var theta=Math.acos(c)*180/Math.PI;
    var contact=intRingMinContact(a,b,dx,dy,dz);
    if(!(contact>=INT_PI_MIN_ATOM))return null;
    if(theta<=INT_PP_ANGLE&&R<=INT_PP_DISTANCE&&contact<=INT_PP_ATOM_CUT){
        var gap=Math.abs(dx*a.nx+dy*a.ny+dz*a.nz);
        if(gap>=INT_PP_MIN_GAP&&gap<=INT_PP_MAX_GAP){
            var lat=Math.sqrt(Math.max(0,R*R-gap*gap));
            if(lat<=INT_PP_LATERAL){
                return {type:'pp',subtype:(lat>=INT_PP_PD_OFFSET?'PD':'P'),d:R,angle:theta};
            }
        }
    }
    if(theta>=INT_T_ANGLE&&contact<=INT_T_ATOM_CUT){
        var hA=Math.abs(dx*a.nx+dy*a.ny+dz*a.nz);
        var hB=Math.abs(dx*b.nx+dy*b.ny+dz*b.nz);
        var latA=Math.sqrt(Math.max(0,R*R-hA*hA));
        var latB=Math.sqrt(Math.max(0,R*R-hB*hB));
        if((hA>=INT_T_HEIGHT&&latA<=INT_T_LATERAL)||(hB>=INT_T_HEIGHT&&latB<=INT_T_LATERAL)){
            return {type:'pp',subtype:'T',d:R,angle:theta};
        }
    }
    return null;
}

function intRingImageLabel(r,t){
    var s=r.shortLabel||r.label;
    if(t[0]||t[1]||t[2])s+=' ['+t[0]+','+t[1]+','+t[2]+']';
    return s;
}

function intCollectMoleculePiPi(rings,rows,trunc){
    if(!rings||!rings.length)return;
    var cell=INT_T_DISTANCE,grid={},i;
    for(i=0;i<rings.length;i++){
        var r=rings[i];
        var key=Math.floor(r.cx/cell)+','+Math.floor(r.cy/cell)+','+Math.floor(r.cz/cell);
        (grid[key]||(grid[key]=[])).push(i);
    }
    var budget={left:INT_PP_WORK,hit:false};
    for(i=0;i<rings.length&&!budget.hit;i++){
        var a=rings[i];
        var gx=Math.floor(a.cx/cell),gy=Math.floor(a.cy/cell),gz=Math.floor(a.cz/cell);
        for(var ox=-1;ox<=1&&!budget.hit;ox++)for(var oy=-1;oy<=1&&!budget.hit;oy++)for(var oz=-1;oz<=1&&!budget.hit;oz++){
            var list=grid[(gx+ox)+','+(gy+oy)+','+(gz+oz)];
            if(!list)continue;
            for(var q=0;q<list.length;q++){
                var j=list[q];
                if(j<=i)continue;
                if(--budget.left<0){budget.hit=true;break}
                var b=rings[j];
                var dx=b.cx-a.cx,dy=b.cy-a.cy,dz=b.cz-a.cz;
                if(dx*dx+dy*dy+dz*dz>INT_T_DISTANCE*INT_T_DISTANCE)continue;
                var row=intPiRow(a,b,dx,dy,dz);
                if(!row)continue;
                row.mol=true;
                row.ringALabel=a.shortLabel;row.ringBLabel=b.shortLabel;
                row.ringAFull=a.label;row.ringBFull=b.label;
                row.cx=a.cx;row.cy=a.cy;row.cz=a.cz;
                row.ax=b.cx;row.ay=b.cy;row.az=b.cz;
                row.selAtoms=a.idx.concat(b.idx);
                row.key='ppm_'+i+'_'+j;
                rows.push(row);
            }
        }
    }
    if(budget.hit&&trunc)trunc.pp=true;
}

// One crystal pi-pi row; the second ring is the image t of ring b. The same
// construction serves the ring's own images (b is a) and a distinct pair.
function intCrystalPiRow(a,b,ti,tj,t,vx,vy,vz,rows){
    var row=intPiRow(a,b,vx,vy,vz);
    if(!row)return;
    row.ringALabel=a.shortLabel;
    row.ringBLabel=intRingImageLabel(b,t);
    row.ringAFull=a.label;row.ringBFull=b.label;
    row.cx=a.cx;row.cy=a.cy;row.cz=a.cz;
    row.ax=a.cx+vx;row.ay=a.cy+vy;row.az=a.cz+vz;
    row.sitesA=[];row.sitesB=[];
    var k;
    for(k=0;k<a.idx.length;k++)row.sitesA.push({bi:a.idx[k],sh:a.shifts[k]});
    for(k=0;k<b.idx.length;k++)row.sitesB.push({bi:b.idx[k],
        sh:[b.shifts[k][0]+t[0],b.shifts[k][1]+t[1],b.shifts[k][2]+t[2]]});
    row.key='ppc_'+ti+'_'+tj+'_'+t.join(',');
    rows.push(row);
}

// Ring centroids are paired through a fractional grid at the centroid cutoff,
// so a cell with thousands of rings no longer runs the all-pairs scan. Every
// pair that owns any image within reach still lands in this bin or a
// neighbour; each image tested is charged to INT_PP_WORK, and the summary says
// when that cap stopped the walk.
function intCollectCrystalPiPi(rings,rows,trunc){
    if(!rings||!rings.length)return;
    var lv=CRY.latticeVectors,widths=crystalPlaneWidths(lv);
    var lim=intImageLimit(widths,INT_T_DISTANCE),i;
    var bins=new Array(rings.length);
    for(i=0;i<rings.length;i++)bins[i]=rings[i].fc;
    var grid=intFracGridBuild(bins,widths,INT_T_DISTANCE);
    var budget={left:INT_PP_WORK,hit:false};
    for(i=0;i<rings.length&&!budget.hit;i++){
        var a=rings[i];
        var self=i;
        // The ring's own images: reachable only in a small cell, and the
        // canonical shift keeps one name for the pair t/-t.
        intForEachPointImage(a.fc,a.fc,lim,INT_T_DISTANCE,function(t,vx,vy,vz){
            if(t[0]===0&&t[1]===0&&t[2]===0)return;
            if(!intShiftCanonical(t))return;
            intCrystalPiRow(a,a,self,self,t,vx,vy,vz,rows);
        },budget);
        if(budget.hit)break;
        intFracGridVisit(grid,a.fc,function(j){
            if(budget.hit||j<=self)return;
            var b=rings[j];
            var fx=b.fc[0]-a.fc[0],fy=b.fc[1]-a.fc[1],fz=b.fc[2]-a.fc[2];
            if(intFracLowerBound(fx,fy,fz,widths)>INT_T_DISTANCE+1e-9)return;
            intForEachPointImage(a.fc,b.fc,lim,INT_T_DISTANCE,function(t,vx,vy,vz){
                if(budget.hit)return;
                intCrystalPiRow(a,b,self,j,t,vx,vy,vz,rows);
            },budget);
        });
    }
    if(budget.hit&&trunc)trunc.pp=true;
}

function intAdjHits(adj,idx,set){
    var list=adj[idx];
    if(!list)return false;
    for(var i=0;i<list.length;i++)if(set[list[i]])return true;
    return false;
}

function intMoleculeHList(atoms,bonds){
    var adj=new Array(atoms.length),i;
    for(i=0;i<atoms.length;i++)adj[i]=[];
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        if(!b||b.atom1==null||b.atom2==null)continue;
        var i1=b.atom1,i2=b.atom2;
        if(i1===i2)continue;
        var a1=atoms[i1],a2=atoms[i2];
        if(!a1||!a2||isGhostAtom(a1)||isGhostAtom(a2))continue;
        adj[i1].push(i2);adj[i2].push(i1);
    }
    var list=[];
    for(i=0;i<atoms.length;i++){
        var a=atoms[i];
        if(!a||a.element!=='H'||isGhostAtom(a))continue;
        var nb=adj[i];
        if(nb.length!==1)continue;
        var c=nb[0],ca=atoms[c];
        if(!ca||ca.element!=='C'||isGhostAtom(ca))continue;
        list.push({h:i,c:c});
    }
    return {list:list,adj:adj};
}

// CH-pi: a carbon-bound H within 3.0 A of a ring centroid, pointing at it
// (C-H...centroid > 120 degrees), at least 1.0 A above the ring plane and
// projecting inside the ring face. An H or C bonded to the ring itself is not
// a CH-pi contact and is skipped.
function intCollectMoleculeChPi(rings,rows){
    if(!rings||!rings.length)return;
    var atoms=MD.atoms||[],bonds=MD.bonds||[];
    var hs=intMoleculeHList(atoms,bonds);
    if(!hs.list.length)return;
    var cell=INT_CHPI_HC,grid={},i;
    for(i=0;i<hs.list.length;i++){
        var H=atoms[hs.list[i].h];
        var key=Math.floor(H.x/cell)+','+Math.floor(H.y/cell)+','+Math.floor(H.z/cell);
        (grid[key]||(grid[key]=[])).push(i);
    }
    for(i=0;i<rings.length;i++){
        var r=rings[i],set={},q;
        for(q=0;q<r.idx.length;q++)set[r.idx[q]]=1;
        var gx=Math.floor(r.cx/cell),gy=Math.floor(r.cy/cell),gz=Math.floor(r.cz/cell);
        for(var ox=-1;ox<=1;ox++)for(var oy=-1;oy<=1;oy++)for(var oz=-1;oz<=1;oz++){
            var list=grid[(gx+ox)+','+(gy+oy)+','+(gz+oz)];
            if(!list)continue;
            for(var li=0;li<list.length;li++){
                var rec=hs.list[list[li]],hIdx=rec.h,cIdx=rec.c;
                if(set[hIdx]||set[cIdx])continue;
                if(intAdjHits(hs.adj,hIdx,set)||intAdjHits(hs.adj,cIdx,set))continue;
                var hh=atoms[hIdx],cc=atoms[cIdx];
                var dx=hh.x-r.cx,dy=hh.y-r.cy,dz=hh.z-r.cz;
                var d2=dx*dx+dy*dy+dz*dz;
                if(!(d2<INT_CHPI_HC*INT_CHPI_HC))continue;
                var d=Math.sqrt(d2);
                var ang=intAngleVectors(cc.x-hh.x,cc.y-hh.y,cc.z-hh.z,-dx,-dy,-dz);
                if(!(ang>INT_CHPI_ANGLE))continue;
                var hpl=Math.abs(dx*r.nx+dy*r.ny+dz*r.nz);
                if(hpl<INT_CHPI_MIN_HEIGHT)continue;
                var lat=Math.sqrt(Math.max(0,d2-hpl*hpl));
                if(lat>r.radius+INT_CHPI_OVERHANG)continue;
                rows.push({type:'chpi',mol:true,d:d,angle:ang,
                    hLabel:atomLabel(hh,hIdx),cLabel:atomLabel(cc,cIdx),
                    ringLabel:r.shortLabel,ringFull:r.label,
                    cx:hh.x,cy:hh.y,cz:hh.z,ax:r.cx,ay:r.cy,az:r.cz,
                    selAtoms:r.idx.concat([cIdx,hIdx]),
                    key:'chpim_'+i+'_'+hIdx});
            }
        }
    }
}

// The crystal CH-pi search fixes the ring in its home cell and walks the H
// images near it; the C-H bond's own integer shift carries the carbon along.
// The H and C are skipped when a base-cell bond already joins either to a ring
// atom (an N-H or a benzylic C-H pointing at its own ring is not CH-pi).
function intCollectCrystalChPi(rings,rows,fracs){
    if(!rings||!rings.length)return;
    var base=CRY.baseAtoms,lv=CRY.latticeVectors,bonds=CRY.baseBonds||[];
    if(!base||!lv)return;
    if(!fracs)fracs=intBaseFracs();
    var widths=crystalPlaneWidths(lv);
    var lim=intImageLimit(widths,INT_CHPI_HC);
    var adj=new Array(base.length),shiftMap={},i;
    for(i=0;i<base.length;i++)adj[i]=[];
    for(i=0;i<bonds.length;i++){
        var b=bonds[i];
        if(!b||b.atom1==null||b.atom2==null)continue;
        var i1=b.atom1,i2=b.atom2;
        if(i1===i2)continue;
        var a1=base[i1],a2=base[i2];
        if(!a1||!a2||isGhostAtom(a1)||isGhostAtom(a2))continue;
        adj[i1].push(i2);adj[i2].push(i1);
        var s=b.shift||ZERO_SHIFT;
        if(shiftMap[i1+'_'+i2]===undefined){
            shiftMap[i1+'_'+i2]=[s[0],s[1],s[2]];
            shiftMap[i2+'_'+i1]=[-s[0],-s[1],-s[2]];
        }
    }
    var hlist=[];
    for(i=0;i<base.length;i++){
        var a=base[i];
        if(!a||a.element!=='H'||isGhostAtom(a)||!neighBaseVisible(a))continue;
        var nb=adj[i];
        if(!nb||nb.length!==1)continue;
        var c=nb[0],ca=base[c];
        if(!ca||ca.element!=='C'||isGhostAtom(ca)||!neighBaseVisible(ca))continue;
        var csh=shiftMap[i+'_'+c];
        if(!csh)continue;
        hlist.push({h:i,c:c,csh:csh});
    }
    if(!hlist.length)return;
    // A fractional-space grid whose cells are the search radius wide along
    // each axis: an H within reach of a ring has to be in the same or an
    // adjacent cell, wrapped through the cell faces. It keeps the scan local
    // in a supercell that carries thousands of ring hydrogens.
    var step=[],cells=[],hgrid={},q;
    for(q=0;q<3;q++){
        step[q]=widths[q]>INT_CHPI_HC?INT_CHPI_HC/widths[q]:1;
        cells[q]=Math.max(1,Math.ceil(1/step[q]-1e-9));
    }
    function hCell(f,k){
        var v=Math.floor(f/step[k])%cells[k];
        return v<0?v+cells[k]:v;
    }
    for(q=0;q<hlist.length;q++){
        var hf=fracs[hlist[q].h];
        var hk=hCell(hf[0],0)+','+hCell(hf[1],1)+','+hCell(hf[2],2);
        (hgrid[hk]||(hgrid[hk]=[])).push(q);
    }
    for(i=0;i<rings.length;i++){
        var r=rings[i],set={};
        for(q=0;q<r.idx.length;q++)set[r.idx[q]]=1;
        var ci=[hCell(r.fc[0],0),hCell(r.fc[1],1),hCell(r.fc[2],2)];
        for(var ox=-1;ox<=1;ox++)for(var oy=-1;oy<=1;oy++)for(var oz=-1;oz<=1;oz++){
            var gx=((ci[0]+ox)%cells[0]+cells[0])%cells[0],
                gy=((ci[1]+oy)%cells[1]+cells[1])%cells[1],
                gz=((ci[2]+oz)%cells[2]+cells[2])%cells[2];
            var cand=hgrid[gx+','+gy+','+gz];
            if(!cand)continue;
            for(var li=0;li<cand.length;li++){
                var rec=hlist[cand[li]];
                if(set[rec.h]||set[rec.c])continue;
                if(intAdjHits(adj,rec.h,set)||intAdjHits(adj,rec.c,set))continue;
                var fh=fracs[rec.h];
                if(intFracLowerBound(fh[0]-r.fc[0],fh[1]-r.fc[1],fh[2]-r.fc[2],widths)>INT_CHPI_HC+1e-9)continue;
                intForEachPointImage(r.fc,fh,lim,INT_CHPI_HC,function(t,vx,vy,vz,d){
                    if(!(d<INT_CHPI_HC))return;
                    var hx=r.cx+vx,hy=r.cy+vy,hz=r.cz+vz;
                    var hb=base[rec.h],cb=base[rec.c];
                    var ux=cb.x-hb.x+rec.csh[0]*lv[0][0]+rec.csh[1]*lv[1][0]+rec.csh[2]*lv[2][0];
                    var uy=cb.y-hb.y+rec.csh[0]*lv[0][1]+rec.csh[1]*lv[1][1]+rec.csh[2]*lv[2][1];
                    var uz=cb.z-hb.z+rec.csh[0]*lv[0][2]+rec.csh[1]*lv[1][2]+rec.csh[2]*lv[2][2];
                    var ang=intAngleVectors(ux,uy,uz,r.cx-hx,r.cy-hy,r.cz-hz);
                    if(!(ang>INT_CHPI_ANGLE))return;
                    var hpl=Math.abs((hx-r.cx)*r.nx+(hy-r.cy)*r.ny+(hz-r.cz)*r.nz);
                    if(hpl<INT_CHPI_MIN_HEIGHT)return;
                    var lat=Math.sqrt(Math.max(0,d*d-hpl*hpl));
                    if(lat>r.radius+INT_CHPI_OVERHANG)return;
                    var csh=[t[0]+rec.csh[0],t[1]+rec.csh[1],t[2]+rec.csh[2]];
                    var row={type:'chpi',d:d,angle:ang,
                        hLabel:atomLabel(hb,rec.h)+intShiftText(t),
                        cLabel:atomLabel(cb,rec.c)+intShiftText(csh),
                        ringLabel:r.shortLabel,ringFull:r.label,
                        cx:hx,cy:hy,cz:hz,ax:r.cx,ay:r.cy,az:r.cz,
                        sitesA:[],sitesB:[
                            {bi:rec.c,sh:csh},
                            {bi:rec.h,sh:[t[0],t[1],t[2]]}],
                        key:'chpic_'+i+'_'+rec.h+'_'+t.join(',')};
                    for(var k=0;k<r.idx.length;k++)row.sitesA.push({bi:r.idx[k],sh:r.shifts[k]});
                    rows.push(row);
                });
            }
        }
    }
}

function intByDistance(a,b){
    if(Math.abs(a.d-b.d)>1e-9)return a.d-b.d;
    if(a.type!==b.type)return a.type<b.type?-1:1;
    if(!!a.mol!==!!b.mol)return a.mol?1:-1;
    if(a.key!=null||b.key!=null){
        var ak=a.key||'',bk=b.key||'';
        if(ak!==bk)return ak<bk?-1:1;
        return 0;
    }
    if(a.mol)return (a.dIdx-b.dIdx)||((a.aIdx||0)-(b.aIdx||0));
    return (a.dBi-b.dBi)||(a.aBi-b.aBi)||
        ((a.aShift[0]-b.aShift[0])||(a.aShift[1]-b.aShift[1])||(a.aShift[2]-b.aShift[2]));
}

// One entry per physical X...A contact: the substituents of X are only the
// reference the sigma-hole angle is measured from, so the most linear one wins.
// Hydrogen bonds and the ring contacts are already one row per contact.
function intDedupeSigma(rows){
    var out=[],map={},i;
    for(i=0;i<rows.length;i++){
        var r=rows[i];
        if(!INT_SIGMA_TYPE[r.type]){out.push(r);continue;}
        var key=r.mol?('m'+r.type+'_'+r.xIdx+'_'+r.aIdx)
                     :('c'+r.type+'_'+r.xBi+'_'+r.aBi+'_'+r.aShift.join(','));
        var at=map[key];
        if(at===undefined){map[key]=out.length;out.push(r);}
        else if(r.angle>out[at].angle){out[at]=r;}
    }
    return out;
}

// The molecule search's fixed data: the donors, the sigma-hole centres and the
// acceptor grid, built once per structure and shared by every type that reads
// it. Each donor then looks only at acceptors in its own or a neighbouring
// cell of a Cartesian grid whose cell is the search radius, instead of every
// atom in the structure.
function intMoleculePrep(){
    var atoms=MD.atoms||[],bonds=MD.bonds||[],i;
    var hbonds=[],sigmas=[],acceptors=[];
    for(i=0;i<bonds.length;i++){
        var bb=bonds[i];
        if(!bb)continue;
        var a=atoms[bb.atom1],b=atoms[bb.atom2];
        if(!a||!b||isGhostAtom(a)||isGhostAtom(b))continue;
        if(b.element==='H'&&a.element!=='H')hbonds.push({d:bb.atom1,h:bb.atom2});
        else if(a.element==='H'&&b.element!=='H')hbonds.push({d:bb.atom2,h:bb.atom1});
        if(INT_SIGMA[a.element])sigmas.push({type:INT_SIGMA[a.element],x:bb.atom1,d:bb.atom2});
        if(INT_SIGMA[b.element])sigmas.push({type:INT_SIGMA[b.element],x:bb.atom2,d:bb.atom1});
    }
    for(i=0;i<atoms.length;i++){
        var at=atoms[i];
        if(isGhostAtom(at)||!INT_ACCEPTORS[at.element])continue;
        acceptors.push(i);
    }
    return {atoms:atoms,hbonds:hbonds,sigmas:sigmas,rings:null,
        accGridHb:intCartGridBuild(atoms,acceptors,INT_HB_DISTANCE),
        accGridSigma:intCartGridBuild(atoms,acceptors,INT_SEARCH_RADIUS)};
}
function intMoleculeHbRows(prep,rows){
    var atoms=prep.atoms;
    for(var i=0;i<prep.hbonds.length;i++){
        var p=prep.hbonds[i],D=atoms[p.d],H=atoms[p.h];
        if(!D||!H)continue;
        intCartGridVisit(prep.accGridHb,H.x,H.y,H.z,function(ai){
            if(ai===p.h||ai===p.d)return;
            var A=atoms[ai];
            var dx=A.x-H.x,dy=A.y-H.y,dz=A.z-H.z;
            var d=Math.sqrt(dx*dx+dy*dy+dz*dz);
            if(!(d<INT_HB_DISTANCE))return;
            var ang=intAngleVectors(D.x-H.x,D.y-H.y,D.z-H.z,dx,dy,dz);
            if(!(ang>INT_HB_ANGLE))return;
            rows.push({type:'hb',mol:true,d:d,angle:ang,dIdx:p.d,hIdx:p.h,aIdx:ai,
                dLabel:atomLabel(D,p.d),hLabel:atomLabel(H,p.h),aLabel:atomLabel(A,ai),
                cx:H.x,cy:H.y,cz:H.z,ax:A.x,ay:A.y,az:A.z});
        });
    }
}
function intMoleculeSigmaRows(prep,type,rows){
    var atoms=prep.atoms;
    for(var i=0;i<prep.sigmas.length;i++){
        var q=prep.sigmas[i];
        if(q.type!==type)continue;
        var X=atoms[q.x],sub=atoms[q.d];
        if(!X||!sub)continue;
        var svx=sub.x-X.x,svy=sub.y-X.y,svz=sub.z-X.z;
        intCartGridVisit(prep.accGridSigma,X.x,X.y,X.z,function(aj){
            if(aj===q.x||aj===q.d)return;
            var A2=atoms[aj];
            var vx=A2.x-X.x,vy=A2.y-X.y,vz=A2.z-X.z;
            var dd=Math.sqrt(vx*vx+vy*vy+vz*vz);
            var cut=(VDW_R[X.element]||DEFAULT_VDW_R)+(VDW_R[A2.element]||DEFAULT_VDW_R)-INT_VDW_SLACK;
            if(!(dd<cut))return;
            var ag=intAngleVectors(svx,svy,svz,vx,vy,vz);
            if(!(ag>INT_SIGMA_ANGLE))return;
            rows.push({type:q.type,mol:true,d:dd,angle:ag,xIdx:q.x,dIdx:q.d,aIdx:aj,
                dLabel:atomLabel(sub,q.d),xLabel:atomLabel(X,q.x),aLabel:atomLabel(A2,aj),
                cx:X.x,cy:X.y,cz:X.z,ax:A2.x,ay:A2.y,az:A2.z});
        });
    }
}

function intBaseFracs(){
    var out=new Array(CRY.baseAtoms.length),i;
    for(i=0;i<out.length;i++)out[i]=neighFractional(CRY.baseAtoms[i]);
    return out;
}

function intImageLimit(widths,radius){
    radius=radius==null?INT_SEARCH_RADIUS:radius;
    var lim=[1,1,1],k;
    for(k=0;k<3;k++){
        var q=widths[k]>1e-9?radius/widths[k]:CRYSTAL_MAX_IMAGE_REPEAT;
        lim[k]=Math.min(CRYSTAL_MAX_IMAGE_REPEAT,Math.max(1,Math.ceil(q+1e-9)));
    }
    return lim;
}

// The crystal search's fixed data. The acceptors are binned in fractional
// space at the search radius, so a donor reads the few sites its own or a
// neighbouring cell holds instead of every acceptor in the cell; the donors
// and sigma-hole centres are read from the cell's bond table once.
function intCrystalPrep(){
    var base=CRY.baseAtoms,lv=CRY.latticeVectors;
    var fracs=intBaseFracs();
    var widths=crystalPlaneWidths(lv);
    var lim=intImageLimit(widths);
    var limHb=intImageLimit(widths,INT_HB_DISTANCE);
    var accList=[],accBins=[],bi;
    for(bi=0;bi<base.length;bi++){
        var ba=base[bi];
        if(!ba||isGhostAtom(ba)||!neighBaseVisible(ba))continue;
        if(!INT_ACCEPTORS[ba.element])continue;
        accList.push(bi);
        accBins.push(fracs[bi]);
    }
    var hbonds=[],sigmas=[],bonds=CRY.baseBonds||[],i;
    for(i=0;i<bonds.length;i++){
        var bb=bonds[i];
        if(!bb||bb.atom1==null||bb.atom2==null)continue;
        var a1=base[bb.atom1],a2=base[bb.atom2];
        if(!a1||!a2||isGhostAtom(a1)||isGhostAtom(a2))continue;
        if(!neighBaseVisible(a1)||!neighBaseVisible(a2))continue;
        var s=bb.shift||ZERO_SHIFT;
        if(a2.element==='H'&&a1.element!=='H')hbonds.push({d:bb.atom1,h:bb.atom2,hShift:[s[0],s[1],s[2]]});
        else if(a1.element==='H'&&a2.element!=='H')hbonds.push({d:bb.atom2,h:bb.atom1,hShift:[-s[0],-s[1],-s[2]]});
        if(INT_SIGMA[a2.element])sigmas.push({type:INT_SIGMA[a2.element],x:bb.atom2,d:bb.atom1,dShift:[s[0],s[1],s[2]]});
        if(INT_SIGMA[a1.element])sigmas.push({type:INT_SIGMA[a1.element],x:bb.atom1,d:bb.atom2,dShift:[-s[0],-s[1],-s[2]]});
    }
    return {base:base,lv:lv,fracs:fracs,widths:widths,lim:lim,limHb:limHb,
        hbonds:hbonds,sigmas:sigmas,accList:accList,rings:null,
        accGridHb:intFracGridBuild(accBins,widths,INT_HB_DISTANCE),
        accGridSigma:intFracGridBuild(accBins,widths,INT_SEARCH_RADIUS)};
}
function intCrystalHbRows(prep,rows){
    var base=prep.base,lv=prep.lv,fracs=prep.fracs;
    for(var i=0;i<prep.hbonds.length;i++){
        var p=prep.hbonds[i],db=base[p.d],hb=base[p.h];
        if(!db||!hb)continue;
        var hs=p.hShift;
        var hx=hb.x+hs[0]*lv[0][0]+hs[1]*lv[1][0]+hs[2]*lv[2][0];
        var hy=hb.y+hs[0]*lv[0][1]+hs[1]*lv[1][1]+hs[2]*lv[2][1];
        var hz=hb.z+hs[0]*lv[0][2]+hs[1]*lv[1][2]+hs[2]*lv[2][2];
        var dvx=db.x-hx,dvy=db.y-hy,dvz=db.z-hz;
        var fc=[fracs[p.h][0]+hs[0],fracs[p.h][1]+hs[1],fracs[p.h][2]+hs[2]];
        intFracGridVisit(prep.accGridHb,fc,function(m){
            var ai=prep.accList[m];
            intForEachPointImage(fc,fracs[ai],prep.limHb,INT_HB_DISTANCE,function(shift,vx,vy,vz,d){
                if(!(d<INT_HB_DISTANCE))return;
                if(ai===p.h&&shift[0]===0&&shift[1]===0&&shift[2]===0)return;
                if(ai===p.d&&shift[0]===0&&shift[1]===0&&shift[2]===0)return;
                var ang=intAngleVectors(dvx,dvy,dvz,vx,vy,vz);
                if(!(ang>INT_HB_ANGLE))return;
                rows.push({type:'hb',d:d,angle:ang,
                    dBi:p.d,hBi:p.h,hShift:[p.hShift[0],p.hShift[1],p.hShift[2]],
                    aBi:ai,aShift:[shift[0],shift[1],shift[2]],
                    dLabel:atomLabel(db,p.d),hLabel:atomLabel(hb,p.h),
                    aLabel:atomLabel(base[ai],ai)+((shift[0]||shift[1]||shift[2])?' ['+shift[0]+','+shift[1]+','+shift[2]+']':''),
                    cx:hx,cy:hy,cz:hz,ax:hx+vx,ay:hy+vy,az:hz+vz});
            });
        });
    }
}
function intCrystalSigmaRows(prep,type,rows){
    var base=prep.base,lv=prep.lv,fracs=prep.fracs;
    for(var i=0;i<prep.sigmas.length;i++){
        var q=prep.sigmas[i];
        if(q.type!==type)continue;
        var xb=base[q.x],sub=base[q.d];
        if(!xb||!sub)continue;
        var ps=q.dShift;
        var svx=sub.x+ps[0]*lv[0][0]+ps[1]*lv[1][0]+ps[2]*lv[2][0]-xb.x;
        var svy=sub.y+ps[0]*lv[0][1]+ps[1]*lv[1][1]+ps[2]*lv[2][1]-xb.y;
        var svz=sub.z+ps[0]*lv[0][2]+ps[1]*lv[1][2]+ps[2]*lv[2][2]-xb.z;
        intFracGridVisit(prep.accGridSigma,fracs[q.x],function(m){
            var ai=prep.accList[m];
            intForEachPointImage(fracs[q.x],fracs[ai],prep.lim,INT_SEARCH_RADIUS,function(shift,vx,vy,vz,d){
                if(ai===q.x&&shift[0]===0&&shift[1]===0&&shift[2]===0)return;
                var el=base[ai].element;
                var cut=(VDW_R[xb.element]||DEFAULT_VDW_R)+(VDW_R[el]||DEFAULT_VDW_R)-INT_VDW_SLACK;
                if(!(d<cut))return;
                var ang=intAngleVectors(svx,svy,svz,vx,vy,vz);
                if(!(ang>INT_SIGMA_ANGLE))return;
                rows.push({type:q.type,d:d,angle:ang,
                    xBi:q.x,dBi:q.d,aBi:ai,aShift:[shift[0],shift[1],shift[2]],
                    dLabel:atomLabel(sub,q.d),xLabel:atomLabel(xb,q.x),
                    aLabel:atomLabel(base[ai],ai)+((shift[0]||shift[1]||shift[2])?' ['+shift[0]+','+shift[1]+','+shift[2]+']':''),
                    cx:xb.x,cy:xb.y,cz:xb.z,ax:xb.x+vx,ay:xb.y+vy,az:xb.z+vz});
            });
        });
    }
}

// The rows of every type asked for, computed once per structure. A checkbox
// that turns a type off and on again re-reads its cached rows; the search runs
// again only for a type that has not been computed yet, or after a scene
// rebuild (which bumps intStructRev).
function intCollect(){
    var crystal=!!(CRY&&CRY_INV&&CRY.baseAtoms&&CRY.baseAtoms.length&&CRY.latticeVectors);
    if(!intCache||intCache.rev!==intStructRev||intCache.crystal!==crystal){
        intCache={rev:intStructRev,crystal:crystal,rows:{},prep:null,trunc:{}};
    }
    var c=intCache;
    if(crystal){
        if(!c.prep)c.prep=intCrystalPrep();
        var cp=c.prep;
        if(intShowTypes.hb&&!c.rows.hb){c.rows.hb=[];intCrystalHbRows(cp,c.rows.hb);}
        if(intShowTypes.xb&&!c.rows.xb){c.rows.xb=[];intCrystalSigmaRows(cp,'xb',c.rows.xb);}
        if(intShowTypes.cb&&!c.rows.cb){c.rows.cb=[];intCrystalSigmaRows(cp,'cb',c.rows.cb);}
        if(intShowTypes.pn&&!c.rows.pn){c.rows.pn=[];intCrystalSigmaRows(cp,'pn',c.rows.pn);}
        if((intShowTypes.pp||intShowTypes.chpi)&&!cp.rings)cp.rings=intBuildCrystalRings(cp.fracs);
        if(intShowTypes.pp&&!c.rows.pp){c.rows.pp=[];intCollectCrystalPiPi(cp.rings,c.rows.pp,c.trunc);}
        if(intShowTypes.chpi&&!c.rows.chpi){c.rows.chpi=[];intCollectCrystalChPi(cp.rings,c.rows.chpi,cp.fracs);}
    }else{
        if(!c.prep)c.prep=intMoleculePrep();
        var mp=c.prep;
        if(intShowTypes.hb&&!c.rows.hb){c.rows.hb=[];intMoleculeHbRows(mp,c.rows.hb);}
        if(intShowTypes.xb&&!c.rows.xb){c.rows.xb=[];intMoleculeSigmaRows(mp,'xb',c.rows.xb);}
        if(intShowTypes.cb&&!c.rows.cb){c.rows.cb=[];intMoleculeSigmaRows(mp,'cb',c.rows.cb);}
        if(intShowTypes.pn&&!c.rows.pn){c.rows.pn=[];intMoleculeSigmaRows(mp,'pn',c.rows.pn);}
        if((intShowTypes.pp||intShowTypes.chpi)&&!mp.rings)mp.rings=intBuildMoleculeRings();
        if(intShowTypes.pp&&!c.rows.pp){c.rows.pp=[];intCollectMoleculePiPi(mp.rings,c.rows.pp,c.trunc);}
        if(intShowTypes.chpi&&!c.rows.chpi){c.rows.chpi=[];intCollectMoleculeChPi(mp.rings,c.rows.chpi);}
    }
    var rows=[],i,j,list;
    for(i=0;i<INT_TYPE_ORDER.length;i++){
        var t=INT_TYPE_ORDER[i];
        if(!intShowTypes[t])continue;
        list=c.rows[t];
        if(!list)continue;
        for(j=0;j<list.length;j++)rows.push(list[j]);
    }
    rows=intDedupeSigma(rows);
    rows.sort(intByDistance);
    return rows;
}

function intRowTypeName(row){
    if(row.type==='pp'){
        if(row.subtype==='T')return '\u03c0-\u03c0 stacking (T-shaped)';
        if(row.subtype==='PD')return '\u03c0-\u03c0 stacking (parallel-displaced)';
        return '\u03c0-\u03c0 stacking (parallel)';
    }
    return INT_TYPE_NAME[row.type]||row.type;
}

function intRowName(row){
    if(row.type==='hb')return row.dLabel+'-'+row.hLabel+'...'+row.aLabel;
    if(row.type==='pp')return row.ringALabel+'...'+row.ringBLabel;
    if(row.type==='chpi')return row.hLabel+'...'+row.ringLabel;
    return row.dLabel+'-'+row.xLabel+'...'+row.aLabel;
}

function intRowTitle(row){
    var extra='';
    if(row.type==='pp'&&row.ringAFull)extra=' (rings '+row.ringAFull+' and '+row.ringBFull+')';
    else if(row.type==='chpi'&&row.ringFull)extra=' (ring '+row.ringFull+')';
    return intRowTypeName(row)+': '+intRowName(row)+', '+row.d.toFixed(3)+' A, '+
        row.angle.toFixed(1)+' deg'+extra+' - click to select its sites';
}

function intDisplayIndexFor(bi,shift){
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a.baseIdx===bi&&a.imgShift&&
           a.imgShift[0]===shift[0]&&a.imgShift[1]===shift[1]&&a.imgShift[2]===shift[2])return i;
    }
    return -1;
}

// Both ends of a crystal interaction on screen, in the same translated cell.
// The interaction is reported for the donor's home cell, but any drawn image of
// that donor is as good: the acceptor is looked for in the cell the donor image
// actually sits in.
function intDisplayPair(row){
    var di=row.xBi!=null?row.xBi:row.dBi;
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a.baseIdx!==di||!a.imgShift)continue;
        var want=[a.imgShift[0]+row.aShift[0],
                  a.imgShift[1]+row.aShift[1],
                  a.imgShift[2]+row.aShift[2]];
        var j=intDisplayIndexFor(row.aBi,want);
        if(j>=0)return [i,j];
    }
    return null;
}

// A crystal ring contact names several atoms at each end, each with its own
// image shift. Anchor the first site-A atom on any image drawn on screen, then
// take every other site in the cell that anchor's image actually sits in.
function intDisplaySites(row){
    var A=row.sitesA,B=row.sitesB;
    if(!A||!B||!A.length||!B.length)return null;
    var map={},i;
    for(i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a.baseIdx==null||!a.imgShift)continue;
        map[a.baseIdx+'_'+a.imgShift[0]+'_'+a.imgShift[1]+'_'+a.imgShift[2]]=i;
    }
    function at(bi,sh){return map[bi+'_'+sh[0]+'_'+sh[1]+'_'+sh[2]];}
    var anchor=A[0];
    for(i=0;i<MD.atoms.length;i++){
        var d=MD.atoms[i];
        if(d.baseIdx!==anchor.bi||!d.imgShift)continue;
        var ox=d.imgShift[0]-anchor.sh[0],
            oy=d.imgShift[1]-anchor.sh[1],
            oz=d.imgShift[2]-anchor.sh[2];
        var picks=[],ok=true,q,j;
        for(q=0;q<A.length;q++){
            j=at(A[q].bi,[A[q].sh[0]+ox,A[q].sh[1]+oy,A[q].sh[2]+oz]);
            if(j===undefined){ok=false;break;}
            picks.push(j);
        }
        if(!ok)continue;
        for(q=0;q<B.length;q++){
            j=at(B[q].bi,[B[q].sh[0]+ox,B[q].sh[1]+oy,B[q].sh[2]+oz]);
            if(j===undefined){ok=false;break;}
            picks.push(j);
        }
        if(ok)return picks;
    }
    return null;
}

function intRowClick(index){
    var row=intRows[index];
    if(!row)return;
    var picks;
    if(row.selAtoms)picks=row.selAtoms;
    else if(row.mol)picks=[row.xIdx!=null?row.xIdx:row.dIdx,row.aIdx];
    else picks=intDisplaySites(row)||intDisplayPair(row);
    if(!picks||!picks.length){
        modeInfoEl.textContent='Interaction '+intRowName(row)+' is outside the drawn supercell - grow the bounds to select it';
        return;
    }
    selectedAtoms=picks;diffSelectedAtoms=[];originalCoords=null;
    highlightSelected();
    var tname=intRowTypeName(row);
    if(row.selAtoms||row.sitesA){
        selInfoEl.textContent='Selected: '+picks.length+' atoms ('+(INT_TYPE_NAME[row.type]||tname)+', '+row.d.toFixed(3)+' A)';
    }else{
        selInfoEl.textContent='Selected: '+formatAtomList(picks)+' ('+tname+', '+row.d.toFixed(3)+' A)';
    }
    modeInfoEl.textContent=tname+' '+intRowName(row)+': '+row.d.toFixed(3)+
        ' A, '+row.angle.toFixed(1)+' deg';
}

function intClearLines(){
    if(!intLines)return;
    if(intLines.parent)moleculeGroup.remove(intLines);
    if(intLines.geometry)intLines.geometry.dispose();
    if(intLines.material)intLines.material.dispose();
    intLines=null;
    needsRender=true;
}

function intDrawLines(shown){
    intClearLines();
    if(!intActive||diffMode||!intShowLines)return;
    var pts=[],cols=[],count=0,i;
    for(i=0;i<shown.length&&i<INT_MAX_ROWS;i++){
        var r=shown[i];
        if(r.cx==null||r.ax==null)continue;
        var col=new THREE.Color(INT_TYPE_COLOR[r.type]);
        pts.push(r.cx-CX,r.cy-CY,r.cz-CZ,r.ax-CX,r.ay-CY,r.az-CZ);
        cols.push(col.r,col.g,col.b,col.r,col.g,col.b);
        count++;
    }
    if(count===0)return;
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(pts,3));
    geo.setAttribute('color',new THREE.Float32BufferAttribute(cols,3));
    var mat=new THREE.LineDashedMaterial({vertexColors:true,dashSize:0.12,gapSize:0.09,transparent:true,opacity:0.95});
    var lines=new THREE.LineSegments(geo,mat);
    lines.computeLineDistances();
    moleculeGroup.add(lines);
    intLines=lines;
    needsRender=true;
}

function intTypeCountText(t,n){
    if(n===1)return INT_TYPE_NAME[t]||t;
    return INT_TYPE_PLURAL[t]||((INT_TYPE_NAME[t]||t)+'s');
}

function intRefresh(){
    if(!intActive)return;
    var sumEl=document.getElementById('int-summary'),listEl=document.getElementById('int-list');
    if(diffMode){
        intRows=[];
        if(sumEl)sumEl.textContent='Interactions are not available in Diff.';
        if(listEl)listEl.innerHTML='<div class="ng-empty">Close the comparison to detect interactions.</div>';
        intClearLines();
        return;
    }
    var all=intCollect();
    var shown=all.length>INT_MAX_ROWS?all.slice(0,INT_MAX_ROWS):all;
    intRows=shown;
    var count={},i;
    for(i=0;i<INT_TYPE_ORDER.length;i++)count[INT_TYPE_ORDER[i]]=0;
    for(i=0;i<all.length;i++)if(count[all[i].type]!=null)count[all[i].type]++;
    var txt;
    if(all.length===0){
        txt='No non-covalent interactions match the geometry criteria';
    }else{
        var parts=[];
        INT_TYPE_ORDER.forEach(function(t){
            if(count[t]>0)parts.push(count[t]+' '+intTypeCountText(t,count[t]));
        });
        txt=all.length+' interaction'+(all.length===1?'':'s')+': '+parts.join(', ');
        if(all.length>shown.length)txt+=' - nearest '+shown.length+' listed';
    }
    if(intShowTypes.pp&&intCache&&intCache.trunc&&intCache.trunc.pp){
        txt+=' - pi-contact search hit its work budget, some pairs were not checked';
    }
    if(sumEl)sumEl.textContent=txt;
    var html='';
    for(i=0;i<shown.length;i++){
        var row=shown[i];
        var tag=INT_TYPE_TAG[row.type]||row.type;
        if(row.subtype)tag+=' '+row.subtype;
        html+='<div class="ng-item" data-row="'+i+'" title="'+escapeHtml(intRowTitle(row))+'">'+
            '<span class="ng-type ng-type-'+row.type+'">'+escapeHtml(tag)+'</span>'+
            '<span class="ng-name">'+escapeHtml(intRowName(row))+'</span>'+
            '<span class="ng-a">'+row.angle.toFixed(0)+'&deg;</span>'+
            '<span class="ng-d">'+row.d.toFixed(3)+' A</span>'+
            '</div>';
    }
    if(listEl)listEl.innerHTML=html||'<div class="ng-empty">No interactions match the criteria.</div>';
    intDrawLines(shown);
}

function intStartSession(){
    intActive=true;
    intRows=[];
    // A fresh session reads the structure in front of it, not a cache left by
    // an earlier one.
    intCache=null;
    var show=document.getElementById('int-show');
    if(show)show.checked=intShowLines;
    INT_TYPE_ORDER.forEach(function(t){
        var el=document.getElementById('int-'+t);
        if(el)el.checked=!!intShowTypes[t];
    });
    intRefresh();
}

function intEndSession(){
    intActive=false;
    intClearLines();
    intRows=[];
    intCache=null;
}

(function(){
    var list=document.getElementById('int-list');
    if(!list)return;
    INT_TYPE_ORDER.forEach(function(t){
        var el=document.getElementById('int-'+t);
        if(!el)return;
        el.addEventListener('change',function(){intShowTypes[t]=this.checked;intRefresh()});
    });
    var show=document.getElementById('int-show');
    if(show)show.addEventListener('change',function(){intShowLines=this.checked;intRefresh()});
    list.addEventListener('click',function(e){
        var node=e.target;
        while(node&&node!==list&&!(node.classList&&node.classList.contains('ng-item')))node=node.parentNode;
        if(!node||node===list)return;
        var idx=parseInt(node.getAttribute('data-row'),10);
        if(idx>=0)intRowClick(idx);
    });
})();

// ===== Lattice planes (hkl) =====
// The cell panel's spacing calculator: two integer triples, the d-spacing of
// the first and the angle between the two plane normals. "Show Plane" draws
// the (hkl) cross-section through the centre of the drawn box - the one plane
// of the family every set of bounds is guaranteed to cut.
function hklReciprocals(lv){
    if(!lv||lv.length<3)return null;
    var a=lv[0],b=lv[1],c=lv[2];
    var bcx=b[1]*c[2]-b[2]*c[1],bcy=b[2]*c[0]-b[0]*c[2],bcz=b[0]*c[1]-b[1]*c[0];
    var det=a[0]*bcx+a[1]*bcy+a[2]*bcz;
    if(Math.abs(det)<1e-9)return null;
    return [
        {x:bcx/det,y:bcy/det,z:bcz/det},
        {x:(c[1]*a[2]-c[2]*a[1])/det,y:(c[2]*a[0]-c[0]*a[2])/det,z:(c[0]*a[1]-c[1]*a[0])/det},
        {x:(a[1]*b[2]-a[2]*b[1])/det,y:(a[2]*b[0]-a[0]*b[2])/det,z:(a[0]*b[1]-a[1]*b[0])/det}
    ];
}
// The reciprocal vector g = h a* + k b* + l c*: its length is 1/d, its
// direction the plane normal. (0 0 0) has neither.
function hklNormal(rec,h,k,l){
    if(!rec)return null;
    var x=h*rec[0].x+k*rec[1].x+l*rec[2].x;
    var y=h*rec[0].y+k*rec[1].y+l*rec[2].y;
    var z=h*rec[0].z+k*rec[1].z+l*rec[2].z;
    var len=Math.sqrt(x*x+y*y+z*z);
    if(!(len>1e-9))return null;
    return {x:x/len,y:y/len,z:z/len,len:len};
}
function hklAngleDeg(n1,n2){
    if(!n1||!n2)return null;
    // A plane is unoriented, so the acute angle between the normals is the
    // angle between the planes.
    var d=Math.abs(n1.x*n2.x+n1.y*n2.y+n1.z*n2.z);
    if(d>1)return 0;
    return Math.acos(d)*180/Math.PI;
}
function hklIndexText(p){return '('+p.h+' '+p.k+' '+p.l+')'}
function hklUpdateReadout(){
    var d1=document.getElementById('hk-d1'),d2=document.getElementById('hk-d2');
    var rec=CRY?hklReciprocals(CRY.latticeVectors):null;
    if(d1){
        if(!CRY)d1.textContent='';
        else if(!rec)d1.textContent='d = - (the cell has no volume)';
        else{
            var n=hklNormal(rec,hklPlane.h,hklPlane.k,hklPlane.l);
            d1.textContent=n?('d'+hklIndexText(hklPlane)+' = '+(1/n.len).toFixed(4)+' A'):
                ('d'+hklIndexText(hklPlane)+' = -');
        }
    }
    if(d2){
        if(!CRY)d2.textContent='';
        else if(!rec)d2.textContent='angle = -';
        else{
            var a1=hklNormal(rec,hklPlane.h,hklPlane.k,hklPlane.l);
            var a2=hklNormal(rec,hklPlane2.h,hklPlane2.k,hklPlane2.l);
            var ang=hklAngleDeg(a1,a2);
            d2.textContent=ang==null?('angle to '+hklIndexText(hklPlane2)+' = -'):
                ('angle to '+hklIndexText(hklPlane2)+' = '+ang.toFixed(2)+' deg');
        }
    }
}
// The plane through the fractional centre of the drawn box, cut down to the
// box: the 12 edges are tested, the crossings ordered around their centroid.
// A box is convex, so its edge crossings already are the polygon.
function hklPlanePolygon(lv,h,k,l){
    var rec=hklReciprocals(lv),nrm=hklNormal(rec,h,k,l);
    if(!nrm)return null;
    var cf=[(boundary.aMin+boundary.aMax)/2,(boundary.bMin+boundary.bMax)/2,(boundary.cMin+boundary.cMax)/2];
    var cx=cf[0]*lv[0][0]+cf[1]*lv[1][0]+cf[2]*lv[2][0];
    var cy=cf[0]*lv[0][1]+cf[1]*lv[1][1]+cf[2]*lv[2][1];
    var cz=cf[0]*lv[0][2]+cf[1]*lv[1][2]+cf[2]*lv[2][2];
    var av=[boundary.aMin,boundary.aMax],bv=[boundary.bMin,boundary.bMax],cv=[boundary.cMin,boundary.cMax];
    var vs=[],i,j,m;
    for(i=0;i<2;i++)for(j=0;j<2;j++)for(m=0;m<2;m++){
        var x=av[i]*lv[0][0]+bv[j]*lv[1][0]+cv[m]*lv[2][0];
        var y=av[i]*lv[0][1]+bv[j]*lv[1][1]+cv[m]*lv[2][1];
        var z=av[i]*lv[0][2]+bv[j]*lv[1][2]+cv[m]*lv[2][2];
        vs.push([x,y,z,nrm.x*(x-cx)+nrm.y*(y-cy)+nrm.z*(z-cz)]);
    }
    var pts=[];
    for(i=0;i<8;i++)if(Math.abs(vs[i][3])<1e-9)pts.push([vs[i][0],vs[i][1],vs[i][2]]);
    var edges=[[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]];
    for(i=0;i<edges.length;i++){
        var p0=vs[edges[i][0]],p1=vs[edges[i][1]];
        if(p0[3]*p1[3]<-1e-18){
            var t=p0[3]/(p0[3]-p1[3]);
            pts.push([p0[0]+t*(p1[0]-p0[0]),p0[1]+t*(p1[1]-p0[1]),p0[2]+t*(p1[2]-p0[2])]);
        }
    }
    var uniq=[];
    for(i=0;i<pts.length;i++){
        var dup=false;
        for(j=0;j<uniq.length;j++){
            var dx=pts[i][0]-uniq[j][0],dy=pts[i][1]-uniq[j][1],dz=pts[i][2]-uniq[j][2];
            if(dx*dx+dy*dy+dz*dz<1e-12){dup=true;break}
        }
        if(!dup)uniq.push(pts[i]);
    }
    if(uniq.length<3)return null;
    var gx=0,gy=0,gz=0;
    for(i=0;i<uniq.length;i++){gx+=uniq[i][0];gy+=uniq[i][1];gz+=uniq[i][2]}
    gx/=uniq.length;gy/=uniq.length;gz/=uniq.length;
    var ux=uniq[0][0]-gx,uy=uniq[0][1]-gy,uz=uniq[0][2]-gz;
    var ul=Math.sqrt(ux*ux+uy*uy+uz*uz);
    if(ul<1e-9)return null;
    ux/=ul;uy/=ul;uz/=ul;
    var vx=nrm.y*uz-nrm.z*uy,vy=nrm.z*ux-nrm.x*uz,vz=nrm.x*uy-nrm.y*ux;
    // Sort by polar angle in the plane's own frame (u, v = n x u), so the fan
    // has no folds whichever way the box was cut.
    var ordered=uniq.slice();
    function ang2(p){return Math.atan2((p[0]-gx)*vx+(p[1]-gy)*vy+(p[2]-gz)*vz,(p[0]-gx)*ux+(p[1]-gy)*uy+(p[2]-gz)*uz)}
    ordered.sort(function(p,q){return ang2(p)-ang2(q)});
    return ordered;
}
function hklClearPlane(){
    if(hklPlaneMesh){moleculeGroup.remove(hklPlaneMesh);disposeMesh(hklPlaneMesh);hklPlaneMesh=null;needsRender=true}
    if(hklPlaneBorder){moleculeGroup.remove(hklPlaneBorder);disposeMesh(hklPlaneBorder);hklPlaneBorder=null;needsRender=true}
}
function hklRefresh(){
    hklUpdateReadout();
    hklClearPlane();
    if(!hklShown||!CRY||diffMode)return;
    var poly=hklPlanePolygon(CRY.latticeVectors,hklPlane.h,hklPlane.k,hklPlane.l);
    if(!poly)return;
    var tri=[],i;
    for(i=1;i<poly.length-1;i++){
        tri.push(poly[0][0]-CX,poly[0][1]-CY,poly[0][2]-CZ,
                 poly[i][0]-CX,poly[i][1]-CY,poly[i][2]-CZ,
                 poly[i+1][0]-CX,poly[i+1][1]-CY,poly[i+1][2]-CZ);
    }
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(tri,3));
    // transparent and depthWrite off: the plane is a reading aid, not matter,
    // so it never hides an atom behind it.
    hklPlaneMesh=new THREE.Mesh(geo,new THREE.MeshBasicMaterial({color:0x66ccff,transparent:true,opacity:0.22,side:THREE.DoubleSide,depthWrite:false}));
    moleculeGroup.add(hklPlaneMesh);
    var loop=[];
    for(i=0;i<poly.length;i++)loop.push(poly[i][0]-CX,poly[i][1]-CY,poly[i][2]-CZ);
    var bgeo=new THREE.BufferGeometry();
    bgeo.setAttribute('position',new THREE.Float32BufferAttribute(loop,3));
    hklPlaneBorder=new THREE.LineLoop(bgeo,new THREE.LineBasicMaterial({color:0x66ccff,transparent:true,opacity:0.9,depthWrite:false}));
    moleculeGroup.add(hklPlaneBorder);
    needsRender=true;
}
(function(){
    var showBtn=document.getElementById('hk-show');
    if(!showBtn)return;
    function readInt(id,fallback){
        var e=document.getElementById(id);
        if(!e)return fallback;
        var v=parseInt(e.value,10);
        return isNaN(v)?fallback:v;
    }
    function sync(){
        hklPlane={h:readInt('hk-h',0),k:readInt('hk-k',0),l:readInt('hk-l',1)};
        hklPlane2={h:readInt('hk2-h',1),k:readInt('hk2-k',0),l:readInt('hk2-l',0)};
        hklRefresh();
    }
    ['hk-h','hk-k','hk-l','hk2-h','hk2-k','hk2-l'].forEach(function(id){
        var e=document.getElementById(id);
        if(e)e.addEventListener('change',sync);
    });
    showBtn.addEventListener('click',function(){
        hklShown=!hklShown;
        this.textContent=hklShown?'Hide Plane':'Show Plane';
        hklRefresh();
        layoutPanels();
    });
    sync();
})();

// ===== Powder X-ray diffraction =====
// The unit cell drawn on screen is a crystal, so it owns a powder pattern the
// way it owns a spacing and a plane. Every triple (h k l) in the reciprocal
// sphere is enumerated - one loop, no symmetry shortcuts - and each acceptor
// is turned into a structure factor over the cell's atoms, so a systematic
// absence (bcc's (100), fcc's (200)) is zero because the atoms cancel, not
// because a rule was applied. Reflections that share a 2theta naturally add up
// in the same bin, which is exactly what a powder does to them, so no
// multiplicity table is needed either. The peaks are then broadened with the
// Caglioti variance of a laboratory diffractometer and normalised to 100.
//
// The form factors are the International Tables for Crystallography, Vol. C
// (1992), Table 6.1.1.4 coefficients for neutral atoms, laid out as
// [a1,b1,a2,b2,a3,b3,a4,b4,c] with f(s) = sum a_i exp(-b_i s^2) + c and
// s = sin(theta)/lambda in 1/A. Elements outside the table are left out of
// the pattern and named in the panel's note; a deuterium reads as hydrogen.
var XRD_FORM_FACTORS={
H:[0.493002,10.5109,0.322912,26.1257,0.140191,3.14236,0.04081,57.7997,0.003038],
He:[0.8734,9.1037,0.6309,3.3568,0.3112,22.9276,0.178,0.9821,0.0064],
Li:[1.1282,3.9546,0.7508,1.0524,0.6175,85.3905,0.4653,168.261,0.0377],
Be:[1.5919,43.6427,1.1278,1.8623,0.5391,103.483,0.7029,0.542,0.0385],
B:[2.0545,23.2185,1.3326,1.021,1.0979,60.3498,0.7068,0.1403,-0.1932],
C:[2.31,20.8439,1.02,10.2075,1.5886,0.5687,0.865,51.6512,0.2156],
N:[12.2126,0.0057,3.1322,9.8933,2.0125,28.9975,1.1663,0.5826,-11.529],
O:[3.0485,13.2771,2.2868,5.7011,1.5463,0.3239,0.867,32.9089,0.2508],
F:[3.5392,10.2825,2.6412,4.2944,1.517,0.2615,1.0243,26.1476,0.2776],
Ne:[3.9553,8.4042,3.1125,3.4262,1.4546,0.2306,1.1251,21.7184,0.3515],
Na:[4.7626,3.285,3.1736,8.8422,1.2674,0.3136,1.1128,129.424,0.676],
Mg:[5.4204,2.8275,2.1735,79.2611,1.2269,0.3808,2.3073,7.1937,0.8584],
Al:[6.4202,3.0387,1.9002,0.7426,1.5936,31.5472,1.9646,85.0886,1.1151],
Si:[6.2915,2.4386,3.0353,32.3337,1.9891,0.6785,1.541,81.6937,1.1407],
P:[6.4345,1.9067,4.1791,27.157,1.78,0.526,1.4908,68.1645,1.1149],
S:[6.9053,1.4679,5.2034,22.2151,1.4379,0.2536,1.5863,56.172,0.8669],
Cl:[11.4604,0.0104,7.1964,1.1662,6.2556,18.5194,1.6455,47.7784,-9.5574],
Ar:[7.4845,0.9072,6.7723,14.8407,0.6539,43.8983,1.6442,33.3929,1.4445],
K:[8.2186,12.7949,7.4398,0.7748,1.0519,213.187,0.8659,41.6841,1.4228],
Ca:[8.6266,10.4421,7.3873,0.6599,1.5899,85.7484,1.0211,178.437,1.3751],
Sc:[9.189,9.0213,7.3679,0.5729,1.6409,136.108,1.468,51.3531,1.3329],
Ti:[9.7595,7.8508,7.3558,0.5,1.6991,35.6338,1.9021,116.105,1.2807],
V:[10.2971,6.8657,7.3511,0.4385,2.0703,26.8938,2.0571,102.478,1.2199],
Cr:[10.6406,6.1038,7.3537,0.392,3.324,20.2626,1.4922,98.7399,1.1832],
Mn:[11.2819,5.3409,7.3573,0.3432,3.0193,17.8674,2.2441,83.7543,1.0896],
Fe:[11.7695,4.7611,7.3573,0.3072,3.5222,15.3535,2.3045,76.8805,1.0369],
Co:[12.2841,4.2791,7.3409,0.2784,4.0034,13.5359,2.3488,71.1692,1.0118],
Ni:[12.8376,3.8785,7.292,0.2565,4.4438,12.1763,2.38,66.3421,1.0341],
Cu:[13.338,3.5828,7.1676,0.247,5.6158,11.3966,1.6735,64.8126,1.191],
Zn:[14.0743,3.2655,7.0318,0.2333,5.1652,10.3163,2.41,58.7097,1.3041],
Ga:[15.2354,3.0669,6.7006,0.2412,4.3591,10.7805,2.9623,61.4135,1.7189],
Ge:[16.0816,2.8509,6.3747,0.2516,3.7068,11.4468,3.683,54.7625,2.1313],
As:[16.6723,2.6345,6.0701,0.2647,3.4313,12.9479,4.2779,47.7972,2.531],
Se:[17.0006,2.4098,5.8196,0.2726,3.9731,15.2372,4.3543,43.8163,2.8409],
Br:[17.1789,2.1723,5.2358,16.5796,5.6377,0.2609,3.9851,41.4328,2.9557],
Kr:[17.3555,1.9384,6.7286,16.5623,5.5493,0.2261,3.5375,39.3972,2.825],
Rb:[17.1784,1.7888,9.6435,17.3151,5.1399,0.2748,1.5292,164.934,3.4873],
Sr:[17.5663,1.5564,9.8184,14.0988,5.422,0.1664,2.6694,132.376,2.5064],
Y:[17.776,1.4029,10.2946,12.8006,5.72629,0.125599,3.26588,104.354,1.91213],
Zr:[17.8765,1.27618,10.948,11.916,5.41732,0.117622,3.65721,87.6627,2.06929],
Nb:[17.6142,1.18865,12.0144,11.766,4.04183,0.204785,3.53346,69.7957,3.75591],
Mo:[3.7025,0.2772,17.2356,1.0958,12.8876,11.004,3.7429,61.6584,4.3875],
Tc:[19.1301,0.864132,11.0948,8.14487,4.64901,21.5707,2.71263,86.8472,5.40428],
Ru:[19.2674,0.80852,12.9182,8.43467,4.86337,24.7997,1.56756,94.2928,5.37874],
Rh:[19.2957,0.751536,14.3501,8.21758,4.73425,25.8749,1.28918,98.6062,5.328],
Pd:[19.3319,0.698655,15.5017,7.98929,5.29537,25.2052,0.605844,76.8986,5.26593],
Ag:[19.2808,0.6446,16.6885,7.4726,4.8045,24.6605,1.0463,99.8156,5.179],
Cd:[19.2214,0.5946,17.6444,6.9089,4.461,24.7008,1.6029,87.4825,5.0694],
In:[19.1624,0.5476,18.5596,6.3776,4.2948,25.8499,2.0396,92.8029,4.9391],
Sn:[19.1889,5.8303,19.1005,0.5031,4.4585,26.8909,2.4663,83.9571,4.7821],
Sb:[19.6418,5.3034,19.0455,0.4607,5.0371,27.9074,2.6827,75.2825,4.5909],
Te:[19.9644,4.81742,19.0138,0.420885,6.14487,28.5284,2.5239,70.8403,4.352],
I:[20.1472,4.347,18.9949,0.3814,7.5138,27.766,2.2735,66.8776,4.0712],
Xe:[20.2933,3.9282,19.0298,0.344,8.9767,26.4659,1.99,64.2658,3.7118],
Cs:[20.3892,3.569,19.1062,0.3107,10.662,24.3879,1.4953,213.904,3.3352],
Ba:[20.3361,3.216,19.297,0.2756,10.888,20.2073,2.6959,167.202,2.7731],
La:[20.578,2.94817,19.599,0.244475,11.3727,18.7726,3.28719,133.124,2.14678],
Ce:[21.1671,2.81219,19.7695,0.226836,11.8513,17.6083,3.33049,127.113,1.86264],
Pr:[22.044,2.77393,19.6697,0.222087,12.3856,16.7669,2.82428,143.644,2.0583],
Nd:[22.6845,2.66248,19.6847,0.210628,12.774,15.885,2.85137,137.903,1.98486],
Pm:[23.3405,2.5627,19.6095,0.202088,13.1235,15.1009,2.87516,132.721,2.02876],
Sm:[24.0042,2.47274,19.4258,0.196451,13.4396,14.3996,2.89604,128.007,2.20963],
Eu:[24.6274,2.3879,19.0886,0.1942,13.7603,13.7546,2.9227,123.174,2.5745],
Gd:[25.0709,2.25341,19.0798,0.181951,13.8518,12.9331,3.54545,101.398,2.4196],
Tb:[25.8976,2.24256,18.2185,0.196143,14.3167,12.6648,2.95354,115.362,3.58324],
Dy:[26.507,2.1802,17.6383,0.202172,14.5596,12.1899,2.96577,111.874,4.29728],
Ho:[26.9049,2.07051,17.294,0.19794,14.5583,11.4407,3.63837,92.6566,4.56796],
Er:[27.6563,2.07356,16.4285,0.223545,14.9779,11.3604,2.98233,105.703,5.92046],
Tm:[28.1819,2.02859,15.8851,0.238849,15.1542,10.9975,2.98706,102.961,6.75621],
Yb:[28.6641,1.9889,15.4345,0.257119,15.3087,10.6647,2.98963,100.417,7.56672],
Lu:[28.9476,1.90182,15.2208,9.98519,15.1,0.261033,3.71601,84.3298,7.97628],
Hf:[29.144,1.83262,15.1726,9.5999,14.7586,0.275116,4.30013,72.029,8.58154],
Ta:[29.2024,1.77333,15.2293,9.37046,14.5135,0.295977,4.76492,63.3644,9.24354],
W:[29.0818,1.72029,15.43,9.2259,14.4327,0.321703,5.11982,57.056,9.8875],
Re:[28.7621,1.67191,15.7189,9.09227,14.5564,0.3505,5.44174,52.0861,10.472],
Os:[28.1894,1.62903,16.155,8.97948,14.9305,0.382661,5.67589,48.1647,11.0005],
Ir:[27.3049,1.59279,16.7296,8.86553,15.6115,0.417916,5.83377,45.0011,11.4722],
Pt:[27.0059,1.51293,17.7639,8.81174,15.7131,0.424593,5.7837,38.6103,11.6883],
Au:[16.8819,0.4611,18.5913,8.6216,25.5582,1.4826,5.86,36.3956,12.0658],
Hg:[20.6809,0.545,19.0417,8.4484,21.6575,1.5729,5.9676,38.3246,12.6089],
Tl:[27.5446,0.65515,19.1584,8.70751,15.538,1.96347,5.52593,45.8149,13.1746],
Pb:[31.0617,0.6902,13.0637,2.3576,18.442,8.618,5.9696,47.2579,13.4118],
Bi:[33.3689,0.704,12.951,2.9238,16.5877,8.7937,6.4692,48.0093,13.5782],
Po:[34.6726,0.700999,15.4733,3.55078,13.1138,9.55642,7.02588,47.0045,13.677],
At:[35.3163,0.68587,19.0211,3.97458,9.49887,11.3824,7.42518,45.4715,13.7108],
Rn:[35.5631,0.6631,21.2816,4.0691,8.0037,14.0422,7.4433,44.2473,13.6905],
Fr:[35.9299,0.646453,23.0547,4.17619,12.1439,23.1052,2.11253,150.645,13.7247],
Ra:[35.763,0.616341,22.9064,3.87135,12.4739,19.9887,3.21097,142.325,13.6211],
Ac:[35.6597,0.589092,23.1032,3.65155,12.5977,18.599,4.08655,117.02,13.5266],
Th:[35.5645,0.563359,23.4219,3.46204,12.7473,17.8309,4.80703,99.1722,13.4314],
Pa:[35.8847,0.547751,23.2948,3.41519,14.1891,16.9235,4.17287,105.251,13.4287],
U:[36.0228,0.5293,23.4128,3.3253,14.9491,16.0927,4.188,100.613,13.3966],
Np:[36.1874,0.511929,23.5964,3.25396,15.6402,15.3622,4.1855,97.4908,13.3573],
Pu:[36.5254,0.499384,23.8083,3.26371,16.7707,14.9455,3.47947,105.98,13.3812],
Am:[36.6706,0.483629,24.0992,3.20647,17.3415,14.3136,3.49331,102.273,13.3592],
Cm:[36.6488,0.465154,24.4096,3.08997,17.399,13.4346,4.21665,88.4834,13.2887],
Bk:[36.7881,0.451018,24.7736,3.04619,17.8919,12.8946,4.23284,86.003,13.2754],
Cf:[36.9185,0.437533,25.1995,3.00775,18.3317,12.4044,4.24391,83.7881,13.2674]
};
// lambda is the value the panel prints (the Kalpha1 line, or the weighted
// mean of the doublet); k1/k2 are the two lines a Kalpha source really
// emits. A doublet entry profiles both at their own angles with the tube's
// 2:1 intensity ratio, so high-angle lines split as the measurement does.
var XRD_RADIATIONS=[
    {id:'CuKa1',label:'Cu Kα1',lambda:1.54056},
    {id:'CuKa',label:'Cu Kα',lambda:1.54184,k1:1.54056,k2:1.54439},
    {id:'MoKa1',label:'Mo Kα1',lambda:0.70930},
    {id:'MoKa',label:'Mo Kα',lambda:0.71073,k1:0.70930,k2:0.71359},
    {id:'CoKa1',label:'Co Kα1',lambda:1.78897},
    {id:'FeKa1',label:'Fe Kα1',lambda:1.93604},
    {id:'CrKa1',label:'Cr Kα1',lambda:2.28970},
    {id:'AgKa1',label:'Ag Kα1',lambda:0.55941}
];
// I(Kalpha2)/I(Kalpha1) for a laboratory tube: the Kalpha2 line carries half
// the photons, so a resolved pair shows the strong and weak components in 2:1.
var XRD_K2_RATIO=0.5;
var XRD_TWO_THETA_CHOICES=[60,90,120,150];
var XRD_2THETA_MIN=2;
var XRD_MAX_INDEX=40;
var XRD_MAX_REFLECTIONS=20000;
var XRD_STEP=0.02;
var XRD_LIST_MAX=30;
// Caglioti FWHM^2 = U tan^2(theta) + V tan(theta) + W, in degrees squared; the
// values describe a routine laboratory powder diffractometer and keep the
// variance positive over 0-180 degrees.
var XRD_CAGLIOTI_U=0.004;
var XRD_CAGLIOTI_V=-0.002;
var XRD_CAGLIOTI_W=0.008;
var xrdOpen=false;
var xrdRadiation='CuKa1';
var xrdTwoThetaMax=90;
var xrdResult=null;
// The experimental pattern is a session object: points, the name of the file
// they came from, the peak list derived from them, and the normalisation
// denominator, all thrown away when the webview reloads. The calculated curve
// is normalised to its own maximum and the experimental one to its own, so
// overlaying one never rescales the other.
var xrdExp=null;
var xrdSelection=null;   // the clicked calculated peak, and its experimental match
var xrdPlot=null;        // geometry of the last draw, so a click can be mapped back
var XRD_EXP_MATCH_WINDOW=1.0;   // deg: an experimental peak further away is not a match
var XRD_CLICK_WINDOW=1.0;       // deg: clicking between peaks selects nothing

function xrdFormFactor(el,s){
    var ff=XRD_FORM_FACTORS[el];
    if(!ff)return null;
    var s2=s*s;
    return ff[0]*Math.exp(-ff[1]*s2)+ff[2]*Math.exp(-ff[3]*s2)+
        ff[4]*Math.exp(-ff[5]*s2)+ff[6]*Math.exp(-ff[7]*s2)+ff[8];
}
// The Debye-Waller factor of one site at (h k l), from its Cartesian tensor.
// The thermal average of exp(i q.u) is exp(-2 pi^2 q U q), with q the
// scattering vector in 1/Angstrom: q = h a* + k b* + l c*. A tensor that is
// not positive definite along q would make the reflection brighter, which no
// vibration can do, so that direction is held at one; a site with no tensor
// has no stated vibration and keeps a factor of one as well.
function xrdDwb(u,h,k,l,rec){
    if(!u)return 1;
    var qx=h*rec[0].x+k*rec[1].x+l*rec[2].x;
    var qy=h*rec[0].y+k*rec[1].y+l*rec[2].y;
    var qz=h*rec[0].z+k*rec[1].z+l*rec[2].z;
    var e=u[0]*qx*qx+u[1]*qy*qy+u[2]*qz*qz+
        2*(u[3]*qx*qy+u[4]*qx*qz+u[5]*qy*qz);
    if(!(e>0))return 1;
    return Math.exp(-2*Math.PI*Math.PI*e);
}
function xrdRadiationById(id){
    for(var i=0;i<XRD_RADIATIONS.length;i++){
        if(XRD_RADIATIONS[i].id===id)return XRD_RADIATIONS[i];
    }
    return XRD_RADIATIONS[0];
}
// Lorentz-polarisation for an unpolarised laboratory beam, in the powder form.
function xrdLp(twoDeg){
    var t=twoDeg*Math.PI/180;
    var sh=Math.sin(t/2),ch=Math.cos(t/2);
    if(!(sh>1e-8)||!(ch>1e-8))return 0;
    return (1+Math.cos(t)*Math.cos(t))/(sh*sh*ch);
}
function xrdCagliotiFwhm(twoDeg){
    var t=Math.tan(twoDeg*Math.PI/360);
    var f2=XRD_CAGLIOTI_U*t*t+XRD_CAGLIOTI_V*t+XRD_CAGLIOTI_W;
    return f2>0?Math.sqrt(f2):0;
}
function xrdReciprocalForm(rec){
    function dot(p,q){return p.x*q.x+p.y*q.y+p.z*q.z}
    return {
        aa:dot(rec[0],rec[0]),bb:dot(rec[1],rec[1]),cc:dot(rec[2],rec[2]),
        ab:2*dot(rec[0],rec[1]),ac:2*dot(rec[0],rec[2]),bc:2*dot(rec[1],rec[2])
    };
}
// The atoms whose electrons scatter: the drawn cell, reduced by the disorder
// filter and without ghosts, each with its fractional coordinates.
function xrdCellSites(){
    if(!CRY||!CRY.baseAtoms||!CRY_INV)return null;
    var sites=[],missing={},i;
    for(i=0;i<CRY.baseAtoms.length;i++){
        var ba=CRY.baseAtoms[i];
        if(isGhostAtom(ba))continue;
        if(!neighBaseVisible(ba))continue;
        var el=ba.element==='D'?'H':ba.element;
        if(!XRD_FORM_FACTORS[el]){missing[el]=true;continue}
        // The same occupancy the formula, the density and every save read: a
        // selected disorder group at full occupancy weighs 1 in the pattern
        // too, not the file's 0.5.
        var occ=crystalOccOf(ba);
        if(!(occ>0))continue;
        var f=neighFractional(ba);
        // The file's displacement tensor, when it has a usable one: the six
        // numbers are Cartesian Angstrom squared, the frame the structure
        // factor wants. A tensor with a non-positive trace cannot describe a
        // vibration and is left out rather than guessed at.
        var u=ba.uAniso,uu=null;
        if(u&&u.length>=6){
            var okU=true,q;
            for(q=0;q<6;q++)if(typeof u[q]!=='number'||!isFinite(u[q]))okU=false;
            if(okU&&u[0]+u[1]+u[2]>0)uu=[u[0],u[1],u[2],u[3],u[4],u[5]];
        }
        sites.push({el:el,f:[f[0],f[1],f[2]],occ:occ,u:uu});
    }
    var miss=[];
    for(var k in missing){
        if(Object.prototype.hasOwnProperty.call(missing,k))miss.push(k);
    }
    miss.sort();
    return {sites:sites,missing:miss};
}
// The pattern is built as a resumable job so a large cell cannot freeze the
// page: enumeration, structure factors, peak merging and profile broadening
// each remember where they stopped, and xrdPatternStep runs until the slice's
// deadline. The sums and the order of the stages are the ones the pattern
// always had; only the chunking is new. xrdPattern() drives the same job to
// completion in one call, which is what the engine tests read.
function xrdProfileAdd(job,two,inten){
    var fwhm=xrdCagliotiFwhm(two);
    if(!(fwhm>0.005))fwhm=0.005;
    if(fwhm>4)fwhm=4;
    var half=2.5*fwhm;
    var n=job.n,step=XRD_STEP,min=XRD_2THETA_MIN,prof=job.prof;
    var i0=Math.max(0,Math.ceil((two-half-min)/step));
    var i1=Math.min(n-1,Math.floor((two+half-min)/step));
    var k4=4*Math.LN2;
    for(var q=i0;q<=i1;q++){
        var dx=(min+q*step-two)/fwhm;
        prof[q]+=inten*Math.exp(-k4*dx*dx);
    }
}
// Pure: sites + lattice + wavelength + limit in, a job out. The index bounds
// are the dual relation |h| = |g.a| <= |g||a|, so a reflection inside the
// reciprocal sphere can never fall outside the loops. A second wavelength
// (lambda2) profiles the Kalpha2 line of every peak beside its Kalpha1 line,
// at the tube's 2:1 ratio; without it the pattern is the single line.
function xrdPatternStart(sites,lv,lambda,twoMax,lambda2){
    if(!sites||sites.length===0||!(lambda>0)||!(twoMax>XRD_2THETA_MIN))return null;
    if(!(lambda2>0)||lambda2===lambda)lambda2=0;
    var rec=hklReciprocals(lv);
    if(!rec)return null;
    var form=xrdReciprocalForm(rec);
    var gmax=2*Math.sin(twoMax*Math.PI/360)/lambda;
    function len(v){return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])}
    var indexCapped=len(lv[0])*gmax>XRD_MAX_INDEX||len(lv[1])*gmax>XRD_MAX_INDEX||len(lv[2])*gmax>XRD_MAX_INDEX;
    var limA=Math.min(XRD_MAX_INDEX,Math.ceil(len(lv[0])*gmax)+1);
    var limB=Math.min(XRD_MAX_INDEX,Math.ceil(len(lv[1])*gmax)+1);
    var limC=Math.min(XRD_MAX_INDEX,Math.ceil(len(lv[2])*gmax)+1);
    return {
        sites:sites,rec:rec,form:form,lambda:lambda,lambda2:lambda2,twoMax:twoMax,
        indexCapped:indexCapped,limA:limA,limB:limB,limC:limC,
        phase:'enum',h:-limA,k:-limB,list:[],
        found:0,capped:false,i:0,peaks:[],
        peakIndex:0,n:0,prof:null
    };
}
function xrdPatternStep(job,deadline){
    for(;;){
        if(job.phase==='enum'){
            var form=job.form;
            while(job.h<=job.limA){
                var h=job.h;
                while(job.k<=job.limB){
                    var k=job.k;
                    for(var l=-job.limC;l<=job.limC;l++){
                        var g2=form.aa*h*h+form.bb*k*k+form.cc*l*l+form.ab*h*k+form.ac*h*l+form.bc*k*l;
                        if(!(g2>1e-12))continue;
                        var g=Math.sqrt(g2);
                        var sinT=job.lambda*g/2;
                        if(sinT>=1)continue;
                        var two=2*Math.asin(sinT)*180/Math.PI;
                        if(two<XRD_2THETA_MIN||two>job.twoMax+1e-9)continue;
                        job.list.push({two:two,h:h,k:k,l:l,s:g/2});
                    }
                    job.k++;
                    if(deadline!==Infinity&&Date.now()>=deadline)return false;
                }
                job.k=-job.limB;
                job.h++;
                if(deadline!==Infinity&&Date.now()>=deadline)return false;
            }
            job.phase='sort';
            continue;
        }
        if(job.phase==='sort'){
            job.list.sort(function(p,q){
                if(p.two!==q.two)return p.two-q.two;
                if(p.h!==q.h)return p.h-q.h;
                if(p.k!==q.k)return p.k-q.k;
                return p.l-q.l;
            });
            job.found=job.list.length;
            job.capped=job.found>XRD_MAX_REFLECTIONS;
            if(job.capped)job.list.length=XRD_MAX_REFLECTIONS;
            job.phase='sf';
            continue;
        }
        if(job.phase==='sf'){
            var list=job.list,sites=job.sites,rec=job.rec,peaks=job.peaks;
            while(job.i<list.length){
                var r=list[job.i],re=0,im=0;
                for(var j=0;j<sites.length;j++){
                    var st=sites[j];
                    var f=xrdFormFactor(st.el,r.s)*st.occ;
                    if(!f)continue;
                    if(st.u)f*=xrdDwb(st.u,r.h,r.k,r.l,rec);
                    var ph=2*Math.PI*(r.h*st.f[0]+r.k*st.f[1]+r.l*st.f[2]);
                    re+=f*Math.cos(ph);im+=f*Math.sin(ph);
                }
                var inten=(re*re+im*im)*xrdLp(r.two);
                var last=peaks.length>0?peaks[peaks.length-1]:null;
                if(last&&Math.abs(r.two-last.two)<=1e-4){
                    last.I+=inten;last.n++;
                    if(inten>last.rawI){last.rawI=inten;last.hkl=[r.h,r.k,r.l]}
                }else{
                    peaks.push({two:r.two,d:1/(2*r.s),hkl:[r.h,r.k,r.l],I:inten,rawI:inten,n:1});
                }
                job.i++;
                if(deadline!==Infinity&&Date.now()>=deadline)return false;
            }
            job.phase='keep';
            continue;
        }
        if(job.phase==='keep'){
            // A systematic absence is zero only in exact arithmetic: the
            // structure factor of a forbidden reflection leaves rounding noise
            // (about 1e-30 of the strongest line), so the cut is relative to
            // the pattern's own peak.
            var pk=job.peaks,imax=0,i;
            for(i=0;i<pk.length;i++)if(pk[i].I>imax)imax=pk[i].I;
            if(!(imax>0))imax=1;
            var keep=[];
            for(i=0;i<pk.length;i++){
                if(!(pk[i].I>imax*1e-12))continue;
                // (h k l) and (-h -k -l) are one plane: label the line by the
                // reading whose first nonzero index is positive, which is how
                // it is written.
                var hh=pk[i].hkl;
                if(hh[0]<0||(hh[0]===0&&(hh[1]<0||(hh[1]===0&&hh[2]<0))))pk[i].hkl=[-hh[0],-hh[1],-hh[2]];
                keep.push(pk[i]);
            }
            job.peaks=keep;
            for(i=0;i<keep.length;i++)keep[i].rel=keep[i].I/imax*100;
            job.n=Math.floor((job.twoMax-XRD_2THETA_MIN)/XRD_STEP+1e-9)+1;
            if(job.n<2)job.n=2;
            job.prof=new Float64Array(job.n);
            job.peakIndex=0;
            job.phase='profile';
            continue;
        }
        if(job.phase==='profile'){
            var w2=job.lambda2?XRD_K2_RATIO/(1+XRD_K2_RATIO):0;
            var w1=1-w2;
            while(job.peakIndex<job.peaks.length){
                var p=job.peaks[job.peakIndex];
                xrdProfileAdd(job,p.two,p.I*w1);
                if(w2>0){
                    var sin2=job.lambda2/(2*p.d);
                    if(sin2<1){
                        var two2=2*Math.asin(sin2)*180/Math.PI;
                        if(two2>=XRD_2THETA_MIN&&two2<=job.twoMax+1e-9)xrdProfileAdd(job,two2,p.I*w2);
                    }
                }
                job.peakIndex++;
                if(deadline!==Infinity&&Date.now()>=deadline)return false;
            }
            job.phase='normalize';
            continue;
        }
        if(job.phase==='normalize'){
            var kk,pmax=0;
            for(kk=0;kk<job.n;kk++)if(job.prof[kk]>pmax)pmax=job.prof[kk];
            if(!(pmax>0))pmax=1;
            for(kk=0;kk<job.n;kk++)job.prof[kk]=job.prof[kk]/pmax*100;
            job.phase='done';
            return true;
        }
        return true;
    }
}
function xrdPatternResult(job){
    return {
        atoms:job.sites.length,peaks:job.peaks,reflections:job.list.length,found:job.found,
        capped:job.capped,indexCapped:job.indexCapped,doublet:job.lambda2>0,
        twoMax:job.twoMax,min:XRD_2THETA_MIN,step:XRD_STEP,values:job.prof
    };
}
function xrdPattern(sites,lv,lambda,twoMax,lambda2){
    var job=xrdPatternStart(sites,lv,lambda,twoMax,lambda2);
    if(!job)return null;
    while(!xrdPatternStep(job,Infinity)){}
    return xrdPatternResult(job);
}
// How far along the slice job is, for the panel's summary line while it runs.
function xrdPatternProgress(job){
    if(job.phase==='enum')return 3;
    if(job.phase==='sort')return 10;
    if(job.phase==='sf'){
        var n=job.list.length||1;
        return 10+Math.floor(job.i/n*70);
    }
    if(job.phase==='keep')return 82;
    if(job.phase==='profile'){
        var pn=job.peaks.length||1;
        return 83+Math.floor(job.peakIndex/pn*16);
    }
    return 99;
}
function xrdJobProgressText(job){
    var txt='Calculating the powder pattern... '+xrdPatternProgress(job)+'%';
    if(job.phase==='sf')txt+=' ('+job.i+' of '+job.list.length+' reflections)';
    else if(job.phase==='profile')txt+=' ('+job.peakIndex+' of '+job.peaks.length+' lines)';
    return txt;
}
// The panel's reading of the cell: the sites, the radiation and a started job.
function xrdPrepareJob(){
    if(!CRY||!CRY.latticeVectors)return {ok:false,message:'Not a crystal structure.'};
    var cell=xrdCellSites();
    if(!cell)return {ok:false,message:'Not a crystal structure.'};
    if(cell.sites.length===0)return {ok:false,message:'No atom in the cell has a scattering factor.',missing:cell.missing};
    var rad=xrdRadiationById(xrdRadiation);
    // A Kalpha entry carries both lines of the doublet; the primary wavelength
    // places the peaks (the table lists the Kalpha1 reading) and the secondary
    // one profiles the weak line beside each of them.
    var job=xrdPatternStart(cell.sites,CRY.latticeVectors,rad.k1||rad.lambda,xrdTwoThetaMax,rad.k2||0);
    if(!job)return {ok:false,message:'The cell has no volume, so no lattice can be built.',missing:cell.missing};
    job.rad=rad;
    job.missing=cell.missing;
    return {ok:true,rad:rad,job:job,missing:cell.missing};
}
function xrdCompute(){
    var prep=xrdPrepareJob();
    if(!prep.ok)return {ok:false,message:prep.message,missing:prep.missing};
    var job=prep.job;
    while(!xrdPatternStep(job,Infinity)){}
    return {ok:true,rad:prep.rad,pattern:xrdPatternResult(job),missing:prep.missing};
}
// The measured pattern is a plain two-column text file: 2theta then intensity,
// one sample a line. Comment lines (#, ;, //), comma separators, extra columns
// and header words are tolerated. Both curves are normalised to their own
// maximum, so importing or clearing a measurement never rescales the
// calculated pattern; the two share only the 2theta axis.
function xrdParseXY(text){
    var lines=String(text||'').split(/\\r?\\n/);
    var pts=[],i;
    for(i=0;i<lines.length;i++){
        var line=lines[i];
        var hash=line.indexOf('#');
        if(hash>=0)line=line.substring(0,hash);
        line=line.replace(/\\/\\/.*$/,'').replace(/,/g,' ').trim();
        if(!line||line.charAt(0)===';'||line.charAt(0)==='!')continue;
        var f=line.split(/\\s+/);
        if(f.length<2)continue;
        var x=Number(f[0]),y=Number(f[1]);
        if(!isFinite(x)||!isFinite(y)||x<0||x>180)continue;
        pts.push([x,y]);
    }
    if(pts.length<3)return null;
    pts.sort(function(a,b){return a[0]-b[0]});
    var clean=[];
    for(i=0;i<pts.length;i++){
        var last=clean.length>0?clean[clean.length-1]:null;
        if(last&&pts[i][0]-last[0]<1e-9){
            if(pts[i][1]>last[1])last[1]=pts[i][1];
            continue;
        }
        clean.push([pts[i][0],pts[i][1]]);
    }
    var imax=-Infinity;
    for(i=0;i<clean.length;i++)if(clean[i][1]>imax)imax=clean[i][1];
    if(!(imax>0))return null;
    var peaks=[];
    for(i=1;i<clean.length-1;i++){
        if(!(clean[i][1]>=clean[i-1][1])||!(clean[i][1]>clean[i+1][1]))continue;
        var v=clean[i][1]/imax*100;
        if(v<2)continue;
        var prev=peaks.length>0?peaks[peaks.length-1]:null;
        if(prev&&clean[i][0]-prev.two<=0.05){
            if(v>prev.I){prev.two=clean[i][0];prev.I=v}
            continue;
        }
        peaks.push({two:clean[i][0],I:v});
    }
    return {points:clean,peaks:peaks,max:imax,from:clean[0][0],to:clean[clean.length-1][0]};
}
function xrdNearestExpPeak(two){
    if(!xrdExp||!xrdExp.peaks.length)return {peak:null,distance:Infinity};
    var best=null,bd=Infinity,i;
    for(i=0;i<xrdExp.peaks.length;i++){
        var d=Math.abs(xrdExp.peaks[i].two-two);
        if(d<bd){bd=d;best=xrdExp.peaks[i]}
    }
    return {peak:bd<=XRD_EXP_MATCH_WINDOW?best:null,distance:bd};
}
function xrdPeakInfoText(){
    var el=document.getElementById('xrd-peak-info');
    if(!el)return;
    if(!xrdSelection){
        el.textContent=xrdExp?
            'Click a calculated peak or a row to read the nearest experimental peak.':
            'Click a calculated peak or a row to mark it; import an experimental .xy to compare positions.';
        return;
    }
    var p=xrdSelection;
    var text='('+p.hkl.join(' ')+') at 2theta = '+p.two.toFixed(3)+' deg, I = '+p.rel.toFixed(1)+' %: ';
    if(!xrdExp){el.textContent=text+'no experimental pattern loaded.';return}
    var near=xrdNearestExpPeak(p.two);
    if(!isFinite(near.distance)){
        text+='the experimental curve has no detected peak to compare with.';
    }else if(near.peak){
        var d=p.two-near.peak.two;
        text+='nearest experimental peak at '+near.peak.two.toFixed(3)+' deg, delta 2theta (calc - exp) = '+
            (d>=0?'+':'')+d.toFixed(3)+' deg.';
        // The calculated position is the Kalpha1 line; a measured maximum is
        // the Kalpha1/Kalpha2 blend at the line's own width, so part of the
        // delta is systematic and grows with 2theta.
        if(xrdResult&&xrdResult.ok&&xrdResult.pattern&&xrdResult.pattern.doublet){
            text+=' The calculated value is the Kα1 line while the experimental maximum is the Kα1/Kα2 blend, '+
                'so part of the delta is systematic (it grows with 2theta and depends on the line width).';
        }
    }else{
        text+='no experimental peak within '+XRD_EXP_MATCH_WINDOW.toFixed(1)+' deg (nearest at '+
            near.distance.toFixed(2)+' deg away).';
    }
    el.textContent=text;
}
function xrdExpInfo(){
    var el=document.getElementById('xrd-exp-info');
    var clearBtn=document.getElementById('xrd-clear-exp');
    if(clearBtn)clearBtn.disabled=!xrdExp;
    if(!el)return;
    if(!xrdExp){
        el.textContent='No experimental pattern loaded: each curve is normalised to its own maximum.';
        return;
    }
    el.textContent='Experimental (grey): '+xrdExp.name+' - '+xrdExp.points.length+' points, '+
        xrdExp.from.toFixed(2)+'-'+xrdExp.to.toFixed(2)+' deg, '+xrdExp.peaks.length+' peak'+
        (xrdExp.peaks.length===1?'':'s')+', normalised to its own maximum.';
}
function xrdMarkSelectedRow(){
    var listEl=document.getElementById('xrd-list');
    if(!listEl)return;
    var rows=listEl.getElementsByClassName('xr-row');
    for(var i=0;i<rows.length;i++){
        var t=parseFloat(rows[i].getAttribute('data-two'));
        var hit=xrdSelection&&Math.abs(t-xrdSelection.two)<1e-9;
        rows[i].className=hit?'xr-row sel':'xr-row';
    }
}
function xrdSelectPeak(peak){
    xrdSelection=peak;
    xrdPeakInfoText();
    xrdMarkSelectedRow();
    xrdDraw(xrdResult&&xrdResult.ok?xrdResult.pattern:null,
        xrdResult&&!xrdResult.ok?(xrdResult.message||''):'');
}
function xrdLoadExperimental(text,name){
    var exp=xrdParseXY(text);
    if(!exp){
        var infoEl=document.getElementById('xrd-exp-info');
        if(infoEl)infoEl.textContent='Could not read a two-column 2theta/intensity pattern from '+name+'.'+
            (xrdExp?' Keeping '+xrdExp.name+'.':'');
        return false;
    }
    exp.name=name;
    xrdExp=exp;
    xrdSelection=null;
    xrdUpdate();
    layoutPanels();
    return true;
}
function xrdClearExperimental(){
    if(!xrdExp)return;
    xrdExp=null;
    xrdSelection=null;
    xrdUpdate();
    layoutPanels();
}
function xrdCanvasClick(ev){
    if(!xrdPlot||!xrdResult||!xrdResult.ok)return;
    var cv=document.getElementById('xrd-canvas');
    if(!cv)return;
    var rect=cv.getBoundingClientRect();
    var px=ev.clientX-rect.left,py=ev.clientY-rect.top;
    var pl=xrdPlot;
    if(px<pl.L||px>pl.L+pl.pw||py<pl.T||py>pl.T+pl.ph)return;
    var two=pl.twoMin+(px-pl.L)/pl.pw*(pl.twoMax-pl.twoMin);
    var peaks=xrdResult.pattern.peaks,best=null,bd=Infinity,i;
    for(i=0;i<peaks.length;i++){
        var d=Math.abs(peaks[i].two-two);
        if(d<bd){bd=d;best=peaks[i]}
    }
    if(!best||bd>XRD_CLICK_WINDOW){
        xrdSelection=null;
        xrdMarkSelectedRow();
        var info=document.getElementById('xrd-peak-info');
        if(info)info.textContent='No calculated peak within '+XRD_CLICK_WINDOW.toFixed(1)+
            ' deg of 2theta = '+two.toFixed(2)+' deg.';
        xrdDraw(xrdResult.pattern,'');
        return;
    }
    xrdSelectPeak(best);
}
// The two overlays, both on the plot rectangle kept by the last draw.
function xrdDrawExperimental(g,plot){
    if(!xrdExp||xrdExp.points.length<2)return;
    var pts=xrdExp.points,span=plot.twoMax-plot.twoMin,i;
    g.save();
    g.beginPath();g.rect(plot.L,plot.T,plot.pw,plot.ph);g.clip();
    g.beginPath();
    for(i=0;i<pts.length;i++){
        var v=pts[i][1]/xrdExp.max*100;
        if(v<0)v=0;else if(v>100)v=100;
        var x=plot.L+(pts[i][0]-plot.twoMin)/span*plot.pw;
        var y=plot.T+plot.ph*(1-v/100);
        if(i===0)g.moveTo(x,y);else g.lineTo(x,y);
    }
    g.strokeStyle='rgba(205,205,205,0.9)';
    g.lineWidth=1.3;
    g.stroke();
    g.restore();
}
function xrdDrawMarkers(g,plot){
    if(!xrdSelection)return;
    var span=plot.twoMax-plot.twoMin;
    var x=plot.L+(xrdSelection.two-plot.twoMin)/span*plot.pw;
    if(x<plot.L||x>plot.L+plot.pw)return;
    var y=plot.T+plot.ph*(1-Math.max(0,Math.min(100,xrdSelection.rel))/100);
    g.save();
    g.setLineDash([3,3]);
    g.lineWidth=1;
    g.strokeStyle='#4ec9b0';
    g.beginPath();g.moveTo(x,plot.T);g.lineTo(x,plot.T+plot.ph);g.stroke();
    g.setLineDash([]);
    g.fillStyle='#4ec9b0';
    g.beginPath();g.arc(x,y,3,0,2*Math.PI);g.fill();
    var near=xrdExp?xrdNearestExpPeak(xrdSelection.two):null;
    if(near&&near.peak){
        var xe=plot.L+(near.peak.two-plot.twoMin)/span*plot.pw;
        if(xe>=plot.L&&xe<=plot.L+plot.pw){
            g.strokeStyle='#b0b0b0';
            g.setLineDash([3,3]);
            g.beginPath();g.moveTo(xe,plot.T);g.lineTo(xe,plot.T+plot.ph);g.stroke();
        }
    }
    g.restore();
}
function xrdDraw(pat,message){
    var cv=document.getElementById('xrd-canvas');
    if(!cv)return;
    var cw=Math.max(220,Math.round(cv.clientWidth||440));
    var ch=150;
    var dpr=Math.min(2,window.devicePixelRatio||1);
    cv.width=Math.round(cw*dpr);cv.height=Math.round(ch*dpr);
    var g=cv.getContext('2d');
    if(!g)return;
    g.setTransform(dpr,0,0,dpr,0,0);
    g.fillStyle='#0d1117';g.fillRect(0,0,cw,ch);
    var L=38,R=10,T=10,B=22;
    var pw=cw-L-R,ph=ch-T-B;
    var twoMax=pat?pat.twoMax:xrdTwoThetaMax;
    var twoMin=pat?pat.min:XRD_2THETA_MIN;
    g.lineWidth=1;
    g.strokeStyle='rgba(128,128,128,0.25)';
    var gi;
    for(gi=0;gi<=4;gi++){
        var yy=Math.round(T+ph*gi/4)+0.5;
        g.beginPath();g.moveTo(L,yy);g.lineTo(L+pw,yy);g.stroke();
    }
    var span=twoMax-twoMin;
    // The geometry a click maps back through: the same rectangle every overlay
    // and the grid use, rebuilt on each draw because the panel can be resized.
    xrdPlot={L:L,T:T,pw:pw,ph:ph,twoMin:twoMin,twoMax:twoMax};
    var tickStep=span>60?20:(span>30?10:5);
    g.fillStyle='#8b949e';
    g.font='10px sans-serif';
    g.textBaseline='middle';g.textAlign='right';
    g.fillText('100',L-4,T+2);g.fillText('50',L-4,T+ph/2);g.fillText('0',L-4,T+ph-2);
    g.textAlign='center';g.textBaseline='top';
    for(var tt=Math.ceil(twoMin/tickStep)*tickStep;tt<=twoMax+1e-9;tt+=tickStep){
        var xx=Math.round(L+(tt-twoMin)/span*pw)+0.5;
        g.beginPath();g.moveTo(xx,T);g.lineTo(xx,T+ph);g.stroke();
        g.fillText(String(tt),xx,T+ph+4);
    }
    g.textAlign='left';g.fillText('I (%)',L+3,T+1);
    g.textAlign='right';g.fillText('2theta (deg)',L+pw,T+1);
    if(!pat||!pat.values||pat.values.length<2){
        xrdDrawExperimental(g,xrdPlot);
        if(message){
            g.fillStyle='#8b949e';g.textAlign='center';g.textBaseline='middle';
            g.font='11px sans-serif';
            g.fillText(message,L+pw/2,T+ph/2);
        }
        return;
    }
    var vals=pat.values,n=vals.length,i;
    function xOf(i){return L+i/(n-1)*pw}
    function yOf(v){return T+ph*(1-v/100)}
    g.beginPath();g.moveTo(xOf(0),yOf(0));
    for(i=0;i<n;i++)g.lineTo(xOf(i),yOf(vals[i]));
    g.lineTo(xOf(n-1),yOf(0));g.closePath();
    g.fillStyle='rgba(78,201,176,0.16)';g.fill();
    g.beginPath();g.moveTo(xOf(0),yOf(vals[0]));
    for(i=1;i<n;i++)g.lineTo(xOf(i),yOf(vals[i]));
    g.strokeStyle='#4ec9b0';g.lineWidth=1.2;g.stroke();
    xrdDrawExperimental(g,xrdPlot);
    xrdDrawMarkers(g,xrdPlot);
}
// The panel's calculation is one job at a time: a new structure, radiation or
// range drops the old job (its next tick sees a stale sequence and stops), and
// a job that takes longer than a slice reports its progress and keeps the page
// responsive. The Recalculate button becomes Cancel while it runs.
var xrdJob=null,xrdJobSeq=0;
var XRD_SLICE_MS=14;

function xrdAbortJob(){
    if(!xrdJob)return false;
    xrdJob=null;
    xrdJobSeq++;
    return true;
}
function xrdSetRunButton(label,title){
    var btn=document.getElementById('xrd-refresh');
    if(btn){btn.textContent=label;btn.title=title}
}
function xrdRenderError(res){
    var sumEl=document.getElementById('xrd-summary');
    var noteEl=document.getElementById('xrd-note');
    var listEl=document.getElementById('xrd-list');
    var expBtn=document.getElementById('xrd-export');
    if(sumEl)sumEl.textContent=res.message||'';
    var nm=[];
    if(res.missing&&res.missing.length)nm.push('No scattering factor for '+res.missing.join(', ')+' in the cell.');
    if(noteEl)noteEl.textContent=nm.join(' ');
    if(listEl)listEl.innerHTML='<div class="xr-empty">No pattern.</div>';
    if(expBtn)expBtn.disabled=true;
    xrdPeakInfoText();
    xrdExpInfo();
    xrdDraw(null,res.message||'');
}
function xrdRenderResult(res){
    var pat=res.pattern,rad=res.rad;
    var sumEl=document.getElementById('xrd-summary');
    var noteEl=document.getElementById('xrd-note');
    var listEl=document.getElementById('xrd-list');
    var expBtn=document.getElementById('xrd-export');
    if(sumEl){
        sumEl.textContent=rad.label+' '+rad.lambda.toFixed(4)+' A - '+pat.peaks.length+' peak'+
            (pat.peaks.length===1?'':'s')+' from '+pat.reflections+' reflection'+
            (pat.reflections===1?'':'s')+' ('+pat.min+'-'+pat.twoMax+' deg)';
    }
    var notes=[];
    if(pat.capped)notes.push('Only the '+XRD_MAX_REFLECTIONS+' lowest-angle reflections were kept (the cell is too large for the rest).');
    if(pat.indexCapped)notes.push('Reflections needing an index above '+XRD_MAX_INDEX+' are not enumerated.');
    if(pat.doublet)notes.push('The curve carries the Kα1/Kα2 doublet at its 2:1 ratio; the peak list is the Kα1 line, '+
        'so a 2theta read off an experimental maximum sits slightly higher at high angle (the Kα2 blend).');
    if(res.missing.length)notes.push('No scattering factor for '+res.missing.join(', ')+' - those atoms are not in the pattern.');
    if(noteEl)noteEl.textContent=notes.join(' ');
    if(expBtn)expBtn.disabled=false;
    var rows=pat.peaks.slice();
    rows.sort(function(a,b){return b.rel-a.rel});
    rows=rows.slice(0,XRD_LIST_MAX);
    rows.sort(function(a,b){return a.two-b.two});
    var html='<div class="xr-head"><span class="xr-c1">2theta</span><span class="xr-c2">d (A)</span><span class="xr-c3">(h k l)</span><span class="xr-c4">I %</span><span class="xr-c5">x</span></div>';
    for(var ri=0;ri<rows.length;ri++){
        var rp=rows[ri];
        html+='<div class="xr-row" data-two="'+rp.two+'" title="Click to compare this line with the experimental curve"><span class="xr-c1">'+rp.two.toFixed(3)+'</span>'+
            '<span class="xr-c2">'+rp.d.toFixed(4)+'</span>'+
            '<span class="xr-c3">('+rp.hkl.join(' ')+')</span>'+
            '<span class="xr-c4">'+rp.rel.toFixed(1)+'</span>'+
            '<span class="xr-c5">'+rp.n+'</span></div>';
    }
    if(rows.length===0)html+='<div class="xr-empty">No reflections in this range.</div>';
    if(listEl)listEl.innerHTML=html;
    xrdPeakInfoText();
    xrdExpInfo();
    xrdDraw(pat,'');
}
function xrdJobTick(seq){
    if(!xrdJob||seq!==xrdJobSeq)return;
    var job=xrdJob;
    var done;
    try{
        done=xrdPatternStep(job,Date.now()+XRD_SLICE_MS);
    }catch(e){
        xrdJob=null;
        xrdJobSeq++;
        xrdSetRunButton('Recalculate','Recalculate after the structure changed');
        xrdResult=null;
        xrdRenderError({message:'The powder calculation stopped: '+
            (e&&e.message?e.message:'unexpected error')});
        if(xrdOpen)layoutPanels();
        return;
    }
    if(!xrdJob||seq!==xrdJobSeq)return;
    if(!done){
        var sum=document.getElementById('xrd-summary');
        if(sum)sum.textContent=xrdJobProgressText(job);
        setTimeout(function(){xrdJobTick(seq)},0);
        return;
    }
    xrdJob=null;
    xrdSetRunButton('Recalculate','Recalculate after the structure changed');
    xrdResult={ok:true,rad:job.rad,pattern:xrdPatternResult(job),missing:job.missing};
    xrdRenderResult(xrdResult);
    if(xrdOpen)layoutPanels();
}
function xrdCancelCompute(){
    if(!xrdAbortJob())return;
    xrdSetRunButton('Recalculate','Recalculate after the structure changed');
    xrdResult=null;
    xrdSelection=null;
    var sumEl=document.getElementById('xrd-summary');
    var noteEl=document.getElementById('xrd-note');
    var listEl=document.getElementById('xrd-list');
    var expBtn=document.getElementById('xrd-export');
    if(sumEl)sumEl.textContent='Calculation cancelled.';
    if(noteEl)noteEl.textContent='';
    if(listEl)listEl.innerHTML='<div class="xr-empty">Calculation cancelled. Recalculate to start again.</div>';
    if(expBtn)expBtn.disabled=true;
    xrdPeakInfoText();
    xrdDraw(null,'Calculation cancelled.');
    if(xrdOpen)layoutPanels();
}
function xrdUpdate(){
    // A new pattern invalidates the marked line: the peak it pointed at may not
    // exist any more.
    xrdAbortJob();
    xrdSelection=null;
    xrdResult=null;
    var expBtn=document.getElementById('xrd-export');
    if(expBtn)expBtn.disabled=true;
    var prep=xrdPrepareJob();
    if(!prep.ok){
        xrdSetRunButton('Recalculate','Recalculate after the structure changed');
        xrdRenderError(prep);
        return;
    }
    var job=prep.job;
    xrdJob=job;
    var seq=++xrdJobSeq;
    xrdSetRunButton('Cancel','Stop the powder calculation');
    var sumEl=document.getElementById('xrd-summary');
    if(sumEl)sumEl.textContent=xrdJobProgressText(job);
    var noteEl=document.getElementById('xrd-note');
    if(noteEl)noteEl.textContent='';
    var listEl=document.getElementById('xrd-list');
    if(listEl)listEl.innerHTML='<div class="xr-empty">Calculating...</div>';
    xrdPeakInfoText();
    xrdExpInfo();
    xrdDraw(null,'Calculating the powder pattern...');
    // The first slice runs now, so a cell that fits in one slice is finished
    // before this function returns and the panel is complete for the layout.
    xrdJobTick(seq);
}
function xrdRefresh(){
    if(!xrdOpen)return;
    var radSel=document.getElementById('xrd-rad');
    if(radSel&&radSel.value)xrdRadiation=radSel.value;
    var maxSel=document.getElementById('xrd-2max');
    if(maxSel){
        var mv=parseInt(maxSel.value,10);
        if(mv>XRD_2THETA_MIN&&mv<=180)xrdTwoThetaMax=mv;
    }
    // A new radiation or range can add or drop a long note above the plot and
    // change the panel height; the layout is told about the new size.
    xrdUpdate();
    layoutPanels();
}
function xrdFillControls(){
    var rad=document.getElementById('xrd-rad');
    var i;
    if(rad&&rad.options.length===0){
        for(i=0;i<XRD_RADIATIONS.length;i++){
            var o=document.createElement('option');
            o.value=XRD_RADIATIONS[i].id;
            o.textContent=XRD_RADIATIONS[i].label+' ('+XRD_RADIATIONS[i].lambda.toFixed(4)+' A)';
            rad.appendChild(o);
        }
    }
    if(rad)rad.value=xrdRadiation;
    var mx=document.getElementById('xrd-2max');
    if(mx&&mx.options.length===0){
        for(i=0;i<XRD_TWO_THETA_CHOICES.length;i++){
            var o2=document.createElement('option');
            o2.value=String(XRD_TWO_THETA_CHOICES[i]);
            o2.textContent=XRD_TWO_THETA_CHOICES[i]+'\u00b0';
            mx.appendChild(o2);
        }
    }
    if(mx)mx.value=String(xrdTwoThetaMax);
}
function xrdOpenPanel(){
    var panel=document.getElementById('xrd-panel');
    if(!panel||!CRY)return;
    // The pattern is of the cell, not of the pair a comparison draws; the
    // right-hand structure has no pattern to show beside it.
    if(diffMode){
        modeInfoEl.textContent='Powder XRD is not available in diff mode';
        return;
    }
    xrdOpen=true;
    panelMinState.xrd=false;
    xrdFillControls();
    panel.style.display='block';
    // A small cell finishes its first slice before the layout measures the
    // panel; a large one reports progress and is laid out again when it ends.
    xrdUpdate();
    syncPanelChips();
    layoutPanels();
}
// The x folds the pattern into its chip; the crystal panel's button and the
// chip both open it again, with the choices it was left on.
function xrdMinimizePanel(){
    if(!xrdOpen)return;
    panelMinState.xrd=true;
    syncPanelChips();
    layoutPanels();
}
function xrdClosePanel(){
    // A running job has no panel to fill; stop it with the panel.
    if(xrdAbortJob())xrdSetRunButton('Recalculate','Recalculate after the structure changed');
    var panel=document.getElementById('xrd-panel');
    if(panel)panel.style.display='none';
    xrdOpen=false;
    panelMinState.xrd=false;
    syncPanelChips();
    layoutPanels();
}
function xrdExport(){
    if(!xrdResult||!xrdResult.ok||!xrdResult.pattern)return;
    var pat=xrdResult.pattern,rad=xrdResult.rad;
    var out=[];
    out.push('# Molecular Viewer - powder X-ray diffraction');
    out.push('# Radiation: '+rad.label+' (lambda = '+rad.lambda.toFixed(5)+' A)');
    if(rad.k2)out.push('# Doublet: Kα1 = '+rad.k1.toFixed(5)+' A, Kα2 = '+rad.k2.toFixed(5)+
        ' A, I(Kα2)/I(Kα1) = '+XRD_K2_RATIO.toFixed(2));
    out.push('# Structure: '+(MD.title||'crystal'));
    out.push('# Range: '+pat.min.toFixed(2)+' - '+(pat.min+pat.step*(pat.values.length-1)).toFixed(2)+' deg, step '+pat.step+' deg');
    out.push('# 2theta(deg) intensity');
    for(var i=0;i<pat.values.length;i++){
        out.push((pat.min+i*pat.step).toFixed(3)+' '+pat.values[i].toFixed(4));
    }
    var base=provBaseName(MD.filePath)||'structure.xyz';
    var stem=base;
    var dot=base.lastIndexOf('.');
    if(dot>0)stem=base.substring(0,dot);
    vscodeApi.postMessage({command:'saveFile',content:out.join('\\n')+'\\n',suggestedName:stem+'.xy',filePath:MD.filePath||''});
}
(function(){
    var openBtn=document.getElementById('xrd-open');
    var panel=document.getElementById('xrd-panel');
    if(!openBtn||!panel)return;
    openBtn.addEventListener('click',xrdOpenPanel);
    var closeBtn=document.getElementById('xrd-close');
    if(closeBtn)closeBtn.addEventListener('click',xrdMinimizePanel);
    var rad=document.getElementById('xrd-rad');
    if(rad)rad.addEventListener('change',xrdRefresh);
    var mx=document.getElementById('xrd-2max');
    if(mx)mx.addEventListener('change',xrdRefresh);
    var reBtn=document.getElementById('xrd-refresh');
    if(reBtn)reBtn.addEventListener('click',function(){
        if(xrdJob)xrdCancelCompute();else xrdRefresh();
    });
    var exBtn=document.getElementById('xrd-export');
    if(exBtn)exBtn.addEventListener('click',xrdExport);
    var importBtn=document.getElementById('xrd-import');
    var fileInput=document.getElementById('xrd-file');
    if(importBtn&&fileInput){
        importBtn.addEventListener('click',function(){fileInput.value='';fileInput.click()});
        fileInput.addEventListener('change',function(){
            var f=fileInput.files&&fileInput.files[0];
            if(!f)return;
            var reader=new FileReader();
            reader.onload=function(){xrdLoadExperimental(String(reader.result||''),f.name)};
            reader.onerror=function(){
                var el=document.getElementById('xrd-exp-info');
                if(el)el.textContent='Could not read '+f.name+'.';
            };
            reader.readAsText(f);
        });
    }
    var clearExp=document.getElementById('xrd-clear-exp');
    if(clearExp)clearExp.addEventListener('click',xrdClearExperimental);
    var cv=document.getElementById('xrd-canvas');
    if(cv)cv.addEventListener('click',xrdCanvasClick);
    var list=document.getElementById('xrd-list');
    if(list)list.addEventListener('click',function(ev){
        var row=ev.target&&ev.target.closest?ev.target.closest('.xr-row'):null;
        if(!row)return;
        var two=parseFloat(row.getAttribute('data-two'));
        var peaks=xrdResult&&xrdResult.ok?xrdResult.pattern.peaks:null;
        if(!peaks||!isFinite(two))return;
        var best=null,bd=Infinity,i;
        for(i=0;i<peaks.length;i++){
            var d=Math.abs(peaks[i].two-two);
            if(d<bd){bd=d;best=peaks[i]}
        }
        if(best)xrdSelectPeak(best);
    });
    xrdFillControls();
})();

// ---------------------------------------------------------------------------
// Voids and cavities (H9)
// ---------------------------------------------------------------------------
// A probe-accessible map of one unit cell. f(p) = min_i(|p - p_i| - (R_i +
// probe)) is the VMD Surf field: p is free where f > 0, i.e. where a probe
// sphere of the chosen radius sits without touching van der Waals matter. The
// map is sampled on a grid spanning the cell at a 0.25 A target spacing, at
// most 128 divisions per axis (128^3 is exactly the voxel budget), and the free
// voxels are split into 6-connected components that wrap through the faces.
//
// A component is "open" when its BFS lift meets an already-visited voxel with
// a different periodic offset: that is a path winding once around the crystal,
// so the void extends for ever and is a channel. Offsets that are merely
// nonzero do not matter: a cage centred on a cell corner is reached through the
// boundary but its lift is a single bounded blob, so it stays a closed cavity
// (drawn as two halves joined through the face, as the crystal repeats it).
//
// Volumes are per unit cell whatever supercell is on screen, and the report is
// read-only: it counts the space between the atoms, it does not touch them.
var CAV_TARGET=0.25;
var CAV_MAX_DIV=128;
// A probe-accessible pocket smaller than this is grid noise, not a cavity the
// reader would call one; the panel states the cut-off.
var CAV_MIN_VOLUME=0.5;
// Drawing a surface walks eight cubes per free voxel; past this many the panel
// reports the volume but declines to draw, rather than freezing the tab.
var CAV_MAX_SURFACE=150000;
var CAV_PROBE_DEFAULT=1.2;
// The deepest-point refinement costs one atom sweep per sample; a cell with
// hundreds of pockets would add up, so only this many of the largest regions
// get the full search. The rest keep their seed voxel's exact radius.
var CAV_DEEP_MAX=512;
var CAV_NB=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
var cavOpen=false,cavCache=null,cavSelected=-1,cavSurfaceMesh=null,cavNote='';

function cavWrap01(v){
    v=v-Math.floor(v);
    if(v>1-1e-6)v=0;
    return v;
}

function cavCellVolume(lv){
    return Math.abs(lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])
                   -lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])
                   +lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]));
}

function cavEdge(lv,k){
    return Math.sqrt(lv[k][0]*lv[k][0]+lv[k][1]*lv[k][1]+lv[k][2]*lv[k][2]);
}

// One division every CAV_TARGET Angstrom, rounded up, with 128 as the ceiling
// on each axis. The grid is a fraction-space lattice, so a skew cell is
// sampled as evenly along its edges as an orthogonal one with the same lengths.
function cavGridDims(lv){
    var d=[];
    for(var k=0;k<3;k++){
        var n=Math.ceil(cavEdge(lv,k)/CAV_TARGET);
        if(!(n>=2))n=2;
        if(n>CAV_MAX_DIV)n=CAV_MAX_DIV;
        d.push(n);
    }
    return d;
}

// The cell as the picture shows it: the ordered sites plus the selected
// disorder alternative, exactly the atoms the other readings work from.
function cavityAtoms(){
    var src=(CRY&&CRY.baseAtoms)||[],out=[];
    for(var i=0;i<src.length;i++){
        var a=src[i];
        if(!a||!neighBaseVisible(a))continue;
        var r=typeof VDW_R[a.element]==='number'?VDW_R[a.element]:DEFAULT_VDW_R;
        var f=neighFractional(a);
        out.push({r:r,
            fx:f[0]-Math.floor(f[0]),fy:f[1]-Math.floor(f[1]),fz:f[2]-Math.floor(f[2])});
    }
    return out;
}

// The clearance field, marked atom by atom: only the voxels inside an atom's
// bounding sphere can be made nearer to that atom, so the whole cell never
// meets the whole atom list. The fractional half-width of the Cartesian sphere
// is exact (it is the norm of an inverse-lattice row), and the atom's images in
// the neighbouring cells are walked too, so a site on a face blocks the space
// on the other side of it. The far field is a finite stand-in: an infinite
// value would turn the isosurface interpolation into NaN.
//
// The image walk is not clamped. In a skew cell the nearest image of a site
// can sit several cells away, and stopping at three would leave voxels free
// that the crystal actually blocks. A probe whose blocking sphere is wider
// than the cell's whole reach blocks every voxel through some image, and that
// case is answered without walking any of them.
var CAV_SLICE_MS=14;
function cavNow(){return Date.now()}
function cavFieldJob(atoms,lv,inv,probe,nx,ny,nz){
    var n=nx*ny*nz,i;
    var field=new Float32Array(n);
    for(i=0;i<n;i++)field[i]=1e6;
    var rx=Math.sqrt(inv[0][0]*inv[0][0]+inv[0][1]*inv[0][1]+inv[0][2]*inv[0][2]);
    var ry=Math.sqrt(inv[1][0]*inv[1][0]+inv[1][1]*inv[1][1]+inv[1][2]*inv[1][2]);
    var rz=Math.sqrt(inv[2][0]*inv[2][0]+inv[2][1]*inv[2][1]+inv[2][2]*inv[2][2]);
    var reach=0.5*(cavEdge(lv,0)+cavEdge(lv,1)+cavEdge(lv,2));
    for(i=0;i<atoms.length;i++){
        if(atoms[i].r+probe>=reach){
            for(i=0;i<n;i++)field[i]=-1e6;
            return {atoms:atoms,lv:lv,probe:probe,nx:nx,ny:ny,nz:nz,
                field:field,rx:rx,ry:ry,rz:rz,ai:atoms.length};
        }
    }
    return {atoms:atoms,lv:lv,probe:probe,nx:nx,ny:ny,nz:nz,
        field:field,rx:rx,ry:ry,rz:rz,ai:0};
}
// Marks a slice of the atom list; returns true when every atom is done. The
// time is checked between atoms, so a cell with thousands of sites still lets
// the page paint progress and answer a Cancel.
function cavFieldChunk(job,deadline){
    var atoms=job.atoms,lv=job.lv,probe=job.probe;
    var nx=job.nx,ny=job.ny,nz=job.nz,field=job.field;
    var rx=job.rx,ry=job.ry,rz=job.rz;
    while(job.ai<atoms.length){
        var at=atoms[job.ai++];
        var R=at.r+probe;
        if(R>0){
            var ex=R*rx,ey=R*ry,ez=R*rz;
            var mx=Math.ceil(ex-1e-9),my=Math.ceil(ey-1e-9),mz=Math.ceil(ez-1e-9);
            for(var sx=-mx;sx<=mx;sx++){
                var cx=at.fx+sx;
                var i0=Math.ceil((cx-ex)*nx-0.5),i1=Math.floor((cx+ex)*nx-0.5);
                if(i0<0)i0=0;if(i1>nx-1)i1=nx-1;
                if(i1<i0)continue;
                for(var sy=-my;sy<=my;sy++){
                    var cy=at.fy+sy;
                    var j0=Math.ceil((cy-ey)*ny-0.5),j1=Math.floor((cy+ey)*ny-0.5);
                    if(j0<0)j0=0;if(j1>ny-1)j1=ny-1;
                    if(j1<j0)continue;
                    for(var sz=-mz;sz<=mz;sz++){
                        var cz=at.fz+sz;
                        var k0=Math.ceil((cz-ez)*nz-0.5),k1=Math.floor((cz+ez)*nz-0.5);
                        if(k0<0)k0=0;if(k1>nz-1)k1=nz-1;
                        if(k1<k0)continue;
                        for(var ii=i0;ii<=i1;ii++){
                            var du=(ii+0.5)/nx-cx;
                            var ux=du*lv[0][0],uy=du*lv[0][1],uz=du*lv[0][2];
                            for(var jj=j0;jj<=j1;jj++){
                                var dv=(jj+0.5)/ny-cy;
                                var vx=ux+dv*lv[1][0],vy=uy+dv*lv[1][1],vz=uz+dv*lv[1][2];
                                for(var kk=k0;kk<=k1;kk++){
                                    var dw=(kk+0.5)/nz-cz;
                                    var dx=vx+dw*lv[2][0],dy=vy+dw*lv[2][1],dz=vz+dw*lv[2][2];
                                    var idx=(kk*ny+jj)*nx+ii;
                                    var d2=dx*dx+dy*dy+dz*dz;
                                    var lim=field[idx]+R;
                                    // sqrt(d2) - R beats the stored value only when
                                    // d2 < (stored + R)^2, so the root is taken on
                                    // the few voxels an atom actually improves.
                                    if(lim>0&&d2<lim*lim)field[idx]=Math.sqrt(d2)-R;
                                }
                            }
                        }
                    }
                }
            }
        }
        if(deadline&&cavNow()>=deadline)return false;
    }
    return true;
}
function cavField(atoms,lv,inv,probe,nx,ny,nz){
    var job=cavFieldJob(atoms,lv,inv,probe,nx,ny,nz);
    while(!cavFieldChunk(job,0)){}
    return job.field;
}

// Six-connected components of the free voxels, with a lift into the universal
// cover. The offsets are exact integers (a face crossing is +/-1 on one axis),
// so the winding test cannot drift. Per component: its voxel count and volume,
// whether any cycle winds, the unwrapped centroid (a real centre for a closed
// cavity, even one split by a cell face) and one interior voxel, which the
// deepest-point search uses when no coarse sample lands inside a small pocket.
function cavComponents(field,nx,ny,nz,voxelVolume){
    var n=nx*ny*nz,nxy=nx*ny;
    var comp=new Int32Array(n);
    var queue=new Int32Array(n);
    var oxA=new Int16Array(n),oyA=new Int16Array(n),ozA=new Int16Array(n);
    var comps=[],id=0;
    for(var s=0;s<n;s++){
        if(comp[s]||!(field[s]>0))continue;
        id++;
        var head=0,tail=0;
        queue[tail++]=s;comp[s]=id;
        var count=0,sumU=0,sumV=0,sumW=0;
        var wrapX=false,wrapY=false,wrapZ=false;
        while(head<tail){
            var cur=queue[head++];
            var x=cur%nx,y=((cur-x)/nx)%ny,z=(cur-x-y*nx)/nxy;
            count++;
            sumU+=x+oxA[cur]*nx;sumV+=y+oyA[cur]*ny;sumW+=z+ozA[cur]*nz;
            for(var d=0;d<6;d++){
                var nb=CAV_NB[d];
                var x2=x+nb[0],y2=y+nb[1],z2=z+nb[2];
                var ox=oxA[cur],oy=oyA[cur],oz=ozA[cur];
                if(x2<0){x2=nx-1;ox--}else if(x2>=nx){x2=0;ox++}
                if(y2<0){y2=ny-1;oy--}else if(y2>=ny){y2=0;oy++}
                if(z2<0){z2=nz-1;oz--}else if(z2>=nz){z2=0;oz++}
                var ni=(z2*ny+y2)*nx+x2;
                if(!(field[ni]>0))continue;
                if(!comp[ni]){
                    comp[ni]=id;queue[tail++]=ni;
                    oxA[ni]=ox;oyA[ni]=oy;ozA[ni]=oz;
                }else if(oxA[ni]!==ox||oyA[ni]!==oy||ozA[ni]!==oz){
                    if(oxA[ni]!==ox)wrapX=true;
                    if(oyA[ni]!==oy)wrapY=true;
                    if(ozA[ni]!==oz)wrapZ=true;
                }
            }
        }
        comps.push({
            id:id,count:count,volume:count*voxelVolume,
            closed:!wrapX&&!wrapY&&!wrapZ,wraps:[wrapX,wrapY,wrapZ],
            center:[cavWrap01(sumU/count/nx),cavWrap01(sumV/count/ny),cavWrap01(sumW/count/nz)],
            seed:s
        });
    }
    return {comp:comp,components:comps};
}

// The radius of the largest sphere that fits at one voxel: the distance to the
// nearest van der Waals surface, through the minimum image. Exact, unlike the
// marked field, which is only meaningful where an atom reached the voxel.
function cavFreeRadius(x,y,z,atoms,lv,nx,ny,nz){
    var u=(x+0.5)/nx,v=(y+0.5)/ny,w=(z+0.5)/nz;
    var best=1e9;
    for(var a=0;a<atoms.length;a++){
        var at=atoms[a];
        var du=u-at.fx,dv=v-at.fy,dw=w-at.fz;
        du-=Math.round(du);dv-=Math.round(dv);dw-=Math.round(dw);
        var dx=du*lv[0][0]+dv*lv[1][0]+dw*lv[2][0];
        var dy=du*lv[0][1]+dv*lv[1][1]+dw*lv[2][1];
        var dz=du*lv[0][2]+dv*lv[1][2]+dw*lv[2][2];
        var d=Math.sqrt(dx*dx+dy*dy+dz*dz)-at.r;
        if(d<best)best=d;
    }
    return best;
}

// The deepest point of every listed region and the radius of the sphere that
// fits there. The marked field cannot answer this - its free voxels all carry
// the far-field stand-in - so the exact free radius is sampled on a coarse
// stride shared by every region, then the best sample of each is refined by
// halving the step down to one voxel and finished with a hill climb, which
// makes the reported point a genuine local maximum rather than just the best
// of a fixed stencil.
//
// The walk is a resumable job: every free-radius sample costs one atom sweep,
// so on a cell with thousands of atoms the caller needs to stop between
// samples, paint its progress and answer a Cancel. The synchronous path runs
// the same states with no deadline, which keeps one implementation.
var CAV_N26=(function(){
    var out=[],oz,oy,ox;
    for(oz=-1;oz<=1;oz++)for(oy=-1;oy<=1;oy++)for(ox=-1;ox<=1;ox++){
        if(ox||oy||oz)out.push([ox,oy,oz]);
    }
    return out;
})();
function cavDeepenJob(comps,allCount,comp,atoms,lv,nx,ny,nz){
    var total=allCount+1,i;
    var bestV=new Float64Array(total);
    var bestX=new Int32Array(total),bestY=new Int32Array(total),bestZ=new Int32Array(total);
    for(i=0;i<total;i++)bestV[i]=-1e9;
    return {comps:comps,total:total,comp:comp,atoms:atoms,lv:lv,
        nx:nx,ny:ny,nz:nz,skip:!comps.length,
        bestV:bestV,bestX:bestX,bestY:bestY,bestZ:bestZ,
        cz:0,cy:0,cx:0,
        ci:0,stage:0,
        bx:0,by:0,bz:0,bv:0,
        hstep:4,nIdx:0,moved:false,guard:0,
        lbv:0,lbx:0,lby:0,lbz:0};
}
function cavDeepenCoarseChunk(job,deadline){
    if(job.skip)return true;
    var comp=job.comp,nx=job.nx,ny=job.ny,nz=job.nz;
    while(job.cz<nz){
        var id=comp[(job.cz*ny+job.cy)*nx+job.cx];
        if(id&&id<job.total){
            var v=cavFreeRadius(job.cx,job.cy,job.cz,job.atoms,job.lv,nx,ny,nz);
            if(v>job.bestV[id]){
                job.bestV[id]=v;
                job.bestX[id]=job.cx;job.bestY[id]=job.cy;job.bestZ[id]=job.cz;
            }
        }
        job.cx+=8;
        if(job.cx>=nx){
            job.cx=0;job.cy+=8;
            if(job.cy>=ny){job.cy=0;job.cz+=8}
        }
        if(deadline&&cavNow()>=deadline)return false;
    }
    return true;
}
function cavDeepenRefineChunk(job,deadline){
    if(job.skip)return true;
    var nx=job.nx,ny=job.ny,nz=job.nz;
    while(job.ci<job.comps.length){
        var c=job.comps[job.ci];
        if(job.stage===0){
            var bx=job.bestX[c.id],by=job.bestY[c.id],bz=job.bestZ[c.id],bv=job.bestV[c.id];
            if(bv<-1e8){
                // No coarse sample fell inside a pocket this small: start from
                // the region's own seed and let the refinement walk to its
                // middle.
                bx=c.seed%nx;by=((c.seed-bx)/nx)%ny;bz=(c.seed-bx-by*nx)/(nx*ny);
                bv=cavFreeRadius(bx,by,bz,job.atoms,job.lv,nx,ny,nz);
            }
            if(job.ci>=CAV_DEEP_MAX){
                c.deepest=[(bx+0.5)/nx,(by+0.5)/ny,(bz+0.5)/nz];
                c.radius=bv;
                job.ci++;
                continue;
            }
            job.bx=bx;job.by=by;job.bz=bz;job.bv=bv;
            job.hstep=4;job.nIdx=0;job.stage=1;
            continue;
        }
        if(job.stage===1){
            // Halving: 26 neighbours at four, then at two voxels, and only
            // the best of each pass is kept.
            if(job.nIdx===0){job.lbv=job.bv;job.lbx=job.bx;job.lby=job.by;job.lbz=job.bz}
            var nb=CAV_N26[job.nIdx++];
            var x2=((job.bx+nb[0]*job.hstep)%nx+nx)%nx;
            var y2=((job.by+nb[1]*job.hstep)%ny+ny)%ny;
            var z2=((job.bz+nb[2]*job.hstep)%nz+nz)%nz;
            var vv=cavFreeRadius(x2,y2,z2,job.atoms,job.lv,nx,ny,nz);
            if(vv>job.lbv){job.lbv=vv;job.lbx=x2;job.lby=y2;job.lbz=z2}
            if(job.nIdx>=CAV_N26.length){
                job.bx=job.lbx;job.by=job.lby;job.bz=job.lbz;job.bv=job.lbv;
                job.nIdx=0;
                if(job.hstep===4)job.hstep=2;
                else{job.stage=2;job.moved=true;job.guard=0}
            }
            if(deadline&&cavNow()>=deadline)return false;
            continue;
        }
        // Hill climb, one 26-neighbour sweep per pass, at most 64 passes.
        var nb2=CAV_N26[job.nIdx++];
        var x3=((job.bx+nb2[0])%nx+nx)%nx;
        var y3=((job.by+nb2[1])%ny+ny)%ny;
        var z3=((job.bz+nb2[2])%nz+nz)%nz;
        var v3=cavFreeRadius(x3,y3,z3,job.atoms,job.lv,nx,ny,nz);
        // Strictly better only: a tie would let the walk shuttle between
        // equal voxels.
        if(v3>job.bv+1e-12){job.bv=v3;job.bx=x3;job.by=y3;job.bz=z3;job.moved=true}
        if(job.nIdx>=CAV_N26.length){
            job.nIdx=0;
            if(job.moved&&job.guard<64){job.guard++;job.moved=false}
            else{
                c.deepest=[(job.bx+0.5)/nx,(job.by+0.5)/ny,(job.bz+0.5)/nz];
                c.radius=job.bv;
                job.ci++;job.stage=0;
            }
        }
        if(deadline&&cavNow()>=deadline)return false;
    }
    return true;
}
// The regions the report lists, with the refinement order: largest volume
// first, so the deepest-point budget lands on the rows a reader is looking at.
function cavGroup(split){
    var cavities=[],channels=[],free=0,i;
    for(i=0;i<split.components.length;i++){
        var c=split.components[i];
        if(c.volume<CAV_MIN_VOLUME)continue;
        free+=c.volume;
        (c.closed?cavities:channels).push(c);
    }
    function byVolume(a,b){return b.volume-a.volume}
    cavities.sort(byVolume);channels.sort(byVolume);
    var refine=cavities.concat(channels);
    refine.sort(byVolume);
    return {cavities:cavities,channels:channels,free:free,refine:refine};
}

// One analysis, from the atoms to the finished cache. The job runs in four
// states (field, splitting, coarse depth, refinement) so the panel can drive
// it a slice at a time; cavityCompute drives the same states without a
// deadline and is what the report's own tests measure.
function cavJobStart(probe){
    if(!CRY||!CRY_INV)return {error:{report:{ok:false,reason:'No unit cell.'}}};
    var lv=CRY.latticeVectors;
    var atoms=cavityAtoms();
    if(atoms.length===0)return {error:{report:{ok:false,reason:'No atoms in the unit cell.'}}};
    var dims=cavGridDims(lv);
    var nx=dims[0],ny=dims[1],nz=dims[2];
    var cellVolume=cavCellVolume(lv);
    var voxelVolume=cellVolume/(nx*ny*nz);
    return {probe:probe,atoms:atoms,lv:lv,nx:nx,ny:ny,nz:nz,
        cellVolume:cellVolume,voxelVolume:voxelVolume,
        stage:0,field:null,fieldJob:null,split:null,group:null,deepen:null};
}
function cavJobStep(job,deadline){
    if(job.stage===0){
        job.fieldJob=cavFieldJob(job.atoms,job.lv,CRY_INV,job.probe,job.nx,job.ny,job.nz);
        job.field=job.fieldJob.field;
        job.stage=1;
    }
    if(job.stage===1){
        if(!cavFieldChunk(job.fieldJob,deadline))return false;
        job.split=cavComponents(job.field,job.nx,job.ny,job.nz,job.voxelVolume);
        job.group=cavGroup(job.split);
        job.stage=2;
    }
    if(job.stage===2){
        if(!job.deepen)job.deepen=cavDeepenJob(job.group.refine,job.split.components.length,
            job.split.comp,job.atoms,job.lv,job.nx,job.ny,job.nz);
        if(!cavDeepenCoarseChunk(job.deepen,deadline))return false;
        job.stage=3;
    }
    if(job.stage===3){
        if(!cavDeepenRefineChunk(job.deepen,deadline))return false;
        job.stage=4;
    }
    return true;
}
function cavJobPercent(job){
    if(job.stage<=1){
        var fj=job.fieldJob;
        return fj&&fj.atoms.length?fj.ai/fj.atoms.length:1;
    }
    if(job.stage===2){
        var dj=job.deepen;
        return dj&&dj.nz>0?dj.cz/dj.nz:0;
    }
    if(job.stage===3&&job.deepen&&job.deepen.comps.length){
        return job.deepen.ci/job.deepen.comps.length;
    }
    return 1;
}
function cavJobProgressText(job){
    if(job.stage<=1)return 'Analysing: marking the atoms ('+Math.round(cavJobPercent(job)*100)+'%)';
    if(job.stage===2)return 'Analysing: sampling the deepest points ('+Math.round(cavJobPercent(job)*100)+'%)';
    if(job.stage===3)return 'Analysing: refining the deepest points ('+Math.round(cavJobPercent(job)*100)+'%)';
    return 'Analysing...';
}
function cavAssemble(job){
    var step=Math.max(cavEdge(job.lv,0)/job.nx,cavEdge(job.lv,1)/job.ny,cavEdge(job.lv,2)/job.nz);
    var g=job.group;
    var report={ok:true,reason:'',probe:job.probe,nx:job.nx,ny:job.ny,nz:job.nz,step:step,
        cellVolume:job.cellVolume,voxelVolume:job.voxelVolume,atoms:job.atoms.length,
        cavities:g.cavities,channels:g.channels,freeVolume:g.free,
        freeFraction:job.cellVolume>0?g.free/job.cellVolume:0};
    return {report:report,field:job.field,comp:job.split.comp,
        components:job.split.components,surfaceId:0,surfacePos:null};
}

function cavityCompute(probe){
    var job=cavJobStart(probe);
    if(job.error)return job.error;
    while(!cavJobStep(job,0)){}
    return cavAssemble(job);
}

// The isosurface f = 0 of one component, as a plain triangle list in absolute
// Cartesian coordinates. Only cubes whose corners include that component (and
// no other free component) are walked, so a point on a list chooses one void
// rather than every void in the cell. The six tetrahedra and the emitter are
// the ones the VMD Surf patches already use; a corner value is -f, which makes
// "free" the inside (negative) half exactly as those patches expect.
function cavitySurfacePositions(cache,id){
    var nx=cache.report.nx,ny=cache.report.ny,nz=cache.report.nz;
    var nxy=nx*ny;
    var comp=cache.comp,field=cache.field,lv=CRY.latticeVectors;
    var cubes=new Set();
    for(var s=0;s<comp.length;s++){
        if(comp[s]!==id)continue;
        var x=s%nx,y=((s-x)/nx)%ny,z=(s-x-y*nx)/nxy;
        for(var dx=-1;dx<=0;dx++){
            var bx=x+dx;if(bx<0)bx=nx-1;
            for(var dy=-1;dy<=0;dy++){
                var by=y+dy;if(by<0)by=ny-1;
                for(var dz=-1;dz<=0;dz++){
                    var bz=z+dz;if(bz<0)bz=nz-1;
                    cubes.add((bz*ny+by)*nx+bx);
                }
            }
        }
    }
    var pos=[];
    var v=new Array(8),px=new Array(8),py=new Array(8),pz=new Array(8);
    var cidx=new Int32Array(8);
    cubes.forEach(function(cube){
        var bx=cube%nx,by=((cube-bx)/nx)%ny,bz=(cube-bx-by*nx)/nxy;
        var hasTarget=false,hasOther=false;
        for(var k=0;k<8;k++){
            var x2=(bx+VMD_CORNER_X[k])%nx;
            var y2=(by+VMD_CORNER_Y[k])%ny;
            var z2=(bz+VMD_CORNER_Z[k])%nz;
            var ci=comp[(z2*ny+y2)*nx+x2];
            cidx[k]=ci;
            if(ci===id)hasTarget=true;
            else if(ci)hasOther=true;
            v[k]=-field[(z2*ny+y2)*nx+x2];
            // The positions are lifted out of the cell: a cube on a face gets
            // its far corner one cell over, or the tetrahedron would span the
            // whole cell backwards.
            var fx=(bx+VMD_CORNER_X[k]+0.5)/nx;
            var fy=(by+VMD_CORNER_Y[k]+0.5)/ny;
            var fz=(bz+VMD_CORNER_Z[k]+0.5)/nz;
            px[k]=fx*lv[0][0]+fy*lv[1][0]+fz*lv[2][0];
            py[k]=fx*lv[0][1]+fy*lv[1][1]+fz*lv[2][1];
            pz[k]=fx*lv[0][2]+fy*lv[1][2]+fz*lv[2][2];
        }
        if(!hasTarget||hasOther)return;
        for(var t=0;t<6;t++){
            var tet=VMD_TETRA[t];
            if(cidx[tet[0]]!==id&&cidx[tet[1]]!==id&&cidx[tet[2]]!==id&&cidx[tet[3]]!==id)continue;
            vmdEmitTet(pos,px,py,pz,v,tet[0],tet[1],tet[2],tet[3]);
        }
    });
    return pos;
}

function cavityClearSurface(){
    if(!cavSurfaceMesh)return;
    if(cavSurfaceMesh.parent)cavSurfaceMesh.parent.remove(cavSurfaceMesh);
    disposeMesh(cavSurfaceMesh);
    cavSurfaceMesh=null;
    needsRender=true;
}

// A live drag moves the atoms while the surface still describes the geometry
// before it: the mesh is taken down and the next rebuildScene puts it back.
function cavityHideSurface(){
    cavityClearSurface();
}

function cavityShowSelected(){
    cavityClearSurface();
    cavNote='';
    var note=document.getElementById('cav-note');
    if(note)note.textContent='';
    if(!cavCache||cavSelected<0||!moleculeGroup)return;
    var item=null,i;
    for(i=0;i<cavCache.components.length;i++){
        if(cavCache.components[i].id===cavSelected){item=cavCache.components[i];break}
    }
    if(!item)return;
    if(item.count>CAV_MAX_SURFACE){
        cavNote='Surface too large to draw ('+item.count+' free voxels); the volume is still exact.';
        if(note)note.textContent=cavNote;
        return;
    }
    var pos=cavCache.surfaceId===cavSelected?cavCache.surfacePos:cavitySurfacePositions(cavCache,cavSelected);
    if(!pos||pos.length===0){
        cavNote='No surface could be extracted for this void.';
        if(note)note.textContent=cavNote;
        return;
    }
    cavCache.surfaceId=cavSelected;cavCache.surfacePos=pos;
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
    geo.computeVertexNormals();
    var mat=new THREE.MeshPhongMaterial({color:0xb18cff,transparent:true,opacity:0.30,side:THREE.DoubleSide,depthWrite:false,shininess:24,specular:0x222222});
    cavSurfaceMesh=new THREE.Mesh(geo,mat);
    // Vertices are absolute, like the field and the mesh the surface was cut
    // from; the group is the one the atoms are centred in.
    cavSurfaceMesh.position.set(-CX,-CY,-CZ);
    moleculeGroup.add(cavSurfaceMesh);
    needsRender=true;
    cavNote='Surface drawn: '+(pos.length/9)+' triangles over '+item.count+' free voxels.';
    if(note)note.textContent=cavNote;
}

function cavityRefresh(){
    cavityClearSurface();
    if(!CRY||diffMode||!cavCache||cavSelected<0)return;
    cavityShowSelected();
}

function cavityInvalidate(){
    if(cavAbortJob())cavSetRunButton('Analyse','Scan the unit cell with this probe');
    cavCache=null;
    cavJobPrev=null;
    cavSelected=-1;
    cavNote='';
    cavityClearSurface();
    if(cavOpen)cavityRenderPanel();
}

function cavityRenderPanel(){
    var sum=document.getElementById('cav-summary');
    var list=document.getElementById('cav-list');
    if(!sum||!list)return;
    if(cavJob){
        sum.textContent=cavJobProgressText(cavJob);
        list.innerHTML='';
        return;
    }
    if(!cavCache){
        sum.textContent=CRY?'Not analysed yet.':'No unit cell.';
        list.innerHTML='';
        return;
    }
    var rep=cavCache.report;
    if(!rep||!rep.ok){
        sum.textContent=(rep&&rep.reason)?rep.reason:'No unit cell.';
        list.innerHTML='';
        return;
    }
    var nc=rep.cavities.length,nch=rep.channels.length;
    sum.innerHTML='<b>'+nc+'</b> closed cave'+(nc===1?'':'s')+' &middot; <b>'+nch+'</b> open channel'+(nch===1?'':'s')+
        '<br><span class="cv-sub">grid '+rep.nx+'&times;'+rep.ny+'&times;'+rep.nz+' ('+rep.step.toFixed(2)+' &Aring;) &middot; probe '+rep.probe.toFixed(2)+' &Aring; &middot; '+
        rep.atoms+' atoms &middot; free '+rep.freeVolume.toFixed(1)+' &Aring;&sup3; ('+(rep.freeFraction*100).toFixed(1)+'% of the cell)</span>';
    var html='<div class="cv-head"><span class="cv-name">void</span><span class="cv-vol">vol/&Aring;&sup3;</span><span class="cv-c">centre (fractional)</span><span class="cv-r">r_free/&Aring;</span></div>';
    function addRows(arr,prefix){
        for(var i=0;i<arr.length;i++){
            var c=arr[i],ctr=c.closed?c.center:c.deepest;
            html+='<div class="cv-row'+(cavSelected===c.id?' sel':'')+'" data-i="'+c.id+'" title="Click to draw or hide this surface">'+
                '<span class="cv-name">'+prefix+' '+(i+1)+'</span>'+
                '<span class="cv-vol">'+c.volume.toFixed(1)+'</span>'+
                '<span class="cv-c">'+ctr[0].toFixed(3)+' '+ctr[1].toFixed(3)+' '+ctr[2].toFixed(3)+'</span>'+
                '<span class="cv-r">'+c.radius.toFixed(2)+'</span></div>';
        }
    }
    addRows(rep.cavities,'Cavity');
    addRows(rep.channels,'Channel');
    if(nc+nch===0)html+='<div class="cv-empty">No probe-accessible void at this radius.</div>';
    list.innerHTML=html;
}

function cavitySelect(id){
    cavSelected=cavSelected===id?-1:id;
    cavityRenderPanel();
    cavityShowSelected();
}

// The analysis is driven a slice at a time: each pass works for CAV_SLICE_MS
// and hands the event loop back, so even a cell with thousands of atoms keeps
// the page alive and the summary line counts up. The Analyse button becomes
// Cancel while the job runs, and a cell edit drops the job with the report.
var cavJob=null,cavJobSeq=0,cavJobPrev=null;

function cavSetRunButton(label,title){
    var run=document.getElementById('cav-run');
    if(run){run.textContent=label;run.title=title}
}

function cavAbortJob(){
    if(!cavJob)return false;
    cavJob=null;
    cavJobSeq++;
    return true;
}

function cavityAnalyse(){
    if(!CRY)return;
    var input=document.getElementById('cav-probe');
    var probe=parseFloat(input?input.value:CAV_PROBE_DEFAULT);
    if(!(probe>=0)||!isFinite(probe))probe=CAV_PROBE_DEFAULT;
    if(probe>20)probe=20;
    if(input)input.value=String(probe);
    // A probe change while a run is going starts the new run from scratch.
    if(cavAbortJob())cavSetRunButton('Analyse','Scan the unit cell with this probe');
    // The old report stays reachable for a cancel, but it is off the panel
    // while the new run works: its rows would draw surfaces the new numbers
    // no longer describe.
    cavJobPrev=cavCache;
    cavCache=null;
    cavSelected=-1;
    cavNote='';
    cavityClearSurface();
    var noteEl=document.getElementById('cav-note');
    if(noteEl)noteEl.textContent='';
    var job=cavJobStart(probe);
    if(job.error){
        cavJobPrev=null;
        cavCache=job.error;
        cavityRenderPanel();
        if(cavOpen){syncPanelChips();layoutPanels();}
        return;
    }
    cavJob=job;
    var seq=++cavJobSeq;
    cavSetRunButton('Cancel','Stop the cavity analysis');
    cavityRenderPanel();
    setTimeout(function(){cavJobTick(seq)},0);
}

function cavJobTick(seq){
    if(!cavJob||seq!==cavJobSeq)return;
    var job=cavJob;
    var done;
    try{
        done=cavJobStep(job,cavNow()+CAV_SLICE_MS);
    }catch(e){
        cavJob=null;
        cavJobSeq++;
        cavJobPrev=null;
        cavSetRunButton('Analyse','Scan the unit cell with this probe');
        cavCache={report:{ok:false,reason:'The analysis stopped: '+
            (e&&e.message?e.message:'unexpected error')}};
        cavityRenderPanel();
        return;
    }
    if(!cavJob||seq!==cavJobSeq)return;
    if(!done){
        var sum=document.getElementById('cav-summary');
        if(sum)sum.textContent=cavJobProgressText(job);
        setTimeout(function(){cavJobTick(seq)},0);
        return;
    }
    cavJob=null;
    cavJobPrev=null;
    cavCache=cavAssemble(job);
    cavSetRunButton('Analyse','Scan the unit cell with this probe');
    cavityRenderPanel();
    var note=document.getElementById('cav-note');
    if(note)note.textContent='';
    if(cavOpen){syncPanelChips();layoutPanels();}
}

function cavityCancelAnalysis(){
    if(!cavAbortJob())return;
    cavSetRunButton('Analyse','Scan the unit cell with this probe');
    if(cavJobPrev){cavCache=cavJobPrev;cavJobPrev=null}
    cavityRenderPanel();
    var note=document.getElementById('cav-note');
    if(note)note.textContent=cavCache?'Analysis cancelled; the report above is the previous run.':'Analysis cancelled.';
    if(cavOpen)layoutPanels();
}

function cavityOpenPanel(){
    if(!CRY)return;
    var panel=document.getElementById('cav-panel');
    if(!panel)return;
    // The scan reads one cell; in a comparison it would cover the right-hand
    // structure and report voids for the left one only.
    if(diffMode){
        modeInfoEl.textContent='Cavity analysis is not available in diff mode';
        return;
    }
    cavOpen=true;
    panelMinState.cav=false;
    panel.style.display='block';
    if(!cavCache&&!cavJob)cavityAnalyse();
    else cavityRenderPanel();
    syncPanelChips();
    layoutPanels();
}

// The x folds the report into its chip; the chip and the crystal panel's
// button both open it again with the probe and the drawn void it had.
function cavityMinimizePanel(){
    if(!cavOpen)return;
    panelMinState.cav=true;
    syncPanelChips();
    layoutPanels();
}

function cavityClosePanel(){
    // A scan still running has no panel to report to; stop it with the report.
    if(cavAbortJob())cavSetRunButton('Analyse','Scan the unit cell with this probe');
    var panel=document.getElementById('cav-panel');
    if(panel)panel.style.display='none';
    cavOpen=false;
    panelMinState.cav=false;
    cavSelected=-1;
    cavityClearSurface();
    syncPanelChips();
    layoutPanels();
}

(function(){
    var openBtn=document.getElementById('cav-open');
    if(!openBtn)return;
    openBtn.addEventListener('click',cavityOpenPanel);
    var closeBtn=document.getElementById('cav-close');
    if(closeBtn)closeBtn.addEventListener('click',cavityMinimizePanel);
    var run=document.getElementById('cav-run');
    if(run)run.addEventListener('click',function(){
        if(cavJob)cavityCancelAnalysis();
        else cavityAnalyse();
    });
    var probe=document.getElementById('cav-probe');
    if(probe)probe.addEventListener('change',cavityAnalyse);
    var clear=document.getElementById('cav-clear');
    if(clear)clear.addEventListener('click',function(){
        cavSelected=-1;
        cavityClearSurface();
        cavityRenderPanel();
        var n=document.getElementById('cav-note');
        if(n)n.textContent='';
    });
    var list=document.getElementById('cav-list');
    if(list)list.addEventListener('click',function(ev){
        var row=ev.target&&ev.target.closest?ev.target.closest('.cv-row'):null;
        if(!row)return;
        var id=parseInt(row.getAttribute('data-i'),10);
        if(!isNaN(id))cavitySelect(id);
    });
})();

function fmtMeasurement(kind,sel,atomsArr){
    if(!sel||sel.length===0)return null;
    if(kind==='bondLength'&&sel.length===2){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]];
        var d=dist(a1,a2);
        return {label:atomLabel(a1,sel[0])+' - '+atomLabel(a2,sel[1]),value:d.toFixed(4)+' A',num:d};
    }
    if(kind==='bondAngle'&&sel.length===3){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]];
        var ang=angle(a1,a2,a3);
        return {label:atomLabel(a1,sel[0])+' - '+atomLabel(a2,sel[1])+' - '+atomLabel(a3,sel[2]),value:ang.toFixed(2)+' deg',num:ang};
    }
    if(kind==='dihedral'&&sel.length===4){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]],a4=atomsArr[sel[3]];
        var dih=dihedral(a1,a2,a3,a4);
        return {label:atomLabel(a1,sel[0])+' - '+atomLabel(a2,sel[1])+' - '+atomLabel(a3,sel[2])+' - '+atomLabel(a4,sel[3]),value:dih.toFixed(2)+' deg',num:dih};
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

function showModal(html,cb,wide){if(vibActive||vibPaused)stopVibration();stopMoveSession();modalEl.className=wide?'modal-wide':'';modalEl.innerHTML=html;modalOverlay.classList.add('show');modalCallback=cb}
function hideModal(){modalOverlay.classList.remove('show');modalCallback=null;modalEl.className=''}

function showBondLengthModal(){
    // A bond drawn as a stub across a cell face has its two atoms a whole cell
    // apart on purpose: the distance on screen is not the bond's length, and
    // writing that distance "back" would move the cell by a displacement the
    // bond does not have. Both ends have to belong to one periodic image.
    if(crystalCrossFaceBondBetween(selectedAtoms[0],selectedAtoms[1])){
        modeInfoEl.textContent='Bond Length: this bond crosses a cell face - both atoms must be in one periodic image (grow the supercell and pick them there)';
        return;
    }
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]];
    var cur=dist(a1,a2);
    // In a crystal the two atoms picked can be two periodic images rather than
    // the two ends of one bond in the cell. The distance on screen is what the
    // dialog edits, but the structure's own distance between those atoms is
    // the shortest periodic image of the displacement, and when the two differ
    // the number being edited is not a bond length at all. The row is read-only:
    // editing semantics (which reject a split bond outright) are unchanged.
    var shortestRow='';
    if(CRY&&CRY_INV){
        var lvS=CRY.latticeVectors,invS=CRY_INV;
        var fdx=invS[0][0]*(a2.x-a1.x)+invS[0][1]*(a2.y-a1.y)+invS[0][2]*(a2.z-a1.z);
        var fdy=invS[1][0]*(a2.x-a1.x)+invS[1][1]*(a2.y-a1.y)+invS[1][2]*(a2.z-a1.z);
        var fdz=invS[2][0]*(a2.x-a1.x)+invS[2][1]*(a2.y-a1.y)+invS[2][2]*(a2.z-a1.z);
        var nearS=crystalShortestFractional(fdx,fdy,fdz,cur,lvS,crystalPlaneWidths(lvS));
        var offImage=Math.abs(nearS.dist-cur)>1e-3;
        shortestRow='<div class="current-val">Shortest image: '+nearS.dist.toFixed(4)+' A'+
            (offImage?' (these two images are not the nearest pair of their atoms)':'')+'</div>';
    }
    var existingBond=MD.bonds.find(function(b){return(b.atom1===selectedAtoms[0]&&b.atom2===selectedAtoms[1])||(b.atom1===selectedAtoms[1]&&b.atom2===selectedAtoms[0])});
    var curOrder=existingBond?existingBond.order:1;
    var hadBondAtOpen=!!existingBond;
    saveOriginal();
    var n1=atomLabel(a1,selectedAtoms[0]),n2=atomLabel(a2,selectedAtoms[1]);
    showModal('<h3>Adjust Bond Length</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(4)+' A, Bond order: '+curOrder+'</div>'+
        shortestRow+
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
            if(existingBond){
                // A crystal's bonds on screen are drawn from its base bonds, so a
                // removal has to reach the cell bond the picture was drawn from;
                // otherwise the next rebuildCrystal() would draw it again.
                if(CRY){
                    var gone=crystalFindBaseBond(crystalBaseBondFromDisplay(existingBond));
                    if(gone)CRY.baseBonds.splice(CRY.baseBonds.indexOf(gone),1);
                }
                MD.bonds=MD.bonds.filter(function(b){return b!==existingBond});existingBond=null;rebuildScene()
            }
        }else{
            if(!existingBond){
                // A ghost has no valence, so this pair cannot be given a bond
                // even from here: the dialog is the one place a bond can be
                // created by hand, and a ghost would otherwise leave with one
                // (which the drawing refuses and a saved connect section would
                // carry into a calculation). The control goes back to None so it
                // never claims a bond that does not exist.
                if(isGhostAtom(a1)||isGhostAtom(a2)){
                    this.value='0';
                    modeInfoEl.textContent='A ghost atom cannot be bonded - no bond was created';
                    return;
                }
                existingBond={atom1:selectedAtoms[0],atom2:selectedAtoms[1],order:newOrder};
                // A bond a crystal did not have lives in the base cell too: its
                // shift is the number of cells between the two images picked.
                if(CRY){
                    var made=crystalBaseBondFromDisplay(existingBond);
                    if(made&&!crystalFindBaseBond(made))CRY.baseBonds.push({atom1:made.atom1,atom2:made.atom2,order:newOrder,shift:made.shift});
                }
                MD.bonds.push(existingBond);
            }else{
                existingBond.order=newOrder;
                if(CRY){
                    var baseBond=crystalFindBaseBond(crystalBaseBondFromDisplay(existingBond));
                    if(baseBond)baseBond.order=newOrder;
                }
            }
            rebuildScene();
        }
    });
    document.getElementById('m-ok').addEventListener('click',function(){
        // The dialog edits the picture from the moment it opens; a crystal saves
        // its base cell, so the previewed geometry is carried there now, before
        // the session baseline this dialog captured is dropped.
        if(CRY)commitCrystalCoords();
        provMarkIfMoved('bondLen');hideModal();originalCoords=null;resetSelection()
    });
    document.getElementById('m-cancel').addEventListener('click',function(){
        if(existingBond&&originalCoords){existingBond.order=curOrder}
        // A bond this dialog created is not part of the structure it opened on,
        // so cancelling takes it out of the picture as well as out of the cell.
        if(existingBond&&!hadBondAtOpen)MD.bonds=MD.bonds.filter(function(b){return b!==existingBond});
        // A bond-order change went into the base cell as it was made; the
        // snapshot this dialog opened with is the state to put back before the
        // snapshot is dropped.
        if(CRY)restoreCrystalBaseFromUndoTop();
        restoreOriginal();rebuildScene();undoStack.pop();updateUndoBtn();hideModal();originalCoords=null;resetSelection()
    });
}

function showBondAngleModal(){
    // Same rule as the bond length: an angle read off a bond drawn across a
    // face is an angle between two different images, not the angle in the cell.
    if(crystalCrossFaceBondBetween(selectedAtoms[0],selectedAtoms[1])||crystalCrossFaceBondBetween(selectedAtoms[1],selectedAtoms[2])){
        modeInfoEl.textContent='Bond Angle: a bond here crosses a cell face - both atoms of each bond must be in one periodic image';
        return;
    }
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]];
    var cur=angle(a1,a2,a3);
    saveOriginal();
    var n1=atomLabel(a1,selectedAtoms[0]),n2=atomLabel(a2,selectedAtoms[1]),n3=atomLabel(a3,selectedAtoms[2]);
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
    document.getElementById('m-ok').addEventListener('click',function(){if(CRY)commitCrystalCoords();provMarkIfMoved('bondAngle');hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){restoreOriginal();rebuildScene();undoStack.pop();updateUndoBtn();hideModal();originalCoords=null;resetSelection()});
}

function showDihedralModal(){
    // And the dihedral: any of its three bonds crossing a face makes the
    // measured angle one between images rather than the cell's own.
    if(crystalCrossFaceBondBetween(selectedAtoms[0],selectedAtoms[1])||crystalCrossFaceBondBetween(selectedAtoms[1],selectedAtoms[2])||crystalCrossFaceBondBetween(selectedAtoms[2],selectedAtoms[3])){
        modeInfoEl.textContent='Dihedral: a bond here crosses a cell face - both atoms of each bond must be in one periodic image';
        return;
    }
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]],a4=MD.atoms[selectedAtoms[3]];
    var cur=dihedral(a1,a2,a3,a4);
    saveOriginal();
    var n1=atomLabel(a1,selectedAtoms[0]),n2=atomLabel(a2,selectedAtoms[1]),n3=atomLabel(a3,selectedAtoms[2]),n4=atomLabel(a4,selectedAtoms[3]);
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
    document.getElementById('m-ok').addEventListener('click',function(){if(CRY)commitCrystalCoords();provMarkIfMoved('dihedral');hideModal();originalCoords=null;resetSelection()});
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
        '<option value="X">X (dummy point)</option>'+
        '<option selected>H</option><option>He</option>'+
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
        '<label style="display:flex;align-items:center;gap:6px;cursor:pointer">'+
        '<input type="checkbox" id="m-ghost" style="flex:0 0 auto;width:auto">'+
        '<span>Ghost atom - no nucleus, no electrons, no bond</span></label>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    // A ghost carries no bond, and a dummy point has no element to bond with:
    // the two controls that describe the new bond are turned off rather than
    // silently ignored, so the dialog says what it will do.
    var ghostSync=function(){
        var g=document.getElementById('m-ghost');
        var el=document.getElementById('m-elem');
        if(!g||!el)return;
        var isG=g.checked||el.value==='X';
        var bo=document.getElementById('m-bond-order');
        if(bo){bo.disabled=isG;bo.title=isG?'A ghost atom is never bonded':''}
    };
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
        ghostSync();
    });
    document.getElementById('m-ghost').addEventListener('change',ghostSync);
    document.getElementById('m-ok').addEventListener('click',function(){
        var el=document.getElementById('m-elem').value;
        var bl=parseFloat(document.getElementById('m-val').value)||1.5;
        var bondOrder=parseFloat(document.getElementById('m-bond-order').value)||1;
        var ghost=document.getElementById('m-ghost').checked||el==='X';
        // Ghosts are refused in a periodic cell: every crystal format this
        // viewer writes (CIF, POSCAR, VESTA, cube) states a species per site
        // and has no way to say "this one has no nucleus", so a ghost would
        // come back as a real atom on the next load.
        if(ghost&&CRY){
            hideModal();resetSelection();
            modeInfoEl.textContent='Ghost atoms are not available for crystal structures';
            return;
        }
        pushUndo();
        provMark(ghost?'addGhostAtoms':'addAtoms',1);
        var dir={x:0,y:0,z:1};
        if(CRY){
            var baseAnchorIdx=anchor.baseIdx!=null?anchor.baseIdx:anchorIdx;
            var baseAnchor=CRY.baseAtoms[baseAnchorIdx];
            // Neighbour directions are read through the image each bond names,
            // not through the raw Cartesian difference: a molecule wrapped
            // across a cell face stores those neighbours a whole cell away, and
            // the raw difference would point the new atom through the wall.
            var nbrs=[];
            CRY.baseBonds.forEach(function(b){
                var s=b.shift||ZERO_SHIFT;
                var other=-1,sign=0;
                if(b.atom1===baseAnchorIdx){other=b.atom2;sign=1}
                else if(b.atom2===baseAnchorIdx){other=b.atom1;sign=-1}
                if(other<0||!CRY.baseAtoms[other])return;
                var v=crystalImageDirection(baseAnchor,CRY.baseAtoms[other],sign,s);
                if(v)nbrs.push(v);
            });
            if(nbrs.length>0)dir=pickFreeDirection(nbrs);
            // The new site is wrapped into the cell — that is where every other
            // base atom lives and where rebuildCrystal() would put it anyway, so
            // the display and the file agree — and the bond records which image
            // of it the anchor is bonded to. Without that shift a site added
            // across a face would be drawn as one long bond straight through the
            // cell instead of as a bond into the next cell.
            var newX=baseAnchor.x+dir.x*bl,newY=baseAnchor.y+dir.y*bl,newZ=baseAnchor.z+dir.z*bl;
            var newShift=[0,0,0];
            if(CRY_INV){
                var invC=CRY_INV,lvC=CRY.latticeVectors;
                var afx=invC[0][0]*baseAnchor.x+invC[0][1]*baseAnchor.y+invC[0][2]*baseAnchor.z;
                var afy=invC[1][0]*baseAnchor.x+invC[1][1]*baseAnchor.y+invC[1][2]*baseAnchor.z;
                var afz=invC[2][0]*baseAnchor.x+invC[2][1]*baseAnchor.y+invC[2][2]*baseAnchor.z;
                var nfx=invC[0][0]*newX+invC[0][1]*newY+invC[0][2]*newZ;
                var nfy=invC[1][0]*newX+invC[1][1]*newY+invC[1][2]*newZ;
                var nfz=invC[2][0]*newX+invC[2][1]*newY+invC[2][2]*newZ;
                // The image of the new site that the bond points at: the one
                // whose position is the spot that was just computed, measured
                // from the cell the anchor's own coordinates sit in.
                newShift=[Math.floor(nfx)-Math.floor(afx),Math.floor(nfy)-Math.floor(afy),Math.floor(nfz)-Math.floor(afz)];
                var wfx=nfx-Math.floor(nfx),wfy=nfy-Math.floor(nfy),wfz=nfz-Math.floor(nfz);
                newX=wfx*lvC[0][0]+wfy*lvC[1][0]+wfz*lvC[2][0];
                newY=wfx*lvC[0][1]+wfy*lvC[1][1]+wfz*lvC[2][1];
                newZ=wfx*lvC[0][2]+wfy*lvC[1][2]+wfz*lvC[2][2];
            }
            var newBaseIdx=CRY.baseAtoms.length;
            CRY.baseAtoms.push({element:el,x:newX,y:newY,z:newZ,index:newBaseIdx,baseIdx:newBaseIdx,occupancy:1});
            CRY.baseBonds.push({atom1:baseAnchorIdx,atom2:newBaseIdx,order:bondOrder,shift:newShift});
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
            var newAtom={element:el,x:anchor.x+dir.x*bl,y:anchor.y+dir.y*bl,z:anchor.z+dir.z*bl,color:MD.atomColors[el]||'#FF1493',index:newIdx};
            if(ghost)newAtom.ghost=true;
            MD.atoms.push(newAtom);
            // The bond is the one thing a ghost does not get: a ghost has no
            // valence, and bonding one to its anchor would make the bond
            // perception, the force field and every saved format disagree about
            // what it is. The anchor being a ghost is refused for the same
            // reason - a point has nothing to bond with.
            if(!ghost&&!isGhostAtom(anchor))MD.bonds.push({atom1:anchorIdx,atom2:newIdx,order:bondOrder});
            rebuildScene();hideModal();originalCoords=null;resetSelection();
            if(ghost){
                modeInfoEl.textContent='Added a ghost atom of '+(el==='X'?'no element (dummy point)':el)+' - not bonded';
            }else if(isGhostAtom(anchor)){
                modeInfoEl.textContent='Added '+el+' - not bonded: the anchor is a ghost atom';
            }
        }
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
}

function showReplaceAtomModal(){
    if(selectedAtoms.length===0)return;
    var names=formatAtomList(selectedAtoms);
    // Pre-ticking the box when everything selected is already a ghost makes the
    // dialog show the state it will keep, so "Replace" on a ghost selection is
    // not a silent un-ghosting.
    var allGhost=selectedAtoms.every(function(i){return isGhostAtom(MD.atoms[i])});
    showModal('<h3>Replace Atom</h3>'+
        '<div class="current-val">Replace '+selectedAtoms.length+' atom(s): '+names+'</div>'+
        '<label>New element:</label><select id="m-elem">'+
        '<option value="X">X (dummy point)</option>'+
        '<option selected>H</option><option>He</option>'+
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
        '<label style="display:flex;align-items:center;gap:6px;cursor:pointer">'+
        '<input type="checkbox" id="m-ghost"'+(allGhost?' checked':'')+' style="flex:0 0 auto;width:auto">'+
        '<span>Ghost - keep each atom\\'s element, drop its bonds</span></label>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Replace</button></div>',null);
    var ghostSync=function(){
        var g=document.getElementById('m-ghost');
        var el=document.getElementById('m-elem');
        if(!g||!el)return;
        // Ticked: the element list has no say - every atom stays the element it
        // already is and merely loses its nucleus. That is the gesture a
        // counterpoise setup needs ("make this whole fragment ghosts"), and it
        // is the only reading under which a mixed selection makes sense.
        el.disabled=g.checked;
        el.title=g.checked?'A ghost keeps each atom\\'s own element':'';
    };
    document.getElementById('m-ghost').addEventListener('change',ghostSync);
    ghostSync();
    document.getElementById('m-ok').addEventListener('click',function(){
        var ghost=document.getElementById('m-ghost').checked;
        var el=document.getElementById('m-elem').value;
        if(ghost&&CRY){
            hideModal();resetSelection();
            modeInfoEl.textContent='Ghost atoms are not available for crystal structures';
            return;
        }
        // Counted before the change so the title records what actually happened:
        // a Replace that re-declares atoms already in the state asked for is not
        // an edit and does not get a clause. Choosing X keeps an atom a ghost
        // whatever the switch says - the dummy element *is* one - so a real
        // conversion needs both the switch off and a real element.
        var n=selectedAtoms.length,changedGhost=0,changedReal=0;
        selectedAtoms.forEach(function(idx){
            var a=MD.atoms[idx];
            if(!a)return;
            var was=isGhostAtom(a);
            if(ghost&&!was)changedGhost++;
            else if(!ghost&&was&&el!=='X')changedReal++;
        });
        pushUndo();
        if(ghost){
            if(changedGhost>0)provMark('ghostAtoms',changedGhost);
        }else{
            provMark('replaceAtoms',n);
            if(changedReal>0)provMark('unghostAtoms',changedReal);
        }
        var becameReal=[];
        selectedAtoms.forEach(function(idx){
            var a=MD.atoms[idx];
            if(!a)return;
            if(ghost){
                a.ghost=true;
            }else{
                if(isGhostAtom(a)&&el!=='X')becameReal.push(idx);
                a.element=el;
                a.ghost=false;
                a.color=MD.atomColors[el]||'#FF1493';
            }
            if(CRY){
                var bi=a.baseIdx!=null?a.baseIdx:idx;
                if(CRY.baseAtoms[bi]){
                    CRY.baseAtoms[bi].element=a.element;
                    CRY.baseAtoms[bi].color=a.color;
                    CRY.baseAtoms[bi].ghost=a.ghost;
                }
            }
        });
        // A ghost has no bonds, so becoming one takes the atom's bonds with it;
        // becoming real again is the reverse - the atom is re-derived from the
        // geometry the same way a load would, and only the bonds that touch a
        // converted atom are added, so no existing bond changes.
        dropGhostBonds();
        if(becameReal.length>0)rebondAtoms(becameReal);
        if(CRY){rebuildCrystal();}
        rebuildScene();hideModal();originalCoords=null;resetSelection();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
}

/** Remove every bond that touches a ghost atom of the structure on screen. */
function dropGhostBonds(){
    if(MD.bonds.length===0)return;
    var kept=MD.bonds.filter(function(b){
        var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
        if(!a1||!a2)return false;
        return !isGhostAtom(a1)&&!isGhostAtom(a2);
    });
    if(kept.length!==MD.bonds.length)MD.bonds=kept;
}

/**
 * Give the named atoms the bonds a fresh perception would derive for them.
 *
 * Used when a ghost is turned back into a real atom: the ghost had no bonds, so
 * the atom is left isolated by the change, and no other pass would reconnect it
 * (the viewer only re-derives bonds for a moved region). The orders come from
 * the same detector a load uses, and only bonds that touch one of the converted
 * atoms and do not exist yet are added — a bond the structure already had is
 * never rewritten by this.
 */
function rebondAtoms(indices){
    if(!indices||indices.length===0)return;
    var want={};
    indices.forEach(function(i){want[i]=1});
    var have={};
    MD.bonds.forEach(function(b){
        if(b.atom1==null||b.atom2==null)return;
        have[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=1;
    });
    var derived;
    try{derived=detectBondsFromAtoms(MD.atoms,MD.charge)}catch(e){return}
    derived.forEach(function(b){
        if(!want[b.atom1]&&!want[b.atom2])return;
        var key=Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2);
        if(have[key])return;
        have[key]=1;
        MD.bonds.push({atom1:Math.min(b.atom1,b.atom2),atom2:Math.max(b.atom1,b.atom2),order:b.order});
    });
}

function showDeleteAtomModal(){
    var idx=selectedAtoms[0];
    var a=MD.atoms[idx];
    if(!a)return;
    var name=atomLabel(a,idx);
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
                index:newIdx,occupancy:ba.occupancy!=null?ba.occupancy:1,baseIdx:newIdx,color:color,
                siteLabel:ba.siteLabel,asymIdx:ba.asymIdx,symOp:ba.symOp,
                constraint:ba.constraint,disorderGroup:ba.disorderGroup,uAniso:ba.uAniso});
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
    for(i=0;i<provEdits.length;i++)edits.push({k:provEdits[i].k,n:provEdits[i].n,d:provEdits[i].d||null});
    for(k in provModes){if(provModes.hasOwnProperty(k))modes[k]={freq:provModes[k].freq,scale:provModes[k].scale}}
    return {edits:edits,merged:provMerged.slice(),modes:modes};
}
function provRestore(p){
    var i,k;
    provEdits=[];provMerged=[];provModes={};
    if(!p)return;
    if(p.edits)for(i=0;i<p.edits.length;i++)provEdits.push({k:p.edits[i].k,n:p.edits[i].n,d:p.edits[i].d||null});
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
    provEdits.push({k:kind,n:n||1,d:null});
}
/** A stateful choice rather than a counted edit: the cell in force, the
 *  disorder alternative on screen, whether its sites count as fully occupied.
 *  Only the latest choice describes the structure, so it replaces the clause
 *  instead of adding a second one ("2 x cell edited" would be wrong). */
function provSet(kind,detail){
    for(var i=0;i<provEdits.length;i++){
        if(provEdits[i].k===kind){
            provEdits[i].n=1;
            provEdits[i].d=detail||null;
            return;
        }
    }
    provEdits.push({k:kind,n:1,d:detail||null});
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
        case 'addGhostAtoms':return n===1?'1 ghost atom added':n+' ghost atoms added';
        case 'ghostAtoms':return n===1?'1 atom made a ghost':n+' atoms made ghosts';
        case 'unghostAtoms':return n===1?'1 ghost made real':n+' ghosts made real';
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
        case 'disorderGroup':return e.d?('disorder group '+e.d+' shown'):'disorder alternative shown';
        case 'occNormalize':return 'disorder sites counted at full occupancy';
        case 'cellEdit':return e.d?('cell edited to '+e.d):'cell edited';
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
function provTitle(fallback,keepOriginal,extraHead){
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
    // A caller that has something the ordinary clauses cannot say — the cell and
    // bounds a cluster was cut from — puts it here, ahead of the edits, so the
    // 190-character cap can only cost the edit list and never the structure's
    // identity (a cluster file that does not say which cell it came from is not
    // traceable at all).
    if(extraHead)head.push(extraHead);
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
    var nGhosts=0;
    MD.atoms.forEach(function(a){if(isGhostAtom(a))nGhosts++});
    // How each format spells a ghost, and how it spells one that stands for a
    // named element. The two programs that have a ghost concept at all take the
    // element with them - Gaussian's C-Bq and ORCA's C: both mean "the basis
    // functions of a carbon, no nucleus" - and every other format here can only
    // say "a point", so it gets the dummy spelling. That is deliberate: a ghost
    // must never come back as a real atom, and X (or Du, q, *, XX) is
    // the one thing those readers agree is not one.
    var xyzEl=function(a){return isGhostAtom(a)?'X':a.element};
    var gjfEl=function(a){return isGhostAtom(a)?(a.element==='X'?'Bq':a.element+'-Bq'):a.element};
    var orcaEl=function(a){return isGhostAtom(a)?(a.element==='X'?'X':a.element+':'):a.element};
    var xyz=MD.atoms.length+'\\n'+saveTitle+'\\n';
    MD.atoms.forEach(function(a){xyz+=xyzEl(a)+'  '+a.x.toFixed(6)+'  '+a.y.toFixed(6)+'  '+a.z.toFixed(6)+'\\n'});
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
        var el=gjfEl(a); while(el.length<2) el+=' ';
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

    var coord='$coord\\n';
    var ANG_TO_BOHR=1.8897259886;
    // TURBOMOLE's dummy atom is q in $coord; it has no way to say which
    // element a ghost stands for, so an elemental ghost is written the same way.
    MD.atoms.forEach(function(a){coord+='  '+(a.x*ANG_TO_BOHR).toFixed(8)+'  '+(a.y*ANG_TO_BOHR).toFixed(8)+'  '+(a.z*ANG_TO_BOHR).toFixed(8)+' '+(isGhostAtom(a)?'q':a.element.toLowerCase())+'\\n'});
    coord+='$end\\n';
    if(chrg!==0)coord+='$chrg '+chrg+'\\n';
    if(mult!==1)coord+='$spin '+(mult-1)/2+'\\n';

    // ORCA's first line is a comment line, so the provenance is one there and
    // the keywords keep their own line below it. ORCA spells a ghost with a
    // colon after the element and a dummy with X.
    var orcaInp='# '+saveTitle+'\\n! B3LYP def2-SVP\\n\\n* xyz '+chrg+' '+mult+'\\n';
    MD.atoms.forEach(function(a){orcaInp+=orcaEl(a)+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+'\\n'});
    orcaInp+='*\\n';

    // Tripos knows one dummy: atom type Du, which is also the only field a
    // MOL2 atom line can carry it in.
    var mol2='@<TRIPOS>MOLECULE\\n'+saveTitle+'\\n'+MD.atoms.length+' '+MD.bonds.length+' 0 0 0\\nSMALL\\nNO_CHARGES\\n\\n';
    mol2+='@<TRIPOS>ATOM\\n';
    MD.atoms.forEach(function(a,i){
        if(isGhostAtom(a)){mol2+=(i+1)+' X'+(i+1)+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+' Du 1 UNK 0.000\\n';return}
        mol2+=(i+1)+' '+a.element+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+' '+a.element+' 1 UNK 0.000\\n';
    });
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
        // MDL's dummy is the * wildcard atom (a ghost has no element to name
        // and no nuclear charge to write, which is exactly what * says).
        var sym=isGhostAtom(a)?'*':a.element;
        var sx=(a.x*10).toFixed(4);
        var sy=(a.y*10).toFixed(4);
        var sz=(a.z*10).toFixed(4);
        while(sx.length<10)sx=' '+sx;
        while(sy.length<10)sy=' '+sy;
        while(sz.length<10)sz=' '+sz;
        mol+=sx+sy+sz+' '+sym+' 0  0  0  0  0  0  0  0  0  0  0  0\\n';
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
    // A ghost is written with the element column X (the standard unknown/dummy
    // symbol) and the residue DUM, so a PDB reader shows a dummy instead of
    // reading a ghost as the atom it stands for.
    var pdb=provPdbTitle(saveTitle);
    MD.atoms.forEach(function(a,i){
        var gh=isGhostAtom(a);
        var serial=(i+1).toString();while(serial.length<5)serial=' '+serial;
        var sym=gh?'X':a.element;
        var name;
        if(sym.length===1){name=' '+sym+'  '}
        else{name=sym+'  '}
        var sx=a.x.toFixed(3);while(sx.length<8)sx=' '+sx;
        var sy=a.y.toFixed(3);while(sy.length<8)sy=' '+sy;
        var sz=a.z.toFixed(3);while(sz.length<8)sz=' '+sz;
        var el=sym.charAt(0).toUpperCase();
        if(sym.length>1)el+=sym.charAt(1).toLowerCase();
        while(el.length<2)el=' '+el;
        pdb+=(gh?'HETATM':'ATOM  ')+serial+name+(gh?'DUM':'UNK')+' A   1    '+sx+sy+sz+'  1.00  0.00          '+el+'\\n';
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
    // MOPAC's dummy atom is XX (X on newer versions); a ghost has no parameters
    // to run with, so it is only a geometry point there.
    MD.atoms.forEach(function(a){
        var sx=a.x.toFixed(5);var sy=a.y.toFixed(5);var sz=a.z.toFixed(5);
        mopac+=(isGhostAtom(a)?'XX':a.element)+' '+sx+' 1 '+sy+' 1 '+sz+' 1\\n';
    });
    mopac+='\\n';

    // ===== CIF =====
    // What a saved CIF states, and what it must not: the file carries the *unit
    // cell*, which is what every edit is written back into and what every
    // rebuild redraws from, with every symmetry operation the original file
    // declared already expanded into it. Writing that expanded cell under the
    // original operation list is the mixed body that expands a second time on
    // the next read — a P2(1)/c file's 147 sites recovering as 294 atoms, and
    // 294 on the read after that — so the export declares P1 with the identity
    // operation and says which body it is in _audit_creation_method. The
    // fractional coordinates are wrapped into [0,1) for the same reason: an
    // atom moved across a face (see Local Crystal Modifications) is the same
    // site as the wrapped one, and a reader that folds coordinates itself must
    // not be handed a value that means a different site than it says.
    var cifContent='';
    if(CRY){
        var cifVol=crystalCellVolume();
        var origSg=CRY.spaceGroup?String(CRY.spaceGroup).split("'").join(''):'';
        cifContent='data_modified\\n';
        cifContent+='_audit_creation_method \\'MV: full unit cell, symmetry expanded';
        if(origSg)cifContent+=' (from '+origSg+')';
        cifContent+=', written as P1\\'\\n';
        // The cell's composition, in Hill order (the same table the info panel
        // shows), one count per site. The occupancy column says how much of
        // each site is filled, normalized to 1 for a group shown as fully
        // occupied (crystalOccOf); All writes the file's values unchanged.
        cifContent+='_chemical_formula_sum \\''+crystalFormulaSum(crystalVisibleBaseAtoms())+'\\'\\n';
        cifContent+='_cell_length_a '+CRY.a.toFixed(6)+'\\n';
        cifContent+='_cell_length_b '+CRY.b.toFixed(6)+'\\n';
        cifContent+='_cell_length_c '+CRY.c.toFixed(6)+'\\n';
        cifContent+='_cell_angle_alpha '+CRY.alpha.toFixed(6)+'\\n';
        cifContent+='_cell_angle_beta '+CRY.beta.toFixed(6)+'\\n';
        cifContent+='_cell_angle_gamma '+CRY.gamma.toFixed(6)+'\\n';
        cifContent+='_cell_volume '+cifVol.toFixed(6)+'\\n';
        cifContent+='_symmetry_space_group_name_H-M \\'P 1\\'\\n';
        cifContent+='loop_\\n_symmetry_equiv_pos_as_xyz\\n\\'x, y, z\\'\\n';
        cifContent+='loop_\\n_atom_site_label\\n_atom_site_type_symbol\\n_atom_site_fract_x\\n_atom_site_fract_y\\n_atom_site_fract_z\\n_atom_site_occupancy\\n_atom_site_site_symmetry_multiplicity\\n';
        // The save follows the disorder selection: the file describes the
        // alternative the picture shows, not every alternative at once.
        var cifBase=crystalVisibleBaseAtoms();
        var cifAtomOf={};
        var cifLabelOf={};
        cifBase.forEach(function(ba,i){
            var f=crystalFracWrapped(ba.x,ba.y,ba.z);
            var label=ba.element+(i+1);
            var occ=crystalOccOf(ba).toFixed(3);
            cifAtomOf[ba.baseIdx!=null?ba.baseIdx:i]=ba;
            cifLabelOf[ba.baseIdx!=null?ba.baseIdx:i]=label;
            // P1 has one operation, so every site is its own multiplicity.
            cifContent+=label+' '+ba.element+' '+f[0].toFixed(6)+' '+f[1].toFixed(6)+' '+f[2].toFixed(6)+' '+occ+' 1\\n';
        });
        // Anisotropic displacement parameters travel with the site. They are
        // already Cartesian Angstrom squared and already rotated into the P1
        // frame, so the six numbers are written as they stand; the loop names
        // only the sites that carry a tensor, which the dictionary allows. A
        // saved file without any aniso data is byte-identical to the old one.
        var cifAnisoRows='';
        cifBase.forEach(function(ba,i){
            var u=ba.uAniso;
            if(!u||u.length<6)return;
            var ok=true;
            for(var k=0;k<6;k++){if(typeof u[k]!=='number'||!isFinite(u[k]))ok=false}
            if(!ok)return;
            cifAnisoRows+=ba.element+(i+1);
            for(var k2=0;k2<6;k2++)cifAnisoRows+=' '+u[k2].toFixed(6);
            cifAnisoRows+='\\n';
        });
        if(cifAnisoRows!==''){
            cifContent+='loop_\\n_atom_site_aniso_label\\n_atom_site_aniso_U_11\\n'+
                '_atom_site_aniso_U_22\\n_atom_site_aniso_U_33\\n_atom_site_aniso_U_12\\n'+
                '_atom_site_aniso_U_13\\n_atom_site_aniso_U_23\\n'+cifAnisoRows;
        }
        // The bond table, in the labels just written. The cell is declared P1,
        // so a bond that points into a neighbouring cell is stated with the one
        // operation and a translation (1_455 is -a): the file carries the
        // same connectivity it was drawn with, and reading it back rebuilds the
        // same base bonds instead of a face cut. A bond whose endpoint is not
        // on screen (a disorder alternative that is not selected) is left out
        // with the site it names.
        var cifBonds=CRY.baseBonds||[];
        var cifAnyValence=false;
        cifBonds.forEach(function(bb){if(bb.valence>0)cifAnyValence=true});
        var cifBondRows='';
        cifBonds.forEach(function(bb){
            var l1=cifLabelOf[bb.atom1],l2=cifLabelOf[bb.atom2];
            var ba1=cifAtomOf[bb.atom1],ba2=cifAtomOf[bb.atom2];
            if(!l1||!l2||!ba1||!ba2)return;
            var s=bb.shift||[0,0,0];
            if(Math.abs(s[0])>4||Math.abs(s[1])>4||Math.abs(s[2])>4)return;
            var code=(s[0]||s[1]||s[2])?'1_'+(s[0]+5)+(s[1]+5)+(s[2]+5):'.';
            var lvC=CRY.latticeVectors;
            var dx=ba2.x+s[0]*lvC[0][0]+s[1]*lvC[1][0]+s[2]*lvC[2][0]-ba1.x;
            var dy=ba2.y+s[0]*lvC[0][1]+s[1]*lvC[1][1]+s[2]*lvC[2][1]-ba1.y;
            var dz=ba2.z+s[0]*lvC[0][2]+s[1]*lvC[1][2]+s[2]*lvC[2][2]-ba1.z;
            var dist=Math.sqrt(dx*dx+dy*dy+dz*dz);
            var ord=(bb.order>0?bb.order:1);
            cifBondRows+=l1+' '+l2+' '+dist.toFixed(4)+' '+code+' '+ord;
            if(cifAnyValence)cifBondRows+=' '+(bb.valence>0?bb.valence:'?');
            cifBondRows+='\\n';
        });
        if(cifBondRows!==''){
            cifContent+='loop_\\n_geom_bond_atom_site_label_1\\n_geom_bond_atom_site_label_2\\n'+
                '_geom_bond_distance\\n_geom_bond_site_symmetry_2\\n_geom_bond_order\\n';
            if(cifAnyValence)cifContent+='_geom_bond_valence\\n';
            cifContent+=cifBondRows;
        }
    }

    // ===== VASP POSCAR =====
    // The old export wrote the *displayed* atoms against the *unit cell's*
    // lattice, so every atom of a supercell sat outside the box it declared
    // (fractions of 2.0 and more) and the file described no structure at all.
    // The atoms and the lattice now come from the same object — the base cell
    // that CIF/VESTA write and that every edit is committed into — and the
    // comment line says so, so a POSCAR of a 2x2x2 view is the unit cell and
    // not a silently mislabelled supercell.
    var vaspContent='';
    if(CRY){
        var lv=CRY.latticeVectors;
        var base=crystalVisibleBaseAtoms();
        var elemOrder=[];
        var elemCounts={};
        var elemAtoms={};
        base.forEach(function(ba){
            if(elemCounts[ba.element]===undefined){
                elemOrder.push(ba.element);elemCounts[ba.element]=0;elemAtoms[ba.element]=[];
            }
            elemCounts[ba.element]++;
            elemAtoms[ba.element].push(ba);
        });
        // POSCAR assigns elements through the order of the symbol line and the
        // counts line, so the coordinates must be written grouped the same way.
        // A base cell in site order (Na, Cl, Na) against a two-name header
        // would relabel the second Na site as Cl on the next read.
        var vaspAtoms=[];
        elemOrder.forEach(function(e){elemAtoms[e].forEach(function(ba){vaspAtoms.push(ba)})});
        // Selective dynamics is written only when some atom carries flags; an
        // atom without any is free, which is T T T.
        var anyConstraint=base.some(function(ba){return !!ba.constraint});
        var vaspTitle=provTitle(MD.title||'Modified structure',true,'unit cell, '+base.length+' atoms (not the displayed supercell)');
        vaspContent=vaspTitle+'\\n';
        vaspContent+='1.0\\n';
        for(var ii=0;ii<3;ii++){
            vaspContent+='  '+lv[ii][0].toFixed(8)+'  '+lv[ii][1].toFixed(8)+'  '+lv[ii][2].toFixed(8)+'\\n';
        }
        vaspContent+=elemOrder.join('  ')+'\\n';
        vaspContent+=elemOrder.map(function(e){return elemCounts[e]}).join('  ')+'\\n';
        if(anyConstraint)vaspContent+='Selective dynamics\\n';
        vaspContent+='Direct\\n';
        vaspAtoms.forEach(function(ba){
            var f=crystalFracWrapped(ba.x,ba.y,ba.z);
            var line='  '+f[0].toFixed(8)+'  '+f[1].toFixed(8)+'  '+f[2].toFixed(8);
            if(anyConstraint){
                var flags=(ba.constraint&&/^[TF] [TF] [TF]$/.test(ba.constraint))?ba.constraint:'T T T';
                line+='  '+flags;
            }
            vaspContent+=line+'\\n';
        });
    }

    // ===== Gaussian cube =====
    // A negative atom count makes a cube file announce orbital records it does
    // not carry, and a 2x2x2 grid announces a volume that is not there either;
    // the file has to be a real one to be read: a positive atom count (which
    // also states the units are Bohr), the unit cell as the grid's box, and a
    // grid the caller can size. The values are zeros because the viewer carries
    // no volumetric data — the file is a geometry + grid container, not a
    // density — and a reader may sample it as such. It is built when the user
    // presses Save, because the grid size is chosen in the dialog.
    var vestaContent='';
    if(CRY){
        var vestaTitle=provTitle(MD.title||'Modified structure',true,'unit cell, '+crystalVisibleBaseAtoms().length+' atoms');
        vestaContent='#VESTA_FORMAT_VERSION 3.5.4\\n\\n\\nCRYSTAL\\n\\nTITLE\\n'+vestaTitle+'\\n\\nGROUP\\n    1   1 P 1\\n\\nSYMOP\\n 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1\\n-1.0 -1.0 -1.0  1 0 0  0 0 0  0 0 0\\n\\nTRANM 0\\n 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1\\n\\nLTRANSL\\n -1\\n 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\\n\\nLORIENT\\n -1   0   0   0   0\\n 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000\\n 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000\\n\\nLMATRIX\\n 1.000000  0.000000  0.000000  0.000000\\n 0.000000  1.000000  0.000000  0.000000\\n 0.000000  0.000000  1.000000  0.000000\\n 0.000000  0.000000  0.000000  1.000000\\n 0.000000  0.000000  0.000000\\n\\nCELLP\\n';
        vestaContent+=' '+CRY.a.toFixed(6)+'  '+CRY.b.toFixed(6)+'  '+CRY.c.toFixed(6)+'  '+CRY.alpha.toFixed(6)+'  '+CRY.beta.toFixed(6)+'  '+CRY.gamma.toFixed(6)+'\\n';
        vestaContent+=' 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\\n\\nSTRUC\\n';
        crystalVisibleBaseAtoms().forEach(function(ba,i){
            var f=crystalFracWrapped(ba.x,ba.y,ba.z);
            var label=ba.element+(i+1);
            var occ=crystalOccOf(ba).toFixed(4);
            var idxStr=''+(i+1);while(idxStr.length<3)idxStr=' '+idxStr;
            var elStr=ba.element;while(elStr.length<3)elStr=elStr+' ';
            var labelStr=label;while(labelStr.length<4)labelStr=labelStr+' ';
            vestaContent+=idxStr+' '+elStr+' '+labelStr+' '+occ+'   '+f[0].toFixed(6)+'   '+f[1].toFixed(6)+'   '+f[2].toFixed(6)+'    2i       1\\n';
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
    // A ghost survives the trip into every format, but not in the same words -
    // and in the formats that have no ghost concept it can only survive as a
    // dummy point. Saying which is which is the difference between a counterpoise
    // input and a file that quietly contains an extra atom.
    var ghostNote=nGhosts>0
        ?'<div class="current-val"><b>'+nGhosts+' ghost atom'+(nGhosts===1?'':'s')+'</b>: '+
            'written as <code>C-Bq</code> in the Gaussian input and <code>C:</code> in the ORCA input '+
            '(element kept), and as a dummy point - <code>X</code> / <code>Du</code> / <code>q</code> / '+
            '<code>*</code> / <code>XX</code> - in the formats that have no ghost atom.</div>'
        :'';
    // A crystal's two families of files do different jobs, and the dialog says
    // so before it writes either: the crystal formats describe the unit cell
    // (the cell the session edits, symmetry expanded, P1), while the cluster
    // formats cut what is on screen out as a molecule — the file a QM job
    // starts from, and the only way to get one cage out of a framework.
    var crystalNote=CRY
        ?'<div class="current-val">Crystal formats (CIF / POSCAR / VESTA / cube) write the <b>unit cell</b> - '+
            'every symmetry operation expanded into one cell, declared P1. '+
            'The cluster formats (XYZ / GJF - displayed cluster) write <b>what is on screen</b>, '+
            'cut to the bounds below, as a molecule.</div>'
        :'';
    showModal('<h3>Save File</h3>'+vdNote+ghostNote+crystalNote+
        '<label>Format:</label><select id="m-fmt">'+
        // The molecular formats keep their old order (plain XYZ is still the
        // default selection), and the cluster formats sit next to them because
        // they are the same question asked of a crystal.
        '<option value="xyz">XYZ (.xyz)</option>'+
        '<option value="gjf">Gaussian Input (.gjf)</option>'+
        (CRY?'<option value="xyz-cluster">XYZ - displayed cluster (.xyz)</option>':'')+
        (CRY?'<option value="gjf-cluster">Gaussian Input - displayed cluster (.gjf)</option>':'')+
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
        (CRY?'<div class="modal-row" id="m-keepmol-row" style="display:none"><label style="min-width:0">Cluster:</label>'+
            '<input type="checkbox" id="m-keepmol" checked style="flex:0 0 auto;width:auto">'+
            '<span class="current-val" style="margin:0">keep whole molecules (untick to cut strict at the bounds)</span></div>':'')+
        (CRY?'<div class="modal-row" id="m-grid-row" style="display:none"><label>Cube grid:</label>'+
            '<input type="number" id="m-grid" value="'+CUBE_GRID_DEFAULT+'" step="1" min="2" max="'+CUBE_GRID_MAX+'">'+
            '<span class="current-val" style="margin:0">points per cell axis (all-zero values)</span></div>':'')+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Save</button></div>',null);
    // The two option rows belong to specific formats, so they appear with them:
    // a cluster control shown over a CIF export would be a switch that does
    // nothing, and the grid shown everywhere would read as a CIF option.
    var fmtSelEl=document.getElementById('m-fmt');
    function syncSaveRows(){
        var kr=document.getElementById('m-keepmol-row');
        var gr=document.getElementById('m-grid-row');
        var f=fmtSelEl?fmtSelEl.value:'';
        if(kr&&kr.style)kr.style.display=(f==='xyz-cluster'||f==='gjf-cluster')?'flex':'none';
        if(gr&&gr.style)gr.style.display=(f==='cube')?'flex':'none';
    }
    if(fmtSelEl)fmtSelEl.addEventListener('change',syncSaveRows);
    syncSaveRows();
    document.getElementById('m-ok').addEventListener('click',function(){
        var fmt=document.getElementById('m-fmt').value;
        var content,ext;
        try{
        switch(fmt){
            case 'gjf':content=gjf;ext='.gjf';break;
            case 'cif':content=cifContent;ext='.cif';break;
            case 'vasp':content=vaspContent;ext='.vasp';break;
            // The grid is the one thing about the cube the dialog lets a user
            // size, so it is read here rather than at open time. A grid whose
            // value block would be unreasonably large is refused before a file
            // is written, with the limit named, instead of being clamped to a
            // size the user did not ask for.
            case 'cube':{
                var gEl=document.getElementById('m-grid');
                var nGrid=gEl?Math.floor(Number(gEl.value)):CUBE_GRID_DEFAULT;
                if(!(nGrid>=2))nGrid=CUBE_GRID_DEFAULT;
                if(nGrid>CUBE_GRID_MAX){
                    alert('Cube grid '+nGrid+' is too large: the limit is '+CUBE_GRID_MAX+' points per axis ('+(CUBE_GRID_MAX*CUBE_GRID_MAX*CUBE_GRID_MAX)+' values). Reduce it and save again.');
                    return;
                }
                content=buildCubeContent(nGrid);ext='.cube';break;
            }
            case 'vesta':content=vestaContent;ext='.vesta';break;
            case 'coord':content=coord;ext='.coord';break;
            case 'inp':content=orcaInp;ext='.inp';break;
            case 'mol2':content=mol2;ext='.mol2';break;
            case 'mol':content=mol;ext='.mol';break;
            case 'pdb':content=pdb;ext='.pdb';break;
            case 'mop':content=mopac;ext='.mop';break;
            case 'xyz-cluster':content=buildClusterXyz(keepWholeMolecules());ext='.xyz';break;
            case 'gjf-cluster':content=buildClusterGjf(keepWholeMolecules());ext='.gjf';break;
            default:content=xyz;ext='.xyz';
        }
        }catch(err){
            alert('Failed to generate save content: '+(err&&err.message?err.message:String(err)));
            return;
        }
        vscodeApi.postMessage({command:'saveFile',content:content,suggestedName:'molecule_modified'+ext,filePath:MD.filePath||''});
        hideModal();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal()});
}

// ===== Crystal export: the unit cell, and the cluster cut out of it =====
//
// Every crystal format in the save dialog answers one of two questions, and
// each answer is written by the functions below instead of a second copy of
// the geometry arithmetic:
//
//   "what is the cell?"   CIF / POSCAR / VESTA / cube — baseAtoms, each site
//                         wrapped into [0,1), the cell's own lattice, and the
//                         declaration (P1, audit line, comment line) that says
//                         the symmetry is already expanded;
//   "what is on screen?"  xyz-cluster / gjf-cluster — the displayed atoms cut
//                         to the supercell bounds, as the molecule a QM job
//                         starts from.
//
// The two are deliberately different files: a cell is periodic and a cluster is
// not, and a format that says one while containing the other is the bug E2/E3
// were about (a POSCAR listing 2x2x2 atoms against a 1x1x1 lattice, a CIF
// listing an expanded cell under the original operation list).

// The cube's grid is the only export setting that changes the file's size, so
// it has an explicit ceiling: 64 points per axis is 262 144 zeros, a few
// megabytes of text. Above that a file is refused rather than clamped.
var CUBE_GRID_DEFAULT=40;
var CUBE_GRID_MAX=64;
var CUBE_BOHR_PER_ANG=1.8897259886;

/** The unit cell's volume (the lattice determinant). */
function crystalCellVolume(){
    if(!CRY||!CRY.latticeVectors)return 0;
    var lv=CRY.latticeVectors;
    var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])
          -lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])
          +lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
    return Math.abs(det);
}
/** The cell's composition in Hill order, spaced the way a CIF writes it. */
function crystalFormulaSum(atoms){
    var f=formulaFromAtoms(atoms||[]);
    if(!f)return '';
    var re=/([A-Z][a-z]?)(\\d+)?/g,m,out='';
    while((m=re.exec(f))!==null){
        if(out)out+=' ';
        out+=m[1]+(m[2]||'');
    }
    return out;
}
/**
 * Fractional coordinates wrapped into [0,1).
 *
 * Wrapping is not cosmetic: an atom moved across a face is the same site as its
 * wrapped position, and a file that lists 1.02 hands a reader a coordinate it
 * has to fold itself — while a value that rounds to 1.000000 at six decimals is
 * written as 0 instead of as the upper face, which is the next cell.
 */
function crystalFracWrapped(x,y,z){
    var inv=CRY_INV;
    var f=[inv[0][0]*x+inv[0][1]*y+inv[0][2]*z,
           inv[1][0]*x+inv[1][1]*y+inv[1][2]*z,
           inv[2][0]*x+inv[2][1]*y+inv[2][2]*z];
    for(var i=0;i<3;i++){
        var w=f[i]-Math.floor(f[i]);
        if(w>=1-5e-7)w=0;
        f[i]=w;
    }
    return f;
}
/**
 * The displayed region as a molecule: atoms plus the bonds among them.
 *
 * keepWhole=true is what the picture already shows — the builder places a
 * finite molecule complete once its centroid is inside the bounds, so its atoms
 * can sit outside the box and still belong to the cage the user cut out.
 * keepWhole=false is the strict reading, atom by atom against the bounds, which
 * is what a fragment cut for a QM job often needs to be.
 *
 * A bond that crosses a cell face (a bond with a shift) is dropped: it is the
 * stub the renderer draws towards an image that is not in the file, and writing
 * it as a bond between two atoms of the cluster would invent a length across
 * the cell.
 */
function crystalCluster(keepWhole){
    var inv=CRY_INV;
    var map={},atoms=[],bonds=[];
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(!keepWhole){
            var fx=inv[0][0]*a.x+inv[0][1]*a.y+inv[0][2]*a.z;
            var fy=inv[1][0]*a.x+inv[1][1]*a.y+inv[1][2]*a.z;
            var fz=inv[2][0]*a.x+inv[2][1]*a.y+inv[2][2]*a.z;
            if(fx<boundary.aMin-1e-9||fx>boundary.aMax+1e-9)continue;
            if(fy<boundary.bMin-1e-9||fy>boundary.bMax+1e-9)continue;
            if(fz<boundary.cMin-1e-9||fz>boundary.cMax+1e-9)continue;
        }
        map[i]=atoms.length;
        atoms.push({element:a.element,ghost:a.ghost,x:a.x,y:a.y,z:a.z,index:atoms.length});
    }
    MD.bonds.forEach(function(b){
        if(b.shift&&(b.shift[0]||b.shift[1]||b.shift[2]))return;
        var a1=map[b.atom1],a2=map[b.atom2];
        if(a1===undefined||a2===undefined)return;
        bonds.push({atom1:a1,atom2:a2,order:b.order});
    });
    return {atoms:atoms,bonds:bonds};
}
/** The cell a cluster was cut from and the bounds that cut it, for its title. */
function crystalClusterTitle(nAtoms,keepWhole){
    var c=CRY;
    var cell='cell '+c.a.toFixed(3)+'/'+c.b.toFixed(3)+'/'+c.c.toFixed(3)+' A '+
        c.alpha.toFixed(2)+'/'+c.beta.toFixed(2)+'/'+c.gamma.toFixed(2);
    var bnd='bounds a['+boundary.aMin+','+boundary.aMax+'] b['+boundary.bMin+','+boundary.bMax+'] c['+boundary.cMin+','+boundary.cMax+']';
    return 'cluster '+cell+', '+bnd+', '+(keepWhole?'whole molecules':'cut at the bounds')+', '+nAtoms+' atoms';
}
/** The cluster-dialog switch: whole molecules unless the user unticked it. */
function keepWholeMolecules(){
    var el=document.getElementById('m-keepmol');
    return el?!!el.checked:true;
}
function buildClusterXyz(keepWhole){
    var cl=crystalCluster(keepWhole);
    var title=provTitle(MD.title||'Modified structure',false,crystalClusterTitle(cl.atoms.length,keepWhole));
    var out=cl.atoms.length+'\\n'+title+'\\n';
    cl.atoms.forEach(function(a){
        out+=(isGhostAtom(a)?'X':a.element)+'  '+a.x.toFixed(6)+'  '+a.y.toFixed(6)+'  '+a.z.toFixed(6)+'\\n';
    });
    return out;
}
function buildClusterGjf(keepWhole){
    var cl=crystalCluster(keepWhole);
    var title=provTitle(MD.title||'Modified structure',false,crystalClusterTitle(cl.atoms.length,keepWhole));
    var chrg=MD.charge||0,mult=MD.multiplicity||1;
    // A crystal has no Gaussian route of its own to preserve, so the cluster
    // gets the same default header a fresh structure gets.
    var out='%chk=cluster.chk\\n%mem=4GB\\n%nproc=4\\n# B3LYP/6-31G(d)\\n\\n'+title+'\\n\\n'+chrg+' '+mult+'\\n';
    cl.atoms.forEach(function(a){
        var el=isGhostAtom(a)?(a.element==='X'?'Bq':a.element+'-Bq'):a.element;
        while(el.length<2)el+=' ';
        var sx=a.x.toFixed(8);while(sx.length<17)sx=' '+sx;
        var sy=a.y.toFixed(8);while(sy.length<17)sy=' '+sy;
        var sz=a.z.toFixed(8);while(sz.length<17)sz=' '+sz;
        out+=' '+el+' '+sx+' '+sy+' '+sz+'\\n';
    });
    out+='\\n';
    // One connect line per atom, each bond written once (on the line of the
    // lower-numbered atom) — the convention the molecular GJF writer uses, so
    // a cluster input and an ordinary input read the same way.
    var conn=[];
    for(var ci=0;ci<cl.atoms.length;ci++)conn.push([]);
    cl.bonds.forEach(function(b){
        if(b.atom1<b.atom2)conn[b.atom1].push(b.atom2+1,b.order.toFixed(1));
        else if(b.atom2<b.atom1)conn[b.atom2].push(b.atom1+1,b.order.toFixed(1));
    });
    for(var ai=0;ai<cl.atoms.length;ai++){
        out+=[ai+1].concat(conn[ai]).join(' ')+'\\n';
    }
    out+='\\n';
    return out;
}
/**
 * A real Gaussian cube of the unit cell: the atoms, the cell as the grid's box,
 * and a grid of zeros the caller sizes.
 *
 * A positive atom count is written (a negative one announces orbital records
 * this file does not carry, and readers then consume the grid as that record),
 * and it also states the units: cube files count in Bohr unless the count is
 * negative, so every coordinate here is converted rather than written in Å
 * under a header that says otherwise. The values are zeros — the viewer has no
 * density to sample — and the header's second line says so, so a consumer can
 * tell "a grid that happens to be flat" from "a density that came out zero".
 */
function buildCubeContent(nGrid){
    var lv=CRY.latticeVectors;
    // The cube follows the disorder selection like every other crystal export,
    // so a QM job starts from the alternative that is on screen.
    var base=crystalVisibleBaseAtoms();
    var title=provTitle(MD.title||'Modified structure',true,'unit cell, '+base.length+' atoms, '+nGrid+'^3 zero grid');
    var out='Molecular Viewer export\\n'+title+'\\n';
    out+=base.length+'  0.000000  0.000000  0.000000\\n';
    for(var ax=0;ax<3;ax++){
        out+=nGrid+'  '+(lv[ax][0]/nGrid*CUBE_BOHR_PER_ANG).toFixed(8)+
            '  '+(lv[ax][1]/nGrid*CUBE_BOHR_PER_ANG).toFixed(8)+
            '  '+(lv[ax][2]/nGrid*CUBE_BOHR_PER_ANG).toFixed(8)+'\\n';
    }
    base.forEach(function(ba){
        var f=crystalFracWrapped(ba.x,ba.y,ba.z);
        var x=(f[0]*lv[0][0]+f[1]*lv[1][0]+f[2]*lv[2][0])*CUBE_BOHR_PER_ANG;
        var y=(f[0]*lv[0][1]+f[1]*lv[1][1]+f[2]*lv[2][1])*CUBE_BOHR_PER_ANG;
        var z=(f[0]*lv[0][2]+f[1]*lv[1][2]+f[2]*lv[2][2])*CUBE_BOHR_PER_ANG;
        out+=(AN[ba.element]||0)+'  0.000000  '+x.toFixed(8)+'  '+y.toFixed(8)+'  '+z.toFixed(8)+'\\n';
    });
    var nVals=nGrid*nGrid*nGrid;
    for(var vi=0;vi<nVals;vi++){
        out+='0.0000e+00';
        if((vi+1)%6===0||vi===nVals-1)out+='\\n';else out+='  ';
    }
    return out;
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

// ===== Display style (ball and stick / space-filling / stick) =====
// The control sits in the View card, but the state it sets does not belong to
// View: the scene build and the bond builder read it, so a structure keeps the
// style it was given when the toolbar moves on to another tool, when the frame
// or the molecule changes, and when an undo rewrites the structure. One function
// keeps the control in step with the state, so a value the page does not
// recognise can never leave the two disagreeing.
var styleSelectEl=document.getElementById('style-select');
function syncStyleSelect(){
    if(styleSelectEl&&styleSelectEl.value!==displayStyle)styleSelectEl.value=displayStyle;
}
// A style *is* the geometry - sphere radii live in the (pooled) sphere geometry
// and a space-filling model draws no sticks - so it is applied by rebuilding the
// scene, the same path a frame switch or an edit takes. That is also what keeps
// everything else in step with it: the selection highlight, the comparison's
// right half, the surface patches and the info panel are all parts of that one
// rebuild.
function setDisplayStyle(style){
    // hasOwnProperty, not a plain lookup: 'constructor' and 'toString' are on
    // every object, and a style is only ever one of the three keys.
    if(!Object.prototype.hasOwnProperty.call(DISPLAY_STYLES,style)||style===displayStyle){syncStyleSelect();return}
    displayStyle=style;
    syncStyleSelect();
    rebuildScene();
    // In Diff mode the right-hand structure is a model of its own, built by
    // buildDiffModel rather than by rebuildScene, so it is restyled here — with
    // the marks it is already wearing, so the comparison does not lose them.
    if(diffMode)rebuildDiffModel();
    // The panels keep clear of the part of the view the structure occupies, and
    // that box is measured from the drawn spheres: a space-filling model is
    // several times the footprint of the ball model it replaced, so the layout
    // is taken again rather than left around a structure that changed size.
    layoutPanels();
}
if(styleSelectEl){
    syncStyleSelect();
    styleSelectEl.addEventListener('change',function(){setDisplayStyle(this.value)});
}

// ===== Thermal ellipsoids: the View card's switch =====
// The switch is offered only when the structure on screen states at least one
// usable ADP; a molecule from a coordinate file has none, and a checked box
// that draws the same spheres would be a claim the file does not support.
// rebuildScene calls this after every build, so deleting the last anisotropic
// site (or removing a disorder alternative that held it) turns the switch off
// together with the picture.
function syncOrtepControl(){
    var chk=document.getElementById('ortep-chk');
    var line=document.getElementById('ortep-line');
    if(!chk)return;
    var any=false;
    for(var i=0;i<MD.atoms.length;i++){
        var a=MD.atoms[i];
        if(a&&isAnisoUsable(a.uAniso)){any=true;break}
    }
    if(!any)ortepEnabled=false;
    // A polyhedron's centre marker is not a ball for an ellipsoid to replace,
    // so the switch is held off while that style is on. The choice itself is
    // kept: switching back to Ball & Stick brings the ellipsoids back.
    var poly=displayStyle==='polyhedra';
    chk.disabled=!any||poly;
    chk.checked=ortepEnabled&&!poly;
    if(line){
        if(any&&!poly)line.classList.remove('off');else line.classList.add('off');
        line.title=!any?'This structure states no anisotropic displacement parameters':
            poly?'Thermal ellipsoids are not drawn in Polyhedra style; the chosen state returns in Ball & Stick':
            'Draw atoms that carry anisotropic displacement parameters as 50% probability ellipsoids';
    }
}
var ortepChkEl=document.getElementById('ortep-chk');
if(ortepChkEl){
    ortepChkEl.addEventListener('change',function(){
        ortepEnabled=!!this.checked;
        rebuildScene();
        layoutPanels();
        modeInfoEl.textContent=ortepEnabled?
            'Thermal ellipsoids: 50% probability surfaces':
            (MODE_INFO[currentMode]||currentMode);
    });
}

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
    // frozen for the drag, so both are a re-placement of existing objects. A
    // rotation bake turns the ring planes with the atoms, so the aromatic sides
    // are re-planned first — the same rule the fragment-rotation drag follows,
    // and a translation bake leaves every direction (and the plan) untouched.
    if(bakeKind==='rot')planAromaticDashes();
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
    var key=(diffMode||CRY)?('x:'+bakeUnavailableReason()):'y';
    if(hint.__hpKey===key)return;
    hint.__hpKey=key;
    if(diffMode||CRY){
        hint.textContent=bakeUnavailableReason()+'.';
    }else{
        hint.innerHTML='Hold <b>Alt</b> and drag: left rotates the structure, right/middle moves it. '+
            'The atom coordinates change, so Save writes the new geometry.';
    }
    // The two hints are not the same height, and the card is laid out from its
    // measured height: a swap without a re-layout left the column budget (and
    // the cards below it) sized for the old text until something else moved.
    layoutPanels();
}
updateBakeHint();

// ===== Picking =====
// One door for "which atom is under the cursor": every click — the selection,
// the three measurements, Replace Atom, Add H's anchor, Move Atoms' grab,
// Delete — and the hover tooltip come through here, so the two can never
// disagree about what is being pointed at.
//
// The spheres are the targets. The stick model, however, draws an atom no wider
// than the tubes meeting at it (that is what the model is), which leaves a joint
// a few pixels across once a structure is large enough to be zoomed out: a
// 0.24 A joint is about 3 px on a 25 A cage, and that is not a target a hand can
// aim at. A press or hover that misses every sphere is therefore re-tested
// against a cone around each atom whose radius is fixed in *pixels*, so the
// precision the pointer needs is the same at every zoom level. It applies to the
// stick model only: the default balls and the van der Waals spheres are their own
// targets, and picking in them is exactly what it has always been.
var PICK_SLOP_PX=6;
var _pickInv=new THREE.Matrix4();
var _pickO=new THREE.Vector3();
var _pickDir=new THREE.Vector3();
function stickPick(atoms,meshes,cx,cy,cz,group){
    if(displayStyle!=='stick'||!group||!atoms||atoms.length===0)return -1;
    var h=canvas.clientHeight||0;
    if(!(h>1))return -1;
    // Radius of the cone per unit of depth: the canvas height spans
    // 2*tan(fov/2) at depth 1, so PICK_SLOP_PX pixels are this fraction of it.
    var perDepth=(2*Math.tan((camera.fov||60)*Math.PI/360)*PICK_SLOP_PX)/h;
    // The ray is taken into the group's own frame, where the atoms are, so the
    // test is the plain distance from a point to a line. The group's transform is
    // a rotation and a translation, which preserve both lengths and the parameter
    // along the ray, so the answer is the one the world-space ray would give.
    group.updateWorldMatrix(true,false);
    _pickInv.copy(group.matrixWorld).invert();
    _pickO.copy(raycaster.ray.origin).applyMatrix4(_pickInv);
    _pickDir.copy(raycaster.ray.direction).transformDirection(_pickInv);
    var n=meshes?Math.min(atoms.length,meshes.length):atoms.length;
    var best=-1,bestT=Infinity;
    for(var i=0;i<n;i++){
        // An atom that is not drawn is not a target, which is exactly what its
        // mesh says (a VMD script hides an atom by giving it sphere scale 0).
        if(meshes&&meshes[i]&&meshes[i].visible===false)continue;
        var a=atoms[i];if(!a)continue;
        var ox=a.x-cx-_pickO.x,oy=a.y-cy-_pickO.y,oz=a.z-cz-_pickO.z;
        var t=ox*_pickDir.x+oy*_pickDir.y+oz*_pickDir.z;
        if(!(t>0)||t>=bestT)continue;
        var qx=ox-_pickDir.x*t,qy=oy-_pickDir.y*t,qz=oz-_pickDir.z*t;
        var lim=t*perDepth;
        if(qx*qx+qy*qy+qz*qz<=lim*lim){best=i;bestT=t}
    }
    return best;
}

// The spheres are instanced now, so the hit test is the analytic one the ray
// always meant: the nearest sphere along the ray, in the group's own frame (a
// rotation and a translation preserve lengths and the ray parameter, exactly
// as stickPick argues). A tessellated SphereGeometry raycast would answer the
// same atom; this answers it without walking a geometry.
function pickAtomSphere(atoms,meshes,cx,cy,cz,group){
    if(!group||!atoms||!meshes||atoms.length===0)return -1;
    group.updateWorldMatrix(true,false);
    _pickInv.copy(group.matrixWorld).invert();
    _pickO.copy(raycaster.ray.origin).applyMatrix4(_pickInv);
    _pickDir.copy(raycaster.ray.direction).transformDirection(_pickInv);
    var n=Math.min(atoms.length,meshes.length);
    var best=-1,bestT=Infinity;
    for(var i=0;i<n;i++){
        var m=meshes[i];
        if(!m||m.visible===false)continue;
        var r=m.userData?m.userData.radius:0;
        if(!(r>0))continue;
        var o=m.position;
        var ox=o.x-_pickO.x,oy=o.y-_pickO.y,oz=o.z-_pickO.z;
        var t=ox*_pickDir.x+oy*_pickDir.y+oz*_pickDir.z;
        var qx=ox-_pickDir.x*t,qy=oy-_pickDir.y*t,qz=oz-_pickDir.z*t;
        var d2=qx*qx+qy*qy+qz*qz;
        if(d2>r*r)continue;
        var dt=Math.sqrt(r*r-d2);
        var t0=t-dt;
        var hit=t0>0?t0:t+dt;
        if(hit<0)continue;
        if(hit<bestT){best=i;bestT=hit}
    }
    return best;
}
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
            var pl=pickAtomSphere(MD.atoms,atomMeshes,CX,CY,CZ,moleculeGroup);
            if(pl>=0)return pl;
            var stL=stickPick(MD.atoms,atomMeshes,CX,CY,CZ,moleculeGroup);
            if(stL>=0)return stL;
        }else{
            diffActiveSide='right';
            mouse.x=((localX-halfW)/(w-halfW))*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=(w-halfW)/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            var pr=pickAtomSphere(diffData.atoms,diffAtomMeshes,diffCX,diffCY,diffCZ,diffMolGroup);
            if(pr>=0)return pr;
            var stR=stickPick(diffData.atoms,diffAtomMeshes,diffCX,diffCY,diffCZ,diffMolGroup);
            if(stR>=0)return stR;
        }
        return-1;
    }
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var hit=pickAtomSphere(MD.atoms,atomMeshes,CX,CY,CZ,moleculeGroup);
    if(hit>=0)return hit;
    return stickPick(MD.atoms,atomMeshes,CX,CY,CZ,moleculeGroup);
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
    // The hover probe and the press ask the same question through the same door,
    // so they can never disagree about which atom is pointed at - and the stick
    // model's wider pick target applies to both.
    var hi=getClickedAtom(e);
    if(hi>=0){
        var diffRight=diffMode&&diffActiveSide==='right';
        var a=diffRight?diffData.atoms[hi]:MD.atoms[hi];
        var tip=atomLabel(a,hi)+(isGhostAtom(a)?' (ghost)':'')+' ('+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4)+')';
        // A crystal atom also says where it sits in the lattice: its fractional
        // coordinates, which periodic image (i,j,k) of its base site this is,
        // the file's own site label and its occupancy. A molecular structure's
        // tooltip reads exactly what it always read — none of these fields
        // exist there, and the comparison's right half is another molecule.
        if(!diffRight&&CRY&&CRY_INV&&a&&a.baseIdx!=null&&a.cellI!=null){
            var tfx=CRY_INV[0][0]*a.x+CRY_INV[0][1]*a.y+CRY_INV[0][2]*a.z;
            var tfy=CRY_INV[1][0]*a.x+CRY_INV[1][1]*a.y+CRY_INV[1][2]*a.z;
            var tfz=CRY_INV[2][0]*a.x+CRY_INV[2][1]*a.y+CRY_INV[2][2]*a.z;
            tip+=' · frac ('+tfx.toFixed(3)+', '+tfy.toFixed(3)+', '+tfz.toFixed(3)+')';
            tip+=' · cell ('+a.cellI+','+a.cellJ+','+a.cellK+')';
            var tba=CRY.baseAtoms[a.baseIdx];
            if(tba&&tba.siteLabel){
                tip+=' · site '+tba.siteLabel;
                // G9: the operation that placed this image, by its position in
                // the file's own table. Only atoms expanded from an asymmetric
                // unit have it; a POSCAR's atoms have no operation to name.
                if(tba.symOp!=null&&CRY.symmetryOps&&CRY.symmetryOps[tba.symOp]!=null){
                    tip+=' (op #'+(tba.symOp+1)+': '+CRY.symmetryOps[tba.symOp]+')';
                }
            }
            if(a.occupancy!=null)tip+=' · occ '+a.occupancy;
            // H7: a bond-valence sum belongs to the base site, not to the
            // image under the pointer: every periodic image has the same
            // neighbours. It is quoted only when the site is a centre the
            // parameter table can evaluate, so other structures read as before.
            tip+=bvsTooltipSuffix(a);
        }
        tooltipEl.textContent=tip;
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
        // The atom handles are records, not scene children, so a box test takes
        // the same local position the instance matrices carry and applies the
        // group's world matrix itself — the number a Mesh.getWorldPosition()
        // used to return.
        if(diffMode){
            var dw=rect.width,dh=rect.height,dhalfW=Math.floor(dw/2);
            diffActiveSide='left';
            camera.aspect=dhalfW/dh;camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);camera.updateMatrixWorld(true);
            moleculeGroup.updateWorldMatrix(true,false);
            atomMeshes.forEach(function(m,i){
                var a=MD.atoms[i];if(!a)return;
                v.set(a.x-CX,a.y-CY,a.z-CZ).applyMatrix4(moleculeGroup.matrixWorld);
                v.project(camera);
                var sx=(v.x*0.5+0.5)*dhalfW+rect.left;
                var sy=(-v.y*0.5+0.5)*dh+rect.top;
                if(sx>=bx&&sx<=bx+bw&&sy>=by&&sy<=by+bh){selectAtom(m.userData.index);}
            });
            if(diffAtomMeshes&&diffAtomMeshes.length>0){
                diffActiveSide='right';
                camera.aspect=(dw-dhalfW)/dh;camera.updateProjectionMatrix();
                camera.position.set(0,0,diffCamDist);camera.updateMatrixWorld(true);
                diffMolGroup.updateWorldMatrix(true,false);
                diffAtomMeshes.forEach(function(m,i){
                    var a=diffData.atoms[i];if(!a)return;
                    v.set(a.x-diffCX,a.y-diffCY,a.z-diffCZ).applyMatrix4(diffMolGroup.matrixWorld);
                    v.project(camera);
                    var sx=(v.x*0.5+0.5)*(dw-dhalfW)+rect.left+dhalfW;
                    var sy=(-v.y*0.5+0.5)*dh+rect.top;
                    if(sx>=bx&&sx<=bx+bw&&sy>=by&&sy<=by+bh){selectAtom(m.userData.index);}
                });
            }
        }else{
            moleculeGroup.updateWorldMatrix(true,false);
            atomMeshes.forEach(function(m,i){
                var a=MD.atoms[i];if(!a)return;
                v.set(a.x-CX,a.y-CY,a.z-CZ).applyMatrix4(moleculeGroup.matrixWorld);
                v.project(camera);
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
    // The powder pattern is a reading aid with a close button of its own: Esc
    // backs out of it before it takes away a tool or a selection.
    if(xrdOpen){xrdClosePanel();return true}
    if(cavOpen){cavityClosePanel();return true}
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
