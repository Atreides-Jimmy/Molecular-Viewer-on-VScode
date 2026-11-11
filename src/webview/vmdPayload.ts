// VMD scene resolution: egg metadata and the drawable payload of one script
// molecule.
import * as vscode from 'vscode';
import * as path from 'path';
import { MolecularData, VmdMolecule, VmdScene } from '../types';
import {
    VmdEgg, parseUnitFile, parseEggFile, matchEggsToSegments, resolveVmdScene,
    resolveSelection, vdwRadius
} from '../parsers/index';
import {
    VmdAtomPayload, VmdBondPayload, VmdSegmentPayload, VmdRepPayload, VmdMoleculePayload
} from './payloads';
import { decodeText } from './loading';

/** Egg metadata read from the directory next to a structure file. */
export interface EggMeta { eggFile?: string; unitFile?: string }

/** What the eggFile chooser needs to score each candidate against the structure. */
export interface EggMatchContext {
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
export async function readEggMetadata(
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

export async function readEggMetadataFromDir(dirPath: string, matchContext?: EggMatchContext): Promise<EggMeta | undefined> {
    try {
        const dir = vscode.Uri.file(dirPath);
        const entries = await vscode.workspace.fs.readDirectory(dir);
        let unitName: string | undefined;
        const eggNames: string[] = [];
        for (const [name, type] of entries) {
            if ((type & vscode.FileType.File) === 0) continue;
            if (/^unitfile$/i.test(name)) unitName = name;
            else if (/^eggfile/i.test(name)) eggNames.push(name);
        }
        if (!unitName && eggNames.length === 0) return undefined;
        const readText = async (name: string): Promise<string | undefined> => {
            try {
                return decodeText(await vscode.workspace.fs.readFile(vscode.Uri.joinPath(dir, name)));
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
export function resolveVmdSource(scriptUri: vscode.Uri, sourceFile: string): vscode.Uri | undefined {
    if (!sourceFile) return undefined;
    if (path.isAbsolute(sourceFile)) return vscode.Uri.file(sourceFile);
    return vscode.Uri.file(path.resolve(path.dirname(scriptUri.fsPath), sourceFile));
}

/** Build the webview payload for one already-loaded VMD molecule. */
export function buildVmdMoleculePayload(
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
        warnings: resolved.warnings,
        // The structure's own state travels with the molecule: a script can
        // load structures of different charge, and switching to one must not
        // leave the previous molecule's state on screen or in a saved file.
        charge: data.charge,
        multiplicity: data.multiplicity,
        hasExplicitBonds: !!data.hasExplicitBonds,
        gjfMeta: data.gjfMeta,
        parseWarnings: data.warnings && data.warnings.length > 0 ? data.warnings.slice() : undefined
    };
}
