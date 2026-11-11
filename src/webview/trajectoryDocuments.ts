// Trajectory documents: pair a topology, flatten frames, and build the payload
// the page navigates and Diff/Import choose from.
import * as vscode from 'vscode';
import * as path from 'path';
import { Atom, Bond, MolecularData, TrajectoryInfo } from '../types';
import {
    ParseFileOptions, LammpsParseOptions, LammpsDumpTrajectory, TrajectoryFormat,
    parseFile, readBinaryTrajectory, lammpsLengthFactor, parseLammpsDumpTrajectory,
    parseXyzTrajectory, isExtXyz, isLammpsDump, parseLammpsInputHints
} from '../parsers/index';
import { ensureBonds } from '../tools/bondDetector';
import { TrajectoryFramePayload, TrajectoryDocument } from './payloads';
import { readDirectoryNames, isLammpsScriptName, readReferencedFile, scanLammpsChain } from './parseSetup';
import { decodeText } from './loading';

/**
 * The companion files of a trajectory.
 *
 * A binary trajectory states no elements, no bonds and — except for XTC/TRR's
 * nanometres — no unit, so everything it is drawn with comes from beside it:
 * the topology GROMACS wrote under the same name, or the LAMMPS script that
 * ran the dump and the data file it read. The scripts are filtered by whether
 * they name this trajectory, because a folder can hold several jobs whose
 * `units` and masses would otherwise rescale or relabel this one.
 */
export interface TrajectoryCompanions {
    units?: string;
    typeElements: { [type: number]: string };
    readData: { name: string; fromDir: vscode.Uri }[];
    /** Scripts that name this trajectory (0 when none does). */
    scripts: number;
    warnings: string[];
}

export interface TrajectoryTopologyChoice {
    data: MolecularData;
    name: string;
    path: string;
}

export interface TrajectoryTopologyLookup {
    topology?: TrajectoryTopologyChoice;
    warnings: string[];
}

export async function collectTrajectoryCompanions(uri: vscode.Uri, fileName: string): Promise<TrajectoryCompanions> {
    const result: TrajectoryCompanions = { typeElements: {}, readData: [], scripts: 0, warnings: [] };
    const dir = vscode.Uri.file(path.dirname(uri.fsPath));
    const index = await readDirectoryNames(dir);
    const scriptFiles: { uri: vscode.Uri; text: string }[] = [];
    // A folder may hold several LAMMPS scripts; reading them in a stable order
    // keeps the hints (units, type elements) the same from one open to the next.
    for (const name of Array.from(index.values()).sort((a, b) => a.localeCompare(b))) {
        if (name.toLowerCase() === fileName.toLowerCase()) continue;
        if (!isLammpsScriptName(name)) continue;
        try {
            const scriptUri = vscode.Uri.joinPath(dir, name);
            scriptFiles.push({
                uri: scriptUri,
                text: decodeText(await vscode.workspace.fs.readFile(scriptUri))
            });
        } catch {
            // An unreadable script contributes no hints.
        }
    }
    const escapedName = fileName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    const namesFile = new RegExp('(^|[^A-Za-z0-9_.-])' + escapedName + '($|[^A-Za-z0-9_.-])', 'i');
    let scripts = scriptFiles.filter(script => namesFile.test(script.text));
    // One script in the folder is that folder's job even when it builds the
    // dump name from a variable the name filter cannot see.
    if (scripts.length === 0 && scriptFiles.length === 1) scripts = scriptFiles;
    result.scripts = scripts.length;

    let units: string | undefined;
    for (const script of scripts) {
        const chain = await scanLammpsChain(script.uri, script.text, dir);
        result.warnings.push(...chain.warnings);
        result.readData.push(...chain.readData);
        for (const item of chain.texts) {
            const hints = parseLammpsInputHints(item.content, item.variables);
            if (!units && hints.units) units = hints.units;
            for (const key of Object.keys(hints.typeElements)) {
                const type = Number(key);
                if (!result.typeElements[type]) result.typeElements[type] = hints.typeElements[type];
            }
        }
    }
    if (units) result.units = units;
    return result;
}

/** A same-atom-count fallback candidate larger than this is not a topology. */
export const MAX_FALLBACK_TOPOLOGY_BYTES = 64 * 1024 * 1024;

/**
 * Choose the structure a trajectory's atom order belongs to.
 *
 * The same-basename files come first because that is the name GROMACS itself
 * uses (`water.xtc` beside `water.gro`), then the `read_data` file the LAMMPS
 * script names, and only when neither exists a coordinate file with the same
 * atom count. Every candidate is parsed through the ordinary readers, so a
 * `.data`, `.gro` or `.pdb` topology means exactly what opening it means.
 */
export async function chooseTrajectoryTopology(
    uri: vscode.Uri,
    fileName: string,
    natoms: number,
    companions: TrajectoryCompanions,
    coordinateTopology: boolean
): Promise<TrajectoryTopologyLookup> {
    const lookup: TrajectoryTopologyLookup = { warnings: [] };
    if (natoms <= 0) return lookup;
    const dir = vscode.Uri.file(path.dirname(uri.fsPath));
    const index = await readDirectoryNames(dir);
    const stem = fileName.replace(/\.(xtc|trr|dcd|lammpstrj|dump|crd|trj)$/i, '');
    const lammpsBase: LammpsParseOptions = {};
    if (companions.units) lammpsBase.units = companions.units;
    if (Object.keys(companions.typeElements).length > 0) lammpsBase.typeElements = companions.typeElements;

    const parseTopology = async (name: string): Promise<MolecularData | undefined> => {
        const actual = index.get(name.toLowerCase());
        if (!actual) return undefined;
        try {
            const text = decodeText(
                await vscode.workspace.fs.readFile(vscode.Uri.joinPath(dir, actual))
            );
            const isData = /\.data$/i.test(actual) || /^data\b/i.test(actual);
            const options: ParseFileOptions = {};
            const lammps: LammpsParseOptions = { ...lammpsBase };
            // Binary dumps (`dump dcd`, `dump xtc`) and this reader's text
            // dumps order atoms by ID; a data file whose Atoms section is not
            // sorted must be put in that order before its elements are paired.
            if (isData) lammps.orderById = true;
            if (Object.keys(lammps).length > 0) options.lammps = lammps;
            const parsed = parseFile(text, actual, options);
            if (parsed.atoms.length !== natoms) {
                lookup.warnings.push('the topology "' + actual + '" holds ' + parsed.atoms.length +
                    ' atoms, not ' + natoms + ': not paired');
                return undefined;
            }
            return parsed;
        } catch {
            return undefined;
        }
    };
    const accept = (data: MolecularData, name: string): TrajectoryTopologyLookup => {
        const actual = index.get(name.toLowerCase()) || name;
        lookup.topology = {
            data,
            name: actual,
            path: vscode.Uri.joinPath(dir, actual).fsPath
        };
        return lookup;
    };

    if (coordinateTopology) {
        for (const ext of ['.gro', '.pdb', '.ent', '.xyz', '.extxyz', '.data']) {
            const data = await parseTopology(stem + ext);
            if (data) return accept(data, stem + ext);
        }
    } else {
        // A LAMMPS dump's rows follow atom IDs. Only a data file states the
        // same order, so a coordinate topology is not paired with it.
        const data = await parseTopology(stem + '.data');
        if (data) return accept(data, stem + '.data');
    }

    for (const target of companions.readData) {
        const file = await readReferencedFile(dir, target.fromDir, target.name);
        if (!file) continue;
        try {
            const options: ParseFileOptions = {
                lammps: { ...lammpsBase, orderById: true }
            };
            const data = parseFile(file.content, file.fileName, options);
            if (data.atoms.length === natoms) return accept(data, file.fileName);
            lookup.warnings.push('the read_data file "' + file.fileName + '" holds ' + data.atoms.length +
                ' atoms, not ' + natoms + ': not paired');
        } catch {
            // A data file this viewer cannot read is not a topology candidate.
        }
    }

    if (!coordinateTopology || natoms <= 0) return lookup;
    // Last resort: exactly one coordinate file in the folder with the same
    // atom count. With several there is no way to tell which order the
    // trajectory used, so none is picked rather than guessing.
    const candidates: string[] = [];
    for (const name of Array.from(index.values()).sort((a, b) => a.localeCompare(b))) {
        if (name.toLowerCase() === fileName.toLowerCase()) continue;
        if (!/\.(gro|pdb|ent|xyz|extxyz)$/i.test(name) && !/\.data$/i.test(name) && !/^data\b/i.test(name)) continue;
        candidates.push(name);
        if (candidates.length > 40) break;
    }
    let match: { data: MolecularData; name: string } | undefined;
    let matchedCount = 0;
    for (const name of candidates) {
        const actual = index.get(name.toLowerCase());
        if (!actual) continue;
        try {
            const candidateUri = vscode.Uri.joinPath(dir, actual);
            // Parsing every same-sized file in a folder can mean reading MD
            // dumps that are larger than the trajectory itself; a coordinate
            // topology is a small file, so an oversized candidate is skipped.
            const info = await vscode.workspace.fs.stat(candidateUri);
            if (info.size > MAX_FALLBACK_TOPOLOGY_BYTES) continue;
            const text = decodeText(await vscode.workspace.fs.readFile(candidateUri));
            const isData = /\.data$/i.test(actual) || /^data\b/i.test(actual);
            const options: ParseFileOptions = {};
            if (isData) options.lammps = { ...lammpsBase, orderById: true };
            const parsed = parseFile(text, actual, options);
            if (parsed.atoms.length !== natoms) continue;
            matchedCount++;
            if (!match) match = { data: parsed, name: actual };
            if (matchedCount > 1) break;
        } catch {
            // Not a readable structure.
        }
    }
    if (matchedCount === 1 && match) {
        lookup.warnings.push('no same-name topology: paired with "' + match.name +
            '" (same atom count); if that file orders its atoms differently, the elements will not match');
        return accept(match.data, match.name);
    }
    if (matchedCount > 1) {
        lookup.warnings.push('several files in the folder hold ' + natoms +
            ' atoms: no topology was paired, so its elements and bonds are not known');
    }
    return lookup;
}

export function formatTrajectoryNumber(value: number): string {
    if (!isFinite(value)) return String(value);
    const rounded = Math.round(value);
    if (Math.abs(value - rounded) < 1e-6) return String(rounded);
    return value.toPrecision(6).replace(/\.?0+$/, '');
}

export function trajectoryStepLabel(
    format: TrajectoryFormat | 'dump',
    frame: { step?: number; time?: number; timestep?: string },
    index: number
): string {
    let label: string;
    if (format === 'dump') {
        label = frame.timestep ? 'Timestep ' + frame.timestep : 'Frame ' + (index + 1);
    } else {
        label = frame.step !== undefined ? 'Step ' + frame.step : 'Frame ' + (index + 1);
    }
    if (frame.time !== undefined && isFinite(frame.time)) {
        // DCD's time column is in the program's own unit — AKMA for
        // CHARMM/NAMD, ps for LAMMPS — so only XTC/TRR's ps is labelled.
        label += ' · t = ' + formatTrajectoryNumber(frame.time) +
            (format === 'xtc' || format === 'trr' ? ' ps' : '');
    }
    return label;
}

export function trajectoryTitle(format: TrajectoryFormat | 'dump'): string {
    if (format === 'xtc') return 'XTC trajectory';
    if (format === 'trr') return 'TRR trajectory';
    if (format === 'dcd') return 'DCD trajectory';
    return 'LAMMPS dump trajectory';
}

/**
 * A trajectory template: the atom list every frame is drawn with.
 *
 * The elements come from the topology when there is one. Without it a binary
 * trajectory has no element information at all; its atoms are shown as `Xx`,
 * the conventional unknown species, rather than `X`, which this viewer treats
 * as a ghost atom and would hide.
 */
export function trajectoryTemplate(
    natoms: number,
    elements: string[] | undefined,
    topology: TrajectoryTopologyChoice | undefined
): { atoms: Atom[]; bonds: Bond[]; hasExplicitBonds: boolean; warnings: string[] } {
    const warnings: string[] = [];
    const atoms: Atom[] = [];
    for (let i = 0; i < natoms; i++) {
        const element = elements && elements[i] ? elements[i] : 'Xx';
        const source = topology ? topology.data.atoms[i] : undefined;
        atoms.push({
            element,
            x: 0, y: 0, z: 0,
            index: i,
            occupancy: source ? source.occupancy : undefined,
            siteLabel: source ? source.siteLabel : undefined
        });
    }
    if (!topology) {
        return { atoms, bonds: [], hasExplicitBonds: false, warnings };
    }
    const bonded = ensureBonds({
        atoms: topology.data.atoms,
        bonds: topology.data.bonds,
        title: topology.data.title,
        hasExplicitBonds: topology.data.hasExplicitBonds
    });
    return {
        atoms,
        bonds: bonded.bonds,
        hasExplicitBonds: bonded.hasExplicitBonds,
        warnings
    };
}

export function finishTrajectoryDocument(
    title: string,
    template: { atoms: Atom[]; bonds: Bond[]; hasExplicitBonds: boolean },
    trajectoryFrames: TrajectoryFramePayload[],
    metadata: TrajectoryInfo,
    warnings: string[]
): TrajectoryDocument {
    if (trajectoryFrames.length > 0) {
        const first = trajectoryFrames[0];
        for (let i = 0; i < template.atoms.length; i++) {
            template.atoms[i].x = first.xyz[i * 3];
            template.atoms[i].y = first.xyz[i * 3 + 1];
            template.atoms[i].z = first.xyz[i * 3 + 2];
        }
    }
    const data: MolecularData = {
        atoms: template.atoms,
        bonds: template.bonds,
        title,
        hasExplicitBonds: template.hasExplicitBonds,
        trajectory: metadata,
        warnings: Array.from(new Set(warnings))
    };
    if (data.warnings && data.warnings.length === 0) delete data.warnings;
    return { data, trajectoryFrames };
}

export function emptyTrajectoryDocument(title: string, warnings: string[]): TrajectoryDocument {
    return {
        data: {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: Array.from(new Set(warnings))
        },
        trajectoryFrames: []
    };
}

export async function loadBinaryTrajectoryDocument(
    uri: vscode.Uri,
    fileName: string,
    bytes: Uint8Array,
    format: TrajectoryFormat
): Promise<TrajectoryDocument> {
    const title = trajectoryTitle(format);
    const companions = await collectTrajectoryCompanions(uri, fileName);
    const warnings = companions.warnings.slice();
    if (companions.units && companions.units.toLowerCase() === 'lj') {
        warnings.push('the LAMMPS run uses units lj: DCD coordinates are in reduced units and are shown unscaled');
    }
    const dcdUnits = companions.units ? lammpsLengthFactor(companions.units) : 1;
    if (format === 'dcd' && companions.units && dcdUnits === undefined) {
        warnings.push('the LAMMPS units "' + companions.units +
            '" are not a length unit this reader knows: DCD coordinates are shown unscaled');
    }
    const lengthScale = format === 'dcd' ? (dcdUnits === undefined ? 1 : dcdUnits) : 10;
    const trajectory = readBinaryTrajectory(bytes, { lengthScale });
    warnings.push(...trajectory.warnings);
    if (trajectory.sampled) {
        warnings.push(trajectory.frameCount + ' frames in the file: showing ' +
            trajectory.frames.length + ' sampled evenly');
    }
    if (trajectory.frames.length === 0) {
        warnings.push('the trajectory holds no frame this reader could decode');
        return emptyTrajectoryDocument(title, warnings);
    }
    let topology: TrajectoryTopologyChoice | undefined;
    if (trajectory.natoms > 0) {
        const lookup = await chooseTrajectoryTopology(uri, fileName, trajectory.natoms, companions, true);
        warnings.push(...lookup.warnings);
        topology = lookup.topology;
    }
    if (!topology) {
        warnings.push('no topology was found for this trajectory: its elements are unknown (shown as Xx) ' +
            'and no bonds are drawn; open or keep the matching .gro/.pdb/.xyz or LAMMPS data file beside it');
    }
    const template = trajectoryTemplate(
        trajectory.natoms,
        topology ? topology.data.atoms.map(a => a.element) : undefined,
        topology
    );
    warnings.push(...template.warnings);
    const frames: TrajectoryFramePayload[] = [];
    for (let k = 0; k < trajectory.frames.length; k++) {
        const frame = trajectory.frames[k];
        if (frame.positions.length !== trajectory.natoms * 3) continue;
        frames.push({
            stepLabel: trajectoryStepLabel(format, frame, k),
            xyz: Array.from(frame.positions),
            cell: frame.cell,
            time: frame.time,
            step: frame.step
        });
    }
    return finishTrajectoryDocument(
        topology ? title + ' — ' + topology.name : title,
        template,
        frames,
        {
            format,
            source: fileName,
            topology: topology ? topology.name : undefined,
            framesInFile: trajectory.frameCount,
            sampled: trajectory.sampled,
            truncated: trajectory.truncated
        },
        warnings
    );
}

export async function loadDumpTrajectoryDocument(
    uri: vscode.Uri,
    fileName: string,
    trajectory: LammpsDumpTrajectory,
    options: LammpsParseOptions | undefined
): Promise<TrajectoryDocument> {
    const title = trajectoryTitle('dump');
    const companions = await collectTrajectoryCompanions(uri, fileName);
    const warnings = companions.warnings.slice();
    if (companions.units && companions.units.toLowerCase() === 'lj') {
        warnings.push('the LAMMPS run uses units lj: dump coordinates are in reduced units and are shown unscaled');
    }
    warnings.push(...trajectory.warnings);
    if (trajectory.frames.length === 0) {
        return emptyTrajectoryDocument(title, warnings);
    }
    if (trajectory.sampled) {
        warnings.push(trajectory.frameCount + ' frames in the file: showing ' +
            trajectory.frames.length + ' sampled evenly');
    }
    const lookup = await chooseTrajectoryTopology(uri, fileName, trajectory.natoms, companions, false);
    warnings.push(...lookup.warnings);
    const topology = lookup.topology;
    const first = trajectory.frames[0];
    const dumpElements = first.sites.map(s => s.element).filter(e => !/^T\d+$/.test(e));
    const elements = topology
        ? topology.data.atoms.map(a => a.element)
        : (dumpElements.length === first.sites.length ? dumpElements : undefined);
    const template = trajectoryTemplate(trajectory.natoms, elements, topology);
    if (!topology && !elements) {
        warnings.push('the dump names its atoms by LAMMPS type and no data file or element list was found: ' +
            'the types are shown as labels and no bonds are drawn');
    } else if (!topology) {
        warnings.push('no topology was found for this dump: its element column names the atoms but no bonds are drawn; ' +
            'keep the matching .data file beside it to draw bonds');
    }
    const frames: TrajectoryFramePayload[] = [];
    for (let k = 0; k < trajectory.frames.length; k++) {
        const frame = trajectory.frames[k];
        const xyz: number[] = [];
        for (const site of frame.sites) {
            const cart = site.cart || [0, 0, 0];
            xyz.push(cart[0], cart[1], cart[2]);
        }
        if (xyz.length !== trajectory.natoms * 3) continue;
        frames.push({
            stepLabel: trajectoryStepLabel('dump', { step: frame.step, timestep: frame.timestep }, k),
            xyz,
            cell: frame.lattice,
            step: frame.step
        });
    }
    return finishTrajectoryDocument(
        topology ? title + ' — ' + topology.name : title,
        template,
        frames,
        {
            format: 'lammps-dump',
            source: fileName,
            topology: topology ? topology.name : undefined,
            framesInFile: trajectory.frameCount,
            sampled: trajectory.sampled,
            truncated: trajectory.truncated
        },
        warnings
    );
}

export interface PickedTextFrame {
    atoms: Atom[];
    title: string;
    label: string;
    /** The frame's cell, when the file states one (a LAMMPS dump). */
    lattice?: number[][];
}

/**
 * Every frame a text trajectory states, for Diff or Import to choose from.
 *
 * Opening such a file gives frame navigation, so silently using its first
 * frame would compare or import a geometry the user never chose. A file that
 * is not a multi-frame trajectory yields an empty list and the caller keeps
 * its ordinary single-structure path.
 */
export function parseTextTrajectoryFrames(
    content: string,
    fileName: string,
    ext: string,
    options: ParseFileOptions | undefined
): PickedTextFrame[] {
    const frames: PickedTextFrame[] = [];
    if ((ext === 'xyz' || ext === 'trj') && !isExtXyz(content)) {
        const trajectory = parseXyzTrajectory(content);
        if (trajectory.frames.length < 2) return [];
        for (let k = 0; k < trajectory.frames.length; k++) {
            const frame = trajectory.frames[k];
            frames.push({
                atoms: frame.atoms,
                title: frame.title,
                label: 'Frame ' + (k + 1) + (frame.title ? ': ' + frame.title : '')
            });
        }
    } else if (isLammpsDump(content)) {
        const trajectory = parseLammpsDumpTrajectory(content, options && options.lammps);
        if (trajectory.frames.length < 2) return [];
        for (const frame of trajectory.frames) {
            const atoms: Atom[] = frame.sites.map((site, i) => {
                const cart = site.cart || [0, 0, 0];
                return { element: site.element, x: cart[0], y: cart[1], z: cart[2], index: i };
            });
            frames.push({
                atoms,
                title: fileName,
                label: 'Frame ' + (frame.index + 1) + ' · Timestep ' + frame.timestep,
                lattice: frame.lattice
            });
        }
    }
    return frames;
}

/**
 * Ask which frame of a text trajectory to use.
 *
 * The picker is deliberately outside the parse's progress notification: a
 * spinner that keeps running while the user decides is not a status, it is a
 * stuck notification. 'cancelled' means the picker was dismissed.
 */
export async function pickTextTrajectoryFrame(
    frames: PickedTextFrame[],
    action: 'compare' | 'import'
): Promise<PickedTextFrame | 'cancelled'> {
    const picked = await vscode.window.showQuickPick(
        frames.map((frame, i) => ({
            label: frame.label,
            description: frame.atoms.length + ' atoms',
            index: i
        })),
        {
            placeHolder: 'Select frame to ' + action + ' (' + frames.length + ' frames found)',
            canPickMany: false
        }
    );
    if (!picked) return 'cancelled';
    return frames[picked.index];
}

/** Cell edge lengths and angles of a lattice, for the diff cell panel. */
export function latticeCellParameters(lattice: number[][]) {
    const length = (v: number[]) => Math.hypot(v[0], v[1], v[2]);
    const angle = (u: number[], v: number[]) => {
        const denom = length(u) * length(v);
        if (!(denom > 0)) return 90;
        const cos = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / denom;
        return Math.acos(Math.max(-1, Math.min(1, cos))) * 180 / Math.PI;
    };
    return {
        a: length(lattice[0]), b: length(lattice[1]), c: length(lattice[2]),
        alpha: angle(lattice[1], lattice[2]),
        beta: angle(lattice[0], lattice[2]),
        gamma: angle(lattice[0], lattice[1])
    };
}
