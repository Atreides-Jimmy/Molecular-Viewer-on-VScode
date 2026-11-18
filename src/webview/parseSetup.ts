// Parse preparation: resolve the files a main input names, expand the include
// chains of LAMMPS/CP2K/ABINIT/ABACUS/ORCA, and build ParseFileOptions.
import * as vscode from 'vscode';
import * as path from 'path';
import { MolecularData } from '../types';
import {
    ParseFileOptions, DeepmdParseOptions, LammpsParseOptions, LammpsHints, TrajectoryFormat,
    detectBinaryTrajectory, looksLikeBinaryTrajectory,
    isLammpsInput, isLammpsData, isLammpsDump,
    isCp2kInput, isAbacusInput, isAbinitInput, scanAbacusInput, scanAbinitXyzfile,
    scanOrcaXyzfile, scanLammpsInput, parseLammpsInputHints, lammpsTypeElements,
    deepmdFirstBox, deepmdBoxFrames, parseDeepmdTypes, parseDeepmdTypeMap
} from '../parsers/index';
import { decodeText, readFileWithProgress } from './loading';

/** A file Diff or Import reads: its bytes, its text (unless binary) and its companion options. */
export interface ToolFileLoad {
    content: Uint8Array;
    binary: TrajectoryFormat | null;
    /** True when the bytes are binary but not a trajectory format this reader knows. */
    binaryUnsupported?: boolean;
    text: string;
    setup: ParseSetup | undefined;
}

/**
 * Read, decode and prepare one file for Diff or Import.
 *
 * The read is the only part that can be reported byte by byte; decoding and
 * `buildParseOptions` follow inside the same notification, so the bar does not
 * close and leave a silent gap before the parser's own notification starts.
 */
export async function readToolFile(
    uri: vscode.Uri,
    fileName: string,
    size: number,
    report?: (fraction: number, message: string) => void,
    isCancelled?: () => boolean
): Promise<ToolFileLoad> {
    const content = report
        ? await readFileWithProgress(uri, size, (fraction, message) => report(fraction * 0.8, message),
            isCancelled || (() => false))
        : await vscode.workspace.fs.readFile(uri);
    if (report) report(0.85, 'Decoding the text');
    const binary = detectBinaryTrajectory(content);
    // `.trj` is a text trajectory for xTB and a binary one for CHARMM/NAMD.
    // A binary file that no magic recognised must not be decoded into mojibake
    // and reported as a parse failure; the caller can name the real reason.
    const binaryUnsupported = !binary &&
        (/\.(xtc|trr|dcd|trj)$/i.test(fileName)) &&
        looksLikeBinaryTrajectory(content);
    const text = binary || binaryUnsupported ? '' : decodeText(content);
    if (report) report(0.9, 'Preparing the parse');
    const setup = binary || binaryUnsupported ? undefined : await buildParseOptions(uri, fileName, content, text);
    if (report) report(1, 'Prepared');
    return { content, binary, binaryUnsupported, text, setup };
}

/** A structure file a main input names, already read from disk. */
export interface StructureReference {
    uri: vscode.Uri;
    content: string;
    fileName: string;
}

/** Everything reading one file's companions produced. */
export interface ParseSetup {
    options: ParseFileOptions;
    /** The structure a main input names (`read_data`, `xyzfile`, `stru_file`, `COORD_FILE_NAME`). */
    reference?: StructureReference;
    /** The input names a structure file, even if that file could not be read. */
    referenceExpected?: boolean;
    /** The text as parsed: include files spliced in when the format has directives. */
    text?: string;
    warnings: string[];
}

/** The file itself, or the same name in another case when the folder allows it. */
export async function fileIfExists(uri: vscode.Uri): Promise<vscode.Uri | undefined> {
    try {
        const stat = await vscode.workspace.fs.stat(uri);
        if ((stat.type & vscode.FileType.File) !== 0) return uri;
        return undefined;
    } catch {
        // The folder may be case-sensitive; try the name case-insensitively.
    }
    try {
        const folder = vscode.Uri.file(path.dirname(uri.fsPath));
        const wanted = path.basename(uri.fsPath).toLowerCase();
        for (const [name, type] of await vscode.workspace.fs.readDirectory(folder)) {
            if ((type & vscode.FileType.File) !== 0 && name.toLowerCase() === wanted) {
                return vscode.Uri.joinPath(folder, name);
            }
        }
    } catch {
        // An unreadable folder has no file to offer.
    }
    return undefined;
}

/**
 * Resolve a path a main input states.
 *
 * The job's working directory is usually the folder of the main input, so a
 * relative path is tried there first and then beside the file that named it
 * (an `include`). A path written for another launch directory falls back to
 * the bare file name beside those two places, which is how the reader has
 * always found a companion; an absolute path is used as written. A value that
 * still holds `${var}` or `@var` is reported as unresolved instead of being
 * searched for literally.
 */
export async function resolveReferencedUri(
    mainDir: vscode.Uri,
    fromDir: vscode.Uri,
    wanted: string
): Promise<vscode.Uri | undefined> {
    const cleaned = wanted.trim().replace(/^['"]|['"]$/g, '');
    if (cleaned === '') return undefined;
    // Anything still holding a LAMMPS/CP2K variable cannot be opened as
    // written; searching for the literal text would name the wrong file.
    if (/\$/.test(cleaned)) return undefined;
    if (/^@[A-Za-z_]/.test(cleaned)) return undefined;
    const normalized = cleaned.replace(/\\/g, '/');
    const bases: string[] = [];
    for (const folder of [mainDir.fsPath, fromDir.fsPath]) {
        if (!bases.includes(folder)) bases.push(folder);
    }
    if (path.isAbsolute(normalized) || /^[A-Za-z]:[\\/]/.test(cleaned)) {
        const direct = await fileIfExists(vscode.Uri.file(normalized));
        if (direct) return direct;
        // A path written on another machine (`/home/u/data.lmp`) is not on
        // this disk, but the same file is often beside the job that names it:
        // fall back to its bare name in the two folders searched below.
        const absoluteBare = path.basename(normalized);
        if (absoluteBare === '') return undefined;
        for (const folder of bases) {
            const hit = await fileIfExists(vscode.Uri.file(path.join(folder, absoluteBare)));
            if (hit) return hit;
        }
        return undefined;
    }
    const names = [normalized];
    const bare = path.basename(normalized);
    if (bare !== normalized) names.push(bare);
    for (const name of names) {
        for (const folder of bases) {
            const hit = await fileIfExists(vscode.Uri.file(path.join(folder, name)));
            if (hit) return hit;
        }
    }
    return undefined;
}

export async function readReferencedFile(
    mainDir: vscode.Uri,
    fromDir: vscode.Uri,
    wanted: string
): Promise<StructureReference | undefined> {
    const uri = await resolveReferencedUri(mainDir, fromDir, wanted);
    if (!uri) return undefined;
    try {
        return {
            uri,
            content: decodeText(await vscode.workspace.fs.readFile(uri)),
            fileName: path.basename(uri.fsPath)
        };
    } catch {
        return undefined;
    }
}

/**
 * Every file an `include` names.
 *
 * LAMMPS expands `*` and `?` in an include name, so `include coords.*` reads
 * every matching file in alphabetical order. A name without wildcards is the
 * ordinary single-file lookup.
 */
export async function resolveIncludeUris(
    mainDir: vscode.Uri,
    fromDir: vscode.Uri,
    wanted: string
): Promise<vscode.Uri[]> {
    const cleaned = wanted.trim().replace(/^['"]|['"]$/g, '');
    if (cleaned === '' || /\$/.test(cleaned) || /^@[A-Za-z_]/.test(cleaned)) return [];
    const normalized = cleaned.replace(/\\/g, '/');
    if (!/[*?]/.test(normalized)) {
        const hit = await resolveReferencedUri(mainDir, fromDir, cleaned);
        return hit ? [hit] : [];
    }
    const folderName = path.dirname(normalized);
    const basePattern = path.basename(normalized);
    const escape = (ch: string): string => ch.replace(/[.+^${}()|[\]\\]/g, '\\$&');
    const matcher = new RegExp('^' + basePattern.split('').map(ch =>
        ch === '*' ? '.*' : ch === '?' ? '.' : escape(ch)).join('') + '$',
        process.platform === 'win32' ? 'i' : '');
    const folders: vscode.Uri[] = [];
    if (path.isAbsolute(normalized) || /^[A-Za-z]:[\\/]/.test(cleaned)) {
        folders.push(vscode.Uri.file(folderName));
    } else {
        for (const base of [mainDir, fromDir]) {
            const folder = vscode.Uri.file(path.join(base.fsPath, folderName));
            if (!folders.some(f => f.fsPath === folder.fsPath)) folders.push(folder);
        }
    }
    const hits: vscode.Uri[] = [];
    for (const folder of folders) {
        try {
            const names = (await vscode.workspace.fs.readDirectory(folder))
                .filter(([, type]) => (type & vscode.FileType.File) !== 0)
                .map(([name]) => name)
                .filter(name => matcher.test(name))
                .sort();
            for (const name of names) hits.push(vscode.Uri.joinPath(folder, name));
        } catch {
            // An unreadable folder contributes no include files.
        }
    }
    return hits;
}

export interface LammpsInputChain {
    /** Each file with the variables visible when it started. */
    texts: { content: string; variables: { [name: string]: string } }[];
    readData: { name: string; fromDir: vscode.Uri }[];
    readRestart?: string;
    warnings: string[];
}

/**
 * Follow a LAMMPS input and its `include` files in command order.
 *
 * A script can split `units`, `mass` and `read_data` across included files, so
 * the chain contributes both the hints and the structure reference. Includes
 * are read relative to the main input's folder first and the including file's
 * folder second; a cycle, an unreadable include and an over-deep chain are
 * warnings, never a hang.
 */
export async function scanLammpsChain(
    mainUri: vscode.Uri,
    mainText: string,
    mainDir: vscode.Uri
): Promise<LammpsInputChain> {
    const chain: LammpsInputChain = { texts: [], readData: [], warnings: [] };
    // LAMMPS executes an `include` every time it appears: the same file may be
    // included twice to define different regions or to re-use a variable block,
    // and the second execution is part of the script's meaning. Only a file on
    // the *current* include stack is a cycle; a later, non-recursive include of
    // the same file is read and scanned again.
    const visiting = new Set<string>();
    const visit = async (
        fileUri: vscode.Uri,
        content: string,
        depth: number,
        inherited: { [name: string]: string }
    ): Promise<{ [name: string]: string }> => {
        const key = fileUri.fsPath.toLowerCase();
        if (visiting.has(key)) {
            chain.warnings.push('the LAMMPS include "' + path.basename(fileUri.fsPath) +
                '" forms a cycle: stopped there');
            return inherited;
        }
        if (depth > 8) {
            chain.warnings.push('the LAMMPS include chain is deeper than 8 files: stopped at ' +
                path.basename(fileUri.fsPath));
            return inherited;
        }
        visiting.add(key);
        try {
            const fromDir = vscode.Uri.file(path.dirname(fileUri.fsPath));
            const scan = scanLammpsInput(content, inherited);
            chain.texts.push({ content, variables: { ...inherited } });
            for (const name of scan.readData) chain.readData.push({ name, fromDir });
            if (scan.readRestart) chain.readRestart = scan.readRestart;
            let variables = { ...inherited, ...scan.variables };
            for (const include of scan.includes) {
                const uris = await resolveIncludeUris(mainDir, fromDir, include);
                if (uris.length === 0) {
                    chain.warnings.push('the LAMMPS include "' + include + '" could not be found');
                    continue;
                }
                for (const uri of uris) {
                    let text: string;
                    try {
                        text = decodeText(await vscode.workspace.fs.readFile(uri));
                    } catch {
                        chain.warnings.push('the LAMMPS include "' + include + '" could not be read');
                        continue;
                    }
                    variables = await visit(uri, text, depth + 1, variables);
                }
            }
            return variables;
        } finally {
            visiting.delete(key);
        }
    };
    await visit(mainUri, mainText, 0, {});
    return chain;
}

/** A directive value without surrounding quotes or a trailing comment. */
export function directiveValue(rest: string): string {
    const value = rest.trim();
    const quoted = /^(['"])(.*?)\1/.exec(value);
    if (quoted) return quoted[2].trim();
    const cut = value.search(/[#!]/);
    return (cut >= 0 ? value.slice(0, cut) : value).trim();
}

/** Substitute CP2K's `$VAR` and `${VAR}` references; unknown ones stay as written. */
export function cp2kSubstitute(text: string, variables: { [name: string]: string }): string {
    return text.replace(/\$\{([A-Za-z_][A-Za-z0-9_]*)\}|\$([A-Za-z_][A-Za-z0-9_]*)/g,
        (whole: string, braced: string, bare: string) => {
            const name = (braced || bare).toLowerCase();
            return variables[name] !== undefined ? variables[name] : whole;
        });
}

/** One `@IF` expression: `==`, `/=`, or a bare value's truthiness. */
export function cp2kConditionTrue(expression: string, variables: { [name: string]: string }): boolean {
    let body = expression.trim();
    if (body.startsWith('(') && body.endsWith(')')) body = body.slice(1, -1).trim();
    const compare = /^(.+?)\s*(==|\/=|!=)\s*(.+?)$/.exec(body);
    if (!compare) {
        const value = cp2kSubstitute(body, variables).trim().replace(/^['"]|['"]$/g, '').toLowerCase();
        return value !== '' && value !== '0' && value !== 'false' && value !== '.false.';
    }
    const left = cp2kSubstitute(compare[1], variables).trim().replace(/^['"]|['"]$/g, '');
    const right = cp2kSubstitute(compare[3], variables).trim().replace(/^['"]|['"]$/g, '');
    const equal = left === right;
    return compare[2] === '==' ? equal : !equal;
}

export interface IncludeExpansion {
    /** The file's text with its include files spliced in. */
    text: string;
    /** Every directory the text came from, the main file's first. */
    dirs: vscode.Uri[];
    warnings: string[];
}

export function addIncludeDir(dirs: vscode.Uri[], uri: vscode.Uri): void {
    const folder = vscode.Uri.file(path.dirname(uri.fsPath));
    if (!dirs.some(d => d.fsPath === folder.fsPath)) dirs.push(folder);
}

/**
 * Expand CP2K's preprocessor enough to find the structure the input names:
 * `@INCLUDE` splices a file in, `@SET`/`$VAR` resolve a path written through a
 * variable, and `@IF`/`@ELSE`/`@ENDIF` hide the branches CP2K would not read.
 * Cycles, unreadable files and over-deep chains warn instead of hanging.
 */
export async function expandCp2kInput(
    mainUri: vscode.Uri,
    mainText: string,
    mainDir: vscode.Uri
): Promise<IncludeExpansion> {
    const warnings: string[] = [];
    const dirs: vscode.Uri[] = [];
    const variables: { [name: string]: string } = {};
    // The files on the current include stack: a repeat is a cycle, while the
    // same file included twice along different branches is ordinary CP2K.
    const visiting = new Set<string>();
    addIncludeDir(dirs, mainUri);

    const expand = async (fileUri: vscode.Uri, content: string, depth: number): Promise<string> => {
        const key = fileUri.fsPath.toLowerCase();
        if (visiting.has(key)) {
            warnings.push('the CP2K include "' + path.basename(fileUri.fsPath) +
                '" forms a cycle: stopped there');
            return '';
        }
        visiting.add(key);
        if (depth > 8) {
            visiting.delete(key);
            warnings.push('the CP2K include chain is deeper than 8 files: stopped at ' +
                path.basename(fileUri.fsPath));
            return content;
        }
        const fromDir = vscode.Uri.file(path.dirname(fileUri.fsPath));
        const out: string[] = [];
        const active: boolean[] = [];
        const isActive = () => active.every(v => v);

        for (const raw of content.split(/\r?\n/)) {
            const trimmed = raw.trim();
            const includeMatch = /^@INCLUDE\s*(?:=\s*)?(.+)$/i.exec(trimmed);
            if (includeMatch) {
                if (!isActive()) { out.push(''); continue; }
                const wantedRaw = directiveValue(includeMatch[1]);
                const wanted = cp2kSubstitute(wantedRaw, variables);
                const target = wanted === '' ? undefined : await resolveReferencedUri(mainDir, fromDir, wanted);
                if (!target) {
                    warnings.push('the CP2K include "' + wantedRaw + '" could not be found');
                    out.push('');
                    continue;
                }
                let includeText: string;
                try {
                    includeText = decodeText(await vscode.workspace.fs.readFile(target));
                } catch {
                    warnings.push('the CP2K include "' + wantedRaw + '" could not be read');
                    out.push('');
                    continue;
                }
                addIncludeDir(dirs, target);
                out.push(await expand(target, includeText, depth + 1));
                continue;
            }
            const setMatch = /^@SET\s+([A-Za-z_][A-Za-z0-9_]*)\s*(?:=|:)?\s*(.+)$/i.exec(trimmed);
            if (setMatch) {
                if (isActive()) {
                    variables[setMatch[1].toLowerCase()] = cp2kSubstitute(directiveValue(setMatch[2]), variables);
                }
                out.push('');
                continue;
            }
            const ifdefMatch = /^@IFDEF\s+([A-Za-z_][A-Za-z0-9_]*)/i.exec(trimmed);
            if (ifdefMatch) {
                active.push(isActive() && variables[ifdefMatch[1].toLowerCase()] !== undefined);
                out.push('');
                continue;
            }
            const ifndefMatch = /^@IFNDEF\s+([A-Za-z_][A-Za-z0-9_]*)/i.exec(trimmed);
            if (ifndefMatch) {
                active.push(isActive() && variables[ifndefMatch[1].toLowerCase()] === undefined);
                out.push('');
                continue;
            }
            const ifMatch = /^@IF\s+(.+)$/i.exec(trimmed);
            if (ifMatch) {
                active.push(isActive() && cp2kConditionTrue(ifMatch[1], variables));
                out.push('');
                continue;
            }
            if (/^@ELSE\b/i.test(trimmed)) {
                if (active.length > 0) active[active.length - 1] = !active[active.length - 1];
                out.push('');
                continue;
            }
            if (/^@ENDIF\b/i.test(trimmed)) {
                if (active.length > 0) active.pop();
                out.push('');
                continue;
            }
            out.push(isActive() ? cp2kSubstitute(raw, variables) : '');
        }
        visiting.delete(key);
        return out.join('\n');
    };

    const text = await expand(mainUri, mainText, 0);
    return { text, dirs, warnings };
}

/** ABINIT's include line: `include "file"` (quotes optional). */
export function abinitIncludeName(line: string): string | undefined {
    const match = /^\s*include\s*(?:=\s*)?(.+)$/i.exec(line);
    if (!match) return undefined;
    return directiveValue(match[1]);
}

/** Splice ABINIT `include` files into one text, with the same cycle guards. */
export async function expandAbinitInput(
    mainUri: vscode.Uri,
    mainText: string,
    mainDir: vscode.Uri
): Promise<IncludeExpansion> {
    const warnings: string[] = [];
    const dirs: vscode.Uri[] = [];
    // A file repeated on the current stack is a cycle; the same file included
    // again along another branch is simply read again, as ABINIT would.
    const visiting = new Set<string>();
    addIncludeDir(dirs, mainUri);

    const expand = async (fileUri: vscode.Uri, content: string, depth: number): Promise<string> => {
        const key = fileUri.fsPath.toLowerCase();
        if (visiting.has(key)) {
            warnings.push('the ABINIT include "' + path.basename(fileUri.fsPath) +
                '" forms a cycle: stopped there');
            return '';
        }
        visiting.add(key);
        if (depth > 8) {
            visiting.delete(key);
            warnings.push('the ABINIT include chain is deeper than 8 files: stopped at ' +
                path.basename(fileUri.fsPath));
            return content;
        }
        const fromDir = vscode.Uri.file(path.dirname(fileUri.fsPath));
        const out: string[] = [];
        for (const line of content.split(/\r?\n/)) {
            const wanted = abinitIncludeName(line);
            if (wanted === undefined) { out.push(line); continue; }
            if (wanted === '') { out.push(''); continue; }
            const target = await resolveReferencedUri(mainDir, fromDir, wanted);
            if (!target) {
                warnings.push('the ABINIT include "' + wanted + '" could not be found');
                out.push('');
                continue;
            }
            let includeText: string;
            try {
                includeText = decodeText(await vscode.workspace.fs.readFile(target));
            } catch {
                warnings.push('the ABINIT include "' + wanted + '" could not be read');
                out.push('');
                continue;
            }
            addIncludeDir(dirs, target);
            out.push(await expand(target, includeText, depth + 1));
        }
        visiting.delete(key);
        return out.join('\n');
    };

    const text = await expand(mainUri, mainText, 0);
    return { text, dirs, warnings };
}

/** Read a named companion from the main folder first, then from each include's folder. */
export async function readReferencedFromDirs(
    mainDir: vscode.Uri,
    dirs: vscode.Uri[],
    wanted: string
): Promise<StructureReference | undefined> {
    for (const fromDir of dirs) {
        const file = await readReferencedFile(mainDir, fromDir, wanted);
        if (file) return file;
    }
    return undefined;
}

/** Join host warnings into the parsed structure without duplicating a line. */
export function applySetupWarnings(data: MolecularData, setup: ParseSetup | undefined): void {
    if (!setup) return;
    let warnings = data.warnings ? data.warnings.slice() : [];
    // The parser only sees the opened script, so an include-only read_data
    // looks like "no read_data"; the host followed the chain and knows better.
    if (setup.referenceExpected && !setup.reference) {
        warnings = warnings.filter(w => !/builds its structure at run time/.test(w));
    }
    for (const warning of setup.warnings) {
        if (!warnings.includes(warning)) warnings.push(warning);
    }
    if (warnings.length > 0) data.warnings = warnings;
}

export async function readDirectoryNames(dir: vscode.Uri): Promise<Map<string, string>> {
    const index = new Map<string, string>();
    try {
        for (const [name, type] of await vscode.workspace.fs.readDirectory(dir)) {
            if ((type & vscode.FileType.File) !== 0) index.set(name.toLowerCase(), name);
        }
    } catch {
        // An unreadable folder simply has no companions.
    }
    return index;
}

/**
 * A file name that can hold a LAMMPS job script.
 *
 * `in.*` and `*.in` are the conventional names, but the same commands are
 * routinely saved as `*.lmp`/`*.lammps` — the trajectory side has always
 * accepted those, and the data/dump side must read the same scripts or a job
 * named `run.lmp` loses its `units`, `atom_style` and type labels.
 */
export function isLammpsScriptName(name: string): boolean {
    return /^in\./i.test(name) || /\.in$/i.test(name) || /\.lmp$/i.test(name) || /\.lammps$/i.test(name);
}

/**
 * The companion files a structure names, read from beside it.
 *
 * A DeepMD dataset is a directory (`type.raw`, `type_map.raw`, `box.raw`, a
 * `nopbc` marker), a LAMMPS dump takes its types and units from the input
 * script that wrote it, and a CP2K input can keep its coordinates in a separate
 * file. The parsers themselves never touch the file system — they receive what
 * the host read — so this is the one place that decides what "beside it" means
 * and how far up it looks.
 */
export async function buildParseOptions(
    uri: vscode.Uri,
    fileName: string,
    bytes: Uint8Array,
    text: string
): Promise<ParseSetup | undefined> {
    const dir = vscode.Uri.file(path.dirname(uri.fsPath));
    const parent = vscode.Uri.file(path.dirname(dir.fsPath));
    const lowerName = fileName.toLowerCase();
    const coordArray = /^coord\.(raw|npy)$/i.test(fileName);
    // A program output or a plain coordinate file names no companion input of
    // the kinds below. Testing every family against a 60 MB text cost seconds
    // before the parser started; the names that can actually carry companions
    // are the LAMMPS inputs/data/dumps, the DeepMD arrays and the program
    // inputs, all excluded here. A `.xyz` that is really a LAMMPS dump is kept
    // by its own name or its first header line.
    const lammpsNamed = /^(dump|data)\b/i.test(fileName) ||
        /^ITEM:\s*TIMESTEP/m.test(text.slice(0, 4096));
    const noCompanion = /\.(log|out|xyz|extxyz|trj|gro|pdb|ent|mol|mol2|cif|vasp|poscar|contcar|vesta|cube|gjf|gjf03|gjf09|gjf16|com|mop|mopac|fchk|fch|tcl|md_dump|xdatcar|outcar|pwi|pwo|abo|molden|coord|xml)$/i.test(lowerName) ||
        /^(POSCAR|CONTCAR|OUTCAR|XDATCAR|STRU|MD_dump)$/i.test(fileName);
    if (noCompanion && !lammpsNamed && !coordArray && !lowerName.endsWith('.npy')) return undefined;
    const lammpsLog = /^LAMMPS\s*\(/im.test(text);
    // A `.lammps`/`.lmp` file that is neither data nor a dump is the job's
    // input even when it holds only one command (a lone `include`).
    const lammpsInput = !lammpsLog && (isLammpsInput(text) ||
        (/\.(lmp|lammps)$/i.test(lowerName) && !isLammpsData(text) && !isLammpsDump(text)));
    // `.dat` is MOPAC's input name but also a common LAMMPS data name, and a
    // renamed data/dump file keeps its content. Reading the sibling scripts is
    // what supplies the units and type labels in those cases.
    const lammpsFamily = lammpsInput || isLammpsData(text) || isLammpsDump(text) ||
        /\.(lammpstrj|dump|lmp|lammps|data)$/i.test(lowerName) ||
        /^dump\b/i.test(fileName) || /^data\b/i.test(fileName);
    const cp2kInput = !lammpsInput && /\.(inp|in)$/i.test(lowerName) &&
        (isCp2kInput(text) || /^\s*&(GLOBAL|FORCE_EVAL|SUBSYS|COORD|CELL)\b/im.test(text) ||
         // A main input can hold nothing but preprocessor directives, with
         // every section in an include; CP2K syntax still decides the file.
         /^\s*@(INCLUDE|SET|IF|ELSE|ENDIF|IFDEF|IFNDEF)\b/im.test(text));
    const abacusInput = !lammpsInput && !cp2kInput &&
        (/^INPUT$/i.test(fileName) || (/\.(inp|in)$/i.test(lowerName) && isAbacusInput(text)));
    const abinitInput = !lammpsInput && !cp2kInput && !abacusInput &&
        (/\.abi$/i.test(lowerName) || (/\.(inp|in)$/i.test(lowerName) && isAbinitInput(text)));
    const orcaReferences = !lammpsInput && !cp2kInput && !abacusInput && !abinitInput &&
        /\.(inp|in)$/i.test(lowerName) ? scanOrcaXyzfile(text) : undefined;
    const orcaInput = orcaReferences !== undefined && orcaReferences.length > 0;
    if (!coordArray && !lammpsFamily && !cp2kInput && !abacusInput && !abinitInput &&
        !orcaInput && !lowerName.endsWith('.npy')) return undefined;

    const warnings: string[] = [];
    let reference: StructureReference | undefined;
    let referenceExpected = false;
    let effectiveText: string | undefined;

    const indexes = new Map<string, Map<string, string>>();
    const indexOf = async (folder: vscode.Uri): Promise<Map<string, string>> => {
        const key = folder.toString();
        let index = indexes.get(key);
        if (!index) {
            index = new Map<string, string>();
            try {
                for (const [name, type] of await vscode.workspace.fs.readDirectory(folder)) {
                    if ((type & vscode.FileType.File) !== 0) index.set(name.toLowerCase(), name);
                }
            } catch {
                // An unreadable folder simply has no companions.
            }
            indexes.set(key, index);
        }
        return index;
    };
    const readCompanion = async (folder: vscode.Uri, wanted: string): Promise<{ name: string; bytes: Uint8Array } | undefined> => {
        const actual = (await indexOf(folder)).get(wanted.toLowerCase());
        if (!actual) return undefined;
        try {
            return { name: actual, bytes: await vscode.workspace.fs.readFile(vscode.Uri.joinPath(folder, actual)) };
        } catch {
            return undefined;
        }
    };

    const options: ParseFileOptions = {};
    if (lowerName.endsWith('.npy')) options.binary = bytes;

    if (coordArray) {
        const deepmd: DeepmdParseOptions = {};
        if (lowerName.endsWith('.npy')) deepmd.binary = bytes;
        // `box` accompanies the arrays (a raw dataset has it in the same folder,
        // an npy one keeps it in `set.NNN`), while the type table describes the
        // whole dataset and sits at its root, one level up for `set.NNN`.
        const roots = [dir];
        if (/^set\.\d+$/i.test(path.basename(dir.fsPath))) roots.push(parent);
        const boxFile = await readCompanion(dir, 'box.raw') || await readCompanion(dir, 'box.npy');
        if (boxFile) {
            const np = boxFile.name.toLowerCase().endsWith('.npy');
            const boxText = np ? '' : decodeText(boxFile.bytes);
            const box = deepmdFirstBox(boxText, np ? boxFile.bytes : undefined);
            if (box) deepmd.box = box;
            // The box holds one row per frame, so its row count is what tells a
            // `(natoms, 3)` coordinate array (a single frame) from `(frames, 3)`
            // when `type.raw` is missing.
            const boxFrames = deepmdBoxFrames(boxText, np ? boxFile.bytes : undefined);
            if (boxFrames > 0) deepmd.framesHint = boxFrames;
        }
        for (const root of roots) {
            if (!deepmd.types) {
                const types = await readCompanion(root, 'type.raw');
                if (types) deepmd.types = parseDeepmdTypes(decodeText(types.bytes));
            }
            if (!deepmd.typeMap) {
                const map = await readCompanion(root, 'type_map.raw');
                if (map) deepmd.typeMap = parseDeepmdTypeMap(decodeText(map.bytes));
            }
            if (!deepmd.nopbc && await readCompanion(root, 'nopbc')) deepmd.nopbc = true;
        }
        options.deepmd = deepmd;
    }

    if (lammpsFamily) {
        const hints: LammpsHints = { typeElements: {}, typeMasses: {} };
        let units: string | undefined;
        let atomStyle: string | undefined;
        const index = await indexOf(dir);
        const chainTexts: { content: string; variables: { [name: string]: string } }[] = [];
        const readData: { name: string; fromDir: vscode.Uri }[] = [];
        let readRestart: string | undefined;

        const addChain = async (scriptUri: vscode.Uri, scriptText: string): Promise<void> => {
            const chain = await scanLammpsChain(scriptUri, scriptText, dir);
            chainTexts.push(...chain.texts);
            readData.push(...chain.readData);
            if (!readRestart && chain.readRestart) readRestart = chain.readRestart;
            warnings.push(...chain.warnings);
        };

        if (lammpsInput) {
            // The opened script is the job; its own include chain describes it,
            // not every other input that happens to sit in the folder.
            await addChain(uri, text);
        } else {
            const scriptFiles: { uri: vscode.Uri; text: string }[] = [];
            for (const name of index.values()) {
                if (name.toLowerCase() === lowerName) continue;
                if (!isLammpsScriptName(name)) continue;
                const file = await readCompanion(dir, name);
                if (file) {
                    scriptFiles.push({
                        uri: vscode.Uri.joinPath(dir, file.name),
                        text: decodeText(file.bytes)
                    });
                }
            }
            // A folder often holds one input per job. With several, only the
            // scripts that name this file describe it: an unrelated job's units
            // and masses would rescale or relabel this geometry.
            const escapedName = fileName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
            const namesFile = new RegExp('(^|[^A-Za-z0-9_.-])' + escapedName + '($|[^A-Za-z0-9_.-])', 'i');
            let scripts = scriptFiles.filter(script => namesFile.test(script.text));
            if (scripts.length === 0 && scriptFiles.length === 1) scripts = scriptFiles;
            for (const script of scripts) await addChain(script.uri, script.text);
        }

        for (const item of chainTexts) {
            const parsed = parseLammpsInputHints(item.content, item.variables, units);
            if (!units && parsed.units) units = parsed.units;
            if (!atomStyle && parsed.atomStyle) atomStyle = parsed.atomStyle;
            for (const key of Object.keys(parsed.typeElements)) {
                const type = Number(key);
                if (!hints.typeElements[type]) hints.typeElements[type] = parsed.typeElements[type];
            }
        }

        // Follow read_data in the same command order LAMMPS executes it. The
        // first readable file is the geometry this editor shows; additional
        // data files are reported instead of being silently dropped.
        if (lammpsInput) {
            for (const target of readData) {
                const file = await readReferencedFile(dir, target.fromDir, target.name);
                if (!file) continue;
                reference = file;
                break;
            }
            referenceExpected = readData.length > 0;
            if (reference) {
                options.externalStructureFile = {
                    content: reference.content,
                    fileName: reference.fileName,
                    fsPath: reference.uri.fsPath
                };
                if (readData.length > 1) {
                    warnings.push(readData.length + ' read_data files in the LAMMPS input: showing ' + reference.fileName);
                }
            } else if (readData.length > 0) {
                warnings.push('the LAMMPS input reads its structure from "' + readData[0].name +
                    '", but that file could not be read');
            } else if (readRestart) {
                warnings.push('the LAMMPS input continues from the restart file "' + readRestart +
                    '": a restart file holds no readable atom list');
            }
        }

        // The masses in the data file name its types too, and following
        // `read_data` reaches the file the script itself used. Without that
        // pointer, several sibling data files describe several jobs: merging
        // all their masses would relabel this dump with another job's types.
        let dataText: string | undefined;
        if (lammpsInput) {
            if (reference) dataText = reference.content;
        } else {
            for (const target of readData) {
                const file = await readReferencedFile(dir, target.fromDir, target.name);
                if (file) { dataText = file.content; break; }
            }
            if (dataText === undefined) {
                const candidates: string[] = [];
                for (const name of index.values()) if (/^data\b/i.test(name)) candidates.push(name);
                if (candidates.length === 1) {
                    const file = await readCompanion(dir, candidates[0]);
                    if (file) dataText = decodeText(file.bytes);
                }
            }
        }
        if (dataText !== undefined) {
            const masses = lammpsTypeElements(dataText, units);
            for (const key of Object.keys(masses)) {
                const type = Number(key);
                if (!hints.typeElements[type]) hints.typeElements[type] = masses[type];
            }
        }
        if (units || atomStyle || Object.keys(hints.typeElements).length > 0) {
            const lammps: LammpsParseOptions = {};
            if (units) lammps.units = units;
            if (atomStyle) lammps.atomStyle = atomStyle;
            if (Object.keys(hints.typeElements).length > 0) lammps.typeElements = hints.typeElements;
            options.lammps = lammps;
        }
    }

    if (cp2kInput) {
        const expanded = await expandCp2kInput(uri, text, dir);
        warnings.push(...expanded.warnings);
        effectiveText = expanded.text;
        // CP2K keeps the last value of a key; a file can repeat TOPOLOGY.
        const coordPattern = /^\s*COORD_FILE_NAME\s+(.+?)\s*$/gim;
        let coordMatch: RegExpExecArray | null;
        let coordLine: string | undefined;
        while ((coordMatch = coordPattern.exec(expanded.text)) !== null) coordLine = coordMatch[1];
        if (coordLine !== undefined) {
            const wanted = directiveValue(coordLine);
            referenceExpected = true;
            const file = wanted === '' ? undefined : await readReferencedFromDirs(dir, expanded.dirs, wanted);
            if (file) {
                reference = file;
                options.externalCoordinateFile = { content: file.content, fileName: file.fileName };
            } else {
                warnings.push('the CP2K input reads its structure from "' + wanted +
                    '", but that file could not be read');
            }
        }
    }

    if (abacusInput) {
        referenceExpected = true;
        const scan = scanAbacusInput(text);
        // ABACUS defaults the structure name to `STRU` when INPUT omits it.
        const wanted = scan.struFile || 'STRU';
        const file = await readReferencedFile(dir, dir, wanted);
        if (file) {
            reference = file;
            options.externalStructureFile = {
                content: file.content,
                fileName: file.fileName,
                fsPath: file.uri.fsPath
            };
        } else {
            warnings.push('the ABACUS INPUT reads its structure from "' + wanted +
                '", but that file could not be read');
        }
    }

    if (abinitInput) {
        const expanded = await expandAbinitInput(uri, text, dir);
        warnings.push(...expanded.warnings);
        effectiveText = expanded.text;
        const scan = scanAbinitXyzfile(expanded.text);
        if (scan.xyzfile) {
            referenceExpected = true;
            const file = await readReferencedFromDirs(dir, expanded.dirs, scan.xyzfile);
            if (file) {
                reference = file;
                options.externalStructureFile = {
                    content: file.content,
                    fileName: file.fileName,
                    fsPath: file.uri.fsPath
                };
            } else {
                warnings.push('the ABINIT input reads its structure from "' + scan.xyzfile +
                    '", but that file could not be read');
            }
        }
    }

    if (orcaInput && orcaReferences) {
        const wanted = orcaReferences.length > 0 ? orcaReferences[0].name : undefined;
        if (wanted !== undefined) {
            referenceExpected = true;
            const file = await readReferencedFile(dir, dir, wanted);
            if (file) {
                reference = file;
                options.externalStructureFile = {
                    content: file.content,
                    fileName: file.fileName,
                    fsPath: file.uri.fsPath
                };
            } else {
                warnings.push('the ORCA input reads its structure from "' + wanted +
                    '", but that file could not be read');
            }
        }
    }

    return { options, warnings, reference, referenceExpected, text: effectiveText };
}
