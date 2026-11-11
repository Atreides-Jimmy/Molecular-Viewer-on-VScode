import { MolecularData, AtomGroup, OptStep, NormalMode, OptThresholds, OrbitalLevels, SolventCavity } from '../types';
import { parseGjf } from './gjfParser';
import { parseXyz, parseXyzTrajectory, XyzFrame, XyzTrajectory } from './xyzParser';
import { parseMol2 } from './mol2Parser';
import { parseMol, isMdlMol } from './molParser';
import { parseGaussianLog, LogFrame, GaussianLogResult, RouteSection, RouteKeyword, LogChargeSet, LogSummary } from './logParser';
import { parseCoord } from './coordParser';
import { parseOrcainp, scanOrcaXyzfile, OrcaExternalStructure, OrcaXyzReference } from './orcaInpParser';
import { parseOrcaOut, isOrcaOutput, OrcaFrame, OrcaLogResult } from './orcaOutParser';
import { parsePdb } from './pdbParser';
import { parseMopac, parseMopacOutput, isMopacOutput } from './mopacParser';
import { parseCif } from './cifParser';
import { parseVasp, parseXdatcar, parseOutcar, parseVasprunXml, isXdatcar, isOutcar, isVasprunXml } from './vaspParser';
import { parseCube } from './cubeParser';
import { parseVesta } from './vestaParser';
import { parseXtbLog, parseXtbOutput, isXtbLog, isXtbOutput, XtbLogResult } from './xtbParser';
import { parseTcl, TclParseResult, resolveVmdScene, VmdResolvedScene, hasScalarChannel, resolveSelection, vdwRadius, parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg } from './vmdVisualization';
import { parseAbacusStru, parseAbacusLog, parseAbacusMdDump, parseAbacusInput, scanAbacusInput, isAbacusLog, isAbacusStru, isAbacusMdDump, isAbacusInput, AbacusExternalStructure, AbacusInputScan } from './abacusParser';
import { parseAbinitInput, parseAbinitOutput, scanAbinitXyzfile, isAbinitOutput, isAbinitInput, AbinitExternalStructure, AbinitXyzfileScan } from './abinitParser';
import { parseCp2kInput, parseCp2kOutput, isCp2kInput, isCp2kOutput, Cp2kParseOptions } from './cp2kParser';
import { parseQeInput, parseQeOutput, isQeInput, isQeOutput } from './qeParser';
import { parseLammpsData, parseLammpsDump, parseLammpsDumpTrajectory, isLammpsData, isLammpsDump, isLammpsInput, isLammpsLog, parseLammpsInputStructure, scanLammpsInput, LammpsParseOptions, LammpsHints, LammpsInputScan, parseLammpsInputHints, lammpsTypeElements, lammpsLengthFactor, LammpsDumpFrame, LammpsDumpTrajectory, LammpsDumpTrajectoryOptions } from './lammpsParser';
import { parseGro, isGro, isGromacsLog } from './gromacsParser';
import { parseExtXyz, isExtXyz, isGpumdInput } from './extXyzParser';
import { parseDeepmdDataset, isDeepmdDatasetFile, isDpgenModelDevi, DeepmdParseOptions, parseDeepmdTypes, parseDeepmdTypeMap, deepmdFirstBox } from './deepmdParser';
import { parsePyscfOutput, parseMolden, isPyscfOutput, isMolden } from './pyscfParser';
import { parseFchk, parseFchkText, readFchkColumn, FchkResult, FchkOrbitals, FchkArrayLayout, FchkOrbitalSpin, FchkBasis } from './fchkParser';
import { detectBinaryTrajectory, readBinaryTrajectory, isTrajectoryFileName, looksLikeBinaryTrajectory, BinaryTrajectory, TrajectoryFrame, TrajectoryFormat, TrajectoryReadOptions } from './trajectoryParser';

export { parseGaussianLog, LogFrame, GaussianLogResult };
export { RouteSection, RouteKeyword };
export { LogChargeSet, LogSummary };
export { parseOrcaOut, isOrcaOutput, OrcaFrame, OrcaLogResult };
export { scanOrcaXyzfile, OrcaExternalStructure, OrcaXyzReference };
export { parseXyzTrajectory, XyzFrame, XyzTrajectory };
export { parseCif };
export { parseMol, isMdlMol };
export { parseVasp, parseXdatcar, parseOutcar, parseVasprunXml, isXdatcar, isOutcar, isVasprunXml };
export { parseCube };
export { parseVesta };
export { parseXtbLog, parseXtbOutput, isXtbLog, isXtbOutput, XtbLogResult };
export { parseMopac, parseMopacOutput, isMopacOutput };
export { parseAbacusStru, parseAbacusLog, parseAbacusMdDump, parseAbacusInput, scanAbacusInput, isAbacusLog, isAbacusStru, isAbacusMdDump, isAbacusInput, AbacusExternalStructure, AbacusInputScan };
export { parseAbinitInput, parseAbinitOutput, scanAbinitXyzfile, isAbinitOutput, isAbinitInput, AbinitExternalStructure, AbinitXyzfileScan };
export { parseCp2kInput, parseCp2kOutput, isCp2kInput, isCp2kOutput, Cp2kParseOptions };
export { parseQeInput, parseQeOutput, isQeInput, isQeOutput };
export { parseLammpsData, parseLammpsDump, parseLammpsDumpTrajectory, isLammpsData, isLammpsDump, isLammpsInput, isLammpsLog, parseLammpsInputStructure, scanLammpsInput, LammpsParseOptions, LammpsHints, LammpsInputScan, parseLammpsInputHints, lammpsTypeElements, lammpsLengthFactor, LammpsDumpFrame, LammpsDumpTrajectory, LammpsDumpTrajectoryOptions };
export { parseGro, isGro, isGromacsLog };
export { parseExtXyz, isExtXyz, isGpumdInput };
export { parseDeepmdDataset, isDeepmdDatasetFile, isDpgenModelDevi, DeepmdParseOptions, parseDeepmdTypes, parseDeepmdTypeMap, deepmdFirstBox };
export { parsePyscfOutput, parseMolden, isPyscfOutput, isMolden };
export { parseFchk, parseFchkText, readFchkColumn, FchkResult, FchkOrbitals, FchkArrayLayout, FchkOrbitalSpin, FchkBasis };
export { detectBinaryTrajectory, readBinaryTrajectory, isTrajectoryFileName, looksLikeBinaryTrajectory, BinaryTrajectory, TrajectoryFrame, TrajectoryFormat, TrajectoryReadOptions };
// Visualization files: the VMD scripts and the bare scalar-carrying coordinate
// files they describe (see src/parsers/vmdVisualization.ts).
export { parseTcl, TclParseResult };
export { resolveVmdScene, VmdResolvedScene };
export { hasScalarChannel };
export { resolveSelection, vdwRadius };
export { parseUnitFile, parseEggFile, matchEggsToSegments, VmdEgg };

export interface LogFileResult {
    frames: LogFrame[] | OrcaFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    // Route-card sections (Gaussian's keywords lines, ORCA's simple-input line),
    // shown in the webview's Route panel. Other formats leave it undefined.
    routes?: RouteSection[];
    // The convergence targets the file itself states, when it states them.
    // The panel draws them as dashed target lines; a program that states none
    // (xtb's gnorm criterion, a trajectory with no table) gets no lines rather
    // than a borrowed scale.
    optThresholds?: OptThresholds;
    // Which log format produced this result. The webview uses it to label the
    // panel; the numeric scales come from `optThresholds`.
    source?: 'gaussian' | 'xtb' | 'orca';
    // The per-atom charge tables the file prints (a Gaussian population
    // analysis) and the headline results of each job. Both come from the file
    // as it states them; a format that prints neither leaves them undefined,
    // which is what keeps the charge and summary panels from appearing on a
    // log that has no such data.
    charges?: LogChargeSet[];
    summary?: LogSummary[];
    // The last eigenvalue dump of the job (energies and occupations, never
    // coefficients) and the solvation cavity it describes. A log states both;
    // only a formatted checkpoint states the coefficients beside them.
    orbitals?: OrbitalLevels;
    cavity?: SolventCavity;
}

export interface ParseFileOptions {
    /** Companion coordinate file referenced by a CP2K `COORD_FILE_NAME`. */
    externalCoordinateFile?: { content: string; fileName: string };
    /**
     * Structure file a main input names (LAMMPS `read_data`, ORCA/ABINIT
     * `xyzfile`, ABACUS `stru_file`). The host resolves the path and reads the
     * bytes; `fsPath` is returned here so the viewer can save and trace the
     * structure beside its real file.
     */
    externalStructureFile?: { content: string; fileName: string; fsPath?: string };
    /** Type/element hints and units read from a LAMMPS input or data file. */
    lammps?: LammpsParseOptions;
    /** Type map, types and box read from the rest of a DeepMD dataset. */
    deepmd?: DeepmdParseOptions;
    /** Raw bytes when the primary file is binary (`.npy`). */
    binary?: Uint8Array;
    /**
     * Internal: the caller has already run the structured-output detectors and
     * they found nothing. `tryAutoParse` then skips its own copy of that scan,
     * which is what keeps a large non-structured output from being walked a
     * second time per `parseFile` call.
     */
    structuredProbed?: boolean;
}

/**
 * Whether the text is a frame log this reader parses as geometry frames
 * (an xTB trajectory, a Gaussian log or an ORCA output).
 *
 * These are the three outputs that share `.log`/`.out` with the structured
 * programs below. Their markers are found with `indexOf` or the first two
 * lines, so the test is cheap enough to run before a 60 MB scan.
 */
export function isFrameLogOutput(content: string): boolean {
    return isXtbLog(content) ||
        content.indexOf('Standard orientation:') >= 0 ||
        content.indexOf('Input orientation:') >= 0 ||
        content.indexOf('Gaussian, Inc.') >= 0 ||
        isOrcaOutput(content);
}

/**
 * The program outputs this reader understands, in detection order. One list
 * feeds both the cheap "is this structured output?" test and the parser, so a
 * caller that probes and then parses cannot check a different set of formats
 * from the one it ends up reading.
 */
interface StructuredReader {
    test: (content: string) => boolean;
    read: (content: string, options?: ParseFileOptions) => MolecularData;
}

/** A recognised program output that carries no coordinates. */
function noStructureOutput(title: string, warning: string): MolecularData {
    return { atoms: [], bonds: [], title, hasExplicitBonds: false, warnings: [warning] };
}

const STRUCTURED_READERS: StructuredReader[] = [
    { test: isQeOutput, read: content => parseQeOutput(content) },
    { test: isCp2kOutput, read: content => parseCp2kOutput(content) },
    { test: isPyscfOutput, read: content => parsePyscfOutput(content) },
    { test: isAbinitOutput, read: content => parseAbinitOutput(content) },
    { test: isVasprunXml, read: content => parseVasprunXml(content) },
    { test: isOutcar, read: content => parseOutcar(content) },
    { test: isLammpsDump, read: (content, options) => parseLammpsDump(content, options && options.lammps) },
    { test: isAbacusMdDump, read: content => parseAbacusMdDump(content) },
    { test: isAbacusLog, read: content => parseAbacusLog(content) },
    { test: isXdatcar, read: content => parseXdatcar(content) },
    // MOPAC and xTB print their final geometry in output that shares `.out`
    // and `.log` with Gaussian, ORCA and the trajectory formats.
    { test: isMopacOutput, read: content => parseMopacOutput(content) },
    { test: isXtbOutput, read: content => parseXtbOutput(content) },
    // Outputs that share `.log`/`.out` with the structure-bearing programs.
    // Naming them here keeps the Gaussian/ORCA/xTB readers from claiming them
    // and replaces an empty canvas with the file that does hold the geometry.
    {
        test: isGromacsLog,
        read: () => noStructureOutput('GROMACS log',
            'a GROMACS log records energies, temperatures and performance, not coordinates: ' +
            'open the .gro file mdrun wrote (the -deffnm output, e.g. water.gro)')
    },
    {
        test: isLammpsLog,
        read: () => noStructureOutput('LAMMPS log',
            'a LAMMPS log records the run and its thermo output, not coordinates: ' +
            'open the data file or a dump file')
    },
    {
        test: isDpgenModelDevi,
        read: () => noStructureOutput('DP-GEN model deviation',
            'model_devi.out is a model-deviation table, not a structure: ' +
            'open the corresponding POSCAR/XYZ frames or the LAMMPS dump')
    },
    {
        test: isGpumdInput,
        read: () => noStructureOutput('GPUMD input',
            'a GPUMD run.in selects the potential and the files to dump, not a structure: ' +
            'open model.xyz or the dump.xyz it writes')
    }
];

/**
 * Whether a `.log`/`.out` file is structured output of a program this reader
 * understands, rather than a Gaussian/ORCA/xTB trajectory log. The host uses
 * this to choose between the frame pipeline and `parseFile`.
 */
export function isStructuredProgramOutput(content: string, onProgress?: (fraction: number) => void): boolean {
    // A frame log is not structured output. Its markers are cheap, so checking
    // them first keeps a 60 MB file from being scanned by every detector below.
    if (isFrameLogOutput(content)) {
        onProgress?.(1);
        return false;
    }
    // Each detector is a full pass over the text, so one step per detector is
    // the only progress a caller can honestly be shown before the result is
    // known. The short-circuit is unchanged: the first match wins.
    for (let i = 0; i < STRUCTURED_READERS.length; i++) {
        const matched = STRUCTURED_READERS[i].test(content);
        onProgress?.((i + 1) / STRUCTURED_READERS.length);
        if (matched) return true;
    }
    return false;
}

/**
 * Detect and read structured program output in one pass.
 *
 * A caller that needs the structure (the open path, Diff and Import) used to
 * ask `isStructuredProgramOutput` and then hand the same text to `parseFile`,
 * whose `.log` branch ran the detector list a second time. This runs the
 * detectors once and returns the parsed result, so a large QE/CP2K/LAMMPS
 * output is scanned once rather than twice.
 */
export function parseStructuredProgramOutput(
    content: string,
    options?: ParseFileOptions,
    onProgress?: (fraction: number) => void
): MolecularData | null {
    if (isFrameLogOutput(content)) {
        onProgress?.(1);
        return null;
    }
    for (let i = 0; i < STRUCTURED_READERS.length; i++) {
        const reader = STRUCTURED_READERS[i];
        if (reader.test(content)) {
            const data = reader.read(content, options);
            onProgress?.(1);
            return data;
        }
        onProgress?.((i + 1) / STRUCTURED_READERS.length);
    }
    return null;
}

/** Recognise structured program output before falling back to a Gaussian/ORCA log. */
function structuredOutput(content: string, options?: ParseFileOptions): MolecularData | null {
    for (const reader of STRUCTURED_READERS) {
        if (reader.test(content)) return reader.read(content, options);
    }
    return null;
}

/**
 * A Gaussian/ORCA/xTB log opened under an input extension (`.in`/`.inp`).
 *
 * The host sends `.log`/`.out` files to the frame pipeline, so a renamed log
 * normally never reaches `parseFile`. When one does, show its first structure
 * with the program's own title instead of letting the ORCA input reader claim
 * the file and return nothing.
 */
function singleFrameLogOutput(content: string): MolecularData | null {
    const isLog = isOrcaOutput(content) || isXtbLog(content) ||
        content.includes('Standard orientation:') || content.includes('Input orientation:');
    return isLog ? parseLogAsSingleFrame(content) : null;
}

export function parseFile(content: string, fileName: string, options?: ParseFileOptions): MolecularData {
    const baseName = fileName.split(/[\\/]/).pop() || fileName;
    const ext = baseName.toLowerCase().split('.').pop() || '';

    // DeepMD datasets and LAMMPS data files are named by convention, not by a
    // unique extension (`coord.raw`, `data.si`), so the name is consulted
    // before the switch.
    if (/^coord\.(raw|npy)$/i.test(baseName)) {
        const deepmd: DeepmdParseOptions = { ...(options && options.deepmd ? options.deepmd : {}) };
        if (options && options.binary && !deepmd.binary) deepmd.binary = options.binary;
        return parseDeepmdDataset(content, baseName, deepmd);
    }
    if (/^data\b/i.test(baseName) && isLammpsData(content)) {
        return parseLammpsData(content, options && options.lammps);
    }
    // ABACUS names its parameter file exactly `INPUT`; that name and the
    // `stru_file` key together leave no other format this could be.
    if (/^INPUT$/i.test(baseName) &&
        (isAbacusInput(content) || (options && options.externalStructureFile))) {
        return parseAbacusInput(content, options && options.externalStructureFile);
    }

    switch (ext) {
        case 'gjf':
        case 'gjf03':
        case 'gjf09':
        case 'gjf16':
        case 'com':
            return parseGjf(content);
        case 'xyz':
        case 'trj':
            // A GPUMD/ASE file carries its cell on the comment line; the plain
            // XYZ reader would throw that cell away. xTB writes its trajectory
            // as `xtb.trj`, which is the same back-to-back XYZ layout.
            if (ext === 'trj') {
                // CHARMM/NAMD also write binary trajectories with this suffix.
                // When the bytes are not one of the formats `readBinaryTrajectory`
                // recognises, the host has already avoided decoding them; this
                // guard is for a caller that handed over the decoded text anyway.
                const head = content.slice(0, 8192);
                const replacements = (head.match(/\uFFFD/g) || []).length;
                if (/[\u0000-\u0008\u000b\u000c\u000e-\u001f]/.test(head) || replacements > 8) {
                    return noStructureOutput('Binary trajectory',
                        'this .trj file is binary and is not an XTC, TRR or DCD trajectory: ' +
                        'export it as DCD/XTC/TRR or as a text trajectory (e.g. .xyz)');
                }
            }
            return isExtXyz(content) ? parseExtXyz(content) : parseXyz(content);
        case 'extxyz':
            return parseExtXyz(content);
        case 'mol2':
            return parseMol2(content);
        case 'mol':
            // Save As writes MDL Mol, so the same reader opens it back. A file
            // whose header was damaged still gets the warning from `parseMol`
            // rather than the empty XYZ structure the fallback produced.
            return parseMol(content);
        case 'log':
        case 'out': {
            // The caller may already have probed the structured detectors. When
            // it has not, run them here; when it has, do not scan the text a
            // second time.
            const structured = options && options.structuredProbed ? null : structuredOutput(content, options);
            if (structured) return structured;
            // An xTB trajectory log has an XYZ block per frame and would
            // otherwise be read by the plain XYZ fallback below, losing the
            // frame's charge/multiplicity. Keep its own frame reader first.
            if (isXtbLog(content)) return parseLogAsSingleFrame(content);
            const probed: ParseFileOptions = options ? { ...options, structuredProbed: true } : { structuredProbed: true };
            // Old MOPAC outputs state no version banner, so the detector that
            // knows their Cartesian table lives in the content sniffers. A
            // Gaussian/ORCA/xTB frame log reaches its own reader there as well;
            // `tryAutoParse` only says nothing (no atoms, no warning) when no
            // reader recognised the file at all, which is when this branch
            // reports the frame reader's own explanation.
            const guessed = tryAutoParse(content, fileName, probed);
            if (guessed.atoms.length > 0 || (guessed.warnings && guessed.warnings.length > 0)) {
                return guessed;
            }
            return parseLogAsSingleFrame(content);
        }
        case 'coord':
            return parseCoord(content);
        case 'inp':
        case 'in': {
            // Programs echo their input into their output. A CP2K output is
            // full of `&SECTION` lines and a QE/ABINIT/PySCF output can look
            // like a namelist, so what wrote the file has to win before the
            // input syntax it repeats does. MOPAC and xTB outputs share these
            // extensions too when a job renamed them.
            const structured = options && options.structuredProbed ? null : structuredOutput(content, options);
            if (structured) return structured;
            const logOutput = singleFrameLogOutput(content);
            if (logOutput) return logOutput;
            // `.inp` is ORCA's input extension, but Gaussian and CP2K users
            // name their inputs that way too. The CP2K section syntax is
            // unambiguous; ORCA's is `* xyz`/`* xyzfile`/`%coords`, and
            // anything else that starts with a route card is Gaussian.
            if (isCp2kInput(content)) return parseCp2kInput(content, options);
            if (isQeInput(content)) return parseQeInput(content);
            if (isAbacusInput(content)) {
                return parseAbacusInput(content, options && options.externalStructureFile);
            }
            if (isAbinitInput(content)) {
                return parseAbinitInput(content, options && options.externalStructureFile);
            }
            if (isLammpsInput(content)) {
                return parseLammpsInputStructure(content, options && options.externalStructureFile, options && options.lammps);
            }
            const orca = parseOrcainp(content, options && options.externalStructureFile);
            if (orca.atoms.length > 0) return orca;
            // `.in` is not Gaussian's extension, and the scripts that do use
            // it start with comments (`#` in GPUMD's run.in), so the content
            // sniffers get the file before the route-card fallback. A renamed
            // `.inp` gets the same treatment: QE, PySCF and ABINIT inputs and
            // outputs are not ORCA inputs merely because ORCA uses the suffix.
            const probed: ParseFileOptions = options ? { ...options, structuredProbed: true } : { structuredProbed: true };
            const guessed = tryAutoParse(content, fileName, probed);
            if (guessed.atoms.length === 0 && (!guessed.warnings || guessed.warnings.length === 0)) {
                guessed.warnings = ['no structure was found in this input file: ' +
                    'check that it contains coordinates or names the file that does'];
            }
            return guessed;
        }
        case 'pdb':
        case 'ent':
            return parsePdb(content);
        case 'mop':
        case 'mopac':
            return isMopacOutput(content) ? parseMopacOutput(content) : parseMopac(content);
        case 'dat': {
            // `.dat` is MOPAC's input name, but LAMMPS data files use it too
            // and MOPAC output can be renamed to it.
            if (isLammpsData(content)) return parseLammpsData(content, options && options.lammps);
            if (isLammpsDump(content)) return parseLammpsDump(content, options && options.lammps);
            if (isMopacOutput(content)) return parseMopacOutput(content);
            return parseMopac(content);
        }
        case 'cif':
            return parseCif(content);
        case 'vasp':
        case 'poscar':
        case 'contcar':
            return parseVasp(content);
        case 'cube':
            return parseCube(content);
        case 'fchk':
        case 'fch': {
            // The streaming reader is used by every path that can open a large
            // checkpoint; this branch exists for a checkpoint that reached
            // `parseFile` as text (a renamed file, a small one). Only the
            // structure travels here — the orbitals, charges and cavity need
            // the dedicated document loader.
            const result = parseFchkText(content);
            return {
                atoms: result.atoms,
                bonds: [],
                title: /^molecule name$/i.test(result.title.trim()) ? fileName : result.title,
                hasExplicitBonds: false,
                charge: result.charge,
                multiplicity: result.multiplicity,
                warnings: result.warnings.length > 0 ? result.warnings.slice() : undefined
            };
        }
        case 'vesta':
            return parseVesta(content);
        case 'abi':
            return parseAbinitInput(content, options && options.externalStructureFile);
        case 'abo':
            return parseAbinitOutput(content);
        case 'stru':
            return parseAbacusStru(content);
        case 'md_dump':
            return parseAbacusMdDump(content);
        case 'gro':
            return parseGro(content);
        case 'molden':
            return parseMolden(content);
        case 'pwi':
            return parseQeInput(content);
        case 'pwo':
            return parseQeOutput(content);
        case 'lammpstrj':
            return parseLammpsDump(content, options && options.lammps);
        case 'data':
            // `system.data` is the usual LAMMPS data name; any other `.data`
            // file is not a structure this viewer knows.
            if (isLammpsData(content)) return parseLammpsData(content, options && options.lammps);
            return tryAutoParse(content, fileName, options);
        case 'lmp':
        case 'lammps':
            if (isLammpsData(content)) return parseLammpsData(content, options && options.lammps);
            if (isLammpsDump(content)) return parseLammpsDump(content, options && options.lammps);
            // A log echoes its input commands; it is not the input file, and
            // there is no LAMMPS log reader to fall through to.
            if (/^LAMMPS\s*\(/im.test(content)) return tryAutoParse(content, fileName, options);
            // The extension names a LAMMPS script; even a one-line `include`
            // main file must not fall through to the coordinate guessers.
            return parseLammpsInputStructure(content, options && options.externalStructureFile, options && options.lammps);
        case 'xdatcar':
            return parseXdatcar(content);
        case 'outcar':
            return parseOutcar(content);
        case 'xml':
            if (isVasprunXml(content)) return parseVasprunXml(content);
            return tryAutoParse(content, fileName, options);
        case 'npy':
        case 'raw':
            // Only a coordinate file is a structure; any other array is not
            // something this viewer can draw.
            return tryAutoParse(content, fileName, options);
        default:
            return tryAutoParse(content, fileName, options);
    }
}

/**
 * Parse a log/trajectory file into the frames and the analysis data the viewer
 * needs.
 *
 * The format is decided by the CONTENT, not by the extension: Gaussian and ORCA
 * are both routinely told to write `.out`, and xTB shares `.log` with Gaussian.
 * Every detector below matches a string only one program writes, so the sniffing
 * cannot pick the wrong parser; the extension is not consulted at all, which
 * means a mis-named file (an ORCA run saved as `.log`, a Gaussian log saved as
 * `.out`) still opens correctly.
 */
export function parseLogFile(content: string, fileName?: string, onProgress?: (fraction: number) => void): LogFileResult {
    if (isOrcaOutput(content)) {
        const result = parseOrcaOut(content, onProgress);
        return {
            frames: result.frames,
            title: result.title,
            optSteps: result.optSteps,
            normalModes: result.normalModes,
            routes: result.routes,
            optThresholds: result.optThresholds,
            source: 'orca'
        };
    }
    if (isXtbLog(content)) {
        const xResult = parseXtbLog(content, onProgress);
        return {
            frames: xResult.frames,
            title: xResult.title,
            optSteps: xResult.optSteps,
            normalModes: xResult.normalModes,
            source: 'xtb'
        };
    }
    const gResult = parseGaussianLog(content, onProgress);
    return {
        frames: gResult.frames,
        title: gResult.title,
        optSteps: gResult.optSteps,
        normalModes: gResult.normalModes,
        routes: gResult.routes,
        source: 'gaussian',
        charges: gResult.charges,
        summary: gResult.summary,
        orbitals: gResult.orbitals,
        cavity: gResult.cavity
    };
}

function parseLogAsSingleFrame(content: string): MolecularData {
    const result = parseLogFile(content);
    if (result.frames.length > 0) {
        const f = result.frames[0];
        return {
            atoms: f.atoms,
            bonds: f.bonds,
            title: f.title,
            hasExplicitBonds: f.hasExplicitBonds,
            charge: f.charge,
            multiplicity: f.multiplicity
        };
    }
    return {
        atoms: [], bonds: [], title: 'Empty', hasExplicitBonds: false,
        warnings: ['no Gaussian, ORCA or xTB structure frame was found in this file: ' +
            'open the output file that holds the geometry']
    };
}

export function tryAutoParse(content: string, fileName: string = '', options?: ParseFileOptions): MolecularData {
    const lines = content.split(/\r?\n/).filter(l => l.trim() !== '');

    // Program outputs echo their own inputs, so what wrote the file has to win
    // over the input syntax it repeats. An ORCA output contains its `* xyz`
    // block; reading that echo as an ORCA input returned no atoms and no
    // warning, so structured output is detected first. A caller that already
    // ran the detectors (`structuredProbed`) must not pay for them again: this
    // is the second scan a `.log`/`.out` used to make inside one `parseFile`.
    if (!options || !options.structuredProbed) {
        const structured = structuredOutput(content, options);
        if (structured) return structured;
    }

    // ORCA output is recognised by ORCA-only markers, not by its echoed input.
    // A renamed output (a `.txt` or an extensionless file) therefore opens
    // with its geometry instead of being read as the input it repeats.
    if (isOrcaOutput(content)) {
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
        return {
            atoms: [], bonds: [], title: result.title || 'ORCA output', hasExplicitBonds: false,
            warnings: ['this is an ORCA output but no geometry frame was found in it: ' +
                'open the output file that holds the coordinates']
        };
    }

    if (/^[ \t]*\$coord\b/im.test(content)) {
        return parseCoord(content);
    }

    if (content.match(/\*\s*xyz/i) || content.match(/\*\s*xyzfile/i)) {
        return parseOrcainp(content, options && options.externalStructureFile);
    }

    // The extended XYZ comment line is the only thing that distinguishes a
    // GPUMD/ASE trajectory from a plain XYZ; it must be tested before the
    // first-line-is-a-count fallback below.
    if (isExtXyz(content)) {
        return parseExtXyz(content);
    }

    if (isLammpsInput(content)) {
        return parseLammpsInputStructure(content, options && options.externalStructureFile, options && options.lammps);
    }

    if (isMolden(content)) {
        return parseMolden(content);
    }

    // GROMACS names its structures `.gro`, but a backup (`#water.gro.1#`) and
    // a renamed file keep the format, not the extension. The fixed-width
    // residue/atom/serial/xyz layout is distinctive enough to recognise.
    if (isGro(content)) {
        return parseGro(content);
    }

    // PDB's record name is followed by the fixed-width serial; matching a bare
    // `ATOM` prefix would claim ABACUS's `ATOMIC_SPECIES` as a PDB file.
    if (content.match(/^(ATOM  |HETATM)/m)) {
        return parsePdb(content);
    }

    // An older MOPAC output may not state a version banner; the program name
    // plus the Cartesian table still identifies it, and the output reader
    // understands the serial number in front of each element.
    if (/^\s*CARTESIAN COORDINATES\b/im.test(content) && /\bMOPAC\b/i.test(content)) {
        return parseMopacOutput(content);
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

    if (isAbacusStru(content)) {
        return parseAbacusStru(content);
    }

    if (isLammpsData(content)) {
        return parseLammpsData(content, options && options.lammps);
    }

    if (isAbacusInput(content)) {
        return parseAbacusInput(content, options && options.externalStructureFile);
    }

    if (isCp2kInput(content)) {
        return parseCp2kInput(content, options);
    }

    if (isQeInput(content)) {
        return parseQeInput(content);
    }

    if (isAbinitInput(content)) {
        return parseAbinitInput(content, options && options.externalStructureFile);
    }

    // A molfile has no extension of its own in some toolchains; the CTfile
    // counts line and `M  END` identify it before the first-line-is-a-count
    // fallback can claim the title line.
    if (isMdlMol(content)) {
        return parseMol(content);
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
