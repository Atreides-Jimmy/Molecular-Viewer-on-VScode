import { MolecularData, AtomGroup, OptStep, NormalMode, OptThresholds } from '../types';
import { parseGjf } from './gjfParser';
import { parseXyz, parseXyzTrajectory, XyzFrame, XyzTrajectory } from './xyzParser';
import { parseMol2 } from './mol2Parser';
import { parseMol, isMdlMol } from './molParser';
import { parseGaussianLog, LogFrame, GaussianLogResult, RouteSection, RouteKeyword } from './logParser';
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
import { detectBinaryTrajectory, readBinaryTrajectory, isTrajectoryFileName, BinaryTrajectory, TrajectoryFrame, TrajectoryFormat, TrajectoryReadOptions } from './trajectoryParser';

export { parseGaussianLog, LogFrame, GaussianLogResult };
export { RouteSection, RouteKeyword };
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
export { detectBinaryTrajectory, readBinaryTrajectory, isTrajectoryFileName, BinaryTrajectory, TrajectoryFrame, TrajectoryFormat, TrajectoryReadOptions };
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
}

/**
 * Whether a `.log`/`.out` file is structured output of a program this reader
 * understands, rather than a Gaussian/ORCA/xTB trajectory log. The host uses
 * this to choose between the frame pipeline and `parseFile`.
 */
export function isStructuredProgramOutput(content: string): boolean {
    return isQeOutput(content) || isCp2kOutput(content) || isPyscfOutput(content) ||
        isAbinitOutput(content) || isVasprunXml(content) || isOutcar(content) ||
        isLammpsDump(content) || isAbacusMdDump(content) || isAbacusLog(content) ||
        isXdatcar(content) || isMopacOutput(content) || isXtbOutput(content) ||
        isGromacsLog(content) || isLammpsLog(content) ||
        isDpgenModelDevi(content) || isGpumdInput(content);
}

/** A recognised program output that carries no coordinates. */
function noStructureOutput(title: string, warning: string): MolecularData {
    return { atoms: [], bonds: [], title, hasExplicitBonds: false, warnings: [warning] };
}

/** Recognise structured program output before falling back to a Gaussian/ORCA log. */
function structuredOutput(content: string, options?: ParseFileOptions): MolecularData | null {
    if (isQeOutput(content)) return parseQeOutput(content);
    if (isCp2kOutput(content)) return parseCp2kOutput(content);
    if (isPyscfOutput(content)) return parsePyscfOutput(content);
    if (isAbinitOutput(content)) return parseAbinitOutput(content);
    if (isVasprunXml(content)) return parseVasprunXml(content);
    if (isOutcar(content)) return parseOutcar(content);
    if (isLammpsDump(content)) return parseLammpsDump(content, options && options.lammps);
    if (isAbacusMdDump(content)) return parseAbacusMdDump(content);
    if (isAbacusLog(content)) return parseAbacusLog(content);
    if (isXdatcar(content)) return parseXdatcar(content);
    // MOPAC and xTB print their final geometry in output that shares `.out`
    // and `.log` with Gaussian, ORCA and the trajectory formats.
    if (isMopacOutput(content)) return parseMopacOutput(content);
    if (isXtbOutput(content)) return parseXtbOutput(content);
    // Outputs that share `.log`/`.out` with the structure-bearing programs.
    // Naming them here keeps the Gaussian/ORCA/xTB readers from claiming them
    // and replaces an empty canvas with the file that does hold the geometry.
    if (isGromacsLog(content)) {
        return noStructureOutput('GROMACS log',
            'a GROMACS log records energies, temperatures and performance, not coordinates: ' +
            'open the .gro file mdrun wrote (the -deffnm output, e.g. water.gro)');
    }
    if (isLammpsLog(content)) {
        return noStructureOutput('LAMMPS log',
            'a LAMMPS log records the run and its thermo output, not coordinates: ' +
            'open the data file or a dump file');
    }
    if (isDpgenModelDevi(content)) {
        return noStructureOutput('DP-GEN model deviation',
            'model_devi.out is a model-deviation table, not a structure: ' +
            'open the corresponding POSCAR/XYZ frames or the LAMMPS dump');
    }
    if (isGpumdInput(content)) {
        return noStructureOutput('GPUMD input',
            'a GPUMD run.in selects the potential and the files to dump, not a structure: ' +
            'open model.xyz or the dump.xyz it writes');
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
            const structured = structuredOutput(content, options);
            return structured || parseLogAsSingleFrame(content);
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
            const structured = structuredOutput(content, options);
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
            const guessed = tryAutoParse(content, fileName, options);
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
export function parseLogFile(content: string, fileName?: string): LogFileResult {
    if (isOrcaOutput(content)) {
        const result = parseOrcaOut(content);
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
        const xResult = parseXtbLog(content);
        return {
            frames: xResult.frames,
            title: xResult.title,
            optSteps: xResult.optSteps,
            normalModes: xResult.normalModes,
            source: 'xtb'
        };
    }
    const gResult = parseGaussianLog(content);
    return {
        frames: gResult.frames,
        title: gResult.title,
        optSteps: gResult.optSteps,
        normalModes: gResult.normalModes,
        routes: gResult.routes,
        source: 'gaussian'
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

    if (content.includes('$coord')) {
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

    const structured = structuredOutput(content, options);
    if (structured) return structured;

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
