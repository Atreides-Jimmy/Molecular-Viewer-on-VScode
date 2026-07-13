import { Atom, Bond, MolecularData, OptStep, NormalMode } from '../types';

const ATOMIC_NUMBER_MAP: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr'
};

export interface LogFrame {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    stepLabel: string;
    charge?: number;
    multiplicity?: number;
}

export interface GaussianLogResult {
    frames: LogFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
}

function skipDashedLines(lines: string[], startIdx: number, count: number): number {
    let i = startIdx;
    let skipped = 0;
    while (i < lines.length && skipped < count) {
        if (lines[i].includes('---')) {
            skipped++;
        }
        i++;
    }
    return i;
}

function parseOptStepsAndModes(content: string, lines: string[]): { optSteps?: OptStep[]; normalModes?: NormalMode[] } {
    const optSteps: OptStep[] = [];
    let pendingEnergy: number | undefined;
    let stepCounter = 0;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];

        // SCF Done:  E(RHF) =  -76.0107469158     A.U. after   10 cycles
        const scfMatch = line.match(/SCF Done:\s+E\([^)]+\)\s*=\s*(-?\d+\.\d+)/);
        if (scfMatch) {
            pendingEnergy = parseFloat(scfMatch[1]);
            continue;
        }

        // Convergence criteria block starts with "Maximum Force" data line
        if (line.includes('Maximum Force') && /-?\d+\.\d/.test(line)) {
            const step: OptStep = { step: stepCounter + 1, energy: pendingEnergy };
            stepCounter++;
            pendingEnergy = undefined;

            // Parse this block: 4 lines (Maximum Force, RMS Force, Maximum Displacement, RMS Displacement)
            for (let k = 0; k < 4 && i < lines.length; k++, i++) {
                const cur = lines[i];
                const valMatch = cur.match(/(-?\d+\.\d+(?:[EDed][-+]?\d+)?)/);
                if (!valMatch) continue;
                const val = parseFloat(valMatch[1].replace(/[EDed]/i, 'e'));
                if (cur.includes('Maximum Force') && !cur.includes('RMS')) {
                    step.maxForce = val;
                } else if (cur.includes('RMS') && cur.includes('Force')) {
                    step.rmsForce = val;
                } else if (cur.includes('Maximum Displacement')) {
                    step.maxDisplacement = val;
                } else if (cur.includes('RMS') && cur.includes('Displacement')) {
                    step.rmsDisplacement = val;
                }
            }
            i--; // compensate for loop increment
            optSteps.push(step);
            continue;
        }
    }

    const normalModes = parseNormalModes(lines);

    return {
        optSteps: optSteps.length > 0 ? optSteps : undefined,
        normalModes: normalModes && normalModes.length > 0 ? normalModes : undefined
    };
}

function parseNormalModes(lines: string[]): NormalMode[] | undefined {
    // Find the LAST "Harmonic frequencies (cm**-1)" header (CalcAll has multiple blocks)
    let headerIdx = -1;
    for (let i = 0; i < lines.length; i++) {
        if (lines[i].includes('Harmonic frequencies') && lines[i].includes('cm**-1')) {
            headerIdx = i;
        }
    }
    if (headerIdx < 0) return undefined;

    const modes: NormalMode[] = [];
    let i = headerIdx + 1;

    // Skip blank lines and column header text until we reach the first frequency block
    while (i < lines.length) {
        const line = lines[i].trim();

        // Detect end of normal modes section
        if (line === '' && modes.length > 0) {
            // Check if next non-blank line is another block or end
            let j = i;
            while (j < lines.length && lines[j].trim() === '') j++;
            if (j >= lines.length) break;
            const nextLine = lines[j].trim();
            // If next line looks like another block header (starts with numbers), continue
            if (!/^\d/.test(nextLine) && !nextLine.includes('Frequencies')) break;
        }

        // Block header line: "1                      2                      3" (1-3 columns)
        const colMatch = line.match(/^(\d+)(?:\s+(\d+))?(?:\s+(\d+))?\s*$/);
        if (!colMatch) {
            i++;
            continue;
        }

        // Next line: symmetry labels "A                      A                      A"
        const symLine = lines[i + 1] || '';
        const symParts = symLine.trim().split(/\s+/).filter(function(s){return s.length>0});

        // Frequencies line: "Frequencies --   1583.7822              3674.5814              3795.6092"
        const freqLine = lines[i + 2] || '';
        if (!freqLine.includes('Frequencies')) {
            i++;
            continue;
        }
        const freqParts = freqLine.split(/Frequencies\s+--/)[1] || '';
        const freqs = freqParts.trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));

        if (freqs.length === 0) { i++; continue; }

        // Reduced masses, force constants, IR intensities
        let redMasses: number[] = [];
        let frcConsts: number[] = [];
        let irIntens: number[] = [];
        let atomStartIdx = i + 3;
        for (let k = i + 3; k < Math.min(lines.length, i + 12); k++) {
            const kl = lines[k];
            if (kl.includes('Red. masses')) {
                const parts = kl.split(/Red\. masses\s+--/);
                if (parts[1]) redMasses = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('Frc consts')) {
                const parts = kl.split(/Frc consts\s+--/);
                if (parts[1]) frcConsts = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('IR Inten')) {
                const parts = kl.split(/IR Inten\s+--/);
                if (parts[1]) irIntens = parts[1].trim().split(/\s+/).map(parseFloat).filter(v => !isNaN(v));
            } else if (kl.includes('Atom') && kl.includes('AN')) {
                atomStartIdx = k + 1;
                break;
            }
        }

        // Parse displacement vectors until blank line or non-atom line
        const nCols = freqs.length;
        const displacementsByMode: number[][][] = Array.from({ length: nCols }, () => []);

        let rowIdx = atomStartIdx;
        while (rowIdx < lines.length) {
            const rl = lines[rowIdx].trim();
            if (rl === '' || rl.includes('---') || rl.includes('Frequencies') || /^(\d+)(?:\s+(\d+))?(?:\s+(\d+))?\s*$/.test(rl)) break;

            const parts = rl.split(/\s+/).map(parseFloat);
            // Format: atomIdx atomicNum X1 Y1 Z1 [X2 Y2 Z2] [X3 Y3 Z3]
            if (parts.length >= 2 + nCols * 3 && !parts.some(p => isNaN(p))) {
                for (let c = 0; c < nCols; c++) {
                    const base = 2 + c * 3;
                    displacementsByMode[c].push([parts[base], parts[base + 1], parts[base + 2]]);
                }
            }
            rowIdx++;
        }

        for (let c = 0; c < freqs.length; c++) {
            const globalIdx = modes.length + 1;
            modes.push({
                index: globalIdx,
                frequency: freqs[c],
                symmetry: symParts[c] || undefined,
                reducedMass: redMasses[c],
                forceConstant: frcConsts[c],
                irIntensity: irIntens[c],
                displacements: displacementsByMode[c] || []
            });
        }

        i = rowIdx;
    }

    return modes.length > 0 ? modes : undefined;
}

export function parseGaussianLog(content: string): GaussianLogResult {
    const lines = content.split(/\r?\n/);
    const frames: LogFrame[] = [];
    let title = '';

    const titleMatch = content.match(/#\s*[A-Za-z]/);
    if (titleMatch) {
        const routeIdx = content.indexOf(titleMatch[0]);
        const afterRoute = content.substring(0, routeIdx);
        const titleLines = afterRoute.split(/\r?\n/).filter(l => l.trim() !== '' && !l.trim().startsWith('%'));
        if (titleLines.length > 0) title = titleLines[titleLines.length - 1].trim();
    }

    const chargeMultMatch = content.match(/Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)/);
    let chargeMultLine = '';
    let logCharge: number | undefined;
    let logMultiplicity: number | undefined;
    if (chargeMultMatch) {
        chargeMultLine = `Charge=${chargeMultMatch[1]} Mult=${chargeMultMatch[2]}`;
        logCharge = parseInt(chargeMultMatch[1], 10);
        logMultiplicity = parseInt(chargeMultMatch[2], 10);
    }

    let optStep = 0;
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i];

        if (line.includes('Standard orientation:') || line.includes('Input orientation:')) {
            const isStandard = line.includes('Standard orientation:');

            i = skipDashedLines(lines, i + 1, 2);

            const atoms: Atom[] = [];
            while (i < lines.length) {
                const coordLine = lines[i].trim();
                if (coordLine === '' || coordLine.includes('---')) break;

                const parts = coordLine.split(/\s+/);
                if (parts.length >= 6) {
                    const atomicNum = parseInt(parts[1], 10);
                    if (isNaN(atomicNum)) { i++; continue; }
                    const element = ATOMIC_NUMBER_MAP[atomicNum] || 'X';
                    const x = parseFloat(parts[3]);
                    const y = parseFloat(parts[4]);
                    const z = parseFloat(parts[5]);

                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ element, x, y, z, index: atoms.length });
                    }
                }
                i++;
            }

            if (atoms.length > 0) {
                const label = isStandard
                    ? `Step ${optStep + 1}${chargeMultLine ? ' (' + chargeMultLine + ')' : ''}`
                    : `Input${chargeMultLine ? ' (' + chargeMultLine + ')' : ''}`;
                frames.push({
                    atoms,
                    bonds: [],
                    title: label,
                    hasExplicitBonds: false,
                    stepLabel: label,
                    charge: logCharge,
                    multiplicity: logMultiplicity
                });
                if (isStandard) optStep++;
            }
        }
    }

    if (frames.length === 0) {
        for (let i = 0; i < lines.length; i++) {
            if (lines[i].includes('Coordinates (Angstroms)')) {
                i = skipDashedLines(lines, i + 1, 2);

                const atoms: Atom[] = [];
                while (i < lines.length) {
                    const coordLine = lines[i].trim();
                    if (coordLine === '' || coordLine.includes('---')) break;

                    const parts = coordLine.split(/\s+/);
                    if (parts.length >= 4) {
                        const element = parts[1].replace(/[0-9]/g, '');
                        const el = element.charAt(0).toUpperCase() + element.slice(1).toLowerCase();
                        const x = parseFloat(parts[2]);
                        const y = parseFloat(parts[3]);
                        const z = parseFloat(parts[4]);

                        if (el && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                            atoms.push({ element: el, x, y, z, index: atoms.length });
                        }
                    }
                    i++;
                }

                if (atoms.length > 0) {
                    frames.push({
                        atoms,
                        bonds: [],
                        title: 'Coordinates',
                        hasExplicitBonds: false,
                        stepLabel: 'Coordinates',
                        charge: logCharge,
                        multiplicity: logMultiplicity
                    });
                }
                break;
            }
        }
    }

    const extra = parseOptStepsAndModes(content, lines);

    return {
        frames,
        title,
        optSteps: extra.optSteps,
        normalModes: extra.normalModes
    };
}
