"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseGaussianLog = parseGaussianLog;
const ATOMIC_NUMBER_MAP = {
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
function skipDashedLines(lines, startIdx, count) {
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
function parseGaussianLog(content) {
    const lines = content.split(/\r?\n/);
    const frames = [];
    let title = '';
    const titleMatch = content.match(/#\s*[A-Za-z]/);
    if (titleMatch) {
        const routeIdx = content.indexOf(titleMatch[0]);
        const afterRoute = content.substring(0, routeIdx);
        const titleLines = afterRoute.split(/\r?\n/).filter(l => l.trim() !== '' && !l.trim().startsWith('%'));
        if (titleLines.length > 0)
            title = titleLines[titleLines.length - 1].trim();
    }
    const chargeMultMatch = content.match(/Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)/);
    let chargeMultLine = '';
    let logCharge;
    let logMultiplicity;
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
            const atoms = [];
            while (i < lines.length) {
                const coordLine = lines[i].trim();
                if (coordLine === '' || coordLine.includes('---'))
                    break;
                const parts = coordLine.split(/\s+/);
                if (parts.length >= 6) {
                    const atomicNum = parseInt(parts[1], 10);
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
                if (isStandard)
                    optStep++;
            }
        }
    }
    if (frames.length === 0) {
        for (let i = 0; i < lines.length; i++) {
            if (lines[i].includes('Coordinates (Angstroms)')) {
                i = skipDashedLines(lines, i + 1, 2);
                const atoms = [];
                while (i < lines.length) {
                    const coordLine = lines[i].trim();
                    if (coordLine === '' || coordLine.includes('---'))
                        break;
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
    return { frames, title };
}
//# sourceMappingURL=logParser.js.map