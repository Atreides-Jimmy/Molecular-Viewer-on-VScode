"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseGjf = parseGjf;
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
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
};
function resolveElement(raw) {
    let cleaned = raw.replace(/\(.*?\)/g, '');
    const asNum = parseInt(cleaned, 10);
    if (!isNaN(asNum) && asNum > 0 && ATOMIC_NUMBER_MAP[asNum]) {
        return ATOMIC_NUMBER_MAP[asNum];
    }
    cleaned = cleaned.replace(/[0-9]/g, '');
    const capitalized = cleaned.charAt(0).toUpperCase() + cleaned.slice(1).toLowerCase();
    if (capitalized.length >= 1 && capitalized.length <= 2 && /[A-Z][a-z]?/.test(capitalized)) {
        return capitalized;
    }
    return cleaned;
}
function isFixedAtomMarker(val) {
    return val === -1 || val === 0 || val === 1;
}
function parseAtomLine(parts) {
    if (parts.length < 4)
        return null;
    const element = resolveElement(parts[0]);
    if (parts.length === 4) {
        const x = parseFloat(parts[1]);
        const y = parseFloat(parts[2]);
        const z = parseFloat(parts[3]);
        if (!isNaN(x) && !isNaN(y) && !isNaN(z))
            return { element, x, y, z };
        return null;
    }
    const coordVals = [];
    const markerIndices = [];
    for (let pi = 1; pi < parts.length; pi++) {
        const v = parseFloat(parts[pi]);
        if (isNaN(v))
            continue;
        if (isFixedAtomMarker(v) && Number.isInteger(v)) {
            markerIndices.push(coordVals.length);
        }
        coordVals.push(v);
    }
    if (coordVals.length === 3) {
        return { element, x: coordVals[0], y: coordVals[1], z: coordVals[2] };
    }
    if (coordVals.length === 4 && markerIndices.length === 1) {
        const coords = coordVals.filter((_, idx) => idx !== markerIndices[0]);
        if (coords.length === 3)
            return { element, x: coords[0], y: coords[1], z: coords[2] };
    }
    if (coordVals.length >= 3) {
        const floatVals = coordVals.filter(v => !isFixedAtomMarker(v) || !Number.isInteger(v));
        if (floatVals.length >= 3) {
            return { element, x: floatVals[floatVals.length - 3], y: floatVals[floatVals.length - 2], z: floatVals[floatVals.length - 1] };
        }
    }
    return null;
}
function parseGjf(content) {
    const lines = content.split(/\r?\n/);
    const atoms = [];
    const bonds = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge;
    let multiplicity;
    const link0Lines = [];
    let routeLine = '';
    const titleLines = [];
    let chargeMultLine = '';
    let afterConnectContent = '';
    let i = 0;
    while (i < lines.length && lines[i].trim() !== '') {
        if (lines[i].trim().startsWith('%')) {
            link0Lines.push(lines[i]);
        }
        else {
            routeLine = lines[i].trim();
        }
        i++;
    }
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }
    while (i < lines.length && lines[i].trim() !== '') {
        titleLines.push(lines[i]);
        i++;
    }
    title = titleLines.map(l => l.trim()).join(' ');
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }
    if (i < lines.length) {
        chargeMultLine = lines[i];
        const cmParts = lines[i].trim().split(/\s+/);
        const chrg = parseInt(cmParts[0], 10);
        const mult = parseInt(cmParts[1], 10);
        if (!isNaN(chrg))
            charge = chrg;
        if (!isNaN(mult))
            multiplicity = mult;
        i++;
    }
    let atomIndex = 0;
    while (i < lines.length) {
        const line = lines[i].trim();
        if (line === '' || line.startsWith('--')) {
            break;
        }
        const parts = line.split(/\s+/);
        const result = parseAtomLine(parts);
        if (result) {
            atoms.push({ element: result.element, x: result.x, y: result.y, z: result.z, index: atomIndex });
            atomIndex++;
        }
        i++;
    }
    let connectStartLine = -1;
    while (i < lines.length && lines[i].trim() === '') {
        i++;
    }
    if (i < lines.length) {
        connectStartLine = i;
        let allNumeric = true;
        let lineCount = 0;
        let maxAtomNum = 0;
        for (let li = i; li < lines.length; li++) {
            const tl = lines[li].trim();
            if (tl === '' || tl.startsWith('--'))
                break;
            lineCount++;
            const tp = tl.split(/\s+/);
            const firstNum = parseInt(tp[0], 10);
            if (isNaN(firstNum)) {
                allNumeric = false;
                break;
            }
            if (firstNum > maxAtomNum)
                maxAtomNum = firstNum;
            for (let ti = 1; ti < tp.length; ti++) {
                if (isNaN(parseFloat(tp[ti]))) {
                    allNumeric = false;
                    break;
                }
            }
            if (!allNumeric)
                break;
        }
        if (allNumeric && lineCount > 0 && maxAtomNum <= atoms.length) {
            hasExplicitBonds = true;
            let connectEndLine = i;
            for (let li = i; li < lines.length; li++) {
                const tl = lines[li].trim();
                if (tl === '' || tl.startsWith('--')) {
                    connectEndLine = li;
                    break;
                }
                const cparts = tl.split(/\s+/);
                const atom1Num = parseInt(cparts[0], 10);
                if (atom1Num > atoms.length) {
                    connectEndLine = li;
                    break;
                }
                parseConnectLine(tl, atoms.length, bonds);
                connectEndLine = li + 1;
            }
            let afterIdx = connectEndLine;
            while (afterIdx < lines.length && lines[afterIdx].trim() === '') {
                afterIdx++;
            }
            if (afterIdx < lines.length) {
                afterConnectContent = lines.slice(afterIdx).join('\n');
            }
        }
        else {
            afterConnectContent = lines.slice(i).join('\n');
        }
    }
    else if (i < lines.length) {
        afterConnectContent = lines.slice(i).join('\n');
    }
    const gjfMeta = {
        link0Lines,
        routeLine,
        titleLines,
        chargeMultLine,
        afterConnectContent
    };
    return { atoms, bonds, title, hasExplicitBonds, gjfMeta, charge, multiplicity };
}
function parseConnectLine(line, totalAtoms, bonds) {
    const cparts = line.split(/\s+/);
    if (cparts.length < 2)
        return;
    const atom1 = parseInt(cparts[0], 10) - 1;
    if (isNaN(atom1) || atom1 < 0 || atom1 >= totalAtoms)
        return;
    let j = 1;
    while (j + 1 < cparts.length) {
        const atom2 = parseInt(cparts[j], 10) - 1;
        const bondOrder = parseFloat(cparts[j + 1]) || 1;
        if (!isNaN(atom2) && atom2 >= 0 && atom2 < totalAtoms && atom1 !== atom2) {
            const exists = bonds.some(b => (b.atom1 === atom1 && b.atom2 === atom2) ||
                (b.atom1 === atom2 && b.atom2 === atom1));
            if (!exists) {
                bonds.push({ atom1, atom2, order: bondOrder });
            }
        }
        j += 2;
    }
    if (cparts.length === 2) {
        const atom2 = parseInt(cparts[1], 10) - 1;
        if (!isNaN(atom2) && atom2 >= 0 && atom2 < totalAtoms && atom1 !== atom2) {
            const exists = bonds.some(b => (b.atom1 === atom1 && b.atom2 === atom2) ||
                (b.atom1 === atom2 && b.atom2 === atom1));
            if (!exists) {
                bonds.push({ atom1, atom2, order: 1 });
            }
        }
    }
}
//# sourceMappingURL=gjfParser.js.map