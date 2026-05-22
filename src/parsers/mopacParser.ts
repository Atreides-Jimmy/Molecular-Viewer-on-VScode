import { Atom, Bond, MolecularData } from '../types';

const ATOMIC_NUMBER_TO_ELEMENT: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn'
};

function parseMopacElement(token: string): string {
    const num = parseInt(token, 10);
    if (!isNaN(num) && num > 0 && ATOMIC_NUMBER_TO_ELEMENT[num]) {
        return ATOMIC_NUMBER_TO_ELEMENT[num];
    }
    const elem = token.charAt(0).toUpperCase() + token.slice(1).toLowerCase();
    if (/^[A-Za-z]{1,2}$/.test(elem)) return elem;
    return 'C';
}

export function parseMopac(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let hasExplicitBonds = false;
    let title = '';
    let charge: number | undefined;
    let multiplicity: number | undefined;
    let atomIndex = 0;

    let inCoords = false;
    let coordLines: string[] = [];
    let foundChargeMult = false;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        if (i === 0) {
            title = line || 'MOPAC Calculation';
            continue;
        }

        if (i === 1) continue;

        if (!foundChargeMult) {
            const cmMatch = line.match(/CHARGE\s*=\s*(-?\d+)/i);
            if (cmMatch) charge = parseInt(cmMatch[1], 10);

            const msMatch = line.match(/(?:MULT|MS|MULTIPlicity)\s*=\s*(\d+)/i);
            if (msMatch) multiplicity = parseInt(msMatch[1], 10);

            if (cmMatch || msMatch) {
                foundChargeMult = true;
            }

            if (line.match(/^\s*[A-Za-z]{1,2}\s+/) || line.match(/^\s*\d+\s+/)) {
                inCoords = true;
            }
        }

        if (inCoords) {
            if (line === '' || line.toLowerCase().startsWith('old') || line.toLowerCase().startsWith('===')) {
                if (coordLines.length > 0) break;
                continue;
            }
            coordLines.push(line);
        }
    }

    for (const cl of coordLines) {
        const parts = cl.split(/\s+/);
        if (parts.length < 4) continue;

        const elem = parseMopacElement(parts[0]);
        const xVal = parseFloat(parts[1]);
        const yVal = parseFloat(parts[3]);
        const zVal = parseFloat(parts[5]);

        if (isNaN(xVal) || isNaN(yVal) || isNaN(zVal)) continue;

        let x: number, y: number, z: number;
        const xLabel = parts.length > 2 ? parts[2] : '';
        const yLabel = parts.length > 4 ? parts[4] : '';
        const zLabel = parts.length > 6 ? parts[6] : '';

        const isAngstrom = xLabel.toUpperCase() === '1' || xLabel.toUpperCase() === 'A' || xLabel === 'Å';
        const isBohr = xLabel.toUpperCase() === '0' || xLabel.toUpperCase() === 'B';

        if (isAngstrom) {
            x = xVal; y = yVal; z = zVal;
        } else if (isBohr) {
            x = xVal * 0.529177249;
            y = yVal * 0.529177249;
            z = zVal * 0.529177249;
        } else {
            x = xVal; y = yVal; z = zVal;
        }

        atoms.push({ element: elem, x, y, z, index: atomIndex });
        atomIndex++;
    }

    if (atoms.length === 0) {
        inCoords = false;
        coordLines = [];
        let pastHeader = false;

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();

            if (line.match(/CARTESIAN COORDINATES/i) || line.match(/ATOM\s+X\s+Y\s+Z/i) || line.match(/FINAL GEOMETRY/i)) {
                pastHeader = true;
                continue;
            }

            if (pastHeader) {
                if (line === '' || line.startsWith('-') || line.startsWith('=')) {
                    if (coordLines.length > 0) break;
                    continue;
                }
                const parts = line.split(/\s+/);
                if (parts.length >= 4) {
                    const elem = parseMopacElement(parts[0]);
                    const x = parseFloat(parts[1]);
                    const y = parseFloat(parts[2]);
                    const z = parseFloat(parts[3]);
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ element: elem, x, y, z, index: atomIndex });
                        atomIndex++;
                    }
                }
            }
        }
    }

    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
