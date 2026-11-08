import { Atom, CrystalData, MolecularData } from '../types';

const BOHR_TO_ANGSTROM = 0.52917721067;

const ELEMENT_BY_NUMBER: { [key: number]: string } = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf'
};

function isValidNumber(v: number): boolean {
    return typeof v === 'number' && isFinite(v) && !isNaN(v);
}

function latticeToCellParams(lattice: number[][]): { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } {
    const [aVec, bVec, cVec] = lattice;
    const a = Math.sqrt(aVec[0] * aVec[0] + aVec[1] * aVec[1] + aVec[2] * aVec[2]);
    const b = Math.sqrt(bVec[0] * bVec[0] + bVec[1] * bVec[1] + bVec[2] * bVec[2]);
    const c = Math.sqrt(cVec[0] * cVec[0] + cVec[1] * cVec[1] + cVec[2] * cVec[2]);
    const cosA = (bVec[0] * cVec[0] + bVec[1] * cVec[1] + bVec[2] * cVec[2]) / (b * c);
    const cosB = (aVec[0] * cVec[0] + aVec[1] * cVec[1] + aVec[2] * cVec[2]) / (a * c);
    const cosG = (aVec[0] * bVec[0] + aVec[1] * bVec[1] + aVec[2] * bVec[2]) / (a * b);
    const rad = 180 / Math.PI;
    return {
        a, b, c,
        alpha: Math.acos(Math.max(-1, Math.min(1, cosA))) * rad,
        beta: Math.acos(Math.max(-1, Math.min(1, cosB))) * rad,
        gamma: Math.acos(Math.max(-1, Math.min(1, cosG))) * rad
    };
}

export function parseCube(content: string): MolecularData {
    const lines = content.split(/\r?\n/);

    if (lines.length < 6) {
        return { atoms: [], bonds: [], title: 'Invalid cube file', hasExplicitBonds: false };
    }

    const title = lines[0].trim() || lines[1].trim() || 'Cube Structure';
    const warnings: string[] = [];

    const parseLine3 = lines[2].trim().split(/\s+/).map(Number);
    // A negative NATOMS is part of the cube convention, not an error: it means
    // the coordinates are Angstrom and one extra "orbital records" line follows
    // the atoms. Reading it as a count would read nothing at all.
    const declaredNatoms = parseLine3[0];
    const headerValid = parseLine3.length >= 4 &&
        isValidNumber(declaredNatoms) && Math.trunc(declaredNatoms) === declaredNatoms &&
        parseLine3.slice(1, 4).every(isValidNumber);
    if (!headerValid) {
        return {
            atoms: [], bonds: [], title: 'Invalid cube file', hasExplicitBonds: false,
            warnings: ['the cube header line 3 is unreadable (atom count or origin): no structure to show']
        };
    }
    const natoms = Math.abs(declaredNatoms);
    if (natoms === 0) {
        return {
            atoms: [], bonds: [], title, hasExplicitBonds: false,
            warnings: ['the cube declares zero atoms: no structure to show']
        };
    }
    const origin = [parseLine3[1], parseLine3[2], parseLine3[3]];

    const voxelLines: number[][] = [];
    const voxelCounts: number[] = [];
    for (let i = 0; i < 3; i++) {
        const parts = lines[3 + i].trim().split(/\s+/).map(Number);
        voxelCounts.push(parts[0]);
        voxelLines.push([parts[1] || 0, parts[2] || 0, parts[3] || 0]);
    }

    const voxelValid = voxelLines.every(v => v.every(isValidNumber)) &&
        voxelLines.some(v => v.some(x => Math.abs(x) > 1e-10));

    // Cube convention: NATOMS > 0 means Bohr, NATOMS < 0 means Angstrom
    let isBohr = declaredNatoms > 0;
    if (voxelValid) {
        const maxVoxel = Math.max(...voxelLines.flat().map(Math.abs));
        const maxCoordAbs = 50;
        if (maxVoxel > maxCoordAbs * 0.5) isBohr = false;
    }

    const lengthFactor = isBohr ? BOHR_TO_ANGSTROM : 1.0;

    let lattice: number[][] | null = null;
    if (voxelValid) {
        lattice = voxelLines.map((v, i) => [
            v[0] * voxelCounts[i] * lengthFactor,
            v[1] * voxelCounts[i] * lengthFactor,
            v[2] * voxelCounts[i] * lengthFactor
        ]);
    }

    const atoms: Atom[] = [];
    let atomIndex = 0;
    const atomStartLine = 6;
    let rowsRead = 0;
    let atomBlockEnd = atomStartLine;

    for (let i = atomStartLine; i < lines.length && rowsRead < natoms; i++) {
        const line = lines[i].trim();
        if (line === '') {
            atomBlockEnd = i + 1;
            continue;
        }

        const parts = line.split(/\s+/);
        // An atom row is `Z charge x y z`; a grid-data row is a list of
        // floating-point values. Requiring an integer Z keeps a truncated atom
        // table from turning grid data into ghost atoms.
        const atomicNumber = Number(parts[0]);
        const isAtomRow = parts.length === 5 &&
            Number.isInteger(atomicNumber) && Math.abs(atomicNumber) <= 118 &&
            parts.slice(1).every(p => isValidNumber(Number(p)));
        if (!isAtomRow) break; // grid data starts here: the atom table was cut short

        rowsRead++;
        atomBlockEnd = i + 1;
        if (atomicNumber <= 0) continue;

        const x = Number(parts[2]) * lengthFactor;
        const y = Number(parts[3]) * lengthFactor;
        const z = Number(parts[4]) * lengthFactor;

        const element = ELEMENT_BY_NUMBER[atomicNumber] || 'X';

        atoms.push({
            element,
            x, y, z,
            index: atomIndex,
            occupancy: 1,
            baseIdx: atomIndex
        });
        atomIndex++;
    }

    if (rowsRead < natoms) {
        warnings.push('the cube declares ' + natoms + ' atom' + (natoms === 1 ? '' : 's') +
            ' but only ' + rowsRead + ' atom row' + (rowsRead === 1 ? ' was' : 's were') +
            ' readable: the file may be truncated');
    }

    if (declaredNatoms < 0 && rowsRead === natoms) {
        // The line after the atom block names how many data grids follow. Only
        // its count is used, and only for the warning; the grids are not read.
        // It is the first non-empty line, not the first numeric line anywhere
        // below the atoms: a truncated file has grid values there instead.
        let records = 1;
        for (let i = atomBlockEnd; i < lines.length; i++) {
            const parts = lines[i].trim().split(/\s+/).map(Number);
            if (parts.length === 0 || !isValidNumber(parts[0])) continue;
            if (Math.trunc(parts[0]) === parts[0] && parts[0] > 0) records = parts[0];
            break;
        }
        warnings.push(`${records} orbital record${records === 1 ? '' : 's'} present - geometry only`);
    }

    if (!lattice && atoms.length > 0) {
        let minX = Infinity, maxX = -Infinity;
        let minY = Infinity, maxY = -Infinity;
        let minZ = Infinity, maxZ = -Infinity;
        atoms.forEach(a => {
            minX = Math.min(minX, a.x); maxX = Math.max(maxX, a.x);
            minY = Math.min(minY, a.y); maxY = Math.max(maxY, a.y);
            minZ = Math.min(minZ, a.z); maxZ = Math.max(maxZ, a.z);
        });
        const pad = 5.0;
        const sx = (maxX - minX) + 2 * pad;
        const sy = (maxY - minY) + 2 * pad;
        const sz = (maxZ - minZ) + 2 * pad;
        const ox = minX - pad;
        const oy = minY - pad;
        const oz = minZ - pad;
        atoms.forEach(a => {
            a.x -= ox; a.y -= oy; a.z -= oz;
        });
        lattice = [[sx, 0, 0], [0, sy, 0], [0, 0, sz]];
    }

    const result: MolecularData = {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false
    };
    if (warnings.length > 0) result.warnings = warnings;

    if (lattice) {
        const cellParams = latticeToCellParams(lattice);
        const crystal: CrystalData = {
            ...cellParams,
            latticeVectors: lattice,
            symmetryOps: ['x, y, z'],
            baseAtoms: atoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
            baseBonds: []
        };
        result.crystal = crystal;
    }

    return result;
}
