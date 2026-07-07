import { Atom, CrystalData, MolecularData } from '../types';

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

function fracToCart(frac: number[], lattice: number[][]): [number, number, number] {
    const x = frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0];
    const y = frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1];
    const z = frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2];
    return [x, y, z];
}

export function parseVasp(content: string): MolecularData {
    const lines = content.split(/\r?\n/);

    if (lines.length < 8) {
        return { atoms: [], bonds: [], title: 'Invalid VASP file', hasExplicitBonds: false };
    }

    const title = lines[0].trim() || 'VASP Structure';

    const scale = parseFloat(lines[1].trim()) || 1.0;

    const lattice: number[][] = [];
    for (let i = 0; i < 3; i++) {
        const parts = lines[2 + i].trim().split(/\s+/).map(Number);
        if (parts.length < 3 || parts.some(isNaN)) {
            return { atoms: [], bonds: [], title: 'Invalid lattice vectors', hasExplicitBonds: false };
        }
        lattice.push([parts[0] * scale, parts[1] * scale, parts[2] * scale]);
    }

    let lineIdx = 5;

    let elementTypes: string[] = [];
    let elementCounts: number[] = [];

    const line6 = lines[lineIdx].trim();
    const line6Tokens = line6.split(/\s+/);

    const isElementTypeLine = line6Tokens.length > 0 &&
        line6Tokens.every(t => /^[A-Za-z][a-z]?$/.test(t)) &&
        !(line6Tokens.length === 1 && /^-?\d/.test(line6));

    if (isElementTypeLine) {
        elementTypes = line6Tokens.map(t => t.charAt(0).toUpperCase() + t.slice(1).toLowerCase());
        lineIdx++;

        const countLine = lines[lineIdx].trim();
        const countTokens = countLine.split(/\s+/).map(Number);
        if (countTokens.length === elementTypes.length && countTokens.every(n => !isNaN(n) && n >= 0)) {
            elementCounts = countTokens;
            lineIdx++;
        }
    }

    let coordLine = lines[lineIdx].trim();
    if (/^s/i.test(coordLine)) {
        lineIdx++;
        coordLine = lines[lineIdx].trim();
    }

    const isFractional = /^d/i.test(coordLine);
    const isCartesian = /^c/i.test(coordLine);
    lineIdx++;

    if (elementTypes.length === 0 || elementCounts.length === 0) {
        const totalAtoms = lines.length - lineIdx;
        elementTypes = ['X'];
        elementCounts = [totalAtoms];
    }

    const atoms: Atom[] = [];
    let atomIndex = 0;
    let elementTypeIdx = 0;
    let elementRemaining = elementCounts[0] || 0;

    for (let i = lineIdx; i < lines.length && elementTypeIdx < elementTypes.length; i++) {
        const line = lines[i].trim();
        if (line === '') continue;

        const parts = line.split(/\s+/);
        const x = parseFloat(parts[0]);
        const y = parseFloat(parts[1]);
        const z = parseFloat(parts[2]);

        if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

        if (elementRemaining <= 0) {
            elementTypeIdx++;
            if (elementTypeIdx >= elementTypes.length) break;
            elementRemaining = elementCounts[elementTypeIdx] || 0;
            if (elementRemaining <= 0) continue;
        }

        let cx: number, cy: number, cz: number;
        if (isCartesian) {
            cx = x * scale;
            cy = y * scale;
            cz = z * scale;
        } else {
            const wrapped = [x, y, z];
            const [fx, fy, fz] = fracToCart(wrapped, lattice);
            cx = fx;
            cy = fy;
            cz = fz;
        }

        atoms.push({
            element: elementTypes[elementTypeIdx],
            x: cx, y: cy, z: cz,
            index: atomIndex,
            occupancy: 1,
            baseIdx: atomIndex
        });
        atomIndex++;
        elementRemaining--;
    }

    const cellParams = latticeToCellParams(lattice);

    const crystal: CrystalData = {
        ...cellParams,
        latticeVectors: lattice,
        symmetryOps: ['x, y, z'],
        baseAtoms: atoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: []
    };

    return {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        crystal
    };
}
