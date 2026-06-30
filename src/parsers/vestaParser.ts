import { Atom, CrystalData, MolecularData } from '../types';

function cellToLatticeVectors(a: number, b: number, c: number, alpha: number, beta: number, gamma: number): number[][] {
    const rad = Math.PI / 180;
    const ar = alpha * rad;
    const br = beta * rad;
    const gr = gamma * rad;

    const cosA = Math.cos(ar);
    const cosB = Math.cos(br);
    const cosG = Math.cos(gr);
    const sinG = Math.sin(gr);

    const aVec = [a, 0, 0];
    const bVec = [b * cosG, b * sinG, 0];
    const cX = c * cosB;
    const cY = c * (cosA - cosB * cosG) / sinG;
    const cZ2 = 1 - cosA * cosA - cosB * cosB - cosG * cosG + 2 * cosA * cosB * cosG;
    const cZ = c * Math.sqrt(Math.max(0, cZ2)) / sinG;
    const cVec = [cX, cY, cZ];

    return [aVec, bVec, cVec];
}

function fracToCart(frac: number[], lattice: number[][]): [number, number, number] {
    const x = frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0];
    const y = frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1];
    const z = frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2];
    return [x, y, z];
}

function wrapFrac(v: number): number {
    return v - Math.floor(v);
}

interface VestaAtom {
    element: string;
    occupancy: number;
    frac: [number, number, number];
}

const SECTION_KEYWORDS = /^(GROUP|SYMOP|TRANM|LTRANSL|LORIENT|LMATRIX|CELLP|STRUC|THERI|ATOMM|BONDM|POLYM|SURFM|FORMM|SECCL|TEXCL|LIGHT|MAPM|VECTM|SBOND|SANG|BDSPL|BDTYP|CMBOND|ORBIT|NETCHG|AATOMM|SATOMM)\b/;

export function parseVesta(content: string): MolecularData {
    const lines = content.split(/\r?\n/);

    let title = 'VESTA Structure';
    let spaceGroup: string | undefined;
    let cellParams: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } | null = null;
    const vestaAtoms: VestaAtom[] = [];

    let i = 0;
    while (i < lines.length) {
        const line = lines[i].trim();

        if (line === 'TITLE') {
            i++;
            const titleParts: string[] = [];
            while (i < lines.length && lines[i].trim() !== '' && !SECTION_KEYWORDS.test(lines[i].trim())) {
                titleParts.push(lines[i].trim());
                i++;
            }
            if (titleParts.length > 0) title = titleParts.join(' ');
            continue;
        }

        if (line === 'GROUP') {
            i++;
            if (i < lines.length) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 3) {
                    spaceGroup = parts.slice(2).join(' ');
                }
            }
            i++;
            continue;
        }

        if (line === 'CELLP') {
            i++;
            if (i < lines.length) {
                const parts = lines[i].trim().split(/\s+/).map(Number);
                if (parts.length >= 6 && parts.every(v => isFinite(v) && !isNaN(v))) {
                    cellParams = {
                        a: parts[0], b: parts[1], c: parts[2],
                        alpha: parts[3], beta: parts[4], gamma: parts[5]
                    };
                }
            }
            i++;
            continue;
        }

        if (line === 'STRUC') {
            i++;
            while (i < lines.length) {
                const strucLine = lines[i].trim();
                if (strucLine === '' || SECTION_KEYWORDS.test(strucLine)) {
                    break;
                }
                const parts = strucLine.split(/\s+/);
                if (parts.length >= 7) {
                    const element = parts[1];
                    const occupancy = parseFloat(parts[3]);
                    const fx = parseFloat(parts[4]);
                    const fy = parseFloat(parts[5]);
                    const fz = parseFloat(parts[6]);

                    if (/^[A-Za-z][A-Za-z]?$/.test(element) && isFinite(fx) && isFinite(fy) && isFinite(fz)) {
                        vestaAtoms.push({
                            element,
                            occupancy: isFinite(occupancy) ? occupancy : 1,
                            frac: [fx, fy, fz]
                        });
                        i += 2;
                        continue;
                    }
                }
                break;
            }
            continue;
        }

        if (line === 'SYMOP') {
            i++;
            while (i < lines.length) {
                const parts = lines[i].trim().split(/\s+/).map(Number);
                if (parts.length >= 3 && parts[0] === -1 && parts[1] === -1 && parts[2] === -1) break;
                i++;
            }
            i++;
            continue;
        }

        i++;
    }

    if (!cellParams) {
        const atoms: Atom[] = vestaAtoms.map((va, idx) => ({
            element: va.element,
            x: va.frac[0], y: va.frac[1], z: va.frac[2],
            index: idx, occupancy: va.occupancy, baseIdx: idx
        }));
        return { atoms, bonds: [], title, hasExplicitBonds: false };
    }

    const lattice = cellToLatticeVectors(cellParams.a, cellParams.b, cellParams.c, cellParams.alpha, cellParams.beta, cellParams.gamma);

    const cartAtoms: Atom[] = vestaAtoms.map((va, idx) => {
        const wf = [wrapFrac(va.frac[0]), wrapFrac(va.frac[1]), wrapFrac(va.frac[2])];
        const [cx, cy, cz] = fracToCart(wf, lattice);
        return {
            element: va.element,
            x: cx, y: cy, z: cz,
            index: idx,
            occupancy: va.occupancy,
            baseIdx: idx
        };
    });

    const crystal: CrystalData = {
        ...cellParams,
        latticeVectors: lattice,
        spaceGroup,
        symmetryOps: ['x, y, z'],
        baseAtoms: cartAtoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: []
    };

    return {
        atoms: cartAtoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        crystal
    };
}
