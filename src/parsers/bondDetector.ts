import { Atom, Bond, MolecularData } from '../types';

const COVALENT_RADII: { [key: string]: number } = {
    H: 0.31, He: 0.28, Li: 1.28, Be: 0.96, B: 0.85, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Ne: 0.58,
    Na: 1.66, Mg: 1.41, Al: 1.21, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02, Ar: 1.06,
    K: 2.03, Ca: 1.76, Sc: 1.70, Ti: 1.60, V: 1.53, Cr: 1.39, Mn: 1.39, Fe: 1.32,
    Co: 1.26, Ni: 1.24, Cu: 1.32, Zn: 1.22, Ga: 1.22, Ge: 1.20, As: 1.19, Se: 1.20,
    Br: 1.20, Kr: 1.16, I: 1.39,
    Rb: 2.20, Sr: 1.95, Y: 1.90, Zr: 1.75, Nb: 1.64, Mo: 1.54, Tc: 1.47, Ru: 1.46,
    Rh: 1.42, Pd: 1.39, Ag: 1.45, Cd: 1.44, In: 1.42, Sn: 1.39, Sb: 1.39, Te: 1.38,
    Xe: 1.40,
    Cs: 2.44, Ba: 2.15, La: 2.07, Ce: 2.04, Pr: 2.03, Nd: 2.01, Pm: 1.99, Sm: 1.98,
    Eu: 1.98, Gd: 1.96, Tb: 1.94, Dy: 1.92, Ho: 1.92, Er: 1.89, Tm: 1.90, Yb: 1.87,
    Lu: 1.87, Hf: 1.75, Ta: 1.70, W: 1.62, Re: 1.51, Os: 1.44, Ir: 1.41, Pt: 1.36,
    Au: 1.36, Hg: 1.32, Tl: 1.45, Pb: 1.46, Bi: 1.48, Po: 1.40, At: 1.50, Rn: 1.50
};

interface BondSpec {
    order: number;
    length: number;
    tolerance: number;
}

function pairKey(el1: string, el2: string): string {
    const e1 = el1.toUpperCase();
    const e2 = el2.toUpperCase();
    return e1 < e2 ? e1 + e2 : e2 + e1;
}

function specKey(el1: string, el2: string): string {
    const e1 = el1.charAt(0).toUpperCase() + el1.slice(1).toLowerCase();
    const e2 = el2.charAt(0).toUpperCase() + el2.slice(1).toLowerCase();
    return e1 < e2 ? e1 + '+' + e2 : e2 + '+' + e1;
}

const BOND_SPECS: { [key: string]: BondSpec[] } = {
    'C+C': [
        { order: 3.0, length: 1.20, tolerance: 0.05 },
        { order: 1.5, length: 1.39, tolerance: 0.05 },
        { order: 2.0, length: 1.38, tolerance: 0.05 },
        { order: 1.0, length: 1.51, tolerance: 0.10 }
    ],
    'C+N': [
        { order: 2.0, length: 1.26, tolerance: 0.05 },
        { order: 1.5, length: 1.36, tolerance: 0.05 },
        { order: 1.0, length: 1.43, tolerance: 0.10 },
        { order: 3.0, length: 1.16, tolerance: 0.06 }
    ],
    'C+O': [
        { order: 2.0, length: 1.24, tolerance: 0.05 },
        { order: 1.0, length: 1.39, tolerance: 0.05 }
    ],
    'N+N': [
        { order: 1.0, length: 1.41, tolerance: 0.10 },
        { order: 2.0, length: 1.25, tolerance: 0.06 },
        { order: 3.0, length: 1.10, tolerance: 0.06 }
    ],
    'N+O': [
        { order: 2.0, length: 1.20, tolerance: 0.06 },
        { order: 1.5, length: 1.30, tolerance: 0.06 },
        { order: 1.0, length: 1.40, tolerance: 0.15 }
    ],
    'O+O': [
        { order: 2.0, length: 1.21, tolerance: 0.06 },
        { order: 1.0, length: 1.48, tolerance: 0.15 }
    ],
    'C+S': [
        { order: 1.5, length: 1.73, tolerance: 0.06 },
        { order: 2.0, length: 1.60, tolerance: 0.10 },
        { order: 1.0, length: 1.82, tolerance: 0.15 }
    ],
    'C+F': [
        { order: 1.0, length: 1.33, tolerance: 0.10 }
    ],
    'C+H': [{ order: 1.0, length: 0.97, tolerance: 0.15 }],
    'N+H': [{ order: 1.0, length: 0.88, tolerance: 0.15 }],
    'O+H': [{ order: 1.0, length: 0.85, tolerance: 0.15 }]
};

const BOND_CUTOFF: { [key: string]: number } = {
    'HH': 0.0, 'CH': 1.3, 'HO': 1.2, 'HN': 1.3,
    'CC': 1.9, 'CO': 1.7, 'CN': 1.7, 'NN': 1.7,
    'NO': 1.8, 'CF': 1.6, 'CS': 2.0
};

const MAX_VALENCE: { [key: string]: number } = {
    H: 1, C: 4, N: 3, O: 2, F: 1, S: 6, P: 5, Cl: 1, Br: 1, I: 1, B: 3
};

function distance(a1: Atom, a2: Atom): number {
    const dx = a1.x - a2.x;
    const dy = a1.y - a2.y;
    const dz = a1.z - a2.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function getBondOrder(el1: string, el2: string, dist: number): number {
    const pk = pairKey(el1, el2);
    const cutoff = BOND_CUTOFF[pk];
    if (cutoff !== undefined) {
        if (dist > cutoff) return 0;
    } else {
        const r1 = COVALENT_RADII[el1.charAt(0).toUpperCase() + el1.slice(1).toLowerCase()] || 1.5;
        const r2 = COVALENT_RADII[el2.charAt(0).toUpperCase() + el2.slice(1).toLowerCase()] || 1.5;
        if (dist > (r1 + r2) + 0.5) return 0;
    }

    const sk = specKey(el1, el2);
    const specs = BOND_SPECS[sk];
    if (specs) {
        for (const spec of specs) {
            if (Math.abs(dist - spec.length) <= spec.tolerance) {
                return spec.order;
            }
        }
        let bestOrder = 1.0;
        let minDiff = Infinity;
        for (const spec of specs) {
            const diff = Math.abs(dist - spec.length);
            if (diff < minDiff) {
                minDiff = diff;
                bestOrder = spec.order;
            }
        }
        return bestOrder;
    }

    const r1 = COVALENT_RADII[el1.charAt(0).toUpperCase() + el1.slice(1).toLowerCase()] || 1.5;
    const r2 = COVALENT_RADII[el2.charAt(0).toUpperCase() + el2.slice(1).toLowerCase()] || 1.5;
    const rSum = r1 + r2;
    const ratio = rSum !== 0 ? dist / rSum : 1.0;
    if (ratio < 0.85) return 3.0;
    if (ratio < 0.90) return 2.0;
    return 1.0;
}

function refineBondOrders(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[]
): void {
    const maxIter = 10;
    for (let iter = 0; iter < maxIter; iter++) {
        let changed = false;
        const valencies = new Map<number, number>();
        for (let i = 0; i < atoms.length; i++) valencies.set(i, 0);
        bondMap.forEach((neighbors, i) => {
            neighbors.forEach((order, _j) => {
                valencies.set(i, (valencies.get(i) || 0) + order);
            });
        });

        const violations: number[] = [];
        valencies.forEach((val, i) => {
            const el = atoms[i].element.charAt(0).toUpperCase() + atoms[i].element.slice(1).toLowerCase();
            const maxV = MAX_VALENCE[el] || 100;
            if (val > maxV + 0.1) violations.push(i);
        });

        if (violations.length === 0) break;

        for (const i of violations) {
            const neighbors = bondMap.get(i);
            if (!neighbors) continue;
            const currentValency = Array.from(neighbors.values()).reduce((s, v) => s + v, 0);
            const el = atoms[i].element.charAt(0).toUpperCase() + atoms[i].element.slice(1).toLowerCase();
            const maxV = MAX_VALENCE[el] || 100;
            if (currentValency <= maxV + 0.1) continue;

            let bestBond: [number, number] | null = null;
            let minLoss = Infinity;
            const el1 = atoms[i].element;

            neighbors.forEach((order, j) => {
                if (order <= 1.0) return;
                let nextOrder: number;
                if (order === 3.0) nextOrder = 2.0;
                else if (order === 2.0) nextOrder = 1.5;
                else if (order === 1.5) nextOrder = 1.0;
                else return;

                const el2 = atoms[j].element;
                const dist = distance(atoms[i], atoms[j]);
                const sk = specKey(el1, el2);
                const specs = BOND_SPECS[sk];
                if (!specs) return;

                let currentIdeal = 0, nextIdeal = 0;
                let minD = Infinity;
                for (const s of specs) {
                    if (s.order === order && Math.abs(dist - s.length) < minD) {
                        minD = Math.abs(dist - s.length);
                        currentIdeal = s.length;
                    }
                }
                minD = Infinity;
                for (const s of specs) {
                    if (s.order === nextOrder && Math.abs(dist - s.length) < minD) {
                        minD = Math.abs(dist - s.length);
                        nextIdeal = s.length;
                    }
                }
                if (currentIdeal === 0 || nextIdeal === 0) return;

                const loss = Math.abs(dist - nextIdeal) - Math.abs(dist - currentIdeal);
                if (loss < minLoss) {
                    minLoss = loss;
                    bestBond = [j, nextOrder];
                }
            });

            if (bestBond) {
                const bj = bestBond[0];
                const bOrder = bestBond[1];
                neighbors.set(bj, bOrder);
                bondMap.get(bj)?.set(i, bOrder);
                changed = true;
            }
        }

        if (!changed) break;
    }
}

function fixBondOrders(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[]
): void {
    bondMap.forEach((neighbors, i) => {
        neighbors.forEach((order, j) => {
            if (i > j) return;
            const el1 = atoms[i].element.toUpperCase();
            const el2 = atoms[j].element.toUpperCase();
            const elements = new Set([el1, el2]);

            if (elements.has('C') && elements.has('O')) {
                if (order !== 1.0 && order !== 2.0) {
                    const newOrder = order < 1.7 ? 1.0 : 2.0;
                    neighbors.set(j, newOrder);
                    bondMap.get(j)?.set(i, newOrder);
                }
            }

            if (el1 === 'BR' || el2 === 'BR') {
                neighbors.set(j, 1.0);
                bondMap.get(j)?.set(i, 1.0);
            }
        });
    });

    bondMap.forEach((neighbors, i) => {
        const el = atoms[i].element.toUpperCase();
        if (el !== 'N') return;
        const neighborList = Array.from(neighbors.entries());
        const numNeighbors = neighborList.length;

        if (numNeighbors === 3) {
            for (const [j, _order] of neighborList) {
                neighbors.set(j, 1.0);
                bondMap.get(j)?.set(i, 1.0);
            }
        } else if (numNeighbors === 2) {
            const [n1, bo1] = neighborList[0];
            const [n2, bo2] = neighborList[1];
            const n1IsH = atoms[n1].element.toUpperCase() === 'H';
            const n2IsH = atoms[n2].element.toUpperCase() === 'H';

            let f1: number, f2: number;
            if (n1IsH) {
                f1 = 1.0; f2 = 2.0;
            } else if (n2IsH) {
                f1 = 2.0; f2 = 1.0;
            } else {
                const d1 = Math.abs(bo1 - 1.5) + Math.abs(bo2 - 1.5);
                const d2 = Math.abs(bo1 - 2.0) + Math.abs(bo2 - 1.0);
                const d3 = Math.abs(bo1 - 1.0) + Math.abs(bo2 - 2.0);
                const best = Math.min(d1, d2, d3);
                if (best === d1) { f1 = 1.5; f2 = 1.5; }
                else if (best === d2) { f1 = 2.0; f2 = 1.0; }
                else { f1 = 1.0; f2 = 2.0; }
            }
            neighbors.set(n1, f1);
            bondMap.get(n1)?.set(i, f1);
            neighbors.set(n2, f2);
            bondMap.get(n2)?.set(i, f2);
        }
    });
}

export function detectBonds(atoms: Atom[]): Bond[] {
    const n = atoms.length;
    const bondMap = new Map<number, Map<number, number>>();

    for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
            const d = distance(atoms[i], atoms[j]);
            const bo = getBondOrder(atoms[i].element, atoms[j].element, d);
            if (bo > 0) {
                if (!bondMap.has(i)) bondMap.set(i, new Map());
                if (!bondMap.has(j)) bondMap.set(j, new Map());
                bondMap.get(i)!.set(j, bo);
                bondMap.get(j)!.set(i, bo);
            }
        }
    }

    fixBondOrders(bondMap, atoms);
    refineBondOrders(bondMap, atoms);

    const bonds: Bond[] = [];
    const seen = new Set<string>();
    bondMap.forEach((neighbors, i) => {
        neighbors.forEach((order, j) => {
            const key = Math.min(i, j) + '-' + Math.max(i, j);
            if (!seen.has(key)) {
                seen.add(key);
                bonds.push({ atom1: Math.min(i, j), atom2: Math.max(i, j), order });
            }
        });
    });

    return bonds;
}

export function ensureBonds(data: MolecularData): MolecularData {
    if (data.hasExplicitBonds && data.bonds.length > 0) {
        return data;
    }

    const bonds = detectBonds(data.atoms);
    return {
        ...data,
        bonds,
        hasExplicitBonds: false
    };
}
