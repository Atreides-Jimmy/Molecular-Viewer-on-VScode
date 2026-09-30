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

// Elements that carry exactly one bond and therefore have no bond order left
// to reassign. A valence fix must never touch the bond these atoms make: a
// hydrogen's single sigma bond IS its chemistry, so "lowering" it produces a
// hydrogen drawn with an empty valence shell instead of a plain single bond.
const TERMINAL_ELEMENTS = new Set([
    'H', 'F', 'Cl', 'Br', 'I',
    'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'
]);

function isTerminalElement(el: string): boolean {
    return TERMINAL_ELEMENTS.has(el.charAt(0).toUpperCase() + el.slice(1).toLowerCase());
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

// Above this many atoms the valence relaxation is skipped. It is an interface
// polish, not a correctness requirement: the bond table and the geometric
// fallback already produced a usable bond list. Structures this large are
// trajectories and protein-scale assemblies, where the extra pass costs more
// than the drawing is worth — a multi-hundred-thousand-atom log would otherwise
// stall the editor while the user waits for a picture.
const MAX_RELAX_ATOMS = 20000;

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

// Constraint-based bond-order assignment for a single atom.
//
// A valence ceiling alone is not enough: a multivalent atom (O, N) carrying one
// extra near neighbour — the common case is a hydroxyl oxygen that also accepts
// a hydrogen bond — starts one unit over its ceiling. The previous version
// treated that as "remove one unit of bond order" and walked every neighbour
// DOWN one step, which is lossy: nothing ever raised an order again, so the
// deficit stayed behind as an over-tight neighbour bond (a hydroxyl oxygen came
// out with a double bond to hydrogen).
//
// The correct question is not "which neighbour do I downgrade" but "which of
// the allowed order combinations puts every neighbour closest to the length it
// was actually measured at". Choosing that combination restores the neighbours
// whose real geometry says single bond (a normal O-H) and pushes the excess
// onto the neighbour whose geometry argues for it (the long contact). Both
// directions of change are available, so the result is stable: re-running it
// leaves the assignment exactly as it is.
//
// Local optimum over the neighbours of ONE centre. The search contributes no
// score of its own, only relative differences, which is what makes it
// reversible.
//
// The cost model is linear rather than a Cartesian sweep. A plan is described
// by how many units of order each free neighbour gains or loses relative to its
// current value; the centre's valence conservation makes the gains and losses
// balance, so the search is over small integer "give one unit to neighbour k"
// moves. Two phases cover it: single transfers, then a balanced pair swap. At
// most a handful of neighbours, so the candidate count stays tiny.
function relaxCenterValence(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[],
    i: number,
    entries: Array<[number, number]>,
    curDev: number,
    curSum: number
): boolean {
    const el = atoms[i].element;
    const maxV = MAX_VALENCE[el.charAt(0).toUpperCase() + el.slice(1).toLowerCase()] || 100;

    // Non-terminal neighbours whose order may be reassigned. A terminal
    // neighbour (hydrogen, halogen) stays a plain single bond: it has nowhere
    // to put extra order and nothing to give away.
    const free: number[] = [];
    let terminalSum = 0;
    for (let k = 0; k < entries.length; k++) {
        const j = entries[k][0];
        if (isTerminalElement(atoms[j].element)) terminalSum += 1.0;
        else free.push(k);
    }
    if (free.length === 0) return false;

    const nf = free.length;

    // Candidate orders per free neighbour, and the deviation each value implies
    // at the measured distance. Everything below is a lookup, so the sweep does
    // no distance maths.
    //
    // A candidate is only admitted when the measured length actually supports
    // it: either the distance lands inside that spec's tolerance, or the spec is
    // the pair's closest fit. Without that filter the search would happily
    // invent multiple bonds from arbitrary long contacts, because a stretching
    // bond's deviation grows monotonically and "any order is an improvement" is
    // false for a distance that is simply too long for a bond.
    const cand: number[][] = [];
    const devs: number[][] = [];
    const curIdx: number[] = [];
    for (let fi = 0; fi < nf; fi++) {
        const j = entries[free[fi]][0];
        const dist = distance(atoms[i], atoms[j]);
        const specs = BOND_SPECS[specKey(el, atoms[j].element)];
        const ladder: number[] = [];
        const add = (o: number) => {
            if (o > 0 && ladder.indexOf(o) < 0) ladder.push(o);
        };
        const devOf: number[] = [];

        if (specs) {
            // Closest-fitting order, then only those orders whose tolerance
            // actually brackets the measured distance.
            let closest = specs[0].order;
            let closestDiff = Infinity;
            for (const s of specs) {
                const diff = Math.abs(dist - s.length);
                if (diff < closestDiff) {
                    closestDiff = diff;
                    closest = s.order;
                }
            }
            add(closest);
            for (const s of specs) {
                if (Math.abs(dist - s.length) <= s.tolerance) add(s.order);
            }
            const cur = entries[free[fi]][1];
            if (ladder.indexOf(cur) < 0) add(cur);
            for (let k = 0; k < ladder.length; k++) {
                let ideal = dist;
                let minD = Infinity;
                for (const s of specs) {
                    if (s.order === ladder[k] && Math.abs(dist - s.length) < minD) {
                        minD = Math.abs(dist - s.length);
                        ideal = s.length;
                    }
                }
                devOf.push(Math.abs(dist - ideal));
            }
        } else {
            // No spec for the pair: the initial geometric fallback is all the
            // evidence there is, so it stays put.
            add(entries[free[fi]][1]);
            devOf.push(0);
        }
        cand.push(ladder);
        devs.push(devOf);
        const ci = ladder.indexOf(entries[free[fi]][1]);
        curIdx.push(ci >= 0 ? ci : 0);
    }

    // Room each neighbour has under its own ceiling, so a promotion is never
    // allowed to push a different atom over its limit. This is what keeps the
    // pass from trading one violation for another.
    const headroom: number[] = [];
    for (let fi = 0; fi < nf; fi++) {
        const j = entries[free[fi]][0];
        const el2 = atoms[j].element.charAt(0).toUpperCase() + atoms[j].element.slice(1).toLowerCase();
        const cap = MAX_VALENCE[el2] || 100;
        const neighbors = bondMap.get(j);
        let own = 0;
        if (neighbors) neighbors.forEach((o, k) => { if (k !== i) own += o; });
        headroom.push(cap - own);
    }

    const chosen = curIdx.slice();
    let bestSigned = 0;
    let bestAbs = 0;
    let found = false;

    // Change neighbour fi to candidate value i1 and fj to i2, then keep the plan
    // if the centre still fits its ceiling, no neighbour is pushed over its own,
    // and the neighbours sit closer to their measured lengths. A deficit at the
    // centre is rejected outright: that is what stops the pass from trading a
    // positive violation for a negative one.
    const consider = (fi: number, i1: number, fj: number, i2: number): void => {
        let sum = terminalSum;
        for (let k = 0; k < nf; k++) sum += (k === fi ? cand[k][i1] : k === fj ? cand[k][i2] : cand[k][curIdx[k]]);
        if (sum < curSum - 0.001) return;
        if (sum > maxV + 0.1) return;

        let signed = 0;
        let abs = 0;
        for (let k = 0; k < nf; k++) {
            const v = k === fi ? cand[k][i1] : k === fj ? cand[k][i2] : cand[k][curIdx[k]];
            const d = k === fi ? devs[k][i1] : k === fj ? devs[k][i2] : devs[k][curIdx[k]];
            const delta = v - cand[k][curIdx[k]];
            if (delta > headroom[k] + 0.001) return;
            signed += d - devs[k][curIdx[k]];
            abs += d;
        }
        // Prefer the plan the geometry argues for; break ties on the total
        // deviation so equal-improvement plans stay deterministic.
        if (!found || signed < bestSigned - 1e-9 || (Math.abs(signed - bestSigned) <= 1e-9 && abs < bestAbs - 1e-9)) {
            bestSigned = signed;
            bestAbs = abs;
            for (let k = 0; k < nf; k++) {
                chosen[k] = k === fi ? i1 : k === fj ? i2 : curIdx[k];
            }
            found = true;
        }
    };

    // Phase 1 — transfer one unit of order from one free neighbour to another.
    for (let fi = 0; fi < nf; fi++) {
        for (let i1 = 0; i1 < cand[fi].length; i1++) {
            const delta = cand[fi][i1] - cand[fi][curIdx[fi]];
            if (delta === 0) continue;
            for (let fj = 0; fj < nf; fj++) {
                if (fj === fi) continue;
                for (let i2 = 0; i2 < cand[fj].length; i2++) {
                    if (Math.abs((cand[fj][i2] - cand[fj][curIdx[fj]]) + delta) > 0.001) continue;
                    consider(fi, i1, fj, i2);
                }
            }
        }
    }

    // Phase 2 — two neighbours move together by the same amount (e.g. two
    // bonds both settling to 1.5), which single transfers cannot express.
    for (let fi = 0; fi < nf; fi++) {
        for (let i1 = 0; i1 < cand[fi].length; i1++) {
            const delta = cand[fi][i1] - cand[fi][curIdx[fi]];
            if (delta === 0) continue;
            for (let fj = fi + 1; fj < nf; fj++) {
                for (let i2 = 0; i2 < cand[fj].length; i2++) {
                    if (Math.abs((cand[fj][i2] - cand[fj][curIdx[fj]]) + delta) > 0.001) continue;
                    consider(fi, i1, fj, i2);
                }
            }
        }
    }

    // The unchanged plan, so a centre that cannot be improved reports false and
    // the outer loop terminates instead of spinning.
    //
    // Rejecting a plan that is not actually better than leaving the atom alone
    // is what makes the pass idempotent.
    if (!found || bestSigned > -1e-9) return false;

    let improved = false;
    for (let fi = 0; fi < nf; fi++) {
        const k = free[fi];
        const next = cand[fi][chosen[fi]];
        if (next === entries[k][1]) continue;
        const j = entries[k][0];
        bondMap.get(i)?.set(j, next);
        bondMap.get(j)?.set(i, next);
        improved = true;
    }
    return improved;
}

// Relax every atom whose valence drifts away from its ceiling. The pass runs to
// a fixed point: each accepted plan can move a neighbour away from its own
// ceiling, so the affected neighbourhoods are re-queued until nothing changes.
// The deviation terms are relative, so no absolute score is accumulated and no
// convergence tolerance has to be guessed.
//
// The iteration cap scales with the worst over-coordination in the structure
// rather than being a fixed constant: a molecule whose oxygen accepts three
// hydrogen bonds needs more rounds than a clean one, and capping at a flat
// number silently left those atoms unresolved. Real structures settle in one or
// two rounds; the bound only exists so a pathological input cannot loop.
function refineBondOrders(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[]
): void {
    let worstOver = 0;
    const valencies = new Map<number, number>();
    for (let i = 0; i < atoms.length; i++) valencies.set(i, 0);
    bondMap.forEach((neighbors, i) => {
        neighbors.forEach((order) => {
            valencies.set(i, (valencies.get(i) || 0) + order);
        });
    });
    valencies.forEach((val, i) => {
        const el = atoms[i].element.charAt(0).toUpperCase() + atoms[i].element.slice(1).toLowerCase();
        const maxV = MAX_VALENCE[el] || 100;
        const over = val - maxV;
        if (over > worstOver) worstOver = over;
    });

    const maxIter = Math.min(64, 4 + Math.ceil(worstOver) * 2);
    for (let iter = 0; iter < maxIter; iter++) {
        let changed = false;

        const current = new Map<number, number>();
        for (let i = 0; i < atoms.length; i++) current.set(i, 0);
        bondMap.forEach((neighbors, i) => {
            neighbors.forEach((order) => {
                current.set(i, (current.get(i) || 0) + order);
            });
        });

        const violations: number[] = [];
        current.forEach((val, i) => {
            const el = atoms[i].element.charAt(0).toUpperCase() + atoms[i].element.slice(1).toLowerCase();
            const maxV = MAX_VALENCE[el] || 100;
            if (val > maxV + 0.1) violations.push(i);
        });
        if (violations.length === 0) break;

        for (const i of violations) {
            const neighbors = bondMap.get(i);
            if (!neighbors) continue;
            const entries = Array.from(neighbors.entries());
            if (entries.length === 0) continue;

            let curSum = 0;
            let curDev = 0;
            for (const [j, order] of entries) {
                curSum += order;
                const dist = distance(atoms[i], atoms[j]);
                const specs = BOND_SPECS[specKey(atoms[i].element, atoms[j].element)];
                let ideal = dist;
                if (specs) {
                    let minD = Infinity;
                    for (const s of specs) {
                        if (s.order === order && Math.abs(dist - s.length) < minD) {
                            minD = Math.abs(dist - s.length);
                            ideal = s.length;
                        }
                    }
                }
                curDev += Math.abs(dist - ideal);
            }

            if (relaxCenterValence(bondMap, atoms, i, entries, curDev, curSum)) {
                changed = true;
            }
        }

        if (!changed) break;
    }
}

// Every spec order a given element pair is allowed to carry, read straight out
// of BOND_SPECS so the set can never drift away from the geometry table.
function allowedOrders(el1: string, el2: string): number[] {
    const specs = BOND_SPECS[specKey(el1, el2)];
    if (!specs) return [1.0];
    const orders: number[] = [];
    for (const s of specs) {
        if (orders.indexOf(s.order) < 0) orders.push(s.order);
    }
    return orders;
}

function nearestAllowedOrder(allowed: number[], order: number): number {
    let best = allowed[0];
    let bestDiff = Infinity;
    for (const a of allowed) {
        const diff = Math.abs(order - a);
        if (diff < bestDiff) {
            bestDiff = diff;
            best = a;
        }
    }
    return best;
}

// Snap every bond onto an order its element pair actually supports. This is a
// table-driven generalisation of the rule that used to be spelled out for a
// handful of pairs by hand: an order the table does not list for that pair is
// meaningless, so it collapses onto the closest listed one.
function fixBondOrders(
    bondMap: Map<number, Map<number, number>>,
    atoms: Atom[]
): void {
    bondMap.forEach((neighbors, i) => {
        neighbors.forEach((order, j) => {
            if (i > j) return;
            const allowed = allowedOrders(atoms[i].element, atoms[j].element);
            if (allowed.indexOf(order) >= 0) return;
            const snapped = nearestAllowedOrder(allowed, order);
            neighbors.set(j, snapped);
            bondMap.get(j)?.set(i, snapped);
        });
    });
}

// Distance above which two atoms of the given (normalized) elements can
// never bond: the pair-specific cutoff when defined, otherwise the covalent
// radii sum + 0.5 (mirrors getBondOrder's early-exit test exactly).
function pairCutoff(e1: string, e2: string): number {
    const c = BOND_CUTOFF[pairKey(e1, e2)];
    if (c !== undefined) return c;
    const r1 = COVALENT_RADII[e1] || 1.5;
    const r2 = COVALENT_RADII[e2] || 1.5;
    return r1 + r2 + 0.5;
}

export function detectBonds(atoms: Atom[]): Bond[] {
    const n = atoms.length;
    const bondMap = new Map<number, Map<number, number>>();
    if (n === 0) return [];

    // --- Spatial hash grid (3 Å cells) instead of the O(N²) pair loop ---
    // A bond requires d <= cutoff(element pair), so per element we compute
    // the maximum possible cutoff against every other element present; the
    // grid search radius floor(maxCutoff/CELL)+1 then safely covers every
    // pair that could bond. Candidate pairs are enumerated in exactly the
    // same (i asc, j asc) order as the full pair loop, so the bondMap
    // insertion sequence — and therefore bond order refinement results —
    // are identical to the previous implementation.
    const CELL = 3;
    const normElements = atoms.map(a => a.element.charAt(0).toUpperCase() + a.element.slice(1).toLowerCase());
    const distinct = Array.from(new Set(normElements));
    const rangeFor = new Map<string, number>();
    for (const e1 of distinct) {
        let maxC = 0;
        for (const e2 of distinct) {
            const c = pairCutoff(e1, e2);
            if (c > maxC) maxC = c;
        }
        rangeFor.set(e1, Math.floor(maxC / CELL) + 1);
    }

    const grid = new Map<number, number[]>();
    const cellX = new Int32Array(n), cellY = new Int32Array(n), cellZ = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        const cx = Math.floor(atoms[i].x / CELL);
        const cy = Math.floor(atoms[i].y / CELL);
        const cz = Math.floor(atoms[i].z / CELL);
        cellX[i] = cx; cellY[i] = cy; cellZ[i] = cz;
        // XOR hash may collide between different cells; collisions only merge
        // buckets (extra candidates filtered by the distance test), never lose one.
        const key = (cx * 73856093) ^ (cy * 19349663) ^ (cz * 83492791);
        let bucket = grid.get(key);
        if (!bucket) { bucket = []; grid.set(key, bucket); }
        bucket.push(i);
    }

    const candidates: number[] = [];
    for (let i = 0; i < n; i++) {
        const range = rangeFor.get(normElements[i])!;
        const cx = cellX[i], cy = cellY[i], cz = cellZ[i];
        candidates.length = 0;
        for (let dx = -range; dx <= range; dx++) {
            for (let dy = -range; dy <= range; dy++) {
                for (let dz = -range; dz <= range; dz++) {
                    const key = ((cx + dx) * 73856093) ^ ((cy + dy) * 19349663) ^ ((cz + dz) * 83492791);
                    const bucket = grid.get(key);
                    if (!bucket) continue;
                    for (let k = 0; k < bucket.length; k++) {
                        const j = bucket[k];
                        if (j > i) candidates.push(j);
                    }
                }
            }
        }
        if (candidates.length === 0) continue;
        candidates.sort((a, b) => a - b);
        let prev = -1;
        for (let ci = 0; ci < candidates.length; ci++) {
            const j = candidates[ci];
            if (j === prev) continue; // hash-collision duplicate
            prev = j;
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
    if (n <= MAX_RELAX_ATOMS) refineBondOrders(bondMap, atoms);

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
