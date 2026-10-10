/**
 * Hydrogen completion ("Add H")
 * =============================
 * Fills the selected atoms whose bond valence falls short of their element's
 * empirical valence with hydrogens, and does nothing else: the structure is
 * read, never written, so every atom that was already there — and every bond
 * between them — stays exactly where it was. The new X–H bonds take their
 * length, and the new H–X–Y angles their value, from the SMILES -> 3D engine's
 * own experience tables (`empiricalBondLength` / `empiricalBondAngle` /
 * `empiricalValence`), so an added bond is the same bond a freshly generated
 * structure would carry instead of a second, drifting set of numbers:
 *
 *   - the X–H length comes from the engine's bond length table (C–H 1.09,
 *     N–H 1.01, O–H 0.96, S–H 1.34, …), clamped to a sane band;
 *   - the H–X–Y and H–X–H angles are projected onto the centre element's
 *     theoretical angle for the centre's hybridization, read off the bond
 *     orders the structure already has the same way the engine's own
 *     hybridization pass reads them (a triple bond -> sp, a double or an
 *     aromatic bond -> sp2, otherwise sp3).
 *
 * Directions are found in two steps: a deterministic greedy start (each
 * hydrogen takes the direction farthest from every existing bond and every
 * hydrogen placed before it), then Gauss-Seidel projection of the angle
 * constraints in which only the new hydrogens ever move. A final pass rotates
 * a hydrogen out of the way when it would sit closer to a non-bonded atom than
 * 0.62 of their van der Waals radii allows — a safety net that sits below any
 * distance real chemistry produces, so a placement that is already sane comes
 * through it untouched.
 *
 * Contract (the one the parsers keep): a structured result, never throws,
 * never mutates its input, no randomness — one structure always yields one
 * answer.
 */

import {
    empiricalBondAngle,
    empiricalBondLength,
    empiricalValence,
    empiricalVdwRadius,
} from './smilesTo3d';

/** New hydrogens one call may add before it stops and says so. */
const MAX_NEW_ATOMS = 5000;
/** Sanity band for an X–H length the table hands out. */
const MIN_XH_LENGTH = 0.70;
const MAX_XH_LENGTH = 2.50;
/** Hybridization codes, the engine's own convention. */
const SP3 = 0;
const SP2 = 1;
const SP = 2;
/** Hydrogens one centre may receive — the empirical table never asks for more. */
const MAX_PER_CENTRE = 6;
/** Angle relaxation. */
const ANGLE_ITERATIONS = 300;
const ANGLE_STEP_DEG = 12.0;
const ANGLE_DAMPING = 0.5;
const ANGLE_TOLERANCE_DEG = 0.02;
/** Direction candidates of the greedy start (deterministic Fibonacci sphere). */
const DIRECTION_CANDIDATES = 64;
/** Non-bonded clash: below this fraction of the van der Waals sum the hydrogen moves away. */
const CONTACT_FACTOR = 0.62;
const AVOID_ROUNDS = 4;
/** Rotation steps tried when a hydrogen is in the way of another atom (degrees). */
const AVOID_STEPS_DEG = [5.0, 10.0, 15.0, 20.0];
/**
 * Floor a placed hydrogen may not come below against any atom it is not bonded
 * to. Shorter than the shortest real bond (H–H, 0.74 Å), so a structure that
 * puts two centres on top of each other cannot smuggle a contact shorter than
 * a chemical bond into the viewer — such a hydrogen is left out and reported.
 */
const MIN_PLACED_CONTACT = 0.80;
/** Cell size of the clash grid; larger than any contact it has to find. */
const CONTACT_CELL = 3.0;
/** Report cap: a pathological selection cannot turn the message into a wall of text. */
const MAX_SKIP_REPORT = 40;

/** One new hydrogen: its centre in the input structure, and where it goes. */
export interface AddedHydrogenAtom {
    element: string;
    x: number;
    y: number;
    z: number;
}

/** One new X–H bond: `atom1` in the input structure, `atom2` = inputCount + its position in `atoms`. */
export interface AddedHydrogenBond {
    atom1: number;
    atom2: number;
    order: number;
}

/** A selected atom that could not be completed, with the numbers behind the verdict. */
export interface SkippedAtomInfo {
    index: number;
    element: string;
    valence: number;
    target: number;
}

export interface AddHydrogensResult {
    /** The new hydrogens, to be appended to the structure in this order. */
    atoms: AddedHydrogenAtom[];
    /** The new X–H bonds, indices already resolved against the extended atom list. */
    bonds: AddedHydrogenBond[];
    /** Centres that received at least one hydrogen. */
    centres: number;
    /** Selected atoms whose valence is already above their empirical valence. */
    overValent: SkippedAtomInfo[];
    /** Selected atoms of elements the empirical table does not cover. */
    unknownValence: number;
    /** Selected entries that were not usable indices or had no usable coordinates. */
    invalid: number;
    /** Hydrogens the structure had no room for (a non-bonded contact would have been shorter than a bond). */
    notPlaced: number;
    /** True when the MAX_NEW_ATOMS budget cut the addition short. */
    truncated: boolean;
    /** One-line summary for the status bar. */
    message: string;
}

/** A hydrogen waiting for its final position; `dir` is a unit vector from its centre. */
interface PendingHydrogen {
    anchor: number;
    dir: number[];
    length: number;
}

/**
 * Add hydrogens to the selected atoms of a structure.
 *
 * `rawAtoms` / `rawBonds` / `rawSelected` are taken as `unknown` on purpose:
 * this runs on whatever the webview sends, so a malformed record degrades the
 * one atom it belongs to instead of failing the call.
 */
export function addHydrogens(rawAtoms: unknown, rawBonds: unknown, rawSelected: unknown): AddHydrogensResult {
    const out: AddHydrogensResult = {
        atoms: [],
        bonds: [],
        centres: 0,
        overValent: [],
        unknownValence: 0,
        invalid: 0,
        notPlaced: 0,
        truncated: false,
        message: '',
    };
    try {
        placeHydrogens(rawAtoms, rawBonds, rawSelected, out);
    } catch (e) {
        const msg = e instanceof Error ? e.message : String(e);
        out.atoms = [];
        out.bonds = [];
        out.centres = 0;
        out.truncated = false;
        out.message = 'Add H failed: ' + msg;
        return out;
    }
    out.message = describe(out);
    return out;
}

// ---------------------------------------------------------------------------
// Placement
// ---------------------------------------------------------------------------

function placeHydrogens(rawAtoms: unknown, rawBonds: unknown, rawSelected: unknown, out: AddHydrogensResult): void {
    const src = Array.isArray(rawAtoms) ? rawAtoms : [];
    const n = src.length;
    const elements: string[] = new Array(n);
    const xs: number[] = new Array(n);
    const ys: number[] = new Array(n);
    const zs: number[] = new Array(n);
    const usable: boolean[] = new Array(n);
    for (let i = 0; i < n; i++) {
        const a = src[i] as { element?: unknown; x?: unknown; y?: unknown; z?: unknown } | null;
        elements[i] = normaliseElement(a ? a.element : '');
        xs[i] = toNumber(a ? a.x : undefined);
        ys[i] = toNumber(a ? a.y : undefined);
        zs[i] = toNumber(a ? a.z : undefined);
        usable[i] = Number.isFinite(xs[i]) && Number.isFinite(ys[i]) && Number.isFinite(zs[i]);
    }

    // Adjacency and total valence from the bonds the caller already has.
    const adj: number[][] = new Array(n);
    const adjOrders: number[][] = new Array(n);
    for (let i = 0; i < n; i++) {
        adj[i] = [];
        adjOrders[i] = [];
    }
    const valence = new Float64Array(n);
    const bondList = Array.isArray(rawBonds) ? rawBonds : [];
    for (let k = 0; k < bondList.length; k++) {
        const b = bondList[k] as { atom1?: unknown; atom2?: unknown; order?: unknown } | null;
        if (!b) {
            continue;
        }
        const i = toIndex(b.atom1, n);
        const j = toIndex(b.atom2, n);
        if (i < 0 || j < 0 || i === j) {
            continue;
        }
        let order = toNumber(b.order);
        if (!Number.isFinite(order) || order <= 0) {
            order = 1;
        }
        if (order > 3) {
            order = 3;
        }
        adj[i].push(j);
        adjOrders[i].push(order);
        adj[j].push(i);
        adjOrders[j].push(order);
        valence[i] += order;
        valence[j] += order;
    }

    const wanted = selectedIndices(rawSelected, n, out);
    const placed: PendingHydrogen[] = [];
    let pending = 0;
    for (let t = 0; t < wanted.length; t++) {
        const idx = wanted[t];
        if (!usable[idx]) {
            out.invalid++;
            continue;
        }
        const el = elements[idx];
        const target = empiricalValence(el);
        if (!(target > 0)) {
            // No empirical valence for this element: completing it would be a guess.
            out.unknownValence++;
            continue;
        }
        const v = valence[idx];
        const count = Math.min(Math.round(target - v), MAX_PER_CENTRE);
        if (count < 1) {
            // Nothing is added when the valence is already at (or above) the
            // empirical value — hydrogens cannot take a bond away.
            if (v > target + 1e-9) {
                out.overValent.push({ index: idx, element: el, valence: v, target: target });
            }
            continue;
        }
        if (pending + count > MAX_NEW_ATOMS) {
            // The budget counts the hydrogens this call would place, not the
            // ones already written out — placement happens after the loop.
            out.truncated = true;
            break;
        }
        // Directions of the bonds that are already there — the hydrogens are
        // the only things that will ever move, so these stay as they are.
        const fixed: number[][] = [];
        for (let k = 0; k < adj[idx].length; k++) {
            const j = adj[idx][k];
            if (!usable[j]) {
                continue;
            }
            const dir = unitVector(xs[idx], ys[idx], zs[idx], xs[j], ys[j], zs[j]);
            if (dir !== null) {
                fixed.push(dir);
            }
        }
        const hybrid = hybridizationOf(adjOrders[idx], neighbourConjugated(idx, adj, adjOrders), el);
        const angle = empiricalBondAngle(el, hybrid);
        const length = clampLength(empiricalBondLength(el, 'H', 1));
        const dirs = initialDirections(fixed, count);
        relaxDirections(dirs, fixed, angle);
        if (hybrid === SP2) {
            snapToPlane(dirs, fixed);
        }
        for (let d = 0; d < dirs.length; d++) {
            placed.push({ anchor: idx, dir: dirs[d], length: length });
        }
        pending += count;
    }
    if (placed.length === 0) {
        return;
    }

    avoidContacts(placed, n, xs, ys, zs, elements, usable);
    const fitting = dropOverlaps(placed, n, xs, ys, zs, elements, usable, out);
    // Centres are counted after the clash filter: an atom whose only hydrogen
    // had to be left out did not receive one.
    const served: { [key: number]: boolean } = {};
    for (let h = 0; h < fitting.length; h++) {
        served[fitting[h].anchor] = true;
    }
    out.centres = Object.keys(served).length;

    for (let h = 0; h < fitting.length; h++) {
        const p = fitting[h];
        out.atoms.push({
            element: 'H',
            x: xs[p.anchor] + p.dir[0] * p.length,
            y: ys[p.anchor] + p.dir[1] * p.length,
            z: zs[p.anchor] + p.dir[2] * p.length,
        });
        out.bonds.push({ atom1: p.anchor, atom2: n + h, order: 1 });
    }
}

/** Selected indices, de-duplicated and range-checked; unusable entries are counted. */
function selectedIndices(rawSelected: unknown, n: number, out: AddHydrogensResult): number[] {
    const wanted: number[] = [];
    const seen: { [key: number]: boolean } = {};
    const list = Array.isArray(rawSelected) ? rawSelected : [];
    for (let t = 0; t < list.length; t++) {
        const idx = toIndex(list[t], n);
        if (idx < 0) {
            out.invalid++;
            continue;
        }
        if (seen[idx] === true) {
            continue;
        }
        seen[idx] = true;
        wanted.push(idx);
    }
    return wanted;
}

/**
 * Hybridization of a centre from the bond orders already around it, read the
 * way the engine's own pass reads a bond list: a triple bond is sp, two
 * double bonds on a two-coordinate centre are sp (cumulated), a double or an
 * aromatic bond is sp2, an amide-like nitrogen next to a multiple bond is sp2,
 * everything else is sp3.
 */
function hybridizationOf(orders: number[], neighbourConjugated: boolean, element: string): number {
    let doubles = 0;
    let unsaturated = false;
    for (let i = 0; i < orders.length; i++) {
        const o = orders[i];
        if (o >= 2.5) {
            return SP;
        }
        if (o >= 1.75) {
            doubles++;
        } else if (o >= 1.25) {
            unsaturated = true;
        }
    }
    if (doubles >= 2 && orders.length === 2) {
        return SP;
    }
    if (doubles >= 1 || unsaturated) {
        return SP2;
    }
    if (neighbourConjugated && (element === 'N' || element === 'P')) {
        return SP2;
    }
    return SP3;
}

/** True when a neighbour carries a multiple bond of its own — an amide/phosphoramide centre. */
function neighbourConjugated(centre: number, adj: number[][], adjOrders: number[][]): boolean {
    for (let k = 0; k < adj[centre].length; k++) {
        const j = adj[centre][k];
        for (let m = 0; m < adj[j].length; m++) {
            if (adj[j][m] === centre) {
                continue;
            }
            if (adjOrders[j][m] >= 1.25) {
                return true;
            }
        }
    }
    return false;
}

/** X–H length inside a band that keeps a corrupt table entry from breaking the view. */
function clampLength(length: number): number {
    if (!Number.isFinite(length) || length <= 0) {
        return 1.09;
    }
    if (length < MIN_XH_LENGTH) {
        return MIN_XH_LENGTH;
    }
    if (length > MAX_XH_LENGTH) {
        return MAX_XH_LENGTH;
    }
    return length;
}

// ---------------------------------------------------------------------------
// Direction search
// ---------------------------------------------------------------------------

/** Greedy start: each hydrogen takes the candidate farthest from every direction already taken. */
function initialDirections(fixed: number[][], count: number): number[][] {
    const out: number[][] = [];
    for (let c = 0; c < count; c++) {
        let best = CANDIDATE_DIRECTIONS[0];
        let bestMin = -1;
        for (let i = 0; i < CANDIDATE_DIRECTIONS.length; i++) {
            const cand = CANDIDATE_DIRECTIONS[i];
            let minAngle = Infinity;
            for (let j = 0; j < fixed.length; j++) {
                const a = angleBetween(cand, fixed[j]);
                if (a < minAngle) {
                    minAngle = a;
                }
            }
            for (let j = 0; j < out.length; j++) {
                const a = angleBetween(cand, out[j]);
                if (a < minAngle) {
                    minAngle = a;
                }
            }
            if (minAngle > bestMin) {
                bestMin = minAngle;
                best = cand;
            }
        }
        out.push([best[0], best[1], best[2]]);
    }
    return out;
}

/**
 * Gauss-Seidel projection of the angle constraints: every pair of directions
 * around a centre — existing bond to new hydrogen, and new hydrogen to new
 * hydrogen — is driven towards the element's theoretical angle. Only the
 * hydrogens move; the bonds that were already there are the frame the
 * placement has to fit into.
 */
function relaxDirections(moving: number[][], fixed: number[][], targetDeg: number): void {
    if (moving.length === 0 || (moving.length === 1 && fixed.length === 0)) {
        return;
    }
    for (let iter = 0; iter < ANGLE_ITERATIONS; iter++) {
        let worst = 0;
        for (let i = 0; i < moving.length; i++) {
            for (let j = 0; j < fixed.length; j++) {
                const step = angleStep(moving[i], fixed[j], targetDeg);
                if (step > worst) {
                    worst = step;
                }
            }
            for (let j = 0; j < moving.length; j++) {
                if (j === i) {
                    continue;
                }
                const step = angleStep(moving[i], moving[j], targetDeg);
                if (step > worst) {
                    worst = step;
                }
            }
        }
        if (worst < ANGLE_TOLERANCE_DEG) {
            return;
        }
    }
}

/**
 * sp2 centres: the angle constraints converge to within their tolerance long
 * before the out-of-plane component does — a coplanar direction is an extremum
 * of every angle around such a centre, so the solver is blind to it, and a
 * residual of a few hundredths of a degree shows up as a hydrogen a tenth of
 * an Angstrom out of plane. Projecting the new directions onto the plane of
 * the bonds already present removes it without touching the angles.
 */
function snapToPlane(dirs: number[][], fixed: number[][]): void {
    if (fixed.length < 2) {
        return;
    }
    const normal = bestFitNormal(fixed);
    if (normal === null) {
        return;
    }
    for (let d = 0; d < dirs.length; d++) {
        const dir = dirs[d];
        const dot = dir[0] * normal[0] + dir[1] * normal[1] + dir[2] * normal[2];
        const px = dir[0] - dot * normal[0];
        const py = dir[1] - dot * normal[1];
        const pz = dir[2] - dot * normal[2];
        const len = Math.sqrt(px * px + py * py + pz * pz);
        if (len > 1e-6) {
            dir[0] = px / len;
            dir[1] = py / len;
            dir[2] = pz / len;
        }
    }
}

/** Normal of the plane the existing bonds span: the cross product of the pair least parallel to each other. */
function bestFitNormal(fixed: number[][]): number[] | null {
    let best: number[] | null = null;
    let bestCross = -1;
    for (let i = 0; i < fixed.length; i++) {
        for (let j = i + 1; j < fixed.length; j++) {
            const n = cross(fixed[i], fixed[j]);
            const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            if (len > bestCross) {
                bestCross = len;
                best = [n[0] / len, n[1] / len, n[2] / len];
            }
        }
    }
    return bestCross > 1e-6 ? best : null;
}

/** One projection step of `dir` against `other`; returns the angle change asked for. */
function angleStep(dir: number[], other: number[], targetDeg: number): number {
    const current = angleBetween(dir, other);
    let step = current - targetDeg;
    if (Math.abs(step) < ANGLE_TOLERANCE_DEG) {
        return 0;
    }
    step *= ANGLE_DAMPING;
    if (step > ANGLE_STEP_DEG) {
        step = ANGLE_STEP_DEG;
    } else if (step < -ANGLE_STEP_DEG) {
        step = -ANGLE_STEP_DEG;
    }
    // Rotating `dir` about (dir x other) by a positive angle turns it towards
    // `other`, which closes the angle — so the signed excess is the step.
    let axis = cross(dir, other);
    const len = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if (len > 1e-8) {
        axis = [axis[0] / len, axis[1] / len, axis[2] / len];
    } else {
        const fallback = perpendicularTo(dir);
        if (fallback === null) {
            return 0;
        }
        axis = fallback;
    }
    const turned = rotateAbout(dir, axis, (step * Math.PI) / 180.0);
    dir[0] = turned[0];
    dir[1] = turned[1];
    dir[2] = turned[2];
    return Math.abs(step);
}

// ---------------------------------------------------------------------------
// Non-bonded clash avoidance
// ---------------------------------------------------------------------------

/** One entry of the clash grid: an existing atom, or a hydrogen already placed. */
interface ContactPoint {
    x: number;
    y: number;
    z: number;
    radius: number;
    atom: number;
    hydrogen: number;
}

/**
 * Rotate a hydrogen around its centre when a non-bonded atom sits closer than
 * CONTACT_FACTOR of their van der Waals sum. Each round tries a few rotation
 * steps away from the offending atom and keeps the one that clears the most
 * space, so a hydrogen between two atoms moves to the larger gap instead of
 * being pushed back and forth. Rotating keeps the X–H length exact, and the
 * threshold is well below any distance real chemistry produces, so ordinary
 * placements are left alone and only genuine overlaps move.
 */
function avoidContacts(placed: PendingHydrogen[], n: number, xs: number[], ys: number[], zs: number[],
                       elements: string[], usable: boolean[]): void {
    const hRadius = empiricalVdwRadius('H');
    for (let round = 0; round < AVOID_ROUNDS; round++) {
        const grid = buildContactGrid(placed, n, xs, ys, zs, elements, usable);
        let moved = false;
        for (let h = 0; h < placed.length; h++) {
            const p = placed[h];
            const clearance = hydrogenClearance(p, h, grid, hRadius, xs, ys, zs);
            if (clearance.ratio >= CONTACT_FACTOR || clearance.blocker === null) {
                continue;
            }
            // Turn the hydrogen away from the blocker: a rotation about the axis
            // perpendicular to both directions opens the angle between them and
            // leaves the X–H length untouched. Both senses and a few sizes are
            // tried, because a hydrogen can sit between two atoms.
            const blocker = clearance.blocker;
            const bdir = unitVector(xs[p.anchor], ys[p.anchor], zs[p.anchor], blocker.x, blocker.y, blocker.z);
            if (bdir === null) {
                continue;
            }
            const axis = rotationAxis(p.dir, bdir);
            if (axis === null) {
                continue;
            }
            const home = [p.dir[0], p.dir[1], p.dir[2]];
            let bestRatio = clearance.ratio;
            let bestDir: number[] | null = null;
            for (let s = 0; s < AVOID_STEPS_DEG.length; s++) {
                for (let sign = -1; sign <= 1; sign += 2) {
                    const turned = rotateAbout(home, axis, (sign * AVOID_STEPS_DEG[s] * Math.PI) / 180.0);
                    p.dir = turned;
                    const candidate = hydrogenClearance(p, h, grid, hRadius, xs, ys, zs).ratio;
                    if (candidate > bestRatio) {
                        bestRatio = candidate;
                        bestDir = turned;
                    }
                }
            }
            if (bestDir !== null) {
                p.dir = bestDir;
                moved = true;
            } else {
                p.dir = home;
            }
        }
        if (!moved) {
            return;
        }
    }
}

/** Closest non-bonded approach of one hydrogen: distance over the van der Waals sum. */
function hydrogenClearance(p: PendingHydrogen, h: number, grid: ContactGrid, hRadius: number,
                           xs: number[], ys: number[], zs: number[]): { ratio: number; blocker: ContactPoint | null } {
    const ax = xs[p.anchor] + p.dir[0] * p.length;
    const ay = ys[p.anchor] + p.dir[1] * p.length;
    const az = zs[p.anchor] + p.dir[2] * p.length;
    const blockers = grid.query(ax, ay, az);
    let ratio = Infinity;
    let blocker: ContactPoint | null = null;
    for (let t = 0; t < blockers.length; t++) {
        const q = blockers[t];
        if (q.hydrogen === h || q.atom === p.anchor) {
            continue;
        }
        const dx = q.x - ax;
        const dy = q.y - ay;
        const dz = q.z - az;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const value = dist / (q.radius + hRadius);
        if (value < ratio) {
            ratio = value;
            blocker = q;
        }
    }
    return { ratio: ratio === Infinity ? CONTACT_FACTOR : ratio, blocker: blocker };
}

/** Axis perpendicular to both directions; null when they are parallel (any perpendicular will do then). */
function rotationAxis(dir: number[], other: number[]): number[] | null {
    const axis = cross(dir, other);
    const len = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if (len > 1e-8) {
        return [axis[0] / len, axis[1] / len, axis[2] / len];
    }
    return perpendicularTo(dir);
}

/**
 * A hydrogen that would still sit closer than MIN_PLACED_CONTACT of an atom it
 * is not bonded to is left out: two centres a fraction of an Angstrom apart
 * (a coordinate defect, not chemistry) leave no room for one, and a viewer
 * should not gain an atom overlap from an edit. Returns the hydrogens that fit.
 */
function dropOverlaps(placed: PendingHydrogen[], n: number, xs: number[], ys: number[], zs: number[],
                      elements: string[], usable: boolean[], out: AddHydrogensResult): PendingHydrogen[] {
    if (placed.length === 0) {
        return placed;
    }
    const grid = buildContactGrid(placed, n, xs, ys, zs, elements, usable);
    const kept: PendingHydrogen[] = [];
    for (let h = 0; h < placed.length; h++) {
        const p = placed[h];
        const ax = xs[p.anchor] + p.dir[0] * p.length;
        const ay = ys[p.anchor] + p.dir[1] * p.length;
        const az = zs[p.anchor] + p.dir[2] * p.length;
        const blockers = grid.query(ax, ay, az);
        let room = true;
        for (let t = 0; t < blockers.length; t++) {
            const q = blockers[t];
            if (q.hydrogen === h || q.atom === p.anchor) {
                continue;
            }
            const dx = q.x - ax;
            const dy = q.y - ay;
            const dz = q.z - az;
            if (Math.sqrt(dx * dx + dy * dy + dz * dz) < MIN_PLACED_CONTACT) {
                room = false;
                break;
            }
        }
        if (room) {
            kept.push(p);
        } else {
            out.notPlaced++;
        }
    }
    return kept;
}

function buildContactGrid(placed: PendingHydrogen[], n: number, xs: number[], ys: number[], zs: number[],
                          elements: string[], usable: boolean[]): ContactGrid {
    const grid = new ContactGrid();
    for (let i = 0; i < n; i++) {
        if (!usable[i]) {
            continue;
        }
        grid.add(xs[i], ys[i], zs[i], empiricalVdwRadius(elements[i]), i, -1);
    }
    for (let h = 0; h < placed.length; h++) {
        const p = placed[h];
        grid.add(
            xs[p.anchor] + p.dir[0] * p.length,
            ys[p.anchor] + p.dir[1] * p.length,
            zs[p.anchor] + p.dir[2] * p.length,
            empiricalVdwRadius('H'),
            -1,
            h
        );
    }
    return grid;
}

/** Spatial hash of the atoms a clash check has to look at (cell size CONTACT_CELL). */
class ContactGrid {
    private cells: { [key: string]: ContactPoint[] } = {};
    private points: ContactPoint[] = [];

    add(x: number, y: number, z: number, radius: number, atom: number, hydrogen: number): void {
        const point: ContactPoint = { x: x, y: y, z: z, radius: radius, atom: atom, hydrogen: hydrogen };
        this.points.push(point);
        const key = this.key(x, y, z);
        const cell = this.cells[key];
        if (cell === undefined) {
            this.cells[key] = [point];
        } else {
            cell.push(point);
        }
    }

    query(x: number, y: number, z: number): ContactPoint[] {
        const cx = Math.floor(x / CONTACT_CELL);
        const cy = Math.floor(y / CONTACT_CELL);
        const cz = Math.floor(z / CONTACT_CELL);
        const out: ContactPoint[] = [];
        for (let i = -1; i <= 1; i++) {
            for (let j = -1; j <= 1; j++) {
                for (let k = -1; k <= 1; k++) {
                    const cell = this.cells[(cx + i) + ',' + (cy + j) + ',' + (cz + k)];
                    if (cell !== undefined) {
                        for (let t = 0; t < cell.length; t++) {
                            out.push(cell[t]);
                        }
                    }
                }
            }
        }
        return out;
    }

    private key(x: number, y: number, z: number): string {
        return Math.floor(x / CONTACT_CELL) + ',' + Math.floor(y / CONTACT_CELL) + ',' + Math.floor(z / CONTACT_CELL);
    }
}

// ---------------------------------------------------------------------------
// Small vector helpers
// ---------------------------------------------------------------------------

/** Unit vector from a to b, or null when the two points coincide. */
function unitVector(ax: number, ay: number, az: number, bx: number, by: number, bz: number): number[] | null {
    const dx = bx - ax;
    const dy = by - ay;
    const dz = bz - az;
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (!(len > 1e-9)) {
        return null;
    }
    return [dx / len, dy / len, dz / len];
}

function angleBetween(a: number[], b: number[]): number {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    if (dot > 1) {
        dot = 1;
    } else if (dot < -1) {
        dot = -1;
    }
    return (Math.acos(dot) * 180.0) / Math.PI;
}

function cross(a: number[], b: number[]): number[] {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
}

/** Any unit vector perpendicular to `v`; deterministic, so results stay reproducible. */
function perpendicularTo(v: number[]): number[] | null {
    const ax = Math.abs(v[0]);
    const ay = Math.abs(v[1]);
    const az = Math.abs(v[2]);
    let ref: number[];
    if (ax <= ay && ax <= az) {
        ref = [1, 0, 0];
    } else if (ay <= az) {
        ref = [0, 1, 0];
    } else {
        ref = [0, 0, 1];
    }
    const axis = cross(v, ref);
    const len = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if (!(len > 1e-9)) {
        return null;
    }
    return [axis[0] / len, axis[1] / len, axis[2] / len];
}

/** Rodrigues rotation of a unit vector about a unit axis. */
function rotateAbout(v: number[], axis: number[], radians: number): number[] {
    const cos = Math.cos(radians);
    const sin = Math.sin(radians);
    const k = cross(axis, v);
    const dot = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
    return [
        v[0] * cos + k[0] * sin + axis[0] * dot * (1 - cos),
        v[1] * cos + k[1] * sin + axis[1] * dot * (1 - cos),
        v[2] * cos + k[2] * sin + axis[2] * dot * (1 - cos),
    ];
}

/** Deterministic Fibonacci sphere, the greedy start's candidate set. */
const CANDIDATE_DIRECTIONS: number[][] = (function (): number[][] {
    const out: number[][] = [];
    const golden = Math.PI * (3.0 - Math.sqrt(5.0));
    for (let i = 0; i < DIRECTION_CANDIDATES; i++) {
        const y = 1.0 - (i / (DIRECTION_CANDIDATES - 1)) * 2.0;
        const r = Math.sqrt(Math.max(0.0, 1.0 - y * y));
        const th = golden * i;
        out.push([Math.cos(th) * r, y, Math.sin(th) * r]);
    }
    return out;
})();

// ---------------------------------------------------------------------------
// Input normalisation and the result summary
// ---------------------------------------------------------------------------

function normaliseElement(value: unknown): string {
    if (typeof value !== 'string' || value.length === 0) {
        return '';
    }
    return value.charAt(0).toUpperCase() + value.slice(1).toLowerCase();
}

function toNumber(value: unknown): number {
    if (typeof value === 'number') {
        return value;
    }
    const parsed = typeof value === 'string' ? Number(value) : NaN;
    return parsed;
}

/** Index from a value the caller sent; -1 when it is not a whole index in range. */
function toIndex(value: unknown, n: number): number {
    const num = typeof value === 'number' ? value : NaN;
    if (!Number.isFinite(num) || Math.floor(num) !== num || num < 0 || num >= n) {
        return -1;
    }
    return num;
}

function describe(out: AddHydrogensResult): string {
    if (out.atoms.length === 0) {
        const reasons: string[] = [];
        if (out.centres === 0 && out.overValent.length === 0 && out.unknownValence === 0 && out.invalid === 0 && out.notPlaced === 0) {
            return 'Add H: nothing to complete in the selection';
        }
        if (out.notPlaced > 0) {
            reasons.push(out.notPlaced + ' H left out (no room at the centre)');
        }
        if (out.overValent.length > 0) {
            reasons.push(out.overValent.length + ' atom(s) already at or above their empirical valence');
        }
        if (out.unknownValence > 0) {
            reasons.push(out.unknownValence + ' atom(s) of elements with no empirical valence');
        }
        if (out.invalid > 0) {
            reasons.push(out.invalid + ' invalid selection/coordinate entr' + (out.invalid === 1 ? 'y' : 'ies'));
        }
        return 'Add H: nothing added (' + reasons.join('; ') + ')';
    }
    let text = 'Add H: added ' + out.atoms.length + ' H to ' + out.centres + ' atom' + (out.centres === 1 ? '' : 's');
    const notes: string[] = [];
    if (out.notPlaced > 0) {
        notes.push(out.notPlaced + ' H left out (no room: a contact would have been shorter than a bond)');
    }
    if (out.overValent.length > 0) {
        notes.push(out.overValent.length + ' skipped (valence above the empirical value)');
    }
    if (out.unknownValence > 0) {
        notes.push(out.unknownValence + ' skipped (no empirical valence)');
    }
    if (out.invalid > 0) {
        notes.push(out.invalid + ' invalid');
    }
    if (out.truncated) {
        notes.push('budget of ' + MAX_NEW_ATOMS + ' H reached');
    }
    if (notes.length > 0) {
        text += ' — ' + notes.join('; ');
    }
    if (out.overValent.length > 0 && out.overValent.length <= MAX_SKIP_REPORT) {
        text += ' [' + out.overValent.map(function (s) {
            return s.element + (s.index + 1) + ' val ' + round1(s.valence) + '/' + s.target;
        }).join(', ') + ']';
    }
    return text;
}

function round1(v: number): number {
    return Math.round(v * 10) / 10;
}