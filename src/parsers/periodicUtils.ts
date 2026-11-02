import { Atom, CrystalData, MolecularData } from '../types';
import { cartToFrac, fracToCart, wrapFrac } from '../tools/crystalGeometry';

/** CODATA 2018 Bohr radius, the value every program below converts with. */
export const BOHR_TO_ANGSTROM = 0.529177210903;

export function bohrToAngstrom(v: number): number {
    return v * BOHR_TO_ANGSTROM;
}

export function angstromToBohr(v: number): number {
    return v / BOHR_TO_ANGSTROM;
}

export interface CellParameters {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
}

/**
 * The six cell parameters a lattice states, in Angstrom and degrees.
 *
 * A reader that has lattice vectors (an ABINIT `rprim`, a LAMMPS triclinic
 * box, a CP2K `A/B/C`) still has to fill `CrystalData.a..gamma`, because every
 * panel — density, cell editor, XRD — reads those fields. The arithmetic is
 * the same for all of them, so it lives here rather than in each reader.
 */
export function latticeToCellParams(lattice: number[][]): CellParameters {
    const [aVec, bVec, cVec] = lattice;
    const norm = (v: number[]) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    const a = norm(aVec);
    const b = norm(bVec);
    const c = norm(cVec);
    const dot = (x: number[], y: number[]) => x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    const safeAcos = (v: number) => Math.acos(Math.max(-1, Math.min(1, v)));
    const rad = 180 / Math.PI;
    return {
        a, b, c,
        alpha: safeAcos(dot(bVec, cVec) / (b * c || 1)) * rad,
        beta: safeAcos(dot(aVec, cVec) / (a * c || 1)) * rad,
        gamma: safeAcos(dot(aVec, bVec) / (a * b || 1)) * rad
    };
}

export interface Site {
    element: string;
    /** Fractional coordinates; wrapped into [0,1) when building a cell. */
    frac?: number[];
    /** Cartesian coordinates in Angstrom, used when `frac` is absent. */
    cart?: number[];
    occupancy?: number;
    constraint?: string;
    label?: string;
    ghost?: boolean;
}

function siteCart(site: Site, lattice: number[][] | null): [number, number, number] {
    if (site.frac) {
        return fracToCart(site.frac, lattice as number[][]);
    }
    const c = site.cart || [0, 0, 0];
    return [c[0], c[1], c[2]];
}

/**
 * A complete cell described by explicit sites and no symmetry.
 *
 * STRU, ABINIT, CP2K, QE, LAMMPS, GRO and the VASP trajectories all list every
 * atom of the cell; the only operation is the identity, and `baseAtoms` *is*
 * the file's atom list. Building that shape in one place keeps the crystal
 * invariants (`baseIdx`, `cellI/J/K`, wrapped fractional positions) identical
 * across readers, which is what the viewer's editing passes rely on.
 */
export function buildPeriodicStructure(
    lattice: number[][],
    sites: Site[],
    title: string,
    warnings?: string[]
): MolecularData {
    const atoms: Atom[] = [];
    for (let i = 0; i < sites.length; i++) {
        const site = sites[i];
        let frac: number[];
        let cart: [number, number, number];
        if (site.frac) {
            frac = [wrapFrac(site.frac[0]), wrapFrac(site.frac[1]), wrapFrac(site.frac[2])];
            cart = fracToCart(frac, lattice);
        } else {
            cart = siteCart(site, lattice);
            frac = cartToFrac(cart, lattice);
        }
        const atom: Atom = {
            element: site.element,
            x: cart[0],
            y: cart[1],
            z: cart[2],
            index: i,
            occupancy: site.occupancy !== undefined && isFinite(site.occupancy) ? site.occupancy : 1,
            baseIdx: i,
            cellI: 0,
            cellJ: 0,
            cellK: 0
        };
        if (site.constraint) atom.constraint = site.constraint;
        if (site.label) atom.siteLabel = site.label;
        if (site.ghost) atom.ghost = true;
        atoms.push(atom);
    }

    const crystal: CrystalData = {
        ...latticeToCellParams(lattice),
        latticeVectors: lattice,
        symmetryOps: ['x, y, z'],
        baseAtoms: atoms.map(a => ({ ...a })),
        baseBonds: []
    };

    const data: MolecularData = {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        crystal
    };
    if (warnings && warnings.length > 0) data.warnings = warnings;
    return data;
}

/** A molecule described by explicit Cartesian sites and no cell. */
export function buildMolecularStructure(
    sites: Site[],
    title: string,
    warnings?: string[]
): MolecularData {
    const atoms: Atom[] = [];
    for (let i = 0; i < sites.length; i++) {
        const site = sites[i];
        const c = site.frac ? [site.frac[0], site.frac[1], site.frac[2]] : (site.cart || [0, 0, 0]);
        const atom: Atom = {
            element: site.element,
            x: c[0],
            y: c[1],
            z: c[2],
            index: i
        };
        if (site.occupancy !== undefined && isFinite(site.occupancy)) atom.occupancy = site.occupancy;
        if (site.label) atom.siteLabel = site.label;
        if (site.ghost) atom.ghost = true;
        atoms.push(atom);
    }
    const data: MolecularData = {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false
    };
    if (warnings && warnings.length > 0) data.warnings = warnings;
    return data;
}

/** Whether a 3x3 lattice has the volume a cell needs to be usable. */
export function isUsableLattice(lattice: number[][]): boolean {
    if (!lattice || lattice.length !== 3) return false;
    for (const row of lattice) {
        if (!row || row.length !== 3 || row.some(v => !isFinite(v))) return false;
    }
    const [a, b, c] = lattice;
    const det = a[0] * (b[1] * c[2] - b[2] * c[1]) -
        a[1] * (b[0] * c[2] - b[2] * c[0]) +
        a[2] * (b[0] * c[1] - b[1] * c[0]);
    return Math.abs(det) > 1e-8;
}

/** `1.0 -> "1.0"` without the `-0` and exponent spellings programs dislike. */
export function formatNumber(v: number, digits: number): string {
    if (!isFinite(v)) return '0';
    if (Math.abs(v) < 0.5 * Math.pow(10, -digits)) v = 0;
    return v.toFixed(digits);
}

/** Trim a title to one line, so a multi-line header cannot break a format. */
export function oneLineTitle(text: string, fallback: string): string {
    const line = String(text || '').split(/\r?\n/).map(s => s.trim()).filter(s => s !== '')[0];
    return line || fallback;
}
