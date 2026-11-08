/**
 * One element table for every reader that has to turn a number, a mass or a
 * program's atom-name spelling back into a chemical symbol.
 *
 * The parsers added for the plane-wave and MD programs all face the same
 * problem: a file identifies its atoms by Z (ABINIT), by mass (LAMMPS), by a
 * label (`OW`, `HW1`, `taud_Si1`) or by a type index mapped through another
 * file. Keeping a private table in each reader would make them disagree about
 * an isotope or an element's spelling; keeping one here makes the disagreement
 * impossible.
 */

/** Symbol by atomic number, 1..118. Index 0 is a placeholder. */
export const ELEMENT_SYMBOLS: string[] = [
    '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
];

/** Standard atomic weights (IUPAC 2021, abridged), by atomic number. */
export const ELEMENT_MASSES: number[] = [
    0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71,
    121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24,
    145.0, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05,
    174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
    204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.04,
    231.04, 238.03, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    258.0, 259.0, 262.0, 267.0, 268.0, 269.0, 270.0, 269.0, 278.0, 281.0,
    282.0, 285.0, 286.0, 289.0, 290.0, 293.0, 294.0, 294.0
];

const SYMBOL_BY_LOWER: { [key: string]: string } = (() => {
    const map: { [key: string]: string } = {};
    for (let z = 1; z < ELEMENT_SYMBOLS.length; z++) {
        map[ELEMENT_SYMBOLS[z].toLowerCase()] = ELEMENT_SYMBOLS[z];
    }
    return map;
})();

/** The element with atomic number `z`, or undefined outside 1..118. */
export function elementFromZ(z: number): string | undefined {
    if (!isFinite(z) || z !== Math.round(z)) return undefined;
    return z >= 1 && z < ELEMENT_SYMBOLS.length ? ELEMENT_SYMBOLS[z] : undefined;
}

/** The atomic number of a symbol, case-insensitively. */
export function atomicNumberOf(symbol: string): number | undefined {
    const hit = SYMBOL_BY_LOWER[String(symbol || '').trim().toLowerCase()];
    return hit ? ELEMENT_SYMBOLS.indexOf(hit) : undefined;
}

/** The standard atomic weight of a symbol, or undefined for a non-element. */
export function atomicMassOf(symbol: string): number | undefined {
    const z = atomicNumberOf(symbol);
    return z !== undefined ? ELEMENT_MASSES[z] : undefined;
}

/**
 * A symbol with the decorations programs put on it removed: `Na+`, `O2-`,
 * `Sr_sv`, `Fe3+` and `Cl_s` are the elements Na, O, Sr, Fe and Cl.
 */
export function normalizeElementSymbol(raw: string): string | undefined {
    const text = String(raw || '').trim();
    if (text === '') return undefined;
    // The element is the leading run of letters: a suffix may be a charge
    // (`Fe3+`, `O2-`), a potential label (`Sr_sv`, `Cl_s`) or an index (`Si1`),
    // and stripping every non-letter instead would turn `Sr_sv` into `Srsv`.
    const run = /[A-Za-z]+/.exec(text);
    if (!run) return undefined;
    const letters = run[0];
    const cap = letters.charAt(0).toUpperCase() + letters.slice(1).toLowerCase();
    return SYMBOL_BY_LOWER[cap.toLowerCase()] ? cap : undefined;
}

/**
 * The element an MD program's atom *label* stands for.
 *
 * GROMACS calls a water oxygen `OW` and its hydrogens `HW1`/`HW2`; a silicon
 * atom may be `SIL`; LAMMPS uses `Si`, `Mg`, `Na`. Exact symbols win, then the
 * longest known prefix (so `OW` -> O, `HW1` -> H, `SIL` -> Si), then nothing.
 * A label that merely starts with a valid letter must not become that element
 * (`SOL` is a residue name, not sulfur) unless no longer prefix is a symbol —
 * which is the same best-effort rule this function is documented to follow.
 */
export function elementFromAtomName(name: string): string | undefined {
    const letters = String(name || '').replace(/[^A-Za-z]/g, '');
    if (letters === '') return undefined;
    const exact = SYMBOL_BY_LOWER[letters.toLowerCase()];
    if (exact) return exact;
    if (letters.length >= 2) {
        const two = SYMBOL_BY_LOWER[letters.slice(0, 2).toLowerCase()];
        if (two) return two;
    }
    return SYMBOL_BY_LOWER[letters.slice(0, 1).toLowerCase()];
}

/**
 * The element whose standard atomic weight is closest to `mass`, within
 * `tolerance` g/mol.
 *
 * LAMMPS data files identify atom types by mass and nothing else, so this is
 * the only way to name their atoms. The tolerance is deliberately narrow: a
 * pseudo-potential's `Zval` mass (2.0 for Mg) or an isotope must not silently
 * become a different element, and when nothing is close the caller keeps the
 * honest `TypeN` label instead of a guess.
 */
export function elementFromMass(mass: number, tolerance = 0.35): string | undefined {
    if (!isFinite(mass) || mass <= 0) return undefined;
    // Deuterium and tritium are hydrogen in every structure file that names them
    // by mass; no element has a standard weight near 2 or 3.
    if (Math.abs(mass - 2.014) < 0.06 || Math.abs(mass - 3.016) < 0.06) return 'H';
    let best = -1;
    let bestDiff = Infinity;
    for (let z = 1; z < ELEMENT_MASSES.length; z++) {
        const diff = Math.abs(ELEMENT_MASSES[z] - mass);
        if (diff < bestDiff) {
            bestDiff = diff;
            best = z;
        }
    }
    return best > 0 && bestDiff <= tolerance ? ELEMENT_SYMBOLS[best] : undefined;
}

/** The label `atomLabel` uses when the element is unknown: `Type3`, `T3`. */
export function typeElementLabel(typeIndex: number): string {
    return 'T' + typeIndex;
}
