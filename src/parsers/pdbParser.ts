import { Atom, Bond, MolecularData } from '../types';
import { cellToLatticeVectors } from '../tools/crystalGeometry';
import { buildPeriodicStructure, CellParameters, isUsableLattice, Site } from './periodicUtils';

const ATOMIC_NUMBER_MAP: { [key: string]: number } = {
    H: 1, D: 1, HE: 2, LI: 3, BE: 4, B: 5, C: 6, N: 7, O: 8, F: 9, NE: 10,
    NA: 11, MG: 12, AL: 13, SI: 14, P: 15, S: 16, CL: 17, AR: 18, K: 19, CA: 20,
    SC: 21, TI: 22, V: 23, CR: 24, MN: 25, FE: 26, CO: 27, NI: 28, CU: 29, ZN: 30,
    GA: 31, GE: 32, AS: 33, SE: 34, BR: 35, KR: 36, RB: 37, SR: 38, Y: 39, ZR: 40,
    NB: 41, MO: 42, TC: 43, RU: 44, RH: 45, PD: 46, AG: 47, CD: 48, IN: 49, SN: 50,
    SB: 51, TE: 52, I: 53, XE: 54, CS: 55, BA: 56, LA: 57, CE: 58, PR: 59, ND: 60,
    PM: 61, SM: 62, EU: 63, GD: 64, TB: 65, DY: 66, HO: 67, ER: 68, TM: 69, YB: 70,
    LU: 71, HF: 72, TA: 73, W: 74, RE: 75, OS: 76, IR: 77, PT: 78, AU: 79, HG: 80,
    TL: 81, PB: 82, BI: 83, PO: 84, AT: 85, RN: 86, FR: 87, RA: 88, AC: 89, TH: 90,
    PA: 91, U: 92, NP: 93, PU: 94, AM: 95, CM: 96, BK: 97, CF: 98
};

/**
 * The element of an atom, from the PDB element column when it is present and
 * from the atom-name field otherwise.
 *
 * PDB alignment is what tells a protein's C-alpha from a calcium ion: a
 * one-letter element's atom name is right-justified (column 13 is blank,
 * ` CA ` in a protein), while a two-letter element's name starts in column 13
 * (`CA  ` for calcium). Trimming the name first loses that distinction and
 * turned every protein `CA`, `CD`, `CE`, `HG` into a metal.
 */
function pdbElement(atomNameField: string, altElement?: string): string {
    if (altElement && altElement.trim() !== '') {
        const e = altElement.trim();
        const upper = e.toUpperCase();
        if (ATOMIC_NUMBER_MAP[upper] !== undefined) {
            return upper.charAt(0) + upper.charAt(1).toLowerCase();
        }
        return e.charAt(0).toUpperCase() + e.slice(1).toLowerCase();
    }
    const field = atomNameField || '';
    const name = field.trim();
    if (!name) return 'C';
    const startsInColumn13 = field.charAt(0) !== ' ' && field.charAt(0) !== '\t';
    const numericPrefix = /^[0-9]/.test(name);
    const letters = name.replace(/^[0-9]+/, '');
    if (!letters) return 'C';

    const two = letters.slice(0, 2);
    let twoSymbol = '';
    const upper2 = two.toUpperCase();
    if (two.length === 2 && ATOMIC_NUMBER_MAP[upper2] !== undefined) {
        twoSymbol = two.charAt(0).toUpperCase() + two.charAt(1).toLowerCase();
    }
    // Only a name that starts in the element column can name a two-letter
    // element; a right-justified name (` CA `) is a single-letter element.
    if (startsInColumn13 && !numericPrefix && twoSymbol) {
        return twoSymbol;
    }
    const first = letters.charAt(0).toUpperCase();
    if (ATOMIC_NUMBER_MAP[first] !== undefined) {
        return first;
    }
    if (twoSymbol) return twoSymbol;
    return first;
}

/**
 * One bond a `CONECT` record states, keyed by its two serials.
 *
 * The record names a partner once per bond, so a partner listed twice is a
 * double bond — the format's only way of stating an order at all. Records are
 * normally written from both ends (`CONECT    1    2    2` beside
 * `CONECT    2    1`), which is why the two directions are counted apart: the
 * order is one direction's count, not the total of both.
 */
interface PdbConect {
    lo: number;
    hi: number;
    loCount: number;
    hiCount: number;
}

/**
 * The unit cell a `CRYST1` record states, and the space-group name beside it.
 *
 * `CRYST1` is the only place a PDB states a cell, and it is fixed-width:
 * a(9.3) in columns 7-15, b in 16-24, c in 25-33, the three angles as 7.2 in
 * 34-54 and the Hermann-Mauguin symbol in 56-66, with the Z value after it.
 * The columns are the authority as long as they agree with the same six
 * numbers read as tokens; a writer that collapsed the padding still leaves
 * those six numbers as the record's first tokens, and that is what is read
 * when a slice disagrees with them — a shifted column picks up its neighbour
 * (`640 20.0` for b) rather than failing. Anything the record cannot state as
 * a cell — a non-positive axis length, an angle outside 0-180 — is rejected
 * here and the caller reports the record instead of building a degenerate
 * lattice from it.
 */
function readCryst1(line: string): { cell: CellParameters; spaceGroup: string } | null {
    const sliceNumber = (from: number, to: number): number | null => {
        const v = parseFloat(line.substring(from, to).trim());
        return isFinite(v) ? v : null;
    };
    const parts = line.trim().split(/\s+/);
    const toks = parts.slice(1, 7).map(t => parseFloat(t));
    const tokenised = toks.length === 6 && toks.every(v => isFinite(v)) ? toks : null;

    const fixed = [sliceNumber(6, 15), sliceNumber(15, 24), sliceNumber(24, 33),
        sliceNumber(33, 40), sliceNumber(40, 47), sliceNumber(47, 54)];
    const fixedOk = fixed.every(v => v !== null);
    // The columns are trusted only when they agree with the same numbers read
    // as tokens. A writer that collapsed the padding shifts every field, and a
    // column slice then picks up the neighbouring value (`640 20.0` for b), so
    // a disagreement means the line is not in the format's columns — the
    // tokens, which cannot be shifted, are what is read then.
    const aligned = fixedOk && (tokenised === null ||
        fixed.every((v, i) => Math.abs((v as number) - tokenised[i]) <= 1e-6 * Math.max(1, Math.abs(tokenised[i]))));

    let values: number[];
    let spaceGroup: string;
    if (aligned) {
        values = fixed as number[];
        spaceGroup = line.substring(55, 66).trim();
    } else if (tokenised !== null) {
        values = tokenised;
        // The trailing Z is a bare integer; a symbol's settings are digits too
        // (`P 21 21 21`), so only a token after at least two others is taken
        // as the Z value.
        const rest = parts.slice(7);
        if (rest.length >= 3 && /^\d+$/.test(rest[rest.length - 1])) rest.pop();
        spaceGroup = rest.join(' ');
    } else {
        return null;
    }
    const [a, b, c, alpha, beta, gamma] = values;
    if (!(a > 0 && b > 0 && c > 0)) return null;
    if (!(alpha > 0 && alpha < 180 && beta > 0 && beta < 180 && gamma > 0 && gamma < 180)) return null;
    return {
        cell: { a, b, c, alpha, beta, gamma },
        spaceGroup: /[A-Za-z]/.test(spaceGroup) ? spaceGroup : ''
    };
}

/**
 * The cell a molecule-only writer emits when it has none: a 1 Angstrom cube.
 *
 * OpenBabel and several converters write `CRYST1 1.000 1.000 1.000 90 90 90 P
 * 1` for a structure that has no cell at all. Taken at face value it folds
 * every atom into a one-Angstrom box, so it is recognised as the placeholder it
 * is — no crystal has an axis shorter than its own bonds.
 */
function isPlaceholderCell(cell: CellParameters): boolean {
    const near = (v: number, want: number, tol: number) => Math.abs(v - want) <= tol;
    return near(cell.a, 1, 1e-3) && near(cell.b, 1, 1e-3) && near(cell.c, 1, 1e-3) &&
        near(cell.alpha, 90, 1e-2) && near(cell.beta, 90, 1e-2) && near(cell.gamma, 90, 1e-2);
}

export function parsePdb(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let hasExplicitBonds = false;
    let title = '';
    let charge: number | undefined;
    let multiplicity: number | undefined;
    let atomIndex = 0;

    const atomSerialMap = new Map<number, number>();
    const seenSerials = new Set<number>();
    // The first `CRYST1` record wins: a file that repeats it (two jobs
    // concatenated) describes the same cell both times, and the atoms are the
    // ones written under the first. A record that is present but unreadable is
    // remembered so the reader can say so rather than silently drawing a
    // molecule.
    let cell: CellParameters | null = null;
    let cellSpaceGroup = '';
    let cellBroken = false;
    // Collected rather than pushed: a `CONECT` record may name an atom that is
    // listed further down the file, and the order of a bond is only known once
    // every record has been seen.
    const conect = new Map<string, PdbConect>();

    for (const line of lines) {
        const recType = line.substring(0, 6).trim();

        if (recType === 'ATOM' || recType === 'HETATM') {
            const serial = parseInt(line.substring(6, 11).trim(), 10);
            // A blank or hybrid-36 serial is not a number. `Set.has(NaN)` is
            // true after the first one, so adding NaN to the set would drop
            // every further atom with such a serial; only real numbers take
            // part in the duplicate check.
            if (!isNaN(serial)) {
                if (seenSerials.has(serial)) continue;
                seenSerials.add(serial);
            }

            const atomNameField = line.substring(12, 16);
            let altElement = '';
            if (line.length > 76) {
                const tail = line.substring(76).trim();
                if (tail.length > 0) {
                    const m = tail.match(/^([A-Za-z]{1,2})/);
                    if (m) altElement = m[1];
                }
            }
            if (!altElement && line.length > 12) {
                const tail = line.replace(/\s+$/, '');
                const lastWord = tail.split(/\s+/).pop() || '';
                if (/^[A-Za-z]{1,2}$/.test(lastWord)) {
                    altElement = lastWord;
                }
            }
            const element = pdbElement(atomNameField, altElement);

            const x = parseFloat(line.substring(30, 38).trim());
            const y = parseFloat(line.substring(38, 46).trim());
            const z = parseFloat(line.substring(46, 54).trim());

            if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                // Occupancy (columns 55-60) and the temperature-factor column
                // (61-66) are read as plain numbers: the visualization toolchain
                // writes an arbitrary per-atom scalar into the latter (unit
                // index, force error, …) rather than a B-factor, and the viewer
                // can colour the structure by it.
                const isDummy = element === 'X' || element === 'Xx' || element === 'XX' ||
                    element === 'Du' || element === 'DUM';
                const atom: Atom = { element: isDummy ? 'X' : element, x, y, z, index: atomIndex };
                // PDB has no ghost atom, so a dummy is written with the element
                // column `X` (the standard "unknown/dummy" symbol) — the same
                // spelling the viewer uses for a bare dummy point.
                if (isDummy) atom.ghost = true;
                const occupancy = parseFloat(line.substring(54, 60).trim());
                if (!isNaN(occupancy) && line.length > 59) {
                    atom.occupancy = occupancy;
                }
                const beta = parseFloat(line.substring(60, 66).trim());
                if (!isNaN(beta) && line.length > 65) {
                    atom.beta = beta;
                }
                if (!isNaN(serial)) atomSerialMap.set(serial, atomIndex);
                atoms.push(atom);
                atomIndex++;
            }
        } else if (recType === 'CONECT') {
            hasExplicitBonds = true;
            const fields: number[] = [];
            // CONECT carries up to four partners per line (columns 7-11, 12-16,
            // 17-21, 22-26) and a bond order is expressed by repeating a serial;
            // the scan therefore runs over every field the line has, not just
            // the four the format reserves.
            for (let pos = 6; pos < Math.min(line.length, 46); pos += 5) {
                const val = parseInt(line.substring(pos, pos + 5).trim(), 10);
                if (!isNaN(val)) fields.push(val);
            }
            if (fields.length >= 2) {
                const anchor = fields[0];
                for (let k = 1; k < fields.length; k++) {
                    const partner = fields[k];
                    if (partner === anchor) continue;
                    const lo = anchor < partner ? anchor : partner;
                    const hi = anchor < partner ? partner : anchor;
                    const key = lo + '_' + hi;
                    let pair = conect.get(key);
                    if (pair === undefined) {
                        pair = { lo: lo, hi: hi, loCount: 0, hiCount: 0 };
                        conect.set(key, pair);
                    }
                    if (anchor === lo) pair.loCount++;
                    else pair.hiCount++;
                }
            }
        } else if (recType === 'TITLE') {
            title += line.substring(10).trim() + ' ';
        } else if (recType === 'CRYST1') {
            if (cell === null && !cellBroken) {
                const read = readCryst1(line);
                if (read === null) cellBroken = true;
                else {
                    cell = read.cell;
                    cellSpaceGroup = read.spaceGroup;
                }
            }
        } else if (recType === 'HEADER' || recType === 'REMARK' || recType === 'COMPND' || recType === 'SOURCE' || recType === 'KEYWDS' || recType === 'EXPDTA' || recType === 'AUTHOR' || recType === 'REVDAT' || recType === 'JRNL' || recType === 'SEQRES' || recType === 'MODRES' || recType === 'DBREF' || recType === 'FORMUL' || recType === 'HELIX' || recType === 'SHEET' || recType === 'TURN' || recType === 'SSBOND' || recType === 'MASTER') {
            // skip
        }
    }

    // ----- the bonds the file stated -----
    // Emitted once every record has been read, so a record that names a later
    // atom still finds it. The order is the larger of the two directions'
    // counts: a channel written from both ends lists each single bond once on
    // each side, and summing would turn every one of them into a double bond.
    // Nothing is invented for a file that states no order — a bond listed once
    // stays a single bond, and that is exactly what the reader knows.
    conect.forEach(pair => {
        const i1 = atomSerialMap.get(pair.lo);
        const i2 = atomSerialMap.get(pair.hi);
        if (i1 === undefined || i2 === undefined || i1 === i2) return;
        let order = pair.loCount > pair.hiCount ? pair.loCount : pair.hiCount;
        if (order < 1) order = 1;
        if (order > 3) order = 3;
        bonds.push({ atom1: Math.min(i1, i2), atom2: Math.max(i1, i2), order });
    });

    title = title.trim() || 'PDB Structure';

    // ----- the cell the file states -----
    // A `CRYST1` record is the file's statement that the atoms below are a
    // unit cell, not a lone molecule: without it every periodic property —
    // periodic bonds, supercells, the cell editor, XRD — is unavailable and a
    // NaCl slab is drawn as 16 disconnected ions. The operations are the
    // identity, though: a PDB model is the coordinates the depositor wrote
    // (for a protein, the asymmetric unit), and expanding it under the stated
    // group would multiply every molecule in the file. The group is kept for
    // the cell panel, which is where the file states it.
    const warnings: string[] = [];
    if (cellBroken) {
        warnings.push('the CRYST1 record does not state a usable cell: the structure is read as a molecule');
    }
    if (cell !== null && atoms.length > 0) {
        const lattice = cellToLatticeVectors(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
        if (isPlaceholderCell(cell)) {
            warnings.push('the CRYST1 record states the 1 A placeholder cell some writers emit for a ' +
                'structure without one: the structure is read as a molecule');
        } else if (isUsableLattice(lattice)) {
            const sites: Site[] = atoms.map(atom => {
                const site: Site = { element: atom.element, cart: [atom.x, atom.y, atom.z] };
                if (atom.occupancy !== undefined) site.occupancy = atom.occupancy;
                if (atom.ghost) site.ghost = true;
                return site;
            });
            const data = buildPeriodicStructure(lattice, sites, title, warnings);
            if (data.crystal) {
                if (cellSpaceGroup !== '') data.crystal.spaceGroup = cellSpaceGroup;
                // `CONECT` bonds are the file's own table, exactly like a CIF's
                // `_geom_bond` loop: they address these cell atoms, so they
                // travel as `fileBonds` and are merged over the perceived set.
                if (bonds.length > 0) {
                    data.crystal.fileBonds = bonds;
                    data.hasExplicitBonds = true;
                }
            }
            return data;
        } else {
            warnings.push('the CRYST1 cell is degenerate: the structure is read as a molecule');
        }
    }

    const result: MolecularData = { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
    if (warnings.length > 0) result.warnings = warnings;
    return result;
}
