import { Atom, Bond, MolecularData } from '../types';
import { cellToLatticeVectors } from '../tools/crystalGeometry';
import { lookupSpaceGroupOps } from '../tools/spaceGroupTable';
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
 * One `ATOM`/`HETATM` record exactly as the file wrote it.
 *
 * Records are collected before any of them becomes an atom, because the two
 * decisions that discard some of them need the whole file: an alternate
 * conformation is a *group* of records of which one survives, and a repeated
 * serial is only a duplicate when the record beside it is the same atom. Doing
 * either while the records stream past would decide on half the evidence.
 */
interface PdbAtomRecord {
    /** Position in the file, so the survivors can be put back into file order. */
    order: number;
    serial: number;
    /** Atom-name field, columns 13-16, alignment included. */
    name: string;
    /** Alternate-location indicator, column 17; empty for the primary conformer. */
    altLoc: string;
    resName: string;
    chain: string;
    resSeq: string;
    iCode: string;
    element: string;
    x: number;
    y: number;
    z: number;
    occupancy?: number;
    beta?: number;
}

/** A `SSBOND` record: the two cysteine residues a disulfide joins. */
interface PdbSsBond {
    chain1: string;
    seq1: string;
    iCode1: string;
    sym1: string;
    chain2: string;
    seq2: string;
    iCode2: string;
    sym2: string;
}

/**
 * One `SSBOND` record, which names its two cysteines by residue rather than by
 * atom: chain, sequence number, insertion code and residue name.
 *
 * The record's symmetry operators say whether the partner is the residue as
 * listed (the identity operator, `555`) or one of its symmetry images. A bond
 * to an image is not a bond this reader can draw — the page's molecule graph has
 * no cell to place that image in — so those records are counted and reported
 * instead of being turned into a bond to the wrong atom.
 *
 * The operators are read as the tokens before the S-S distance rather than from
 * columns 44-46 and 51-53. The blocks in the wild do not line up with those
 * columns — the format's own printed example puts them twenty columns further
 * right — and a column that picks up the wrong characters turns a real
 * disulfide into a "symmetry bond" that is then thrown away.
 */
function readSsBond(line: string): PdbSsBond {
    const tail = line.substring(36).trim().split(/\s+/).filter(t => t !== '');
    // The record's last field is the S-S distance in Angstrom.
    if (tail.length > 0 && /^[0-9.]+$/.test(tail[tail.length - 1])) tail.pop();
    const sym = tail.slice(-2);
    return {
        chain1: line.substring(15, 16).trim(),
        seq1: line.substring(17, 21).trim(),
        iCode1: line.substring(21, 22).trim(),
        sym1: sym.length > 0 ? sym[0] : '',
        chain2: line.substring(29, 30).trim(),
        seq2: line.substring(31, 35).trim(),
        iCode2: line.substring(35, 36).trim(),
        sym2: sym.length > 1 ? sym[1] : ''
    };
}

/**
 * Whether a symmetry code means "the residue as it is listed".
 *
 * The code is three digits — an operation number and two translations — and the
 * identity is `555`, spelled `1555` when the operation is given as well: no
 * rotation, no translation. Anything with another digit names a real image.
 */
function isIdentitySymmetryCode(code: string): boolean {
    if (code === '') return true;
    return /^[15]+$/.test(code);
}

/**
 * Which record of an alternate-conformation group survives.
 *
 * The occupancies the file states decide, because that is what the depositor
 * refined against; a tie goes to the blank altLoc (the record the file marked
 * as the primary conformer) and then to the one written first. A group with a
 * single altLoc — the usual case, every atom blank — never reaches here.
 */
function chooseAltLoc(list: PdbAtomRecord[]): string {
    const occupancyOf = (r: PdbAtomRecord) => r.occupancy !== undefined ? r.occupancy : 1;
    let best = list[0];
    for (const rec of list) {
        const occ = occupancyOf(rec);
        const bestOcc = occupancyOf(best);
        if (occ > bestOcc + 1e-9) {
            best = rec;
        } else if (Math.abs(occ - bestOcc) <= 1e-9 && best.altLoc !== '' && rec.altLoc === '') {
            best = rec;
        }
    }
    return best.altLoc;
}

/**
 * Reduce every group of alternate conformations to one conformer.
 *
 * A file refined with alternate conformations lists the same atom two or three
 * times (`CA` as altLoc `A` and as altLoc `B`), each with its own occupancy and
 * position. Read as-is, every such atom appeared twice in the structure, the
 * bond table drew both copies, and a formula or an electron count doubled. The
 * group is identified the way the format identifies it — chain, residue
 * number, insertion code and atom name — and only the chosen conformer's
 * records survive.
 */
function selectAlternateConformations(records: PdbAtomRecord[]): {
    kept: PdbAtomRecord[];
    dropped: number;
    sites: number;
} {
    const groups = new Map<string, PdbAtomRecord[]>();
    const order: string[] = [];
    for (const rec of records) {
        // The residue name is deliberately not part of the key: alternate
        // *residues* (a serine and a cysteine refined in the same place) share
        // the backbone atoms' names and must collapse the same way.
        const key = rec.chain + '|' + rec.resSeq + '|' + rec.iCode + '|' + rec.name;
        const list = groups.get(key);
        if (list) list.push(rec);
        else {
            groups.set(key, [rec]);
            order.push(key);
        }
    }

    const kept: PdbAtomRecord[] = [];
    let dropped = 0;
    let sites = 0;
    for (const key of order) {
        const list = groups.get(key) as PdbAtomRecord[];
        const distinct: { [alt: string]: true } = {};
        for (const rec of list) distinct[rec.altLoc] = true;
        if (Object.keys(distinct).length <= 1) {
            for (const rec of list) kept.push(rec);
            continue;
        }
        sites++;
        const chosen = chooseAltLoc(list);
        for (const rec of list) {
            if (rec.altLoc === chosen) kept.push(rec);
            else dropped++;
        }
    }
    // Every index this reader produces is file order, so the survivors go back
    // into it before anything else looks at them.
    kept.sort((a, b) => a.order - b.order);
    return { kept, dropped, sites };
}

/** Whether two records describe the very same atom, serial aside. */
function isSameRecord(a: PdbAtomRecord, b: PdbAtomRecord): boolean {
    return a.name === b.name && a.element === b.element && a.chain === b.chain &&
        a.resSeq === b.resSeq && a.iCode === b.iCode &&
        Math.abs(a.x - b.x) < 1e-3 && Math.abs(a.y - b.y) < 1e-3 &&
        Math.abs(a.z - b.z) < 1e-3;
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
        const cell: CellParameters = {
            a: tokenised[0], b: tokenised[1], c: tokenised[2],
            alpha: tokenised[3], beta: tokenised[4], gamma: tokenised[5]
        };
        spaceGroup = spaceGroupFromTokens(parts.slice(7), cell);
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
 * The space group a collapsed `CRYST1` states, with its trailing Z removed.
 *
 * The record ends with Z, the number of formula units, so a bare last token
 * looks like a value to drop — but a symbol carries numbers of its own
 * (`P 21 21 21`, `P 1 21/c 1`), and dropping the last one turned the
 * orthorhombic symbol into `P 21 21`, which names no group at all and loses the
 * setting every later tool would need. The built-in table settles it: the
 * trailing token is read as Z only when removing it leaves a symbol the table
 * knows and the untrimmed line is not one.
 */
function spaceGroupFromTokens(rest: string[], cell: CellParameters): string {
    if (rest.length === 0) return '';
    const full = rest.join(' ');
    if (rest.length < 2) return full;
    if (!/^\d+$/.test(rest[rest.length - 1])) return full;
    const trimmed = rest.slice(0, rest.length - 1).join(' ');
    if (!/[A-Za-z]/.test(trimmed)) return full;
    if (lookupSpaceGroupOps(full, cell) !== null) return full;
    if (lookupSpaceGroupOps(trimmed, cell) !== null) return trimmed;
    // Neither reading is a group this table knows. Keeping every token is the
    // safer failure: a stray Z in a caption is cosmetic, a symbol that lost its
    // setting is not.
    return full;
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
    const records: PdbAtomRecord[] = [];
    const conect = new Map<string, PdbConect>();
    const ssBondRecords: PdbSsBond[] = [];
    let hasExplicitBonds = false;
    let title = '';
    const charge: number | undefined = undefined;
    const multiplicity: number | undefined = undefined;
    let recordOrder = 0;
    // `MODEL` counts the models a deposition holds. An NMR ensemble writes each
    // model's atoms under the same serials, so reading them all folded every
    // conformer of the ensemble into one structure; reading them as separate
    // frames needs the trajectory machinery, which a PDB does not feed. Model 1
    // is the structure, and saying so is the rest.
    let modelCount = 0;
    let currentModel = 0;
    let atomsPastModel1 = 0;

    // The first `CRYST1` record wins: a file that repeats it (two jobs
    // concatenated) describes the same cell both times, and the atoms are the
    // ones written under the first. A record that is present but unreadable is
    // remembered so the reader can say so rather than silently drawing a
    // molecule.
    let cell: CellParameters | null = null;
    let cellSpaceGroup = '';
    let cellBroken = false;

    for (const line of lines) {
        const recType = line.substring(0, 6).trim();

        if (recType === 'ATOM' || recType === 'HETATM') {
            if (currentModel > 1) {
                atomsPastModel1++;
                continue;
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
            if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

            const rec: PdbAtomRecord = {
                order: recordOrder++,
                serial: parseInt(line.substring(6, 11).trim(), 10),
                name: atomNameField,
                altLoc: line.length > 16 ? line.substring(16, 17).trim() : '',
                resName: line.substring(17, 20).trim(),
                chain: line.substring(21, 22).trim(),
                resSeq: line.substring(22, 26).trim(),
                iCode: line.substring(26, 27).trim(),
                element,
                x, y, z
            };
            // Occupancy (columns 55-60) and the temperature-factor column
            // (61-66) are read as plain numbers: the visualization toolchain
            // writes an arbitrary per-atom scalar into the latter (unit
            // index, force error, …) rather than a B-factor, and the viewer
            // can colour the structure by it.
            const occupancy = parseFloat(line.substring(54, 60).trim());
            if (!isNaN(occupancy) && line.length > 59) rec.occupancy = occupancy;
            const beta = parseFloat(line.substring(60, 66).trim());
            if (!isNaN(beta) && line.length > 65) rec.beta = beta;
            records.push(rec);
        } else if (recType === 'MODEL') {
            modelCount++;
            currentModel++;
        } else if (recType === 'ENDMDL') {
            // Nothing to close: the model number stays as it is, so the atoms
            // that follow without a further `MODEL` still belong to it.
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
        } else if (recType === 'SSBOND') {
            ssBondRecords.push(readSsBond(line));
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
        } else if (recType === 'HEADER' || recType === 'REMARK' || recType === 'COMPND' || recType === 'SOURCE' || recType === 'KEYWDS' || recType === 'EXPDTA' || recType === 'AUTHOR' || recType === 'REVDAT' || recType === 'JRNL' || recType === 'SEQRES' || recType === 'MODRES' || recType === 'DBREF' || recType === 'FORMUL' || recType === 'HELIX' || recType === 'SHEET' || recType === 'TURN' || recType === 'TER' || recType === 'MASTER') {
            // skip. `TER` is the chain break and carries no atom of its own; it
            // is named here because it is also where a writer is free to
            // restart its atom numbering, which the serial check below allows
            // for instead of dropping the atoms after the restart.
        }
    }

    const warnings: string[] = [];
    if (modelCount > 1) {
        warnings.push('the file holds ' + modelCount + ' models (an NMR ensemble or a multi-model ' +
            'deposition): model 1 is shown, ' + atomsPastModel1 + ' atom record' +
            (atomsPastModel1 === 1 ? '' : 's') + ' of the later models are not');
    }

    // ----- alternate conformations, then duplicate serials -----
    const alt = selectAlternateConformations(records);
    const atoms: Atom[] = [];
    const atomSerialMap = new Map<number, number>();
    // Which atom indices share a cysteine residue, so an `SSBOND` record can
    // name the two sulfur atoms it joins.
    const sulphurByResidue = new Map<string, number[]>();
    const firstOfSerial = new Map<number, PdbAtomRecord>();
    let duplicateRecords = 0;
    let repeatedSerials = 0;
    let atomIndex = 0;

    for (const rec of alt.kept) {
        // A blank or hybrid-36 serial is not a number. `Set.has(NaN)` is true
        // after the first one, so adding NaN to the set would drop every
        // further atom with such a serial; only real numbers take part in the
        // duplicate check.
        if (!isNaN(rec.serial)) {
            const first = firstOfSerial.get(rec.serial);
            if (first === undefined) {
                firstOfSerial.set(rec.serial, rec);
            } else if (isSameRecord(first, rec)) {
                // The same atom written twice — a duplicated record, not a
                // numbering restart.
                duplicateRecords++;
                continue;
            } else {
                // A serial used twice for different atoms: the numbering was
                // restarted (per chain, per model, or by concatenating two
                // files). The atom is real and is kept; only the bond table
                // has to pick one of them, and it picks the first.
                repeatedSerials++;
            }
        }

        const isDummy = rec.element === 'X' || rec.element === 'Xx' || rec.element === 'XX' ||
            rec.element === 'Du' || rec.element === 'DUM';
        const atom: Atom = {
            element: isDummy ? 'X' : rec.element,
            x: rec.x, y: rec.y, z: rec.z,
            index: atomIndex
        };
        // PDB has no ghost atom, so a dummy is written with the element
        // column `X` (the standard "unknown/dummy" symbol) — the same
        // spelling the viewer uses for a bare dummy point.
        if (isDummy) atom.ghost = true;
        if (rec.occupancy !== undefined) atom.occupancy = rec.occupancy;
        if (rec.beta !== undefined) atom.beta = rec.beta;
        if (!isNaN(rec.serial) && !atomSerialMap.has(rec.serial)) atomSerialMap.set(rec.serial, atomIndex);
        if (rec.element === 'S') {
            const key = rec.chain + '|' + rec.resSeq + '|' + rec.iCode;
            const list = sulphurByResidue.get(key);
            // A disulfide joins the two SG atoms, so one named SG comes first:
            // a residue with any other sulfur still resolves to its own SG.
            if (list === undefined) sulphurByResidue.set(key, [atomIndex]);
            else if (/^\s*SG/i.test(rec.name)) list.unshift(atomIndex);
            else list.push(atomIndex);
        }
        atoms.push(atom);
        atomIndex++;
    }

    if (alt.dropped > 0) {
        warnings.push('alternate conformations: ' + alt.dropped + ' atom record' +
            (alt.dropped === 1 ? '' : 's') + ' at ' + alt.sites + ' site' + (alt.sites === 1 ? '' : 's') +
            ' were dropped, keeping the highest-occupancy conformer of each');
    }
    if (duplicateRecords > 0) {
        warnings.push(duplicateRecords + ' atom record' + (duplicateRecords === 1 ? '' : 's') +
            ' repeat an earlier record exactly (same serial, name and position): skipped');
    }
    if (repeatedSerials > 0) {
        warnings.push(repeatedSerials + ' atom record' + (repeatedSerials === 1 ? '' : 's') +
            ' reuse a serial number for a different atom (the numbering restarted): every one is kept, ' +
            'and CONECT records bond the first of each serial');
    }

    // ----- the bonds the file stated -----
    // Emitted once every record has been read, so a record that names a later
    // atom still finds it. The order is the larger of the two directions'
    // counts: a channel written from both ends lists each single bond once on
    // each side, and summing would turn every one of them into a double bond.
    // Nothing is invented for a file that states no order — a bond listed once
    // stays a single bond, and that is exactly what the reader knows.
    const conectBonds: Bond[] = [];
    conect.forEach(pair => {
        const i1 = atomSerialMap.get(pair.lo);
        const i2 = atomSerialMap.get(pair.hi);
        if (i1 === undefined || i2 === undefined || i1 === i2) return;
        let order = pair.loCount > pair.hiCount ? pair.loCount : pair.hiCount;
        if (order < 1) order = 1;
        if (order > 3) order = 3;
        conectBonds.push({ atom1: Math.min(i1, i2), atom2: Math.max(i1, i2), order });
    });

    // ----- the disulfides the file stated -----
    // `SSBOND` names two cysteines by residue; the bond is between their SG
    // atoms, which is what the record's own length column measures. It is a
    // declaration beside the `CONECT` table rather than part of it, so it is
    // kept apart and merged over the perceived bonds — a file with SSBOND and
    // no CONECT still gets geometry perception for everything else.
    const statedBonds: Bond[] = [];
    let ssbondUnresolved = 0;
    let ssbondSymmetry = 0;
    for (const ss of ssBondRecords) {
        if (!isIdentitySymmetryCode(ss.sym1) || !isIdentitySymmetryCode(ss.sym2)) {
            ssbondSymmetry++;
            continue;
        }
        const findSulphur = (chain: string, seq: string, iCode: string): number | undefined => {
            const list = sulphurByResidue.get(chain + '|' + seq + '|' + iCode);
            if (list === undefined || list.length === 0) return undefined;
            return list[0];
        };
        const i1 = findSulphur(ss.chain1, ss.seq1, ss.iCode1);
        const i2 = findSulphur(ss.chain2, ss.seq2, ss.iCode2);
        if (i1 === undefined || i2 === undefined || i1 === i2) {
            ssbondUnresolved++;
            continue;
        }
        statedBonds.push({ atom1: Math.min(i1, i2), atom2: Math.max(i1, i2), order: 1 });
    }
    if (ssbondUnresolved > 0) {
        warnings.push('SSBOND: ' + ssbondUnresolved + ' record' + (ssbondUnresolved === 1 ? '' : 's') +
            ' name a residue whose sulfur atom is not in the atom list: no bond is drawn for ' +
            (ssbondUnresolved === 1 ? 'it' : 'them'));
    }
    if (ssbondSymmetry > 0) {
        warnings.push('SSBOND: ' + ssbondSymmetry + ' record' + (ssbondSymmetry === 1 ? '' : 's') +
            ' join a residue to one of its symmetry images: that bond is not part of the molecule ' +
            'graph and is left out');
    }

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
                // An `SSBOND` record states a bond of the same kind and joins
                // the same table.
                const fileBonds = conectBonds.concat(statedBonds);
                if (fileBonds.length > 0) {
                    data.crystal.fileBonds = fileBonds;
                    data.hasExplicitBonds = conectBonds.length > 0;
                }
            }
            return data;
        } else {
            warnings.push('the CRYST1 cell is degenerate: the structure is read as a molecule');
        }
    }

    const result: MolecularData = { atoms, bonds: conectBonds, title, hasExplicitBonds, charge, multiplicity };
    // A molecule's bonds are perceived from its geometry; an SSBOND record is
    // the one bond the file states beside that, so it travels as a stated bond
    // for the loader to merge over the perceived set instead of replacing it.
    if (statedBonds.length > 0) result.statedBonds = statedBonds;
    if (warnings.length > 0) result.warnings = warnings;
    return result;
}
