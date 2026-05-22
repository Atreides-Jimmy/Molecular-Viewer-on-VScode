"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parsePdb = parsePdb;
const ATOMIC_NUMBER_MAP = {
    H: 1, D: 1, HE: 2, LI: 3, BE: 4, B: 5, C: 6, N: 7, O: 8, F: 9, NE: 10,
    NA: 11, MG: 12, AL: 13, SI: 14, P: 15, S: 16, CL: 17, AR: 18, K: 19, CA: 20,
    SC: 21, TI: 22, V: 23, CR: 24, MN: 25, FE: 26, CO: 27, NI: 28, CU: 29, ZN: 30,
    GA: 31, GE: 32, AS: 33, SE: 34, BR: 35, KR: 36, RB: 37, SR: 38, Y: 39, ZR: 40,
    NB: 41, MO: 42, TC: 43, RU: 44, RH: 45, PD: 46, AG: 47, CD: 48, IN: 49, SN: 50,
    SB: 51, TE: 52, I: 53, XE: 54, CS: 55, BA: 56, LA: 57, CE: 58, PR: 59, ND: 60,
    PM: 61, SM: 62, EU: 63, GD: 64, TB: 65, DY: 66, HO: 67, ER: 68, TM: 69, YB: 70,
    LU: 71, HF: 72, TA: 73, W: 74, RE: 75, OS: 76, IR: 77, PT: 78, AU: 79, HG: 80,
    TL: 81, PB: 82, BI: 83, PO: 84, AT: 85, RN: 86
};
function pdbElement(atomName, altElement) {
    if (altElement && altElement.trim() !== '') {
        const e = altElement.trim();
        const upper = e.toUpperCase();
        if (ATOMIC_NUMBER_MAP[upper] !== undefined) {
            return upper.charAt(0) + upper.charAt(1).toLowerCase();
        }
        return e.charAt(0).toUpperCase() + e.slice(1).toLowerCase();
    }
    let name = atomName.trim();
    const match = name.match(/^([A-Za-z]+)/);
    if (!match)
        return 'C';
    let elem = match[1];
    if (elem.length >= 2) {
        const first = elem.charAt(0).toUpperCase();
        const second = elem.charAt(1).toLowerCase();
        const upper2 = elem.substring(0, 2).toUpperCase();
        if (ATOMIC_NUMBER_MAP[upper2] !== undefined) {
            return first + second;
        }
        if (ATOMIC_NUMBER_MAP[first] !== undefined) {
            return first;
        }
    }
    if (elem.length === 1) {
        const upper = elem.toUpperCase();
        if (ATOMIC_NUMBER_MAP[upper] !== undefined) {
            return upper;
        }
    }
    return elem.charAt(0).toUpperCase();
}
function parsePdb(content) {
    const lines = content.split(/\r?\n/);
    const atoms = [];
    const bonds = [];
    let hasExplicitBonds = false;
    let title = '';
    let charge;
    let multiplicity;
    let atomIndex = 0;
    const atomSerialMap = new Map();
    const seenSerials = new Set();
    for (const line of lines) {
        const recType = line.substring(0, 6).trim();
        if (recType === 'ATOM' || recType === 'HETATM') {
            const serial = parseInt(line.substring(6, 11).trim(), 10);
            if (seenSerials.has(serial))
                continue;
            seenSerials.add(serial);
            const atomName = line.substring(12, 16).trim();
            let altElement = '';
            if (line.length > 76) {
                const tail = line.substring(76).trim();
                if (tail.length > 0) {
                    const m = tail.match(/^([A-Za-z]{1,2})/);
                    if (m)
                        altElement = m[1];
                }
            }
            if (!altElement && line.length > 12) {
                const tail = line.replace(/\s+$/, '');
                const lastWord = tail.split(/\s+/).pop() || '';
                if (/^[A-Za-z]{1,2}$/.test(lastWord)) {
                    altElement = lastWord;
                }
            }
            const element = pdbElement(atomName, altElement);
            const x = parseFloat(line.substring(30, 38).trim());
            const y = parseFloat(line.substring(38, 46).trim());
            const z = parseFloat(line.substring(46, 54).trim());
            if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                atomSerialMap.set(serial, atomIndex);
                atoms.push({ element, x, y, z, index: atomIndex });
                atomIndex++;
            }
        }
        else if (recType === 'CONECT') {
            hasExplicitBonds = true;
            const fields = [];
            for (let pos = 6; pos < Math.min(line.length, 31); pos += 5) {
                const val = parseInt(line.substring(pos, pos + 5).trim(), 10);
                if (!isNaN(val))
                    fields.push(val);
            }
            if (fields.length >= 2) {
                const a1Serial = fields[0];
                for (let k = 1; k < fields.length; k++) {
                    const a2Serial = fields[k];
                    const a1 = atomSerialMap.get(a1Serial);
                    const a2 = atomSerialMap.get(a2Serial);
                    if (a1 !== undefined && a2 !== undefined && a1 < a2) {
                        const exists = bonds.some(b => (b.atom1 === a1 && b.atom2 === a2) || (b.atom1 === a2 && b.atom2 === a1));
                        if (!exists)
                            bonds.push({ atom1: a1, atom2: a2, order: 1 });
                    }
                }
            }
        }
        else if (recType === 'TITLE') {
            title += line.substring(10).trim() + ' ';
        }
        else if (recType === 'CRYST1' || recType === 'HEADER' || recType === 'REMARK' || recType === 'COMPND' || recType === 'SOURCE' || recType === 'KEYWDS' || recType === 'EXPDTA' || recType === 'AUTHOR' || recType === 'REVDAT' || recType === 'JRNL' || recType === 'SEQRES' || recType === 'MODRES' || recType === 'DBREF' || recType === 'FORMUL' || recType === 'HELIX' || recType === 'SHEET' || recType === 'TURN' || recType === 'SSBOND' || recType === 'MASTER') {
            // skip
        }
    }
    title = title.trim() || 'PDB Structure';
    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
//# sourceMappingURL=pdbParser.js.map