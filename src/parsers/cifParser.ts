import { Atom, CrystalData, MolecularData } from '../types';

interface CifLoop {
    keys: string[];
    rows: string[][];
}

interface CifData {
    items: { [key: string]: string };
    loops: CifLoop[];
}

function stripEsd(token: string): string {
    return token.replace(/\(\d+\)$/, '');
}

function parseNumber(token: string): number | null {
    if (token === '?' || token === '.' || token === '') return null;
    const cleaned = stripEsd(token).replace(/[()]/g, '');
    const v = parseFloat(cleaned);
    return isNaN(v) ? null : v;
}

function tokenize(content: string): CifData {
    const data: CifData = { items: {}, loops: [] };
    const lines = content.split(/\r?\n/);
    let i = 0;
    let inSemicolon = false;
    let semicolonBuffer: string[] = [];
    let semicolonKey: string | null = null;

    function flushSemicolon() {
        if (semicolonKey !== null) {
            data.items[semicolonKey] = semicolonBuffer.join('\n').trim();
        }
        inSemicolon = false;
        semicolonBuffer = [];
        semicolonKey = null;
    }

    function readValue(rest: string): { value: string; consumed: number } {
        const trimmed = rest.trimStart();
        if (trimmed.startsWith("'")) {
            const end = trimmed.indexOf("'", 1);
            if (end >= 0) return { value: trimmed.substring(1, end), consumed: (rest.length - trimmed.length) + end + 1 };
        }
        if (trimmed.startsWith('"')) {
            const end = trimmed.indexOf('"', 1);
            if (end >= 0) return { value: trimmed.substring(1, end), consumed: (rest.length - trimmed.length) + end + 1 };
        }
        const parts = trimmed.split(/\s+/);
        return { value: parts[0], consumed: rest.length - trimmed.length + parts[0].length };
    }

    while (i < lines.length) {
        const rawLine = lines[i];
        const line = rawLine.trim();

        if (inSemicolon) {
            if (line.startsWith(';')) {
                flushSemicolon();
                i++;
                continue;
            }
            semicolonBuffer.push(rawLine);
            i++;
            continue;
        }

        if (line === '' || line.startsWith('#')) {
            i++;
            continue;
        }

        if (line.startsWith(';')) {
            inSemicolon = true;
            semicolonBuffer = [];
            const after = line.substring(1).trim();
            if (after !== '') semicolonBuffer.push(after);
            i++;
            continue;
        }

        if (line.startsWith('data_')) {
            i++;
            continue;
        }

        if (line.startsWith('save_')) {
            i++;
            continue;
        }

        if (line.startsWith('loop_')) {
            const loop: CifLoop = { keys: [], rows: [] };
            i++;
            let startedRows = false;
            while (i < lines.length) {
                const l = lines[i].trim();
                if (l === '' || l.startsWith('#')) { i++; continue; }
                if (l.startsWith('_')) {
                    if (startedRows) break;
                    const parts = l.split(/\s+/);
                    if (parts.length > 1) {
                        // Key with inline value = scalar item, not part of loop
                        break;
                    }
                    loop.keys.push(parts[0]);
                    i++;
                } else if (l.startsWith('loop_') || l.startsWith('data_') || l.startsWith('save_')) {
                    break;
                } else {
                    startedRows = true;
                    const row: string[] = [];
                    let lineRest = lines[i];
                    for (let k = 0; k < loop.keys.length; k++) {
                        if (lineRest.trim() === '') break;
                        if (lineRest.trim().startsWith(';')) {
                            const buf: string[] = [];
                            const firstPart = lineRest.trim().substring(1);
                            if (firstPart !== '') buf.push(firstPart);
                            i++;
                            while (i < lines.length && !lines[i].trim().startsWith(';')) {
                                buf.push(lines[i]);
                                i++;
                            }
                            row.push(buf.join('\n').trim());
                            lineRest = '';
                            if (i < lines.length) {
                                lineRest = lines[i].replace(/^\s*;\s*/, '');
                            }
                        } else {
                            const { value, consumed } = readValue(lineRest);
                            row.push(value);
                            lineRest = lineRest.substring(consumed);
                        }
                    }
                    if (row.length > 0) loop.rows.push(row);
                    i++;
                }
            }
            if (loop.keys.length > 0) {
                if (loop.rows.length === 0 || loop.rows[0].length === loop.keys.length) {
                    data.loops.push(loop);
                }
            }
            continue;
        }

        if (line.startsWith('_')) {
            const match = line.match(/^(_\S+)\s+(.*)$/);
            if (match) {
                let value = match[2].trim();
                if (value.startsWith(';')) {
                    inSemicolon = true;
                    semicolonKey = match[1];
                    semicolonBuffer = [];
                    const after = value.substring(1);
                    if (after.trim() !== '') semicolonBuffer.push(after);
                    i++;
                    continue;
                }
                if (value.startsWith("'") && value.endsWith("'") && value.length >= 2) {
                    value = value.substring(1, value.length - 1);
                } else if (value.startsWith('"') && value.endsWith('"') && value.length >= 2) {
                    value = value.substring(1, value.length - 1);
                } else {
                    value = value.split(/\s+/)[0];
                }
                data.items[match[1]] = value;
            } else {
                const key = line.split(/\s+/)[0];
                data.items[key] = '';
            }
            i++;
            continue;
        }

        i++;
    }

    flushSemicolon();
    return data;
}

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

interface SymOp {
    R: number[][];
    t: number[];
}

function parseSymOp(s: string): SymOp {
    const parts = s.split(',').map(p => p.trim());
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const t = [0, 0, 0];

    for (let idx = 0; idx < 3; idx++) {
        if (idx >= parts.length) continue;
        const expr = parts[idx].replace(/\s+/g, '').toLowerCase();
        let sign = 1;
        let k = 0;
        while (k < expr.length) {
            const ch = expr[k];
            if (ch === '+') { sign = 1; k++; continue; }
            if (ch === '-') { sign = -1; k++; continue; }
            if (ch === 'x' || ch === 'y' || ch === 'z') {
                const col = ch === 'x' ? 0 : ch === 'y' ? 1 : 2;
                R[idx][col] += sign;
                sign = 1;
                k++;
                continue;
            }
            const numMatch = expr.substring(k).match(/^(\d+\/\d+|\d+\.?\d*)/);
            if (numMatch) {
                let v: number;
                if (numMatch[1].includes('/')) {
                    const [p, q] = numMatch[1].split('/').map(Number);
                    v = p / q;
                } else {
                    v = parseFloat(numMatch[1]);
                }
                if (k > 0 && (expr[k - 1] === 'x' || expr[k - 1] === 'y' || expr[k - 1] === 'z')) {
                    t[idx] += sign * v;
                } else {
                    t[idx] += sign * v;
                }
                sign = 1;
                k += numMatch[1].length;
                continue;
            }
            k++;
        }
    }

    return { R, t };
}

function applySymOp(op: SymOp, frac: number[]): number[] {
    const result = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        result[i] = op.R[i][0] * frac[0] + op.R[i][1] * frac[1] + op.R[i][2] * frac[2] + op.t[i];
    }
    return result;
}

function wrapFrac(v: number): number {
    let r = v - Math.floor(v);
    if (r >= 1 - 1e-6) r -= 1;
    if (r < 0) r += 1;
    return r;
}

function findLoop(data: CifData, keySubstr: string): CifLoop | null {
    for (const loop of data.loops) {
        for (const key of loop.keys) {
            if (key.includes(keySubstr)) return loop;
        }
    }
    return null;
}

function getLoopColumn(loop: CifLoop, keySubstr: string): number {
    for (let i = 0; i < loop.keys.length; i++) {
        if (loop.keys[i].includes(keySubstr)) return i;
    }
    return -1;
}

export function parseCif(content: string): MolecularData {
    const data = tokenize(content);

    const a = parseNumber(data.items['_cell_length_a']) || 1;
    const b = parseNumber(data.items['_cell_length_b']) || 1;
    const c = parseNumber(data.items['_cell_length_c']) || 1;
    const alpha = parseNumber(data.items['_cell_angle_alpha']) ?? 90;
    const beta = parseNumber(data.items['_cell_angle_beta']) ?? 90;
    const gamma = parseNumber(data.items['_cell_angle_gamma']) ?? 90;
    const spaceGroup = data.items['_symmetry_space_group_name_H-M'] || undefined;

    const lattice = cellToLatticeVectors(a, b, c, alpha, beta, gamma);

    const symmetryOps: string[] = [];
    const symLoop = findLoop(data, '_symmetry_equiv_pos_as_xyz');
    if (symLoop) {
        const col = getLoopColumn(symLoop, '_symmetry_equiv_pos_as_xyz');
        if (col >= 0) {
            for (const row of symLoop.rows) {
                if (row[col]) symmetryOps.push(row[col]);
            }
        }
    }
    if (symmetryOps.length === 0) {
        symmetryOps.push('x, y, z');
    }

    const parsedSymOps = symmetryOps.map(parseSymOp);

    const atomLoop = findLoop(data, 'fract_x') || findLoop(data, 'Cartn_x');
    if (!atomLoop) {
        const empty: MolecularData = {
            atoms: [], bonds: [], title: 'No atom sites found in CIF',
            hasExplicitBonds: false,
            crystal: {
                a, b, c, alpha, beta, gamma,
                latticeVectors: lattice, spaceGroup, symmetryOps,
                baseAtoms: [], baseBonds: []
            }
        };
        return empty;
    }

    const useFract = getLoopColumn(atomLoop, 'fract_x') >= 0;
    const labelCol = getLoopColumn(atomLoop, 'label');
    const typeCol = getLoopColumn(atomLoop, 'type_symbol');
    const xCol = getLoopColumn(atomLoop, useFract ? 'fract_x' : 'Cartn_x');
    const yCol = getLoopColumn(atomLoop, useFract ? 'fract_y' : 'Cartn_y');
    const zCol = getLoopColumn(atomLoop, useFract ? 'fract_z' : 'Cartn_z');
    const occCol = getLoopColumn(atomLoop, 'occupancy');

    interface RawAtom {
        element: string;
        label: string;
        frac: number[];
        occupancy: number;
    }
    const rawAtoms: RawAtom[] = [];

    const cols = [labelCol, typeCol, xCol, yCol, zCol].filter(c => c >= 0);
    const minCols = cols.length > 0 ? Math.max(...cols) + 1 : 0;

    for (const row of atomLoop.rows) {
        if (row.length < minCols) continue;
        const label = labelCol >= 0 ? row[labelCol] : '';
        const typeSymRaw = typeCol >= 0 ? row[typeCol] : (label ? label.replace(/[0-9].*$/, '') : '');
        const elementRaw = typeSymRaw || (label ? label.replace(/[0-9].*$/, '') : 'C');
        const element = elementRaw.charAt(0).toUpperCase() + elementRaw.slice(1).toLowerCase();
        const xv = parseNumber(row[xCol]);
        const yv = parseNumber(row[yCol]);
        const zv = parseNumber(row[zCol]);
        if (xv === null || yv === null || zv === null) continue;
        const occupancy = occCol >= 0 ? (parseNumber(row[occCol]) ?? 1) : 1;

        const frac = useFract ? [xv, yv, zv] :
            cartToFrac([xv, yv, zv], lattice);
        rawAtoms.push({ element, label, frac, occupancy });
    }

    const baseAtoms: Atom[] = [];
    let atomIndex = 0;
    for (const raw of rawAtoms) {
        for (const op of parsedSymOps) {
            const transformed = applySymOp(op, raw.frac);
            const wrapped = transformed.map(wrapFrac);
            const [cx, cy, cz] = fracToCart(wrapped, lattice);

            let duplicate = false;
            for (const existing of baseAtoms) {
                const dx = existing.x - cx;
                const dy = existing.y - cy;
                const dz = existing.z - cz;
                if (dx * dx + dy * dy + dz * dz < 0.25) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) continue;

            baseAtoms.push({
                element: raw.element,
                x: cx, y: cy, z: cz,
                index: atomIndex++,
                occupancy: raw.occupancy,
                baseIdx: baseAtoms.length
            });
        }
    }

    const atoms: Atom[] = baseAtoms.map(a => ({
        element: a.element,
        x: a.x, y: a.y, z: a.z,
        index: a.index,
        occupancy: a.occupancy,
        baseIdx: a.index,
        cellI: 0, cellJ: 0, cellK: 0
    }));

    const crystal: CrystalData = {
        a, b, c, alpha, beta, gamma,
        latticeVectors: lattice,
        spaceGroup,
        symmetryOps,
        baseAtoms: baseAtoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: []
    };

    const title = data.items['_chemical_name_common'] ||
                  data.items['_chemical_formula_sum'] ||
                  'CIF Structure';

    return {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        crystal
    };
}

function cartToFrac(cart: number[], lattice: number[][]): number[] {
    const [aVec, bVec, cVec] = lattice;
    const det = aVec[0] * (bVec[1] * cVec[2] - bVec[2] * cVec[1]) -
                aVec[1] * (bVec[0] * cVec[2] - bVec[2] * cVec[0]) +
                aVec[2] * (bVec[0] * cVec[1] - bVec[1] * cVec[0]);
    if (Math.abs(det) < 1e-10) return [0, 0, 0];
    const invDet = 1 / det;
    const inv = [
        [
            (bVec[1] * cVec[2] - bVec[2] * cVec[1]) * invDet,
            (aVec[2] * cVec[1] - aVec[1] * cVec[2]) * invDet,
            (aVec[1] * bVec[2] - aVec[2] * bVec[1]) * invDet
        ],
        [
            (bVec[2] * cVec[0] - bVec[0] * cVec[2]) * invDet,
            (aVec[0] * cVec[2] - aVec[2] * cVec[0]) * invDet,
            (aVec[2] * bVec[0] - aVec[0] * bVec[2]) * invDet
        ],
        [
            (bVec[0] * cVec[1] - bVec[1] * cVec[0]) * invDet,
            (aVec[1] * cVec[0] - aVec[0] * cVec[1]) * invDet,
            (aVec[0] * bVec[1] - aVec[1] * bVec[0]) * invDet
        ]
    ];
    return [
        inv[0][0] * cart[0] + inv[0][1] * cart[1] + inv[0][2] * cart[2],
        inv[1][0] * cart[0] + inv[1][1] * cart[1] + inv[1][2] * cart[2],
        inv[2][0] * cart[0] + inv[2][1] * cart[1] + inv[2][2] * cart[2]
    ];
}

export function cartesianToFractional(cart: number[], lattice: number[][]): number[] {
    return cartToFrac(cart, lattice);
}

export function fractionalToCartesian(frac: number[], lattice: number[][]): [number, number, number] {
    return fracToCart(frac, lattice);
}
