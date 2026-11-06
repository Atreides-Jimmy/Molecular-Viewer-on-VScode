import { Atom, CrystalData, MolecularData } from '../types';
import { normalizeElementSymbol } from '../tools/periodicTable';
import { isUsableLattice } from './periodicUtils';

function latticeToCellParams(lattice: number[][]): { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } {
    const [aVec, bVec, cVec] = lattice;
    const a = Math.sqrt(aVec[0] * aVec[0] + aVec[1] * aVec[1] + aVec[2] * aVec[2]);
    const b = Math.sqrt(bVec[0] * bVec[0] + bVec[1] * bVec[1] + bVec[2] * bVec[2]);
    const c = Math.sqrt(cVec[0] * cVec[0] + cVec[1] * cVec[1] + cVec[2] * cVec[2]);

    const cosA = (bVec[0] * cVec[0] + bVec[1] * cVec[1] + bVec[2] * cVec[2]) / (b * c);
    const cosB = (aVec[0] * cVec[0] + aVec[1] * cVec[1] + aVec[2] * cVec[2]) / (a * c);
    const cosG = (aVec[0] * bVec[0] + aVec[1] * bVec[1] + aVec[2] * bVec[2]) / (a * b);

    const rad = 180 / Math.PI;
    return {
        a, b, c,
        alpha: Math.acos(Math.max(-1, Math.min(1, cosA))) * rad,
        beta: Math.acos(Math.max(-1, Math.min(1, cosB))) * rad,
        gamma: Math.acos(Math.max(-1, Math.min(1, cosG))) * rad
    };
}

function fracToCart(frac: number[], lattice: number[][]): [number, number, number] {
    const x = frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0];
    const y = frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1];
    const z = frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2];
    return [x, y, z];
}

function crystalDataFor(lattice: number[][], atoms: Atom[]): CrystalData {
    return {
        ...latticeToCellParams(lattice),
        latticeVectors: lattice,
        symmetryOps: ['x, y, z'],
        baseAtoms: atoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: []
    };
}

function molecularDataFor(atoms: Atom[], title: string, warnings: string[]): MolecularData {
    const data: MolecularData = { atoms, bonds: [], title, hasExplicitBonds: false };
    if (warnings.length > 0) data.warnings = warnings;
    return data;
}

function isElementSymbolLine(tokens: string[]): boolean {
    // A dummy species (`X`, `Bq`, `Du`) is a legal name for a site with no
    // element; requiring a real symbol would read the counts line below it as
    // the symbol line and lose the ghost's count.
    return tokens.length > 0 && tokens.every(t =>
        normalizeElementSymbol(t) !== undefined || /^(x|bq|du|dummy|\*)$/i.test(t));
}

export function parseVasp(content: string): MolecularData {
    const lines = content.split(/\r?\n/);

    if (lines.length < 8) {
        return {
            atoms: [], bonds: [], title: 'Invalid VASP file', hasExplicitBonds: false,
            warnings: ['too few lines for a POSCAR: no structure to show']
        };
    }

    const title = lines[0].trim() || 'VASP Structure';
    const warnings: string[] = [];

    // VASP accepts one scale for the whole cell or three per-axis factors.
    // Reading the first number of a three-factor line scaled all three axes by
    // it, so a cell stretched differently along x and z came out cubic.
    const scaleTokens = lines[1].trim().split(/\s+/).map(Number);
    const perAxis = scaleTokens.length >= 3 &&
        scaleTokens.slice(0, 3).every(v => isFinite(v) && v > 0);
    let scaleX: number, scaleY: number, scaleZ: number;
    let scale = 1.0;
    if (perAxis) {
        scaleX = scaleTokens[0];
        scaleY = scaleTokens[1];
        scaleZ = scaleTokens[2];
    } else {
        const rawScale = parseFloat(lines[1].trim());
        scale = isFinite(rawScale) && rawScale !== 0 ? rawScale : 1.0;
        scaleX = scaleY = scaleZ = scale;
    }

    const rawLattice: number[][] = [];
    for (let i = 0; i < 3; i++) {
        const parts = lines[2 + i].trim().split(/\s+/).map(Number);
        if (parts.length < 3 || parts.some(isNaN)) {
            return {
                atoms: [], bonds: [], title: 'Invalid lattice vectors', hasExplicitBonds: false,
                warnings: ['unreadable lattice vectors: no structure to show']
            };
        }
        rawLattice.push([parts[0], parts[1], parts[2]]);
    }

    // A negative scale is VASP's other spelling of a scale: the target cell
    // volume in cubic Angstrom, reached by the factor (V/V0)^(1/3). Applying
    // the number itself would turn a requested 100 A^3 cell into a lattice
    // mirrored through the origin and scaled by 100. The factor applies to
    // Cartesian coordinates too, exactly as a positive scale does.
    if (!perAxis && scale < 0) {
        const target = -scale;
        const v0 = Math.abs(
            rawLattice[0][0] * (rawLattice[1][1] * rawLattice[2][2] - rawLattice[1][2] * rawLattice[2][1]) -
            rawLattice[0][1] * (rawLattice[1][0] * rawLattice[2][2] - rawLattice[1][2] * rawLattice[2][0]) +
            rawLattice[0][2] * (rawLattice[1][0] * rawLattice[2][1] - rawLattice[1][1] * rawLattice[2][0])
        );
        if (v0 > 1e-8) {
            scale = Math.cbrt(target / v0);
        } else {
            warnings.push('negative scale factor on a zero-volume cell: read as 1.0');
            scale = 1;
        }
        scaleX = scaleY = scaleZ = scale;
    }
    const lattice = rawLattice.map((v, i) => v.map(c => c * (i === 0 ? scaleX : i === 1 ? scaleY : scaleZ)));
    if (!isUsableLattice(lattice)) {
        // The scale block above may already have reported a zero-volume target;
        // replacing that list with the degenerate-cell line would lose the
        // reason the file could not be read.
        return {
            atoms: [], bonds: [], title: 'Invalid lattice vectors', hasExplicitBonds: false,
            warnings: warnings.concat(['the POSCAR lattice is degenerate: no structure to show'])
        };
    }

    let lineIdx = 5;

    let elementTypes: string[] = [];
    let elementCounts: number[] = [];

    const line6 = lines[lineIdx] !== undefined ? lines[lineIdx].trim() : '';
    const line6Tokens = line6.split(/\s+/).filter(t => t !== '');

    // A symbol line may spell a species `Si`, `Si_sv` or `H.75`; the symbol
    // the table names is what the viewer can draw, so normalize each token.
    const isElementTypeLine = isElementSymbolLine(line6Tokens) &&
        !(line6Tokens.length === 1 && /^-?\d/.test(line6));

    let symbolCountsMissing = 0;
    // The site total the file states, when it states one. Zero means the
    // counts line was unreadable, and the fallback below reads every remaining
    // row; there is then no declared total to compare the parsed sites with.
    let declaredTotal = 0;
    if (isElementTypeLine) {
        elementTypes = line6Tokens.map(t => normalizeElementSymbol(t) || 'X');
        lineIdx++;

        const countLine = lines[lineIdx] !== undefined ? lines[lineIdx].trim() : '';
        const countTokens = countLine.split(/\s+/).filter(t => t !== '').map(Number);
        if (countTokens.length === elementTypes.length && countTokens.every(n => !isNaN(n) && n >= 0)) {
            elementCounts = countTokens;
            declaredTotal = countTokens.reduce((s, n) => s + n, 0);
            lineIdx++;
        } else {
            warnings.push('atom counts line missing or unreadable: every site is read as the dummy element X');
            elementTypes = [];
        }
    } else {
        // VASP 4 wrote a counts line without symbol names. The counts still
        // say how many sites follow, so the layout is recoverable; only the
        // elements are unknown, and saying so is the point.
        const countTokens = line6Tokens.map(Number);
        if (countTokens.length > 0 &&
            countTokens.every(n => !isNaN(n) && n >= 0 && n % 1 === 0) &&
            countTokens.some(n => n > 0)) {
            symbolCountsMissing = countTokens.reduce((s, n) => s + n, 0);
            declaredTotal = symbolCountsMissing;
            lineIdx++;
        }
        warnings.push('no element symbol line in the file: every site is read as the dummy element X');
    }

    let coordLine = lines[lineIdx] !== undefined ? lines[lineIdx].trim() : '';
    let selective = false;
    if (/^s/i.test(coordLine)) {
        selective = true;
        lineIdx++;
        coordLine = lines[lineIdx] !== undefined ? lines[lineIdx].trim() : '';
    }

    const isFractional = /^d/i.test(coordLine);
    const isCartesian = /^c/i.test(coordLine);
    lineIdx++;

    if (elementTypes.length === 0 || elementCounts.length === 0) {
        const totalAtoms = symbolCountsMissing > 0
            ? symbolCountsMissing
            : Math.max(0, lines.length - lineIdx);
        elementTypes = ['X'];
        elementCounts = [totalAtoms];
    }

    const atoms: Atom[] = [];
    let atomIndex = 0;
    let elementTypeIdx = 0;
    let elementRemaining = elementCounts[0] || 0;

    for (let i = lineIdx; i < lines.length && elementTypeIdx < elementTypes.length; i++) {
        const line = lines[i].trim();
        if (line === '') continue;

        const parts = line.split(/\s+/);
        if (parts.length < 3) continue;
        const x = parseFloat(parts[0]);
        const y = parseFloat(parts[1]);
        const z = parseFloat(parts[2]);

        if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

        // A type with zero sites is legal in a POSCAR, and the row that is
        // being read belongs to the next type that has one. Advancing only on
        // the next iteration would drop that row (a `1 0 1` cell lost the
        // second element's atom), so the zero-count types are stepped over
        // before the coordinates are assigned.
        while (elementRemaining <= 0) {
            elementTypeIdx++;
            if (elementTypeIdx >= elementTypes.length) break;
            elementRemaining = elementCounts[elementTypeIdx] || 0;
        }
        if (elementTypeIdx >= elementTypes.length) break;

        let cx: number, cy: number, cz: number;
        if (isCartesian) {
            cx = x * scaleX;
            cy = y * scaleY;
            cz = z * scaleZ;
        } else {
            const wrapped = [x, y, z];
            const [fx, fy, fz] = fracToCart(wrapped, lattice);
            cx = fx;
            cy = fy;
            cz = fz;
        }

        // Selective dynamics keeps its three flags on the atom so a save can
        // write the column back; a file without the line leaves them unset.
        let constraint: string | undefined;
        if (selective && parts.length >= 6) {
            const f1 = parts[3].charAt(0).toUpperCase();
            const f2 = parts[4].charAt(0).toUpperCase();
            const f3 = parts[5].charAt(0).toUpperCase();
            if ((f1 === 'T' || f1 === 'F') && (f2 === 'T' || f2 === 'F') && (f3 === 'T' || f3 === 'F')) {
                constraint = f1 + ' ' + f2 + ' ' + f3;
            }
        }

        // Occupancy follows the selective flags (VASP 5) or the coordinates
        // (VASP 4); reading it as a flag would drop a partially occupied site
        // to full occupancy without a trace.
        let occupancy = 1;
        const occIndex = selective ? 6 : 3;
        if (parts.length > occIndex) {
            const occ = parseFloat(parts[occIndex]);
            if (isFinite(occ) && occ > 0 && occ <= 1.0001) occupancy = occ;
        }

        atoms.push({
            element: elementTypes[elementTypeIdx],
            x: cx, y: cy, z: cz,
            index: atomIndex,
            occupancy,
            baseIdx: atomIndex,
            constraint
        });
        atomIndex++;
        elementRemaining--;
    }

    // A truncated coordinate table is a fact about the file, not a reason to
    // show the part that happened to fit without comment: the counts line is
    // the file's own statement of how many sites it holds.
    if (declaredTotal > atoms.length) {
        warnings.push('the POSCAR declares ' + declaredTotal + ' sites but only ' + atoms.length +
            ' coordinate row' + (atoms.length === 1 ? ' was' : 's were') + ' readable');
    }

    const cellParams = latticeToCellParams(lattice);

    const crystal: CrystalData = {
        ...cellParams,
        latticeVectors: lattice,
        symmetryOps: ['x, y, z'],
        baseAtoms: atoms.map(a => ({ ...a, cellI: 0, cellJ: 0, cellK: 0 })),
        baseBonds: []
    };

    const data: MolecularData = {
        atoms,
        bonds: [],
        title,
        hasExplicitBonds: false,
        crystal
    };
    if (warnings.length > 0) data.warnings = warnings;
    return data;
}

/**
 * XDATCAR: the header is a POSCAR, and each `Direct configuration=` block is
 * one frame. The viewer's trajectory machinery only moves atoms, not the cell,
 * so the first frame is shown and the frame count is reported instead.
 */
export function parseXdatcar(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];
    let configIndex = -1;
    let mode = 'Direct';
    for (let i = 0; i < lines.length; i++) {
        const m = /^\s*(Direct|Cartesian)\s+configuration\s*=/i.exec(lines[i]);
        if (!m) continue;
        configIndex = i;
        mode = /^cartesian/i.test(m[1]) ? 'Cartesian' : 'Direct';
        break;
    }
    if (configIndex < 0) {
        return {
            atoms: [], bonds: [], title: 'VASP XDATCAR', hasExplicitBonds: false,
            warnings: ['no configuration block in the XDATCAR: no structure to show']
        };
    }

    const rows: string[] = [];
    for (let i = configIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line === '') continue;
        if (/^\s*(Direct|Cartesian)\s+configuration\s*=/i.test(line)) break;
        const parts = line.split(/\s+/);
        if (parts.length < 3 || !isFinite(parseFloat(parts[0])) || !isFinite(parseFloat(parts[1])) || !isFinite(parseFloat(parts[2]))) break;
        rows.push(line);
    }
    if (rows.length === 0) {
        return {
            atoms: [], bonds: [], title: 'VASP XDATCAR', hasExplicitBonds: false,
            warnings: ['the first XDATCAR configuration has no coordinates: no structure to show']
        };
    }

    const header = lines.slice(0, configIndex);
    const selective = rows.length > 0 && rows.every(r => {
        const p = r.split(/\s+/);
        return p.length >= 6 && /^[TF]$/i.test(p[3]) && /^[TF]$/i.test(p[4]) && /^[TF]$/i.test(p[5]);
    });
    // The XDATCAR header is a POSCAR header and may already declare Selective
    // dynamics. A second declaration would be consumed as the coordinate mode,
    // so the line is only added when the header does not state one.
    const headerSelective = header.some(l => /^\s*Selective\s+dynamics\s*$/i.test(l));
    const virtual = header.concat(selective && !headerSelective ? ['Selective dynamics'] : [], [mode], rows).join('\n');
    const data = parseVasp(virtual);
    const allWarnings = warnings.concat(data.warnings || []);
    const frames = (content.match(/(?:Direct|Cartesian)\s+configuration\s*=/gi) || []).length;
    if (frames > 1) allWarnings.push(frames + ' configurations in the XDATCAR: showing the first');
    if (allWarnings.length > 0) data.warnings = allWarnings;
    if (data.title === 'VASP Structure') data.title = 'VASP XDATCAR';
    return data;
}

function outcarElements(lines: string[]): string[] {
    for (const line of lines) {
        const m = /^\s*elements:\s+(.+)$/i.exec(line);
        if (!m) continue;
        const symbols = m[1].trim().split(/\s+/).map(t => normalizeElementSymbol(t) || t);
        if (symbols.length > 0) return symbols;
    }
    const fromTitel: string[] = [];
    for (const line of lines) {
        // `TITEL = PAW_PBE Sr_sv 04Feb2008`: the functional comes first, the
        // element second, the date third. Capturing the third token would name
        // the element after its dataset date.
        const m = /^\s*TITEL\s*=\s*\S+\s+(\S+)/.exec(line);
        if (m) {
            const symbol = normalizeElementSymbol(m[1]) || m[1];
            if (fromTitel[fromTitel.length - 1] !== symbol) fromTitel.push(symbol);
        }
    }
    if (fromTitel.length > 0) return fromTitel;
    const fromVrhfin: string[] = [];
    for (const line of lines) {
        const m = /VRHFIN\s*=\s*([A-Za-z]+)\s*:/.exec(line);
        if (m) fromVrhfin.push(normalizeElementSymbol(m[1]) || m[1]);
    }
    return fromVrhfin;
}

function outcarIonCounts(lines: string[]): number[] {
    let counts: number[] = [];
    for (const line of lines) {
        const m = /ions per type\s*=\s*([\d\s]+)/i.exec(line);
        if (!m) continue;
        counts = m[1].trim().split(/\s+/).map(Number).filter(v => isFinite(v));
    }
    return counts;
}

/**
 * OUTCAR: the cell is the last `direct lattice vectors` group and the geometry
 * is the last POSITION/TOTAL-FORCE table. Both are Cartesian Angstrom. The
 * table itself names no elements; they come from `elements:` or the POTCAR
 * title lines, matched to the `ions per type` counts.
 */
export function parseOutcar(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];

    let lattice: number[][] | null = null;
    for (let i = 0; i < lines.length; i++) {
        if (!/direct lattice vectors/i.test(lines[i])) continue;
        const rows: number[][] = [];
        for (let j = i + 1; j < lines.length && rows.length < 3; j++) {
            const line = lines[j];
            if (line === undefined) break;
            const parts = line.trim().split(/\s+/);
            if (parts.length < 3) break;
            const x = parseFloat(parts[0]);
            const y = parseFloat(parts[1]);
            const z = parseFloat(parts[2]);
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
            rows.push([x, y, z]);
        }
        if (rows.length === 3) lattice = rows;
    }

    let positions: number[][] = [];
    for (let i = 0; i < lines.length; i++) {
        if (!/^\s*POSITION\s+TOTAL-FORCE/i.test(lines[i])) continue;
        const rows: number[][] = [];
        for (let j = i + 1; j < lines.length; j++) {
            const line = lines[j];
            if (line === undefined) break;
            const text = line.trim();
            if (text === '') { if (rows.length > 0) break; continue; }
            if (/^-{3,}/.test(text)) { if (rows.length > 0) break; continue; }
            const parts = text.split(/\s+/);
            if (parts.length < 3) break;
            const x = parseFloat(parts[0]);
            const y = parseFloat(parts[1]);
            const z = parseFloat(parts[2]);
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) break;
            rows.push([x, y, z]);
        }
        if (rows.length > 0) positions = rows;
    }

    if (positions.length === 0) {
        return {
            atoms: [], bonds: [], title: 'VASP OUTCAR', hasExplicitBonds: false,
            warnings: ['the OUTCAR states no POSITION table: no structure to show']
        };
    }

    const species = outcarElements(lines);
    const counts = outcarIonCounts(lines);
    const types: string[] = [];
    for (let t = 0; t < species.length; t++) {
        const n = counts[t] !== undefined ? counts[t] : 0;
        for (let k = 0; k < n; k++) types.push(species[t]);
    }
    if (types.length !== positions.length) {
        warnings.push('the OUTCAR species/count lines do not match the POSITION table: sites are read as the dummy element X');
    }
    const titleMatch = /^\s*POSCAR:\s*(.+)$/m.exec(content);
    const title = titleMatch ? titleMatch[1].trim() : 'VASP OUTCAR';

    const atoms: Atom[] = positions.map((p, i) => ({
        element: types[i] || 'X',
        x: p[0], y: p[1], z: p[2],
        index: i,
        occupancy: 1,
        baseIdx: i
    }));
    if (lattice && isUsableLattice(lattice)) {
        const data: MolecularData = {
            atoms, bonds: [], title, hasExplicitBonds: false,
            crystal: crystalDataFor(lattice, atoms)
        };
        if (warnings.length > 0) data.warnings = warnings;
        return data;
    }
    warnings.push(lattice
        ? 'the OUTCAR direct lattice vectors are degenerate: the geometry is read as a molecule'
        : 'the OUTCAR states no direct lattice vectors: the geometry is read as a molecule');
    return molecularDataFor(atoms, title, warnings);
}

function xmlVarrays(block: string, name: string): number[][] {
    // vasprun.xml is written by hand-formatted code that puts a space before
    // the closing angle bracket (`<varray name="basis" >`); requiring the
    // quote to be followed by `>` directly lost every basis and position.
    const varray = new RegExp('<varray\\s+name\\s*=\\s*"' + name + '"\\s*>([\\s\\S]*?)</varray\\s*>', 'i').exec(block);
    if (!varray) return [];
    const rows: number[][] = [];
    const re = /<v>([\s\S]*?)<\/v>/gi;
    let m: RegExpExecArray | null;
    while ((m = re.exec(varray[1])) !== null) {
        const values = m[1].trim().split(/\s+/).map(Number).filter(v => isFinite(v));
        if (values.length >= 3) rows.push(values.slice(0, 3));
    }
    return rows;
}

/** vasprun.xml: the last `<structure>` carries the final basis and positions. */
export function parseVasprunXml(content: string): MolecularData {
    const warnings: string[] = [];
    if (!/<modeling/i.test(content)) {
        return {
            atoms: [], bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
            warnings: ['the file is not a vasprun.xml: no structure to show']
        };
    }

    const elements: string[] = [];
    const atomArray = /<array\s+name\s*=\s*"atoms"\s*>([\s\S]*?)<\/array\s*>/i.exec(content);
    if (atomArray) {
        const re = /<rc>([\s\S]*?)<\/rc>/gi;
        let m: RegExpExecArray | null;
        while ((m = re.exec(atomArray[1])) !== null) {
            const first = /<c>\s*([A-Za-z]{1,3})\s*<\/c>/.exec(m[1]);
            if (first) elements.push(normalizeElementSymbol(first[1]) || first[1]);
        }
    }

    let structure: string | null = null;
    const structureRe = /<structure[^>]*>([\s\S]*?)<\/structure>/gi;
    let match: RegExpExecArray | null;
    while ((match = structureRe.exec(content)) !== null) structure = match[1];
    if (structure === null) {
        return {
            atoms: [], bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
            warnings: ['the vasprun.xml has no <structure> block: no structure to show']
        };
    }

    const basis = xmlVarrays(structure, 'basis');
    const frac = xmlVarrays(structure, 'positions');
    if (basis.length < 3) {
        return {
            atoms: [], bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
            warnings: ['the vasprun.xml structure states no basis: no structure to show']
        };
    }
    if (frac.length === 0) {
        return {
            atoms: [], bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
            warnings: ['the vasprun.xml structure states no positions: no structure to show']
        };
    }
    const lattice = basis.slice(0, 3);
    if (!isUsableLattice(lattice)) {
        return {
            atoms: [], bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
            warnings: ['the vasprun.xml basis is degenerate: no structure to show']
        };
    }
    if (elements.length !== frac.length) {
        warnings.push('the vasprun.xml atom list does not match its positions: unmatched sites are read as the dummy element X');
    }

    const atoms: Atom[] = frac.map((f, i) => {
        const [x, y, z] = fracToCart(f, lattice);
        return {
            element: elements[i] || 'X',
            x, y, z,
            index: i,
            occupancy: 1,
            baseIdx: i
        };
    });
    const data: MolecularData = {
        atoms, bonds: [], title: 'VASP vasprun.xml', hasExplicitBonds: false,
        crystal: crystalDataFor(lattice, atoms)
    };
    if (warnings.length > 0) data.warnings = warnings;
    return data;
}

export function isXdatcar(content: string): boolean {
    return /(?:Direct|Cartesian)\s+configuration\s*=/i.test(content) &&
        content.split(/\r?\n/).length > 7;
}

export function isOutcar(content: string): boolean {
    return /^\s*vasp\.\d/im.test(content) &&
        /direct lattice vectors/i.test(content) &&
        /POSITION\s+TOTAL-FORCE/i.test(content);
}

export function isVasprunXml(content: string): boolean {
    return /<modeling[\s>]/i.test(content) &&
        /<generator>[\s\S]*?<i\s+name="program"/i.test(content) &&
        /<atominfo>/i.test(content);
}
