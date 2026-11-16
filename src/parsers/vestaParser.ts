import { Atom, AsymAtom, CrystalData, MolecularData } from '../types';
import { normalizeElementSymbol } from '../tools/periodicTable';
import { lookupSpaceGroupOps } from '../tools/spaceGroupTable';
import {
    applySymOp,
    cellToLatticeVectors,
    createNearIndex,
    EXACT_MERGE,
    fracToCart,
    SymOp,
    symOpKey,
    symOpToString,
    wrapFrac
} from '../tools/crystalGeometry';

/**
 * VESTA reader.
 *
 * A VESTA file states the same three things a CIF does — a cell, a list of
 * symmetry operations written as numbers, and the asymmetric-unit sites — in a
 * fixed-order, column-oriented layout. `CELLP` carries a/b/c/alpha/beta/gamma,
 * `SYMOP` carries each operation as a translation and a rotation matrix, and
 * `STRUC` lists each site with its label, occupancy and fractional coordinates.
 *
 * The expansion itself is not re-invented here: the cell arithmetic, the
 * coincidence tolerance and its grid live in `crystalGeometry` and are the same
 * code the CIF reader runs. That is deliberate. When the two readers kept
 * private copies they disagreed about how close two images must be to count as
 * one atom, and the same structure read through the two extensions came back
 * with different atom counts — the CIF reader folding only images within 0.1 A
 * while this one folded everything within 0.5 A, silently deleting split sites.
 *
 * Two things the file states were previously thrown away and are now kept:
 *
 *  - The operations themselves. `SYMOP` was applied but `crystal.symmetryOps`
 *    was hard-coded to `x, y, z`, so a save-back or a diff could not tell how
 *    the cell was generated. Each numeric operation is serialized to the CIF
 *    spelling (`-x+1/2, y, -z+1/2`), which is what the rest of the pipeline
 *    already speaks.
 *  - The asymmetric unit. `STRUC` *is* the asymmetric unit, so `asymAtoms`
 *    records each listed site with the operations that actually placed an
 *    image of it — the same traceability the CIF reader provides.
 *
 * Policy follows the CIF reader: a file without a cell comes back empty with a
 * warning instead of being drawn in an invented 1 A cube, and every non-trivial
 * decision (a merge above round-off, an unreadable operation list) is reported
 * through `warnings`.
 */

interface VestaAtom {
    element: string;
    label: string;
    occupancy: number;
    frac: [number, number, number];
}

// Every section header VESTA writes. An unrecognised header does not end the
// STRUC table, so the section's own rows are read as sites (SITET and SHAPE
// rows carry a label and enough numbers to look like one). The list includes
// the sections of the current format that carry no structure data as well, so
// this reader stops at all of them rather than just the ones it consumes.
const SECTION_KEYWORDS = /^(GROUP|SYMOP|TRANM|LTRANSL|LORIENT|LMATRIX|CELLP|CELLT|STRUC|THERI|THERT|ATOMM|ATOMP|ATOMT|BONDM|BONDP|POLYM|POLYP|SURFM|FORMM|FORMP|SECCL|SECTP|TEXCL|TEX3P|LIGHT\d*|MAPM|VECTM|VECTR|VECTT|SBOND|SANG|BDSPL|BDTYP|CMBOND|ORBIT|NETCHG|AATOMM|SATOMM|BOUND|SITET|SHAPE|HKLPM|HKLPP|ISURF|LBLAT|LBLSP|DLATM|DLBND|DLPLY|PLN2D|SPLAN|STYLE|CONTR|SCENE|BKGRC|UCOLP|GRID|OBJECT)\b/;

/** Every site merged away by an inexact match is worth one line, but no more. */
const MAX_MERGE_WARNINGS = 6;

/**
 * The element a `STRUC` row's name column stands for.
 *
 * The column is meant to hold a bare symbol, but files also carry an oxidation
 * state (`Fe3+`, `O2-`), an isotope spelling (`D`, `T`) or one of VESTA's
 * dummy names (`X`, `Bq`, `Du`). Left alone, `Fe3+` is no element the tables
 * know and the site would be dropped; the dummy names are positions, not
 * atoms, and the viewer already draws the dummy element `X` as a ghost.
 */
function vestaElement(raw: string): string | undefined {
    const text = String(raw || '').trim();
    if (text === '') return undefined;
    const symbol = normalizeElementSymbol(text);
    if (symbol) return symbol;
    if (/^[DT]$/i.test(text)) return text.toUpperCase();
    if (/^(x|bq|du|dummy|\*)$/i.test(text)) return 'X';
    return undefined;
}

/**
 * Whether a `GROUP` row's third field is a Hermann-Mauguin symbol.
 *
 * The field is meant to hold one (`P 1`, `P 21/c`, `P a -3`, `R -3 c:r`), so a
 * name has to start with a lattice letter and consist of nothing but the
 * letters, digits and separators the symbols are written with. Anything else is
 * a damaged or truncated row, and keeping it would put a string in the cell
 * panel that names no group at all.
 */
function isHermannMauguin(name: string): boolean {
    const s = name.trim();
    if (s === '') return false;
    if (!/^[PABCIFR]\b/i.test(s) && !/^[PABCIFR][a-z0-9\-\/]/i.test(s)) return false;
    return /^[A-Za-z0-9\s\-\/:._]+$/.test(s);
}

export function parseVesta(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];

    let title = 'VESTA Structure';
    let spaceGroup: string | undefined;
    // The number the `GROUP` row states beside the name, and a name the row
    // carried that is not one — both are checked against the table below.
    let groupNumber = 0;
    let badGroupName = '';
    let cellParams: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } | null = null;
    const vestaAtoms: VestaAtom[] = [];
    const vestaSymOps: SymOp[] = [];
    let symopSectionSeen = false;

    let i = 0;
    while (i < lines.length) {
        const line = lines[i].trim();

        if (line === 'TITLE') {
            i++;
            const titleParts: string[] = [];
            while (i < lines.length && lines[i].trim() !== '' && !SECTION_KEYWORDS.test(lines[i].trim())) {
                titleParts.push(lines[i].trim());
                i++;
            }
            if (titleParts.length > 0) title = titleParts.join(' ');
            continue;
        }

        if (line === 'GROUP') {
            i++;
            if (i < lines.length) {
                const parts = lines[i].trim().split(/\s+/);
                // The row is `number  Z  symbol`: the space-group number, the
                // number of formula units, then the Hermann-Mauguin symbol,
                // which may carry spaces of its own (`P a -3`, `P 21/c`).
                const number = parseInt(parts[0], 10);
                const symbol = parts.slice(2).join(' ');
                groupNumber = isFinite(number) && number >= 1 && number <= 230 ? number : 0;
                // A name is only kept when it is one: a truncated or damaged
                // row used to become a group string of whatever was left, and
                // the panel then named a space group the file never stated.
                if (isHermannMauguin(symbol)) {
                    spaceGroup = symbol;
                } else if (parts.length >= 3) {
                    badGroupName = symbol;
                }
            }
            i++;
            continue;
        }

        if (line === 'CELLP') {
            i++;
            if (i < lines.length) {
                // The first line holds the six parameters; the second holds
                // their standard uncertainties and is not part of the cell.
                const parts = lines[i].trim().split(/\s+/).map(Number);
                if (parts.length >= 6 && parts.slice(0, 6).every(v => isFinite(v) && !isNaN(v))) {
                    cellParams = {
                        a: parts[0], b: parts[1], c: parts[2],
                        alpha: parts[3], beta: parts[4], gamma: parts[5]
                    };
                }
            }
            i++;
            continue;
        }

        if (line === 'STRUC') {
            i++;
            let skippedRows = 0;
            while (i < lines.length) {
                const strucLine = lines[i].trim();
                // A blank line is whitespace, not a terminator: VESTA writes
                // one between the table and the next section, and a file that
                // leaves one between two sites must not lose everything below
                // it. The section ends at the next keyword, not at an empty
                // line.
                if (strucLine === '') {
                    i++;
                    continue;
                }
                if (SECTION_KEYWORDS.test(strucLine)) break;
                const parts = strucLine.split(/\s+/);
                // Each site may be followed by a detail line of displacement
                // parameters. It is numeric throughout, and treating it as a
                // damaged site row would put one warning per atom in the panel.
                if (parts.every(p => p !== '' && isFinite(Number(p)))) {
                    i++;
                    continue;
                }
                const element = parts.length >= 7 ? vestaElement(parts[1]) : undefined;
                const fx = parts.length >= 7 ? parseFloat(parts[4]) : NaN;
                const fy = parts.length >= 7 ? parseFloat(parts[5]) : NaN;
                const fz = parts.length >= 7 ? parseFloat(parts[6]) : NaN;
                if (element !== undefined && isFinite(fx) && isFinite(fy) && isFinite(fz)) {
                    const occupancy = parseFloat(parts[3]);
                    vestaAtoms.push({
                        element,
                        label: parts[2],
                        occupancy: isFinite(occupancy) ? occupancy : 1,
                        frac: [fx, fy, fz]
                    });
                } else {
                    // One damaged row is skipped, not treated as the end of the
                    // table: a single unreadable line must not hide every site
                    // below it. The count is reported once, after the section.
                    skippedRows++;
                }
                i++;
            }
            if (skippedRows > 0) {
                warnings.push(skippedRows + ' STRUC row' + (skippedRows === 1 ? '' : 's') +
                    ' skipped: not a site with an element and three coordinates');
            }
            continue;
        }

        if (line === 'SYMOP') {
            symopSectionSeen = true;
            i++;
            while (i < lines.length) {
                const symopLine = lines[i].trim();
                if (symopLine === '') {
                    i++;
                    continue;
                }
                const parts = symopLine.split(/\s+/).map(Number);
                if (parts.length >= 3 && parts[0] === -1 && parts[1] === -1 && parts[2] === -1) {
                    i++;
                    break;
                }
                // A row holds t (3), R (9) and VESTA's operation-kind column;
                // the trailing column is optional in files written by other
                // tools. The first non-numeric line is the next section, and
                // walking into it would read its numeric rows as operations.
                if (parts.some(v => !isFinite(v) || isNaN(v))) break;
                if (parts.length >= 12) {
                    vestaSymOps.push({
                        t: [parts[0], parts[1], parts[2]],
                        R: [
                            [parts[3], parts[4], parts[5]],
                            [parts[6], parts[7], parts[8]],
                            [parts[9], parts[10], parts[11]]
                        ]
                    });
                }
                i++;
            }
            continue;
        }

        i++;
    }

    // ----- cell -----
    // No cell means no periodicity and no Cartesian frame. Inventing one (the
    // old code treated the fractional coordinates as Angstrom) would draw a
    // structure the file never described, so the structure comes back empty and
    // the warning says why — the same policy the CIF reader follows.
    if (!cellParams) {
        warnings.push('no unit cell in the file (CELLP missing): no structure to show');
        return { atoms: [], bonds: [], title, hasExplicitBonds: false, warnings };
    }

    const lattice = cellToLatticeVectors(cellParams.a, cellParams.b, cellParams.c, cellParams.alpha, cellParams.beta, cellParams.gamma);

    // ----- symmetry -----
    // Duplicate operations would expand every site twice; compare the numbers,
    // not the text, so a table with two spellings of the same operation cannot
    // double the cell.
    const parsedSymOps: SymOp[] = [];
    const seenOps: { [key: string]: true } = {};
    for (const op of vestaSymOps) {
        const key = symOpKey(op);
        if (seenOps[key]) continue;
        seenOps[key] = true;
        parsedSymOps.push(op);
    }

    if (parsedSymOps.length === 0) {
        parsedSymOps.push({ R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] });
        const declared = spaceGroup ? ' (declared space group ' + spaceGroup + ')' : '';
        if (symopSectionSeen) {
            warnings.push('no readable symmetry operations in the file: read as P1' + declared);
        } else {
            warnings.push('no symmetry operations in the file: read as P1' + declared);
        }
    }
    const symmetryOps = parsedSymOps.map(symOpToString);

    // ----- expand -----
    const baseAtoms: Atom[] = [];
    const asymAtoms: AsymAtom[] = [];
    const nearMisses: { label: string; dist: number; otherLabel?: string; crossSite?: boolean }[] = [];
    const index = createNearIndex();
    // Which listed site each placed atom came from. It is what separates a
    // special position folding an image onto itself (normal and silent) from
    // two different sites landing on the same point — a row, its occupancy and
    // its label would otherwise vanish without a trace.
    const siteOf: number[] = [];

    for (let siteIdx = 0; siteIdx < vestaAtoms.length; siteIdx++) {
        const site = vestaAtoms[siteIdx];
        const hitOps: number[] = [];
        for (let opIdx = 0; opIdx < parsedSymOps.length; opIdx++) {
            const transformed = applySymOp(parsedSymOps[opIdx], site.frac);
            const wrapped = transformed.map(wrapFrac);
            const [cx, cy, cz] = fracToCart(wrapped, lattice);

            const near = index.find(cx, cy, cz);
            if (near) {
                // Landing on an atom of the *same* site is what symmetry does
                // to a site on a special position. Landing on another site
                // means two rows the file lists separately are being folded
                // into one, and the second row's occupancy and position leave
                // the cell — that is reported whether the overlap is exact or
                // merely close, and the message names both rows.
                const otherSite = siteOf[near.index];
                const crossSite = otherSite !== undefined && otherSite !== siteIdx;
                if (crossSite) {
                    const other = vestaAtoms[otherSite];
                    nearMisses.push({
                        label: site.label || site.element + (siteIdx + 1),
                        dist: near.dist,
                        crossSite: true,
                        otherLabel: other.label || other.element + (otherSite + 1)
                    });
                } else if (near.dist > EXACT_MERGE) {
                    nearMisses.push({ label: site.label || site.element + (siteIdx + 1), dist: near.dist });
                }
                continue;
            }

            index.add(cx, cy, cz);
            siteOf.push(siteIdx);
            baseAtoms.push({
                element: site.element,
                x: cx, y: cy, z: cz,
                index: baseAtoms.length,
                occupancy: site.occupancy,
                baseIdx: baseAtoms.length,
                siteLabel: site.label || undefined,
                asymIdx: asymAtoms.length,
                symOp: opIdx,
                cellI: 0, cellJ: 0, cellK: 0
            });
            hitOps.push(opIdx);
        }

        const siteFrac: [number, number, number] = [wrapFrac(site.frac[0]), wrapFrac(site.frac[1]), wrapFrac(site.frac[2])];
        const [sx, sy, sz] = fracToCart(siteFrac, lattice);
        asymAtoms.push({
            element: site.element,
            label: site.label,
            x: sx, y: sy, z: sz,
            frac: siteFrac,
            index: asymAtoms.length,
            occupancy: site.occupancy,
            symmetryOps: hitOps
        });
    }

    if (vestaAtoms.length === 0) {
        warnings.push('no atom sites in the file');
    }
    for (let k = 0; k < nearMisses.length && k < MAX_MERGE_WARNINGS; k++) {
        const m = nearMisses[k];
        if (m.crossSite) {
            warnings.push('site ' + m.label + ': its image ' + m.dist.toFixed(3) +
                ' A from site ' + m.otherLabel + ' was merged into it; ' +
                'the second site and its occupancy are not in the unit-cell list');
        } else {
            warnings.push('site ' + m.label + ': an image ' + m.dist.toFixed(3) +
                ' A from an existing atom was merged into it');
        }
    }
    if (nearMisses.length > MAX_MERGE_WARNINGS) {
        warnings.push('and ' + (nearMisses.length - MAX_MERGE_WARNINGS) +
            ' more images merged into an atom less than 0.1 A away');
    }

    // ----- the space group the file names -----
    // `GROUP` states the number and the symbol, and the operations above are
    // what the picture is actually built from. A name the table does not know,
    // or one whose number disagrees with its own symbol, describes a group the
    // file contradicts itself about; drawing is unaffected (the file's own
    // operations win) but the panel must not present the name as fact.
    if (badGroupName !== '') {
        warnings.push('the GROUP record does not state a space-group symbol (' + badGroupName +
            '): it is not shown as one');
    }
    if (spaceGroup !== undefined) {
        const hit = lookupSpaceGroupOps(spaceGroup, cellParams);
        if (hit === null) {
            warnings.push('the GROUP record names ' + spaceGroup +
                ', which the built-in space-group table does not know: the cell is drawn from the file\'s own operations');
        } else if (groupNumber > 0 && hit.number !== groupNumber) {
            warnings.push('the GROUP record names ' + spaceGroup + ' (No. ' + hit.number +
                ') beside the number ' + groupNumber + ': the cell is drawn from the file\'s own operations');
        } else if (parsedSymOps.length !== hit.ops.length) {
            warnings.push('the GROUP record names ' + spaceGroup + ', whose ' + hit.ops.length +
                ' operation' + (hit.ops.length === 1 ? '' : 's') + ' expand a cell, but the file lists ' +
                parsedSymOps.length + ': the file\'s operations are used');
        }
    }

    const crystal: CrystalData = {
        ...cellParams,
        latticeVectors: lattice,
        spaceGroup,
        symmetryOps,
        // The operations are the file's own table, so a save-back or a diff can
        // say where the cell came from — the same provenance a CIF carries.
        symmetrySource: parsedSymOps.length > 0 && vestaSymOps.length > 0 ? 'file' : undefined,
        baseAtoms,
        baseBonds: [],
        asymAtoms
    };

    return {
        atoms: baseAtoms.map(a => ({ ...a })),
        bonds: [],
        title,
        hasExplicitBonds: false,
        warnings,
        crystal
    };
}
