import { Atom, AsymAtom, CrystalData, MolecularData } from '../types';
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

const SECTION_KEYWORDS = /^(GROUP|SYMOP|TRANM|LTRANSL|LORIENT|LMATRIX|CELLP|STRUC|THERI|ATOMM|BONDM|POLYM|SURFM|FORMM|SECCL|TEXCL|LIGHT|MAPM|VECTM|SBOND|SANG|BDSPL|BDTYP|CMBOND|ORBIT|NETCHG|AATOMM|SATOMM)\b/;

/** Every site merged away by an inexact match is worth one line, but no more. */
const MAX_MERGE_WARNINGS = 6;

export function parseVesta(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const warnings: string[] = [];

    let title = 'VESTA Structure';
    let spaceGroup: string | undefined;
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
                if (parts.length >= 3) {
                    spaceGroup = parts.slice(2).join(' ');
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
            while (i < lines.length) {
                const strucLine = lines[i].trim();
                if (strucLine === '' || SECTION_KEYWORDS.test(strucLine)) {
                    break;
                }
                const parts = strucLine.split(/\s+/);
                if (parts.length >= 7) {
                    const element = parts[1];
                    const label = parts[2];
                    const occupancy = parseFloat(parts[3]);
                    const fx = parseFloat(parts[4]);
                    const fy = parseFloat(parts[5]);
                    const fz = parseFloat(parts[6]);

                    if (/^[A-Za-z][A-Za-z]?$/.test(element) && isFinite(fx) && isFinite(fy) && isFinite(fz)) {
                        vestaAtoms.push({
                            element,
                            label,
                            occupancy: isFinite(occupancy) ? occupancy : 1,
                            frac: [fx, fy, fz]
                        });
                        i++;
                        // VESTA writes a second, detail line per atom (its
                        // displacement parameters). Skip it, but only when the
                        // next line really is a detail line — a file with the
                        // detail line omitted must not lose its next atom.
                        if (i < lines.length) {
                            const nextLine = lines[i].trim();
                            if (nextLine === '' || SECTION_KEYWORDS.test(nextLine)) continue;
                            const nextParts = nextLine.split(/\s+/);
                            const nextElem = nextParts.length >= 2 ? nextParts[1] : '';
                            const nextFx = nextParts.length >= 5 ? parseFloat(nextParts[4]) : NaN;
                            if (!(nextParts.length >= 7 && /^[A-Za-z][A-Za-z]?$/.test(nextElem) && isFinite(nextFx))) {
                                i++;
                            }
                        }
                        continue;
                    }
                }
                break;
            }
            continue;
        }

        if (line === 'SYMOP') {
            symopSectionSeen = true;
            i++;
            while (i < lines.length) {
                const parts = lines[i].trim().split(/\s+/).map(Number);
                if (parts.length >= 3 && parts[0] === -1 && parts[1] === -1 && parts[2] === -1) break;
                // t (3), R (9) and VESTA's trailing operation-kind column.
                if (parts.length >= 13 && parts.every(v => isFinite(v) && !isNaN(v))) {
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
            i++;
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
    const nearMisses: { label: string; dist: number }[] = [];
    const index = createNearIndex();

    for (let siteIdx = 0; siteIdx < vestaAtoms.length; siteIdx++) {
        const site = vestaAtoms[siteIdx];
        const hitOps: number[] = [];
        for (let opIdx = 0; opIdx < parsedSymOps.length; opIdx++) {
            const transformed = applySymOp(parsedSymOps[opIdx], site.frac);
            const wrapped = transformed.map(wrapFrac);
            const [cx, cy, cz] = fracToCart(wrapped, lattice);

            const near = index.find(cx, cy, cz);
            if (near) {
                // Landing exactly on an existing atom is what symmetry does to
                // a site on a special position. Landing *near* one means two
                // sites the file lists separately are being folded together,
                // which is the case the tolerance exists for and the only one
                // the user needs to hear about.
                if (near.dist > EXACT_MERGE) {
                    nearMisses.push({ label: site.label || site.element + (siteIdx + 1), dist: near.dist });
                }
                continue;
            }

            index.add(cx, cy, cz);
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
        warnings.push('site ' + m.label + ': an image ' + m.dist.toFixed(3) +
            ' A from an existing atom was merged into it');
    }
    if (nearMisses.length > MAX_MERGE_WARNINGS) {
        warnings.push('and ' + (nearMisses.length - MAX_MERGE_WARNINGS) +
            ' more images merged into an atom less than 0.1 A away');
    }

    const crystal: CrystalData = {
        ...cellParams,
        latticeVectors: lattice,
        spaceGroup,
        symmetryOps,
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
