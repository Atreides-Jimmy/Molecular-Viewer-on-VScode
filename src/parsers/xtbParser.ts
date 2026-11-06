import { Atom, MolecularData, OptStep, NormalMode } from '../types';
import { normalizeElementSymbol } from '../tools/periodicTable';
import { LogFrame } from './logParser';
import { resolveElement } from './xyzParser';

export interface XtbLogResult {
    frames: LogFrame[];
    title: string;
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
}

/**
 * Sniff whether a .log file is an xtb trajectory (extended multi-frame XYZ
 * with `xtb:` marker in the comment line) rather than a Gaussian log.
 *
 * Heuristic: first line is a positive integer (atom count) AND second line
 * contains the `xtb:` program marker. Gaussian logs never match this pattern
 * (they start with a route card or banner text, not a bare integer).
 */
export function isXtbLog(content: string): boolean {
    // Only the first two lines are needed; slicing avoids allocating the full
    // line array just for the sniff (the parser splits the content again).
    const firstNl = content.indexOf('\n');
    if (firstNl < 0) return false;
    const count = parseInt(content.slice(0, firstNl).trim(), 10);
    if (isNaN(count) || count <= 0 || count > 1000000) return false;
    const rest = content.slice(firstNl + 1);
    const secondNl = rest.indexOf('\n');
    const second = (secondNl < 0 ? rest : rest.slice(0, secondNl)).trim();
    return /\bxtb:\s*\S+/i.test(second);
}

/**
 * Parse an xtb optimization trajectory (extended multi-frame XYZ).
 *
 * Each frame is:
 *   <atomCount>
 *   energy: <Eh> gnorm: <grad_norm> xtb: <version> (<commit>)
 *   <element> <x> <y> <z>
 *   ...
 *
 * - Coordinates are in Angstrom.
 * - Element symbols are lowercase in the file (resolved to proper case).
 * - The last frame may be an appended "final geometry" with energy=0 and
 *   gnorm=0 as placeholders; it is labeled "Final" and excluded from optSteps.
 * - gnorm is mapped to rmsForce so the convergence panel renders a force curve.
 */
export function parseXtbLog(content: string): XtbLogResult {
    const lines = content.split(/\r?\n/);
    const frames: LogFrame[] = [];
    const optSteps: OptStep[] = [];
    let stepCounter = 0;

    let i = 0;
    while (i < lines.length) {
        const atomCount = parseInt(lines[i].trim(), 10);
        if (isNaN(atomCount) || atomCount <= 0) { i++; continue; }

        if (i + 1 >= lines.length) break;
        const comment = lines[i + 1].trim();

        // Parse energy and gnorm from the comment line
        const energyMatch = comment.match(/energy:\s*(-?\d+\.\d+)/i);
        const gnormMatch = comment.match(/gnorm:\s*(-?\d+\.\d+)/i);
        const energy = energyMatch ? parseFloat(energyMatch[1]) : undefined;
        const gnorm = gnormMatch ? parseFloat(gnormMatch[1]) : undefined;

        // Parse atom lines
        const atoms: Atom[] = [];
        let atomIndex = 0;
        const atomStart = i + 2;
        const atomEnd = Math.min(atomStart + atomCount, lines.length);
        for (let j = atomStart; j < atomEnd; j++) {
            const line = lines[j].trim();
            if (line === '') continue;
            const parts = line.split(/\s+/);
            if (parts.length >= 4) {
                const element = resolveElement(parts[0]);
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (element && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                    atoms.push({ element, x, y, z, index: atomIndex });
                    atomIndex++;
                }
            }
        }

        if (atoms.length === 0) { i++; continue; }

        // xtb appends a final-geometry frame with energy=0 AND gnorm=0 as
        // placeholders (the real final energy is in the preceding frame).
        // A molecular total energy is never exactly 0 Hartree, so this is safe.
        const isPlaceholder = energy === 0 && gnorm === 0;

        let stepLabel: string;
        if (isPlaceholder) {
            stepLabel = 'Final';
        } else {
            stepCounter++;
            stepLabel = `Step ${stepCounter}`;
        }

        frames.push({
            atoms,
            bonds: [],
            title: comment,
            hasExplicitBonds: false,
            stepLabel,
            // xtb prints one energy per frame, so every frame maps 1:1 onto an
            // optStep — the webview uses it to highlight the point of the frame
            // being viewed. [2026-09-12 | Atreides-Jimmy]
            energy,
        });

        if (!isPlaceholder) {
            optSteps.push({
                step: stepCounter,
                energy,
                rmsForce: gnorm,
            });
        }

        i = atomEnd;
    }

    return {
        frames,
        title: 'xtb optimization',
        optSteps: optSteps.length > 0 ? optSteps : undefined,
        normalModes: undefined,
    };
}

function xtbOutputAtom(line: string): Atom | undefined {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) return undefined;
    const element = normalizeElementSymbol(parts[0]);
    if (!element) return undefined;
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return undefined;
    return { element, x, y, z, index: 0 };
}

/**
 * Whether the text is xTB's standard output (the run log) rather than the
 * multi-frame trajectory xTB writes separately.
 *
 * The run log ends with a `final structure:` block in XYZ layout; the xTB
 * banner in the header tells it apart from a CREST or a Gaussian output that
 * happens to contain a similarly labelled section.
 */
export function isXtbOutput(content: string): boolean {
    // A trajectory starts with the atom count, never with the banner.
    if (isXtbLog(content)) return false;
    if (!/^\s*\*?\s*final structure:\s*$/im.test(content)) return false;
    const head = content.slice(0, 120000);
    return /\bxtb\b/i.test(head) && /\bversion\b/i.test(head);
}

/**
 * Read the final geometry from an xTB run log.
 *
 * The block layout is
 *
 *   ===============
 *    final structure:
 *   ===============
 *   6
 *   xtb: 6.7.1 (edcfbbe)
 *   O  0.0000  0.0000  0.1173
 *   ...
 *
 * The line after the atom count is xTB's comment, not a coordinate row; some
 * builds add another version line or leave the comment empty, so rows are
 * recognised by their shape instead of by counting lines.
 */
export function parseXtbOutput(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    let atoms: Atom[] = [];
    for (let i = 0; i < lines.length; i++) {
        if (!/^\s*\*?\s*final structure:\s*$/i.test(lines[i])) continue;
        // Find the atom count after the heading and its ruler.
        let countAt = i + 1;
        while (countAt < lines.length && !/^\s*\d+\s*$/.test(lines[countAt])) {
            if (countAt - i > 8) { countAt = lines.length; break; }
            countAt++;
        }
        if (countAt >= lines.length) continue;
        const count = parseInt(lines[countAt].trim(), 10);
        if (!(count > 0) || count > 1000000) { i = countAt; continue; }

        const block: Atom[] = [];
        let skipped = 0;
        let k = countAt + 1;
        while (k < lines.length && block.length < count) {
            const line = lines[k].trim();
            k++;
            if (line === '' || /^[=*-]+$/.test(line)) continue;
            const atom = xtbOutputAtom(line);
            if (!atom) {
                // The comment and version lines live between the count and the
                // first row; the table itself ends at the first non-row.
                if (block.length === 0 && ++skipped <= 5) continue;
                break;
            }
            atom.index = block.length;
            block.push(atom);
        }
        if (block.length === count) atoms = block;
        i = k - 1;
    }

    if (atoms.length === 0) {
        return {
            atoms: [], bonds: [], title: 'xtb output', hasExplicitBonds: false,
            warnings: ['the xtb output has no readable final structure: ' +
                'open the trajectory (xtbopt.log or xtb.trj) for the individual frames']
        };
    }
    return { atoms, bonds: [], title: 'xtb output', hasExplicitBonds: false };
}
