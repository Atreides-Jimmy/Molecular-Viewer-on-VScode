import { Atom, OptStep, NormalMode } from '../types';
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
    const lines = content.split(/\r?\n/);
    if (lines.length < 2) return false;
    const count = parseInt(lines[0].trim(), 10);
    if (isNaN(count) || count <= 0 || count > 1000000) return false;
    return /\bxtb:\s*\S+/i.test(lines[1]);
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
