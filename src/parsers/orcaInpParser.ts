import { Atom, Bond, MolecularData } from '../types';
import { ghostFromLabel } from '../tools/ghostAtoms';
import { parseXyzTrajectory } from './xyzParser';

/** A structure file an ORCA input names with `* xyzfile`. */
export interface OrcaExternalStructure {
    content: string;
    fileName: string;
}

/** One `* xyzfile` reference, in file order. */
export interface OrcaXyzReference {
    name: string;
    charge: number;
    multiplicity: number;
}

/**
 * One atom label of an ORCA coordinate block.
 *
 * ORCA spells a dummy `DA` / `X` / `Xx` and a ghost `<element>:` (a colon right
 * after the symbol), so the ghost reading comes first; everything else keeps
 * the old rule — a one- or two-letter element symbol.
 */
function orcaAtomLabel(raw: string): { element: string; ghost: boolean } | null {
    const ghost = ghostFromLabel(raw);
    if (ghost) return { element: ghost.element, ghost: true };
    if (!/^[A-Za-z]{1,2}$/.test(raw)) return null;
    return { element: raw.charAt(0).toUpperCase() + raw.slice(1).toLowerCase(), ghost: false };
}

/** Remove an ORCA `#` comment, but not a `#` inside a quoted file name. */
function stripOrcaComment(line: string): string {
    let quote = '';
    for (let i = 0; i < line.length; i++) {
        const ch = line.charAt(i);
        if (quote !== '') {
            if (ch === quote) quote = '';
            continue;
        }
        if (ch === '"' || ch === "'") {
            quote = ch;
            continue;
        }
        if (ch === '#') return line.slice(0, i);
    }
    return line;
}

/** Split one ORCA line into words, keeping a quoted file name as one word. */
function orcaWords(line: string): string[] {
    const matches = stripOrcaComment(line).match(/"[^"]*"|'[^']*'|\S+/g);
    if (!matches) return [];
    return matches.map(word => {
        if (word.length >= 2 && ((word.charAt(0) === '"' && word.endsWith('"')) ||
            (word.charAt(0) === "'" && word.endsWith("'")))) {
            return word.slice(1, -1);
        }
        return word;
    });
}

/**
 * The geometry keyword and its arguments on a `*` line. ORCA accepts both
 * `* xyzfile 0 1 coords.xyz` and the unspaced `*xyzfile` spelling.
 */
function orcaStarLine(words: string[]): { keyword: string; args: string[] } | null {
    if (words.length === 0 || words[0].charAt(0) !== '*') return null;
    if (words[0].length > 1) {
        return { keyword: words[0].slice(1).toLowerCase(), args: words.slice(1) };
    }
    if (words.length < 2) return { keyword: '', args: [] };
    return { keyword: words[1].toLowerCase(), args: words.slice(2) };
}

/** Every `* xyzfile <charge> <multiplicity> <name>` reference in an input. */
export function scanOrcaXyzfile(content: string): OrcaXyzReference[] {
    const refs: OrcaXyzReference[] = [];
    for (const raw of content.split(/\r?\n/)) {
        const star = orcaStarLine(orcaWords(raw));
        if (!star || star.keyword !== 'xyzfile') continue;
        const name = star.args.length >= 3 ? star.args[2] : '';
        if (name === '') continue;
        refs.push({
            name,
            charge: parseInt(star.args[0], 10) || 0,
            multiplicity: parseInt(star.args[1], 10) || 1
        });
    }
    return refs;
}

export function parseOrcainp(content: string, external?: OrcaExternalStructure): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge = 0;
    let multiplicity = 1;

    const references = scanOrcaXyzfile(content);

    let inXyz = false;
    let inCoordsBlock = false;
    let inCoordsBohr = false;
    let atomIndex = 0;

    for (let i = 0; i < lines.length; i++) {
        const rawLine = lines[i];
        const line = rawLine.trim();

        if (line.startsWith('#')) continue;

        if (line.startsWith('*')) {
            const star = orcaStarLine(orcaWords(rawLine));
            if (star && (star.keyword === 'xyz' || star.keyword === 'xyzfile')) {
                // `xyzfile` has no inline rows: its block is the path in args.
                inXyz = star.keyword === 'xyz';
                if (star.args.length >= 2) {
                    charge = parseInt(star.args[0], 10) || 0;
                    multiplicity = parseInt(star.args[1], 10) || 1;
                }
                continue;
            }
            inXyz = false;
            continue;
        }

        // A `%coords` block after an unterminated inline block must not be
        // read twice: the inline block ends when a new `%` section starts.
        if (line.startsWith('%')) inXyz = false;

        if (inXyz && line !== '') {
            const parts = line.split(/\s+/);
            if (parts.length >= 4) {
                const label = orcaAtomLabel(parts[0]);
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (label && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                    const atom: Atom = { element: label.element, x, y, z, index: atomIndex };
                    if (label.ghost) atom.ghost = true;
                    atoms.push(atom);
                    atomIndex++;
                }
            }
        }

        if (line.toLowerCase().startsWith('%coords')) {
            inCoordsBlock = true;
            inCoordsBohr = false;
            continue;
        }

        if (inCoordsBlock) {
            if (line.toLowerCase() === 'end') {
                inCoordsBlock = false;
                inCoordsBohr = false;
                continue;
            }
            const lowerLine = line.toLowerCase();
            if (lowerLine.startsWith('xyz')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 3) {
                    charge = parseInt(parts[1], 10) || 0;
                    multiplicity = parseInt(parts[2], 10) || 1;
                }
                continue;
            }
            if (lowerLine.startsWith('units')) {
                if (lowerLine.includes('bohr') || lowerLine.includes('angs')) {
                    inCoordsBohr = lowerLine.includes('bohr');
                }
                continue;
            }
            if (line !== '' && !line.startsWith('%')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 4) {
                    const label = orcaAtomLabel(parts[0]);
                    let x = parseFloat(parts[1]);
                    let y = parseFloat(parts[2]);
                    let z = parseFloat(parts[3]);
                    if (label && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        if (inCoordsBohr) {
                            const BOHR = 0.529177249;
                            x *= BOHR; y *= BOHR; z *= BOHR;
                        }
                        const atom: Atom = { element: label.element, x, y, z, index: atomIndex };
                        if (label.ghost) atom.ghost = true;
                        atoms.push(atom);
                        atomIndex++;
                    }
                }
            }
        }
    }

    // Only the first xyzfile is the geometry this input starts from; a later
    // job in the same file would be a different structure.
    const reference = references.length > 0 ? references[0] : undefined;
    if (reference) {
        if (external) {
            const trajectory = parseXyzTrajectory(external.content);
            const frame = trajectory.frames.length > 0 ? trajectory.frames[0] : undefined;
            if (frame && frame.atoms.length > 0) {
                // The named file is the geometry; a `%coords` block read above
                // was only a fallback and must not add a second copy.
                atoms.length = 0;
                atomIndex = 0;
                for (const atom of frame.atoms) {
                    const copied: Atom = {
                        element: atom.element,
                        x: atom.x,
                        y: atom.y,
                        z: atom.z,
                        index: atomIndex++
                    };
                    if (atom.ghost) copied.ghost = true;
                    atoms.push(copied);
                }
                if (title === '' && frame.title !== '') title = frame.title;
            } else {
                return {
                    atoms: [], bonds: [], title, hasExplicitBonds: false,
                    charge, multiplicity,
                    warnings: ['xyzfile ' + reference.name + ': the coordinate file has no atoms']
                };
            }
        } else {
            return {
                atoms: [], bonds: [], title, hasExplicitBonds: false,
                charge, multiplicity,
                warnings: ['the ORCA input reads its structure from "' + reference.name +
                    '", but that file could not be read']
            };
        }
    }

    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
