import { Atom, Bond, isGhostAtom } from '../types';

/**
 * Ghost atoms: one spelling in memory, many in the files.
 *
 * A ghost is a position without a nucleus, electrons or bonds — the marker a
 * counterpoise calculation puts where the other fragment's atoms are, and the
 * point a geometry is defined against. No two programs agree on how to write
 * one, so the parsers all funnel their atom label through here and the model
 * keeps a single shape: `element` names the atom the ghost stands for (a ghost
 * carbon keeps `C`, because Gaussian's `C-Bq` needs it) plus `ghost: true`, and
 * a bare dummy point is `element: 'X'`.
 *
 * The spellings recognised, and where each comes from:
 *
 *   `C-Bq`   Gaussian — a ghost with the basis functions of a carbon
 *   `Bq`     Gaussian, and the pseudo-atom older tools write for a ghost
 *   `C:`     ORCA — the same ghost, marked with a trailing colon
 *   `DA`     ORCA `/` MOPAC `/` ADF — dummy atom (`X`, `Xx`, `DA`, `XX`)
 *   `Du`     Tripos MOL2 atom type for a dummy
 *   `q`      TURBOMOLE `$coord` dummy (handled by the coord parser only: `Q`
 *            means a point charge in ORCA input, so it is not a global alias)
 *
 * Normalising on the way in is what makes a structure survive a round trip: a
 * file's `Du` becomes `X` here, is written back as `Du` to a MOL2 and as `Bq` to
 * a Gaussian input, and is still the same ghost in the viewer.
 */

/** Labels that mean "a dummy point" in at least one program, upper-cased. */
const DUMMY_LABELS: { [key: string]: number } = {
    X: 1, XX: 1, DA: 1, DU: 1, DUM: 1, DUMMY: 1, BQ: 1
};

/** `C-Bq`, `C_bq`, `C(Bq)`, `C-BQ` … — a ghost of a named element. */
const ELEMENT_GHOST_RE = /^([A-Za-z]{1,2})\s*(?:[-_]\s*bq|\(bq\))$/i;

/** `C:` — ORCA's ghost marker, after the element symbol. */
const COLON_GHOST_RE = /^([A-Za-z]{1,2})\s*:$/;

function cleanLabel(raw: string): string {
    let s = String(raw == null ? '' : raw).trim();
    // A quoted label (CIF-style) may arrive with its quotes.
    if (s.length > 1 && s.charAt(0) === s.charAt(s.length - 1) && (s.charAt(0) === '"' || s.charAt(0) === "'")) {
        s = s.substring(1, s.length - 1).trim();
    }
    return s;
}

/** Capitalise an element symbol the way the model stores it (`cl` → `Cl`). */
export function ghostElementSymbol(raw: string): string {
    const s = cleanLabel(raw);
    if (s.length === 0) return s;
    return s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
}

/**
 * Read one atom label as a ghost, or return null when it is an ordinary atom.
 *
 * The answer is `{ element, ghost: true }`: the element the ghost stands for,
 * or `X` for a bare dummy point. Returning null is the common case and leaves
 * the label to the parser's own element resolution, so nothing about how a
 * normal atom is read changes.
 */
export function ghostFromLabel(raw: string): { element: string; ghost: true } | null {
    const label = cleanLabel(raw);
    if (label.length === 0) return null;

    const elementGhost = ELEMENT_GHOST_RE.exec(label);
    if (elementGhost) return { element: ghostElementSymbol(elementGhost[1]), ghost: true };

    const colonGhost = COLON_GHOST_RE.exec(label);
    if (colonGhost) return { element: ghostElementSymbol(colonGhost[1]), ghost: true };

    const upper = label.toUpperCase();
    if (DUMMY_LABELS[upper] === 1) return { element: 'X', ghost: true };

    return null;
}

/**
 * Drop every bond that touches a ghost.
 *
 * A ghost has no valence: a file's own bond table (PDB `CONECT`, MOL2 `BOND`,
 * a Gaussian `connect` block) can still name one, and carrying such a bond into
 * the viewer would draw a stick into a point that has no atom — and would make
 * the MM force field and the valence solver treat a nucleus-less point as a
 * bonded centre. Called on every load path through `ensureBonds`, so no parser
 * has to remember it.
 */
export function stripGhostBonds(atoms: Atom[], bonds: Bond[]): Bond[] {
    let hasGhost = false;
    for (let i = 0; i < atoms.length; i++) {
        if (isGhostAtom(atoms[i])) { hasGhost = true; break; }
    }
    if (!hasGhost) return bonds;
    return bonds.filter(function (b) {
        const a1 = atoms[b.atom1];
        const a2 = atoms[b.atom2];
        if (!a1 || !a2) return false;
        return !isGhostAtom(a1) && !isGhostAtom(a2);
    });
}
