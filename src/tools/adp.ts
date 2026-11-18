import { applySymOp, cartToFrac, CELL_VOLUME_EPSILON, fracToCart, SymOp } from './crystalGeometry';

/**
 * Anisotropic displacement parameters (ADPs), as a CIF states them and as the
 * viewer stores them.
 *
 * The six numbers are:
 *
 *     [U11, U22, U33, U12, U13, U23]
 *
 * A CIF states them in the **crystal-axis basis**: the displacement is
 * `u = x a + y b + z c`, and U_ij is the mean product of the components along
 * the (non-orthogonal) cell axes, in Angstrom squared. They are *not* Cartesian,
 * and the two only coincide when the cell is orthogonal. The dictionary's own
 * conversion is the one cctbx implements: scale the tensor by the reciprocal
 * cell lengths (`u_star_ii = r_i^2 u_cif_ii`, `u_star_ij = r_i r_j u_cif_ij`),
 * then rotate by the lattice matrix (`U_cart = M^T u_star M`, M's rows being
 * the lattice vectors in the viewer's frame).
 *
 * The viewer stores every U as the Cartesian tensor, because that is where the
 * ellipsoid, the symmetry rotation and the renderer all work. `uCifAsCart`
 * converts a file's numbers on read; `uCartAsCif` is the inverse a CIF writer
 * needs, or a non-orthogonal cell would save a tensor that means a different
 * ellipsoid when read back.
 *
 * A site is listed once but a cell holds every symmetry image of it, and an
 * ADP is a second-rank tensor: the image produced by an operation R is
 * `R U R^T`, not `U`. A general-position site under a rotation therefore has
 * one orientation per image, which is what makes an ORTEP drawing of a cubic
 * group readable. `rotateAnisoTensor` is that transform, done in Cartesian
 * space because that is where the stored tensor lives.
 */

/** `[U11, U22, U33, U12, U13, U23]`, Cartesian Angstrom squared. */
export type AnisoTensor = [number, number, number, number, number, number];

/**
 * The Cartesian rotation a fractional-coordinate operation induces.
 *
 * The page's frame maps fractional `f` to Cartesian `c = f · M`, with the
 * lattice vectors as the rows of `M` (see `fracToCart`). An operation acts as
 * `f' = R f` (column form), so `c' = c · M⁻¹ R M`: the matrix returned here,
 * applied to a column Cartesian vector, is `Mᵀ R (M⁻¹)ᵀ`. It is computed by
 * pushing the three unit Cartesian axes through the operation rather than by
 * multiplying matrices, which keeps the oblique-cell case in one place — the
 * same `cartToFrac`/`fracToCart` pair every other periodic calculation uses.
 */
export function symOpCartesianRotation(op: SymOp, lattice: number[][]): number[][] {
    const R: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let col = 0; col < 3; col++) {
        const unit = [col === 0 ? 1 : 0, col === 1 ? 1 : 0, col === 2 ? 1 : 0];
        const frac = cartToFrac(unit, lattice);
        const mapped = applySymOp({ R: op.R, t: [0, 0, 0] }, frac);
        const cart = fracToCart(mapped, lattice);
        for (let row = 0; row < 3; row++) R[row][col] = cart[row];
    }
    return R;
}

/** `U' = R U Rᵀ` for a symmetric tensor stored as its six independent values. */
export function rotateAnisoTensor(u: AnisoTensor, R: number[][]): AnisoTensor {
    const U = [
        [u[0], u[3], u[4]],
        [u[3], u[1], u[5]],
        [u[4], u[5], u[2]]
    ];
    const RU: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            RU[i][j] = R[i][0] * U[0][j] + R[i][1] * U[1][j] + R[i][2] * U[2][j];
        }
    }
    const out: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            out[i][j] = RU[i][0] * R[j][0] + RU[i][1] * R[j][1] + RU[i][2] * R[j][2];
        }
    }
    // Symmetrise the off-diagonals: round-off from the two products would
    // otherwise make U12 and U21 differ in the last bits, and a reader of the
    // six-number form must see one number per pair.
    return [
        out[0][0],
        out[1][1],
        out[2][2],
        0.5 * (out[0][1] + out[1][0]),
        0.5 * (out[0][2] + out[2][0]),
        0.5 * (out[1][2] + out[2][1])
    ];
}

// ===== Crystal-axis <-> Cartesian =====

function cross(a: number[], b: number[]): number[] {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
}

function dot(a: number[], b: number[]): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * The reciprocal vectors and their lengths, or null for a degenerate cell.
 *
 * The rows of this matrix are `a*`, `b*`, `c*`: `a* = (b x c)/V` and cyclic.
 * They are exactly the map from a Cartesian vector to its fractional
 * components (`f_i = a*_i . r`), and their lengths are the reciprocal cell
 * parameters the CIF-to-Cartesian scaling uses. Both facts come from the same
 * volume, which is why they are computed together: a cell whose volume is
 * indistinguishable from zero has no reciprocal basis to convert with.
 */
function reciprocalVectors(lattice: number[][]): { rows: number[][]; lengths: number[] } | null {
    if (!lattice || lattice.length < 3) return null;
    const a = lattice[0], b = lattice[1], c = lattice[2];
    const bc = cross(b, c);
    const volume = dot(a, bc);
    if (!isFinite(volume) || Math.abs(volume) < CELL_VOLUME_EPSILON) return null;

    const rows = [bc, cross(c, a), cross(a, b)].map(v => v.map(x => x / volume));
    const lengths = rows.map(v => Math.sqrt(dot(v, v)));
    if (!isFinite(lengths[0]) || !isFinite(lengths[1]) || !isFinite(lengths[2])) return null;
    if (!(lengths[0] > 0 && lengths[1] > 0 && lengths[2] > 0)) return null;
    return { rows, lengths };
}

/** M is stored as its three rows; the Cartesian transform needs M^T. */
function transpose3(m: number[][]): number[][] {
    return [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]]
    ];
}

/**
 * A CIF's crystal-axis U_ij as the Cartesian tensor the viewer stores.
 *
 * Returns null for a cell with no usable reciprocal basis; the caller keeps
 * the file's numbers and says so rather than drawing an ellipsoid derived
 * from a division by zero.
 */
export function uCifAsCart(u: AnisoTensor, lattice: number[][]): AnisoTensor | null {
    const rec = reciprocalVectors(lattice);
    if (rec === null) return null;
    const r = rec.lengths;
    const uStar: AnisoTensor = [
        u[0] * r[0] * r[0],
        u[1] * r[1] * r[1],
        u[2] * r[2] * r[2],
        u[3] * r[0] * r[1],
        u[4] * r[0] * r[2],
        u[5] * r[1] * r[2]
    ];
    // U_cart = M^T u_star M, and rotateAnisoTensor(X, M^T) is M^T X (M^T)^T.
    return rotateAnisoTensor(uStar, transpose3(lattice));
}

/**
 * The inverse: a Cartesian tensor as the crystal-axis U_ij a CIF states.
 *
 * The reciprocal rows are the fractional-coordinate map, so the star tensor is
 * `N U N^T` with `N = [a*; b*; c*]`; dividing by the length products returns to
 * the file's own scaling. Without this a saved CIF of a monoclinic or
 * hexagonal cell would carry a tensor rotated relative to the one it was read
 * from.
 */
export function uCartAsCif(u: AnisoTensor, lattice: number[][]): AnisoTensor | null {
    const rec = reciprocalVectors(lattice);
    if (rec === null) return null;
    const r = rec.lengths;
    const uStar = rotateAnisoTensor(u, rec.rows);
    return [
        uStar[0] / (r[0] * r[0]),
        uStar[1] / (r[1] * r[1]),
        uStar[2] / (r[2] * r[2]),
        uStar[3] / (r[0] * r[1]),
        uStar[4] / (r[0] * r[2]),
        uStar[5] / (r[1] * r[2])
    ];
}
