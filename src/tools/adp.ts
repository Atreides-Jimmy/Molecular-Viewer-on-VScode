import { applySymOp, cartToFrac, fracToCart, SymOp } from './crystalGeometry';

/**
 * Anisotropic displacement parameters (ADPs), as a CIF states them.
 *
 * The six numbers are the Cartesian tensor in Angstrom squared:
 *
 *     [U11, U22, U33, U12, U13, U23]
 *
 * The CIF dictionary defines the U_ij in the standard orthogonal frame (a
 * along x, b in the xy plane), which is the frame `cellToLatticeVectors`
 * builds, so a file's numbers need no conversion before the viewer draws them.
 *
 * A site is listed once but a cell holds every symmetry image of it, and an
 * ADP is a second-rank tensor: the image produced by an operation R is
 * `R U R^T`, not `U`. A general-position site under a rotation therefore has
 * one orientation per image, which is what makes an ORTEP drawing of a cubic
 * group readable. `rotateAnisoTensor` is that transform, done in Cartesian
 * space because that is where the tensor is stated.
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
