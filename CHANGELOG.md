# Change Log

## [Unreleased]

### Added

- **Vibration pause/resume** — The Stop button in the Normal Modes panel has been replaced with a Pause/Resume toggle. Clicking Pause halts the animation at the current displacement (atoms stay in place); the button then becomes Resume, which continues the animation from the same phase. Clicking a different frequency row or closing the panel fully stops the vibration as before
- **Vibration mode highlight** — The currently playing (or paused) frequency row is highlighted in amber/yellow to clearly indicate which mode is being demonstrated

### Fixed

- **GJF fixed-atom marker parsing** — Atom lines like `C  -1  -4.67083500  -0.24081200  0.53486900` (where `-1` is a fixed-atom marker between the element and coordinates) were incorrectly parsed, treating `-1` as the x coordinate. The parser now distinguishes markers (plain integers 0/-1/1 without decimal points) from coordinates (which always have decimal points in Gaussian GJF format), correctly handling markers in both positions: after the element symbol or after the three coordinates
- **Gaussian log duplicate frames** — Optimization log files that output both "Input orientation" and "Standard orientation" for each step (e.g., `tetramethylammonium.log`) produced double the expected frames. "Input orientation" of step N+1 shows the same molecular structure as "Standard orientation" of step N, but in a different coordinate system (original input vs principal axes), making them visually redundant in the 3D viewer. The parser now keeps only the first "Input" frame (initial geometry) and skips all subsequent "Input" frames when "Standard" frames are present, giving a clean sequence: initial geometry → step 1 → step 2 → ...
- **Convergence panel auto-open** — The convergence and normal modes panels no longer auto-open when loading a file with optimization/frequency data; instead, the reopen buttons (📈 Convergence / 🎵 Modes) are shown for the user to click when ready
- **Panel scrollbar layout** — Fixed two scrollbar issues in the convergence and normal modes panels: (1) The normal modes panel scrollbar could exceed the panel's visual boundary when many frequencies were listed, because the body had a fixed `max-height:480px` independent of the panel's actual height. Both panels now use flexbox layout (`display:flex;flex-direction:column` on the panel, `flex:1 1 auto;min-height:0` on the body) so the scrollable body always fills the available panel height. (2) The convergence panel showed a horizontal scrollbar before the vertical scrollbar appeared, because `overflow-y:auto` caused the vertical scrollbar to appear late (after the canvas was already drawn at the wider pre-scrollbar width), triggering horizontal overflow. Changed to `overflow-y:scroll` so the scrollbar space is always reserved, and canvas width is measured from `parentElement.clientWidth` (which accounts for the scrollbar) instead of `canvas.clientWidth`

## [0.9.4] - 2026-07-07

### Added

- **Optimization convergence panel** — Gaussian structure optimization log files (`.log`) now display a collapsible convergence panel showing energy, force (max + RMS), and displacement (max + RMS) curves versus optimization step, mirroring GaussView's convergence plot. Each chart is Canvas-rendered with gridlines, threshold lines (dashed amber), data points, and a legend. The panel includes a vertical resize handle (drag up/down to adjust height) and a close/reopen toggle (📈 Convergence button) to reclaim screen space. Energy is extracted from `SCF Done` lines; force/displacement values from the `Maximum Force` / `RMS Force` / `Maximum Displacement` / `RMS Displacement` criteria block after each optimization step
- **Normal mode vibration playback** — Gaussian log files containing a frequency calculation (`freq` keyword or `CalcAll`) now display a collapsible Normal Modes panel listing every harmonic frequency (cm⁻¹) with symmetry label and a ▶ play button. Clicking a frequency animates the molecular vibration by modulating each atom's position with a sine wave along its mass-weighted displacement vector (amplitude normalized to 0.35 Å max, 60 FPS via `requestAnimationFrame` with in-place mesh updates). Vibration playback automatically jumps to the last optimization frame (where the frequency calculation was performed) before animating. Clicking the playing mode again or the ⏹ Stop button stops the animation and restores equilibrium positions. Vibration auto-stops on frame switch, undo, diff entry, or opening any edit modal to prevent stale-position conflicts. Imaginary frequencies (transition states) are marked with an asterisk (`*`)
- **Frequency data parsing** — `logParser.ts` now extracts the full `Harmonic frequencies` block from Gaussian logs: frequencies, symmetry labels, reduced masses, force constants, IR intensities, and per-atom Cartesian displacement vectors. Handles 1–3 columns per frequency block (last block may have fewer modes). For `CalcAll` calculations (where a Hessian is computed at every optimization step), the parser uses the LAST frequency block rather than the first. Displacement vectors are stored as `[atomIndex][xyz]` arrays for vibration animation

### Fixed

- **Vibration animation lag** — The vibration playback was stuttering because every frame disposed and recreated all bond meshes (CylinderGeometry + MeshPhongMaterial + Mesh per bond half), causing 100–300 GPU object allocations/deallocations per frame. Fixed by storing bond metadata (bond index, half, perpendicular offset, original length) on each mesh's `userData` during creation, then updating mesh `position`/`quaternion`/`scale.y` in-place during vibration — zero geometry creation or disposal per frame. Additionally switched from `setInterval` (12.5–20 FPS) to `requestAnimationFrame` (60 FPS) with time-based phase calculation, giving smooth animation with a consistent 3-second cycle period
- **Convergence panel scrollbar** — The convergence panel required scrolling horizontally before the vertical scrollbar appeared, because the Canvas chart had a fixed pixel width (`canvas.style.width=w+'px'`) that didn't shrink when the vertical scrollbar appeared, causing ~15px horizontal overflow. Fixed by using responsive `canvas.style.width='100%'` (automatically adjusts when scrollbar appears) and adding `scrollbar-gutter:stable` to `.opt-body`/`.freq-body` to reserve scrollbar space from the start, preventing layout shift entirely
- **Toolbar and panel overflow on narrow views** — When the VS Code sidebar was open or the editor area was narrow, toolbar buttons (especially frame navigation controls) overflowed beyond the visible area and became inaccessible. Fixed by adding `overflow-x:auto` with a thin 3px scrollbar to the toolbar, `flex-shrink:0` on buttons/separators to prevent compression, and a `@media (max-width:600px)` rule that reduces button padding and font size. Additionally, all floating panels (`#opt-panel`, `#freq-panel`, `#select-panel`, `#crystal-panel`, `#mol-info`) now have `max-width:calc(100% - 16px)` to prevent overflow, the modal dialog uses `min-width:min(320px,calc(100vw - 32px))` to clamp to the viewport, and `#frame-info` truncates long step labels with `text-overflow:ellipsis`

## [0.9.3] - 2026-06-30

### Added

- **VESTA format support** — New parser for VESTA structure files (`.vesta`); reads CELLP (cell parameters), STRUC (fractional coordinates + occupancy), and GROUP (space group) sections; converts fractional coordinates to Cartesian via lattice vectors; supports fractional occupancy for disordered sites; wraps negative fractional coordinates to [0, 1) range
- **Minimum Image Convention (MIC) bond detection** — Crystal bond detection now uses the minimum image convention: for each pair of base atoms, the fractional coordinate difference is wrapped to [-0.5, 0.5] before computing the Cartesian distance, ensuring bonds are always detected at their shortest periodic distance. Each base bond stores a `shift` vector (lattice translation to the nearest image cell)
- **Box Select** — New independent Box Select toolbar button; in this mode, hold the left mouse button and drag a rectangle, then release to select all atoms inside the box via screen-space projection (atoms are projected to screen coordinates and tested against the rectangle); selections accumulate with the same yellow highlighting as the Select Atoms mode, and clear on mode switch
- **Molecular formula display** — The info panel now shows the Hill-system molecular formula (when carbon is present, C first then H then remaining elements alphabetically; without carbon, all elements alphabetically); for crystal files, an additional Unit Cell composition line shows the full unit-cell formula based on the base atoms
- **Undo button** — New Undo toolbar button reverts the most recent edit/delete operation; snapshots (atoms, bonds, baseAtoms, baseBonds deep-copied via JSON) are pushed automatically when an edit modal opens (bond length/angle/dihedral/add atom) and before delete-atom/remove-disorder execute; the stack is capped at 50 entries and resets on file reload
- **VESTA format output** — Save As now supports VESTA (`.vesta`) format for crystal structures; outputs P1 symmetry with identity SYMOP, CELLP cell parameters, and STRUC entries with fractional coordinates converted from Cartesian via the inverse lattice matrix; available only when crystal data is present
- **Select Atoms inline panel** — The Select Atoms toolbar button now toggles an inline panel (top-right corner) instead of a modal dialog; the panel stays visible while in Select Atoms mode, accepts indices/ranges/symbols at any time, and hides automatically when switching to another mode; direct atom clicking for selection still works alongside the panel
- **Keyboard arrow key rotation** — Arrow keys (↑/↓/←/→) now rotate the molecule; left/right rotate around the Y axis, up/down around the X axis; rotation speed is 0.03 rad/frame; keys are ignored when focus is in an input/textarea/select element; works in both normal and diff modes
- **Translate Group / Translate Atom** — Bond Length, Bond Angle, and Dihedral adjustment modals now include a Move dropdown with two options: Translate Group (default, moves the atom plus all connected atoms via BFS) or Translate Atom (moves only the selected atom), matching GaussView's translate group / translate atom behavior
- **Replace Atom** — New toolbar button to replace selected atom(s) with a different element while preserving coordinates; click atoms to select (toggle), click the Replace Atom button again to open the element picker, then confirm to update element and color for all selected atoms (works with crystal base atoms too)
- **Delete keyboard shortcut** — Press the Delete key to quickly delete atoms; if atoms are already selected (via Select Atoms or Box Select), opens the batch delete confirmation; otherwise enters Delete Atom mode; ignored when typing in input fields and in diff mode

### Fixed

- **Cross-boundary bond rendering** — Bonds crossing periodic cell boundaries (e.g. atom at x=0.05 bonded to atom at x=0.95 across the a=1 boundary) were previously rendered as long bonds spanning the entire view. Now they are rendered as short bonds extending outside the view boundary via virtual atom positions, with both forward and reverse split bonds drawn at each edge of the supercell. The `createBond()` function applies the lattice translation vector to compute the virtual position; `rebuildCrystal()` performs a bidirectional propagation pass (forward + reverse) to ensure split bonds appear on all edges of the supercell
- **Cube file bond detection** — Cube files with positive NATOMS (Bohr units) were incorrectly treated as Angstrom, causing atom coordinates to be ~1.89× too large and all interatomic distances to exceed bond length thresholds. Fixed by using the Cube convention: NATOMS > 0 means Bohr, NATOMS < 0 means Angstrom; the voxel-based heuristic is now only used as a secondary check when voxel data is valid
- **VASP/Cube/POSCAR right-click menu missing** — The explorer/context and editor/title/context menus did not include `.vasp`, `.cube`, or bare `POSCAR`/`CONTCAR` filenames in their `when` clauses. Added these extensions and a case-insensitive regex match for bare POSCAR/CONTCAR filenames
- **Crystal bond detection on initial load** — When opening crystal files (CIF/VASP/Cube/VESTA), bonds were not detected on initial load because `rebuildCrystal()` was not called before `rebuildScene()`. Additionally, the extension's `ensureBonds()` (which uses plain distance without MIC) was previously applied to crystal structures and its results copied to `baseBonds`, overriding the webview's MIC detection. Fixed by: (1) skipping `ensureBonds()` for crystal files in the extension, and (2) calling `rebuildCrystal()` before `rebuildScene()` on initial load when crystal data is present
- **GPU memory leak / rendering performance regression** — `rebuildScene()` and `updateScenePositions()` rebuilt meshes without disposing the old geometry/material/texture, causing GPU resources to accumulate within a session and making both molecular and crystal files progressively laggy. A unified `disposeMesh()` helper now releases the old resources before removal. In addition, `animate()` was switched to dirty-flag rendering driven by a global `needsRender` flag: `requestAnimationFrame` still ticks, but rendering only proceeds when `needsRender` is set (on interaction, scene rebuild, highlight, zoom, resize); idle frames no longer re-render, eliminating continuous GPU usage
- **VESTA symmetry expansion** — VESTA files with symmetry operations (SYMOP section) were previously parsed but the symmetry operations were discarded, leaving only the asymmetric unit atoms in `baseAtoms`. This caused MIC bond detection to find phantom long bonds between distant asymmetric-unit atoms, resulting in stretched bonds spanning the cell, atom-bond misalignment, duplicate ghosting, and overall structural distortion. Fixed by parsing SYMOP into `{t, R}` pairs and expanding each atom through all symmetry operations (with `wrapFrac` and 0.25 Å² dedup), mirroring the CIF symmetry expansion logic
- **Rotate/zoom/pan stuttering for crystal files** — The `updateAxesIndicator()` function rebuilt its SVG via `innerHTML` string concatenation on every `mousemove` event (~200 Hz), causing DOM thrash and frame-by-frame stuttering when dragging crystal structures. Fixed by moving the call from `updateTransform()` to `animate()` after the `needsRender=false` check, throttling it to the animation frame rate (~60 Hz) and only when a render is actually scheduled
- **Add Atom overlap on repeated additions** — Adding multiple atoms to the same anchor produced identical directions every time (negated average of existing neighbor vectors), causing all new atoms and bonds to overlap. Fixed by replacing the deterministic negated-average with a max-min-angle candidate selection over a Fibonacci sphere of 50 unit directions: the new atom is placed at the direction that maximizes the minimum angular separation from all existing neighbors, naturally producing linear/trigonal/tetrahedral/octahedral coordination geometry
- **Add Atom missing undo snapshot** — The Add Atom OK handler modified `CRY.baseAtoms`/`MD.atoms` without first calling `pushUndo()`, so the Undo button could not revert atom additions. Fixed by calling `pushUndo()` at the start of the OK handler
- **Crystal bond rendering anomaly (non-orthogonal lattices)** — VESTA and CIF files with non-orthogonal lattices (α/β/γ ≠ 90°) displayed numerous long bonds piercing the cell boundary with dangling ends. Root cause: `computeCrystalInv()` computed A⁻¹ (inverse of the lattice matrix A, where lattice vectors are rows), but `cartToFrac` requires (Aᵀ)⁻¹ = (A⁻¹)ᵀ because `fracToCart` computes cart = Aᵀ·frac. The missing transpose caused ~2 Å errors in fractional coordinates, wrong MIC shift vectors, normal bonds between distant atoms (up to 17.9 Å), and dangling split bonds. Fixed by swapping the 6 off-diagonal elements in both `computeCrystalInv()` (webview) and `cartToFrac()` (cifParser.ts). The bug was invisible for orthogonal lattices where A⁻¹ = (A⁻¹)ᵀ
- **Crystal bond rendering anomaly (atoms outside unit cell)** — VESTA files containing atoms with Cartesian coordinates outside the unit cell (fractional coordinates outside [0, 1)) still showed a few long bonds even after the transpose fix. Root cause: `detectCrystalBaseBonds()` computed MIC shift vectors from unwrapped fractional coordinates, while `rebuildCrystal()` wrapped atoms to [0, 1) before placing them in cells. This inconsistency produced shift vectors pointing to non-existent cells, creating split bonds with virtual endpoints ~40 Å away. Fixed by wrapping fractional coordinates to [0, 1) in `detectCrystalBaseBonds()` before computing differences, making shift vectors consistent with the wrapped positions
- **Bond length slider drift after Add Atom** — Adjusting the bond length of a newly added bond via slider input produced visible drift between the entered value and the displayed bond length (e.g. entering 0.9 showed a very short bond, re-checking showed 0.4). Root cause: `applyBondLength()` computed `curLen` from `MD.atoms` (current/modified coordinates) instead of `originalCoords` (snapshot at modal open), causing floating-point error accumulation across slider input events. Fixed by using `originalCoords` consistently for `curLen`, matching the pattern already used in `applyBondAngle()` and `applyDihedral()`
- **Supercell boundary mirror images** — Setting non-integer supercell boundaries (e.g. a=-0.4 to 0.4) generated mirror images: atoms at fractional coordinate 0.6 appeared wrapped at -0.4 on the opposite side, and deleting one image deleted both. Root cause: `rebuildCrystal()` used `Math.floor(aMin)`/`Math.ceil(aMax)` for cell iteration range, producing cells -1/0/1 for a -0.4 to 0.4 range. Fixed by switching to `Math.ceil(aMin)`/`Math.floor(aMax)`, which only iterates cell 0 for the same range. Integer boundaries (0/1, -1/1) are unaffected since ceil/floor produce the same result as floor/ceil
- **Finite molecule truncation at supercell boundary** — Finite molecules (molecules not extending periodically) straddling the supercell boundary were truncated, showing broken bonds at the edge. Fixed by adding a finite-molecule completion pass in `rebuildCrystal()`: connected components in `baseBonds` where all bonds have zero shift vectors are identified as finite; if any atom of a finite molecule is within the boundary, all atoms of that molecule are shown (extending beyond the boundary). Infinite periodic chains are still truncated at the boundary via split bonds
- **Select Atoms panel layout compression** — Clicking the Select Atoms button showed an input panel that caused the molecular view to narrow and the molecule to be horizontally compressed. Root cause: the `#select-panel` div was a sibling of `#container` in the body flex layout, causing it to affect layout calculations. Fixed by moving the panel inside `#container` so its `position:absolute` is relative to the container, preventing any impact on the molecular view layout
- **Select Atoms toggle selection** — In Select Atoms mode, clicking an already-highlighted atom did nothing instead of deselecting it. Fixed by adding toggle logic in `selectAtom()`: clicking a selected atom in selectAtoms (or replaceAtom) mode now removes it from the selection, and clicking again re-selects it
- **GJF save bond corruption** — After editing and saving a GJF file, the connect section was completely scrambled: bonds were auto-detected (geometric) instead of preserving explicit bonds, and stale connect data from `afterConnectContent` was appended, producing extra lines referencing non-existent atoms (e.g. 164 connect lines for 156 atoms). Root cause: (1) `gjfParser.ts` validation rejected the ENTIRE connect section when any atom reference exceeded the atom count, storing it all as `afterConnectContent`; (2) `ensureBonds()` then replaced bonds with geometric auto-detected bonds (different from explicit bonds, especially long-range bonds); (3) `parseConnectLine` stored bonds without normalization (atom1 = line atom, not min). Fixed by: relaxing parser validation to parse valid lines and skip invalid ones; normalizing bonds to `{atom1: min, atom2: max}` in `parseConnectLine`; filtering connect-like lines from `afterConnectContent` on save to prevent stale data from being appended
- **Translate Group scaling instead of rigid translation** — When adjusting bond length with Translate Group mode, all bonds in the connected group shortened/stretched synchronously instead of the group moving as a rigid unit. Root cause: `applyBondLength()` used `scale = targetLen/curLen` and applied it to all atoms relative to the fixed point (scaling), instead of computing a translation vector. Fixed by computing the moved atom's new position, deriving a translation vector `t = newPos - oldPos`, and applying the same translation to all movable atoms, preserving their relative positions (rigid translation)
- **Supercell boundary sudden appearance** — Adjusting the a/b/c lower bound from 0 toward -1 showed no new content until -1.0, then suddenly displayed the entire adjacent cell. Root cause: `Math.ceil(aMin)`/`Math.floor(aMax)` skipped boundary cells for non-integer values (e.g. ceil(-0.99)=0, only cell 0 iterated). Fixed by reverting to `Math.floor(aMin)`/`Math.ceil(aMax)`, which iterates the boundary cell and lets the fractional position filter progressively reveal atoms as the boundary expands. The prior mirror-image concern is moot since `propagateToAllCells()` was removed (deletions no longer propagate to periodic images)
- **Initial load shows Select panel and keyboard unresponsive** — On first opening a file, the Select Atoms panel was briefly visible and arrow key rotation did not respond until another toolbar button was clicked. Root cause: no explicit focus management on initialization — the browser focus could be on an input element, causing the keyboard event filter to ignore arrow keys. Fixed by setting `container.tabIndex=0` and calling `container.focus()` on init, and explicitly hiding the select panel before first render
- **Select panel layout impact** — The Select Atoms input panel still affected the molecular view height when shown/hidden. Root cause: the panel was inside `#container` with `position:absolute`, which could influence flex layout. Fixed by moving the panel to body level with `position:fixed;bottom:40px;left:8px`, making it a truly floating overlay that never affects the molecular view dimensions; repositioned from top-right to bottom-left to avoid overlapping the molecular info panel
- **Gold focus border on initial open** — On first opening a file, a gold/yellow outline appeared around the 3D view container, then never reappeared after switching modes. Root cause: `container.focus()` on init applied the browser default focus outline. Fixed by adding `outline:none` to the `#container` CSS rule
- **Toggle selection not working** — Clicking an already-highlighted atom in Select Atoms mode did not deselect it. Root cause: the toggle logic added in a prior fix was lost. Fixed by re-applying the toggle block in `selectAtom()`: clicking a selected atom in selectAtoms or replaceAtom mode removes it from the selection, and clicking again re-selects it
- **Selection cleared between Select Atoms and Box Select** — Switching between Select Atoms and Box Select modes cleared all highlighted atoms, forcing the user to re-select. Fixed by modifying `setMode()` to preserve `selectedAtoms` when transitioning between selectAtoms↔boxSelect; selection is only cleared when switching to any other mode
- **Replace Atom button missing** — The Replace Atom toolbar button was not visible, preventing access to the feature. Root cause: the button HTML and `MODE_INFO` entry were lost from a prior edit. Fixed by re-adding the `<button data-mode="replaceAtom">` toolbar button and the corresponding `MODE_INFO` entry
- **Translate Group fixed atoms drifting visually** — During bond length/angle/dihedral adjustment with Translate Group mode, the fixed atoms/part appeared to shift in the view even though their coordinates didn't change. Root cause: `updateScenePositions()` recalculated the geometric center of all atoms on every call; when the movable group translated, the center shifted, causing fixed atoms to appear displaced relative to the new center. Fixed by adding a `keepCenter` parameter to `updateScenePositions()`; the three live-preview functions (`applyBondLength`, `applyBondAngle`, `applyDihedral`) now call `updateScenePositions(true)` to reuse the existing center, keeping fixed atoms visually anchored in place
- **abc boundary negative direction still empty** — Adjusting the a/b/c lower bound from 0 toward -1 still showed no new content until exactly -1.0. Root cause: the `Math.floor(aMin)`/`Math.ceil(aMax)` fix from a prior release had regressed back to `Math.ceil(aMin)`/`Math.floor(aMax)`, which skipped boundary cells for negative non-integer values (e.g. ceil(-0.99)=0). Fixed by restoring `Math.floor(aMin)`/`Math.ceil(aMax)` so the boundary cell is iterated and the fractional position filter progressively reveals atoms as the boundary expands in the negative direction
- **Select panel still affecting view width** — The Select Atoms panel still influenced the molecular view width and appeared at the wrong position despite multiple relocation attempts. Root cause: all CSS was defined as a long inline `style` attribute on the div, which was unreliable in the webview rendering environment. Fixed by moving all CSS to a `#select-panel` rule in the `<style>` block (with child rules for `.sp-title`, `input[type=text]`, `.sp-btns`), removing all inline styles from the HTML, matching the proven pattern used by `#crystal-panel` and `#diff-panel` which have zero layout impact
- **Diff label overlapping molecular info** — In diff mode, the "Original: ..." label (`#diff-label`) was positioned at `top:8px;left:8px`, identical to the molecular info panel (`#mol-info`), causing the two to overlap and obscure the molecular formula, charge, and electron count. Fixed by repositioning `#diff-label` to `left:25%;transform:translateX(-50%)` (top-center of the left viewport half) and `#diff-label-right` to `left:75%;transform:translateX(-50%)` (top-center of the right viewport half), keeping both labels visible without overlapping any corner panels
- **Selection info bar overflow covering toolbar** — When selecting a large number of atoms (e.g. via Box Select), the selection info text in the status bar (`#selection-info`) grew to multiple lines, overflowing the 24px-tall `#status-bar` and covering the toolbar buttons above. Fixed by adding `overflow:hidden` to `#status-bar` and `white-space:nowrap;overflow:hidden;text-overflow:ellipsis;flex:1;min-width:0` to `#selection-info` for CSS-level truncation; additionally, a `formatAtomList()` helper now limits the displayed atom list to the first 15 atoms followed by `... (N atoms)` instead of listing all selected atoms
- **Floating panel collision detection system** — In diff mode, the diff labels (`#diff-label` at 25% / `#diff-label-right` at 75%) still overlapped the molecular info panel and crystal panel in narrow viewports, since the fixed percentage positions didn't account for other panels' widths. Fixed by implementing a general `layoutPanels()` collision avoidance system: all floating panels (mol-info, crystal-panel, diff-panel, select-panel, diff-label, diff-label-right, axes-indicator) are assigned a priority order; on each mode switch, diff enter/exit, panel toggle, mouseup (rotation/pan/box-select completion), wheel zoom, or window resize, the system resets inline overrides, then for each visible panel (in priority order) checks overlap with all higher-priority panels via `getBoundingClientRect()` and resolves conflicts by pushing the lower-priority panel down/up; if no space is available, it progressively shrinks the panel (reducing max-width with CSS ellipsis truncation, then font size down to 8px); if the current panel can't be shrunk further, it shrinks the conflicting higher-priority panel; multiple passes handle cascading overlaps. The system also computes the molecule's screen-space bounding box via `getMoleculeScreenBox()` (projects the 3D bounding box corners through the camera, handling both single-viewport and split-viewport diff modes with camera state save/restore) and treats it as an immovable collision object — panels that overlap the molecule are pushed towards the nearest screen edge (for top-center panels like diff labels) or progressively shrunk, ensuring the 3D structure view is never obscured. Diff labels also gained `max-width:40%;white-space:nowrap;overflow:hidden;text-overflow:ellipsis` CSS for graceful text truncation when shrunk
- **Box select selects all atoms after diff collision** — After entering diff mode and triggering the panel collision detection (when the molecule overlapped a diff label), subsequent Box Select operations selected every atom in the view instead of only those inside the dragged rectangle. Root cause: `getMoleculeScreenBox()` saved the camera state, then inside a `try` block called `camera.updateProjectionMatrix()` (which updates both `projectionMatrix` AND `projectionMatrixInverse`), but the `finally` block only restored `projectionMatrix` via `camera.projectionMatrix.copy(sProj)`, leaving `projectionMatrixInverse` with the diff viewport's value. Since `Raycaster.setFromCamera()` and `Vector3.project()` ultimately depend on the camera matrices being consistent, the stale `projectionMatrixInverse` corrupted subsequent projections. Additionally, a throttled `layoutPanels()` call (300ms) inside the `animate()` loop ran `getMoleculeScreenBox()` during active mouse interactions, re-corrupting the camera mid-drag. Fixed by: (1) replacing `camera.projectionMatrix.copy(sProj)` with `camera.updateProjectionMatrix()` in the `finally` block so both matrices are recomputed from the restored `aspect`, removing the `sProj` variable entirely; (2) removing the throttled call from `animate()` and instead invoking `layoutPanels()` at interaction endpoints only — mouseup (after rotation/pan/box-select completes) and wheel (after zoom) — so the camera is never mutated during an in-progress drag
- **Box select atom mismatch in diff mode** — In diff mode, Box Select selected atoms that did not visually match the dragged rectangle — the projected atom positions felt offset or wrong relative to what was displayed. Root cause: the Box Select `mouseup` handler was the only interaction handler that did NOT account for the split-viewport diff layout. It projected `atomMeshes` (left-viewport atoms) using whatever camera state was left over from the last hover/click (which could be the right-viewport camera with `aspect=(w-halfW)/h` and `position=(0,0,diffCamDist)`), and mapped NDC coordinates to the full canvas width (`rect.width`) instead of the left viewport width (`halfW`). This produced incorrect screen coordinates, so atoms appeared inside the box when they shouldn't (or vice versa). The handler also ignored `diffAtomMeshes` entirely, making right-viewport box selection impossible. Fixed by adding a dedicated `diffMode` branch that mirrors the proven pattern from `getClickedAtom()` and the hover handler: for the left viewport, set up `camera.aspect=halfW/h`, `camera.position=(0,0,camDist)`, `updateProjectionMatrix()` + `updateMatrixWorld(true)`, project `atomMeshes` with NDC mapped to `[rect.left, rect.left+halfW]`, set `diffActiveSide='left'`, and call `selectAtom()` for hits; for the right viewport, set up `camera.aspect=(w-halfW)/h`, `camera.position=(0,0,diffCamDist)`, project `diffAtomMeshes` with NDC mapped to `[rect.left+halfW, rect.right]`, set `diffActiveSide='right'`, and call `selectAtom()` for hits. A box spanning both viewports now correctly selects atoms from both sides using each viewport's correct camera and screen mapping
- **Del key not working in diff mode** — The Delete key was explicitly disabled in diff mode via `if(e.key==='Delete'&&!diffMode)`, preventing atom deletion during structural comparison. Fixed by removing the `!diffMode` restriction so Del works in all modes
- **Diff highlights stale after structural edits** — After deleting, adding, or replacing atoms in diff mode, the diff highlights and mapping became stale (referenced old atom indices), causing incorrect or missing highlights. Fixed by creating `refreshDiff()` which re-detects bonds from the current `MD.atoms`, re-computes the atom mapping via `findAtomMapping()`, and calls `recomputeDiff()` to update highlights and the diff panel; `rebuildScene()` now auto-calls `refreshDiff()` at the end when in diff mode, covering all structural modification paths (delete, add, replace, undo, frame switch)
- **`checkSelectionComplete` missing modes in diff** — The diff-mode branch of `checkSelectionComplete` only handled `bondLength`, `bondAngle`, `dihedral`, and `selectAtoms` modes; `deleteAtom`, `addAtom`, `boxSelect`, and `replaceAtom` modes were missing, causing the corresponding modals/info to not appear in diff mode. Fixed by adding all missing mode handlers
- **exitDiff GPU memory leak** — `exitDiff()` removed `diffPivot` from the scene but did not dispose the geometry/material of its child meshes (diff atoms and bonds), causing GPU resource accumulation across repeated diff enter/exit cycles. Fixed by traversing `diffPivot` and disposing all geometry/material (and material.map) before removal
- **refreshDiff skeleton-differ branch incomplete** — When `refreshDiff()` detected that skeletons differed after an edit (mapping failed), it cleared left-side highlights but not right-side (`diffAtomMeshes`, `diffBondMeshes`) highlights, leaving stale orange highlights on the right viewport. Fixed by clearing all diff highlight flags and resetting emissive/color on both sides in the else branch
- **CRY.symmetryOps crash on CIF save** — `doSave()` accessed `CRY.symmetryOps.forEach(...)` without a null check; crystal data loaded from some formats could have `symmetryOps` undefined, causing a crash on CIF export. Fixed with `(CRY.symmetryOps||['x, y, z'])` fallback
- **Replace Atom / Add Atom double undo snapshot** — `showReplaceAtomModal` and `showAddAtomModal` called `saveOriginal()` (which internally calls `pushUndo()`) at the top, then their OK handlers called `pushUndo()` again, pushing two undo snapshots for one operation. Fixed by removing the `saveOriginal()` call (neither modal uses `originalCoords` for live preview) so only one `pushUndo()` fires in the OK handler
- **Delete/Batch Delete undo pollution on cancel** — `showDeleteAtomModal` and `showBatchDeleteModal` called `pushUndo()` before showing the modal; if the user cancelled, the undo stack retained a stale snapshot that wasted an undo press. Fixed by moving `pushUndo()` into the OK click handler so the snapshot is only pushed when the deletion is confirmed
- **Slider modal undo pollution on cancel** — Bond Length / Bond Angle / Dihedral modals called `saveOriginal()` (→ `pushUndo()`) at the top; cancelling after moving the slider left a stale undo snapshot even though `restoreOriginal()` reverted the coordinates. Fixed by calling `undoStack.pop(); updateUndoBtn()` in the Cancel handlers to remove the stale snapshot
- **Frame switch undo confusion** — Switching frames replaced `MD.atoms` with the new frame's atoms but left the undo stack populated with snapshots from the previous frame; pressing Undo after a frame switch would restore the wrong frame's state. Fixed by clearing the undo stack and resetting selection on frame switch
- **Selection cleared when entering Replace Atom** — Switching from Select Atoms or Box Select to Replace Atom mode cleared all selected atoms, breaking the workflow of selecting atoms then switching to Replace Atom to replace them. Fixed by adding `replaceAtom` to the `preserveSel` set in `setMode()` so selection is preserved across selectAtoms↔boxSelect↔replaceAtom transitions
- **Select panel Clear button didn't sync diff side** — The Select Atoms panel's Clear button only cleared `selectedAtoms` (left side), not `diffSelectedAtoms` (right side), leaving right-viewport atoms highlighted after Clear in diff mode. Fixed by clearing both arrays
- **formatAtomList crash on stale indices** — `formatAtomList` accessed `atoms[i].element` without bounds checking; stale selected indices after structural edits could cause an undefined access. Fixed with a null check (`a?a.element:'?'`)
- **VASP element detection rejected C/S/D elements** — `isElementTypeLine` had a guard `!/^[DCSdc]/i.test(line6)` intended to avoid misidentifying coordinate keywords (Direct/Cartesian/Selective) as element lines, but it also rejected legitimate element lines starting with C, S, or D — including Cs, Cd, Ca, Cl, Co, Cr, Cu, Sr, Sc, Si, Sn, Sm, Sb, Dy. Fixed by replacing the broad prefix guard with a narrow check that only rejects single-token numeric lines; the existing `every(t => /^[A-Za-z][a-z]?$/)` check already correctly distinguishes element lines from coordinate lines
- **GJF fixed-atom marker misdetected for integer coordinates** — `parseAtomLine` flagged any numeric value of 0, 1, or -1 as a fixed-atom marker, but integer coordinates (e.g. `C 0 1 0 0`) were misidentified, stripping legitimate coordinate values. Fixed by using column-count-based detection: for 4-column lines, all values are coordinates; for 5+ column lines, the first 3 numeric values after the element are coordinates and trailing integers are markers
- **VESTA STRUC parser assumed 2 lines per atom** — The STRUC block parser unconditionally did `i += 2` after each atom, assuming a detail line always followed; VESTA files with single-line STRUC entries skipped legitimate atom lines. Fixed by peeking at the next line and only skipping it if it does NOT look like a valid atom line (wrong column count or non-element symbol)
- **PDB CONECT only parsed 4 connections** — The CONECT record parser read columns 7–31 (1 center + 4 bonded atoms), missing columns 32–46 which can hold 3 additional bonded atoms per the PDB specification. Fixed by extending the parse range to column 46, allowing up to 7 bonded atoms per CONECT line
- **ORCA A.U. coordinate block never parsed** — The ORCA `.out` parser's A.U. (Bohr) coordinate block required `parts.length >= 7` and read `parts[5], parts[6], parts[7]`, but ORCA's A.U. format has only 6 columns (index, element, nuclear charge, x, y, z), so `parts[7]` was always undefined and the block was silently skipped. Fixed by finding the element symbol (first non-numeric token) and taking the last 3 numeric values as coordinates, making the parser robust to column count variations
- **TCL parser didn't support quoted file paths** — `mol new` path extraction used `(\S+)` which split paths containing spaces (e.g. `mol new "C:\My Folder\file.xyz"`). Fixed by matching quoted strings (`"..."` or `'...'`) first, falling back to `\S+` for unquoted paths
- **openCustomDocument no error handling** — `parseFile()` and `parseLogFile()` could throw on malformed input, crashing the editor open flow with an unhandled exception. Fixed by wrapping the parse dispatch in a try/catch that returns a placeholder `MolecularData` with the error message in the title
- **`.mol`/`.sdf` registered without parser** — The package.json customEditors selector, context menus, and open/diff dialogs registered `.mol` and `.sdf` extensions, but no MDL Mol parser existed, causing these files to fall through to `tryAutoParse` → `parseXyz` (garbage output). Fixed by removing `.mol`/`.sdf` from all registration points (the `.mol` save format is preserved since a writer exists)
- **MOL2 `am` bond order treated as double** — MOL2 `am` (amide) bonds were mapped to order 2 (double), but amide bonds should be order 1 (single with partial double character from resonance, not an actual double bond). Fixed by mapping `am` → 1; also added explicit `du` (dummy) → 1 and `ar` (aromatic) → 1.5 handling
- **ORCA inp `%coords` units ignored** — The ORCA inp parser skipped `units` lines in `%coords` blocks without parsing them; files using `units bohr` had coordinates interpreted as Angstrom (~1.89× too small). Fixed by tracking `inCoordsBohr` and converting Bohr → Angstrom when the `units bohr` directive is present
- **coordParser `$chrg=1` not parsed** — Turbomole `$chrg` and `$spin`/`$mult` directives with `=` syntax (e.g. `$chrg=1`) were not parsed because the code used `line.split(/\s+/)[1]` which returns undefined for `=`-delimited values. Fixed by first trying a `=\s*(-?\d+)` regex match, falling back to whitespace splitting
- **PDB missing actinide elements** — The PDB `ATOMIC_NUMBER_MAP` stopped at Rn (86), so actinide-containing PDB files (Th, U, Pu, etc.) had their elements misidentified (e.g. "Th" → "T"). Fixed by extending the map through Cf (98)
- **Webview bond detector missing covalent radii** — The webview's `detectBondsFromAtoms` and `detectCrystalBaseBonds` covalent radii map stopped at I (53), using a 1.5 Å fallback for Rb–Rn. This caused inaccurate bond detection for metals, lanthanides, and heavy elements. Fixed by extending the CR2 map to match the backend `bondDetector.ts` COVALENT_RADII (through Rn, 86)

### Changed

- **File association** — `package.json` customEditors selector and explorer/editor context menus now register `.vesta` extension; context menus also match bare `POSCAR`/`CONTCAR` filenames (case-insensitive) via regex
- **Crystal bond detection architecture** — Crystal bonds are now detected exclusively in the webview using MIC via `detectCrystalBaseBonds()`, not in the extension via `ensureBonds()`; this ensures periodic boundary conditions are correctly handled for all crystal formats (CIF, VASP, Cube, VESTA)
- **Save robustness** — `doSave()` now guards against empty atom arrays and NaN coordinates at the entry point, checks for a degenerate crystal lattice (zero determinant) before generating crystal formats, wraps the entire format-generation block in a try/catch that surfaces an informative alert on failure, and reuses the global `CRY_INV` inverse matrix for the CIF and VASP branches instead of recomputing the 3×3 adjugate matrix in each branch
- **Batch delete from selection** — When atoms are already highlighted (via Select Atoms or Box Select), clicking the Delete Atom button now opens a batch delete confirmation showing all selected atom names; confirming removes all selected atoms at once using the `oldToNew` remapping pattern (crystal) or descending-splice with ascending reindex (non-crystal). The original click-to-delete flow (enter delete mode, click one atom) is preserved when no selection exists or when already in delete mode
- **Adaptive level-of-detail** — Sphere segment count now scales with atom count: ≤200 atoms uses 32×24 (full quality, unchanged), 201–1000 uses 20×16, >1000 uses 12×10. This keeps small molecules crisp while ensuring large crystal structures (e.g. VESTA-expanded MOFs with 1000+ atoms) remain interactive
- **Crystal modifications no longer propagated** — Atom edits (bond length, angle, dihedral adjustments) in crystal structures now modify only the specific atom being edited, not all periodic images of the same base atom. The `propagateToAllCells()` function has been removed; each periodic image is independently editable. This matches the expectation that crystal file modifications are local, not periodic

## [0.9.2] - 2026-06-30

### Added

- **CIF crystal structure support** — New parser for Crystallographic Information Files (`.cif`); reads cell lengths/angles, space group symmetry operators, and atom site fractional/Cartesian coordinates; applies symmetry operations to generate the complete unit cell; supports fractional occupancy with partial coloration rendering; outputs fractional coordinates in minimal CIF format on save
- **VASP POSCAR/CONTCAR support** — New parser for VASP structure files (`.vasp`, `POSCAR`, `CONTCAR`); handles scale factor, lattice vectors, element types/counts, selective dynamics flag, and both Direct (fractional) and Cartesian coordinate systems; converts fractional coordinates to Cartesian via lattice vectors; outputs Direct (fractional) coordinates on save
- **Gaussian Cube support** — New parser for Gaussian Cube volumetric files (`.cube`); reads voxel vectors and atomic positions; auto-detects Bohr/Angstrom units; robustly handles corrupted voxel vectors (e.g. `1.#INF00`) by inferring an orthogonal lattice from atomic coordinate ranges with padding; supports atom-only cube files without volumetric data
- **Crystal structure viewer** — Comprehensive crystal visualization system: initial unit cell display with cell wireframe; rotating axes indicator (bottom-left corner) that synchronizes with the main view; supercell boundary controls with 6 input boxes (a/b/c min/max) supporting fractional values (step 0.1) for partial cell selection; periodic boundary modification propagation via `propagateToAllCells()` so atom edits apply to all periodic images; "Remove Disorder" button to clean up disordered sites (occupancy < 0.5 removed, > 0.5 set to 1, = 0.5 kept)
- **Save as VASP/Cube** — Export modified crystal structures in VASP POSCAR format (Direct fractional coordinates grouped by element) or Gaussian Cube format (with minimal 2×2×2 volumetric grid); both available only when crystal data is present

### Fixed

- **Crystal deformation bug when adjusting a/b/c boundaries** — Adjusting supercell a/b/c boundary values caused the crystal structure to appear stretched/distorted. Root cause: the 3×3 inverse matrix used for Cartesian→fractional coordinate conversion contained 3 incorrect elements in the adjugate matrix formula (affects `rebuildCrystal()`, CIF save in `doSave()`, and `cartToFrac()` in cifParser.ts). The incorrect elements were `inv[0][1]` (used `lv[0][0]` instead of `lv[0][1]`), `inv[0][2]` (completely wrong formula), and `inv[1][2]` (completely wrong formula). Fixed by applying the correct adjugate (cofactor transpose) formula in all three locations. The bug was previously unnoticed because CIF parsing uses the forward fractional→Cartesian transform, and diagonal lattices happen to produce correct results despite the matrix errors

### Changed

- **Save dialog filters** — Save As dialog now includes VASP POSCAR and Gaussian Cube format filters alongside all existing formats; diff file picker now accepts VASP and Cube files
- **File association** — `package.json` customEditors selector and explorer context menu now register `.vasp`, `POSCAR`/`CONTCAR` (case-insensitive), and `.cube` extensions

## [0.9.1] - 2026-06-23/28

### Added

- **Adjustable diff threshold slider** — A percentage slider (0–20%, step 0.5%, default 2%) in the diff results panel lets you interactively adjust the sensitivity for what counts as a "difference" in bond lengths, bond angles, and dihedral angles; the diff list updates in real time as you drag the slider

### Fixed

- **Select Atoms space-separated parsing** — Space-separated atom input (e.g. `1 2 3`) only selected the first atom because regex escape sequences (`\s`, `\d`) inside the webview template literal were stripped during evaluation, causing the split regex to match the letter `s` instead of whitespace. Range matching (e.g. `3-10`) was also broken. Fixed by double-escaping backslashes (`\\s`, `\\d`) so the browser receives the correct regexes

### Changed

- **GJF save format** — Atom coordinates in saved GJF files now use 8 decimal places (zero-padded) with GaussView-style alignment (element left-justified to 2 chars, each coordinate right-justified to 17 chars), replacing the previous 6-decimal fixed-separator format that produced misaligned columns
- **Skeleton matching ignores bond order** — The diff skeleton isomorphism check now considers only connectivity (whether atoms are bonded), not bond order. Previously, if corresponding atoms in the two molecules had different bond orders (e.g. single vs. double), the skeleton would be incorrectly considered different. Now all bonds are treated as single bonds for skeleton matching purposes; bond order differences do not prevent conformation diff comparison

### Improved

- **Select Atoms accumulation** — Clicking the Select Atoms toolbar button while already in select mode now reopens the input dialog (previously required switching to another mode first). Selections accumulate across multiple dialog sessions within the same mode — atoms stay highlighted until switching to another mode (View, Add, etc.). Previously each dialog session replaced the previous selection

## [0.9.0] - 2026-05-22

### Added

- **Molecular structure diff** — New "Diff" toolbar button that compares the currently open structure against another molecular file (any supported format, no need for matching formats); if the second file is an optimization LOG/OUT, a QuickPick lets you select which frame to compare
- **Skeleton comparison** — Uses graph isomorphism via BFS-based atom matching (element + degree + neighbor-element signature) to determine if two molecules share the same connectivity skeleton; if skeletons differ, shows a clear "skeletons differ" message with atom/bond counts
- **Conformation diff** — When skeletons match, computes bond length differences (threshold 0.05 Å), bond angle differences (threshold 2°), and dihedral angle differences (threshold 5°); displays a sorted diff panel listing all differences with values from both structures
- **Side-by-side rendering** — In diff mode, the 3D view splits into left (original) and right (comparison) viewports using Three.js scissor/viewport; both molecules rotate synchronously; differing atoms and bonds are highlighted in orange on both sides
- **Diff info panel** — Floating panel showing all differences grouped by type (bonds, angles, dihedrals), sorted by magnitude, with atom labels and values from both structures; closable via × button
- **Diff exit** — Clicking the Diff button again exits diff mode and restores the single-molecule view

## [0.8.2] - 2026-05-22

### Fixed

- **Gaussian LOG external program coordinate filtering** — When Gaussian calls external programs during optimization (e.g., ONIOM), the LOG file may contain coordinate blocks where atomic numbers are written as element symbols (e.g., `C` instead of `6`) instead of numeric atomic numbers. Previously, these lines were parsed with `parseInt` returning `NaN`, causing atoms to be misidentified as element `X`. Now, lines with non-numeric atomic numbers in `Standard orientation:` and `Input orientation:` blocks are skipped, ensuring only valid Gaussian-native coordinate entries are read

## [0.8.1] - 2026-05-22

### Added

- **VMD TCL script support** — New parser for VMD TCL visualization scripts (`.tcl`); reads `mol new <filename>` to locate the referenced molecular structure file (supports relative and absolute paths); parses `mol color ColorID N` + `mol selection "index ..."` pairs to extract atom color group assignments; maps VMD ColorID (0-32) to hex colors (blue, red, gray, orange, yellow, tan, silver, green, white, pink, cyan, purple, etc.); atoms in each group are rendered with their assigned VMD color instead of the default CPK element color
- **TCL file integration** — Opening a `.tcl` file automatically resolves and loads the referenced molecular file, applies the color groups, and displays the structure with per-group coloring; supports all molecular file formats as source (PDB, XYZ, GJF, etc.); works with both relative and absolute file paths; if the source file cannot be found, shows an informative error message
- **AtomGroup data type** — New `AtomGroup` interface in types.ts with `colorId`, `color`, and `indices` fields; `MolecularData` now has an optional `atomGroups` field for color group information

### Changed

- **Package name** — Changed from `molecular-viewer` to `Molecular-Viewer` (capitalized) to allow re-publishing after accidental deletion of the original extension
- **.gitignore** — Added `test/` folder to git ignore list

## [0.7.2] - 2026-04-15

### Added

- **PDB file support** — New parser for Protein Data Bank files (`.pdb`, `.ent`); reads ATOM and HETATM records with fixed-column PDB format parsing; element identification from columns 77-78 (element field) or inferred from atom name with 2-letter element detection (e.g. FE→Fe, CL→Cl); CONECT records for explicit bond connectivity; handles duplicate atom serial numbers; title from TITLE records
- **MOPAC file support** — New parser for MOPAC input/output files (`.mop`, `.mopac`, `.dat`); reads internal coordinate format (`ELEM x 1 y 1 z 1`); supports both element symbols and atomic numbers; auto-detects Å/Bohr units from coordinate flags (1/A=Å, 0/B=Bohr); extracts CHARGE and MS/MULT keywords from keyword line; falls back to CARTESIAN COORDINATES output blocks for MOPAC output files
- **Save as PDB** — Export molecular structure in PDB format with ATOM records and CONECT records for bond connectivity
- **Save as MOPAC Input** — Export in MOPAC input format with PM7 method, CHARGE and MS keywords, and internal coordinates in Å

## [0.7.1] - 2026-04-15

### Changed

- **Completely rewritten bond order detection algorithm** — Replaced simple covalent radii ratio-based estimation with a sophisticated multi-stage approach:
  1. **Element-pair-specific bond length specifications** — BOND_SPECS table defines expected bond lengths and tolerances for 11 common atom pairs (C-C, C-N, C-O, N-N, N-O, O-O, C-S, C-F, C-H, N-H, O-H) with multiple bond orders per pair
  2. **Pair-specific distance cutoffs** — BOND_CUTOFF dictionary provides maximum bond distances for 11 common atom pairs (e.g., CC: 1.9Å, CO: 1.7Å, CH: 1.3Å); other pairs use covalent radii sum + 0.5Å
  3. **Best-match fallback** — For BOND_SPECS pairs where no spec matches within tolerance, the closest spec is chosen by distance
  4. **Post-processing fixes** — C-O bonds forced to 1.0 or 2.0; Br bonds forced to single; N with 3 neighbors → all single; N with 2 neighbors → smart assignment (1.5+1.5, 2+1, or 1+2)
  5. **Valence-based refinement** — Iterative algorithm (max 10 iterations) that reduces bond orders when atoms exceed their maximum valence (H:1, C:4, N:3, O:2, F:1, S:6, P:5, Cl:1, Br:1, I:1, B:3); prefers reducing the bond whose ideal length is closest to the actual distance
  6. **Fallback ratio** — For atom pairs not in BOND_SPECS, uses covalent radii ratio (ratio < 0.85 → triple, < 0.90 → double, else single)

## [0.6.3] - 2026-04-15

### Added

- **Save as Turbomole .coord** — Export molecular structure in Turbomole coordinate format with automatic Å→Bohr conversion; writes `$coord`, `$end`, `$chrg`, and `$spin` sections
- **Save as ORCA Input (.inp)** — Export in ORCA input format with `* xyz CHARGE MULT ... *` block; uses charge and multiplicity from the original file
- **Save as MOL2 (.mol2)** — Export in Tripos MOL2 format with `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, and `@<TRIPOS>BOND` sections; bond orders mapped correctly (1→1, 1.5→ar, 2→2, 3→3)
- **Save as MDL Mol (.mol)** — Export in MDL Mol V2000 format with atom block (coordinates in Å×10) and bond block; aromatic bonds represented as type 4

## [0.6.2] - 2026-04-15

### Added

- **Turbomole .coord file support** — New parser for Turbomole coordinate files (`.coord`); reads `$coord` section with automatic Bohr → Ångström conversion; reads `$chrg` for charge and `$spin`/`$mult` for spin multiplicity
- **ORCA input file support** — New parser for ORCA input files (`.inp`); reads `* xyz CHARGE MULT ... *` coordinate blocks and `%coords` blocks; supports both inline xyz and xyzfile reference formats; extracts charge and multiplicity
- **ORCA output file support** — New parser for ORCA output files (`.out`); reads `CARTESIAN COORDINATES (ANGSTROEM)` blocks; supports multi-frame optimization trajectory with frame navigation; extracts charge and multiplicity from `Total Charge` and `Multiplicity` lines; falls back to `CARTESIAN COORDINATES (A.U.)` with Bohr → Å conversion
- **Molecular info display** — Shows atom count, charge, electron count, and spin multiplicity in the top-left corner of the 3D view; electron count calculated from atomic numbers minus charge; info updates automatically when atoms are added/deleted
- **Charge and multiplicity extraction** — GJF, LOG, coord, ORCA input, and ORCA output parsers now extract and return charge and multiplicity information

### Changed

- **File extension handling** — `.out` files are now dispatched to ORCA output parser (previously treated as Gaussian LOG); `.log` files remain as Gaussian LOG; auto-detection also checks for ORCA-specific markers
- **`parseLogFile` signature** — Now accepts optional `fileName` parameter to distinguish between Gaussian LOG and ORCA output formats

## [0.6.1] - 2026-04-15

### Fixed

- **Rotation direction bug** — When rotating the molecule to its back side, the rotation direction would reverse (up became down, left became right). Root cause: Euler angle rotation with 'YXZ' order suffers from gimbal lock near ±90°. Fixed by replacing Euler-based rotation with quaternion-based rotation (`THREE.Quaternion`), which has no gimbal lock and always rotates consistently in screen-space directions
- **Rotation sensitivity** — Increased rotation sensitivity from 0.005 to 0.008 per pixel of mouse movement, making rotation feel more responsive

### Changed

- **GJF connect section output** — When saving as GJF, all atom lines in the connect section are now written, including atoms with no bonds (e.g. `2` alone on a line) and atoms whose bond info was already listed in other atoms' lines. Previously, lines with only an atom number (no bond pairs) were skipped

## [0.6.0] - 2026-04-15

### Added

- **Bond order 0 (Remove bond)** — Bond Length modal now includes a "None (0) - Remove bond" option; selecting it deletes the bond between the two atoms
- **Select Atoms feature** — New toolbar button "Select Atoms" opens an input dialog where you can enter atom indices (1-based), ranges (e.g. `3-10`), or element symbols (e.g. `C H`), separated by spaces or commas; selected atoms are highlighted in yellow
- **Stay in editing mode** — After completing an edit operation (OK/Cancel), the viewer now stays in the current editing mode instead of switching back to View mode, allowing continuous adjustments
- **GJF file structure preservation** — When saving as GJF, the original file's Link 0 lines (`%chk`, `%mem`, `%nproc`), route line, title, charge/multiplicity, and any content after the connect section (e.g. mixed basis set info) are preserved; only the coordinate and connect sections are updated

## [0.5.0] - 2026-04-15

### Added

- **MOL2 file support** — New parser for Tripos MOL2 format (`.mol2`), reads `@<TRIPOS>ATOM` and `@<TRIPOS>BOND` sections with bond order support (including aromatic `ar` → 1.5)
- **Gaussian LOG file support** — New parser for Gaussian output files (`.log`, `.out`), reads `Standard orientation:` and `Input orientation:` coordinate blocks
- **Optimization trajectory navigation** — When opening a LOG file with multiple structures, ◀ Prev / Next ▶ buttons appear in the toolbar to step through optimization frames; frame counter shows current step and label
- **Jump to frame** — Direct input field to jump to a specific frame number in optimization trajectory
- **Auto play** — Automatically cycle through all optimization frames with 500ms interval
- **Fixed atom notation in GJF** — Enhanced GJF parser to handle coordinates with fixed atom markers like `C  -1  -7.678  -1.467  1.374` or `C  -7.678  -1.467  1.374  -1` by extracting the last 3 numeric values as coordinates
- **Aromatic bond display** — Aromatic bonds (order 1.5) now render as one solid line + one dashed line, distinguishing them from double bonds

### Changed

- **GJF connect section parsing** — No longer searches for explicit `connect` label; instead directly checks content after the blank line following coordinates. Lines are validated as all-numeric with max atom number ≤ total atoms; single-number lines (atoms with no additional bonds) are allowed
- **Bond order estimation thresholds relaxed** — Tightened double/triple bond thresholds to reduce false positives: triple ≤ 0.78 (was 0.80), double ≤ 0.88 (was 0.90), single > 0.88

### Fixed

- **Bond order display bug** — Bond orders > 1 from file connect sections were not displayed correctly in 3D (all showed as single bonds). Root cause: GJF parser used `Math.round()` which could alter bond orders, and `createBond` used strict `===` comparison that failed for float values. Fixed by preserving original float bond orders and using range-based comparison (`ord < 1.25` for single, `ord < 1.75` for aromatic, `ord < 2.5` for double, `ord < 3.5` for triple)
- **GJF connect section not detected** — Lines with only an atom number (no bond pairs, e.g. `43`) were incorrectly rejected by the `tp.length < 2` check, causing the entire connect section to be skipped. Fixed by allowing single-number lines in validation

## [0.4.0] - 2026-04-15

### Added

- **Extended element selection** — Add Atom dialog now includes 70+ elements (periods 1-6, common transition metals and lanthanides) instead of just 10
- **Bond order selection when adding atom** — Choose single, aromatic (1.5), double, or triple bond when adding a new atom
- **Bond order editing** — Bond Length modal now shows current bond order and allows changing it (single ↔ aromatic ↔ double ↔ triple); changes are reflected immediately in 3D display
- **GJF connect section in saved files** — When saving as GJF format, bond connectivity information (atom indices + bond orders) is now correctly written in the connect section, keeping atom numbering and bond data consistent after add/delete/edit operations
- **Default save path** — Save As dialog now defaults to the directory of the currently opened file instead of an arbitrary location

### Fixed

- **Variable declaration order** — Moved `selectedAtoms` and other state variables before `rebuildScene()` call to fix `indexOf` undefined error that caused black screen

## \[0.3.0] - 2026-04-14

### Fixed

- **Black screen root cause (final fix): Three.js now inlined into HTML** — Instead of loading Three.js as an external file (which failed due to CDN unreachability and local file loading issues in webviews), the library is now read at runtime and embedded directly into the HTML. This eliminates ALL script loading failures regardless of network, CSP, or webview configuration
- **CSP simplified** — `script-src` no longer needs external domains or `webview.cspSource`; only `'nonce-xxx'` is required since Three.js is inline
- **Visible loading indicator** — Shows "Loading 3D Viewer..." text while initializing, so users know the webview is working
- **Comprehensive error handling** — Added `try/catch` around entire initialization; any JavaScript error now displays a visible red error message instead of a silent black screen
- **Loading state management** — Loading indicator is hidden once Three.js renderer is created; error messages properly hide the loading indicator

## \[0.2.2] - 2026-04-14

### Fixed

- **Black screen root cause: Three.js CDN unreachable** — Bundled Three.js r128 locally in `media/three.min.js` instead of loading from cdnjs.cloudflare.com, which is often blocked or slow in China and corporate networks
- **Local resource loading** — Changed script loading to use `webview.asWebviewUri()` for reliable local file access, compatible with both local and Remote-SSH scenarios
- **CSP updated** — Replaced `https://cdnjs.cloudflare.com` in Content Security Policy with `webview.cspSource` for proper local resource authorization

## \[0.2.1] - 2026-04-14

### Fixed

- **Black screen on startup** — Changed layout from absolute positioning to CSS flexbox for the toolbar, status bar, and 3D container, ensuring the container always has correct dimensions
- **Container dimension fallback** — Added fallback to `window.innerWidth/innerHeight` when `container.clientWidth/clientHeight` is 0, preventing Three.js renderer from creating a 0×0 canvas
- **Three.js load check** — Added check for `THREE` undefined with visible error message when CDN fails to load
- **`acquireVsCodeApi()`** **multiple call bug** — Moved API acquisition to top-level (once per session) instead of inside `doSave()`, preventing error on second Save As click
- **CSP inline style violation** — Replaced inline `style="background:#c33"` on Delete button with CSS class `.mbtn-danger`, complying with Content Security Policy
- **Error display element** — Added `#error-msg` element for showing runtime errors to users instead of silent black screen

## \[0.2.0] - 2026-04-14

### Added

- **Molecular Editing Toolbar** — Mode-based toolbar with View, Bond Length, Bond Angle, Dihedral, Add Atom, Delete Atom, Save As, and Reset View buttons
- **Bond Length Adjustment** — Click 2 atoms to select, modal shows current bond length, choose which atom to fix/move, adjust via numeric input or slider with real-time 3D preview
- **Bond Angle Adjustment** — Click 3 atoms (2nd is central), modal shows current angle, fix/move either side, real-time 3D preview
- **Dihedral Angle Adjustment** — Click 4 atoms, modal shows current dihedral, fix/move either side, real-time 3D preview
- **Add Atom** — Click anchor atom, choose element (H/C/N/O/F/P/S/Cl/Br/I), set bond length, auto-calculated direction based on existing bonds
- **Delete Atom** — Click atom, confirm deletion; automatically re-indexes atoms and bonds
- **Save As** — Export modified structure as XYZ or Gaussian GJF format via VS Code save dialog (original file is never modified)
- **Atom Selection Highlighting** — Selected atoms glow yellow (emissive) for clear visual feedback
- **Status Bar** — Shows current editing mode and selected atoms
- **Cancel/Undo** — Cancel button in edit modals restores original coordinates
- **BFS Fragment Detection** — `getMovable()` uses breadth-first search to correctly identify which atoms move when adjusting geometry
- **Rodrigues Rotation** — `rotAroundAxis()` implements Rodrigues' rotation formula for accurate rotation around arbitrary axes
- **Real-time 3D Preview** — Slider and input changes immediately update the 3D molecular view

## \[0.1.0] - 2026-04-14

### Added

- 3D ball-and-stick molecular rendering using Three.js
- Gaussian `.gjf` / `.gjf03` / `.gjf09` / `.gjf16` / `.com` file parser
  - Reads Link 0 commands, route section, title, charge/multiplicity
  - Supports atomic numbers and element symbols (e.g., `6` → `C`)
  - Reads `connect` section with bond orders
- XYZ format parser with atomic number support
- Automatic bond detection using covalent radii (118 elements) + 0.45 Å tolerance
- Dual-colored bonds (half atom1 color, half atom2 color)
- Visual bond order distinction (single, double, triple)
- Interactive mouse controls:
  - Left drag: Rotate around molecule center
  - Scroll: Zoom
  - Middle/Right drag: Pan
- Atom hover tooltip showing element name and coordinates
- Reset View button
- Touch support (single-finger rotate, pinch-to-zoom)
- Custom editor integration (open .gjf/.xyz files directly)
- Explorer context menu integration
- Command palette integration
- Remote-SSH compatibility
- CPK atom coloring scheme
