# Molecular Viewer 3D

## English

**English** | [简体中文](#简体中文)

A VS Code / Trae extension for visualizing and editing molecular structures in 3D, designed for computational chemists working with Gaussian, ORCA, and other quantum chemistry software — especially on remote servers where GUI tools like GaussView are unavailable.

### Features

- **3D Rendering in Three Display Styles** — Atoms are spheres with CPK coloring and bonds are dual-colored cylinders, drawn in the model you pick: **Ball & Stick** (the default — a sphere of half the covalent radius per atom), **Space-filling** (each atom at its van der Waals radius, with no sticks) or **Stick** (the sticks stay and every sphere shrinks to the thickness of the sticks meeting at it, so the model reads as tubes without the beads). The **Display** selector is the second section of the View card at the bottom-right. What it sets is a property of the view, not of a tool: leaving View for Bond Length, Move Atoms or MM, switching frame or molecule, undoing an edit or running an optimisation all keep the style you chose, and the View card is only where you change it. Gaps are filled where two sticks meet, so a bend in a stick model shows no notch; a ghost atom stays the small marker it is in every style, and an atom always stays clickable — the stick model's joint is only a pixel or two across on a large structure, so a click that misses every sphere is resolved in screen space rather than being lost
- **Bond Order Support** — Visual distinction for single (1 line), aromatic (1 solid + 1 dashed), double (2 lines), and triple (3 lines) bonds. An aromatic bond's dashed line is planned per connected aromatic system rather than per bond, so consecutive aromatic bonds keep their dashes on one side: a ring shows them all inside (as it is drawn on paper) and a chain or each arm of a branch shows one continuous band, whichever way the ring happens to be stored in the file or turned in space
- **Auto Bond Detection** — When files lack explicit connectivity, bonds are derived from the coordinates by a single pass that applies two checks, each covering what the other cannot. The **geometric** check reads each length against element-pair-specific bond length specifications (C-C, C-N, C-O, etc.) with tolerances and pair-specific distance cutoffs, placing every bond on the order the measured length actually supports (candidate orders are limited to a spec tolerance or the pair's closest fit, so a long non-bonded contact is never promoted), pins terminal elements (H, halogens, noble gases) to single bonds, and falls back to a covalent-radii ratio for pairs with no entry. A **valence-conservation** check is then applied to any connected component that still holds an atom off its allowed valence: every atom is given a target valence from its element's permitted set (N 3 or 4, O 2/1/3, S 6/3/2/1), unsaturated atoms are paired by a maximum-cardinality matching until each atom reaches its target, and an assignment is accepted only when the formal charges it implies sum to the total charge the file states — so a nitro group comes out as one N=O plus one N–O and an N⁺ instead of a five-valent nitrogen, and an ammonium or a carboxylate is told apart from its neutral form. Geometric candidates are tightened first to the band an atom's degree and its neighbours' spare capacity allow; aromatic rings the second solver re-derived as alternating single/double bonds are restored to order 1.5 when that leaves every atom's valence unchanged; and a structure whose geometry is already self-consistent is left exactly as the geometric pass produced it. Every re-detection runs this one module — loading a file, both sides of a Diff comparison, and the Move Atoms release (which re-derives only the connected components containing the moved atoms, so a move lands on the same answer a fresh load of the moved coordinates would give) — and the webview receives it from the extension host rather than keeping a second copy of it. **The default is the geometric pass, unchanged** — the valence solver is opt-in, because it re-derives a whole connected component and would visibly re-assign bonds that were never in doubt. What the default *does* add cannot change a bond: **formal charges** for every atom, derived from the finished orders, plus a report of the atoms whose valence is not one their element may take and a count of how many components a solver could settle. Whichever mode runs, the result is held to a measured monotonicity property against the detector it replaces — no structure may end up with more atoms off their allowed valence, and every bond whose length fit gets worse must sit in a component the old result already had a valence problem in — backed by a runtime safety valve that rolls back a solver answer which fails to lower the component's total valence gap
- **GJF Connect Section** — Reads explicit bond information from coordinate section in GJF files, including bond orders (1.0, 1.5, 2.0, 3.0)
- **Molecular Info Display** — Shows the Hill-system molecular formula, atom count, charge, electron count, and spin multiplicity in the top-left corner of the 3D view; for crystal files, an additional Unit Cell composition line shows the full unit-cell formula. When a parser had to work around something in the file — a CIF with no readable symmetry operations, an incomplete unit cell, two sites folded together — the panel ends with one amber line naming the first of them and how many others there are, so a structure that is not what the file literally says is never shown without saying so. A structure that read cleanly shows nothing extra, and the panel reads exactly as before
- **Two-row toolbar** — The tool buttons are laid out over two rows so that the ones at the end are reachable without dragging a scrollbar: the editing modes and tools (View, Bond Length, Bond Angle, Dihedral, Rotate Group, Move Atoms, MM Optimize, Add/Delete/Replace Atom, Add H, Import Structure) on the first row, the selection pair and the structure and view actions plus the trajectory and molecule navigation (Select Atoms, Box Select, Diff, Undo, Save As, Reset View, ◀ n/N ▶, molecule selector) on the second. Each row scrolls horizontally on its own, and only when its own content does not fit — so a narrow panel still reaches every button, and a loaded trajectory never pushes the actions out of sight again
- **Interactive Mouse Control**:
  - Left drag → Rotate around molecule center
  - Scroll → Zoom in/out
  - Middle/Right drag → Pan (cursor-locked: the grabbed point stays exactly under the cursor, 1:1 cursor speed)
  - **Alt** + left drag / right (or middle) drag → instead of only turning or panning the view, write that rotation/translation **into the atom coordinates** (see below)
  - Hover atom → Show element name + coordinates
  - Arrow keys (↑/↓/←/→) → Rotate molecule (left/right around Y axis, up/down around X axis; ignored when typing in input fields; with **Auto Rotate** on, holding a key spins the view up and releasing it keeps the speed — see below)
  - Delete key → If atoms are selected (via Select Atoms / Box Select), opens batch delete confirmation; otherwise enters Delete Atom mode
- **Auto Rotate** — A small panel at the bottom-right of the view appears while the **View** mode is active and carries the **Auto Rotate** toggle. With it on, the view can be *thrown* from either input, and what it is given is angular momentum, so it keeps turning about a fixed world axis instead of stopping when the input does:
  - **A left-drag that rotates the molecule, released**: the throw's direction becomes the spin axis and its speed the spin rate (distance over time) — drag farther or flick faster to spin faster. Only the last 120 ms of the drag count, so a drag brought to a standstill before the button comes up does not spin.
  - **An arrow key, held**: it spins the view up like an accelerator rather than turning it at a fixed rate — **the longer the hold, the faster** — and releasing the key hands over the speed it reached. The axes are the same ones a drag uses (←/→ on Y, ↑/↓ on X), and two keys together give the diagonal.
  
  The spin rate is capped and has a floor, so a wild flick cannot blur the view and a near-zero throw reads as a stop. **The next left press stops it; a right/middle press pans and deliberately leaves it running.** Reset View, leaving Diff mode, leaving View mode and switching the toggle off also end it, and the first arrow key of a hold takes over from any spin already running. View-only: it rotates the camera, never the coordinates, so no edit or playback is affected, and in Diff mode the half you spun or keyed is the half that spins
- **Alt-drag writes the drag into the coordinates** — In **View** mode, holding **Alt** while you press turns the drag into an edit: **left** rotates the structure, **right/middle** translates it, and it is the atoms that move rather than the camera — **Save As** then writes the geometry you dragged it to, in any format. Orient a structure for a figure, or bring several of them into a common frame, without leaving the viewer. The picture is exactly the one a plain drag would have shown: the rotation is written about the structure's own centre (its centroid, which a rotation about it leaves fixed) and the translation along the same screen direction the pan uses, both expressed in the structure's own frame, so nothing is turned twice and nothing jumps on release. One undo step per drag — an Alt-click that moves nothing leaves the stack alone — and no bond is re-derived, because a rigid transform changes no distance. **Reset View** afterwards returns the camera to the file's axes: the orientation the drag produced now lives in the coordinates, not in the view. Only the press decides, so releasing Alt mid-drag keeps the drag an edit. Where the coordinates are tied to something else — **Diff** mode (two structures, two view states) and **crystal** structures (the cell) — the drag stays view-only and the panel says so; a running vibration or frame auto-play is ended first, and Auto Rotate is not launched from such a drag. It is stated in its own **Edit Coordinates** card, which sits directly above the **Auto Rotate** card — both are separate boxes in the panel layout, so they stack instead of overlapping each other or the Diff panel below them
- **Touch Support** — Single-finger rotate, pinch-to-zoom
- **Orientation Axes Indicator** — A small axes gizmo in the bottom-left corner is shown for every structure: molecular files display the Cartesian X/Y/Z axes (X red, Y green, Z blue) of the initial coordinate frame, tracking every view rotation so you always know which direction the molecule faces; crystal files display the lattice a/b/c axes. Reset View returns the indicator to the initial orientation
- **Remote-SSH Compatible** — Works seamlessly when editing files on remote Linux servers via VS Code/Trae Remote-SSH
- **VMD TCL Script Visualization** — Opening a VMD visualization script reproduces its display recipe: representations (CPK/VDW/Licorice/Lines/Surf) with their relative sphere scales, bond radii and probe radii; colour modes `Beta` (per-atom scalar on the R→G→B gradient, with the script's `mol scaleminmax` range or the data's own), `ColorID` (the reference renderer's 33-colour palette at measured values) and `Name`; `mol selection` expressions resolved against the loaded atoms; molecules from repeated `mol new` selectable one at a time; and `Surf` patches extracted as real solvent-accessible surfaces. Visualizations that ship as a bare PDB instead of a script are coloured by their own per-atom scalar automatically. See [VMD TCL Script Visualization](#vmd-tcl-script-visualization)
- **Performance** — Dirty-flag rendering: the scene only re-renders on interaction or structural change (not every frame), so idle GPU usage is eliminated; meshes are properly disposed (geometry, material, texture) on rebuild to prevent GPU memory leaks that previously caused progressive lag within a session; adaptive level-of-detail reduces sphere segment count for large structures (>200 atoms: 20×16, >1000 atoms: 12×10) to maintain interactive framerate; the axes indicator SVG is rebuilt at most once per animation frame (~60 Hz) instead of on every mousemove event; automatic bond detection uses a 3 Å spatial hash grid (O(N·neighbors) instead of O(N²) pair tests — ~18× faster at 2000 atoms, ~130× on sparse systems) with byte-identical results, applied to initial detection and to the post-release bond re-detection in Move Atoms mode, whose re-judgement is further narrowed to the affected region (the moved atoms and their bonded neighbors) instead of the whole molecule; Move Atoms drags update bond meshes incrementally — only bonds touching the moving atoms are rebuilt each frame, and the release-time refresh likewise rebuilds only the affected region's bond meshes (surviving meshes are kept and re-indexed into the new bond list), instead of rebuilding every bond mesh in the structure; Gaussian optimization logs are parsed in a single pass (geometry frames, convergence steps and frequency tables used to be three separate full-file scans), and ORCA .out / MOPAC .out files likewise parse in a single pass (their fallback geometry-block scans were merged into the main traversal)
- **Large Structures** — The scene's GPU objects are shared rather than owned per drawn thing, and neither the bonds nor the atoms are one object each any more. Every bond cylinder (a bond is two to eight cylinders — one per bond half, plus a dashed line of six for an aromatic bond) is drawn by a **single instanced mesh** with the radius and length in the instance matrix and the colour in the instance colour buffer. **The atom spheres are instanced too**: one shared unit sphere (tessellation chosen from the atom count) draws every atom of an appearance, scaled per instance by its drawn radius, with the element colour in the instance colour buffer; plain, selected and difference-marked atoms are **three instanced meshes of that appearance**, so a highlight change is two matrix writes and never adds an object or a draw call. A structure of N atoms therefore costs a fixed handful of three.js objects, geometries and materials rather than N of each — measured against the real WebGL context, CCDC 142619 at 2×2×2 (**2352 atoms**) renders in **14 draw calls**, where it used to be one per atom, and a selection repaint stays at 15. Scalar colouring is the only part that still grows with the data, because every distinct scalar is its own instance colour (the pool is a bounded cache). Picking became the analytic ray–sphere test the ray always meant — the nearest sphere along the ray in the group's own frame — and answers exactly what it answered before (every atom reachable in all three display styles, a click selecting the atom the hover resolved). A bond's transform is read off its two atoms by one routine that the scene build, a Move Atoms drag, a streamed MM chunk and a normal-mode playback all call, so a bond can never disagree with the atoms it spans. The viewer also handles the graphics context itself: a context lost mid-session is reported and the scene rebuilt when it comes back, the context is released when the page goes away, and a renderer that cannot be created explains that the window's GPU process is gone (and that reloading the window is what brings it back) instead of only repeating three.js

#### Molecular Editing

- **Bond Length Adjustment** — Select 2 atoms, view current bond length, choose which atom to fix, choose Translate Group (rigidly move atom + connected subtree, preserving relative positions) or Translate Atom (move only the selected atom), adjust via numeric input or slider with real-time 3D preview
- **Bond Angle Adjustment** — Select 3 atoms (2nd is the vertex), view current angle, fix/move either side, Translate Group/Atom option, real-time preview
- **Dihedral Angle Adjustment** — Select 4 atoms, view current dihedral, fix/move either side, Translate Group/Atom option, real-time preview
- **Bond Order Editing** — Change bond order (none / single / aromatic 1.5 / double / triple) in the Bond Length modal; selecting "None (0)" removes the bond; changes reflected immediately in 3D display
- **Add Atom** — Click anchor atom, choose element (70+ elements), set bond length and bond order, direction auto-calculated for polyhedral coordination with existing bonds (maximum angular separation from existing neighbors, naturally producing linear/trigonal/tetrahedral/octahedral geometry); repeated additions to the same anchor fan out instead of overlapping
- **Delete Atom** — Click atom and confirm; atoms and bonds are automatically re-indexed; or delete all currently highlighted atoms at once by selecting them first via Select Atoms / Box Select, then clicking Delete
- **Replace Atom** — Click the Replace Atom toolbar button, then click one or more atoms in the 3D view to select them (click again to deselect); click the Replace Atom button again to open the element picker (70+ elements), then confirm to replace all selected atoms with the chosen element while preserving their coordinates; works with crystal base atoms too
- **Ghost Atoms (counterpoise markers)** — An atom can be marked a **ghost**: a position with no nucleus, no electrons and no bonds — the marker a counterpoise (BSSE) calculation puts where the other fragment's atoms are, and the point a geometry is defined against. The **Add Atom** dialog gains a *Ghost atom* switch and an `X (dummy point)` element; **Replace Atom** gains a *Ghost* switch that turns the whole selection into ghosts **of their own elements** in one step, which is the counterpoise gesture (select one monomer, make it ghosts, save the Gaussian input) and the only reading that makes sense for a mixed selection. A ghost keeps the element it stands for — `C-Bq` in a Gaussian input, `C:` in an ORCA one — or is a bare dummy point (`X`); it is drawn as a small, faint, grey-blue marker, never as the atom it stands for, and it is labelled `C*3` / `X5` wherever the viewer names an atom. It is left out of everything that needs a nucleus: **bond perception never bonds it** (a ghost laid on top of a real carbon adds no sticks, and the real atoms' bonds are unchanged), the **formula and electron count** skip it (the info panel lists the ghosts separately instead), **Add H** and **Molecular Mechanics** skip it (a ghost is never moved and exerts no force, so a ghost parked 1 Å from an atom cannot drag it), and the **Bond Length** dialog refuses to create a bond to one. Turning a ghost back into a real atom re-derives that atom's bonds from the geometry — and only that atom's — in the same undo step. Ghosts are refused for crystal structures, whose formats (CIF, POSCAR, VESTA, cube) have no way to say "this site has no nucleus"
- **Ghost atoms survive every format** — Save As writes each format's own spelling of a ghost, and the parsers read all of them back, so a ghost survives a round trip: `C-Bq` / `Bq` (Gaussian input), `C:` / `X` (ORCA input), `X` (XYZ), `X` + `HETATM`/`DUM` (PDB), `Du` (MOL2), `*` (MDL Mol), `q` (`$coord`), `XX` (MOPAC). The formats that have no ghost concept are given the dummy spelling rather than the atom's — a ghost must never come back as a real atom — and the Save dialog says so before the file is written. A structure without ghosts saves byte-identically to before, and the saved provenance line records a ghosting (`2 atoms made ghosts`, `1 ghost atom added`, `1 ghost made real`)
- **Add Hydrogens (H Completion)** — The **Add H** toolbar button (next to Replace Atom) fills the selected atoms' missing valences with hydrogens. Each selected atom whose bond valence is below its element's empirical valence receives exactly the difference, at the X–H length and the H–X–Y angles the plugin's SMILES → 3D engine already uses — S completes to 2 and P to 3, a C–H bond is 1.09 Å, N–H 1.01, O–H 0.96, S–H 1.34, P–H 1.42, and the angles are the element's theoretical ones (C 109.47° sp3 / 120° sp2 / 180° sp; N 106.5°, O 104.5°, S 98°) — with the centre's hybridization read off the bond orders it already has (triple → sp, double or aromatic → sp2, amide-like N → sp2, else sp3), so a ring carbon's hydrogen lands in the ring plane at 120°, a methyl group comes out tetrahedral and an alkyne's hydrogen linear. Click atoms in the mode (or carry a selection in from Select Atoms / Box Select / Replace Atom) and click the button again to place them: **the atoms that were already there never move** — every hydrogen is placed around its centre, the whole addition is one Undo step, re-running it adds nothing, and an atom already at or above its empirical valence is reported instead of completed. A hydrogen the structure has no room for (a non-bonded contact that would come out shorter than a chemical bond) is left out and reported rather than allowed to overlap an atom. Not available for crystal structures or in Diff mode
- **Molecular Mechanics Optimization** — The **MM Optimize** toolbar button (immediately right of Move Atoms) opens a panel listing the built-in force field's parameters, **a selection control for every settable parameter** (the non-bonded / torsion / improper / electrostatics / hydrogen-bond / metal-coordination / polarisation term switches, the convergence force, step limit, time budget, step size and the non-bonded cutoff of 5–12 Å, default 8), **three precision presets that fill those controls in one click** — *Low* (quick pre-optimisation: bonds and angles only, every distance- and dihedral-dependent term off, 0.4 kcal/mol/Å, 120 steps, 1.5 s, 0.30 Å steps, 6 Å cutoff), *Medium* (the module's own defaults — the full geometric field: Lennard-Jones, the light Coulomb term and the hydrogen bond) and *High* (careful: everything including metal coordination and polarisation, 0.01 kcal/mol/Å, 2000 steps, 20 s, 0.05 Å steps, 10 Å cutoff) — with any hand-edited value turning the choice into "Custom", the atoms a run would move, and a **Start Optimization** button. Each control opens on the module's own default, the presets carry the module's own values, and every value is validated again host-side, so a choice can never put a run into an unreasonable state. Select atoms by clicking in the mode (toggle), or carry a selection in from Select Atoms / Box Select / Replace Atom / Add H / Rotate Group / Move Atoms — the selection survives the switch both ways and the panel follows it. **This is a restrained optimisation, not an isolated one**: only the selected atoms move, but every fixed atom is a full participant — its bonds, angles, torsions, impropers and non-bonded contacts all pull on the selection, and only a term whose atoms are all fixed is dropped (it cannot change). Pressing the button relaxes exactly those atoms and **every other atom comes back bit-identical**, so a fragment can be settled against the rest of a structure, or a whole distorted structure can be relaxed, without any risk to what was not selected. The panel then reports energy before → after (kcal/mol), the number of steps and why it stopped (converged / energy converged / step or time budget / no further gain), the final max and RMS force, the largest displacement, and how many of each term were live; the status bar carries the same summary. **A run is a session the host steps in chunks, so the panel watches it work**: the structure follows each chunk live (spheres and only the bonds that touch the moved atoms are redrawn), two canvases draw the energy and the max/RMS force against the step as the run goes, and a progress line names the step, the energy, its change and the forces. **One Pause / Resume toggle button** (the same single-button pattern the vibration panel uses) stops the run where it stands and continues it again — the session is kept, so **Resume** carries on from exactly that geometry, step size and conjugate direction, and the time spent paused does not count against the wall-clock budget; the label flips on the click itself, so the panel can never be left in a state it cannot leave. **Restore initial** writes the geometry the run started from back, as one Undo step. The geometry the run itself reaches is one Undo step too, taken in front of the first coordinate that actually changes — so a run that reports its starting geometry unchanged leaves the Undo stack exactly as it found it. Leaving the mode for anything else — another tool, a frame, a molecule, an Undo, an Import, a comparison — **interrupts** the run, so it can never keep writing coordinates into a structure the user has moved on from, and a late answer from a replaced run is ignored (every message carries the run's id). Closing the panel with its × only hides it: a small **MM Optimize** button next to the Convergence / Modes / Route buttons reopens it, and it doubles as the run's state while one is working ("running – step N" / "paused"). The whole run is **one Undo step**, and afterwards the bonds of the moved region are re-judged by the same perception module a fresh load of the moved coordinates would use (the Move Atoms release path), so a relaxation can only ever produce the geometry a reload would agree with. The force field is part of the plugin (`src/tools/mmOptimizer.ts`, zero dependencies): harmonic bond stretching whose reference lengths are the same table the SMILES → 3D engine and Add H use (`r₀` per element pair and bond order, aromatic included), harmonic angle bending towards the centre element's theoretical sp/sp²/sp³ angle (force-free beyond 30° of a *compressed* angle, so a strained ring is left as drawn, but up to 60° of a *widened* one, which is a drawing to correct), three-fold staggered torsions on sp³–sp³ bonds and two-fold planar ones where an sp² centre or a double/aromatic bond is involved (no term where hydrogen defines the four-atom set), an out-of-plane restraint at three-coordinate sp² centres, and a Lennard-Jones non-bonded term with σ from the van der Waals radii, 1–4 pairs at half strength, switched off smoothly between 6 and 8 Å, plus a 1-3 clash guard that re-admits a pair as pure repulsion if it comes inside 0.45 of its van der Waals sum. Beyond the bonded and Lennard-Jones terms the field carries the **light** additions the higher precision levels switch in: a **Coulomb** term whose partial charges come from the Pauling electronegativity difference (weighted by bond order, neutralised per structure, so no formal charges are needed) with a distance-dependent dielectric εᵣ = r and the same pair list, taper and 1-4 halving; a **hydrogen bond** that blends a 10-12 well over the H···A Lennard-Jones pair inside ≈1.6–2.8 Å and a D–H···A window (donors and acceptors N/O/F); a **metal-coordination** term that bends a metal's ligand–metal–ligand angles toward the ideal geometry of its coordination number (linear, trigonal, tetrahedral, square planar for the d⁸ metals, and cis/trans for five or six ligands) and leaves the metal–ligand distances to the bond term; and a **damped induced-dipole polarisation** (∝ 1/r⁴, calibrated so α = 1.45 Å³, q = 0.4 e at 2 Å gives ≈ −0.5 kcal/mol). The per-element tables cover the **main-group elements of periods 4 and 5** (Ga Ge Rb Sr In Sn Sb and their group analogues) alongside the first rows, added by one documented group-analogue rule, so no existing value changed and the heavier elements no longer fall through to the table default. The optimiser is a preconditioned conjugate gradient in which **every accepted step must lower the energy**, so a run is monotone: whatever it returns — converged, out of steps, or out of time — is never worse than what it was given. The run is bounded by 3000 atoms, 500 steps and 5 seconds. Not available for crystal structures (the force field has no image atoms) or in Diff mode (there are two structures to choose from). Verified with finite-difference gradient checks (worst relative error 1.3e-8), a restrained-optimisation proof (a lone mobile hydrogen lands at the table length from a fixed carbon, all three H–C–H angles to frozen hydrogens come out exactly 109.47°, and a terminal carbon feels the torsion of a frozen backbone), a 19-file corpus in which every structure improves, and a stub-driven run of the real client script. The streaming, preset and cutoff work is verified on top of that: the answers of the build *before* the session refactor were dumped for 16 structures × 8 option sets and every field (coordinates included) is **byte-identical** after it; a chunk-size sweep (1 / 3 / 7 / 1000 steps per chunk) and a pause/resume in the middle of each run reproduce the single call exactly; a whole streamed run driven through the compiled extension code — mmProgress → mmResult — is byte-identical to the module's single call, with Pause stopping the pump, Resume continuing it and Cancel (or a disposing panel) stopping it for good; and the new client code was driven against a stub DOM (start → live coordinates → pause → resume → result → restore, one undo snapshot at the first write, a stale run id that cannot write anything, the navigation cancel). The cutoff gradients still match central finite differences inside every taper window (5, 6, 8, 10, 12 Å; worst relative error 2.3e-7), the default cutoff reproduces the previous behaviour bit for bit, and the three presets were measured against each other (medium 51 steps to max|F| 0.027 kcal/mol/Å, high 81 steps to 0.0093, low 29 steps to 0.377; medium identical to the no-options answer). **A run is steerable while the rest of the viewer keeps working**, and the three ways the panel's Pause toggle, its curves and its state could disagree with what the host is actually doing are closed: a chunk already in flight when Pause is clicked is **recorded** (its steps join the curves, so a paused run's energy and force history is complete) while its coordinates are held back — so the geometry the user paused on is the geometry they keep looking at, until the host's own pause reply brings the snapshot it actually stopped at; the host's pump is **bound to the run it was scheduled for**, so a chunk chain left over from a run that was cancelled and replaced can no longer wake up and step the newer one; and a pause or resume for a run the host no longer owns is **answered** (as an ending) instead of ignored, so the button cannot sit in a state nothing will confirm. The curves are drawn into the canvas's own box and redrawn whenever the layout reflows, so a window resize, a panel move or a closed-then-reopened panel leaves them filling the panel instead of clipped at their right edge. Closing the panel with its × leaves the run going with the reopen button showing its state, changing a parameter or clicking new atoms mid-run rebuilds the panel around the running session without disturbing it, starting a run is refused while an **Add H** request is still in flight, and an Add H answer that lands later ends the run it would otherwise be appending atoms underneath. Nothing else that writes the structure can fight a run either: starting a run stops a **frame auto-play** first, and a **normal-mode playback** ends the run before it plays (the mirror of a run ending a playback when it starts) — both before the run captures the geometry **Restore initial** will put back, so the two can never own the coordinates at once. While a run is live the panel names **both** atom sets — the run's own, which its buttons act on, and the view's selection, which is the next run's — and its chunk is a fixed amount of *work* rather than a fixed pass count, so a Pause lands promptly on a large selection too.
- **Import Structure** — Click the Import Structure toolbar button (next to Add/Delete/Replace Atom); a dialog first asks for both the **number of copies** (typed into its input box, pre-filled with `1`) and the **source** — **From File** or **From SMILES** — so the count is settled before the structure. *From File* loads a second structure from any supported file (XYZ, GJF, MOL2, PDB, log/out trajectories, crystal formats, …) into the current view. Multi-frame files (optimization logs) open a frame picker; crystal files (CIF/VASP/Cube/VESTA) contribute their unit-cell atoms as a molecule. *From SMILES* opens an input box that validates the string **as you type** — showing the parser's own error message (with the offending position) while the input is not a valid SMILES — and builds the 3D structure when you press Enter: a locally generated conformer with explicit hydrogens, correct bond orders, aromatic ring perception, cis/trans double bonds and tetrahedral stereochemistry (`@`/`@@`), produced by the bundled zero-dependency SMILES → 3D engine (`src/tools/smilesTo3d.ts`, a faithful port of the reference `test/smiles_to_3d.py`), so no Python and no external cheminformatics toolkit is required. Both sources then behave identically: the imported structure is placed beside the existing one with a bounding-sphere separation test (both structures padded by 2.5 Å + a 2 Å gap), so it never overlaps; importing into an empty view centers the structure at the origin. The import can be repeated any number of times with the same or different sources — each copy is stacked further along +X at one uniform stride, and a multi-copy import is laid out and placed in a single pass, so the copy count does not multiply the placement work. Bonds merge smartly: if neither side has explicit bonds, connectivity is re-detected for the combined structure; if either side has authored bonds (MOL2/PDB/GJF connect, or any SMILES import), existing bonds are preserved and only the imported part's bonds are added (aromatic bonds carry order 1.5, so they render as 1 solid + 1 dashed line). Each import is one undo step, the camera zooms out to fit, and saving exports the combined structure
- **Select Atoms** — Click the Select Atoms toolbar button to show a floating panel (bottom-left corner) where you can type atom indices (1-based), ranges (e.g. `3-10`), or element symbols (e.g. `C H`) to highlight specific atoms in yellow; press Enter or click Select to apply, Clear to reset; you can also directly click atoms in the 3D view to select them; clicking a highlighted atom again deselects it (toggle); selections accumulate while in Select Atoms mode; the panel is a truly floating overlay that does not affect the molecular view size; the panel toggles with the toolbar button and hides automatically when switching to another mode; selections are preserved when switching between Select Atoms and Box Select
- **Box Select** — Independent toolbar button; in this mode, hold the left mouse button and drag a rectangle, then release to select all atoms inside the box (atoms are projected to screen coordinates and tested against the rectangle); selections accumulate with the same yellow highlighting as Select Atoms, and are preserved when switching between Box Select and Select Atoms; switching to any other mode clears the selection
- **Undo** — Toolbar button reverts the most recent edit/delete operation; snapshots are captured automatically when an edit modal opens (including Add Atom) and before delete/remove-disorder/batch-delete execute; the stack is capped at 50 entries and resets on file reload
- **Save As** — Export modified structure in 14 formats: XYZ, Gaussian GJF (preserving original file structure), CIF (fractional coordinates), VASP POSCAR (Direct fractional coordinates), Gaussian Cube (volumetric grid), VESTA (P1 symmetry, fractional coordinates), Turbomole Coord (Å→Bohr conversion), ORCA Input, MOL2 (with bond orders), MDL Mol, PDB (with CONECT records), MOPAC Input, and — for crystal structures — **XYZ / GJF − displayed cluster**, which write the region on screen as a molecule (the cage, pore or fragment a QM job starts from) instead of as a cell: the title line records the source file, the **cell parameters**, the **supercell bounds** and the atom count, so the file says which region of which structure it holds; a **Keep whole molecules** switch keeps the whole-molecule behaviour of the display (a molecule whose centroid is inside the bounds is written complete) or, unticked, cuts strictly at the bounds (atoms outside are dropped and bonds crossing a cell face are not written); GJF output preserves original Link 0, route, title card (written below the provenance line), charge/mult, and post-connect content; GJF atom coordinates use GaussView-style 8-decimal aligned format; connect section includes all atom lines; CIF/VASP/Cube/VESTA formats available only when crystal data is present
- **The saved title line records the structure's history** — Every format that owns a title or comment line (XYZ comment, GJF title card, MOL2 molecule name, MDL title line, PDB `TITLE` record, MOPAC title, ORCA comment, POSCAR/CUBE/VESTA headers) is written with a line of provenance, e.g. `MV: from opt.log, frame 12 of 40, Step 12; 2 atoms deleted, translated, displaced along mode 7 (38.2 cm-1) +0.45 A`: the file it came from, the frame of an optimisation trajectory (or the molecule of a VMD script), the edits made here — atoms added/deleted/replaced, hydrogens added, translated, rotated, bond length/angle/dihedral adjusted, MM-optimized, disorder removed, files merged in, supercell bounds — and any manual normal-mode displacement with its mode number, frequency and amplitude. The record is written as the edits happen and rides in the **Undo** snapshot, so undoing an edit takes its clause off the title with the geometry; a GJF keeps its own title card under the provenance line, and saving a file the viewer produced never stacks provenance on provenance. Formats with no title line (Turbomole `$coord`, CIF) are written exactly as before. The file's own title is appended after a `|` only when it describes the structure on screen — a trajectory's title belongs to the frame the file opens on, so saving frame 79 of a log does not end with frame 1's comment. The wording is short and ASCII, and the line is capped at 190 characters
- **Continuous Editing** — After completing an edit, the viewer stays in the current editing mode for repeated adjustments
- **Cancel/Undo** — Cancel button restores original coordinates before confirming edits
- **Rotate Group Around Axis** — A "Rotate Group" toolbar mode (grouped with View / Bond Length / Bond Angle / Dihedral) lets you select a group of atoms (by clicking in the 3D view or typing atom indices/element symbols, including the `ElementN` shorthand like `H581`) and define a rotation axis from 2+ axis atoms, then rotate the group 0–360° around that axis. The angle can be adjusted by dragging the smooth slider or by typing a value in the numeric input box and pressing Enter; both stay in sync. The axis is visualized as an orange line through the molecule. Axis atoms are validated chemically: any bonded pair among them must be a single bond (double/triple bonds are rejected, as they cannot rotate without breaking). The first and last axis atoms define the direction; intermediate atoms are included for validation only. Works for multi-atom collinear axes (e.g. `C32 C33 C43 C53`) or a simple two-atom bond axis (e.g. `C33 C43`). The axis atoms stay fixed while the rest of the selected group rotates around the line. Each angle adjustment is one undo step (Undo reverts a single adjustment and updates the slider display), while the panel's Reset button reverts the whole rotation session at once
- **Move Atoms** — A "Move Atoms" toolbar mode (next to Rotate Group) for translating and rotating a selected group of atoms. Select atoms by clicking in the mode (toggle), by typing in the Select Atoms panel, or by carrying over a selection from Select Atoms / Box Select / Replace Atom / Add H / Rotate Group. While in the mode, plain controls act on the whole view and Ctrl-modified controls act on the selected fragment: **left-drag / arrow keys** rotate the view; **right-drag** pans the view (cursor-locked); **Ctrl+left-drag / Ctrl+arrow keys** rotate the fragment around its own center, in the same drag direction as view rotation regardless of orientation; **Ctrl+right-drag** translates the fragment — the movement plane passes through the atom under the cursor, so that atom stays exactly under the cursor while dragging (when starting from empty space, the fragment-center plane is used and the fragment follows the cursor at 1:1 speed). Release Ctrl, the mouse button, or all arrow keys to finish. While dragging, existing bonds stretch with the atoms and only the bond meshes touching the moving atoms are rebuilt (synchronously every frame, so bonds never lag behind); connectivity and bond orders are re-detected once when the session ends (on release), re-judging only the bonds that touch the moved atoms plus their neighbors' bonds — the cost scales with the moved region rather than the molecule size, with results identical to a full re-detection (a safety check falls back to the full detector for degenerate pre-existing valence violations); the release refresh also rebuilds only the affected region's bond meshes instead of the whole scene (surviving meshes keep their GPU resources and are re-indexed into the new bond list, and a full rebuild runs only when the safety fallback trips). Files with explicit bond blocks such as MOL2/PDB (and GJF files with a connect section) keep their authored bonds and orders — only bonds whose endpoints are pulled beyond the bonding cutoff are dropped on release (deletion only: no re-ordering, and no new bonds are added when a fragment approaches). Each drag or key-hold is one undo step. Not available for crystal structures

#### Optimization Trajectory Navigation (LOG files)

- **Frame Stepping** — ◀ Prev / Next ▶ buttons to step through optimization frames
- **Jump to Frame** — Direct input field to jump to a specific frame number
- **Auto Play** — Plays the trajectory through all frames at 350 ms per step; playback stops on the last structure (the converged geometry) instead of wrapping back to the first one, and pressing Play there runs it again from the start
- **Route Keywords Panel** — Gaussian logs and ORCA outputs show a third collapsible panel button (⚙️ Route) alongside Convergence and Modes, displaying the parse result of the run's keywords: for every `--Link1--` job section, the raw keywords line plus the parsed keyword breakdown, and for ORCA the `!` simple-input line echoed into the output (recovered from the `INPUT FILE` block, joined when it wraps), one keyword per token, with the `opt` badge on runs that optimize. Gaussian: (print level, method/basis, each keyword with its option list, an `opt` badge on sections that request an optimization). Tokenization is paren- and quote-aware (`scrf=(SMD, Solvent=Acetone)` keeps its option group, `external="..."` quoted values stay intact, `iop(5/17=17)` parses as keyword + option), wrapped continuation lines are joined (`... opt` / `=(calcall) freq` folds the options into `opt`), and full-width characters are normalized. The same parse drives the Convergence panel visibility, so the two always agree. ORCA's simple input has no `key=value` option syntax (options live in `%block` directives), so its tokens are listed as keywords without options
- **Convergence Panel** — Gaussian, ORCA and xtb structure optimization runs display a collapsible panel (📈 Convergence) showing energy, force (max + RMS), and displacement (max + RMS) curves versus optimization step, mirroring GaussView's convergence plot; includes threshold lines (dashed amber) that are **the targets the file itself states**: an ORCA output carries them in its `Convergence Tolerances` block (`TolMAXG`/`TolRMSG`/`TolMAXD`/`TolRMSD`, with the table's own `Tolerance` column as the fallback when that block was not printed), so a run told to converge tighter shows the line it was actually held to; a Gaussian log states them only as the same table's Tolerance column, so its documented defaults stand in (force 0.000450/0.000300, displacement 0.001800/0.001200); a log from a program that states neither (xtb's `gnorm` criterion, a trajectory file) draws the curves with no borrowed threshold scale — plus a vertical resize handle to adjust panel height, and a close/reopen toggle. The panel appears only when the run actually optimizes — for Gaussian, single-point jobs with the `force` keyword (or `freq` jobs) print the same convergence table without optimizing, so the parser checks every `#...` route section for the `opt` keyword (`opt`, `opt=...`, full spelling `optimization`, `freq=calcall`; case-insensitive whole-token match) and shows no Convergence button when no route has it; ORCA prints its table only while it relaxes, so a frequency-only run shows none. Force and displacement charts use a logarithmic Y axis — axis ticks show the log₁₀ values with a rotated `log₁₀` axis label — so the convergence-threshold lines stay readable even when early steps dwarf the final values (energy keeps a linear axis since it may be negative). Energies are read from `SCF Done` lines, or from the optimizer's `Energy= ... NIter=` line in `opt=external` runs (e.g. Gaussian + xTB/XO combined jobs, which print no `SCF Done`), and from ORCA's `FINAL SINGLE POINT ENERGY` (each step's energy is the one printed for the geometry it relaxes — ORCA prints the geometry first, so the parser binds the energy back to the frame waiting for it). For xtb trajectories, the gradient norm (`gnorm`) from each step's comment line is plotted as the RMS Force curve. The Energy axis prints one shared integer base in the chart title and only the changing fraction on the ticks (e.g. `-6912 +` with ticks `.98 1.18 … 1.81`); the number of decimals adapts to the data so the five gridlines are always distinct and monotonic, even for a run whose span crosses an integer boundary or shrinks to a few milli-Hartree. Hovering a data point shows a tooltip with the step index, the point's value and (on the Force/Displacement charts) that step's energy, each printed with enough decimals to separate neighbouring steps
- **Normal Mode Vibration** — Gaussian logs with a frequency calculation (`freq` / `CalcAll`) and ORCA outputs with a `VIBRATIONAL FREQUENCIES` / `NORMAL MODES` pair display a collapsible Normal Modes panel (🎵 Modes) listing all harmonic frequencies (cm⁻¹) with symmetry labels; ORCA prints all 3N modes, so the 3N−6 translations and rotations it lists as `0.00 cm**-1` are dropped (a mode is only dropped when its frequency *and* its displacement vector are zero, so a genuinely soft mode survives), and its imaginary modes keep the negative sign ORCA prints them with, which the panel marks with `*` exactly as it marks a Gaussian one; clicking a frequency plays the vibration animation (60 FPS via `requestAnimationFrame`, amplitude normalized to 0.35 Å, 3-second cycle); imaginary frequencies (transition states) marked with `*`; the animation ends whenever something else takes over the structure — a frame switch, a molecule switch, Undo, Import, Diff, an Add H answer landing, an edit modal opening, or a Molecular Mechanics run starting — and starting one ends a frame auto-play rather than being stopped by its next tick (exactly one thing writes the coordinates at a time). A mode's row now carries **two separate selections**: the circle at the end of the row picks that mode for the displacement slider, while the row itself (or its ▶) plays the vibration — playing never re-aims the slider, and picking never starts a playback. With a mode picked, **Manual Displacement** writes it into the coordinates instead of animating it: drag the slider — a proportion of that mode from **−1.00 to +1.00**, where 1.00 moves that mode's most-moving atom by 1 Å and every other atom by its own share (weights and phases kept exactly as calculated), with the distance it comes to printed beside the fraction — and the structure is set to *equilibrium + Σ scaleₖ × modeₖ* on every move. **A mode is added to the distortion, not swapped in**: the block lists every component (`#2 +0.50 · #1 +0.30`), every contributing row is marked, dragging a component back to 0 takes that mode out, **Reset** clears all of them, and the largest displacement of the whole combination is shown next to the list — it is not the sum of the parts, because modes can partially cancel. Every write is rebuilt from the equilibrium snapshot, never accumulated, so dragging back and forth cannot drift and one slider value always means one geometry; merely picking a mode never rewrites anything. The session survives everything that does not replace or rewrite the structure — rotating the view, measuring, closing and reopening the panel, saving — and **Save As** writes the displaced geometry in any format (the dialog says which mode and how far). Playing any mode afterwards, that one included, swings from the **undisplaced** equilibrium and puts the displaced geometry back when it stops, so the animation always shows the mode the calculation produced and costs nothing of the adjustment; **Reset** returns exactly to the equilibrium and takes the session's undo step with it. Every coordinate edit closes the session first and keeps the geometry it produced (the edit starts from the displaced structure), while a frame switch, a molecule switch or an Undo drops it; a live Molecular Mechanics run refuses the slider (it owns those coordinates) and a Move Atoms session is committed before the slider takes over

#### Structure Diff

- **Diff Button** — Compare the current structure against another molecular file (any supported format, no need for matching formats)
- **Frame Selection** — If the comparison file is an optimization LOG/OUT, a QuickPick lets you choose which frame to compare
- **Skeleton Check** — Graph isomorphism via recursive backtracking atom matching (element + degree + neighbor signature with most-constrained-variable heuristic); uses the plugin's own bond detection algorithm (not the file's original bond orders) to ensure consistent skeleton comparison even when files specify different bond orders for the same connectivity; matching considers **connectivity only** (bonded or not), not bond order — a single bond in one molecule and a double bond in the other at the same position does not prevent skeleton matching
- **Side-by-Side View** — Left viewport shows the original molecule, right shows the comparison; each side has **independent** rotation, pan, and zoom controls. Left-drag rotates, scroll zooms and middle/right-drag pans **the half you grabbed** — panning is cursor-locked (the grabbed point stays exactly under the cursor) and is scoped to that viewport, so dragging the right structure never moves the left one and vice versa
- **Correct Aspect Ratio** — Each viewport uses its own camera aspect ratio matching the half-width, preventing horizontal distortion
- **Conformation Diff** — When skeletons match, highlights bond length, bond angle, and dihedral differences in orange on both sides; full list of all differences displayed (no truncation); each difference shows the percentage change relative to the average value
- **Grouped Coordinate Statistics** — When the skeletons match, a collapsible **Internal Coordinate Statistics** block sits above the difference list and summarises the whole comparison: for bond lengths (Å), bond angles (°) and dihedrals (°) it gives **n, rms, maxd (max |Δ|), mean and std** for every element-composition group (e.g. `C-C`, `C-C-H`, `C-C-C-H`) plus a combined **All** row. The statistics cover every compared coordinate, not just the ones above the threshold slider, so they stay put while you drag the slider (the slider only changes what is highlighted and itemised). Values are signed differences (left − right), a torsion's difference folded into [−180, 180] so a chain and its reverse cannot disagree; rms = √mean(Δ²), maxd = max|Δ|, mean = mean(Δ), std = population standard deviation. Group keys are canonical, so one coordinate always lands in the same group whatever the atom numbering is. Click the block header to collapse or expand it
- **Adjustable Threshold** — A percentage slider (0–20%) in the diff results panel lets you interactively adjust the sensitivity for what counts as a "difference"; all three types use the same relative difference formula (Δ/avg × 100%); for dihedrals near 0° where the average is very small, a fallback ensures robustness; the 3D highlights and difference list update in real-time as you drag the slider
- **Diff Info Panel** — Floating panel listing all differences sorted by magnitude, with atom labels and values from both structures; can be closed and reopened via "Show Results" button; diff results are always accessible and not overwritten by measurement info
- **Auto-Mirror Selection** — When selecting atoms for bond length/angle/dihedral measurement on one side, the corresponding atoms are automatically selected on the other side via the atom mapping, and both values plus the Δ difference are displayed simultaneously
- **Crystal Diff** — A crystal is compared as a **cell**, not as the two atom lists the screen and the file happen to hold: the displayed supercell (2352 atoms at 2×2×2) and the file's wrapped unit cell (294) are two drawings of the same object, so the comparison first checks the **six cell parameters** (a/b/c within 0.01 Å, α/β/γ within 0.05°) and reports **Unit Cells Differ** with both parameter sets when they really differ, then matches the **base cell's periodic graph** (the same skeleton matching, with every bond's image shift honoured) and reports **Skeletons Differ** with the base-cell counts when the graph differs. A file compared with itself is therefore **Structures are identical** at any supercell size. Bond lengths, angles and dihedrals are measured **through each bond's own periodic image**, so a bond across a cell face is its real length rather than a whole cell; on the left every displayed image of a differing site wears the orange mark, on the right the matching site does, and the auto-mirror selection still pairs the two sides. A crystal and a molecule are explicitly **not comparable** rather than silently "skeletons differ", and the comparison re-runs in place after an edit or a display-style change

#### Crystal Structure Viewing

- **Unit Cell Display** — When opening CIF, VASP, Cube, or VESTA files, the crystal unit cell is rendered with a wireframe outline showing the lattice boundaries; when the displayed region spans several cells the frame is drawn as a grid (interior cell-facet lines fainter than the outline), and the three edges meeting at the origin carry **a**/**b**/**c** labels in the axes colours
- **Supercell Boundary Controls** — 6 input boxes (a/b/c min/max) let you define a custom region of the crystal to display; supports fractional values (step 0.1) for partial cell selection; the view is a cutout of the crystal, not a deformation of the lattice; non-integer boundaries (e.g. -0.4 to 0.4) do not generate mirror images — atoms outside the range are simply not shown. **Molecules stay whole**: each finite molecule (one that does not extend periodically) is first unwrapped along the bond network into a single connected image, and copies are then placed by centroid — every image whose centroid falls inside the boundary draws the whole molecule, so a molecule straddling the boundary is shown completely (extending beyond it, as before) instead of being sliced by a cell face, while truly periodic networks (diamond, polymer chains) are still truncated at the boundary via split bonds. **One-click presets**: 2×2×2 / 3×3×3 / 4×4×4 buttons, plus **Reset** (back to 0–1) and **Align Origin** (keep the span, move the minimum to 0), sit next to the inputs. **The view does not jump**: a boundary change keeps the previous rotation centre, so the atoms you were looking at hold their place on screen, and the current atom selection is remapped onto the rebuilt structure by base-atom index and cell image — what was selected stays selected. **A size guard**: a boundary change that would exceed 20,000 displayed atoms is rejected with a message and the inputs restored, instead of freezing the editor
- **Atom Identity & Cell Info** — hovering an atom of a crystal shows, besides its Cartesian coordinates, its fractional coordinates, which periodic image it is (`cell (i,j,k)`), and — where present — its asymmetric-unit site label and its occupancy, so two images of the same base atom can be told apart; the info panel gains the unit-cell **Volume**, the molecules-per-cell count **Z** (counted from the placed finite molecules and isolated atoms; the line is omitted for purely periodic networks, where it would not be meaningful) and the crystal **Density** (occupancy-weighted, computed from atomic masses)
- **Local Crystal Modifications** — A geometry edit made on a displayed crystal atom (bond length, bond angle, dihedral, rotate group) is **written back into the unit cell**: the displacement of the image you edited becomes the fractional displacement of the base site it was drawn from, so the change survives a Save As to CIF/VESTA and appears in every periodic copy after the next rebuild, and one Undo step takes it back (the undo also restores the supercell bounds the edit was made under). A site that is on screen several times takes the edit **once**, not once per copy. Bonds drawn as stubs **across a cell face** are refused by the dialogs — the two atoms on screen are a whole cell apart on purpose, so the distance between them is not the bond length; grow the supercell and pick an image with both ends present. **Add Atom follows the same rule**: the new site is wrapped into the cell and its bond records which neighbouring image it points at, so a long bond is never drawn through the cell
- **Minimum Image Convention Bonds** — Crystal bonds are found by an exhaustive minimum image search, not by rounding the fractional difference: a repeat is only worth trying while the spacing between the planes it moves a point across is narrower than the cutoff, so an ordinary cell enumerates the classic 27 images and a narrow or skewed one widens to 5×5×5. Bond orders come from the **same perception module the molecular path runs** — the valence relaxation and the aromatic 1.5 included — so a molecule in a cell is drawn with the bonds it would have as an XYZ, and a molecule straddling a face is unwrapped and perceived whole instead of as two fragments. Bonds crossing cell boundaries are rendered as short split bonds extending outside the view (not as long bonds spanning across the view)
- **Rotating Axes Indicator** — For crystal structures the bottom-left axes gizmo shows the lattice a/b/c axes, rotating synchronously with the main 3D view (see the general Orientation Axes Indicator feature above, which covers molecular structures)
- **Fractional Occupancy Rendering** — Atoms with partial occupancy (disordered crystal sites) are rendered with partial coloration proportional to their occupancy value
- **Remove Disorder** — One-click button to clean up disordered sites: occupancy < 0.5 removed, > 0.5 set to 1, = 0.5 kept
- **Crystal Save** — Export modified crystal structures as CIF (fractional coordinates), VASP POSCAR (Direct fractional coordinates grouped by element), Gaussian Cube (volumetric grid), or VESTA (P1 symmetry, fractional coordinates). All four write the **unit cell** — the `baseAtoms` the session edits and every rebuild redraws from — and each file says so: the CIF declares `_symmetry_space_group_name_H-M 'P 1'` with the single identity operation and an `_audit_creation_method` line that names the body (the symmetry-expanded full cell, not the asymmetric unit), so reading the file back and saving it again keeps the atom count instead of expanding it a second time; its coordinates are wrapped into [0,1) and the cell carries `_cell_volume`, `_chemical_formula_sum` and `_atom_site_site_symmetry_multiplicity`. The POSCAR's comment line names the cell it lists, and its atoms are the cell's own (a 2×2×2 view still saves the unit cell, so no atom sits outside the declared box), the Cube is a real cube — a positive atom count (Bohr coordinates), the cell as the grid's box, and a grid of the size chosen in the dialog (default 40 points per axis, all-zero values) — instead of the negative count and 2×2×2 stub that readers took for orbital data, and the VESTA coordinates are wrapped the same way. The CIF and VESTA files carry the unit cell **as edited** (see Local Crystal Modifications), so a geometry changed on screen is the geometry in the file; the region on screen is exported as a molecule through the **cluster formats** of Save As

#### VMD TCL Script Visualization

Opening a `.tcl` script renders the structure it loads **the way the script asks for it**, instead of applying the viewer's own defaults:

- **The script is a display recipe** — Every `mol new` it performs becomes a molecule (with a ◀ n/N ▶ selector in the toolbar when there is more than one, as in the force-error series scripts that load one PDB per optimisation round), and each representation is reproduced: style and arguments (`mol representation CPK 0.97 0.10`, `Surf 1.0 1.0`, `VDW`, `Licorice`, `Lines`), colour mode (`mol color Beta` / `ColorID n` / `Name` / `Element`), selection, material and any `mol scaleminmax` colour-scale range. `mol modstyle` / `modcolor` / `modselect` / `modmaterial` modify an existing representation and `mol delrep` removes one — including the default representation `mol new` creates, so representation indices address what the script's author meant.
- **Selections are resolved against the loaded atoms** — `index` lists and `to` ranges (0-based; out-of-range indices are dropped), `element` / `name` / `type` symbol lists, `not` negation, and `all` / `none`. A selection that needs residue metadata the viewer does not keep (e.g. `resname BIB`) is widened to the whole molecule **and reported** on the status bar and in the info panel, never applied silently.
- **Sizes stay ball-and-stick, with the script's relative scale** — A script's representation arguments are relative to the reference renderer's own defaults, where `CPK 1.0` means van der Waals spheres (a space-filling model). Applying them as absolute radii would turn every structure into one solid blob and hide the bonds and the segmentation, so the viewer applies the script's **ratios** to its own geometry: sphere radius = the viewer's own radius × the script's sphere scale (an ordinary script comes out the size of an ordinary structure, `CPK 0.7` comes out smaller, `CPK 0.1` tiny), and bond radius is rescaled the same way against the reference default of 0.3 (`CPK … 0.3` is the usual stick thickness, `0.1` a thinner one, `Lines` a thin line). Anything physical is *not* scaled down: a `Surf` probe radius is an Å value, so the surface sits that far outside the van der Waals surface — which is what makes `show_eggs.tcl` read as coloured envelopes around the atoms. Atoms and bonds that no representation covers are not drawn at all, so a partial selection shows a partial view — unless the **Display** style is Space-filling or Stick, which are models of the whole structure rather than of one representation and therefore draw every atom (the script's colours and its surfaces still apply)
- **Colour modes** — `Beta` paints the per-atom scalar on the reference renderer's RGB ramp as it is actually rendered: red → green → blue with the midpoint at 0.5, but with the two passive channels on a 15/255 floor and the active channel ramping at twice the nominal rate, so the low end is `#FF0F0F`, the midpoint `#0FFF0F` and the high end `#0F0FFF` rather than the pure primaries its scale nominally names. The range comes from the script's `mol scaleminmax` when it writes one (a written range is treated as an unconditional statement that the column is a scalar channel, so a constant column still gets painted — `force_error_delta.3.pdb` is all zeros on a −10…+10 scale and renders on the middle of the ramp) and from the representation's actual data range otherwise (a script's hardcoded `0 – 100` scale really does squeeze 0.04–3.92 force errors into the red end — that is what the script does, and the info panel shows the range in use). `ColorID` uses the reference renderer's 33-entry palette at its measured values, including script overrides from `color change rgb`. `Name` / `Element` uses the reference renderer's element colours — carbon is the soft cyan `#40BFBF`, not grey.
- **`Surf` representations are real surfaces** — A `mol representation Surf <probe> <isovalue>` patch becomes the solvent-accessible isosurface of exactly the atoms it selects, extracted by marching tetrahedra at a 0.45 Å grid. Patches are drawn translucent: a script like `show_eggs.tcl` stacks fourteen of them, and opaque ones would hide everything but the outermost. They ride inside the molecule group, so they rotate and pan with the structure exactly as the atoms do. A fourteen-patch, 72-atom scene builds ~164 000 triangles in ~50 ms.
- **An egg list** — When the structure sits next to the toolchain's egg metadata (an `eggFile*` naming each egg and the unit numbers it wraps — `U12F1 : 1,10-15` — plus a `unitfile` with every unit's atoms), an **Eggs** button appears and a panel lists every egg the script draws: each row carries the egg's colour as a left accent bar, followed by its name and the units it contains (`units 1, 10–15`). All eggs are shown initially; unchecking a row hides exactly that egg's surface patch, checking brings it back, and a **Hide all / Show all** button toggles every egg in one click (completing a partial selection rather than resetting it). The choices survive scene rebuilds (edits, undo) and are reset when the patch list changes (switching molecules or frames). Metadata files are optional — without them the patches render unlabelled. When the directory holds **several** `eggFile*` files (leftovers from other settings or rounds), usually only one describes the structure: the viewer scores every candidate against the script's surface patches and asks which to use, each entry showing how many patches it matches (`14/14 patches match` vs `0/14 patches match`), so the right file is unmistakable; cancelling, or picking a file that matches nothing, simply leaves the patches unlabelled. The question is asked once per (directory, structure) pair, so a multi-molecule script does not repeat it.
- **Bare PDB visualizations are coloured by their own scalar** — Not every picture ships as a script: `A4` and the per-round files of `A5`/`A6` are PDB files whose 61-66 column *is* the visualization (a per-atom force error, a round-to-round difference, a unit index). A file is recognised as such when it carries nothing but coordinate records (no deposited-structure header) and its column is non-constant, and it then opens already painted on the same gradient over its own range, with the range on the status bar. A constant column (the exports written with 0.00) colours nothing, and a deposited PDB with a header record is left in its element colours. The **Color by Beta** toolbar button appears whenever such data is present, starts active for the automatic case and toggles it off again.
- **The viewer keeps its own background** — A script's `color Display Background white` is a setting of the reference renderer and is deliberately not applied here: a white canvas would clash with the editor theme and wash out white atoms (which is what `Name` colouring gives hydrogen). The directive is still parsed and shown in the scene model; `color scale method` and `color change rgb` overrides are applied.

A script also participates in **Structure Diff**: the comparison side is the molecule the script draws, coloured the way the script colours it.

All of it lives in one module, `src/parsers/vmdVisualization.ts` — the reference renderer's palettes and radii, the resolution of a parsed script into drawing properties, the script parser itself and the recognition of a bare scalar file — so the visualization-file concern has a single home (the parser index re-exports its API).

#### Supported File Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Gaussian Input | `.gjf`, `.gjf03`, `.gjf09`, `.gjf16`, `.com` | Reads Link 0, route, title, charge/mult, coordinates, connect section; supports fixed atom notation (`C -1 x y z` or `C x y z -1`) |
| XYZ | `.xyz` | Standard XYZ format with atom count header. A file holding several frames back to back (ORCA's `<job>_trj.xyz`, and the trajectory writers that copy its layout) is opened as a trajectory, with frame navigation and — when every comment line states an energy — the convergence chart; a frame that was cut off by a killed job is dropped rather than shown truncated |
| MOL2 | `.mol2` | Tripos MOL2 format; reads `@<TRIPOS>ATOM` and `@<TRIPOS>BOND` sections with bond order support (aromatic `ar` → 1.5) |
| Gaussian LOG | `.log` | Reads `Standard orientation:` / `Input orientation:` blocks; skips lines where atomic numbers are written as element symbols (e.g. from ONIOM external program calls); supports multi-frame optimization trajectory; extracts SCF energy and convergence criteria (force/displacement) for the convergence panel (shown only when a route-card section requests an optimization — `force`/`freq` single-point routes print the same tables but show no panel); parses `Harmonic frequencies` block for normal mode vibration playback. xtb trajectories (extended multi-frame XYZ with `xtb:` comment marker) are auto-detected and parsed as optimization trajectories with energy + gradient norm convergence |
| ORCA Input | `.inp` | Reads `* xyz CHARGE MULT ... *` blocks and `%coords` blocks; supports xyz and xyzfile coordinate formats; parses `units bohr/angs` for coordinate conversion |
| ORCA Output | `.out` | Reads `CARTESIAN COORDINATES (ANGSTROEM)` blocks, falling back to the `(A.U.)` blocks when a run suppressed the Ångström ones; supports the multi-frame optimization trajectory (each cycle plus the final energy evaluation, labelled `Cycle 1` … `Final`); extracts charge and multiplicity from the echoed input and the SCF settings; reads the `|Geometry convergence|` table, the run's own convergence tolerances, `VIBRATIONAL FREQUENCIES` + `NORMAL MODES` + `IR SPECTRUM`, and the `!` simple-input line for the Route panel. Format is decided by content, so an ORCA run saved as `.log` (or a Gaussian log saved as `.out`) still opens correctly |
| Turbomole Coord | `.coord` | Reads `$coord` section (Bohr → Å conversion), `$chrg` and `$spin`/`$mult` for charge and multiplicity; supports both `$chrg 1` and `$chrg=1` syntax |
| PDB | `.pdb`, `.ent` | Reads ATOM/HETATM records with fixed-column parsing; element from columns 77-78 or atom name; occupancy from columns 55-60 and the temperature-factor column (61-66) as a per-atom scalar, which is painted automatically when the file is a bare coordinate file with a non-constant column; CONECT records for explicit bonds (up to 8 connections per line); handles duplicate serials; ATOMIC_NUMBER_MAP extended to Cf (98) |
| MOPAC | `.mop`, `.mopac`, `.dat` | Reads MOPAC input format with `ELEM x 1 y 1 z 1` internal coordinates; supports atomic numbers or element symbols; auto-detects Å/Bohr units; extracts CHARGE and MS keywords; falls back to `CARTESIAN COORDINATES` output blocks |
| VMD TCL | `.tcl` | VMD visualization script, read as a display recipe rather than as a data file: every `mol new` (relative/absolute/quoted/braced paths, any supported structure format) becomes a selectable molecule, and each `mol representation` / `mol color` / `mol selection` / `mol material` / `mol addrep` becomes a representation — styles CPK, VDW, Licorice, Lines and Surf with their numeric arguments (sphere scale, bond radius, probe radius — applied as ratios to the viewer's own ball-and-stick geometry), colour modes `Beta`, `ColorID` and `Name`/`Element`, `mol scaleminmax` ranges, `mol modstyle`/`modcolor`/`modselect`/`modmaterial`, `color scale method` and `color change rgb`; a script's `color Display Background` is parsed but not applied (the viewer keeps its own background). Sidecar metadata (`eggFile*` + `unitfile` next to the structure) adds an egg list for toggling the surface patches |
| CIF | `.cif` | Crystallographic Information File; reads cell lengths/angles, the space group and the symmetry operations, and fractional/Cartesian atom sites; applies the operations to build the full unit cell; supports fractional occupancy with partial coloration. **All three dictionary generations are read** — `_symmetry_equiv_pos_as_xyz` (CIF 1.1), `_space_group_symop_operation_xyz` (CIF 1.2/core) and `_space_group_symop.operation_xyz` (CIF 2) are one tag in three spellings, the same for `_symmetry_space_group_name_H-M` / `_space_group_name_H-M_alt`, and a CIF 2 file's dotted `_atom_site.fract_x` is the same column as `_atom_site_fract_x` — so a modern file is no longer read as P1 with a quarter of its atoms. Operations may be a loop, a one-row item or a semicolon block, and a repeated operation is listed once. A file that states no operations, or states ones that cannot be parsed, is read as P1 **and says so** in the info panel instead of quietly shrinking the cell. Site elements are read from `type_symbol` with the oxidation state stripped (`Fe3+` → `Fe`, `O2-` → `O`, `Cu(II)` → `Cu`), `?`/`.` are treated as absent everywhere (a missing space group is not a group called "?" and is not written back on save), and an incomplete unit cell comes back as an empty structure with a warning rather than a 1 Å cube. Symmetry images are folded together only when they are the same point (within 0.1 Å), so two sites listed 0.3 Å apart stay two atoms; a near-coincident merge is reported with the site label and the distance. The file's own asymmetric unit is kept alongside the expanded cell |
| VASP POSCAR | `.vasp`, `POSCAR`, `CONTCAR` | VASP structure file; handles scale factor, lattice vectors, element types/counts, selective dynamics, and both Direct (fractional) and Cartesian coordinates; robust element-type detection (no false rejection of Cs/Cd/Ca/Cl etc.) |
| Gaussian Cube | `.cube` | Gaussian volumetric file; reads voxel vectors and atomic positions; auto-detects Bohr/Angstrom units; robustly handles corrupted voxel vectors by inferring an orthogonal lattice from atomic coordinate ranges |
| VESTA | `.vesta` | VESTA structure file; reads CELLP (cell parameters), STRUC (fractional coordinates + occupancy), GROUP (space group), and SYMOP (symmetry operations) sections; applies symmetry operations to expand the full unit cell; supports fractional occupancy for disordered sites |

### Installation

#### From VSIX (Recommended)

1. Download the latest `.vsix` file from [Releases](https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode/releases)
2. In VS Code / Trae, press `Ctrl+Shift+P`
3. Type `Extensions: Install from VSIX...`
4. Select the downloaded `.vsix` file
5. **For Remote-SSH**: Make sure to install the extension **on the remote server** (choose "Install on Remote" when prompted)

#### From Source

```bash
git clone https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode.git
cd Molecular-Viewer-on-VScode
npm install
npm run compile
npx vsce package --no-dependencies
# Then install the generated .vsix file
```

### Usage

#### Opening Molecular Files

1. **Right-click** a supported file in the Explorer → **Molecular Viewer: Open 3D Viewer**
2. **Command Palette** (`Ctrl+Shift+P`) → `Molecular Viewer: Open 3D Viewer`
3. **Custom Editor** — Double-click a supported file and select "Molecular 3D Viewer"

#### Importing a Structure (File or SMILES)

1. Click **Import Structure** in the toolbar (works in View mode; not available in Diff mode or for crystal structures)
2. One dialog asks for both things: type the **number of copies** into its input box (pre-filled with `1`; any whole number from 1 to 50 — a count that is not one re-opens the dialog with your input kept and the reason in its title), then pick where the structure comes from. The two entries below the box stay visible while you type:
   - **From File** — pick any supported file. Multi-frame optimization logs open a frame picker; crystal files contribute their unit-cell atoms as a molecule
   - **From SMILES** — type a SMILES string (e.g. `CCO`, `c1ccccc1`, `CC(=O)Oc1ccccc1C(=O)O`, `OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O`). The input box checks the string on every keystroke and keeps showing the parser's error message until it is valid; press **Enter** to generate the 3D structure
3. A structure large enough that the requested copies would add more than 20 000 atoms is trimmed to fit, with a warning naming the new count; cancelling any prompt (Esc) cancels the whole import
4. Every copy is added beside the existing structure without overlapping it (bounding-sphere separation) and copies are spaced at one uniform stride, so they stack along +X and never overlap each other; the whole import is a single **Undo** step and **Save As** exports everything together

SMILES coverage and limits — the generator reads the OpenSMILES organic subset plus the common extensions: branches, `%nn` ring closures, aromaticity, explicit bond symbols, dot-disconnected fragments, isotopes, charges, and `@`/`@@` tetrahedral configuration. Planar aromatic systems are emitted flat: ring substituents (hydrogens included) lie in the ring plane and fused rings — naphthalene, phenanthrene, the purine core of caffeine — are coplanar with each other; an aromatic ring fused to a saturated ring (indane, tetralin) likewise gets both of its benzylic bonds placed in the aromatic plane, as closely as the rigid saturated ring allows. Tetrahedral chirality is honoured; a double bond written without `/` or `\` defaults to **trans**; the conformation comes from local geometric rules rather than a conformer search, so a long flexible chain may fold tightly; isotopes are parsed but do not influence geometry; organometallic geometry is not special-cased.

#### Set as Default Viewer

Add to your `settings.json`:

```json
{
  "workbench.editorAssociations": {
    "*.gjf": "molecularViewer.editor",
    "*.gjf03": "molecularViewer.editor",
    "*.gjf09": "molecularViewer.editor",
    "*.gjf16": "molecularViewer.editor",
    "*.com": "molecularViewer.editor",
    "*.xyz": "molecularViewer.editor",
    "*.mol2": "molecularViewer.editor",
    "*.log": "molecularViewer.editor",
    "*.out": "molecularViewer.editor",
    "*.coord": "molecularViewer.editor",
    "*.inp": "molecularViewer.editor",
    "*.pdb": "molecularViewer.editor",
    "*.ent": "molecularViewer.editor",
    "*.mop": "molecularViewer.editor",
    "*.mopac": "molecularViewer.editor",
    "*.dat": "molecularViewer.editor",
    "*.tcl": "molecularViewer.editor",
    "*.cif": "molecularViewer.editor",
    "*.vasp": "molecularViewer.editor",
    "*.{POSCAR,CONTCAR,poscar,contcar}": "molecularViewer.editor",
    "*.cube": "molecularViewer.editor",
    "*.vesta": "molecularViewer.editor"
  }
}
```

The list mirrors the extensions the viewer registers; `.dat` is MOPAC's alternate input/output extension. Drop any line you do not want to divert from your normal editor.

#### Controls

| Action | Effect |
|--------|--------|
| Left mouse drag | Rotate molecule around its center |
| Mouse scroll | Zoom in / out |
| Middle / Right mouse drag | Pan view (cursor-locked; in Diff mode pans the structure on the grabbed half only) |
| Hover over atom | Show element + coordinates tooltip |
| Arrow keys | Rotate molecule (←/→ around Y, ↑/↓ around X); with Auto Rotate on, a held key spins the view up instead — longer hold, faster spin |
| `Esc` | Back out of one level, and only one: close the open dialog through its own **Cancel** (a previewed geometry is put back), abandon a box-select rectangle, leave the current tool for View, or clear the selection — a second press can only take back a smaller thing |
| `Ctrl+Z` (`Cmd+Z` on macOS) | Undo, same as the Undo button. Left alone while a dialog is open (it previews geometry on top of an undo step of its own) and inside a text field (that is the field's text undo) |
| `Ctrl+S` (`Cmd+S` on macOS) | Open the Save dialog, same as **Save As**; with that dialog already up, press it again to save the format already chosen |
| `Ctrl+A` (`Cmd+A` on macOS) | Select every atom on screen — in Diff mode, both halves. Nothing else in the view is selectable text, so the key can only mean this (a text field keeps its own select-all) |
| `Ctrl+C` (`Cmd+C` on macOS) | Copy the selected atoms and the bonds among them to the viewer's own clipboard |
| `Ctrl+X` (`Cmd+X` on macOS) | Cut: copy, then delete through the same path as Delete, as one undo step |
| `Ctrl+V` (`Cmd+V` on macOS) | Paste the clipboard beside the structure — the same 7 Å bounding-sphere spacing Import uses — and select what arrived, as one undo step. Not available for crystal structures or in Diff mode |
| ◀ ▶ (molecule selector) | Switch between the molecules a VMD TCL script creates (`mol new`), shown only when a script defines more than one |
| Color by Beta | Toggle scalar colouring for structures whose PDB columns 61-66 carry a scalar — each atom's sphere and the bond halves that meet it take that atom's scalar colour, so a force-error picture reads as one gradient rather than gradient spheres between element-coloured bonds; starts active when such a file was coloured automatically at load, and is hidden when a script is driving the view |
| Auto Rotate button | In the panel that appears in View mode (bottom-right): toggle flick-to-spin. With it on, a left-drag that rotates the molecule keeps turning after release, at the direction and speed it was released with (farther/faster drag → faster spin), and a held arrow key spins the view up — the longer the hold, the faster — handing over the reached speed on release; the next left press stops it and a right/middle press pans while leaving it running |
| Reset View button | Return to default view |

These shortcuts are handled by the viewer itself — nothing is added to the editor's keybindings, so no VS Code or Trae shortcut is taken over — and they work whenever the 3D view has focus. Keystrokes typed into a field still belong to the field, Esc excepted. The view is not a document: apart from the form fields, nothing on the page is selectable text (Ctrl+A therefore always means Select All of the structure, and a drag that ends over a panel leaves no selection behind) — only the load-error message stays selectable, so a failure can be copied. The fragment clipboard is the viewer's own, not the system one: copying never replaces text you copied in another application, and nothing is ever read back from it. A copy is mirrored into the panel's local storage on a best-effort basis, so a fragment copied in one panel can be pasted in another.

#### Editing Workflow

1. Click a toolbar button to enter an editing mode (e.g., **Bond Length**)
2. Click atoms in the 3D view to select them (selected atoms glow yellow)
3. A modal dialog appears showing the current value
4. Choose which atoms to **fix** vs. **move** using the dropdown
5. Adjust the value using the **numeric input** or **slider** — the 3D view updates in real-time
6. Click **OK** to confirm the change, or **Cancel** to revert
7. Use **Save As** to export the modified structure to a new file

### Architecture

```
┌─────────────────────┐          ┌──────────────────────────┐
│   Local (Windows)   │   SSH    │   Remote (Linux Server)  │
│                     │ ───────> │                          │
│  Trae IDE (UI)      │          │  Trae Server (Extension) │
│  ├─ Webview 3D      │ <─────── │  ├─ Parse .gjf/.xyz/.log │
│  ├─ Three.js (inline)│  data   │  ├─ Bond detection       │
│  ├─ Editing UI      │          │  ├─ Save file (VS Code)  │
│  └─ Mouse events    │          │  └─ Return molecule data  │
└─────────────────────┘          └──────────────────────────┘
```

The extension runs on the **remote side** (reading files, parsing, saving), while the Webview renders on the **local side** (Three.js inlined into HTML, mouse interaction, editing UI).

### Project Structure

```
molecular-viewer/
├── src/
│   ├── extension.ts           # Extension entry point
│   ├── types.ts               # TypeScript type definitions
│   ├── parsers/               # File-format readers: text in, molecule out
│   │   ├── index.ts           # Parser dispatcher (auto-detect format)
│   │   ├── gjfParser.ts       # Gaussian .gjf parser (connect section, fixed atoms)
│   │   ├── xyzParser.ts       # XYZ parser; multi-frame XYZ (trajectories)
│   │   ├── mol2Parser.ts      # Tripos MOL2 format parser
│   │   ├── logParser.ts       # Gaussian LOG parser (optimization trajectory)
│   │   ├── xtbParser.ts       # xtb LOG parser (optimization trajectory, content-sniffed)
│   │   ├── coordParser.ts     # Turbomole .coord parser (Bohr → Å)
│   │   ├── orcaInpParser.ts   # ORCA input .inp parser
│   │   ├── orcaOutParser.ts   # ORCA output .out parser (trajectory, convergence, frequencies)
│   │   ├── pdbParser.ts       # PDB format parser (ATOM/HETATM/CONECT)
│   │   ├── mopacParser.ts     # MOPAC input/output parser
│   │   ├── vmdVisualization.ts # Visualization files: VMD TCL scripts (molecules, representations, selections) + palettes/radii + representation → colour, size, surface patches + bare scalar-file recognition
│   │   ├── cifParser.ts       # CIF crystal structure parser (symmetry, occupancy)
│   │   ├── vaspParser.ts      # VASP POSCAR/CONTCAR parser (scale, Direct/Cartesian)
│   │   ├── cubeParser.ts      # Gaussian Cube parser (voxel vectors, robust lattice)
│   │   └── vestaParser.ts     # VESTA structure parser (CELLP, STRUC, SYMOP)
│   ├── tools/                 # Chemistry engines: molecule in, molecule out (no file I/O)
│   │   ├── bondPerception.ts  # Bond detection and order perception — pair tables, order lookup, valence solver, aromaticity; serialised into the webview as `BP` so the page runs the same module the host does
│   │   ├── bondDetector.ts    # Thin entry point over bondPerception (`ensureBonds`), the single call the host makes
│   │   ├── smilesTo3d.ts      # SMILES → 3D coordinates (zero-dependency, port of test/smiles_to_3d.py)
│   │   ├── addHydrogens.ts    # Add H: empirical-valence completion from the engine's tables, frozen skeleton
│   │   └── mmOptimizer.ts     # MM Optimize: built-in force field (bond/angle/torsion/improper/LJ + light electrostatics/H-bond/metal-coordination/polarisation), period-4/5 element tables, settable cutoff, low/medium/high presets, monotone optimizer as a steppable session (pause/resume/history)
│   └── webview/
│       └── molecularViewer.ts # Custom editor + Three.js webview + editing
├── dist/                      # Compiled JavaScript (pre-built)
├── test/                      # Sample molecular files
├── media/
│   └── three.min.js           # Three.js r128 (bundled locally)
├── package.json
├── tsconfig.json
└── LICENSE
```

### Development

```bash
# Install dependencies
npm install

# Compile TypeScript
npm run compile

# Watch for changes
npm run watch

# Type check without emitting
npm run lint

# Package as .vsix
npm run package
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full release history.

### License

[MIT](LICENSE)

---

## 简体中文

[English](#english) | **简体中文**

一个用于在 3D 中可视化与编辑分子结构的 VS Code / Trae 插件，面向使用 Gaussian、ORCA 等量子化学软件的科研人员——尤其是在无法使用 GaussView 之类图形界面的远程服务器上工作的场景。

### 功能特性

- **三种显示样式** — 原子以 CPK 配色的球体、化学键以双色圆柱绘制，并可选模型：**球棍**（默认，原子为半个共价半径的球）、**空间填充**（每个原子取范德华半径，不画键）与 **stick**（保留键，球的半径缩到该原子上所连键的粗细，整体读作无珠的管）。**Display** 选择器位于右下角 View 卡片的第二节。它设定的是**视图**的属性而不是某个工具的：切到 Bond Length、Move Atoms、MM，切换帧或分子、撤销一次编辑、跑一次优化，都不会退回球棍模型——View 卡片只是修改它的地方。两根键相接处由同样的粗细补上，因此 stick 模型的折角没有缺口；幽灵原子在任何样式下都保持它那个小标记；原子也始终点得中——大结构里 stick 的接头只有一两个像素宽，因此落空的点击会在屏幕空间里重新判定，而不是丢失
- **键级支持** — 单键（1 条线）、芳香键（1 实线 + 1 虚线）、双键（2 条线）、三键（3 条线）在显示上可区分。芳香键的虚线按「一个连通的芳香体系」统一规划而不是逐键决定，因此相邻的芳香键虚线落在同一侧：芳环的虚线全部朝环内（与纸面画法一致），直链与叉链的每条臂各是一条连续的带——无论文件里键的存储顺序、还是结构在空间中的朝向如何
- **自动键感知** — 文件不携带显式连接关系时，由坐标推导化学键：一趟遍历施加两道检查，各自覆盖对方做不到的部分。**几何**检查把每条键长对照元素对专属的键长规格（C–C、C–N、C–O 等）在容差与元素对专属距离截断之内比对，让每条键落在其实测键长真正支撑的键级上（候选键级被限制在规格容差之内、或取该元素对最贴合的值，因此一条较长的非键接触永远不会被提升为键），把端基元素（H、卤素、稀有气体）钉在单键上，对表中没有的元素对回退到共价半径比值。随后对**仍存在偏离允许价的原子**的连通分量施加**价态守恒**检查：每个原子按元素允许价集合（N 3 或 4、O 2/1/3、S 6/3/2/1）取得目标价，未饱和的原子以最大基数匹配两两配对直至各自达到目标价，只有当该指派所隐含的形式电荷之和等于文件声明的总电荷时才被接受——于是硝基得到 1 个 N=O 加 1 个 N–O 和 1 个 N⁺，而不是一个五价氮；铵离子或羧酸根也与其中性形式区分开来。几何候选先被收紧到「该原子的度数与其邻居剩余容量所允许的带」；被第二个求解器重推成交替单双键的芳香环，在「每个原子的价都不改变」的前提下恢复为 1.5 级；几何本身已经自洽的结构则完全保持几何检查所产出的结果。每一次重推都走这同一个模块——载入文件、Diff 比较的两侧、以及 Move Atoms 释放时（只重推包含被移动原子的连通分量，因此一次移动得到的结果与直接载入移动后坐标所得的结果完全一致）——并且由扩展宿主下发给 webview，而不是在页面里保留第二份实现。**默认即几何检查，行为不变**——价态求解器是选择性启用的，因为它会重推整个连通分量，会明显改写从未存疑的键。默认路径新增的内容不可能改变任何一条键：由最终键级推导出的**每个原子的形式电荷**，加上一份「价不在其元素允许集合内」的原子清单，以及一个「求解器可解决的份量」计数。无论运行哪种模式，结果都要通过一项可度量的单调性检验——对照它所替代的检测器：任何结构都不得出现更多偏离允许价的原子，而任何「键长拟合变差」的键都必须位于旧结果本就存在价态问题的分量之内——并由一个运行时安全阀兜底：求解器给出的答案若未能降低该分量的总价态缺口，就被回滚

- **GJF Connect 段** — 读取 GJF 文件坐标段中的显式键信息，含键级（1.0、1.5、2.0、3.0）
- **分子信息显示** — 在 3D 视图左上角显示 Hill 体系分子式、原子数、电荷、电子数与自旋多重度；晶体文件额外显示一行「晶胞组成」，给出整个晶胞的分子式。当解析器不得不绕开文件里的某个问题时——CIF 里读不到任何对称操作、晶胞不完整、两个位点被并成一个——面板末尾会多出一行琥珀色文字，给出其中第一条以及还有几条，因此**「屏幕上这个结构与文件字面不完全一致」这件事不会被悄悄吞掉**。解析器没遇到问题的结构一概不加，面板与以前逐字相同
- **两行工具栏** — 工具按钮排成两行，末尾的按钮不必再拖滚动条才点得到：编辑模式与工具（View、Bond Length、Bond Angle、Dihedral、Rotate Group、Move Atoms、MM Optimize、Add/Delete/Replace Atom、Add H、Import Structure）在第一行，选择工具以及结构/视图操作与轨迹、分子导航（Select Atoms、Box Select、Diff、Undo、Save As、Reset View、◀ n/N ▶、分子选择器）在第二行。两行各自独立横向滚动，且只在自身内容放不下时才出现滚动条——因此面板再窄也能点到每一个按钮，载入轨迹后也不会再把操作按钮挤到看不见的地方
- **交互式鼠标控制**：
  - 左键拖拽 → 绕分子中心旋转
  - 滚轮 → 缩放
  - 中键/右键拖拽 → 平移（光标锁定：抓取点始终精确停在光标下，光标 1:1 速度）
  - 按住 Alt + 左键 / 右键（中键）拖拽 → 不只是转/平移视图，而是把这次旋转/平移**写进原子坐标**（见下文）
  - 悬停原子 → 显示元素名 + 坐标
  - 方向键（↑/↓/←/→）→ 旋转分子（左右绕 Y 轴、上下绕 X 轴；在输入框中打字时忽略；开启 **Auto Rotate** 后，按住方向键变成加速自转、松开保留转速，见下文）
  - Delete 键 → 若已有选中原子（经「选择原子 / 框选」），弹出批量删除确认；否则进入「删除原子」模式
- **自动旋转（Auto Rotate）** — **View** 模式下，视图右下角会出现一个小面板，其中放着 **Auto Rotate** 开关（避开左下角坐标轴指示器所在的位置），工具栏不会多出一个按钮。开启后，两种输入都能把视角「甩」出去，交给它的是角动量，因此输入停下来后它会绕一条固定的世界轴继续自转：
  - **按住左键旋转后松开**：这次甩动的方向即自转轴、速度即自转速率（拖动距离 ÷ 时间）——拖得越远、甩得越快，转得越快；且只统计**松手前最后 120 ms**，所以拖到停住再松手不会转。
  - **按住方向键**：此时方向键像油门——视角被「加速转起来」，而不再是以固定速率转动，**按住越久转速越快**；松开该键即把已达成的转速交给自由旋转继续。轴向与拖拽完全一致（←/→ 绕 Y、↑/↓ 绕 X），同时按两个键则得到对角轴。

  自转速率有上限、也有下限（低于下限按停止处理），因此再野的一甩也不会把画面甩糊。**下一次按住左键即停止；按住右键/中键是平移，转动不会停。** Reset View、退出 Diff 模式、**离开 View 模式**、以及关掉该开关也都会结束自转；一次按住中的**第一个**方向键会接管已有自转（之后的按键与按住产生的重复 keydown 故意不接管，否则每次重复都会把加速重新从基数开始，永远加不上去）。它只操作视角（与拖拽写的是同一个旋转四元数），从不改动坐标，因此不会与任何编辑或回放冲突；Diff 模式下你甩/按键的是哪一侧就转哪一侧。Move Atoms 里 Ctrl+方向键旋转的是选中片段，与这里无关
- **Alt 拖拽把这次拖拽写进坐标** — 在 **View** 模式下按住 **Alt** 再拖拽，拖拽就变成一次编辑：**左键旋转**结构、**右键/中键平移**结构，动的是原子而不是相机，因此随后的**另存为**（任意格式）写出的就是拖到的那套几何——不用离开查看器就能把结构摆成做图要的朝向，或把几个结构摆进同一个坐标系。画面与普通拖拽**完全一致**：旋转绕结构自身的中心（质心，绕它转动时它不动）写入，平移沿与视图平移相同的屏幕方向写入，两者都换算到结构自身的坐标系，所以既不会旋转两次，也不会在松手时跳一下。每次拖拽是单步 Undo（按住 Alt 点一下而没有位移则不会在栈上留下任何东西），且**不重推任何键**——刚体变换不改变任何距离。此后的 **Reset View** 把相机复位到文件的坐标轴：这次拖出朝向已经写在坐标里，而不是写在视图里。只有按下那一刻的态度算数，中途松开 Alt 仍然继续编辑。坐标被别的东西绑住的地方——**Diff** 模式（两个结构、两套视图状态）与**晶体结构**（晶胞）——拖拽保持纯视图行为，面板里会说明；开始拖拽会先结束正在播放的振动或帧自动播放，Auto Rotate 也不会由这种拖拽甩出去。这项功能写在它自己的 **Edit Coordinates** 卡片里，就排在 **Auto Rotate** 卡片正上方——两张卡片在面板布局里各有一个碰撞箱，因此会上下堆叠，而不是互相压住或压住下面的 Diff 面板
- **触摸支持** — 单指旋转，双指捏合缩放
- **方向坐标轴指示器** — 左下角为每个结构显示一个小坐标轴：分子文件显示初始坐标系的笛卡尔 X/Y/Z 轴（X 红、Y 绿、Z 蓝），并跟随每一次视图旋转，因此始终知道分子朝向哪一边；晶体文件显示晶格 a/b/c 轴。Reset View 会把指示器复位到初始朝向
- **远程 SSH 兼容** — 通过 VS Code/Trae Remote-SSH 编辑远程 Linux 服务器上的文件时无缝工作
- **VMD TCL 脚本可视化** — 打开 VMD 可视化脚本时复现它的显示「配方」：表示方式（CPK/VDW/Licorice/Lines/Surf）及其相对球半径、键半径与探针半径；配色模式 `Beta`（逐原子标量走 R→G→B 渐变，取脚本 `mol scaleminmax` 指定的范围或数据自身范围）、`ColorID`（参考渲染器按实测值取用的 33 色板）与 `Name`；`mol selection` 表达式针对已载入原子求解；重复 `mol new` 产生的多分子可逐个切换；`Surf` 片被提取为真实的溶剂可及表面。以裸 PDB 而非脚本形式提供的可视化，会按其自带的逐原子标量自动上色。参见 [VMD TCL 脚本可视化](#vmd-tcl-脚本可视化)
- **性能** — 脏标记渲染：场景只在交互或结构变化时重绘（而非每帧），因此空闲时 GPU 占用为零；每次重建都正确释放网格的几何体、材质与贴图，杜绝此前会导致同一会话内逐渐卡顿的 GPU 显存泄漏；自适应细节层次会为大结构降低球体分段数（>200 原子：20×16，>1000 原子：12×10）以维持交互帧率；坐标轴指示器 SVG 每帧最多重建一次（约 60 Hz），而不是每次 mousemove 都重建；自动键感知使用 3 Å 空间哈希网格（O(N·邻居数) 取代 O(N²) 的两两比较——2000 原子约快 18×，稀疏体系约快 130×）且结果逐字节一致，用于初始检测以及 Move Atoms 模式释放后的键重推，而后者的重判范围又被进一步收窄到受影响区域（被移动原子及其成键邻居）而非整个分子；Move Atoms 拖拽以增量方式更新键网格——每帧只重建与移动原子相关的键，释放时的刷新同样只重建受影响区域的键网格（存活的网格被保留并重新索引进新的键表），而不是重建结构中的每一个键网格；Gaussian 优化日志改为单趟解析（几何帧、收敛步与频率表原本是三次整文件扫描），ORCA .out / MOPAC .out 同样单趟解析（其回退的几何块扫描并入主遍历）
- **大结构** — 场景的 GPU 对象是共享的，而不是每个被绘制对象各持一份，而且键与原子都不再是「一个原子一个对象」。每一段键圆柱（一条键由 2 到 8 段圆柱组成——两个半键各一段，芳香键另加一条六段的虚线）都由**同一个实例化网格**绘制，半径与长度写在实例矩阵里、颜色写在实例颜色缓冲里。**原子球同样实例化**：一个共享的单位球（细分程度按原子数选取）绘制同一外观的所有原子，每个实例按实际绘制半径缩放，元素颜色走实例颜色缓冲；普通、选中、差异标记三种状态各是一张实例化网格，因此改一次高亮只是两次矩阵写入，**不会新增对象或 draw call**。一个 N 原子的结构因此只占固定的一小批 three.js 对象、几何体与材质，而不是各 N 份——在真实 WebGL 上下文上实测，CCDC 142619 的 2×2×2（**2352 原子**）一帧**只有 14 个 draw call**（此前是每个原子一个），选中重绘仍是 15 个。标量配色是唯一仍会随数据增长的部分——每个不同标量就是自己的实例颜色（资源池是有上限的缓存）。拾取改成了射线本来就表达的解析射线–球求交（在结构自身坐标系里取沿射线最近的球），结果与以前逐点一致（三种显示样式下每个原子都可指向，点击选中的正是悬停解析出的那个原子）。键的变换由**同一个函数**从其两端原子算出，场景构建、Move Atoms 拖拽、MM 分块流式推进与简正模式播放都调用它，因此键绝不会与它跨越的两个原子不一致。查看器同时自行处理图形上下文：会话中途丢失的上下文会被上报，并在它恢复后重建场景；页面真正销毁时释放上下文；而完全无法创建渲染器时会说明窗口的 GPU 进程已失效（并指出重新加载窗口即可恢复），而不是只把 three.js 的原文抛给用户

#### 分子编辑

- **键长调节** — 选择 2 个原子，查看当前键长，选择固定哪一个原子，选择「平移基团」（刚性地移动该原子及其连通子树，保持相对位置）或「平移原子」（只移动所选原子），通过数值输入或滑块调节，3D 实时预览
- **键角调节** — 选择 3 个原子（第 2 个为顶点），查看当前键角，固定/移动任一侧，可选平移基团/原子，实时预览
- **二面角调节** — 选择 4 个原子，查看当前二面角，固定/移动任一侧，可选平移基团/原子，实时预览
- **键级编辑** — 在键长弹窗中修改键级（无 / 单键 / 芳香 1.5 / 双键 / 三键）；选择「None (0)」即删除该键；改动立即反映到 3D 显示
- **添加原子** — 点击锚点原子，选择元素（70+ 种），设定键长与键级，方向会针对与已有键构成的多面体配位自动计算（与已有邻居的角分离最大，自然给出直线/三角/四面体/八面体构型）；对同一锚点的重复添加会扇形展开而不是互相重叠
- **删除原子** — 点击原子并确认；原子与键会自动重新编号；也可以先用「选择原子 / 框选」选中，再点击 Delete 一次性删除所有当前高亮的原子
- **替换原子** — 点击工具栏的「Replace Atom」按钮，然后在 3D 视图中点击一个或多个原子以选中它们（再次点击取消选中）；再次点击「Replace Atom」按钮打开元素选择器（70+ 种元素），确认后把所有选中的原子替换为指定元素并保留其坐标；对晶体的基础原子同样有效
- **幽灵原子（counterpoise 标记）** — 原子可被标记为**幽灵原子**：一个有位置、却没有原子核、没有电子、也没有任何键的点——即 counterpoise（BSSE）计算中把「另一片段的原子」放在原处的标记，也是几何定义所参照的点。**添加原子**弹窗新增 *Ghost atom* 开关与 `X (dummy point)` 元素；**替换原子**弹窗新增 *Ghost* 开关，一次把整个选区变成**各自元素**的幽灵原子——这正是 counterpoise 的手势（选中一个单体、整体设为幽灵、另存为 Gaussian 输入），也是混合选区下唯一说得通的读法。幽灵原子保留其所代表的元素（Gaussian 输入写作 `C-Bq`，ORCA 输入写作 `C:`），或是一个纯占位点（`X`）；它被画成一个半透明的小灰蓝标记，而不是它所代表的原子，凡页面给出原子标签处都记作 `C*3` / `X5`。凡需要原子核的地方都把它排除在外：**成键感知绝不给它连键**（把幽灵原子叠在真实碳原子上不会长出任何棍，真实原子的键也一字不变），**分子式与电子数**把它排除（信息面板另列一行统计幽灵数），**补氢**与**分子力学**跳过它（幽灵原子不会被移动，也不产生任何力，因此放在某原子 1 Å 处的幽灵拉不动该原子），**键长**弹窗也拒绝给它建键。把幽灵变回真实原子时，只按几何重新推导该原子的键（且只推导它），与变更同属一步 undo。晶体结构拒绝幽灵原子——CIF、POSCAR、VESTA、cube 都无法表达「这个位点没有原子核」
- **幽灵原子在每种格式中都能写出与读回** — 另存为会按每种格式自己的写法写出幽灵原子，各解析器也都能读回，因此幽灵原子可以往返：`C-Bq` / `Bq`（Gaussian 输入）、`C:` / `X`（ORCA 输入）、`X`（XYZ）、`X` + `HETATM`/`DUM`（PDB）、`Du`（MOL2）、`*`（MDL Mol）、`q`（`$coord`）、`XX`（MOPAC）。没有幽灵概念的格式一律写成占位标记而不是元素符号——幽灵绝不能变回真实原子——保存对话框会在写文件之前说明这一点。没有幽灵原子的结构导出结果与改动前逐字节一致，且来历行会记录变更（`2 atoms made ghosts`、`1 ghost atom added`、`1 ghost made real`）
- **补氢（H 补齐）** — 工具栏的 **Add H** 按钮（位于「Replace Atom」右侧）会为所选原子补齐缺失的价。每个成键价低于其元素经验价的原子，恰好获得其差额，采用插件 SMILES → 3D 引擎已经在用的 X–H 键长与 H–X–Y 键角——S 补到 2、P 补到 3，C–H 键 1.09 Å、N–H 1.01、O–H 0.96、S–H 1.34、P–H 1.42，角度取该元素的理论值（C 109.47° sp3 / 120° sp2 / 180° sp；N 106.5°、O 104.5°、S 98°）——中心的杂化态由其已有键级读出（三键 → sp，双键或芳香 → sp2，类酰胺 N → sp2，其余 sp3），因此环上碳的氢落在环平面内并呈 120°，甲基呈四面体、炔氢呈直线。在该模式下点击原子（或从「选择原子 / 框选 / 替换原子」带入一个选择），再点击按钮完成放置：**原本存在的原子绝不会移动**——每个氢都围绕其中心放置，整次添加是单步 Undo，重复运行不会新增任何东西，而已经达到或超过其经验价的原子会被报告而不是被补齐。结构中放不下的氢（非键接触会短于一条化学键）会被略过并报告，而不是允许它与某个原子重叠。晶体结构与 Diff 模式下不可用
- **分子力学优化（MM Optimize）** — 工具栏的 **MM Optimize** 按钮（紧邻「Move Atoms」右侧）会打开一个面板，列出内置力场的参数、**每一个可设置参数的下拉/选择控件**（非键 / 扭转 / 异常 / 静电 / 氢键 / 金属配位 / 极化各项开关，收敛力、步数上限、时间预算、步长，以及 5–12 Å 的非键截断，默认 8）、**一键填好全部控件的三档精度预设**——*Low*（快速预优化：只算键与角，所有依赖距离与二面角的项关闭，0.4 kcal/mol/Å，120 步，1.5 s，0.30 Å 步长，6 Å 截断）、*Medium*（模块自身默认值——完整几何力场：Lennard-Jones、轻量 Coulomb 项与氢键）与 *High*（精细：包含金属配位与极化的全部项，0.01 kcal/mol/Å，2000 步，20 s，0.05 Å 步长，10 Å 截断）——任何手工改动过的值都会让选择变为「Custom」，面板还会列出本次运行会移动的原子，以及一个 **Start Optimization** 按钮。每个控件的初始值就是模块自身的默认值，预设携带的也是模块自身的取值，且每个值都会在宿主侧再次校验，因此任何选择都不可能把一次运行带入不合理状态。在该模式下点击原子即可选择（可切换），或从「选择原子 / 框选 / 替换原子 / Add H / Rotate Group / Move Atoms」带入一个选择——选择在双向切换中都保留，面板会跟随它。**这是一次「受约束的松弛」，不是「孤立的松弛」**：只有被选中的原子移动，但每一个固定原子都是完全的参与者——它的键、角、二面角、平面外约束与非键接触都在拉扯选区，只有当某个项涉及的原子全部固定时才会被丢弃（它不可能改变）。按下按钮后恰好只有这些原子被松弛，而**其余每个原子都逐位返回原值**，因此既可以让一个片段在结构其余部分的存在下安顿下来，也可以松弛一个整体扭曲的结构，而完全不必担心未被选中的部分。面板随后报告能量 前 → 后（kcal/mol）、步数与停止原因（收敛 / 能量收敛 / 步数或时间预算 / 再无增益）、最终的 max 与 RMS 力、最大位移，以及每个项有多少个处于活动状态；状态栏给出同样的摘要。**一次运行是宿主分块推进的会话，面板看着它工作**：结构随每个分块实时跟进（只重绘球体，以及与移动原子相关的键），两块画布随运行进程绘制能量与 max/RMS 力对步数的曲线，一行进度文字给出步数、能量、其变化量与受力。**一个 Pause / Resume 切换按钮**（与振动面板相同的单按钮模式）就地停下并继续运行——会话被保留，因此 **Resume** 从完全相同的几何、步长与共轭方向继续，暂停期间的时间不计入墙钟预算；标签在点击的那一刻翻转，因此面板绝不会停留在它无法离开的状态。**Restore initial** 把运行开始时的几何写回，作为单步 Undo。运行本身达到的几何也是单步 Undo，且快照在第一个真正发生变化的坐标之前拍摄——因此一次「报告其起始几何未变」的运行会让 Undo 栈保持原样。因为别的原因离开该模式——切换工具、帧、分子、Undo、Import、比较——都会**中断**运行，因此它绝不会继续往用户已经离开的结构里写入坐标，而被替换掉的运行的迟到答复会被忽略（每条消息都带该运行的 id）。用 × 关闭面板只是把它隐藏：Convergence / Modes / Route 按钮旁边的小 **MM Optimize** 按钮可重新打开它，并在有运行进行时兼任其状态显示（「running – step N」/「paused」）。整次运行是**单步 Undo**，之后被移动区域的键会由「载入移动后坐标时同一个感知模块」重判（即 Move Atoms 释放路径），因此一次松弛只会产生「重新载入也会认可」的几何。力场属于插件本体（`src/tools/mmOptimizer.ts`，零依赖）：谐振键伸缩，其参考键长取自 SMILES → 3D 引擎与 Add H 使用的同一张表（按元素对与键级给出 `r₀`，含芳香）；谐振键角弯曲，朝向中心元素的理论 sp/sp²/sp³ 角（对*被压缩*的角超过 30° 后无作用力，因此张力环保持原样，而对*被张开*的角可一直拉到 60°，因为张得太开的角是需要纠正的画法）；sp³–sp³ 键上三重重叠式扭转，涉及 sp² 中心或双键/芳香键处为二重平面式（当四原子集合由氢定义时没有该项）；三配位 sp² 中心有平面外约束；Lennard-Jones 非键项，σ 取自范德华半径，1–4 对折半，在 6 到 8 Å 之间平滑关闭，另有一道 1-3 冲突保护：若一对原子接近到其范德华和的 0.45 以内，就以纯排斥重新纳入。除成键项与 Lennard-Jones 项之外，力场还带有更高精度档位才启用的**轻量**补充：**Coulomb** 项，其分电荷来自 Pauling 电负性差（按键级加权、按结构整体中和，因此无需形式电荷），配合与距离相关的介电常数 εᵣ = r 以及同一套配对表、渐变与 1–4 折半；**氢键**项，在 ≈1.6–2.8 Å 的 H···A 对上叠加一个 10-12 势阱，并施加 D–H···A 角窗口（给体与受体为 N/O/F）；**金属配位**项，把金属的配体–金属–配体角弯向其配位数的理想几何（直线型、三角型、四面体型，d⁸ 金属为平面四方型，五或六个配体时用 90° 与 180° 双极小值的顺/反式项），而把金属–配体距离留给键项；以及一个**阻尼诱导偶极极化**项（∝ 1/r⁴，标定为 α = 1.45 Å³、q = 0.4 e、r = 2 Å 时约 −0.5 kcal/mol）。逐元素表覆盖**第 4、5 周期的主族元素**（Ga Ge Rb Sr In Sn Sb 及其同族类比者）以及第一、二周期的元素，新增项都遵循一条成文规则——同族取最近表项的值——因此没有任何既有值被改动，重元素也不再落到表的默认值上。优化器是预处理共轭梯度，且**每一步被接受的移动都必须降低能量**，因此一次运行是单调的：它返回的任何结果——收敛、步数用尽或时间用尽——都不会比给它的输入更差。单次运行受 3000 原子、500 步与 5 秒限制。晶体结构（力场没有镜像原子）与 Diff 模式（有两个结构需要选择一侧）下不可用。验证包括：梯度与中心差分的核对（最差相对误差 1.3e-8）、受约束松弛的证明（一个孤立的可动氢会落在与固定碳之间的表键长上，三个 H–C–H 角与冻结氢恰好给出 109.47°，一个端基碳能感受到被冻结骨架的扭转）、一个 19 文件语料库中每个结构都有改善，以及对真实客户端脚本的桩驱动运行。流式、预设与截断的工作在此基础上另做验证：会话重构**之前**那一版的答案被导出为 16 个结构 × 8 组选项，重构后每个字段（含坐标）都**逐字节一致**；分块大小扫描（每块 1 / 3 / 7 / 1000 步）以及在其中每次运行中途做一次暂停/恢复，都精确复现单次调用；整段流式运行驱动编译后的扩展代码——mmProgress → mmResult——与模块的单次调用逐字节一致，且 Pause 确实停住泵、Resume 继续它、Cancel（或面板被销毁）永久停止它；新客户端代码以桩 DOM 驱动（start → 实时坐标 → pause → resume → result → restore，首次写入时的一次 undo 快照，一个不会写入任何东西的过期 run id，以及导航取消）。截断梯度在每个渐变窗口（5、6、8、10、12 Å；最差相对误差 2.3e-7）内仍与中心差分吻合，默认截断逐位复现先前行为，三档预设也相互测过（medium 51 步到 max|F| 0.027 kcal/mol/Å，high 81 步到 0.0093，low 29 步到 0.377；medium 与不给选项的答案完全一致）。**一次运行在进行中仍可被操控，而查看器的其余部分继续可用**，并且「面板的 Pause 切换、它的曲线、它的状态」可能与宿主实际行为不一致的三条路径都被堵住：点击 Pause 时已在途的分块会被**记录**（其步数进入曲线，因此暂停运行的能与力历史是完整的）而坐标被扣下——于是用户按下暂停时看到的几何就是他一直看的几何，直到宿主自己的暂停答复带来它真正停下的快照；宿主的泵**绑定到它被调度时所属的那次运行**，因此一次被取消并被替换的运行所遗留的分块链再也不能醒来推进新的运行；对宿主已不再拥有的运行发出的暂停或恢复会被**答复**（作为一次结束）而不是被忽略，因此按钮不会停在没有任何东西会来确认的状态。曲线绘制在画布自身的盒子内，并在布局回流时重绘，因此窗口缩放、面板移动或「关掉再打开」都不会让它们被右边缘裁掉。用 × 关闭面板后运行继续，重开按钮显示其状态；运行中途改参数或点选新原子会在不打扰运行的前提下围绕会话重建面板；在 **Add H** 请求仍在途时启动运行会被拒绝，而稍后落地的 Add H 答复会先结束它、否则就会在运行之下追加原子。其他任何写结构的东西也斗不过一次运行：启动运行会先停掉**帧自动播放**，而**简正模式播放**会在它播放前先结束运行（与「运行启动时结束播放」互为镜像）——两者都在运行捕获 **Restore initial** 将要写回的几何**之前**完成，因此二者绝不可能同时拥有坐标。运行进行中时，面板会同时标注**两套原子集合**——运行自身的那一套（其按钮作用的对象）与视图的当前选择（下一次运行的那一套）——并且它的分块是一份固定的*工作量*而非固定的趟数，因此在大选区上点击 Pause 也能及时落地。
- **导入结构** — 点击工具栏的「Import Structure」按钮（位于 添加/删除/替换原子 旁边）；一个对话框先同时询问**副本数量**（在输入框中键入，预填 `1`）与**来源**——**从文件** 或 **从 SMILES**——因此数量在结构之前就确定下来。*从文件* 会把第二个结构从任一受支持文件（XYZ、GJF、MOL2、PDB、log/out 轨迹、晶体格式……）载入当前视图。多帧文件（优化日志）会打开帧选择器；晶体文件（CIF/VASP/Cube/VESTA）以分子形式贡献其晶胞原子。*从 SMILES* 打开一个输入框，**边输入边校验**——在字符串不是合法 SMILES 时显示解析器自己的报错信息（含出错位置）——并在按 Enter 时构建 3D 结构：本地生成带显式氢的构象，具有正确的键级、芳香环感知、顺/反双键与四面体立体化学（`@`/`@@`），由内置的零依赖 SMILES → 3D 引擎（`src/tools/smilesTo3d.ts`，参考实现 `test/smiles_to_3d.py` 的忠实移植）产生，因此不需要 Python，也不需要任何外部化学信息学工具包。两种来源此后行为一致：导入的结构以包围球分离检测（两侧各加 2.5 Å 再加 2 Å 间隙）被放置在已有结构旁边，因此绝不重叠；导入到空视图时结构会以原点居中。导入可用相同或不同来源重复任意次——每份副本沿 +X 以统一步距进一步堆叠，且多副本导入在一趟中完成布局与放置，因此副本数量不会让放置工作量成倍增加。键的合并很聪明：若两侧都没有显式键，则对合并后的结构重新感知连接性；若任一侧有自带键（MOL2/PDB/GJF connect，或任何 SMILES 导入），则保留已有键，只补充导入部分的键（芳香键携带 1.5 级，因此渲染为 1 实线 + 1 虚线）。每次导入是单步 undo，相机自动缩放到合适视野，保存时导出合并后的结构
- **选择原子** — 点击工具栏的「Select Atoms」按钮会显示一个浮动面板（左下角），可在其中键入原子编号（从 1 开始）、区间（例如 `3-10`）或元素符号（例如 `C H`）来以黄色高亮特定原子；按 Enter 或点击 Select 生效，Clear 重置；也可以直接在 3D 视图中点击原子来选择；再次点击已高亮的原子会取消选中（可切换）；在「选择原子」模式下选择会累积；该面板是真正浮动的覆盖层，不影响分子视图尺寸；面板随工具栏按钮切换，切换到其他模式时自动隐藏；在「选择原子」与「框选」之间切换时选择会保留
- **框选** — 独立的工具栏按钮；在该模式下按住左键拖出一个矩形，松开即选中框内所有原子（原子被投影到屏幕坐标并与矩形比较）；选择以与「选择原子」相同的黄色高亮累积，且在「框选」与「选择原子」之间切换时保留；切换到任何其他模式都会清除选择
- **撤销** — 工具栏按钮回退最近一次编辑/删除操作；快照在编辑弹窗打开时（含「添加原子」）以及在删除/清除无序/批量删除执行之前自动拍摄；栈上限 50 条，文件重新载入时重置
- **另存为** — 以 14 种格式导出修改后的结构：XYZ、Gaussian GJF（保留原文件结构）、CIF（分数坐标）、VASP POSCAR（Direct 分数坐标）、Gaussian Cube（含体积网格）、VESTA（P1 对称，分数坐标）、Turbomole Coord（Å→Bohr 换算）、ORCA Input、MOL2（含键级）、MDL Mol、PDB（含 CONECT 记录）、MOPAC Input，以及——仅对晶体结构——**XYZ / GJF − displayed cluster（当前显示团簇）**：把屏幕上的区域当作分子写出（QM 计算要的笼、孔或片段），而不是当成一个晶胞。标题行记录来源文件、**晶胞参数**、**超胞 boundary** 与原子数，因此文件能说清自己取自哪个结构的哪一块；**Keep whole molecules** 开关保留显示时的整分子规则（质心落在 boundary 内的分子完整写出），取消勾选则在 boundary 处严格切割（界外原子丢弃，跨晶面的键不写）。GJF 输出保留原 Link 0、route、标题（写在出处行之下）、电荷/多重度以及 connect 段之后的内容；GJF 原子坐标采用 GaussView 风格的 8 位小数对齐格式；connect 段包含所有原子行；CIF/VASP/Cube/VESTA 仅在存在晶体数据时可用
- **保存的标题行记录结构来历** — 凡是有标题行/注释行的格式（XYZ 注释行、GJF 标题段、MOL2 分子名、MDL 标题行、PDB `TITLE` 记录、MOPAC 标题行、ORCA 注释行，以及 POSCAR/CUBE/VESTA 的头部），写出时都带一行来历，例如 `MV: from opt.log, frame 12 of 40, Step 12; 2 atoms deleted, translated, displaced along mode 7 (38.2 cm-1) +0.45 A`：来自哪个文件、来自构型优化文件的第几个结构（VMD 脚本则写第几个分子）、此后在这里做过什么——增删改原子、补氢、平移、旋转、调整键长/键角/二面角、MM 优化、清除无序、合并进来的文件、超胞范围——以及沿简正模式的位移（模式编号、频率与幅度）。这些信息在每一步操作发生时记录，并且随 **Undo** 快照一起回滚：撤掉某次编辑，标题里对应的那句也一起消失；GJF 会把原有标题段原样保留在来历行之下，重复保存不会层层叠加。没有标题行的格式（Turbomole `$coord`、CIF）保持原样；文件**原有的标题**只有在它描述的就是屏幕上这个结构时才以 `|` 追加在后面——轨迹的标题属于文件打开时的那一帧，所以保存 79 帧日志的第 79 帧不会以第 1 帧的注释收尾；措辞简短且为 ASCII，整行上限 190 字符
- **连续编辑** — 完成一次编辑后，查看器保持当前编辑模式以便反复调节
- **取消/撤销** — Cancel 按钮在确认编辑前恢复原始坐标
- **绕轴旋转基团** — 工具栏的「Rotate Group」模式（与 View / Bond Length / Bond Angle / Dihedral 归为一组）允许选择一个原子组（在 3D 视图中点击，或键入原子编号/元素符号，含 `H581` 这样的 `元素+编号` 简写），并由 2 个以上轴原子定义一条旋转轴，然后让该原子组绕轴旋转 0–360°。角度可通过拖动平滑滑块或在数值输入框中键入数值并按 Enter 来调节，两者保持同步。轴以穿过分子的橙线可视化。轴原子会接受化学校验：其中任何成键对都必须是单键（双键/三键被拒绝，因为它们无法在不破坏的前提下旋转）。首尾两个轴原子定义方向；中间原子仅参与校验。适用于多原子共线轴（例如 `C32 C33 C43 C53`）或简单的两原子键轴（例如 `C33 C43`）。轴原子保持固定，所选组的其余部分绕该直线旋转。每次角度调节是单步 undo（Undo 回退一次调节并更新滑块显示），而面板的 Reset 按钮一次性回退整个旋转会话
- **移动原子** — 工具栏的「Move Atoms」模式（紧邻「Rotate Group」）用于平移与旋转所选原子组。可在该模式中点击选择原子（可切换）、在「选择原子」面板中键入，或从「选择原子 / 框选 / 替换原子 / Add H / Rotate Group」带入选择。在该模式下，普通操作作用于整个视图，Ctrl 修饰的操作作用于所选片段：**左键拖拽 / 方向键**旋转视图；**右键拖拽**平移视图（光标锁定）；**Ctrl+左键拖拽 / Ctrl+方向键**让片段绕其自身中心旋转，拖拽方向与视图旋转一致而不受朝向影响；**Ctrl+右键拖拽**平移片段——移动平面穿过光标下的那个原子，因此该原子在拖拽中始终精确停在光标下（从空白处开始时使用片段中心平面，片段以 1:1 速度跟随光标）。松开 Ctrl、鼠标键或全部方向键即结束。拖拽过程中，已有键随原子拉伸，且只有与移动原子相关的键网格被重建（每帧同步，因此键绝不落后于原子）；会话结束时（释放时）重新感知一次连接性与键级，只重判「触及被移动原子的键」及其邻居的键——代价随被移动区域而非分子大小增长，且结果与完整重感知一致（对退化性的既有价态违规有一道安全检查会回退到完整检测器）；释放刷新同样只重建受影响区域的键网格而不是整个场景（存活的网格保留其 GPU 资源并重新索引进新的键表，只有在安全检查触发时才做完整重建）。带显式键块的文件如 MOL2/PDB（以及带 connect 段的 GJF）保留其自带键与键级——只有在释放时端点被拉到成键截断之外的那些键才会被删除（只删不增：不重新排序，片段靠近时也不新增键）。每次拖拽或按住方向键是单步 undo。晶体结构下不可用

#### 优化轨迹浏览（LOG 文件）

- **逐帧步进** — ◀ Prev / Next ▶ 按钮在优化帧之间步进
- **跳转帧** — 直接输入帧号跳转
- **自动播放** — 以 350 ms 每帧的间隔依次播放全部帧；播到最后一个结构（收敛后的几何）即停止，不再跳回第一帧循环，在该帧再点播放则从第一帧重新开始
- **路由关键字面板** — Gaussian 日志与 ORCA 输出都会与 Convergence、Modes 并列显示第三个可折叠面板按钮（⚙️ Route），展示该次运行关键字的解析结果：对每一个 `--Link1--` 作业段，给出原始关键字行以及解析后的关键字分解（打印级别、方法/基组、每个关键字及其选项列表，对要求优化的段加 `opt` 徽标）；ORCA 则给出被回显进输出的 `!` 简单输入行（从 `INPUT FILE` 段恢复，折行时拼接），每个词元算一个关键字，优化类运行带 `opt` 徽标。分词器能正确处理括号与引号（`scrf=(SMD, Solvent=Acetone)` 保留其选项组，`external="..."` 的引号值保持完整，`iop(5/17=17)` 解析为关键字 + 选项），续行会被拼接（`... opt` / `=(calcall) freq` 把选项折进 `opt`），全角字符会被规范化。同一份解析结果也驱动 Convergence 面板的显隐，因此两者始终一致。ORCA 的简单输入没有 `key=value` 选项语法（选项写在 `%block` 指令里），因此它的词元作为无选项的关键字列出
- **收敛面板** — Gaussian、ORCA 与 xtb 的结构优化运行会显示一个可折叠面板（📈 Convergence），给出能量、力（max + RMS）与位移（max + RMS）随优化步数的曲线，对应 GaussView 的收敛图；阈值线（琥珀色虚线）画的是**文件自己写明的目标**：ORCA 输出在 `Convergence Tolerances` 段里给出它们（`TolMAXG`/`TolRMSG`/`TolMAXD`/`TolRMSD`，该段没打印时退回用收敛表自身的 `Tolerance` 列），因此被要求收得更紧的运行显示的就是它真正被考核的那条线；Gaussian 日志只在同一张表的 Tolerance 列里给出它们，故以其成文默认值代替（力 0.000450/0.000300，位移 0.001800/0.001200）；两者都没有的程序（xtb 的 `gnorm` 判据、轨迹文件）只画曲线，不借用阈值尺度——另有竖直拖拽把手可调整面板高度，以及关闭/重开切换。只有该次运行真的在优化时面板才出现——对 Gaussian 而言，带 `force` 关键字的单点作业（或 `freq` 作业）会打印同样的收敛表但并不优化，因此解析器会检查每一个 `#...` 路由段是否含 `opt` 关键字（`opt`、`opt=...`、全拼 `optimization`、`freq=calcall`；大小写不敏感的整词匹配），只要没有任何路由段含它就完全不显示 Convergence 按钮；ORCA 只在松弛循环中打印该表，因此纯频率运行不会出现面板。力与位移图使用对数 Y 轴——刻度显示 log₁₀ 值并带一个旋转的 `log₁₀` 轴标签——这样即使早先的步数把最终值压得很小，收敛阈值线仍然可读（能量保持线性轴，因为它可能为负）。能量取自 `SCF Done` 行，或 `opt=external` 运行中优化器的 `Energy= ... NIter=` 行（例如 Gaussian + xTB/XO 联合作业，它们不打印 `SCF Done`），以及 ORCA 的 `FINAL SINGLE POINT ENERGY`（每一步的能量就是针对它所松弛的那个几何打印的那一个——ORCA 先打印几何，因此解析器把能量反向绑定到正在等待它的那一帧）。对 xtb 轨迹，把每一步注释行中的梯度范数（`gnorm`）绘制为 RMS Force 曲线。能量轴在图表标题中打印一个共同的整数基数，刻度上只显示变化的小数部分（例如标题为 `-6912 +`，刻度为 `.98 1.18 … 1.81`）；小数位数随数据自适应，使五条网格线始终可区分且单调，即便某次运行的范围跨越整数边界、或缩小到几毫哈特里。悬停数据点会显示提示框，给出步序号、该点的值，以及在力/位移图上给出该步的能量，各自打印足够的小数位以区分相邻步
- **简正模式振动** — 带频率计算（`freq` / `CalcAll`）的 Gaussian 日志，以及带 `VIBRATIONAL FREQUENCIES` / `NORMAL MODES` 成对输出的 ORCA 结果，都会显示一个可折叠的 Normal Modes 面板（🎵 Modes），列出所有简谐频率（cm⁻¹）及其对称性标签；ORCA 会打印全部 3N 个模式，因此它按 `0.00 cm**-1` 列出的 3N−6 个平动与转动模式会被剔除（只有频率**与**位移向量同时为零才剔除，所以一个频率恰好舍入为 0.00 的软模式仍会保留），虚频则保留 ORCA 打印的负号，面板用 `*` 标出它，与标注 Gaussian 虚频的方式完全一致；点击某个频率即播放振动动画（通过 `requestAnimationFrame` 达到 60 FPS，振幅归一化到 0.35 Å，周期 3 秒）；虚频（过渡态）以 `*` 标出；当别的东西接管结构时动画即结束——切换帧、切换分子、Undo、Import、Diff、Add H 答复落地、编辑弹窗打开，或分子力学运行启动——而启动振动动画会结束帧自动播放，而不是被它的下一个 tick 停掉（同一时刻只有一样东西在写坐标）。每个模式行现在带**两种彼此独立的选择**：行尾的小圆圈把该模式指定给位移滑条，行本身（或它的 ▶）播放振动动画——播放不会再改写滑条指向的模式，选取模式也不会启动动画。指定模式后，**Manual Displacement** 把它真正写进坐标：拖动滑条——量程是该模式的**比例系数，从 −1.00 到 +1.00**，1.00 表示"该模式中位移最大的那个原子走 1 Å"、其余原子按模式自身的相对权重同步移动（权重与相位原样保留），系数旁边同时显示它对应的最大位移距离——每一步都把结构设为 *平衡构型 + Σ 系数ₖ × 模式ₖ*。**新模式是叠加进畸变、而不是替换上一个**：面板列出每一个分量（`#2 +0.50 · #1 +0.30`），每个有贡献的行都带标记，把某个分量的滑条拖回 0 即把它移出，**Reset** 清空全部；列表旁给出整个组合的最大位移——它不等于各分量之和，因为模式之间可能部分抵消。每次写入都从平衡快照重算、绝不累加，所以来回拖动不会漂移，同一个滑条值永远对应同一个几何；仅"选取"某个模式不会改写任何东西。只要不替换或重写结构，会话就一直保留——旋转视图、量键长键角、关掉再打开面板、保存——**Save As** 会以任意格式写出被位移后的几何（对话框会注明是哪个模式、多大幅度）。此后播放任何模式（包括它自己）都以**未位移的平衡构型**为平衡位置摆动，停止时把位移后的几何写回来，所以动画始终显示计算给出的那个模式，也不会丢掉已有的调整；**Reset** 精确回到平衡构型，并连同该会话的 undo 快照一起撤掉。任何坐标编辑都会先关闭会话并保留它产生的几何（编辑从被位移后的结构开始），而切换帧、切换分子、Undo 会丢弃会话；分子力学运行进行中拒绝滑条（那批坐标归它所有），Move Atoms 会话则在滑条接管前先提交

#### 结构对比（Structure Diff）

- **Diff 按钮** — 把当前结构与另一个分子文件比较（任意受支持格式，无需格式一致）
- **帧选择** — 若比较文件是优化 LOG/OUT，用 QuickPick 选择要比较的帧
- **骨架校验** — 通过递归回溯的原子匹配做图同构（元素 + 度数 + 邻居签名，配合最受约束变量启发式）；使用插件自身的键感知算法（而非文件原有的键级），以确保即使文件对同一连接性给出不同键级，骨架比较仍然一致；匹配只考虑**连接性**（是否成键），不考虑键级——同一位置在一个分子中是单键、在另一个中是双键，不影响骨架匹配
- **并排视图** — 左侧视口显示原始分子，右侧显示比较对象；每一侧都有**独立**的旋转、平移与缩放控制。左键拖拽旋转、滚轮缩放、中/右键拖拽平移**你抓住的那一半**——平移是光标锁定的（抓取点精确停在光标下）且限定在该视口内，因此拖动右侧结构绝不会移动左侧，反之亦然
- **正确的宽高比** — 每个视口使用与半宽匹配的独立相机宽高比，避免横向变形
- **构象差异** — 骨架匹配时，在两侧以橙色高亮键长、键角与二面角的差异；完整列出全部差异（不截断）；每个差异显示相对于平均值的百分比变化
- **分组内坐标统计** — 骨架匹配时，差异列表上方会显示一个可折叠的 **Internal Coordinate Statistics** 区块，汇总整个比较：对键长（Å）、键角（°）与二面角（°），为每个元素组成分组（例如 `C-C`、`C-C-H`、`C-C-C-H`）给出 **n、rms、maxd（最大 |Δ|）、mean 与 std**，另加一行合并的 **All**。统计覆盖所有被比较的坐标，而不仅是超过阈值滑块的那些，因此在拖动滑块时统计值保持不动（滑块只改变高亮与逐条列出的内容）。数值是有符号差（左 − 右），二面角的差折入 [−180, 180]，因此一条链与其反向不会互相矛盾；rms = √mean(Δ²)，maxd = max|Δ|，mean = mean(Δ)，std = 总体标准差。分组键是规范化的，因此无论原子编号如何，同一个坐标总是落在同一个分组。点击区块标题即可折叠或展开
- **可调阈值** — 差异结果面板中的百分比滑块（0–20%）让你交互式调整「算不算差异」的灵敏度；三种类型使用同一个相对差异公式（Δ/avg × 100%）；对平均值非常小的接近 0° 的二面角有一条回退分支以保证稳健性；拖动滑块时 3D 高亮与差异列表实时更新
- **差异信息面板** — 浮动面板按幅度排序列出全部差异，附原子标签与两侧结构中的数值；可通过「Show Results」按钮关闭与重开；差异结果始终可访问，不会被测量信息覆盖
- **自动镜像选择** — 在一侧选择原子做键长/键角/二面角测量时，另一侧通过原子映射自动选中对应原子，并同时显示两侧数值与 Δ 差值
- **晶体对比** — 晶体按**晶胞**比较，而不是按屏幕上与文件里恰好各自持有的两份原子表：显示的 2×2×2 超胞（2352 原子）与文件的单胞（294 原子）只是同一个对象的两种画法，因此比较先核对**六个晶胞参数**（a/b/c 相差 ≤ 0.01 Å、α/β/γ ≤ 0.05°；真的不同时以 **Unit Cells Differ** 列出两侧参数），再匹配**基础晶胞的周期性图**（同一套骨架匹配，尊重每根键的周期像；图不同时以 **Skeletons Differ** 给出基础晶胞的原子/键数）。因此同一个文件与自身比较，在任何超胞尺寸下都是 **Structures are identical**。键长、键角与二面角都**沿各自键的周期像**测量，跨晶面的键得到的是真实键长而不是一个晶胞；左侧同一基础位点的每个显示像都带橙色标记，右侧标记对应位点，自动镜像选择仍然把两侧配对。晶体与分子比较会明确提示**不可比较**，而不是含糊地报「骨架不同」；编辑或切换显示样式后比较在原位重算

#### 晶体结构查看

- **晶胞显示** — 打开 CIF、VASP、Cube 或 VESTA 文件时，晶胞以线框轮廓渲染，显示晶格边界；当显示区域跨多个晶胞时，线框按胞细分为网格（内部晶胞面线比外框更淡），且在原点相交的三条棱带 **a**/**b**/**c** 标签（颜色与坐标轴指示器一致）
- **超胞边界控制** — 6 个输入框（a/b/c 的最小/最大值）用于定义要显示的晶体区域；支持小数取值（步长 0.1）以选取部分晶胞；该视图是晶体的「剖切」，而不是对晶格的形变；非整数边界（例如 -0.4 到 0.4）不会产生镜像——范围之外的原子直接不显示。**分子保持完整**：每个有限分子（非周期性延伸者）先沿键网络解包（unwrap）为同一像内的连通整体，再按质心放置副本——质心落在边界内的每个像都会画出整个分子，因此跨越边界的分子完整显示（并照旧延伸到边界之外），而不是被晶面切断；真正周期性的网络（金刚石、聚合链）仍在边界处通过分割键截断。**一键预设**：输入框旁提供 2×2×2 / 3×3×3 / 4×4×4 按钮，以及 **Reset**（回到 0–1）与 **Align Origin**（保持跨度、把最小值归 0）。**视图不跳动**：修改边界时保留原先的旋转中心，因此正在看的原子在屏幕上的位置不变；当前原子选择按「基础原子索引 + 像」重映射到重建后的结构——原来选中的仍是那些原子。**规模保护**：边界改动若使显示原子数超过 20,000，将被拒绝并提示、输入框回退原值，而不是卡死编辑器
- **原子身份与晶胞信息** — 悬停晶体原子时，除笛卡尔坐标外还显示其分数坐标、所处的周期像（`cell (i,j,k)`），以及存在时的不对称单元位点标签与占据率，从而区分同一基础原子的不同像；信息面板新增晶胞**体积**（Volume）、每胞分子数 **Z**（由放置的有限分子与孤立原子计数；纯周期网络无此行）与晶体**密度**（按占据加权、由原子质量表计算）
- **晶体局部修改** — 在显示的晶体原子上做的几何编辑（键长、键角、二面角、分组旋转）会**回写到晶胞**：被编辑的那个像的位移会换算成它所在基础位点的分数位移，因此改动在另存为 CIF/VESTA 后依然存在，并在下一次重建后出现在每个周期副本上；一步 Undo 即可撤回（undo 还会恢复编辑时所处的超胞边界）。同一基础位点在屏幕上出现多次时，编辑只**作用一次**，而不是每个副本各一次。被画成**跨晶面短桩**的键会被弹窗拒绝——屏幕上两个原子本来就相隔整个晶胞，它们之间的距离不是键长；扩大超胞、选一张两端都在同一像里的原子即可编辑。**添加原子遵循同一规则**：新位点先 wrap 到晶胞内，其键记录它指向哪个相邻像，因此永远不会画出横穿晶胞的长键
- **最小镜像约定成键** — 晶体键按穷举最小镜像搜索判定，而不是四舍五入分数差：只有当某个方向重复一次所跨越的晶面间距小于截断距离时才需要尝试该重复，因此普通晶胞枚举经典的 27 个像，窄或斜的晶胞自动扩展到 5×5×5。键级来自**分子路径所用的同一套感知模块**（含价态松弛与芳香 1.5），所以晶胞里的分子画出的键与它作为 XYZ 时一致；跨越晶面的分子先解包再整体判定，而不是被切成两段。跨越晶胞边界的键仍绘制为延伸到视图之外的短分割键（而不是横跨整个视图的长键）
- **随动坐标轴指示器** — 晶体结构的左下角坐标轴显示晶格 a/b/c 轴，与主 3D 视图同步旋转（见上文通用的「方向坐标轴指示器」，它覆盖分子结构）
- **分数占据渲染** — 具有部分占据（晶体无序位点）的原子按与其占据值成比例的部分着色渲染
- **清除无序** — 一键按钮清理无序位点：占据 < 0.5 删除，> 0.5 置为 1，= 0.5 保留
- **晶体保存** — 把修改后的晶体结构导出为 CIF（分数坐标）、VASP POSCAR（Direct 分数坐标，按元素分组）、Gaussian Cube（体积网格）或 VESTA（P1 对称，分数坐标）。四种格式写的都是**晶胞**——即会话中编辑、每次重建都据此重绘的 `baseAtoms`——且每个文件都把这个事实说出来：CIF 以 `_symmetry_space_group_name_H-M 'P 1'` 加唯一恒等操作声明，并在 `_audit_creation_method` 里写明所写的是「对称已展开的整胞」而不是不对称单元，因此把文件读回再存一次，原子数不会翻倍膨胀；分数坐标 wrap 到 [0,1)，并补上 `_cell_volume`、`_chemical_formula_sum` 与 `_atom_site_site_symmetry_multiplicity`。POSCAR 的注释行写明它列的是晶胞，原子也取自晶胞本身（2×2×2 的视图保存下来仍是单胞，不会有原子落在声明盒子之外）；Cube 是真正的 cube——正的原子数（坐标为 Bohr）、以晶胞为网格盒子、网格点数在对话框中可选（默认每轴 40，数值全零）——不再是会被读成轨道数据的负原子数 + 2×2×2 空壳；VESTA 的坐标同样 wrap。CIF 与 VESTA 携带的是**编辑后**的晶胞（见「晶体局部修改」），屏幕上改过的几何就是文件里的几何；屏幕上的区域则通过另存为的**团簇格式**导出为分子

#### VMD TCL 脚本可视化

打开一个 `.tcl` 脚本时，它载入的结构会**按脚本要求的方式**渲染，而不是套用查看器自己的默认值：

- **脚本就是一套显示配方** — 它执行的每一个 `mol new` 都成为一个分子（当不止一个时，工具栏出现 ◀ n/N ▶ 选择器，正如力误差系列脚本每轮优化载入一个 PDB 的情形），并且每一种表示都被复现：样式与参数（`mol representation CPK 0.97 0.10`、`Surf 1.0 1.0`、`VDW`、`Licorice`、`Lines`）、颜色模式（`mol color Beta` / `ColorID n` / `Name` / `Element`）、选择、材质以及任何 `mol scaleminmax` 色标范围。`mol modstyle` / `modcolor` / `modselect` / `modmaterial` 修改已有表示，`mol delrep` 删除一个——包括 `mol new` 创建的默认表示，因此表示的编号索引指向脚本作者所指向的东西。
- **选择针对已载入的原子求解** — `index` 列表与 `to` 区间（0 起；越界索引被丢弃）、`element` / `name` / `type` 符号列表、`not` 取反，以及 `all` / `none`。需要查看器并不保留的残基元数据的选择（例如 `resname BIB`）会被放宽为整个分子并**予以报告**（状态栏与信息面板），绝不静默套用。
- **尺寸保持球棍风格，并采用脚本的相对比例** — 脚本的表示参数是相对于参考渲染器自身默认值的，那里 `CPK 1.0` 意味着范德华球（空间填充模型）。若把它们当作绝对半径套用，每个结构都会变成一整块实心球，把键与分段全部遮住，因此查看器把脚本的**比例**套用到自身几何上：球半径 = 查看器自身半径 × 脚本球缩放（普通脚本得到普通结构的尺寸，`CPK 0.7` 更小，`CPK 0.1` 极小），键半径同样以参考默认值 0.3 为基准重标（`CPK … 0.3` 是通常的棍粗细，`0.1` 更细，`Lines` 是一条细线）。任何「物理量」都*不*缩小：`Surf` 探针半径是 Å 值，因此表面位于范德华表面之外那么远——这正是 `show_eggs.tcl` 读起来像是原子周围的彩色包络的原因。没有任何表示覆盖的原子与键完全不绘制，因此部分选择给出部分视图——除非把 **Display** 选为空间填充或 stick：那是「整个结构」的模型而不是某一条表示的模型，因此会画出每个原子（脚本的颜色与曲面仍然生效）。
- **颜色模式** — `Beta` 把逐原子标量按参考渲染器实际渲染时的 RGB 斜坡上色：红 → 绿 → 蓝，中点位于 0.5，但两个非活动通道有 15/255 的底，活动通道以标称速率的两倍爬升，因此低端是 `#FF0F0F`、中点是 `#0FFF0F`、高端是 `#0F0FFF`，而不是其色标名义上给出的纯原色。范围在脚本写出了 `mol scaleminmax` 时取自脚本（写出范围被当作「该列是标量通道」的无条件声明，因此常量列也会被上色——`force_error_delta.3.pdb` 在 −10…+10 色标上全为零，就渲染在斜坡中段），否则取自该表示的实际数据范围（脚本写死的 `0 – 100` 色标确实会把 0.04–3.92 的力误差挤到红色端——那就是脚本做的事，信息面板会显示实际使用的范围）。`ColorID` 使用参考渲染器 33 项色板的实测取值，包括 `color change rgb` 带来的脚本覆盖。`Name` / `Element` 使用参考渲染器的元素颜色——碳是柔和的青色 `#40BFBF`，不是灰色。
- **`Surf` 表示是真实曲面** — 一个 `mol representation Surf <probe> <isovalue>` 片会变成它所选原子的溶剂可及等值面，以 marching tetrahedra 在 0.45 Å 网格上提取。曲面片以半透明绘制：像 `show_eggs.tcl` 这样的脚本会叠十四层，不透明的话除了最外层什么都看不见。它们随分子组一起运动，因此与原子完全一样地旋转与平移。一个十四片、72 原子的场景约在 50 ms 内构建约 164 000 个三角形。
- **蛋（egg）列表** — 当结构旁边放有该工具链的蛋元数据（一个 `eggFile*` 命名每个蛋及其包裹的单元编号——`U12F1 : 1,10-15`——外加一个含全部单元原子的 `unitfile`）时，会出现一个 **Eggs** 按钮，面板列出脚本绘制的每一个蛋：每行以该蛋的颜色作为左侧强调条，随后是它的名称与所含单元（`units 1, 10–15`）。初始显示所有蛋；取消勾选某行即隐藏恰好那一个蛋的曲面片，重新勾选则恢复，而 **Hide all / Show all** 按钮一键切换所有蛋（是补全部分选择，而不是重置它）。这些选择会在场景重建（编辑、undo）后保留，并在曲面片列表变化时（切换分子或帧）重置。元数据文件是可选的——没有它们时曲面片就不带标签渲染。当目录里存在**多个** `eggFile*` 文件（其他设置或轮次留下的）时，通常只有一份描述该结构：查看器会把每个候选与脚本的曲面片打分并询问使用哪一个，每个条目显示它匹配了多少片（`14/14 patches match` 对 `0/14 patches match`），因此正确文件一目了然；取消选择、或选了一个什么都不匹配的文件，只是让曲面片不带标签。该问题对每个（目录，结构）组合只问一次，因此多分子脚本不会重复提问。
- **裸 PDB 可视化按自身标量上色** — 并非每张图都以脚本形式提供：`A4` 以及 `A5`/`A6` 的逐轮文件是 PDB 文件，其 61-66 列*就是*那张可视化（逐原子力误差、轮间差异、单元编号）。当一个文件除坐标记录之外什么都没有（没有沉积结构头部），且该列非常量时，它就被识别为这类文件，并以同一渐变按其自身范围直接上色打开，范围显示在状态栏上。常量列（导出时写成 0.00 的那些）不上色，带头部记录的沉积 PDB 保持其元素颜色。只要存在这类数据，工具栏的 **Color by Beta** 按钮就会出现，在自动情形下初始即为启用，并可再次点击关闭。
- **查看器保留自己的背景** — 脚本的 `color Display Background white` 是参考渲染器的设置，这里有意不套用：白色画布会与编辑器主题冲突，并把白色原子冲掉（而 `Name` 配色正好把氢着成白色）。该指令仍会被解析并出现在场景模型中；`color scale method` 与 `color change rgb` 覆盖则会被套用。

脚本也会参与**结构对比**：比较一侧就是脚本绘制的那个分子，按脚本的方式上色。

以上全部位于同一个模块 `src/parsers/vmdVisualization.ts`——参考渲染器的色板与半径、把解析后的脚本解析成绘制属性的过程、脚本解析器本身，以及裸标量文件的识别——因此「可视化文件」这件事只有一个归属（解析器索引导出它的 API）。

#### 支持的文件格式

| 格式 | 扩展名 | 说明 |
|--------|-----------|-------|
| Gaussian Input | `.gjf`、`.gjf03`、`.gjf09`、`.gjf16`、`.com` | 读取 Link 0、route、标题、电荷/多重度、坐标、connect 段；支持固定原子写法（`C -1 x y z` 或 `C x y z -1`） |
| XYZ | `.xyz` | 标准 XYZ 格式，带原子数头行。首尾相接存放多帧的文件（ORCA 的 `<job>_trj.xyz`，以及沿用其布局的各类轨迹写出工具）会作为轨迹打开，可用帧导航；当每帧的注释行都写明能量时，还给出收敛曲线；被作业中断截断的那一帧会被丢弃，而不是显示成残缺的分子 |
| MOL2 | `.mol2` | Tripos MOL2 格式；读取 `@<TRIPOS>ATOM` 与 `@<TRIPOS>BOND` 段，支持键级（芳香 `ar` → 1.5） |
| Gaussian LOG | `.log` | 读取 `Standard orientation:` / `Input orientation:` 块；跳过把原子序数写成元素符号的行（例如来自 ONIOM 外部程序调用的输出）；支持多帧优化轨迹；提取 SCF 能量与收敛判据（力/位移）供收敛面板使用（仅在路由卡段要求优化时才显示——`force`/`freq` 单点路由会打印同样的表但不显示面板）；解析 `Harmonic frequencies` 块用于简正模式振动播放。xtb 轨迹（带 `xtb:` 注释标记的扩展多帧 XYZ）会被自动识别并按优化轨迹解析，给出能量 + 梯度范数收敛信息 |
| ORCA Input | `.inp` | 读取 `* xyz CHARGE MULT ... *` 块与 `%coords` 块；支持 xyz 与 xyzfile 坐标格式；解析 `units bohr/angs` 以换算坐标 |
| ORCA Output | `.out` | 读取 `CARTESIAN COORDINATES (ANGSTROEM)` 块，若该轮用 `%output` 抑制了它则回退到 `(A.U.)` 块；支持多帧优化轨迹（每个循环一帧，外加末尾的最终能量评估，标记为 `Cycle 1` … `Final`）；从回显输入与 SCF 设置中提取电荷与多重度；读取 `\|Geometry convergence\|` 表、该次运行自带的收敛判据、`VIBRATIONAL FREQUENCIES` + `NORMAL MODES` + `IR SPECTRUM`，以及 `!` 简单输入行供 Route 面板显示。格式按内容判定，因此存成 `.log` 的 ORCA 输出（或存成 `.out` 的 Gaussian 日志）同样能正确打开 |
| Turbomole Coord | `.coord` | 读取 `$coord` 段（Bohr → Å 换算）、`$chrg` 与 `$spin`/`$mult` 给出电荷与多重度；同时支持 `$chrg 1` 与 `$chrg=1` 语法 |
| PDB | `.pdb`、`.ent` | 以固定列方式读取 ATOM/HETATM 记录；元素取自 77-78 列或原子名；占据度取自 55-60 列，温度因子列（61-66）作为逐原子标量，当文件是裸坐标文件且该列非常量时会自动上色；CONECT 记录给出显式键（每行最多 8 个连接）；处理重复序号；ATOMIC_NUMBER_MAP 扩展到 Cf (98) |
| MOPAC | `.mop`、`.mopac`、`.dat` | 读取使用 `ELEM x 1 y 1 z 1` 内坐标的 MOPAC 输入格式；支持原子序数或元素符号；自动识别 Å/Bohr 单位；提取 CHARGE 与 MS 关键字；回退到 `CARTESIAN COORDINATES` 输出块 |
| VMD TCL | `.tcl` | VMD 可视化脚本，按「显示配方」而非数据文件读取：每一个 `mol new`（相对/绝对/带引号/带花括号路径，任意受支持的结构格式）都成为一个可切换的分子，每一个 `mol representation` / `mol color` / `mol selection` / `mol material` / `mol addrep` 都成为一个表示——含 CPK、VDW、Licorice、Lines 与 Surf 样式及其数值参数（球缩放、键半径、探针半径——作为比例套用到查看器自身的球棍几何上），颜色模式 `Beta`、`ColorID` 与 `Name`/`Element`，`mol scaleminmax` 范围，`mol modstyle`/`modcolor`/`modselect`/`modmaterial`，`color scale method` 与 `color change rgb`；脚本的 `color Display Background` 会被解析但不套用（查看器保留自己的背景）。旁车元数据（结构旁的 `eggFile*` + `unitfile`）会添加一个蛋列表，用于切换曲面片 |
| CIF | `.cif` | 晶体学信息文件；读取晶胞长度/角度、空间群与对称操作、分数/笛卡尔原子位点；套用对称操作构建完整晶胞；支持分数占据并做部分着色。**三代字典的标签都读**——`_symmetry_equiv_pos_as_xyz`（CIF 1.1）、`_space_group_symop_operation_xyz`（CIF 1.2/核心字典）与 `_space_group_symop.operation_xyz`（CIF 2）是同一个标签的三种写法，空间群名 `_symmetry_space_group_name_H-M` / `_space_group_name_H-M_alt` 同理，CIF 2 文件里的点号写法 `_atom_site.fract_x` 与 `_atom_site_fract_x` 也是同一列——因此现代文件不再被当成 P1 读成只剩四分之一的原子。对称操作可以写成 loop、单条 item 或分号文本块，重复出现的操作只列一次。文件没有声明任何对称操作、或声明的操作无法解析时，按 P1 读取**并在信息面板里说明**，而不是悄悄把晶胞缩小。元素取自 `type_symbol` 并剥离氧化态（`Fe3+` → `Fe`、`O2-` → `O`、`Cu(II)` → `Cu`）；`?`/`.` 一律视为缺失（缺失的空间群不是名叫「?」的空间群，保存时也不会写回去）；晶胞不完整时返回空结构并给出警告，而不是造一个 1 Å 的立方胞。对称像只有在就是同一个点（0.1 Å 以内）时才合并，因此相距 0.3 Å 的两个位点仍是两个原子；近乎重合的合并会连同位点标签与距离一起报告。文件自带的不对称单元与展开后的晶胞一并保留 |
| VASP POSCAR | `.vasp`、`POSCAR`、`CONTCAR` | VASP 结构文件；处理缩放因子、晶格矢量、元素种类/数量、选择性动力学，以及 Direct（分数）与 Cartesian 两种坐标；稳健的元素类型识别（不会误拒 Cs/Cd/Ca/Cl 等） |
| Gaussian Cube | `.cube` | Gaussian 体积文件；读取体素矢量与原子位置；自动识别 Bohr/Å 单位；对损坏的体素矢量能通过原子坐标范围推断正交晶格，稳健处理 |
| VESTA | `.vesta` | VESTA 结构文件；读取 CELLP（晶胞参数）、STRUC（分数坐标 + 占据）、GROUP（空间群）与 SYMOP（对称操作）段；套用对称操作展开完整晶胞；支持无序位点的分数占据 |

### 安装

#### 从 VSIX 安装（推荐）

1. 从 [Releases](https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode/releases) 下载最新的 `.vsix` 文件
2. 在 VS Code / Trae 中按 `Ctrl+Shift+P`
3. 输入 `Extensions: Install from VSIX...`
4. 选择下载的 `.vsix` 文件
5. **远程 SSH 用户**：请务必在**远程服务器上**安装该扩展（弹出提示时选择「Install on Remote」）

#### 从源码构建

```bash
git clone https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode.git
cd Molecular-Viewer-on-VScode
npm install
npm run compile
npx vsce package --no-dependencies
# 然后安装生成的 .vsix 文件
```

### 使用

#### 打开分子文件

1. 在资源管理器中**右键点击**受支持的文件 → **Molecular Viewer: Open 3D Viewer**
2. **命令面板**（`Ctrl+Shift+P`）→ `Molecular Viewer: Open 3D Viewer`
3. **自定义编辑器** — 双击受支持的文件并选择「Molecular 3D Viewer」

#### 导入结构（文件或 SMILES）

1. 点击工具栏中的 **Import Structure**（在 View 模式下可用；Diff 模式或晶体结构下不可用）
2. 一个对话框同时询问两件事：在输入框中键入**副本数量**（预填 `1`；1 到 50 的任意整数——不是 1 的取值会带着你的输入重新打开对话框，原因写在标题里），然后选择结构来源。输入框下方的两个条目在你键入时保持可见：
   - **From File** — 选择任意受支持的文件。多帧优化日志会打开帧选择器；晶体文件以分子形式贡献其晶胞原子
   - **From SMILES** — 键入一个 SMILES 字符串（例如 `CCO`、`c1ccccc1`、`CC(=O)Oc1ccccc1C(=O)O`、`OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O`）。输入框在每次按键时校验该字符串，并在合法之前持续显示解析器的报错信息；按 **Enter** 生成 3D 结构
3. 若结构大到「所请求副本数会使原子数超过 20 000」，会被裁剪到合适规模，并以警告说明新的数量；取消任何提示（Esc）都会取消整个导入
4. 每份副本都被添加在已有结构旁边且不与之重叠（包围球分离），副本之间以统一步距间隔，因此沿 +X 堆叠、彼此绝不重叠；整个导入是单步 **Undo**，**Save As** 会一起导出全部内容

SMILES 覆盖范围与限制——生成器读取 OpenSMILES 有机子集及常见扩展：分支、`%nn` 环闭合、芳香性、显式键符号、点号分隔的片段、同位素、电荷，以及 `@`/`@@` 四面体构型。平面芳香体系按平面输出：环上的取代基（含氢）位于环平面内，稠环——萘、菲、咖啡因的嘌呤母核——彼此共平面；芳香环与饱和环稠合时（茚满、四氢萘）同样把它的两条苄位键放到芳香平面内，只要刚性饱和环允许如此接近。四面体手性被遵守；未写 `/` 或 `\` 的双键默认取 **trans**；构象来自局部几何规则而非构象搜索，因此一条长的柔性链可能折叠得比较紧；同位素会被解析但不影响几何；有机金属几何不做特殊处理。

#### 设为默认查看器

在 `settings.json` 中添加：

```json
{
  "workbench.editorAssociations": {
    "*.gjf": "molecularViewer.editor",
    "*.gjf03": "molecularViewer.editor",
    "*.gjf09": "molecularViewer.editor",
    "*.gjf16": "molecularViewer.editor",
    "*.com": "molecularViewer.editor",
    "*.xyz": "molecularViewer.editor",
    "*.mol2": "molecularViewer.editor",
    "*.log": "molecularViewer.editor",
    "*.out": "molecularViewer.editor",
    "*.coord": "molecularViewer.editor",
    "*.inp": "molecularViewer.editor",
    "*.pdb": "molecularViewer.editor",
    "*.ent": "molecularViewer.editor",
    "*.mop": "molecularViewer.editor",
    "*.mopac": "molecularViewer.editor",
    "*.dat": "molecularViewer.editor",
    "*.tcl": "molecularViewer.editor",
    "*.cif": "molecularViewer.editor",
    "*.vasp": "molecularViewer.editor",
    "*.{POSCAR,CONTCAR,poscar,contcar}": "molecularViewer.editor",
    "*.cube": "molecularViewer.editor",
    "*.vesta": "molecularViewer.editor"
  }
}
```

该列表与查看器注册的扩展名一一对应；`.dat` 是 MOPAC 的备用输入/输出扩展名。不想从常规编辑器改走的行可以自行删掉。

#### 控制

| 操作 | 效果 |
|--------|--------|
| 左键拖拽 | 绕分子中心旋转分子 |
| 滚轮 | 放大 / 缩小 |
| 中键 / 右键拖拽 | 平移视图（光标锁定；Diff 模式下只平移你抓住的那一侧结构） |
| 悬停原子 | 显示元素 + 坐标提示 |
| 方向键 | 旋转分子（←/→ 绕 Y，↑/↓ 绕 X）；开启 Auto Rotate 后改为按住即加速自转——按得越久转得越快 |
| `Esc` | 只后退一层：用对话框自己的 **Cancel** 关掉它（预览过的几何会被还原）、放弃正在拖的框选矩形、从当前工具退回 View、或取消选择——再按一次只可能撤销更小的一件事 |
| `Ctrl+Z`（macOS 为 `Cmd+Z`） | 撤销，等同 Undo 按钮。对话框打开时、以及光标在输入框中时不会被接管（对话框本身压着一步 undo；输入框里那是浏览器的文本撤销） |
| `Ctrl+S`（macOS 为 `Cmd+S`） | 打开保存对话框，等同 **Save As**；对话框已经打开时再按一次，即按当前所选格式保存 |
| `Ctrl+A`（macOS 为 `Cmd+A`） | 选中屏幕上全部原子（Diff 模式下左右两侧都选）。页面其余部分都不是可选文本，所以这个键只可能是这个含义（输入框内仍是输入框自己的全选） |
| `Ctrl+C`（macOS 为 `Cmd+C`） | 把所选原子及**它们之间的键**复制到查看器自己的剪贴板 |
| `Ctrl+X`（macOS 为 `Cmd+X`） | 剪切：先复制，再走与 Delete 相同的删除路径，整体算一步 undo |
| `Ctrl+V`（macOS 为 `Cmd+V`） | 把剪贴板内容粘到结构旁边——与导入相同的 7 Å 包围球避让间距——并选中新到的原子，算一步 undo。晶体结构与 Diff 模式下不可用 |
| ◀ ▶（分子选择器） | 在 VMD TCL 脚本创建的分子（`mol new`）之间切换，仅当脚本定义了多于一个时显示 |
| Color by Beta | 为 PDB 61-66 列携带标量的结构切换标量着色——每个原子的球体以及与其相连的半键都取该原子的标量颜色，力误差图因此呈现为一条连续的渐变，而不是渐变球体之间夹着元素色的键；当此类文件在载入时被自动上色则初始为启用；当视图由脚本驱动时隐藏 |
| Auto Rotate 按钮 | 位于 View 模式下出现的小面板中：切换「甩动即自转」。开启后，按住左键旋转并松开，分子会以松开时的方向与速度继续自转（拖得越远/越快，转得越快）；下一次按住左键即停止，而中键/右键拖拽是平移、不会停 |
| Reset View 按钮 | 返回默认视图 |

这些快捷键由查看器自己处理——不向编辑器贡献任何 keybinding，因此不会占用 VS Code / Trae 的任何快捷键——只要 3D 视图有焦点就可用。输入框中键入的按键仍归输入框，Esc 除外。视图不是文档：除输入框之外，页面上没有任何可选文本（因此 Ctrl+A 永远是「全选结构」，拖拽划过面板也不会留下选中痕迹）——只有加载失败的报错文字仍可选，方便复制反馈。片段剪贴板是查看器自己的，不是系统剪贴板：复制不会顶掉你在别的程序里复制的文字，也从不读取它；副本会尽力镜像到面板的 localStorage，因此在一个面板里复制的片段可以在另一个面板里粘贴。

#### 编辑流程

1. 点击工具栏按钮进入某种编辑模式（例如 **Bond Length**）
2. 在 3D 视图中点击原子来选择它们（选中的原子发黄光）
3. 弹出模态对话框显示当前值
4. 用下拉框选择哪些原子**固定**、哪些**移动**
5. 用**数值输入**或**滑块**调节取值——3D 视图实时更新
6. 点击 **OK** 确认改动，或 **Cancel** 还原
7. 用 **Save As** 把修改后的结构导出为新文件

### 架构

```
┌─────────────────────┐          ┌──────────────────────────┐
│   Local (Windows)   │   SSH    │   Remote (Linux Server)  │
│                     │ ───────> │                          │
│  Trae IDE (UI)      │          │  Trae Server (Extension) │
│  ├─ Webview 3D      │ <─────── │  ├─ Parse .gjf/.xyz/.log │
│  ├─ Three.js (inline)│  data   │  ├─ Bond detection       │
│  ├─ Editing UI      │          │  ├─ Save file (VS Code)  │
│  └─ Mouse events    │          │  └─ Return molecule data  │
└─────────────────────┘          └──────────────────────────┘
```

扩展运行在**远程一侧**（读取文件、解析、保存），而 Webview 在**本地一侧**渲染（Three.js 内联进 HTML、鼠标交互、编辑界面）。

### 项目结构

```
molecular-viewer/
├── src/
│   ├── extension.ts           # 扩展入口
│   ├── types.ts               # TypeScript 类型定义
│   ├── parsers/               # 文件格式读取器：文本进，分子出
│   │   ├── index.ts           # 解析器分发（自动识别格式）
│   │   ├── gjfParser.ts       # Gaussian .gjf 解析器（connect 段、固定原子）
│   │   ├── xyzParser.ts       # XYZ 解析器；多帧 XYZ（轨迹）
│   │   ├── mol2Parser.ts      # Tripos MOL2 格式解析器
│   │   ├── logParser.ts       # Gaussian LOG 解析器（优化轨迹）
│   │   ├── xtbParser.ts       # xtb LOG 解析器（优化轨迹，按内容识别）
│   │   ├── coordParser.ts     # Turbomole .coord 解析器（Bohr → Å）
│   │   ├── orcaInpParser.ts   # ORCA 输入 .inp 解析器
│   │   ├── orcaOutParser.ts   # ORCA 输出 .out 解析器（优化轨迹、收敛、频率）
│   │   ├── pdbParser.ts       # PDB 格式解析器（ATOM/HETATM/CONECT）
│   │   ├── mopacParser.ts     # MOPAC 输入/输出解析器
│   │   ├── vmdVisualization.ts # 可视化文件：VMD TCL 脚本（分子、表示、选择）+ 色板/半径 + 表示 → 颜色、尺寸、曲面片 + 裸标量文件识别
│   │   ├── cifParser.ts       # CIF 晶体结构解析器（对称性、占据）
│   │   ├── vaspParser.ts      # VASP POSCAR/CONTCAR 解析器（缩放、Direct/Cartesian）
│   │   ├── cubeParser.ts      # Gaussian Cube 解析器（体素矢量、稳健晶格）
│   │   └── vestaParser.ts     # VESTA 结构解析器（CELLP、STRUC、SYMOP）
│   ├── tools/                 # 化学引擎：分子进，分子出（不做文件 I/O）
│   │   ├── bondPerception.ts  # 键检测与键级感知——元素对表、键级查找、价态求解器、芳香性；以 `BP` 序列化进 webview，使页面与宿主运行同一份模块
│   │   ├── bondDetector.ts    # bondPerception 之上的薄入口（`ensureBonds`），宿主唯一的调用点
│   │   ├── smilesTo3d.ts      # SMILES → 3D 坐标（零依赖，test/smiles_to_3d.py 的移植）
│   │   ├── addHydrogens.ts    # 补氢：依据引擎的表做经验价补齐，骨架冻结
│   │   └── mmOptimizer.ts     # MM Optimize：内置力场（键/角/扭转/异常/LJ + 轻量静电/氢键/金属配位/极化）、第 4/5 周期元素表、可设截断、low/medium/high 预设，单调优化器以可分块的会话形式提供（pause/resume/history）
│   └── webview/
│       └── molecularViewer.ts # 自定义编辑器 + Three.js webview + 编辑
├── dist/                      # 编译后的 JavaScript（预构建）
├── test/                      # 示例分子文件
├── media/
│   └── three.min.js           # Three.js r128（本地打包）
├── package.json
├── tsconfig.json
└── LICENSE
```

### 开发

```bash
# 安装依赖
npm install

# 编译 TypeScript
npm run compile

# 监听改动
npm run watch

# 只做类型检查、不产出
npm run lint

# 打包为 .vsix
npm run package
```

### 贡献

1. Fork 本仓库
2. 创建特性分支（`git checkout -b feature/amazing-feature`）
3. 提交改动（`git commit -m 'Add amazing feature'`）
4. 推送分支（`git push origin feature/amazing-feature`）
5. 发起 Pull Request

### 更新日志

完整发布历史见 [CHANGELOG.md](CHANGELOG.md)。

### 许可证

[MIT](LICENSE)
