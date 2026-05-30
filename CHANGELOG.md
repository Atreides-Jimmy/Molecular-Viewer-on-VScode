# Change Log

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
