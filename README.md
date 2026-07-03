# Molecular Viewer

A VS Code / Trae extension for visualizing and editing molecular structures in 3D, designed for computational chemists working with Gaussian, ORCA, and other quantum chemistry software — especially on remote servers where GUI tools like GaussView are unavailable.

## Features

- **3D Ball-and-Stick Rendering** — Atoms rendered as spheres (scaled by covalent radius) with CPK coloring; bonds rendered as dual-colored cylinders
- **Bond Order Support** — Visual distinction for single (1 line), aromatic (1 solid + 1 dashed), double (2 lines), and triple (3 lines) bonds
- **Auto Bond Detection** — When files lack explicit connectivity, bonds are automatically detected using element-pair-specific bond length specifications (C-C, C-N, C-O, etc.) with tolerances; pair-specific distance cutoffs; valence-based bond order refinement (max valence constraints for H/C/N/O/F/S/P/Cl/Br/I/B); post-processing fixes for N and C-O bonds; fallback to covalent radii ratio for unspecified pairs
- **GJF Connect Section** — Reads explicit bond information from coordinate section in GJF files, including bond orders (1.0, 1.5, 2.0, 3.0)
- **Molecular Info Display** — Shows the Hill-system molecular formula, atom count, charge, electron count, and spin multiplicity in the top-left corner of the 3D view; for crystal files, an additional Unit Cell composition line shows the full unit-cell formula
- **Interactive Mouse Control**:
  - Left drag → Rotate around molecule center
  - Scroll → Zoom in/out
  - Middle/Right drag → Pan
  - Hover atom → Show element name + coordinates
  - Arrow keys (↑/↓/←/→) → Rotate molecule (left/right around Y axis, up/down around X axis; ignored when typing in input fields)
  - Delete key → If atoms are selected (via Select Atoms / Box Select), opens batch delete confirmation; otherwise enters Delete Atom mode
- **Touch Support** — Single-finger rotate, pinch-to-zoom
- **Remote-SSH Compatible** — Works seamlessly when editing files on remote Linux servers via VS Code/Trae Remote-SSH
- **Performance** — Dirty-flag rendering: the scene only re-renders on interaction or structural change (not every frame), so idle GPU usage is eliminated; meshes are properly disposed (geometry, material, texture) on rebuild to prevent GPU memory leaks that previously caused progressive lag within a session; adaptive level-of-detail reduces sphere segment count for large structures (>200 atoms: 20×16, >1000 atoms: 12×10) to maintain interactive framerate; the axes indicator SVG is rebuilt at most once per animation frame (~60 Hz) instead of on every mousemove event

### Molecular Editing

- **Bond Length Adjustment** — Select 2 atoms, view current bond length, choose which atom to fix, choose Translate Group (rigidly move atom + connected subtree, preserving relative positions) or Translate Atom (move only the selected atom), adjust via numeric input or slider with real-time 3D preview
- **Bond Angle Adjustment** — Select 3 atoms (2nd is the vertex), view current angle, fix/move either side, Translate Group/Atom option, real-time preview
- **Dihedral Angle Adjustment** — Select 4 atoms, view current dihedral, fix/move either side, Translate Group/Atom option, real-time preview
- **Bond Order Editing** — Change bond order (none / single / aromatic 1.5 / double / triple) in the Bond Length modal; selecting "None (0)" removes the bond; changes reflected immediately in 3D display
- **Add Atom** — Click anchor atom, choose element (70+ elements), set bond length and bond order, direction auto-calculated for polyhedral coordination with existing bonds (maximum angular separation from existing neighbors, naturally producing linear/trigonal/tetrahedral/octahedral geometry); repeated additions to the same anchor fan out instead of overlapping
- **Delete Atom** — Click atom and confirm; atoms and bonds are automatically re-indexed; or delete all currently highlighted atoms at once by selecting them first via Select Atoms / Box Select, then clicking Delete
- **Replace Atom** — Click the Replace Atom toolbar button, then click one or more atoms in the 3D view to select them (click again to deselect); click the Replace Atom button again to open the element picker (70+ elements), then confirm to replace all selected atoms with the chosen element while preserving their coordinates; works with crystal base atoms too
- **Select Atoms** — Click the Select Atoms toolbar button to show a floating panel (bottom-left corner) where you can type atom indices (1-based), ranges (e.g. `3-10`), or element symbols (e.g. `C H`) to highlight specific atoms in yellow; press Enter or click Select to apply, Clear to reset; you can also directly click atoms in the 3D view to select them; clicking a highlighted atom again deselects it (toggle); selections accumulate while in Select Atoms mode; the panel is a truly floating overlay that does not affect the molecular view size; the panel toggles with the toolbar button and hides automatically when switching to another mode; selections are preserved when switching between Select Atoms and Box Select
- **Box Select** — Independent toolbar button; in this mode, hold the left mouse button and drag a rectangle, then release to select all atoms inside the box (atoms are projected to screen coordinates and tested against the rectangle); selections accumulate with the same yellow highlighting as Select Atoms, and are preserved when switching between Box Select and Select Atoms; switching to any other mode clears the selection
- **Undo** — Toolbar button reverts the most recent edit/delete operation; snapshots are captured automatically when an edit modal opens (including Add Atom) and before delete/remove-disorder/batch-delete execute; the stack is capped at 50 entries and resets on file reload
- **Save As** — Export modified structure in 12 formats: XYZ, Gaussian GJF (preserving original file structure), CIF (fractional coordinates), VASP POSCAR (Direct fractional coordinates), Gaussian Cube (with volumetric grid), VESTA (P1 symmetry, fractional coordinates), Turbomole Coord (Å→Bohr conversion), ORCA Input, MOL2 (with bond orders), MDL Mol, PDB (with CONECT records), or MOPAC Input; GJF output preserves original Link 0, route, title, charge/mult, and post-connect content; GJF atom coordinates use GaussView-style 8-decimal aligned format; connect section includes all atom lines; CIF/VASP/Cube/VESTA formats available only when crystal data is present
- **Continuous Editing** — After completing an edit, the viewer stays in the current editing mode for repeated adjustments
- **Cancel/Undo** — Cancel button restores original coordinates before confirming edits

### Optimization Trajectory Navigation (LOG files)

- **Frame Stepping** — ◀ Prev / Next ▶ buttons to step through optimization frames
- **Jump to Frame** — Direct input field to jump to a specific frame number
- **Auto Play** — Automatically cycle through all frames with 500ms interval

### Structure Diff

- **Diff Button** — Compare the current structure against another molecular file (any supported format, no need for matching formats)
- **Frame Selection** — If the comparison file is an optimization LOG/OUT, a QuickPick lets you choose which frame to compare
- **Skeleton Check** — Graph isomorphism via recursive backtracking atom matching (element + degree + neighbor signature with most-constrained-variable heuristic); uses the plugin's own bond detection algorithm (not the file's original bond orders) to ensure consistent skeleton comparison even when files specify different bond orders for the same connectivity; matching considers **connectivity only** (bonded or not), not bond order — a single bond in one molecule and a double bond in the other at the same position does not prevent skeleton matching
- **Side-by-Side View** — Left viewport shows the original molecule, right shows the comparison; each side has **independent** rotation, pan, and zoom controls
- **Correct Aspect Ratio** — Each viewport uses its own camera aspect ratio matching the half-width, preventing horizontal distortion
- **Conformation Diff** — When skeletons match, highlights bond length, bond angle, and dihedral differences in orange on both sides; full list of all differences displayed (no truncation); each difference shows the percentage change relative to the average value
- **Adjustable Threshold** — A percentage slider (0–20%) in the diff results panel lets you interactively adjust the sensitivity for what counts as a "difference"; all three types use the same relative difference formula (Δ/avg × 100%); for dihedrals near 0° where the average is very small, a fallback ensures robustness; the 3D highlights and difference list update in real-time as you drag the slider
- **Diff Info Panel** — Floating panel listing all differences sorted by magnitude, with atom labels and values from both structures; can be closed and reopened via "Show Results" button; diff results are always accessible and not overwritten by measurement info
- **Auto-Mirror Selection** — When selecting atoms for bond length/angle/dihedral measurement on one side, the corresponding atoms are automatically selected on the other side via the atom mapping, and both values plus the Δ difference are displayed simultaneously

### Crystal Structure Viewing

- **Unit Cell Display** — When opening CIF, VASP, Cube, or VESTA files, the crystal unit cell is rendered with a wireframe outline showing the lattice boundaries
- **Supercell Boundary Controls** — 6 input boxes (a/b/c min/max) let you define a custom region of the crystal to display; supports fractional values (step 0.1) for partial cell selection; the view is a cutout of the crystal, not a deformation of the lattice; non-integer boundaries (e.g. -0.4 to 0.4) do not generate mirror images — atoms outside the range are simply not shown; finite molecules (not extending periodically) straddling the boundary are shown completely, extending beyond the boundary to avoid broken bonds; infinite periodic chains are truncated at the boundary via split bonds
- **Local Crystal Modifications** — Atom edits (bond length, angle, dihedral) in crystal structures modify only the specific atom being edited, not all periodic images of the same base atom; each periodic image is independently editable
- **Minimum Image Convention Bonds** — Crystal bonds are detected using the minimum image convention (MIC): interatomic distances are computed at their shortest periodic distance, and bonds crossing cell boundaries are rendered as short split bonds extending outside the view (not as long bonds spanning across the view)
- **Rotating Axes Indicator** — A small axes gizmo in the bottom-left corner rotates synchronously with the main 3D view, providing spatial orientation reference
- **Fractional Occupancy Rendering** — Atoms with partial occupancy (disordered crystal sites) are rendered with partial coloration proportional to their occupancy value
- **Remove Disorder** — One-click button to clean up disordered sites: occupancy < 0.5 removed, > 0.5 set to 1, = 0.5 kept
- **Crystal Save** — Export modified crystal structures as CIF (fractional coordinates, minimal format), VASP POSCAR (Direct fractional coordinates grouped by element), Gaussian Cube (with minimal volumetric grid), or VESTA (P1 symmetry, fractional coordinates)

### Supported File Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Gaussian Input | `.gjf`, `.gjf03`, `.gjf09`, `.gjf16`, `.com` | Reads Link 0, route, title, charge/mult, coordinates, connect section; supports fixed atom notation (`C -1 x y z` or `C x y z -1`) |
| XYZ | `.xyz` | Standard XYZ format with atom count header |
| MOL2 | `.mol2` | Tripos MOL2 format; reads `@<TRIPOS>ATOM` and `@<TRIPOS>BOND` sections with bond order support (aromatic `ar` → 1.5) |
| Gaussian LOG | `.log` | Reads `Standard orientation:` / `Input orientation:` blocks; skips lines where atomic numbers are written as element symbols (e.g. from ONIOM external program calls); supports multi-frame optimization trajectory |
| ORCA Input | `.inp` | Reads `* xyz CHARGE MULT ... *` blocks and `%coords` blocks; supports xyz and xyzfile coordinate formats |
| ORCA Output | `.out` | Reads `CARTESIAN COORDINATES (ANGSTROEM)` blocks; supports multi-frame optimization trajectory; extracts charge and multiplicity |
| Turbomole Coord | `.coord` | Reads `$coord` section (Bohr → Å conversion), `$chrg` and `$spin`/`$mult` for charge and multiplicity |
| PDB | `.pdb`, `.ent` | Reads ATOM/HETATM records with fixed-column parsing; element from columns 77-78 or atom name; CONECT records for explicit bonds; handles duplicate serials |
| MOPAC | `.mop`, `.mopac`, `.dat` | Reads MOPAC input format with `ELEM x 1 y 1 z 1` internal coordinates; supports atomic numbers or element symbols; auto-detects Å/Bohr units; extracts CHARGE and MS keywords; falls back to `CARTESIAN COORDINATES` output blocks |
| VMD TCL | `.tcl` | VMD visualization script; reads `mol new <file>` to load the referenced molecular file; parses `mol color ColorID N` + `mol selection "index ..."` to assign per-group atom colors using VMD's 33-color palette; supports relative and absolute file paths |
| CIF | `.cif` | Crystallographic Information File; reads cell lengths/angles, space group symmetry operators, and fractional/Cartesian atom sites; applies symmetry operations to build the full unit cell; supports fractional occupancy with partial coloration |
| VASP POSCAR | `.vasp`, `POSCAR`, `CONTCAR` | VASP structure file; handles scale factor, lattice vectors, element types/counts, selective dynamics, and both Direct (fractional) and Cartesian coordinates |
| Gaussian Cube | `.cube` | Gaussian volumetric file; reads voxel vectors and atomic positions; auto-detects Bohr/Angstrom units; robustly handles corrupted voxel vectors by inferring an orthogonal lattice from atomic coordinate ranges |
| VESTA | `.vesta` | VESTA structure file; reads CELLP (cell parameters), STRUC (fractional coordinates + occupancy), GROUP (space group), and SYMOP (symmetry operations) sections; applies symmetry operations to expand the full unit cell; supports fractional occupancy for disordered sites |
| MDL Mol | `.mol` | Basic support |
| SDF | `.sdf` | Basic support |

## Installation

### From VSIX (Recommended)

1. Download the latest `.vsix` file from [Releases](https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode/releases)
2. In VS Code / Trae, press `Ctrl+Shift+P`
3. Type `Extensions: Install from VSIX...`
4. Select the downloaded `.vsix` file
5. **For Remote-SSH**: Make sure to install the extension **on the remote server** (choose "Install on Remote" when prompted)

### From Source

```bash
git clone https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode.git
cd molecular-viewer
npm install
npm run compile
npx vsce package --no-dependencies
# Then install the generated .vsix file
```

## Usage

### Opening Molecular Files

1. **Right-click** a supported file in the Explorer → **Molecular Viewer: Open 3D Viewer**
2. **Command Palette** (`Ctrl+Shift+P`) → `Molecular Viewer: Open 3D Viewer`
3. **Custom Editor** — Double-click a supported file and select "Molecular 3D Viewer"

### Set as Default Viewer

Add to your `settings.json`:

```json
{
  "workbench.editorAssociations": {
    "*.gjf": "molecularViewer.editor",
    "*.xyz": "molecularViewer.editor",
    "*.com": "molecularViewer.editor",
    "*.mol2": "molecularViewer.editor",
    "*.log": "molecularViewer.editor",
    "*.out": "molecularViewer.editor",
    "*.coord": "molecularViewer.editor",
    "*.inp": "molecularViewer.editor",
    "*.pdb": "molecularViewer.editor",
    "*.ent": "molecularViewer.editor",
    "*.mop": "molecularViewer.editor",
    "*.tcl": "molecularViewer.editor",
    "*.cif": "molecularViewer.editor",
    "*.vasp": "molecularViewer.editor",
    "*.cube": "molecularViewer.editor"
  }
}
```

### Controls

| Action | Effect |
|--------|--------|
| Left mouse drag | Rotate molecule around its center |
| Mouse scroll | Zoom in / out |
| Middle / Right mouse drag | Pan view |
| Hover over atom | Show element + coordinates tooltip |
| Arrow keys | Rotate molecule (←/→ around Y, ↑/↓ around X) |
| Reset View button | Return to default view |

### Editing Workflow

1. Click a toolbar button to enter an editing mode (e.g., **Bond Length**)
2. Click atoms in the 3D view to select them (selected atoms glow yellow)
3. A modal dialog appears showing the current value
4. Choose which atoms to **fix** vs. **move** using the dropdown
5. Adjust the value using the **numeric input** or **slider** — the 3D view updates in real-time
6. Click **OK** to confirm the change, or **Cancel** to revert
7. Use **Save As** to export the modified structure to a new file

## Architecture

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

## Project Structure

```
molecular-viewer/
├── src/
│   ├── extension.ts           # Extension entry point
│   ├── types.ts               # TypeScript type definitions
│   ├── parsers/
│   │   ├── index.ts           # Parser dispatcher (auto-detect format)
│   │   ├── gjfParser.ts       # Gaussian .gjf parser (connect section, fixed atoms)
│   │   ├── xyzParser.ts       # XYZ format parser
│   │   ├── mol2Parser.ts      # Tripos MOL2 format parser
│   │   ├── logParser.ts       # Gaussian LOG parser (optimization trajectory)
│   │   ├── coordParser.ts     # Turbomole .coord parser (Bohr → Å)
│   │   ├── orcaInpParser.ts   # ORCA input .inp parser
│   │   ├── orcaOutParser.ts   # ORCA output .out parser (optimization trajectory)
│   │   ├── pdbParser.ts       # PDB format parser (ATOM/HETATM/CONECT)
│   │   ├── mopacParser.ts     # MOPAC input/output parser
│   │   ├── tclParser.ts       # VMD TCL script parser (color groups)
│   │   ├── cifParser.ts       # CIF crystal structure parser (symmetry, occupancy)
│   │   ├── vaspParser.ts      # VASP POSCAR/CONTCAR parser (scale, Direct/Cartesian)
│   │   ├── cubeParser.ts      # Gaussian Cube parser (voxel vectors, robust lattice)
│   │   ├── vestaParser.ts     # VESTA structure parser (CELLP, STRUC, SYMOP)
│   │   └── bondDetector.ts    # Covalent radii bond detection + order estimation
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

## Development

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

## Roadmap

- [ ] MOL/SDF full parser with explicit bond info
- [ ] Multiple display styles (wireframe, space-filling, licorice)
- [ ] Vibration animation from frequency calculations
- [ ] Export as PNG/SVG
- [ ] Undo/redo history for edits

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

[MIT](LICENSE)
