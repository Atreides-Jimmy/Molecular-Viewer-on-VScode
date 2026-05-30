# Molecular Viewer

A VS Code / Trae extension for visualizing and editing molecular structures in 3D, designed for computational chemists working with Gaussian, ORCA, and other quantum chemistry software — especially on remote servers where GUI tools like GaussView are unavailable.

## Features

- **3D Ball-and-Stick Rendering** — Atoms rendered as spheres (scaled by covalent radius) with CPK coloring; bonds rendered as dual-colored cylinders
- **Bond Order Support** — Visual distinction for single (1 line), aromatic (1 solid + 1 dashed), double (2 lines), and triple (3 lines) bonds
- **Auto Bond Detection** — When files lack explicit connectivity, bonds are automatically detected using element-pair-specific bond length specifications (C-C, C-N, C-O, etc.) with tolerances; pair-specific distance cutoffs; valence-based bond order refinement (max valence constraints for H/C/N/O/F/S/P/Cl/Br/I/B); post-processing fixes for N and C-O bonds; fallback to covalent radii ratio for unspecified pairs
- **GJF Connect Section** — Reads explicit bond information from coordinate section in GJF files, including bond orders (1.0, 1.5, 2.0, 3.0)
- **Molecular Info Display** — Shows atom count, charge, electron count, and spin multiplicity in the top-left corner of the 3D view
- **Interactive Mouse Control**:
  - Left drag → Rotate around molecule center
  - Scroll → Zoom in/out
  - Middle/Right drag → Pan
  - Hover atom → Show element name + coordinates
- **Touch Support** — Single-finger rotate, pinch-to-zoom
- **Remote-SSH Compatible** — Works seamlessly when editing files on remote Linux servers via VS Code/Trae Remote-SSH

### Molecular Editing

- **Bond Length Adjustment** — Select 2 atoms, view current bond length, choose which atom to fix, adjust via numeric input or slider with real-time 3D preview
- **Bond Angle Adjustment** — Select 3 atoms (2nd is the vertex), view current angle, fix/move either side, real-time preview
- **Dihedral Angle Adjustment** — Select 4 atoms, view current dihedral, fix/move either side, real-time preview
- **Bond Order Editing** — Change bond order (none / single / aromatic 1.5 / double / triple) in the Bond Length modal; selecting "None (0)" removes the bond; changes reflected immediately in 3D display
- **Add Atom** — Click anchor atom, choose element (70+ elements), set bond length and bond order, direction auto-calculated from existing bonds
- **Delete Atom** — Click atom and confirm; atoms and bonds are automatically re-indexed
- **Select Atoms** — Input atom indices (1-based), ranges (e.g. `3-10`), or element symbols (e.g. `C H`) to highlight specific atoms in yellow
- **Save As** — Export modified structure in 8 formats: XYZ, Gaussian GJF (preserving original file structure), Turbomole Coord (Å→Bohr conversion), ORCA Input, MOL2 (with bond orders), MDL Mol, PDB (with CONECT records), or MOPAC Input; GJF output preserves original Link 0, route, title, charge/mult, and post-connect content; connect section includes all atom lines
- **Continuous Editing** — After completing an edit, the viewer stays in the current editing mode for repeated adjustments
- **Cancel/Undo** — Cancel button restores original coordinates before confirming edits

### Optimization Trajectory Navigation (LOG files)

- **Frame Stepping** — ◀ Prev / Next ▶ buttons to step through optimization frames
- **Jump to Frame** — Direct input field to jump to a specific frame number
- **Auto Play** — Automatically cycle through all frames with 500ms interval

### Supported File Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Gaussian Input | `.gjf`, `.gjf03`, `.gjf09`, `.gjf16`, `.com` | Reads Link 0, route, title, charge/mult, coordinates, connect section; supports fixed atom notation (`C -1 x y z` or `C x y z -1`) |
| XYZ | `.xyz` | Standard XYZ format with atom count header |
| MOL2 | `.mol2` | Tripos MOL2 format; reads `@<TRIPOS>ATOM` and `@<TRIPOS>BOND` sections with bond order support (aromatic `ar` → 1.5) |
| Gaussian LOG | `.log` | Reads `Standard orientation:` / `Input orientation:` blocks; supports multi-frame optimization trajectory |
| ORCA Input | `.inp` | Reads `* xyz CHARGE MULT ... *` blocks and `%coords` blocks; supports xyz and xyzfile coordinate formats |
| ORCA Output | `.out` | Reads `CARTESIAN COORDINATES (ANGSTROEM)` blocks; supports multi-frame optimization trajectory; extracts charge and multiplicity |
| Turbomole Coord | `.coord` | Reads `$coord` section (Bohr → Å conversion), `$chrg` and `$spin`/`$mult` for charge and multiplicity |
| PDB | `.pdb`, `.ent` | Reads ATOM/HETATM records with fixed-column parsing; element from columns 77-78 or atom name; CONECT records for explicit bonds; handles duplicate serials |
| MOPAC | `.mop`, `.mopac`, `.dat` | Reads MOPAC input format with `ELEM x 1 y 1 z 1` internal coordinates; supports atomic numbers or element symbols; auto-detects Å/Bohr units; extracts CHARGE and MS keywords; falls back to `CARTESIAN COORDINATES` output blocks |
| VMD TCL | `.tcl` | VMD visualization script; reads `mol new <file>` to load the referenced molecular file; parses `mol color ColorID N` + `mol selection "index ..."` to assign per-group atom colors using VMD's 33-color palette; supports relative and absolute file paths |
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
    "*.tcl": "molecularViewer.editor"
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
- [ ] CIF crystal structure support
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
