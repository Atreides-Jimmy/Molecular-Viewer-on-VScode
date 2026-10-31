# Molecular Viewer 3D

**English** | [简体中文](#简体中文)

**Jump to**

- **Features** — [overview](#features) · [molecular editing](#molecular-editing) · [trajectory and vibrations](#trajectory-and-vibrations) · [structure diff](#structure-diff) · [crystal structures](#crystal-structures) · [VMD TCL scripts](#vmd-tcl-scripts) · [supported formats](#supported-file-formats)
- **Getting started** — [installation](#installation) · [usage](#usage) · [opening files](#opening-files) · [importing a structure](#importing-a-structure) · [controls](#controls)
- **More** — [development](#development) · [contributing](#contributing) · [changelog](#changelog) · [license](#license)

## English

A VS Code / Trae extension that brings a GaussView-style 3D molecular viewer and editor into your IDE. It opens Gaussian, ORCA, MOPAC, xtb and VMD output alongside common structure formats, renders and edits them interactively, follows optimisation trajectories, compares structures and inspects crystals — including files on remote servers over Remote-SSH.

### Features

- **Four display styles** — **Ball & Stick** (default), **Space-filling**, **Stick** and **Polyhedra**. Polyhedra draws a translucent coordination hull around every metal or metalloid centre with at least four non-centre neighbours and completes the shell across crystal faces; centres with too few neighbours, a planar shell or an all-metal shell stay ball-and-stick. Polyhedra is available for crystal pages only and cannot be used in Diff mode. Choose the style under **Display** in the **View** card: it belongs to the view, so edits, frame changes, undo and optimisation keep it. Ghost atoms stay small markers and every atom remains clickable in every style.
- **Bond orders** — Single, aromatic (one solid + one dashed line), double and triple bonds are drawn differently. Aromatic dashes are planned per aromatic system, so they stay inside a ring or form one continuous band along a chain.
- **Automatic bond detection** — When a file carries no connectivity, bonds and bond orders are inferred from coordinates using element-pair bond lengths and a valence/charge check, so groups such as nitro, ammonium and carboxylate are drawn correctly and formal charges are derived. Atoms outside their usual valence are reported. Bonds stated by GJF connect sections, MOL2 files and PDB CONECT records are preserved, and connectivity is re-derived after atoms are moved.
- **Molecular info panel** — The top-left corner shows the Hill formula, atom count, charge, electron count and spin multiplicity; crystal files add the full unit-cell composition. If the parser had to work around something in the file, an amber line names the issue instead of hiding it.
- **Non-covalent interactions** — The **Interactions** mode lists H-bonds (H···A < 2.5 Å, D–H···A > 120°), halogen/chalcogen/pnictogen contacts, π–π stacking (parallel, parallel-displaced or T-shaped) and CH–π contacts, nearest first. Each type has its own checkbox, count and dashed guide; click a row to select both partners. Crystal searches include periodic images. Escape or View closes the mode.
- **Two-row toolbar** — Editing tools occupy the first row; contact tools, selection tools and view/action controls the second. Each row scrolls horizontally on its own, so every button is reachable even in a narrow panel.
- **Collapsible cards** — Every floating card can be hidden with **×** and leaves a named chip at its screen edge; clicking the chip restores the card at the same place and size. The layout narrows and stacks cards so they never cover one another or leave the viewport.
- **Mouse and touch** — Drag to rotate, scroll to zoom, middle/right-drag to pan with the grabbed point locked under the cursor, and hover an atom for its element and coordinates. Arrow keys rotate, Delete deletes the selection (or enters Delete Atom mode), and touch supports one-finger rotation and pinch zoom. See [Controls](#controls) for the full list.
- **Auto Rotate** — In View mode, turn on **Auto Rotate** in the bottom-right panel to throw the view. A rotating drag keeps spinning in the direction and at the speed it was released with (only the last 120 ms count), and a held arrow key accelerates the spin. The next left press stops it; a right/middle press pans without stopping. It turns the camera only, never the coordinates.
- **Alt-drag edits the coordinates** — In View mode, hold **Alt** while dragging: left-drag rotates the structure and right/middle-drag translates it. The atoms move rather than the camera, so **Save As** writes the orientation you dragged to. One undo step per drag, and no bonds are re-derived because a rigid transform changes no distance. In Diff mode and on crystals the drag stays view-only.
- **Orientation axes** — A gizmo at bottom-left shows Cartesian X/Y/Z axes for molecules and lattice a/b/c axes for crystals, following every rotation. **Reset View** restores the initial orientation.
- **Remote-SSH** — Parsing and saving run on the remote host while the 3D view runs locally, so the extension works unchanged on a remote Linux server.
- **VMD TCL scripts** — Opening a VMD script reproduces its representations, colours, selections and surfaces, including multiple molecules and optional egg metadata. See [VMD TCL scripts](#vmd-tcl-scripts).
- **Performance with large structures** — The scene redraws only on interaction or change, and GPU resources are released on rebuild. Detail adapts to structure size, and spheres and bonds are drawn with instancing, so structures with thousands of atoms stay responsive. Move Atoms rebuilds only the bonds touching the selection, and long logs are parsed in a single pass.

#### Molecular Editing

- **Bond length** — Select two atoms, fix one side, choose **Translate Group** (move the atom and its connected subtree) or **Translate Atom**, then drag the slider or type a value; the 3D view previews in real time. The same dialog also edits bond order (none / single / aromatic / double / triple — choosing "None" deletes the bond).
- **Bond angle** — Select three atoms (the second is the vertex), fix or move either side, and adjust with the slider or numeric input; preview is live.
- **Dihedral** — Select four atoms and adjust the dihedral the same way.
- **Add atom** — Click an anchor atom, pick one of 70+ elements, set the bond length and order. The direction maximizes separation from existing neighbours, producing linear, trigonal, tetrahedral or octahedral geometry; repeated additions fan out.
- **Delete atom** — Click an atom and confirm; indices are updated automatically. To delete several atoms, select them first and press Delete.
- **Replace atom** — Click **Replace Atom**, click one or more atoms (click again to deselect), then click the button again to pick a new element from 70+. Coordinates are preserved, and crystal base atoms can be replaced too.
- **Ghost atoms** — Mark a position with no nucleus, electrons or bonds — the counterpoise (BSSE) marker. Add one with the Add Atom dialog's Ghost switch or the `X (dummy point)` element; the Replace Atom dialog's Ghost switch turns a whole selection into ghosts of their own elements. Ghosts are excluded from bonding, formula, electron count, Add H and MM, are drawn as faint markers, and round-trip through every format that can express them (Gaussian `C-Bq`, ORCA `C:`, XYZ/PDB dummy markers and others). Crystal structures cannot contain ghosts.
- **Add H** — **Add H** fills the selected atoms' missing valences with hydrogens, using element-specific bond lengths, angles and the centre's hybridization. Existing atoms never move, re-running adds nothing, hydrogens that would clash are skipped and reported, and the whole addition is one undo step. Not available for crystals or in Diff.
- **MM Optimize** — The **MM Optimize** panel runs the built-in, dependency-free force field. Set every parameter yourself or use **Low / Medium / High** presets, choose the atoms to relax and start. Only the selected atoms move, but fixed atoms still exert bonded and non-bonded forces. The run streams live coordinates with energy and force plots, reports the energy change, steps, stop reason and final forces, and offers **Pause/Resume** and **Restore initial**. A run is one undo step and is interrupted when you switch tool, frame or molecule; limits are 3000 atoms, 500 steps or 5 s. Not available for crystals or in Diff.
- **Import structure** — Bring in a second structure from a file or a SMILES string, optionally as multiple copies; see [Importing a structure](#importing-a-structure).
- **Select Atoms** — Type atom numbers (1-based), ranges (`3-10`) or elements (`C H`) in the floating panel, or click atoms directly; selected atoms glow yellow. Clicking a selected atom deselects it, and selections accumulate.
- **Box Select** — Drag a rectangle to select every atom inside it. The selection carries over between Box Select and Select Atoms and into most editing tools.
- **Clipboard and batch delete** — **Ctrl+C / Ctrl+X / Ctrl+V** copy, cut and paste fragments with the viewer's own clipboard; **Ctrl+A** selects all atoms; **Delete** removes the selection. Pasting is not available for crystals or in Diff.
- **Undo** — Every edit is one undo step, the stack holds 50 entries, and it resets when the file is reloaded. **Cancel** in a dialog restores the coordinates it started from.
- **Save As** — Export to XYZ, Gaussian GJF, ORCA Input, MOL2, MDL Mol, PDB, MOPAC, Turbomole Coord, CIF, VASP POSCAR, Gaussian Cube or VESTA. GJF keeps its Link 0, route, original title, charge/multiplicity and connect section. For crystals, two extra "displayed cluster" formats (XYZ/GJF) write the region on screen as a molecule, with a **Keep whole molecules** switch. Crystal formats are available whenever a crystal is loaded.
- **Saved files record their history** — Title or comment lines receive a short provenance line such as `MV: from opt.log, frame 12 of 40; 2 atoms deleted, translated`. It follows edits and undo, is not stacked when a viewer-produced file is saved again, and formats without a title line are written as before.
- **Rotate Group** — Select a group and 2+ axis atoms, then rotate 0–360° around the orange axis with a slider or numeric input. Axis atoms stay fixed, bonded axis atoms must be joined by single bonds, each adjustment is one undo step, and **Reset** reverts the whole session.
- **Move Atoms** — Plain left-drag or arrow keys rotate the view and right-drag pans it; **Ctrl+left-drag** or **Ctrl+arrow keys** rotate the selected fragment, and **Ctrl+right-drag** translates it. Bonds stretch with the atoms and are re-judged on release; authored MOL2/PDB/GJF bonds are kept unless pulled beyond the cutoff. One undo step per drag. Not available for crystals.
- **Continuous editing** — After confirming an edit, the current mode stays active so you can make the next adjustment immediately.

#### Trajectory and Vibrations

- **Frame navigation** — ◀/▶ steps through optimisation frames, a jump field goes to any frame, and **Auto Play** runs at 350 ms per frame, stopping on the last (converged) geometry; pressing Play there starts again.
- **Route keywords** — Gaussian logs and ORCA outputs gain a collapsible **Route** panel that breaks the run's keywords into tokens and options, with an `opt` badge on optimisation jobs. The same parse decides whether the Convergence panel appears, so single-point `force`/`freq` jobs show no convergence curve.
- **Convergence** — The **Convergence** panel plots energy, max/RMS force and max/RMS displacement against step, with threshold lines taken from the file (ORCA tolerances or the Gaussian defaults). Force and displacement use log axes; hovering a point shows its step, value and energy. Gaussian `SCF Done`, external-optimizer energies, ORCA cycles and xtb gradient norms are supported.
- **Normal modes and manual displacement** — Frequency runs add a **Modes** panel listing harmonic frequencies and symmetries; imaginary modes are marked `*` and ORCA's zero translational/rotational modes are removed. Click a row (or ▶) to animate a mode. The circle at the row end selects a mode for **Manual Displacement**: set a scale from −1.00 to +1.00, combine several modes, and the structure becomes equilibrium + Σ(scale × mode). **Reset** returns to equilibrium and **Save As** writes the displaced geometry; playback always starts from the undisplaced equilibrium. Edits that replace the structure end the session.

#### Structure Diff

- **Any format** — Compare the current structure with any supported file; multi-frame logs ask which frame to use.
- **Skeleton check** — Structures match by connectivity (graph isomorphism using element, degree and neighbour signature), ignoring bond order, so different bond orders at the same connectivity still match.
- **Side-by-side view** — Each half rotates, pans and zooms independently and keeps its own aspect ratio.
- **Differences** — Matching skeletons highlight bond length, angle and dihedral differences in orange on both sides, with percentage changes and a complete difference list. A 0–20% threshold slider controls what is highlighted; the results panel can be closed and reopened with **Show Results**.
- **Statistics** — A collapsible **Internal Coordinate Statistics** block gives n, RMS, max |Δ|, mean and standard deviation per element-composition group and overall. It covers every compared coordinate, not only those above the threshold.
- **Auto-mirror selection** — Measuring a bond, angle or dihedral on one side selects the mapped atoms on the other and shows both values plus Δ.
- **Crystal Diff** — Crystals are compared as unit cells: cell parameters first, then the periodic base-cell graph, so a file compared with itself is identical at any supercell size. Distances are measured through periodic images, and every displayed copy of a differing site is marked. A crystal and a molecule are reported as not comparable, and Polyhedra is refused during Diff.

#### Crystal Structures

- **Unit cell** — CIF, VASP, Cube and VESTA files open with a wireframe cell (grid lines when the view spans several cells) and a/b/c edge labels.
- **Supercell bounds** — Six inputs define the displayed region in fractional coordinates, with 2×2×2, 3×3×3 and 4×4×4 presets plus **Reset** (0–1) and **Align Origin**. The region is a cut-out of the crystal, not a deformation; atoms outside the range are simply not shown, and finite molecules are shown whole even across a boundary. The view keeps its centre and selection when bounds change, and a change that would exceed 30,000 atoms is refused with a message.
- **Atom identity and cell info** — Hovering shows Cartesian and fractional coordinates, the periodic image `cell (i,j,k)`, site label and occupancy. The info panel adds cell volume, Z and density plus a collapsible asymmetric-unit table; expanded images name their site and symmetry operation.
- **Local edits** — Bond-length, angle, dihedral and Rotate Group edits are written back into the unit cell as fractional displacements of the base site, appear in every periodic copy and survive crystal export. Each site is edited once. Bonds drawn as stubs across a cell face cannot be edited — grow the supercell instead. Add Atom wraps new sites into the cell, and the Bond Length dialog reports the shortest periodic image.
- **Bonds and symmetry** — Bonds use a minimum-image search (27 images normally, more for narrow or skewed cells) and the same perception module as molecules, including aromatic orders, and are drawn as split stubs at cell faces. A CIF's own `_geom_bond` table is read with its symmetry codes and stated order, and unusable rows are reported. CIFs that state only a space-group name use a built-in table of all 230 groups (530 settings); unreadable operations are reported as P1.
- **Neighbours** — Lists contacts around a clicked atom within an adjustable radius (3.5 Å default, up to 12 Å), skipping bonded atoms unless enabled. Rows name the periodic image, clicking selects both ends, and **Select listed** takes the whole shell. Escape or View clears the dashed guides.
- **Lattice planes** — Two h k l triples give d-spacing and the angle between normals; **Show Plane** draws the plane through the centre of the box. Hidden in Diff.
- **Powder XRD** — Calculates a pattern for Cu/Mo/Co/Fe/Cr/Ag Kα radiation and 2θ limits of 60/90/120/150°, with real systematic absences, enumerated multiplicities, Debye–Waller damping and the Kα1/Kα2 doublet. The panel plots the curve, lists the strongest lines and notes elements outside the form-factor table; it runs in the background with progress and Cancel, and exports as `.xy`. Not available in Diff.
- **Cavities** — Probes the unit cell with an adjustable probe radius on a 0.25 Å grid, reporting open channels versus closed cavities, per-cell volume, centre and largest fitting sphere. Click a row to draw the void as a translucent surface; very large voids report their volume without drawing. Read-only and not available in Diff.
- **Edit Cell** — Change a, b, c, α, β, γ while sites keep their fractional coordinates; invalid values or a cell above the atom ceiling are refused with the reason shown. Every export uses the edited cell.
- **Disorder groups** — Choose **All** or a disorder group; the selected group can be counted at full occupancy in the formula, density and exports. Switching back to All is lossless.
- **Occupancy and cleanup** — Partially occupied atoms are drawn partially coloured. **Remove Disorder** deletes sites below 0.5 occupancy, sets those above 0.5 to 1 and keeps 0.5.
- **Crystal save** — CIF, VASP POSCAR, Cube and VESTA all write the unit cell rather than the supercell. CIF/VESTA declare P1 with wrapped fractional coordinates; POSCAR writes the cell; Cube writes a real cubic grid with a chosen number of points. Use the cluster XYZ/GJF formats to export the displayed region as a molecule.

#### VMD TCL Scripts

Opening a `.tcl` file reproduces the script's display recipe instead of applying the viewer's defaults.

- **Representations** — Every `mol new` becomes a selectable molecule. CPK, VDW, Licorice, Lines and Surf styles are reproduced with their sphere scales, bond radii and probe radii; `mol modstyle`, `mol modcolor`, `mol modselect`, `mol modmaterial` and `mol delrep` are honoured, including the default representation created by `mol new`.
- **Selections** — `index`, `to` ranges, `element`/`name`/`type`, `not`, `all` and `none` are resolved against the loaded atoms. Selections that need residue metadata are widened to the whole molecule and reported.
- **Sizing** — Script sizes are applied as ratios to the viewer's own ball-and-stick geometry, so `CPK 0.97` looks like an ordinary structure and `CPK 0.1` stays tiny; probe radii remain real Å values. Atoms not covered by any representation are hidden unless Display is Space-filling or Stick.
- **Colours** — `Beta` maps per-atom scalars onto the reference red → green → blue ramp using the script's `mol scaleminmax` range or the data's own range. `ColorID` uses the 33-colour reference palette including `color change rgb` overrides; `Name`/`Element` uses the reference element colours.
- **Surf surfaces** — `Surf` representations become real solvent-accessible isosurfaces of the selected atoms, drawn translucent and moving with the structure.
- **Egg list** — With `eggFile*` and `unitfile` metadata present, an **Eggs** panel lists each egg with its colour and units and can hide or show individual eggs. If several metadata files exist, they are scored against the script's patches and the best match is offered.
- **Bare PDB visualizations** — A header-less PDB whose column 61–66 is non-constant opens already coloured by that scalar, with the range on the status bar. **Color by Beta** toggles the colouring; constant columns are ignored.
- **Background** — The script's white background is deliberately not applied; the viewer keeps its theme background. `color scale method` and `color change rgb` are applied.

A script also works as the comparison side of **Structure Diff**, using the script's colours.

#### Supported File Formats

| Format | Extensions | What is read |
|--------|------------|--------------|
| Gaussian input | `.gjf` `.gjf03` `.gjf09` `.gjf16` `.com` | Link 0, route, title, charge/multiplicity, coordinates, connect section, fixed-atom notation |
| XYZ | `.xyz` | Standard XYZ; multi-frame trajectories (including ORCA `_trj.xyz`) with frame navigation and convergence when every frame has an energy; a truncated final frame is dropped |
| MOL2 | `.mol2` | ATOM and BOND sections, aromatic `ar` bonds |
| Gaussian log | `.log` | Optimisation trajectories, SCF energies, convergence, harmonic frequencies; xtb trajectories are detected automatically |
| ORCA input | `.inp` | `* xyz … *` and `%coords` blocks, xyz/xyzfile data, Bohr/Å units |
| ORCA output | `.out` | Cartesian coordinates in Å or a.u., cycles and final geometry, charge/multiplicity, convergence tolerances, frequencies/IR, route line; format detected by content |
| Turbomole | `.coord` | `$coord` (Bohr → Å), `$chrg`, `$spin`/`$mult` |
| PDB | `.pdb` `.ent` | ATOM/HETATM, element, occupancy, B-factor scalar, CONECT, duplicate serials |
| MOPAC | `.mop` `.mopac` `.dat` | Internal coordinates, Å/Bohr auto-detection, CHARGE and MS, Cartesian fallback |
| VMD TCL | `.tcl` | Display recipe with molecules, representations, colours, selections, surfaces and egg metadata |
| CIF | `.cif` | Cell, space group and symmetry operations, fractional/Cartesian sites, occupancy, asymmetric unit, `_geom_bond`; CIF 1.1, 1.2 and 2 tags |
| VASP POSCAR | `.vasp` `POSCAR` `CONTCAR` | Scale (negative = target volume), lattice, element counts, selective dynamics, Direct/Cartesian |
| Gaussian Cube | `.cube` | Voxel vectors and atoms, Bohr/Å auto-detection, robust lattice fallback |
| VESTA | `.vesta` | CELLP, STRUC, GROUP and SYMOP; symmetry expansion and occupancy, with the same logic as CIF |

### Installation

#### From VSIX (Recommended)

1. Download the latest `.vsix` from [Releases](https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode/releases).
2. In VS Code or Trae, press `Ctrl+Shift+P`.
3. Run `Extensions: Install from VSIX...` and select the downloaded file.
4. For Remote-SSH, install the extension **on the remote server** when prompted.

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

#### Opening Files

1. Right-click a supported file in the Explorer and choose **Molecular Viewer: Open 3D Viewer**.
2. Or open the Command Palette (`Ctrl+Shift+P`) and run `Molecular Viewer: Open 3D Viewer`.
3. Or double-click a supported file and choose **Molecular 3D Viewer** as the editor.

#### Importing a Structure

1. Click **Import Structure** in the toolbar (View mode; not available for crystals or in Diff).
2. Enter the number of copies (1–50) and choose **From File** or **From SMILES**. SMILES input is validated as you type with the parser's error message; press Enter to build the 3D structure.
3. A request that would exceed 20,000 atoms is trimmed with a warning.
4. Copies are placed beside the existing structure without overlap and spaced along +X. The whole import is one Undo step, and Save As exports everything together.

SMILES support covers the organic subset plus branches, `%nn` ring closures, aromaticity, explicit bonds, disconnected fragments, isotopes, charges and `@`/`@@` tetrahedral configuration. Aromatic systems are planar, unspecified double bonds default to trans, and the conformation comes from local rules rather than a conformer search; organometallic geometry is not special-cased.

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

The list mirrors the extensions the viewer registers; drop any line you do not want to divert from your normal editor.

#### Controls

| Action | Effect |
|--------|--------|
| Left-drag | Rotate around the molecule centre |
| Scroll | Zoom in / out |
| Middle/right-drag | Pan (cursor-locked; in Diff mode only the grabbed half) |
| Hover atom | Show element and coordinates |
| Arrow keys | Rotate (←/→ around Y, ↑/↓ around X); with Auto Rotate on, a held key spins the view up |
| `Esc` | One step back: cancel a dialog or box select, leave the current tool for View, or clear the selection |
| `Ctrl+Z` / `Cmd+Z` | Undo (ignored while a dialog or text field has focus) |
| `Ctrl+S` / `Cmd+S` | Open Save As; press again to save the chosen format |
| `Ctrl+A` / `Cmd+A` | Select every atom (both halves in Diff) |
| `Ctrl+C` / `Ctrl+X` / `Ctrl+V` | Copy, cut and paste fragments using the viewer's clipboard |
| `Delete` | Delete the selection, or enter Delete Atom mode when nothing is selected |
| ◀ ▶ (molecule selector) | Switch between molecules created by a VMD script |
| Color by Beta | Toggle scalar colouring for PDB scalar files |
| Auto Rotate | Toggle flick-to-spin in View mode |
| Reset View | Return to the initial view orientation |

Shortcuts work when the 3D view has focus and do not override editor keybindings; text fields keep their own keys. Copy and paste use the viewer's own clipboard, so they never replace text you copied elsewhere.

#### Editing Workflow

1. Click a toolbar button to enter an editing mode (for example **Bond Length**).
2. Click atoms in the 3D view; selected atoms glow yellow.
3. The dialog shows the current value.
4. Choose which atoms stay fixed and which move.
5. Adjust the value with the numeric input or slider; the 3D view previews in real time.
6. Click **OK** to confirm or **Cancel** to revert.
7. Use **Save As** to export the modified structure.

### Development

```bash
npm install       # Install dependencies
npm run compile   # Compile TypeScript
npm run watch     # Compile on change
npm run lint      # Type-check without emitting
npm run package   # Package as .vsix
```

### Contributing

Fork the repository, create a feature branch, commit your changes, push the branch and open a Pull Request. Please keep changes focused and describe what they change.

### Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full release history.

### License

[MIT](LICENSE)

---

## 简体中文

**简体中文** | [English](#english)

一个把 GaussView 风格的 3D 分子查看与编辑功能带进 VS Code / Trae 的插件。它可以打开 Gaussian、ORCA、MOPAC、xtb、VMD 的输出以及常见结构格式，交互式渲染与编辑，浏览优化轨迹、比较结构并查看晶体——包括通过 Remote-SSH 编辑远程服务器上的文件。

**目录**

- **功能** — [功能特性](#功能特性) · [分子编辑](#分子编辑) · [轨迹与振动](#轨迹与振动) · [结构对比](#结构对比) · [晶体结构](#晶体结构) · [VMD TCL 脚本](#vmd-tcl-脚本) · [支持的文件格式](#支持的文件格式)
- **快速上手** — [安装](#安装) · [使用](#使用) · [打开文件](#打开文件) · [导入结构](#导入结构) · [操作与快捷键](#操作与快捷键)
- **其他** — [开发](#开发) · [贡献](#贡献) · [更新日志](#更新日志) · [许可证](#许可证)

### 功能特性

- **四种显示样式** — **球棍**（默认）、**空间填充**、**Stick** 与**多面体**。多面体会为每个至少有四个非中心邻居的金属或类金属中心画出半透明配位多面体，并跨越晶面补全配位壳；邻居不足、壳层共面或全为金属的中心仍按球棍绘制。多面体只在晶体页面提供，Diff 模式下不可用。在 **View** 卡片的 **Display** 一节选择样式：样式属于视图，编辑、切换帧、撤销与优化都不会改变它。幽灵原子始终保持小标记，任何样式下原子都点得中。
- **键级** — 单键、芳香键（1 实线 + 1 虚线）、双键和三键显示不同。芳香虚线按芳香体系统一规划，环内全部朝内，链上则保持连续的一条边。
- **自动键感知** — 文件没有显式连接关系时，按键长与元素对规格推断化学键与键级，并用价态/电荷检查正确处理硝基、铵离子、羧酸根等结构，同时推导形式电荷；价态异常的原子会被报告。GJF connect 段、MOL2 与 PDB 中写明的键会保留；移动原子后会重新判定连接性。
- **分子信息面板** — 左上角显示 Hill 分子式、原子数、电荷、电子数与自旋多重度；晶体文件额外显示整个晶胞的组成。解析器若不得不绕开文件中的问题，会用琥珀色行说明，而不是静默隐藏。
- **非共价相互作用** — **Interactions** 模式由近及远列出氢键（H···A < 2.5 Å，D–H···A > 120°）、卤键/硫族键/氮族键接触、π–π 堆积（平行、平行错位或 T 型）与 CH–π 接触。每类各有复选框、计数与虚线；点击一行会选中两端原子。晶体搜索包含周期像。Escape 或 View 关闭模式。
- **两行工具栏** — 编辑工具在第一行，接触、选择与视图/操作工具在第二行；两行各自横向滚动，窄面板中也能点到每个按钮。
- **可折叠卡片** — 每张浮动卡片都能用 **×** 收起，并在屏幕边缘留下一个带名字的标签；点击标签即可在原位置、原尺寸恢复。布局会自动收窄与堆叠卡片，使其不会互相遮挡或跑出窗口。
- **鼠标与触摸** — 拖拽旋转、滚轮缩放、中/右键平移（抓取点始终锁定在光标下），悬停原子显示元素与坐标。方向键旋转，Delete 删除选区（无选区时进入删除原子模式），触摸支持单指旋转与双指捏合缩放。完整列表见 [操作与快捷键](#操作与快捷键)。
- **自动旋转** — 在 View 模式右下角面板中开启 **Auto Rotate**，即可把视角「甩」出去：旋转拖拽松手后，视角按松手时的方向与速度继续自转（只统计最后 120 ms），按住方向键则会加速自转。下一次左键按下停止，右键/中键平移不会打断。它只旋转相机，不改动坐标。
- **Alt 拖拽写入坐标** — 在 View 模式下按住 **Alt** 拖拽：左键旋转结构，右键/中键平移结构。移动的是原子而不是相机，因此 **另存为** 写出的就是拖到的朝向。每次拖拽一步 Undo，且刚体变换不改变距离，所以不会重推化学键。Diff 模式与晶体结构中仍保持纯视图拖拽。
- **方向坐标轴** — 左下角显示分子的笛卡尔 X/Y/Z 轴或晶体的晶格 a/b/c 轴，并跟随视图旋转；**Reset View** 恢复初始朝向。
- **远程 SSH** — 解析与保存运行在远程主机，3D 视图在本地渲染，因此在远程 Linux 服务器上编辑文件也能直接使用。
- **VMD TCL 脚本** — 打开 VMD 脚本会复现它的表示、颜色、选择与曲面，并支持多分子与可选的 egg 元数据。详见 [VMD TCL 脚本](#vmd-tcl-脚本)。
- **大结构性能** — 场景只在交互或变化时重绘，空闲时不占用 GPU；重建时释放显存。细节层次随结构规模调整，球与键采用实例化绘制，数千原子的结构仍可流畅操作。Move Atoms 只重建受影响的键，长日志单趟解析。

#### 分子编辑

- **键长** — 选择 2 个原子，固定一侧，选择**平移基团**（连同连通子树一起移动）或**平移原子**，再拖动滑块或输入数值；3D 视图实时预览。同一个对话框还可修改键级（无 / 单键 / 芳香 / 双键 / 三键，选择「无」会删除该键）。
- **键角** — 选择 3 个原子（第 2 个为顶点），固定或移动任一侧，用滑块或数值输入调节，实时预览。
- **二面角** — 选择 4 个原子，用同样的方式调节。
- **添加原子** — 点击锚点原子，选择元素（70+ 种），设定键长与键级。方向会最大化与已有邻居的夹角，自然得到直线、三角、四面体或八面体构型；对同一锚点重复添加会扇形展开。
- **删除原子** — 点击原子并确认，原子与键会自动重新编号。要一次删除多个原子，先选中它们再按 Delete。
- **替换原子** — 点击 **Replace Atom**，在 3D 视图中点击一个或多个原子（再次点击取消选中），然后再次点击按钮，从 70+ 种元素中选择新元素。坐标保持不变，晶体的基础原子也可以替换。
- **幽灵原子** — 一个只有位置、没有原子核、电子和化学键的标记，用于 counterpoise（BSSE）计算。可在添加原子弹窗中打开 Ghost 开关或选择 `X (dummy point)` 元素；替换原子弹窗的 Ghost 开关可以把整个选区变成各自元素的幽灵原子。幽灵原子不参与成键、分子式、电子数、补氢与 MM，绘制为淡色小标记，并能在各格式中往返（Gaussian `C-Bq`、ORCA `C:`、XYZ/PDB 占位标记等）。晶体结构中不允许幽灵原子。
- **补氢** — **Add H** 为所选原子补齐缺失的价，采用元素专属的键长、键角与中心杂化态。已有原子绝不会移动，重复运行不会新增，放不下（会与别的原子重叠）的氢会被跳过并报告，整次添加是一步 Undo。晶体结构与 Diff 模式下不可用。
- **MM 优化** — **MM Optimize** 面板运行内置、零依赖的力场。可以逐项设置参数，或一键选择 **Low / Medium / High** 三档预设，再选择要松弛的原子并开始。只有选中的原子会移动，但固定原子仍通过键、角、二面角与非键作用施力。运行中实时更新坐标，并绘制能量与受力曲线；结束后报告能量变化、步数、停止原因与最终受力，支持 **Pause/Resume** 与 **Restore initial**。整次运行是一步 Undo，切换工具、帧或分子会中断运行；上限为 3000 原子、500 步或 5 秒。晶体结构与 Diff 模式下不可用。
- **导入结构** — 从文件或 SMILES 字符串导入第二个结构，可选择副本数量；详见 [导入结构](#导入结构)。
- **选择原子** — 在浮动面板中输入原子编号（从 1 开始）、范围（如 `3-10`）或元素（如 `C H`），也可以直接在 3D 视图中点击原子；选中的原子显示为黄色。再次点击已选原子会取消选中，选择会累积。
- **框选** — 拖出一个矩形，选中框内的所有原子。选择会在「框选」与「选择原子」之间保留，并可带入大多数编辑工具。
- **剪贴板与批量删除** — **Ctrl+C / Ctrl+X / Ctrl+V** 使用查看器自己的剪贴板复制、剪切与粘贴片段；**Ctrl+A** 全选原子；**Delete** 删除选区。晶体结构与 Diff 模式下不可粘贴。
- **撤销** — 每次编辑是一步 Undo，栈上限 50 条，重新载入文件时重置。对话框中的 **Cancel** 会恢复到打开对话框时的坐标。
- **另存为** — 可导出为 XYZ、Gaussian GJF、ORCA Input、MOL2、MDL Mol、PDB、MOPAC、Turbomole Coord、CIF、VASP POSCAR、Gaussian Cube 或 VESTA。GJF 会保留 Link 0、route、原标题、电荷/多重度与 connect 段。对晶体结构，额外的两个「当前显示团簇」格式（XYZ/GJF）会把屏幕上的区域当作分子写出，并提供 **Keep whole molecules** 开关。载入晶体数据时才会提供晶体格式。
- **保存的文件记录来历** — 标题行或注释行会写入一行简短的来历，例如 `MV: from opt.log, frame 12 of 40; 2 atoms deleted, translated`。它跟随编辑与 Undo，重复保存不会层层叠加；没有标题行的格式保持原样。
- **绕轴旋转基团** — 选择一个原子组和 2 个以上轴原子，然后用滑块或数值输入让原子组绕橙色轴线旋转 0–360°。轴原子保持不动，成键的轴原子之间必须是单键；每次调节是一步 Undo，面板的 **Reset** 一次性回退整个旋转会话。
- **移动原子** — 普通左键拖拽或方向键旋转视图，右键拖拽平移视图；**Ctrl+左键拖拽**或 **Ctrl+方向键**旋转选中片段，**Ctrl+右键拖拽**平移片段。拖拽时键随原子拉伸，松开后重新判定连接性；MOL2/PDB/GJF 中写明的键会保留，除非被拉到成键范围之外。每次拖拽一步 Undo。晶体结构下不可用。
- **连续编辑** — 确认一次编辑后，当前模式保持激活，可以立即进行下一次调整。

#### 轨迹与振动

- **帧浏览** — ◀/▶ 在优化帧之间步进，跳转输入框可直达指定帧；**Auto Play** 以每帧 350 ms 播放，在最后一个（收敛）结构处停止，在该帧再次点击播放则从头开始。
- **路由关键字** — Gaussian 日志与 ORCA 输出会多出一个可折叠的 **Route** 面板，把运行的关键字拆成词元与选项，优化作业带 `opt` 徽标。同一份解析结果也决定 Convergence 面板是否出现，因此单点的 `force`/`freq` 作业不会显示收敛曲线。
- **收敛面板** — **Convergence** 面板绘制能量、最大/RMS 力与最大/RMS 位移随步数的变化，阈值线取自文件（ORCA 的收敛容差或 Gaussian 的默认值）。力与位移使用对数轴，悬停数据点可查看步数、数值与能量。支持 Gaussian 的 `SCF Done`、外部优化器能量、ORCA 的每个循环以及 xtb 的梯度范数。
- **简正模式与手动位移** — 频率计算会添加 **Modes** 面板，列出谐振频率与对称性；虚频以 `*` 标记，ORCA 的零平动/转动模式会被剔除。点击行（或 ▶）播放振动动画。行尾的小圆圈把模式指定给 **Manual Displacement**：设置 −1.00 到 +1.00 的比例、可叠加多个模式，结构会变为「平衡构型 + Σ(比例 × 模式)」。**Reset** 回到平衡构型，**Save As** 可写出位移后的几何；播放始终从无位移的平衡构型开始。会替换结构的编辑会结束该会话。

#### 结构对比

- **任意格式** — 把当前结构与任一受支持文件比较；多帧日志会询问使用哪一帧。
- **骨架校验** — 按连接性做图同构匹配（元素、度数、邻居签名），不考虑键级，因此连接性相同但键级写得不一致的分子仍能匹配。
- **并排视图** — 左右两半可各自独立旋转、平移与缩放，并保持各自的宽高比。
- **差异显示** — 骨架匹配时，两侧以橙色高亮键长、键角与二面角差异，显示百分比变化和完整差异列表。0–20% 的阈值滑块控制高亮内容；结果面板可用 **Show Results** 关闭与重开。
- **分组统计** — 可折叠的 **Internal Coordinate Statistics** 区块按元素组成分组，给出 n、RMS、最大 |Δ|、均值与标准差，并附一行总体统计。统计覆盖所有比较过的坐标，而不只是超过阈值的那些。
- **自动镜像选择** — 在一侧测量键长、键角或二面角时，另一侧会自动选中对应原子，并同时显示两侧数值与 Δ。
- **晶体对比** — 晶体按晶胞比较：先核对六个晶胞参数，再匹配基础晶胞的周期性图，因此同一文件在任何超胞尺寸下都与自身相同。距离沿周期像测量，差异位点的每个显示副本都会标记。晶体与分子会明确提示不可比较，Diff 模式下也不能使用多面体。

#### 晶体结构

- **晶胞显示** — 打开 CIF、VASP、Cube 或 VESTA 文件时，晶胞以线框显示（显示多个晶胞时内部晶面线更淡），原点相交的三条棱标注 a/b/c。
- **超胞边界** — 6 个输入框以分数坐标定义显示区域，并提供 2×2×2、3×3×3、4×4×4 预设以及 **Reset**（0–1）与 **Align Origin**。显示区域是晶体的剖切，不是晶格的形变；范围之外的原子不显示，有限分子即使跨越边界也完整显示。修改边界时视图保持中心与选择不变；若显示原子数将超过 30,000，会拒绝修改并提示。
- **原子身份与晶胞信息** — 悬停原子可看到笛卡尔坐标、分数坐标、周期像 `cell (i,j,k)`、位点标签与占据率。信息面板额外显示晶胞体积、Z 与密度，以及可折叠的不对称单元表；展开的像会注明其位点与对称操作。
- **晶体局部修改** — 键长、键角、二面角与绕轴旋转会以基础位点的分数位移写回晶胞，出现在每个周期副本中，并在导出 CIF/VESTA 后继续存在。同一基础位点只编辑一次。被画成跨晶面短桩的键不能编辑，应扩大超胞；添加原子会把新位点 wrap 回晶胞，键长对话框会报告最短周期像。
- **成键与对称性** — 晶体键使用最小镜像搜索（普通晶胞枚举 27 个像，窄或斜的晶胞会自动扩大），并采用与分子相同的感知模块（含芳香键级），跨晶面的键绘制为延伸到视图外的短分割键。CIF 自带的 `_geom_bond` 表会连同对称码与写明的键级一起读取，无法使用的行会被报告。只写空间群名的 CIF 使用内置的 230 个空间群（530 个设置）表；无法读取的对称操作会按 P1 报告。
- **邻居** — 列出点击原子周围的接触，半径可调（默认 3.5 Å，最大 12 Å），默认跳过已成键的原子（可用复选框恢复）。每行注明周期像，点击一行会选中两端，**Select listed** 选中整个壳层。Escape 或 View 清除虚线。
- **晶面** — 两组 h k l 给出 d 间距与法线夹角；**Show Plane** 在显示盒子中心画出该晶面。Diff 中不显示。
- **粉末 XRD** — 按 Cu/Mo/Co/Fe/Cr/Ag Kα 辐射与 60/90/120/150° 的 2θ 上限计算衍射图，包含真实的系统消光、枚举出的多重度、Debye–Waller 衰减与 Kα1/Kα2 双线。面板绘制曲线、列出最强线并标注形状因子表之外的元素；计算在后台分片进行，带进度与 Cancel，并可把曲线导出为 `.xy`。Diff 模式下不可用。
- **空腔** — 用可调探针半径在 0.25 Å 网格上扫描晶胞，区分开放通道与闭合空腔，并给出单胞体积、中心与可容纳的最大球。点击一行可把空腔绘制为半透明曲面；体素过多时只报告体积。只读，Diff 模式下不可用。
- **编辑晶胞** — 修改 a、b、c、α、β、γ，所有位点保持分数坐标；非法数值或会超过原子上限的晶胞会被拒绝并说明原因。所有导出都使用修改后的晶胞。
- **无序组** — 在 **All** 与某个无序组之间切换，所选组可在分子式、密度与导出中按满占据计算。切回 All 不会丢失任何信息。
- **占据与清理** — 部分占据的原子按占据率部分着色。**Remove Disorder** 删除占据率低于 0.5 的位点，把高于 0.5 的置为 1，等于 0.5 的保留。
- **晶体保存** — CIF、VASP POSCAR、Cube 与 VESTA 都写出晶胞而不是超胞。CIF/VESTA 声明 P1 并写入 wrap 后的分数坐标；POSCAR 写出晶胞；Cube 写出真正的立方网格（网格点数可选）。要把屏幕上的区域导出为分子，请使用 XYZ/GJF 的「当前显示团簇」格式。

#### VMD TCL 脚本

打开 `.tcl` 文件时，会按脚本要求的方式渲染，而不是套用查看器的默认值。

- **表示方式** — 每个 `mol new` 都会成为一个可切换的分子。CPK、VDW、Licorice、Lines 与 Surf 样式连同球缩放、键半径与探针半径都会被复现；`mol modstyle`、`mol color`、`mol selection`、`mol material` 与 `mol delrep` 都会生效，包括 `mol new` 创建的默认表示。
- **选择** — `index`、`to` 区间、`element`/`name`/`type`、`not`、`all` 与 `none` 会针对已载入原子求解。需要残基信息的选择会被放宽为整个分子并明确报告。
- **尺寸** — 脚本中的尺寸按比例套用到查看器自己的球棍几何上：`CPK 0.97` 看起来是普通结构，`CPK 0.1` 则很小；探针半径仍按真实 Å 值使用。未被任何表示覆盖的原子不会绘制，除非 Display 选择空间填充或 Stick。
- **颜色** — `Beta` 把逐原子标量映射到参考渲染器的红 → 绿 → 蓝渐变，范围取自脚本的 `mol scaleminmax` 或数据自身范围。`ColorID` 使用 33 色参考调色板，并应用 `color change rgb` 覆盖；`Name`/`Element` 使用参考渲染器的元素颜色。
- **Surf 曲面** — `Surf` 表示会变成所选原子的真实溶剂可及等值面，半透明绘制，并随结构一起旋转与平移。
- **Egg 列表** — 当存在 `eggFile*` 与 `unitfile` 元数据时，**Eggs** 面板会列出每个蛋及其颜色和包含的单元，并可单独隐藏或显示。若存在多个元数据文件，会按与曲面片的匹配数评分并推荐最合适的一份。
- **裸 PDB 可视化** — 没有头部记录且 61–66 列非常量的 PDB 会按该列标量自动上色，范围显示在状态栏。**Color by Beta** 可切换着色；常量列不着色。
- **背景** — 脚本的白色背景有意不套用，查看器保留自己的主题背景；`color scale method` 与 `color change rgb` 会生效。

脚本也可以作为**结构对比**的一侧，并按脚本的颜色显示。

#### 支持的文件格式

| 格式 | 扩展名 | 读取内容 |
|------|--------|----------|
| Gaussian 输入 | `.gjf` `.gjf03` `.gjf09` `.gjf16` `.com` | Link 0、route、标题、电荷/多重度、坐标、connect 段、固定原子写法 |
| XYZ | `.xyz` | 标准 XYZ；多帧轨迹（含 ORCA `_trj.xyz`）支持帧导航，每帧都带能量时显示收敛信息；被截断的末帧会被丢弃 |
| MOL2 | `.mol2` | ATOM 与 BOND 段、芳香 `ar` 键 |
| Gaussian 日志 | `.log` | 优化轨迹、SCF 能量、收敛、谐振频率；xtb 轨迹会自动识别 |
| ORCA 输入 | `.inp` | `* xyz … *` 与 `%coords` 块、xyz/xyzfile 数据、Bohr/Å 单位 |
| ORCA 输出 | `.out` | Å 或 a.u. 的笛卡尔坐标、每个循环与最终几何、电荷/多重度、收敛容差、频率/IR、route 行；按内容识别格式 |
| Turbomole | `.coord` | `$coord`（Bohr → Å）、`$chrg`、`$spin`/`$mult` |
| PDB | `.pdb` `.ent` | ATOM/HETATM、元素、占据率、B 因子标量、CONECT、重复序号 |
| MOPAC | `.mop` `.mopac` `.dat` | 内坐标、Å/Bohr 自动识别、CHARGE 与 MS、笛卡尔坐标回退 |
| VMD TCL | `.tcl` | 显示配方：分子、表示、颜色、选择、曲面与 egg 元数据 |
| CIF | `.cif` | 晶胞、空间群与对称操作、分数/笛卡尔位点、占据率、不对称单元、`_geom_bond`；支持 CIF 1.1、1.2 与 2 标签 |
| VASP POSCAR | `.vasp` `POSCAR` `CONTCAR` | 缩放因子（负值表示目标体积）、晶格、元素数量、选择性动力学、Direct/Cartesian |
| Gaussian Cube | `.cube` | 体素矢量与原子、Bohr/Å 自动识别、晶格损坏时的稳健回退 |
| VESTA | `.vesta` | CELLP、STRUC、GROUP 与 SYMOP；对称展开与占据率，逻辑与 CIF 一致 |

### 安装

#### 从 VSIX 安装（推荐）

1. 从 [Releases](https://github.com/Atreides-Jimmy/Molecular-Viewer-on-VScode/releases) 下载最新的 `.vsix` 文件。
2. 在 VS Code 或 Trae 中按 `Ctrl+Shift+P`。
3. 运行 `Extensions: Install from VSIX...` 并选择下载的文件。
4. 远程 SSH 用户请务必在**远程服务器上**安装该扩展（弹出提示时选择「Install on Remote」）。

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

#### 打开文件

1. 在资源管理器中右键点击受支持的文件，选择 **Molecular Viewer: Open 3D Viewer**。
2. 或按 `Ctrl+Shift+P` 打开命令面板，运行 `Molecular Viewer: Open 3D Viewer`。
3. 或双击受支持的文件，选择 **Molecular 3D Viewer** 作为编辑器。

#### 导入结构

1. 在工具栏点击 **Import Structure**（View 模式；晶体结构与 Diff 模式下不可用）。
2. 输入副本数量（1–50），并选择 **From File** 或 **From SMILES**。SMILES 输入会边输入边校验，并显示解析器的错误信息；按 Enter 构建 3D 结构。
3. 若请求会超过 20,000 原子，会自动裁剪并给出提示。
4. 副本会放在已有结构旁边、互不重叠，并沿 +X 方向等距排开。整次导入是一步 Undo，Save As 会导出全部结构。

SMILES 支持有机子集以及分支、`%nn` 环闭合、芳香性、显式键、断开片段、同位素、电荷和 `@`/`@@` 四面体构型。芳香体系保持平面，未指定构型的双键默认为反式，构象由局部规则生成而非构象搜索；有机金属几何不做特殊处理。

#### 设为默认查看器

在 `settings.json` 中加入：

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

该列表与查看器注册的扩展名一致；不需要用查看器打开的行可以删掉。

#### 操作与快捷键

| 操作 | 效果 |
|------|------|
| 左键拖拽 | 绕分子中心旋转 |
| 滚轮 | 缩放 |
| 中键/右键拖拽 | 平移（抓取点锁定在光标下；Diff 模式下只平移抓到的那一半） |
| 悬停原子 | 显示元素与坐标 |
| 方向键 | 旋转（←/→ 绕 Y 轴，↑/↓ 绕 X 轴）；开启 Auto Rotate 后按住方向键会加速自转 |
| `Esc` | 后退一步：取消对话框或框选、从当前工具回到 View、或清除选择 |
| `Ctrl+Z` / `Cmd+Z` | 撤销（对话框或文本输入框获得焦点时不触发） |
| `Ctrl+S` / `Cmd+S` | 打开另存为；对话框已打开时再按一次保存当前格式 |
| `Ctrl+A` / `Cmd+A` | 全选原子（Diff 模式下两侧全选） |
| `Ctrl+C` / `Ctrl+X` / `Ctrl+V` | 使用查看器剪贴板复制、剪切、粘贴片段 |
| `Delete` | 删除选区；没有选区时进入删除原子模式 |
| ◀ ▶（分子选择器） | 在 VMD 脚本创建的多个分子之间切换 |
| Color by Beta | 切换 PDB 标量文件的逐原子着色 |
| Auto Rotate | 在 View 模式下切换甩动自转 |
| Reset View | 恢复初始视图朝向 |

快捷键只在 3D 视图获得焦点时生效，不会占用编辑器的快捷键；文本输入框保留自己的按键。复制粘贴使用查看器自己的剪贴板，不会替换你在别处复制的文本。

#### 编辑流程

1. 点击工具栏按钮进入编辑模式（例如 **Bond Length**）。
2. 在 3D 视图中点击原子；选中的原子显示为黄色。
3. 对话框显示当前数值。
4. 选择哪些原子固定、哪些移动。
5. 用数值输入或滑块调节，3D 视图实时预览。
6. 点击 **OK** 确认，或 **Cancel** 还原。
7. 用 **Save As** 导出修改后的结构。

### 开发

```bash
npm install       # 安装依赖
npm run compile   # 编译 TypeScript
npm run watch     # 监听改动并编译
npm run lint      # 只做类型检查
npm run package   # 打包为 .vsix
```

### 贡献

Fork 仓库，创建功能分支，提交改动，推送分支并发起 Pull Request。请保持改动聚焦，并说明改动内容。

### 更新日志

完整发布历史见 [CHANGELOG.md](CHANGELOG.md)。

### 许可证

[MIT](LICENSE)
