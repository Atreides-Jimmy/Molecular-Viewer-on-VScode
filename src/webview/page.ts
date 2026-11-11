// Build the webview page: inline the shipped assets (media/three.min.js,
// media/viewer.css, media/viewer.js) with a nonce, add this document's payload
// as window.__MV_HOST__, and hand it to the panel.
import * as vscode from 'vscode';
import { Atom, AsymAtom, MolecularData, NormalMode, OptStep, OptThresholds, OrbitalLevels, SolventCavity } from '../types';
import { FchkOrbitals, LogChargeSet, LogFrame, LogSummary, OrcaFrame, RouteSection } from '../parsers/index';
import { BOND_PERCEPTION_SOURCE } from '../tools/bondPerception';
import { BVS_SOURCE } from '../tools/bvs';
import { describeForceField } from '../tools/mmOptimizer';
import { ATOM_COLORS, TrajectoryFramePayload, VmdPayload } from './payloads';
import { decodeText, yieldToHost } from './loading';

/**
 * Everything the page needs to be built, gathered from the document.
 *
 * The viewer page is an asset pair (media/viewer.css and media/viewer.js) the
 * host assembles with this document's payload: keeping the builder out of the
 * provider means the provider only opens documents, and the page only renders
 * one.
 */
export interface ViewerPageOptions {
    data: MolecularData;
    frames?: (LogFrame | OrcaFrame)[];
    optSteps?: OptStep[];
    normalModes?: NormalMode[];
    logSource?: string;
    routes?: RouteSection[];
    vmdPayload?: VmdPayload;
    autoBeta?: boolean;
    optThresholds?: OptThresholds;
    trajectoryFrames?: TrajectoryFramePayload[];
    charges?: LogChargeSet[];
    summary?: LogSummary[];
    orbitalLevels?: OrbitalLevels;
    orbitals?: FchkOrbitals;
    cavity?: SolventCavity;
    onProgress?: (fraction: number, message?: string) => void;
}

/**
 * The payload shape of one base-cell atom. The page keeps its own copy of the
 * fields a tooltip and a save need, so the mapping lives here once and is used
 * by both the file's expansion and a space-group alternative's.
 */
function crystalBasePayload(a: Atom) {
    return {
        element: a.element, x: a.x, y: a.y, z: a.z,
        index: a.index, occupancy: a.occupancy, baseIdx: a.baseIdx,
        siteLabel: a.siteLabel || undefined,
        asymIdx: a.asymIdx, symOp: a.symOp,
        constraint: a.constraint, disorderGroup: a.disorderGroup,
        uAniso: a.uAniso
    };
}

/** The payload shape of one asymmetric-unit site (label, frac, ops). */
function crystalAsymPayload(a: AsymAtom) {
    return {
        element: a.element, label: a.label, x: a.x, y: a.y, z: a.z,
        frac: a.frac, index: a.index, occupancy: a.occupancy,
        symmetryOps: a.symmetryOps, disorderGroup: a.disorderGroup,
        uAniso: a.uAniso
    };
}

function getNonce(): string {
    let text = '';
    const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    for (let i = 0; i < 32; i++) {
        text += possible.charAt(Math.floor(Math.random() * possible.length));
    }
    return text;
}

/**
 * A page for a webview whose renderer could not be read.
 *
 * It is static HTML: no script, no data payload, nothing that could fail a
 * second time. The message names the missing asset and the reason, so a broken
 * installation can be diagnosed from the panel itself.
 */
function webviewLoadErrorHtml(nonce: string, message: string): string {
    const safe = message.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
    return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'nonce-${nonce}';">
<title>Molecular Viewer</title>
<style nonce="${nonce}">
html,body{width:100%;height:100%;margin:0}
body{display:flex;align-items:center;justify-content:center;background:var(--vscode-editor-background,#1e1e1e);color:var(--vscode-editor-foreground,#ccc);font-family:var(--vscode-font-family,sans-serif)}
.box{max-width:640px;padding:20px 24px;border:1px solid var(--vscode-panel-border,#444);border-radius:4px;background:var(--vscode-editor-background,#1e1e1e)}
h2{font-size:15px;margin:0 0 10px}
p{font-size:12px;line-height:1.6;margin:0 0 10px;color:var(--vscode-descriptionForeground,#a9a9a9)}
pre{font-size:11px;line-height:1.5;white-space:pre-wrap;word-break:break-word;margin:0;padding:8px;border-radius:3px;background:rgba(128,128,128,0.12);user-select:text}
</style>
</head>
<body>
<div class="box">
<h2>The 3D viewer could not start</h2>
<p>The extension could not load the renderer it needs. Reinstall the extension (or restore its <b>media</b> folder), then reopen this file.</p>
<pre>${safe}</pre>
</div>
</body>
</html>`;
}

export async function renderViewerPage(extensionUri: vscode.Uri, options: ViewerPageOptions): Promise<string> {
    const {
        data, frames = [], optSteps, normalModes, logSource, routes, vmdPayload,
        autoBeta = false, optThresholds, trajectoryFrames, charges, summary,
        orbitalLevels, orbitals, cavity, onProgress
    } = options;
    const nonce = getNonce();
    onProgress?.(0, 'Reading the renderer');

    // The renderer and the page's own assets are shipped beside each other.
    // An installation that lost media/ (a partial update, a copy of the folder
    // without it) must say so on the page that opens rather than leave the
    // generic "Unknown error" panel behind and nothing else; the message names
    // the file that could not be read.
    let threeJsContent: string;
    let viewerCss: string;
    let viewerJs: string;
    const readAsset = async (name: string): Promise<string> => {
        try {
            return decodeText(await vscode.workspace.fs.readFile(
                vscode.Uri.joinPath(extensionUri, 'media', name)
            ));
        } catch (e: any) {
            throw new Error('media/' + name + ': ' + (e && e.message ? e.message : e));
        }
    };
    try {
        [threeJsContent, viewerCss, viewerJs] = await Promise.all([
            readAsset('three.min.js'),
            readAsset('viewer.css'),
            readAsset('viewer.js')
        ]);
    } catch (e: any) {
        return webviewLoadErrorHtml(nonce,
            'Could not read the viewer assets: ' + (e && e.message ? e.message : e));
    }
    onProgress?.(0.05, 'Preparing the structure');

    const atomColors = ATOM_COLORS;

    const atomData = data.atoms.map(a => {
        let color = atomColors[a.element] || '#FF1493';
        return {
            element: a.element, x: a.x, y: a.y, z: a.z,
            color: color,
            beta: a.beta,
            occupancy: a.occupancy,
            ghost: a.ghost === true ? true : undefined,
            baseIdx: a.baseIdx,
            cellI: a.cellI, cellJ: a.cellJ, cellK: a.cellK,
            uAniso: a.uAniso
        };
    });

    if (data.atomGroups && data.atomGroups.length > 0) {
        const atomColorOverride: { [key: number]: string } = {};
        for (const group of data.atomGroups) {
            for (const idx of group.indices) {
                atomColorOverride[idx] = group.color;
            }
        }
        for (let i = 0; i < atomData.length; i++) {
            if (atomColorOverride[i]) {
                atomData[i].color = atomColorOverride[i];
            }
        }
    }

    const bondData = data.bonds.map(b => ({
        atom1: b.atom1, atom2: b.atom2, order: b.order
    }));

    type FramePayloadData = {
        atoms: { element: string; x: number; y: number; z: number; color: string; ghost?: boolean }[];
        bonds: { atom1: number; atom2: number; order: number }[];
        stepLabel: string;
        charge?: number;
        multiplicity?: number;
        energy?: number;
    };
    const framesData: FramePayloadData[] = [];
    for (let fi = 0; fi < frames.length; fi++) {
        const f = frames[fi];
        framesData.push({
            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: atomColors[a.element] || '#FF1493', ghost: a.ghost === true ? true : undefined })),
            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
            stepLabel: f.stepLabel,
            // The electronic state the geometry was computed at. A log can run
            // several jobs with different values (`--Link1--`, an ORCA file
            // with more than one `* xyz` block), so every frame keeps its own;
            // the page adopts them on a switch and a save writes them out.
            charge: f.charge,
            multiplicity: f.multiplicity,
            // Single-point energy of this frame's structure (opt logs only), used
            // by the convergence charts to highlight the selected step and to
            // report its exact value. [2026-09-12 | Atreides-Jimmy]
            energy: (f as { energy?: number }).energy
        });
        // A 647-frame log is 650 000 atom objects; serialising them in one
        // synchronous loop holds the host for seconds. Only the view that
        // asked for a notification yields here, so a small page keeps the
        // old single-pass speed.
        if (onProgress && ((fi & 15) === 0 || fi === frames.length - 1)) {
            onProgress(0.05 + 0.45 * ((fi + 1) / frames.length),
                'Building frame ' + (fi + 1) + ' of ' + frames.length);
            await yieldToHost();
        }
    }

    const atomGroupsData = data.atomGroups ? data.atomGroups.map(g => ({
        colorId: g.colorId, color: g.color, indices: g.indices
    })) : [];

    const crystalData = data.crystal ? {
        a: data.crystal.a, b: data.crystal.b, c: data.crystal.c,
        alpha: data.crystal.alpha, beta: data.crystal.beta, gamma: data.crystal.gamma,
        latticeVectors: data.crystal.latticeVectors,
        spaceGroup: data.crystal.spaceGroup || '',
        symmetryOps: data.crystal.symmetryOps,
        symmetrySource: data.crystal.symmetrySource,
        symmetryCheck: data.crystal.symmetryCheck ? {
            declared: data.crystal.symmetryCheck.declared,
            symbol: data.crystal.symmetryCheck.symbol,
            settingLabel: data.crystal.symmetryCheck.settingLabel,
            number: data.crystal.symmetryCheck.number,
            fileOps: data.crystal.symmetryCheck.fileOps,
            tableOps: data.crystal.symmetryCheck.tableOps,
            fileAtoms: data.crystal.symmetryCheck.fileAtoms,
            tableAtoms: data.crystal.symmetryCheck.tableAtoms
        } : undefined,
        symmetryAlternative: data.crystal.symmetryAlternative ? {
            ops: data.crystal.symmetryAlternative.ops,
            symbol: data.crystal.symmetryAlternative.symbol,
            settingLabel: data.crystal.symmetryAlternative.settingLabel,
            baseAtoms: data.crystal.symmetryAlternative.baseAtoms.map(crystalBasePayload),
            asymAtoms: data.crystal.symmetryAlternative.asymAtoms.map(crystalAsymPayload)
        } : undefined,
        baseAtoms: data.crystal.baseAtoms.map(crystalBasePayload),
        baseBonds: data.crystal.baseBonds.map(b => ({
            atom1: b.atom1, atom2: b.atom2, order: b.order, shift: b.shift, valence: b.valence
        })),
        fileBonds: (data.crystal.fileBonds || []).map(b => ({
            atom1: b.atom1, atom2: b.atom2, order: b.order, shift: b.shift, valence: b.valence
        })),
        asymAtoms: (data.crystal.asymAtoms || []).map(crystalAsymPayload)
    } : null;

    const optStepsData = optSteps ? optSteps.map(s => ({
        step: s.step, energy: s.energy, maxForce: s.maxForce, rmsForce: s.rmsForce,
        maxDisplacement: s.maxDisplacement, rmsDisplacement: s.rmsDisplacement
    })) : null;

    const normalModesData = normalModes ? normalModes.map(m => ({
        index: m.index, frequency: m.frequency, symmetry: m.symmetry || '',
        reducedMass: m.reducedMass, forceConstant: m.forceConstant, irIntensity: m.irIntensity,
        displacements: m.displacements
    })) : null;

    const routesData = routes ? routes.map(r => ({
        raw: r.raw, hasOpt: r.hasOpt,
        keywords: r.keywords.map(k => ({ name: k.name, options: k.options }))
    })) : null;

    // The population analyses and the job results, as the file states them.
    // A charge value is `null` at an atom the table does not cover, which
    // the page reads as "no value" and never as zero.
    const chargesData = charges && charges.length > 0 ? charges.map(c => ({
        label: c.label,
        values: c.values,
        elements: c.elements,
        covered: c.covered,
        complete: c.complete,
        frame: c.frame
    })) : null;

    const summaryData = summary && summary.length > 0 ? summary.map(s => ({
        job: s.job,
        nAtoms: s.nAtoms,
        nBasis: s.nBasis,
        charge: s.charge,
        multiplicity: s.multiplicity,
        energies: s.energies.map(e => ({ label: e.label, value: e.value })),
        scfCycles: s.scfCycles,
        rmsGradient: s.rmsGradient,
        rmsGradientLabel: s.rmsGradientLabel,
        dipole: s.dipole ? { x: s.dipole.x, y: s.dipole.y, z: s.dipole.z, total: s.dipole.total } : null,
        pointGroup: s.pointGroup,
        cpuTime: s.cpuTime,
        elapsedTime: s.elapsedTime,
        termination: s.termination
    })) : null;

    const vmdData = vmdPayload && vmdPayload.molecules.length > 0 ? vmdPayload : null;

    // --- Orbitals and the solvation cavity ---
    //
    // The level diagram needs the energies, the occupations and nothing
    // else, so that is what travels: 1673 (or 4680) numbers, not the
    // coefficient block beside them. `canDraw` says whether a surface is
    // possible at all — a log states no coefficients, and a basis this
    // engine does not evaluate says so in `reason` rather than failing
    // when the button is pressed.
    let orbitalsData: { [key: string]: unknown } | null = null;
    if (orbitalLevels && orbitalLevels.spins.length > 0) {
        const drawable = !!orbitals && orbitals.basis.shellType.length > 0 &&
            orbitals.spins.some(s => !!s.layout && s.layout.uniform);
        let reason = '';
        if (!orbitals) {
            reason = 'This log prints orbital energies but not the coefficients: a surface needs the ' +
                'formatted checkpoint (.fchk) of the same job.';
        } else if (orbitals.basis.maxL > 3) {
            reason = 'This basis uses g functions (angular momentum ' + orbitals.basis.maxL +
                '); this viewer evaluates s, p, d and f only.';
        } else if (!orbitals.path) {
            reason = 'This checkpoint was read without a local file path, so its coefficient block ' +
                'cannot be read back for a surface; the levels above are its own.';
        } else if (!drawable) {
            reason = 'This checkpoint does not carry the orbital coefficients in a layout this viewer can read.';
        }
        orbitalsData = {
            closedShell: orbitalLevels.closedShell,
            nBasis: orbitals ? orbitals.nBasis : 0,
            model: orbitals ? orbitals.model : '',
            canDraw: drawable,
            reason: reason,
            spins: orbitalLevels.spins.map(s => ({
                label: s.label,
                energies: s.energies,
                occupations: s.occupations,
                occupied: s.occupied
            }))
        };
    }

    const cavityData = cavity ? {
        source: cavity.source,
        model: cavity.model || '',
        type: cavity.type || '',
        solvent: cavity.solvent || '',
        alpha: cavity.alpha,
        area: cavity.area,
        volume: cavity.volume,
        sphereCount: cavity.sphereCount !== undefined ? cavity.sphereCount : (cavity.spheres ? cavity.spheres.length : 0),
        spheres: cavity.spheres ? cavity.spheres.map(s => ({
            x: Math.round(s.x * 1e4) / 1e4,
            y: Math.round(s.y * 1e4) / 1e4,
            z: Math.round(s.z * 1e4) / 1e4,
            r: Math.round(s.r * 1e4) / 1e4
        })) : null
    } : null;

    const optThresholdsData = optThresholds ? {
        maxForce: optThresholds.maxForce, rmsForce: optThresholds.rmsForce,
        maxDisplacement: optThresholds.maxDisplacement, rmsDisplacement: optThresholds.rmsDisplacement
    } : null;

    // A trajectory frame carries no atom objects: the elements, colours and
    // bonds are the same in every frame, so only the flat positions and the
    // frame's own cell travel. The page pairs them with the template it
    // already has, which keeps a 500-frame trajectory a few numbers per
    // frame instead of a duplicate structure in every frame.
    const trajectoryFramesData = trajectoryFrames && trajectoryFrames.length > 0
        ? trajectoryFrames.map(f => ({
            stepLabel: f.stepLabel,
            xyz: f.xyz,
            cell: f.cell,
            time: f.time,
            step: f.step
        }))
        : null;

    const trajectoryData = data.trajectory ? {
        format: data.trajectory.format,
        source: data.trajectory.source,
        topology: data.trajectory.topology,
        framesInFile: data.trajectory.framesInFile,
        sampled: data.trajectory.sampled,
        truncated: data.trajectory.truncated
    } : null;

    onProgress?.(0.55, 'Building the page payload');
    // The payload is written into a <script> block, so a title or comment
    // that contains `</script>` would close the block early and leave the
    // page blank (and let the file inject markup). `<` is escaped to its
    // JSON spelling, which is valid JSON and stops the HTML parser from
    // seeing a tag; U+2028/U+2029 are escaped because they are line
    // terminators to a JavaScript parser.
    const jsonData = JSON.stringify({ atoms: atomData, bonds: bondData, title: data.title, atomColors: atomColors, filePath: data.filePath || '', frames: framesData, gjfMeta: data.gjfMeta || null, charge: data.charge, multiplicity: data.multiplicity, atomGroups: atomGroupsData, crystal: crystalData, optSteps: optStepsData, normalModes: normalModesData, optSource: logSource || null, optThresholds: optThresholdsData, routes: routesData, vmd: vmdData, charges: chargesData, summary: summaryData, orbitals: orbitalsData, cavity: cavityData, autoBeta: !!(autoBeta && !vmdData), hasExplicitBonds: !!data.hasExplicitBonds, warnings: data.warnings || [], mmParams: describeForceField(), trajectory: trajectoryData, trajectoryFrames: trajectoryFramesData })
        .replace(/</g, '\\u003c')
        .replace(/\u2028/g, '\\u2028')
        .replace(/\u2029/g, '\\u2029');
    onProgress?.(0.90, 'Assembling the page');

    // The bond perception module, as source. Both sides run this text —
    // the host calls the compiled function directly, the webview evaluates
    // the serialisation below — so the viewport's re-detection can never
    // drift away from what the parser produced.
    const bondPerceptionSource = BOND_PERCEPTION_SOURCE;

    // The bond-valence module, as source, for the same reason: the report
    // the page shows computes its sums with the module this host ships,
    // not a second copy kept in the template. It carries a reading aid
    // only — nothing it returns is written back to the structure.
    const bvsSource = BVS_SOURCE;

    onProgress?.(0.95, 'Assembling the page');

    return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'nonce-${nonce}'; style-src 'nonce-${nonce}';">
<title>Molecular Viewer</title>
<style nonce="${nonce}">
${viewerCss.trimEnd()}
</style>
</head>
<body>
<div id="toolbar">
<div class="trow" id="toolbar-modes">
<button class="tbtn active" data-mode="view">View</button>
<button class="tbtn" data-mode="bondLength">Bond Length</button>
<button class="tbtn" data-mode="bondAngle">Bond Angle</button>
<button class="tbtn" data-mode="dihedral">Dihedral</button>
<button class="tbtn" data-mode="rotateGroup">Rotate Group</button>
<button class="tbtn" data-mode="moveAtoms">Move Atoms</button>
<button class="tbtn" data-mode="mmOpt">MM Optimize</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="addAtom">Add Atom</button>
<button class="tbtn" data-mode="deleteAtom">Delete Atom</button>
<button class="tbtn" data-mode="replaceAtom">Replace Atom</button>
<button class="tbtn" data-mode="addH">Add H</button>
<button class="tbtn" id="import-btn">Import Structure</button>
</div>
<div class="trow" id="toolbar-actions">
<button class="tbtn" data-mode="neighbours">Neighbours</button>
<button class="tbtn" data-mode="interactions">Interactions</button>
<button class="tbtn" id="diff-btn">Diff</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="selectAtoms">Select Atoms</button>
<button class="tbtn" data-mode="boxSelect">Box Select</button>
<div class="tsep"></div>
<button class="tbtn" id="beta-btn">Color by Beta</button>
<button class="tbtn" id="egg-btn">Eggs</button>
<button class="tbtn" id="undo-btn" disabled>Undo</button>
<button class="tbtn" id="save-btn">Save As</button>
<button class="tbtn" id="snap-btn">PNG</button>
<button class="tbtn" id="reset-btn">Reset View</button>
<div class="tsep" id="frame-sep"></div>
<div id="frame-nav">
<button class="tbtn" id="prev-frame">◀</button>
<input type="number" id="frame-num" min="1" value="1">
<span id="frame-info">1/1</span>
<button class="tbtn" id="next-frame">▶</button>
<button class="tbtn" id="auto-play">⏵ Play</button>
</div>
<div class="tsep hidden" id="mol-sep"></div>
<div id="mol-nav">
<button class="tbtn" id="prev-mol">◀</button>
<span id="mol-label">1/1</span>
<button class="tbtn" id="next-mol">▶</button>
</div>
</div>
</div>
<div id="status-bar"><span id="mode-info">View Mode</span><span id="selection-info"></span></div>
<div id="container"><div id="loading">Loading 3D Viewer...</div><div id="mol-info"></div><div id="mm-panel"></div><div id="mm-reopen">MM Optimize</div><div id="egg-panel"></div><div id="diff-label"></div><div id="diff-label-right"></div><div id="diff-panel"></div><div id="diff-reopen">📊 Show Results</div><div id="axes-indicator"></div><div id="crystal-panel"><div class="pn-head"><h4>Supercell Bounds</h4><button class="pn-close" id="crystal-close" type="button" title="Minimize the crystal panel">&times;</button></div><div class="bnd-row bnd-axis"><label for="bnd-a-min">a</label><input type="number" id="bnd-a-min" value="0" step="0.1" title="lower bound along a" aria-label="lower bound along a"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-a-max" value="1" step="0.1" title="upper bound along a" aria-label="upper bound along a"></div><div class="bnd-row bnd-axis"><label for="bnd-b-min">b</label><input type="number" id="bnd-b-min" value="0" step="0.1" title="lower bound along b" aria-label="lower bound along b"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-b-max" value="1" step="0.1" title="upper bound along b" aria-label="upper bound along b"></div><div class="bnd-row bnd-axis"><label for="bnd-c-min">c</label><input type="number" id="bnd-c-min" value="0" step="0.1" title="lower bound along c" aria-label="lower bound along c"><span class="bnd-dash">&ndash;</span><input type="number" id="bnd-c-max" value="1" step="0.1" title="upper bound along c" aria-label="upper bound along c"></div><div class="bnd-presets"><button class="bnd-mini" id="bnd-p2" title="Eight unit cells, 0..2 along every axis">2×2×2</button><button class="bnd-mini" id="bnd-p3" title="Twenty-seven unit cells, 0..3 along every axis">3×3×3</button><button class="bnd-mini" id="bnd-p4" title="Sixty-four unit cells, 0..4 along every axis">4×4×4</button></div><div class="bnd-presets"><button class="bnd-mini" id="bnd-reset" title="Back to the single unit cell, 0..1 along every axis">Reset</button><button class="bnd-mini" id="bnd-origin" title="Keep the spans, move the low corner onto the lattice origin">Align Origin</button></div><div class="bnd-cell" id="bnd-cell-line"></div><div class="bnd-sym" id="bnd-sym-note"></div><div class="bnd-row"><label>Disorder</label><select id="bnd-disorder" title="Show the ordered part plus one disorder group; All shows every alternative the cell lists"><option value="">All</option></select></div><label class="bnd-occ off" id="bnd-occ-label" for="bnd-occ-normalize" title="When one disorder group is selected, count its sites at occupancy 1 in the formula, density and saved files, so the cell describes the ordered variant on screen. All keeps every occupancy the file states."><input type="checkbox" id="bnd-occ-normalize" checked disabled>Selected group fully occupied</label><div class="hk-title">Lattice Planes (h k l)</div><div class="bnd-row hk-row"><label>(h k l)</label><input type="number" id="hk-h" value="0" step="1" title="h"><input type="number" id="hk-k" value="0" step="1" title="k"><input type="number" id="hk-l" value="1" step="1" title="l"></div><div class="hk-val" id="hk-d1"></div><div class="bnd-row hk-row"><label>2nd</label><input type="number" id="hk2-h" value="1" step="1" title="h of the second plane"><input type="number" id="hk2-k" value="0" step="1" title="k of the second plane"><input type="number" id="hk2-l" value="0" step="1" title="l of the second plane"></div><div class="hk-val" id="hk-d2"></div><button id="hk-show" class="bnd-btn" title="Draw the (hkl) plane through the centre of the box">Show Plane</button><button id="xrd-open" class="bnd-btn" title="Calculate the powder X-ray diffraction pattern of the unit cell">Powder XRD Pattern...</button><button id="bvs-open" class="bnd-btn" title="Bond-valence sums per site: how well the coordination around each cation and anion matches its formal valence (read-only)">Bond Valence (BVS)...</button><button id="cav-open" class="bnd-btn" title="Probe the unit cell for closed cavities and open channels and report their volume (read-only)">Cavities...</button><button id="bnd-edit-cell" class="bnd-btn" title="Change the cell parameters; every site keeps its fractional coordinates">Edit Cell...</button><button id="bnd-remove-disorder" class="bnd-btn">Remove Disorder &lt;0.5</button></div><div id="view-panel"><div class="pn-head"><span class="vp-title">Auto Rotate</span><button class="pn-close" id="view-close" type="button" title="Minimize the View controls">&times;</button></div><button class="tbtn" id="spin-btn">Off</button><div class="vp-hint">Drag and release to spin, or hold an arrow key &mdash; the longer the hold, the faster. The next left press stops it.</div><div class="vp-title vp-sub">Display</div><select id="style-select" title="How the structure is drawn: ball and stick (the default), van der Waals spheres, sticks without the atom spheres, or coordination polyhedra for a crystal"><option value="ball">Ball &amp; Stick</option><option value="spacefill">Space-filling</option><option value="stick">Stick</option>${(function(){try{return JSON.parse(jsonData).crystal?'<option value="polyhedra">Polyhedra</option>':'';}catch(e){return '';}}())}</select><label class="vp-check off" id="ortep-line" title="Draw atoms that carry anisotropic displacement parameters as 50% probability ellipsoids"><input type="checkbox" id="ortep-chk" disabled> Thermal ellipsoids (50%)</label></div><div id="bake-panel"><div class="pn-head"><span class="vp-title">Edit Coordinates</span><button class="pn-close" id="bake-close" type="button" title="Minimize the Edit Coordinates card">&times;</button></div><div class="vp-hint" id="vp-bake-hint">Hold <b>Alt</b> and drag: left rotates the structure, right/middle moves it. The atom coordinates change, so Save writes the new geometry.</div></div><div id="select-panel"><div class="pn-head"><span class="sp-title">Select Atoms</span><button class="pn-close" id="select-close" type="button" title="Minimize the selection panel">&times;</button></div><input type="text" id="sel-panel-input" placeholder="e.g. 1 3-5 C H 8"><div class="sp-btns"><button class="mbtn mbtn-ok" id="sel-panel-ok">Select</button><button class="mbtn mbtn-cancel" id="sel-panel-clear">Clear</button></div></div><div id="neigh-panel"><div class="pn-head"><span class="ng-title">Neighbours</span><button class="pn-close" id="neigh-close" type="button" title="Minimize the neighbour list">&times;</button></div><div class="ng-line"><label for="ng-radius">Radius (A)</label><input type="number" id="ng-radius" value="3.5" min="0.1" max="12" step="0.1"><label class="ng-chk" for="ng-skip"><input type="checkbox" id="ng-skip" checked>skip bonded</label></div><div class="ng-line"><label>Anchor</label><span class="ng-anchor" id="ng-anchor">click an atom</span></div><div class="ng-summary" id="ng-summary"></div><div class="ng-list" id="ng-list"></div><div class="ng-btns"><button id="ng-select">Select listed</button></div><div class="ng-hint" id="ng-hint">Click a row to select that atom. Crystal images outside the drawn cell are listed with their cell offset.</div></div><div id="int-panel"><div class="pn-head"><span class="ng-title">Non-covalent Interactions</span><button class="pn-close" id="int-close" type="button" title="Minimize the interaction list">&times;</button></div><div class="ng-line ng-types"><label class="ng-chk" for="int-hb"><input type="checkbox" id="int-hb" checked>H-bond</label><label class="ng-chk" for="int-xb"><input type="checkbox" id="int-xb" checked>halogen</label><label class="ng-chk" for="int-cb"><input type="checkbox" id="int-cb" checked>chalcogen</label><label class="ng-chk" for="int-pn"><input type="checkbox" id="int-pn" checked>pnictogen</label><label class="ng-chk" for="int-pp"><input type="checkbox" id="int-pp" checked>&pi;-&pi;</label><label class="ng-chk" for="int-chpi"><input type="checkbox" id="int-chpi" checked>CH-&pi;</label></div><div class="ng-line"><label class="ng-chk" for="int-show"><input type="checkbox" id="int-show" checked>draw dashed lines</label></div><div class="ng-summary" id="int-summary"></div><div class="ng-list" id="int-list"></div><div class="ng-hint" id="int-hint">H-bond: d(H&middot;&middot;&middot;A) &lt; 2.5 &Aring; and D-H&middot;&middot;&middot;A &gt; 120&deg;. Halogen, chalcogen and pnictogen bonds: d(X&middot;&middot;&middot;A) &lt; vdW sum &minus; 0.1 &Aring; and D-X&middot;&middot;&middot;A &gt; 130&deg;; acceptors N, O, S, F, Cl, Br, I. &pi;-&pi;: aromatic 5&ndash;7 rings with a plane angle &le; 30&deg;, centroid distance &le; 4.5 &Aring; and a 2.8&ndash;3.8 &Aring; face gap (parallel/offset), or a &ge; 60&deg; T-shape with one ring edge 3.0&ndash;5.5 &Aring; over the other face. CH-&pi;: H&middot;&middot;&middot;centroid &lt; 3.0 &Aring;, C-H&middot;&middot;&middot;centroid &gt; 120&deg; and the H at least 1.0 &Aring; above the ring face.</div></div><div id="xrd-panel"><div class="xr-title"><span>Powder XRD</span><button class="xr-close" id="xrd-close" title="Close (Esc)">&times;</button></div><div class="xr-line"><label for="xrd-rad">Radiation</label><select id="xrd-rad" title="The X-ray wavelength the pattern is calculated for"></select><label for="xrd-2max">2&theta; max</label><select id="xrd-2max" class="xr-narrow" title="The highest diffraction angle in the pattern"></select></div><div class="xr-summary" id="xrd-summary"></div><canvas id="xrd-canvas" width="440" height="150"></canvas><div class="xr-note" id="xrd-peak-info"></div><div class="xr-note" id="xrd-exp-info"></div><div class="xr-note" id="xrd-note"></div><div class="xr-list" id="xrd-list"></div><div class="xr-btns"><button id="xrd-refresh" class="xr-secondary" title="Recalculate after the structure changed">Recalculate</button><button id="xrd-export" title="Save the drawn curve as a two-column text file">Export .xy</button></div><div class="xr-btns"><button id="xrd-import" class="xr-secondary" title="Overlay a measured two-column 2theta/intensity pattern">Import experimental .xy...</button><button id="xrd-clear-exp" class="xr-secondary" disabled title="Remove the experimental curve">Clear exp.</button></div><input type="file" id="xrd-file" accept=".xy,.txt,.csv,text/plain"><div class="xr-hint">Pattern of the unit cell on screen (the supercell bounds do not change it) with Cromer-Mann form factors, Lorenz-polarisation and a Caglioti Gaussian profile. Click a calculated peak to read its position against the experimental curve.</div></div><div id="cav-panel"><div class="cv-title"><span>Voids and Cavities</span><button class="cv-close" id="cav-close" title="Close (Esc)">&times;</button></div><div class="cv-line"><label for="cav-probe">Probe radius</label><input type="number" id="cav-probe" value="1.2" min="0" max="20" step="0.1" title="Radius of the spherical probe rolled through the cell"><span class="cv-unit">&Aring;</span><button id="cav-run" title="Scan the unit cell with this probe">Analyse</button><button id="cav-clear" class="cv-secondary" title="Remove the drawn surface">Clear</button></div><div class="cv-summary" id="cav-summary">Not analysed yet.</div><div class="cv-note" id="cav-note"></div><div class="cv-list" id="cav-list"></div><div class="cv-hint">A point is free when a probe of this radius fits there without touching an atom (van der Waals radii, periodic images included). The free space is split into components: one whose walk winds through a cell face is an open channel, otherwise it is a closed cavity. Click a row to draw that void as a translucent surface. Volumes are per unit cell, so a supercell shows the same numbers; pockets under 0.5 &Aring;&sup3; are ignored. r_free is the largest sphere that fits at the void's deepest point.</div></div><div id="rotate-panel"><div class="pn-head"><span class="rp-title">Rotate Group Around Axis</span><button class="pn-close" id="rotate-close" type="button" title="Minimize the rotate panel">&times;</button></div><div class="rp-field"><label>Group atoms (click atoms or type, e.g. H581 O380 C44 ...):</label><input type="text" id="rp-group-input" placeholder="e.g. H581 O380 C44 C43 C42 ..."></div><div class="rp-field"><label>Axis atoms (2+ atoms, single-bonded chain):</label><input type="text" id="rp-axis-input" placeholder="e.g. C32 C33 C43 C53  or  C33 C43"></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-apply">Apply &amp; Rotate</button><button class="mbtn mbtn-cancel" id="rp-clear">Clear</button></div><div class="rp-error" id="rp-error"></div><div id="rp-slider-wrap"><div class="rp-slider-row"><label>Angle:</label><input type="range" id="rp-slider" min="0" max="360" step="0.5" value="0"><input type="number" id="rp-angle-input" min="0" max="360" step="0.5" value="0" title="Type an angle (0-360) and press Enter"><span class="rp-angle-val" id="rp-angle-val">0°</span></div><div class="rp-btns"><button class="mbtn mbtn-ok" id="rp-done">Done</button><button class="mbtn mbtn-danger" id="rp-reset">Reset</button></div></div><div class="rp-hint">Axis atoms are kept fixed; the rest of the selected group rotates around the line through the axis atoms (0–360°).</div></div><div id="charge-labels"></div><div id="charge-panel"></div><div id="summary-panel"></div><div id="opt-panel"></div><div id="opt-reopen">📈 Convergence</div><div id="freq-panel"></div><div id="freq-reopen">🎵 Modes</div><div id="route-panel"></div><div id="route-reopen">⚙️ Route</div><div id="charge-reopen" class="pn-chip" title="Atomic charges, solvent cavity and molecular orbitals">⚛️ Properties</div><div id="summary-reopen" class="pn-chip" title="Open the job summary">📋 Summary</div><div id="crystal-reopen" class="pn-chip" title="Open the crystal panel">Crystal</div><div id="egg-reopen" class="pn-chip" title="Open the egg list">Eggs</div><div id="view-reopen" class="pn-chip" title="Open the View controls">View</div><div id="bake-reopen" class="pn-chip" title="Open the Edit Coordinates card">Coordinates</div><div id="select-reopen" class="pn-chip" title="Open the Select Atoms panel">Select Atoms</div><div id="neigh-reopen" class="pn-chip" title="Open the Neighbours panel">Neighbours</div><div id="int-reopen" class="pn-chip" title="Open the Non-covalent Interactions panel">Interactions</div><div id="rotate-reopen" class="pn-chip" title="Open the Rotate Group panel">Rotate Group</div><div id="xrd-reopen" class="pn-chip" title="Open the powder XRD panel">Powder XRD</div><div id="cav-reopen" class="pn-chip" title="Open the cavity and void analysis">Cavities</div></div>
<div id="error-msg"></div>
<div id="empty-msg"></div>
<div id="atom-tooltip"></div>
<div id="box-select-rect"></div>
<div id="modal-overlay"><div id="modal"></div></div>
<script nonce="${nonce}">
window.__MV_HOST__={data:${jsonData},bondPerception:${bondPerceptionSource},bvs:${bvsSource}};
</script>
<script nonce="${nonce}">
${threeJsContent}
</script>
<script nonce="${nonce}">
${viewerJs.trimEnd()}
</script>
</body>
</html>`;
}
