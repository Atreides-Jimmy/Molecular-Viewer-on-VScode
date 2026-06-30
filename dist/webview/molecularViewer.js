"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || (function () {
    var ownKeys = function(o) {
        ownKeys = Object.getOwnPropertyNames || function (o) {
            var ar = [];
            for (var k in o) if (Object.prototype.hasOwnProperty.call(o, k)) ar[ar.length] = k;
            return ar;
        };
        return ownKeys(o);
    };
    return function (mod) {
        if (mod && mod.__esModule) return mod;
        var result = {};
        if (mod != null) for (var k = ownKeys(mod), i = 0; i < k.length; i++) if (k[i] !== "default") __createBinding(result, mod, k[i]);
        __setModuleDefault(result, mod);
        return result;
    };
})();
Object.defineProperty(exports, "__esModule", { value: true });
exports.MolecularViewerProvider = void 0;
const vscode = __importStar(require("vscode"));
const path = __importStar(require("path"));
const index_1 = require("../parsers/index");
const bondDetector_1 = require("../parsers/bondDetector");
class MolecularViewerProvider {
    constructor(context) {
        this.context = context;
    }
    async openCustomDocument(uri, _openContext, _token) {
        const content = await vscode.workspace.fs.readFile(uri);
        const textContent = new TextDecoder().decode(content);
        const fileName = uri.path.split('/').pop() || 'unknown.xyz';
        const ext = fileName.toLowerCase().split('.').pop() || '';
        let data;
        let frames = [];
        let atomGroups;
        if (ext === 'tcl') {
            const tclResult = (0, index_1.parseTcl)(textContent);
            let sourceUri;
            if (tclResult.sourceFile) {
                if (path.isAbsolute(tclResult.sourceFile)) {
                    sourceUri = vscode.Uri.file(tclResult.sourceFile);
                }
                else {
                    const tclDir = path.dirname(uri.fsPath);
                    const resolvedPath = path.resolve(tclDir, tclResult.sourceFile);
                    sourceUri = vscode.Uri.file(resolvedPath);
                }
            }
            if (sourceUri) {
                try {
                    const sourceContent = await vscode.workspace.fs.readFile(sourceUri);
                    const sourceText = new TextDecoder().decode(sourceContent);
                    const sourceFileName = sourceUri.path.split('/').pop() || 'unknown.xyz';
                    const sourceExt = sourceFileName.toLowerCase().split('.').pop() || '';
                    if (sourceExt === 'log' || sourceExt === 'out') {
                        const logResult = (0, index_1.parseLogFile)(sourceText, sourceFileName);
                        frames = logResult.frames;
                        if (frames.length > 0) {
                            data = (0, bondDetector_1.ensureBonds)({
                                atoms: frames[0].atoms,
                                bonds: frames[0].bonds,
                                title: frames[0].title,
                                hasExplicitBonds: frames[0].hasExplicitBonds,
                                charge: frames[0].charge,
                                multiplicity: frames[0].multiplicity
                            });
                        }
                        else {
                            data = { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false };
                        }
                    }
                    else {
                        data = (0, index_1.parseFile)(sourceText, sourceFileName);
                        if (data.crystal) {
                            data.bonds = [];
                            data.crystal.baseBonds = [];
                        }
                        else {
                            data = (0, bondDetector_1.ensureBonds)(data);
                        }
                    }
                    data.filePath = sourceUri.fsPath;
                }
                catch {
                    data = { atoms: [], bonds: [], title: 'Failed to load source file: ' + tclResult.sourceFile, hasExplicitBonds: false };
                }
            }
            else {
                data = { atoms: [], bonds: [], title: 'No source file specified in TCL', hasExplicitBonds: false };
            }
            atomGroups = tclResult.groups.length > 0 ? tclResult.groups : undefined;
            data.atomGroups = atomGroups;
        }
        else if (ext === 'log' || ext === 'out') {
            const logResult = (0, index_1.parseLogFile)(textContent, fileName);
            frames = logResult.frames;
            if (frames.length > 0) {
                data = (0, bondDetector_1.ensureBonds)({
                    atoms: frames[0].atoms,
                    bonds: frames[0].bonds,
                    title: frames[0].title,
                    hasExplicitBonds: frames[0].hasExplicitBonds,
                    charge: frames[0].charge,
                    multiplicity: frames[0].multiplicity
                });
            }
            else {
                data = { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false };
            }
        }
        else {
            data = (0, index_1.parseFile)(textContent, fileName);
            if (data.crystal) {
                data.bonds = [];
                data.crystal.baseBonds = [];
            }
            else {
                data = (0, bondDetector_1.ensureBonds)(data);
            }
        }
        if (ext !== 'tcl') {
            data.filePath = uri.fsPath;
        }
        return new MolecularDocument(uri, data, frames);
    }
    async resolveCustomEditor(document, webviewPanel, _token) {
        webviewPanel.webview.options = {
            enableScripts: true,
        };
        webviewPanel.webview.html = await this.getHtmlForWebview(webviewPanel.webview, document.data, document.frames);
        webviewPanel.webview.onDidReceiveMessage(async (message) => {
            switch (message.command) {
                case 'saveFile':
                    try {
                        const srcPath = message.filePath || '';
                        const srcDir = srcPath ? srcPath.substring(0, srcPath.replace(/\\/g, '/').lastIndexOf('/')) : '';
                        const defaultName = message.suggestedName || 'molecule.xyz';
                        const defaultUri = srcDir ? vscode.Uri.file(srcDir + '/' + defaultName) : vscode.Uri.file(defaultName);
                        const uri = await vscode.window.showSaveDialog({
                            defaultUri: defaultUri,
                            filters: {
                                'XYZ Files': ['xyz'],
                                'Gaussian Input': ['gjf'],
                                'CIF Files': ['cif'],
                                'VASP POSCAR': ['vasp', 'poscar', 'contcar'],
                                'Gaussian Cube': ['cube'],
                                'Turbomole Coord': ['coord'],
                                'ORCA Input': ['inp'],
                                'MOL2': ['mol2'],
                                'MDL Mol': ['mol'],
                                'PDB': ['pdb'],
                                'MOPAC': ['mop'],
                                'All Files': ['*']
                            }
                        });
                        if (uri) {
                            await vscode.workspace.fs.writeFile(uri, new TextEncoder().encode(message.content));
                            vscode.window.showInformationMessage('Saved: ' + uri.fsPath);
                        }
                    }
                    catch (e) {
                        vscode.window.showErrorMessage('Save failed: ' + (e.message || e));
                    }
                    break;
                case 'diffFile':
                    try {
                        const result = await vscode.window.showOpenDialog({
                            canSelectMany: false,
                            openLabel: 'Select File to Compare',
                            filters: {
                                'Molecular Files': ['gjf', 'xyz', 'mol', 'sdf', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol2', 'log', 'out', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'tcl', 'cif', 'vasp', 'poscar', 'contcar', 'cube', 'vesta'],
                                'All Files': ['*']
                            }
                        });
                        if (!result || result.length === 0) {
                            webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                            break;
                        }
                        const diffUri = result[0];
                        const diffContent = await vscode.workspace.fs.readFile(diffUri);
                        const diffText = new TextDecoder().decode(diffContent);
                        const diffFileName = diffUri.path.split('/').pop() || 'unknown.xyz';
                        const diffExt = diffFileName.toLowerCase().split('.').pop() || '';
                        let diffData;
                        let diffFrames = [];
                        if (diffExt === 'tcl') {
                            const tclResult = (0, index_1.parseTcl)(diffText);
                            let sourceUri;
                            if (tclResult.sourceFile) {
                                if (path.isAbsolute(tclResult.sourceFile)) {
                                    sourceUri = vscode.Uri.file(tclResult.sourceFile);
                                }
                                else {
                                    const tclDir = path.dirname(diffUri.fsPath);
                                    sourceUri = vscode.Uri.file(path.resolve(tclDir, tclResult.sourceFile));
                                }
                            }
                            if (sourceUri) {
                                const srcContent = await vscode.workspace.fs.readFile(sourceUri);
                                const srcText = new TextDecoder().decode(srcContent);
                                const srcFileName = sourceUri.path.split('/').pop() || 'unknown.xyz';
                                const srcExt = srcFileName.toLowerCase().split('.').pop() || '';
                                if (srcExt === 'log' || srcExt === 'out') {
                                    const logResult = (0, index_1.parseLogFile)(srcText, srcFileName);
                                    diffFrames = logResult.frames;
                                    if (diffFrames.length > 0) {
                                        diffData = (0, bondDetector_1.ensureBonds)({
                                            atoms: diffFrames[0].atoms, bonds: diffFrames[0].bonds,
                                            title: diffFrames[0].title, hasExplicitBonds: diffFrames[0].hasExplicitBonds,
                                            charge: diffFrames[0].charge, multiplicity: diffFrames[0].multiplicity
                                        });
                                    }
                                    else {
                                        diffData = { atoms: [], bonds: [], title: 'No structures', hasExplicitBonds: false };
                                    }
                                }
                                else {
                                    diffData = (0, index_1.parseFile)(srcText, srcFileName);
                                    if (diffData.crystal) {
                                        diffData.bonds = [];
                                        diffData.crystal.baseBonds = [];
                                    }
                                    else {
                                        diffData = (0, bondDetector_1.ensureBonds)(diffData);
                                    }
                                }
                                diffData.filePath = sourceUri.fsPath;
                            }
                            else {
                                diffData = { atoms: [], bonds: [], title: 'No source', hasExplicitBonds: false };
                            }
                            diffData.atomGroups = tclResult.groups.length > 0 ? tclResult.groups : undefined;
                        }
                        else if (diffExt === 'log' || diffExt === 'out') {
                            const logResult = (0, index_1.parseLogFile)(diffText, diffFileName);
                            diffFrames = logResult.frames;
                            if (diffFrames.length > 0) {
                                diffData = (0, bondDetector_1.ensureBonds)({
                                    atoms: diffFrames[0].atoms, bonds: diffFrames[0].bonds,
                                    title: diffFrames[0].title, hasExplicitBonds: diffFrames[0].hasExplicitBonds,
                                    charge: diffFrames[0].charge, multiplicity: diffFrames[0].multiplicity
                                });
                            }
                            else {
                                diffData = { atoms: [], bonds: [], title: 'No structures', hasExplicitBonds: false };
                            }
                        }
                        else {
                            diffData = (0, index_1.parseFile)(diffText, diffFileName);
                            if (diffData.crystal) {
                                diffData.bonds = [];
                                diffData.crystal.baseBonds = [];
                            }
                            else {
                                diffData = (0, bondDetector_1.ensureBonds)(diffData);
                            }
                            diffData.filePath = diffUri.fsPath;
                        }
                        let selectedFrameIdx = 0;
                        if (diffFrames.length > 1) {
                            const items = diffFrames.map((f, i) => ({
                                label: `Frame ${i + 1}: ${f.stepLabel || ''}`,
                                description: `${f.atoms.length} atoms`,
                                index: i
                            }));
                            const picked = await vscode.window.showQuickPick(items, {
                                placeHolder: `Select frame to compare (${diffFrames.length} frames found)`,
                                canPickMany: false
                            });
                            if (!picked) {
                                webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                                break;
                            }
                            selectedFrameIdx = picked.index;
                            const sf = diffFrames[selectedFrameIdx];
                            diffData = (0, bondDetector_1.ensureBonds)({
                                atoms: sf.atoms, bonds: sf.bonds, title: sf.title,
                                hasExplicitBonds: sf.hasExplicitBonds, charge: sf.charge, multiplicity: sf.multiplicity
                            });
                        }
                        const diffAtomColors = {
                            H: '#FFFFFF', C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050',
                            S: '#FFFF30', Cl: '#1FF01F', Br: '#A62929', I: '#940094', P: '#FF8000',
                            Na: '#AB5CF2', Mg: '#8AFF00', K: '#8F40D4', Ca: '#3DFF00', Fe: '#E06633',
                            Cu: '#C88033', Zn: '#7D80B0', Ag: '#C0C0C0', Au: '#FFD123', Si: '#F0C8A0',
                            B: '#FFB5B5', Li: '#CC80FF', Be: '#C2FF00', Al: '#BFA6A6'
                        };
                        const diffAtoms = diffData.atoms.map(a => ({
                            element: a.element, x: a.x, y: a.y, z: a.z,
                            color: diffAtomColors[a.element] || '#FF1493'
                        }));
                        if (diffData.atomGroups && diffData.atomGroups.length > 0) {
                            const override = {};
                            for (const g of diffData.atomGroups)
                                for (const idx of g.indices)
                                    override[idx] = g.color;
                            for (let i = 0; i < diffAtoms.length; i++)
                                if (override[i])
                                    diffAtoms[i].color = override[i];
                        }
                        const diffBonds = diffData.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order }));
                        const diffFramesData = diffFrames.map(f => ({
                            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: diffAtomColors[a.element] || '#FF1493' })),
                            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
                            stepLabel: f.stepLabel
                        }));
                        webviewPanel.webview.postMessage({
                            command: 'diffResult',
                            cancelled: false,
                            fileName: diffFileName,
                            atoms: diffAtoms,
                            bonds: diffBonds,
                            frames: diffFramesData,
                            selectedFrame: selectedFrameIdx,
                            title: diffData.title
                        });
                    }
                    catch (e) {
                        vscode.window.showErrorMessage('Diff failed: ' + (e.message || e));
                        webviewPanel.webview.postMessage({ command: 'diffResult', cancelled: true });
                    }
                    break;
                case 'info':
                    vscode.window.showInformationMessage(message.text);
                    break;
                case 'error':
                    vscode.window.showErrorMessage(message.text);
                    break;
            }
        });
    }
    async getHtmlForWebview(webview, data, frames = []) {
        const nonce = getNonce();
        const threeJsBytes = await vscode.workspace.fs.readFile(vscode.Uri.joinPath(this.context.extensionUri, 'media', 'three.min.js'));
        const threeJsContent = new TextDecoder().decode(threeJsBytes);
        const atomColors = {
            H: '#FFFFFF', He: '#D9FFFF', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5',
            C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050', Ne: '#B3E3F5',
            Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6', Si: '#F0C8A0', P: '#FF8000',
            S: '#FFFF30', Cl: '#1FF01F', Ar: '#80D1E3', K: '#8F40D4', Ca: '#3DFF00',
            Sc: '#E6E6E6', Ti: '#BFC2C7', V: '#A6A6AB', Cr: '#8A99C7', Mn: '#9C7AC7',
            Fe: '#E06633', Co: '#F090A0', Ni: '#50D050', Cu: '#C88033', Zn: '#7D80B0',
            Ga: '#C28F8F', Ge: '#668F8F', As: '#BD80E3', Se: '#FFA100', Br: '#A62929',
            Kr: '#5CB8D1', Rb: '#702EB0', Sr: '#00FF00', Y: '#94FFFF', Zr: '#94E0E0',
            Nb: '#73C2C9', Mo: '#54B5B5', Tc: '#3B9E9E', Ru: '#248F8F', Rh: '#0A7D8C',
            Pd: '#006985', Ag: '#C0C0C0', Cd: '#FFD98F', In: '#A67573', Sn: '#668080',
            Sb: '#9E63B5', Te: '#D47A00', I: '#940094', Xe: '#429EB0', Cs: '#57178F',
            Ba: '#00C900', La: '#70D4FF', Ce: '#FFFFC7', Pr: '#D9FFC7', Nd: '#C7FFC7',
            Pm: '#A3FFC7', Sm: '#8FFFC7', Eu: '#61FFC7', Gd: '#45FFC7', Tb: '#30FFC7',
            Dy: '#1FFFC7', Ho: '#00FF9C', Er: '#00E675', Tm: '#00D452', Yb: '#00BF38',
            Lu: '#00AB24', Hf: '#4DC2FF', Ta: '#4DA6FF', W: '#2194D6', Re: '#267DAB',
            Os: '#266696', Ir: '#175487', Pt: '#D0D0E0', Au: '#FFD123', Hg: '#B8B8D0',
            Tl: '#A6544D', Pb: '#575961', Bi: '#9E4FB5', Po: '#AB5C00', At: '#754F45',
            Rn: '#428296', Fr: '#420066', Ra: '#007D00', Ac: '#70ABFA', Th: '#00BAFF',
            Pa: '#00A1FF', U: '#008FFF', Np: '#0080FF', Pu: '#006BFF', Am: '#545CF2',
            Cm: '#785CE3', Bk: '#8A4FE3', Cf: '#A136D4', Es: '#B31FD4', Fm: '#B31FBA',
            Md: '#B30DA6', No: '#BD0D87', Lr: '#C70066', Rf: '#CC0059', Db: '#D9004F',
            Sg: '#E00045', Bh: '#E6002E', Hs: '#EB0026'
        };
        const atomData = data.atoms.map(a => {
            let color = atomColors[a.element] || '#FF1493';
            return {
                element: a.element, x: a.x, y: a.y, z: a.z,
                color: color,
                occupancy: a.occupancy,
                baseIdx: a.baseIdx,
                cellI: a.cellI, cellJ: a.cellJ, cellK: a.cellK
            };
        });
        if (data.atomGroups && data.atomGroups.length > 0) {
            const atomColorOverride = {};
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
        const framesData = frames.map(f => ({
            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: atomColors[a.element] || '#FF1493' })),
            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
            stepLabel: f.stepLabel
        }));
        const atomGroupsData = data.atomGroups ? data.atomGroups.map(g => ({
            colorId: g.colorId, color: g.color, indices: g.indices
        })) : [];
        const crystalData = data.crystal ? {
            a: data.crystal.a, b: data.crystal.b, c: data.crystal.c,
            alpha: data.crystal.alpha, beta: data.crystal.beta, gamma: data.crystal.gamma,
            latticeVectors: data.crystal.latticeVectors,
            spaceGroup: data.crystal.spaceGroup || '',
            symmetryOps: data.crystal.symmetryOps,
            baseAtoms: data.crystal.baseAtoms.map(a => ({
                element: a.element, x: a.x, y: a.y, z: a.z,
                index: a.index, occupancy: a.occupancy, baseIdx: a.baseIdx
            })),
            baseBonds: data.crystal.baseBonds.map(b => ({
                atom1: b.atom1, atom2: b.atom2, order: b.order
            }))
        } : null;
        const jsonData = JSON.stringify({ atoms: atomData, bonds: bondData, title: data.title, atomColors: atomColors, filePath: data.filePath || '', frames: framesData, gjfMeta: data.gjfMeta || null, charge: data.charge, multiplicity: data.multiplicity, atomGroups: atomGroupsData, crystal: crystalData });
        return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'nonce-${nonce}'; style-src 'nonce-${nonce}';">
<title>Molecular Viewer</title>
<style nonce="${nonce}">
*{margin:0;padding:0;box-sizing:border-box}
html{width:100%;height:100%;overflow:hidden}
body{width:100%;height:100%;overflow:hidden;display:flex;flex-direction:column;background:var(--vscode-editor-background,#1e1e1e);font-family:var(--vscode-font-family,sans-serif);color:var(--vscode-editor-foreground,#ccc)}
#toolbar{height:36px;flex-shrink:0;background:var(--vscode-editor-background,#1e1e1e);border-bottom:1px solid var(--vscode-panel-border,#444);display:flex;align-items:center;padding:0 8px;z-index:20;gap:2px}
.tbtn{background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff);border:1px solid var(--vscode-panel-border,#444);padding:3px 10px;border-radius:3px;cursor:pointer;font-size:11px;white-space:nowrap}
.tbtn:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
.tbtn.active{background:var(--vscode-button-background,#0e639c);border-color:var(--vscode-button-background,#0e639c)}
.tsep{width:1px;height:20px;background:var(--vscode-panel-border,#444);margin:0 4px}
#status-bar{height:24px;flex-shrink:0;background:var(--vscode-statusBar-background,#007acc);color:var(--vscode-statusBar-foreground,#fff);display:flex;align-items:center;padding:0 10px;font-size:11px;z-index:20;gap:12px}
#container{flex:1;position:relative;overflow:hidden;min-height:0}
#mol-info{position:absolute;top:8px;left:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.55);padding:6px 10px;border-radius:4px;z-index:25;pointer-events:none;line-height:1.6}
#axes-indicator{position:absolute;bottom:12px;left:12px;width:90px;height:90px;z-index:25;pointer-events:none}
#axes-indicator svg{width:100%;height:100%}
#crystal-panel{position:absolute;top:8px;right:8px;background:rgba(0,0,0,0.6);border:1px solid var(--vscode-panel-border,#444);border-radius:4px;padding:8px 10px;font-size:11px;color:var(--vscode-editor-foreground,#ccc);z-index:25;display:none;min-width:180px}
#crystal-panel h4{margin:0 0 6px 0;font-size:11px;font-weight:600}
#crystal-panel .bnd-row{display:flex;align-items:center;gap:4px;margin:3px 0}
#crystal-panel .bnd-row label{width:42px;text-align:right}
#crystal-panel .bnd-row input{width:48px;background:var(--vscode-input-background,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border:1px solid var(--vscode-input-border,#444);border-radius:2px;padding:1px 3px;font-size:11px;text-align:center}
#crystal-panel .bnd-btn{margin-top:6px;width:100%;background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border:none;border-radius:3px;padding:4px 6px;font-size:11px;cursor:pointer}
#crystal-panel .bnd-btn:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
canvas{display:block}
#atom-tooltip{position:absolute;display:none;color:var(--vscode-editor-foreground,#ccc);font-size:12px;background:var(--vscode-editor-background,#1e1e1e);padding:4px 8px;border-radius:3px;border:1px solid var(--vscode-panel-border,#444);pointer-events:none;z-index:30}
#modal-overlay{position:fixed;top:0;left:0;right:0;bottom:0;background:rgba(0,0,0,0.5);z-index:100;display:none;align-items:center;justify-content:center}
#modal-overlay.show{display:flex}
#modal{background:var(--vscode-editor-background,#1e1e1e);border:1px solid var(--vscode-panel-border,#444);border-radius:6px;padding:16px 20px;min-width:320px;max-width:420px;box-shadow:0 8px 30px rgba(0,0,0,0.5)}
#modal h3{font-size:14px;margin-bottom:10px;color:var(--vscode-editor-foreground,#ccc)}
#modal label{font-size:12px;display:block;margin:6px 0 2px}
#modal input[type=number],#modal select{width:100%;padding:4px 8px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:12px}
#modal input[type=range]{width:100%;margin:4px 0}
#modal .modal-row{display:flex;gap:8px;align-items:center;margin:4px 0}
#modal .modal-row label{margin:0;white-space:nowrap;min-width:60px}
#modal .modal-row input,#modal .modal-row select{flex:1}
#modal .modal-btns{display:flex;justify-content:flex-end;gap:8px;margin-top:14px}
#modal .mbtn{padding:5px 16px;border-radius:3px;cursor:pointer;font-size:12px;border:1px solid var(--vscode-panel-border,#444)}
#modal .mbtn-ok{background:var(--vscode-button-background,#0e639c);color:var(--vscode-button-foreground,#fff);border-color:var(--vscode-button-background,#0e639c)}
#modal .mbtn-ok:hover{background:var(--vscode-button-hoverBackground,#1177bb)}
#modal .mbtn-cancel{background:var(--vscode-button-secondaryBackground,#3a3d41);color:var(--vscode-button-secondaryForeground,#fff)}
#modal .mbtn-cancel:hover{background:var(--vscode-button-secondaryHoverBackground,#45494e)}
#modal .mbtn-danger{background:#c33;border-color:#c33;color:#fff}
#modal .current-val{font-size:13px;color:var(--vscode-descriptionForeground,#999);margin-bottom:6px}
#loading{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:var(--vscode-editor-foreground,#ccc);font-size:14px}
#error-msg{display:none;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);color:#f66;padding:20px;font-size:13px;text-align:center;max-width:80%;z-index:50}
.hidden{display:none!important}
#frame-nav{display:none;align-items:center;gap:2px}
#frame-nav.show{display:flex}
#frame-num{width:40px;padding:2px 4px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:11px;text-align:center}
#frame-info{color:var(--vscode-statusBar-foreground,#fff);font-size:11px;padding:0 4px;white-space:nowrap}
#auto-play.playing{background:var(--vscode-button-background,#0e639c);border-color:var(--vscode-button-background,#0e639c)}
#diff-panel{display:none;position:absolute;bottom:40px;right:8px;color:var(--vscode-editor-foreground,#ccc);font-size:11px;background:rgba(0,0,0,0.65);padding:8px 12px;border-radius:4px;z-index:25;max-width:40%;max-height:60%;overflow-y:auto;pointer-events:auto;line-height:1.5}
#diff-panel.show{display:block}
#diff-panel h4{font-size:12px;margin-bottom:4px;color:var(--vscode-textLink-foreground,#3794ff)}
#diff-panel .diff-row{margin:2px 0}
#diff-panel .diff-close{position:absolute;top:4px;right:8px;cursor:pointer;color:var(--vscode-descriptionForeground,#999);font-size:14px}
#diff-panel .diff-close:hover{color:var(--vscode-errorForeground,#f66)}
#diff-reopen{display:none;position:absolute;bottom:40px;right:8px;z-index:26;padding:4px 10px;border-radius:4px;border:1px solid var(--vscode-button-border,#555);background:rgba(0,0,0,0.7);color:var(--vscode-editor-foreground,#ccc);font-size:11px;cursor:pointer;pointer-events:auto}
#diff-reopen.show{display:block}
#diff-reopen:hover{background:var(--vscode-button-background,#0e639c)}
#diff-label,#diff-label-right{display:none;position:absolute;top:8px;color:var(--vscode-editor-foreground,#ccc);font-size:13px;font-weight:bold;background:rgba(0,0,0,0.55);padding:4px 10px;border-radius:4px;z-index:25;pointer-events:none}
#diff-label.show,#diff-label-right.show{display:block}
#diff-label{left:8px}
#diff-label-right{right:8px}
</style>
</head>
<body>
<div id="toolbar">
<button class="tbtn active" data-mode="view">View</button>
<button class="tbtn" data-mode="bondLength">Bond Length</button>
<button class="tbtn" data-mode="bondAngle">Bond Angle</button>
<button class="tbtn" data-mode="dihedral">Dihedral</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="addAtom">Add Atom</button>
<button class="tbtn" data-mode="deleteAtom">Delete Atom</button>
<div class="tsep"></div>
<button class="tbtn" data-mode="selectAtoms">Select Atoms</button>
<div class="tsep"></div>
<button class="tbtn" id="diff-btn">Diff</button>
<button class="tbtn" id="save-btn">Save As</button>
<button class="tbtn" id="reset-btn">Reset View</button>
<div class="tsep" id="frame-sep"></div>
<div id="frame-nav">
<button class="tbtn" id="prev-frame">◀</button>
<input type="number" id="frame-num" min="1" value="1">
<span id="frame-info">1/1</span>
<button class="tbtn" id="next-frame">▶</button>
<button class="tbtn" id="auto-play">⏵ Play</button>
</div>
</div>
<div id="status-bar"><span id="mode-info">View Mode</span><span id="selection-info"></span></div>
<div id="container"><div id="loading">Loading 3D Viewer...</div><div id="mol-info"></div><div id="diff-label"></div><div id="diff-label-right"></div><div id="diff-panel"></div><div id="diff-reopen">📊 Show Results</div><div id="axes-indicator"></div><div id="crystal-panel"><h4>Supercell Bounds</h4><div class="bnd-row"><label>a min</label><input type="number" id="bnd-a-min" value="0" step="0.1"></div><div class="bnd-row"><label>a max</label><input type="number" id="bnd-a-max" value="1" step="0.1"></div><div class="bnd-row"><label>b min</label><input type="number" id="bnd-b-min" value="0" step="0.1"></div><div class="bnd-row"><label>b max</label><input type="number" id="bnd-b-max" value="1" step="0.1"></div><div class="bnd-row"><label>c min</label><input type="number" id="bnd-c-min" value="0" step="0.1"></div><div class="bnd-row"><label>c max</label><input type="number" id="bnd-c-max" value="1" step="0.1"></div><button id="bnd-remove-disorder" class="bnd-btn">Remove Disorder &lt;0.5</button></div></div>
<div id="error-msg"></div>
<div id="atom-tooltip"></div>
<div id="modal-overlay"><div id="modal"></div></div>
<script nonce="${nonce}">
${threeJsContent}
</script>
<script nonce="${nonce}">
(function(){
try{
var MD=${jsonData};
var AN={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Ne:10,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Sc:21,Ti:22,V:23,Cr:24,Mn:25,Fe:26,Co:27,Ni:28,Cu:29,Zn:30,Ga:31,Ge:32,As:33,Se:34,Br:35,Kr:36,Rb:37,Sr:38,Y:39,Zr:40,Nb:41,Mo:42,Tc:43,Ru:44,Rh:45,Pd:46,Ag:47,Cd:48,In:49,Sn:50,Sb:51,Te:52,I:53,Xe:54,Cs:55,Ba:56,La:57,Ce:58,Pr:59,Nd:60,Pm:61,Sm:62,Eu:63,Gd:64,Tb:65,Dy:66,Ho:67,Er:68,Tm:69,Yb:70,Lu:71,Hf:72,Ta:73,W:74,Re:75,Os:76,Ir:77,Pt:78,Au:79,Hg:80,Tl:81,Pb:82,Bi:83,Po:84,At:85,Rn:86};
var CRY=MD.crystal||null;
var boundary={aMin:0,aMax:1,bMin:0,bMax:1,cMin:0,cMax:1};
var occTextureCache={};
var CRY_INV=null;
function computeCrystalInv(){
    if(!CRY){CRY_INV=null;return}
    var lv=CRY.latticeVectors;
    var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
    if(Math.abs(det)<1e-12){CRY_INV=null;return}
    var invDet=1/det;
    CRY_INV=[
        [(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])*invDet,(lv[0][2]*lv[2][1]-lv[0][1]*lv[2][2])*invDet,(lv[0][1]*lv[1][2]-lv[0][2]*lv[1][1])*invDet],
        [(lv[1][2]*lv[2][0]-lv[1][0]*lv[2][2])*invDet,(lv[0][0]*lv[2][2]-lv[0][2]*lv[2][0])*invDet,(lv[0][2]*lv[1][0]-lv[0][0]*lv[1][2])*invDet],
        [(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0])*invDet,(lv[0][1]*lv[2][0]-lv[0][0]*lv[2][1])*invDet,(lv[0][0]*lv[1][1]-lv[0][1]*lv[1][0])*invDet]
    ];
}
computeCrystalInv();
function updateMolInfo(){
    var infoEl=document.getElementById('mol-info');
    if(!infoEl)return;
    var nAtoms=MD.atoms.length;
    var chrg=MD.charge!=null?MD.charge:'-';
    var mult=MD.multiplicity!=null?MD.multiplicity:'-';
    var nElectrons=0;
    MD.atoms.forEach(function(a){nElectrons+=(AN[a.element]||0)});
    if(typeof chrg==='number')nElectrons-=chrg;
    var html='Atoms: '+nAtoms+'<br>Charge: '+chrg+'<br>Electrons: '+nElectrons+'<br>Multiplicity: '+mult;
    if(CRY){
        var sa=boundary.aMax-boundary.aMin,sb=boundary.bMax-boundary.bMin,sc=boundary.cMax-boundary.cMin;
        html+='<br>Cell: '+CRY.a.toFixed(3)+' × '+CRY.b.toFixed(3)+' × '+CRY.c.toFixed(3);
        html+='<br>Angles: '+CRY.alpha.toFixed(1)+'°, '+CRY.beta.toFixed(1)+'°, '+CRY.gamma.toFixed(1)+'°';
        if(CRY.spaceGroup)html+='<br>SG: '+CRY.spaceGroup;
        html+='<br>Supercell: '+sa.toFixed(2)+'×'+sb.toFixed(2)+'×'+sc.toFixed(2);
    }
    infoEl.innerHTML=html;
}
updateMolInfo();
var container=document.getElementById('container');
var loadingEl=document.getElementById('loading');
var errorEl=document.getElementById('error-msg');
var tooltipEl=document.getElementById('atom-tooltip');
var modeInfoEl=document.getElementById('mode-info');
var selInfoEl=document.getElementById('selection-info');
var modalOverlay=document.getElementById('modal-overlay');
var modalEl=document.getElementById('modal');
var vscodeApi=acquireVsCodeApi();

function showError(msg){if(loadingEl)loadingEl.style.display='none';errorEl.style.display='block';errorEl.textContent=msg}

if(typeof THREE==='undefined'){showError('Three.js library failed to load. Please reinstall the extension.');return}

var cw=container.clientWidth||window.innerWidth;
var ch=container.clientHeight||(window.innerHeight-60);
if(ch<1)ch=window.innerHeight-60;
if(cw<1)cw=window.innerWidth;

var scene=new THREE.Scene();
scene.background=new THREE.Color(0x1e1e1e);
var camera=new THREE.PerspectiveCamera(60,cw/ch,0.1,1000);
var renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setSize(cw,ch);
renderer.setPixelRatio(window.devicePixelRatio);
container.appendChild(renderer.domElement);
if(loadingEl)loadingEl.style.display='none';

scene.add(new THREE.AmbientLight(0x404040,1.5));
var dl1=new THREE.DirectionalLight(0xffffff,0.8);dl1.position.set(5,10,7);scene.add(dl1);
var dl2=new THREE.DirectionalLight(0xffffff,0.4);dl2.position.set(-5,-3,-5);scene.add(dl2);

var pivotGroup=new THREE.Group();scene.add(pivotGroup);
var moleculeGroup=new THREE.Group();pivotGroup.add(moleculeGroup);

var CX=0,CY=0,CZ=0;
MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length;

var CR={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.84,C:0.76,N:0.71,O:0.66,F:0.57,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,K:2.03,Ca:1.76,Fe:1.32,Cu:1.32,Zn:1.22,Br:1.20,I:1.39};
function getR(el){return(CR[el]||1.50)*0.5}

function getOccupancyTexture(color,occupancy){
    var key=color+'_'+occupancy.toFixed(3);
    if(occTextureCache[key])return occTextureCache[key];
    var canvas=document.createElement('canvas');
    canvas.width=4;canvas.height=256;
    var ctx=canvas.getContext('2d');
    var boundY=Math.round((1-occupancy)*255);
    ctx.fillStyle='#FFFFFF';
    ctx.fillRect(0,0,4,boundY);
    ctx.fillStyle=color;
    ctx.fillRect(0,boundY,4,256-boundY);
    var tex=new THREE.CanvasTexture(canvas);
    tex.needsUpdate=true;
    occTextureCache[key]=tex;
    return tex;
}

function detectCrystalBaseBonds(){
    if(!CRY||!CRY_INV)return;
    var baseAtoms=CRY.baseAtoms;
    var lv=CRY.latticeVectors;
    var inv=CRY_INV;
    var n=baseAtoms.length;
    var CR2={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.85,C:0.76,N:0.71,O:0.66,F:0.57,Ne:0.58,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,Ar:1.06,K:2.03,Ca:1.76,Sc:1.70,Ti:1.60,V:1.53,Cr:1.39,Mn:1.39,Fe:1.32,Co:1.26,Ni:1.24,Cu:1.32,Zn:1.22,Ga:1.22,Ge:1.20,As:1.19,Se:1.20,Br:1.20,Kr:1.16,I:1.39};
    var BS={'C+C':[{o:3,l:1.20,t:0.05},{o:1.5,l:1.39,t:0.05},{o:2,l:1.38,t:0.05},{o:1,l:1.51,t:0.10}],'C+N':[{o:2,l:1.26,t:0.05},{o:1.5,l:1.36,t:0.05},{o:1,l:1.43,t:0.10},{o:3,l:1.16,t:0.06}],'C+O':[{o:2,l:1.24,t:0.05},{o:1,l:1.39,t:0.05}],'N+N':[{o:1,l:1.41,t:0.10},{o:2,l:1.25,t:0.06},{o:3,l:1.10,t:0.06}],'N+O':[{o:2,l:1.20,t:0.06},{o:1.5,l:1.30,t:0.06},{o:1,l:1.40,t:0.15}],'O+O':[{o:2,l:1.21,t:0.06},{o:1,l:1.48,t:0.15}],'C+S':[{o:1.5,l:1.73,t:0.06},{o:2,l:1.60,t:0.10},{o:1,l:1.82,t:0.15}],'C+F':[{o:1,l:1.33,t:0.10}],'C+H':[{o:1,l:0.97,t:0.15}],'N+H':[{o:1,l:0.88,t:0.15}],'O+H':[{o:1,l:0.85,t:0.15}]};
    var BC={HH:0,CH:1.3,HO:1.2,HN:1.3,CC:1.9,CO:1.7,CN:1.7,NN:1.7,NO:1.8,CF:1.6,CS:2.0};
    function pk(e1,e2){e1=e1.toUpperCase();e2=e2.toUpperCase();return e1<e2?e1+e2:e2+e1}
    function sk(e1,e2){e1=e1.charAt(0).toUpperCase()+e1.slice(1).toLowerCase();e2=e2.charAt(0).toUpperCase()+e2.slice(1).toLowerCase();return e1<e2?e1+'+'+e2:e2+'+'+e1}
    function gbo(el1,el2,d){
        var p=pk(el1,el2);var co=BC[p];
        if(co!==undefined){if(d>co)return 0}else{var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;if(d>(r1+r2)+0.5)return 0}
        var s=sk(el1,el2);var sp=BS[s];
        if(sp){for(var k=0;k<sp.length;k++){if(Math.abs(d-sp[k].l)<=sp[k].t)return sp[k].o}var bo=1,md=Infinity;for(var k=0;k<sp.length;k++){var df=Math.abs(d-sp[k].l);if(df<md){md=df;bo=sp[k].o}}return bo}
        var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;var rs=r1+r2;var ratio=rs?d/rs:1;
        if(ratio<0.85)return 3;if(ratio<0.90)return 2;return 1;
    }
    function cartToFrac(a){
        return[inv[0][0]*a.x+inv[0][1]*a.y+inv[0][2]*a.z,
               inv[1][0]*a.x+inv[1][1]*a.y+inv[1][2]*a.z,
               inv[2][0]*a.x+inv[2][1]*a.y+inv[2][2]*a.z];
    }
    var bonds=[];
    for(var i=0;i<n;i++){
        var f1=cartToFrac(baseAtoms[i]);
        for(var j=i+1;j<n;j++){
            var f2=cartToFrac(baseAtoms[j]);
            var fdx=f2[0]-f1[0],fdy=f2[1]-f1[1],fdz=f2[2]-f1[2];
            var sx=Math.round(fdx),sy=Math.round(fdy),sz=Math.round(fdz);
            var wfx=fdx-sx,wfy=fdy-sy,wfz=fdz-sz;
            var dx=wfx*lv[0][0]+wfy*lv[1][0]+wfz*lv[2][0];
            var dy=wfx*lv[0][1]+wfy*lv[1][1]+wfz*lv[2][1];
            var dz=wfx*lv[0][2]+wfy*lv[1][2]+wfz*lv[2][2];
            var dist=Math.sqrt(dx*dx+dy*dy+dz*dz);
            var bo=gbo(baseAtoms[i].element,baseAtoms[j].element,dist);
            if(bo>0){
                bonds.push({atom1:i,atom2:j,order:bo,shift:[-sx,-sy,-sz]});
            }
        }
    }
    CRY.baseBonds=bonds;
}

function rebuildCrystal(){
    if(!CRY)return;
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    var baseAtoms=CRY.baseAtoms;
    var baseBonds=CRY.baseBonds;

    if(!CRY_INV){MD.atoms=[];MD.bonds=[];updateMolInfo();return}
    var inv=CRY_INV;

    if(baseBonds.length===0||!baseBonds[0].shift){
        detectCrystalBaseBonds();
        baseBonds=CRY.baseBonds;
    }

    var iMin=Math.floor(aMin),iMax=Math.ceil(aMax);
    var jMin=Math.floor(bMin),jMax=Math.ceil(bMax);
    var kMin=Math.floor(cMin),kMax=Math.ceil(cMax);
    var EPS=1e-9;

    var atoms=[];
    var bonds=[];
    var idx=0;
    var cellMap={};

    baseAtoms.forEach(function(ba,bi){
        var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
        var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
        var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
        fx-=Math.floor(fx);fy-=Math.floor(fy);fz-=Math.floor(fz);
        for(var ci=iMin;ci<=iMax;ci++){
            for(var cj=jMin;cj<=jMax;cj++){
                for(var ck=kMin;ck<=kMax;ck++){
                    var ax=fx+ci,ay=fy+cj,az=fz+ck;
                    if(ax<aMin-EPS||ax>aMax+EPS)continue;
                    if(ay<bMin-EPS||ay>bMax+EPS)continue;
                    if(az<cMin-EPS||az>cMax+EPS)continue;
                    var ox=ax*lv[0][0]+ay*lv[1][0]+az*lv[2][0];
                    var oy=ax*lv[0][1]+ay*lv[1][1]+az*lv[2][1];
                    var oz=ax*lv[0][2]+ay*lv[1][2]+az*lv[2][2];
                    var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
                    atoms.push({
                        element:ba.element,
                        x:ox,y:oy,z:oz,
                        index:idx,
                        occupancy:ba.occupancy,
                        color:color,
                        baseIdx:ba.baseIdx!=null?ba.baseIdx:bi,
                        cellI:ci,cellJ:cj,cellK:ck
                    });
                    cellMap[bi+'_'+ci+'_'+cj+'_'+ck]=idx;
                    idx++;
                }
            }
        }
    });

    var seenPairs={};
    var seenSplits={};
    baseBonds.forEach(function(bb){
        var b1=bb.atom1,b2=bb.atom2;
        var s=bb.shift||[0,0,0];
        var hasShift=s[0]||s[1]||s[2];
        for(var ci1=iMin;ci1<=iMax;ci1++){
            for(var cj1=jMin;cj1<=jMax;cj1++){
                for(var ck1=kMin;ck1<=kMax;ck1++){
                    var idx1=cellMap[b1+'_'+ci1+'_'+cj1+'_'+ck1];
                    if(idx1==null)continue;
                    var ci2=ci1+s[0],cj2=cj1+s[1],ck2=ck1+s[2];
                    var idx2=cellMap[b2+'_'+ci2+'_'+cj2+'_'+ck2];
                    if(idx2!=null){
                        var key=idx1<idx2?idx1+'_'+idx2:idx2+'_'+idx1;
                        if(!seenPairs[key]){
                            seenPairs[key]=true;
                            bonds.push({atom1:idx1,atom2:idx2,order:bb.order});
                        }
                    }else{
                        var idx2s=cellMap[b2+'_'+ci1+'_'+cj1+'_'+ck1];
                        if(idx2s!=null){
                            var skey='F_'+idx1+'_'+idx2s+'_'+s[0]+'_'+s[1]+'_'+s[2];
                            if(!seenSplits[skey]){
                                seenSplits[skey]=true;
                                bonds.push({atom1:idx1,atom2:idx2s,order:bb.order,shift:s});
                            }
                        }
                    }
                }
            }
        }
        if(hasShift){
            var rs=[-s[0],-s[1],-s[2]];
            for(var ci2=iMin;ci2<=iMax;ci2++){
                for(var cj2=jMin;cj2<=jMax;cj2++){
                    for(var ck2=kMin;ck2<=kMax;ck2++){
                        var idx2r=cellMap[b2+'_'+ci2+'_'+cj2+'_'+ck2];
                        if(idx2r==null)continue;
                        var ri1=ci2+rs[0],rj1=cj2+rs[1],rk1=ck2+rs[2];
                        var ridx1=cellMap[b1+'_'+ri1+'_'+rj1+'_'+rk1];
                        if(ridx1==null){
                            var idx1r=cellMap[b1+'_'+ci2+'_'+cj2+'_'+ck2];
                            if(idx1r!=null){
                                var skey='R_'+idx2r+'_'+idx1r+'_'+rs[0]+'_'+rs[1]+'_'+rs[2];
                                if(!seenSplits[skey]){
                                    seenSplits[skey]=true;
                                    bonds.push({atom1:idx2r,atom2:idx1r,order:bb.order,shift:rs});
                                }
                            }
                        }
                    }
                }
            }
        }
    });

    MD.atoms=atoms;
    MD.bonds=bonds;
    updateMolInfo();
}

function removeDisorder(){
    if(!CRY)return;
    var newBaseAtoms=[];
    var oldToNew={};
    CRY.baseAtoms.forEach(function(ba,i){
        var occ=ba.occupancy!=null?ba.occupancy:1;
        if(occ<0.5)return;
        var newIdx=newBaseAtoms.length;
        oldToNew[i]=newIdx;
        var newOcc=occ>0.5?1:occ;
        var color=ba.color||((MD.atomColors&&MD.atomColors[ba.element])||'#FF1493');
        newBaseAtoms.push({
            element:ba.element,x:ba.x,y:ba.y,z:ba.z,
            index:newIdx,occupancy:newOcc,baseIdx:newIdx,color:color
        });
    });
    CRY.baseBonds=CRY.baseBonds.filter(function(b){
        return oldToNew[b.atom1]!=null&&oldToNew[b.atom2]!=null;
    }).map(function(b){
        return{atom1:oldToNew[b.atom1],atom2:oldToNew[b.atom2],order:b.order,shift:b.shift};
    });
    CRY.baseAtoms=newBaseAtoms;
    rebuildCrystal();
    rebuildScene();
}

var cellWireframe=null;
function buildCellWireframe(){
    if(cellWireframe){moleculeGroup.remove(cellWireframe);cellWireframe=null}
    if(!CRY)return;
    var lv=CRY.latticeVectors;
    var aMin=boundary.aMin,aMax=boundary.aMax;
    var bMin=boundary.bMin,bMax=boundary.bMax;
    var cMin=boundary.cMin,cMax=boundary.cMax;
    var corners=[];
    for(var i=0;i<2;i++)for(var j=0;j<2;j++)for(var k=0;k<2;k++){
        var ai=i?aMax:aMin,aj=j?bMax:bMin,ak=k?cMax:cMin;
        var x=ai*lv[0][0]+aj*lv[1][0]+ak*lv[2][0];
        var y=ai*lv[0][1]+aj*lv[1][1]+ak*lv[2][1];
        var z=ai*lv[0][2]+aj*lv[1][2]+ak*lv[2][2];
        corners.push(new THREE.Vector3(x-CX,y-CY,z-CZ));
    }
    var edges=[[0,4],[1,5],[2,6],[3,7],[0,2],[1,3],[4,6],[5,7],[0,1],[2,3],[4,5],[6,7]];
    var positions=[];
    edges.forEach(function(e){
        positions.push(corners[e[0]].x,corners[e[0]].y,corners[e[0]].z);
        positions.push(corners[e[1]].x,corners[e[1]].y,corners[e[1]].z);
    });
    var geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
    var mat=new THREE.LineBasicMaterial({color:0x88aaff});
    cellWireframe=new THREE.LineSegments(geo,mat);
    moleculeGroup.add(cellWireframe);
}

function updateAxesIndicator(){
    if(!CRY)return;
    var el=document.getElementById('axes-indicator');
    if(!el)return;
    var lv=CRY.latticeVectors;
    var vectors=[
        {dir:new THREE.Vector3(lv[0][0],lv[0][1],lv[0][2]).normalize(),color:'#FF4444',label:'a'},
        {dir:new THREE.Vector3(lv[1][0],lv[1][1],lv[1][2]).normalize(),color:'#44FF44',label:'b'},
        {dir:new THREE.Vector3(lv[2][0],lv[2][1],lv[2][2]).normalize(),color:'#4488FF',label:'c'}
    ];
    vectors.forEach(function(v){v.dir.applyQuaternion(rotQuat)});
    vectors.sort(function(a,b){return a.dir.z-b.dir.z});
    var cx=45,cy=45,scale=30;
    var svg='<svg viewBox="0 0 90 90">';
    vectors.forEach(function(v){
        var x2=cx+v.dir.x*scale;
        var y2=cy-v.dir.y*scale;
        var depth=(v.dir.z+1)/2;
        var opacity=0.4+depth*0.6;
        var width=1.5+depth*1.5;
        svg+='<line x1="'+cx+'" y1="'+cy+'" x2="'+x2.toFixed(1)+'" y2="'+y2.toFixed(1)+'" stroke="'+v.color+'" stroke-width="'+width.toFixed(1)+'" stroke-linecap="round" opacity="'+opacity.toFixed(2)+'"/>';
        svg+='<circle cx="'+x2.toFixed(1)+'" cy="'+y2.toFixed(1)+'" r="3" fill="'+v.color+'" opacity="'+opacity.toFixed(2)+'"/>';
        svg+='<text x="'+(x2+4).toFixed(1)+'" y="'+(y2+4).toFixed(1)+'" fill="'+v.color+'" font-size="10" font-weight="bold" opacity="'+opacity.toFixed(2)+'">'+v.label+'</text>';
    });
    svg+='<circle cx="'+cx+'" cy="'+cy+'" r="2" fill="#888"/>';
    svg+='</svg>';
    el.innerHTML=svg;
}

var atomMeshes=[];
var bondMeshes=[];

function rebuildScene(){
    while(moleculeGroup.children.length>0)moleculeGroup.remove(moleculeGroup.children[0]);
    atomMeshes.length=0;
    bondMeshes.length=0;
    cellWireframe=null;
    CX=0;CY=0;CZ=0;
    MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
    if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    MD.atoms.forEach(function(a,i){
        a.index=i;
        var r=getR(a.element);
        var g=new THREE.SphereGeometry(r,32,24);
        var occ=a.occupancy!=null?a.occupancy:1;
        var m;
        if(occ<0.999){
            var tex=getOccupancyTexture(a.color,occ);
            m=new THREE.MeshPhongMaterial({color:0xffffff,map:tex,shininess:80,specular:0x444444});
        }else{
            m=new THREE.MeshPhongMaterial({color:new THREE.Color(a.color),shininess:80,specular:0x444444});
        }
        var mesh=new THREE.Mesh(g,m);
        mesh.position.set(a.x-CX,a.y-CY,a.z-CZ);
        mesh.userData={element:a.element,index:i};
        moleculeGroup.add(mesh);
        atomMeshes.push(mesh);
    });
    MD.bonds.forEach(function(b){createBond(b)});
    if(CRY)buildCellWireframe();
    highlightSelected();
}
function updateScenePositions(){
    CX=0;CY=0;CZ=0;
    MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
    if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    atomMeshes.forEach(function(m,i){var a=MD.atoms[i];if(a)m.position.set(a.x-CX,a.y-CY,a.z-CZ)});
    for(var i=bondMeshes.length-1;i>=0;i--){moleculeGroup.remove(bondMeshes[i])}
    bondMeshes.length=0;
    MD.bonds.forEach(function(b){createBond(b)});
    if(CRY)buildCellWireframe();
    highlightSelected();
}

function getPerp(dir){
    var up=Math.abs(dir.y)<0.99?new THREE.Vector3(0,1,0):new THREE.Vector3(1,0,0);
    return new THREE.Vector3().crossVectors(dir,up).normalize();
}

function createBond(b){
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2)return;
    var a2x=a2.x,a2y=a2.y,a2z=a2.z;
    if(b.shift&&CRY){
        var lv=CRY.latticeVectors;
        a2x+=b.shift[0]*lv[0][0]+b.shift[1]*lv[1][0]+b.shift[2]*lv[2][0];
        a2y+=b.shift[0]*lv[0][1]+b.shift[1]*lv[1][1]+b.shift[2]*lv[2][1];
        a2z+=b.shift[0]*lv[0][2]+b.shift[1]*lv[1][2]+b.shift[2]*lv[2][2];
    }
    var s=new THREE.Vector3(a1.x-CX,a1.y-CY,a1.z-CZ);
    var e=new THREE.Vector3(a2x-CX,a2y-CY,a2z-CZ);
    var d=new THREE.Vector3().subVectors(e,s);
    var l=d.length();
    var mp=new THREE.Vector3().addVectors(s,e).multiplyScalar(0.5);
    var br=0.12,ord=b.order||1;
    var c1=new THREE.Color(a1.color),c2=new THREE.Color(a2.color);
    if(ord<1.25){hBond(s,mp,d,l/2,br,c1);hBond(mp,e,d,l/2,br,c2)}
    else if(ord<1.75){var off=0.10,p=getPerp(d).multiplyScalar(off);
        hBond(s,mp,d,l/2,br,c1);hBond(mp,e,d,l/2,br,c2);
        hDashedBond(s.clone().add(p),e.clone().add(p),d,l,br*0.7,c1,6);
    }else if(ord<2.5){var off=0.12,p=getPerp(d).multiplyScalar(off);
        hBond(s.clone().add(p),mp.clone().add(p),d,l/2,br*0.6,c1);hBond(mp.clone().add(p),e.clone().add(p),d,l/2,br*0.6,c2);
        hBond(s.clone().sub(p),mp.clone().sub(p),d,l/2,br*0.6,c1);hBond(mp.clone().sub(p),e.clone().sub(p),d,l/2,br*0.6,c2);
    }else if(ord<3.5){var off=0.15,p=getPerp(d).multiplyScalar(off);
        hBond(s,mp,d,l/2,br*0.45,c1);hBond(mp,e,d,l/2,br*0.45,c2);
        hBond(s.clone().add(p),mp.clone().add(p),d,l/2,br*0.45,c1);hBond(mp.clone().add(p),e.clone().add(p),d,l/2,br*0.45,c2);
        hBond(s.clone().sub(p),mp.clone().sub(p),d,l/2,br*0.45,c1);hBond(mp.clone().sub(p),e.clone().sub(p),d,l/2,br*0.45,c2);
    }else{hBond(s,mp,d,l/2,br,c1);hBond(mp,e,d,l/2,br,c2)}
}

function hDashedBond(s,e,d,hl,r,c,dashes){
    var seg=hl/dashes, gap=seg*0.35, dashLen=seg-gap;
    var dir=d.clone().normalize();
    for(var k=0;k<dashes;k++){
        var t0=k*seg+gap*0.5;
        var t1=t0+dashLen;
        if(t1>hl)t1=hl;
        var ds=s.clone().add(dir.clone().multiplyScalar(t0));
        var de=s.clone().add(dir.clone().multiplyScalar(t1));
        var dm=new THREE.Vector3().addVectors(ds,de).multiplyScalar(0.5);
        var dl=t1-t0;
        if(dl<0.001)continue;
        var g=new THREE.CylinderGeometry(r,r,dl,6,1);
        var m=new THREE.MeshPhongMaterial({color:c,shininess:40,specular:0x222222});
        var mesh=new THREE.Mesh(g,m);
        mesh.position.copy(dm);
        var axis=new THREE.Vector3(0,1,0);
        mesh.quaternion.copy(new THREE.Quaternion().setFromUnitVectors(axis,dir));
        moleculeGroup.add(mesh);
        bondMeshes.push(mesh);
    }
}
function hBond(s,e,d,hl,r,c){
    var g=new THREE.CylinderGeometry(r,r,hl,8,1);
    var m=new THREE.MeshPhongMaterial({color:c,shininess:40,specular:0x222222});
    var mesh=new THREE.Mesh(g,m);
    var mid=new THREE.Vector3().addVectors(s,e).multiplyScalar(0.5);
    mesh.position.copy(mid);
    var axis=new THREE.Vector3(0,1,0);
    mesh.quaternion.copy(new THREE.Quaternion().setFromUnitVectors(axis,d.clone().normalize()));
    moleculeGroup.add(mesh);
    bondMeshes.push(mesh);
}

var isRot=false,isPan=false,prevM={x:0,y:0},panX=0,panY=0,camDist=10;
var rotQuat=new THREE.Quaternion();
var currentMode='view';
var selectedAtoms=[];
var originalCoords=null;
var modalCallback=null;

if(CRY){rebuildCrystal();}
rebuildScene();

if(CRY){
    var crystalPanel=document.getElementById('crystal-panel');
    if(crystalPanel)crystalPanel.style.display='block';
    var axesEl=document.getElementById('axes-indicator');
    if(axesEl)axesEl.style.display='block';
    updateAxesIndicator();
    ['bnd-a-min','bnd-a-max','bnd-b-min','bnd-b-max','bnd-c-min','bnd-c-max'].forEach(function(id){
        var inp=document.getElementById(id);
        if(inp){
            inp.addEventListener('change',function(){
                var v=parseFloat(inp.value);
                if(isNaN(v))v=0;
                v=Math.max(-5,Math.min(5,v));
                inp.value=v;
                if(id==='bnd-a-min')boundary.aMin=v;
                else if(id==='bnd-a-max')boundary.aMax=v;
                else if(id==='bnd-b-min')boundary.bMin=v;
                else if(id==='bnd-b-max')boundary.bMax=v;
                else if(id==='bnd-c-min')boundary.cMin=v;
                else if(id==='bnd-c-max')boundary.cMax=v;
                if(boundary.aMax<boundary.aMin)boundary.aMax=boundary.aMin;
                if(boundary.bMax<boundary.bMin)boundary.bMax=boundary.bMin;
                if(boundary.cMax<boundary.cMin)boundary.cMax=boundary.cMin;
                rebuildCrystal();
                rebuildScene();
            });
        }
    });
    var rmDisorderBtn=document.getElementById('bnd-remove-disorder');
    if(rmDisorderBtn){
        rmDisorderBtn.addEventListener('click',function(){
            if(!CRY)return;
            var removedCount=0;
            CRY.baseAtoms.forEach(function(ba){
                var occ=ba.occupancy!=null?ba.occupancy:1;
                if(occ<0.5)removedCount++;
            });
            if(removedCount===0){
                showModal('<h3>No Disorder</h3><div class="current-val">No atoms with occupancy &lt; 0.5 found.</div><div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">OK</button></div>',null);
                document.getElementById('m-cancel').addEventListener('click',hideModal);
                return;
            }
            showModal('<h3>Remove Disorder</h3>'+
                '<div class="current-val">Remove '+removedCount+' atom(s) with occupancy &lt; 0.5?<br>Atoms with occupancy &gt; 0.5 will be set to 1.0.<br>Atoms with occupancy = 0.5 will be kept unchanged.</div>'+
                '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok mbtn-danger" id="m-ok">Remove</button></div>',null);
            document.getElementById('m-ok').addEventListener('click',function(){hideModal();removeDisorder()});
            document.getElementById('m-cancel').addEventListener('click',hideModal);
        });
    }
}

var maxD=0;
MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
var initCam=maxD*2.5+5;
camera.position.set(0,0,initCam);camera.lookAt(0,0,0);
camDist=initCam;

var MODE_INFO={view:'View Mode',bondLength:'Bond Length - Click 2 atoms',bondAngle:'Bond Angle - Click 3 atoms (central 2nd)',dihedral:'Dihedral - Click 4 atoms',addAtom:'Add Atom - Click anchor atom',deleteAtom:'Delete Atom - Click atom to delete',selectAtoms:'Select Atoms - Input indices or element symbols'};

function setMode(m){
    if(currentMode===m){
        if(m==='selectAtoms')showSelectAtomsModal();
        return;
    }
    currentMode=m;selectedAtoms=[];diffSelectedAtoms=[];originalCoords=null;hideModal();highlightSelected();
    modeInfoEl.textContent=MODE_INFO[m]||m;
    selInfoEl.textContent='';
    document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.classList.toggle('active',b.dataset.mode===m)});
    if(m==='selectAtoms')showSelectAtomsModal();
}

function resetSelection(){
    selectedAtoms=[];diffSelectedAtoms=[];originalCoords=null;highlightSelected();
    selInfoEl.textContent='';
}

document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.addEventListener('click',function(){setMode(this.dataset.mode)})});
document.getElementById('reset-btn').addEventListener('click',function(){
    rotQuat.identity();panX=0;panY=0;
    if(diffMode){
        camDist=camDiffInit();diffCamDist=camDist;
        diffRotQuat.identity();diffPanX=0;diffPanY=0;
    }else{
        camDist=initCam;camera.position.set(0,0,camDist);
    }
    updateTransform();
});
document.getElementById('save-btn').addEventListener('click',doSave);
document.getElementById('diff-btn').addEventListener('click',function(){
    if(diffMode){
        exitDiff();
    }else{
        modeInfoEl.textContent='Diff: selecting file...';
        vscodeApi.postMessage({command:'diffFile'});
    }
});

var diffMode=false;
var diffData=null;
var diffPivot=null;
var diffMolGroup=null;
var diffAtomMeshes=[];
var diffBondMeshes=[];
var diffCX=0,diffCY=0,diffCZ=0;
var diffRotQuat=new THREE.Quaternion();
var diffThresholdPct=2;
var diffDetBonds1=null;
var diffDetBonds2=null;
var diffPanX=0,diffPanY=0;
var diffCamDist=10;
var diffTransformSide='left';
var diffActiveSide='left';
var diffSelectedAtoms=[];
var diffMapping=null;
var diffReverseMapping=null;
var diffPanelEl=document.getElementById('diff-panel');
var diffReopenEl=document.getElementById('diff-reopen');
var diffPanelHTML='';
var diffResultsHTML='';
var diffLabelLeft=document.getElementById('diff-label');
var diffLabelRight=document.getElementById('diff-label-right');

window.addEventListener('message',function(event){
    var msg=event.data;
    if(msg.command==='diffResult'){
        if(msg.cancelled){
            modeInfoEl.textContent='View Mode';
            return;
        }
        startDiff(msg);
    }
});

function buildAdjacency(atoms,bonds){
    var adj=[];
    for(var i=0;i<atoms.length;i++)adj[i]=[];
    bonds.forEach(function(b){
        if(adj[b.atom1])adj[b.atom1].push({to:b.atom2,order:b.order});
        if(adj[b.atom2])adj[b.atom2].push({to:b.atom1,order:b.order});
    });
    return adj;
}

function atomSignature(i,atoms,adj){
    var el=atoms[i].element;
    var deg=adj[i].length;
    var nbEls=adj[i].map(function(e){return atoms[e.to].element}).sort().join(',');
    return el+'|'+deg+'|'+nbEls;
}

function findAtomMapping(atoms1,bonds1,atoms2,bonds2){
    var n1=atoms1.length,n2=atoms2.length;
    if(n1!==n2)return null;

    var adj1=buildAdjacency(atoms1,bonds1);
    var adj2=buildAdjacency(atoms2,bonds2);

    var elCount1={},elCount2={};
    atoms1.forEach(function(a){elCount1[a.element]=(elCount1[a.element]||0)+1});
    atoms2.forEach(function(a){elCount2[a.element]=(elCount2[a.element]||0)+1});
    var els1=Object.keys(elCount1).sort(),els2=Object.keys(elCount2).sort();
    if(els1.length!==els2.length)return null;
    for(var k=0;k<els1.length;k++){
        if(els1[k]!==els2[k]||elCount1[els1[k]]!==elCount2[els2[k]])return null;
    }

    var deg1={},deg2={};
    for(var i=0;i<n1;i++){var d=adj1[i].length;deg1[d]=(deg1[d]||0)+1}
    for(var i=0;i<n2;i++){var d=adj2[i].length;deg2[d]=(deg2[d]||0)+1}
    var dk1=Object.keys(deg1).sort(),dk2=Object.keys(deg2).sort();
    if(dk1.length!==dk2.length)return null;
    for(var k=0;k<dk1.length;k++){if(dk1[k]!==dk2[k]||deg1[dk1[k]]!==deg2[dk2[k]])return null}

    // Precompute signatures and candidate lists
    var sig1=[],sig2=[];
    for(var i=0;i<n1;i++)sig1[i]=atomSignature(i,atoms1,adj1);
    for(var j=0;j<n2;j++)sig2[j]=atomSignature(j,atoms2,adj2);
    var candidates=[];
    for(var i=0;i<n1;i++){
        candidates[i]=[];
        for(var j=0;j<n2;j++){
            if(sig1[i]===sig2[j])candidates[i].push(j);
        }
        if(candidates[i].length===0)return null;
    }

    var map=new Array(n1).fill(-1);
    var rmap=new Array(n2).fill(-1);

    function consistent(a1,a2){
        if(atoms1[a1].element!==atoms2[a2].element)return false;
        var nb2set={};
        adj2[a2].forEach(function(e){nb2set[e.to]=true});
        for(var k=0;k<adj1[a1].length;k++){
            var nb1=adj1[a1][k].to;
            if(map[nb1]!==-1){
                if(nb2set[map[nb1]]===undefined)return false;
            }
        }
        var nb1set={};
        adj1[a1].forEach(function(e){nb1set[e.to]=true});
        for(var k=0;k<adj2[a2].length;k++){
            var nb2=adj2[a2][k].to;
            if(rmap[nb2]!==-1){
                if(nb1set[rmap[nb2]]===undefined)return false;
            }
        }
        return true;
    }

    function nextUnassigned(){
        var best=-1,bestKey=null;
        for(var i=0;i<n1;i++){
            if(map[i]!==-1)continue;
            var mappedNbrs=0;
            for(var k=0;k<adj1[i].length;k++){
                if(map[adj1[i][k].to]!==-1)mappedNbrs++;
            }
            var avail=0;
            for(var c=0;c<candidates[i].length;c++){
                if(rmap[candidates[i][c]]===-1)avail++;
            }
            var key=(mappedNbrs>0?1000000:0)-avail*1000+mappedNbrs;
            if(bestKey===null||key>bestKey){bestKey=key;best=i}
        }
        return best;
    }

    var iterCount=0;
    var maxIter=500000;

    function backtrack(depth){
        if(iterCount++>maxIter)return false;
        if(depth===n1)return true;
        var a1=nextUnassigned();
        if(a1===-1)return false;
        for(var c=0;c<candidates[a1].length;c++){
            var a2=candidates[a1][c];
            if(rmap[a2]!==-1)continue;
            if(!consistent(a1,a2))continue;
            map[a1]=a2;rmap[a2]=a1;
            if(backtrack(depth+1))return true;
            map[a1]=-1;rmap[a2]=-1;
        }
        return false;
    }

    if(backtrack(0))return map;
    return null;
}

function dist3d(a,b){var dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;return Math.sqrt(dx*dx+dy*dy+dz*dz)}
function angle3d(a,b,c){
    var v1={x:a.x-b.x,y:a.y-b.y,z:a.z-b.z},v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var d1=Math.sqrt(v1.x*v1.x+v1.y*v1.y+v1.z*v1.z),d2=Math.sqrt(v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(d1*d2);
    return Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI;
}
function dihedral3d(a,b,c,d){
    var v1={x:b.x-a.x,y:b.y-a.y,z:b.z-a.z};
    var v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var v3={x:d.x-c.x,y:d.y-c.y,z:d.z-c.z};
    var n1={x:v1.y*v2.z-v1.z*v2.y,y:v1.z*v2.x-v1.x*v2.z,z:v1.x*v2.y-v1.y*v2.x};
    var n2={x:v2.y*v3.z-v2.z*v3.y,y:v2.z*v3.x-v2.x*v3.z,z:v2.x*v3.y-v2.y*v3.x};
    var d1=Math.sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z),d2=Math.sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(n1.x*n2.x+n1.y*n2.y+n1.z*n2.z)/(d1*d2);
    var sign=v1.x*n2.x+v1.y*n2.y+v1.z*n2.z;
    var ang=Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI;
    return sign<0?-ang:ang;
}

function computeConformationDiffs(atoms1,bonds1,atoms2,bonds2,mapping,threshPct){
    var diffs={bonds:[],angles:[],dihedrals:[]};
    var tp=threshPct||2;

    var bondSet2={};
    bonds2.forEach(function(b){bondSet2[Math.min(b.atom1,b.atom2)+'-'+Math.max(b.atom1,b.atom2)]=b.order});
    bonds1.forEach(function(b){
        var a1m=Math.min(mapping[b.atom1],mapping[b.atom2]);
        var a2m=Math.max(mapping[b.atom1],mapping[b.atom2]);
        var key=a1m+'-'+a2m;
        if(bondSet2[key]!==undefined){
            var d1=dist3d(atoms1[b.atom1],atoms1[b.atom2]);
            var d2=dist3d(atoms2[mapping[b.atom1]],atoms2[mapping[b.atom2]]);
            var dd=Math.abs(d1-d2);
            var avg=(d1+d2)/2;
            var pct=avg>0.01?dd/avg*100:0;
            if(pct>tp){
                diffs.bonds.push({i1:b.atom1,i2:b.atom2,d1:d1,d2:d2,diff:dd,pct:pct});
            }
        }
    });

    var adj1=buildAdjacency(atoms1,bonds1);
    var angleSeen={};
    for(var j=0;j<atoms1.length;j++){
        var nbs=adj1[j].map(function(e){return e.to});
        for(var k=0;k<nbs.length;k++){
            for(var m=k+1;m<nbs.length;m++){
                var i1=nbs[k],i3=nbs[m];
                var key=Math.min(i1,j)+'-'+j+'-'+Math.max(i3,j);
                if(angleSeen[key])continue;
                angleSeen[key]=true;
                var a1=angle3d(atoms1[i1],atoms1[j],atoms1[i3]);
                var a2=angle3d(atoms2[mapping[i1]],atoms2[mapping[j]],atoms2[mapping[i3]]);
                var da=Math.abs(a1-a2);
                var avgA=(a1+a2)/2;
                var pctA=avgA>1.0?da/avgA*100:da*100;
                if(pctA>tp){
                    diffs.angles.push({i1:i1,i2:j,i3:i3,a1:a1,a2:a2,diff:da,pct:pctA});
                }
            }
        }
    }

    var dihedralSeen={};
    bonds1.forEach(function(b){
        var a1=b.atom1,a2=b.atom2;
        [[a1,a2],[a2,a1]].forEach(function(pair){
            var start=pair[0],mid=pair[1];
            adj1[mid].forEach(function(e2){
                var end=e2.to;
                if(end===start)return;
                adj1[start].forEach(function(e1){
                    var s0=e1.to;
                    if(s0===mid)return;
                    var key=[s0,start,mid,end].sort().join('-');
                    if(dihedralSeen[key])return;
                    dihedralSeen[key]=true;
                    var d1=dihedral3d(atoms1[s0],atoms1[start],atoms1[mid],atoms1[end]);
                    var d2=dihedral3d(atoms2[mapping[s0]],atoms2[mapping[start]],atoms2[mapping[mid]],atoms2[mapping[end]]);
                    var dd=Math.abs(d1-d2);
                    if(dd>180)dd=360-dd;
                    var avgD=(Math.abs(d1)+Math.abs(d2))/2;
                    var pctD=avgD>1.0?dd/avgD*100:dd*10;
                    if(pctD>tp){
                        diffs.dihedrals.push({i1:s0,i2:start,i3:mid,i4:end,d1:d1,d2:d2,diff:dd,pct:pctD});
                    }
                });
            });
        });
    });

    return diffs;
}

function startDiff(msg){
    diffData={atoms:msg.atoms,bonds:msg.bonds,title:msg.title||msg.fileName,fileName:msg.fileName};

    var detBonds1=detectBondsFromAtoms(MD.atoms);
    var detBonds2=detectBondsFromAtoms(diffData.atoms);
    diffDetBonds1=detBonds1;
    diffDetBonds2=detBonds2;

    var mapping=findAtomMapping(MD.atoms,detBonds1,diffData.atoms,detBonds2);

    if(!mapping){
        diffMapping=null;diffReverseMapping=null;
        showDiffPanel('<span class="diff-close">×</span>'+
            '<h4>Skeletons Differ</h4>'+
            '<div class="diff-row">Molecular skeletons are different — cannot compare.</div>'+
            '<div class="diff-row">Left: '+MD.atoms.length+' atoms, '+detBonds1.length+' bonds (detected)</div>'+
            '<div class="diff-row">Right: '+diffData.atoms.length+' atoms, '+detBonds2.length+' bonds (detected)</div>'+
            '<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Click Diff again to exit.</div>',true);
        diffPanelEl.classList.add('show');
        diffMode=true;
        diffLabelLeft.textContent='Original: '+(MD.title||'molecule');
        diffLabelLeft.classList.add('show');
        diffLabelRight.textContent='Diff: '+diffData.fileName+' (skeleton differs)';
        diffLabelRight.classList.add('show');
        modeInfoEl.textContent='Diff Mode (skeletons differ)';
        enterDiffRender();
        return;
    }

    var diffs=computeConformationDiffs(MD.atoms,detBonds1,diffData.atoms,detBonds2,mapping,diffThresholdPct);
    var totalDiff=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;

    diffMapping=mapping;
    diffReverseMapping=new Array(diffData.atoms.length).fill(-1);
    for(var mi=0;mi<mapping.length;mi++){diffReverseMapping[mapping[mi]]=mi}

    var html=buildDiffResultsHTML(diffs,totalDiff);
    showDiffPanel(html,true);
    wireDiffSlider();
    diffMode=true;
    diffLabelLeft.textContent='Original: '+(MD.title||'molecule');
    diffLabelLeft.classList.add('show');
    diffLabelRight.textContent='Diff: '+diffData.fileName;
    diffLabelRight.classList.add('show');
    modeInfoEl.textContent='Diff Mode ('+totalDiff+' differences)';

    enterDiffRender(diffs,mapping);
}

function buildDiffResultsHTML(diffs,totalDiff){
    var html='<span class="diff-close">×</span>';
    html+='<h4>Diff Results</h4>';
    html+='<div class="diff-row" style="margin-bottom:4px">';
    html+='<label style="font-size:11px">Threshold: <b><span id="thresh-val">'+diffThresholdPct.toFixed(1)+'</span>%</b></label>';
    html+='<input type="range" id="thresh-slider" min="0" max="20" step="0.5" value="'+diffThresholdPct+'" style="width:100%;accent-color:var(--vscode-textLink-foreground,#3794ff)">';
    html+='<div style="font-size:10px;color:var(--vscode-descriptionForeground,#999)">Relative difference: Δ/avg × 100%</div>';
    html+='</div>';
    html+='<div id="diff-results-list">';
    html+=buildDiffListHTML(diffs,totalDiff);
    html+='</div>';
    return html;
}

function buildDiffListHTML(diffs,totalDiff){
    var html='';
    if(totalDiff===0){
        html+='<div class="diff-row" style="color:#4ec9b0">Structures are identical (same conformation within thresholds).</div>';
    }else{
        html+='<div class="diff-row">Mapping: '+MD.atoms.length+' atoms matched. '+totalDiff+' differences found.</div>';
        if(diffs.bonds.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Length Differences ('+diffs.bonds.length+')</b></div>';
            diffs.bonds.sort(function(a,b){return b.pct-a.pct});
            diffs.bonds.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+
                    ': '+d.d1.toFixed(3)+' vs '+d.d2.toFixed(3)+' Å (Δ='+d.diff.toFixed(3)+', '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.angles.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Bond Angle Differences ('+diffs.angles.length+')</b></div>';
            diffs.angles.sort(function(a,b){return b.pct-a.pct});
            diffs.angles.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+'-'+MD.atoms[d.i3].element+(d.i3+1)+
                    ': '+d.a1.toFixed(1)+'° vs '+d.a2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+'°, '+d.pct.toFixed(2)+'%)</div>';
            });
        }
        if(diffs.dihedrals.length>0){
            html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Dihedral Differences ('+diffs.dihedrals.length+')</b></div>';
            diffs.dihedrals.sort(function(a,b){return b.pct-a.pct});
            diffs.dihedrals.forEach(function(d){
                html+='<div class="diff-row">  '+MD.atoms[d.i1].element+(d.i1+1)+'-'+MD.atoms[d.i2].element+(d.i2+1)+'-'+MD.atoms[d.i3].element+(d.i3+1)+'-'+MD.atoms[d.i4].element+(d.i4+1)+
                    ': '+d.d1.toFixed(1)+'° vs '+d.d2.toFixed(1)+'° (Δ='+d.diff.toFixed(1)+'°, '+d.pct.toFixed(2)+'%)</div>';
            });
        }
    }
    return html;
}

function wireDiffSlider(){
    var slider=diffPanelEl.querySelector('#thresh-slider');
    if(slider){
        slider.addEventListener('input',function(){
            diffThresholdPct=parseFloat(slider.value);
            var valEl=diffPanelEl.querySelector('#thresh-val');
            if(valEl)valEl.textContent=diffThresholdPct.toFixed(1);
            recomputeDiff();
        });
    }
}

function recomputeDiff(){
    if(!diffMapping||!diffDetBonds1||!diffDetBonds2)return;
    var diffs=computeConformationDiffs(MD.atoms,diffDetBonds1,diffData.atoms,diffDetBonds2,diffMapping,diffThresholdPct);
    var totalDiff=diffs.bonds.length+diffs.angles.length+diffs.dihedrals.length;
    var listEl=diffPanelEl.querySelector('#diff-results-list');
    if(listEl){
        listEl.innerHTML=buildDiffListHTML(diffs,totalDiff);
    }
    var fullHTML=buildDiffResultsHTML(diffs,totalDiff);
    diffPanelHTML=fullHTML;
    diffResultsHTML=fullHTML;
    modeInfoEl.textContent='Diff Mode ('+totalDiff+' differences)';
    updateDiffHighlights(diffs);
}

function updateDiffHighlights(diffs){
    var highlightAtoms={};
    var highlightBonds={};
    diffs.bonds.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;
        highlightBonds[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true;
    });
    diffs.angles.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;
    });
    diffs.dihedrals.forEach(function(d){
        highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;highlightAtoms[d.i4]=true;
    });

    atomMeshes.forEach(function(m,i){
        m.userData.diffHi=highlightAtoms[i]?true:false;
    });

    diffAtomMeshes.forEach(function(m,i){
        var origI=diffReverseMapping?diffReverseMapping[i]:-1;
        var isHi=origI>=0&&highlightAtoms[origI];
        m.userData.diffHi=isHi?true:false;
        if(m.userData.origColor){
            m.material.color.copy(m.userData.origColor);
        }
    });

    diffBondMeshes.forEach(function(mesh){
        if(mesh.userData&&mesh.userData.bondKey){
            var isHi=highlightBonds[mesh.userData.bondKey];
            if(isHi){
                mesh.material.color.set(0xff6600);
                mesh.material.emissive=new THREE.Color(0xff6600);
                mesh.material.emissiveIntensity=0.3;
            }else{
                mesh.material.color.copy(mesh.userData.origColor);
                mesh.material.emissive=new THREE.Color(0x000000);
                mesh.material.emissiveIntensity=0;
            }
        }
    });

    highlightSelected();
}

function showDiffPanel(html,isResults){
    diffPanelHTML=html;
    if(isResults)diffResultsHTML=html;
    diffPanelEl.innerHTML=html;
    diffPanelEl.classList.add('show');
    diffReopenEl.classList.remove('show');
    var closeBtn=diffPanelEl.querySelector('.diff-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            diffPanelEl.classList.remove('show');
            diffReopenEl.classList.add('show');
        });
    }
    wireDiffSlider();
}

diffReopenEl.addEventListener('click',function(){
    var html=diffResultsHTML||diffPanelHTML;
    diffPanelEl.innerHTML=html;
    diffPanelEl.classList.add('show');
    diffReopenEl.classList.remove('show');
    var closeBtn=diffPanelEl.querySelector('.diff-close');
    if(closeBtn){
        closeBtn.addEventListener('click',function(){
            diffPanelEl.classList.remove('show');
            diffReopenEl.classList.add('show');
        });
    }
    wireDiffSlider();
});

function enterDiffRender(diffs,mapping){
    diffPivot=new THREE.Group();scene.add(diffPivot);
    diffMolGroup=new THREE.Group();diffPivot.add(diffMolGroup);

    diffCX=0;diffCY=0;diffCZ=0;
    diffData.atoms.forEach(function(a){diffCX+=a.x;diffCY+=a.y;diffCZ+=a.z});
    if(diffData.atoms.length>0){diffCX/=diffData.atoms.length;diffCY/=diffData.atoms.length;diffCZ/=diffData.atoms.length}

    var highlightAtoms={};
    var highlightBonds={};
    if(diffs&&mapping){
        diffs.bonds.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;
            highlightBonds[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true;
        });
        diffs.angles.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;
        });
        diffs.dihedrals.forEach(function(d){
            highlightAtoms[d.i1]=true;highlightAtoms[d.i2]=true;highlightAtoms[d.i3]=true;highlightAtoms[d.i4]=true;
        });
    }

    diffData.atoms.forEach(function(a,i){
        var r=getR(a.element);
        var g=new THREE.SphereGeometry(r,32,24);
        var isHi=highlightAtoms[mapping?mapping.indexOf(i):-1];
        var origColor=new THREE.Color(a.color);
        var col=isHi?new THREE.Color(0xff6600):origColor;
        var m=new THREE.MeshPhongMaterial({color:col,shininess:80,specular:0x444444});
        if(isHi){m.emissive=new THREE.Color(0xff6600);m.emissiveIntensity=0.4}
        var mesh=new THREE.Mesh(g,m);
        mesh.position.set(a.x-diffCX,a.y-diffCY,a.z-diffCZ);
        mesh.userData={element:a.element,index:i,diffHi:isHi?true:false,origColor:origColor};
        diffMolGroup.add(mesh);
        diffAtomMeshes.push(mesh);
    });

    diffData.bonds.forEach(function(b){
        var a1=diffData.atoms[b.atom1],a2=diffData.atoms[b.atom2];
        if(!a1||!a2)return;
        var s=new THREE.Vector3(a1.x-diffCX,a1.y-diffCY,a1.z-diffCZ);
        var e=new THREE.Vector3(a2.x-diffCX,a2.y-diffCY,a2.z-diffCZ);
        var d=new THREE.Vector3().subVectors(e,s);
        var l=d.length();
        var mp=new THREE.Vector3().addVectors(s,e).multiplyScalar(0.5);
        var br=0.12;
        var origI1=mapping?mapping.indexOf(b.atom1):-1;
        var origI2=mapping?mapping.indexOf(b.atom2):-1;
        var bondKey=Math.min(origI1,origI2)+'-'+Math.max(origI1,origI2);
        var isHi=highlightBonds[bondKey];
        var c1Orig=new THREE.Color(a1.color);
        var c2Orig=new THREE.Color(a2.color);
        var c1=isHi?new THREE.Color(0xff6600):c1Orig;
        var c2=isHi?new THREE.Color(0xff6600):c2Orig;
        hBondDiff(s,mp,d,l/2,br,c1);
        hBondDiff(mp,e,d,l/2,br,c2);
        if(diffBondMeshes.length>=2){
            diffBondMeshes[diffBondMeshes.length-2].userData={origColor:c1Orig,bondKey:bondKey};
            diffBondMeshes[diffBondMeshes.length-1].userData={origColor:c2Orig,bondKey:bondKey};
        }
    });

    if(diffs&&mapping){
        var leftHi={};
        diffs.bonds.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true});
        diffs.angles.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true;leftHi[d.i3]=true});
        diffs.dihedrals.forEach(function(d){leftHi[d.i1]=true;leftHi[d.i2]=true;leftHi[d.i3]=true;leftHi[d.i4]=true});
        atomMeshes.forEach(function(m,i){
            m.userData.diffHi=leftHi[i]?true:false;
            if(leftHi[i]){
                m.material.emissive=new THREE.Color(0xff6600);
                m.material.emissiveIntensity=0.4;
            }
        });
        var leftBondHi={};
        diffs.bonds.forEach(function(d){leftBondHi[Math.min(d.i1,d.i2)+'-'+Math.max(d.i1,d.i2)]=true});
        bondMeshes.forEach(function(mesh){
            mesh.material.emissive=new THREE.Color(0x000000);
        });
    }

    var allAtoms=MD.atoms.concat(diffData.atoms);
    var maxD2=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD2)maxD2=dd});
    diffData.atoms.forEach(function(a){var dx=a.x-diffCX,dy=a.y-diffCY,dz=a.z-diffCZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD2)maxD2=dd});
    var newCam=maxD2*2.5+5;
    camDist=newCam;diffCamDist=newCam;
    camera.position.set(0,0,camDist);
    panX=0;panY=0;
    diffPanX=0;diffPanY=0;
    diffRotQuat.copy(rotQuat);
    diffTransformSide='left';
    updateTransform();
}

function hBondDiff(s,e,d,hl,r,c){
    var g=new THREE.CylinderGeometry(r,r,hl,8,1);
    var m=new THREE.MeshPhongMaterial({color:c,shininess:40,specular:0x222222});
    var mesh=new THREE.Mesh(g,m);
    var mid=new THREE.Vector3().addVectors(s,e).multiplyScalar(0.5);
    mesh.position.copy(mid);
    var axis=new THREE.Vector3(0,1,0);
    mesh.quaternion.copy(new THREE.Quaternion().setFromUnitVectors(axis,d.clone().normalize()));
    diffMolGroup.add(mesh);
    diffBondMeshes.push(mesh);
}

function exitDiff(){
    diffMode=false;
    if(diffPivot){scene.remove(diffPivot);diffPivot=null;diffMolGroup=null}
    diffAtomMeshes=[];diffBondMeshes=[];
    diffSelectedAtoms=[];
    diffMapping=null;diffReverseMapping=null;
    diffDetBonds1=null;diffDetBonds2=null;
    diffActiveSide='left';
    diffTransformSide='left';
    diffRotQuat.identity();
    diffPanX=0;diffPanY=0;diffCamDist=10;
    diffPanelEl.classList.remove('show');
    diffPanelEl.innerHTML='';
    diffReopenEl.classList.remove('show');
    diffPanelHTML='';
    diffResultsHTML='';
    diffLabelLeft.classList.remove('show');
    diffLabelRight.classList.remove('show');
    modeInfoEl.textContent='View Mode';
    atomMeshes.forEach(function(m){
        m.userData.diffHi=false;
        m.material.emissive=new THREE.Color(0x000000);
        m.material.emissiveIntensity=0;
    });
    highlightSelected();
    pivotGroup.visible=true;
    var w=container.clientWidth||window.innerWidth;
    var h=container.clientHeight||(window.innerHeight-60);
    if(w<1)w=window.innerWidth;
    if(h<1)h=window.innerHeight-60;
    renderer.setScissorTest(false);
    renderer.setViewport(0,0,w,h);
    renderer.setScissor(0,0,w,h);
    renderer.autoClear=true;
    camera.aspect=w/h;
    camera.updateProjectionMatrix();
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    camDist=maxD*2.5+5;
    camera.position.set(0,0,camDist);
    updateTransform();
}

function camDiffInit(){
    var maxD=0;
    MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    if(diffData){
        diffData.atoms.forEach(function(a){var dx=a.x-diffCX,dy=a.y-diffCY,dz=a.z-diffCZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
    }
    return maxD*2.5+5;
}

var currentFrame=0;
var totalFrames=MD.frames?MD.frames.length:0;
var autoPlayTimer=null;
var isAutoPlaying=false;
var frameNavEl=document.getElementById('frame-nav');
var frameSepEl=document.getElementById('frame-sep');
if(totalFrames>1){
    frameNavEl.classList.add('show');
    frameSepEl.classList.remove('hidden');
    updateFrameInfo();
}else{
    frameNavEl.classList.remove('show');
    frameSepEl.classList.add('hidden');
}
function updateFrameInfo(){
    var el=document.getElementById('frame-info');
    var numEl=document.getElementById('frame-num');
    if(el)el.textContent='/'+totalFrames+(MD.frames[currentFrame]?' - '+MD.frames[currentFrame].stepLabel:'');
    if(numEl)numEl.value=currentFrame+1;
}
function switchFrame(idx){
    if(idx<0||idx>=totalFrames)return;
    currentFrame=idx;
    var f=MD.frames[idx];
    MD.atoms=f.atoms.map(function(a,i){a.index=i;return a});
    MD.bonds=f.bonds||[];
    if(MD.bonds.length===0){
        MD.bonds=detectBondsFromAtoms(MD.atoms);
    }
    rebuildScene();
    updateFrameInfo();
}
function detectBondsFromAtoms(atoms){
    var CR2={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.85,C:0.76,N:0.71,O:0.66,F:0.57,Ne:0.58,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,Ar:1.06,K:2.03,Ca:1.76,Sc:1.70,Ti:1.60,V:1.53,Cr:1.39,Mn:1.39,Fe:1.32,Co:1.26,Ni:1.24,Cu:1.32,Zn:1.22,Ga:1.22,Ge:1.20,As:1.19,Se:1.20,Br:1.20,Kr:1.16,I:1.39};
    var BS={'C+C':[{o:3,l:1.20,t:0.05},{o:1.5,l:1.39,t:0.05},{o:2,l:1.38,t:0.05},{o:1,l:1.51,t:0.10}],'C+N':[{o:2,l:1.26,t:0.05},{o:1.5,l:1.36,t:0.05},{o:1,l:1.43,t:0.10},{o:3,l:1.16,t:0.06}],'C+O':[{o:2,l:1.24,t:0.05},{o:1,l:1.39,t:0.05}],'N+N':[{o:1,l:1.41,t:0.10},{o:2,l:1.25,t:0.06},{o:3,l:1.10,t:0.06}],'N+O':[{o:2,l:1.20,t:0.06},{o:1.5,l:1.30,t:0.06},{o:1,l:1.40,t:0.15}],'O+O':[{o:2,l:1.21,t:0.06},{o:1,l:1.48,t:0.15}],'C+S':[{o:1.5,l:1.73,t:0.06},{o:2,l:1.60,t:0.10},{o:1,l:1.82,t:0.15}],'C+F':[{o:1,l:1.33,t:0.10}],'C+H':[{o:1,l:0.97,t:0.15}],'N+H':[{o:1,l:0.88,t:0.15}],'O+H':[{o:1,l:0.85,t:0.15}]};
    var BC={HH:0,CH:1.3,HO:1.2,HN:1.3,CC:1.9,CO:1.7,CN:1.7,NN:1.7,NO:1.8,CF:1.6,CS:2.0};
    var MV={H:1,C:4,N:3,O:2,F:1,S:6,P:5,Cl:1,Br:1,I:1,B:3};
    function pk(e1,e2){e1=e1.toUpperCase();e2=e2.toUpperCase();return e1<e2?e1+e2:e2+e1}
    function sk(e1,e2){e1=e1.charAt(0).toUpperCase()+e1.slice(1).toLowerCase();e2=e2.charAt(0).toUpperCase()+e2.slice(1).toLowerCase();return e1<e2?e1+'+'+e2:e2+'+'+e1}
    function gbo(el1,el2,d){
        var p=pk(el1,el2);var co=BC[p];
        if(co!==undefined){if(d>co)return 0}else{var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;if(d>(r1+r2)+0.5)return 0}
        var s=sk(el1,el2);var sp=BS[s];
        if(sp){for(var k=0;k<sp.length;k++){if(Math.abs(d-sp[k].l)<=sp[k].t)return sp[k].o}var bo=1,md=Infinity;for(var k=0;k<sp.length;k++){var df=Math.abs(d-sp[k].l);if(df<md){md=df;bo=sp[k].o}}return bo}
        var r1=CR2[el1.charAt(0).toUpperCase()+el1.slice(1).toLowerCase()]||1.5;var r2=CR2[el2.charAt(0).toUpperCase()+el2.slice(1).toLowerCase()]||1.5;var rs=r1+r2;var ratio=rs?d/rs:1;
        if(ratio<0.85)return 3;if(ratio<0.90)return 2;return 1;
    }
    var n=atoms.length;var bm=new Map();
    for(var i=0;i<n;i++){for(var j=i+1;j<n;j++){
        var dx=atoms[i].x-atoms[j].x,dy=atoms[i].y-atoms[j].y,dz=atoms[i].z-atoms[j].z;
        var d=Math.sqrt(dx*dx+dy*dy+dz*dz);var bo=gbo(atoms[i].element,atoms[j].element,d);
        if(bo>0){if(!bm.has(i))bm.set(i,new Map());if(!bm.has(j))bm.set(j,new Map());bm.get(i).set(j,bo);bm.get(j).set(i,bo)}
    }}
    bm.forEach(function(nb,i){nb.forEach(function(o,j){if(i>j)return;var e1=atoms[i].element.toUpperCase(),e2=atoms[j].element.toUpperCase();var els=new Set([e1,e2]);
        if(els.has('C')&&els.has('O')&&o!==1&&o!==2){var no=o<1.7?1:2;nb.set(j,no);bm.get(j).set(i,no)}
        if(e1==='BR'||e2==='BR'){nb.set(j,1);bm.get(j).set(i,1)}
    })});
    bm.forEach(function(nb,i){var el=atoms[i].element.toUpperCase();if(el!=='N')return;var nl=Array.from(nb.entries());var nn=nl.length;
        if(nn===3){for(var k=0;k<nl.length;k++){nb.set(nl[k][0],1);bm.get(nl[k][0]).set(i,1)}}
        else if(nn===2){var n1=nl[0][0],bo1=nl[0][1],n2=nl[1][0],bo2=nl[1][1];var n1H=atoms[n1].element.toUpperCase()==='H',n2H=atoms[n2].element.toUpperCase()==='H';var f1,f2;
            if(n1H){f1=1;f2=2}else if(n2H){f1=2;f2=1}else{var d1=Math.abs(bo1-1.5)+Math.abs(bo2-1.5),d2=Math.abs(bo1-2)+Math.abs(bo2-1),d3=Math.abs(bo1-1)+Math.abs(bo2-2);var best=Math.min(d1,d2,d3);if(best===d1){f1=1.5;f2=1.5}else if(best===d2){f1=2;f2=1}else{f1=1;f2=2}}
            nb.set(n1,f1);bm.get(n1).set(i,f1);nb.set(n2,f2);bm.get(n2).set(i,f2)}
    });
    for(var iter=0;iter<10;iter++){var changed=false;var val=new Map();for(var i=0;i<n;i++)val.set(i,0);
        bm.forEach(function(nb,i){nb.forEach(function(o){val.set(i,(val.get(i)||0)+o)})});
        var viols=[];val.forEach(function(v,i){var el=atoms[i].element.charAt(0).toUpperCase()+atoms[i].element.slice(1).toLowerCase();var mv=MV[el]||100;if(v>mv+0.1)viols.push(i)});
        if(viols.length===0)break;
        for(var vi=0;vi<viols.length;vi++){var i=viols[vi];var nb=bm.get(i);if(!nb)continue;var cv=Array.from(nb.values()).reduce(function(s,v){return s+v},0);var el=atoms[i].element.charAt(0).toUpperCase()+atoms[i].element.slice(1).toLowerCase();var mv=MV[el]||100;if(cv<=mv+0.1)continue;
            var bb=null,ml=Infinity;nb.forEach(function(o,j){if(o<=1)return;var no;if(o===3)no=2;else if(o===2)no=1.5;else if(o===1.5)no=1;else return;
                var el2=atoms[j].element;var dd=Math.sqrt(Math.pow(atoms[i].x-atoms[j].x,2)+Math.pow(atoms[i].y-atoms[j].y,2)+Math.pow(atoms[i].z-atoms[j].z,2));var s=sk(atoms[i].element,el2);var sp=BS[s];if(!sp)return;
                var ci=0,ni=0,md=Infinity;for(var k=0;k<sp.length;k++){if(sp[k].o===o&&Math.abs(dd-sp[k].l)<md){md=Math.abs(dd-sp[k].l);ci=sp[k].l}}md=Infinity;for(var k=0;k<sp.length;k++){if(sp[k].o===no&&Math.abs(dd-sp[k].l)<md){md=Math.abs(dd-sp[k].l);ni=sp[k].l}}if(!ci||!ni)return;
                var loss=Math.abs(dd-ni)-Math.abs(dd-ci);if(loss<ml){ml=loss;bb=[j,no]}
            });if(bb){nb.set(bb[0],bb[1]);bm.get(bb[0]).set(i,bb[1]);changed=true}
        }if(!changed)break;
    }
    var bonds=[];var seen=new Set();bm.forEach(function(nb,i){nb.forEach(function(o,j){var key=Math.min(i,j)+'-'+Math.max(i,j);if(!seen.has(key)){seen.add(key);bonds.push({atom1:Math.min(i,j),atom2:Math.max(i,j),order:o})}})});
    return bonds;
}
document.getElementById('prev-frame').addEventListener('click',function(){if(currentFrame>0)switchFrame(currentFrame-1)});
document.getElementById('next-frame').addEventListener('click',function(){if(currentFrame<totalFrames-1)switchFrame(currentFrame+1)});
document.getElementById('frame-num').addEventListener('change',function(){var n=parseInt(this.value);if(!isNaN(n)&&n>=1&&n<=totalFrames)switchFrame(n-1)});
document.getElementById('auto-play').addEventListener('click',function(){
    if(isAutoPlaying){
        clearInterval(autoPlayTimer);autoPlayTimer=null;isAutoPlaying=false;
        this.textContent='⏵ Play';this.classList.remove('playing');
    }else{
        isAutoPlaying=true;this.textContent='⏸ Stop';this.classList.add('playing');
        autoPlayTimer=setInterval(function(){
            var next=currentFrame+1;if(next>=totalFrames)next=0;
            switchFrame(next);
        },500);
    }
});

function highlightSelected(){
    atomMeshes.forEach(function(m,i){
        var sel=selectedAtoms.indexOf(i)>=0;
        if(sel){m.material.emissive=new THREE.Color(0xffff00);m.material.emissiveIntensity=0.6}
        else if(m.userData.diffHi){m.material.emissive=new THREE.Color(0xff6600);m.material.emissiveIntensity=0.4}
        else{m.material.emissive=new THREE.Color(0x000000);m.material.emissiveIntensity=0}
    });
    if(diffMode){
        diffAtomMeshes.forEach(function(m,i){
            var sel=diffSelectedAtoms.indexOf(i)>=0;
            if(sel){m.material.emissive=new THREE.Color(0xffff00);m.material.emissiveIntensity=0.6}
            else if(m.userData.diffHi){m.material.emissive=new THREE.Color(0xff6600);m.material.emissiveIntensity=0.4}
            else{m.material.emissive=new THREE.Color(0x000000);m.material.emissiveIntensity=0}
        });
    }
}

function selectAtom(idx){
    if(diffMode&&diffActiveSide==='right'){
        var need=diffSelectedAtoms.length>=requiredCount();
        if(need){diffSelectedAtoms=[idx]}else if(diffSelectedAtoms.indexOf(idx)>=0){return}else{diffSelectedAtoms.push(idx)}
        if(diffMapping){
            var needL=selectedAtoms.length>=requiredCount();
            if(needL){selectedAtoms=[]}
            diffSelectedAtoms.forEach(function(i){
                var li=diffReverseMapping[i];
                if(li>=0&&selectedAtoms.indexOf(li)<0)selectedAtoms.push(li);
            });
        }
        highlightSelected();
        updateDiffSelInfo();
        checkSelectionComplete();
        return;
    }
    if(diffMode){
        var needL=selectedAtoms.length>=requiredCount();
        if(needL){selectedAtoms=[idx]}else if(selectedAtoms.indexOf(idx)>=0){return}else{selectedAtoms.push(idx)}
        if(diffMapping){
            var needR=diffSelectedAtoms.length>=requiredCount();
            if(needR){diffSelectedAtoms=[]}
            selectedAtoms.forEach(function(i){
                var ri=diffMapping[i];
                if(ri>=0&&diffSelectedAtoms.indexOf(ri)<0)diffSelectedAtoms.push(ri);
            });
        }
        highlightSelected();
        updateDiffSelInfo();
        checkSelectionComplete();
        return;
    }
    if(selectedAtoms.indexOf(idx)>=0)return;
    selectedAtoms.push(idx);
    highlightSelected();
    var names=selectedAtoms.map(function(i){return MD.atoms[i].element+(i+1)}).join(', ');
    selInfoEl.textContent='Selected: '+names;
    checkSelectionComplete();
}

function updateDiffSelInfo(){
    var lNames=selectedAtoms.map(function(i){return MD.atoms[i].element+(i+1)}).join(', ');
    var rNames=diffSelectedAtoms.map(function(i){return diffData.atoms[i].element+(i+1)}).join(', ');
    selInfoEl.textContent='[Left] '+(lNames||'-')+'   |   [Right] '+(rNames||'-');
}

function requiredCount(){
    if(currentMode==='bondLength')return 2;
    if(currentMode==='bondAngle')return 3;
    if(currentMode==='dihedral')return 4;
    if(currentMode==='addAtom'||currentMode==='deleteAtom')return 1;
    return 999;
}

function checkSelectionComplete(){
    if(diffMode){
        if(currentMode==='bondLength'){
            if(selectedAtoms.length===2||diffSelectedAtoms.length===2)showDiffMeasurement('bondLength');
        }else if(currentMode==='bondAngle'){
            if(selectedAtoms.length===3||diffSelectedAtoms.length===3)showDiffMeasurement('bondAngle');
        }else if(currentMode==='dihedral'){
            if(selectedAtoms.length===4||diffSelectedAtoms.length===4)showDiffMeasurement('dihedral');
        }else if(currentMode==='selectAtoms'){
            updateDiffSelInfo();
        }
        return;
    }
    if(currentMode==='bondLength'&&selectedAtoms.length===2)showBondLengthModal();
    else if(currentMode==='bondAngle'&&selectedAtoms.length===3)showBondAngleModal();
    else if(currentMode==='dihedral'&&selectedAtoms.length===4)showDihedralModal();
    else if(currentMode==='addAtom'&&selectedAtoms.length===1)showAddAtomModal();
    else if(currentMode==='deleteAtom'&&selectedAtoms.length===1)showDeleteAtomModal();
    else if(currentMode==='selectAtoms'){
        var names=selectedAtoms.map(function(i){return MD.atoms[i].element+(i+1)}).join(', ');
        selInfoEl.textContent='Selected: '+names;
    }
}

function fmtMeasurement(kind,sel,atomsArr){
    if(!sel||sel.length===0)return null;
    if(kind==='bondLength'&&sel.length===2){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]];
        var d=dist(a1,a2);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1),value:d.toFixed(4)+' A',num:d};
    }
    if(kind==='bondAngle'&&sel.length===3){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]];
        var ang=angle(a1,a2,a3);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1)+' - '+a3.element+(sel[2]+1),value:ang.toFixed(2)+' deg',num:ang};
    }
    if(kind==='dihedral'&&sel.length===4){
        var a1=atomsArr[sel[0]],a2=atomsArr[sel[1]],a3=atomsArr[sel[2]],a4=atomsArr[sel[3]];
        var dih=dihedral(a1,a2,a3,a4);
        return {label:a1.element+(sel[0]+1)+' - '+a2.element+(sel[1]+1)+' - '+a3.element+(sel[2]+1)+' - '+a4.element+(sel[3]+1),value:dih.toFixed(2)+' deg',num:dih};
    }
    return null;
}

function showDiffMeasurement(kind){
    var kindLabel=kind==='bondLength'?'Bond Length':kind==='bondAngle'?'Bond Angle':'Dihedral Angle';
    var left=fmtMeasurement(kind,selectedAtoms,MD.atoms);
    var right=fmtMeasurement(kind,diffSelectedAtoms,diffData.atoms);
    var html='<span class="diff-close">×</span><h4>'+kindLabel+' (Diff Mode)</h4>';
    html+='<div class="diff-row" style="color:#4ec9b0"><b>Left: '+(MD.title||'original')+'</b></div>';
    if(left){
        html+='<div class="diff-row">  '+left.label+': '+left.value+'</div>';
    }else{
        html+='<div class="diff-row" style="color:var(--vscode-descriptionForeground,#999)">  Select '+(kind==='bondLength'?2:kind==='bondAngle'?3:4)+' atoms on left...</div>';
    }
    html+='<div class="diff-row" style="margin-top:6px;color:#f0c674"><b>Right: '+diffData.fileName+'</b></div>';
    if(right){
        html+='<div class="diff-row">  '+right.label+': '+right.value+'</div>';
    }else{
        html+='<div class="diff-row" style="color:var(--vscode-descriptionForeground,#999)">  Select '+(kind==='bondLength'?2:kind==='bondAngle'?3:4)+' atoms on right...</div>';
    }
    if(left&&right){
        var diff=Math.abs(left.num-right.num);
        if(kind==='dihedral'&&diff>180)diff=360-diff;
        html+='<div class="diff-row" style="margin-top:6px;color:#ff6600"><b>Δ = '+diff.toFixed(4)+(kind==='bondLength'?' A':' deg')+'</b></div>';
    }
    html+='<div class="diff-row" style="margin-top:6px;color:var(--vscode-descriptionForeground,#999)">Selecting atoms on one side auto-selects the corresponding atoms on the other.</div>';
    showDiffPanel(html);
}

function dist(a,b){var dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;return Math.sqrt(dx*dx+dy*dy+dz*dz)}
function angle(a,b,c){var v1={x:a.x-b.x,y:a.y-b.y,z:a.z-b.z},v2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z};
    var d1=Math.sqrt(v1.x*v1.x+v1.y*v1.y+v1.z*v1.z),d2=Math.sqrt(v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
    if(d1<1e-10||d2<1e-10)return 0;
    var dot=(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(d1*d2);
    return Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI}

function dihedral(a,b,c,d){
    var b1={x:b.x-a.x,y:b.y-a.y,z:b.z-a.z},b2={x:c.x-b.x,y:c.y-b.y,z:c.z-b.z},b3={x:d.x-c.x,y:d.y-c.y,z:d.z-c.z};
    function cross(u,v){return{x:u.y*v.z-u.z*v.y,y:u.z*v.x-u.x*v.z,z:u.x*v.y-u.y*v.x}}
    function dot(u,v){return u.x*v.x+u.y*v.y+u.z*v.z}
    function norm(v){var l=Math.sqrt(v.x*v.x+v.y*v.y+v.z*v.z);return l<1e-10?{x:0,y:0,z:0}:{x:v.x/l,y:v.y/l,z:v.z/l}}
    var n1=cross(b1,b2),n2=cross(b2,b3);
    var m=cross(norm(n1),norm(b2));
    var x=dot(n1,n2),y=dot(m,n2);
    return-Math.atan2(y,x)*180/Math.PI;
}

function getMovable(fixedSet,startIdx){
    var adj=[];for(var i=0;i<MD.atoms.length;i++)adj[i]=[];
    MD.bonds.forEach(function(b){adj[b.atom1].push(b.atom2);adj[b.atom2].push(b.atom1)});
    var visited=new Set(fixedSet);visited.add(startIdx);
    var queue=[startIdx],result=[startIdx];
    while(queue.length>0){var cur=queue.shift();adj[cur].forEach(function(nb){
        if(!visited.has(nb)){visited.add(nb);queue.push(nb);result.push(nb)}})}
    return result;
}

function rotAroundAxis(px,py,pz,ox,oy,oz,dx,dy,dz,angle){
    var c=Math.cos(angle),s=Math.sin(angle);
    var x=px-ox,y=py-oy,z=pz-oz;
    var kx=dx,ky=dy,kz=dz;
    var l=Math.sqrt(kx*kx+ky*ky+kz*kz);if(l<1e-10)return{x:px,y:py,z:pz};
    kx/=l;ky/=l;kz/=l;
    var dot=x*kx+y*ky+z*kz;
    var rx=x*c+(ky*z-kz*y)*s+kx*dot*(1-c);
    var ry=y*c+(kz*x-kx*z)*s+ky*dot*(1-c);
    var rz=z*c+(kx*y-ky*x)*s+kz*dot*(1-c);
    return{x:rx+ox,y:ry+oy,z:rz+oz};
}

function saveOriginal(){originalCoords=MD.atoms.map(function(a){return{x:a.x,y:a.y,z:a.z}})}
function restoreOriginal(){if(!originalCoords)return;originalCoords.forEach(function(c,i){MD.atoms[i].x=c.x;MD.atoms[i].y=c.y;MD.atoms[i].z=c.z})}

function propagateToAllCells(){
    if(!CRY||!originalCoords)return;
    var modified=[];
    MD.atoms.forEach(function(a,i){
        var oc=originalCoords[i];
        if(oc&&(Math.abs(a.x-oc.x)>1e-12||Math.abs(a.y-oc.y)>1e-12||Math.abs(a.z-oc.z)>1e-12)){
            modified.push(i);
        }
    });
    var modifiedSet=new Set(modified);
    modified.forEach(function(idx){
        var a=MD.atoms[idx];
        if(a.baseIdx==null)return;
        var oc=originalCoords[idx];
        var dx=a.x-oc.x,dy=a.y-oc.y,dz=a.z-oc.z;
        MD.atoms.forEach(function(other,i){
            if(modifiedSet.has(i))return;
            if(other.baseIdx===a.baseIdx){
                var oc2=originalCoords[i];
                other.x=oc2.x+dx;
                other.y=oc2.y+dy;
                other.z=oc2.z+dz;
                modifiedSet.add(i);
            }
        });
    });
}

function applyBondLength(targetLen,fixFirst){
    var i1=selectedAtoms[0],i2=selectedAtoms[1];
    var a1=MD.atoms[i1],a2=MD.atoms[i2];
    var dx=a2.x-a1.x,dy=a2.y-a1.y,dz=a2.z-a1.z;
    var curLen=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(curLen<1e-10)return;
    var nx=dx/curLen,ny=dy/curLen,nz=dz/curLen;
    if(fixFirst){
        var fixedSet=new Set([i1]);var movable=getMovable(fixedSet,i2);
        movable.forEach(function(idx){var a=MD.atoms[idx];
            var ox=originalCoords[idx].x,oy=originalCoords[idx].y,oz=originalCoords[idx].z;
            var vx=ox-a1.x,vy=oy-a1.y,vz=oz-a1.z;
            var proj=vx*nx+vy*ny+vz*nz;
            var scale=targetLen/curLen;
            a.x=a1.x+vx*scale;a.y=a1.y+vy*scale;a.z=a1.z+vz*scale;
        });
    }else{
        var fixedSet=new Set([i2]);var movable=getMovable(fixedSet,i1);
        movable.forEach(function(idx){var a=MD.atoms[idx];
            var ox=originalCoords[idx].x,oy=originalCoords[idx].y,oz=originalCoords[idx].z;
            var vx=ox-a2.x,vy=oy-a2.y,vz=oz-a2.z;
            var nx2=a1.x-a2.x,ny2=a1.y-a2.y,nz2=a1.z-a2.z;
            var curLen2=Math.sqrt(nx2*nx2+ny2*ny2+nz2*nz2);
            if(curLen2<1e-10)return;
            nx2/=curLen2;ny2/=curLen2;nz2/=curLen2;
            var scale=targetLen/curLen2;
            a.x=a2.x+vx*scale;a.y=a2.y+vy*scale;a.z=a2.z+vz*scale;
        });
    }
    if(CRY)propagateToAllCells();
    updateScenePositions();
}

function applyBondAngle(targetDeg,fixFirstTwo){
    var i1=selectedAtoms[0],i2=selectedAtoms[1],i3=selectedAtoms[2];
    var a1=MD.atoms[i1],a2=MD.atoms[i2],a3=MD.atoms[i3];
    var curDeg=angle(originalCoords[i1],originalCoords[i2],originalCoords[i3]);
    var delta=(targetDeg-curDeg)*Math.PI/180;
    if(Math.abs(delta)<1e-10)return;
    var v1={x:originalCoords[i1].x-originalCoords[i2].x,y:originalCoords[i1].y-originalCoords[i2].y,z:originalCoords[i1].z-originalCoords[i2].z};
    var v2={x:originalCoords[i3].x-originalCoords[i2].x,y:originalCoords[i3].y-originalCoords[i2].y,z:originalCoords[i3].z-originalCoords[i2].z};
    var cx=v1.y*v2.z-v1.z*v2.y,cy=v1.z*v2.x-v1.x*v2.z,cz=v1.x*v2.y-v1.y*v2.x;
    var cl=Math.sqrt(cx*cx+cy*cy+cz*cz);if(cl<1e-10)return;
    cx/=cl;cy/=cl;cz/=cl;
    var ox=originalCoords[i2].x,oy=originalCoords[i2].y,oz=originalCoords[i2].z;
    if(fixFirstTwo){
        var fixedSet=new Set([i1,i2]);var movable=getMovable(fixedSet,i3);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,cx,cy,cz,delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }else{
        var fixedSet=new Set([i2,i3]);var movable=getMovable(fixedSet,i1);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,cx,cy,cz,-delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }
    if(CRY)propagateToAllCells();
    updateScenePositions();
}

function applyDihedral(targetDeg,fixFirstThree){
    var i1=selectedAtoms[0],i2=selectedAtoms[1],i3=selectedAtoms[2],i4=selectedAtoms[3];
    var curDeg=dihedral(originalCoords[i1],originalCoords[i2],originalCoords[i3],originalCoords[i4]);
    var delta=(targetDeg-curDeg)*Math.PI/180;
    if(Math.abs(delta)<1e-10)return;
    var ax=originalCoords[i3].x-originalCoords[i2].x,ay=originalCoords[i3].y-originalCoords[i2].y,az=originalCoords[i3].z-originalCoords[i2].z;
    var al=Math.sqrt(ax*ax+ay*ay+az*az);if(al<1e-10)return;
    ax/=al;ay/=al;az/=al;
    var ox=originalCoords[i2].x,oy=originalCoords[i2].y,oz=originalCoords[i2].z;
    if(fixFirstThree){
        var fixedSet=new Set([i1,i2,i3]);var movable=getMovable(fixedSet,i4);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,ax,ay,az,delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }else{
        var fixedSet=new Set([i2,i3,i4]);var movable=getMovable(fixedSet,i1);
        movable.forEach(function(idx){
            var oc=originalCoords[idx];
            var r=rotAroundAxis(oc.x,oc.y,oc.z,ox,oy,oz,ax,ay,az,-delta);
            MD.atoms[idx].x=r.x;MD.atoms[idx].y=r.y;MD.atoms[idx].z=r.z;
        });
    }
    if(CRY)propagateToAllCells();
    updateScenePositions();
}

function showModal(html,cb){modalEl.innerHTML=html;modalOverlay.classList.add('show');modalCallback=cb}
function hideModal(){modalOverlay.classList.remove('show');modalCallback=null}

function showBondLengthModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]];
    var cur=dist(a1,a2);
    var existingBond=MD.bonds.find(function(b){return(b.atom1===selectedAtoms[0]&&b.atom2===selectedAtoms[1])||(b.atom1===selectedAtoms[1]&&b.atom2===selectedAtoms[0])});
    var curOrder=existingBond?existingBond.order:1;
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1);
    showModal('<h3>Adjust Bond Length</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(4)+' A, Bond order: '+curOrder+'</div>'+
        '<label>Fix atom:</label><select id="m-fix"><option value="1">Fix '+n1+' (move '+n2+')</option><option value="2">Fix '+n2+' (move '+n1+')</option></select>'+
        '<label>Bond order:</label><select id="m-order"><option value="0"'+(curOrder===0?' selected':'')+'>None (0) - Remove bond</option><option value="1"'+(curOrder===1?' selected':'')+'>Single (1.0)</option><option value="1.5"'+(curOrder===1.5?' selected':'')+'>Aromatic (1.5)</option><option value="2"'+(curOrder===2?' selected':'')+'>Double (2.0)</option><option value="3"'+(curOrder===3?' selected':'')+'>Triple (3.0)</option></select>'+
        '<label>Target length (A):</label><input type="number" id="m-val" value="'+cur.toFixed(4)+'" step="0.01" min="0.3" max="6">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(4)+'" min="0.3" max="6" step="0.01">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix'),orderEl=document.getElementById('m-order');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyBondLength(parseFloat(this.value),fixEl.value==='1')});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyBondLength(parseFloat(this.value),fixEl.value==='1')});
    fixEl.addEventListener('change',function(){applyBondLength(parseFloat(valEl.value),this.value==='1')});
    orderEl.addEventListener('change',function(){
        var newOrder=parseFloat(this.value);
        if(newOrder===0){
            if(existingBond){MD.bonds=MD.bonds.filter(function(b){return b!==existingBond});existingBond=null;rebuildScene()}
        }else{
            if(!existingBond){
                existingBond={atom1:selectedAtoms[0],atom2:selectedAtoms[1],order:newOrder};
                MD.bonds.push(existingBond);
            }else{existingBond.order=newOrder}
            rebuildScene();
        }
    });
    document.getElementById('m-ok').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){
        if(existingBond&&originalCoords){existingBond.order=curOrder}
        restoreOriginal();rebuildScene();hideModal();originalCoords=null;resetSelection()
    });
}

function showBondAngleModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]];
    var cur=angle(a1,a2,a3);
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1),n3=a3.element+(selectedAtoms[2]+1);
    showModal('<h3>Adjust Bond Angle</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(2)+' deg</div>'+
        '<label>Fix side:</label><select id="m-fix"><option value="1">Fix '+n1+'-'+n2+' (move '+n3+')</option><option value="2">Fix '+n2+'-'+n3+' (move '+n1+')</option></select>'+
        '<label>Target angle (deg):</label><input type="number" id="m-val" value="'+cur.toFixed(2)+'" step="0.5" min="5" max="175">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(2)+'" min="5" max="175" step="0.5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyBondAngle(parseFloat(this.value),fixEl.value==='1')});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyBondAngle(parseFloat(this.value),fixEl.value==='1')});
    fixEl.addEventListener('change',function(){applyBondAngle(parseFloat(valEl.value),this.value==='1')});
    document.getElementById('m-ok').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){restoreOriginal();rebuildScene();hideModal();originalCoords=null;resetSelection()});
}

function showDihedralModal(){
    var a1=MD.atoms[selectedAtoms[0]],a2=MD.atoms[selectedAtoms[1]],a3=MD.atoms[selectedAtoms[2]],a4=MD.atoms[selectedAtoms[3]];
    var cur=dihedral(a1,a2,a3,a4);
    saveOriginal();
    var n1=a1.element+(selectedAtoms[0]+1),n2=a2.element+(selectedAtoms[1]+1),n3=a3.element+(selectedAtoms[2]+1),n4=a4.element+(selectedAtoms[3]+1);
    showModal('<h3>Adjust Dihedral Angle</h3>'+
        '<div class="current-val">Current: '+cur.toFixed(2)+' deg</div>'+
        '<label>Fix side:</label><select id="m-fix"><option value="1">Fix '+n1+'-'+n2+'-'+n3+' (move '+n4+')</option><option value="2">Fix '+n2+'-'+n3+'-'+n4+' (move '+n1+')</option></select>'+
        '<label>Target dihedral (deg):</label><input type="number" id="m-val" value="'+cur.toFixed(2)+'" step="1" min="-180" max="180">'+
        '<input type="range" id="m-slider" value="'+cur.toFixed(2)+'" min="-180" max="180" step="0.5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    var valEl=document.getElementById('m-val'),sliderEl=document.getElementById('m-slider'),fixEl=document.getElementById('m-fix');
    sliderEl.addEventListener('input',function(){valEl.value=this.value;applyDihedral(parseFloat(this.value),fixEl.value==='1')});
    valEl.addEventListener('input',function(){sliderEl.value=this.value;applyDihedral(parseFloat(this.value),fixEl.value==='1')});
    fixEl.addEventListener('change',function(){applyDihedral(parseFloat(valEl.value),this.value==='1')});
    document.getElementById('m-ok').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
    document.getElementById('m-cancel').addEventListener('click',function(){restoreOriginal();rebuildScene();hideModal();originalCoords=null;resetSelection()});
}

function showAddAtomModal(){
    var anchorIdx=selectedAtoms[0];
    var anchor=MD.atoms[anchorIdx];
    saveOriginal();
    showModal('<h3>Add Atom</h3>'+
        '<label>Element:</label><select id="m-elem">'+
        '<option>H</option><option>He</option>'+
        '<option>Li</option><option>Be</option><option>B</option><option>C</option><option>N</option><option>O</option><option>F</option><option>Ne</option>'+
        '<option>Na</option><option>Mg</option><option>Al</option><option>Si</option><option>P</option><option>S</option><option>Cl</option><option>Ar</option>'+
        '<option>K</option><option>Ca</option><option>Sc</option><option>Ti</option><option>V</option><option>Cr</option><option>Mn</option><option>Fe</option>'+
        '<option>Co</option><option>Ni</option><option>Cu</option><option>Zn</option><option>Ga</option><option>Ge</option><option>As</option><option>Se</option>'+
        '<option>Br</option><option>Kr</option><option>Rb</option><option>Sr</option><option>Y</option><option>Zr</option><option>Nb</option><option>Mo</option>'+
        '<option>Ru</option><option>Rh</option><option>Pd</option><option>Ag</option><option>Cd</option><option>In</option><option>Sn</option><option>Sb</option>'+
        '<option>Te</option><option>I</option><option>Xe</option><option>Cs</option><option>Ba</option><option>La</option><option>Ce</option><option>Pr</option>'+
        '<option>Nd</option><option>Sm</option><option>Eu</option><option>Gd</option><option>Tb</option><option>Dy</option><option>Ho</option><option>Er</option>'+
        '<option>Tm</option><option>Yb</option><option>Lu</option><option>Hf</option><option>Ta</option><option>W</option><option>Re</option><option>Os</option>'+
        '<option>Ir</option><option>Pt</option><option>Au</option><option>Hg</option><option>Tl</option><option>Pb</option><option>Bi</option>'+
        '</select>'+
        '<label>Bond order:</label><select id="m-bond-order"><option value="1">Single (1.0)</option><option value="1.5">Aromatic (1.5)</option><option value="2">Double (2.0)</option><option value="3">Triple (3.0)</option></select>'+
        '<label>Bond length (A):</label><input type="number" id="m-val" value="1.09" step="0.01" min="0.3" max="5">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">OK</button></div>',null);
    document.getElementById('m-elem').addEventListener('change',function(){
        var defaults={H:1.09,He:1.30,Li:2.12,Be:1.67,B:1.58,C:1.54,N:1.47,O:1.43,F:1.36,Ne:1.35,
            Na:1.88,Mg:1.63,Al:1.84,Si:1.87,P:1.80,S:1.82,Cl:1.77,Ar:1.74,
            K:2.34,Ca:2.00,Sc:1.88,Ti:1.87,V:1.79,Cr:1.79,Mn:1.83,Fe:1.80,
            Co:1.78,Ni:1.73,Cu:1.84,Zn:1.88,Ga:1.87,Ge:1.88,As:1.87,Se:1.90,
            Br:1.94,Kr:1.90,Rb:2.53,Sr:2.15,Y:2.12,Zr:2.06,Nb:2.04,Mo:2.08,
            Ru:2.07,Rh:2.09,Pd:2.05,Ag:2.10,Cd:2.07,In:2.10,Sn:2.17,Sb:2.12,
            Te:2.14,I:2.14,Xe:2.16,Cs:2.65,Ba:2.22,La:2.32,Ce:2.30,Pr:2.31,
            Nd:2.30,Sm:2.29,Eu:2.29,Gd:2.29,Tb:2.28,Dy:2.28,Ho:2.27,Er:2.27,
            Tm:2.26,Yb:2.26,Lu:2.26,Hf:2.23,Ta:2.22,W:2.18,Re:2.20,Os:2.19,
            Ir:2.16,Pt:2.13,Au:2.14,Hg:2.14,Tl:2.20,Pb:2.22,Bi:2.23};
        document.getElementById('m-val').value=defaults[this.value]||1.5;
    });
    document.getElementById('m-ok').addEventListener('click',function(){
        var el=document.getElementById('m-elem').value;
        var bl=parseFloat(document.getElementById('m-val').value)||1.5;
        var bondOrder=parseFloat(document.getElementById('m-bond-order').value)||1;
        var dir={x:0,y:0,z:1};
        if(CRY){
            var baseAnchorIdx=anchor.baseIdx!=null?anchor.baseIdx:anchorIdx;
            var baseAnchor=CRY.baseAtoms[baseAnchorIdx];
            var bondedBase=[];
            CRY.baseBonds.forEach(function(b){
                if(b.atom1===baseAnchorIdx)bondedBase.push(b.atom2);
                if(b.atom2===baseAnchorIdx)bondedBase.push(b.atom1);
            });
            if(bondedBase.length>0){
                var avg={x:0,y:0,z:0};
                bondedBase.forEach(function(bi){avg.x+=CRY.baseAtoms[bi].x-baseAnchor.x;avg.y+=CRY.baseAtoms[bi].y-baseAnchor.y;avg.z+=CRY.baseAtoms[bi].z-baseAnchor.z});
                var al=Math.sqrt(avg.x*avg.x+avg.y*avg.y+avg.z*avg.z);
                if(al>1e-10){dir={x:-avg.x/al,y:-avg.y/al,z:-avg.z/al}}
            }
            var newBaseIdx=CRY.baseAtoms.length;
            CRY.baseAtoms.push({element:el,x:baseAnchor.x+dir.x*bl,y:baseAnchor.y+dir.y*bl,z:baseAnchor.z+dir.z*bl,index:newBaseIdx,baseIdx:newBaseIdx,occupancy:1});
            CRY.baseBonds.push({atom1:baseAnchorIdx,atom2:newBaseIdx,order:bondOrder});
            rebuildCrystal();rebuildScene();hideModal();originalCoords=null;resetSelection();
        }else{
            var bonded=[];
            MD.bonds.forEach(function(b){
                if(b.atom1===anchorIdx)bonded.push(b.atom2);
                if(b.atom2===anchorIdx)bonded.push(b.atom1);
            });
            if(bonded.length>0){
                var avg={x:0,y:0,z:0};
                bonded.forEach(function(bi){avg.x+=MD.atoms[bi].x-anchor.x;avg.y+=MD.atoms[bi].y-anchor.y;avg.z+=MD.atoms[bi].z-anchor.z});
                var al=Math.sqrt(avg.x*avg.x+avg.y*avg.y+avg.z*avg.z);
                if(al>1e-10){dir={x:-avg.x/al,y:-avg.y/al,z:-avg.z/al}}
            }
            var newIdx=MD.atoms.length;
            MD.atoms.push({element:el,x:anchor.x+dir.x*bl,y:anchor.y+dir.y*bl,z:anchor.z+dir.z*bl,color:MD.atomColors[el]||'#FF1493',index:newIdx});
            MD.bonds.push({atom1:anchorIdx,atom2:newIdx,order:bondOrder});
            rebuildScene();hideModal();originalCoords=null;resetSelection();
        }
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();originalCoords=null;resetSelection()});
}

function showDeleteAtomModal(){
    var idx=selectedAtoms[0];
    var a=MD.atoms[idx];
    var name=a.element+(idx+1);
    showModal('<h3>Delete Atom</h3>'+
        '<div class="current-val">Delete '+name+'?</div>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok mbtn-danger" id="m-ok">Delete</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        if(CRY){
            var baseIdx=a.baseIdx!=null?a.baseIdx:idx;
            CRY.baseAtoms.splice(baseIdx,1);
            CRY.baseAtoms.forEach(function(ba,i){ba.index=i;ba.baseIdx=i});
            CRY.baseBonds=CRY.baseBonds.filter(function(b){return b.atom1!==baseIdx&&b.atom2!==baseIdx}).map(function(b){
                return{atom1:b.atom1>baseIdx?b.atom1-1:b.atom1,atom2:b.atom2>baseIdx?b.atom2-1:b.atom2,order:b.order,shift:b.shift};
            });
            rebuildCrystal();rebuildScene();hideModal();resetSelection();
        }else{
            MD.atoms.splice(idx,1);
            MD.atoms.forEach(function(a,i){a.index=i});
            MD.bonds=MD.bonds.filter(function(b){return b.atom1!==idx&&b.atom2!==idx}).map(function(b){
                return{atom1:b.atom1>idx?b.atom1-1:b.atom1,atom2:b.atom2>idx?b.atom2-1:b.atom2,order:b.order};
            });
            rebuildScene();hideModal();resetSelection();
        }
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal();resetSelection()});
}

function doSave(){
    var xyz=MD.atoms.length+'\\n'+(MD.title||'Modified structure')+'\\n';
    MD.atoms.forEach(function(a){xyz+=a.element+'  '+a.x.toFixed(6)+'  '+a.y.toFixed(6)+'  '+a.z.toFixed(6)+'\\n'});
    var gjf='';
    var meta=MD.gjfMeta;
    if(meta){
        meta.link0Lines.forEach(function(l){gjf+=l+'\\n'});
        gjf+=meta.routeLine+'\\n\\n';
        meta.titleLines.forEach(function(l){gjf+=l+'\\n'});
        gjf+='\\n'+meta.chargeMultLine+'\\n';
    }else{
        gjf='%chk=molecule.chk\\n%mem=4GB\\n%nproc=4\\n# B3LYP/6-31G(d)\\n\\n'+(MD.title||'Modified structure')+'\\n\\n0 1\\n';
    }
    MD.atoms.forEach(function(a){
        var el=a.element; while(el.length<2) el+=' ';
        var sx=a.x.toFixed(8); while(sx.length<17) sx=' '+sx;
        var sy=a.y.toFixed(8); while(sy.length<17) sy=' '+sy;
        var sz=a.z.toFixed(8); while(sz.length<17) sz=' '+sz;
        gjf+=' '+el+' '+sx+' '+sy+' '+sz+'\\n';
    });
    gjf+='\\n';
    MD.atoms.forEach(function(a,i){
        var parts=[i+1];
        MD.bonds.forEach(function(b){
            if(b.atom1===i)parts.push(b.atom2+1,b.order.toFixed(1));
        });
        gjf+=parts.join(' ')+'\\n';
    });
    gjf+='\\n';
    if(meta&&meta.afterConnectContent){gjf+=meta.afterConnectContent+'\\n'}

    var chrg=MD.charge||0;
    var mult=MD.multiplicity||1;
    var AN2={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Ne:10,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Sc:21,Ti:22,V:23,Cr:24,Mn:25,Fe:26,Co:27,Ni:28,Cu:29,Zn:30,Ga:31,Ge:32,As:33,Se:34,Br:35,Kr:36,Rb:37,Sr:38,Y:39,Zr:40,Nb:41,Mo:42,Tc:43,Ru:44,Rh:45,Pd:46,Ag:47,Cd:48,In:49,Sn:50,Sb:51,Te:52,I:53,Xe:54,Cs:55,Ba:56,La:57,Ce:58,Pr:59,Nd:60,Pm:61,Sm:62,Eu:63,Gd:64,Tb:65,Dy:66,Ho:67,Er:68,Tm:69,Yb:70,Lu:71,Hf:72,Ta:73,W:74,Re:75,Os:76,Ir:77,Pt:78,Au:79,Hg:80,Tl:81,Pb:82,Bi:83,Po:84,At:85,Rn:86};

    var coord='$coord\\n';
    var ANG_TO_BOHR=1.8897259886;
    MD.atoms.forEach(function(a){coord+='  '+(a.x*ANG_TO_BOHR).toFixed(8)+'  '+(a.y*ANG_TO_BOHR).toFixed(8)+'  '+(a.z*ANG_TO_BOHR).toFixed(8)+' '+a.element.toLowerCase()+'\\n'});
    coord+='$end\\n';
    if(chrg!==0)coord+='$chrg '+chrg+'\\n';
    if(mult!==1)coord+='$spin '+(mult-1)/2+'\\n';

    var orcaInp='! B3LYP def2-SVP\\n\\n* xyz '+chrg+' '+mult+'\\n';
    MD.atoms.forEach(function(a){orcaInp+=a.element+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+'\\n'});
    orcaInp+='*\\n';

    var mol2='@<TRIPOS>MOLECULE\\n'+(MD.title||'Modified structure')+'\\n'+MD.atoms.length+' '+MD.bonds.length+' 0 0 0\\nSMALL\\nNO_CHARGES\\n\\n';
    mol2+='@<TRIPOS>ATOM\\n';
    MD.atoms.forEach(function(a,i){mol2+=(i+1)+' '+a.element+' '+a.x.toFixed(6)+' '+a.y.toFixed(6)+' '+a.z.toFixed(6)+' '+a.element+' 1 UNK 0.000\\n'});
    mol2+='@<TRIPOS>BOND\\n';
    var bondIdx=1;
    MD.bonds.forEach(function(b){
        var bt='1';
        if(b.order>=2.5)bt='3';
        else if(b.order>=1.75)bt='2';
        else if(b.order>=1.25)bt='ar';
        mol2+=bondIdx+' '+(b.atom1+1)+' '+(b.atom2+1)+' '+bt+'\\n';
        bondIdx++;
    });

    var mol='\\n '+MD.atoms.length+'  '+MD.bonds.length+'  0  0  0  0  0  0  0  0999 V2000\\n';
    MD.atoms.forEach(function(a){
        var z=AN2[a.element]||0;
        var sx=(a.x*10).toFixed(4);
        var sy=(a.y*10).toFixed(4);
        var sz=(a.z*10).toFixed(4);
        while(sx.length<10)sx=' '+sx;
        while(sy.length<10)sy=' '+sy;
        while(sz.length<10)sz=' '+sz;
        mol+=sx+sy+sz+' '+a.element+' 0  0  0  0  0  0  0  0  0  0  0  0\\n';
    });
    MD.bonds.forEach(function(b){
        var bt=1;
        if(b.order>=2.5)bt=3;
        else if(b.order>=1.75)bt=2;
        else if(b.order>=1.25)bt=4;
        var a1=b.atom1+1,a2=b.atom2+1;
        var s1=''+a1,s2=''+a2,s3=''+bt;
        while(s1.length<3)s1=' '+s1;
        while(s2.length<3)s2=' '+s2;
        while(s3.length<3)s3=' '+s3;
        mol+=s1+s2+s3+'  0  0  0  0\\n';
    });
    mol+='M  END\\n';

    var pdb='';
    MD.atoms.forEach(function(a,i){
        var serial=(i+1).toString();while(serial.length<5)serial=' '+serial;
        var name;
        if(a.element.length===1){name=' '+a.element+'  '}
        else{name=a.element+'  '}
        var sx=a.x.toFixed(3);while(sx.length<8)sx=' '+sx;
        var sy=a.y.toFixed(3);while(sy.length<8)sy=' '+sy;
        var sz=a.z.toFixed(3);while(sz.length<8)sz=' '+sz;
        var el=a.element.charAt(0).toUpperCase();
        if(a.element.length>1)el+=a.element.charAt(1).toLowerCase();
        while(el.length<2)el=' '+el;
        pdb+='ATOM  '+serial+name+'UNK A   1    '+sx+sy+sz+'  1.00  0.00          '+el+'\\n';
    });
    MD.bonds.forEach(function(b){
        var a1=(b.atom1+1).toString();while(a1.length<5)a1=' '+a1;
        var a2=(b.atom2+1).toString();while(a2.length<5)a2=' '+a2;
        pdb+='CONECT'+a1+a2+'\\n';
    });
    pdb+='END\\n';

    var mopac=(MD.title||'Modified structure')+'\\n';
    var chrgStr=chrg!==0?' CHARGE='+chrg:'';
    var multStr=mult!==1?' MS='+mult:'';
    mopac+='PM7'+chrgStr+multStr+'\\n\\n';
    MD.atoms.forEach(function(a){
        var sx=a.x.toFixed(5);var sy=a.y.toFixed(5);var sz=a.z.toFixed(5);
        mopac+=a.element+' '+sx+' 1 '+sy+' 1 '+sz+' 1\\n';
    });
    mopac+='\\n';

    var cifContent='';
    if(CRY){
        var lv=CRY.latticeVectors;
        var det=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
        var invDet=1/det;
        var inv=[
            [(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])*invDet,(lv[0][2]*lv[2][1]-lv[0][1]*lv[2][2])*invDet,(lv[0][1]*lv[1][2]-lv[0][2]*lv[1][1])*invDet],
            [(lv[1][2]*lv[2][0]-lv[1][0]*lv[2][2])*invDet,(lv[0][0]*lv[2][2]-lv[0][2]*lv[2][0])*invDet,(lv[0][2]*lv[1][0]-lv[0][0]*lv[1][2])*invDet],
            [(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0])*invDet,(lv[0][1]*lv[2][0]-lv[0][0]*lv[2][1])*invDet,(lv[0][0]*lv[1][1]-lv[0][1]*lv[1][0])*invDet]
        ];
        cifContent='data_modified\\n';
        cifContent+='_cell_length_a '+CRY.a.toFixed(6)+'\\n';
        cifContent+='_cell_length_b '+CRY.b.toFixed(6)+'\\n';
        cifContent+='_cell_length_c '+CRY.c.toFixed(6)+'\\n';
        cifContent+='_cell_angle_alpha '+CRY.alpha.toFixed(6)+'\\n';
        cifContent+='_cell_angle_beta '+CRY.beta.toFixed(6)+'\\n';
        cifContent+='_cell_angle_gamma '+CRY.gamma.toFixed(6)+'\\n';
        if(CRY.spaceGroup)cifContent+='_symmetry_space_group_name_H-M \\''+CRY.spaceGroup+'\\'\\n';
        cifContent+='loop_\\n_symmetry_equiv_pos_as_xyz\\n';
        CRY.symmetryOps.forEach(function(op){cifContent+='\\''+op+'\\'\\n'});
        cifContent+='loop_\\n_atom_site_label\\n_atom_site_type_symbol\\n_atom_site_fract_x\\n_atom_site_fract_y\\n_atom_site_fract_z\\n_atom_site_occupancy\\n';
        CRY.baseAtoms.forEach(function(ba,i){
            var fx=inv[0][0]*ba.x+inv[0][1]*ba.y+inv[0][2]*ba.z;
            var fy=inv[1][0]*ba.x+inv[1][1]*ba.y+inv[1][2]*ba.z;
            var fz=inv[2][0]*ba.x+inv[2][1]*ba.y+inv[2][2]*ba.z;
            var label=ba.element+(i+1);
            var occ=(ba.occupancy!=null?ba.occupancy:1).toFixed(3);
            cifContent+=label+' '+ba.element+' '+fx.toFixed(6)+' '+fy.toFixed(6)+' '+fz.toFixed(6)+' '+occ+'\\n';
        });
    }

    var vaspContent='';
    if(CRY){
        var lv=CRY.latticeVectors;
        var vDet=lv[0][0]*(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])-lv[0][1]*(lv[1][0]*lv[2][2]-lv[1][2]*lv[2][0])+lv[0][2]*(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0]);
        var vInvDet=1/vDet;
        var vInv=[
            [(lv[1][1]*lv[2][2]-lv[1][2]*lv[2][1])*vInvDet,(lv[0][2]*lv[2][1]-lv[0][1]*lv[2][2])*vInvDet,(lv[0][1]*lv[1][2]-lv[0][2]*lv[1][1])*vInvDet],
            [(lv[1][2]*lv[2][0]-lv[1][0]*lv[2][2])*vInvDet,(lv[0][0]*lv[2][2]-lv[0][2]*lv[2][0])*vInvDet,(lv[0][2]*lv[1][0]-lv[0][0]*lv[1][2])*vInvDet],
            [(lv[1][0]*lv[2][1]-lv[1][1]*lv[2][0])*vInvDet,(lv[0][1]*lv[2][0]-lv[0][0]*lv[2][1])*vInvDet,(lv[0][0]*lv[1][1]-lv[0][1]*lv[1][0])*vInvDet]
        ];
        var elemOrder=[];
        var elemCounts={};
        MD.atoms.forEach(function(a){
            if(elemCounts[a.element]===undefined){elemOrder.push(a.element);elemCounts[a.element]=0}
            elemCounts[a.element]++;
        });
        vaspContent=(MD.title||'Modified structure')+'\\n';
        vaspContent+='1.0\\n';
        for(var ii=0;ii<3;ii++){
            vaspContent+='  '+lv[ii][0].toFixed(8)+'  '+lv[ii][1].toFixed(8)+'  '+lv[ii][2].toFixed(8)+'\\n';
        }
        vaspContent+=elemOrder.join('  ')+'\\n';
        vaspContent+=elemOrder.map(function(e){return elemCounts[e]}).join('  ')+'\\n';
        vaspContent+='Direct\\n';
        MD.atoms.forEach(function(a){
            var fx=vInv[0][0]*a.x+vInv[0][1]*a.y+vInv[0][2]*a.z;
            var fy=vInv[1][0]*a.x+vInv[1][1]*a.y+vInv[1][2]*a.z;
            var fz=vInv[2][0]*a.x+vInv[2][1]*a.y+vInv[2][2]*a.z;
            vaspContent+='  '+fx.toFixed(8)+'  '+fy.toFixed(8)+'  '+fz.toFixed(8)+'\\n';
        });
    }

    var cubeContent='';
    if(CRY){
        var clv=CRY.latticeVectors;
        var nGrid=2;
        cubeContent='Molecular Viewer export\\n';
        cubeContent+='Generated from '+(MD.title||'structure')+'\\n';
        cubeContent+='-'+MD.atoms.length+'  0.000000  0.000000  0.000000\\n';
        cubeContent+=nGrid+'  '+(clv[0][0]/nGrid).toFixed(8)+'  '+(clv[0][1]/nGrid).toFixed(8)+'  '+(clv[0][2]/nGrid).toFixed(8)+'\\n';
        cubeContent+=nGrid+'  '+(clv[1][0]/nGrid).toFixed(8)+'  '+(clv[1][1]/nGrid).toFixed(8)+'  '+(clv[1][2]/nGrid).toFixed(8)+'\\n';
        cubeContent+=nGrid+'  '+(clv[2][0]/nGrid).toFixed(8)+'  '+(clv[2][1]/nGrid).toFixed(8)+'  '+(clv[2][2]/nGrid).toFixed(8)+'\\n';
        MD.atoms.forEach(function(a){
            var z=AN2[a.element]||0;
            cubeContent+=z+'  0.000000  '+a.x.toFixed(8)+'  '+a.y.toFixed(8)+'  '+a.z.toFixed(8)+'\\n';
        });
        var nVals=nGrid*nGrid*nGrid;
        for(var vi=0;vi<nVals;vi++){
            cubeContent+='0.0000e+00';
            if((vi+1)%6===0||vi===nVals-1){cubeContent+='\\n'}else{cubeContent+='  '}
        }
    }

    showModal('<h3>Save File</h3>'+
        '<label>Format:</label><select id="m-fmt">'+
        '<option value="xyz">XYZ (.xyz)</option>'+
        '<option value="gjf">Gaussian Input (.gjf)</option>'+
        (CRY?'<option value="cif">CIF (.cif)</option>':'')+
        (CRY?'<option value="vasp">VASP POSCAR (.vasp)</option>':'')+
        (CRY?'<option value="cube">Gaussian Cube (.cube)</option>':'')+
        '<option value="coord">Turbomole Coord (.coord)</option>'+
        '<option value="inp">ORCA Input (.inp)</option>'+
        '<option value="mol2">MOL2 (.mol2)</option>'+
        '<option value="mol">MDL Mol (.mol)</option>'+
        '<option value="pdb">PDB (.pdb)</option>'+
        '<option value="mop">MOPAC Input (.mop)</option>'+
        '</select>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Save</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        var fmt=document.getElementById('m-fmt').value;
        var content,ext;
        switch(fmt){
            case 'gjf':content=gjf;ext='.gjf';break;
            case 'cif':content=cifContent;ext='.cif';break;
            case 'vasp':content=vaspContent;ext='.vasp';break;
            case 'cube':content=cubeContent;ext='.cube';break;
            case 'coord':content=coord;ext='.coord';break;
            case 'inp':content=orcaInp;ext='.inp';break;
            case 'mol2':content=mol2;ext='.mol2';break;
            case 'mol':content=mol;ext='.mol';break;
            case 'pdb':content=pdb;ext='.pdb';break;
            case 'mop':content=mopac;ext='.mop';break;
            default:content=xyz;ext='.xyz';
        }
        vscodeApi.postMessage({command:'saveFile',content:content,suggestedName:'molecule_modified'+ext,filePath:MD.filePath||''});
        hideModal();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal()});
}

function showSelectAtomsModal(){
    showModal('<h3>Select Atoms</h3>'+
        '<div class="current-val">Enter indices (1-based), ranges (e.g. 3-10), or element symbols. Separate with spaces or commas.</div>'+
        '<input type="text" id="m-sel-input" placeholder="e.g. 1 3-5 C H 8" style="width:100%;padding:6px 8px;background:var(--vscode-input-background,#3c3c3c);border:1px solid var(--vscode-input-border,#3c3c3c);color:var(--vscode-input-foreground,#ccc);border-radius:3px;font-size:12px">'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Close</button><button class="mbtn mbtn-ok" id="m-ok">Select</button></div>',null);
    var inputEl=document.getElementById('m-sel-input');
    inputEl.focus();
    inputEl.addEventListener('keydown',function(e){if(e.key==='Enter'){document.getElementById('m-ok').click()}});
    document.getElementById('m-ok').addEventListener('click',function(){
        var input=inputEl.value.trim();
        if(!input){hideModal();return}
        var tokens=input.split(/[\\s,;]+/);
        tokens.forEach(function(tok){
            if(!tok)return;
            var rangeMatch=tok.match(/^(\\d+)-(\\d+)$/);
            if(rangeMatch){
                var start=parseInt(rangeMatch[1],10);
                var end=parseInt(rangeMatch[2],10);
                if(!isNaN(start)&&!isNaN(end)){
                    for(var k=start;k<=end;k++){
                        var idx=k-1;
                        if(idx>=0&&idx<MD.atoms.length&&selectedAtoms.indexOf(idx)<0)selectedAtoms.push(idx);
                    }
                }
                return;
            }
            var num=parseInt(tok,10);
            if(!isNaN(num)&&num>0){
                var idx2=num-1;
                if(idx2<MD.atoms.length&&selectedAtoms.indexOf(idx2)<0)selectedAtoms.push(idx2);
                return;
            }
            var el=tok.charAt(0).toUpperCase()+tok.slice(1).toLowerCase();
            MD.atoms.forEach(function(a,i){
                if(a.element===el&&selectedAtoms.indexOf(i)<0)selectedAtoms.push(i);
            });
        });
        highlightSelected();
        var names=selectedAtoms.map(function(i){return MD.atoms[i].element+(i+1)}).join(', ');
        selInfoEl.textContent='Selected: '+names+' ('+selectedAtoms.length+' atoms)';
        hideModal();
    });
    document.getElementById('m-cancel').addEventListener('click',function(){hideModal()});
}

function updateTransform(){
    if(diffMode){
        moleculeGroup.quaternion.copy(rotQuat);
        pivotGroup.position.set(panX,panY,0);
        if(diffMolGroup){
            diffMolGroup.quaternion.copy(diffRotQuat);
            diffPivot.position.set(diffPanX,diffPanY,0);
        }
    }else{
        moleculeGroup.quaternion.copy(rotQuat);
        pivotGroup.position.set(panX,panY,0);
    }
    if(CRY)updateAxesIndicator();
}

var canvas=renderer.domElement;
var raycaster=new THREE.Raycaster();
var mouse=new THREE.Vector2();

function getClickedAtom(e){
    var rect=canvas.getBoundingClientRect();
    if(diffMode){
        var w=rect.width;
        var h=rect.height;
        var halfW=Math.floor(w/2);
        var localX=e.clientX-rect.left;
        if(localX<halfW){
            diffActiveSide='left';
            mouse.x=(localX/halfW)*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=halfW/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            var hitsL=raycaster.intersectObjects(atomMeshes);
            if(hitsL.length>0)return hitsL[0].object.userData.index;
        }else{
            diffActiveSide='right';
            mouse.x=((localX-halfW)/(w-halfW))*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=(w-halfW)/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            var hitsR=raycaster.intersectObjects(diffAtomMeshes);
            if(hitsR.length>0)return hitsR[0].object.userData.index;
        }
        return-1;
    }
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var hits=raycaster.intersectObjects(atomMeshes);
    if(hits.length>0)return hits[0].object.userData.index;
    return-1;
}

canvas.addEventListener('mousedown',function(e){
    if(currentMode!=='view'&&e.button===0){
        var idx=getClickedAtom(e);
        if(idx>=0){selectAtom(idx);e.preventDefault();return}
    }
    if(e.button===0)isRot=true;
    else if(e.button===1||e.button===2)isPan=true;
    if(diffMode){
        var rect=canvas.getBoundingClientRect();
        var halfW=Math.floor(rect.width/2);
        diffTransformSide=(e.clientX-rect.left)<halfW?'left':'right';
    }
    prevM={x:e.clientX,y:e.clientY};
    e.preventDefault();
});

canvas.addEventListener('mousemove',function(e){
    var dm={x:e.clientX-prevM.x,y:e.clientY-prevM.y};
    if(diffMode){
        var side=diffTransformSide;
        var rq=side==='right'?diffRotQuat:rotQuat;
        if(isRot){
            var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
            var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
            rq.premultiply(qx);rq.premultiply(qy);rq.normalize();
            updateTransform()
        }
        if(isPan){
            if(side==='right'){diffPanX+=dm.x*0.01*(diffCamDist/20);diffPanY-=dm.y*0.01*(diffCamDist/20)}
            else{panX+=dm.x*0.01*(camDist/20);panY-=dm.y*0.01*(camDist/20)}
            updateTransform()
        }
    }else{
        if(isRot){
            var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
            var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
            rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();
            updateTransform()
        }
        if(isPan){panX+=dm.x*0.01*(camDist/20);panY-=dm.y*0.01*(camDist/20);updateTransform()}
    }
    prevM={x:e.clientX,y:e.clientY};
    var rect=canvas.getBoundingClientRect();
    if(diffMode){
        var w=rect.width;
        var h=rect.height;
        var halfW=Math.floor(w/2);
        var localX=e.clientX-rect.left;
        var hitsD;
        if(localX<halfW){
            mouse.x=(localX/halfW)*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=halfW/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,camDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            hitsD=raycaster.intersectObjects(atomMeshes);
        }else{
            mouse.x=((localX-halfW)/(w-halfW))*2-1;
            mouse.y=-((e.clientY-rect.top)/h)*2+1;
            camera.aspect=(w-halfW)/h;
            camera.updateProjectionMatrix();
            camera.position.set(0,0,diffCamDist);
            camera.updateMatrixWorld();
            raycaster.setFromCamera(mouse,camera);
            hitsD=raycaster.intersectObjects(diffAtomMeshes);
        }
        if(hitsD&&hitsD.length>0){
            var o=hitsD[0].object,i=o.userData.index;
            var a=diffMode&&o.userData.diffHi!==undefined&&localX>=halfW?diffData.atoms[i]:MD.atoms[i];
            tooltipEl.textContent=a.element+(i+1)+' ('+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4)+')';
            tooltipEl.style.display='block';
            tooltipEl.style.left=(e.clientX-container.getBoundingClientRect().left+15)+'px';
            tooltipEl.style.top=(e.clientY-container.getBoundingClientRect().top-10)+'px';
        }else{tooltipEl.style.display='none'}
        return;
    }
    mouse.x=((e.clientX-rect.left)/rect.width)*2-1;
    mouse.y=-((e.clientY-rect.top)/rect.height)*2+1;
    raycaster.setFromCamera(mouse,camera);
    var hits=raycaster.intersectObjects(atomMeshes);
    if(hits.length>0){var o=hits[0].object,i=o.userData.index,a=MD.atoms[i];
        tooltipEl.textContent=a.element+(i+1)+' ('+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4)+')';
        tooltipEl.style.display='block';
        tooltipEl.style.left=(e.clientX-container.getBoundingClientRect().left+15)+'px';
        tooltipEl.style.top=(e.clientY-container.getBoundingClientRect().top-10)+'px';
    }else{tooltipEl.style.display='none'}
});

canvas.addEventListener('mouseup',function(){isRot=false;isPan=false});
canvas.addEventListener('mouseleave',function(){isRot=false;isPan=false;tooltipEl.style.display='none'});
canvas.addEventListener('wheel',function(e){e.preventDefault();
    if(diffMode){
        var rect=canvas.getBoundingClientRect();
        var halfW=Math.floor(rect.width/2);
        var side=(e.clientX-rect.left)<halfW?'left':'right';
        if(side==='right'){
            diffCamDist*=e.deltaY>0?1.1:0.9;
            diffCamDist=Math.max(1,Math.min(500,diffCamDist));
        }else{
            camDist*=e.deltaY>0?1.1:0.9;
            camDist=Math.max(1,Math.min(500,camDist));
        }
    }else{
        camDist*=e.deltaY>0?1.1:0.9;
        camDist=Math.max(1,Math.min(500,camDist));
        camera.position.z=camDist;
    }
},{passive:false});
canvas.addEventListener('contextmenu',function(e){e.preventDefault()});

var touchSD=0;
canvas.addEventListener('touchstart',function(e){e.preventDefault();
    if(e.touches.length===1){
        isRot=true;prevM={x:e.touches[0].clientX,y:e.touches[0].clientY};
        if(diffMode){var rect=canvas.getBoundingClientRect();var halfW=Math.floor(rect.width/2);diffTransformSide=(e.touches[0].clientX-rect.left)<halfW?'left':'right'}
    }
    else if(e.touches.length===2){isRot=false;var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;touchSD=Math.sqrt(dx*dx+dy*dy)}
},{passive:false});
canvas.addEventListener('touchmove',function(e){e.preventDefault();
    if(e.touches.length===1&&isRot){
        var dm={x:e.touches[0].clientX-prevM.x,y:e.touches[0].clientY-prevM.y};
        var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
        var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
        if(diffMode){
            var rq=diffTransformSide==='right'?diffRotQuat:rotQuat;
            rq.premultiply(qx);rq.premultiply(qy);rq.normalize();
        }else{
            rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();
        }
        updateTransform();
        prevM={x:e.touches[0].clientX,y:e.touches[0].clientY}
    }
    else if(e.touches.length===2){var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY,d=Math.sqrt(dx*dx+dy*dy);
        if(touchSD>0){
            if(diffMode){
                var sc=touchSD/d;
                if(diffTransformSide==='right'){diffCamDist*=sc;diffCamDist=Math.max(1,Math.min(500,diffCamDist))}
                else{camDist*=sc;camDist=Math.max(1,Math.min(500,camDist))}
            }else{
                camDist*=touchSD/d;camDist=Math.max(1,Math.min(500,camDist));camera.position.z=camDist;
            }
        }
        touchSD=d}
},{passive:false});
canvas.addEventListener('touchend',function(e){isRot=false;if(e.touches.length<2)touchSD=0});

window.addEventListener('resize',function(){var rw=container.clientWidth||window.innerWidth;var rh=container.clientHeight||(window.innerHeight-60);if(rw<1)rw=window.innerWidth;if(rh<1)rh=window.innerHeight-60;camera.aspect=rw/rh;camera.updateProjectionMatrix();renderer.setSize(rw,rh)});

function animate(){
    requestAnimationFrame(animate);
    var w=container.clientWidth||window.innerWidth;
    var h=container.clientHeight||(window.innerHeight-60);
    if(w<1)w=window.innerWidth;
    if(h<1)h=window.innerHeight-60;
    if(diffMode){
        var halfW=Math.floor(w/2);
        renderer.autoClear=false;
        renderer.setScissorTest(true);
        renderer.setViewport(0,0,halfW,h);
        renderer.setScissor(0,0,halfW,h);
        camera.aspect=halfW/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,camDist);
        renderer.clear();
        if(diffPivot)diffPivot.visible=false;
        pivotGroup.visible=true;
        renderer.render(scene,camera);
        renderer.setViewport(halfW,0,w-halfW,h);
        renderer.setScissor(halfW,0,w-halfW,h);
        camera.aspect=(w-halfW)/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,diffCamDist);
        renderer.clear();
        if(diffPivot)diffPivot.visible=true;
        pivotGroup.visible=false;
        renderer.render(scene,camera);
        renderer.setScissorTest(false);
        renderer.autoClear=true;
        pivotGroup.visible=true;
        if(diffPivot)diffPivot.visible=true;
    }else{
        camera.aspect=w/h;
        camera.updateProjectionMatrix();
        camera.position.set(0,0,camDist);
        renderer.render(scene,camera);
    }
}
animate();
}catch(e){var el=document.getElementById('error-msg');var ll=document.getElementById('loading');if(ll)ll.style.display='none';if(el){el.style.display='block';el.textContent='Error: '+e.message}}
})();
</script>
</body>
</html>`;
    }
}
exports.MolecularViewerProvider = MolecularViewerProvider;
class MolecularDocument {
    constructor(uri, data, frames = []) {
        this.uri = uri;
        this.data = data;
        this.frames = frames;
    }
    dispose() { }
}
function getNonce() {
    let text = '';
    const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    for (let i = 0; i < 32; i++) {
        text += possible.charAt(Math.floor(Math.random() * possible.length));
    }
    return text;
}
//# sourceMappingURL=molecularViewer.js.map