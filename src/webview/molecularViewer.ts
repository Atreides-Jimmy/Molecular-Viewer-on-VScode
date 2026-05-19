import * as vscode from 'vscode';
import { parseFile, parseLogFile, LogFrame, OrcaFrame } from '../parsers/index';
import { ensureBonds } from '../parsers/bondDetector';
import { MolecularData } from '../types';

export class MolecularViewerProvider implements vscode.CustomReadonlyEditorProvider<MolecularDocument> {
    constructor(private readonly context: vscode.ExtensionContext) {}

    async openCustomDocument(
        uri: vscode.Uri,
        _openContext: vscode.CustomDocumentOpenContext,
        _token: vscode.CancellationToken
    ): Promise<MolecularDocument> {
        const content = await vscode.workspace.fs.readFile(uri);
        const textContent = new TextDecoder().decode(content);
        const fileName = uri.path.split('/').pop() || 'unknown.xyz';

        const ext = fileName.toLowerCase().split('.').pop() || '';
        let data: MolecularData;
        let frames: (LogFrame | OrcaFrame)[] = [];

        if (ext === 'log' || ext === 'out') {
            const logResult = parseLogFile(textContent, fileName);
            frames = logResult.frames;
            if (frames.length > 0) {
                data = ensureBonds({
                    atoms: frames[0].atoms,
                    bonds: frames[0].bonds,
                    title: frames[0].title,
                    hasExplicitBonds: frames[0].hasExplicitBonds,
                    charge: frames[0].charge,
                    multiplicity: frames[0].multiplicity
                });
            } else {
                data = { atoms: [], bonds: [], title: 'No structures found', hasExplicitBonds: false };
            }
        } else {
            data = parseFile(textContent, fileName);
            data = ensureBonds(data);
        }

        data.filePath = uri.fsPath;

        return new MolecularDocument(uri, data, frames);
    }

    async resolveCustomEditor(
        document: MolecularDocument,
        webviewPanel: vscode.WebviewPanel,
        _token: vscode.CancellationToken
    ): Promise<void> {
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
                                'All Files': ['*']
                            }
                        });
                        if (uri) {
                            await vscode.workspace.fs.writeFile(uri, new TextEncoder().encode(message.content));
                            vscode.window.showInformationMessage('Saved: ' + uri.fsPath);
                        }
                    } catch (e: any) {
                        vscode.window.showErrorMessage('Save failed: ' + (e.message || e));
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

    private async getHtmlForWebview(webview: vscode.Webview, data: MolecularData, frames: (LogFrame | OrcaFrame)[] = []): Promise<string> {
        const nonce = getNonce();

        const threeJsBytes = await vscode.workspace.fs.readFile(
            vscode.Uri.joinPath(this.context.extensionUri, 'media', 'three.min.js')
        );
        const threeJsContent = new TextDecoder().decode(threeJsBytes);

        const atomColors: { [key: string]: string } = {
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

        const atomData = data.atoms.map(a => ({
            element: a.element, x: a.x, y: a.y, z: a.z,
            color: atomColors[a.element] || '#FF1493'
        }));

        const bondData = data.bonds.map(b => ({
            atom1: b.atom1, atom2: b.atom2, order: b.order
        }));

        const framesData = frames.map(f => ({
            atoms: f.atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z, color: atomColors[a.element] || '#FF1493' })),
            bonds: f.bonds.map(b => ({ atom1: b.atom1, atom2: b.atom2, order: b.order })),
            stepLabel: f.stepLabel
        }));

        const jsonData = JSON.stringify({ atoms: atomData, bonds: bondData, title: data.title, atomColors: atomColors, filePath: data.filePath || '', frames: framesData, gjfMeta: data.gjfMeta || null, charge: data.charge, multiplicity: data.multiplicity });

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
<div id="container"><div id="loading">Loading 3D Viewer...</div><div id="mol-info"></div></div>
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
function updateMolInfo(){
    var infoEl=document.getElementById('mol-info');
    if(!infoEl)return;
    var nAtoms=MD.atoms.length;
    var chrg=MD.charge!=null?MD.charge:'-';
    var mult=MD.multiplicity!=null?MD.multiplicity:'-';
    var nElectrons=0;
    MD.atoms.forEach(function(a){nElectrons+=(AN[a.element]||0)});
    if(typeof chrg==='number')nElectrons-=chrg;
    infoEl.innerHTML='Atoms: '+nAtoms+'<br>Charge: '+chrg+'<br>Electrons: '+nElectrons+'<br>Multiplicity: '+mult;
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

var atomMeshes=[];
var bondMeshes=[];

function rebuildScene(){
    while(moleculeGroup.children.length>0)moleculeGroup.remove(moleculeGroup.children[0]);
    atomMeshes.length=0;
    bondMeshes.length=0;
    CX=0;CY=0;CZ=0;
    MD.atoms.forEach(function(a){CX+=a.x;CY+=a.y;CZ+=a.z});
    if(MD.atoms.length>0){CX/=MD.atoms.length;CY/=MD.atoms.length;CZ/=MD.atoms.length}
    MD.atoms.forEach(function(a,i){
        a.index=i;
        var r=getR(a.element);
        var g=new THREE.SphereGeometry(r,32,24);
        var m=new THREE.MeshPhongMaterial({color:new THREE.Color(a.color),shininess:80,specular:0x444444});
        var mesh=new THREE.Mesh(g,m);
        mesh.position.set(a.x-CX,a.y-CY,a.z-CZ);
        mesh.userData={element:a.element,index:i};
        moleculeGroup.add(mesh);
        atomMeshes.push(mesh);
    });
    MD.bonds.forEach(function(b){createBond(b)});
    highlightSelected();
}

function getPerp(dir){
    var up=Math.abs(dir.y)<0.99?new THREE.Vector3(0,1,0):new THREE.Vector3(1,0,0);
    return new THREE.Vector3().crossVectors(dir,up).normalize();
}

function createBond(b){
    var a1=MD.atoms[b.atom1],a2=MD.atoms[b.atom2];
    if(!a1||!a2)return;
    var s=new THREE.Vector3(a1.x-CX,a1.y-CY,a1.z-CZ);
    var e=new THREE.Vector3(a2.x-CX,a2.y-CY,a2.z-CZ);
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

rebuildScene();

var maxD=0;
MD.atoms.forEach(function(a){var dx=a.x-CX,dy=a.y-CY,dz=a.z-CZ,dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd>maxD)maxD=dd});
var initCam=maxD*2.5+5;
camera.position.set(0,0,initCam);camera.lookAt(0,0,0);
camDist=initCam;

var MODE_INFO={view:'View Mode',bondLength:'Bond Length - Click 2 atoms',bondAngle:'Bond Angle - Click 3 atoms (central 2nd)',dihedral:'Dihedral - Click 4 atoms',addAtom:'Add Atom - Click anchor atom',deleteAtom:'Delete Atom - Click atom to delete',selectAtoms:'Select Atoms - Input indices or element symbols'};

function setMode(m){
    if(currentMode===m)return;
    currentMode=m;selectedAtoms=[];originalCoords=null;hideModal();highlightSelected();
    modeInfoEl.textContent=MODE_INFO[m]||m;
    selInfoEl.textContent='';
    document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.classList.toggle('active',b.dataset.mode===m)});
    if(m==='selectAtoms')showSelectAtomsModal();
}

function resetSelection(){
    selectedAtoms=[];originalCoords=null;highlightSelected();
    selInfoEl.textContent='';
}

document.querySelectorAll('.tbtn[data-mode]').forEach(function(b){b.addEventListener('click',function(){setMode(this.dataset.mode)})});
document.getElementById('reset-btn').addEventListener('click',function(){rotQuat.identity();panX=0;panY=0;camDist=initCam;camera.position.set(0,0,camDist);updateTransform()});
document.getElementById('save-btn').addEventListener('click',doSave);

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
    var CR={H:0.31,He:0.28,Li:1.28,Be:0.96,B:0.84,C:0.76,N:0.71,O:0.66,F:0.57,Na:1.66,Mg:1.41,Al:1.21,Si:1.11,P:1.07,S:1.05,Cl:1.02,K:2.03,Ca:1.76,Fe:1.32,Cu:1.32,Zn:1.22,Br:1.20,I:1.39};
    var tol=0.45;
    var bonds=[];
    for(var i=0;i<atoms.length;i++){
        for(var j=i+1;j<atoms.length;j++){
            var dx=atoms[i].x-atoms[j].x,dy=atoms[i].y-atoms[j].y,dz=atoms[i].z-atoms[j].z;
            var d=Math.sqrt(dx*dx+dy*dy+dz*dz);
            var r1=CR[atoms[i].element]||1.5,r2=CR[atoms[j].element]||1.5;
            var maxD=r1+r2+tol;
            if(d<=maxD){
                var ratio=d/(r1+r2);
                var order=1;
                if(ratio<=0.78)order=3;
                else if(ratio<=0.88)order=2;
                bonds.push({atom1:i,atom2:j,order:order});
            }
        }
    }
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
        else{m.material.emissive=new THREE.Color(0x000000);m.material.emissiveIntensity=0}
    });
}

function selectAtom(idx){
    if(selectedAtoms.indexOf(idx)>=0)return;
    selectedAtoms.push(idx);
    highlightSelected();
    var names=selectedAtoms.map(function(i){return MD.atoms[i].element+(i+1)}).join(', ');
    selInfoEl.textContent='Selected: '+names;
    checkSelectionComplete();
}

function checkSelectionComplete(){
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
    rebuildScene();
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
    rebuildScene();
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
    rebuildScene();
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
        '<input type="range" id="m-slider" value="'+cur.toFixed(2)+'" min="-180" max="180" step="1">'+
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
        MD.atoms.splice(idx,1);
        MD.atoms.forEach(function(a,i){a.index=i});
        MD.bonds=MD.bonds.filter(function(b){return b.atom1!==idx&&b.atom2!==idx}).map(function(b){
            return{atom1:b.atom1>idx?b.atom1-1:b.atom1,atom2:b.atom2>idx?b.atom2-1:b.atom2,order:b.order};
        });
        rebuildScene();hideModal();resetSelection();
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
    MD.atoms.forEach(function(a){gjf+=' '+a.element+'   '+a.x.toFixed(6)+'   '+a.y.toFixed(6)+'   '+a.z.toFixed(6)+'\\n'});
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

    showModal('<h3>Save File</h3>'+
        '<label>Format:</label><select id="m-fmt">'+
        '<option value="xyz">XYZ (.xyz)</option>'+
        '<option value="gjf">Gaussian Input (.gjf)</option>'+
        '<option value="coord">Turbomole Coord (.coord)</option>'+
        '<option value="inp">ORCA Input (.inp)</option>'+
        '<option value="mol2">MOL2 (.mol2)</option>'+
        '<option value="mol">MDL Mol (.mol)</option>'+
        '</select>'+
        '<div class="modal-btns"><button class="mbtn mbtn-cancel" id="m-cancel">Cancel</button><button class="mbtn mbtn-ok" id="m-ok">Save</button></div>',null);
    document.getElementById('m-ok').addEventListener('click',function(){
        var fmt=document.getElementById('m-fmt').value;
        var content,ext;
        switch(fmt){
            case 'gjf':content=gjf;ext='.gjf';break;
            case 'coord':content=coord;ext='.coord';break;
            case 'inp':content=orcaInp;ext='.inp';break;
            case 'mol2':content=mol2;ext='.mol2';break;
            case 'mol':content=mol;ext='.mol';break;
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
        selectedAtoms=[];
        var tokens=input.split(/[\s,;]+/);
        tokens.forEach(function(tok){
            if(!tok)return;
            var rangeMatch=tok.match(/^(\d+)-(\d+)$/);
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
    moleculeGroup.quaternion.copy(rotQuat);
    pivotGroup.position.set(panX,panY,0);
}

var canvas=renderer.domElement;
var raycaster=new THREE.Raycaster();
var mouse=new THREE.Vector2();

function getClickedAtom(e){
    var rect=canvas.getBoundingClientRect();
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
    prevM={x:e.clientX,y:e.clientY};
    e.preventDefault();
});

canvas.addEventListener('mousemove',function(e){
    var dm={x:e.clientX-prevM.x,y:e.clientY-prevM.y};
    if(isRot){
        var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);
        var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);
        rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();
        updateTransform()
    }
    if(isPan){panX+=dm.x*0.01*(camDist/20);panY-=dm.y*0.01*(camDist/20);updateTransform()}
    prevM={x:e.clientX,y:e.clientY};
    var rect=canvas.getBoundingClientRect();
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
canvas.addEventListener('wheel',function(e){e.preventDefault();camDist*=e.deltaY>0?1.1:0.9;camDist=Math.max(1,Math.min(500,camDist));camera.position.z=camDist},{passive:false});
canvas.addEventListener('contextmenu',function(e){e.preventDefault()});

var touchSD=0;
canvas.addEventListener('touchstart',function(e){e.preventDefault();
    if(e.touches.length===1){isRot=true;prevM={x:e.touches[0].clientX,y:e.touches[0].clientY}}
    else if(e.touches.length===2){isRot=false;var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;touchSD=Math.sqrt(dx*dx+dy*dy)}
},{passive:false});
canvas.addEventListener('touchmove',function(e){e.preventDefault();
    if(e.touches.length===1&&isRot){var dm={x:e.touches[0].clientX-prevM.x,y:e.touches[0].clientY-prevM.y};var qx=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0,1,0),dm.x*0.008);var qy=new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1,0,0),dm.y*0.008);rotQuat.premultiply(qx);rotQuat.premultiply(qy);rotQuat.normalize();updateTransform();prevM={x:e.touches[0].clientX,y:e.touches[0].clientY}}
    else if(e.touches.length===2){var dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY,d=Math.sqrt(dx*dx+dy*dy);
        if(touchSD>0){camDist*=touchSD/d;camDist=Math.max(1,Math.min(500,camDist));camera.position.z=camDist}touchSD=d}
},{passive:false});
canvas.addEventListener('touchend',function(e){isRot=false;if(e.touches.length<2)touchSD=0});

window.addEventListener('resize',function(){var rw=container.clientWidth||window.innerWidth;var rh=container.clientHeight||(window.innerHeight-60);if(rw<1)rw=window.innerWidth;if(rh<1)rh=window.innerHeight-60;camera.aspect=rw/rh;camera.updateProjectionMatrix();renderer.setSize(rw,rh)});

function animate(){requestAnimationFrame(animate);renderer.render(scene,camera)}
animate();
}catch(e){var el=document.getElementById('error-msg');var ll=document.getElementById('loading');if(ll)ll.style.display='none';if(el){el.style.display='block';el.textContent='Error: '+e.message}}
})();
</script>
</body>
</html>`;
    }
}

class MolecularDocument implements vscode.CustomDocument {
    constructor(
        public readonly uri: vscode.Uri,
        public readonly data: MolecularData,
        public readonly frames: (LogFrame | OrcaFrame)[] = []
    ) {}

    dispose(): void {}
}

function getNonce(): string {
    let text = '';
    const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    for (let i = 0; i < 32; i++) {
        text += possible.charAt(Math.floor(Math.random() * possible.length));
    }
    return text;
}
