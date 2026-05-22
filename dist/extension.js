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
exports.activate = activate;
exports.deactivate = deactivate;
const vscode = __importStar(require("vscode"));
const molecularViewer_1 = require("./webview/molecularViewer");
function activate(context) {
    const provider = new molecularViewer_1.MolecularViewerProvider(context);
    context.subscriptions.push(vscode.window.registerCustomEditorProvider('molecularViewer.editor', provider, {
        webviewOptions: {
            retainContextWhenHidden: true,
        },
        supportsMultipleEditorsPerDocument: false,
    }));
    context.subscriptions.push(vscode.commands.registerCommand('molecularViewer.openViewer', async (uri) => {
        let fileUri;
        if (uri) {
            fileUri = uri;
        }
        else {
            const activeEditor = vscode.window.activeTextEditor;
            if (activeEditor) {
                fileUri = activeEditor.document.uri;
            }
            else {
                const result = await vscode.window.showOpenDialog({
                    canSelectMany: false,
                    openLabel: 'Select Molecular File',
                    filters: {
                        'Molecular Files': ['gjf', 'xyz', 'mol', 'sdf', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol2', 'log', 'out', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat'],
                        'All Files': ['*']
                    }
                });
                if (!result || result.length === 0) {
                    return;
                }
                fileUri = result[0];
            }
        }
        await vscode.commands.executeCommand('vscode.openWith', fileUri, 'molecularViewer.editor');
    }));
}
function deactivate() { }
//# sourceMappingURL=extension.js.map