import * as vscode from 'vscode';
import { MolecularViewerProvider } from './webview/molecularViewer';
import { parseFile } from './parsers/index';
import { ensureBonds } from './parsers/bondDetector';
import { MolecularData } from './types';

export function activate(context: vscode.ExtensionContext) {
    const provider = new MolecularViewerProvider(context);

    context.subscriptions.push(
        vscode.window.registerCustomEditorProvider(
            'molecularViewer.editor',
            provider,
            {
                webviewOptions: {
                    retainContextWhenHidden: true,
                },
                supportsMultipleEditorsPerDocument: false,
            }
        )
    );

    context.subscriptions.push(
        vscode.commands.registerCommand('molecularViewer.openViewer', async (uri?: vscode.Uri) => {
            let fileUri: vscode.Uri;

            if (uri) {
                fileUri = uri;
            } else {
                const activeEditor = vscode.window.activeTextEditor;
                if (activeEditor) {
                    fileUri = activeEditor.document.uri;
                } else {
                    const result = await vscode.window.showOpenDialog({
                        canSelectMany: false,
                        openLabel: 'Select Molecular File',
                        filters: {
                            'Molecular Files': ['gjf', 'xyz', 'mol', 'sdf', 'gjf03', 'gjf09', 'gjf16', 'com', 'mol2', 'log', 'out', 'coord', 'inp', 'pdb', 'ent', 'mop', 'mopac', 'dat', 'tcl'],
                            'All Files': ['*']
                        }
                    });
                    if (!result || result.length === 0) {
                        return;
                    }
                    fileUri = result[0];
                }
            }

            await vscode.commands.executeCommand(
                'vscode.openWith',
                fileUri,
                'molecularViewer.editor'
            );
        })
    );
}

export function deactivate() {}
