// The Import dialog: source and copy count, and the one import message the
// page receives.
import * as vscode from 'vscode';
import { ImportAtomPayload, ImportBondPayload } from './payloads';

/** 一次导入的份数上限, 与新增原子总数上限同时生效。 */
export const IMPORT_COPY_LIMIT = 50;
export const IMPORT_ATOM_LIMIT = 20000;

/** 对话框里的来源选项: 用自己的字段而不是标签文本做判定。 */
export interface ImportSourceItem extends vscode.QuickPickItem {
    source: 'file' | 'smiles';
}

/** 用户选定的来源与份数。 */
export interface ImportChoice {
    source: 'file' | 'smiles';
    copies: number;
}

/** 解析对话框输入框里的份数; 非法返回 null (空、非数字、0、负数、超上限)。 */
export function parseCopyCount(text: string): number | null {
    const trimmed = text.trim();
    if (!/^[0-9]+$/.test(trimmed)) {
        return null;
    }
    const count = parseInt(trimmed, 10);
    if (count < 1 || count > IMPORT_COPY_LIMIT) {
        return null;
    }
    return count;
}

/**
 * 弹出导入对话框: **输入框填份数, 列表里点来源** —— 先定份数, 再定结构。
 *
 * 输入框是这个对话框本来就有的那个 (快速选择总会在列表上方渲染一个输入框);
 * 两个来源选项都带 alwaysShow, 因此往里打数字不会把选项过滤掉, 它就从"没有用
 * 的筛选框"变成了"份数输入框"。
 *
 * 份数不合法时带着原输入重新弹出, 原因写在标题里 —— 不弹模态提示, 免得打断,
 * 用户也可以直接把数字改掉。返回 undefined 表示取消。
 */
export async function askImportChoice(): Promise<ImportChoice | undefined> {
    let preset = '1';
    let complaint = '';
    for (;;) {
        const items: ImportSourceItem[] = [
            { label: 'From File', description: 'pick a structure file (.xyz, .gjf, .mol2, .log, ...)', source: 'file', alwaysShow: true },
            { label: 'From SMILES', description: 'type a SMILES string; coordinates are generated locally', source: 'smiles', alwaysShow: true }
        ];
        const quickPick = vscode.window.createQuickPick<ImportSourceItem>();
        quickPick.title = complaint === '' ? 'Import structure' : 'Import structure — ' + complaint;
        quickPick.placeholder = 'Copies to import (1-' + IMPORT_COPY_LIMIT + '), then pick a source';
        quickPick.value = preset;
        quickPick.ignoreFocusOut = true;
        quickPick.matchOnDescription = false;
        quickPick.items = items;
        quickPick.activeItems = [items[0]];
        const picked = await new Promise<{ item: ImportSourceItem; text: string } | undefined>((resolve) => {
            let accepted = false;
            quickPick.onDidAccept(() => {
                accepted = true;
                const chosen = quickPick.selectedItems.length > 0
                    ? quickPick.selectedItems[0]
                    : (quickPick.activeItems.length > 0 ? quickPick.activeItems[0] : items[0]);
                resolve({ item: chosen, text: quickPick.value });
                quickPick.hide();
            });
            quickPick.onDidHide(() => {
                quickPick.dispose();
                if (!accepted) {
                    resolve(undefined);          // Esc / 被其它界面挤掉 = 取消
                }
            });
            quickPick.show();
        });
        if (picked === undefined) {
            return undefined;
        }
        const copies = parseCopyCount(picked.text);
        if (copies !== null) {
            return { source: picked.item.source, copies: copies };
        }
        preset = picked.text;                    // 保留原输入, 直接改数字即可
        complaint = 'copies must be a whole number from 1 to ' + IMPORT_COPY_LIMIT;
    }
}

/** 相邻两份之间的间隔参数: 与客户端摆放使用的缓冲/间隙保持一致。 */
export const IMPORT_PLACEMENT_BUFFER = 2.5;
export const IMPORT_PLACEMENT_GAP = 2.0;

/**
 * 把 N 份拷贝拼成一个整体。份与份之间用**同一个**间距 (
 * 两份各自的包围球半径 + 缓冲 + 间隙), 以共同质心为中心均匀排开 —— 于是每两份
 * 之间的最小距离都与"单独导入两次"所保证的一样, 不会重叠。
 *
 * 单份时原样返回 (不复制数组), 保证"只导入一份"的行为与以前逐位相同。
 */
export function buildImportBlock(atoms: ImportAtomPayload[], bonds: ImportBondPayload[],
                          copies: number): { atoms: ImportAtomPayload[]; bonds: ImportBondPayload[] } {
    if (copies <= 1 || atoms.length === 0) {
        return { atoms: atoms, bonds: bonds };
    }
    // 单份的包围球 (与客户端摆放判定用的是同一个量)
    let cx = 0.0;
    let cy = 0.0;
    let cz = 0.0;
    for (let i = 0; i < atoms.length; i++) {
        cx += atoms[i].x;
        cy += atoms[i].y;
        cz += atoms[i].z;
    }
    cx /= atoms.length;
    cy /= atoms.length;
    cz /= atoms.length;
    let radius = 0.0;
    for (let i = 0; i < atoms.length; i++) {
        const dx = atoms[i].x - cx;
        const dy = atoms[i].y - cy;
        const dz = atoms[i].z - cz;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > radius) {
            radius = d;
        }
    }
    const stride = 2.0 * (radius + IMPORT_PLACEMENT_BUFFER) + IMPORT_PLACEMENT_GAP;
    const outAtoms: ImportAtomPayload[] = [];
    const outBonds: ImportBondPayload[] = [];
    for (let c = 0; c < copies; c++) {
        const shift = (c - (copies - 1) / 2.0) * stride;
        const base = c * atoms.length;
        for (let i = 0; i < atoms.length; i++) {
            const a = atoms[i];
            outAtoms.push({ element: a.element, x: a.x + shift, y: a.y, z: a.z, color: a.color });
        }
        for (let k = 0; k < bonds.length; k++) {
            const b = bonds[k];
            outBonds.push({ atom1: base + b.atom1, atom2: base + b.atom2, order: b.order });
        }
    }
    return { atoms: outAtoms, bonds: outBonds };
}

/**
 * 按份数拼好后**一次性**发出 (只发一条导入消息)。
 *
 * 位置计算因此只发生一次: 客户端对"整个块"做一次包围球避让摆放 (第 1 份的位置
 * 由这次计算确定), 份与份之间则是上面那个等间距, 客户端还会一次成键、一次场景
 * 重建, 撤销也只要一步 —— 逐份重算位置、逐份重建场景的开销全部省掉。
 *
 * 份数已在对话框里校验; 这里再按"新增原子总数"上限截断 (此时才知道结构有多大)
 * 并告知用户, 免得一次导入把视图撑到无法交互。
 */
export function postImportCopies(webview: vscode.Webview, fileName: string,
                          atoms: ImportAtomPayload[], bonds: ImportBondPayload[],
                          hasExplicitBonds: boolean, copies: number): void {
    const affordable = Math.max(1, Math.floor(IMPORT_ATOM_LIMIT / Math.max(1, atoms.length)));
    let count = copies;
    if (count > affordable) {
        count = affordable;
        vscode.window.showWarningMessage('Import: ' + copies + ' copies of a ' + atoms.length +
            '-atom structure would add ' + (copies * atoms.length) + ' atoms, so this import is limited to ' +
            count + ' copies.');
    }
    const block = buildImportBlock(atoms, bonds, count);
    webview.postMessage({
        command: 'importResult',
        cancelled: false,
        fileName: count > 1 ? fileName + ' x ' + count : fileName,
        atoms: block.atoms,
        bonds: block.bonds,
        hasExplicitBonds: hasExplicitBonds
    });
}
