import { AtomGroup } from '../types';

const VMD_COLORS: { [id: number]: string } = {
    0: '#0000FF',
    1: '#FF0000',
    2: '#808080',
    3: '#FF6600',
    4: '#FFFF00',
    5: '#D2B48C',
    6: '#C0C0C0',
    7: '#00FF00',
    8: '#FFFFFF',
    9: '#FF69B4',
    10: '#00FFFF',
    11: '#800080',
    12: '#00FF00',
    13: '#E0B0FF',
    14: '#CC7722',
    15: '#5F9EA0',
    16: '#000000',
    17: '#CCCC00',
    18: '#999900',
    19: '#00CC00',
    20: '#009900',
    21: '#00CCCC',
    22: '#009999',
    23: '#0000CC',
    24: '#000099',
    25: '#CC00CC',
    26: '#990099',
    27: '#FF00FF',
    28: '#CC0099',
    29: '#CC0000',
    30: '#990000',
    31: '#CC6600',
    32: '#996600'
};

export function getVmdColor(colorId: number): string {
    return VMD_COLORS[colorId] || '#FF1493';
}

export interface TclParseResult {
    sourceFile: string;
    sourceType: string;
    groups: AtomGroup[];
}

function parseSelectionIndices(selection: string): number[] {
    const indices: number[] = [];
    const trimmed = selection.trim();

    if (trimmed === 'all') {
        return [];
    }

    const indexMatch = trimmed.match(/index\s+([\d\s]+)/);
    if (indexMatch) {
        const nums = indexMatch[1].trim().split(/\s+/);
        for (const n of nums) {
            const val = parseInt(n, 10);
            if (!isNaN(val)) {
                indices.push(val);
            }
        }
        return indices;
    }

    const indexRangeMatch = trimmed.match(/index\s+(\d+)\s+to\s+(\d+)/);
    if (indexRangeMatch) {
        const start = parseInt(indexRangeMatch[1], 10);
        const end = parseInt(indexRangeMatch[2], 10);
        if (!isNaN(start) && !isNaN(end)) {
            for (let i = start; i <= end; i++) {
                indices.push(i);
            }
        }
        return indices;
    }

    return indices;
}

export function parseTcl(content: string): TclParseResult {
    const lines = content.split(/\r?\n/);
    let sourceFile = '';
    let sourceType = '';
    const groups: AtomGroup[] = [];

    let currentColorId: number | null = null;
    let currentSelection = '';

    for (const line of lines) {
        const trimmed = line.trim();

        const molNewMatch = trimmed.match(/^mol\s+new\s+(\S+)/);
        if (molNewMatch) {
            sourceFile = molNewMatch[1];
        }

        const typeMatch = trimmed.match(/^mol\s+new\s+\S+\s+type\s+(\S+)/);
        if (typeMatch) {
            sourceType = typeMatch[1];
        }

        const colorMatch = trimmed.match(/^mol\s+color\s+ColorID\s+(\d+)/);
        if (colorMatch) {
            currentColorId = parseInt(colorMatch[1], 10);
        }

        const colorBetaMatch = trimmed.match(/^mol\s+color\s+Beta/);
        if (colorBetaMatch) {
            currentColorId = null;
        }

        const selMatch = trimmed.match(/^mol\s+selection\s+\{(.+?)\}/);
        if (!selMatch) {
            const selMatch2 = trimmed.match(/^mol\s+selection\s+"(.+?)"/);
            if (selMatch2) {
                currentSelection = selMatch2[1];
            }
        } else {
            currentSelection = selMatch[1];
        }

        if (trimmed === 'mol addrep top') {
            if (currentColorId !== null && currentSelection) {
                const indices = parseSelectionIndices(currentSelection);
                if (indices.length > 0) {
                    groups.push({
                        colorId: currentColorId,
                        color: getVmdColor(currentColorId),
                        indices: indices
                    });
                }
            }
            currentColorId = null;
            currentSelection = '';
        }
    }

    return { sourceFile, sourceType, groups };
}
