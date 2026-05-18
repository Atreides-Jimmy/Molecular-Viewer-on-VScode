"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseCoord = parseCoord;
const BOHR_TO_ANG = 0.529177249;
function parseCoord(content) {
    const lines = content.split(/\r?\n/);
    const atoms = [];
    const bonds = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge = 0;
    let multiplicity = 1;
    let inCoord = false;
    let atomIndex = 0;
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.toLowerCase() === '$coord' || line.toLowerCase() === '$coord ') {
            inCoord = true;
            continue;
        }
        if (line.startsWith('$') && inCoord) {
            inCoord = false;
        }
        if (inCoord && line !== '' && !line.startsWith('#')) {
            const parts = line.split(/\s+/);
            if (parts.length >= 4) {
                const x = parseFloat(parts[0]);
                const y = parseFloat(parts[1]);
                const z = parseFloat(parts[2]);
                const elem = parts[3];
                if (!isNaN(x) && !isNaN(y) && !isNaN(z) && /^[A-Za-z]{1,2}$/.test(elem)) {
                    atoms.push({
                        element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                        x: x * BOHR_TO_ANG,
                        y: y * BOHR_TO_ANG,
                        z: z * BOHR_TO_ANG,
                        index: atomIndex
                    });
                    atomIndex++;
                }
            }
        }
        if (line.toLowerCase().startsWith('$chrg')) {
            const chrgVal = parseInt(line.split(/\s+/)[1], 10);
            if (!isNaN(chrgVal))
                charge = chrgVal;
        }
        if (line.toLowerCase().startsWith('$spin') || line.toLowerCase().startsWith('$mult')) {
            const spinVal = parseInt(line.split(/\s+/)[1], 10);
            if (!isNaN(spinVal))
                multiplicity = spinVal;
        }
    }
    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
//# sourceMappingURL=coordParser.js.map