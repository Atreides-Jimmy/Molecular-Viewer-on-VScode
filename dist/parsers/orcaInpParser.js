"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseOrcainp = parseOrcainp;
function parseOrcainp(content) {
    const lines = content.split(/\r?\n/);
    const atoms = [];
    const bonds = [];
    let title = '';
    let hasExplicitBonds = false;
    let charge = 0;
    let multiplicity = 1;
    let inXyz = false;
    let inCoordsBlock = false;
    let atomIndex = 0;
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.startsWith('#'))
            continue;
        if (line.startsWith('*')) {
            const lowerLine = line.toLowerCase();
            if (lowerLine.match(/\*\s*xyz/) || lowerLine.match(/\*\s*xyzfile/)) {
                inXyz = true;
                const parts = line.split(/\s+/);
                for (let p = 0; p < parts.length; p++) {
                    if (parts[p].toLowerCase() === 'xyz' || parts[p].toLowerCase() === 'xyzfile') {
                        if (p + 2 < parts.length) {
                            charge = parseInt(parts[p + 1], 10) || 0;
                            multiplicity = parseInt(parts[p + 2], 10) || 1;
                        }
                    }
                }
                continue;
            }
            else if (inXyz) {
                inXyz = false;
                continue;
            }
        }
        if (inXyz && line !== '') {
            const parts = line.split(/\s+/);
            if (parts.length >= 4) {
                const elem = parts[0];
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (!isNaN(x) && !isNaN(y) && !isNaN(z) && /^[A-Za-z]{1,2}$/.test(elem)) {
                    atoms.push({
                        element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                        x, y, z, index: atomIndex
                    });
                    atomIndex++;
                }
            }
        }
        if (line.toLowerCase().startsWith('%coords')) {
            inCoordsBlock = true;
            continue;
        }
        if (inCoordsBlock) {
            if (line.toLowerCase() === 'end') {
                inCoordsBlock = false;
                continue;
            }
            const lowerLine = line.toLowerCase();
            if (lowerLine.startsWith('xyz')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 3) {
                    charge = parseInt(parts[1], 10) || 0;
                    multiplicity = parseInt(parts[2], 10) || 1;
                }
                continue;
            }
            if (lowerLine.startsWith('units'))
                continue;
            if (line !== '' && !line.startsWith('%')) {
                const parts = line.split(/\s+/);
                if (parts.length >= 4) {
                    const elem = parts[0];
                    const x = parseFloat(parts[1]);
                    const y = parseFloat(parts[2]);
                    const z = parseFloat(parts[3]);
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z) && /^[A-Za-z]{1,2}$/.test(elem)) {
                        atoms.push({
                            element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                            x, y, z, index: atomIndex
                        });
                        atomIndex++;
                    }
                }
            }
        }
    }
    return { atoms, bonds, title, hasExplicitBonds, charge, multiplicity };
}
//# sourceMappingURL=orcaInpParser.js.map