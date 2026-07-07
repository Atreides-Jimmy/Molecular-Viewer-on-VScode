import { Atom, Bond, MolecularData } from '../types';

export function parseMol2(content: string): MolecularData {
    const lines = content.split(/\r?\n/);
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    let title = '';
    let hasExplicitBonds = false;

    let i = 0;
    let atomCount = 0;
    let bondCount = 0;

    while (i < lines.length) {
        const line = lines[i].trim();

        if (line.startsWith('@<TRIPOS>MOLECULE')) {
            i++;
            if (i < lines.length) {
                title = lines[i].trim();
            }
        } else if (line.startsWith('@<TRIPOS>ATOM')) {
            i++;
            let atomIndex = 0;
            while (i < lines.length && lines[i].trim() !== '' && !lines[i].trim().startsWith('@<TRIPOS>')) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 6) {
                    const elementRaw = parts[5].split('.')[0];
                    const element = elementRaw.charAt(0).toUpperCase() + elementRaw.slice(1).toLowerCase();
                    const x = parseFloat(parts[2]);
                    const y = parseFloat(parts[3]);
                    const z = parseFloat(parts[4]);

                    if (element && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ element, x, y, z, index: atomIndex });
                        atomIndex++;
                    }
                }
                i++;
            }
            atomCount = atomIndex;
            continue;
        } else if (line.startsWith('@<TRIPOS>BOND')) {
            hasExplicitBonds = true;
            i++;
            while (i < lines.length && lines[i].trim() !== '' && !lines[i].trim().startsWith('@<TRIPOS>')) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    const atom1 = parseInt(parts[1], 10) - 1;
                    const atom2 = parseInt(parts[2], 10) - 1;
                    const bondTypeStr = parts[3].toLowerCase();
                    let order = 1;
                    if (bondTypeStr === '2') {
                        order = 2;
                    } else if (bondTypeStr === 'ar') {
                        order = 1.5;
                    } else if (bondTypeStr === 'am') {
                        order = 1;
                    } else if (bondTypeStr === '3') {
                        order = 3;
                    } else if (bondTypeStr === 'du') {
                        order = 1;
                    } else {
                        const parsed = parseFloat(bondTypeStr);
                        if (!isNaN(parsed)) order = parsed;
                    }

                    if (!isNaN(atom1) && !isNaN(atom2) && atom1 >= 0 && atom2 >= 0 && atom1 !== atom2) {
                        const exists = bonds.some(b =>
                            (b.atom1 === atom1 && b.atom2 === atom2) ||
                            (b.atom1 === atom2 && b.atom2 === atom1)
                        );
                        if (!exists) {
                            bonds.push({ atom1, atom2, order });
                        }
                    }
                }
                i++;
            }
            bondCount = bonds.length;
            continue;
        }

        i++;
    }

    return { atoms, bonds, title, hasExplicitBonds };
}
