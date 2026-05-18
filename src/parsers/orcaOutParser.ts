import { Atom, Bond, MolecularData } from '../types';

export interface OrcaFrame {
    atoms: Atom[];
    bonds: Bond[];
    title: string;
    hasExplicitBonds: boolean;
    stepLabel: string;
    charge?: number;
    multiplicity?: number;
}

export function parseOrcaOut(content: string): { frames: OrcaFrame[], title: string, charge?: number, multiplicity?: number } {
    const lines = content.split(/\r?\n/);
    const frames: OrcaFrame[] = [];
    let title = '';
    let charge: number | undefined;
    let multiplicity: number | undefined;
    let stepCount = 0;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        if (line.includes('Total Charge') && line.includes('Charge')) {
            const m = line.match(/Charge\s*\.\.\.\s*(-?\d+)/i);
            if (m) charge = parseInt(m[1], 10);
        }

        if (line.includes('Multiplicity') && line.includes('Mult')) {
            const m = line.match(/Mult\s*\.\.\.\s*(\d+)/i);
            if (m) multiplicity = parseInt(m[1], 10);
        }

        if (line.includes('CARTESIAN COORDINATES (ANGSTROEM)')) {
            const atoms: Atom[] = [];
            let j = i + 2;
            while (j < lines.length) {
                const tl = lines[j].trim();
                if (tl === '' || tl.startsWith('-')) break;
                const parts = tl.split(/\s+/);
                if (parts.length >= 4) {
                    const elem = parts[0];
                    const x = parseFloat(parts[1]);
                    const y = parseFloat(parts[2]);
                    const z = parseFloat(parts[3]);
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z) && /^[A-Za-z]{1,2}$/.test(elem)) {
                        atoms.push({
                            element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                            x, y, z, index: atoms.length
                        });
                    }
                }
                j++;
            }
            if (atoms.length > 0) {
                stepCount++;
                frames.push({
                    atoms,
                    bonds: [],
                    title,
                    hasExplicitBonds: false,
                    stepLabel: 'Step ' + stepCount,
                    charge,
                    multiplicity
                });
            }
        }
    }

    if (frames.length === 0) {
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line.includes('CARTESIAN COORDINATES (A.U.)')) {
                const atoms: Atom[] = [];
                let j = i + 2;
                while (j < lines.length) {
                    const tl = lines[j].trim();
                    if (tl === '' || tl.startsWith('-')) break;
                    const parts = tl.split(/\s+/);
                    if (parts.length >= 7) {
                        const elem = parts[1];
                        const x = parseFloat(parts[5]) * 0.529177249;
                        const y = parseFloat(parts[6]) * 0.529177249;
                        const z = parseFloat(parts[7]) * 0.529177249;
                        if (!isNaN(x) && !isNaN(y) && !isNaN(z) && /^[A-Za-z]{1,2}$/.test(elem)) {
                            atoms.push({
                                element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                                x, y, z, index: atoms.length
                            });
                        }
                    }
                    j++;
                }
                if (atoms.length > 0) {
                    stepCount++;
                    frames.push({
                        atoms,
                        bonds: [],
                        title,
                        hasExplicitBonds: false,
                        stepLabel: 'Step ' + stepCount,
                        charge,
                        multiplicity
                    });
                }
            }
        }
    }

    return { frames, title, charge, multiplicity };
}
