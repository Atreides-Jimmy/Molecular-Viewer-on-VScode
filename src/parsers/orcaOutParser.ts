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
    // Single pass: charge/multiplicity, Angstrom geometry blocks and A.U.
    // geometry blocks are all detected in ONE traversal (the A.U. blocks used
    // to be a second full scan). Consumed atom rows never match the other
    // detectors, so the merged scan yields exactly the previous results; the
    // A.U. frames are used only when no Angstrom frame exists.
    const angstromFrames: OrcaFrame[] = [];
    const auFrames: OrcaFrame[] = [];
    let title = '';
    let charge: number | undefined;
    let multiplicity: number | undefined;

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
                angstromFrames.push({
                    atoms,
                    bonds: [],
                    title,
                    hasExplicitBonds: false,
                    stepLabel: 'Step ' + (angstromFrames.length + 1),
                    charge,
                    multiplicity
                });
            }
            i = j - 1; // resume after the consumed atom block
            continue;
        }

        if (line.includes('CARTESIAN COORDINATES (A.U.)')) {
            const atoms: Atom[] = [];
            let j = i + 2;
            while (j < lines.length) {
                const tl = lines[j].trim();
                if (tl === '' || tl.startsWith('-')) break;
                const parts = tl.split(/\s+/);
                // ORCA A.U. format: "index elem nuclear_charge x y z" (6 cols)
                // but some variants have extra columns. Find the element symbol
                // (first non-numeric token) and take the last 3 numbers as coords.
                let elem = '';
                const nums: number[] = [];
                for (const p of parts) {
                    if (/^[A-Za-z]{1,2}$/.test(p) && !elem) {
                        elem = p;
                    } else {
                        const n = parseFloat(p);
                        if (!isNaN(n)) nums.push(n);
                    }
                }
                if (elem && nums.length >= 3) {
                    const x = nums[nums.length - 3] * 0.529177249;
                    const y = nums[nums.length - 2] * 0.529177249;
                    const z = nums[nums.length - 1] * 0.529177249;
                    atoms.push({
                        element: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                        x, y, z, index: atoms.length
                    });
                }
                j++;
            }
            if (atoms.length > 0) {
                auFrames.push({
                    atoms,
                    bonds: [],
                    title,
                    hasExplicitBonds: false,
                    stepLabel: 'Step ' + (auFrames.length + 1),
                    charge,
                    multiplicity
                });
            }
            i = j - 1; // resume after the consumed atom block
            continue;
        }
    }

    const frames = angstromFrames.length > 0 ? angstromFrames : auFrames;
    return { frames, title, charge, multiplicity };
}
