"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseVesta = exports.parseCube = exports.parseVasp = exports.parseCif = exports.parseTcl = exports.parseOrcaOut = exports.parseGaussianLog = void 0;
exports.parseFile = parseFile;
exports.parseLogFile = parseLogFile;
const gjfParser_1 = require("./gjfParser");
const xyzParser_1 = require("./xyzParser");
const mol2Parser_1 = require("./mol2Parser");
const logParser_1 = require("./logParser");
Object.defineProperty(exports, "parseGaussianLog", { enumerable: true, get: function () { return logParser_1.parseGaussianLog; } });
const coordParser_1 = require("./coordParser");
const orcaInpParser_1 = require("./orcaInpParser");
const orcaOutParser_1 = require("./orcaOutParser");
Object.defineProperty(exports, "parseOrcaOut", { enumerable: true, get: function () { return orcaOutParser_1.parseOrcaOut; } });
const pdbParser_1 = require("./pdbParser");
const mopacParser_1 = require("./mopacParser");
const tclParser_1 = require("./tclParser");
Object.defineProperty(exports, "parseTcl", { enumerable: true, get: function () { return tclParser_1.parseTcl; } });
const cifParser_1 = require("./cifParser");
Object.defineProperty(exports, "parseCif", { enumerable: true, get: function () { return cifParser_1.parseCif; } });
const vaspParser_1 = require("./vaspParser");
Object.defineProperty(exports, "parseVasp", { enumerable: true, get: function () { return vaspParser_1.parseVasp; } });
const cubeParser_1 = require("./cubeParser");
Object.defineProperty(exports, "parseCube", { enumerable: true, get: function () { return cubeParser_1.parseCube; } });
const vestaParser_1 = require("./vestaParser");
Object.defineProperty(exports, "parseVesta", { enumerable: true, get: function () { return vestaParser_1.parseVesta; } });
function parseFile(content, fileName) {
    const ext = fileName.toLowerCase().split('.').pop() || '';
    switch (ext) {
        case 'gjf':
        case 'gjf03':
        case 'gjf09':
        case 'gjf16':
        case 'com':
            return (0, gjfParser_1.parseGjf)(content);
        case 'xyz':
            return (0, xyzParser_1.parseXyz)(content);
        case 'mol2':
            return (0, mol2Parser_1.parseMol2)(content);
        case 'log':
        case 'out':
            return parseLogAsSingleFrame(content);
        case 'coord':
            return (0, coordParser_1.parseCoord)(content);
        case 'inp':
            return (0, orcaInpParser_1.parseOrcainp)(content);
        case 'pdb':
        case 'ent':
            return (0, pdbParser_1.parsePdb)(content);
        case 'mop':
        case 'mopac':
        case 'dat':
            return (0, mopacParser_1.parseMopac)(content);
        case 'cif':
            return (0, cifParser_1.parseCif)(content);
        case 'vasp':
        case 'poscar':
        case 'contcar':
            return (0, vaspParser_1.parseVasp)(content);
        case 'cube':
            return (0, cubeParser_1.parseCube)(content);
        case 'vesta':
            return (0, vestaParser_1.parseVesta)(content);
        default:
            return tryAutoParse(content, fileName);
    }
}
function parseLogFile(content, fileName) {
    const ext = (fileName || '').toLowerCase().split('.').pop() || '';
    if (ext === 'out') {
        const result = (0, orcaOutParser_1.parseOrcaOut)(content);
        return { frames: result.frames, title: result.title };
    }
    return (0, logParser_1.parseGaussianLog)(content);
}
function parseLogAsSingleFrame(content) {
    const result = (0, logParser_1.parseGaussianLog)(content);
    if (result.frames.length > 0) {
        return {
            atoms: result.frames[0].atoms,
            bonds: result.frames[0].bonds,
            title: result.frames[0].title,
            hasExplicitBonds: result.frames[0].hasExplicitBonds
        };
    }
    return { atoms: [], bonds: [], title: 'Empty', hasExplicitBonds: false };
}
function tryAutoParse(content, fileName = '') {
    const lines = content.split(/\r?\n/).filter(l => l.trim() !== '');
    if (content.includes('$coord')) {
        return (0, coordParser_1.parseCoord)(content);
    }
    if (content.match(/\*\s*xyz/i) || content.match(/\*\s*xyzfile/i)) {
        return (0, orcaInpParser_1.parseOrcainp)(content);
    }
    if (content.match(/^(ATOM|HETATM)/m)) {
        return (0, pdbParser_1.parsePdb)(content);
    }
    if (content.match(/CARTESIAN COORDINATES/i) && content.match(/MOPAC/i)) {
        return (0, mopacParser_1.parseMopac)(content);
    }
    if (content.match(/^data_/m) && content.includes('_cell_length_a')) {
        return (0, cifParser_1.parseCif)(content);
    }
    if (fileName.match(/\.vasp$/i) || fileName.match(/poscar/i) || fileName.match(/contcar/i)) {
        return (0, vaspParser_1.parseVasp)(content);
    }
    if (fileName.match(/\.cube$/i)) {
        return (0, cubeParser_1.parseCube)(content);
    }
    if (fileName.match(/\.vesta$/i) || content.includes('#VESTA_FORMAT_VERSION')) {
        return (0, vestaParser_1.parseVesta)(content);
    }
    if (content.includes('CARTESIAN COORDINATES (ANGSTROEM)')) {
        const result = (0, orcaOutParser_1.parseOrcaOut)(content);
        if (result.frames.length > 0) {
            return {
                atoms: result.frames[0].atoms,
                bonds: result.frames[0].bonds,
                title: result.frames[0].title,
                hasExplicitBonds: result.frames[0].hasExplicitBonds,
                charge: result.charge,
                multiplicity: result.multiplicity
            };
        }
    }
    if (lines.length > 0) {
        const firstLine = lines[0].trim();
        const possibleCount = parseInt(firstLine, 10);
        if (!isNaN(possibleCount) && possibleCount > 0 && possibleCount < 100000) {
            return (0, xyzParser_1.parseXyz)(content);
        }
    }
    if (content.includes('@<TRIPOS>')) {
        return (0, mol2Parser_1.parseMol2)(content);
    }
    if (content.includes('Standard orientation:') || content.includes('Input orientation:')) {
        return parseLogAsSingleFrame(content);
    }
    if (content.includes('--Link1--') || content.match(/#\s*[A-Za-z]/)) {
        return (0, gjfParser_1.parseGjf)(content);
    }
    return (0, xyzParser_1.parseXyz)(content);
}
//# sourceMappingURL=index.js.map