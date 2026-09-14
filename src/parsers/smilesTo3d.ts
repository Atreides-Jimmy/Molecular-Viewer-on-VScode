/**
 * SMILES -> 3D 分子坐标生成器 (零第三方依赖)
 * ================================================================
 *
 * 输入一条 SMILES 字符串, 输出带显式氢的三维坐标与完整键级信息。
 * 只用标准 Math, 不引入任何 npm 包。
 *
 * 契约:
 * 1. 健壮: 任何失败路径都返回结构化错误码 + 中文说明, 不抛异常、不打印、
 *    不写文件; 所有入口都做边界检查 (非法字符、括号/环号不配对、价态异常、
 *    输入过长、原子数超上限)。
 * 2. 全面: 支持 OpenSMILES 有机子集与常用扩展 (分支、环闭合 %nn、芳香性、
 *    显式键、点断开、同位素、电荷、手性、多片段)。
 * 3. 确定: 同一输入恒得同一坐标 (全程无随机数)。
 *
 * 分层结构:
 *   1. 元素数据表          ELEMENT_MASS / ELEMENT_VALENCE / VDW / COV
 *   2. 常量与错误码         ERR_*
 *   3. 向量与角度工具       vAdd / vSub / vScale / vDot / vCross / vLen / vDist / vNorm / modPositive
 *   4. SMILES 解析器        parseSmiles
 *   5. 价键与隐式氢         computeImplicitHydrogens / checkValence / expandHydrogens
 *   6. 环感知              findRings
 *   7. 杂化与几何参数       atomHybridization / bondAngleAt
 *   8. 三维坐标生成         generateCoordinates
 *   9. 结构松弛            relax / enforceDihedralTargets
 *  10. 手性修正            fixChirality
 *  11. 公开 API            analyzeSmiles / smilesTo3dSafe / validateSmiles / toXyz
 *
 * 坐标系约定: 右手系, 单位埃 (Angstrom)。
 * 算法说明: 生成树 + 二面角扫描给初值, 再用"全距离约束松弛"(1-2 键长 /
 * 1-3 键角 / 非键排斥) 迭代收敛, 顺反与共轭二面角另用整支刚性旋转精确落位。
 * 全部约束都表达为距离约束, 因此求解器只需一种投影算子, 实现短且稳定。
 */

// =====================================================================
// 1. 元素数据表
// =====================================================================

/** 原子量 (IUPAC 2021 相对原子质量, 取常用值) */
const ELEMENT_MASS: { [key: string]: number } = {
    H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81,
    C: 12.011, N: 14.007, O: 15.999, F: 18.998, Ne: 20.180,
    Na: 22.990, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974,
    S: 32.06, Cl: 35.45, Ar: 39.948, K: 39.098, Ca: 40.078,
    Sc: 44.956, Ti: 47.867, V: 50.942, Cr: 51.996, Mn: 54.938,
    Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38,
    Ga: 69.723, Ge: 72.630, As: 74.922, Se: 78.971, Br: 79.904,
    Kr: 83.798, Rb: 85.468, Sr: 87.62, Y: 88.906, Zr: 91.224,
    Nb: 92.906, Mo: 95.95, Tc: 98.0, Ru: 101.07, Rh: 102.91,
    Pd: 106.42, Ag: 107.87, Cd: 112.41, In: 114.82, Sn: 118.71,
    Sb: 121.76, Te: 127.60, I: 126.90, Xe: 131.29, Cs: 132.91,
    Ba: 137.33, La: 138.91, Ce: 140.12, Pr: 140.91, Nd: 144.24,
    Pm: 145.0, Sm: 150.36, Eu: 151.96, Gd: 157.25, Tb: 158.93,
    Dy: 162.50, Ho: 164.93, Er: 167.26, Tm: 168.93, Yb: 173.05,
    Lu: 174.97, Hf: 178.49, Ta: 180.95, W: 183.84, Re: 186.21,
    Os: 190.23, Ir: 192.22, Pt: 195.08, Au: 196.97, Hg: 200.59,
    Tl: 204.38, Pb: 207.2, Bi: 208.98, Po: 209.0, At: 210.0,
    Rn: 222.0, Fr: 223.0, Ra: 226.0, Ac: 227.0, Th: 232.04,
    Pa: 231.04, U: 238.03, Np: 237.0, Pu: 244.0, Am: 243.0,
    Cm: 247.0, Bk: 247.0, Cf: 251.0, Es: 252.0, Fm: 257.0,
    Md: 258.0, No: 259.0, Lr: 262.0, Rf: 267.0, Db: 268.0,
    Sg: 269.0, Bh: 270.0, Hs: 269.0, Mt: 278.0, Ds: 281.0,
    Rg: 281.0, Cn: 285.0, Nh: 286.0, Fl: 289.0, Mc: 289.0,
    Lv: 293.0, Ts: 293.0, Og: 294.0,
};

/** 默认价数 (用于推断隐式氢; 0 表示"不自动补氢") */
const ELEMENT_VALENCE: { [key: string]: number } = {
    H: 1, B: 3, C: 4, N: 3, O: 2, F: 1,
    Si: 4, P: 3, S: 2, Cl: 1, Ge: 4, As: 3,
    Se: 2, Br: 1, Sn: 4, Sb: 3, Te: 2, I: 1,
};

/** 元素可接受的最大价数 */
const ELEMENT_MAX_VALENCE: { [key: string]: number } = {
    H: 1, B: 4, C: 4, N: 4, O: 3, F: 1,
    Si: 6, P: 6, S: 6, Cl: 7, Ge: 4, As: 5,
    Se: 6, Br: 7, Sn: 4, Sb: 5, Te: 6, I: 7,
};

/** 范德华半径 (Bondi 1964 + Batsanov 2001 补充), 单位 Å */
const VDW_RADIUS: { [key: string]: number } = {
    H: 1.20, He: 1.40, Li: 1.82, Be: 1.53, B: 1.92,
    C: 1.70, N: 1.55, O: 1.52, F: 1.47, Ne: 1.54,
    Na: 2.27, Mg: 1.73, Al: 1.84, Si: 2.10, P: 1.80,
    S: 1.80, Cl: 1.75, Ar: 1.88, K: 2.75, Ca: 2.31,
    Ti: 2.11, Cr: 2.06, Mn: 2.05, Fe: 2.04, Co: 2.00,
    Ni: 1.97, Cu: 1.96, Zn: 2.01, Ga: 1.87, Ge: 2.11,
    As: 1.85, Se: 1.90, Br: 1.85, Kr: 2.02, Rb: 3.03,
    Sr: 2.49, Zr: 2.23, Nb: 2.18, Mo: 2.17, Ru: 2.07,
    Rh: 2.02, Pd: 2.05, Ag: 2.03, Cd: 2.18, In: 1.93,
    Sn: 2.17, Sb: 2.06, Te: 2.06, I: 1.98, Xe: 2.16,
    Cs: 3.43, Ba: 2.68, Pt: 2.13, Au: 2.14, Hg: 2.23,
    Tl: 1.96, Pb: 2.02, Bi: 2.07, U: 1.86,
};

/** 共价半径 (Cordero 2008), 单位 Å */
const COVALENT_RADIUS: { [key: string]: number } = {
    H: 0.31, He: 0.28, Li: 1.28, Be: 0.96, B: 0.84,
    C: 0.76, N: 0.71, O: 0.66, F: 0.57, Ne: 0.58,
    Na: 1.66, Mg: 1.41, Al: 1.21, Si: 1.11, P: 1.07,
    S: 1.05, Cl: 1.02, Ar: 1.06, K: 2.03, Ca: 1.76,
    Sc: 1.70, Ti: 1.60, V: 1.53, Cr: 1.39, Mn: 1.39,
    Fe: 1.32, Co: 1.26, Ni: 1.24, Cu: 1.32, Zn: 1.22,
    Ga: 1.22, Ge: 1.20, As: 1.19, Se: 1.20, Br: 1.20,
    Kr: 1.16, Rb: 2.20, Sr: 1.95, Y: 1.90, Zr: 1.75,
    Nb: 1.64, Mo: 1.54, Ru: 1.46, Rh: 1.42, Pd: 1.39,
    Ag: 1.45, Cd: 1.44, In: 1.42, Sn: 1.39, Sb: 1.39,
    Te: 1.38, I: 1.39, Xe: 1.40, Cs: 2.44, Ba: 2.15,
    La: 2.07, Ce: 2.04, Hf: 1.75, Ta: 1.70, W: 1.62,
    Re: 1.51, Os: 1.44, Ir: 1.41, Pt: 1.36, Au: 1.36,
    Hg: 1.32, Tl: 1.45, Pb: 1.46, Bi: 1.48, U: 1.96,
};

/**
 * 显式键长表 (Å), 优先于共价半径估算。
 * 键为 "元素A|元素B|键级" 字符串; 查表时两个方向都会试。
 */
const BOND_LENGTH_OVERRIDE: { [key: string]: number } = (function () {
    const push = (tbl: { [key: string]: number }, a: string, b: string, order: number, len: number) => {
        tbl[a + '|' + b + '|' + order] = len;
    };
    const t: { [key: string]: number } = {};
    push(t, 'C', 'C', 1, 1.54); push(t, 'C', 'C', 2, 1.34); push(t, 'C', 'C', 3, 1.20);
    push(t, 'C', 'N', 1, 1.47); push(t, 'C', 'N', 2, 1.28); push(t, 'C', 'N', 3, 1.16);
    push(t, 'C', 'O', 1, 1.43); push(t, 'C', 'O', 2, 1.21); push(t, 'C', 'O', 3, 1.13);
    push(t, 'C', 'S', 1, 1.81); push(t, 'C', 'S', 2, 1.60); push(t, 'C', 'F', 1, 1.35);
    push(t, 'C', 'Cl', 1, 1.77); push(t, 'C', 'Br', 1, 1.94); push(t, 'C', 'I', 1, 2.14);
    push(t, 'C', 'H', 1, 1.09); push(t, 'C', 'P', 1, 1.84); push(t, 'C', 'Si', 1, 1.87);
    push(t, 'C', 'B', 1, 1.56); push(t, 'C', 'Se', 1, 1.97); push(t, 'C', 'Te', 1, 2.07);
    push(t, 'N', 'N', 1, 1.45); push(t, 'N', 'N', 2, 1.25); push(t, 'N', 'N', 3, 1.10);
    push(t, 'N', 'O', 1, 1.40); push(t, 'N', 'O', 2, 1.22); push(t, 'N', 'H', 1, 1.01);
    push(t, 'N', 'S', 1, 1.68); push(t, 'N', 'P', 1, 1.70); push(t, 'N', 'B', 1, 1.44);
    push(t, 'O', 'O', 1, 1.48); push(t, 'O', 'O', 2, 1.21); push(t, 'O', 'H', 1, 0.96);
    push(t, 'O', 'S', 2, 1.44); push(t, 'O', 'P', 1, 1.60); push(t, 'O', 'P', 2, 1.48);
    push(t, 'O', 'B', 1, 1.36); push(t, 'O', 'Si', 1, 1.63); push(t, 'O', 'Se', 1, 1.61);
    push(t, 'S', 'S', 1, 2.05); push(t, 'S', 'S', 2, 1.89); push(t, 'S', 'H', 1, 1.34);
    push(t, 'S', 'P', 1, 2.10); push(t, 'S', 'C', 2, 1.60); push(t, 'S', 'N', 2, 1.60);
    push(t, 'P', 'H', 1, 1.42); push(t, 'P', 'P', 1, 2.22); push(t, 'P', 'F', 1, 1.55);
    push(t, 'Se', 'H', 1, 1.46); push(t, 'Se', 'Se', 1, 2.34); push(t, 'Te', 'H', 1, 1.69);
    push(t, 'Si', 'H', 1, 1.48); push(t, 'Si', 'O', 1, 1.63); push(t, 'Si', 'F', 1, 1.57);
    push(t, 'Si', 'Cl', 1, 2.03); push(t, 'B', 'H', 1, 1.19); push(t, 'B', 'N', 1, 1.44);
    push(t, 'F', 'F', 1, 1.42); push(t, 'Cl', 'Cl', 1, 1.99); push(t, 'Cl', 'H', 1, 1.27);
    push(t, 'Br', 'Br', 1, 2.28); push(t, 'Br', 'H', 1, 1.41); push(t, 'I', 'I', 1, 2.67);
    push(t, 'I', 'H', 1, 1.61); push(t, 'F', 'H', 1, 0.92);
    return t;
})();

/** 芳香键的键长收缩系数 (相对同元素单键) */
const AROMATIC_BOND_SCALE = 0.90;

/** 键级对共价半径和的收缩系数 (查表未命中时使用) */
const BOND_ORDER_SCALE: { [key: number]: number } = { 1: 1.00, 2: 0.87, 3: 0.78 };

const DEFAULT_BOND_LENGTH = 1.50;
const DEFAULT_VDW_RADIUS = 1.80;
const DEFAULT_COVALENT_RADIUS = 1.50;

/** OpenSMILES 允许不加方括号的有机子集 (保留常数, 判定逻辑用 AROMATIC_ALIAS) */
const ORGANIC_SUBSET = [
    'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I',
    'b', 'c', 'n', 'o', 'p', 's',
];
void ORGANIC_SUBSET;

/** 可以小写书写 (芳香) 的元素, 映射到标准元素符号 */
const AROMATIC_ALIAS: { [key: string]: string } = {
    b: 'B', c: 'C', n: 'N', o: 'O', p: 'P', s: 'S',
    se: 'Se', te: 'Te', as: 'As', si: 'Si',
};

/** 芳香原子在计算隐式氢时的目标价数 */
const AROMATIC_VALENCE: { [key: string]: number } = {
    B: 3, C: 4, N: 3, O: 2, P: 3, S: 2,
    As: 3, Se: 2, Te: 2,
};

/** 芳香键在统计价数时的贡献 (1.5) 的两倍整数值 */
const AROMATIC_BOND_ORDER_X2 = 3;

/** 元素理论杂化键角表: 元素 -> [sp3, sp2, sp] */
const BOND_ANGLE_TABLE: { [key: string]: number[] } = {
    C: [109.47, 120.0, 180.0],
    N: [106.5, 120.0, 180.0],
    O: [104.5, 120.0, 180.0],
    S: [98.0, 106.0, 180.0],
    P: [100.0, 106.0, 180.0],
    Si: [109.5, 120.0, 180.0],
    B: [109.5, 120.0, 180.0],
    Se: [96.0, 106.0, 180.0],
    Te: [94.0, 102.0, 180.0],
    As: [98.0, 106.0, 180.0],
    H: [109.47, 120.0, 180.0],
};
const DEFAULT_BOND_ANGLES = [109.47, 120.0, 180.0];

/** 芳香环内角覆盖 (按环大小) */
const AROMATIC_RING_ANGLE: { [key: number]: number } = { 3: 60.0, 4: 90.0, 5: 108.0, 6: 120.0, 7: 128.6, 8: 135.0 };

/** 小环 (非芳香) 的环内角 */
const RING_ANGLE: { [key: number]: number } = { 3: 60.0, 4: 90.0, 5: 108.0, 6: 120.0, 7: 128.6, 8: 135.0, 9: 140.0 };
void RING_ANGLE;

/** 元素原子量; 未知元素返回 0 */
function elementMass(symbol: string): number {
    const v = ELEMENT_MASS[symbol];
    return v === undefined ? 0.0 : v;
}

/** 元素范德华半径 (表中缺失时用默认值)。 */
function vdwRadius(symbol: string): number {
    const v = VDW_RADIUS[symbol];
    return v === undefined ? DEFAULT_VDW_RADIUS : v;
}

/** 元素共价半径 (表中缺失时用默认值)。 */
function covalentRadius(symbol: string): number {
    const v = COVALENT_RADIUS[symbol];
    return v === undefined ? DEFAULT_COVALENT_RADIUS : v;
}

/** 元素默认价数 (用于推断隐式氢; 0 表示不自动补氢)。 */
function baseValence(symbol: string): number {
    const v = ELEMENT_VALENCE[symbol];
    return v === undefined ? 0 : v;
}

/** 元素可接受的最大价数; 未收录元素退回默认价数。 */
function maxValence(symbol: string): number {
    const v = ELEMENT_MAX_VALENCE[symbol];
    if (v !== undefined) {
        return v;
    }
    return baseValence(symbol);
}

/**
 * 键长 (Å)。order: 1/2/3; 0 或 undefined 表示芳香键。
 * 查找顺序: 显式表 -> 共价半径和 x 键级系数。
 */
function bondLength(a: string, b: string, order?: number): number {
    const aromatic = (order === undefined || order === null || order === 0);
    const keyOrder = aromatic ? 1 : (order as number);
    let length: number;
    const k1 = a + '|' + b + '|' + keyOrder;
    const k2 = b + '|' + a + '|' + keyOrder;
    if (BOND_LENGTH_OVERRIDE[k1] !== undefined) {
        length = BOND_LENGTH_OVERRIDE[k1];
    } else if (BOND_LENGTH_OVERRIDE[k2] !== undefined) {
        length = BOND_LENGTH_OVERRIDE[k2];
    } else {
        const scale = BOND_ORDER_SCALE[keyOrder] === undefined ? 1.0 : BOND_ORDER_SCALE[keyOrder];
        length = (covalentRadius(a) + covalentRadius(b)) * scale;
        if (length <= 0.0) {
            length = DEFAULT_BOND_LENGTH;
        }
    }
    if (aromatic) {
        length = length * AROMATIC_BOND_SCALE;
    }
    return length;
}

/** 元素的理论杂化键角三元组 (sp3, sp2, sp); 未收录元素用通用值。 */
function bondAnglesFor(symbol: string): number[] {
    const v = BOND_ANGLE_TABLE[symbol];
    return v === undefined ? DEFAULT_BOND_ANGLES : v;
}

// =====================================================================
// 2. 常量与错误码
// =====================================================================

const ERR_EMPTY = 'empty_input';                 // 输入为空 / 非字符串
const ERR_CHAR = 'invalid_character';            // 非法字符
const ERR_BRACKET = 'invalid_bracket_atom';      // [] 内容非法
const ERR_ELEMENT = 'unknown_element';           // 未知元素符号
const ERR_PAREN = 'unbalanced_parenthesis';      // 括号不匹配
const ERR_BRACKET_UNCLOSED = 'unclosed_bracket'; // 方括号未闭合
const ERR_RING_OPEN = 'unclosed_ring';           // 环号未闭合
const ERR_RING_FORMAT = 'invalid_ring_number';   // 环号格式错误
const ERR_RING_DUP = 'duplicate_ring_bond';      // 环闭合到已成键的原子
const ERR_RING_SELF = 'ring_closure_on_self';    // 环闭合到自身
const ERR_NO_ATOM = 'no_atoms';                  // 没有原子
const ERR_VALENCE = 'valence_error';             // 价态错误
const ERR_TOO_LARGE = 'too_many_atoms';          // 超出规模上限

const MAX_ATOMS_HARD = 2000;        // 硬上限: 超过直接拒绝
const MAX_INPUT_LEN = 100000;       // 输入字符串长度上限
const MAX_CONSTRAINT_STEP = 0.6;    // 单次距离投影允许修正的最大偏差 (Å)
const RELAX_FACTOR = 0.5;           // 距离投影的松弛因子

/** 默认选项: 原子数上限、松弛轮数与阈值、是否加氢与居中、非键排斥强度、二面角扫描步长、片段间距。 */
export interface SmilesOptions {
    max_atoms: number;
    relax_iterations: number;
    relax_tolerance: number;
    optimize: boolean;
    stereo: boolean;
    add_hydrogens: boolean;
    center: boolean;
    contact_weight: number;
    dihedral_step: number;
    fragment_gap: number;
}

export const SMILES_DEFAULT_OPTIONS: SmilesOptions = {
    max_atoms: MAX_ATOMS_HARD,
    relax_iterations: 500,
    relax_tolerance: 0.002,
    optimize: true,
    stereo: true,
    add_hydrogens: true,
    center: true,
    contact_weight: 0.45,
    dihedral_step: 15.0,
    fragment_gap: 4.0,
};

/**
 * 合并用户选项与默认值 (返回新对象, 不修改默认值)。
 * 只接受已知字段, 未知字段被忽略。
 */
function makeOptions(overrides?: Partial<SmilesOptions>): SmilesOptions {
    const opts: SmilesOptions = {
        max_atoms: SMILES_DEFAULT_OPTIONS.max_atoms,
        relax_iterations: SMILES_DEFAULT_OPTIONS.relax_iterations,
        relax_tolerance: SMILES_DEFAULT_OPTIONS.relax_tolerance,
        optimize: SMILES_DEFAULT_OPTIONS.optimize,
        stereo: SMILES_DEFAULT_OPTIONS.stereo,
        add_hydrogens: SMILES_DEFAULT_OPTIONS.add_hydrogens,
        center: SMILES_DEFAULT_OPTIONS.center,
        contact_weight: SMILES_DEFAULT_OPTIONS.contact_weight,
        dihedral_step: SMILES_DEFAULT_OPTIONS.dihedral_step,
        fragment_gap: SMILES_DEFAULT_OPTIONS.fragment_gap,
    };
    if (overrides) {
        const keys = Object.keys(overrides) as (keyof SmilesOptions)[];
        for (let i = 0; i < keys.length; i++) {
            const key = keys[i];
            const value = overrides[key];
            if (value === undefined) {
                continue;
            }
            if (typeof value === typeof opts[key]) {
                (opts as unknown as { [k: string]: unknown })[key as string] = value;
            }
        }
    }
    return opts;
}

/** 解析错误信息模板 (中文, 与错误码一一对应) */
const ERROR_TEXT: { [key: string]: string } = {
    [ERR_EMPTY]: '输入为空',
    [ERR_CHAR]: '非法字符',
    [ERR_BRACKET]: '方括号原子内容非法',
    [ERR_ELEMENT]: '未知元素符号',
    [ERR_PAREN]: '括号不匹配',
    [ERR_BRACKET_UNCLOSED]: '方括号未闭合',
    [ERR_RING_OPEN]: '存在未闭合的环',
    [ERR_RING_FORMAT]: '环编号格式错误',
    [ERR_RING_DUP]: '环闭合连接的两个原子之间已存在键',
    [ERR_RING_SELF]: '环闭合不能连接到自身',
    [ERR_NO_ATOM]: '没有解析出任何原子',
    [ERR_VALENCE]: '原子价态异常',
    [ERR_TOO_LARGE]: '原子数超过上限',
};

// =====================================================================
// 3. 向量工具 (数组形式 [x, y, z])
// =====================================================================

type Vec3 = number[];

/**
 * 角度取模: 把 x 回绕到 [0, m), 结果恒非负 (modPositive(-30, 360) == 330)。
 * 不能直接用 JS 的 %: 它保留被除数的符号, 负角度会拿到负余数。
 */
function modPositive(x: number, m: number): number {
    const r = x % m;
    return r < 0 ? r + m : r;
}

/** 向量加法 a + b。 */
function vAdd(a: Vec3, b: Vec3): Vec3 {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

/** 向量减法 a - b。 */
function vSub(a: Vec3, b: Vec3): Vec3 {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

/** 向量数乘 a * s。 */
function vScale(a: Vec3, s: number): Vec3 {
    return [a[0] * s, a[1] * s, a[2] * s];
}

/** 向量点积 a · b。 */
function vDot(a: Vec3, b: Vec3): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/** 向量叉积 a × b。 */
function vCross(a: Vec3, b: Vec3): Vec3 {
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]];
}

/** 向量长度 |a|。 */
function vLen(a: Vec3): number {
    return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/** 两点距离 |a - b|。 */
function vDist(a: Vec3, b: Vec3): number {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    const dz = a[2] - b[2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/** 单位向量; 零向量退化为 +X 方向, 避免除零。 */
function vNorm(a: Vec3): Vec3 {
    const n = vLen(a);
    if (n < 1e-12) {
        return [1.0, 0.0, 0.0];
    }
    return [a[0] / n, a[1] / n, a[2] / n];
}

/** 计算 p0-p1-p2-p3 二面角, 返回 (-180, 180] 度 */
function dihedralAngle(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3): number {
    const b0 = vSub(p1, p0);
    const b1 = vSub(p2, p1);
    const b2 = vSub(p3, p2);
    const n1 = vCross(b0, b1);
    const n2 = vCross(b1, b2);
    const length1 = vLen(n1);
    const length2 = vLen(n2);
    if (length1 < 1e-9 || length2 < 1e-9) {
        return 0.0;
    }
    const m = vCross(vNorm(n1), vNorm(b1));
    let x = vDot(n1, n2) / (length1 * length2);
    if (x > 1.0) {
        x = 1.0;
    } else if (x < -1.0) {
        x = -1.0;
    }
    let ang = Math.acos(x);
    if (vDot(m, n2) < 0.0) {
        ang = -ang;
    }
    return ang * 180.0 / Math.PI;
}

/** 三点夹角 (度) */
function angleBetween(a: Vec3, center: Vec3, b: Vec3): number {
    const u = vSub(a, center);
    const v = vSub(b, center);
    const lu = vLen(u);
    const lv = vLen(v);
    if (lu < 1e-9 || lv < 1e-9) {
        return 0.0;
    }
    let x = vDot(u, v) / (lu * lv);
    if (x > 1.0) {
        x = 1.0;
    } else if (x < -1.0) {
        x = -1.0;
    }
    return Math.acos(x) * 180.0 / Math.PI;
}
void angleBetween;

/** 绕经过 center、方向为 axis 的轴旋转 point (Rodrigues 公式) */
function rotateAboutAxis(point: Vec3, axis: Vec3, angleDeg: number, center: Vec3): Vec3 {
    const rad = angleDeg * Math.PI / 180.0;
    const c = Math.cos(rad);
    const s = Math.sin(rad);
    const k = vNorm(axis);
    const p = vSub(point, center);
    const term1 = vScale(p, c);
    const term2 = vScale(vCross(k, p), s);
    const term3 = vScale(k, vDot(k, p) * (1.0 - c));
    return vAdd(center, vAdd(vAdd(term1, term2), term3));
}

// =====================================================================
// 4. SMILES 解析器
// =====================================================================
//
// 输出分子图 mol 的结构:
//   mol.atoms[i] = { element, aromatic, charge, isotope, hcount, bracket,
//                    chirality, stereo_h_first };  hcount = -1 表示待推断
//   mol.bonds[k] = { a, b, order, aromatic, direction }
//   mol.adj[i]   = [k, ...];  与原子 i 相连的键索引, 按键入顺序排列
//
// 关键设计: 邻接只保存"键索引", 每条键在 bonds 中只出现一次, 因此从结构上
// 排除了"邻接表与键表两套数据可能不一致"以及重复加键的问题。

interface SmilesAtom {
    element: string;
    aromatic: boolean;
    charge: number;
    isotope: number;
    hcount: number;
    bracket: boolean;
    chirality: string;
    stereo_h_first: boolean;
}

interface SmilesBond {
    a: number;
    b: number;
    order: number;
    aromatic: boolean;
    direction: string;
}

interface Mol {
    atoms: SmilesAtom[];
    bonds: SmilesBond[];
    adj: number[][];
    heavy_count?: number;
}

const SINGLE_BOND = 1;
const DOUBLE_BOND = 2;
const TRIPLE_BOND = 3;
const AROMATIC_BOND = 0;

/** 允许不加方括号书写的两字母元素符号 */
const BARE_TWO_LETTER: { [key: string]: boolean } = { Cl: true, Br: true, Si: true, Se: true, Te: true };

/** 新建空的分子图 (atoms / bonds / adj)。 */
function newMol(): Mol {
    return { atoms: [], bonds: [], adj: [] };
}

/** 追加一个原子并扩展邻接表, 返回其索引。 */
function addAtom(mol: Mol, element: string, aromatic: boolean, charge: number, isotope: number,
                hcount: number, bracket: boolean, chirality: string): number {
    const atom: SmilesAtom = {
        element: element,
        aromatic: aromatic,
        charge: charge,
        isotope: isotope,
        hcount: hcount,
        bracket: bracket,
        chirality: chirality,
        stereo_h_first: false,
    };
    mol.atoms.push(atom);
    mol.adj.push([]);
    return mol.atoms.length - 1;
}

/** 追加一条键, 返回其索引; 两端邻接表各记录该键索引一次。 */
function addBond(mol: Mol, a: number, b: number, order: number, aromatic: boolean, direction: string): number {
    const bond: SmilesBond = {
        a: a, b: b,
        order: order,
        aromatic: aromatic,
        direction: direction,
    };
    mol.bonds.push(bond);
    const index = mol.bonds.length - 1;
    mol.adj[a].push(index);
    mol.adj[b].push(index);
    return index;
}

/** 若 a-b 之间已有键则返回键索引, 否则返回 -1 */
function bondBetween(mol: Mol, a: number, b: number): number {
    const list = mol.adj[a];
    for (let i = 0; i < list.length; i++) {
        const k = list[i];
        const bond = mol.bonds[k];
        if ((bond.a === a && bond.b === b) || (bond.a === b && bond.b === a)) {
            return k;
        }
    }
    return -1;
}

interface ParseResult {
    ok: boolean;
    error: string;
    error_pos: number;
    error_text: string;
    mol: Mol | null;
}

/** 填写失败结果 (错误码 + 出错位置 + 中文说明), 返回同一个结果对象。 */
function fail(result: ParseResult, code: string, pos: number, detail?: string): ParseResult {
    result.ok = false;
    result.error = code;
    result.error_pos = pos;
    const base = ERROR_TEXT[code] === undefined ? code : ERROR_TEXT[code];
    const tail = detail === undefined ? '' : detail;
    if (pos >= 0) {
        result.error_text = '位置 ' + pos + ': ' + base + tail;
    } else {
        result.error_text = base + tail;
    }
    return result;
}

/** ASCII 字母判定 (刻意不用 Unicode 判定) */
function isAlpha(ch: string): boolean {
    return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z');
}

/** ASCII 小写字母判定。 */
function isLower(ch: string): boolean {
    return ch >= 'a' && ch <= 'z';
}

/** ASCII 数字判定 */
function isDigit(ch: string): boolean {
    return ch >= '0' && ch <= '9';
}

interface AtomFields {
    element: string;
    aromatic: boolean;
    charge: number;
    isotope: number;
    hcount: number;
    chirality: string;
}

type AtomParseOutcome = { ok: true; fields: AtomFields; next: number } | { ok: false; code: string };

/**
 * 解析 [ ] 原子。语法 (OpenSMILES):
 *   '[' [同位素] 元素 [手性] [H数量] [电荷] [':原子类'] ']'
 */
function parseBracketAtom(text: string, pos: number): AtomParseOutcome {
    const n = text.length;
    let i = pos + 1;
    let isotope = 0;
    let element = '';
    let chirality = '';
    let hcount = 0;
    let charge = 0;

    // --- 同位素 (可省略) ---
    const start = i;
    while (i < n && isDigit(text[i])) {
        i += 1;
    }
    if (i > start) {
        isotope = parseInt(text.substring(start, i), 10);
    }

    // --- 元素符号 ---
    if (i >= n) {
        return { ok: false, code: ERR_BRACKET_UNCLOSED };
    }
    const ch = text[i];
    if (ch === '*') {
        element = '*';
        i += 1;
    } else if (isAlpha(ch)) {
        const first = text[i];
        i += 1;
        element = '';
        // SMILES 规定: 两字母元素符号的第二个字母必须小写 (Cl / Na / Si / Se)。
        // 因此 "[NH4+]" 中的 "NH" 不是元素 Nh, 而是 N 加 4 个氢。
        if (i < n && isLower(text[i])) {
            const two = canonicalElement(first + text[i]);
            if (two.length === 2 && ELEMENT_MASS[two] !== undefined) {
                element = two;
                i += 1;
            }
        }
        if (!element) {
            element = canonicalElement(first);
        }
        if (ELEMENT_MASS[element] === undefined && element !== '*') {
            return { ok: false, code: ERR_ELEMENT };
        }
    } else {
        return { ok: false, code: ERR_BRACKET };
    }

    // --- 手性 @ / @@ ---
    if (i < n && text[i] === '@') {
        i += 1;
        if (i < n && text[i] === '@') {
            chirality = '@@';
            i += 1;
        } else {
            chirality = '@';
        }
    }

    // --- 氢数量 H / Hn ---
    if (i < n && text[i] === 'H') {
        i += 1;
        if (i < n && isDigit(text[i])) {
            hcount = parseInt(text[i], 10);
            i += 1;
        } else {
            hcount = 1;
        }
    }

    // --- 电荷 + / - / ++ / -- / +n / -n ---
    if (i < n && (text[i] === '+' || text[i] === '-')) {
        const sign = text[i] === '+' ? 1 : -1;
        const symbol = text[i];
        i += 1;
        if (i < n && isDigit(text[i])) {
            charge = sign * parseInt(text[i], 10);
            i += 1;
        } else {
            charge = sign;
            while (i < n && text[i] === symbol) {
                charge += sign;
                i += 1;
            }
        }
    }

    // --- 原子类 :n (仅记录, 不影响几何) ---
    if (i < n && text[i] === ':') {
        i += 1;
        const clsStart = i;
        while (i < n && isDigit(text[i])) {
            i += 1;
        }
        if (i === clsStart) {
            return { ok: false, code: ERR_BRACKET };
        }
    }

    if (i >= n) {
        return { ok: false, code: ERR_BRACKET_UNCLOSED };
    }
    if (text[i] !== ']') {
        return { ok: false, code: ERR_BRACKET };
    }

    const isAromatic = (AROMATIC_VALENCE[element] !== undefined) && sourceWasLowercase(text, pos + 1, i);
    const fields: AtomFields = {
        element: element,
        aromatic: isAromatic,
        charge: charge,
        isotope: isotope,
        hcount: hcount,
        chirality: chirality,
    };
    return { ok: true, fields: fields, next: i + 1 };
}

/** 检查 [ ] 内的元素符号是否以小写书写 (用于判定芳香性) */
function sourceWasLowercase(text: string, start: number, end: number): boolean {
    for (let k = start; k < end; k++) {
        if (isAlpha(text[k])) {
            return isLower(text[k]);
        }
    }
    return false;
}

/**
 * 把元素符号归一化为标准写法 (大小写敏感, 不猜测)。
 */
function canonicalElement(symbol: string): string {
    if (ELEMENT_MASS[symbol] !== undefined) {
        return symbol;
    }
    const lowered = symbol.toLowerCase();
    if (AROMATIC_ALIAS[lowered] !== undefined) {
        return AROMATIC_ALIAS[lowered];
    }
    if (symbol.length === 1 && ELEMENT_MASS[symbol.toUpperCase()] !== undefined) {
        return symbol.toUpperCase();
    }
    return symbol;
}

/** 解析不加方括号的原子 */
function parseOrganicAtom(text: string, pos: number): AtomParseOutcome {
    const ch = text[pos];
    if (ch === '*') {
        return {
            ok: true,
            next: pos + 1,
            fields: {
                element: '*', aromatic: false, charge: 0, isotope: 0,
                hcount: -1, chirality: '',
            },
        };
    }

    if (isLower(ch)) {
        // 芳香单字母: b c n o p s
        if (AROMATIC_ALIAS[ch] === undefined || AROMATIC_ALIAS[ch].length !== 1) {
            return { ok: false, code: ERR_CHAR };
        }
        return {
            ok: true,
            next: pos + 1,
            fields: {
                element: AROMATIC_ALIAS[ch], aromatic: true, charge: 0,
                isotope: 0, hcount: -1, chirality: '',
            },
        };
    }

    // 大写: 先尝试两字母
    if (pos + 1 < text.length) {
        const two = ch + text[pos + 1];
        if (BARE_TWO_LETTER[two] === true) {
            if (ELEMENT_MASS[two] === undefined) {
                return { ok: false, code: ERR_ELEMENT };
            }
            return {
                ok: true,
                next: pos + 2,
                fields: {
                    element: two, aromatic: false, charge: 0,
                    isotope: 0, hcount: -1, chirality: '',
                },
            };
        }
    }

    if (ELEMENT_MASS[ch] === undefined) {
        return { ok: false, code: ERR_ELEMENT };
    }
    return {
        ok: true,
        next: pos + 1,
        fields: {
            element: ch, aromatic: false, charge: 0, isotope: 0,
            hcount: -1, chirality: '',
        },
    };
}

/**
 * 解析 SMILES。
 * 返回 { ok, error, error_pos, error_text, mol }。不抛异常、不打印。
 */
export function parseSmiles(smiles: unknown): ParseResult {
    const result: ParseResult = { ok: false, error: '', error_pos: -1, error_text: '', mol: null };

    if (typeof smiles !== 'string') {
        return fail(result, ERR_EMPTY, -1);
    }

    // 去掉所有空白 (含换行): 部分来源会把长 SMILES 折行
    let text = '';
    for (let p = 0; p < smiles.length; p++) {
        const c = smiles[p];
        if (!/\s/.test(c)) {
            text += c;
        }
    }
    if (text.length > MAX_INPUT_LEN) {
        return fail(result, ERR_TOO_LARGE, -1, ' (输入长度 ' + text.length + ')');
    }
    if (!text) {
        return fail(result, ERR_EMPTY, -1);
    }

    const mol = newMol();
    let i = 0;
    const n = text.length;
    let prev = -1;                  // 当前链尾原子
    let pendingOrder = SINGLE_BOND;
    let pendingArom = false;
    let pendingDir = '';
    let pendingSet = false;         // 是否已显式给出键符号
    const branchStack: number[] = [];
    const ringOpen: { [key: number]: [number, number, boolean] } = {};   // 环号 -> (原子索引, 键级, 是否芳香)

    /** 把 current 原子接上前一个原子 */
    const add = (current: number, needBond: boolean): void => {
        if (prev >= 0 && needBond) {
            let order = pendingOrder;
            let aromatic = pendingArom;
            if (aromatic && !(mol.atoms[prev].aromatic && mol.atoms[current].aromatic)) {
                // ':' 用于非芳香原子时退化为单键
                aromatic = false;
                order = SINGLE_BOND;
            }
            addBond(mol, prev, current, order, aromatic, pendingDir);
        }
        prev = current;
    };

    while (i < n) {
        const ch = text[i];

        // ---------- 分支 ----------
        if (ch === '(') {
            if (prev < 0) {
                return fail(result, ERR_PAREN, i, " ('(' 之前没有原子)");
            }
            branchStack.push(prev);
            i += 1;
            continue;
        }

        if (ch === ')') {
            if (branchStack.length === 0) {
                return fail(result, ERR_PAREN, i, " (多余的 ')')");
            }
            prev = branchStack.pop() as number;
            pendingOrder = SINGLE_BOND;
            pendingArom = false;
            pendingDir = '';
            pendingSet = false;
            i += 1;
            continue;
        }

        // ---------- 断开 ----------
        if (ch === '.') {
            if (prev < 0) {
                return fail(result, ERR_CHAR, i, " (多余的 '.')");
            }
            prev = -1;
            pendingOrder = SINGLE_BOND;
            pendingArom = false;
            pendingDir = '';
            pendingSet = false;
            i += 1;
            continue;
        }

        // ---------- 键符号 ----------
        if (ch === '-' || ch === '=' || ch === '#' || ch === ':' || ch === '$' || ch === '/' || ch === '\\') {
            if (prev < 0) {
                return fail(result, ERR_CHAR, i, " (键符号 '" + ch + "' 之前没有原子)");
            }
            if (pendingSet) {
                return fail(result, ERR_CHAR, i, ' (重复的键符号)');
            }
            if (ch === '-') {
                pendingOrder = SINGLE_BOND;
            } else if (ch === '=') {
                pendingOrder = DOUBLE_BOND;
            } else if (ch === '#') {
                pendingOrder = TRIPLE_BOND;
            } else if (ch === ':') {
                pendingOrder = AROMATIC_BOND;
                pendingArom = true;
            } else {
                pendingOrder = SINGLE_BOND;
                pendingDir = ch;
            }
            pendingSet = true;
            i += 1;
            continue;
        }

        // ---------- 环闭合 ----------
        if (isDigit(ch) || ch === '%') {
            let ringId: number;
            if (ch === '%') {
                if (i + 2 >= n || !isDigit(text[i + 1]) || !isDigit(text[i + 2])) {
                    return fail(result, ERR_RING_FORMAT, i, " ('%' 后需要两位数字)");
                }
                ringId = parseInt(text.substring(i + 1, i + 3), 10);
                i += 3;
            } else {
                ringId = parseInt(ch, 10);
                i += 1;
            }
            if (ringId === 0) {
                return fail(result, ERR_RING_FORMAT, i - 1, ' (环编号不能为 0)');
            }

            if (prev < 0) {
                return fail(result, ERR_RING_FORMAT, i - 1, ' (环编号之前没有原子)');
            }

            if (ringOpen[ringId] !== undefined) {
                const entry = ringOpen[ringId];
                delete ringOpen[ringId];
                const other = entry[0];
                const order0 = entry[1];
                const arom0 = entry[2];
                if (other === prev) {
                    return fail(result, ERR_RING_SELF, i - 1);
                }
                if (bondBetween(mol, other, prev) >= 0) {
                    return fail(result, ERR_RING_DUP, i - 1);
                }
                let order = pendingSet ? pendingOrder : order0;
                let arom = pendingArom || arom0;
                if (arom && !(mol.atoms[other].aromatic && mol.atoms[prev].aromatic)) {
                    arom = false;
                    order = SINGLE_BOND;
                }
                addBond(mol, other, prev, order, arom, pendingDir);
            } else {
                // 注意: 环编号在闭合之后可以再次使用 (如 c1ccccc1C1CC1),
                // 因此这里只检查"当前是否已打开", 不记录历史。
                ringOpen[ringId] = [prev, pendingOrder, pendingArom];
            }
            pendingOrder = SINGLE_BOND;
            pendingArom = false;
            pendingDir = '';
            pendingSet = false;
            continue;
        }

        // ---------- 方括号原子 ----------
        if (ch === '[') {
            const outcome = parseBracketAtom(text, i);
            if (!outcome.ok) {
                return fail(result, outcome.code, i);
            }
            const f = outcome.fields;
            const idx = addAtom(mol, f.element, f.aromatic, f.charge, f.isotope,
                                f.hcount, true, f.chirality);
            if (f.chirality && f.hcount > 0) {
                mol.atoms[idx].stereo_h_first = true;
            }
            add(idx, true);
            i = outcome.next;
            pendingOrder = SINGLE_BOND;
            pendingArom = false;
            pendingDir = '';
            pendingSet = false;
            continue;
        }

        // ---------- 普通原子 ----------
        if (isAlpha(ch) || ch === '*') {
            const outcome = parseOrganicAtom(text, i);
            if (!outcome.ok) {
                return fail(result, outcome.code, i);
            }
            const f = outcome.fields;
            const idx = addAtom(mol, f.element, f.aromatic, f.charge, f.isotope,
                                f.hcount, false, f.chirality);
            add(idx, true);
            i = outcome.next;
            pendingOrder = SINGLE_BOND;
            pendingArom = false;
            pendingDir = '';
            pendingSet = false;
            continue;
        }

        // ---------- 其它: 非法字符 ----------
        return fail(result, ERR_CHAR, i, " ('" + ch + "')");
    }

    if (branchStack.length > 0) {
        return fail(result, ERR_PAREN, -1, " (未闭合的 '(' 共 " + branchStack.length + ' 个)');
    }
    const openKeys = Object.keys(ringOpen).map((k) => parseInt(k, 10)).sort((a, b) => a - b);
    if (openKeys.length > 0) {
        // 环号列表按 [1] / [1, 2] 的形式打印, 便于一眼看出还有哪些环号没闭合
        return fail(result, ERR_RING_OPEN, -1, ' (环编号 [' + openKeys.join(', ') + '])');
    }
    if (pendingSet) {
        return fail(result, ERR_CHAR, -1, ' (末尾多余的键符号)');
    }
    if (prev < 0) {
        return fail(result, ERR_CHAR, -1, " (末尾多余的 '.')");
    }
    if (mol.atoms.length === 0) {
        return fail(result, ERR_NO_ATOM, -1);
    }

    // 统一处理芳香性: 芳香原子之间的键若未显式指定, 视为芳香键
    markAromaticBonds(mol);

    result.ok = true;
    result.error = '';
    result.mol = mol;
    return result;
}

/**
 * 把"芳香原子-芳香原子"之间的单键升级为芳香键 (order=0)。
 */
function markAromaticBonds(mol: Mol): void {
    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        const a = mol.atoms[bond.a];
        const b = mol.atoms[bond.b];
        if (a.aromatic && b.aromatic && bond.order === SINGLE_BOND) {
            bond.order = AROMATIC_BOND;
            bond.aromatic = true;
        }
    }
}

/** 键级的 2 倍整数表示: 芳香键为 3 (即 1.5*2), 其余为 order*2 */
function bondOrderX2(bond: SmilesBond): number {
    if (bond.aromatic) {
        return AROMATIC_BOND_ORDER_X2;
    }
    return bond.order * 2;
}

/** 键级的浮点表示 (芳香键为 1.5), 用于几何计算 */
function orderValue(bond: SmilesBond): number {
    if (bond.aromatic) {
        return 1.5;
    }
    return bond.order;
}
void orderValue;

// =====================================================================
// 5. 价键与隐式氢
// =====================================================================

/**
 * 为每个原子计算 hcount: 目标价数 - 已用键级 = 氢数。
 * 芳香键按 1.5 计数 (用整数 3 表示 1.5*2, 避免浮点误差)。
 */
function computeImplicitHydrogens(mol: Mol): void {
    const atoms = mol.atoms;
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        if (atom.bracket) {
            continue;                      // 方括号原子: 氢数已显式给出
        }
        let target: number;
        if (atom.aromatic) {
            const av = AROMATIC_VALENCE[atom.element];
            target = av === undefined ? baseValence(atom.element) : av;
        } else {
            target = baseValence(atom.element);
        }
        let usedX2 = 0;
        const list = mol.adj[i];
        for (let t = 0; t < list.length; t++) {
            usedX2 += bondOrderX2(mol.bonds[list[t]]);
        }
        const diff = target * 2 - usedX2;
        // 向下取整的"四舍五入": floor(x + 0.5), 对负值也安全
        const h = Math.floor(diff / 2.0 + 0.5);
        atom.hcount = h > 0 ? h : 0;
    }
}

/**
 * 检查价态。返回错误信息列表 (空列表表示通过)。
 * 计数规则: 价数 = σ 键数 (芳香键计 1, 其余按键级) + 已连接氢数。
 */
function checkValence(mol: Mol): string[] {
    const problems: string[] = [];
    const atoms = mol.atoms;
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        const element = atom.element;
        if (element === '*') {
            continue;
        }
        let allowed = maxValence(element);
        if (allowed <= 0) {
            continue;                      // 表中未定义最大价数的元素 (金属等) 放行
        }
        if (atom.bracket && allowed >= 5) {
            allowed += 1;                  // 显式书写的超价物种 (如 [Cl](=O)(=O)=O)
        }
        let usedX2 = 0;
        const list = mol.adj[i];
        for (let t = 0; t < list.length; t++) {
            const bond = mol.bonds[list[t]];
            if (bond.aromatic) {
                usedX2 += 2;               // 芳香键按单键计
            } else {
                usedX2 += bond.order * 2;
            }
        }
        const hydrogen = atom.hcount;
        if (hydrogen > 0) {
            usedX2 += hydrogen * 2;
        }
        if (usedX2 > allowed * 2) {
            problems.push('原子 ' + i + ' (' + element + ') 价态异常: ' + (usedX2 / 2.0) + ' > ' + allowed);
        }
    }
    return problems;
}

/**
 * 把隐式氢变成真实原子, 返回扩充后的新分子图。
 * 氢原子追加在重原子之后。
 */
function expandHydrogens(mol: Mol): Mol {
    const out: Mol = {
        atoms: [],
        bonds: [],
        adj: [],
        heavy_count: mol.atoms.length,
    };
    for (let i = 0; i < mol.atoms.length; i++) {
        const atom = mol.atoms[i];
        const newAtom: SmilesAtom = {
            element: atom.element, aromatic: atom.aromatic, charge: atom.charge,
            isotope: atom.isotope, hcount: atom.hcount, bracket: atom.bracket,
            chirality: atom.chirality, stereo_h_first: atom.stereo_h_first,
        };
        out.atoms.push(newAtom);
        out.adj.push([]);
    }
    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        out.bonds.push({ a: bond.a, b: bond.b, order: bond.order, aromatic: bond.aromatic, direction: bond.direction });
        const index = out.bonds.length - 1;
        out.adj[bond.a].push(index);
        out.adj[bond.b].push(index);
    }

    for (let i = 0; i < mol.atoms.length; i++) {
        const count = mol.atoms[i].hcount;
        for (let t = 0; t < count; t++) {
            const hidx = out.atoms.length;
            out.atoms.push({
                element: 'H',
                aromatic: false,
                charge: 0,
                isotope: 0,
                hcount: 0,
                bracket: false,
                chirality: '',
                stereo_h_first: false,
            });
            out.adj.push([]);
            out.bonds.push({ a: i, b: hidx, order: SINGLE_BOND, aromatic: false, direction: '' });
            const bi = out.bonds.length - 1;
            out.adj[i].push(bi);
            out.adj[hidx].push(bi);
        }
    }
    return out;
}

// =====================================================================
// 6. 环感知
// =====================================================================

/**
 * 反复摘除度为 1 的原子; 剩下的就是处在环上的原子。
 * 返回 removed 布尔列表。
 */
function pruneAcyclicAtoms(mol: Mol): boolean[] {
    const n = mol.atoms.length;
    const degree: number[] = new Array<number>(n);
    for (let i = 0; i < n; i++) {
        degree[i] = mol.adj[i].length;
    }
    const removed: boolean[] = new Array<boolean>(n);
    for (let i = 0; i < n; i++) {
        removed[i] = false;
    }
    const stack: number[] = [];
    for (let i = 0; i < n; i++) {
        if (degree[i] <= 1) {
            stack.push(i);
        }
    }
    while (stack.length > 0) {
        const cur = stack.pop() as number;
        if (removed[cur]) {
            continue;
        }
        removed[cur] = true;
        const list = mol.adj[cur];
        for (let t = 0; t < list.length; t++) {
            const bond = mol.bonds[list[t]];
            const other = bond.b === cur ? bond.a : bond.b;
            if (!removed[other]) {
                degree[other] -= 1;
                if (degree[other] <= 1) {
                    stack.push(other);
                }
            }
        }
    }
    return removed;
}

/** 避开 skipBond, 求 src 到 dst 的最短原子路径 (含两端) */
function shortestPath(mol: Mol, src: number, dst: number, skipBond: number): number[] | null {
    const n = mol.atoms.length;
    const prev: number[] = new Array<number>(n);
    const seen: boolean[] = new Array<boolean>(n);
    for (let i = 0; i < n; i++) {
        prev[i] = -1;
        seen[i] = false;
    }
    seen[src] = true;
    const queue: number[] = [src];
    let head = 0;
    let found = (src === dst);
    while (head < queue.length && !found) {
        const cur = queue[head];
        head += 1;
        const list = mol.adj[cur];
        for (let t = 0; t < list.length; t++) {
            const k = list[t];
            if (k === skipBond) {
                continue;
            }
            const bond = mol.bonds[k];
            const nxt = bond.b === cur ? bond.a : bond.b;
            if (!seen[nxt]) {
                seen[nxt] = true;
                prev[nxt] = cur;
                if (nxt === dst) {
                    found = true;
                    break;
                }
                queue.push(nxt);
            }
        }
    }
    if (!found) {
        return null;
    }
    const path: number[] = [];
    let cur = dst;
    while (cur !== -1) {
        path.push(cur);
        if (cur === src) {
            break;
        }
        cur = prev[cur];
    }
    if (path.length === 0 || path[path.length - 1] !== src) {
        return null;
    }
    path.reverse();
    return path;
}

/**
 * 求环集合 (最小环优先)。返回原子索引列表的列表, 每个环按顺序排列。
 */
function findRings(mol: Mol): number[][] {
    const removed = pruneAcyclicAtoms(mol);
    const rings: number[][] = [];
    const seen: { [key: string]: boolean } = {};
    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        const a = bond.a;
        const b = bond.b;
        if (removed[a] || removed[b]) {
            continue;                      // 该键不在任何环上
        }
        const path = shortestPath(mol, a, b, k);
        if (path === null || path.length < 3) {
            continue;
        }
        const key = path.slice().sort((x, y) => x - y).join(',');
        if (seen[key] === true) {
            continue;
        }
        seen[key] = true;
        rings.push(path);
    }
    return rings;
}

interface RingMaps {
    ring_edge: { [key: string]: number };
    ring_size: number[];
}

/**
 * 构建环查询表:
 *   ring_edge["min|max"] = 环索引 (取最小的那个)
 *   ring_size[i]         = 环 i 的原子数
 */
function buildRingMaps(mol: Mol, rings: number[][]): RingMaps {
    void mol;
    const ringEdge: { [key: string]: number } = {};
    const ringSize: number[] = [];
    for (let r = 0; r < rings.length; r++) {
        ringSize.push(rings[r].length);
    }
    // 按环从小到大处理, 使小环优先占据 ring_edge
    const order: number[] = [];
    for (let r = 0; r < rings.length; r++) {
        order.push(r);
    }
    order.sort((x, y) => rings[x].length - rings[y].length);
    for (let oi = 0; oi < order.length; oi++) {
        const r = order[oi];
        const ring = rings[r];
        const count = ring.length;
        for (let i = 0; i < count; i++) {
            const a = ring[i];
            const b = ring[(i + 1) % count];
            const key = a < b ? (a + '|' + b) : (b + '|' + a);
            if (ringEdge[key] === undefined) {
                ringEdge[key] = r;
            }
        }
    }
    return { ring_edge: ringEdge, ring_size: ringSize };
}

/** ring_edge / targets 的键: 两端索引较小者在前, 形如 "a|b"。 */
function ringKey(a: number, b: number): string {
    return a < b ? (a + '|' + b) : (b + '|' + a);
}

// =====================================================================
// 7. 杂化、键角与扭转目标
// =====================================================================

const SP3 = 0;
const SP2 = 1;
const SP = 2;

const HYBRID_NAMES = ['sp3', 'sp2', 'sp'];

/** 可能形成酰胺/共轭型的氮: 邻接原子带双键时按 sp2 处理 (平面) */
const AMIDE_LIKE_ELEMENTS: { [key: string]: boolean } = { N: true, P: true };

/**
 * 判断原子杂化类型, 返回 SP3 / SP2 / SP。
 */
function atomHybridization(mol: Mol, index: number): number {
    const atom = mol.atoms[index];
    const bonds = mol.adj[index];
    const nbonds = bonds.length;
    let aromaticCount = 0;
    let doubleCount = 0;
    let neighborConjugated = false;
    for (let bi = 0; bi < bonds.length; bi++) {
        const k = bonds[bi];
        const bond = mol.bonds[k];
        if (bond.aromatic) {
            aromaticCount += 1;
            continue;
        }
        if (bond.order >= 3) {
            return SP;
        }
        if (bond.order === 2) {
            doubleCount += 1;
        }
        const other = bond.b === index ? bond.a : bond.b;
        const list2 = mol.adj[other];
        for (let t = 0; t < list2.length; t++) {
            const k2 = list2[t];
            const bond2 = mol.bonds[k2];
            if (k2 !== k && (bond2.aromatic || bond2.order >= 2)) {
                neighborConjugated = true;
                break;
            }
        }
    }
    if (doubleCount >= 2 && nbonds === 2) {
        return SP;
    }
    if (doubleCount >= 1 || aromaticCount >= 1) {
        return SP2;
    }
    if (neighborConjugated && AMIDE_LIKE_ELEMENTS[atom.element] === true) {
        return SP2;
    }
    // 带负电的氧/氮 (如羧酸根) 也按 sp2 处理
    if (atom.charge < 0 && (atom.element === 'O' || atom.element === 'N' || atom.element === 'S')) {
        return SP2;
    }
    return SP3;
}

/** 小环的环内角。环太大时返回 null (用杂化角) */
function ringAngleFor(size: number, aromaticLike: boolean): number | null {
    if (size === 3) {
        return 60.0;
    }
    if (size === 4) {
        return 90.0;
    }
    if (size === 5) {
        return aromaticLike ? AROMATIC_RING_ANGLE[5] : 104.0;
    }
    return null;
}

/**
 * 原子 center 上 nb1-center-nb2 的理论键角 (度)。
 */
function bondAngleAt(mol: Mol, center: number, nb1: number, nb2: number,
                     ringEdge: { [key: string]: number }, ringSize: number[], hybrid: number): number {
    const r1 = ringEdge[ringKey(center, nb1)] === undefined ? -1 : ringEdge[ringKey(center, nb1)];
    const r2 = ringEdge[ringKey(center, nb2)] === undefined ? -1 : ringEdge[ringKey(center, nb2)];
    if (r1 >= 0 && r1 === r2 && r1 < ringSize.length) {
        const size = ringSize[r1];
        const atoms = mol.atoms;
        let aromaticLike = atoms[center].aromatic;
        if (!aromaticLike) {
            aromaticLike = hybrid === SP2;
        }
        const angle = ringAngleFor(size, aromaticLike);
        if (angle !== null) {
            return angle;
        }
    }
    const table = bondAnglesFor(mol.atoms[center].element);
    return table[hybrid];
}

/** 返回原子 index 的全部邻居原子索引。 */
function neighborAtoms(mol: Mol, index: number): number[] {
    const out: number[] = [];
    const list = mol.adj[index];
    for (let t = 0; t < list.length; t++) {
        const bond = mol.bonds[list[t]];
        out.push(bond.b === index ? bond.a : bond.b);
    }
    return out;
}

/**
 * 为目标扭转角建立查表。
 * 返回 { 键索引: { "a端取代基|b端取代基": 目标二面角(度) } }
 */
function buildTorsionTargets(mol: Mol, rings: number[][], ringEdge: { [key: string]: number }): { [key: number]: { [key: string]: number } } {
    const targets: { [key: number]: { [key: string]: number } } = {};
    const atoms = mol.atoms;
    const hybridCache: number[] = [];
    for (let i = 0; i < atoms.length; i++) {
        hybridCache[i] = atomHybridization(mol, i);
    }

    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        if (bond.aromatic) {
            continue;
        }
        const a = bond.a;
        const b = bond.b;
        // sp 原子 (累积双键中心、腈基等) 的取代基与键轴共线, 二面角无定义。
        // 若强行按"默认反式"给约束, 求解器会把数值噪声当成偏差, 把取代基
        // 甩得满处乱转 —— 实测丙二烯会因此被反复扭弯。因此这类键一律跳过。
        if (hybridCache[a] === SP || hybridCache[b] === SP) {
            continue;
        }
        if (ringEdge[ringKey(a, b)] !== undefined) {
            continue;                      // 环内键交给环
        }
        const subsA = neighborAtoms(mol, a);
        const subsB = neighborAtoms(mol, b);
        for (let x = subsA.length - 1; x >= 0; x--) {
            if (subsA[x] === b) {
                subsA.splice(x, 1);
            }
        }
        for (let y = subsB.length - 1; y >= 0; y--) {
            if (subsB[y] === a) {
                subsB.splice(y, 1);
            }
        }
        if (subsA.length === 0 || subsB.length === 0) {
            continue;
        }

        // 参照取代基 refA/refB: 双键时优先取带 / \ 方向记号者, 否则取第一个;
        // 表格循环里还要用它们判定"同侧/异侧", 因此声明在 if 之外。
        let base: number;
        let refA = -1;
        let refB = -1;
        if (bond.order === 2) {
            for (let t = 0; t < subsA.length; t++) {
                if (bondDirection(mol, a, subsA[t])) {
                    refA = subsA[t];
                    break;
                }
            }
            for (let t = 0; t < subsB.length; t++) {
                if (bondDirection(mol, b, subsB[t])) {
                    refB = subsB[t];
                    break;
                }
            }
            if (refA < 0) {
                refA = subsA[0];
            }
            if (refB < 0) {
                refB = subsB[0];
            }
            const dirA = bondDirection(mol, a, refA);
            const dirB = bondDirection(mol, b, refB);
            if (dirA && dirB) {
                base = dirA === dirB ? 180.0 : 0.0;
            } else {
                base = 180.0;              // 未说明时默认反式
            }
        } else if (hybridCache[a] === SP2 && hybridCache[b] === SP2) {
            refA = subsA[0];
            refB = subsB[0];
            base = 180.0;
        } else {
            continue;
        }

        const table: { [key: string]: number } = {};
        for (let t1 = 0; t1 < subsA.length; t1++) {
            const x = subsA[t1];
            for (let t2 = 0; t2 < subsB.length; t2++) {
                const y = subsB[t2];
                const sameSide = ((x === refA) === (y === refB));
                if (sameSide) {
                    table[x + '|' + y] = base;
                } else {
                    table[x + '|' + y] = Math.abs(base - 180.0) < 1e-6 ? 0.0 : 180.0;
                }
            }
        }
        targets[k] = table;
    }

    // --- 全 sp2 环的平面性 ---
    for (let ri = 0; ri < rings.length; ri++) {
        const ring = rings[ri];
        const count = ring.length;
        if (count < 4) {
            continue;
        }
        let planar = true;
        for (let t = 0; t < ring.length; t++) {
            if (hybridCache[ring[t]] !== SP2) {
                planar = false;
                break;
            }
        }
        if (!planar) {
            continue;
        }
        for (let position = 0; position < count; position++) {
            const b = ring[position];
            const c = ring[(position + 1) % count];
            const a = ring[(position - 1 + count) % count];
            const d = ring[(position + 2) % count];
            if (a === c || d === b || a === d) {
                continue;
            }
            const k = findBond(mol, b, c);
            if (k < 0) {
                continue;
            }
            let table = targets[k];
            if (table === undefined) {
                table = {};
                targets[k] = table;
            }
            if (mol.bonds[k].a === b) {
                table[a + '|' + d] = 0.0;
            } else {
                table[d + '|' + a] = 0.0;
            }
        }
    }
    return targets;
}

/** 返回 a-b 键上的方向记号 ('/', '\\' 或 '') */
function bondDirection(mol: Mol, a: number, b: number): string {
    const list = mol.adj[a];
    for (let t = 0; t < list.length; t++) {
        const bond = mol.bonds[list[t]];
        if ((bond.a === a && bond.b === b) || (bond.a === b && bond.b === a)) {
            return bond.direction;
        }
    }
    return '';
}

/** 返回连接 a-b 的键索引, 不存在返回 -1 */
function findBond(mol: Mol, a: number, b: number): number {
    const list = mol.adj[a];
    for (let t = 0; t < list.length; t++) {
        const bond = mol.bonds[list[t]];
        if ((bond.a === a && bond.b === b) || (bond.a === b && bond.b === a)) {
            return list[t];
        }
    }
    return -1;
}

/**
 * 查询 dihedral(p0,p1,p2,p3) 的目标值; 无目标返回 null。
 * 中心键为 p1-p2。
 */
function lookupTorsionTarget(mol: Mol, targets: { [key: number]: { [key: string]: number } },
                             p0: number, p1: number, p2: number, p3: number): number | null {
    const k = findBond(mol, p1, p2);
    if (k < 0) {
        return null;
    }
    const table = targets[k];
    if (table === undefined || Object.keys(table).length === 0) {
        return null;
    }
    const bond = mol.bonds[k];
    const key = bond.a === p1 ? (p0 + '|' + p3) : (p3 + '|' + p0);
    const v = table[key];
    return v === undefined ? null : v;
}

// =====================================================================
// 8. 三维坐标生成
// =====================================================================

/**
 * 均匀网格空间索引。用于把 O(n^2) 的近邻搜索降到 O(n)。
 * 只增不改: 摆放阶段原子逐个加入, 结构松弛阶段每轮重建。
 */
class SpatialGrid {
    private cell: number;
    private buckets: { [key: string]: number[] };

    constructor(cell: number) {
        this.cell = cell;
        this.buckets = {};
    }

    /** 把坐标映射为网格单元键 floor(x/cell), floor(y/cell), floor(z/cell)。 */
    keyOf(p: Vec3): string {
        const c = this.cell;
        return Math.floor(p[0] / c) + ',' + Math.floor(p[1] / c) + ',' + Math.floor(p[2] / c);
    }

    /** 把原子索引加入它所在的那个单元。 */
    add(index: number, p: Vec3): void {
        const k = this.keyOf(p);
        const bucket = this.buckets[k];
        if (bucket === undefined) {
            this.buckets[k] = [index];
        } else {
            bucket.push(index);
        }
    }

    /** 返回与 p 所在格相距 radius 格以内的所有已加入原子索引 */
    query(p: Vec3, radius?: number): number[] {
        const r = radius === undefined ? 1 : radius;
        const parts = this.keyOf(p).split(',');
        const kx = parseInt(parts[0], 10);
        const ky = parseInt(parts[1], 10);
        const kz = parseInt(parts[2], 10);
        const out: number[] = [];
        for (let dx = -r; dx <= r; dx++) {
            for (let dy = -r; dy <= r; dy++) {
                for (let dz = -r; dz <= r; dz++) {
                    const bucket = this.buckets[(kx + dx) + ',' + (ky + dy) + ',' + (kz + dz)];
                    if (bucket !== undefined) {
                        for (let t = 0; t < bucket.length; t++) {
                            out.push(bucket[t]);
                        }
                    }
                }
            }
        }
        return out;
    }
}

/** 把无序原子对压成一个整数键 (i<j) */
function pairKey(i: number, j: number, n: number): number {
    if (i > j) {
        const t = i;
        i = j;
        j = t;
    }
    return i * n + j;
}

/**
 * 放置新原子 X, 精确满足三个几何量:
 *   |X - origin| = length, 角 X-origin-ref1 = angleDeg,
 *   二面角 X-origin-ref1-ref2 = dihedralDeg
 * ref2 为 null 时, 参考平面任意选取。
 */
function placeAtom(origin: Vec3, ref1: Vec3, ref2: Vec3 | null, length: number,
                   angleDeg: number, dihedralDeg: number): Vec3 {
    const axis = vNorm(vSub(ref1, origin));
    let u: Vec3;
    if (ref2 === null) {
        const aux: Vec3 = Math.abs(axis[2]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0];
        u = vNorm(vCross(axis, aux));
    } else {
        const d = vSub(ref2, ref1);
        const perp = vSub(d, vScale(axis, vDot(d, axis)));
        if (vLen(perp) < 1e-8) {
            const aux: Vec3 = Math.abs(axis[2]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0];
            u = vNorm(vCross(axis, aux));
        } else {
            u = vNorm(perp);
        }
    }
    const w = vCross(axis, u);
    const rad = angleDeg * Math.PI / 180.0;
    const phi = dihedralDeg * Math.PI / 180.0;
    const perpDir = vAdd(vScale(u, Math.cos(phi)), vScale(w, Math.sin(phi)));
    const direction = vAdd(vScale(axis, Math.cos(rad)), vScale(perpDir, Math.sin(rad)));
    return vAdd(origin, vScale(direction, length));
}

/** 返回连通片段列表, 每个片段是原子索引列表 (按出现顺序) */
function findFragments(mol: Mol): number[][] {
    const n = mol.atoms.length;
    const seen: boolean[] = new Array<boolean>(n);
    for (let i = 0; i < n; i++) {
        seen[i] = false;
    }
    const fragments: number[][] = [];
    for (let start = 0; start < n; start++) {
        if (seen[start]) {
            continue;
        }
        seen[start] = true;
        const queue: number[] = [start];
        let head = 0;
        const group: number[] = [start];
        while (head < queue.length) {
            const cur = queue[head];
            head += 1;
            const list = mol.adj[cur];
            for (let t = 0; t < list.length; t++) {
                const bond = mol.bonds[list[t]];
                const nxt = bond.b === cur ? bond.a : bond.b;
                if (!seen[nxt]) {
                    seen[nxt] = true;
                    group.push(nxt);
                    queue.push(nxt);
                }
            }
        }
        fragments.push(group);
    }
    return fragments;
}

/** 返回与 start 拓扑距离 <= depth 的原子集合 (含 start) */
function nearAtoms(mol: Mol, start: number, depth: number): { [key: number]: boolean } {
    const seen: { [key: number]: boolean } = {};
    seen[start] = true;
    let frontier: number[] = [start];
    for (let d = 0; d < depth; d++) {
        const nxt: number[] = [];
        for (let f = 0; f < frontier.length; f++) {
            const cur = frontier[f];
            const list = mol.adj[cur];
            for (let t = 0; t < list.length; t++) {
                const bond = mol.bonds[list[t]];
                const other = bond.b === cur ? bond.a : bond.b;
                if (seen[other] !== true) {
                    seen[other] = true;
                    nxt.push(other);
                }
            }
        }
        frontier = nxt;
        if (frontier.length === 0) {
            break;
        }
    }
    return seen;
}

interface PlacementCtx {
    mol: Mol;
    elements: string[];
    radii: number[];
    grid: SpatialGrid;
}

/**
 * 候选位置的代价 = 成环键长度偏差 + 严重重叠 + 非键排斥。
 */
function placementEnergy(ctx: PlacementCtx, pos: Vec3, curr: number, newIndex: number,
                         coords: Vec3[], placed: boolean[], near: { [key: number]: boolean }): number {
    let energy = 0.0;
    const elements = ctx.elements;
    const radii = ctx.radii;
    const mol = ctx.mol;
    const rNew = radii[newIndex];
    const list = mol.adj[newIndex];
    for (let t = 0; t < list.length; t++) {
        const bond = mol.bonds[list[t]];
        const other = bond.b === newIndex ? bond.a : bond.b;
        if (other === curr || other === newIndex || !placed[other]) {
            continue;
        }
        const target = bondLength(elements[newIndex], elements[other], bond.order);
        const dd = vDist(pos, coords[other]);
        energy += 20.0 * (dd - target) * (dd - target);
    }
    const candidates = ctx.grid.query(pos);
    for (let t = 0; t < candidates.length; t++) {
        const idx = candidates[t];
        if ((!placed[idx]) || idx === newIndex) {
            continue;
        }
        const dd = vDist(pos, coords[idx]);
        if (dd < 1e-6) {
            energy += 1000.0;
            continue;
        }
        const radiiSum = rNew + radii[idx];
        if (dd < 0.5 * radiiSum) {
            energy += 100.0;
            continue;
        }
        if (near[idx] === true) {
            continue;
        }
        const contact = 0.78 * radiiSum;
        if (dd < contact) {
            const ratio = contact / dd;
            const ratio2 = ratio * ratio;
            energy += ratio2 * ratio2 * ratio2;
        }
    }
    return energy;
}

/**
 * 扫描二面角, 取代价最小者。先粗扫再在最优值附近细化。
 */
function scanDihedral(ctx: PlacementCtx, curr: number, newIndex: number, ref1: number,
                      ref2: Vec3 | null, length: number, angle: number, coords: Vec3[],
                      placed: boolean[], near: { [key: number]: boolean }, step: number,
                      target: number | null): number {
    let bestAngle = 0.0;
    let bestEnergy: number | null = null;
    let dihedral = 0.0;
    while (dihedral < 360.0) {
        const pos = placeAtom(coords[curr], coords[ref1], ref2, length, angle, dihedral);
        let energy = placementEnergy(ctx, pos, curr, newIndex, coords, placed, near);
        if (target !== null) {
            const delta = modPositive(dihedral - target + 180.0, 360.0) - 180.0;
            energy += 1.5 * (delta / 180.0) * (delta / 180.0);
        }
        if (bestEnergy === null || energy < bestEnergy) {
            bestEnergy = energy;
            bestAngle = dihedral;
        }
        dihedral += step;
    }
    const base = bestAngle;
    dihedral = base - step;
    while (dihedral <= base + step + 1e-9) {
        const pos = placeAtom(coords[curr], coords[ref1], ref2, length, angle, dihedral);
        let energy = placementEnergy(ctx, pos, curr, newIndex, coords, placed, near);
        if (target !== null) {
            const delta = modPositive(dihedral - target + 180.0, 360.0) - 180.0;
            energy += 1.5 * (delta / 180.0) * (delta / 180.0);
        }
        if (energy < (bestEnergy as number)) {
            bestEnergy = energy;
            bestAngle = dihedral;
        }
        dihedral += step / 3.0;
    }
    return bestAngle;
}

/**
 * 生成初始三维坐标。
 */
function generateCoordinates(mol: Mol, targets: { [key: number]: { [key: string]: number } },
                             ringEdge: { [key: string]: number }, ringSize: number[],
                             opts: SmilesOptions): Vec3[] {
    const n = mol.atoms.length;
    const ctx: PlacementCtx = {
        mol: mol,
        elements: [],
        radii: [],
        grid: new SpatialGrid(5.0),
    };
    for (let i = 0; i < n; i++) {
        ctx.elements.push(mol.atoms[i].element);
        ctx.radii.push(vdwRadius(mol.atoms[i].element));
    }
    const elements = ctx.elements;
    const coords: Vec3[] = [];
    for (let i = 0; i < n; i++) {
        coords.push([0.0, 0.0, 0.0]);
    }
    const placed: boolean[] = new Array<boolean>(n);
    const parent: number[] = new Array<number>(n);
    const hybridization: number[] = new Array<number>(n);
    for (let i = 0; i < n; i++) {
        placed[i] = false;
        parent[i] = -1;
        hybridization[i] = atomHybridization(mol, i);
    }

    let step = opts.dihedral_step;
    if (step <= 0.0) {
        step = 15.0;
    }

    const fragments = findFragments(mol);
    let cursorX = 0.0;                 // 已排布片段占用的 +x 边界

    for (let gi = 0; gi < fragments.length; gi++) {
        const group = fragments[gi];
        // 每个片段单独用一张空间网格, 因为片段整体会被平移, 旧网格会失效
        const grid = new SpatialGrid(5.0);
        ctx.grid = grid;
        const heavy: number[] = [];
        let hydrogens: number[] = [];
        for (let t = 0; t < group.length; t++) {
            const i = group[t];
            if (elements[i] === 'H') {
                hydrogens.push(i);
            } else {
                heavy.push(i);
            }
        }
        if (heavy.length === 0) {
            for (let t = 0; t < hydrogens.length; t++) {
                heavy.push(hydrogens[t]);
            }
            hydrogens = [];
        }

        let root = heavy[0];
        let bestDegree = -1;
        for (let t = 0; t < heavy.length; t++) {
            const i = heavy[t];
            const degree = mol.adj[i].length;
            if (degree > bestDegree) {
                bestDegree = degree;
                root = i;
            }
        }

        // 重原子生成树: 用深度优先而非广度优先 —— 广度优先会在环上留下很长的缺口,
        // 把最后几个原子摆到环的另一侧 (实测苯环会出现 5 Å 的"键"); 深度优先沿环一路
        // 走到底, 只留一个缺口, 而缺口两端此时都已摆好, 由成环键约束闭合成正多边形。
        const order: number[] = [root];
        const seen: { [key: number]: boolean } = {};
        seen[root] = true;
        const stack: number[] = [root];
        while (stack.length > 0) {
            const cur = stack.pop() as number;
            const children: number[] = [];
            const list = mol.adj[cur];
            for (let t = 0; t < list.length; t++) {
                const bond = mol.bonds[list[t]];
                const nxt = bond.b === cur ? bond.a : bond.b;
                if (seen[nxt] === true || elements[nxt] === 'H') {
                    continue;
                }
                seen[nxt] = true;
                parent[nxt] = cur;
                children.push(nxt);
            }
            for (let t = 0; t < children.length; t++) {
                order.push(children[t]);
            }
            let position = children.length - 1;
            while (position >= 0) {
                stack.push(children[position]);
                position -= 1;
            }
        }
        // 追加氢 (以及任何漏掉的原子)
        for (let t = 0; t < group.length; t++) {
            const i = group[t];
            if (seen[i] !== true) {
                parent[i] = -1;
                order.push(i);
                seen[i] = true;
            }
        }

        const start = order[0];
        coords[start] = [0.0, 0.0, 0.0];
        placed[start] = true;
        grid.add(start, coords[start]);

        for (let oi = 1; oi < order.length; oi++) {
            const index = order[oi];
            const placedNeighbors: number[] = [];
            const list = mol.adj[index];
            for (let t = 0; t < list.length; t++) {
                const bond = mol.bonds[list[t]];
                const other = bond.b === index ? bond.a : bond.b;
                if (placed[other]) {
                    placedNeighbors.push(other);
                }
            }

            let anchor = parent[index];
            if (anchor < 0 || !placed[anchor]) {
                anchor = placedNeighbors.length > 0 ? placedNeighbors[0] : -1;
            }
            if (anchor < 0) {
                // 兜底: 与已放置部分不连通, 单独摆一个位置
                coords[index] = [1.5 + 1.5 * order.length, 0.0, 0.0];
                placed[index] = true;
                grid.add(index, coords[index]);
                continue;
            }

            const kBond = findBond(mol, anchor, index);
            let length: number;
            if (kBond < 0) {
                length = 1.5;
            } else {
                length = bondLength(elements[anchor], elements[index], mol.bonds[kBond].order);
            }

            // 选参考原子: 旋转轴必须过 anchor 并指向另一个已摆放的原子, 故取父原子作
            // ref1 (旋转轴方向), 祖父作 ref2 (二面角零点); anchor 是片段根时改用 anchor
            // 其它已摆放的邻居。ref2 仍取不到时, 再找 ref1 上另一个已摆放的邻居 —— 顺反
            // 异构的目标角正是靠这个参考点才能查到。
            let ref1 = -1;
            let ref2 = -1;
            const grand = parent[anchor];
            if (grand >= 0 && placed[grand] && grand !== index) {
                ref1 = grand;
                const great = parent[grand];
                if (great >= 0 && placed[great] && great !== index) {
                    ref2 = great;
                }
            }
            if (ref1 < 0) {
                // anchor 是片段的根: 用 anchor 其它已放置的邻居作参考
                const anchorPlaced: number[] = [];
                const listA = mol.adj[anchor];
                for (let t = 0; t < listA.length; t++) {
                    const bond = mol.bonds[listA[t]];
                    const other = bond.b === anchor ? bond.a : bond.b;
                    if (other !== index && placed[other]) {
                        anchorPlaced.push(other);
                    }
                }
                if (anchorPlaced.length > 0) {
                    ref1 = anchorPlaced[0];
                }
                if (anchorPlaced.length > 1) {
                    ref2 = anchorPlaced[1];
                }
            }

            // ref2 兜底: 若取不到曾祖父, 就用 ref1 上另一个已放置的邻居。
            if (ref1 >= 0 && ref2 < 0) {
                const listR = mol.adj[ref1];
                for (let t = 0; t < listR.length; t++) {
                    const bond = mol.bonds[listR[t]];
                    const other = bond.b === ref1 ? bond.a : bond.b;
                    if (other !== anchor && other !== index && placed[other]) {
                        ref2 = other;
                        break;
                    }
                }
            }

            if (ref1 < 0) {
                coords[index] = [coords[anchor][0] + length, coords[anchor][1], coords[anchor][2]];
                placed[index] = true;
                grid.add(index, coords[index]);
                continue;
            }

            const angle = bondAngleAt(mol, anchor, index, ref1, ringEdge, ringSize, hybridization[anchor]);
            const ref2Coord: Vec3 | null = ref2 >= 0 ? coords[ref2] : null;

            let target: number | null = null;
            if (ref2 >= 0) {
                target = lookupTorsionTarget(mol, targets, index, anchor, ref1, ref2);
            }

            // 若该原子还有"已放置的成环键伙伴", 成环优先: 改为扫描。
            let ringPartnerPlaced = false;
            const listI = mol.adj[index];
            for (let t = 0; t < listI.length; t++) {
                const bond = mol.bonds[listI[t]];
                const other = bond.b === index ? bond.a : bond.b;
                if (other !== anchor && other !== index && placed[other]) {
                    ringPartnerPlaced = true;
                    break;
                }
            }

            let dihedral: number;
            if (target !== null && !ringPartnerPlaced) {
                dihedral = target;
            } else {
                // 只把 1-2、1-3 排除在立体排斥之外; 1-4 必须计入。
                const near = nearAtoms(mol, index, 2);
                dihedral = scanDihedral(ctx, anchor, index, ref1, ref2Coord,
                                        length, angle, coords, placed, near,
                                        step, target);
            }

            coords[index] = placeAtom(coords[anchor], coords[ref1], ref2Coord, length, angle, dihedral);
            placed[index] = true;
            grid.add(index, coords[index]);
        }

        // 把本片段整体平移: 令其 x 最小端正好落在 cursor_x 上, 再留出间隙。
        let minX: number | null = null;
        let maxX: number | null = null;
        for (let t = 0; t < group.length; t++) {
            const x = coords[group[t]][0];
            if (minX === null || x < minX) {
                minX = x;
            }
            if (maxX === null || x > maxX) {
                maxX = x;
            }
        }
        if (minX === null) {
            continue;
        }
        const shift = cursorX - minX;
        if (shift !== 0.0) {
            for (let t = 0; t < group.length; t++) {
                coords[group[t]][0] += shift;
            }
        }
        cursorX = (maxX as number) + shift + opts.fragment_gap;
    }

    for (let i = 0; i < n; i++) {
        if (!placed[i]) {
            coords[i] = [cursorX, 0.0, 0.0];
            cursorX += opts.fragment_gap;
        }
    }
    return coords;
}

// =====================================================================
// 9. 结构松弛
// =====================================================================
//
// 所有几何要求都统一表达为"两点距离约束", 因此只需要一种投影算子:
//   键长 (1-2) -> 目标 = 标准键长
//   键角 (1-3) -> 目标 = 由键长与理论键角算出的第三边距离
//   非键(>1-4)-> 下限 = 范德华接触距离 (只推不拉)
// 扭转角不在这里: 它不能用"1-4 距离"来表达 (平面处该距离对二面角不敏感),
// 由 enforceDihedralTargets 用整支刚性旋转精确求解。

/** 返回 a-b 两原子之间的键级; 无键时按单键处理。 */
function bondOrderOf(mol: Mol, a: number, b: number): number {
    const k = findBond(mol, a, b);
    if (k < 0) {
        return 1;
    }
    return mol.bonds[k].order;
}

/**
 * 构建"距离型"约束 (键长与键角), 返回 [i, j, 目标距离, 刚度] 或
 * [i, j, 0, 刚度, 1, 中心原子] (直线型中心的实时目标)。
 */
function buildConstraints(mol: Mol, targets: { [key: number]: { [key: string]: number } },
                          ringEdge: { [key: string]: number }, ringSize: number[]): number[][] {
    void targets;
    const elements: string[] = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        elements.push(mol.atoms[i].element);
    }
    const n = mol.atoms.length;
    const constraints: number[][] = [];
    const hybridization: number[] = new Array<number>(n);
    for (let i = 0; i < n; i++) {
        hybridization[i] = atomHybridization(mol, i);
    }

    // --- 1-2 键长 ---
    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        const a = bond.a;
        const b = bond.b;
        const target = bondLength(elements[a], elements[b], bond.order);
        constraints.push([a, b, target, 1.0]);
    }

    // --- 1-3 键角 ---
    // 对直线型中心 (sp 且只有两个邻居) 用"实时目标": 目标距离取当前两根键的
    // 实测长度之和。180 度附近角度对距离极其敏感: 目标距离差 0.002 Å 就会让键角差 8 度,
    // 若用理想键长算出固定目标, 实测丙二烯只能到 168 度; 改为实时目标后可精确到 180 度。
    const neighborLists: number[][] = [];
    for (let i = 0; i < n; i++) {
        neighborLists.push(neighborAtoms(mol, i));
    }
    for (let c = 0; c < n; c++) {
        const neighbors = neighborLists[c];
        const degree = neighbors.length;
        if (degree < 2) {
            continue;
        }
        let stiffness = degree <= 2 ? 0.75 : 0.45;
        const linear = (degree === 2 && hybridization[c] === SP);
        if (linear) {
            stiffness = 1.0;
        }
        for (let x = 0; x < degree; x++) {
            const a = neighbors[x];
            for (let y = x + 1; y < degree; y++) {
                const b = neighbors[y];
                if (linear) {
                    constraints.push([a, b, 0.0, stiffness, 1, c]);
                    continue;
                }
                const angle = bondAngleAt(mol, c, a, b, ringEdge, ringSize, hybridization[c]);
                const la = bondLength(elements[c], elements[a], bondOrderOf(mol, c, a));
                const lb = bondLength(elements[c], elements[b], bondOrderOf(mol, c, b));
                const cosv = Math.cos(angle * Math.PI / 180.0);
                const target = Math.sqrt(la * la + lb * lb - 2.0 * la * lb * cosv);
                constraints.push([a, b, target, stiffness]);
            }
        }
    }
    return constraints;
}

/**
 * 挑出必须保持平面的环: 环上所有原子都是 sp2。
 */
function buildPlanarRings(mol: Mol, rings: number[][]): number[][] {
    const planar: number[][] = [];
    for (let ri = 0; ri < rings.length; ri++) {
        const ring = rings[ri];
        if (ring.length < 4) {
            continue;
        }
        let flat = true;
        for (let t = 0; t < ring.length; t++) {
            if (atomHybridization(mol, ring[t]) !== SP2) {
                flat = false;
                break;
            }
        }
        if (flat) {
            planar.push(ring);
        }
    }
    return planar;
}

/**
 * 把环上的原子投影到最佳拟合平面, 使其变平 (Newell 法求法向)。
 * 返回最大位移量。
 */
function projectRingToPlane(coords: Vec3[], ring: number[]): number {
    const count = ring.length;
    let nx = 0.0;
    let ny = 0.0;
    let nz = 0.0;
    for (let i = 0; i < count; i++) {
        const p = coords[ring[i]];
        const q = coords[ring[(i + 1) % count]];
        nx += (p[1] - q[1]) * (p[2] + q[2]);
        ny += (p[2] - q[2]) * (p[0] + q[0]);
        nz += (p[0] - q[0]) * (p[1] + q[1]);
    }
    const length = Math.sqrt(nx * nx + ny * ny + nz * nz);
    if (length < 1e-9) {
        return 0.0;
    }
    const unit = [nx / length, ny / length, nz / length];
    let cx = 0.0;
    let cy = 0.0;
    let cz = 0.0;
    for (let t = 0; t < ring.length; t++) {
        cx += coords[ring[t]][0];
        cy += coords[ring[t]][1];
        cz += coords[ring[t]][2];
    }
    cx /= count;
    cy /= count;
    cz /= count;
    let worst = 0.0;
    for (let t = 0; t < ring.length; t++) {
        const p = coords[ring[t]];
        const offset = vDot(vSub(p, [cx, cy, cz]), unit);
        p[0] -= unit[0] * offset;
        p[1] -= unit[1] * offset;
        p[2] -= unit[2] * offset;
        if (Math.abs(offset) > worst) {
            worst = Math.abs(offset);
        }
    }
    return worst;
}

/**
 * 挑出需要保持局部平面的 sp2 中心: 恰好 3 个邻居的 sp2 原子 (跳过 sp2 环上的原子)。
 */
function buildSp2Centers(mol: Mol, planarRings: number[][]): number[][] {
    const inPlanarRing: { [key: number]: boolean } = {};
    for (let ri = 0; ri < planarRings.length; ri++) {
        const ring = planarRings[ri];
        for (let t = 0; t < ring.length; t++) {
            inPlanarRing[ring[t]] = true;
        }
    }
    const centers: number[][] = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        if (inPlanarRing[i] === true) {
            continue;
        }
        const neighbors = neighborAtoms(mol, i);
        if (neighbors.length !== 3) {
            continue;
        }
        if (atomHybridization(mol, i) !== SP2) {
            continue;
        }
        centers.push([i, neighbors[0], neighbors[1], neighbors[2]]);
    }
    return centers;
}

/** 把 center 投影到 n1、n2、n3 所在平面上, 返回位移量 (Å) */
function projectAtomToPlane(coords: Vec3[], center: number, n1: number, n2: number, n3: number): number {
    const p1 = coords[n1];
    const p2 = coords[n2];
    const p3 = coords[n3];
    const normal = vCross(vSub(p2, p1), vSub(p3, p1));
    const length = vLen(normal);
    if (length < 1e-9) {
        return 0.0;
    }
    const unit = [normal[0] / length, normal[1] / length, normal[2] / length];
    const p = coords[center];
    const offset = vDot(vSub(p, p1), unit);
    p[0] -= unit[0] * offset;
    p[1] -= unit[1] * offset;
    p[2] -= unit[2] * offset;
    return Math.abs(offset);
}

interface Exclusions {
    excluded: { [key: number]: boolean };
    tight: { [key: number]: boolean };
}

/**
 * 构建非键作用的排除表: excluded (拓扑距离 <= 3), tight (拓扑距离 == 4)。
 */
function buildExclusions(mol: Mol): Exclusions {
    const n = mol.atoms.length;
    const excluded: { [key: number]: boolean } = {};
    const tight: { [key: number]: boolean } = {};
    for (let i = 0; i < n; i++) {
        const adj1 = mol.adj[i];
        for (let a1 = 0; a1 < adj1.length; a1++) {
            const bond1 = mol.bonds[adj1[a1]];
            const j = bond1.b === i ? bond1.a : bond1.b;
            excluded[pairKey(i, j, n)] = true;
            const adj2 = mol.adj[j];
            for (let a2 = 0; a2 < adj2.length; a2++) {
                const bond2 = mol.bonds[adj2[a2]];
                const m = bond2.b === j ? bond2.a : bond2.b;
                if (m === i) {
                    continue;
                }
                excluded[pairKey(i, m, n)] = true;
                const adj3 = mol.adj[m];
                for (let a3 = 0; a3 < adj3.length; a3++) {
                    const bond3 = mol.bonds[adj3[a3]];
                    const p = bond3.b === m ? bond3.a : bond3.b;
                    if (p === j || p === i) {
                        continue;
                    }
                    excluded[pairKey(i, p, n)] = true;
                    const adj4 = mol.adj[p];
                    for (let a4 = 0; a4 < adj4.length; a4++) {
                        const bond4 = mol.bonds[adj4[a4]];
                        const q = bond4.b === p ? bond4.a : bond4.b;
                        if (q === m || q === j || q === i) {
                            continue;
                        }
                        const key = pairKey(i, q, n);
                        if (excluded[key] !== true) {
                            tight[key] = true;
                        }
                    }
                }
            }
        }
    }
    return { excluded: excluded, tight: tight };
}

/**
 * 按当前位置收集可能发生范德华接触的原子对 (Verlet 邻居表)。
 * 每项为 [i, j, 接触系数]。
 */
function buildPairList(mol: Mol, coords: Vec3[], excluded: { [key: number]: boolean },
                       tight: { [key: number]: boolean }, cutoff: number): number[][] {
    const n = mol.atoms.length;
    const grid = new SpatialGrid(cutoff);
    for (let i = 0; i < n; i++) {
        grid.add(i, coords[i]);
    }
    const pairs: number[][] = [];
    const seen: { [key: number]: boolean } = {};
    for (let i = 0; i < n; i++) {
        const candidates = grid.query(coords[i]);
        for (let t = 0; t < candidates.length; t++) {
            const j = candidates[t];
            if (j === i) {
                continue;
            }
            const key = pairKey(i, j, n);
            if (seen[key] === true || excluded[key] === true) {
                continue;
            }
            seen[key] = true;
            if (tight[key] === true) {
                pairs.push([i, j, 0.62]);
            } else {
                pairs.push([i, j, 0.80]);
            }
        }
    }
    return pairs;
}

/**
 * 迭代投影求解全部约束。返回实际迭代轮数。
 */
function relax(mol: Mol, coords: Vec3[], constraints: number[][], planarRings: number[][],
               sp2Centers: number[][], excluded: { [key: number]: boolean },
               tight: { [key: number]: boolean }, opts: SmilesOptions): number {
    const n = coords.length;
    if (n < 2) {
        return 0;
    }
    const radii: number[] = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        radii.push(vdwRadius(mol.atoms[i].element));
    }
    const contactWeight = opts.contact_weight;
    let maxIterations = opts.relax_iterations;
    if (n > 200) {
        // 大分子每轮代价高, 而局部几何早已收敛, 因此限制轮数控制耗时
        if (maxIterations > 250) {
            maxIterations = 250;
        }
    }
    const tolerance = opts.relax_tolerance;
    const cutoff = 5.0;

    let pairs = buildPairList(mol, coords, excluded, tight, cutoff);
    const rebuildEvery = 25;
    let iterations = 0;

    for (let iteration = 0; iteration < maxIterations; iteration++) {
        iterations = iteration + 1;
        let maxMove = 0.0;

        // 距离约束 (高斯-赛德尔逐个投影)
        for (let ci = 0; ci < constraints.length; ci++) {
            const item = constraints[ci];
            const i = item[0];
            const j = item[1];
            let target = item[2];
            const stiffness = item[3];
            if (item.length > 4) {
                // 实时目标 (直线型中心): 目标 = 中心到两端的实测键长之和
                const center = item[5];
                target = vDist(coords[i], coords[center]) + vDist(coords[center], coords[j]);
            }
            const pi = coords[i];
            const pj = coords[j];
            const dx = pj[0] - pi[0];
            const dy = pj[1] - pi[1];
            const dz = pj[2] - pi[2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (dist < 1e-9) {
                continue;
            }
            // 限制单次修正幅度: 两点几乎重合时未加限幅会造成周期性震荡而无法收敛。
            let gap = dist - target;
            if (gap > MAX_CONSTRAINT_STEP) {
                gap = MAX_CONSTRAINT_STEP;
            } else if (gap < -MAX_CONSTRAINT_STEP) {
                gap = -MAX_CONSTRAINT_STEP;
            }
            const correction = gap / dist * RELAX_FACTOR * stiffness;
            const mx = dx * correction;
            const my = dy * correction;
            const mz = dz * correction;
            pi[0] += mx;
            pi[1] += my;
            pi[2] += mz;
            pj[0] -= mx;
            pj[1] -= my;
            pj[2] -= mz;
            const move = Math.abs(mx) + Math.abs(my) + Math.abs(mz);
            if (move > maxMove) {
                maxMove = move;
            }
        }

        // 环平面性: 把 sp2 环投影到最佳拟合平面。
        for (let ri = 0; ri < planarRings.length; ri++) {
            const moved = projectRingToPlane(coords, planarRings[ri]);
            if (moved > maxMove) {
                maxMove = moved;
            }
        }

        // 局部 sp2 平面性: 把 sp2 中心压回三个邻居所在平面
        for (let si = 0; si < sp2Centers.length; si++) {
            const item = sp2Centers[si];
            const moved = projectAtomToPlane(coords, item[0], item[1], item[2], item[3]);
            if (moved > maxMove) {
                maxMove = moved;
            }
        }

        // 非键排斥 (只推不拉, 避免分子塌缩)
        for (let pi2 = 0; pi2 < pairs.length; pi2++) {
            const pair = pairs[pi2];
            const i = pair[0];
            const j = pair[1];
            const contact = pair[2] * (radii[i] + radii[j]);
            const pi = coords[i];
            const pj = coords[j];
            const dx = pj[0] - pi[0];
            const dy = pj[1] - pi[1];
            const dz = pj[2] - pi[2];
            const dist2 = dx * dx + dy * dy + dz * dz;
            if (dist2 > contact * contact || dist2 < 1e-12) {
                continue;
            }
            const dist = Math.sqrt(dist2);
            const correction = (dist - contact) / dist * RELAX_FACTOR * contactWeight;
            const mx = dx * correction;
            const my = dy * correction;
            const mz = dz * correction;
            pi[0] += mx;
            pi[1] += my;
            pi[2] += mz;
            pj[0] -= mx;
            pj[1] -= my;
            pj[2] -= mz;
            const move = Math.abs(mx) + Math.abs(my) + Math.abs(mz);
            if (move > maxMove) {
                maxMove = move;
            }
        }

        if (iteration > 0 && iteration % rebuildEvery === 0) {
            pairs = buildPairList(mol, coords, excluded, tight, cutoff);
        }

        if (maxMove < tolerance) {
            break;
        }
    }
    return iterations;
}

/**
 * 把顺反/共轭二面角精确落到目标值 (整支刚性旋转)。返回修正的键数。
 */
function enforceDihedralTargets(mol: Mol, coords: Vec3[],
                                targets: { [key: number]: { [key: string]: number } },
                                ringEdge: { [key: string]: number }): number {
    const elements: string[] = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        elements.push(mol.atoms[i].element);
    }
    let fixed = 0;
    const keys = Object.keys(targets).map((k) => parseInt(k, 10)).sort((a, b) => a - b);
    for (let ki = 0; ki < keys.length; ki++) {
        const k = keys[ki];
        const bond = mol.bonds[k];
        const a = bond.a;
        const b = bond.b;
        if (ringEdge[ringKey(a, b)] !== undefined) {
            continue;                      // 环内键交给环几何
        }
        const table = targets[k];
        // 选参照: 优先两端都不是氢的那一对
        let bestWeight = -1;
        let bestX = -1;
        let bestY = -1;
        let bestValue = 0.0;
        const pairKeys = Object.keys(table);
        for (let pi = 0; pi < pairKeys.length; pi++) {
            const parts = pairKeys[pi].split('|');
            const x = parseInt(parts[0], 10);
            const y = parseInt(parts[1], 10);
            let weight = 0;
            if (elements[x] === 'H') {
                weight += 1;
            }
            if (elements[y] === 'H') {
                weight += 1;
            }
            if (bestWeight < 0 || weight < bestWeight) {
                bestWeight = weight;
                bestX = x;
                bestY = y;
                bestValue = table[pairKeys[pi]];
            }
        }
        if (bestWeight < 0) {
            continue;
        }
        const x = bestX;
        const y = bestY;
        const wanted = bestValue;
        const current = dihedralAngle(coords[x], coords[a], coords[b], coords[y]);
        let delta = wanted - current;
        while (delta > 180.0) {
            delta -= 360.0;
        }
        while (delta < -180.0) {
            delta += 360.0;
        }
        if (Math.abs(delta) < 0.05) {
            continue;
        }
        const axis = vSub(coords[b], coords[a]);
        // 符号约定同 relax: 绕 a->b 轴正转会使二面角变小
        const step = current - wanted;
        // 把 b 端的所有分支一起转, 保持 b 上取代基之间的相对构型
        const listB = mol.adj[b];
        for (let t = 0; t < listB.length; t++) {
            const bond2 = mol.bonds[listB[t]];
            const other = bond2.b === b ? bond2.a : bond2.b;
            if (other === a) {
                continue;
            }
            const group = collectBranch(mol, b, other);
            for (let gi = 0; gi < group.length; gi++) {
                const index = group[gi];
                coords[index] = rotateAboutAxis(coords[index], axis, step, coords[b]);
            }
        }
        fixed += 1;
    }
    return fixed;
}

/** 把分子平移到几何中心位于原点 */
function centerCoordinates(coords: Vec3[]): void {
    const n = coords.length;
    if (n === 0) {
        return;
    }
    let sx = 0.0;
    let sy = 0.0;
    let sz = 0.0;
    for (let i = 0; i < n; i++) {
        sx += coords[i][0];
        sy += coords[i][1];
        sz += coords[i][2];
    }
    const cx = sx / n;
    const cy = sy / n;
    const cz = sz / n;
    for (let i = 0; i < n; i++) {
        coords[i][0] -= cx;
        coords[i][1] -= cy;
        coords[i][2] -= cz;
    }
}

// =====================================================================
// 10. 手性修正
// =====================================================================

/**
 * 四面体中心的取向符号。
 * +1 表示其余三个按逆时针排列 (对应 @@), -1 表示顺时针 (@), 0 表示无法判别。
 */
function chiralParity(coords: Vec3[], center: number, order: number[]): number {
    const c = coords[center];
    const p0 = coords[order[0]];
    const p1 = coords[order[1]];
    const p2 = coords[order[2]];
    const p3 = coords[order[3]];
    const normal = vCross(vSub(p2, p1), vSub(p3, p1));
    if (vLen(normal) < 1e-9) {
        return 0;
    }
    const projection = vDot(normal, vSub(p0, c));
    if (Math.abs(projection) < 1e-12) {
        return 0;
    }
    return projection > 0.0 ? 1 : -1;
}

/** 收集去掉 center 之后 start 所在的连通分支 (不含 center) */
function collectBranch(mol: Mol, center: number, start: number): number[] {
    const seen: { [key: number]: boolean } = {};
    seen[center] = true;
    seen[start] = true;
    const group: number[] = [start];
    const queue: number[] = [start];
    let head = 0;
    while (head < queue.length) {
        const cur = queue[head];
        head += 1;
        const list = mol.adj[cur];
        for (let t = 0; t < list.length; t++) {
            const bond = mol.bonds[list[t]];
            const nxt = bond.b === cur ? bond.a : bond.b;
            if (seen[nxt] !== true) {
                seen[nxt] = true;
                group.push(nxt);
                queue.push(nxt);
            }
        }
    }
    return group;
}

/**
 * 返回手性中心的取代基顺序; 无法处理时返回 null。
 */
function tetrahedralOrder(mol: Mol, center: number): number[] | null {
    const neighbors = neighborAtoms(mol, center);
    if (neighbors.length !== 4) {
        return null;
    }
    if (!mol.atoms[center].stereo_h_first) {
        return neighbors;
    }
    let hy = -1;
    for (let t = 0; t < neighbors.length; t++) {
        if (mol.atoms[neighbors[t]].element === 'H') {
            hy = neighbors[t];
            break;
        }
    }
    if (hy < 0) {
        return neighbors;
    }
    const ordered: number[] = [hy];
    for (let t = 0; t < neighbors.length; t++) {
        if (neighbors[t] !== hy) {
            ordered.push(neighbors[t]);
        }
    }
    return ordered;
}

/**
 * 在中心 center 上交换两个取代基 (及其整个分支) 的位置。
 * 绕轴 w = normalize(d_a + d_b) 把两个分支各旋转 180 度。
 */
function swapBranches(mol: Mol, coords: Vec3[], center: number, idxA: number, idxB: number): boolean {
    const c = coords[center];
    const directionA = vNorm(vSub(coords[idxA], c));
    const directionB = vNorm(vSub(coords[idxB], c));
    const axis = vAdd(directionA, directionB);
    if (vLen(axis) < 1e-6) {
        return false;                 // 两个取代基完全反向, 无法交换
    }
    const branchA = collectBranch(mol, center, idxA);
    const branchB = collectBranch(mol, center, idxB);
    for (let t = 0; t < branchA.length; t++) {
        coords[branchA[t]] = rotateAboutAxis(coords[branchA[t]], axis, 180.0, c);
    }
    for (let t = 0; t < branchB.length; t++) {
        coords[branchB[t]] = rotateAboutAxis(coords[branchB[t]], axis, 180.0, c);
    }
    return true;
}

/**
 * 按 @/@@ 修正四面体手性中心。返回被修正的中心个数。
 */
function fixChirality(mol: Mol, coords: Vec3[], opts: SmilesOptions): number {
    void opts;
    let fixed = 0;
    for (let center = 0; center < mol.atoms.length; center++) {
        const mark = mol.atoms[center].chirality;
        if (!mark) {
            continue;
        }
        const order = tetrahedralOrder(mol, center);
        if (order === null) {
            continue;
        }
        const want = mark === '@@' ? 1 : -1;
        if (chiralParity(coords, center, order) === want) {
            continue;
        }
        // 优先交换原子数最少的两个分支, 减少被移动的原子数
        const candidates: number[][] = [];
        for (let t = 0; t < order.length; t++) {
            candidates.push([collectBranch(mol, center, order[t]).length, order[t]]);
        }
        candidates.sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));
        if (candidates.length < 2) {
            continue;
        }
        const idxA = candidates[0][1];
        const idxB = candidates[1][1];
        if (swapBranches(mol, coords, center, idxA, idxB)) {
            if (chiralParity(coords, center, order) === want) {
                fixed += 1;
            }
        }
    }
    return fixed;
}

// =====================================================================
// 11. 公开 API
// =====================================================================

/** 单原子信息 (原子明细表的一项)。 */
export interface SmilesAtomDetail {
    index: number;
    element: string;
    aromatic: boolean;
    charge: number;
    isotope: number;
    implicitHydrogens: number;
    chirality: string;
    hybridization: string;
    x: number;
    y: number;
    z: number;
}

/** 单键信息 (键明细表的一项)。 */
export interface SmilesBondDetail {
    a: number;
    b: number;
    order: number;
    aromatic: boolean;
    length: number;
    idealLength: number;
}

/** 几何自检结果 */
export interface SmilesQuality {
    maxBondDeviation: number;
    worstBond: number;
    minNonbondedRatio: number;
    clashes: number;
}

/** analyze() 的分子信息 */
export interface SmilesInfo {
    smiles: string;
    numAtoms: number;
    numHeavyAtoms: number;
    numBonds: number;
    numRings: number;
    formula: string;
    molecularWeight: number;
    totalCharge: number;
    /** 兼容旧接口: [索引, 元素] 序列, coordinates 按下标对齐 */
    atoms: Array<[number, string]>;
    coordinates: number[][];
    atomDetails: SmilesAtomDetail[];
    bondDetails: SmilesBondDetail[];
    rings: number[][];
    fragments: number[][];
    stereoCenters: { atom: number; mark: string; parity: number }[];
    relaxIterations: number;
    quality: SmilesQuality;
    options: SmilesOptions;
}

/** analyze() 的返回 */
export interface SmilesAnalysisResult {
    ok: boolean;
    error: string;
    errorPos: number;
    errorText: string;
    smiles: string;
    info: SmilesInfo | null;
}

/**
 * 完整流程: 解析 -> 验价 -> 环感知 -> 加氢 -> 生成坐标 -> 松弛。
 * 不做任何 I/O: 不打印、不写文件、不抛异常。
 */
export function analyzeSmiles(smiles: unknown, options?: Partial<SmilesOptions>): SmilesAnalysisResult {
    const opts = makeOptions(options);
    const result: SmilesAnalysisResult = {
        ok: false,
        error: '',
        errorPos: -1,
        errorText: '',
        smiles: typeof smiles === 'string' ? smiles : '',
        info: null,
    };

    const parsed = parseSmiles(smiles);
    if (!parsed.ok) {
        result.error = parsed.error;
        result.errorPos = parsed.error_pos;
        result.errorText = parsed.error_text;
        return result;
    }

    const heavy = parsed.mol as Mol;
    computeImplicitHydrogens(heavy);

    const problems = checkValence(heavy);
    if (problems.length > 0) {
        result.error = ERR_VALENCE;
        result.errorPos = -1;
        result.errorText = problems.join('; ');
        return result;
    }

    // 环与扭转目标在"重原子图"上计算即可 (氢不可能成环)
    const rings = findRings(heavy);
    const maps = buildRingMaps(heavy, rings);
    const ringEdge = maps.ring_edge;
    const ringSize = maps.ring_size;

    let mol: Mol;
    if (opts.add_hydrogens) {
        mol = expandHydrogens(heavy);
    } else {
        // 保留 hcount (即隐式氢数目) 供调用方参考, 但不生成氢原子坐标
        mol = heavy;
    }

    if (mol.atoms.length > opts.max_atoms) {
        result.error = ERR_TOO_LARGE;
        result.errorText = '原子数 ' + mol.atoms.length + ' 超过上限 ' + opts.max_atoms;
        return result;
    }

    const targets = buildTorsionTargets(mol, rings, ringEdge);
    const coords = generateCoordinates(mol, targets, ringEdge, ringSize, opts);

    let iterations = 0;
    if (opts.optimize && coords.length > 1) {
        const constraints = buildConstraints(mol, targets, ringEdge, ringSize);
        const planarRings = buildPlanarRings(mol, rings);
        const sp2Centers = buildSp2Centers(mol, planarRings);
        const exclusions = buildExclusions(mol);
        const excluded = exclusions.excluded;
        const tight = exclusions.tight;
        // 先落位顺反/共轭二面角, 再松弛。
        enforceDihedralTargets(mol, coords, targets, ringEdge);
        iterations = relax(mol, coords, constraints, planarRings, sp2Centers,
                           excluded, tight, opts);
        // 松弛过程可能小幅扰动二面角, 再精确落一次
        enforceDihedralTargets(mol, coords, targets, ringEdge);
        // 碰撞消解: 用加大的非键排斥再压几轮, 专门拆开残存的近接触。
        const relief = makeOptions(opts);
        relief.relax_iterations = coords.length > 200 ? 60 : 150;
        relief.contact_weight = 2.0;
        iterations += relax(mol, coords, constraints, planarRings, sp2Centers,
                            excluded, tight, relief);
        enforceDihedralTargets(mol, coords, targets, ringEdge);
        if (opts.stereo) {
            const stereoOpts = makeOptions(opts);
            stereoOpts.relax_iterations = 60;
            const fixed = fixChirality(mol, coords, stereoOpts);
            if (fixed > 0) {
                iterations += relax(mol, coords, constraints, planarRings,
                                    sp2Centers, excluded, tight, stereoOpts);
                enforceDihedralTargets(mol, coords, targets, ringEdge);
            }
        }
    }

    if (opts.center) {
        centerCoordinates(coords);
    }

    result.info = buildInfo(typeof smiles === 'string' ? smiles : '', mol, heavy, coords, rings, iterations, opts);
    result.ok = true;
    return result;
}

/** 组装对外暴露的分子信息 (全部为 JSON 可序列化的基础类型) */
function buildInfo(smiles: string, mol: Mol, heavy: Mol, coords: Vec3[], rings: number[][],
                   iterations: number, opts: SmilesOptions): SmilesInfo {
    const atoms = mol.atoms;
    const elements: string[] = [];
    for (let i = 0; i < atoms.length; i++) {
        elements.push(atoms[i].element);
    }
    const n = atoms.length;

    const atomList: SmilesAtomDetail[] = [];
    for (let i = 0; i < n; i++) {
        const atom = atoms[i];
        const p = coords[i];
        atomList.push({
            index: i,
            element: atom.element,
            aromatic: atom.aromatic,
            charge: atom.charge,
            isotope: atom.isotope,
            implicitHydrogens: atom.hcount,
            chirality: atom.chirality,
            hybridization: HYBRID_NAMES[atomHybridization(mol, i)],
            x: p[0], y: p[1], z: p[2],
        });
    }

    const bondList: SmilesBondDetail[] = [];
    for (let k = 0; k < mol.bonds.length; k++) {
        const bond = mol.bonds[k];
        const a = bond.a;
        const b = bond.b;
        const ideal = bondLength(elements[a], elements[b], bond.order);
        bondList.push({
            a: a, b: b,
            order: bond.aromatic ? 1 : bond.order,
            aromatic: bond.aromatic,
            length: vDist(coords[a], coords[b]),
            idealLength: ideal,
        });
    }

    const ringList: number[][] = [];
    for (let r = 0; r < rings.length; r++) {
        ringList.push(rings[r].slice());
    }

    const formula = hillFormula(elements);
    let mass = 0.0;
    for (let i = 0; i < elements.length; i++) {
        mass += elementMass(elements[i]);
    }
    let totalCharge = 0;
    for (let i = 0; i < atoms.length; i++) {
        totalCharge += atoms[i].charge;
    }

    const stereoCenters: { atom: number; mark: string; parity: number }[] = [];
    for (let i = 0; i < n; i++) {
        if (atoms[i].chirality) {
            const order = tetrahedralOrder(mol, i);
            const parity = order !== null ? chiralParity(coords, i, order) : 0;
            stereoCenters.push({
                atom: i,
                mark: atoms[i].chirality,
                parity: parity,
            });
        }
    }

    const quality = measureQuality(mol, coords, bondList);

    // 兼容旧接口: atoms 为 [索引, 元素] 序列, coordinates 按下标对齐
    const atomPairs: Array<[number, string]> = [];
    for (let i = 0; i < n; i++) {
        atomPairs.push([i, elements[i]]);
    }
    const coordinates: Vec3[] = [];
    for (let i = 0; i < coords.length; i++) {
        coordinates.push([coords[i][0], coords[i][1], coords[i][2]]);
    }

    return {
        smiles: smiles,
        numAtoms: n,
        numHeavyAtoms: heavy.atoms.length,
        numBonds: bondList.length,
        numRings: rings.length,
        formula: formula,
        molecularWeight: mass,
        totalCharge: totalCharge,
        atoms: atomPairs,
        coordinates: coordinates,
        atomDetails: atomList,
        bondDetails: bondList,
        rings: ringList,
        fragments: findFragments(mol),
        stereoCenters: stereoCenters,
        relaxIterations: iterations,
        quality: quality,
        options: opts,
    };
}

/**
 * 几何自检: 键长偏差、最近非键接触、重叠原子数。
 * 冲突只统计拓扑距离 > 4 的原子对。
 */
function measureQuality(mol: Mol, coords: Vec3[], bondList: SmilesBondDetail[]): SmilesQuality {
    let maxBondDeviation = 0.0;
    let worstBond = -1;
    for (let k = 0; k < bondList.length; k++) {
        const item = bondList[k];
        const deviation = Math.abs(item.length - item.idealLength);
        if (deviation > maxBondDeviation) {
            maxBondDeviation = deviation;
            worstBond = k;
        }
    }

    const n = coords.length;
    const exclusions = buildExclusions(mol);
    const excluded = exclusions.excluded;
    const tight = exclusions.tight;
    const grid = new SpatialGrid(5.0);
    for (let i = 0; i < n; i++) {
        grid.add(i, coords[i]);
    }
    let minRatio: number | null = null;
    let closePairs = 0;
    for (let i = 0; i < n; i++) {
        const candidates = grid.query(coords[i]);
        for (let t = 0; t < candidates.length; t++) {
            const j = candidates[t];
            if (j <= i) {
                continue;
            }
            const key = pairKey(i, j, n);
            if (excluded[key] === true || tight[key] === true) {
                continue;
            }
            const d = vDist(coords[i], coords[j]);
            const contact = vdwRadius(mol.atoms[i].element) + vdwRadius(mol.atoms[j].element);
            if (contact <= 0.0) {
                continue;
            }
            const ratio = d / contact;
            if (minRatio === null || ratio < minRatio) {
                minRatio = ratio;
            }
            if (ratio < 0.65) {
                closePairs += 1;
            }
        }
    }
    return {
        maxBondDeviation: maxBondDeviation,
        worstBond: worstBond,
        // 最近非键接触距离 / 范德华半径和; -1 表示分子太小, 没有可统计的原子对
        minNonbondedRatio: minRatio !== null ? minRatio : -1.0,
        clashes: closePairs,
    };
}

/** Hill 排序分子式: C 优先, 然后 H, 其余按字母序 */
function hillFormula(elements: string[]): string {
    const counts: { [key: string]: number } = {};
    for (let i = 0; i < elements.length; i++) {
        let symbol = elements[i];
        if (symbol === '*') {
            symbol = 'R';
        }
        counts[symbol] = (counts[symbol] === undefined ? 0 : counts[symbol]) + 1;
    }
    const parts: Array<[string, number]> = [];
    if (counts['C'] !== undefined) {
        parts.push(['C', counts['C']]);
        delete counts['C'];
        if (counts['H'] !== undefined) {
            parts.push(['H', counts['H']]);
            delete counts['H'];
        }
    } else {
        if (counts['H'] !== undefined) {
            parts.push(['H', counts['H']]);
            delete counts['H'];
        }
    }
    const rest = Object.keys(counts).sort();
    for (let i = 0; i < rest.length; i++) {
        parts.push([rest[i], counts[rest[i]]]);
    }
    let text = '';
    for (let i = 0; i < parts.length; i++) {
        text += parts[i][0];
        if (parts[i][1] > 1) {
            text += String(parts[i][1]);
        }
    }
    return text;
}

/** smilesTo3dSafe 的返回: 永不抛出, 宿主只需检查 ok。 */
export interface SmilesTo3dOutcome {
    ok: boolean;
    coordinates: number[][] | null;
    info: SmilesInfo | null;
    errorText: string;
}

/**
 * 推荐接口: 生成三维坐标, 任何失败都通过返回值表达, 永不抛异常。
 */
export function smilesTo3dSafe(smiles: unknown, options?: Partial<SmilesOptions>): SmilesTo3dOutcome {
    try {
        const result = analyzeSmiles(smiles, options);
        if (!result.ok || result.info === null) {
            return { ok: false, coordinates: null, info: null, errorText: result.errorText };
        }
        return { ok: true, coordinates: result.info.coordinates, info: result.info, errorText: '' };
    } catch (e) {
        // 兜底: 保证调用方不会因内部异常而中断 (本模块契约是永不抛出、不打印)
        const msg = e instanceof Error ? e.message : String(e);
        return { ok: false, coordinates: null, info: null, errorText: '内部错误: ' + msg };
    }
}

/**
 * 轻量校验 (不生成坐标): 解析 + 隐式氢 + 价态检查。
 * 适合在输入框的 validateInput 回调里逐键调用。
 */
export function validateSmiles(smiles: unknown): { ok: boolean; errorText: string } {
    if (typeof smiles !== 'string') {
        return { ok: false, errorText: ERROR_TEXT[ERR_EMPTY] };
    }
    const parsed = parseSmiles(smiles);
    if (!parsed.ok) {
        return { ok: false, errorText: parsed.error_text };
    }
    const mol = parsed.mol as Mol;
    computeImplicitHydrogens(mol);
    const problems = checkValence(mol);
    if (problems.length > 0) {
        return { ok: false, errorText: problems.join('; ') };
    }
    return { ok: true, errorText: '' };
}

/** %-2s: 左对齐宽度 2。 */
function padRight2(s: string): string {
    return s.length >= 2 ? s : s + ' ';
}

/** %12.6f: 右对齐宽度 12, 6 位小数。 */
function fmtF12(v: number): string {
    let s = v.toFixed(6);
    while (s.length < 12) {
        s = ' ' + s;
    }
    return s;
}

/** 把结果输出为 XYZ 文本 (不写文件) */
export function toXyz(info: SmilesInfo, comment?: string): string {
    const coordinates = info.coordinates;
    const atoms = info.atoms;
    const lines: string[] = [];
    lines.push(String(atoms.length));
    const title = comment ? comment : (info.smiles || '');
    lines.push(title);
    for (let i = 0; i < atoms.length; i++) {
        const pair = atoms[i];
        const p = coordinates[pair[0]];
        const symbol = pair[1] !== '*' ? pair[1] : 'X';
        lines.push(padRight2(symbol) + ' ' + fmtF12(p[0]) + ' ' + fmtF12(p[1]) + ' ' + fmtF12(p[2]));
    }
    return lines.join('\n') + '\n';
}
