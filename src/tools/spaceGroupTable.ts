/**
 * A space-group name turned into the operations it names.
 *
 * A CIF is supposed to list its symmetry operations, but a fair number state
 * only a space-group name and expect the reader to know the group. Read as P1,
 * such a file shows one asymmetric unit as if it were the whole cell, which is
 * wrong in a way that looks plausible. This table closes that gap.
 *
 * The data is International Tables for Crystallography Vol. B (2001), Table
 * A1.4.2.7: one Hall symbol (Hall 1981) for every setting of all 230 space
 * groups, 530 settings in all, including both origins where a group has two
 * and both the hexagonal and rhombohedral axes of the seven R groups. A Hall
 * symbol is a compact list of generators, so the operations are produced by
 * closing the group rather than stored one by one: the whole table costs a few
 * kilobytes instead of the tens of thousands of operation strings an explicit
 * list would need.
 *
 * All arithmetic is exact. A translation is carried in units of 1/24, the
 * least common multiple of every denominator a space group needs, so folding
 * an image back into the cell compares integers and cannot turn a third into
 * 0.3333333333. The operations leave through `symOpToString`, the same helper
 * the file readers use, so an operation read from this table is spelled
 * exactly like the one the same structure's CIF would state.
 *
 * The CIF reader consults this table only when the file has no symmetry tag at
 * all. A file that states operations is the authority on its own structure;
 * a file whose operation tag is present but unreadable is reported as P1, not
 * quietly given some other group's operations.
 */

import { SymOp, symOpToString } from './crystalGeometry';

/**
 * `number|setting|Hermann-Mauguin|Hall symbol`, in the order of the IT table.
 * `setting` is the table's setting code (`b1`, `cab`, `1`, `2`, `h`, `r`, or
 * empty for the standard setting); the HM field may repeat an origin or axis
 * choice as a `:1`/`:2`/`:h`/`:r` suffix, which is stripped before matching.
 */
const SG_HALL_DATA: string[] = [
    '1||P 1|p 1',
    '2||P -1|-p 1',
    '3|b|P 1 2 1|p 2y',
    '3|c|P 1 1 2|p 2',
    '3|a|P 2 1 1|p 2x',
    '4|b|P 1 21 1|p 2yb',
    '4|c|P 1 1 21|p 2c',
    '4|a|P 21 1 1|p 2xa',
    '5|b1|C 1 2 1|c 2y',
    '5|b2|A 1 2 1|a 2y',
    '5|b3|I 1 2 1|i 2y',
    '5|c1|A 1 1 2|a 2',
    '5|c2|B 1 1 2|b 2',
    '5|c3|I 1 1 2|i 2',
    '5|a1|B 2 1 1|b 2x',
    '5|a2|C 2 1 1|c 2x',
    '5|a3|I 2 1 1|i 2x',
    '6|b|P 1 m 1|p -2y',
    '6|c|P 1 1 m|p -2',
    '6|a|P m 1 1|p -2x',
    '7|b1|P 1 c 1|p -2yc',
    '7|b2|P 1 n 1|p -2yac',
    '7|b3|P 1 a 1|p -2ya',
    '7|c1|P 1 1 a|p -2a',
    '7|c2|P 1 1 n|p -2ab',
    '7|c3|P 1 1 b|p -2b',
    '7|a1|P b 1 1|p -2xb',
    '7|a2|P n 1 1|p -2xbc',
    '7|a3|P c 1 1|p -2xc',
    '8|b1|C 1 m 1|c -2y',
    '8|b2|A 1 m 1|a -2y',
    '8|b3|I 1 m 1|i -2y',
    '8|c1|A 1 1 m|a -2',
    '8|c2|B 1 1 m|b -2',
    '8|c3|I 1 1 m|i -2',
    '8|a1|B m 1 1|b -2x',
    '8|a2|C m 1 1|c -2x',
    '8|a3|I m 1 1|i -2x',
    '9|b1|C 1 c 1|c -2yc',
    '9|b2|A 1 n 1|a -2yab',
    '9|b3|I 1 a 1|i -2ya',
    '9|-b1|A 1 a 1|a -2ya',
    '9|-b2|C 1 n 1|c -2yac',
    '9|-b3|I 1 c 1|i -2yc',
    '9|c1|A 1 1 a|a -2a',
    '9|c2|B 1 1 n|b -2ab',
    '9|c3|I 1 1 b|i -2b',
    '9|-c1|B 1 1 b|b -2b',
    '9|-c2|A 1 1 n|a -2ab',
    '9|-c3|I 1 1 a|i -2a',
    '9|a1|B b 1 1|b -2xb',
    '9|a2|C n 1 1|c -2xac',
    '9|a3|I c 1 1|i -2xc',
    '9|-a1|C c 1 1|c -2xc',
    '9|-a2|B n 1 1|b -2xab',
    '9|-a3|I b 1 1|i -2xb',
    '10|b|P 1 2/m 1|-p 2y',
    '10|c|P 1 1 2/m|-p 2',
    '10|a|P 2/m 1 1|-p 2x',
    '11|b|P 1 21/m 1|-p 2yb',
    '11|c|P 1 1 21/m|-p 2c',
    '11|a|P 21/m 1 1|-p 2xa',
    '12|b1|C 1 2/m 1|-c 2y',
    '12|b2|A 1 2/m 1|-a 2y',
    '12|b3|I 1 2/m 1|-i 2y',
    '12|c1|A 1 1 2/m|-a 2',
    '12|c2|B 1 1 2/m|-b 2',
    '12|c3|I 1 1 2/m|-i 2',
    '12|a1|B 2/m 1 1|-b 2x',
    '12|a2|C 2/m 1 1|-c 2x',
    '12|a3|I 2/m 1 1|-i 2x',
    '13|b1|P 1 2/c 1|-p 2yc',
    '13|b2|P 1 2/n 1|-p 2yac',
    '13|b3|P 1 2/a 1|-p 2ya',
    '13|c1|P 1 1 2/a|-p 2a',
    '13|c2|P 1 1 2/n|-p 2ab',
    '13|c3|P 1 1 2/b|-p 2b',
    '13|a1|P 2/b 1 1|-p 2xb',
    '13|a2|P 2/n 1 1|-p 2xbc',
    '13|a3|P 2/c 1 1|-p 2xc',
    '14|b1|P 1 21/c 1|-p 2ybc',
    '14|b2|P 1 21/n 1|-p 2yn',
    '14|b3|P 1 21/a 1|-p 2yab',
    '14|c1|P 1 1 21/a|-p 2ac',
    '14|c2|P 1 1 21/n|-p 2n',
    '14|c3|P 1 1 21/b|-p 2bc',
    '14|a1|P 21/b 1 1|-p 2xab',
    '14|a2|P 21/n 1 1|-p 2xn',
    '14|a3|P 21/c 1 1|-p 2xac',
    '15|b1|C 1 2/c 1|-c 2yc',
    '15|b2|A 1 2/n 1|-a 2yab',
    '15|b3|I 1 2/a 1|-i 2ya',
    '15|-b1|A 1 2/a 1|-a 2ya',
    '15|-b2|C 1 2/n 1|-c 2yac',
    '15|-b3|I 1 2/c 1|-i 2yc',
    '15|c1|A 1 1 2/a|-a 2a',
    '15|c2|B 1 1 2/n|-b 2ab',
    '15|c3|I 1 1 2/b|-i 2b',
    '15|-c1|B 1 1 2/b|-b 2b',
    '15|-c2|A 1 1 2/n|-a 2ab',
    '15|-c3|I 1 1 2/a|-i 2a',
    '15|a1|B 2/b 1 1|-b 2xb',
    '15|a2|C 2/n 1 1|-c 2xac',
    '15|a3|I 2/c 1 1|-i 2xc',
    '15|-a1|C 2/c 1 1|-c 2xc',
    '15|-a2|B 2/n 1 1|-b 2xab',
    '15|-a3|I 2/b 1 1|-i 2xb',
    '16||P 2 2 2|p 2 2',
    '17||P 2 2 21|p 2c 2',
    '17|cab|P 21 2 2|p 2a 2a',
    '17|bca|P 2 21 2|p 2 2b',
    '18||P 21 21 2|p 2 2ab',
    '18|cab|P 2 21 21|p 2bc 2',
    '18|bca|P 21 2 21|p 2ac 2ac',
    '19||P 21 21 21|p 2ac 2ab',
    '20||C 2 2 21|c 2c 2',
    '20|cab|A 21 2 2|a 2a 2a',
    '20|bca|B 2 21 2|b 2 2b',
    '21||C 2 2 2|c 2 2',
    '21|cab|A 2 2 2|a 2 2',
    '21|bca|B 2 2 2|b 2 2',
    '22||F 2 2 2|f 2 2',
    '23||I 2 2 2|i 2 2',
    '24||I 21 21 21|i 2b 2c',
    '25||P m m 2|p 2 -2',
    '25|cab|P 2 m m|p -2 2',
    '25|bca|P m 2 m|p -2 -2',
    '26||P m c 21|p 2c -2',
    '26|ba-c|P c m 21|p 2c -2c',
    '26|cab|P 21 m a|p -2a 2a',
    '26|-cba|P 21 a m|p -2 2a',
    '26|bca|P b 21 m|p -2 -2b',
    '26|a-cb|P m 21 b|p -2b -2',
    '27||P c c 2|p 2 -2c',
    '27|cab|P 2 a a|p -2a 2',
    '27|bca|P b 2 b|p -2b -2b',
    '28||P m a 2|p 2 -2a',
    '28|ba-c|P b m 2|p 2 -2b',
    '28|cab|P 2 m b|p -2b 2',
    '28|-cba|P 2 c m|p -2c 2',
    '28|bca|P c 2 m|p -2c -2c',
    '28|a-cb|P m 2 a|p -2a -2a',
    '29||P c a 21|p 2c -2ac',
    '29|ba-c|P b c 21|p 2c -2b',
    '29|cab|P 21 a b|p -2b 2a',
    '29|-cba|P 21 c a|p -2ac 2a',
    '29|bca|P c 21 b|p -2bc -2c',
    '29|a-cb|P b 21 a|p -2a -2ab',
    '30||P n c 2|p 2 -2bc',
    '30|ba-c|P c n 2|p 2 -2ac',
    '30|cab|P 2 n a|p -2ac 2',
    '30|-cba|P 2 a n|p -2ab 2',
    '30|bca|P b 2 n|p -2ab -2ab',
    '30|a-cb|P n 2 b|p -2bc -2bc',
    '31||P m n 21|p 2ac -2',
    '31|ba-c|P n m 21|p 2bc -2bc',
    '31|cab|P 21 m n|p -2ab 2ab',
    '31|-cba|P 21 n m|p -2 2ac',
    '31|bca|P n 21 m|p -2 -2bc',
    '31|a-cb|P m 21 n|p -2ab -2',
    '32||P b a 2|p 2 -2ab',
    '32|cab|P 2 c b|p -2bc 2',
    '32|bca|P c 2 a|p -2ac -2ac',
    '33||P n a 21|p 2c -2n',
    '33|ba-c|P b n 21|p 2c -2ab',
    '33|cab|P 21 n b|p -2bc 2a',
    '33|-cba|P 21 c n|p -2n 2a',
    '33|bca|P c 21 n|p -2n -2ac',
    '33|a-cb|P n 21 a|p -2ac -2n',
    '34||P n n 2|p 2 -2n',
    '34|cab|P 2 n n|p -2n 2',
    '34|bca|P n 2 n|p -2n -2n',
    '35||C m m 2|c 2 -2',
    '35|cab|A 2 m m|a -2 2',
    '35|bca|B m 2 m|b -2 -2',
    '36||C m c 21|c 2c -2',
    '36|ba-c|C c m 21|c 2c -2c',
    '36|cab|A 21 m a|a -2a 2a',
    '36|-cba|A 21 a m|a -2 2a',
    '36|bca|B b 21 m|b -2 -2b',
    '36|a-cb|B m 21 b|b -2b -2',
    '37||C c c 2|c 2 -2c',
    '37|cab|A 2 a a|a -2a 2',
    '37|bca|B b 2 b|b -2b -2b',
    '38||A m m 2|a 2 -2',
    '38|ba-c|B m m 2|b 2 -2',
    '38|cab|B 2 m m|b -2 2',
    '38|-cba|C 2 m m|c -2 2',
    '38|bca|C m 2 m|c -2 -2',
    '38|a-cb|A m 2 m|a -2 -2',
    '39||A b m 2|a 2 -2b',
    '39|ba-c|B m a 2|b 2 -2a',
    '39|cab|B 2 c m|b -2a 2',
    '39|-cba|C 2 m b|c -2a 2',
    '39|bca|C m 2 a|c -2a -2a',
    '39|a-cb|A c 2 m|a -2b -2b',
    '40||A m a 2|a 2 -2a',
    '40|ba-c|B b m 2|b 2 -2b',
    '40|cab|B 2 m b|b -2b 2',
    '40|-cba|C 2 c m|c -2c 2',
    '40|bca|C c 2 m|c -2c -2c',
    '40|a-cb|A m 2 a|a -2a -2a',
    '41||A b a 2|a 2 -2ab',
    '41|ba-c|B b a 2|b 2 -2ab',
    '41|cab|B 2 c b|b -2ab 2',
    '41|-cba|C 2 c b|c -2ac 2',
    '41|bca|C c 2 a|c -2ac -2ac',
    '41|a-cb|A c 2 a|a -2ab -2ab',
    '42||F m m 2|f 2 -2',
    '42|cab|F 2 m m|f -2 2',
    '42|bca|F m 2 m|f -2 -2',
    '43||F d d 2|f 2 -2d',
    '43|cab|F 2 d d|f -2d 2',
    '43|bca|F d 2 d|f -2d -2d',
    '44||I m m 2|i 2 -2',
    '44|cab|I 2 m m|i -2 2',
    '44|bca|I m 2 m|i -2 -2',
    '45||I b a 2|i 2 -2c',
    '45|cab|I 2 c b|i -2a 2',
    '45|bca|I c 2 a|i -2b -2b',
    '46||I m a 2|i 2 -2a',
    '46|ba-c|I b m 2|i 2 -2b',
    '46|cab|I 2 m b|i -2b 2',
    '46|-cba|I 2 c m|i -2c 2',
    '46|bca|I c 2 m|i -2c -2c',
    '46|a-cb|I m 2 a|i -2a -2a',
    '47||P m m m|-p 2 2',
    '48|1|P n n n:1|p 2 2 -1n',
    '48|2|P n n n:2|-p 2ab 2bc',
    '49||P c c m|-p 2 2c',
    '49|cab|P m a a|-p 2a 2',
    '49|bca|P b m b|-p 2b 2b',
    '50|1|P b a n:1|p 2 2 -1ab',
    '50|2|P b a n:2|-p 2ab 2b',
    '50|1cab|P n c b:1|p 2 2 -1bc',
    '50|2cab|P n c b:2|-p 2b 2bc',
    '50|1bca|P c n a:1|p 2 2 -1ac',
    '50|2bca|P c n a:2|-p 2a 2c',
    '51||P m m a|-p 2a 2a',
    '51|ba-c|P m m b|-p 2b 2',
    '51|cab|P b m m|-p 2 2b',
    '51|-cba|P c m m|-p 2c 2c',
    '51|bca|P m c m|-p 2c 2',
    '51|a-cb|P m a m|-p 2 2a',
    '52||P n n a|-p 2a 2bc',
    '52|ba-c|P n n b|-p 2b 2n',
    '52|cab|P b n n|-p 2n 2b',
    '52|-cba|P c n n|-p 2ab 2c',
    '52|bca|P n c n|-p 2ab 2n',
    '52|a-cb|P n a n|-p 2n 2bc',
    '53||P m n a|-p 2ac 2',
    '53|ba-c|P n m b|-p 2bc 2bc',
    '53|cab|P b m n|-p 2ab 2ab',
    '53|-cba|P c n m|-p 2 2ac',
    '53|bca|P n c m|-p 2 2bc',
    '53|a-cb|P m a n|-p 2ab 2',
    '54||P c c a|-p 2a 2ac',
    '54|ba-c|P c c b|-p 2b 2c',
    '54|cab|P b a a|-p 2a 2b',
    '54|-cba|P c a a|-p 2ac 2c',
    '54|bca|P b c b|-p 2bc 2b',
    '54|a-cb|P b a b|-p 2b 2ab',
    '55||P b a m|-p 2 2ab',
    '55|cab|P m c b|-p 2bc 2',
    '55|bca|P c m a|-p 2ac 2ac',
    '56||P c c n|-p 2ab 2ac',
    '56|cab|P n a a|-p 2ac 2bc',
    '56|bca|P b n b|-p 2bc 2ab',
    '57||P b c m|-p 2c 2b',
    '57|ba-c|P c a m|-p 2c 2ac',
    '57|cab|P m c a|-p 2ac 2a',
    '57|-cba|P m a b|-p 2b 2a',
    '57|bca|P b m a|-p 2a 2ab',
    '57|a-cb|P c m b|-p 2bc 2c',
    '58||P n n m|-p 2 2n',
    '58|cab|P m n n|-p 2n 2',
    '58|bca|P n m n|-p 2n 2n',
    '59|1|P m m n:1|p 2 2ab -1ab',
    '59|2|P m m n:2|-p 2ab 2a',
    '59|1cab|P n m m:1|p 2bc 2 -1bc',
    '59|2cab|P n m m:2|-p 2c 2bc',
    '59|1bca|P m n m:1|p 2ac 2ac -1ac',
    '59|2bca|P m n m:2|-p 2c 2a',
    '60||P b c n|-p 2n 2ab',
    '60|ba-c|P c a n|-p 2n 2c',
    '60|cab|P n c a|-p 2a 2n',
    '60|-cba|P n a b|-p 2bc 2n',
    '60|bca|P b n a|-p 2ac 2b',
    '60|a-cb|P c n b|-p 2b 2ac',
    '61||P b c a|-p 2ac 2ab',
    '61|ba-c|P c a b|-p 2bc 2ac',
    '62||P n m a|-p 2ac 2n',
    '62|ba-c|P m n b|-p 2bc 2a',
    '62|cab|P b n m|-p 2c 2ab',
    '62|-cba|P c m n|-p 2n 2ac',
    '62|bca|P m c n|-p 2n 2a',
    '62|a-cb|P n a m|-p 2c 2n',
    '63||C m c m|-c 2c 2',
    '63|ba-c|C c m m|-c 2c 2c',
    '63|cab|A m m a|-a 2a 2a',
    '63|-cba|A m a m|-a 2 2a',
    '63|bca|B b m m|-b 2 2b',
    '63|a-cb|B m m b|-b 2b 2',
    '64||C m c a|-c 2ac 2',
    '64|ba-c|C c m b|-c 2ac 2ac',
    '64|cab|A b m a|-a 2ab 2ab',
    '64|-cba|A c a m|-a 2 2ab',
    '64|bca|B b c m|-b 2 2ab',
    '64|a-cb|B m a b|-b 2ab 2',
    '65||C m m m|-c 2 2',
    '65|cab|A m m m|-a 2 2',
    '65|bca|B m m m|-b 2 2',
    '66||C c c m|-c 2 2c',
    '66|cab|A m a a|-a 2a 2',
    '66|bca|B b m b|-b 2b 2b',
    '67||C m m a|-c 2a 2',
    '67|ba-c|C m m b|-c 2a 2a',
    '67|cab|A b m m|-a 2b 2b',
    '67|-cba|A c m m|-a 2 2b',
    '67|bca|B m c m|-b 2 2a',
    '67|a-cb|B m a m|-b 2a 2',
    '68|1|C c c a:1|c 2 2 -1ac',
    '68|2|C c c a:2|-c 2a 2ac',
    '68|1ba-c|C c c b:1|c 2 2 -1ac',
    '68|2ba-c|C c c b:2|-c 2a 2c',
    '68|1cab|A b a a:1|a 2 2 -1ab',
    '68|2cab|A b a a:2|-a 2a 2b',
    '68|1-cba|A c a a:1|a 2 2 -1ab',
    '68|2-cba|A c a a:2|-a 2ab 2b',
    '68|1bca|B b c b:1|b 2 2 -1ab',
    '68|2bca|B b c b:2|-b 2ab 2b',
    '68|1a-cb|B b a b:1|b 2 2 -1ab',
    '68|2a-cb|B b a b:2|-b 2b 2ab',
    '69||F m m m|-f 2 2',
    '70|1|F d d d:1|f 2 2 -1d',
    '70|2|F d d d:2|-f 2uv 2vw',
    '71||I m m m|-i 2 2',
    '72||I b a m|-i 2 2c',
    '72|cab|I m c b|-i 2a 2',
    '72|bca|I c m a|-i 2b 2b',
    '73||I b c a|-i 2b 2c',
    '73|ba-c|I c a b|-i 2a 2b',
    '74||I m m a|-i 2b 2',
    '74|ba-c|I m m b|-i 2a 2a',
    '74|cab|I b m m|-i 2c 2c',
    '74|-cba|I c m m|-i 2 2b',
    '74|bca|I m c m|-i 2 2a',
    '74|a-cb|I m a m|-i 2c 2',
    '75||P 4|p 4',
    '76||P 41|p 4w',
    '77||P 42|p 4c',
    '78||P 43|p 4cw',
    '79||I 4|i 4',
    '80||I 41|i 4bw',
    '81||P -4|p -4',
    '82||I -4|i -4',
    '83||P 4/m|-p 4',
    '84||P 42/m|-p 4c',
    '85|1|P 4/n:1|p 4ab -1ab',
    '85|2|P 4/n:2|-p 4a',
    '86|1|P 42/n:1|p 4n -1n',
    '86|2|P 42/n:2|-p 4bc',
    '87||I 4/m|-i 4',
    '88|1|I 41/a:1|i 4bw -1bw',
    '88|2|I 41/a:2|-i 4ad',
    '89||P 4 2 2|p 4 2',
    '90||P 4 21 2|p 4ab 2ab',
    '91||P 41 2 2|p 4w 2c',
    '92||P 41 21 2|p 4abw 2nw',
    '93||P 42 2 2|p 4c 2',
    '94||P 42 21 2|p 4n 2n',
    '95||P 43 2 2|p 4cw 2c',
    '96||P 43 21 2|p 4nw 2abw',
    '97||I 4 2 2|i 4 2',
    '98||I 41 2 2|i 4bw 2bw',
    '99||P 4 m m|p 4 -2',
    '100||P 4 b m|p 4 -2ab',
    '101||P 42 c m|p 4c -2c',
    '102||P 42 n m|p 4n -2n',
    '103||P 4 c c|p 4 -2c',
    '104||P 4 n c|p 4 -2n',
    '105||P 42 m c|p 4c -2',
    '106||P 42 b c|p 4c -2ab',
    '107||I 4 m m|i 4 -2',
    '108||I 4 c m|i 4 -2c',
    '109||I 41 m d|i 4bw -2',
    '110||I 41 c d|i 4bw -2c',
    '111||P -4 2 m|p -4 2',
    '112||P -4 2 c|p -4 2c',
    '113||P -4 21 m|p -4 2ab',
    '114||P -4 21 c|p -4 2n',
    '115||P -4 m 2|p -4 -2',
    '116||P -4 c 2|p -4 -2c',
    '117||P -4 b 2|p -4 -2ab',
    '118||P -4 n 2|p -4 -2n',
    '119||I -4 m 2|i -4 -2',
    '120||I -4 c 2|i -4 -2c',
    '121||I -4 2 m|i -4 2',
    '122||I -4 2 d|i -4 2bw',
    '123||P 4/m m m|-p 4 2',
    '124||P 4/m c c|-p 4 2c',
    '125|1|P 4/n b m:1|p 4 2 -1ab',
    '125|2|P 4/n b m:2|-p 4a 2b',
    '126|1|P 4/n n c:1|p 4 2 -1n',
    '126|2|P 4/n n c:2|-p 4a 2bc',
    '127||P 4/m b m|-p 4 2ab',
    '128||P 4/m n c|-p 4 2n',
    '129|1|P 4/n m m:1|p 4ab 2ab -1ab',
    '129|2|P 4/n m m:2|-p 4a 2a',
    '130|1|P 4/n c c:1|p 4ab 2n -1ab',
    '130|2|P 4/n c c:2|-p 4a 2ac',
    '131||P 42/m m c|-p 4c 2',
    '132||P 42/m c m|-p 4c 2c',
    '133|1|P 42/n b c:1|p 4n 2c -1n',
    '133|2|P 42/n b c:2|-p 4ac 2b',
    '134|1|P 42/n n m:1|p 4n 2 -1n',
    '134|2|P 42/n n m:2|-p 4ac 2bc',
    '135||P 42/m b c|-p 4c 2ab',
    '136||P 42/m n m|-p 4n 2n',
    '137|1|P 42/n m c:1|p 4n 2n -1n',
    '137|2|P 42/n m c:2|-p 4ac 2a',
    '138|1|P 42/n c m:1|p 4n 2ab -1n',
    '138|2|P 42/n c m:2|-p 4ac 2ac',
    '139||I 4/m m m|-i 4 2',
    '140||I 4/m c m|-i 4 2c',
    '141|1|I 41/a m d:1|i 4bw 2bw -1bw',
    '141|2|I 41/a m d:2|-i 4bd 2',
    '142|1|I 41/a c d:1|i 4bw 2aw -1bw',
    '142|2|I 41/a c d:2|-i 4bd 2c',
    '143||P 3|p 3',
    '144||P 31|p 31',
    '145||P 32|p 32',
    '146|h|R 3:h|r 3',
    '146|r|R 3:r|p 3*',
    '147||P -3|-p 3',
    '148|h|R -3:h|-r 3',
    '148|r|R -3:r|-p 3*',
    '149||P 3 1 2|p 3 2',
    '150||P 3 2 1|p 3 2"',
    '151||P 31 1 2|p 31 2 (0 0 4)',
    '152||P 31 2 1|p 31 2"',
    '153||P 32 1 2|p 32 2 (0 0 2)',
    '154||P 32 2 1|p 32 2"',
    '155|h|R 3 2:h|r 3 2"',
    '155|r|R 3 2:r|p 3* 2',
    '156||P 3 m 1|p 3 -2"',
    '157||P 3 1 m|p 3 -2',
    '158||P 3 c 1|p 3 -2"c',
    '159||P 3 1 c|p 3 -2c',
    '160|h|R 3 m:h|r 3 -2"',
    '160|r|R 3 m:r|p 3* -2',
    '161|h|R 3 c:h|r 3 -2"c',
    '161|r|R 3 c:r|p 3* -2n',
    '162||P -3 1 m|-p 3 2',
    '163||P -3 1 c|-p 3 2c',
    '164||P -3 m 1|-p 3 2"',
    '165||P -3 c 1|-p 3 2"c',
    '166|h|R -3 m:h|-r 3 2"',
    '166|r|R -3 m:r|-p 3* 2',
    '167|h|R -3 c:h|-r 3 2"c',
    '167|r|R -3 c:r|-p 3* 2n',
    '168||P 6|p 6',
    '169||P 61|p 61',
    '170||P 65|p 65',
    '171||P 62|p 62',
    '172||P 64|p 64',
    '173||P 63|p 6c',
    '174||P -6|p -6',
    '175||P 6/m|-p 6',
    '176||P 63/m|-p 6c',
    '177||P 6 2 2|p 6 2',
    '178||P 61 2 2|p 61 2 (0 0 5)',
    '179||P 65 2 2|p 65 2 (0 0 1)',
    '180||P 62 2 2|p 62 2 (0 0 4)',
    '181||P 64 2 2|p 64 2 (0 0 2)',
    '182||P 63 2 2|p 6c 2c',
    '183||P 6 m m|p 6 -2',
    '184||P 6 c c|p 6 -2c',
    '185||P 63 c m|p 6c -2',
    '186||P 63 m c|p 6c -2c',
    '187||P -6 m 2|p -6 2',
    '188||P -6 c 2|p -6c 2',
    '189||P -6 2 m|p -6 -2',
    '190||P -6 2 c|p -6c -2c',
    '191||P 6/m m m|-p 6 2',
    '192||P 6/m c c|-p 6 2c',
    '193||P 63/m c m|-p 6c 2',
    '194||P 63/m m c|-p 6c 2c',
    '195||P 2 3|p 2 2 3',
    '196||F 2 3|f 2 2 3',
    '197||I 2 3|i 2 2 3',
    '198||P 21 3|p 2ac 2ab 3',
    '199||I 21 3|i 2b 2c 3',
    '200||P m -3|-p 2 2 3',
    '201|1|P n -3:1|p 2 2 3 -1n',
    '201|2|P n -3:2|-p 2ab 2bc 3',
    '202||F m -3|-f 2 2 3',
    '203|1|F d -3:1|f 2 2 3 -1d',
    '203|2|F d -3:2|-f 2uv 2vw 3',
    '204||I m -3|-i 2 2 3',
    '205||P a -3|-p 2ac 2ab 3',
    '206||I a -3|-i 2b 2c 3',
    '207||P 4 3 2|p 4 2 3',
    '208||P 42 3 2|p 4n 2 3',
    '209||F 4 3 2|f 4 2 3',
    '210||F 41 3 2|f 4d 2 3',
    '211||I 4 3 2|i 4 2 3',
    '212||P 43 3 2|p 4acd 2ab 3',
    '213||P 41 3 2|p 4bd 2ab 3',
    '214||I 41 3 2|i 4bd 2c 3',
    '215||P -4 3 m|p -4 2 3',
    '216||F -4 3 m|f -4 2 3',
    '217||I -4 3 m|i -4 2 3',
    '218||P -4 3 n|p -4n 2 3',
    '219||F -4 3 c|f -4a 2 3',
    '220||I -4 3 d|i -4bd 2c 3',
    '221||P m -3 m|-p 4 2 3',
    '222|1|P n -3 n:1|p 4 2 3 -1n',
    '222|2|P n -3 n:2|-p 4a 2bc 3',
    '223||P m -3 n|-p 4n 2 3',
    '224|1|P n -3 m:1|p 4n 2 3 -1n',
    '224|2|P n -3 m:2|-p 4bc 2bc 3',
    '225||F m -3 m|-f 4 2 3',
    '226||F m -3 c|-f 4a 2 3',
    '227|1|F d -3 m:1|f 4d 2 3 -1d',
    '227|2|F d -3 m:2|-f 4vw 2vw 3',
    '228|1|F d -3 c:1|f 4d 2 3 -1ad',
    '228|2|F d -3 c:2|-f 4ud 2vw 3',
    '229||I m -3 m|-i 4 2 3',
    '230||I a -3 d|-i 4bd 2c 3',
];

// ===== Exact integer operations =====

const DEN = 24;

interface IntOp {
    R: number[][];
    t: number[];
}

const I3: number[][] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/**
 * The principal rotation matrices, by axis and order, in the Hall convention
 * (right-handed, with the first non-zero entry positive). x, y and z are the
 * cell axes; the matrices act on fractional column vectors.
 */
const PRIN: { [axis: string]: { [order: number]: number[][] } } = {
    x: {
        1: I3,
        2: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
        3: [[1, 0, 0], [0, 0, -1], [0, 1, -1]],
        4: [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
        6: [[1, 0, 0], [0, 1, -1], [0, 1, 0]]
    },
    y: {
        1: I3,
        2: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
        3: [[-1, 0, 1], [0, 1, 0], [-1, 0, 0]],
        4: [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
        6: [[0, 0, 1], [0, 1, 0], [-1, 0, 1]]
    },
    z: {
        1: I3,
        2: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        3: [[0, -1, 0], [1, -1, 0], [0, 0, 1]],
        4: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        6: [[1, -1, 0], [1, 0, 0], [0, 0, 1]]
    }
};

/** Two-fold rotations about face diagonals, written `'` and `"` in Hall symbols. */
const DIAG: { [axis: string]: number[][] } = {
    'x\'': [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
    'x"': [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    'y\'': [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],
    'y"': [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
    'z\'': [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
    'z"': [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
};

const DIAG_DIR: { [axis: string]: number[] } = {
    'x\'': [0, 1, -1], 'x"': [0, 1, 1],
    'y\'': [1, 0, -1], 'y"': [1, 0, 1],
    'z\'': [1, -1, 0], 'z"': [1, 1, 0]
};

/** The three-fold axis about the body diagonal of the rhombohedral setting. */
const BODY3: number[][] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];

/** Hall translation letters, in 24ths. */
const LETTER_T: { [letter: string]: number[] } = {
    a: [12, 0, 0], b: [0, 12, 0], c: [0, 0, 12],
    n: [12, 12, 12], u: [6, 0, 0], v: [0, 6, 0], w: [0, 0, 6],
    d: [6, 6, 6]
};

/** Lattice centring vectors, in 24ths. */
const CENTER: { [letter: string]: number[][] } = {
    p: [[0, 0, 0]],
    a: [[0, 0, 0], [0, 12, 12]],
    b: [[0, 0, 0], [12, 0, 12]],
    c: [[0, 0, 0], [12, 12, 0]],
    i: [[0, 0, 0], [12, 12, 12]],
    f: [[0, 0, 0], [0, 12, 12], [12, 0, 12], [12, 12, 0]],
    r: [[0, 0, 0], [16, 8, 8], [8, 16, 16]]
};

function mod24(v: number): number {
    return ((v % DEN) + DEN) % DEN;
}

function cloneMatrix(m: number[][]): number[][] {
    return [m[0].slice(), m[1].slice(), m[2].slice()];
}

function negMatrix(m: number[][]): number[][] {
    return m.map(row => row.map(v => -v));
}

function opMul(a: IntOp, b: IntOp): IntOp {
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            R[i][j] = a.R[i][0] * b.R[0][j] + a.R[i][1] * b.R[1][j] + a.R[i][2] * b.R[2][j];
        }
    }
    const t = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        t[i] = mod24(a.R[i][0] * b.t[0] + a.R[i][1] * b.t[1] + a.R[i][2] * b.t[2] + a.t[i]);
    }
    return { R, t };
}

function opKey(op: IntOp): string {
    return op.R[0].join(',') + ';' + op.R[1].join(',') + ';' + op.R[2].join(',') + '|' + op.t.join(',');
}

interface Generator {
    op: IntOp;
    order: number;
    axis: string;
}

/**
 * One Hall generator: an order (`4`), an optional axis (`x`, `'`, `*`), an
 * optional improper sign (`-4`) and translation letters or fractions (`2_1`
 * is written `21`, `4 2vw`, `-4ud`). The axis a bare generator turns about is
 * part of the Hall grammar: the first is z, the second is x after a 2- or
 * 4-fold and the x-y diagonal after a 3- or 6-fold, and the third is the body
 * diagonal.
 */
function parseGenerator(token: string, prev: Generator | null, index: number): Generator {
    let i = 0;
    let improper = false;
    if (token.charAt(i) === '-') {
        improper = true;
        i++;
    }
    const orderChar = token.charAt(i);
    if (!/[1-6]/.test(orderChar)) throw new Error('no rotation order in ' + token);
    const order = parseInt(orderChar, 10);
    i++;

    let axisChars = '';
    while (i < token.length && 'xyz\'"*'.indexOf(token.charAt(i)) >= 0) {
        axisChars += token.charAt(i);
        i++;
    }
    const transChars = token.substring(i);

    let axis: string;
    if (axisChars === '') {
        if (order === 1) axis = 'z';
        else if (index === 0) axis = 'z';
        else if (index === 1) {
            const p = prev === null ? 1 : prev.order;
            if (p === 2 || p === 4) axis = 'x';
            else if (p === 3 || p === 6) axis = 'z\'';
            else throw new Error('no default second axis after order ' + p + ' in ' + token);
        } else if (index === 2) {
            axis = '*';
        } else {
            throw new Error('too many generators in ' + token);
        }
    } else if (axisChars === 'x' || axisChars === 'y' || axisChars === 'z' || axisChars === '*') {
        axis = axisChars;
    } else if (axisChars === '\'' || axisChars === '"') {
        if (prev === null || index === 0) throw new Error('diagonal axis first in ' + token);
        const base = prev.axis.charAt(0);
        if (base !== 'x' && base !== 'y' && base !== 'z') throw new Error('diagonal axis after ' + prev.axis);
        axis = base + axisChars;
    } else {
        throw new Error('unknown axis "' + axisChars + '" in ' + token);
    }

    let R: number[][];
    if (axis === '*') {
        if (order !== 3) throw new Error('body diagonal with order ' + order);
        R = cloneMatrix(BODY3);
    } else if (axis.length === 1) {
        const table = PRIN[axis];
        if (!table || !table[order]) throw new Error('order ' + order + ' about ' + axis);
        R = cloneMatrix(table[order]);
    } else {
        if (order !== 2) throw new Error('diagonal axis with order ' + order + ' in ' + token);
        R = cloneMatrix(DIAG[axis]);
    }
    if (improper) R = negMatrix(R);

    const dir = axis === '*' ? [1, 1, 1]
        : axis.length === 1
            ? (axis === 'x' ? [1, 0, 0] : axis === 'y' ? [0, 1, 0] : [0, 0, 1])
            : DIAG_DIR[axis];

    const t = [0, 0, 0];
    for (let k = 0; k < transChars.length; k++) {
        const ch = transChars.charAt(k);
        const letter = LETTER_T[ch];
        if (letter) {
            t[0] += letter[0];
            t[1] += letter[1];
            t[2] += letter[2];
        } else if (ch >= '1' && ch <= '5') {
            const step = parseInt(ch, 10) * DEN / order;
            if (step !== Math.round(step)) throw new Error('fraction ' + ch + '/' + order + ' in ' + token);
            t[0] += dir[0] * step;
            t[1] += dir[1] * step;
            t[2] += dir[2] * step;
        } else {
            throw new Error('unknown translation "' + ch + '" in ' + token);
        }
    }
    return { op: { R, t: [mod24(t[0]), mod24(t[1]), mod24(t[2])] }, order, axis };
}

/**
 * The operations a Hall symbol generates, in 24ths, identity first.
 *
 * A leading `-` on the lattice letter adds the inversion; an origin shift in
 * parentheses (twelfths of a cell edge, used by a few settings) re-origins the
 * whole group with `t' = t + (I - R)v`, which is what makes those settings'
 * operations agree with the coordinates the same file lists.
 */
function decodeHall(hall: string): IntOp[] {
    let text = hall.toLowerCase().replace(/_/g, ' ').trim();
    let shift = [0, 0, 0];
    const pm = /\(([^)]*)\)/.exec(text);
    if (pm) {
        const parts = pm[1].trim().split(/[\s,]+/).map(Number);
        if (parts.length !== 3 || parts.some(v => !Number.isInteger(v))) {
            throw new Error('origin shift "' + pm[1] + '"');
        }
        shift = parts.map(v => v * DEN / 12);
        text = text.replace(pm[0], ' ').trim();
    }

    const tokens = text.split(/\s+/).filter(t => t !== '');
    let latticeToken = tokens.length > 0 ? tokens.shift() as string : '';
    let centrosymmetric = false;
    if (latticeToken.charAt(0) === '-') {
        centrosymmetric = true;
        latticeToken = latticeToken.substring(1);
    }
    const lattice = latticeToken.toLowerCase();
    if (!CENTER[lattice]) throw new Error('unknown lattice "' + latticeToken + '"');

    const generators: IntOp[] = [];
    let prev: Generator | null = null;
    for (let idx = 0; idx < tokens.length; idx++) {
        const g = parseGenerator(tokens[idx], prev, idx);
        generators.push(g.op);
        prev = g;
    }

    const all: IntOp[] = [];
    const seen: { [key: string]: true } = {};
    function add(op: IntOp): boolean {
        const key = opKey(op);
        if (seen[key]) return false;
        seen[key] = true;
        all.push(op);
        return true;
    }

    add({ R: cloneMatrix(I3), t: [0, 0, 0] });
    const centerings = CENTER[lattice];
    for (const v of centerings) add({ R: cloneMatrix(I3), t: v.slice() });
    if (centrosymmetric) {
        for (const v of centerings) add({ R: negMatrix(I3), t: v.slice() });
    }
    for (const g of generators) {
        for (const v of centerings) {
            add({ R: cloneMatrix(g.R), t: [mod24(g.t[0] + v[0]), mod24(g.t[1] + v[1]), mod24(g.t[2] + v[2])] });
        }
    }

    const base = all.slice();
    for (let changed = true; changed;) {
        changed = false;
        const current = all.slice();
        for (const a of base) {
            for (const b of current) {
                if (add(opMul(a, b))) changed = true;
            }
        }
    }

    if (shift[0] || shift[1] || shift[2]) {
        for (const op of all) {
            for (let i = 0; i < 3; i++) {
                const d = shift[i] - (op.R[i][0] * shift[0] + op.R[i][1] * shift[1] + op.R[i][2] * shift[2]);
                op.t[i] = mod24(op.t[i] + d);
            }
        }
    }
    return all;
}

// ===== Names =====

interface SgRow {
    number: number;
    setting: string;
    symbol: string;
    hall: string;
    keys: string[];
}

function stripSettingSuffix(hm: string): string {
    const m = /\s*:\s*([12hr])\s*$/i.exec(hm);
    return m ? hm.substring(0, m.index).trim() : hm.trim();
}

function symbolTokens(hm: string): string[] {
    return hm.toLowerCase()
        .replace(/_/g, '')
        .split(/\s+/)
        .filter(t => t !== '');
}

function symbolKey(hm: string): string {
    return symbolTokens(hm).join('');
}

const ROWS: SgRow[] = [];
const INDEX = new Map<string, number[]>();

function addIndex(key: string, rowIndex: number): void {
    if (key === '') return;
    const list = INDEX.get(key);
    if (!list) {
        INDEX.set(key, [rowIndex]);
        return;
    }
    if (list.indexOf(rowIndex) < 0) list.push(rowIndex);
}

for (const line of SG_HALL_DATA) {
    const parts = line.split('|');
    if (parts.length !== 4) continue;
    const number = parseInt(parts[0], 10);
    if (!(number >= 1 && number <= 230)) continue;
    const setting = parts[1].trim();
    const hm = parts[2].trim();
    const hall = parts[3].trim();
    if (hm === '' || hall === '') continue;

    const base = stripSettingSuffix(hm);
    const tokens = symbolTokens(base);
    if (tokens.length < 2) continue;
    const keys = [symbolKey(base)];
    // Monoclinic HM symbols are written full (P 1 21/c 1); the short form the
    // rest of the world writes (P 21/c) is the same group with the redundant
    // 1 positions dropped. Only the unique-axis-b rows get that alias: a short
    // monoclinic name means the standard b setting, and the b position is the
    // only one whose reduction cannot collide with another setting's full
    // symbol (both P 1 21/c 1 and P 21/c 1 1 reduce to P 21/c).
    if (number >= 3 && number <= 15 && setting.charAt(0) === 'b') {
        const positions = tokens.slice(1);
        const real = positions.filter(t => t !== '1');
        if (real.length === 1 && positions.length > 1) keys.push(tokens[0] + real[0]);
    }

    ROWS.push({ number, setting, symbol: base, hall, keys });
    for (const key of keys) addIndex(key, ROWS.length - 1);
}

interface ParsedName {
    setting: string;
    keys: string[];
}

/**
 * The keys a file's spelling of a name may match, most specific first.
 *
 * Besides the spelling as written, two crystallographic reductions are
 * offered: a monoclinic full symbol can lose its redundant 1 positions, and a
 * full symbol of any system can lose the `2/` (or `2_1/`) that a short symbol
 * omits (P 4/m 2/m 2/m is P 4/m m m; F 4_1/d -3 2/m is F d -3 m). The
 * reductions are only fallbacks — an exact hit is always tried first, so a
 * name that already pins a setting is never reinterpreted.
 */
function derivedKeys(tokens: string[]): string[] {
    const out: string[] = [];
    if (tokens.length < 2) return out;
    const lattice = tokens[0];
    if (!/^[pabcifr]$/.test(lattice)) return out;
    const positions = tokens.slice(1);

    const real = positions.filter(t => t !== '1');
    if (real.length === 1 && positions.length > real.length) {
        out.push(lattice + real[0]);
    }

    if (positions.length === 3) {
        const reduced = positions.slice();
        const dropTwo = (t: string): string => {
            const m = /^2[0-9]?\/(.+)$/.exec(t);
            return m ? m[1] : t;
        };
        reduced[1] = dropTwo(reduced[1]);
        reduced[2] = dropTwo(reduced[2]);
        // The first position of a cubic full symbol carries the same 4/n or
        // 4_1/n; the presence of a -3 or 3 in the second position is what says
        // this is a cubic symbol rather than a tetragonal one (P 4/m m m must
        // not reduce to the orthorhombic P m m m).
        if (reduced[1] === '-3' || reduced[1] === '3') {
            const m = /^4[0-9]?\/(.+)$/.exec(reduced[0]);
            if (m) reduced[0] = m[1];
        }
        out.push(lattice + reduced.join(''));
    }
    return out;
}

/**
 * Parse a file's space-group name into lookup keys.
 *
 * Files write the same group as `P 21/c`, `P21/c`, `P 2₁/c` or `P 1 21/c 1`,
 * and some append the setting as `:2`, `:H` or an ICSD `Z`. All of those are
 * normalised here: unicode minus and overbars, subscripts, quotes, braces,
 * underscores and spacing disappear, and the setting suffix is separated.
 */
function parseName(raw: string): ParsedName | null {
    let s = String(raw === undefined || raw === null ? '' : raw);
    s = s.replace(/[\u2212\u2013\u2014]/g, '-');
    s = s.replace(/([0-9])[\u0304\u0305]/g, '-$1');
    s = s.replace(/[\u0304\u0305]/g, '-');
    s = s.replace(/[\u2080-\u2089]/g, ch => String(ch.charCodeAt(0) - 0x2080));
    s = s.replace(/\(\s*origin\s+choice\s*([12])\s*\)/gi, ':$1');
    s = s.replace(/\([^)]*\)/g, ' ');
    s = s.replace(/["'{}]/g, ' ');
    s = s.replace(/_/g, '');
    s = s.replace(/\s+/g, ' ').trim();
    s = s.toLowerCase();
    if (s === '') return null;

    let setting = '';
    const sm = /:\s*([12hr])\s*$/i.exec(s);
    if (sm) {
        setting = sm[1].toLowerCase();
        s = s.substring(0, sm.index).trim();
    } else {
        // ICSD spells origin choice 2 as a trailing Z.
        const zm = /\s+z$/i.exec(s);
        if (zm) {
            setting = '2';
            s = s.substring(0, zm.index).trim();
        }
    }

    let tokens = s.split(/\s+/).filter(t => t !== '');
    // `Fm-3m`: one token, but the leading lattice letter is always separable.
    if (tokens.length === 1 && /^[pabcifr]/.test(tokens[0])) {
        tokens = [tokens[0].charAt(0), tokens[0].substring(1)];
    }
    if (tokens.length < 2) return null;

    const keys: string[] = [tokens.join('')];
    for (const key of derivedKeys(tokens)) {
        if (keys.indexOf(key) < 0) keys.push(key);
    }
    return { setting, keys };
}

// ===== Choosing a setting =====

/**
 * Whether a cell is the primitive rhombohedral cell of an R group rather than
 * its hexagonal one. The distinction matters because the two settings use
 * different coordinates for the same crystal; the table carries both.
 */
function looksRhombohedral(cell: CellShape | undefined): boolean {
    if (!cell) return false;
    const lens = [cell.a, cell.b, cell.c];
    if (!(lens[0] > 0 && lens[1] > 0 && lens[2] > 0)) return false;
    const maxLen = Math.max(lens[0], lens[1], lens[2]);
    const minLen = Math.min(lens[0], lens[1], lens[2]);
    if ((maxLen - minLen) / maxLen > 0.02) return false;
    const angles = [cell.alpha, cell.beta, cell.gamma];
    if (Math.max(angles[0], angles[1], angles[2]) - Math.min(angles[0], angles[1], angles[2]) > 0.5) return false;
    return Math.abs(cell.alpha - 90) > 0.5;
}

/**
 * The origin-choice digit of a setting code. Orthorhombic settings of a group
 * with two origins concatenate the two (`2cab`), so the digit is the leading
 * character, not the whole code.
 */
function originChoiceOf(setting: string): string {
    const ch = setting.charAt(0);
    return ch === '1' || ch === '2' ? ch : '';
}

/** Whether a setting code answers a `:1`, `:2`, `:h` or `:r` request. */
function matchesSetting(rowSetting: string, wanted: string): boolean {
    if (rowSetting === wanted) return true;
    if (wanted === '1' || wanted === '2') return originChoiceOf(rowSetting) === wanted;
    return false;
}

function pickRow(rows: SgRow[], cell: CellShape | undefined): SgRow {
    if (rows.length === 1) return rows[0];

    const hex = rows.filter(r => r.setting === 'h');
    const rhom = rows.filter(r => r.setting === 'r');
    if (hex.length > 0 && rhom.length > 0) {
        return looksRhombohedral(cell) ? rhom[0] : hex[0];
    }

    // Origin choice 2 is the setting the International Tables and the CIF
    // dictionaries recommend; choice 1 is kept for files that ask for it.
    const origin2 = rows.filter(r => originChoiceOf(r.setting) === '2');
    if (origin2.length > 0) return origin2[0];

    // Monoclinic: unique axis b, cell choice 1; orthorhombic: the first row
    // (abc). That is the table's own order, so the first row is the default.
    return rows[0];
}

function settingLabelFor(row: SgRow, rows: SgRow[]): string {
    const origin = originChoiceOf(row.setting);
    const hasOriginChoice = rows.some(r => originChoiceOf(r.setting) === '1') &&
        rows.some(r => originChoiceOf(r.setting) === '2');
    if (hasOriginChoice) return 'origin choice ' + origin;
    if (row.setting === 'h') return 'hexagonal axes';
    if (row.setting === 'r') return 'rhombohedral axes';
    if (row.number >= 3 && row.number <= 15) {
        if (rows.length === 1) return '';
        const axis = row.setting.charAt(0);
        if (axis === 'a' || axis === 'b' || axis === 'c') {
            const choice = row.setting.substring(1);
            return 'unique axis ' + axis + (choice === '' ? '' : ', cell choice ' + choice);
        }
    }
    if (rows.length > 1 && row.setting !== '') return 'setting ' + row.setting;
    return '';
}

const OPERATION_CACHE = new Map<string, SymOp[]>();

function operationsFor(row: SgRow): SymOp[] {
    const cacheKey = row.number + '|' + row.setting + '|' + row.hall;
    const cached = OPERATION_CACHE.get(cacheKey);
    if (cached) return cached;
    const ops = decodeHall(row.hall).map(op => ({
        R: op.R.map(r => r.slice()),
        t: op.t.map(v => v / DEN)
    }));
    OPERATION_CACHE.set(cacheKey, ops);
    return ops;
}

// ===== Public surface =====

export interface CellShape {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
}

export interface SpaceGroupOps {
    /** The operations, in CIF spelling (`-x+1/2, y, -z+1/2`), identity first. */
    ops: string[];
    /** The space-group number, 1-230. */
    number: number;
    /** The table's Hermann-Mauguin symbol for the chosen setting. */
    symbol: string;
    /** The setting code (`b1`, `2`, `h`, `r`, `cab`, or empty). */
    setting: string;
    /** A human description of an origin or axis choice, when there is one. */
    settingLabel: string;
    /** Whether the name left a choice of setting open and a convention filled it. */
    assumed: boolean;
    /** Whether origin choice 2 was assumed because the file did not state one. */
    assumedOrigin: boolean;
}

/**
 * The operations a declared space-group name stands for, or null when the
 * name is not one this table knows. The `cell` is only consulted to choose
 * between the hexagonal and rhombohedral settings of an R group; everything
 * else is decided by the name's own settings suffix and the table's defaults.
 */
export function lookupSpaceGroupOps(name: string, cell?: CellShape): SpaceGroupOps | null {
    const parsed = parseName(name);
    if (!parsed) return null;

    let rows: SgRow[] | null = null;
    for (const key of parsed.keys) {
        const hit = INDEX.get(key);
        if (!hit || hit.length === 0) continue;
        const candidates = hit.map(i => ROWS[i]);
        // A key shared by different space groups is not enough to identify
        // one; the less specific keys that follow are tried instead.
        const number = candidates[0].number;
        if (candidates.some(r => r.number !== number)) continue;
        rows = candidates;
        break;
    }
    if (!rows) return null;

    if (parsed.setting !== '') {
        const filtered = rows.filter(r => matchesSetting(r.setting, parsed.setting));
        if (filtered.length > 0) rows = filtered;
    }

    const chosen = pickRow(rows, cell);
    const assumedOrigin = rows.length > 1 && originChoiceOf(chosen.setting) === '2' &&
        rows.some(r => originChoiceOf(r.setting) === '1');

    return {
        ops: operationsFor(chosen).map(symOpToString),
        number: chosen.number,
        symbol: chosen.symbol,
        setting: chosen.setting,
        settingLabel: settingLabelFor(chosen, rows),
        assumed: rows.length > 1,
        assumedOrigin
    };
}
