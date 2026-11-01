# Field-level descriptions of the 27 coupling-rule ingredients

This directory answers a referee request for "generator matrices, field conventions, ingredient-to-cell mappings and
exact rank-check outputs" for the coupling cells of the p10 manuscript (Sec. V, Table IV), including ingredients
that the manuscript specifies only by a candidate file or a construction description, such as the bordered
`[[16,0,8]]_4` state and the exceptional AME factors.

The frozen coupling certificates (`../certificates/AME_n_d.json`) contain each ingredient as a stabilizer state in
the verifier's canonical form: generator rows over a prime field `F_p`. Those rows are what was certified, and
they are authoritative. They are split out as `../ingredients/<id>.json` (27 distinct states), and
`../cells.json` maps every cell to its certificate and its two ingredients. For each of the 27 ingredients this
directory gives the object one level up, at field level. That is a Hermitian self-dual code `[I | A]` over
`F_{q^2}`, a classical code `C` whose minimal-support state is used, or a graph `Gamma`. Each comes with an
explicit field convention, structural parameters, the claimed property, and two independent checks:

* `check_fieldlevel.py` checks the claimed property **at field level**, using its own `GF(Q)` arithmetic built
  only from the printed modulus. The checks are exact rank conditions over `F_Q` on all column selections, plus a
  rebuild of the printed matrix from the structural parameters.
* `convert.py` maps each field-level object to canonical `F_p` rows with the documented converter (that of
  `qlib/ame.py`, which produced the certified rows). It confirms that these rows span **exactly** the certified
  `F_p` row space, by comparing reduced row echelon forms (`MATCH`). It also confirms the ingredient-to-cell
  mapping against the frozen certificates.

Result: 27/27 ingredients PASS the field-level check and 27/27 MATCH the certified rows. All 28 certificates
(26 new cells plus the 2 calibrations AME(5,10) and AME(11,10)) map to the stated ingredients, and every cell
satisfies the conditions of the coupling lemma.

## Files

| file | content |
|---|---|
| `<id>.json` (27) | field-level description of ingredient `<id>` (format below) |
| `check_fieldlevel.py` | field-level verification, independent of any conversion (Python 3.9 stdlib, own field arithmetic) |
| `convert.py` | converter field level -> canonical `F_p` rows, `MATCH` test against `../ingredients/<id>.json`, cell-mapping test against `../certificates/`, negative controls |
| `derive.py` | provenance: re-derives every `<id>.json` from the certified rows by inverting the converter, and asserts the structure stated in the construction notes. `--check` writes nothing and compares with the stored files |
| `run.sh` | runs `convert.py`, `check_fieldlevel.py` and `derive.py --check`, diffs each output with `expected/`, exits nonzero on any failure or difference |
| `expected/*.out` | the recorded outputs (`convert.out`, `check_fieldlevel.out`, `derive_check.out`) |

All scripts use only the Python 3.9 standard library. They import nothing from outside this directory
(`derive.py` imports `convert.py`), and they read only `../ingredients`, `../cells.json` and `../certificates`.
Runtime is about 3 s in total and memory about 15 MB.

## Conventions

### Field elements

`GF(Q)`, `Q = p^m`. A field element is an int `e = c_0 + c_1 p + ... + c_{m-1} p^{m-1}` (`0 <= c_i < p`), meaning
`c_0 + c_1 x + ... + c_{m-1} x^{m-1}` in `F_p[x]/(f)`. The modulus `f` is printed in every file as `field.modulus`,
the monic coefficient list `[c_0, ..., c_m]`, and as `field.modulus_str`. For a prime field (`m = 1`,
`modulus: null`) the elements are the ints `0..p-1` with arithmetic mod `p`. The moduli are those of the code
that produced the certified rows (`qlib` `GF`): the lexicographically first primitive monic polynomial, with the
coefficient vector `(c_0, ..., c_{m-1})` read as the base-`p` number `sum c_i p^i`, smallest first. Hence `x` (the int
`p`) is a primitive element in every non-prime field here. The fields used:

| field | modulus | used by |
|---|---|---|
| `F_16` | `x^4 + x + 1` (same as Table III of the paper) | `[[14,0,7]]_4`, `[[16,0,8]]_4`, RS AME(12,16) |
| `F_25` | `x^2 + x + 2` | AME(12,5), AME(14,5) |
| `F_27` | `x^3 + 2x + 1` | RS AME(12,27) |
| `F_49` | `x^2 + x + 3` | AME(12,7), AME(14,7) |
| `F_64` | `x^6 + x + 1` | AME(12,8) |
| `F_81` | `x^4 + x + 2` (the HZ field `F_3(g)`, `g^4 + g + 2 = 0`, with `g = x`) | AME(12,9), AME(14,9), AME(16,9) |
| `F_121` | `x^2 + x + 7` | AME(14,11), AME(16,11) |
| `F_169` | `x^2 + x + 2` | AME(16,13) |
| `F_p`, p prime | (none) | `F_2`: wheel, dodecacode. `F_3`: Golay. `F_5`, `F_11`, `F_13`, ... `F_31`: RS |

For Hermitian codes, `q` is the local dimension, `Q = q^2`, and `conj(z) = z^q`. The subfield `F_q` is the set
`{z : z^q = z}` of ints of `F_Q`.

### Kinds of field-level object

* **`hermitian_systematic`**: `A` (`n/2 x n/2` over `F_Q`, `Q = q^2`). The code is `C = rowspace [I | A]`,
  Hermitian self-dual (`<u, v> = sum u_i v_i^q`), which holds iff `A conj(A)^T = -I` (`-I = I` in characteristic 2).
  `C` defines an `F_q`-linear stabilizer state of `n` parties of dimension `q` [KKKS06, GBR04]. The state is
  `t`-uniform iff `d(C) >= t+1`, iff every `n-t` columns of `[I | A]` have rank `n/2`. A codeword supported on a
  `t`-set `T` exists iff the columns outside `T` have rank `< n/2`. AME(`n`, `q`) means `t = n/2`, that is, `C` is
  MDS `[n, n/2, n/2+1]`.
* **`css_minimal_support`**: `G` over `F_q`, `C = rowspace G`. The state is `|C> = |C|^{-1/2} sum_{c in C} |c>`,
  stabilized by `X(c)` (`c in C`) and `Z(h)` (`h in C^perp`, Euclidean dual). It is `t`-uniform iff `d(C) >= t+1`
  (every `n-t` columns of `G` have rank `dim C`) and `d(C^perp) >= t+1` (every `t` columns of `G` are
  independent). For `dim C = n/2` and `t = n/2` both conditions say that every `n/2` columns are independent.
* **`fq_graph`**: `Gamma`, symmetric with zero diagonal over `F_q` (here `q = 2`). The state is the graph state
  (eq. (1) of the paper), stabilized by `X(e_i) Z(Gamma_i)`. It is `t`-uniform iff `Gamma[A, A^c]` has full row
  rank `t` for every set `A` of `t` parties. This is Lemma 2 of the paper, stated there for `t = floor(n/2)`.

The dodecacode is given as an `fq_graph`, because its certified rows are literally `[I | Gamma]` with `Gamma` the
circulant graph on `Z_12` with connection set `{3, 5, 6, 7, 9}`. The file also gives the additive cyclic GF(4)
description (`additive_gf4`): the cyclic generator is `X_0 Z_3 Z_5 Z_6 Z_7 Z_9` (`x = 0x001`, `z = 0x2e8`), that is,
the GF(4) word `1 0 0 w 0 w w w 0 w 0 0` under `X(a)Z(b) -> a + b w`.

### Structural parameters (`structure`)

`check_fieldlevel.py` rebuilds the printed matrix from these parameters and requires equality.

| `type` | rule |
|---|---|
| `circulant` | `A[i][j] = first_row[(j - i) mod c]` |
| `bordered_circulant` | `A[0][0] = alpha`, `A[0][j] = beta`, `A[j][0] = gamma` (`1 <= j <= c`), `A[1+i][1+j] = core_first_row[(j - i) mod c]`. Optionally `core_exponents` with `core_first_row[t] = g^core_exponents[t]` (HZ parent) and `corner_relation` `alpha = -conj(sum core) beta / conj(beta)` (group circulants of `grp_orbits3.py`) |
| `block_circulant` | `A[3r + x][3s + y] = blocks_first_rows[r][s][(y - x) mod 3]` |
| `hz_kernel_submatrix` | `A = diag(Lambda) K[rows][cols] diag(D)`. `K` is the Hontarenko-Zyczkowski kernel `A_9` (formula in the file, from `work/mds/mdslib.hz_matrix`). `K` is rebuilt from the formula, compared with the printed `kernel_matrix`, and checked to satisfy `K conj(K)^T = -I`. Normalisation `D[0] = 1` |
| `explicit_labels` | the printed source matrix (Bevins-Bidav `A_5`). The label `e` means `(e mod 5) + (e div 5) r`, with `r` a root of the printed minimal polynomial `x^2 + x + 2`. Here `r = x` (int 5), so the labels equal the ints |
| `reed_solomon`, `doubly_extended_reed_solomon` | `G[i][j] = x_j^i` (`0^0 = 1`), `i < k`, for the evaluation points `x_j` (field ints `0, 1, ..., min(n,q)-1`). For `n = q+1`, an extra column `e_{k-1}` (point at infinity). `G_systematic` is the RREF of the same code |
| `systematic_bordered_circulant` | `G = [I_6 \| P]`, with `P` a bordered circulant (the Paley matrix of the extended ternary Golay code) |
| `wheel` | 5-cycle `0-1-2-3-4-0` plus hub `5` |

### Converter (field level -> canonical rows), `convert.py`

Canonical form: `n` parties, local dimension `d = p^{m'}`, `N = m' n` qudits of dimension `p`, and party `j` owns qudits
`j m' ... j m' + m' - 1`. There are `N` rows `(a | b)` in `F_p^{2N}` for `X(a) Z(b)`. The converter is copied from
`qlib/ame.py` without behavioural change. Its constants are printed in each file's `converter` block, and
`convert.py` asserts that they are the ones it uses.

* `canonical_from_fq(S, F, q)`: `S` holds rows `(a | b)` over `F_q`. For each row and each `beta_l` of the basis
  `beta = (1, g, ..., g^{k-1})` of `F_q / F_p`, where `g` is the smallest-int primitive element of `F_q`, `k = log_p q`,
  and `beta = (1)` when `q` is prime, it writes `beta_l (a | b)` with `X` coordinates in `beta` and `Z` coordinates in
  the trace-dual basis `gamma` (`Tr_{F_q/F_p}(beta_i gamma_j) = delta_ij`).
* `canonical_from_hermitian(G, F_Q, q)`: `w` is the smallest int with `{w, w^q}` an `F_q`-basis of `F_Q` (normal
  basis), and `xi` the smallest int outside `F_q`. Each row `g` of `G = [I | A]` gives the two `F_q` rows `(a | b)`
  of `lambda g`, `lambda in {1, xi}`, with `lambda g = w a + w^q b` coordinatewise. Then `canonical_from_fq`.
* `mds_state_rows(G, F)` (CSS): `S = [G | 0]` plus `[0 | H]`, `H` a basis of `C^perp`. Then `canonical_from_fq`.
* `canonical_from_graph(Gamma, p)`: rows `[e_i | Gamma_i]`.

`derive.py` inverts this. It maps digits back to `a_j = sum x_{j,i} beta_i`, `b_j = sum z_{j,i} gamma_i`. In the
Hermitian case it forms `v = w a + w^q b`, takes the `F_Q` span, and brings it to the RREF `[I | A]` (unique). In
the CSS case the `X` parts span `C` and the `Z` parts span `C^perp`, and both are checked. In the graph case it
takes the RREF over `F_p`. It then asserts the structure given in the construction notes. That is how the
`<id>.json` files were produced.

### JSON format

`id`, `n`, `local_dim`, `kind`, `field` (`Q`, `p`, `m`, `modulus`, `modulus_str`, `encoding`, and for Hermitian
codes `q` and `conjugation`), then the matrix (`A`, `G` or `Gamma`, int arrays in the encoding above) and
`structure`. Then `claimed_property` (`statement`, `uniformity_k`, `n`, `d`), `code` (`params` `[n, n/2, distance]_Q`,
`distance`, `duality`), `stabilizer_code` (`[[n, 0, distance]]_d`), `used_as`, `used_in` (cells), `converter`,
`source` (`script` as recorded in the certificate, `builder` = the script that builds the matrix, `note` = the
certificate's construction note, `board_posts`), `canonical_file`, and optional `remarks`.

## The 27 ingredients

| id | kind | field | field-level code | property | structure | used in AME(n,d) | source (builder; board posts) |
|---|---|---|---|---|---|---|---|
| `a_14_0_7_d4_hermitian_circulant` | hermitian_systematic | F_16 = F_2[x]/(x^4 + x + 1) | [14,7,7]_16 | 6-uniform | circulant, first row (0,1,2,15,6,2,7) | (13,20), (13,28), (13,36), (13,44), (13,52) | work/couple/g14_4.py 0 (hit 0 of work/couple/near14_4.py); #597 |
| `a_16_0_8_d4_bordered_circulant` | hermitian_systematic | F_16 = F_2[x]/(x^4 + x + 1) | [16,8,8]_16 | 7-uniform | bordered Z_7 circulant: alpha=2, beta=4, gamma=4, core (2,1,1,0,2,8,5) | (15,36), (15,44), (15,52) | work/couple/g16_4.py 0 (hit 0 of work/couple/near_bord.py 4 7 8); #616, #644, #653 |
| `a_dodecacode_n12_d2` | fq_graph | F_2 | [[12,0,6]]_2 | 5-uniform | circulant graph on Z_12, connection set {3,5,6,7,9} | (11,10), (11,14), (11,18), (11,22), (11,26), (11,34), (11,38), (11,46), (11,54), (11,58), (11,62) | work/couple/lib.py dodeca(); #586 |
| `a_ternary_golay_n12_d3` | css_minimal_support | F_3 | [12,6,6]_3 (Euclidean self-dual) | 5-uniform | [I_6 \| P], P bordered Paley (core (0,1,2,2,1)) | (11,15), (11,21), (11,24), (11,33), (11,39), (11,48), (11,51), (11,57) | work/couple/lib.py golay12(); #586 |
| `a_wheel_ame6_2` | fq_graph | F_2 | [[6,0,4]]_2 | AME(6,2) | wheel: 5-cycle + hub | (5,10) | work/couple/lib.py wheel6(); #586 |
| `b_ame12_11_rs` | css_minimal_support | F_11 | [12,6,7]_11 | AME(12,11) | doubly extended RS, points 0..10 + inf | (11,22), (11,33) | work/couple/lib.py ame_rs() (qlib rs_extended); #586 |
| `b_ame12_13_rs` | css_minimal_support | F_13 | [12,6,7]_13 | AME(12,13) | RS, points 0..11 | (11,26), (11,39) | lib.py ame_rs(); #586 |
| `b_ame12_16_rs` | css_minimal_support | F_16 = F_2[x]/(x^4 + x + 1) | [12,6,7]_16 | AME(12,16) | RS, points 0..11 | (11,48) | lib.py ame_rs(); #586 |
| `b_ame12_17_rs` | css_minimal_support | F_17 | [12,6,7]_17 | AME(12,17) | RS, points 0..11 | (11,34), (11,51) | lib.py ame_rs(); #586 |
| `b_ame12_19_rs` | css_minimal_support | F_19 | [12,6,7]_19 | AME(12,19) | RS, points 0..11 | (11,38), (11,57) | lib.py ame_rs(); #586 |
| `b_ame12_23_rs` | css_minimal_support | F_23 | [12,6,7]_23 | AME(12,23) | RS, points 0..11 | (11,46) | lib.py ame_rs(); #586 |
| `b_ame12_27_rs` | css_minimal_support | F_27 = F_3[x]/(x^3 + 2x + 1) | [12,6,7]_27 | AME(12,27) | RS, points 0..11 | (11,54) | lib.py ame_rs(); #586 |
| `b_ame12_29_rs` | css_minimal_support | F_29 | [12,6,7]_29 | AME(12,29) | RS, points 0..11 | (11,58) | lib.py ame_rs(); #586 |
| `b_ame12_31_rs` | css_minimal_support | F_31 | [12,6,7]_31 | AME(12,31) | RS, points 0..11 | (11,62) | lib.py ame_rs(); #586 |
| `b_ame12_5_bevins_bidav` | hermitian_systematic | F_25 = F_5[x]/(x^2 + x + 2) | [12,6,7]_25 | AME(12,5) | explicit (Bevins-Bidav A_5, printed labels) | (11,10), (11,15) | work/couple/lib.py ame12_5(); #597 |
| `b_ame12_7_circulant` | hermitian_systematic | F_49 = F_7[x]/(x^2 + x + 3) | [12,6,7]_49 | AME(12,7) | Z_6 circulant, first row (14,45,1,13,34,23) | (11,14), (11,21) | work/couple/run_gorb.py 7 6 0 (work/mds/grp_orbits3.py, first hit); #597 |
| `b_ame12_8_circulant` | hermitian_systematic | F_64 = F_2[x]/(x^6 + x + 1) | [12,6,7]_64 | AME(12,8) | 2x2 blocks of 3x3 circulants, first rows [[(45,46,45),(28,6,57)],[(42,51,43),(24,24,60)]] | (11,24) | work/couple/ame12_8.py 1; #610, #616 |
| `b_ame12_9_circulant` | hermitian_systematic | F_81 = F_3[x]/(x^4 + x + 2) | [12,6,7]_81 | AME(12,9) | bordered Z_5 circulant: alpha=13, beta=gamma=5, core (2,13,80,41,24) | (11,18) | run_gorb.py 9 5 1 (grp_orbits3.py, first hit); #597 |
| `b_ame14_11_circulant` | hermitian_systematic | F_121 = F_11[x]/(x^2 + x + 7) | [14,7,8]_121 | AME(14,11) | bordered Z_6 circulant: alpha=109, beta=gamma=27, core (23,51,79,34,40,41) | (13,44) | run_gorb.py 11 6 1 (grp_orbits3.py, first hit); #597 |
| `b_ame14_13_rs` | css_minimal_support | F_13 | [14,7,8]_13 | AME(14,13) | doubly extended RS, points 0..12 + inf | (13,52) | lib.py ame_rs(); #597 |
| `b_ame14_5_circulant` | hermitian_systematic | F_25 = F_5[x]/(x^2 + x + 2) | [14,7,8]_25 | AME(14,5) | bordered Z_6 circulant: alpha=1, beta=gamma=8, core (23,11,16,9,3,4) | (13,20) | run_gorb.py 5 6 1 (grp_orbits3.py, first hit); #597 |
| `b_ame14_7_circulant` | hermitian_systematic | F_49 = F_7[x]/(x^2 + x + 3) | [14,7,8]_49 | AME(14,7) | bordered Z_6 circulant: alpha=8, beta=gamma=2, core (4,36,20,38,2,26) | (13,28) | run_gorb.py 7 6 1 (grp_orbits3.py, first hit); #597 |
| `b_ame14_9_hz` | hermitian_systematic | F_81 = F_3[x]/(x^4 + x + 2) | [14,7,8]_81 | AME(14,9) | diag(Lambda) K[R][C] diag(D), K = HZ kernel A_9, R = C = (0,1,2,3,5,6,7) | (13,36) | work/couple/hzsub.py 9 7; #622, #624 |
| `b_ame16_11_circulant` | hermitian_systematic | F_121 = F_11[x]/(x^2 + x + 7) | [16,8,9]_121 | AME(16,11) | bordered Z_7 circulant: alpha=92, beta=gamma=1, core (58,4,21,30,109,118,38) | (15,44) | work/couple/ame16_11.py 0 (hit 0 of f121z7.py); #653 |
| `b_ame16_13_circulant` | hermitian_systematic | F_169 = F_13[x]/(x^2 + x + 2) | [16,8,9]_169 | AME(16,13) | bordered Z_7 circulant: alpha=12, beta=gamma=3, core (2,110,22,29,144,164,76) | (15,52) | run_gorb.py 13 7 1 (grp_orbits3.py, first hit); #614, #616 |
| `b_ame16_9_hz` | hermitian_systematic | F_81 = F_3[x]/(x^4 + x + 2) | [16,8,9]_81 | AME(16,9) | HZ parent B_7(1,1;[29,21,37,2,40,57,76]_g): bordered Z_7 circulant, alpha=beta=gamma=1, core (53,68,80,9,2,19,41) = g^(29,21,37,2,40,57,76), g = x = 3 | (15,36) | work/couple/ame16_9.py (arXiv:2609.11796 suppl.); #637, #644 |
| `b_ame6_5_rs` | css_minimal_support | F_5 | [6,3,4]_5 | AME(6,5) | doubly extended RS, points 0..4 + inf | (5,10) | lib.py ame_rs(); #586 |

Board posts are in `board/board.db` (table `posts`). The builder files named in the JSONs
(`work/couple/gorb/grp3_*.json`, `work/couple/ing/*.json`, `work/couple/near14_4_result.json`,
`work/couple/nearb_q4_c7_d8.json`) were cross-checked when these files were made: each one gives the certified row
space, and the stored hits are the printed parameters. That cross-check needs the repository and is not part of
`run.sh`.

## Cells (ingredient-to-cell mapping)

`../cells.json` lists 28 certificates: the 26 cells of Table IV and the 2 calibrations AME(5,10) and AME(11,10)
(MATCHES_RECORD). Cell AME(`2k+1`, `ab`) = coupling of `G_a` (`k`-uniform, `2k+2` parties, dimension `a`) with
`G_b` (AME(`2k+2`, `b`)), `a <= b`.

| cell | k | verdict | a | G_a | b | G_b | certificate |
|---|---|---|---|---|---|---|---|
| AME(5,10) | 2 | MATCHES_RECORD | 2 | `a_wheel_ame6_2` | 5 | `b_ame6_5_rs` | `certificates/AME_5_10.json` |
| AME(11,10) | 5 | MATCHES_RECORD | 2 | `a_dodecacode_n12_d2` | 5 | `b_ame12_5_bevins_bidav` | `certificates/AME_11_10.json` |
| AME(11,14) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 7 | `b_ame12_7_circulant` | `certificates/AME_11_14.json` |
| AME(11,15) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 5 | `b_ame12_5_bevins_bidav` | `certificates/AME_11_15.json` |
| AME(11,18) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 9 | `b_ame12_9_circulant` | `certificates/AME_11_18.json` |
| AME(11,21) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 7 | `b_ame12_7_circulant` | `certificates/AME_11_21.json` |
| AME(11,22) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 11 | `b_ame12_11_rs` | `certificates/AME_11_22.json` |
| AME(11,24) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 8 | `b_ame12_8_circulant` | `certificates/AME_11_24.json` |
| AME(11,26) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 13 | `b_ame12_13_rs` | `certificates/AME_11_26.json` |
| AME(11,33) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 11 | `b_ame12_11_rs` | `certificates/AME_11_33.json` |
| AME(11,34) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 17 | `b_ame12_17_rs` | `certificates/AME_11_34.json` |
| AME(11,38) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 19 | `b_ame12_19_rs` | `certificates/AME_11_38.json` |
| AME(11,39) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 13 | `b_ame12_13_rs` | `certificates/AME_11_39.json` |
| AME(11,46) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 23 | `b_ame12_23_rs` | `certificates/AME_11_46.json` |
| AME(11,48) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 16 | `b_ame12_16_rs` | `certificates/AME_11_48.json` |
| AME(11,51) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 17 | `b_ame12_17_rs` | `certificates/AME_11_51.json` |
| AME(11,54) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 27 | `b_ame12_27_rs` | `certificates/AME_11_54.json` |
| AME(11,57) | 5 | IMPROVEMENT | 3 | `a_ternary_golay_n12_d3` | 19 | `b_ame12_19_rs` | `certificates/AME_11_57.json` |
| AME(11,58) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 29 | `b_ame12_29_rs` | `certificates/AME_11_58.json` |
| AME(11,62) | 5 | IMPROVEMENT | 2 | `a_dodecacode_n12_d2` | 31 | `b_ame12_31_rs` | `certificates/AME_11_62.json` |
| AME(13,20) | 6 | IMPROVEMENT | 4 | `a_14_0_7_d4_hermitian_circulant` | 5 | `b_ame14_5_circulant` | `certificates/AME_13_20.json` |
| AME(13,28) | 6 | IMPROVEMENT | 4 | `a_14_0_7_d4_hermitian_circulant` | 7 | `b_ame14_7_circulant` | `certificates/AME_13_28.json` |
| AME(13,36) | 6 | IMPROVEMENT | 4 | `a_14_0_7_d4_hermitian_circulant` | 9 | `b_ame14_9_hz` | `certificates/AME_13_36.json` |
| AME(13,44) | 6 | IMPROVEMENT | 4 | `a_14_0_7_d4_hermitian_circulant` | 11 | `b_ame14_11_circulant` | `certificates/AME_13_44.json` |
| AME(13,52) | 6 | IMPROVEMENT | 4 | `a_14_0_7_d4_hermitian_circulant` | 13 | `b_ame14_13_rs` | `certificates/AME_13_52.json` |
| AME(15,36) | 7 | IMPROVEMENT | 4 | `a_16_0_8_d4_bordered_circulant` | 9 | `b_ame16_9_hz` | `certificates/AME_15_36.json` |
| AME(15,44) | 7 | IMPROVEMENT | 4 | `a_16_0_8_d4_bordered_circulant` | 11 | `b_ame16_11_circulant` | `certificates/AME_15_44.json` |
| AME(15,52) | 7 | IMPROVEMENT | 4 | `a_16_0_8_d4_bordered_circulant` | 13 | `b_ame16_13_circulant` | `certificates/AME_15_52.json` |

## Running

```bash
cd release/coupling/fieldlevel
sh run.sh                         # about 3 s; prints ALL OK and exits 0, or FAILED and exits 1
python3 check_fieldlevel.py       # field-level checks, negative controls, per-cell summary
python3 convert.py                # MATCH per ingredient, cell mapping, negative controls
python3 derive.py --check         # stored JSON == re-derivation from the certified rows
```

`run.sh` uses `$PYTHON` if set, otherwise `python3`, with `-B` so that nothing is written here. It works on any copy
of `release/coupling/`, because it needs only `fieldlevel/`, `ingredients/`, `cells.json` and `certificates/`. The
full expected outputs are in `expected/`. The summary lines are:

```
check_fieldlevel.py:  ingredients: 27 PASS, 0 FAIL      (negative controls: 4 x FAIL as expected)   cells: 28 OK, 0 FAIL
convert.py:           ingredients: 27 MATCH, 0 MISMATCH  cells: 28 OK, 0 FAIL   (negative controls: 5 x MISMATCH as expected)
derive.py --check:    27 x "same"
```

Exact rank-check numbers per ingredient, from `expected/check_fieldlevel.out`:

| ingredient | condition | column sets (or vertex sets) checked | of full rank |
|---|---|---|---|
| AME(n, b), n = 6 / 12 / 14 / 16 (1 / 13 / 5 / 3 ingredients) | every n/2 columns of `[I\|A]` resp. `G` independent | 20 / 924 / 3432 / 12870 | all |
| `[[14,0,7]]_4` circulant, 6-uniform | every 8 columns of `[I\|A]` rank 7 | 3003 | 3003 |
| (exactness, d = 7) | 7-column sets of rank 7 | 3432 | 3292 (140 deficient = 70 bad 7\|7 cuts) |
| `[[16,0,8]]_4` bordered circulant, 7-uniform | every 9 columns of `[I\|A]` rank 8 | 11440 | 11440 |
| (exactness, d = 8) | 8-column sets of rank 8 | 12870 | 12296 (574 deficient = 287 bad 8\|8 cuts) |
| dodecacode, 5-uniform | `Gamma[A, A^c]` rank 5, \|A\| = 5 | 792 | 792 |
| (exactness) | \|A\| = 6 | 924 | 540 |
| ternary Golay, 5-uniform | 7 columns of G rank 6; 5 columns independent | 792; 792 | 792; 792 |
| (exactness, d = 6) | 6-column sets of rank 6 | 924 | 792 |
| wheel AME(6,2) | `Gamma[A, A^c]` rank 3, \|A\| = 3 | 20 | 20 |

In addition, every Hermitian ingredient satisfies `A conj(A)^T = -I`, and every printed matrix is rebuilt
exactly from its structural parameters. For AME(16,9) that includes `g^(29,21,37,2,40,57,76)` with `g = x` of
order 80. For AME(14,9) it includes the HZ kernel `A_9` from its formula, unitary, and the factorisation
`A = diag(Lambda) K[R][C] diag(D)`.

The deficient counts for `[[14,0,7]]_4` and `[[16,0,8]]_4` reproduce the numbers on the board: 70 of 1716 bad
7|7 cuts (#597) and 287 of 6435 bad 8|8 cuts (#616). A bad cut is rank-deficient on both sides.

## Remarks and limits

* **Structure versus notes.** Every structural statement in the construction notes was confirmed on the certified
  rows by `derive.py`, which asserts it, and by `check_fieldlevel.py`, which rebuilds it.
  * `[[14,0,7]]_4` is the circulant `(0,1,2,15,6,2,7)` printed in the paper.
  * The note of `[[16,0,8]]_4` prints `[2, 4, 4, 2, 1, 1, 0, 2, 8, 5]`. This list is `[alpha, beta, gamma, core_first_row]`
    in the format of `near_bord.py`. Its row border and column border are both 4, and the corner is 2.
  * The `b` ingredients from `grp_orbits3.py` have the symmetric border `beta = gamma`.
  * AME(16,11) has border 1, as on the board (#653).
  * The HZ parent exponents are exactly `[29,21,37,2,40,57,76]`.
  * The Bevins-Bidav labels coincide with the ints of this file.
  * The RS ingredients are qlib's `rs_extended`. That is the doubly extended code exactly when `n = q+1`: AME(6,5),
    AME(12,11), AME(14,13).

  No ingredient differs from its note.
* **`source.script`.** The ingredient files record the script of the first certificate in which the state occurs.
  For `a_16_0_8_d4_bordered_circulant` this is `work/couple/ame16_9.py` (certificate AME(15,36)). The matrix itself
  is built by `work/couple/g16_4.py 0`, from hit 0 of `near_bord.py 4 7 8`; `source.builder` records this.
* **Normalisation of AME(14,9).** `Lambda` and `D` are determined only up to `Lambda -> c Lambda`, `D -> D/c`. The
  file uses `D[0] = 1`. `hzsub.py` may have used another representative. Either way `A` is the same.
* **AME(12,8)** is printed as a block circulant. The Fourier parametrisation used in the search (`ame12_8.py`: each
  of the 3 Fourier modes is a 2x2 unitary matrix over `F_64`) is not reproduced; it is not needed to define or check
  the code.
* **What the checks rely on.** `check_fieldlevel.py` uses the standard equivalences stated above: Hermitian
  self-dual `F_{q^2}`-linear code <-> `F_q`-linear stabilizer state, `t`-uniform <-> distance `>= t+1`, CSS
  minimal-support states, and Lemma 2 for graph states. `convert.py` ties each field-level object to the certified
  `F_p` rows, which the frozen verifier checked directly (Lemma 1 on all cuts, coupling type). Neither script
  re-proves the coupling lemma (Lemma 3 of the paper). The cell summary only checks its hypotheses: `n = 2k+1`,
  ingredients on `2k+2` parties, `G_a` `k`-uniform, `G_b` AME, `a <= b`, `d = ab`.
