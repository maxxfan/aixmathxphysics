# Coupling certificates (paper Sec. V, Lemma 3, Table IV)

**Claim.** The 26 cells AME(11,d) for d = 14, 15, 18, 21, 22, 24, 26, 33, 34, 38, 39, 46, 48, 51, 54, 57, 58, 62,
AME(13,d) for d = 20, 28, 36, 44, 52, and AME(15,d) for d = 36, 44, 52 exist. Each follows from Lemma 3 (coupling): if
G_a is a k-uniform state of 2k+2 parties of local dimension a, G_b is an AME(2k+2, b) state and a <= b, then
psi = sum_(t<a) (<t|_last G_a) (x) (<t|_last G_b) is AME(2k+1, ab). The lemma holds for arbitrary states; its proof is
in the paper (and in `verifier/verify.py` `check_coupling` of the repository). A certificate is a pair of explicit
stabilizer states; what must be checked is that G_a is k-uniform and G_b is AME(2k+2, b). Two further certificates,
AME(5,10) and AME(11,10), are calibrations for cells that were already known.

## Files

| path | content |
|---|---|
| `certificates/AME_n_d.json` | the 28 frozen coupling certificates exactly as submitted to the project verifier: `{"type": "coupling", "k", "a", "b"}`, a and b stabilizer states in canonical F_p form |
| `ingredients/<id>.json` | the 27 distinct ingredient states (identical row spaces are identified), canonical F_p form, with `used_in` |
| `fieldlevel/` | each ingredient at field level: field, modulus, encoding, matrix (A of [I\|A] for Hermitian codes, Gamma for graphs, G for CSS states), structural parameters, code parameters, the claimed property, sources; a field-level checker, the converter to canonical form and a MATCH test (see `fieldlevel/README.md`) |
| `cells.json` | the cell -> certificate -> ingredient map (also below) |
| `check_mapping.py` | checks that each certificate's a and b span the same F_p row spaces as the named ingredient files, and the shapes (n = 2k+1, both ingredients on 2k+2 parties, d = ab, a <= b) |
| `negative/` | four negative controls (a regrouped ingredient a; a and b exchanged; two valid ingredients with a = 11 > b = 7, so that only a <= b fails; b replaced by a GHZ state): the checker must reject each |
| `run.sh` | replay; `expected/` holds the expected outputs |

**Canonical form** (what is certified): n parties of local dimension d = p^m; N = m*n qudits of prime dimension p, party
j owns qudits j*m .. j*m+m-1; N generator rows (a | b) in F_p^(2N) for the Pauli operators X(a)Z(b); rows are base-p
digit strings or integer lists. For F_q-linear and F_(q^2)-Hermitian constructions the X coordinates are taken in the
polynomial basis of GF(q) over F_p and the Z coordinates in its trace-dual basis (`fieldlevel/README.md` gives the
constants per ingredient).

**Exact checks.** `../checker/ame_check.py` (standalone, shares no code with the project verifier) checks isotropy,
rank N and, for every set B of n-k parties, that the N x 2m(n-k) column block of B has rank N over F_p (Lemma 1 of the
paper); an exhaustive C enumerator (`../checker/cutrank.c`) plus a seeded Python sample, or `--python-only` for a fully
exhaustive Python run. The counts it prints for every ingredient are in `expected/coupling_check_C.txt`.

## Replay

```bash
./run.sh          # about 15 s: mapping, 28 certificates, negative controls, field-level checks
./run.sh full     # also every certificate exhaustively in pure Python (about 70 s)
```

## Cell -> ingredients

| cell | k | ingredient a (k-uniform) | ingredient b (AME(2k+2, b)) | certificate |
|---|---|---|---|---|
| AME(5,10) (calibration, known cell) | 2 | `a_wheel_ame6_2` (dim 2) | `b_ame6_5_rs` (dim 5) | `certificates/AME_5_10.json` |
| AME(11,10) (calibration, known cell) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_5_bevins_bidav` (dim 5) | `certificates/AME_11_10.json` |
| AME(11,14) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_7_circulant` (dim 7) | `certificates/AME_11_14.json` |
| AME(11,15) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_5_bevins_bidav` (dim 5) | `certificates/AME_11_15.json` |
| AME(11,18) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_9_circulant` (dim 9) | `certificates/AME_11_18.json` |
| AME(11,21) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_7_circulant` (dim 7) | `certificates/AME_11_21.json` |
| AME(11,22) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_11_rs` (dim 11) | `certificates/AME_11_22.json` |
| AME(11,24) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_8_circulant` (dim 8) | `certificates/AME_11_24.json` |
| AME(11,26) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_13_rs` (dim 13) | `certificates/AME_11_26.json` |
| AME(11,33) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_11_rs` (dim 11) | `certificates/AME_11_33.json` |
| AME(11,34) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_17_rs` (dim 17) | `certificates/AME_11_34.json` |
| AME(11,38) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_19_rs` (dim 19) | `certificates/AME_11_38.json` |
| AME(11,39) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_13_rs` (dim 13) | `certificates/AME_11_39.json` |
| AME(11,46) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_23_rs` (dim 23) | `certificates/AME_11_46.json` |
| AME(11,48) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_16_rs` (dim 16) | `certificates/AME_11_48.json` |
| AME(11,51) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_17_rs` (dim 17) | `certificates/AME_11_51.json` |
| AME(11,54) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_27_rs` (dim 27) | `certificates/AME_11_54.json` |
| AME(11,57) | 5 | `a_ternary_golay_n12_d3` (dim 3) | `b_ame12_19_rs` (dim 19) | `certificates/AME_11_57.json` |
| AME(11,58) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_29_rs` (dim 29) | `certificates/AME_11_58.json` |
| AME(11,62) | 5 | `a_dodecacode_n12_d2` (dim 2) | `b_ame12_31_rs` (dim 31) | `certificates/AME_11_62.json` |
| AME(13,20) | 6 | `a_14_0_7_d4_hermitian_circulant` (dim 4) | `b_ame14_5_circulant` (dim 5) | `certificates/AME_13_20.json` |
| AME(13,28) | 6 | `a_14_0_7_d4_hermitian_circulant` (dim 4) | `b_ame14_7_circulant` (dim 7) | `certificates/AME_13_28.json` |
| AME(13,36) | 6 | `a_14_0_7_d4_hermitian_circulant` (dim 4) | `b_ame14_9_hz` (dim 9) | `certificates/AME_13_36.json` |
| AME(13,44) | 6 | `a_14_0_7_d4_hermitian_circulant` (dim 4) | `b_ame14_11_circulant` (dim 11) | `certificates/AME_13_44.json` |
| AME(13,52) | 6 | `a_14_0_7_d4_hermitian_circulant` (dim 4) | `b_ame14_13_rs` (dim 13) | `certificates/AME_13_52.json` |
| AME(15,36) | 7 | `a_16_0_8_d4_bordered_circulant` (dim 4) | `b_ame16_9_hz` (dim 9) | `certificates/AME_15_36.json` |
| AME(15,44) | 7 | `a_16_0_8_d4_bordered_circulant` (dim 4) | `b_ame16_11_circulant` (dim 11) | `certificates/AME_15_44.json` |
| AME(15,52) | 7 | `a_16_0_8_d4_bordered_circulant` (dim 4) | `b_ame16_13_circulant` (dim 13) | `certificates/AME_15_52.json` |

## Ingredients (details in `fieldlevel/README.md`)

The k-uniform ingredients are the dodecacode [[12,0,6]]_2 (a = 2, n = 11 cells), the extended ternary Golay state
[[12,0,6]]_3 (a = 3), a 6-uniform [[14,0,7]]_4 from the Hermitian self-dual circulant code [I | circ(0,1,2,15,6,2,7)]
over F_16 (n = 13 cells), and a 7-uniform [[16,0,8]]_4 from a bordered circulant Hermitian self-dual [16,8,8]_16 code
(corner 2, borders 4, core first row (2,1,1,0,2,8,5); n = 15 cells). The AME ingredients are Reed-Solomon minimal-support
states (AME(12,b) for b = 11, 13, 16, 17, 19, 23, 27, 29, 31; AME(14,13); AME(6,5)) and Hermitian self-dual MDS codes:
AME(12,5) (Bevins-Bidav), AME(12,7) (Z_6 circulant), AME(12,8) (Z_3-block circulant), AME(12,9) (bordered Z_5 circulant),
AME(14,5/7/11) (bordered Z_6 circulants), AME(14,9) (a submatrix of the Hontarenko-Zyczkowski kernel), AME(16,9) (the
Hontarenko-Zyczkowski parent M_(16,9)), AME(16,11) and AME(16,13) (bordered Z_7 circulants).
