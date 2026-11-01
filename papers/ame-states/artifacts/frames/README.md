# AME(8,4) frame searches: replication package

This folder contains everything needed to replay the exhaustive frame searches of Sec. VI B: the reduction and why it is
exhaustive, all ten gauge representatives (from both implementations, with the correspondence between their class numbers),
the enumeration code of BOTH independent implementations, the positive controls, and SHA-256 "count hashes" of every result.
`./run.sh` recompiles and reruns everything from source and compares each result with the recorded one.

## 1. Claim

**There is no F_2-linear (qubit-pair) stabilizer AME(8,4) state.** Equivalently: no 16-dimensional isotropic subspace
S of (F_2^4)^8, with the pairwise symplectic form, has every nonzero element supported on at least 5 of the 8 parties.
Equivalently again: there is no additive (F_2-linear) [[8,0,5]]_4 stabilizer code. Each party is two qubits, i.e. F_2^4 with a
nondegenerate alternating form. All such forms are GL(4,2)-equivalent, so the F_4 Weyl-Heisenberg (F_4-linear) case is
included. No symmetry is assumed.

The result is computer-assisted. It **confirms** Ball, Moreno and Simoens, *Stabiliser codes over fields of even order*,
arXiv:2401.06618 (v2, abstract and Sec. 5.2: no [[7,1,4]]_4 and no [[8,0,5]]_4). They work in the same model (additive
stabilisers, F_2^h per party, h = 2) and prove it by a different, geometric method. Our search is an independent computational
confirmation of their result. It is not the first proof.

Not covered: Z_4-labelled stabilizers and non-stabilizer states. The project treats the Z_4 case separately; non-stabilizer
AME(8,4) states remain open.

## 2. Quick start

```
./run.sh              # 4 parallel processes (JOBS=4, default)
JOBS=8 ./run.sh       # more parallelism if you have the cores
```
The whole replay is about 24 CPU-minutes of work. With JOBS=4 on an otherwise idle machine, expect roughly 7 minutes (an
estimate). The measured times, on a 10-core Apple M4 that was heavily loaded by other jobs, were 12.6 and 14.2 minutes with
JOBS=4; see Sec. 7. JOBS=1 should take about 25 minutes (estimate).

Requirements: bash, a C compiler invoked as `cc` (or set `CC`), and `python3` >= 3.9 (or set `PYTHON`). Python uses only
the standard library. There is no network access and no third-party package. Every process stays below about 130 MB of RAM
(the largest is `classcount.py` at about 127 MB; the C searches use a few MB).

`run.sh` changes into its own directory. It creates `build/` (binaries, `cc -O2`) and `out/` (logs, raw outputs, canonical
results, `out/SUMMARY.txt`), prints one `PASS`/`FAIL` line per check, and ends with `OVERALL: PASS` or `OVERALL: FAIL`.
It exits nonzero if any check fails. Rerunning deletes `build/` and `out/` first.

## 3. The reduction and why it is exhaustive

The full argument is in [PROOF.md](PROOF.md), copied unchanged from the project. The auditor's section-by-section review and
independent re-derivation are in [REVIEW.md](REVIEW.md). Here is a complete summary.

**Notation.** Parties 0..7. A = {0,1,2,3}, B = {4,5,6,7}. Each party has P = F_2^4 with coordinates (x1,x2,z1,z2) and the
form J(a,b) = x·z' + z·x'. W_l (l = 0..3) is the copy of P belonging to B-party 4+l. A_j (j = 0..3) is the copy of P
belonging to A-party j. Sp(4,2) is the isometry group of (P,J), of order 720.

1. **Normal form.** Let S be a stabilizer AME(8,4), i.e. a 16-dimensional isotropic subspace with no nonzero element of
   party weight <= 4. Restriction S -> P^A has kernel {elements supported on the 4 B-parties} = 0, and 16 = dim P^A, so it is
   a bijection. Hence S = {(u, Mu)} with M: P^A -> P^B linear and bijective. Isotropy says that M is an isometry.
2. **Frames.** Put V_j := M(A_j). Then V_0..V_3 are mutually orthogonal, nondegenerate 4-dimensional subspaces with
   P^B = V_0 + V_1 + V_2 + V_3. Conversely, every such frame together with any isometries E_j: P -> V_j defines an isotropic
   S. Different choices of E_j differ by local Cliffords on the A-parties. So, up to local Cliffords on A, S *is* its frame.
3. **AME condition (PROOF.md Lemma 2; Lemma 4 of the paper).** S is AME(8,4) iff V_J ∩ W_L = 0 for all J, L ⊂ {0..3} with |J| + |L| = 4, where
   V_J = ⊕_{j∈J} V_j and W_L = ⊕_{l∈L} W_l. The reason: v ∈ V_J ∩ W_L exactly when the element (M^{-1}v, v) has A-support
   inside J and B-support inside L. There are 68 such conditions; the duality (J,L) -> (J^c,L^c) reduces them to 34
   independent ones.
4. **Graph form.** The conditions with (|J|,|L|) = (1,3) say that every projection π_l|V_j is bijective. Therefore
   V_j = {(w, θ_{j,1}w, θ_{j,2}w, θ_{j,3}w)} with θ_{j,l} ∈ GL(4,2). Orthogonality V_i ⊥ V_j reads
   Σ_l θ_{i,l}^T J θ_{j,l} = 0, which is **linear** in θ_j.
5. **Gauge.** The group G = (Sp(4,2)^4 ⋊ S_4 on P^B: local Cliffords and permutations of the B-parties) × (S_4 relabelling
   the V_j) preserves ω, orthogonality, nondegeneracy and the family of conditions in item 3. Every solution frame contains
   some member, which can be relabelled V_0. It therefore suffices to fix V_0 to one representative of each G-class and
   enumerate **all ordered** completions (V_1, V_2, V_3).
6. **Classification of V_0.** Let Q_l be the pullback of J along π_l|V_0, l = 0..3. Each Q_l is one of the 28 nondegenerate
   alternating forms on F_2^4, and Σ Q_l must be nondegenerate. If two tuples (Q_0..Q_3) lie in the same orbit under
   GL(4,2) (simultaneous congruence) × S_4 (permuting the tuple), the corresponding V_0 are G-equivalent: D_l := π'_{σ(l)} h
   (π_l|V_0)^{-1} is an isometry. Only this direction is needed for exhaustiveness (REVIEW.md 5.18). There are 268,800
   admissible tuples out of 28^4 = 614,656, and they fall into **10 orbits** of sizes 80640, 60480, 40320, 20160 (three times),
   13440, 10080, 2016 and 1344. The orbit sizes sum to 268,800; `v0reps.py` asserts this, so no V_0 is missed. The second
   implementation derives the same 10 classes from its own code, and `classcount.py` recounts the orbit sizes with all 20160
   elements of GL(4,2).
7. **Enumeration.** For each of the 10 representatives, all ordered (V_1, V_2, V_3) are enumerated. Only **necessary**
   conditions are used for pruning, and every condition of item 3 is tested at some level. A frame surviving all levels is
   therefore an AME frame, and a frame that is pruned violates a condition. Both implementations find **0** surviving frames
   in every class.

**Why this is exhaustive.** Every stabilizer AME(8,4) would give a frame (items 1-2) satisfying all 68 conditions (item 3).
After a gauge transformation (item 5), its V_0 equals one of the 10 representatives (item 6), and its (V_1, V_2, V_3) is
then one of the enumerated completions (item 7). The enumeration reports none.

### The two implementations

| | impl1 (`impl1_classify/`, author "classify") | impl2 (`impl2_auditor/`, the auditor, written independently) |
|---|---|---|
| V_0 classes | `v0reps.py`: flood fill of the 268,800 4-tuples under the 12 elementary transvections of GL(4,2) plus S_4; representative normalised to Q_0 = J | `vclasses8.py`: triples (Q_1,Q_2,Q_3) with Q_0 = J fixed (9600 admissible); orbits under all 720 elements of Sp(4,2), S_3 and the base-party swap; also the 72 orbits under Sp(4,2) alone |
| parametrisation | V_j as graphs, matrices θ_{j,l} ∈ GL(4,2) | V_j as 4-dimensional subspaces of F_2^16, given by a basis (16-bit vectors) |
| level j = 1, 2 | Gaussian elimination of the linear orthogonality system: an affine space of dimension 32 (j = 1) or 16 (j = 2) for every node. All 2^32 or 2^16 points are visited in Gray-code order and tested for θ invertible, V_j nondegenerate, and every condition with j ∈ J ⊂ {0..j} | candidate vectors in (V_0 + ... + V_{j-1})^⊥ with all 4 nibbles nonzero and party weight >= 3 against every element of the earlier V_i; every 4-dimensional subspace inside the candidate set is listed by depth-first search over lexicographically least bases (`enumtest.c` proves this enumerator exact); then a nondegeneracy test |
| level 3 | V_3 is forced by a square linear system (at most one solution); all conditions with 3 ∈ J are tested | V_3 = (V_0+V_1+V_2)^⊥ is forced; its single-element and pair conditions are checked directly. The (3,1) conditions are the V_3 single conditions by duality |
| condition test | rank of stacked basis vectors masked to W_{L^c} | party weights of vectors (number of nonzero nibbles) |

The two implementations share the reduction (items 1-6), which the auditor re-derived in writing (REVIEW.md 5.12 and 5.13).
They share no code. Their level-2 definitions differ: impl1 level 2 includes the four conditions with J = {0,1,2}, impl2 level 2
does not. The level-2 numbers are therefore compared using impl1 with `JMASKS=124`, which drops exactly those four conditions.

## 4. The ten gauge representatives

The two class lists are numbered differently. impl1 sorts by decreasing orbit size, breaking ties by discovery order. impl2
sorts by the smallest triple code (a_1·28 + a_2)·28 + a_3 in the orbit. Orbit sizes also use different units: impl1 counts
4-tuples (Q_0..Q_3), impl2 counts triples with Q_0 = J, so impl1 = 28 × impl2.

The correspondence below is established by a gauge **invariant**, not by matching counts. `tools/crosswalk.py`, written for
this release and sharing no code with either implementation, computes the GL(4,2) × S_4 congruence class of (Q_0..Q_3) for
every representative. It agrees with the count-based matching recorded during the audit.

Packing conventions:
* impl1 stores a 4×4 matrix over F_2 as an integer: entry (r,c) is bit 4r+c, and V_0 = {(w, ψ_1 w, ψ_2 w, ψ_3 w)}.
* impl2 stores V_0 as 4 basis vectors. Each is 16 bits: party l is nibble l (bits 4l..4l+3), and within a nibble bits 0..3
  are (x1, x2, z1, z2).
* The stabiliser of V_0 in Sp(4,2)^4 ⋊ S_4 has order 483840 / (impl1 orbit size).

| impl1 class | impl1 orbit (4-tuples) | stabiliser | ψ_1 ψ_2 ψ_3 (packed) | impl2 rep | impl2 orbit (triples) | impl2 V_0 basis (packed) | impl2 Sp-only reps in class |
|---|---|---|---|---|---|---|---|
| 0 | 80640 | 6 | 5673 4804 5160 | 5 | 2880 | 4369 9282 33316 27016 | 12 |
| 1 | 60480 | 8 | 5673 4684 4740 | 2 | 2160 | 4369 17474 8996 39048 | 6 |
| 2 | 40320 | 12 | 5673 5160 4685 | 4 | 1440 | 4369 17474 29220 35208 | 4 |
| 3 | 20160 | 24 | 4680 5673 4804 | 1 | 720 | 4369 9282 33316 59528 | 12 |
| 4 | 20160 | 24 | 4680 5673 5160 | 0 | 720 | 4369 17474 8740 39048 | 6 |
| 5 | 20160 | 24 | 5673 4804 4748 | 6 | 720 | 4369 16962 14372 56968 | 12 |
| 6 | 13440 | 36 | 5673 5160 4681 | 3 | 480 | 4369 17474 25124 35208 | 4 |
| 7 | 10080 | 48 | 4680 4804 5770 | 9 | 360 | 4369 17474 29220 35224 | 6 |
| 8 | 2016 | 240 | 4680 4804 4748 | 8 | 72 | 4369 17474 25124 35224 | 6 |
| 9 | 1344 | 360 | 4680 4680 4804 | 7 | 48 | 4369 17474 8740 39320 | 4 |

The same ten representatives written out bit by bit (`expected/representatives.txt`, regenerated by run.sh as
`out/canon/representatives.txt`). A vector is written as `W_0|W_1|W_2|W_3`, each nibble as x1 x2 z1 z2. A matrix ψ is written
row by row, rows separated by `/`.

```
impl1 class 0 (orbit 80640 tuples = 2880 x 28; Sp(4,2)^4 x| S_4 stabiliser order 6)  <->  impl2 rep 5 (orbit 2880 triples)
  impl1 psi_1, psi_2, psi_3 = 1001/0100/0110/1000  0010/0011/0100/1000  0001/0100/0010/1000  (packed 5673 4804 5160)
  impl1 V_0 basis (w = e_0..e_3):  1000|1001|0001|0001  0100|0110|0010|0100  0010|0010|1100|0010  0001|1000|0100|1000
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0100  0010|0100|0100|0001  0001|0001|1001|0110  (packed 4369 9282 33316 27016)
impl1 class 1 (orbit 60480 tuples = 2160 x 28; Sp(4,2)^4 x| S_4 stabiliser order 8)  <->  impl2 rep 2 (orbit 2160 triples)
  impl1 psi_1, psi_2, psi_3 = 1001/0100/0110/1000  0011/0010/0100/1000  0010/0001/0100/1000  (packed 5673 4684 4740)
  impl1 V_0 basis (w = e_0..e_3):  1000|1001|0001|0001  0100|0110|0010|0010  0010|0010|1100|1000  0001|1000|1000|0100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|1100|0100  0001|0001|0001|1001  (packed 4369 17474 8996 39048)
impl1 class 2 (orbit 40320 tuples = 1440 x 28; Sp(4,2)^4 x| S_4 stabiliser order 12)  <->  impl2 rep 4 (orbit 1440 triples)
  impl1 psi_1, psi_2, psi_3 = 1001/0100/0110/1000  0001/0100/0010/1000  1011/0010/0100/1000  (packed 5673 5160 4685)
  impl1 V_0 basis (w = e_0..e_3):  1000|1001|0001|1001  0100|0110|0100|0010  0010|0010|0010|1100  0001|1000|1000|1000
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|1110  0001|0001|1001|0001  (packed 4369 17474 29220 35208)
impl1 class 3 (orbit 20160 tuples = 720 x 28; Sp(4,2)^4 x| S_4 stabiliser order 24)  <->  impl2 rep 1 (orbit 720 triples)
  impl1 psi_1, psi_2, psi_3 = 0001/0010/0100/1000  1001/0100/0110/1000  0010/0011/0100/1000  (packed 4680 5673 4804)
  impl1 V_0 basis (w = e_0..e_3):  1000|0001|1001|0001  0100|0010|0110|0010  0010|0100|0010|1100  0001|1000|1000|0100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0100  0010|0100|0100|0001  0001|0001|0001|0111  (packed 4369 9282 33316 59528)
impl1 class 4 (orbit 20160 tuples = 720 x 28; Sp(4,2)^4 x| S_4 stabiliser order 24)  <->  impl2 rep 0 (orbit 720 triples)
  impl1 psi_1, psi_2, psi_3 = 0001/0010/0100/1000  1001/0100/0110/1000  0001/0100/0010/1000  (packed 4680 5673 5160)
  impl1 V_0 basis (w = e_0..e_3):  1000|0001|1001|0001  0100|0010|0110|0100  0010|0100|0010|0010  0001|1000|1000|1000
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|0100  0001|0001|0001|1001  (packed 4369 17474 8740 39048)
impl1 class 5 (orbit 20160 tuples = 720 x 28; Sp(4,2)^4 x| S_4 stabiliser order 24)  <->  impl2 rep 6 (orbit 720 triples)
  impl1 psi_1, psi_2, psi_3 = 1001/0100/0110/1000  0010/0011/0100/1000  0011/0001/0100/1000  (packed 5673 4804 4748)
  impl1 V_0 basis (w = e_0..e_3):  1000|1001|0001|0001  0100|0110|0010|0010  0010|0010|1100|1000  0001|1000|0100|1100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0100|0010  0010|0100|0001|1100  0001|0001|0111|1011  (packed 4369 16962 14372 56968)
impl1 class 6 (orbit 13440 tuples = 480 x 28; Sp(4,2)^4 x| S_4 stabiliser order 36)  <->  impl2 rep 3 (orbit 480 triples)
  impl1 psi_1, psi_2, psi_3 = 1001/0100/0110/1000  0001/0100/0010/1000  1001/0010/0100/1000  (packed 5673 5160 4681)
  impl1 V_0 basis (w = e_0..e_3):  1000|1001|0001|1001  0100|0110|0100|0010  0010|0010|0010|0100  0001|1000|1000|1000
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|0110  0001|0001|1001|0001  (packed 4369 17474 25124 35208)
impl1 class 7 (orbit 10080 tuples = 360 x 28; Sp(4,2)^4 x| S_4 stabiliser order 48)  <->  impl2 rep 9 (orbit 360 triples)
  impl1 psi_1, psi_2, psi_3 = 0001/0010/0100/1000  0010/0011/0100/1000  0101/0001/0110/1000  (packed 4680 4804 5770)
  impl1 V_0 basis (w = e_0..e_3):  1000|0001|0001|0001  0100|0010|0010|1010  0010|0100|1100|0010  0001|1000|0100|1100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|1110  0001|1001|1001|0001  (packed 4369 17474 29220 35224)
impl1 class 8 (orbit 2016 tuples = 72 x 28; Sp(4,2)^4 x| S_4 stabiliser order 240)  <->  impl2 rep 8 (orbit 72 triples)
  impl1 psi_1, psi_2, psi_3 = 0001/0010/0100/1000  0010/0011/0100/1000  0011/0001/0100/1000  (packed 4680 4804 4748)
  impl1 V_0 basis (w = e_0..e_3):  1000|0001|0001|0001  0100|0010|0010|0010  0010|0100|1100|1000  0001|1000|0100|1100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|0110  0001|1001|1001|0001  (packed 4369 17474 25124 35224)
impl1 class 9 (orbit 1344 tuples = 48 x 28; Sp(4,2)^4 x| S_4 stabiliser order 360)  <->  impl2 rep 7 (orbit 48 triples)
  impl1 psi_1, psi_2, psi_3 = 0001/0010/0100/1000  0001/0010/0100/1000  0010/0011/0100/1000  (packed 4680 4680 4804)
  impl1 V_0 basis (w = e_0..e_3):  1000|0001|0001|0001  0100|0010|0010|0010  0010|0100|0100|1100  0001|1000|1000|0100
  impl2 V_0 basis:                 1000|1000|1000|1000  0100|0010|0010|0010  0010|0100|0100|0100  0001|1001|1001|1001  (packed 4369 17474 8740 39320)
```

The 72 Sp(4,2)-only representatives of impl2 are in `original_logs/impl2_auditor/v0reps8_sp.txt`, and each one's class is
listed in `expected/crosswalk.txt`. Searching all 72 removes any reliance on the party-permutation part of the gauge (S_4 on
the B-parties and the base-party swap). impl2 does this, and it reproduces the class counts exactly.

## 5. Expected results (all reproduced by run.sh, both implementations)

Full search at minimum party weight 5 (AME), per class. The columns are:
* #V_1: valid first completions. impl1 `nodes[1]`, impl2 `n_V1`; they are equal.
* level 2: triples (V_0,V_1,V_2) satisfying every condition except the four with J = {0,1,2}. impl1 frames3 `JMASKS=124`
  `nodes[2]`, impl2 `n_V2`.
* #V_2 (impl1, t = 5): impl1 level 2 with all conditions.
* solutions: impl1 `nsol`, impl2 `n_SOL`.

| impl1 class | impl2 rep | #V_1 | level 2 | #V_2 (impl1, t = 5) | solutions impl1 / impl2 |
|---|---|---|---|---|---|
| 0 | 5 | 26556 | 702 | 0 | 0 / 0 |
| 1 | 2 | 24128 | 176 | 0 | 0 / 0 |
| 2 | 4 | 24128 | 416 | 0 | 0 / 0 |
| 3 | 1 | 28960 | 1776 | 0 | 0 / 0 |
| 4 | 0 | 24488 | 712 | 0 | 0 / 0 |
| 5 | 6 | 28960 | 3640 | 0 | 0 / 0 |
| 6 | 3 | 24488 | 2444 | 0 | 0 / 0 |
| 7 | 9 | 28960 | 1712 | 0 | 0 / 0 |
| 8 | 8 | 28656 | 8912 | 0 | 0 / 0 |
| 9 | 7 | 28656 | 3744 | 0 | 0 / 0 |
| total | | 267,980 | 24,234 | 0 | 0 / 0 |

More detail on these numbers:
* At level 1, impl1's affine solution space has dimension 32. At level 2 it has dimension 16 for every one of the 267,980
  V_1, so no degenerate linear system is skipped (histograms `L1{32:1} L2{16:#V_1}` in the check `i1_search_t5_frames3`).
* In impl2, for all 24,234 level-2 triples, the forced V_3 passes its pair conditions and fails its single-element conditions
  (check `i2_v3diag_full`).
* The 72 Sp-only representatives give 72 × 0 solutions, with per-representative counts equal to those of their class
  (`i2_search_sp72`, `crosswalk`).
* Gauge invariance: 3 random Sp(4,2)^4 × S_4 transforms of each impl2 representative give identical counts (`i2_gauge`).

**Positive controls.** The same code finds solutions wherever they exist:

| check | what | expected |
|---|---|---|
| `i1_ctrl_n4` | impl1, n = 4 (m = 2, t = 3) | 1 class (336 tuples), 1 completion. Number of labelled stabilizer AME(4,4) = 8640 · 720² = **4,478,976,000** |
| `i1_ctrl_n6` | impl1, n = 6 (m = 3, t = 4) | 5 classes with **112, 64, 288, 368, 2272** completions. Σ (orbit/28)·nsol = 54400, so N(6) = 54400 · 720^5 = **10,525,951,918,080,000,000** |
| `i2_ctrl_counts_n4_n6` | impl2 `vcount.py`, separate code | N(4) = 4,478,976,000. n = 6: 11 Sp-orbits with t = [2272, 288, 368, 288, 288, 64, 112, 368, 112, 368, 112], Σ = 54400, the same N(6) |
| `cross_n4_n6` | both | equal N(4), N(6) and the same set of n = 6 completion counts |
| `i1_ctrl_m4t4` | impl1, n = 8 with minimum weight 4 instead of 5, class 9 | frames found at once (stops after 20: nodes 1, 2, 1700, 20) |
| `i1_ctrl_relaxed` | impl1 frames3 `JMASKS=64479` (drops the six (2,2) conditions with J = {0,2} and their duals {1,3}), class 9 | 60 leaves (nodes 1, 24, 60, 60) |
| `i1_framecheck_*` | brute force (`tools/framecheck.py`) on every frame found above | all frames orthogonal and nondegenerate. n = 6: all 3104 have party-weight enumerator [1,0,0,0,225,1080,2790], i.e. are AME(6,4). n = 8 controls: minimum weight exactly 4. Relaxed leaves: impl1's (J,L) condition evaluator agrees with an exhaustive word check on all 68 conditions for all 60 |
| `i2_ctrl_relaxed` | impl2 mode 2 (thresholds of minimum weight 4), rep 0, first 400 V_1 | n_V2 = 16,619,964, n_SOL = 8,215,834 |
| `i2_vbrute_relaxed` | brute-force minimum party weight of the first 50 of those | {4: 50} |
| `i2_nearmiss_rep0`, `i2_vbrute_nearmiss` | the first 30 near-misses of rep 0 in mode 0 | minimum weight 4. A_4 = 30 (18 cases), 6 (6), 12 (2), 18 (4) |
| `i2_enumtest` | exactness of impl2's subspace enumerator | 651, 11811, 200787 = Gaussian binomials [6,4]_2, [7,4]_2, [8,4]_2. Random candidate sets: DFS count = brute force / 20160 |
| `i2_classcount` | orbit sizes recounted with all of GL(4,2) × S_h | h = 2: [336]. h = 3: [5040, 2520, 1260, 1008, 28]. h = 4: the ten sizes above |

**Count hashes.** `expected/HASHES.txt` lists the SHA-256 of each canonical result file. run.sh recomputes the hash of every
replayed file and compares. The hashes of the main results are:

```
729bc3e16c1160f07eb1651c0d1d6c874c7c540475499ac3c385e82d255d5562  i1_search_t5_frames3.txt   (impl1 full search, 10 classes)
307598d8e1904038529329b2a7e569e43ad250549da9fab9b40660bbb2a2a66a  i1_search_t5_frames.txt    (impl1, first version frames.c)
8cdfb6be037de6cf8c08b06a21d27f49541098dadffc67a0d63f2d978d70efa6  i1_level2_m124.txt         (impl1 level-2 counts, total 24234)
cec37ce119168472cfc261fd1c74a6973d27d3405d9e6ca44277d62b4e1273c7  i2_search_full.txt         (impl2 full search, 10 classes)
d072a3ca114af954c9a7bec47bbf13d4dc80aae635d970a4fc7cdf93440bdfec  i2_search_sp72.txt         (impl2, 72 Sp-only classes)
5be137babec28d759a30ccd4cf3c03551f93295843017b11f5d16edda279fa89  crosswalk.txt              (class correspondence + agreement)
```

## 6. What run.sh does

1. It compiles `frames.c`, `frames3.c` (impl1), `vsearch8.c`, `vcount.c` and `enumtest.c` (impl2) into `build/`.
2. It re-derives the V_0 representatives with each implementation's own code: `v0reps.py 2|3|4` (impl1) and `vclasses8.py`
   (impl2: the 10 full-gauge and 72 Sp-only representatives). The searches then use these freshly derived files, not the
   recorded ones.
3. It runs the fast controls: impl1 n = 4, n = 6, n = 8 t = 4 and the relaxed n = 8 run; impl2 `vcount.py`.
4. It runs 125 search tasks in parallel under `nice -n 10`:
   * impl1 `frames3` and `frames` at t = 5 for the 10 classes;
   * impl1 `frames3 JMASKS=124` for the 10 classes;
   * impl2 `vsearch8` for the 10 + 72 representatives;
   * the impl2 gauge test for the 10 representatives;
   * impl2 relaxed control, `classcount.py 4` and `enumtest`.
5. It brute-forces the impl2 relaxed solutions and near-misses (`vbrute8.py`).
6. It canonicalises every result (`tools/checks.py replay` -> `out/canon/`, using `tools/canon.py`) and compares each file
   byte for byte with `expected/` and with `expected/HASHES.txt`. It also regenerates `expected/` from `original_logs/` into
   `out/expected_regen/` and checks that the two are identical. This shows that the expected values are mechanically derived
   from the recorded logs.

`python3 tools/checks.py list` prints every check with the source of its expected value. Most expected files are
canonicalised copies of the original logs (`original_logs/`). A few runs were never logged to a file. For those, the
expected values are the numbers documented in PROOF.md / REVIEW.md, written as constants in `tools/checks.py` together with
their source:
* `i1_classes_summary`
* `i2_classes_summary`
* `i2_ctrl_counts_n4_n6`
* `i2_classcount`
* `i2_enumtest`
* `i2_vbrute_relaxed`
* `i2_vbrute_nearmiss`

## 7. Runtimes

Measured on a fresh copy of this folder outside the repository on
2026-10-31. Hardware: Apple M4 (10 cores), macOS, Apple clang 21, Python 3.9.6. The machine was shared with other jobs (load
average 73 at the start and 34 at the end), so wall times are pessimistic. The C programs report CPU time (`time=`), which is
less sensitive to the load.

| step | this replay | original logs (2026-10-16) |
|---|---|---|
| impl1 frames3 t = 5, one class (CPU) | 28.0-30.1 s | 25.1-28.8 s |
| impl1 frames.c t = 5, one class (CPU) | 27.2-29.2 s | 21.0-25.6 s |
| impl1 frames3 JMASKS=124, one class (CPU) | 28.3-30.5 s | 20.8-22.3 s |
| impl2 vsearch8, one representative (wall, loaded) | 7.4-9.0 s | 3.1-3.9 s |
| impl2 relaxed control, 400 V_1 (wall, loaded) | 130 s | 34 s |
| class derivation v0reps.py 2,3,4 / vclasses8.py (wall) | 15 s / 11 s | -- |
| **whole run.sh, JOBS=4 (wall)** | **757 s = 12.6 min, OVERALL: PASS (38/38 checks)**; user CPU 1447 s, system 20 s | -- |
| first test run (scratch copy, same load) | 851 s; 37/38 (the line-count issue of Sec. 9, since fixed) | -- |

## 8. Files

```
README.md            this file
PROOF.md             the completeness argument (runs/classify/PROOF.md of the project, unchanged below a release note)
REVIEW.md            the auditor's review (results/argument.md 5.12-5.14, 5.18, 5.19; closure note), verbatim excerpts
run.sh               the replay
impl1_classify/      frames.c, frames3.c (search; frames3 = frames + diagnostics and the JMASKS/PATTERN/DROPJ control hooks),
                     v0reps.py + f2m.py (V_0 classes), run.py (nice + log wrapper), run_m124_all.py (JMASKS=124 driver),
                     RESULT_TABLE.md (the original result table)
impl2_auditor/       vclasses8.py (V_0 classes), vsearch8.c (search), run_v8.py (driver), gauge_test8.py, run_v8c.py (relaxed
                     control), vbrute8.py (brute-force weights), vcount.py + vcount.c (n = 4, 6 counts), classcount.py
                     (orbit-size recount), enumtest.c (enumerator exactness)
tools/               written for this release: checks.py (list of checks), canon.py (normalisers), crosswalk.py (class
                     correspondence by invariant), framecheck.py (brute-force frame checker), hashes.py
expected/            canonical expected results, one file per check, and HASHES.txt
original_logs/       the recorded logs and output files from the project (impl1_classify/, impl2_auditor/); repro.log,
                     repro_all_s8.log and c_m124.log are earlier reruns with identical counts, kept for provenance only
```

### Changes relative to the project's files (provenance)

Byte-identical copies:
* impl1: `frames.c`, `frames3.c`, `v0reps.py`, `f2m.py`, `RESULT_TABLE.md`.
* impl2: `vsearch8.c`, `vcount.c`, `enumtest.c`, `classcount.py`.

Adapted files. Each starts with a `RELEASE NOTE` describing the change.
* `impl1_classify/run.py`: the log header now also records `JMASKS`/`LOGL2`. The original recorded only DROPJ/PATTERN, so its logs show
  `env={}` for JMASKS runs.
* `impl1_classify/run_m124_all.py`: binary, class file and log are arguments, plus an optional class list.
* `impl2_auditor/vclasses8.py` and `impl2_auditor/vbrute8.py`: ported from numpy to the standard library.
  * `vclasses8.py`: same algorithm; output files byte-identical to the originals. Checked by run.sh and, during preparation,
    against the original numpy script, whose stdout is identical too.
  * `vbrute8.py`: same formula, same results on the recorded near-misses. One extra output line gives the (minimum weight,
    multiplicity) histogram.
* `impl2_auditor/run_v8.py`, `gauge_test8.py`, `run_v8c.py`, `vcount.py`: binaries and files are passed as arguments instead of being
  compiled or read next to the script. `run_v8.py` writes only the RESULT line into its log and puts the diagnostic 'SOL'
  lines in a side file. `run_v8c.py` streams the 8.2 million solution lines instead of holding them in memory (about 0.6 GB in
  the original).
* `original_logs/`: copies of the project's files. In 6 logs (`final_t5.log`, `m124_all.log`, `ame84_A.log`, `ame84_B.log`,
  `repro.log`, `repro_all_s8.log`) the absolute path prefix of the author's machine was removed, so the paths read relative to
  `problems/p10-ame-states/`. No other edit was made; all other files are byte-identical.

## 9. What is not reproduced exactly, and other caveats

* **One recorded number cannot be reproduced: the raw count of 'SOL' output lines in impl2's relaxed control.**
  * Recorded (`run_v8c.log`): `#SOL lines 8215834`, equal to n_SOL = 8,215,834.
  * Replay: `#SOL lines 8215864`, i.e. 30 more, while n_SOL = 8,215,834 and n_V1, n_V1_with_V2, n_V2 are identical.
  * Cause: the current `vsearch8.c` (see the next item) prints, in every mode, the first 30 near-misses with the same `SOL`
    prefix (`if (n_pair_only <= 30) printf("SOL ...")` in `cb_v2`). In mode 2 these are relaxed near-misses. The binary of the
    original run was compiled before that edit and printed only the solutions.
  * The check `i2_ctrl_relaxed` therefore compares the search counts (n_V1, n_V1_with_V2, n_V2, n_SOL) and not the raw line
    count. The difference was observed in the replay and is not hidden: the raw count is in `out/logs/i2_relaxed.log`.
  * All 50 solution lines that `vbrute8.py` checks have minimum party weight exactly 4. A near-miss line would have weight
    <= 3, so none of the 50 is a near-miss.
* **vsearch8.c was edited after the main impl2 runs.** The logs `run_v8_0_0/1_4/5_9.log` (22:47-22:48 on 2026-10-16) predate
  the last modification of `vsearch8.c` (22:53). That edit added the V_3 diagnostic counters and the printing of the first 30
  near-misses as 'SOL' lines; the earlier version was not kept. The release uses the current file. All recorded counts
  (n_V1, n_V1_with_V2, n_V2, n_SOL, cand1) are reproduced exactly with it, and so are the diagnostic counters, which match
  `run_v8h.log` (written by `vsearch8h.c` = vsearch8.c + histograms).
* **The 50 relaxed solutions were not kept.** `run_v8c.py` wrote them to `work/auditor/v8c_solutions.txt`, which
  `run_v8n.py` later overwrote with the 30 near-misses of rep 0. For the relaxed control the replay therefore compares the
  recorded counts (`run_v8c.log`) and the documented brute-force result ({4: 50}), not the solution lines. The 30 near-misses
  are compared line by line.
* **The gauge test seed was not recorded.** The replay uses seeds 1000 + rep, so its random transforms differ from the
  original ones. Only the counts and the agreement (30/30 identical) are compared.
* **enumtest uses the C library's rand().** On a different platform the random candidate sets differ. Only the Gaussian
  binomial counts and the equality DFS = brute force / 20160 are compared.
* Some expected values are documented numbers rather than logs (listed in Sec. 6).
* **Not part of this release, because they are not needed for exhaustiveness:**
  * classify's numpy recounts (`recount_level1.py`, `recount_level2.py`, `classcheck.py`, `pattern_check.py` on the
    `frames2.c` control);
  * the n = 6 Monte Carlo (`mc6.c`: 574 hits in 4·10^7 vs 584 ± 24 predicted);
  * the SAT encodings (inconclusive at n = 8, see PROOF.md §10);
  * the near-miss census programs (`vsearch8h.c`, `vsearch8h2.c`).

  They are described in PROOF.md and REVIEW.md and are available in the project repository. `tools/framecheck.py` (new)
  replaces the brute-force verifications of found frames for which classify used numpy.
