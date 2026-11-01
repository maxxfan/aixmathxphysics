> **Release note (release/frames, 2026-10-31).** This is `runs/classify/PROOF.md` of the project repository, copied without
> edits below this box.  Paths in it are relative to `problems/p10-ame-states/`; in this release the code of `work/classify/`
> is in `impl1_classify/` (frames.c, frames3.c, v0reps.py, f2m.py, run.py, run_m124_all.py), the auditor's code in
> `impl2_auditor/`, and the logs cited (final_t5.log, ame84_*.log, m124_all.log, cal.log, relax_leaf.log, repro_all_s8.log,
> run_v8_*.log, gauge_test8.log, ...) in `original_logs/`.  Section 8 items 2-6c that use numpy tools (frame2stab.py,
> pattern_check.py, recount_level1.py, recount_level2.py, classcheck.py) and the Monte Carlo (mc6.c) and SAT runs (section 10)
> are NOT replayed by run.sh; their role is covered in the replay by tools/framecheck.py (brute-force check of every frame
> found by the controls, including the 68-condition pattern) and by the second implementation.  See README.md.

# Exhaustive classification of F_2-linear stabilizer AME(8,4) states (classify, 2026-10-16)

**Claim (computer-assisted). Status 2026-10-16 23:00: AUDITED — the auditor re-derived every section (results/argument.md 5.13, board #327) and
reproduced the whole search with an independent implementation (work/auditor/vsearch8.c: different enumeration order, own V_0 classification;
board #322-#324, #329): same 10 classes, identical valid-V_1 counts in all classes, 0 survivors, same n=4 and n=6 counts.** There is **no** stabilizer state on 8 ququarts (ququart = 2 qubits; S ⊂ F_2^32 isotropic for the
pairwise symplectic form, dim S = 16) in which every nonzero element has support on at least 5 parties. Equivalently there is
no F_2-linear (qubit-pair) stabilizer AME(8,4). Nothing is assumed about symmetry. Z_4 stabilizers are excluded separately
(pairs #191(2), audited #199); non-stabilizer states and non-group error bases are NOT covered.

Everything below is elementary linear algebra over F_2; the only computation is the enumeration in Step 7.

## 0. Notation
Parties 0..7. A = {0,1,2,3}, B = {4,5,6,7}. Each party has the space P = F_2^4 with coordinates (x1,x2,z1,z2) and the
alternating form J(a,b) = x·z' + z·x' (so J pairs x1↔z1, x2↔z2). The Pauli space is (P^8, ω = ⊕J). In code J is the 4×4
matrix with ones at (0,2),(1,3),(2,0),(3,1); 4×4 matrices over F_2 are packed in 16 bits, entry (r,c) = bit 4r+c.
Sp(4,2) = isometry group of (P,J), order 720; GL(4,2) has order 20160. Over F_2, J^{-1} = J, and a 4×4 alternating matrix is
nondegenerate iff it is invertible (28 of the 64 alternating matrices).

## 1. Normal form (information set)
Let S be a stabilizer AME(8,4), i.e. S is a 16-dim isotropic subspace of P^8 (maximal, since ω is nondegenerate) with no nonzero
element of party-support ≤ 4.
* The restriction map S → P^A (coordinates of parties 0..3) has kernel = elements supported on B (4 parties), which is 0.
  dim S = 16 = dim P^A, so it is a bijection; the same for B. Hence S = {(u, Mu) : u ∈ P^A} with M: P^A → P^B a linear bijection.
* Isotropy: ω_A(u,u') + ω_B(Mu,Mu') = 0, and in characteristic 2 this says M is an **isometry** (P^A, ω_A) → (P^B, ω_B).
* Write W_l (l = 0..3) for the party-(4+l) copy of P inside P^B, π_l for the coordinate projection onto W_l, and
  A_j (j = 0..3) for party j's copy of P inside P^A. Put **V_j := M(A_j)**, 4-dim subspaces of P^B. Because the A_j are
  mutually orthogonal and nondegenerate and M is an isometry and a bijection:
  **V_0,…,V_3 are mutually ω-orthogonal, nondegenerate, and P^B = V_0 ⊕ V_1 ⊕ V_2 ⊕ V_3.**
* Conversely, given such a decomposition and any isometries E_j: (P,J) → (V_j, ω|V_j) (they exist because each V_j is a
  nondegenerate 4-dim alternating space; they are unique up to Sp(4,2)), M|A_j := E_j defines an isotropic S of dim 16.
  Different choices of E_j differ by the local Clifford group Sp(4,2) of party j. **So S up to local Cliffords on the A-parties is
  exactly the unordered/ordered frame (V_0,…,V_3), and the E_j never matter.**

## 2. The AME condition in frame language
For J ⊂ {0..3} put V_J = ⊕_{j∈J} V_j, and for L ⊂ {0..3} put W_L = ⊕_{l∈L} W_l.
A nonzero v ∈ P^B corresponds to the stabilizer element s_v = (M^{-1}v, v). Its A-support is {j : v_j ≠ 0} (v = Σ v_j with
v_j ∈ V_j), its B-support is {l : π_l v ≠ 0}. So  v ∈ V_J ∩ W_L ⇔ A-support ⊂ J and B-support ⊂ L.
**Lemma 2.** S is AME(8,4) ⇔ (C): V_J ∩ W_L = 0 for all J, L ⊂ {0..3} with |J| + |L| = 4.
Proof. S has a nonzero element of weight ≤ 4 iff some nonzero v has |A-supp| + |B-supp| ≤ 4; enlarge J ⊇ A-supp and
L ⊇ B-supp to sizes summing to exactly 4 (possible since each can grow up to size 4). Conversely a nonzero v in V_J ∩ W_L
with |J|+|L| = 4 has weight ≤ 4. ∎
The conditions with |J| or |L| = 0 are automatic (V_J ∩ 0, resp. v = 0 for J = ∅). Remaining pairs: (|J|,|L|) = (1,3): 16;
(2,2): 36; (3,1): 16 — 68 conditions. (Under (J,L) ↦ (J^c,L^c), V_J ∩ W_L = 0 ⇔ V_{J^c} ∩ W_{L^c} = 0 because
(V_J ∩ W_L)^⊥ = V_{J^c} + W_{L^c}; so there are 34 independent conditions. The code does not use this.)
Since dim V_J = 4|J|, V_J ∩ W_L = 0 ⇔ the projection of V_J onto W_{L^c} (the coordinates of the other parties) is injective.

## 3. Gauge
Let G = (Sp(4,2) ≀ S_4) acting on P^B (local Cliffords of the B-parties and permutations of the B-parties) and S_4 acting by
renumbering the V_j (permuting A-parties). G preserves ω, orthogonality, nondegeneracy and the family of conditions (C)
(permutations permute the pairs (J,L)). A frame is a solution iff its G-image is, so **it suffices to find all solutions with V_0
a prescribed representative of its G-class, and the remaining V_1,V_2,V_3 arbitrary (ordered).** (Every solution frame has a
member V_j; renumber it to V_0.) We therefore need all G-classes of V_0 that can occur, and then ALL completions.

## 4. Graph form
(C) with (|J|,|L|) = (1,3) says V_j ∩ W_{L} = 0 for every 3-set L, i.e. π_l|V_j is injective, hence bijective, for every l.
In particular every V_j is a graph over W_0: **V_j = {(w, θ_{j,1}w, θ_{j,2}w, θ_{j,3}w) : w ∈ W_0}** with θ_{j,l} ∈ GL(4,2)
(θ_{j,0} := I). The ω-form of V_j in the basis w is Q_j = J + Σ_{l=1}^{3} θ_{j,l}^T J θ_{j,l} (using ω((w,…),(w',…)) =
w^T J w' + Σ_l (θ_l w)^T J (θ_l w')); V_j nondegenerate ⇔ Q_j invertible.
Orthogonality V_i ⊥ V_j (i ≠ j) ⇔ **Σ_{l=0}^{3} θ_{i,l}^T J θ_{j,l} = 0** (4×4 matrix identity; l = 0 term is J) — this is **linear in
θ_j = (θ_{j,1},θ_{j,2},θ_{j,3}) ∈ F_2^48** for fixed θ_i.

## 5. Classification of V_0 up to G
Let V_0 satisfy the (1,3) conditions. Each π_l|V_0: V_0 → W_l is a bijection, so Q_l := (π_l|V_0)^* J is a nondegenerate alternating
form on the 4-dim space V_0 (28 choices in coordinates), and Σ_l Q_l = ω|V_0 must be nondegenerate. Claim: the G-class of V_0 is
the orbit of the 4-tuple (Q_0,…,Q_3) under GL(V_0) (congruence) × S_4 (permuting the tuple).
(⇒) a gauge element (D_l, σ) with D_l ∈ Sp(W_l → W_σ(l)) maps V_0 to V_0' via the linear iso h = π'_{σ(l)}^{-1} D_l π_l|V_0 (independent of l),
and h^*Q'_{σ(l)} = Q_l. (⇐) if h: V_0 → V_0' and σ satisfy h^*Q'_{σ(l)} = Q_l then D_l := π'_{σ(l)} h (π_l|V_0)^{-1}: W_l → W_{σ(l)} is an
isometry and (D_l) maps V_0 onto V_0'. ∎
`work/classify/v0reps.py 4` enumerates all 28^4 tuples with nondegenerate sum (268,800 of them), generates the orbits under
the 12 elementary transvections generating GL(4,2) and the adjacent transpositions of S_4 by flood fill, asserts that the orbit
sizes add up to 268,800, normalises each orbit representative to Q_0 = J by a g ∈ GL(4,2) (coordinates w ∈ W_0) and chooses ψ_l with
ψ_l^T J ψ_l = Q'_l. **Result: 10 classes** (V_0 = {(w,ψ_1 w,ψ_2 w,ψ_3 w)}), orbit sizes (in tuples) 80640, 60480, 40320,
20160, 20160, 20160, 13440, 10080, 2016, 1344; file `work/classify/v0_m4.txt` (ψ values in packed form):
| class | orbit | ψ_1 ψ_2 ψ_3 |
|---|---|---|
| 0 | 80640 | 5673 4804 5160 |
| 1 | 60480 | 5673 4684 4740 |
| 2 | 40320 | 5673 5160 4685 |
| 3 | 20160 | 4680 5673 4804 |
| 4 | 20160 | 4680 5673 5160 |
| 5 | 20160 | 5673 4804 4748 |
| 6 | 13440 | 5673 5160 4681 |
| 7 | 10080 | 4680 4804 5770 |
| 8 | 2016 | 4680 4804 4748 |
| 9 | 1344 | 4680 4680 4804 |
(Same procedure: n = 4 → 1 class of 336 tuples; n = 6 → 5 classes of 5040, 2520, 1260, 1008, 28 tuples, 9856 total.)

## 6. Enumeration (what `work/classify/frames.c` does; same logic in frames2.c/frames3.c plus diagnostics)
Fix a class representative V_0 (data ψ_l). For j = 1, 2, 3 in turn, given V_0..V_{j-1}:
1. Solve the linear system Σ_l θ_{i,l}^T J θ_{j,l} = J-term for all i < j (16j equations in 48 unknowns) by Gaussian elimination:
   an affine space of dimension f_j (observed: f_1 = 32, f_2 = 16, f_3 = 0 for every node).
2. Enumerate ALL 2^{f_j} points in Gray-code order. Test: θ_{j,l} ∈ GL(4,2) for l = 1,2,3; Q_j invertible; and **every condition of (C)
   whose J contains j and is contained in {0..j}** (rank test of the stacked basis vectors masked to W_{L^c}).
3. Recurse on every point that passes; at j = 3 the system is square and has at most one solution.
**Nothing is lost**: each test in 2 is a necessary condition for a frame to be a solution (§§2, 4), applied to a frame prefix, and the
union over j of the conditions tested is all of (C); a leaf that survives all four levels satisfies (C), orthogonality and nondegeneracy,
hence is an AME frame (and is then converted to stabilizer rows and checked by rankcheck).
The search is exhaustive over: all (V_1,V_2,V_3) for each V_0 class (hence, by §§3,5, over all AME frames).

## 7. Result
Per class (all 10 classes, `work/classify/final_t5.log`; the first run `ame84_*.log` with frames.c gave identical node counts):
see the table in `runs/classify/RESULT_TABLE.md`. **Level-1 survivors V_1: 24128–28960 per class; level-2 survivors (V_0,V_1,V_2): 0 for
every class; hence 0 AME frames.** Time ≈ 25 s per class.
Reproduction pass (shift 8, 2026-10-16 23:52-23:56): `frames3.c` recompiled to `frames3_s8` and rerun on all ten classes (`work/classify/repro_all_s8.log`);
`cmp_repro.py` shows all ten (orbit size, node counts, nsol) tuples identical to `final_t5.log`; total valid V_1 = 267,980, total V_2 = 0, nsol = 0.

## 8. Controls and cross-checks (all done by the same author; the auditor's independent implementation is the next step)
1. n = 4: 1 class (336 tuples). V_0 count = 336/28·720 = 8640; every V_0 has V_1 = V_0^⊥ valid ⇒ 8640·720² = 4,478,976,000 AME(4,4) F_2-stabilizer
   states, equal to the pairs first-moment note's "exact 4.478e9" (work/pairs/moment_pairs.py).
2. n = 6, t = 4 (AME(6,4)): all 5 classes have completions (112, 64, 288, 368, 2272); frames converted to stabilizers pass rankcheck
   (bad_cuts 0/10) and have exact party-weight histogram [1,0,0,0,225,1080,2790] (all 4096 words). Total predicted number of AME(6,4)
   stabilizers = 54400·720^5 = 1.0526e19 = 1.461e-5 of all 7.2e23 Lagrangians; Monte Carlo with uniformly random Lagrangians
   (`mc6.c`, sampler asserted isotropic, seed 3): **40,000,000 samples, 574 hits, fraction 1.435e-5 vs predicted 1.461e-5 (expected 584 ± 24)**;
   a 1,000,000-sample pilot (seed 2) gave 15 hits (expected 14.6). This validates, end to end at m = 3, the frame theory (Lemma 2, graph form,
   orthogonality equations), the V_0 class reduction with its orbit weights, and the enumeration (the weights 54400 = Σ_c (orbit_c/28)·nsol_c enter the prediction).
   (A first version of mc6.c copied pivot rows too early and found 0 hits in 2·10^7 samples; fixed before the real run and isotropy is now asserted.)
3. n = 8, t = 4 (graph form forced): frames found at the same code path; all have exact min party weight 4 (2^16-word enumeration, frame2stab.py).
4. C condition evaluator vs exhaustive word check: for 40 n=8 frames, the validity of every one of the 68 conditions (J,L) agrees with
   "exists a nonzero stabilizer word with A-support ⊂ J, B-support ⊂ L" (pattern_check.py).
5. Independent recount of level 1 in numpy with a different algorithm (φ_1, φ_2 ∈ GL(4,2) scanned, φ_3 solved, Schur-complement conditions;
   recount_level1.py): class 9: 28656 = frames.c.
6. Independent recount of level 2 (V_2 = section + Hom(W_0, K), K = ker π_0|U with U = (V_0+V_1)^⊥, Python-int rank tests; recount_level2.py)
   against frames3 under relaxed condition sets so that completions exist: J = {0,2}: 135 = 135; J = {1,2}: 169 = 169; J = {0,1,2}: 1546 = 1546
   (25 level-1 nodes each, class 9, per-node counts all equal). The full level-2 test is the conjunction of these three families.
6b. Gauge-invariance of the class reduction (classcheck.py): 4 random graph-form V_0 (ψ_l uniform in GL(4,2)) were assigned to classes 0, 4, 0, 4 by the
   flood-fill table, and the independent numpy recount of valid V_1 gave 26556, 24488, 26556, 24488 = the representatives' counts from frames.c.
6c. Leaf-level control at t = 5 (frames3 with JMASKS=64479: all conditions except the six (2,2) conditions with J = {0,2} and their duals J = {1,3}, class 9):
   60 leaves found at once; all 60 pass the isotropy assertion and have exact min party weight 4; the verifier's rankcheck gives bad_cuts 1-4 of 35
   (≤ 6 possible) — so levels 1-3 (including the forced level-3 solve) deliver genuine stabilizer states, and pattern_check agrees with the exhaustive
   word check on all 68 conditions for all 60.
7. Free-dimension histograms from the instrumented run: level 1: 32 (all), level 2: 16 (all V_1) — the generic dimensions, no degenerate systems.
8. Expectation (shift 1; **withdrawn as evidence after auditor review, 2026-10-17 00:20, argument.md 5.18**): random-model count of level-2 survivors per class
   ≈ 2^-19·28,000 ≈ 0.05, total ≈ 0.5. This heuristic is NOT used in the argument: the conditional counts do not follow the independence model. 24,234 triples
   (V_0,V_1,V_2) pass every condition except the four with J = {0,1,2} (my rerun with JMASKS=124 reproduces the auditor's total, per-class 702, 176, 416, 1776,
   712, 3640, 2444, 1712, 8912, 3744; `work/classify/m124_all.log`); independence for the four omitted conditions predicts 24234·0.3076^4 ≈ 220 survivors, observed 0, and
   V_3 always contains 3..15 nonzero vectors meeting at most three B-parties (auditor, argument.md 5.14(c)). So the null is a pattern we cannot explain; the theorem rests
   on the enumeration alone.

## 9. What would break the argument (checklist for the auditor)
* §1 (information set) and Lemma 2 — pure algebra.
* §5: orbit claim (⇒/⇐) and the flood-fill (generators of GL(4,2): elementary transvections generate; S_4: adjacent transpositions).
* §4: formula for orthogonality (convention ω(x,y) = x^T J y, matrices acting on column vectors) and the Gaussian elimination in frames.c.
* §6 step 2: that tests are necessary conditions only (no heuristic pruning anywhere).
The independent implementation should NOT use the V_0 class table or frames.c; re-derive the classes (e.g. directly as Sp(4,2)^4 ⋊ S_4 orbits
of V_0 on a different parametrisation) and enumerate V_1, V_2 by its own route.

## 10. SAT route (shift 2, 2026-10-16 23:05-23:40): calibration only, NO UNSAT certificate
`work/classify/sat_ame.py` encodes the normal form of §1 directly (no frames, no classification): M ∈ F_2^{16×16}; isotropy ω_B(Me_a, Me_b) = ω_A(e_a, e_b)
(120 XORs of 16 ANDs); every party-aligned k×k block M[R,C] (k = 1..3, or only k ≤ 2 with option `dual`, justified by the (J,L) ↔ (J^c,L^c) duality of §2) gets
an inverse variable matrix Y with M[R,C]·Y = I. Options: `v0=c` fixes the columns of A_0 to the class-c representative (§5, `v0cols.py`), `gaugeA`
fixes the local Cliffords of A-parties 1..3 completely (X_j = M[W_0,A_j] ranges over one representative of each of the 28 right cosets X·Sp(4,2),
invariant X J X^T). Solutions are converted to stabilizer rows and checked with exact weight enumeration + rankcheck (`check_sat.py`).
* Positive: n = 4 (0.0 s, bad_cuts 0/3), n = 6 (0.5 s unfixed; 5/5 V_0 classes with and without gaugeA, ≤ 0.8 s; `dual` with k = 1 only: 0.2 s, bad_cuts 0/10, hist
  [1,0,0,0,225,1080,2790] = the AME(6,4) enumerator), n = 8 relaxed (the six (2,2) blocks C ∈ {(0,2),(1,3)} dropped, V_0 = class 9): 27 s (120 s with gaugeA), solution
  has party-weight histogram [1,0,0,0,30,720,4800,20880,39105] and bad_cuts 1/35, identical to the relaxed leaves of frames3 (§8.6c).
* Negative / inconclusive: n = 8 full (68 blocks): unfixed CaDiCaL: no answer in 5 min; V_0 = class 9 fixed: CaDiCaL TIMEOUT 120 s; V_0 = class 9 and gaugeA: CaDiCaL 1.9.5 TIMEOUT 600 s, Kissat 4.0.4 TIMEOUT 600 s, CryptoMiniSat 5.17 (native XOR clauses, `dual`) TIMEOUT 270 s; Kissat with `dual` + gaugeA: TIMEOUT 400 s;
  **positive control failure: n = 10 (known to exist), `dual` encoding (125 blocks, 113k vars), unfixed CaDiCaL: TIMEOUT 330 s; CryptoMiniSat with native XORs: TIMEOUT 400 s**, so CDCL cannot find an AME(10,4) F_2-stabilizer
  either and the n = 8 timeouts carry no weight. XOR-heavy instances with a large residual symmetry group are the known weak spot of CDCL (the role brief anticipated this).
Conclusion: the SAT encoding reproduces the n ≤ 6 and relaxed-n=8 facts of §8 (a third, independent encoding of the normal form agrees with the frame
machinery on every satisfiable control), but it is not a proof of the n = 8 nonexistence. The nonexistence rests on §§1-7 and the two exhaustive C implementations.

## 11. Prior literature (added 2026-10-16 23:50 after auditor alert #358) -- this is a CONFIRMATION, not a new theorem
Ball, Moreno, Simoens, *Stabiliser codes over fields of even order*, arXiv:2401.06618 (v1 12 Jan 2024, v2 9 Sep 2024) state that there is no [[7,1,4]]_4
and no [[8,0,5]]_4 stabiliser code (Pauli group labelled by F_{2^h}, stabiliser F_2-linear = our model with h = 2), obtained by a different method (quantum sets in
symplectic polar spaces, a computer classification of [[6,2,3]]_4, Sage + GAP/FinInG code). Our §§1-7 (V-picture, 10 V_0 classes, 0 frames) and the auditor's
independent re-implementation are therefore an independent computational confirmation of their [[8,0,5]]_4 statement, not a first proof. We did NOT run their code and
did NOT check their proof; the sentence "no stabiliser AME(8,4) with parties F_2^4" is claimed by us only as far as §§1-8 and the audit establish it.
Remark on [[7,1,4]]_4 (not used anywhere above): a pure [[7,1,4]]_4 code is MDS (n-k = 2(d-1)), its stabiliser has minimum weight 5, and purifying with one
reference ququart gives a stabiliser state whose nonzero elements all have weight >= 5, i.e. an AME(8,4); conversely projecting one party of an AME(8,4)
stabiliser state onto |0,0> gives a [[7,1,4]]_4. For stabiliser codes the two nonexistence statements are therefore equivalent PROVIDED every [[7,1,4]]_4 stabiliser
code is pure (I recall that quantum MDS codes are pure, Rains 1999, but I did not re-check that; this remark is ours and not audited). What remains open for AME(8,4)
is non-stabiliser states only.

**Confirmation of the citation (auditor, board #380, 23:52; the auditor fetched and read arXiv:2401.06618v2, work/auditor/bms/, results/audits.json key
literature_bms_2401_06618_fulltext; I did not read the paper myself).** (1) The abstract states that there is no [[7,1,4]]_4 code and no [[8,0,5]]_4 code and that
the latter "rules out the existence of a stabilised absolutely maximum entangled state on 8 quqits with local dimension 4" (their [17] = Huber-Wyderka table); Sec. 5.2
derives [[8,0,5]]_4 from [[7,1,4]]_4 by Lemma 49 ([[n,k,d]]_{2^h} gives [[n-1,k+1,d-1]]_{2^h}). (2) Their model (Defs. 1/3, Thms 17/21) is the F_2-linear
(additive) stabiliser inside the trace-symplectic dual, i.e. exactly our two qubits per party, not only F_4-linear; Z/4Z-labelled Paulis are a different model. (3) Purity of
[[7,1,4]]_4: their Thm 9 (quantum Singleton bound): if k = n-2d+2 the minimum nonzero weight of S is n-d+2 and Q is pure, so the equivalence remark above holds
(the "provided pure" proviso is discharged by that theorem, again via the auditor's reading). (4) Their proof: 341 six-solid configurations in PG(7,2), at most 1311
polarity solutions ([[6,2,3]]_4), no 7th solid. The citation for the final write-up is therefore: abstract + Sec. 5.2 of arXiv:2401.06618v2; our §§1-8 are the
independent computational confirmation.
