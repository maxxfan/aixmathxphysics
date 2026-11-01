<!-- Copied verbatim from results/ame84_f2linear_closure.md of the p10 repository for the Appendix-B release (sector table: who proved it, who re-checked it). -->

# F_2-linear stabilizer AME(8,4): consolidated status of the symmetry classification (auditor, 2026-10-16 04:15)

**Statement under audit.** An F_2-linear stabilizer AME(8,4) (ququart = 2 qubits, S an isotropic F_2-subspace of
F_2^32 of dimension 16, no nonzero element supported on <= 4 parties) has a TRIVIAL automorphism group in
GL(4,2) wr S_8 (equal to its symplectic automorphism group in Sp(4,2) wr S_8 by pairs #191(1)).
This is a negative result about symmetric ansaetze, NOT a nonexistence proof: a stabilizer AME(8,4) with trivial
automorphism group is not excluded. No candidate was submitted and nothing here is a win.

## Reduction to prime order
A nontrivial automorphism has a power of prime order l in {2,3,5,7} (|GL(4,2) wr S_8| = 20160^8 * 8!; primes
2,3,5,7 only). So it suffices to exclude each prime order. Z_{p^2} (non-elementary) stabilizers are excluded by the
Z_4 lemma (pairs #191(2), additive #197; audited #199), so the F_2-linear qubit-pair family is the only Pauli-type
family at d = 4.

## Sector table (who proved it, who re-checked it with independent code)

| sector | claim | first derivation | independent re-check | auditor files / argument.md |
|---|---|---|---|---|
| local (permutation-free) automorphism, any odd or even prime, fixed points allowed | none | auditor #118/#119 (corrected #132) | consistent with pairs #129 use (not independently re-derived by pairs); non-symplectic case reduced to it by pairs #191(1) (audited #199) | `fixpt_count.py`, 5.3 |
| uniform fpf order 3, F_4-structure H, S types | no state | pairs (dfs4.c) | auditor `hdfs.c` (H, #87), `hsdfs.c` (S, #99): independent DFS, calibrated n=4,6,8 | 5, 5.2 |
| uniform M type, mixed H/S/M | excluded by chain-ring / Singleton argument | auditor 5.2 (#103/#106) and 5.3 closure | single-author (the l = 3 case of 5.3 covers all fpf order-3 local maps) | 5.2, 5.3 |
| involution, t = 4 swaps, f = 0 | none | pairs #129 (37 AME(4,4) classes) | auditor: AME(4,4) stabilizers = ONE orbit under Sp(4,2)^4 x S_4 (`ame44_classes.py`), 250 Hermitian + ~146 random codes, 0 solutions (`inv4*.py`) | 5.4, 5.5 |
| involution, t = 3 | none | pairs #161 (8640 K DFS) | auditor `t3k0audit.py` (12 form pairs = one orbit, 720 K up to gauge) + `t3lift.py` (own lift DFS: 0 solutions for all 10 orbit reps; THR = 4 control 1,433,150 solutions) | 5.7 |
| involution, t = 2 | none | pairs #163 | auditor `t2audit.py`: all 200787 four-dim subspaces of F_2^8, exactly 432 with block weight >= 3, all with 2-dim radical | 5.7 |
| involution, t = 1 | none | pairs #165 (counting) | auditor #169 re-derived by hand (counting forces all twist ranks r equal; r = 0, 1, 2 each impossible; r = 2 via the OA(16,6,4,2) Rao bound) | #169 |
| 7-cycle (+ fixed party) | none | auditor `c7sym.py` (all 5917 F_8-subspaces of F_8^4; best min party weight 4) | pairs #175 not needed; counting variant also in 5.4 | 5.4 |
| 5-cycle with order-5 twists | none | pairs #175 (`c5twist.py`, 331,500 candidates) | auditor `c5audit.py`, `c5graph.py`: same histogram {2:4500, 3:54000, 4:273000}; 333,125 = (q+1)(q^3+1)(q^5+1) generators of H(5,16) | 5.6 |
| single 3-cycle, order-3 twists | none | pairs #175 (counting) | auditor `b3audit.py`: 0 of 243 patterns survive | 5.6 |
| (3,3,1,1), twists h6 = h7 = 1 | none | pairs #138 | auditor `z3_audit.py` (#146): 77760 = |GU(4,2)| equivariant isometries, same low-weight word counts | 5.5 |
| (3,3,1,1), both twists with 2-dim fixed space (a = 2) | none | auditor #168 (`z3twist2.py`, `lagsearch2.c`: 1152 T, 0 Lagrangians; positive control `lag_control.py`) | **CONFIRMED by pairs #238 (04:24, own code `a2_audit.py` / `a2_stage3.py`, no auditor code read)**: 270 X_A lines (one orbit under Sp_A x C(h6) x C(h7); my 3 classes of 90 under Sp_A alone), 384 T completions by an exhaustive DFS (88,863 nodes), |Vs| = 2928 of 4096 for every T, 0 Lagrangians for all 384 T, THR = 4 control finds L for 40/40. Earlier note, now superseded: pairs had re-derived only the parity part by hand (#191(3)). Within-author cross-check of the last stage: the generic unstaged DFS `lagsearch.c` (5,735,899 nodes per T, no use of the dim-2 structure) gives 0 Lagrangians on T = 0-6, 384-386, 768-770 (all three X_A classes), agreeing with the staged `lagsearch2.c` (3565 nodes per T, all 1152 T). First structural step re-derived by hand (04:20): X_A = F_4-line {x, Omega x, Omega^2 x}, x = (a0,a1,a2; c; d) with a_i nonzero distinct nibbles summing to 0 (weight 3), c in N6 minus 0 (fixed to one basis vector by the F_4 scalar), d in N7 minus 0 (3 keys); omega(x, Omega x) = omega(a0,a1) + 1 + 1 (c, d each contribute 1: N_j is a nondegenerate F_4-line), so isotropy iff omega(a0,a1) = 0: 15 x 6 = 90 pairs per key, transitive under Sp(4,2) (stabilizer of u acts as S_3 on u-perp/u and the transvection by u swaps within a coset), total 3 x 90 = 270 = the 270 X_A / 3 classes of 90 of #168. The upstream T-completion count (891 Lagrangians of the Hermitian quotient, 384 of them with T-alone weight >= 5) is still single-author | 5.5(c) |
| (3,3,1,1), mixed twist types (r6, r7) = (2,0), (1,0), (2,1), ... | none | auditor #168 and pairs #191(3) by parity / weight argument | two hand derivations that agree | #168, #191 |
| (3,3,1,1), both twists fixed-point-free (sector A) | none | pairs #186 (3 T_A orbits, 391,910,400 pairs, 0 valid g) | auditor #202: `z3twist.py`, `z3twist5.py`, `z3twist6.py` (233280 X_A, 7 pair classes, completions 24/24/24/18/16/16/16, f-test 0 hits; direct-S control 21/21 agree) | 5.8 |

Order-3 elements of GL(4,2) wr S_8 have k = 1 or 2 three-cycles (n = 8 leaves at most two); k = 1 is the "single
3-cycle" row, k = 2 is the (3,3,1,1) rows, k = 0 is the local row.

Positive controls (a sector search that returns solutions when the weight threshold is relaxed from 5 to 4) exist for
every computational row except sector A step 3 of the auditor's own pipeline, where the substitute control is the
direct 2^16-word verification of 21 unreduced pairs plus the near-miss (bad_cuts 3, min weight 4).

## Open or single-author items (what a further reviewer should look at)
1. CLOSED at 04:24: the `a = 2` sector (#168) now has a second, independent implementation (pairs #238, all
   numbers agree). Every computational sector therefore has two independent implementations; what remains
   single-author are the structural/counting lemmas (5.3, 5.2 for mixed H/S/M types), each reviewed by hand by the
   other agent. The structure it relies on
   (X_A of F_4-dim exactly 1, 270 solutions, 3 classes, 384 completions each) is easy to compare with a count.
2. Theorem 5.3 (trivial local automorphism group) is single-author (auditor) and has no independent re-derivation by pairs, only use; pairs #129/#191 use it, and its counting
   inputs (Singleton with unequal alphabets, OA(256,8,4,4) via Bush, binary A(8,5) = 4) are classical.
3. The sector list is by prime order; elements of composite order are covered because a power has prime order, but
   nobody has independently checked that the list of cycle types of a prime-order element of GL(4,2) wr S_8 is
   exhaustive. For l = 3: 3-cycles k = 1, 2 with fixed-party twists (k = 2 is (3,3,1,1) above; k = 1 is the single
   3-cycle). l = 5: one 5-cycle. l = 7: one 7-cycle. l = 2: t = 1..4 swaps. Pure local elements: 5.3. This is the
   complete list for n = 8.
4. Heuristic context (pairs #198): first moment N_sol/|G| ~ 2^-16.7 in the random-Lagrangian model, so the only
   surviving escape is a trivial-Aut Lagrangian with probability about 2^-17 per class. This is evidence, not proof;
   pairs themselves flag it as weak (runs/pairs/WRITEUP.md, #240): the same model gives about 2^-132 orbits at n = 10
   although AME(10,4) exists, so the real evidence is the closure of every nontrivial-Aut sector, not the moment.

## Scope of "stabilizer"
- Z_4 Weyl-Heisenberg stabilizers: empty (pairs #191(2), additive #197, audited #199).
- Any group-type (Klappenecker-Roetteler) nice error basis, abelian or NONABELIAN, with |K| = 4^8: the index group of
  every party must be abelian and equal, hence Z_2^4 (this note) or Z_4^2 (empty) (pairs #208, audited #224: bijection
  pi_S for every 4-set S, party-5 homomorphism argument, cancellation in finite abelian groups).
- Not covered: non-group unitary error bases (Hadamard / Latin-square types), |K| < 4^8, non-stabilizer states.

## PHASE 2 UPDATE (auditor, 2026-10-16 ~23:00): exhaustive frame search, two independent implementations
classify #315/#320 and the auditor (#318, #319, #322-#324) independently searched the V-picture of a stabilizer AME(8,4)
(S = {(u, Mu)}; V_j = M(A_j) four mutually orthogonal nondegenerate 4-dim subspaces of B with V_J + W_L = B for |J|+|L| = 4;
reduction re-derived in argument.md 5.12). Both find **no solution in any of the 10 V_0 gauge classes** (auditor also in all 72
Sp(4,2)-only classes), with agreeing intermediate numbers (V_0 class orbit sizes, 24128-28960 valid V_1 per class, n=4 count
4478976000, n=6 per-class completions {112, 64, 288, 368, 2272}). Auditor code: `work/auditor/{vcount,vclasses8,vsearch8,run_v8,
gauge_test8,vbrute8,enumtest}`; controls: enumerator exact, gauge invariance 30/30, relaxed-threshold run finds millions (brute-force
verified min party weight 4), near-misses have A_4 >= 6. The written review of `runs/classify/PROOF.md` (argument.md 5.13, post #327,
no gap found) completes the audit: **no F_2-linear (2 qubits per party) stabilizer AME(8,4) exists, with or without symmetry**
(computer-assisted, two independent implementations); together with the Z_4 and group-NEB exclusions below, no stabilizer-type
AME(8,4) remains. This does not touch non-stabilizer states, so the cell AME(8,4) stays "open" in status.json. The symmetry
classification above is superseded for AME(8,4) (kept as an independent cross-check of the table's structure).

**PRIOR LITERATURE (auditor, 2026-10-16 23:45; supersedes the "appears to be new" wording of literature.md sec 7):** the same
nonexistence is already published. S. Ball, E. Moreno, R. Simoens, "Stabiliser codes over fields of even order", arXiv:2401.06618
(v1 12 Jan 2024, v2 9 Sep 2024; the authors' GitHub README says IEEE Trans. Inf. Theory 2024, not confirmed by us) prove that
there is no [[7,1,4]]_4 and no [[8,0,5]]_4 stabiliser code, and state that this "rules out the existence of a stabilised
absolutely maximum entangled state on 8 quqits with local dimension 4" (citing the Huber-Wyderka table). Their model is ours:
Pauli group labelled by F_{2^h} (h = 2: two qubits per party, Thm 21 bijection with binary [[2n,2k]]_2 codes), additive = F_2-linear,
NOT F_4-linear. Method: computer-assisted geometry (7 solids of PG(11,2) with symplectic polarities; classify all [[6,2,3]]_4
codes, 341 classes, at most 1311 polarity solutions; every extension to a 7th solid contradicts full rank of X_2; Lemma 49 projects
[[8,0,5]]_4 to [[7,1,4]]_4; Sage code plus a GAP/FinInG reimplementation giving the same results). They also note that Z/4Z-type
codes with these parameters would imply binary [[8,0,5]]_2, which does not exist. So our result (two independent V-picture
implementations + PROOF.md review) is an INDEPENDENT CONFIRMATION by a different route, not a new theorem; any write-up must cite
2401.06618 and say "confirmed". The cell AME(8,4) stays open for non-stabilizer states only.
*Update 2026-10-17 00:05 (auditor shift 20):* the statements above were re-verified against the FULL TEXT of v2 (PDF text-extracted,
work/auditor/bms/; audits.json key literature_bms_2401_06618_fulltext): the abstract states both nonexistences, Sec 5.2 derives
[[8,0,5]]_4 from [[7,1,4]]_4 by Lemma 49, and [[7,1,4]]_4 is pure because it meets the quantum Singleton bound (Thm 9). Their
"stabiliser code" is the joint eigenspace of an abelian subgroup of the F_q-labelled Pauli group (additive codes, Thm 17), i.e. every
2-qubits-per-party stabilizer state. Recommended citation: Sec 5.2 and the abstract of arXiv:2401.06618v2.

*Update 2026-10-17 00:30 (auditor shift 22):* classify's drop-in paper subsection (`runs/classify/paper_subsection_classification.tex`) was reviewed (argument.md 5.18, board #413/#416/#419): numbers recounted with new code, mathematics correct; the sentence calling the null result 'typical, no structural explanation' must go (n = 6 calibration of the independence heuristic: observed 0.056 vs predicted 0.029; n = 8: 0 of 24 234 vs ~220 predicted). Classify (#414) is cross-checking my 24 234 level-2 triples with frames3 JMASKS=124; per-class numbers posted in #419.

## What is NOT covered by any of the above
- (STALE before phase 2; superseded by the PHASE 2 UPDATE above: the V-picture frame search IS exhaustive over all F_2-linear
  stabilizer states, with or without symmetry.) Stabilizer states with trivial automorphism group (generic Lagrangians of F_2^32)
  were thought to have no exhaustive method (Lagrangian Grassmannian of size prod (2^i + 1), i = 1..16, about 2^136).
- A third, method-independent argument is NOT available: the quadratic-form-refined MacWilliams LP, the weight-5 sub-code form
  decomposition and the near-miss analysis (argument.md 5.14, #338/#340) give no obstruction; the closure rests on two exhaustive
  V-picture implementations (classify #325, auditor #322/#327) plus the written PROOF.md review (a SAT cross-check by classify, #336, is in progress at 23:15).
- Non-stabilizer AME(8,4) states (numerics, #205: pure-permutation symmetric floors f = 0.0364, 0.0395 are saddles; generic floor 0.047).
- The semantics of codetables.de q = 4 bound 4 for [[8,0,d]]_4 (undecided; literature.md).
