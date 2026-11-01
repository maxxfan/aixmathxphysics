# Ladder strategies: raw dimension-versus-value data

`ladders.tsv` lists every rung of the strategy ladders discussed in Sec. V C of the paper (Fig. 1). The
columns are:

| column | meaning |
|---|---|
| `inequality` | catalog name |
| `construction` | bulk-cell insertion (one Schmidt pair, with its coupling block, duplicated in the bulk of a pair chain), ancilla doubling (`A ⊗ 1_m`, `B ⊗ 1_m`, `psi ⊗ Phi_m`, perturbed, re-optimized), partial doubling (only the K Schmidt levels of smallest weight), or a plain seesaw / Schmidt-form iteration |
| `rung` | the step: number of Schmidt pairs, or the doubling stage |
| `local_dim` | local dimension of the strategy (for certified rungs, the dimension of the certificate) |
| `active_rank` | Schmidt rank after compression, where recorded |
| `value` | violation |
| `status` | `certified`: an exact strategy, checked by the frozen verifier on the board and again by `replay/replay.py` (value rounded down to 1e-10). `numerical`: a floating-point seesaw value copied from the agents' notebooks (`notebooks/`). These rows are not certified. |
| `source` | board submission number, or the notebook entry |
| `board_file`, `certificate`, `sha256` | for certified rungs: the frozen submission and its copy in this release |

## J73 (bulk-cell insertion)

From 6 to 16 pairs (local dimension 12 to 32) the gains are

| pairs k | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| gain Δ_k | 2.885e-4 | 1.494e-4 | 7.92e-5 | 4.26e-5 | 2.33e-5 | 1.27e-5 | 7.1e-6 | 3.9e-6 | 2.1e-6 | 1.2e-6 |

Consecutive ratios lie between 0.52 and 0.59, and a least-squares fit of log Δ_k gives 0.545. Summing
the fitted geometric tail beyond k = 16 gives ≈ 1.1223884. The rungs at k = 10 (1.1223367271) and k = 16
(1.1223870369, the paper's lower bound) are certified. The certified upper bound is 1.1224104368.

## J62, J48, J30 and the A list

* J62: 3, 4, 5, 6 and 8 pairs, then several pairs inserted at once up to dimension 32 (0.7500885544, certified).
* J48: the chain peaks at 6 pairs (0.7516219860, Pál and Vértesi's R^12 optimum). Further insertions
  lower the value (7 and 8 pairs). This is a negative result.
* J30: doubling the R^7 optimum twice reproduces Pál and Vértesi's R^17 value 0.4363842240 and
  saturates there. Partial doubling gives nothing higher.
* A68, A89, A62, A82: the doubling sequences behind the four raised A-list lower bounds.

## What these data do and do not show

The ladders are evidence, not proofs. Three claims in the paper rest on them, and all three are conjectural:

1. *The J73 and J62 values converge geometrically.* This holds only for the rungs listed, within one family
   of constructions. A different family, or a strategy outside any ladder, may do better in a fixed finite
   dimension.
2. *The extrapolated limits (≈ 1.1223884 for J73, ≈ 0.7500886 for J62).* These come from a fitted
   geometric tail. They are neither lower nor upper bounds on the quantum maximum.
3. *The maximum is approached only in infinite dimension.* This is an analogy with I3322, for which it has
   been proved. Our data do not establish it.

What a proof would need:

* A global upper bound matching the ladder limit: an SOS certificate in the projector algebra with value
  at most the limit plus the desired tolerance. The present certified gap for J73 is 2.3e-5, so this needs
  a much stronger relaxation, or a certificate whose zero set forces the ladder structure.
* A proof of non-attainment: an argument that no finite-dimensional strategy reaches the supremum, as
  was given for I3322 (Pauwels, arXiv:2608.29734; Coladangelo, arXiv:2609.06038; Douglas,
  arXiv:2609.05555). Nothing of this kind has been done for J73 or J62.
* For the value of a finite-dimensional maximum (for example "the best strategy in dimension 32"): a
  dimension-restricted upper bound, which the NPA hierarchy does not provide.

## Notebooks

`notebooks/lower1_NOTES.md` and `notebooks/lower2_NOTES.md` are the working notebooks of the two strategy
agents, verbatim. The numerical rows of `ladders.tsv` come from them (search for the values). They also
record dead ends: seeds, ranks and constructions that found nothing higher.
