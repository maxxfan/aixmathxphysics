# F_2-linear stabilizer AME(8,4): exhaustive frame search, per V_0 class (t = 5, all conditions V_J ∩ W_L = 0, |J|+|L| = 4)

Binary `work/classify/frames3` (= frames.c + diagnostics, no env vars set), log `work/classify/final_t5.log`; the first run with frames.c
(`ame84_*.log`) has identical node counts. nodes = (V_0, #V_1 valid, #(V_0,V_1,V_2) valid, #(V_0..V_3) valid).

| class | V_0 tuple-orbit | ψ_1 ψ_2 ψ_3 | #V_1 | #V_2 | #V_3 | free dim L1 / L2 | time |
|---|---|---|---|---|---|---|---|
| 0 | 80640 | 5673 4804 5160 | 26556 | 0 | 0 | 32 / 16 | 28 s |
| 1 | 60480 | 5673 4684 4740 | 24128 | 0 | 0 | 32 / 16 | 28 s |
| 2 | 40320 | 5673 5160 4685 | 24128 | 0 | 0 | 32 / 16 | 28 s |
| 3 | 20160 | 4680 5673 4804 | 28960 | 0 | 0 | 32 / 16 | 27 s |
| 4 | 20160 | 4680 5673 5160 | 24488 | 0 | 0 | 32 / 16 | 27 s |
| 5 | 20160 | 5673 4804 4748 | 28960 | 0 | 0 | 32 / 16 | 29 s |
| 6 | 13440 | 5673 5160 4681 | 24488 | 0 | 0 | 32 / 16 | 27 s |
| 7 | 10080 | 4680 4804 5770 | 28960 | 0 | 0 | 32 / 16 | 28 s |
| 8 | 2016 | 4680 4804 4748 | 28656 | 0 | 0 | 32 / 16 | 27 s |
| 9 | 1344 | 4680 4680 4804 | 28656 | 0 | 0 | 32 / 16 | 25 s |

Every V_1 had free dimension exactly 16 at level 2 (histograms `L2{16:#V_1}` in the log), so no degenerate linear system was skipped.
Total ≈ 4.6 min CPU. **0 AME frames, i.e. no F_2-linear stabilizer AME(8,4).**

## Level 2 without the four J = {0,1,2} conditions (added shift 18, 00:25; `JMASKS=124`, log `work/classify/m124_all.log`, driver `work/classify/run_m124_all.py`)
Same binary source (`frames3_s8`), all conditions except the four with J = {0,1,2} (= the V_3 single-party conditions by duality), ten classes, 22 s each:
#(V_0,V_1,V_2) = 702, 176, 416, 1776, 712, 3640, 2444, 1712, 8912, 3744 (classes 0..9), **total 24,234**, and #V_3 = 0 in every class (V_3 forced, then tested).
Auditor's independent count (argument.md 5.14(c)/5.18): 24,234 triples, per-class range 176..8912 -- identical total. The independence model for the four omitted
conditions (0.3076^4 = 0.9 % each triple) would leave ~220 survivors, observed 0: the V_3 failure is a pattern, not a typical null (auditor, 5.18).

## Nearest misses (by-product, relaxed run)
`JMASKS=64479` (impose everything except the six (2,2) conditions with A-set {0,2} and their duals {1,3}) on class 9 yields leaves at once
(`work/classify/relax_leaf.log`, `sol_relax_m4.txt`): genuine stabilizer states, exact party-weight histogram e.g. [1,0,0,0,30,720,4800,20880,39105],
min party weight 4, **bad_cuts = 1 of 35** (verifier's rankcheck). In that state the single bad cut has F_2 rank deficiency 4 (30 = 2·15 words of
party weight 4, supported on the cut's two 4-sets), so its total cut deficit in the pairs measure (sum of n - rank) is 4, not better than pairs'
anneal best of 2; it is just the cheapest way to get leaves at this level. Frames with bad_cuts 1,2,3,4 appear in the first 12. No file submitted (not a win).
