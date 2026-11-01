#!/usr/bin/env python3
"""Audit of pairs #175 (B): single 3-cycle (0 1 2), fixed parties 3..7 with order-3 twists h_j (or 1).
Omega = sigma*h symplectic of order 3, V = V^Omega (+) V_nf, S = S^Omega (+) S_nf, both Lagrangian.
Per fixed party j: u_j = dim Fix(h_j) in {4 (h_j = 1), 2 (2-dim fixed space), 0 (fpf)}, w_j = 4 - u_j.
  V^Omega = {(a,a,a; f_3..f_7)}  : dim 4 + U,     U = sum u_j ; S^Omega Lagrangian: (4+U)/2
  V_nf    = {cycle part (dim 8)} + sum W_j : dim 8 + W, W = sum w_j = 20 - U ; S_nf Lagrangian: (8+W)/2
Necessary counting conditions (all words of S have party weight >= 5 on 8 parties):
  (i)  S_c = {s in S^Omega : a = 0} has dim >= (4+U)/2 - 4; its words live on the 5 fixed parties, so have
       weight exactly 5 (all five fixed components nonzero) and S_c injects into U_j for each j.
  (ii) T_c = {s in S_nf : cycle part 0} has dim >= (8+W)/2 - 8; same argument: injects into W_j for each j.
Enumerate all u in {0,2,4}^5 and report which survive both conditions.
"""
import itertools
surv = []
tot = 0
for u in itertools.product((0, 2, 4), repeat=5):
    tot += 1
    U = sum(u)
    w = [4 - x for x in u]
    W = sum(w)
    dimSO = (4 + U) // 2
    dimSnf = (8 + W) // 2
    assert dimSO + dimSnf == 16
    lbSc = max(0, dimSO - 4)
    lbTc = max(0, dimSnf - 8)
    ok1 = lbSc <= min(u)   # injects into every Fix_j (also trivially ok if lb = 0)
    ok2 = lbTc <= min(w)
    if ok1 and ok2:
        surv.append(u)
print(f"{tot} twist patterns u in {{0,2,4}}^5; survivors of (i)+(ii): {len(surv)}")
for s in sorted(set(tuple(sorted(x)) for x in surv)):
    print("  survivor (sorted):", s)
# What if the injectivity bound is relaxed to 'only needs some party'?  (sanity: the argument really needs all)
