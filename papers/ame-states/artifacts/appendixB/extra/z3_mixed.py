#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Table VI row l = 3, sigma = (3^2, 1^2), "mixed" twist types.  Omega = sigma * (h_6, h_7), sigma = (0 1 2)(3 4 5), h_j of
order dividing 3 on the fixed parties 6, 7.  Hand derivation: auditor #168 and pairs #191(3), audited #199.
F_4-dimensions (Omega acts on V_nf as an F_4-structure): the cycle parts A_nf, B_nf have dim 4 each, the twist part
N_j = ker(1 + h_j + h_j^2) of party j has dim r_j = 0 (h_j = 1), 1 (h_j of order 3 with a 2-dim fixed space), 2 (h_j
fixed-point-free).  T = S cap V_nf is an F_4-Lagrangian of V_nf (dim 8 + r_6 + r_7), so dim T = (8 + r_6 + r_7)/2 must
be an integer (parity).  A subspace X of V_nf supported on 4 parties (A_nf + N_j or B_nf + N_j, parties {0,1,2,j} or
{3,4,5,j}) meets T in dim >= dim T + dim X - dim V_nf; a nonzero intersection is a stabilizer element of party weight
<= 4, impossible.  This script runs over all 9 (r_6, r_7) and prints which ones die by parity, which by weight, and which
are left for the computational rows of Table VI.
"""
LEFT = {(0, 0): "row h = 1 (77760 isometries)", (1, 1): "row 'h with 2-dim fixed space' (a = 2)",
        (2, 2): "row 'h fixed-point-free' (sector A)"}
for r6 in (0, 1, 2):
    for r7 in (0, 1, 2):
        dimV = 8 + r6 + r7
        if dimV % 2:
            print("(r6,r7) = (%d,%d): dim V_nf = %d odd -> no F_4-Lagrangian: excluded by parity" % (r6, r7, dimV))
            continue
        dimT = dimV // 2
        lbs = {"A+N6": dimT + 4 + r6 - dimV, "A+N7": dimT + 4 + r7 - dimV,
               "B+N6": dimT + 4 + r6 - dimV, "B+N7": dimT + 4 + r7 - dimV}
        worst = max(lbs.values())
        if worst > 0:
            which = [k for k, v in lbs.items() if v > 0]
            print("(r6,r7) = (%d,%d): dim T = %d, dim(T cap X) >= %d for X = %s (4 parties): excluded by weight"
                  % (r6, r7, dimT, worst, ",".join(which)))
        else:
            print("(r6,r7) = (%d,%d): dim T = %d, no forced 4-party word -> %s" % (r6, r7, dimT, LEFT[(r6, r7)]))
