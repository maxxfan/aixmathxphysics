#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Table VI row l = 2, sigma = (2, 1^6): one swapped pair of parties {0, 1}, six fixed parties j = 2..7 carrying local
involutions l_j with r_j = rank(1 + l_j) in {0, 1, 2}.  The exclusion is a hand derivation (pairs #165, re-derived by the
auditor #169); this script runs its counting chain over ALL 3^6 rank patterns and performs the two exhaustive code checks
it ends with.  Notation: Omega^2 = 1, N = 1 + Omega, K = S^Omega = ker N|_S, d = dim K, dim NS = 16 - d, NS subset K.
  (i)   d >= 8                                   (rank-nullity: NS subset K, dim NS = 16 - d)
  (ii)  d <= 12 - r_a - r_b   for all fixed a != b (K -> (w, f_a, f_b) injective; Fix(l_j) has dim 4 - r_j)
  (iii) dim NS_0 >= 12 - d,  NS_0 = NS cap {w = 0}  (dim NS - 4)
  (iv)  dim NS_0 <= r_c + r_e for all fixed c != e  (NS_0 lives on the 6 fixed blocks in (+) im N_j, weight >= 5, so it
        injects into any two blocks)
  (v)   dim NS_0 <= dim NS = 16 - d.
A pattern survives (i)-(v) iff some integers d, k = dim NS_0 satisfy all of them.  Then
  r = 0: NS lies in {(w, w)} on the swapped pair (weight <= 2), so NS = 0 and d = 16, contradicting d <= 12;
  r = 1: NS_0 is a binary code of length 6 (im N_j is 1-dim), dim >= 2, weight >= 5: needs A(6,5) >= 4, but A(6,5) = 2;
  r = 2: NS_0 is an F_2-linear code of F_2-dim 4 on 6 blocks of alphabet F_2^2 with block weight >= 5, i.e. an additive
         OA(16,6,4,2) of index 1 = additive MDS code with (n, k) = (6, 2): none exists (enumeration below; Rao bound
         16 < 1 + 6*3 = 19; Bose-Bush: at most s + 1 = 5 factors).
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classical_codes as cc   # noqa: E402  (release-written enumerations, same folder)

surv = []
for r in itertools.product((0, 1, 2), repeat=6):
    pairs = list(itertools.combinations(range(6), 2))
    maxpair = max(r[a] + r[b] for a, b in pairs)
    minpair = min(r[a] + r[b] for a, b in pairs)
    ok = False
    for d in range(8, 17):                        # (i)
        if d > 12 - maxpair:                      # (ii)
            continue
        for k in range(0, 17 - d):                # (v)
            if k >= 12 - d and k <= minpair:      # (iii), (iv)
                ok = True
    if ok:
        surv.append(r)
print("rank patterns (r_2..r_7) in {0,1,2}^6: 729; surviving the counting chain (i)-(v): %d" % len(surv))
for r in surv:
    print("  survivor r = %s  (d = %d, dim NS_0 = %d)" % ("".join(map(str, r)), 12 - 2 * r[0], 2 * r[0]))
print("r = 0: NS inside {(w,w)} on the swapped pair has party weight <= 2 < 5 => NS = 0 => d = 16 > 12: excluded")
print("r = 1: binary [6,2,>=5] code needed; A(6,5) = %d < 4: excluded" % cc.A_nd(6, 5))
print("r = 2: additive MDS (n,k) = (6,2) over F_2^2 [OA(16,6,4,2), index 1]: %d codes (control (5,2): %d): excluded"
      % (cc.count_additive_mds(6, 2), cc.count_additive_mds(5, 2)))
