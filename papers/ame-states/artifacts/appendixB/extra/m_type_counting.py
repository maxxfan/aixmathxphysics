#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Prop. prop:loc, case (0,4)^8, type M ("Type M makes S a module over a non-commutative chain ring of order 16.  A Singleton
argument on its socle again forces an [8,4,5]_4 MDS code, and a similar count excludes M parties mixed with H or S
parties").  The argument is the auditor's (argument.md 5.2, board #103/#106); it is a hand derivation.  This script runs
its counting over all cases.  Facts used (checked by auditor/f4struct.py): on an M party the algebra A = <w, w*> has 16
elements, radical N with |N| = 4, A/N = F_4, F_2^4 = A free of rank 1; N V_j has 4 elements.

S is an A-module (when R = F_2[Omega, Omega*] = A acts with N = 0 on the a_O = 8 - a_M non-M parties; the other case,
R = A x F_4, splits S as a product across (M | others), impossible because the smaller side has <= 4 parties and a pure
reduced state).  S = A^k1 + (A/N)^k2, 2 k1 + k2 = 8 (F_4-dims).
  (i)  T = N S = N^k1 is an F_4-linear code of F_4-dim k1 on the a_M M-parties (alphabet N V_j = F_4), min weight >= 5:
       k1 > 0 needs a_M >= 5 and k1 <= a_M - 4 (Singleton);
  (ii) S_0 = {s : N s = 0} has 4^(8-k1) elements, alphabets of size 4 on M parties and 16 on the others, min weight >= 5,
       so it injects into any 4 parties: 4^(8-k1) <= product of the 4 smallest alphabets.
  uniform M (a_M = 8): T has F_4-dim 8 - k1 inside N^8 = F_4^8 (it is S cap N V^8 there), weight >= 5, so 8 - k1 <= 4,
       k1 = 4, k2 = 0 and T is an [8,4,5]_4 MDS code -- none exists (classical_codes.py (1)).
"""
print("a_M = number of M-type parties (a_O = 8 - a_M parties of type H or S), k1 = rank of the free part")
for aM in range(1, 8):
    surv = []
    for k1 in range(0, 5):
        small = sorted([4] * aM + [16] * (8 - aM))[:4]
        prod = 1
        for x in small:
            prod *= x
        ok_ii = 4 ** (8 - k1) <= prod
        ok_i = (k1 == 0) or (aM >= 5 and k1 <= aM - 4)
        if ok_i and ok_ii:
            surv.append(k1)
    print("  a_M = %d (mixed with %d H/S parties): k1 surviving (i)+(ii): %s" % (aM, 8 - aM, surv if surv else "none"))
print("  a_M = 8 (uniform M): 8 - k1 <= 4 forces k1 = 4: T is an [8,4,5]_4 MDS code, which does not exist")
print("  a_M = 0: H and S parties only; mixed H/S splits S as a product (hand), uniform H / S are the exhaustive searches")
