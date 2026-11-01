#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Table VI row l = 5, sigma = (5, 1^3): "Counting forces two twisted fixed parties".  Hand derivation: pairs #175 (C),
re-derived by the auditor #184 / argument.md 5.6.  Omega = sigma * h, sigma = (0 1 2 3 4), fixed parties 5, 6, 7 with
h_j of order dividing 5; an order-5 element of Sp(4,2) is fixed-point-free (F_2[h_j] = F_16; see group_orders.py), so each
fixed party is either twisted (fixed space 0) or untwisted (h_j = 1).  n5 = #twisted, n1 = 3 - n5 = #untwisted.
  (a) V^Omega = {(c,c,c,c,c; f_j on the untwisted parties)} has F_2-dim 4 + 4 n1; S^Omega is a Lagrangian of it,
      dim 2 + 2 n1.  Its words with cycle constant c = 0 are supported on the n1 <= 3 untwisted parties (weight <= 3 < 5),
      so S^Omega injects into c in F_2^4: 2 + 2 n1 <= 4.
  (b) V_nf has F_16-dim 4 + n5 (cycle part F_16^4, one F_16 per twisted party) and S_nf is an F_16-Lagrangian of it
      (Hermitian form), so 4 + n5 is even.
The script enumerates n5 = 0..3; the survivor is the case computed exhaustively by c5twist.py / c5graph.py.
"""
for n5 in range(4):
    n1 = 3 - n5
    a = 2 + 2 * n1 <= 4
    b = (4 + n5) % 2 == 0
    verdict = "SURVIVES (computational row: S_nf is a generator of H(5,16))" if a and b else "excluded"
    why = [] if a else ["(a) dim S^Omega = %d > 4" % (2 + 2 * n1)]
    if not b:
        why.append("(b) F_16-dim of V_nf = %d odd" % (4 + n5))
    print("n5 = %d, n1 = %d: %s %s" % (n5, n1, verdict, "; ".join(why)))
q = 4
print("generators of H(5,q^2), q = 4: (q+1)(q^3+1)(q^5+1) = %d = 1105*300 + 325*5 = %d"
      % ((q + 1) * (q ** 3 + 1) * (q ** 5 + 1), 1105 * 300 + 325 * 5))
