#!/usr/bin/env python3
"""Independent re-check (stdlib only) of the two ingredients of the Lean-verified AME(9,10) state of
facebookresearch/atlas-lean, Atlas/FC/OpenQuantumProblem35_ame_9_10_open_solution.lean (first commit 2026-08-28),
as compiled in niketp03/atlas-fc-verified, AtlasFCSolutions/Oqp35.lean (commit 15e4b3a, 2026-09-14), read 2026-11-10.

The state: each party is Z_2 x Z_5; amplitude(u, v) is proportional to
    sum over a in Z_2^2 of (-1)^(phase_C9(u) + label(a).u) * zeta_5^(phase_G(q5Logical(a), v)),
phase_C9 the graph phase of the 9-cycle, label(a) = a0*l1 + a1*l2 with l1, l2 the two check vectors below, G the
symmetric circulant on Z_10 over Z_5 with c = (3,1,1,0,0), and q5Logical(a) = a0 + 2 a1 in {0,1,2,3} the value of
vertex 0 of G.  In the language of the general coupling lemma (paper, Lemma 3): k = 4, a = 2, K = 4, b = 5;
|A_t> = Z^label(t) |C_9>, |G_b> = |G>.

Checks: (1) |G> is AME(10,5): G[A^c, A] is invertible mod 5 for all 252 five-sets A.
        (2) the uniform mixture of the four |A_t> has maximally mixed 4-party marginals: every nonzero x in Z_2^9
            with l1.x = l2.x = 0 gives a stabilizer element X^x Z^(C9 x) of party weight >= 5.
        (3) brute force on the 4 + 4 parties of a few random cuts is not attempted; (1) and (2) with Lemma 3 suffice.
"""
import itertools


def det_mod(M, p):
    M = [r[:] for r in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d = d * M[c][c] % p
        inv = pow(M[c][c], p - 2, p)
        for i in range(c + 1, n):
            f = M[i][c] * inv % p
            M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    return d % p


def q5(i, j):
    e = abs(i - j)
    e = min(e, 10 - e)
    return 3 if e == 1 else (1 if e in (2, 3) else 0)


bad = 0
for A in itertools.combinations(range(10), 5):
    Ac = [i for i in range(10) if i not in A]
    if det_mod([[q5(i, j) for j in A] for i in Ac], 5) == 0:
        bad += 1
print("AME(10,5) circulant c=(3,1,1,0,0): 252 cuts,", bad, "singular")

l1 = [1 if i in (0, 1, 3, 4, 6, 7) else 0 for i in range(9)]
l2 = [1 if i in (0, 2, 3, 5, 6, 8) else 0 for i in range(9)]
minw, count = 99, 0
for x in itertools.product(range(2), repeat=9):
    if not any(x) or sum(a * b for a, b in zip(l1, x)) % 2 or sum(a * b for a, b in zip(l2, x)) % 2:
        continue
    count += 1
    z = [(x[(i - 1) % 9] + x[(i + 1) % 9]) % 2 for i in range(9)]
    minw = min(minw, sum(1 for i in range(9) if x[i] or z[i]))
print("9-cycle code: nonzero stabilizer elements", count, "minimum party weight", minw)
ok = bad == 0 and minw >= 5 and count == 2 ** 7 - 1
print("AME(9,10) INGREDIENTS OK" if ok else "FAILED")
