#!/usr/bin/env python3
"""Partial independent check of pairs #161 (involution sector t=3 swaps + 2 fixed parties): the K_0 count.

K = S^Omega = N S (8-dim, w-part in (F_2^4)^3, f-part in L_1+L_2).  K_0 = K cap {f=0} = {(x, Mx, Nx)} with M,N in GL(4,2)
and, since Q(k,k') = sum_p <w_p,w_p'> = b(f,f') vanishes on K_0:  M^T J M + N^T J N = J.
Then K = K_0^{perp_Q} (dim 8 = dim K) is determined by K_0, automatically MDS over the 3 blocks when M,N invertible, and
f: K -> L_1+L_2 is an isometry K_0^perp/K_0 -> (L_1+L_2, b_6+b_7): |Sp(4,2)| = 720 choices.
Gauge: local symplectic maps g1,g2,g3 on the three pair positions: (M,N) -> (g2 M g1^-1, g3 N g1^-1), so K_0 classes are the
orbits of Sp(J) (g1) on ordered pairs (F_M, F_N) of nondegenerate alternating forms with F_M + F_N = J.
pairs claim: 12 K_0 classes, 720 K each, 8640 K.
"""
import itertools
from collections import Counter

# 4x4 matrices over F_2 as tuples of 4 row bitmasks
def mat_mul(A, B):
    # (AB)[i] = xor of B[k] for k in bits of A[i]
    out = []
    for i in range(4):
        r = 0
        for k in range(4):
            if A[i] >> k & 1: r ^= B[k]
        out.append(r)
    return tuple(out)
def transpose(A):
    return tuple(sum(((A[j] >> i) & 1) << j for j in range(4)) for i in range(4))
J = (0b0010, 0b0001, 0b1000, 0b0100)  # J[0] has bit1, J[1] bit0, J[2] bit3, J[3] bit2
def det_nonzero(A):
    rows = list(A); r = 0
    for bit in range(4):
        piv = None
        for i in range(r, 4):
            if rows[i] >> bit & 1: piv = i; break
        if piv is None: return False
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(4):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return True
def add(A, B): return tuple(a ^ b for a, b in zip(A, B))

# all alternating 4x4 matrices: symmetric, zero diagonal
alts = []
for m in range(64):
    bits = [(m >> k) & 1 for k in range(6)]
    pos = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    A = [[0] * 4 for _ in range(4)]
    for (i, j), v in zip(pos, bits):
        A[i][j] = A[j][i] = v
    alts.append(tuple(sum(A[i][j] << j for j in range(4)) for i in range(4)))
nondeg = [A for A in alts if det_nonzero(A)]
print("nondegenerate alternating forms:", len(nondeg), "(expect 28)")
assert J in nondeg

# Sp(J) = {g invertible: g^T J g = J}
sp = []
for code in range(1 << 16):
    g = tuple((code >> (4 * i)) & 15 for i in range(4))
    if not det_nonzero(g): continue
    if mat_mul(mat_mul(transpose(g), J), g) == J:
        sp.append(g)
print("|Sp(4,2)| =", len(sp), "(expect 720)")

pairs = [(F, add(F, J)) for F in nondeg if det_nonzero(add(F, J))]
print("ordered pairs (F_M, F_N) with F_M + F_N = J, both nondegenerate:", len(pairs))
pairset = set(pairs)
seen = set(); orbits = []
for p in pairs:
    if p in seen: continue
    orb = set()
    for g in sp:
        gT = transpose(g)
        q = (mat_mul(mat_mul(gT, p[0]), g), mat_mul(mat_mul(gT, p[1]), g))
        assert q in pairset
        orb.add(q)
    seen |= orb
    orbits.append(len(orb))
print("Sp(4,2)-orbits on these pairs:", len(orbits), "sizes", sorted(orbits))
print("=> K_0 classes:", len(orbits), "; times 720 isometries f:", len(orbits) * 720, "K  (pairs: 12 and 8640)")
