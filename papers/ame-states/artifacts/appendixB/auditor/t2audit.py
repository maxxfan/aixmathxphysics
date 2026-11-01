#!/usr/bin/env python3
"""Independent re-run of the computational part of pairs #163 (involution sector t=2 swaps, 4 fixed parties).

Claim (iv): the 4-dim subspaces Z of L = (F_2^2)^4 (4 blocks, b = sum of the standard alternating forms on the
blocks) with minimum BLOCK weight >= 3 (additive (4,16,3)_4 MDS codes) number 432, and all have a 2-dim radical for
b|Z, so none is b-nondegenerate; a relaxed control (block weight >= 2) finds 42264 nondegenerate ones.

Independent method: enumerate ALL [8 choose 4]_2 = 200787 four-dim subspaces of F_2^8 through reduced row echelon
forms (pivot sets x free entries), compute block weights of the 15 nonzero vectors and the rank of the b-Gram matrix.
"""
import itertools, sys
from collections import Counter

def bw(v):  # number of nonzero 2-bit blocks of an 8-bit vector
    return sum(1 for j in range(4) if (v >> (2 * j)) & 3)
def b(u, v):
    t = 0
    for j in range(4):
        uu = (u >> (2 * j)) & 3; vv = (v >> (2 * j)) & 3
        t ^= ((uu & 1) & (vv >> 1)) ^ ((uu >> 1) & (vv & 1))
    return t
def rank2(rows):
    rows = list(rows); r = 0
    for bit in range(8):
        piv = None
        for i in range(r, len(rows)):
            if rows[i] >> bit & 1: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return r

n_sub = 0
stats = {thr: Counter() for thr in (2, 3)}
for piv in itertools.combinations(range(8), 4):
    # RREF with pivot columns piv (bit index = column): row i has 1 at piv[i], 0 at other pivots, free elsewhere to the right of piv[i]
    free = []
    for i, p in enumerate(piv):
        free.append([c for c in range(p + 1, 8) if c not in piv])
    nf = sum(len(f) for f in free)
    for mask in range(1 << nf):
        rows = []; k = 0
        for i, p in enumerate(piv):
            r = 1 << p
            for c in free[i]:
                if mask >> k & 1: r |= 1 << c
                k += 1
            rows.append(r)
        n_sub += 1
        # all 15 nonzero vectors
        mn = 99
        for m in range(1, 16):
            v = 0
            for i in range(4):
                if m >> i & 1: v ^= rows[i]
            w = bw(v)
            if w < mn: mn = w
        if mn >= 2:
            G = [[b(rows[i], rows[j]) for j in range(4)] for i in range(4)]
            rk = rank2([sum(G[i][j] << j for j in range(4)) for i in range(4)])
            for thr in (2, 3):
                if mn >= thr:
                    stats[thr][rk] += 1
print("total 4-dim subspaces enumerated:", n_sub, "(Gaussian binomial [8,4]_2 = 200787)")
for thr in (3, 2):
    tot = sum(stats[thr].values())
    print(f"min block weight >= {thr}: {tot} subspaces; rank of b|Z histogram {dict(sorted(stats[thr].items()))}  "
          f"-> radical dims {dict(sorted((4 - r, c) for r, c in stats[thr].items()))}")
