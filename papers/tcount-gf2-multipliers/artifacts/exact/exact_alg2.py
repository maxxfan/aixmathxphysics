#!/usr/bin/env python3
"""auditor: exhaustive T-count of the 2-dimensional F_2-algebras F_2[x]/(P), deg P = 2, by enumerating the whole coset (2^22 elements, N = 6).
Same method as exact_n2.py (which is the GF(4) case) but P is a command-line argument, so it also covers F_2[y]/(y^2), F_2 x F_2 = F_2[x]/(x^2+x).
usage: exact_alg2.py P0 P1   (P = x^2 + P1 x + P0, e.g. `1 1` = GF(4), `0 0` = F2[y]/(y^2), `0 1` = F2 x F2)
The tensor is S(a,b,c) = c . (a b) in the basis 1, x; T-count is invariant under change of basis and of the pairing (paper Sec. II)."""
import itertools
import sys

import numpy as np

p0, p1 = int(sys.argv[1]), int(sys.argv[2])
n, N = 2, 6
red = [p0, p1]  # x^2 = p0 + p1 x
# x^e mod P for e = 0..2
pw = [[1, 0], [0, 1], red[:]]
tgt = {}
for i in range(n):
    for j in range(n):
        for k in range(n):
            if pw[i + j][k]:
                tgt[tuple(sorted((i, n + j, 2 * n + k)))] = 1
mons = list(itertools.combinations_with_replacement(range(N), 3))
t = np.array([tgt.get(mm, 0) for mm in mons], dtype=np.uint8)
cols = list(range(1, 1 << N))
A = np.zeros((len(mons), len(cols)), dtype=np.uint8)
for c, u in enumerate(cols):
    bits = [(u >> i) & 1 for i in range(N)]
    for r, (i, j, k) in enumerate(mons):
        A[r, c] = bits[i] & bits[j] & bits[k]
M = np.concatenate([A, t[:, None]], axis=1).copy()
rows, ncols = A.shape
piv, r = [], 0
for c in range(ncols):
    p = next((rr for rr in range(r, rows) if M[rr, c]), None)
    if p is None:
        continue
    M[[r, p]] = M[[p, r]]
    for rr in range(rows):
        if rr != r and M[rr, c]:
            M[rr] ^= M[r]
    piv.append(c)
    r += 1
    if r == rows:
        break
for rr in range(r, rows):
    assert M[rr, -1] == 0, "inconsistent system"
free = [c for c in range(ncols) if c not in set(piv)]
print(f"algebra x^2 = {p0} + {p1} x: N={N} columns={ncols} rank={len(piv)} kernel dim={len(free)}")
s0 = np.zeros(ncols, dtype=np.uint8)
for ri, c in enumerate(piv):
    s0[c] = M[ri, -1]


def pack(v):
    x = 0
    for c in np.nonzero(v)[0]:
        x |= 1 << int(c)
    return x


arr = np.array([pack(s0)], dtype=np.uint64)
for f in free:
    v = np.zeros(ncols, dtype=np.uint8)
    v[f] = 1
    for ri, c in enumerate(piv):
        v[c] = M[ri, f]
    arr = np.concatenate([arr, arr ^ np.uint64(pack(v))])
x = arr.copy()
x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
w = ((x * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.int64)
mn = int(w.min())
hist = np.bincount(w)
print(f"coset size {len(arr)}, min weight {mn}, minimisers {int((w == mn).sum())}, histogram",
      {i: int(h) for i, h in enumerate(hist) if h and i <= mn + 4})
s = arr[int(np.argmin(w))]
sel = [cols[c] for c in range(ncols) if (int(s) >> c) & 1]
chk = {}
for u in sel:
    bits = [(u >> i) & 1 for i in range(N)]
    for mm in mons:
        if bits[mm[0]] & bits[mm[1]] & bits[mm[2]]:
            chk[mm] = chk.get(mm, 0) ^ 1
print("minimiser re-verified against the target monomials:", all(chk.get(mm, 0) == tgt.get(mm, 0) for mm in mons), "| factors:", len(sel))
