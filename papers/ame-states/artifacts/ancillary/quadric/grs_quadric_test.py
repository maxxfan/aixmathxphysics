#!/usr/bin/env python3
"""GRS test for the shortened [17,8,10]_256 code of the extend AME(18,16) code (auditor, phase 3).
Generator [I_9 | A] (A circulant, c = argv[1]); shorten at information coordinate 0 -> G' = [I_8 | A[1:, :]] (8 x 17).
Ball (1907.04391, Sec 6, citing Glynn Thm 3.3 and Ball-Pepe Thm 6): for an MDS code of dimension k and length >= 2k+1, the code is GRS iff the
space of quadrics (in the k coordinates) vanishing on the n columns of a generator matrix has dimension C(k-1,2).
Here k = 8, n = 17 = 2k+1: GRS <=> dim = 21 <=> rank of the 17 x 36 evaluation matrix = 36 - 21 = 15; generic points give rank 17.
usage: grs_quadric_test.py c0,...,c8 [row_to_remove]"""
import sys, itertools
MOD = 0x11D
EXP = [0] * 512
LOG = [0] * 256
x = 1
for i in range(255):
    EXP[i] = x
    LOG[x] = i
    x <<= 1
    if x & 256:
        x ^= MOD
for i in range(255, 512):
    EXP[i] = EXP[i - 255]


def mul(a, b):
    return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]


def inv(a):
    return EXP[255 - LOG[a]]


def rank(M):
    M = [r[:] for r in M]
    rk = 0
    for c in range(len(M[0])):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        iv = inv(M[rk][c])
        M[rk] = [mul(iv, v) for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [a ^ mul(f, b) for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


c = [int(t) for t in sys.argv[1].split(",")]
m = len(c)                       # 9
drop = int(sys.argv[2]) if len(sys.argv) > 2 else 0
A = [[c[(j - i) % m] for j in range(m)] for i in range(m)]
rows = [i for i in range(m) if i != drop]
k = len(rows)                    # 8
cols = []
for t in range(k):               # identity columns
    cols.append([1 if s == t else 0 for s in range(k)])
for j in range(m):               # A' columns
    cols.append([A[i][j] for i in rows])
n = len(cols)
mons = list(itertools.combinations_with_replacement(range(k), 2))
E = [[mul(col[i], col[j]) for (i, j) in mons] for col in cols]
r = rank(E)
print(f"k={k} n={n} monomials={len(mons)} rank of evaluation matrix={r}; dim of quadrics through the columns = {len(mons) - r}; "
      f"GRS would give {len(mons) - (len(mons) - (k - 1) * (k - 2) // 2)} (= C(k-1,2) = {(k - 1) * (k - 2) // 2}) i.e. rank {len(mons) - (k - 1) * (k - 2) // 2}; generic: rank {min(n, len(mons))}")
# control: columns of a genuine GRS code = points (1, t, t^2, ..., t^{k-1}) with distinct t in GF(256) plus the point at infinity
pts = [[1 if s == 0 else EXP[(LOG[t] * s) % 255] if t else 0 for s in range(k)] for t in range(1, n)]
pts.append([0] * (k - 1) + [1])
Ec = [[mul(col[i], col[j]) for (i, j) in mons] for col in pts]
rc = rank(Ec)
print(f"control GRS (rational normal curve, n={len(pts)}): rank {rc}, quadric dimension {len(mons) - rc} (expected {(k - 1) * (k - 2) // 2})")
print("VERDICT:", "GRS-equivalent" if r == len(mons) - (k - 1) * (k - 2) // 2 else "NOT GRS (shortened code has no rational-normal-curve structure)")
