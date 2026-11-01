#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Calibration count quoted in Prop. prop:loc: "34 992 of the 77 760 unitary 4x4 matrices over F_4" are H-type AME(4,4)
matrices, i.e. M M^dagger = I (dagger = transpose + Frobenius x -> x^2) and every party block minor nonsingular
(party = 2 consecutive coordinates; 2x2 blocks M[party r, party c] and the full 4x4 matrix).  The auditor's hform.py
counts the 34 992 with a pruned DFS; this script enumerates ALL unitary matrices first (no pruning), checks that their
number is |GU(4,2)| = 2^6 * 3 * 3 * 9 * 15 = 77 760, and then filters.
"""
import itertools

MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [0, 1, 3, 2]
CONJ = [0, 1, 3, 2]


def herm(u, v):
    s = 0
    for a, b in zip(u, v):
        s ^= MUL[a][CONJ[b]]
    return s


def nonsingular(M):
    M = [row[:] for row in M]
    n = len(M)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c]), None)
        if p is None:
            return False
        M[c], M[p] = M[p], M[c]
        iv = INV[M[c][c]]
        for r in range(c + 1, n):
            if M[r][c]:
                f = MUL[M[r][c]][iv]
                M[r] = [x ^ MUL[f][y] for x, y in zip(M[r], M[c])]
    return True


vecs = list(itertools.product(range(4), repeat=4))
units = [v for v in vecs if herm(v, v) == 1]
unitary = []


def dfs(rows):
    if len(rows) == 4:
        unitary.append([list(r) for r in rows])
        return
    for v in units:
        if all(herm(v, r) == 0 for r in rows):
            rows.append(v)
            dfs(rows)
            rows.pop()


dfs([])
gu = 2 ** 6
for i in range(1, 5):
    gu *= 2 ** i - (-1) ** i
good = 0
for M in unitary:
    ok = all(nonsingular([M[2 * r + i][2 * c:2 * c + 2] for i in range(2)]) for r in range(2) for c in range(2)) \
        and nonsingular(M)
    good += ok
print("unit vectors of F_4^4 (h(v,v) = 1): %d" % len(units))
print("unitary 4x4 matrices over F_4: %d   |GU(4,2)| = 2^6 * prod_{i=1..4} (2^i - (-1)^i) = %d" % (len(unitary), gu))
print("with all party block minors nonsingular (H-type AME(4,4) in systematic form): %d" % good)
print("normalised count 34992 / 18^3 = %s" % (good / 18 ** 3))
