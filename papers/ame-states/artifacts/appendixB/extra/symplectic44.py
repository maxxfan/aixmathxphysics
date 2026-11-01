#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Brute-force calibration of the auditor's S-type search (auditor/hsdfs.c) at n = 4, the S-type analogue of the
"34 992 of the 77 760 unitary matrices" calibration of the H type.  S-type AME(4,4) in systematic form: S = rowspace
[I_4 | M] over F_4 with M J M^T = J, J = diag([[0,1],[1,0]], [[0,1],[1,0]]) (the F_4-bilinear alternating form of
hsdfs.c), and every party block minor of M nonsingular (party = coordinates 2j, 2j+1).  hsdfs.c reports N_norm = 2 with
#solutions = N_norm * 60^h * 60^(h-1) = 2 * 60^3 = 432 000 at h = 2.  This script enumerates all of Sp(4,4) (rows forming
hyperbolic frames; |Sp(4,4)| = 4^4 (4^2 - 1)(4^4 - 1) = 979 200) and filters by the minors.
"""
import itertools

MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [0, 1, 3, 2]


def form(a, b):
    s = 0
    for c in (0, 2):
        s ^= MUL[a[c]][b[c + 1]] ^ MUL[a[c + 1]][b[c]]
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


J = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
vecs = [list(v) for v in itertools.product(range(4), repeat=4) if any(v)]
count = good = 0
for r0 in vecs:
    for r1 in vecs:
        if form(r0, r1) != 1:
            continue
        rest = [v for v in vecs if form(r0, v) == 0 and form(r1, v) == 0]
        for r2 in rest:
            for r3 in rest:
                if form(r2, r3) != 1:
                    continue
                M = [r0, r1, r2, r3]
                count += 1
                if all(nonsingular([M[2 * r + i][2 * c:2 * c + 2] for i in range(2)]) for r in range(2) for c in range(2)) \
                        and nonsingular(M):
                    good += 1
print("matrices M over F_4 with M J M^T = J: %d   |Sp(4,4)| = 4^4 (4^2-1)(4^4-1) = %d" % (count, 4 ** 4 * 15 * 255))
print("with all party block minors nonsingular (S-type AME(4,4) in systematic form): %d" % good)
print("normalised count %d / 60^3 = %s   (auditor/hsdfs.c 2: N_norm = 2)" % (good, good / 60 ** 3))
