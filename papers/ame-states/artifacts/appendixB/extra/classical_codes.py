#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Exhaustive enumerations for the classical coding facts that the proofs of Appendix B quote from the literature.  Each
nonexistence count comes with a positive control at nearby parameters where the object exists.

 (1) No F_4-linear [8,4,5]_4 MDS code (used for (2,2)^8 in Prop. prop:loc and the M-type socle argument; MS77).
     A systematic generator [I_4 | A] is MDS iff every square submatrix of A is nonsingular; rows and columns of A can be
     scaled so that the first row and column are all 1.  Control: [6,3,4]_4 (the hexacode) exists.
 (2) Binary codes of length 8 and minimum distance 5 have at most 4 words, A(8,5) = 4 (used for (1,3)^8 at l = 7 and for
     rank pattern (1)^8 at l = 2); so a binary LINEAR [8,k,5] code has k <= 2.  Also A(6,5) = 2 (row (2,1^6), r = 1:
     a binary [6,2,>=5] code would need 4 words).  Maximum-clique search on the distance graph.
 (3) No F_2-linear (additive) MDS code over the alphabet F_2^2 with (n, k) = (8, 4), i.e. no additive OA(256,8,4,4) of
     index 1 (rank pattern (2)^8 at l = 2; Bush's bound gives at most t + 1 = 5 factors), and none with (n, k) = (6, 2),
     i.e. no additive OA(16,6,4,2) (row (2,1^6), r = 2; Rao: 16 < 1 + 6*3 = 19; Bose-Bush: at most s + 1 = 5 factors).
     An additive code of F_2-dimension 2k in (F_2^2)^n with minimum party weight n - k + 1 projects bijectively onto any
     k parties, so it is the row space of [I_2k | M] with M a k x (n-k) array of 2 x 2 binary blocks, and it is MDS iff
     every square block submatrix of M is nonsingular over F_2.  Local GL(2,2) changes of basis normalise the first block
     row and column to the identity.  Controls: (n, k) = (6, 3) (additive hexacode) and (5, 2) exist.
     The proofs only need F_2-linear codes here: (Omega - 1)S and NS_0 are F_2-subspaces.
"""
import itertools
import sys

# ---------------------------------------------------------------- F_4 (0, 1, w = 2, w^2 = 3; addition = xor)
MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [0, 1, 3, 2]


def det_f4(M):
    M = [row[:] for row in M]
    n = len(M)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c]), None)
        if p is None:
            return 0
        M[c], M[p] = M[p], M[c]
        iv = INV[M[c][c]]
        for r in range(c + 1, n):
            if M[r][c]:
                f = MUL[M[r][c]][iv]
                M[r] = [x ^ MUL[f][y] for x, y in zip(M[r], M[c])]
    return 1


def all_minors_f4(A, k, m):
    for s in range(1, min(k, m) + 1):
        for R in itertools.combinations(range(k), s):
            for C in itertools.combinations(range(m), s):
                if not det_f4([[A[r][c] for c in C] for r in R]):
                    return False
    return True


def count_f4_mds(n, k):
    """number of k x (n-k) matrices A over F_4 with first row and column all 1 and all square submatrices nonsingular"""
    m = n - k
    cnt = 0
    for vals in itertools.product((1, 2, 3), repeat=(k - 1) * (m - 1)):
        A = [[1] * m] + [[1] + list(vals[i * (m - 1):(i + 1) * (m - 1)]) for i in range(k - 1)]
        if all_minors_f4(A, k, m):
            cnt += 1
    return cnt


# ---------------------------------------------------------------- binary codes: A(n, d) by maximum clique
def A_nd(n, d):
    words = [w for w in range(1, 1 << n) if bin(w).count("1") >= d]      # WLOG the code contains 0
    best = [1]

    def grow(clique, cands):
        if len(clique) + 1 > best[0]:
            best[0] = len(clique) + 1
        for i, w in enumerate(cands):
            if len(clique) + 1 + (len(cands) - i) <= best[0]:
                return
            grow(clique + [w], [v for v in cands[i + 1:] if bin(v ^ w).count("1") >= d])
    grow([], words)
    return best[0]


def max_linear_dim(n, d):
    best = 0
    for k in range(1, n + 1):
        found = False
        for gens in itertools.combinations(range(1, 1 << n), k):
            span = {0}
            for g in gens:
                span |= {x ^ g for x in span}
            if len(span) == 1 << k and min(bin(x).count("1") for x in span if x) >= d:
                found = True
                break
        if not found:
            break
        best = k
    return best


# ---------------------------------------------------------------- additive (F_2-linear) codes over F_2^2: block minors
GL22 = []
for a, b, c, dd in itertools.product((0, 1), repeat=4):
    if (a * dd + b * c) % 2:
        GL22.append(((a, b), (c, dd)))
I2 = ((1, 0), (0, 1))


def rank_f2(rows):
    rows = list(rows)
    r = 0
    for bit in reversed(range(64)):
        p = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def block_nonsingular(M, R, C):
    rows = []
    for r in R:
        for i in range(2):
            v = 0
            for j, c in enumerate(C):
                blk = M[r][c]
                v |= (blk[i][0] << (2 * j)) | (blk[i][1] << (2 * j + 1))
            rows.append(v)
    return rank_f2(rows) == 2 * len(R)


def count_additive_mds(n, k):
    """normalised k x (n-k) block arrays (first block row/column = I) with every square block minor nonsingular"""
    m = n - k
    M = [[None] * m for _ in range(k)]
    for c in range(m):
        M[0][c] = I2
    for r in range(k):
        M[r][0] = I2
    cells = [(r, c) for r in range(1, k) for c in range(1, m)]
    count = [0]

    def ok_new(r, c):
        # every square block minor whose last row is r and last column is c, using only filled cells
        for s in range(1, min(r, c) + 2):
            for R in itertools.combinations(range(r), s - 1):
                for C in itertools.combinations(range(c), s - 1):
                    if not block_nonsingular(M, R + (r,), C + (c,)):
                        return False
        return True

    for r in range(k):           # minors among the normalised first row / column
        for c in range(m):
            if (r == 0 or c == 0) and not ok_new(r, c):
                return 0

    def dfs(i):
        if i == len(cells):
            count[0] += 1
            return
        r, c = cells[i]
        for g in GL22:
            M[r][c] = g
            if ok_new(r, c):
                dfs(i + 1)
        M[r][c] = None
    dfs(0)
    return count[0]


if __name__ == "__main__":
    print("(1) F_4-linear MDS, normalised systematic [I_k | A] (first row/column of A all 1):")
    print("    [8,4,5]_4: %d matrices A" % count_f4_mds(8, 4))
    print("    control [6,3,4]_4 (hexacode): %d matrices A" % count_f4_mds(6, 3))
    print("(2) binary codes, maximum size by clique search:")
    print("    A(8,5) = %d   max dim of a binary linear [8,k,>=5] code: %d" % (A_nd(8, 5), max_linear_dim(8, 5)))
    print("    A(6,5) = %d   (a binary [6,2,>=5] code needs 4 words)" % A_nd(6, 5))
    print("(3) F_2-linear MDS codes over F_2^2, normalised block arrays [I | M]:")
    print("    (n,k) = (8,4)  [OA(256,8,4,4), index 1]: %d" % count_additive_mds(8, 4))
    print("    (n,k) = (6,2)  [OA(16,6,4,2), index 1]: %d" % count_additive_mds(6, 2))
    print("    control (n,k) = (6,3) [additive hexacode]: %d" % count_additive_mds(6, 3))
    print("    control (n,k) = (5,2): %d" % count_additive_mds(5, 2))
    print("    bounds: Bush (t >= s, index 1): k <= t + 1 = 5 < 8 factors; Rao for OA(16,6,4,2): 1 + 6*3 = 19 > 16;"
          " Bose-Bush (index 1, t = 2): k <= s + 1 = 5 < 6 factors")
