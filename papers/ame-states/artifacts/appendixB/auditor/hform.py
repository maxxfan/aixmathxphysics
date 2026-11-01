"""Independent normal form for H-type AME(n,4) stabilizer states (auditor, calibrates @pairs' Hermitian-h DFS).
S = rowspace [I_n | M] over F_4 (qubit-level coordinates, party = 2 coordinates), parties 0..n/2-1 are the information set.
Hermitian self-duality  <=>  M M^dagger = I  (char 2).   AME(n,4) <=> for every k >= 1 and party sets R (rows), C (cols), |R|=|C|=k,
the 2k x 2k submatrix M[R,C] (party blocks of 2 coordinates) is nonsingular over F_4.
F_4 = {0,1,2=w,3=w^2}: add = xor, conj(x) = x^2.
usage: python3 work/auditor/hform.py n [seconds] [seed]
  n = 4: enumerate all unitary M (77760) and count those satisfying the minors;  n = 6: randomised DFS, count solutions found in the time.
"""
import sys, time, random, itertools
MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [0, 1, 3, 2]
CONJ = [0, 1, 3, 2]


def herm(a, b):
    s = 0
    for x, y in zip(a, b):
        s ^= MUL[x][CONJ[y]]
    return s


def rank(rows):
    rows = [list(r) for r in rows]
    rk, ncol = 0, len(rows[0]) if rows else 0
    for c in range(ncol):
        p = None
        for i in range(rk, len(rows)):
            if rows[i][c]:
                p = i
                break
        if p is None:
            continue
        rows[rk], rows[p] = rows[p], rows[rk]
        iv = INV[rows[rk][c]]
        rows[rk] = [MUL[iv][x] for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]
                rows[i] = [x ^ MUL[f][y] for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def blocks_ok(M, upto, np_):
    """check all minors M[R,C] with R subset of parties 0..upto containing upto, C any parties, |R|=|C|."""
    for k in range(1, upto + 2):
        for Rr in itertools.combinations(range(upto), k - 1):
            R = Rr + (upto,)
            rows = [2 * a + t for a in R for t in (0, 1)]
            for C in itertools.combinations(range(np_), k):
                cols = [2 * b + t for b in C for t in (0, 1)]
                sub = [[M[r][c] for c in cols] for r in rows]
                if rank(sub) < 2 * k:
                    return False
    return True


def main():
    n = int(sys.argv[1])
    secs = float(sys.argv[2]) if len(sys.argv) > 2 else 60
    rnd = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    half = n // 2          # number of info parties = number of other parties
    dim = n                # F_4 coordinates in M
    units = [v for v in itertools.product(range(4), repeat=dim) if herm(v, v) == 1]
    t0 = time.time()
    sols = []
    count_unitary = 0
    nodes = 0

    def dfs(rows):
        nonlocal count_unitary, nodes
        if time.time() - t0 > secs:
            return True
        i = len(rows)
        if i == dim:
            sols.append([list(r) for r in rows])
            return n == 6 and len(sols) >= 1 and False
        cands = units[:]
        if n > 4:
            rnd.shuffle(cands)
        for v in cands:
            if all(herm(v, r) == 0 for r in rows):
                rows.append(v)
                nodes += 1
                if i % 2 == 1:      # a block row (parties) completed
                    if not blocks_ok(rows, i // 2, half):
                        rows.pop()
                        continue
                stop = dfs(rows)
                rows.pop()
                if stop:
                    return True
                if n == 6 and sols:
                    return True
        return False
    dfs([])
    print(f"n={n}: solutions found {len(sols)} (time {time.time()-t0:.0f}s, nodes {nodes}, limit {secs}s)")
    if sols:
        M = sols[0]
        print("first M =", M)
        # verify unitary and all minors again from scratch
        ok = all(herm(M[i], M[j]) == (1 if i == j else 0) for i in range(dim) for j in range(dim))
        ok2 = all(blocks_ok(M, a, half) for a in range(half))
        print("re-check unitary:", ok, " all party minors:", ok2)


if __name__ == "__main__":
    main()
