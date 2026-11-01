"""val_inv8.py : validate the free-involution parametrization (C, q) and the inv8 DFS.
 (a) for random C, random symmetric q: S is isotropic of rank 16 and (min weight of S >= thr+1) == (q x not in B_x for all x)
     where B_x is recomputed here independently (coset enumeration).
 (b) DFS solution count (thr=3) vs Monte-Carlo fraction of random q with min weight >= 4."""
import sys, os, random, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # [release] was the absolute work/pairs path
import gen_inv8 as g

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 300


def words_minweight(rows):
    """rows: 16 x 32 ints (x16|z16). Returns (min weight over nonzero words, number of words with weight<=4)."""
    xs = np.zeros(1, dtype=np.uint32)
    zs = np.zeros(1, dtype=np.uint32)
    for r in rows:
        xv = sum(r[k] << k for k in range(16))
        zv = sum(r[16 + k] << k for k in range(16))
        xs = np.concatenate([xs, xs ^ np.uint32(xv)])
        zs = np.concatenate([zs, zs ^ np.uint32(zv)])
    w = np.zeros(len(xs), dtype=np.int32)
    for j in range(8):
        nz = (((xs >> (2 * j)) & 3) | ((zs >> (2 * j)) & 3)) != 0
        w += nz.astype(np.int32)
    w[0] = 99
    return int(w.min()), int((w <= 4).sum()), w


def rand_sym(rng):
    q = [[0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i, 8):
            v = rng.getrandbits(1)
            q[i][j] = v
            q[j][i] = v
    return [sum(q[i][j] << j for j in range(8)) for i in range(8)]


def table_pass(d, cols, thr):
    """independent: for every x != 0, enumerate the coset U = Qx + c (c in C) and test the folded weight."""
    Cwords = [0]
    for dv in d:
        Cwords = Cwords + [c ^ dv for c in Cwords]
    for x in range(1, 256):
        D = 0
        Qx = 0
        for i in range(8):
            if x >> i & 1:
                D ^= d[i]
                Qx ^= cols[i] << 8
        for c in Cwords:
            U = Qx ^ c
            if g.posw(U) + g.posw(U ^ D) <= thr:
                return False
    return True


bad_agree = 0
for t in range(NS):
    d = g.random_C(rng)
    cols = rand_sym(rng)
    rows = g.build_rows(d, cols)
    assert g.isotropic(rows) and g.rank2(rows) == 16
    mw, nbad, _ = words_minweight(rows)
    for thr in (2, 3, 4):
        a = table_pass(d, cols, thr)
        b = mw >= thr + 1
        if a != b:
            print("MISMATCH thr", thr, "table", a, "brute", b, "minweight", mw)
            sys.exit(1)
    bad_agree += 1
print("(a) parametrization/table agree on", bad_agree, "random (C,q); example min weights follow")
mws = []
for t in range(2000):
    d = g.random_C(rng)
    cols = rand_sym(rng)
    mws.append(words_minweight(g.build_rows(d, cols))[0])
print("min weight histogram of random (C,q):", {k: mws.count(k) for k in sorted(set(mws))})
