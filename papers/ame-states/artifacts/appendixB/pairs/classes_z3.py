"""classes_z3.py : reps of AME(4,4) graph forms (4 positions = A, B, f1, f2) modulo GL(2,2)^4 only (no position permutations).
Writes work/pairs/z3_reps.txt : one line per rep, 8 hex ints d_k = (e_k | Gamma_k) in W layout (bits 0..7 x, bits 8..15 z)."""
import os, itertools, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

HERE = os.environ.get("APPB_WORK", ".")  # [release] was the absolute path of work/pairs; z3_reps.txt is written to the working directory


def rank4(m):
    rows = [(m >> (4 * i)) & 15 for i in range(4)]
    r = 0
    for col in range(4):
        piv = None
        for i in range(r, 4):
            if rows[i] >> col & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(4):
            if i != r and rows[i] >> col & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


inv4 = np.array([rank4(m) == 4 for m in range(65536)], dtype=bool)
tr = np.array([(n & 1) | (((n >> 2) & 1) << 1) | (((n >> 1) & 1) << 2) | (((n >> 3) & 1) << 3) for n in range(16)], dtype=np.uint8)
g = np.arange(1 << 24, dtype=np.uint32)
nib = [((g >> (4 * t)) & 15).astype(np.uint8) for t in range(6)]
del g


def f(n, a):
    return (n >> (2 * a)) & 3


def enc(r0, r1, r2, r3):
    return r0.astype(np.uint32) | (r1.astype(np.uint32) << 4) | (r2.astype(np.uint32) << 8) | (r3.astype(np.uint32) << 12)


n0, n1, n2, n3, n4, n5 = nib
ok = inv4[enc(f(n1, 0) | f(n2, 0) << 2, f(n1, 1) | f(n2, 1) << 2, f(n3, 0) | f(n4, 0) << 2, f(n3, 1) | f(n4, 1) << 2)]
t3, t4, t5 = tr[n3], tr[n4], tr[n5]
ok &= inv4[enc(f(n0, 0) | f(n2, 0) << 2, f(n0, 1) | f(n2, 1) << 2, f(t3, 0) | f(n5, 0) << 2, f(t3, 1) | f(n5, 1) << 2)]
ok &= inv4[enc(f(n0, 0) | f(n1, 0) << 2, f(n0, 1) | f(n1, 1) << 2, f(t4, 0) | f(t5, 0) << 2, f(t4, 1) | f(t5, 1) << 2)]
del nib, n0, n1, n2, n3, n4, n5, t3, t4, t5
print("AME(4,4) graph forms:", int(ok.sum()), flush=True)

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
GL = [np.array(m, dtype=np.int64) for m in ([[1, 0], [0, 1]], [[0, 1], [1, 0]], [[1, 1], [0, 1]], [[1, 0], [1, 1]], [[0, 1], [1, 1]], [[1, 1], [1, 0]])]
Bs = np.array([[GL[i0], GL[i1], GL[i2], GL[i3]] for i0, i1, i2, i3 in itertools.product(range(6), repeat=4)])


def blocks_of(gv):
    G = np.zeros((4, 4, 2, 2), dtype=np.int64)
    for t, (p, r) in enumerate(PAIRS):
        nb = (gv >> (4 * t)) & 15
        blk = np.array([[nb & 1, (nb >> 1) & 1], [(nb >> 2) & 1, (nb >> 3) & 1]], dtype=np.int64)
        G[p, r] = blk
        G[r, p] = blk.T
    return G


def orbit(gv):
    G = blocks_of(gv)
    code = np.zeros(Bs.shape[0], dtype=np.int64)
    for t, (p, r) in enumerate(PAIRS):
        new = np.einsum("nab,bc,ndc->nad", Bs[:, p], G[p, r], Bs[:, r]) % 2
        nibv = new[:, 0, 0] | (new[:, 0, 1] << 1) | (new[:, 1, 0] << 2) | (new[:, 1, 1] << 3)
        code |= nibv << (4 * t)
    return np.unique(code)


visited = np.zeros(1 << 24, dtype=bool)
reps = []
t0 = time.time()
for gv in np.nonzero(ok)[0]:
    if visited[gv]:
        continue
    orb = orbit(int(gv))
    assert ok[orb].all()
    visited[orb] = True
    reps.append((int(gv), len(orb)))
assert sum(s for _, s in reps) == int(ok.sum())
print("orbits under GL(2,2)^4:", len(reps), "(%.1fs)" % (time.time() - t0), flush=True)

with open(HERE + "/z3_reps.txt", "w") as fo:
    for gv, sz in reps:
        G = blocks_of(gv)
        d = []
        for k in range(8):
            z = 0
            for j in range(8):
                if k // 2 != j // 2 and G[k // 2, j // 2, k % 2, j % 2]:
                    z |= 1 << j
            d.append((1 << k) | (z << 8))
        fo.write(" ".join("%x" % v for v in d) + "\n")
print("wrote", len(reps), "reps", flush=True)
