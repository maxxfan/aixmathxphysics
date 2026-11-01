"""classes_inv8.py : enumerate ALL AME(4,4) graph forms C (4 ququart positions = qubit pairs, 24 cross bits), reduce modulo
H = GL(2,2)^4 x S_4 (graph-form-preserving local Cliffords and position permutations), run inv8 (free-involution AME(8,4) DFS) on
one representative per H-orbit.  Every Lagrangian C is LC-equivalent (diagonal local Clifford on (a, a') keeps sigma-symmetry) to a
graph form, so this covers the whole free-involution sector.  Log: work/pairs/inv8_classes.log"""
import os, sys, itertools, subprocess, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

# [release] was ROOT/HERE = absolute work/pairs path; log now in $APPB_WORK (default: cwd), inv8 binary from $BUILD
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = open(os.path.join(os.environ.get("APPB_WORK", "."), "inv8_classes.log"), "a")
INV8 = os.path.join(os.environ["BUILD"], "inv8") if os.environ.get("BUILD") else os.path.join(HERE, "inv8")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


# ---- invertibility table for 4x4 matrices over F_2 (rows r0..r3 as 4-bit fields) ----
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
P("AME(4,4) graph forms among 2^24:", int(ok.sum()))

# ---- orbits under H ----
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
GL = [np.array(m, dtype=np.int64) for m in ([[1, 0], [0, 1]], [[0, 1], [1, 0]], [[1, 1], [0, 1]], [[1, 0], [1, 1]], [[0, 1], [1, 1]], [[1, 1], [1, 0]])]
Bs = np.array([[GL[i0], GL[i1], GL[i2], GL[i3]] for i0, i1, i2, i3 in itertools.product(range(6), repeat=4)])  # (1296,4,2,2)


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
    imgs = []
    for perm in itertools.permutations(range(4)):
        Gp = np.zeros_like(G)
        for p in range(4):
            for r in range(4):
                if p != r:
                    Gp[p, r] = G[perm[p], perm[r]]
        code = np.zeros(Bs.shape[0], dtype=np.int64)
        for t, (p, r) in enumerate(PAIRS):
            new = np.einsum("nab,bc,ndc->nad", Bs[:, p], Gp[p, r], Bs[:, r]) % 2  # B_p G B_r^T
            nibv = new[:, 0, 0] | (new[:, 0, 1] << 1) | (new[:, 1, 0] << 2) | (new[:, 1, 1] << 3)
            code |= nibv << (4 * t)
        imgs.append(code)
    return np.unique(np.concatenate(imgs))


visited = np.zeros(1 << 24, dtype=bool)
reps = []
cand_idx = np.nonzero(ok)[0]
t0 = time.time()
for gv in cand_idx:
    if visited[gv]:
        continue
    orb = orbit(int(gv))
    assert ok[orb].all()
    visited[orb] = True
    reps.append((int(gv), len(orb)))
P("orbits of AME(4,4) graph forms under GL(2,2)^4 x S_4:", len(reps), "orbit sizes sum", sum(s for _, s in reps),
  "(%.1fs)" % (time.time() - t0))
assert sum(s for _, s in reps) == int(ok.sum())


def basis_from(gv):
    G = blocks_of(gv)
    d = []
    for k in range(8):
        z = 0
        for j in range(8):
            if G[k // 2, j // 2, k % 2, j % 2] and (k // 2 != j // 2):
                z |= 1 << j
        d.append((1 << k) | (z << 8))
    return d


total = 0
nsol_total = 0
for idx, (gv, sz) in enumerate(reps):
    d = basis_from(gv)
    t1 = time.time()
    pr = subprocess.run([INV8, "5", "1000000000000", "4"], input=" ".join("%x" % v for v in d) + "\n",
                        capture_output=True, text=True)
    line = [l for l in pr.stdout.splitlines() if l.startswith("NSOL")][0]
    sols = [l for l in pr.stdout.splitlines() if l.startswith("SOL")]
    ns = int(line.split()[1])
    nsol_total += ns
    P("class", idx, "gv=%06x" % gv, "orbit", sz, "d=" + " ".join("%x" % v for v in d), "::", line, "(%.1fs)" % (time.time() - t1))
    for s in sols:
        P("   ", s)
P("DONE classes", len(reps), "total NSOL", nsol_total)
