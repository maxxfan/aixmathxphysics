#!/usr/bin/env python3
"""Independent audit of pairs #175 (C): stabilizer AME(8,4) with an automorphism
Omega = (5-cycle on parties 0-4) x (order-5 twists h5,h6 on parties 5,6) x (identity on party 7).

Claim to test: there is NO Omega-invariant Lagrangian S_nf of V_nf = ker Phi_5(Omega) (F_2-dim 24,
F_16-dim 6) all of whose nonzero words have party weight >= 5.

Method (shares no code with work/pairs/c5twist.py and does not use the 'graph over w-perp' reduction):
  1. build Omega on F_2^32 (party j = bits 4j..4j+3, B_4(u,v)=u0v1+u1v0+u2v3+u3v2), V_nf = Im(Omega+1);
  2. enumerate all 2^24 vectors of V_nf with numpy, compute party weights;
  3. an F_16-line through x is {p(Omega)x}: keep the lines that are totally isotropic and all of whose
     15 nonzero vectors have weight >= T (T = 5 is the claim, T = 4,3,2 are controls);
  4. a good Lagrangian is spanned by three mutually orthogonal good lines; search for them.

Usage: c5audit.py [T] [h5idx] [h6idx]      (T = weight threshold, default 5)
"""
import os, sys, time, itertools
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

T = int(sys.argv[1]) if len(sys.argv) > 1 else 5
H5IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
H6IDX = int(sys.argv[3]) if len(sys.argv) > 3 else 0

# ---------- Sp(4,2) ----------
def mat_apply(M, v):  # M list of 4 row masks, v 4-bit int
    r = 0
    for i in range(4):
        r |= (bin(M[i] & v).count("1") & 1) << i
    return r

def B4(u, v):
    return ((u & 1) & (v >> 1 & 1)) ^ ((u >> 1 & 1) & (v & 1)) ^ ((u >> 2 & 1) & (v >> 3 & 1)) ^ ((u >> 3 & 1) & (v >> 2 & 1))

sp = []
for code in range(1 << 16):
    M = [(code >> (4 * i)) & 15 for i in range(4)]
    img = [mat_apply(M, v) for v in range(16)]
    if len(set(img)) != 16:
        continue
    ok = True
    for u in (1, 2, 4, 8):
        for v in (1, 2, 4, 8):
            if B4(img[u], img[v]) != B4(u, v):
                ok = False
    if ok:
        sp.append(tuple(img))
assert len(sp) == 720, len(sp)

def compose(f, g):  # f after g
    return tuple(f[g[v]] for v in range(16))
ident = tuple(range(16))
def order(f):
    k, g = 1, f
    while g != ident:
        g = compose(f, g); k += 1
    return k
ord5 = [f for f in sp if order(f) == 5]
assert len(ord5) == 144, len(ord5)
h5 = ord5[H5IDX % 144]
h6 = ord5[H6IDX % 144]

# ---------- Omega on 32-bit ints ----------
def omega_int(x):
    nib = [(x >> (4 * j)) & 15 for j in range(8)]
    y = [0] * 8
    for i in range(5):
        y[(i + 1) % 5] = nib[i]
    y[5] = h5[nib[5]]
    y[6] = h6[nib[6]]
    y[7] = nib[7]
    r = 0
    for j in range(8):
        r |= y[j] << (4 * j)
    return r

def bform(x, y):
    t = 0
    for j in range(8):
        t ^= B4((x >> (4 * j)) & 15, (y >> (4 * j)) & 15)
    return t

# Omega symplectic and of order 5
for _ in range(20):
    x = int(np.random.default_rng(_).integers(0, 1 << 32)); y = int(np.random.default_rng(100 + _).integers(0, 1 << 32))
    assert bform(omega_int(x), omega_int(y)) == bform(x, y)
    z = x
    for _k in range(5):
        z = omega_int(z)
    assert z == x

# ---------- V_nf = Im(Omega + 1), basis by elimination ----------
def reduce_basis(vecs):
    basis = []  # list of (pivotbit, vec)
    for v in vecs:
        for p, b in basis:
            if v >> p & 1:
                v ^= b
        if v:
            p = v.bit_length() - 1
            # keep reduced echelon: eliminate p from others
            basis = [(q, (b ^ v) if (b >> p & 1) else b) for q, b in basis]
            basis.append((p, v))
    return [b for _, b in basis]

imgs = [omega_int(1 << k) ^ (1 << k) for k in range(32)]
bnf = reduce_basis(imgs)
assert len(bnf) == 24, len(bnf)
# V_nf is killed by Phi_5(Omega) = 1+Omega+...+Omega^4
for b in bnf:
    z, s = b, 0
    for _k in range(5):
        s ^= z; z = omega_int(z)
    assert s == 0
t0 = time.time()

def build(vecs):  # all 2^len XOR combinations, index bit k <-> vecs[k]
    arr = np.zeros(1, dtype=np.uint32)
    for v in vecs:
        arr = np.concatenate([arr, arr ^ np.uint32(v)])
    return arr

om = [bnf]
for m in range(1, 4):
    om.append([omega_int(v) for v in om[-1]])
X = [build(om[m]) for m in range(4)]
N = X[0].shape[0]
assert N == 1 << 24

M1 = np.uint32(0x11111111); M5 = np.uint32(0x55555555)

def wt(v):
    y = v | (v >> np.uint32(1))
    y = y | (y >> np.uint32(2))
    y = y & M1
    return ((y * M1) >> np.uint32(28)).astype(np.uint8)

def par32(v):
    v = v ^ (v >> np.uint32(16)); v = v ^ (v >> np.uint32(8)); v = v ^ (v >> np.uint32(4))
    v = v ^ (v >> np.uint32(2)); v = v ^ (v >> np.uint32(1))
    return v & np.uint32(1)

def pc(v):  # popcount mod 2 only needed
    return par32(v)

def swap(v):
    return ((v & M5) << np.uint32(1)) | ((v >> np.uint32(1)) & M5)

def Bar(x, y):
    return pc(x & swap(y))

# isotropy of the F_16-line of x: B(x, Omega^k x) = 0 for k = 1,2
iso = (Bar(X[0], X[1]) == 0) & (Bar(X[0], X[2]) == 0)
nz = X[0] != 0
print(f"[{time.time()-t0:5.1f}s] V_nf built; isotropic nonzero vectors = {int((iso & nz).sum())} "
      f"(expected 15 * 279825 = {15*279825})", flush=True)

# all 15 nonzero F_16-multiples must have weight >= T
good = np.ones(N, dtype=bool)
minw = np.full(N, 99, dtype=np.uint8)
for c in range(1, 16):
    v = np.zeros(N, dtype=np.uint32)
    for m in range(4):
        if c >> m & 1:
            v ^= X[m]
    w = wt(v)
    minw = np.minimum(minw, w)
    good &= (w >= T)
good[0] = False
cand = good & iso & nz
print(f"[{time.time()-t0:5.1f}s] T={T}: good isotropic vectors {int(cand.sum())}  -> lines {int(cand.sum())//15} "
      f"(lines of all 7 weights: min-weight histogram of isotropic lines below)", flush=True)
hist = np.bincount(minw[iso & nz], minlength=9)
print("  min-weight histogram over isotropic nonzero vectors (value = min over 15 multiples):",
      {i: int(hist[i]) // 15 for i in range(len(hist)) if hist[i]}, "lines", flush=True)

idx = np.nonzero(cand)[0]
vecs = X[0][idx].astype(np.uint64)
# one representative per line: keep the numerically smallest of the 15 multiples
mult = []
seen = set()
reps = []
vlist = [int(v) for v in vecs]
vset = set(vlist)
for v in vlist:
    if v in seen:
        continue
    line = set()
    z = [v, omega_int(v), omega_int(omega_int(v)), omega_int(omega_int(omega_int(v)))]
    for c in range(1, 16):
        s = 0
        for m in range(4):
            if c >> m & 1:
                s ^= z[m]
        line.add(s)
    assert line <= vset, "multiples of a good line must be good"
    seen |= line
    reps.append((v, z))
print(f"[{time.time()-t0:5.1f}s] distinct good isotropic lines: {len(reps)}", flush=True)
if not reps:
    print("RESULT: no good isotropic line at all => no good Lagrangian"); sys.exit(0)

# orthogonality graph between lines (generator x_i, basis z_i[0..3] of the line)
L = len(reps)
if L > 40000:
    print("too many lines for the quadratic orthogonality pass; stop"); sys.exit(1)
Zs = np.array([[z for z in r[1]] for r in reps], dtype=np.uint32)  # L x 4
adj = np.ones((L, L), dtype=bool)
for m in range(4):
    # B(x_i, Omega^m x_j) = 0 for all m (use basis of line j: Zs[:, m])
    for i in range(L):
        pass
xi = Zs[:, 0]
for m in range(4):
    yj = Zs[:, m]
    for i in range(L):
        adj[i] &= (Bar(np.full(L, xi[i], dtype=np.uint32), yj) == 0)
print(f"[{time.time()-t0:5.1f}s] orthogonality computed, mean degree {adj.sum(1).mean():.1f}", flush=True)

def span_good(gens, thr):
    """gens: list of F_2 vectors (ints) spanning F_2-space; True iff all nonzero words weight >= thr and dim = len."""
    arr = build(gens)
    if len(set(arr.tolist())) != len(arr):
        return False
    return int(wt(arr[1:]).min()) >= thr

def line_basis(i):
    return [int(z) for z in Zs[i]]

best_dim = 1
found = []
cnt2 = 0
for i in range(L):
    for j in np.nonzero(adj[i])[0]:
        if j <= i:
            continue
        g2 = line_basis(i) + line_basis(int(j))
        if not span_good(g2, T):
            continue
        cnt2 += 1
        best_dim = max(best_dim, 2)
        for k in np.nonzero(adj[i] & adj[j])[0]:
            if k <= j:
                continue
            g3 = g2 + line_basis(int(k))
            if span_good(g3, T):
                found.append((i, int(j), int(k)))
                best_dim = 3
print(f"[{time.time()-t0:5.1f}s] T={T}: good isotropic 2-dim (F_16) subspaces from line pairs: {cnt2}; "
      f"good Lagrangians (3-dim) found: {len(found)}; best dim {best_dim}", flush=True)
