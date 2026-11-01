#!/usr/bin/env python3
"""Independent audit of pairs #175 (C), second route: enumerate the 'graph' family in MY coordinates and
histogram the minimum party weight of S_nf.

Model: Omega = (5-cycle on parties 0-4) x h5 (party 5) x h6 (party 6) x 1 (party 7); party j = bits 4j..4j+3;
B_4(u,v) = u0v1+u1v0+u2v3+u3v2.  C = cycle part of V_nf (F_2-dim 16, F_16-dim 4), W = W5+W6 (F_2-dim 8).
A Lagrangian L of V_nf = C + W (F_16-dim 3) has L_0 = L cap C of F_16-dim >= 1 (projection to W has dim <= 2).
  * dim L_0 = 1: L_0 = <w> (w an isotropic F_16-line of C, 1105 of them), L injects into C (kernel lives on
    parties 5,6, weight <= 2), image = w-perp (dim 3), L = {a + phi(a)}, phi: w-perp -> W with kernel <w>, i.e.
    phi is given by images u1,u2 in W of an F_16-basis e1,e2 of (w-perp)/<w>; L isotropic iff
    B(e_i,Om^m e_j) = B(u_i,Om^m u_j) for all i,j,m.  (Expected: 300 solutions per w = |U(2,4)|.)
  * dim L_0 = 2: L = L_0 + <(0,u)> contains a weight-2 word (u has both coordinates nonzero) -> bad.
The script solves the isotropy equations by brute force over u1,u2 in W (256 x 256) and takes the minimum
party weight over all 4095 nonzero words of every L, with no use of F_16 arithmetic or of unitary groups.

Usage: c5graph.py [h5idx] [h6idx] [wstart wend]
"""
import os, sys, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

H5IDX = int(sys.argv[1]) if len(sys.argv) > 1 else 0
H6IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
LOG = os.path.join(os.environ.get("APPB_WORK", "."), f"c5graph_{H5IDX}_{H6IDX}.log")  # [release] was an absolute path in work/auditor/
logf = open(LOG, "w")
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); logf.write(s + "\n"); logf.flush()

def mat_apply(M, v):
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
    if all(B4(img[u], img[v]) == B4(u, v) for u in (1, 2, 4, 8) for v in (1, 2, 4, 8)):
        sp.append(tuple(img))
assert len(sp) == 720
ident = tuple(range(16))
def compose(f, g): return tuple(f[g[v]] for v in range(16))
def order(f):
    k, g = 1, f
    while g != ident:
        g = compose(f, g); k += 1
    return k
ord5 = [f for f in sp if order(f) == 5]
assert len(ord5) == 144
h5 = ord5[H5IDX % 144]; h6 = ord5[H6IDX % 144]

def omega_int(x):
    nib = [(x >> (4 * j)) & 15 for j in range(8)]
    y = [0] * 8
    for i in range(5):
        y[(i + 1) % 5] = nib[i]
    y[5] = h5[nib[5]]; y[6] = h6[nib[6]]; y[7] = nib[7]
    return sum(y[j] << (4 * j) for j in range(8))
def bform(x, y):
    t = 0
    for j in range(8):
        t ^= B4((x >> (4 * j)) & 15, (y >> (4 * j)) & 15)
    return t

def reduce_basis(vecs):
    basis = []
    for v in vecs:
        for p, b in basis:
            if v >> p & 1: v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, (b ^ v) if (b >> p & 1) else b) for q, b in basis]
            basis.append((p, v))
    return [b for _, b in basis]

def build(vecs):
    arr = np.zeros(1, dtype=np.uint32)
    for v in vecs:
        arr = np.concatenate([arr, arr ^ np.uint32(v)])
    return arr
M1 = np.uint32(0x11111111); M5 = np.uint32(0x55555555)
def wt_arr(v):
    y = v | (v >> np.uint32(1)); y = y | (y >> np.uint32(2)); y = y & M1
    return ((y * M1) >> np.uint32(28)).astype(np.int16)
def par32(v):
    v = v ^ (v >> np.uint32(16)); v = v ^ (v >> np.uint32(8)); v = v ^ (v >> np.uint32(4))
    v = v ^ (v >> np.uint32(2)); v = v ^ (v >> np.uint32(1))
    return v & np.uint32(1)
def swap(v): return ((v & M5) << np.uint32(1)) | ((v >> np.uint32(1)) & M5)
def Bar(x, y): return par32(x & swap(y))
def cyc_wt(x):  # number of nonzero nibbles among parties 0..4 (python int)
    return sum(1 for j in range(5) if (x >> (4 * j)) & 15)

def om_pow(x, m):
    for _ in range(m): x = omega_int(x)
    return x
def f16mult(x, c):  # c in 0..15 acts as sum c_m Omega^m
    s = 0
    for m in range(4):
        if c >> m & 1: s ^= om_pow(x, m)
    return s

# cycle space C and W
Cb = reduce_basis([omega_int(1 << k) ^ (1 << k) for k in range(20)])
assert len(Cb) == 16
Wb = reduce_basis([omega_int(1 << k) ^ (1 << k) for k in range(20, 28)])
assert len(Wb) == 8
log(f"h5={h5} h6={h6}; dim C={len(Cb)} dim W={len(Wb)}")

# isotropic F16-lines of C
Cv = build(Cb)
Cim = [build([om_pow(b, m) for b in Cb]) for m in range(4)]
iso = (Bar(Cv, Cim[1]) == 0) & (Bar(Cv, Cim[2]) == 0) & (Cv != 0)
isovec = Cv[iso]
log("isotropic nonzero vectors in C:", len(isovec), "-> lines", len(isovec) // 15, "(expected 1105)")
seen = set(); wlines = []
for x in isovec.tolist():
    if x in seen: continue
    mult = {f16mult(x, c) for c in range(1, 16)}
    assert len(mult) == 15
    seen |= mult
    wlines.append(x)
assert len(wlines) == 1105, len(wlines)

# W-side tables
Wv = build(Wb)  # 256 vectors, party5 / party6 nibbles
Wim = [build([om_pow(b, m) for b in Wb]) for m in range(4)]
Mu = np.zeros((16, 256), dtype=np.uint32)  # Mu[c][i] = f16mult(Wv[i], c)
for c in range(16):
    acc = np.zeros(256, dtype=np.uint32)
    for m in range(4):
        if c >> m & 1: acc ^= Wim[m]
    Mu[c] = acc
Wi = [np.array(Wim[m]) for m in range(4)]

hist = {}
bad_dim_count = 0
nsol_hist = {}
total_L = 0
t0 = time.time()
best = []
for wi, w in enumerate(wlines):
    wm = [om_pow(w, m) for m in range(4)]
    # w-perp in C: kernel of x -> (B(x, Om^m w))_m on Cb
    piv = {}; ker = []
    for b in Cb:
        f = sum(bform(b, wm[m]) << m for m in range(4)); v = b
        while f:
            hb = f.bit_length() - 1
            if hb in piv:
                f ^= piv[hb][0]; v ^= piv[hb][1]
            else:
                piv[hb] = (f, v); break
        if f == 0:
            ker.append(v)
    assert len(ker) == 12, len(ker)
    Wspan = {f16mult(w, c) for c in range(16)}
    def span_with(base, e):
        return {a ^ f16mult(e, c) for a in base for c in range(16)}
    e1 = next(v for v in ker if v not in Wspan)
    S1 = span_with(Wspan, e1)
    e2 = next(v for v in ker if v not in S1)
    for e in (e1, e2):
        for m in range(4): assert bform(e, wm[m]) == 0
    # cycle-weight table [aw][c1][c2]
    wcT = np.zeros((16, 16, 16), dtype=np.int16)
    for aw in range(16):
        a = f16mult(w, aw)
        for c1 in range(16):
            x1 = a ^ f16mult(e1, c1)
            for c2 in range(16):
                wcT[aw, c1, c2] = cyc_wt(x1 ^ f16mult(e2, c2))
    # isometry targets
    t11 = [bform(e1, om_pow(e1, m)) for m in range(4)]
    t22 = [bform(e2, om_pow(e2, m)) for m in range(4)]
    t12 = [bform(e1, om_pow(e2, m)) for m in range(4)]
    ok1 = np.ones(256, dtype=bool); ok2 = np.ones(256, dtype=bool)
    for m in range(1, 4):
        ok1 &= (Bar(Wv, Wi[m]) == t11[m]); ok2 &= (Bar(Wv, Wi[m]) == t22[m])
    i1 = np.nonzero(ok1)[0]; i2 = np.nonzero(ok2)[0]
    okp = np.ones((len(i1), len(i2)), dtype=bool)
    for m in range(4):
        okp &= (Bar(Wv[i1][:, None], Wi[m][i2][None, :]) == t12[m])
    a1, a2 = np.nonzero(okp)
    U1 = i1[a1]; U2 = i2[a2]
    nsol = len(U1)
    nsol_hist[nsol] = nsol_hist.get(nsol, 0) + 1
    total_L += nsol
    if nsol == 0: continue
    # words: (aw, c1, c2): cycle weight + party5/6 flags of c1*u1 + c2*u2
    Phi = Mu[:, U1][:, None, :] ^ Mu[:, U2][None, :, :]   # (16,16,S)
    f = ((Phi >> np.uint32(20)) & np.uint32(15) != 0).astype(np.int16) + (((Phi >> np.uint32(24)) & np.uint32(15)) != 0).astype(np.int16)
    # total weight[aw,c1,c2,s] = wcT[aw,c1,c2] + f[c1,c2,s]; exclude (0,0,0)
    tot = wcT[:, :, :, None] + f[None, :, :, :]
    tot[0, 0, 0, :] = 99
    mw = tot.reshape(-1, nsol).min(axis=0)
    for v in mw.tolist():
        hist[v] = hist.get(v, 0) + 1
    if wi % 100 == 0:
        log(f"[{time.time()-t0:6.1f}s] w {wi}/1105 nsol {nsol} cumulative L {total_L} hist {dict(sorted(hist.items()))}")
log("solutions-per-w histogram (expect {300: 1105}):", nsol_hist)
log("total graph-type Lagrangians:", total_L, "(expected 1105*300 = 331500)")
log("min-party-weight histogram over all graph-type L (pairs: {2: 4500, 3: 54000, 4: 273000}):", dict(sorted(hist.items())))
log("max over L of min weight:", max(hist))
