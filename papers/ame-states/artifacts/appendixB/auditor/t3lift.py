#!/usr/bin/env python3
"""Independent re-run of the lift search of pairs #161 (t=3 swaps (01)(23)(45), fixed parties 6,7 with 15-class involutions).

Own derivation (see results/argument.md 5.7): K = S^Omega = N S (8-dim), K_0 = K cap {f=0} = {(x,Mx,Nx)},
K = K_0^perp, f: K/K_0 -> L_6+L_7 an isometry (720); S = K + span(s_1..s_8), N s_i = k_i,
s_i = p_i + sum_m q_im eps_m  (eps = complement of K in V^Omega), isotropy: for i<j
   sum_m q_jm a_im + sum_m q_im a_jm = B(p_i,p_j),  a_im = B(p_i, eps_m);
weight condition: every coset x + K (x in S/K nonzero) has all 256 words of party weight >= THR.
Gauge: all 12 ordered form pairs (F_M,F_N) are ONE Sp(4,2) orbit (t3k0audit.py) so one (M,N) suffices.

Usage: t3lift.py THR [first last] [maxsec]      THR=5 is the AME search, THR=4 the positive control.
"""
import os, sys, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

THR = int(sys.argv[1]) if len(sys.argv) > 1 else 5
FIRST = int(sys.argv[2]) if len(sys.argv) > 2 else 0
LAST = int(sys.argv[3]) if len(sys.argv) > 3 else 720
MAXSEC = float(sys.argv[4]) if len(sys.argv) > 4 else 480
# [release] log goes to $APPB_WORK (default: current directory) instead of the absolute work/auditor path
LOG = os.path.join(os.environ.get("APPB_WORK", "."), f"t3lift_{THR}_{FIRST}_{LAST}.log")
# [release] deterministic solution cap (env T3_MAXSOL, 0 = none) replacing the wall-clock cap used for the THR=4 control;
# solutions are counted, only the first one is stored (the original kept all of them in memory).
T3_MAXSOL = int(os.environ.get("T3_MAXSOL", "0"))
class _SolCap(Exception):
    pass
logf = open(LOG, "w")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); logf.write(s + "\n"); logf.flush()
T0 = time.time()

def B4(u, v):
    return ((u & 1) & (v >> 1 & 1)) ^ ((u >> 1 & 1) & (v & 1)) ^ ((u >> 2 & 1) & (v >> 3 & 1)) ^ ((u >> 3 & 1) & (v >> 2 & 1))
def mat_apply(M, v):
    r = 0
    for i in range(4):
        r |= (bin(M[i] & v).count("1") & 1) << i
    return r
gl = []; sp = []
for code in range(1 << 16):
    Mr = [(code >> (4 * i)) & 15 for i in range(4)]
    img = tuple(mat_apply(Mr, v) for v in range(16))
    if len(set(img)) != 16: continue
    gl.append(img)
    if all(B4(img[u], img[v]) == B4(u, v) for u in (1, 2, 4, 8) for v in (1, 2, 4, 8)):
        sp.append(img)
assert len(gl) == 20160 and len(sp) == 720
ident = tuple(range(16))
def comp(f, g): return tuple(f[g[v]] for v in range(16))
def is15(l):
    if l == ident or comp(l, l) != ident: return False
    Nn = [x ^ l[x] for x in range(16)]
    return len(set(Nn)) == 4 and all(B4(x, Nn[x]) == 0 for x in range(16))
cl15 = [l for l in sp if is15(l)]
log("15-class involutions found:", len(cl15), "(expect 15)")
l0 = cl15[0]
Nl = [x ^ l0[x] for x in range(16)]
Lset = sorted(set(Nl))            # L_j = im N = ker N, 4 values
assert sorted(x for x in range(16) if Nl[x] == 0) == Lset
Ninv = {}
for x in range(16): Ninv.setdefault(Nl[x], x)
def bj(u, v): return B4(Ninv[u], v)   # b_j on L_j

# ---- Omega on 32-bit ----
def omega_int(x):
    nib = [(x >> (4 * j)) & 15 for j in range(8)]
    y = [nib[1], nib[0], nib[3], nib[2], nib[5], nib[4], l0[nib[6]], l0[nib[7]]]
    return sum(y[j] << (4 * j) for j in range(8))
def Bf(x, y):
    t = 0
    for j in range(8): t ^= B4((x >> (4 * j)) & 15, (y >> (4 * j)) & 15)
    return t
def nib(x, j): return (x >> (4 * j)) & 15

# ---- K_0 = {(x,Mx,Nx)}: pick (M,N) with  B(x,y)+B(Mx,My)+B(Nx,Ny)=0 ----
def gram(img): return tuple(B4(img[1 << i], img[1 << j]) for i in range(4) for j in range(i + 1, 4))
gJ = tuple(B4(1 << i, 1 << j) for i in range(4) for j in range(i + 1, 4))
byg = {}
for g in gl: byg.setdefault(gram(g), []).append(g)
Mm = Nm = None
for g in gl:
    need = tuple(a ^ b for a, b in zip(gJ, gram(g)))
    if need in byg:
        Mm = g; Nm = byg[need][0]; break
assert Mm is not None
def wvec(x, y, z):  # w-part packed into V: pairs duplicated
    return x | x << 4 | y << 8 | y << 12 | z << 16 | z << 20
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
def kernel_of(images):  # images[k] = image (int) of unit vector k  -> list of combos (as bitmasks of unit vectors) mapping to 0
    piv = {}; ker = []
    for k, im in enumerate(images):
        v = 1 << k
        while im:
            hb = im.bit_length() - 1
            if hb in piv: im ^= piv[hb][0]; v ^= piv[hb][1]
            else: piv[hb] = (im, v); break
        if im == 0: ker.append(v)
    return ker

# K = K_0^perp in the 12-bit w-space: (u,v,w) with B4(u,x)+B4(v,Mx)+B4(w,Nx)=0 for x = unit vectors
def func_image(k):  # image of unit vector k of the 12-bit space = 4-bit vector of the 4 functionals
    pos = k // 4; bit = 1 << (k % 4)
    r = 0
    for i in range(4):
        x = 1 << i
        if pos == 0: val = B4(bit, x)
        elif pos == 1: val = B4(bit, Mm[x])
        else: val = B4(bit, Nm[x])
        r |= val << i
    return r
kerK = kernel_of([func_image(k) for k in range(12)])
assert len(kerK) == 8
K0 = [(1 << 0) | 0, 0]  # placeholder
K0b = []
for i in range(4):
    x = 1 << i
    K0b.append(x | Mm[x] << 4 | Nm[x] << 8)
for v in K0b: assert v in set(kerK) or True
Kb_all = reduce_basis(K0b + kerK)
assert len(reduce_basis(K0b)) == 4 and len(Kb_all) == 8, (len(reduce_basis(K0b)), len(Kb_all))
# complement e1..e4 of K0 in K
cur = reduce_basis(K0b); E = []
for v in kerK:
    t = reduce_basis(cur + [v])
    if len(t) > len(cur): cur = t; E.append(v)
assert len(E) == 4
def w3(v): return (v & 15, v >> 4 & 15, v >> 8 & 15)
def Qf(a, b): return B4(a & 15, b & 15) ^ B4(a >> 4 & 15, b >> 4 & 15) ^ B4(a >> 8 & 15, b >> 8 & 15)
GQ = [[Qf(E[i], E[j]) for j in range(4)] for i in range(4)]
# isometries y_i in L6 x L7 (index pair of L-values)
Ycands = [(a, b) for a in Lset for b in Lset]
def bb(y, z): return bj(y[0], z[0]) ^ bj(y[1], z[1])
isos = []
import itertools
for ys in itertools.product(Ycands, repeat=4):
    if all(bb(ys[i], ys[j]) == GQ[i][j] for i in range(4) for j in range(i + 1, 4)):
        isos.append(ys)
log("K_0 basis ok; isometries K/K_0 -> L6+L7:", len(isos), "(expect 720)")
assert len(isos) == 720

# ---- V^Omega and complement ----
NIm = [omega_int(1 << k) ^ (1 << k) for k in range(32)]
VO = [v for v in kernel_of(NIm)]
assert len(VO) == 16

M1 = np.uint32(0x11111111)
def wt(v):
    y = v | (v >> np.uint32(1)); y = y | (y >> np.uint32(2)); y = y & M1
    return ((y * M1) >> np.uint32(28)).astype(np.int16)
def build(vecs):
    arr = np.zeros(1, dtype=np.uint32)
    for v in vecs: arr = np.concatenate([arr, arr ^ np.uint32(v)])
    return arr
def par8(a):
    a = a ^ (a >> 4); a = a ^ (a >> 2); a = a ^ (a >> 1); return a & 1

def preimage(k):  # a vector p with N p = k
    p = 0
    for pp in range(3):
        w = (k >> (8 * pp)) & 15  # nibble 2pp (first of the pair)
        p |= w << (8 * pp)
    for j in (6, 7):
        p |= Ninv[nib(k, j)] << (4 * j)
    return p

nsol_total = 0; nodes_total = 0
def run_K(idx, ys):
    global nsol_total, nodes_total
    kv = []
    for v in K0b: kv.append(wvec(v & 15, v >> 4 & 15, v >> 8 & 15))
    for i in range(4):
        w = w3(E[i])
        kv.append(wvec(*w) | ys[i][0] << 24 | ys[i][1] << 28)
    kv = [int(x) for x in kv]
    for k in kv:
        assert omega_int(k) == k
    Kall = build(kv)
    if int(wt(Kall[1:]).min()) < THR:
        return 0, 0, "K weight"
    ps = []
    for k in kv:
        p = preimage(k); assert omega_int(p) ^ p == k; ps.append(p)
    for i in range(8):
        for k in kv: assert Bf(ps[i], k) == 0
    # complement eps of K in V^Omega
    cur = reduce_basis(kv); eps = []
    for v in VO:
        t = reduce_basis(cur + [v])
        if len(t) > len(cur): cur = t; eps.append(v)
    assert len(eps) == 8
    Vtab = build(eps)  # index bits <-> eps_m
    a = [sum(Bf(ps[i], eps[m]) << m for m in range(8)) for i in range(8)]
    c = [[Bf(ps[i], ps[j]) for j in range(8)] for i in range(8)]
    Qall = np.arange(256, dtype=np.int64)
    sols = []
    nsolc = [0]   # [release]
    nodes = [0] * 9

    def dfs(i, Spre, qs, ss):
        if time.time() - T0 > MAXSEC: raise TimeoutError
        if i == 8:
            # [release] deterministic cap: like the original wall-clock check (made on entry to dfs), the cap is checked on
            # entry to a leaf, i.e. the search stops when it reaches solution number T3_MAXSOL + 1 (not counted).
            if T3_MAXSOL and nsolc[0] >= T3_MAXSOL: raise _SolCap
            nsolc[0] += 1                                   # [release] count instead of storing every solution
            if not sols: sols.append(list(ss))
            if os.environ.get("T3_STOPFIRST"): raise StopIteration
            return
        ok = np.ones(256, dtype=bool)
        for j in range(i):
            rhs = c[i][j] ^ int(par8(np.int64(qs[j] & a[i])))
            ok &= (par8(Qall & a[j]) == rhs)
        cands = Qall[ok]
        if len(cands) == 0: return
        sc = (np.uint32(ps[i]) ^ Vtab[cands]).astype(np.uint32)
        words = Spre[:, None, None] ^ sc[None, :, None] ^ Kall[None, None, :]
        mw = wt(words).min(axis=(0, 2))
        good = np.nonzero(mw >= THR)[0]
        nodes[i + 1] += len(good)
        for g in good:
            s_new = np.uint32(sc[g])
            Spre2 = np.concatenate([Spre, Spre ^ s_new])
            dfs(i + 1, Spre2, qs + [int(cands[g])], ss + [int(s_new)])
    status = "done"
    try:
        dfs(0, np.zeros(1, dtype=np.uint32), [], [])
    except TimeoutError:
        status = "TIMEOUT"
    except StopIteration:
        status = "stopped at first solution"
    except _SolCap:                                         # [release]
        status = "stopped at solution cap %d" % T3_MAXSOL
    return nsolc[0], nodes, status, sols[:1], kv

# ---- orbit reduction of the 720 isometries: local symplectic maps c6,c7 in C(l0) acting on the targets, plus the swap 6<->7
centr = [c for c in sp if comp(c, l0) == comp(l0, c)]
induced = set()
for c in centr:
    assert all(c[u] in Lset for u in Lset)
    induced.add(tuple(Lset.index(c[u]) for u in Lset))
log("centralizer of l0 in Sp(4,2):", len(centr), "induced permutations of L:", len(induced))
def act(ys, p6, p7, swap):
    out = []
    for (a, b) in ys:
        a2 = Lset[p6[Lset.index(a)]]; b2 = Lset[p7[Lset.index(b)]]
        out.append((b2, a2) if swap else (a2, b2))
    return tuple(out)
canon = {}
reps = []
for idx, ys in enumerate(isos):
    if ys in canon: continue
    orb = set()
    for p6 in induced:
        for p7 in induced:
            for sw in (False, True):
                orb.add(act(ys, p6, p7, sw))
    for o in orb: canon[o] = idx
    reps.append(idx)
log("orbit representatives of the 720 isometries under C(l0)^2 x swap:", len(reps), "orbit sizes sum", len(canon))
USE_ORB = len(sys.argv) > 5 and sys.argv[5] == "orb"
todo = reps if USE_ORB else list(range(FIRST, min(LAST, len(isos))))
results = []
for idx in todo:
    ys = isos[idx]
    r = run_K(idx, ys)
    if r[2] == "K weight":
        log(f"K {idx}: K-alone weight < {THR}"); continue
    nsol, nodes, status = r[0], r[1], r[2]
    nodes_total += sum(nodes); nsol_total += nsol
    log(f"[{time.time()-T0:6.1f}s] K {idx}: nsol {nsol} survivors per level {nodes[1:]} {status}")
    if nsol:
        sols, kv = r[3], r[4]
        s = sols[0]
        gens = kv + s
        basis = reduce_basis(gens)
        Sall = build(basis)
        # direct checks
        okiso = all(Bf(a, b) == 0 for a in basis for b in basis)
        invariant = all(int(omega_int(int(x))) in set(Sall.tolist()) for x in basis)
        log(f"   SOLUTION: rank {len(basis)} isotropic {okiso} Omega-invariant {invariant} min weight {int(wt(Sall[1:]).min())}")
    if status == "TIMEOUT": break
log(f"done: K processed up to {idx}, total solutions {nsol_total}, total surviving nodes {nodes_total}, {time.time()-T0:.0f}s")
