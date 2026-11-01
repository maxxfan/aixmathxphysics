"""z3a_stage1.py: sigma=(012)(345), fixed parties 6,7 with fpf order-3 twists (auditor's case (A), a=0), F_2-linear stabilizer AME(8,4).
Stage 1: enumerate T_A = S_nf cap {B-part=0} = {(a(f); f, M f) : f in F_4^2}  (F_4-dim 2 inside V_A + V_6 + V_7, Hermitian isotropic,
every nonzero word has cycle weight 3 on A and f6,f7 != 0), count them (auditor: 233280) and compute orbits under
G = Sp(4,2)_A x GU(2,2)_6 x GU(2,2)_7 (order 720*18*18 = 233280).
Model: GF(4)={0,1,2=w,3=w^2}; V_A = F_4^4, y_i(a) = (Tr(a_m w^i))_m, i=0,1,2; Hermitian H_A(a,b) = a1 b2^2 + a2 b1^2 + a3 b4^2 + a4 b3^2
(sum_i B4(y_i(a),y_i(b)) = Tr H_A(a,b)); V_6 = V_7 = F_4^2 with standard Hermitian <f,g> = f1 g1^2 + f2 g2^2.
Usage: python z3a_stage1.py [orbits]
"""
import os, sys, itertools
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

MT = np.zeros((4, 4), dtype=np.uint8)
EXP4 = [1, 2, 3]
LOG4 = {1: 0, 2: 1, 3: 2}
for a in range(1, 4):
    for b in range(1, 4):
        MT[a, b] = EXP4[(LOG4[a] + LOG4[b]) % 3]
CONJ = np.array([0, 1, 3, 2], dtype=np.uint8)
TR = np.array([a ^ int(CONJ[a]) for a in range(4)], dtype=np.uint8)  # in {0,1}
assert set(TR.tolist()) == {0, 1}
W = [1, 2, 3]  # w^0, w^1, w^2

IDX = np.arange(256)
A4 = np.stack([IDX & 3, (IDX >> 2) & 3, (IDX >> 4) & 3, (IDX >> 6) & 3], axis=1).astype(np.uint8)
Y = np.zeros((256, 3), dtype=np.uint8)
for i in range(3):
    v = np.zeros(256, dtype=np.uint8)
    for m in range(4):
        v |= (TR[MT[A4[:, m], W[i]]].astype(np.uint8) << m)
    Y[:, i] = v
WT = (Y != 0).sum(axis=1)


def aidx(a):
    a = np.asarray(a)
    return a[..., 0].astype(np.int64) + 4 * a[..., 1] + 16 * a[..., 2] + 64 * a[..., 3]


def H(a, b):
    return MT[a[..., 0], CONJ[b[..., 1]]] ^ MT[a[..., 1], CONJ[b[..., 0]]] ^ MT[a[..., 2], CONJ[b[..., 3]]] ^ MT[a[..., 3], CONJ[b[..., 2]]]


def smul(t, a):
    return MT[t, a]


# form check
rng = np.random.default_rng(0)
for _ in range(500):
    ia, ib = rng.integers(0, 256, 2)
    s = 0
    for i in range(3):
        u, v = int(Y[ia, i]), int(Y[ib, i])
        u1, u2, u3, u4 = [(u >> k) & 1 for k in range(4)]
        v1, v2, v3, v4 = [(v >> k) & 1 for k in range(4)]
        s ^= (u1 & v2) ^ (u2 & v1) ^ (u3 & v4) ^ (u4 & v3)
    assert s == TR[H(A4[ia], A4[ib])]
    za = smul(2, A4[ia])
    assert [int(Y[aidx(za), i]) for i in range(3)] == [int(Y[ia, (i + 1) % 3]) for i in range(3)]
print("form checks ok", flush=True)

# GL(2,4)
mats = []
for m11, m12, m21, m22 in itertools.product(range(4), repeat=4):
    det = MT[m11, m22] ^ MT[m12, m21]
    if det != 0:
        mats.append((m11, m12, m21, m22))
print("GL(2,4):", len(mats))


def norm(x):
    return int(MT[x, CONJ[x]])


def inner(u, v):  # <u,v> = u1 conj v1 + u2 conj v2
    return int(MT[u[0], CONJ[v[0]]] ^ MT[u[1], CONJ[v[1]]])


# Unitary group GU(2,2): M with <Me_i,Me_j> = delta_ij
def cols(M):
    return ((M[0], M[2]), (M[1], M[3]))  # (Me1, Me2) with M = (m11,m12,m21,m22)


GU = [M for M in mats if inner(cols(M)[0], cols(M)[0]) == 1 and inner(cols(M)[1], cols(M)[1]) == 1 and inner(cols(M)[0], cols(M)[1]) == 0]
print("GU(2,2):", len(GU), "(expect 18)")

# Sp(4,2) acting on a-coordinates as F_2 matrices g (4x4) : a'_m = sum g_mn a_n
J = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
SP = []
for bits in range(1 << 16):
    g = np.array([(bits >> k) & 1 for k in range(16)]).reshape(4, 4)
    if ((g.T @ J @ g) % 2 == J).all():
        SP.append(g.astype(np.uint8))
print("Sp(4,2):", len(SP), "(expect 720)", flush=True)


def apply_g(g, a):  # a: (...,4) uint8 over F_4; g F_2 matrix
    out = np.zeros_like(a)
    for m in range(4):
        for n in range(4):
            if g[m, n]:
                out[..., m] ^= a[..., n]
    return out


# enumerate X_A
f_proj = [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3)]
sols = []
for M in mats:
    Me1, Me2 = cols(M)
    G11 = 1 ^ inner(Me1, Me1)
    G22 = 1 ^ inner(Me2, Me2)
    G12 = inner(Me1, Me2)  # <e1,e2>=0
    c1 = np.nonzero((WT == 3) & (H(A4, A4) == G11))[0]
    c2 = np.nonzero((WT == 3) & (H(A4, A4) == G22))[0]
    # H(alpha,beta) = G12
    pr = H(A4[c1][:, None, :], A4[c2][None, :, :])
    ii, jj = np.nonzero(pr == G12)
    for i, j in zip(ii, jj):
        al, be = A4[c1[i]], A4[c2[j]]
        ok = True
        for (f1, f2) in f_proj[2:]:
            v = smul(f1, al) ^ smul(f2, be)
            if WT[aidx(v)] != 3:
                ok = False
                break
        if ok:
            sols.append((M, int(c1[i]), int(c2[j])))
print("X_A solutions:", len(sols), "(auditor: 233280)", flush=True)
from collections import Counter
print("per-M counts histogram:", sorted(Counter(Counter(s[0] for s in sols).values()).items()))
if len(sys.argv) > 1 and sys.argv[1] == "orbits":
    # orbits under G: (g, U6, U7): M -> U7 M U6^{-1}, a -> g a
    def mmul(P, Q):
        return (int(MT[P[0], Q[0]] ^ MT[P[1], Q[2]]), int(MT[P[0], Q[1]] ^ MT[P[1], Q[3]]), int(MT[P[2], Q[0]] ^ MT[P[3], Q[2]]), int(MT[P[2], Q[1]] ^ MT[P[3], Q[3]]))

    def minv(P):
        det = int(MT[P[0], P[3]] ^ MT[P[1], P[2]])
        di = EXP4[(-LOG4[det]) % 3]
        return (int(MT[di, P[3]]), int(MT[di, P[1]]), int(MT[di, P[2]]), int(MT[di, P[0]]))

    GUinv = {U: minv(U) for U in GU}
    sset = {s: k for k, s in enumerate(sols)}
    parent = list(range(len(sols)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    # generators: a few random elements of Sp, GU (union-find with all single-factor generators suffices if we take full factor groups)
    gens = []
    gsp = [SP[k] for k in rng.choice(len(SP), 6, replace=False)]
    gu = [GU[k] for k in rng.choice(len(GU), 5, replace=False)]
    for si, (M, ia, ib) in enumerate(sols):
        al, be = A4[ia], A4[ib]
        for g in gsp:
            s2 = (M, int(aidx(apply_g(g, al))), int(aidx(apply_g(g, be))))
            parent[find(si)] = find(sset[s2])
        for U in gu:
            # U6 acts on f: new M = M U6^{-1}; with a'(f'') = a(U6^{-1} f'') : alpha' , beta' = columns of U6^{-1}
            Ui = GUinv[U]
            # new a(e1') = a(Ui e1) = Ui11*al + Ui21*be ; a(e2') = Ui12*al + Ui22*be
            al2 = smul(Ui[0], al) ^ smul(Ui[2], be)
            be2 = smul(Ui[1], al) ^ smul(Ui[3], be)
            s2 = (mmul(M, Ui), int(aidx(al2)), int(aidx(be2)))
            parent[find(si)] = find(sset[s2])
            # U7 acts on M only
            s3 = (mmul(U, M), ia, ib)
            parent[find(si)] = find(sset[s3])
    roots = Counter(find(i) for i in range(len(sols)))
    print("orbits under random subgroup generators:", len(roots), sorted(roots.values())[:20])
    import json
    reps = [dict(M=list(sols[r][0]), alpha=sols[r][1], beta=sols[r][2], size=roots[r]) for r in roots]
    json.dump(reps, open(os.path.join(os.environ.get("APPB_WORK", "."), "z3a_reps.json"), "w"))  # [release] was the script dir
    print("reps saved", reps)
