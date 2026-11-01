"""Exact count N_tot of F_2-linear stabilizer AME(4,4) codes (party = F_2^4 symplectic, Lagrangian C in (F_2^4)^4, min party
weight 3) and size of the automorphism group of one representative in G = Sp(4,2)^4 x| S_4.
If N_tot == |G| / |Aut(C0)| then there is exactly ONE orbit, so a single C0 suffices in the sigma-involution family
(pairs #129 uses 37 classes under the smaller group GL(2,2)^4 x S_4 of graph forms).
C = graph of isometry M: (P0+P1) -> (P2+P3), blocks M = [[A,B],[C,D]] (A,B: ->P2; C,D: ->P3), all four blocks in GL(4,2).
Own derivation, own code.  Usage: python3 ame44_classes.py
"""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import itertools, sys, collections
import numpy as np

# 4x4 matrices over F_2 as numpy uint8
J = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], dtype=np.uint8)

def key(m):
    return int(sum(int(v) << i for i, v in enumerate(m.reshape(-1))))

def mm(a, b):
    return (a.astype(np.int64) @ b.astype(np.int64) % 2).astype(np.uint8)

def is_inv(m):
    # rank over F_2
    rows = [int(sum(int(m[i, j]) << j for j in range(4))) for i in range(4)]
    r = 0
    for bit in range(4):
        piv = None
        for k in range(r, 4):
            if (rows[k] >> bit) & 1:
                piv = k; break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for k in range(4):
            if k != r and (rows[k] >> bit) & 1:
                rows[k] ^= rows[r]
        r += 1
    return r == 4

def inv(m):
    # inverse by brute force via adjugate-free Gauss
    a = np.concatenate([m.astype(np.uint8), np.eye(4, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(4):
        piv = None
        for k in range(r, 4):
            if a[k, c]:
                piv = k; break
        assert piv is not None
        a[[r, piv]] = a[[piv, r]]
        for k in range(4):
            if k != r and a[k, c]:
                a[k] ^= a[r]
        r += 1
    return a[:, 4:].copy()

GL = []
for bits in range(1 << 16):
    m = np.array([(bits >> i) & 1 for i in range(16)], dtype=np.uint8).reshape(4, 4)
    if is_inv(m):
        GL.append(m)
print("|GL(4,2)| =", len(GL))
Q = {}   # A^T J A key -> list of indices
Qa = []
for idx, A in enumerate(GL):
    q = mm(mm(A.T, J), A)
    kq = key(q)
    Qa.append(kq)
    Q.setdefault(kq, []).append(idx)
SP = [A for A, kq in zip(GL, Qa) if kq == key(J)]
print("|Sp(4,2)| =", len(SP))
SPk = {key(g): i for i, g in enumerate(SP)}
GLk = {key(g): i for i, g in enumerate(GL)}
SPa = np.array(SP, dtype=np.int64)          # 720 x 4 x 4

# right coset reps of GL / Sp: A ~ A g
seen = set()
reps = []
for i, A in enumerate(GL):
    if key(A) in seen:
        continue
    reps.append(A)
    for g in SP:
        seen.add(key(mm(A, g)))
print("coset reps:", len(reps))

# symplectic form on 8-dim target E-space = P2 + P3 : omega8(u, v) = u^T (J + J) v blockwise
J8 = np.zeros((8, 8), dtype=np.uint8)
J8[:4, :4] = J; J8[4:, 4:] = J

def omega8(u, v):
    return int((u.astype(np.int64) @ J8.astype(np.int64) @ v.astype(np.int64)) % 2)

all8 = [np.array([(b >> i) & 1 for i in range(8)], dtype=np.uint8) for b in range(256)]

def sym_basis(vecs):
    """find u1,v1,u2,v2 in vecs (list of 8-vectors, a 4-dim symplectic subspace) with omega(u1,v1)=1, omega(u2,v2)=1, rest 0,
    matching the standard form J (pairs (x1,z1),(x2,z2))."""
    nz = [v for v in vecs if v.any()]
    for u1 in nz:
        for v1 in nz:
            if omega8(u1, v1) != 1: continue
            rest = [w for w in nz if omega8(w, u1) == 0 and omega8(w, v1) == 0]
            for u2 in rest:
                for v2 in rest:
                    if omega8(u2, v2) == 1:
                        return [u1, v1, u2, v2]
    raise RuntimeError

def n_of_A(A, want_example=False):
    total = 0
    ex = None
    kA = Qa[GLk[key(A)]]
    need = key(J) ^ kA          # Q(C) must equal J + A^T J A  (all as 16-bit keys; xor of keys == F_2 sum)
    cs = Q.get(need, [])
    for ci in cs:
        C = GL[ci]
        # E = {(A x, C x)}: columns of 8x4 matrix
        Em = np.concatenate([A, C], axis=0)      # 8 x 4
        Evecs = [(Em.astype(np.int64) @ np.array([(b >> i) & 1 for i in range(4)]) % 2).astype(np.uint8) for b in range(16)]
        # perp
        perp = [v for v in all8 if all(omega8(v, e) == 0 for e in Evecs[1:5] + [Evecs[1], Evecs[2], Evecs[4], Evecs[8]])]
        # E is spanned by Evecs[1],[2],[4],[8]
        perp = [v for v in all8 if all(omega8(v, e) == 0 for e in (Evecs[1], Evecs[2], Evecs[4], Evecs[8]))]
        assert len(perp) == 16
        U = np.array(sym_basis(perp), dtype=np.int64).T       # 8 x 4 ; columns u1,v1,u2,v2
        # candidates (B;D) = U g, g in Sp(4,2); need B, D invertible
        # U g : (8x4)(4x4)
        UG = (np.einsum('ij,njk->nik', U, SPa) % 2).astype(np.uint8)   # 720 x 8 x 4
        kb = (UG[:, :4, :].reshape(720, 16).astype(np.int64) * WGT).sum(1)
        kd = (UG[:, 4:, :].reshape(720, 16).astype(np.int64) * WGT).sum(1)
        okm = ISGL[kb] & ISGL[kd]
        total += int(okm.sum())
        if want_example and ex is None and okm.any():
            n = int(np.argmax(okm))
            ex = (A.copy(), UG[n, :4, :].copy(), C.copy(), UG[n, 4:, :].copy())
    return total, ex

WGT = (1 << np.arange(16)).astype(np.int64)
ISGL = np.zeros(1 << 16, dtype=bool)
for A in GL:
    ISGL[key(A)] = True
counts = []
ex_all = None
for r, A in enumerate(reps):
    t, ex = n_of_A(A, want_example=(ex_all is None))
    if ex is not None and ex_all is None:
        ex_all = ex
    counts.append(t)
    print("rep", r, "n(A) =", t, flush=True)
N_tot = 720 * sum(counts)
print("N_tot =", N_tot)

# --- representative C0 and its automorphism group in Sp(4,2)^4 x| S_4
A0, B0, C0, D0 = ex_all
def code_words(A, B, C, D):
    """C = {(x0,x1,x2,x3)}: x2 = A x0 + B x1, x3 = C x0 + D x1 ; return set of 16-bit ints (nibble per party)"""
    W = set()
    for x0 in range(16):
        v0 = np.array([(x0 >> i) & 1 for i in range(4)], dtype=np.int64)
        for x1 in range(16):
            v1 = np.array([(x1 >> i) & 1 for i in range(4)], dtype=np.int64)
            x2 = (A.astype(np.int64) @ v0 + B.astype(np.int64) @ v1) % 2
            x3 = (C.astype(np.int64) @ v0 + D.astype(np.int64) @ v1) % 2
            n2 = int(sum(int(b) << i for i, b in enumerate(x2)))
            n3 = int(sum(int(b) << i for i, b in enumerate(x3)))
            W.add(x0 | (x1 << 4) | (n2 << 8) | (n3 << 12))
    return W

W0 = code_words(A0, B0, C0, D0)
assert len(W0) == 256

def pw(w):
    return sum(1 for i in range(4) if (w >> (4 * i)) & 15)

mw = min(pw(w) for w in W0 if w)
print("C0 min party weight:", mw)
assert mw == 3
# Lagrangian check
def symp_nib(u, v):
    x1, z1, x2, z2 = u & 1, (u >> 1) & 1, (u >> 2) & 1, (u >> 3) & 1
    X1, Z1, X2, Z2 = v & 1, (v >> 1) & 1, (v >> 2) & 1, (v >> 3) & 1
    return (x1 & Z1) ^ (z1 & X1) ^ (x2 & Z2) ^ (z2 & X2)
import random
random.seed(5)
Wl = sorted(W0)
for _ in range(3000):
    a = random.choice(Wl); b = random.choice(Wl)
    s = 0
    for i in range(4):
        s ^= symp_nib((a >> 4 * i) & 15, (b >> 4 * i) & 15)
    assert s == 0
print("C0 isotropic (3000 random pairs), 256 words, min weight 3: ok")

# automorphisms: (pi, g0..g3): w' with w'_{pi(i)} = g_i w_i ; count those with image set == W0
SPlist = SP
def apply_nibble_matrix(g, nib):
    v = np.array([(nib >> i) & 1 for i in range(4)], dtype=np.int64)
    r = (g.astype(np.int64) @ v) % 2
    return int(sum(int(b) << i for i, b in enumerate(r)))

# Precompute action tables: for each g in SP (720) a table nib->nib
tab = []
for g in SP:
    tab.append([apply_nibble_matrix(g, n) for n in range(16)])

def perm_words(W, pi):
    out = set()
    for w in W:
        nw = 0
        for i in range(4):
            nw |= ((w >> 4 * i) & 15) << (4 * pi[i])
        out.add(nw)
    return out

def find_local_maps(W1, W2):
    """all (g0..g3) in Sp(4,2)^4 with g(W1) = W2, W1, W2 AME(4,4) codes"""
    # brute force with propagation: choose g0 (720), g1 (720), g2,g3 determined? Instead search: for each g0,g1 compute image of the
    # sub-code words (x0,x1,*,*): W1 words have x2,x3 functions of (x0,x1); likewise W2.
    # build dicts: for W1: (x0,x1) -> (x2,x3) ; W2 same
    d1 = {}
    for w in W1:
        d1[(w & 15, (w >> 4) & 15)] = ((w >> 8) & 15, (w >> 12) & 15)
    d2 = {}
    for w in W2:
        d2[(w & 15, (w >> 4) & 15)] = ((w >> 8) & 15, (w >> 12) & 15)
    assert len(d1) == 256 and len(d2) == 256
    sols = []
    # g2 determined by g0 via x1=0 words: (x0,0) -> (a(x0), c(x0)) ; need g2 a1(x0) = a2(g0 x0), g3 c1(x0)=c2(g0 x0)
    # enumerate g0 and find g2 as the unique element (if any) of SP with this property; same for g1 with x0 = 0 words.
    def find_g(f1, f2, gsrc):
        """find all g in SP with g(f1(x)) = f2(gsrc x) for all x in 0..15 ; f1,f2 : list nib->nib"""
        res = []
        for gi, t in enumerate(tab):
            if all(t[f1[x]] == f2[gsrc[x]] for x in range(16)):
                res.append(gi)
        return res
    a1 = [d1[(x, 0)][0] for x in range(16)]; c1 = [d1[(x, 0)][1] for x in range(16)]
    b1 = [d1[(0, x)][0] for x in range(16)]; dd1 = [d1[(0, x)][1] for x in range(16)]
    a2 = [d2[(x, 0)][0] for x in range(16)]; c2 = [d2[(x, 0)][1] for x in range(16)]
    b2 = [d2[(0, x)][0] for x in range(16)]; dd2 = [d2[(0, x)][1] for x in range(16)]
    for g0i, t0 in enumerate(tab):
        g2s = find_g(a1, a2, t0)
        if not g2s: continue
        g3s = find_g(c1, c2, t0)
        if not g3s: continue
        for g1i, t1 in enumerate(tab):
            # g2 must also satisfy b-part with g1
            for g2i in g2s:
                if not all(tab[g2i][b1[x]] == b2[t1[x]] for x in range(16)): continue
                for g3i in g3s:
                    if not all(tab[g3i][dd1[x]] == dd2[t1[x]] for x in range(16)): continue
                    sols.append((g0i, g1i, g2i, g3i))
    return sols

# pre-filter g1: to avoid 720*720 python loop cost, it's fine (518400 per perm x checks lazily) -> use early exit
aut_count = 0
per_perm = {}
for pi in itertools.permutations(range(4)):
    Wp = perm_words(W0, pi)
    # Wp must have 256 words and be graph over (x0,x1); (true since all 2-party projections of an AME(4,4) are bijective)
    s = find_local_maps(Wp, W0)
    # verify each found solution completely
    good = 0
    for (g0, g1, g2, g3) in s:
        ok = True
        for w in Wp:
            nw = (tab[g0][w & 15]) | (tab[g1][(w >> 4) & 15] << 4) | (tab[g2][(w >> 8) & 15] << 8) | (tab[g3][(w >> 12) & 15] << 12)
            if nw not in W0:
                ok = False; break
        good += ok
    per_perm[pi] = good
    aut_count += good
    print("perm", pi, "local maps:", good, flush=True)
print("|Aut(C0)| in Sp(4,2)^4 x S_4 =", aut_count)
G = 720 ** 4 * 24
print("|G| =", G, " |G|/|Aut| =", G / aut_count if aut_count else None)
print("N_tot =", N_tot, " single orbit:", G == N_tot * aut_count)
