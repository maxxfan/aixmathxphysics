"""Independent audit of pairs #138: order-3 party permutation sigma=(0 1 2)(3 4 5), parties 6,7 fixed with trivial twist.
S_nf = {(a, phi a)} with a in W = {(y0,y1,y2): y0+y1+y2=0}, y_i in F_2^4 (2 qubits), phi sigma-equivariant isometry W->W.
Claim: for every such phi, S_nf has a nonzero word of party weight <= 4.
Own code: no sharing with work/pairs.  Usage: python3 z3_audit.py
"""
import sys, collections

def symp4(u, v):
    # qubit1 = bits 0 (x),1 (z); qubit2 = bits 2 (x),3 (z)
    x1, z1, x2, z2 = u & 1, (u >> 1) & 1, (u >> 2) & 1, (u >> 3) & 1
    X1, Z1, X2, Z2 = v & 1, (v >> 1) & 1, (v >> 2) & 1, (v >> 3) & 1
    return (x1 & Z1) ^ (z1 & X1) ^ (x2 & Z2) ^ (z2 & X2)

SY = [[symp4(u, v) for v in range(16)] for u in range(16)]

def par(a):
    return bin(a).count("1") & 1  # [release] was a.bit_count() & 1 (Python >= 3.10 only)

def enc(y0, y1, y2):
    return y0 | (y1 << 4) | (y2 << 8)

def dec(a):
    return a & 15, (a >> 4) & 15, (a >> 8) & 15

def om(a, b):
    a0, a1, a2 = dec(a)
    b0, b1, b2 = dec(b)
    return SY[a0][b0] ^ SY[a1][b1] ^ SY[a2][b2]

def sig(a):
    y0, y1, y2 = dec(a)
    return enc(y2, y0, y1)   # party 0->1->2->0

def wt(a):
    return sum(1 for y in dec(a) if y)

W = [enc(y0, y1, y0 ^ y1) for y0 in range(16) for y1 in range(16)]
Wset = set(W)
assert len(Wset) == 256 and all(sig(a) in Wset for a in W)
# nondegenerate?
rad = [a for a in W if all(om(a, b) == 0 for b in W)]
assert rad == [0], rad
# weights
cnt = collections.Counter(wt(a) for a in W)
print("weights in W:", dict(cnt))
P2 = [a for a in W if wt(a) == 2]
print("weight-2 vectors:", len(P2))

# F_4-basis: e_i with {e_i, sigma e_i} F_2-independent
def span_add(basis, v):
    # basis: list of (pivot,vec) reduced; returns True if independent and adds
    for p, b in basis:
        if v ^ b < v:
            v ^= b
    if v == 0:
        return False
    basis.append((v.bit_length() - 1, v))
    basis.sort(reverse=True)
    return True

def reduce(basis, v):
    for p, b in basis:
        if (v >> p) & 1:
            v ^= b
    return v

basis = []
E = []
for a in W:
    if a == 0:
        continue
    trial = list(basis)
    # use independent test with gaussian elimination
    def indep(bs, v):
        v2 = v
        for p, b in bs:
            if (v2 >> p) & 1:
                v2 ^= b
        return v2
    v = indep(basis, a)
    if v == 0:
        continue
    basis.append((v.bit_length() - 1, v)); basis.sort(reverse=True)
    sv = indep(basis, sig(a))
    if sv == 0:
        # sigma a in span: undo
        basis = [x for x in basis if x[1] != v]
        continue
    basis.append((sv.bit_length() - 1, sv)); basis.sort(reverse=True)
    E.append(a)
    if len(E) == 4:
        break
assert len(E) == 4
B = []
for e in E:
    B.append(e); B.append(sig(e))
# check F_2-independence of the 8
def rank(vs):
    bs = []
    for v in vs:
        for p, b in bs:
            if (v >> p) & 1:
                v ^= b
        if v:
            bs.append((v.bit_length() - 1, v)); bs.sort(reverse=True)
    return len(bs)
assert rank(B) == 8
print("F4 basis E =", E)

# gram constraints
G = [[om(B[k], B[l]) for l in range(8)] for k in range(8)]

# DFS over images f_i of e_i (images of sigma e_i are sigma f_i)
sols = []
def dfs(i, imgs):
    if i == 4:
        sols.append(list(imgs))
        return
    for f in W:
        if f == 0:
            continue
        ok = True
        cand = list(imgs) + [f]
        # new basis vectors indices 2i, 2i+1 : images f, sigma f
        newvecs = [f, sig(f)]
        allv = []
        for j in range(i):
            allv.append(imgs[j]); allv.append(sig(imgs[j]))
        # constraints between new and old, and among new
        for a_idx, nv in enumerate(newvecs):
            k = 2 * i + a_idx
            for l in range(2 * i):
                if om(nv, allv[l]) != G[k][l]:
                    ok = False; break
            if not ok: break
        if not ok: continue
        if om(newvecs[0], newvecs[1]) != G[2 * i][2 * i + 1]:
            continue
        imgs.append(f)
        dfs(i + 1, imgs)
        imgs.pop()

dfs(0, [])
print("number of sigma-equivariant isometries:", len(sols), "(|GU(4,2)| = 77760)")

# build phi on W by F_2-linearity
def build_phi(imgs):
    # images of the 8 basis vectors
    im8 = []
    for f in imgs:
        im8.append(f); im8.append(sig(f))
    # express each a in W in basis B: solve via table of span
    return im8

# table: coordinates of every a in W wrt B
coord = {0: 0}
# generate by combination
span = {0: 0}
for k, b in enumerate(B):
    new = {}
    for v, c in span.items():
        new[v ^ b] = c | (1 << k)
    span.update(new)
assert len(span) == 256

hist = collections.Counter()
hist_w = collections.Counter()
minw_all = 99
for imgs in sols:
    im8 = build_phi(imgs)
    phi = {}
    for v, c in span.items():
        r = 0
        for k in range(8):
            if (c >> k) & 1:
                r ^= im8[k]
        phi[v] = r
    # sanity: isometry on all W, equivariant
    # (spot check a few)
    nlow = 0
    minw = 99
    for a in W:
        if a == 0: continue
        w = wt(a) + wt(phi[a])
        if w < minw: minw = w
        if w <= 4: nlow += 1
    hist[nlow] += 1
    minw_all = min(minw_all, minw)
    hist_w[minw] += 1
print("phi count by #words of party weight <=4:", dict(sorted(hist.items())))
print("min party weight over all phi:", minw_all, " distribution of minw:", dict(hist_w))
# full-isotropy and weight check on one random solution embedded in 8 parties
import random
random.seed(1)
imgs = random.choice(sols)
im8 = build_phi(imgs)
phi = {}
for v, c in span.items():
    r = 0
    for k in range(8):
        if (c >> k) & 1:
            r ^= im8[k]
    phi[v] = r
for _ in range(2000):
    a = random.choice(W); b = random.choice(W)
    assert om(phi[a], phi[b]) == om(a, b)
    assert phi[sig(a)] == sig(phi[a])
    assert phi[a ^ b] == phi[a] ^ phi[b]
print("spot checks (isometry, equivariance, linearity) ok")
