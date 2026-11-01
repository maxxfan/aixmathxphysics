"""AME(8,4) F_2-linear stabilizer with Omega = sigma*h, sigma=(0 1 2)(3 4 5), order-3 local twists h6,h7 on the fixed parties 6,7.
S = S^Omega + T, T = S cap V_nf (F_4-linear, Omega = mult by omega).  Step 1: structure X_A = T cap (A+N6+N7), X_B (dim r over F_4,
r = dim_F4 N6, N6 = ker(1+h+h^2) on party 6): count the possible (X_A, X_B) pairs.  Usage: z3twist.py CASE   (CASE = fpf | nonfpf)
Own code (no sharing with work/pairs).
"""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import sys, itertools, collections, time
import numpy as np

case = sys.argv[1] if len(sys.argv) > 1 else "fpf"

def symp4(u, v):
    x1, z1, x2, z2 = u & 1, (u >> 1) & 1, (u >> 2) & 1, (u >> 3) & 1
    X1, Z1, X2, Z2 = v & 1, (v >> 1) & 1, (v >> 2) & 1, (v >> 3) & 1
    return (x1 & Z1) ^ (z1 & X1) ^ (x2 & Z2) ^ (z2 & X2)

SY = [[symp4(u, v) for v in range(16)] for u in range(16)]
std = [1, 2, 4, 8]
SP = []
for m in itertools.product(range(1, 16), repeat=4):
    if all(SY[m[i]][m[j]] == SY[std[i]][std[j]] for i in range(4) for j in range(i + 1, 4)):
        tab = []
        for nib in range(16):
            r = 0
            for i in range(4):
                if (nib >> i) & 1:
                    r ^= m[i]
            tab.append(r)
        if len(set(tab)) == 16:
            SP.append(tab)
assert len(SP) == 720, len(SP)
ID = list(range(16))
def comp(g, h):  # g after h
    return [g[h[n]] for n in range(16)]
ord3 = []
for g in SP:
    if g != ID and comp(g, comp(g, g)) == ID:
        nfix = sum(1 for n in range(16) if g[n] == n)
        ord3.append((nfix, g))
print("order-3 elements:", collections.Counter(n for n, g in ord3))
want = 1 if case == "fpf" else 4
h6 = next(g for n, g in ord3 if n == want)
h7 = h6   # conjugate classes are equivalent by local conjugation; same rep is WLOG
mats = [ID] * 6 + [h6, h7]
perm = [1, 2, 0, 4, 5, 3, 6, 7]

def Om(v):
    out = 0
    for j in range(8):
        out |= mats[j][(v >> (4 * j)) & 15] << (4 * perm[j])
    return out

# numpy helpers
M32 = np.uint32
try:
    popc = np.bitwise_count
except AttributeError:
    PC = np.array([bin(i).count("1") for i in range(65536)], dtype=np.uint8)
    def popc(x):
        x = x.astype(np.uint32)
        return PC[x & 0xFFFF] + PC[x >> 16]

def omg(u, v):
    u = np.asarray(u, dtype=np.uint32); v = np.asarray(v, dtype=np.uint32)
    sv = ((v & M32(0x55555555)) << M32(1)) | ((v >> M32(1)) & M32(0x55555555))
    return (popc(u & sv) & 1).astype(np.uint8)

def pw(v):
    v = np.asarray(v, dtype=np.uint32)
    y = (v | (v >> M32(1)) | (v >> M32(2)) | (v >> M32(3))) & M32(0x11111111)
    return popc(y).astype(np.int64)

def OmN(v):
    v = np.asarray(v, dtype=np.uint32)
    out = np.zeros_like(v)
    for j in range(8):
        t = np.array(mats[j], dtype=np.uint32)
        out |= t[(v >> M32(4 * j)) & M32(15)] << M32(4 * perm[j])
    return out

# sanity: Omega order 3, symplectic
import random
random.seed(3)
for _ in range(200):
    u = random.getrandbits(32); w = random.getrandbits(32)
    assert Om(Om(Om(u))) == u
    assert int(omg(Om(u), Om(w))) == int(omg(u, w))

# nf pieces
VA = [y0 | (y1 << 4) | ((y0 ^ y1) << 8) for y0 in range(16) for y1 in range(16)]
wt3 = [a for a in VA if ((a & 15) and ((a >> 4) & 15) and ((a >> 8) & 15))]
assert len(wt3) == 210
N6 = [c for c in range(16) if c ^ h6[c] ^ h6[h6[c]] == 0]
N7 = N6
r = {16: 2, 4: 1}[len(N6)]
print("case", case, "|N6| =", len(N6), "r =", r)
def f4line(c):
    return {0, c, h6[c], h6[h6[c]]}
cb = []
span = {0}
for c in N6:
    if c and c not in span:
        cb.append(c)
        span |= f4line(c)
        if len(cb) == r:
            break
assert len(cb) == r
print("c basis:", cb)

# candidate x_i = a | c_i<<24 | d<<28 with a weight-3, d in N7 nonzero; isotropy omega(x, Om x) = 0
cands = []
for i in range(r):
    L = np.array([a | (cb[i] << 24) | (d << 28) for a in wt3 for d in N7 if d], dtype=np.uint32)
    ok = omg(L, OmN(L)) == 0
    cands.append(L[ok])
    print("cand", i, len(L), "->", len(cands[-1]))

def mulF4(k, x):
    """k in 0..3 : 0 ->0, 1 -> x, 2 -> Omega x, 3 -> Omega^2 x = x + Omega x"""
    x = np.asarray(x, dtype=np.uint32)
    if k == 0: return np.zeros_like(x)
    if k == 1: return x
    ox = OmN(x)
    if k == 2: return ox
    return x ^ ox

def XA_solutions():
    if r == 1:
        return [(x,) for x in cands[0]]
    L1, L2 = cands
    sols = []
    O2 = OmN(L2)
    for x1 in L1:
        x1 = M32(x1)
        ox1 = OmN(np.array([x1]))[0]
        c1 = omg(np.full(len(L2), x1), L2) == 0
        c2 = omg(np.full(len(L2), x1), O2) == 0
        c3 = omg(np.full(len(L2), ox1), L2) == 0
        ok = c1 & c2 & c3
        idx = np.nonzero(ok)[0]
        for j in idx:
            x2 = L2[j]
            good = True
            for l1 in range(4):
                for l2 in range(4):
                    if l1 == 0 and l2 == 0: continue
                    v = mulF4(l1, np.array([x1]))[0] ^ mulF4(l2, np.array([x2]))[0]
                    if pw(np.array([v]))[0] < 5:
                        good = False; break
                if not good: break
            if good:
                sols.append((x1, x2))
    return sols

t0 = time.time()
SA = XA_solutions()
print("X_A solutions:", len(SA), "time %.1fs" % (time.time() - t0))
import pickle
pickle.dump({"case": case, "r": r, "cb": cb, "SA": [tuple(int(v) for v in s) for s in SA]}, open(os.path.join(os.environ.get("APPB_WORK", "."), "z3twist_SA_%s.pkl" % case), "wb"))  # [release] output to the working dir
# d-keys
def dpart(x):
    return int(x) >> 28
keys = collections.Counter(tuple(dpart(x) for x in s) for s in SA)
print("distinct d-keys:", len(keys), " solutions per key (min/max):", min(keys.values()) if keys else None, max(keys.values()) if keys else None)

# cross compatibility between A-solution (shifted copy B): omega(xA_i, xB_j) only through c and d parts (a, b in different parties)
def shiftB(x):
    x = int(x)
    a = x & 0xFFF
    return (x & ~0xFFF) | (a << 12)

def basisvecs(s):
    out = []
    for x in s:
        x = M32(x)
        out.append(int(x)); out.append(int(OmN(np.array([x]))[0]))
    return out

keylist = list(keys.keys())
def compat(kA, kB):
    # build only c,d parts of the vectors
    vA = []
    for i in range(r):
        x = (cb[i] << 24) | (kA[i] << 28)
        vA.append(x); vA.append(int(OmN(np.array([x], dtype=np.uint32))[0]))
    vB = []
    for i in range(r):
        x = (cb[i] << 24) | (kB[i] << 28)
        vB.append(x); vB.append(int(OmN(np.array([x], dtype=np.uint32))[0]))
    for u in vA:
        for w in vB:
            if int(omg(u, w)): return False
    return True

# all (kA, kB) compatible pairs among keys
comp_pairs = 0
total_pairs = 0
partner = collections.defaultdict(list)
for kA in keylist:
    for kB in keylist:
        if compat(kA, kB):
            partner[kA].append(kB)
            comp_pairs += 1
            total_pairs += keys[kA] * keys[kB]
print("compatible key pairs:", comp_pairs, " total (X_A, X_B) pairs:", total_pairs)
print("partners per key (distribution):", collections.Counter(len(v) for v in partner.values()))
