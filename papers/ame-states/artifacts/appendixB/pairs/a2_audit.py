"""a2_audit.py : independent re-implementation (pairs) of the auditor's a=2 sector (#168):
sigma=(012)(345), fixed parties 6,7 with order-3 twists h having a 2-dim fixed space.
S = L (+) T, T = S cap V_nf  (F_4-Lagrangian, 5-dim over F_4, all words weight >= 5), L = S cap V^Omega (6-dim Lagrangian in 12-dim V^Omega).
Stage 1: X_A lines (T cap {B-part=0}), orbits.  Stage 2: X_B compatible, completions T.  Stage 3: L search for each T.
usage: a2_audit.py [stage_limit]   (prints counts; compare with auditor: 270 lines / 3 classes of 90 / 384 completions per rep pair / Vs=2928)
"""
import os, sys, itertools, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
import numpy as np

XM = 0
ZM = 0
for j in range(8):
    XM |= 3 << (4 * j)
    ZM |= 12 << (4 * j)

def parity32(x):
    x = x.copy() if isinstance(x, np.ndarray) else np.uint32(x)
    x ^= x >> np.uint32(16); x ^= x >> np.uint32(8); x ^= x >> np.uint32(4)
    x ^= x >> np.uint32(2); x ^= x >> np.uint32(1)
    return x & np.uint32(1)

def omega(u, v):
    """symplectic form between mask u (int or array) and mask v (int or array)"""
    u = np.asarray(u, dtype=np.uint32); v = np.asarray(v, dtype=np.uint32)
    a = ((u & np.uint32(XM)) << np.uint32(2)) & (v & np.uint32(ZM))
    b = ((v & np.uint32(XM)) << np.uint32(2)) & (u & np.uint32(ZM))
    return parity32(a ^ b)

def wt(v):
    v = np.asarray(v, dtype=np.uint32)
    w = np.zeros(v.shape, dtype=np.int32)
    for j in range(8):
        w += (((v >> np.uint32(4 * j)) & np.uint32(15)) != 0).astype(np.int32)
    return w

# ---- Sp(4,2) on nibbles
def nib_omega(u, v):
    return (bin((u & 3) & ((v >> 2) & 3)).count("1") + bin(((u >> 2) & 3) & (v & 3)).count("1")) & 1

def transvection(w):
    return tuple((u ^ w) if nib_omega(u, w) else u for u in range(16))

def compose(g, h):  # g after h
    return tuple(g[h[u]] for u in range(16))

IDN = tuple(range(16))
def sp42():
    gens = [transvection(w) for w in range(1, 16)]
    seen = {IDN}; frontier = [IDN]
    while frontier:
        nf = []
        for g in frontier:
            for t in gens:
                x = compose(t, g)
                if x not in seen:
                    seen.add(x); nf.append(x)
        frontier = nf
    return sorted(seen), gens

SP, TRANSV = sp42()
assert len(SP) == 720
def order(g):
    k, x = 1, g
    while x != IDN:
        x = compose(g, x); k += 1
    return k
def nfixed(g):
    return sum(1 for u in range(16) if g[u] == u)

# pick h of order 3 with 4 fixed nibbles (2-dim fixed space)
H = [g for g in SP if order(g) == 3 and nfixed(g) == 4]
print("order-3 elements with 2-dim fixed space:", len(H), " fpf order-3:", sum(1 for g in SP if order(g) == 3 and nfixed(g) == 1))
h = H[0]
CH = [g for g in SP if compose(g, h) == compose(h, g)]
print("|C(h)| =", len(CH))

# ---- Omega on masks
def party_map(parts):
    """parts: dict j -> (dest j', table). returns function on uint32 array"""
    def f(v):
        v = np.asarray(v, dtype=np.uint32)
        out = np.zeros(v.shape, dtype=np.uint32)
        for j, (jd, tab) in parts.items():
            nib = ((v >> np.uint32(4 * j)) & np.uint32(15)).astype(np.int64)
            out |= (np.array(tab, dtype=np.uint32)[nib] << np.uint32(4 * jd))
        return out
    return f

idt = IDN
OMEGA = party_map({0: (1, idt), 1: (2, idt), 2: (0, idt), 3: (4, idt), 4: (5, idt), 5: (3, idt), 6: (6, h), 7: (7, h)})
def Om(v):
    return OMEGA(v)

def gauss_basis(vecs):
    basis = []  # list of (pivot, vec)
    for v in vecs:
        v = int(v)
        for p, b in basis:
            if (v >> p) & 1:
                v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, (b ^ v) if (b >> p) & 1 else b) for q, b in basis]
            basis.append((p, v))
    return [b for _, b in basis]

def span(basis):
    a = np.zeros(1, dtype=np.uint32)
    for b in basis:
        a = np.concatenate([a, a ^ np.uint32(b)])
    return a

unit = np.array([1 << i for i in range(32)], dtype=np.uint32)
Vnf_basis = gauss_basis((unit ^ Om(unit)).tolist())
print("dim V_nf =", len(Vnf_basis))
VNF = span(Vnf_basis)
assert np.all(Om(Om(VNF)) ^ Om(VNF) ^ VNF == 0)
# V^Omega
Vfix_basis = []
for i in range(32):
    pass
# fixed space = {v : Om v = v}: use sum over orbit  v + Om v + Om^2 v  (image of norm map)
Vfix_basis = gauss_basis((unit ^ Om(unit) ^ Om(Om(unit))).tolist())
print("dim V^Omega =", len(Vfix_basis))
VFIX = span(Vfix_basis)
assert np.all(Om(VFIX) == VFIX)

def line_key(v):
    v = int(v); a = int(Om(np.uint32(v))); b = int(Om(np.uint32(a)))
    return min(v, a, b)

def cls_wt(v):
    return int(wt(np.uint32(v)))

# ---- Stage 1: X_A lines
Bmask = 0
for j in (3, 4, 5): Bmask |= 15 << (4 * j)
Amask = 0
for j in (0, 1, 2): Amask |= 15 << (4 * j)
def lines_with_zero(maskzero):
    sel = VNF[((VNF & np.uint32(maskzero)) == 0) & (wt(VNF) >= 5)]
    sel = sel[omega(sel, Om(sel)) == 0]
    return sorted(set(line_key(v) for v in sel.tolist()))
XA = lines_with_zero(Bmask)
XB = lines_with_zero(Amask)
print("X_A lines:", len(XA), " X_B lines:", len(XB), "(auditor: 270)")

def apply_local(v, parties, g):
    v = int(v)
    for j in parties:
        n = (v >> (4 * j)) & 15
        v = (v & ~(15 << (4 * j))) | (g[n] << (4 * j))
    return v

def orbits(items, gens):
    items = set(items); seen = set(); orbs = []
    for it in sorted(items):
        if it in seen: continue
        orb = {it}; fr = [it]
        while fr:
            nf = []
            for x in fr:
                for g in gens:
                    y = line_key(g(x))
                    assert y in items, "group maps outside set"
                    if y not in orb:
                        orb.add(y); nf.append(y)
            fr = nf
        seen |= orb; orbs.append(sorted(orb))
    return orbs

gensA = [(lambda v, t=t: apply_local(v, (0, 1, 2), t)) for t in TRANSV]
gensB = [(lambda v, t=t: apply_local(v, (3, 4, 5), t)) for t in TRANSV]
gensC6 = [(lambda v, t=t: apply_local(v, (6,), t)) for t in CH]
gensC7 = [(lambda v, t=t: apply_local(v, (7,), t)) for t in CH]
orbA = orbits(XA, gensA + gensC6 + gensC7)
print("orbits of X_A under Sp_A x C6 x C7:", [len(o) for o in orbA], "(auditor: 3 classes of 90)")
orbB = orbits(XB, gensB + gensC6 + gensC7)
print("orbits of X_B:", [len(o) for o in orbB])

def herm_orth(u, v):
    """line u, vector v: omega(u,v)=0 and omega(u,Om v)=0"""
    u = np.uint32(u)
    return (omega(u, v) == 0) & (omega(u, Om(v)) == 0)

def compat(xa, xb):
    return bool(herm_orth(xa, np.uint32(xb)))

# ---- Stage 2: for each X_A orbit rep, find compatible X_B classes under Sp_B alone
orbB_SpB = orbits(XB, gensB)
print("X_B orbits under Sp_B alone:", [len(o) for o in orbB_SpB])

def canon(basis):
    return tuple(sorted(gauss_basis(basis)))

def complete_T(xa, xb, thr=5):
    """all F_4-Lagrangians T of V_nf containing lines xa, xb with all words weight >= thr.
    Binary branching (contains line of w0 / excludes its coset class): each T is reached once."""
    base = gauss_basis([xa, int(Om(np.uint32(xa))), xb, int(Om(np.uint32(xb)))])
    cur = span(base)
    pool = VNF
    for b in base:
        pool = pool[omega(np.uint32(b), pool) == 0]
    ok = np.ones(len(pool), dtype=bool)
    for s in cur:
        ok &= (wt(pool ^ s) >= thr)
    pool = pool[ok & ~np.isin(pool, cur)]
    pool = pool[omega(pool, Om(pool)) == 0]
    sols = []
    nodes = [0]
    def rec(basis, pool):
        nodes[0] += 1
        r = (10 - len(basis)) // 2
        if r == 0:
            sols.append(canon(basis)); return
        if len(pool) < (4 ** r - 1):   # need at least (4^r-1) nonzero vectors (3 per line)
            return
        w = int(pool[0])
        ow = int(Om(np.uint32(w)))
        nb = gauss_basis(basis + [w, ow])
        newspan = span(nb)
        # branch 1: include line of w0
        p = pool
        p = p[(omega(np.uint32(w), p) == 0) & (omega(np.uint32(ow), p) == 0)]
        p = p[~np.isin(p, newspan)]
        okp = np.ones(len(p), dtype=bool)
        curw = span(basis)
        newwords = newspan[~np.isin(newspan, curw)]
        for n in newwords:
            okp &= (wt(p ^ n) >= thr)
        # isotropy of pool vectors with their Omega-images already ensured
        if (wt(newwords) >= thr).all():
            rec(nb, p[okp])
        # branch 2: exclude the whole coset class
        rec(basis, pool[~np.isin(pool, newspan)])
    rec(base, pool)
    print("   nodes visited", nodes[0], " initial pool", len(pool), flush=True)
    return set(sols)

def complete_T_old(xa, xb, thr=5):
    start = [xa, int(Om(np.uint32(xa))), xb, int(Om(np.uint32(xb)))]
    sols = set()
    def words(basis):
        return span(basis)
    def rec(basis, depth):
        cur = span(basis)
        if len(basis) == 10:
            sols.add(canon(basis)); return
        # candidates: V_nf vectors orthogonal (Hermitian) to the whole current span
        cand = VNF
        mask = np.ones(len(cand), dtype=bool)
        for b in basis:
            mask &= (omega(np.uint32(b), cand) == 0)
        cand = cand[mask]
        # not in current span: weight test over coset: all w + s, w + Om w + s must have wt >= thr
        # test w itself, Om w ; all three nonzero F_4 multiples of w plus every s in cur
        keep = []
        okw = np.ones(len(cand), dtype=bool)
        for s in cur:
            okw &= (wt(cand ^ s) >= thr) | ((cand ^ s) == 0) & False
        # cand not in cur  (cand ^ 0 weights>=thr for s=0 already excludes zero; elements of cur have wt>=thr but may lie in cur)
        incur = np.isin(cand, cur)
        okw &= ~incur
        cand = cand[okw]
        # isotropy of w with Om w
        cand = cand[omega(cand, Om(cand)) == 0]
        # test Om w cosets too
        okw2 = np.ones(len(cand), dtype=bool)
        ow = Om(cand)
        for s in cur:
            okw2 &= (wt(ow ^ s) >= thr)
        cand = cand[okw2]
        # cross term: w + Om w
        okw3 = np.ones(len(cand), dtype=bool)
        sw = cand ^ Om(cand)
        for s in cur:
            okw3 &= (wt(sw ^ s) >= thr)
        cand = cand[okw3]
        seen_lines = set()
        for w in cand.tolist():
            lk = line_key(w)
            # canonical: only extend by lines whose coset class not seen at this node
            nb = gauss_basis(basis + [w, int(Om(np.uint32(w)))])
            key = tuple(sorted(nb))
            if key in seen_lines: continue
            seen_lines.add(key)
            rec(nb, depth + 1)
    rec(gauss_basis(start), 2)
    return sols

if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else 99
    t0 = time.time()
    reps = []
    for oa in orbA:
        xa = oa[0]
        comp = [orb for orb in orbB_SpB if compat(xa, orb[0])]
        # compatibility must be Sp_B invariant
        for orb in comp:
            assert all(compat(xa, x) for x in orb)
        print("X_A rep", xa, "orbit size", len(oa), ": compatible X_B Sp_B-classes sizes", [len(o) for o in comp])
        for orb in comp:
            reps.append((xa, orb[0]))
    print("rep pairs:", len(reps))
    allT = []
    for (xa, xb) in reps[:lim]:
        Ts = complete_T(xa, xb)
        print("rep pair", xa, xb, ": completions T with all words wt>=5:", len(Ts), "t=%.0fs" % (time.time() - t0), flush=True)
        allT.append(((xa, xb), sorted(Ts)))
    import pickle
    pickle.dump(allT, open(os.path.join(os.environ.get("APPB_WORK", "."), "a2_allT.pkl"), "wb"))  # [release] was work/pairs/a2_allT.pkl
