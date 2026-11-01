"""Sector (A): sigma=(0 1 2)(3 4 5), both order-3 twists h6,h7 FIXED-POINT-FREE.  Usage: z3twist5.py fpf [maxsec] [startclass]
Step 1 (imported from z3twist): all X_A solutions (233280) with their d-keys.
Step 2 (new): orbits of A-solutions under Sp(4,2)_A (SA reps), then the pair classes (O_A, O_B) of compatible orbits under the
   centralizers C(h6) x C(h7) (order 18 x 18, F_4-linear on N6, N7), which act simultaneously on X_A and X_B.
Step 3: for each class rep BFS over F_4-Lagrangian completions T (all 4096 words of weight >= 5), then for each T count f in Sp(4,2)
   with S^Omega = {(a,a,a,f(a),f(a),f(a),0,0)} such that every s+t has weight >= 5.  Any (T,f) hit is an AME(8,4) CANDIDATE
   (stored in work/auditor/z3twist5_hits.txt) -> must be rebuilt as 16-qubit rows and checked with ame.bad_cuts.
Usage: z3twist5.py fpf [maxsec] [startclass] [endclass]
Log: work/auditor/z3twist5_<startclass>.log  (two parallel runs need different startclass values)
"""
import sys
sys.argv = [sys.argv[0], "fpf"] + sys.argv[1:]
MAXSEC = float(sys.argv[2]) if len(sys.argv) > 2 else 800
START = int(sys.argv[3]) if len(sys.argv) > 3 else 0
END = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9     # classes START..END-1 (parallel runs: use two ranges)
from z3twist import *
T_START = time.time()
LOGF = open(os.path.join(os.environ.get("APPB_WORK", "."), "z3twist5_%d.log" % START), "w")  # [release] output to the working dir
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()
log("step 1 done: X_A solutions", len(SA), "keys", len(keys), "comp key pairs", comp_pairs, "total pairs", total_pairs)

def actA(x, g):
    x = int(x); a = x & 0xFFF; na = 0
    for j in range(3):
        na |= g[(a >> (4 * j)) & 15] << (4 * j)
    return (x & ~0xFFF) | na

by_key = collections.defaultdict(list)
for s in SA:
    by_key[tuple(dpart(x) for x in s)].append(tuple(int(x) for x in s))
orbid = {}; orb_rep = []; orb_key = []
for k, sols in sorted(by_key.items()):
    for s in sols:
        if s in orbid: continue
        oid = len(orb_rep); orb_rep.append(s); orb_key.append(k)
        for g in SP:
            orbid[tuple(actA(x, g) for x in s)] = oid
log("A-orbits under Sp(4,2)_A:", len(orb_rep), " covering", len(orbid), "solutions (expect", len(SA), ")")
assert len(orbid) == len(SA)

C6 = [g for g in SP if comp(g, h6) == comp(h6, g)]
C7 = C6
log("centralizer of h6 in Sp(4,2):", len(C6))
def act67(x, g6, g7):
    x = int(x)
    return (x & 0x00FFFFFF) | (g6[(x >> 24) & 15] << 24) | (g7[(x >> 28) & 15] << 28)
def cpart(x): return (int(x) >> 24) & 15
def f4comb(l1, l2, x1, x2):
    v1 = int(mulF4(l1, np.array([x1], dtype=np.uint32))[0]); v2 = int(mulF4(l2, np.array([x2], dtype=np.uint32))[0])
    return v1 ^ v2
def normalize(y1, y2):
    """y1,y2 span an F_4-plane with c-parts a basis of N6; return the basis (z1,z2) with c-parts cb[0], cb[1]"""
    out = []
    for target in cb:
        found = None
        for l1 in range(4):
            for l2 in range(4):
                if l1 == 0 and l2 == 0: continue
                v = f4comb(l1, l2, y1, y2)
                if cpart(v) == target:
                    found = v; break
            if found is not None: break
        assert found is not None
        out.append(found)
    return tuple(out)
# action on orbits
t1 = time.time()
def gen_set(G):
    gens = []; closure = {tuple(ID)}
    for g in G:
        if tuple(g) in closure: continue
        gens.append(g)
        changed = True
        while changed:
            changed = False
            for a in list(closure):
                for b in gens:
                    c = tuple(comp(list(a), b))
                    if c not in closure:
                        closure.add(c); changed = True
    assert len(closure) == len(G)
    return gens
gens6 = gen_set(C6); gens7 = gen_set(C7)
log("generators of C6, C7:", len(gens6), len(gens7))
act = {}
for tag, (g6, g7) in [(("6", i), (g, ID)) for i, g in enumerate(gens6)] + [(("7", i), (ID, g)) for i, g in enumerate(gens7)]:
    mp = []
    for oid, s in enumerate(orb_rep):
        y = normalize(act67(s[0], g6, g7), act67(s[1], g6, g7))
        assert y in orbid, "image leaves solution set"
        mp.append(orbid[y])
    act[tag] = mp
log("action of the generators of C6 x C7 on A-orbits computed: %.0fs" % (time.time() - t1))
# distinct permutations
perms = {tuple(v) for v in act.values()}
log("distinct permutations of the A-orbits:", len(perms))
# compatible orbit pairs
partner_orb = collections.defaultdict(list)
for kA, kBs in partner.items():
    pass
orbs_by_key = collections.defaultdict(list)
for oid, k in enumerate(orb_key): orbs_by_key[k].append(oid)
pairs = []
for kA in keylist:
    for kB in partner[kA]:
        for oa in orbs_by_key[kA]:
            for ob in orbs_by_key[kB]:
                pairs.append((oa, ob))
log("orbit pairs (O_A,O_B) compatible:", len(pairs))
pset = {p: i for i, p in enumerate(pairs)}
parent = list(range(len(pairs)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for pm in perms:
    for i, (oa, ob) in enumerate(pairs):
        j = pset[(pm[oa], pm[ob])]
        ri, rj = find(i), find(j)
        if ri != rj: parent[max(ri, rj)] = min(ri, rj)
classes = collections.defaultdict(list)
for i in range(len(pairs)): classes[find(i)].append(i)
creps = sorted(classes.keys())
log("pair classes under C6 x C7 (x Sp_A x Sp_B): %d  (class sizes %s)" % (len(creps), sorted(len(classes[c]) for c in creps)))

# ---- step 3 ----
def shiftB_(x):
    x = int(x); a = x & 0xFFF
    return (x & ~0xFFF) | (a << 12)
def words_of(basis):
    ws = np.zeros(1, dtype=np.uint32)
    for b in basis: ws = np.concatenate([ws, ws ^ M32(b)])
    return ws
def rref_insert(basis, v):
    for p, b in sorted(basis.items(), reverse=True):
        if (v >> p) & 1: v ^= b
    if v == 0: return None
    p = v.bit_length() - 1
    nb = {}
    for q, b in basis.items():
        if (b >> p) & 1: b ^= v
        nb[q] = b
    nb[p] = v
    return nb
def reduce_arr(basis, arr):
    for p, b in sorted(basis.items(), reverse=True):
        arr = np.where(((arr >> M32(p)) & M32(1)) == 1, arr ^ M32(b), arr)
    return arr
VB = [a << 12 for a in VA]
Vnf = (np.array(VA, dtype=np.uint32)[:, None, None, None] | np.array(VB, dtype=np.uint32)[None, :, None, None]
       | (np.array(N6, dtype=np.uint32) << M32(24))[None, None, :, None]
       | (np.array(N7, dtype=np.uint32) << M32(28))[None, None, None, :]).reshape(-1)
log("|V_nf| =", len(Vnf))

def complete(XA, XB):
    xs = []
    for x in list(XA) + [shiftB_(y) for y in XB]:
        xs.append(int(x)); xs.append(int(OmN(np.array([x], dtype=np.uint32))[0]))
    basis = {}
    for x in xs:
        basis = rref_insert(basis, x)
        assert basis is not None
    W = Vnf
    for x in xs:
        W = W[omg(W, np.full(len(W), x, dtype=np.uint32)) == 0]
    xw = words_of(xs)
    keep = np.ones(len(W), dtype=bool)
    for w in xw: keep &= pw(W ^ w) >= 5
    W = W[keep]
    W = W[omg(W, OmN(W)) == 0]
    states = {tuple(sorted(basis.values())): basis}
    while True:
        new_states = {}
        for key, B in states.items():
            Uw = words_of(list(B.values()))
            Wr = W
            for u in B.values(): Wr = Wr[omg(Wr, np.full(len(Wr), u, dtype=np.uint32)) == 0]
            Wr = np.unique(reduce_arr(B, Wr)); Wr = Wr[Wr != 0]
            if len(Wr) == 0: continue
            ok = np.ones(len(Wr), dtype=bool)
            for w in Uw: ok &= pw(Wr ^ w) >= 5
            for t in Wr[ok]:
                t = int(t)
                nb = rref_insert(B, t)
                if nb is None: continue
                ot = int(OmN(np.array([t], dtype=np.uint32))[0])
                nb = rref_insert(nb, ot)
                if nb is None: continue
                kk = tuple(sorted(nb.values()))
                if kk not in new_states: new_states[kk] = nb
        if not new_states: break
        states = new_states
    return states

HITS = os.path.join(os.environ.get("APPB_WORK", "."), "z3twist5_hits.txt")  # [release] output to the working dir
FULL = (len(Vnf).bit_length() - 1) // 2
tot_T = 0; tot_hit = 0
for ci, c in enumerate(creps):
    if ci < START: continue
    if ci >= END: break
    if time.time() - T_START > MAXSEC:
        log("TIME CAP reached before class", ci, "- resume with startclass", ci); break
    oa, ob = pairs[c]
    t2 = time.time()
    st = complete(orb_rep[oa], orb_rep[ob])
    Ts = [b for b in st.values() if len(b) == FULL]
    nhit = 0
    for B in Ts:
        Tw = words_of(list(B.values()))
        good = np.zeros((16, 16), dtype=bool)
        for a in range(1, 16):
            for b in range(1, 16):
                s = a | (a << 4) | (a << 8) | (b << 12) | (b << 16) | (b << 20)
                good[a, b] = int(pw(Tw ^ M32(s)).min()) >= 5
        for f in SP:
            if all(good[a, f[a]] for a in range(1, 16)):
                nhit += 1
                open(HITS, "a").write("class %d T=%s f=%s\n" % (ci, sorted(B.values()), f))
    tot_T += len(Ts); tot_hit += nhit
    log("[%6.0fs] class %d/%d (%d pairs) keys %s %s: T found %d, (T,f) hits %d, %.0fs" % (
        time.time() - T_START, ci, len(creps), len(classes[c]), orb_key[oa], orb_key[ob], len(Ts), nhit, time.time() - t2))
log("DONE: classes processed up to", ci, " total T", tot_T, " total (T,f) hits", tot_hit)
