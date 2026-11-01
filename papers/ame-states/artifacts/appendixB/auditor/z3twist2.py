"""Step 2 of z3twist: orbit reps of (X_A, X_B) under Sp(4,2)_A x Sp(4,2)_B, then BFS over F_4-Lagrangian completions T with
min party weight >= 5 (T alone).  Usage: z3twist2.py CASE"""
import sys
from z3twist import *      # re-runs step 1 (fast); gives SA, keys, partner, cb, r, N6, N7, VA, mats, Om, OmN, omg, pw, SP ...

t_start = time.time()

def actA(x, g, shift=0):
    x = int(x)
    a = (x >> shift) & 0xFFF
    na = 0
    for j in range(3):
        na |= g[(a >> (4 * j)) & 15] << (4 * j)
    return (x & ~(0xFFF << shift)) | (na << shift)

# group solutions by key, compute orbits of A-solutions under Sp(4,2) acting on a-parts
by_key = collections.defaultdict(list)
for s in SA:
    by_key[tuple(dpart(x) for x in s)].append(tuple(int(x) for x in s))
orbit_reps = {}
for k, sols in by_key.items():
    seen = set()
    reps = []
    sset = set(sols)
    for s in sols:
        if s in seen: continue
        reps.append(s)
        orb = set()
        for g in SP:
            im = tuple(actA(x, g) for x in s)
            assert im in sset, "orbit leaves solution set"
            orb.add(im)
        seen |= orb
    orbit_reps[k] = reps
    print("key", k, ":", len(sols), "solutions,", len(reps), "orbit(s) under Sp(4,2)_A; orbit sizes", [sum(1 for _ in [0]) for _ in reps][:0])
# orbit sizes
for k, sols in by_key.items():
    sset = set(sols); seen = set(); sizes = []
    for s in sols:
        if s in seen: continue
        orb = {tuple(actA(x, g) for x in s) for g in SP}
        sizes.append(len(orb)); seen |= orb
    print("key", k, "orbit sizes", sizes)

def words_of(basis):
    ws = np.zeros(1, dtype=np.uint32)
    for b in basis:
        ws = np.concatenate([ws, ws ^ M32(b)])
    return ws

def rref_insert(basis, v):
    """basis: dict pivot->vec (fully reduced). returns new dict or None if dependent"""
    for p, b in sorted(basis.items(), reverse=True):
        if (v >> p) & 1:
            v ^= b
    if v == 0:
        return None
    p = v.bit_length() - 1
    nb = {}
    for q, b in basis.items():
        if (b >> p) & 1:
            b ^= v
        nb[q] = b
    nb[p] = v
    return nb

def reduce_arr(basis, arr):
    for p, b in sorted(basis.items(), reverse=True):
        arr = np.where(((arr >> M32(p)) & M32(1)) == 1, arr ^ M32(b), arr)
    return arr

VB = [a << 12 for a in VA]
Vnf = (np.array(VA, dtype=np.uint32)[:, None, None, None]
       | np.array(VB, dtype=np.uint32)[None, :, None, None]
       | (np.array(N6, dtype=np.uint32) << M32(24))[None, None, :, None]
       | (np.array(N7, dtype=np.uint32) << M32(28))[None, None, None, :]).reshape(-1)
print("|V_nf| =", len(Vnf))

def complete(XA, XB, verbose=True):
    """BFS over Omega-invariant Lagrangian completions T of X = X_A + X_B in V_nf with all words of weight >= 5.
    returns list of full bases"""
    xs = []
    for x in list(XA) + [shiftB(y) for y in XB]:
        xs.append(int(x)); xs.append(int(OmN(np.array([x], dtype=np.uint32))[0]))
    basis = {}
    for x in xs:
        basis = rref_insert(basis, x)
        assert basis is not None
    # Cd: vectors of V_nf orthogonal to X, weight-5 safe against all X words, isotropic omega(t, Om t) = 0
    W = Vnf
    for x in xs:
        W = W[omg(W, np.full(len(W), x, dtype=np.uint32)) == 0]
    xw = words_of(xs)
    keep = np.ones(len(W), dtype=bool)
    for w in xw:
        keep &= pw(W ^ w) >= 5
    W = W[keep]
    W = W[omg(W, OmN(W)) == 0]
    if verbose: print("   |Cd| =", len(W), " (X words", len(xw), ")")
    states = {tuple(sorted(basis.values())): basis}
    target_f4 = (len(Vnf).bit_length() - 1) // 4 * 0   # unused
    nlevel = 0
    while True:
        new_states = {}
        for key, B in states.items():
            Uw = words_of(list(B.values()))
            Wr = W
            for u in B.values():
                Wr = Wr[omg(Wr, np.full(len(Wr), u, dtype=np.uint32)) == 0]
            Wr = np.unique(reduce_arr(B, Wr))
            Wr = Wr[Wr != 0]
            if len(Wr) == 0: continue
            ok = np.ones(len(Wr), dtype=bool)
            for w in Uw:
                ok &= pw(Wr ^ w) >= 5
            for t in Wr[ok]:
                t = int(t)
                nb = rref_insert(B, t)
                if nb is None: continue
                ot = int(OmN(np.array([t], dtype=np.uint32))[0])
                nb = rref_insert(nb, ot)
                if nb is None: continue
                kk = tuple(sorted(nb.values()))
                if kk not in new_states:
                    new_states[kk] = nb
        if not new_states:
            break
        states = new_states
        nlevel += 1
        if verbose: print("   level", nlevel, ": subspaces", len(states), " F2 dim", len(next(iter(states.values()))), flush=True)
    return states

results = {}
total_T = 0
for kA, reps in orbit_reps.items():
    for kB in partner[kA]:
        # B solutions: shifted copies of A-solutions with key kB ; take orbit reps of that key as well
        for sA in reps:
            for sB in orbit_reps[kB]:
                print("pair keys", kA, kB, flush=True)
                st = complete(sA, sB)
                dims = collections.Counter(len(b) for b in st.values())
                print("   final subspaces by F2-dim:", dict(dims), " time %.0fs" % (time.time() - t_start), flush=True)
                results[(kA, kB, sA, sB)] = st
full_dim = (len(Vnf).bit_length() - 1) // 2
nfull = sum(1 for st in results.values() for b in st.values() if len(b) == full_dim)
print("full-dimension T found (dim %d over F_2):" % full_dim, nfull)
