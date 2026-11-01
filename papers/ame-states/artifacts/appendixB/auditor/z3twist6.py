"""Controls for the sector-A closure of z3twist5.py (fpf): (i) direct verification of found T and of the S = S^Omega + T construction
(isotropic, Omega-invariant, rank 16, min party weight by 2^16 enumeration == the table prediction); (ii) orbit-reduction check: random
UNREDUCED pairs (X_A, X_B) must give the same T-count as the representative of their class; (iii) near-miss statistics (max number of the 15
nonzero a that pass for any f) and ame.bad_cuts of the best near-miss.   Usage: z3twist6.py   (runs z3twist5 step 1-3 first, ~3 min)"""
import sys, random
sys.argv = [sys.argv[0], "1e9", "0"]  # [release] was "30" (a wall-clock cap in z3twist5's own loop); 1e9 = no cap, so the output is machine-independent
from z3twist5 import *
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))  # [release] lib/ame.py (copy of qlib/ame.py)
import ame
random.seed(6)
L6 = open(os.path.join(os.environ.get("APPB_WORK", "."), "z3twist6.log"), "w")  # [release] output to the working dir
def log6(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); L6.write(s + "\n"); L6.flush()

def to_row(v):
    A = [0] * 16; B = [0] * 16
    for j in range(8):
        nib = (v >> (4 * j)) & 15
        A[2 * j] = nib & 1; B[2 * j] = (nib >> 1) & 1
        A[2 * j + 1] = (nib >> 2) & 1; B[2 * j + 1] = (nib >> 3) & 1
    return A + B
def s_of(a, b): return a | (a << 4) | (a << 8) | (b << 12) | (b << 16) | (b << 20)
def table(Tw):
    good = np.zeros((16, 16), dtype=bool)
    mn = np.zeros((16, 16), dtype=np.int64)
    for a in range(1, 16):
        for b in range(1, 16):
            m = int(pw(Tw ^ M32(s_of(a, b))).min()); mn[a, b] = m; good[a, b] = m >= 5
    return good, mn

class_T = {}
best_overall = (0, None)
# (i) + (iii) on every class representative
for ci, c in enumerate(creps):
    oa, ob = pairs[c]
    st = complete(orb_rep[oa], orb_rep[ob])
    Ts = [b for b in st.values() if len(b) == FULL]
    class_T[ci] = len(Ts)
    for ti, B in enumerate(Ts):
        vecs = list(B.values())
        Tw = words_of(vecs)
        # T direct checks
        assert int(pw(Tw[1:]).min()) >= 5
        assert all(int(omg(np.array([u], dtype=np.uint32), np.array(vecs, dtype=np.uint32)).sum()) == 0 for u in vecs)
        assert all(int(OmN(np.array([u], dtype=np.uint32))[0]) in set(Tw.tolist()) for u in vecs)
        good, mn = table(Tw)
        # best f
        best = (-1, None)
        for f in SP:
            sc = sum(1 for a in range(1, 16) if good[a, f[a]])
            if sc > best[0]: best = (sc, f)
        if best[0] > best_overall[0]: best_overall = (best[0], (ci, ti, vecs, best[1]))
        # control (i): build S for the best f and for a random f, compare min weight by full enumeration with table prediction
        for f in (best[1], random.choice(SP)):
            Sv = vecs + [s_of(a, f[a]) for a in (1, 2, 4, 8)]
            rows = [to_row(v) for v in Sv]
            Sall = words_of(Sv)
            assert len(set(Sall.tolist())) == 65536
            direct = int(pw(Sall[1:]).min())
            pred = min(int(pw(Tw[1:]).min()), min(int(mn[a, f[a]]) for a in range(1, 16)))
            assert direct == pred, (direct, pred)
            assert all(sum(r[i] * s[16 + i] + r[16 + i] * s[i] for i in range(16)) % 2 == 0 for r in rows for s in rows)
log6("(i) all T verified (weight >= 5, isotropic, Omega-invariant); S = S^Omega + T: 65536 distinct words, direct min weight == table prediction, isotropic")
log6("T counts per class:", class_T)
sc, (ci, ti, vecs, f) = best_overall
log6("(iii) best near-miss over all T,f: %d of 15 nonzero a pass (class %d, T #%d)" % (sc, ci, ti))
Sv = vecs + [s_of(a, f[a]) for a in (1, 2, 4, 8)]
rows = [to_row(v) for v in Sv]
bad, cuts = ame.bad_cuts(rows, 2, 2, 8)
log6("   ame.bad_cuts of that near-miss: %s of 35, direct min weight %d" % (bad, int(pw(words_of(Sv)[1:]).min())))

# (ii) orbit reduction control: unreduced random pairs
orbit_class = {}
for ci, c in enumerate(creps):
    for i in classes[c]: orbit_class[pairs[i]] = ci
nchk = 0
for trial in range(21):
    ci_target = trial % 7
    cands = [pairs[i] for i in classes[creps[ci_target]]]
    oa, ob = random.choice(cands)
    # random members of the orbits
    gA = random.choice(SP); gB = random.choice(SP)
    sA = tuple(actA(x, gA) for x in orb_rep[oa]); sB = tuple(actA(x, gB) for x in orb_rep[ob])
    assert orbid[sA] == oa and orbid[sB] == ob
    # also random C6 x C7 image
    g6 = random.choice(C6); g7 = random.choice(C7)
    sA = normalize(act67(sA[0], g6, g7), act67(sA[1], g6, g7)); sB = normalize(act67(sB[0], g6, g7), act67(sB[1], g6, g7))
    st = complete(sA, sB)
    Ts = [b for b in st.values() if len(b) == FULL]
    nh = 0
    for B in Ts:
        Tw = words_of(list(B.values()))
        good, mn = table(Tw)
        for f in SP:
            if all(good[a, f[a]] for a in range(1, 16)): nh += 1
    ok = (len(Ts) == class_T[ci_target]) and nh == 0
    nchk += ok
    log6("(ii) trial %d class %d: unreduced pair gives T=%d (rep %d), hits %d  %s" % (trial, ci_target, len(Ts), class_T[ci_target], nh, "OK" if ok else "MISMATCH"))
log6("(ii) %d of 21 unreduced pairs agree with their class representative" % nchk)
