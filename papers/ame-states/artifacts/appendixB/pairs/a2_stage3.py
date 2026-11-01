"""a2_stage3.py : stage 3 of the independent a=2 audit.  For each T (from a2_allT.pkl) search a 6-dim Lagrangian L of V^Omega
with every nonzero word s of L satisfying wt(s ^ t) >= 5 for all t in T.  usage: a2_stage3.py [first] [last]
Structure used (re-derived): L cap {A-part=0} and L cap {B-part=0} each contain a 2-dim subspace (dim L=6, codim 4), so
enumerate all pairs (Lx0, Ly0) of 2-dim allowed subspaces, ortho to each other, and extend D=Lx0+Ly0 (isotropic, dim 4) by two vectors of D^perp.
"""
import os, sys, pickle, time  # [release] os for the working-dir path
import numpy as np
sys.argv_saved = sys.argv
import a2_audit as A  # noqa  (runs stage 1; harmless)

first = int(sys.argv_saved[1]) if len(sys.argv_saved) > 1 else 0
last = int(sys.argv_saved[2]) if len(sys.argv_saved) > 2 else 10**9
THR = int(sys.argv_saved[3]) if len(sys.argv_saved) > 3 else 5   # THR=4: positive control (L words may have weight 4 against T)

VF = A.VFIX            # 4096 masks, VF[i^j] = VF[i]^VF[j]
N = len(VF)
Amask = A.Amask; Bmask = A.Bmask
OM = np.zeros((N, N), dtype=np.uint8)
for i in range(N):
    OM[i] = A.omega(VF[i], VF).astype(np.uint8)
assert (OM == OM.T).all()
Apart = (VF & np.uint32(Amask)) == 0
Bpart = (VF & np.uint32(Bmask)) == 0
idx = np.arange(N)

def run_T(Tbasis):
    Tw = A.span(list(Tbasis))
    wts = np.zeros((N, len(Tw)), dtype=np.int8)
    allowed = np.ones(N, dtype=bool)
    for k in range(0, N, 512):
        blk = VF[k:k + 512, None] ^ Tw[None, :]
        allowed[k:k + 512] = (A.wt(blk) >= THR).all(axis=1)
    allowed[0] = False
    nallowed = int(allowed.sum())
    X0 = idx[allowed & Apart]    # words with A-part = 0
    Y0 = idx[allowed & Bpart]
    def planes(S):
        Sset = np.zeros(N, dtype=bool); Sset[S] = True
        out = []
        for a in S:
            for b in S:
                if b > a and Sset[a ^ b] and (a ^ b) > b:
                    out.append((int(a), int(b)))
        return out
    PX = planes(X0); PY = planes(Y0)
    nL = 0; npairs = 0; nD = 0
    found = None
    if PX and PY:
        pyA = np.array([p[0] for p in PY]); pyB = np.array([p[1] for p in PY])
        for (a, b) in PX:
            ok = (OM[a][pyA] == 0) & (OM[a][pyB] == 0) & (OM[b][pyA] == 0) & (OM[b][pyB] == 0)
            for j in np.nonzero(ok)[0]:
                c, d = int(pyA[j]), int(pyB[j])
                npairs += 1
                D = [0, a, b, a ^ b]
                D = [x ^ y for x in D for y in (0, c, d, c ^ d)]
                D = np.array(D)
                if not allowed[D[D != 0]].all():   # mixed words Lx0+Ly0 must be allowed too (bug fix 04:21)
                    continue
                nD += 1
                perp = (OM[a] == 0) & (OM[b] == 0) & (OM[c] == 0) & (OM[d] == 0)
                cos = np.ones(N, dtype=bool)
                for x in D:
                    cos &= allowed[idx ^ x] | (idx ^ x == 0) & False
                # coset c0+D allowed: all elements nonzero & allowed -> exclude cosets meeting 0 (those are D itself)
                cos &= ~np.isin(idx, D)
                P = idx[perp & cos]
                # reps modulo D
                reps = {}
                for p in P:
                    key = min(int(p) ^ int(x) for x in D)
                    reps.setdefault(key, int(p))
                R = list(reps.values())
                Pset = set(int(x) for x in P)
                for i1 in range(len(R)):
                    for i2 in range(i1 + 1, len(R)):
                        p, q = R[i1], R[i2]
                        if OM[p][q] == 0 and (p ^ q) in Pset:
                            found = (a, b, c, d, p, q)
                            break
                    if found: break
                if found: break
            if found: break
    if found:   # direct verification: S = L + T, all 2^16 words, min party weight; isotropy
        gens = [int(VF[x]) for x in found] + [int(x) for x in Tbasis]
        words = A.span(gens)
        print("   FOUND: |span|=%d min nonzero party weight=%d" % (len(set(words.tolist())), int(A.wt(words[1:]).min())), flush=True)
    return nallowed, len(PX), len(PY), (npairs, nD), found

if __name__ == "__main__":
    allT = pickle.load(open(os.path.join(os.environ.get("APPB_WORK", "."), "a2_allT.pkl"), "rb"))  # [release] was work/pairs/a2_allT.pkl
    t0 = time.time(); k = 0; nfound = 0
    for (pair, Ts) in allT:
        for T in Ts:
            if k >= first and k < last:
                na, nx, ny, npr, f = run_T(T)
                nfound += 1 if f else 0
                if k % 20 == 0 or f:
                    print("T#%d allowed=%d planesX=%d planesY=%d orthopairs=%d goodD=%d found=%s t=%.0fs" % (k, na, nx, ny, npr[0], npr[1], f, time.time() - t0), flush=True)
            k += 1
    print("DONE Ts processed", min(k, last) - first, " Lagrangians found:", nfound, flush=True)
