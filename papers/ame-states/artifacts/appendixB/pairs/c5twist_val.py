"""c5twist_val.py: validation of c5twist.py's weight formula and construction on random candidates (isotropy, rank 16,
direct 2^16-word min party weight == formula min, ame.bad_cuts). Usage: python c5twist_val.py NTRIALS SEED"""
import os, sys
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))  # [release] ame.py lives in ../lib
import numpy as np
import c5twist as C
import ame

ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
hw = C.Hc(C.A4, C.A4)
iso = [i for i in range(1, 65536) if hw[i] == 0]
PSI = np.array([(u, v) for u in range(16) for v in range(16)], dtype=np.uint8)
NORM = C.MT[PSI[:, 0], C.CONJ[PSI[:, 0]]] ^ C.MT[PSI[:, 1], C.CONJ[PSI[:, 1]]]
done = 0
while done < ntr:
    iw = int(rng.choice(iso))
    w = C.A4[iw]
    hv = C.Hc(C.A4, np.broadcast_to(w, C.A4.shape))
    Aidx = np.nonzero(hv == 0)[0]
    perm = rng.permutation(len(Aidx))
    def inspan(v, vs):
        # v in F_16-span of vs (list of arrays)?
        import itertools
        if len(vs) == 1:
            return any(all(int(C.MT[t, int(z)]) == int(y) for z, y in zip(vs[0], v)) for t in range(16))
        for t in range(16):
            tv = C.MT[t, vs[-1]]
            if inspan(np.array([int(p) ^ int(q) for p, q in zip(v, tv)], dtype=np.uint8), vs[:-1]):
                return True
        return False
    a1 = next(C.A4[Aidx[i]] for i in perm if not inspan(C.A4[Aidx[i]], [w]))
    a2 = next(C.A4[Aidx[i]] for i in perm if not inspan(C.A4[Aidx[i]], [w]) and not inspan(C.A4[Aidx[i]], [w, a1]) and not inspan(C.A4[Aidx[i]], [a1]) and not inspan(C.A4[Aidx[i]], [w, a1]))
    # make sure a2 not in span(w,a1)
    G11 = int(C.Hc(a1, a1)); G22 = int(C.Hc(a2, a2)); G12 = int(C.Hc(a1, a2))
    P1 = np.nonzero(NORM == G11)[0]
    P2 = np.nonzero(NORM == G22)[0]
    pr = C.MT[PSI[P1][:, None, 0], C.CONJ[PSI[P2][None, :, 0]]] ^ C.MT[PSI[P1][:, None, 1], C.CONJ[PSI[P2][None, :, 1]]]
    ii, jj = np.nonzero(pr == G12)
    if len(ii) == 0:
        continue
    k = int(rng.integers(len(ii)))
    p1 = PSI[P1[ii[k]]]; p2 = PSI[P2[jj[k]]]
    # formula min over all (t,c) != 0
    T = np.array([(t0, t1, t2) for t0 in range(16) for t1 in range(16) for t2 in range(16)], dtype=np.uint8)
    av = C.MT[T[:, 0][:, None], w[None, :]] ^ C.MT[T[:, 1][:, None], a1[None, :]] ^ C.MT[T[:, 2][:, None], a2[None, :]]
    ia = C.aindex(av)
    u5 = C.MT[T[:, 1], p1[0]] ^ C.MT[T[:, 2], p2[0]]
    u6 = C.MT[T[:, 1], p1[1]] ^ C.MT[T[:, 2], p2[1]]
    nzu = (u5 != 0).astype(np.int16) + (u6 != 0)
    wf = C.WTC[ia].astype(np.int16) + nzu[:, None] + (np.arange(16)[None, :] != 0)
    nz = (T.sum(axis=1) != 0)
    mn = min(int(wf[nz].min()), int(wf[~nz][:, 1:].min()))
    rows = C.build_S_rows([int(z) for z in w], [int(z) for z in a1], [int(z) for z in a2], [int(z) for z in p1], [int(z) for z in p2])
    isoq, rk = C.gf2_rank_isotropic(rows)
    R = np.array(rows, dtype=np.uint8)
    allw = np.zeros((1, 32), dtype=np.uint8)
    for r in R:
        allw = np.concatenate([allw, allw ^ r[None, :]], axis=0)
    X, Z = allw[:, :16], allw[:, 16:]
    party = ((X[:, 0::2] | X[:, 1::2] | Z[:, 0::2] | Z[:, 1::2]) != 0).sum(axis=1)
    pm = int(party[1:].min())
    print("trial", done, "isotropic", isoq, "rank", rk, "direct min party wt", pm, "formula", mn, "MATCH" if pm == mn else "MISMATCH", "bad_cuts", ame.bad_cuts(rows, 2, 2, 8), flush=True)
    done += 1
