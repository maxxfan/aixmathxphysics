"""z3a_val.py: validation of z3a_full.py: (1) total number of isotropic lines in T^perp/T without weight filter = 27,
(2) build the 16-qubit rows of S = S^Omega(g) + L for random (T_B, line, g), check isotropic/rank 16/ direct min party weight == formula.
Usage: python z3a_val.py REP NTB SEED"""
import os, sys, json, itertools
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import z3a_full as Z
from z3a_full import *

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))  # [release] lib/ame.py (copy of qlib/ame.py)
import ame

reps = json.load(open(os.path.join(os.environ.get("APPB_WORK", "."), "z3a_reps.json")))  # [release] written by z3a_stage1.py orbits in the working dir
rep = reps[int(sys.argv[1])]
ntb = int(sys.argv[2])
rng = np.random.default_rng(int(sys.argv[3]))

# party basis for F_4^2 with Tr<b_k,b_l> = J
V16 = [(u0, u1) for u0 in range(4) for u1 in range(4)]


def trin(u, v):
    return int(Z.TR[Z.inner2(np.array(u, dtype=np.uint8), np.array(v, dtype=np.uint8))])


Jt = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
PB = None
for b in itertools.product(V16, repeat=4):
    if all(trin(b[k], b[l]) == Jt[k][l] for k in range(4) for l in range(4)):
        PB = b
        break
assert PB is not None
COORD = {}
for t in range(16):
    u = np.zeros(2, dtype=np.uint8)
    for k in range(4):
        if (t >> k) & 1:
            u = u ^ np.array(PB[k], dtype=np.uint8)
    COORD[tuple(int(z) for z in u)] = t
assert len(COORD) == 16


def locs_row(locs):
    X = [0] * 16
    Zr = [0] * 16
    for j, t in enumerate(locs):
        t1, t2, t3, t4 = [(t >> k) & 1 for k in range(4)]
        X[2 * j], Zr[2 * j], X[2 * j + 1], Zr[2 * j + 1] = t1, t2, t3, t4
    return X + Zr


def rows_of(Lbasis, g):
    rows = []
    for v in Lbasis:
        for s in (1, 2):
            w = MT[s, np.asarray(v)]
            a = Z.aidx(w[0:4]); b = Z.aidx(w[4:8])
            locs = [int(Y[a, i]) for i in range(3)] + [int(Y[b, i]) for i in range(3)] + [COORD[tuple(int(z) for z in w[8:10])], COORD[tuple(int(z) for z in w[10:12])]]
            rows.append(locs_row(locs))
    for k in range(4):
        c = 1 << k
        d = int(GC[g, c])
        rows.append(locs_row([c] * 3 + [d] * 3 + [0, 0]))
    return rows


def check_rows(rows):
    M = np.array(rows, dtype=np.uint8)
    X, Zm = M[:, :16], M[:, 16:]
    iso = not ((X @ Zm.T + Zm @ X.T) % 2).any()
    A = [int("".join(map(str, r)), 2) for r in rows]
    basis = []
    for v in A:
        for bb in basis:
            v = min(v, v ^ bb)
        if v:
            basis.append(v)
    allw = np.zeros((1, 32), dtype=np.uint8)
    for r in M:
        allw = np.concatenate([allw, allw ^ r[None, :]], axis=0)
    X, Zm = allw[:, :16], allw[:, 16:]
    party = ((X[:, 0::2] | X[:, 1::2] | Zm[:, 0::2] | Zm[:, 1::2]) != 0).sum(axis=1)
    return iso, len(basis), int(party[1:].min())


M = tuple(rep['M'])
al, be = A4[rep['alpha']], A4[rep['beta']]
Me1, Me2 = cols(M)
zero4 = np.zeros(4, dtype=np.uint8)
tA = [np.concatenate([al, zero4, [1, 0], Me1]).astype(np.uint8), np.concatenate([be, zero4, [0, 1], Me2]).astype(np.uint8)]
Mp = minv(mdagger(M))
Mp1, Mp2 = cols(Mp)
TBs = enumA(Mp)
for trial in range(ntb):
    ig, id_ = TBs[int(rng.integers(len(TBs)))]
    ga, de = A4[ig], A4[id_]
    tB = [np.concatenate([zero4, ga, [1, 0], Mp1]).astype(np.uint8), np.concatenate([zero4, de, [0, 1], Mp2]).astype(np.uint8)]
    T = tA + tB
    rows, piv = rref(T)
    P = perp_of(T)
    V = P.copy()
    for r_, p_ in zip(rows, piv):
        V ^= MT[V[:, p_][:, None], r_[None, :]]
    key = np.zeros(len(V), dtype=np.int64)
    for k in range(12):
        key = key * 4 + V[:, k].astype(np.int64)
    uk, inv = np.unique(key, return_inverse=True)
    first = np.zeros(len(uk), dtype=np.int64)
    first[inv[::-1]] = np.arange(len(V))[::-1]
    R = V[first]
    iso = (H12(R, R) == 0)
    nzm = (uk != 0)
    ii = np.nonzero(iso & nzm)[0]
    O = (H12(R[ii][:, None, :], R[ii][None, :, :]) == 0)
    I, J = np.nonzero(np.triu(O, 1))
    allpos = []
    for a_, b_ in zip(I, J):
        pos = []
        for (x, y) in ALL15:
            cv = MT[x, R[ii[a_]]] ^ MT[y, R[ii[b_]]]
            kk = 0
            for k in range(12):
                kk = kk * 4 + int(cv[k])
            pos.append(int(np.searchsorted(uk, kk)))
        if len(set(pos)) == 15:
            allpos.append(tuple(sorted(pos)))
    lines = sorted(set(allpos))
    print("trial", trial, "isotropic lines (no weight filter):", len(lines), "(expect 27)", flush=True)
    # build a random line, random g, check against formula
    ln = lines[int(rng.integers(len(lines)))]
    u1 = R[ln[0]]
    u2 = next(R[q] for q in ln[1:] if H12(u1, R[q]) == 0 and not any((MT[s, u1] == R[q]).all() for s in range(1, 4)))
    basis = T + [u1, u2]
    g = int(rng.integers(720))
    rws = rows_of(basis, g)
    isoq, rk, pm = check_rows(rws)
    words = np.zeros((4096, 12), dtype=np.uint8)
    for k_, b_ in enumerate(basis):
        words ^= MT[COEF6[:, k_][:, None], np.asarray(b_)[None, :]]
    wa = Z.aidx(words[:, 0:4]); wb = Z.aidx(words[:, 4:8])
    base = words[:, 8:10].any(axis=1).astype(np.int16) + words[:, 10:12].any(axis=1).astype(np.int16)
    mn = int((WT12(words)[1:]).min())
    for c in range(1, 16):
        d = int(GC[g, c])
        mn = min(mn, int((WTC[wa, c] + WTC[wb, d] + base).min()))
    print("   isotropic", isoq, "rank", rk, "direct min party weight", pm, "formula", mn, "MATCH" if pm == mn else "MISMATCH", "bad_cuts", ame.bad_cuts(rws, 2, 2, 8), flush=True)
