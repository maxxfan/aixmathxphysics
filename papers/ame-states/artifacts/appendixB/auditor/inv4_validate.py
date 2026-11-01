"""Validate inv4.py's reduction: for a random AME(4,4) graph code C and RANDOM symmetric beta, the stabilizer S built by inv4.build_rows must be isotropic of rank 16,
and its direct minimum party weight (all 2^16 words) must equal min over w != 0 of the minimum over x in C of wt8(y_w + x, y_w + x + w) computed with the same y_w.
Also: S is AME  <=>  every (w, beta(w,.)) is 'good'; checked on beta solutions of the RELAXED threshold 4 (so that solutions exist): bad_cuts is reported.
usage: python3 work/auditor/inv4_validate.py SEED"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "lib"))   # [release] was ROOT (qlib package)
import numpy as np
import inv4, c7sym as C7
import ame   # [release] was: from qlib import ame

rng = random.Random(int(sys.argv[1]))
while True:
    Gam = [[0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1, 8):
            Gam[i][j] = Gam[j][i] = rng.randint(0, 1)
    g = inv4.graph_code(Gam)
    Call = inv4.code_elements(g)
    if inv4.party_weights4(Call[1:]).min() >= 3:
        break
h = [inv4.hvec(j) for j in range(8)]
for trial in range(5):
    M = [[0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i, 8):
            M[i][j] = M[j][i] = rng.randint(0, 1)
    beta = [sum(M[i][j] << j for j in range(8)) for i in range(8)]
    rows = inv4.build_rows(g, Call, beta)
    iso, rk = C7.isotropic(rows), C7.rank2(rows)
    R = np.array(rows, dtype=np.int64)
    coeffs = ((np.arange(65536)[:, None] >> np.arange(16)[None, :]) & 1).astype(np.int64)
    allw = (coeffs @ R) % 2
    x = allw[:, :16].reshape(-1, 8, 2); z = allw[:, 16:].reshape(-1, 8, 2)
    pw = ((x | z).sum(axis=2) > 0).sum(axis=1)
    direct_min = int(pw[1:].min())
    # formula: min over w != 0 (coefficient vector wi), x in C
    best = 99
    for wi in range(1, 256):
        w = int(Call[wi])
        lam = 0
        for i in range(8):
            if (wi >> i) & 1:
                lam ^= beta[i]
        y = 0
        for j in range(8):
            if (lam >> j) & 1:
                y ^= h[j]
        Z = y ^ Call
        wt = sum(((Z >> (4 * t)) & 15) != 0 for t in range(4)) + sum((((Z ^ w) >> (4 * t)) & 15) != 0 for t in range(4))
        best = min(best, int(wt.min()))
    # also w = 0 part: 2*wt4(x) >= 6
    bad, cuts = ame.bad_cuts(rows, 2, 2, 8)
    print("trial", trial, "isotropic", iso, "rank", rk, "direct min party weight", direct_min, "formula min (w != 0 part)", best, "bad_cuts", bad)
