"""Auditor: stabilizer AME(8,4) with a PURE involution sigma = (a_1 b_1)(a_2 b_2)(a_3 b_3)(a_4 b_4) (4 transpositions, no fixed party; counting forces t = 4).
Structure (results/argument.md 5.4): with N = sigma - 1, S cap V^sigma = {(x,x): x in C}, N S = {(w,w): w in C}, C a Lagrangian in F_2^16 (4 parties of F_2^4) with party weight >= 3,
i.e. C = stabilizer AME(4,4); S = {(y_w + x, y_w + x + w)}, x in C, where y_w is any vector with B4(y_w, v) = beta(w, v) for all v in C and beta is a SYMMETRIC bilinear
form on C.  weight(s) = sum_t [z_t != 0] + [(z+w)_t != 0] with z = y_w + x, must be >= 5 for all w != 0, x in C.
This script: for a given AME(4,4) code C (random 8-vertex graph state with the 4x4 cut blocks invertible, qubits paired (2t,2t+1) = party t) enumerate ALL symmetric beta by DFS
(rows 0..7, 2^(8-k) choices at level k, pruning with the weight condition for all w with top bit k) and count solutions; a solution is verified with ame.bad_cuts.
usage: python3 work/auditor/inv4.py SEED NSAMPLES [--limit L]
"""
import os, sys, random, itertools, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# [release] the copy of qlib/ame.py lives in ../lib; c7sym.py (imported below) uses "from qlib import ame",
# so register a 'qlib' package alias pointing to the same module.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import ame
import types as _types
if "qlib" not in sys.modules:
    _q = _types.ModuleType("qlib"); _q.__path__ = []; _q.ame = ame
    sys.modules["qlib"] = _q; sys.modules["qlib.ame"] = ame
import c7sym as C7


def nib(v, t):
    return (v >> (4 * t)) & 15


def swapJ(v):
    out = 0
    for t in range(4):
        u = nib(v, t)
        u2 = ((u & 1) << 1) | ((u >> 1) & 1) | (((u >> 2) & 1) << 3) | (((u >> 3) & 1) << 2)
        out |= u2 << (4 * t)
    return out


def B4(y, v):
    return bin(y & swapJ(v)).count("1") & 1


def graph_code(Gam):
    """8 rows of the graph state [I | Gam] as packed ints (party t = qubits 2t, 2t+1; nibble bits x_{2t}, z_{2t}, x_{2t+1}, z_{2t+1})."""
    rows = []
    for i in range(8):
        v = 0
        for q in range(8):
            t, k = divmod(q, 2)
            xbit = 1 if q == i else 0
            zbit = Gam[i][q]
            v |= (xbit << (2 * k)) << (4 * t)
            v |= (zbit << (2 * k + 1)) << (4 * t)
        rows.append(v)
    return rows


def hvec(j):
    t, k = divmod(j, 2)
    return (1 << (2 * k + 1)) << (4 * t)


def code_elements(g):
    out = np.zeros(256, dtype=np.int64)
    for c in range(256):
        v = 0
        for i in range(8):
            if (c >> i) & 1:
                v ^= g[i]
        out[c] = v
    return out


def party_weights4(arr):
    return sum(((arr >> (4 * t)) & 15) != 0 for t in range(4))


def good_table(g, Call):
    """good[w_idx][lam] = (min over c in C of wt8(y_lam + c, y_lam + c + w) >= 5)."""
    h = [hvec(j) for j in range(8)]
    for j in range(8):
        for i in range(8):
            assert B4(h[j], g[i]) == (1 if i == j else 0)
    good = np.zeros((256, 256), dtype=bool)
    Wn = np.stack([(Call >> (4 * t)) & 15 for t in range(4)], axis=-1)   # (256w, 4)
    for lam in range(256):
        y = 0
        for j in range(8):
            if (lam >> j) & 1:
                y ^= h[j]
        Z = y ^ Call                                                     # (256c,)
        Zn = np.stack([(Z >> (4 * t)) & 15 for t in range(4)], axis=-1)  # (256c, 4)
        wt1 = (Zn != 0).sum(axis=1)                                      # (256c,)
        Z2 = Zn[None, :, :] ^ Wn[:, None, :]                             # (256w,256c,4)
        wt = wt1[None, :] + (Z2 != 0).sum(axis=2)
        good[:, lam] = wt.min(axis=1) >= 5
    good[0, :] = True
    return good


def dfs(good, limit):
    sols = []
    beta = [0] * 8      # beta[i] = row i as int (bit j = beta_ij)
    nodes = [0]

    def rec(k):
        if len(sols) >= limit:
            return
        low = 0
        for j in range(k):
            low |= ((beta[j] >> k) & 1) << j
        for free in range(1 << (8 - k)):
            row = low | (free << k)
            beta[k] = row
            ok = True
            # all w with top bit k: w = (1<<k) | s, s subset of lower k bits
            for s in range(1 << k):
                w = (1 << k) | s
                lam = 0
                ww = w
                i = 0
                while ww:
                    if ww & 1:
                        lam ^= beta[i]
                    ww >>= 1
                    i += 1
                if not good[w, lam]:
                    ok = False
                    break
            nodes[0] += 1
            if not ok:
                continue
            if k == 7:
                sols.append(list(beta))
                if len(sols) >= limit:
                    return
            else:
                rec(k + 1)
        beta[k] = 0

    rec(0)
    return sols, nodes[0]


def build_rows(g, Call, beta):
    h = [hvec(j) for j in range(8)]
    words = []
    for i in range(8):
        lam = beta[i]
        y = 0
        for j in range(8):
            if (lam >> j) & 1:
                y ^= h[j]
        words.append((y, y ^ g[i]))
    for i in range(8):
        words.append((g[i], g[i]))
    ws = []
    for a, b in words:
        ws.append([nib(a, t) for t in range(4)] + [nib(b, t) for t in range(4)])
    return C7.to_rows(ws, None)


def main():
    seed = int(sys.argv[1])
    nsamp = int(sys.argv[2])
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 1
    rng = random.Random(seed)
    done = 0
    t0 = time.time()
    tries = 0
    while done < nsamp:
        tries += 1
        Gam = [[0] * 8 for _ in range(8)]
        for i in range(8):
            for j in range(i + 1, 8):
                Gam[i][j] = Gam[j][i] = rng.randint(0, 1)
        g = graph_code(Gam)
        Call = code_elements(g)
        if party_weights4(Call[1:]).min() < 3:
            continue
        good = good_table(g, Call)
        sols, nodes = dfs(good, limit)
        rec = {"sample": done, "tries": tries, "n_good_pairs": int(good[1:].sum()), "frac_good": round(float(good[1:].mean()), 5), "dfs_nodes": nodes, "solutions": len(sols), "sec": round(time.time() - t0, 1)}
        if sols:
            rows = build_rows(g, Call, sols[0])
            bad, cuts = ame.bad_cuts(rows, 2, 2, 8)
            rec["bad_cuts_of_first_solution"] = bad
            rec["iso"] = C7.isotropic(rows)
            rec["rank"] = C7.rank2(rows)
            rec["graph"] = Gam
            rec["beta"] = sols[0]
        print(json.dumps(rec), flush=True)
        done += 1


if __name__ == "__main__":
    main()
