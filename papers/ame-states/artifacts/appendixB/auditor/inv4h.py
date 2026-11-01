"""Involution-symmetric AME(8,4) search for STRUCTURED AME(4,4) codes C: F_16-linear Hermitian self-dual MDS [4,2,3]_16 codes (qlib canonical_from_hermitian).
Same reduction as inv4.py (C = S cap V^sigma = N S, beta symmetric bilinear on C), DFS in C (inv4dfs).  usage: python3 work/auditor/inv4h.py NSAMPLES SEED"""
import os, sys, random, json, time, subprocess, itertools
os.environ.setdefault("OMP_NUM_THREADS", "2")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "lib"))   # [release] was ROOT (qlib package)
import numpy as np
import inv4
import ame   # [release] was: from qlib import ame
WORKDIR = os.environ.get("APPB_WORK", ".")   # [release] logs / .bin / cands go here instead of work/auditor
DFSBIN = os.path.join(os.environ["BUILD"], "inv4dfs") if os.environ.get("BUILD") else os.path.join(HERE, "inv4dfs")  # [release]

nsamp = int(sys.argv[1]); seed = int(sys.argv[2])
rng = random.Random(seed)
F = ame.GF(16)
q = 4
log = open(os.path.join(WORKDIR, f"inv4h_{seed}.log"), "a")   # [release]


def norm(x):
    return F.pow(x, q + 1)


def hermitian_mds_codes():
    """all (a,b,c,d) giving G = [[1,0,a,b],[0,1,c,d]] Hermitian self-dual (rows) with all entries and det nonzero (MDS [4,2,3])."""
    res = []
    minus1 = F.neg(1)
    for a, b in itertools.product(range(1, 16), repeat=2):
        if F.add(F.add(1, norm(a)), norm(b)) != 0:
            continue
        for c, d in itertools.product(range(1, 16), repeat=2):
            if F.add(F.add(1, norm(c)), norm(d)) != 0:
                continue
            if F.add(F.mul(a, F.pow(c, q)), F.mul(b, F.pow(d, q))) != 0:
                continue
            det = F.sub(F.mul(a, d), F.mul(b, c))
            if det == 0:
                continue
            # MDS: all 2x2 minors of [I | A] nonzero <=> entries of A and det nonzero (A = [[a,b],[c,d]] since [I|A] 2 x 4)
            # here G = [I2 | A] with A = [[a,b],[c,d]]
            res.append((a, b, c, d))
    return res


codes = hermitian_mds_codes()
log.write(json.dumps({"n_hermitian_mds_codes": len(codes)}) + "\n"); log.flush()
rng.shuffle(codes)
done = 0
t0 = time.time()
for (a, b, c, d) in codes:
    if done >= nsamp:
        break
    G = [[1, 0, a, b], [0, 1, c, d]]
    rows = ame.canonical_from_hermitian(G, F, q)       # 8 rows (x(8)|z(8))
    g = []
    for r in rows:
        v = 0
        for qb in range(8):
            t, k = divmod(qb, 2)
            v |= (r[qb] << (2 * k)) << (4 * t)
            v |= (r[8 + qb] << (2 * k + 1)) << (4 * t)
        g.append(v)
    # check Lagrangian and weights
    assert all(inv4.B4(g[i], g[j]) == 0 for i in range(8) for j in range(8))
    Call = inv4.code_elements(g)
    assert len(set(Call.tolist())) == 256
    if inv4.party_weights4(Call[1:]).min() < 3:
        continue
    # dual basis h_j : B(h_j, g_i) = delta_ij, by Gaussian elimination over F_2 on the 8 x 16 matrix M (row i = swapJ(g_i))
    M = [inv4.swapJ(gi) for gi in g]
    h = []
    for j in range(8):
        # solve parity(y & M[i]) = delta_ij  for y in F_2^16  (brute force over reduced system)
        # elimination on augmented rows
        rowsA = [(M[i], 1 if i == j else 0) for i in range(8)]
        piv = []
        red = []
        for (m, rhs) in rowsA:
            for (pm, pr, pc) in red:
                if (m >> pc) & 1:
                    m ^= pm; rhs ^= pr
            if m == 0:
                assert rhs == 0
                continue
            pc = m.bit_length() - 1
            red = [((x ^ m) if (x >> pc) & 1 else x, (r ^ rhs) if (x >> pc) & 1 else r, cc) for (x, r, cc) in red]
            red.append((m, rhs, pc))
        y = 0
        for (pm, pr, pc) in red:
            if pr:
                y |= 1 << pc
        assert all((bin(y & M[i]).count("1") & 1) == (1 if i == j else 0) for i in range(8))
        h.append(y)
    good = np.zeros((256, 256), dtype=bool)
    Wn = np.stack([(Call >> (4 * t)) & 15 for t in range(4)], axis=-1)
    for lam in range(256):
        y = 0
        for j in range(8):
            if (lam >> j) & 1:
                y ^= h[j]
        Z = y ^ Call
        Zn = np.stack([(Z >> (4 * t)) & 15 for t in range(4)], axis=-1)
        wt1 = (Zn != 0).sum(axis=1)
        Z2 = Zn[None, :, :] ^ Wn[:, None, :]
        wt = wt1[None, :] + (Z2 != 0).sum(axis=2)
        good[:, lam] = wt.min(axis=1) >= 5
    good[0, :] = True
    fn = os.path.join(WORKDIR, f"goodh_{seed}.bin")   # [release]
    open(fn, "wb").write(good.astype(np.uint8).tobytes())
    out = subprocess.run([DFSBIN, fn, "1"], capture_output=True, text=True, timeout=600).stdout.split()
    nodes, nsol = int(out[0]), int(out[1])
    rec = {"sample": done, "abcd": [a, b, c, d], "frac_good": round(float(good[1:].mean()), 4), "nodes": nodes, "sols": nsol, "t": round(time.time() - t0, 1)}
    if nsol:
        beta = [int(x) for x in out[2:10]]
        words = []
        for i in range(8):
            y = 0
            for j in range(8):
                if (beta[i] >> j) & 1:
                    y ^= h[j]
            words.append((y, y ^ g[i]))
        for i in range(8):
            words.append((g[i], g[i]))
        import c7sym as C7
        ws = [[inv4.nib(x, t) for t in range(4)] + [inv4.nib(y2, t) for t in range(4)] for x, y2 in words]
        rws = C7.to_rows(ws, None)
        bad, cuts = ame.bad_cuts(rws, 2, 2, 8)
        rec.update({"bad_cuts": bad, "beta": beta})
        if bad == 0:
            os.makedirs(os.path.join(WORKDIR, "cands"), exist_ok=True); ame.write_stabilizer(os.path.join(WORKDIR, "cands", f"inv4h_hit_{seed}.json"), 8, 4, 2, rws, script=f"work/auditor/inv4h.py {nsamp} {seed}", note="involution-symmetric AME(8,4), AUDITOR search")
    log.write(json.dumps(rec) + "\n"); log.flush()
    done += 1
log.write(json.dumps({"DONE": done, "elapsed": round(time.time() - t0, 1)}) + "\n")
log.close()
