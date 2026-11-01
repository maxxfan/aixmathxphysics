#!/usr/bin/env python3
"""Column-generation 'pricing' of extra NPA words from a saved L3 interior-point solution (work/upper/*.npz).

Idea: the saved z (multipliers of the non-identity moment classes) gives the primal moment matrix Gamma_ij = -z_class(ij)
(Gamma_00 = 1) of the L3 optimum.  A new word u cuts the L3 optimum off iff the known part of its moment row cannot be
completed to a PSD matrix.  For u let K = {i in W : class(W_i^dag u) is an existing L3 moment}; c_i = moment(W_i^dag u),
d = moment(u^dag u) (always known for words with na<=3, nb<=3).  Then
        schur(u) = d - c_K^T (Gamma_KK + eps I)^-1 c_K      (< 0  <=>  Gamma cannot be extended: u is a useful word)
Report the most negative schur(u) (violation = -schur).  Read-only: uses work/upper/*.npz, writes nothing but stdout.

    python3 price_words.py J97 work/upper/J97_L3.npz [--shapes 2,2;3,1;1,3;2,1;1,2;4,0;0,4;3,0;0,3] [--top 40] [--eps 1e-9]
"""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import argparse
import itertools
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "qlib"))
import bell   # noqa: E402
import certs  # noqa: E402


def key(w):
    return "".join(str(i + 1) for i in w[0]) + "|" + "".join(str(i + 1) for i in w[1])


def alt_words(m, n):
    return [a for a in itertools.product(range(m), repeat=n) if all(a[i] != a[i + 1] for i in range(n - 1))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target")
    ap.add_argument("npz")
    ap.add_argument("--level", default="3")
    ap.add_argument("--shapes", default="2,2;3,1;1,3;2,1;1,2;4,0;0,4")
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--eps", type=float, default=1e-9)
    ap.add_argument("--chosen", default="", help="SA|SB spec (ipm.py --extra format) to report the scores of")
    ap.add_argument("--extra", default="", help="words already in the saved solution (ipm.py --extra format), for re-pricing")
    ap.add_argument("--emit", type=int, default=0, help="print an ipm.py --extra spec of the N best words")
    ap.add_argument("--lowest", type=int, default=0, help="also list the N words with smallest X_ii of the saved solution")
    ap.add_argument("--dropwords", default="", help="exact words dropped from L3 in the saved solution (ipm4.py --dropwords format)")
    a = ap.parse_args()
    c = bell.load(a.target)
    W = bell.words(c, a.level)
    for part in filter(None, a.extra.split(";")):
        sa, sb = part.split("|")
        for x in (sa.split(",") if sa else [""]):
            for yb in (sb.split(",") if sb else [""]):
                w = (tuple(int(ch) - 1 for ch in x), tuple(int(ch) - 1 for ch in yb))
                if w not in set(W):
                    W.append(w)
    dw = set()
    for part in filter(None, a.dropwords.split(";")):
        sa, sb = part.split("|")
        for x in (sa.split(",") if sa else [""]):
            for yb in (sb.split(",") if sb else [""]):
                dw.add((tuple(int(ch) - 1 for ch in x), tuple(int(ch) - 1 for ch in yb)))
    W = [w for w in W if w not in dw]
    n = len(W)
    cl, inv = certs.classes(W)
    one = cl[0, 0]
    keep = np.array([k for k in range(len(inv)) if k != one])
    z = np.load(a.npz)["z"]
    assert len(z) == len(keep), (len(z), len(keep))
    y = np.zeros(len(inv))
    y[one] = 1.0
    y[keep] = -z
    G = y[cl]
    ev = np.linalg.eigvalsh(G)
    print(f"{a.target} L{a.level}: n={n} classes={len(inv)}  Gamma min eig {ev[0]:.3e}  rank(>1e-7) {int((ev > 1e-7).sum())}")
    t = certs.targets(c, inv)
    print(f"check bell value b.z - b0 = {float(t[keep] @ z) - c['b0']:.10f}")
    idx = {m: k for k, m in enumerate(inv)}
    seen = set(W)
    shapes = [tuple(int(x) for x in s.split(",")) for s in a.shapes.split(";")]
    cands = []
    for na, nb in shapes:
        for aw in alt_words(c["mA"], na):
            for bw in alt_words(c["mB"], nb):
                w = (aw, bw)
                if w not in seen:
                    cands.append(w)

    def score(u):
        cvec = np.full(n, np.nan)
        for i in range(n):
            k = idx.get(bell.canon_real(bell.product(W[i], u)))
            if k is not None:
                cvec[i] = y[k]
        K = np.where(~np.isnan(cvec))[0]
        kd = idx.get(bell.canon_real(bell.product(u, u)))
        if kd is None:
            return None
        d = y[kd]
        GK = G[np.ix_(K, K)] + a.eps * np.eye(len(K))
        cK = cvec[K]
        s = d - cK @ np.linalg.solve(GK, cK)
        return s, len(K)

    res = []
    for u in cands:
        r = score(u)
        if r is not None:
            res.append((r[0], r[1], u))
    res.sort(key=lambda x: x[0])
    print(f"{len(cands)} candidate words, {len(res)} with known d; most violating (schur<0):")
    for s, k, u in res[:a.top]:
        print(f"  {key(u):>14s}  shape {len(u[0])},{len(u[1])}  |K|={k:3d}  schur {s:+.4e}")
    by = {}
    for s, k, u in res:
        by.setdefault((len(u[0]), len(u[1])), []).append(s)
    for sh, v in sorted(by.items()):
        v = np.array(v)
        print(f"shape {sh}: n={len(v)} min {v.min():+.3e} median {np.median(v):+.3e} #(schur<-1e-6) {int((v < -1e-6).sum())}")
    if a.emit:
        neg = [u for s, k, u in res if s < 0][:a.emit]
        print("EMIT --extra '" + ";".join("".join(str(i + 1) for i in u[0]) + "|" + "".join(str(i + 1) for i in u[1]) for u in neg) + "'")
    if a.lowest:
        d = np.diag(np.load(a.npz)["X"])
        assert len(d) == n, (len(d), n)
        order = np.argsort(d)
        print("LOWEST X_ii (ipm4.py --dropwords candidates): " + " ".join(f"{key(W[i])}:{d[i]:.3f}" for i in order[:a.lowest] if i > 0))
    if a.chosen:
        from_spec = []
        for part in filter(None, a.chosen.split(";")):
            sa, sb = part.split("|")
            for x in (sa.split(",") if sa else [""]):
                for yb in (sb.split(",") if sb else [""]):
                    from_spec.append((tuple(int(ch) - 1 for ch in x), tuple(int(ch) - 1 for ch in yb)))
        print("chosen words (e.g. upper's set):")
        rank = {u: r for r, (_, _, u) in enumerate(res)}
        for u in from_spec:
            r = score(u)
            print(f"  {key(u):>14s}  schur {r[0]:+.4e}  rank {rank.get(u, -1)}/{len(res)}")


if __name__ == "__main__":
    main()
