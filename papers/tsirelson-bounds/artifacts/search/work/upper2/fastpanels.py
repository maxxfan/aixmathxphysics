"""fastpanels.py: drop-in replacement of ipm_half_u.schur_panels exploiting that S_k X has only |rows(S_k)| nonzero rows,
so  Zinv S_k X = Zinv[:, R_k] @ (S_k X)[R_k, :]   (cost n^2 |R_k| instead of 2 n^3 per constraint; sum_k |R_k| <= n^2).
    import fastpanels; fastpanels.install()      # patches ipm_half_u.schur_panels
Self-test:  python fastpanels.py TARGET [LEVEL] [EXTRA]   (compares against the dense formation on a small problem)
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "qlib"))
sys.path.insert(0, os.path.join(ROOT, "work", "upper"))
import ipm_half_u as hu  # noqa: E402

_CACHE = {}


def _rows(P):
    """per chunk: list over constraints of the unique row indices touched by S_k (computed once per Prob)."""
    key = id(P)
    if key in _CACHE:
        return _CACHE[key]
    n = P.n
    out = []
    for l0, l1, G in P.chunks:
        c = l1 - l0
        nz = np.diff(G.indptr).reshape(c, n) > 0
        out.append([np.flatnonzero(nz[k]) for k in range(c)])
    _CACHE.clear()
    _CACHE[key] = out
    return out


def schur_panels_fast(P, X, Zinv):
    n, m = P.n, P.m
    panels, bounds = [], []
    ch = P.chunks
    R = _rows(P)
    GROUP = hu.GROUP
    for g in range(0, len(ch), GROUP):
        grp = ch[g:g + GROUP]
        L0, L1 = grp[0][0], grp[-1][1]
        PmS = P.Pm[L0:]
        pan = np.empty((m - L0, L1 - L0), order="F")
        for gi, (l0, l1, G) in enumerate(grp):
            c = l1 - l0
            H = (G @ X).reshape(c, n, n)
            V = np.zeros((c, n, n))
            for k, r in enumerate(R[g + gi]):
                if len(r) == 0:
                    continue
                V[k] = Zinv[:, r] @ H[k][r]
            pan[:, l0 - L0:l1 - L0] = PmS @ V.reshape(c, n * n).T
        panels.append(pan)
        bounds.append((L0, L1))
    return panels, bounds


def install():
    hu.schur_panels = schur_panels_fast


if __name__ == "__main__":
    os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
    import bell
    import ipm
    c = bell.load(sys.argv[1])
    level = sys.argv[2] if len(sys.argv) > 2 else "3"
    W = list(bell.words(c, level))
    ex = sys.argv[3] if len(sys.argv) > 3 else ""
    if ex.startswith("@"):
        ex = [l[6:].strip() for l in open(ex[1:]) if l.startswith("EXTRA ")][0]
    for w in ipm.parse_extra(ex):
        if w not in set(W):
            W.append(w)
    P = ipm.Prob(c, W, chunk=4)
    hu.GROUP = 32
    n = P.n
    rng = np.random.default_rng(1)
    A = rng.standard_normal((n, n)); X = A @ A.T / n + np.eye(n)
    B = rng.standard_normal((n, n)); Zinv = B @ B.T / n + np.eye(n)
    t0 = time.time(); p1, b1 = hu.schur_panels(P, X, Zinv); t1 = time.time() - t0
    t0 = time.time(); p2, b2 = schur_panels_fast(P, X, Zinv); t2 = time.time() - t0
    err = max(np.abs(a - b).max() for a, b in zip(p1, p2))
    sc = max(np.abs(a).max() for a in p1)
    print(f"n={n} m={P.m}: dense {t1:.2f}s fast {t2:.2f}s speedup {t1 / t2:.1f}x max|diff| {err:.2e} (scale {sc:.2e})")
