"""ipm_half_v.py: ipm_half_u.py + iterative refinement of every Schur solve against the exact operator
M v = A(X At(v) Zinv) (cheap, O(n^3)), and diagonal-shift retry if the panel Cholesky breaks down.
Goal: keep pinf ~1e-10 near the end (ipm_half_u lost accuracy: pinf 4e-9..5e-8 at mu<1e-8, making N*tau of the certificate 1e-7..1e-6).
    ipm_half_v.py --target J97 --level 3 --extra "SA|SB;..." --tol 3e-10 --out cand.json --save x.npz [--nref 3]
"""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import argparse
import sys
import time
from fractions import Fraction

import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "qlib"))
sys.path.insert(0, os.path.join(ROOT, "work", "upper"))
import bell    # noqa: E402
import certs   # noqa: E402
import ipm     # noqa: E402
import mkcand  # noqa: E402
import ipm_half_u as hu  # noqa: E402


def factor(P, X, Zinv):
    """panels + Cholesky; on breakdown retry with diagonal shifts relative to max diag(M)."""
    for rel in (0.0, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10):
        panels, bounds = hu.schur_panels(P, X, Zinv)
        if rel > 0:
            dmax = max(np.diag(p[:p.shape[1]]).max() for p in panels)
            for p in panels:
                c = p.shape[1]
                idx = np.arange(c)
                p[idx, idx] += rel * dmax
        try:
            hu.panel_chol(panels, bounds)
            return panels, bounds, rel
        except np.linalg.LinAlgError:
            del panels
    raise np.linalg.LinAlgError("breakdown even with shift 1e-10")


T_START = time.time()


def solve_ipm(P, tol=1e-10, maxit=60, nref=3, log=print):
    n, m, b = P.n, P.m, P.b
    C = np.zeros((n, n)); C[0, 0] = 1.0
    normA = np.sqrt(np.bincount(P.kidx.ravel(), minlength=m + 1)[:m].astype(float))
    xi = max(10.0, np.sqrt(n), n * np.max((1 + np.abs(b)) / (1 + normA)))
    eta = max(10.0, np.sqrt(n), normA.max(), 1.0)
    X = xi * np.eye(n); Z = eta * np.eye(n); z = np.zeros(m)
    best = None
    for it in range(maxit):
        t0 = time.time()
        rp = b - P.A(X)
        Rd = C - Z - P.At(z)
        pobj = X[0, 0]; dobj = b @ z
        mu = np.sum(X * Z) / n
        relgap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj))
        pinf = np.linalg.norm(rp) / (1 + np.linalg.norm(b)); dinf = np.linalg.norm(Rd) / 2.0
        log(f"it {it:2d} pobj {pobj:.12f} dobj {dobj:.12f} gap {relgap:.2e} pinf {pinf:.2e} dinf {dinf:.2e} mu {mu:.2e} "
            f"rss {ipm.rss_mb():.0f}MB", flush=True)
        score = max(relgap, pinf, dinf)
        if best is None or score < best[0]:
            best = (score, X.copy(), z.copy(), Z.copy(), dict(pobj=pobj, dobj=dobj, mu=mu, pinf=pinf, dinf=dinf, iters=it))
        if score < tol:
            break
        if time.time() - T_START > float(os.environ.get("IPM_WALL", "1e9")):
            log("WALL limit reached; returning best iterate"); break
        try:
            Lz = np.linalg.cholesky(Z)
            Zinv = sla.cho_solve((Lz, True), np.eye(n), check_finite=False)
            Zinv = (Zinv + Zinv.T) / 2
            panels, bounds, rel = factor(P, X, Zinv)

            def msolve(r):
                y = hu.panel_solve(panels, bounds, r)
                for _ in range(nref):
                    res = r - P.A(X @ P.At(y) @ Zinv)
                    y = y + hu.panel_solve(panels, bounds, res)
                return y

            XRZ = P.A(X @ Rd @ Zinv)
            AZ = P.A(Zinv)
            dz = msolve(XRZ + b)
            dZ = Rd - P.At(dz)
            dX = -X - X @ dZ @ Zinv; dX = (dX + dX.T) / 2
            ap = min(1.0, 0.95 * ipm.max_step(X, dX)); ad = min(1.0, 0.95 * ipm.max_step(Z, dZ))
            mu_aff = np.sum((X + ap * dX) * (Z + ad * dZ)) / n
            sigma = min(1.0, (mu_aff / mu) ** 3) if mu > 0 else 0.0
            mut = sigma * mu
            corr = dX @ dZ @ Zinv
            dz = msolve(XRZ - mut * AZ + b + P.A(corr))
            dZ = Rd - P.At(dz)
            dX = mut * Zinv - X - X @ dZ @ Zinv - corr; dX = (dX + dX.T) / 2
            gam = 0.9 + 0.09 * min(ap, ad)
            ap = min(1.0, gam * ipm.max_step(X, dX)); ad = min(1.0, gam * ipm.max_step(Z, dZ))
        except np.linalg.LinAlgError:
            log("numerical breakdown (Cholesky); returning best iterate"); break
        X = X + ap * dX; z = z + ad * dz; Z = Z + ad * dZ
        Z = (Z + Z.T) / 2
        del panels
        log(f"      sigma {sigma:.2e} ap {ap:.3f} ad {ad:.3f} shift {rel:.0e} time {time.time() - t0:.1f}s", flush=True)
    _, X, z, Z, info = best
    return X, z, Z, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--level", default="3")
    ap.add_argument("--extra", default="")
    ap.add_argument("--dropwords", default="", help="L3 words to drop, same SA|SB format as --extra (e.g. '234|;|234')")
    ap.add_argument("--tol", type=float, default=2e-10)
    ap.add_argument("--maxit", type=int, default=50)
    ap.add_argument("--chunk", type=int, default=64)
    ap.add_argument("--nref", type=int, default=3)
    ap.add_argument("--out")
    ap.add_argument("--bits", type=int, default=46)
    ap.add_argument("--save")
    ap.add_argument("--group", type=int, default=hu.GROUP)
    a = ap.parse_args()
    hu.GROUP = a.group
    c = bell.load(a.target)
    dw = set(ipm.parse_extra(a.dropwords))
    W = [w for i, w in enumerate(bell.words(c, a.level)) if i == 0 or w not in dw]
    seen = set(W)
    for w in ipm.parse_extra(a.extra):
        if w not in seen:
            seen.add(w); W.append(w)
    t0 = time.time()
    P = ipm.Prob(c, W, chunk=a.chunk)
    print(f"target {a.target} level {a.level}: n={P.n} m={P.m} (half-M = {4 * P.m ** 2 / 1e6:.0f} MB) setup {time.time() - t0:.1f}s", flush=True)
    X, z, Z, info = solve_ipm(P, tol=a.tol, maxit=a.maxit, nref=a.nref)
    ev = np.linalg.eigvalsh(X)
    print(f"RESULT target={a.target} level={a.level} n={P.n} X00-b0={X[0, 0] - c['b0']:.12f} dual={info['dobj'] - c['b0']:.12f} "
          f"minEig(X)={ev[0]:.3e} pinf={info['pinf']:.2e} time={time.time() - t0:.0f}s maxrss={ipm.rss_mb():.0f}MB", flush=True)
    if a.save:
        np.savez(a.save, X=X, z=z)
    if a.out:
        cand, md = mkcand.make(c, W, X, bits=a.bits, cls=(P.cl, P.inv), author="upper",
                               construction=f"NPA level {a.level} + extra words [{a.extra}], upper's IPM (HKM) with half-memory panel Cholesky (work/upper/ipm_half_v.py, refined Schur solves), rounded, class residual spread evenly",
                               script=f"work/upper/ipm_half_v.py --target {a.target} --level {a.level} --extra \"{a.extra}\" --dropwords \"{a.dropwords}\"")
        certs.save(cand, a.out)
        print(f"saved {a.out} (max class-spread correction {md:.3e}, lambda-b0={float(Fraction(cand['lambda'])) - c['b0']:.12f})", flush=True)


if __name__ == "__main__":
    main()
