"""smooth.py for large n (k = 7, n ~ 10^6).

Same max-min search, with three changes that keep memory and time near O(n log n):
  * extra points (exact local minima of earlier passes) are sparse rows that interpolate the
    oversampled grid with a Gaussian-windowed sinc instead of dense O(n) rows;
  * local minima are Newton-refined with q', q'' from their own oversampled grids, interpolated the
    same way, instead of dense cos/sin sums;
  * the k - 1 transforms of a value/gradient evaluation run as one batched (threaded) DCT.
The best chain is saved after every pass, so a long run can be stopped and certified.

usage: python smooth_big.py K N [--m0 100] [--rho 1.05] [--R 8] [--tag T]
"""
import argparse
import os
import pathlib
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy import fft as sfft
from scipy.optimize import minimize

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "qlib"))
import fggs  # noqa: E402
import lp_basis  # noqa: E402
from interp import Interp  # noqa: E402
from smooth import omega  # noqa: E402

WORKERS = min(8, os.cpu_count() or 1)


class SmoothBig:
    def __init__(self, k, n, m0, rho, R, W=24):
        self.k, self.n, self.R = k, n, R
        ev, od = fggs.half_parity(n)
        self.par = lambda t: ev if t % 2 == 1 else od
        mult = fggs.circle_mult(n)
        self.fixed = {0: fggs.fejer_values(n)[od], k - 1: np.ones(self.par(k - 1).sum())}
        self.Phi = {t: lp_basis.basis(n, self.par(t), m0, rho) for t in range(1, k - 1)}
        sizes = [self.Phi[t].shape[1] for t in range(1, k - 1)]
        self.offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
        self.nv = int(self.offs[-1])
        self.crow = {t: self.Phi[t].T @ mult[self.par(t)] for t in range(1, k - 1)}
        N = n * R
        self.N = N
        self.winv = 1.0 / omega(np.pi * np.arange(N + 1) / N, n)
        self.W1 = np.full(n + 1, 2.0)
        self.W1[[0, n]] = 1.0
        self.WR = np.full(N + 1, 2.0)
        self.WR[[0, N]] = 1.0
        self.fixmask = np.zeros((k - 1, N + 1), bool)
        self.fixmask[0, R * np.nonzero(self.par(0))[0]] = True
        self.fixmask[k - 2, R * np.nonzero(self.par(k - 1))[0]] = True
        # 1/omega with +inf at the fixed samples (their values are positive), so q * winv_m masks
        # them without a boolean scatter on every evaluation
        self.winv_m = np.tile(self.winv, (k - 1, 1))
        self.winv_m[self.fixmask] = np.inf
        self.interp = Interp(n, R, W)
        self.extra = {t: [sp.csr_matrix((0, N + 1)), np.zeros(0), np.zeros(0)] for t in range(1, k)}

    def add_extra(self, t, thetas):
        B0, w0, th0 = self.extra[t]
        self.extra[t] = [sp.vstack([B0, self.interp.matrix(thetas)]).tocsr(),
                         np.r_[w0, 1.0 / omega(thetas, self.n)], np.r_[th0, thetas]]

    def n_extra(self):
        return sum(len(e[1]) for e in self.extra.values())

    def us(self, a):
        u = [self.fixed[0]]
        for t in range(1, self.k - 1):
            u.append(self.Phi[t] @ a[self.offs[t - 1]:self.offs[t]])
        u.append(self.fixed[self.k - 1])
        return u

    def fine(self, V):
        """(k-1, n+1) half-grid values -> (k-1, N+1) values on the R-times finer grid"""
        n = self.n
        C = sfft.dct(V, type=1, axis=-1, workers=WORKERS) / (2 * n)
        X = np.zeros((V.shape[0], self.N + 1))
        X[:, :n] = C[:, :n]
        return sfft.dct(X, type=1, axis=-1, workers=WORKERS)

    def fine_adj(self, G):
        n = self.n
        Z = self.WR * sfft.dct(G / self.WR, type=1, axis=-1, workers=WORKERS)
        Y = np.zeros((G.shape[0], n + 1))
        Y[:, :n] = Z[:, :n]
        return self.W1 * sfft.dct(Y / self.W1, type=1, axis=-1, workers=WORKERS) / (2 * n)

    def value_grad(self, a, mu, cap=np.inf):
        """softmin_mu of min(q/omega, cap) over grid + extra points: with a finite cap, points
        above it stop pulling (feasibility mode), so the chain settles while dips are fixed"""
        k = self.k
        V = np.array(fggs.chain_values(self.us(a), self.n))
        Q = self.fine(V)
        S = Q * self.winv_m
        X = [(B @ Q[t - 1]) * w for t, (B, w, _) in self.extra.items()]
        smin = min(S.min(), min((x.min() for x in X if len(x)), default=np.inf), cap)
        over = S > cap if cap < np.inf else None
        if over is not None:
            S[over] = cap
        np.subtract(S, smin, out=S)
        np.multiply(S, -1.0 / mu, out=S)
        np.exp(S, out=S)
        EX = [np.exp(-(np.minimum(x, cap) - smin) / mu) for x in X]
        Z = S.sum() + sum(e.sum() for e in EX)
        F = smin - mu * np.log(Z)
        if over is not None:
            S[over] = 0.0
            EX = [np.where(x > cap, 0.0, e) for x, e in zip(X, EX)]
        S *= self.winv / Z                                   # dF/dQ from the grid
        for t, (B, w, _) in self.extra.items():
            if len(w):
                S[t - 1] += B.T @ (EX[t - 1] * w / Z)
        GV = self.fine_adj(S)
        g = np.zeros(self.nv)
        for t in range(1, k - 1):
            p = self.par(t)
            gu = GV[t - 1][p] + GV[t][p]                     # u_t enters q_t and q_{t+1}
            gt = self.Phi[t].T @ gu
            r = self.crow[t]
            g[self.offs[t - 1]:self.offs[t]] = gt - (gt @ r) / (r @ r) * r
        return F, g, smin

    def init_guess(self):
        """u_t = U^(1 - t/(k-1)) normalized, U = Fejer envelope; fitted to the basis"""
        n, k = self.n, self.k
        mult = fggs.circle_mult(n)
        a = np.zeros(self.nv)
        for t in range(1, k - 1):
            m = np.nonzero(self.par(t))[0]
            th = np.pi * m / n
            U = 1.0 / (n * np.sin(th / 2) ** 2 + np.pi ** 2 / (4 * n))
            u = U ** (1 - t / (k - 1))
            u *= n / (mult[m] @ u)
            P = self.Phi[t]
            sol = np.linalg.solve((P.T @ P).toarray(), P.T @ u)
            r = self.crow[t]
            sol *= n / (r @ sol)
            a[self.offs[t - 1]:self.offs[t]] = sol
        return a


def local_minima(c, n, Rc=16, W=32, keep_quantile=1.0, slack=np.inf, newton=5, chunk=200000):
    """Newton-refined local minima of q (coefficients c, c[n] = 0) on [0, pi] whose parabolic
    estimate is within `slack` of the lowest ones; returns (theta, q/omega).  The default refines
    every grid minimum: at large theta omega ~ 1/n is below the parabolic estimate's error
    (~1e-5 at Rc = 16), so a filter on the estimate drops real negative minima (a band of 20 at
    k = 6, n = 90000 was missed that way)."""
    N = n * Rc
    vals = fggs.values_from_coeffs(c, Rc)
    h = np.pi / N
    i = np.arange(1, N)
    ismin = (vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1])
    d1 = vals[i + 1] - vals[i - 1]
    d2 = vals[i + 1] - 2 * vals[i] + vals[i - 1]
    est = (vals[i] - d1 * d1 / (8 * np.maximum(d2, 1e-300))) / omega(i * h, n)
    cand = i[ismin & (est < np.quantile(est[ismin], keep_quantile) + slack)]
    j = np.arange(n)
    x = np.zeros(N + 1)
    x[:n] = j * j * c[:n]
    g2 = -sfft.dct(x, type=1, workers=WORKERS)               # q''
    y = np.zeros(N - 1)
    y[:n - 1] = j[1:] * c[1:n]
    g1 = np.zeros(N + 1)
    g1[1:N] = -sfft.dst(y, type=1, workers=WORKERS)          # q'
    I = Interp(n, Rc, W)
    ths, qvs = [], []
    for s in range(0, len(cand), chunk):
        th = cand[s:s + chunk] * h
        for _ in range(newton):
            q1 = I.matrix(th, odd=True) @ g1
            q2 = I.matrix(th) @ g2
            ok = q2 > 0
            step = np.clip(q1 / np.where(ok, q2, 1.0), -2 * h, 2 * h)
            th = np.where(ok, np.clip(th - step, 0, np.pi), th)
        ths.append(th)
        qvs.append(I.matrix(th) @ vals)
    th = np.concatenate(ths) if ths else np.zeros(0)
    qv = np.concatenate(qvs) if qvs else np.zeros(0)
    return th, qv / omega(th, n)


def run(k, n, m0=100, rho=1.05, R=8, mus=(1.0, 0.3, 0.1, 0.03, 0.01, 0.003), maxiter=3000,
        early_maxiter=1500, polish_maxiter=3000, outer=24, polish_mus=(0.01, 0.003), add_below=0.3,
        save=None, init=None, param="lin", cap=np.inf, resume=None, verbose=True):
    """param="lin": a = sc * x with a fixed diagonal scaling sc = |a_init| and projected mean
    constraints.  param="log": every free block is a = n e^x / (r . e^x), so a > 0, the mean
    constraint holds identically and steps are relative (the large-n chains span many decades)."""
    S = SmoothBig(k, n, m0, rho, R)
    a = S.init_guess() if init is None else init
    if resume:                       # chain and extra points of a previous run's last pass
        z = np.load(resume)
        a = z["a"]
        for t in range(1, k):
            if len(z[f"extra{t}"]):
                S.add_extra(t, z[f"extra{t}"])
        mus = polish_mus
    if verbose:
        print(f"k={k} n={n}: {S.nv} variables, grid {S.N + 1} points per q, {WORKERS} fft workers, "
              f"param {param}", flush=True)
    t0 = time.time()
    best = None
    blocks = [slice(S.offs[t - 1], S.offs[t]) for t in range(1, k - 1)]
    if param == "log":
        a = np.maximum(a, 1e-12 * np.abs(a).max())

        def to_a(x):
            a = np.empty_like(x)
            for t, sl in enumerate(blocks, start=1):
                e = np.exp(x[sl] - x[sl].max())
                a[sl] = n * e / (S.crow[t] @ e)
            return a

        def to_x(a):
            return np.log(a)

        def fg(x, mu):
            a = to_a(x)
            F, g, _ = S.value_grad(a, mu, cap)
            gx = np.empty_like(g)
            for t, sl in enumerate(blocks, start=1):
                gx[sl] = a[sl] * (g[sl] - (g[sl] @ a[sl]) / n * S.crow[t])
            return -F, -gx
    else:
        sc = np.maximum(np.abs(a), 1e-6 * np.abs(a).max())

        def to_a(x):
            return sc * x

        def to_x(a):
            return a / sc

        def fg(x, mu):
            F, g, _ = S.value_grad(sc * x, mu, cap)
            g = sc * g
            for t, sl in enumerate(blocks, start=1):
                r = S.crow[t] * sc[sl]
                g[sl] -= (g[sl] @ r) / (r @ r) * r
            return -F, -g

    for o in range(outer):
        for mu in (mus if o == 0 else polish_mus):
            res = minimize(fg, to_x(a), args=(mu,), jac=True, method="L-BFGS-B",
                           options={"maxiter": (polish_maxiter if o else maxiter if mu < 0.1
                                                else early_maxiter), "maxcor": 30})
            a = to_a(res.x)
            F, g, smin = S.value_grad(a, mu)
            if verbose:
                print(f"  [{o}] mu={mu:g}: softmin {F:.5g}, min over grid+extra {smin:.5g}, "
                      f"{res.nit} its ({time.time() - t0:.0f}s)", flush=True)
        worst, added = np.inf, 0
        for t, v in enumerate(fggs.chain_values(S.us(a), n), start=1):
            c = fggs.coeffs_from_values(v)
            c[n] = 0.0
            th, rel = local_minima(c, n)
            worst = min(worst, rel.min())
            new = th[rel < min(smin - 1e-4 * max(1.0, abs(smin)), add_below)]
            if len(new):
                S.add_extra(t, new)
                added += len(new)
        if best is None or worst > best[1]:
            best = (a.copy(), worst, smin)
            if save:
                np.savez(save, k=k, n=n, beta=smin, worst=worst, a=a, m0=m0, rho=rho,
                         **{f"u{i}": v for i, v in enumerate(S.us(a))})
        if save:
            np.savez(str(save).replace(".npz", "_last.npz"), k=k, n=n, a=a, worst=worst,
                     **{f"extra{t}": S.extra[t][2] for t in range(1, k)})
        if verbose:
            print(f"  [{o}] true min q/omega {worst:.5g}; extra points +{added} "
                  f"(total {S.n_extra()}); best {best[1]:.5g} ({time.time() - t0:.0f}s)", flush=True)
        if added == 0:
            break
    if verbose:
        print(f"k={k} n={n} nv={S.nv}: true min q/omega {best[1]:.5g} (best pass; "
              f"{time.time() - t0:.0f}s)", flush=True)
    return S, best[0], best[2], best[1]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--m0", type=int, default=100)
    ap.add_argument("--rho", type=float, default=1.05)
    ap.add_argument("--R", type=int, default=8)
    ap.add_argument("--maxiter", type=int, default=3000)
    ap.add_argument("--early", type=int, default=1500)
    ap.add_argument("--polish", type=int, default=3000)
    ap.add_argument("--add-below", type=float, default=0.3)
    ap.add_argument("--param", default="lin", choices=["lin", "log"])
    ap.add_argument("--mus", default="1,0.3,0.1,0.03,0.01,0.003")
    ap.add_argument("--outer", type=int, default=24)
    ap.add_argument("--cap", type=float, default=np.inf)
    ap.add_argument("--resume", default=None, help="_last.npz of a previous run")
    ap.add_argument("--init", default=None, help="npz whose 'a' starts the full mu continuation")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    run(a.k, a.n, a.m0, a.rho, a.R, maxiter=a.maxiter, early_maxiter=a.early,
        polish_maxiter=a.polish, add_below=a.add_below, param=a.param, outer=a.outer, cap=a.cap,
        resume=a.resume, init=np.load(a.init)["a"] if a.init else None,
        mus=tuple(float(x) for x in a.mus.split(",")), save=f"out/big_k{a.k}_n{a.n}{a.tag}.npz")
