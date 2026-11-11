"""Smooth max-min search for the FGGS chain at large n (no LP, no cutting planes).

Maximize softmin_mu over t = 1..k-1 and an R-times oversampled grid of q_t(theta) / omega(theta),
with every free half-vector in the reduced basis of lp_basis.py (samples for m <= m0, log-theta
hats beyond).  All q_t values come from two DCT-I transforms, the gradient from their adjoints,
and the mean-1 constraints are kept by projecting the gradient.  L-BFGS with continuation in mu.

usage: python smooth.py K N [--m0 100] [--rho 1.05] [--R 8] [--init FILE.npz]
"""
import argparse
import pathlib
import sys
import time

import numpy as np
from scipy import fft as sfft
from scipy.optimize import minimize

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "qlib"))
import fggs  # noqa: E402
import lp_basis  # noqa: E402


def omega(th, n):
    s2 = np.sin(np.asarray(th) / 2) ** 2
    return np.minimum(1.0, 1.0 / np.maximum(n * s2, 1e-300))


class Smooth:
    def __init__(self, k, n, m0, rho, R):
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
        self.th = np.pi * np.arange(N + 1) / N
        self.winv = 1.0 / omega(self.th, n)
        self.W1 = np.full(n + 1, 2.0)
        self.W1[[0, n]] = 1.0
        self.WR = np.full(N + 1, 2.0)
        self.WR[[0, N]] = 1.0
        # grid points that are fixed samples (u_0 in q_1, u_{k-1} = 1 in q_{k-1}) do not depend on
        # the variables; the u_0 tail sits exactly at q/omega = 1 and would swamp the softmin
        self.fixmask = np.zeros((k - 1, N + 1), bool)
        self.fixmask[0, R * np.nonzero(self.par(0))[0]] = True
        self.fixmask[k - 2, R * np.nonzero(self.par(k - 1))[0]] = True
        self.chain = lp_basis.Chain(k, n, m0, rho)        # same basis and offsets; rows_for()
        self.extra = {t: (np.zeros((0, self.nv)), np.zeros(0), np.zeros(0), np.zeros(0))
                      for t in range(1, k)}

    def add_extra(self, t, thetas):
        G, h = self.chain.rows_for(t, thetas)
        G0, h0, w0, th0 = self.extra[t]
        self.extra[t] = (np.vstack([G0, G[:, :-1]]), np.r_[h0, h], np.r_[w0, 1.0 / omega(thetas, self.n)],
                         np.r_[th0, thetas])

    def us(self, a):
        u = [self.fixed[0]]
        for t in range(1, self.k - 1):
            u.append(self.Phi[t] @ a[self.offs[t - 1]:self.offs[t]])
        u.append(self.fixed[self.k - 1])
        return u

    def fine(self, v):
        n = self.n
        c = sfft.dct(v, type=1) / (2 * n)
        x = np.zeros(n * self.R + 1)
        x[:n] = c[:n]
        return sfft.dct(x, type=1)

    def fine_adj(self, g):
        n = self.n
        z = self.WR * sfft.dct(g / self.WR, type=1)          # A_R^T g
        y = np.zeros(n + 1)
        y[:n] = z[:n]
        return self.W1 * sfft.dct(y / self.W1, type=1) / (2 * n)  # A_1^T y / 2n

    def value_grad(self, a, mu):
        u = self.us(a)
        vs = fggs.chain_values(u, self.n)
        S = np.array([self.fine(v) * self.winv for v in vs])    # (k-1) x (N+1)
        S[self.fixmask] = np.inf
        X = [(G @ a + h) * w for (G, h, w, _) in (self.extra[t] for t in range(1, self.k))]
        smin = min(S.min(), min((x.min() for x in X if len(x)), default=np.inf))
        E = np.exp(-(S - smin) / mu)
        EX = [np.exp(-(x - smin) / mu) for x in X]
        Z = E.sum() + sum(e.sum() for e in EX)
        F = smin - mu * np.log(Z)
        P = E / Z
        gx = np.zeros(self.nv)
        for t in range(1, self.k):
            G, h, w, _ = self.extra[t]
            if len(h):
                gx += G.T @ (EX[t - 1] / Z * w)
        gu = {t: np.zeros(self.par(t).sum()) for t in range(1, self.k - 1)}
        for t in range(1, self.k):
            gv = self.fine_adj(P[t - 1] * self.winv)
            for s in (t - 1, t):
                if s in gu:
                    gu[s] += gv[self.par(s)]
        g = np.zeros(self.nv)
        for t in range(1, self.k - 1):
            gt = self.Phi[t].T @ gu[t] + gx[self.offs[t - 1]:self.offs[t]]
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
            sol = np.linalg.lstsq(P.toarray() if P.shape[0] * P.shape[1] < 4e7 else P.toarray(), u,
                                  rcond=None)[0]
            r = self.crow[t]
            sol *= n / (r @ sol)
            a[self.offs[t - 1]:self.offs[t]] = sol
        return a


def warm_start(S, path):
    """basis coefficients for S from a solved chain at another n: each free u_t is interpolated in
    (log theta, log u), renormalized to mean 1, and fitted to S's basis by least squares"""
    z = np.load(path)
    n0, n = int(z["n"]), S.n
    ev0, od0 = fggs.half_parity(n0)
    par0 = lambda t: ev0 if t % 2 == 1 else od0
    mult = fggs.circle_mult(n)
    a = np.zeros(S.nv)
    for t in range(1, S.k - 1):
        m0_ = np.nonzero(par0(t))[0]
        th0 = np.pi * np.maximum(m0_, 0.5) / n0
        lu = np.log(np.maximum(z[f"u{t}"], 1e-300))
        m = np.nonzero(S.par(t))[0]
        th = np.pi * np.maximum(m, 0.5) / n
        # near the peak the profile scales with n (Fejer peak ~ n); far away it is a function of theta
        u = np.exp(np.interp(np.log(th), np.log(th0), lu))
        u *= n / (mult[m] @ u)
        P = S.Phi[t]
        sol = np.linalg.lstsq(P.toarray(), u, rcond=None)[0]
        r = S.crow[t]
        sol *= n / (r @ sol)
        a[S.offs[t - 1]:S.offs[t]] = sol
    return a


def run(k, n, m0=100, rho=1.05, R=8, mus=(1.0, 0.3, 0.1, 0.03, 0.01, 0.003), maxiter=3000,
        init=None, verbose=True, save=None, outer=16, polish_mus=(0.01, 0.003), add_below=0.3):
    S = Smooth(k, n, m0, rho, R)
    if isinstance(init, str):
        a = warm_start(S, init)
        mus = (0.03, 0.01, 0.003)
    else:
        a = S.init_guess() if init is None else init
    if verbose:
        print(f"k={k} n={n}: {S.nv} variables, grid {n * R + 1} points per q", flush=True)
    t0 = time.time()
    worst = -np.inf
    best = None
    # diagonal scaling: optimize b with a = sc * b (variables span ~8 orders of magnitude)
    sc = np.maximum(np.abs(a), 1e-6 * np.abs(a).max())

    def fg(b, mu):
        F, g, _ = S.value_grad(sc * b, mu)
        g = sc * g
        for t in range(1, k - 1):                 # project onto the scaled mean constraint
            sl = slice(S.offs[t - 1], S.offs[t])
            r = S.crow[t] * sc[sl]
            g[sl] -= (g[sl] @ r) / (r @ r) * r
        return -F, -g

    for o in range(outer):
        for mu in (mus if o == 0 else polish_mus):
            res = minimize(fg, a / sc, args=(mu,), jac=True, method="L-BFGS-B",
                           options={"maxiter": maxiter, "maxcor": 30})
            a = sc * res.x
            F, g, smin = S.value_grad(a, mu)
            if verbose:
                print(f"  [{o}] mu={mu:g}: softmin {F:.5g}, min over grid+extra {smin:.5g}, "
                      f"{res.nit} its ({time.time() - t0:.0f}s)", flush=True)
        u = S.us(a)
        worst, added = np.inf, 0
        for t, v in enumerate(fggs.chain_values(u, n), start=1):
            c = fggs.coeffs_from_values(v)
            c[n] = 0.0
            th, rel = local_minima_rel(c, n)
            worst = min(worst, rel.min())
            new = th[rel < min(smin - 1e-4 * max(1.0, abs(smin)), add_below)]
            if len(new):
                S.add_extra(t, new)
                added += len(new)
        if best is None or worst > best[1]:
            best = (a.copy(), worst, smin)
        if verbose:
            print(f"  [{o}] true min q/omega {worst:.5g}; extra points +{added} "
                  f"(total {sum(len(e[1]) for e in S.extra.values())})", flush=True)
        if added == 0:
            break
    if verbose:
        print(f"k={k} n={n} nv={S.nv}: true min q/omega {best[1]:.5g} (best pass; "
              f"{time.time() - t0:.0f}s)", flush=True)
    a, worst, smin = best
    if save:
        np.savez(save, k=k, n=n, beta=smin, worst=worst, a=a, m0=m0, rho=rho,
                 **{f"u{i}": v for i, v in enumerate(S.us(a))})
    return S, a, smin, worst


def local_minima_rel(c, n, R=16, keep_quantile=0.05, slack=0.5):
    """Newton-refined local minima of q (coefficients c) on [0, pi] whose parabolic estimate is
    within `slack` of the lowest ones; returns (theta, q/omega)"""
    vals = fggs.values_from_coeffs(c, R)
    th_all = np.pi * np.arange(len(vals)) / (n * R)
    i = np.arange(1, len(vals) - 1)
    ismin = (vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1])
    d1 = vals[i + 1] - vals[i - 1]
    d2 = vals[i + 1] - 2 * vals[i] + vals[i - 1]
    est = (vals[i] - d1 * d1 / (8 * np.maximum(d2, 1e-300))) / omega(th_all[i], n)
    cand = i[ismin & (est < np.quantile(est[ismin], keep_quantile) + slack)]
    th = th_all[cand]
    jj = np.arange(1, n)
    cc = c[1:n]
    for _ in range(4):
        d1 = np.concatenate([-2 * (np.sin(th[a:a + 256, None] * jj) @ (jj * cc)) for a in range(0, len(th), 256)])
        d2 = np.concatenate([-2 * (np.cos(th[a:a + 256, None] * jj) @ (jj * jj * cc)) for a in range(0, len(th), 256)])
        ok = d2 > 0
        th = np.where(ok, np.clip(th - d1 / np.where(ok, d2, 1.0), 0, np.pi), th)
    qv = np.concatenate([c[0] + 2 * np.cos(th[a:a + 256, None] * jj) @ cc for a in range(0, len(th), 256)])
    return th, qv / omega(th, n)


def true_min_rel(c, n, R=16):
    vals = fggs.values_from_coeffs(c, R)
    th_all = np.pi * np.arange(len(vals)) / (n * R)
    i = np.arange(1, len(vals) - 1)
    ismin = (vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1])
    rel = vals / omega(th_all, n)
    d1 = vals[i + 1] - vals[i - 1]
    d2 = vals[i + 1] - 2 * vals[i] + vals[i - 1]
    est = (vals[i] - d1 * d1 / (8 * np.maximum(d2, 1e-300))) / omega(th_all[i], n)
    cand = i[ismin & (est < np.quantile(est[ismin], 0.02) + 0.5)]
    th = th_all[cand]
    jj = np.arange(1, n)
    cc = c[1:n]
    for _ in range(4):
        d1 = np.concatenate([-2 * (np.sin(th[a:a + 256, None] * jj) @ (jj * cc)) for a in range(0, len(th), 256)])
        d2 = np.concatenate([-2 * (np.cos(th[a:a + 256, None] * jj) @ (jj * jj * cc)) for a in range(0, len(th), 256)])
        ok = d2 > 0
        th = np.where(ok, np.clip(th - d1 / np.where(ok, d2, 1.0), 0, np.pi), th)
    qv = np.concatenate([c[0] + 2 * np.cos(th[a:a + 256, None] * jj) @ cc for a in range(0, len(th), 256)])
    return min(rel.min(), (qv / omega(th, n)).min())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--m0", type=int, default=100)
    ap.add_argument("--rho", type=float, default=1.05)
    ap.add_argument("--R", type=int, default=8)
    ap.add_argument("--maxiter", type=int, default=3000)
    ap.add_argument("--init", default=None, help="warm start from a solved chain (npz)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    run(a.k, a.n, a.m0, a.rho, a.R, maxiter=a.maxiter, init=a.init,
        save=f"out/smooth_k{a.k}_n{a.n}{a.tag}.npz")
