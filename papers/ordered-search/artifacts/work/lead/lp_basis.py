"""Reduced-basis cutting-plane LP for the FGGS chain at large n.

Each free half-vector u_t (t = 1..k-2) is u_t = Phi_t a_t with Phi_t = identity on the first `m0`
parity points plus piecewise-linear hats in log(theta) (ratio `rho`) beyond.  Constraints:
  * q_t(theta_m) = u sample >= beta at every half-grid point (sparse rows),
  * q_t(theta) >= beta at cut angles (dense rows, Dirichlet closed form),
  * full-circle mean of every u_t = 1.
Cuts are the local minima of the current q_t on a fine FFT grid, refined by Newton.

usage: python lp_basis.py K N [--m0 200] [--rho 1.03] [--iters 40] [--init FILE.npz]
"""
import argparse
import pathlib
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "qlib"))
import fggs  # noqa: E402


def basis(n, par, m0, rho):
    """sparse Phi (parity points x nb)"""
    m = np.nonzero(par)[0]
    th = np.pi * m / n
    near = m <= m0
    nn = near.sum()
    rows, cols, vals = list(np.nonzero(near)[0]), list(range(nn)), [1.0] * nn
    far = np.nonzero(~near)[0]
    if len(far):
        lo, hi = np.log(th[far[0]]), np.log(np.pi)
        nk = max(2, int(np.ceil((hi - lo) / np.log(rho))) + 1)
        knots = np.linspace(lo, hi, nk)
        x = np.log(th[far])
        pos = np.clip(np.searchsorted(knots, x, side="right") - 1, 0, nk - 2)
        f = (x - knots[pos]) / (knots[pos + 1] - knots[pos])
        rows += list(far) + list(far)
        cols += list(nn + pos) + list(nn + pos + 1)
        vals += list(1 - f) + list(f)
        nb = nn + nk
    else:
        nb = nn
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(m), nb))


def omega(th, n):
    """margin weight: q_t >= beta * omega; omega <= every fixed sample, so beta is not capped at 1/n"""
    s2 = np.sin(np.asarray(th) / 2) ** 2
    return OMEGA_SCALE * np.minimum(1.0, 1.0 / np.maximum(n * s2, 1e-300))


OMEGA_SCALE = 0.1   # the fixed samples cap beta at 1 / OMEGA_SCALE; keep the cap out of reach


def dirichlet(x, n):
    """D(x) = sum_{j=1}^{n-1} cos(j x)"""
    s = np.sin(x / 2)
    small = np.abs(s) < 1e-9
    out = np.where(small, n - 1.0, np.sin((n - 0.5) * x) / (2 * np.where(small, 1.0, s)) - 0.5)
    return out


def eval_rows(n, thetas, parmask):
    """dense rows r(theta) restricted to the half-grid points in parmask: q(theta) = r . v"""
    m = np.nonzero(parmask)[0]
    phi = np.pi * m / n
    w = fggs.circle_mult(n)[m]
    T = thetas[:, None]
    return (1.0 + dirichlet(T - phi, n) + dirichlet(T + phi, n)) * w / (2 * n)


def local_minima(c, R=16, thresh=np.inf):
    n = len(c) - 1
    vals = fggs.values_from_coeffs(c, R)
    i = np.arange(1, len(vals) - 1)
    mins = i[(vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1]) & (vals[i] < thresh)]
    th = np.pi * mins / (n * R)
    jj = np.arange(1, n)
    cc = c[1:n]
    for _ in range(3):
        if len(th) == 0:
            break
        d1 = np.zeros(len(th))
        d2 = np.zeros(len(th))
        for a in range(0, len(th), 256):
            s = np.sin(np.outer(th[a:a + 256], jj))
            co = np.cos(np.outer(th[a:a + 256], jj))
            d1[a:a + 256] = -2 * (s @ (jj * cc))
            d2[a:a + 256] = -2 * (co @ (jj * jj * cc))
        ok = d2 > 0
        th = np.where(ok, np.clip(th - d1 / np.where(ok, d2, 1.0), 0, np.pi), th)
    return th, vals.min()


def binof(th, n, m0, ratio=1.002):
    """bin index: one bin per grid cell for m < 4 m0, log bins beyond"""
    x = th * n / np.pi
    lim = 4 * m0
    return np.where(x < lim, np.floor(x), lim + np.floor(np.log(np.maximum(x, lim) / lim) / np.log(ratio)))


class Chain:
    def __init__(self, k, n, m0, rho):
        self.k, self.n, self.m0 = k, n, m0
        ev, od = fggs.half_parity(n)
        self.par = lambda t: ev if t % 2 == 1 else od
        self.mult = fggs.circle_mult(n)
        F = fggs.fejer_values(n)
        self.fixed = {0: F[od], k - 1: np.ones(self.par(k - 1).sum())}
        self.Phi = {t: basis(n, self.par(t), m0, rho) for t in range(1, k - 1)}
        sizes = [self.Phi[t].shape[1] for t in range(1, k - 1)]
        self.offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
        self.nv = int(self.offs[-1]) + 1
        self.cuts = {t: np.zeros(0) for t in range(1, k)}
        self.cutrows = {t: (np.zeros((0, self.nv)), np.zeros(0)) for t in range(1, k)}

    def us(self, x):
        u = [self.fixed[0]]
        for t in range(1, self.k - 1):
            u.append(self.Phi[t] @ x[self.offs[t - 1]:self.offs[t]])
        u.append(self.fixed[self.k - 1])
        return u

    def rows_for(self, t, thetas):
        """G, h with q_t(theta) = G[:, :-1] x + h (beta column left zero)"""
        G = np.zeros((len(thetas), self.nv))
        h = np.zeros(len(thetas))
        for a in range(0, len(thetas), 512):
            th = thetas[a:a + 512]
            for s in (t - 1, t):
                E = eval_rows(self.n, th, self.par(s))
                if s in self.fixed:
                    h[a:a + 512] += E @ self.fixed[s]
                else:
                    G[a:a + 512, self.offs[s - 1]:self.offs[s]] += (self.Phi[s].T @ E.T).T
        return G, h

    def add_cuts(self, t, thetas):
        if len(thetas) == 0:
            return
        G, h = self.rows_for(t, thetas)
        G0, h0 = self.cutrows[t]
        self.cutrows[t] = (np.vstack([G0, G]), np.concatenate([h0, h]))
        self.cuts[t] = np.concatenate([self.cuts[t], thetas])

    def solve_lp(self, beta_cap=None):
        k, nv = self.k, self.nv
        blocks, rhs = [], []
        for t in range(1, k):
            # sample rows: q_t equals u_t at its parity points (once per free u_t)
            for s in ((t,) if t < k - 1 else ()):
                P = self.Phi[s]
                A = sp.hstack([sp.csr_matrix((P.shape[0], self.offs[s - 1])), -P,
                               sp.csr_matrix((P.shape[0], nv - 1 - self.offs[s])),
                               sp.csr_matrix(omega(np.pi * np.nonzero(self.par(s))[0] / self.n,
                                                   self.n)[:, None])]).tocsr()
                blocks.append(A)
                rhs.append(np.zeros(P.shape[0]))
            G, h = self.cutrows[t]
            if len(h):
                A = -G.copy()
                A[:, -1] = omega(self.cuts[t], self.n)
                blocks.append(sp.csr_matrix(A))
                rhs.append(h)
        Aub = sp.vstack(blocks).tocsr()
        bub = np.concatenate(rhs)
        Aeq = np.zeros((k - 2, nv))
        for t in range(1, k - 1):
            Aeq[t - 1, self.offs[t - 1]:self.offs[t]] = self.Phi[t].T @ self.mult[self.par(t)]
        beq = np.full(k - 2, float(self.n))
        c = np.zeros(nv)
        c[-1] = -1.0
        bounds = [(0, None)] * (nv - 1) + [(None, beta_cap)]
        if getattr(self, "center", None) is not None:   # box trust region around the last iterate
            a = self.center[:-1]
            d = self.trust * np.maximum(np.abs(a), 1e-3 * np.abs(a).max())
            bounds = [(max(0.0, lo), hi) for lo, hi in zip(a - d, a + d)] + [(None, beta_cap)]
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds,
                      method=getattr(self, "method", "highs"))
        if res.status != 0:
            raise RuntimeError(res.message)
        return res.x, -res.fun, Aub.shape

    def separate(self, x, beta, tol, R=16, maxcuts=2000, screen=0.02):
        """add cuts at local minima with q/omega < beta - tol; returns (worst q/omega, #added)"""
        worst, added = np.inf, 0
        n = self.n
        for t, v in enumerate(fggs.chain_values(self.us(x), n), start=1):
            c = fggs.coeffs_from_values(v)
            c[n] = 0.0
            vals = fggs.values_from_coeffs(c, R)
            th_all = np.pi * np.arange(len(vals)) / (n * R)
            rel = vals / omega(th_all, n)
            i = np.arange(1, len(vals) - 1)
            # parabolic vertex through each grid minimum estimates the dip between samples
            d1 = vals[i + 1] - vals[i - 1]
            d2 = vals[i + 1] - 2 * vals[i] + vals[i - 1]
            est = vals[i] - d1 * d1 / (8 * np.maximum(d2, 1e-300))
            est_rel = est / omega(th_all[i], n)
            ismin = (vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1])
            cand = i[ismin & (np.minimum(rel[i], est_rel) < beta + screen)]
            worst = min(worst, rel.min())
            if len(cand) == 0:
                continue
            th = th_all[cand]
            jj = np.arange(1, n)
            cc = c[1:n]
            for _ in range(3):
                d1 = np.concatenate([-2 * (np.sin(th[a:a + 256, None] * jj) @ (jj * cc))
                                     for a in range(0, len(th), 256)])
                d2 = np.concatenate([-2 * (np.cos(th[a:a + 256, None] * jj) @ (jj * jj * cc))
                                     for a in range(0, len(th), 256)])
                ok = d2 > 0
                th = np.where(ok, np.clip(th - d1 / np.where(ok, d2, 1.0), 0, np.pi), th)
            qv = np.concatenate([c[0] + 2 * np.cos(th[a:a + 256, None] * jj) @ cc
                                 for a in range(0, len(th), 256)])
            relq = qv / omega(th, n)
            worst = min(worst, relq.min())
            bad = np.nonzero(relq < beta - tol)[0]
            if len(bad):
                # most violated minimum per bin: one sample per grid cell near the peak,
                # log-spaced (ratio 1.002) beyond, so cuts cover the whole envelope
                b = binof(th[bad], n, self.m0)
                order = np.lexsort((relq[bad], b))
                first = np.r_[True, b[order][1:] != b[order][:-1]]
                bad = bad[order[first]]
                if len(bad) > maxcuts:
                    bad = bad[np.argsort(relq[bad])[:maxcuts]]
            self.add_cuts(t, th[bad])
            added += len(bad)
        return worst, added

    def prune(self, x, beta, keep):
        """drop cut rows whose relative slack exceeds `keep`"""
        dropped = 0
        for t in self.cuts:
            G, h = self.cutrows[t]
            if not len(h):
                continue
            rel = (G[:, :-1] @ x[:-1] + h) / omega(self.cuts[t], self.n) - beta
            m = rel < keep
            dropped += (~m).sum()
            self.cutrows[t] = (G[m], h[m])
            self.cuts[t] = self.cuts[t][m]
        return dropped


def run(k, n, m0=200, rho=1.03, iters=60, tol=1e-6, beta_cap=None, verbose=True, save=None,
        method="highs", keep=None, stop_frac=0.95, seed_far=2000, maxcuts=6000, trust=None):
    ch = Chain(k, n, m0, rho)
    ch.method = method
    near = np.pi * (np.arange(0, min(n, 4 * m0)) + 0.5) / n
    # far field: log-spaced bins, each with 4 cuts spread over one oscillation period 2 pi / n
    nb = max(2, int(np.log(n / (4 * m0)) / np.log(1.01))) if n > 4 * m0 else 0
    centres = np.exp(np.linspace(np.log(near[-1]), np.log(np.pi), nb)) if nb else np.zeros(0)
    far = (centres[:, None] + np.arange(4) * np.pi / (2 * n)).ravel()
    seed = np.concatenate([near, far[far <= np.pi]])
    for t in range(1, k):
        ch.add_cuts(t, seed)
    t0 = time.time()
    best = None
    for it in range(iters):
        x, beta, shape = ch.solve_lp(beta_cap)
        worst, added = ch.separate(x, beta, tol, maxcuts=maxcuts)
        if trust:
            ch.trust = trust
            ch.center = x if best is None or worst > best[1] else best[0]
        dropped = ch.prune(x, beta, keep) if keep and it % 5 == 4 else 0
        if verbose:
            print(f"  it {it}: LP beta={beta:.6g} true min~{worst:.6g} cuts+={added} -{dropped} "
                  f"rows={shape} ({time.time() - t0:.0f}s)", flush=True)
        if best is None or worst > best[1]:
            best = (x, worst, beta)
        if added == 0 or (beta > 0 and worst >= stop_frac * beta):
            break
    x, worst, beta = best
    if save:
        u = ch.us(x)
        np.savez(save, k=k, n=n, beta=beta, worst=worst, m0=m0, rho=rho,
                 **{f"u{i}": v for i, v in enumerate(u)})
    if verbose:
        print(f"k={k} n={n} nv={ch.nv}: LP beta={beta:.6g}, best true min {worst:.6g} "
              f"({time.time() - t0:.0f}s)")
    return ch, x, beta, worst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--m0", type=int, default=200)
    ap.add_argument("--rho", type=float, default=1.03)
    ap.add_argument("--iters", type=int, default=60)
    ap.add_argument("--keep", type=float, default=None)
    ap.add_argument("--trust", type=float, default=None)
    ap.add_argument("--cap", type=float, default=None)
    a = ap.parse_args()
    run(a.k, a.n, a.m0, a.rho, a.iters, beta_cap=a.cap, keep=a.keep, trust=a.trust,
        save=f"out/basis_k{a.k}_n{a.n}.npz")
