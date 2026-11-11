"""Cutting-plane LP for the FGGS chain (small n): maximize beta = min_t min_theta q_t(theta).

Rows are point evaluations q_t(theta) >= beta at a growing set of angles: a coarse grid plus the
local minima found on a very fine FFT grid after each solve (refined by Newton).
usage: python lp_cut.py K N
"""
import sys
import pathlib
import time
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "qlib"))
import fggs  # noqa: E402


def eval_rows(n, thetas):
    """rows r(theta) with q(theta) = r . v for half-grid values v (degree < n, c[n] dropped)"""
    j = np.arange(n + 1)
    m = np.arange(n + 1)
    e = 2.0 * np.cos(np.outer(thetas, j))
    e[:, 0] = 1.0
    e[:, n] = 0.0
    w = fggs.circle_mult(n)
    C = np.cos(np.pi * np.outer(j, m) / n)   # (n+1) x (n+1)
    return (e @ C) * w / (2 * n)


def local_minima(c, R=64):
    """angles of local minima of q (coeffs c[0..n]) on [0, pi], refined by Newton"""
    n = len(c) - 1
    vals = fggs.values_from_coeffs(c, R)
    i = np.arange(1, len(vals) - 1)
    mins = i[(vals[i] <= vals[i - 1]) & (vals[i] <= vals[i + 1])]
    th = list(np.pi * mins / (n * R))
    if vals[0] <= vals[1]:
        th.append(0.0)
    if vals[-1] <= vals[-2]:
        th.append(np.pi)
    th = np.array(th)
    jj = np.arange(1, n)
    cc = c[1:n]
    for _ in range(4):
        s = np.sin(np.outer(th, jj))
        co = np.cos(np.outer(th, jj))
        d1 = -2 * (s * (jj * cc)).sum(1)
        d2 = -2 * (co * (jj * jj * cc)).sum(1)
        step = np.where(d2 > 0, d1 / np.where(d2 > 0, d2, 1), 0.0)
        th = np.clip(th - step, 0, np.pi)
    qv = c[0] + 2 * (np.cos(np.outer(th, jj)) * cc).sum(1)
    return th, qv


def solve(k, n, R0=4, iters=30, tol=1e-7, verbose=True):
    ev, od = fggs.half_parity(n)
    mult = fggs.circle_mult(n)
    par = lambda t: ev if t % 2 == 1 else od
    F = fggs.fejer_values(n)
    fixed = {0: F[od], k - 1: np.ones(par(k - 1).sum())}
    sizes = [par(t).sum() for t in range(1, k - 1)]
    offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    nv = offs[-1] + 1
    pts = {t: list(np.pi * np.arange(n * R0 + 1) / (n * R0)) for t in range(1, k)}
    Aeq = np.zeros((k - 2, nv))
    for t in range(1, k - 1):
        Aeq[t - 1, offs[t - 1]:offs[t]] = mult[par(t)]
    beq = np.full(k - 2, float(n))
    cobj = np.zeros(nv)
    cobj[-1] = -1.0
    bounds = [(0, None)] * (nv - 1) + [(None, None)]
    t0 = time.time()
    for it in range(iters):
        rows, rhs = [], []
        for t in range(1, k):
            E = eval_rows(n, np.array(pts[t]))
            G = np.zeros((E.shape[0], nv))
            h = np.zeros(E.shape[0])
            for s, ps in ((t - 1, par(t - 1)), (t, par(t))):
                if s in fixed:
                    h += E[:, ps] @ fixed[s]
                else:
                    G[:, offs[s - 1]:offs[s]] -= E[:, ps]
            G[:, -1] = 1.0
            rows.append(G)
            rhs.append(h)
        res = linprog(cobj, A_ub=np.vstack(rows), b_ub=np.concatenate(rhs), A_eq=Aeq, b_eq=beq,
                      bounds=bounds, method="highs")
        beta = -res.fun
        u = [fixed[0]] + [res.x[offs[i]:offs[i + 1]] for i in range(k - 2)] + [fixed[k - 1]]
        worst, added = np.inf, 0
        for t, v in enumerate(fggs.chain_values(u, n), start=1):
            c = fggs.coeffs_from_values(v)
            c[n] = 0.0
            th, qv = local_minima(c)
            worst = min(worst, qv.min())
            bad = th[qv < beta - tol]
            pts[t].extend(bad)
            added += len(bad)
        if verbose:
            print(f"  it {it}: LP beta={beta:.6g} true min={worst:.6g} added={added} "
                  f"rows={sum(len(p) for p in pts.values())} ({time.time() - t0:.0f}s)", flush=True)
        if added == 0:
            break
    if verbose:
        print(f"k={k} n={n}: beta*={beta:.6g} (true min {worst:.6g})")
    return beta, worst, u


def save(path, k, n, u, beta):
    np.savez(path, k=k, n=n, beta=beta, **{f"u{i}": v for i, v in enumerate(u)})


if __name__ == "__main__":
    K, N = int(sys.argv[1]), int(sys.argv[2])
    b, w, u = solve(K, N, iters=int(sys.argv[3]) if len(sys.argv) > 3 else 30)
    save(f"lp_k{K}_n{N}.npz", K, N, u, w)
