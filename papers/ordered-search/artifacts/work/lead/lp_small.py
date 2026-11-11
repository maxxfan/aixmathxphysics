"""Dense grid LP for the FGGS chain (small n): maximize beta = min_t min_grid q_t.

usage: python lp_small.py K N [R]
beta* >= 0 on the grid is necessary for an exact K-query algorithm on N elements.
"""
import sys
import pathlib
import time
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "qlib"))
import fggs  # noqa: E402


def eval_matrix(n, R):
    """dense (nR+1) x (n+1) map: half-grid values -> fine values"""
    I = np.eye(n + 1)
    cols = []
    for m in range(n + 1):
        c = fggs.coeffs_from_values(I[m])
        c[n] = 0.0  # degree < n is enforced by the mean constraints; drop the Nyquist term
        cols.append(fggs.values_from_coeffs(c, R))
    return np.array(cols).T


def solve(k, n, R=8, verbose=True):
    ev, od = fggs.half_parity(n)
    mult = fggs.circle_mult(n)
    A = eval_matrix(n, R)
    # Nyquist term: with c[n] dropped, the map is only exact when both halves have mean 1, which
    # the constraints below enforce.
    par = lambda t: ev if t % 2 == 1 else od  # parity points of u_t (u_0 on odd)
    F = fggs.fejer_values(n)
    u0 = F[od]
    ukm1 = np.ones(par(k - 1).sum())
    sizes = [par(t).sum() for t in range(1, k - 1)]
    offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    nv = offs[-1] + 1  # + beta
    rows, rhs = [], []
    for t in range(1, k):
        Ma = A[:, par(t - 1)]
        Mb = A[:, par(t)]
        G = np.zeros((A.shape[0], nv))
        h = np.zeros(A.shape[0])
        # constraint: -(Ma u_{t-1} + Mb u_t) + beta <= 0
        if t - 1 == 0:
            h += Ma @ u0
        else:
            G[:, offs[t - 2]:offs[t - 1]] -= Ma
        if t == k - 1:
            h += Mb @ ukm1
        else:
            G[:, offs[t - 1]:offs[t]] -= Mb
        G[:, -1] = 1.0
        rows.append(G)
        rhs.append(h)
    Aub = np.vstack(rows)
    bub = np.concatenate(rhs)
    Aeq = np.zeros((k - 2, nv))
    for t in range(1, k - 1):
        Aeq[t - 1, offs[t - 1]:offs[t]] = mult[par(t)]
    beq = np.full(k - 2, float(n))
    cobj = np.zeros(nv)
    cobj[-1] = -1.0
    bounds = [(0, None)] * (nv - 1) + [(None, None)]
    t0 = time.time()
    res = linprog(cobj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs")
    if verbose:
        print(f"k={k} n={n} R={R}: status={res.status} beta*={-res.fun if res.status == 0 else None}"
              f"  ({time.time() - t0:.1f}s, {Aub.shape})")
    us = [u0] + [res.x[offs[i]:offs[i + 1]] for i in range(k - 2)] + [ukm1] if res.status == 0 else None
    return res, us


if __name__ == "__main__":
    k, n = int(sys.argv[1]), int(sys.argv[2])
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    solve(k, n, R)
