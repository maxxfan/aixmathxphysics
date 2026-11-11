#!/usr/bin/env python3
"""Gradient optimisation of permutation-symmetric rank-r inputs for the depolarizing channel
(search tool; certification is verifier/verify.py).

I_c(V, p) and its exact gradient in V (2 x (n+1) complex, Dicke basis):
  f(Y) = sum_k [H(sigma_k) - H(tau_k)]  (nats, homogeneous of degree 1 in Y = V V^dag)
  df   = sum_k < log tau_k - log sigma_k (x) 1_R , L_k(dY) >,   L_k = s_k Phi_d Tr^{2k}
  L_k^dag = s_k E^{2k} Phi_d                                    (Phi_d self-adjoint)
  I_c  = f / Z,  grad_V I_c = (2 G V - 2 I_c V) / Z  with G = sum_k L_k^dag(W_k).

usage:
  depopt.py calib                       re-optimise KG2026 state A at n=45
  depopt.py grow N0 N1 [--seed FILE]    grow a state from N0 to N1, optimising at each n
"""
import argparse
import json
import math
import pathlib
import sys
import time

import numpy as np
from scipy.optimize import brentq, minimize

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = pathlib.Path(__file__).resolve().parent / "states"
LN2 = math.log(2)


def tr1(T):
    N = T.shape[-1] - 1
    c = np.arange(N)
    s0, s1 = np.sqrt(N - c), np.sqrt(c + 1)
    return (np.multiply.outer(s0, s0) * T[..., :N, :N] + np.multiply.outer(s1, s1) * T[..., 1:, 1:]) / N


def emb1(Z):
    N = Z.shape[-1]
    a = np.arange(N + 1)
    s0, s1 = np.sqrt(N - a), np.sqrt(a)
    out = np.zeros(Z.shape[:-2] + (N + 1, N + 1), dtype=Z.dtype)
    out[..., :N, :N] += np.multiply.outer(s0[:N], s0[:N]) * Z
    out[..., 1:, 1:] += np.multiply.outer(s1[1:], s1[1:]) * Z
    return out / N


class Depol:
    def __init__(self, n, p):
        self.n, self.p = n, p
        A = 2 * p * (1 - 2 * p)
        beta = 1 - 4 * p
        self.K = n // 2
        self.cs, self.sk = [], []
        for k in range(self.K + 1):
            d = n - 2 * k
            mj = math.comb(n, k) - (math.comb(n, k - 1) if k else 0)
            # log-space for large n
            lc = [math.lgamma(d + 1) - math.lgamma(s + 1) - math.lgamma(d - s + 1)
                  + (d - s) * math.log(beta) + s * math.log(2 * p) for s in range(d + 1)]
            self.cs.append(np.exp(lc))
            self.sk.append(math.exp(math.log(mj) + k * math.log(A)))

    def phi(self, chain, start, d, cs):
        R = cs[d] * chain[start + d]
        for s in range(d - 1, -1, -1):
            R = cs[s] * chain[start + s] + emb1(R)
        return R

    def blocks(self, Y):
        chain = [Y]
        for _ in range(self.n):
            chain.append(tr1(chain[-1]))
        out = []
        for k in range(self.K + 1):
            d = self.n - 2 * k
            out.append(self.sk[k] * self.phi(chain, 2 * k, d, self.cs[k]))
        return out

    def value_grad(self, V, grad=True):
        r = V.shape[0]
        Y = np.einsum("ia,jb->ijab", V, V.conj())
        Z = float(np.real(np.einsum("iiaa->", Y)))
        f = 0.0
        G = np.zeros_like(Y) if grad else None
        for k, M in enumerate(self.blocks(Y)):
            d1 = M.shape[-1]
            tau = M.transpose(0, 2, 1, 3).reshape(r * d1, r * d1)
            tau = (tau + tau.conj().T) / 2
            sig = np.einsum("iiab->ab", M)
            sig = (sig + sig.conj().T) / 2
            ws, Us = np.linalg.eigh(sig)
            wt, Ut = np.linalg.eigh(tau)
            ps, pt = ws[ws > 0], wt[wt > 0]
            f += float(-(ps * np.log(ps)).sum() + (pt * np.log(pt)).sum())
            if not grad:
                continue
            ls = (Us * np.log(np.maximum(ws, 1e-300))) @ Us.conj().T
            lt = (Ut * np.log(np.maximum(wt, 1e-300))) @ Ut.conj().T
            W = lt.reshape(r, d1, r, d1).transpose(0, 2, 1, 3).copy()
            for q in range(r):
                W[q, q] -= ls
            # L_k^dag(W) = s_k E^{2k} Phi_d(W)
            d = d1 - 1
            ch = [W]
            for _ in range(d):
                ch.append(tr1(ch[-1]))
            U = self.sk[k] * self.phi(ch, 0, d, self.cs[k])
            for _ in range(2 * k):
                U = emb1(U)
            G += U
        I = f / Z
        if not grad:
            return I, None
        g = np.einsum("ijab,jb->ia", G, V)
        return I, (2 * g - 2 * I * V) / Z


def ic_bits(V, p):
    return Depol(V.shape[1] - 1, p).value_grad(V, grad=False)[0] / LN2


def threshold(V, lo=0.060, hi=0.0685, margin=0.0):
    """largest p with I_c(V, p) = margin bits (I_c is nonincreasing in p)"""
    f = lambda p: ic_bits(V, p) - margin
    if f(lo) <= 0:
        return None
    while f(hi) > 0:
        hi += 0.002
    return brentq(f, lo, hi, xtol=1e-10)


def pack(V):
    return np.concatenate([V.real.ravel(), V.imag.ravel()])


def unpack(x, r, n1):
    h = len(x) // 2
    return (x[:h] + 1j * x[h:]).reshape(r, n1)


def optimise(V, p, iters=400, verbose=False):
    r, n1 = V.shape
    model = Depol(n1 - 1, p)

    def fun(x):
        Vx = unpack(x, r, n1)
        I, g = model.value_grad(Vx)
        return -I, -pack(g)

    res = minimize(fun, pack(V), jac=True, method="L-BFGS-B", options={"maxiter": iters, "gtol": 1e-14, "ftol": 1e-16})
    Vn = unpack(res.x, r, n1)
    Vn /= np.linalg.norm(Vn)
    if verbose:
        print(f"   opt p={p:.6f}: I_c={-res.fun / LN2:.4e} bits, {res.nit} its, {res.message}")
    return Vn


def climb(V, rounds=6, iters=400, verbose=True):
    """alternate: optimise I_c at p slightly above the current threshold, recompute threshold"""
    th = threshold(V)
    for _ in range(rounds):
        W = optimise(V, th + 2e-5, iters, verbose)
        t2 = threshold(W)
        if verbose:
            print(f"  threshold {th:.8f} -> {t2}")
        if t2 is None or t2 <= th:
            break
        V, th, gain = W, t2, t2 - th
        if gain < 1e-8:
            break
    return V, th


def grow(V, phi=(1.0, 0.0)):
    """Sym^n -> Sym^{n+1}: rows v -> P_sym(v (x) phi), normalised"""
    al, be = phi
    r, n1 = V.shape
    n = n1 - 1
    a = np.arange(n + 2)
    W = np.zeros((r, n + 2), dtype=complex)
    W[:, :n1] += al * np.sqrt((n + 1 - a[:n1]) / (n + 1)) * V
    W[:, 1:] += be * np.sqrt(a[1:] / (n + 1)) * V
    return W / np.linalg.norm(W)


def save(V, n, th, tag):
    OUT.mkdir(exist_ok=True)
    f = OUT / f"dep_n{n}_{tag}.json"
    f.write_text(json.dumps({"n": n, "threshold_float": th, "rows": [[[z.real, z.imag] for z in row] for row in V]}))
    return f


def load_A():
    raw = (ROOT / "data" / "sources" / "kg2026" / "A.c128").read_bytes()
    return np.frombuffer(raw, dtype="<c16").reshape(2, 46).copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["calib", "grow", "check"])
    ap.add_argument("n0", type=int, nargs="?")
    ap.add_argument("n1", type=int, nargs="?")
    ap.add_argument("--seed")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--tag", default="grow")
    ap.add_argument("--iters", type=int, default=400)
    a = ap.parse_args()
    if a.mode == "check":  # finite-difference gradient check
        rng = np.random.default_rng(0)
        V = rng.normal(size=(2, 9)) + 1j * rng.normal(size=(2, 9))
        m = Depol(8, 0.06)
        I, g = m.value_grad(V)
        x = pack(V)
        gx = pack(g)
        for i in rng.choice(len(x), 6, replace=False):
            e = np.zeros_like(x)
            e[i] = 1e-6
            fd = (m.value_grad(unpack(x + e, 2, 9), False)[0] - m.value_grad(unpack(x - e, 2, 9), False)[0]) / 2e-6
            print(f"coord {i}: analytic {gx[i]:+.8e}  fd {fd:+.8e}")
        return
    if a.mode == "calib":
        V = load_A()
        V /= np.linalg.norm(V)
        t0 = time.time()
        print(f"A: threshold {threshold(V):.8f}")
        V, th = climb(V, a.rounds)
        print(f"calibrated n=45 threshold {th:.8f} ({time.time() - t0:.0f}s) -> {save(V, 45, th, 'calib')}")
        return
    if a.seed:
        d = json.loads(pathlib.Path(a.seed).read_text())
        V = np.array([[complex(x, y) for x, y in row] for row in d["rows"]])
    else:
        V = load_A()
        V /= np.linalg.norm(V)
    n = V.shape[1] - 1
    assert n == a.n0, f"seed has n={n}"
    best = None
    while n < a.n1:
        t0 = time.time()
        for _ in range(a.stride):
            cands = []
            for phi in [(1.0, 0.0), (0.0, 1.0)]:
                W = grow(V, phi)
                th = threshold(W)
                cands.append((th if th else 0, W))
            V = max(cands, key=lambda c: c[0])[1]
            n += 1
        V, th = climb(V, a.rounds, iters=a.iters, verbose=False)
        f = save(V, n, th, a.tag)
        best = th if best is None else max(best, th)
        print(f"n={n}: threshold {th:.8f}  ({time.time() - t0:.0f}s)  {f.name}", flush=True)


if __name__ == "__main__":
    main()
