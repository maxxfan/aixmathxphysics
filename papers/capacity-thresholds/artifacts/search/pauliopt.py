#!/usr/bin/env python3
"""Gradient optimisation of permutation-symmetric rank-r inputs for any qubit Pauli channel
(search tool; certification is verifier/verify.py).  Same block formula as qlib/icfloat.py
(path 1), vectorised, with the exact gradient:

  block_k(Y) = m_k D_d^{-1/2} DP_d( corr^k( D_n^{-1/2} Y D_n^{-1/2} ) ) D_d^{-1/2}
  DP_d is self-adjoint (the coefficient tensor of (x^dag Z y)^d is symmetric under swapping
  input and output monomials), corr^dag = conv, so
  L_k^dag(W) = m_k D_n^{-1/2} conv^k( DP_d( D_d^{-1/2} W D_d^{-1/2} ) ) D_n^{-1/2}.

usage:
  pauliopt.py check                                  finite-difference gradient check
  pauliopt.py climb FAMILY SEED.json                 optimise a seed (rows) at its n
  pauliopt.py grow FAMILY N0 N1 SEED.json            grow n -> n+1 and re-optimise
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import brentq, minimize

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import depopt  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent / "states"
LN2 = math.log(2)
FAM = {
    "depolarizing": lambda p: (1 - 3 * p, p, p, p),
    "two_pauli": lambda p: (1 - 2 * p, p, 0.0, p),
    "independent_xz": lambda p: ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p)),
    # the same channels conjugated by exp(-i pi/4 X) per qubit (Y <-> Z): p_X = p_Y, so the channel
    # is z-covariant and the optimal input is real in the Dicke basis (yframe.py maps states across)
    "two_pauli@y": lambda p: (1 - 2 * p, p, p, 0.0),
    "independent_xz@y": lambda p: ((1 - p) ** 2, p * (1 - p), p * (1 - p), p * p),
}
RANGE = {"depolarizing": (0.06, 0.0685), "two_pauli": (0.10, 0.125), "independent_xz": (0.10, 0.125),
         "two_pauli@y": (0.10, 0.125), "independent_xz@y": (0.10, 0.125)}
DA = [(0, 0), (1, 1), (2, 2), (0, 2), (2, 0)]


class Pauli:
    def __init__(self, n, probs):
        self.n = n
        pI, pX, pY, pZ = probs
        cA, cB, cC, cE = pI + pZ, pX + pY, pI - pZ, pX - pY
        self.det = [cA * cB, cA * cA + cB * cB - cC * cC - cE * cE, cA * cB, -cC * cE, -cC * cE]
        # (da, db) -> [(dm, dm', w)]
        self.kern = [((0, 0), [(0, 0, cA), (1, 1, cB)]), ((1, 1), [(0, 0, cB), (1, 1, cA)]),
                     ((0, 1), [(0, 1, cC), (1, 0, cE)]), ((1, 0), [(1, 0, cC), (0, 1, cE)])]
        self.isq = [1 / np.sqrt(np.array([math.comb(d, a) for a in range(d + 1)], dtype=float)) for d in range(n + 1)]
        self.mk = [math.comb(n, k) - (math.comb(n, k - 1) if k else 0) for k in range(n // 2 + 1)]

    def corr(self, Y):
        m = Y.shape[-1] - 2
        out = np.zeros(Y.shape[:-2] + (m, m), dtype=Y.dtype)
        for (al, be), w in zip(DA, self.det):
            if w:
                out += w * Y[..., al:al + m, be:be + m]
        return out

    def conv(self, W):
        m = W.shape[-1]
        out = np.zeros(W.shape[:-2] + (m + 2, m + 2), dtype=W.dtype)
        for (al, be), w in zip(DA, self.det):
            if w:
                out[..., al:al + m, be:be + m] += w * W
        return out

    def dp(self, X):
        """sum_{a,b} H^{(d)}[a,b,m,m'] X[a,b], X of shape (..., d+1, d+1)"""
        d = X.shape[-1] - 1
        lead = X.shape[:-2]
        # level s = d: G[al, be, 0, 0] = X[al, be]
        G = X.reshape(lead + (d + 1, d + 1, 1, 1))
        for s in range(d - 1, -1, -1):
            L = d - s  # current inner size; new inner size L+1
            Gn = np.zeros(lead + (s + 1, s + 1, L + 1, L + 1), dtype=X.dtype)
            for (da, db), outs in self.kern:
                src = G[..., da:da + s + 1, db:db + s + 1, :, :]
                for dm, dm2, w in outs:
                    if w:
                        Gn[..., dm:dm + L, dm2:dm2 + L] += w * src
            G = Gn
        return G[..., 0, 0, :, :]

    def blocks(self, Y):
        n = self.n
        Yh = Y * np.multiply.outer(self.isq[n], self.isq[n])
        out = []
        for k in range(n // 2 + 1):
            if k:
                Yh = self.corr(Yh)
            d = n - 2 * k
            out.append(self.mk[k] * self.dp(Yh) * np.multiply.outer(self.isq[d], self.isq[d]))
        return out

    def adjoint(self, k, W):
        n, d = self.n, self.n - 2 * k
        U = self.dp(W * np.multiply.outer(self.isq[d], self.isq[d]))
        for _ in range(k):
            U = self.conv(U)
        return self.mk[k] * U * np.multiply.outer(self.isq[n], self.isq[n])

    def value_grad(self, V, grad=True):
        r = V.shape[0]
        Y = np.einsum("ia,jb->ijab", V, V.conj())
        Z = float(np.real(np.einsum("iiaa->", Y)))
        f, G = 0.0, (np.zeros_like(Y) if grad else None)
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
            if grad:
                ls = (Us * np.log(np.maximum(ws, 1e-300))) @ Us.conj().T
                lt = (Ut * np.log(np.maximum(wt, 1e-300))) @ Ut.conj().T
                W = lt.reshape(r, d1, r, d1).transpose(0, 2, 1, 3).copy()
                for q in range(r):
                    W[q, q] -= ls
                G += self.adjoint(k, W)
        I = f / Z
        if not grad:
            return I, None
        g = np.einsum("ijab,jb->ia", G, V)
        return I, (2 * g - 2 * I * V) / Z


def model(fam, n, p):
    if fam == "depolarizing":
        return depopt.Depol(n, p)
    return Pauli(n, FAM[fam](p))


def ic_bits(fam, V, p):
    return model(fam, V.shape[1] - 1, p).value_grad(V, grad=False)[0] / LN2


def threshold(fam, V, margin=0.0):
    lo, hi = RANGE[fam]
    f = lambda p: ic_bits(fam, V, p) - margin
    if f(lo) <= 0:
        return None
    while f(hi) > 0:
        hi += 0.003
    return brentq(f, lo, hi, xtol=1e-10)


def optimise(fam, V, p, iters=400):
    r, n1 = V.shape
    m = model(fam, n1 - 1, p)

    def fun(x):
        I, g = m.value_grad(depopt.unpack(x, r, n1))
        return -I, -depopt.pack(g)

    res = minimize(fun, depopt.pack(V), jac=True, method="L-BFGS-B",
                   options={"maxiter": iters, "gtol": 1e-14, "ftol": 1e-16})
    W = depopt.unpack(res.x, r, n1)
    return W / np.linalg.norm(W)


def climb(fam, V, rounds=6, iters=400, step=None, verbose=True):
    step = step or (2e-5 if fam == "depolarizing" else 1e-4)
    th = threshold(fam, V)
    for _ in range(rounds):
        W = optimise(fam, V, th + step, iters)
        t2 = threshold(fam, W)
        if verbose:
            print(f"  {fam} threshold {th:.8f} -> {t2}", flush=True)
        if t2 is None or t2 <= th:
            break
        V, th, gain = W, t2, t2 - th
        if gain < 1e-9:
            break
    return V, th


def load_rows(path):
    d = json.loads(pathlib.Path(path).read_text())
    if "rows" in d:
        return np.array([[complex(x, y) for x, y in row] for row in d["rows"]])
    import base64
    raw = base64.b64decode(d["payload_c128_b64"])
    return np.frombuffer(raw, dtype="<c16").reshape(d.get("rank", 2), -1).copy()


def save(fam, V, th, tag):
    OUT.mkdir(exist_ok=True)
    n = V.shape[1] - 1
    f = OUT / f"{fam}_n{n}_{tag}.json"
    f.write_text(json.dumps({"target": fam, "n": n, "threshold_float": th,
                             "rows": [[[z.real, z.imag] for z in row] for row in V]}))
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["check", "climb", "grow"])
    ap.add_argument("args", nargs="*")
    ap.add_argument("--rounds", type=int, default=6)
    a = ap.parse_args()
    if a.mode == "check":
        rng = np.random.default_rng(1)
        V = rng.normal(size=(2, 8)) + 1j * rng.normal(size=(2, 8))
        for fam in ("two_pauli", "independent_xz", "depolarizing"):
            m = Pauli(7, FAM[fam](0.09))
            I, g = m.value_grad(V)
            I2 = depopt.Depol(7, 0.09).value_grad(V, False)[0] if fam == "depolarizing" else None
            x, gx = depopt.pack(V), depopt.pack(g)
            errs = []
            for i in rng.choice(len(x), 5, replace=False):
                e = np.zeros_like(x)
                e[i] = 1e-6
                fd = (m.value_grad(depopt.unpack(x + e, 2, 8), False)[0]
                      - m.value_grad(depopt.unpack(x - e, 2, 8), False)[0]) / 2e-6
                errs.append(abs(fd - gx[i]))
            print(fam, "max |fd - grad| =", max(errs), "" if I2 is None else f"(dep paths differ by {abs(I - I2):.1e})")
        return
    fam = a.args[0]
    if a.mode == "climb":
        V = load_rows(a.args[1])
        V /= np.linalg.norm(V)
        print(f"start threshold {threshold(fam, V)}")
        V, th = climb(fam, V, a.rounds)
        print(f"{fam} n={V.shape[1] - 1}: {th:.8f} -> {save(fam, V, th, 'climb')}")
        return
    n0, n1, seed = int(a.args[1]), int(a.args[2]), a.args[3]
    V = load_rows(seed)
    V /= np.linalg.norm(V)
    assert V.shape[1] - 1 == n0
    n = n0
    while n < n1:
        cands = []
        for phi in [(1.0, 0.0), (0.0, 1.0), (2 ** -0.5, 2 ** -0.5)]:
            W = depopt.grow(V, phi)
            th = threshold(fam, W)
            cands.append((th or 0.0, W))
        V = max(cands, key=lambda c: c[0])[1]
        n += 1
        V, th = climb(fam, V, a.rounds, verbose=False)
        print(f"{fam} n={n}: threshold {th:.8f}  {save(fam, V, th, 'grow').name}", flush=True)


if __name__ == "__main__":
    main()
