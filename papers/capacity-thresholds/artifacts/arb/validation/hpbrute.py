"""High-precision full 2^n density-matrix reference for I_c (validation of arb_ic.py only).

mpmath at `dps` digits: builds the Dicke vectors explicitly, applies the Pauli channel qubit by
qubit as the mixture of Pauli conjugations (on the system of rho and of the purification), and
takes eigenvalues with mpmath's own Hermitian eigensolver.  Shares no code with arb_ic.py's
channel (no Schur blocks, no generator) or with its spectral certification (no Gershgorin).
Its error is ~10^-(dps-5), far below arb_ic's enclosure widths, so containment is a real test.
"""
from fractions import Fraction
from math import comb

import mpmath as mp
import numpy as np


def _h(ev):
    s = mp.mpf(0)
    for x in ev:
        x = mp.re(x)
        if x > 0:
            s -= x * mp.log(x)
    return s / mp.log(2)


def _channel(T, n, r, probs, cplx):
    dim = (1 << n) * r
    idx = np.arange(dim)
    x = idx // r
    q = idx % r
    pI, pX, pY, pZ = [mp.mpf(Fraction(v).numerator) / Fraction(v).denominator for v in probs]
    I_ = mp.mpc(0, 1)
    for k in range(n):
        bit = (x >> k) & 1
        perm = ((x ^ (1 << k)) * r + q)
        Tf = T[np.ix_(perm, perm)]
        out = T * pI
        if pX != 0:
            out = out + Tf * pX
        if pY != 0:
            f = np.array([(-I_ if b == 0 else I_) for b in bit], dtype=object)
            fc = np.array([mp.conj(v) for v in f], dtype=object)
            out = out + (f[:, None] * Tf * fc[None, :]) * pY
        if pZ != 0:
            z = np.array([1 if b == 0 else -1 for b in bit], dtype=object)
            out = out + (z[:, None] * T * z[None, :]) * pZ
        T = out
    return T


def ic_hp(rows, probs, dps=50):
    """rows: list of rows of (re, im) python floats (Dicke coefficients); probs: Fractions."""
    mp.mp.dps = dps
    r = len(rows)
    n = len(rows[0]) - 1
    N = 1 << n
    cplx = any(im != 0 for row in rows for (_, im) in row)
    Z = sum(Fraction(re) ** 2 + Fraction(im) ** 2 for row in rows for (re, im) in row)
    Zm = mp.mpf(Z.numerator) / Z.denominator
    w = [bin(x).count('1') for x in range(N)]
    vec = np.empty((r, N), dtype=object)
    for qq in range(r):
        for x in range(N):
            re, im = rows[qq][w[x]]
            vec[qq, x] = mp.mpc(re, im) / mp.sqrt(comb(n, w[x]))
    sig = np.empty((N, N), dtype=object)
    for a in range(N):
        for b in range(N):
            s = mp.mpc(0)
            for qq in range(r):
                s += vec[qq, a] * mp.conj(vec[qq, b])
            sig[a, b] = s / Zm
    psi = np.empty(N * r, dtype=object)
    for x in range(N):
        for qq in range(r):
            psi[x * r + qq] = vec[qq, x]
    tau = np.empty((N * r, N * r), dtype=object)
    for a in range(N * r):
        for b in range(N * r):
            tau[a, b] = psi[a] * mp.conj(psi[b]) / Zm
    sig = _channel(sig, n, 1, probs, cplx)
    tau = _channel(tau, n, r, probs, cplx)
    out = []
    for M in (sig, tau):
        d = M.shape[0]
        Mm = mp.matrix(d, d)
        for i in range(d):
            for j in range(d):
                Mm[i, j] = M[i, j]
        if cplx:
            ev = mp.eighe(Mm, eigvals_only=True)
        else:
            Mr = mp.matrix(d, d)
            for i in range(d):
                for j in range(d):
                    Mr[i, j] = mp.re(Mm[i, j])
            ev = mp.eigsy(Mr, eigvals_only=True)
        out.append(_h(ev))
    return out[0] - out[1]
