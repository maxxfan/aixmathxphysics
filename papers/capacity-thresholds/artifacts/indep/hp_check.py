#!/usr/bin/env python3
"""High-precision re-run of indep_ic's pipeline (Python decimal, 45 digits) to measure the float64
error.  Everything is recomputed in decimal: Lindblad rates, Clebsch-Gordan factors,
multiplicity ratios, Poisson weights, the uniformisation series and (for selected blocks) the
eigenvalues (mpmath eigh at 40 digits).  Plain frame (no rotation), so the input doubles are
used exactly.

    hp_check.py FILE [P] [--blocks k1,k2,...]
"""
import argparse
import math
import os
import sys
import time
from decimal import Decimal as Dm, getcontext
from fractions import Fraction

import mpmath
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import indep_ic as I  # noqa: E402

getcontext().prec = 45


def D(fr):
    fr = Fraction(fr)
    return Dm(fr.numerator) / Dm(fr.denominator)


def hp_rates(probs):
    pI, pX, pY, pZ = probs
    lx = pI + pX - pY - pZ
    ly = pI - pX + pY - pZ
    lz = pI - pX - pY + pZ
    Lx, Ly, Lz = (-(D(l).ln()) / 2 for l in (lx, ly, lz))
    g = [(Ly + Lz - Lx) / 2, (Lx + Lz - Ly) / 2, (Lx + Ly - Lz) / 2]
    big = max(abs(x) for x in g)
    return [Dm(0) if abs(x) < big * Dm('1e-30') else x for x in g]


def cg(j2pp, j2, m2, s2):
    den = 2 * (j2pp + 1)
    if j2 == j2pp + 1:
        num = j2pp + s2 * m2 + 1
        return (Dm(max(num, 0)) / den).sqrt()
    num = j2pp - s2 * m2 + 1
    v = (Dm(max(num, 0)) / den).sqrt()
    return -v if s2 == 1 else v


IDX = {1: 0, -1: 1}
PA = {
    'X': [[0, 1], [1, 0]],
    'Y': [[0, -1j], [1j, 0]],
    'Z': [[1, 0], [0, -1]],
}


def hp_superop(L, P):
    """dict (row, col) -> Decimal for rho -> sum_i P_i rho P_i on weighted Schur blocks."""
    n = L.n
    A = B = PA[P]
    out = {}
    for kp, j2p in enumerate(L.j2):
        dp = j2p + 1
        for j2pp in (j2p - 1, j2p + 1):
            if j2pp < 0 or j2pp > n - 1:
                continue
            mpp = I.mult(n - 1, j2pp)
            for j2in in (j2pp - 1, j2pp + 1):
                if j2in < 0 or j2in > n:
                    continue
                kin = (n - j2in) // 2
                din = j2in + 1
                F = D(Fraction(n * mpp, I.mult(n, j2in)))
                for s2 in (1, -1):
                    for t2 in (1, -1):
                        Ast = A[IDX[s2]][IDX[t2]]
                        if Ast == 0:
                            continue
                        for s2p in (1, -1):
                            for t2p in (1, -1):
                                Bts = B[IDX[t2p]][IDX[s2p]]
                                if Bts == 0:
                                    continue
                                ph = complex(Ast) * complex(Bts)
                                assert ph.imag == 0
                                ph = int(ph.real)
                                er = (j2in - j2p + s2 - t2) // 2
                                ec = (j2in - j2p + s2p - t2p) // 2
                                rowf = []
                                for a in range(dp):
                                    m2 = j2p - 2 * a
                                    ai = a + er
                                    if ai < 0 or ai > j2in:
                                        rowf.append(None)
                                        continue
                                    rowf.append(cg(j2pp, j2p, m2, s2) * cg(j2pp, j2in, m2 - s2 + t2, t2))
                                colf = []
                                for b in range(dp):
                                    m2 = j2p - 2 * b
                                    bi = b + ec
                                    if bi < 0 or bi > j2in:
                                        colf.append(None)
                                        continue
                                    colf.append(cg(j2pp, j2in, m2 - s2p + t2p, t2p) * cg(j2pp, j2p, m2, s2p))
                                for a in range(dp):
                                    ra = rowf[a]
                                    if ra is None or ra == 0:
                                        continue
                                    fa = F * ra * ph
                                    row0 = L.off[kp] + a * dp
                                    col0 = L.off[kin] + (a + er) * din
                                    for b in range(dp):
                                        cb = colf[b]
                                        if cb is None or cb == 0:
                                            continue
                                        key = (row0 + b, col0 + b + ec)
                                        out[key] = out.get(key, Dm(0)) + fa * cb
    return out


def to_rows(N, coef_dicts):
    """sum of (scalar, dict) -> list of rows [(col, coef), ...]."""
    acc = {}
    for s, dct in coef_dicts:
        for key, v in dct.items():
            acc[key] = acc.get(key, Dm(0)) + s * v
    rows = [[] for _ in range(N)]
    for (r, c), v in acc.items():
        if v != 0:
            rows[r].append((c, v))
    return rows


def matvec(rows, X):
    """X: list of C columns, each a list of N Decimals."""
    out = []
    for xc in X:
        out.append([sum((v * xc[c] for c, v in row), Dm(0)) for row in rows])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('p', nargs='?', default=None)
    ap.add_argument('--blocks', default='')
    args = ap.parse_args()
    target, pstr, V, note = I.load_candidate(args.file)
    if args.p:
        pstr = args.p
    p = Fraction(pstr)
    r, n1 = V.shape
    n = n1 - 1
    L = I.Layout(n)
    probs = I.pauli_probs(target, p)
    t0 = time.time()
    g = hp_rates(probs)
    names = ['X', 'Y', 'Z']
    S = {}
    for nm, gg in zip(names, g):
        if gg != 0:
            S[nm] = hp_superop(L, nm)
    t1 = time.time()
    print(f'{target} n={n} r={r} p={pstr} rates {[float(x) for x in g]}  superops {t1 - t0:.0f}s '
          f'nnz {sum(len(s) for s in S.values())}', flush=True)
    # input (exact doubles), normalised in decimal
    Zs = sum(Dm(float(x.real)) ** 2 + Dm(float(x.imag)) ** 2 for x in V.ravel())
    sq = Zs.sqrt()
    vre = [[Dm(float(V[q, a].real)) / sq for a in range(n1)] for q in range(r)]
    vim = [[Dm(float(V[q, a].imag)) / sq for a in range(n1)] for q in range(r)]
    N = L.N
    X = []  # columns: for comp c=(q,q'), real then imag; |v_q><v_q'| entries v_q[a] conj(v_q'[b])
    for q in range(r):
        for qq in range(r):
            re = [Dm(0)] * N
            im = [Dm(0)] * N
            d = L.d[0]
            for a in range(d):
                for b in range(d):
                    xr, xi = vre[q][a], vim[q][a]
                    yr, yi = vre[qq][b], -vim[qq][b]
                    re[a * d + b] = xr * yr - xi * yi
                    im[a * d + b] = xr * yi + xi * yr
            X.append(re)
            X.append(im)
    pos = [(nm, gg) for nm, gg in zip(names, g) if gg > 0]
    neg = [(nm, gg) for nm, gg in zip(names, g) if gg < 0]
    if neg:
        shift = sum(gg * n for nm, gg in neg)
        rows = to_rows(N, [(gg, S[nm]) for nm, gg in neg])
        for i in range(N):
            rows[i].append((i, -shift))
        bound = sum(-2 * gg * n for nm, gg in neg)
        acc = [list(c) for c in X]
        T = X
        k = 0
        while True:
            k += 1
            T = matvec(rows, T)
            T = [[x / k for x in c] for c in T]
            for ci in range(len(acc)):
                acc[ci] = [a + b for a, b in zip(acc[ci], T[ci])]
            mx = max(abs(x) for c in T for x in c)
            if k > bound and mx < Dm('1e-44'):
                break
        X = acc
        print(f'  negative-rate Taylor: {k} terms  {time.time() - t1:.0f}s', flush=True)
    Gam = sum(gg for nm, gg in pos)
    mu = n * Gam
    rows = to_rows(N, [(gg / (n * Gam), S[nm]) for nm, gg in pos])
    w = (-mu).exp()
    acc = [[w * x for x in c] for c in X]
    T = X
    k = 0
    cum = w
    while True:
        k += 1
        w = w * mu / k
        T = matvec(rows, T)
        for ci in range(len(acc)):
            acc[ci] = [a + w * b for a, b in zip(acc[ci], T[ci])]
        cum += w
        if k > mu and 1 - cum < Dm('1e-42'):
            break
    t2 = time.time()
    print(f'  uniformisation: mu={float(mu):.4f} {k + 1} terms  {t2 - t1:.0f}s', flush=True)
    # float64 pipeline for comparison
    res = I.evaluate(target, p, V, method='uniform', frame='plain')
    Xf = res['X']
    # assemble hp blocks, compare
    C = r * r
    hpX = np.zeros((N, C), dtype=complex)
    for c in range(C):
        hpX[:, c] = np.array([float(x) for x in acc[2 * c]]) + 1j * np.array([float(x) for x in acc[2 * c + 1]])
    worst_rel = 0.0
    for k_, d in enumerate(L.d):
        sl = slice(L.off[k_], L.off[k_ + 1])
        tr = sum(abs(hpX[sl, q * r + q].reshape(d, d).trace()) for q in range(r))
        dif = np.max(np.abs(Xf[sl] - hpX[sl]))
        worst_rel = max(worst_rel, dif / max(tr, 1e-300))
    ic_f64, _ = I.block_entropies(L, Xf, r)
    ic_hpround, _ = I.block_entropies(L, hpX, r)
    print(f'  float64 I_c            {ic_f64:.15e}')
    print(f'  hp blocks -> f64 eig   {ic_hpround:.15e}   diff {ic_f64 - ic_hpround:+.2e}')
    print(f'  max |X_f64 - X_hp| / tr(block) over blocks: {worst_rel:.2e}')
    # high-precision eigenvalues for selected blocks
    blocks = [int(x) for x in args.blocks.split(',') if x != '']
    mpmath.mp.dps = 40
    for kb in blocks:
        d = L.d[kb]
        off = L.off[kb]
        # tau[(a,q),(b,q')] = X[(a,b),(q,q')]
        tau = mpmath.matrix(d * r, d * r)
        sig = mpmath.matrix(d, d)
        for a in range(d):
            for b in range(d):
                for q in range(r):
                    for qq in range(r):
                        c = q * r + qq
                        val = mpmath.mpc(mpmath.mpf(str(acc[2 * c][off + a * d + b])),
                                         mpmath.mpf(str(acc[2 * c + 1][off + a * d + b])))
                        tau[a * r + q, b * r + qq] = val
                        if q == qq:
                            sig[a, b] += val
        # hermitise (the series is Hermitian to ~1e-45)
        tau = (tau + tau.H) / 2
        sig = (sig + sig.H) / 2
        te = time.time()
        es = mpmath.eigh(sig, eigvals_only=True)
        et = mpmath.eigh(tau, eigvals_only=True)

        def Hmp(ev):
            s = mpmath.mpf(0)
            for x in ev:
                if x > 0:
                    s -= x * mpmath.log(x, 2)
            return s
        hs, ht = Hmp(es), Hmp(et)
        # float64 entropies of the float64 block
        blk = Xf[off:off + d * d].reshape(d, d, r, r)
        tauf = blk.transpose(0, 2, 1, 3).reshape(d * r, d * r)
        sigf = sum(blk[:, :, q, q] for q in range(r))
        hsf = I._H(np.linalg.eigvalsh((sigf + sigf.conj().T) / 2))
        htf = I._H(np.linalg.eigvalsh((tauf + tauf.conj().T) / 2))
        print(f'  block k={kb} (2j={L.j2[kb]}, dim {d}/{d * r}): weight {float(mpmath.re(sum(sig[i, i] for i in range(d)))):.4e}  '
              f'min eig sigma {mpmath.nstr(min(es), 3)} tau {mpmath.nstr(min(et), 3)}  '
              f'dH_sigma(f64-hp) {float(hsf - hs):+.2e}  dH_tau {float(htf - ht):+.2e}  '
              f'd(Hs-Ht) {float((hsf - htf) - (hs - ht)):+.2e}  ({time.time() - te:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
