"""rigeig: rigorous enclosure of  sum_i h(lambda_i),  h(x) = -x log2 x,  over the eigenvalues of a
Hermitian positive semidefinite matrix known only as an Arb ball matrix H (plus a spectral-norm
perturbation of size <= slack).

Method
------
1. Hint: float64 eigendecomposition of mid(H) (LAPACK via numpy).  Not trusted.
2. Refinement (heuristic, any error only widens the result): Newton-Schulz re-orthonormalisation
   X <- X (3I - X*X)/2, first-order Jacobi corrections X <- X (I + E), E_ij = s_ij/(l_j - l_i)
   for separated pairs of S = X* A X, and re-diagonalisation of clusters of close Rayleigh
   quotients (float eigh of the shifted and rescaled cluster block, which resolves the next
   scale of a multi-scale spectrum).  All products in Arb at the working precision.
3. Certification (rigorous):
   * B = X* H X and G = X* X - I in Arb ball arithmetic (X exact), eta >= ||G||_2 (max row sum).
   * Gershgorin on B: discs [b_ii -+ sum_{j!=i} |b_ij|]; a connected component of k discs holds
     exactly k eigenvalues of every Hermitian matrix in the ball B.
   * Ostrowski: lambda_k(X* A X) = theta_k lambda_k(A), theta_k in [1-eta, 1+eta]; so the
     eigenvalues of A with the indices of a component lie in the component interval divided by
     theta.  Weyl: the slack widens each by +-slack.  PSD: lower ends clipped at 0 (an interval
     entirely below 0 raises, it would contradict positivity).
   * h is concave with maximum 1/(e ln 2) at 1/e: per component k * [min h(ends), max h(I)],
     evaluated with Arb's rigorous log.
A block whose certified trace is below `tiny` is bounded without any eigenvalues:
sum h(l_i) in [h(t), t log2(d/t)] (subadditivity and concavity of h, t <= 1).
"""
import math
import time

import numpy as np
from flint import acb, acb_mat, arb, arb_mat, ctx

_INV_E = None
_LN2 = None


def _consts():
    global _INV_E, _LN2
    if _INV_E is None or _INV_E.rel_accuracy_bits() < ctx.prec - 8:
        _INV_E = (-arb(1)).exp()
        _LN2 = arb(2).log()
    return _INV_E, _LN2


def h_exact(x):
    """h(x) = -x log2 x for an exact arb x >= 0 (ball result)."""
    inv_e, ln2 = _consts()
    if x == 0:
        return arb(0)
    return -x * x.log() / ln2


def h_bounds(a, b):
    """rigorous [lower, upper] (arb exact endpoints) of h over [a, b], 0 <= a <= b exact arbs."""
    inv_e, ln2 = _consts()
    ha = h_exact(a)
    hb = h_exact(b)
    lo = ha.lower() if ha.lower() < hb.lower() else hb.lower()
    if b < inv_e:
        hi = hb.upper()
    elif a > inv_e:
        hi = ha.upper()
    else:
        hi = (inv_e / ln2).upper()
    return lo, hi


def _MT(cplx):
    return acb_mat if cplx else arb_mat


def to_np(M, cplx):
    d = M.nrows()
    e = M.entries()
    if cplx:
        return np.array([complex(x) for x in e], dtype=np.complex128).reshape(d, M.ncols())
    return np.array([float(x) for x in e], dtype=np.float64).reshape(d, M.ncols())


def from_np(a, cplx):
    return _MT(cplx)(a.tolist())


def ct(X, cplx):
    return X.conjugate().transpose() if cplx else X.transpose()


def eye(d, cplx):
    M = _MT(cplx)(d, d)
    for i in range(d):
        M[i, i] = 1
    return M


def _diag_re(M, cplx, i):
    x = M[i, i]
    return x.real if cplx else x


def _clusters(lam_hi, rad):
    """connected components of overlapping intervals [l_i - rad_i, l_i + rad_i]."""
    d = len(lam_hi)
    order = np.argsort(lam_hi)
    cl = np.empty(d, dtype=np.int64)
    cid = 0
    cur_hi = None
    for idx in order:
        lo_i = lam_hi[idx] - rad[idx]
        if cur_hi is not None and lo_i <= cur_hi:
            cur_hi = max(cur_hi, lam_hi[idx] + rad[idx])
        else:
            cid += 1
            cur_hi = lam_hi[idx] + rad[idx]
        cl[idx] = cid
    return cl


def refine(Am, cplx, maxit=10, verbose=False):
    """heuristic high-precision eigenvectors of the exact Hermitian matrix Am (mid matrix).

    Each iteration: Newton-Schulz re-orthonormalisation X <- X (3I - X*X)/2 until X*X = I to
    working precision, S = X* A X, first-order Jacobi correction
    X <- X (I + E), E_ij = s_ij/(l_j - l_i) for pairs in different clusters, and float
    re-diagonalisation of every cluster block of S (shifted by its mean and rescaled), which
    resolves the next scale of a multi-scale spectrum.  Only a hint: certify() trusts nothing
    here."""
    d = Am.nrows()
    prec = ctx.prec
    A_np = to_np(Am, cplx)
    A_np = (A_np + A_np.conj().T) / 2
    w, Q = np.linalg.eigh(A_np)
    normA = float(np.max(np.abs(w))) if d else 0.0
    info = dict(iters=0, offmax=[], clusters=[], ns=0)
    if d == 1 or normA == 0.0:
        return from_np(Q, cplx), info
    floor = normA * d * 2.0 ** (-prec + 6)
    ofloor = d * 2.0 ** (-prec + 6)
    I = eye(d, cplx)
    X = from_np(Q, cplx)
    for it in range(maxit):
        # Newton-Schulz until X*X = I to working precision (quadratic: a defect below 1e-20
        # is squared below the floor by one more step)
        for _ in range(8):
            G = (ct(X, cplx) * X).mid()
            defect = float(np.abs(to_np(G - I, cplx)).max())
            if defect <= ofloor:
                break
            X = ((X * (I * 3 - G)) * arb(0.5)).mid()
            info['ns'] += 1
            if defect < 1e-20:
                break
        S = (ct(X, cplx) * (Am * X).mid()).mid()
        S_np = to_np(S, cplx)
        lam = [_diag_re(S, cplx, i) for i in range(d)]
        lam_hi = np.array([float(x) for x in lam])
        lam_lo = np.array([float(x - arb(float(x))) for x in lam])
        F = S_np.copy()
        np.fill_diagonal(F, 0)
        absF = np.abs(F)
        offmax = float(absF.max())
        info['offmax'].append(offmax)
        info['iters'] = it + 1
        if offmax <= floor:
            break
        rho = np.sqrt((absF ** 2).sum(axis=1))
        cl = _clusters(lam_hi, 4 * rho + floor)
        diff = (lam_hi[None, :] - lam_hi[:, None]) + (lam_lo[None, :] - lam_lo[:, None])  # l_j - l_i
        same = cl[:, None] == cl[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            E = np.where(same, 0.0, F / np.where(same, 1.0, diff))
        X = (X + (X * from_np(E, cplx)).mid()).mid()
        # cluster re-diagonalisation (on the pre-update S; the error is second order)
        V = np.eye(d, dtype=np.complex128 if cplx else np.float64)
        nrot = 0
        sizes = []
        for c in np.unique(cl):
            idx = np.nonzero(cl == c)[0]
            if len(idx) < 2:
                continue
            sizes.append(len(idx))
            mu = sum((lam[i] for i in idx), arb(0)) / len(idx)
            sub = S_np[np.ix_(idx, idx)].copy()
            for a_, i in enumerate(idx):
                sub[a_, a_] = float(lam[i] - mu)
            scale = float(np.max(np.abs(sub)))
            if scale <= floor:
                continue
            sub = (sub + sub.conj().T) / (2 * scale)
            _, Vc = np.linalg.eigh(sub)
            V[np.ix_(idx, idx)] = Vc
            nrot += 1
        info['clusters'].append((len(sizes), max(sizes) if sizes else 0, nrot))
        if nrot:
            X = (X * from_np(V, cplx)).mid()
        if verbose:
            print(f'    it {it}: offmax {offmax:.2e} floor {floor:.1e} clusters {info["clusters"][-1]}')
    else:
        # iteration budget exhausted after an update: re-orthonormalise before certification
        for _ in range(8):
            G = (ct(X, cplx) * X).mid()
            defect = float(np.abs(to_np(G - I, cplx)).max())
            if defect <= ofloor:
                break
            X = ((X * (I * 3 - G)) * arb(0.5)).mid()
            info['ns'] += 1
        info['maxit_reached'] = True
    return X, info


def certify(H, X, cplx, slack):
    """rigorous (lo, hi) of sum h(lambda_i) for every Hermitian A in the ball H, perturbed by
    any Hermitian E with ||E||_2 <= slack (arb, exact), assuming A + E is PSD."""
    d = H.nrows()
    Xh = ct(X, cplx)
    B = Xh * H * X
    G = Xh * X - eye(d, cplx)
    eb = B.entries()
    eg = G.entries()
    # eta = ||G||_2 <= max row abs sum (G Hermitian)
    eta = arb(0)
    for i in range(d):
        s = arb(0)
        for j in range(d):
            s += abs(eg[i * d + j]).upper()
        su = s.upper()
        if su > eta:
            eta = su
    if not eta < 0.25:
        raise RuntimeError('eigenvector basis too far from orthonormal: eta = %s' % eta)
    ivals = []
    for i in range(d):
        s = arb(0)
        for j in range(d):
            if j != i:
                s += abs(eb[i * d + j]).upper()
        c = eb[i * d + i].real if cplx else eb[i * d + i]
        ivals.append(((c - s).lower(), (c + s).upper()))
    ivals.sort(key=lambda t: t[0])
    comps = []
    for lo, hi in ivals:
        if comps and lo <= comps[-1][1]:
            c = comps[-1]
            comps[-1] = [c[0], hi if hi > c[1] else c[1], c[2] + 1]
        else:
            comps.append([lo, hi, 1])
    S_lo = arb(0)
    S_hi = arb(0)
    one_m = 1 - eta
    one_p = 1 + eta
    nneg = 0
    worst_width = 0.0
    for lo, hi, k in comps:
        a = (lo / one_m) if lo < 0 else (lo / one_p)
        a = (a - slack).lower()
        b = (hi / one_p) if hi < 0 else (hi / one_m)
        b = (b + slack).upper()
        if b < 0:
            raise RuntimeError('VERIFIER_ALARM: eigenvalue interval entirely negative [%s, %s]' % (a, b))
        if a < 0:
            nneg += k
            a = arb(0)
        hl, hh = h_bounds(a, b)
        S_lo += hl * k
        S_hi += hh * k
        worst_width = max(worst_width, float(b - a))
    info = dict(d=d, eta=float(eta), ncomp=len(comps), nclip=nneg,
                maxcomp=max(c[2] for c in comps), worst_interval=worst_width)
    return S_lo.lower(), S_hi.upper(), info


def trace_bounds(H, cplx, slack):
    d = H.nrows()
    t = arb(0)
    for i in range(d):
        x = H[i, i]
        t += x.real if cplx else x
    return (t - slack * d).lower(), (t + slack * d).upper()


def entropy_enclosure(H, cplx, slack, tiny=1e-30, verbose=False):
    """rigorous (lo, hi) of sum_i h(lambda_i) [bits] for the PSD matrix in ball H (+- slack)."""
    t0 = time.time()
    d = H.nrows()
    tlo, thi = trace_bounds(H, cplx, slack)
    if thi < tiny:
        # sum h in [h(t), t log2(d/t)] for 0 <= t <= 1
        if thi <= 0:
            return arb(0), arb(0), dict(d=d, mode='zero', time=time.time() - t0)
        inv_e, ln2 = _consts()
        hi = (thi * (arb(d) / thi).log() / ln2).upper()
        return arb(0), hi, dict(d=d, mode='tiny-trace', trace_hi=float(thi), time=time.time() - t0)
    X, rinfo = refine(H.mid(), cplx, verbose=verbose)
    lo, hi, cinfo = certify(H, X, cplx, slack)
    cinfo.update(mode='eig', refine_iters=rinfo['iters'],
                 offmax=rinfo['offmax'][-1] if rinfo['offmax'] else 0.0,
                 time=time.time() - t0)
    return lo, hi, cinfo
