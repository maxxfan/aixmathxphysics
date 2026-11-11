"""Turn a float FGGS chain into an exact certificate for verifier/verify.py.

1. Exact rounding that keeps every root condition: q_1..q_{k-2} are rounded step by step and each
   increment is projected (in integers) onto the subspace allowed at that step; q_{k-1} is then
   *solved* from q_{k-2} and q_k = 1 (it is determined by them), all over D = 2 n 2^b.
2. Nonnegativity: for each q_t, gamma = margin/2 and h = minimum-phase spectral factor of
   q_t - gamma (cepstrum, then Wilson's Newton iteration; where that stalls, cepstral steps from
   h = 1 with spike-corrected cepstral coefficients and a power-series exponential), rounded to
   2^-bh.  The verifier needs the residual q_t - |h|^2 to be diagonally dominant.
"""
import numpy as np

import fggs
from interp import Interp, grid_cos, grid_sin


def exact_chain(cs, n, b=56):
    """cs = float coefficient vectors c_1..c_{k-1} (length >= n).  Returns (D, [C_1..C_{k-1}])."""
    k = len(cs) + 1
    D0 = n << b
    prev = [D0 - j * (D0 // n) for j in range(n)]
    out = []
    for t in range(1, k - 1):
        s = 1 if t % 2 == 0 else -1
        target = [int(round(float(x) * D0)) for x in cs[t - 1][:n]]
        d = [target[j] - prev[j] for j in range(n)]
        d[0] = 0
        for j in range(1, n // 2 + 1):
            jj = n - j
            if j == jj:
                if s == 1:
                    d[j] = 0
            elif s == -1:   # d_j = d_{n-j}
                m = (d[j] + d[jj]) // 2
                d[j] = d[jj] = m
            else:           # d_j = -d_{n-j}
                m = (d[j] - d[jj]) // 2
                d[j], d[jj] = m, -m
        cur = [prev[j] + d[j] for j in range(n)]
        out.append(cur)
        prev = cur
    # q_{k-1}: x - prev in V_{k-1} (sign s1), x - 1 in V_k (sign -s1); doubled to stay integral
    s1 = 1 if (k - 1) % 2 == 0 else -1
    out = [[2 * v for v in C] for C in out]
    prev = [2 * v for v in prev]
    D = 2 * D0
    x = [0] * n
    x[0] = D
    for j in range(1, n // 2 + 1):
        jj = n - j
        if j == jj:
            x[j] = prev[j] if s1 == 1 else 0
            continue
        alpha = prev[j] + s1 * prev[jj]        # x_j + s1 x_{n-j} = alpha and x_j = s1 x_{n-j}
        assert alpha % 2 == 0
        x[j] = alpha // 2
        x[jj] = s1 * (alpha // 2)
    out.append(x)
    return D, out


def eval_grid(c, M):
    """q(theta_r) = c0 + 2 sum_j c_j cos(j theta_r), theta_r = 2 pi r / M, r = 0..M/2 (q is even,
    so the half grid is the whole grid; all transforms below are real-symmetric half spectra)"""
    x = np.zeros(M)
    x[: len(c)] = c
    x[1: len(c)] *= 2
    return np.fft.rfft(x).real


def autocorr(h):
    m = len(h)
    L = 1 << int(np.ceil(np.log2(2 * m)))
    f = np.fft.rfft(h, L)
    return np.fft.irfft(f * np.conj(f), L)[:m]


def causal(X, M):
    """{X}_+ for a real even X given on the half grid (length M/2 + 1): zero lag halved, positive
    lags kept.  Returns the half spectrum (the full one is conjugate-symmetric)."""
    g = np.fft.irfft(X, M)
    g[0] *= 0.5
    g[M // 2:] = 0.0
    return np.fft.rfft(g)


def spectral_factor(c, gamma, tol, M0=None, Mmax=1 << 26, wilson=12, verbose=False):
    """minimum-phase h (length n) with |h|^2 ~ q - gamma to max coefficient error <= tol.

    Cepstrum on an M-point grid, then Wilson's Newton step H <- H {P/|H|^2 + 1}_+ with h truncated
    to degree < n after each step (the step is exact Newton once the grid resolves 1/|H|^2).  M is
    doubled until the error is below tol.  Returns (h, err, M)."""
    n = len(c)
    target = np.asarray(c, float).copy()
    target[0] -= gamma
    M = M0 or 1 << int(np.ceil(np.log2(64 * n)))
    best = None
    while M <= Mmax:
        p = eval_grid(c, M) - gamma
        if p.min() <= 0:
            raise ValueError(f"q - gamma not positive on the grid (min {p.min():.3e})")
        Hf = np.exp(causal(np.log(p), M))
        prev = np.inf
        for _ in range(wilson):
            h = np.fft.irfft(Hf, M)[:n]
            err = float(np.max(np.abs(autocorr(h) - target)))
            if best is None or err < best[1]:
                best = (h, err, M)
            if err <= tol or err > 0.9 * prev:
                break
            prev = err
            Hf = np.fft.rfft(h, M)
            Hf = Hf * causal(p / np.abs(Hf) ** 2 + 1.0, M)
        if verbose:
            print(f"    M=2^{int(np.log2(M))}: err {best[1]:.2e} (tol {tol:.1e})", flush=True)
        if best[1] <= tol:
            break
        del p, Hf
        M *= 2
    return best


def chain_coeffs(us, n):
    """float coefficient vectors c_1..c_{k-1} (length n) from half-vectors u_0..u_{k-1}"""
    out = []
    for v in fggs.chain_values(us, n):
        c = fggs.coeffs_from_values(v)
        out.append(c[:n])
    return out


def _chi(x, a, b):
    """1 for |x| <= a, 0 for |x| >= b, C-infinity in between"""
    t = np.clip((np.abs(x) - a) / (b - a), 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        e0 = np.where(t > 0, np.exp(-1.0 / np.maximum(t, 1e-300)), 0.0)
        e1 = np.where(t < 1, np.exp(-1.0 / np.maximum(1.0 - t, 1e-300)), 0.0)
    return e1 / (e0 + e1)


def find_spikes(c, gamma, M, Rf=16, W=24, dmax=3.5, edge=32):
    """local minima of q (coefficients c) where the roots of q - gamma come within dmax grid
    spacings (2 pi / M) of the unit circle: 1 - |root| ~ d = sqrt(2 (q0 - gamma) / q'').  There
    1/|h|^2 is a spike too narrow for the M-point grid.  Returns (theta, d), theta in (0, pi)."""
    n = len(c)
    N = n * Rf
    hs = np.pi / N
    j = np.arange(n)
    v = grid_cos(c, N)
    g1 = -grid_sin(j * c, N)                 # q'
    g2 = grid_cos(-(j * j) * c, N)           # q''
    i = np.arange(1, N)
    cand = i[(v[i] <= v[i - 1]) & (v[i] <= v[i + 1])]
    Delta = 2 * np.pi / M
    dcrude = np.sqrt(2 * np.maximum(v[cand] - gamma, 0) / np.maximum(g2[cand], 1e-300))
    cand = cand[dcrude < dmax * Delta + 2 * hs]
    th = cand * hs
    I = Interp(n, Rf, W)
    for _ in range(6):
        q1 = I.matrix(th, odd=True) @ g1
        q2 = I.matrix(th) @ g2
        ok = q2 > 0
        th = np.where(ok, np.clip(th - np.clip(q1 / np.where(ok, q2, 1.0), -hs, hs), 0, np.pi), th)
    A = I.matrix(th)
    q0, q2 = A @ v, A @ g2
    d = np.sqrt(2 * np.maximum(q0 - gamma, 1e-300) / np.maximum(q2, 1e-300))
    keep = (d < dmax * Delta) & (th > edge * Delta) & (th < np.pi - edge * Delta)
    return th[keep], d[keep]


def _graded_nodes(ds, a, npan, npt):
    """per-spike Gauss nodes and weights on [-a, a]: a panel [0, d/2], then npan - 1 geometric
    panels from d/2 to a, mirrored.  Shape (S, 2 npan npt)."""
    gx, gw = np.polynomial.legendre.leggauss(npt)
    lo = np.minimum(ds / 2, a / 4)
    ratio = (a / lo) ** (1.0 / (npan - 1))
    e = np.concatenate([np.zeros((len(ds), 1)), lo[:, None] * ratio[:, None] ** np.arange(npan)], axis=1)
    e[:, -1] = a
    mid, half = (e[:, 1:] + e[:, :-1]) / 2, (e[:, 1:] - e[:, :-1]) / 2
    x = (mid[:, :, None] + half[:, :, None] * gx).reshape(len(ds), -1)
    w = (half[:, :, None] * gw).reshape(len(ds), -1)
    return np.concatenate([-x[:, ::-1], x], axis=1), np.concatenate([w[:, ::-1], w], axis=1)


def spike_coeffs(hf, r, n, M, Fg, th_s, d_s, Rf=16, W=24, P=26, a_mult=8, b_mult=32, ncore=96,
                 npanel=6, nper=12, chunk=2000, core="tan", npan=10, npt=10, fun=None):
    """Correction to the M-point trapezoid coefficients F_k (k < n) of F = r/|H|^2 from the
    spikes at th_s: with a smooth window chi_s (1 on |x| <= a, 0 beyond b), the integral of
    F chi_s e^{-ik theta} is done by a Gauss rule on the core and Gauss panels on the
    transition, the grid's own share of F chi_s is subtracted, and both are folded into all k at
    once through moments in n (theta - th_s) and the type-1 sums of `Interp`.
    core="tan": the map x = d tan(phi), exact for a Lorentzian of width d; core="graded":
    geometric panels from d/2 to a, for integrands with features of any width in between.  With
    fun, the integrand is fun(r/|H|^2) and Fg must hold fun of the grid values."""
    Delta = 2 * np.pi / M
    a, b = a_mult * Delta, b_mult * Delta
    N = n * Rf
    reH = (grid_cos(hf, N) + hf[0]) / 2
    imH = grid_sin(hf, N) / 2
    rg = grid_cos(r, N)
    I = Interp(n, Rf, W)
    gx, gw = np.polynomial.legendre.leggauss(ncore)
    px, pw = np.polynomial.legendre.leggauss(nper)
    e = np.linspace(a, b, npanel + 1)
    tx = ((e[:-1, None] + e[1:, None]) / 2 + (e[1:, None] - e[:-1, None]) / 2 * px).ravel()
    tw = ((e[1:, None] - e[:-1, None]) / 2 * pw).ravel()
    tx, tw = np.r_[-tx[::-1], tx], np.r_[tw[::-1], tw]
    tw = tw * _chi(tx, a, b)
    L = int(np.ceil(2 * b / Delta)) + 2
    kn = np.arange(n) / n
    out = np.zeros(n)
    for c0 in range(0, len(th_s), chunk):
        ts, ds = th_s[c0:c0 + chunk], d_s[c0:c0 + chunk]
        S = len(ts)
        if core == "tan":
            phi0 = np.arctan(a / ds)
            phi = phi0[:, None] * gx
            cx, cw = ds[:, None] * np.tan(phi), phi0[:, None] * gw * ds[:, None] / np.cos(phi) ** 2
        else:
            cx, cw = _graded_nodes(ds, a, npan, npt)
        X = np.concatenate([cx, np.broadcast_to(tx, (S, len(tx)))], axis=1)
        Wt = np.concatenate([cw, np.broadcast_to(tw, (S, len(tx)))], axis=1)
        flat = (ts[:, None] + X).ravel()
        Ae, Ao = I.matrix(flat), I.matrix(flat, odd=True)
        Fv = (Ae @ rg) / ((Ae @ reH) ** 2 + (Ao @ imH) ** 2)
        if fun is not None:
            Fv = fun(Fv)
        Y = Wt * Fv.reshape(S, -1)
        mm = np.ceil((ts - b) / Delta).astype(np.int64)[:, None] + np.arange(L)
        xm = mm * Delta - ts[:, None]
        Ym = -Delta * Fg[np.clip(mm, 0, M // 2)] * _chi(xm, a, b)
        tN, tM = n * X, n * xm
        pN, pM = np.ones_like(tN), np.ones_like(tM)
        mom = np.empty((P, S))
        for p in range(P):
            mom[p] = (Y * pN).sum(1) + (Ym * pM).sum(1)
            pN *= tN
            pM *= tM
        Ac, As = I.matrix(ts), I.matrix(ts, odd=True)
        kp, fact = np.ones(n), 1.0
        for p in range(P):
            if p % 2 == 0:
                term = I.cos_sums(ts, mom[p], Ac) * (1 if p % 4 == 0 else -1)
            else:
                term = I.sin_sums(ts, mom[p], As) * (-1 if p % 4 == 1 else 1)
            out += kp / fact * term
            kp = kp * kn
            fact *= p + 1
    return out / np.pi


def exact_residual(Ct, D, G, hi, H):
    """r_j = (q - gamma - |h|^2)_j for j < n as floats, from the exact integers (FLINT product):
    q = Ct / D, gamma = G / D, h = hi / H"""
    from flint import fmpz_poly
    n, m = len(Ct), len(hi)
    sq = fmpz_poly(hi) * fmpz_poly(hi[::-1])
    DH2 = D * H * H
    return np.array([((Ct[j] - (G if j == 0 else 0)) * H * H - D * (int(sq[m - 1 + j]) if j < m else 0)) / DH2
                     for j in range(n)])


def _smul(a, b, m):
    """(a b) mod z^m for real power series"""
    a, b = a[:m], b[:m]
    L = 1 << int(np.ceil(np.log2(len(a) + len(b))))
    return np.fft.irfft(np.fft.rfft(a, L) * np.fft.rfft(b, L), L)[:m]


def _sinv(a, m):
    """1/a mod z^m (Newton: b <- b (2 - a b))"""
    b, L = np.array([1.0 / a[0]]), 1
    while L < m:
        L = min(2 * L, m)
        e = -_smul(a, b, L)
        e[0] += 2.0
        b = _smul(b, e, L)
    return b


def _slog(e, m):
    """log e mod z^m for e_0 = 1"""
    q = _smul(np.arange(1, m) * e[1:m], _sinv(e, m - 1), m - 1)
    return np.r_[0.0, q / np.arange(1, m)]


def _sexp(a, m):
    """exp a mod z^m (Newton: e <- e (1 + a - log e))"""
    e, L = np.array([np.exp(a[0])]), 1
    while L < m:
        L = min(2 * L, m)
        ee = np.r_[e, np.zeros(L - len(e))]
        f = a[:L] - _slog(ee / ee[0], L)
        f[0] = a[0] - np.log(ee[0]) + 1.0
        e = _smul(ee, f, L)
    return e


def cepstral_step(Ct, D, gamma, hi, H, M, spikes=None):
    """Multiplicative step h <- h exp({G}_+), G = log((q - gamma)/|h|^2) = log1p(r/|h|^2).

    For minimum-phase h, log h* - log h = {G}_+, so one exact step gives the factor h* of q - gamma;
    and exp of a power series mod z^n needs only its first n terms, so only G_0..G_{n-1} are
    needed: the M-point grid plus `spike_coeffs` corrections (graded rule, fun = log1p).  Wilson's
    step 1 + {r/|h|^2}_+ is its linearisation; it cannot move a root from ~50 d to d from the
    circle in one step, and at (7, 850000) it stalled where h lacked such near-zeros.  Returns the
    new integer factor."""
    n = len(Ct)
    hi = list(hi) + [0] * (n - len(hi))
    r = exact_residual(Ct, D, int(round(gamma * D)), hi, H)
    hf = np.array(hi, dtype=float) / H
    Hf = np.fft.rfft(hf, M)
    Gg = np.log1p(eval_grid(r, M) / np.abs(Hf) ** 2)
    del Hf
    Gk = np.fft.irfft(Gg, M)[:n]
    if spikes is not None and len(spikes[0]):
        Gk += spike_coeffs(hf, r, n, M, Gg, *spikes, core="graded", fun=np.log1p, chunk=1000)
    del Gg
    Gk[0] *= 0.5                                        # {.}_+
    E = _sexp(Gk, n)
    E[0] -= 1.0
    delta = _smul(hf, E, n)                             # h (exp{G}_+ - 1), added in integers
    return [v + int(round(d * H)) for v, d in zip(hi, delta)]


def cepstral_factor(Ct, D, gamma, H, M, tol, spikes=None, steps=5, verbose=False):
    """Factor of q - gamma by cepstral steps from h = 1 (minimum phase), while each step gains at
    least 10x.  Returns the integer h with the smallest residual and that max |r_j|."""
    n = len(Ct)
    G = int(round(gamma * D))
    hi = [H] + [0] * (n - 1)
    best = (hi, float(np.max(np.abs(exact_residual(Ct, D, G, hi, H)))))
    for s in range(steps):
        hi = cepstral_step(Ct, D, gamma, hi, H, M, spikes)
        err = float(np.max(np.abs(exact_residual(Ct, D, G, hi, H))))
        if verbose:
            print(f"    cepstral step {s + 1}: max residual {err:.2e} (tol {tol:.1e})", flush=True)
        prev = best[1]
        if err < prev:
            best = (hi, err)
        if err <= tol or err > 0.1 * prev:
            break
    return best


def refine(Ct, D, gamma, h, H, M, tol, iters=6, spikes=None, verbose=False):
    """Mixed-precision Newton refinement of an integer spectral factor.

    The residual r = (q - gamma) - |h|^2 is computed exactly (FLINT products of the integer h);
    the correction is Wilson's Newton step d = {r/|h|^2}_+ h evaluated in floating point.  Since r
    is small, float errors in the step are relative to r, so the residual keeps shrinking below
    the float64 floor of |h|^2.  The Fourier coefficients of r/|h|^2 come from an M-point grid,
    plus `spike_coeffs` corrections at `spikes` = (theta, d) where 1/|h|^2 is narrower than the
    grid.  Stops early once the residual has grown twice in a row (the step diverges at spikes
    the grid cannot resolve; at (7, 775000) it then kept doubling).  Returns the integer h with
    the smallest residual seen, and that max |r_j|."""
    n = len(Ct)
    G = int(round(gamma * D))
    hi = [int(round(x * H)) for x in h]
    best, prev, rises = None, np.inf, 0
    for it in range(iters + 1):
        r = exact_residual(Ct, D, G, hi, H)
        err = float(np.max(np.abs(r)))
        if verbose:
            print(f"    refine {it}: max residual {err:.2e} (tol {tol:.1e})", flush=True)
        if best is None or err < best[1]:
            best = (hi, err)
        rises = rises + 1 if err > prev else 0
        prev = err
        if err <= tol or it == iters or rises == 2:
            break
        hf = np.array(hi, dtype=float) / H
        Hf = np.fft.rfft(hf, M)
        Fg = eval_grid(r, M) / np.abs(Hf) ** 2
        del Hf
        Fk = np.fft.irfft(Fg, M)[:n]                    # trapezoid coefficients of r/|h|^2
        if spikes is not None and len(spikes[0]):
            Fk += spike_coeffs(hf, r, n, M, Fg, *spikes)
        del Fg
        Fk[0] *= 0.5                                    # {.}_+
        Lc = 1 << int(np.ceil(np.log2(2 * n)))
        delta = np.fft.irfft(np.fft.rfft(hf, Lc) * np.fft.rfft(Fk, Lc), Lc)[:n]
        hi = [v + int(round(d * H)) for v, d in zip(hi, delta)] + [0] * (n - len(hi))
    return best


def make_certificate(cs, n, b=56, bh=60, gamma_frac=0.5, meta=None, iters=20, Mretry=1 << 27,
                     verbose=False):
    """cs: float coefficient vectors c_1..c_{k-1}.  Returns (certificate dict, diagnostics).
    If refinement does not converge, the factor is recomputed on a grid twice as fine (up to
    Mretry) and refined again: at (7, 700000) the first Newton step from the 2^26 cepstrum of q_6
    (error 1.2e-5) overshot at the spikes, while from 2^27 (error 2.6e-11) it converged in 2.
    At (7, 775000) q_6 needed 9 steps at 2^27 (linear phase), hence iters=20.  If that still
    fails, the factor is rebuilt from h = 1 by `cepstral_factor` on the finest grid and refined:
    at (7, 850000) q_6 (closest root 2.8e-9) stalled at 1.3e-7 under refinement at 2^27, and two
    cepstral steps reached 5.7e-17."""
    k = len(cs) + 1
    D, C = exact_chain(cs, n, b)
    H = 1 << bh
    hs, diag = [], []
    for t, Ct in enumerate(C, start=1):
        c = np.array([v / D for v in Ct])     # exact -> float (rel. error 1e-16)
        qmin = float(eval_grid(c, 1 << int(np.ceil(np.log2(64 * n)))).min())
        gamma = gamma_frac * qmin
        budget = gamma / (2 * n)              # diagonal dominance needs sum_j |r_j| <~ gamma / 2
        if verbose:
            print(f"  q_{t}: grid min {qmin:.3e}, gamma {gamma:.3e}, budget {budget:.1e}", flush=True)
        h, err, M = spectral_factor(c, gamma, tol=0.05 * budget, verbose=verbose)
        hi = None
        while err > 0.05 * budget:
            spikes = find_spikes(c, gamma, M)
            if verbose:
                print(f"    {len(spikes[0])} spikes narrower than the grid "
                      f"(closest root distance {spikes[1].min() if len(spikes[1]) else np.inf:.2e})",
                      flush=True)
            hi, err = refine(Ct, D, gamma, h, H, M, tol=0.05 * budget, iters=iters, spikes=spikes,
                             verbose=verbose)
            if err <= 0.05 * budget or 2 * M > Mretry:
                break
            h, err, M = spectral_factor(c, gamma, tol=0.05 * budget, M0=2 * M, Mmax=2 * M,
                                        verbose=verbose)
            hi = None
        if err > 0.05 * budget:
            spikes = find_spikes(c, gamma, M)
            hi, err = cepstral_factor(Ct, D, gamma, H, M, 0.05 * budget, spikes, verbose=verbose)
            if err > 0.05 * budget:
                hi, err = refine(Ct, D, gamma, np.array(hi, dtype=float) / H, H, M, tol=0.05 * budget,
                                 iters=iters, spikes=spikes, verbose=verbose)
        if hi is None:
            hi = [int(round(x * H)) for x in h]
        hs.append(hi)
        diag.append({"t": t, "grid_min": qmin, "gamma": gamma, "factor_err": err,
                     "budget": budget, "M": M})
    cert = {"type": "fggs_chain", "k": k, "n": n, "den": D, "q": C, "h_den": H, "h": hs}
    if meta:
        cert["meta"] = meta
    return cert, diag
