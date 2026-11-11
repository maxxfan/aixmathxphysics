#!/usr/bin/env python3
"""kernel_ic: a second evaluation of the coherent information I_c of a p08 input.  NOT RIGOROUS.

Purpose: an evaluation of the inputs certified by results/rigorous/arb_ic.py (release: arb/)
that shares no code and no channel representation with it.  arb_ic.py integrates a master
equation on the Schur blocks (Clebsch-Gordan superoperator, Poisson series) and certifies spectra
with Gershgorin discs.  This file uses the generating function of the blocks (the "kernel"
det(Z)^k Sym^{2j}(Z)), built by a recurrence in the spin with positive terms only, and takes
eigenvalues from LAPACK (float mode) or mpmath (mp mode).  It imports numpy and mpmath only, and
nothing from verifier/, results/rigorous/ or results/indep/.

It returns a number, not an enclosure.  In float mode the error is that of double precision
(observed: 1e-15 to 1e-13 bits, see validate_kernel.py and run_headline.py); in mp mode it is that
of the chosen number of digits.  Use it as corroboration of a certified sign, never as a proof.

Usage:
    kernel_ic.py FILE [--p P] [--mp DIGITS]

Mathematics (derived for this file)
-----------------------------------
Input: rows v_q in Sym^n(C^2) in the Dicke basis |D_a> (a = number of ones), rho = sum_q
|v_q><v_q| / Z, purification sum_q |v_q>|q> / sqrt Z.  The output-reference state is
sum_{q,q'} N^{(x)n}(|v_q><v_q'|) (x) |q><q'| / Z.

1. z-covariant frame.  A Pauli channel with p_X = p_Y acts on matrix units as
       |0><0| -> (1-f)|0><0| + f|1><1|,   |1><1| -> f|0><0| + (1-f)|1><1|,
       |0><1| -> c |0><1|,                |1><0| -> c |1><0|,
   with f = p_X + p_Y and c = p_I - p_Z.  The depolarizing channel has this form.  For the
   independent X-Z and two-Pauli channels, W = (1 - iX)/sqrt 2 maps Y -> Z and Z -> -Y, so
   I_c(rho, N) = I_c(W rho W^+, N') with N' the Pauli channel with p_Y and p_Z exchanged, which
   has p_X = p_Y.  On Sym^n,
       <D_b| W^{(x)n} |D_a> = 2^{-n/2} sqrt(C(n,b)/C(n,a)) (-i)^{b-a} K(a,b),
       K(a,b) = sum_j (-1)^j C(n-b,j) C(b,a-j)      (an exact integer),
   from  W(|0> + t|1>) = ((1 - it)|0> + (t - i)|1>)/sqrt 2.

2. Generating function.  (|0> + x|1>)^{(x)n} = sum_a x^a sqrt(C(n,a)) |D_a>, so
       N^{(x)n}(|D_a><D_b|) = [x^a y^b] Z^{(x)n} / sqrt(C(n,a) C(n,b)),
       Z = N((|0> + x|1>)(<0| + y<1|)) = [[(1-f) + f s, c y], [c x, f + (1-f) s]],  s = x y.
   On the Schur sector with spin j = n/2 - k (d = 2j, multiplicity m_k = C(n,k) - C(n,k-1)),
   Z^{(x)n} acts as det(Z)^k Sym^d(Z) (x) 1, and
       sqrt(C(d,al) C(d,be)) Sym^d(Z)[al,be] = [u^al t^be] (Z00 + Z01 t + Z10 u + Z11 u t)^d.
   For al - be = D >= 0 the right side is x^D P_d[al,be](s) with a polynomial P of degree
   <= d - D, and  det Z = e0 + e1 s + e0 s^2,  e0 = f(1-f),  e1 = (1-f)^2 + f^2 - c^2 >= 0.

3. Recurrence (all coefficients >= 0).  With A = (1-f) + f s and B = f + (1-f) s,
       P_d[al,be] = A P_{d-1}[al,be] + c s P_{d-1}[al,be-1] + c P_{d-1}[al-1,be] + B P_{d-1}[al-1,be-1]
   (for D = 0 the third term is c s P_{d-1}[al,al-1], by the symmetry al <-> be), P_0[0,0] = 1.

4. Blocks.  The weighted block (trace = probability of the sector) of N^{(x)n}(X), X = sum
   X_ab |D_a><D_b|, is, for al - be = D >= 0,
       B_k[al,be] = m_k / sqrt(C(d,al) C(d,be)) * sum_m P_d[al,be]_m u^{(k,D)}_m,
       u^{(k,D)}_m = sum_b X_{b+D,b} / sqrt(C(n,b+D) C(n,b)) * [s^{b-m}] det(Z)^k,
   and u^{(k+1)}_m = e0 u^{(k)}_m + e1 u^{(k)}_{m+1} + e0 u^{(k)}_{m+2}.  Entries with al < be
   follow from Hermiticity.  tau_k is assembled from the r^2 operators |v_q><v_q'|, sigma_k is
   its partial trace over q, and I_c = sum_k [H(sigma_k) - H(tau_k)] in bits.

Scaling: P is carried divided by 2^d and det by det(1) = 1 - c^2, so nothing overflows for
n <= 1000; the factors are restored exactly in the block prefactor.
"""
import argparse
import base64
import json
import math
import sys
import time
from fractions import Fraction

import numpy as np

try:
    import mpmath
except Exception:                                    # float mode does not need mpmath
    mpmath = None


# ------------------------------------------------------------------------------ input

def _num(x):
    if isinstance(x, bool):
        raise ValueError('boolean is not a number')
    if isinstance(x, str):
        t = x.strip()
        return float.fromhex(t) if t.lower().lstrip('+-').startswith('0x') else float(t)
    return float(x)


def load(path):
    """-> (target, p as string, rows: list of lists of python complex; the exact doubles)."""
    with open(path) as fh:
        d = json.load(fh)
    if 'rows' in d:
        rows = []
        for row in d['rows']:
            rr = []
            for z in row:
                if isinstance(z, (list, tuple)):
                    rr.append(complex(_num(z[0]), _num(z[1])))
                else:
                    rr.append(complex(_num(z), 0.0))
            rows.append(rr)
    else:
        raw = np.frombuffer(base64.b64decode(d['payload_c128_b64']), dtype='<c16')
        rows = [[complex(z) for z in row] for row in raw.reshape(int(d['rank']), -1)]
    if len({len(r) for r in rows}) != 1:
        raise ValueError('rows of unequal length')
    return d['target'], d['p'], rows


def pauli_probs(target, p):
    p = Fraction(p)
    if target == 'depolarizing':
        return (1 - 3 * p, p, p, p)
    if target == 'independent_xz':
        return ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p))
    if target == 'two_pauli':
        return (1 - 2 * p, p, Fraction(0), p)
    raise ValueError('unknown target %r' % (target,))


# ------------------------------------------------------------------------------ arithmetic

class Float:
    name = 'float64'

    def real(self, q):
        q = Fraction(q)
        return q.numerator / q.denominator

    def sqrt_int(self, m):
        return math.sqrt(m)

    def zeros(self, shape, cplx=False):
        return np.zeros(shape, dtype=np.complex128 if cplx else np.float64)

    def amp(self, z):
        return complex(z)

    def entropy(self, M):
        """sum -x log2 x over the eigenvalues of the Hermitian matrix M."""
        if M.shape[0] == 0:
            return 0.0
        if not M.imag.any():
            M = M.real
        w = np.linalg.eigvalsh((M + M.conj().T) / 2)
        w = w[w > 0]
        return float(-(w * np.log2(w)).sum())

    def trace(self, M):
        return float(np.trace(M).real)

    def conj(self, a):
        return np.conj(a)

    def out(self, x):
        return float(x)


class MP:
    def __init__(self, dps):
        if mpmath is None:
            raise RuntimeError('mp mode needs mpmath')
        mpmath.mp.dps = dps
        self.name = 'mpmath-%d' % dps
        self.zero = mpmath.mpf(0)

    def real(self, q):
        q = Fraction(q)
        return mpmath.mpf(q.numerator) / q.denominator

    def sqrt_int(self, m):
        return mpmath.sqrt(mpmath.mpf(m))

    def zeros(self, shape, cplx=False):
        a = np.empty(shape, dtype=object)
        a.fill(mpmath.mpc(0) if cplx else mpmath.mpf(0))
        return a

    def amp(self, z):
        return mpmath.mpc(z.real, z.imag)             # doubles are converted exactly

    def entropy(self, M):
        d = M.shape[0]
        if d == 0:
            return mpmath.mpf(0)
        A = mpmath.matrix(d, d)
        for i in range(d):
            for j in range(d):
                A[i, j] = (M[i, j] + mpmath.conj(M[j, i])) / 2
        ev = mpmath.eigh(A, eigvals_only=True)
        s = mpmath.mpf(0)
        for x in ev:
            x = mpmath.re(x)
            if x > 0:
                s -= x * mpmath.log(x)
        return s / mpmath.log(2)

    def trace(self, M):
        return sum((mpmath.re(M[i, i]) for i in range(M.shape[0])), mpmath.mpf(0))

    def conj(self, a):
        out = np.empty(a.shape, dtype=object)
        for i, x in enumerate(a):
            out[i] = mpmath.conj(x)
        return out

    def out(self, x):
        return x


# ------------------------------------------------------------------------------ frame change

def krawtchouk(n):
    """K[a][b] = sum_j (-1)^j C(n-b,j) C(b,a-j), exact integers."""
    C = [[math.comb(m, j) for j in range(m + 1)] for m in range(n + 1)]
    K = [[0] * (n + 1) for _ in range(n + 1)]
    for b in range(n + 1):
        Cnb, Cb = C[n - b], C[b]
        for a in range(n + 1):
            t = 0
            for j in range(max(0, a - b), min(a, n - b) + 1):
                v = Cnb[j] * Cb[a - j]
                t += -v if j & 1 else v
            K[a][b] = t
    return K


def rotate_rows(rows, ar):
    """W^{(x)n} v for every row, W = (1 - iX)/sqrt 2."""
    n = len(rows[0]) - 1
    K = krawtchouk(n)
    binom = [math.comb(n, a) for a in range(n + 1)]
    root = [ar.sqrt_int(b) for b in binom]
    if isinstance(ar, Float):
        scale = lambda kint, b, a: float(kint) * 2.0 ** (-n / 2) * root[b] / root[a]                 # noqa: E731
        unit = [1, -1j, -1, 1j]
    else:
        two = mpmath.mpf(2) ** (-mpmath.mpf(n) / 2)
        scale = lambda kint, b, a: mpmath.mpf(kint) * two * root[b] / root[a]                       # noqa: E731
        unit = [mpmath.mpc(1), mpmath.mpc(0, -1), mpmath.mpc(-1), mpmath.mpc(0, 1)]
    out = []
    for row in rows:
        amps = [ar.amp(z) for z in row]
        new = []
        for b in range(n + 1):
            s = 0
            for a in range(n + 1):
                if row[a] != 0:
                    s = s + unit[(b - a) % 4] * scale(K[a][b], b, a) * amps[a]
            new.append(s)
        out.append(new)
    return out


# ------------------------------------------------------------------------------ evaluation

def evaluate(target, p, rows, mp_digits=None, verbose=False):
    """I_c in bits (float, or mpf in mp mode) and a dict of diagnostics."""
    t0 = time.time()
    ar = MP(mp_digits) if mp_digits else Float()
    n = len(rows[0]) - 1
    r = len(rows)
    pI, pX, pY, pZ = pauli_probs(target, p)
    if target == 'depolarizing':
        V = [[ar.amp(z) for z in row] for row in rows]
    else:
        pY, pZ = pZ, pY
        V = rotate_rows(rows, ar)
    if pX != pY:
        raise ValueError('channel is not z-covariant in this frame')
    fq, cq = pX + pY, pI - pZ
    e0q = fq * (1 - fq)
    e1q = (1 - fq) ** 2 + fq ** 2 - cq ** 2
    det1q = 2 * e0q + e1q                             # det Z at s = 1, = 1 - c^2
    assert e1q >= 0 and det1q > 0
    half = Fraction(1, 2)
    a0, a1, cc = ar.real((1 - fq) * half), ar.real(fq * half), ar.real(cq * half)
    g0, g1 = ar.real(e0q / det1q), ar.real(e1q / det1q)
    # exact normalisation Z = sum |v|^2 of the given doubles (the frame change is unitary)
    Zq = Fraction(0)
    for row in rows:
        for z in row:
            Zq += Fraction(z.real) ** 2 + Fraction(z.imag) ** 2
    Zinv = ar.real(1 / Zq)
    binom_n = [math.comb(n, a) for a in range(n + 1)]
    rootn = [ar.sqrt_int(b) for b in binom_n]
    pairs = [(q, q2) for q in range(r) for q2 in range(r)]
    npair = len(pairs)
    cplx = True
    nk = n // 2 + 1
    # u[k][D]: array (d - D + 1, npair), d = n - 2k
    conj = (lambda z: z.conjugate()) if isinstance(ar, Float) else mpmath.conj
    u_prev = []
    for D in range(n + 1):
        w = ar.zeros((n - D + 1, npair), cplx)
        for ip, (q, q2) in enumerate(pairs):
            for b in range(n - D + 1):
                w[b, ip] = V[q][b + D] * conj(V[q2][b]) * Zinv / (rootn[b + D] * rootn[b])
        u_prev.append(w)
    U = [u_prev]
    for k in range(1, nk):
        d = n - 2 * k
        cur = []
        for D in range(d + 1):
            L = d - D + 1
            w = u_prev[D]
            cur.append(g0 * w[:L] + g1 * w[1:L + 1] + g0 * w[2:L + 2])
        U.append(cur)
        u_prev = cur
    t_u = time.time() - t0
    # P tables, old/new buffers; P[D] has shape (n + 1 - D, n + 1 - D), indices [be, m]
    old = [ar.zeros((n + 1 - D, n + 1 - D)) for D in range(n + 1)]
    new = [ar.zeros((n + 1 - D, n + 1 - D)) for D in range(n + 1)]
    old[0][0, 0] = ar.real(1)
    Ic = ar.real(0)
    total_trace = ar.real(0)
    Hs = Ht = ar.real(0)

    def emit(d, P):
        nonlocal Ic, total_trace, Hs, Ht
        k = (n - d) // 2
        mk = math.comb(n, k) - (math.comb(n, k - 1) if k else 0)
        # prefactor m_k 2^d det(1)^k / sqrt(C(d,al) C(d,be)), kept as row and column factors
        if isinstance(ar, Float):
            scal = float(mk) * 2.0 ** d * float(det1q) ** k
        else:
            scal = mpmath.mpf(mk) * mpmath.mpf(2) ** d * ar.real(det1q) ** k
        rd = np.array([1 / ar.sqrt_int(math.comb(d, a)) for a in range(d + 1)],
                      dtype=np.float64 if isinstance(ar, Float) else object)
        dim = d + 1
        tau = ar.zeros((r * dim, r * dim), cplx)
        for D in range(d + 1):
            L = d - D + 1
            Y = P[D][:L, :L].dot(U[k][D][:L])          # (be, pair)
            be = np.arange(L)
            al = be + D
            pref = scal * rd[al] * rd[be]
            for ip, (q, q2) in enumerate(pairs):
                vals = Y[:, ip] * pref
                tau[q * dim + al, q2 * dim + be] = vals
                if D:
                    tau[q2 * dim + be, q * dim + al] = ar.conj(vals)
        sig = ar.zeros((dim, dim), cplx)
        for q in range(r):
            sig = sig + tau[q * dim:(q + 1) * dim, q * dim:(q + 1) * dim]
        hs, ht = ar.entropy(sig), ar.entropy(tau)
        Hs, Ht = Hs + hs, Ht + ht
        Ic = Ic + (hs - ht)
        total_trace = total_trace + ar.trace(sig)
        if verbose:
            print(f'  block k={k} d={d}: tr {ar.out(ar.trace(sig))}  H(sigma)-H(tau) {ar.out(hs - ht)}',
                  flush=True)

    if n % 2 == 0:
        emit(0, old)
    for d in range(1, n + 1):
        for D in range(d + 1):
            s = d - D + 1
            N = new[D]
            O = old[D]
            N[:s, :s] = a0 * O[:s, :s]
            N[:s, 1:s] += a1 * O[:s, :s - 1]
            N[1:s, :s] += a1 * O[:s - 1, :s]
            N[1:s, 1:s] += a0 * O[:s - 1, :s - 1]
            if D + 1 <= d - 1 and s > 1:
                N[1:s, 1:s] += (2 * cc if D == 0 else cc) * old[D + 1][:s - 1, :s - 1]
            if D >= 1:
                N[:s, :s] += cc * old[D - 1][:s, :s]
        old, new = new, old
        if (n - d) % 2 == 0:
            emit(d, old)
    info = dict(n=n, rank=r, mode=ar.name, f=float(fq), c=float(cq), trace=ar.out(total_trace),
                H_out=ar.out(Hs), H_joint=ar.out(Ht), seconds=time.time() - t0, seconds_u=t_u)
    return ar.out(Ic), info


def fmt(x, digits=None):
    if isinstance(x, float):
        return '%.12e' % x
    return mpmath.nstr(x, digits or 25)


def main(argv=None):
    ap = argparse.ArgumentParser(description='second, non-rigorous evaluation of I_c (kernel formula)')
    ap.add_argument('file')
    ap.add_argument('--p', default=None, help='exact p overriding the file ("a/b" or decimal)')
    ap.add_argument('--mp', type=int, default=None, metavar='DIGITS', help='mpmath at this many digits')
    ap.add_argument('--tag', default='input')
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args(argv)
    target, pstr, rows = load(a.file)
    if a.p:
        pstr = a.p
    val, info = evaluate(target, Fraction(pstr), rows, a.mp, a.verbose)
    print(f'# {a.file}: {target} n={info["n"]} rank {info["rank"]} p={pstr}  mode {info["mode"]}  '
          f'trace-1 = {float(info["trace"] - 1):.1e}  {info["seconds"]:.1f} s')
    print(f'RESULT {a.tag} {target} n={info["n"]} p={pstr} {info["mode"]} I_c={fmt(val)}')
    return val


if __name__ == '__main__':
    main()
