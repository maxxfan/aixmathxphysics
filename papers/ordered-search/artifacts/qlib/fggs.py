"""FGGS polynomial chain for exact translation-invariant ordered search (float helpers).

Theorem (FGGS99; Thm 8 of arXiv:2503.21090): an exact k-query translation-invariant algorithm
for n-element ordered search exists iff there are symmetric Laurent polynomials q_0..q_k of degree
< n, nonnegative on the unit circle, with q_0 = F_n (Fejer kernel), q_k = 1, constant term 1, and
q_t = q_{t-1} on the n roots of z^n = (-1)^t.

Coefficient form: q(theta) = c[0] + 2 sum_{j=1}^{n-1} c[j] cos(j theta).

Value form: sample on the 2n-grid theta_m = pi m / n, half grid m = 0..n (symmetric).  Roots of
z^n = 1 are the even m, roots of z^n = -1 the odd m.  So q_t and q_{t-1} share the odd samples when
t is odd and the even samples when t is even, and the chain is a path of half-vectors

    u_0 = odd samples of F_n,  u_1 = even samples of q_1,  u_2 = odd samples of q_2, ...,
    u_{k-1} = 1,

where q_t interleaves u_{t-1} and u_t, and u_t sits on the even points for odd t and on the odd
points for even t (see `chain_values`).
Each u has full-circle mean 1 (that is c[0] = 1 and c[n] = 0).
"""
import numpy as np
from scipy import fft as sfft


def half_parity(n):
    m = np.arange(n + 1)
    return (m % 2 == 0), (m % 2 == 1)


def circle_mult(n):
    """multiplicity of half-grid point m on the full 2n-circle"""
    w = np.full(n + 1, 2.0)
    w[0] = 1.0
    w[n] = 1.0
    return w


def fejer_values(n):
    """F_n on the half grid m = 0..n"""
    m = np.arange(n + 1)
    v = np.zeros(n + 1)
    v[0] = n
    odd = m % 2 == 1
    v[odd] = 1.0 / (n * np.sin(np.pi * m[odd] / (2 * n)) ** 2)
    return v


def coeffs_from_values(v):
    """half-grid values (length n+1) -> c[0..n] (c[n] must vanish for degree < n)"""
    n = len(v) - 1
    return sfft.dct(v, type=1) / (2 * n)


def values_from_coeffs(c, R=1):
    """c[0..n] -> q at theta = pi r/(n R), r = 0..nR"""
    n = len(c) - 1
    N = n * R
    x = np.zeros(N + 1)
    x[: n + 1] = c
    # dct1: y_r = x_0 + (-1)^r x_N + 2 sum_{j=1}^{N-1} x_j cos(pi j r / N)
    return sfft.dct(x, type=1)


def chain_values(u, n):
    """u = [u_0, ..., u_{k-1}] half-vectors on their parity points (u_0 odd points).
    Returns list of full half-grid value vectors for q_1..q_{k-1}.
    q_t interleaves u_{t-1} and u_t; u_t lives on even points if t odd, odd points if t even."""
    ev, od = half_parity(n)
    out = []
    for t in range(1, len(u)):
        v = np.empty(n + 1)
        a, b = u[t - 1], u[t]
        if t % 2 == 1:      # u_{t-1} on odd, u_t on even
            v[od] = a
            v[ev] = b
        else:
            v[ev] = a
            v[od] = b
        out.append(v)
    return out
