"""Gaussian-windowed sinc interpolation of trig polynomials from an oversampled grid, and its
transpose (a type-1 nonuniform cosine/sine transform).

A real trig polynomial of degree < n sampled at theta_r = pi r / N, r = 0..N (N = R n) is
evaluated at arbitrary theta in [0, pi] from the 2W + 1 nearest samples; relative error
~ exp(-W pi (1 - 1/R) / 2) (3.9e-14 at R = 8, W = 24).  Because the rule reproduces cos(k .) and
sin(k .) for k < n exactly (to that error), its transpose gives sum_s y_s cos(k theta_s) and
sum_s y_s sin(k theta_s) for all k < n from one sparse product and one DCT/DST (`cos_sums`,
`sin_sums`).
"""
import numpy as np
import scipy.sparse as sp
from scipy import fft as sfft


class Interp:
    def __init__(self, n, R, W=24):
        self.n, self.N, self.W = n, n * R, W
        self.a = np.pi * (1 - 1 / R) / (2 * W)

    def matrix(self, th, odd=False):
        N, W = self.N, self.W
        x = np.asarray(th) * N / np.pi
        r = np.rint(x).astype(np.int64)[:, None] + np.arange(-W, W + 1)
        y = x[:, None] - r
        w = np.sinc(y) * np.exp(-self.a * y * y)
        lo, hi = r < 0, r > N
        if odd:
            w = np.where(lo | hi, -w, w)
        r = np.where(lo, -r, np.where(hi, 2 * N - r, r))
        rows = np.repeat(np.arange(len(x)), 2 * W + 1)
        return sp.csr_matrix((w.ravel(), (rows, r.ravel())), shape=(len(x), N + 1))

    def cos_sums(self, th, y, A=None):
        """sum_s y_s cos(k th_s), k = 0..n-1 (A: a cached matrix(th))"""
        g = (self.matrix(th) if A is None else A).T @ y
        g[1:-1] *= 0.5
        return sfft.dct(g, type=1)[:self.n]

    def sin_sums(self, th, y, A=None):
        """sum_s y_s sin(k th_s), k = 0..n-1 (A: a cached matrix(th, odd=True))"""
        g = (self.matrix(th, odd=True) if A is None else A).T @ y
        out = np.zeros(self.n)
        out[1:] = 0.5 * sfft.dst(g[1:-1], type=1)[:self.n - 1]
        return out


def grid_cos(c, N):
    """c_0 + 2 sum_{j>=1} c_j cos(j theta) at theta_r = pi r / N, r = 0..N"""
    x = np.zeros(N + 1)
    x[:len(c)] = c
    return sfft.dct(x, type=1)


def grid_sin(c, N):
    """2 sum_{j>=1} c_j sin(j theta) at theta_r = pi r / N"""
    y = np.zeros(N - 1)
    y[:len(c) - 1] = c[1:]
    out = np.zeros(N + 1)
    out[1:N] = sfft.dst(y, type=1)
    return out
