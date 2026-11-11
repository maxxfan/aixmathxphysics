#!/usr/bin/env python3
"""Unit test of rigeig.entropy_enclosure on matrices with a known spectrum.

A = H diag(lam) H^T with H a product of exact rational Householder reflections (exactly
orthogonal, so the spectrum of A is exactly lam).  The spectra mimic the hard cases: exactly
degenerate eigenvalues, near-degenerate clusters (gaps 1e-25 .. 1e-40), a geometric tail down to
1e-60, exact zeros, and a large dynamic range.  The exact sum h(lam) must lie in the enclosure.
"""
import os
import sys
from fractions import Fraction

import numpy as np
from flint import acb_mat, arb, arb_mat, ctx, fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rigeig as RE   # noqa: E402


def householder_product(d, k, rng, cplx):
    """exact rational orthogonal (real) matrix: product of k Householder reflections."""
    M = arb_mat(d, d)
    for i in range(d):
        M[i, i] = 1
    for _ in range(k):
        v = [Fraction(int(x), 7) for x in rng.integers(-9, 10, size=d)]
        if all(x == 0 for x in v):
            v[0] = Fraction(1)
        nv = sum(x * x for x in v)
        H = arb_mat(d, d)
        for i in range(d):
            for j in range(d):
                h = (1 if i == j else 0) - 2 * v[i] * v[j] / nv
                H[i, j] = arb(fmpq(h.numerator, h.denominator))
        M = H * M
    return M


def spectrum(d, rng, kind):
    lam = []
    if kind == 'gauss':
        # Gaussian envelope like the record states: lam_i ~ exp(-c i^2)
        for i in range(d):
            lam.append(Fraction(1, 2 ** (i * i // 6 + 1)))
    elif kind == 'clusters':
        base = [Fraction(1, 3), Fraction(1, 7), Fraction(1, 10 ** 5)]
        for b in base:
            lam += [b, b, b + Fraction(1, 10 ** 30), b + Fraction(1, 10 ** 40)]
        while len(lam) < d - 4:
            lam.append(Fraction(1, 10 ** (10 + 2 * len(lam))))
        lam += [Fraction(0)] * (d - len(lam))
    else:
        lam = [Fraction(int(x), 1000) for x in rng.integers(0, 1000, size=d)]
    s = sum(lam)
    return [x / s for x in lam]


def exact_entropy(lam):
    ln2 = arb(2).log()
    t = arb(0)
    for x in lam:
        if x > 0:
            a = arb(fmpq(x.numerator, x.denominator))
            t += -a * a.log() / ln2
    return t


def main():
    ctx.prec = 128
    rng = np.random.default_rng(7)
    ok = True
    for d, kind in ((20, 'gauss'), (40, 'gauss'), (20, 'clusters'), (30, 'clusters'), (25, 'random')):
        lam = spectrum(d, rng, kind)
        Q = householder_product(d, 4, rng, False)
        D = arb_mat(d, d)
        for i, x in enumerate(lam):
            D[i, i] = arb(fmpq(x.numerator, x.denominator))
        A = Q * D * Q.transpose()          # ball matrix containing the exact A
        for cplx in (False, True):
            H = acb_mat(A) if cplx else A
            lo, hi, info = RE.entropy_enclosure(H, cplx, arb(0))
            ex = exact_entropy(lam)
            inside = lo <= ex and ex <= hi
            ok &= bool(inside)
            print(f'{"ok  " if inside else "FAIL"} d {d:3d} {kind:9s} {"complex" if cplx else "real   "} '
                  f'exact {ex.str(25, radius=False)}  enclosure width {float(hi - lo):.1e}  '
                  f'components {info.get("ncomp")} (largest {info.get("maxcomp")}), clipped {info.get("nclip")}')
    print('ALL OK' if ok else 'FAILURES')


if __name__ == '__main__':
    main()
