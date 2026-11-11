#!/usr/bin/env python3
"""Validation of kernel_ic.py against a full 2^n density-matrix computation written here.

  float : random complex inputs, n = 1..7, ranks 1-3, the three families at two noise levels each.
  mp    : the same comparison with kernel_ic in mpmath (30 digits) for n <= 5, which measures the
          error of the double-precision brute force rather than of kernel_ic.
  frame : the frame change W^{(x)n} = exp(-i pi/2 J_x) of kernel_ic (exact Krawtchouk integers)
          against a matrix exponential of J_x computed from its eigendecomposition, n <= 80.

Prints one line per case and a final line "kernel validation: N cases, F failures, ...".
"""
import math
import sys
import os
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kernel_ic as K       # noqa: E402

PAULI = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]


def brute(target, p, rows):
    """I_c in bits from the 2^n r dimensional output-reference density matrix."""
    n = len(rows[0]) - 1
    r = len(rows)
    probs = [float(x) for x in K.pauli_probs(target, p)]
    dim = 2 ** n
    psi = np.zeros((dim, r), dtype=complex)
    for x in range(dim):
        a = bin(x).count('1')
        for q in range(r):
            psi[x, q] = rows[q][a] / math.sqrt(math.comb(n, a))
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi.reshape(-1), psi.reshape(-1).conj()).reshape([2] * n + [r] + [2] * n + [r])
    for i in range(n):
        acc = 0
        for P, w in zip(PAULI, probs):
            if w == 0:
                continue
            t = np.tensordot(P, rho, axes=([1], [i]))
            t = np.moveaxis(t, 0, i)
            t = np.tensordot(t, P.conj().T, axes=([n + 1 + i], [0]))
            t = np.moveaxis(t, -1, n + 1 + i)
            acc = acc + w * t
        rho = acc
    joint = rho.reshape(dim * r, dim * r)
    out = np.einsum('xqyq->xy', rho.reshape(dim, r, dim, r))

    def H(M):
        w = np.linalg.eigvalsh((M + M.conj().T) / 2)
        w = w[w > 1e-300]
        return float(-(w * np.log2(w)).sum())
    return H(out) - H(joint)


def main():
    rng = np.random.default_rng(20261110)
    PS = {'depolarizing': ['3/50', '27/400'], 'independent_xz': ['11/100', '0.1277'],
          'two_pauli': ['23/200', '0.1268']}
    ncase = nfail = 0
    worst = 0.0
    print('== float: kernel_ic against brute force ==')
    for n in range(1, 8):
        for r in (1, 2, 3):
            for target in PS:
                for ps in PS[target]:
                    rows = [[complex(a, b) for a, b in zip(rng.normal(size=n + 1), rng.normal(size=n + 1))]
                            for _ in range(r)]
                    ref = brute(target, Fraction(ps), rows)
                    val, info = K.evaluate(target, Fraction(ps), rows)
                    dd = abs(val - ref)
                    ok = dd < 5e-13 and abs(info['trace'] - 1) < 1e-12
                    ncase += 1
                    nfail += not ok
                    worst = max(worst, dd)
                    print(f'{"ok  " if ok else "FAIL"} {target:15s} n {n} r {r} p {ps:7s} kernel {val:+.15f} '
                          f'brute {ref:+.15f} diff {dd:.1e}')
    print(f'float: worst |kernel - brute| = {worst:.2e}')
    print('== mp (30 digits): kernel_ic against brute force ==')
    worst_mp = 0.0
    for n in range(2, 6):
        for target in PS:
            ps = PS[target][1]
            rows = [[complex(a, b) for a, b in zip(rng.normal(size=n + 1), rng.normal(size=n + 1))]
                    for _ in range(2)]
            ref = brute(target, Fraction(ps), rows)
            val, info = K.evaluate(target, Fraction(ps), rows, mp_digits=30)
            fl, _ = K.evaluate(target, Fraction(ps), rows)
            dd = abs(float(val) - ref)
            df = abs(float(val) - fl)
            ok = dd < 5e-13 and df < 5e-13
            ncase += 1
            nfail += not ok
            worst_mp = max(worst_mp, df)
            print(f'{"ok  " if ok else "FAIL"} {target:15s} n {n} r 2 p {ps:7s} mp {K.fmt(val, 20)} '
                  f'|mp - brute| {dd:.1e} |mp - float| {df:.1e}')
    print(f'mp: worst |mp - float kernel| = {worst_mp:.2e}')
    print('== frame change: Krawtchouk integers against exp(-i pi/2 J_x) ==')
    worst_rot = 0.0
    for n in (5, 20, 48, 80):
        v = rng.normal(size=n + 1) + 1j * rng.normal(size=n + 1)
        got = np.array(K.rotate_rows([list(v)], K.Float())[0])
        a = np.arange(n)
        off = 0.5 * np.sqrt((a + 1) * (n - a))
        J = np.diag(off, 1) + np.diag(off, -1)
        w, Q = np.linalg.eigh(J)
        m = -n / 2 + np.arange(n + 1)                       # the exact spectrum of J_x
        ref = Q @ (np.exp(-1j * np.pi / 2 * m) * (Q.T @ v))
        dd = float(np.abs(got - ref).max())
        ok = dd < 1e-11 and abs(np.linalg.norm(got) - np.linalg.norm(v)) < 1e-12
        ncase += 1
        nfail += not ok
        worst_rot = max(worst_rot, dd)
        print(f'{"ok  " if ok else "FAIL"} n {n}: max |difference| {dd:.1e}')
    print(f'kernel validation: {ncase} cases, {nfail} failures, worst float {worst:.1e}, '
          f'worst mp-float {worst_mp:.1e}, worst frame {worst_rot:.1e}')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
