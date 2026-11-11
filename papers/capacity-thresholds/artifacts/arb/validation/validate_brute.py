#!/usr/bin/env python3
"""Validation (a): arb_ic.py against full 2^n brute force.

  float  : results/indep/brute.py (float64, error ~1e-15) for n = 2..8, all three families,
           ranks 1-3, three p each, random complex states (plus real states for depolarizing).
           Check: distance from the float value to the enclosure <= 1e-12, width tiny.
  hp     : hpbrute.py (mpmath, 50 digits, own eigensolver) for n <= 5 (rank 2) / n <= 4 (rank 3).
           Check: the 50-digit value lies inside the enclosure.

usage: validate_brute.py [--float-only | --hp-only] [--nmax-float 8]
"""
import argparse
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'indep'))
import arb_ic                 # noqa: E402
import rigcore as RC          # noqa: E402
import hpbrute                # noqa: E402
import brute                  # noqa: E402  (results/indep/brute.py, float64)
from flint import arb         # noqa: E402

TD = os.path.join(HERE, 'testdata')
os.makedirs(TD, exist_ok=True)
PS = {'depolarizing': ['0.02', '0.0625', '27/400'],
      'independent_xz': ['0.05', '0.1183', '0.13'],
      'two_pauli': ['0.04', '0.118', '0.15']}


def make(target, n, r, seed, real=False):
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(r, n + 1))
    if not real:
        V = V + 1j * rng.normal(size=(r, n + 1))
    rows = [[[float(np.real(z)), float(np.imag(z))] for z in row] for row in V]
    fn = os.path.join(TD, f'{target}_n{n}_r{r}_s{seed}{"_real" if real else ""}.json')
    with open(fn, 'w') as f:
        json.dump(dict(target=target, p=PS[target][0], rows=rows, note='validation state'), f)
    return fn, V, [[tuple(e) for e in row] for row in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--float-only', action='store_true')
    ap.add_argument('--hp-only', action='store_true')
    ap.add_argument('--nmax-float', type=int, default=8)
    ap.add_argument('--hp-extra', action='store_true',
                    help='only the larger 50-digit cases: n = 6 (ranks 2, 3) and n = 7 (rank 2)')
    args = ap.parse_args()
    worst_float = 0.0
    worst_width = 0.0
    nfail = 0
    ntests = 0
    t0 = time.time()
    if not args.hp_only and not args.hp_extra:
        print('== float brute force (results/indep/brute.py) ==', flush=True)
        for n in range(2, args.nmax_float + 1):
            for target in ('depolarizing', 'independent_xz', 'two_pauli'):
                for r in (1, 2, 3):
                    variants = [False, True] if target == 'depolarizing' else [False]
                    for real in variants:
                        seed = 1000 * n + 10 * r + (5 if real else 0)
                        fn, V, rows = make(target, n, r, seed, real)
                        for pstr in PS[target]:
                            R = arb_ic.evaluate(fn, p_override=pstr, nproc=1)
                            probs = [float(x) for x in RC.pauli_probs(target, Fraction(pstr))]
                            b = brute.ic_brute(V, probs)[0]
                            dist = max(float(R['lo'] - b), float(b - R['hi']), 0.0)
                            worst_float = max(worst_float, dist)
                            worst_width = max(worst_width, R['width'])
                            ok = dist <= 1e-12 and R['width'] < 1e-20
                            ntests += 1
                            nfail += (not ok)
                            print(f'{"ok  " if ok else "FAIL"} {target:15s} n {n} r {r} {"real" if real else "cplx"} '
                                  f'p {pstr:7s} I_c [{R["lo"].str(17, radius=False)}, {R["hi"].str(17, radius=False)}] '
                                  f'w {R["width"]:.1e}  float {b:+.16e}  dist {dist:.1e}', flush=True)
        print(f'float: {ntests} tests, {nfail} failures, max distance to enclosure {worst_float:.2e}, '
              f'max width {worst_width:.2e}, {time.time() - t0:.0f} s', flush=True)
    if not args.float_only:
        print('== 50-digit brute force (hpbrute.py, mpmath) ==', flush=True)
        nhp = 0
        nhpfail = 0
        worst_hp_width = 0.0
        cases = []
        for target in ('depolarizing', 'independent_xz', 'two_pauli'):
            for n in (2, 3, 4, 5):
                cases.append((target, n, 2, False))
            for n in (2, 3, 4):
                cases.append((target, n, 3, False))
            cases.append((target, 3, 1, False))
        cases.append(('depolarizing', 5, 2, True))
        cases.append(('depolarizing', 6, 2, True))
        if args.hp_extra:
            cases = [(t, n, r, False) for t in ('depolarizing', 'independent_xz', 'two_pauli')
                     for (n, r) in ((6, 2), (6, 3), (7, 2))]
        for target, n, r, real in cases:
            seed = 77000 + 100 * n + r + (50 if real else 0)
            fn, V, rows = make(target, n, r, seed, real)
            pstr = PS[target][1]
            t1 = time.time()
            R = arb_ic.evaluate(fn, p_override=pstr, nproc=1)
            t2 = time.time()
            v = hpbrute.ic_hp(rows, RC.pauli_probs(target, Fraction(pstr)), dps=50)
            t3 = time.time()
            vs = hpbrute.mp.nstr(v, 45)
            inside = R['lo'] <= arb(vs) and arb(vs) <= R['hi']
            # arb(vs) is a ball of radius ~1e-45 around the 45-digit decimal
            inside = bool(inside) or bool(R['lo'].overlaps(arb(vs)) or R['hi'].overlaps(arb(vs)))
            nhp += 1
            nhpfail += (not inside)
            worst_hp_width = max(worst_hp_width, R['width'])
            print(f'{"ok  " if inside else "FAIL"} {target:15s} n {n} r {r} {"real" if real else "cplx"} p {pstr}: '
                  f'enclosure [{R["lo"].str(30, radius=False)}, {R["hi"].str(30, radius=False)}] '
                  f'(w {R["width"]:.1e}); 50-digit {hpbrute.mp.nstr(v, 32)}  '
                  f'({t2 - t1:.1f} s / {t3 - t2:.1f} s)', flush=True)
        print(f'hp: {nhp} tests, {nhpfail} outside the enclosure, max width {worst_hp_width:.2e}', flush=True)


if __name__ == '__main__':
    main()
