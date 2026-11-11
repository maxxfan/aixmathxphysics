#!/usr/bin/env python3
"""Test of the series-truncation bound of arb_ic.py (added 2026-11-10 for the release; the
evaluator is unchanged).

arb_ic.py truncates the Poisson series (and, for two_pauli, the Taylor series of the negative-rate
factor) where a rigorous tail bound falls below `--tail` (default 1e-45), and widens every
eigenvalue interval by the resulting trace-norm bound.  With the default the truncation error is
far below the ball radii, so the validation against brute force does not exercise the bound.
Here the series are cut much earlier (tail targets 1e-30 ... 1e-6), so that the truncation error
dominates, and the enclosure must still contain the true value:

  (a) small inputs (n = 4, 5; all three families, ranks 2 and 3): the 50-digit brute-force value
      of hpbrute.py must lie in the enclosure for every tail target;
  (b) the final two-Pauli and X-Z board records (n = 80): the enclosure at a loose tail target
      must contain the 128-bit enclosure obtained with the default.

Prints one line per case and "tail test: N cases, F failures".
usage: test_tail.py [--small-only] [--nproc 4]
"""
import argparse
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))                 # release layout: arb/validation/
import arb_ic                 # noqa: E402
import rigcore as RC          # noqa: E402
import hpbrute                # noqa: E402
from flint import arb         # noqa: E402

TD = os.path.join(HERE, 'testdata')
PS = {'depolarizing': '0.0625', 'independent_xz': '0.1183', 'two_pauli': '0.118'}
TAILS = (1e-45, 1e-30, 1e-20, 1e-12, 1e-6)
BOARD = {      # final n = 80 board records and their default-tail enclosures (x 1e-30 bits)
    '20261027-124857_lead_585247a2e5ee8280.json': (4688004925413416462978, 4688004925413416462986),
    '20261027-124857_lead_2a327ea83c2b10f7.json': (1605068705609815603826, 1605068705609815603827),
}


def find_board(name):
    for parts in (('..', '..', 'certificates', 'board'), ('..', '..', 'candidates')):
        path = os.path.join(HERE, *parts, name)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--small-only', action='store_true')
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    ncase = nfail = 0
    print('== (a) 50-digit brute force inside the enclosure at every tail target ==', flush=True)
    for target in ('depolarizing', 'independent_xz', 'two_pauli'):
        for n, r in ((5, 2), (4, 3)):
            fn = os.path.join(TD, f'{target}_n{n}_r{r}_s{77000 + 100 * n + r}.json')
            with open(fn) as f:
                rows = [[tuple(e) for e in row] for row in json.load(f)['rows']]
            pstr = PS[target]
            v = hpbrute.ic_hp(rows, RC.pauli_probs(target, Fraction(pstr)), dps=50)
            vb = arb(hpbrute.mp.nstr(v, 45))
            for tail in TAILS:
                R = arb_ic.evaluate(fn, p_override=pstr, nproc=1, tail=tail)
                inside = bool(R['lo'] <= vb and vb <= R['hi']) or bool(R['lo'].overlaps(vb) or R['hi'].overlaps(vb))
                ncase += 1
                nfail += not inside
                print(f'{"ok  " if inside else "FAIL"} {target:15s} n {n} r {r} p {pstr:7s} tail {tail:.0e}: '
                      f'{R["poisson_terms"]:3d} Poisson, {R["taylor_terms"]:3d} Taylor terms, bound {R["tail_bound"]:.1e}, '
                      f'width {R["width"]:.1e}, 50-digit value inside: {inside}', flush=True)
    if not a.small_only:
        print('== (b) n = 80 board records: loose-tail enclosure contains the default one ==', flush=True)
        for name, (lo30, hi30) in BOARD.items():
            fn = find_board(name)
            for tail in (1e-20, 1e-10):
                R = arb_ic.evaluate(fn, nproc=a.nproc, tail=tail)
                scale = arb(10) ** 30
                contains = bool(R['lo'] * scale <= lo30) and bool(R['hi'] * scale >= hi30)
                ncase += 1
                nfail += not contains
                print(f'{"ok  " if contains else "FAIL"} {R["target"]:15s} n {R["n"]} p {R["p"]} tail {tail:.0e}: '
                      f'{R["poisson_terms"]:3d} Poisson, {R["taylor_terms"]:3d} Taylor terms, bound {R["tail_bound"]:.1e}, '
                      f'I_c in [{R["lo"].str(8, radius=False)}, {R["hi"].str(8, radius=False)}], '
                      f'contains the default enclosure: {contains}', flush=True)
    print(f'tail test: {ncase} cases, {nfail} failures')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
