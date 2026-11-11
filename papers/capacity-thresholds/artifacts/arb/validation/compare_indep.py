#!/usr/bin/env python3
"""Validation (c)/(d): run results/indep/indep_ic.py (float64 auditor tool) on a candidate at a
given p and print its value next to an arb_ic.py JSON result (if given).

usage: compare_indep.py FILE [--p P] [--arb-json RESULT.json] [--fast]
"""
import argparse
import json
import os
import sys
from fractions import Fraction

for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'indep'))
import indep_ic   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('--p', default=None)
    ap.add_argument('--arb-json', default=None)
    ap.add_argument('--fast', action='store_true')
    a = ap.parse_args()
    target, pstr, V, note = indep_ic.load_candidate(a.file)
    if a.p:
        pstr = a.p
    res = indep_ic.evaluate(target, Fraction(pstr), V, probes=not a.fast)
    ic = res['Ic']
    print(f'indep_ic  {a.file}  p {pstr}:  I_c = {ic:.15e} bits  ({res["times"]["total"]:.1f} s)')
    if not a.fast:
        d = res['diag']
        print(f'          |evd-evr| {abs(d["eig_alt"] - ic):.1e}  |probe| {abs(d["probe"] - ic):.1e}  '
              f'first-order bound {d["bound1"]:.1e}  floor {d["floor"]:.1e}')
    if a.arb_json:
        with open(a.arb_json) as f:
            R = json.load(f)
        from flint import arb
        lo, hi = arb(R['lo']), arb(R['hi'])
        mid = (lo + hi) / 2
        print(f'arb_ic    [{R["lo"]}, {R["hi"]}]  sign {R["sign"]}')
        print(f'          indep - arb_ic midpoint = {float(arb(repr(ic)) - mid):+.3e} bits; '
              f'inside: {bool(lo <= arb(repr(ic)) <= hi)}')


if __name__ == '__main__':
    main()
