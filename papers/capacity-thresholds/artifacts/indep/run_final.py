#!/usr/bin/env python3
"""indep_ic.py (double precision, master-equation representation; not rigorous) on the final
X-Z, two-Pauli and n = 240 depolarizing board records.  Added on 2026-11-10 for the replay
script of the release; indep_ic.py itself is unchanged.

Prints  RESULT <tag> <family> n=<n> p=<p> float64 I_c=<value>  lines.
"""
import os
import sys
from fractions import Fraction

for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import indep_ic   # noqa: E402

RUNS = [
    ('rec28_xz_n80', '20261027-124857_lead_2a327ea83c2b10f7.json'),
    ('rec27_2p_n80', '20261027-124857_lead_585247a2e5ee8280.json'),
    ('rec20_dep_n240', '20261027-053852_lead_9e078e22e51cd14b.json'),
]


def find(name):
    for parts in (('..', 'certificates', 'board'), ('..', '..', 'candidates')):
        path = os.path.join(HERE, *parts, name)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(name)


def main():
    for tag, name in RUNS:
        target, pstr, V, note = indep_ic.load_candidate(find(name))
        res = indep_ic.evaluate(target, Fraction(pstr), V, probes=False)
        print(f'RESULT {tag} {target} n={V.shape[1] - 1} p={pstr} float64 I_c={res["Ic"]:.12e}', flush=True)


if __name__ == '__main__':
    main()
