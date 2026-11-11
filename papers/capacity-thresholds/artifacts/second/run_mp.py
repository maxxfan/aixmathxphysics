#!/usr/bin/env python3
"""kernel_ic.py in mpmath (30 digits; not rigorous, but with 20 more digits than a double) on the
inputs that are small enough for a pure-Python multiprecision eigensolver:

  the two 45-copy states of Krohn-Grimberghe (arXiv:2608.15870) at their four certified p,
  our re-optimised 45-copy input at p = 6497/100000 (Krohn-Grimberghe's exact verifier gives
  I_c in [5.34543241111018806874e-9, 5.34543241111018807061e-9]),
  the final X-Z and two-Pauli board records (n = 80),
  the X-Z record padded to n = 120 at p = 0.1259 (above the frozen verifier's size limit for
  this family; the verifier's enclosure there contains zero).

The n = 45 runs take about 30 s each, the n = 80 runs about 5 minutes each, the n = 120 run
17 minutes; the cost grows like n^4.  The headline inputs (n = 320, 480) are out of reach of
this mode; run_headline.py evaluates them in double precision.

usage: run_mp.py [tag ...]
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kernel_ic as K            # noqa: E402
from inputs import find          # noqa: E402

RUNS = [
    ('kgA_pos', 'kg2026_A_n45.json', '16239/250000'),
    ('kgA_neg', 'kg2026_A_n45.json', '203/3125'),
    ('kgR_pos', 'kg2026_R_n45.json', '4057/62500'),
    ('kgR_neg', 'kg2026_R_n45.json', '16229/250000'),
    ('ours_n45', 'depolarizing_n45_0.06497.json', '6497/100000'),
    ('rec28_xz_n80', '20261027-124857_lead_2a327ea83c2b10f7.json', None),
    ('rec27_2p_n80', '20261027-124857_lead_585247a2e5ee8280.json', None),
    ('xz_n120', 'independent_xz_n120_0.1259_board_path1.json', None),
]


def main():
    only = set(sys.argv[1:])
    for tag, name, p_over in RUNS:
        if only and tag not in only:
            continue
        target, pstr, rows = K.load(find(name))
        if p_over:
            pstr = p_over
        val, info = K.evaluate(target, Fraction(pstr), rows, mp_digits=30)
        print(f'# {name}: trace-1 = {float(info["trace"] - 1):.1e}, {info["seconds"]:.0f} s', flush=True)
        print(f'RESULT {tag} {target} n={info["n"]} p={pstr} {info["mode"]} I_c={K.fmt(val, 25)}', flush=True)


if __name__ == '__main__':
    main()
