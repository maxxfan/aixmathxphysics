#!/usr/bin/env python3
"""Second evaluation (double precision, kernel formula; not rigorous) of the three headline
inputs: at the certified p, where I_c must be positive and close to the value certified by
arb_ic.py, and at p + 1.1e-5, where the same input has I_c < 0 (the float zero crossing lies
1.0e-5 above the certified p).  Also the final board records.  About 5 minutes, one process,
1.5 GB for n = 480.

Prints  RESULT <tag> <family> n=<n> p=<p> float64 I_c=<value>  lines (checked by
check_outputs.py second second_headline against expected/second.json).
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kernel_ic as K            # noqa: E402
from inputs import find          # noqa: E402

RUNS = [
    ('rec28_xz_n80', '20261027-124857_lead_2a327ea83c2b10f7.json', None),
    ('rec27_2p_n80', '20261027-124857_lead_585247a2e5ee8280.json', None),
    ('rec20_dep_n240', '20261027-053852_lead_9e078e22e51cd14b.json', None),
    ('xz_n240', 'independent_xz_n240_0.12737956.json', None),
    ('headline_xz_n320', 'independent_xz_n320_0.12770564.json', None),
    ('headline_xz_n320_above', 'independent_xz_n320_0.12770564.json', '0.12771664'),
    ('headline_2p_n320', 'two_pauli_n320_0.12677472.json', None),
    ('headline_2p_n320_above', 'two_pauli_n320_0.12677472.json', '0.12678572'),
    ('rec30_dep_n480', '20261028-024810_lead_5f09c5df10a5804f.json', None),
    ('headline_dep_n480', 'depolarizing_n480_0.06688201.json', None),
    ('headline_dep_n480_above', 'depolarizing_n480_0.06688201.json', '0.06689301'),
]


def main():
    only = set(sys.argv[1:])
    for tag, name, p_over in RUNS:
        if only and tag not in only:
            continue
        target, pstr, rows = K.load(find(name))
        if p_over:
            pstr = p_over
        val, info = K.evaluate(target, Fraction(pstr), rows)
        print(f'# {name}: trace-1 = {info["trace"] - 1:.1e}, H(out) = {info["H_out"]:.12f} bits, '
              f'{info["seconds"]:.1f} s', flush=True)
        print(f'RESULT {tag} {target} n={info["n"]} p={pstr} {info["mode"]} I_c={K.fmt(val)}', flush=True)


if __name__ == '__main__':
    main()
