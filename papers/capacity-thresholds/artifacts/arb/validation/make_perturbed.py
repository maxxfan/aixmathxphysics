#!/usr/bin/env python3
"""Negative control (d): write copies of the target candidates with one amplitude perturbed.

The largest-modulus entry of row 0 gets Re += DELTA (the resulting double is what both tools
read).  Output: testdata/pert_<tag>.json.  Evaluate with arb_ic.py and results/indep/indep_ic.py.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
DELTA = 1e-4
FILES = {
    'dep_n240': 'certificates/board/20261027-053852_lead_9e078e22e51cd14b.json',
    'xz_n80': 'certificates/board/20261027-124857_lead_2a327ea83c2b10f7.json',
    '2p_n80': 'certificates/board/20261027-124857_lead_585247a2e5ee8280.json',
}


def num(x):
    if isinstance(x, str):
        s = x.strip()
        return float.fromhex(s) if s.lower().lstrip('+-').startswith('0x') else float(s)
    return float(x)


def main():
    os.makedirs(os.path.join(HERE, 'testdata'), exist_ok=True)
    for tag, path in FILES.items():
        with open(os.path.join(ROOT, path)) as f:
            d = json.load(f)
        row = d['rows'][0]
        mods = []
        for a, x in enumerate(row):
            re, im = (num(x[0]), num(x[1])) if isinstance(x, list) else (num(x), 0.0)
            mods.append((re * re + im * im, a, re, im))
        _, a, re, im = max(mods)
        new = re + DELTA
        row[a] = [new, im] if isinstance(row[a], list) else new
        d['note'] = f'NEGATIVE CONTROL: {path} with Re v_0[{a}] {re!r} -> {new!r}'
        out = os.path.join(HERE, 'testdata', f'pert_{tag}.json')
        with open(out, 'w') as f:
            json.dump(d, f)
        print(out, d['note'])


if __name__ == '__main__':
    main()
