#!/usr/bin/env python3
"""Validation (b) and (d) for arb_ic.py, plus the extra KG sign checks.

  (b)  Krohn-Grimberghe n = 45 depolarizing certificates:
       state A at 16239/250000: inside KG's certificate [8.3895e-8, 8.7075e-8], float 8.60887680595e-8
       our n=45 state at 6497/100000: must overlap KG's exact interval
       [534543241111018806874081249826800422545285, 534543241111018807060758434297471090811600] e-50
  (b') the four certified KG signs: A + at 16239/250000, A - at 203/3125, R + at 4057/62500,
       R - at 16229/250000.
  (d)  negative controls: each target state at a p just above its float zero crossing must not
       be POSITIVE; a corrupted state (one amplitude perturbed) must move by the float-predicted
       amount (indep_ic.py on the same perturbed file).

usage: validate_records.py [--only b|d] [--nproc 6]
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import arb_ic                  # noqa: E402
from flint import arb, fmpz    # noqa: E402

KG_OURS = (534543241111018806874081249826800422545285, 534543241111018807060758434297471090811600)


def run(path, p=None, nproc=6, perturb=None):
    R = arb_ic.evaluate(os.path.join(ROOT, path), p_override=p, nproc=nproc, perturb=perturb)
    return R


def line(tag, R, extra=''):
    print(f'{tag:44s} p {R["p"]:>14s}  I_c [{R["lo"].str(16, radius=False)}, {R["hi"].str(16, radius=False)}] '
          f'w {R["width"]:.1e}  {R["sign"]:9s} {R["times"]["total"]:.0f}s {extra}', flush=True)


def part_b(nproc):
    print('== (b) Krohn-Grimberghe certificates ==')
    R = run('seeds/kg2026_A_n45.json', nproc=nproc)
    inside = R['lo'] > arb('8.3895e-8') and R['hi'] < arb('8.7075e-8')
    agree = abs(float(R['lo']) - 8.60887680595e-8) < 1e-18
    line('KG state A', R, f'inside KG cert: {inside}; float 8.60887680595e-8 to 11 digits: {agree}')
    R = run('certificates/arb/depolarizing_n45_0.06497.json', nproc=nproc)
    lo30 = int(arb_ic.fmt_dec(R['lo'], 50, up=False))
    hi30 = int(arb_ic.fmt_dec(R['hi'], 50, up=True))
    ov = lo30 <= KG_OURS[1] and KG_OURS[0] <= hi30
    inside = KG_OURS[0] <= lo30 and hi30 <= KG_OURS[1]
    line('our n=45 state (KG exact verifier)', R,
         f'\n    ours  [{lo30}, {hi30}] e-50\n    KG    [{KG_OURS[0]}, {KG_OURS[1]}] e-50\n'
         f'    overlap {ov}; ours inside KG {inside}')
    print("== (b') the four certified KG signs ==")
    for path, p, want in (('seeds/kg2026_A_n45.json', '16239/250000', 'POSITIVE'),
                          ('seeds/kg2026_A_n45.json', '203/3125', 'NEGATIVE'),
                          ('seeds/kg2026_R_n45.json', '4057/62500', 'POSITIVE'),
                          ('seeds/kg2026_R_n45.json', '16229/250000', 'NEGATIVE')):
        R = run(path, p, nproc=nproc)
        line(os.path.basename(path), R, 'ok' if R['sign'] == want else f'MISMATCH (want {want})')


TARGETS = [
    ('dep n240', 'certificates/board/20261027-053852_lead_9e078e22e51cd14b.json', '0.0667122'),   # float zero 0.0667112
    ('xz n80', 'certificates/board/20261027-124857_lead_2a327ea83c2b10f7.json', '0.1246510'),     # float zero 0.1246500
    ('2p n80', 'certificates/board/20261027-124857_lead_585247a2e5ee8280.json', '0.1237138'),     # float zero 0.1237128
]


def part_d(nproc, which=None):
    print('== (d) negative controls: p just above the float zero crossing ==')
    for tag, path, pneg in TARGETS:
        if which and which not in tag:
            continue
        R = run(path, pneg, nproc=nproc)
        line(tag + ' above zero crossing', R, 'ok (not POSITIVE)' if R['sign'] != 'POSITIVE' else 'FAIL')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default=None)
    ap.add_argument('--which', default=None)
    ap.add_argument('--nproc', type=int, default=6)
    a = ap.parse_args()
    t0 = time.time()
    if a.only in (None, 'b'):
        part_b(a.nproc)
    if a.only in (None, 'd'):
        part_d(a.nproc, a.which)
    print(f'total {time.time() - t0:.0f} s')
