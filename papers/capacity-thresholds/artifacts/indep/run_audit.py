#!/usr/bin/env python3
"""Audit driver: evaluates the board candidates, negative controls at p + 1e-6, and a
monotonicity grid, with the independent evaluator indep_ic.py.  Prints a log."""
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import indep_ic as I  # noqa: E402

PROB = os.path.dirname(HERE)
CANDS = [
    ('certificates/board/20261026-210916_lead_481ce3870c2d4e21.json', (1.5821e-9, 2.1906e-9)),
    ('certificates/board/20261026-230622_lead_32025a4e90f6200f.json', (2.9994e-9, 4.4447e-9)),
    ('certificates/board/20261026-210926_lead_c9b944ed9a6f8a4a.json', (9.9207e-9, 1.0197e-8)),
    ('certificates/board/20261026-210955_lead_ad42e025717074f4.json', (9.8282e-9, 1.0237e-8)),
    ('certificates/arb/depolarizing_n45_0.06497.json', (5.3454324111e-9, 5.3454324111e-9)),
]


def run(path, p, methods):
    target, pstr, V, note = I.load_candidate(os.path.join(PROB, path))
    out = {}
    for m, fr in methods:
        t0 = time.time()
        res = I.evaluate(target, Fraction(p), V, method=m, frame=fr)
        out[(m, fr)] = (res['Ic'], time.time() - t0, res['diag'])
    return target, V.shape, out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    for path, encl in CANDS:
        target, pstr, V, note = I.load_candidate(os.path.join(PROB, path))
        n = V.shape[1] - 1
        p0 = Fraction(pstr)
        if which in ('all', 'main'):
            meths = [('uniform', 'auto'), ('expm', 'auto')]
            if target != 'depolarizing':
                meths += [('dense', 'rot'), ('uniform', 'plain'), ('expm', 'plain')]
            elif n <= 60:
                meths += [('dense', 'plain')]
            _, _, out = run(path, pstr, meths)
            vals = [v[0] for v in out.values()]
            spread = max(vals) - min(vals)
            Ic, tt, dg = out[('uniform', 'auto')]
            inside = encl[0] <= Ic <= encl[1]
            print(f'{target:15s} n={n:3d} p={pstr:12s} I_c={Ic:.10e}  spread({len(vals)} variants)={spread:.1e}  '
                  f'floor={dg["floor"]:.1e}  enclosure=[{encl[0]:.4e},{encl[1]:.4e}] inside={inside}  t={tt:.1f}s',
                  flush=True)
            for k, v in out.items():
                print(f'      {k[0]:8s} {k[1]:6s} {v[0]:.13e}  ({v[1]:.1f} s)', flush=True)
        if which == 'mono' and path.startswith('certificates/board/20261026-210916'):
            grid = ['0.0640', '0.0650', '0.0655', '0.0660', '0.0662', '0.06628', '0.06629494',
                    '0.06630', '0.0664', '0.0666', '0.0670', '0.0680']
            prev = None
            for ps in grid:
                res = I.evaluate(target, Fraction(ps), V, method='uniform')
                mono = '' if prev is None else ('decreasing' if res['Ic'] < prev else 'NOT DECREASING')
                print(f'   MONO {target} n={n} p={ps:11s} I_c={res["Ic"]:+.10e}  {mono}', flush=True)
                prev = res['Ic']
        if which == 'cross':
            # zero crossing of I_c(p) for this fixed state (secant on 1e-10 decimal grid)
            f = lambda ps: I.evaluate(target, Fraction(ps), V, method='uniform')['Ic']  # noqa: E731
            pa = p0
            pb = p0 + (Fraction(1, 10 ** 4) if target == 'depolarizing' and n > 100 else Fraction(1, 10 ** 6))
            fa, fb = f(pa), f(pb)
            for _ in range(8):
                pc = pb - fb * (pb - pa) / (fb - fa)
                pc = Fraction(round(pc * 10 ** 10), 10 ** 10)
                if pc in (pa, pb):
                    break
                pa, fa, pb, fb = pb, fb, pc, f(pc)
                if abs(fb) < 1e-14:
                    break
            slope = (fb - fa) / float(pb - pa)
            print(f'  CROSS {target:15s} n={n:3d} claimed p={float(p0):.10f}  zero crossing p*={float(pb):.10f} '
                  f'(I_c there {fb:+.1e})  p*-p={float(pb - p0):.2e}  slope dI_c/dp={slope:.3e} bits', flush=True)
            for dp in (Fraction(-1, 10 ** 6), Fraction(1, 10 ** 6)):
                ps = pb + dp
                v = f(ps)
                print(f'        p*{"+" if dp > 0 else "-"}1e-6 = {float(ps):.10f}: I_c = {v:+.4e}', flush=True)
        if which in ('all', 'neg'):
            pn = p0 + Fraction(1, 10 ** 6)
            _, _, out = run(path, str(pn), [('uniform', 'auto'), ('expm', 'auto')])
            a = out[('uniform', 'auto')][0]
            b = out[('expm', 'auto')][0]
            print(f'   NEG {target:15s} n={n:3d} p+1e-6={float(pn):.10f}  I_c={a:.6e} (expm {b:.6e})  '
                  f'sign={"NEG" if a < 0 and b < 0 else "POS?"}', flush=True)


if __name__ == '__main__':
    main()
