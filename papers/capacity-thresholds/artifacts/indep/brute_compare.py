"""Compare indep_ic (3 exponential methods x 2 frames) with the 2^n brute force, n<=8, ranks 1-3, all families."""
import sys, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
import indep_ic as I
import brute as B

rng = np.random.default_rng(7)
worst = 0.0
for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    for r in [1, 2, 3]:
        for target, p in [('depolarizing', Fraction(3, 50)), ('independent_xz', Fraction(11, 100)),
                          ('two_pauli', Fraction(23, 200))]:
            V = rng.normal(size=(r, n + 1)) + 1j * rng.normal(size=(r, n + 1))
            probs = [float(x) for x in I.pauli_probs(target, p)]
            ib = B.ic_brute(V, probs)[0]
            line = [f'n={n} r={r} {target:15s} brute {ib:+.15f}']
            for method in ['uniform', 'expm', 'dense']:
                for frame in ['plain', 'rot']:
                    if method == 'dense' and frame == 'plain' and target != 'depolarizing':
                        continue
                    res = I.evaluate(target, p, V, method=method, frame=frame)
                    dd = res['Ic'] - ib
                    worst = max(worst, abs(dd))
                    line.append(f'{method[0]}{frame[0]} {dd:+.1e}')
            print('  '.join(line), flush=True)
print('worst', worst)
