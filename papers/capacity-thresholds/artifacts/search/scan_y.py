"""Nelder-Mead over (a0, sfac, q) for the mod-3 ansatz in the Y-pole frame (real arithmetic):
usage scan_y.py FAMILY N A0 SFAC Q [MAXFEV]   (saves {fam}_n{N}_y_nm.json in the target frame)"""
import sys
sys.path.insert(0, '.')
import numpy as np, pauliopt, yframe
from ansatz import ansatz
from scipy.optimize import minimize
fam, n = sys.argv[1], int(sys.argv[2])
x0 = [float(t) for t in sys.argv[3:6]]
maxfev = int(sys.argv[6]) if len(sys.argv) > 6 else 70
best = [0, None]
def f(x):
    a0, sf, q = x
    if not (2 < a0 < n - 2 and 0.3 < sf < 5 and 0.05 < q < 0.95):
        return 1.0
    V = ansatz(n, a0, sf * np.sqrt(a0), q).real
    th = pauliopt.threshold(fam + "@y", V, margin=1e-8) or 0.0
    if th > best[0]:
        best[:] = [th, (a0, sf, q)]
        print(f"{fam} n={n} a0={a0:.3f} sf={sf:.4f} q={q:.4f}: {th:.8f}", flush=True)
        pauliopt.save(fam, yframe.from_frame(V), th, "y_nm")
    return -th
minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-8, "maxfev": maxfev,
         "initial_simplex": [x0, [x0[0] + 4, x0[1], x0[2]], [x0[0], x0[1] + 0.2, x0[2]], [x0[0], x0[1], x0[2] + 0.04]]})
print("best", best, flush=True)
