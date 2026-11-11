"""Nelder-Mead over the mod-3 ansatz parameters (a0, sfac, q) at fixed n: usage scan_nm.py N A0 SFAC Q"""
import sys
sys.path.insert(0, '.')
import numpy as np, depopt
from ansatz import ansatz
from scipy.optimize import minimize
n = int(sys.argv[1])
MARGIN = float(sys.argv[5]) if len(sys.argv) > 5 else 2e-9
x0 = [float(t) for t in sys.argv[2:5]]
best = [0, None]
def f(x):
    a0, sf, q = x
    if not (2 < a0 < n - 2 and 0.3 < sf < 5 and 0.05 < q < 0.95):
        return 1.0
    th = depopt.threshold(ansatz(n, a0, sf * np.sqrt(a0), q).real, lo=0.05, margin=MARGIN) or 0.0
    if th > best[0]:
        best[:] = [th, (a0, sf, q)]
        print(f"n={n} a0={a0:.5f} sf={sf:.6f} q={q:.6f}: {th:.8f} (margin {MARGIN})", flush=True)
        depopt.save(ansatz(n, a0, sf * np.sqrt(a0), q), n, th, "nm_margin")
    return -th
minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-8, "initial_simplex": [x0, [x0[0] + 6, x0[1], x0[2]], [x0[0], x0[1] + 0.3, x0[2]], [x0[0], x0[1], x0[2] + 0.05]], "maxfev": 80})
print("best", best, flush=True)
