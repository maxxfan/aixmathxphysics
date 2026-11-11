"""Nelder-Mead over (a0, sfac, q) for the Y-pole mod-3 ansatz, any family:
usage scan_nm2.py FAMILY N A0 SFAC Q"""
import sys
sys.path.insert(0, '.')
import numpy as np, pauliopt
from ansatz import ansatz, rotate
from scipy.optimize import minimize
fam, n = sys.argv[1], int(sys.argv[2])
x0 = [float(t) for t in sys.argv[3:6]]
best = [0, None]
def f(x):
    a0, sf, q = x
    if not (2 < a0 < n - 2 and 0.3 < sf < 5 and 0.05 < q < 0.95):
        return 1.0
    V = rotate(ansatz(n, a0, sf * np.sqrt(a0), q), "x", np.pi / 2)
    th = pauliopt.threshold(fam, V, margin=1e-8) or 0.0
    if th > best[0]:
        best[:] = [th, (a0, sf, q)]
        print(f"{fam} n={n} a0={a0:.3f} sf={sf:.4f} q={q:.4f}: {th:.8f}", flush=True)
        pauliopt.save(fam, V, th, "ypole_nm")
    return -th
minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-8, "maxfev": 70,
         "initial_simplex": [x0, [x0[0] + 4, x0[1], x0[2]], [x0[0], x0[1] + 0.3, x0[2]], [x0[0], x0[1], x0[2] + 0.05]]})
print("best", best, flush=True)
