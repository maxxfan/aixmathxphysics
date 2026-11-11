"""margin-aware gradient climb for X-Z / 2-Pauli in the Y-pole frame (real arithmetic):
usage climb_y.py FAMILY STATE.json MARGIN ROUNDS ITERS TAG
STATE is in the target frame (as saved by pauliopt.save); saves {fam}_n{n}_{TAG}.json, also in the target frame."""
import sys
sys.path.insert(0, '.')
import numpy as np, pauliopt, yframe
from scipy.optimize import minimize
fam, f, margin, rounds, iters, tag = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
fy = fam + "@y"
V = yframe.to_frame(pauliopt.load_rows(f))
V /= np.linalg.norm(V)
n = V.shape[1] - 1
th = pauliopt.threshold(fy, V, margin=margin)
print(f"{fam} n={n} start threshold(margin {margin}) {th:.8f}", flush=True)
for r in range(rounds):
    m = pauliopt.model(fy, n, th + 1e-5)
    def fun(x):
        I, g = m.value_grad(x.reshape(V.shape))
        return -I, -g.ravel()
    res = minimize(fun, V.ravel(), jac=True, method="L-BFGS-B", options={"maxiter": iters, "gtol": 1e-14, "ftol": 1e-16})
    W = res.x.reshape(V.shape); W /= np.linalg.norm(W)
    t2 = pauliopt.threshold(fy, W, margin=margin)
    print(f"  round {r}: {th:.8f} -> {t2}  ({res.nit} its)", flush=True)
    if t2 is None or t2 <= th:
        break
    V, th = W, t2
    pauliopt.save(fam, yframe.from_frame(V), th, tag)
print(f"final {fam} n={n}: {th:.8f}", flush=True)
