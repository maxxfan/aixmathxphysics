"""margin-aware gradient climb in real arithmetic (depolarizing):
usage climb_margin.py STATE.json MARGIN ROUNDS ITERS TAG"""
import sys, json
sys.path.insert(0, '.')
import numpy as np, depopt
from scipy.optimize import minimize
f, margin, rounds, iters, tag = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
d = json.load(open(f))
V = np.array([[complex(x, y) for x, y in r] for r in d["rows"]])
assert np.abs(V.imag).max() == 0, "real states only"
V = V.real / np.linalg.norm(V.real)
n = V.shape[1] - 1
th = depopt.threshold(V, lo=0.05, margin=margin)
print(f"n={n} start threshold(margin {margin}) {th:.8f}", flush=True)
for r in range(rounds):
    m = depopt.Depol(n, th + 1e-5)
    def fun(x):
        I, g = m.value_grad(x.reshape(V.shape))
        return -I, -g.ravel()
    res = minimize(fun, V.ravel(), jac=True, method="L-BFGS-B", options={"maxiter": iters, "gtol": 1e-14, "ftol": 1e-16})
    W = res.x.reshape(V.shape); W /= np.linalg.norm(W)
    t2 = depopt.threshold(W, lo=0.05, margin=margin)
    print(f"  round {r}: {th:.8f} -> {t2}", flush=True)
    if t2 is None or t2 <= th:
        break
    V, th = W, t2
    depopt.save(V + 0j, n, th, tag)
print(f"final n={n}: {th:.8f}", flush=True)
