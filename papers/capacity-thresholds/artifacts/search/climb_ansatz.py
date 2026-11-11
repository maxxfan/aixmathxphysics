"""climb from the mod-3 Gaussian ansatz: usage climb_ansatz.py N A0 SFAC Q ROUNDS ITERS"""
import sys
sys.path.insert(0, '.')
import numpy as np, depopt
from ansatz import ansatz
n, a0, sf, q, rounds, iters = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
V = ansatz(n, a0, sf * np.sqrt(a0), q)
th = depopt.threshold(V, lo=0.05)
print(f"n={n} ansatz a0={a0} s={sf}sqrt(a0) q={q}: {th}", flush=True)
for r in range(rounds):
    W = depopt.optimise(V, th + 3e-5, iters)
    t2 = depopt.threshold(W, lo=0.05)
    print(f"  round {r}: {th:.8f} -> {t2}", flush=True)
    if t2 is None or t2 <= th:
        break
    V, th = W, t2
    depopt.save(V, n, th, f"ans{r}")
print(f"final n={n}: {th:.8f}", flush=True)
