"""time one depolarizing I_c evaluation (and value+grad) of the mod-3 ansatz at several n:
usage time_n.py N [N ...]"""
import sys
import time
sys.path.insert(0, '.')
import numpy as np, depopt
from ansatz import ansatz

for n in map(int, sys.argv[1:]):
    V = ansatz(n, 36.3, 1.324 * np.sqrt(36.3), 0.5106).real
    t = time.time()
    I = depopt.ic_bits(V, 0.0665)
    t1 = time.time() - t
    print(f"n={n}: I_c(0.0665) = {I:.6e} bits, value {t1:.1f} s", flush=True)
