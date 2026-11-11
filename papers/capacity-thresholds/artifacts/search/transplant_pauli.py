"""Transplant an X-Z / 2-Pauli state to more copies: zero-pad its Dicke amplitudes in the Y-pole frame
(where the state is real and lives on low levels), map back, and find its float zero crossing and
margin threshold with the O(n^4) evaluator of results/indep (the pauliopt model is O(n^6)).
usage transplant_pauli.py FAMILY STATE.json N [N ...]   (saves {fam}_n{N}_transplant.json)"""
import pathlib
import sys
sys.path.insert(0, '.')
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "results" / "indep"))
import numpy as np, pauliopt, yframe
from fractions import Fraction
from scipy.optimize import brentq
import indep_ic

fam, f = sys.argv[1], sys.argv[2]
W0 = yframe.to_frame(pauliopt.load_rows(f))
n0 = W0.shape[1] - 1
MARGIN = 1e-8


def ic(V, p):
    return indep_ic.evaluate(fam, Fraction(p).limit_denominator(10 ** 12), V, probes=False)["Ic"]


for N in map(int, sys.argv[3:]):
    W = np.zeros((W0.shape[0], N + 1))
    W[:, :n0 + 1] = W0
    V = yframe.from_frame(W / np.linalg.norm(W))
    lo, hi = 0.120, 0.1255
    while ic(V, hi) > 0:
        hi += 0.001
    p0 = brentq(lambda p: ic(V, p), lo, hi, xtol=1e-9)
    pm = brentq(lambda p: ic(V, p) - MARGIN, lo, p0, xtol=1e-9)
    slope = (ic(V, p0) - ic(V, p0 - 1e-5)) / 1e-5
    print(f"{fam} n={n0} -> {N}: zero crossing {p0:.8f}, margin-{MARGIN:g} threshold {pm:.8f}, "
          f"dI_c/dp {slope:.3e}", flush=True)
    pauliopt.save(fam, V, pm, "transplant")
