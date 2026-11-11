"""Transplant a depolarizing state to more copies by zero-padding its Dicke amplitudes.

The record states live on low Dicke levels (a <~ 80), where Sym^n looks like a bosonic Fock space
(Holstein-Primakoff; the Gaussian limit of arXiv:2609.39747), so the same amplitudes should stay
near-optimal at larger n.  Compares with the mod-3 ansatz at the n=240 optimum parameters.
usage transplant.py STATE.json N [MARGIN]   (saves dep_n{N}_transplant.json)"""
import sys
sys.path.insert(0, '.')
import numpy as np, depopt, pauliopt
from ansatz import ansatz

V0 = pauliopt.load_rows(sys.argv[1])
N = int(sys.argv[2])
margin = float(sys.argv[3]) if len(sys.argv) > 3 else 4e-9
V = np.zeros((V0.shape[0], N + 1), dtype=complex)
V[:, :V0.shape[1]] = V0
V = V.real / np.linalg.norm(V)
th = depopt.threshold(V, lo=0.0655, margin=margin)
print(f"transplant n={V0.shape[1] - 1} -> {N}: threshold(margin {margin:g}) {th:.8f}", flush=True)
depopt.save(V + 0j, N, th, "transplant")
A = ansatz(N, 36.30612, 1.324115 * np.sqrt(36.30612), 0.510607).real
tha = depopt.threshold(A, lo=0.0655, margin=margin)
print(f"ansatz (n=240 parameters) at n={N}: threshold(margin {margin:g}) {tha:.8f}", flush=True)
