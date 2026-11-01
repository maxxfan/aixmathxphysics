"""python3 work/lower1/mkcand3.py PKL OUT.json [--den 12] [--construction ...] [--script ...]
Schmidt-basis state pkl (pvit/pvdbl) -> compress to active rank, polish with pvit.iterate, exact rational candidate (real)."""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import sys, pickle, argparse
sys.path.insert(0, "qlib"); sys.path.insert(0, "work/lower1")
import numpy as np
import bell, certs, pvit

ap = argparse.ArgumentParser()
ap.add_argument("pkl"); ap.add_argument("out")
ap.add_argument("--den", type=int, default=12)
ap.add_argument("--construction", default="Schmidt-basis iteration + ancilla doubling ladder")
ap.add_argument("--script", default="work/lower1/pvdbl.py")
a = ap.parse_args()
d = pickle.load(open(a.pkl, "rb"))
c = bell.load(d["c"])
A, B, lam = pvit.from_state(d["st"])
r = int((np.abs(lam) > 1e-7).sum())
idx = np.argsort(-np.abs(lam))[:r]
A = A[:, idx][:, :, idx]; B = B[:, idx][:, :, idx]; lam = lam[idx] / np.linalg.norm(lam[idx])
A = np.array([bell._pos_projector(P - 0.5 * np.eye(r)) for P in A]); B = np.array([bell._pos_projector(P - 0.5 * np.eye(r)) for P in B])
v, A, B, lam = pvit.iterate(c, A, B, lam, 200000, tol=1e-17)
st = pvit.to_state(A, B, lam)
print("active rank", r, "polished value %.12f" % bell.value(c, st), "rank A", [int(round(np.trace(P))) for P in A], "B", [int(round(np.trace(P))) for P in B], "lam", np.round(np.abs(lam), 4))
cand = certs.strategy_candidate(c, st, field="real", den=10 ** a.den, author="lower1", construction=a.construction, script=a.script)
certs.save(cand, a.out)
print("wrote", a.out, "dA", cand["dA"], "dB", cand["dB"])
