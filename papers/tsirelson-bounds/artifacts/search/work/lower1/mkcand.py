"""python3 mkcand.py PKL OUT.json FIELD [--den 12] [--script-note ...]
Polish a pickled strategy with extra seesaw iterations and write an exact-rational candidate."""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import sys, pickle, argparse
sys.path.insert(0, "qlib")
import numpy as np
import bell, certs

ap = argparse.ArgumentParser()
ap.add_argument("pkl"); ap.add_argument("out"); ap.add_argument("field")
ap.add_argument("--den", type=int, default=12)
ap.add_argument("--construction", default="seesaw")
ap.add_argument("--script", default="work/lower1/ss.py")
a = ap.parse_args()
d = pickle.load(open(a.pkl, "rb"))
c = bell.load(d["c"])
st = d["st"]
A, B = st["A"], st["B"]
dA, dB = A[0].shape[0], B[0].shape[0]
I_A, I_B = np.eye(dA), np.eye(dB)
old = -1
for it in range(20000):
    W = bell.bell_operator(c, A, B)
    w, V = np.linalg.eigh(W)
    R = V[:, -1].reshape(dA, dB)
    for x in range(c["mA"]):
        Ox = c["bA"][x] * I_B + sum(c["bAB"][x, y] * B[y] for y in range(c["mB"]))
        A[x] = bell._pos_projector(R @ Ox.T @ R.conj().T)
    for y in range(c["mB"]):
        Oy = c["bB"][y] * I_A + sum(c["bAB"][x, y] * A[x] for x in range(c["mA"]))
        B[y] = bell._pos_projector(R.T @ Oy.T @ R.conj())
    if w[-1] - old < 1e-16:
        break
    old = w[-1]
st = {"A": A, "B": B, "psi": np.linalg.eigh(bell.bell_operator(c, A, B))[1][:, -1]}
v = bell.value(c, st)
print("polished value %.12f after %d iters; rank A" % (v, it), [int(round(np.trace(P).real)) for P in A], "B", [int(round(np.trace(P).real)) for P in B])
cand = certs.strategy_candidate(c, st, field=a.field, den=10 ** a.den,
                                author="lower1", construction=a.construction, script=a.script)
certs.save(cand, a.out)
print("wrote", a.out)
