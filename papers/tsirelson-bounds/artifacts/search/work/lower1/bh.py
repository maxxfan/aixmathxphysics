"""Basin hopping on seesaw optima.  python3 bh.py TARGET dA dB FIELD SEED --tmax s [--sigma s] [--T t]
Start: best of a few random seesaw runs.  Kick: perturb every projector's operator by sigma*GUE then
re-seesaw; Metropolis-like acceptance at temperature T (violation units)."""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import sys, time, argparse, pickle
sys.path.insert(0, "qlib")
import numpy as np
import bell

ap = argparse.ArgumentParser()
ap.add_argument("target"); ap.add_argument("dA", type=int); ap.add_argument("dB", type=int)
ap.add_argument("field"); ap.add_argument("seed", type=int)
ap.add_argument("--tmax", type=float, default=300)
ap.add_argument("--sigma", type=float, default=0.4)
ap.add_argument("--T", type=float, default=0.003)
ap.add_argument("--iters", type=int, default=1500)
ap.add_argument("--save", default=None)
ap.add_argument("--init", default=None)
a = ap.parse_args()
c = bell.load(a.target)
rec = float(bell.record(a.target)["lower"])
cplx = a.field == "complex"
rng = np.random.default_rng(a.seed)

def rh(d):
    X = rng.standard_normal((d, d)) + (1j * rng.standard_normal((d, d)) if cplx else 0)
    return (X + X.conj().T) / 2

def run_seesaw(A, B, iters=a.iters, tol=1e-13):
    old = -np.inf
    for it in range(iters):
        W = bell.bell_operator(c, A, B)
        w, V = np.linalg.eigh(W)
        psi = V[:, -1]
        R = psi.reshape(A[0].shape[0], B[0].shape[0])
        for x in range(c["mA"]):
            Ox = c["bA"][x] * np.eye(R.shape[1]) + sum(c["bAB"][x, y] * B[y] for y in range(c["mB"]))
            A[x] = bell._pos_projector(R @ Ox.T @ R.conj().T)
        for y in range(c["mB"]):
            Oy = c["bB"][y] * np.eye(R.shape[0]) + sum(c["bAB"][x, y] * A[x] for x in range(c["mA"]))
            B[y] = bell._pos_projector(R.T @ Oy.T @ R.conj())
        if w[-1] - old < tol:
            break
        old = w[-1]
    st = {"A": A, "B": B, "psi": np.linalg.eigh(bell.bell_operator(c, A, B))[1][:, -1]}
    return bell.value(c, st), st

if a.init:
    d0 = pickle.load(open(a.init, "rb"))
    st = d0["st"]
    st = {"A": [np.array(P, dtype=complex if cplx else float) for P in st["A"]],
          "B": [np.array(P, dtype=complex if cplx else float) for P in st["B"]], "psi": st["psi"]}
    v = bell.value(c, st)
else:
    v, st = bell.seesaw(c, a.dA, a.field, restarts=20, iters=400, seed=a.seed, dB=a.dB)
cur, cur_st = v, st
best, best_st = v, st
print("start", v, flush=True)
t0 = time.time(); n = 0; acc = 0
seen = {}
while time.time() - t0 < a.tmax:
    n += 1
    # kick: project (P + sigma*H) for a random subset of projectors
    A = [p.copy() for p in cur_st["A"]]; B = [p.copy() for p in cur_st["B"]]
    for x in range(len(A)):
        if rng.random() < 0.6:
            A[x] = bell._pos_projector(A[x] - 0.5 * np.eye(a.dA) + a.sigma * rh(a.dA))
    for y in range(len(B)):
        if rng.random() < 0.6:
            B[y] = bell._pos_projector(B[y] - 0.5 * np.eye(a.dB) + a.sigma * rh(a.dB))
    v, st = run_seesaw(A, B)
    seen[round(v, 6)] = seen.get(round(v, 6), 0) + 1
    if v > cur or rng.random() < np.exp((v - cur) / a.T):
        cur, cur_st = v, st; acc += 1
    if v > best:
        best, best_st = v, st
        print(f"  new best {best:.10f} at hop {n}  ({time.time()-t0:.0f}s) delta vs PV {best-rec:.3e}", flush=True)
top = sorted(seen.items(), reverse=True)[:8]
print(f"{a.target} bh ({a.dA},{a.dB}) {a.field} sigma={a.sigma} T={a.T}: {n} hops acc {acc}; best {best:.10f}; delta vs PV {best-rec:.3e}")
print("top distinct:", top)
if a.save:
    pickle.dump({"c": a.target, "st": best_st, "v": best}, open(a.save, "wb"))
