"""Fast ancilla-doubling ladder in the Schmidt basis (real):
    python3 work/lower1/pvdbl.py TARGET PKL M NTRIALS SEED [--sigma s0 s1] [--tmax s] [--save out.pkl] [--dbl K] [--iters N]
Compress PKL to active rank r, embed A (x) I_M, B (x) I_M, lam (x) 1/sqrt(M)  (n = r*M), or with --dbl K double only the K
smallest-weight levels (n = r+K); add sigma*GUE noise to the projector generators, iterate (pvit.iterate).  Saves best."""
import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
import sys, time, argparse, pickle
sys.path.insert(0, "qlib"); sys.path.insert(0, "work/lower1")
import numpy as np
import bell
import pvit

ap = argparse.ArgumentParser()
ap.add_argument("target"); ap.add_argument("pkl"); ap.add_argument("m", type=int)
ap.add_argument("ntrials", type=int); ap.add_argument("seed", type=int)
ap.add_argument("--sigma", type=float, nargs=2, default=[0.03, 0.5])
ap.add_argument("--tmax", type=float, default=300)
ap.add_argument("--iters", type=int, default=3000)
ap.add_argument("--dbl", type=int, default=0)
ap.add_argument("--save", default=None)
a = ap.parse_args()
c = bell.load(a.target)
rec = float(bell.record(a.target)["lower"])
rng = np.random.default_rng(a.seed)
st = pickle.load(open(a.pkl, "rb"))["st"]
A, B, lam = pvit.from_state(st)
r = int((np.abs(lam) > 1e-5).sum())
idx = np.argsort(-np.abs(lam))[:r]
A = A[:, idx][:, :, idx]; B = B[:, idx][:, :, idx]; lam = lam[idx]
lam = lam / np.linalg.norm(lam)
A = np.array([bell._pos_projector(P - 0.5 * np.eye(r)) for P in A]); B = np.array([bell._pos_projector(P - 0.5 * np.eye(r)) for P in B])
v0, A, B, lam = pvit.iterate(c, A, B, lam, 20000)
print(f"active rank {r}, base value {v0:.10f} (vs PV rec {v0-rec:+.3e}); lam", np.round(np.abs(lam), 4), flush=True)
order = np.argsort(-np.abs(lam)); A = A[:, order][:, :, order]; B = B[:, order][:, :, order]; lam = lam[order]

if a.dbl:
    K = a.dbl
    mult = [1] * (r - K) + [2] * K
    n = sum(mult)
    V = np.zeros((n, r)); off = 0
    for i, mi in enumerate(mult):
        V[off:off + mi, i] = 1.0 / np.sqrt(mi); off += mi
else:
    n = r * a.m
    mult = [a.m] * r
    V = np.kron(np.eye(r), np.ones((a.m, 1)) / np.sqrt(a.m))
emb = (lambda P: np.kron(P, np.eye(a.m))) if not a.dbl else (lambda P: V @ P @ V.T)
# Schmidt vector of the embedded state: copy k of level i gets lam_i / sqrt(mult_i)
lam_n = np.concatenate([np.full(mi, lam[i] / np.sqrt(mi)) for i, mi in enumerate(mult)])
lam_n = lam_n / np.linalg.norm(lam_n)

def gue(n):
    X = rng.standard_normal((n, n))
    return (X + X.T) / 2

best, bst = v0, (A, B, lam)
vals = []
t0 = time.time(); t = 0
while t < a.ntrials and time.time() - t0 < a.tmax:
    t += 1
    sg = float(np.exp(rng.uniform(np.log(a.sigma[0]), np.log(a.sigma[1]))))
    An = np.array([bell._pos_projector(emb(P) - 0.5 * np.eye(n) + sg * gue(n)) for P in A])
    Bn = np.array([bell._pos_projector(emb(P) - 0.5 * np.eye(n) + sg * gue(n)) for P in B])
    lam0 = lam_n * (1 + 0.05 * rng.standard_normal(n))
    try:
        v, A2, B2, l2 = pvit.iterate(c, An, Bn, lam0 / np.linalg.norm(lam0), a.iters)
    except np.linalg.LinAlgError:
        continue
    vals.append(v)
    if v > best + 1e-10:
        best, bst = v, (A2, B2, l2)
        print(f"  trial {t} sigma {sg:.3f}: new best {best:.10f} (vs PV rec {best-rec:+.3e}) ({time.time()-t0:.0f}s)", flush=True)
        if a.save:
            pickle.dump({"c": a.target, "st": pvit.to_state(*bst), "v": best}, open(a.save, "wb"))
vals = np.array(vals)
print(f"{a.target} r={r} m={a.m} dbl={a.dbl} -> n={n}: {t} trials {time.time()-t0:.0f}s; best {best:.10f} (base {v0:.10f}, vs PV rec {best-rec:+.3e}); top:",
      " ".join(f"{x:.7f}" for x in np.unique(np.round(vals, 7))[::-1][:8]), flush=True)
