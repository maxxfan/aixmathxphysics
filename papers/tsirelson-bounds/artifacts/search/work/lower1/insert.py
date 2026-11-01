"""Bulk-cell insertion: python3 insert.py TARGET PKL M SIGMAS KS SEED [--save out.pkl]
Go to Schmidt basis of a pair-chain strategy (r active dims = r/2 pairs), duplicate bulk cell k
(pair k + coupling to pair k+1) M times, keep the edge pair at the end, add sigma*GUE noise, seesaw.
KS, SIGMAS comma-separated lists."""
import os
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import sys, time, argparse, pickle
sys.path.insert(0, "qlib")
import numpy as np
import bell

ap = argparse.ArgumentParser()
ap.add_argument("target"); ap.add_argument("pkl"); ap.add_argument("m", type=int)
ap.add_argument("sigmas"); ap.add_argument("ks"); ap.add_argument("seed", type=int)
ap.add_argument("--ntry", type=int, default=3)
ap.add_argument("--iters", type=int, default=8000)
ap.add_argument("--save", default=None)
ap.add_argument("--log", default=None)
a = ap.parse_args()
c = bell.load(a.target)
rng = np.random.default_rng(a.seed)
d0 = pickle.load(open(a.pkl, "rb")); st = d0["st"]
A0 = [np.array(P, float) for P in st["A"]]; B0 = [np.array(P, float) for P in st["B"]]
dA = A0[0].shape[0]
R = np.asarray(st["psi"]).reshape(dA, dA); R /= np.linalg.norm(R)
U, s, Vh = np.linalg.svd(R)
r = int((s > 1e-5).sum()); npairs = r // 2
Ar = [(U.T @ P @ U)[:r, :r] for P in A0]
Br = [(Vh @ P @ Vh.T)[:r, :r] for P in B0]
print("pairs", npairs, "value", bell.value(c, st), flush=True)

def run_seesaw(A, B, iters=a.iters, tol=1e-13):
    old = -np.inf
    n1 = A[0].shape[0]
    for it in range(iters):
        W = bell.bell_operator(c, A, B)
        w, V = np.linalg.eigh(W)
        Rr = V[:, -1].reshape(n1, n1)
        for x in range(c["mA"]):
            Ox = c["bA"][x] * np.eye(n1) + sum(c["bAB"][x, y] * B[y] for y in range(c["mB"]))
            A[x] = bell._pos_projector(Rr @ Ox.T @ Rr.conj().T)
        for y in range(c["mB"]):
            Oy = c["bB"][y] * np.eye(n1) + sum(c["bAB"][x, y] * A[x] for x in range(c["mA"]))
            B[y] = bell._pos_projector(Rr.T @ Oy.T @ Rr.conj())
        if w[-1] - old < tol:
            break
        old = w[-1]
    stt = {"A": A, "B": B, "psi": np.linalg.eigh(bell.bell_operator(c, A, B))[1][:, -1]}
    return bell.value(c, stt), stt

def blk(P, i, j):
    return P[2 * i:2 * i + 2, 2 * j:2 * j + 2]

def build(P, k, m):
    """new pair list: old pairs 0..k, then m copies of pair k, then old k+1..; couplings follow cells."""
    order = list(range(k + 1)) + [k] * m + list(range(k + 1, npairs))
    N = len(order)
    M = np.zeros((2 * N, 2 * N))
    for j in range(N):
        M[2 * j:2 * j + 2, 2 * j:2 * j + 2] = blk(P, order[j], order[j])
        if j + 1 < N:
            # coupling to the right: cell index = order[j] unless this is the last copy of k, then old (k,k+1)
            oj = order[j]
            if order[j + 1] == oj:         # between two copies of k: use bulk coupling (k,k+1)
                C = blk(P, oj, oj + 1)
                CT = blk(P, oj + 1, oj)
            else:
                C = blk(P, oj, order[j + 1]); CT = blk(P, order[j + 1], oj)
            M[2 * j:2 * j + 2, 2 * j + 2:2 * j + 4] = C
            M[2 * j + 2:2 * j + 4, 2 * j:2 * j + 2] = CT
    return M

def rh(n):
    X = rng.standard_normal((n, n)); return (X + X.T) / 2

best = -1; best_st = None; rows = []; best_full = -1
t0 = time.time()
for k in [int(x) for x in a.ks.split(",")]:
    for sg in [float(x) for x in a.sigmas.split(",")]:
        for tr in range(a.ntry if sg > 0 else 1):
            nn = r + 2 * a.m
            A = []; B = []
            for P in Ar:
                M = build(P, k, a.m)
                A.append(bell._pos_projector(M - 0.5 * np.eye(nn) + (sg * rh(nn) if sg > 0 else 0)))
            for P in Br:
                M = build(P, k, a.m)
                B.append(bell._pos_projector(M - 0.5 * np.eye(nn) + (sg * rh(nn) if sg > 0 else 0)))
            v, stt = run_seesaw(A, B)
            sv = np.linalg.svd(np.asarray(stt["psi"]).reshape(nn, nn), compute_uv=False)
            act = int((sv > 1e-4).sum())
            rows.append((v, k, sg, act))
            print(f"k={k} sigma={sg} try {tr}: {v:.10f} active {act}  ({time.time()-t0:.0f}s)", flush=True)
            if v > best:
                best, best_st = v, stt
                if a.save:
                    pickle.dump({"c": a.target, "st": best_st, "v": best}, open(a.save, "wb"))
            if act == nn and v > best_full:
                best_full = v
                if a.save:
                    pickle.dump({"c": a.target, "st": stt, "v": v}, open(a.save + ".full.pkl", "wb"))
print(f"best {best:.10f}")
if a.log:
    open(a.log, "a").write(f"{a.target} insert m={a.m} ks={a.ks} sigmas={a.sigmas} seed {a.seed}: best {best:.10f}\n")
