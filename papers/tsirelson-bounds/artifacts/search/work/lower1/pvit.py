"""PV 2010 Schmidt-basis iteration (n x n matrices only, ~30x faster than the full seesaw), real field.
    python3 work/lower1/pvit.py TARGET n RESTARTS SEED [--tmax s] [--iters N] [--save out.pkl] [--minrank r]
State: psi = sum_i lam_i |ii>, A_x, B_y real symmetric projectors n x n (Schmidt basis).
Saves best as pkl {"c","st":{"A","B","psi" (n*n vector)},"v"}.  Importable: iterate(c, A, B, lam)."""
import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
import sys, time, argparse, pickle
sys.path.insert(0, "qlib")
import numpy as np
import bell


def pos_proj_batch(X):
    X = (X + np.swapaxes(X, -1, -2)) / 2
    w, V = np.linalg.eigh(X)
    Vp = V * (w > 0)[..., None, :]
    return np.einsum("...ik,...jk->...ij", Vp, V)


def kmat(c, A, B):
    n = A.shape[-1]
    I = np.eye(n)
    K = -c["b0"] * I + np.einsum("x,xij->ij", c["bA"], A) * I + np.einsum("y,yij->ij", c["bB"], B) * I
    K = K + np.einsum("xy,xij,yij->ij", c["bAB"], A, B)
    return K


def iterate(c, A, B, lam, iters=3000, tol=1e-13):
    """A: (mA,n,n), B: (mB,n,n), lam: (n,). Returns value, A, B, lam."""
    n = len(lam)
    I = np.eye(n)
    old = -np.inf
    val = old
    for it in range(iters):
        L = np.outer(lam, lam)
        XA = (c["bA"][:, None, None] * I + np.einsum("xy,yij->xij", c["bAB"], B)) * L
        A = pos_proj_batch(XA)
        XB = (c["bB"][:, None, None] * I + np.einsum("xy,xij->yij", c["bAB"], A)) * L
        B = pos_proj_batch(XB)
        K = kmat(c, A, B)
        w, V = np.linalg.eigh((K + K.T) / 2)
        val = w[-1]
        lam = V[:, -1]
        if val - old < tol:
            break
        old = val
    return val, A, B, lam


def to_state(A, B, lam):
    n = len(lam)
    psi = np.zeros(n * n)
    psi[np.arange(n) * n + np.arange(n)] = lam
    return {"A": [a.copy() for a in A], "B": [b.copy() for b in B], "psi": psi}


def from_state(st):
    """full state -> Schmidt basis arrays (A,B,lam) real, active rank kept at full n."""
    A0 = [np.array(P, float) for P in st["A"]]; B0 = [np.array(P, float) for P in st["B"]]
    n1, n2 = A0[0].shape[0], B0[0].shape[0]
    R = np.asarray(st["psi"], float).reshape(n1, n2)
    R = R / np.linalg.norm(R)
    U, s, Vh = np.linalg.svd(R)
    k = min(n1, n2)
    A = np.array([(U.T @ P @ U)[:k, :k] for P in A0]); B = np.array([(Vh @ P @ Vh.T)[:k, :k] for P in B0])
    return A, B, s[:k]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("target"); ap.add_argument("n", type=int); ap.add_argument("restarts", type=int)
    ap.add_argument("seed", type=int)
    ap.add_argument("--tmax", type=float, default=120)
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--minrank", type=int, default=0)
    ap.add_argument("--save", default=None)
    a = ap.parse_args()
    c = bell.load(a.target)
    rec = float(bell.record(a.target)["lower"])
    rng = np.random.default_rng(a.seed)
    n = a.n
    best, bst = -1, None
    vals = []
    t0 = time.time()
    for k in range(a.restarts):
        Bm = []
        for y in range(c["mB"]):
            r = rng.integers(1, n)
            D = np.diag((np.arange(n) < r).astype(float))
            Q = np.linalg.qr(rng.standard_normal((n, n)))[0]
            Bm.append(Q @ D @ Q.T)
        lam = rng.uniform(0.2, 1.0, n); lam /= np.linalg.norm(lam)
        A0 = np.zeros((c["mA"], n, n))
        val, A, B, lam2 = iterate(c, A0, np.array(Bm), lam, a.iters)
        vals.append(val)
        eff = int((lam2 ** 2 > 1e-6).sum())
        if val > best and eff >= a.minrank:
            best, bst = val, (A, B, lam2)
        if time.time() - t0 > a.tmax:
            break
    vals = np.array(vals)
    u = np.unique(np.round(vals, 7))[::-1][:6]
    A, B, lam2 = bst
    st = to_state(A, B, lam2)
    print(f"{a.target} PV n={a.n}: {len(vals)} restarts {time.time()-t0:.0f}s best {best:.10f} (vs PV rec {best-rec:+.3e}) fullvalue {bell.value(c, st):.10f}; top", " ".join(f"{x:.7f}" for x in u), flush=True)
    print("lam", np.round(lam2, 4), flush=True)
    if a.save:
        pickle.dump({"c": a.target, "st": st, "v": best}, open(a.save, "wb"))
