"""Numerical helpers for p05 (numpy; run with the repo .venv python).  Not part of the verifier.

    import sys; sys.path.insert(0, "qlib"); import bell
    c = bell.load("J73")                       # coefficients as numpy arrays
    v, st = bell.seesaw(c, d=4, field="real", restarts=200, seed=1)
    bell.value(c, st)                          # recompute the violation of a strategy
    ub = bell.npa_upper(c, level="1+ab")       # numerical NPA bound (cvxpy + SCS/Clarabel)

A strategy is a dict {"A": [P_1..P_mA], "B": [Q_1..Q_mB], "psi": vector of length dA*dB}
with P_x, Q_y orthogonal projectors (outcome "1") and psi indexed a*dB + b.
Conventions follow data/make_records.py (Pal-Vertesi 2009, Eq. (2)).
"""
import itertools
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalog.json"
RECORDS = ROOT / "data" / "records.json"


def load(name):
    c = json.loads(CATALOG.read_text())["inequalities"][name]
    return {"name": name, "bA": np.array(c["bA"], float), "bB": np.array(c["bB"], float),
            "bAB": np.array(c["bAB"], float), "b0": c["b0"], "mA": c["mA"], "mB": c["mB"]}


def record(name):
    return json.loads(RECORDS.read_text())["records"][name]


def local_bound(c):
    best = -np.inf
    for a in itertools.product((0, 1), repeat=c["mA"]):
        a = np.array(a, float)
        s = c["bB"] + a @ c["bAB"]
        best = max(best, a @ c["bA"] + np.clip(s, 0, None).sum())
    return best


def bell_operator(c, A, B):
    dA, dB = A[0].shape[0], B[0].shape[0]
    IA, IB = np.eye(dA), np.eye(dB)
    W = sum(c["bA"][x] * np.kron(A[x], IB) for x in range(c["mA"]))
    W = W + sum(c["bB"][y] * np.kron(IA, B[y]) for y in range(c["mB"]))
    W = W + sum(c["bAB"][x, y] * np.kron(A[x], B[y]) for x in range(c["mA"]) for y in range(c["mB"]))
    return W


def value(c, st):
    psi = np.asarray(st["psi"])
    psi = psi / np.linalg.norm(psi)
    W = bell_operator(c, st["A"], st["B"])
    return float(np.real(np.vdot(psi, W @ psi))) - c["b0"]


def _pos_projector(M):
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    V = V[:, w > 0]
    return V @ V.conj().T


def _rand_proj(d, rng, field):
    X = rng.standard_normal((d, d))
    if field == "complex":
        X = X + 1j * rng.standard_normal((d, d))
    H = X + X.conj().T
    return _pos_projector(H)


def seesaw(c, d, field="real", restarts=50, iters=300, seed=0, tol=1e-13, dB=None):
    """Alternating optimization: state <- top eigenvector; each projector <- positive part of its
    effective operator.  Returns (best violation, strategy)."""
    rng = np.random.default_rng(seed)
    dA, dB = d, dB or d
    best, best_st = -np.inf, None
    for _ in range(restarts):
        A = [_rand_proj(dA, rng, field) for _ in range(c["mA"])]
        B = [_rand_proj(dB, rng, field) for _ in range(c["mB"])]
        old = -np.inf
        for _ in range(iters):
            W = bell_operator(c, A, B)
            w, V = np.linalg.eigh(W)
            psi = V[:, -1]
            R = psi.reshape(dA, dB)            # psi = sum R[a,b] |a>|b>
            # <psi| A (x) O |psi> = Tr(A R O^T R^dag);  <psi| O (x) B |psi> = Tr(B R^T O^T conj(R))
            for x in range(c["mA"]):
                Ox = c["bA"][x] * np.eye(dB) + sum(c["bAB"][x, y] * B[y] for y in range(c["mB"]))
                A[x] = _pos_projector(R @ Ox.T @ R.conj().T)
            for y in range(c["mB"]):
                Oy = c["bB"][y] * np.eye(dA) + sum(c["bAB"][x, y] * A[x] for x in range(c["mA"]))
                B[y] = _pos_projector(R.T @ Oy.T @ R.conj())
            val = w[-1]
            if val - old < tol:
                break
            old = val
        st = {"A": A, "B": B, "psi": np.linalg.eigh(bell_operator(c, A, B))[1][:, -1]}
        v = value(c, st)
        if v > best:
            best, best_st = v, st
    return best, best_st


# ----------------------------------------------------------------------------- NPA (numerical)

def reduce_word(w):
    """w: tuple of (party, setting); party 0 = Alice, 1 = Bob.  Alice letters commute with Bob's;
    projectors are idempotent.  Normal form: (alice tuple, bob tuple) with no equal neighbours."""
    a = [s for p, s in w if p == 0]
    b = [s for p, s in w if p == 1]
    return (_collapse(a), _collapse(b))


def _collapse(seq):
    out = []
    for s in seq:
        if out and out[-1] == s:
            continue
        out.append(s)
    return tuple(out)


def words(c, level):
    """Monomial list for an NPA level: an int L (all words of length <= L) or a string like
    "1+ab", "2+aab", "2+aab+abb", "3+aaab+aabb+abbb"."""
    mA, mB = c["mA"], c["mB"]

    def all_words(na, nb):
        out = []
        for a in itertools.product(range(mA), repeat=na):
            if any(a[i] == a[i + 1] for i in range(na - 1)):
                continue
            for b in itertools.product(range(mB), repeat=nb):
                if any(b[i] == b[i + 1] for i in range(nb - 1)):
                    continue
                out.append((a, b))
        return out
    if isinstance(level, int) or str(level).isdigit():
        L = int(level)
        shapes = [(i, j) for i in range(L + 1) for j in range(L + 1 - i)]
    else:
        base, *extra = str(level).split("+")
        L = int(base)
        shapes = [(i, j) for i in range(L + 1) for j in range(L + 1 - i)]
        for e in extra:
            shapes.append((e.count("a"), e.count("b")))
    seen, out = set(), []
    for na, nb in shapes:
        for w in all_words(na, nb):
            if w not in seen:
                seen.add(w)
                out.append(w)
    return out


def product(u, v):
    """(u^dagger) v in normal form; u, v are (alice tuple, bob tuple)."""
    return (_collapse(tuple(reversed(u[0])) + v[0]), _collapse(tuple(reversed(u[1])) + v[1]))


def canon_real(m):
    """Real moment matrices: identify m with its adjoint (reverse both tuples)."""
    r = (tuple(reversed(m[0])), tuple(reversed(m[1])))
    return min(m, r)


def npa_upper(c, level=2, solver="SCS", verbose=False, eps=1e-9, return_data=False):
    """Numerical NPA bound (real moment matrix).  Returns the optimal value (violation)."""
    import cvxpy as cp
    W = words(c, level)
    n = len(W)
    index = {}
    pos = np.empty((n, n), dtype=int)
    for i in range(n):
        for j in range(i, n):
            m = canon_real(product(W[i], W[j]))
            k = index.setdefault(m, len(index))
            pos[i, j] = pos[j, i] = k
    y = cp.Variable(len(index))
    G = cp.Variable((n, n), symmetric=True)
    cons = [G >> 0]
    # tie entries to moments
    for i in range(n):
        for j in range(i, n):
            cons.append(G[i, j] == y[pos[i, j]])
    one = index[((), ())]
    cons.append(y[one] == 1)
    obj = 0
    for x in range(c["mA"]):
        obj += c["bA"][x] * y[index[canon_real(((x,), ()))]]
    for yy in range(c["mB"]):
        obj += c["bB"][yy] * y[index[canon_real(((), (yy,)))]]
    for x in range(c["mA"]):
        for yy in range(c["mB"]):
            obj += c["bAB"][x, yy] * y[index[canon_real(((x,), (yy,)))]]
    prob = cp.Problem(cp.Maximize(obj), cons)
    kw = {"eps": eps, "max_iters": 200000} if solver == "SCS" else {}
    prob.solve(solver=solver, verbose=verbose, **kw)
    if return_data:
        return prob.value - c["b0"], {"words": W, "index": index, "pos": pos, "problem": prob, "G": G, "y": y}
    return prob.value - c["b0"]
