"""Build candidate files for the p05 verifier (numpy + cvxpy; run with the repo .venv python).

    import sys; sys.path.insert(0, "qlib"); import bell, certs
    c = bell.load("J97")
    v, st = bell.seesaw(c, 2, "real", restarts=50)
    cand = certs.strategy_candidate(c, st, field="real")          # exact rational strategy
    certs.save(cand, "work/me/J97_lower.json")

    lam, Y, W = certs.solve_dual(c, level=2)                        # numerical NPA dual
    cand = certs.sos_candidate(c, W, Y, eps=1e-7)                  # rounded, identity fixed exactly
    certs.save(cand, "work/me/J97_upper_L2.json")

Then check it with the stdlib verifier:  python3 verifier/verify.py work/me/J97_upper_L2.json
"""
import json
from fractions import Fraction

import numpy as np

import bell


# ----------------------------------------------------------------------------- strategies

def _span(P, tol=1e-6):
    w, V = np.linalg.eigh((P + P.conj().T) / 2)
    return [V[:, k] for k in range(len(w)) if w[k] > 0.5]


def _q(x, den):
    return str(Fraction(int(round(float(x) * den)), den))


def _entry(z, den, field):
    if field == "real":
        return _q(np.real(z), den)
    return [_q(np.real(z), den), _q(np.imag(z), den)]


def strategy_candidate(c, st, field=None, den=10 ** 12, **meta):
    """Rationalize a numerical strategy: projectors -> rounded eigenvectors spanning them (the
    verifier rebuilds the exact projector from the span), psi -> rounded entries."""
    A, B, psi = st["A"], st["B"], np.asarray(st["psi"])
    if field is None:
        field = "complex" if any(np.iscomplexobj(M) and np.abs(np.imag(M)).max() > 1e-12 for M in A + B + [psi]) else "real"
    dA, dB = A[0].shape[0], B[0].shape[0]
    psi = psi / np.linalg.norm(psi)
    if field == "real":
        psi = np.real(psi * np.exp(-1j * np.angle(psi[np.argmax(np.abs(psi))]))) if np.iscomplexobj(psi) else psi
    cand = {"type": "strategy", "inequality": c["name"], "field": field, "dA": dA, "dB": dB,
            "A": [[[_entry(z, den, field) for z in v] for v in _span(P)] for P in A],
            "B": [[[_entry(z, den, field) for z in v] for v in _span(Q)] for Q in B],
            "psi": [_entry(z, den, field) for z in psi]}
    cand.update(meta)
    return cand


# ----------------------------------------------------------------------------- SOS / NPA dual

def word_str(w):
    return "".join(f"A{x + 1}" for x in w[0]) + "".join(f"B{y + 1}" for y in w[1])


def classes(W):
    """Map each pair (i, j) to its adjoint-merged reduced word; returns (cls index matrix, class list)."""
    n = len(W)
    idx, cl = {}, np.empty((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            m = bell.canon_real(bell.product(W[i], W[j]))
            cl[i, j] = idx.setdefault(m, len(idx))
    inv = [None] * len(idx)
    for m, k in idx.items():
        inv[k] = m
    return cl, inv


def targets(c, inv):
    """t(M) for each merged class: sum over m in M of the coefficient of m in  -W  (identity excluded)."""
    t = np.zeros(len(inv))
    for k, m in enumerate(inv):
        a, b = m
        if len(a) == 1 and len(b) == 0:
            t[k] = -c["bA"][a[0]]
        elif len(a) == 0 and len(b) == 1:
            t[k] = -c["bB"][b[0]]
        elif len(a) == 1 and len(b) == 1:
            t[k] = -c["bAB"][a[0], b[0]]
    return t


def solve_dual(c, level=2, solver="CLARABEL", verbose=False, **kw):
    if solver == "CLARABEL":
        kw = dict(dict(tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=400), **kw)
    """min lambda s.t. Y >= 0 and sum_{(i,j) in M} Y_ij = t(M) for every class M != 1, lambda = Y_00.
    Returns (lambda - b0, Y, W)."""
    import cvxpy as cp
    import scipy.sparse as sp
    W = bell.words(c, level)
    if W[0] != ((), ()):
        raise ValueError("word list must start with the identity")
    n = len(W)
    cl, inv = classes(W)
    one = cl[0, 0]
    t = targets(c, inv)
    rows = cl.reshape(-1)
    A = sp.csr_matrix((np.ones(n * n), (rows, np.arange(n * n))), shape=(len(inv), n * n))
    keep = np.array([k for k in range(len(inv)) if k != one])
    Y = cp.Variable((n, n), PSD=True)
    prob = cp.Problem(cp.Minimize(Y[0, 0]), [A[keep] @ cp.vec(Y, order="C") == t[keep]])
    prob.solve(solver=solver, verbose=verbose, **kw)
    return prob.value - c["b0"], np.array(Y.value), W


def sos_candidate(c, W, Y, eps=1e-8, den=2 ** 44, **meta):
    """Round Y + eps*I to integers over 2*den, then repair the identity exactly by moving each class
    residual onto one representative entry (a diagonal one if the class has one)."""
    n = len(W)
    D2 = 2 * den
    Yr = np.rint((np.asarray(Y) + eps * np.eye(n)) * den).astype(object) * 2
    Yr = [[int(Yr[i, j]) if i <= j else 0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            Yr[i][j] = Yr[j][i]
    cl, inv = classes(W)
    t = targets(c, inv)
    one = cl[0, 0]
    members = [[] for _ in inv]
    for i in range(n):
        for j in range(n):
            members[cl[i, j]].append((i, j))
    for k, mem in enumerate(members):
        if k == one:
            continue
        s = sum(Yr[i][j] for i, j in mem)
        want = int(round(t[k])) * D2
        r = want - s
        if r == 0:
            continue
        diag = [(i, j) for i, j in mem if i == j]
        if diag:
            i, _ = diag[0]
            Yr[i][i] += r
        else:
            i, j = min(mem, key=lambda p: (min(p), max(p)))
            if r % 2:
                raise ArithmeticError("odd residual; increase den")
            Yr[i][j] += r // 2
            Yr[j][i] += r // 2
    lam = Fraction(Yr[0][0], D2)
    cand = {"type": "sos", "inequality": c["name"], "words": [word_str(w) for w in W],
            "Y": {"denominator": str(D2), "upper": [[str(Yr[i][j]) for j in range(i, n)] for i in range(n)]},
            "lambda": str(lam)}
    cand.update(meta)
    return cand


def sos_candidate_auto(c, W, Y, eps_list=None, den=2 ** 44, **meta):
    """Try increasing eps until the rounded, repaired matrix is numerically positive definite with a
    small margin; returns (candidate, eps, min eigenvalue).  The verifier then needs tau = 0."""
    n = len(W)
    for eps in eps_list or [1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5]:
        cand = sos_candidate(c, W, Y, eps=eps, den=den, **meta)
        D2 = int(cand["Y"]["denominator"])
        M = np.zeros((n, n))
        for i, row in enumerate(cand["Y"]["upper"]):
            for k, e in enumerate(row):
                M[i, i + k] = M[i + k, i] = int(e) / D2
        lmin = np.linalg.eigvalsh(M)[0]
        if lmin > 1e-12 * max(1.0, np.abs(M).max()) * n:
            return cand, eps, lmin
    return cand, eps, lmin


def save(cand, path):
    with open(path, "w") as fh:
        json.dump(cand, fh)
    return path
