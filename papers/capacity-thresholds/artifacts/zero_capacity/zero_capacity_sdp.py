"""Where Q = 0 is forced by antidegradability, for each target family.

A channel is antidegradable iff its Choi state rho_AB has an extension rho_AB1B2 symmetric in the two
outputs. For each family we solve  max t  s.t.  W - t*1 >= 0,  Tr_B2 W = Tr_B1 W = rho  (W on A B1 B2)
and bisect p on t* >= -tol. The 2-Pauli Choi state has a zero eigenvalue (p_Y = 0), so any extension is
singular and t* is 0 on the antidegradable side; hence the tolerance.

Compares with the closed form sqrt(p_I) <= sqrt(p_X) + sqrt(p_Y) + sqrt(p_Z).
Needs cvxpy with CLARABEL (repo .venv).  Output (2026-10-27):
  depolarizing    SDP 0.083333  closed form 0.083333   (1/12)
  independent_xz  SDP 0.146446  closed form 0.146447   ((1 - 1/sqrt 2)/2)
  two_pauli       SDP 0.166545  closed form 0.166667   (1/6; SDP side shifted by the 1e-7 tolerance)
"""
import cvxpy as cp
import numpy as np

I2 = np.eye(2)
PAULI = (I2, np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
PHI = np.array([1, 0, 0, 1]) / np.sqrt(2)
FAM = {
    "depolarizing": lambda p: (1 - 3 * p, p, p, p),
    "independent_xz": lambda p: ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p)),
    "two_pauli": lambda p: (1 - 2 * p, p, 0.0, p),
}
SWAP = np.zeros((8, 8))  # A B1 B2 -> A B2 B1
for a in range(2):
    for b in range(2):
        for c in range(2):
            SWAP[a * 4 + c * 2 + b, a * 4 + b * 2 + c] = 1


def choi(probs):
    out = np.zeros((4, 4), complex)
    for w, P in zip(probs, PAULI):
        v = np.kron(I2, P) @ PHI
        out += w * np.outer(v, v.conj())
    return out


def min_eig_of_best_extension(probs):
    rho = choi(probs)
    W, t = cp.Variable((8, 8), hermitian=True), cp.Variable()
    cons = [W - t * np.eye(8) >> 0,
            cp.partial_trace(W, [2, 2, 2], axis=2) == rho,
            cp.partial_trace(SWAP @ W @ SWAP.T, [2, 2, 2], axis=2) == rho]
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver="CLARABEL")
    assert prob.status == "optimal", prob.status
    return t.value


def bisect(test, lo, hi, steps=40):
    for _ in range(steps):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if test(mid) else (mid, hi)
    return hi


if __name__ == "__main__":
    for fam, f in FAM.items():
        sdp = bisect(lambda p: min_eig_of_best_extension(f(p)) >= -1e-7, 0.05, 0.3, 30)
        closed = bisect(lambda p: np.sqrt(f(p)[0]) <= sum(np.sqrt(x) for x in f(p)[1:]), 0.0, 0.3, 60)
        print(f"{fam:15s} SDP {sdp:.6f}  closed form {closed:.6f}")
