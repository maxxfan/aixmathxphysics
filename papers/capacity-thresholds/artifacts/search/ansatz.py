"""mod-3 Gaussian ansatz (structure of the Gaussian-limit optimum of arXiv:2609.39747): row r is
supported on Dicke levels a = r (mod 3), envelope exp(-(a-a0)^2/(4 s^2)), weights q, 1-q."""
import numpy as np


def ansatz(n, a0, s, q, mod=3):
    a = np.arange(n + 1)
    env = np.exp(-(a - a0) ** 2 / (4 * s * s))
    V = np.zeros((2, n + 1))
    for r in range(2):
        v = np.where(a % mod == r, env, 0.0)
        V[r] = v / np.linalg.norm(v) * np.sqrt(q if r == 0 else 1 - q)
    return V + 0j


def spin_ops(n):
    """J_x, J_y, J_z on Sym^n in the Dicke basis (a = number of ones, m = n/2 - a)"""
    a = np.arange(n + 1)
    Jz = np.diag(n / 2 - a).astype(complex)
    Jm = np.zeros((n + 1, n + 1), dtype=complex)  # lowers m: |a> -> |a+1>
    for k in range(n):
        Jm[k + 1, k] = np.sqrt((k + 1) * (n - k))
    Jp = Jm.T.conj()
    return (Jp + Jm) / 2, (Jp - Jm) / 2j, Jz


def rotate(V, axis, theta):
    from scipy.linalg import expm
    n = V.shape[1] - 1
    J = spin_ops(n)["xyz".index(axis)]
    U = expm(-1j * theta * J)
    return V @ U.T
