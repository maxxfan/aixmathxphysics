#!/usr/bin/env python3
"""Full 2^n density-matrix reference for I_c (auditor, independent of verifier/, qlib/, work/lead/).

Builds the n-qubit Dicke vectors explicitly, applies the Pauli channel qubit by qubit as a
mixture of Pauli conjugations on the system (and on the system of the purification), and takes
entropies of the full output matrices.  Only usable for n <= 10 or so.
"""
import numpy as np
from math import comb

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def dicke_matrix(n):
    """(2^n, n+1) matrix whose column a is |D^n_a> (a ones, equal weights)."""
    N = 1 << n
    w = np.array([bin(x).count('1') for x in range(N)])
    D = np.zeros((N, n + 1))
    for a in range(n + 1):
        D[w == a, a] = 1.0 / np.sqrt(comb(n, a))
    return D


def apply_pauli_channel(T, n, rdim, probs):
    """Apply the qubit Pauli channel with probs (pI,pX,pY,pZ) to every one of the n qubits of
    the operator T on (C^2)^{(x)n} (x) C^rdim (system first, qubit 0 most significant)."""
    Ps = [PAULI['I'], PAULI['X'], PAULI['Y'], PAULI['Z']]
    dim = (1 << n) * rdim
    T = T.reshape(dim, dim)
    for k in range(n):
        left = 1 << k
        right = (1 << (n - k - 1)) * rdim
        T6 = T.reshape(left, 2, right, left, 2, right)
        out = np.zeros_like(T6)
        for p, P in zip(probs, Ps):
            if p == 0:
                continue
            out += p * np.einsum('as,isjktl,bt->iajkbl', P, T6, P.conj(), optimize=True)
        T = out.reshape(dim, dim)
    return T


def entropy_bits(M):
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    ev = ev[ev > 1e-300]
    return float(-(ev * np.log2(ev)).sum()), ev


def ic_brute(V, probs):
    """V: (r, n+1) complex Dicke coefficients (unnormalised).  probs: floats (pI,pX,pY,pZ)."""
    V = np.asarray(V, dtype=complex)
    r, n1 = V.shape
    n = n1 - 1
    Z = float(np.sum(np.abs(V) ** 2))
    D = dicke_matrix(n)
    Vf = (D @ V.T) / np.sqrt(Z)            # (2^n, r): column q is |v_q>
    sigma = Vf @ Vf.conj().T
    psi = Vf.reshape(-1)                   # system (x) reference, index x*r+q
    tau = np.outer(psi, psi.conj())
    sig_out = apply_pauli_channel(sigma, n, 1, probs)
    tau_out = apply_pauli_channel(tau, n, r, probs)
    S_B, evB = entropy_bits(sig_out)
    S_BR, evBR = entropy_bits(tau_out)
    return S_B - S_BR, evB, evBR


if __name__ == '__main__':
    import sys, json
    rng = np.random.default_rng(1)
    n, r = 5, 2
    V = rng.normal(size=(r, n + 1)) + 1j * rng.normal(size=(r, n + 1))
    p = 0.05
    print(ic_brute(V, (1 - 3 * p, p, p, p))[0])
