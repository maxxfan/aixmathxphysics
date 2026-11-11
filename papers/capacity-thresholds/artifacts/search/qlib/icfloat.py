#!/usr/bin/env python3
"""Fast floating-point coherent information of permutation-symmetric inputs to n uses of a qubit
Pauli channel (SEARCH TOOL ONLY: nothing here is rigorous; the verifier decides).

Input: rows V (r x (n+1) complex), Dicke basis |D_a> (a = number of ones), rho = sum_r |v_r><v_r|.
Channel: Pauli probabilities (pI, pX, pY, pZ); helpers depol(p), twopauli(p), indep(p).

Block formula (Schur-Weyl; see PROBLEM.md, "Block formula"): the output of N^{(x)n} on a symmetric
input decomposes over two-row sectors j = n/2, n/2-1, ... with multiplicity
m_j = C(n,k) - C(n,k-1), k = n/2 - j, and on the (2j+1)-dimensional spin-j space it is

    sigma_j = coefficient extraction of det(Z)^k (x^dag Z y)^{2j},   Z = N(|u><w|),

i.e. "contract k pairs of qubits with the channel's determinant form, then apply
P_sym N^{(x)2j} P_sym".  For the depolarizing channel this is
    sigma_j = (2p(1-2p))^k  Phi_{2j}( Tr_{2k} rho ),
    Phi_d(Y) = sum_s C(d,s) (1-4p)^{d-s} (2p)^s  E^s Tr^s (Y)        (E = Tr_1^dagger),
evaluated by Horner in O(d^3).  General Pauli channels use the O(d^5) dynamic programme below.

Coherent information (bits) of the normalized state:
    I_c = (1/Z) sum_j [ H(m_j sigma_j) - H(m_j tau_j) ],   Z = sum_r |v_r|^2,
tau_j the joint (output x reference) block, H(X) = -sum eig log2 eig.

usage (library):  from icfloat import coh_info, depol;  coh_info(V, depol(0.0649))
"""
import math

import numpy as np


def depol(p):
    return (1 - 3 * p, p, p, p)


def twopauli(p):
    return (1 - 2 * p, p, 0.0, p)


def indep(p):
    return ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p))


CHANNELS = {"depolarizing": depol, "two_pauli": twopauli, "independent_xz": indep}


def mult(n, k):
    """dimension of the S_n irrep (n-k, k)"""
    return math.comb(n, k) - (math.comb(n, k - 1) if k > 0 else 0)


def tr1(T):
    """partial trace of one qubit, Dicke basis, degree N -> N-1; T has shape (..., N+1, N+1)"""
    N = T.shape[-1] - 1
    c = np.arange(N)
    s0 = np.sqrt(N - c)
    s1 = np.sqrt(c + 1)
    return (s0[:, None] * s0[None, :] * T[..., :N, :N] + s1[:, None] * s1[None, :] * T[..., 1:, 1:]) / N


def emb1(Z):
    """E = Tr_1^dagger: Sym^{N-1} -> Sym^N, Z |-> P_sym (Z (x) I) P_sym"""
    N = Z.shape[-1]
    a = np.arange(N + 1)
    out = np.zeros(Z.shape[:-2] + (N + 1, N + 1), dtype=Z.dtype)
    s0 = np.sqrt(N - a)
    s1 = np.sqrt(a)
    out[..., :N, :N] += (s0[:N, None] * s0[None, :N]) * Z
    out[..., 1:, 1:] += (s1[1:, None] * s1[None, 1:]) * Z
    return out / N


def entropy_bits(M):
    w = np.linalg.eigvalsh((M + M.conj().T) / 2)
    w = w[w > 0]
    return float(-(w * np.log2(w)).sum())


def joint(Y):
    """Y[r, r', a, b] -> (r(d+1)) x (r(d+1)) matrix"""
    r, _, d1, _ = Y.shape
    return Y.transpose(0, 2, 1, 3).reshape(r * d1, r * d1)


def blocks_depol(V, p):
    """yield (j, m_j, M) with M[r, r', :, :] = m_j * sigma-block of |v_r><v_r'| (depolarizing)"""
    V = np.asarray(V, dtype=complex)
    n = V.shape[1] - 1
    Y = np.einsum("ia,jb->ijab", V, V.conj())
    chain = [Y]
    for _ in range(n):
        chain.append(tr1(chain[-1]))
    A = 2 * p * (1 - 2 * p)
    beta = 1 - 4 * p
    for k in range(n // 2 + 1):
        d = n - 2 * k
        mj = mult(n, k)
        cs = [math.comb(d, s) * beta ** (d - s) * (2 * p) ** s for s in range(d + 1)]
        R = cs[d] * chain[2 * k + d]
        for s in range(d - 1, -1, -1):
            R = cs[s] * chain[2 * k + s] + emb1(R)
        yield d / 2, mj, mj * A ** k * R


def blocks_pauli(V, probs):
    """general Pauli channel (pI, pX, pY, pZ): O(d^5) per block; same output as blocks_depol"""
    V = np.asarray(V, dtype=complex)
    n = V.shape[1] - 1
    pI, pX, pY, pZ = probs
    cA, cB, cC, cE = pI + pZ, pX + pY, pI - pZ, pX - pY
    sq = np.sqrt(np.array([math.comb(n, a) for a in range(n + 1)], dtype=float))
    Y = np.einsum("ia,jb->ijab", V / sq, (V / sq).conj())  # monomial-basis coefficients
    det = {(0, 0): cA * cB, (2, 2): cA * cB, (1, 1): cA * cA + cB * cB - cC * cC - cE * cE,
           (0, 2): -cC * cE, (2, 0): -cC * cE}
    kern = {(0, 0): {(0, 0): cA, (1, 1): cB}, (1, 1): {(0, 0): cB, (1, 1): cA},
            (0, 1): {(0, 1): cC, (1, 0): cE}, (1, 0): {(1, 0): cC, (0, 1): cE}}
    Yk = Y
    for k in range(n // 2 + 1):
        if k > 0:
            m = Yk.shape[-1] - 2
            nxt = np.zeros(Yk.shape[:2] + (m, m), dtype=complex)
            for (al, be), w in det.items():
                if w:
                    nxt += w * Yk[..., al:al + m, be:be + m]
            Yk = nxt
        d = n - 2 * k
        # G[(alpha,beta)] arrays of shape (r, r, d-s+1, d-s+1); start s = d
        G = {(al, be): Yk[:, :, al:al + 1, be:be + 1] for al in range(d + 1) for be in range(d + 1)}
        for s in range(d - 1, -1, -1):
            L = d - s + 1
            newG = {}
            for al in range(s + 1):
                for be in range(s + 1):
                    acc = np.zeros(Yk.shape[:2] + (L, L), dtype=complex)
                    for (da, db), kd in kern.items():
                        X = G[(al + da, be + db)]
                        for (dm, dm2), w in kd.items():
                            if w:
                                acc[..., dm:dm + L - 1, dm2:dm2 + L - 1] += w * X
                    newG[(al, be)] = acc
            G = newG
        out = G[(0, 0)]
        cd = np.sqrt(np.array([math.comb(d, m) for m in range(d + 1)], dtype=float))
        out = out / (cd[:, None] * cd[None, :])  # coefficient sums -> Dicke basis
        mj = mult(n, k)
        yield d / 2, mj, mj * out


def coh_info(V, probs, method=None):
    """I_c in bits per block of n channel uses (not per channel use)"""
    V = np.asarray(V, dtype=complex)
    Z = float((np.abs(V) ** 2).sum())
    pI, pX, pY, pZ = probs
    if method is None:
        method = "depol" if (pX == pY == pZ) else "pauli"
    gen = blocks_depol(V, pX) if method == "depol" else blocks_pauli(V, probs)
    tot = 0.0
    for _, _, M in gen:
        sig = np.einsum("iiab->ab", M)
        tot += entropy_bits(sig) - entropy_bits(joint(M))
    return tot / Z


def brute(V, probs):
    """independent brute force on the full 2^n space (n <= 7)"""
    V = np.asarray(V, dtype=complex)
    r, n1 = V.shape
    n = n1 - 1
    dim = 2 ** n
    pops = np.array([bin(x).count("1") for x in range(dim)])
    psi = np.zeros((dim, r), dtype=complex)
    for x in range(dim):
        a = pops[x]
        psi[x, :] = V[:, a] / math.sqrt(math.comb(n, a))
    rho = np.einsum("xi,yj->xiyj", psi, psi.conj()).reshape(dim * r, dim * r)
    P = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    for q in range(n):
        new = np.zeros_like(rho)
        for pr, Pm in zip(probs, P):
            if pr == 0:
                continue
            op = np.kron(np.kron(np.eye(2 ** q), Pm), np.eye(2 ** (n - q - 1)))
            op = np.kron(op, np.eye(r))
            new += pr * op @ rho @ op.conj().T
        rho = new
    rhoB = np.einsum("xiyi->xy", rho.reshape(dim, r, dim, r))
    Z = np.trace(rhoB).real
    return (entropy_bits(rhoB / Z) - entropy_bits(rho / Z))


def load_rows(path):
    """read rows from a raw payload (.c128: little-endian complex128, shape (r, n+1), r given by
    the file size and n) or JSON {"rows": [[[re, im], ...], ...]}"""
    import json
    import pathlib
    p = pathlib.Path(path)
    if p.suffix == ".json":
        d = json.loads(p.read_text())
        return np.array([[complex(x, y) for x, y in row] for row in d["rows"]])
    raw = np.frombuffer(p.read_bytes(), dtype="<c16")
    return raw.reshape(2, -1)
