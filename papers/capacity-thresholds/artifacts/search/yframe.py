"""Y-pole frame for the X-Z and 2-Pauli targets.

I_c(N, U^{(x)n} rho U^{(x)n dag}) = I_c(U^dag N U, rho), and with U = exp(-i pi/4 X) per qubit
(Dicke basis: exp(-i pi/2 J_x)) both channels become z-covariant (p_X = p_Y, family "<fam>@y" in
pauliopt.FAM).  Searching there keeps amplitudes real (about 2.5x faster); candidates for
verifier/verify.py are mapped back with from_frame.
"""
import numpy as np

from ansatz import rotate


def to_frame(V):
    W = rotate(V, "x", -np.pi / 2)
    if np.abs(W.imag).max() > 1e-9 * np.abs(W).max():
        raise ValueError("state is not real in the Y-pole frame")
    return W.real.copy()


def from_frame(W):
    return rotate(W, "x", np.pi / 2)
