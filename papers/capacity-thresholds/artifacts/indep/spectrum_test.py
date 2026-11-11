#!/usr/bin/env python3
"""Structural test of the Schur-block superoperators S_P = sum_i P_i . P_i at moderate n.

On S_n-invariant operators the Pauli-weight superoperator W = (3n - S_X - S_Y - S_Z)/4 has
eigenvalue w (= Pauli weight) with multiplicity (w+1)(w+2)/2 (symmetric monomials in X,Y,Z of
degree w), and S_Z has eigenvalue n-2u (u = number of X/Y letters) with multiplicity
(u+1)(n-u+1).  Both are checked per conserved Delta = m-m' sector with dense eigenvalues.
"""
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep_ic as I  # noqa: E402


def check(n):
    L = I.Layout(n)
    SX = I.build_superop(L, I.PAULI['X'], I.PAULI['X'])
    SY = I.build_superop(L, I.PAULI['Y'], I.PAULI['Y'])
    SZ = I.build_superop(L, I.PAULI['Z'], I.PAULI['Z'])
    Se = I.build_superop(L, I.E01, I.E10) + I.build_superop(L, I.E10, I.E01)
    W = (3 * n * np.eye(1) - 0) if False else None
    import scipy.sparse as sp
    Wm = (3 * n * sp.identity(L.N) - SX - SY - SZ) / 4
    delta = L.delta_of()
    cW, cZ = Counter(), Counter()
    errW = errZ = 0.0
    for D in np.unique(delta):
        idx = np.nonzero(delta == D)[0]
        for M, cnt, name in ((Wm, cW, 'W'), (SZ, cZ, 'Z')):
            # Delta-conservation of W and S_Z
            if abs(M[idx]).sum() - abs(M[idx][:, idx]).sum() > 1e-12:
                raise SystemExit(f'{name} leaks out of Delta sector')
            ev = np.linalg.eigvals(M[idx][:, idx].toarray())
            rr = np.round(ev.real)
            e = max(np.max(np.abs(ev - rr)), 0.0)
            if name == 'W':
                errW = max(errW, e)
            else:
                errZ = max(errZ, e)
            cnt.update(rr.astype(int).tolist())
    expW = Counter({w: (w + 1) * (w + 2) // 2 for w in range(n + 1)})
    expZ = Counter({n - 2 * u: (u + 1) * (n - u + 1) for u in range(n + 1)})
    # S_X + S_Y from Pauli form vs from matrix units
    dxy = abs((SX + SY) - 2 * Se).max()
    print(f'n={n} N={L.N}: W spectrum {"OK" if cW == expW else "MISMATCH"} (max |ev-int| {errW:.1e}); '
          f'S_Z spectrum {"OK" if cZ == expZ else "MISMATCH"} ({errZ:.1e}); |S_X+S_Y-2S_e| {dxy:.1e}')
    return cW == expW and cZ == expZ


if __name__ == '__main__':
    ok = True
    for n in [int(x) for x in (sys.argv[1:] or ['2', '3', '7', '12', '20', '25'])]:
        ok &= check(n)
    print('ALL OK' if ok else 'FAIL')
