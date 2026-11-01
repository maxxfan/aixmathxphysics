"""Validate the 7-cycle derivation: for random lines L1, the party-weight histogram of the full 2^16-word stabilizer built from rows
(S_nf basis + graph of f = id) must equal the histogram predicted by the formulas used in c7sym.criterion:
  weight(s + (c^7; f(c))) = #{i: s_i != c} + 1 (c != 0),  weight(s) (c = 0)."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import c7sym as C
import ame  # [release] was: from qlib import ame (ame.py lives in ../lib, path set by c7sym)

rng = random.Random(3)
lines = list(C.rref_subspaces(2))
for trial in range(6):
    L1 = rng.choice(lines)
    a_el = C.span(L1)
    mask = np.ones(len(C.ALL4), dtype=bool)
    for v in L1:
        mask &= (C.ip(v, C.ALL4) == 0)
    b_el = C.ALL4[mask]
    basis = []
    for b in b_el:
        if not b.any():
            continue
        cand = basis + [list(map(int, b))]
        if len(np.unique(C.span(cand), axis=0)) == 8 ** len(cand):
            basis = cand
        if len(basis) == 2:
            break
    words = C.gen_words(L1, basis)
    rows = np.array(C.to_rows(words, None), dtype=np.int64)       # (16, 32)
    # direct enumeration of all 2^16 words
    coeffs = ((np.arange(65536)[:, None] >> np.arange(16)[None, :]) & 1).astype(np.int64)
    allw = (coeffs @ rows) % 2                                     # (65536, 32)
    x = allw[:, :16].reshape(-1, 8, 2)
    z = allw[:, 16:].reshape(-1, 8, 2)
    party_nz = ((x | z).sum(axis=2) > 0).sum(axis=1)
    h_direct = np.bincount(party_nz, minlength=9)
    # predicted
    sym = C.words_of(a_el, b_el)                                    # (4096, 7)
    h_pred = np.zeros(9, dtype=np.int64)
    h_pred += np.bincount((sym != 0).sum(axis=1), minlength=9)      # c = 0
    for c in range(1, 16):
        h_pred += np.bincount((sym != c).sum(axis=1) + 1, minlength=9)
    print("trial", trial, "rank", C.rank2(rows.tolist()), "iso", C.isotropic(rows.tolist()), "direct", h_direct.tolist(), "pred", h_pred.tolist(), "match", bool((h_direct == h_pred).all()))
