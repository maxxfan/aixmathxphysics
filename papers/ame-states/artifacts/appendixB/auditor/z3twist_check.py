"""Sanity check of z3twist2 output: random full completions T are Omega-invariant, isotropic, of dimension 10 (case nonfpf) / 12 (fpf),
inside V_nf, with direct min party weight >= 5.  Usage: z3twist_check.py CASE"""
import sys
from z3twist2 import *
import random
random.seed(11)
allT = [(key, kk, b) for key, st in results.items() for kk, b in st.items() if len(b) == full_dim]
print("number of T:", len(allT))
cnt = 0
for key, kk, b in random.sample(allT, min(12, len(allT))):
    basis = list(b.values())
    ws = words_of(basis)
    assert len(set(int(w) for w in ws)) == len(ws) == 2 ** full_dim
    minw = pw(ws[1:]).min()
    # isotropic
    iso = all(int(omg(np.uint32(u), np.uint32(v))) == 0 for u in basis for v in basis)
    # invariant
    wset = set(int(w) for w in ws)
    inv = all(int(OmN(np.array([w], dtype=np.uint32))[0]) in wset for w in basis)
    # inside V_nf: 1 + Omega + Omega^2 = 0 on basis
    nf = all((int(w) ^ int(OmN(np.array([w], dtype=np.uint32))[0]) ^ int(OmN(OmN(np.array([w], dtype=np.uint32)))[0])) == 0 for w in basis)
    print("T: dim", len(basis), "min weight", minw, "isotropic", iso, "Omega-invariant", inv, "in V_nf", nf)
    assert minw >= 5 and iso and inv and nf
    cnt += 1
print("checked", cnt, "random T: all fine")
