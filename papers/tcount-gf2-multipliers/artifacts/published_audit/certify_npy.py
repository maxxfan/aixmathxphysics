#!/usr/bin/env python3
"""Map a released parity matrix (rows = phase gadgets, columns = qubits) to the polynomial basis with
qlib/basis_map.py and run the frozen verifier on it.  usage: certify_npy.py FILE.npy [OUT.json]"""
import json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "qlib")); sys.path.insert(0, str(ROOT / "verifier"))
import basis_map, verify as V
M = np.load(sys.argv[1], allow_pickle=False).astype(np.uint8)
n = M.shape[1] // 3
F = [sum(int(b) << i for i, b in enumerate(r)) for r in M]
for seed in range(1, 9):
    try:
        P, mapped = basis_map.standardize(F, n, seed=seed)
        break
    except Exception as e:
        P = None
if P is None:
    print(pathlib.Path(sys.argv[1]).name, "T=%d: no isotopy to GF(2^%d) found" % (len(F), n)); sys.exit(1)
cand = {"n": n, "modulus": format(P, "b"), "factors": [format(u, f"0{3*n}b")[::-1] for u in mapped],
        "author": "published", "construction": "released matrix " + pathlib.Path(sys.argv[1]).name + " mapped by qlib/basis_map.py"}
rep = V.verify_candidate(cand, quick=n > 8)
print(pathlib.Path(sys.argv[1]).name, "->", rep["summary"])
if len(sys.argv) > 2:
    pathlib.Path(sys.argv[2]).write_text(json.dumps(cand) + "\n")
