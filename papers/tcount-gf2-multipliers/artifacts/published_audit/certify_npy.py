#!/usr/bin/env python3
"""Map a released parity matrix (rows = phase gadgets, columns = qubits) to the polynomial basis with
qlib/basis_map.py and run the frozen verifier on it.

usage: python3 certify_npy.py FILE.npy [OUT.json]      (needs numpy)

The script works in two layouts:
  * in the release bundle: published_audit/certify_npy.py, with the mapping code in
    published_audit/qlib/ and the frozen verifier in verifier/ (which reads data/records.json);
  * in the project repository: tools/published_audit/certify_npy.py, with qlib/ and verifier/ at the
    problem root.
The verifier compiles verifier/phasecheck next to its source on first use (n <= 8 uses it for the
circuit check).  The release checker ignores that build product.  The line "pins" of the verifier
output reports MODIFIED: phasecheck on other machines, because the pinned hash of the compiled
binary is machine-specific; the pins of verify.py, phasecheck.c and records.json must match.
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
for QLIB, ROOT in ((HERE / "qlib", HERE.parent), (HERE.parents[1] / "qlib", HERE.parents[1])):
    if (QLIB / "basis_map.py").exists() and (ROOT / "verifier" / "verify.py").exists():
        break
else:
    sys.exit("certify_npy.py: qlib/basis_map.py or verifier/verify.py not found")
sys.path.insert(0, str(QLIB))
sys.path.insert(0, str(ROOT / "verifier"))
import basis_map  # noqa: E402
import verify as V  # noqa: E402

M = np.load(sys.argv[1], allow_pickle=False).astype(np.uint8)
n = M.shape[1] // 3
F = [sum(int(b) << i for i, b in enumerate(r)) for r in M]
P = None
for seed in range(1, 9):
    try:
        P, mapped = basis_map.standardize(F, n, seed=seed)
        break
    except Exception:
        P = None
name = pathlib.Path(sys.argv[1]).name
if P is None:
    print(name, "T=%d: no isotopy to GF(2^%d) found" % (len(F), n))
    sys.exit(1)
cand = {"n": n, "modulus": format(P, "b"), "factors": [format(u, f"0{3 * n}b")[::-1] for u in mapped],
        "author": "published", "construction": "released matrix " + name + " mapped by qlib/basis_map.py"}
rep = V.verify_candidate(cand, quick=n > 8)
print(name, "->", rep["summary"])
pins = rep.get("pins", "")
src_ok = pins == "ok" or pins == "MODIFIED: phasecheck"
print(name, "   tensor check:", rep.get("tensor_check"), "| pins:", pins,
      "(sources and record table match)" if src_ok else "(SOURCE OR RECORD TABLE DIFFERS FROM THE PINS)")
if len(sys.argv) > 2:
    pathlib.Path(sys.argv[2]).write_text(json.dumps(cand) + "\n")
sys.exit(0 if rep.get("tensor_check") and src_ok else 1)
