"""Negative controls: perturb one hint entry in block 0 and rewrite the certificate canonically."""
import json
import sys

src, dst, which = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(src))
entry = d["blocks"][0]["output"][which]["entries"][0]
entry[0] = str(int(entry[0]) + 2)  # stays odd, so still canonical dyadic
if which == "center_dyadic":
    pass  # entry 0 is a diagonal entry (real), so Hermiticity is preserved
open(dst, "wb").write(json.dumps(d, allow_nan=False, ensure_ascii=True, separators=(",", ":"),
                                 sort_keys=False).encode("ascii") + b"\n")
print("perturbed", which, "entry (0,0) of block 0 output ->", entry)
