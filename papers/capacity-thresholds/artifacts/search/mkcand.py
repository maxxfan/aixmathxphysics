#!/usr/bin/env python3
"""Turn an optimised state into a candidate at a safe rational p just below its float threshold.

usage: mkcand.py STATE.json [--target FAMILY] [--margin BITS]
p is the largest multiple of 1e-8 with I_c(p) >= margin bits (float screen); the candidate goes
to cands/<target>_n<n>_<p>.json for verifier/verify.py.
"""
import argparse
import json
import pathlib
import sys
from decimal import Decimal

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pauliopt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("--target")
    ap.add_argument("--margin", type=float, default=2e-9)
    a = ap.parse_args()
    d = json.loads(pathlib.Path(a.state).read_text())
    fam = a.target or d.get("target", "depolarizing")
    V = pauliopt.load_rows(a.state)
    th = pauliopt.threshold(fam, V, margin=a.margin)  # I_c(th) = margin; I_c nonincreasing in p
    q = int(th * 1e8)
    while pauliopt.ic_bits(fam, V, q / 1e8) < a.margin:
        q -= 1
    p = str(Decimal(q) / Decimal(10 ** 8))
    n = V.shape[1] - 1
    cand = {"target": fam, "p": p, "rows": [[[z.real, z.imag] for z in row] for row in V],
            "note": f"n={n} rank-{V.shape[0]} symmetric state, float threshold {th:.10f}; "
                    f"L-BFGS with the analytic gradient (work/lead/depopt.py, pauliopt.py) from {pathlib.Path(a.state).name}"}
    out = HERE / "cands" / f"{fam}_n{n}_{p}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(cand))
    print(f"{out}  p={p}  I_c(float)={pauliopt.ic_bits(fam, V, q / 1e8):.3e} bits  threshold {th:.10f}")


if __name__ == "__main__":
    main()
