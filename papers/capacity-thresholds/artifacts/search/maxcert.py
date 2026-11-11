#!/usr/bin/env python3
"""Highest p at which verifier/verify.py certifies a saved state.

The large-n states are flat in p (|dI_c/dp| ~ 1e-4 bits per unit p at n >= 64), so a fixed float
margin can waste up to ~7e-5 in p. This bisects p between a certified point (the threshold of
I_c - MARGIN) and the float zero crossing, running the frozen verifier at every step.

usage: maxcert.py STATE.json [--target FAMILY] [--margin BITS] [--tol P] [--min-lower BITS]

--min-lower L asks for the highest p whose certified lower bound on I_c is at least L bits (a
conservative point for publication) instead of merely positive.
"""
import argparse
import json
import pathlib
import subprocess
import sys
from decimal import Decimal

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pauliopt  # noqa: E402
import yframe  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
VERIFY = HERE.parents[1] / "verifier" / "verify.py"


def verify(fam, V, q, note, min_lower=0.0):
    p = str(Decimal(q) / Decimal(10 ** 8))
    out = HERE / "cands" / f"{fam}_n{V.shape[1] - 1}_{p}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"target": fam, "p": p, "note": note,
                               "rows": [[[z.real, z.imag] for z in row] for row in V]}))
    res = subprocess.run([sys.executable, str(VERIFY), "--json", str(out)], capture_output=True, text=True)
    rep = json.loads(res.stdout)
    print(f"  {rep['summary']}  ({rep['seconds']} s)", flush=True)
    return rep["verdict"] == "IMPROVEMENT" and float(rep["I_c_lo"]) >= min_lower, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("--target")
    ap.add_argument("--margin", type=float, default=4e-9)
    ap.add_argument("--tol", type=float, default=5e-7)
    ap.add_argument("--min-lower", type=float, default=0.0)
    a = ap.parse_args()
    fam = a.target or json.loads(pathlib.Path(a.state).read_text()).get("target", "depolarizing")
    V = pauliopt.load_rows(a.state)
    V /= np.linalg.norm(V)
    F, W = (fam, V) if fam == "depolarizing" else (fam + "@y", yframe.to_frame(V))
    p0 = pauliopt.threshold(F, W)
    slope = (pauliopt.ic_bits(F, W, p0) - pauliopt.ic_bits(F, W, p0 - 1e-5)) / 1e-5
    pm = pauliopt.threshold(F, W, margin=a.margin)
    print(f"{fam} n={V.shape[1] - 1}: float zero crossing {p0:.10f}, dI_c/dp {slope:.3e} bits, "
          f"margin-{a.margin:g} threshold {pm:.10f}", flush=True)
    note = (f"n={V.shape[1] - 1} rank-{V.shape[0]} symmetric state from {pathlib.Path(a.state).name}; "
            f"float zero crossing {p0:.10f}; p chosen by bisection against verifier/verify.py (work/lead/maxcert.py)")
    lo, hi = int(pm * 1e8), int(p0 * 1e8)
    ok, best = verify(fam, V, lo, note, a.min_lower)
    while not ok:
        lo -= max(hi - lo, 100)
        ok, best = verify(fam, V, lo, note, a.min_lower)
    while hi - lo > a.tol * 1e8:
        mid = (lo + hi) // 2
        ok, f = verify(fam, V, mid, note, a.min_lower)
        if ok:
            lo, best = mid, f
        else:
            hi = mid
    print(f"best certified p = {lo / 1e8:.8f}: {best}", flush=True)


if __name__ == "__main__":
    main()
