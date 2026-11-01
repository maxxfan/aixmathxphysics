#!/usr/bin/env python3
"""Reproducible check of the A19 case study (standard library only; run from anywhere).

1. The inequality: coefficients from manifest.json (Collins-Gisin form), compared with the source file
   data/naceur2025_A19.txt (our W is minus their polynomial), local bound 0 by enumerating all 2^4 * 2^5
   deterministic strategies.
2. Four certificates, each checked exactly by ../../replay/exactcheck.py:
     complex qubit strategy            (the published optimum, rounded: 0.6226300)
     real dimension-4 strategy          (rank-2 projectors)
     SOS certificate, N = 202 words     (NPA level 2 + aab + abb)
     SOS certificate, N = 318 words     (NPA level 3)
   For the SOS certificates, positive semidefiniteness of Y + tau I is proven at the stated shift by the
   exact-integer residual (and the bound lambda + N tau is recomputed exactly).
3. The comparison with the published value 0.6226300 (Pal-Vertesi 2009, Table II; same value in the
   benchmark data of Naceur et al.).
"""
import itertools
import json
import math
import pathlib
import sys
from fractions import Fraction as Q

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "replay"))
import exactcheck as X  # noqa: E402

PUBLISHED = Q("0.6226300")


def naceur_polynomial(path):
    """{(x, y): coefficient} of the polynomial in data/naceur2025_A19.txt (x, y 1-based; 0 = no letter)."""
    line = path.read_text().splitlines()[2].replace(" ", "")
    poly = {}
    for term in line.replace("-", "+-").split("+"):
        if not term:
            continue
        c, *letters = term.split("*")
        x = y = 0
        for l in letters:
            if l[0] == "X":
                x = int(l[1:])
            else:
                y = int(l[1:])
        poly[(x, y)] = poly.get((x, y), 0) + int(c)
    return poly


def main():
    manifest = json.loads((ROOT / "manifest.json").read_text())
    e = next(x for x in manifest["inequalities"] if x["name"] == "A19")
    co = e["coefficients"]
    cat = {"A19": co}
    print("A19 in Collins-Gisin form (violation = <W> - b0, b0 = %d):" % co["b0"])
    print("  bA  =", co["bA"])
    print("  bB  =", co["bB"])
    print("  bAB =", co["bAB"])
    # 1a. the source file, with the opposite sign
    src = naceur_polynomial(HERE / "data" / "naceur2025_A19.txt")
    ours = {}
    for x in range(co["mA"]):
        if co["bA"][x]:
            ours[(x + 1, 0)] = -co["bA"][x]
        for y in range(co["mB"]):
            if co["bAB"][x][y]:
                ours[(x + 1, y + 1)] = -co["bAB"][x][y]
    for y in range(co["mB"]):
        if co["bB"][y]:
            ours[(0, y + 1)] = -co["bB"][y]
    print("  equals minus the polynomial of data/naceur2025_A19.txt:", ours == src)
    # 1b. local bound
    best = max(sum(co["bA"][x] * a[x] for x in range(co["mA"])) + sum(co["bB"][y] * b[y] for y in range(co["mB"]))
               + sum(co["bAB"][x][y] * a[x] * b[y] for x in range(co["mA"]) for y in range(co["mB"]))
               for a in itertools.product((0, 1), repeat=co["mA"]) for b in itertools.product((0, 1), repeat=co["mB"]))
    print("  local bound (max over 512 deterministic strategies):", best)
    ok = ours == src and best == co["b0"] == 0

    reports = json.loads((HERE / "data" / "board_reports.json").read_text())
    print("\nexact certificates:")
    vals = {}
    for rel, label in (("case_studies/A19/data/complex_qubit_lower.json", "complex qubit strategy"),
                       ("certificates/A19/lower.json", "real d=4 strategy, rank-2 projectors"),
                       ("case_studies/A19/data/upper_N202.json", "SOS, NPA level 2 + aab + abb"),
                       ("certificates/A19/upper.json", "SOS, NPA level 3")):
        cand = json.loads((ROOT / rel).read_text())
        rep = reports[rel]["report"]
        if cand["type"] == "strategy":
            out = X.check(cand, cat)
            v = Q(out["violation_exact"])
            vals[label] = v
            print(f"  {label:40s} d={cand['dA']} {cand['field']:7s} violation = {float(v):.12f}  (exact, "
                  f"{len(out['violation_exact'])} characters)   board #{reports[rel]['board_id']}")
            ok &= X.rational(rep["violation"]) <= v
        else:
            tau = Q(rep["sigma"]) + Q(rep["residual_bound"])
            out = X.check(cand, cat, rep["sigma"], tau)
            bound = Q(out["lambda_exact"]) + out["n_words"] * tau - co["b0"]
            vals[label] = bound
            print(f"  {label:40s} N={out['n_words']:3d} lambda = {float(Q(out['lambda_exact'])):.12f}  "
                  f"tau_stated = {float(tau):.3e}  tau_exact = {out['tau_exact_float']:.3e}  "
                  f"bound = {float(bound):.12f}   board #{reports[rel]['board_id']}")
            ok &= out["tau_ok"] and out.get("bareiss_at_stated_tau", True) and bound <= X.rational(rep["bound"])
    lo, up = vals["real d=4 strategy, rank-2 projectors"], vals["SOS, NPA level 3"]
    cq = vals["complex qubit strategy"]
    down = math.floor(lo * 10 ** 10) / 10 ** 10
    upr = math.ceil(up * 10 ** 10) / 10 ** 10
    print(f"\ncertified window      [{down:.10f}, {upr:.10f}] (rounded outward)  exact width {float(up - lo):.2e}")
    print(f"published value       {float(PUBLISHED):.7f}  (Pal-Vertesi 2009 Table II: C^2, ranks 1111/11111, NPA level 2a;"
          f" data/pv2009_table_row.tex; also data/naceur2025_maximal_violations_A19.txt)")
    print(f"complex qubit value   {float(cq):.10f}  (rounds to {float(cq):.7f}: the published strategy value)")
    print(f"lower - published     {float(lo - PUBLISHED):.3e}  (> 0: an exact strategy beats the published maximum"
          f" and the published NPA upper bound)")
    print(f"lower - complex qubit {float(lo - cq):.3e}")
    ok &= lo - PUBLISHED > Q(4, 10 ** 7) and up >= lo
    print("\nA19 case study:", "REPRODUCED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
