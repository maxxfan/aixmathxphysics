"""Independent shadow-inequality audit of data/status.json (auditor, shift 38; cross-check of extend #695).

Shadow inequality in purity form (Huber-Eltschka-Siewert-Guehne 2018, any local dimension D):
    W(S) = sum_{T subset [n]} (-1)^{|S cap T|} Tr(rho_T^2)  >= 0   for every S subset [n].
For an AME(n, D) state Tr(rho_T^2) = D^{-min(|T|, n-|T|)}, so W depends on s = |S| only:
    W(s) = sum_{t=0}^{n} sum_{j} (-1)^j C(s,j) C(n-s,t-j) D^{-min(t,n-t)}   (exact Fractions).
Usage: python3 shadow_check.py            -> audit ../../snapshot/status_phase3.json (or $STATUS_JSON): every 'exists'/'verified' cell must pass, every cell whose source
                                             mentions Shadow must fail; list of violations printed, plus new exclusions among open cells.
       python3 shadow_check.py n D [...]  -> print pass/fail and (min W, s) for the given cells.
"""
import json, sys
from fractions import Fraction
from math import comb


_KRAW = {}


def kraw(n):
    """K[s][t] = coefficient of z^t in (1-z)^s (1+z)^(n-s): sum_j (-1)^j C(s,j) C(n-s,t-j) (integers, cached per n)."""
    if n in _KRAW:
        return _KRAW[n]
    K = []
    for s in range(n + 1):
        poly = [1]
        for _ in range(s):  # multiply by (1 - z)
            new = poly + [0]
            for i, c in enumerate(poly):
                new[i + 1] -= c
            poly = new
        for _ in range(n - s):  # multiply by (1 + z)
            new = poly + [0]
            for i, c in enumerate(poly):
                new[i + 1] += c
            poly = new
        K.append(poly)
    _KRAW[n] = K
    return K


def shadow(n, D):
    h = n // 2
    pur = [D ** (h - min(t, n - t)) for t in range(n + 1)]  # purities scaled by D^h (integers)
    K = kraw(n)
    worst = None
    for s in range(n + 1):
        W = sum(K[s][t] * pur[t] for t in range(n + 1))
        if worst is None or W < worst[0]:
            worst = (W, s)
    return worst[0] >= 0, (Fraction(worst[0], D ** h), worst[1])


def main():
    if len(sys.argv) > 2:
        a = sys.argv[1:]
        for i in range(0, len(a), 2):
            n, D = int(a[i]), int(a[i + 1])
            ok, (w, s) = shadow(n, D)
            print(f"AME({n},{D}): shadow {'pass' if ok else 'FAIL'}  min W = {float(w):.4g} at |S| = {s}")
        return
    import os
    path = os.environ.get("STATUS_JSON") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "snapshot", "status_phase3.json")
    st = json.load(open(path))
    cells = st["cells"] if isinstance(st, dict) and "cells" in st else st
    items = list(cells.values()) if isinstance(cells, dict) else cells
    counts = {}
    bad_exists, shadow_src_pass, new_excl, nonex_pass = [], [], [], []
    for c in items:
        n, d, status = c["n"], c["d"], c["status"]
        if n < 2:
            continue
        ok, _ = shadow(n, d)
        counts[(status, ok)] = counts.get((status, ok), 0) + 1
        src = c.get("source", "")
        if status == "exists" and not ok:
            bad_exists.append((n, d, src))
        if "Shadow" in src and ok:
            shadow_src_pass.append((n, d, src))
        if status == "open" and not ok:
            new_excl.append((n, d, src))
        if status == "nonexistent" and ok:
            nonex_pass.append((n, d, src))
    print("counts (status, shadow ok):", {f"{k[0]}/{'pass' if k[1] else 'FAIL'}": v for k, v in sorted(counts.items(), key=str)})
    print("EXISTS cells failing the shadow inequality (would be a contradiction):", bad_exists)
    print("cells whose source says Shadow but my test passes:", shadow_src_pass)
    print("OPEN cells failing the shadow inequality (would be new non-existence):", new_excl)
    print("nonexistent cells that pass shadow (excluded by other arguments, e.g. Scott/Rains/Euler):", nonex_pass[:40], len(nonex_pass))


main()
