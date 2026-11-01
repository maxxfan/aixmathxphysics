#!/usr/bin/env python3
"""Check the ingredient-to-cell mapping (stdlib): for every cell in cells.json, the ingredients embedded in its
certificate span exactly the same F_p row spaces as the named files in ingredients/, the certificate's k, n and
dimensions agree with the cell, and the cell is AME(2k+1, a*b). Prints one line per cell."""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
DIG = "0123456789abcdefghijklmnopqrstuvwxyz"


def rref(obj):
    p = int(obj["p"])
    M = [[DIG.index(c) for c in r] if isinstance(r, str) else [int(v) for v in r] for r in obj["generators"]]
    rk = 0
    for col in range(len(M[0])):
        pv = next((i for i in range(rk, len(M)) if M[i][col] % p), None)
        if pv is None:
            continue
        M[rk], M[pv] = M[pv], M[rk]
        iv = pow(M[rk][col], p - 2, p)
        M[rk] = [(v * iv) % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(u - f * v) % p for u, v in zip(M[i], M[rk])]
        rk += 1
    return (p, int(obj["n"]), int(obj["d"]), tuple(tuple(r) for r in M[:rk]))


cells = json.loads((HERE / "cells.json").read_text())["cells"]
ok = True
for c in cells:
    cert = json.loads((HERE / c["certificate"]).read_text())
    ia = json.loads((HERE / "ingredients" / f"{c['a']}.json").read_text())
    ib = json.loads((HERE / "ingredients" / f"{c['b']}.json").read_text())
    same_a, same_b = rref(cert["a"]) == rref(ia), rref(cert["b"]) == rref(ib)
    k = int(cert["k"])
    shape = (int(cert["a"]["n"]) == int(cert["b"]["n"]) == 2 * k + 2 and c["n"] == 2 * k + 1
             and c["d"] == int(cert["a"]["d"]) * int(cert["b"]["d"]) and int(cert["a"]["d"]) <= int(cert["b"]["d"]))
    good = same_a and same_b and shape
    ok &= good
    print(f"{c['cell']:11s} k={k} a={c['a']:34s} (d={cert['a']['d']}) b={c['b']:26s} (d={cert['b']['d']})  "
          f"a-match={'yes' if same_a else 'NO'} b-match={'yes' if same_b else 'NO'} shape={'ok' if shape else 'BAD'}")
print("MAPPING", "PASS" if ok else "FAIL", f"({len(cells)} cells)")
sys.exit(0 if ok else 1)
