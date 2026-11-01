#!/usr/bin/env python3
"""Build negative controls for the coupling checker (deterministic, stdlib). The checker must FAIL on each.

  neg_regrouped_a.json  AME(15,44) with one F_2 coordinate of party 0 exchanged with party 1 in ingredient a
                        (same stabilizer group, regrouped: isotropy and rank survive, 7-uniformity does not)
  neg_a_gt_b.json       AME(11,22) with the roles of a and b exchanged (violates a <= b and the uniformity orders)
  neg_a_gt_b_only.json  a = the AME(12,11) ingredient of AME(11,22) (6-uniform, hence 5-uniform), b = the AME(12,7)
                        ingredient of AME(11,14): both ingredients pass, only a <= b fails (a = 11 > b = 7)
  neg_b_not_ame.json    AME(11,22) with ingredient b replaced by the 12-party GHZ state over F_11 (X^(x12) and
                        Z_i Z_(i+1)^-1): isotropic, rank 12, only 1-uniform, so not AME(12,11)
"""
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DIG = "0123456789abcdefghijklmnopqrstuvwxyz"


def rows_of(obj):
    return [[DIG.index(c) for c in r] if isinstance(r, str) else list(r) for r in obj["generators"]]


def swap_cols(obj, i, j):
    """Exchange qudit i and qudit j (both X and Z columns) in every row."""
    rows = rows_of(obj)
    N = len(rows)
    for r in rows:
        r[i], r[j] = r[j], r[i]
        r[N + i], r[N + j] = r[N + j], r[N + i]
    out = dict(obj)
    out["generators"] = rows
    return out


def ghz(n, p):
    """GHZ stabilizer state on n parties of prime dimension p: X^(xn) and Z_i Z_(i+1)^(-1)."""
    rows = [[1] * n + [0] * n]
    for i in range(n - 1):
        z = [0] * n
        z[i], z[i + 1] = 1, p - 1
        rows.append([0] * n + z)
    return {"type": "stabilizer", "n": n, "d": p, "p": p, "generators": rows, "note": f"GHZ state, {n} parties, F_{p}"}


c = json.loads((HERE / "certificates" / "AME_15_44.json").read_text())
c["a"] = swap_cols(c["a"], 1, 2)  # party 0 owns qudits 0,1; party 1 owns 2,3: exchange qudit 1 with qudit 2
(HERE / "negative" / "neg_regrouped_a.json").write_text(json.dumps(c) + "\n")
c = json.loads((HERE / "certificates" / "AME_11_22.json").read_text())
c["a"], c["b"] = c["b"], c["a"]
(HERE / "negative" / "neg_a_gt_b.json").write_text(json.dumps(c) + "\n")
c = json.loads((HERE / "certificates" / "AME_11_22.json").read_text())
c["a"] = c["b"]
c["b"] = json.loads((HERE / "certificates" / "AME_11_14.json").read_text())["b"]
c["note"] = "negative control: a = AME(12,11), b = AME(12,7); only a <= b fails"
(HERE / "negative" / "neg_a_gt_b_only.json").write_text(json.dumps(c) + "\n")
c = json.loads((HERE / "certificates" / "AME_11_22.json").read_text())
c["b"] = ghz(12, 11)
(HERE / "negative" / "neg_b_not_ame.json").write_text(json.dumps(c) + "\n")
print("wrote 4 negative controls")
