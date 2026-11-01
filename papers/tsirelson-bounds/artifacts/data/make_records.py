#!/usr/bin/env python3
"""Build data/catalog.json, data/records.json and data/targets.json from the frozen TeX snapshots
in data/sources/ (Pal-Vertesi 2009 tables; Brunner-Gisin 2008 appendix).  Stdlib only.

  python3 data/make_records.py            # rebuild the three JSON files
  python3 data/make_records.py --check    # rebuild in memory and compare with the files on disk

Conventions (Pal-Vertesi, Eq. (2)): Alice has projectors A_1..A_mA, Bob B_1..B_mB (outcome "1").
A Bell expression is  sum_x bA[x] p(A_x) + sum_y bB[y] p(B_y) + sum_xy bAB[x][y] p(A_x B_y) <= b0,
with p(A_x) = <A_x (x) 1>, p(A_x B_y) = <A_x (x) B_y>.  The "violation" of a strategy is its value
minus b0 (b0 = 0 for every inequality used as a target).
"""
import json
import pathlib
import re
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "sources"
DELTA = "0.000001"  # minimum improvement over a published value (and alarm tolerance)


def parse_pv_tables():
    tex = (SRC / "pv2009_tables.tex").read_text()
    tabs = dict((re.search(r"\\label\{([^}]*)\}", t).group(1), t)
                for t in re.findall(r"\\begin\{table\*\}.*?\\end\{table\*\}", tex, re.S))
    # Table I: J^k_4422 coefficient strings (single-digit integers, signs inline)
    J = {}
    for k, code in re.findall(r"\$J\^\{(\d+)\}_\{4422\}\$&([-0-9]+)", tabs["table:newineqs"]):
        toks = [int(t) for t in re.findall(r"-?\d", code)]
        if len(toks) != 24:
            raise ValueError(f"J{k}: {len(toks)} coefficients")
        J[int(k)] = {"bA": toks[0:4], "bB": toks[4:8],
                     "bAB": [toks[8 + 4 * i: 12 + 4 * i] for i in range(4)]}
    if sorted(J) != list(range(1, 130)):
        raise ValueError("Table I incomplete")
    # Tables II and III: name & violation & H & D_A & D_B & level
    vals = {}
    cell = re.compile(r"\$([^$]+)\$&([0-9.]+)&\$\\(RR|CC)\^(\d)\$&(\d+)&(\d+)&(\w+)")
    for key in ("table:b5", "table:vt"):
        for name, v, fld, dim, da, db, lev in cell.findall(tabs[key]):
            vals[pv_name(name)] = {"violation": v, "field": "real" if fld == "RR" else "complex",
                                   "dim": int(dim), "proj_ranks_A": da, "proj_ranks_B": db, "level": lev}
    # Table IV: gaps (upper - lower) at increasing NPA levels for the 20 open cases
    gaps = {}
    for line in tabs["table:highlowdiff"].splitlines():
        m = re.match(r"\$([^$]+)\$((?:&[0-9.]*)+)\\\\", line.strip())
        if m:
            g = [x for x in m.group(2).split("&")[1:] if x]
            gaps[pv_name(m.group(1))] = g
    return J, vals, gaps


def parse_pv2010():
    """Pal-Vertesi 2010 (PRA 82, 022116; arXiv:1006.3032) Table I: updated windows for the 19 open cases.
    Row: name & lower & H & new & (upper - lower) & level."""
    tex = (SRC / "pv2010_table1.tex").read_text()
    out = {}
    row = re.compile(r"\$([^$]+)\$&([0-9.]+)&\$\\(RR|CC)\^\{(\d+)\}\$&(yes|no)&([0-9.]+)&\$([^$]*)\$")
    for name, lo, fld, dim, new, gap, lev in row.findall(tex):
        out[pv_name(name)] = {"lower": lo, "gap": gap, "field": "real" if fld == "RR" else "complex",
                              "dim": int(dim), "new_in_2010": new == "yes", "level": lev.replace("'", "")}
    if len(out) != 19:
        raise ValueError(f"PV2010 table: {len(out)} rows")
    return out


def parse_alist():
    """Ito-Imai-Avis A-list from the Naceur-Wang-Magron-Acin repository (MIT; data/sources/naceur2025_alist/):
    each A<k>.txt has m=, n= and a polynomial f in projectors X1..Xm (Alice), Y1..Yn (Bob) to be MINIMIZED;
    the Pal-Vertesi form (classical value 0) is -f."""
    out = {}
    d = SRC / "naceur2025_alist"
    term = re.compile(r"([+-]?\d+)\*((?:[XY]\d+)(?:\*[XY]\d+)?)")
    for f in sorted(d.glob("A*.txt"), key=lambda q: int(q.stem[1:])):
        lines = f.read_text().split("\n")
        mA, mB = int(lines[0].split("=")[1]), int(lines[1].split("=")[1])
        poly = lines[2].replace(" ", "")
        if term.sub("", poly):
            raise ValueError(f"{f.name}: unparsed objective {term.sub('', poly)!r}")
        bA, bB, bAB = [0] * mA, [0] * mB, [[0] * mB for _ in range(mA)]
        for c, mono in term.findall(poly):
            c = -int(c)  # PV form = -f
            ops = mono.split("*")
            if len(ops) == 1:
                k = int(ops[0][1:]) - 1
                if ops[0][0] == "X":
                    bA[k] += c
                else:
                    bB[k] += c
            else:
                x, y = sorted(ops)  # X before Y
                if x[0] != "X" or y[0] != "Y":
                    raise ValueError(f"{f.name}: unexpected monomial {mono}")
                bAB[int(x[1:]) - 1][int(y[1:]) - 1] += c
        if local_max(bA, bB, bAB) != 0:  # A10 is stored in PV form already: take the sign with local bound 0
            bA, bB, bAB = [-x for x in bA], [-x for x in bB], [[-x for x in r] for r in bAB]
            if local_max(bA, bB, bAB) != 0:
                raise ValueError(f"{f.name}: neither sign has local bound 0")
        out[f.stem] = {"bA": bA, "bB": bB, "bAB": bAB, "b0": 0}
    return out


def local_max(bA, bB, bAB):
    best = None
    for a in range(1 << len(bA)):
        av = [(a >> x) & 1 for x in range(len(bA))]
        v = sum(bA[x] * av[x] for x in range(len(bA)))
        v += sum(max(0, bB[y] + sum(bAB[x][y] * av[x] for x in range(len(bA)))) for y in range(len(bB)))
        best = v if best is None or v > best else best
    return best


def parse_alist_bounds():
    """Order-3 and order-4 NPA values for the open A cases from maximal_violations_exact.txt (numerical,
    published in the repository, not in the paper): name, -exact, -bound3, 3, gap3, t3, -bound4, 4, gap4, t4."""
    tex = (SRC / "naceur2025_alist" / "maximal_violations_exact.txt").read_text()
    out = {}
    pat = re.compile(r"(A\d+) (-?[\d.]+) (-?[\d.]+) 3 ([\d.]+) ([\d.]+) (-?[\d.]+) 4 ([\d.]+) ([\d.]+)")
    for name, ex, b3, g3, t3, b4, g4, t4 in pat.findall(tex):
        out[name] = {"exact": ex.lstrip("-"), "order3": b3.lstrip("-"), "order4": b4.lstrip("-"), "gap4": g4}
    return out


def pv_name(s):
    s = s.replace(" ", "")
    m = re.fullmatch(r"J\^\{(\d+)\}_\{4422\}", s)
    if m:
        return f"J{m.group(1)}"
    m = re.fullmatch(r"I\^\{?(\d+)\}?_\{(\d{4})\}", s)
    if m:
        return f"I{m.group(2)}_{m.group(1)}"
    m = re.fullmatch(r"I_\{(\d{4})\}", s)
    if m:
        return f"I{m.group(1)}"
    m = re.fullmatch(r"(A|AS|AII)_\{?(\d+)\}?", s)
    if m:
        return f"{m.group(1)}{m.group(2)}"
    raise ValueError(f"unknown name {s!r}")


def parse_bg_appendix():
    tex = (SRC / "bg2008_appendix.tex").read_text()
    out = {}
    pat = re.compile(r"((?:I|A|AS|AII)_\{[^}]*\}(?:\^\{?\d+\}?)?)\s*\\equiv\s*\\begin\{array\}\{[^}]*\}"
                     r"(.*?)\\end\{array\}\s*\\leq\s*(-?\d+)", re.S)
    for label, body, b0 in pat.findall(tex):
        rows = [r.strip() for r in body.replace("\\hline", "").split("\\\\") if r.strip()]
        head = [int(x) for x in rows[0].split("&")[1:]]
        bA, bAB = [], []
        for r in rows[1:]:
            xs = [int(x) for x in r.split("&")]
            bA.append(xs[0])
            bAB.append(xs[1:])
        if any(len(r) != len(head) for r in bAB):
            raise ValueError(f"{label}: ragged table")
        name = bg_name(label)
        out[name] = {"bA": bA, "bB": head, "bAB": bAB, "b0": int(b0)}
    return out


def bg_name(label):
    s = label.replace(" ", "")
    m = re.fullmatch(r"I_\{4422\}\^\{?(\d+)\}?", s)
    if m:
        return f"I4422_{m.group(1)}"
    m = re.fullmatch(r"(A|AS|AII)_\{(\d+)\}", s)
    if m:
        return f"{m.group(1)}{m.group(2)}"
    raise ValueError(label)


def build():
    J, vals, gaps = parse_pv_tables()
    bg = parse_bg_appendix()
    cat = {
        "CHSH": {"bA": [-1, 0], "bB": [-1, 0], "bAB": [[1, 1], [1, -1]], "b0": 0,
                 "source": "Pal-Vertesi 2009, Eq. (3)"},
        "I3322": {"bA": [-2, -1, 0], "bB": [-1, 0, 0], "bAB": [[1, 1, 1], [1, 1, -1], [1, -1, 0]], "b0": 0,
                  "source": "Pal-Vertesi 2009, Eq. (5)"},
    }
    for name, c in bg.items():
        cat[name] = dict(c, source="Brunner-Gisin 2008, Appendix A")
    for k, c in J.items():
        cat[f"J{k}"] = dict(c, b0=0, source="Pal-Vertesi 2009, Table I")
    for name, c in parse_alist().items():
        if name in cat or name in ("A2", "A3", "A4", "A7"):
            continue  # CHSH, I3322, I3422_2, I4422_1 under other names; A5, A6 from Brunner-Gisin
        cat[name] = dict(c, source="Ito-Imai-Avis A-list via Naceur et al., github.com/nininaceur/CertifiedQuantumBounds "
                                   "(MIT), data/sources/naceur2025_alist")
    for name, c in cat.items():
        c["mA"], c["mB"] = len(c["bA"]), len(c["bB"])

    recs = {}
    for name in cat:
        if name == "CHSH":
            recs[name] = {"lower": "0.2071068", "upper": "0.2071068", "status": "closed",
                          "note": "Tsirelson: (sqrt2-1)/2", "source": "Tsirelson 1980"}
            continue
        v = vals.get(name)
        if v is None:
            continue
        r = {"lower": v["violation"], "pv_field": v["field"], "pv_dim": v["dim"],
             "pv_proj_ranks": [v["proj_ranks_A"], v["proj_ranks_B"]], "pv_level": v["level"],
             "source": "Pal-Vertesi 2009, Tables II-IV"}
        if v["level"] == "no":
            g = gaps[name]
            r["gaps_by_level"] = g
            r["upper"] = dec(Fraction(v["violation"]) + Fraction(g[-1]))
            r["status"] = "open"
        else:
            r["upper"] = v["violation"]
            r["status"] = "matched"  # PV's NPA bound met their strategy (numerically)
        recs[name] = r
    pv10 = parse_pv2010()
    for name, u in pv10.items():
        if name not in recs:
            continue  # A-list rows: no coefficients yet
        r = recs[name]
        r["pv2009"] = {"lower": r["lower"], "upper": r["upper"]}
        r.update(lower=u["lower"], upper=dec(Fraction(u["lower"]) + Fraction(u["gap"])), pv_field=u["field"],
                 pv_dim=u["dim"], pv_level=u["level"], source="Pal-Vertesi 2010, Table I (supersedes 2009)")
        r["status"] = "open" if Fraction(u["gap"]) > 0 else "closed_numerically"
        if name == "J97":
            r["note"] = ("PV2010: an L3+aabb run (Borchers' code, reduced accuracy) gave a value below the lower "
                         "bound; published reliable upper bound is L3")
    # Naceur et al. repository: numerical NPA order-4 values for the open A cases (lower than PV 2010's)
    for name, nb in parse_alist_bounds().items():
        r = recs.get(name)
        if r is None:
            continue
        if Fraction(nb["order4"]) < Fraction(r["upper"]):
            r.setdefault("pv2010", {"lower": r["lower"], "upper": r["upper"]})
            r["upper"] = nb["order4"]
            r["naceur2025"] = {"order3": nb["order3"], "order4": nb["order4"],
                               "note": "numerical NPA order 4 from the repository's maximal_violations_exact.txt (not in the paper)"}
            r["source"] = r["source"] + "; upper: Naceur et al. 2025 repository, NPA order 4"
            r["pv_level"] = "4 (Naceur et al.)"
            if Fraction(r["upper"]) - Fraction(r["lower"]) <= 0:
                r["status"] = "closed_numerically"
    # phase-1 overlay: our own certified bests (frozen in data/sources/p05_phase1.json) become the bar
    ph1 = json.loads((SRC / "p05_phase1.json").read_text())["best"]
    for key, b in ph1.items():
        name, side = key.split(":")
        r = recs[name]
        better = (Fraction(b["value"]) > Fraction(r["lower"])) if side == "lower" else (Fraction(b["value"]) < Fraction(r["upper"]))
        if name == "J97" and side == "lower":
            better = True  # exact strategy value; PV's printed 1.1583626 is rounded up
        if better:
            r.setdefault("pv2010", {"lower": r["lower"], "upper": r["upper"]})
            r[side] = b["value"]
            r.setdefault("phase1", {})[side] = {"value": b["value"], "file": b["file"], "detail": b["detail"]}
            r["source"] = "Pal-Vertesi 2010, Table I, improved by p05 phase 1 (data/sources/p05_phase1.json)"
    for name, r in recs.items():
        if r["status"] == "open" and Fraction(r["upper"]) - Fraction(r["lower"]) < Fraction(DELTA):
            r["status"] = "closed"
            r["note"] = (r.get("note", "") + " Closed by p05 phase 1: certified window narrower than delta.").strip()
    if recs["I3322"]["status"] == "open":  # superseded: Pal-Vertesi 2010 (arXiv:1006.3032) and the 2026 proofs
        recs["I3322"].update(status="superseded", lower="0.2508754", upper="0.2508754",
                             note="NPA Q4=Q5 match the infinite-dimensional value 0.25087538451 to 1e-16 "
                                  "(Pal-Vertesi 2010); excluded from the targets")
    open_ids = sorted((n for n, r in recs.items() if r["status"] == "open"), key=sort_key)
    targets = []
    for n in open_ids:
        r = recs[n]
        gap = Fraction(r["upper"]) - Fraction(r["lower"])
        for side in ("lower", "upper"):
            targets.append({
                "id": f"{n}:{side}", "band": n, "inequality": n, "side": side,
                "record": r[side], "score": float(gap),
                "display": (f"{'raise' if side == 'lower' else 'lower'} the {side} bound; current window "
                            f"[{r['lower']}, {r['upper']}] gap {dec(gap)}; PV strategy {r['pv_field']}^{r['pv_dim']}, upper at {r['pv_level']}"),
            })
    meta = {"source": "Pal & Vertesi, PRA 79, 022120 (2009), arXiv:0810.1615v4, updated by Pal & Vertesi, PRA 82, 022116 "
                      "(2010), arXiv:1006.3032 Table I; coefficients from PV2009 Table I and Brunner & Gisin, PLA 372, 3162 (2008)",
            "snapshot": "2026-10-11 (2010 overlay 23:55; phase-1 overlay 2026-10-12 23:00)", "delta": DELTA,
            "note": "A-list cases (A14, A21, A47, A62, A64, A68, A80, A82, A84, A89 open in PV2010; A65, A67 closed there) "
                    "need the Ito-Imai-Avis coefficients (web list offline at setup); not targets yet. I4422_20 and J34 "
                    "are closed numerically in PV2010 (gap 0 at L3 + 1/4 aabb)."}
    return ({"meta": meta, "inequalities": cat},
            {"meta": meta, "records": recs},
            {"meta": meta, "targets": targets})


def dec(fr):
    return f"{float(fr):.7f}"


def sort_key(n):
    m = re.match(r"([A-Za-z]+)(\d*)_?(\d*)", n)
    return (m.group(1), int(m.group(2) or 0), int(m.group(3) or 0))


def main():
    files = dict(zip(("catalog.json", "records.json", "targets.json"), build()))
    if "--check" in sys.argv:
        bad = [f for f, d in files.items() if json.loads((HERE / f).read_text()) != json.loads(json.dumps(d))]
        print("ok" if not bad else f"MISMATCH: {bad}")
        sys.exit(1 if bad else 0)
    for f, d in files.items():
        (HERE / f).write_text(json.dumps(d, indent=1) + "\n")
    c = files["catalog.json"]["inequalities"]
    r = files["records.json"]["records"]
    print(f"catalog: {len(c)} inequalities; records: {len(r)} "
          f"({sum(x['status'] == 'open' for x in r.values())} open); targets: {len(files['targets.json']['targets'])}")


if __name__ == "__main__":
    main()
