#!/usr/bin/env python3
"""Build the self-contained artifact release in release/ from the repository.

Copies the frozen data, the pinned verifier with its tests, every certificate behind Appendix A (byte-identical
to the frozen board submission), the intermediate ladder strategies, the A19 case study inputs and the search
programs; writes manifest.json (exact values, computed with the release's independent checker
replay/exactcheck.py and cross-checked against the verifier's board reports), appendix_A.tex/.tsv, table_I.tsv,
ladders/ladders.csv and SHA256SUMS.  The hand-written parts of the release (README.md, replay/, generate/,
case_studies/A19/README.md and check_A19.py, ladders/README.md) are kept.

usage: python3 results/build_release.py [--jobs 8]
Run after results/build_table.py, results/build_summary.py and paper/make_tables.py."""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import pathlib
import shutil
import sqlite3
import subprocess
import sys
from fractions import Fraction

P = pathlib.Path(__file__).resolve().parents[1]
R = P / "release"
sys.path.insert(0, str(R / "replay"))
import exactcheck as X  # noqa: E402

GENERATED = ["data", "verifier", "certificates", "search", "ladders/certificates", "ladders/notebooks",
             "case_studies/A19/data"]
SEARCH = ["qlib/bell.py", "qlib/certs.py", "work/upper/ipm_half_v.py", "work/upper/price2.py",
          "work/upper2/fastpanels.py", "work/certify/batch.py", "work/certify/ipm_half_v.py",
          "work/certify/price2.py", "work/lower1/bh.py", "work/lower1/insert.py", "work/lower1/pvit.py",
          "work/lower1/pvdbl.py", "work/lower1/mkcand.py", "work/lower1/mkcand3.py"]
VERIFIER = ["verify.py", "qexact.py", "rescheck.c", "PINS.json", "freeze.sh", "tests/test_verifier.py"]


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def dec_down(q):
    return f"{math.floor(q * 10 ** 10) / 10 ** 10:.10f}"


def dec_up(q):
    return f"{math.ceil(q * 10 ** 10) / 10 ** 10:.10f}"


def family(n):
    if n[0] == "J":
        return "Pal-Vertesi J list (4422)"
    if n.startswith("I4422"):
        return "Brunner-Gisin I4422 list"
    if n[0] == "A" and n[1:].isdigit():
        return "Ito-Imai-Avis A list (up to 5522)"
    return "CHSH, I3322, AS1, AS2, AII1, AII2"


def order_key(n):  # the order of Appendix A (paper/make_tables.py)
    import re
    fam = {"C": 0, "I": 1, "A": 3, "J": 2}.get(n[0], 4)
    if n.startswith(("AS", "AII")):
        fam = 4
    if n == "I3322":
        fam = 0
    num = int(re.findall(r"\d+", n)[-1]) if re.findall(r"\d+", n) else 0
    return (fam, n.rstrip("0123456789"), num)


def exact_lower(args):
    name, path = args
    cat = json.loads((P / "data" / "catalog.json").read_text())["inequalities"]
    return name, X.check(json.loads((P / path).read_text()), cat)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    cat = json.loads((P / "data" / "catalog.json").read_text())["inequalities"]
    recs = json.loads((P / "data" / "records.json").read_text())["records"]
    table = json.loads((P / "results" / "certified_table.json").read_text())
    windows = json.loads((P / "results" / "windows.json").read_text())
    con = sqlite3.connect(P / "board" / "board.db")
    board = {f: (i, json.loads(r)) for i, f, r in con.execute("SELECT id, file, report FROM submissions ORDER BY id")}

    for p in R.rglob("*"):  # the repository freezes data/ and verifier/ read-only (files and directories)
        p.chmod(0o755 if p.is_dir() else 0o644)
    for d in GENERATED:
        shutil.rmtree(R / d, ignore_errors=True)
    # ---- data and verifier
    shutil.copytree(P / "data", R / "data", ignore=shutil.ignore_patterns("__pycache__"))
    for f in VERIFIER:
        (R / "verifier" / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(P / "verifier" / f, R / "verifier" / f)
    shutil.copytree(P / "verifier" / "tests" / "fixtures", R / "verifier" / "tests" / "fixtures")
    for f in SEARCH:
        (R / "search" / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(P / f, R / "search" / f)
    shutil.copy2(__file__, R / "generate" / "build_release.py")
    for p in R.rglob("*"):  # the copies stay writable
        p.chmod(0o755 if p.is_dir() or p.suffix == ".sh" else 0o644)

    # ---- exact values of all strategies (independent checker), in parallel
    lowers = [(n, r["lower_file"]) for n, r in table["table"].items() if r["lower_file"]]
    with mp.Pool(a.jobs) as pool:
        exact = dict(pool.map(exact_lower, lowers, chunksize=1))

    # ---- certificates and manifest
    entries = []
    for row, name in enumerate(sorted(table["table"], key=order_key), 1):
        t, rec, ineq = table["table"][name], recs.get(name, {}), cat[name]
        pub = rec.get("pv2010", {}).get("lower", rec.get("lower"))
        if name == "CHSH":
            pub = "0.2071068"
        e = {"row": row, "name": name, "family": family(name),
             "coefficients": {k: ineq[k] for k in ("mA", "mB", "bA", "bB", "bAB", "b0")},
             "catalog_source": ineq.get("source"), "published": pub,
             "reference": {k: rec.get(k) for k in ("lower", "upper", "status", "source")} if rec else None,
             "window": {"lower": t["lower"], "upper": t["upper"]}, "lower": None, "upper": None}
        if t["lower"] and t["upper"]:
            e["window"]["width"] = str(Fraction(t["upper"]) - Fraction(t["lower"]))
        d = R / "certificates" / name
        if t["lower_file"]:
            d.mkdir(parents=True, exist_ok=True)
            shutil.copy2(P / t["lower_file"], d / "lower.json")
            v = Fraction(exact[name]["violation_exact"])
            bid, rep = board.get(t["lower_file"], (None, {}))
            if rep:
                assert rep["violation"] == dec_down(v), (name, rep["violation"], dec_down(v))
                if rep.get("violation_exact"):
                    assert Fraction(rep["violation_exact"]) == v, name
            assert dec_down(v) == t["lower"], (name, t["lower"], dec_down(v))
            c = json.loads((P / t["lower_file"]).read_text())
            e["lower"] = {"file": f"certificates/{name}/lower.json", "sha256": sha256(d / "lower.json"),
                          "board_file": t["lower_file"], "board_id": bid,
                          "field": c.get("field", "real"), "dims": [c["dA"], c["dB"]],
                          "ranks": [[len(s) for s in c["A"]], [len(s) for s in c["B"]]],
                          "violation_exact": str(v), "violation_down10": dec_down(v)}
        if t["upper_file"]:
            d.mkdir(parents=True, exist_ok=True)
            shutil.copy2(P / t["upper_file"], d / "upper.json")
            bid, rep = board[t["upper_file"]]
            c = json.loads((P / t["upper_file"]).read_text())
            lam, n, _, _ = X.sos_lambda(c, ineq)
            tau = Fraction(rep["sigma"]) + Fraction(rep["residual_bound"])
            bound = lam + n * tau - Fraction(ineq["b0"])
            assert n == rep["n_words"] and dec_up(lam) == rep["lambda"], name
            assert dec_up(bound) == rep["bound"] == t["upper"], (name, rep["bound"], dec_up(bound), t["upper"])
            if rep.get("bound_exact"):
                assert Fraction(rep["bound_exact"]) == bound, name
            e["upper"] = {"file": f"certificates/{name}/upper.json", "sha256": sha256(d / "upper.json"),
                          "board_file": t["upper_file"], "board_id": bid, "n_words": n,
                          "lambda_exact": str(lam), "sigma_hex": float(rep["sigma"]).hex(), "sigma": rep["sigma"],
                          "residual_hex": float(rep["residual_bound"]).hex(), "tau_stated_exact": str(tau),
                          "tau_stated": float(tau), "penalty": float(n * tau), "bound_exact": str(bound),
                          "bound_up10": dec_up(bound),
                          "verifier_exact_pd_2tau": rep.get("exact_pd_crosscheck")}
        entries.append(e)
    counts = table["counts"]
    manifest = {
        "release": "tsirelson-bounds-v2",
        "paper": "Certified Tsirelson bounds for 229 tight bipartite Bell inequalities",
        "convention": "violation = <psi|W|psi> - b0 with W = sum_x bA[x] A_x + sum_y bB[y] B_y + "
                      "sum_xy bAB[x][y] A_x B_y (Collins-Gisin form, projectors A_x, B_y); local bound 0",
        "lower_bound": "exact strategy: projectors V (V^dag V)^-1 V^dag from the given spanning vectors, state psi; "
                       "violation_exact is the exact rational value",
        "upper_bound": "sum_ij Y_ij w_i^* w_j = lambda 1 - W in the universal projector algebra and "
                       "Y + tau_stated I >= 0 imply <W> <= lambda + n_words * tau_stated for every commuting-"
                       "operator strategy; bound_exact = lambda + n_words * tau_stated - b0, tau_stated = sigma + "
                       "residual (binary64 values given exactly in hex)",
        "printed_values": "window.lower = violation rounded down to 1e-10, window.upper = bound rounded up to 1e-10",
        "counts": {"inequalities": len(entries), "lower": sum(1 for e in entries if e["lower"]),
                   "upper": sum(1 for e in entries if e["upper"]), "both": counts["both"],
                   "width_below_1e-7": counts["lt1e-7"], "width_below_1e-8": counts["lt1e-8"]},
        "inequalities": entries}
    (R / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    subprocess.run([sys.executable, str(R / "generate" / "make_appendix.py")], check=True)
    if (R / "appendix_A.tex").read_text() != (P / "paper" / "tab_catalog.tex").read_text():
        sys.exit("appendix_A.tex differs from paper/tab_catalog.tex")

    # ---- Table I (open cases)
    lines = ["inequality\tliterature_lower\tliterature_upper\tlower\tupper\tlower_certified\tupper_certified"]
    tol = Fraction(1, 10 ** 7)
    for n in ["J97", "I4422_20", "J34", "A65", "J73", "J62", "J30", "J48", "A68", "A89", "A62", "A82", "A84",
              "A21", "A14", "A47", "A64", "A67", "A80"]:
        w = windows[n]
        lo0, up0 = w["literature"]
        cl = w["certified_lower"] is not None and Fraction(w["certified_lower"]) > Fraction(lo0) - tol
        cu = w["certified_upper"] is not None and Fraction(w["certified_upper"]) < Fraction(up0) + tol
        lines.append(f"{n}\t{lo0}\t{up0}\t{w['now'][0]}\t{w['now'][1]}\t{'yes' if cl else 'no'}\t{'yes' if cu else 'no'}")
    (R / "table_I.tsv").write_text("\n".join(lines) + "\n")

    # ---- ladders: certified intermediate rungs
    lad = (R / "ladders" / "ladders_source.tsv").read_text().splitlines()
    out = [lad[0] + "\tcertificate\tsha256"]
    (R / "ladders" / "certificates").mkdir(parents=True, exist_ok=True)
    for line in lad[1:]:
        f = line.split("\t")
        board_file = f[-1]
        cert, h = "", ""
        if board_file:
            bid, rep = board[board_file]
            assert rep["violation"] == f[5] and rep["dims"] == [int(f[3])] * 2, (line, rep["violation"], rep["dims"])
            if board_file == table["table"][f[0]]["lower_file"]:
                cert = f"certificates/{f[0]}/lower.json"
            else:
                cert = f"ladders/certificates/{f[0]}_d{f[3]}.json"
                shutil.copy2(P / board_file, R / cert)
            h = sha256(R / cert)
        out.append(line + f"\t{cert}\t{h}")
    (R / "ladders" / "ladders.tsv").write_text("\n".join(out) + "\n")
    (R / "ladders" / "notebooks").mkdir(parents=True, exist_ok=True)
    for nb in ("lower1", "lower2"):
        shutil.copy2(P / "runs" / nb / "NOTES.md", R / "ladders" / "notebooks" / f"{nb}_NOTES.md")

    # ---- A19 case study inputs
    d = R / "case_studies" / "A19" / "data"
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy2(P / "data" / "sources" / "naceur2025_alist" / "A19.txt", d / "naceur2025_A19.txt")
    for f, out in (("candidates/20261013-003441_certify_bf15394840d8bfac.json", "complex_qubit_lower.json"),
                   ("candidates/20261013-005508_certify_a64acd3b91866969.json", "upper_N202.json")):
        assert board[f][1]["inequality"] == "A19"
        shutil.copy2(P / f, d / out)
    reps = {}
    for f, out in ((table["table"]["A19"]["lower_file"], "certificates/A19/lower.json"),
                   (table["table"]["A19"]["upper_file"], "certificates/A19/upper.json"),
                   ("candidates/20261013-003441_certify_bf15394840d8bfac.json", "case_studies/A19/data/complex_qubit_lower.json"),
                   ("candidates/20261013-005508_certify_a64acd3b91866969.json", "case_studies/A19/data/upper_N202.json")):
        bid, rep = board[f]
        reps[out] = {"board_id": bid, "board_file": f, "report": {k: v for k, v in rep.items() if k not in ("sources",)}}
    (d / "board_reports.json").write_text(json.dumps(reps, indent=1) + "\n")
    src = (P / "data" / "sources" / "pv2009_tables.tex").read_text().splitlines()
    (d / "pv2009_table_row.tex").write_text("\n".join(l for l in src if "A_{19}" in l) + "\n")
    mx = (P / "data" / "sources" / "naceur2025_alist" / "maximal_violations_exact.txt").read_text()
    i = mx.index("A19 ")
    (d / "naceur2025_maximal_violations_A19.txt").write_text(mx[i:mx.index(" A20 ", i)] + "\n")

    # ---- SHA256SUMS (permissions normalised first: the frozen submissions are read-only in the repository)
    for p in R.rglob("*"):
        p.chmod(0o755 if p.is_dir() or p.suffix == ".sh" else 0o644)
    files = sorted(p for p in R.rglob("*") if p.is_file() and p.name != "SHA256SUMS"
                   and "__pycache__" not in p.parts and p.name != ".DS_Store")
    (R / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(R).as_posix()}\n" for p in files))
    print(f"release: {len(entries)} inequalities, {manifest['counts']}, {len(files)} files")


if __name__ == "__main__":
    main()
