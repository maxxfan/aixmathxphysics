#!/usr/bin/env python3
"""Build data/status.json (existence status of every AME(n,d) cell, 2 <= d <= 11) and data/targets.json.

Sources (frozen in data/sources/):
  huber_wyderka_table_2026-06-08.html   F. Huber and N. Wyderka, "Table of AME states" (huberfe.github.io/ame,
                                         raw table.html from github.com/huberfe/huberfe.github.io, gh-pages,
                                         page says "updated: 08.06.2026"; fetched 2026-10-16)
  codetables_k0_2026-10-10.json          codetables.de [[n,0,d]]_q bounds for q <= 8 (F_q-linear stabilizer codes)
Overlay (2026 preprints, read 2026-10-16; phase-1 results and AME(9,6) added after the phase-1 run; see PROBLEM.md):
  arXiv:2608.01011  Shi, Zhang, Zhao, Li: AME(7,d) exists for every d >= 3
  arXiv:2608.05781  Bevins, Bidav: AME(12,5), AME(18,11), AME(18,13) (+ projections AME(17,11), AME(17,13))
  arXiv:2609.11796  Hontarenko, Zyczkowski: for every odd prime power q >= 3: AME(2q,q), AME(2q-1,q),
                    AME(4k+2,q), AME(4k+1,q) for 1 <= k <= (q-1)/2, AME(q+3,q); Table 1 parents AME(16,9),
                    AME(15,9), AME(20,13), AME(19,13)
Closure rules applied to the union (both are standard):
  product   AME(n,d1) and AME(n,d2) => AME(n,d1*d2)                 (tensor product, party-wise)
  project   AME(2k,d) => AME(2k-1,d)                                (project one party onto a basis state)
            hence AME(2k-1,d) nonexistent => AME(2k,d) nonexistent

Blank cells of the Huber-Wyderka table are unknown (open). Status values: exists, nonexistent, open.
Run: python3 data/build_status.py   (stdlib only; writes data/status.json and data/targets.json)
"""
import html
import json
import math
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "sources"
DMAX = 11      # the Huber-Wyderka table covers 2 <= d <= 11
DEXT = 64      # phase 3: cells 12 <= d <= 64, 2 <= n <= NEXT, from published families only ("beyond_table")
NEXT = 40


def scott_ok(n, d):
    """Scott's bound: AME(n,d) needs n <= 2(d^2-1) (n even), n <= 2d(d+1)-1 (n odd)."""
    return n <= (2 * (d * d - 1) if n % 2 == 0 else 2 * d * (d + 1) - 1)


def parse_hw():
    src = (SRC / "huber_wyderka_table_2026-06-08.html").read_text(encoding="utf-8")
    line = next(l for l in src.split("\n") if '"exists":[' in l)
    start = line.find("{")
    depth, end, ins, esc = 0, None, False, False
    for i in range(start, len(line)):
        ch = line[i]
        if ins:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                ins = False
            continue
        if ch == '"':
            ins = True
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    doc = json.loads(line[start:end])
    found = []

    def walk(o):
        if isinstance(o, dict):
            if "exists" in o and "ident" in o:
                found.append(o)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(doc)
    data = found[0]
    cells = {}
    for i in range(len(data["ident"])):
        n, d = int(data["n"][i]), int(data["D"][i])
        refs = re.sub(r"<[^>]+>", " ", html.unescape(data["refs"][i]))
        cells[(n, d)] = {"hw": data["exists"][i], "const": data["const"][i], "comments": data["comments"][i],
                         "refs": " ".join(refs.split())}
    return cells


def is_prime_power(d):
    for p in range(2, d + 1):
        if d % p == 0:
            while d % p == 0:
                d //= p
            return d == 1
    return False


def overlay():
    ex = {}
    for d in range(3, DEXT + 1):
        ex[(7, d)] = "arXiv:2608.01011 (Shi, Zhang, Zhao, Li 2026): AME(7,d) for all d >= 3"
    for nd in [(12, 5), (18, 11), (17, 11), (18, 13), (17, 13)]:
        ex[nd] = "arXiv:2608.05781 (Bevins, Bidav 2026): Hermitian self-dual MDS codes (and projection)"
    for d in range(3, DEXT + 1):  # orthogonal Latin squares of every order d != 2, 6 (Bose-Shrikhande-Parker 1960); d = 6: golden AME(4,6)
        ex.setdefault((4, d), "AME(4,d) <=> a pair of orthogonal Latin squares of order d (exists for d != 2, 6); AME(4,6): Rather et al. PRL 128, 080507 (2022)")
    for q in range(DMAX + 1, DEXT + 1):  # beyond the table: minimal-support AME from extended Reed-Solomon codes
        if is_prime_power(q):
            for n in range(2, q + 2):
                ex.setdefault((n, q), "MDS (extended Reed-Solomon) code [n, n//2] over GF(q), n <= q+1: minimal-support AME")
    for q in range(3, DEXT + 1, 2):
        if not is_prime_power(q):
            continue
        tag = "arXiv:2609.11796 (Hontarenko, Zyczkowski 2026)"
        ex[(2 * q, q)] = tag + ": AME(2q,q), odd prime power q"
        ex[(2 * q - 1, q)] = tag + ": AME(2q-1,q) by projection"
        for k in range(1, (q - 1) // 2 + 1):
            ex.setdefault((4 * k + 2, q), tag + ": AME(4k+2,q)")
            ex.setdefault((4 * k + 1, q), tag + ": AME(4k+1,q)")
        ex.setdefault((q + 3, q), tag + ": AME(q+3,q)")
    ex[(9, 6)] = ("K. Kitamura, github.com/KitaKen1/ame-9-6-lean (Lean 4 proof, 2026-09-22; cited in arXiv:2610.17170); "
                  "non-stabilizer, re-checked numerically by the p10 auditor (board #58)")
    ex[(15, 8)] = "p10 phase 1 (this work, 2026-10-16): GF(8)-linear symmetric circulant graph state, board #11, audited"
    ex[(17, 8)] = "p10 phase 1 (this work, 2026-10-16): GF(8)-linear symmetric circulant graph state, board #13, audited"
    ex[(11, 10)] = ("K. Kitamura, github.com/KitaKen1/ame-11-10-lean (Lean 4); coupling of the dodecacode with AME(12,5), "
                    "re-derived by the p10 auditor (argument.md 5.11, board #317)")
    ex[(20, 13)] = "arXiv:2609.11796 Table 1: explicit parent AME(20,13)"
    ex[(19, 13)] = "arXiv:2609.11796 Table 1: projection AME(19,13)"
    ex[(16, 9)] = "arXiv:2609.11796 Table 1: explicit parent AME(16,9)"
    ex[(15, 9)] = "arXiv:2609.11796 Table 1: projection AME(15,9)"
    return ex


def main():
    hw = parse_hw()
    ct = json.loads((SRC / "codetables_k0_2026-10-10.json").read_text())["bounds"]
    ov = overlay()
    NMAX = {d: (2 * d * (d + 1) - 1) for d in range(2, DMAX + 1)}
    st = {}
    for d in range(2, DMAX + 1):
        for n in range(2, NMAX[d] + 1):
            c = hw.get((n, d))
            if not scott_ok(n, d):
                s, src = "nonexistent", "Scott bound"
            elif c is None or c["hw"] == "Unknown":
                s, src = "open", "Huber-Wyderka 2026-06-08: " + ("unknown" if c else "blank")
            elif c["hw"] == "Yes":
                s, src = "exists", "Huber-Wyderka 2026-06-08 (constructions: " + (c["const"] or "not listed") + ")"
            else:
                s, src = "nonexistent", "Huber-Wyderka 2026-06-08: " + c["hw"]
            if (n, d) in ov and s != "exists":
                if s == "nonexistent":
                    raise SystemExit(f"overlay contradicts table at AME({n},{d})")
                s, src = "exists", ov[(n, d)]
            st[(n, d)] = [s, src]
    for d in range(DMAX + 1, DEXT + 1):  # beyond the table: published families only
        for n in range(2, NEXT + 1):
            if not scott_ok(n, d):
                st[(n, d)] = ["nonexistent", "Scott bound"]
            elif (n, d) in ov:
                st[(n, d)] = ["exists", ov[(n, d)]]
            else:
                st[(n, d)] = ["open", "beyond the Huber-Wyderka table (d <= 11); not covered by the published families in data/build_status.py"]
    changed = True
    while changed:  # closure: products and even->odd projection
        changed = False
        for (n, d), (s, _) in list(st.items()):
            if s == "open":
                for d1 in range(2, d):
                    if d % d1 == 0 and st.get((n, d1), ["?"])[0] == "exists" and st.get((n, d // d1), ["?"])[0] == "exists":
                        st[(n, d)] = ["exists", f"product AME({n},{d1}) x AME({n},{d // d1})"]
                        changed = True
                        break
            if s == "open" and n % 2 == 1 and st.get((n + 1, d), ["?"])[0] == "exists":
                st[(n, d)] = ["exists", f"projection of AME({n + 1},{d})"]
                changed = True
            if s == "open" and n % 2 == 0 and st.get((n - 1, d), ["?"])[0] == "nonexistent":
                st[(n, d)] = ["nonexistent", f"AME({n - 1},{d}) does not exist (projection)"]
                changed = True
    cells = []
    for (n, d), (s, src) in sorted(st.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        row = {"n": n, "d": d, "status": s, "source": src}
        if d > DMAX:
            row["beyond_table"] = True
        if (n, d) == (8, 4):
            row["stabilizer"] = ("excluded for every group-type nice error basis: no F_2-linear [[8,0,5]]_4 "
                                 "(Ball, Moreno, Simoens, arXiv:2401.06618; p10 phase 2 classify/auditor) and Z_4^2 reduces to AME(8,2)")
        if str(d) in ct and str(n) in ct[str(d)]:
            lo, hi = ct[str(d)][str(n)]
            need = n // 2 + 1
            row["codetables"] = [lo, hi]
            row["stabilizer_Fq_linear"] = "excluded" if hi < need else ("exists" if lo >= need else "open")
            if row["stabilizer_Fq_linear"] == "exists" and s != "exists":
                print(f"warning: codetables has [[{n},0,{lo}]]_{d} but status is {s}")
        cells.append(row)
    meta = {"built_by": "data/build_status.py", "sources": ["sources/huber_wyderka_table_2026-06-08.html",
            "sources/codetables_k0_2026-10-10.json", "arXiv:2608.01011", "arXiv:2608.05781", "arXiv:2609.11796"],
            "dmax": DMAX, "dext": DEXT, "next": NEXT,
            "note": ("status of AME(n,d) for 2 <= d <= 11 and n up to Scott's bound (Huber-Wyderka + overlay); for 12 <= d <= 64 "
                     "and n <= 40 from published families only (RS n <= q+1, Shi n = 7, Hontarenko-Zyczkowski odd q, "
                     "Bevins-Bidav, Kitamura, product and projection closure), flagged beyond_table; cells outside are not in the snapshot")}
    (HERE / "status.json").write_text(json.dumps({"meta": meta, "cells": cells}, indent=0) + "\n")
    by = {(c["n"], c["d"]): c for c in cells}
    print("open cells (n <= 24):")
    for d in range(2, DEXT + 1):
        o = [c for c in cells if c["d"] == d and c["status"] == "open" and c["n"] <= 24]
        if not o:
            continue
        print(f"  d={d}: " + ", ".join(f"{c['n']}" + ("x" if c.get("stabilizer_Fq_linear") == "excluded" else "") for c in o))
    write_targets(by)


def write_targets(by):
    """Phase-3 targets (2026-10-18): the beyond-table frontier (q = 16 circulants; composite d by coupling), then the
    remaining in-table cells at low priority. Any open cell of data/status.json is a win."""
    T = []

    def add(n, d, score, band, why):
        c = by[(n, d)]
        if c["status"] != "open":
            print(f"skip target AME({n},{d}): status {c['status']}")
            return
        ctb = f"; codetables [[{n},0,d]]_{d} in [{c['codetables'][0]},{c['codetables'][1]}]" if "codetables" in c else ""
        T.append({"id": f"AME({n},{d})", "n": n, "d": d, "band": band, "score": score,
                  "display": f"AME({n},{d}){' [beyond table]' if c.get('beyond_table') else ''}: {why}{ctb}"})

    # q = 16: the even-q analogue of AME(q+3,q); symmetric circulants reach n = 19 (coordinator probe: 5 hits in 6e6 samples)
    add(19, 16, 6.0, "q16", "GF(16) symmetric circulant graph states (lambda 9.8); coordinator probe found it, classify the solutions")
    add(21, 16, 5.0, "q16", "lambda 33.8 for plain circulants; multiplier-symmetric circulants (Z_21^x orbits) predict ~15 hits")
    add(23, 16, 2.0, "q16", "lambda 118; needs large symmetry")
    add(18, 16, 2.0, "q16", "even n: lambda 87; Hermitian self-dual [18,9,10]_256 or bordered/Paley-type")
    add(20, 16, 1.0, "q16", "even n")
    for d in (13, 17, 19):
        add(23, d, 1.0, f"q{d}", "next cell beyond the Hontarenko-Zyczkowski families; large symmetry needed")
    # composite d by the coupling rule (argument.md 5.11): k-uniform (2k+2)-party G_a (dim a) + AME(2k+2, b), a <= b -> AME(2k+1, ab)
    for n, d, sc, why in [
        (11, 14, 5.0, "coupling: dodecacode (5-uniform 12 qubits) + AME(12,7)"),
        (11, 18, 4.0, "coupling: dodecacode + AME(12,9)"),
        (11, 22, 4.0, "coupling: dodecacode + AME(12,11)"),
        (11, 26, 3.0, "coupling: dodecacode + AME(12,13)"),
        (13, 20, 4.0, "coupling: 6-uniform 14-ququart [[14,0,7]]_4 + AME(14,5)"),
        (13, 28, 3.0, "coupling: [[14,0,7]]_4 + AME(14,7)"),
        (13, 15, 3.0, "coupling: 6-uniform 14-qutrit [[14,0,7]]_3 (if it exists) + AME(14,5)"),
        (13, 21, 3.0, "coupling: [[14,0,7]]_3 + AME(14,7)"),
        (11, 12, 3.0, "coupling: 5-uniform 12-qutrit [[12,0,6]]_3 + AME(12,4)? (AME(12,4) does not exist) -> needs another route"),
        (8, 6, 2.0, "non-stabilizer only; Kitamura-type quadratic family excluded (board #59)"),
    ]:
        add(n, d, sc, "coupling" if "coupling" in why else "nonstab", why)
    # in-table remainder (phase 1-2 closed most structured families; new ideas only)
    for n, d, sc, why in [
        (8, 4, 3.0, "non-stabilizer only: every group-type stabilizer is excluded (Ball-Moreno-Simoens 2024; phase 2)"),
        (16, 8, 2.0, "F_64-linear needs |Aut| <= 4 (involution lemma chain); F_8 circulants excluded (parity lemma)"),
        (19, 8, 2.0, "F_8 circulants exhausted (0 of 8^9)"),
        (20, 11, 2.0, "announced in arXiv:2610.17170 (forthcoming revision of 2609.11796); no construction public"),
        (19, 11, 2.0, "projection of AME(20,11)"),
        (16, 7, 1.0, "beyond 2q"), (16, 5, 1.0, "beyond 2q"), (14, 4, 1.0, "q = 4"), (13, 4, 1.0, "q = 4"),
        (15, 4, 1.0, "q = 4"), (18, 8, 1.0, "q = 8"), (20, 8, 1.0, "q = 8"), (19, 9, 1.0, "beyond 2q"), (20, 9, 1.0, "beyond 2q"),
    ]:
        add(n, d, sc, f"q{d}" if d != 4 or n != 8 else "marquee", why)
    meta = {"source": "data/status.json (Huber-Wyderka table 2026-06-08 + 2026 overlay; d > 11 from published families); "
                      "built by data/build_status.py",
            "snapshot": "2026-10-18 (phase 3: beyond-table frontier added; AME(11,10) Kitamura; stabilizer AME(8,4) excluded)",
            "win": "IMPROVEMENT = an exactly verified AME(n,d) state for ANY cell whose status is open in data/status.json "
                   "(within the verifier limits); beyond-table cells (d > 11) are wins pending the auditor's novelty check",
            "score_key": "bad_cuts (0 = AME)"}
    (HERE / "targets.json").write_text(json.dumps({"meta": meta, "targets": T}, indent=1) + "\n")
    print(f"{len(T)} targets written")


if __name__ == "__main__":
    main()
