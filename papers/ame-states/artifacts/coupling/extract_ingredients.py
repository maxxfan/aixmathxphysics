#!/usr/bin/env python3
"""Maintainer script (run inside the repository, not part of the replay): copy the coupling certificates from
candidates/ into certificates/, split out the distinct ingredient states into ingredients/<id>.json, and write
cells.json (cell -> certificate -> ingredient ids). Distinct ingredients are identified by their canonical row space
(RREF over F_p), so the same state used in several certificates gets one id.

usage: python3 release/coupling/extract_ingredients.py
"""
import json
import pathlib
import shutil
import sys

P = pathlib.Path(__file__).resolve().parents[2]
OUT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(P / "verifier"))
import verify  # noqa: E402  (only rref_mod_p and parse_row, to identify equal row spaces)

cells = json.loads((P / "paper" / "phase3_cells.json").read_text())["cells"]
coupling = [c for c in cells if c["type"] == "coupling"]

NAMES = {  # (n, d, note-prefix) -> readable id
    "dodecacode": "a_dodecacode_n12_d2",
    "extended ternary Golay": "a_ternary_golay_n12_d3",
    "[[14,0,7]]_4": "a_14_0_7_d4_hermitian_circulant",
    "[[16,0,8]]_4": "a_16_0_8_d4_bordered_circulant",
    "AME(6,2) wheel": "a_wheel_ame6_2",
}


def ing_id(ing):
    note = ing.get("note", "")
    for key, name in NAMES.items():
        if key in note:
            return name
    n, d = ing["n"], ing["d"]
    tag = note.split(" ")[1] if note.startswith("AME(") and len(note.split(" ")) > 1 else "x"
    kind = "rs" if " RS" in note or note.endswith("RS") else (
        "bevins_bidav" if "Bevins" in note else "hz" if "HZ" in note else "circulant")
    return f"b_ame{n}_{d}_{kind}"


def rref_key(ing):
    p = ing["p"]
    rows = [verify.parse_row(r, p, len(ing["generators"][0]) if isinstance(ing["generators"][0], list) else len(ing["generators"][0]), i)
            for i, r in enumerate(ing["generators"])]
    R, _ = verify.rref_mod_p(rows, p)
    return json.dumps(R)


ingredients, mapping = {}, []
for c in sorted(coupling, key=lambda x: (x["n"], x["d"])):
    src = P / c["file"]
    cand = json.loads(src.read_text())
    dst = OUT / "certificates" / f"AME_{c['n']}_{c['d']}.json"
    shutil.copyfile(src, dst)
    ids = {}
    for role in ("a", "b"):
        ing = cand[role]
        key = rref_key(ing)
        name = ing_id(ing)
        prev = ingredients.get(name)
        if prev and prev["key"] != key:
            name = name + "_v2"
        if name not in ingredients:
            ingredients[name] = {"key": key, "obj": ing, "used_in": []}
        ingredients[name]["used_in"].append(f"AME({c['n']},{c['d']})")
        ids[role] = name
    mapping.append({"cell": f"AME({c['n']},{c['d']})", "n": c["n"], "d": c["d"], "k": cand["k"],
                    "verdict_at_submission": c["verdict"], "certificate": f"certificates/{dst.name}",
                    "board_candidate": c["file"], "a": ids["a"], "b": ids["b"],
                    "a_dim": cand["a"]["d"], "b_dim": cand["b"]["d"]})
for name, rec in ingredients.items():
    obj = dict(rec["obj"])
    obj["id"] = name
    obj["used_in"] = rec["used_in"]
    (OUT / "ingredients" / f"{name}.json").write_text(json.dumps(obj) + "\n")
(OUT / "cells.json").write_text(json.dumps({"cells": mapping}, indent=1) + "\n")
print(len(mapping), "certificates,", len(ingredients), "distinct ingredients")
for name, rec in sorted(ingredients.items()):
    o = rec["obj"]
    print(f"  {name:40s} n={o['n']:2d} d={o['d']:2d} p={o['p']:2d}  used in {len(rec['used_in'])}  | {o.get('note','')[:90]}")
