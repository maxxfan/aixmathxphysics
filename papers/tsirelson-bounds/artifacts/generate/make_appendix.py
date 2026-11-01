#!/usr/bin/env python3
"""Regenerate the paper's Appendix A table (appendix_A.tex, LaTeX rows, two entries per row) and its plain-text
form (appendix_A.tsv) from manifest.json.  Standard library only.

usage:  python3 generate/make_appendix.py            write appendix_A.tex and appendix_A.tsv
        python3 generate/make_appendix.py --check    compare with the files in the release (exit 1 on a difference)
"""
import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sci(x):
    if x is None:
        return "--"
    if x == 0:
        return "0"
    e = int(f"{x:.1e}".split("e")[1])
    m = x / 10 ** e
    return f"${m:.1f}\\times10^{{{e}}}$"


def plain(n):
    return "I" + n.split("_")[1] if n.startswith("I4422_") else n


def render(manifest):
    cells, tsv = [], ["row\tinequality\tpublished\tlower\tupper\twidth"]
    for e in manifest["inequalities"]:
        w = e["window"]
        width = None
        if w["lower"] and w["upper"]:
            width = Fraction(w["upper"]) - Fraction(w["lower"])
        cells.append(f"{plain(e['name'])} & {e['published']} & {w['lower'] or '--'} & {w['upper'] or '--'} & "
                     f"{sci(float(width)) if width is not None else '--'}")
        tsv.append(f"{e['row']}\t{e['name']}\t{e['published']}\t{w['lower'] or ''}\t{w['upper'] or ''}\t"
                   f"{width if width is not None else ''}")
    half = (len(cells) + 1) // 2
    rows = []
    for i in range(half):
        right = cells[half + i] if half + i < len(cells) else "& & & &"
        rows.append(f"{cells[i]} & {right} \\\\")
    return "\n".join(rows) + "\n", "\n".join(tsv) + "\n"


def main():
    manifest = json.loads((ROOT / "manifest.json").read_text())
    tex, tsv = render(manifest)
    if "--check" in sys.argv:
        ok = (ROOT / "appendix_A.tex").read_text() == tex and (ROOT / "appendix_A.tsv").read_text() == tsv
        print("appendix_A.tex and appendix_A.tsv " + ("match manifest.json" if ok else "DIFFER from manifest.json"))
        return 0 if ok else 1
    (ROOT / "appendix_A.tex").write_text(tex)
    (ROOT / "appendix_A.tsv").write_text(tsv)
    print(f"wrote appendix_A.tex and appendix_A.tsv ({len(manifest['inequalities'])} inequalities)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
