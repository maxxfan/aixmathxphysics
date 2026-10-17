#!/usr/bin/env python3
"""Negative controls for check_release.py.  Each case copies the whole bundle, applies one corruption,
re-signs MANIFEST.sha256 (so that only the mathematical and structural checks can catch it, except in
the 'manifest' case) and requires the checker to reject it; the unchanged copy must pass.

usage: python3 -I negative_controls.py [RELEASE_DIR]   (works in a temporary directory)
"""
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile


def resign(d):
    files = sorted(p for p in d.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    (d / "MANIFEST.sha256").write_text("".join(
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(d)}\n" for p in files))


def qasm(d, n):
    return next((d / "circuits").glob(f"gf2_{n}_T*.qasm"))


def record(d, n):
    return next((d / "records").glob(f"gf2_{n}_T*.json"))


def edit_lines(path, fn):
    path.write_text("\n".join(fn(path.read_text().splitlines())) + "\n")


def first(lines, prefix):
    return next(k for k, line in enumerate(lines) if line.startswith(prefix))


def last(lines, prefix):
    return max(k for k, line in enumerate(lines) if line.startswith(prefix))


def drop(prefix, where=first):
    return lambda L: L[:where(L, prefix)] + L[where(L, prefix) + 1:]


def t_to_tdg(L):
    i = first(L, "t q[")
    return L[:i] + ["tdg" + L[i][1:]] + L[i + 1:]


def header_inject(L):
    i = first(L, "include")
    return L[:i] + [L[i] + " x q[0];"] + L[i + 1:]


def other_include(L):
    i = first(L, "include")
    return L[:i] + ['include "other.inc";'] + L[i + 1:]


def json_edit(path, fn):
    r = json.loads(path.read_text())
    fn(r)
    path.write_text(json.dumps(r))


def flip_bit(r):
    s = r["factors"][3]
    r["factors"][3] = s[:5] + ("1" if s[5] == "0" else "0") + s[6:]


def fullwidth(L):
    i = first(L, "t q[")
    j = L[i].index("[") + 1
    return L[:i] + [L[i][:j] + chr(0xFF10 + int(L[i][j])) + L[i][j + 1:]] + L[i + 1:]


def table_bar(d, n, t):
    p = d / "data" / "records.json"
    p.chmod(0o644)
    tab = json.loads(p.read_text())
    tab["records"][str(n)]["tcount"] = t
    p.write_text(json.dumps(tab))


def empty(d):
    for x in ("records", "circuits", "blocks"):
        shutil.rmtree(d / x)
        (d / x).mkdir()


def block(d, stem):
    return d / "blocks" / f"{stem}.json"


CASES = [
    ("unchanged control", lambda d: None, 0),
    ("drop one cz from the Clifford correction", lambda d: edit_lines(qasm(d, 7), drop("cz ", last)), 1),
    ("change one t into tdg", lambda d: edit_lines(qasm(d, 32), t_to_tdg), 1),
    ("drop one cx from a ladder", lambda d: edit_lines(qasm(d, 7), drop("cx ")), 1),
    ("flip one bit of one factor", lambda d: json_edit(record(d, 7), flip_bit), 1),
    ("swap in another irreducible modulus (n=7)", lambda d: json_edit(record(d, 7), lambda r: r.update(modulus="10001001")), 1),
    ("statement hidden on the include line", lambda d: edit_lines(qasm(d, 7), header_inject), 1),
    ("include another gate library", lambda d: edit_lines(qasm(d, 7), other_include), 1),
    ("remove one record", lambda d: (record(d, 13).unlink(), qasm(d, 13).unlink()), 1),
    ("empty bundle (no records, no blocks)", lambda d: empty(d), 1),
    ("block T-count field wrong", lambda d: json_edit(block(d, "block_ring_e2"), lambda r: r.update(tcount=3)), 1),
    ("ring F2[x]/(x^3+1) relabelled GF(8)", lambda d: json_edit(block(d, "block_alpha_n3"),
                                                               lambda r: r.update(modulus="1001")), 1),
    ("form-feed separator hides a gate in a comment", lambda d: edit_lines(qasm(d, 7), lambda L: L[:4] + ["// x\x0ch q[10];"] + L[4:]), 1),
    ("fullwidth digits in a qubit index", lambda d: edit_lines(qasm(d, 7), fullwidth), 1),
    ("gate_counts field wrong", lambda d: json_edit(record(d, 7), lambda r: r["gate_counts"].update(CNOT=1)), 1),
    ("clifford_correction field wrong", lambda d: json_edit(record(d, 7), lambda r: r["clifford_correction"].update(CZ=[])), 1),
    ("duplicate record for the same n", lambda d: (shutil.copy(record(d, 7), d / "records" / "gf2_7_T97b.json"),
                                                   shutil.copy(qasm(d, 7), d / "circuits" / "gf2_7_T97b.qasm")), 1),
    ("modulus with a sign ('-100101')", lambda d: json_edit(record(d, 5), lambda r: r.update(modulus="-100101")), 1),
    ("factor with an underscore", lambda d: json_edit(record(d, 7), lambda r: r["factors"].__setitem__(0, r["factors"][0][:3] + "_" + r["factors"][0][4:])), 1),
    ("leading zero in a qubit index", lambda d: edit_lines(qasm(d, 7), lambda L: L[:-1] + [L[-1].replace("q[", "q[0", 1)]), 1),
    ("modulus_poly text wrong", lambda d: json_edit(record(d, 5), lambda r: r.update(modulus_poly="x^5 + x + 1")), 1),
    ("comparison value changed", lambda d: json_edit(record(d, 6), lambda r: r["comparison"].update(reported=79)), 1),
    ("tcount as a float", lambda d: json_edit(record(d, 7), lambda r: r.update(tcount=97.0)), 1),
    ("record bar changed in both record and table", lambda d: (json_edit(record(d, 5), lambda r: r["published_record"].update(tcount=49)),
                                                              table_bar(d, 5, 49)), 1),
    ("file missing from the manifest", "unlisted", 1),
    ("file changed without updating the manifest", None, 1),
]


def main():
    src = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parent)
    bad = 0
    with tempfile.TemporaryDirectory() as tmp:
        for k, (name, mutate, want) in enumerate(CASES):
            d = pathlib.Path(tmp) / f"case{k}"
            shutil.copytree(src, d)
            if mutate is None and want:
                p = qasm(d, 7)
                p.write_text(p.read_text() + "// tampered\n")      # no re-sign
            elif mutate == "unlisted":
                resign(d)
                (d / "extra.txt").write_text("not in the manifest\n")
            else:
                if mutate:
                    mutate(d)
                resign(d)
            r = subprocess.run([sys.executable, "-I", str(d / "check_release.py"), str(d)], capture_output=True, text=True)
            got = r.returncode
            ok = (got != 0) == bool(want)
            bad += not ok
            print(f"  {name:<45} checker exit {got}: {'as expected' if ok else 'UNEXPECTED'}")
    print("NEGATIVE CONTROLS OK" if not bad else f"{bad} NEGATIVE CONTROL(S) FAILED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
