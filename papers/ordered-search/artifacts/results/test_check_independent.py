#!/usr/bin/env python3
"""Tests for results/check_independent.py (needs python-flint).

  * the two small real certificates (verifier/tests/fixtures) are accepted with the margins the
    verifier reports;
  * malformed or tampered inputs are rejected with `Rejected`, never accepted and never with
    another exception;
  * the command line rejects tampered, truncated and non-JSON files with exit status 1, also
    under `python -O` (acceptance must not depend on assert statements).

usage: .venv/bin/python results/test_check_independent.py   (exit status 1 on any failure)
"""
import copy
import gzip
import json
import pathlib
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "verifier"))
import check_independent as ci  # noqa: E402
import verify  # noqa: E402

FIX = {name: json.loads((ROOT / "verifier" / "tests" / "fixtures" / f"{name}.json").read_text())
       for name in ("cert_k3_n56", "cert_k4_n605")}
failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


# 1. real certificates
for name, c in FIX.items():
    res = ci.check(copy.deepcopy(c))
    ref = verify.verify_candidate(copy.deepcopy(c))["margins"]
    check(res["k"] == c["k"] and res["n"] == c["n"] and not res["beats_prior_explicit"],
          f"{name}: accepted, k and n reproduced, does not beat (5, 7265)")
    check(all(abs(a - b) <= 1e-12 * abs(b) for a, b in zip(res["margins"], ref)),
          f"{name}: margins agree with the verifier {res['margins']}")


# 2. malformed and tampered inputs (built from the (4, 605) certificate)
def tamper(f):
    c = copy.deepcopy(FIX["cert_k4_n605"])
    out = f(c)
    return c if out is None else out


def setq(c, t, j, v):
    c["q"][t][j] = v


CASES = [
    ("not a dict", lambda c: [c]),
    ("wrong type", lambda c: c.update(type="other")),
    ("missing q", lambda c: c.pop("q") and None),
    ("missing h_den", lambda c: c.pop("h_den") and None),
    ("k is bool", lambda c: c.update(k=True)),
    ("k is float", lambda c: c.update(k=4.0)),
    ("k = 1", lambda c: c.update(k=1, q=[], h=[])),
    ("n mismatch", lambda c: c.update(n=121)),
    ("den not divisible by n", lambda c: c.update(den=c["den"] + 1)),
    ("den negative", lambda c: c.update(den=-c["den"])),
    ("h_den zero", lambda c: c.update(h_den=0)),
    ("h_den negative", lambda c: c.update(h_den=-c["h_den"])),
    ("den is float", lambda c: c.update(den=float(c["den"]))),
    ("q not a list", lambda c: c.update(q="q")),
    ("q too few rows", lambda c: c.update(q=c["q"][:-1])),
    ("q row too short", lambda c: c["q"][0].pop() and None),
    ("q entry float", lambda c: setq(c, 1, 3, float(c["q"][1][3]))),
    ("q entry bool", lambda c: setq(c, 1, 3, True)),
    ("q entry string", lambda c: setq(c, 1, 3, "7")),
    ("h empty", lambda c: c["h"].__setitem__(0, [])),
    ("h too long", lambda c: c["h"][0].extend([0, 0])),
    ("h entry float", lambda c: c["h"][2].__setitem__(5, 0.5)),
    ("constant term", lambda c: setq(c, 0, 0, c["q"][0][0] + 1)),
    ("root condition (+1 on one coefficient)", lambda c: setq(c, 1, 7, c["q"][1][7] + 1)),
    ("positivity (h scaled up 0.1%)", lambda c: c["h"].__setitem__(0, [v + v // 1000 for v in c["h"][0]])),
    ("positivity (h dropped)", lambda c: c["h"].__setitem__(2, [0])),
]
for name, f in CASES:
    bad = tamper(f)
    try:
        ci.check(bad)
        check(False, f"rejected: {name} (was ACCEPTED)")
    except ci.Rejected as e:
        check(True, f"rejected: {name} -> {e}")
    except Exception as e:  # noqa: BLE001 -- any other exception is a checker bug
        check(False, f"rejected: {name} (raised {type(e).__name__}: {e})")


# 3. command line, with and without -O
with tempfile.TemporaryDirectory() as tmp:
    tmp = pathlib.Path(tmp)
    good = tmp / "good.json"
    good.write_text(json.dumps(FIX["cert_k3_n56"]))
    files = {"root condition": tamper(lambda c: setq(c, 1, 7, c["q"][1][7] + 1)),
             "positivity": tamper(lambda c: c["h"].__setitem__(2, [0])),
             "den is float": tamper(lambda c: c.update(den=float(c["den"])))}
    paths = {}
    for name, c in files.items():
        p = tmp / f"{name.replace(' ', '_')}.json"
        p.write_text(json.dumps(c))
        paths[name] = p
    trunc = tmp / "truncated.json.gz"
    trunc.write_bytes(gzip.compress(json.dumps(FIX["cert_k3_n56"]).encode())[:-40])
    paths["truncated gzip"] = trunc
    junk = tmp / "junk.json"
    junk.write_text("{not json")
    paths["not JSON"] = junk
    for flags in ([], ["-O"]):
        tag = " ".join(flags) or "no flags"
        run = subprocess.run([sys.executable, *flags, str(HERE / "check_independent.py"), str(good)],
                             capture_output=True, text=True)
        check(run.returncode == 0 and "REJECTED" not in run.stdout, f"cli ({tag}): real certificate accepted")
        for name, p in paths.items():
            run = subprocess.run([sys.executable, *flags, str(HERE / "check_independent.py"), str(p)],
                                 capture_output=True, text=True)
            check(run.returncode == 1 and "REJECTED" in run.stdout and not run.stderr,
                  f"cli ({tag}): {name} rejected with exit status 1")

print(f"{len(failures)} failures")
sys.exit(1 if failures else 0)
