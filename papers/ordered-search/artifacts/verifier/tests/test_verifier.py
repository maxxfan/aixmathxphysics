#!/usr/bin/env python3
"""Tests for the p09 verifier (stdlib only; run before freeze.sh).

  * real certificates reproduce known thresholds: 3 queries for n = 56 (CLP07) and 4 queries for
    n = 605 (CLP07), both BELOW_RECORD and not beating the prior explicit (5, 7265);
  * every verdict path: IMPROVEMENT / MATCHES_RECORD (record moved by monkeypatching),
    VERIFIER_ALARM (published limit lowered), SCREEN_ONLY, INVALID (17 malformed or tampered inputs);
  * independent cross-checks: Kronecker autocorrelation against an O(n^2) loop on random signed
    vectors; the exact root conditions against complex evaluation of q_t - q_{t-1} at the roots;
    the certified polynomials against a float evaluation on a fine grid.
"""
import cmath
import copy
import json
import math
import pathlib
import random
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import verify  # noqa: E402

FIX = {name: json.loads((HERE / "fixtures" / f"{name}.json").read_text())
       for name in ("cert_k3_n56", "cert_k4_n605")}
failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


def verdict(c, **kw):
    return verify.verify_candidate(c, **kw)


# 1. known thresholds
for name, (k, n) in (("cert_k3_n56", (3, 56)), ("cert_k4_n605", (4, 605))):
    r = verdict(FIX[name])
    check(r["verdict"] == "BELOW_RECORD" and r["k"] == k and r["n"] == n,
          f"{name}: {r['summary']}")
    check(all(m > 0 for m in r["margins"]), f"{name}: positive certified margins {r['margins']}")
    check(r["target"] == f"k{k}", f"{name}: target {r['target']}")
    check(r["beats_prior_explicit"] is False, f"{name}: does not beat the prior explicit (5, 7265)")

# 2. verdict paths via the record / published limits
saved_rec, saved_pub = dict(verify.RECORD), dict(verify.PUBLISHED_MAX)
try:
    verify.RECORD.update(k=4, n=605)
    check(verdict(FIX["cert_k4_n605"])["verdict"] == "MATCHES_RECORD", "MATCHES_RECORD when record = (4, 605)")
    verify.RECORD.update(k=4, n=600)
    check(verdict(FIX["cert_k4_n605"])["verdict"] == "IMPROVEMENT", "IMPROVEMENT when record = (4, 600)")
    verify.RECORD.update(k=4, n=606)
    check(verdict(FIX["cert_k4_n605"])["verdict"] == "BELOW_RECORD", "BELOW_RECORD when record = (4, 606)")
    verify.RECORD.clear()
    verify.RECORD.update(saved_rec)
    verify.PUBLISHED_MAX[3] = (50, "test")
    check(verdict(FIX["cert_k3_n56"])["verdict"] == "VERIFIER_ALARM", "VERIFIER_ALARM above a published limit")
finally:
    verify.RECORD.clear()
    verify.RECORD.update(saved_rec)
    verify.PUBLISHED_MAX.clear()
    verify.PUBLISHED_MAX.update(saved_pub)
check(verdict(FIX["cert_k3_n56"], quick=True)["verdict"] == "SCREEN_ONLY", "SCREEN_ONLY with quick=True")

# 3. malformed and tampered inputs
base = FIX["cert_k4_n605"]


def tamper(f):
    c = copy.deepcopy(base)
    f(c)
    return c


bad = {
    "not a dict": [1, 2],
    "wrong type": tamper(lambda c: c.update(type="other")),
    "k = 1": tamper(lambda c: c.update(k=1)),
    "k is bool": tamper(lambda c: c.update(k=True)),
    "n mismatch": tamper(lambda c: c.update(n=121)),
    "den not divisible by n": tamper(lambda c: c.update(den=c["den"] + 1)),
    "negative h_den": tamper(lambda c: c.update(h_den=-c["h_den"])),
    "missing q": tamper(lambda c: c.pop("q")),
    "q too short": tamper(lambda c: c["q"][0].pop()),
    "q float entry": tamper(lambda c: c["q"][1].__setitem__(3, 1.5)),
    "h too long": tamper(lambda c: c["h"][0].extend([0, 0])),
    "constant term": tamper(lambda c: c["q"][0].__setitem__(0, c["q"][0][0] + 1)),
    "root condition (+1 on one coefficient)": tamper(lambda c: c["q"][1].__setitem__(7, c["q"][1][7] + 1)),
    "root condition on q_{k-1}": tamper(lambda c: c["q"][2].__setitem__(10, c["q"][2][10] - 2)),
    "consistent pair change breaks the next step": tamper(
        lambda c: (c["q"][0].__setitem__(5, c["q"][0][5] + 3), c["q"][0].__setitem__(600, c["q"][0][600] + 3))),
    "positivity (h scaled up 0.1%)": tamper(lambda c: c["h"].__setitem__(0, [v + v // 1000 for v in c["h"][0]])),
    "positivity (h dropped)": tamper(lambda c: c["h"].__setitem__(2, [0])),
}
for name, c in bad.items():
    r = verdict(c)
    check(r["verdict"] == "INVALID", f"INVALID: {name} -> {r.get('reason', r['verdict'])[:90]}")

# 4. independent cross-checks
rng = random.Random(7)
for trial in range(30):
    m = rng.randint(1, 40)
    h = [rng.randint(-2 ** 70, 2 ** 70) for _ in range(m)]
    brute = [sum(h[i] * h[i + j] for i in range(m - j)) for j in range(m)]
    if verify.autocorrelation(h) != brute:
        check(False, f"autocorrelation trial {trial}")
        break
else:
    check(True, "Kronecker autocorrelation == O(n^2) loop on 30 random signed vectors")

for name in FIX:
    c = FIX[name]
    k, n, D = c["k"], c["n"], c["den"]
    chain = [[(n - j) / n for j in range(n)]] + [[v / D for v in q] for q in c["q"]] + [[1.0] + [0.0] * (n - 1)]
    worst = 0.0
    for t in range(1, k + 1):
        s = (-1) ** t
        for r in range(0, n, max(1, n // 25)):
            z = cmath.exp(1j * (math.pi * (2 * r + (0 if s == 1 else 1)) / n))
            d = [chain[t][j] - chain[t - 1][j] for j in range(n)]
            val = d[0] + sum(d[j] * (z ** j + z ** (-j)) for j in range(1, n))
            worst = max(worst, abs(val))
    check(worst < 1e-9, f"{name}: q_t - q_(t-1) vanishes at sampled roots (max {worst:.1e})")
    fine = 8 * n
    mins = []
    for q in chain[1:-1]:
        mins.append(min(q[0] + 2 * sum(q[j] * math.cos(j * math.pi * r / fine) for j in range(1, n))
                        for r in range(0, fine + 1, max(1, fine // 2000))))
    check(min(mins) > 0, f"{name}: float evaluation of q_1..q_(k-1) on a grid is positive (min {min(mins):.2e})")

print(f"\n{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
