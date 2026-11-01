#!/usr/bin/env python3
"""p05 verifier tests (stdlib only).

  * data: every catalog inequality is a facet with local bound b0; data/*.json rebuild from sources
  * known records reproduce: exact strategies at the Pal-Vertesi values, SOS certificates at the
    known Tsirelson values (fixtures generated with qlib/certs.py)
  * every verdict path (records injected where no real improvement exists)
  * invalid input is rejected and never crashes
  * independent cross-checks on random instances: Bell values against an explicit Kronecker
    computation written here; the SOS identity against a string-rewriting implementation;
    the PSD shift against exact eigenvalue constructions (Cayley-rational orthogonal matrices)

run: python3 verifier/tests/test_verifier.py [--trials 30]
"""
import argparse
import copy
import itertools
import json
import pathlib
import random
import subprocess
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))
import qexact as qx  # noqa: E402
import verify as V  # noqa: E402

FIX = HERE / "fixtures"
FAIL = []


def check(name, cond, info=""):
    print(("ok   " if cond else "FAIL ") + name + (f"  [{info}]" if info and not cond else ""))
    if not cond:
        FAIL.append(name)


def fx(name):
    return json.loads((FIX / name).read_text())


def run(cand):
    return V.verify_candidate(cand, lock=False)


# ----------------------------------------------------------------------------- data

def test_data():
    cat = json.loads(V.CATALOG.read_text())["inequalities"]
    bad = []
    for name, c in cat.items():
        mA, mB = c["mA"], c["mB"]
        best, sat = None, []
        for a in itertools.product((0, 1), repeat=mA):
            for b in itertools.product((0, 1), repeat=mB):
                v = (sum(c["bA"][x] * a[x] for x in range(mA)) + sum(c["bB"][y] * b[y] for y in range(mB))
                     + sum(c["bAB"][x][y] * a[x] * b[y] for x in range(mA) for y in range(mB)))
                p = list(a) + list(b) + [a[x] * b[y] for x in range(mA) for y in range(mB)] + [1]
                if best is None or v > best:
                    best, sat = v, [p]
                elif v == best:
                    sat.append(p)
        if best != c["b0"] or rank(sat) != mA + mB + mA * mB:
            bad.append(name)
    check(f"catalog: {len(cat)} inequalities are facets with local bound b0", not bad, bad[:5])
    r = subprocess.run([sys.executable, str(ROOT / "data" / "make_records.py"), "--check"],
                       capture_output=True, text=True)
    check("data/*.json rebuild exactly from data/sources (make_records.py --check)", r.returncode == 0, r.stdout + r.stderr)
    recs = json.loads(V.RECORDS.read_text())["records"]
    opened = sorted(n for n, x in recs.items() if x["status"] == "open")
    check("open targets: A14 A21 A47 A62 A64 A68 A82 A84 A89 (A-list, PV 2010 / Naceur order 4) and J30 J48 J62 J73 (J97 closed)",
          opened == ["A14", "A21", "A47", "A62", "A64", "A68", "A82", "A84", "A89", "J30", "J48", "J62", "J73"]
          and recs["J97"]["status"] == "closed", opened)
    check("A-list: 84 entries from data/sources/naceur2025_alist; A14 = PV 2010 lower 0.4759513, upper Naceur order 4 0.4768133",
          sum(1 for k in json.loads(V.CATALOG.read_text())["inequalities"] if k[0] == "A" and k[1:].isdigit()) == 84
          and recs["A14"]["lower"] == "0.4759513" and recs["A14"]["upper"] == "0.4768133")
    check("overlays: PV2010 values kept (J73 1.1223170, J30 0.4363842); phase-1 bests are the bar; I4422_20, J34 closed numerically",
          recs["J73"]["pv2010"]["lower"] == "1.1223170" and recs["J30"]["lower"] == "0.4363842"
          and recs["J73"]["lower"] == "1.1223870369" and recs["J97"]["upper"] == "1.1583625755"
          and recs["I4422_20"]["status"] == recs["J34"]["status"] == "closed_numerically")
    bad_win = [n for n, x in recs.items() if Fraction(x["upper"]) < Fraction(x["lower"])]
    check("every record window is consistent (lower <= upper)", not bad_win, bad_win)


def rank(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    r = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


# ----------------------------------------------------------------------------- known records

def test_fixtures():
    cases = [("CHSH_lower.json", "MATCHES_RECORD", 0.2071067), ("J97_lower.json", "MATCHES_RECORD", 1.1583625),
             ("J73_lower.json", "BELOW_RECORD", 1.1204591), ("I4422_20_lower.json", "MATCHES_RECORD", 0.4676793)]
    for f, want, approx in cases:
        rep = run(fx(f))
        check(f"strategy {f}: {want} at {approx}", rep["verdict"] == want and abs(rep.get("value_float", 0) - approx) < 2e-7,
              rep.get("summary"))
    for f, want, lo, hi in [("CHSH_upper_1.json", "MATCHES_RECORD", 0.20710678118, 0.2071068),
                            ("J9_upper_2.json", "MATCHES_RECORD", 0.5, 0.5000001),
                            ("J73_upper_2.json", "BELOW_RECORD", 1.15345, 1.15346)]:
        rep = run(fx(f))
        v = rep.get("value_float", -1)
        check(f"sos {f}: {want}, bound in [{lo}, {hi}]", rep["verdict"] == want and lo <= v <= hi, rep.get("summary"))
        check(f"sos {f}: exact PD cross-check ran", rep.get("exact_pd_crosscheck") is True, rep.get("exact_pd_crosscheck"))
    # Tsirelson: the certified CHSH bound is >= (sqrt2-1)/2 exactly
    rep = run(fx("CHSH_upper_1.json"))
    b = Fraction(rep["bound_exact"])
    check("CHSH certified bound >= (sqrt(2)-1)/2 (exact: (2b+1)^2 >= 2)", (2 * b + 1) ** 2 >= 2)
    s = run(fx("CHSH_lower.json"))
    v = Fraction(s["violation_exact"])
    check("CHSH exact strategy value <= (sqrt(2)-1)/2 (exact: (2v+1)^2 <= 2)", (2 * v + 1) ** 2 <= 2)


# ----------------------------------------------------------------------------- verdict paths

class Inject:
    """Temporarily replace the record of one inequality."""

    def __init__(self, name, **rec):
        self.name, self.rec = name, rec

    def __enter__(self):
        self.orig = V.load_problem

        def fake(n):
            ineq, rec = self.orig(n)
            if n == self.name:
                rec = dict(rec or {}, **self.rec)
            return ineq, rec
        V.load_problem = fake

    def __exit__(self, *a):
        V.load_problem = self.orig


def test_paths():
    s = fx("J97_lower.json")
    with Inject("J97", lower="1.1580000", upper="1.1590000", status="open"):
        rep = run(s)
    check("lower IMPROVEMENT path (injected record)", rep["verdict"] == "IMPROVEMENT", rep.get("summary"))
    check("IMPROVEMENT report has target, gain, summary", rep.get("target") == "J97:lower" and rep.get("gain", 0) > 0
          and "IMPROVEMENT" in rep.get("summary", ""))
    with Inject("J97", lower="1.1580000", upper="1.1582000", status="open"):
        rep = run(s)
    check("lower VERIFIER_ALARM when value exceeds the published upper bound", rep["verdict"] == "VERIFIER_ALARM", rep.get("summary"))
    with Inject("J97", lower="1.1580000", upper="1.1590000", status="matched"):
        rep = run(s)
    check("no IMPROVEMENT on a non-open record", rep["verdict"] != "IMPROVEMENT", rep.get("summary"))
    with Inject("J97", lower="1.1600000", upper="1.1700000", status="open"):
        rep = run(s)
    check("lower BELOW_RECORD path", rep["verdict"] == "BELOW_RECORD", rep.get("summary"))
    u = fx("J73_upper_2.json")
    with Inject("J73", lower="1.1204592", upper="1.1600000", status="open"):
        rep = run(u)
    check("upper IMPROVEMENT path (injected record)", rep["verdict"] == "IMPROVEMENT" and rep["target"] == "J73:upper",
          rep.get("summary"))
    with Inject("J73", lower="1.1600000", upper="1.1700000", status="open"):
        rep = run(u)
    check("upper VERIFIER_ALARM when bound is below a published strategy", rep["verdict"] == "VERIFIER_ALARM", rep.get("summary"))
    rep = run(u)
    check("upper BELOW_RECORD path (real record)", rep["verdict"] == "BELOW_RECORD", rep.get("summary"))
    # a weak random strategy is BELOW_RECORD
    w = copy.deepcopy(s)
    w["psi"] = ["1", "0", "0", "1"]
    rep = run(w)
    check("weak strategy BELOW_RECORD", rep["verdict"] == "BELOW_RECORD", rep.get("summary"))
    check("report carries pins and hash", "pins" in rep and rep.get("hash"))
    big = json.loads((FIX / "J73_d10.json").read_text())
    rep = run(big)
    check("d=10 strategy uses the Bob-side re-check and verifies", rep.get("recheck") == "bob-side"
          and rep["verdict"] in ("BELOW_RECORD", "MATCHES_RECORD"), rep.get("summary"))
    check("identical content gives identical hash; metadata ignored",
          run(dict(s, author="x", script="y"))["hash"] == run(s)["hash"])


# ----------------------------------------------------------------------------- invalid input

def test_invalid():
    s, u = fx("J97_lower.json"), fx("J73_upper_2.json")
    bad = {
        "not a dict": [1, 2],
        "no type": {"inequality": "J97"},
        "unknown inequality": dict(s, inequality="J999"),
        "float entries": dict(s, psi=[0.5, 0.0, 0.0, 0.5]),
        "wrong psi length": dict(s, psi=["1", "0"]),
        "dim too large": dict(s, dA=99),
        "too many Alice projectors": dict(s, A=s["A"] + [[]]),
        "dependent spanning vectors": dict(s, A=[[["1", "0"], ["2", "0"]]] + s["A"][1:]),
        "spanning vector wrong length": dict(s, A=[[["1", "0", "0"]]] + s["A"][1:]),
        "complex entry in real field": dict(s, psi=[["1", "1"], "0", "0", "1"]),
        "zero state": dict(s, psi=["0", "0", "0", "0"]),
        "bad field": dict(s, field="quaternion"),
        "garbage number": dict(s, psi=["1/0x", "0", "0", "1"]),
        "huge exponent": dict(s, psi=["1e999999", "0", "0", "1"]),
        "sos: bad word": dict(u, words=["", "A9"] + u["words"][2:]),
        "sos: duplicate word": dict(u, words=["", "A1", "A1A1"] + u["words"][3:]),
        "sos: wrong row length": dict(u, Y={"denominator": u["Y"]["denominator"], "upper": u["Y"]["upper"][:-1]}),
        "sos: non-integer entry": dict(u, Y={"denominator": u["Y"]["denominator"],
                                             "upper": [["1/2"] + r[1:] if i == 0 else r for i, r in enumerate(u["Y"]["upper"])]}),
        "sos: wrong lambda": dict(u, **{"lambda": "1"}),
        "sos: no Y": {k: v for k, v in u.items() if k != "Y"},
        "sos: negative denominator": dict(u, Y={"denominator": "-5", "upper": u["Y"]["upper"]}),
    }
    for name, cand in bad.items():
        try:
            rep = run(cand)
            ok = rep["verdict"] == "INVALID"
        except Exception as e:  # must never raise
            ok, rep = False, {"reason": repr(e)}
        check(f"INVALID: {name}", ok, rep.get("verdict", "") + " " + str(rep.get("reason", "")))
    # tampering: change one off-diagonal numerator -> identity fails
    t = copy.deepcopy(u)
    t["Y"]["upper"][1][3] = str(int(t["Y"]["upper"][1][3]) + 1)
    rep = run(t)
    check("tampered SOS entry -> INVALID (identity fails)", rep["verdict"] == "INVALID" and "identity" in rep.get("reason", ""),
          rep.get("reason"))
    # an indefinite Y that satisfies the identity is not accepted cheaply: penalty makes it useless
    t = copy.deepcopy(u)
    d = int(t["Y"]["denominator"])
    # add -0.5 to Y_11 and +0.5 to Y_00? Y_00 is lambda; changing Y_11 alone breaks identity, so shift
    # both entries of class 'A1' consistently: Y_{0,A1} += 1/4 (twice) and Y_{A1,A1} -= 1/2
    iA1 = u["words"].index("A1")
    t["Y"]["upper"][0][iA1] = str(int(t["Y"]["upper"][0][iA1]) + d // 4)
    t["Y"]["upper"][iA1][0] = str(int(t["Y"]["upper"][iA1][0]) - d // 2)
    rep = run(t)
    check("identity-preserving indefinite perturbation is caught by the PSD shift (bound gets worse)",
          rep["verdict"] in ("BELOW_RECORD", "INVALID") and rep.get("value_float", 99) > 1.16, rep.get("summary"))


# ----------------------------------------------------------------------------- independent checks

def kron_value(ineq, A, B, psi, dA, dB, field):
    """Independent: explicit full matrices with complex Fractions as pairs."""
    def cm(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def ca(a, b):
        return (a[0] + b[0], a[1] + b[1])

    def to_pair(z):
        return (z.re, z.im) if field == "complex" else (z, Fraction(0))
    n = dA * dB
    Wm = [[(Fraction(0), Fraction(0))] * n for _ in range(n)]
    Ia = [[Fraction(int(i == j)) for j in range(dA)] for i in range(dA)]
    Ib = [[Fraction(int(i == j)) for j in range(dB)] for i in range(dB)]
    terms = []
    for x in range(ineq["mA"]):
        terms.append((ineq["bA"][x], A[x], None))
        for y in range(ineq["mB"]):
            terms.append((ineq["bAB"][x][y], A[x], B[y]))
    for y in range(ineq["mB"]):
        terms.append((ineq["bB"][y], None, B[y]))
    for c, P, Q in terms:
        if not c:
            continue
        for a1 in range(dA):
            for a2 in range(dA):
                pa = to_pair(P[a1][a2]) if P is not None else (Ia[a1][a2], Fraction(0))
                if pa == (0, 0):
                    continue
                for b1 in range(dB):
                    for b2 in range(dB):
                        qb = to_pair(Q[b1][b2]) if Q is not None else (Ib[b1][b2], Fraction(0))
                        if qb == (0, 0):
                            continue
                        pr = cm(pa, qb)
                        r, s_ = a1 * dB + b1, a2 * dB + b2
                        Wm[r][s_] = ca(Wm[r][s_], (c * pr[0], c * pr[1]))
    ps = [to_pair(z) for z in psi]
    num = (Fraction(0), Fraction(0))
    den = Fraction(0)
    for r in range(n):
        cr = (ps[r][0], -ps[r][1])
        den += ps[r][0] ** 2 + ps[r][1] ** 2
        for s_ in range(n):
            num = ca(num, cm(cr, cm(Wm[r][s_], ps[s_])))
    assert num[1] == 0
    return num[0] / den


def test_random_values(trials, rng):
    cat = json.loads(V.CATALOG.read_text())["inequalities"]
    names = ["CHSH", "I3322", "J73", "I4422_20", "J30"]
    agree = 0
    for t in range(trials):
        name = rng.choice(names)
        ineq = cat[name]
        field = rng.choice(["real", "complex"])
        F = qx.Field(field)
        dA, dB = rng.randint(1, 3), rng.randint(1, 3)

        def rnum():
            return Fraction(rng.randint(-9, 9), rng.randint(1, 5))

        def rent():
            return rnum() if field == "real" else qx.GQ(rnum(), rnum())

        def rproj(d):
            r = rng.randint(0, d)
            while True:
                vecs = [[rent() for _ in range(d)] for _ in range(r)]
                try:
                    return qx.projector_from_span(F, vecs, d)
                except qx.Invalid:
                    continue
        A = [rproj(dA) for _ in range(ineq["mA"])]
        B = [rproj(dB) for _ in range(ineq["mB"])]
        psi = [rent() for _ in range(dA * dB)]
        if all(F.is_zero(z) for z in psi):
            psi[0] = F.one
        ok = all(qx.is_projector(F, P) for P in A + B)
        v1 = qx.strategy_value(F, ineq, A, B, psi, dA, dB)
        v2 = qx.strategy_value_kron(F, ineq, A, B, psi, dA, dB)
        v3 = kron_value(ineq, A, B, psi, dA, dB, field)
        v4 = qx.strategy_value_bob(F, ineq, A, B, psi, dA, dB)
        agree += ok and v1 == v2 == v3 == v4
    check(f"random strategies: exact projectors and four Bell-value computations agree ({agree}/{trials})", agree == trials)


def word_to_str(m):
    return "".join(f"a{x}" for x in m[0]) + "".join(f"b{y}" for y in m[1])


def rewrite_product(u, v):
    """Independent normal form of u^dag v via strings: letters 'aX'/'bY'; commute a's left; drop repeats."""
    toks = [("a", x) for x in reversed(u[0])] + [("b", y) for y in reversed(u[1])] + \
           [("a", x) for x in v[0]] + [("b", y) for y in v[1]]
    a = [k for p, k in toks if p == "a"]
    b = [k for p, k in toks if p == "b"]
    changed = True
    while changed:
        changed = False
        for seq in (a, b):
            for i in range(len(seq) - 1):
                if seq[i] == seq[i + 1]:
                    del seq[i]
                    changed = True
                    break
    return word_to_str((tuple(a), tuple(b)))


def test_random_sos_identity(trials, rng):
    agree = 0
    for t in range(trials):
        mA, mB = rng.randint(2, 4), rng.randint(2, 4)
        allw = set()
        while len(allw) < rng.randint(3, 12):
            na, nb = rng.randint(0, 3), rng.randint(0, 2)
            w = "".join(f"A{rng.randint(1, mA)}" for _ in range(na)) + "".join(f"B{rng.randint(1, mB)}" for _ in range(nb))
            allw.add(qx.parse_word(w, mA, mB))
        words = sorted(allw)
        n = len(words)
        Y = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                Y[i][j] = Y[j][i] = rng.randint(-5, 5)
        got = {word_to_str(m): c for m, c in qx.sos_identity(words, Y).items() if c}
        want = {}
        for i in range(n):
            for j in range(n):
                if Y[i][j]:
                    k = rewrite_product(words[i], words[j])
                    want[k] = want.get(k, 0) + Y[i][j]
        want = {k: v for k, v in want.items() if v}
        agree += got == want
    check(f"random SOS identities match an independent rewriting implementation ({agree}/{trials})", agree == trials)


def cayley(n, rng):
    """Exact rational orthogonal matrix Q = (I - S)(I + S)^{-1}, S skew-symmetric."""
    S = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = Fraction(rng.randint(-3, 3), rng.randint(1, 4))
            S[i][j], S[j][i] = v, -v
    F = qx.Field("real")
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    ImS = [[I[i][j] - S[i][j] for j in range(n)] for i in range(n)]
    IpS = [[I[i][j] + S[i][j] for j in range(n)] for i in range(n)]
    return qx.mat_mul(F, ImS, qx.mat_inv(F, IpS))


def test_psd(trials, rng):
    F = qx.Field("real")
    sound, tight = 0, 0
    for t in range(trials):
        n = rng.randint(2, 9)
        Q = cayley(n, rng)
        lam = [Fraction(rng.randint(1, 50), 10) for _ in range(n)]
        neg = rng.random() < 0.5
        if neg:
            lam[0] = -Fraction(rng.randint(1, 1000), 10 ** rng.randint(3, 9))
        else:
            lam[0] = Fraction(rng.randint(0, 3), 10 ** rng.randint(6, 12))
        D = [[lam[i] if i == j else Fraction(0) for j in range(n)] for i in range(n)]
        M = qx.mat_mul(F, qx.mat_mul(F, Q, D), qx.transpose(Q))
        den = 1
        for row in M:
            for e in row:
                den = den * e.denominator // __import__("math").gcd(den, e.denominator)
        Ynum = [[int(e * den) for e in row] for row in M]
        best, _ = V.psd_shift(Ynum, den)
        tau = best[0]
        lmin = min(lam)
        sound += tau >= -lmin  # tau must cover the most negative eigenvalue
        tight += tau <= max(0, -lmin) * Fraction(101, 100) + Fraction(1, 10 ** 8)
    check(f"PSD shift is sound on exact-eigenvalue matrices ({sound}/{trials})", sound == trials)
    check(f"PSD shift is tight within 1% + 1e-8 ({tight}/{trials})", tight == trials)
    # exact PD test sanity
    check("exact Bareiss PD test: [[2,1],[1,2]] yes, [[1,2],[2,1]] no, [[0]] no",
          qx.ldl_pd([[2, 1], [1, 2]]) and not qx.ldl_pd([[1, 2], [2, 1]]) and not qx.ldl_pd([[0]]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=30)
    a = ap.parse_args()
    rng = random.Random(20261011)
    test_data()
    test_fixtures()
    test_paths()
    test_invalid()
    test_random_values(a.trials, rng)
    test_random_sos_identity(a.trials * 3, rng)
    test_psd(a.trials, rng)
    print(f"\n{len(FAIL)} failure(s)" if FAIL else "\nall tests passed")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
