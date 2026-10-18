#!/usr/bin/env python3
"""Tests for the p04 verifier (stdlib only).  usage: test_verifier.py [--trials N]

1. hull3 against a brute-force facet enumeration on random rational point sets (incl. degenerate ones)
2. hull3.validate rejects open, doubled, flipped and non-convex surfaces
3. exact shrinking factors: octahedron 1/3, cube 1, and the published algebraic values for the
   BellPolytopes.jl 16- and 46-vertex polyhedra (DIB+23 supplemental) after exact rationalization
4. accum.c against a pure-Python accumulation (strategies, signed permutations, transpose)
5. every verdict path of verify_candidate, and INVALID on malformed input
"""
import argparse
import itertools
import json
import math
import pathlib
import random
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import hull3  # noqa: E402
import verify  # noqa: E402

ROOT = HERE.parent.parent
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg)


def hom(p):
    """Fraction triple -> homogeneous integer point."""
    den = 1
    for c in p:
        den = den * c.denominator // math.gcd(den, c.denominator)
    return (int(p[0] * den), int(p[1] * den), int(p[2] * den), den)


def rational_unit(v, maxden=10 ** 15):
    """Exact rational unit vector near the float unit vector v (inverse stereographic projection)."""
    x, y, z = v
    flip = z < 0
    if flip:
        z = -z
    a = Fraction(x / (1 + z)).limit_denominator(maxden)
    b = Fraction(y / (1 + z)).limit_denominator(maxden)
    s = 1 + a * a + b * b
    p = (2 * a / s, 2 * b / s, (1 - a * a - b * b) / s)
    return (p[0], p[1], -p[2]) if flip else p


# ----------------------------------------------------------------------------- 1. hull vs brute force

def test_hull_random(trials):
    rng = random.Random(1)
    for t in range(trials):
        n = rng.randint(5, 14)
        kind = t % 3
        pts = []
        while len(pts) < n:
            if kind == 0:
                p = (rng.randint(-6, 6), rng.randint(-6, 6), rng.randint(-6, 6), rng.randint(1, 4))
            elif kind == 1:   # coarse grid: many coplanar / collinear points
                p = (rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-2, 2), 1)
            else:             # centrally symmetric sets, as the verifier builds them
                q = (rng.randint(-5, 5), rng.randint(-5, 5), rng.randint(-5, 5), rng.randint(1, 3))
                pts.append((-q[0], -q[1], -q[2], q[3]))
                p = q
            pts.append(p)
        try:
            F = hull3.hull(pts, seed=t)
        except hull3.HullError as e:
            check(any(w in str(e) for w in ("coplanar", "collinear", "coincide")), f"hull trial {t}: unexpected {e}")
            continue
        Pf = hull3.to_fracs(pts)
        planes = [hull3._plane(Pf, i, j, k) for (i, j, k) in F]
        got = hull3.facet_planes_normalized(planes)
        want = hull3.brute_facet_planes(pts)
        check(got == want, f"hull trial {t}: facet planes differ from brute force")
        if kind == 2:  # centrally symmetric and full-dimensional: the full validation must pass
            try:
                hull3.validate(pts, F)
            except hull3.HullError as e:
                check(False, f"hull trial {t}: validate rejected a correct hull: {e}")


# ----------------------------------------------------------------------------- 2. validate rejects

def test_validate_rejects():
    cube = [(x, y, z, 1) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    F = hull3.hull(cube)
    hull3.validate(cube, F)
    for name, bad in [("open", F[1:]),
                      ("doubled", F + F),
                      ("flipped", [(i, k, j) for (i, j, k) in F])]:
        try:
            hull3.validate(cube, bad)
            check(False, f"validate accepted a {name} surface")
        except hull3.HullError:
            pass
    # a single facet with flipped orientation (det < 0) breaks both closedness and positive orientation
    bad = [(F[0][0], F[0][2], F[0][1])] + F[1:]
    try:
        hull3.validate(cube, bad)
        check(False, "validate accepted a surface with one flipped facet")
    except hull3.HullError:
        pass
    # a point outside a facet plane by an exact margin far below any floating-point resolution (v1.1: exact tests)
    tiny = [(10 ** 40 + 1, 0, 0, 10 ** 40)]
    try:
        hull3.validate(cube + tiny, F)
        check(False, "validate accepted a point outside a facet by 1e-40")
    except hull3.HullError:
        pass
    # a closed surface that is not the hull: octahedron facets on points where one is pushed outward
    octa = [(1, 0, 0, 1), (-1, 0, 0, 1), (0, 1, 0, 1), (0, -1, 0, 1), (0, 0, 1, 1), (0, 0, -1, 1)]
    Fo = hull3.hull(octa)
    pushed = list(octa) + [(2, 2, 2, 1)]
    try:
        hull3.validate(pushed, Fo)
        check(False, "validate accepted a surface with a point outside")
    except hull3.HullError:
        pass


# ----------------------------------------------------------------------------- 3. shrinking factors

def test_eta():
    e2, _ = verify.eta2_of([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    check(e2 == Fraction(1, 3), f"octahedron eta^2 = {e2}")
    e2, _ = verify.eta2_of([(1, 1, 1, 1), (1, 1, -1, 1), (1, -1, 1, 1), (-1, 1, 1, 1)])
    check(e2 == 1, f"cube eta^2 = {e2}")
    seeds = json.loads((ROOT / "seeds" / "bellpolytopes_polyhedra.json").read_text())["polyhedra"]
    # published exact values (DIB+23 supplemental, Appendix A): largest roots of
    #   1 - 40X + 530X^2 - 2480X^3 + 2245X^4  (m = 16)  and  81 - 3024X + 16794X^2 - 29640X^3 + 16045X^4  (m = 46)
    for m, want in (("16", 0.8511948093773854), ("46", 0.9440979216955705)):
        vecs = [hom(rational_unit(v)) for v in seeds[m]]
        e2, nf = verify.eta2_of(vecs)
        check(abs(float(e2) - want) < 1e-11, f"pol{m}: eta^2 {float(e2)!r} vs published {want!r}")


# ----------------------------------------------------------------------------- 4. accumulation

def py_accum(mA, mB, strategies, transforms, tr):
    N0 = [[0] * mB for _ in range(mA)]
    for (W, a, b) in strategies:
        for i in range(mA):
            si = 1 if (a >> i) & 1 else -1
            for j in range(mB):
                N0[i][j] += W * si * (1 if (b >> j) & 1 else -1)
    S0 = [[0] * mB for _ in range(mA)]
    for (pa, sa, pb, sb) in transforms:
        for i in range(mA):
            for j in range(mB):
                sg = (1 if (sa >> i) & 1 else -1) * (1 if (sb >> j) & 1 else -1)
                S0[i][j] += sg * N0[pa[i]][pb[j]]
    if tr:
        S0 = [[S0[i][j] + S0[j][i] for j in range(mB)] for i in range(mA)]
    return N0, S0


def test_accum(trials):
    rng = random.Random(2)
    for t in range(trials):
        mA = rng.randint(1, 70)
        tr = t % 2 == 1
        mB = mA if tr else rng.randint(1, 70)
        K = rng.randint(1, 30)
        strategies = [(rng.randint(1, 2 ** 40), rng.getrandbits(mA), rng.getrandbits(mB)) for _ in range(K)]
        G = rng.randint(1, 5)
        transforms = []
        for _ in range(G):
            pa = list(range(mA))
            pb = list(range(mB))
            rng.shuffle(pa)
            rng.shuffle(pb)
            transforms.append((pa, rng.getrandbits(mA), pb, rng.getrandbits(mB)))
        N0, S = verify.run_accum(mA, mB, strategies, transforms, int(tr), 60)
        N0p, Sp = py_accum(mA, mB, strategies, transforms, tr)
        check(list(N0) == [x for r in N0p for x in r], f"accum trial {t}: N0 differs")
        check(list(S) == [x for r in Sp for x in r], f"accum trial {t}: S differs")
        verify.spot_check(mA, mB, strategies, transforms, tr, N0, S, seed=t)


# ----------------------------------------------------------------------------- 5. verdicts

def octa_candidate(**kw):
    c = {"type": "local_model", "target": "KG3:upper", "alice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
         "bob": "alice", "scale": "1", "weight_den": 8,
         "strategies": [[1, format(s, "x"), format(s, "x")] for s in range(8)]}
    c.update(kw)
    return c


def run(c):
    return verify.verify_candidate(c, lock=False)


def test_verdicts():
    r = run(octa_candidate())
    check(r["verdict"] == "BELOW_RECORD" and r["kg3_upper_str"] == "3.0000000000", f"octahedron: {r.get('summary')}")
    # the same model via transforms: one strategy (+,+,+) for both, symmetrized over the 8 sign flips
    flips = [[[0, 1, 2], format(s, "x"), [0, 1, 2], format(s, "x")] for s in range(8)]
    r = run(octa_candidate(weight_den=1, strategies=[[1, "7", "7"]], transforms=flips))
    check(r["verdict"] == "BELOW_RECORD" and r["kg3_upper_str"] == "3.0000000000", f"transforms: {r.get('summary')}")
    # transpose: strategies (s, t) and their transposes average to the identity as well
    r = run(octa_candidate(transpose=True))
    check(r["verdict"] == "BELOW_RECORD", f"transpose: {r.get('summary')}")
    # scale below 1 leaves slack: residual 0.01*I has norm 0.01*sqrt(3) < slack 0.02 with weight 0.98
    r = run(octa_candidate(scale="99/100", weight_den=800, strategies=[[98, format(s, "x"), format(s, "x")] for s in range(8)]))
    check(r["verdict"] == "BELOW_RECORD", f"slack: {r.get('summary')} {r.get('reason')}")
    # residual too big: target 1.01*I with full weight
    r = run(octa_candidate(scale="101/100"))
    check(r["verdict"] == "INVALID" and "not certified local" in r["reason"], f"residual: {r.get('reason')}")
    # weight above one
    r = run(octa_candidate(weight_den=7))
    check(r["verdict"] == "INVALID" and "exceeds 1" in r["reason"], f"weight: {r.get('reason')}")
    # a non-local target that would give a fake bound: CHSH-like target with full weight on 1 strategy
    r = run(octa_candidate(weight_den=1, strategies=[[1, "7", "7"]]))
    check(r["verdict"] == "INVALID", "single strategy cannot give the identity")
    # malformed inputs
    for name, c in [("float coord", octa_candidate(alice=[[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])),
                    ("2d vectors", octa_candidate(alice=[[1, 0], [0, 1]])),
                    ("coplanar", octa_candidate(alice=[[1, 0, 0], [0, 1, 0], [1, 1, 0]])),
                    ("bad hex", octa_candidate(strategies=[[8, "xyz", "1"]])),
                    ("bits beyond m", octa_candidate(strategies=[[8, "f", "1"]])),
                    ("zero weight", octa_candidate(strategies=[[0, "1", "1"]])),
                    ("bad perm", octa_candidate(transforms=[[[0, 0, 2], "7", [0, 1, 2], "7"]])),
                    ("negative scale", octa_candidate(scale="-1")),
                    ("bad type", octa_candidate(type="bell")),
                    ("transpose mismatch", octa_candidate(bob=[[1, 0, 0], [0, 1, 0]], transpose=True)),
                    ("not a dict", [1, 2, 3])]:
        r = run(c)
        check(r["verdict"] == "INVALID", f"{name}: expected INVALID, got {r['verdict']}")
    # verdict paths that need other records: patch the record reader
    orig = verify.records
    try:
        verify.records = lambda: (Fraction(7, 2), Fraction(1, 1), Fraction(1, 10 ** 6))
        check(run(octa_candidate())["verdict"] == "IMPROVEMENT", "IMPROVEMENT path")
        verify.records = lambda: (Fraction(3), Fraction(1), Fraction(1, 10 ** 6))
        check(run(octa_candidate())["verdict"] == "MATCHES_RECORD", "MATCHES_RECORD path")
        verify.records = lambda: (Fraction(4), Fraction(31, 10), Fraction(1, 10 ** 6))
        check(run(octa_candidate())["verdict"] == "VERIFIER_ALARM", "VERIFIER_ALARM path")
    finally:
        verify.records = orig
    # records file sanity
    up, low, d = verify.records()
    check(up == Fraction("1.4494177506") and low == Fraction("1.4367003433") and d == Fraction(1, 10 ** 6),
          "records.json values")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=40)
    a = ap.parse_args()
    test_hull_random(a.trials)
    test_validate_rejects()
    test_eta()
    test_accum(a.trials)
    test_verdicts()
    print("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} FAILURES")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
