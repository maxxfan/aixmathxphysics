#!/usr/bin/env python3
"""Standalone checker for the K_G(3) <= 1.4494177506 certificate (Python 3 standard library only).

usage:  python3 -I check_release.py [BUNDLE_DIR] [--cert FILE] [--no-verifier]

It shares no code with the project verifier (verifier/verify.py, verifier/hull3.py, verifier/accum.c) and uses
different algorithms: packed big-integer accumulation instead of C int64, gift wrapping instead of an incremental
hull, residual entries rounded up to 2^-80 instead of 2^-64. Every test is exact integer or rational arithmetic;
floating point is used only to PROPOSE hull candidates, never to decide anything.

What is checked (paper: Theorem 1, Lemma 1, Appendix A):
  0. the SHA-256 manifest of the bundle;
  1. the certificate's format: rational 3-vectors u_i (and w_j), rational scale lam > 0, integer weights n_k >= 1 with
     sum n_k <= D, sign masks, signed permutations, transposition flag;
  2. the local model L = 1/(G(1+tr)) sum_h [h.N (+ (h.N)^T)] / D with N = sum_k n_k s_k t_k^T, accumulated exactly;
  3. sum_k n_k/D + ||lam U V^T - L||_F <= 1, with every |E_ij| rounded UP to a multiple of 2^-80;
  4. for conv{+-u_i} (and conv{+-w_j}): a facet list F such that
       (a) every directed edge occurs once and its reverse once (closed oriented surface),
       (b) every facet (a,b,c) has det(a,b,c) > 0 (origin strictly inside its plane),
       (c) generic rays from the origin each meet exactly one facet (radial degree 1, so the facets cover all
           directions; then conv >= eta*Ball with eta = min distance of a facet plane from the origin),
       (d) every point lies on the closed inner side of every facet plane (sharpness: eta is the exact inradius),
       (e) Euler characteristic 2;
  5. the bound K_G(3) <= 1/(lam eta_A eta_B), rounded up to 10 digits, equals the claimed value, lies below the
     published 1.4546336567 (Designolle et al. 2023) and above the proven lower bound 1.4367003433.
  6. (optional) the frozen project verifier, run as a separate program, reports the same bound.
Exit code 0 iff everything passes.
"""
import hashlib
import json
import math
import os
import pathlib
import random
import subprocess
import sys
import time
from fractions import Fraction

CLAIM = "1.4494177506"
PUBLISHED = Fraction("1.4546336567")
PROVEN_LOWER = Fraction("1.4367003433")
LANE = 64                 # bits per lane in the packed accumulation (sum of weights <= 2^50)
PREC = 80                 # residual entries rounded up to multiples of 2^-PREC


class Fail(Exception):
    pass


def log(msg):
    print(msg, flush=True)


# ----------------------------------------------------------------------------- 0. manifest

def check_manifest(bundle):
    man = bundle / "MANIFEST.sha256"
    if not man.exists():
        raise Fail("MANIFEST.sha256 missing")
    n = 0
    for line in man.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        name = name.strip().lstrip("*")
        f = bundle / name
        if not f.is_file():
            raise Fail(f"manifest lists a missing file: {name}")
        if hashlib.sha256(f.read_bytes()).hexdigest() != digest:
            raise Fail(f"SHA-256 mismatch: {name}")
        n += 1
    return n


# ----------------------------------------------------------------------------- 1. parsing

def rat(x, what):
    if isinstance(x, bool) or isinstance(x, float):
        raise Fail(f"{what}: not an exact rational")
    if isinstance(x, int):
        return Fraction(x)
    if not isinstance(x, str) or len(x) > 1000:
        raise Fail(f"{what}: not an exact rational string")
    try:
        return Fraction(x)
    except (ValueError, ZeroDivisionError):
        raise Fail(f"{what}: unparsable rational")


def vectors(v, what):
    if not isinstance(v, list) or not v:
        raise Fail(f"{what}: need a non-empty list of 3-vectors")
    out = []
    for i, row in enumerate(v):
        if not isinstance(row, list) or len(row) != 3:
            raise Fail(f"{what}[{i}]: not a 3-vector")
        q = [rat(c, f"{what}[{i}]") for c in row]
        den = 1
        for c in q:
            den = den * c.denominator // math.gcd(den, c.denominator)
        out.append((int(q[0] * den), int(q[1] * den), int(q[2] * den), den))   # homogeneous (X, Y, Z, W), W > 0
    return out


def intval(x, what, lo=1):
    if isinstance(x, bool):
        raise Fail(f"{what}: not an integer")
    if isinstance(x, str) and x.isdigit():
        x = int(x)
    if not isinstance(x, int) or x < lo:
        raise Fail(f"{what}: must be an integer >= {lo}")
    return x


def mask(h, m, what):
    if not isinstance(h, str):
        raise Fail(f"{what}: sign mask must be a hex string")
    try:
        x = int(h, 16)
    except ValueError:
        raise Fail(f"{what}: not hexadecimal")
    if x < 0 or x >> m:
        raise Fail(f"{what}: bits beyond the number of settings")
    return x


def parse(cert):
    if cert.get("type") != "local_model":
        raise Fail("type must be local_model")
    U = vectors(cert.get("alice"), "alice")
    V = U if cert.get("bob") == "alice" else vectors(cert.get("bob"), "bob")
    mA, mB = len(U), len(V)
    lam = rat(cert.get("scale"), "scale")
    if lam <= 0:
        raise Fail("scale must be positive")
    D = intval(cert.get("weight_den"), "weight_den")
    strat = cert.get("strategies")
    if not isinstance(strat, list) or not strat:
        raise Fail("strategies: need a non-empty list")
    S = []
    for k, s in enumerate(strat):
        if not (isinstance(s, list) and len(s) == 3):
            raise Fail(f"strategies[{k}]: [weight, mask_A, mask_B]")
        S.append((intval(s[0], f"strategies[{k}].weight"), mask(s[1], mA, f"strategies[{k}].A"),
                  mask(s[2], mB, f"strategies[{k}].B")))
    wsum = sum(w for w, _, _ in S)
    if wsum > D:
        raise Fail("total weight exceeds 1")
    tr = cert.get("transpose", False)
    if not isinstance(tr, bool) or (tr and mA != mB):
        raise Fail("transpose flag invalid")
    H = []
    raw = cert.get("transforms")
    if raw is None:
        H = [(list(range(mA)), [1] * mA, list(range(mB)), [1] * mB)]
    else:
        for i, t in enumerate(raw):
            if not (isinstance(t, list) and len(t) == 4):
                raise Fail(f"transforms[{i}]: [perm_A, signs_A, perm_B, signs_B]")
            pa, sa, pb, sb = t
            for p, m in ((pa, mA), (pb, mB)):
                if not (isinstance(p, list) and len(p) == m and sorted(p) == list(range(m))):
                    raise Fail(f"transforms[{i}]: not a permutation")
            sa, sb = mask(sa, mA, f"transforms[{i}].sA"), mask(sb, mB, f"transforms[{i}].sB")
            H.append((pa, [1 if (sa >> x) & 1 else -1 for x in range(mA)],
                      pb, [1 if (sb >> x) & 1 else -1 for x in range(mB)]))
    if not H:
        raise Fail("transforms: empty list")
    return dict(U=U, V=V, mA=mA, mB=mB, lam=lam, D=D, S=S, wsum=wsum, tr=tr, H=H, same=cert.get("bob") == "alice")


# ----------------------------------------------------------------------------- 2. accumulation

def accumulate(c):
    """N[i][j] = sum_k n_k s_k(i) t_k(j), exactly, by packing Bob's index into 64-bit lanes of big integers."""
    mA, mB, S = c["mA"], c["mB"], c["S"]
    if c["wsum"] >= 1 << (LANE - 1):
        raise Fail("weights too large for the lane width")
    Q = []
    for (w, a, b) in S:
        wb = w.to_bytes(LANE // 8, "little")
        buf = bytearray(LANE // 8 * mB)
        x = b
        while x:
            low = x & -x
            j = low.bit_length() - 1
            buf[8 * j:8 * j + 8] = wb
            x ^= low
        Q.append(int.from_bytes(bytes(buf), "little"))
    Ct = sum(Q)                                   # lanes: sum of n_k over k with t_k(j) = +1
    rows = [[] for _ in range(mA)]
    wplus = [0] * mA
    for k, (w, a, b) in enumerate(S):
        x = a
        while x:
            low = x & -x
            i = low.bit_length() - 1
            rows[i].append(k)
            wplus[i] += w
            x ^= low
    nb = LANE // 8 * mB + 8

    def lanes(z):
        bts = z.to_bytes(nb, "little")
        return [int.from_bytes(bts[8 * j:8 * j + 8], "little") for j in range(mB)]
    ct = lanes(Ct)
    W = c["wsum"]
    N = []
    for i in range(mA):
        A = lanes(sum(Q[k] for k in rows[i])) if rows[i] else [0] * mB
        base = W - 2 * wplus[i]
        N.append([4 * A[j] - 2 * ct[j] + base for j in range(mB)])   # 4(++) - 2 W+_i - 2 Ct_j + W
    return N


def symmetrize(c, N):
    mA, mB = c["mA"], c["mB"]
    L = [[0] * mB for _ in range(mA)]
    for (pa, sa, pb, sb) in c["H"]:
        for i in range(mA):
            src = N[pa[i]]
            s = sa[i]
            row = L[i]
            L[i] = [x + s * sb[j] * src[pb[j]] for j, x in enumerate(row)]
    if c["tr"]:
        L = [[L[i][j] + L[j][i] for j in range(mB)] for i in range(mA)]
    return L, len(c["H"]) * (2 if c["tr"] else 1) * c["D"]          # model = L / Qden


def residual_ok(c, L, Qden):
    """Exact test sum n_k/D + ||lam U V^T - L/Qden||_F <= 1, with |E_ij| rounded up to multiples of 2^-PREC."""
    ln, ld = c["lam"].numerator, c["lam"].denominator
    total = 0
    for i, (X, Y, Z, al) in enumerate(c["U"]):
        row = L[i]
        a1 = ln * Qden
        a2 = ld * al
        for j, (X2, Y2, Z2, be) in enumerate(c["V"]):
            num = a1 * (X * X2 + Y * Y2 + Z * Z2) - a2 * be * row[j]
            if num:
                den = a2 * be * Qden
                ub = -((-abs(num) << PREC) // den)
                total += ub * ub
    slack = c["D"] - c["wsum"]
    ok = total * c["D"] ** 2 <= (slack * slack) << (2 * PREC)
    return ok, math.sqrt(total) / 2.0 ** PREC, slack / c["D"]


# ----------------------------------------------------------------------------- 4. hull (gift wrapping)

def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def plane(P, i, j, k):
    """Integer (n1, n2, n3, n4) with value(p) = n . (X, Y, Z, W)_p = det4[P_i; P_j; P_k; p]; value(origin) = det3."""
    rows = [P[i], P[j], P[k]]

    def minor(cols):
        r = [[row[c] for c in cols] for row in rows]
        return det3(r[0], r[1], r[2])
    return (-minor((1, 2, 3)), minor((0, 2, 3)), -minor((0, 1, 3)), minor((0, 1, 2)))


def val(pl, p):
    return pl[0] * p[0] + pl[1] * p[1] + pl[2] * p[2] + pl[3] * p[3]


def gift_wrap(P, k_near=36):
    n = len(P)
    fl = [(p[0] / p[3], p[1] / p[3], p[2] / p[3]) for p in P]      # floats only to propose candidates
    near = []
    for i in range(n):
        xi = fl[i]
        d = sorted(range(n), key=lambda j: (fl[j][0] - xi[0]) ** 2 + (fl[j][1] - xi[1]) ** 2 + (fl[j][2] - xi[2]) ** 2)
        near.append([j for j in d[1:k_near + 1]])

    def is_facet(i, j, k):
        pl = plane(P, i, j, k)
        return pl[3] > 0 and all(val(pl, p) >= 0 for p in P)
    first = None
    for j in near[0][:16]:
        for k in near[0][:16]:
            if j != k and det3(P[0][:3], P[j][:3], P[k][:3]) > 0 and is_facet(0, j, k):
                first = (0, j, k)
                break
        if first:
            break
    if first is None:
        raise Fail("hull: no initial facet found")
    facets, owner, todo = [], {}, []

    def add(f):
        i, j, k = f
        for e in ((i, j), (j, k), (k, i)):
            if e in owner:
                raise Fail("hull: directed edge used twice during wrapping")
            owner[e] = len(facets)
        facets.append(f)
        for e in ((i, j), (j, k), (k, i)):
            if (e[1], e[0]) not in owner:
                todo.append(e)
    add(first)
    while todo:
        a, b = todo.pop()
        if (b, a) in owner:
            continue
        cand = [x for x in dict.fromkeys(near[a] + near[b]) if x not in (a, b)]
        d = cand[0]
        for p in cand[1:]:
            if val(plane(P, b, a, d), P[p]) < 0:
                d = p
        add((b, a, d))
        if len(facets) > 4 * n:
            raise Fail("hull: wrapping does not terminate")
    return facets


def hull_check(P, F, rng):
    """Exact checks (a)-(e); returns eta^2 = min over facets of squared plane distance."""
    n = len(P)
    de = {}
    for (i, j, k) in F:
        for e in ((i, j), (j, k), (k, i)):
            de[e] = de.get(e, 0) + 1
    if any(v != 1 for v in de.values()) or any((b, a) not in de for (a, b) in de):
        raise Fail("hull (a): not a closed oriented surface")
    planes = [plane(P, i, j, k) for (i, j, k) in F]
    if any(pl[3] <= 0 for pl in planes):
        raise Fail("hull (b): a facet with det <= 0")
    rays = 0
    tries = 0
    while rays < 5:
        tries += 1
        if tries > 200:
            raise Fail("hull (c): no generic ray found")
        d = (rng.randint(-10 ** 15, 10 ** 15), rng.randint(-10 ** 15, 10 ** 15), rng.randint(-10 ** 15, 10 ** 15))
        hits, generic = 0, True
        for (i, j, k) in F:
            a, b, c = P[i][:3], P[j][:3], P[k][:3]                  # positive multiples of the points: same cones
            x, y, z = det3(d, b, c), det3(a, d, c), det3(a, b, d)   # barycentric numerators (det3(a,b,c) > 0)
            if min(x, y, z) >= 0:
                if 0 in (x, y, z):
                    generic = False
                    break
                hits += 1
        if not generic:
            continue
        if hits != 1:
            raise Fail(f"hull (c): a generic ray meets {hits} facets")
        rays += 1
    for pl in planes:                                               # (d) exact, all pairs
        for p in P:
            if val(pl, p) < 0:
                raise Fail("hull (d): a point outside a facet plane")
    if len({x for f in F for x in f}) - len(de) // 2 + len(F) != 2:
        raise Fail("hull (e): Euler characteristic is not 2")
    best = None
    for (i, j, k) in F:
        a = [Fraction(t, P[i][3]) for t in P[i][:3]]
        b = [Fraction(t, P[j][3]) for t in P[j][:3]]
        c = [Fraction(t, P[k][3]) for t in P[k][:3]]
        u = [b[t] - a[t] for t in range(3)]
        v = [c[t] - a[t] for t in range(3)]
        nrm = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        h = nrm[0] * a[0] + nrm[1] * a[1] + nrm[2] * a[2]
        r = h * h / (nrm[0] ** 2 + nrm[1] ** 2 + nrm[2] ** 2)
        best = r if best is None or r < best else best
    return best


def eta2(vecs, rng):
    pts, seen = [], set()
    for (X, Y, Z, W) in vecs + [(-X, -Y, -Z, W) for (X, Y, Z, W) in vecs]:
        g = math.gcd(math.gcd(math.gcd(abs(X), abs(Y)), abs(Z)), W)
        key = (X // g, Y // g, Z // g, W // g)
        if key not in seen:
            seen.add(key)
            pts.append(key)
    F = gift_wrap(pts)
    return hull_check(pts, F, rng), len(F)


# ----------------------------------------------------------------------------- driver

def up10(q):
    t = q * 10 ** 20
    r = math.isqrt(t.numerator // t.denominator)
    while Fraction(r * r) < t:
        r += 1
    s = str(r)
    return f"{s[:-10]}.{s[-10:]}"


def check(cert_path, seed=2026):
    """Run checks 1-5 on one certificate; return a dict, raise Fail on any failure."""
    t0 = time.time()
    c = parse(json.loads(pathlib.Path(cert_path).read_text()))
    N = accumulate(c)
    L, Qden = symmetrize(c, N)
    ok, resid, slack = residual_ok(c, L, Qden)
    if not ok:
        raise Fail(f"residual: weight {1 - slack:.12f} + ||E||_F {resid:.3e} > 1")
    rng = random.Random(seed)
    eA, fA = eta2(c["U"], rng)
    eB, fB = (eA, fA) if c["same"] else eta2(c["V"], rng)
    K2 = 1 / (c["lam"] ** 2 * eA * eB)
    return dict(K=up10(K2), K2=K2, eta2_A=eA, eta2_B=eB, facets=(fA, fB), m=(c["mA"], c["mB"]), atoms=len(c["S"]),
                G=len(c["H"]), weight=1 - slack, resid=resid, lam=c["lam"], seconds=round(time.time() - t0, 1))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    bundle = pathlib.Path(args[0] if args else pathlib.Path(__file__).resolve().parent).resolve()
    cert = bundle / "certificate" / "KG3_upper_1.4494177506_rad721.json"
    if "--cert" in sys.argv:
        cert = pathlib.Path(sys.argv[sys.argv.index("--cert") + 1]).resolve()
    fails = 0

    def step(name, fn):
        nonlocal fails
        try:
            log(f"PASS  {name}: {fn()}")
        except Fail as e:
            fails += 1
            log(f"FAIL  {name}: {e}")

    step("0 manifest", lambda: f"{check_manifest(bundle)} files match MANIFEST.sha256")
    res = {}

    def core():
        res.update(check(cert))
        return (f"m = {res['m']}, {res['atoms']} orbit atoms x {res['G']} transforms (+ transposition), "
                f"weight {res['weight']:.12f}, ||E||_F <= {res['resid']:.3e}, facets {res['facets']}, "
                f"{res['seconds']} s")
    step("1-4 format, exact accumulation, residual, hull (a)-(e)", core)
    if res:
        def bound():
            lam, e = res["lam"], res["eta2_A"]
            e15 = math.floor(e * 10 ** 15)
            msg = (f"lambda = {lam.numerator}/{lam.denominator}; eta_A^2 >= 0.{str(e15)[-15:]} (exact rational with "
                   f"{len(str(e.denominator))}-digit denominator); K_G(3) <= {res['K']}")
            if res["K"] != CLAIM:
                raise Fail(msg + f" != claimed {CLAIM}")
            if not res["K2"] < PUBLISHED ** 2:
                raise Fail(msg + " is not below the published 1.4546336567")
            if not res["K2"] > PROVEN_LOWER ** 2:
                raise Fail(msg + " contradicts the proven lower bound 1.4367003433")
            return msg + f" (claimed {CLAIM}; published 1.4546336567; proven lower bound 1.4367003433)"
        step("5 bound", bound)
    if "--no-verifier" not in sys.argv and (bundle / "verifier" / "verify.py").exists():
        def verifier():
            r = subprocess.run([sys.executable, "-I", str(bundle / "verifier" / "verify.py"), str(cert), "--json",
                                "--no-lock"], capture_output=True, text=True, timeout=1800)
            try:
                rep = json.loads(r.stdout)
            except json.JSONDecodeError:
                raise Fail(f"verifier did not run: {r.stderr.strip()[-300:]}")
            if rep.get("kg3_upper_str") != CLAIM or rep.get("verdict") in ("INVALID", "VERIFIER_ALARM"):
                raise Fail(f"verifier says {rep.get('verdict')} {rep.get('kg3_upper_str')}: {rep.get('reason')}")
            return f"{rep['verifier']}: {rep['verdict']}, K_G(3) <= {rep['kg3_upper_str']} ({rep['seconds']} s)"
        step("6 project verifier (separate code)", verifier)
    log("ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
