#!/usr/bin/env python3
"""Verifier for p04: certified upper bounds on the Grothendieck constant K_G(3) (Table 8 item 4).

Candidate type `local_model` (format in PROBLEM.md). It certifies

    K_G(3) <= 1 / (lam * eta_A * eta_B)

from four exact objects:
  * rational vectors u_1..u_mA (Alice) and w_1..w_mB (Bob) in Q^3, and a rational scale lam > 0;
    the target correlation matrix is T_ij = lam * u_i . w_j;
  * a local model: positive integer weights W_k over a denominator D, deterministic strategies
    s_k in {+-1}^mA, t_k in {+-1}^mB, an optional list of G signed index permutations h, and an
    optional transpose flag, defining
        L = 1/(G (1+tr)) sum_h [ h.N0 (+ its transpose) ] / D,   N0 = sum_k W_k s_k t_k^T,
        (h.X)_ij = sgA_h(i) sgB_h(j) X[piA_h(i), piB_h(j)]
    (every term is W_k/D times a deterministic strategy, so L has total weight sum_k W_k / D);
  * the residual E = T - L, which must satisfy   sum_k W_k / D + ||E||_F <= 1.
    The Frobenius unit ball lies in the local polytope of full correlators (E_{s,t}<M, s t^T>^2 =
    ||M||_F^2 for uniform s, t, so max <M, s t^T> >= ||M||_F for every M), hence T is local, with
    zero marginals after pairing every strategy with its negation;
  * eta_A^2 = squared inradius of conv{+-u_i}, eta_B^2 likewise, computed by an exact 3D hull
    (hull3.py) that is then validated independently: closed oriented surface, every facet positively oriented
    (origin strictly inside its plane), a generic ray meets exactly one facet, so the radial projection has degree 1
    and covers every direction (this alone certifies conv{+-u_i} >= eta*Ball, see hull3.py), and, for sharpness, every
    point inside every facet plane (v1.1: all pairs tested exactly, no floating-point screen).
For unit vectors a, b: eta_A a is a convex combination of the +-u_i and eta_B b of the +-w_j, so the
model gives correlations lam*eta_A*eta_B (a.b) with zero marginals for ALL projective measurements;
hence lam*eta_A*eta_B * sum M_ij a_i.b_j <= L(M) for every Bell matrix M and K_G(3) <= 1/(lam eta_A eta_B)
(equivalently the two-qubit Werner state is local for projective measurements at visibility
lam*eta_A*eta_B; Acin-Gisin-Toner 2006).

Exactness: vectors, scale and weights are exact rationals; N0 and the symmetrization are exact int64
sums in accum.c (inputs capped so |entries| <= 2^62), spot-checked entrywise in Python; the residual
norm is bounded above with every |E_ij| rounded up to a multiple of 2^-64; the final comparison is
exact rational arithmetic on squares.  No floating-point value enters any test (v1.1; v1.0 used a float screen
to skip point/facet tests for points far inside a facet, which never changed a verdict).

Verdicts (records in data/records.json, delta = 1e-6):
  IMPROVEMENT     certified bound <= record_upper - delta
  MATCHES_RECORD  within delta of the record
  BELOW_RECORD    a valid but weaker bound
  VERIFIER_ALARM  certified bound below the proven lower bound 1.4367003433 (contradiction)
  INVALID         malformed input, or the decomposition / hull does not certify the claim

usage: verify.py CANDIDATE.json [--json]
"""
import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import os
import pathlib
import random
import struct
import subprocess
import sys
import tempfile
import time
from array import array
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import hull3  # noqa: E402

RECORDS = ROOT / "data" / "records.json"
ACC_SRC, ACC_BIN = HERE / "accum.c", HERE / "accum"
LOCKFILE = HERE / ".heavy.lock"
PINNED = ["verify.py", "hull3.py", "accum.c"]
TARGET = "KG3:upper"

MAX_M = 2000              # vectors per party
MAX_K = 5_000_000         # strategies
MAX_G = 2048              # symmetrizing transforms
MAX_D = 2 ** 50           # weight denominator
MAX_DIGITS = 400          # digits in any numerator / denominator
PREC = 64                 # residual entries rounded up to multiples of 2^-PREC


class Invalid(Exception):
    pass


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()[:16]


def source_hashes():
    h = {f: sha(HERE / f) for f in PINNED}
    h["records.json"] = sha(RECORDS)
    if ACC_BIN.exists():
        h["accum"] = sha(ACC_BIN)
    return h


def check_pins():
    pins = HERE / "PINS.json"
    if not pins.exists():
        return "unpinned"
    want = json.loads(pins.read_text())
    cur = source_hashes()
    bad = [f for f in set(want) | set(cur) if want.get(f) != cur.get(f)]
    return "ok" if not bad else "MODIFIED: " + ",".join(sorted(bad))


def ensure_accum():
    if not ACC_BIN.exists() or ACC_BIN.stat().st_mtime < ACC_SRC.stat().st_mtime:
        subprocess.run(["cc", "-O2", "-o", str(ACC_BIN), str(ACC_SRC)], check=True)


@contextlib.contextmanager
def heavy_lock(enabled=True):
    if not enabled:
        yield
        return
    with open(LOCKFILE, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def canonical_hash(cand):
    keep = {k: v for k, v in cand.items() if k not in ("author", "construction", "script", "note", "notes",
                                                         "claimed", "created")}
    return hashlib.sha256(json.dumps(keep, sort_keys=True).encode()).hexdigest()[:16]


def records():
    r = json.loads(RECORDS.read_text())
    return Fraction(r["upper"]["value"]), Fraction(r["lower"]["value"]), Fraction(r["meta"]["delta"])


# ----------------------------------------------------------------------------- parsing

def parse_int(x, what, lo=None, hi=None):
    if isinstance(x, bool):
        raise Invalid(f"{what}: booleans are not integers")
    if isinstance(x, str):
        if not x.strip().lstrip("-").isdigit() or len(x) > MAX_DIGITS:
            raise Invalid(f"{what}: not an integer string")
        x = int(x)
    if not isinstance(x, int):
        raise Invalid(f"{what}: must be an integer (or a string of digits)")
    if (lo is not None and x < lo) or (hi is not None and x > hi):
        raise Invalid(f"{what}: {x} outside [{lo}, {hi}]")
    return x


def parse_rat(x, what):
    if isinstance(x, bool):
        raise Invalid(f"{what}: booleans are not numbers")
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        if len(x) > 2 * MAX_DIGITS + 8:
            raise Invalid(f"{what}: too long")
        try:
            q = Fraction(x.strip())
        except (ValueError, ZeroDivisionError):
            raise Invalid(f"{what}: not a rational string")
        if len(str(abs(q.numerator))) > MAX_DIGITS or len(str(q.denominator)) > MAX_DIGITS:
            raise Invalid(f"{what}: more than {MAX_DIGITS} digits")
        return q
    raise Invalid(f"{what}: give an integer or an exact rational string like \"3/7\" (no JSON floats)")


def parse_vectors(v, what):
    if not isinstance(v, list) or not 1 <= len(v) <= MAX_M:
        raise Invalid(f"{what}: need a list of 1..{MAX_M} vectors")
    out = []
    for i, row in enumerate(v):
        if not isinstance(row, list) or len(row) != 3:
            raise Invalid(f"{what}[{i}]: each vector has 3 rational coordinates")
        q = [parse_rat(c, f"{what}[{i}]") for c in row]
        den = 1
        for c in q:
            den = den * c.denominator // math.gcd(den, c.denominator)
        out.append((int(q[0] * den), int(q[1] * den), int(q[2] * den), den))
    return out


def parse_signs(h, m, what):
    if not isinstance(h, str) or not 1 <= len(h) <= m // 4 + 3:
        raise Invalid(f"{what}: sign vector must be a hex string of at most {m // 4 + 1} digits")
    try:
        x = int(h, 16)
    except ValueError:
        raise Invalid(f"{what}: not hexadecimal")
    if x < 0 or x >> m:
        raise Invalid(f"{what}: has bits beyond index {m - 1}")
    return x


def parse_transform(t, mA, mB, i):
    if not (isinstance(t, list) and len(t) == 4):
        raise Invalid(f"transforms[{i}]: [permA, signsA_hex, permB, signsB_hex]")
    pa, sa, pb, sb = t
    for p, m, nm in ((pa, mA, "permA"), (pb, mB, "permB")):
        if not (isinstance(p, list) and len(p) == m and all(isinstance(x, int) and not isinstance(x, bool)
                                                              for x in p) and sorted(p) == list(range(m))):
            raise Invalid(f"transforms[{i}].{nm}: must be a permutation of 0..{m - 1}")
    return pa, parse_signs(sa, mA, f"transforms[{i}].signsA"), pb, parse_signs(sb, mB, f"transforms[{i}].signsB")


# ----------------------------------------------------------------------------- pieces

def eta2_of(vecs):
    """Exact squared inradius of conv{+-v}: hull, independent validation, min facet distance^2."""
    pts = list(vecs) + [(-X, -Y, -Z, W) for (X, Y, Z, W) in vecs]
    # drop exact duplicates (same rational point) so the hull sees each point once
    seen, uniq = set(), []
    for (X, Y, Z, W) in pts:
        g = math.gcd(math.gcd(math.gcd(abs(X), abs(Y)), abs(Z)), W)
        key = (X // g, Y // g, Z // g, W // g)
        if key not in seen:
            seen.add(key)
            uniq.append(key)
    try:
        F = hull3.hull(uniq)
        planes = hull3.validate(uniq, F)
    except hull3.HullError as e:
        raise Invalid(f"convex hull of +-vectors: {e}")
    return hull3.inradius2(planes), len(F)


def run_accum(mA, mB, strategies, transforms, tr, timeout):
    ensure_accum()
    wa, wb = (mA + 63) // 64, (mB + 63) // 64
    with tempfile.TemporaryDirectory(prefix="p04acc") as td:
        fin, fout = pathlib.Path(td) / "in.bin", pathlib.Path(td) / "out.bin"
        with open(fin, "wb") as f:
            f.write(struct.pack("<5q", mA, mB, len(strategies), len(transforms), tr))
            for (W, a, b) in strategies:
                f.write(struct.pack("<q", W))
                f.write(a.to_bytes(8 * wa, "little"))
                f.write(b.to_bytes(8 * wb, "little"))
            for (pa, sa, pb, sb) in transforms:
                f.write(array("i", pa).tobytes())
                f.write(bytes((1 if (sa >> i) & 1 else 255) for i in range(mA)))
                f.write(array("i", pb).tobytes())
                f.write(bytes((1 if (sb >> j) & 1 else 255) for j in range(mB)))
        if sys.byteorder != "little":
            raise Invalid("verifier requires a little-endian machine")
        r = subprocess.run([str(ACC_BIN), str(fin), str(fout)], capture_output=True, timeout=timeout)
        if r.returncode != 0:
            raise Invalid(f"accumulation failed: {r.stderr.decode()[:200]}")
        raw = fout.read_bytes()
    n = mA * mB
    if len(raw) != 16 * n:
        raise Invalid("accumulation output has the wrong size")
    N0, S = array("q"), array("q")
    N0.frombytes(raw[:8 * n])
    S.frombytes(raw[8 * n:])
    return N0, S


def spot_check(mA, mB, strategies, transforms, tr, N0, S, seed):
    """Recompute random entries of N0 (from the strategies) and of S (from N0) in pure Python."""
    rng = random.Random(seed)
    K = len(strategies)
    nspots = max(4, min(24, 3_000_000 // max(K, 1)))
    for _ in range(nspots):
        i, j = rng.randrange(mA), rng.randrange(mB)
        v = 0
        for (W, a, b) in strategies:
            v += W if ((a >> i) & 1) == ((b >> j) & 1) else -W
        if v != N0[i * mB + j]:
            raise Invalid(f"internal: accumulation spot check failed at N0[{i}][{j}] (please report)")
    for _ in range(24):
        i, j = rng.randrange(mA), rng.randrange(mB)

        def sym(i, j):
            v = 0
            for (pa, sa, pb, sb) in transforms:
                x = N0[pa[i] * mB + pb[j]]
                v += x if ((sa >> i) & 1) == ((sb >> j) & 1) else -x
            return v
        want = sym(i, j) + (sym(j, i) if tr else 0)
        if want != S[i * mB + j]:
            raise Invalid(f"internal: symmetrization spot check failed at S[{i}][{j}] (please report)")


def residual_bound2(A, B, lam, S, Q, mB):
    """Return (sum ub_ij^2) where ub_ij = ceil(|E_ij| 2^PREC), E = lam*Gram(A,B) - S/Q, exactly."""
    ln, ld = lam.numerator, lam.denominator
    total = 0
    for i, (X, Y, Z, al) in enumerate(A):
        row = i * mB
        lnQ = ln * Q
        ldal = ld * al
        for j, (U, V, Wc, be) in enumerate(B):
            num = lnQ * (X * U + Y * V + Z * Wc) - ldal * be * S[row + j]
            if num:
                den = ldal * be * Q
                ub = -((-abs(num) << PREC) // den)
                total += ub * ub
    return total


def dec_up_sqrt(q, digits=10):
    """Decimal string >= sqrt(q) (q a positive Fraction)."""
    t = q * 10 ** (2 * digits)
    r = math.isqrt(t.numerator // t.denominator)
    while Fraction(r * r) < t:
        r += 1
    s = str(r).rjust(digits + 1, "0")
    return f"{s[:-digits]}.{s[-digits:]}"


def dec_down_sqrt(q, digits=10):
    t = q * 10 ** (2 * digits)
    r = math.isqrt(t.numerator // t.denominator)
    s = str(r).rjust(digits + 1, "0")
    return f"{s[:-digits]}.{s[-digits:]}"


# ----------------------------------------------------------------------------- main check

def check_local_model(cand, timeout=1800, seed=1):
    t0 = time.time()
    if cand.get("type", "local_model") != "local_model":
        raise Invalid("type must be 'local_model'")
    if cand.get("target", TARGET) != TARGET:
        raise Invalid(f"target must be {TARGET!r}")
    A = parse_vectors(cand.get("alice"), "alice")
    same = cand.get("bob") == "alice"
    B = A if same else parse_vectors(cand.get("bob"), "bob")
    mA, mB = len(A), len(B)
    lam = parse_rat(cand.get("scale"), "scale")
    if lam <= 0:
        raise Invalid("scale must be positive")
    D = parse_int(cand.get("weight_den"), "weight_den", 1, MAX_D)
    raw = cand.get("strategies")
    if not isinstance(raw, list) or not 1 <= len(raw) <= MAX_K:
        raise Invalid(f"strategies: need a list of 1..{MAX_K} [W, signsA_hex, signsB_hex]")
    strategies, sumW = [], 0
    for k, s in enumerate(raw):
        if not (isinstance(s, list) and len(s) == 3):
            raise Invalid(f"strategies[{k}]: [W, signsA_hex, signsB_hex]")
        W = parse_int(s[0], f"strategies[{k}].W", 1, MAX_D)
        strategies.append((W, parse_signs(s[1], mA, f"strategies[{k}].A"), parse_signs(s[2], mB, f"strategies[{k}].B")))
        sumW += W
    if sumW > D:
        raise Invalid(f"total weight {sumW}/{D} exceeds 1")
    tr = cand.get("transpose", False)
    if not isinstance(tr, bool):
        raise Invalid("transpose must be true or false")
    if tr and mA != mB:
        raise Invalid("transpose needs mA == mB")
    rt = cand.get("transforms")
    if rt is None:
        full = (1 << mA) - 1, (1 << mB) - 1
        transforms = [(list(range(mA)), full[0], list(range(mB)), full[1])]
    else:
        if not isinstance(rt, list) or not 1 <= len(rt) <= MAX_G:
            raise Invalid(f"transforms: need a list of 1..{MAX_G} signed permutations (or omit it)")
        transforms = [parse_transform(t, mA, mB, i) for i, t in enumerate(rt)]
    G = len(transforms)
    if G * (2 if tr else 1) > 4096:
        raise Invalid("too many transforms")

    etaA2, fA = eta2_of(A)
    etaB2, fB = (etaA2, fA) if same else eta2_of(B)
    t_hull = time.time() - t0

    N0, S = run_accum(mA, mB, strategies, transforms, int(tr), timeout)
    spot_check(mA, mB, strategies, transforms, tr, N0, S, seed)
    Q = D * G * (2 if tr else 1)
    sub2 = residual_bound2(A, B, lam, S, Q, mB)
    slack = D - sumW
    certified = sub2 * D * D <= (slack * slack) << (2 * PREC)
    resid_ub = math.sqrt(sub2) / 2.0 ** PREC
    weight = sumW / D
    info = dict(mA=mA, mB=mB, K=len(strategies), G=G, transpose=tr, weight=weight, resid_ub=resid_ub,
                margin=1.0 - weight - resid_ub, facets_A=fA, facets_B=fB,
                eta_A=math.sqrt(float(etaA2)), eta_B=math.sqrt(float(etaB2)), scale=float(lam),
                seconds=round(time.time() - t0, 1), hull_seconds=round(t_hull, 1))
    if not certified:
        raise Invalid(f"not certified local: weight {weight:.12f} + residual norm {resid_ub:.3e} > 1 "
                      f"(scale the target down, e.g. lam -> lam/(1+||E||))")

    v2 = lam * lam * etaA2 * etaB2      # squared certified Werner visibility
    U2 = 1 / v2                          # squared certified upper bound on K_G(3)
    rec_up, rec_low, delta = records()
    U_str, v_str = dec_up_sqrt(U2), dec_down_sqrt(v2)
    rec_s, low_s = f"{float(rec_up):.10f}", f"{float(rec_low):.10f}"
    out = dict(target=TARGET, kg3_upper=float(U_str), kg3_upper_str=U_str, werner_v_str=v_str, **info)
    if U2 < rec_low * rec_low:
        out.update(verdict="VERIFIER_ALARM",
                   reason=f"certified K_G(3) <= {U_str} is below the proven lower bound {low_s}",
                   summary=f"VERIFIER_ALARM: K_G(3) <= {U_str} contradicts K_G(3) >= {low_s}")
    elif U2 <= (rec_up - delta) ** 2:
        out.update(verdict="IMPROVEMENT",
                   summary=f"IMPROVEMENT: K_G(3) <= {U_str} (record {rec_s}); Werner v_c >= {v_str} "
                           f"[mA={mA}, mB={mB}, K={len(strategies)}, G={G}]")
    elif U2 <= (rec_up + delta) ** 2:
        out.update(verdict="MATCHES_RECORD", summary=f"MATCHES_RECORD: K_G(3) <= {U_str} (record {rec_s})")
    else:
        out.update(verdict="BELOW_RECORD",
                   summary=f"BELOW_RECORD: K_G(3) <= {U_str} (record {rec_s}); Werner v_c >= {v_str} "
                           f"[mA={mA}, mB={mB}, eta_A={info['eta_A']:.6f}, eta_B={info['eta_B']:.6f}, "
                           f"scale={float(lam):.6f}]")
    return out


def check_candidate(cand, timeout=1800, seed=1, **_):
    if not isinstance(cand, dict):
        raise Invalid("candidate must be a JSON object")
    return check_local_model(cand, timeout=timeout, seed=seed)


def verify_candidate(cand, exact=False, quick=False, threads=8, timeout=1800, budget=0, seed=1,
                     lock=True, **kw):
    rep = {"verifier": "p04-grothendieck-kg3/1.1", "pins": check_pins(), "sources": source_hashes()}
    try:
        rep["hash"] = canonical_hash(cand) if isinstance(cand, dict) else None
        with heavy_lock(lock):
            rep.update(check_candidate(cand, timeout=timeout, seed=seed))
    except Invalid as e:
        rep.update(verdict="INVALID", reason=str(e), target=TARGET, summary=f"INVALID: {e}")
    except subprocess.TimeoutExpired:
        rep.update(verdict="UNDECIDED", reason="timeout", target=TARGET, summary="UNDECIDED: timeout")
    except Exception as e:  # never let bad input crash the board
        rep.update(verdict="INVALID", reason=f"{type(e).__name__}: {e}", target=TARGET,
                   summary=f"INVALID: {type(e).__name__}: {e}")
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--no-lock", action="store_true")
    a = ap.parse_args()
    try:
        cand = json.loads(pathlib.Path(a.candidate).read_text())
    except (OSError, json.JSONDecodeError) as e:
        print(json.dumps({"verdict": "INVALID", "reason": str(e)}))
        sys.exit(1)
    rep = verify_candidate(cand, lock=not a.no_lock)
    print(json.dumps(rep, indent=1) if a.json else rep.get("summary", rep["verdict"]))
    sys.exit(0 if rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") else 1)


if __name__ == "__main__":
    main()
