#!/usr/bin/env python3
"""Verifier for p03: T-count of GF(2^n) multiplication (Table 8 item 3).

Candidate JSON:
  {"n": 7, "modulus": "11110001",              # irreducible P, coefficients x^n .. x^0
   "factors": ["0101...", ...],                # length 3n; char i = variable i
                                               # (a_0..a_{n-1}, b_0..b_{n-1}, c_0..c_{n-1})
   "author": "...", "construction": "...", "script": "..."}

The multiplier |a,b,c> -> |a,b,c+ab> is H_c D H_c with D = (-1)^{c.(a*b mod P)}.  Each factor u
is one T gate (a phase gadget exp(i pi/4 (u.x mod 2))); the T-count is the number of factors.
The T-count does not depend on P or on the basis (those changes are CNOT circuits).

Checks (both exact):
  1. tensor   sum_u u^{x3} = S_P over GF(2): every variable and every pair of variables lies in
              an even number of factors, and every triple lies in an odd number iff it is a
              monomial a_i b_j c_k of c.(a*b mod P)            (Heyfron-Campbell criterion)
  2. circuit  (n <= 8) phasecheck.c, written independently: 4f(x) - sum_u (u.x) mod 8 is a
              diagonal Clifford phase for every input x               (direct certificate)
Disagreement between 1 and 2, or a T-count below a proven lower bound, is VERIFIER_ALARM.

Verdicts: IMPROVEMENT (T-count < record), MATCHES_RECORD, BELOW_RECORD (worse than record),
NO_TABLE_ENTRY, INVALID, VERIFIER_ALARM.  usage: verify.py CANDIDATE.json [--json]
"""
import argparse
import hashlib
import itertools
import json
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RECORDS = ROOT / "data" / "records.json"
PHASE_SRC, PHASE_BIN = HERE / "phasecheck.c", HERE / "phasecheck"
PINNED = ["verify.py", "phasecheck.c"]
MAX_PHASE_N = 8


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()[:16]


def source_hashes():
    h = {f: sha(HERE / f) for f in PINNED}
    h["records.json"] = sha(RECORDS)
    if PHASE_BIN.exists():
        h["phasecheck"] = sha(PHASE_BIN)
    return h


def ensure_phasecheck():
    if not PHASE_BIN.exists() or PHASE_BIN.stat().st_mtime < PHASE_SRC.stat().st_mtime:
        subprocess.run(["cc", "-O3", "-o", str(PHASE_BIN), str(PHASE_SRC)], check=True)


def check_pins():
    pins = HERE / "PINS.json"
    if not pins.exists():
        return "unpinned"
    ensure_phasecheck()
    want, cur = json.loads(pins.read_text()), source_hashes()
    bad = [k for k in sorted(set(want) | set(cur)) if want.get(k) != cur.get(k)]
    return "ok" if not bad else "MODIFIED: " + ",".join(bad)


class Invalid(Exception):
    pass


# ------------------------------------------------------------------ GF(2)[x]

def pmod(a, p):
    dp = p.bit_length() - 1
    while a and a.bit_length() - 1 >= dp:
        a ^= p << (a.bit_length() - 1 - dp)
    return a


def pmulmod(a, b, p):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> (p.bit_length() - 1) & 1:
            a ^= p
    return r


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def irreducible(p):
    """Rabin's test: P of degree n is irreducible iff x^(2^n) = x mod P and
    gcd(x^(2^(n/q)) - x, P) = 1 for every prime q dividing n.  Polynomial time."""
    n = p.bit_length() - 1
    if n < 1:
        return False
    if n == 1:
        return True
    if not p & 1:
        return False
    def frob(k):  # x^(2^k) mod P
        y = 2
        for _ in range(k):
            y = pmulmod(y, y, p)
        return y
    if frob(n) != 2:
        return False
    primes = [q for q in range(2, n + 1) if n % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    return all(pgcd(p, frob(n // q) ^ 2) == 1 for q in primes)


def target(n, p):
    """Odd entries of the signature tensor, keyed by support: only triples (a_i, b_j, c_k)."""
    out = set()
    for i in range(n):
        for j in range(n):
            v = pmod(1 << (i + j), p)
            for k in range(n):
                if v >> k & 1:
                    out.add((i, n + j, 2 * n + k))
    return out


def tensor_rows(factors):
    """rows[(i, j)] (i <= j) = XOR of the factors whose support contains i and j.
    Bit k of rows[(i, j)] is the entry (i, j, k) of sum u^{x3} mod 2, so these rows determine the
    whole symmetric tensor; cost sum_u w(w+1)/2 big-int XORs (O(w^2) instead of O(w^3))."""
    rows = {}
    for u in factors:
        s = [i for i in range(u.bit_length()) if u >> i & 1]
        for x in range(len(s)):
            for y in range(x, len(s)):
                key = (s[x], s[y])
                rows[key] = rows.get(key, 0) ^ u
    return rows


def target_rows(n, p):
    rows = {}
    for (i, j, k) in target(n, p):
        for (x, y, z) in ((i, j, k), (i, k, j), (j, k, i)):
            rows[(x, y)] = rows.get((x, y), 0) | (1 << z)
    return rows


def tensor_mismatch(factors, n, p):
    """None if sum u^{x3} == S_P, else up to three mismatching sorted index triples."""
    got, want = tensor_rows(factors), target_rows(n, p)
    bad = []
    for key in sorted(set(got) | set(want)):
        diff = got.get(key, 0) ^ want.get(key, 0)
        while diff and len(bad) < 3:
            k = diff.bit_length() - 1
            bad.append(sorted({key[0], key[1], k}))
            diff ^= 1 << k
        if len(bad) >= 3:
            break
    return bad or None


def cube_sum(factors):
    """Reference implementation (enumerates triples; O(w^3)); used by the tests as a cross-check."""
    odd = set()
    for u in factors:
        s = [i for i in range(u.bit_length()) if u >> i & 1]
        for r in (1, 2, 3):
            for t in itertools.combinations(s, r):
                odd.symmetric_difference_update((t,))
    return odd


# ------------------------------------------------------------------ verify

def load(cand):
    n = cand.get("n") if isinstance(cand, dict) else None
    if isinstance(n, bool) or not isinstance(n, int):
        raise Invalid("n must be an integer")
    if not 2 <= n <= 64:
        raise Invalid("n must be in 2..64")
    mod = str(cand.get("modulus", "")).strip()
    if len(mod) != n + 1 or set(mod) - {"0", "1"} or mod[0] != "1":
        raise Invalid("modulus must be a 0/1 string of length n+1 with leading 1 (x^n .. x^0)")
    p = int(mod, 2)
    if not irreducible(p):
        raise Invalid(f"modulus {mod} is not irreducible over GF(2)")
    fs = cand.get("factors")
    if not isinstance(fs, list) or not fs:
        raise Invalid("factors must be a non-empty list")
    if len(fs) > 20000:
        raise Invalid("more than 20000 factors")
    N = 3 * n
    factors = []
    for i, s in enumerate(fs):
        if not isinstance(s, str) or len(s) != N or set(s) - {"0", "1"}:
            raise Invalid(f"factor {i} must be a 0/1 string of length 3n = {N}")
        u = int(s[::-1], 2)
        if u == 0:
            raise Invalid(f"factor {i} is zero")
        factors.append(u)
    return n, p, mod, factors


def phase_check(n, p, factors, timeout):
    ensure_phasecheck()
    inp = f"{n} {p} {len(factors)}\n" + " ".join(map(str, factors)) + "\n"
    r = subprocess.run([str(PHASE_BIN)], input=inp, capture_output=True, text=True, timeout=timeout)
    return json.loads(r.stdout.strip().splitlines()[-1])


def verify_candidate(cand, exact=False, quick=False, threads=8, timeout=1800, budget=0, seed=1,
                     lock=True, records_path=None, **kw):
    """records_path is for tests only; the board always uses data/records.json."""
    t0 = time.time()
    rep = {"verifier": "p03-tcount/2", "pins": check_pins(), "sources": source_hashes()}
    try:
        n, p, mod, factors = load(cand)
    except Invalid as e:
        rep.update(verdict="INVALID", reason=str(e), seconds=round(time.time() - t0, 3))
        return rep
    m = len(factors)
    rep.update(n=n, modulus=mod, tcount=m, target=f"gf2^{n}",
               hash=hashlib.sha256(f"{n}:{p}:".encode() + b",".join(
                   str(u).encode() for u in sorted(factors))).hexdigest()[:16])
    if len(set(factors)) < m:
        rep["note"] = "duplicate factors (each pair cancels in the tensor and only wastes T gates)"

    bad = tensor_mismatch(factors, n, p)
    tensor_ok = bad is None
    if not tensor_ok:
        rep["tensor_mismatch"] = bad
    rep["tensor_check"] = tensor_ok

    phase = None
    if n <= MAX_PHASE_N and not quick:
        phase = phase_check(n, p, factors, timeout)
        rep["circuit_check"] = phase
        if phase.get("ok") != tensor_ok:
            rep.update(verdict="VERIFIER_ALARM", seconds=round(time.time() - t0, 3),
                       reason="tensor check and independent circuit check disagree")
            rep["summary"] = f"VERIFIER_ALARM: checks disagree on GF(2^{n}) T-count {m}"
            return rep
    else:
        rep["circuit_check"] = {"ok": None, "reason": f"skipped (n > {MAX_PHASE_N} or --quick)"}

    if not tensor_ok:
        rep.update(verdict="INVALID", reason="factors do not reproduce the multiplier's signature tensor",
                   seconds=round(time.time() - t0, 3))
        rep["summary"] = f"INVALID: GF(2^{n}) candidate with {m} factors fails the tensor check"
        return rep

    recs = json.loads(pathlib.Path(records_path or RECORDS).read_text())
    rec = recs["records"].get(str(n))
    low = recs.get("proven_lower", {}).get(str(n))
    rep["record"] = rec
    if low and m < low["tcount"]:
        rep.update(verdict="VERIFIER_ALARM",
                   reason=f"T-count {m} is below the proven lower bound {low['tcount']}")
    elif rec is None:
        rep["verdict"] = "NO_TABLE_ENTRY"
    elif m < rec["tcount"]:
        rep["verdict"] = "IMPROVEMENT"
    elif m == rec["tcount"]:
        rep["verdict"] = "MATCHES_RECORD"
    else:
        rep["verdict"] = "BELOW_RECORD"
    recs_txt = f"record {rec['tcount']} ({rec['by']})" if rec else "no record"
    checks = "tensor+circuit" if phase else "tensor"
    rep["summary"] = f"{rep['verdict']}: GF(2^{n}) multiplier T-count {m} vs {recs_txt}; checks: {checks}"
    rep["seconds"] = round(time.time() - t0, 3)
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidate")
    ap.add_argument("--quick", action="store_true", help="tensor check only")
    ap.add_argument("--exact", action="store_true", help="(accepted for interface compatibility)")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    rep = verify_candidate(json.loads(pathlib.Path(a.candidate).read_text()), quick=a.quick)
    if a.json:
        print(json.dumps(rep, indent=1))
    else:
        print(rep.get("summary") or f"{rep['verdict']}: {rep.get('reason', '')}")
        for k in ("reason", "note", "tensor_mismatch", "circuit_check", "pins", "seconds"):
            if k in rep:
                print(f"{k}: {rep[k]}")
    sys.exit(0 if rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") else 1)


if __name__ == "__main__":
    main()
