#!/usr/bin/env python3
"""Verifier for p09: exact quantum ordered search (Table 8 item 9, query complexity of a specific
function).  Standard library only.

Theorem (Farhi-Goldstone-Gutmann-Sipser 1999; Childs-Landahl-Parrilo 2007; Thm 8 of
arXiv:2503.21090): an exact k-query translation-invariant quantum algorithm for n-element ordered
search exists iff there are symmetric Laurent polynomials q_0, ..., q_k of degree < n, nonnegative
on the unit circle, with

    q_0 = F_n (Fejer kernel, coefficients 1 - |j|/n),   q_k = 1,   constant term of every q_t = 1,
    q_t(z) = q_{t-1}(z) whenever z^n = (-1)^t           (t = 1..k).

Running a k-query n-element algorithm recursively searches N elements with k*ceil(log_n N) queries,
so the asymptotic constant is k / log2(n).  Explicit finite algorithms (snapshot 2026-10-29):
k = 5, n = 7265 (0.389810 log2 N, arXiv:2503.21090, LP + numerical nonnegativity check) and
k = 6, n = 90000 (0.364572 log2 N, Wu et al., ICML 2026, "strong numerical evidence", not
certified; they certify N*(6) < 94000).  Carolan-Childs (arXiv:2609.34293) prove the asymptotic
(1/pi) ln n + o(log n) without explicit finite instances.  The record is the best explicit
finite claim, (6, 90000); `beats_prior_explicit` reports n^5 > 7265^k.

Candidate (JSON):
    {"type": "fggs_chain", "k": K, "n": N,
     "den": D,      q_t coefficient j is q[t-1][j] / D   (t = 1..K-1, j = 0..N-1; D divisible by N)
     "q": [[...N ints...], ...K-1 lists...],
     "h_den": H,    h_t(z) = sum_j h[t-1][j] / H * z^j   (at most N ints each)
     "h": [[...], ...K-1 lists...]}

with q_t(theta) = q_t[0] + 2 sum_{j>=1} q_t[j] cos(j theta).  Everything is checked in exact
integer arithmetic:
  1. q_0 = F_n and q_K = 1 are fixed; every q_t has constant term 1 (q[t-1][0] = D);
  2. the root conditions: q_t - q_{t-1} vanishes on z^n = s (s = (-1)^t) iff, for j = 1..n-1,
         d_j + s d_{n-j} = 0,   d = coefficients of q_t - q_{t-1}
     (on those roots z^{-j} = s z^{n-j}, and 1, z, ..., z^{n-1} are independent there);
  3. nonnegativity of q_1..q_{K-1}: the residual r_t = q_t - |h_t|^2 (|h|^2 is the
     autocorrelation of h, computed exactly by Kronecker substitution) must be diagonally
     dominant, r_t[0] >= 2 sum_{j>=1} |r_t[j]|; then q_t(theta) >= |h_t(e^{i theta})|^2 >= 0.

Verdicts:
  IMPROVEMENT     certified and n^6 > 90000^k       (k / log2 n below 0.364572...)
  MATCHES_RECORD  certified and n^6 = 90000^k       (a rigorous certificate of the numerical record)
  BELOW_RECORD    certified and n^6 < 90000^k
  INVALID         malformed, or a condition fails
  VERIFIER_ALARM  certified beyond a proven or published bound: k < ln(n)/pi - 2 (Hoyer-Neerbek-
                  Shi exact lower bound, weakened), or n above the published feasibility limits
                  (k=2: 6, k=4: 605 [CLP07], k=5: 7265 [arXiv:2503.21090], k=6: 93999 [Wu et
                  al., certified N*(6) < 94000]) -- both mean a bug
                  here or an error in the literature
  SCREEN_ONLY     --quick: steps 1-2 only (no positivity, no verdict)

usage: verify.py CANDIDATE.json[.gz] [--json] [--quick]
"""
import argparse
import gzip
import hashlib
import json
import math
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
VERSION = "p09-verify-1.1"
PINNED = ["verify.py"]
RECORD = {"k": 6, "n": 90000, "source": "Wu, Liu, Gao, Su, Li, Ge, Ye, ICML 2026 (numerical, not certified)"}
PRIOR_EXPLICIT = {"k": 5, "n": 7265, "source": "arXiv:2503.21090 (Carolan, Childs, Kovacs-Deak, Schaeffer 2025)"}
# largest n with a known exact k-query algorithm whose infeasibility at n+1 is published
PUBLISHED_MAX = {2: (6, "FGGS99/CLP07"), 4: (605, "CLP07"), 5: (7265, "arXiv:2503.21090"),
                 6: (93999, "Wu et al., ICML 2026")}


class Invalid(Exception):
    pass


class Alarm(Exception):
    pass


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()[:16]


def source_hashes():
    return {f: sha256(HERE / f) for f in PINNED}


def check_pins():
    pins = HERE / "PINS.json"
    if not pins.exists():
        return "unpinned"
    want = json.loads(pins.read_text())
    bad = [f for f, h in source_hashes().items() if want.get(f) != h]
    return "ok" if not bad else "MODIFIED: " + ",".join(bad)


def canonical_hash(cand):
    return hashlib.sha256(json.dumps(cand, sort_keys=True).encode()).hexdigest()[:16]


def _int(x, what):
    if isinstance(x, bool) or not isinstance(x, int):
        raise Invalid(f"{what} must be an integer")
    return x


def _intlist(x, what, maxlen, exact=False):
    if not isinstance(x, list):
        raise Invalid(f"{what} must be a list")
    if (exact and len(x) != maxlen) or len(x) > maxlen or not x:
        raise Invalid(f"{what} has length {len(x)} (expected {'exactly ' if exact else 'at most '}{maxlen})")
    for v in x:
        if isinstance(v, bool) or not isinstance(v, int):
            raise Invalid(f"{what} must contain integers only")
    return x


def autocorrelation(h):
    """a[j] = sum_i h[i] h[i+j], j = 0..len(h)-1, exactly (Kronecker substitution)."""
    m = len(h)
    bound = max(abs(v) for v in h)
    # |a_j| <= m * bound^2 < 2^(W-1); W a multiple of 8 for byte slicing
    W = (2 * bound.bit_length() + m.bit_length() + 2 + 7) // 8 * 8
    X = 1 << W

    nb = W // 8
    off = 1 << (W - 2)                       # |h_i| < off, so h_i + off is a digit in [0, X)
    ones = lambda L: (X ** L - 1) // (X - 1)  # sum_{i<L} X^i

    def pack(seq):
        # sum seq[i] X^i, built in linear time from nonnegative digits
        raw = b"".join((v + off).to_bytes(nb, "little") for v in seq)
        return int.from_bytes(raw, "little") - off * ones(len(seq))

    prod = pack(h) * pack(h[::-1])           # coefficient m-1+j is a[j] (and a[-j] mirrors it)
    L = 2 * m - 1
    half = 1 << (W - 1)
    bias = half * ones(L)                    # every digit shifted into [0, X)
    raw = (prod + bias).to_bytes(L * W // 8 + 1, "little")
    out = []
    for l in range(m - 1, L):
        out.append(int.from_bytes(raw[l * nb:(l + 1) * nb], "little") - half)
    return out


def check_candidate(cand, quick=False, **_):
    if not isinstance(cand, dict):
        raise Invalid("candidate must be a JSON object")
    if cand.get("type") != "fggs_chain":
        raise Invalid("type must be 'fggs_chain'")
    k = _int(cand.get("k"), "k")
    n = _int(cand.get("n"), "n")
    if k < 2 or n < 2:
        raise Invalid("need k >= 2 and n >= 2")
    if n > 10 ** 7:
        raise Invalid("n too large for this verifier")
    D = _int(cand.get("den"), "den")
    H = _int(cand.get("h_den"), "h_den")
    if D <= 0 or H <= 0:
        raise Invalid("denominators must be positive")
    if D % n:
        raise Invalid("den must be divisible by n (so that F_n is representable)")
    q = cand.get("q")
    hs = cand.get("h")
    if not isinstance(q, list) or len(q) != k - 1 or not isinstance(hs, list) or len(hs) != k - 1:
        raise Invalid("q and h must each hold k-1 lists")
    q = [_intlist(v, f"q[{i}]", n, exact=True) for i, v in enumerate(q)]
    hs = [_intlist(v, f"h[{i}]", n) for i, v in enumerate(hs)]
    target = f"k{k}"
    # 1-2: the chain, exactly (common denominator D)
    step = D // n
    chain = [[D - j * step for j in range(n)]] + q + [[D] + [0] * (n - 1)]
    for t in range(1, k + 1):
        a, b = chain[t - 1], chain[t]
        if b[0] != a[0]:
            raise Invalid(f"constant term of q_{t} is not 1")
        s = 1 if t % 2 == 0 else -1
        for j in range(1, n):
            if (b[j] - a[j]) + s * (b[n - j] - a[n - j]) != 0:
                raise Invalid(f"q_{t} != q_{t-1} on z^n = {s} (coefficient pair {j}, {n - j})")
    out = {"target": target, "k": k, "n": n}
    if quick:
        out.update(verdict="SCREEN_ONLY", summary=f"SCREEN_ONLY: chain conditions hold for k={k}, n={n}")
        return out
    # 3: positivity certificates
    margins = []
    H2 = H * H
    for t in range(1, k):
        c = chain[t]
        ac = autocorrelation(hs[t - 1]) + [0] * (n - len(hs[t - 1]))
        r = [c[j] * H2 - ac[j] * D for j in range(n)]       # residual * D * H^2
        off = 2 * sum(abs(v) for v in r[1:])
        if r[0] < off:
            raise Invalid(f"q_{t}: residual q - |h|^2 not diagonally dominant "
                          f"(margin {float((r[0] - off) / (D * H2)):.3e})")
        margins.append(float((r[0] - off) / (D * H2)))
    out["margins"] = margins
    # verdict
    const = k / math.log2(n)
    rk, rn = RECORD["k"], RECORD["n"]
    lhs, rhs = n ** rk, rn ** k
    if k < math.log(n) / math.pi - 2:
        raise Alarm(f"k={k} queries for n={n} beats the proven exact lower bound ln(n)/pi")
    if k in PUBLISHED_MAX and n > PUBLISHED_MAX[k][0]:
        raise Alarm(f"k={k}, n={n} exceeds the published feasibility limit "
                    f"{PUBLISHED_MAX[k][0]} ({PUBLISHED_MAX[k][1]})")
    rec = rk / math.log2(rn)
    if lhs > rhs:
        verdict = "IMPROVEMENT"
    elif lhs == rhs:
        verdict = "MATCHES_RECORD"
    else:
        verdict = "BELOW_RECORD"
    out["beats_prior_explicit"] = n ** PRIOR_EXPLICIT["k"] > PRIOR_EXPLICIT["n"] ** k
    out.update(verdict=verdict, constant=const,
               summary=f"{verdict}: exact ordered search with {k} queries for n={n}: "
                       f"{const:.6f} log2 N (record {rec:.6f}, k={rk}, n={rn}); "
                       f"min certified margin {min(margins):.2e}")
    return out


def verify_candidate(cand, exact=False, quick=False, threads=8, timeout=1800, budget=0, seed=1,
                     lock=True, **kw):
    t0 = time.time()
    rep = {"verifier": VERSION, "pins": check_pins(), "sources": source_hashes()}
    try:
        rep["hash"] = canonical_hash(cand) if isinstance(cand, dict) else None
        rep.update(check_candidate(cand, quick=quick))
    except Invalid as e:
        rep.update(verdict="INVALID", reason=str(e))
    except Alarm as e:
        rep.update(verdict="VERIFIER_ALARM", reason=str(e))
    except Exception as e:  # never let bad input crash the board
        rep.update(verdict="INVALID", reason=f"{type(e).__name__}: {e}")
    rep["seconds"] = round(time.time() - t0, 2)
    if "summary" not in rep:
        rep["summary"] = f"{rep['verdict']}: {rep.get('reason', '')}"[:300]
    return rep


def load(path):
    p = pathlib.Path(path)
    data = gzip.decompress(p.read_bytes()) if p.suffix == ".gz" else p.read_bytes()
    return json.loads(data)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--quick", action="store_true", help="chain conditions only (no verdict)")
    a = ap.parse_args()
    rep = verify_candidate(load(a.candidate), quick=a.quick)
    print(json.dumps(rep, indent=1) if a.json else rep["summary"])
    sys.exit(0 if rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") else 1)


if __name__ == "__main__":
    main()
