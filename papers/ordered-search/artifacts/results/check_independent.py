#!/usr/bin/env python3
"""Independent checker for p09 certificates (python-flint; shares no code with verifier/verify.py).

Different route for each condition:
  * root conditions: P(z) = z^(n-1) (q_t - q_(t-1))(z) is an integer polynomial of degree <= 2n-2,
    and q_t = q_(t-1) on every root of z^n = s iff (z^n - s) divides P  (FLINT divmod);
  * |h|^2: FLINT product h(z) * z^(m-1) h(1/z);
  * positivity: q_t - |h_t|^2 = r_t with r_t[0] >= 2 sum_(j>=1) |r_t[j]|, in fmpq arithmetic.
Also re-derives the constant k / log2 n and the comparison n^5 vs 7265^k (beats the prior explicit
5-query algorithm).

Every acceptance condition raises `Rejected` explicitly (none relies on `assert`, so `python -O`
checks the same conditions). Numbers must be exact JSON integers (not floats or booleans) and both
denominators must be positive. This program checks certificates; it is not a service for hostile
input: it does not limit the size of a file, so a crafted file can exhaust memory or time.

usage: check_independent.py CERT.json[.gz] [...]   (exit status 1 if any certificate is rejected)
"""
import gzip
import json
import math
import pathlib
import sys
import zlib

from flint import fmpq, fmpz, fmpz_poly


class Rejected(Exception):
    pass


def require(cond, msg):
    if not cond:
        raise Rejected(msg)


def exact_int(x):
    return type(x) is int                  # excludes bool, float and int subclasses


def int_list(x, name, lo, hi):
    require(type(x) is list, f"{name} must be a list")
    require(lo <= len(x) <= hi, f"{name} has length {len(x)} (expected {lo}..{hi})")
    require(all(exact_int(v) for v in x), f"{name} must contain exact integers only")
    return x


def load(path):
    p = pathlib.Path(path)
    raw = p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix == ".gz" else raw)


def laurent_times_zn1(coef, n):
    """symmetric Laurent sum_{|j|<n} coef[|j|] z^j, times z^(n-1), as integer coefficient list"""
    return [coef[n - 1 - i] for i in range(n - 1)] + [coef[0]] + [coef[j] for j in range(1, n)]


def check(cert):
    require(type(cert) is dict, "certificate must be a JSON object")
    for key in ("type", "k", "n", "den", "h_den", "q", "h"):
        require(key in cert, f"missing key {key!r}")
    require(cert["type"] == "fggs_chain", "type must be 'fggs_chain'")
    k, n, D, H = cert["k"], cert["n"], cert["den"], cert["h_den"]
    require(all(exact_int(v) for v in (k, n, D, H)), "k, n, den and h_den must be exact integers")
    require(k >= 2 and n >= 2, "need k >= 2 and n >= 2")
    require(D > 0 and H > 0, "denominators must be positive")
    require(D % n == 0, "den must be divisible by n")
    require(type(cert["q"]) is list and len(cert["q"]) == k - 1, "q must hold k-1 lists")
    require(type(cert["h"]) is list and len(cert["h"]) == k - 1, "h must hold k-1 lists")
    q = [int_list(row, f"q[{i}]", n, n) for i, row in enumerate(cert["q"])]
    hs = [int_list(row, f"h[{i}]", 1, n) for i, row in enumerate(cert["h"])]
    chain = [[D * (n - j) // n for j in range(n)]] + q + [[D] + [0] * (n - 1)]
    for t, row in enumerate(chain):
        require(row[0] == D, f"constant term of q_{t} is not 1")
    zn = {1: fmpz_poly([-1] + [0] * (n - 1) + [1]), -1: fmpz_poly([1] + [0] * (n - 1) + [1])}
    for t in range(1, k + 1):
        s = (-1) ** t
        d = [chain[t][j] - chain[t - 1][j] for j in range(n)]
        P = fmpz_poly(laurent_times_zn1(d, n))
        _, r = divmod(P, zn[s])
        require(r == 0, f"root condition fails at t={t}")
    margins = []
    for t in range(1, k):
        h = hs[t - 1]
        m = len(h)
        sq = fmpz_poly(h) * fmpz_poly(h[::-1])          # coefficient m-1+j is sum_i h_i h_(i+j)
        a = [int(sq[m - 1 + j]) if j < m else 0 for j in range(n)]
        r = [fmpq(chain[t][j], D) - fmpq(a[j], H * H) for j in range(n)]
        off = 2 * sum((abs(x) for x in r[1:]), fmpq(0))
        require(r[0] >= off, f"q_{t} residual not diagonally dominant")
        margins.append(float(r[0] - off))
    beats = fmpz(n) ** 5 > fmpz(7265) ** k
    return {"k": k, "n": n, "constant": k / math.log2(n), "beats_prior_explicit": bool(beats),
            "margins": margins}


def main(paths):
    bad = 0
    for path in paths:
        try:
            print(path, json.dumps(check(load(path))))
        except (Rejected, OSError, EOFError, ValueError, zlib.error) as e:  # unreadable, bad gzip or JSON
            bad += 1
            print(path, "REJECTED:", e)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
