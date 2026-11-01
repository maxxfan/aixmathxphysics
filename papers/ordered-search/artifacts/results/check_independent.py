#!/usr/bin/env python3
"""Independent checker for p09 certificates (python-flint; shares no code with verifier/verify.py).

Different route for each condition:
  * root conditions: P(z) = z^(n-1) (q_t - q_(t-1))(z) is an integer polynomial of degree <= 2n-2,
    and q_t = q_(t-1) on every root of z^n = s iff (z^n - s) divides P  (FLINT divmod);
  * |h|^2: FLINT product h(z) * z^(m-1) h(1/z);
  * positivity: q_t - |h_t|^2 = r_t with r_t[0] >= 2 sum_(j>=1) |r_t[j]|, in fmpq arithmetic.
Also re-derives the constant k / log2 n and the comparison n^5 vs 7265^k.

usage: check_independent.py CERT.json[.gz] [...]
"""
import gzip
import json
import math
import pathlib
import sys

from flint import fmpq, fmpz, fmpz_poly


def load(path):
    p = pathlib.Path(path)
    raw = p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix == ".gz" else raw)


def laurent_times_zn1(coef, n):
    """symmetric Laurent sum_{|j|<n} coef[|j|] z^j, times z^(n-1), as integer coefficient list"""
    return [coef[n - 1 - i] for i in range(n - 1)] + [coef[0]] + [coef[j] for j in range(1, n)]


def check(cert):
    assert cert["type"] == "fggs_chain"
    k, n, D, H = cert["k"], cert["n"], cert["den"], cert["h_den"]
    assert k >= 2 and n >= 2 and D > 0 and H > 0 and D % n == 0
    assert len(cert["q"]) == k - 1 and len(cert["h"]) == k - 1
    chain = [[D * (n - j) // n for j in range(n)]] + cert["q"] + [[D] + [0] * (n - 1)]
    for row in chain:
        assert len(row) == n and all(type(v) is int for v in row)
        assert row[0] == D, "constant term"
    zn = {1: fmpz_poly([-1] + [0] * (n - 1) + [1]), -1: fmpz_poly([1] + [0] * (n - 1) + [1])}
    for t in range(1, k + 1):
        s = (-1) ** t
        d = [chain[t][j] - chain[t - 1][j] for j in range(n)]
        P = fmpz_poly(laurent_times_zn1(d, n))
        q, r = divmod(P, zn[s])
        assert r == 0, f"root condition fails at t={t}"
    margins = []
    for t in range(1, k):
        h = cert["h"][t - 1]
        assert 1 <= len(h) <= n and all(type(v) is int for v in h)
        m = len(h)
        sq = fmpz_poly(h) * fmpz_poly(h[::-1])          # coefficient m-1+j is sum_i h_i h_(i+j)
        a = [int(sq[m - 1 + j]) if j < m else 0 for j in range(n)]
        r = [fmpq(chain[t][j], D) - fmpq(a[j], H * H) for j in range(n)]
        off = 2 * sum((abs(x) for x in r[1:]), fmpq(0))
        assert r[0] >= off, f"q_{t} residual not diagonally dominant"
        margins.append(float(r[0] - off))
    beats = fmpz(n) ** 5 > fmpz(7265) ** k
    return {"k": k, "n": n, "constant": k / math.log2(n), "beats_record": bool(beats),
            "margins": margins}


if __name__ == "__main__":
    bad = 0
    for path in sys.argv[1:]:
        try:
            print(path, json.dumps(check(load(path))))
        except AssertionError as e:
            bad += 1
            print(path, "REJECTED:", e)
    sys.exit(1 if bad else 0)
