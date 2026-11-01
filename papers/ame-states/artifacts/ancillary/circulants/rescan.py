#!/usr/bin/env python3
"""Exhaustive symmetric-circulant scans of the paper (Sec. III, Table II), self-contained (stdlib + cc).

The C kernel circ.c (the additive agent's search code, unchanged) enumerates every symmetric circulant
Gamma_ij = c_{min(|i-j|, n-|i-j|)} over GF(q) with first nonzero coefficient 1, i.e. (q^floor(n/2) - 1)/(q - 1)
vectors c, and tests Lemma 2 (Gamma[A, A^c] of full row rank floor(n/2) for every set A of floor(n/2) parties; by the
dihedral symmetry of a symmetric circulant it tests one representative per dihedral orbit of sets).

This script writes the field tables itself (prime fields mod p; GF(4) = F_2[x]/(x^2+x+1), GF(8) = F_2[x]/(x^3+x+1),
GF(9) = F_3[x]/(x^2+x+2); an element of GF(p^2) or GF(2^3) is the integer sum c_i p^i), runs the kernel, and for
the two cells with hits (n = 15, 17 over GF(8)):
  * re-checks every hit (full) or 24 seeded hits (quick) in pure Python with Lemma 2 over ALL sets A (no symmetry);
  * groups the hits into orbits under multipliers i -> u i (u a unit mod n), the Frobenius c -> c^2 and scaling
    c -> lambda c (absorbed by normalizing the first nonzero coefficient to 1); all three preserve the AME property.
usage: rescan.py hits [quick|full]       the n = 15, 17 scans (336 and 48 hits; 28 and 4 orbits)
       rescan.py zero SPEC...            zero-hit scans, SPEC = n,q (e.g. 16,8 19,8 16,7 19,7 16,5 19,5 13,4 15,4)
"""
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))


def field(q):
    """(add, mul) tables for GF(q), q in {2,3,4,5,7,8,9,11,13}."""
    if q in (2, 3, 5, 7, 11, 13):
        return [[(a + b) % q for b in range(q)] for a in range(q)], [[(a * b) % q for b in range(q)] for a in range(q)]
    p, m, poly = {4: (2, 2, [1, 1]), 8: (2, 3, [1, 1, 0]), 9: (3, 2, [2, 1])}[q]  # x^m = -(poly . (1, x, ...))

    def dig(a):
        return [(a // p ** i) % p for i in range(m)]

    def num(ds):
        return sum((d % p) * p ** i for i, d in enumerate(ds))

    def mul(a, b):
        A, B = dig(a), dig(b)
        prod = [0] * (2 * m - 1)
        for i, x in enumerate(A):
            for j, y in enumerate(B):
                prod[i + j] += x * y
        for deg in range(2 * m - 2, m - 1, -1):  # reduce x^deg using x^m = -sum poly_i x^i
            c = prod[deg] % p
            prod[deg] = 0
            for i, pc in enumerate(poly):
                prod[deg - m + i] -= c * pc
        return num(prod[:m])
    add = [[num([x + y for x, y in zip(dig(a), dig(b))]) for b in range(q)] for a in range(q)]
    return add, [[mul(a, b) for b in range(q)] for a in range(q)]


def build():
    td = tempfile.mkdtemp(prefix="circ_")
    exe = os.path.join(td, "circ")
    subprocess.run(["cc", "-O2", "-o", exe, os.path.join(HERE, "circ.c")], check=True)
    return td, exe


def run_scan(exe, td, n, q):
    add, mul = field(q)
    tab = os.path.join(td, f"gf{q}.txt")
    with open(tab, "w") as fh:
        fh.write(f"{q}\n")
        for r in add:
            fh.write(" ".join(map(str, r)) + "\n")
        for r in mul:
            fh.write(" ".join(map(str, r)) + "\n")
    t0 = time.time()
    out = subprocess.run([exe, tab, str(n), "0", "0", str(q ** (n // 2)), "1"], capture_output=True, text=True, check=True).stdout
    hits = [tuple(int(v) for v in l.split()[1:]) for l in out.split("\n") if l.startswith("HIT")]
    done = next(l for l in out.split("\n") if l.startswith("DONE"))
    return hits, int(done.split("tested=")[1].split()[0]), time.time() - t0, (add, mul)


def lemma2_all_sets(c, n, q, add, mul):
    inv = {a: next(b for b in range(1, q) if mul[a][b] == 1) for a in range(1, q)}
    neg = {a: next(b for b in range(q) if add[a][b] == 0) for a in range(q)}
    G = [[c[min(abs(i - j), n - abs(i - j)) - 1] if i != j else 0 for j in range(n)] for i in range(n)]
    k = n // 2
    for A in itertools.combinations(range(n), k):
        Ac = [j for j in range(n) if j not in A]
        M = [[G[i][j] for j in Ac] for i in A]
        rk = 0
        for col in range(len(Ac)):
            pv = next((r for r in range(rk, k) if M[r][col]), None)
            if pv is None:
                continue
            M[rk], M[pv] = M[pv], M[rk]
            iv = inv[M[rk][col]]
            M[rk] = [mul[iv][v] for v in M[rk]]
            for r in range(k):
                if r != rk and M[r][col]:
                    f = M[r][col]
                    M[r] = [add[u][neg[mul[f][v]]] for u, v in zip(M[r], M[rk])]
            rk += 1
        if rk < k:
            return False
    return True


def orbit_key(c, n, q, mul):
    h = n // 2
    inv = {a: next(b for b in range(1, q) if mul[a][b] == 1) for a in range(1, q)}
    keys = []
    for u in range(1, n):
        if gcd(u, n) != 1:
            continue
        cu = []
        for k in range(1, h + 1):
            m = (u * k) % n
            cu.append(c[min(m, n - m) - 1])
        for f in range(3):  # Frobenius x -> x^2 has order 3 on GF(8)
            cf = cu
            for _ in range(f):
                cf = [mul[x][x] for x in cf]
            lead = next(x for x in cf if x)
            keys.append(tuple(mul[inv[lead]][x] for x in cf))
    return min(keys)


def main():
    mode = sys.argv[1]
    td, exe = build()
    if mode == "hits":
        full = len(sys.argv) > 2 and sys.argv[2] == "full"
        res = {}
        for n in (15, 17):
            hits, tested, secs, (add, mul) = run_scan(exe, td, n, 8)
            rnd = random.Random(n)
            sample = hits if full else rnd.sample(hits, min(24, len(hits)))
            bad = sum(1 for c in sample if not lemma2_all_sets(c, n, 8, add, mul))
            orbits = {}
            for c in hits:
                orbits.setdefault(orbit_key(c, n, 8, mul), []).append(c)
            sizes = sorted({len(v) for v in orbits.values()})
            print(f"n={n} q=8: circulants tested {tested}, AME hits {len(hits)}, orbits {len(orbits)} (orbit sizes {sizes}); "
                  f"python re-check over all sets: {len(sample)} hits checked, {bad} failed")
            res[n] = {"tested": tested, "hits": sorted(hits), "orbit_representatives": sorted(orbits)}
        with open(os.path.join(HERE, "out", "hits.json"), "w") as fh:
            json.dump(res, fh)
    elif mode == "zero":
        for spec in sys.argv[2:]:
            n, q = map(int, spec.split(","))
            hits, tested, secs, _ = run_scan(exe, td, n, q)
            print(f"n={n} q={q}: circulants tested {tested}, AME hits {len(hits)}")
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
