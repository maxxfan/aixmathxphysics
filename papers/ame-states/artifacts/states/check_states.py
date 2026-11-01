#!/usr/bin/env python3
"""Field-level check of the six explicit AME states (paper Theorems 1 and 2), stdlib only, own GF(2^m) arithmetic.

For each entry of printed.json:
  graph (F_q circulant):  Lemma 2 -- Gamma[A, A^c] has full row rank floor(n/2) over GF(q) for every set A of
                          floor(n/2) parties. Gamma is circulant, so rank(Gamma[A+s, A^c+s]) = rank(Gamma[A, A^c]); we test one
                          representative per rotation orbit (the lexicographically smallest rotation) and report orbit
                          counts and the number of sets they cover.
  hermitian (F_256):      A conj(A)^T = I (characteristic 2, conj x -> x^16), so [I | A] generates a Hermitian self-dual
                          code, and every choice of n/2 columns of [I | A] is invertible (MDS, distance n/2 + 1), which is
                          AME(n, 16) by the Hermitian -> stabilizer map. All C(n, n/2) column selections are tested.
  twist (AME(21,16)):     c_{8j mod 21} = c_j^4 for all j.
Then (optional, --convert) it rebuilds the canonical F_2 rows from the printed data with the documented converter
qlib_ame.py (the code that produced the certified candidates) and compares the F_2 row spaces with candidates/.
"""
import itertools
import json
import pathlib
import sys
from math import comb

HERE = pathlib.Path(__file__).resolve().parent
MOD = {8: 0b1011, 16: 0b10011, 256: 0b100011101}


class GF2m:
    def __init__(self, q):
        self.q, self.mod = q, MOD[q]
        self.exp, self.log = [0] * (2 * q), [0] * q
        x = 1
        for i in range(q - 1):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & q:
                x ^= self.mod
        for i in range(q - 1, 2 * q):
            self.exp[i] = self.exp[i - (q - 1)]
        assert len(set(self.exp[:q - 1])) == q - 1, "modulus is not primitive"

    def mul(self, a, b):
        return 0 if a == 0 or b == 0 else self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        return self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)]

    def pw(self, a, e):
        return 0 if a == 0 else self.exp[(self.log[a] * e) % (self.q - 1)]


def rank(F, M):
    M = [list(r) for r in M]
    rk = 0
    for c in range(len(M[0])):
        pv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if pv is None:
            continue
        M[rk], M[pv] = M[pv], M[rk]
        iv = F.inv(M[rk][c])
        M[rk] = [F.mul(iv, v) for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [u ^ F.mul(f, v) for u, v in zip(M[i], M[rk])]
        rk += 1
    return rk


def check_graph(st):
    F, n, c = GF2m(st["q"]), st["n"], st["c"]
    G = [[c[min(abs(i - j), n - abs(i - j)) - 1] if i != j else 0 for j in range(n)] for i in range(n)]
    k = n // 2
    total = reps = covered = bad = 0
    for A in itertools.combinations(range(n), k):
        total += 1
        rots = [tuple(sorted((x + s) % n for x in A)) for s in range(n)]
        if A != min(rots):
            continue
        reps += 1
        covered += len(set(rots))
        Ac = [j for j in range(n) if j not in A]
        if rank(F, [[G[i][j] for j in Ac] for i in A]) < k:
            bad += 1
    assert covered == total == comb(n, k)
    line = (f"{st['cell']}: graph over GF({st['q']}), {total} sets of {k} parties in {reps} rotation orbits; "
            f"rank-deficient orbits {bad}")
    ok = bad == 0
    if "twist" in st:
        u, e = st["twist"]["u"], st["twist"]["frobenius_power"]
        cf = lambda x: c[min(x % n, n - x % n) - 1] if x % n else 0  # noqa: E731
        tw = all(cf(u * x) == F.pw(cf(x), e) for x in range(1, n))
        line += f"; twist c_(8j) = c_j^4: {'holds' if tw else 'FAILS'}"
        ok &= tw
    return ok, line


def circ(r):
    m = len(r)
    return [[r[(j - i) % m] for j in range(m)] for i in range(m)]


def check_hermitian(st):
    F, n = GF2m(256), st["n"]
    b = st["blocks"]
    if "A" in b:
        A = circ(b["A"])
    else:
        C00, C01, C11 = circ(b["C00"]), circ(b["C01"]), circ(b["C11"])
        C10 = [list(r) for r in zip(*C01)]
        A = [C00[i] + C01[i] for i in range(5)] + [C10[i] + C11[i] for i in range(5)]
    h = n // 2
    conjA = [[F.pw(x, 16) for x in row] for row in A]
    unitary = all(
        (lambda s: s == (1 if i == j else 0))(
            __import__("functools").reduce(lambda acc, t: acc ^ F.mul(A[i][t], conjA[j][t]), range(h), 0))
        for i in range(h) for j in range(h))
    Gm = [[1 if j == i else 0 for j in range(h)] + A[i] for i in range(h)]
    total = bad = 0
    for S in itertools.combinations(range(n), h):
        total += 1
        if rank(F, [[Gm[i][j] for j in S] for i in range(h)]) < h:
            bad += 1
    line = (f"{st['cell']}: Hermitian [I|A] over GF(256): A conj(A)^T = I {'holds' if unitary else 'FAILS'}; "
            f"{total} selections of {h} columns, singular {bad}")
    return unitary and bad == 0, line


def convert_check(st):
    sys.path.insert(0, str(HERE))
    import qlib_ame as ame
    n = st["n"]

    def rowspace(rows):
        basis = {}
        for r in rows:
            v = int("".join(map(str, r)), 2)
            while v:
                h = v.bit_length() - 1
                if h in basis:
                    v ^= basis[h]
                else:
                    basis[h] = v
                    break
        # full reduction
        for h in sorted(basis):
            for g in basis:
                if g != h and (basis[g] >> h) & 1:
                    basis[g] ^= basis[h]
        return sorted(basis.values())
    if st["kind"] == "graph":
        F = ame.GF(st["q"])
        c = st["c"]
        G = [[c[min(abs(i - j), n - abs(i - j)) - 1] if i != j else 0 for j in range(n)] for i in range(n)]
        rows = ame.canonical_from_fq([[int(k == i) for k in range(n)] + G[i] for i in range(n)], F)
    else:
        F = ame.GF(256)
        b = st["blocks"]
        if "A" in b:
            A = circ(b["A"])
        else:
            C00, C01, C11 = circ(b["C00"]), circ(b["C01"]), circ(b["C11"])
            C10 = [list(r) for r in zip(*C01)]
            A = [C00[i] + C01[i] for i in range(5)] + [C10[i] + C11[i] for i in range(5)]
        h = n // 2
        rows = ame.canonical_from_hermitian([[int(j == i) for j in range(h)] + A[i] for i in range(h)], F, 16)
    cand = json.loads((HERE / st["candidate"]).read_text())
    crow = [[int(ch) for ch in r] if isinstance(r, str) else r for r in cand["generators"]]
    same = rowspace(rows) == rowspace(crow)
    return same, f"{st['cell']}: printed data -> canonical F_2 rows (qlib_ame.py) {'equal' if same else 'DIFFER'} to {st['candidate']}"


def main():
    data = json.loads((HERE / "printed.json").read_text())
    ok = True
    for st in data["states"]:
        good, line = (check_graph if st["kind"] == "graph" else check_hermitian)(st)
        print(("PASS  " if good else "FAIL  ") + line, flush=True)
        ok &= good
        if "--convert" in sys.argv:
            good, line = convert_check(st)
            print(("PASS  " if good else "FAIL  ") + line, flush=True)
            ok &= good
    print("STATES FIELD-LEVEL", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
