#!/usr/bin/env python3
"""Exhaustive minimum T-count of the three 2-dimensional F2-algebras (Python 3 standard library only).

Algebras A = F2[x]/(x^2 + p1 x + p0):
    GF(4)        x^2 = x + 1      (p0, p1) = (1, 1)
    F2[y]/(y^2)  x^2 = 0          (p0, p1) = (0, 0)
    F2 x F2      x^2 = x          (p0, p1) = (0, 1)
The diagonal gate is D = (-1)^f with f(a, b, c) = c . (a b) on N = 6 variables
(a_0, a_1, b_0, b_1, c_0, c_1), basis 1, x.

Model (paper, Sec. II): a set S of nonzero vectors u in F2^6 implements D, up to a diagonal Clifford
gate, if and only if every variable and every pair of variables lies in the support of an even
number of the u, and a triple lies in an odd number exactly when it is a monomial of f (Lemma 1).
A repeated vector cancels in pairs, so a minimum solution is a set.  Writing w in F2^63 for the
indicator vector of S (one coordinate per nonzero u), the condition is the linear system
    sum_{u contains s} w_u = [s is a monomial of f]    for the 41 sets s of 1, 2 or 3 variables.
Its solutions form one coset w0 + C of a code C of dimension 63 - 41 = 22.  The T-count of A is
the minimum Hamming weight in the coset.  This program visits all 2^22 elements by Gray code and
prints the complete weight distribution, so the minimum is exact within the model.  It is not a
lower bound for circuits with intermediate Hadamard gates, measurements or feed-forward.

For GF(4) the minimum 17 was proved by Ruiz et al. (AlphaTensor-Quantum, Z3 solver) and obtained
again by Mazumder (exact coset decoding); this run is an independent confirmation.

With a release directory as argument the program also checks that the released blocks
blocks/block_alpha_n2.json (GF(4), 17) and blocks/block_ring_e2.json (F2[y]/(y^2), 15) lie in the
right coset and attain the minimum.

usage: python3 -I exact/enumerate_small.py [RELEASE_DIR]      (about 15 s)
"""
import json
import pathlib
import sys

N = 6
L = (1 << N) - 1  # coordinates: nonzero u = 1..63, stored at bit u-1

ALGEBRAS = [  # label, (p0, p1), released block or None
    ("GF(4) = F2[x]/(x^2+x+1)", (1, 1), "block_alpha_n2.json"),
    ("F2[y]/(y^2)", (0, 0), "block_ring_e2.json"),
    ("F2 x F2 = F2[x]/(x^2+x)", (0, 1), None),
]


def monomials(p0, p1):
    """Variable sets {a_i, b_j, c_k} (as 6-bit masks) with coefficient 1 in f = c . (a b)."""
    power = [0b01, 0b10, p0 | (p1 << 1)]  # x^0, x^1, x^2 in the basis 1, x
    out = set()
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if power[i + j] >> k & 1:
                    out.add((1 << i) | (1 << (2 + j)) | (1 << (4 + k)))
    return out


def system(p0, p1):
    """Rows (mask over the 63 coordinates, right-hand side) for all sets of 1..3 variables."""
    tgt = monomials(p0, p1)
    rows = []
    for s in range(1, 1 << N):
        if bin(s).count("1") > 3:
            continue
        mask = 0
        for u in range(1, L + 1):
            if u & s == s:
                mask |= 1 << (u - 1)
        rows.append((mask, 1 if s in tgt else 0))
    return rows


def solve(rows):
    """Reduced echelon form; returns (particular solution, basis of the kernel)."""
    red = []  # (pivot column, mask, rhs)
    for mask, rhs in rows:
        for c, m, r in red:
            if mask >> c & 1:
                mask ^= m
                rhs ^= r
        if not mask:
            if rhs:
                raise ValueError("inconsistent system")
            continue
        c = (mask & -mask).bit_length() - 1
        for i, (c2, m2, r2) in enumerate(red):
            if m2 >> c & 1:
                red[i] = (c2, m2 ^ mask, r2 ^ rhs)
        red.append((c, mask, rhs))
    pivots = {c for c, _, _ in red}
    w0 = 0
    for c, _, r in red:
        if r:
            w0 |= 1 << c
    basis = []
    for f in range(L):
        if f in pivots:
            continue
        v = 1 << f
        for c, m, _ in red:
            if m >> f & 1:
                v |= 1 << c
        basis.append(v)
    return w0, basis


def satisfies(w, rows):
    return all(bin(w & mask).count("1") & 1 == rhs for mask, rhs in rows)


def enumerate_coset(w0, basis):
    hist = [0] * (L + 1)
    cur = w0
    best = w0
    bw = bin(w0).count("1")
    hist[bw] += 1
    for g in range(1, 1 << len(basis)):
        cur ^= basis[(g & -g).bit_length() - 1]
        w = bin(cur).count("1")
        hist[w] += 1
        if w < bw:
            bw, best = w, cur
    return hist, best


def block_word(path):
    r = json.loads(path.read_text())
    w = 0
    for s in r["factors"]:
        u = int(s[::-1], 2)  # character i of the string = variable i
        w ^= 1 << (u - 1)
    return w, len(r["factors"])


def main():
    root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None
    fails = 0
    for label, (p0, p1), block in ALGEBRAS:
        rows = system(p0, p1)
        w0, basis = solve(rows)
        ok0 = satisfies(w0, rows) and all(satisfies(v, [(m, 0) for m, _ in rows]) for v in basis)
        hist, best = enumerate_coset(w0, basis)
        mn = min(w for w in range(L + 1) if hist[w])
        okb = satisfies(best, rows) and bin(best).count("1") == mn
        fails += not (ok0 and okb and sum(hist) == 1 << len(basis))
        print(f"{label}: N = {N}, {len(rows)} equations, code dimension {len(basis)}, "
              f"coset size 2^{len(basis)} = {sum(hist)}")
        print(f"  minimum weight {mn}, attained by {hist[mn]} decompositions")
        print("  weight distribution: " + ", ".join(f"{w}: {hist[w]}" for w in range(L + 1) if hist[w]))
        print("  a minimizer (vectors as integers, bit i = variable i): "
              + " ".join(str(u) for u in range(1, L + 1) if best >> (u - 1) & 1)
              + ("   [re-checked against the 41 equations: ok]" if okb else "   [RE-CHECK FAILED]"))
        if root is not None and block is not None:
            wb, m = block_word(root / "blocks" / block)
            good = satisfies(wb, rows) and bin(wb).count("1") == m == mn
            fails += not good
            print(f"  released block blocks/{block}: {m} vectors, in the coset and of minimum weight: "
                  f"{'ok' if good else 'FAIL'}")
    print("EXACT ENUMERATION OK" if not fails else f"{fails} EXACT ENUMERATION CHECK(S) FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
