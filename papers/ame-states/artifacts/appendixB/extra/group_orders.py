#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Group-order facts used by Propositions prop:loc and prop:aut:
  * |GL(4,2)| = 20160 = 2^6 * 3^2 * 5 * 7, so a nontrivial local automorphism has a power of prime order 2, 3, 5 or 7;
  * |Sp(4,2)| = 720 (= |S_6|) and the element orders of Sp(4,2) are 1, 2, 3, 4, 5, 6 only: no element of order 7, so
    the fixed-party twist h is trivial for a 7-cycle (Prop. prop:aut (b));
  * |GL(4,2) wr S_8| = 20160^8 * 8! has prime divisors 2, 3, 5, 7 only (reduction to prime order);
  * the order-3 elements of Sp(4,2): 40 fixed-point-free and 40 with a 2-dim fixed space; order-5 elements are
    fixed-point-free (used in the l = 3 and l = 5 rows of Table VI).
Brute force over all 2^16 4x4 binary matrices; the symplectic form is B(x, y) = x^T J y with J = [[0, I], [I, 0]].
"""
from collections import Counter


def apply(M, v):                     # M = 4 row bitmasks, v = 4-bit vector
    r = 0
    for i in range(4):
        r |= (bin(M[i] & v).count("1") & 1) << i
    return r


def B(u, v):
    return ((u & 1) & (v >> 2 & 1)) ^ ((u >> 1 & 1) & (v >> 3 & 1)) ^ ((u >> 2 & 1) & (v & 1)) ^ ((u >> 3 & 1) & (v >> 1 & 1))


def compose(f, g):
    return tuple(f[g[x]] for x in range(16))


ident = tuple(range(16))
gl, sp = [], []
for code in range(1 << 16):
    M = [(code >> (4 * i)) & 15 for i in range(4)]
    img = tuple(apply(M, v) for v in range(16))
    if len(set(img)) != 16:
        continue
    gl.append(img)
    if all(B(img[u], img[v]) == B(u, v) for u in (1, 2, 4, 8) for v in (1, 2, 4, 8)):
        sp.append(img)


def order(f):
    k, g = 1, f
    while g != ident:
        g = compose(f, g)
        k += 1
    return k


def factor(n):
    out, p = {}, 2
    while n > 1:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    return out


def fmt(f):
    return " * ".join("%d^%d" % (p, e) if e > 1 else "%d" % p for p, e in sorted(f.items()))


ord_gl = Counter(order(g) for g in gl)
ord_sp = Counter(order(g) for g in sp)
print("|GL(4,2)| = %d = %s" % (len(gl), fmt(factor(len(gl)))))
print("element orders in GL(4,2): %s" % dict(sorted(ord_gl.items())))
print("|Sp(4,2)| = %d = %s" % (len(sp), fmt(factor(len(sp)))))
print("element orders in Sp(4,2): %s" % dict(sorted(ord_sp.items())))
print("Sp(4,2) has an element of order 7: %s" % (7 in ord_sp))
fix = lambda g: sum(1 for x in range(16) if g[x] == x)
o3 = Counter(fix(g) for g in sp if order(g) == 3)
o5 = Counter(fix(g) for g in sp if order(g) == 5)
print("order-3 elements of Sp(4,2) by number of fixed vectors: %s  (1 = fixed-point-free, 4 = 2-dim fixed space)"
      % dict(sorted(o3.items())))
print("order-5 elements of Sp(4,2) by number of fixed vectors: %s" % dict(sorted(o5.items())))
fw = factor(len(gl) ** 8 * 40320)
print("prime divisors of |GL(4,2) wr S_8| = 20160^8 * 8!: %s" % sorted(fw))

# ---- Prop. prop:aut (c): the sector list.  Cycle types of permutations of prime order l of the 8 parties, and the
#      possible twists h_j (h_j^l = 1) on the parties fixed by sigma, classified by fixed-space dimension (order 3) or by
#      rank(1 + h) (order 2).
import itertools as _it
types = set()
for perm in _it.permutations(range(8)):
    seen, cyc = set(), []
    for i in range(8):
        if i not in seen:
            j, L = i, 0
            while j not in seen:
                seen.add(j); j = perm[j]; L += 1
            cyc.append(L)
    lens = sorted(cyc, reverse=True)
    o = 1
    for L in lens:
        o = o * L // __import__("math").gcd(o, L)
    if o in (2, 3, 5, 7):
        types.add((o, tuple(lens)))
for o in (2, 3, 5, 7):
    ts = sorted((t for oo, t in types if oo == o), reverse=True)
    print("prime order %d permutations of 8 parties, cycle types: %s"
          % (o, ", ".join("(" + ",".join(map(str, t)) + ")" for t in ts)))
inv = [g for g in sp if order(g) == 2]
rank1h = lambda g: 4 - (fix(g).bit_length() - 1)          # rank(1 + h) = 4 - dim ker(1 + h), |ker| = number of fixed vectors
print("involutions of Sp(4,2) by rank of 1 + h: %s" % dict(sorted(Counter(rank1h(g) for g in inv).items())))
print("order-3 twists: h = 1, h with 2-dim fixed space (%d), h fixed-point-free (%d); order-5 twists: fixed-point-free (%d)"
      % (o3[4], o3[1], o5[1]))
