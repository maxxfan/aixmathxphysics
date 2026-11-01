#!/usr/bin/env python3
"""[release] SUPPLEMENTARY script written for the Appendix-B release (NOT part of the original run; stdlib only).

Table VI row l = 5, sigma = (0 1 2 3 4), order-5 twists h5, h6 on parties 5, 6, party 7 untwisted.  S_nf is an
Omega-invariant (F_16-linear) Lagrangian L of V_nf = C + W (C = cycle part, F_16-dim 4; W = W5 + W6, F_16-dim 2), i.e. a
generator of the Hermitian polar space H(5,16).  L_0 = L cap C has F_16-dim 1 or 2 (L projects to W, dim <= 2; C has no
isotropic 3-space).  The auditor's c5graph.py and the pairs' c5twist.py enumerate the 331,500 Lagrangians with dim L_0 = 1
("graph type", 1105 x 300).  The remaining ones (dim L_0 = 2) were excluded BY HAND in #184 / argument.md 5.6:
"L_0 one of the (q+1)(q^3+1) = 325 Lagrangians of C, plus one of the q+1 = 5 isotropic points of W, which has both
coordinates nonzero -- these contain the weight-2 word (0, u_5, u_6); 331 500 + 325*5 = 333 125 = (q+1)(q^3+1)(q^5+1)".
This script makes that step executable in the auditor's coordinates (same Omega as c5graph.py with h5idx = h6idx = 0, or the
indices given on the command line):
  (a) enumerate the isotropic F_16-lines of C (expect 1105) and all totally isotropic F_16-planes of C (expect 325);
  (b) enumerate the isotropic F_16-points of W (expect 5) and check that each has both party coordinates nonzero;
  (c) build all 325 x 5 Lagrangians L = L_0 + <u>, check isotropy, Omega-invariance and F_2-dim 12, and compute the minimum
      party weight of all 4095 nonzero words of every L (expect max over L of min weight = 2);
  (d) print the count identity 1105*300 + 325*5 = 333125 = (q+1)(q^3+1)(q^5+1), q = 4.
Usage: c5_nongraph.py [h5idx] [h6idx]
"""
import sys

H5IDX = int(sys.argv[1]) if len(sys.argv) > 1 else 0
H6IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0


def popc(x):
    return bin(x).count("1")


def mat_apply(M, v):
    r = 0
    for i in range(4):
        r |= (popc(M[i] & v) & 1) << i
    return r


def B4(u, v):
    return ((u & 1) & (v >> 1 & 1)) ^ ((u >> 1 & 1) & (v & 1)) ^ ((u >> 2 & 1) & (v >> 3 & 1)) ^ ((u >> 3 & 1) & (v >> 2 & 1))


# Sp(4,2) in the same enumeration order as c5graph.py / c5audit.py, so that the h-indices mean the same elements
sp = []
for code in range(1 << 16):
    M = [(code >> (4 * i)) & 15 for i in range(4)]
    img = [mat_apply(M, v) for v in range(16)]
    if len(set(img)) != 16:
        continue
    if all(B4(img[u], img[v]) == B4(u, v) for u in (1, 2, 4, 8) for v in (1, 2, 4, 8)):
        sp.append(tuple(img))
assert len(sp) == 720
ident = tuple(range(16))


def compose(f, g):
    return tuple(f[g[v]] for v in range(16))


def order(f):
    k, g = 1, f
    while g != ident:
        g = compose(f, g)
        k += 1
    return k


ord5 = [f for f in sp if order(f) == 5]
assert len(ord5) == 144
h5, h6 = ord5[H5IDX % 144], ord5[H6IDX % 144]


def omega(x):
    nib = [(x >> (4 * j)) & 15 for j in range(8)]
    y = [0] * 8
    for i in range(5):
        y[(i + 1) % 5] = nib[i]
    y[5] = h5[nib[5]]
    y[6] = h6[nib[6]]
    y[7] = nib[7]
    return sum(y[j] << (4 * j) for j in range(8))


M5 = 0x55555555


def bform(x, y):   # sum over parties of B4: parity of x & swap(y) (swap = exchange the two bits of each qubit pair 0<->1, 2<->3)
    sy = ((y & M5) << 1) | ((y >> 1) & M5)
    return popc(x & sy) & 1


def om_pow(x, m):
    for _ in range(m):
        x = omega(x)
    return x


def f16mult(x, c):  # c in 0..15 acts as sum_m c_m Omega^m (F_2[Omega] = F_16 on V_nf)
    s = 0
    for m in range(4):
        if c >> m & 1:
            s ^= om_pow(x, m)
    return s


def reduce_basis(vecs):
    basis = []
    for v in vecs:
        for p, b in basis:
            if v >> p & 1:
                v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, (b ^ v) if (b >> p & 1) else b) for q, b in basis]
            basis.append((p, v))
    return [b for _, b in basis]


def all_combos(basis):
    out = [0]
    for b in basis:
        out = out + [x ^ b for x in out]
    return out


def party_weight(x):
    return sum(1 for j in range(8) if (x >> (4 * j)) & 15)


def line_iso(x):   # the F_16-line of x is totally isotropic iff B(x, Omega x) = B(x, Omega^2 x) = 0
    return bform(x, omega(x)) == 0 and bform(x, om_pow(x, 2)) == 0


def h_orth(x, y):  # F_16-Hermitian orthogonality of the lines of x and y: B(x, Omega^m y) = 0, m = 0..3
    return all(bform(x, om_pow(y, m)) == 0 for m in range(4))


# ---- (a) C, isotropic lines, totally isotropic planes
Cb = reduce_basis([omega(1 << k) ^ (1 << k) for k in range(20)])
Wb = reduce_basis([omega(1 << k) ^ (1 << k) for k in range(20, 28)])
assert len(Cb) == 16 and len(Wb) == 8
Cv = all_combos(Cb)
line_of = {}
lines = []
for x in Cv:
    if x == 0 or x in line_of or not line_iso(x):
        continue
    mult = [f16mult(x, c) for c in range(1, 16)]
    assert len(set(mult)) == 15
    for v in mult:
        line_of[v] = len(lines)
    lines.append(x)
print("C: F_2-dim %d; isotropic nonzero vectors %d -> isotropic F_16-lines %d (expected 1105 = (q^4-1)(q^3+1)/(q^2-1))"
      % (len(Cb), len(line_of), len(lines)))
planes = set()
for i, x in enumerate(lines):
    covered = set()
    for j, y in enumerate(lines):
        if j == i or j in covered or not h_orth(x, y):
            continue
        pl = set()
        for a in range(16):
            ax = f16mult(x, a)
            for b in range(16):
                if a == 0 and b == 0:
                    continue
                v = ax ^ f16mult(y, b)
                assert v in line_of, "a vector of a plane spanned by two orthogonal isotropic lines is not isotropic"
                pl.add(line_of[v])
        assert len(pl) == 17
        covered |= pl
        planes.add(frozenset(pl))
planes = sorted(planes, key=lambda s: sorted(s))
print("C: totally isotropic F_16-planes (Lagrangians of C) %d (expected 325 = (q+1)(q^3+1)); lines per plane 17, planes per line %s"
      % (len(planes), sorted({sum(1 for p in planes if i in p) for i in range(len(lines))})))

# ---- (b) W, isotropic points
Wv = all_combos(Wb)
wpts = []
seenw = set()
for u in Wv:
    if u == 0 or u in seenw or not line_iso(u):
        continue
    mult = {f16mult(u, c) for c in range(1, 16)}
    seenw |= mult
    wpts.append(u)
both = sum(1 for u in wpts if all(((v >> 20) & 15) and ((v >> 24) & 15) for v in {f16mult(u, c) for c in range(1, 16)}))
print("W: F_2-dim %d; isotropic F_16-points %d (expected q+1 = 5); with both party coordinates (5 and 6) nonzero: %d"
      % (len(Wb), len(wpts), both))

# ---- (c) the 325 x 5 Lagrangians L = L_0 + <u>
hist = {}
nbad = 0
for pl in planes:
    i, j = sorted(pl)[:2]
    x, y = lines[i], lines[j]
    L0 = reduce_basis([om_pow(x, m) for m in range(4)] + [om_pow(y, m) for m in range(4)])
    assert len(L0) == 8
    for u in wpts:
        Lb = reduce_basis(L0 + [om_pow(u, m) for m in range(4)])
        words = all_combos(Lb)
        wset = set(words)
        ok = (len(Lb) == 12 and all(bform(a, b) == 0 for a in Lb for b in Lb)   # Lagrangian of V_nf (F_2-dim 24)
              and all(omega(b) in wset for b in Lb))                            # Omega-invariant
        if not ok:
            nbad += 1
        mw = min(party_weight(v) for v in words if v)
        hist[mw] = hist.get(mw, 0) + 1
print("non-graph Lagrangians L_0 + <u>: %d (expected 325 * 5 = 1625); failing isotropy/rank-12/Omega-invariance: %d"
      % (sum(hist.values()), nbad))
print("min-party-weight histogram over these L:", dict(sorted(hist.items())), " max over L of min weight:", max(hist))
# ---- (d) count identity
q = 4
print("count: graph type 1105*300 = %d (c5graph.py / c5twist.py) + non-graph %d = %d; (q+1)(q^3+1)(q^5+1) = %d"
      % (1105 * 300, len(planes) * len(wpts), 1105 * 300 + len(planes) * len(wpts), (q + 1) * (q ** 3 + 1) * (q ** 5 + 1)))
