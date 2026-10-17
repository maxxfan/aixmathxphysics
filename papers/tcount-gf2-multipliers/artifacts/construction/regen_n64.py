#!/usr/bin/env python3
"""Independent regeneration of the raw n=64 network (33 GF(16) blocks of 39 T gates) from this bundle alone:
construction/curve64.json and blocks/block_alpha_n4.json.  It pulls the 39 block factors back through
L_k = (lambda_k, lambda_k, gamma_k), with gamma_k(c)_i = parity(c AND omega_k theta^i), checks the result with its
own tensor test, and compares sha256(json.dumps(sorted factor ints, bit i = variable i)) with
curve64.json["raw_network_sha256"].  Written by the project's auditor agent; Python 3 standard library only.
usage: python3 -I construction/regen_n64.py [BUNDLE_DIR]"""
import hashlib
import json
import os
import random
import sys

BUNDLE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
d = json.load(open(os.path.join(BUNDLE, "construction", "curve64.json")))
blk = json.load(open(os.path.join(BUNDLE, "blocks", "block_alpha_n4.json")))
assert blk["modulus"] == d["K_modulus"]
Vm = int(d["V_modulus"], 2)
Km = int(d["K_modulus"], 2)
N = 64


def vmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> N & 1:
            a ^= Vm
    return r


theta = int(d["theta"], 16)
thp = [1]
for _ in range(3):
    thp.append(vmul(thp[-1], theta))


def kV(c):  # K element (4-bit int) -> V
    v = 0
    for i in range(4):
        if c >> i & 1:
            v ^= thp[i]
    return v


def kmul(a, b):
    r = 0
    for i in range(4):
        if b >> i & 1:
            r ^= a << i
    for i in range(7, 3, -1):
        if r >> i & 1:
            r ^= Km << (i - 4)
    return r


beta = [int(b, 16) for b in d["beta"]]
cols = [vmul(thp[r], beta[s]) for s in range(16) for r in range(4)]  # col index s*4+r
# solve: a = XOR_{idx in S} cols[idx]; find matrix Inv with coords(a) = bits; Gaussian elimination
piv = {}  # leading bit -> (vec, comb)
for idx, v in enumerate(cols):
    comb = 1 << idx
    while v:
        h = v.bit_length() - 1
        if h in piv:
            v ^= piv[h][0]
            comb ^= piv[h][1]
        else:
            piv[h] = (v, comb)
            break
    else:
        raise SystemExit("cols dependent")


def coords(a):  # 64-bit int of V -> 64-bit int of coefficients w.r.t. cols (bit s*4+r)
    comb = 0
    while a:
        h = a.bit_length() - 1
        a ^= piv[h][0]
        comb ^= piv[h][1]
    return comb


# linear maps on the 64 bits of a: alpha_s (K coord) ; lambda'_k(a) = sum_s alpha_s * lam[k][s]
# matrix Lk[i] = 64-bit mask over a-bits j such that bit i of lambda'_k(e_j) = 1
ecoords = [coords(1 << j) for j in range(N)]  # coefficient vector of basis vector e_j


def lam_of_coords(k, cv):
    out = 0
    for s in range(16):
        al = (cv >> (4 * s)) & 15
        out ^= kmul(al, d["lambda"][k][s])
    return out


Lk = []
for k in range(33):
    rows = [0] * 4
    for j in range(N):
        v = lam_of_coords(k, ecoords[j])
        for i in range(4):
            if v >> i & 1:
                rows[i] |= 1 << j
    Lk.append(rows)

# gamma_k(c)_i = parity(c & (omega_k * theta^i))
omega = [int(o, 16) for o in d["omega"]]
Gk = []
for k in range(33):
    Gk.append([vmul(omega[k], thp[i]) for i in range(4)])

# sanity: the identity a*b = sum_k lam_k(a) lam_k(b) omega_k on random inputs
rng = random.Random(5)
for _ in range(20):
    a, b = rng.getrandbits(N), rng.getrandbits(N)
    acc = 0
    for k in range(33):
        la = sum(((bin(Lk[k][i] & a).count("1") & 1) << i) for i in range(4))
        lb = sum(((bin(Lk[k][i] & b).count("1") & 1) << i) for i in range(4))
        acc ^= vmul(omega[k], kV(kmul(la, lb)))
    assert acc == vmul(a, b), "identity fails"
print("identity a*b = sum lam lam omega OK on 20 random pairs")

# block factors: string chars: a_0..a_3, b_0..b_3, c_0..c_3 (character i = variable i)
bf = [int(s[::-1], 2) for s in blk["factors"]]  # bit i = variable i
assert len(bf) == 39
out = {}
for k in range(33):
    for u in bf:
        ux, uy, uz = u & 15, (u >> 4) & 15, (u >> 8) & 15
        fa = fb = fc = 0
        for i in range(4):
            if ux >> i & 1:
                fa ^= Lk[k][i]
            if uy >> i & 1:
                fb ^= Lk[k][i]
            if uz >> i & 1:
                fc ^= Gk[k][i]
        vec = fa | (fb << N) | (fc << (2 * N))
        if vec:
            out[vec] = out.get(vec, 0) ^ 1
F = sorted(v for v, m in out.items() if m)
print("pulled-back factors: %d (before cancellation %d)" % (len(F), 33 * 39))

# own tensor check (singles, pairs even; triples odd exactly at monomials a_i b_j c_k of c.(ab))
NN = 3 * N
# target: for i,j: x^(i+j) mod Vm bits k -> triple (i, N+j, 2N+k)
single = [0] * NN
pair = {}
trip = {}
for u in F:
    s = [i for i in range(NN) if u >> i & 1]
    for x in range(len(s)):
        single[s[x]] ^= 1
        for y in range(x + 1, len(s)):
            key = (s[x], s[y])
            pair[key] = pair.get(key, 0) ^ 1
            above = u & ~((1 << (s[y] + 1)) - 1)
            if above:
                trip[key] = trip.get(key, 0) ^ above
bad_single = sum(single)
bad_pair = sum(1 for v in pair.values() if v)
want = {}
xp = [1]
for e in range(1, 2 * N):
    v = xp[-1] << 1
    if v >> N & 1:
        v ^= Vm
    xp.append(v)
for i in range(N):
    for j in range(N):
        v = xp[i + j]
        m = 0
        for kk in range(N):
            if v >> kk & 1:
                m |= 1 << (2 * N + kk)
        if m:
            want[(i, N + j)] = want.get((i, N + j), 0) | m
keys = set(trip) | set(want)
bad_trip = sum(1 for k in keys if trip.get(k, 0) != want.get(k, 0))
print("singles odd: %d, pairs odd: %d, triple mismatches: %d" % (bad_single, bad_pair, bad_trip))
print("VALID decomposition of GF(2^64) multiplication" if not (bad_single or bad_pair or bad_trip) else "INVALID")

h = hashlib.sha256(json.dumps(F).encode()).hexdigest()
ok = not (bad_single or bad_pair or bad_trip) and len(F) == d["raw_network_tcount"] and h == d["raw_network_sha256"]
print("raw network: %d factors, sha256 %s... %s" % (len(F), h[:16], "matches curve64.json" if h == d["raw_network_sha256"] else "DIFFERS"))
print("REGENERATION OK" if ok else "REGENERATION FAILED")
sys.exit(0 if ok else 1)
