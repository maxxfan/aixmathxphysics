"""crosswalk.py (release/frames, written for this release; standard library only; shares no code with impl1 or impl2).

Matches the V_0 class representatives of the two implementations by a gauge INVARIANT, not by their counts, and then checks that
matched classes have equal orbit sizes and equal search counts.

Invariant (PROOF.md section 5, argument.md 5.12): for V_0 with basis b_0..b_3 (16-bit vectors, party l = nibble l, coordinates
(x1,x2,z1,z2) = bits 0..3 of a nibble, symplectic form J pairing bit 0<->2 and 1<->3), let Q_l be the Gram matrix of J on the
projections of b_0..b_3 to party l (l = 0..3).  The class of V_0 under Sp(4,2)^4 x| S_4 (local Cliffords and permutations of the
four B-parties) is the orbit of the tuple (Q_0,..,Q_3) under GL(4,2) (simultaneous congruence, change of basis of V_0) x S_4
(permutation of the tuple).  Canonical form used here: the lexicographically smallest sorted tuple of 6-bit form codes over all
20160 elements of GL(4,2).  For the Sp(4,2)-only representatives of impl2 the same invariant gives their full-gauge class.

usage:
  crosswalk.py check V0_M4 V0REPS8 V0REPS8_SP I1_T5.canon I1_M124.canon I2_10.canon I2_SP72.canon
  crosswalk.py table V0_M4 V0REPS8          (human-readable list of the ten representatives of both implementations)
"""
import itertools
import re
import sys

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def jform(x, y):
    """symplectic form on one party F_2^4, coordinates (x1,x2,z1,z2): x1 z1' + z1 x1' + x2 z2' + z2 x2'."""
    s = ((y >> 2) & 3) | ((y & 3) << 2)
    return bin(x & s).count("1") & 1


def gram_code(vecs):
    """6-bit code of the alternating Gram matrix of 4 vectors of F_2^4 under J."""
    return sum(jform(vecs[i], vecs[j]) << k for k, (i, j) in enumerate(PAIRS))


def code_bil(code):
    """table B[a][b] = Q(a,b) for the alternating form with this code."""
    Q = [[0] * 4 for _ in range(4)]
    for k, (i, j) in enumerate(PAIRS):
        if code >> k & 1:
            Q[i][j] = Q[j][i] = 1
    B = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for b in range(16):
            B[a][b] = sum(Q[i][j] for i in range(4) for j in range(4) if (a >> i & 1) and (b >> j & 1)) & 1
    return B


def rank_cols(cols):
    basis = []
    for x in cols:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
            basis.sort(reverse=True)
    return len(basis)


def nondeg(code):
    # Pfaffian of the alternating 4x4 matrix: a01 a23 + a02 a13 + a03 a12
    a = [(code >> k) & 1 for k in range(6)]
    return (a[0] & a[5]) ^ (a[1] & a[4]) ^ (a[2] & a[3])


GL = [c for c in itertools.product(range(1, 16), repeat=4) if rank_cols(list(c)) == 4]
assert len(GL) == 20160
CODES = [c for c in range(64) if nondeg(c)]
assert len(CODES) == 28
BIL = {c: code_bil(c) for c in range(64)}
# congruence action: ACT[g][code] = code of g^T Q g, i.e. entries Q(g e_i, g e_j)
ACT = [{c: sum(BIL[c][g[i]][g[j]] << k for k, (i, j) in enumerate(PAIRS)) for c in CODES} for g in GL]


def invariant(basis):
    qs = [gram_code([(b >> (4 * l)) & 15 for b in basis]) for l in range(4)]
    assert all(nondeg(q) for q in qs), "a projection of V_0 is not an isomorphism"
    assert nondeg(qs[0] ^ qs[1] ^ qs[2] ^ qs[3]), "V_0 degenerate"
    return min(tuple(sorted(A[q] for q in qs)) for A in ACT)


def col(th, w):
    return sum(((th >> (4 * s + w)) & 1) << s for s in range(4))


def i1_basis(psis):
    I4 = sum(1 << (5 * i) for i in range(4))
    th = [I4] + psis
    return [sum(col(th[l], w) << (4 * l) for l in range(4)) for w in range(4)]


def read_rows(f):
    return [[int(x) for x in l.split()] for l in open(f) if l.strip()]


def parse_canon_frames(f):
    out = {}
    for l in open(f):
        m = re.match(r"class=(\d+) orbit=(\d+) nodes=([\d,]+) nsol=(\d+)", l)
        if m:
            out[int(m.group(1))] = (int(m.group(2)), [int(x) for x in m.group(3).split(",")], int(m.group(4)))
    return out


def parse_canon_v8(f):
    out = {}
    for l in open(f):
        m = re.match(r"rep=(\d+) orbit=(\d+) V0=([\d,]+) mode=\d+ cand1=\S+ n_V1=(\d+) n_V1_with_V2=\d+ n_V2=(\d+) n_SOL=(\d+)", l)
        if m:
            out[int(m.group(1))] = (int(m.group(2)), [int(x) for x in m.group(3).split(",")], int(m.group(4)), int(m.group(5)), int(m.group(6)))
    return out


def check(a):
    v0m4, r8, r8sp, i1t5, i1m124, i210, i2sp = a
    i1 = read_rows(v0m4)
    inv1 = [invariant(i1_basis(r[1:])) for r in i1]
    assert len(set(inv1)) == len(inv1) == 10, "impl1 representatives not pairwise inequivalent"
    where = {v: c for c, v in enumerate(inv1)}
    t5, m124 = parse_canon_frames(i1t5), parse_canon_frames(i1m124)
    v10, vsp = parse_canon_v8(i210), parse_canon_v8(i2sp)
    lines = ["# impl2 full-gauge representative -> impl1 class, matched by the GL(4,2) x S_4 congruence class of (Q_0..Q_3)",
             "# columns: impl2 rep, impl1 class, orbit (impl2 x 28 vs impl1), #V_1 (impl2 n_V1 vs impl1 nodes[1] at t=5),",
             "#          level 2 (impl2 n_V2 vs impl1 nodes[2] with JMASKS=124), solutions (impl2 n_SOL, impl1 nsol at t=5)"]
    ok = True
    used = []
    for k, row in enumerate(read_rows(r8)):
        c = where.get(invariant(row[1:]))
        used.append(c)
        if c is None:
            lines.append("i2rep=%d -> NO MATCH" % k); ok = False; continue
        o2, _, nv1, nv2, ns = v10[k]
        o1, n1, s1 = t5[c]
        n2 = m124[c][1]
        good = (o2 * 28 == o1 == i1[c][0] and nv1 == n1[1] and nv2 == n2[2] and ns == 0 and s1 == 0 and n1[2] == 0)
        ok &= good
        lines.append("i2rep=%d -> i1class=%d orbit=%d*28=%d/%d n_V1=%d/%d level2=%d/%d solutions=%d/%d %s" % (
            k, c, o2, o2 * 28, o1, nv1, n1[1], nv2, n2[2], ns, s1, "OK" if good else "MISMATCH"))
    bij = sorted(x for x in used if x is not None) == list(range(10))
    ok &= bij
    lines.append("bijection_between_class_lists=%s" % ("yes" if bij else "NO"))
    lines.append("# impl2 Sp(4,2)-only representatives (72) -> impl1 class (same invariant)")
    per = {c: [] for c in range(10)}
    for k, row in enumerate(read_rows(r8sp)):
        c = where.get(invariant(row[1:]))
        if c is None:
            lines.append("sp=%d -> NO MATCH" % k); ok = False; continue
        o2, _, nv1, nv2, ns = vsp[k]
        good = nv1 == t5[c][1][1] and nv2 == m124[c][1][2] and ns == 0
        ok &= good
        per[c].append(o2)
        lines.append("sp=%d -> i1class=%d orbit=%d n_V1=%d n_V2=%d n_SOL=%d %s" % (k, c, o2, nv1, nv2, ns, "OK" if good else "MISMATCH"))
    for c in range(10):
        good = sum(per[c]) * 28 == i1[c][0]
        ok &= good
        lines.append("i1class=%d sp_reps=%d sp_orbit_sum=%d*28=%d/%d %s" % (c, len(per[c]), sum(per[c]), sum(per[c]) * 28, i1[c][0], "OK" if good else "MISMATCH"))
    lines.append("ALL_CONSISTENT=%s" % ("yes" if ok else "NO"))
    return lines


def bits4(x):
    return "".join(str((x >> s) & 1) for s in range(4))


def table(a):
    v0m4, r8 = a
    i1 = read_rows(v0m4)
    inv1 = [invariant(i1_basis(r[1:])) for r in i1]
    where = {v: c for c, v in enumerate(inv1)}
    i2 = read_rows(r8)
    to2 = {where[invariant(r[1:])]: k for k, r in enumerate(i2)}
    out = ["Vectors: 4 nibbles W_0|W_1|W_2|W_3, each written as its coordinates x1 x2 z1 z2 (= bits 0..3).",
           "Matrices psi: rows r = 0..3 written as entries (r,0)(r,1)(r,2)(r,3), separated by '/'; packed integer = sum entry(r,c) 2^(4r+c).",
           ""]
    for c, r in enumerate(i1):
        k = to2[c]
        mats = ["/".join(bits4((p >> (4 * rr)) & 15) for rr in range(4)) for p in r[1:]]
        b1 = i1_basis(r[1:])
        out.append("impl1 class %d (orbit %d tuples = %d x 28; Sp(4,2)^4 x| S_4 stabiliser order %d)  <->  impl2 rep %d (orbit %d triples)" % (
            c, r[0], r[0] // 28, 20160 * 24 // r[0], k, i2[k][0]))
        out.append("  impl1 psi_1, psi_2, psi_3 = %s  (packed %s)" % ("  ".join(mats), " ".join(map(str, r[1:]))))
        out.append("  impl1 V_0 basis (w = e_0..e_3):  %s" % "  ".join("|".join(bits4((b >> (4 * l)) & 15) for l in range(4)) for b in b1))
        out.append("  impl2 V_0 basis:                 %s  (packed %s)" % (
            "  ".join("|".join(bits4((b >> (4 * l)) & 15) for l in range(4)) for b in i2[k][1:]), " ".join(map(str, i2[k][1:]))))
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    res = check(sys.argv[2:]) if mode == "check" else table(sys.argv[2:])
    sys.stdout.write("\n".join(res) + "\n")
    if mode == "check" and res[-1] != "ALL_CONSISTENT=yes":
        sys.exit(1)
