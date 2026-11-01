"""Auditor phase 2: independent count of labelled stabilizer AME(4,4) and AME(6,4) in the V-picture of @classify #315.

S = {(u, M u)}: A = parties 0..h-1, B = parties h..2h-1, M symplectic; V_j = M(A_j) are h mutually orthogonal
nondegenerate 4-dim subspaces with B = sum V_j; AME <=> V_J cap W_L = 0 for |J|+|L| = h.
Number of labelled stabilizer AME(2h,4) states (= number of valid M, = number of graph Lagrangians)
  N(2h) = 720^h * #{ordered valid (V_0..V_{h-1})}.
n=4 control: N(4) must equal 4478976000 (auditor ame44_classes.py, |G|/1440).
Run: python3 vcount.py   (compiles vcount.c with cc)

RELEASE NOTE (release/frames): copied from work/auditor/vcount.py; the only change is that the compiled counter is passed as an
argument (run.sh compiles vcount.c into build/) instead of being compiled next to this file.
usage: python3 vcount.py VCOUNT_BINARY
"""
import itertools, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def sw(v):
    return ((v >> 2) & 0x333) | ((v & 0x333) << 2)


def om(u, v):
    return bin(u & sw(v)).count("1") & 1


def rank_rows(rows, nbits):
    rows = list(rows)
    rk = 0
    for bit in range(nbits):
        p = None
        for i in range(rk, len(rows)):
            if rows[i] >> bit & 1:
                p = i
                break
        if p is None:
            continue
        rows[rk], rows[p] = rows[p], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i] >> bit & 1:
                rows[i] ^= rows[rk]
        rk += 1
    return rk


def nonsing(rows):
    return rank_rows(rows, 4) == 4 and len(rows) == 4


def qform(cols):
    """alternating form Q(e_i,e_j) = omega(phi e_i, phi e_j) as 4 row-nibbles (bit j of row i)"""
    q = [0] * 4
    for i in range(4):
        for j in range(4):
            if om(cols[i], cols[j]):
                q[i] |= 1 << j
    return tuple(q)


def main():
    # GL(4,2): columns = images of e_i
    gl = []
    for cols in itertools.product(range(1, 16), repeat=4):
        if nonsing(list(cols)):
            gl.append(cols)
    assert len(gl) == 20160, len(gl)
    J = qform((1, 2, 4, 8))
    forms = {}
    for cols in gl:
        forms.setdefault(qform(cols), []).append(cols)
    assert len(forms) == 28 and all(len(v) == 720 for v in forms.values())
    sp = forms[J]  # phi with phi^T J phi = J; the Sp(4,2) elements
    print("GL(4,2) 20160, nondegenerate alternating forms 28, Sp(4,2)", len(sp))

    def addq(*qs):
        return tuple(a ^ b ^ (c if False else 0) for a, b in zip(qs[0], qs[1])) if len(qs) == 2 else tuple(
            qs[0][i] ^ qs[1][i] ^ qs[2][i] for i in range(4))

    # ---- n = 4 control
    nvalid = 0
    for Q, phis in forms.items():
        phi = phis[0]
        S = tuple(J[i] ^ Q[i] for i in range(4))
        if not nonsing(list(S)):
            continue
        # V_0 = {(x, phi x)} basis: e_i + (col_i << 4); V_1 = V_0^perp inside 8 bits
        v0 = [(1 << i) | (phi[i] << 4) for i in range(4)]
        perp = [x for x in range(256) if all(om(x, b) == 0 for b in v0)]
        assert len(perp) == 16
        # basis of perp
        basis, span = [], {0}
        for x in perp:
            if x not in span:
                basis.append(x)
                span |= {s ^ x for s in span}
        assert len(basis) == 4
        ok = all(nonsing([(b >> (4 * p)) & 15 for b in basis]) for p in range(2))
        if ok:
            nvalid += 1
    T4 = 720 * nvalid
    print("n=4: valid Q forms %d of 28; ordered (V0,V1) pairs %d; N(4) = 720^2 * T4 = %d (control 4478976000)" % (
        nvalid, T4, 720 ** 2 * T4))

    # ---- n = 6
    def gact(g_cols, Q):
        """Q -> g^T Q g : (g^T Q g)(e_i,e_j) = Q(g e_i, g e_j)"""
        def bil(a, b):
            qb = 0
            for i in range(4):
                if bin(Q[i] & b).count("1") & 1:
                    qb |= 1 << i
            return bin(a & qb).count("1") & 1
        out = [0] * 4
        for i in range(4):
            for j in range(4):
                if bil(g_cols[i], g_cols[j]):
                    out[i] |= 1 << j
        return tuple(out)

    flist = list(forms.keys())
    pairs = []
    for Q4 in flist:
        for Q5 in flist:
            S = tuple(J[i] ^ Q4[i] ^ Q5[i] for i in range(4))
            if nonsing(list(S)):
                pairs.append((Q4, Q5))
    print("n=6: (Q4,Q5) with V_0 nondegenerate:", len(pairs), "of", 28 * 28)
    pset = set(pairs)
    seen, orbits = set(), []
    for pr in pairs:
        if pr in seen:
            continue
        orb = {(gact(g, pr[0]), gact(g, pr[1])) for g in sp}
        assert orb <= pset
        seen |= orb
        orbits.append((pr, len(orb)))
    print("Sp(4,2)-orbits on admissible pairs:", len(orbits), "sizes", sorted(o[1] for o in orbits))
    # build V_0 and R for each representative
    lines = [str(len(orbits))]
    for (Q4, Q5), _ in orbits:
        phi4, phi5 = forms[Q4][0], forms[Q5][0]
        v0 = [(1 << i) | (phi4[i] << 4) | (phi5[i] << 8) for i in range(4)]
        perp = [x for x in range(4096) if all(om(x, b) == 0 for b in v0)]
        assert len(perp) == 256
        basis, span = [], {0}
        for x in perp:
            if x not in span:
                basis.append(x)
                span |= {s ^ x for s in span}
        assert len(basis) == 8
        lines.append(" ".join(map(str, basis)))
    cexe = sys.argv[1]
    res = subprocess.run([cexe], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    ts = [tuple(map(int, l.split())) for l in res.stdout.strip().splitlines()]
    assert all(tot == 200787 for _, tot in ts), ts[:3]
    print("t per orbit rep (valid V_1) :", [t for t, _ in ts])
    S_sum = sum(orb * t for (_, orb), (t, _) in zip(orbits, ts))
    T6 = 720 ** 2 * S_sum
    N6 = 720 ** 3 * T6
    print("sum over pairs of t = %d ; ordered triples T(6) = 720^2 * that = %d ; N(6) labelled stabilizer AME(6,4) = %d = 2^%.3f" % (
        S_sum, T6, N6, __import__("math").log2(N6) if N6 else 0))
    g6 = 720 ** 6
    print("N(6)/720^6 = %.4f  (number of Sp(4,2)^6 gauge classes if all stabilizers trivial; fewer when Aut nontrivial)" % (N6 / g6))
    print("N(6)/(720^6 * 6!) = %.5f" % (N6 / (g6 * 720)))


main()
