"""framecheck.py (release/frames, written for this release; standard library only).

Independent brute-force check of the frames printed by impl1 (frames.c / frames3.c solution files).  It replaces classify's
numpy tools frame2stab.py / pattern_check.py (PROOF.md section 8, items 2-4 and 6c), which are not part of the release.

Input line format (frames.c print_frame): for i = 0..m-1 the matrices th[i][1..m-1] (packed 4x4 over F_2, entry (r,c) = bit 4r+c),
groups separated by '|', optionally followed by 'P: <bits>' (frames3 with PATTERN=1: validity of every condition (J,L) with
|J|+|L| = m, in the loop order jm = 1..2^m-1, lm = 1..2^m-1).  V_i = {(w, th[i][1] w, ..., th[i][m-1] w)}, party l = nibble l.

For every frame it checks from scratch:
  * every th[i][l] invertible, each V_i nondegenerate for omega(u,v) = parity(u & sw(v)) (sw swaps bits 0<->2, 1<->3 of each
    nibble, i.e. J pairs x1<->z1, x2<->z2), V_i orthogonal to V_j for i != j, V_0 + ... + V_{m-1} = F_2^{4m};
  * the party-weight enumerator of the stabilizer S = {(u, Mu)} by brute force over all 16^m tuples (m_0..m_{m-1}), m_j in V_j:
    weight = #{j : m_j != 0} + #{nonzero nibbles of sum m_j}  (argument.md 5.12: independent of the isometries M|A_j);
  * if a P string is present: for each (J,L) the bit must equal 'no nonzero word has A-support in J and B-support in L'.
usage: framecheck.py M SOLFILE
"""
import sys
from collections import Counter


def sw(v):
    m = 0
    for k in range(8):
        m |= 0x3 << (4 * k)
    return ((v >> 2) & m) | ((v & m) << 2)


def om(u, v):
    return bin(u & sw(v)).count("1") & 1


def rank(vecs):
    basis = []
    for x in vecs:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
            basis.sort(reverse=True)
    return len(basis)


def inv4(a):
    rows = [(a >> (4 * r)) & 15 for r in range(4)]
    return rank(rows) == 4


def col(th, w):
    return sum(((th >> (4 * s + w)) & 1) << s for s in range(4))


def span(b):
    out = []
    for msk in range(16):
        x = 0
        for i in range(4):
            if msk >> i & 1:
                x ^= b[i]
        out.append(x)
    return out


def main():
    m = int(sys.argv[1])
    fname = sys.argv[2]
    nb = 4 * m
    NZ = []
    for x in range(1 << nb):
        msk = 0
        for l in range(m):
            if (x >> (4 * l)) & 15:
                msk |= 1 << l
        NZ.append(msk)
    PC = [bin(k).count("1") for k in range(1 << m)]
    I4 = sum(1 << (5 * i) for i in range(4))
    conds = [(jm, lm) for jm in range(1, 1 << m) for lm in range(1, 1 << m) if PC[jm] + PC[lm] == m]
    nframes = nvalid = npat = nagree = 0
    minw = Counter()
    enums = Counter()
    for line in open(fname):
        line = line.strip()
        if not line:
            continue
        P = None
        if "P:" in line:
            line, P = line.split("P:")
            P = P.strip()
        groups = [g.split() for g in line.split("|")]
        groups = [g for g in groups if g]
        assert len(groups) == m and all(len(g) == m - 1 for g in groups), line
        th = [[I4] + [int(x) for x in g] for g in groups]
        nframes += 1
        vec = [[sum(col(th[i][l], w) << (4 * l) for l in range(m)) for w in range(4)] for i in range(m)]
        ok = all(inv4(th[i][l]) for i in range(m) for l in range(m))
        for i in range(m):
            gram = [sum(om(vec[i][a], vec[i][b]) << b for b in range(4)) for a in range(4)]
            ok = ok and rank(gram) == 4
            for j in range(i + 1, m):
                ok = ok and all(om(a, b) == 0 for a in vec[i] for b in vec[j])
        ok = ok and rank([v for vi in vec for v in vi]) == nb
        nvalid += ok
        E = [span(vi) for vi in vec]
        hist = [0] * (2 * m + 1)
        pairs = set()
        # all 16^m tuples, iteratively: (partial sum, A-support mask)
        cur = [(0, 0)]
        for i in range(m):
            Ei = E[i]
            cur = [(s ^ Ei[k], a | ((k != 0) << i)) for (s, a) in cur for k in range(16)]
        for s, a in cur:
            b = NZ[s]
            hist[PC[a] + PC[b]] += 1
            if a:
                pairs.add((a, b))
        assert hist[0] == 1 and sum(hist) == 1 << nb
        mw = next(w for w in range(1, 2 * m + 1) if hist[w])
        minw[mw] += 1
        enums[",".join(map(str, hist))] += 1
        if P is not None:
            npat += 1
            assert len(P) == len(conds), (len(P), len(conds))
            mine = "".join("0" if any((a & ~jm) == 0 and (b & ~lm) == 0 for (a, b) in pairs) else "1" for jm, lm in conds)
            nagree += mine == P
    print("frames=%d m=%d" % (nframes, m))
    print("valid_frames=%d (theta invertible, V_i nondegenerate, pairwise orthogonal, direct sum)" % nvalid)
    print("min_party_weight_histogram=%s" % dict(sorted(minw.items())))
    for e, c in sorted(enums.items()):
        print("party_weight_enumerator=[%s] frames=%d" % (e, c))
    print("pattern_lines=%d pattern_agree=%d (%d conditions each)" % (npat, nagree, len(conds)))


main()
