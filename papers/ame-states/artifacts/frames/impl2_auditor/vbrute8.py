"""Independent brute-force check of V-quadruples (V_0..V_3) at n=8, no matrix M needed.
S = {(u, Mu)}: element with u = sum u_j (u_j in A_j), Mu = sum m_j, m_j in V_j (m_j != 0 iff u_j != 0, M_j bijective).
Party weight = |J| + wt(sum_j m_j), J = {j : m_j != 0}.  min over the 16^4 - 1 tuples = min party weight of S (any isometries M_j).
AME(8,4) <=> min >= 5.   Reads lines 'SOL V0=a,b,c,d V1=.. V2=..', recomputes V_3 = (V0+V1+V2)^perp.

RELEASE NOTE (release/frames): standard-library port of work/auditor/vbrute8.py (the original used numpy to evaluate the
16^4 grid).  Same formula, same V_3 (the greedy basis of the perp space scanned in increasing order; the perp space is now obtained
from a null-space basis and sorted, which yields the same list of 16 vectors as the original scan of all 65536 x).  Changes:
the solution file is an argument (the original read v8c_solutions.txt next to itself), and one extra output line gives the
histogram of (min party weight, number of words of that weight), i.e. A_4 for the near-misses (argument.md 5.14(c)).
usage: vbrute8.py SOLFILE [N]   (check first N solutions; default 20)
"""
import sys


def sw(v):
    return ((v >> 2) & 0x3333) | ((v & 0x3333) << 2)


def om(u, v):
    return bin(u & sw(v)).count("1") & 1


def span(b):
    out = []
    for m in range(16):
        x = 0
        for i in range(4):
            if m >> i & 1:
                x ^= b[i]
        out.append(x)
    return out


WT = [((x | (x >> 1) | (x >> 2) | (x >> 3)) & 0x1111) for x in range(65536)]
WT = [(n & 1) + ((n >> 4) & 1) + ((n >> 8) & 1) + ((n >> 12) & 1) for n in WT]


def min_party_weight(Vs):
    E = [span(b) for b in Vs]
    best, cnt = 99, 0
    E0, E1, E2, E3 = E
    for a in range(16):
        xa, ja = E0[a], (a != 0)
        for b in range(16):
            xab, jab = xa ^ E1[b], ja + (b != 0)
            for c in range(16):
                xabc, jabc = xab ^ E2[c], jab + (c != 0)
                for d in range(16):
                    if a == 0 and b == 0 and c == 0 and d == 0:
                        continue
                    w = jabc + (d != 0) + WT[xabc ^ E3[d]]
                    if w < best:
                        best, cnt = w, 1
                    elif w == best:
                        cnt += 1
    return best, cnt


def basis_of(vecs):
    basis, sp = [], {0}
    for x in vecs:
        if x not in sp:
            basis.append(x)
            sp |= {s ^ x for s in sp}
    return basis


def perp_space(gens):
    """all x in F_2^16 with om(x, g) = 0 for every g, in increasing order (null space of the functionals x -> <x, sw(g)>)."""
    rows = [sw(g) for g in gens]
    # reduced row echelon form over F_2, pivot = highest bit
    piv = []
    for r in rows:
        for p, pr in piv:
            if r >> p & 1:
                r ^= pr
        if r:
            p = r.bit_length() - 1
            piv = [(q, (qr ^ r) if (qr >> p & 1) else qr) for q, qr in piv]
            piv.append((p, r))
    pivcols = {p for p, _ in piv}
    free = [c for c in range(16) if c not in pivcols]
    null = []
    for f in free:
        x = 1 << f
        for p, pr in piv:
            if pr >> f & 1:
                x |= 1 << p
        null.append(x)
    out = sorted(span(null)) if len(null) == 4 else None
    assert out is not None and all(om(x, g) == 0 for x in out for g in gens)
    return out


def main():
    fname = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    lines = [l for l in open(fname) if l.startswith("SOL")][:n]
    hist, hist2 = {}, {}
    for l in lines:
        parts = l.split()
        V = [list(map(int, p.split("=")[1].split(","))) for p in parts[1:4]]
        allv = [x for b in V for x in span(b)]
        sp3 = basis_of(allv)
        assert len(sp3) == 12, len(sp3)  # V0+V1+V2 direct sum
        perp = perp_space(sp3)
        V3 = basis_of(perp)
        assert len(V3) == 4
        mw, cnt = min_party_weight(V + [V3])
        hist[mw] = hist.get(mw, 0) + 1
        hist2[(mw, cnt)] = hist2.get((mw, cnt), 0) + 1
    print("checked", len(lines), "relaxed solutions; histogram of brute-force min party weight:", hist)
    print("histogram of (min party weight, number of words of that weight):", dict(sorted(hist2.items())))


main()
