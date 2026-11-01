"""Auditor phase 2: classes of V_0 for the n=8 V-picture (classify #315), independent count.

V_0 = {(x, phi5 x, phi6 x, phi7 x)} in B = W4+W5+W6+W7 with phi_l in GL(4,2); only Q_l = phi_l^T J phi_l matters (28 forms each),
Q_4 = J.  Admissible: J+Q5+Q6+Q7 nondegenerate (V_0 nondegenerate).  Equivalence: B-gauge Sp(4,2)^4 (already quotiented by using Q),
simultaneous GL(4,2)-congruence fixing Q_4 = J (i.e. Sp(4,2)), S_3 on parties 5,6,7, and the swap of the base party 4 <-> 5 (renormalised).
Output: number of admissible triples, number of orbits under Sp only and under Sp x S_4 (full party permutations of W), with orbit sizes.

RELEASE NOTE (release/frames, 2026-10-31): this is a standard-library port of work/auditor/vclasses8.py, which used numpy only
for vectorised bookkeeping (action tables, orbit labelling by iterated minimum propagation).  The mathematics, the group
generators (all 720 Sp(4,2) elements, the 6 permutations of parties 5,6,7, the base swap 4 <-> 5), the choice of the form
representatives phi (first GL(4,2) element in itertools.product order with that form), the class representative (the smallest
triple code (a5*28+a6)*28+a7 in the orbit) and the output order (ascending representative code) are unchanged; orbits are
computed as connected components by breadth-first search instead of numpy minimum propagation (same partition, same minima).
The output files are byte-identical to the original v0reps8.txt / v0reps8_sp.txt (checked by run.sh).
Run: python3 vclasses8.py [OUTDIR]     (writes OUTDIR/v0reps8.txt and OUTDIR/v0reps8_sp.txt; default OUTDIR = current directory)
"""
import itertools
import os
import sys
from collections import deque


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
    return len(rows) == 4 and rank_rows(rows, 4) == 4


def qform(cols):
    q = [0] * 4
    for i in range(4):
        for j in range(4):
            if om(cols[i], cols[j]):
                q[i] |= 1 << j
    return tuple(q)


def apply(cols, x):
    y = 0
    for i in range(4):
        if x >> i & 1:
            y ^= cols[i]
    return y


def inverse(cols):
    inv = {apply(cols, x): x for x in range(16)}
    return tuple(inv[1 << i] for i in range(4))


def bil(Q, a, b):
    qb = 0
    for i in range(4):
        if bin(Q[i] & b).count("1") & 1:
            qb |= 1 << i
    return bin(a & qb).count("1") & 1


def gact(cols, Q):
    """Q(g e_i, g e_j) = g^T Q g"""
    out = [0] * 4
    for i in range(4):
        for j in range(4):
            if bil(Q, cols[i], cols[j]):
                out[i] |= 1 << j
    return tuple(out)


def main(outdir):
    gl = [c for c in itertools.product(range(1, 16), repeat=4) if nonsing(list(c))]
    assert len(gl) == 20160
    J = qform((1, 2, 4, 8))
    forms = {}
    for c in gl:
        forms.setdefault(qform(c), []).append(c)
    flist = sorted(forms)
    idx = {Q: i for i, Q in enumerate(flist)}
    assert len(flist) == 28
    sp = forms[J]
    assert len(sp) == 720
    # action tables
    act_sp = [[idx[gact(g, Q)] for Q in flist] for g in sp]  # 720 x 28
    phi = {Q: forms[Q][0] for Q in flist}
    # renormalisation maps for the swap move: psi = phi_Q^{-1}; act(psi, Q')
    ren = [[0] * 28 for _ in range(28)]  # ren[a][b] = idx of act(phi_a^{-1}, Q_b)
    for a, Qa in enumerate(flist):
        psi = inverse(phi[Qa])
        for b, Qb in enumerate(flist):
            ren[a][b] = idx[gact(psi, Qb)]
    iJ = idx[J]

    # triples (a5,a6,a7) in 28^3, code = (a5*28 + a6)*28 + a7
    N = 28 ** 3

    def code(a, b, c):
        return (a * 28 + b) * 28 + c

    def dec(t):
        return t // 784, (t // 28) % 28, t % 28

    adm = [False] * N
    for t in range(N):
        a5, a6, a7 = dec(t)
        s = tuple(J[i] ^ flist[a5][i] ^ flist[a6][i] ^ flist[a7][i] for i in range(4))
        adm[t] = nonsing(list(s))
    print("admissible triples:", sum(adm), "of", N)

    perms = list(itertools.permutations(range(3)))

    def sp_image(g, t):
        A = act_sp[g]
        a5, a6, a7 = dec(t)
        return code(A[a5], A[a6], A[a7])

    def perm_image(p, t):
        cols_ = dec(t)
        return code(cols_[p[0]], cols_[p[1]], cols_[p[2]])

    def swap_image(t):
        # base party 4 <-> 5:  new triple = (act(psi,J), act(psi,Q6), act(psi,Q7)), psi = phi_{Q5}^{-1}
        a5, a6, a7 = dec(t)
        return code(ren[a5][iJ], ren[a5][a6], ren[a5][a7])

    # sanity: the maps preserve admissibility (same subset of maps as the original: 5 Sp elements, the 6 perms, the swap)
    for t in range(N):
        imgs = [sp_image(g, t) for g in range(5)] + [perm_image(p, t) for p in perms] + [swap_image(t)]
        for u in imgs:
            assert adm[u] == adm[t]

    def neighbours_sp(t):
        return [sp_image(g, t) for g in range(720)]

    def neighbours_all(t):
        return neighbours_sp(t) + [perm_image(p, t) for p in perms] + [swap_image(t)]

    def orbits(neigh):
        """connected components of the admissible triples under the given maps; returns list of (min code, size), ascending.
        Every map is a bijection of a finite set and the generated group contains its inverses, so components = orbits."""
        seen = [False] * N
        out = []
        for t in range(N):  # ascending: the first unseen element of a component is its minimum
            if not adm[t] or seen[t]:
                continue
            seen[t] = True
            q = deque([t])
            size = 0
            while q:
                x = q.popleft()
                size += 1
                for y in neigh(x):
                    if not seen[y]:
                        seen[y] = True
                        q.append(y)
            out.append((t, size))
        return out

    orb_sp = orbits(neighbours_sp)
    orb_all = orbits(neighbours_all)
    for name, orb in (("Sp(4,2) only", orb_sp), ("Sp x S_3 x swap (= Sp(4,2)^4 x S_4 on V_0)", orb_all)):
        print(name, ": orbits", len(orb), "sizes", sorted(c for _, c in orb))

    def write(fname, orb):
        with open(os.path.join(outdir, fname), "w") as fh:
            for t, cnt in orb:
                a5, a6, a7 = dec(t)
                p5, p6, p7 = phi[flist[a5]], phi[flist[a6]], phi[flist[a7]]
                basis = [(1 << i) | (p5[i] << 4) | (p6[i] << 8) | (p7[i] << 12) for i in range(4)]
                fh.write("%d %d %d %d %d\n" % (cnt, *basis))

    # one representative V_0 per full-gauge class: 4 basis vectors, 16-bit (party p = bits 4p..4p+3; W_4 first)
    write("v0reps8.txt", orb_all)
    print("wrote v0reps8.txt (orbit size, 4 basis vectors) for", len(orb_all), "classes")
    # also the 72 classes under Sp(4,2) congruence only (no S_4 / swap assumption)
    write("v0reps8_sp.txt", orb_sp)
    print("wrote v0reps8_sp.txt for", len(orb_sp), "Sp-only classes")


main(sys.argv[1] if len(sys.argv) > 1 else ".")
