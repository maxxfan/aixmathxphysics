"""c5twist.py: exhaustive search of F_2-linear stabilizer AME(8,4) with an automorphism Omega = (5-cycle on parties 0..4) with
order-5 fixed-point-free local twists on parties 5,6 and trivial twist on party 7 (the only 5-cycle-with-twists sector not
excluded by counting, see runs/pairs/NOTES.md shift 5).

Structure (derived): S = S^Omega + S_nf, S^Omega = {(c,c,c,c,c ; 0,0 ; c)} (gauge g=1 on party 7),
S_nf = F_16-Hermitian Lagrangian L of V_nf = F_16^4 (cycle part, hyperbolic Hermitian form H_c) + F_16 + F_16 (parties 5,6, norm forms),
L = graph over A = w^perp (w isotropic point) of an isometry A/<w> -> F_16^2 (Phi=(phi5,phi6)).
Enumerates 1105 isotropic w x (isometries), tests exact party weights of all 2^16 words (prune: c=0 words first).
Usage: python c5twist.py [MAXW] [NSAMPLE_VALID]
"""
import os, sys
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np
# [release] was ROOT = <absolute repo path>, sys.path.insert(0, ROOT + "/qlib"); ame.py now lives in ../lib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

# ---- GF(16), x^4+x+1, g = 2
EXP = [0] * 30
LOG = [0] * 16
x = 1
for i in range(15):
    EXP[i] = x
    LOG[x] = i
    x <<= 1
    if x & 16:
        x ^= 0b10011
for i in range(15, 30):
    EXP[i] = EXP[i - 15]
MT = np.zeros((16, 16), dtype=np.uint8)
for a in range(16):
    for b in range(16):
        MT[a, b] = 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]


def mul(a, b):
    return MT[a, b]


def pw(a, k):
    return 0 if a == 0 else EXP[(LOG[a] * k) % 15]


CONJ = np.array([pw(a, 4) for a in range(16)], dtype=np.uint8)
TR = np.array([a ^ pw(a, 2) ^ pw(a, 4) ^ pw(a, 8) for a in range(16)], dtype=np.uint8)
assert set(TR.tolist()) == {0, 1}
ZETA = EXP[3]
ZP = [1]
for i in range(1, 5):
    ZP.append(int(MT[ZP[-1], ZETA]))
assert int(MT[ZP[4], ZETA]) == 1

# ---- all a in F_16^4 as index a1 + 16 a2 + 256 a3 + 4096 a4
IDX = np.arange(65536)
A4 = np.stack([IDX & 15, (IDX >> 4) & 15, (IDX >> 8) & 15, (IDX >> 12) & 15], axis=1).astype(np.uint8)
# x_i(a) = bits (Tr(a_m zeta^i))_m  as int with bit m
XI = np.zeros((65536, 5), dtype=np.uint8)
for i in range(5):
    v = np.zeros(65536, dtype=np.uint8)
    for m in range(4):
        v |= (TR[MT[A4[:, m], ZP[i]]].astype(np.uint8) << m)
    XI[:, i] = v
WTCYC = (XI != 0).sum(axis=1).astype(np.uint8)
# wtc[a][c] = #{i : x_i(a) != c}
WTC = np.zeros((65536, 16), dtype=np.uint8)
for c in range(16):
    WTC[:, c] = (XI != c).sum(axis=1)


def Hc(a, b):
    """Hermitian form on F_16^4 (arrays (...,4))."""
    r = MT[a[..., 0], CONJ[b[..., 1]]] ^ MT[a[..., 1], CONJ[b[..., 0]]] ^ MT[a[..., 2], CONJ[b[..., 3]]] ^ MT[a[..., 3], CONJ[b[..., 2]]]
    return r


def aindex(a):
    return a[..., 0].astype(np.int64) + 16 * a[..., 1].astype(np.int64) + 256 * a[..., 2].astype(np.int64) + 4096 * a[..., 3].astype(np.int64)


def check_form():
    """sum_i B4(x_i(a), x_i(b)) == Tr(Hc(a,b)) for random a,b"""
    rng = np.random.default_rng(1)
    for _ in range(2000):
        ia, ib = rng.integers(0, 65536, 2)
        s = 0
        for i in range(5):
            u, v = int(XI[ia, i]), int(XI[ib, i])
            u1, u2, u3, u4 = [(u >> k) & 1 for k in range(4)]
            v1, v2, v3, v4 = [(v >> k) & 1 for k in range(4)]
            s ^= (u1 & v2) ^ (u2 & v1) ^ (u3 & v4) ^ (u4 & v3)
        assert s == TR[Hc(A4[ia], A4[ib])], (ia, ib)
    # zeta action: x_i(zeta a) = x_{i+1}(a)
    for _ in range(200):
        ia = int(rng.integers(0, 65536))
        za = MT[A4[ia], ZETA]
        ja = int(aindex(za))
        assert [int(XI[ja, i]) for i in range(5)] == [int(XI[ia, (i + 1) % 5]) for i in range(5)]


# ---- coordinates of F_16 elements for a party (basis with Tr(b_k conj(b_l)) = J)
def party_basis():
    J = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
    import itertools
    for b in itertools.permutations(range(1, 16), 4):
        ok = True
        for k in range(4):
            for l in range(k, 4):
                if int(TR[MT[b[k], CONJ[b[l]]]]) != J[k][l]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return b
    raise RuntimeError


PB = party_basis()
COORD = {}
for t in range(16):
    u = 0
    for k in range(4):
        if (t >> k) & 1:
            u ^= PB[k]
    COORD[u] = t
assert len(COORD) == 16  # basis => bijection


def rows_from(gens):
    """gens: list of tuples (a (4 F16 elems) , u5, u6, c) -> F_2 rows [X(16)|Z(16)] (party j qubits 2j,2j+1; local coords t1..t4 = X0 Z0 X1 Z1)"""
    rows = []
    for (a, u5, u6, c) in gens:
        ia = int(aindex(np.array(a, dtype=np.uint8)))
        locs = [int(XI[ia, i]) ^ c_shift for i, c_shift in zip(range(5), [0] * 5)]
        # S^Omega part is added separately by caller via c (c^5 on cycle, c on party 7)
        locs = [int(XI[ia, i]) for i in range(5)] + [COORD[u5], COORD[u6], 0]
        rows.append(locs)
    return rows


def locs_to_row(locs):
    X = [0] * 16
    Z = [0] * 16
    for j, t in enumerate(locs):
        t1, t2, t3, t4 = [(t >> k) & 1 for k in range(4)]
        X[2 * j], Z[2 * j], X[2 * j + 1], Z[2 * j + 1] = t1, t2, t3, t4
    return X + Z


def build_S_rows(wvec, a1, a2, psi1, psi2):
    """12 F_2-gens of L (F_16-basis (w;0,0),(a1;psi1),(a2;psi2) times basis 1,g,g^2,g^3) + 4 gens of S^Omega."""
    gens = []
    base = [(wvec, 0, 0), (a1, psi1[0], psi1[1]), (a2, psi2[0], psi2[1])]
    rows = []
    for (a, u5, u6) in base:
        for k in range(4):
            xi = EXP[k]
            aa = [int(MT[xi, int(t)]) for t in a]
            ia = int(aindex(np.array(aa, dtype=np.uint8)))
            locs = [int(XI[ia, i]) for i in range(5)] + [COORD[int(MT[xi, u5])], COORD[int(MT[xi, u6])], 0]
            rows.append(locs_to_row(locs))
    for k in range(4):
        c = 1 << k
        rows.append(locs_to_row([c] * 5 + [0, 0, c]))
    return rows


def gf2_rank_isotropic(rows):
    import numpy as _np
    M = _np.array(rows, dtype=_np.uint8)
    N = 16
    X, Z = M[:, :N], M[:, N:]
    G = (X @ Z.T + Z @ X.T) % 2
    iso = not G.any()
    # rank
    A = [int("".join(map(str, r)), 2) for r in rows]
    rk = 0
    basis = []
    for v in A:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            rk += 1
    return iso, rk


def main():
    import time
    t0 = time.time()
    maxw = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
    nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    check_form()
    print("form checks ok; party basis", PB, flush=True)
    # isotropic projective points
    hw = Hc(A4, A4)
    iso_pts = []
    seen = set()
    for i in range(1, 65536):
        if hw[i] != 0:
            continue
        # normal form: scale so that first nonzero coord = 1
        a = A4[i]
        k = int(np.nonzero(a)[0][0])
        inv = EXP[(15 - LOG[int(a[k])]) % 15]
        s = MT[inv, a]
        j = int(aindex(s))
        if j == i:
            iso_pts.append(i)
    print("isotropic projective points:", len(iso_pts), "(expect 1105)", flush=True)
    stats = dict(w=0, cand=0, pass_c0=0, pass_all=0)
    best = []
    samples = []
    minw_hist = {}
    # psi enumeration tables
    PSI = np.array([(u, v) for u in range(16) for v in range(16)], dtype=np.uint8)  # 256 x 2
    NORM = MT[PSI[:, 0], CONJ[PSI[:, 0]]] ^ MT[PSI[:, 1], CONJ[PSI[:, 1]]]
    for wi, iw in enumerate(iso_pts):
        if wi >= maxw:
            break
        w = A4[iw]
        # A = w^perp
        hv = Hc(A4, np.broadcast_to(w, A4.shape))
        Aidx = np.nonzero(hv == 0)[0]
        assert len(Aidx) == 4096
        # basis completion
        def span_of(vs):
            out = {tuple([0, 0, 0, 0])}
            for v in vs:
                new = set()
                for t in range(16):
                    tv = tuple(int(MT[t, int(z)]) for z in v)
                    for o in out:
                        new.add(tuple(p ^ q for p, q in zip(o, tv)))
                out = new
            return out
        sp = span_of([w])
        a1 = next(A4[i] for i in Aidx if tuple(int(z) for z in A4[i]) not in sp)
        sp = span_of([w, a1])
        a2 = next(A4[i] for i in Aidx if tuple(int(z) for z in A4[i]) not in sp)
        G11 = int(Hc(a1, a1)); G22 = int(Hc(a2, a2)); G12 = int(Hc(a1, a2))
        # candidate psi1: norm G11 ; psi2: norm G22
        P1 = np.nonzero(NORM == G11)[0]
        P2 = np.nonzero(NORM == G22)[0]
        # psi1 . conj(psi2) = G12
        pr = MT[PSI[P1][:, None, 0], CONJ[PSI[P2][None, :, 0]]] ^ MT[PSI[P1][:, None, 1], CONJ[PSI[P2][None, :, 1]]]
        ii, jj = np.nonzero(pr == G12)
        K = len(ii)
        stats['w'] += 1
        stats['cand'] += K
        if K == 0:
            continue
        psi1 = PSI[P1[ii]]  # K x 2
        psi2 = PSI[P2[jj]]
        # all 4096 coefficient triples (t0,t1,t2)
        T = np.array([(t0, t1, t2) for t0 in range(16) for t1 in range(16) for t2 in range(16)], dtype=np.uint8)
        # a = t0 w + t1 a1 + t2 a2
        av = MT[T[:, 0][:, None], w[None, :]] ^ MT[T[:, 1][:, None], a1[None, :]] ^ MT[T[:, 2][:, None], a2[None, :]]
        ia = aindex(av)
        wcy = WTCYC[ia].astype(np.int16)  # 4096
        # u5 = t1 psi1[0] + t2 psi2[0]
        u5 = MT[T[:, 1][None, :], psi1[:, 0][:, None]] ^ MT[T[:, 2][None, :], psi2[:, 0][:, None]]  # K x 4096
        u6 = MT[T[:, 1][None, :], psi1[:, 1][:, None]] ^ MT[T[:, 2][None, :], psi2[:, 1][:, None]]
        wt0 = wcy[None, :] + (u5 != 0) + (u6 != 0)
        nz = (T.sum(axis=1) != 0)
        wt0 = np.where(nz[None, :], wt0, 99)
        mn0 = wt0.min(axis=1)
        for v in mn0:
            minw_hist[int(v)] = minw_hist.get(int(v), 0) + 1
        ok0 = np.nonzero(mn0 >= 5)[0]
        stats['pass_c0'] += len(ok0)
        for k in ok0:
            nzu = (u5[k] != 0).astype(np.int16) + (u6[k] != 0)
            wfull = WTC[ia][:, 1:] + nzu[:, None] + 1  # 4096 x 15 (c != 0)
            mnf = int(wfull.min())
            if mnf >= 5:
                stats['pass_all'] += 1
                best.append((wi, int(k)))
                print("SOLUTION candidate", wi, k, flush=True)
                rows = build_S_rows([int(z) for z in w], [int(z) for z in a1], [int(z) for z in a2], [int(z) for z in psi1[k]], [int(z) for z in psi2[k]])
                import ame
                print(gf2_rank_isotropic(rows), ame.bad_cuts(rows, 2, 2, 8), flush=True)
            else:
                minw_hist[('full', mnf)] = minw_hist.get(('full', mnf), 0) + 1
                if len(samples) < nsamp:
                    samples.append((wi, int(k), [int(z) for z in w], [int(z) for z in a1], [int(z) for z in a2], [int(z) for z in psi1[k]], [int(z) for z in psi2[k]], mnf))
        if wi % 100 == 0:
            print("w", wi, stats, round(time.time() - t0), "s", flush=True)
    print("DONE", stats, "hist", sorted(minw_hist.items(), key=lambda kv: str(kv[0])), round(time.time() - t0), "s", flush=True)
    # validation of the formula on a few candidates (those passing c=0 test but failing overall, plus no pass)
    import ame
    for (wi, k, w, a1, a2, p1, p2, mnf) in samples:
        rows = build_S_rows(w, a1, a2, p1, p2)
        iso, rk = gf2_rank_isotropic(rows)
        # direct min weight over 2^16 words
        R = np.array(rows, dtype=np.uint8)
        words = np.zeros((1, 32), dtype=np.uint8)
        allw = np.zeros((1, 32), dtype=np.uint8)
        for r in R:
            allw = np.concatenate([allw, allw ^ r[None, :]], axis=0)
        X, Z = allw[:, :16], allw[:, 16:]
        party = ((X[:, 0::2] | X[:, 1::2] | Z[:, 0::2] | Z[:, 1::2]) != 0).sum(axis=1)
        pm = int(party[1:].min())
        print("VALIDATE cand", wi, k, "isotropic", iso, "rank", rk, "direct min party weight", pm, "formula", mnf, "bad_cuts", ame.bad_cuts(rows, 2, 2, 8), flush=True)


if __name__ == "__main__":
    main()
