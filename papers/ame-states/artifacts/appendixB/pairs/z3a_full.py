"""z3a_full.py: stages 2-3 of the sigma=(012)(345), fpf order-3 twists on parties 6,7 sector (auditor's case (A)).
Reads work/pairs/z3a_reps.json (3 orbit reps of T_A from z3a_stage1.py orbits). For each rep: enumerate all T_B (M'=(M^dagger)^-1, b from Gram/weight),
T = T_A + T_B (F_4-dim 4), T^perp/T (F_4^4), good isotropic cosets (all 256 coset elements of weight >= 5), isotropic good lines = completions L of S_nf,
then S^Omega = {(c^3,(gc)^3)} for all g in Sp(4,2): S = S^Omega + L must have min party weight >= 5.
Usage: python z3a_full.py [REP_INDEX|all] [MAX_TB]
"""
import os, sys, json, itertools, time
os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MT = np.zeros((4, 4), dtype=np.uint8)
EXP4 = [1, 2, 3]
LOG4 = {1: 0, 2: 1, 3: 2}
for a in range(1, 4):
    for b in range(1, 4):
        MT[a, b] = EXP4[(LOG4[a] + LOG4[b]) % 3]
CONJ = np.array([0, 1, 3, 2], dtype=np.uint8)
TR = np.array([a ^ int(CONJ[a]) for a in range(4)], dtype=np.uint8)
W = [1, 2, 3]
IDX = np.arange(256)
A4 = np.stack([IDX & 3, (IDX >> 2) & 3, (IDX >> 4) & 3, (IDX >> 6) & 3], axis=1).astype(np.uint8)
Y = np.zeros((256, 3), dtype=np.uint8)
for i in range(3):
    v = np.zeros(256, dtype=np.uint8)
    for m in range(4):
        v |= (TR[MT[A4[:, m], W[i]]].astype(np.uint8) << m)
    Y[:, i] = v
WT = (Y != 0).sum(axis=1).astype(np.int16)
WTC = np.zeros((256, 16), dtype=np.int16)  # #{i : y_i(a) != c}
for c in range(16):
    WTC[:, c] = (Y != c).sum(axis=1)


def aidx(a):
    a = np.asarray(a)
    return a[..., 0].astype(np.int64) + 4 * a[..., 1].astype(np.int64) + 16 * a[..., 2].astype(np.int64) + 64 * a[..., 3].astype(np.int64)


def H4(a, b):
    return MT[a[..., 0], CONJ[b[..., 1]]] ^ MT[a[..., 1], CONJ[b[..., 0]]] ^ MT[a[..., 2], CONJ[b[..., 3]]] ^ MT[a[..., 3], CONJ[b[..., 2]]]


def inner2(u, v):
    return MT[u[..., 0], CONJ[v[..., 0]]] ^ MT[u[..., 1], CONJ[v[..., 1]]]


def H12(u, v):
    return H4(u[..., 0:4], v[..., 0:4]) ^ H4(u[..., 4:8], v[..., 4:8]) ^ inner2(u[..., 8:10], v[..., 8:10]) ^ inner2(u[..., 10:12], v[..., 10:12])


def WT12(v):
    return WT[aidx(v[..., 0:4])] + WT[aidx(v[..., 4:8])] + (v[..., 8:10].any(axis=-1)) + (v[..., 10:12].any(axis=-1))


def inner(u, v):
    return int(MT[u[0], CONJ[v[0]]] ^ MT[u[1], CONJ[v[1]]])


def mmul(P, Q):
    return (int(MT[P[0], Q[0]] ^ MT[P[1], Q[2]]), int(MT[P[0], Q[1]] ^ MT[P[1], Q[3]]), int(MT[P[2], Q[0]] ^ MT[P[3], Q[2]]), int(MT[P[2], Q[1]] ^ MT[P[3], Q[3]]))


def minv(P):
    det = int(MT[P[0], P[3]] ^ MT[P[1], P[2]])
    di = EXP4[(-LOG4[det]) % 3]
    return (int(MT[di, P[3]]), int(MT[di, P[1]]), int(MT[di, P[2]]), int(MT[di, P[0]]))


def mdagger(P):
    return (int(CONJ[P[0]]), int(CONJ[P[2]]), int(CONJ[P[1]]), int(CONJ[P[3]]))


def cols(M):
    return ((M[0], M[2]), (M[1], M[3]))


def enumA(M):
    """all (alpha, beta) index pairs with Gram H(a_i,a_j) = delta_ij + <M e_i, M e_j> and cycle weight 3 on all nonzero a(f)."""
    Me1, Me2 = cols(M)
    G11 = 1 ^ inner(Me1, Me1)
    G22 = 1 ^ inner(Me2, Me2)
    G12 = inner(Me1, Me2)
    hh = H4(A4, A4)
    c1 = np.nonzero((WT == 3) & (hh == G11))[0]
    c2 = np.nonzero((WT == 3) & (hh == G22))[0]
    pr = H4(A4[c1][:, None, :], A4[c2][None, :, :])
    ii, jj = np.nonzero(pr == G12)
    out = []
    for i, j in zip(ii, jj):
        al, be = A4[c1[i]], A4[c2[j]]
        ok = True
        for (f1, f2) in [(1, 1), (1, 2), (1, 3)]:
            v = MT[f1, al] ^ MT[f2, be]
            if WT[aidx(v)] != 3:
                ok = False
                break
        if ok:
            out.append((int(c1[i]), int(c2[j])))
    return out


# Sp(4,2) as F_2 matrices; action on 4-bit c: bit m of (gc) = sum_n g[m,n] c_n
Jm = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
SP = []
for bits in range(1 << 16):
    g = np.array([(bits >> k) & 1 for k in range(16)]).reshape(4, 4)
    if ((g.T @ Jm @ g) % 2 == Jm).all():
        SP.append(g)
assert len(SP) == 720
GC = np.zeros((720, 16), dtype=np.int64)  # (g c) as 4-bit int
for k, g in enumerate(SP):
    for c in range(16):
        cv = np.array([(c >> m) & 1 for m in range(4)])
        d = (g @ cv) % 2
        GC[k, c] = sum(int(d[m]) << m for m in range(4))


def rref(vecs):
    """RREF over F_4 of list of 12-vectors -> (rows, pivots)."""
    rows = [np.array(v, dtype=np.uint8) for v in vecs]
    piv = []
    r = 0
    for c in range(12):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = EXP4[(-LOG4[int(rows[r][c])]) % 3]
        rows[r] = MT[inv, rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                rows[i] = rows[i] ^ MT[rows[i][c], rows[r]]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows[:r], piv


def perp_of(T):
    """all vectors v of F_4^12 of the T^perp form, T = [tA1,tA2,tB1,tB2]; return array (N,12). Uses the join over (f6,f7)."""
    # enumerate keys (f6,f7) in F_4^4; a-part constraint from T_A (first two), b-part from T_B (last two)
    fa = np.array(list(itertools.product(range(4), repeat=4)), dtype=np.uint8)  # (256,4) = (f6_0,f6_1,f7_0,f7_1)
    out = []
    for key in fa:
        # a candidates: H12(t, v) = 0 for t in T_A with v = (a,0,key)
        va = np.zeros((256, 12), dtype=np.uint8)
        va[:, 0:4] = A4
        va[:, 8:12] = key
        ma = np.ones(256, dtype=bool)
        for t in T[:2]:
            ma &= (H12(np.broadcast_to(t, va.shape), va) == 0)
        vb = np.zeros((256, 12), dtype=np.uint8)
        vb[:, 4:8] = A4
        vb[:, 8:12] = key
        mb = np.ones(256, dtype=bool)
        for t in T[2:]:
            mb &= (H12(np.broadcast_to(t, vb.shape), vb) == 0)
        ai = np.nonzero(ma)[0]
        bi = np.nonzero(mb)[0]
        if len(ai) == 0 or len(bi) == 0:
            continue
        blk = np.zeros((len(ai) * len(bi), 12), dtype=np.uint8)
        blk[:, 0:4] = np.repeat(A4[ai], len(bi), axis=0)
        blk[:, 4:8] = np.tile(A4[bi], (len(ai), 1))
        blk[:, 8:12] = key
        out.append(blk)
    return np.concatenate(out, axis=0)


def process(rep, maxtb, stats, found):
    M = tuple(rep['M'])
    al, be = A4[rep['alpha']], A4[rep['beta']]
    Me1, Me2 = cols(M)
    zero4 = np.zeros(4, dtype=np.uint8)
    tA = [np.concatenate([al, zero4, [1, 0], Me1]).astype(np.uint8), np.concatenate([be, zero4, [0, 1], Me2]).astype(np.uint8)]
    Mp = minv(mdagger(M))
    Mp1, Mp2 = cols(Mp)
    TBs = enumA(Mp)
    stats['TB'] = stats.get('TB', 0) + len(TBs)
    print("rep", rep, "M'", Mp, "#T_B:", len(TBs), flush=True)
    t0 = time.time()
    for n, (ig, id_) in enumerate(TBs):
        if n >= maxtb:
            break
        ga, de = A4[ig], A4[id_]
        tB = [np.concatenate([zero4, ga, [1, 0], Mp1]).astype(np.uint8), np.concatenate([zero4, de, [0, 1], Mp2]).astype(np.uint8)]
        T = tA + tB
        # sanity: T isotropic
        for x in T:
            for y in T:
                assert H12(x, y) == 0, "T not isotropic"
        # words of T: weight >= 5 (T alone)
        rows, piv = rref(T)
        assert len(rows) == 4
        P = perp_of(T)
        assert len(P) == 65536, len(P)
        # reduce mod T
        V = P.copy()
        for r_, p_ in zip(rows, piv):
            coef = V[:, p_]
            V ^= MT[coef[:, None], r_[None, :]]
        key = np.zeros(len(V), dtype=np.int64)
        for k in range(12):
            key = key * 4 + V[:, k].astype(np.int64)
        uk, inv = np.unique(key, return_inverse=True)
        assert len(uk) == 256
        wts = WT12(P).astype(np.int16)
        mn = np.full(len(uk), 99, dtype=np.int16)
        np.minimum.at(mn, inv, wts)
        # coset reps
        first = np.zeros(len(uk), dtype=np.int64)
        first[inv[::-1]] = np.arange(len(V))[::-1]
        reps_v = V[first]  # canonical reps (reduced)
        iso = (H12(reps_v, reps_v) == 0)
        nz = (key[first] != 0)
        stats['iso_cosets'] = stats.get('iso_cosets', 0) + int((iso & nz).sum())
        good = np.nonzero(iso & nz & (mn >= 5))[0]
        stats['good_cosets'] = stats.get('good_cosets', 0) + len(good)
        if n == 0:
            print("first T_B: isotropic nonzero cosets", int((iso & nz).sum()), "(expect 135 = 45 points x 3)", "good (min wt>=5)", len(good), flush=True)
        if len(good) < 3:
            continue
        # lines: pairs of orthogonal good cosets, closure
        gv = reps_v[good]
        gk = {int(key[first][g]): i for i, g in enumerate(good)}
        lines = set()
        for i in range(len(good)):
            for j in range(i + 1, len(good)):
                if H12(gv[i], gv[j]) != 0:
                    continue
                ok = True
                members = []
                for (x, y) in [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3)]:
                    v = MT[x, gv[i]] ^ MT[y, gv[j]]
                    for r_, p_ in zip(rows, piv):
                        v = v ^ MT[v[p_], r_]
                    kk = 0
                    for k in range(12):
                        kk = kk * 4 + int(v[k])
                    if kk not in gk:
                        ok = False
                        break
                    members.append(gk[kk])
                if ok:
                    lines.add(frozenset(members))
        stats['lines'] = stats.get('lines', 0) + len(lines)
        if not lines:
            continue
        for ln in lines:
            ms = sorted(ln)
            u1, u2 = gv[ms[0]], gv[ms[1]]
            basis = T + [u1, u2]
            # all 4096 words
            coefs = np.array(list(itertools.product(range(4), repeat=6)), dtype=np.uint8)
            words = np.zeros((4096, 12), dtype=np.uint8)
            for k_, b_ in enumerate(basis):
                words ^= MT[coefs[:, k_][:, None], np.asarray(b_)[None, :]]
            wa = aidx(words[:, 0:4]); wb = aidx(words[:, 4:8])
            nz6 = words[:, 8:10].any(axis=1).astype(np.int16)
            nz7 = words[:, 10:12].any(axis=1).astype(np.int16)
            base = nz6 + nz7
            assert (WT12(words)[1:].min() >= 5), "L has a light word"
            stats['L_ok'] = stats.get('L_ok', 0) + 1
            # all (g, c != 0)
            best = 99
            for c in range(1, 16):
                WA = WTC[wa, c]  # 4096
                dc = GC[:, c]  # 720
                WBm = WTC[wb[None, :], dc[:, None]]  # 720 x 4096
                tot = WA[None, :] + WBm + base[None, :]
                mm = tot.min(axis=1)  # per g
                best = min(best, int(mm.min()))
                if (mm >= 5).any():
                    pass
                if c == 1:
                    okg = mm >= 5
                else:
                    okg &= (mm >= 5)
            stats['min_over_c_g_best'] = min(stats.get('min_over_c_g_best', 99), best)
            if okg.any():
                found.append((rep, n, [list(map(int, b_)) for b_ in basis], np.nonzero(okg)[0].tolist()))
                print("FOUND", found[-1], flush=True)
        if n % 200 == 0:
            print("  T_B", n, stats, round(time.time() - t0), "s", flush=True)


COEF6 = np.array(list(itertools.product(range(4), repeat=6)), dtype=np.uint8)
PROJ = [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3)]
ALL15 = [(x, y) for x in range(4) for y in range(4) if (x, y) != (0, 0)]


def process2(rep, maxtb, stats, found, tstart=0):
    """faster version of process: vectorized line finding; weight test via the (c, d=gc) table ok[c][d]."""
    M = tuple(rep['M'])
    al, be = A4[rep['alpha']], A4[rep['beta']]
    Me1, Me2 = cols(M)
    zero4 = np.zeros(4, dtype=np.uint8)
    tA = [np.concatenate([al, zero4, [1, 0], Me1]).astype(np.uint8), np.concatenate([be, zero4, [0, 1], Me2]).astype(np.uint8)]
    Mp = minv(mdagger(M))
    Mp1, Mp2 = cols(Mp)
    TBs = enumA(Mp)
    stats['TB'] = len(TBs)
    print("rep", rep, "M'", Mp, "#T_B:", len(TBs), flush=True)
    t0 = time.time()
    nlines_all = 0
    for n, (ig, id_) in enumerate(TBs):
        if n < tstart:
            continue
        if n >= maxtb:
            break
        ga, de = A4[ig], A4[id_]
        tB = [np.concatenate([zero4, ga, [1, 0], Mp1]).astype(np.uint8), np.concatenate([zero4, de, [0, 1], Mp2]).astype(np.uint8)]
        T = tA + tB
        if n == 0 or n == tstart:
            for x in T:
                for y in T:
                    assert H12(x, y) == 0, "T not isotropic"
        rows, piv = rref(T)
        assert len(rows) == 4
        P = perp_of(T)
        V = P.copy()
        for r_, p_ in zip(rows, piv):
            V ^= MT[V[:, p_][:, None], r_[None, :]]
        key = np.zeros(len(V), dtype=np.int64)
        for k in range(12):
            key = key * 4 + V[:, k].astype(np.int64)
        uk, inv = np.unique(key, return_inverse=True)
        assert len(uk) == 256
        wts = WT12(P).astype(np.int16)
        mn = np.full(len(uk), 99, dtype=np.int16)
        np.minimum.at(mn, inv, wts)
        first = np.zeros(len(uk), dtype=np.int64)
        first[inv[::-1]] = np.arange(len(V))[::-1]
        R = V[first]  # 256 reduced reps, sorted by key; R is an F_4-linear subspace (pivot cols zero)
        iso = (H12(R, R) == 0)
        nzm = (uk != 0)
        goodm = iso & nzm & (mn >= 5)
        gi = np.nonzero(goodm)[0]
        stats['good_cosets'] = stats.get('good_cosets', 0) + len(gi)
        if n == tstart:
            print("first T_B: isotropic nonzero cosets", int((iso & nzm).sum()), "(expect 135)", "good", len(gi), flush=True)
        if len(gi) < 3:
            continue
        O = (H12(R[gi][:, None, :], R[gi][None, :, :]) == 0)
        I, J = np.nonzero(np.triu(O, 1))
        if len(I) == 0:
            continue
        # 15 combos per pair
        pos_all = []
        okpair = np.ones(len(I), dtype=bool)
        for (x, y) in ALL15:
            cv = MT[x, R[gi][I]] ^ MT[y, R[gi][J]]
            kk = np.zeros(len(cv), dtype=np.int64)
            for k in range(12):
                kk = kk * 4 + cv[:, k].astype(np.int64)
            ps = np.searchsorted(uk, kk)
            ps = np.minimum(ps, len(uk) - 1)
            okpair &= (uk[ps] == kk) & goodm[ps]
            pos_all.append(ps)
        if not okpair.any():
            continue
        PA = np.sort(np.stack(pos_all, axis=1)[okpair], axis=1)
        lines = np.unique(PA, axis=0)
        lines = np.array([ln for ln in lines if len(set(ln.tolist())) == 15], dtype=np.int64).reshape(-1, 15)
        stats['lines'] = stats.get('lines', 0) + len(lines)
        if n == tstart:
            print("first T_B: weight-good lines", len(lines), "(<= 27)", flush=True)
        for ln in lines:
            # basis: two elements of the line (non-proportional): take first two distinct-projective members
            u1 = R[ln[0]]
            u2 = None
            for q in ln[1:]:
                cand = R[q]
                if H12(u1, cand) == 0 and not any((MT[s, u1] == cand).all() for s in range(1, 4)):
                    u2 = cand
                    break
            basis = T + [u1, u2]
            words = np.zeros((4096, 12), dtype=np.uint8)
            for k_, b_ in enumerate(basis):
                words ^= MT[COEF6[:, k_][:, None], np.asarray(b_)[None, :]]
            assert WT12(words)[1:].min() >= 5, "light word in L"
            stats['L_ok'] = stats.get('L_ok', 0) + 1
            wa = aidx(words[:, 0:4]); wb = aidx(words[:, 4:8])
            base = words[:, 8:10].any(axis=1).astype(np.int16) + words[:, 10:12].any(axis=1).astype(np.int16)
            okcd = np.zeros((16, 16), dtype=bool)
            for c in range(1, 16):
                WA = WTC[wa, c] + base
                for d in range(1, 16):
                    okcd[c, d] = (WA + WTC[wb, d]).min() >= 5
            okg = np.ones(720, dtype=bool)
            for c in range(1, 16):
                okg &= okcd[c, GC[:, c]]
            stats['cd_ok_max'] = max(stats.get('cd_ok_max', 0), int(okcd.sum()))
            stats['g_ok_max'] = max(stats.get('g_ok_max', 0), int(okg.sum()))
            if okg.any():
                found.append((rep, n, [list(map(int, b_)) for b_ in basis], np.nonzero(okg)[0].tolist()))
                print("FOUND", found[-1], flush=True)
        if n % 100 == 0:
            print("  T_B", n, stats, round(time.time() - t0), "s", flush=True)


def main():
    reps = json.load(open(os.path.join(os.environ.get("APPB_WORK", "."), "z3a_reps.json")))  # [release] written by z3a_stage1.py orbits in the working dir
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    maxtb = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    found = []
    for ri, rep in enumerate(reps):
        if which != "all" and int(which) != ri:
            continue
        stats = {}
        process2(rep, maxtb, stats, found, int(sys.argv[3]) if len(sys.argv) > 3 else 0)
        print("REP", ri, "DONE", stats, "found", len(found), flush=True)


if __name__ == "__main__":
    main()
