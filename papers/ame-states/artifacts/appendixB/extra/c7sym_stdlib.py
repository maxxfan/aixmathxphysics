"""[release] stdlib port of c7sym.py (written for the Appendix-B release; NOT part of the original run).

Same mathematics, same enumeration and the same printed summary as auditor/c7sym.py, without numpy, so that the l = 7 row
(single implementation in the original run) can also be replayed on a plain Python >= 3.9:
  sigma = 7-cycle on parties 0..6, party 7 fixed; S_nf = L1 + L1^perp (L1 an F_8-subspace of F_8^4), words
  x_i = (Tr(a_m z^i))_m + (Tr(b_m z^-i))_m, i = 0..6 (a in L1, b in L1^perp); weight(s + (c^7; f(c))) = #{i : s_i != c} + 1.
  For every L1 (r = dim 0..4): (i) all nonzero s have >= 5 nonzero symbols; (ii) no s has >= 4 symbols equal to one c != 0;
  min party weight of S = min(min_{s != 0} #nonzero symbols, 8 - max_{s, c != 0} #{i : s_i = c}).
Differences from c7sym.py (implementation only): L1^perp is computed by Gaussian elimination over F_8 instead of filtering
all 4096 vectors; words are packed 28-bit integers; the per-subspace scan stops as soon as its result can no longer change
the printed numbers (a word of weight <= 4 decides (i), and the minimum weight cannot exceed the best found so far).  The
r = 2 cross-check samples the SAME 40 subspaces as c7sym.py (identical random.Random(7) call sequence) and rebuilds each as a
16-qubit stabilizer checked with the verifier's rank test (lib/ame.bad_cuts -> build/rankcheck).
usage: python3 auditor/c7sym_stdlib.py [--r R] [--crosscheck K]
"""
import os, sys, itertools, json, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import ame

F = ame.GF(8)
zeta = F.exp[1]
MUL = [[F.mul(a, b) for b in range(8)] for a in range(8)]
INV = [0] + [next(b for b in range(1, 8) if MUL[a][b] == 1) for a in range(1, 8)]


def tr(x):
    t = x
    y = x
    for _ in range(2):
        y = F.mul(y, y)
        t = F.add(t, y)
    assert t in (0, 1)
    return t


TR1 = [[tr(F.mul(a, F.pow(zeta, i))) for i in range(7)] for a in range(8)]          # Tr(a z^i)
TR2 = [[tr(F.mul(a, F.pow(zeta, (-i) % 7))) for i in range(7)] for a in range(8)]   # Tr(a z^-i)
ALL4 = list(itertools.product(range(8), repeat=4))                                 # same order as c7sym.ALL4


def vidx(v):
    return ((v[0] * 8 + v[1]) * 8 + v[2]) * 8 + v[3]


def packed(v, TR):
    """7 party symbols (bit m of symbol i = Tr(v_m z^{+-i})) packed as nibbles of a 28-bit int."""
    w = 0
    for i in range(7):
        s = 0
        for m in range(4):
            s |= TR[v[m]][i] << m
        w |= s << (4 * i)
    return w


SA = [packed(v, TR1) for v in ALL4]
SB = [packed(v, TR2) for v in ALL4]
NZ16 = [sum(1 for k in range(4) if (x >> (4 * k)) & 15) for x in range(1 << 16)]


def nnz(w):
    return NZ16[w & 0xFFFF] + NZ16[w >> 16]


def maxeq(w):
    cnt = [0] * 16
    for i in range(7):
        cnt[(w >> (4 * i)) & 15] += 1
    cnt[0] = 0
    return max(cnt)


def rref_subspaces(r):   # verbatim from c7sym.py
    for piv in itertools.combinations(range(4), r):
        free_pos = []
        for k, p in enumerate(piv):
            free_pos.append([c for c in range(p + 1, 4) if c not in piv])
        for choice in itertools.product(*[list(itertools.product(range(8), repeat=len(fp))) for fp in free_pos]):
            rows = []
            for k, p in enumerate(piv):
                row = [0] * 4
                row[p] = 1
                for c, v in zip(free_pos[k], choice[k]):
                    row[c] = v
                rows.append(row)
            yield rows


def span(rows):
    out = [(0, 0, 0, 0)]
    for v in rows:
        out = [tuple(o[k] ^ MUL[c][v[k]] for k in range(4)) for c in range(8) for o in out]
    return out


def ip(a, b):
    return MUL[a[0]][b[1]] ^ MUL[a[1]][b[0]] ^ MUL[a[2]][b[3]] ^ MUL[a[3]][b[2]]


def perp_basis(L1):
    """basis of {b : <v,b> = 0 for v in L1}: nullspace of the rows (v1, v0, v3, v2) over F_8."""
    A = [[v[1], v[0], v[3], v[2]] for v in L1]
    piv_cols = []
    r = 0
    for col in range(4):
        p = next((i for i in range(r, len(A)) if A[i][col]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        iv = INV[A[r][col]]
        A[r] = [MUL[iv][x] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][col]:
                f = A[i][col]
                A[i] = [x ^ MUL[f][y] for x, y in zip(A[i], A[r])]
        piv_cols.append(col)
        r += 1
    free = [c for c in range(4) if c not in piv_cols]
    basis = []
    for fc in free:
        b = [0] * 4
        b[fc] = 1
        for i, pc in enumerate(piv_cols):
            b[pc] = A[i][fc]          # char 2: -x = x
        basis.append(b)
    return basis


# ---- verbatim (pure Python) helpers of c7sym.py for the cross-check
def gen_words(L1rows, L2rows):
    out = []
    for v in L1rows:
        for e in (1, zeta, F.mul(zeta, zeta)):
            a = [F.mul(e, c) for c in v]
            w = [sum(TR1[a[m]][i] << m for m in range(4)) for i in range(7)] + [0]
            out.append(w)
    for v in L2rows:
        for e in (1, zeta, F.mul(zeta, zeta)):
            b = [F.mul(e, c) for c in v]
            w = [sum(TR2[b[m]][i] << m for m in range(4)) for i in range(7)] + [0]
            out.append(w)
    for m in range(4):
        c = 1 << m
        out.append([c] * 7 + [c])
    return out


def to_rows(words):
    rows = []
    for w in words:
        x = [0] * 16
        z = [0] * 16
        for j, s in enumerate(w):
            u = [(s >> m) & 1 for m in range(4)]
            x[2 * j], z[2 * j], x[2 * j + 1], z[2 * j + 1] = u[0], u[1], u[2], u[3]
        rows.append(x + z)
    return rows


def rank2(rows):
    rows = [int("".join(map(str, r)), 2) for r in rows]
    rk = 0
    for bit in range(31, -1, -1):
        piv = None
        for i in range(rk, len(rows)):
            if (rows[i] >> bit) & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and (rows[i] >> bit) & 1:
                rows[i] ^= rows[rk]
        rk += 1
    return rk


def isotropic(rows):
    for a in rows:
        for b in rows:
            s = 0
            for k in range(16):
                s ^= (a[k] & b[16 + k]) ^ (a[16 + k] & b[k])
            if s:
                return False
    return True


def scan(a_idx, b_idx, best):
    """returns (ok_i, ok_ii or None, minw upper bound); minw is exact whenever it exceeds `best`."""
    ok_i = True
    m = 99           # running min of S's party weight
    nz_min = 99
    me_max = 0
    for ia in a_idx:
        wa = SA[ia]
        for ib in b_idx:
            w = wa ^ SB[ib]
            if w:
                k = nnz(w)
                if k < nz_min:
                    nz_min = k
                    if k < 5:
                        ok_i = False
            if m > best:                 # only then can this subspace still raise the best value: need exact minw
                e = maxeq(w)
                if e > me_max:
                    me_max = e
            m = min(nz_min, 8 - me_max)
            if not ok_i and m <= best:
                return False, None, m
    ok_ii = (me_max <= 3) if ok_i else None
    return ok_i, ok_ii, m


def main():
    rsel = None
    xcheck = 40
    if "--r" in sys.argv:
        rsel = int(sys.argv[sys.argv.index("--r") + 1])
    if "--crosscheck" in sys.argv:
        xcheck = int(sys.argv[sys.argv.index("--crosscheck") + 1])
    rng = random.Random(7)
    for r in ([rsel] if rsel is not None else [0, 1, 2, 3, 4]):
        n_lines = n_i = n_both = 0
        best_minw = 0
        checked = []
        for L1 in rref_subspaces(r) if r > 0 else [[]]:
            n_lines += 1
            a_idx = [vidx(v) for v in span(L1)] if r > 0 else [0]
            pb = perp_basis(L1)
            b_idx = sorted(vidx(v) for v in span(pb))
            assert len(b_idx) == 8 ** (4 - r) and len(set(b_idx)) == len(b_idx)
            ok_i, ok_ii, minw = scan(a_idx, b_idx, best_minw)
            if ok_i:     # never happens; keep the exact semantics of c7sym.criterion anyway
                ok_ii = all(maxeq(SA[ia] ^ SB[ib]) <= 3 for ia in a_idx for ib in b_idx)
            best_minw = max(best_minw, minw)
            n_i += ok_i
            n_both += bool(ok_i and ok_ii)
            if r == 2 and len(checked) < xcheck and (rng.random() < 0.02 or bool(ok_i and ok_ii) or (ok_i and not ok_ii and rng.random() < 0.05)):
                B = [ALL4[i] for i in b_idx]          # ALL4 order, as c7sym's b_el
                basis = []
                for b in B:
                    if len(basis) == 2:
                        break
                    if not any(b):
                        continue
                    cand = basis + [list(b)]
                    if len(set(span(cand))) == 8 ** len(cand):
                        basis = cand
                rows = to_rows(gen_words(L1, basis))
                ok_rank = rank2(rows) == 16
                ok_iso = isotropic(rows)
                bad, cuts = ame.bad_cuts(rows, 2, 2, 8) if (ok_rank and ok_iso) else (-1, -1)
                checked.append((bool(ok_i), bool(ok_ii), ok_rank, ok_iso, bad))
        summary = {"lines": n_lines, "pass_weight5": n_i, "pass_both": n_both,
                   "best_min_party_weight_of_S_over_all_lines": best_minw, "examples": [],
                   "crosscheck(ok_i,ok_ii,rank16,isotropic,bad_cuts)": checked}
        print(json.dumps({r: summary}), flush=True)


if __name__ == "__main__":
    main()
