"""Auditor: exhaustive classification of stabilizer AME(8,4) (party = 2 qubits = F_2^4, symplectic form u1v2+u2v1+u3v4+u4v3)
with a PURE block-permutation symmetry sigma = 7-cycle on parties 0..6 (party 7 fixed).  Derivation (results/argument.md 5.4):
 V = V^sigma + V_nf, S = S^sigma + S_nf (l = 7 odd, sigma symplectic => orthogonal sum, both Lagrangian).
 V_nf = V1 + V2, V1 = {x_i = (Tr(a_m z^i))_m, a in F_8^4}, V2 = {x_i = (Tr(b_m z^-i))_m}, B(x_a, y_b) = Tr(<a,b>_F8), V1, V2 isotropic.
 S_nf = L1 + L1^perp (L1 an F_8-subspace of F_8^4, perp for the F_8 symplectic form); S^sigma = {(c,c,...,c ; f(c))}, f in Sp(4,2)
 (elements of S^sigma with c = 0 would have weight 1).  weight(s + (c^7;f(c))) = #{i: s_i != c} + 1  for c != 0, so
 AME  <=>  (i) every nonzero s in S_nf has weight >= 5 on the 7-cycle and (ii) for every s in S_nf, c != 0: #{i : s_i = c} <= 3.
 The criterion does NOT involve f.  dim_F8 L1 = 2 is forced by (i) (a plane/point in PG(3,8) always meets PG(3,2) for r = 3 or by duality r = 1); the script
 enumerates ALL r in 0..4 anyway.  For passing lines the full stabilizer is written (f = id) and checked with ame.bad_cuts.
usage: python3 work/auditor/c7sym.py [--r R] [--crosscheck K]
"""
import os, sys, itertools, json, random
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
# [release] was: ROOT = <problem dir>; sys.path.insert(0, ROOT); from qlib import ame   (ame.py now lives in ../lib)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import ame

F = ame.GF(8)
zeta = F.exp[1]
MUL = np.array([[F.mul(a, b) for b in range(8)] for a in range(8)], dtype=np.int64)
def tr(x):
    t = x
    y = x
    for _ in range(2):
        y = F.mul(y, y)
        t = F.add(t, y)
    assert t in (0, 1)
    return t
TR1 = np.array([[tr(F.mul(a, F.pow(zeta, i))) for i in range(7)] for a in range(8)], dtype=np.int64)       # Tr(a z^i)
TR2 = np.array([[tr(F.mul(a, F.pow(zeta, (-i) % 7))) for i in range(7)] for a in range(8)], dtype=np.int64)  # Tr(a z^-i)

ALL4 = np.array(list(itertools.product(range(8), repeat=4)), dtype=np.int64)  # (4096, 4)


def ip(a, B):
    """<a,b> = a1 b2 + a2 b1 + a3 b4 + a4 b3 over F_8 for one a and an array B (k,4)."""
    return MUL[a[0], B[:, 1]] ^ MUL[a[1], B[:, 0]] ^ MUL[a[2], B[:, 3]] ^ MUL[a[3], B[:, 2]]


def span(rows):
    """all F_8-combinations of the given row vectors, (8^r, 4)."""
    r = len(rows)
    out = np.zeros((1, 4), dtype=np.int64)
    for v in rows:
        v = np.array(v, dtype=np.int64)
        new = []
        for c in range(8):
            new.append(out ^ MUL[c, v][None, :])
        out = np.concatenate(new, axis=0)
    return out


def rref_subspaces(r):
    """all r-dim subspaces of F_8^4 as lists of r rows in RREF."""
    for piv in itertools.combinations(range(4), r):
        free_pos = []   # per row the free columns: > pivot and not a pivot column
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


def words_of(a_el, b_el):
    """S_nf words as party symbols (64 x 64 x 7) ints 0..15 (bit m = u_{m+1}); a_el, b_el (k,4) element lists."""
    ka, kb = len(a_el), len(b_el)
    Abits = np.stack([TR1[a_el[:, m]] for m in range(4)], axis=-1)   # (ka,7,4)
    Bbits = np.stack([TR2[b_el[:, m]] for m in range(4)], axis=-1)   # (kb,7,4)
    W = Abits[:, None, :, :] ^ Bbits[None, :, :, :]                    # (ka,kb,7,4)
    sym = W[..., 0] | (W[..., 1] << 1) | (W[..., 2] << 2) | (W[..., 3] << 3)
    return sym.reshape(ka * kb, 7)


def criterion(sym):
    """sym: (K,7) symbols of ALL words of S_nf (including 0). returns (ok_i, ok_ii)."""
    nz = (sym != 0).sum(axis=1)
    nonzero_word = (sym != 0).any(axis=1)
    ok_i = bool((nz[nonzero_word] >= 5).all())
    ok_ii = True
    for c in range(1, 16):
        if ((sym == c).sum(axis=1) > 3).any():
            ok_ii = False
            break
    return ok_i, ok_ii


def to_rows(words, f_images):
    """16 generator words, each a list of 8 party symbols -> canonical rows (x|z), party j owns qubits 2j, 2j+1."""
    rows = []
    for w in words:
        x = [0] * 16
        z = [0] * 16
        for j, s in enumerate(w):
            u = [(s >> m) & 1 for m in range(4)]
            x[2 * j], z[2 * j], x[2 * j + 1], z[2 * j + 1] = u[0], u[1], u[2], u[3]
        rows.append(x + z)
    return rows


def gen_words(L1rows, L2rows):
    """12 basis words of S_nf (7 symbols + party 7 = 0) and 4 words of S^sigma (f = id)."""
    out = []
    for v in L1rows:
        for e in (1, zeta, F.mul(zeta, zeta)):
            a = [F.mul(e, c) for c in v]
            w = [sum(TR1[a[m], i] << m for m in range(4)) for i in range(7)] + [0]
            out.append(w)
    for v in L2rows:
        for e in (1, zeta, F.mul(zeta, zeta)):
            b = [F.mul(e, c) for c in v]
            w = [sum(TR2[b[m], i] << m for m in range(4)) for i in range(7)] + [0]
            out.append(w)
    for m in range(4):
        c = 1 << m
        out.append([c] * 7 + [c])
    return out


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


def main():
    rsel = None
    xcheck = 40
    if "--r" in sys.argv:
        rsel = int(sys.argv[sys.argv.index("--r") + 1])
    if "--crosscheck" in sys.argv:
        xcheck = int(sys.argv[sys.argv.index("--crosscheck") + 1])
    rng = random.Random(7)
    summary = {}
    for r in ([rsel] if rsel is not None else [0, 1, 2, 3, 4]):
        n_lines = n_i = n_both = 0
        best_minw = 0
        examples = []
        checked = []
        for L1 in rref_subspaces(r) if r > 0 else [[]]:
            n_lines += 1
            a_el = span(L1) if r > 0 else np.zeros((1, 4), dtype=np.int64)
            mask = np.ones(len(ALL4), dtype=bool)
            for v in L1:
                mask &= (ip(v, ALL4) == 0)
            b_el = ALL4[mask]
            assert len(b_el) == 8 ** (4 - r), (len(b_el), r)
            sym = words_of(a_el, b_el)
            ok_i, ok_ii = criterion(sym)
            nzw = (sym != 0).sum(axis=1)
            minw = int(nzw[(sym != 0).any(axis=1)].min()) if (sym != 0).any() else 99
            for c in range(1, 16):
                minw = min(minw, int(((sym != c).sum(axis=1) + 1).min()))
            best_minw = max(best_minw, minw)
            n_i += ok_i
            n_both += (ok_i and ok_ii)
            if ok_i and ok_ii and len(examples) < 3:
                examples.append(L1)
            # cross-check a sample against the exact C rank check (only r = 2 words are meaningful for a dim-16 S)
            if r == 2 and len(checked) < xcheck and (rng.random() < 0.02 or (ok_i and ok_ii) or (ok_i and not ok_ii and rng.random() < 0.05)):
                L2 = []
                # basis of L1-perp: pick 2 independent vectors from b_el
                B = b_el
                basis = []
                for b in B:
                    if len(basis) == 2:
                        break
                    if not b.any():
                        continue
                    cand = basis + [list(map(int, b))]
                    sp = np.unique(span(cand), axis=0)
                    if len(sp) == 8 ** len(cand):
                        basis = cand
                L2 = basis
                words = gen_words(L1, L2)
                rows = to_rows(words, None)
                ok_rank = rank2(rows) == 16
                ok_iso = isotropic(rows)
                bad, cuts = ame.bad_cuts(rows, 2, 2, 8) if (ok_rank and ok_iso) else (-1, -1)
                checked.append((ok_i, ok_ii, ok_rank, ok_iso, bad))
        summary[r] = {"lines": n_lines, "pass_weight5": n_i, "pass_both": n_both, "best_min_party_weight_of_S_over_all_lines": best_minw, "examples": [list(map(list, e)) for e in examples], "crosscheck(ok_i,ok_ii,rank16,isotropic,bad_cuts)": checked}
        print(json.dumps({r: summary[r]}), flush=True)


if __name__ == "__main__":
    main()
