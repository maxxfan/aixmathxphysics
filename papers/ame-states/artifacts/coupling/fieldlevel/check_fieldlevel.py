#!/usr/bin/env python3
"""Field-level check of the 27 coupling-rule ingredients, independent of any conversion to F_p stabilizer rows.

For every fieldlevel/<id>.json this script builds GF(Q) from the modulus printed in the file (own arithmetic, no
primitive-polynomial search, no code shared with convert.py), rebuilds the matrix from the structural parameters
(circulant first rows, bordered-circulant parameters, evaluation points, kernel formula) and compares it with the
printed matrix, and verifies the claimed property by exact rank conditions:

  hermitian_systematic  C = rowspace [I | A] over F_Q, Q = q^2, conj(z) = z^q.
      Hermitian self-duality: A conj(A)^T = -I.  The F_q-linear stabilizer state of C is t-uniform iff
      d(C) >= t+1 iff every n-t columns of [I | A] have rank n/2 (a word supported on a t-set T exists iff the
      columns outside T are dependent on the row side).  AME(n, q): t = n/2, every n/2 columns independent (MDS).
  css_minimal_support   |C> = sum_{c in C} |c>, C = rowspace G over F_q, stabilizers X(c), c in C, Z(h), h in C^perp.
      t-uniform iff d(C) >= t+1 (every n-t columns of G have rank dim C) and d(C^perp) >= t+1 (every t columns of
      G are linearly independent).  AME(n, q), dim C = n/2: both conditions say every n/2 columns have rank n/2.
  fq_graph              Gamma symmetric with zero diagonal over F_q; t-uniform iff Gamma[A, A^c] has full row rank t
      for every set A of t parties (Lemma 2 of the paper, stated there for t = floor(n/2)).

For the k-uniform ingredients that are not AME, the exact distance is also shown: some (n-t-1)-column set is
rank-deficient, so d = t+1 exactly.  Finally every cell of ../cells.json is listed with the verdicts of its two
ingredients and the arithmetic conditions of the coupling lemma (n = 2k+1, ingredients on 2k+2 parties, G_a
k-uniform, G_b AME, a <= b, d = ab).

Self-contained: Python 3.9 standard library only.  usage: python3 check_fieldlevel.py [id ...]
"""
import itertools
import json
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


# ----------------------------------------------------------------------------------------------------------------
# GF(p^m) from an explicit modulus (independent implementation)
# ----------------------------------------------------------------------------------------------------------------
class Field:
    """GF(p^m) = F_p[x]/(f), f = modulus (coefficients c_0..c_m, monic); the int e = sum c_i p^i is sum c_i x^i.
    For m = 1 (modulus None) this is F_p on the ints 0..p-1."""

    def __init__(self, p, m, modulus):
        self.p, self.m, self.Q = p, m, p ** m
        Q = self.Q
        if m == 1:
            assert modulus is None
            self.mul_t = [[(a * b) % p for b in range(Q)] for a in range(Q)]
            self.add_t = [[(a + b) % p for b in range(Q)] for a in range(Q)]
        else:
            f = list(modulus)
            assert len(f) == m + 1 and f[-1] == 1, "modulus must be monic of degree m"
            self.f = f
            dig = [self.digits(a) for a in range(Q)]
            self.add_t = [[self.undigits([(x + y) % p for x, y in zip(dig[a], dig[b])]) for b in range(Q)] for a in range(Q)]
            self.mul_t = [[self.undigits(self.polymulmod(dig[a], dig[b])) for b in range(Q)] for a in range(Q)]
        self.neg_t = [next(b for b in range(Q) if self.add_t[a][b] == 0) for a in range(Q)]
        self.sub_t = [[self.add_t[a][self.neg_t[b]] for b in range(Q)] for a in range(Q)]
        # field test: every nonzero element has a multiplicative inverse (=> modulus irreducible, F is a field)
        self.inv_t = [0] * Q
        for a in range(1, Q):
            inv = [b for b in range(1, Q) if self.mul_t[a][b] == 1]
            if len(inv) != 1:
                raise ValueError("modulus is not irreducible: %d has no inverse" % a)
            self.inv_t[a] = inv[0]

    def digits(self, a):
        return [(a // self.p ** i) % self.p for i in range(self.m)]

    def undigits(self, ds):
        return sum(d * self.p ** i for i, d in enumerate(ds))

    def polymulmod(self, u, v):
        p, m, f = self.p, self.m, self.f
        prod = [0] * (2 * m - 1)
        for i, x in enumerate(u):
            if x:
                for j, y in enumerate(v):
                    prod[i + j] = (prod[i + j] + x * y) % p
        for d in range(2 * m - 2, m - 1, -1):   # reduce x^d = x^(d-m) * (-(f_0 + .. + f_{m-1} x^{m-1}))
            c = prod[d]
            if c:
                prod[d] = 0
                for i in range(m):
                    prod[d - m + i] = (prod[d - m + i] - c * f[i]) % p
        return prod[:m]

    def pow(self, a, e):
        r = 1
        for _ in range(e):
            r = self.mul_t[r][a]
        return r

    def fast_pow(self, a, e):
        r, b = 1, a
        while e:
            if e & 1:
                r = self.mul_t[r][b]
            b = self.mul_t[b][b]
            e >>= 1
        return r

    def order(self, a):
        o, y = 1, a
        while y != 1:
            y = self.mul_t[y][a]
            o += 1
        return o


def rank(M, F):
    """rank over F of the matrix M (list of rows), forward elimination"""
    M = [r[:] for r in M]
    if not M:
        return 0
    rows, cols = len(M), len(M[0])
    mul, sub, inv = F.mul_t, F.sub_t, F.inv_t
    r = 0
    for c in range(cols):
        pr = None
        for i in range(r, rows):
            if M[i][c]:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = inv[M[r][c]]
        piv = [mul[iv][x] for x in M[r]]
        M[r] = piv
        for i in range(r + 1, rows):
            f = M[i][c]
            if f:
                mf = mul[f]
                Mi = M[i]
                M[i] = [sub[x][mf[y]] for x, y in zip(Mi, piv)]
        r += 1
        if r == rows:
            break
    return r


def rref_key(M, F):
    """reduced row echelon form (for row-space equality)"""
    M = [r[:] for r in M]
    rows, cols = len(M), len(M[0])
    r = 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if M[i][c]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = F.inv_t[M[r][c]]
        M[r] = [F.mul_t[iv][x] for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [F.sub_t[x][F.mul_t[f][y]] for x, y in zip(M[i], M[r])]
        r += 1
    return M[:r]


def columns(M, S):
    return [[row[j] for j in S] for row in M]


def count_full_rank_column_sets(G, F, size, need):
    """(sets checked, sets whose column block G[:, S] has rank `need`) over all |S| = size"""
    n = len(G[0])
    tot = ok = 0
    for S in itertools.combinations(range(n), size):
        tot += 1
        if rank(columns(G, S), F) == need:
            ok += 1
    return tot, ok


# ----------------------------------------------------------------------------------------------------------------
# structure rebuilds
# ----------------------------------------------------------------------------------------------------------------
def circulant(first_row):
    c = len(first_row)
    return [[first_row[(j - i) % c] for j in range(c)] for i in range(c)]


def bordered_circulant(st):
    c = st["c"]
    core = circulant(st["core_first_row"])
    assert len(core) == c
    A = [[st["corner_alpha"]] + [st["row_border_beta"]] * c]
    for i in range(c):
        A.append([st["column_border_gamma"]] + core[i])
    return A


def block_circulant(st):
    b = st["block_size"]
    blocks = st["blocks_first_rows"]
    nb = len(blocks)
    A = [[0] * (nb * b) for _ in range(nb * b)]
    for r in range(nb):
        for s in range(nb):
            for x in range(b):
                for y in range(b):
                    A[b * r + x][b * s + y] = blocks[r][s][(y - x) % b]
    return A


def hz_kernel(F, q, g):
    """Hontarenko-Zyczkowski kernel A_q (as printed in the JSON's structure.kernel)"""
    gam = F.fast_pow(g, q + 1)
    two = 2 % F.p
    j = next(j for j in range(q - 1) if F.fast_pow(gam, j) == two)
    b = F.fast_pow(g, j)
    bb = F.fast_pow(b, q)
    xs = [F.fast_pow(gam, i) for i in range(q - 1)]
    K = [[0] * q for _ in range(q)]
    K[0][0] = F.neg_t[F.mul_t[b][F.inv_t[bb]]]
    for c in range(q - 1):
        K[0][c + 1] = b
        K[c + 1][0] = b
    for a, x in enumerate(xs):
        for c, y in enumerate(xs):
            num = F.mul_t[F.mul_t[two][y]][F.add_t[y][F.mul_t[g][x]]]
            den = F.sub_t[F.mul_t[gam][F.mul_t[x][x]]][F.mul_t[y][y]]
            K[a + 1][c + 1] = F.mul_t[num][F.inv_t[den]]
    return K, {"g": g, "gam": gam, "b": b, "j": j}


def poly_str(cs):
    """[c_0, .., c_d] -> 'x^d + .. + c_0'"""
    terms = []
    for i in range(len(cs) - 1, -1, -1):
        cf = cs[i]
        if not cf:
            continue
        mon = "" if i == 0 else ("x" if i == 1 else "x^%d" % i)
        terms.append(("%d" % cf) if i == 0 else (mon if cf == 1 else "%d*%s" % (cf, mon)))
    return " + ".join(terms)


def herm_gram(A, F, q):
    """A conj(A)^T"""
    k = len(A)
    conj = [F.fast_pow(x, q) for x in range(F.Q)]
    out = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            s = 0
            for t in range(len(A[i])):
                s = F.add_t[s][F.mul_t[A[i][t]][conj[A[j][t]]]]
            out[i][j] = s
    return out


def check_structure(fl, F, lines):
    """rebuild the printed matrix from the structural parameters; returns True if equal"""
    st = fl.get("structure", {})
    t = st.get("type")
    kind = fl["kind"]
    key = {"hermitian_systematic": "A", "css_minimal_support": "G", "fq_graph": "Gamma"}[kind]
    M = fl[key]
    ok = True
    if t == "circulant":
        R = circulant(st["first_row"])
        desc = "circulant, first row %s" % st["first_row"]
    elif t == "bordered_circulant":
        R = bordered_circulant(st)
        desc = "bordered circulant, alpha=%d beta=%d gamma=%d core first row %s" % (
            st["corner_alpha"], st["row_border_beta"], st["column_border_gamma"], st["core_first_row"])
        if "core_exponents" in st:     # Hontarenko-Zyczkowski parent: core = g^exponents
            g = st["g"]
            okg = F.order(g) == F.Q - 1 and [F.fast_pow(g, e) for e in st["core_exponents"]] == st["core_first_row"]
            ok = ok and okg
            desc += "; core = g^%s with g = %d of order %d: %s" % (st["core_exponents"], g, F.order(g),
                                                                  "yes" if okg else "NO")
        if "corner_relation" in st:    # alpha = -conj(F0) beta / conj(beta), F0 = sum of the core first row
            q = fl["q"]
            F0 = 0
            for x in st["core_first_row"]:
                F0 = F.add_t[F0][x]
            b = st["row_border_beta"]
            al = F.neg_t[F.mul_t[F.mul_t[F.fast_pow(F0, q)][b]][F.inv_t[F.fast_pow(b, q)]]]
            okc = al == st["corner_alpha"]
            ok = ok and okc
            desc += "; alpha = -conj(sum core) beta/conj(beta): %s" % ("yes" if okc else "NO")
    elif t == "block_circulant":
        R = block_circulant(st)
        desc = "2x2 blocks of %dx%d circulants, first rows %s" % (st["block_size"], st["block_size"], st["blocks_first_rows"])
    elif t == "hz_kernel_submatrix":
        q = fl["q"]
        K, kp = hz_kernel(F, q, st["kernel_params"]["g"])
        okK = K == st["kernel_matrix"] and kp == st["kernel_params"]
        gram = herm_gram(K, F, q)
        mone = F.neg_t[1]
        okU = all(gram[i][j] == (mone if i == j else 0) for i in range(q) for j in range(q))
        P = [[K[i][j] for j in st["cols"]] for i in st["rows"]]
        lam, D = st["Lambda"], st["D"]
        R = [[F.mul_t[F.mul_t[lam[i]][P[i][j]]][D[j]] for j in range(len(D))] for i in range(len(lam))]
        ok = ok and okK and okU
        desc = ("Hontarenko-Zyczkowski kernel A_%d rebuilt from its formula (g=%d): %s; kernel K conj(K)^T = -I: %s; "
                "A = diag(Lambda) K[rows %s][cols %s] diag(D)" % (q, kp["g"], "equal to the printed kernel" if okK else
                                                                  "DIFFERENT", "yes" if okU else "NO", st["rows"], st["cols"]))
    elif t == "explicit_labels":
        # printed labels e -> (e mod p) + (e div p) r, r a root of the printed minimal polynomial
        r, mp = st["r"], st["r_minimal_polynomial"]
        val = 0
        for i, cf in enumerate(mp):
            val = F.add_t[val][F.mul_t[cf % F.p][F.fast_pow(r, i)]]
        okr = val == 0
        ok = ok and okr
        R = [[F.add_t[e % F.p][F.mul_t[(e // F.p) % F.p][r]] for e in row] for row in st["source_labels"]]
        desc = ("explicit printed matrix (%s); labels e -> (e mod %d) + (e div %d) r with r = %d a root of %s: %s"
                % (st.get("description", ""), F.p, F.p, r, poly_str(mp), "yes" if okr else "NO"))
    elif t in ("reed_solomon", "doubly_extended_reed_solomon"):
        k = st["k"]
        pts = st["evaluation_points"]
        R = []
        for i in range(k):
            row = [F.fast_pow(x, i) if x else (1 if i == 0 else 0) for x in pts]
            if st["point_at_infinity"]:
                row.append(1 if i == k - 1 else 0)
            R.append(row)
        desc = "%s [%d,%d] over F_%d, evaluation points %s%s" % (
            "doubly extended RS" if st["point_at_infinity"] else "RS", len(R[0]), k, F.Q, pts,
            " + infinity" if st["point_at_infinity"] else "")
        if "G_systematic" in fl:
            oks = rref_key(fl["G_systematic"], F) == rref_key(M, F)
            ok = ok and oks
            desc += "; G_systematic spans the same code: %s" % ("yes" if oks else "NO")
    elif t == "systematic_bordered_circulant":
        P = bordered_circulant(st)
        k = len(P)
        R = [[1 if i == j else 0 for j in range(k)] + P[i] for i in range(k)]
        desc = "[I_%d | P], P bordered circulant alpha=%d beta=%d gamma=%d core first row %s" % (
            k, st["corner_alpha"], st["row_border_beta"], st["column_border_gamma"], st["core_first_row"])
    elif t == "wheel":
        n = fl["n"]
        R = [[0] * n for _ in range(n)]
        for i in range(n - 1):
            j = (i + 1) % (n - 1)
            R[i][j] = R[j][i] = 1
            R[i][n - 1] = R[n - 1][i] = 1
        desc = "wheel graph: %d-cycle plus hub %d" % (n - 1, n - 1)
    else:
        raise ValueError("unknown structure type %r" % t)
    same = R == M
    ok = ok and same
    lines.append("  structure: %s" % desc)
    lines.append("  printed %s rebuilt from these parameters: %s" % (key, "equal" if same else "DIFFERENT"))
    return ok


# ----------------------------------------------------------------------------------------------------------------
# per-kind property checks
# ----------------------------------------------------------------------------------------------------------------
def check_one(fl):
    lines = []
    fd = fl["field"]
    F = Field(fd["p"], fd["m"], fd["modulus"])
    assert F.Q == fd["Q"]
    n, kind = fl["n"], fl["kind"]
    t = fl["claimed_property"]["uniformity_k"]
    if F.m == 1:
        finfo = "field F_%d (prime field): all %d nonzero elements invertible" % (F.Q, F.Q - 1)
    else:
        finfo = "field F_%d = F_%d[x]/(%s): all %d nonzero elements invertible (modulus irreducible)" % (
            F.Q, F.p, fd["modulus_str"], F.Q - 1)
    ok = True
    if fd.get("x_is_primitive"):
        ox = F.order(F.p)          # the int p is the element x
        finfo += ", x primitive (order %d)" % ox
        ok = ok and ox == F.Q - 1
    lines.append("  " + finfo)
    ok = check_structure(fl, F, lines) and ok

    if kind == "hermitian_systematic":
        q = fl["q"]
        assert q * q == F.Q
        A = fl["A"]
        k = len(A)
        assert 2 * k == n and all(len(r) == k for r in A)
        gram = herm_gram(A, F, q)
        mone = F.neg_t[1]
        herm = all(gram[i][j] == (mone if i == j else 0) for i in range(k) for j in range(k))
        lines.append("  Hermitian self-duality A conj(A)^T = -I (conj z = z^%d, %dx%d): %s" % (q, k, k, "yes" if herm else "NO"))
        ok = ok and herm
        G = [[1 if i == j else 0 for j in range(k)] + list(A[i]) for i in range(k)]
        dim = k
        tot, good = count_full_rank_column_sets(G, F, n - t, dim)
        lines.append("  d(C) >= %d (state %d-uniform): column sets of size n-%d = %d of [I|A]: %d checked, %d of rank %d"
                     % (t + 1, t, t, n - t, tot, good, dim))
        ok = ok and tot == good == math.comb(n, n - t)
        verdict_obj = "[%d,%d,%d]_%d Hermitian self-dual" % (n, k, fl["code"]["distance"], F.Q)
    elif kind == "css_minimal_support":
        G = fl["G"]
        dim = rank(G, F)
        lines.append("  C = rowspace G over F_%d: dim C = %d" % (F.Q, dim))
        if "duality" in fl.get("code", {}) and "self-dual" in fl["code"]["duality"]:
            sd = all(_dot(u, v, F) == 0 for u in G for v in G) and 2 * dim == n
            lines.append("  Euclidean self-duality G G^T = 0 and dim C = n/2: %s" % ("yes" if sd else "NO"))
            ok = ok and sd
        tot1, good1 = count_full_rank_column_sets(G, F, n - t, dim)
        lines.append("  d(C) >= %d: column sets of size n-%d = %d of G: %d checked, %d of rank %d"
                     % (t + 1, t, n - t, tot1, good1, dim))
        ok = ok and tot1 == good1
        if n - t != t:
            tot2, good2 = count_full_rank_column_sets(G, F, t, t)
            lines.append("  d(C^perp) >= %d: column sets of size %d of G: %d checked, %d of rank %d (independent)"
                         % (t + 1, t, tot2, good2, t))
            ok = ok and tot2 == good2
        else:
            lines.append("  d(C^perp) >= %d: same condition (every %d columns of G independent, dim C = n/2)" % (t + 1, t))
        verdict_obj = "C = [%d,%d,%d]_%d" % (n, dim, fl["code"]["distance"], F.Q)
    elif kind == "fq_graph":
        Gam = fl["Gamma"]
        sym = all(Gam[i][j] == Gam[j][i] for i in range(n) for j in range(n)) and all(Gam[i][i] == 0 for i in range(n))
        lines.append("  Gamma symmetric with zero diagonal: %s" % ("yes" if sym else "NO"))
        ok = ok and sym
        tot = good = 0
        for S in itertools.combinations(range(n), t):
            Sc = [j for j in range(n) if j not in S]
            tot += 1
            good += rank([[Gam[i][j] for j in Sc] for i in S], F) == t
        lines.append("  Lemma 2: Gamma[A, A^c] of full row rank %d for |A| = %d: %d sets checked, %d of full row rank"
                     % (t, t, tot, good))
        ok = ok and tot == good == math.comb(n, t)
        verdict_obj = "graph state over F_%d" % F.Q
    else:
        raise ValueError(kind)

    # exact distance for the k-uniform (non-AME) ingredients: some (t+1)-set carries a stabilizer
    if t < n // 2:
        if kind == "fq_graph":
            tot = good = 0
            for S in itertools.combinations(range(n), t + 1):
                Sc = [j for j in range(n) if j not in S]
                tot += 1
                good += rank([[Gam[i][j] for j in Sc] for i in S], F) == t + 1
            lines.append("  exactness: |A| = %d: %d sets, %d of full row rank, %d not => not %d-uniform, distance exactly %d"
                         % (t + 1, tot, good, tot - good, t + 1, t + 1))
        else:
            tot, good = count_full_rank_column_sets(G, F, n - t - 1, dim)
            lines.append("  exactness: column sets of size %d: %d checked, %d of rank %d, %d rank-deficient => d(C) = %d"
                         % (n - t - 1, tot, good, dim, tot - good, t + 1))
        ok = ok and good < tot and fl["code"]["distance"] == t + 1
    elif kind != "fq_graph":
        ok = ok and fl["code"]["distance"] == n // 2 + 1   # Singleton bound: MDS distance is exact
    claim = fl["claimed_property"]["statement"]
    lines.append("  VERDICT: %s  (%s => %s, %s)" % ("PASS" if ok else "FAIL", verdict_obj, fl["stabilizer_code"], claim))
    return ok, lines


def _dot(u, v, F):
    s = 0
    for x, y in zip(u, v):
        s = F.add_t[s][F.mul_t[x][y]]
    return s


def main(argv):
    ids = argv[1:] or sorted(f.stem for f in HERE.glob("*.json"))
    print("check_fieldlevel.py: field-level verification of the coupling ingredients (exact ranks over F_Q)")
    verdict = {}
    for i in ids:
        fl = json.loads((HERE / (i + ".json")).read_text())
        fd = fl["field"]
        par = fl["code"]["params"]
        print("%s  [%s, F_%d%s, %s, claim: %s]" % (i, fl["kind"], fd["Q"], "" if fd["m"] == 1 else " mod " + fd["modulus_str"],
                                                  par, fl["claimed_property"]["statement"]))
        ok, lines = check_one(fl)
        for ln in lines:
            print(ln)
        verdict[i] = (ok, fl)
    npass = sum(ok for ok, _ in verdict.values())
    print("ingredients: %d PASS, %d FAIL" % (npass, len(verdict) - npass))

    # negative controls: altered data and over-claims must FAIL (the checks are not vacuous)
    nctrl = 0
    if len(argv) == 1:
        print()
        print("negative controls (expected FAIL)")
        controls = []
        fl = json.loads((HERE / "b_ame12_7_circulant.json").read_text())
        fl["A"][0][0] = (fl["A"][0][0] + 1) % fl["field"]["Q"]
        fl["structure"]["first_row"][0] = fl["A"][0][0]
        fl["A"] = circulant(fl["structure"]["first_row"])
        controls.append(("b_ame12_7_circulant, first-row entry 0 changed by +1 (as an int)", fl))
        fl = json.loads((HERE / "a_14_0_7_d4_hermitian_circulant.json").read_text())
        fl["claimed_property"]["uniformity_k"] = 7
        fl["code"]["distance"] = 8
        controls.append(("a_14_0_7_d4_hermitian_circulant claimed 7-uniform (AME(14,4))", fl))
        fl = json.loads((HERE / "a_dodecacode_n12_d2.json").read_text())
        fl["claimed_property"]["uniformity_k"] = 6
        controls.append(("a_dodecacode_n12_d2 claimed 6-uniform (AME(12,2))", fl))
        fl = json.loads((HERE / "b_ame16_9_hz.json").read_text())
        fl["structure"]["core_exponents"][0] += 1
        controls.append(("b_ame16_9_hz with core exponent 29 -> 30 (matrix as printed)", fl))
        for label, fl in controls:
            ok, lines = check_one(fl)
            nctrl += not ok
            failed = [ln.strip() for ln in lines if ("NO" in ln or "DIFFERENT" in ln or " checked, " in ln)]
            print("%-66s %s" % (label, "FAIL (as expected)" if not ok else "PASS (CONTROL FAILED)"))
            for ln in failed:
                print("    " + ln)
        allctrl = nctrl == len(controls)
    else:
        allctrl = True

    allok = npass == len(verdict) and allctrl
    cpath = HERE.parent / "cells.json"
    if len(argv) == 1 and cpath.exists():
        print()
        print("cells (coupling lemma: G_a k-uniform on 2k+2 parties, G_b AME(2k+2,b), a <= b => AME(2k+1, ab))")
        cells = json.loads(cpath.read_text())["cells"]
        nok = 0
        for c in cells:
            oa, fa = verdict[c["a"]]
            ob, fb = verdict[c["b"]]
            k = c["k"]
            a, b = fa["local_dim"], fb["local_dim"]
            ta, tb = fa["claimed_property"]["uniformity_k"], fb["claimed_property"]["uniformity_k"]
            cond = (c["n"] == 2 * k + 1 and fa["n"] == fb["n"] == 2 * k + 2 and ta >= k and tb == k + 1 and a <= b
                    and a * b == c["d"])
            ok = oa and ob and cond
            nok += ok
            print("%-11s %-17s G_a=%-32s %d-uniform %s | G_b=%-23s AME(%d,%d) %s | a=%d, b=%d, k=%d | %s"
                  % (c["cell"], "(" + c["verdict_at_submission"].lower() + ")", c["a"], ta, "PASS" if oa else "FAIL",
                     c["b"], fb["n"], b, "PASS" if ob else "FAIL", a, b, k, "OK" if ok else "FAIL"))
        print("cells: %d OK, %d FAIL" % (nok, len(cells) - nok))
        allok = allok and nok == len(cells)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
