#!/usr/bin/env python3
"""Field-level description -> canonical F_p stabilizer rows, and comparison with the certified rows.

For every fieldlevel/<id>.json this script rebuilds the canonical stabilizer generator rows (the form the frozen
verifier checked) from the field-level object, using the documented converter below, and compares the F_p row space
with the certified rows in ../ingredients/<id>.json (equal reduced row echelon forms over F_p => MATCH).
It then checks the ingredient-to-cell mapping: for every cell in ../cells.json, the a- and b-ingredients of the frozen
certificate ../certificates/AME_n_d.json span the same F_p row spaces as the ingredient files named in cells.json.

Self-contained: Python 3.9 standard library only. The finite-field code and the converters are copied (with minor
adaptations, no behavioural change) from the repository's qlib/ame.py, which generated the certified rows.

Conventions (identical to qlib/ame.py)
  GF(Q), Q = p^m: elements are ints 0..Q-1; base-p digit i of the int is the coefficient of x^i, where x is a root of
  the LEXICOGRAPHICALLY FIRST PRIMITIVE monic polynomial of degree m over F_p (coefficient vector (c_0..c_{m-1}) read
  as the base-p number sum c_i p^i, smallest first, c_0 != 0). For m = 1 this is F_p with the usual integers.
  Canonical form: n parties, local dimension d = p^m_loc, N = m_loc * n qudits of dimension p, party j owns qudits
  j*m_loc .. j*m_loc + m_loc - 1; N rows (a | b) in F_p^(2N) for X(a) Z(b).
  canonical_from_fq(S, F, q): S = rows (a | b) over F_q (F_q inside F). For each row and each beta_l of the basis
      beta = (1, g, .., g^(k-1)) of F_q/F_p (g = smallest-int primitive element of the subfield F_q, k = log_p q;
      beta = (1) if q is prime), the row beta_l * (a | b) is written with X coordinates in the basis beta and Z
      coordinates in the trace-dual basis gamma (Tr_{F_q/F_p}(beta_i gamma_j) = delta_ij).
  canonical_from_hermitian(G, F2, q): G (rows over F2 = GF(q^2)) Hermitian self-orthogonal. With w = the smallest
      int such that {w, w^q} is an F_q-basis of F_{q^2} (normal basis) and xi = the smallest int outside F_q, each
      row g gives the two F_q-rows (a | b) of lam * g, lam in (1, xi), where lam*g = w*a + w^q*b coordinatewise;
      then canonical_from_fq(S, F2, q).
  mds_state_rows(G, F): state sum_{c in rowspace G}|c>: S = [G | 0] (X part) + [0 | H] (Z part), H = basis of the
      Euclidean dual C^perp (standard nullspace from the RREF of G); then canonical_from_fq(S, F, F.q).
  canonical_from_graph(Gamma, p): rows [e_i | Gamma_i] over F_p (graph state, one qudit per party here).

usage: python3 convert.py [id ...]     (default: every fieldlevel/*.json)
"""
import itertools
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
COUPLING = HERE.parent
DIG = "0123456789abcdefghijklmnopqrstuvwxyz"


# ----------------------------------------------------------------------------------------------------------------
# finite fields (copy of qlib/ame.py GF)
# ----------------------------------------------------------------------------------------------------------------
def _factor_prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            m, x = 0, q
            while x % p == 0:
                x //= p
                m += 1
            if x != 1:
                raise ValueError("%d is not a prime power" % q)
            return p, m
    raise ValueError(q)


def _prime_factors(x):
    f, k = [], 2
    while k * k <= x:
        if x % k == 0:
            f.append(k)
            while x % k == 0:
                x //= k
        k += 1
    return f + ([x] if x > 1 else [])


class GF:
    """GF(p^m) as in qlib: ints 0..q-1, digit i (base p) = coefficient of x^i, x a root of the lexicographically first
    primitive polynomial poly = x^m + c_{m-1} x^{m-1} + .. + c_0 (self.poly = [c_0, .., c_{m-1}])."""

    def __init__(self, q):
        self.q = q
        self.p, self.m = _factor_prime_power(q)
        p, m = self.p, self.m
        if m == 1:
            self.poly = None
            g = 1 if p == 2 else next(g for g in range(2, p) if all(pow(g, (p - 1) // f, p) != 1 for f in _prime_factors(p - 1)))
            self.exp = [pow(g, i, p) for i in range(p - 1)]
        else:
            self.poly = self._primitive_poly()
            self.exp = []
            x = 1
            for _ in range(q - 1):
                self.exp.append(x)
                x = self._mulx(x)
        self.log = {v: i for i, v in enumerate(self.exp)}
        if len(self.log) != q - 1:
            raise RuntimeError("not primitive")
        self._add = [[self._add_slow(a, b) for b in range(q)] for a in range(q)]
        self._neg = [next(b for b in range(q) if self._add[a][b] == 0) for a in range(q)]

    def modulus(self):
        """monic modulus coefficients [c_0, .., c_{m-1}, 1] (None for a prime field)"""
        return None if self.poly is None else list(self.poly) + [1]

    def _digits(self, a):
        return [(a // self.p ** i) % self.p for i in range(self.m)]

    def _from_digits(self, ds):
        return sum((d % self.p) * self.p ** i for i, d in enumerate(ds))

    def _add_slow(self, a, b):
        return self._from_digits([x + y for x, y in zip(self._digits(a), self._digits(b))])

    def _mulx(self, a):
        ds = [0] + self._digits(a)
        top = ds.pop()
        return self._from_digits([d - top * c for d, c in zip(ds, self.poly)])

    def _primitive_poly(self):
        p, m, q = self.p, self.m, self.q
        for code in range(p ** m):
            cs = [(code // p ** i) % p for i in range(m)]
            if cs[0] == 0:
                continue
            self.poly = cs
            x, seen, ok = 1, set(), True
            for _ in range(q - 1):
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = self._mulx(x)
            if ok and x == 1 and len(seen) == q - 1:
                return cs
        raise RuntimeError("no primitive polynomial")

    def add(self, a, b):
        return self._add[a][b]

    def sub(self, a, b):
        return self._add[a][self._neg[b]]

    def neg(self, a):
        return self._neg[a]

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[(self.log[a] + self.log[b]) % (self.q - 1)]

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError
        return self.exp[(-self.log[a]) % (self.q - 1)]

    def pow(self, a, e):
        if a == 0:
            return 0 if e else 1
        return self.exp[(self.log[a] * e) % (self.q - 1)]

    def trace_to_prime(self, a, levels=None):
        levels = levels or self.m
        s, x = 0, a
        for _ in range(levels):
            s = self.add(s, x)
            x = self.pow(x, self.p)
        return s

    def subfield(self, q_sub):
        return [a for a in range(self.q) if self.pow(a, q_sub) == a]

    def dot(self, u, v):
        s = 0
        for x, y in zip(u, v):
            if x and y:
                s = self.add(s, self.mul(x, y))
        return s


_FIELDS = {}


def field(q):
    if q not in _FIELDS:
        _FIELDS[q] = GF(q)
    return _FIELDS[q]


# ----------------------------------------------------------------------------------------------------------------
# converters (copy of qlib/ame.py)
# ----------------------------------------------------------------------------------------------------------------
def _solve_coords(F, basis):
    p = F.p
    table = {}
    for coeffs in itertools.product(range(p), repeat=len(basis)):
        v = 0
        for c, b in zip(coeffs, basis):
            v = F.add(v, F.mul(c, b))
        table[v] = list(coeffs)
    if len(table) != p ** len(basis):
        raise ValueError("basis is not independent over F_p")
    return table


def subfield_bases(F, q_sub):
    """basis beta of GF(q_sub)/F_p (inside F) and its trace-dual basis gamma"""
    p = F.p
    k = _factor_prime_power(q_sub)[1]
    sub = F.subfield(q_sub)
    gen = next(a for a in sub if a and all(F.pow(a, (q_sub - 1) // f) != 1 for f in _prime_factors(q_sub - 1))) if q_sub > 2 else 1
    beta = [F.pow(gen, i) for i in range(k)]
    T = [[F.trace_to_prime(F.mul(bi, bj), k) for bj in beta] for bi in beta]
    M = [row[:] + [1 if i == j else 0 for j in range(k)] for i, row in enumerate(T)]
    for c in range(k):
        piv = next(i for i in range(c, k) if M[i][c] % p)
        M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p - 2, p)
        M[c] = [x * iv % p for x in M[c]]
        for i in range(k):
            if i != c and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    Tinv = [row[k:] for row in M]
    gamma = []
    for j in range(k):
        g = 0
        for i in range(k):
            g = F.add(g, F.mul(Tinv[i][j], beta[i]))
        gamma.append(g)
    for i in range(k):
        for j in range(k):
            assert F.trace_to_prime(F.mul(beta[i], gamma[j]), k) == (1 if i == j else 0)
    return beta, gamma


def canonical_from_fq(S, F, q=None):
    q = q or F.q
    beta, gamma = subfield_bases(F, q)
    cx, cz = _solve_coords(F, beta), _solve_coords(F, gamma)
    n = len(S[0]) // 2
    out = []
    for row in S:
        for bl in beta:
            r = [F.mul(bl, x) for x in row]
            out.append([c for x in r[:n] for c in cx[x]] + [c for x in r[n:] for c in cz[x]])
    return out


def normal_element(F2, q):
    """smallest int w with {w, w^q} an F_q-basis of F_{q^2} (qlib's choice)"""
    sub = F2.subfield(q)
    return next(b for b in range(1, F2.q)
                if len({F2.add(F2.mul(x, b), F2.mul(y, F2.pow(b, q))) for x in sub for y in sub}) == q * q)


def canonical_from_hermitian(G, F2, q):
    sub = F2.subfield(q)
    beta = normal_element(F2, q)
    bq = F2.pow(beta, q)
    ab = {F2.add(F2.mul(x, beta), F2.mul(y, bq)): (x, y) for x in sub for y in sub}
    xi = next(x for x in range(F2.q) if x not in sub)
    S = []
    for g in G:
        for lam in (1, xi):
            v = [F2.mul(lam, x) for x in g]
            S.append([ab[x][0] for x in v] + [ab[x][1] for x in v])
    return canonical_from_fq(S, F2, q)


def canonical_from_graph(Gam, p):
    N = len(Gam)
    return [[1 if k == i else 0 for k in range(N)] + [Gam[i][k] % p for k in range(N)] for i in range(N)]


def rs_extended(F, n, k):
    """qlib's [n, k, n-k+1] doubly-extended RS generator: points = the ints 0..min(n,q)-1 as field elements,
    row i = (x^i)_x with 0^0 = 1, plus the column e_{k-1} (point at infinity) when n = q+1"""
    q = F.q
    if n > q + 1:
        raise ValueError("n <= q+1")
    pts = list(range(min(n, q)))
    G = []
    for i in range(k):
        row = [F.pow(x, i) if x else (1 if i == 0 else 0) for x in pts]
        if n == q + 1:
            row.append(1 if i == k - 1 else 0)
        G.append(row)
    return G


def nullspace(G, F):
    k, n = len(G), len(G[0])
    M = [r[:] for r in G]
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, k) if M[i][c]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = F.inv(M[r][c])
        M[r] = [F.mul(iv, x) for x in M[r]]
        for i in range(k):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    H = []
    for f in (j for j in range(n) if j not in piv):
        h = [0] * n
        h[f] = 1
        for i, pc in enumerate(piv):
            h[pc] = F.neg(M[i][f])
        H.append(h)
    return H


def mds_state_rows(G, F):
    n = len(G[0])
    H = nullspace(G, F)
    S = [list(g) + [0] * n for g in G] + [[0] * n + list(h) for h in H]
    return canonical_from_fq(S, F)


# ----------------------------------------------------------------------------------------------------------------
# F_p linear algebra
# ----------------------------------------------------------------------------------------------------------------
def rref_mod_p(rows, p):
    M = [[x % p for x in r] for r in rows]
    piv, r = [], 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = pow(M[r][c], p - 2, p)
        M[r] = [x * iv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    return M[:r], piv


def parse_rows(obj):
    gens = obj["generators"]
    if gens and isinstance(gens[0], str):
        return [[DIG.index(ch) for ch in r] for r in gens]
    return [list(map(int, r)) for r in gens]


def symplectic_isotropic(rows, p):
    N = len(rows[0]) // 2
    for u in rows:
        for v in rows:
            s = sum(u[i] * v[N + i] - u[N + i] * v[i] for i in range(N)) % p
            if s:
                return False
    return True


# ----------------------------------------------------------------------------------------------------------------
# field-level JSON -> canonical rows
# ----------------------------------------------------------------------------------------------------------------
def check_field_convention(fl, F):
    """the JSON's stated field must be qlib's field (same modulus), otherwise the integer labels mean something else"""
    fd = fl["field"]
    assert fd["Q"] == F.q and fd["p"] == F.p and fd["m"] == F.m, (fd, F.q)
    assert fd["modulus"] == F.modulus(), ("modulus differs from qlib's lexicographically first primitive polynomial",
                                          fd["modulus"], F.modulus())


def to_canonical(fl):
    """returns (p, m_loc, rows) for the field-level object fl"""
    kind = fl["kind"]
    fd = fl["field"]
    F = field(fd["Q"])
    check_field_convention(fl, F)
    n = fl["n"]
    cv = fl.get("converter", {})
    if kind == "hermitian_systematic":
        q = fl["q"]
        k = len(fl["A"])
        G = [[1 if i == j else 0 for j in range(k)] + list(fl["A"][i]) for i in range(k)]
        rows = canonical_from_hermitian(G, F, q)
        p, m_loc = _factor_prime_power(q)
        # the converter constants printed in the JSON are the ones used here
        beta, gamma = subfield_bases(F, q)
        w = normal_element(F, q)
        xi = next(x for x in range(F.q) if x not in set(F.subfield(q)))
        assert (cv["subfield_basis_beta"], cv["trace_dual_basis_gamma"], cv["normal_element_w"], cv["w_to_the_q"],
                cv["xi"]) == (beta, gamma, w, F.pow(w, q), xi), "converter constants in the JSON differ"
    elif kind == "css_minimal_support":
        rows = mds_state_rows([list(r) for r in fl["G"]], F)
        p, m_loc = F.p, F.m
        assert (cv["subfield_basis_beta"], cv["trace_dual_basis_gamma"]) == subfield_bases(F, F.q), \
            "converter constants in the JSON differ"
    elif kind == "fq_graph":
        assert F.m == 1, "graph converter implemented for prime fields (all graph ingredients here are qubit states)"
        rows = canonical_from_graph(fl["Gamma"], F.p)
        p, m_loc = F.p, 1
    else:
        raise ValueError(kind)
    assert len(rows) == m_loc * n
    return p, m_loc, rows


def compare(fl):
    ing = json.loads((COUPLING / "ingredients" / (fl["id"] + ".json")).read_text())
    p, m_loc, rows = to_canonical(fl)
    ok_meta = ing["n"] == fl["n"] and ing["d"] == p ** m_loc == fl["local_dim"] and ing["p"] == p
    cert = parse_rows(ing)
    R1, _ = rref_mod_p(rows, p)
    R2, _ = rref_mod_p(cert, p)
    N = m_loc * fl["n"]
    iso = symplectic_isotropic(rows, p)
    match = ok_meta and R1 == R2 and len(R1) == N and iso
    return match, dict(p=p, m=m_loc, N=N, rank_built=len(R1), rank_cert=len(R2), isotropic=iso, meta=ok_meta)


def main(argv):
    ids = argv[1:] or sorted(f.stem for f in HERE.glob("*.json"))
    print("convert.py: field-level object -> canonical F_p rows (qlib conventions) vs certified rows ../ingredients/<id>.json")
    nmatch = 0
    for i in ids:
        fl = json.loads((HERE / (i + ".json")).read_text())
        match, info = compare(fl)
        nmatch += match
        print("%-34s %-20s p=%-2d m=%d N=%-2d rank(built)=%-2d rank(certified)=%-2d isotropic=%s  %s" % (
            i, fl["kind"], info["p"], info["m"], info["N"], info["rank_built"], info["rank_cert"],
            info["isotropic"], "MATCH" if match else "MISMATCH"))
    print("ingredients: %d MATCH, %d MISMATCH" % (nmatch, len(ids) - nmatch))

    # ingredient-to-cell mapping: certificate a/b rows == ingredient rows (F_p row spaces)
    cells = json.loads((COUPLING / "cells.json").read_text())["cells"]
    print()
    print("cell mapping: frozen certificate's a/b == ingredient files named in cells.json (equal F_p row spaces)")
    rcache = {}

    def ing_rref(name):
        if name not in rcache:
            ob = json.loads((COUPLING / "ingredients" / (name + ".json")).read_text())
            rcache[name] = rref_mod_p(parse_rows(ob), ob["p"])[0]
        return rcache[name]

    ncell_ok = 0
    for c in cells:
        cert = json.loads((COUPLING / c["certificate"]).read_text())
        ok = cert["type"] == "coupling" and cert["k"] == c["k"]
        res = []
        for role in ("a", "b"):
            ob = cert[role]
            same = rref_mod_p(parse_rows(ob), ob["p"])[0] == ing_rref(c[role])
            ok = ok and same
            res.append("%s=%s %s" % (role, c[role], "same" if same else "DIFFERENT"))
        n_par = cert["a"]["n"]
        ok = ok and n_par == cert["b"]["n"] == 2 * c["k"] + 2 == c["n"] + 1 and cert["a"]["d"] * cert["b"]["d"] == c["d"] \
            and cert["a"]["d"] <= cert["b"]["d"]
        ncell_ok += ok
        print("%-11s k=%d %-10s %s  %s" % (c["cell"], c["k"], "(%dx%d)" % (cert["a"]["d"], cert["b"]["d"]), "  ".join(res),
                                           "OK" if ok else "FAIL"))
    print("cells: %d OK, %d FAIL" % (ncell_ok, len(cells) - ncell_ok))

    # negative controls: one altered entry must change the row space (the comparison is not vacuous)
    print()
    print("negative controls (one entry altered; expected MISMATCH)")
    nctrl = 0
    for cid, key in (("a_14_0_7_d4_hermitian_circulant", "A"), ("b_ame16_9_hz", "A"), ("b_ame12_16_rs", "G"),
                     ("a_ternary_golay_n12_d3", "G"), ("a_dodecacode_n12_d2", "Gamma")):
        path = HERE / (cid + ".json")
        if not path.exists():
            continue
        fl = json.loads(path.read_text())
        F = field(fl["field"]["Q"])
        M = [list(r) for r in fl[key]]
        if key == "Gamma":   # flip a symmetric pair (stays a graph state)
            M[0][1] = M[1][0] = (M[0][1] + 1) % F.p
        else:
            M[0][-1] = F.add(M[0][-1], 1)
        fl[key] = M
        match, info = compare(fl)
        nctrl += not match
        print("%-34s %s[0][%s] altered: %s" % (cid, key, "1" if key == "Gamma" else "-1",
                                               "MISMATCH (as expected)" if not match else "MATCH (CONTROL FAILED)"))
    return 0 if nmatch == len(ids) and ncell_ok == len(cells) and nctrl == 5 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
