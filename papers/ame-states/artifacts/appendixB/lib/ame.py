"""Helpers for p10 (open AME cases). Not trusted: the frozen verifier checks only the canonical form.

Canonical stabilizer form (what verifier/verify.py checks): N = m*n qudits of prime dimension p, party j owns qudits
j*m .. j*m+m-1, d = p^m; N rows (a | b) in F_p^(2N) for X(a)Z(b), isotropic for a.b' - b.a', rank N.
AME(n,d) <=> every set of ceil(n/2) parties has a full-rank (N x 2m*ceil(n/2)) column block.

Contents
  GF(q)                         finite field q = p^m <= 256 (ints 0..q-1 = polynomial-basis digits base p)
  canonical_from_fq(S, F)       F_q-linear symplectic code (n rows (a|b) over GF(q)) -> canonical F_p rows
  canonical_from_hermitian(G, F2, q)
                                F_{q^2}-linear Hermitian self-orthogonal code (n/2 x n over GF(q^2)) -> canonical rows
  canonical_from_graph(Gam, p)  qudit graph state over F_p (any grouping m): rows [e_i | Gam_i]
  mds_state_rows(G, F)          sum_{c in C}|c> for C = rowspace(G) over GF(q) (k = floor(n/2)); AME iff C is MDS
  rs_extended(F, n, k)          doubly-extended Reed-Solomon generator (n <= q+1)
  project_last_party(rows, p, m, n)   AME(2k,d) -> AME(2k-1,d) (projection onto |0..0> of the last party)
  bad_cuts(rows, p, m, n)       exact count via the verifier's compiled rankcheck (fast; same as the board)
  hermitian_self_orthogonal(G, F2, q), is_mds_systematic(A, F)   quick checks for code searches
  graph_state_phases(Gam, p)    qudit graph state amplitudes zeta_p^(sum_{i<j} G_ij x_i x_j + sum_i G_ii x_i^2/2 ...)
                                (p = 2: (-1)^(x^T G_upper x); odd p: zeta_p^(x^T G x / 2)), for the state type
  write_stabilizer(path, ...), write_state(path, ...)   candidate JSON writers
"""
import itertools
import json
import math
import os
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
# release/appendixB: rankcheck.c is compiled by run.sh into build/ (override with $AME_RANKCHECK)
RANKCHECK = pathlib.Path(os.environ.get("AME_RANKCHECK", str(ROOT / "build" / "rankcheck")))
DIG = "0123456789abcdefghijklmnopqrstuvwxyz"


def _is_prime(x):
    return x >= 2 and all(x % k for k in range(2, int(x ** 0.5) + 1))


def _factor_prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            m, x = 0, q
            while x % p == 0:
                x //= p
                m += 1
            if x != 1:
                raise ValueError(f"{q} is not a prime power")
            return p, m
    raise ValueError(q)


class GF:
    """GF(p^m) with log/exp tables. Elements are ints 0..q-1; digit i (base p) = coefficient of alpha^i,
    where alpha is a root of the lexicographically first primitive polynomial of degree m over F_p.
    The prime field F_p is {0..p-1}."""

    def __init__(self, q):
        self.q = q
        self.p, self.m = _factor_prime_power(q)
        p, m = self.p, self.m
        if q > 256:
            raise ValueError("GF(q): q <= 256 (add table is q x q)")
        if m == 1:
            self.poly = None
            g = 1 if p == 2 else next(g for g in range(2, p) if all(pow(g, (p - 1) // f, p) != 1 for f in self._pf(p - 1)))
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

    @staticmethod
    def _pf(x):
        f, k = [], 2
        while k * k <= x:
            if x % k == 0:
                f.append(k)
                while x % k == 0:
                    x //= k
            k += 1
        return f + ([x] if x > 1 else [])

    def _digits(self, a):
        return [(a // self.p ** i) % self.p for i in range(self.m)]

    def _from_digits(self, ds):
        return sum((d % self.p) * self.p ** i for i, d in enumerate(ds))

    def _add_slow(self, a, b):
        return self._from_digits([x + y for x, y in zip(self._digits(a), self._digits(b))])

    def _mulx(self, a):
        ds = [0] + self._digits(a)
        top = ds.pop()
        # alpha^m = -sum c_i alpha^i, poly = x^m + c_{m-1} x^{m-1} + ... + c_0
        return self._from_digits([d - top * c for d, c in zip(ds, self.poly)])

    def _primitive_poly(self):
        p, m, q = self.p, self.m, self.q
        for code in range(p ** m):
            cs = [(code // p ** i) % p for i in range(m)]
            if cs[0] == 0:
                continue
            self.poly = cs
            x, seen = 1, set()
            ok = True
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

    def frob(self, a, k=1):
        return self.pow(a, self.p ** k)

    def trace_to_prime(self, a, levels=None):
        """Tr_{GF(p^levels)/F_p}(a) for a in the subfield GF(p^levels) (default the whole field)."""
        levels = levels or self.m
        s, x = 0, a
        for _ in range(levels):
            s = self.add(s, x)
            x = self.frob(x)
        return s

    def subfield(self, q_sub):
        """Elements of the subfield GF(q_sub) (q_sub = p^k, k | m), as ints of this field."""
        return [a for a in range(self.q) if self.pow(a, q_sub) == a]

    def dot(self, u, v):
        s = 0
        for x, y in zip(u, v):
            if x and y:
                s = self.add(s, self.mul(x, y))
        return s


def _solve_coords(F, basis, levels):
    """Map every element of span_Fp(basis) to its coordinate vector (dict)."""
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
    """A basis beta of GF(q_sub) over F_p (inside F) and its trace-dual basis gamma (Tr(beta_i gamma_j) = delta_ij)."""
    p = F.p
    k = _factor_prime_power(q_sub)[1]
    sub = F.subfield(q_sub)
    # beta: powers of a generator of the subfield
    gen = next(a for a in sub if a and all(F.pow(a, (q_sub - 1) // f) != 1 for f in GF._pf(q_sub - 1))) if q_sub > 2 else 1
    beta = [F.pow(gen, i) for i in range(k)]
    T = [[F.trace_to_prime(F.mul(bi, bj), k) for bj in beta] for bi in beta]
    # invert T over F_p
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
    """F_q-linear stabilizer code: n rows (a | b) over GF(q) (as elements of F, which may be GF(q) itself or an
    extension containing it). Returns the m*n canonical rows over F_p (X in basis beta, Z in the trace-dual gamma)."""
    q = q or F.q
    beta, gamma = subfield_bases(F, q)
    cx, cz = _solve_coords(F, beta, None), _solve_coords(F, gamma, None)
    n = len(S[0]) // 2
    out = []
    for row in S:
        for bl in beta:
            r = [F.mul(bl, x) for x in row]
            out.append([c for x in r[:n] for c in cx[x]] + [c for x in r[n:] for c in cz[x]])
    return out


def canonical_from_hermitian(G, F2, q):
    """Hermitian self-orthogonal F_{q^2}-linear code (rows of G over F2 = GF(q^2)) -> canonical rows.
    Uses v = beta*a + beta^q*b with {beta, beta^q} a normal basis of F_{q^2}/F_q; then the trace-alternating form
    equals tr_{q/p}(b.a' - a.b'), so the image is a symplectic self-orthogonal F_q-linear code."""
    sub = F2.subfield(q)
    beta = next(b for b in range(1, F2.q) if len({F2.add(F2.mul(x, b), F2.mul(y, F2.pow(b, q))) for x in sub for y in sub}) == q * q)
    bq = F2.pow(beta, q)
    ab = {F2.add(F2.mul(x, beta), F2.mul(y, bq)): (x, y) for x in sub for y in sub}
    xi = next(x for x in range(F2.q) if x not in sub)  # F_q-basis {1, xi} of F_{q^2}
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
    """Generator of the [n, k, n-k+1]_q doubly-extended RS code over F (n <= q+1)."""
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


def _nullspace(G, F):
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
    """sum_{c in rowspace(G)} |c> over GF(q): X(c) for c in C, Z(h) for h in C^perp. Canonical rows."""
    n = len(G[0])
    H = _nullspace(G, F)
    S = [list(g) + [0] * n for g in G] + [[0] * n + list(h) for h in H]
    return canonical_from_fq(S, F)


def hermitian_self_orthogonal(G, F2, q):
    return all(F2.dot(u, [F2.pow(x, q) for x in v]) == 0 for u in G for v in G)


def is_mds_systematic(A, F):
    """[I | A] is MDS iff every square submatrix of A is nonsingular (brute force; fine for A up to ~8x8)."""
    k, r = len(A), len(A[0])
    for s in range(1, min(k, r) + 1):
        for R in itertools.combinations(range(k), s):
            for C in itertools.combinations(range(r), s):
                if _det([[A[i][j] for j in C] for i in R], F) == 0:
                    return False
    return True


def _det(M, F):
    M = [r[:] for r in M]
    n, det = len(M), 1
    for c in range(n):
        pr = next((i for i in range(c, n) if M[i][c]), None)
        if pr is None:
            return 0
        if pr != c:
            M[c], M[pr] = M[pr], M[c]
            det = F.neg(det)
        det = F.mul(det, M[c][c])
        iv = F.inv(M[c][c])
        for i in range(c + 1, n):
            if M[i][c]:
                f = F.mul(M[i][c], iv)
                M[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(M[i], M[c])]
    return det


def bad_cuts(rows, p, m, n):
    """Exact (bad_cuts, cuts) with the verifier's compiled rankcheck. bad_cuts = 0 <=> AME (rows must be a valid
    isotropic rank-N set; check that separately or just run the verifier)."""
    N, t = m * n, (n + 1) // 2
    inp = f"{p} {N} {n} {m} {t}\n" + "\n".join(" ".join(map(str, r)) for r in rows) + "\n"
    pr = subprocess.run([str(RANKCHECK)], input=inp, capture_output=True, text=True, check=True)
    head = pr.stdout.split()
    nbad = int(head[3])
    return (nbad // 2 if n % 2 == 0 else nbad), (math.comb(n, t) // 2 if n % 2 == 0 else math.comb(n, t))


def graph_state_phases(Gam, p):
    """Amplitude exponents (root_order 2 for p = 2, else p) of the qudit graph state on N = len(Gam) qudits,
    qudit 0 most significant; party grouping is by consecutive qudits, so the index matches the state type
    with d = p^m. For p = 2 off-diagonal entries only; for odd p, exponent = sum_{i<j} G_ij x_i x_j."""
    N = len(Gam)
    out = []
    for x in range(p ** N):
        digs = [(x // p ** (N - 1 - i)) % p for i in range(N)]
        e = sum(Gam[i][j] * digs[i] * digs[j] for i in range(N) for j in range(i + 1, N))
        out.append(e % p)
    return out


def project_last_party(rows, p, m, n):
    """Project party n-1 onto |0...0> (Z eigenstate): keep {lam.S : X part of lam.S on the last party = 0} and delete
    that party's columns. AME(2k,d) -> AME(2k-1,d). Returns canonical rows for n-1 parties (N - m rows)."""
    N = m * n
    xcols = [(n - 1) * m + k for k in range(m)]
    # nullspace of lam -> lam . S[:, xcols] over F_p
    M = [[rows[i][c] % p for i in range(N)] for c in xcols]  # m x N
    piv, r = [], 0
    M = [row[:] for row in M]
    for c in range(N):
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
    lams = []
    for f in (j for j in range(N) if j not in piv):
        lam = [0] * N
        lam[f] = 1
        for i, pc in enumerate(piv):
            lam[pc] = (-M[i][f]) % p
        lams.append(lam)
    keep = [k for k in range(2 * N) if not ((n - 1) * m <= k < N or N + (n - 1) * m <= k < 2 * N)]
    out = []
    for lam in lams:
        v = [sum(l * rows[i][k] for i, l in enumerate(lam) if l) % p for k in range(2 * N)]
        out.append([v[k] for k in keep])
    return out


def write_stabilizer(path, n, d, p, rows, script="", note=""):
    rows = ["".join(DIG[v] for v in r) for r in rows] if p <= 36 else rows
    obj = {"type": "stabilizer", "n": n, "d": d, "p": p, "generators": rows, "script": script, "note": note}
    pathlib.Path(path).write_text(json.dumps(obj) + "\n")
    return path


def write_state(path, n, d, root_order, norm, phases=None, amplitudes=None, script="", note=""):
    obj = {"type": "state", "n": n, "d": d, "root_order": root_order, "norm": norm, "script": script, "note": note}
    if phases is not None:
        obj["phases"] = phases
    else:
        obj["amplitudes"] = amplitudes
    pathlib.Path(path).write_text(json.dumps(obj) + "\n")
    return path
