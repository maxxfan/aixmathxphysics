"""Map a decomposition of an *equivalent* tensor (e.g. AlphaTensor-Quantum's, in its circuit's own
qubit basis) onto the standard GF(2^n) multiplier target of gf2mult.py.

Method (all GF(2) linear algebra; the result is self-checked by gf2mult.cube_sum_matches):
 1. Registers.  For the standard form the a-, b- and c-registers are n-dim subspaces U with
    S(U,U,.) = 0; their nonzero points x are exactly those where the alternating form
    B_x = S(x,.,.) has the minimal rank 2n.  A random (2n+1)-dim subspace meets every n-dim
    subspace of GF(2)^{3n}, so enumerating it finds points of all three registers; the
    register through x is the radical of B_x.
 2. Isotopy.  In register coordinates the form is t(a,b,c) = <M(a,b), c>.  Kaplansky's trick
    (u o v = M(R^-1 u, L^-1 v) with R = M(., e_b), L = M(e_a, .)) gives a unital algebra on
    C*, which is the field; a generator's minimal polynomial is P and its powers give the
    polynomial basis.
 3. Factors u (linear forms u.x) transform as u -> M^T u under x_old = M x_new.
"""
import random
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gf2mult as G  # noqa: E402


# ---------------------------------------------------------------- GF(2) linear algebra (int rows)

def rank(rows):
    return G.gf2_rank(list(rows))


def solve_basis(vectors):
    """Echelon basis (list of ints) of span(vectors)."""
    basis = []
    for v in vectors:
        for b in basis:
            if v ^ b < v:
                v ^= b
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return basis


def kernel(rows, ncols):
    """Basis of {x : <row, x> = 0 for all rows} in GF(2)^ncols."""
    piv = []
    R = []
    for r in rows:
        for (p, rr) in zip(piv, R):
            if r >> p & 1:
                r ^= rr
        if r:
            p = r.bit_length() - 1
            for i in range(len(R)):
                if R[i] >> p & 1:
                    R[i] ^= r
            piv.append(p)
            R.append(r)
    free = [c for c in range(ncols) if c not in set(piv)]
    out = []
    for f in free:
        x = 1 << f
        for (p, rr) in zip(piv, R):
            if rr >> f & 1:
                x |= 1 << p
        out.append(x)
    return out


def mat_inverse(cols, n):
    """Inverse of the n x n matrix whose column j is the int cols[j]; returns columns."""
    aug = [(sum(((cols[j] >> i) & 1) << j for j in range(n))) | (1 << (n + i)) for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if aug[r] >> c & 1), None)
        if p is None:
            raise ValueError("singular")
        aug[c], aug[p] = aug[p], aug[c]
        for r in range(n):
            if r != c and aug[r] >> c & 1:
                aug[r] ^= aug[c]
    rows_inv = [aug[i] >> n for i in range(n)]  # row i of inverse
    return [sum(((rows_inv[i] >> j) & 1) << i for i in range(n)) for j in range(n)]


def mat_vec(cols, v):
    out = 0
    j = 0
    while v:
        if v & 1:
            out ^= cols[j]
        v >>= 1
        j += 1
    return out


# ---------------------------------------------------------------- tensor helpers

def slice_rows(triples, N):
    """S_gamma as alternating matrices: S[g][a] = bitmask of b with {a,b,g} a triple."""
    S = [[0] * N for _ in range(N)]
    for t in triples:
        for g in t:
            a, b = [w for w in t if w != g]
            S[g][a] ^= 1 << b
            S[g][b] ^= 1 << a
    return S


def form_B(S, x, N):
    B = [0] * N
    g = 0
    while x:
        if x & 1:
            for a in range(N):
                B[a] ^= S[g][a]
        x >>= 1
        g += 1
    return B


def trilinear(triples, x, y, z):
    s = 0
    for (i, j, k) in triples:
        xi, xj, xk = x >> i & 1, x >> j & 1, x >> k & 1
        yi, yj, yk = y >> i & 1, y >> j & 1, y >> k & 1
        zi, zj, zk = z >> i & 1, z >> j & 1, z >> k & 1
        s ^= (xi & yj & zk) ^ (xi & yk & zj) ^ (xj & yi & zk) ^ (xj & yk & zi) ^ (xk & yi & zj) ^ (xk & yj & zi)
    return s


def registers(triples, n, rng):
    N = 3 * n
    S = slice_rows(triples, N)
    for _attempt in range(20):
        W = solve_basis([rng.getrandbits(N) for _ in range(2 * n + 1)])
        if len(W) < 2 * n + 1:
            continue
        found = {}
        x, Bx = 0, [0] * N
        for i in range(1, 1 << len(W)):  # Gray code over W
            bit = (i & -i).bit_length() - 1
            x ^= W[bit]
            Bw = form_B(S, W[bit], N)
            Bx = [p ^ q for p, q in zip(Bx, Bw)]
            if rank(Bx) == 2 * n:
                rad = tuple(sorted(solve_basis(kernel(Bx, N))))
                if len(rad) == n:
                    found[rad] = True
        regs = list(found)
        if len(regs) >= 3:
            # choose three registers spanning everything
            for i in range(len(regs)):
                for j in range(i + 1, len(regs)):
                    for k in range(j + 1, len(regs)):
                        if rank(list(regs[i]) + list(regs[j]) + list(regs[k])) == N:
                            return [list(regs[i]), list(regs[j]), list(regs[k])]
    raise RuntimeError("registers not found")


def standardize(factors, n, seed=1):
    """Return (P, mapped_factors) with mapped factors decomposing the standard target for P."""
    rng = random.Random(seed)
    N = 3 * n
    tri = G.triples_of(factors)
    if any(len(t) < 3 for t in G.odd_entries(factors)):
        raise ValueError("tensor has nonzero diagonal entries; not a pure cubic form")
    UA, UB, UC = registers(tri, n, rng)
    # dual pairing for c: coordinates of M(a,b) via functionals c -> t(a,b,c) on basis UC
    def M(a, b):  # vector in GF(2)^n indexed by UC basis
        return sum(trilinear(tri, a, b, UC[k]) << k for k in range(n))

    def comb(basis, coords):
        v = 0
        for k in range(n):
            if coords >> k & 1:
                v ^= basis[k]
        return v

    for _ in range(200):
        ea = comb(UA, rng.getrandbits(n) or 1)
        eb = comb(UB, rng.getrandbits(n) or 1)
        Rcols = [M(UA[i], eb) for i in range(n)]      # R: A-coords -> C*
        Lcols = [M(ea, UB[j]) for j in range(n)]      # L: B-coords -> C*
        try:
            Rinv, Linv = mat_inverse(Rcols, n), mat_inverse(Lcols, n)
        except ValueError:
            continue
        e = M(ea, eb)

        def mult(u, v):  # unital product on C* (n-bit vectors)
            return M(comb(UA, mat_vec(Rinv, u)), comb(UB, mat_vec(Linv, v)))

        # find generator g with 1, g, ..., g^{n-1} independent; P from g^n
        for _g in range(200):
            g = rng.getrandbits(n)
            pw = [e]
            for _k in range(n):
                pw.append(mult(pw[-1], g))
            if rank(pw[:n]) < n:
                continue
            Pcols = mat_inverse(pw[:n], n)        # C*-vector -> coords in basis {g^k}
            gn = mat_vec(Pcols, pw[n])            # g^n = sum c_k g^k
            P = (1 << n) | gn
            if not G.is_irreducible(P):
                continue
            # new coordinates: alpha = coords of R(a) in g-basis, beta = coords of L(b)
            # product alpha*beta mod P = coords of R(a) o L(b) = coords of M(a, b)
            # old point x = sum a_i UA_i + b_j UB_j + c_k UC_k  (register coordinates)
            # Build x_old = Mbig x_new with x_new = (alpha, beta, gamma), gamma dual to g-basis.
            Ra = [mat_vec(Pcols, Rcols[i]) for i in range(n)]   # alpha = Ra-matrix * a_coords
            Lb = [mat_vec(Pcols, Lcols[j]) for j in range(n)]
            A_of_alpha = mat_inverse(Ra, n)                      # a_coords = A_of_alpha * alpha
            B_of_beta = mat_inverse(Lb, n)
            # c: t = <M(a,b), c_coords>_UC = <Pinv * coords_g(M), c> ; with M = Pg * m_g (Pg = pw cols)
            # <Pg m, c> = <m, Pg^T c> -> gamma = Pg^T c  =>  c = (Pg^T)^{-1} gamma
            PgT_cols = [sum(((pw[r] >> k) & 1) << r for r in range(n)) for k in range(n)]
            C_of_gamma = mat_inverse(PgT_cols, n)
            cols = []
            for i in range(n):
                cols.append(comb(UA, A_of_alpha[i]))
            for j in range(n):
                cols.append(comb(UB, B_of_beta[j]))
            for k in range(n):
                cols.append(comb(UC, C_of_gamma[k]))
            # factor u (form u.x_old) -> u' with u'_m = u . cols[m]
            mapped = [sum(((bin(u & cols[m]).count("1")) & 1) << m for m in range(N)) for u in factors]
            if G.cube_sum_matches(mapped, n, P):
                return P, mapped
    raise RuntimeError("isotopy to a field not found")
