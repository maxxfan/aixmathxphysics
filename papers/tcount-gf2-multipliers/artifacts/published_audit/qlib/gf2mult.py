"""GF(2^n) multiplication as an alternating trilinear form over GF(2), and Waring decompositions.

Variables: a_0..a_{n-1} -> 0..n-1, b_j -> n+j, c_k -> 2n+k (N = 3n).  The multiplier
|a,b,c> -> |a,b,c+ab> equals H_c D H_c with D = (-1)^{f}, f(a,b,c) = c . (a*b mod P).
Its signature tensor S is the zero-diagonal symmetric tensor with S_{ijk} = 1 exactly on the
triples (a_i, b_j, c_k) where x^k appears in x^{i+j} mod P.  A set of nonzero vectors
u_1..u_m in GF(2)^N with sum_r u_r^{x3} = S (mod 2) gives a Clifford+T implementation with
T-count m (Heyfron-Campbell signature-tensor theorem); m is the quantity to minimize.
Factors are Python ints (bit i = variable i).  Stdlib only.
"""
import itertools


def poly_mod(a, p):
    dp = p.bit_length() - 1
    while a and a.bit_length() - 1 >= dp:
        a ^= p << (a.bit_length() - 1 - dp)
    return a


def is_irreducible(p):
    n = p.bit_length() - 1
    if n < 1 or not p & 1:
        return n == 1
    for d in range(2, 1 << (n // 2 + 1)):
        if 1 <= d.bit_length() - 1 <= n // 2 and poly_mod(p, d) == 0:
            return False
    return True


def irreducibles(n):
    return [p for p in range(1 << n, 1 << (n + 1)) if is_irreducible(p)]


def default_modulus(n):
    """Lowest-weight, then smallest, irreducible polynomial of degree n (trinomial if one exists)."""
    return min(irreducibles(n), key=lambda p: (bin(p).count("1"), p))


def gf_mul(a, b, p):
    n = p.bit_length() - 1
    r = 0
    for i in range(n):
        if b >> i & 1:
            r ^= a << i
    return poly_mod(r, p)


def target_triples(n, p):
    """Set of sorted index triples (i, n+j, 2n+k) where the cubic form has a monomial."""
    S = set()
    for i in range(n):
        for j in range(n):
            v = gf_mul(1 << i, 1 << j, p)
            for k in range(n):
                if v >> k & 1:
                    S.add((i, n + j, 2 * n + k))
    return S


def odd_entries(factors):
    """Odd entries of sum u^{x3}, keyed by support: singles (i,), pairs (i,j), triples (i,j,k)."""
    odd = set()
    for u in factors:
        s = [i for i in range(u.bit_length()) if u >> i & 1]
        for k in (1, 2, 3):
            for t in itertools.combinations(s, k):
                if t in odd:
                    odd.remove(t)
                else:
                    odd.add(t)
    return odd


def cube_sum_matches(factors, n, p):
    """True iff sum u^{x3} equals the multiplier's signature tensor (all entries, incl. diagonal)."""
    return odd_entries(factors) == target_triples(n, p)


def triples_of(factors):
    return {t for t in odd_entries(factors) if len(t) == 3}


def apply_linear(factors, L):
    """Map each factor u to L^T-image: new_u_j = sum_i u_i L[i][j]  (L given as list of row ints)."""
    out = []
    for u in factors:
        v = 0
        for i, row in enumerate(L):
            if u >> i & 1:
                v ^= row
        out.append(v)
    return out


def slice_rank_profile(triples, N):
    """GL(N)-invariant fingerprint: histogram of rank of the alternating form S(u,.,.) over u."""
    hist = {}
    by_var = {i: [] for i in range(N)}
    for t in triples:
        for x in t:
            by_var[x].append(t)
    for u in range(1, 1 << N):
        B = [0] * N  # bilinear form rows
        for i in range(N):
            if u >> i & 1:
                for (x, y, z) in by_var[i]:
                    a, b = [w for w in (x, y, z) if w != i]
                    B[a] ^= 1 << b
                    B[b] ^= 1 << a
        r = gf2_rank(B)
        hist[r] = hist.get(r, 0) + 1
    return dict(sorted(hist.items()))


def gf2_rank(rows):
    rows = [r for r in rows if r]
    rank = 0
    while rows:
        piv = max(rows)
        hb = piv.bit_length() - 1
        rows = [r ^ piv if r >> hb & 1 else r for r in rows if r != piv]
        rows = [r for r in rows if r]
        rank += 1
    return rank
