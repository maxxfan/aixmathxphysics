"""f2m.py : 4x4 matrices over F_2 packed in 16 bits (entry (r,c) = bit 4r+c), symplectic form J on (x1,x2,z1,z2).
Shared helpers for the frame classification of F_2-linear stabilizer AME(2m,4) (see runs/classify/PROOF.md)."""


def mul(a, b):
    """4x4 product over F_2 (packed)."""
    c = 0
    for r in range(4):
        row = (a >> (4 * r)) & 15
        acc = 0
        for k in range(4):
            if (row >> k) & 1:
                acc ^= (b >> (4 * k)) & 15
        c |= acc << (4 * r)
    return c


def tr(a):
    c = 0
    for r in range(4):
        for k in range(4):
            if (a >> (4 * r + k)) & 1:
                c |= 1 << (4 * k + r)
    return c


def rank(a):
    rows = [(a >> (4 * r)) & 15 for r in range(4)]
    rk = 0
    for bit in range(4):
        piv = next((i for i in range(rk, 4) if (rows[i] >> bit) & 1), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(4):
            if i != rk and (rows[i] >> bit) & 1:
                rows[i] ^= rows[rk]
        rk += 1
    return rk


I4 = sum(1 << (5 * i) for i in range(4))  # identity
# J[x1][z1]=J[x2][z2]=J[z1][x1]=J[z2][x2]=1 in coordinate order (x1,x2,z1,z2)
J = 0
for (r, c) in [(0, 2), (1, 3), (2, 0), (3, 1)]:
    J |= 1 << (4 * r + c)


def form(a):
    """theta^T J theta: pullback of J along a."""
    return mul(mul(tr(a), J), a)


def gl_list():
    return [a for a in range(1 << 16) if rank(a) == 4]


def inverse(a):
    for b in GL:
        if mul(a, b) == I4:
            return b
    raise ValueError


GL = gl_list()


def alt_nondeg_forms():
    """All nondegenerate alternating forms (symmetric, zero diagonal) on F_2^4 packed: 28 of them."""
    out = []
    for bits in range(64):
        pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        a = 0
        for k, (r, c) in enumerate(pairs):
            if (bits >> k) & 1:
                a |= 1 << (4 * r + c)
                a |= 1 << (4 * c + r)
        if rank(a) == 4:
            out.append(a)
    return out


NA = alt_nondeg_forms()
assert len(NA) == 28 and J in NA


def sp_list():
    return [a for a in GL if form(a) == J]


SP = sp_list()
assert len(SP) == 720


def congr(g, q):
    """g^T q g"""
    return mul(mul(tr(g), q), g)


def psi_for_form(q):
    """some psi in GL(4,2) with psi^T J psi = q."""
    for a in GL:
        if form(a) == q:
            return a
    raise ValueError
