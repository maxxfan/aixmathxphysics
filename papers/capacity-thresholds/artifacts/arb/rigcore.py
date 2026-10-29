"""rigcore: exact inputs and rigorous (Arb ball) channel application for p08.

Written for results/rigorous/arb_ic.py.  Shares no code with verifier/ (or with results/indep/,
whose mathematics it follows: weighted Schur blocks, Lindblad generator, spin-1/2 Clebsch-Gordan
superoperator, Poisson uniformisation).

Conventions
-----------
* Dicke index a = number of ones; |0> = spin up (s2 = +1), |1> = spin down (s2 = -1).
* Schur blocks k = 0..n//2, spin j = n/2 - k, j2 = 2j = n - 2k, dimension d = j2 + 1, row index
  a = j - m (a = 0 is m = j).  An S_n-invariant operator X = (+)_j X_j (x) 1_{m_j} is stored by
  its weighted blocks Xt_j = m_j X_j (so tr Xt_j is the probability of sector j and
  ||X||_1 = sum_j ||Xt_j||_1).
* A Pauli channel with PTM eigenvalues lx = ly (z-covariant: depolarizing in the plain frame,
  independent_xz and two_pauli after the per-qubit frame change U = exp(-i pi/4 X), which swaps
  p_Y and p_Z) is exp(G),  G = g_xy (S_X + S_Y - 2n) + g_z (S_Z - n),  S_P = sum_i P_i . P_i.
  S_X + S_Y = 2 sum_i (E01_i . E10_i + E10_i . E01_i)  (E01 = |0><1|).  Both conserve
  Delta = b - a, so every Delta-diagonal ("sector") of every block evolves independently.
* exp(G) = exp(G_+) exp(G_-) (the S_P commute).  G_+ (positive rates) is applied by Poisson
  uniformisation  exp(G_+) = sum_k Pois(k; mu) M^k  with M = (g_xy (S_X+S_Y) + g_z^+ S_Z)/mu a
  convex combination of Pauli conjugations (CPTP, so ||M^k Y||_1 <= ||Y||_1 and the tail is
  bounded by Pois(> K) ||Y||_1).  G_- = g_z (S_Z - n) with g_z < 0 (two_pauli only) is applied by
  its Taylor series, tail bounded with ||S_Z - n||_{1->1} <= 2n.

Rigour: every number below is an Arb ball (python-flint) or an exact integer/rational; the
superoperator coefficients are exact square roots of rationals evaluated in Arb; the Poisson
weights and rates use Arb's rigorous exp/log.  The ball iteration contains the exact truncated
series; the truncation tails are bounded in trace norm and returned separately (eps_tail) to be
added as a Weyl perturbation to every block spectrum.
"""
import base64
import json
from fractions import Fraction
from math import comb

import numpy as np
from flint import acb, arb, ctx, fmpq, fmpz

# ----------------------------------------------------------------------------- input


def _parse_double(x):
    if isinstance(x, bool):
        raise ValueError('boolean is not a number')
    if isinstance(x, str):
        s = x.strip()
        if s.lower().lstrip('+-').startswith('0x'):
            return float.fromhex(s)
        return float(s)
    if isinstance(x, (int, float)):
        return float(x)
    raise ValueError('bad number %r' % (x,))


def load_candidate(path):
    """returns dict(target, p (Fraction), pstr, rows: list of list of (re, im) python floats)."""
    with open(path) as f:
        d = json.load(f)
    target = d['target']
    pstr = d['p']
    if not isinstance(pstr, str):
        raise ValueError('p must be a string')
    if 'rows' in d:
        rows = []
        for row in d['rows']:
            rr = []
            for x in row:
                if isinstance(x, (list, tuple)):
                    if len(x) != 2:
                        raise ValueError('complex entry must be [re, im]')
                    rr.append((_parse_double(x[0]), _parse_double(x[1])))
                else:
                    rr.append((_parse_double(x), 0.0))
            rows.append(rr)
    else:
        raw = base64.b64decode(d['payload_c128_b64'])
        V = np.frombuffer(raw, dtype='<c16')
        r = int(d['rank'])
        V = V.reshape(r, -1)
        rows = [[(float(z.real), float(z.imag)) for z in row] for row in V]
    n1 = len(rows[0])
    if any(len(row) != n1 for row in rows):
        raise ValueError('rows of unequal length')
    return dict(target=target, pstr=pstr, p=Fraction(pstr), rows=rows, note=d.get('note', ''))


def pauli_probs(target, p):
    p = Fraction(p)
    if target == 'depolarizing':
        return (1 - 3 * p, p, p, p)
    if target == 'independent_xz':
        return ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p))
    if target == 'two_pauli':
        return (1 - 2 * p, p, Fraction(0), p)
    raise ValueError('unknown target ' + target)


def A(q):
    """exact Fraction -> arb ball."""
    q = Fraction(q)
    return arb(fmpq(q.numerator, q.denominator))


# ----------------------------------------------------------------------------- Schur layout


def mult(N, j2):
    """multiplicity of spin j = j2/2 in N qubits."""
    if j2 < 0 or j2 > N or (N - j2) % 2:
        return 0
    k = (N - j2) // 2
    return comb(N, k) - (comb(N, k - 1) if k >= 1 else 0)


# single-qubit operators as {(row s2, col t2): value}, s2 = +1 <-> |0>
OP_E01 = {(+1, -1): 1}
OP_E10 = {(-1, +1): 1}
OP_Z = {(+1, +1): 1, (-1, -1): -1}

_SQRT_CACHE = {}


def _sqrt_q(num, den):
    key = (num, den)
    v = _SQRT_CACHE.get(key)
    if v is None:
        v = arb(fmpq(num, den)).sqrt()
        _SQRT_CACHE[key] = v
    return v


def cg(j2pp, j2, m2, s2):
    """Condon-Shortley <j'', m - s; 1/2, s | j, m> (twice-values), exact zero outside range."""
    if abs(m2) > j2 or abs(m2 - s2) > j2pp or j2pp < 0:
        return None
    den = 2 * (j2pp + 1)
    if j2 == j2pp + 1:
        num = j2pp + s2 * m2 + 1
        if num <= 0:
            return None
        return _sqrt_q(num, den)
    if j2 == j2pp - 1:
        num = j2pp - s2 * m2 + 1
        if num <= 0:
            return None
        v = _sqrt_q(num, den)
        return -v if s2 > 0 else v
    raise ValueError('cg: |j - j\'\'| != 1/2')


ZERO = arb(0)


def rowfac(j2pp, j2p, j2in, op, d):
    """R[a] = sum_{(s,t)} op[s,t] <j'',m-s;1/2 s|j',m> <j'',m-s;1/2 t|j_in,m-s+t>, m2 = j2p - 2a."""
    out = np.empty(d, dtype=object)
    for a in range(d):
        m2 = j2p - 2 * a
        acc = None
        for (s2, t2), val in op.items():
            c1 = cg(j2pp, j2p, m2, s2)
            if c1 is None:
                continue
            c2 = cg(j2pp, j2in, m2 - s2 + t2, t2)
            if c2 is None:
                continue
            term = c1 * c2
            if val != 1:
                term = term * val
            acc = term if acc is None else acc + term
        out[a] = ZERO if acc is None else acc
    return out


def colfac(j2pp, j2p, j2in, op, d):
    """C[b] = sum_{(t',s')} op[t',s'] <j'',m'-s';1/2 t'|j_in,m'-s'+t'> <j'',m'-s';1/2 s'|j',m'>."""
    out = np.empty(d, dtype=object)
    for b in range(d):
        m2 = j2p - 2 * b
        acc = None
        for (t2, s2), val in op.items():
            c1 = cg(j2pp, j2p, m2, s2)
            if c1 is None:
                continue
            c2 = cg(j2pp, j2in, m2 - s2 + t2, t2)
            if c2 is None:
                continue
            term = c1 * c2
            if val != 1:
                term = term * val
            acc = term if acc is None else acc + term
        out[b] = ZERO if acc is None else acc
    return out


def _er(op):
    ers = {s2 - t2 for (s2, t2) in op}
    assert len(ers) == 1
    return ers.pop()


class Stencil:
    """Superoperator  Y -> sum_ops w_op sum_i A_i Y B_i  (+ shift * Y) on weighted Schur blocks,
    stored as per-(block, delta, op) rank-one factor lists, to be specialised per Delta-sector.

    ops: list of (A, B, weight arb).  key = (delta, er) with j2_in = j2_out + 2 delta,
    a_in = a + delta + er/2, b_in = b + delta + er/2."""

    def __init__(self, n, ops, shift=None):
        self.n = n
        self.nb = n // 2 + 1
        self.shift = shift
        self.terms = {}   # (kp, key) -> list of (Fw arb, R array, C array)
        for kp in range(self.nb):
            j2p = n - 2 * kp
            d = j2p + 1
            for delta in (-1, 0, 1):
                j2in = j2p + 2 * delta
                if j2in < 0 or j2in > n:
                    continue
                for j2pp in (j2p - 1, j2p + 1):
                    if j2pp < 0 or j2pp > n - 1 or abs(j2pp - j2in) != 1:
                        continue
                    mpp = mult(n - 1, j2pp)
                    if mpp == 0:
                        continue
                    F = arb(fmpq(n * mpp, mult(n, j2in)))
                    for opA, opB, w in ops:
                        er = _er(opA)
                        ec = -_er(opB)        # opB keys are (row t', col s'): ec = s' - t'
                        assert er == ec, 'operator pair is not Delta-conserving'
                        R = rowfac(j2pp, j2p, j2in, opA, d)
                        C = colfac(j2pp, j2p, j2in, opB, d)
                        key = (delta, er)
                        self.terms.setdefault((kp, key), []).append((F * w, R, C))
        self.keys = sorted({key for (_, key) in self.terms})

    def sector(self, D):
        """flat layout and (SRC, COEF) arrays for sector Delta = D >= 0."""
        n = self.n
        nb = (n - D) // 2 + 1                    # blocks with d > D
        dims = [n - 2 * k + 1 - D for k in range(nb)]
        off = [0]
        for x in dims:
            off.append(off[-1] + x)
        L = off[-1]
        out = []
        for key in self.keys:
            delta, er = key
            SRC = np.full(L, L, dtype=np.int64)  # L = dummy zero row
            COEF = np.full(L, ZERO, dtype=object)
            for kp in range(nb):
                lst = self.terms.get((kp, key))
                if not lst:
                    continue
                kin = kp - delta
                dkp = dims[kp]
                a = np.arange(dkp)
                ain = a + delta + er // 2
                if 0 <= kin < nb:
                    valid = (ain >= 0) & (ain < dims[kin])
                else:
                    valid = np.zeros(dkp, dtype=bool)
                c = None
                for Fw, R, C in lst:
                    t = Fw * (R[:dkp] * C[D:D + dkp])
                    c = t if c is None else c + t
                # masked entries must have an exactly zero coefficient (CG out of range)
                for idx in np.nonzero(~valid)[0]:
                    if not c[idx].is_zero():
                        raise AssertionError('nonzero coefficient on an out-of-range source')
                idxv = np.nonzero(valid)[0]
                SRC[off[kp] + idxv] = off[kin] + ain[idxv]
                COEF[off[kp] + idxv] = c[idxv]
            if self.shift is not None and key == (0, 0):
                # identity part of (S - shift): every entry maps to itself (all valid for (0,0))
                assert np.all(SRC == np.arange(L)), 'identity key must be complete'
                COEF = COEF + self.shift
            out.append((SRC, COEF))
        if self.shift is not None and (0, 0) not in self.keys:
            out.append((np.arange(L, dtype=np.int64), np.full(L, self.shift, dtype=object)))
        return dict(nb=nb, dims=dims, off=off, L=L, ops=out)


# ----------------------------------------------------------------------------- rates, frame


class Channel:
    def __init__(self, target, p, n):
        self.target = target
        self.p = Fraction(p)
        self.n = n
        probs = pauli_probs(target, self.p)
        if target == 'depolarizing':
            self.frame = 'plain'
            pf = probs
        else:
            self.frame = 'rot'
            pI, pX, pY, pZ = probs
            pf = (pI, pX, pZ, pY)          # U Y U^+ = Z, U Z U^+ = -Y for U = exp(-i pi/4 X)
        self.probs = probs
        self.probs_frame = pf
        pI, pX, pY, pZ = pf
        lx = pI + pX - pY - pZ
        ly = pI - pX + pY - pZ
        lz = pI - pX - pY + pZ
        if not (lx > 0 and ly > 0 and lz > 0):
            raise ValueError('non-positive PTM eigenvalue')
        if lx != ly:
            raise ValueError('channel is not z-covariant in this frame (l_x != l_y)')
        self.lam = (lx, ly, lz)
        Lx = -A(lx).log() / 2
        Lz = -A(lz).log() / 2
        self.g_xy = Lz / 2                         # g_X = g_Y = (L_y + L_z - L_x)/2
        if lz == lx * ly:
            self.g_z_sign = 0
            self.g_z = arb(0)
        else:
            self.g_z_sign = 1 if lz > lx * ly else -1   # g_z = L_x - L_z/2 > 0 iff lx^2 < lz
            self.g_z = Lx - Lz / 2
            if self.g_z_sign > 0:
                assert self.g_z > 0
            else:
                assert self.g_z < 0
        gpos = 2 * self.g_xy + (self.g_z if self.g_z_sign > 0 else 0)
        self.mu = n * gpos
        w_xy = 2 * self.g_xy / self.mu
        ops = [(OP_E01, OP_E10, w_xy), (OP_E10, OP_E01, w_xy)]
        if self.g_z_sign > 0:
            ops.append((OP_Z, OP_Z, self.g_z / self.mu))
        self.ops_pos = ops
        self.ops_neg = [(OP_Z, OP_Z, arb(1))] if self.g_z_sign < 0 else None


def poisson_weights(mu, tail_target):
    """w_0..w_K (arb) and a rigorous upper bound on P(Pois(mu) > K) <= tail_target."""
    w = [(-mu).exp()]
    k = 0
    while True:
        k += 1
        w.append(w[-1] * mu / k)
        K = k
        # tail after K: sum_{i > K} w_i <= w_{K+1} / (1 - mu/(K+2))
        if K + 2 > 2 * mu.upper():
            wn = w[-1] * mu / (K + 1)
            tail = wn / (1 - mu / (K + 2))
            if tail.upper() < tail_target:
                return w, tail.upper()
        if k > 100000:
            raise RuntimeError('Poisson series too long')


def taylor_terms(x, tail_target):
    """number of Taylor terms K with sum_{k>K} x^k/k! <= tail_target (x arb upper bound)."""
    t = arb(1)
    k = 0
    while True:
        k += 1
        t = t * x / k
        if k + 2 > 2 * x.upper():
            tn = t * x / (k + 1)
            tail = tn / (1 - x / (k + 2))
            if tail.upper() < tail_target:
                return k, tail.upper()
        if k > 100000:
            raise RuntimeError('Taylor series too long')


def krawtchouk_rotation(n):
    """D[b][a] = <D_b| U^{(x)n} |D_a>, U = exp(-i pi/4 X) = (1 - iX)/sqrt 2, as acb balls.

    <y|U^{(x)n}|x> = 2^{-n/2} (-i)^{d_H(x,y)}  =>
    D[b,a] = 2^{-n/2} sqrt(C(n,a)/C(n,b)) (-i)^{a+b} sum_c (-1)^c C(a,c) C(n-a,b-c)."""
    two = arb(2)
    s = two ** (-arb(n) / 2) if n % 2 else arb(fmpq(1, 2 ** (n // 2)))
    ipow = [(1, 0), (0, -1), (-1, 0), (0, 1)]       # (-i)^k
    D = [[None] * (n + 1) for _ in range(n + 1)]
    for a in range(n + 1):
        for b in range(n + 1):
            K = 0
            for c in range(max(0, b - (n - a)), min(a, b) + 1):
                t = comb(a, c) * comb(n - a, b - c)
                K += -t if c % 2 else t
            re, im = ipow[(a + b) % 4]
            mag = s * arb(fmpq(comb(n, a), comb(n, b))).sqrt() * K
            D[b][a] = acb(mag * re, mag * im)
    return D


def input_rows(rows, frame):
    """rows of python floats -> list of rows of acb (exact for plain frame)."""
    V = [[acb(arb(re), arb(im)) for (re, im) in row] for row in rows]
    if frame == 'plain':
        return V
    n = len(rows[0]) - 1
    D = krawtchouk_rotation(n)
    out = []
    for row in V:
        nr = []
        for b in range(n + 1):
            s = acb(0)
            for a in range(n + 1):
                if row[a] != 0:
                    s += D[b][a] * row[a]
            nr.append(s)
        out.append(nr)
    return out


def exact_Z(rows):
    """Z = sum_q ||v_q||^2, exact rational (the rows are exact doubles)."""
    Z = Fraction(0)
    for row in rows:
        for re, im in row:
            Z += Fraction(re) ** 2 + Fraction(im) ** 2
    return Z


# ----------------------------------------------------------------------------- sector iteration


def apply_stencil(sec, Y, Lact):
    """sum over keys of COEF[:Lact] * Y[SRC[:Lact]]  (Y has a trailing zero row)."""
    acc = None
    for SRC, COEF in sec['ops']:
        t = COEF[:Lact, None] * Y[SRC[:Lact]]
        acc = t if acc is None else acc + t
    return acc


def run_sector(D, ctxd):
    """Evolve sector Delta = D of all columns.  ctxd: dict with stencils, input, weights.
    returns (off, list of per-column object arrays of arb, stats)."""
    ctx.prec = ctxd['prec']
    st_pos = ctxd['stencil_pos']
    sec = st_pos.sector(D)
    L, off, nb = sec['L'], sec['off'], sec['nb']
    cols = ctxd['cols']                 # list of (q, q2, part)
    V = ctxd['V']                       # list of rows of acb (frame)
    Zinv = ctxd['Zinv']                 # arb
    ncol = len(cols)
    d0 = sec['dims'][0]                 # = n + 1 - D
    Y = np.full((L + 1, ncol), ZERO, dtype=object)
    for ci, (q, q2, part) in enumerate(cols):
        for a in range(d0):
            z = V[q][a] * V[q2][a + D].conjugate() * Zinv
            Y[a, ci] = z.real if part == 0 else z.imag
    active = 0                           # highest block index that may be nonzero
    stats = {}
    if ctxd['taylor'] is not None:
        st_neg, gz, KT = ctxd['taylor']
        secn = st_neg.sector(D)
        T = Y
        S = Y[:L].copy()
        for k in range(1, KT + 1):
            active = min(active + 1, nb - 1)
            La = off[active + 1]
            Tn = np.full((L + 1, ncol), ZERO, dtype=object)
            Tn[:La] = apply_stencil(secn, T, La) * (gz / k)
            T = Tn
            S[:La] = S[:La] + T[:La]
        Y = np.full((L + 1, ncol), ZERO, dtype=object)
        Y[:L] = S
        stats['taylor_terms'] = KT
    w = ctxd['weights']
    acc = Y[:L] * w[0]
    for k in range(1, len(w)):
        active = min(active + 1, nb - 1)
        La = off[active + 1]
        Yn = np.full((L + 1, ncol), ZERO, dtype=object)
        Yn[:La] = apply_stencil(sec, Y, La)
        Y = Yn
        acc[:La] = acc[:La] + Y[:La] * w[k]
    stats['poisson_terms'] = len(w)
    return off, acc, stats


def ser_arb(x):
    m, e = x.mid().man_exp()
    rm, re = x.rad().man_exp()
    return (int(m), int(e), int(rm), int(re))


def deser_arb(t):
    m, e, rm, re = t
    return arb((m, e), (rm, re))


def ser_array(xs):
    """object array of arb -> (mantissas: object array of python ints, exps, rad mantissas,
    rad exps: int64 arrays); exact."""
    L = len(xs)
    man = np.empty(L, dtype=object)
    ex = np.empty(L, dtype=np.int64)
    rman = np.empty(L, dtype=np.int64)
    rex = np.empty(L, dtype=np.int64)
    for i, x in enumerate(xs):
        m, e = x.mid().man_exp()
        rm, re = x.rad().man_exp()
        man[i] = int(m)
        ex[i] = int(e)
        rman[i] = int(rm)
        rex[i] = int(re)
    return (man, ex, rman, rex)


def deser_at(arrs, i):
    man, ex, rman, rex = arrs
    return arb((man[i], int(ex[i])), (int(rman[i]), int(rex[i])))
