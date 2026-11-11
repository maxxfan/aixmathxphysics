#!/usr/bin/env python3
"""Independent evaluator of the coherent information I_c(N^{(x)n}, rho) for p08 candidates.

Usage:
    indep_ic.py FILE [P] [--no-check] [--fast] [--method uniform|expm|dense] [--frame auto|rot|plain]

    FILE  candidate JSON ("rows" or "payload_c128_b64"+"rank", Dicke basis, see PROBLEM.md)
    P     optional exact p ("a/b" or decimal) overriding the file's p
    Default: uniformisation + an expm_multiply re-run + eigen-stage probes (n=240 rank 2: ~2 min,
    ~3 GB; rank 4: ~6-10 min, ~5.5 GB).  --fast: one method, no probes (n=240 rank 2: ~30 s).

    Error lines: 'empirical' = 10 x the largest of |evd - evr| (two LAPACK eigensolvers),
    |random-rotation probe| (blocks conjugated by a random unitary before diagonalising) and
    |uniform - expm|; 'first-order bound' and 'worst-case floor' charge every eigenvalue a
    LAPACK-sized error (2 sqrt(d) u ||A||, resp. 4 d u ||A|| for the unresolved ones).  The sign
    line uses max(bound, floor) + empirical.  Not a rigorous enclosure (that is verifier/'s job).

Auditor code, written without reading verifier/, qlib/ or work/lead/.

Method (a different route from the coherent-state block formula det^k Sym^{2j}):
  * The input operators |v_q><v_q'| commute with S_n, as does the channel, so every operator in
    play is  (+)_j  X_j (x) 1_{m_j}  on the Schur-Weyl decomposition.  We store the weighted blocks
    Xt_j = m_j X_j (Dicke-type basis |j,m>, m = j..-j), whose traces are sector probabilities.
  * A Pauli channel with all PTM eigenvalues positive is exp(G) with the Lindblad generator
    G = sum_i sum_P gamma_P (P_i . P_i - id), gamma_P from the eigenvalues (two_pauli has one
    negative gamma; that factor is a commuting 'anti-dephasing' handled separately).
  * S_P = sum_i P_i . P_i acting on S_n-invariant operators is built in closed form by splitting
    off the last qubit: spin-1/2 Clebsch-Gordan coefficients, j -> j'' = j +- 1/2 -> j' and the
    branching multiplicities m^{(n-1)}_{j''} (a PIQS-style Dicke-space Liouvillian).
  * exp(G) is applied by uniformisation: exp(n Gamma (M - 1)) = sum_k Pois(k; n Gamma) M^k with
    M = (1/(n Gamma)) sum_{gamma_P>0} gamma_P S_P a trace-preserving CP map (all terms are images
    of the input under a channel, so no cancellation between terms).  The negative-rate factor
    (two_pauli) is a Taylor series applied first.  Cross-check: scipy expm_multiply on G.
  * independent_xz / two_pauli are evaluated in the rotated frame U = exp(-i pi/4 X) per qubit,
    where p_Y <-> p_Z and the channel becomes z-covariant (generator conserves m-m'); the
    rotation on Sym^n is computed exactly (Gaussian-integer coefficients, mpmath rounding).
    --frame plain evaluates in the original frame instead (Delta-mixing generator).
  * I_c = sum_j [H(Xt_j^B) - H(Xt_j^{BR})] in bits (the log m_j terms cancel because the traces
    agree), eigenvalues by LAPACK on each (2j+1) and (2j+1)r block.
"""
import argparse
import base64
import json
import math
import sys
import time
from fractions import Fraction
from math import comb

import mpmath
import numpy as np
import scipy.sparse as sp

U_ROUND = 2.0 ** -53
LN2 = math.log(2.0)

# ----------------------------------------------------------------------------- input


def _num(x):
    if isinstance(x, str):
        s = x.strip()
        if s.lower().startswith(('0x', '-0x', '+0x')):
            return float.fromhex(s)
        return float(s)
    return float(x)


def load_candidate(path):
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
                    rr.append(complex(_num(x[0]), _num(x[1])))
                else:
                    rr.append(complex(_num(x), 0.0))
            rows.append(rr)
        V = np.array(rows, dtype=np.complex128)
    else:
        raw = base64.b64decode(d['payload_c128_b64'])
        V = np.frombuffer(raw, dtype='<c16').astype(np.complex128)
        V = V.reshape(int(d['rank']), -1)
    return target, pstr, V, d.get('note', '')


def pauli_probs(target, p):
    p = Fraction(p)
    if target == 'depolarizing':
        return (1 - 3 * p, p, p, p)
    if target == 'independent_xz':
        return ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p))
    if target == 'two_pauli':
        return (1 - 2 * p, p, Fraction(0), p)
    raise ValueError('unknown target ' + target)


def lindblad_rates(probs):
    """gamma_X, gamma_Y, gamma_Z (floats) with exp(sum_P gamma_P (P.P - 1)) = the channel."""
    pI, pX, pY, pZ = probs
    lx = pI + pX - pY - pZ
    ly = pI - pX + pY - pZ
    lz = pI - pX - pY + pZ
    if min(lx, ly, lz) <= 0:
        raise ValueError('channel has a non-positive PTM eigenvalue; not exp of a generator')
    with mpmath.workdps(50):
        Lx = -mpmath.log(mpmath.mpf(lx.numerator) / lx.denominator) / 2
        Ly = -mpmath.log(mpmath.mpf(ly.numerator) / ly.denominator) / 2
        Lz = -mpmath.log(mpmath.mpf(lz.numerator) / lz.denominator) / 2
        gX = (Ly + Lz - Lx) / 2
        gY = (Lx + Lz - Ly) / 2
        gZ = (Lx + Ly - Lz) / 2
        big = max(abs(gX), abs(gY), abs(gZ))
        # exact zeros (e.g. gamma_Y of independent_xz) come out as ~1e-50 residues: snap them
        return tuple(0.0 if abs(g) < 1e-30 * big else float(g) for g in (gX, gY, gZ))


# ----------------------------------------------------------------------------- Schur layout


def mult(N, j2):
    """multiplicity of spin j = j2/2 in N qubits."""
    k = (N - j2) // 2
    if k < 0 or (N - j2) % 2:
        return 0
    return comb(N, k) - (comb(N, k - 1) if k >= 1 else 0)


class Layout:
    def __init__(self, n):
        self.n = n
        self.j2 = [n - 2 * k for k in range(n // 2 + 1)]
        self.d = [j2 + 1 for j2 in self.j2]
        self.off = [0]
        for d in self.d:
            self.off.append(self.off[-1] + d * d)
        self.N = self.off[-1]
        self.m = [mult(n, j2) for j2 in self.j2]

    def delta_of(self):
        """Delta = b - a ( = m - m') for every variable."""
        out = np.empty(self.N, dtype=np.int64)
        for k, d in enumerate(self.d):
            a = np.arange(d)
            out[self.off[k]:self.off[k + 1]] = (a[None, :] - a[:, None]).ravel()
        return out


def _cg(j2pp, j2, m2, s2):
    """<j'', m-s; 1/2, s | j, m> (Condon-Shortley), all angular momenta as twice-values."""
    den = 2.0 * (j2pp + 1)
    if j2 == j2pp + 1:
        num = j2pp + s2 * m2 + 1
        return np.sqrt(np.maximum(num, 0) / den)
    if j2 == j2pp - 1:
        num = j2pp - s2 * m2 + 1
        return -s2 * np.sqrt(np.maximum(num, 0) / den)
    raise ValueError


_IDX = {1: 0, -1: 1}   # s2=+1 <-> |0> (spin up), s2=-1 <-> |1>


def build_superop(L, A, B):
    """Sparse matrix of  rho -> sum_i A_i rho B_i  on weighted Schur blocks Xt_j = m_j X_j.

    Derivation: Y_{j'}[m,m'] = (n/m_{j'}) sum_alpha <j' m alpha| A_n rho B_n |j' m' alpha>; split
    alpha = (j'', alpha'') over the first n-1 qubits, |j' m (j'' a'')> = sum_s CG |j'', m-s, a''>|s>,
    and the same for rho's basis; the alpha'' sum gives m^{(n-1)}_{j''}.
    """
    n = L.n
    rows, cols, vals = [], [], []
    for kp, j2p in enumerate(L.j2):
        dp = j2p + 1
        a = np.arange(dp)
        m2 = j2p - 2 * a
        acc = {}
        for j2pp in (j2p - 1, j2p + 1):
            if j2pp < 0 or j2pp > n - 1:
                continue
            mpp = mult(n - 1, j2pp)
            for j2in in (j2pp - 1, j2pp + 1):
                if j2in < 0 or j2in > n:
                    continue
                kin = (n - j2in) // 2
                F = float(Fraction(n * mpp, mult(n, j2in)))
                for s2 in (1, -1):
                    for t2 in (1, -1):
                        Ast = A[_IDX[s2], _IDX[t2]]
                        if Ast == 0:
                            continue
                        rowfac = _cg(j2pp, j2p, m2, s2) * _cg(j2pp, j2in, m2 - s2 + t2, t2) * Ast
                        for s2p in (1, -1):
                            for t2p in (1, -1):
                                Bts = B[_IDX[t2p], _IDX[s2p]]
                                if Bts == 0:
                                    continue
                                colfac = Bts * _cg(j2pp, j2in, m2 - s2p + t2p, t2p) * _cg(j2pp, j2p, m2, s2p)
                                key = (kin, s2 - t2, s2p - t2p)
                                term = F * np.outer(rowfac, colfac)
                                if key in acc:
                                    acc[key] = acc[key] + term
                                else:
                                    acc[key] = term
        for (kin, er, ec), C in acc.items():
            j2in = L.j2[kin]
            din = j2in + 1
            a_in = a + (j2in - j2p + er) // 2
            b_in = a + (j2in - j2p + ec) // 2
            va = (a_in >= 0) & (a_in <= j2in)
            vb = (b_in >= 0) & (b_in <= j2in)
            Cv = C[np.ix_(va, vb)]
            ra = a[va]
            rb = a[vb]
            R = L.off[kp] + ra[:, None] * dp + rb[None, :]
            Cc = L.off[kin] + a_in[va][:, None] * din + b_in[vb][None, :]
            nz = Cv != 0
            rows.append(R[nz].astype(np.int64))
            cols.append(Cc[nz].astype(np.int64))
            vals.append(Cv[nz])
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    if np.iscomplexobj(vals):
        if np.max(np.abs(vals.imag), initial=0.0) > 0:
            raise ValueError('complex superoperator coefficients')
        vals = vals.real
    M = sp.coo_matrix((vals, (rows, cols)), shape=(L.N, L.N)).tocsr()
    M.sum_duplicates()
    M.eliminate_zeros()
    return M


PAULI = {
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}
E01 = np.array([[0, 1], [0, 0]], dtype=complex)
E10 = np.array([[0, 0], [1, 0]], dtype=complex)


def pauli_superops(L, rates):
    """list of (gamma, S_P) with S_P = sum_i P_i . P_i (X+Y merged when gamma_X == gamma_Y)."""
    gX, gY, gZ = rates
    out = []
    if gX == gY:
        if gX != 0:
            Se = build_superop(L, E01, E10) + build_superop(L, E10, E01)
            out.append((gX, 2.0 * Se, 'X+Y'))
    else:
        if gX != 0:
            out.append((gX, build_superop(L, PAULI['X'], PAULI['X']), 'X'))
        if gY != 0:
            out.append((gY, build_superop(L, PAULI['Y'], PAULI['Y']), 'Y'))
    if gZ != 0:
        out.append((gZ, build_superop(L, PAULI['Z'], PAULI['Z']), 'Z'))
    # count: 'X+Y' carries two Paulis with the same rate
    return out


def _npaulis(name):
    return 2 if name == 'X+Y' else 1


# ----------------------------------------------------------------------------- frame rotation


def rotation_sym(n):
    """D[b,a] = <D_b| U^{(x)n} |D_a>, U = exp(-i pi/4 X), as mpmath matrix (exact integers, 60 digits)."""
    # [s^a t^b] (1 - i s - i t + s t)^n = sum_{n3} n!/(n0! n1! n2! n3!) (-i)^{n1+n2}
    fact = [1] * (n + 1)
    for i in range(1, n + 1):
        fact[i] = fact[i - 1] * i
    ipow = [(1, 0), (0, -1), (-1, 0), (0, 1)]   # (-i)^k
    with mpmath.workdps(60):
        c = mpmath.sqrt(mpmath.mpf(1) / 2) ** n
        D = [[None] * (n + 1) for _ in range(n + 1)]
        for a in range(n + 1):
            for b in range(n + 1):
                re = 0
                im = 0
                for n3 in range(0, min(a, b) + 1):
                    n1 = a - n3
                    n2 = b - n3
                    n0 = n - n1 - n2 - n3
                    if n0 < 0:
                        continue
                    mn = fact[n] // (fact[n0] * fact[n1] * fact[n2] * fact[n3])
                    pr, pi = ipow[(n1 + n2) % 4]
                    re += mn * pr
                    im += mn * pi
                s = c / mpmath.sqrt(mpmath.mpf(comb(n, a) * comb(n, b)))
                D[b][a] = mpmath.mpc(re * s, im * s)
    return D


def rotate_rows(V):
    """apply U^{(x)n} to each row of V (Dicke coefficients), in mpmath, rounded to complex128."""
    r, n1 = V.shape
    n = n1 - 1
    D = rotation_sym(n)
    out = np.zeros_like(V)
    with mpmath.workdps(60):
        for q in range(r):
            vq = [mpmath.mpc(float(V[q, a].real), float(V[q, a].imag)) for a in range(n1)]
            for b in range(n1):
                s = mpmath.mpc(0)
                for a in range(n1):
                    s += D[b][a] * vq[a]
                out[q, b] = complex(s)
    return out


# ----------------------------------------------------------------------------- channel application


def input_matrix(L, V):
    """N x r^2 complex: column q*r+q' = |v_q><v_q'| / Z in block j=n/2."""
    r, n1 = V.shape
    Z = float(np.sum(np.abs(V) ** 2))
    Vn = V / math.sqrt(Z)
    X = np.zeros((L.N, r * r), dtype=np.complex128)
    d = L.d[0]
    for q in range(r):
        for qq in range(r):
            X[0:d * d, q * r + qq] = np.outer(Vn[q], Vn[qq].conj()).ravel()
    return X, Z


def _poisson_weights(mu, tail=1e-20):
    """Poisson(mu) weights w_0..w_K with tail mass < tail (mpmath, then float)."""
    with mpmath.workdps(40):
        mu = mpmath.mpf(mu)
        w = []
        k = 0
        term = mpmath.exp(-mu)
        cum = mpmath.mpf(0)
        while True:
            w.append(term)
            cum += term
            if k > mu and 1 - cum < tail:
                break
            k += 1
            term = term * mu / k
        return [float(x) for x in w], float(1 - cum)


def _realify(X):
    return np.concatenate([X.real, X.imag], axis=1)


def _complexify(Y):
    h = Y.shape[1] // 2
    return Y[:, :h] + 1j * Y[:, h:]


def apply_uniform(L, terms, X0, verbose=False):
    """exp(G) X0 with G = sum gamma_P (S_P - n_P n).  Positive rates: uniformisation;
    negative rates: Taylor series applied first (it commutes with the rest)."""
    n = L.n
    Y = _realify(X0)
    pos = [(g, S, nm) for g, S, nm in terms if g > 0]
    neg = [(g, S, nm) for g, S, nm in terms if g < 0]
    info = {}
    if neg:
        # Gm = sum gamma (S - n*np) with gamma<0: spectrum in [0, 2 |gamma| n np]
        Gm = None
        shift = 0.0
        for g, S, nm in neg:
            Gm = g * S if Gm is None else Gm + g * S
            shift += g * n * _npaulis(nm)
        Gm = (Gm - shift * sp.identity(L.N, format='csr')).tocsr()
        bound = sum(-2 * g * n * _npaulis(nm) for g, S, nm in neg)
        acc = Y.copy()
        T = Y.copy()
        k = 0
        while True:
            k += 1
            T = (Gm @ T) / k
            acc += T
            if k > bound and np.max(np.abs(T)) < 1e-19 * max(1.0, np.max(np.abs(acc))):
                break
            if k > 2000:
                raise RuntimeError('negative-rate Taylor series did not converge')
        Y = acc
        info['neg_taylor_terms'] = k
        info['neg_growth_bound'] = math.exp(bound)
    Gam = sum(g * _npaulis(nm) for g, S, nm in pos)
    M = None
    for g, S, nm in pos:
        M = (g / (n * Gam)) * S if M is None else M + (g / (n * Gam)) * S
    M = M.tocsr()
    mu = n * Gam
    w, tail = _poisson_weights(mu)
    acc = w[0] * Y
    T = Y
    for k in range(1, len(w)):
        T = M @ T
        acc += w[k] * T
    info['poisson_mu'] = mu
    info['poisson_terms'] = len(w)
    info['poisson_tail'] = tail
    return _complexify(acc), info


def generator(L, terms):
    n = L.n
    G = None
    shift = 0.0
    for g, S, nm in terms:
        G = g * S if G is None else G + g * S
        shift += g * n * _npaulis(nm)
    G = (G - shift * sp.identity(L.N, format='csr')).tocsr()
    return G


def apply_expm_multiply(L, terms, X0):
    from scipy.sparse.linalg import expm_multiply
    G = generator(L, terms)
    Y = expm_multiply(G, _realify(X0))
    return _complexify(Y), {}


def apply_dense(L, terms, X0):
    """z-covariant generators only: dense expm on each Delta sector that the input touches."""
    from scipy.linalg import expm
    G = generator(L, terms)
    delta = L.delta_of()
    out = np.zeros_like(X0)
    used = np.nonzero(np.any(X0 != 0, axis=1))[0]
    for D in np.unique(delta[used]):
        idx = np.nonzero(delta == D)[0]
        Gd = G[idx][:, idx]
        # z-covariance check: no coupling out of the sector
        leak = abs(G[idx]).sum() - abs(Gd).sum()
        if leak > 0:
            raise ValueError('generator is not Delta-conserving; use --method uniform/expm')
        E = expm(Gd.toarray())
        out[idx] = E @ X0[idx]
    return out, {}


# ----------------------------------------------------------------------------- entropies


def _H(ev):
    pos = ev[ev > 0]
    return float(-(pos * np.log2(pos)).sum())


def _h(x):
    return 0.0 if x <= 0 else -x * math.log2(x)


def _first_order(ev, dim):
    """Entropy uncertainty if every eigenvalue may be off by eps = 2 sqrt(dim) u ||A|| (a typical
    backward-error size for Householder/QR on this block plus the block's own float64 error).
    Resolved eigenvalues: eps |log2 lam + 1/ln2|; unresolved ones (|lam| <= eps): h(2 eps)."""
    lam_max = float(np.max(np.abs(ev), initial=0.0))
    if lam_max == 0:
        return 0.0
    eps = 2.0 * math.sqrt(dim) * U_ROUND * lam_max
    big = ev[ev > eps]
    s = float(np.sum(eps * np.abs(np.log2(big) + 1.0 / LN2)))
    s += int(np.sum(ev <= eps)) * _h(2 * eps)
    return s


def block_entropies(L, X, r, want_detail=False, probe_seed=12345, probes=True):
    """I_c = sum_k H(sigma_k) - H(tau_k); returns I_c and diagnostics.

    Diagnostics: eig_alt  = same with LAPACK MRRR (evr) instead of divide-and-conquer (evd);
                 probe    = same after a random unitary similarity of every block (float64), i.e.
                            the entropies' sensitivity to perturbations of rounding size;
                 bound1   = first-order sum of _first_order over all blocks (semi-conservative);
                 floor    = h(4 d u ||A||) for every eigenvalue below 4 d u ||A|| (worst-case
                            LAPACK constant; counts only the unresolved eigenvalues)."""
    import scipy.linalg as sla
    rng = np.random.default_rng(probe_seed)
    total = 0.0
    total_alt = 0.0
    total_probe = 0.0
    floor = 0.0
    bound1 = 0.0
    trace_B = 0.0
    herm = 0.0
    negmin = 0.0
    detail = []

    def rot(A):
        m = A.shape[0]
        Q, R = np.linalg.qr(rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))
        B = Q @ A @ Q.conj().T
        return (B + B.conj().T) / 2

    for k, d in enumerate(L.d):
        blk = X[L.off[k]:L.off[k + 1]].reshape(d, d, r, r)
        tau = blk.transpose(0, 2, 1, 3).reshape(d * r, d * r)
        sig = np.zeros((d, d), dtype=complex)
        for q in range(r):
            sig += blk[:, :, q, q]
        herm = max(herm, float(np.max(np.abs(tau - tau.conj().T), initial=0.0)))
        tau = (tau + tau.conj().T) / 2
        sig = (sig + sig.conj().T) / 2
        es = np.linalg.eigvalsh(sig)
        et = np.linalg.eigvalsh(tau)
        hs, ht = _H(es), _H(et)
        total += hs - ht
        # second LAPACK driver (MRRR) as an eigen-stage noise probe
        if probes:
            es2 = sla.eigvalsh(sig, driver='evr')
            et2 = sla.eigvalsh(tau, driver='evr')
            total_alt += _H(es2) - _H(et2)
            total_probe += _H(np.linalg.eigvalsh(rot(sig))) - _H(np.linalg.eigvalsh(rot(tau)))
        else:
            total_alt += hs - ht
            total_probe += hs - ht
        trace_B += float(np.trace(sig).real)
        for ev, dim in ((es, d), (et, d * r)):
            lam_max = max(float(np.max(np.abs(ev), initial=0.0)), 0.0)
            if lam_max == 0:
                continue
            bound1 += _first_order(ev, dim)
            eps = 4.0 * dim * U_ROUND * lam_max
            small = np.abs(ev) < eps
            negmin = min(negmin, float(np.min(ev)) / lam_max)
            if np.any(small):
                floor += int(small.sum()) * _h(eps)
        if want_detail:
            detail.append((k, L.j2[k], float(np.trace(sig).real), hs, ht))
    diag = dict(trace=trace_B, herm_defect=herm, eig_alt=total_alt, probe=total_probe,
                bound1=bound1, floor=floor, negmin_rel=negmin)
    if want_detail:
        diag['detail'] = detail
    return total, diag


# ----------------------------------------------------------------------------- driver


def evaluate(target, p, V, method='uniform', frame='auto', verbose=False, probes=True):
    """returns dict with I_c (bits) and diagnostics."""
    t0 = time.time()
    probs = pauli_probs(target, p)
    r, n1 = V.shape
    n = n1 - 1
    rot = False
    if frame == 'auto':
        rot = target in ('independent_xz', 'two_pauli')
    elif frame == 'rot':
        rot = True
    if rot:
        pI, pX, pY, pZ = probs
        probs_f = (pI, pX, pZ, pY)        # U Y U^dag = Z, U Z U^dag = -Y for U = exp(-i pi/4 X)
        Vf = rotate_rows(V)
    else:
        probs_f = probs
        Vf = V
    rates = lindblad_rates(probs_f)
    L = Layout(n)
    t1 = time.time()
    terms = pauli_superops(L, rates)
    t2 = time.time()
    X0, Z = input_matrix(L, Vf)
    if method == 'uniform':
        X, info = apply_uniform(L, terms, X0)
    elif method == 'expm':
        X, info = apply_expm_multiply(L, terms, X0)
    elif method == 'dense':
        X, info = apply_dense(L, terms, X0)
    else:
        raise ValueError(method)
    t3 = time.time()
    ic, diag = block_entropies(L, X, r, probes=probes)
    t4 = time.time()
    return dict(n=n, rank=r, target=target, p=str(p), frame='rot' if rot else 'plain',
                method=method, rates=rates, Ic=ic, diag=diag, info=info, X=X, L=L,
                times=dict(setup=t1 - t0, superop=t2 - t1, channel=t3 - t2, entropy=t4 - t3,
                           total=t4 - t0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('file')
    ap.add_argument('p', nargs='?', default=None)
    ap.add_argument('--method', default='uniform', choices=['uniform', 'expm', 'dense'])
    ap.add_argument('--frame', default='auto', choices=['auto', 'rot', 'plain'])
    ap.add_argument('--fast', action='store_true',
                    help='skip the second method and the eigen-stage probes (about 4x faster)')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--check', dest='check', action='store_true', default=None,
                   help='also run a second exponential method, scipy expm_multiply (default: on)')
    g.add_argument('--no-check', dest='check', action='store_false')
    args = ap.parse_args(argv)

    target, pstr, V, note = load_candidate(args.file)
    if args.p is not None:
        pstr = args.p
    p = Fraction(pstr)
    r, n1 = V.shape
    n = n1 - 1
    res = evaluate(target, p, V, method=args.method, frame=args.frame, probes=not args.fast)
    check = args.check if args.check is not None else (not args.fast)
    d = res['diag']
    print(f'file     {args.file}')
    print(f'family   {target}   n {n}   rank {r}   p {pstr} (= {float(p):.12g})   frame {res["frame"]}')
    print(f'rates    gamma_X,Y,Z = {res["rates"][0]:.15g} {res["rates"][1]:.15g} {res["rates"][2]:.15g}')
    print(f'method   {res["method"]}  {res["info"]}')
    print(f'I_c      {res["Ic"]:.12e} bits')
    eig_noise = abs(d['eig_alt'] - res['Ic'])
    probe = abs(d['probe'] - res['Ic'])
    print(f'diag     trace-1 {d["trace"] - 1:.1e}  herm-defect {d["herm_defect"]:.1e}  '
          f'min-eig/max {d["negmin_rel"]:.1e}')
    if not args.fast:
        print(f'         |evd - evr| {eig_noise:.1e}   |random-rotation probe| {probe:.1e}')
    emp = [eig_noise, probe] if not args.fast else []
    if check:
        alt = 'expm' if res['method'] != 'expm' else 'uniform'
        res2 = evaluate(target, p, V, method=alt, frame=args.frame, probes=False)
        diffm = abs(res2['Ic'] - res['Ic'])
        print(f'check    {alt}: I_c {res2["Ic"]:.12e}   |diff| {diffm:.1e}   ({res2["times"]["total"]:.1f} s)')
        emp.append(diffm)
    err_emp = 10.0 * max(emp) if emp else 0.0
    err_b1 = d['bound1']
    if emp:
        print(f'error    empirical {err_emp:.1e} bits (10 x max of the differences above)')
    else:
        print('error    empirical: not measured (--fast)')
    print(f'         first-order bound {err_b1:.1e} bits (every eigenvalue off by 2 sqrt(dim) u ||block||)')
    print(f'         worst-case floor  {d["floor"]:.1e} bits (unresolved eigenvalues, 4 dim u ||block||)')
    err = max(err_b1, d['floor']) + err_emp
    sign = 'POSITIVE' if res['Ic'] - err > 0 else ('NEGATIVE' if res['Ic'] + err < 0 else 'UNRESOLVED')
    print(f'sign     {sign}   (I_c -+ (max(bound, floor) + empirical) = '
          f'[{res["Ic"] - err:.4e}, {res["Ic"] + err:.4e}])')
    t = res['times']
    print(f'time     total {t["total"]:.1f} s (superop {t["superop"]:.1f}, channel {t["channel"]:.1f}, '
          f'entropy {t["entropy"]:.1f})')
    return res


if __name__ == '__main__':
    main()
