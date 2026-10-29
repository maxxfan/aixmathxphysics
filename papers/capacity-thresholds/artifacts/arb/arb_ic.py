#!/usr/bin/env python3
"""Rigorous (Arb ball arithmetic) enclosure of the coherent information I_c for p08 candidates.

Usage:
    arb_ic.py CANDIDATE.json [--p P] [--prec BITS] [--nproc N] [--json OUT] [--verbose]

Prints a rigorous enclosure [lo, hi] of I_c (bits) and the sign (POSITIVE if lo > 0, NEGATIVE if
hi < 0, else UNDECIDED).  Independent of verifier/ (no shared code, different block formula):
see rigcore.py (channel) and rigeig.py (spectra) for the method and the rigour argument.

  Inputs:   p exact rational; amplitudes are the exact doubles; Z = sum ||v_q||^2 exact.
  Channel:  exp(G) on weighted Schur blocks, Delta-sector by Delta-sector, in Arb balls
            (CG coefficients = exact sqrt of rationals; Poisson weights from Arb exp/log);
            truncation tails bounded in trace norm -> Weyl slack on every block spectrum.
  Spectra:  float hint -> Arb refinement -> Ostrowski + Gershgorin certification.
  Entropy:  -x log2 x on every eigenvalue interval with Arb logs, concavity, max 1/(e ln 2).
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction

for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np                       # noqa: E402
from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq, fmpz   # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rigcore as RC                     # noqa: E402
import rigeig as RE                      # noqa: E402

CTX = {}       # inherited by forked workers
SECT = {}


def _sector_job(D):
    t0 = time.time()
    off, acc, stats = RC.run_sector(D, CTX)
    ncol = acc.shape[1]
    ser = [RC.ser_array(acc[:, ci]) for ci in range(ncol)]
    stats['time'] = time.time() - t0
    return D, off, ser, stats


def _entry(colmap, q, q2, k, a, b, cplx):
    """X^{(q,q2)}_k [a, b] as arb/acb ball."""
    if b >= a:
        D = b - a
        off, ser = SECT[D]
        pos = off[k] + a
        re = RC.deser_at(ser[colmap[(q, q2, 0)]], pos)
        if not cplx:
            return re
        im = RC.deser_at(ser[colmap[(q, q2, 1)]], pos)
        return acb(re, im)
    D = a - b
    off, ser = SECT[D]
    pos = off[k] + b
    re = RC.deser_at(ser[colmap[(q2, q, 0)]], pos)
    if not cplx:
        return re
    im = RC.deser_at(ser[colmap[(q2, q, 1)]], pos)
    return acb(re, -im)


def _block_job(k):
    t0 = time.time()
    ctx.prec = CTX['prec']
    n, r, cplx = CTX['n'], CTX['r'], CTX['cplx']
    colmap = CTX['colmap']
    slack = CTX['slack']
    tiny = CTX['tiny']
    d = n - 2 * k + 1
    # trace of sigma_k from the Delta = 0 sector
    off0, ser0 = SECT[0]
    tr = arb(0)
    for q in range(r):
        ci = colmap[(q, q, 0)]
        for a in range(d):
            tr += RC.deser_at(ser0[ci], off0[k] + a)
    out = dict(k=k, d=d, trace=RC.ser_arb(tr))
    thi = (tr + slack * (r * d)).upper()
    if thi < tiny:
        inv_e, ln2 = RE._consts()
        for name, dim in (('sigma', d), ('tau', r * d)):
            hi = (thi * (arb(dim) / thi).log() / ln2).upper() if thi > 0 else arb(0)
            out[name] = (RC.ser_arb(arb(0)), RC.ser_arb(hi), dict(d=dim, mode='tiny-trace'))
        out['time'] = time.time() - t0
        return out
    MT = acb_mat if cplx else arb_mat
    X = {}
    for q in range(r):
        for q2 in range(r):
            X[(q, q2)] = [[_entry(colmap, q, q2, k, a, b, cplx) for b in range(d)] for a in range(d)]
    sig = MT(d, d)
    for a in range(d):
        for b in range(d):
            s = X[(0, 0)][a][b]
            for q in range(1, r):
                s = s + X[(q, q)][a][b]
            sig[a, b] = s
    tau = MT(r * d, r * d)
    for q in range(r):
        for q2 in range(r):
            blk = X[(q, q2)]
            for a in range(d):
                row = blk[a]
                for b in range(d):
                    tau[q * d + a, q2 * d + b] = row[b]
    del X
    for name, H in (('sigma', sig), ('tau', tau)):
        lo, hi, info = RE.entropy_enclosure(H, cplx, slack * r, tiny=tiny)
        out[name] = (RC.ser_arb(lo), RC.ser_arb(hi), info)
    out['time'] = time.time() - t0
    return out


def fmt_dec(x, digits=30, up=False):
    """exact arb x -> decimal string rounded outward at 10^-digits."""
    s = x * fmpz(10) ** digits
    v = s.upper().ceil() if up else s.lower().floor()
    m = v.unique_fmpz()
    if m is None:
        raise ValueError('rounding not exact')
    return '%d' % int(m)


def evaluate(path, p_override=None, prec=128, nproc=6, tail=1e-45, tiny=1e-30, verbose=False,
             perturb=None):
    t_start = time.time()
    ctx.prec = prec
    cand = RC.load_candidate(path)
    target = cand['target']
    pstr = p_override if p_override is not None else cand['pstr']
    p = Fraction(pstr)
    if not (0 < p <= Fraction(1, 4)):
        raise ValueError('p out of range')
    rows = cand['rows']
    if perturb is not None:
        q, a, delta = perturb
        re, im = rows[q][a]
        rows[q][a] = (re + delta, im)
    r = len(rows)
    n = len(rows[0]) - 1
    ch = RC.Channel(target, p, n)
    Zq = RC.exact_Z(rows)
    Zinv = RC.A(1 / Zq)
    V = RC.input_rows(rows, ch.frame)
    cplx = ch.frame == 'rot' or any(im != 0.0 for row in rows for (_, im) in row)
    cols = []
    for q in range(r):
        for q2 in range(r):
            for part in ((0, 1) if cplx else (0,)):
                cols.append((q, q2, part))
    colmap = {c: i for i, c in enumerate(cols)}
    t0 = time.time()
    st_pos = RC.Stencil(n, ch.ops_pos)
    t_stencil = time.time() - t0
    weights, tailP = RC.poisson_weights(ch.mu, tail)
    taylor = None
    growth = arb(1)
    tailT = arb(0)
    KT = 0
    if ch.ops_neg is not None:
        st_neg = RC.Stencil(n, ch.ops_neg, shift=arb(-n))
        x = (2 * n * abs(ch.g_z)).upper()
        KT, tailT = RC.taylor_terms(x, tail)
        growth = x.exp().upper()
        taylor = (st_neg, ch.g_z, KT)
    # trace-norm bound of (exact output - exact truncated series), per column (||X0||_1 <= 1)
    eps = (tailT + tailP * growth).upper()
    CTX.clear()
    CTX.update(prec=prec, stencil_pos=st_pos, cols=cols, colmap=colmap, V=V, Zinv=Zinv,
               weights=weights, taylor=taylor, n=n, r=r, cplx=cplx, slack=eps, tiny=tiny)
    hdr = dict(file=path, target=target, n=n, rank=r, p=pstr, p_float=float(p), frame=ch.frame,
               prec=prec, complex=cplx, mu=float(ch.mu), poisson_terms=len(weights),
               taylor_terms=KT, tail_bound=float(eps), g_xy=float(ch.g_xy), g_z=float(ch.g_z),
               stencil_time=t_stencil)
    if verbose:
        print('setup', hdr, flush=True)
    mpc = mp.get_context('fork')
    # ---------------- phase 1: channel, sector by sector
    t1 = time.time()
    sectors = sorted(range(n + 1), key=lambda D: -(n + 1 - D) ** 2)
    SECT.clear()
    stats1 = {}
    if nproc > 1:
        with mpc.Pool(nproc) as pool:
            for D, off, ser, st in pool.imap_unordered(_sector_job, sectors):
                SECT[D] = (off, ser)
                stats1[D] = st
    else:
        for D in sectors:
            D, off, ser, st = _sector_job(D)
            SECT[D] = (off, ser)
            stats1[D] = st
    t_channel = time.time() - t1
    if verbose:
        print(f'channel phase {t_channel:.1f} s', flush=True)
    # ---------------- phase 2: block spectra
    t2 = time.time()
    blocks = list(range(n // 2 + 1))
    res = {}
    if nproc > 1:
        with mpc.Pool(nproc) as pool:
            for out in pool.imap_unordered(_block_job, blocks):
                res[out['k']] = out
                if verbose:
                    print('  block', out['k'], 'd', out['d'], 'time %.1f' % out['time'],
                          out['sigma'][2].get('mode'), out['tau'][2].get('mode'), flush=True)
    else:
        for k in blocks:
            out = _block_job(k)
            res[k] = out
            if verbose:
                print('  block', k, out['time'], out['sigma'][2], out['tau'][2], flush=True)
    t_eig = time.time() - t2
    lo = arb(0)
    hi = arb(0)
    trace = arb(0)
    nclip = 0
    maxeta = 0.0
    for k in blocks:
        o = res[k]
        slo, shi, si = o['sigma']
        tlo, thi, ti = o['tau']
        slo, shi, tlo, thi = (RC.deser_arb(x) for x in (slo, shi, tlo, thi))
        lo += slo - thi
        hi += shi - tlo
        trace += RC.deser_arb(o['trace'])
        for info in (si, ti):
            nclip += info.get('nclip', 0)
            maxeta = max(maxeta, info.get('eta', 0.0))
    lo = lo.lower()
    hi = hi.upper()
    sign = 'POSITIVE' if lo > 0 else ('NEGATIVE' if hi < 0 else 'UNDECIDED')
    trace_ok = (trace + arb(0, (eps * r * (n + 1)))).contains(1)
    result = dict(hdr, lo=lo, hi=hi, sign=sign, width=float(hi - lo), trace=trace,
                  trace_contains_1=bool(trace_ok), nclip=nclip, max_eta=maxeta,
                  times=dict(channel=t_channel, eig=t_eig, total=time.time() - t_start),
                  blocks=res)
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('file')
    ap.add_argument('--p', default=None, help='exact p overriding the file ("a/b" or decimal)')
    ap.add_argument('--prec', type=int, default=128)
    ap.add_argument('--nproc', type=int, default=6)
    ap.add_argument('--tail', type=float, default=1e-45, help='series truncation tail target')
    ap.add_argument('--tiny', type=float, default=1e-30,
                    help='blocks with certified trace below this are bounded by their trace')
    ap.add_argument('--perturb', default=None,
                    help='negative control: "q,a,delta" adds delta to Re v_q[a]')
    ap.add_argument('--json', default=None)
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args(argv)
    if args.prec < 128:
        print('warning: precision below 128 bits', file=sys.stderr)
    pert = None
    if args.perturb:
        q, a, dl = args.perturb.split(',')
        pert = (int(q), int(a), float(dl))
    R = evaluate(args.file, args.p, args.prec, args.nproc, args.tail, args.tiny, args.verbose, pert)
    lo, hi = R['lo'], R['hi']
    print(f'file      {R["file"]}')
    print(f'family    {R["target"]}  n {R["n"]}  rank {R["rank"]}  p {R["p"]} (= {R["p_float"]:.12g})  '
          f'frame {R["frame"]}  {"complex" if R["complex"] else "real"}  prec {R["prec"]} bits')
    if pert:
        print(f'PERTURBED Re v_{pert[0]}[{pert[1]}] += {pert[2]!r}  (negative control)')
    print(f'series    Poisson mu {R["mu"]:.6f}, {R["poisson_terms"]} terms; Taylor {R["taylor_terms"]} terms; '
          f'trace-norm tail <= {R["tail_bound"]:.2e}')
    print(f'trace     sum_k tr sigma_k = {R["trace"].str(20, radius=True)}  (contains 1: {R["trace_contains_1"]})')
    print(f'spectra   max eta {R["max_eta"]:.1e}; eigenvalue intervals clipped at 0: {R["nclip"]}')
    print(f'I_c       [{lo.str(16, radius=False)}, {hi.str(16, radius=False)}] bits  (width {R["width"]:.2e})')
    print(f'I_c exact [{fmt_dec(lo, 30, up=False)}, {fmt_dec(hi, 30, up=True)}] x 1e-30 bits (outward)')
    print(f'sign      {R["sign"]}')
    t = R['times']
    print(f'time      total {t["total"]:.1f} s (channel {t["channel"]:.1f}, spectra {t["eig"]:.1f}; '
          f'{args.nproc} processes)')
    if args.json:
        out = {k: v for k, v in R.items() if k not in ('blocks',)}
        out['lo_e30'] = fmt_dec(lo, 30, up=False)
        out['hi_e30'] = fmt_dec(hi, 30, up=True)
        out['lo'] = lo.str(25, radius=False)
        out['hi'] = hi.str(25, radius=False)
        out['trace'] = R['trace'].str(25, radius=True)
        out['blocks'] = {str(k): dict(d=o['d'], time=o['time'],
                                      sigma=[RC.deser_arb(o['sigma'][0]).str(20, radius=False),
                                             RC.deser_arb(o['sigma'][1]).str(20, radius=False),
                                             o['sigma'][2]],
                                      tau=[RC.deser_arb(o['tau'][0]).str(20, radius=False),
                                           RC.deser_arb(o['tau'][1]).str(20, radius=False),
                                           o['tau'][2]])
                         for k, o in R['blocks'].items()}
        with open(args.json, 'w') as f:
            json.dump(out, f, indent=1, default=str)
    return R


if __name__ == '__main__':
    main()
