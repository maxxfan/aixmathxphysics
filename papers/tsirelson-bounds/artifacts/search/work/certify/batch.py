#!/usr/bin/env python3
"""certify pipeline for one catalog inequality (or a list).

  nice -n 10 .venv/bin/python work/certify/batch.py NAME [NAME ...] [--part lower|upper|both]
        [--restarts 20] [--seed 1] [--den 14] [--maxwords 150] [--submit] [--verify]

lower: seesaw in PV's (field, dim) from records.json (fixed seeds), rationalised with den 10**DEN.
upper: NPA dual (Clarabel) at PV's level, escalating through LEVELS until the estimated window < GAP;
       Y = Y' + delta*1 with Y' >= 0 solved for the shifted class targets, delta chosen by a sweep.
Candidates go to work/certify/cands/, one JSON log line per run to work/certify/log.jsonl.
"""
import argparse
import json
import os
import subprocess
import sys
import time

os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
sys.path.insert(0, "qlib")
import numpy as np
import scipy.sparse as sp

import bell
import certs

WORK = "work/certify"
BOARD = "<quantum-information>/core/board/board.py"
LEVELS = ["1", "1+ab", "2", "2+aab", "2+abb", "2+aab+abb", "3", "3+aabb"]
PVLEV = {"1": "1", "1a": "1+ab", "2": "2", "2a": "2+aab", "2b": "2+aab+abb", "3": "3", "L3": "3", "L3+aabb": "3+aabb",
         "no": "3", None: "1+ab"}


def log(rec):
    rec["t"] = time.strftime("%H:%M:%S")
    with open(f"{WORK}/log.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


# ------------------------------------------------------------------------------------ lower
def do_lower(name, c, rec, restarts, seed, den):
    field = rec.get("pv_field") or "real"
    d = rec.get("pv_dim") or 2
    t0 = time.time()
    best, st = bell.seesaw(c, d, field, restarts=restarts, iters=3000, seed=seed, tol=1e-15)
    cand = certs.strategy_candidate(c, st, field=field, den=10 ** den, author="certify",
                                    construction=f"seesaw d={d} {field} restarts={restarts} seed={seed}",
                                    script=f"work/certify/batch.py {name} --part lower --seed {seed}")
    path = f"{WORK}/cands/{name}_lower.json"
    certs.save(cand, path)
    log(dict(name=name, part="lower", value=best, d=d, field=field, path=path, secs=round(time.time() - t0, 1)))
    return best, path


# ------------------------------------------------------------------------------------ upper
def dual_shifted(c, W, delta, solver="CLARABEL", tol=1e-11):
    """min Y'_00 s.t. Y' >= 0, A(Y') = t - delta*A(1); returns (Y = Y' + delta*1, lambda - b0)."""
    import cvxpy as cp
    n = len(W)
    cl, inv = certs.classes(W)
    one = cl[0, 0]
    t = certs.targets(c, inv)
    rows = cl.reshape(-1)
    A = sp.csr_matrix((np.ones(n * n), (rows, np.arange(n * n))), shape=(len(inv), n * n))
    dcount = np.bincount(np.diag(cl), minlength=len(inv)).astype(float)
    keep = np.array([k for k in range(len(inv)) if k != one])
    Yv = cp.Variable((n, n), PSD=True)
    prob = cp.Problem(cp.Minimize(Yv[0, 0]), [A[keep] @ cp.vec(Yv, order="C") == (t - delta * dcount)[keep]])
    prob.solve(solver=solver, tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol, max_iter=500)
    Y = np.array(Yv.value)
    Y = (Y + Y.T) / 2 + delta * np.eye(n)
    return Y, Y[0, 0] - c["b0"]


def est_bound(cand, n, b0):
    D2 = int(cand["Y"]["denominator"])
    M = np.zeros((n, n))
    for i, row in enumerate(cand["Y"]["upper"]):
        for k, e in enumerate(row):
            M[i, i + k] = M[i + k, i] = int(e) / D2
    lmin = np.linalg.eigvalsh(M)[0]
    lam = int(cand["Y"]["upper"][0][0]) / D2
    return lam - b0 + n * max(0.0, -lmin) * 1.5, lmin


def upper_at(name, c, level, deltas, ref=None, gap=1.2e-7):
    W = bell.words(c, int(level) if level.isdigit() else level)
    n = len(W)
    best = None
    for delta in deltas:
        if best is not None and ref is not None and best[0] - ref < gap:
            break
        Y, lam = dual_shifted(c, W, delta)
        cand = certs.sos_candidate(c, W, Y, eps=0.0, den=2 ** 46, author="certify",
                                   construction=f"NPA dual level {level} (N={n}), Clarabel, shift {delta:g}",
                                   script=f"work/certify/batch.py {name} --part upper")
        ub, lmin = est_bound(cand, n, c["b0"])
        if best is None or ub < best[0]:
            best = (ub, cand, delta, lmin, lam)
    return n, best


def do_upper(name, c, rec, lower, maxwords, gap, deltas, start=None):
    lev0 = start or PVLEV.get(rec.get("pv_level"), "1+ab")
    i0 = LEVELS.index(lev0)
    out = None
    for level in LEVELS[i0:]:
        W = bell.words(c, int(level) if level.isdigit() else level)
        if len(W) > maxwords:
            log(dict(name=name, part="upper", level=level, n=len(W), note="skipped: too many words"))
            break
        t0 = time.time()
        ref0 = float(rec["lower"]) if rec.get("lower") else lower
        n, (ub, cand, delta, lmin, lam) = upper_at(name, c, level, deltas, ref0, gap)
        path = f"{WORK}/cands/{name}_upper_{level.replace('+', '')}.json"
        certs.save(cand, path)
        ref = float(rec["lower"]) if rec.get("lower") else lower   # PV's 7-digit value is the target
        win = None if ref is None else ub - ref
        log(dict(name=name, part="upper", level=level, n=n, lam=lam, est_ub=ub, delta=delta, lmin=lmin,
                 est_window=win, path=path, secs=round(time.time() - t0, 1)))
        out = (ub, path, level)
        if win is not None and win < gap:
            break
    return out


# ------------------------------------------------------------------------------------ verify / submit
def verify(path):
    t0 = time.time()
    p = subprocess.run([sys.executable, "verifier/verify.py", path, "--json"], capture_output=True, text=True)
    try:
        rep = json.loads(p.stdout)
    except Exception:
        rep = {"verdict": "ERROR", "stderr": p.stderr[-500:], "stdout": p.stdout[-300:]}
    keep = {k: rep.get(k) for k in ("verdict", "side", "violation", "bound", "penalty", "tau", "n_words", "summary")}
    keep["secs"] = round(time.time() - t0, 1)
    return keep


def submit(path):
    p = subprocess.run(["python3", BOARD, "submit", path, "--as", "certify", "--quiet"], capture_output=True, text=True)
    return (p.stdout + p.stderr)[-600:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--part", default="both")
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--den", type=int, default=14)
    ap.add_argument("--maxwords", type=int, default=150)
    ap.add_argument("--gap", type=float, default=1.2e-7)   # vs PV's 7-digit value (rounding +-5e-8)
    ap.add_argument("--start", default=None, help="force starting level string")
    ap.add_argument("--deltas", default="1e-11")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--submit", action="store_true")
    a = ap.parse_args()
    os.makedirs(f"{WORK}/cands", exist_ok=True)
    deltas = [float(x) for x in a.deltas.split(",")]
    for name in a.names:
        try:
            c = bell.load(name)
            rec = bell.record(name)
            lower, lpath, upath = None, None, None
            if a.part in ("lower", "both"):
                lower, lpath = do_lower(name, c, rec, a.restarts, a.seed, a.den)
            if a.part in ("upper", "both"):
                res = do_upper(name, c, rec, lower, a.maxwords, a.gap, deltas, a.start)
                upath = res[1] if res else None
            for p in (lpath, upath):
                if p and (a.verify or a.submit):
                    rep = verify(p)
                    log(dict(name=name, verify=p, **rep))
                    if a.submit and rep.get("verdict") not in ("INVALID", "VERIFIER_ALARM", "ERROR"):
                        log(dict(name=name, submit=p, out=submit(p)))
        except Exception as e:  # keep the batch going
            log(dict(name=name, error=repr(e)))


if __name__ == "__main__":
    main()
