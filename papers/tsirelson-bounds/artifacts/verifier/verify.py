#!/usr/bin/env python3
"""Verifier for p05: Tsirelson bounds of specific Bell inequalities (Table 8 item 5).

Two candidate types, both checked exactly (see PROBLEM.md for the formats):

  strategy   a finite-dimensional quantum strategy: for each setting a projector given by exact
             (rational or Gaussian-rational) spanning vectors, and an exact state vector.  The
             projector is V (V^dag V)^{-1} V^dag, computed exactly, so it is a projector by
             construction.  The violation <psi|W|psi>/<psi|psi> - b0 is an exact rational, computed
             twice (Alice-side reduced-state formula; then an entrywise Kronecker sum for
             dA*dB <= 64, else the Bob-side reduced formula) -> a certified LOWER bound.

  sos        a sum-of-squares (NPA dual) certificate: words w_1..w_N in the projectors and an exact
             symmetric matrix Y with  sum_ij Y_ij w_i^dag w_j = lambda*1 - W  as an identity in the
             universal algebra (projectors idempotent, Alice commutes with Bob), checked
             coefficient by coefficient in the basis of reduced words.  Y + tau*1 >= 0 is proved
             rigorously (rescheck.c: floating Cholesky L of Y + sigma*1, interval enclosure of the
             residual F = Y + sigma*1 - L L^T, tau = sigma + bound on ||F||_F); for small N
             Y + 2*tau*1 is also checked exactly (fraction-free elimination).  Since the moment matrix of any state has
             trace <= N (each w_i has norm <= 1), every quantum strategy (commuting operators, any
             dimension) has  <W> <= lambda + N*tau  -> a certified UPPER bound.

Records come from data/records.json (Pal-Vertesi 2009).  With delta = 1e-6:
  strategy  IMPROVEMENT if violation >= record_lower + delta; VERIFIER_ALARM if it exceeds the
            published upper bound by more than delta (contradiction).
  sos       IMPROVEMENT if bound <= record_upper - delta; VERIFIER_ALARM if it is below the published
            lower bound (an explicit strategy) by more than delta.
  MATCHES_RECORD within delta of the record, BELOW_RECORD otherwise; INVALID for malformed input.

usage: verify.py CANDIDATE.json [--json]
"""
import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import pathlib
import struct
import subprocess
import sys
import time
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import qexact as qx  # noqa: E402

CATALOG = ROOT / "data" / "catalog.json"
RECORDS = ROOT / "data" / "records.json"
RES_SRC, RES_BIN = HERE / "rescheck.c", HERE / "rescheck"
LOCKFILE = HERE / ".heavy.lock"
PINNED = ["verify.py", "qexact.py", "rescheck.c"]
DELTA = Fraction(1, 10 ** 6)
MAX_DIM = 32          # local dimension cap for strategies (PV 2010 used up to 25)
MAX_WORDS = 2000      # SOS basis size cap
KRON_MAX = 64         # entrywise Kronecker re-check up to dA*dB = 64, Bob-side reduced formula above
EXACT_PD_MAX_N = 120  # also run the exact PD test up to this size
SIGMA_LADDER = [0.0] + [2.0 ** -k for k in range(60, -3, -2)]  # Cholesky shifts tried in rescheck
SIGMA_BISECT = 8       # geometric refinement steps between the last failing and first good shift


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()[:16]


def source_hashes():
    h = {f: sha(HERE / f) for f in PINNED}
    h["catalog.json"] = sha(CATALOG)
    h["records.json"] = sha(RECORDS)
    if RES_BIN.exists():
        h["rescheck"] = sha(RES_BIN)
    return h


def ensure_rescheck():
    if not RES_BIN.exists() or RES_BIN.stat().st_mtime < RES_SRC.stat().st_mtime:
        subprocess.run(["cc", "-O2", "-ffp-contract=off", "-fno-fast-math", "-o", str(RES_BIN),
                        str(RES_SRC), "-lm"], check=True)


def check_pins():
    pins = HERE / "PINS.json"
    if not pins.exists():
        return "unpinned"
    want = json.loads(pins.read_text())
    cur = source_hashes()
    bad = [f for f in set(want) | set(cur) if want.get(f) != cur.get(f)]
    return "ok" if not bad else "MODIFIED: " + ",".join(sorted(bad))


@contextlib.contextmanager
def heavy_lock(enabled=True):
    if not enabled:
        yield
        return
    with open(LOCKFILE, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def canonical_hash(cand):
    keep = {k: v for k, v in cand.items() if k not in ("author", "construction", "script", "note", "notes",
                                                         "claimed", "claimed_value", "level", "created")}
    return hashlib.sha256(json.dumps(keep, sort_keys=True).encode()).hexdigest()[:16]


def dec_down(q, digits=10):
    s = 10 ** digits
    return f"{math.floor(q * s) / s:.{digits}f}"


def dec_up(q, digits=10):
    s = 10 ** digits
    return f"{math.ceil(q * s) / s:.{digits}f}"


def load_problem(name):
    cat = json.loads(CATALOG.read_text())["inequalities"]
    recs = json.loads(RECORDS.read_text())["records"]
    if not isinstance(name, str) or name not in cat:
        raise qx.Invalid(f"unknown inequality {name!r} (see data/catalog.json)")
    return cat[name], recs.get(name)


# ----------------------------------------------------------------------------- strategy

def check_strategy(cand):
    ineq, rec = load_problem(cand.get("inequality"))
    F = qx.Field(cand.get("field", "real"))
    dA, dB = cand.get("dA"), cand.get("dB")
    if not (isinstance(dA, int) and isinstance(dB, int) and 1 <= dA <= MAX_DIM and 1 <= dB <= MAX_DIM):
        raise qx.Invalid(f"dA, dB must be integers in 1..{MAX_DIM}")
    As, Bs, psi = cand.get("A"), cand.get("B"), cand.get("psi")
    if not (isinstance(As, list) and len(As) == ineq["mA"] and isinstance(Bs, list) and len(Bs) == ineq["mB"]):
        raise qx.Invalid(f"need {ineq['mA']} Alice and {ineq['mB']} Bob projectors (lists of spanning vectors)")
    if not (isinstance(psi, list) and len(psi) == dA * dB):
        raise qx.Invalid(f"psi must have dA*dB = {dA * dB} entries")

    def proj(spec, d):
        if not isinstance(spec, list) or len(spec) > d or any(not isinstance(v, list) for v in spec):
            raise qx.Invalid("each projector is a list of at most d spanning vectors")
        return qx.projector_from_span(F, [[F.parse(e) for e in v] for v in spec], d)
    A = [proj(s, dA) for s in As]
    B = [proj(s, dB) for s in Bs]
    for P in A + B:
        if not qx.is_projector(F, P):
            return {"verdict": "VERIFIER_ALARM", "reason": "constructed operator is not a projector"}
    ps = [F.parse(e) for e in psi]
    v1 = qx.strategy_value(F, ineq, A, B, ps, dA, dB)
    if dA * dB <= KRON_MAX:
        v2, how = qx.strategy_value_kron(F, ineq, A, B, ps, dA, dB), "kron"
    else:
        v2, how = qx.strategy_value_bob(F, ineq, A, B, ps, dA, dB), "bob-side"
    if v1 != v2:
        return {"verdict": "VERIFIER_ALARM", "reason": f"the two exact evaluations disagree ({how})"}
    viol = v1 - ineq["b0"]
    out = {"type": "strategy", "inequality": cand["inequality"], "side": "lower",
           "target": f"{cand['inequality']}:lower", "dims": [dA, dB], "field": F.kind,
           "ranks": [[len(s) for s in As], [len(s) for s in Bs]], "recheck": how,
           "violation": dec_down(viol), "violation_exact": str(viol) if len(str(viol)) < 400 else None,
           "value_float": float(viol)}
    return judge(out, viol, rec, "lower")


# ----------------------------------------------------------------------------- SOS certificate

def interval(q):
    """Smallest double interval containing the rational q."""
    f = float(q)  # correctly rounded
    if math.isinf(f) or math.isnan(f):
        raise qx.Invalid("matrix entry out of double range")
    fq = Fraction(f)
    if fq == q:
        return f, f
    if fq > q:
        return math.nextafter(f, -math.inf), f
    return f, math.nextafter(f, math.inf)


def packed_lower(Ynum, den):
    """Lower triangle of Ynum/den as enclosing double intervals (binary, for rescheck)."""
    n = len(Ynum)
    buf = bytearray()
    for i in range(n):
        row = Ynum[i]
        for j in range(i + 1):
            lo, hi = interval(Fraction(row[j], den))
            buf += struct.pack("dd", lo, hi)
    return bytes(buf)


def rescheck(n, packed, sigma):
    """Run rescheck.c with shift sigma; returns f (a float) or None if the Cholesky failed."""
    ensure_rescheck()
    p = subprocess.run([str(RES_BIN)], input=struct.pack("i", n) + struct.pack("d", sigma) + packed,
                       capture_output=True, timeout=3600)
    line = p.stdout.decode().strip()
    if line.startswith("OK "):
        f = float(line.split()[1])
        if math.isfinite(f) and f >= 0:
            return f
    return None


def psd_shift(Ynum, den):
    """Smallest tau found with a rigorous proof that Y + tau*1 >= 0 (tau = sigma + f, see rescheck.c)."""
    n = len(Ynum)
    packed = packed_lower(Ynum, den)
    runs, best, prev = 0, None, None
    for sg in SIGMA_LADDER:
        f = rescheck(n, packed, sg)
        runs += 1
        if f is not None:
            best = (Fraction(sg) + Fraction(f), sg, f)
            break
        prev = sg
    if best is None:
        return None, runs
    if prev is not None and best[1] > 0:
        lo, hi = (prev if prev > 0 else best[1] / 1024), best[1]
        for _ in range(SIGMA_BISECT):
            mid = math.sqrt(lo * hi)
            if not lo < mid < hi:
                break
            f = rescheck(n, packed, mid)
            runs += 1
            if f is None:
                lo = mid
            else:
                hi = mid
                cand = Fraction(mid) + Fraction(f)
                if cand < best[0]:
                    best = (cand, mid, f)
    return best, runs


def check_sos(cand, lock=True):
    ineq, rec = load_problem(cand.get("inequality"))
    mA, mB = ineq["mA"], ineq["mB"]
    ws = cand.get("words")
    if not isinstance(ws, list) or not 1 <= len(ws) <= MAX_WORDS:
        raise qx.Invalid(f"words must be a list of 1..{MAX_WORDS} strings")
    words = [qx.parse_word(w, mA, mB) for w in ws]
    if len(set(words)) != len(words):
        raise qx.Invalid("duplicate words (after normal form)")
    n = len(words)
    Yspec = cand.get("Y")
    if not isinstance(Yspec, dict):
        raise qx.Invalid("Y must be {'denominator': D, 'upper': [...]} or {'denominator': D, 'entries': [[i,j,num],...]}")
    den = qx.parse_q(Yspec.get("denominator", 1))
    if den.denominator != 1 or den <= 0:
        raise qx.Invalid("Y.denominator must be a positive integer")
    den = int(den)
    Ynum = [[0] * n for _ in range(n)]
    if "upper" in Yspec:
        up = Yspec["upper"]
        if not isinstance(up, list) or len(up) != n:
            raise qx.Invalid("Y.upper must have one row per word")
        for i, row in enumerate(up):
            if not isinstance(row, list) or len(row) != n - i:
                raise qx.Invalid(f"Y.upper row {i} must have {n - i} entries (j = i..N-1)")
            for k, e in enumerate(row):
                v = qx.parse_q(e)
                if v.denominator != 1:
                    raise qx.Invalid("Y entries are integer numerators over Y.denominator")
                Ynum[i][i + k] = Ynum[i + k][i] = int(v)
    elif "entries" in Yspec:
        seen = set()
        for t in Yspec["entries"]:
            if not (isinstance(t, list) and len(t) == 3 and all(isinstance(x, int) for x in t[:2])):
                raise qx.Invalid("Y.entries items are [i, j, numerator]")
            i, j = t[0], t[1]
            if not (0 <= i < n and 0 <= j < n) or (min(i, j), max(i, j)) in seen:
                raise qx.Invalid(f"bad or repeated entry index {i},{j}")
            seen.add((min(i, j), max(i, j)))
            v = qx.parse_q(t[2])
            if v.denominator != 1:
                raise qx.Invalid("Y entries are integer numerators over Y.denominator")
            Ynum[i][j] = Ynum[j][i] = int(v)
    else:
        raise qx.Invalid("Y needs 'upper' or 'entries'")

    # 1. the identity sum_ij Y_ij w_i^dag w_j = lambda - W, coefficient by coefficient (times den)
    coef = qx.sos_identity(words, Ynum)
    target = {m: -c * den for m, c in qx.bell_polynomial(ineq).items()}
    lam_num = coef.get(((), ()), 0)
    lam = Fraction(lam_num, den)
    if "lambda" in cand and qx.parse_q(cand["lambda"]) != lam:
        raise qx.Invalid(f"lambda {cand['lambda']} does not match the identity coefficient {lam}")
    bad = []
    for m in set(coef) | set(target):
        if m == ((), ()):
            continue
        if coef.get(m, 0) != target.get(m, 0):
            bad.append(m)
    if bad:
        m = sorted(bad)[0]
        return {"verdict": "INVALID", "reason": f"SOS identity fails on {len(bad)} monomial(s), e.g. "
                f"{fmt_word(m)}: have {Fraction(coef.get(m, 0), den)}, need {Fraction(target.get(m, 0), den)}",
                "type": "sos", "inequality": cand["inequality"], "target": f"{cand['inequality']}:upper"}

    # 2. a rigorous tau with Y + tau*1 >= 0 (Cholesky + interval residual bound, rescheck.c)
    with heavy_lock(lock):
        best, runs = psd_shift(Ynum, den)
    if best is None:
        return {"verdict": "INVALID", "reason": f"Y + sigma*1 has no floating Cholesky factor for sigma <= {SIGMA_LADDER[-1]}",
                "type": "sos", "inequality": cand["inequality"], "target": f"{cand['inequality']}:upper"}
    tau, sigma, fres = best
    exact_pd = None
    if n <= EXACT_PD_MAX_N and tau > 0:
        # independent exact check that Y + 2*tau*1 is positive definite (integer Bareiss elimination)
        tq = 2 * tau * den
        scale = tq.denominator
        M = [[Ynum[i][j] * scale + (tq.numerator if i == j else 0) for j in range(n)] for i in range(n)]
        exact_pd = qx.ldl_pd(M)
        if not exact_pd:
            return {"verdict": "VERIFIER_ALARM", "reason": "residual bound and exact elimination disagree",
                    "type": "sos", "inequality": cand["inequality"], "target": f"{cand['inequality']}:upper"}
    bound = lam + n * tau - ineq["b0"]
    out = {"type": "sos", "inequality": cand["inequality"], "side": "upper",
           "target": f"{cand['inequality']}:upper", "n_words": n, "tau": float(tau), "penalty": float(n * tau),
           "lambda": dec_up(lam), "bound": dec_up(bound),
           "bound_exact": str(bound) if len(str(bound)) < 400 else None, "value_float": float(bound),
           "sigma": sigma, "residual_bound": fres, "psd_runs": runs, "exact_pd_crosscheck": exact_pd}
    return judge(out, bound, rec, "upper")


def fmt_word(m):
    return "".join(f"A{x + 1}" for x in m[0]) + "".join(f"B{y + 1}" for y in m[1]) or "1"


# ----------------------------------------------------------------------------- verdict

def judge(out, val, rec, side):
    name = out["inequality"]
    if rec is None:
        out.update(verdict="NO_TABLE_ENTRY", summary=f"NO_TABLE_ENTRY: {name} {side} {out.get('violation') or out.get('bound')}")
        return out
    lo, hi = Fraction(rec["lower"]), Fraction(rec["upper"])
    out["record"] = {"lower": rec["lower"], "upper": rec["upper"], "status": rec["status"]}
    if side == "lower":
        gain = val - lo
        out["gain"] = float(gain)
        if val > hi + DELTA:
            v, why = "VERIFIER_ALARM", f"strategy value exceeds the published NPA upper bound {rec['upper']} by more than 1e-6"
        elif gain >= DELTA and rec["status"] == "open":
            v, why = "IMPROVEMENT", None
        elif gain > -DELTA:
            v, why = "MATCHES_RECORD", None
        else:
            v, why = "BELOW_RECORD", None
        shown = out["violation"]
    else:
        gain = hi - val
        out["gain"] = float(gain)
        if val < lo - DELTA:
            v, why = "VERIFIER_ALARM", f"certified upper bound is below the published strategy value {rec['lower']} by more than 1e-6"
        elif gain >= DELTA and rec["status"] == "open":
            v, why = "IMPROVEMENT", None
        elif gain > -DELTA:
            v, why = "MATCHES_RECORD", None
        else:
            v, why = "BELOW_RECORD", None
        shown = out["bound"]
    out["verdict"] = v
    if why:
        out["reason"] = why
    rel = ">=" if side == "lower" else "<="
    out["summary"] = (f"{v}: {name} quantum max {rel} {shown} (record {side} {rec[side]}, window "
                      f"[{rec['lower']}, {rec['upper']}], gain {out['gain']:+.2e})")
    return out


def check_candidate(cand, lock=True, **opts):
    if not isinstance(cand, dict):
        raise qx.Invalid("candidate must be a JSON object")
    kind = cand.get("type")
    if kind == "strategy":
        return check_strategy(cand)
    if kind == "sos":
        return check_sos(cand, lock=lock)
    raise qx.Invalid("type must be 'strategy' or 'sos'")


def verify_candidate(cand, exact=False, quick=False, threads=8, timeout=1800, budget=0, seed=1,
                     lock=True, **kw):
    t0 = time.time()
    rep = {"verifier": "p05-tsirelson-bounds/1", "pins": check_pins(), "sources": source_hashes()}
    try:
        rep["hash"] = canonical_hash(cand) if isinstance(cand, dict) else None
        rep.update(check_candidate(cand, lock=lock))
    except qx.Invalid as e:
        rep.update(verdict="INVALID", reason=str(e))
    except Exception as e:  # never let bad input crash the board
        rep.update(verdict="INVALID", reason=f"{type(e).__name__}: {e}")
    rep["seconds"] = round(time.time() - t0, 2)
    if "summary" not in rep:
        rep["summary"] = f"{rep['verdict']}: {rep.get('reason', '')}"[:300]
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--quick", action="store_true", help="accepted for compatibility; checks are always exact")
    a = ap.parse_args()
    rep = verify_candidate(json.loads(pathlib.Path(a.candidate).read_text()))
    print(json.dumps(rep, indent=1) if a.json else rep["summary"])
    sys.exit(0 if rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") else 1)


if __name__ == "__main__":
    main()
