#!/usr/bin/env python3
"""Verifier for p08: noise thresholds for positive quantum capacity of qubit Pauli channels
(Table 8 item 8).

A candidate is an explicit input: r <= 4 rows v_1..v_r over the Dicke basis of Sym^n(C^2)
(rho = sum_q |v_q><v_q|, normalised by the verifier), a channel family (target) and a rational
noise level p.  The verifier certifies the sign of the coherent information

    I_c = S(N^{(x)n}(rho)) - S((N^{(x)n} (x) id_R)(|Psi><Psi|)),   |Psi> = sum_q |v_q>|q>_R,

which, when positive, proves Q(N_p) >= I_c / n > 0.  Method (PROBLEM.md, "Block formula"):
Schur-Weyl duality splits both states into sectors j = n/2 - k with multiplicity
m_j = C(n,k) - C(n,k-1); on the spin-j space the output block is the coefficient extraction of
det(Z)^k (x^dag Z y)^{2j}, Z = N(|u><w|).  Then

    I_c = (1/Z) sum_k [ H(m_j sigma_k) - H(m_j tau_k) ]          (nats; H(X) = -tr X ln X)

with Z = sum_q |v_q|^2, sigma_k the output block and tau_k the joint output-reference block.

Rigour:
  * exact rational inputs (p, the rows as the doubles given, all combinatorial scalars);
  * icball.c builds every block in ball arithmetic (centre + radius, round-to-nearest with
    explicit error terms) and encloses every eigenvalue: a floating eigendecomposition is only a
    hint; its residual and orthogonality defect are bounded with error-free (Dot2) dot products,
    and Weyl + polar-factor bounds give |lambda_i - centre_i| <= delta for the sorted spectra;
  * entropies are bounded here with decimal arithmetic (ln is correctly rounded in Python's
    decimal module) and directed rounding; the sign of I_c is decided by the enclosure.
  * depolarizing inputs with n <= RECHECK_N are evaluated twice, by the Horner path (path 0) and by
    the general Pauli dynamic programme (path 1): both enclosures must agree.
Alarms: trace of the blocks inconsistent with the exact input trace, a block eigenvalue certified
negative, the two paths disagreeing, or I_c > 0 certified where the channel is known to have zero
capacity (data/records.json: zero_capacity_from).

Verdicts (data/records.json: match_from <= bar per target):
  IMPROVEMENT     I_c > 0 certified at p >= bar
  MATCHES_RECORD  I_c > 0 certified at match_from <= p < bar
  BELOW_RECORD    I_c > 0 certified at p < match_from, or I_c < 0 certified
  UNDECIDED       the enclosure contains 0 (move p down, or improve conditioning)
  SCREEN_ONLY     --quick: same computation, no verdict
  INVALID / VERIFIER_ALARM

usage: verify.py CANDIDATE.json [--json] [--quick]
"""
import argparse
import base64
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
from decimal import ROUND_CEILING, ROUND_FLOOR, Context, Decimal, localcontext
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RECORDS = ROOT / "data" / "records.json"
SRC, BIN = HERE / "icball.c", HERE / "icball"
LOCKFILE = HERE / ".heavy.lock"
PINNED = ["verify.py", "icball.c"]
VERSION = "p08-capacity-thresholds/2"  # /2 (2026-10-27): depolarizing cap 240 -> 480

MAX_RANK = 4
NMAX = {0: 480, 1: 80}   # depolarizing Horner path / general Pauli DP
RECHECK_N = 48             # depolarizing: also run path 1 and require agreement
MAX_DEN = 10 ** 15         # denominator of p
DPREC = 40                 # decimal digits for entropy bounds
FAMILIES = ("depolarizing", "two_pauli", "independent_xz")


class Invalid(Exception):
    pass


class Alarm(Exception):
    pass


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()[:16]


def source_hashes():
    h = {f: sha(HERE / f) for f in PINNED}
    h["records.json"] = sha(RECORDS)
    if BIN.exists():
        h["icball"] = sha(BIN)
    return h


def check_pins():
    pins = HERE / "PINS.json"
    if not pins.exists():
        return "unpinned"
    want = json.loads(pins.read_text())
    cur = source_hashes()
    bad = [f for f in set(want) | set(cur) if want.get(f) != cur.get(f)]
    return "ok" if not bad else "MODIFIED: " + ",".join(sorted(bad))


def ensure_binary():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(["cc", "-O2", "-ffp-contract=off", "-fno-fast-math", "-o", str(BIN), str(SRC), "-lm"],
                       check=True)


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
    keep = {k: v for k, v in cand.items() if k in ("target", "p", "rows", "payload_c128_b64", "rank")}
    return hashlib.sha256(json.dumps(keep, sort_keys=True).encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------- exact inputs

def parse_p(x):
    if isinstance(x, bool) or not isinstance(x, (str, int)):
        raise Invalid("p must be a string 'a/b' or a decimal string (exact)")
    try:
        p = Fraction(x)
    except (ValueError, ZeroDivisionError):
        raise Invalid(f"cannot parse p={x!r}")
    if not (0 < p <= Fraction(1, 4)):
        raise Invalid("need 0 < p <= 1/4")
    if p.denominator > MAX_DEN:
        raise Invalid(f"denominator of p exceeds {MAX_DEN}")
    return p


def channel_probs(fam, p):
    if fam == "depolarizing":
        return (1 - 3 * p, p, p, p)
    if fam == "two_pauli":
        return (1 - 2 * p, p, Fraction(0), p)
    if fam == "independent_xz":
        return ((1 - p) ** 2, p * (1 - p), p * p, p * (1 - p))
    raise Invalid(f"unknown target {fam!r}; one of {FAMILIES}")


def num(x):
    if isinstance(x, bool):
        raise Invalid("booleans are not amplitudes")
    if isinstance(x, (int, float)):
        v = float(x)
    elif isinstance(x, str):
        v = float.fromhex(x) if x.strip().lower().lstrip("+-").startswith("0x") else float(x)
    else:
        raise Invalid("amplitudes must be numbers or numeric strings")
    if not math.isfinite(v):
        raise Invalid("non-finite amplitude")
    return v


def parse_rows(cand):
    if "payload_c128_b64" in cand:
        raw = base64.b64decode(cand["payload_c128_b64"], validate=True)
        r = cand.get("rank")
        if not isinstance(r, int) or r < 1 or len(raw) % (16 * r):
            raise Invalid("payload_c128_b64 needs an integer 'rank' dividing the payload")
        vals = struct.unpack("<%dd" % (len(raw) // 8), raw)
        m = len(vals) // (2 * r)
        rows = [[(vals[2 * (q * m + a)], vals[2 * (q * m + a) + 1]) for a in range(m)] for q in range(r)]
        for row in rows:
            for re, im in row:
                if not (math.isfinite(re) and math.isfinite(im)):
                    raise Invalid("non-finite amplitude")
    else:
        R = cand.get("rows")
        if not isinstance(R, list) or not R:
            raise Invalid("need 'rows': r lists of n+1 [re, im] pairs (Dicke basis)")
        rows = []
        for row in R:
            if not isinstance(row, list):
                raise Invalid("each row must be a list of [re, im]")
            out = []
            for z in row:
                if isinstance(z, (list, tuple)) and len(z) == 2:
                    out.append((num(z[0]), num(z[1])))
                else:
                    out.append((num(z), 0.0))
            rows.append(out)
    r = len(rows)
    if not 1 <= r <= MAX_RANK:
        raise Invalid(f"rank must be 1..{MAX_RANK}")
    n1 = len(rows[0])
    if n1 < 2 or any(len(row) != n1 for row in rows):
        raise Invalid("all rows need the same length n+1 >= 2")
    Z = sum(Fraction(re) ** 2 + Fraction(im) ** 2 for row in rows for re, im in row)
    if Z == 0:
        raise Invalid("zero state")
    return rows, n1 - 1, Z


# ----------------------------------------------------------------------------- balls

def ball(x):
    """double ball (mid, rad) around the exact rational x"""
    x = Fraction(x)
    mid = float(x)  # correctly rounded
    if not math.isfinite(mid):
        raise Invalid("scalar out of double range (n too large)")
    err = abs(x - Fraction(mid))
    rad = 0.0 if err == 0 else math.nextafter(float(err), math.inf)
    if not math.isfinite(rad):
        raise Invalid("scalar radius out of range")
    return mid, rad


def invsqrt_ball(N):
    with localcontext(Context(prec=50)):
        v = Decimal(1) / Decimal(N).sqrt()
    mid = float(v)
    return mid, math.ldexp(abs(mid), -50)  # >= conversion error + decimal error


def hx(x):
    return float(x).hex()


def mult(n, k):
    return math.comb(n, k) - (math.comb(n, k - 1) if k else 0)


def build_input(rows, n, path, probs):
    pI, pX, pY, pZ = probs
    cA, cB, cC, cE = pI + pZ, pX + pY, pI - pZ, pX - pY
    K = n // 2
    lines = [f"{n} {len(rows)} {path}"]

    def put(b):
        lines.append(f"{hx(b[0])} {hx(b[1])}")

    if path == 0:
        # c_s = C(d,s) cC^(d-s) cB^s exactly, as integers over a common power of Q
        P, Q = cB.numerator * cC.denominator, cB.denominator * cC.denominator
        Cn = cC.numerator * cB.denominator
        # cB = P/Q, cC = Cn/Q
        for k in range(K + 1):
            d = n - 2 * k
            Qd = Q ** d
            pc = [1] * (d + 1)
            for s in range(1, d + 1):
                pc[s] = pc[s - 1] * Cn
            pb = 1
            for s in range(d + 1):
                put(ball(Fraction(math.comb(d, s) * pc[d - s] * pb, Qd)))
                pb *= P
        AB = cA * cB
        for k in range(K + 1):
            put(ball(mult(n, k) * AB ** k))
    else:
        for c in (cA, cB, cC, cE):
            put(ball(c))
        for c in (cA * cB, cA * cA + cB * cB - cC * cC - cE * cE, cA * cB, -cC * cE, -cC * cE):
            put(ball(c))
        for a in range(n + 1):
            put(invsqrt_ball(math.comb(n, a)))
        for k in range(K + 1):
            d = n - 2 * k
            for m in range(d + 1):
                put(invsqrt_ball(math.comb(d, m)))
        for k in range(K + 1):
            put(ball(mult(n, k)))
    for row in rows:
        for re, im in row:
            lines.append(f"{re.hex()} {im.hex()}")
    return "\n".join(lines) + "\n"


def run_engine(inp, timeout):
    ensure_binary()
    proc = subprocess.run([str(BIN)], input=inp, capture_output=True, text=True, timeout=timeout)
    out = proc.stdout.split("\n")
    if proc.returncode != 0:
        msg = next((ln for ln in out if ln.startswith("FAIL")), proc.stderr.strip()[:200])
        raise RuntimeError(f"engine failed: {msg}")
    blocks, cur, i = [], None, 0
    while i < len(out):
        ln = out[i].split()
        i += 1
        if not ln:
            continue
        if ln[0] == "BLOCK":
            cur = {"k": int(ln[1]), "d": int(ln[2])}
            blocks.append(cur)
        elif ln[0] == "TR":
            cur["tr" + ln[1]] = (float.fromhex(ln[2]), float.fromhex(ln[3]))
        elif ln[0] == "EIG":
            m, emb, delta = int(ln[2]), int(ln[3]), float.fromhex(ln[4])
            vals = [float.fromhex(out[i + t]) for t in range(m)]
            i += m
            cur["eig" + ln[1]] = (emb, delta, vals)
        elif ln[0] == "END":
            break
    return blocks


# ----------------------------------------------------------------------------- entropy bounds

_LO = Context(prec=DPREC, rounding=ROUND_FLOOR)
_HI = Context(prec=DPREC, rounding=ROUND_CEILING)
_RN = Context(prec=DPREC)
_INV_E = _RN.divide(Decimal(1), _RN.exp(Decimal(1)))
INV_E_LO, INV_E_HI = _INV_E.next_minus(_RN).next_minus(_RN), _INV_E.next_plus(_RN).next_plus(_RN)


def h_bounds(x):
    """rigorous [lo, hi] for -x ln x at the exact double x >= 0"""
    if x == 0:
        return Decimal(0), Decimal(0)
    X = Decimal(x)
    L = X.ln(_RN)  # correctly rounded (half-even): true ln within one ulp
    Llo, Lhi = L.next_minus(_RN), L.next_plus(_RN)
    return _LO.minus(_HI.multiply(X, Lhi)), _HI.minus(_LO.multiply(X, Llo))


def entropy_bounds(eig, where):
    """[lo, hi] for H = -sum lambda ln lambda (nats) from sorted centres and a uniform radius"""
    emb, delta, vals = eig
    lo, hi = Decimal(0), Decimal(0)
    for c in vals:
        if c + delta < 0:
            u = math.nextafter(c + delta, math.inf)
            if u < 0:
                raise Alarm(f"{where}: eigenvalue certified negative ({c:.3e} +- {delta:.3e}) in a PSD block")
        a = max(0.0, math.nextafter(c - delta, -math.inf))
        b = max(0.0, math.nextafter(c + delta, math.inf))
        ha, hb = h_bounds(a), h_bounds(b)
        mn = min(ha[0], hb[0])
        if Decimal(b) < INV_E_LO:
            mx = hb[1]
        elif Decimal(a) > INV_E_HI:
            mx = ha[1]
        else:
            mx = INV_E_HI
        lo, hi = _LO.add(lo, mn), _HI.add(hi, mx)
    if emb:
        lo, hi = _LO.divide(lo, 2), _HI.divide(hi, 2)
    return lo, hi


def evaluate(rows, n, Z, probs, path, timeout):
    blocks = run_engine(build_input(rows, n, path, probs), timeout)
    if len(blocks) != n // 2 + 1:
        raise RuntimeError("engine output incomplete")
    lo, hi = Decimal(0), Decimal(0)
    trs_m, trs_r = Fraction(0), Fraction(0)
    for B in blocks:
        (sm, sr), (jm, jr) = B["trS"], B["trJ"]
        if abs(Fraction(sm) - Fraction(jm)) > Fraction(sr) + Fraction(jr):
            raise Alarm(f"block k={B['k']}: trace(sigma) and trace(tau) balls disjoint")
        trs_m += Fraction(sm)
        trs_r += Fraction(sr)
        Slo, Shi = entropy_bounds(B["eigS"], f"k={B['k']} sigma")
        Jlo, Jhi = entropy_bounds(B["eigJ"], f"k={B['k']} tau")
        lo = _LO.add(lo, _LO.subtract(Slo, Jhi))
        hi = _HI.add(hi, _HI.subtract(Shi, Jlo))
    if abs(trs_m - Z) > trs_r:
        raise Alarm(f"total output trace {float(trs_m)} +- {float(trs_r):.2e} != input trace {float(Z)}")
    # I_c in bits = nats / (Z ln 2)
    zlo = _LO.divide(Decimal(Z.numerator), Decimal(Z.denominator))
    zhi = _HI.divide(Decimal(Z.numerator), Decimal(Z.denominator))
    l2 = Decimal(2).ln(_RN)
    l2lo, l2hi = l2.next_minus(_RN), l2.next_plus(_RN)
    den_lo, den_hi = _LO.multiply(zlo, l2lo), _HI.multiply(zhi, l2hi)
    blo = _LO.divide(lo, den_hi if lo >= 0 else den_lo)
    bhi = _HI.divide(hi, den_lo if hi >= 0 else den_hi)
    delta_max = max(max(B["eigS"][1], B["eigJ"][1]) for B in blocks)
    return blo, bhi, delta_max


# ----------------------------------------------------------------------------- decision

def check_candidate(cand, quick=False, timeout=1800, lock=True):
    if not isinstance(cand, dict):
        raise Invalid("candidate must be a JSON object")
    fam = cand.get("target")
    if fam not in FAMILIES:
        raise Invalid(f"target must be one of {FAMILIES}")
    p = parse_p(cand.get("p"))
    rows, n, Z = parse_rows(cand)
    probs = channel_probs(fam, p)
    path = 0 if fam == "depolarizing" else 1
    if n > NMAX[path]:
        raise Invalid(f"n={n} exceeds the verifier cap {NMAX[path]} for this target")
    cell = json.loads(RECORDS.read_text())["cells"][fam]
    t0 = time.time()
    with heavy_lock(lock and not quick):
        lo, hi, dmax = evaluate(rows, n, Z, probs, path, timeout)
        rechecked = None
        if path == 0 and n <= RECHECK_N and not quick:
            lo2, hi2, _ = evaluate(rows, n, Z, probs, 1, timeout)
            if lo2 > hi or lo > hi2:
                raise Alarm(f"paths disagree: [{lo:.6e}, {hi:.6e}] vs [{lo2:.6e}, {hi2:.6e}]")
            if (lo > 0) != (lo2 > 0) or (hi < 0) != (hi2 < 0):
                rechecked = "path-1 recheck could not confirm the sign"
                lo, hi = max(lo, lo2), min(hi, hi2)  # intersection is still rigorous
            else:
                rechecked = "path 1 agrees"
    secs = time.time() - t0
    pd = f"{float(p):.9f}".rstrip("0")
    enc = f"I_c in [{lo:.4e}, {hi:.4e}] bits"
    info = {"target": fam, "p": str(p), "n": n, "rank": len(rows), "I_c_lo": str(lo), "I_c_hi": str(hi),
            "delta_max": dmax, "engine_seconds": round(secs, 2)}
    if rechecked:
        info["recheck"] = rechecked
    tail = f"{fam} n={n} r={len(rows)} p={p} ({pd}): {enc}"
    if quick:
        return dict(info, verdict="SCREEN_ONLY", summary=f"SCREEN_ONLY: {tail}")
    if lo > 0:
        if p >= Fraction(cell["zero_capacity_from"]):
            raise Alarm(f"I_c > 0 certified at p={p} where Q = 0 ({cell['zero_capacity_reason']})")
        info["certified_p"] = float(p)  # board score: best p with certified I_c > 0
        bar, mfrom = Fraction(cell["bar"]), Fraction(cell["match_from"])
        if p >= bar:
            v = "IMPROVEMENT"
            note = f"Q > 0 certified at p={pd} >= bar {float(bar)}"
        elif p >= mfrom:
            v = "MATCHES_RECORD"
            note = f"Q > 0 certified at p={pd} (record window [{float(mfrom)}, {float(bar)}))"
        else:
            v = "BELOW_RECORD"
            note = f"Q > 0 certified at p={pd} < {float(mfrom)}"
        return dict(info, verdict=v, summary=f"{v}: {note}; {tail}")
    if hi < 0:
        return dict(info, verdict="BELOW_RECORD", summary=f"BELOW_RECORD: I_c < 0 certified; {tail}")
    return dict(info, verdict="UNDECIDED", reason="enclosure contains 0",
                summary=f"UNDECIDED: enclosure contains 0; {tail}")


def verify_candidate(cand, exact=False, quick=False, threads=8, timeout=1800, budget=0, seed=1,
                     lock=True, **kw):
    t0 = time.time()
    rep = {"verifier": VERSION, "pins": check_pins(), "sources": source_hashes()}
    try:
        rep["hash"] = canonical_hash(cand) if isinstance(cand, dict) else None
        rep.update(check_candidate(cand, quick=quick, timeout=timeout, lock=lock))
    except Invalid as e:
        rep.update(verdict="INVALID", reason=str(e))
    except Alarm as e:
        rep.update(verdict="VERIFIER_ALARM", reason=str(e))
    except subprocess.TimeoutExpired:
        rep.update(verdict="UNDECIDED", reason=f"engine timeout after {timeout}s")
    except RuntimeError as e:  # engine could not certify (e.g. eigensolver), not a bad candidate
        rep.update(verdict="UNDECIDED", reason=str(e))
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
    ap.add_argument("--quick", action="store_true", help="screen only (no verdict)")
    a = ap.parse_args()
    rep = verify_candidate(json.loads(pathlib.Path(a.candidate).read_text()), quick=a.quick)
    print(json.dumps(rep, indent=1) if a.json else rep["summary"])
    sys.exit(0 if rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") else 1)


if __name__ == "__main__":
    main()
