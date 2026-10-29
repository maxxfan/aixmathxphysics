#!/usr/bin/env python3
"""p08 verifier tests (run before any search; freeze.sh runs them).

1. malformed candidates -> INVALID
2. an independent brute force (full 2^n density matrices, pure-Python Jacobi) lies inside the
   certified enclosure for random small inputs, every channel family, both engine paths
3. the four certified endpoints of arXiv:2608.15870 reproduce (A, R at n = 45)
4. every verdict path: MATCHES_RECORD, BELOW_RECORD (positive and negative), IMPROVEMENT and
   VERIFIER_ALARM (doctored records), UNDECIDED (rank one: I_c = 0 exactly), SCREEN_ONLY,
   alarms for path disagreement, a negative eigenvalue and an inconsistent trace
5. the arXiv:2605.09138 states (if data/fetch_public_states.py has been run)

run: python3 verifier/tests/test_verifier.py [--fast]
"""
import base64
import cmath
import json
import math
import pathlib
import random
import sys
import tempfile
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))
import verify  # noqa: E402

FAST = "--fast" in sys.argv
FAIL = []


def check(name, cond, info=""):
    print(("ok   " if cond else "FAIL ") + name + (f"  [{info}]" if info and not cond else ""))
    if not cond:
        FAIL.append(name)


# ----------------------------------------------------------------------------- brute force

def jacobi_eigs(A):
    """eigenvalues of a real symmetric matrix (list of lists), cyclic Jacobi"""
    n = len(A)
    A = [row[:] for row in A]
    for _ in range(100):
        off = sum(A[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-30:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(A[p][q]) < 1e-300:
                    continue
                th = (A[q][q] - A[p][p]) / (2 * A[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p], A[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k], A[q][k] = c * apk - s * aqk, s * apk + c * aqk
    return [A[i][i] for i in range(n)]


def herm_entropy_nats(M):
    n = len(M)
    E = [[0.0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            z = (M[i][j] + M[j][i].conjugate()) / 2
            E[i][j] = E[i + n][j + n] = z.real
            E[i][j + n] = -z.imag
            E[i + n][j] = z.imag
    return sum(-x * math.log(x) for x in jacobi_eigs(E) if x > 1e-300) / 2


def brute_ic_bits(rows, probs):
    r, n = len(rows), len(rows[0]) - 1
    dim = 2 ** n
    pops = [bin(x).count("1") for x in range(dim)]
    psi = [[rows[q][pops[x]] / math.sqrt(math.comb(n, pops[x])) for q in range(r)] for x in range(dim)]
    D = dim * r
    rho = [[psi[x][q] * psi[y][qq].conjugate() for y in range(dim) for qq in range(r)]
           for x in range(dim) for q in range(r)]
    pI, pX, pY, pZ = [float(t) for t in probs]
    for i in range(n):
        bit = 1 << i
        new = [[pI * rho[a][b] for b in range(D)] for a in range(D)]
        for pr, kind in ((pX, "X"), (pY, "Y"), (pZ, "Z")):
            if pr == 0:
                continue
            img, ph = [0] * D, [0j] * D
            for x in range(dim):
                b = 1 if x & bit else 0
                if kind == "X":
                    tx, f = x ^ bit, 1
                elif kind == "Z":
                    tx, f = x, (-1 if b else 1)
                else:
                    tx, f = x ^ bit, (-1j if b else 1j)
                for q in range(r):
                    img[x * r + q], ph[x * r + q] = tx * r + q, f
            for a in range(D):
                ia, fa = img[a], ph[a]
                for b in range(D):
                    new[ia][img[b]] += pr * fa * ph[b].conjugate() * rho[a][b]
        rho = new
    rhoB = [[sum(rho[x * r + q][y * r + q] for q in range(r)) for y in range(dim)] for x in range(dim)]
    Z = sum(rhoB[x][x].real for x in range(dim))
    rhoB = [[z / Z for z in row] for row in rhoB]
    rho = [[z / Z for z in row] for row in rho]
    return (herm_entropy_nats(rhoB) - herm_entropy_nats(rho)) / math.log(2)


# ----------------------------------------------------------------------------- helpers

def payload_cand(path, p, target="depolarizing"):
    raw = pathlib.Path(path).read_bytes()
    return {"target": target, "p": p, "rank": 2, "payload_c128_b64": base64.b64encode(raw).decode()}


def rand_rows(rng, n, r, real=False):
    return [[complex(rng.gauss(0, 1), 0 if real else rng.gauss(0, 1)) for _ in range(n + 1)] for _ in range(r)]


def as_json_rows(rows):
    return [[[z.real, z.imag] for z in row] for row in rows]


def with_records(mod):
    """context: a doctored copy of data/records.json"""
    class C:
        def __enter__(self):
            d = json.loads(verify.RECORDS.read_text())
            mod(d)
            self.tmp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
            json.dump(d, self.tmp)
            self.tmp.close()
            self.old = verify.RECORDS
            verify.RECORDS = pathlib.Path(self.tmp.name)

        def __exit__(self, *a):
            verify.RECORDS = self.old
            pathlib.Path(self.tmp.name).unlink()
    return C()


# ----------------------------------------------------------------------------- 1. invalid input

KG = ROOT / "data" / "sources" / "kg2026"
A, R = KG / "A.c128", KG / "R.c128"
good = payload_cand(A, "16239/250000")
bad = [
    ("garbage", {"garbage": True}),
    ("not a dict", [1, 2]),
    ("unknown target", dict(good, target="amplitude_damping")),
    ("p as float", dict(good, p=0.064956)),
    ("p out of range", dict(good, p="3/10")),
    ("p zero", dict(good, p="0")),
    ("p huge denominator", dict(good, p="1/" + "9" * 20)),
    ("ragged rows", {"target": "depolarizing", "p": "1/16", "rows": [[[1, 0], [0, 0]], [[1, 0]]]}),
    ("rank 5", {"target": "depolarizing", "p": "1/16", "rows": [[[1, 0], [0, 1]]] * 5}),
    ("zero state", {"target": "depolarizing", "p": "1/16", "rows": [[[0, 0], [0, 0]]]}),
    ("nan amplitude", {"target": "depolarizing", "p": "1/16", "rows": [[["nan", 0], [1, 0]]]}),
    ("n too large", {"target": "two_pauli", "p": "1/10", "rows": [[[1, 0]] * 82]}),
    ("bad payload rank", dict(good, rank=3)),
]
for name, c in bad:
    rep = verify.verify_candidate(c)
    check(f"INVALID: {name}", rep["verdict"] == "INVALID", rep.get("reason"))
check("report has pins field", "pins" in verify.verify_candidate(good, quick=True))

# ----------------------------------------------------------------------------- 2. brute force

rng = random.Random(20261025)
cases = [(n, r, fam) for n in (1, 2, 3) for r in (2, 3) for fam in verify.FAMILIES]
cases += [(4, 2, "depolarizing"), (4, 2, "two_pauli")] if not FAST else []
worst = 0.0
for i, (n, r, fam) in enumerate(cases):
    rows = rand_rows(rng, n, r, real=(i % 4 == 0))
    p = Fraction(rng.randint(30, 200), 1000)
    probs = verify.channel_probs(fam, p)
    vrows = [[(z.real, z.imag) for z in row] for row in rows]
    Z = sum(Fraction(a) ** 2 + Fraction(b) ** 2 for row in vrows for a, b in row)
    bf = brute_ic_bits(rows, probs)
    for path in ((0, 1) if fam == "depolarizing" else (1,)):
        lo, hi, _ = verify.evaluate(vrows, n, Z, probs, path, 600)
        inside = float(lo) - 1e-9 <= bf <= float(hi) + 1e-9
        worst = max(worst, float(hi - lo))
        check(f"brute force inside enclosure: {fam} n={n} r={r} p={p} path {path}",
              inside and float(hi - lo) < 1e-9, f"bf={bf} enc=[{lo}, {hi}]")
print(f"     (widest enclosure {worst:.2e} bits)")

# ----------------------------------------------------------------------------- 3. arXiv:2608.15870

def sign_case(path, p, want_verdict, want_sign, label):
    rep = verify.verify_candidate(payload_cand(path, p))
    lo, hi = Fraction(rep.get("I_c_lo", "0")), Fraction(rep.get("I_c_hi", "0"))
    sign = 1 if lo > 0 else (-1 if hi < 0 else 0)
    check(f"KG2026 {label} at {p}: {want_verdict}, sign {want_sign:+d}",
          rep["verdict"] == want_verdict and sign == want_sign and rep.get("recheck") == "path 1 agrees",
          rep["summary"])
    return rep


repA = sign_case(A, "16239/250000", "MATCHES_RECORD", +1, "A")
lo = float(repA.get("I_c_lo", 0))
check("KG2026 A enclosure matches the paper (8.4e-8 .. 8.7e-8 bits)", 8.4e-8 <= lo <= 8.7e-8, repA["summary"])
sign_case(A, "203/3125", "BELOW_RECORD", -1, "A")
if not FAST:
    sign_case(R, "4057/62500", "BELOW_RECORD", +1, "R")
    sign_case(R, "16229/250000", "BELOW_RECORD", -1, "R")

# ----------------------------------------------------------------------------- 4. verdict paths

with with_records(lambda d: d["cells"]["depolarizing"].update(bar="16239/250000", match_from="1/16")):
    rep = verify.verify_candidate(good)
    check("IMPROVEMENT when the bar is below a certified point (doctored records)",
          rep["verdict"] == "IMPROVEMENT", rep["summary"])
with with_records(lambda d: d["cells"]["depolarizing"].update(zero_capacity_from="1/20")):
    rep = verify.verify_candidate(good)
    check("VERIFIER_ALARM when I_c > 0 is certified where Q = 0 (doctored bound)",
          rep["verdict"] == "VERIFIER_ALARM", rep["summary"])

rows1 = as_json_rows(rand_rows(rng, 6, 1))
rep = verify.verify_candidate({"target": "depolarizing", "p": "1/16", "rows": rows1})
check("UNDECIDED for a pure input (I_c = 0 exactly)", rep["verdict"] == "UNDECIDED", rep["summary"])
rep = verify.verify_candidate(good, quick=True)
check("SCREEN_ONLY in quick mode", rep["verdict"] == "SCREEN_ONLY", rep["summary"])

orig_eval = verify.evaluate


def shifted(rows, n, Z, probs, path, timeout):
    lo, hi, dm = orig_eval(rows, n, Z, probs, path, timeout)
    return (lo + 1, hi + 1, dm) if path == 1 else (lo, hi, dm)


verify.evaluate = shifted
rep = verify.verify_candidate(good)
verify.evaluate = orig_eval
check("VERIFIER_ALARM when the two engine paths disagree", rep["verdict"] == "VERIFIER_ALARM", rep["summary"])

try:
    verify.entropy_bounds((0, 1e-12, [0.5, -1e-6]), "test")
    check("VERIFIER_ALARM on a certified negative eigenvalue", False)
except verify.Alarm:
    check("VERIFIER_ALARM on a certified negative eigenvalue", True)

orig_run = verify.run_engine


def bad_trace(inp, timeout):
    blocks = orig_run(inp, timeout)
    m, r = blocks[0]["trS"]
    blocks[0]["trS"] = (m + 1e-3, r)
    blocks[0]["trJ"] = (m + 1e-3, r)
    return blocks


verify.run_engine = bad_trace
rep = verify.verify_candidate({"target": "two_pauli", "p": "1/10", "rows": as_json_rows(rand_rows(rng, 5, 2))})
verify.run_engine = orig_run
check("VERIFIER_ALARM on an output trace that differs from the input trace",
      rep["verdict"] == "VERIFIER_ALARM", rep["summary"])

h1 = verify.verify_candidate(dict(good, note="x"), quick=True)["hash"]
h2 = verify.verify_candidate(dict(good, note="y"), quick=True)["hash"]
check("hash ignores notes", h1 == h2)

# ----------------------------------------------------------------------------- 5. arXiv:2605.09138

cache = ROOT / "data" / "cache" / "akllsss2026"
if (cache / "xz30.c128").exists():
    for key, target, p_ok, p_bar in (("xz30", "independent_xz", "0.118370", "0.118372"),
                                     ("tp30", "two_pauli", "0.118066", "0.118068")):
        rep = verify.verify_candidate(payload_cand(cache / f"{key}.c128", p_ok, target))
        check(f"AKLLSSS2026 {key} certified at {p_ok}: MATCHES_RECORD", rep["verdict"] == "MATCHES_RECORD",
              rep["summary"])
        rep = verify.verify_candidate(payload_cand(cache / f"{key}.c128", p_bar, target))
        check(f"AKLLSSS2026 {key} negative at bar {p_bar}", rep["verdict"] == "BELOW_RECORD"
              and Fraction(rep["I_c_hi"]) < 0, rep["summary"])
    check("AKLLSSS2026 dep45 is KG2026 R (byte-identical)", (cache / "dep45.c128").read_bytes() == R.read_bytes())
else:
    print("skip AKLLSSS2026 states (run data/fetch_public_states.py)")

print(f"\n{len(FAIL)} failure(s)")
sys.exit(1 if FAIL else 0)
