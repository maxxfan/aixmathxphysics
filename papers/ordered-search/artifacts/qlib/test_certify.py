#!/usr/bin/env python3
"""Tests for the spectral-factor parts of qlib/certify.py (needs python-flint).

  * power series: log(exp a) = a, exp(a + b) = exp(a) exp(b), a (1/a) = 1;
  * spike corrections: the tan-map and graded core rules agree, and refine is unchanged by
    factoring out exact_residual (it still reaches the tolerance);
  * cepstral_factor from h = 1 reaches the tolerance for q_1 of (4, 605) on a 2^16 grid and q_2
    of (5, 5000) on 2^18, where 1/|h|^2 has spikes narrower than the grid;
  * a (4, 605) certificate whose factors all come from cepstral_factor passes the verifier.

The exact chains are read from the committed certificates in results/certificates.
usage: .venv/bin/python qlib/test_certify.py   (exit status 1 on any failure)
"""
import gzip
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "verifier"))
import certify  # noqa: E402
import verify  # noqa: E402

failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


def load(name):
    p = ROOT / "results" / "certificates" / name
    raw = p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix == ".gz" else raw)


def setup(cert, t, Mexp):
    n, D = cert["n"], cert["den"]
    Ct = cert["q"][t - 1]
    c = np.array([v / D for v in Ct])
    qmin = float(certify.eval_grid(c, 1 << int(np.ceil(np.log2(64 * n)))).min())
    gamma = 0.5 * qmin
    tol = 0.05 * gamma / (2 * n)
    M = 1 << Mexp
    return Ct, D, c, gamma, tol, M, certify.find_spikes(c, gamma, M)


# 1. power series
rng = np.random.default_rng(0)
a = rng.standard_normal(300) * 0.3 / np.arange(1, 301)
b = rng.standard_normal(300) * 0.3 / np.arange(1, 301)
ea = certify._sexp(a, 300)
check(np.abs(certify._slog(ea / ea[0], 300)[1:] - a[1:]).max() < 1e-14 and abs(ea[0] - np.exp(a[0])) < 1e-14,
      "log(exp a) = a")
check(np.abs(certify._sexp(a + b, 300) - certify._smul(ea, certify._sexp(b, 300), 300)).max() < 1e-14,
      "exp(a + b) = exp(a) exp(b)")
check(np.abs(certify._smul(ea, certify._sinv(ea, 300), 300) - np.r_[1.0, np.zeros(299)]).max() < 1e-14,
      "a (1/a) = 1")

# 2. spike rules against a fine graded reference, at a factor that is not yet converged (the 2^17
#    cepstrum of q_1, residual 3.6e-8), and refine from that factor
k4 = load("cert_k4_n605.json")
H = 1 << 60
Ct, D, c, gamma, tol, M, sp = setup(k4, 1, 17)
n = len(Ct)
G = int(round(gamma * D))
h0 = certify.spectral_factor(c, gamma, tol=1e-30, M0=M, Mmax=M)[0]
hi0 = [int(round(x * H)) for x in h0]
r = certify.exact_residual(Ct, D, G, hi0, H)
hf = np.array(hi0, dtype=float) / H
Fg = certify.eval_grid(r, M) / np.abs(np.fft.rfft(hf, M)) ** 2
ref = certify.spike_coeffs(hf, r, n, M, Fg, *sp, core="graded", npan=20, npt=16)
for core in ("tan", "graded"):
    s = certify.spike_coeffs(hf, r, n, M, Fg, *sp, core=core)
    rel = np.abs(s - ref).max() / np.abs(ref).max()
    check(len(sp[0]) > 0 and rel <= 1e-6, f"(4, 605) q_1, {len(sp[0])} spikes at 2^17: {core} rule within {rel:.1e} "
          f"of the fine reference")
hr, er = certify.refine(Ct, D, gamma, h0, H, M, tol, iters=6, spikes=sp)
check(np.abs(certify.exact_residual(Ct, D, G, hi0, H)).max() > tol >= er,
      f"refine from the 2^17 cepstrum reaches the tolerance ({er:.1e} <= {tol:.1e})")

# 3. cepstral steps from h = 1; on the 2^16 grid refine from the cepstrum fails for (4, 605) q_1
Ct, D, c, gamma, tol, M, sp = setup(k4, 1, 16)
h16 = certify.spectral_factor(c, gamma, tol=1e-30, M0=M, Mmax=M)[0]
e16 = certify.refine(Ct, D, gamma, h16, H, M, tol, iters=6, spikes=sp)[1]
check(e16 > tol, f"(4, 605) q_1 at 2^16: cepstrum + refine stays at {e16:.1e} > {tol:.1e} (the case for cepstral_factor)")
for name, t, Mexp in (("cert_k4_n605.json", 1, 16), ("cert_k5_n5000.json.gz", 2, 18)):
    cert = k4 if name == "cert_k4_n605.json" else load(name)
    Ct, D, c, gamma, tol, M, sp = setup(cert, t, Mexp)
    hi, err = certify.cepstral_factor(Ct, D, gamma, H, M, tol, sp)
    check(err <= tol, f"({cert['k']}, {cert['n']}) q_{t}: cepstral_factor at 2^{Mexp} ({len(sp[0])} spikes, closest "
          f"{sp[1].min():.1e} vs grid spacing {2 * np.pi / M:.1e}) reaches {err:.1e} <= {tol:.1e}")

# 4. a certificate built from cepstral factors
cert = {key: k4[key] for key in ("type", "k", "n", "den", "q", "h_den")}
cert["h"] = []
for t in range(1, k4["k"]):
    Ct, D, c, gamma, tol, M, sp = setup(k4, t, 16)
    cert["h"].append(certify.cepstral_factor(Ct, D, gamma, H, M, tol, sp)[0])
rep = verify.verify_candidate(cert)
check(rep.get("margins") is not None and min(rep["margins"]) > 0,
      f"(4, 605) with cepstral factors: verifier accepts, min margin "
      f"{min(rep['margins']) if rep.get('margins') else None}")

print(f"{len(failures)} failures")
sys.exit(1 if failures else 0)
