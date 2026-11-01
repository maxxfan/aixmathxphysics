#!/usr/bin/env python3
"""Audit p09 certificates: frozen verifier, independent FLINT checker, float spot check.

For each certificate: (1) verifier/verify.py verdict and pins; (2) results/check_independent.py
(python-flint, different route for every condition); (3) a float evaluation of q_1..q_{k-1} on a
16x oversampled grid (sanity only).  Appends one record per certificate to results/audits.json.

usage: .venv/bin/python results/audit.py CERT.json.gz [...]
"""
import gzip
import hashlib
import json
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "verifier"))
sys.path.insert(0, str(HERE))
import check_independent  # noqa: E402
import verify  # noqa: E402


def float_min(cert):
    k, n, D = cert["k"], cert["n"], cert["den"]
    out = []
    for q in cert["q"]:
        c = np.array([v / D for v in q])
        M = 1 << int(np.ceil(np.log2(32 * n)))
        x = np.zeros(M)
        x[:n] = c
        x[1:n] *= 2
        out.append(float(np.real(np.fft.fft(x)).min()))
    return out


def main():
    audits = json.loads((HERE / "audits.json").read_text())
    for path in sys.argv[1:]:
        p = pathlib.Path(path)
        raw = p.read_bytes()
        cert = json.loads(gzip.decompress(raw) if p.suffix == ".gz" else raw)
        t0 = time.time()
        rep = verify.verify_candidate(cert)
        t1 = time.time()
        ind = check_independent.check(cert)
        t2 = time.time()
        fm = float_min(cert)
        rec = {"file": str(p.relative_to(ROOT)) if p.is_absolute() else str(p),
               "sha256": hashlib.sha256(raw).hexdigest()[:16],
               "k": cert["k"], "n": cert["n"],
               "verifier": {"version": rep.get("verifier"), "pins": rep["pins"], "verdict": rep["verdict"],
                            "beats_prior_explicit": rep.get("beats_prior_explicit"),
                            "margins": rep.get("margins"), "seconds": round(t1 - t0, 1)},
               "independent": {"constant": ind["constant"], "beats_7265": ind["beats_prior_explicit"],
                               "margins": ind["margins"], "seconds": round(t2 - t1, 1)},
               "float_grid_min": fm,
               "agree": rep["verdict"] not in ("INVALID", "VERIFIER_ALARM") and
                        all(abs(a - b) <= 1e-12 * max(1.0, abs(a)) for a, b in zip(rep["margins"], ind["margins"])),
               "date": time.strftime("%Y-%m-%d %H:%M")}
        audits["audits"][f"k{cert['k']}_n{cert['n']}"] = rec
        audits["events"].append(f"{rec['date']} audited {rec['file']}: {rep['summary']}")
        print(json.dumps(rec, indent=1))
    (HERE / "audits.json").write_text(json.dumps(audits, indent=1) + "\n")


if __name__ == "__main__":
    main()
