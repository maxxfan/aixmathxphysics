"""Exact certificate from a saved chain (npz with u0..u_{k-1}), then verifier + independent check.

usage: python make_cert.py out/basis_k6_n45000.npz [--out ../../candidates/NAME.json.gz]
"""
import argparse
import gzip
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "qlib"))
sys.path.insert(0, str(ROOT / "verifier"))
import certify  # noqa: E402
import verify  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("npz")
ap.add_argument("--out")
ap.add_argument("--gamma", type=float, default=0.5)
a = ap.parse_args()
z = np.load(a.npz)
k, n = int(z["k"]), int(z["n"])
u = [z[f"u{i}"] for i in range(k)]
t0 = time.time()
cs = certify.chain_coeffs(u, n)
name = pathlib.Path(a.npz).name
solver = {"big": "smooth_big.py", "smooth": "smooth.py", "basis": "lp_basis.py"}.get(name.split("_")[0], "?")
meta = {"source": f"work/lead/{solver} chain {name} (beta {float(z['beta']):.4g}, "
                  f"true relative min {float(z['worst']):.4g}); qlib/certify.py",
        "made": time.strftime("%Y-%m-%d %H:%M")}
cert, diag = certify.make_certificate(cs, n, gamma_frac=a.gamma, meta=meta, verbose=True)
for d in diag:
    print(d, flush=True)
out = pathlib.Path(a.out or ROOT / "candidates" / f"cert_k{k}_n{n}.json.gz")
out.write_bytes(gzip.compress(json.dumps(cert, separators=(",", ":")).encode()))
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB) in {time.time() - t0:.0f}s", flush=True)
t1 = time.time()
rep = verify.verify_candidate(cert)
print(f"verifier ({time.time() - t1:.1f}s): {rep['summary']}  margins={rep.get('margins')}", flush=True)
r = subprocess.run([sys.executable, str(ROOT / "results" / "check_independent.py"), str(out)],
                   capture_output=True, text=True)
print("independent:", r.stdout.strip(), r.stderr.strip()[-300:])
