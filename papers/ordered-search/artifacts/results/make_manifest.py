#!/usr/bin/env python3
"""Write MANIFEST.sha256: the SHA-256 of every released checker, test and data file and of every
certificate, with paths relative to the p09 folder.  Check with `shasum -a 256 -c MANIFEST.sha256`
(or `sha256sum -c`) from that folder.  Run again after any change to a listed file.

usage: python3 results/make_manifest.py
"""
import hashlib
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
FILES = ["verifier/verify.py", "verifier/PINS.json", "verifier/freeze.sh", "verifier/tests/test_verifier.py",
         "verifier/tests/fixtures/*.json", "results/check_independent.py", "results/test_check_independent.py",
         "results/audit.py", "results/audits.json", "results/make_manifest.py",
         "results/certificates/SHA256SUMS", "results/certificates/cert_*.json*", "data/*.json"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


lines = []
for pat in FILES:
    hits = sorted(ROOT.glob(pat))
    if not hits:
        raise SystemExit(f"no file matches {pat}")
    lines += [f"{sha256(p)}  {p.relative_to(ROOT).as_posix()}" for p in hits]
(ROOT / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
print(f"MANIFEST.sha256: {len(lines)} files")
