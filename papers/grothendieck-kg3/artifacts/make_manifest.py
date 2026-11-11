#!/usr/bin/env python3
"""Write MANIFEST.sha256 for this bundle: SHA-256 of every regular file except the manifest itself, build products
(verifier/accum, __pycache__, lock files) and this script.   usage: python3 -I make_manifest.py"""
import hashlib
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
SKIP = {"MANIFEST.sha256", "make_manifest.py", "build_all_certificates.py", "accum", ".heavy.lock", ".DS_Store"}
lines = []
for f in sorted(HERE.rglob("*")):
    if f.is_file() and f.name not in SKIP and "__pycache__" not in f.parts:
        lines.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(HERE).as_posix()}")
(HERE / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
print(f"{len(lines)} files in MANIFEST.sha256")
