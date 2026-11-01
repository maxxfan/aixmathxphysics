#!/usr/bin/env bash
# Build the PSD residual checker, run the tests, pin hashes (sources, binary, data), make read-only.
# Run by a human after any reviewed change, or after cloning (the binary hash is machine-specific).
set -euo pipefail
cd "$(dirname "$0")/.."
chmod u+w verifier verifier/*.py verifier/*.c verifier/PINS.json data data/*.json 2>/dev/null || true
rm -f verifier/rescheck
cc -O2 -ffp-contract=off -fno-fast-math -o verifier/rescheck verifier/rescheck.c -lm
python3 data/make_records.py --check
python3 verifier/tests/test_verifier.py --trials 40
python3 - <<'PY'
import json, sys
sys.path.insert(0, "verifier")
import verify
json.dump(verify.source_hashes(), open("verifier/PINS.json", "w"), indent=1)
print("pinned:", verify.source_hashes())
PY
chmod a-w verifier/*.py verifier/*.c verifier/rescheck verifier/PINS.json data/*.json data/sources/* data/make_records.py
echo "verifier frozen"
