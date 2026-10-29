#!/usr/bin/env bash
# Build the ball-arithmetic engine, run the tests, pin hashes (sources, binary, records), make
# read-only.  Run by a human after any reviewed change, or after cloning (the binary hash is
# machine-specific).  Optional: python3 data/fetch_public_states.py first, so the tests also cover
# the arXiv:2605.09138 states.
set -euo pipefail
cd "$(dirname "$0")/.."
chmod u+w verifier verifier/*.py verifier/*.c verifier/PINS.json data data/*.json 2>/dev/null || true
rm -f verifier/icball
cc -O2 -ffp-contract=off -fno-fast-math -o verifier/icball verifier/icball.c -lm
python3 verifier/tests/test_verifier.py
python3 - <<'PY'
import json, sys
sys.path.insert(0, "verifier")
import verify
json.dump(verify.source_hashes(), open("verifier/PINS.json", "w"), indent=1)
print("pinned:", verify.source_hashes())
PY
chmod a-w verifier/*.py verifier/*.c verifier/icball verifier/PINS.json data/*.json data/sources/*/*
echo "verifier frozen"
