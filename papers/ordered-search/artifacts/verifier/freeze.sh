#!/usr/bin/env bash
# Run the verifier tests, pin source hashes, and make the verifier read-only.
# Run by a human after any reviewed change to the verifier.
set -euo pipefail
cd "$(dirname "$0")/.."
chmod u+w verifier verifier/*.py 2>/dev/null || true
chmod u+w verifier/PINS.json 2>/dev/null || true
python3 verifier/tests/test_verifier.py
python3 - <<'PY'
import json, sys
sys.path.insert(0, "verifier")
import verify
json.dump(verify.source_hashes(), open("verifier/PINS.json", "w"), indent=1)
print("pinned:", verify.source_hashes())
PY
chmod a-w verifier/*.py verifier/PINS.json
echo "verifier frozen"
