#!/usr/bin/env bash
# Table VI row l = 5, (5,1^3): "Counting forces two twisted fixed parties" (pairs #175, auditor #184), made executable by
# extra/c5_counting.py (written for the release).
set -euo pipefail
"$PY" "$ROOT/extra/c5_counting.py" | tee /dev/stderr
