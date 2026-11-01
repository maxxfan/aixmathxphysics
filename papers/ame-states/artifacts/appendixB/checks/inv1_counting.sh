#!/usr/bin/env bash
# Table VI row l = 2, (2,1^6): hand derivation (pairs #165, auditor #169) made executable by extra/inv_t1.py (written for the
# release): counting chain over all 3^6 rank patterns, then the r = 1 and r = 2 code checks.
set -euo pipefail
"$PY" "$ROOT/extra/inv_t1.py" | tee /dev/stderr
