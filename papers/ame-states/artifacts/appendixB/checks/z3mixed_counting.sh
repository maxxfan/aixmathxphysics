#!/usr/bin/env bash
# Table VI row l = 3, (3^2,1^2), mixed twist types: hand derivation (auditor #168, pairs #191(3)) made executable by
# extra/z3_mixed.py (written for the release): parity of dim V_nf and forced 4-party words, all 9 (r6, r7).
set -euo pipefail
"$PY" "$ROOT/extra/z3_mixed.py" | tee /dev/stderr
