#!/usr/bin/env bash
# Table VI row l = 2, (2^2,1^4), PAIRS implementation: chk_t2f4.py (all subspaces by dimension, 432 with block weight >= 3,
# all with a radical of size 4) and enum_t2f4.py (second enumeration via hyperbolic-plane pairs: 0 symplectic Z_A at
# block weight >= 3; relaxed control at block weight >= 2: 42264).
set -euo pipefail
"$PY" "$ROOT/pairs/chk_t2f4.py" | tee /dev/stderr
"$PY" "$ROOT/pairs/enum_t2f4.py" 4 1000000000 3 2>&1 | tee /dev/stderr | grep -E '^(allowed|hyperbolic|candidate|ordered|TOTAL)' | sed -E 's/ \([0-9]+s\)$//'
"$PY" "$ROOT/pairs/enum_t2f4.py" 4 0 2 2>&1 | tee /dev/stderr | grep -E '^(allowed|hyperbolic|candidate|ordered)'
