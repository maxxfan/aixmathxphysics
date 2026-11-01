#!/usr/bin/env bash
# pairs, full (board #161): enum_t3f2.py 1 1 1 runs the lift DFS on all 8640 K: nsol_total 0. stdlib + C, ~4-6 min.
set -euo pipefail
"$PY" "$ROOT/pairs/enum_t3f2.py" 1 1 1 | tee /dev/stderr | sed -E 's/ \([0-9]+s\)$//'
