#!/usr/bin/env bash
# Validations of the pairs t = 4 parametrization and DFS (board #129): val_inv8.py 7 120 (table test == brute-force 2^16-word
# min weight at thresholds 2,3,4), val2_inv8.py 4 (Monte-Carlo node counts per level vs the exact inv8 counts printed below),
# dump_check.py (the 96 depth-7 survivors of class 0, completed both ways, all have min weight 4).  numpy, ~35 s.
set -euo pipefail
"$PYNP" "$ROOT/pairs/val_inv8.py" 7 120 2>&1 | tee /dev/stderr | grep -E '^\(a\)|^min weight histogram|MISMATCH'
"$PY" "$ROOT/pairs/val2_inv8.py" 4 2>&1 | tee /dev/stderr | grep -E '^level'
echo "exact inv8 nodes for the same C (THR=4):"
echo "4401 d802 8104 6208 8210 8820 b40 3680" | "$BUILD/inv8" 0 1000000000000 4 | tee /dev/stderr | grep '^NSOL'
"$PYNP" "$ROOT/pairs/dump_check.py" 2>&1 | tee /dev/stderr | grep -E '^level-7|^\(min weight|^best near miss'
