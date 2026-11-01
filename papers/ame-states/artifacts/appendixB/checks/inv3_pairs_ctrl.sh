#!/usr/bin/env bash
# pairs controls for the involution sectors (board #161): (i) 45-class twists (0,0) and mixed (0,1) give no candidate K;
# (ii) THR=3 distance-4 control on the same K list: solutions are isotropic, rank 16, Omega-invariant, min weight 4
# (bad_cuts 8 / 2 of 35); (iii) ctrl_fpf6.py 1 40: qdfs finds AME(6,4) states (bad_cuts 0 of 10). stdlib + C, ~15 s.
set -euo pipefail
"$PY" "$ROOT/pairs/enum_t3f2.py" 0 0 | tee /dev/stderr | sed -E 's/ \([0-9]+s\)$//'
"$PY" "$ROOT/pairs/enum_t3f2.py" 0 1 | tee /dev/stderr | sed -E 's/ \([0-9]+s\)$//'
"$PY" "$ROOT/pairs/enum_t3f2.py" 1 1 1 6 3 | tee /dev/stderr | sed -E 's/ \([0-9]+s\)$//'
"$PY" "$ROOT/pairs/ctrl_fpf6.py" 1 40 | tee /dev/stderr | sed -E 's/ \([0-9.]+s\)$//'
