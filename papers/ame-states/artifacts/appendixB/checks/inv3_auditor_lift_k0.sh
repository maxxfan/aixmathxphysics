#!/usr/bin/env bash
# Auditor lift DFS (t3lift.py, numpy), quick subset: unreduced K 0, 1, 2 at threshold 5: 0 solutions, per-level survivor
# profiles == work/auditor/t3lift_5_0_3.log.  ~30 s.
set -euo pipefail
"$PYNP" "$ROOT/auditor/t3lift.py" 5 0 3 1e9 | tee /dev/stderr | sed -E 's/^\[ *[0-9.]+s\] //; s/, [0-9]+s$//'
