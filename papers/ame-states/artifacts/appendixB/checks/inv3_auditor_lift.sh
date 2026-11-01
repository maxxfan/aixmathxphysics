#!/usr/bin/env bash
# Auditor lift DFS, full: the 720 isometries reduce to 10 orbit representatives under C(l0)^2 x swap(6,7); threshold 5:
# 0 solutions for all 10, 1,133,400 surviving nodes (== work/auditor/t3lift_5_0_720.log). numpy, ~3-4 min.
set -euo pipefail
"$PYNP" "$ROOT/auditor/t3lift.py" 5 0 720 1e9 orb | tee /dev/stderr | sed -E 's/^\[ *[0-9.]+s\] //; s/, [0-9]+s$//'
