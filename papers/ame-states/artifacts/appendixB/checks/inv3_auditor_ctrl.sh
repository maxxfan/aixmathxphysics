#!/usr/bin/env bash
# Positive control of t3lift.py as in the paper ("at threshold 4 it finds 1.4x10^6"): the original run (t3lift_4_0_2.log)
# stopped at a 200 s WALL-CLOCK cap with 1,433,150 solutions; here the same DFS is stopped deterministically by the solution
# cap T3_MAXSOL = 1433150 (checked on leaf entry, like the original time check), which reproduces the original count and
# per-level profile. numpy, ~3-5 min, < 50 MB.
set -euo pipefail
T3_MAXSOL=1433150 "$PYNP" "$ROOT/auditor/t3lift.py" 4 0 1 1e9 | tee /dev/stderr | sed -E 's/^\[ *[0-9.]+s\] //; s/, [0-9]+s$//'
