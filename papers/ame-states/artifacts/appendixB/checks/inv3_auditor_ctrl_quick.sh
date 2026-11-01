#!/usr/bin/env bash
# Positive control of t3lift.py: threshold 4 (distance-4 states) on K 0 finds solutions; deterministic cap of 1000 solutions
# (T3_MAXSOL); the first solution is rebuilt directly: rank 16, isotropic, Omega-invariant, min weight 4. numpy, ~2 s.
set -euo pipefail
T3_MAXSOL=1000 "$PYNP" "$ROOT/auditor/t3lift.py" 4 0 1 1e9 | tee /dev/stderr | sed -E 's/^\[ *[0-9.]+s\] //; s/, [0-9]+s$//'
