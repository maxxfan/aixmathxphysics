#!/usr/bin/env bash
# Prop. prop:loc, uniform S type (F_4-linear, F_4-bilinear alternating form): auditor/hsdfs.c, exhaustive DFS over
# M J M^T = J with all party minors nonsingular.  h = 2, 3 (calibration), h = 4 (n = 8: N_norm = 0).
set -euo pipefail
for h in 2 3 4; do "$BUILD/hsdfs" $h 2>&1; done | tee /dev/stderr | grep -vE '^rc=|time='
