#!/usr/bin/env bash
# Prop. prop:loc, uniform H type (Hermitian self-dual [16,8]_4 with paired coordinates): auditor/hdfs.c, exhaustive DFS over
# M M^dagger = I with all party minors nonsingular.  h = 2, 3 (n = 4, 6: calibration), h = 4 (n = 8: N_norm = 0).
set -euo pipefail
for h in 2 3 4; do "$BUILD/hdfs" $h 2>&1; done | tee /dev/stderr | grep -vE '^rc=|time='
