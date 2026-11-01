#!/usr/bin/env bash
# Prop. prop:loc calibration "34 992 of the 77 760 unitary 4x4 matrices over F_4": auditor/hform.py 4 (pruned DFS, original)
# and extra/unitary44.py (all unitary matrices first, then filter; written for the release).
set -euo pipefail
"$PY" "$ROOT/auditor/hform.py" 4 2>&1 | tee /dev/stderr | grep -E '^(n=|re-check)' | sed -E 's/ \(time [0-9]+s, /, (/; s/, limit [0-9.]+s\)/)/'
"$PY" "$ROOT/extra/unitary44.py" | tee /dev/stderr
