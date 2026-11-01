#!/usr/bin/env bash
# Positive control for the t = 4 DFS (ADDED FOR THE RELEASE, pairs code only): inv8 with its relaxed threshold THR=3 on the
# class-0 representative finds solutions (capped at 1000); each printed solution is rebuilt and has min party weight 4.
set -euo pipefail
"$PY" "$ROOT/extra/inv8_thr3_ctrl.py" 1000 | tee /dev/stderr
