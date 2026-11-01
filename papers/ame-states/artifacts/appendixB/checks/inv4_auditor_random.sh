#!/usr/bin/env bash
# Auditor t = 4, full: 146 random AME(4,4) graph-state codes (seed 21, == work/auditor/inv4c_21.log), all 2^36 beta: 0 solutions.
# numpy, ~8 min.
set -euo pipefail
"$PYNP" "$ROOT/auditor/inv4c.py" 21 146 1e9
sed -E 's/, "t": [0-9.]+//; s/"elapsed": [0-9.]+, //; s/, "elapsed": [0-9.]+//' inv4c_21.log
