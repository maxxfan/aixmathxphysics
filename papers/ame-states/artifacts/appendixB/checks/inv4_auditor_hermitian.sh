#!/usr/bin/env bash
# Auditor t = 4, full: all 250 F_16-Hermitian MDS [4,2,3] codes (every AME(4,4) code that is F_16-linear up to coordinate
# order), DFS over all 2^36 symmetric beta: 0 solutions each (inv4h.py 250 6 == work/auditor/inv4h_6.log). numpy, ~2.5 min.
set -euo pipefail
"$PYNP" "$ROOT/auditor/inv4h.py" 250 6
sed -E 's/, "t": [0-9.]+//; s/"elapsed": [0-9.]+, //; s/, "elapsed": [0-9.]+//' inv4h_6.log
