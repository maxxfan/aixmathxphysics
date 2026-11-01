#!/usr/bin/env bash
# Prop. prop:loc calibration of the S-type search at n = 4 by brute force (extra/symplectic44.py, written for the release):
# all 979200 matrices of Sp(4,4), of which N_norm * 60^3 = 2 * 216000 = 432000 have nonsingular party minors.
set -euo pipefail
"$PY" "$ROOT/extra/symplectic44.py" | tee /dev/stderr
