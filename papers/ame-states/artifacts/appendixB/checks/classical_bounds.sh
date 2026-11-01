#!/usr/bin/env bash
# Classical facts quoted in Prop. prop:loc and Table VI row (2,1^6), made executable (extra/classical_codes.py, written for
# the release): no [8,4,5]_4 linear MDS code, A(8,5) = 4, A(6,5) = 2, no additive OA(256,8,4,4) / OA(16,6,4,2) of index 1,
# with positive controls (hexacode, (5,2)).
set -euo pipefail
"$PY" "$ROOT/extra/classical_codes.py" | tee /dev/stderr
