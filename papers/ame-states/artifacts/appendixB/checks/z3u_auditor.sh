#!/usr/bin/env bash
# Table VI, row l = 3, sigma = (3^2,1^2), h = 1 -- auditor's independent code (board #146, argument.md 5.5(b)).
# z3_audit.py (stdlib): DFS over images of an F_4-basis finds the 77760 = |GU(4,2)| sigma-equivariant isometries phi;
# every graph {(a, phi a)} has minimum party weight exactly 4, with 9 / 18 / 45 words of weight <= 4 for 32400 / 43200 / 2160 phi.
set -euo pipefail
"$PY" "$ROOT/auditor/z3_audit.py" 2>&1 | tee /dev/stderr
