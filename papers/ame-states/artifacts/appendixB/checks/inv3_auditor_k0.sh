#!/usr/bin/env bash
# Table VI row l=2, sigma=(2^3,1^2): the 12 ordered pairs of nondegenerate alternating forms (F_M, F_N) with F_M + F_N = J form
# ONE Sp(4,2)-orbit, so pairs' 12 x 720 = 8640 K are 720 up to gauge (auditor t3k0audit.py, stdlib, < 1 s).
set -euo pipefail
"$PY" "$ROOT/auditor/t3k0audit.py" | tee /dev/stderr
