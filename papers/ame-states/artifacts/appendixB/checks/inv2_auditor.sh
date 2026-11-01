#!/usr/bin/env bash
# Table VI row l = 2, (2^2,1^4): auditor/t2audit.py enumerates all 200787 four-dim subspaces of (F_2^2)^4: 432 have block
# weight >= 3, all with b-rank 2 (degenerate); control block weight >= 2: 42264 nondegenerate.
set -euo pipefail
"$PY" "$ROOT/auditor/t2audit.py" | tee /dev/stderr
