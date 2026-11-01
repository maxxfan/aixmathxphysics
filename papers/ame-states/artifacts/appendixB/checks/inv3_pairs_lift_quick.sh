#!/usr/bin/env bash
# Table VI row l=2, sigma=(2^3,1^2), pairs implementation (enum_t3f2.py + invtw.py + qdfs.c, stdlib + C), quick subset:
# all 12 K_0 and 8640 K are enumerated (720 per K_0); the lift DFS is run on the first 720 K: NSOL 0.  ~30 s.
set -euo pipefail
"$PY" "$ROOT/pairs/enum_t3f2.py" 1 1 1 720 | tee /dev/stderr | sed -E 's/ \([0-9]+s\)$//'
