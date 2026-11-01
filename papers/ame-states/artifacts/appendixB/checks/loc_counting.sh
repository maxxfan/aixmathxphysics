#!/usr/bin/env bash
# Prop. prop:loc, counting step (odd l and l = 2): the auditor's fixpt_count.py (argument.md 5.3 and addendum, board #118,
# #119) for n = 8 (the claim) and n = 4, 6 (sanity: AME(4,4), AME(6,4) exist, so the counting must not exclude them).
set -euo pipefail
for n in 8 4 6; do "$PY" "$ROOT/auditor/fixpt_count.py" $n; done | tee /dev/stderr
