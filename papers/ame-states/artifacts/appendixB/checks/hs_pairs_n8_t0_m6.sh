#!/usr/bin/env bash
# Prop. prop:loc, uniform H (Hermitian, type 0) at n = 8, PAIRS implementation pairs/dfs4b.c, row-0 weight m = 6
# (m = 4..7 covers every F_4 graph state: a generator X_v Z_N(v) has weight 1 + deg v >= 5).  Expect nsol = 0.
set -euo pipefail
"$BUILD/dfs4b" 8 0 6 0 1 1000000000 | tee /dev/stderr | grep '^DONE' | sed -E 's/ time=[0-9.]+s//'
