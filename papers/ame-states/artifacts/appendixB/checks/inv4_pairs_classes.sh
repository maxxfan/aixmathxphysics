#!/usr/bin/env bash
# Table VI row l=2, sigma=(2^4), pairs implementation: all 514944 AME(4,4) graph forms, 37 orbits under GL(2,2)^4 x S_4,
# inv8 DFS over all 2^36 symmetric forms q for each representative: NSOL 0 (pairs classes_inv8.py + inv8.c; numpy, ~40 s, 0.5 GB).
set -euo pipefail
"$PYNP" "$ROOT/pairs/classes_inv8.py" > raw.txt
cat raw.txt >&2
sed -E 's/ \([0-9.]+s\)$//' raw.txt
