#!/usr/bin/env bash
# Table VI, row l = 5 -- auditor's c5graph.py for a second choice of the order-5 twists (h5, h6) = indices (5, 101) of the
# 144 order-5 elements of Sp(4,2) (argument.md 5.6: "same result for two different pairs (h5, h6)"; log c5graph_5_101.log).
set -euo pipefail
"$PYNP" "$ROOT/auditor/c5graph.py" 5 101 > c5graph.out 2>&1 || { cat c5graph.out >&2; exit 1; }
cat c5graph.out >&2
grep -vE '^\[' c5graph.out
