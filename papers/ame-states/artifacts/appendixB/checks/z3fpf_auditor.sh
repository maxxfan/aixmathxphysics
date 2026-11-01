#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h fixed-point-free" (sector A), AUDITOR implementation (argument.md 5.8, board #202).
# z3twist6.py imports z3twist5.py (which imports z3twist.py fpf), so one run prints:
#   step 1 (z3twist.py): 233280 X_A solutions, 180 d-keys, 391,910,400 (X_A, X_B) pairs;
#   step 2 (z3twist5.py): 324 Sp(4,2)_A-orbits, 756 compatible orbit pairs, 7 classes under C(h6) x C(h7);
#   step 3 (z3twist5.py): completions T per class 24,24,24,18,16,16,16 (138 T), (T,f) hits 0 for all 720 f in Sp(4,2);
#   controls (z3twist6.py): (i) every T and S = S^Omega + T re-verified directly (2^16 words), (ii) 21 unreduced pairs
#   reproduce their class's T-count with 0 hits, (iii) best near-miss: 14 of 15 a pass, ame.bad_cuts 3 of 35.
set -euo pipefail
"$PYNP" "$ROOT/auditor/z3twist6.py" > z6.out 2>&1 || { cat z6.out >&2; exit 1; }
cat z6.out >&2
sed -E -e 's/ +time [0-9.]+s$//' \
       -e 's/^\[ *[0-9]+s\] //' \
       -e 's/, [0-9]+s$//' \
       -e 's/(computed): [0-9]+s$/\1/' z6.out
