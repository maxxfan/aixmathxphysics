#!/usr/bin/env bash
# Table VI, row l = 3, sigma = (3^2,1^2), h = 1 -- pairs' full pipeline as originally run (work/pairs/z3_scan.log):
# classes_z3.py (numpy): 514944 AME(4,4) graph forms, 414 orbits under GL(2,2)^4 (z3_reps.txt);
# scan_z3.py + scan_z3.c with those 414 representatives: 77760 phi, 0 pass S_nf alone, NSOL = 0.
set -euo pipefail
"$PYNP" "$ROOT/pairs/classes_z3.py" > classes.out 2>&1 || { cat classes.out >&2; exit 1; }
cat classes.out >&2
sed -E 's/ \([0-9.]+s\)$//' classes.out
"$PY" "$ROOT/pairs/scan_z3.py" > scan.out 2>&1 || { cat scan.out >&2; exit 1; }
cat scan.out >&2
sed -E 's/ \([0-9.]+s\)$//' scan.out
