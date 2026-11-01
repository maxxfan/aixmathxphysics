#!/usr/bin/env bash
# Table VI, row l = 3, sigma = (3^2,1^2), h = 1 -- pairs' implementation (board #138), numpy-free variant.
# scan_z3.py (stdlib DFS over F_4-bases): 77760 anti-isometries phi; scan_z3.c: number of phi whose graph S_nf alone has
# minimum party weight >= 5 (= 0).  The AME(4,4) representatives used for the (unneeded) S^sigma stage come from
# classes_z3.py (numpy); here an EMPTY representative list is used (R = 0), which does not affect the S_nf-alone count.
# The run with the 414 representatives, as in the original log, is check z3u_pairs_reps.  z3_check.py: independent
# Python recount of the low-weight words of every graph.
set -euo pipefail
: > empty_reps.txt
Z3_REPS=empty_reps.txt "$PY" "$ROOT/pairs/scan_z3.py" > scan.out 2>&1 || { cat scan.out >&2; exit 1; }
cat scan.out >&2
grep -vE '^progress' scan.out | sed -E 's/ \([0-9.]+s\)$//'
"$PY" "$ROOT/pairs/z3_check.py" 2>&1 | tee /dev/stderr
