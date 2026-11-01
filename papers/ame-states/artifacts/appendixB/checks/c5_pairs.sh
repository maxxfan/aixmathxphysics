#!/usr/bin/env bash
# Table VI, row l = 5, sigma = (5,1^3) with order-5 twists -- pairs' implementation (board #175).
# c5twist.py: exhaustive scan of the graph family (1105 isotropic points w of H(3,16) x 300 isometries = 331,500
# Lagrangians S_nf); prints the histogram of the minimum party weight of S_nf alone.  Then c5twist_val.py: 8 random
# candidates rebuilt as 16-qubit stabilizers; direct 2^16-word minimum party weight must equal the formula.
set -euo pipefail
"$PYNP" "$ROOT/pairs/c5twist.py" > c5twist.out 2>&1 || { cat c5twist.out >&2; exit 1; }
cat c5twist.out >&2
grep -E '^(form checks ok|isotropic projective points|DONE|SOLUTION|VALIDATE)' c5twist.out | sed -E 's/ [0-9]+ s$//'
"$PYNP" "$ROOT/pairs/c5twist_val.py" 8 > val.out 2>&1 || { cat val.out >&2; exit 1; }
cat val.out >&2
# canonical summary (robust to the numpy random stream: only properties that every trial must have)
awk '/^trial/ {n++; if ($0 ~ /isotropic True rank 16 /) ok++; if ($0 ~ / MATCH /) m++;
       for (i = 1; i <= NF; i++) if ($i == "wt") w = $(i+1);
       if (w <= 4) le4++;
       b = $0; sub(/.*bad_cuts \(/, "", b); split(b, a, ","); if (a[1] + 0 > 0) bad++ }
     END { print "c5twist_val trials: " n;
           print "  isotropic and rank 16: " ok;
           print "  direct 2^16-word min party weight == formula: " m;
           print "  min party weight <= 4: " le4;
           print "  bad_cuts > 0 (not AME): " bad }' val.out
