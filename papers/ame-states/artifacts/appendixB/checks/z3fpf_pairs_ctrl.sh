#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h fixed-point-free" (sector A), PAIRS implementation, controls (board #186):
#  (a) positive control: the same pipeline with weight threshold 4 (z3a_full_thr.py, Z3A_THR=4) on representative 1,
#      first T_B, finds valid g for its weight-good lines ("hundreds of valid g per line");
#  (b) validation z3a_val.py 1 4 3: for 4 random T_B, T^perp/T has exactly 27 isotropic lines (no weight filter), and the
#      16-qubit rows built for a random line and g are isotropic, of rank 16, with directly enumerated (2^16 words) minimum
#      party weight equal to the formula used by the search.  Only RNG-independent facts are compared; the sampled weights
#      and ame.bad_cuts values (original: weights 3,3,4,3, bad_cuts 8,26,6,29 of 35) are in the raw log.
set -euo pipefail
cp "$ROOT/pairs/z3a_reps.json" .
Z3A_THR=4 "$PYNP" -u "$ROOT/pairs/z3a_full_thr.py" 1 1 > thr.out 2>&1 || { cat thr.out >&2; exit 1; }
cut -c1-400 thr.out >&2
echo "== Z3A_THR=4 z3a_full_thr.py 1 1 (positive control, representative 1, first T_B)"
grep -E '^(rep |first T_B)' thr.out
"$PY" - thr.out <<'EOF'
import re, sys
found = [l for l in open(sys.argv[1]) if l.startswith("FOUND")]
counts = []
for l in found:
    m = re.search(r"\]\], \[([0-9, ]*)\]\)\s*$", l)
    counts.append(len(m.group(1).split(",")))
for k, c in enumerate(counts[:11]):
    print("FOUND line %d: number of valid g = %d" % (k, c))
print("positive control finds valid g:", len(counts) > 0 and min(counts) > 0)
EOF
"$PYNP" -u "$ROOT/pairs/z3a_val.py" 1 4 3 > val.out 2>&1 || { cat val.out >&2; exit 1; }
cat val.out >&2
echo "== z3a_val.py 1 4 3 (validation)"
grep -E '^trial ' val.out
grep -E '^   isotropic ' val.out | sed -E 's/direct min party weight [0-9]+ formula [0-9]+ //; s/ bad_cuts .*$//'
