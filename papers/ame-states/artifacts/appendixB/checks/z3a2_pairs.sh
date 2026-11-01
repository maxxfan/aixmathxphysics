#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h with 2-dim fixed space" (a = 2 sector), PAIRS implementation (board #238; written without
# reading the auditor's code): a2_audit.py stages 1-2: 270 X_A lines, one orbit under Sp_A x C(h6) x C(h7) (the auditor's
# 3 classes of 90 under Sp_A alone), X_B: 3 Sp_B-orbits of 90, exactly one compatible with the X_A representative, and the
# exhaustive binary-branching DFS (88,863 nodes) giving 384 completions T; a2_stage3.py: |allowed| = 2928 of 4096 for every T,
# 576 + 576 planes, 55,296 orthogonal pairs, 144 good D, 0 Lagrangians for all 384 T; positive control a2_stage3.py with
# threshold 4 on the first 40 T: 40/40 find L, S = L + T has 65536 words and min party weight 4.
set -euo pipefail
"$PYNP" -u "$ROOT/pairs/a2_audit.py" 1 > a2.out 2>&1 || { cat a2.out >&2; exit 1; }
cat a2.out >&2
echo "== a2_audit.py 1 (stages 1-2)"
sed -E 's/ t=[0-9]+s$//' a2.out
"$PYNP" -u "$ROOT/pairs/a2_stage3.py" 0 9999 > s3.out 2>&1 || { cat s3.out >&2; exit 1; }
cat s3.out >&2
echo "== a2_stage3.py 0 9999 (stage 3; printed every 20th T, run time removed)"
grep -E '^T#' s3.out | sed -E 's/^T#[0-9]+ //; s/ t=[0-9]+s$//' | sort | uniq -c | awk '{c=$1; $1=""; print c " x" $0}'
grep -E '^DONE' s3.out
"$PYNP" -u "$ROOT/pairs/a2_stage3.py" 0 40 4 > s3c.out 2>&1 || { cat s3c.out >&2; exit 1; }
cat s3c.out >&2
echo "== a2_stage3.py 0 40 4 (positive control, threshold 4 for L against T)"
grep -E '^   FOUND' s3c.out | sort | uniq -c | awk '{c=$1; $1=""; print c " x" $0}'
grep -E '^DONE' s3c.out
