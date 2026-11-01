#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h with 2-dim fixed space" (a = 2), AUDITOR within-author cross-checks
# (ame84_f2linear_closure.md, a = 2 row; board #168, #241):
#  (1) z3twist_check.py nonfpf: 12 random completions T re-verified directly (dim 10, min party weight 5, isotropic,
#      Omega-invariant, inside V_nf);
#  (2) the generic, unstaged lagsearch.c (no use of the dim-2 structure; 5,735,899 nodes per T) on T = 0-6, 384-386,
#      768-770 (all three X_A classes) finds 0 Lagrangians, agreeing with the staged lagsearch2.c of check z3a2_auditor.
set -euo pipefail
"$PYNP" "$ROOT/auditor/z3twist_check.py" nonfpf > check.out 2>&1 || { cat check.out >&2; exit 1; }
cat check.out >&2
echo "== z3twist_check.py nonfpf"
grep -E '^(number of T|T: dim|checked)' check.out
"$PYNP" "$ROOT/auditor/z3twist4.py" nonfpf > z3twist4.out 2>&1 || { cat z3twist4.out >&2; exit 1; }
cat z3twist4.out >&2
echo "== lagsearch (generic DFS) on T 0-6, 384-386, 768-770"
for r in "0 6" "384 386" "768 770"; do "$BUILD/lagsearch" z3twist_Vs_nonfpf.bin $r; done | tee /dev/stderr
