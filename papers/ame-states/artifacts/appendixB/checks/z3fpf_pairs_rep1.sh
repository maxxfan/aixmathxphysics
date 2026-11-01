#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h fixed-point-free" (sector A), PAIRS implementation, stages 2-3 for T_A orbit
# representative 1 of pairs/z3a_reps.json (board #186, work/pairs/z3a_results.log): all T_B compatible with T_A, T = T_A + T_B,
# the 27 isotropic lines of T^perp/T filtered to weight-good completions L, then every g in Sp(4,2) tested via the
# (c, d = gc) table: g_ok_max = 0 (no valid g). Input: the stage-1 representatives (regenerated and compared in z3fpf_pairs_stage1).
set -euo pipefail
cp "$ROOT/pairs/z3a_reps.json" .
"$PYNP" -u "$ROOT/pairs/z3a_full.py" 1 > full.out 2>&1 || { cat full.out >&2; exit 1; }
cat full.out >&2
grep -vE '^  T_B ' full.out
