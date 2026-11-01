#!/usr/bin/env bash
# Table VI row "l = 3, (3^2,1^2), h fixed-point-free" (sector A), PAIRS implementation, stage 1 (board #186):
# z3a_stage1.py enumerates T_A = S_nf cap {B-part = 0}: 233280 solutions (per-M histogram 108 M with 720, 72 M with 2160)
# and their orbits under Sp(4,2)_A x GU(2,2)_6 x GU(2,2)_7 (order 233280): 3 orbits of sizes 38880 / 77760 / 116640.
# It writes the 3 orbit representatives to z3a_reps.json; this check also confirms that they equal the representatives
# shipped in pairs/z3a_reps.json (copy of the original work/pairs/z3a_reps.json), which the checks z3fpf_pairs_rep{0,1,2}
# and z3fpf_pairs_ctrl use as input.
set -euo pipefail
"$PYNP" "$ROOT/pairs/z3a_stage1.py" orbits > st1.out 2>&1 || { cat st1.out >&2; exit 1; }
cat st1.out >&2
grep -E '^(form checks|GL\(2,4\)|GU\(2,2\)|Sp\(4,2\)|X_A solutions|per-M counts|orbits under)' st1.out
"$PY" - "$ROOT/pairs/z3a_reps.json" z3a_reps.json <<'EOF'
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
print("orbit representatives identical to the shipped pairs/z3a_reps.json:", a == b)
EOF
