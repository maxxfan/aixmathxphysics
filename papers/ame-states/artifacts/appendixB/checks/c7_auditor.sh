#!/usr/bin/env bash
# Table VI, row l = 7, sigma = (7,1) -- auditor only (single implementation; board #132(a), argument.md 5.4).
# c7sym.py: all F_8-subspaces L1 of F_8^4 (r = 0..4: 1 + 585 + 4745 + 585 + 1 = 5917); S_nf = L1 + L1^perp; none has
# min party weight >= 5, best minimum party weight 4.  40 sampled r = 2 lines rebuilt as 16-qubit stabilizers and checked
# with the verifier's rank test (ame.bad_cuts, 7..35 bad cuts of 35).  c7dbg.py: for 6 random lines the directly
# enumerated 2^16-word party-weight histogram equals the histogram predicted by the formulas used in c7sym.
set -euo pipefail
"$PYNP" "$ROOT/auditor/c7sym.py" > c7sym.out 2>&1 || { cat c7sym.out >&2; exit 1; }
cat c7sym.out >&2
"$PY" - c7sym.out <<'PYEOF'
import json, sys
tot = 0
for line in open(sys.argv[1]):
    line = line.strip()
    if not line.startswith("{"):
        continue
    for r, d in json.loads(line).items():
        tot += d["lines"]
        print("r=%s subspaces %d pass_weight5 %d pass_both %d best_min_party_weight %d"
              % (r, d["lines"], d["pass_weight5"], d["pass_both"], d["best_min_party_weight_of_S_over_all_lines"]))
        cc = d["crosscheck(ok_i,ok_ii,rank16,isotropic,bad_cuts)"]
        if cc:
            print("r=%s crosscheck: %d sampled lines, rank 16 and isotropic: %d, bad_cuts min %d max %d (of 35), formula-pass: %d"
                  % (r, len(cc), sum(1 for x in cc if x[2] and x[3]), min(x[4] for x in cc), max(x[4] for x in cc),
                     sum(1 for x in cc if x[0] and x[1])))
print("total subspaces", tot)
PYEOF
"$PYNP" "$ROOT/auditor/c7dbg.py" > c7dbg.out 2>&1 || { cat c7dbg.out >&2; exit 1; }
cat c7dbg.out >&2
grep -E '^trial' c7dbg.out
