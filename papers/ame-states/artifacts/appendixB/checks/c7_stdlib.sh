#!/usr/bin/env bash
# Table VI, row l = 7, sigma = (7,1): numpy-free replay.  extra/c7sym_stdlib.py is a stdlib PORT of the auditor's c7sym.py
# written for this release (same enumeration of all 5917 F_8-subspaces, same printed summary, same 40 sampled r = 2
# cross-checks with the verifier's rank test).  It is NOT a second independent implementation of the original run.
set -euo pipefail
"$PY" "$ROOT/extra/c7sym_stdlib.py" > c7std.out 2>&1 || { cat c7std.out >&2; exit 1; }
cat c7std.out >&2
"$PY" - c7std.out <<'PYEOF'
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
            print("r=%s crosscheck bad_cuts sequence: %s" % (r, " ".join(str(x[4]) for x in cc)))
print("total subspaces", tot)
PYEOF
