#!/usr/bin/env bash
# Table VI, row l = 5 -- auditor's independent implementation (board #184, argument.md 5.6; shares no code with pairs).
# (1) c5audit.py 5: own F_2^32 model of Omega; enumerates all 2^24 vectors of V_nf = Im(Omega+1) and checks the model
#     (15 * 279825 isotropic vectors = isotropic points of H(5,16)); 4155 isotropic F_16-lines contain a word of weight <= 4.
#     Its quadratic line-triple search is NOT run to completion (it stops: "too many lines"); the exhaustive part is (2).
# (2) c5graph.py 0 0: the graph family in the auditor's coordinates, brute force over (u1,u2) in W^2 per w:
#     300 solutions for each of 1105 w, 331,500 Lagrangians, min-weight histogram {2: 4500, 3: 54000, 4: 273000}.
set -euo pipefail
rc=0; "$PYNP" "$ROOT/auditor/c5audit.py" 5 > c5audit.out 2>&1 || rc=$?
cat c5audit.out >&2
# exit code 1 is expected here (the script stops before its quadratic pass); anything else is a failure
if [ "$rc" -ne 0 ] && [ "$rc" -ne 1 ]; then exit "$rc"; fi
sed -E 's/^\[ *[0-9.]+s\] //' c5audit.out
echo "c5audit exit code: $rc"
"$PYNP" "$ROOT/auditor/c5graph.py" 0 0 > c5graph.out 2>&1 || { cat c5graph.out >&2; exit 1; }
cat c5graph.out >&2
grep -vE '^\[' c5graph.out
