#!/usr/bin/env bash
# Replay the six explicit AME states (paper Theorems 1 and 2): AME(15,8), AME(17,8), AME(18,16), AME(19,16), AME(20,16), AME(21,16).
#   1. field level: printed coefficients -> Lemma 2 cut ranks (graphs) or unitarity + all column selections (Hermitian codes)
#   2. printed coefficients -> canonical F_2 rows with the documented converter (qlib_ame.py) == certified candidate files
#   3. canonical level: the standalone checker (../checker) on the candidate files, every balanced cut
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p out
fail=0
check() { if diff -q "$1" "$2" >/dev/null; then echo "PASS  $3"; else echo "FAIL  $3 (diff $1 $2)"; fail=1; fi; }
python3 check_states.py --convert > out/fieldlevel.txt || true
check out/fieldlevel.txt expected/fieldlevel.txt "six states: field-level checks and printed data == candidates (about 35 s)"
python3 ../checker/ame_check.py state candidates/*.json > out/canonical_check.txt || true
check out/canonical_check.txt expected/canonical_check.txt "six states: canonical F_2 check of every balanced cut (about 20 s)"
if [ "$fail" = 0 ]; then echo "STATES: PASS"; else echo "STATES: FAIL"; fi
exit $fail
