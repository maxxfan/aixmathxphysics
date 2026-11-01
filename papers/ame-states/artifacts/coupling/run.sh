#!/usr/bin/env bash
# Replay the coupling certificates (paper Sec. V, Table IV). Usage: ./run.sh [quick|full]
#   quick (default): mapping check, all 28 certificates with the C enumerator + seeded Python sample, negative controls,
#                    field-level descriptions of the 27 ingredients (fieldlevel/run.sh, about 5 s)
#   full: additionally every certificate exhaustively in pure Python (about 70 s)
set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-quick}"
mkdir -p out
fail=0
check() { if diff -q "$1" "$2" >/dev/null; then echo "PASS  $3"; else echo "FAIL  $3 (diff $1 $2)"; fail=1; fi; }

python3 check_mapping.py > out/mapping.txt || true
check out/mapping.txt expected/mapping.txt "ingredient-to-cell mapping (28 certificates, 27 ingredients)"

python3 ../checker/ame_check.py coupling certificates/*.json > out/coupling_check_C.txt || true
check out/coupling_check_C.txt expected/coupling_check_C.txt "28 coupling certificates: both ingredients exactly checked (C + Python sample)"

: > out/negative_controls.txt
neg_ok=1
for f in negative/*.json; do
  if python3 ../checker/ame_check.py coupling "$f" >> out/negative_controls.txt 2>&1; then neg_ok=0; fi
done
check out/negative_controls.txt expected/negative_controls.txt "negative controls (each must FAIL)"
if [ "$neg_ok" = 1 ]; then echo "PASS  every negative control was rejected"; else echo "FAIL  a negative control was accepted"; fail=1; fi

if fieldlevel/run.sh > out/fieldlevel_run.txt 2>&1; then echo "PASS  field-level descriptions of the 27 ingredients (fieldlevel/run.sh: field checks, conversion == certified rows)"; else echo "FAIL  fieldlevel/run.sh (see out/fieldlevel_run.txt)"; fail=1; fi

if [ "$MODE" = "full" ]; then
  python3 ../checker/ame_check.py coupling --python-only certificates/*.json > out/coupling_check_python.txt || true
  check out/coupling_check_python.txt expected/coupling_check_python.txt "28 coupling certificates, pure-Python exhaustive"
fi
if [ "$fail" = 0 ]; then echo "COUPLING: PASS"; else echo "COUPLING: FAIL"; fi
exit $fail
