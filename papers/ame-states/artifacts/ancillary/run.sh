#!/usr/bin/env bash
# Replay the ancillary claims of the paper. Usage: ./run.sh [quick|full]
#   circulants/  Sec. III and Table II: exhaustive scans (336 and 48 AME hits over GF(8), 28 and 4 orbits), first-moment columns
#   quadric/     Sec. IV: the shortened [17,8,10]_256 code of AME(18,16) is not GRS (19 vs 21 quadrics; control on a normal rational curve)
#   shadow/      Sec. IV: the shadow inequality (purity form) excludes no cell that the snapshot leaves uncovered
#   codetables/  Sec. III: the [[15,0,8]]_8 and [[17,0,9]]_8 codes shown on codetables.de (pages cached 2026-11-10 in
#                ../snapshot/sources/) are symmetric circulant graph codes from the scan lists; the second equals the
#                state of Theorem 1 up to local Cliffords and a party permutation, the first is not (an invariant
#                separates the four classes of n = 15 scan orbits)
#   ame910/      Sec. V: the two ingredients of the Lean-verified AME(9,10) state (an instance of Lemma 3 with K = 4)
# quick: about 4 minutes. full: also re-checks every hit in pure Python and runs the zero-hit scans of Table II
#   except (19,9) and (19,11) (about 10 and 60 minutes; commands in README.md, outputs in expected/).
set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-quick}"
mkdir -p circulants/out quadric/out shadow/out codetables/out ame910/out
fail=0
check() { if diff -q "$1" "$2" >/dev/null; then echo "PASS  $3"; else echo "FAIL  $3 (diff $1 $2)"; fail=1; fi; }

python3 circulants/lambda_table.py > circulants/out/lambda_table.txt
check circulants/out/lambda_table.txt circulants/expected/lambda_table.txt "Table II first-moment columns (exact, denominator q^n+1)"

python3 circulants/rescan.py hits "$MODE" > circulants/out/hits_$MODE.txt
check circulants/out/hits_$MODE.txt circulants/expected/hits_$MODE.txt "GF(8) symmetric circulants: 336 hits / 28 orbits (n=15), 48 hits / 4 orbits (n=17)"
python3 - <<'PY' || fail=1
import json, sys
h = json.load(open("circulants/out/hits.json")); e = json.load(open("circulants/expected/hits.json"))
ok = all(sorted(map(tuple, h[n]["orbit_representatives"])) == sorted(map(tuple, e[n]["orbit_representatives"])) and
         sorted(map(tuple, h[n]["hits"])) == sorted(map(tuple, e[n]["hits"])) for n in ("15", "17"))
print(("PASS" if ok else "FAIL") + "  hit lists and orbit representatives identical to expected/hits.json")
sys.exit(0 if ok else 1)
PY

if [ "$MODE" = "full" ]; then
  python3 circulants/rescan.py zero 16,8 16,7 16,5 19,5 13,4 15,4 19,7 19,8 > circulants/out/zero_full.txt
  check circulants/out/zero_full.txt circulants/expected/zero_full.txt "zero-hit scans of Table II (all rows except (19,9), (19,11))"
fi

( python3 quadric/grs_quadric_test.py 92,74,70,173,117,59,208,212,182; python3 quadric/grs_quadric_test.py 92,74,70,173,117,59,208,212,182 4 ) > quadric/out/quadric.txt
check quadric/out/quadric.txt quadric/expected/quadric.txt "AME(18,16): shortened code not GRS (rank 17 -> 19 quadrics vs 21; RNC control 21)"

python3 shadow/shadow_check.py > shadow/out/shadow.txt
check shadow/out/shadow.txt shadow/expected/shadow.txt "shadow inequality over the phase-3 snapshot: no uncovered cell excluded; known exclusions reproduced"
python3 shadow/shadow_check.py 8 4 22 16 23 16 4 2 8 3 > shadow/out/shadow_cells.txt
check shadow/out/shadow_cells.txt shadow/expected/shadow_cells.txt "shadow inequality: AME(8,4), AME(22,16), AME(23,16) pass; controls AME(4,2), AME(8,3) fail"

python3 codetables/compare_codetables.py ../snapshot/sources/codetables_2026-11-10 circulants/expected/hits.json codetables/out/comparison.json > codetables/out/comparison.txt || fail=1
check codetables/out/comparison.txt codetables/expected/comparison.txt "codetables.de [[15,0,8]]_8 and [[17,0,9]]_8 (cached 2026-11-10): in the scan lists; n=17 equivalent to Theorem 1; 4 and 1 classes"
python3 codetables/lc_invariant.py codetables/out/comparison.json circulants/expected/hits.json > codetables/out/lc_invariant.txt || fail=1
check codetables/out/lc_invariant.txt codetables/expected/lc_invariant.txt "n=15: the 4 classes have pairwise different invariants; the code-table code is not locally Clifford equivalent to Theorem 1"

python3 ame910/check_ame910.py > ame910/out/ame910.txt || fail=1
check ame910/out/ame910.txt ame910/expected/ame910.txt "AME(9,10) (Lean proof of 2026-08-28): both ingredients of the K = 4 coupling re-checked"

if [ "$fail" = 0 ]; then echo "ANCILLARY: PASS"; else echo "ANCILLARY: FAIL"; fi
exit $fail
