#!/usr/bin/env bash
# run.sh -- replay of the F_2-linear stabilizer AME(8,4) frame searches (see README.md).
#
# Compiles both implementations, re-derives the V_0 class representatives with each implementation's own code, runs the full
# search over all ten classes with BOTH implementations, runs the positive controls, normalises every result to a canonical text
# (out/canon/<check>.txt), compares it with expected/<check>.txt and with the SHA-256 in expected/HASHES.txt, prints PASS/FAIL per
# check and an overall verdict.  Exit status 0 iff every check passes.  No network access; needs bash, a C compiler (cc) and
# python3 >= 3.9 (standard library only).
#
# Environment: JOBS (parallel processes, default 4), PYTHON (default python3), CC (default cc).
# Work: about 24 CPU-minutes in total (Apple M4), so the wall time is roughly 24/JOBS minutes on an otherwise idle machine plus
# about 1 minute of sequential steps.  Measured: 14 minutes with JOBS=4 on a 10-core Apple M4 that was heavily loaded by other jobs.
set -euo pipefail
cd "$(dirname "$0")"
HERE="$(pwd)"
PY="${PYTHON:-python3}"
CC="${CC:-cc}"
JOBS="${JOBS:-4}"
NICE="nice -n 10"
export PYTHONDONTWRITEBYTECODE=1   # keep the folder clean (no __pycache__)

# ----------------------------------------------------------------------------------------------------------- one parallel task
task() {
  local name="$1" a="$2"
  case "$name" in
    i1_frames3_t5)  # impl1 full search, one class (frames3.c, all conditions, t = 5)
      "$PY" impl1_classify/run.py "out/logs/i1_frames3_t5_c$a.log" build/frames3 4 5 out/work/v0_m4.txt "$a" -1 "out/work/sol_frames3_t5_c$a.txt" ;;
    i1_frames_t5)   # the same with the first version frames.c
      "$PY" impl1_classify/run.py "out/logs/i1_frames_t5_c$a.log" build/frames 4 5 out/work/v0_m4.txt "$a" -1 "out/work/sol_frames_t5_c$a.txt" ;;
    i1_m124)        # impl1 level-2 counts without the four J = {0,1,2} conditions
      $NICE "$PY" impl1_classify/run_m124_all.py build/frames3 out/work/v0_m4.txt "out/logs/i1_frames3_m124_c$a.log" "$a" ;;
    i2_full)        # impl2 full search, one full-gauge representative
      $NICE "$PY" impl2_auditor/run_v8.py build/vsearch8 out/work/v0reps8.txt "$a" "$a" 1200 0 "out/logs/i2_vsearch8_full_r$a.log" ;;
    i2_sp)          # impl2 full search, one Sp(4,2)-only representative
      $NICE "$PY" impl2_auditor/run_v8.py build/vsearch8 out/work/v0reps8_sp.txt "$a" "$a" 1200 0 "out/logs/i2_vsearch8_sp_r$a.log" ;;
    i2_gauge)       # impl2 gauge-invariance test, one representative, 3 random gauge transforms (seed 1000 + rep)
      $NICE "$PY" impl2_auditor/gauge_test8.py build/vsearch8 out/work/v0reps8.txt 3 "$((1000 + a))" "out/logs/i2_gauge_r$a.log" "$a" "$a" > /dev/null ;;
    i2_relaxed)     # impl2 positive control: relaxed thresholds (min party weight 4), rep 0, first 400 V_1
      $NICE "$PY" impl2_auditor/run_v8c.py build/vsearch8 out/work/v0reps8.txt 0 400 1800 out/logs/i2_relaxed.log out/work/v8c_relaxed_solutions.txt > /dev/null ;;
    i2_classcount)  # impl2 independent orbit-size recount, h = 2, 3, 4
      $NICE "$PY" impl2_auditor/classcount.py 4 > out/logs/i2_classcount.log ;;
    i2_enumtest)    # impl2 exactness test of the subspace enumerator
      $NICE build/enumtest > out/logs/i2_enumtest.log ;;
    *) echo "unknown task $name" >&2; return 2 ;;
  esac
}

if [[ "${1:-}" == "__task" ]]; then
  shift
  task "$@"
  exit $?
fi

T0=$(date +%s)
stamp() { echo "[$(( $(date +%s) - T0 )) s] $*"; }

"$PY" -c 'import sys; assert sys.version_info >= (3, 9), "python >= 3.9 required"; print("python", sys.version.split()[0])'
rm -rf build out
mkdir -p build out/logs out/work out/canon

# ------------------------------------------------------------------------------------------------------------------ compile
stamp "compiling (build/)"
$CC -O2 -o build/frames   impl1_classify/frames.c
$CC -O2 -o build/frames3  impl1_classify/frames3.c
$CC -O2 -o build/vsearch8 impl2_auditor/vsearch8.c
$CC -O2 -o build/vcount   impl2_auditor/vcount.c
$CC -O2 -I impl2_auditor -o build/enumtest impl2_auditor/enumtest.c 2> out/logs/enumtest_compile_warnings.log

# ------------------------------------------------------------------------------------- V_0 classes, each with its own code
stamp "impl1: V_0 classes (v0reps.py 2, 3, 4)"
for m in 2 3 4; do
  $NICE "$PY" impl1_classify/v0reps.py "$m" "out/work/v0_m$m.txt" > "out/logs/i1_v0reps_m$m.log"
done
stamp "impl2: V_0 classes (vclasses8.py)"
$NICE "$PY" impl2_auditor/vclasses8.py out/work > out/logs/i2_vclasses8.log

# ------------------------------------------------------------------------------------------------- fast positive controls
stamp "impl1 controls: n = 4, n = 6, n = 8 t = 4, relaxed n = 8"
"$PY" impl1_classify/run.py out/logs/i1_ctrl.log build/frames 2 3 out/work/v0_m2.txt -1 -1 out/work/sol_m2.txt
"$PY" impl1_classify/run.py out/logs/i1_ctrl.log build/frames 3 4 out/work/v0_m3.txt -1 -1 out/work/sol_m3.txt
"$PY" impl1_classify/run.py out/logs/i1_ctrl.log build/frames 4 4 out/work/v0_m4.txt 9 20 out/work/sol_m4_t4.txt
"$PY" impl1_classify/run.py out/logs/i1_ctrl_relaxed.log JMASKS=64479 PATTERN=1 build/frames3 4 5 out/work/v0_m4.txt 9 60 out/work/sol_relax_m4.txt
stamp "impl2 control: n = 4 and n = 6 counts (vcount.py)"
$NICE "$PY" impl2_auditor/vcount.py build/vcount > out/logs/i2_vcount.log

# --------------------------------------------------------------------------------------------- the searches (in parallel)
TASKS=(i2_relaxed - i2_classcount - i2_enumtest -)
for c in 0 1 2 3 4 5 6 7 8 9; do TASKS+=(i1_frames3_t5 "$c"); done
for c in 0 1 2 3 4 5 6 7 8 9; do TASKS+=(i1_frames_t5 "$c"); done
for c in 0 1 2 3 4 5 6 7 8 9; do TASKS+=(i1_m124 "$c"); done
for k in 0 1 2 3 4 5 6 7 8 9; do TASKS+=(i2_gauge "$k"); done
for k in 0 1 2 3 4 5 6 7 8 9; do TASKS+=(i2_full "$k"); done
for k in $(seq 0 71); do TASKS+=(i2_sp "$k"); done
stamp "running $(( ${#TASKS[@]} / 2 )) search tasks with JOBS=$JOBS (about 23 CPU-minutes in total)"
TASKFAIL=0
printf '%s %s\n' "${TASKS[@]}" | xargs -P "$JOBS" -n 2 bash "$HERE/run.sh" __task || TASKFAIL=1
stamp "searches finished (task failures: $TASKFAIL)"

# ----------------------------------------------------------------------------------- brute-force checks of found objects
stamp "brute-force checks of the frames / solutions found by the controls"
"$PY" impl2_auditor/vbrute8.py out/work/v8c_relaxed_solutions.txt 50 > out/logs/i2_vbrute_relaxed.log || true
"$PY" impl2_auditor/vbrute8.py out/logs/i2_vsearch8_full_r0.log.SOLlines_0.txt 30 > out/logs/i2_vbrute_nearmiss.log || true

# ------------------------------------------------------------------------------------------- canonicalise and compare
stamp "canonicalising (out/canon/) and comparing with expected/"
"$PY" tools/checks.py replay out out/canon
"$PY" tools/checks.py expected out/expected_regen

NPASS=0; NFAIL=0
: > out/SUMMARY.txt
report() { echo "$1" | tee -a out/SUMMARY.txt; }
if "$PY" tools/hashes.py verify expected/HASHES.txt expected > out/logs/expected_integrity.log; then
  report "PASS  expected_integrity        expected/*.txt match expected/HASHES.txt"; NPASS=$((NPASS + 1))
else
  report "FAIL  expected_integrity        expected/*.txt do not match expected/HASHES.txt (see out/logs/expected_integrity.log)"; NFAIL=$((NFAIL + 1))
fi
REGEN_OK=1
[[ "$(ls expected | grep -v '^HASHES.txt$')" == "$(ls out/expected_regen)" ]] || REGEN_OK=0
for f in expected/*.txt; do
  n="$(basename "$f")"
  [[ "$n" == "HASHES.txt" ]] && continue
  cmp -s "$f" "out/expected_regen/$n" || REGEN_OK=0
done
if [[ $REGEN_OK == 1 ]]; then
  report "PASS  expected_provenance       expected/ is regenerated identically from original_logs/ by tools/checks.py"; NPASS=$((NPASS + 1))
else
  report "FAIL  expected_provenance       expected/ differs from what tools/checks.py derives from original_logs/"; NFAIL=$((NFAIL + 1))
fi
"$PY" tools/hashes.py verify expected/HASHES.txt out/canon > out/logs/replay_hashes.log || true
while read -r h n; do
  [[ -z "$h" || "$h" == \#* ]] && continue
  name="${n%.txt}"
  if cmp -s "expected/$n" "out/canon/$n" && grep -q "^OK   $n " out/logs/replay_hashes.log; then
    report "$(printf 'PASS  %-27s sha256 %s' "$name" "${h:0:16}")"; NPASS=$((NPASS + 1))
  else
    report "$(printf 'FAIL  %-27s (diff expected/%s out/canon/%s)' "$name" "$n" "$n")"; NFAIL=$((NFAIL + 1))
    diff "expected/$n" "out/canon/$n" | head -8 | sed 's/^/        /' || true
  fi
done < expected/HASHES.txt
if [[ $TASKFAIL != 0 ]]; then
  report "FAIL  search_tasks                at least one search task exited with an error (see out/logs/)"; NFAIL=$((NFAIL + 1))
fi

T1=$(date +%s)
echo
report "checks passed: $NPASS, failed: $NFAIL, wall time $((T1 - T0)) s (JOBS=$JOBS)"
if [[ $NFAIL == 0 ]]; then
  report "OVERALL: PASS -- both implementations reproduce every recorded count; 0 AME(8,4) frames in all ten V_0 classes."
  exit 0
else
  report "OVERALL: FAIL"
  exit 1
fi
