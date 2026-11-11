#!/usr/bin/env bash
# Replay the release.  No network access is needed.
#
# Usage:  ./replay.sh [quick|extra|board|headline|full]    (default: quick)
#
#   quick     checksums, verifier build + tests + board record #20, the validation suite of the
#             128-bit evaluator (with a test of its series-tail bound), its re-certification of
#             the final X-Z and two-Pauli board records, the second (kernel) evaluator on the
#             three headline inputs in double precision, and the floating-point master-equation
#             evaluator.  6 to 8 minutes (4 processes).
#   extra     the frozen verifier on the other three final board records and on the depolarizing
#             headline input (which it cannot decide), the negative controls, the points with
#             margin and the n = 120 / n = 240 X-Z inputs with the 128-bit evaluator, and the
#             second evaluator at 30 digits on the n = 45, 80 and 120 inputs.  About 4 hours with
#             NPROC=2 on a shared machine.
#   board     the 128-bit evaluator on the 17 board records that are not final records; each
#             enclosure must lie inside the frozen verifier's.  2.7 hours with NPROC=2.
#   headline  the three headline values with the 128-bit evaluator, one after the other:
#             depolarizing n = 480 (2.9 to 3.4 h), X-Z n = 320 (1.3 to 2.2 h), two-Pauli n = 320
#             (1.7 to 2.9 h); times of the runs of October and November 2026 with the process
#             counts below.
#   full      quick, extra, board, headline.
#
# Requirements: bash, a C compiler (cc), and Python with the packages of requirements.txt (the
# reference run used Python 3.14.3; python-flint 0.9.0 needs a recent Python).  Choose it with
# PYTHON=/path/to/python (default python3).  NPROC overrides the number of processes of every
# 128-bit run (defaults: 3 for depolarizing n = 480 and two-Pauli n = 320, 4 otherwise; a
# depolarizing n = 480 process needs about 2.5 GB).  ONLY="step1 step2" runs only the named steps.
#
# Outputs go to out/ (logs, one file per step) and out/SUMMARY.txt.  Each step is compared with
# expected/ by check_outputs.py.  Exit status: 0 if every step passed, 1 otherwise.
set -uo pipefail
cd "$(dirname "$0")"
MODE="${1:-quick}"
PY="${PYTHON:-python3}"
OUT=out
mkdir -p "$OUT"
SUMMARY="$OUT/SUMMARY.txt"
status=0
B=certificates/board
A=certificates/arb
export PYTHONPATH="$PWD/arb:$PWD/arb/validation:$PWD/indep:$PWD/second${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1

say() { echo "$*" | tee -a "$SUMMARY"; }

want() {            # want STEP: true if the step is selected
  [ -z "${ONLY:-}" ] && return 0
  case " $ONLY " in *" $1 "*) return 0 ;; esac
  return 1
}

finish() {          # finish STEP SECONDS CHECK-COMMAND...
  local step="$1" secs="$2"; shift 2
  local line
  line="$("$@" 2>&1)"; local rc=$?
  [ $rc -ne 0 ] && status=1
  say "$(printf '%-22s %6ss  %s' "$step" "$secs" "$line")"
}

arb_run() {         # arb_run STEP EXPECTED-NAME INPUT NPROC [extra arb_ic.py arguments]
  local step="$1" name="$2" input="$3" np="${NPROC:-$4}"; shift 4
  want "$step" || return 0
  local t0=$(date +%s)
  "$PY" arb/arb_ic.py "$input" --nproc "$np" --json "$OUT/$step.json" "$@" > "$OUT/$step.log" 2>&1
  finish "$step" $(( $(date +%s) - t0 )) "$PY" check_outputs.py arb "$name" "$OUT/$step.log"
}

verifier_run() {    # verifier_run STEP EXPECTED-NAME INPUT
  local step="$1" name="$2" input="$3"
  want "$step" || return 0
  local t0=$(date +%s)
  "$PY" verifier/verify.py --json "$input" > "$OUT/$step.json" 2> "$OUT/$step.err"
  finish "$step" $(( $(date +%s) - t0 )) "$PY" check_outputs.py verifier "$name" "$OUT/$step.json"
}

text_run() {        # text_run STEP EXPECTED-NAME COMMAND...
  local step="$1" name="$2"; shift 2
  want "$step" || return 0
  local t0=$(date +%s)
  "$@" > "$OUT/$step.log" 2>&1
  finish "$step" $(( $(date +%s) - t0 )) "$PY" check_outputs.py text "$name" "$OUT/$step.log"
}

second_run() {      # second_run STEP EXPECTED-NAME COMMAND...
  local step="$1" name="$2"; shift 2
  want "$step" || return 0
  local t0=$(date +%s)
  "$@" > "$OUT/$step.log" 2>&1
  finish "$step" $(( $(date +%s) - t0 )) "$PY" check_outputs.py second "$name" "$OUT/$step.log"
}

part_quick() {
  if want sums; then
    if command -v shasum >/dev/null; then SHA="shasum -a 256"; else SHA="sha256sum"; fi
    if $SHA -c SHA256SUMS --quiet > "$OUT/sums.log" 2>&1; then
      say "sums                        PASS SHA256SUMS: every released file unchanged"
    else
      say "sums                        FAIL SHA256SUMS (see $OUT/sums.log)"; status=1
    fi
  fi
  if want env; then
    "$PY" check_outputs.py env 2>&1 | tee "$OUT/env.log" | sed 's/^/env                         /' | tee -a "$SUMMARY"
    (cc --version 2>&1 | head -1) | sed 's/^/env                         /' | tee -a "$SUMMARY"
  fi
  if want verifier_build; then
    # the engine is built with the flags of verifier/verify.py and verifier/freeze.sh; the pinned
    # hashes of the three sources are compared in every verifier step (check_outputs.py)
    rm -f verifier/icball
    if cc -O2 -ffp-contract=off -fno-fast-math -o verifier/icball verifier/icball.c -lm > "$OUT/verifier_build.log" 2>&1; then
      say "verifier_build              PASS engine built from verifier/icball.c"
    else
      say "verifier_build              FAIL cc (see $OUT/verifier_build.log)"; status=1
    fi
  fi
  text_run verifier_tests verifier_tests "$PY" verifier/tests/test_verifier.py
  verifier_run verifier_rec20 rec20 $B/20261027-053852_lead_9e078e22e51cd14b.json
  text_run arb_test_eig arb_test_eig "$PY" arb/validation/test_eig.py
  text_run arb_validate_brute arb_validate_brute "$PY" arb/validation/validate_brute.py
  text_run arb_test_tail_small arb_test_tail_small "$PY" arb/validation/test_tail.py --small-only
  text_run arb_validate_kg arb_validate_kg "$PY" arb/validation/validate_records.py --only b --nproc "${NPROC:-4}"
  arb_run arb_rec28 final_xz_n80 $B/20261027-124857_lead_2a327ea83c2b10f7.json 4
  arb_run arb_rec27 final_2p_n80 $B/20261027-124857_lead_585247a2e5ee8280.json 4
  text_run second_validate second_validate "$PY" second/validate_kernel.py
  second_run second_headline second_headline "$PY" second/run_headline.py
  second_run indep_final indep_final "$PY" indep/run_final.py
}

part_extra() {
  verifier_run verifier_rec28 rec28 $B/20261027-124857_lead_2a327ea83c2b10f7.json
  verifier_run verifier_rec27 rec27 $B/20261027-124857_lead_585247a2e5ee8280.json
  verifier_run verifier_rec30 rec30 $B/20261028-024810_lead_5f09c5df10a5804f.json
  verifier_run verifier_dep480_headline dep480_headline_undecided $A/depolarizing_n480_0.06688201.json
  text_run arb_validate_brute_extra arb_validate_brute_extra "$PY" arb/validation/validate_brute.py --hp-extra
  text_run arb_test_tail arb_test_tail "$PY" arb/validation/test_tail.py --nproc "${NPROC:-4}"
  text_run arb_negative_p arb_negative_p "$PY" arb/validation/validate_records.py --only d --nproc "${NPROC:-4}"
  if want make_perturbed; then
    # the perturbed inputs are regenerated from the board records and must carry the same family,
    # p and amplitudes as the shipped ones (the free-text note names the source by its path and is
    # not compared); the shipped files are then restored, so the release tree stays unchanged
    mkdir -p "$OUT/pert_shipped" "$OUT/pert_regenerated"
    cp arb/validation/testdata/pert_*.json "$OUT/pert_shipped/"
    "$PY" arb/validation/make_perturbed.py > "$OUT/make_perturbed.log" 2>&1
    cp arb/validation/testdata/pert_*.json "$OUT/pert_regenerated/"
    cp "$OUT"/pert_shipped/pert_*.json arb/validation/testdata/
    for t in xz_n80 2p_n80 dep_n240; do
      finish "make_perturbed_$t" 0 "$PY" check_outputs.py pert "$OUT/pert_shipped/pert_$t.json" "$OUT/pert_regenerated/pert_$t.json"
    done
  fi
  arb_run arb_pert_xz_n80 pert_xz_n80 arb/validation/testdata/pert_xz_n80.json 4
  arb_run arb_pert_2p_n80 pert_2p_n80 arb/validation/testdata/pert_2p_n80.json 4
  arb_run arb_pert_dep_n240 pert_dep_n240 arb/validation/testdata/pert_dep_n240.json 4
  arb_run arb_margin_xz_n80 cons_independent_xz_n80_0.12463577 $A/independent_xz_n80_0.12463577.json 4
  arb_run arb_margin_2p_n80 cons_two_pauli_n80_0.12366746 $A/two_pauli_n80_0.12366746.json 4
  arb_run arb_margin_dep_n240 cons_depolarizing_n240_0.06666451 $A/depolarizing_n240_0.06666451.json 4
  arb_run arb_rec20 final_dep_n240 $B/20261027-053852_lead_9e078e22e51cd14b.json 4
  arb_run arb_xz_n120 arb_independent_xz_n120_0.1259_board_path1 $A/independent_xz_n120_0.1259_board_path1.json 4
  arb_run arb_xz_n240 arb_independent_xz_n240_0.12737956 $A/independent_xz_n240_0.12737956.json 4
  second_run second_mp second_mp "$PY" second/run_mp.py
}

part_board() {
  # the 17 board records that are not among the four final ones (those are in quick and extra, and
  # record #30 is in arb/logs/): 128-bit evaluator, each enclosure must lie inside the verifier's
  arb_run arb_board_02 board_02 $B/20261026-210916_lead_481ce3870c2d4e21.json 4
  arb_run arb_board_03 board_03 $B/20261026-210926_lead_c9b944ed9a6f8a4a.json 4
  arb_run arb_board_04 board_04 $B/20261026-210955_lead_ad42e025717074f4.json 4
  arb_run arb_board_06 board_06 $B/20261026-230622_lead_32025a4e90f6200f.json 4
  arb_run arb_board_09 board_09 $B/20261026-231833_lead_cc6a44c7cf5414bb.json 4
  arb_run arb_board_10 board_10 $B/20261026-232813_lead_9c563111d96087e5.json 4
  arb_run arb_board_11 board_11 $B/20261027-010137_lead_5a8d8b8895df26cb.json 4
  arb_run arb_board_12 board_12 $B/20261027-010209_lead_ad2fcf9c940ceb22.json 4
  arb_run arb_board_13 board_13 $B/20261027-010328_lead_713c6cf1f278172f.json 4
  arb_run arb_board_15 board_15 $B/20261027-033842_lead_38eff8cefeab3c5d.json 4
  arb_run arb_board_16 board_16 $B/20261027-045237_lead_673b3e61d692d1a4.json 4
  arb_run arb_board_17 board_17 $B/20261027-045237_lead_4675821f074f5b91.json 4
  arb_run arb_board_18 board_18 $B/20261027-050500_lead_ba11e08853949e86.json 4
  arb_run arb_board_21 board_21 $B/20261027-074257_lead_d1c9e883b77efe7c.json 4
  arb_run arb_board_22 board_22 $B/20261027-074257_lead_834edf5841fffd54.json 4
  arb_run arb_board_24 board_24 $B/20261027-101825_lead_46b0d3127de5b834.json 4
  arb_run arb_board_25 board_25 $B/20261027-101825_lead_6b2ddb012ce45bff.json 4
}

part_headline() {
  arb_run headline_dep_n480 arb_depolarizing_n480_0.06688201 $A/depolarizing_n480_0.06688201.json 3
  arb_run headline_xz_n320 arb_independent_xz_n320_0.12770564 $A/independent_xz_n320_0.12770564.json 4
  arb_run headline_2p_n320 arb_two_pauli_n320_0.12677472 $A/two_pauli_n320_0.12677472.json 3
}

say "=== replay $MODE, $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
case "$MODE" in
  quick)    part_quick ;;
  extra)    part_extra ;;
  board)    part_board ;;
  headline) part_headline ;;
  full)     part_quick; part_extra; part_board; part_headline ;;
  *) echo "usage: ./replay.sh [quick|extra|board|headline|full]"; exit 2 ;;
esac
if [ $status -eq 0 ]; then say "RELEASE ($MODE): PASS"; else say "RELEASE ($MODE): FAIL"; fi
exit $status
