#!/usr/bin/env bash
# Replay the whole release. No network access is needed.
# Requirements: bash, a C compiler (cc), Python 3.9+. 24 of the 60 checks in appendixB/ also need numpy: they use
#   $PYTHON_NUMPY, or python3 if it imports numpy; without numpy they are skipped and the verdict is INCOMPLETE.
# Usage: ./replay.sh [quick|full]
#   quick: about 40 CPU-minutes; JOBS=4 ./replay.sh quick took 14.6 minutes of wall time (JOBS=1: about 45 minutes)
#   full:  about 2.5 CPU-hours, about 1 hour with JOBS=4
#   frames/ and appendixB/ run JOBS checks in parallel (frames/ defaults to 4, appendixB/ to 1)
# Each part can also be run on its own: states/run.sh, coupling/run.sh, frames/run.sh, appendixB/run.sh, ancillary/run.sh.
# Exit status: 0 all PASS, 1 a FAIL, 2 no FAIL but appendixB skipped its numpy checks.
set -uo pipefail
cd "$(dirname "$0")"
MODE="${1:-quick}"
status=0
incomplete=0
summary=""
if command -v shasum >/dev/null; then SHA="shasum -a 256"; else SHA="sha256sum"; fi
if $SHA -c SHA256SUMS --quiet; then summary+="PASS  SHA256SUMS (every released file unchanged)\n"; else summary+="FAIL  SHA256SUMS\n"; status=1; fi
for part in states coupling frames appendixB ancillary; do
  if [ -x "$part/run.sh" ]; then
    echo "=== $part ($MODE) ==="
    start=$(date +%s)
    "$part/run.sh" "$MODE"; rc=$?
    case $rc in
      0) r=PASS ;;
      2) if [ "$part" = appendixB ]; then r="INCOMPLETE (numpy checks skipped)"; incomplete=1; else r=FAIL; status=1; fi ;;
      *) r=FAIL; status=1 ;;
    esac
    summary+="$r  $part ($(( $(date +%s) - start )) s)\n"
  else
    summary+="FAIL  $part/run.sh missing\n"; status=1
  fi
done
echo "=== summary ==="
printf "%b" "$summary"
if [ $status != 0 ]; then echo "RELEASE: FAIL"; exit 1; fi
if [ $incomplete != 0 ]; then echo "RELEASE: INCOMPLETE (no FAIL; install numpy or set PYTHON_NUMPY to run every check)"; exit 2; fi
echo "RELEASE: PASS"
exit 0
