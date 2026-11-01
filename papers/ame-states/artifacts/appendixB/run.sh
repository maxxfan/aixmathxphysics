#!/usr/bin/env bash
# Replay of the computations behind Appendix B ("Symmetry sectors (cross-check)") of the p10 AME manuscript.
#
#   ./run.sh            same as ./run.sh quick
#   ./run.sh quick      every check of tier "quick" (target: < 20 min on one laptop core)
#   ./run.sh full       every check (quick + full tiers)
#   ./run.sh list       list the checks (id, tier, needs, Appendix-B claim)
#   ./run.sh ID [ID..]  run the named checks only
#   ./run.sh clean      remove build/, out/, work/
#
# Environment:
#   PYTHON=python3          interpreter for the stdlib-only scripts (Python >= 3.9)
#   PYTHON_NUMPY=...        interpreter with numpy >= 1.20 for the scripts that need numpy
#                           (default: $PYTHON if it can import numpy; otherwise those checks are SKIPPED)
#   CC=cc                   C compiler (compiled with -O2 into build/)
#   JOBS=1                  number of checks run in parallel (each check is a single process, < 1.5 GB)
#
# Every check writes its canonical output to out/ID.txt (raw log: out/ID.log) and is compared, after dropping
# '#' comment lines and blank lines, with expected/ID.txt. Exit status: 0 = all PASS, 1 = some FAIL,
# 2 = no FAIL but some checks SKIPPED (missing numpy).  No network access is used.
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(pwd -P)"
BUILD="$ROOT/build"
OUT="$ROOT/out"
WORK="$ROOT/work"
LIST="$ROOT/checks/LIST"
export ROOT BUILD OUT WORK

PY="${PYTHON:-python3}"
if [ -n "${PYTHON_NUMPY:-}" ]; then
    PYNP="$PYTHON_NUMPY"
elif "$PY" -c 'import numpy' >/dev/null 2>&1; then
    PYNP="$PY"
else
    PYNP=""
fi
CC="${CC:-cc}"
JOBS="${JOBS:-1}"
export PY PYNP CC
export AME_RANKCHECK="$BUILD/rankcheck"
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

# ---------------------------------------------------------------- registry helpers
# checks/LIST: one check per line, fields separated by '|':  ID | tier (quick/full) | needs (py/np/c) | claim
list_ids() {   # $1 = quick | full
    grep -v '^#' "$LIST" | sed '/^[[:space:]]*$/d' | while IFS='|' read -r id tier needs claim; do
        id="$(echo "$id" | tr -d ' ')"; tier="$(echo "$tier" | tr -d ' ')"
        if [ "$1" = full ] || [ "$tier" = quick ]; then echo "$id"; fi
    done
}
field() {      # field ID N
    grep -v '^#' "$LIST" | awk -F'|' -v id="$1" -v n="$2" '{g=$1; gsub(/ /,"",g); if (g==id) {f=$n; sub(/^ +/,"",f); sub(/ +$/,"",f); print f}}'
}

# ---------------------------------------------------------------- build
C_SOURCES="lib/rankcheck.c auditor/hdfs.c auditor/hsdfs.c auditor/inv4dfs.c auditor/lagsearch.c auditor/lagsearch2.c
pairs/dfs4.c pairs/dfs4b.c pairs/inv8.c pairs/qdfs.c pairs/scan_z3.c"
build() {
    mkdir -p "$BUILD"
    for src in $C_SOURCES; do
        exe="$BUILD/$(basename "$src" .c)"
        if [ ! -x "$exe" ] || [ "$src" -nt "$exe" ]; then
            "$CC" -O2 -o "$exe" "$src" -lm
        fi
    done
}

# ---------------------------------------------------------------- one check
run_one() {
    id="$1"
    if [ ! -f "$ROOT/checks/$id.sh" ]; then echo "unknown check: $id" >&2; return 1; fi
    needs="$(field "$id" 3)"
    mkdir -p "$OUT"
    rm -f "$OUT/$id.txt" "$OUT/$id.log" "$OUT/$id.diff" "$OUT/$id.status"
    if [ "$needs" = np ] && [ -z "$PYNP" ]; then
        echo "SKIP $id 0 (numpy not available: set PYTHON_NUMPY)" > "$OUT/$id.status"
        printf '%-6s %-28s %s\n' SKIP "$id" "(needs numpy; set PYTHON_NUMPY=/path/to/python)"
        return 0
    fi
    rm -rf "$WORK/$id"; mkdir -p "$WORK/$id"
    t0=$(date +%s)
    rc=0
    ( cd "$WORK/$id" && nice -n 10 bash "$ROOT/checks/$id.sh" ) > "$OUT/$id.txt" 2> "$OUT/$id.log" || rc=$?
    t1=$(date +%s); dt=$((t1 - t0))
    if [ "$rc" -ne 0 ]; then
        st=FAIL; why="(exit code $rc; see out/$id.log)"
    elif diff <(grep -v '^#' "$ROOT/expected/$id.txt" | sed '/^[[:space:]]*$/d') \
              <(grep -v '^#' "$OUT/$id.txt" | sed '/^[[:space:]]*$/d') > "$OUT/$id.diff"; then
        st=PASS; why=""; rm -f "$OUT/$id.diff"
    else
        st=FAIL; why="(output differs from expected/$id.txt; see out/$id.diff)"
    fi
    echo "$st $id $dt $why" > "$OUT/$id.status"
    printf '%-6s %-28s %6ss %s\n' "$st" "$id" "$dt" "$why"
}

# ---------------------------------------------------------------- driver
run_many() {
    ids="$*"
    build
    n=$(echo $ids | wc -w | tr -d ' ')
    echo "Appendix B replay: $n checks, JOBS=$JOBS, PYTHON=$PY ($("$PY" -V 2>&1)), PYTHON_NUMPY=${PYNP:-none}"
    echo "------------------------------------------------------------------"
    T0=$(date +%s)
    if [ "$JOBS" -gt 1 ]; then
        printf '%s\n' $ids | xargs -P "$JOBS" -n 1 "$ROOT/run.sh" --one
    else
        for id in $ids; do run_one "$id"; done
    fi
    T1=$(date +%s)
    npass=0; nfail=0; nskip=0
    for id in $ids; do
        st="$(cut -d' ' -f1 "$OUT/$id.status" 2>/dev/null || echo FAIL)"
        case "$st" in PASS) npass=$((npass+1));; SKIP) nskip=$((nskip+1));; *) nfail=$((nfail+1));; esac
    done
    echo "------------------------------------------------------------------"
    echo "PASS $npass  FAIL $nfail  SKIP $nskip  (wall time $((T1 - T0)) s)"
    if [ "$nfail" -gt 0 ]; then
        echo "OVERALL: FAIL"
        for id in $ids; do grep -q '^FAIL' "$OUT/$id.status" 2>/dev/null && echo "  failed: $id"; done || true
        return 1
    elif [ "$nskip" -gt 0 ]; then
        echo "OVERALL: INCOMPLETE (no failures; $nskip checks skipped because numpy is missing)"
        return 2
    fi
    echo "OVERALL: PASS"
    return 0
}

mode="${1:-quick}"
case "$mode" in
    --one) run_one "$2" ;;
    list)
        printf '%-28s %-6s %-5s %s\n' ID TIER NEEDS CLAIM
        grep -v '^#' "$LIST" | sed '/^[[:space:]]*$/d' | while IFS='|' read -r id tier needs claim; do
            printf '%-28s %-6s %-5s %s\n' "$(echo $id)" "$(echo $tier)" "$(echo $needs)" "$(echo $claim)"
        done ;;
    clean) rm -rf "$BUILD" "$OUT" "$WORK" ;;
    build) build ;;
    quick) run_many $(list_ids quick) ;;
    full) run_many $(list_ids full) ;;
    -h|--help|help) sed -n '2,22p' "$0" ;;
    *) run_many "$@" ;;
esac
