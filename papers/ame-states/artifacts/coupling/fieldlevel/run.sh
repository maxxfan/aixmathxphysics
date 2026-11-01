#!/bin/sh
# Run the field-level checks of the coupling ingredients and compare with the recorded outputs in expected/.
# Exit status 0 only if every script exits 0 and every output is identical to expected/<name>.out.
# usage: sh run.sh        (PYTHON=/path/to/python3 sh run.sh to choose the interpreter; Python 3.9+, stdlib only)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE" || exit 2
PY=${PYTHON:-python3}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/fieldlevel.XXXXXX") || exit 2
trap 'rm -rf "$TMP"' EXIT
status=0

run() {   # name, then the command
    name=$1
    shift
    "$@" > "$TMP/$name.out" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "$name: exit status $rc"
        status=1
    fi
    if diff -u "expected/$name.out" "$TMP/$name.out" > "$TMP/$name.diff"; then
        echo "$name: exit $rc, output identical to expected/$name.out ($(wc -l < "$TMP/$name.out" | tr -d ' ') lines)"
    else
        echo "$name: OUTPUT DIFFERS from expected/$name.out:"
        cat "$TMP/$name.diff"
        status=1
    fi
}

run convert "$PY" -B convert.py
run check_fieldlevel "$PY" -B check_fieldlevel.py
run derive_check "$PY" -B derive.py --check

grep -h -e "^ingredients:" -e "^cells:" "$TMP/convert.out" "$TMP/check_fieldlevel.out"
if [ $status -eq 0 ]; then
    echo "ALL OK"
else
    echo "FAILED"
fi
exit $status
