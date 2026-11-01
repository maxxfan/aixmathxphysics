#!/usr/bin/env bash
# Prop. prop:loc: "calibrated ... against random sampling at intermediate depths".  Monte-Carlo estimates of the normalised
# counts (auditor hdfs.c / hsdfs.c 'sample' and 'prefix' modes, splitmix64, seed 1) compared with the exact DFS counts.
# Canonical output: whether |estimate - exact| <= 4 sigma (the estimates themselves are printed to the log only).
set -euo pipefail
chk() {  # chk EXACT label cmd...
    exact="$1"; label="$2"; shift 2
    line="$("$@" 2>&1 | tail -1)"; echo "$line" >&2
    "$PY" - "$exact" "$label" "$line" <<'PYEOF'
import re, sys
exact, label, line = float(sys.argv[1]), sys.argv[2], sys.argv[3]
m = re.search(r"= ([0-9.]+) \(\+- ([0-9.]+)\)\s*$", line)
est, sd = float(m.group(1)), float(m.group(2))
print("%-38s exact %-8g within 4 sigma: %s" % (label, exact, "yes" if abs(est - exact) <= 4 * sd else "NO"))
PYEOF
}
chk 6      "H h=2 N_norm (sample)"            "$BUILD/hdfs" 2 sample 20000 1
chk 2074   "H h=3 N_norm (sample)"            "$BUILD/hdfs" 3 sample 400000 1
chk 64     "H h=3 starts (prefix P=0)"        "$BUILD/hdfs" 3 prefix 400000 1 0
chk 2074   "H h=3 party-1 nodes (prefix P=1)" "$BUILD/hdfs" 3 prefix 400000 1 1
chk 624    "H h=4 starts (prefix P=0)"        "$BUILD/hdfs" 4 prefix 1000000 1 0
chk 778584 "H h=4 party-1 nodes (prefix P=1)" "$BUILD/hdfs" 4 prefix 2000000 1 1
chk 2      "S h=2 N_norm (sample)"            "$BUILD/hsdfs" 2 sample 300000 1
chk 380    "S h=3 N_norm (sample)"            "$BUILD/hsdfs" 3 sample 1000000 1
chk 20     "S h=4 starts (prefix P=0)"        "$BUILD/hsdfs" 4 prefix 1000000 1 0
chk 18720  "S h=4 party-1 nodes (prefix P=1)" "$BUILD/hsdfs" 4 prefix 10000000 1 1
