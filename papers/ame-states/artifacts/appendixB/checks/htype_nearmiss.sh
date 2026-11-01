#!/usr/bin/env bash
# Prop. prop:loc: "The best H-type state has one bad cut out of 35."  auditor/hdfs.c with the k = 3 minors dropped at
# party 2 (KMAX = 2, STOPAT = 2): the unique completion of each of the 322608 near-solutions fails 2, 4 or 6 of the 70
# party 4-sets (= 1, 2, 3 of the 35 cuts), never 0.
set -euo pipefail
"$BUILD/hdfs" 4 - 2 2 2>&1 | tee /dev/stderr | grep -vE '^rc=|time='
