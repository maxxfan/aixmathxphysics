#!/usr/bin/env bash
# Table VI, row l = 5: the 1625 = 325 x 5 generators of H(5,16) that are NOT of graph type (dim L cap C = 2), excluded by hand
# in #184 / argument.md 5.6 because they contain the weight-2 word (0; u5, u6).  extra/c5_nongraph.py is a
# SUPPLEMENTARY stdlib script written for this release (not part of the original run); it enumerates them in the auditor's
# coordinates and checks the count identity 331500 + 1625 = 333125 = (q+1)(q^3+1)(q^5+1).
set -euo pipefail
"$PY" "$ROOT/extra/c5_nongraph.py" 0 0 2>&1 | tee /dev/stderr
