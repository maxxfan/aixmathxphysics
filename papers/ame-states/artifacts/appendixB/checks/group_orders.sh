#!/usr/bin/env bash
# Group facts (Prop. prop:loc first line; Prop. prop:aut (b) "Sp(4,2) = S_6 has no element of order 7"):
# auditor/sp42_classes.py (original) and extra/group_orders.py (written for the release).
set -euo pipefail
"$PY" "$ROOT/auditor/sp42_classes.py" | tee /dev/stderr
"$PY" "$ROOT/extra/group_orders.py" | tee /dev/stderr
