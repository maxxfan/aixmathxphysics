#!/usr/bin/env bash
# Prop. prop:loc, counting step, second enumeration written for the release (extra/loc_types.py): all ORDERED type
# assignments, explicit minimum over party subsets.  Must leave the same patterns as auditor/fixpt_count.py.
set -euo pipefail
"$PY" "$ROOT/extra/loc_types.py" | tee /dev/stderr
