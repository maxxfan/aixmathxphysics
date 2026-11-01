#!/usr/bin/env bash
# Table VI row l = 3, (3,1^5), any twists: auditor/b3audit.py, counting over all 3^5 = 243 twist patterns: none survives.
# (pairs #175 derived it by hand; the auditor's script is the enumeration of that argument.)
set -euo pipefail
"$PY" "$ROOT/auditor/b3audit.py" | tee /dev/stderr
