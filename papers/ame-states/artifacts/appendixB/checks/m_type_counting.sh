#!/usr/bin/env bash
# Prop. prop:loc, (0,4)^8 with M-type parties (uniform M, M mixed with H/S): the auditor's chain-ring / Singleton counting
# (argument.md 5.2, hand derivation) made executable by extra/m_type_counting.py (written for the release).
set -euo pipefail
"$PY" "$ROOT/extra/m_type_counting.py" | tee /dev/stderr
