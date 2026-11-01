#!/usr/bin/env bash
# Table VI row l=2, sigma=(2^4), auditor implementation (inv4dfs.c, own reduction C + symmetric beta), quick subset:
# 3 random AME(4,4) graph codes (seed 21, = first 3 lines of work/auditor/inv4c_21.log), the first 20 of the 250
# F_16-Hermitian MDS codes (seed 6, = first 20 lines of inv4h_6.log), all 2^36 beta each: 0 solutions;
# reduction validation inv4_validate.py 5 (random beta: S isotropic, rank 16, direct 2^16-word min weight == formula). numpy.
set -euo pipefail
"$PYNP" "$ROOT/auditor/inv4c.py" 21 3 1e9
sed -E 's/, "t": [0-9.]+//; s/"elapsed": [0-9.]+, //; s/, "elapsed": [0-9.]+//' inv4c_21.log
"$PYNP" "$ROOT/auditor/inv4h.py" 20 6
sed -E 's/, "t": [0-9.]+//; s/"elapsed": [0-9.]+, //; s/, "elapsed": [0-9.]+//' inv4h_6.log
"$PYNP" "$ROOT/auditor/inv4_validate.py" 5
