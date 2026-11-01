#!/usr/bin/env bash
# Table VI, row "l = 3, sigma = (3^2,1^2), h with 2-dim fixed space" (a = 2 sector), AUDITOR implementation
# (argument.md 5.5(c), board #168): z3twist.py/z3twist2.py (imported by z3twist4.py) enumerate the 270 X_A lines,
# 3 key classes of 90 (Sp(4,2)_A transitive), the 3 representative pairs and their 384 completions T each (1152 T),
# export Vs = {s in V^Omega : wt(s+t) >= 5 for all t in T} (2928 of 4096 for every T); lagsearch2.c (staged DFS)
# finds 0 Lagrangians L of V^Omega inside Vs for all 1152 T; lag_control.py is the positive control
# (planted Lagrangian found; one word removed -> 0).
set -euo pipefail
"$PYNP" "$ROOT/auditor/z3twist4.py" nonfpf > z3twist4.out 2>&1 || { cat z3twist4.out >&2; exit 1; }
cat z3twist4.out >&2
echo "== z3twist4.py nonfpf (steps 1-2 and export)"
grep -vE '^exported ' z3twist4.out | sed -E 's/ +time [0-9.]+s$//'
echo "== sha256 of the exported Vs file (original work/auditor/z3twist_Vs_nonfpf.bin)"
"$PY" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' z3twist_Vs_nonfpf.bin
echo "== lagsearch2 z3twist_Vs_nonfpf.bin 0 1151 (count x identical per-T lines, T index dropped)"
"$BUILD/lagsearch2" z3twist_Vs_nonfpf.bin 0 1151 > lag2.out
cat lag2.out >&2
sed -E 's/^T [0-9]+ //' lag2.out | sort | uniq -c | awk '{c=$1; $1=""; print c " x" $0}'
echo "== lag_control.py (positive control)"
"$PY" "$ROOT/auditor/lag_control.py" | tee /dev/stderr
