#!/usr/bin/env bash
# Table VI row l=2, sigma=(2^4): all F_2-linear stabilizer AME(4,4) codes form ONE orbit of G = Sp(4,2)^4 x| S_4:
# N_tot = 4,478,976,000 = |G| / |Aut(C0)| with |Aut(C0)| = 1440.  Auditor's ame44_classes.py (numpy), ~30 s.
set -euo pipefail
"$PYNP" "$ROOT/auditor/ame44_classes.py" > raw.txt
cat raw.txt >&2
grep -vE '^(rep|perm) ' raw.txt
# per-coset-representative counts n(A) and per-permutation local-map counts, summarised
awk '/^rep /{c[$NF]++} END{for (k in c) printf "coset reps with n(A) = %s: %d\n", k, c[k]}' raw.txt | sort
awk '/^perm /{c[$NF]++} END{for (k in c) printf "position permutations with %s local maps: %d\n", k, c[k]}' raw.txt
