#!/usr/bin/env bash
# Prop. prop:loc, uniform H and S types, PAIRS implementation (pairs/dfs4b.c: F_4 graph states with Hermitian (type 0 = H) or
# symmetric (type 1 = S) adjacency, row 0 normalised to m ones, m = n/2..n-1): calibration runs at n = 4, 6 (both types)
# and n = 7 (type 0), which must find solutions.  The n = 8 runs are the full-tier checks hs_pairs_n8_t*_m*.
set -euo pipefail
run() { "$BUILD/dfs4b" "$@" 1000000000 | tee /dev/stderr | grep '^DONE' | sed -E 's/ time=[0-9.]+s//'; }
for t in 0 1; do for n in 4 6; do for m in $(seq $((n / 2)) $((n - 1))); do run $n $t $m 0 1; done; done; done
for m in 3 4 5 6; do run 7 0 $m 0 1; done
# dfs4.c (the same search, minors tested in the opposite order) must give identical lines at n = 4, 6:
same=yes
for t in 0 1; do for n in 4 6; do for m in $(seq $((n / 2)) $((n - 1))); do
    a="$("$BUILD/dfs4" $n $t $m 0 1 1000000000 | grep '^DONE' | sed -E 's/ time=[0-9.]+s//')"
    b="$("$BUILD/dfs4b" $n $t $m 0 1 1000000000 | grep '^DONE' | sed -E 's/ time=[0-9.]+s//')"
    [ "$a" = "$b" ] || same=no
done; done; done
echo "dfs4.c and dfs4b.c give identical DONE lines at n = 4, 6: $same"
