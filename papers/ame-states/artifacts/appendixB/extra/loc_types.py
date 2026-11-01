#!/usr/bin/env python3
"""[release] Written for the Appendix-B release (NOT part of the original run; stdlib only).

Second, independently coded enumeration of the counting step of Proposition prop:loc (the original is the auditor's
auditor/fixpt_count.py, argument.md 5.3).  Setting: S in V = (F_2^4)^n, dim S = 2n, every nonzero element of party weight
>= n/2+1 (AME), Omega = (+) omega_j a LOCAL linear map of prime order l.

  odd l : S = S^Omega (+) S_nf, S^Omega in (+) U_j (U_j = ker(omega_j - 1), dim u_j), S_nf in (+) W_j (dim w_j = 4 - u_j),
          S_nf != 0 (otherwise Omega acts trivially on S, hence Omega = 1);
  l = 2 : N = Omega - 1, N^2 = 0, rank r_j in {0,1,2}; dim S = dim(S cap ker N) + dim(N S), S cap ker N in (+) ker N_j
          (dim 4 - r_j), N S in (+) im N_j (dim r_j), N S != 0.

Both summands are subcodes of S, so every nonzero word has party weight >= dmin = n/2+1 and the projection onto ANY
n - dmin + 1 parties is injective: dim <= min over (n-dmin+1)-subsets P of sum_{j in P} alphabet_j.  Unlike fixpt_count.py
(multisets, "sum of the smallest alphabets"), this script runs over ALL ordered type assignments and takes the minimum over
all party subsets explicitly.  An assignment survives iff bound(A) + bound(B) >= 2n and bound(B) > 0.
Usage: loc_types.py [n]   (default: 8, then the sanity cases n = 4 and n = 6, where AME(4,4) and AME(6,4) exist)
"""
import itertools
import sys

TYPES = {  # l: per-party (alphabet dim of summand A, alphabet dim of summand B)
    2: {"r=0": (4, 0), "r=1": (3, 1), "r=2": (2, 2)},
    3: {"(4,0)": (4, 0), "(2,2)": (2, 2), "(0,4)": (0, 4)},
    5: {"(4,0)": (4, 0), "(0,4)": (0, 4)},
    7: {"(4,0)": (4, 0), "(1,3)": (1, 3)},
}


def bound(alph, n, dmin):
    k = n - dmin + 1                      # number of parties onto which a subcode projects injectively
    return min(sum(alph[j] for j in P) for P in itertools.combinations(range(n), k))


def run(n):
    dmin = n // 2 + 1
    print("n = %d, min party weight %d, dim S = %d" % (n, dmin, 2 * n))
    for l, types in TYPES.items():
        names = list(types)
        surv = {}
        total = 0
        for assign in itertools.product(names, repeat=n):
            total += 1
            a = [types[t][0] for t in assign]
            b = [types[t][1] for t in assign]
            ba, bb = bound(a, n, dmin), bound(b, n, dmin)
            if bb > 0 and ba + bb >= 2 * n:
                key = tuple(sorted(set(assign)))
                surv.setdefault(key, set()).add((ba, bb))
        desc = "; ".join("%s^%d (bounds A+B = %s)" % ("+".join(k), n, ",".join("%d+%d" % x for x in sorted(v)))
                         for k, v in sorted(surv.items())) if surv else "none"
        print("  l = %d: %d ordered assignments, survivors: %s" % (l, total, desc))


if __name__ == "__main__":
    ns = [int(sys.argv[1])] if len(sys.argv) > 1 else [8, 4, 6]
    for n in ns:
        run(n)
