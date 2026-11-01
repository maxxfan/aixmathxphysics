"""Auditor: counting constraints for a stabilizer AME(n,4) (party = F_2^4, S = F_2-subspace of dim 2n, every nonzero
element of party weight >= dm = n/2+1, projection of S onto any <= dm-1 parties is onto) invariant under a LOCAL LINEAR map
Omega = sum omega_j, omega_j in GL(4,2), Omega nontrivial of prime order l in {2,3,5,7}; NO symplecticity assumed.
Odd l: V_j = U_j (fixed, F_2-dim u_j) + W_j (non-fixed, w_j = 4 - u_j); S = S^Omega + S_nf, S^Omega on parties with u_j > 0,
  S_nf on parties with w_j > 0 (S_nf != 0 because Omega|S = 1 would force Omega = 1).
   l = 3: (u,w) in (4,0),(2,2),(0,4);  l = 5: (4,0),(0,4);  l = 7: (4,0),(1,3)
l = 2: N = Omega - 1, N^2 = 0, rank rho_j in 0,1,2; N S in S, S = (S cap ker N) + ... : dim S = dim(S cap ker N) + dim(N S),
  S cap ker N lives in alphabets ker nu_j (dim 4 - rho_j) on all parties, N S lives in alphabets Im nu_j (dim rho_j) on parties with rho_j > 0.
Every summand is a subcode of S, so each nonzero element has party weight >= dm and dim <= the sum of the (m - dm + 1) smallest
alphabet dimensions over its m parties (0 if m < dm).  Print the type assignments with total bound >= 2n and nonzero second summand.
usage: python3 work/auditor/fixpt_count.py [n]      (default 8; try 4 and 6 as sanity checks: AME(4,4), AME(6,4) exist)
"""
import itertools
import sys

n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
dm = n // 2 + 1
target = 2 * n

# each type = (dim of summand A alphabet, dim of summand B alphabet); summand B must be nonzero
TYPES = {
    2: [(4, 0), (3, 1), (2, 2)],   # (dim ker N, dim Im N)
    3: [(4, 0), (2, 2), (0, 4)],
    5: [(4, 0), (0, 4)],
    7: [(4, 0), (1, 3)],
}


def single_bound(dims):
    dims = sorted(d for d in dims if d > 0)
    m = len(dims)
    return 0 if m < dm else sum(dims[: m - (dm - 1)])


print("n =", n, " dmin =", dm, " dim S =", target)
for l, types in TYPES.items():
    seen = set()
    for assign in itertools.combinations_with_replacement(range(len(types)), n):
        a = [types[t][0] for t in assign]
        b = [types[t][1] for t in assign]
        ba, bb = single_bound(a), single_bound(b)
        if bb > 0 and ba + bb >= target:
            key = tuple(sorted(types[t] for t in assign))
            if key not in seen:
                seen.add(key)
                print("  l =", l, "feasible by Singleton counting:", key, " max dims (A,B) =", ba, bb)
    if not seen:
        print("  l =", l, "none")
