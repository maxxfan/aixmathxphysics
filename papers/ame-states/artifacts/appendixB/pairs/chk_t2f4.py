"""chk_t2f4.py: independent check: all 4-dim subspaces of F_2^8 = (F_2^2)^4 with min block weight >= 3 (additive MDS (4,16,3)_4),
classified by the rank of the Gram matrix of b = sum_j b_j (b_j symplectic form on each 2-dim block: b_j(u,v) = u0 v1 + u1 v0)."""
import sys
allowed = set()
for f in range(1, 256):
    if sum(1 for j in range(4) if (f >> (2 * j)) & 3) >= 3:
        allowed.add(f)


def b(f, g):
    r = 0
    for j in range(4):
        u0, u1 = (f >> (2 * j)) & 1, (f >> (2 * j + 1)) & 1
        v0, v1 = (g >> (2 * j)) & 1, (g >> (2 * j + 1)) & 1
        r ^= (u0 & v1) ^ (u1 & v0)
    return r


level = {frozenset([0])}
for d in range(4):
    nxt = set()
    for S in level:
        for v in allowed:
            if v in S:
                continue
            if all((v ^ s) in allowed or (v ^ s) == 0 for s in S):
                nxt.add(frozenset(S | {v ^ s for s in S}))
    level = nxt
    print("dim", d + 1, "subspaces", len(level), flush=True)


def rank_gram(S):
    # radical dimension: elements orthogonal to all of S inside S
    rad = [z for z in S if all(b(z, w) == 0 for w in S)]
    return len(rad)


from collections import Counter
cnt = Counter(rank_gram(S) for S in level)
print("radical sizes (1=nondegenerate):", dict(cnt))
