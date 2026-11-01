"""val2_inv8.py : check inv8's DFS node counts per level against Monte-Carlo (independent coset enumeration) for the
first C of seed 1 in gen_inv8 (nodes 168 20128 343544 3842216 3377472 119040 96 0 at thr=4)."""
import sys, random, subprocess
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # [release] was the absolute work/pairs path
import gen_inv8 as g

d = [int(x, 16) for x in "4401 d802 8104 6208 8210 8820 b40 3680".split()]
Cwords = [0]
for dv in d:
    Cwords = Cwords + [c ^ dv for c in Cwords]
THR = int(sys.argv[1]) if len(sys.argv) > 1 else 4
rng = random.Random(5)
dp = {}
for x in range(1, 256):
    D = 0
    for i in range(8):
        if x >> i & 1:
            D ^= d[i]
    dp[x] = D


def ok_x(x, cols):
    Qx = 0
    for i in range(8):
        if x >> i & 1:
            Qx ^= cols[i] << 8
    D = dp[x]
    for c in Cwords:
        U = Qx ^ c
        if g.posw(U) + g.posw(U ^ D) <= THR:
            return False
    return True


def sample(k):
    """random columns 0..k-1: q_ij symmetric for i,j<k, q_ij free for i<k<=j (stored in col_i bit j and mirrored)."""
    q = [[0] * 8 for _ in range(8)]
    for i in range(k):
        for j in range(i, 8):
            v = rng.getrandbits(1)
            q[i][j] = v
            if j < k:
                q[j][i] = v
    return [sum(q[i][j] << j for j in range(8)) for i in range(k)] + [0] * (8 - k)


for k in (1, 2, 3, 4):
    total = 2 ** (k * (k + 1) // 2 + k * (8 - k))
    ns = 3000
    ok = 0
    for _ in range(ns):
        cols = sample(k)
        if all(ok_x(x, cols) for x in range(1, 1 << k)):
            ok += 1
    print("level", k, "MC fraction", ok / ns, "-> est nodes", ok / ns * total, flush=True)
