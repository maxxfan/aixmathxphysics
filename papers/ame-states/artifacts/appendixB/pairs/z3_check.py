"""z3_check.py : independent Python check of the Z_3 (3,3,1,1) sector: minimum party weight of S_nf (graph of phi) over ALL 77760 phi."""
import collections
import os
HERE = os.environ.get("APPB_WORK", ".")  # [release] was the absolute path of work/pairs; reads z3_in.txt written by scan_z3.py
lines = open(HERE + "/z3_in.txt").read().splitlines()
R, M = map(int, lines[0].split())
phis = [[int(x, 16) for x in l.split()] for l in lines[1 + R:1 + R + M]]
pm = [(3 << (2 * j)) | (3 << (16 + 2 * j)) for j in range(8)]
hist = collections.Counter()
nlow = collections.Counter()
for g in phis:
    w = 0
    mn = 99
    cnt4 = 0
    for i in range(1, 256):
        t = (i & -i).bit_length() - 1
        w ^= g[t]
        k = sum(1 for m in pm if w & m)
        if k < mn:
            mn = k
        if k <= 4:
            cnt4 += 1
    hist[mn] += 1
    nlow[cnt4] += 1
print("M =", M, "min party weight of S_nf over all phi:", dict(hist))
print("number of words of weight <= 4 in S_nf (value: #phi):", dict(sorted(nlow.items())))
