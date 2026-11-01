"""dump_check.py [CLASSLINE_D...] : independent check of the inv8 DFS bottom level. For one C (default: class 0 of inv8_classes.log),
dump the level-7 survivors (cols 0..6 chosen), complete both values of q_77, build S and compute exact min weight and the number
of nonzero words of weight <= 4 by brute force (numpy over all 2^16 words)."""
import sys, subprocess, collections
import numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # [release] was the absolute work/pairs path
import gen_inv8 as g

d = [int(x, 16) for x in (sys.argv[1:] if len(sys.argv) > 1 else "b001 6002 a004 5008 910 720 a40 580".split())]


def words_stats(rows):
    xs = np.zeros(1, dtype=np.uint32)
    zs = np.zeros(1, dtype=np.uint32)
    for r in rows:
        xv = sum(r[k] << k for k in range(16))
        zv = sum(r[16 + k] << k for k in range(16))
        xs = np.concatenate([xs, xs ^ np.uint32(xv)])
        zs = np.concatenate([zs, zs ^ np.uint32(zv)])
    w = np.zeros(len(xs), dtype=np.int32)
    for j in range(8):
        w += ((((xs >> (2 * j)) & 3) | ((zs >> (2 * j)) & 3)) != 0).astype(np.int32)
    w[0] = 99
    return int(w.min()), collections.Counter(w[w < 99].tolist())


pr = subprocess.run([g.BIN, "0", "1000000000000", "4", "1"], input=" ".join("%x" % v for v in d) + "\n", capture_output=True, text=True)
l7 = [l.split()[1:] for l in pr.stdout.splitlines() if l.startswith("L7")]
print("level-7 survivors:", len(l7), [l for l in pr.stdout.splitlines() if l.startswith("NSOL")])
res = collections.Counter()
best = None
for cols7 in l7:
    c = [int(x, 16) for x in cols7]
    for q77 in (0, 1):
        # column 7 = bit 7 of each earlier column (symmetric) plus q77 at bit 7
        c7 = sum(((c[j] >> 7) & 1) << j for j in range(7)) | (q77 << 7)
        cols = c + [c7]
        # symmetry check
        for i in range(8):
            for j in range(8):
                assert (cols[i] >> j & 1) == (cols[j] >> i & 1)
        rows = g.build_rows(d, cols)
        assert g.isotropic(rows) and g.rank2(rows) == 16
        mw, cnt = words_stats(rows)
        n_le4 = sum(v for k, v in cnt.items() if k <= 4)
        res[(mw, n_le4)] += 1
        if best is None or n_le4 < best[0]:
            best = (n_le4, mw, cols)
print("(min weight, #words of weight<=4) over the 2*%d completions:" % len(l7), dict(res))
print("best near miss (n words wt<=4, minwt, cols):", best)
