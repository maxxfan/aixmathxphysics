"""inv8_thr3_ctrl.py [MAXSOL] [d0 .. d7]

[release] Positive control ADDED FOR THE RELEASE (not part of the original swarm run; uses pairs' own code only).
pairs' inv8 DFS (free-involution sector t = 4, fold C + symmetric 8x8 form q) has a built-in relaxed threshold:
argv[3] = THR, a word is bad iff its folded weight is <= THR (THR = 4: AME(8,4); THR = 3: distance-4 control).
This script runs inv8 with THR = 3 on the class-0 representative C of classes_inv8.py (default basis below), capped at
MAXSOL solutions (the DFS order is deterministic, so the capped count and the printed solutions are reproducible), and
rebuilds every printed solution q as a 16-qubit stabilizer with gen_inv8.build_rows, checking it directly:
isotropic, rank 16, minimum party weight over all 2^16 nonzero words = 4 (so not AME), ame.bad_cuts > 0.
The same C with THR = 4 must give NSOL 0 (classes_inv8.py)."""
import os, sys, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pairs"))   # [release] lives in extra/, uses pairs/gen_inv8.py
import gen_inv8 as g
import ame

MAXSOL = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
d = [int(x, 16) for x in (sys.argv[2:10] if len(sys.argv) > 9 else "b001 6002 a004 5008 910 720 a40 580".split())]
inp = " ".join("%x" % v for v in d) + "\n"


def run(thr, maxsol):
    pr = subprocess.run([g.BIN, "3", str(maxsol), str(thr)], input=inp, capture_output=True, text=True, check=True)
    lines = pr.stdout.strip().splitlines()
    return [l for l in lines if l.startswith("NSOL")][0], [l for l in lines if l.startswith("SOL")]


def min_party_weight(rows):
    """direct enumeration of all 2^16 words; party p = qubits 2p, 2p+1 of the 16-qubit rows [x16 | z16]."""
    gens = []
    for r in rows:
        v = 0
        # pack: bits 4p..4p+3 = (x_2p, x_2p+1, z_2p, z_2p+1)
        for p in range(8):
            nb = r[2 * p] | (r[2 * p + 1] << 1) | (r[16 + 2 * p] << 2) | (r[16 + 2 * p + 1] << 3)
            v |= nb << (4 * p)
        gens.append(v)
    span = [0]
    for v in gens:
        span += [v ^ s for s in span]
    best, nbad = 99, 0
    for w in span[1:]:
        c = sum(1 for p in range(8) if (w >> (4 * p)) & 15)
        best = min(best, c)
        nbad += c <= 4
    return best, nbad, len(set(span))


line4, _ = run(4, 10 ** 12)
print("THR=4 (AME search) on C:", " ".join(line4.split()[:2]))
line3, sols = run(3, MAXSOL)
print("THR=3 (distance-4 control) on C, cap %d:" % MAXSOL, line3)
for s in sols:
    cols = [int(x, 16) for x in s.split()[1:]]
    rows = g.build_rows(d, cols)
    iso = g.isotropic(rows)
    rk = g.rank2(rows)
    mw, nbad, nwords = min_party_weight(rows)
    b, tot = ame.bad_cuts(rows, 2, 2, 8)
    print("  solution q =", " ".join(s.split()[1:]), ": isotropic", iso, "rank", rk, "distinct words", nwords,
          "min party weight", mw, "words of weight <= 4:", nbad, "bad_cuts", b, "of", tot)
