"""ctrl_fpf6.py SEED NTRIALS: positive control of qdfs/invtw on n=6 (AME(6,4), THR=3), Omega = 3 swaps, random K = (C,C)."""
import sys, random, time
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))  # [release] was the absolute work/pairs path
from invtw import *
import ame

n = 6
om = Omega(n, [(0, 3), (1, 4), (2, 5)], {})
st = Setup(om)
print("dim ker N", len(st.Vom))
rng = random.Random(int(sys.argv[1]))
ntr = int(sys.argv[2])
Vb = []  # basis of V^Omega: (w,w)
for p in range(3):
    for a in range(4):
        Vb.append((1 << a) << (4 * p) | ((1 << a) << (4 * (p + 3))))


def randvec():
    v = 0
    for b in Vb:
        if rng.getrandbits(1):
            v ^= b
    return v


found = 0
for t in range(ntr):
    # random K: isotropic, 6-dim, all nonzero weight >= 4 (THR=3)
    for attempt in range(1000):
        K = []
        span = [0]
        ok = True
        while len(K) < 6:
            for tries in range(200):
                v = randvec()
                if v in span:
                    continue
                if any(st.beta(v, k) for k in K):
                    continue
                if st.beta(v, v):
                    continue
                if all(wt(v ^ s, n) >= 4 for s in span):
                    break
            else:
                ok = False
                break
            K.append(v)
            span = span + [v ^ s for s in span]
        if ok:
            break
    t0 = time.time()
    r = run_qdfs(st, K, 3)
    if r is None:
        print("trial", t, "invalid K")
        continue
    nsol, sols, info = r
    print("trial", t, "nsol", nsol, "(%.1fs)" % (time.time() - t0), flush=True)
    for gens in sols:
        rows = to_rows(gens, n)
        b, tot = ame.bad_cuts(rows, 2, 2, n)
        print("   SOL bad_cuts", b, "of", tot, flush=True)
        found += 1
    if found >= 3:
        break
