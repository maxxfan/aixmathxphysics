"""gen_inv8.py SEED NTRIALS [MAXSOL] : sample AME(4,4) Lagrangians C (graph form on 4 ququart positions), run inv8 (DFS over the
36-bit symmetric form q), and verify any solution S as a 16-qubit stabilizer state with ame.bad_cuts.
Writes candidate to work/pairs/cands/inv8_<seed>_<trial>.json when bad_cuts == 0."""
import sys, os, random, subprocess, time
# [release] paths: ../lib holds the copy of qlib/ame.py; inv8 binary from $BUILD (fallback: compile into cwd);
# candidate files go to $APPB_WORK/cands (default: ./cands) instead of work/pairs/cands.
HERE0 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE0, "..", "lib"))
import ame

HERE = os.environ.get("APPB_WORK", ".")
BIN = os.path.join(os.environ.get("BUILD", ""), "inv8")
if not os.environ.get("BUILD") or not os.path.exists(BIN):
    BIN = os.path.abspath("inv8")
    if not os.path.exists(BIN):
        subprocess.run([os.environ.get("CC", "cc"), "-O2", "-o", BIN, os.path.join(HERE0, "inv8.c")], check=True)


def posw(w):
    c = 0
    for p in range(4):
        if w & ((3 << (2 * p)) | (3 << (8 + 2 * p))):
            c += 1
    return c


def random_C(rng):
    """Return the 8 basis vectors d_k = (e_k | Gamma_k) of a Lagrangian C with all nonzero words of position-weight >= 3."""
    while True:
        G = [[0] * 8 for _ in range(8)]
        for p in range(4):
            for r in range(p + 1, 4):
                for a in range(2):
                    for b in range(2):
                        v = rng.getrandbits(1)
                        G[2 * p + a][2 * r + b] = v
                        G[2 * r + b][2 * p + a] = v
        d = []
        for k in range(8):
            z = sum(G[k][j] << j for j in range(8))
            d.append((1 << k) | (z << 8))
        ok = True
        for x in range(1, 256):
            w = 0
            for i in range(8):
                if x >> i & 1:
                    w ^= d[i]
            if posw(w) < 3:
                ok = False
                break
        if ok:
            return d


def build_rows(d, cols):
    """S generators as 16-qubit rows [x(16) | z(16)]; party p (p<4) = W position p first copy, party p+4 = second copy."""
    gens = []
    for i in range(8):
        ui = cols[i] << 8  # pure-Z vector with z = q_i
        gens.append((ui, ui ^ d[i]))
        gens.append((d[i], d[i]))
    rows = []
    for (a, b) in gens:
        x = [(a >> t) & 1 for t in range(8)] + [(b >> t) & 1 for t in range(8)]
        z = [(a >> (8 + t)) & 1 for t in range(8)] + [(b >> (8 + t)) & 1 for t in range(8)]
        rows.append(x + z)
    return rows


def isotropic(rows, N=16):
    for r in rows:
        for s in rows:
            t = 0
            for k in range(N):
                t ^= (r[k] & s[N + k]) ^ (r[N + k] & s[k])
            if t:
                return False
    return True


def rank2(rows):
    rs = [int("".join(map(str, r)), 2) for r in rows]
    rk = 0
    basis = []
    for v in rs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            rk += 1
    return rk


if __name__ == "__main__":
    seed, ntr = int(sys.argv[1]), int(sys.argv[2])
    maxsol = sys.argv[3] if len(sys.argv) > 3 else "5"
    maxcount = sys.argv[4] if len(sys.argv) > 4 else "1000000000000"
    rng = random.Random(seed)
    for t in range(ntr):
        d = random_C(rng)
        t0 = time.time()
        inp = " ".join("%x" % v for v in d) + "\n"
        pr = subprocess.run([BIN, maxsol, maxcount], input=inp, capture_output=True, text=True)
        lines = pr.stdout.strip().splitlines()
        summary = [l for l in lines if l.startswith("NSOL")]
        print("trial", t, "C=", " ".join("%x" % v for v in d), "::", summary, "(%.1fs)" % (time.time() - t0), flush=True)
        for l in lines:
            if l.startswith("SOL"):
                cols = [int(x, 16) for x in l.split()[1:]]
                rows = build_rows(d, cols)
                iso = isotropic(rows)
                rk = rank2(rows)
                b, tot = ame.bad_cuts(rows, 2, 2, 8) if (iso and rk == 16) else (-1, -1)
                print("   SOL cols", l.split()[1:], "isotropic", iso, "rank", rk, "bad_cuts", b, "of", tot, flush=True)
                if b == 0:
                    os.makedirs(os.path.join(HERE, "cands"), exist_ok=True)  # [release]
                    out = HERE + "/cands/inv8_%d_%d.json" % (seed, t)
                    script = "work/pairs/gen_inv8.py %d %d" % (seed, ntr)
                    ame.write_stabilizer(out, 8, 4, 2, rows, script=script, note="free-involution F2 stabilizer AME(8,4): fold C + symmetric q")
                    print("   WROTE", out, flush=True)
                    P = ame.project_last_party(rows, 2, 2, 8)
                    b2, _ = ame.bad_cuts(P, 2, 2, 7)
                    print("   projected bad_cuts", b2, flush=True)
                    ame.write_stabilizer(out.replace(".json", "_proj.json"), 7, 4, 2, P, script=script, note="projection of AME(8,4)")
                    sys.exit(0)
