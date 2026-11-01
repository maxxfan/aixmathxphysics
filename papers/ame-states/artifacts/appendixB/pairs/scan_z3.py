"""scan_z3.py : sector sigma = (0 1 2)(3 4 5) (fixed parties 6,7, no local twists) of F_2-linear stabilizer AME(8,4).
S = S^sigma (+) S_nf, S^sigma = lift of an AME(4,4) state on positions (A,B,f1,f2) (414 reps mod GL(2,2)^4, z3_reps.txt),
S_nf = graph of an F_4-anti-isometry phi: V_A -> V_B (enumerated by DFS; expect |GU(4,2)| = 77760).
Writes z3_in.txt, compiles + runs scan_z3.c, log work/pairs/z3_scan.log."""
import os, sys, subprocess, time
# [release] was HERE = <absolute path of work/pairs>; generated files (z3_in.txt, z3_scan.log) now go to the working directory
# ($APPB_WORK, default "."), the AME(4,4) reps file is $Z3_REPS (default ./z3_reps.txt, written by classes_z3.py) and the C
# helper is $BUILD/scan_z3 (compiled by run.sh) or, if $BUILD is unset, compiled here from the copy next to this script.
SRC = os.path.dirname(os.path.abspath(__file__))
HERE = os.environ.get("APPB_WORK", ".")
REPS = os.environ.get("Z3_REPS", os.path.join(HERE, "z3_reps.txt"))
LOG = open(HERE + "/z3_scan.log", "a")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def place(j, t):
    return 1 << (2 * j + t) if t < 2 else 1 << (16 + 2 * j + (t - 2))


PI = {0: 1, 1: 2, 2: 0, 3: 4, 4: 5, 5: 3, 6: 6, 7: 7}


def sigma(w):
    out = 0
    for q in range(16):
        j, s = q // 2, q % 2
        q2 = 2 * PI[j] + s
        if w >> q & 1:
            out |= 1 << q2
        if w >> (16 + q) & 1:
            out |= 1 << (16 + q2)
    return out


def form(a, b):
    return bin(((a & 0xFFFF) & (b >> 16)) ^ ((a >> 16) & (b & 0xFFFF))).count("1") & 1  # [release] was int.bit_count() (3.10+)


v = [place(0, t) ^ place(2, t) for t in range(4)]
sv = [sigma(x) for x in v]
bgens = []
for t in range(4):
    b = place(3, t) ^ place(5, t)
    bgens += [b, sigma(b)]
VB = [0]
for g in bgens:
    VB = VB + [x ^ g for x in VB]
assert len(set(VB)) == 256
sVB = {w: sigma(w) for w in VB}

res = []


def dfs(t, W):
    if t == 4:
        res.append(tuple(W))
        return
    for w in VB:
        sw = sVB[w]
        if form(w, sw) != form(v[t], sv[t]):
            continue
        ok = True
        for s in range(t):
            if (form(W[s], w) != form(v[s], v[t]) or form(W[s], sw) != form(v[s], sv[t]) or
                    form(sVB[W[s]], w) != form(sv[s], v[t])):
                ok = False
                break
        if ok:
            dfs(t + 1, W + [w])


t0 = time.time()
dfs(0, [])
P("anti-isometries phi found:", len(res), "(expect 77760)", "(%.1fs)" % (time.time() - t0))

# lifts of C_1 reps
# [release] next line: was open(HERE + "/z3_reps.txt")
reps = [[int(x, 16) for x in line.split()] for line in open(REPS) if line.strip()]
PARTIES = {0: (0, 1, 2), 1: (3, 4, 5), 2: (6,), 3: (7,)}


def lift(wv):
    out = 0
    for p in range(4):
        for s in range(2):
            xb = wv >> (2 * p + s) & 1
            zb = wv >> (8 + 2 * p + s) & 1
            for j in PARTIES[p]:
                if xb:
                    out |= 1 << (2 * j + s)
                if zb:
                    out |= 1 << (16 + 2 * j + s)
    return out


lifted = [[lift(x) for x in d] for d in reps]
# sanity: lifted gens isotropic
for g in lifted[:20]:
    for a in g:
        for b in g:
            assert form(a, b) == 0
with open(HERE + "/z3_in.txt", "w") as fo:
    fo.write("%d %d\n" % (len(lifted), len(res)))
    for g in lifted:
        fo.write(" ".join("%x" % x for x in g) + "\n")
    for W in res:
        gens = []
        for t in range(4):
            gens += [v[t] ^ W[t], sv[t] ^ sVB[W[t]]]
        fo.write(" ".join("%x" % x for x in gens) + "\n")
P("wrote z3_in.txt: R=%d M=%d" % (len(lifted), len(res)))
# isotropy sanity of S_nf gens
W = res[0]
gens = []
for t in range(4):
    gens += [v[t] ^ W[t], sv[t] ^ sVB[W[t]]]
assert all(form(a, b) == 0 for a in gens for b in gens), "S_nf not isotropic"
# [release] was: compile HERE/scan_z3.c into HERE and run it
if os.environ.get("BUILD") and os.path.exists(os.path.join(os.environ["BUILD"], "scan_z3")):
    EXE = os.path.join(os.environ["BUILD"], "scan_z3")
else:
    EXE = os.path.join(HERE, "scan_z3")
    subprocess.run(["cc", "-O2", "-o", EXE, os.path.join(SRC, "scan_z3.c")], check=True)
pr = subprocess.Popen([EXE, HERE + "/z3_in.txt"], stdout=subprocess.PIPE, text=True)
for line in pr.stdout:
    P(line.rstrip())
pr.wait()
