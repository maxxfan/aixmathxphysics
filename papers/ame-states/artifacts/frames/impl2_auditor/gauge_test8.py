"""Gauge-invariance test of vsearch8: transform each V_0 class rep by random local Sp(4,2)^4 and a random permutation of the 4 W-parties;
n_V1 and n_V2 (and n_SOL) must be IDENTICAL to the untransformed rep.

RELEASE NOTE (release/frames): adapted from work/auditor/gauge_test8.py (which compiled vsearch8.c to vsearch8g next to itself,
read v0reps8.txt and appended to gauge_test8.log in its own directory, and ran all ten reps with one seed).  Changes: binary,
representative file and log are arguments, and an optional rep range lets run.sh run the reps in parallel (one seed per call).
The transformation, the comparison (first '|'-field of the RESULT line) and the log format are unchanged.
usage: gauge_test8.py EXE REPFILE TRIALS_PER_REP SEED LOG [FIRST LAST]"""
import itertools, random, subprocess, sys


def sw(v):
    return ((v >> 2) & 0x3333) | ((v & 0x3333) << 2)


def om(u, v):
    return bin(u & sw(v)).count("1") & 1


def apply(cols, x):
    y = 0
    for i in range(4):
        if x >> i & 1:
            y ^= cols[i]
    return y


def is_inv(cols):
    return len({apply(cols, x) for x in range(16)}) == 16


def is_sp(cols):
    return all(om(apply(cols, 1 << i), apply(cols, 1 << j)) == om(1 << i, 1 << j) for i in range(4) for j in range(4))


def run(exe, b):
    r = subprocess.run([exe] + [str(x) for x in b] + ["0"], capture_output=True, text=True, timeout=300)
    return [l for l in r.stdout.splitlines() if l.startswith("RESULT")][0].split("|")[0].strip()


def main():
    exe, repfile = sys.argv[1], sys.argv[2]
    trials, seed, logname = int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    random.seed(seed)
    sp = [c for c in itertools.product(range(1, 16), repeat=4) if is_inv(c) and is_sp(c)]
    assert len(sp) == 720, len(sp)
    reps = [list(map(int, l.split())) for l in open(repfile) if l.strip()]
    first = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    last = int(sys.argv[7]) if len(sys.argv) > 7 else len(reps) - 1
    bad = 0
    out = open(logname, "a")
    for k in range(first, last + 1):
        rep = reps[k]
        b = rep[1:]
        base = run(exe, b)
        for t in range(trials):
            perm = list(range(4))
            random.shuffle(perm)
            gs = [random.choice(sp) for _ in range(4)]
            nb = []
            for v in b:
                w = 0
                for p in range(4):
                    nib = (v >> (4 * p)) & 15
                    w |= apply(gs[p], nib) << (4 * perm[p])
                nb.append(w)
            res = run(exe, nb)
            ok = res == base
            bad += not ok
            out.write("rep %d trial %d perm %s: %s | base %s | %s\n" % (k, t, perm, res, base, "OK" if ok else "MISMATCH"))
            out.flush()
        print("rep", k, "base", base)
    print("mismatches:", bad)
    out.write("mismatches: %d\n" % bad)


main()
