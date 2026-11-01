"""Relaxed-threshold positive control for vsearch8 (mode 2: single wt >= 3, pair wt >= 2, V_3 likewise = min party weight 4 style).

RELEASE NOTE (release/frames): adapted from work/auditor/run_v8c.py.  The original compiled vsearch8.c next to itself, read
v0reps8.txt from its own directory and held the complete stdout in memory (8.2 million 'SOL' lines, about 0.6 GB, for MAXV1 = 400).
This version takes binary, representative file and output paths as arguments and streams stdout line by line (constant memory),
keeping the first 50 'SOL' lines (the solutions vbrute8.py checks) and counting the rest.  Log line format unchanged.
usage: run_v8c.py EXE REPFILE REP MAXV1 TIMEOUT LOG SOLFILE
"""
import subprocess, sys, time

exe, repfile = sys.argv[1], sys.argv[2]
rep, maxv1, tmo = int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
logname, solfile = sys.argv[6], sys.argv[7]
reps = [list(map(int, l.split())) for l in open(repfile) if l.strip()]
b = reps[rep][1:]
t0 = time.time()
res, sols, nsol = [], [], 0
p = subprocess.Popen([exe] + [str(x) for x in b] + ["2", maxv1], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
for line in p.stdout:
    if line.startswith("SOL"):
        nsol += 1
        if len(sols) < 50:
            sols.append(line.rstrip("\n"))
    elif line.startswith("RESULT"):
        res.append(line.strip())
    if time.time() - t0 > tmo:
        p.kill()
        break
rc = p.wait()
if time.time() - t0 > tmo:
    msg = "rep %d relaxed mode2: TIMEOUT %d s" % (rep, tmo)
else:
    msg = "rep %d relaxed mode2 maxV1 %s: %s ; #SOL lines %d ; %.1fs" % (rep, maxv1, res, nsol, time.time() - t0)
    open(solfile, "w").write("\n".join(sols) + "\n")
open(logname, "a").write(msg + "\n")
print(msg)
