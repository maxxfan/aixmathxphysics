"""Run vsearch8 over the V_0 class representatives (v0reps8.txt or v0reps8_sp.txt).

RELEASE NOTE (release/frames): adapted from work/auditor/run_v8.py, which compiled vsearch8.c next to itself and read/wrote
fixed files in its own directory.  Changes: (1) the binary, the representative file and the log are arguments (run.sh compiles
the binary into build/); (2) the log line keeps the original format 'rep K orbit_size N V0 [...] mode M: <RESULT line> [<stderr>] T s'
but contains only the RESULT line of stdout: the current vsearch8.c also prints, in mode 0, up to 30 diagnostic lines starting
with 'SOL' for the first near-misses (V_0..V_2 and the V_3 pair conditions pass, a V_3 single-element condition fails; these are
NOT solutions; genuine solutions are counted by n_SOL and would be printed as 'SOL' lines as well).  All 'SOL' lines are written
to LOG.SOLlines_K.txt instead (the original kept the 30 near-misses of rep 0 in work/auditor/v8c_solutions.txt via run_v8n.py).
usage: run_v8.py EXE REPFILE FIRST LAST TIMEOUT_PER_REP MODE LOG
"""
import subprocess, sys, time

exe, repfile = sys.argv[1], sys.argv[2]
first, last, tmo = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
mode, logname = sys.argv[6], sys.argv[7]
reps = [list(map(int, l.split())) for l in open(repfile) if l.strip()]
log = open(logname, "a")
for k in range(first, min(last, len(reps) - 1) + 1):
    cnt, b = reps[k][0], reps[k][1:]
    t0 = time.time()
    try:
        r = subprocess.run([exe] + [str(x) for x in b] + [mode], capture_output=True, text=True, timeout=tmo)
        lines = r.stdout.strip().splitlines()
        out = " | ".join(l for l in lines if l.startswith("RESULT"))
        near = [l for l in lines if l.startswith("SOL")]
        if near:
            with open("%s.SOLlines_%d.txt" % (logname, k), "w") as fh:
                fh.write("\n".join(near) + "\n")
        err = r.stderr.strip()
        if r.returncode != 0:
            out += " EXIT %d" % r.returncode
    except subprocess.TimeoutExpired:
        out = "TIMEOUT after %d s" % tmo
        err = ""
    log.write("rep %d orbit_size %d V0 %s mode %s: %s [%s] %.1fs\n" % (k, cnt, b, mode, out, err, time.time() - t0))
    log.flush()
