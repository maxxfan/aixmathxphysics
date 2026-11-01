"""Run frames3 with JMASKS=124 (all conditions except the four with J={0,1,2}) on the V_0 classes, sequentially, nice -n 10.

RELEASE NOTE (release/frames): adapted from work/classify/run_m124_all.py, which had the binary (work/classify/frames3_s8),
the class file and the log hard-wired.  They are now arguments, and an optional list of class indices allows run.sh to run the
classes in parallel (one log per class; run.sh's canonicaliser recomputes the SUMMARY totals from the per-class lines).
usage: run_m124_all.py FRAMES3_BINARY V0FILE LOG [CLASS ...]      (default classes: 0..9)
"""
import os, re, subprocess, sys

exe, v0file, log = sys.argv[1], sys.argv[2], sys.argv[3]
classes = [int(c) for c in sys.argv[4:]] or list(range(10))
env = dict(os.environ)
env["JMASKS"] = "124"
open(log, "w").close()
tot1 = tot2 = tot3 = 0
for c in classes:
    cmd = ["nice", "-n", "10", exe, "4", "5", v0file, str(c), "-1", "-"]
    with open(log, "a") as f:
        f.write("# " + " ".join(cmd) + "  env=JMASKS=124\n")
        f.flush()
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, text=True).stdout
        f.write(out)
    m = re.search(r"nodes (\d+) (\d+) (\d+) (\d+)", out)
    if m:
        tot1 += int(m.group(2)); tot2 += int(m.group(3)); tot3 += int(m.group(4))
with open(log, "a") as f:
    f.write(f"# SUMMARY sum nodes[1]={tot1} sum nodes[2]={tot2} sum nodes[3]={tot3}\n")
