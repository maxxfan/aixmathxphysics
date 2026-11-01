"""run.py LOG [KEY=VAL ...] cmd... : run a command under nice -n 10 (leading KEY=VAL args become environment variables),
append stdout/stderr to LOG, return its exit status.  (binaries cannot be launched directly in don't-ask mode)

RELEASE NOTE (release/frames): copied from work/classify/run.py.  Only change: the '# ... env=' header line now also records
JMASKS and LOGL2 (the original recorded only DROPJ and PATTERN, so its logs show env={} for JMASKS runs)."""
import os
import subprocess
import sys

log = sys.argv[1]
args = sys.argv[2:]
env = dict(os.environ)
while args and "=" in args[0] and not args[0].startswith("/"):
    k, v = args.pop(0).split("=", 1)
    env[k] = v
cmd = ["nice", "-n", "10"] + args
with open(log, "a") as f:
    f.write("# " + " ".join(cmd) + "  env=" + str({k: env[k] for k in ("DROPJ", "PATTERN", "JMASKS", "LOGL2") if k in env}) + "\n")
    f.flush()
    rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    f.write(f"# exit {rc}\n")
sys.exit(rc)
