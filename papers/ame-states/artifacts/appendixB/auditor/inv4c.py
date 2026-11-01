"""Driver: random AME(4,4) graph-state codes C (inv4.py machinery), good table in numpy, DFS over all symmetric beta in C (inv4dfs).
usage: python3 work/auditor/inv4c.py SEED NSAMPLES MAXSECONDS   (run from the problem directory)   -> work/auditor/inv4c_SEED.log
Any solution is rebuilt as a 16x32 stabilizer, checked with ame.bad_cuts, and written to work/auditor/cands/inv4_hit_SEED.json (NOT submitted)."""
import os, sys, random, json, time, subprocess
os.environ.setdefault("OMP_NUM_THREADS", "2")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))   # [release] was ROOT (qlib package)
import numpy as np
import inv4
import ame   # [release] was: from qlib import ame
WORKDIR = os.environ.get("APPB_WORK", ".")   # [release] logs / .bin / cands go here instead of work/auditor
DFSBIN = os.path.join(os.environ["BUILD"], "inv4dfs") if os.environ.get("BUILD") else os.path.join(HERE, "inv4dfs")  # [release]

seed = int(sys.argv[1]); nsamp = int(sys.argv[2]); maxsec = float(sys.argv[3])
rng = random.Random(seed)
log = open(os.path.join(WORKDIR, f"inv4c_{seed}.log"), "a")   # [release]
t0 = time.time()
done = 0
classes = {}
while done < nsamp and time.time() - t0 < maxsec:
    Gam = [[0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1, 8):
            Gam[i][j] = Gam[j][i] = rng.randint(0, 1)
    g = inv4.graph_code(Gam)
    Call = inv4.code_elements(g)
    if inv4.party_weights4(Call[1:]).min() < 3:
        continue
    good = inv4.good_table(g, Call)
    fn = os.path.join(WORKDIR, f"good_{seed}.bin")   # [release]
    open(fn, "wb").write(good.astype(np.uint8).tobytes())
    out = subprocess.run([DFSBIN, fn, "1"], capture_output=True, text=True, timeout=600).stdout.split()
    nodes, nsol = int(out[0]), int(out[1])
    wd = np.bincount(inv4.party_weights4(Call[1:]), minlength=5).tolist()
    classes[tuple(wd)] = classes.get(tuple(wd), 0) + 1
    rec = {"sample": done, "wdist": wd, "frac_good": round(float(good[1:].mean()), 4), "nodes": nodes, "sols": nsol, "t": round(time.time() - t0, 1)}
    if nsol:
        beta = [int(x) for x in out[2:10]]
        rows = inv4.build_rows(g, Call, beta)
        bad, cuts = ame.bad_cuts(rows, 2, 2, 8)
        rec.update({"bad_cuts": bad, "beta": beta, "graph": Gam})
        if bad == 0:
            os.makedirs(os.path.join(WORKDIR, "cands"), exist_ok=True); ame.write_stabilizer(os.path.join(WORKDIR, "cands", f"inv4_hit_{seed}.json"), 8, 4, 2, rows, script=f"work/auditor/inv4c.py {seed} {nsamp} {maxsec}", note="involution-symmetric AME(8,4), AUDITOR search")
    log.write(json.dumps(rec) + "\n"); log.flush()
    done += 1
log.write(json.dumps({"DONE": done, "elapsed": round(time.time() - t0, 1), "weight_distribution_classes": {str(k): v for k, v in classes.items()}}) + "\n")
log.close()
