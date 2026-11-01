"""canon.py (release/frames, written for this release; standard library only).

Normalises the raw output of the two implementations to one canonical text per check: only the counts (class index, orbit size,
node counts, numbers of solutions, ...) are kept; run times, absolute paths and command lines are dropped.  The same functions are
applied to the ORIGINAL logs (tools/make_expected.py -> expected/) and to the fresh replay (run.sh -> out/canon/), so a check passes
iff the replay reproduces the recorded numbers exactly.

usage: canon.py KIND [options] FILE ...   (prints the canonical text on stdout)
  frames [--hist] [--count M] [--select "M T"] LOG...   impl1 frames/frames3 class lines
  v8 LOG...                                             impl2 run_v8 log lines (core counts)
  v8diag LOG...                                         impl2 V_3 diagnostic counters (current vsearch8.c, or vsearch8h run_v8h.log)
  gauge LOG...                                          impl2 gauge_test8 log
  v8relaxed LOG                                         impl2 run_v8c log
  reps FILE                                             a class-representative file (impl1 v0_m*.txt or impl2 v0reps8*.txt)
  sols FILE                                             a solution file (line count and SHA-256 of the whitespace-normalised lines)
  vcount LOG | classcount LOG | enumtest LOG | vclasses LOG | v0reps LOG    stdout of the auditor's / classify's count scripts
  vbrute [--minonly] LOG                                stdout of impl2 vbrute8.py
  text FILE...                                          non-empty lines with trailing blanks removed (framecheck / crosswalk output)
"""
import hashlib
import re
import sys
from fractions import Fraction

FR = re.compile(r"^class (\d+) \(tuple-orbit (\d+)\): nodes ([\d ]+?)\s+nsol=(\d+)\s+time=[\d.]+s(?:\s+free-dim histogram:(.*))?$")


def parse_frames(lines, select=None):
    """class lines of frames/frames3 output; with select='M T', only the sections whose '# ...' header runs 'frames* M T ...'."""
    out = {}
    active = select is None
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#"):
            if select is not None:
                m = re.search(r"/?frames3?(?:_s8)?\s+(\d+)\s+(\d+)\s", line + " ")
                if m:
                    active = ("%s %s" % (m.group(1), m.group(2))) == select
            continue
        m = FR.match(line.strip())
        if not m or not active:
            continue
        c = int(m.group(1))
        nodes = [int(x) for x in m.group(3).split()]
        hist = None
        if m.group(5) is not None:
            hist = " ".join(re.sub(r"\s+", " ", h).replace(" }", "}").replace("{ ", "{") for h in re.findall(r"L\d+\{[^}]*\}", m.group(5)))
        rec = (int(m.group(2)), nodes, int(m.group(4)), hist)
        if c in out and out[c] != rec:
            raise SystemExit("conflicting records for class %d: %s vs %s" % (c, out[c], rec))
        out[c] = rec
    return out


def canon_frames(files, hist=False, count=None, select=None):
    recs = {}
    for f in files:
        for c, r in parse_frames(open(f), select).items():
            if c in recs and recs[c] != r:
                raise SystemExit("conflicting records for class %d" % c)
            recs[c] = r
    lines = []
    tot = None
    nsol_tot = 0
    for c in sorted(recs):
        orb, nodes, nsol, h = recs[c]
        s = "class=%d orbit=%d nodes=%s nsol=%d" % (c, orb, ",".join(map(str, nodes)), nsol)
        if hist:
            s += " hist=%s" % (h if h is not None else "MISSING")
        lines.append(s)
        tot = nodes[:] if tot is None else [a + b for a, b in zip(tot, nodes)]
        nsol_tot += nsol
    lines.append("TOTAL classes=%d nodes=%s nsol=%d" % (len(recs), ",".join(map(str, tot or [])), nsol_tot))
    if count is not None:
        m = count
        # PROOF.md 8.1/8.2: #(Q-tuples with Q_0 = J) in class c = orbit_c / 28, each gives 720^(m-1) graph-form V_0 with these forms,
        # each V_0 has nsol_c ordered completions; labelled states N = 720^m * #ordered frames.
        w = sum(Fraction(recs[c][0], 28) * recs[c][2] for c in recs)
        assert w.denominator == 1
        frames_ = int(w) * 720 ** (m - 1)
        lines.append("DERIVED sum_c(orbit_c/28*nsol_c)=%d ordered_frames=%d N_labelled_stabilizer_states=%d" % (int(w), frames_, frames_ * 720 ** m))
    return lines


V8 = re.compile(r"^rep (\d+) orbit_size (\d+) V0 \[([\d, ]+)\] mode (\d+): (.*)$")


def res_fields(s):
    m = re.search(r"RESULT n_V1=(\d+) n_V1_with_V2=(\d+) n_V2=(\d+) n_SOL=(\d+)", s)
    return None if m is None else tuple(int(x) for x in m.groups())


def canon_v8(files):
    recs = {}
    for f in files:
        for line in open(f):
            m = V8.match(line.strip())
            if not m:
                continue
            k, orb, v0, mode, rest = int(m.group(1)), int(m.group(2)), m.group(3).replace(" ", ""), m.group(4), m.group(5)
            r = res_fields(rest)
            c1 = re.search(r"cand1 = (\d+)", rest)
            if r is None:
                rec = (orb, v0, mode, "NO RESULT: " + rest[:80])
            else:
                rec = (orb, v0, mode, "cand1=%s n_V1=%d n_V1_with_V2=%d n_V2=%d n_SOL=%d" % (((c1.group(1) if c1 else "?"),) + r))
            if k in recs and recs[k] != rec:
                raise SystemExit("conflicting records for rep %d" % k)
            recs[k] = rec
    lines = []
    tv1 = tv2 = ts = to = 0
    for k in sorted(recs):
        orb, v0, mode, s = recs[k]
        lines.append("rep=%d orbit=%d V0=%s mode=%s %s" % (k, orb, v0, mode, s))
        to += orb
        r = re.findall(r"n_V1=(\d+) n_V1_with_V2=\d+ n_V2=(\d+) n_SOL=(\d+)", s)
        if r:
            tv1 += int(r[0][0]); tv2 += int(r[0][1]); ts += int(r[0][2])
    lines.append("TOTAL reps=%d orbit_sum=%d n_V1=%d n_V2=%d n_SOL=%d" % (len(recs), to, tv1, tv2, ts))
    return lines


def canon_v8diag(files):
    recs = {}
    D = re.compile(r"V3 singles ok (\d+), V3 pairs ok (\d+), singles-ok-but-pairs-fail (\d+), pairs-ok-but-singles-fail (\d+)")
    for f in files:
        for line in open(f):
            m = re.match(r"^rep (\d+)[ :]", line)
            d = D.search(line)
            r = res_fields(line)
            if not (m and d and r):
                continue
            k = int(m.group(1))
            recs[k] = "rep=%d n_V2=%d n_SOL=%d V3_singles_ok=%s V3_pairs_ok=%s singles_ok_pairs_fail=%s pairs_ok_singles_fail=%s" % (
                (k, r[2], r[3]) + d.groups())
    return [recs[k] for k in sorted(recs)] + ["TOTAL reps=%d" % len(recs)]


def canon_gauge(files):
    per = {}
    mism = []
    G = re.compile(r"^rep (\d+) trial (\d+) perm \[[\d, ]+\]: (RESULT[^|]*)\| base (RESULT[^|]*)\| (OK|MISMATCH)")
    for f in files:
        for line in open(f):
            m = G.match(line.strip())
            if m:
                k = int(m.group(1))
                base = res_fields(m.group(4))
                got = res_fields(m.group(3))
                d = per.setdefault(k, {"base": base, "trials": 0, "ok": 0})
                if d["base"] != base:
                    raise SystemExit("conflicting base for rep %d" % k)
                d["trials"] += 1
                d["ok"] += (m.group(5) == "OK" and got == base)
            m2 = re.match(r"^mismatches: (\d+)", line.strip())
            if m2:
                mism.append(int(m2.group(1)))
    lines = []
    for k in sorted(per):
        d = per[k]
        lines.append("rep=%d base=n_V1:%d,n_V1_with_V2:%d,n_V2:%d,n_SOL:%d trials=%d all_trials_identical=%s" % (
            (k,) + d["base"] + (d["trials"], "yes" if d["ok"] == d["trials"] else "NO")))
    lines.append("TOTAL reps=%d trials=%d mismatches=%d" % (len(per), sum(d["trials"] for d in per.values()), sum(mism)))
    return lines


def canon_v8relaxed(files):
    """Counts of the relaxed control.  The raw number of 'SOL' output lines ('#SOL lines' in the log) is deliberately NOT part of
    the canonical text: the current vsearch8.c (edited after the original run) also prints its first 30 near-misses with the
    prefix 'SOL', so the replay prints n_SOL + 30 such lines where the original binary printed n_SOL (see README.md, sec. 9)."""
    lines = []
    for f in files:
        for line in open(f):
            m = re.match(r"^rep (\d+) relaxed mode2 maxV1 (\d+): \['RESULT \(capped at (\d+) V_1\) n_V1=(\d+) n_V1_with_V2=(\d+) n_V2=(\d+) n_SOL=(\d+)'\] ; #SOL lines (\d+)", line.strip())
            if m:
                lines.append("rep=%s mode=2 maxV1=%s capped_at=%s n_V1=%s n_V1_with_V2=%s n_V2=%s n_SOL=%s" % m.groups()[:7])
    return lines


def canon_reps(files):
    lines = []
    for f in files:
        rows = [l.split() for l in open(f) if l.strip()]
        for r in rows:
            lines.append(" ".join(r))
        sizes = sorted(int(r[0]) for r in rows)
        lines.append("TOTAL classes=%d orbit_sum=%d sizes=%s" % (len(rows), sum(sizes), sizes))
    return lines


def canon_sols(files):
    lines = []
    for f in files:
        rows = [" ".join(l.split()) for l in open(f) if l.strip()]
        h = hashlib.sha256(("\n".join(rows) + "\n").encode()).hexdigest()
        lines.append("lines=%d sha256=%s" % (len(rows), h))
    return lines


def canon_vcount(files):
    txt = "".join(open(f).read() for f in files)
    g = lambda pat: re.search(pat, txt)
    lines = []
    m = g(r"n=4: valid Q forms (\d+) of 28; ordered \(V0,V1\) pairs (\d+); N\(4\) = 720\^2 \* T4 = (\d+)")
    lines.append("n4 admissible_forms=%s/28 ordered_pairs=%s N4=%s" % m.groups() if m else "n4 MISSING")
    m = g(r"Sp\(4,2\)-orbits on admissible pairs: (\d+)")
    lines.append("n6 sp_orbits=%s" % m.group(1) if m else "n6 sp_orbits MISSING")
    m = g(r"t per orbit rep \(valid V_1\) : (\[[\d, ]+\])")
    lines.append("n6 t_per_orbit_rep=%s" % m.group(1) if m else "n6 t MISSING")
    m = g(r"sum over pairs of t = (\d+) ; ordered triples T\(6\) = 720\^2 \* that = (\d+) ; N\(6\) labelled stabilizer AME\(6,4\) = (\d+)")
    lines.append("n6 sum_orb_t=%s ordered_triples=%s N6=%s" % m.groups() if m else "n6 MISSING")
    return lines


def canon_classcount(files):
    txt = "".join(open(f).read() for f in files)
    lines = []
    for h, tup, tot in re.findall(r"h=(\d+): tuples with nondegenerate sum: (\d+) \(of (\d+)\)", txt):
        m = re.search(r"h=%s: tuples.*?\n\s+orbits under GL\(4,2\) x S_%s: (\d+), sizes (\[[\d, ]+\]), total (\d+)" % (h, h), txt)
        lines.append("h=%s tuples=%s/%s orbits=%s sizes=%s total=%s" % ((h, tup, tot) + (m.groups() if m else ("?", "?", "?"))))
    return lines


def canon_enumtest(files):
    txt = "".join(open(f).read() for f in files)
    lines = ["U=%s subspaces=%s" % m for m in re.findall(r"U=(\d+): enumerated (\d+) subspaces", txt)]
    trials = re.findall(r"DFS (\d+) subspaces ; brute-force ordered bases (\d+) / 20160", txt)
    agree = all(int(b) == 20160 * int(a) for a, b in trials)
    lines.append("random_trials=%d dfs_equals_bruteforce_over_20160=%s" % (len(trials), "yes" if agree and trials else "NO"))
    return lines


def canon_vclasses(files):
    txt = "".join(open(f).read() for f in files)
    lines = []
    m = re.search(r"admissible triples: (\d+) of (\d+)", txt)
    lines.append("admissible_triples=%s/%s" % m.groups() if m else "MISSING")
    m = re.search(r"Sp\(4,2\) only : orbits (\d+) sizes", txt)
    lines.append("sp_only_orbits=%s" % m.group(1) if m else "MISSING")
    m = re.search(r"= Sp\(4,2\)\^4 x S_4 on V_0\) : orbits (\d+) sizes (\[[\d, ]+\])", txt)
    lines.append("full_gauge_orbits=%s sizes=%s" % m.groups() if m else "MISSING")
    return lines


def canon_v0reps(files):
    txt = "".join(open(f).read() for f in files)
    m = re.search(r"m=(\d+): (\d+) classes, (\d+) valid tuples", txt)
    return ["m=%s classes=%s valid_tuples=%s" % m.groups() if m else "MISSING"]


def canon_vbrute(files, minonly=False):
    txt = "".join(open(f).read() for f in files)
    m = re.search(r"checked (\d+) relaxed solutions; histogram of brute-force min party weight: (\{[^}]*\})", txt)
    lines = ["checked=%s min_party_weight_histogram=%s" % m.groups() if m else "MISSING"]
    if not minonly:
        m = re.search(r"histogram of \(min party weight, number of words of that weight\): (\{[^}]*\})", txt)
        lines.append("min_weight_and_multiplicity_histogram=%s" % m.group(1) if m else "MISSING")
    return lines


def main(argv):
    kind = argv[0]
    args = argv[1:]
    opts = {"hist": False, "count": None, "select": None}
    files = []
    minonly = False
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--hist":
            opts["hist"] = True
        elif a == "--count":
            opts["count"] = int(args[i + 1]); i += 1
        elif a == "--select":
            opts["select"] = args[i + 1]; i += 1
        elif a == "--minonly":
            minonly = True
        else:
            files.append(a)
        i += 1
    if kind == "frames":
        out = canon_frames(files, **opts)
    elif kind == "vbrute":
        out = canon_vbrute(files, minonly)
    elif kind == "text":
        out = [l.rstrip() for f in files for l in open(f) if l.strip()]
    else:
        out = {"v8": canon_v8, "v8diag": canon_v8diag, "gauge": canon_gauge, "v8relaxed": canon_v8relaxed, "reps": canon_reps,
               "sols": canon_sols, "vcount": canon_vcount, "classcount": canon_classcount, "enumtest": canon_enumtest,
               "vclasses": canon_vclasses, "v0reps": canon_v0reps}[kind](files)
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
