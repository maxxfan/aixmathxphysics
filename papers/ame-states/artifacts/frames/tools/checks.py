"""checks.py (release/frames, written for this release; standard library only): the list of checks.

Every check has a name, a one-line description, and two producers of the same canonical text:
  * 'expected': from the ORIGINAL records (original_logs/, copies of the project's logs and output files) or, where no log
    file exists, from numbers documented in PROOF.md / REVIEW.md (the source is named in SOURCES below);
  * 'replay'  : from the fresh run (out/logs, out/work) produced by run.sh.
A check passes iff the two texts are byte-identical (and the replay text has the SHA-256 listed in expected/HASHES.txt).

usage:
  checks.py expected DESTDIR          write expected/<name>.txt for every check (run once by the release author; run.sh re-runs it
                                      into out/expected_regen and compares, to show that expected/ is derived from original_logs/)
  checks.py replay OUTROOT DESTDIR    write the canonical replay texts (OUTROOT = run.sh's out/ directory)
  checks.py list                      print name, description and source of every check
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import canon  # noqa: E402

O1 = os.path.join(ROOT, "original_logs", "impl1_classify")
O2 = os.path.join(ROOT, "original_logs", "impl2_auditor")
PY = sys.executable


def tool(script, *args):
    r = subprocess.run([PY, os.path.join(HERE, script)] + list(args), capture_output=True, text=True)
    return [l.rstrip() for l in r.stdout.splitlines() if l.strip()] + ([] if r.returncode == 0 else ["EXIT %d %s" % (r.returncode, r.stderr.strip()[-200:])])


def need(*paths):
    miss = [p for p in paths if not os.path.exists(p)]
    return ["MISSING %s" % os.path.relpath(p, ROOT) for p in miss]


def guarded(paths, fn):
    m = need(*paths)
    return m if m else fn()


# ----------------------------------------------------------------------------------------------------------------- definitions
CHECKS = []  # (name, description, source of expected, expected producer, replay producer(out))


def check(name, desc, source):
    def deco(pair):
        CHECKS.append((name, desc, source, pair[0], pair[1]))
        return pair
    return deco


def L(out, *p):
    return os.path.join(out, "logs", *p)


def W(out, *p):
    return os.path.join(out, "work", *p)


# ---- impl1: V_0 classes
for m in (2, 3, 4):
    check("i1_classes_m%d" % m, "impl1 v0reps.py %d: V_0 gauge classes for n = %d (orbit size, psi_1..psi_%d)" % (m, 2 * m, m - 1),
          "original_logs/impl1_classify/v0_m%d.txt" % m)((
              (lambda m=m: canon.canon_reps([os.path.join(O1, "v0_m%d.txt" % m)])),
              (lambda out, m=m: guarded([W(out, "v0_m%d.txt" % m)], lambda: canon.canon_reps([W(out, "v0_m%d.txt" % m)])))))

check("i1_classes_summary", "impl1 v0reps.py stdout: number of classes and of admissible form tuples, n = 4, 6, 8",
      "PROOF.md sec. 5 (1 class / 336 tuples; 5 classes / 9856; 10 classes / 268800)")((
          lambda: ["m=2 classes=1 valid_tuples=336", "m=3 classes=5 valid_tuples=9856", "m=4 classes=10 valid_tuples=268800"],
          lambda out: [l for m in (2, 3, 4) for l in guarded([L(out, "i1_v0reps_m%d.log" % m)], lambda m=m: canon.canon_v0reps([L(out, "i1_v0reps_m%d.log" % m)]))]))

# ---- impl1: the search
check("i1_search_t5_frames3", "impl1 frames3.c, n = 8, t = 5 (all 68 conditions), all ten classes: nodes (1, #V_1, #V_2, #V_3), nsol, free-dim histograms",
      "original_logs/impl1_classify/final_t5.log")((
          lambda: canon.canon_frames([os.path.join(O1, "final_t5.log")], hist=True),
          lambda out: guarded([L(out, "i1_frames3_t5_c%d.log" % c) for c in range(10)],
                              lambda: canon.canon_frames([L(out, "i1_frames3_t5_c%d.log" % c) for c in range(10)], hist=True))))

check("i1_search_t5_frames", "impl1 frames.c (first version, no diagnostics), n = 8, t = 5, all ten classes: node counts",
      "original_logs/impl1_classify/ame84_{c0,c1,A,B,c9}.log")((
          lambda: canon.canon_frames([os.path.join(O1, f) for f in ("ame84_c0.log", "ame84_c1.log", "ame84_A.log", "ame84_B.log", "ame84_c9.log")]),
          lambda out: guarded([L(out, "i1_frames_t5_c%d.log" % c) for c in range(10)],
                              lambda: canon.canon_frames([L(out, "i1_frames_t5_c%d.log" % c) for c in range(10)]))))

check("i1_level2_m124", "impl1 frames3.c with JMASKS=124 (every condition except the four with J = {0,1,2}): level-2 survivors per class (total 24234)",
      "original_logs/impl1_classify/m124_all.log")((
          lambda: canon.canon_frames([os.path.join(O1, "m124_all.log")], hist=True),
          lambda out: guarded([L(out, "i1_frames3_m124_c%d.log" % c) for c in range(10)],
                              lambda: canon.canon_frames([L(out, "i1_frames3_m124_c%d.log" % c) for c in range(10)], hist=True))))

# ---- impl1: positive controls
for name, sel, cnt, desc in (("i1_ctrl_n4", "2 3", 2, "positive control n = 4 (frames.c m = 2, t = 3, AME(4,4)): completions and the derived number of labelled AME(4,4) stabilizers"),
                             ("i1_ctrl_n6", "3 4", 3, "positive control n = 6 (frames.c m = 3, t = 4, AME(6,4)): completions per class (112, 64, 288, 368, 2272) and N(6) = 54400 * 720^5"),
                             ("i1_ctrl_m4t4", "4 4", None, "positive control n = 8, t = 4 (min weight 4 instead of 5), class 9, stop after 20 frames")):
    check(name, desc, "original_logs/impl1_classify/cal.log")((
        (lambda sel=sel, cnt=cnt: canon.canon_frames([os.path.join(O1, "cal.log")], count=cnt, select=sel)),
        (lambda out, sel=sel, cnt=cnt: guarded([L(out, "i1_ctrl.log")], lambda: canon.canon_frames([L(out, "i1_ctrl.log")], count=cnt, select=sel)))))

check("i1_ctrl_relaxed", "positive control: frames3.c JMASKS=64479 (drop the (2,2) conditions with J = {0,2}, {1,3}), class 9, stop after 60 leaves",
      "original_logs/impl1_classify/relax_leaf.log")((
          lambda: canon.canon_frames([os.path.join(O1, "relax_leaf.log")], hist=True),
          lambda out: guarded([L(out, "i1_ctrl_relaxed.log")], lambda: canon.canon_frames([L(out, "i1_ctrl_relaxed.log")], hist=True))))

SOLS = (("n4", "sol_m2.txt", 2), ("n6", "sol_m3.txt", 3), ("m4t4", "sol_m4_t4.txt", 4), ("relaxed", "sol_relax_m4.txt", 4))
for tag, f, m in SOLS:
    check("i1_solfile_%s" % tag, "impl1 control %s: the list of frames found (exact file content, line count + SHA-256)" % tag,
          "original_logs/impl1_classify/%s" % f)((
              (lambda f=f: canon.canon_sols([os.path.join(O1, f)])),
              (lambda out, f=f: guarded([W(out, f)], lambda: canon.canon_sols([W(out, f)])))))
    check("i1_framecheck_%s" % tag, "tools/framecheck.py on the %s frames: orthogonality, nondegeneracy, brute-force party-weight enumerator%s" % (
        tag, ", (J,L) pattern vs exhaustive word check" if tag == "relaxed" else ""),
          "tools/framecheck.py run on original_logs/impl1_classify/%s (agrees with PROOF.md 8.2, 8.3, 8.6c)" % f)((
              (lambda f=f, m=m: tool("framecheck.py", str(m), os.path.join(O1, f))),
              (lambda out, f=f, m=m: guarded([W(out, f)], lambda: tool("framecheck.py", str(m), W(out, f))))))

# ---- impl2: V_0 classes
check("i2_classes_full", "impl2 vclasses8.py: the ten full-gauge V_0 representatives (orbit size in triples, 4 basis vectors)",
      "original_logs/impl2_auditor/v0reps8.txt")((
          lambda: canon.canon_reps([os.path.join(O2, "v0reps8.txt")]),
          lambda out: guarded([W(out, "v0reps8.txt")], lambda: canon.canon_reps([W(out, "v0reps8.txt")]))))
check("i2_classes_sp", "impl2 vclasses8.py: the 72 Sp(4,2)-only V_0 representatives (no party permutations used)",
      "original_logs/impl2_auditor/v0reps8_sp.txt")((
          lambda: canon.canon_reps([os.path.join(O2, "v0reps8_sp.txt")]),
          lambda out: guarded([W(out, "v0reps8_sp.txt")], lambda: canon.canon_reps([W(out, "v0reps8_sp.txt")]))))
check("i2_classes_summary", "impl2 vclasses8.py stdout: admissible triples, number of orbits and orbit sizes",
      "REVIEW.md (argument.md 5.12): 9600 of 28^3 admissible, 72 Sp-orbits, 10 classes of sizes 48 ... 2880")((
          lambda: ["admissible_triples=9600/21952", "sp_only_orbits=72", "full_gauge_orbits=10 sizes=[48, 72, 360, 480, 720, 720, 720, 1440, 2160, 2880]"],
          lambda out: guarded([L(out, "i2_vclasses8.log")], lambda: canon.canon_vclasses([L(out, "i2_vclasses8.log")]))))

# ---- impl2: the search
check("i2_search_full", "impl2 vsearch8.c mode 0, ten full-gauge representatives: n_V1, n_V1_with_V2, n_V2, n_SOL",
      "original_logs/impl2_auditor/run_v8_{0_0,1_4,5_9}.log")((
          lambda: canon.canon_v8([os.path.join(O2, f) for f in ("run_v8_0_0.log", "run_v8_1_4.log", "run_v8_5_9.log")]),
          lambda out: guarded([L(out, "i2_vsearch8_full_r%d.log" % k) for k in range(10)],
                              lambda: canon.canon_v8([L(out, "i2_vsearch8_full_r%d.log" % k) for k in range(10)]))))
check("i2_search_sp72", "impl2 vsearch8.c mode 0 on all 72 Sp(4,2)-only representatives",
      "original_logs/impl2_auditor/run_v8_v0reps8_sp_{0_35,36_71}.log")((
          lambda: canon.canon_v8([os.path.join(O2, f) for f in ("run_v8_v0reps8_sp_0_35.log", "run_v8_v0reps8_sp_36_71.log")]),
          lambda out: guarded([L(out, "i2_vsearch8_sp_r%d.log" % k) for k in range(72)],
                              lambda: canon.canon_v8([L(out, "i2_vsearch8_sp_r%d.log" % k) for k in range(72)]))))
check("i2_v3diag_full", "impl2 V_3 diagnostics (ten classes): V_3 never passes its single-element conditions, always passes its pair conditions",
      "original_logs/impl2_auditor/run_v8h.log (vsearch8h.c = vsearch8.c + histograms, same counters)")((
          lambda: canon.canon_v8diag([os.path.join(O2, "run_v8h.log")]),
          lambda out: guarded([L(out, "i2_vsearch8_full_r%d.log" % k) for k in range(10)],
                              lambda: canon.canon_v8diag([L(out, "i2_vsearch8_full_r%d.log" % k) for k in range(10)]))))
check("i2_gauge", "impl2 gauge_test8.py: 3 random Sp(4,2)^4 x S_4 transforms of each representative give identical counts",
      "original_logs/impl2_auditor/gauge_test8.log (the random transforms differ; only counts and agreement are compared)")((
          lambda: canon.canon_gauge([os.path.join(O2, "gauge_test8.log")]),
          lambda out: guarded([L(out, "i2_gauge_r%d.log" % k) for k in range(10)],
                              lambda: canon.canon_gauge([L(out, "i2_gauge_r%d.log" % k) for k in range(10)]))))

# ---- impl2: positive controls
check("i2_ctrl_relaxed", "impl2 positive control: vsearch8 mode 2 (min party weight 4 thresholds), rep 0, first 400 V_1: millions of solutions "
      "(counts n_V1, n_V1_with_V2, n_V2, n_SOL; the raw number of 'SOL' output lines is not compared, see README sec. 9)",
      "original_logs/impl2_auditor/run_v8c.log")((
          lambda: canon.canon_v8relaxed([os.path.join(O2, "run_v8c.log")]),
          lambda out: guarded([L(out, "i2_relaxed.log")], lambda: canon.canon_v8relaxed([L(out, "i2_relaxed.log")]))))
check("i2_vbrute_relaxed", "impl2 vbrute8.py: brute-force min party weight of the first 50 relaxed solutions",
      "REVIEW.md (argument.md 5.12): '50 of them brute-force verified with min party weight exactly 4'")((
          lambda: ["checked=50 min_party_weight_histogram={4: 50}"],
          lambda out: guarded([L(out, "i2_vbrute_relaxed.log")], lambda: canon.canon_vbrute([L(out, "i2_vbrute_relaxed.log")], minonly=True))))
check("i2_nearmiss_rep0", "impl2 mode 0, rep 0: the first 30 near-misses ('SOL' lines with n_SOL = 0), exact content",
      "original_logs/impl2_auditor/v8c_solutions.txt (written by run_v8n.py)")((
          lambda: canon.canon_sols([os.path.join(O2, "v8c_solutions.txt")]),
          lambda out: guarded([L(out, "i2_vsearch8_full_r0.log.SOLlines_0.txt")], lambda: canon.canon_sols([L(out, "i2_vsearch8_full_r0.log.SOLlines_0.txt")]))))
check("i2_vbrute_nearmiss", "impl2 vbrute8.py on the 30 near-misses: min party weight 4 and A_4 in {6, 12, 18, 30}",
      "REVIEW.md (argument.md 5.14(c)): 18 x A_4 = 30, 6 x 6, 2 x 12, 4 x 18")((
          lambda: ["checked=30 min_party_weight_histogram={4: 30}", "min_weight_and_multiplicity_histogram={(4, 6): 6, (4, 12): 2, (4, 18): 4, (4, 30): 18}"],
          lambda out: guarded([L(out, "i2_vbrute_nearmiss.log")], lambda: canon.canon_vbrute([L(out, "i2_vbrute_nearmiss.log")]))))
check("i2_ctrl_counts_n4_n6", "impl2 vcount.py: exact numbers of labelled stabilizer AME(4,4) and AME(6,4) states in the V-picture",
      "REVIEW.md (argument.md 5.12: 12 of 28 forms, 8640 pairs, N(4) = 4478976000, N(6) = 720^5 * 54400; 5.18: the 11 t values)")((
          lambda: ["n4 admissible_forms=12/28 ordered_pairs=8640 N4=4478976000", "n6 sp_orbits=11",
                   "n6 t_per_orbit_rep=[2272, 288, 368, 288, 288, 64, 112, 368, 112, 368, 112]",
                   "n6 sum_orb_t=54400 ordered_triples=%d N6=%d" % (720 ** 2 * 54400, 720 ** 5 * 54400)],
          lambda out: guarded([L(out, "i2_vcount.log")], lambda: canon.canon_vcount([L(out, "i2_vcount.log")]))))
check("i2_classcount", "impl2 classcount.py: orbit sizes of form tuples under the full GL(4,2) x S_h (explicit congruence, h = 2, 3, 4)",
      "REVIEW.md (argument.md 5.18) / audits.json paper_subsection_classify_review")((
          lambda: ["h=2 tuples=336/784 orbits=1 sizes=[336] total=336",
                   "h=3 tuples=9856/21952 orbits=5 sizes=[5040, 2520, 1260, 1008, 28] total=9856",
                   "h=4 tuples=268800/614656 orbits=10 sizes=[80640, 60480, 40320, 20160, 20160, 20160, 13440, 10080, 2016, 1344] total=268800"],
          lambda out: guarded([L(out, "i2_classcount.log")], lambda: canon.canon_classcount([L(out, "i2_classcount.log")]))))
check("i2_enumtest", "impl2 enumtest.c: the subspace enumerator of vsearch8.c is exact (Gaussian binomials; random candidate sets vs brute force)",
      "REVIEW.md (argument.md 5.12): 651, 11811, 200787, exact agreement with brute force")((
          lambda: ["U=6 subspaces=651", "U=7 subspaces=11811", "U=8 subspaces=200787", "random_trials=3 dfs_equals_bruteforce_over_20160=yes"],
          lambda out: guarded([L(out, "i2_enumtest.log")], lambda: canon.canon_enumtest([L(out, "i2_enumtest.log")]))))


# ---- cross checks (computed from the canonical texts above)
def cross_n4_n6(i1n4, i1n6, i2c):
    def grab(lines, pat):
        for l in lines:
            m = re.search(pat, l)
            if m:
                return m.group(1)
        return "?"
    a4 = grab(i1n4, r"N_labelled_stabilizer_states=(\d+)")
    a6 = grab(i1n6, r"N_labelled_stabilizer_states=(\d+)")
    b4 = grab(i2c, r"N4=(\d+)")
    b6 = grab(i2c, r"N6=(\d+)")
    c1 = sorted(int(x) for x in re.findall(r"nsol=(\d+)", "\n".join(l for l in i1n6 if l.startswith("class="))))
    t = grab(i2c, r"t_per_orbit_rep=(\[[\d, ]+\])")
    c2 = sorted(set(int(x) for x in re.findall(r"\d+", t))) if t != "?" else []
    return ["N4 impl1=%s impl2=%s equal=%s" % (a4, b4, "yes" if a4 == b4 != "?" else "NO"),
            "N6 impl1=%s impl2=%s equal=%s" % (a6, b6, "yes" if a6 == b6 != "?" else "NO"),
            "n6 completions per class impl1 (sorted)=%s, distinct per-Sp-orbit counts impl2=%s equal=%s" % (c1, c2, "yes" if c1 == c2 and c1 else "NO")]


DERIVED = ("cross_n4_n6", "crosswalk", "representatives")
DESC_DERIVED = {
    "cross_n4_n6": ("both implementations give the same N(4), N(6) and the same n = 6 completion counts", "computed from the expected texts i1_ctrl_n4, i1_ctrl_n6, i2_ctrl_counts_n4_n6"),
    "crosswalk": ("class correspondence by a gauge invariant (tools/crosswalk.py) and equality of orbit sizes, #V_1, level-2 counts and 0 solutions per matched pair, also for the 72 Sp-only representatives",
                  "tools/crosswalk.py run on the original representative files and the expected texts"),
    "representatives": ("the ten V_0 representatives of both implementations, explicit, with their correspondence", "tools/crosswalk.py table on the original representative files"),
}


def write(dest, name, lines):
    with open(os.path.join(dest, name + ".txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


def produce(mode, dest, out=None):
    os.makedirs(dest, exist_ok=True)
    texts = {}
    for name, desc, src, fexp, frep in CHECKS:
        try:
            lines = fexp() if mode == "expected" else frep(out)
        except SystemExit as e:  # canon raises on conflicting records
            lines = ["ERROR %s" % e]
        except Exception as e:  # missing or malformed input: record it, the comparison then fails
            lines = ["ERROR %s: %s" % (type(e).__name__, e)]
        texts[name] = lines
        write(dest, name, lines)
    write(dest, "cross_n4_n6", cross_n4_n6(texts["i1_ctrl_n4"], texts["i1_ctrl_n6"], texts["i2_ctrl_counts_n4_n6"]))
    if mode == "expected":
        reps = (os.path.join(O1, "v0_m4.txt"), os.path.join(O2, "v0reps8.txt"), os.path.join(O2, "v0reps8_sp.txt"))
    else:
        reps = (W(out, "v0_m4.txt"), W(out, "v0reps8.txt"), W(out, "v0reps8_sp.txt"))
    cf = [os.path.join(dest, n + ".txt") for n in ("i1_search_t5_frames3", "i1_level2_m124", "i2_search_full", "i2_search_sp72")]
    write(dest, "crosswalk", guarded(list(reps), lambda: tool("crosswalk.py", "check", *(list(reps) + cf))))
    write(dest, "representatives", guarded(list(reps[:2]), lambda: tool("crosswalk.py", "table", *reps[:2])))


def main(a):
    if a[0] == "list":
        for name, desc, src, _, _ in CHECKS:
            print("%-24s %s\n%-24s   expected from: %s" % (name, desc, "", src))
        for name in DERIVED:
            print("%-24s %s\n%-24s   expected from: %s" % (name, DESC_DERIVED[name][0], "", DESC_DERIVED[name][1]))
    elif a[0] == "names":
        print("\n".join([c[0] for c in CHECKS] + list(DERIVED)))
    elif a[0] == "expected":
        produce("expected", a[1])
    elif a[0] == "replay":
        produce("replay", a[2], out=a[1])


if __name__ == "__main__":
    main(sys.argv[1:])
