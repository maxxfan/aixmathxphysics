#!/usr/bin/env python3
"""Exact replay of every window advertised in the paper, from this release alone.

usage:  python3 replay/replay.py [--jobs 8] [--only A19,J97] [--tests] [--controls] [--out replay_report.json]

Needs Python >= 3.9 and a C compiler (`cc`) for the verifier's residual kernel; no third-party packages.
Run it from any directory; it reads and writes nothing outside the release except the --out report.

 1. SHA256SUMS: every file of the release has the listed hash, and no certificate is unlisted.
 2. Pins: the verifier sources and frozen data have the hashes in verifier/PINS.json.
 3. Floating-point environment (replay/fpenv.c, built with the verifier's flags), plus the negative
    control built with -ffp-contract=fast.
 4. For each of the 239 inequalities in manifest.json, every certificate is checked twice:
      a. by the frozen verifier (verifier/verify.py; exact strategy value computed two ways; exact SOS
         identity; rigorous shift tau from the residual kernel verifier/rescheck.c), and
      b. by replay/exactcheck.py, written independently: its own exact strategy evaluation, its own word
         reduction and identity check, and positive semidefiniteness of Y + tau I proven AT THE STATED
         SHIFT in exact integer arithmetic (plus exact Bareiss elimination for N <= 120).
    and compared with the manifest: exact violation / lambda / N / shift, and the printed window
    (printed lower <= exact violation, exact bound <= printed upper).
 5. Aggregates: the counts in the paper (229 two-sided, 223 below 1e-7, 146 below 1e-8), no certified lower
    bound above a certified upper bound, Appendix A (appendix_A.tex, the LaTeX of the paper) and Table I
    (table_I.tsv) regenerated from the manifest, the ladder certificates, and the A19 case study.
 6. --tests runs the verifier's own test suite; --controls runs the negative controls (replay/controls.py),
    each of which must be rejected.
Exit status 0 iff every check passes.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import pathlib
import platform
import subprocess
import sys
import tempfile
import time
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verifier"))
sys.path.insert(0, str(ROOT / "replay"))
sys.path.insert(0, str(ROOT / "generate"))
import exactcheck as X  # noqa: E402
import verify as V      # noqa: E402

UNLISTED_OK = ("__pycache__", ".heavy.lock", "rescheck", ".DS_Store")


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def manifest():
    return json.loads((ROOT / "manifest.json").read_text())


def catalog():
    return json.loads((ROOT / "data" / "catalog.json").read_text())["inequalities"]


# ------------------------------------------------------------------ 1-3

def check_sums():
    listed, bad = {}, []
    for line in (ROOT / "SHA256SUMS").read_text().splitlines():
        h, name = line.split("  ", 1)
        listed[name] = h
        p = ROOT / name
        if not p.is_file():
            bad.append(f"missing {name}")
        elif sha256(p) != h:
            bad.append(f"hash mismatch {name}")
    extra = [p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file()
             and p.relative_to(ROOT).as_posix() not in listed and p.name != "SHA256SUMS"
             and not any(part in UNLISTED_OK for part in p.relative_to(ROOT).parts)]
    bad += [f"unlisted certificate {x}" for x in extra if x.startswith(("certificates/", "ladders/certificates/"))]
    return {"files": len(listed), "problems": bad, "unlisted_other": [x for x in extra if not x.startswith("certificates/")]}


def check_pins():
    want = json.loads((ROOT / "verifier" / "PINS.json").read_text())
    cur = {f: V.sha(V.HERE / f) for f in V.PINNED}
    cur["catalog.json"], cur["records.json"] = V.sha(V.CATALOG), V.sha(V.RECORDS)
    bad = [f for f in cur if want.get(f) != cur[f]]
    return {"pinned": sorted(cur), "problems": [f"modified {f}" for f in bad],
            "note": "the binary verifier/rescheck is built locally from rescheck.c; its hash is machine-specific"}


def check_fpenv():
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        for tag, flags in (("verifier_flags", ["-O2", "-ffp-contract=off", "-fno-fast-math"]),
                           ("negative_control_contract_fast", ["-O2", "-ffp-contract=fast"])):
            exe = os.path.join(tmp, "fpenv_" + tag)
            subprocess.run(["cc", *flags, "-o", exe, str(ROOT / "replay" / "fpenv.c"), "-lm"], check=True)
            p = subprocess.run([exe], capture_output=True, text=True)
            out[tag] = {"flags": " ".join(flags), "output": p.stdout.split("\n")[:-1], "ok": p.returncode == 0}
    cc = subprocess.run(["cc", "--version"], capture_output=True, text=True).stdout.splitlines()
    out["compiler"] = cc[0] if cc else "?"
    out["machine"] = f"{platform.system()} {platform.machine()} python {platform.python_version()}"
    out["problems"] = [] if out["verifier_flags"]["ok"] else ["floating-point environment check failed"]
    return out


# ------------------------------------------------------------------ 4. per inequality

def compare_lower(e, cand, rep, ex):
    """Problems (list of str) comparing a strategy's two exact evaluations with the manifest and the window."""
    m, bad = e["lower"], []
    if rep.get("verdict") in ("INVALID", "VERIFIER_ALARM"):
        bad.append(f"verifier: {rep.get('verdict')} {rep.get('reason', '')}")
        return bad
    if rep.get("violation") != m["violation_down10"]:
        bad.append(f"verifier violation {rep.get('violation')} != manifest {m['violation_down10']}")
    if rep.get("violation_exact") and Q(rep["violation_exact"]) != Q(m["violation_exact"]):
        bad.append("verifier exact violation differs from manifest")
    if Q(ex["violation_exact"]) != Q(m["violation_exact"]):
        bad.append("independent exact violation differs from manifest")
    if not Q(e["window"]["lower"]) <= Q(ex["violation_exact"]):
        bad.append(f"printed lower {e['window']['lower']} exceeds the exact violation")
    if [cand["dA"], cand["dB"]] != m["dims"] or cand.get("field", "real") != m["field"]:
        bad.append("dimension or field differs from manifest")
    return bad


def compare_upper(e, cand, rep, ex, b0):
    m, bad = e["upper"], []
    if rep.get("verdict") in ("INVALID", "VERIFIER_ALARM"):
        bad.append(f"verifier: {rep.get('verdict')} {rep.get('reason', '')}")
        return bad, {}
    lam, n, tau = Q(m["lambda_exact"]), m["n_words"], Q(m["tau_stated_exact"])
    if tau != Q(float.fromhex(m["sigma_hex"])) + Q(float.fromhex(m["residual_hex"])):
        bad.append("manifest tau_stated is not sigma + residual")
    if Q(m["bound_exact"]) != lam + n * tau - b0:
        bad.append("manifest bound_exact is not lambda + N tau - b0")
    if rep.get("n_words") != n:
        bad.append("verifier N differs from manifest")
    if not (lam <= Q(rep["lambda"]) < lam + Q(1, 10 ** 10)):
        bad.append("verifier lambda differs from manifest")
    if Q(ex["lambda_exact"]) != lam:
        bad.append("independent lambda differs from manifest")
    if not Q(rep["bound"]) <= Q(e["window"]["upper"]):
        bad.append(f"verifier bound {rep['bound']} exceeds printed upper {e['window']['upper']}")
    if not ex.get("tau_ok"):
        bad.append(f"exact shift {ex.get('tau_exact_float')} exceeds stated shift {float(tau)}")
    if ex.get("bareiss_at_stated_tau") is False:
        bad.append("Bareiss: Y + tau_stated I is not positive definite")
    if not Q(m["bound_exact"]) <= Q(e["window"]["upper"]):
        bad.append("manifest bound exceeds printed upper")
    if not Q(ex["bound_exact_replay"]) <= Q(e["window"]["upper"]):
        bad.append("independent bound lambda + N tau_exact - b0 exceeds printed upper")
    same = rep.get("sigma") == float.fromhex(m["sigma_hex"]) and rep.get("residual_bound") == float.fromhex(m["residual_hex"])
    sigma = Q(float.fromhex(m["sigma_hex"]))
    f_stated, f_exact = tau - sigma, Q(ex["tau_exact"]) - sigma if "tau_exact" in ex else None
    return bad, {"verifier_shift_bit_identical": same, "tau_stated": float(tau), "tau_exact": ex.get("tau_exact_float"),
                 "f_exact_over_f_stated": float(f_exact / f_stated) if f_exact is not None and f_stated else None,
                 "bareiss_at_stated_tau": ex.get("bareiss_at_stated_tau")}


def replay_one(e):
    """Replay all certificates of one manifest entry."""
    cat = catalog()
    b0 = Q(cat[e["name"]]["b0"])
    out = {"name": e["name"], "problems": [], "seconds": {}}
    for side in ("lower", "upper"):
        m = e[side]
        if m is None:
            continue
        t0 = time.time()
        path = ROOT / m["file"]
        if sha256(path) != m["sha256"]:
            out["problems"].append(f"{side}: sha256 differs from manifest")
            continue
        cand = json.loads(path.read_text())
        if cand.get("inequality") != e["name"]:
            out["problems"].append(f"{side}: certificate is for {cand.get('inequality')}")
            continue
        try:
            rep = V.check_candidate(cand, lock=False)
        except Exception as err:  # noqa: BLE001
            rep = {"verdict": "INVALID", "reason": f"{type(err).__name__}: {err}"}
        try:
            if side == "lower":
                ex = X.check(cand, cat)
                bad = compare_lower(e, cand, rep, ex)
            else:
                ex = X.check(cand, cat, float.fromhex(m["sigma_hex"]), Q(m["tau_stated_exact"]))
                bad, info = compare_upper(e, cand, rep, ex, b0)
                out["upper_info"] = info
        except X.Bad as err:
            bad = [f"independent checker rejects: {err}"]
        out["problems"] += [f"{side}: {b}" for b in bad]
        out["seconds"][side] = round(time.time() - t0, 2)
    for side in ("lower", "upper"):
        if e["window"][side] and e[side] is None:
            out["problems"].append(f"printed {side} bound without a certificate")
    return out


# ------------------------------------------------------------------ 5. aggregates

def check_aggregates(man):
    bad = []
    es = man["inequalities"]
    both = [e for e in es if e["window"]["lower"] and e["window"]["upper"]]
    widths = [Q(e["window"]["upper"]) - Q(e["window"]["lower"]) for e in both]
    counts = {"inequalities": len(es), "both": len(both), "width_below_1e-7": sum(w < Q(1, 10**7) for w in widths),
              "width_below_1e-8": sum(w < Q(1, 10**8) for w in widths), "negative_widths": sum(w < 0 for w in widths)}
    for k in ("inequalities", "both", "width_below_1e-7", "width_below_1e-8"):
        if counts[k] != man["counts"][k]:
            bad.append(f"count {k}: {counts[k]} != manifest {man['counts'][k]}")
    if (counts["inequalities"], counts["both"], counts["width_below_1e-7"], counts["width_below_1e-8"]) != (239, 229, 223, 146):
        bad.append("counts differ from the paper (239, 229, 223, 146)")
    if counts["negative_widths"]:
        bad.append("a certified lower bound exceeds a certified upper bound")
    import make_appendix
    tex, tsv = make_appendix.render(man)
    if (ROOT / "appendix_A.tex").read_text() != tex or (ROOT / "appendix_A.tsv").read_text() != tsv:
        bad.append("appendix_A.tex/.tsv do not match the manifest")
    byname = {e["name"]: e for e in es}
    rows = (ROOT / "table_I.tsv").read_text().splitlines()[1:]
    for r in rows:
        n, lo0, up0, lo, up, cl, cu = r.split("\t")
        if cl == "yes" and lo != byname[n]["window"]["lower"]:
            bad.append(f"Table I {n}: lower {lo} is not the certified {byname[n]['window']['lower']}")
        if cu == "yes" and up != byname[n]["window"]["upper"]:
            bad.append(f"Table I {n}: upper {up} is not the certified {byname[n]['window']['upper']}")
        if cl == cu == "yes" and Q(lo) > Q(up):  # mixed rows pair a certified value with a rounded literature one
            bad.append(f"Table I {n}: certified lower exceeds certified upper")
    return {"counts": counts, "table_I_rows": len(rows), "problems": bad}


def check_ladders():
    cat, bad, n = catalog(), [], 0
    lines = (ROOT / "ladders" / "ladders.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    for line in lines[1:]:
        r = dict(zip(head, line.split("\t")))
        if r["status"] != "certified":
            continue
        n += 1
        p = ROOT / r["certificate"]
        if sha256(p) != r["sha256"]:
            bad.append(f"{r['inequality']} d={r['local_dim']}: sha256 mismatch")
            continue
        cand = json.loads(p.read_text())
        v = Q(X.check(cand, cat)["violation_exact"])
        if not (Q(r["value"]) <= v < Q(r["value"]) + Q(1, 10 ** 10)) or cand["dA"] != int(r["local_dim"]):
            bad.append(f"{r['inequality']} d={r['local_dim']}: value {float(v)!r} or dimension differs")
    return {"certified_rungs": n, "problems": bad}


def check_A19():
    p = subprocess.run([sys.executable, str(ROOT / "case_studies" / "A19" / "check_A19.py")], capture_output=True, text=True)
    return {"output_tail": p.stdout.splitlines()[-3:], "problems": [] if p.returncode == 0 else ["A19 case study failed"]}


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--only", help="comma-separated inequality names")
    ap.add_argument("--tests", action="store_true", help="also run verifier/tests/test_verifier.py")
    ap.add_argument("--controls", action="store_true", help="also run the negative controls")
    ap.add_argument("--out", default="replay_report.json")
    a = ap.parse_args()
    t0 = time.time()
    man = manifest()
    report = {"release": man["release"], "root": str(ROOT)}
    for key, fn in (("sha256sums", check_sums), ("pins", check_pins), ("fpenv", check_fpenv)):
        report[key] = fn()
        print(f"[{key}] {'ok' if not report[key]['problems'] else report[key]['problems']}", flush=True)
    V.ensure_rescheck()  # build the kernel once, before the workers start
    es = man["inequalities"]
    if a.only:
        keep = set(a.only.split(","))
        es = [e for e in es if e["name"] in keep]
    es = sorted(es, key=lambda e: -((e["upper"] or {}).get("n_words", 0) + 10 * max((e["lower"] or {}).get("dims", [0]))))
    with mp.get_context("spawn").Pool(a.jobs) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(replay_one, es), 1):
            res.append(r)
            if r["problems"]:
                print(f"  FAIL {r['name']}: {r['problems']}", flush=True)
            if i % 40 == 0 or i == len(es):
                print(f"  {i}/{len(es)} inequalities replayed ({time.time() - t0:.0f} s)", flush=True)
    res.sort(key=lambda r: [e["name"] for e in man["inequalities"]].index(r["name"]))
    n_cert = sum(len(r["seconds"]) for r in res)
    fails = [r for r in res if r["problems"]]
    ups = [r["upper_info"] for r in res if "upper_info" in r]
    report["certificates"] = {
        "inequalities": len(res), "certificates_replayed": n_cert, "failed": [r["name"] for r in fails],
        "uppers_with_bareiss_at_stated_tau": sum(1 for u in ups if u["bareiss_at_stated_tau"]),
        "uppers_verifier_shift_bit_identical": sum(1 for u in ups if u["verifier_shift_bit_identical"]),
        "max_tau_exact_over_tau_stated": max((u["tau_exact"] / u["tau_stated"] for u in ups if u["tau_stated"]), default=None),
        "max_f_exact_over_f_stated": max((u["f_exact_over_f_stated"] for u in ups if u["f_exact_over_f_stated"] is not None), default=None),
        "problems": [f"{r['name']}: {p}" for r in fails for p in r["problems"]], "details": res}
    print(f"[certificates] {n_cert} certificates of {len(res)} inequalities: "
          f"{'all pass' if not fails else str(len(fails)) + ' inequalities FAIL'}; "
          f"Bareiss at stated tau {report['certificates']['uppers_with_bareiss_at_stated_tau']}/{len(ups)}; "
          f"verifier shift bit-identical {report['certificates']['uppers_verifier_shift_bit_identical']}/{len(ups)}; "
          f"max tau_exact/tau_stated {report['certificates']['max_tau_exact_over_tau_stated']:.6f}, "
          f"max f_exact/f_stated {report['certificates']['max_f_exact_over_f_stated']:.3f}", flush=True)
    if not a.only:
        for key, fn in (("aggregates", lambda: check_aggregates(man)), ("ladders", check_ladders), ("A19", check_A19)):
            report[key] = fn()
            print(f"[{key}] {'ok' if not report[key]['problems'] else report[key]['problems']}"
                  + (f" {report[key]['counts']}" if key == "aggregates" else ""), flush=True)
    if a.tests:
        p = subprocess.run([sys.executable, str(ROOT / "verifier" / "tests" / "test_verifier.py"), "--trials", "10"],
                           capture_output=True, text=True)
        report["verifier_tests"] = {"returncode": p.returncode, "tail": p.stdout.splitlines()[-5:],
                                    "problems": [] if p.returncode == 0 else ["verifier tests failed"]}
        print(f"[verifier_tests] {'ok' if p.returncode == 0 else 'FAIL'}: {report['verifier_tests']['tail'][-1:]}", flush=True)
    if a.controls:
        import controls
        report["controls"] = controls.run_all(man)
        print(f"[controls] {report['controls']['rejected']}/{report['controls']['total']} rejected as required"
              + (f"; NOT rejected: {report['controls']['problems']}" if report["controls"]["problems"] else ""), flush=True)
    problems = [p for k, v in report.items() if isinstance(v, dict) for p in v.get("problems", [])]
    report["seconds"] = round(time.time() - t0, 1)
    report["verdict"] = "PASS" if not problems else "FAIL"
    pathlib.Path(a.out).write_text(json.dumps(report, indent=1) + "\n")
    print(f"REPLAY {report['verdict']} in {report['seconds']} s (report: {a.out})")
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
