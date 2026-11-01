#!/usr/bin/env python3
"""Standalone exact checker for the p10 AME release (Python 3.9+ stdlib, plus a C compiler for cutrank.c).

It shares no code with the project verifier. It checks stabilizer states given in the canonical F_p form:

  {"type": "stabilizer", "n": n, "d": d, "p": p, "generators": [...]}
  n parties of local dimension d = p^m; N = m*n qudits of prime dimension p, party j owns qudits j*m .. j*m+m-1;
  N rows (a | b) in F_p^(2N) for the Pauli operators X(a)Z(b); rows are digit strings (base p, p <= 36) or int lists.

and coupling certificates

  {"type": "coupling", "k": k, "a": <stabilizer, 2k+2 parties, dim a>, "b": <stabilizer, 2k+2 parties, dim b>}.

Mathematics (paper, Lemmas 1 and 3):
  * The rows define a pure stabilizer state iff they are isotropic for a.b' - b.a' (mod p) and have rank N.
  * Every set of k parties is maximally mixed (k-uniformity) iff no nonzero stabilizer element is supported on k
    parties iff, for every set B of n-k parties, the N x 2m(n-k) block of the X and Z columns of B has rank N.
    AME(n,d) is the case k = floor(n/2).
  * Coupling (Lemma 3): if a is k-uniform on 2k+2 parties of dimension a, b is AME(2k+2, b) and a <= b, the
    contracted state is AME(2k+1, ab).

usage:
  ame_check.py state FILE [FILE...]          AME check of stabilizer files
  ame_check.py uniform K FILE               K-uniformity check
  ame_check.py coupling FILE [FILE...]       coupling certificates (checks both ingredients)
Options: --python-only (no C; exhaustive in Python), --sample S (Python re-check of S seeded random sets, default 200).
Exit status 0 iff every check passes.
"""
import itertools
import json
import os
import pathlib
import random
import subprocess
import sys
import tempfile
from math import comb

HERE = pathlib.Path(__file__).resolve().parent
DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz"
_EXE = None


def is_prime(x):
    return x >= 2 and all(x % k for k in range(2, int(x ** 0.5) + 1))


def load_rows(obj):
    p, n, d = int(obj["p"]), int(obj["n"]), int(obj["d"])
    if not is_prime(p):
        raise ValueError(f"p = {p} is not prime")
    m, x = 0, d
    while x % p == 0 and x > 1:
        x //= p
        m += 1
    if x != 1 or m == 0:
        raise ValueError(f"d = {d} is not a power of p = {p}")
    N = m * n
    rows = []
    for r in obj["generators"]:
        if isinstance(r, str):
            row = [DIGITS.index(ch) for ch in r.strip().lower()]
        else:
            row = [int(v) for v in r]
        if len(row) != 2 * N or any(not (0 <= v < p) for v in row):
            raise ValueError("bad generator row")
        rows.append(row)
    if len(rows) != N:
        raise ValueError(f"expected {N} rows, got {len(rows)}")
    return p, n, d, m, N, rows


def rank_mod_p(mat, p):
    M = [list(r) for r in mat]
    rk, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        pv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if pv is None:
            continue
        M[rk], M[pv] = M[pv], M[rk]
        iv = pow(M[rk][c], p - 2, p)
        M[rk] = [(v * iv) % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(u - f * v) % p for u, v in zip(M[i], M[rk])]
        rk += 1
    return rk


def block_cols(B, m, N):
    cols = []
    for j in B:
        cols += [j * m + t for t in range(m)]
        cols += [N + j * m + t for t in range(m)]
    return cols


def isotropic(rows, p, N):
    for i in range(N):
        for j in range(i + 1, N):
            s = sum(rows[i][t] * rows[j][N + t] - rows[i][N + t] * rows[j][t] for t in range(N)) % p
            if s:
                return False
    return True


def c_count(rows, p, n, m, t):
    global _EXE
    if _EXE is None:
        td = tempfile.mkdtemp(prefix="amecheck_")
        _EXE = os.path.join(td, "cutrank")
        subprocess.run(["cc", "-O2", "-o", _EXE, str(HERE / "cutrank.c")], check=True)
    N = len(rows)
    inp = f"{p} {N} {n} {m} {t}\n" + "\n".join(" ".join(map(str, r)) for r in rows) + "\n"
    out = subprocess.run([_EXE], input=inp, capture_output=True, text=True, check=True).stdout.split()
    i = out.index("SETS")
    return int(out[i + 1]), int(out[i + 3])


def uniformity(obj, k, python_only=False, sample=200, label=""):
    p, n, d, m, N, rows = load_rows(obj)
    res = {"label": label, "n": n, "d": d, "p": p, "m": m, "N": N, "k": k}
    res["isotropic"] = isotropic(rows, p, N)
    res["rank"] = rank_mod_p(rows, p)
    t = n - k
    total = comb(n, t)
    sets = list(itertools.combinations(range(n), t))
    if python_only:
        deficient = sum(1 for B in sets if rank_mod_p([[r[c] for c in block_cols(B, m, N)] for r in rows], p) < N)
        res["sets"], res["deficient"], res["method"] = total, deficient, "python exhaustive"
    else:
        s, deficient = c_count(rows, p, n, m, t)
        assert s == total
        rnd = random.Random(1000003 * n + d)
        pick = rnd.sample(sets, min(sample, total))
        py_def = sum(1 for B in pick if rank_mod_p([[r[c] for c in block_cols(B, m, N)] for r in rows], p) < N)
        res["sets"], res["deficient"] = total, deficient
        res["method"] = f"C exhaustive + python sample {len(pick)} (deficient in sample: {py_def})"
        if deficient == 0 and py_def:
            res["disagreement"] = True
    res["ok"] = res["isotropic"] and res["rank"] == N and res["deficient"] == 0 and not res.get("disagreement")
    return res


def show(r):
    print(f"  {r['label']:44s} n={r['n']:2d} d={r['d']:2d} (F_{r['p']}, m={r['m']}, N={r['N']}) k={r['k']}: "
          f"isotropic={'yes' if r['isotropic'] else 'NO'} rank={r['rank']}/{r['N']} "
          f"sets of {r['n'] - r['k']} parties: {r['sets']}, full rank {r['sets'] - r['deficient']}, deficient {r['deficient']} "
          f"[{r['method']}] -> {'PASS' if r['ok'] else 'FAIL'}")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    py = "--python-only" in sys.argv
    sample = 200
    if "--sample" in sys.argv:
        sample = int(sys.argv[sys.argv.index("--sample") + 1])
        args = [a for a in args if a != str(sample)]
    mode, files = args[0], args[1:]
    ok = True
    if mode == "uniform":
        k, files = int(files[0]), files[1:]
    for f in files:
        obj = json.loads(pathlib.Path(f).read_text())
        name = pathlib.Path(f).name
        if mode == "state":
            r = uniformity(obj, int(obj["n"]) // 2, py, sample, name)
            show(r)
            ok &= r["ok"]
        elif mode == "uniform":
            r = uniformity(obj, k, py, sample, name)
            show(r)
            ok &= r["ok"]
        elif mode == "coupling":
            kk = int(obj["k"])
            a, b = obj["a"], obj["b"]
            na, nb, da, db = int(a["n"]), int(b["n"]), int(a["d"]), int(b["d"])
            good = na == nb == 2 * kk + 2 and da <= db
            ra = uniformity(a, kk, py, sample, f"{name} : a ({kk}-uniform)")
            rb = uniformity(b, kk + 1, py, sample, f"{name} : b (AME({nb},{db}))")
            show(ra)
            show(rb)
            if not good:
                print(f"  {name}: shape check FAIL (parties of a and b: {na}, {nb}, need {2 * kk + 2}; "
                      f"a = {da} <= b = {db}: {'yes' if da <= db else 'no'})")
            cell_ok = good and ra["ok"] and rb["ok"]
            print(f"  => {name}: AME({2 * kk + 1},{da * db}) by coupling (k={kk}, a={da} <= b={db}): {'PASS' if cell_ok else 'FAIL'}")
            ok &= cell_ok
        else:
            sys.exit("mode must be state, uniform or coupling")
    print("OVERALL", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
