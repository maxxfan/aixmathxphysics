#!/usr/bin/env python3
"""A local-Clifford invariant of the circulant AME(15,8) graph codes (stdlib only, own GF(8) arithmetic).

For a stabilizer group S (an F_2-space) and a set B of parties let S_B = {s in S : supp(s) inside B}.  For an AME(15,8)
stabilizer state and |B| = 9, dim_F2 S_B = 9, and for each i in B the subspace S_{B - i} has dimension 3.  For i, j, l
in B the number d = dim_F2 (S_{B-i} + S_{B-j} + S_{B-l}) is unchanged by local symplectic maps (they are F_2-linear and
preserve supports), and the whole collection is permuted by party permutations.  Hence

    profile(S) = multiset over the 5,005 sets B of  #{triples {i,j,l} in B : d < 9}

is the same for two states that are equal up to local Clifford operations and a permutation of the parties.  For an
F_8-linear code S_B is an F_8-space of dimension 3, the S_{B-i} are nine points of the projective plane PG(2,8), and
d < 9 means that three of them are collinear.

The script computes the profile for every orbit representative of the exhaustive scan (28 for n = 15), for the
vector printed in Theorem 1 and for the graph form of the code shown on codetables.de, and checks that it is constant
on the classes found by compare_codetables.py (a consistency test of the invariant and of the classes) and whether it
separates the classes.  Different profiles prove that two states are NOT equal up to local Clifford operations and
party permutations.  Equal profiles would prove nothing.

usage: python3 paper/lc_invariant.py          (reads results/codetables_comparison_2026-11-10.json, paper/circ_hits.json;
                                               writes results/lc_invariant_2026-11-10.json; about 40 s)
       python3 lc_invariant.py COMPARISON_JSON HITS_JSON [OUT_JSON]      (as used in the release, ancillary/codetables/)
"""
import itertools
import json
import pathlib
import sys
from collections import Counter

P = pathlib.Path(__file__).resolve().parents[1]
EXP = [1]
for _ in range(6):
    x = EXP[-1] << 1
    if x & 8:
        x ^= 0b1011
    EXP.append(x)
LOG = {v: i for i, v in enumerate(EXP)}
MUL = [[0 if a == 0 or b == 0 else EXP[(LOG[a] + LOG[b]) % 7] for b in range(8)] for a in range(8)]
INV = [0] + [EXP[(-LOG[a]) % 7] for a in range(1, 8)]


def kernel(M, ncols):
    """basis of {y : y M = 0} for the rows-by-ncols matrix M over GF(8) (left kernel), by elimination on [M | I]"""
    r = len(M)
    A = [list(M[i]) + [1 if j == i else 0 for j in range(r)] for i in range(r)]
    rk = 0
    for c in range(ncols):
        p = next((i for i in range(rk, r) if A[i][c]), None)
        if p is None:
            continue
        A[rk], A[p] = A[p], A[rk]
        iv = INV[A[rk][c]]
        A[rk] = [MUL[iv][x] for x in A[rk]]
        for i in range(r):
            if i != rk and A[i][c]:
                f = A[i][c]
                A[i] = [x ^ MUL[f][y] for x, y in zip(A[i], A[rk])]
        rk += 1
    return [row[ncols:] for row in A[rk:]]


def det3(a, b, c):
    m = MUL
    return (m[a[0]][m[b[1]][c[2]] ^ m[b[2]][c[1]]] ^ m[a[1]][m[b[0]][c[2]] ^ m[b[2]][c[0]]]
            ^ m[a[2]][m[b[0]][c[1]] ^ m[b[1]][c[0]]])


def profile(c, n=15):
    h = n // 2
    G = [[0 if i == j else c[min((j - i) % n, (i - j) % n) - 1] for j in range(n)] for i in range(n)]
    prof = Counter()
    for B in itertools.combinations(range(n), h + 2):
        Bc = [j for j in range(n) if j not in B]
        K = kernel([[G[i][j] for j in Bc] for i in B], len(Bc))   # y in F_8^B with y G[B, B^c] = 0
        assert len(K) == 3, "dim S_B"
        pts = []
        for idx, i in enumerate(B):
            # evaluation of the basis words at party i: (x_i, (x G)_i)
            ev = [[y[idx], 0] for y in K]
            for t, y in enumerate(K):
                acc = 0
                for jdx, j in enumerate(B):
                    if y[jdx]:
                        acc ^= MUL[y[jdx]][G[j][i]]
                ev[t][1] = acc
            ker = kernel(ev, 2)          # coefficient vectors (in the basis K) of the words vanishing at party i
            assert len(ker) == 1, "dim S_{B-i}"
            pts.append(ker[0])
        prof[sum(1 for a, b, d in itertools.combinations(pts, 3) if det3(a, b, d) == 0)] += 1
    return dict(sorted(prof.items()))


def canon(c, n=15):
    """lexicographically least normalized member of the orbit under multipliers, Frobenius and scaling"""
    import math
    h = n // 2
    best = None
    for u in range(1, n):
        if math.gcd(u, n) != 1:
            continue
        v = [0] * h
        for i in range(1, h + 1):
            j = (u * i) % n
            v[min(j, n - j) - 1] = c[i - 1]
        for _ in range(3):
            v = [MUL[x][x] for x in v]
            iv = INV[next(x for x in v if x)]
            w = tuple(MUL[iv][x] for x in v)
            if best is None or w < best:
                best = w
    return best


def main(argv):
    cmp_path = P / "results" / "codetables_comparison_2026-11-10.json"
    hits_path = P / "paper" / "circ_hits.json"
    out_path = P / "results" / "lc_invariant_2026-11-10.json"
    if len(argv) >= 3:
        cmp_path, hits_path = pathlib.Path(argv[1]), pathlib.Path(argv[2])
        out_path = pathlib.Path(argv[3]) if len(argv) > 3 else None
    cmp_ = json.loads(cmp_path.read_text())["cells"][0]
    assert cmp_["cell"] == "[[15,0]]_8"
    classes = cmp_["scan_orbit_classes_under_uniform_local_cliffords"]
    hits = json.loads(hits_path.read_text())["15"]
    reps = sorted(canon(tuple(r)) for r in hits.get("orbit_reps", hits.get("orbit_representatives")))
    assert len(reps) == 28 and len(set(reps)) == 28
    thm = (1, 0, 4, 6, 3, 0, 0)
    table_forms = [tuple(e["c"]) for e in cmp_["graph_forms"]]
    table_orbits = sorted(reps.index(canon(c)) + 1 for c in table_forms)
    thm_orbit = reps.index(canon(thm)) + 1
    assert thm_orbit == cmp_["theorem1_scan_orbit_number"]
    prof = {o: profile(reps[o - 1]) for o in range(1, 29)}
    p_thm = profile(thm)
    p_tab = [profile(c) for c in table_forms]
    cls_prof = []
    ok_const = True
    for cl in classes:
        ps = [prof[o] for o in cl]
        ok_const &= all(x == ps[0] for x in ps)
        cls_prof.append(ps[0])
    ok_distinct = all(cls_prof[i] != cls_prof[j] for i in range(len(classes)) for j in range(i))
    cl_thm = next(i for i, cl in enumerate(classes) if thm_orbit in cl)
    cl_tab = next(i for i, cl in enumerate(classes) if table_orbits[0] in cl)
    ok_members = p_thm == cls_prof[cl_thm] and all(x == cls_prof[cl_tab] for x in p_tab) and \
        all(o in classes[cl_tab] for o in table_orbits)
    lines = ["invariant: for each of the 5,005 nine-sets B, the number of collinear triples among the nine points S_{B-i} of PG(2,8)",
             "profile = {number of collinear triples: number of sets B}"]
    for i, cl in enumerate(classes):
        tags = (" <- Theorem 1 (scan orbit %d)" % thm_orbit if i == cl_thm else "") + \
               (" <- codetables.de [[15,0,8]]_8 (scan orbits %s)" % table_orbits if i == cl_tab else "")
        lines.append(f"class {i + 1}, scan orbits {cl}: {cls_prof[i]}{tags}")
    lines.append(f"profile constant on each class (7 orbits each): {ok_const}")
    lines.append(f"Theorem 1 and all graph forms of the table code have the profile of their class: {ok_members}")
    lines.append(f"the {len(classes)} class profiles are pairwise different: {ok_distinct}")
    sep = cls_prof[cl_thm] != cls_prof[cl_tab]
    lines.append("the codetables.de [[15,0,8]]_8 code and the state of Theorem 1 are "
                 + ("NOT equal up to local Clifford operations and party permutations (profiles differ)" if sep
                    else "not separated by this invariant"))
    for line in lines:
        print(line)
    if out_path is not None:
        out_path.write_text(json.dumps({
            "invariant": "multiset over the 9-subsets B of the number of triples i,j,l in B with dim_F2(S_{B-i}+S_{B-j}+S_{B-l}) < 9",
            "classes": classes, "class_profiles": cls_prof, "theorem1_class": cl_thm + 1, "table_class": cl_tab + 1,
            "theorem1_scan_orbit": thm_orbit, "table_scan_orbits": table_orbits,
            "constant_on_classes": ok_const, "members_match": ok_members, "classes_pairwise_distinct": ok_distinct,
            "table_code_inequivalent_to_theorem1_under_local_cliffords_and_permutations": sep}, indent=1) + "\n")
    good = ok_const and ok_members
    print("INVARIANT DONE" if good else "INVARIANT INCONSISTENT")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
