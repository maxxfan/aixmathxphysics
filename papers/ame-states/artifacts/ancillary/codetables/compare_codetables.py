#!/usr/bin/env python3
"""Compare the [[15,0,8]]_8 and [[17,0,9]]_8 codes shown on codetables.de (pages fetched 2026-11-10, both marked
"last modified: 2026-10-21") with the symmetric circulant graph states of Theorem 1.  Stdlib only; shares no code with
the verifier, qlib or check_constructions.py.  The cached pages are read as data only.

For each cached page data/sources/codetables_2026-11-10/qecc_q64_n{15,17}_k0.html the script
  1. parses the printed stabilizer matrix [X|Z] (entries 0, 1, w^k; one row per F_2-generator);
  2. decides the minimal polynomial of w from the data: for each of the two primitive cubics it takes the F_2-span of
     the rows and tests dimension 3n and isotropy for the trace-symplectic form Tr(a.b' - a'.b);
  3. tests whether the span is closed under multiplication by w (F_8-linear) and computes a reduced F_8-basis;
  4. applies each of the 504 elements T of SL(2,8) as the same local Clifford (a|b) -> (p a + q b | r a + s b) on
     every party (this keeps the code cyclic); whenever the X part becomes an information set the code is
     {(x | x G)} with G symmetric and circulant (asserted), the constant diagonal of G is removed by the local
     Clifford (a|b) -> (a|b + g a), and c = (G_{0,1}, ..., G_{0,floor(n/2)}) is read off;
  5. re-checks, with its own GF(8) elimination, that G[A, A^c] has full row rank for every floor(n/2)-set A
     (all 6,435 resp. 24,310 sets), i.e. that the code is AME(n,8);
  6. computes the orbit of each such c under multipliers i -> u i, the Frobenius map and scaling (canonical form:
     first nonzero coefficient 1, lexicographically least over the orbit) and compares it with the orbit of the
     coefficient vector printed in Theorem 1 and with the orbit representatives of the exhaustive scan
     (paper/circ_hits.json).

What the result does and does not say.  If one of the graph forms of the table's code lies in the orbit of Theorem 1's
vector, the two states are equal up to a party permutation and local Clifford unitaries (a rigorous positive).  If
none does, nothing follows from this script: the states may still be equivalent under non-uniform local unitaries, and
the party-weight enumerator cannot separate them, because it is the same for every AME(n,8) stabilizer state.  The
companion script lc_invariant.py settles that case for n = 15 with a finer invariant (the four classes found here are
pairwise inequivalent under local Clifford operations and party permutations).

usage: python3 paper/compare_codetables.py        (writes results/codetables_comparison_2026-11-10.json; about 40 s)
       python3 compare_codetables.py PAGES_DIR HITS_JSON [OUT_JSON]     (as used in the release, ancillary/codetables/)
"""
import html
import itertools
import json
import math
import pathlib
import re
import sys

P = pathlib.Path(__file__).resolve().parents[1]
SRC = P / "data" / "sources" / "codetables_2026-11-10"
HITS = P / "paper" / "circ_hits.json"   # the release passes ancillary/circulants/expected/hits.json
THEOREM1 = {15: (1, 0, 4, 6, 3, 0, 0), 17: (1, 2, 2, 1, 2, 5, 5, 2)}  # integers: bit i = coefficient of alpha^i


class GF8:
    def __init__(self, modulus):
        self.exp = [1]
        for _ in range(6):
            x = self.exp[-1] << 1
            if x & 8:
                x ^= modulus
            self.exp.append(x)
        assert sorted(self.exp) == list(range(1, 8))
        self.log = {v: i for i, v in enumerate(self.exp)}

    def mul(self, a, b):
        return 0 if a == 0 or b == 0 else self.exp[(self.log[a] + self.log[b]) % 7]

    def inv(self, a):
        return self.exp[(-self.log[a]) % 7]

    def tr(self, a):
        t = a ^ self.mul(a, a) ^ self.mul(self.mul(a, a), self.mul(a, a))
        assert t in (0, 1)
        return t


def parse_page(path):
    text = path.read_text(errors="replace")
    head = re.search(r"Bounds on \[\[(\d+),0\]\]<SUB>(\d+)</SUB>.*?lower bound:</TD><TD>(\d+).*?upper bound:</TD><TD>(\d+)",
                     text, re.S)
    n, q, lo, hi = (int(x) for x in head.groups())
    ctype = re.search(r"Construction type: (.*?)<P>", text).group(1)
    modified = re.search(r"last modified: (\S+)", text).group(1)
    poly = html.unescape(re.search(r"generating polynomial (.*)", text).group(1)).strip()
    rows = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith("[") and "|" in line and line.endswith("]"):
            xs, zs = line[1:-1].split("|")
            rows.append((xs.split(), zs.split()))
    assert all(len(x) == n and len(z) == n for x, z in rows), "row length"
    return {"n": n, "q": q, "lower": lo, "upper": hi, "type": ctype, "modified": modified, "poly": poly, "rows": rows}


def entry(tok, F):
    if tok == "0":
        return 0
    if tok == "1":
        return 1
    m = re.fullmatch(r"w(?:\^(\d+))?", tok)
    assert m, tok
    return F.exp[int(m.group(1) or 1) % 7]


def f2_span_basis(vecs):
    """vecs: integers (bit vectors). Returns a reduced basis as dict pivot -> vector."""
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                break
    return piv


def in_span(piv, v):
    while v:
        h = v.bit_length() - 1
        if h not in piv:
            return False
        v ^= piv[h]
    return True


def pack(vec):
    """tuple of GF(8) elements -> integer (3 bits per element)"""
    r = 0
    for i, a in enumerate(vec):
        r |= a << (3 * i)
    return r


def rank_gf8(M, F):
    M = [list(r) for r in M]
    rk = 0
    ncol = len(M[0])
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        iv = F.inv(M[rk][c])
        M[rk] = [F.mul(iv, x) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [x ^ F.mul(f, y) for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def rref_gf8(M, F):
    M = [list(r) for r in M]
    rk = 0
    pivots = []
    for c in range(len(M[0])):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        iv = F.inv(M[rk][c])
        M[rk] = [F.mul(iv, x) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [x ^ F.mul(f, y) for x, y in zip(M[i], M[rk])]
        pivots.append(c)
        rk += 1
    return M[:rk], pivots


def canon_orbit(c, n, F):
    """orbit of c under multipliers, Frobenius and scaling, as normalized tuples; returns (least, sorted orbit)"""
    h = n // 2

    def norm(v):
        lead = next(x for x in v if x)
        iv = F.inv(lead)
        return tuple(F.mul(iv, x) for x in v)
    orbit = set()
    for u in range(1, n):
        if math.gcd(u, n) != 1:
            continue
        v = [0] * h
        for i in range(1, h + 1):
            j = (u * i) % n
            j = min(j, n - j)
            v[j - 1] = c[i - 1]
        for _ in range(3):
            v = [F.mul(x, x) for x in v]
            orbit.add(norm(v))
    return min(orbit), sorted(orbit)


def graph_forms(R, n, F):
    """R: F_8-basis (n rows, 2n columns) of a cyclic self-dual code.  Apply the same T in SL(2,8) on every party,
    (a|b) -> (p a + q b | r a + s b); when the X part becomes an information set, read off the symmetric circulant G.
    Returns ({least orbit representative: info}, number of T, number of T giving a graph form)."""
    h = n // 2
    reached = {}
    n_T = n_graph = 0
    for p_, q_, r_, s_ in itertools.product(range(8), repeat=4):
        if F.mul(p_, s_) ^ F.mul(q_, r_) != 1:
            continue
        n_T += 1
        M = [[F.mul(p_, row[i]) ^ F.mul(q_, row[n + i]) for i in range(n)] +
             [F.mul(r_, row[i]) ^ F.mul(s_, row[n + i]) for i in range(n)] for row in R]
        R2, piv2 = rref_gf8(M, F)
        if piv2 != list(range(n)):
            continue
        n_graph += 1
        G = [row[n:] for row in R2]
        assert all(G[i][j] == G[j][i] for i in range(n) for j in range(n)), "G symmetric"
        assert all(G[i][j] == G[0][(j - i) % n] for i in range(n) for j in range(n)), "G circulant"
        c = tuple(G[0][d] for d in range(1, h + 1))  # the constant diagonal G_00 is removed by (a|b) -> (a|b + G_00 a)
        least, orb = canon_orbit(c, n, F)
        reached.setdefault(least, {"T": [p_, q_, r_, s_], "c": list(c), "diagonal": G[0][0], "orbit_size": len(orb)})
    return reached, n_T, n_graph


def graph_rows(c, n):
    return [[1 if j == i else 0 for j in range(n)] +
            [0 if i == j else c[min((j - i) % n, (i - j) % n) - 1] for j in range(n)] for i in range(n)]


def scan_classes(n, F, rep_canon):
    """classes of the scan orbits under the uniform local Cliffords: connected components of 'reached'"""
    parent = list(range(len(rep_canon)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, rep in enumerate(rep_canon):
        reached, _, _ = graph_forms(graph_rows(rep, n), n, F)
        assert rep in reached and all(r in rep_canon for r in reached), "closed inside the scan"
        for r in reached:
            parent[find(rep_canon.index(r))] = find(i)
    cls = {}
    for i in range(len(rep_canon)):
        cls.setdefault(find(i), []).append(i + 1)
    return sorted(cls.values())


def explicit_map_to_theorem1(R, info, n, F):
    """For the positive case: find a multiplier u, a Frobenius power k and a scalar mu such that
    T, then (a|b) -> (a|b + g a), then party i -> u i, then entrywise x -> x^(2^k), then (a|b) -> (mu a | mu^-1 b)
    maps every basis row into {(x | x Gamma)} with Gamma the circulant of Theorem 1.  Checked row by row."""
    p_, q_, r_, s_ = info["T"]
    g = info["diagonal"]
    target = graph_rows(THEOREM1[n], n)
    Gt = [row[n:] for row in target]
    rows = []
    for row in R:
        a = [F.mul(p_, row[i]) ^ F.mul(q_, row[n + i]) for i in range(n)]
        b = [F.mul(r_, row[i]) ^ F.mul(s_, row[n + i]) for i in range(n)]
        b = [y ^ F.mul(g, x) for x, y in zip(a, b)]
        rows.append((a, b))
    for u in range(1, n):
        if math.gcd(u, n) != 1:
            continue
        for k in range(3):
            for mu in range(1, 8):
                ok = True
                for a, b in rows:
                    a2, b2 = [0] * n, [0] * n
                    for i in range(n):
                        x, y = a[i], b[i]
                        for _ in range(k):
                            x, y = F.mul(x, x), F.mul(y, y)
                        a2[(u * i) % n], b2[(u * i) % n] = F.mul(mu, x), F.mul(F.inv(mu), y)
                    img = [0] * n
                    for i in range(n):
                        if a2[i]:
                            for j in range(n):
                                img[j] ^= F.mul(a2[i], Gt[i][j])
                    if img != b2:
                        ok = False
                        break
                if ok:
                    return {"T": info["T"], "diagonal_removed": g, "multiplier": u, "frobenius_power": k, "mu": mu}
    return None


def analyse(path):
    pg = parse_page(path)
    n = pg["n"]
    out = {"page": path.name, "cell": f"[[{n},0]]_8", "table_parameter_q": pg["q"], "lower": pg["lower"], "upper": pg["upper"],
           "construction_type": pg["type"], "last_modified": pg["modified"], "generating_polynomial": pg["poly"],
           "rows_printed": len(pg["rows"])}
    good = []
    for modulus in (0b1011, 0b1101):
        F = GF8(modulus)
        rows = [([entry(t, F) for t in x], [entry(t, F) for t in z]) for x, z in pg["rows"]]
        piv = f2_span_basis(pack(x + z) for x, z in rows)
        iso = all((sum(F.tr(F.mul(a, d) ^ F.mul(b, c)) for a, b, c, d in zip(x1, z1, x2, z2)) % 2) == 0
                  for (x1, z1), (x2, z2) in itertools.combinations(rows, 2))
        w = F.exp[1]
        lin = all(in_span(piv, pack([F.mul(w, a) for a in x + z])) for x, z in rows)
        info = {"modulus": "w^3+w+1" if modulus == 0b1011 else "w^3+w^2+1", "f2_dim": len(piv), "isotropic": iso,
                "f8_linear": lin}
        out.setdefault("field_test", []).append(info)
        if len(piv) == 3 * n and iso:
            good.append((modulus, F, rows))
    assert len(good) == 1, "the data should single out the minimal polynomial of w"
    modulus, F, rows = good[0]
    out["w_minimal_polynomial"] = "w^3+w+1" if modulus == 0b1011 else "w^3+w^2+1"
    lin = [i for i in out["field_test"] if i["modulus"] == out["w_minimal_polynomial"]][0]["f8_linear"]
    out["f8_linear"] = lin
    if not lin:
        out["conclusion"] = "not F_8-linear in the printed coordinates; no comparison made"
        return out
    assert modulus == 0b1011, "w = alpha is assumed below (the data single out w^3 = w + 1)"
    R, pivots = rref_gf8([x + z for x, z in rows], F)
    assert len(R) == n, "F_8-dimension"
    out["x_part_information_set_as_printed"] = pivots == list(range(n))
    h = n // 2
    hits = json.loads(HITS.read_text())[str(n)]
    reps = hits.get("orbit_reps", hits.get("orbit_representatives"))
    rep_canon = sorted(set(canon_orbit(tuple(r), n, F)[0] for r in reps))
    hits["orbits"] = hits.get("orbits", len(reps))
    assert len(rep_canon) == hits["orbits"] == len(reps), "orbit representatives of the scan"
    thm, orb_t = canon_orbit(THEOREM1[n], n, F)
    assert thm in rep_canon
    reached, n_T, n_graph = graph_forms(R, n, F)
    assert n_T == 504
    out["uniform_local_cliffords_tried"] = n_T
    out["of_which_give_graph_form"] = n_graph
    out["scan_orbits"] = hits["orbits"]
    lst = []
    for least in sorted(reached):
        info = reached[least]
        c = tuple(info["c"])
        Gm = [[0 if i == j else c[min((j - i) % n, (i - j) % n) - 1] for j in range(n)] for i in range(n)]
        bad = total = 0
        for A in itertools.combinations(range(n), h):
            Ac = [j for j in range(n) if j not in A]
            total += 1
            if rank_gf8([[Gm[i][j] for j in Ac] for i in A], F) != h:
                bad += 1
        lst.append({"orbit_least": list(least), "scan_orbit_number": rep_canon.index(least) + 1 if least in rep_canon else None,
                    "example_T": info["T"], "c": info["c"], "c_powers": ["0" if x == 0 else f"w^{F.log[x]}" for x in c],
                    "diagonal_removed": info["diagonal"], "orbit_size": info["orbit_size"],
                    "cuts_checked": total, "bad_cuts": bad})
    out["graph_forms"] = lst
    out["is_AME"] = bool(lst) and all(e["bad_cuts"] == 0 for e in lst)
    out["all_graph_forms_in_scan"] = all(e["scan_orbit_number"] is not None for e in lst)
    out["orbit_least_theorem1"] = list(thm)
    out["theorem1_scan_orbit_number"] = rep_canon.index(thm) + 1
    out["equivalent_to_theorem1_by_these_maps"] = thm in reached
    classes = scan_classes(n, F, rep_canon)
    out["scan_orbit_classes_under_uniform_local_cliffords"] = classes
    out["number_of_classes"] = len(classes)
    if thm in reached:
        out["explicit_map_to_theorem1"] = explicit_map_to_theorem1(R, reached[thm], n, F)
        assert out["explicit_map_to_theorem1"] is not None, "explicit map"
    if thm in reached:
        out["conclusion"] = (
            f"the code shown on codetables.de is F_8-linear and cyclic; a local Clifford (the same T = {reached[thm]['T']} in "
            "SL(2,8) on every party, then removal of a constant diagonal) brings it to a symmetric circulant graph code whose "
            "coefficient vector lies in the orbit (multipliers, Frobenius, scaling) of the vector printed in Theorem 1; the composite map (explicit_map_to_theorem1) was checked row by row. "
            f"Hence it is equal to the AME({n},8) state of Theorem 1 up to a party permutation and local Clifford unitaries.")
    else:
        out["conclusion"] = (
            f"the code shown on codetables.de is F_8-linear and cyclic, and AME({n},8); the uniform local Cliffords bring it to "
            f"symmetric circulant graph codes in scan orbit(s) {[e['scan_orbit_number'] for e in lst]}, none of which is the "
            f"orbit of the vector printed in Theorem 1 (scan orbit {rep_canon.index(thm) + 1}). Whether the two states are "
            "equivalent under general (non-uniform) local unitaries and permutations is not decided by this test "
            "(see lc_invariant.py, which separates the classes).")
    return out


def summary(res):
    lines = []
    for c in res["cells"]:
        lines.append(f"{c['cell']}: codetables.de lower {c['lower']} upper {c['upper']}, construction type {c['construction_type']}, "
                     f"last modified {c['last_modified']}")
        lines.append(f"  printed matrix: {c['rows_printed']} rows; w has minimal polynomial {c['w_minimal_polynomial']}; "
                     f"F_8-linear {c['f8_linear']}")
        lines.append(f"  uniform local Cliffords giving a graph form: {c['of_which_give_graph_form']} of "
                     f"{c['uniform_local_cliffords_tried']}; scan orbits reached: "
                     f"{[e['scan_orbit_number'] for e in c['graph_forms']]} of {c['scan_orbits']}; "
                     f"cuts checked per form {c['graph_forms'][0]['cuts_checked']}, bad cuts {sum(e['bad_cuts'] for e in c['graph_forms'])}")
        lines.append(f"  AME: {c['is_AME']}; Theorem 1 is scan orbit {c['theorem1_scan_orbit_number']}; "
                     f"equivalent to Theorem 1 by these maps: {c['equivalent_to_theorem1_by_these_maps']}"
                     + (f"; explicit map {c['explicit_map_to_theorem1']}" if c.get("explicit_map_to_theorem1") else ""))
        lines.append(f"  classes of the {c['scan_orbits']} scan orbits under uniform local Cliffords: {c['number_of_classes']} "
                     f"{c['scan_orbit_classes_under_uniform_local_cliffords']}")
    return lines


def main(argv):
    """usage: compare_codetables.py [pages_dir hits.json [out.json]]"""
    global SRC, HITS
    out = P / "results" / "codetables_comparison_2026-11-10.json"
    if len(argv) >= 3:
        SRC, HITS = pathlib.Path(argv[1]), pathlib.Path(argv[2])
        out = pathlib.Path(argv[3]) if len(argv) > 3 else None
    res = {"fetched": "2026-11-10", "cells": [analyse(SRC / f"qecc_q64_n{n}_k0.html") for n in (15, 17)]}
    if out is not None:
        out.write_text(json.dumps(res, indent=1) + "\n")
    for line in summary(res):
        print(line)
    ok = all(c.get("is_AME") and c.get("all_graph_forms_in_scan") for c in res["cells"])
    print("COMPARISON DONE" if ok else "COMPARISON INCOMPLETE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
