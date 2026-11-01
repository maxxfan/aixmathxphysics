#!/usr/bin/env python3
"""Provenance of the field-level JSON files: re-derive every field-level object FROM THE CERTIFIED canonical rows
(../ingredients/<id>.json) by inverting the documented converter of convert.py, identify the structure stated in the
construction note (circulant first rows, bordered-circulant parameters, evaluation points, ...), and write
fieldlevel/<id>.json.  With --check nothing is written; the re-derived objects are compared with the stored JSON files.

Inversion (conventions of convert.py / qlib):
  canonical row (x | z) over F_p, party j digits x_{j,0..m-1}, z_{j,0..m-1}
     -> a_j = sum_i x_{j,i} beta_i,  b_j = sum_i z_{j,i} gamma_i   (beta basis of F_q/F_p, gamma its trace dual)
  Hermitian: v_j = w a_j + w^q b_j in F_{q^2} (w the normal element); the F_{q^2}-span of the v's is C; RREF -> [I | A].
  CSS (sum_{c in C}|c>): the X parts a span C, the Z parts b span C^perp; RREF of the a's -> systematic generator.
  graph: RREF over F_p of the rows is [I | Gamma].
Requires Python 3.9 stdlib and convert.py (same directory).  usage: python3 derive.py [--check]
"""
import json
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import convert as C  # noqa: E402

ING = HERE.parent / "ingredients"


# ----------------------------------------------------------------------------------------------------------------
# linear algebra over GF(Q) (convert.GF)
# ----------------------------------------------------------------------------------------------------------------
def rref(rows, F):
    M = [r[:] for r in rows]
    piv, r = [], 0
    ncol = len(M[0])
    for c in range(ncol):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = F.inv(M[r][c])
        M[r] = [F.mul(iv, x) for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    return M[:r], piv


def to_fq_pairs(rows, F, q, n, m):
    """canonical F_p rows -> list of (a, b) over F_q (elements of F)"""
    beta, gamma = C.subfield_bases(F, q)
    N = n * m
    out = []
    for r in rows:
        a, b = [], []
        for j in range(n):
            x = z = 0
            for i in range(m):
                x = F.add(x, F.mul(r[j * m + i] % F.p, beta[i]))
                z = F.add(z, F.mul(r[N + j * m + i] % F.p, gamma[i]))
            a.append(x)
            b.append(z)
        out.append((a, b))
    return out, beta, gamma


def recover_hermitian(ing, Q):
    F = C.field(Q)
    q = ing["d"]
    p, m = C._factor_prime_power(q)
    n = ing["n"]
    rows = C.parse_rows(ing)
    pairs, beta, gamma = to_fq_pairs(rows, F, q, n, m)
    w = C.normal_element(F, q)
    wq = F.pow(w, q)
    V = [[F.add(F.mul(w, x), F.mul(wq, y)) for x, y in zip(a, b)] for a, b in pairs]
    R, piv = rref(V, F)
    assert len(R) == n // 2, ("F_{q^2}-rank", len(R))
    assert piv == list(range(n // 2)), ("first n/2 columns are not an information set", piv)
    A = [r[n // 2:] for r in R]
    xi = next(x for x in range(F.q) if x not in set(F.subfield(q)))
    conv = {"function": "canonical_from_hermitian", "subfield_basis_beta": beta, "trace_dual_basis_gamma": gamma,
            "normal_element_w": w, "w_to_the_q": wq, "xi": xi}
    return F, A, conv


def recover_css(ing, q):
    F = C.field(q)
    n = ing["n"]
    rows = C.parse_rows(ing)
    pairs, beta, gamma = to_fq_pairs(rows, F, q, n, F.m)
    RX, pivx = rref([a for a, b in pairs], F)
    RZ, pivz = rref([b for a, b in pairs], F)
    # X part spans C, Z part spans C^perp: check orthogonality and dimensions
    assert len(RX) + len(RZ) == n, (len(RX), len(RZ))
    assert all(F.dot(u, v) == 0 for u in RX for v in RZ)
    # and the space is the direct sum: X-only and Z-only rows are in it (dimension count over F_p)
    conv = {"function": "mds_state_rows", "subfield_basis_beta": beta, "trace_dual_basis_gamma": gamma}
    return F, RX, pivx, conv


def recover_graph(ing):
    p = ing["p"]
    rows = C.parse_rows(ing)
    R, piv = C.rref_mod_p(rows, p)
    n = ing["n"]
    assert piv[:n] == list(range(n)) and len(R) == n
    Gam = [r[n:] for r in R]
    assert all(Gam[i][j] == Gam[j][i] for i in range(n) for j in range(n)) and all(Gam[i][i] == 0 for i in range(n))
    return Gam


# ----------------------------------------------------------------------------------------------------------------
# structure identification
# ----------------------------------------------------------------------------------------------------------------
def is_circulant(A):
    n = len(A)
    return all(A[i][j] == A[0][(j - i) % n] for i in range(n) for j in range(n))


def bordered(A):
    """A = [[alpha, beta..beta], [gamma, core circulant]]: returns params or None"""
    k = len(A)
    alpha, beta, gamma = A[0][0], A[0][1], A[1][0]
    if any(A[0][j] != beta for j in range(1, k)) or any(A[i][0] != gamma for i in range(1, k)):
        return None
    core = [r[1:] for r in A[1:]]
    if not is_circulant(core):
        return None
    return {"type": "bordered_circulant", "c": k - 1, "corner_alpha": alpha, "row_border_beta": beta,
            "column_border_gamma": gamma, "core_first_row": core[0],
            "rule": "A[0][0] = alpha; A[0][j] = beta and A[j][0] = gamma (1 <= j <= c); "
                    "A[1+i][1+j] = core_first_row[(j - i) mod c] (0 <= i, j < c)"}


def block_circulant(A, b):
    k = len(A)
    nb = k // b
    blocks = []
    for r in range(nb):
        row = []
        for s in range(nb):
            Bk = [[A[b * r + x][b * s + y] for y in range(b)] for x in range(b)]
            if not is_circulant(Bk):
                return None
            row.append(Bk[0])
        blocks.append(row)
    return {"type": "block_circulant", "block_size": b, "blocks_first_rows": blocks,
            "rule": "A[b*r + x][b*s + y] = blocks_first_rows[r][s][(y - x) mod b]  (b = block_size)"}


def dlog(F, g, x):
    y, e = 1, 0
    while y != x:
        y = F.mul(y, g)
        e += 1
        if e > F.q:
            raise ValueError
    return e


def hz_kernel(F, q):
    """Hontarenko-Zyczkowski kernel A_q (q x q) exactly as work/mds/mdslib.hz_matrix (rows {0} u F_q^*, columns
    {inf} u F_q^*), with g = x (the primitive element, int p) and gam = g^(q+1)."""
    g = F.exp[1]
    gam = F.pow(g, q + 1)
    two = 2 % F.p
    j = next(j for j in range(q - 1) if F.pow(gam, j) == two)
    b = F.pow(g, j)
    bb = F.pow(b, q)
    xs = [F.pow(gam, i) for i in range(q - 1)]
    K = [[0] * q for _ in range(q)]
    K[0][0] = F.neg(F.mul(b, F.inv(bb)))
    for c in range(q - 1):
        K[0][c + 1] = b
        K[c + 1][0] = b
    for a, x in enumerate(xs):
        for c, y in enumerate(xs):
            num = F.mul(F.mul(two, y), F.add(y, F.mul(g, x)))
            den = F.sub(F.mul(gam, F.mul(x, x)), F.mul(y, y))
            K[a + 1][c + 1] = F.mul(num, F.inv(den))
    return K, {"g": g, "gam": gam, "b": b, "j_with_gam^j=2": j}


def field_block(F, q=None):
    fd = {"Q": F.q, "p": F.p, "m": F.m, "modulus": F.modulus()}
    if F.m == 1:
        fd["modulus_str"] = None
        fd["encoding"] = "prime field F_%d: the ints 0..%d with arithmetic mod %d" % (F.p, F.p - 1, F.p)
    else:
        terms = []
        for i, c in reversed(list(enumerate(F.modulus()))):
            if c == 0:
                continue
            mon = "1" if i == 0 else ("x" if i == 1 else "x^%d" % i)
            terms.append(mon if (c == 1 and i > 0) else ("%d" % c if i == 0 else "%d*%s" % (c, mon)))
        fd["modulus_str"] = " + ".join(terms)
        fd["encoding"] = ("F_%d = F_%d[x]/(%s); the int e = sum_{i<%d} c_i %d^i (0 <= c_i < %d) is the element "
                          "sum_i c_i x^i" % (F.q, F.p, fd["modulus_str"], F.m, F.p, F.p))
        fd["x_is_primitive"] = True
        fd["modulus_choice"] = "lexicographically first primitive monic polynomial of degree m over F_p (qlib GF)"
    if q is not None:
        fd["q"] = q
        fd["conjugation"] = "z -> z^%d" % q
    return fd


# ----------------------------------------------------------------------------------------------------------------
# per-ingredient specification (construction notes, sources) and derivation
# ----------------------------------------------------------------------------------------------------------------
BOARD = "board/board.db posts"


def derive(iid):
    ing = json.loads((ING / (iid + ".json")).read_text())
    n, d = ing["n"], ing["d"]
    out = {"id": iid, "n": n, "local_dim": d}
    notes = []
    k_unif = n // 2
    src = {"script": ing["script"], "note": ing["note"]}

    if iid in ("a_wheel_ame6_2", "a_dodecacode_n12_d2"):
        Gam = recover_graph(ing)
        F = C.field(d)
        out["kind"] = "fq_graph"
        out["field"] = field_block(F)
        out["Gamma"] = Gam
        if iid == "a_wheel_ame6_2":
            W = [[0] * 6 for _ in range(6)]
            for i in range(5):
                j = (i + 1) % 5
                W[i][j] = W[j][i] = 1
                W[i][5] = W[5][i] = 1
            assert Gam == W, "not the wheel"
            out["structure"] = {"type": "wheel", "rule": "5-cycle 0-1-2-3-4-0 plus hub 5 adjacent to all: "
                                "Gamma[i][i+1 mod 5] = Gamma[i][5] = 1 (symmetric), else 0"}
            out["claimed_property"] = {"statement": "AME(6,2)", "uniformity_k": 3}
            out["stabilizer_code"] = "[[6,0,4]]_2"
            out["code"] = {"params": "[[6,0,4]]_2", "distance": 4}
            src.update({"builder": "work/couple/lib.py wheel6()", "board_posts": [586]})
            out["used_as"] = "G_a (k = 2, needs 2-uniform; it is AME(6,2), i.e. 3-uniform)"
        else:
            assert is_circulant(Gam)
            z = 0x2e8
            fr = [(z >> i) & 1 for i in range(12)]
            assert Gam[0] == fr, (Gam[0], fr)
            out["structure"] = {"type": "circulant", "first_row": Gam[0],
                                "connection_set": [i for i in range(12) if fr[i]],
                                "rule": "Gamma[i][j] = first_row[(j - i) mod 12] (circulant graph on Z_12)"}
            # GF(4) additive cyclic description: X -> 1, Z -> omega (omega^2 = omega + 1), Y -> omega^2
            out["additive_gf4"] = {
                "description": "the dodecacode: additive cyclic (12, 2^12, 6) self-dual code over GF(4) = {0, 1, w, w^2} "
                               "(w^2 = w + 1). Pauli map X(a)Z(b) -> a + b w per qubit. Generated over F_2 by the 12 "
                               "cyclic shifts of one word (the shifts are the 12 graph-state stabilizers X_i Z_{N(i)}).",
                "cyclic_generator_pauli": "X_0 Z_3 Z_5 Z_6 Z_7 Z_9  (x = 0x001, z = 0x2e8 as 12-bit masks, bit i = qubit i)",
                "cyclic_generator_gf4": " ".join("1" if i == 0 else ("w" if fr[i] else "0") for i in range(12))}
            out["claimed_property"] = {"statement": "5-uniform", "uniformity_k": 5}
            out["stabilizer_code"] = "[[12,0,6]]_2"
            out["code"] = {"params": "[[12,0,6]]_2 (additive (12, 2^12, 6) over GF(4))", "distance": 6}
            src.update({"builder": "work/couple/lib.py dodeca()", "board_posts": [586]})
            out["used_as"] = "G_a (k = 5, needs 5-uniform)"
            k_unif = 5
        out["converter"] = {"function": "canonical_from_graph", "rule": "rows [e_i | Gamma_i] over F_p"}

    elif iid == "a_ternary_golay_n12_d3" or iid.endswith("_rs"):
        F, RX, pivx, conv = recover_css(ing, d)
        out["kind"] = "css_minimal_support"
        out["field"] = field_block(F)
        out["state"] = "|C> = |C|^(-1/2) sum_{c in C} |c>, stabilized by X(c) (c in C) and Z(h) (h in C^perp, Euclidean dual)"
        if iid == "a_ternary_golay_n12_d3":
            Ag = [[0, 1, 1, 1, 1, 1], [1, 0, 1, 2, 2, 1], [1, 1, 0, 1, 2, 2], [1, 2, 1, 0, 1, 2], [1, 2, 2, 1, 0, 1],
                  [1, 1, 2, 2, 1, 0]]
            G = [[1 if i == j else 0 for j in range(6)] + Ag[i] for i in range(6)]
            assert rref(G, F)[0] == RX, "not the stated Golay generator"
            out["G"] = G
            st = bordered(Ag)
            assert st is not None
            st["type"] = "systematic_bordered_circulant"
            st["rule"] = st["rule"].replace("A[", "P[")
            st["description"] = "G = [I_6 | P], P = bordered circulant (Paley matrix, quadratic residues mod 5)"
            st["core_rule"] = "core_first_row[j] = 0 (j = 0), 1 (j a nonzero square mod 5), 2 (j a non-square)"
            out["structure"] = st
            out["claimed_property"] = {"statement": "5-uniform", "uniformity_k": 5}
            out["code"] = {"params": "[12,6,6]_3", "field_order": 3, "length": 12, "dimension": 6, "distance": 6,
                           "duality": "Euclidean self-dual (extended ternary Golay code), so C^perp = C"}
            out["stabilizer_code"] = "[[12,0,6]]_3"
            src.update({"builder": "work/couple/lib.py golay12()", "board_posts": [586]})
            out["used_as"] = "G_a (k = 5, needs 5-uniform)"
            k_unif = 5
        else:
            G = C.rs_extended(F, n, n // 2)
            assert rref(G, F)[0] == RX, "not the stated RS code"
            out["G"] = G
            pts = list(range(min(n, F.q)))
            out["structure"] = {"type": "doubly_extended_reed_solomon" if n == F.q + 1 else "reed_solomon",
                                "k": n // 2, "evaluation_points": pts, "point_at_infinity": n == F.q + 1,
                                "rule": "G[i][j] = x_j^i for the evaluation point x_j (field ints; 0^0 = 1), i < k"
                                        + ("; last column = e_{k-1} (point at infinity)" if n == F.q + 1 else "")}
            out["G_systematic"] = RX
            out["claimed_property"] = {"statement": "AME(%d,%d)" % (n, d), "uniformity_k": n // 2}
            out["code"] = {"params": "[%d,%d,%d]_%d" % (n, n // 2, n // 2 + 1, d), "field_order": d, "length": n,
                           "dimension": n // 2, "distance": n // 2 + 1, "duality": "C^perp is again MDS [%d,%d,%d]" % (n, n // 2, n // 2 + 1)}
            out["stabilizer_code"] = "[[%d,0,%d]]_%d" % (n, n // 2 + 1, d)
            src.update({"builder": "work/couple/lib.py ame_rs() -> qlib rs_extended, mds_state_rows",
                        "board_posts": [586] if n == 12 else ([597] if n == 14 else [586])})
            out["used_as"] = "G_b (AME(%d,%d))" % (n, d)
        out["converter"] = conv

    else:  # Hermitian
        Q = d * d
        F, A, conv = recover_hermitian(ing, Q)
        out["kind"] = "hermitian_systematic"
        out["field"] = field_block(F, d)
        out["q"] = d
        out["generator"] = "[I_%d | A] over F_%d, rows = basis of the Hermitian self-dual code C" % (n // 2, Q)
        out["A"] = A
        out["converter"] = conv
        dist = n // 2 + 1
        st = None
        if iid == "a_14_0_7_d4_hermitian_circulant":
            assert is_circulant(A) and A[0] == [0, 1, 2, 15, 6, 2, 7], A[0]
            st = {"type": "circulant", "description": "7 x 7 circulant over F_16", "first_row": A[0],
                  "rule": "A[i][j] = first_row[(j - i) mod 7]",
                  "search": "near14_4.py: 20475 unitary 7 x 7 circulants over F_16 (A conj(A)^T = I), 840 with d = 7; "
                            "this is hit 0 (work/couple/near14_4_result.json)"}
            dist, k_unif = 7, 6
            src.update({"builder": "work/couple/g14_4.py 0 (hit 0 of work/couple/near14_4.py)", "board_posts": [597]})
        elif iid == "a_16_0_8_d4_bordered_circulant":
            st = bordered(A)
            assert st is not None
            got = [st["corner_alpha"], st["row_border_beta"], st["column_border_gamma"]] + st["core_first_row"]
            assert got == [2, 4, 4, 2, 1, 1, 0, 2, 8, 5], got
            dist, k_unif = 8, 7
            st["description"] = "bordered Z_7 circulant over F_16 (row border beta, column border gamma)"
            st["search"] = ("near_bord.py 4 7 8: 1228500 unitary bordered circulants (A conj(A)^T = -I), first 20 hits "
                            "with d >= 8 stored in work/couple/nearb_q4_c7_d8.json; this is hit 0")
            notes.append("The note's list [2, 4, 4, 2, 1, 1, 0, 2, 8, 5] is [alpha, beta, gamma, core_first_row] "
                         "(near_bord.py hit format).")
            src.update({"builder": "work/couple/g16_4.py 0 (hit 0 of work/couple/near_bord.py 4 7 8, "
                                   "work/couple/nearb_q4_c7_d8.json)", "board_posts": [616, 644, 653]})
        elif iid == "b_ame12_5_bevins_bidav":
            A5 = [[1, 1, 1, 1, 1, 2], [23, 1, 2, 3, 5, 3], [19, 7, 1, 22, 17, 12], [9, 15, 11, 2, 18, 15],
                  [1, 4, 22, 9, 2, 1], [6, 7, 18, 15, 10, 2]]
            # Bevins-Bidav labels e = e0 + 5 e1 <-> e0 + e1*r, r a root of x^2 + x + 2; here r = x (int 5), so labels coincide
            assert F.add(F.add(F.mul(5, 5), 5), 2) == 0
            assert A == A5, A
            st = {"type": "explicit_labels", "description": "Bevins-Bidav A_5 (arXiv:2608.05781), the printed 6 x 6 matrix",
                  "source_labels": A5, "r": 5, "r_minimal_polynomial": [2, 1, 1],
                  "rule": "a printed label e (0 <= e < 25) is the element (e mod 5) + (e div 5) r, r a root of the "
                          "printed minimal polynomial x^2 + x + 2 (transcription: work/couple/lib.py ame12_5, "
                          "work/auditor/bb_kernel.py); r = x = int 5 here, because x^2 + x + 2 is also the modulus, "
                          "so the printed labels and the ints of this file coincide"}
            src.update({"builder": "work/couple/lib.py ame12_5() (data in results/literature.md)", "board_posts": [597]})
        elif iid == "b_ame12_7_circulant":
            assert is_circulant(A)
            st = {"type": "circulant", "description": "Z_6 group circulant", "first_row": A[0],
                  "rule": "A[i][j] = first_row[(j - i) mod 6]"}
            src.update({"builder": "work/couple/run_gorb.py 7 6 0 -> work/mds/grp_orbits3.py (first hit, "
                                   "work/couple/gorb/grp3_q7_Z6_b0_n12.json)", "board_posts": [597]})
        elif iid == "b_ame12_8_circulant":
            st = block_circulant(A, 3)
            assert st is not None
            st["description"] = "Z_3-block circulant: 2 x 2 blocks, each a 3 x 3 circulant"
            src.update({"builder": "work/couple/ame12_8.py 1 (random Fourier sampling of unitary Z_3-block circulants)",
                        "board_posts": [610, 616]})
        elif iid in ("b_ame12_9_circulant", "b_ame14_5_circulant", "b_ame14_7_circulant", "b_ame14_11_circulant",
                     "b_ame16_13_circulant"):
            st = bordered(A)
            assert st is not None and st["row_border_beta"] == st["column_border_gamma"]
            c = st["c"]
            # grp_orbits3 border: alpha = -conj(F0) b / conj(b), F0 = sum of core first row
            F0 = 0
            for x in st["core_first_row"]:
                F0 = F.add(F0, x)
            b = st["row_border_beta"]
            alpha = F.neg(F.mul(F.mul(F.pow(F0, d), b), F.inv(F.pow(b, d))))
            assert alpha == st["corner_alpha"]
            st["description"] = "bordered Z_%d group circulant, symmetric border (beta = gamma)" % c
            st["corner_relation"] = "alpha = -conj(F0) * beta / conj(beta), F0 = sum(core_first_row) (holds)"
            gq = {9: "9 5 1", 5: "5 6 1", 7: "7 6 1", 11: "11 6 1", 13: "13 7 1"}[d]
            fn = {9: "grp3_q9_Z5_b1_n12.json", 5: "grp3_q5_Z6_b1_n14.json", 7: "grp3_q7_Z6_b1_n14.json",
                  11: "grp3_q11_Z6_b1_n14.json", 13: "grp3_q13_Z7_b1_n16.json"}[d]
            src.update({"builder": "work/couple/run_gorb.py %s -> work/mds/grp_orbits3.py (first hit, work/couple/gorb/%s)"
                                   % (gq, fn), "board_posts": [597] if d != 13 else [614, 616]})
        elif iid == "b_ame16_11_circulant":
            st = bordered(A)
            assert st is not None and st["row_border_beta"] == st["column_border_gamma"] == 1
            st["description"] = "bordered Z_7 circulant, border beta = gamma = 1"
            st["search"] = ("f121z7.py: the core is diagonalised over F_{11^6} (7 does not divide 120); hit 0 = "
                            "(m = 0, e1 = 16469); 84 of the 1.6e7 unitary family members are MDS")
            src.update({"builder": "work/couple/ame16_11.py 0 (hit 0 of work/couple/f121z7.py)", "board_posts": [653]})
        elif iid == "b_ame16_9_hz":
            st = bordered(A)
            assert st is not None
            g = 3  # x, a root of x^4 + x + 2 of order 80 (the modulus of F_81 here)
            assert F.modulus() == [2, 1, 0, 0, 1] and F.exp[1] == g
            exps = [dlog(F, g, x) for x in st["core_first_row"]]
            assert st["corner_alpha"] == st["row_border_beta"] == st["column_border_gamma"] == 1
            assert exps == [29, 21, 37, 2, 40, 57, 76], exps
            st["description"] = "Hontarenko-Zyczkowski parent M_{16,9} = B_7(1, 1; [29,21,37,2,40,57,76]_g)"
            st["g"] = g
            st["core_exponents"] = exps
            st["exponent_rule"] = "core_first_row[t] = g^core_exponents[t], g = x = int 3 (root of x^4 + x + 2, order 80)"
            src.update({"builder": "work/couple/ame16_9.py (HZ parent, arXiv:2609.11796 supplement)",
                        "board_posts": [637, 644], "transcription": "work/auditor/hz_parent_16_9.txt (auditor, #637)"})
        elif iid == "b_ame14_9_hz":
            K, kp = hz_kernel(F, d)
            Rr = Cc = [0, 1, 2, 3, 5, 6, 7]
            P = [[K[i][j] for j in Cc] for i in Rr]
            assert all(P[i][j] for i in range(7) for j in range(7))
            # factor A = diag(lam) P diag(D) with the normalisation D_0 = 1
            ratio = [[F.mul(A[i][j], F.inv(P[i][j])) for j in range(7)] for i in range(7)]
            lam = [ratio[i][0] for i in range(7)]
            D = [F.mul(ratio[0][j], F.inv(lam[0])) for j in range(7)]
            assert all(ratio[i][j] == F.mul(lam[i], D[j]) for i in range(7) for j in range(7)), "not Lambda P D"
            st = {"type": "hz_kernel_submatrix",
                  "description": "rescaled 7 x 7 submatrix of the Hontarenko-Zyczkowski kernel A_9: "
                                 "A = diag(Lambda) P diag(D), P = K[rows][cols]",
                  "kernel": "K = A_9 (9 x 9), rows {0} u F_9^*, columns {inf} u F_9^*: K[0][0] = -b/conj(b), "
                            "K[0][c+1] = K[c+1][0] = b, K[a+1][c+1] = 2y(y + g x)/(gam x^2 - y^2) with x = gam^a, "
                            "y = gam^c, g = x (int 3, primitive), gam = g^(q+1) = g^10, b = g^j with gam^j = 2 "
                            "(work/mds/mdslib.hz_matrix)",
                  "kernel_params": {"g": kp["g"], "gam": kp["gam"], "b": kp["b"], "j": kp["j_with_gam^j=2"]},
                  "kernel_matrix": K, "rows": Rr, "cols": Cc, "Lambda": lam, "D": D,
                  "normalisation": "D[0] = 1 (Lambda and D are unique up to Lambda -> c Lambda, D -> D / c)"}
            src.update({"builder": "work/couple/hzsub.py 9 7 (work/couple/ing/hz_q9_m7.json), loaded by batch13b.py",
                        "board_posts": [622, 624]})
        else:
            raise KeyError(iid)
        out["structure"] = st
        if iid.startswith("a_"):
            out["claimed_property"] = {"statement": "%d-uniform" % k_unif, "uniformity_k": k_unif}
            out["used_as"] = "G_a (k = %d, needs %d-uniform)" % (k_unif, k_unif)
        else:
            out["claimed_property"] = {"statement": "AME(%d,%d)" % (n, d), "uniformity_k": n // 2}
            out["used_as"] = "G_b (AME(%d,%d))" % (n, d)
        out["code"] = {"params": "[%d,%d,%d]_%d" % (n, n // 2, dist, Q), "field_order": Q, "length": n,
                       "dimension": n // 2, "distance": dist,
                       "duality": "Hermitian self-dual: <u,v> = sum_i u_i v_i^%d" % d + (" (MDS)" if dist == n // 2 + 1 else "")}
        out["stabilizer_code"] = "[[%d,0,%d]]_%d" % (n, dist, d)
    out["claimed_property"].update({"n": n, "d": d})
    out["used_in"] = ing["used_in"]
    out["source"] = src
    out["canonical_file"] = "../ingredients/%s.json" % iid
    if notes:
        out["remarks"] = notes
    return out


def dumps(obj):
    s = json.dumps(obj, indent=1)
    # put innermost integer lists on one line
    return re.sub(r"\[\n\s*(-?\d+(?:,\n\s*-?\d+)*)\n\s*\]",
                  lambda mo: "[" + ", ".join(x.strip() for x in mo.group(1).split(",")) + "]", s) + "\n"


def main(argv):
    check = "--check" in argv
    ids = sorted(f.stem for f in ING.glob("*.json"))
    bad = 0
    for iid in ids:
        obj = derive(iid)
        path = HERE / (iid + ".json")
        if check:
            same = path.exists() and json.loads(path.read_text()) == json.loads(dumps(obj))
            bad += not same
            print("%-34s %s" % (iid, "same" if same else "DIFFERENT"))
        else:
            path.write_text(dumps(obj))
            print("wrote", path.name, obj["kind"], obj.get("structure", {}).get("type", ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
