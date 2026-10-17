#!/usr/bin/env python3
"""Standalone checker for the GF(2^n) T-count release (Python 3 standard library only).

It shares no code with the project verifier.  For every record it checks:
  1. the modulus P is irreducible over GF(2) (gcd(x^(2^i) - x, P) = 1 for i <= n/2);
  2. the factor list satisfies the signature-tensor criterion: every variable and every pair of
     variables lies in an even number of factors, and a triple lies in an odd number exactly when
     it is a monomial a_i b_j c_k of f = c . (a*b mod P);
  3. the OpenQASM circuit, gate by gate, is EXACTLY the multiplier for every input.  Between the
     two Hadamard layers on c the circuit only uses cx, t, s, sdg, z, cz.  Each wire carries a
     linear form of the inputs (tracked through the cx gates); a phase gate with multiplicity m on
     a wire with form L contributes m * (L.x mod 2), which over Z_8 equals
     m * (sum_i x_i - 2 sum_{i<j} x_i x_j + 4 sum_{i<j<k} x_i x_j x_k) over the support of L.
     The checker sums these polynomials (and the cz terms) and requires the result to equal
     4 f(x) mod 8 coefficient by coefficient, and every wire to end with its own input.  This
     proves that the diagonal part is (-1)^f, so the circuit maps |a,b,c> to |a,b,c + ab>;
  4. as an independent cross-check, it simulates the circuit on basis states (bit-sliced): all
     2^(3n) inputs when 3n <= 24, otherwise 2048 random inputs.
It also checks the small-algebra blocks (criterion 2 for F2[y]/(m)) and the genus-2 construction
data in construction/curve64.json (curve, place Q, and the bilinear identity on random inputs).

usage: python3 -I check_release.py [RELEASE_DIR]
"""
import hashlib
import json
import pathlib
import random
import re
import sys

RNG = random.Random(2026)

# What the paper claims; the checker fails if any of these is missing or different.
EXPECTED_RECORDS = {5: 57, 6: 69, 7: 97, 8: 113, 10: 151, 11: 179, 12: 191, 13: 219, 14: 229, 15: 253, 16: 267,
                    32: 571, 64: 1255}
# (reported value, reported value invalid?, best checked published value) as printed in the paper, Table I
EXPECTED_COMPARISON = {5: (49, True, 59), 6: (77, False, 77), 7: (101, False, 104), 8: (123, False, 123),
                       10: (157, False, 157), 11: (215, True, None), 12: (243, True, None), 13: (293, True, None),
                       14: (319, True, None), 15: (363, True, None), 16: (383, True, None), 32: (1177, True, None),
                       64: (3729, True, None)}
# record bars of the frozen verifier (data/records.json, pinned in verifier/PINS.json) at release time
EXPECTED_BARS = {5: 59, 6: 77, 7: 101, 8: 123, 10: 157, 11: 215, 12: 243, 13: 293, 14: 319, 15: 363, 16: 383,
                 32: 1177, 64: 3729}
EXPECTED_BLOCKS = {  # algebra label -> (T-count, modulus pattern: exact bits, or 'irreducible:k')
    "GF(4)": (17, "irreducible:2"), "GF(8)": (29, "irreducible:3"), "GF(16)": (39, "irreducible:4"),
    "F2[y]/(y^2)": (15, "100"), "F2[y]/(y^3)": (25, "1000"),
    "F2[y]/((y^2+y+1)^2) = GF(4)[e]/(e^2)": (35, "10101")}


# ------------------------------------------------------------------ GF(2)[x]

def pmod(a, p):
    dp = p.bit_length() - 1
    while a and a.bit_length() - 1 >= dp:
        a ^= p << (a.bit_length() - 1 - dp)
    return a


def pmul(a, b, p):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a = pmod(a << 1, p)
    return r


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def is_irreducible(p):
    n = p.bit_length() - 1
    if n < 1:
        return False
    y = 2  # x
    for _ in range(n // 2):
        y = pmul(y, y, p)  # x^(2^i) mod p
        if pgcd(p, y ^ 2) != 1:
            return False
    return True


def bits(s):  # factor string, character i = variable i
    if not isinstance(s, str) or not re.fullmatch(r"[01]+", s):
        raise ValueError("factor must be a 0/1 string")
    return int(s[::-1], 2)


def poly_str(b):
    n = len(b) - 1
    terms = []
    for i, ch in enumerate(b):
        if ch == "1":
            e = n - i
            terms.append("1" if e == 0 else ("x" if e == 1 else f"x^{e}"))
    return " + ".join(terms)


def well_formed(r, dim):
    """Exact types and formats: ints are ints, modulus '1[01]*' of length dim+1, factors '[01]{3 dim}'."""
    ints = all(isinstance(r.get(k), int) and not isinstance(r.get(k), bool) for k in ("tcount",))
    mod = r.get("modulus")
    fac = r.get("factors")
    return (ints and isinstance(mod, str) and re.fullmatch(r"1[01]*", mod) is not None and len(mod) == dim + 1
            and isinstance(fac, list) and all(isinstance(x, str) and re.fullmatch("[01]{%d}" % (3 * dim), x) for x in fac)
            and r.get("modulus_poly") == poly_str(mod))


def target_rows(n, p):
    """rows[(i, j)] = mask of k > j such that x_i x_j x_k is a monomial of f (i < j < k)."""
    rows = {}
    for i in range(n):
        for j in range(n):
            v = pmod(1 << (i + j), p)
            for k in range(n):
                if v >> k & 1:
                    key = (i, n + j)
                    rows[key] = rows.get(key, 0) | (1 << (2 * n + k))
    return rows


# ------------------------------------------------------------------ check 2: tensor criterion

def tensor_ok(n, p, factors):
    N = 3 * n
    single = [0] * N
    pair = {}
    trip = {}
    for u in factors:
        s = [i for i in range(N) if u >> i & 1]
        for x in range(len(s)):
            single[s[x]] ^= 1
            for y in range(x + 1, len(s)):
                k = (s[x], s[y])
                pair[k] = pair.get(k, 0) ^ 1
                above = u & ~((1 << (s[y] + 1)) - 1)
                if above:
                    trip[k] = trip.get(k, 0) ^ above
    if any(single) or any(pair.values()):
        return False
    want = target_rows(n, p)
    keys = set(trip) | set(want)
    return all(trip.get(k, 0) == want.get(k, 0) for k in keys)


# ------------------------------------------------------------------ check 3: exact circuit identity

STMT = re.compile(r"^(h|cx|t|tdg|s|sdg|z|cz) q\[(0|[1-9][0-9]*)\](?:,q\[(0|[1-9][0-9]*)\])?;$")
HEX64 = re.compile(r"^([0-9a-f]{64})  ([A-Za-z0-9_./^-]+)$")
BAD_CHARS = re.compile(r"[^\x20-\x7e\n]")


def parse_qasm(text, n):
    """Strict reader: the first three statements must be exactly the standard header; every other
    non-comment line must be ONE statement of the form 'g q[i];' or 'g q[i],q[j];'."""
    N = 3 * n
    if BAD_CHARS.search(text):
        raise ValueError("file contains characters outside printable ASCII and newline")
    header = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{N}];"]
    lines = []
    for raw in text.split("\n"):
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        lines.append(line)
    if lines[:3] != header:
        raise ValueError("header must be exactly: " + " | ".join(header))
    gates = []
    for line in lines[3:]:
        m = STMT.match(line)
        if not m:
            raise ValueError("not a single allowed statement: " + line[:80])
        name = m.group(1)
        qs = [int(m.group(2))] + ([int(m.group(3))] if m.group(3) is not None else [])
        if (name in ("cx", "cz")) != (len(qs) == 2):
            raise ValueError("wrong arity: " + line)
        if len(set(qs)) != len(qs) or any(not 0 <= q < N for q in qs):
            raise ValueError("bad qubit index: " + line)
        gates.append((name, qs))
    return gates


def split_layers(gates, n):
    """Return the diagonal part between the two H layers on c, or raise."""
    hc = sorted(2 * n + k for k in range(n))
    for layer in (gates[:n], gates[-n:]):
        if any(g != "h" for g, _ in layer) or sorted(qs[0] for _, qs in layer) != hc:
            raise ValueError("circuit must start and end with h on every c wire")
    mid = gates[n:-n]
    allowed = {"cx": 2, "cz": 2, "t": 1, "tdg": 1, "s": 1, "sdg": 1, "z": 1}
    for g, qs in mid:
        if allowed.get(g) != len(qs):
            raise ValueError(f"gate {g} not allowed between the H layers")
    return mid


PHASE = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}


def circuit_identity(n, p, mid):
    """Exact Z_8 phase polynomial of the diagonal part; True iff it equals 4f and the wires return."""
    N = 3 * n
    form = [1 << i for i in range(N)]
    lin = [0] * N          # coefficient of x_i (mod 8)
    pair = {}              # coefficient of x_i x_j (mod 8)
    trip = {}              # parity of the cubic coefficient / 4
    tforms = []
    for g, qs in mid:
        if g == "cx":
            c, t = qs
            form[t] ^= form[c]
            continue
        if g == "cz":
            L, M = form[qs[0]], form[qs[1]]
            # 4 (L.x)(M.x) mod 8 = 4 * sum_{i in L, j in M} x_i x_j  (x_i^2 = x_i)
            for i in range(N):
                if not L >> i & 1:
                    continue
                for j in range(N):
                    if M >> j & 1:
                        if i == j:
                            lin[i] = (lin[i] + 4) % 8
                        else:
                            k = (min(i, j), max(i, j))
                            pair[k] = (pair.get(k, 0) + 4) % 8
            continue
        m = PHASE[g]
        L = form[qs[0]]
        if g in ("t", "tdg"):
            tforms.append(L)
        s = [i for i in range(N) if L >> i & 1]
        for x in range(len(s)):
            lin[s[x]] = (lin[s[x]] + m) % 8
            for y in range(x + 1, len(s)):
                k = (s[x], s[y])
                pair[k] = (pair.get(k, 0) - 2 * m) % 8
                if m & 1:
                    above = L & ~((1 << (s[y] + 1)) - 1)
                    if above:
                        trip[k] = trip.get(k, 0) ^ above
    wires_ok = all(form[i] == 1 << i for i in range(N))
    want = target_rows(n, p)
    cubic_ok = all(trip.get(k, 0) == want.get(k, 0) for k in set(trip) | set(want))
    ok = wires_ok and not any(lin) and not any(pair.values()) and cubic_ok
    return ok, tforms


# ------------------------------------------------------------------ check 4: bit-sliced simulation

def var_slices(N, exhaustive, width):
    if exhaustive:
        size = 1 << N
        sl = []
        for i in range(N):
            per = 1 << (i + 1)
            pat = ((1 << (1 << i)) - 1) << (1 << i)
            s, ln = pat, per
            while ln < size:
                s |= s << ln
                ln *= 2
            sl.append(s)
        return sl, (1 << size) - 1
    return [RNG.getrandbits(width) for _ in range(N)], (1 << width) - 1


def simulate(n, p, mid, exhaustive):
    N = 3 * n
    xs, full = var_slices(N, exhaustive, 2048)
    w = list(xs)
    p0 = p1 = p2 = 0

    def add(v, plane):
        nonlocal p0, p1, p2
        if plane == 0:
            c = p0 & v
            p0 ^= v
            v = c
        if plane <= 1:
            c = p1 & v
            p1 ^= v
            v = c
        p2 ^= v

    for g, qs in mid:
        if g == "cx":
            w[qs[1]] ^= w[qs[0]]
        elif g == "cz":
            add(w[qs[0]] & w[qs[1]], 2)
        else:
            m, v = PHASE[g], w[qs[0]]
            for b in range(3):
                if m >> b & 1:
                    add(v, b)
    f = 0
    for i in range(n):
        for j in range(n):
            ab = xs[i] & xs[n + j]
            v = pmod(1 << (i + j), p)
            for k in range(n):
                if v >> k & 1:
                    f ^= ab & xs[2 * n + k]
    return w == xs and p0 == 0 and p1 == 0 and p2 == f & full


def metadata_ok(n, F, r, gates):
    """gate_counts must equal the counts of the QASM statements, and clifford_correction must equal the
    correction determined by the factors (Lemma 1): S_i^(-N_i/2) and CZ where N_ij = 2 mod 4."""
    if gates is None:
        return False
    N = 3 * n
    names = [g for g, _ in gates]
    want = {"qubits": N, "T": names.count("t") + names.count("tdg"), "CNOT": names.count("cx"),
            "S": names.count("s"), "Sdg": names.count("sdg"), "Z": names.count("z"), "CZ": names.count("cz"),
            "H": names.count("h")}
    if r.get("gate_counts") != want:
        return False
    lin = [0] * N
    pair = {}
    for u in F:
        s = [i for i in range(N) if u >> i & 1]
        for x in range(len(s)):
            lin[s[x]] += 1
            for y in range(x + 1, len(s)):
                k = (s[x], s[y])
                pair[k] = (pair.get(k, 0) + 1) & 3
    sp = {str(i): ((-lin[i]) % 8) // 2 for i in range(N) if (-lin[i]) % 8}
    cz = sorted([list(k) for k, m in pair.items() if m == 2])
    cc = r.get("clifford_correction", {})
    return cc.get("S_power") == sp and sorted(cc.get("CZ", [])) == cz


# ------------------------------------------------------------------ curve data

def curve_ok(d):
    K = int(d["K_modulus"], 2)
    Vm = int(d["V_modulus"], 2)
    kmul = lambda a, b: pmul(a, b, K)
    vmul = lambda a, b: pmul(a, b, Vm)

    def kpow(a, e):
        r = 1
        for _ in range(e):
            r = kmul(r, a)
        return r

    def vpow2k(a, k):
        for _ in range(k):
            a = vmul(a, a)
        return a

    pts = [(x, y) for x in range(16) for y in range(16) if kmul(y, y) ^ y == kpow(x, 5)]
    if len(pts) != 32 or not is_irreducible(K) or not is_irreducible(Vm) or Vm.bit_length() != 65:
        return False, "curve point count or moduli"
    theta = int(d["theta"], 16)
    # theta is a root of K's modulus in V
    acc, tp = 0, 1
    for i in range(K.bit_length()):
        if K >> i & 1:
            acc ^= tp
        tp = vmul(tp, theta)
    if acc:
        return False, "theta is not a root of the K modulus"
    xq, yq = int(d["x_Q"], 16), int(d["y_Q"], 16)
    x5 = vmul(vmul(vmul(xq, xq), vmul(xq, xq)), xq)
    if vmul(yq, yq) ^ yq != x5:
        return False, "Q is not on the curve"
    if vpow2k(xq, 32) == xq and vpow2k(yq, 32) == yq:
        return False, "Q does not have degree 16"
    kv_basis = [1]
    for _ in range(3):
        kv_basis.append(vmul(kv_basis[-1], theta))

    def kV(c):
        v = 0
        for i in range(4):
            if c >> i & 1:
                v ^= kv_basis[i]
        return v

    beta = [int(b, 16) for b in d["beta"]]
    cols = [vmul(kv_basis[r], beta[s]) for s in range(16) for r in range(4)]  # index s*4 + r
    # Gaussian elimination: express any a in V in the basis cols
    piv = []  # (vector, combination)
    for idx, v in enumerate(cols):
        comb = 1 << idx
        for pv, pc in piv:
            if v ^ pv < v:
                v ^= pv
                comb ^= pc
        if not v:
            return False, "beta does not give a basis of V"
        piv.append((v, comb))
        piv.sort(reverse=True)

    def coords(a):
        comb = 0
        for pv, pc in piv:
            if a ^ pv < a:
                a ^= pv
                comb ^= pc
        assert a == 0
        return [(comb >> (4 * s)) & 15 for s in range(16)]

    lam = d["lambda"]
    om = [int(o, 16) for o in d["omega"]]
    if len(lam) != 33 or len(om) != 33:
        return False, "expected 33 functionals"
    for _ in range(20):
        a, b = RNG.getrandbits(64), RNG.getrandbits(64)
        ca, cb = coords(a), coords(b)
        tot = 0
        for k in range(33):
            la = lb = 0
            for s in range(16):
                la ^= kmul(ca[s], lam[k][s])
                lb ^= kmul(cb[s], lam[k][s])
            tot ^= vmul(kV(kmul(la, lb)), om[k])
        if tot != vmul(a, b):
            return False, "bilinear identity fails"
    return True, "33 rational points, Q of degree 16 on the curve, identity a*b = sum_k l_k(a) l_k(b) w_k on 20 random pairs"


# ------------------------------------------------------------------ main

def main():
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parent)
    fails = 0
    man = root / "MANIFEST.sha256"
    actual = {str(q.relative_to(root)) for q in root.rglob("*") if q.is_file() and q.name != "MANIFEST.sha256"
              and "__pycache__" not in q.parts}
    problems = []
    listed = {}
    if not man.exists():
        problems.append("missing")
    else:
        for ln in man.read_text().split("\n"):
            if not ln:
                continue
            m = HEX64.match(ln)
            if not m or ".." in m.group(2).split("/") or m.group(2).startswith("/") or m.group(2) in listed:
                problems.append("malformed line")
                break
            listed[m.group(2)] = m.group(1)
        if set(listed) != actual:
            problems.append(f"file set differs (unlisted {sorted(actual - set(listed))[:3]}, absent {sorted(set(listed) - actual)[:3]})")
        for f, h in listed.items():
            if f in actual and hashlib.sha256((root / f).read_bytes()).hexdigest() != h:
                problems.append("hash mismatch " + f)
    fails += bool(problems)
    print(f"manifest: {'ok' if not problems else '; '.join(problems[:4])}")
    seen = {}
    print("records (n, T-count; reported value, best checked published value): irreducible | tensor criterion | exact circuit identity | simulation")
    recs_table = json.loads((root / "data" / "records.json").read_text())["records"]
    def nkey(q):
        try:
            v = json.loads(q.read_text()).get("n")
            return v if isinstance(v, int) and not isinstance(v, bool) else 10 ** 6
        except Exception:
            return 10 ** 6
    for f in sorted((root / "records").glob("*.json"), key=nkey):
        r = json.loads(f.read_text())
        n = r.get("n")
        if not (isinstance(n, int) and not isinstance(n, bool) and 2 <= n <= 64 and well_formed(r, n)):
            print(f"  {f.name}: malformed record (types, modulus or factor format, modulus_poly)")
            fails += 1
            continue
        if not isinstance(r.get("comparison"), dict) or not isinstance(r.get("clifford_correction"), dict) \
                or not isinstance(r.get("published_record"), dict) or not isinstance(r.get("gate_counts"), dict):
            print(f"  {f.name}: malformed record (field types)")
            fails += 1
            continue
        cmp_ = r.get("comparison", {})
        want = EXPECTED_COMPARISON.get(n)
        if want is None or (cmp_.get("reported"), cmp_.get("reported_invalid"), cmp_.get("checked")) != want \
                or r.get("published_record") != recs_table.get(str(n)) \
                or (r.get("published_record") or {}).get("tcount") != EXPECTED_BARS.get(n):
            print(f"  GF(2^{n}): comparison / published_record fields differ from the paper's table")
            fails += 1
            continue
        p = int(r["modulus"], 2)
        F = [bits(s) for s in r["factors"]]
        irr = is_irreducible(p) and p.bit_length() - 1 == n
        ten = len(set(F)) == len(F) and all(F) and tensor_ok(n, p, F)
        ex = 3 * n <= 24
        try:
            gates = parse_qasm((root / "circuits" / (f.stem + ".qasm")).read_text(), n)
            mid = split_layers(gates, n)
            ident, tforms = circuit_identity(n, p, mid)
            same = sorted(tforms) == sorted(F) and len(tforms) == r["tcount"]
            sim = simulate(n, p, mid, ex)
        except (ValueError, OSError) as e:
            print(f"  GF(2^{n}) circuit rejected: {e}")
            ident = same = sim = False
        meta = metadata_ok(n, F, r, gates if ident is not False or same else None)
        ok = irr and ten and ident and same and sim and len(F) == r["tcount"] and meta and n not in seen
        seen[n] = r["tcount"]
        fails += not ok
        cmp = r.get("comparison", {})
        rv = f"{cmp.get('reported', '-')}{'*' if cmp.get('reported_invalid') else ''}"
        print(f"  GF(2^{n:<2})  T={r['tcount']:<5} ({rv:>6}, {str(cmp.get('checked') or '-'):>4})  "
              f"{'ok' if irr else 'FAIL'} | {'ok' if ten else 'FAIL'} | {'ok' if ident and same else 'FAIL'} | "
              f"{'ok' if sim else 'FAIL'} ({'all 2^' + str(3 * n) if ex else '2048 random'} inputs)")
    if seen != EXPECTED_RECORDS:
        fails += 1
        print(f"  RECORD SET MISMATCH: expected {EXPECTED_RECORDS}, found {seen}")
    print("blocks (algebra, T-count): factor count, distinct, modulus, tensor criterion")
    bseen = {}
    for f in sorted((root / "blocks").glob("*.json")):
        r = json.loads(f.read_text())
        d = r.get("dimension")
        if not (isinstance(d, int) and not isinstance(d, bool) and 1 <= d <= 8 and well_formed(r, d)):
            print(f"  {f.name}: malformed block")
            fails += 1
            continue
        F = [bits(s) for s in r["factors"]]
        want_t, pat = EXPECTED_BLOCKS.get(r["algebra"], (None, None))
        p = int(r["modulus"], 2)
        if pat is None:
            mod_ok = False
        elif pat.startswith("irreducible:"):
            mod_ok = is_irreducible(p) and p.bit_length() - 1 == int(pat.split(":")[1])
        else:
            mod_ok = r["modulus"] == pat
        ok = (len(F) == r["tcount"] == want_t and len(set(F)) == len(F) and all(F) and mod_ok
              and p.bit_length() - 1 == r["dimension"] and tensor_ok(r["dimension"], p, F))
        bseen[r["algebra"]] = r["tcount"]
        fails += not ok
        print(f"  {r['algebra']:<40} T={r['tcount']:<4} {'ok' if ok else 'FAIL'}")
    if bseen != {k: v[0] for k, v in EXPECTED_BLOCKS.items()}:
        fails += 1
        print("  BLOCK SET MISMATCH")
    cp = root / "construction" / "curve64.json"
    ok, msg = curve_ok(json.loads(cp.read_text())) if cp.exists() else (False, "missing")
    fails += not ok
    print(f"genus-2 construction data: {'ok' if ok else 'FAIL'} ({msg})")
    print("ALL CHECKS PASSED" if not fails else f"{fails} CHECK(S) FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
