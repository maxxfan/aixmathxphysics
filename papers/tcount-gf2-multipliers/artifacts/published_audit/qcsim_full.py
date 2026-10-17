#!/usr/bin/env python3
"""Exact sparse state-vector simulation of a .qc circuit (H, X/CNOT/Toffoli 'tof', Z/CZ/CCZ 'Z',
T/T* 'T'/'T*', S/S* 'S'/'S*', 'P' = S) on basis inputs, to test what a GF(2^n) multiplier circuit computes.

usage: qcsim_full.py FILE.qc n MODULUS_BITS     (registers a0..,b0..,c0.. and optional ancillas)
Reports: out-of-place correctness (c = 0, ancillas = 0), in-place correctness (all c, ancillas = 0),
whether ancillas return to 0, and the T-count of the file.
"""
import cmath
import sys

W = cmath.exp(1j * cmath.pi / 4)


def parse(path):
    vs, ins, gates = [], [], []
    for line in open(path):
        t = line.split()
        if not t or t[0] in ("BEGIN", "END"):
            continue
        if t[0] == ".v":
            vs = t[1:]
            continue
        if t[0] == ".i":
            ins = t[1:]
            continue
        if t[0].startswith("."):
            continue
        gates.append((t[0], t[1:]))
    return vs, ins, gates


def simulate(vs, gates, x0):
    idx = {v: i for i, v in enumerate(vs)}
    st = {x0: 1 + 0j}
    for g, args in gates:
        qs = [idx[a] for a in args]
        if g == "H":
            (q,) = qs
            new = {}
            for x, a in st.items():
                b = x >> q & 1
                x0_, x1_ = x & ~(1 << q), x | (1 << q)
                new[x0_] = new.get(x0_, 0) + a / 2 ** .5
                new[x1_] = new.get(x1_, 0) + (-a if b else a) / 2 ** .5
            st = {x: a for x, a in new.items() if abs(a) > 1e-9}
        elif g in ("tof", "X", "cnot"):
            *ctl, t = qs
            st2 = {}
            for x, a in st.items():
                if all(x >> c & 1 for c in ctl):
                    x ^= 1 << t
                st2[x] = st2.get(x, 0) + a
            st = st2
        else:
            ph = {"Z": -1, "T": W, "T*": W.conjugate(), "S": 1j, "P": 1j, "S*": -1j}[g]
            st = {x: (a * ph if all(x >> q & 1 for q in qs) else a) for x, a in st.items()}
    return st


def fmul(a, b, P, n):
    r = 0
    for i in range(n):
        if a >> i & 1:
            r ^= b << i
    for d in range(2 * n - 2, n - 1, -1):
        if r >> d & 1:
            r ^= P << (d - n)
    return r


def main():
    path, n, P = sys.argv[1], int(sys.argv[2]), int(sys.argv[3], 2)
    vs, ins, gates = parse(path)
    idx = {v: i for i, v in enumerate(vs)}
    A = [idx[f"a{i}"] for i in range(n)]
    B = [idx[f"b{i}"] for i in range(n)]
    C = [idx[f"c{i}"] for i in range(n)]
    anc = [i for i, v in enumerate(vs) if i not in A + B + C]
    tcount = sum(1 for g, _ in gates if g in ("T", "T*"))
    ccz = sum(1 for g, a in gates if g == "Z" and len(a) == 3)
    print(f"{path.split('/')[-1]}: {len(vs)} qubits ({len(anc)} ancillas), inputs .i = {' '.join(ins)}; "
          f"T gates {tcount}, CCZ {ccz}")

    def enc(a, b, c):
        x = 0
        for i in range(n):
            x |= (a >> i & 1) << A[i] | (b >> i & 1) << B[i] | (c >> i & 1) << C[i]
        return x

    out_ok = in_ok = anc_ok = tot_in = 0
    phase_issue = 0
    for a in range(1 << n):
        for b in range(1 << n):
            for c in range(1 << n):
                st = simulate(vs, gates, enc(a, b, c))
                big = [(x, amp) for x, amp in st.items() if abs(amp) > 1e-6]
                want = enc(a, b, c ^ fmul(a, b, P, n))
                good = len(big) == 1 and big[0][0] == want
                if len(big) == 1 and any(big[0][0] >> q & 1 for q in anc):
                    pass
                else:
                    anc_ok += 1
                if c == 0:
                    out_ok += good
                in_ok += good
                tot_in += 1
    print(f"  out-of-place (c=0, ancillas 0): correct on {out_ok}/{1 << 2 * n}")
    print(f"  in-place (every c, ancillas 0): correct on {in_ok}/{tot_in}")


if __name__ == "__main__":
    main()
