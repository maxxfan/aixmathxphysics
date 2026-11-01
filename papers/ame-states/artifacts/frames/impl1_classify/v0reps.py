"""v0reps.py m : gauge classes of the first frame member V_0 for n = 2m (m = 2, 3, 4).
V_0 = {(w, psi_1 w, ..., psi_{m-1} w)} in W_0 + ... + W_{m-1}, every psi_l in GL(4,2).  Left gauge Sp(4,2)^m x S_m acts;
a class is an orbit of m-tuples (Q_0..Q_{m-1}) of nondegenerate alternating forms on V_0 (Q_l = pullback of J by pi_l),
with sum nondegenerate, under GL(V_0) (congruence) x S_m (permuting the tuple).  Output: one line per class
 'orbit_size psi_1 ... psi_{m-1}' (psi packed 16-bit matrices, psi^T J psi = g^T Q_l g with Q_0 normalised to J by g)."""
import itertools
import sys
from f2m import *


def main(m, out):
    idx = {q: i for i, q in enumerate(NA)}
    # generators of GL(4,2): the 12 elementary transvections I + E_{ij}
    gens = []
    for i in range(4):
        for j in range(4):
            if i != j:
                gens.append(I4 | (1 << (4 * i + j)))
    gact = [[idx[congr(g, q)] for q in NA] for g in gens]

    def valid(t):
        s = 0
        for i in t:
            s ^= NA[i]
        return rank(s) == 4

    seen = set()
    classes = []
    for t0 in itertools.product(range(28), repeat=m):
        if t0 in seen or not valid(t0):
            continue
        orb = {t0}
        stack = [t0]
        while stack:
            t = stack.pop()
            nxt = []
            for ga in gact:
                nxt.append(tuple(ga[i] for i in t))
            for k in range(m - 1):
                u = list(t)
                u[k], u[k + 1] = u[k + 1], u[k]
                nxt.append(tuple(u))
            u = list(t)
            u = u[1:] + u[:1]
            nxt.append(tuple(u))
            for u in nxt:
                if u not in orb:
                    orb.add(u)
                    stack.append(u)
        seen |= orb
        classes.append((len(orb), t0))
    total = sum(c[0] for c in classes)
    ntup = sum(1 for t in itertools.product(range(28), repeat=m) if valid(t))
    assert total == ntup, (total, ntup)
    psi_of = {}
    for a in GL:
        psi_of.setdefault(form(a), a)
    lines = []
    for size, t in sorted(classes, key=lambda c: -c[0]):
        q = [NA[i] for i in t]
        g = next(g for g in GL if congr(g, q[0]) == J)
        qs = [congr(g, x) for x in q]
        assert qs[0] == J
        psis = [psi_of[x] for x in qs[1:]]
        stab = 20160 * __import__("math").factorial(m) // size
        lines.append(f"{size} " + " ".join(str(p) for p in psis))
        print(f"class size {size} |stab|={stab} psis={psis}")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"m={m}: {len(classes)} classes, {ntup} valid tuples, written to {out}")


if __name__ == "__main__":
    m = int(sys.argv[1])
    main(m, sys.argv[2] if len(sys.argv) > 2 else f"v0_m{m}.txt")
