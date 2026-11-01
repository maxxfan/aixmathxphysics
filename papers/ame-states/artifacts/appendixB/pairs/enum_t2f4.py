"""enum_t2f4.py [THR] [LIMIT]: AME(8,4) F_2-linear stabilizer invariant under Omega = swaps (0,4)(1,5) + rank-2 local involutions
(15-class, <x,lx>=0) on fixed parties 2,3,6,7.  Derivation: S^Omega = N S = K (dim 8), K -> L = sum_j L_j (8-dim, symplectic b = sum_j b_j) bijective,
w_1 = A f, w_2 = B f; isotropy forces b = A*J + B*J, ker A = Z_A, ker B = Z_B, L = Z_A (+) Z_B b-orthogonal, both 4-dim symplectic with every
nonzero element of block weight >= 3 (else K word of weight <= 4); A|Z_B, B|Z_A isometries (unique up to local Clifford gauge).
Enumerates all such Z_A (no symmetry reduction), runs qdfs on each K."""
import sys, time, itertools
import os  # [release]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # [release] was the absolute work/pairs path
from invtw import *
import ame

THR = int(sys.argv[1]) if len(sys.argv) > 1 else 4
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
n = 8
G = sp42()
cls = inv_classes(G)
reps = [c[0] for rk, c in cls if rk == 2]
l = list(reps[1])
fixed = [2, 3, 6, 7]
om = Omega(n, [(0, 4), (1, 5)], {j: l for j in fixed})
st = Setup(om)
assert len(st.Vom) == 16, len(st.Vom)


def Lbasis(lm):
    la = LinAlg()
    out = []
    for v in range(16):
        w = apply_loc(lm, v) ^ v
        if la.insert(w, 0):
            out.append(w)
    return out


Lb = Lbasis(l)
lb = []
for j in fixed:
    for w in Lb:
        lb.append(w << (4 * j))


def vec(f):
    v = 0
    for i in range(8):
        if f >> i & 1:
            v ^= lb[i]
    return v


def bform(f, g):
    return st.beta(vec(f), vec(g))


B = [[bform(f, g) for g in range(256)] for f in range(256)]
assert all(B[f][f] == 0 for f in range(256)), "b not alternating"
assert all(B[f][g] == B[g][f] for f in range(256) for g in range(256))


def bw(f):
    return sum(1 for j in range(4) if (f >> (2 * j)) & 3)


MINBW = int(sys.argv[3]) if len(sys.argv) > 3 else 3
allowed = {f for f in range(1, 256) if bw(f) >= MINBW}
print("allowed vectors", len(allowed), flush=True)
planes = []
for a in allowed:
    for a2 in allowed:
        if a2 > a and B[a][a2] == 1 and (a ^ a2) in allowed:
            planes.append((a, a2))
print("hyperbolic planes", len(planes), flush=True)
Zs = set()
for (a, a2) in planes:
    for (c, c2) in planes:
        if c < a:
            continue
        if B[a][c] or B[a][c2] or B[a2][c] or B[a2][c2]:
            continue
        # span of a,a2,c,c2
        sp_ = [0]
        for g in (a, a2, c, c2):
            sp_ += [g ^ s for s in sp_]
        if len(set(sp_)) != 16:
            continue
        if all(x in allowed for x in sp_ if x):
            Zs.add(frozenset(sp_))
print("candidate Z_A (symplectic, block weight>=3):", len(Zs), flush=True)


def orth(Z):
    return frozenset(f for f in range(256) if all(B[f][z] == 0 for z in Z))


def symp_basis(Z):
    Zl = sorted(Z)
    a = next(z for z in Zl if z)
    a2 = next(z for z in Zl if B[a][z])
    rest = [z for z in Zl if B[a][z] == 0 and B[a2][z] == 0]
    c = next(z for z in rest if z)
    c2 = next(z for z in rest if B[c][z])
    return [a, a2, c, c2]


NIB = [1, 4, 2, 8]  # e0,e2,e1,e3 : <e0,e2>=1, <e1,e3>=1


def vec_pair(p, nib):
    return (nib << (4 * p)) | (nib << (4 * (p + 4)))


t0 = time.time()
nK = nrun = nsol_tot = 0
pairsZ = []
for ZA in Zs:
    ZB = orth(ZA)
    if len(ZB) != 16 or not all(x in allowed for x in ZB if x):
        continue
    pairsZ.append((ZA, ZB))
print("ordered (Z_A,Z_B) pairs", len(pairsZ), "(unordered", len(pairsZ) // 2, ")", flush=True)
for (ZA, ZB) in pairsZ:
    bA, bB = symp_basis(ZA), symp_basis(ZB)
    Kb = []
    for u, nib in zip(bA, NIB):
        Kb.append(vec(u) | vec_pair(1, nib))   # Z_A -> w_2 (B isometry), A kills Z_A
    for u, nib in zip(bB, NIB):
        Kb.append(vec(u) | vec_pair(0, nib))   # Z_B -> w_1
    span = [0]
    for k in Kb:
        span += [k ^ s for s in span]
    if len(set(span)) != 256 or min(wt(v) for v in span[1:]) < 5:
        print("K invalid weights", min(wt(v) for v in span[1:]), flush=True)
        continue
    nK += 1
    if nrun >= LIMIT:
        continue
    nrun += 1
    r = run_qdfs(st, Kb, THR, 3, 3 if THR < 4 else 10**18)
    if r is None:
        print("INVALID K (isotropy)", flush=True)
        continue
    ns, sols, info = r
    nsol_tot += ns
    if nrun <= 3 or nrun % 200 == 0:
        print("  K#", nrun, info, "(%.0fs)" % (time.time() - t0), flush=True)
    if ns:
        print("SOLUTIONS", ns, "K=", " ".join("%x" % v for v in Kb), flush=True)
        for gens in sols:
            rows = to_rows(gens, n)
            la = LinAlg()
            rk = sum(1 for g in gens if la.insert(g, 0))
            iso = all(sform(a, b_, n) == 0 for a in gens for b_ in gens)
            inv = all(la.solve(om.apply(g)) is not None for g in gens)
            sp2 = [0]
            for g in gens:
                sp2 += [g ^ s_ for s_ in sp2]
            mw = min(wt(v) for v in sp2[1:])
            b, tot = ame.bad_cuts(rows, 2, 2, n)
            print("   SOL iso", iso, "rank", rk, "inv", inv, "minw", mw, "bad_cuts", b, "of", tot, "gens", " ".join("%x" % v for v in gens), flush=True)
print("TOTAL pairs", len(pairsZ), "K valid", nK, "ran", nrun, "nsol_total", nsol_tot, "(%.0fs)" % (time.time() - t0), flush=True)
