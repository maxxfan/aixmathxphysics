"""enum_t3f2.py CA CB [RUN] : AME(8,4) F_2-linear stabilizer invariant under Omega = 3 swaps (0,3)(1,4)(2,5) + rank-2 local
involutions l_6 (class CA), l_7 (class CB) on fixed parties (classes: 0 = 45-class, 1 = 15-class of rank-2 involutions in Sp(4,2)).
K = S^Omega = K_0 + span(k_1..k_4): K_0 = {(x,Mx,Nx;0,0)} (M,N coset reps of GL(4,2)/Sp(4,2), J+F_M+F_N=0), k_i = (t_i; phi_i) lifts.
Counts all K (isotropic, weight>=5), optionally runs qdfs (RUN=1) on each.  Output: log lines + solutions."""
import sys, time, itertools
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))  # [release] was the absolute work/pairs path
from invtw import *
import ame

CA, CB = int(sys.argv[1]), int(sys.argv[2])
RUN = int(sys.argv[3]) if len(sys.argv) > 3 else 0
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
THR = int(sys.argv[5]) if len(sys.argv) > 5 else 4
n = 8
G = sp42()
cls = inv_classes(G)
reps = [c[0] for rk, c in cls if rk == 2]
lA, lB = list(reps[CA]), list(reps[CB])
om = Omega(n, [(0, 3), (1, 4), (2, 5)], {6: lA, 7: lB})
st = Setup(om)
assert len(st.Vom) == 16, len(st.Vom)

# GL(4,2) coset reps by form
forms = {}
for M in itertools.product(range(1, 16), repeat=4):
    if len({apply_loc(list(M), v) for v in range(16)}) != 16:
        continue
    f = tuple(sp_nib(M[i], M[j]) for i in range(4) for j in range(i + 1, 4))
    if f not in forms:
        forms[f] = list(M)
Jf = tuple(sp_nib(1 << i, 1 << j) for i in range(4) for j in range(i + 1, 4))
print("forms", len(forms), "J", Jf, flush=True)


def vec_pair(p, nib):  # nibble at party p and p+3
    return (nib << (4 * p)) | (nib << (4 * (p + 3)))


def Lbasis(l):
    la = LinAlg()
    out = []
    for v in range(16):
        w = apply_loc(l, v) ^ v
        if la.insert(w, 0):
            out.append(w)
    return out


L6, L7 = Lbasis(lA), Lbasis(lB)
phis = [v << 24 for v in L6] + [v << 28 for v in L7]
print("L6", L6, "L7", L7, flush=True)
# check V^Omega = span of pair vecs + L6,L7
t0 = time.time()
nK0 = nK = nrun = nsol_tot = 0
sol_found = []
for fM, M in forms.items():
    for fN, N_ in forms.items():
        if any(Jf[i] ^ fM[i] ^ fN[i] for i in range(6)):
            continue
        K0 = [vec_pair(0, 1 << a) | vec_pair(1, apply_loc(M, 1 << a)) | vec_pair(2, apply_loc(N_, 1 << a)) for a in range(4)]
        K0span = [0]
        for k in K0:
            K0span += [k ^ s for s in K0span]
        nK0 += 1
        # candidates per phi: t in 12 bits -> vector, canonical mod K0, beta-orthogonal to K0, beta(k,k)=0
        cands = []
        for phi in phis:
            cl = []
            for tt in range(1 << 12):
                tv = vec_pair(0, tt & 15) | vec_pair(1, (tt >> 4) & 15) | vec_pair(2, (tt >> 8) & 15)
                if tv != min(tv ^ s for s in K0span):
                    continue
                k = tv | phi
                if st.beta(k, k) or any(st.beta(k, d) for d in K0):
                    continue
                if any(wt(k ^ s) < 5 for s in K0span):
                    continue
                cl.append(k)
            cands.append(cl)
        # DFS over the 4 lifts
        def rec(i, chosen, span):
            global nK, nrun, nsol_tot
            if i == 4:
                nK += 1
                if RUN and nrun < LIMIT:
                    nrun += 1
                    Kb = K0 + chosen
                    r = run_qdfs(st, Kb, THR, 3, 3 if THR < 4 else 10**18)
                    if r is None:
                        print("INVALID K?!", flush=True)
                        return
                    ns, sols, info = r
                    nsol_tot += ns
                    if nrun <= 3 or nrun % 1000 == 0:
                        print('  K#', nrun, info, flush=True)
                    if ns:
                        print("SOLUTIONS", ns, "K=", " ".join("%x" % v for v in Kb), flush=True)
                        for gens in sols:
                            rows = to_rows(gens, n)
                            la = LinAlg()
                            rk = sum(1 for g in gens if la.insert(g, 0))
                            iso = all(sform(a, b_, n) == 0 for a in gens for b_ in gens)
                            inv = all(la.solve(om.apply(g)) is not None for g in gens)
                            span = [0]
                            for g in gens:
                                span += [g ^ s_ for s_ in span]
                            mw = min(wt(v) for v in span[1:])
                            print("   check iso", iso, "rank", rk, "Omega-inv", inv, "minweight", mw, flush=True)
                            b, tot = ame.bad_cuts(rows, 2, 2, n)
                            print("   SOL bad_cuts", b, "of", tot, "gens", " ".join("%x" % v for v in gens), flush=True)
                            sol_found.append(gens)
                return
            for k in cands[i]:
                if any(st.beta(k, c) for c in chosen):
                    continue
                # weights of k + span
                newspan = [k ^ s for s in span]
                if any(wt(v) < 5 for v in newspan):
                    continue
                rec(i + 1, chosen + [k], span + newspan)
        before = nK
        rec(0, [], K0span)
        print("K0 forms", fM, fN, "cands", [len(c) for c in cands], "K count", nK - before, "(%.0fs)" % (time.time() - t0), flush=True)
print("TOTAL K0", nK0, "K", nK, "ran", nrun, "nsol_total", nsol_tot, "(%.0fs)" % (time.time() - t0), flush=True)
