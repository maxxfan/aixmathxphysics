"""auditor: local-symmetry facts for the H-type reduction of AME(8,4).
Sp(4,2) = automorphisms of the 4-dim symplectic F_2-space of one ququart party (Clifford group mod Paulis).
Check: |Sp(4,2)| = 720; elements of order 3: classes, which are fixed-point-free; centralizer order of an fpf order-3 element;
for fpf g: F_2^4 becomes F_4^2 (omega acts as g); the symplectic form is Tr of a Hermitian form (checked: for every nonzero x,
B(x, g x) pattern), centralizer = GU(2,2) of order 18.
usage: python3 work/auditor/sp42_classes.py
"""
import itertools

J = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]  # form B(x,y) = x^T J y  (x=(a1,a2|b1,b2))


def mm(A, B):
    return tuple(tuple(sum(A[i][k] & B[k][j] for k in range(4)) & 1 for j in range(4)) for i in range(4))


def tr(A):
    return tuple(tuple(A[j][i] for j in range(4)) for i in range(4))


ID = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
Jt = tuple(tuple(r) for r in J)
G = []
for bits in range(1 << 16):
    A = tuple(tuple((bits >> (4 * i + j)) & 1 for j in range(4)) for i in range(4))
    if mm(mm(tr(A), Jt), A) == Jt:
        G.append(A)
print("|Sp(4,2)| =", len(G))
idx = {A: i for i, A in enumerate(G)}


def inv(A):
    for B in G:
        if mm(A, B) == ID:
            return B


def order(A):
    P, k = A, 1
    while P != ID:
        P = mm(P, A)
        k += 1
    return k


def vec_fixed(A):
    n = 0
    for v in itertools.product((0, 1), repeat=4):
        w = tuple(sum(A[i][j] * v[j] for j in range(4)) & 1 for i in range(4))
        if w == v:
            n += 1
    return n  # includes 0


o3 = [A for A in G if A != ID and order(A) == 3]
print("elements of order 3:", len(o3))
fpf = [A for A in o3 if vec_fixed(A) == 1]
nonfpf = [A for A in o3 if vec_fixed(A) != 1]
print("  fixed-point-free:", len(fpf), " with fixed 2-space:", len(nonfpf))
# conjugacy classes among fpf
inv_cache = {}
def conj_class(A):
    S = set()
    for g in G:
        gi = inv_cache.get(g)
        if gi is None:
            gi = inv(g)
            inv_cache[g] = gi
        S.add(mm(mm(g, A), gi))
    return S
rem = set(fpf)
classes = []
while rem:
    A = next(iter(rem))
    C = conj_class(A)
    classes.append(len(C))
    rem -= C
print("  conjugacy class sizes of fpf order-3 elements:", classes)
A = fpf[0]
cent = [g for g in G if mm(g, A) == mm(A, g)]
print("  centralizer order of an fpf order-3 element:", len(cent), "(GU(2,2) has order 18)")
# F_4-structure: g^2+g+1=0
S = tuple(tuple((mm(A, A)[i][j] + A[i][j] + ID[i][j]) & 1 for j in range(4)) for i in range(4))
print("  g^2+g+1 == 0:", all(S[i][j] == 0 for i in range(4) for j in range(4)))
