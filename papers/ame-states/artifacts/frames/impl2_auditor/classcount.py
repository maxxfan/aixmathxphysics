"""Auditor shift 22: re-count the numbers quoted in classify's drop-in subsection (paper_subsection_classification.tex).
Alternating 4x4 F_2 matrix = 6 bits (a01,a02,a03,a12,a13,a23); nondegenerate iff Pfaffian a01*a23+a02*a13+a03*a12 = 1.
Counts tuples (Q_0..Q_{h-1}) of nondegenerate forms with nondegenerate sum, plus orbit sizes under GL(4,2) x S_h
(orbit computation by explicit congruence Q -> g^T Q g, own code, full GL(4,2) by enumeration)."""
import itertools, sys

def pf(x):
    a01, a02, a03, a12, a13, a23 = [(x >> i) & 1 for i in range(6)]
    return (a01 & a23) ^ (a02 & a13) ^ (a03 & a12)

def to_mat(x):
    M = [[0] * 4 for _ in range(4)]
    idx = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for i, (r, c) in enumerate(idx):
        if (x >> i) & 1:
            M[r][c] = M[c][r] = 1
    return M

def from_mat(M):
    idx = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    return sum(M[r][c] << i for i, (r, c) in enumerate(idx))

def matmul(A, B):
    return [[sum(A[i][k] & B[k][j] for k in range(4)) & 1 for j in range(4)] for i in range(4)]

def tr(A):
    return [[A[j][i] for j in range(4)] for i in range(4)]

def det_f2(M):
    rows = [sum(M[i][j] << j for j in range(4)) for i in range(4)]
    r = 0
    for c in range(4):
        p = None
        for i in range(r, 4):
            if (rows[i] >> c) & 1:
                p = i
                break
        if p is None:
            return 0
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(4):
            if i != r and (rows[i] >> c) & 1:
                rows[i] ^= rows[r]
        r += 1
    return 1

GL = []
for bits in range(1 << 16):
    M = [[(bits >> (4 * i + j)) & 1 for j in range(4)] for i in range(4)]
    if det_f2(M):
        GL.append(M)
assert len(GL) == 20160, len(GL)

# congruence action on the 6-bit codes: table[g][x] = from_mat(g^T Q_x g); nondegenerate forms only (alternating stays alternating)
ND = [x for x in range(64) if pf(x)]
assert len(ND) == 28
act = []
for g in GL:
    gt = tr(g)
    act.append({x: from_mat(matmul(matmul(gt, to_mat(x)), g)) for x in ND})

def run(h):
    tuples = [t for t in itertools.product(ND, repeat=h) if pf(_xor(t))]
    S = set(tuples)
    print("h=%d: tuples with nondegenerate sum: %d (of %d)" % (h, len(tuples), 28 ** h))
    seen = set()
    sizes = []
    for t in tuples:
        if t in seen:
            continue
        orb = set()
        for a in act:
            u = tuple(a[x] for x in t)
            for perm in itertools.permutations(range(h)):
                orb.add(tuple(u[p] for p in perm))
        assert orb <= S
        seen |= orb
        sizes.append(len(orb))
    print("   orbits under GL(4,2) x S_%d: %d, sizes %s, total %d" % (h, len(sizes), sorted(sizes, reverse=True), sum(sizes)))

def _xor(t):
    r = 0
    for x in t:
        r ^= x
    return r

for h in (2, 3):
    run(h)
if len(sys.argv) > 1 and sys.argv[1] == "4":
    run(4)
