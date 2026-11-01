"""invtw.py -- involution-with-twist sectors for F_2-linear stabilizer AME(n,4) (2 qubits/party, nibble per party).
Vectors: int, nibble p (bits 4p..4p+3) = (x0,x1,z0,z1) of party p.  Omega = swaps of party pairs + local involutions on fixed parties.
Library: sp42(), involution_classes(), Omega class (apply, N, kerN basis, preimage), beta, build K, call qdfs.
"""
import sys, os, subprocess, itertools, random
# [release] paths: script directory, ../lib (copy of qlib/ame.py), qdfs binary from $BUILD (run.sh compiles pairs/qdfs.c);
# fallback: compile qdfs.c into the current working directory.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
BIN = os.path.join(os.environ.get("BUILD", ""), "qdfs")
if not os.environ.get("BUILD") or not os.path.exists(BIN):
    BIN = os.path.abspath("qdfs")
    if not os.path.exists(BIN):
        subprocess.run([os.environ.get("CC", "cc"), "-O2", "-o", BIN, os.path.join(HERE, "qdfs.c")], check=True)


def sp_nib(a, b):
    return (bin((a & 3) & (b >> 2)).count("1") + bin((a >> 2) & (b & 3)).count("1")) & 1


def apply_loc(M, v):  # M = list of 4 column images
    r = 0
    for i in range(4):
        if v >> i & 1:
            r ^= M[i]
    return r


def mat_mul(A, B):  # (A*B)(v) = A(B(v))
    return [apply_loc(A, B[i]) for i in range(4)]


def sp42():
    out = []
    E = [[sp_nib(1 << i, 1 << j) for j in range(4)] for i in range(4)]
    for M in itertools.product(range(1, 16), repeat=4):
        if any(sp_nib(M[i], M[j]) != E[i][j] for i in range(4) for j in range(i + 1, 4)):
            continue
        if len({apply_loc(list(M), v) for v in range(16)}) == 16:
            out.append(list(M))
    return out


def inv_classes(G):
    ID = [1, 2, 4, 8]
    invs = [M for M in G if mat_mul(M, M) == ID and M != ID]
    inv = lambda M: next(H for H in G if mat_mul(H, M) == ID)
    seen, classes = set(), []
    for M in invs:
        if tuple(M) in seen:
            continue
        cl = {tuple(mat_mul(mat_mul(g, M), inv(g))) for g in G}
        seen |= cl
        rk = len({apply_loc(M, v) ^ v for v in range(16)}).bit_length() - 1
        classes.append((rk, sorted(cl)))
    return classes


class Omega:
    def __init__(self, n, swaps, fixed):  # swaps: list of (p,q); fixed: dict party -> local matrix
        self.n, self.swaps, self.fixed = n, swaps, fixed
        self.NB = 4 * n

    def apply(self, v):
        r = 0
        nb = [(v >> (4 * p)) & 15 for p in range(self.n)]
        out = list(nb)
        for (p, q) in self.swaps:
            out[p], out[q] = nb[q], nb[p]
        for p, M in self.fixed.items():
            out[p] = apply_loc(M, nb[p])
        for p in range(self.n):
            r |= out[p] << (4 * p)
        return r

    def N(self, v):
        return v ^ self.apply(v)


def wt(v, n=8):
    return sum(1 for p in range(n) if (v >> (4 * p)) & 15)


def sform(a, b, n):
    return sum(sp_nib((a >> (4 * p)) & 15, (b >> (4 * p)) & 15) for p in range(n)) & 1


class LinAlg:
    """F_2 linear basis with combination tracking: insert (key, payload) pairs."""
    def __init__(self):
        self.piv = {}

    def insert(self, key, payload):
        while key:
            h = key.bit_length() - 1
            if h in self.piv:
                k2, p2 = self.piv[h]
                key ^= k2
                payload ^= p2
            else:
                self.piv[h] = (key, payload)
                return True
        return False

    def solve(self, key):  # returns payload with XOR keys == key or None
        payload = 0
        while key:
            h = key.bit_length() - 1
            if h not in self.piv:
                return None
            k2, p2 = self.piv[h]
            key ^= k2
            payload ^= p2
        return payload


class Setup:
    def __init__(self, om):
        self.om = om
        n = om.n
        # preimage solver: N(e_i) -> e_i
        self.pre = LinAlg()
        self.ker = []
        for i in range(4 * n):
            e = 1 << i
            img = om.N(e)
            # reduce img with payload e; if img reduces to 0 then combination in kernel
            key, pay = img, e
            while key:
                h = key.bit_length() - 1
                if h in self.pre.piv:
                    k2, p2 = self.pre.piv[h]
                    key ^= k2
                    pay ^= p2
                else:
                    self.pre.piv[h] = (key, pay)
                    break
            if key == 0:
                self.ker.append(pay)
        self.Vom = self.ker  # basis of ker N (= V^Omega); assumed = im N (checked below by dimension)

    def preimage(self, k):
        r = self.pre.solve(k)
        return r

    def beta(self, k, kp):
        return sform(self.preimage(k), kp, self.om.n)


def run_qdfs(setup, Kb, THR, maxprint=3, maxsol=10**18):
    """Kb: basis (len n) of K inside V^Omega. Returns (nsol, solutions(list of rows), info) or None if K invalid."""
    n = setup.om.n
    kd = len(Kb)
    p = [setup.preimage(k) for k in Kb]
    assert all(setup.om.N(p[i]) == Kb[i] for i in range(kd))
    for i in range(kd):
        for j in range(i, kd):
            if sform(p[i], Kb[j], n):
                return None
    a = [sum(sform(p[i], p[j], n) << j for j in range(kd)) for i in range(kd)]
    # eps_m in V^Omega with <p_i, eps_m> = delta_im
    la = LinAlg()
    for b in setup.Vom:
        key = sum(sform(p[i], b, n) << i for i in range(kd))
        la.insert(key, b)
    eps = []
    for m in range(kd):
        e = la.solve(1 << m)
        if e is None:
            return None
        eps.append(e)
    inp = "%d %d %d\n" % (n, kd, THR) + " ".join("%x" % v for v in Kb) + "\n" + " ".join("%x" % v for v in p) + "\n" + \
          " ".join("%x" % v for v in eps) + "\n" + " ".join("%x" % v for v in a) + "\n"
    pr = subprocess.run([BIN, str(maxprint), str(maxsol)], input=inp, capture_output=True, text=True)
    lines = pr.stdout.strip().splitlines()
    nsol = int([l for l in lines if l.startswith("NSOL")][0].split()[1])
    sols = []
    for l in lines:
        if l.startswith("SOL"):
            rowsq = [int(x, 16) for x in l.split()[1:]]
            gens = []
            for i in range(kd):
                u = p[i]
                for m in range(kd):
                    if rowsq[i] >> m & 1:
                        u ^= eps[m]
                gens.append(u)
            gens += list(Kb)
            sols.append(gens)
    return nsol, sols, lines[-1]


def to_rows(gens, n):
    """32-bit nibble vectors -> rows [x(2n) | z(2n)] for ame (qubit 2p+b of party p)."""
    rows = []
    for v in gens:
        x, z = [], []
        for p in range(n):
            nb = (v >> (4 * p)) & 15
            x += [nb & 1, (nb >> 1) & 1]
            z += [(nb >> 2) & 1, (nb >> 3) & 1]
        rows.append(x + z)
    return rows
