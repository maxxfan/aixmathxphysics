"""Export for the C Lagrangian search: for each full completion T (case nonfpf), the set Vs of V^Omega elements s with pw(s+t) >= 5 for all t in T,
coded in 12 bits: bits 0-3 a nibble, 4-7 b nibble, 8-9 coefficients of F6 basis (p6,q6), 10-11 of F7 basis. omega = pair-swap formula.
Writes work/auditor/z3twist_Vs_CASE.bin (int32 nT, then per T: uint16 N, N x uint16) and z3twist_T_CASE.txt (T bases).  Usage: z3twist4.py CASE"""
import sys
from z3twist2 import *
import struct

F6 = [c for c in range(16) if h6[c] == c]
F7 = [c for c in range(16) if h7[c] == c]
assert len(F6) == 4 and len(F7) == 4
p6, q6 = F6[1], F6[2]
p7, q7 = F7[1], F7[2]
assert SY[p6][q6] == 1 and SY[p7][q7] == 1
def vec_of_code(code):
    a = code & 15; b = (code >> 4) & 15
    c = (p6 if (code >> 8) & 1 else 0) ^ (q6 if (code >> 9) & 1 else 0)
    d = (p7 if (code >> 10) & 1 else 0) ^ (q7 if (code >> 11) & 1 else 0)
    ar = lambda x: x | (x << 4) | (x << 8)
    return ar(a) | (ar(b) << 12) | (c << 24) | (d << 28)
codes = np.arange(4096, dtype=np.uint32)
VOm = np.array([vec_of_code(int(c)) for c in codes], dtype=np.uint32)
# check omega on codes == omega on vectors for random pairs
def omc(u, v):
    sv = ((v & 0x555) << 1) | ((v >> 1) & 0x555)
    return bin(u & sv).count("1") & 1
random.seed(7)
for _ in range(500):
    u = random.randrange(4096); v = random.randrange(4096)
    assert omc(u, v) == int(omg(np.uint32(vec_of_code(u)), np.uint32(vec_of_code(v)))), (u, v)
print("code-omega consistent with vector-omega")
# also each code vector is Omega-fixed
chk = OmN(VOm[:200])
assert (chk == VOm[:200]).all()

path = os.path.join(os.environ.get("APPB_WORK", "."), "")  # [release] output to the working dir
Ts = [(key, kk, Tb) for key, st in results.items() for kk, Tb in st.items() if len(Tb) == full_dim]
out = open(path + "z3twist_Vs_%s.bin" % case, "wb")
out.write(struct.pack("<i", len(Ts)))
tx = open(path + "z3twist_T_%s.txt" % case, "w")
for i, (key, kk, Tb) in enumerate(Ts):
    Tw = words_of(list(Tb.values()))
    ok = np.ones(4096, dtype=bool)
    for t in Tw:
        ok &= pw(VOm ^ t) >= 5
    cs = codes[ok].astype(np.uint16)
    out.write(struct.pack("<H", len(cs)))
    out.write(cs.tobytes())
    tx.write("%d %s %s\n" % (i, key[:2], " ".join(str(v) for v in Tb.values())))
    if i % 100 == 0:
        print("exported", i, len(cs), flush=True)
out.close(); tx.close()
print("done", len(Ts))
