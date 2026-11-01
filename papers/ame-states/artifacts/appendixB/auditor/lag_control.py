"""Positive control for lagsearch2: build a Lagrangian L of F_2^12 (omega = pair-swap form) with the staged structure
(L cap {a=0} dim 2 inside b,c,d != 0; L cap {b=0} dim 2), write Vs = L words (and 200 random extra codes) as a one-T .bin, run lagsearch2.
Expect: Lagrangians >= 1 (and the control Vs without L's words -> 0)."""
import struct, random, subprocess, os
base = os.path.dirname(os.path.abspath(__file__))
def om(u, v):
    sv = ((v & 0x555) << 1) | ((v >> 1) & 0x555)
    return bin(u & sv).count("1") & 1
random.seed(5)
e = [(1 << 4) | (1 << 8) | (1 << 10), (4 << 4) | (2 << 8) | (2 << 10), 1 | (1 << 8) | (1 << 10), 4 | (2 << 8) | (2 << 10)]
for i in range(4):
    for j in range(4):
        assert om(e[i], e[j]) == 0, (i, j)
def span(vs):
    s = {0}
    for v in vs:
        s |= {x ^ v for x in s}
    return s
S4 = span(e)
assert len(S4) == 16
P = [v for v in range(4096) if all(om(v, x) == 0 for x in e)]
assert len(P) == 256
while True:
    e5 = random.choice([v for v in P if v not in S4])
    S5 = span(e + [e5])
    cand = [v for v in P if om(v, e5) == 0 and v not in S5]
    if cand:
        e6 = random.choice(cand); break
L = span(e + [e5, e6])
assert len(L) == 64
assert all(om(x, y) == 0 for x in L for y in L)
def run(vs, name):
    path = os.path.join(os.environ.get("APPB_WORK", "."), "ctl_%s.bin" % name)  # [release] output to the working dir
    with open(path, "wb") as f:
        f.write(struct.pack("<i", 1))
        vs = sorted(set(vs))
        f.write(struct.pack("<H", len(vs)))
        f.write(struct.pack("<%dH" % len(vs), *vs))
    # [release] binary from build/ (compiled by run.sh) instead of the script dir
    out = subprocess.run([os.path.join(os.environ.get("BUILD", os.path.join(base, "..", "build")), "lagsearch2"), path, "0", "0"], capture_output=True, text=True).stdout
    print(name, out.strip())
nz = [x for x in L if x]
run(nz, "L_only")
extra = [random.randrange(1, 4096) for _ in range(200)]
run(nz + extra, "L_plus_random")
run([x for x in nz if x != nz[7]] + extra, "L_minus_one_word_plus_random")
