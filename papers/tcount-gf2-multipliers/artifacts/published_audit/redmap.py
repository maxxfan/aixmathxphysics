import sys, numpy as np, pathlib
def cols(path):
    M = np.load(path, allow_pickle=False).astype(np.uint8)
    n = M.shape[1] // 3
    F = [int("".join(map(str, r[::-1])), 2) for r in M]
    N = 3 * n; trip = {}
    for u in F:
        s = [i for i in range(N) if u >> i & 1]
        for x in range(len(s)):
            for y in range(x + 1, len(s)):
                k = (s[x], s[y]); above = u & ~((1 << (s[y] + 1)) - 1)
                if above: trip[k] = trip.get(k, 0) ^ above
    by = {}
    for (i, j), m in trip.items():
        if i < n <= j < 2 * n:
            by[i + j - n] = (m >> (2 * n))
    return n, [by.get(s, 0) for s in range(2 * n - 1)]
for p in sys.argv[1:]:
    n, c = cols(p)
    print(pathlib.Path(p).name, "x^s ->", [format(v, f"0{n}b")[::-1] for v in c])
