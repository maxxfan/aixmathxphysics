"""For each released matrix: is its signature tensor exactly GF(2^n) multiplication?
Registers assumed a = 0..n-1, b = n..2n-1, c = 2n..3n-1 (checked: every odd entry must be such a triple)."""
import sys, re, pathlib, random
import numpy as np

def poly_from_name(name, n):
    m = re.match(r"gf2_(\d+)_(\d+)$", name)
    if not m:
        return None
    digits = m.group(2)
    # degrees are concatenated (e.g. 520 = 5 2 0, 84320, 1612310 = 16 12 3 1 0, 644310 = 64 4 3 1 0)
    def parse(s, top):
        if not s: return []
        for k in (2, 1):
            if len(s) >= k and int(s[:k]) <= top and (k == 1 or s[0] != "0"):
                rest = parse(s[k:], int(s[:k]) - 1)
                if rest is not None and (not rest or rest[-1] == 0 or True):
                    if int(s[:k]) < top + 1:
                        return [int(s[:k])] + rest
        return None
    if not digits.startswith(str(n)):
        return None
    d = parse(digits[len(str(n)):], n - 1)
    if d is None or d[-1] != 0:
        return None
    P = 1 << n
    for e in d: P |= 1 << e
    return P

def pmod(a, P):
    dp = P.bit_length() - 1
    while a and a.bit_length() - 1 >= dp:
        a ^= P << (a.bit_length() - 1 - dp)
    return a

def analyze(path):
    M = np.load(path, allow_pickle=False).astype(np.uint8)
    T, N = M.shape
    if N % 3: return f"{path.name}: N={N} not divisible by 3 (ancillas?)"
    n = N // 3
    F = [int("".join(map(str, row[::-1])), 2) for row in M]
    # tensor parities via pair rows
    single = 0; pair = {}; trip = {}
    for u in F:
        s = [i for i in range(N) if u >> i & 1]
        for x in range(len(s)):
            single ^= 1 << s[x]
            for y in range(x + 1, len(s)):
                k = (s[x], s[y]); pair[k] = pair.get(k, 0) ^ 1
                above = u & ~((1 << (s[y] + 1)) - 1)
                if above: trip[k] = trip.get(k, 0) ^ above
    bad_low = bool(single) or any(pair.values())
    odd = [(i, j, k) for (i, j), m in trip.items() for k in range(N) if m >> k & 1]
    structured = all(i < n <= j < 2 * n <= k for (i, j, k) in odd)
    if not structured:
        return f"{path.name}: T={T} n={n}: odd triples not of (a,b,c) form ({len(odd)} triples), low-order even={not bad_low}"
    # bilinear map columns: col[i][j] = c-vector of a_i b_j
    col = {}
    for (i, j, k) in odd:
        col[(i, j - n)] = col.get((i, j - n), 0) ^ (1 << (k - 2 * n))
    P = poly_from_name(path.stem, n)
    exact = P is not None and all(col.get((i, j), 0) == pmod(1 << (i + j), P) for i in range(n) for j in range(n))
    # depends only on i+j?  (polynomial product followed by a linear map on degrees)
    by_sum = {}
    hankel = True
    for i in range(n):
        for j in range(n):
            v = col.get((i, j), 0)
            if by_sum.setdefault(i + j, v) != v: hankel = False
    # field test: for every a != 0 the map b -> M(a,b) must be invertible (exhaustive for n <= 16)
    def Ma_rank(a):
        rows = []
        for j in range(n):
            v = 0
            for i in range(n):
                if a >> i & 1: v ^= col.get((i, j), 0)
            rows.append(v)
        r = 0; basis = []
        for v in rows:
            for b in basis:
                v = min(v, v ^ b)
            if v: basis.append(v); r += 1
        return r
    if n <= 16:
        sing = sum(1 for a in range(1, 1 << n) if Ma_rank(a) < n)
        field = f"singular a: {sing}/{(1 << n) - 1}"
    else:
        rng = random.Random(1); tr = 4000
        sing = sum(1 for _ in range(tr) if Ma_rank(rng.getrandbits(n) or 1) < n)
        field = f"singular a among {tr} random: {sing}"
    return (f"{path.parent.name[:24]:<24} {path.name:<24} T={T:<6} n={n:<3} P={'x^'+str(n)+'..' if P else '?':<6} "
            f"low-order-even={not bad_low} equals-mult-mod-P={exact if P is not None else 'n/a'} depends-on-i+j={hankel} {field}")

if __name__ == "__main__":
  for p in sys.argv[1:]:
    try:
        print(analyze(pathlib.Path(p)))
    except Exception as e:
        print(p, "ERROR", e)
