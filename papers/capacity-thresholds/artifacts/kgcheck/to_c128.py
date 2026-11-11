"""Write a p08 rank-2 n=45 candidate's rows as a raw little-endian complex128 payload (2 x 46, row-major)."""
import json
import struct
import sys

cand = json.load(open(sys.argv[1]))
rows = cand["rows"]
assert len(rows) == 2 and all(len(r) == 46 for r in rows)
vals = []
for row in rows:
    for z in row:
        assert isinstance(z, list) and len(z) == 2 and all(type(x) is float for x in z), z
        vals += [z[0], z[1]]
data = struct.pack("<184d", *vals)
open(sys.argv[2], "wb").write(data)
print(len(data), "bytes")
