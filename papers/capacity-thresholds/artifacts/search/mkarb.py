"""Write a candidate file (verifier format) for a saved state at p = float zero crossing - OFFSET, rounded
down to 1e-8, for certification by results/rigorous/arb_ic.py beyond the board verifier's caps.
usage mkarb.py STATE.json ZERO_CROSSING [OFFSET]   (writes arbcands/{target}_n{n}_{p}.json)"""
import json
import pathlib
import sys
from decimal import Decimal, ROUND_FLOOR

d = json.loads(pathlib.Path(sys.argv[1]).read_text())
p0 = Decimal(sys.argv[2])
off = Decimal(sys.argv[3]) if len(sys.argv) > 3 else Decimal("1e-5")
p = (p0 - off).quantize(Decimal("1e-8"), rounding=ROUND_FLOOR)
fam, n = d.get("target", "depolarizing"), len(d["rows"][0]) - 1
out = pathlib.Path(__file__).resolve().parent / "arbcands" / f"{fam}_n{n}_{p}.json"
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps({"target": fam, "p": str(p), "rows": d["rows"],
                           "note": f"{pathlib.Path(sys.argv[1]).name}; float zero crossing {p0}; p = crossing - {off}"}))
print(out)
