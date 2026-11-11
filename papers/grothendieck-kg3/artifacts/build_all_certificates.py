#!/usr/bin/env python3
"""Maintainer script: fill all_certificates/ from the project's candidates/ directory (not part of the bundle).

usage: python3 release/build_all_certificates.py [--no-sanitize]        (from the problem directory; then run
       python3 -I release/make_manifest.py)

For every accepted certificate in results/wins.json it writes all_certificates/<run file name>.json.xz (LZMA; the
41 files take about 25 MB instead of 223 MB) and one entry in all_certificates/INDEX.json. The certificate of the
main theorem is not duplicated: its entry points to certificate/KG3_upper_1.4494177506_rad721.json, which is left
untouched.

Sanitizing. Some certificates carry, in the free-text metadata field "script", absolute paths of the machine the
search ran on. Unless --no-sanitize is given, the directory prefix up to and including the problem directory is
removed from that field, by a byte-level replacement of that one JSON string. Nothing else in the file changes;
this is checked (the parsed objects agree outside "script", and the verifier's canonical hash, which ignores the
metadata fields, is unchanged). INDEX.json records for every file the SHA-256 of the JSON as shipped, the 16-digit
SHA-256 prefix of the file as it was submitted during the run (the one in the file name and on the board), the
canonical hash, and whether "script" was edited. The theorem's certificate is never edited by this script.
"""
import hashlib
import json
import lzma
import math
import pathlib
import re
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
P = HERE.parent
OUT = HERE / "all_certificates"
FINAL = "certificate/KG3_upper_1.4494177506_rad721.json"
sys.path.insert(0, str(P / "verifier"))
import verify  # noqa: E402  (only canonical_hash is used)

# direction sets, keyed by (m, eta^2 to 9 digits) as in paper/make_tables.py
SETS = {
    (406, "0.993522758"): "pol406", (406, "0.993612508"): "pol406_1/4", (541, "0.995337676"): "pol541",
    (541, "0.995920093"): "rad541", (541, "0.996231172"): "rad541'", (606, "0.995904004"): "geo606",
    (721, "0.996560718"): "geo721", (721, "0.996691977"): "rad721",
}
ABS = re.compile(r"/\S*?/p04-grothendieck-kg3/")


def down(q, d):
    s = 10 ** d
    v = math.floor(Fraction(q) * s)
    return f"{v // s}.{str(v % s).rjust(d, '0')}"


def sanitize(raw):
    """Return (bytes, edited): `raw` with absolute directory prefixes removed from the "script" string only."""
    cand = json.loads(raw)
    old = cand.get("script")
    if not isinstance(old, str) or not ABS.search(old):
        return raw, False
    new = ABS.sub("", old)
    for enc in (json.dumps(old), json.dumps(old, ensure_ascii=False)):
        if raw.count(enc.encode()) == 1:
            out = raw.replace(enc.encode(), json.dumps(new).encode())
            break
    else:
        raise SystemExit("script string not found exactly once in the raw file")
    c2 = json.loads(out)
    if c2.pop("script") != new or {k: v for k, v in cand.items() if k != "script"} != c2:
        raise SystemExit("sanitizing changed more than the script field")
    if verify.canonical_hash(cand) != verify.canonical_hash(json.loads(out)):
        raise SystemExit("sanitizing changed the canonical hash")
    return out, True


def main():
    do_sanitize = "--no-sanitize" not in sys.argv
    wins = json.loads((P / "results" / "wins.json").read_text())["wins"]
    best = {}
    for w in wins:
        name = SETS[(w["mA"], f"{float(w['eta2_A']):.9f}")]
        if name not in best or Fraction(w["kg3_upper"]) < Fraction(best[name]["kg3_upper"]):
            best[name] = w
    table2 = {w["post"] for w in best.values()}
    final_post = min(wins, key=lambda w: Fraction(w["kg3_upper"]))["post"]
    OUT.mkdir(exist_ok=True)
    for f in OUT.glob("*.json.xz"):
        f.unlink()
    entries, total = [], 0
    for w in wins:
        src = P / w["file"]
        raw = src.read_bytes()
        run16 = hashlib.sha256(raw).hexdigest()[:16]
        if run16 != w["sha256_16"] or run16 not in src.name:
            raise SystemExit(f"{src.name}: file differs from the one recorded in wins.json")
        cand = json.loads(raw)
        if w["post"] == final_post:
            rel = FINAL
            if (HERE / FINAL).read_bytes() != raw:
                raise SystemExit("the theorem's certificate differs from its candidates/ original")
            shipped, edited = raw, False
        else:
            shipped, edited = sanitize(raw) if do_sanitize else (raw, False)
            rel = f"all_certificates/{src.stem}.json.xz"
            blob = lzma.compress(shipped, preset=6)
            if lzma.decompress(blob) != shipped:
                raise SystemExit("LZMA round trip failed")
            (HERE / rel).write_bytes(blob)
            total += len(blob)
        entries.append({
            "post": w["post"], "time": w["time"], "agent": w["agent"],
            "set": SETS[(w["mA"], f"{float(w['eta2_A']):.9f}")], "m": w["mA"], "atoms": w["K"],
            "table2": w["post"] in table2, "file": rel,
            "sha256": hashlib.sha256(shipped).hexdigest(), "run_sha256_16": run16,
            "canonical_hash": verify.canonical_hash(cand), "script_sanitized": edited,
            "lambda_down7": down(Fraction(cand["scale"]), 7), "werner_v": w["werner_v"], "kg3_upper": w["kg3_upper"],
        })
        print(f"post {w['post']:>3}  {entries[-1]['set']:<10} {rel}{'  (script sanitized)' if edited else ''}", flush=True)
    index = {
        "about": "All certificates accepted by the frozen verifier during the two runs (paper, Appendix B), in order "
                 "of submission. file: path in this bundle (LZMA-compressed JSON, or plain JSON for the theorem's "
                 "certificate). sha256: of the uncompressed JSON as shipped. run_sha256_16: SHA-256 prefix of the "
                 "file as submitted during the run. canonical_hash: the verifier's hash of the certificate without "
                 "its metadata fields. script_sanitized: an absolute directory prefix was removed from the "
                 "free-text field 'script' (nothing else changed). lambda_down7 and werner_v are rounded down, "
                 "kg3_upper is rounded up. table2: the best certificate of its direction set (paper, Table II). "
                 "Check with: python3 -I check_release.py . --all",
        "certificates": entries,
    }
    (OUT / "INDEX.json").write_text(json.dumps(index, indent=1) + "\n")
    print(f"{len(entries)} entries, {sum(e['script_sanitized'] for e in entries)} sanitized, "
          f"{sum(e['table2'] for e in entries)} in Table II, {total / 1e6:.1f} MB compressed")


if __name__ == "__main__":
    main()
