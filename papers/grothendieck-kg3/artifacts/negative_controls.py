#!/usr/bin/env python3
"""Negative controls for check_release.py: corrupt copies of the certificate (and of the bundle) and require the
checker to reject every one.   usage: python3 -I negative_controls.py [BUNDLE_DIR]      (about 3-4 minutes)

Cheap controls fail in parsing; the others change the local model, the target or the shrinking data only slightly
and must fail the exact residual test (the certificate's margin is about 1e-9).
"""
import copy
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_release as cr  # noqa: E402

bundle = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE
CERT = bundle / "certificate" / "KG3_upper_1.4494177506_rad721.json"
base = json.loads(CERT.read_text())
heavy = max(range(len(base["strategies"])), key=lambda k: int(base["strategies"][k][0]))


def frac(s):
    q = Fraction(s)
    return f"{q.numerator}/{q.denominator}"


def c_scale(f):
    def g(c):
        c["scale"] = frac(Fraction(c["scale"]) * f)
    return g


def c_float(c):
    c["alice"][0][0] = float(Fraction(c["alice"][0][0]))


def c_half_den(c):
    c["weight_den"] = str(int(c["weight_den"]) // 2)


def c_mask_beyond(c):
    c["strategies"][0][1] = format(int(c["strategies"][0][1], 16) | (1 << len(c["alice"])), "x")


def c_not_perm(c):
    c["transforms"][1][0][0] = c["transforms"][1][0][1]


def c_drop_vector(c):
    c["alice"] = c["alice"][:-1]
    c["bob"] = "alice"


def c_flip_heavy(c):
    s = c["strategies"][heavy]
    s[1] = format(int(s[1], 16) ^ 1, "x")


def c_drop_atoms(c):
    order = sorted(range(len(c["strategies"])), key=lambda k: -int(c["strategies"][k][0]))
    drop = set(order[:100])
    c["strategies"] = [s for k, s in enumerate(c["strategies"]) if k not in drop]


def c_swap_perm(c):
    p = c["transforms"][5][0]
    p[3], p[4] = p[4], p[3]


def c_flip_sign(c):
    t = c["transforms"][7]
    t[1] = format(int(t[1], 16) ^ (1 << 10), "x")


def c_move_vector(c):
    c["alice"][0] = [frac(Fraction(x) * Fraction(1000001, 1000000)) for x in c["alice"][0]]


def c_no_transpose(c):
    c["transpose"] = False


CONTROLS = [
    ("scale = 0", lambda c: c.__setitem__("scale", "0")),
    ("a JSON float coordinate", c_float),
    ("denominator halved (total weight > 1)", c_half_den),
    ("sign mask with a bit beyond m", c_mask_beyond),
    ("a transform that is not a permutation", c_not_perm),
    ("one direction removed", c_drop_vector),
    ("scale raised by a factor 1 + 1e-9", c_scale(Fraction(1000000001, 1000000000))),
    ("scale lowered by a factor 1 - 1e-6", c_scale(Fraction(999999, 1000000))),
    ("one sign of the heaviest atom flipped", c_flip_heavy),
    ("the 100 heaviest atoms removed", c_drop_atoms),
    ("two entries of one permutation swapped", c_swap_perm),
    ("one transform sign flipped", c_flip_sign),
    ("one direction lengthened by 1e-6", c_move_vector),
    ("transposition flag cleared", c_no_transpose),
]


def main():
    bad = 0
    with tempfile.TemporaryDirectory(prefix="kg3neg") as td:
        for name, f in CONTROLS:
            c = copy.deepcopy(base)
            f(c)
            p = pathlib.Path(td) / "cert.json"
            p.write_text(json.dumps(c))
            try:
                cr.check(p)
                bad += 1
                print(f"NOT REJECTED  {name}", flush=True)
            except cr.Fail as e:
                print(f"rejected      {name}: {e}", flush=True)
            except Exception as e:  # malformed input may also surface as a Python error; that is a rejection
                print(f"rejected      {name}: {type(e).__name__}: {e}", flush=True)
        # bundle tampering: one changed byte in a listed file must fail the manifest check
        tb = pathlib.Path(td) / "bundle"
        shutil.copytree(bundle, tb, ignore=shutil.ignore_patterns("__pycache__", "accum"))
        readme = tb / "README.md"
        readme.write_bytes(readme.read_bytes() + b" ")
        try:
            cr.check_manifest(tb)
            bad += 1
            print("NOT REJECTED  README.md changed after the manifest", flush=True)
        except cr.Fail as e:
            print(f"rejected      README.md changed after the manifest: {e}", flush=True)
    n = len(CONTROLS) + 1
    print(f"{n - bad} of {n} corruptions rejected" + ("" if bad == 0 else f"; {bad} NOT REJECTED"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
