"""hashes.py (release/frames; standard library only): SHA-256 'count hashes' of the canonical result files.

  hashes.py write DIR > HASHES.txt        one line '<sha256>  <name>' per DIR/*.txt (sorted, HASHES.txt itself excluded)
  hashes.py verify HASHES.txt DIR         exit 0 iff every listed file exists in DIR with exactly that hash (prints one line per file)
  hashes.py one FILE                      print the hash of one file
"""
import hashlib
import os
import sys


def h(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main(a):
    if a[0] == "write":
        for n in sorted(os.listdir(a[1])):
            if n.endswith(".txt") and n != "HASHES.txt":
                print("%s  %s" % (h(os.path.join(a[1], n)), n))
    elif a[0] == "verify":
        bad = 0
        for line in open(a[1]):
            if not line.strip() or line.startswith("#"):
                continue
            want, name = line.split()
            p = os.path.join(a[2], name)
            got = h(p) if os.path.exists(p) else "MISSING"
            ok = got == want
            bad += not ok
            print("%s %s %s" % ("OK  " if ok else "BAD ", name, got))
        sys.exit(1 if bad else 0)
    elif a[0] == "one":
        print(h(a[1]))


main(sys.argv[1:])
