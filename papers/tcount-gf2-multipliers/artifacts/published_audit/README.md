# Audit of published values

`AUDIT.md` explains the method and the results. `output.txt` is the raw output of the scripts.

The audited matrices are Fisher et al.'s released data (MIT license). They are not redistributed
here. Download them from

* finals: https://github.com/DanilkaFish/gigaevo-core-vartodd/tree/d9e411d3e8abbe64f5b8aa7f3a30f72a497975af/evolution_results
* starts: https://github.com/DanilkaFish/VarTodd/tree/409788270f9b46fe8deb95fb849d5ee5e4073d51/data/init_npy

and compare them with `SHA256SUMS` (the file names there use `_` for `^`).

## Scripts

All scripts load the matrices with `allow_pickle=False`. They need numpy, except `qcsim_full.py`.
Run them from the bundle root.

| script | what it does |
|---|---|
| `audit_npy.py FILE.npy ...` | Negative test. Prints, for each matrix, whether singles and pairs are even, whether the odd triples have the (a, b, c) form, and the number of nonzero a for which b ↦ M(a, b) is singular (all a for n ≤ 16, 4000 random a otherwise). A GF(2^n) multiplication tensor has none, in any basis. In `output.txt` this is the column `singular a`; the header line explains the others. |
| `certify_npy.py FILE.npy [OUT.json]` | Positive test. Maps the matrix to the polynomial basis with `qlib/basis_map.py` (shipped here, with `qlib/gf2mult.py`) and runs the frozen verifier `../verifier/verify.py` on the result. This is the check behind "the Khoruzhii-start finals pass the verifier". For n = 10 the mapping takes a few minutes. |
| `redmap.py FILE.npy ...` | Prints the image of each power x^s under the released map (used for the examples in `AUDIT.md`). |
| `qcsim_full.py FILE.qc n MODULUS_BITS` | Simulates a `.qc` circuit (used on Vandaele's public circuits). |

Example:

```bash
python3 published_audit/audit_npy.py   path/to/gf2_5_49.npy          # singular a: 12/31
python3 published_audit/certify_npy.py path/to/gf2_8_khor_123.npy    # MATCHES_RECORD ... 123
```

`certify_npy.py` runs the frozen verifier, which compiles `verifier/phasecheck` next to its
source (a C compiler `cc` is needed). `check_release.py` ignores that build product. On a machine
other than the one that froze the verifier, the verifier reports `pins: MODIFIED: phasecheck`,
because the pinned hash of the compiled binary is machine-specific. The pins of `verify.py`,
`phasecheck.c` and `records.json` must match; `certify_npy.py` checks this and says so.
