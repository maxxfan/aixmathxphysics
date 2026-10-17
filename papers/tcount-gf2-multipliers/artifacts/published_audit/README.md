# Audit of published values

`AUDIT.md` explains the method and results. The audited matrices are Fisher et al.'s released data (MIT license), not redistributed here. Download them from
* finals: https://github.com/DanilkaFish/gigaevo-core-vartodd/tree/d9e411d3e8abbe64f5b8aa7f3a30f72a497975af/evolution_results
* starts: https://github.com/DanilkaFish/VarTodd/tree/409788270f9b46fe8deb95fb849d5ee5e4073d51/data/init_npy
and compare with `SHA256SUMS` (file names there use '_' for '^'). In `output.txt` the decisive column is 'singular a' (a GF(2^n) multiplication tensor has none in any basis); the header line explains the others. Then run `python3 -I audit_npy.py FILE.npy` (numpy, loads with allow_pickle=False). `certify_npy.py` needs the project's qlib/basis_map.py and verifier. `qcsim_full.py` simulates .qc circuits (used on Vandaele's public circuits).
