# Record snapshots (what "known", "nonexistent" and "not covered" mean in the paper)

The verifier decided each verdict against a frozen snapshot of the existence status of AME(n,d). Three versions were used.
Each is pinned in the verifier by the first 16 hex digits of its SHA-256.

| file | used in | SHA-256 prefix (verifier pin) | content |
|---|---|---|---|
| `status_phase1.json` | phase 1 (2026-10-16 00:59-04:58) | `e07e43e502967e93` | 2 <= d <= 11: the Huber-Wyderka table + the 2026 overlay below + closure |
| `status_phase2.json` | phase 2 (2026-10-16/17) | `3880ea0daf6f66cc` | as phase 1, plus AME(15,8), AME(17,8) (this work, phase 1) and AME(9,6) (Kitamura) |
| `status_phase3.json` | phase 3 (2026-10-18/19) | `2473d4938ca44314` | as phase 2, plus AME(11,10) (Kitamura), and cells 12 <= d <= 64, n <= 40 from published families only (`beyond_table: true`) |

`build_status.py` is the program that built them (phase-3 version; it documents every rule). `sources/` holds the raw inputs:
the Huber-Wyderka table data `huber_wyderka_table_2026-06-08.html` (the page's `table.html`, retrieved 2026-10-16 from
github.com/huberfe/huberfe.github.io, branch gh-pages; live at https://huberfe.github.io/ame.html, page "updated:
08.06.2026") and `codetables_k0_2026-10-10.json` (Grassl's [[n,0,d]]_q bounds for q <= 8, retrieved 2026-10-10).

## Status of every source

| source | kind | used for |
|---|---|---|
| Huber-Wyderka, Table of AME states (online, updated 8 June 2026) | curated online table | all cells with d <= 11; blank cells are taken as unknown |
| Reed-Solomon / MDS codes (minimal-support AME for n <= q+1, q a prime power) | classical theorem | d > 11 |
| Orthogonal Latin squares of every order d != 2, 6 (Bose, Shrikhande, Parker 1960); golden AME(4,6) (Rather et al., PRL 128, 080507 (2022)) | refereed | AME(4,d), d > 11 |
| Scott bound (PRA 69, 052330 (2004)), shadow inequalities (J. Phys. A 51, 175301 (2018)), AME(7,2) (PRL 118, 200502 (2017)) | refereed | nonexistence (through the table) |
| Shi, Zhang, Zhao, Li, arXiv:2608.01011 (AME(7,d), d >= 3) | preprint, unrefereed | n = 7 |
| Bevins, Bidav, arXiv:2608.05781 (AME(12,5), (18,11), (18,13) and projections) | preprint, unrefereed | those cells |
| Hontarenko, Zyczkowski, arXiv:2609.11796 (AME(2q,q), AME(2q-1,q), AME(4k+1/4k+2,q), AME(q+3,q) for odd prime powers q; Table 1 parents) | preprint, unrefereed | those cells |
| Kitamura, github.com/KitaKen1/ame-9-6-lean and ame-11-10-lean (Lean 4 proofs) | repositories with machine-checked proofs, unrefereed | AME(9,6) (phase 2, 3), AME(11,10) (phase 3) |
| Hontarenko, Zyczkowski, arXiv:2610.17170 Sec. 7: "a forthcoming revision will also include an extension to stabilizer AME(N,q) states for all 2 <= N <= 2q" | **announcement** (no construction, no q-range) | **not used** to mark any cell; mentioned in the paper (Table I caption) |
| Ball, Moreno, Simoens, arXiv:2401.06618 (no [[8,0,5]]_4 stabiliser code) | preprint (a journal version is claimed by the authors' repository, not confirmed) | annotation of AME(8,4) in phase 3 (the cell stays unknown: non-stabilizer states are not covered) |
| this work (AME(15,8), AME(17,8)) | phase-1 results | phases 2 and 3 |

Closure rules: tensor products on the same parties (AME(n,d1) and AME(n,d2) give AME(n,d1 d2)); projection
AME(2k,d) -> AME(2k-1,d); and, contrapositively, AME(2k-1,d) nonexistent gives AME(2k,d) nonexistent.

**Wording.** A cell "not covered by the snapshot" (status `open` in the files) is one that none of these sources gives.
This is a statement about this snapshot, not about the whole literature. Priority of every new cell needed its own
literature search. `literature_log_2026-10-19.md` is the auditing agent's dated log of that search, copied unchanged
from `results/literature.md` of the project repository as of the phase-3 commit (last entry 2026-10-19). It records
which version of each preprint was read (arXiv:2608.01011v1, 2608.05781v2, 2609.11796v1, 2610.17170v1, 2401.06618v2)
and when; the paper cites these versions. It is a working log (agent shifts, board post numbers, repository paths), not
a curated bibliography; the paper states the limits of the search.

Check the pins: `shasum -a 256 status_phase*.json` (the first 16 hex digits are the pins above).
