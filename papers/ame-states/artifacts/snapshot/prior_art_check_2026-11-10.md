# Prior-art check of 2026-11-10 (answer to the exam of 2026-11-09, findings B1, M1, M3, M4, M5)

Everything below was fetched and read on 2026-11-10. "Read" means the full text unless a line says otherwise. Nothing
here was obtained from the authors or maintainers; nobody was contacted. Raw pages: `data/sources/codetables_2026-11-10/`
(with `INDEX.json`). Scripts: `paper/compare_codetables.py`, `paper/lc_invariant.py`, `work/exam20261109/`. In the release the same files are
`snapshot/sources/codetables_2026-11-10/`, `ancillary/codetables/compare_codetables.py`,
`ancillary/codetables/lc_invariant.py` and `ancillary/ame910/check_ame910.py`, and this file is `snapshot/prior_art_check_2026-11-10.md`.

## 1. codetables.de (M. Grassl, bounds on additive quantum codes)

The query parameter `q` of `QECC.php` is the size of the field of the additive code, so the local dimension is its
square root (`q=64`: local dimension 8). 41 pages were fetched with a plain GET; `INDEX.json` lists bounds,
construction type, "last modified" date and SHA-256 of each.

| cell | 10 Oct 2026 snapshot (`codetables_k0_2026-10-10.json`) | page on 10 Nov 2026 | construction shown | last modified |
|---|---|---|---|---|
| [[15,0]]_8 | 7-8 | **8-8** | type "Gra"; cyclic code of length 15 over GF(2^6), explicit generating polynomial and 45-row stabilizer matrix | **2026-10-21** |
| [[17,0]]_8 | 8-9 | **9-9** | type "Gajewski"; cyclic code of length 17, explicit polynomial and 51-row matrix | **2026-10-21** |
| [[16,1]]_8 | (not in the k = 0 snapshot) | 8-8 | | 2026-10-21 |
| [[21,0]]_8 | 9-11 | **10-11** | type "Gra"; cyclic code of length 21 with an explicit polynomial (not AME: AME(21,8) needs 11) | 2026-10-21 |
| [[16,0]]_8 | 8-9 | 8-9 (AME(16,8) undecided there) | stored matrix | 2015-02-09 |
| [[18,0]]_8, [[19,0]]_8, [[20,0]]_8 | 9-10, 9-10, 10-11 | unchanged | | 2015-02-09 |
| [[14,0]]_4 | 7-8 | 7-8 | explicit quasi-cyclic [[14,0,7]]_4 with stabilizer matrix | 2010-05-11 |
| [[16,0]]_4 | 8-8 | 8-8 | explicit quasi-cyclic [[16,0,8]]_4 | 2010-05-11 |
| [[13,0]]_4, [[15,0]]_4 | 6-7, 7-8 | unchanged (AME(13,4), AME(15,4) undecided there) | | 2010 |
| [[18,0]]_4 | 8-10 | 8-10 | | 2010-05-11 |
| [[14,0]]_3, [[16,0]]_3 | 6-7, 7-8 | 6-7, 7-8 | | 2020-11-16 |
| local dimension 9, 11, 13, 16 (`q` = 81, 121, 169, 256) | - | "currently not available" | | |

Consequences.
* **AME(15,8) and AME(17,8)** (findings B1). Both cells were undecided in the snapshot of 10 October 2026 and are blank
  in the Huber-Wyderka table (updated 8 June 2026; unchanged on 10 Nov per the exam). Our constructions were certified on
  2026-10-16 (candidate files `candidates/20261016-010218_additive_...json`, `20261016-010302_additive_...json`; repo
  commit cdcd1ea of 2026-10-16) and the manuscript containing both coefficient vectors was added to the project's public
  site in commit 9c46a78 (commit time 2026-10-16 23:05 -0700; the push time is not recorded in the repository). The
  code tables show codes with these parameters with the date 2026-10-21. **We do not know the origin of the two table
  entries** and make no statement about it.
* **Comparison of the codes** (`paper/compare_codetables.py`, output `results/codetables_comparison_2026-11-10.json`;
  own GF(8) arithmetic; the pages are read as data). Verified facts:
  * Both printed stabilizer matrices span F_2-spaces of dimension 3n that are isotropic and closed under
    multiplication by w for w^3 = w + 1 (and neither for w^3 = w^2 + 1): the codes are F_8-linear and cyclic.
  * Exchanging X and Z on every party puts the table's [[15,0,8]]_8 code in the form {(x | x G)} with G the symmetric
    circulant with c = (w^2, w^3, 0, w^3, w^4, w, 1) and zero diagonal, and the table's [[17,0,9]]_8 code in the same
    form with c = (1, w, w^6, 1, w^6, w, w, w). These are exactly the coefficients of the generating polynomials shown on
    the pages. Both are AME (all 6,435 resp. 24,310 cuts re-checked) and both vectors are in the solution lists of our
    exhaustive scans (scan orbit 20 of 28 for n = 15, scan orbit 2 of 4 for n = 17; orbits numbered by their
    lexicographically least normalized member).
  * **n = 17: the table's code is equal to the state of Theorem 1 up to a party permutation and local Clifford
    unitaries.** Explicit map, checked row by row: the local Clifford (a|b) -> (a + w^3 b | b) on every party (T = [1,3;0,1]
    in the integer labels), removal of the constant diagonal w^4 (label 6), the multiplier i -> 6i, and the scaling
    (a|b) -> (mu a | mu^-1 b) with mu = w^4 (label 6).
  * n = 15: the same-T-on-every-party local Cliffords (all 504 elements of SL(2,8)) take the table's code to scan
    orbits 8, 9, 12, 15, 16, 20, 24. The vector printed in Theorem 1 is in scan orbit 4. Under these maps the 28 orbits
    of n = 15 solutions form 4 classes of 7 orbits each, and the 4 orbits of n = 17 solutions form one class.
  * **n = 15: the table's code and the state of Theorem 1 are not equal up to local Clifford operations and party
    permutations** (`paper/lc_invariant.py`, output `results/lc_invariant_2026-11-10.json`). Invariant: for a stabilizer
    group S and a set B of 9 parties, S_B (elements supported in B) has F_2-dimension 9 and the nine S_{B-i} have
    dimension 3; the number of triples i, j, l with dim(S_{B-i} + S_{B-j} + S_{B-l}) < 9 is preserved by local
    symplectic maps, so the multiset of these numbers over the 5,005 sets B is the same for locally Clifford equivalent
    states (up to permutations). For the F_8-linear codes here it counts collinear triples among nine points of PG(2,8).
    The profile is constant on each of the 4 classes (all 28 orbit representatives computed) and the 4 class profiles are
    pairwise different. Hence the 336 circulant AME(15,8) solutions form exactly 4 classes under local Clifford
    operations and party permutations, the state of Theorem 1 is in one (scan orbits 4, 5, 11, 17, 21, 25, 26) and the
    table's code in another (8, 9, 12, 15, 16, 20, 24). The test does not address local unitaries outside the Clifford
    group. Single implementation, written on 2026-11-10; the consistency checks are that the profile is constant on
    the classes found independently by the uniform-Clifford computation.
  * All 48 circulant AME(17,8) solutions are one state up to local Cliffords and party permutations.
* **[[14,0,7]]_4 and [[16,0,8]]_4** (finding M1): explicit codes with these parameters have been in the table since
  2010-05-11. Our circulant and bordered-circulant codes are other explicit codes with the same parameters, not new
  parameters. Whether they are equivalent to the table's codes was not tested.
* **Local dimension 16**: the table has no entries (`q=256` "currently not available"), so it says nothing about
  AME(18..21,16). Local dimension 11 and 13 (ingredients AME(16,11), AME(16,13)): not available either.

## 2. AME(9,10): an earlier machine-checked coupling (finding M3)

* Source read line by line: `AtlasFCSolutions/Oqp35.lean` in github.com/niketp03/atlas-fc-verified at commit 15e4b3a
  (1,244 lines; commit date 2026-09-14, GitHub API). Its header and the repository README say the proof is
  `Atlas/FC/OpenQuantumProblem35_ame_9_10_open_solution.lean` of github.com/facebookresearch/atlas-lean (branch
  `prepare-atlas-v2`, commit 223d8bc), with local copies of the Formal Conjectures definitions removed. The GitHub API
  dates the first commit of that file in atlas-lean to **2026-08-28** (commit 8a68ff3). The atlas-lean original itself
  was not read. Kitamura's `ame-11-10-lean` (commits of 2026-09-22) cites the atlas-fc-verified file as "an earlier
  AME(9,10) formalization".
* What is proved: `ame_9_10_open : answer(True) <-> ExistsAME 9 10` for the Formal Conjectures statement
  (google-deepmind/formal-conjectures at 990d209a, `OpenQuantumProblems/35.lean`: a normalized state of 9 parties of
  dimension 10 whose reduction to floor(9/2) = 4 parties is maximally mixed after every permutation). Axioms reported by
  the repository: propext, Classical.choice, Quot.sound. We did not compile it.
* The state: each party is Z_2 x Z_5; the amplitude is (2^11 5^9)^(-1/2) times the sum over a in Z_2^2 of
  (-1)^(graph phase of the 9-cycle + label(a).u) zeta_5^(graph phase of G on (q5Logical(a), v)), where label(a) =
  a_0 l_1 + a_1 l_2 with l_1 = 110110110, l_2 = 101101101, G is the symmetric circulant on Z_10 over Z_5 with
  c = (3,1,1,0,0), and q5Logical(a) = a_0 + 2 a_1 is the value of vertex 0 of G.
* Relation to our Lemma 3: it is a coupling with a label space of dimension K = 4: four states Z^label(a)|C_9> of nine
  qubits whose uniform mixture has maximally mixed 4-party marginals (the common stabilizer has minimum weight 5; Lean
  lemma `q2Cycle_subgroup_minWeight`), glued to the AME(10,5) graph state of G with one party restricted to four of its
  five values. It is not an instance of the lemma as stated in revisions <= 4.1 (K = a = 2 would need a 4-uniform state
  of 10 qubits, which does not exist), but it is an instance of the general form now stated as Lemma 3 (labels
  t < K <= b, unit vectors |A_t> whose uniform mixture is k-uniform). It predates Kitamura's AME(11,10).
* Independent re-check of both ingredients (`work/exam20261109/check_ame910.py`): G[A^c, A] is invertible mod 5 for all
  252 five-sets; the 127 nonzero common stabilizers of the four 9-qubit states have minimum weight 5.
* Effect on the paper: AME(9,10) is known. `data/status.json` (pinned, not changed) has it as `open`;
  `paper/make_tables.py` now applies a dated correction (AME(9,10) exists, and by closure AME(9,50) = AME(9,5) x AME(9,10))
  before counting. Table I counts: 119 existing / 905 not covered (was 118 / 906); beyond-table 948 / 1,119 (was 947 /
  1,120). None of the 30 phase-3 cells is affected (no n = 9 cell is claimed).

## 3. Bibliographic data confirmed through Crossref and the arXiv API (finding M4)

| reference | confirmed data |
|---|---|
| Ball, Moreno, Simoens, "Stabilizer codes over fields of even order" | IEEE Trans. Inf. Theory 71(5), 3707-3718 (2025), DOI 10.1109/TIT.2024.3454480; arXiv:2401.06618 (v2 of 2024-09-09) |
| Ball, Zhang, "Error-correcting codes and AME states for mixed-dimensional Hilbert spaces" | Phys. Rev. A 113, 012432 (2026), DOI 10.1103/yqck-92fb; arXiv:2510.17231 |
| Shen, Chen, "AME states in tripartite heterogeneous systems" | Quantum Inf. Process. 20, 111 (2021), DOI 10.1007/s11128-021-03034-y; arXiv:2001.08821 |
| Zhao, Ma | arXiv:2602.06377v1 "Hermitian self-dual generalized Reed-Solomon codes" (2026-02-06). Crossref lists, by the same two authors, "Complete classification and construction of Hermitian self-dual generalized Reed-Solomon codes", Finite Fields Appl. 119, 102952 (issue dated March 2027), DOI 10.1016/j.ffa.2026.102952. The titles differ; that the article is the published version of the preprint is likely but was not checked against the article's text. |
| Zhao, Ma, Yan, Sun, "Hermitian self-dual MDS twisted generalized Reed-Solomon codes" | Entropy 28(11), 1169 (published 2026-10-22), DOI 10.3390/e28111169 (Crossref) |
| AME(16,11) | the Huber-Wyderka snapshot credits AME(16,11) (and (15,11), (14,11), (12,9), (14,9), (13,9), (11,9)) to Grassl-Roetteler, arXiv:1502.05267 |

## 4. Full-text reading for the q = 16 cells and the lemma (finding M5)

| source | what was read | result |
|---|---|---|
| Zhao, Han, Ma, Yan, Sun, arXiv:2508.09687v1 (13 Aug 2025), "Hermitian self-dual twisted generalized Reed-Solomon codes" (29 pages) | abstract, all definitions, theorem, corollary and construction statements of Secs. 3-5 and the conclusion (PDF text); proofs not checked | General criteria (Theorems 1, 5, 7) without explicit codes, and five explicit constructions over F_{q^2}: I (n divides q-1), II (q = 2^s, n = q), III (n = 2 sqrt(q)), IV (q odd, n = q+1), V (q odd, n divides q-1; the only one shown to be MDS). Every explicit construction has n <= q+1, i.e. n <= 17 for q = 16 and n <= 9 for q = 8. It gives no Hermitian self-dual MDS code of length 18 or 20 over GF(256). |
| Zhao, Ma, Yan, Sun, Entropy 28(11), 1169 (22 Oct 2026) | **not read**: mdpi.com answers HTTP 403 to our fetches. Title and authors (Crossref) match the preprint above except for one author (Han) | Not checked. That the article contains no construction beyond the preprint is an assumption. |
| Hu, Ji, arXiv:2609.09943v1 (9 Sep 2026), "Quantum MDS codes from complements of unions of finite-field subsets" (37 pages) | abstract, introduction, the standing assumption of Sec. 2 and the statement headers of all theorems (PDF text); proofs not read | "Throughout this paper, q is an odd prime power"; GRS codes only. Not applicable to q = 8 or 16 (and GRS codes are covered by the Ball-Vilar bound in any case). |
| Hontarenko, Zyczkowski, arXiv:2609.11796 | v1 is still the only version (arXiv API, 2026-11-10). Full text (PDF text) searched for even q, characteristic two, q = 8, q = 16; Table I and Sec. S6 read | Odd prime powers only; Sec. S6 and Lemma S4: the direct characteristic-two analogue is obstructed for q = 2^r >= 8. No q = 8 or q = 16 state. |
| Hontarenko, Zyczkowski, arXiv:2610.17170 | v1 is still the only version (arXiv API, 2026-11-10). Sec. 7.4 | The sentence "A forthcoming revision will also include an extension to stabilizer AME(N, q) states for all 2 <= N <= 2q" is unchanged; no range of q, no construction. |
| Ball, Vilar, arXiv:2106.10180 | full text (project copy, Theorem 4.3), as before | bound used in Sec. IV is quoted correctly (confirmed by the exam) |
| Wan, Zhu, arXiv:2302.06169 (found by a web search on 2026-11-10) | abstract only | quantum MDS codes with d > q/2+1 from Hermitian self-orthogonal GRS codes; GRS, hence subject to the Ball-Vilar bound. Not read in full. |
| Shen-Chen; Raissi et al. (related gluing constructions) | abstracts only, as before | unchanged |
| web searches (2026-11-10) for AME(18,16), AME(20,16), Hermitian self-dual MDS codes of length 18, 20 over GF(256), AME(15,8), AME(17,8) | result pages | no source found for the q = 16 cells; for q = 8 only the code-table entries above |

Verdict. Nothing read gives a Hermitian self-dual MDS code over GF(256) of length 18 or 20, a quantum MDS code
[[19,1,10]]_16 or [[21,1,11]]_16, or an AME(n,16) state with 18 <= n <= 21. The statement remains "not found in what we
read", with these gaps: the Entropy article (not reachable), Wan-Zhu, Shen-Chen and Raissi et al. (abstracts only), and
the announced revision of arXiv:2609.11796 (not yet posted).
