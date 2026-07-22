# Reproducibility artifact

Ancillary files for: *Properly colored paths in the union of two Hamiltonian paths:
the O(√n) collapse, the transit floor, and the failure of block lower bounds*
(Lizardi, July 2026).

On arXiv these files ship in the standard ancillary directory `anc/`; in the source
repository they live in `code/`. The two are the same files — the paper refers to them
as the ancillary artifact and cites individual logs by name.

Environment used to produce `logs/`: Python 3.14, Windows 10, single core; no external
Python dependencies for the core verifiers (standard library only). One optional
script, `cpsat_rho.py`, uses Google OR-Tools (`pip install ortools`) to certify
optimality of the staircase-chain values by an independent constraint-programming
model. The C++ tools compile with `g++ -O2` (`rho_tool.cpp` additionally uses
`-fopenmp`).

## Claim-to-log map (paper section "Computational verification")

The item numbers match the numbered list of verified statements in that section.

| paper item | claim | script | log |
|---|---|---|---|
| 1 | anchors of the half-bound section: ρ(σ_m)=2m+3 (m=1..4); ρ(σ*)=10 at n=18; ρ_min(7)=5, ρ_min(8)=6 | `pc.py`, `search_small.py` | `validate2.log`, `search_small.log` |
| 2 | classification of admissible visits of H_{T_b} is exactly (V1)–(V5), even b ≤ 14; M(T_b)=b for odd b ≤ 11 | `verify_M.py` (minimal independent verifier, written directly from the admissible-visit/transit-number definition, no shared code) | `verify_M.log` |
| 3 | M(T_b)=3 for every even 6 ≤ b ≤ 100; M(T_b)=b for odd b ∈ {7,9,11,13} | `check_all_even.py`, `validate2.py` | `check_all_even.log`, `validate2.log` |
| 4 | min over S_b of M(π) = 3 for 4 ≤ b ≤ 9 (full enumeration) | `search_small.py` | `search_small.log` |
| 5 | explicit lower-bound path (Prop. "Explicit long path") valid at 7 points incl. (20,20) | `verify_path.py` | `verify_path.log` |
| 6 | exact ρ(σ_{a,b}) values up to n = 720 (with reachable-state counts) | `endtoend.py`, `bigtable.py` | `endtoend.log`, `bigtable.log` |
| 7 | joint (b(σ), ρ(σ)) distribution over all of S_n, 7 ≤ n ≤ 12 (C++; 479,001,600 permutations at n = 12) | `cpp/rho_tool.cpp` | `exh7.txt` … `exh12.txt` |
| 8 | landscape probe (exploratory) | `probe.py` | (exploratory; rerun to reproduce) |
| 9 | transit-floor construction (Theorem "Floor lower bound"): the designated matching-union witness is a valid admissible visit of size ≥ 3, validated from the definition — exhaustive all π ∈ S_b for 3 ≤ b ≤ 9, samples to b = 50, 0 failures | `verify_floor.py` | `verify_floor.log` |
| 10 | odd gadget T′_b (Appendix "The odd gadget"): inventory/highway/crossing claims for odd b ∈ {7,…,21}, M(T′_b)=3, with positive and negative controls; transit spectrum for all π ∈ S_b, b ≤ 8 (μ(2)=2, μ(b)=3 for 3 ≤ b ≤ 8, value 4 unattained at odd b) | `verify_oddb.py`, `verify_spectrum.py` | `verify_oddb.log`, `verify_spectrum.log` |
| 11 | the n = 14 threshold: σ† block-free with ρ(σ†) = 13; all 831,283,558 block-free σ ∈ S_13 properly-colored Hamiltonian; exactly 16 block-free failures at n = 14 | `verify_n14.py`; `cpp/bfcensus.cpp` (n = 13 census, hours of CPU); `n14_doublebad.py` + `word_bruteforce.py` (the exact n = 14 count) | `verify_n14.log`, `logs/blockfree13.txt`, `n14_doublebad.log` |
| 12 | staircase chains: θ_s block-free, junction structure, parity lemma (s ≤ 33, k ≤ 6); exact ρ(Θ_{9,k}) = 9k for k ≤ 3 (reference solver); certified-optimal ρ = 36, 38, 40, 42 at k = 4..7 and 60, 62 for s = 15 (CP-SAT, independent formalism) | `verify_stair.py`; `cpsat_rho.py` (needs `ortools`) | `verify_stair.log`, `cpsat_stair.txt` |
| App. A | trap-chain exact values ρ(σ_C) = 2k+2(L−1) (cross-check at 6 points) | `check_trap.py` | `check_trap.log` |
| open-problems §(Lovász side) | G(ab, σ_{a,b}) has a path with ≥ a(b−3) = n − 3a matching edges (explicit snake, edge-by-edge check at 6 points up to n = 600) | `verify_f_linear.py` | `verify_f_linear.log` |

## How to reproduce

```
python verify_M.py        # ~40 s   independent verifier: classification + M values
python verify_path.py     # <1 s    explicit lower-bound path, edge-by-edge
python validate2.py       # ~5 s    solver cross-validation + odd-b parity + two-visit total
python check_trap.py      # ~30 s   trap-chain exact values (Appendix A)
python endtoend.py        # ~1 min  exact rho of inflations vs bound, n ≤ 80
python check_all_even.py  # ~1 min  M(T_b)=3 for all even b in [6,100]
python search_small.py    # ~2 min  rho_min(7,8) + exhaustive min M over S_b, b ≤ 9
python bigtable.py        # ~1 min  exact rho table up to n = 720
python verify_n14.py      # <1 s    the n = 14 block-free witness, exact rho
python verify_stair.py    # ~1 min  staircase structure + exact rho, k ≤ 3
python verify_oddb.py     # ~2 min  odd-gadget structural confirmation + controls
python n14_doublebad.py   # ~2 min  the exact count of 16 failures at n = 14
python verify_spectrum.py # ~20 min transit spectrum, all of S_b for b ≤ 8
python verify_floor.py    # ~30 min transit-floor construction, exhaustive + samples
python cpsat_rho.py "..." # optional (ortools): certify a rho value optimal
python probe.py           # ~5 min  landscape probe (exploratory)
```

The n = 13 block-free census (`cpp/bfcensus.cpp`, compile `g++ -O2`, run
`bfcensus 13`) takes hours of CPU time; its result line is `logs/blockfree13.txt`.

`MANIFEST.sha256` lists checksums of all code and logs.

## Code inventory

- `pc.py` — reference solvers: ρ(σ), constrained PC-path search, transit number M
  (plain recursive).
- `family.py` — memoized iterative (state-dedup) solvers; the gadget family T_b.
- `inflate.py` — the inflation σ = τ[π_1..π_a]; the two-visit total (`M2`).
- `endtoend.py` — exact longest-PC-path search with state dedup (used at scale).
- `verify_M.py` — **minimal independent verifier** (~100 lines, no shared code):
  brute-force enumeration of all admissible visits, directly from the definition.
- `verify_path.py` — mechanical check of the explicit lower-bound path
  (Proposition "Explicit long path: the family is Θ(a+b)").
- `verify_floor.py` — implements the transit-floor proof's designated matching-union
  construction and validates every produced witness from the definition of an
  admissible visit (independent of the construction logic).
- `verify_oddb.py` — independent brute-force confirmation of the odd-gadget
  classification (Appendix "The odd gadget"), with positive controls (reproduces the
  even-gadget values) and negative controls (the M = 2 floor at b = 2).
- `verify_spectrum.py` — exhaustive M(π) over all of S_b, b ≤ 8: the transit floor
  at small sizes and the spectrum (the odd-b gap at the value 4).
- `verify_n14.py` — the n = 14 threshold witness σ†: block-freeness and exact ρ = 13.
- `n14_doublebad.py`, `word_bruteforce.py` — the exhaustive n = 14 census: exactly 16
  block-free permutations are not properly-colored Hamiltonian (via the matching
  dichotomy on doubled-pair words; `word_bruteforce.py` derives the unique bad word).
- `verify_stair.py` — the staircase pattern θ_s and chains Θ_{s,k}: block-freeness,
  junction structure, the parity lemma (edge by edge), and exact ρ for small chains.
- `cpsat_rho.py` — independent CP-SAT (arc-MILP, MTZ) model certifying ρ values
  optimal; a formalism unrelated to the state-space search. Requires `ortools`.
- `verify_f_linear.py` — explicit matching-edge-rich snake in G(ab, σ_{a,b}),
  edge-by-edge (the Lovász-side remark in the open-problems section).
- `cpp/rho_tool.cpp` — independent C++ implementation from a parallel verification
  effort; source of the exhaustive n ≤ 12 statistics (`exh*.txt`) and of cross-checks
  of the exact values at (a,b) ∈ {(6,8),(5,10),(4,12)} and the trap-family values.
- `cpp/bfcensus.cpp` — the block-free census at odd n via the matching dichotomy
  (exhaustive pruned enumeration; prints any non-Hamiltonian σ found).

## Algorithm note

The exact search explores the state graph on triples (visited set, current vertex,
last color), each state visited once; no pruning. Exactness is by construction; cost is
the number of reachable states, which the logs report per instance (57,066 at n = 80;
5,012,186 at n = 600; 7,273,586 at n = 720 for σ_{a,b}; 1,938 for the n = 14
witness). On unstructured permutations the state space is exponential and the method
was used only for n ≤ 19 (plus the structured staircase chains through n = 27).
