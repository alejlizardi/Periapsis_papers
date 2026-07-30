# Reproducibility artifact

Ancillary files for: *Properly colored paths in the union of two Hamiltonian paths:
the O(√n) collapse, the transit floor, and the failure of block lower bounds*
(Lizardi, July 2026).

These files accompany the paper as its ancillary artifact; in the source repository
they live in `code/`. The paper cites individual logs by name.

Environment used to produce `logs/`: Python 3.14, Windows 10, single core; no external
Python dependencies for the core verifiers (standard library only). Two optional
scripts, `cpsat_rho3.py` and `fmu_cpsat.py`, use Google OR-Tools (`pip install
ortools`) to certify optimality of the three-color witnesses and the exact max-μ
values by an independent constraint-programming model. The C++ tools compile with
`g++ -O2` (`rho_tool.cpp` additionally uses `-fopenmp`).

## Claim-to-log map (paper section "Computational verification")

The item numbers match the numbered list of verified statements in that section.

| paper item | claim | script | log |
|---|---|---|---|
| 1 | census anchors: ρ_min(7)=5, ρ_min(8)=6 by full enumeration of S_7, S_8 (Python solvers, independent of the C++ censuses of item 7); see the ρ_min census appendix | `pc.py`, `search_small.py` | `search_small.log` |
| 2 | classification of admissible visits of H_{T_b} is exactly (V1)–(V5), even b ≤ 14; M(T_b)=b for odd b ≤ 11 | `verify_M.py` (minimal independent verifier, written directly from the admissible-visit/transit-number definition, no shared code) | `verify_M.log` |
| 3 | M(T_b)=3 for every even 6 ≤ b ≤ 100; M(T_b)=b for odd b ∈ {7,9,11,13} | `check_all_even.py`, `validate2.py` | `check_all_even.log`, `validate2.log` |
| 4 | min over S_b of M(π) = 3 for 4 ≤ b ≤ 9 (full enumeration) | `search_small.py` | `search_small.log` |
| 5 | explicit lower-bound path (Prop. "Explicit long path") valid at 7 points incl. (20,20) | `verify_path.py` | `verify_path.log` |
| 6 | exact ρ(σ_{a,b}) values up to n = 720 (with reachable-state counts) | `endtoend.py`, `bigtable.py` | `endtoend.log`, `bigtable.log` |
| 7 | joint (b(σ), ρ(σ)) distribution over all of S_n, 7 ≤ n ≤ 12 (C++; 479,001,600 permutations at n = 12) — source of the ρ_min census appendix and the deficient-set counts \|D(n)\| | `cpp/rho_tool.cpp` | `exh7.txt` … `exh12.txt` |
| 8 | landscape probe (exploratory) | `probe.py` | (exploratory; rerun to reproduce) |
| 9 | transit-floor construction (Theorem "Floor lower bound"): the designated matching-union witness is a valid admissible visit of size ≥ 3, validated from the definition — exhaustive all π ∈ S_b for 3 ≤ b ≤ 9, samples to b = 50, 0 failures | `verify_floor.py` | `verify_floor.log` |
| 10 | odd gadget T′_b (Appendix "The odd gadget"): inventory/highway/crossing claims for odd b ∈ {7,…,21}, M(T′_b)=3, with positive and negative controls; transit spectrum for all π ∈ S_b, b ≤ 8 (μ(2)=2, μ(b)=3 for 3 ≤ b ≤ 8, value 4 unattained at odd b) | `verify_oddb.py`, `verify_spectrum.py` | `verify_oddb.log`, `verify_spectrum.log` |
| App. A | trap-chain exact values ρ(σ_C) = 2k+2(L−1) (cross-check at 6 points) | `check_trap.py` | `check_trap.log` |
| open-problems §(matching side) | G(ab, σ_{a,b}) has a path with ≥ a(b−3) = n − 3a matching edges (explicit snake, edge-by-edge check at 6 points up to n = 600) | `verify_f_linear.py` | `verify_f_linear.log` |
| open-problems §(more colors) | ρ_min^(3)(n) = n for 4 ≤ n ≤ 10; ρ_min^(3)(11) = 10, unique deficient triple up to symmetry (coset-triangle census over the deficient set D(n)); all 12 witnesses independently certified by CP-SAT | `k3.cpp` (census; compile `g++ -O2`); `cpsat_rho3.py` (needs `ortools`) | (rerun to reproduce; census outputs archived in the research repo) |
| open-problems §(matching side) | exact f(n): f = n for n ≤ 6, f(7) = 6, f(8) = 7, f(9) = 8; max_π μ ∈ {n−1, n} for all σ, n ≤ 9, with {max μ = n−1} = D(n); exact max_π μ on σ_{a,6} (a = 2..5), σ_{3,8}, and the trap chains (one edge lost per piece beyond the second) | `fmu.cpp` (bitmask DP; validated vs direct search on all S_n, n ≤ 7, and vs the conversion inequality pointwise); `fmu_cpsat.py` (CP-SAT, optimality certified; needs `ortools`); `trap_snake.py` (edge-by-edge snake verifier, sizes to n = 5201) | (rerun to reproduce) |

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
python verify_oddb.py     # ~2 min  odd-gadget structural confirmation + controls
python verify_spectrum.py # ~20 min transit spectrum, all of S_b for b ≤ 8
python verify_floor.py    # ~30 min transit-floor construction, exhaustive + samples
python trap_snake.py      # ~1 min  edge-by-edge chain-snake verifier (max-mu)
python cpsat_rho3.py      # optional (ortools): certify the 12 three-color witnesses
python fmu_cpsat.py       # optional (ortools): certify exact max-mu values
python probe.py           # ~5 min  landscape probe (exploratory)
```

The three-color census (`k3.cpp`) and the max-μ bitmask DP (`fmu.cpp`) compile with
`g++ -O2`; their outputs back the k-color and matching-content items of the
open-problems section.

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
- `verify_f_linear.py` — explicit matching-edge-rich snake in G(ab, σ_{a,b}),
  edge-by-edge (the matching-side remark in the open-problems section).
- `k3.cpp` — the three-color census: deficient-set extraction and coset-triangle
  reduction behind the k-color hierarchy item (compile `g++ -O2`).
- `cpsat_rho3.py` — independent CP-SAT certification of the twelve n = 11 three-color
  witnesses (arc-MTZ three-color model). Requires `ortools`.
- `fmu.cpp` — exact bitmask dynamic program for max_π μ(π) over G(n,σ); validated
  against direct search on all of S_n, n ≤ 7, and against the conversion inequality
  pointwise (compile `g++ -O2`).
- `fmu_cpsat.py` — integer-programming formulation for max_π μ on structured
  instances, optimality certified. Requires `ortools`.
- `trap_snake.py` — edge-by-edge verifier for the chain snake (matching content of
  the trap-chain family, sizes to n = 5201).
- `cpp/rho_tool.cpp` — independently written C++ implementation; source of the
  exhaustive n ≤ 12 statistics (`exh*.txt`) and of cross-checks
  of the exact values at (a,b) ∈ {(6,8),(5,10),(4,12)} and the trap-family values.

## Algorithm note

The exact search explores the state graph on triples (visited set, current vertex,
last color), each state visited once; no pruning. Exactness is by construction; cost is
the number of reachable states, which the logs report per instance (57,066 at n = 80;
5,012,186 at n = 600; 7,273,586 at n = 720 for σ_{a,b}). On unstructured permutations
the state space is exponential and the method was used only for n ≤ 19.
