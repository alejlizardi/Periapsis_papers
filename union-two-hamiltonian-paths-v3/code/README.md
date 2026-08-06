# Reproducibility artifact

Ancillary files for: *How long a properly colored path must the union of two
Hamiltonian paths contain?* (Lizardi, July 2026).

These files accompany the paper as its ancillary artifact; in the source repository
they live in `code/`. The paper cites individual logs by name.

Environment used to produce `logs/`: Python 3.14, Windows 11, single core; no external
Python dependencies for the core verifiers (standard library only). Two optional
scripts, `cpsat_rho3.py` and `fmu_cpsat.py`, use Google OR-Tools (`pip install
ortools`) to certify optimality of the three-color witnesses and the exact max-μ
values by an independent constraint-programming model. The C++ tools compile with
`g++ -O2` (`rho_tool.cpp` additionally uses `-fopenmp`).

## Claim-to-log map (paper section "Computational verification")

The item numbers match the numbered list of verified statements in that section.

| paper item | claim | script | log |
|---|---|---|---|
| 1 | min over S_b of M(π) = 3 for 4 ≤ b ≤ 9 (full enumeration) | `search_small.py` | `search_small.log` |
| 2 | exact ρ(σ_{a,b}) values up to n = 720 (with reachable-state counts) | `endtoend.py`, `bigtable.py` | `endtoend.log`, `bigtable.log` |
| 3 | exhaustive censuses: exact ρ(σ) over all of S_n, 7 ≤ n ≤ 12 (C++; 479,001,600 permutations at n = 12) — source of the ρ_min census appendix and the deficient-set counts \|D(n)\|; anchors ρ_min(7)=5, ρ_min(8)=6 independently reproduced by the Python solvers | `cpp/rho_tool.cpp`; anchors: `pc.py`, `search_small.py` | `exh7.txt` … `exh12.txt`; `search_small.log` |
| 4 | ρ_min(13..16) + all minimizers (threshold censuses, dual-enumerator at n = 13, 14) and \|D(13..16)\| (orbit-weighted counting censuses); every minimizer certified two ways, every count certified two-sided | `census1316/` (see `census1316/RESULTS.md`) | `census1316/cert_rho_n13..16.log`, `census1316/cert_D_n13..16.log` |
| 5 | odd gadget T′_b (Appendix "The odd gadget"): inventory/highway/crossing claims for odd b ∈ {7,…,21}, M(T′_b)=3, with positive and negative controls; transit spectrum for all π ∈ S_b, b ≤ 8 (μ(2)=2, μ(b)=3 for 3 ≤ b ≤ 8, value 4 unattained at odd b) | `verify_oddb.py`, `verify_spectrum.py` | `verify_oddb.log`, `verify_spectrum.log` |
| open-problems §(matching side) | G(ab, σ_{a,b}) has a path with ≥ a(b−3) = n − 3a matching edges (explicit snake, edge-by-edge check at 6 points up to n = 600) | `verify_f_linear.py` | `verify_f_linear.log` |
| open-problems §(more colors) | ρ_min^(3)(n) = n for 4 ≤ n ≤ 10; ρ_min^(3)(11) = 10, unique deficient triple up to symmetry (coset-triangle census over the deficient set D(n)); all 12 witnesses independently certified by CP-SAT | `k3.cpp` (census; compile `g++ -O2`); `cpsat_rho3.py` (needs `ortools`) | (rerun to reproduce; census outputs archived in the research repo) |
| open-problems §(matching side) | exact f(n): f = n for n ≤ 6, f(7) = 6, f(8) = 7, f(9) = 8; max_π μ ∈ {n−1, n} for all σ, n ≤ 9, with {max μ = n−1} = D(n); exact max_π μ on σ_{a,6} (a = 2..6) and σ_{3,8} (one edge lost per gadget copy beyond the second) | `fmu.cpp` (bitmask DP; validated vs direct search on all S_n, n ≤ 7, and vs the conversion inequality pointwise); `fmu_cpsat.py` (CP-SAT, optimality certified; needs `ortools`) | (rerun to reproduce) |

## How to reproduce

```
python verify_M.py        # ~1 s   independent verifier: classification + M values
python validate2.py       # ~1 s   solver cross-validation + odd-b parity + two-visit total
python endtoend.py        # ~1 s   exact rho of inflations vs bound, n ≤ 80
python search_small.py    # ~40 s  rho_min(7,8) + exhaustive min M over S_b, b ≤ 9
python bigtable.py        # ~40 s  exact rho table up to n = 720
python verify_oddb.py     # ~1 s   odd-gadget structural confirmation + controls
python verify_spectrum.py # ~4 min transit spectrum, all of S_b for b ≤ 8
python cpsat_rho3.py      # optional (ortools): certify the 12 three-color witnesses
python fmu_cpsat.py       # optional (ortools): certify exact max-mu values
```

The three-color census (`k3.cpp`) and the max-μ bitmask DP (`fmu.cpp`) compile with
`g++ -O2`; their outputs back the k-color and matching-content items of the
open-problems section.

## The n = 13–16 censuses (`census1316/`)

Self-contained directory backing paper item 4: engines (`census.cpp`,
`census_gpu.py`, `count_gpu.py`), witness lists, and certificates for
ρ_min(13..16), all minimizers, and |D(13..16)|. See `census1316/RESULTS.md` for
the values, the file map, and how each certificate re-verifies. Re-running the
censuses themselves needs a CUDA GPU and hours-to-days; re-verifying the
certificates takes minutes on CPU (`ortools` required), and the exact-once
coverage check (`burnside_canon.py`) runs in a second with no GPU at all.

`MANIFEST.sha256` lists checksums of all code and logs.

## Code inventory

- `pc.py` — reference solvers: ρ(σ), constrained PC-path search, transit number M
  (plain recursive).
- `family.py` — memoized iterative (state-dedup) solvers; the gadget family T_b.
- `inflate.py` — the inflation σ = τ[π_1..π_a]; the two-visit total (`M2`).
- `endtoend.py` — exact longest-PC-path search with state dedup (used at scale).
- `verify_M.py` — **minimal independent verifier** (~100 lines, no shared code):
  brute-force enumeration of all admissible visits, directly from the definition.
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
- `cpp/rho_tool.cpp` — independently written C++ implementation; source of the
  exhaustive n ≤ 12 statistics (`exh*.txt`) and of cross-checks
  of the exact values at (a,b) ∈ {(6,8),(5,10),(4,12)}.

## Algorithm note

The exact search explores the state graph on triples (visited set, current vertex,
last color), each state visited once; no pruning. Exactness is by construction; cost is
the number of reachable states, which the logs report per instance (57,066 at n = 80;
5,012,186 at n = 600; 7,273,586 at n = 720 for σ_{a,b}). On unstructured permutations
the state space is exponential and the method was used only for n ≤ 19.
