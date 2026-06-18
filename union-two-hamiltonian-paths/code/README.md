# Reproducibility artifact

Ancillary files for: *Properly colored paths in the union of two Hamiltonian paths:
from a refuted half-bound to an O(√n) collapse* (Lizardi, June 2026).

On arXiv these files ship in the standard ancillary directory `anc/`; in the source
repository they live in `code/`. The two are the same files — the paper refers to them
as the ancillary artifact and cites individual logs by name.

Environment used to produce `logs/`: Python 3.14, Windows 11, single core; no external
Python dependencies (standard library only). Total wall-clock to regenerate all Python
logs: about 5 minutes. The C++ tool (`cpp/rho_tool.cpp`) compiles with
`g++ -O2 -fopenmp`; it was used for the exhaustive S_n enumerations (n ≤ 12) and for
independent cross-checks.

## Claim-to-log map (paper Section 9, "Computational verification")

The item numbers below match the numbered list of verified statements in Section 9 of
the paper.

| paper item | claim | script | log |
|---|---|---|---|
| 9.1 | anchors of §4: ρ(σ_m)=2m+3 (m=1..4); ρ(σ*)=10 at n=18; ρ_min(7)=5, ρ_min(8)=6 | `pc.py`, `search_small.py` | `validate2.log`, `search_small.log` |
| 9.2 | classification of admissible visits of H_{T_b} is exactly (V1)–(V5), even b ≤ 14; M(T_b)=b for odd b ≤ 11 | `verify_M.py` (minimal independent verifier, written directly from the admissible-visit/transit-number definition, no shared code) | `verify_M.log` |
| 9.3 | M(T_b)=3 for every even 6 ≤ b ≤ 100; M(T_b)=b for odd b ∈ {7,9,11,13} | `check_all_even.py`, `validate2.py` | `check_all_even.log`, `validate2.log` |
| 9.4 | min over S_b of M(π) = 3 for 4 ≤ b ≤ 9 (full enumeration) | `search_small.py` | `search_small.log` |
| 9.5 | explicit lower-bound path (Prop. "Explicit long path") valid at 7 points incl. (20,20) | `verify_path.py` | `verify_path.log` |
| 9.6 | exact ρ(σ_{a,b}) values up to n = 720 (with reachable-state counts) | `endtoend.py`, `bigtable.py` | `endtoend.log`, `bigtable.log` |
| 9.7 | joint (b(σ), ρ(σ)) distribution over all of S_n, 7 ≤ n ≤ 12 (C++; 479,001,600 permutations at n = 12) | `cpp/rho_tool.cpp` | `exh7.txt` … `exh12.txt` |
| 9.8 | landscape probe (exploratory) | `probe.py` | (exploratory; rerun to reproduce) |
| App. A | trap-chain exact values ρ(σ_C) = 2k+2(L−1) (cross-check at 6 points) | `check_trap.py` | `check_trap.log` |
| §10 (Lovász side) | G(ab, σ_{a,b}) has a path with ≥ a(b−3) = n − 3a matching edges (explicit snake, edge-by-edge check at 6 points up to n = 600) | `verify_f_linear.py` | `verify_f_linear.log` |

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
python probe.py           # ~5 min  landscape probe (exploratory)
```

`MANIFEST.sha256` lists checksums of all code and logs.

## Code inventory

- `pc.py` — reference solvers: ρ(σ) and constrained PC-path search (plain recursive).
- `family.py` — memoized iterative (state-dedup) solvers; the gadget family T_b.
- `inflate.py` — the inflation σ = τ[π_1..π_a]; the two-visit total (`M2`).
- `endtoend.py` — exact longest-PC-path search with state dedup (used at scale).
- `verify_M.py` — **minimal independent verifier** (~100 lines, no shared code):
  brute-force enumeration of all admissible visits, directly from the definition.
- `verify_path.py` — mechanical check of the explicit lower-bound path
  (Proposition "Explicit long path: the family is Θ(a+b)").
- `verify_f_linear.py` — explicit matching-edge-rich snake in G(ab, σ_{a,b}),
  edge-by-edge (the Lovász-side remark in §9).
- `cpp/rho_tool.cpp` — independent C++ implementation from a parallel verification
  effort; source of the exhaustive n ≤ 12 statistics (`exh*.txt`) and of cross-checks
  of the exact values at (a,b) ∈ {(6,8),(5,10),(4,12)} and the trap-family values.

## Algorithm note

The exact search explores the state graph on triples (visited set, current vertex,
last color), each state visited once; no pruning. Exactness is by construction; cost is
the number of reachable states, which the logs report per instance (57,066 at n = 80;
5,012,186 at n = 600; 7,273,586 at n = 720 for σ_{a,b}). On unstructured permutations
the state space is exponential and the method was used only for n ≤ 19.
