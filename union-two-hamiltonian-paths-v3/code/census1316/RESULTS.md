# The n = 13–16 censuses

Two exhaustive computations per n, sharing one symmetry framework (the eight
ρ-preserving symmetries of H_σ: reversal, value-complementation, inversion):

1. **Threshold census** (ρ_min and all minimizers): enumerate canonical
   representatives, record every σ with ρ(σ) ≤ K. Engines: `census.cpp` (CPU,
   sharded) and `census_gpu.py` (GPU, Lehmer-rank cursor); full dual-enumerator
   agreement at n = 13, 14.
2. **Counting census** (|D(n)| = #{σ : ρ(σ) < n}): one deficiency decision per
   canonical representative, accumulated with exact orbit weights (never reps × 8).
   Engine: `count_gpu.py` (GPU; `census.cpp` count modes for the CPU cross-checks);
   dual counts at n = 13, 14; exact reproduction of the exhaustive |D(7..12)|.

## Results

| n  | ρ_min | minimizers (canonical) | \|D(n)\|          | \|D(n)\|/n! |
|----|-------|------------------------|-------------------|-------------|
| 13 | 9     | 290 (62)               | 26,310,368        | 0.4225%     |
| 14 | 10    | 2,894 (476)            | 281,638,656       | 0.3231%     |
| 15 | 9     | 8 (3)                  | 6,454,817,248     | 0.4936%     |
| 16 | 10    | 240 (47)               | 81,620,723,456    | 0.3901%     |

Minimizer counts are orbit-weighted (canonical representatives in parentheses).
The n = 16 counting census completed 2026-08-06 (546,151 s of GPU, 9 resume events,
cursor closing at exactly 16! = 20,922,789,888,000).

## Certificates

- `cert_rho_n13..16.log` — **every** minimizer certified two ways: `pc.py`
  (exhaustive DFS reference, no shared code with the engines) and `cpsat_rho.py`
  (CP-SAT arc-MTZ integer programming, OPTIMAL required). 62/62, 476/476, 3/3,
  47/47 ALL PASS. Inputs: `witnesses_n13..16.txt` (SURV lines: `SURV rho b orbit
  [perm]`). Harness: `cert_witnesses.py --file <witnesses>`.
- `count_n13..16.txt` — the counting runs' RESULT records (the |D(n)| counters,
  with canonical-representative and deep-split auxiliaries), emitted by
  `count_gpu.py`.
- `cert_D_n13..16.log` — the counting censuses certified two-sided on random
  samples from the runs' deterministic reservoirs: claimed-deficient σ
  confirmed by `pc.py` AND CP-SAT (OPTIMAL, ρ < n); claimed-non-deficient σ
  confirmed by explicit properly colored Hamiltonian paths. 100 + 100 per n,
  ALL PASS. Harness: `cert_twosided.py`.
- `coverage_certificate.log` + `burnside_canon.py` — an enumeration-free check
  that each census visited every rank **exactly once**. Each run reports
  `canon`, the number of canonical representatives it met while scanning all of
  [0, n!); that must equal the number of orbits of S_n under the eight
  symmetries, which Burnside's lemma gives in closed form. A skipped range would
  make `canon` too small and a double-counted one too large. The reported and
  exact values agree to the digit at every n from 10 to 16 (and the formula is
  brute-force validated for n ≤ 8). At n = 16 both are 2,615,361,578,344.

## Notes

- Compile: `g++ -O2 -fopenmp census.cpp`. The GPU engines need CUDA (numba);
  runtimes range from minutes (n = 13) to days (n = 16, ~190 h on one consumer
  GPU). The certificates re-verify in minutes on CPU (`ortools` required for the
  CP-SAT side) — reproducing the *certificates* does not require re-running the
  censuses.
- `pc.py` here is the same reference solver as in `code/`; it is duplicated so
  this directory is self-contained.
