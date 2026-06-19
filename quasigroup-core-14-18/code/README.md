# Reproducibility artifact

Ancillary files for: *Recursively differentiable quasigroups of orders 14 and 18: the
last open cases of the Couselo–González–Markov–Nechaev conjecture* (Lizardi, June 2026).

On arXiv these files ship in the standard ancillary directory `anc/`; in the source
repository they live in `code/`. The two are the same files.

Environment used: Python 3.14, Windows 11; no external Python dependencies (standard
library only). The verifier runs in well under a second.

## Contents

| File | What it is |
|------|------------|
| `rdiff_verify.py` | The standalone verifier (paper Section 7). Re-checks the tables from the definitions. |
| `L14.json` | The order-14 base Cayley table $L_{14}$ (Table 1 of the paper). |
| `L18.json` | The order-18 base Cayley table $L_{18}$ (Table 3 of the paper). |
| `search/search.cpp` | The prescribed-automorphism search that produced the two tables (paper Section 4). Provided for reproducibility; not needed to verify the results. |
| `LICENSE` | MIT license for the code. |

## What the verifier checks (paper Section 7)

For each order $n \in \{14, 18\}$, `rdiff_verify.py`:

1. checks the base table $L$ is a **Latin square** (a quasigroup): every row and column a
   permutation of $\{0,\dots,n-1\}$;
2. forms the first recursive derivative two ways — the closed form
   $D[a][b] = L[b][L[a][b]]$ and the unrolled recurrence
   $a*_kb=(a*_{k-2}b)\cdot(a*_{k-1}b)$ — and confirms they agree;
3. checks the derivative $D$ is a **Latin square** — i.e. the core is again a quasigroup,
   which is the definition of recursive differentiability;
4. runs the **wrong-reading control**: of the four candidate "core" operations
   $b\cdot(a\cdot b),\ a\cdot(a\cdot b),\ (a\cdot b)\cdot b,\ (a\cdot b)\cdot a$, only the
   correct reading $b\cdot(a\cdot b)$ yields a Latin derivative — so each table is
   certified against the exact definition, not a weaker variant.

The JSON tables are the same integers printed in the paper; one can equally type the
printed tables in by hand and run the same four checks.

## How to reproduce

```bash
python rdiff_verify.py      # < 1 s; prints PASS for both orders, exit code 0
```

Expected output ends with:

```
ALL CHECKS PASSED: recursively differentiable quasigroups of orders 14 and 18 are certified.
```

## Reproducing the search (optional)

The two tables were found by the prescribed-automorphism search in `search/search.cpp`.
It is not needed to verify the results, but to regenerate a witness:

```bash
cd search
g++ -O3 -march=native -std=c++17 search.cpp -o search
./search 14 "cyc:0,1,2,3,4,5,6|7,8,9,10,11,12,13" 600000 1
./search 18 "cyc:0,1,2,3,4,5,6,7,8|9,10,11,12,13,14,15,16,17" 3600000 1
```

Any solution found is a valid witness; the exact table depends on the random seed.

`../MANIFEST.sha256` lists checksums of the code and data files.
