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
| `rdiff_verify.py` | The standalone verifier (Section 7 of the paper). Shares no code with the search that produced the tables; re-checks them from the definitions. |
| `L14.json` | The order-14 base Cayley table $L_{14}$, exactly as printed in Section 5 of the paper. |
| `L18.json` | The order-18 base Cayley table $L_{18}$, exactly as printed in Section 6. |
| `LICENSE` | MIT license for the code. |

## What the verifier checks (paper Section 7)

For each order $n \in \{14, 18\}$, `rdiff_verify.py`:

1. checks the base table $L$ is a **Latin square** (a quasigroup): every row and column a
   permutation of $\{0,\dots,n-1\}$;
2. forms the first recursive derivative **two independent ways** and confirms they agree
   — (a) the closed form $D[a][b] = L[b][L[a][b]]$ (paper eq. (2.4)), and (b) by
   unrolling the CGMN recurrence $a*_kb=(a*_{k-2}b)\cdot(a*_{k-1}b)$ step by step (paper
   eq. (1.1));
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

`../MANIFEST.sha256` lists checksums of the code and data files.
