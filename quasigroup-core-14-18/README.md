# Recursively differentiable quasigroups of orders 14 and 18

A single self-contained paper settling the **last two open cases** of the
Couselo–González–Markov–Nechaev (CGMN) conjecture on recursively differentiable
quasigroups, by explicit construction.

A quasigroup $(\Omega,\cdot)$ is *recursively differentiable* if its first recursive
derivative $x*y := y\cdot(x\cdot y)$ is again a quasigroup (a Latin square). CGMN (1998)
conjectured that one exists for every order $q\notin\{2,6\}$ and verified all orders
except $q\in\{14,18,26,42\}$. The case $q=42$ was settled in 2008 and $q=26$ in 2026,
leaving $q=14$ and $q=18$. This paper:

1. exhibits an explicit recursively differentiable quasigroup of **order 14** and one of
   **order 18** (full Cayley tables and their recursive derivatives are printed in the
   paper), each invariant under a fixed-point-free automorphism with two regular orbits;
2. records the elementary arithmetic reason the standard perfect-Mendelsohn-design
   constructions skip exactly these two orders
   ($14,18\not\equiv 0,1 \bmod 4$ and $\not\equiv 0,1 \bmod 5$), explaining why they
   remained open and motivating the prescribed-automorphism search used here;
3. concludes the CGMN conjecture: a recursively differentiable quasigroup of order $q$
   exists **iff** $q\notin\{2,6\}$.

> ⚠️ **Status: unrefereed preprint, AI-assisted. Read before relying on this.**
> The constructions and computations were AI-assisted by Claude under the Periapsis
> research structure; no journal referee has yet checked the paper. It is written to be
> *checkable without trusting the author or the assistant*: both results are explicit
> finite Cayley tables, printed in full, whose recursive-differentiability is an
> $O(n^2)$ check reproduced by the standalone verifier in `code/`. A dual adversarial
> review (a hostile math referee and a format referee) was run prior to release.
> One open item, stated in the paper's honest-status note: the **novelty** claim (no
> prior explicit order-14/18 example) rests on a literature survey that reached the
> controlling sources but did *not* include a first-hand reading of the original 1998
> CGMN paper (Russian); a reader with access should confirm that before treating the
> result as a priority claim rather than a clean construction.

## Reproduce the computational claims

Self-contained, Python standard library only:

```bash
cd code
python rdiff_verify.py     # < 1 s; certifies both orders, exit code 0
```

The verifier reads the two Cayley tables (the same integers printed in the paper, also
stored as `code/L14.json`, `code/L18.json`), checks each is a Latin square, forms the
recursive derivative two independent ways (closed form and recurrence-unrolling), checks
the derivative is also a Latin square, and runs the wrong-reading control. See
`code/README.md` for the claim-to-check map and `MANIFEST.sha256` for file hashes.

## Contents

| Path | What it is |
|------|------------|
| `quasigroup_core_14_18.tex` | The paper (LaTeX source; compile with pdfLaTeX). |
| `quasigroup_core_14_18.pdf` | The rendered paper (8 pp). |
| `assets/` | Periapsis logos: title-page lockup and bottom-left corner glyph. |
| `code/` | Reproducibility artifact: standalone verifier + the two Cayley tables (JSON). |
| `requirements.txt` | Python dependencies (none — standard library only). |
| `MANIFEST.sha256` | Checksums of the code and data files. |

## Building the PDF

Compiles with **pdfLaTeX** and standard packages (`amsmath`, `amsthm`, `booktabs`,
`array`, `graphicx`, `fancyhdr`, `eso-pic`, `hyperref`). Two passes resolve
cross-references:

```bash
pdflatex quasigroup_core_14_18.tex && pdflatex quasigroup_core_14_18.tex
```

## License

- **Paper text and figures** (`.tex`/`.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT.
- The Periapsis logos in `assets/` are trademarks of Periapsis and are not covered by the
  above; they may not be reused to imply endorsement.

## Citation

> A. Lizardi. *Recursively differentiable quasigroups of orders 14 and 18: the last open
> cases of the Couselo–González–Markov–Nechaev conjecture.* Working paper, Periapsis,
> 2026.
