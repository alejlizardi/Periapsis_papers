# A proof of the conjectured closed form for OEIS A181280

A single self-contained paper proving the Kauers–Koutschan **Conjecture 20** closed form for
the integer sequence **OEIS A181280**. The sequence counts the $4\times n$ matrices $M$ over
$\mathbb{F}_2$ whose four rows (as bit strings) are lexicographically strictly increasing and
whose Gram matrix $G = MM^{\top} \bmod 2$ has its four rows lexicographically strictly
decreasing. Kauers and Koutschan (arXiv:2303.02793) guessed a closed form for $a(n)$ and left
it open (status "O"); the OEIS still marks it "(conjectured)". This paper proves it.

The proof:

1. both defining conditions are decided by reading $M$ **column by column**, so $a(n)$ is the
   number of accepting length-$n$ walks of an explicit **1242-state transfer matrix** $T$, with
   $a(n) = \mathbf u\,T^n\,\mathbf e_0$;
2. the **Krylov dimension** of $T$ at $\mathbf e_0$ is $16$ — an **a-priori** bound (independent
   of the conjectured formula) on the order of a constant-coefficient recurrence for $a(n)$;
3. within that bound, a nonsingular order-11 Hankel matrix (lower bound) together with a
   **residual-propagation** argument — the order-≤16 a-priori recurrence annihilates the residual
   of the candidate, so finitely many zero residuals force all of them (upper bound) — pins the
   minimal recurrence to order **11**, characteristic polynomial $(x-16)(x-8)^2(x+8)(x-4)^3(x+4)^2(x-2)^2$;
4. the conjectured closed form is the unique solution of that recurrence with the matching
   initial data (Stanley, *Enumerative Combinatorics* I, Thm 4.1.1), so it equals $a(n)$ for
   all $n \ge 4$.

The argument is a **finite, terminating certificate**, not a match of finitely many catalog
terms — the order bound comes from the problem's structure, not from the formula. The method
(read the array into an automaton; observe the count is C-finite for fixed row-count; solve the
recurrence) is **entirely standard**; the paper claims no methodological novelty, only the
specific result, which it could not find in the literature.

> ⚠️ **Status: unrefereed preprint, AI-assisted. Read before relying on this.**
> The proof and computations were AI-assisted by Claude under the Periapsis research structure;
> no journal referee has yet checked the paper. It is written to be *checkable without trusting
> the author or the assistant*: the order bound is derived from Cayley–Hamilton, the
> uniqueness is a one-line induction, the recurrence/closed-form equivalence is Stanley's
> Theorem 4.1.1 (pinned, with the convention quoted), and every computational claim is
> reproducible from the artifact in `code/` — including an independent recount by a
> structurally different (Walsh–Hadamard/Fourier) machine and two independent brute forces. A
> dual adversarial review (a hostile math referee and a format referee) was run prior to
> release; the math referee's verdict was "attack failed." The single step a human should still
> read is the modeling of $a(n)$ by the automaton; it is corroborated but not human-audited.

## Reproduce the computational claims

All claims are reproducible from the self-contained artifact in `code/` (Python standard
library + SymPy). See `code/README.md` for the full claim-to-driver map and the reproduce
order, and `code/MANIFEST.sha256` for file hashes. The load-bearing chain:

```bash
cd code
python -m pip install -r ../requirements.txt   # installs sympy
python transfer.py             # 1242-state count vs OEIS, n<=12
python proof_matrix.py         # Krylov dim 16; minimal order 11; char poly; annihilation
python machine_fourier_coupled.py   # independent 8-state Fourier recount (writes a_terms_fourier.json)
python derive_closed_form.py        # independent Berlekamp-Massey + symbolic C20 diff = 0
python recurrence.py           # a(n) == Conjecture 20 for n = 4..60
```

On arXiv the `code/` directory ships as the standard ancillary directory `anc/`.

## Contents

| Path | What it is |
|------|------------|
| `a181280_binary_matrix_gram.tex` | The paper (LaTeX source; compile with pdfLaTeX). |
| `a181280_binary_matrix_gram.pdf` | The rendered paper (7 pp). |
| `assets/` | Periapsis logos: title-page lockup and bottom-left corner glyph. |
| `code/` | The reproducibility artifact: exact-arithmetic verifiers + an independent cross-derivation + logs map. |
| `requirements.txt` | Python dependencies (SymPy). |

## Building the PDF

Compiles with **pdfLaTeX** and standard packages (`amsmath`, `amsthm`, `booktabs`,
`graphicx`, `fancyhdr`, `eso-pic`, `hyperref`). Two passes resolve cross-references:

```bash
pdflatex a181280_binary_matrix_gram.tex && pdflatex a181280_binary_matrix_gram.tex
```

## License

- **Paper text and figures** (`.tex`/`.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT.
- The Periapsis logos in `assets/` are trademarks of Periapsis and are not covered by the
  above licenses; they may not be reused to imply endorsement.

## Citation

> A. Lizardi. *A proof of the conjectured closed form for OEIS A181280.* Working paper,
> Periapsis, 2026.
