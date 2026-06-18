# Properly colored paths in the union of two Hamiltonian paths

A single self-contained paper on the longest *properly colored* path forced in the union
of two Hamiltonian paths on a common vertex set. Writing $\rho(\sigma)$ for the longest
properly colored path in the two-edge-colored multigraph $H_\sigma$ (blue = natural
order, red = value order of a permutation $\sigma$), and
$\rho_{\min}(n) = \min_\sigma \rho(\sigma)$, the paper:

1. introduces $H_\sigma$, the **Conversion Lemma** linking $\rho$ to the matching-edge
   parameter behind the Norin–Steiner–Thomassé–Wollan bound for Lovász's problem, and
   the Lovász motivation;
2. states a natural **half-bound** $\rho(\sigma) \ge \lceil (n+3)/2 \rceil$ (true on a
   block-reversal family and through $n \le 14$) and **refutes** it with an explicit
   $\sigma^*$ at $n = 18$ where $\rho = 10 < 11$;
3. proves the much stronger
   $$\rho_{\min}(n) = O(\sqrt n), \qquad \rho_{\min}(n) \le 8\sqrt n + 20 \ (n \ge 100),$$
   via an **inflation calculus**, the **transit number** $M(\pi)$, and a complete
   classification of the admissible visits of a gadget family ($M(T_b)=3$ for even $b$,
   from a parity obstruction);
4. gives a second, independent **trap-chain** construction (appendix);
5. collects the lower-bound discussion and open problems (even $\rho_{\min}(n)\to\infty$
   is open).

So $\rho_{\min}(n)/n \to 0$: there is **no** constant $c>0$ forcing a properly colored
path on $cn$ vertices. This answers, in the negative, the relaxation route toward
Lovász's problem introduced here — and **only** that route, not the underlying
matching-edge question (see the paper's scope remark).

This paper merges and supersedes two earlier working notes by the author (an
extremal-conjecture note and a separate $O(\sqrt n)$ note). It is the single source of
truth for the result; the two notes are not needed and are not cited.

> ⚠️ **Status: unrefereed preprint, AI-assisted. Read before relying on this.**
> The proofs and computations were AI-assisted by Claude under the Periapsis research
> structure; no journal referee has yet checked the paper. It is written to be
> *checkable without trusting the author or the assistant*: the central case analysis is
> a complete, locally verifiable classification; every computational claim is
> reproducible from the artifact in `code/` (including a minimal independent verifier
> written directly from the definitions, about 100 lines); and the main theorem has a
> second, independent proof in an appendix. A dual adversarial review (a hostile math
> referee and a format referee) was run prior to release. The proofs are elementary and
> self-contained; the computations are exhaustive within their stated ranges and are
> independent checks, not ingredients of the proofs.

## Reproduce the computational claims

All computational claims are reproducible from the self-contained artifact in `code/`
(Python standard library only, with an optional parallel C++ solver). See
`code/README.md` for the full claim-to-log map and driver list, and
`code/MANIFEST.sha256` for file hashes. The load-bearing checks:

```bash
cd code
python -m pip install -r ../requirements.txt   # (none required; stdlib only)
python verify_M.py        # the minimal independent verifier (<100 lines, from the definitions)
python verify_path.py     # edge-by-edge check of the explicit construction
python endtoend.py        # end-to-end reproduction of the bound
```

The exhaustive-search and table logs used in the paper are in `code/logs/`. On arXiv the
`code/` directory ships as the standard ancillary directory `anc/`.

## Contents

| Path | What it is |
|------|------------|
| `union_two_hamiltonian_paths.tex` | The paper (LaTeX source; compile with pdfLaTeX). |
| `union_two_hamiltonian_paths.pdf` | The rendered paper (22 pp). |
| `assets/` | Periapsis logos: title-page lockup and bottom-left corner glyph (`mark.png` is retained but no longer used in the footer). |
| `code/` | The reproducibility artifact: Python verifiers + C++ solver + logs + manifest. |
| `requirements.txt` | Python dependencies for the verifiers (empty — standard library only). |

## Building the PDF

Compiles with **pdfLaTeX** and standard packages (`amsmath`, `amsthm`, `booktabs`,
`graphicx`, `fancyhdr`, `eso-pic`, `hyperref`). Two passes resolve cross-references:

```bash
pdflatex union_two_hamiltonian_paths.tex && pdflatex union_two_hamiltonian_paths.tex
```

## License

- **Paper text and figures** (`union_two_hamiltonian_paths.tex`/`.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT.
- The Periapsis logos in `assets/` are trademarks of Periapsis and are not covered by
  the above licenses; they may not be reused to imply endorsement.

## Citation

> A. Lizardi. *Properly colored paths in the union of two Hamiltonian paths: from a
> refuted half-bound to an O(√n) collapse.* Working paper, Periapsis, 2026.
