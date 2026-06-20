# Properly colored paths in the union of two Hamiltonian paths can have length O(√n)

A proof that the shortest guaranteed properly colored path in the union of two
Hamiltonian paths on a common vertex set is only $O(\sqrt n)$ — so its ratio to $n$
tends to $0$. Writing $\rho(\sigma)$ for the longest properly colored path in the
two-edge-colored multigraph $H_\sigma$ (blue = natural order, red = value order of a
permutation $\sigma$), the paper proves
$$\rho_{\min}(n) = \min_\sigma \rho(\sigma) = O(\sqrt n), \qquad \rho_{\min}(n) \le 8\sqrt n + 20 \ (n \ge 100).$$
This answers, in the negative, the main question of the author's earlier note: there is
**no** constant $c>0$ forcing a properly colored path on $cn$ vertices. It closes that
particular route toward Lovász's problem on Hamiltonian paths in vertex-transitive
graphs (not the underlying matching-edge question — see the paper's scope remark).

The engine of the proof is an **inflation/substitution** operation on permutations with
a per-permutation invariant, the **transit number** $M(\pi)$, and a complete
classification of the admissible terminal-to-terminal visits of an explicit gadget
family (the five admissible visits of $T_b$, with $M(T_b)=3$ independent of $b$). The
classification turns on a parity obstruction in a long run of doubled edges.

> ⚠️ **Status: unrefereed preprint, AI-assisted. Read before relying on this.**
> The proofs and computations were AI-assisted by Claude under the Periapsis research
> structure; no journal referee has yet checked the paper. It is written to be
> *checkable without trusting the author or the assistant*: the central case analysis is
> a complete, locally verifiable classification; every computational claim is reproducible
> from the artifact in `code/` (including a minimal independent verifier written directly
> from the definitions, under 100 lines); and the main theorem has a second, independent
> proof in an appendix. A hostile adversarial review was run prior to release.

## Reproduce the computational claims

All computational claims are reproducible from the self-contained artifact in `code/`
(Python, with an optional parallel C++ solver). See `code/README.md` for the full driver
list and `code/MANIFEST.sha256` for file hashes. The load-bearing checks:

```bash
cd code
python verify_M.py        # the minimal independent verifier (<100 lines, from the definitions)
python verify_path.py     # edge-by-edge check of the construction
python endtoend.py        # end-to-end reproduction of the bound
```

The exhaustive-search and table logs used in the paper are in `code/logs/`.

## Contents

| Path | What it is |
|------|------------|
| `cstar_paper.tex` | The paper (LaTeX source; compile with pdfLaTeX). |
| `cstar_paper.pdf` | The rendered paper (19 pp). |
| `assets/` | Periapsis logos (title-page lockup; footer mark; bottom-left corner glyph). |
| `code/` | The reproducibility artifact: Python verifiers + C++ solver + logs + manifest. |

## Building the PDF

Compiles with **pdfLaTeX** and standard packages (`amsmath`, `amsthm`, `graphicx`,
`fancyhdr`, `eso-pic`, `hyperref`). Two passes resolve cross-references and the
bibliography:

```bash
pdflatex cstar_paper.tex && pdflatex cstar_paper.tex
```

## License

- **Paper text and figures** (`cstar_paper.tex`, `cstar_paper.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT (see `code/README.md` if it states otherwise).
- The Periapsis logos in `assets/` are trademarks of Periapsis and are not covered by the
  above licenses; they may not be reused to imply endorsement.

## Citation

> A. Lizardi. *Properly colored paths in the union of two Hamiltonian paths can have
> length O(√n).* Working paper, Periapsis, 2026.
> https://github.com/alejlizardi/Periapsis_papers
