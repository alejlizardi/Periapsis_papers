# Properly colored paths in the union of two Hamiltonian paths

A single self-contained paper on the longest *properly colored* path forced in the
union of two Hamiltonian paths on a common vertex set. Writing $\rho(\sigma)$ for the
longest properly colored path in the two-edge-colored multigraph $H_\sigma$ (blue =
natural order, red = value order of a permutation $\sigma$), and
$\rho_{\min}(n) = \min_\sigma \rho(\sigma)$, the paper:

1. introduces $H_\sigma$, the **Conversion Lemma** linking $\rho$ to the matching-edge
   parameter behind the Norin–Steiner–Thomassé–Wollan bound for Lovász's problem, and
   the Lovász motivation;
2. states a natural **half-bound** $\rho(\sigma) \ge \lceil (n+3)/2 \rceil$ (true on a
   block-reversal family and through $n \le 12$) and **refutes** it with an explicit
   $\sigma^*$ at $n = 18$ where $\rho = 10 < 11$;
3. proves the much stronger
   $$\rho_{\min}(n) = O(\sqrt n), \qquad \rho_{\min}(n) \le 8\sqrt n + 20 \ (n \ge 100),$$
   via an **inflation calculus**, the **transit number** $M(\pi)$, and a complete
   classification of the admissible visits of a gadget family ($M(T_b)=3$ for even $b$,
   from a parity obstruction);
4. proves the **transit floor**: $\mu(b) = \min_{\pi \in S_b} M(\pi) = 3$ for every
   $b \ge 3$ — the lower bound by an elementary matching-union argument, the upper
   bound by a second, odd-size gadget family $T'_b$ whose admissible visits are also
   classified completely (appendix);
5. proves that **block statistics give no lower bound**: block-free permutations
   ($b(\sigma) = n$) built as chains of a modular staircase pattern have
   $\rho = \Theta(\sqrt n)$, so $\inf_\sigma \rho(\sigma)/b(\sigma) = 0$; the
   inequality $\rho \ge b(\sigma)$, exhaustively true through $n = 12$, first fails at
   $n = 14$ (exact census: 16 block-free failures); the same chains give
   $\rho_{\min}(N) \le 4\sqrt2\,\sqrt N - 8$ on an infinite set of $N$;
6. gives a second, independent **trap-chain** construction for the main theorem
   (appendix), and poses the paper's **main open problem**: does
   $\rho_{\min}(n) \to \infty$? Any lower bound of order $\sqrt n$ would be best
   possible, and no lower bound growing with $n$ is known.

This version extends the June-2026 draft (`union-two-hamiltonian-paths/`) with items
4 and 5, obtained in a follow-up research campaign on that draft's own open
questions; it supersedes that folder, which is kept unchanged for reference.

> ⚠️ **Status: unrefereed preprint, AI-assisted. Read before relying on this.**
> The proofs and computations were AI-assisted by Claude under the Periapsis research
> structure; no journal referee has yet checked the paper, and the new results (items
> 4–5) have not yet been through an external novelty/prior-art check (owed before
> submission). It is written to be *checkable without trusting the author or the
> assistant*: the central case analyses are complete, locally verifiable
> classifications and elementary matching arguments; every computational claim is
> reproducible from the artifact in `code/` (including minimal independent verifiers
> written directly from the definitions, and CP-SAT values certified optimal by a
> solver unrelated to the reference search); and the main theorem has a second,
> independent proof in an appendix. The proofs are elementary and self-contained; the
> computations are exhaustive or certified within their stated ranges and are
> independent checks, not ingredients of the proofs (except where a result is
> explicitly computational: the $n = 18$ and $n = 14$ counterexamples and the census
> counts).

## Reproduce the computational claims

All computational claims are reproducible from the self-contained artifact in `code/`
(Python standard library only; optional: `ortools` for the CP-SAT certifier, a C++
compiler for the census/enumeration tools). See `code/README.md` for the full
claim-to-log map and driver list, and `code/MANIFEST.sha256` for file hashes. The
load-bearing checks:

```bash
cd code
python verify_M.py        # minimal independent verifier (<100 lines, from the definitions)
python verify_path.py     # edge-by-edge check of the explicit construction
python verify_floor.py    # the transit-floor witness, validated from the definitions
python verify_n14.py      # the n=14 block-free threshold witness
python verify_stair.py    # the staircase chains: structure + exact small values
python endtoend.py        # end-to-end reproduction of the main bound
```

The exhaustive-search and table logs used in the paper are in `code/logs/`. On arXiv
the `code/` directory ships as the standard ancillary directory `anc/`.

## Contents

| Path | What it is |
|------|------------|
| `union_two_hamiltonian_paths_v2.tex` | The paper (LaTeX source; compiles with pdfLaTeX or Tectonic). |
| `union_two_hamiltonian_paths_v2.pdf` | The rendered paper (35 pp). |
| `assets/` | Periapsis logos: title-page lockup (`lockup.png`; `glyph_square.png`/`mark.png` retained but unused in this version). |
| `code/` | The reproducibility artifact: Python verifiers + C++ tools + logs + manifest. |
| `requirements.txt` | Python dependencies for the verifiers (stdlib only; `ortools` optional). |

## Building the PDF

Compiles with **pdfLaTeX** (or Tectonic) and standard packages (`amsmath`, `amsthm`,
`booktabs`, `graphicx`, `fancyhdr`, `hyperref`). Two passes resolve cross-references:

```bash
pdflatex union_two_hamiltonian_paths_v2.tex && pdflatex union_two_hamiltonian_paths_v2.tex
```

## License

- **Paper text and figures** (`.tex`/`.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT.
- The Periapsis logos in `assets/` are trademarks of Periapsis and are not covered by
  the above licenses; they may not be reused to imply endorsement.

## Citation

> A. Lizardi. *Properly colored paths in the union of two Hamiltonian paths: the
> O(√n) collapse, the transit floor, and the failure of block lower bounds.* Working
> paper, Periapsis, 2026.
