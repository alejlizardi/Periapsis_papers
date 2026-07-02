# Periapsis papers

Research notes and papers from **Periapsis**, published openly with the code needed to
reproduce and verify each result.

Each subdirectory is a self-contained paper: LaTeX source, the rendered PDF, and a standalone
verification script a reviewer can run.

## Papers

| Directory | Result | Verify |
|-----------|--------|--------|
| [`hexagonal-hardin-family/`](hexagonal-hardin-family/) | Unified order-polynomial proof for the side-2/3/4 hexagonal Hardin arrays: OEIS [A216938](https://oeis.org/A216938), [A216939](https://oeis.org/A216939), and [A216940](https://oeis.org/A216940); includes Kauers-Koutschan Conjecture 12 as the side-4 case. Supersedes the side-4-only draft in `archive/`. | `python hexagonal-hardin-family/code/verify_family.py` |
| [`a269021-kk-conjecture-17/`](a269021-kk-conjecture-17/) | Proof of the conjectured order-4/degree-21 recurrence for OEIS [A269021](https://oeis.org/A269021) (permutations of [2n] with an increasing subsequence of length n) — Kauers-Koutschan Conjecture 17 — via the Borodin-Okounkov-Olshanski discrete Bessel kernel and a Wilf-Zeilberger certificate. Unrefereed; not yet read by any human expert. | `python a269021-kk-conjecture-17/code/verify_a269021.py` |

## Archived

Earlier standalone drafts that were later incorporated into the more comprehensive papers above.
Kept for reference; not maintained.

| Directory | Result | Superseded by |
|-----------|--------|---------------|
| [`archive/a216940-hardin-hexagon/`](archive/a216940-hardin-hexagon/) | Earlier side-4-only draft for OEIS [A216940](https://oeis.org/A216940). | `hexagonal-hardin-family/` |
| [`archive/cstar-properly-colored-paths/`](archive/cstar-properly-colored-paths/) | Properly colored paths in the union of two Hamiltonian paths can have length O(sqrt(n)); rho_min(n)/n -> 0, answering the author's earlier question in the negative. | `union-two-hamiltonian-paths/` |

## Philosophy

These are released as **working drafts under review** unless stated otherwise. Where a result
is computer-assisted, the repository ships the verification code so the claim can be checked
independently rather than taken on trust. Caveats and the review status are stated in each
paper's own README.

## License

Per-paper. Each directory states its own license; in general, paper text is CC BY 4.0 and
code is MIT. The Periapsis logos are trademarks and are excluded from those licenses.
