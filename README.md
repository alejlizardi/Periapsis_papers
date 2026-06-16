# Periapsis papers

Research notes and papers from **Periapsis**, published openly with the code needed to
reproduce and verify each result.

Each subdirectory is a self-contained paper: LaTeX source, the rendered PDF, and a standalone
verification script a reviewer can run.

## Papers

| Directory | Result | Verify |
|-----------|--------|--------|
| [`hexagonal-hardin-family/`](hexagonal-hardin-family/) | Unified order-polynomial proof for the side-2/3/4 hexagonal Hardin arrays: OEIS [A216938](https://oeis.org/A216938), [A216939](https://oeis.org/A216939), and [A216940](https://oeis.org/A216940); includes Kauers-Koutschan Conjecture 12 as the side-4 case. Supersedes the side-4-only draft below. | `python hexagonal-hardin-family/code/verify_family.py` |
| [`a216940-hardin-hexagon/`](a216940-hardin-hexagon/) | Earlier side-4-only draft for OEIS [A216940](https://oeis.org/A216940), now superseded by `hexagonal-hardin-family/`. | `python a216940-hardin-hexagon/code/verify.py` |
| [`cstar-properly-colored-paths/`](cstar-properly-colored-paths/) | Properly colored paths in the union of two Hamiltonian paths can have length O(sqrt(n)); rho_min(n)/n -> 0, answering the author's earlier question in the negative. | `python cstar-properly-colored-paths/code/verify_M.py` |

## Philosophy

These are released as **working drafts under review** unless stated otherwise. Where a result
is computer-assisted, the repository ships the verification code so the claim can be checked
independently rather than taken on trust. Caveats and the review status are stated in each
paper's own README.

## License

Per-paper. Each directory states its own license; in general, paper text is CC BY 4.0 and
code is MIT. The Periapsis logos are trademarks and are excluded from those licenses.
