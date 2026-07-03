# Permutations of [2n] with an increasing subsequence of length n: a proof of the Kauers–Koutschan recurrence

Unrefereed, computer-assisted preprint.
Proves the conjectured order-4, degree-21 recurrence for **OEIS A269021** — Conjecture 17 of
Kauers–Koutschan (arXiv:2303.02793) — via an exact single-sum representation

    a(n) = (2n)!^2 · Σ_{ρ=n}^{2n} (−1)^{ρ−n} C(2ρ−2, ρ−n) / ((ρ!)^2 (2n−ρ)!)

derived from the Robinson–Schensted correspondence, the Borodin–Okounkov–Olshanski discrete
Bessel kernel (one-point function of the poissonized Plancherel measure), and the observation
that a partition of 2n has at most one row of length ≥ n. Zeilberger's algorithm returns exactly the
Kauers–Koutschan operator as a telescoper for this sum; the certificate identities are proved
deterministically (exact evaluation beyond explicit degree bounds). Also yields (to our
knowledge) the first D-finiteness proof for this coupled diagonal and a closed formula for the upper half of the
LIS distribution on S_{2n}.

## Files
- `a269021_conjecture17.tex` — the manuscript (pdfLaTeX, two passes).
- `a269021_conjecture17.pdf` — compiled PDF (11 pp).
- `code/verify_a269021.py` — **standalone reproducibility verifier**. One command, Python 3
  stdlib only, no external data files, exact arithmetic throughout (no floats). Re-derives a(n)
  by brute force over all (2n)! permutations (n ≤ 3) and by hook lengths (n ≤ 10), checks the
  representation, checks that the conjectured operator annihilates it (n ≤ 52), proves the two
  WZ-certificate identities on integer grids exceeding their structural degree bounds, checks
  telescoping at every integer point for n = 5..40, and checks the equivalence with the paper's
  ã-normalized form. Exit 0 iff all 12 checks pass (~10 s).
- `code/certificate_xhat.json` — the WZ certificate numerator X̂(n,ρ) as machine-readable data
  (identical to the `XHAT` constant embedded in the verifier; the script remains standalone).
- `requirements.txt` — records that no third-party Python packages are required.

SHA-256 of `code/verify_a269021.py` as of this revision:
`2dee0cdc613f6d96beaa10f6c2f9f915f3cd1f42e03c7cf25724adda13538042`
- `assets/` — Periapsis logo (title lockup + page-corner mark).

## Reproduce
```
python code/verify_a269021.py       # expect "A269021 VERIFICATION: ALL 12 CHECKS PASS"
```
Compile the paper:
```
pdflatex a269021_conjecture17.tex && pdflatex a269021_conjecture17.tex
```

## Honest status
- **Computer-assisted; unrefereed.** The certificate was discovered by Maxima's `zeilberger`
  package (discovery only); every step of the proof is verified by the shipped script. The
  research phase included an independent automated audit that re-ran all computations plus
  brute-force cross-checks, and a subsequent cross-model skeptical review (different AI vendor)
  that attempted to break the derivation and the boundary treatment and could not; its revision
  requests are incorporated in this version.
- Two ingredients are **cited, not re-derived**: the Robinson–Schensted correspondence
  (Schensted 1961) and the Borodin–Okounkov–Olshanski one-point formula with its kernel series
  (arXiv:math/9905032v2, Theorem 2 + Proposition 2.6, eq. (2.14); statements quoted in the paper and
  cross-checked numerically in exact arithmetic).
- **Novelty:** the ingredients are classical (all available by 2000) and the representation is a
  natural consequence of BOO; we found no prior statement of it or of the resulting proof
  (OEIS entries A269021/A269042/A047874, arXiv, and the LIS-distribution literature were
  searched), but make **no priority claim** — pointers to prior appearances are welcome.
- **Author:** Alejandro Lizardi (independent researcher; recently B.S., Georgia Institute of
  Technology). Correspondence: alejlizardi05@gmail.com.
- AI-assisted by Claude (Anthropic) under the Periapsis research architecture; the verification
  is designed so that no step rests on the assistant's (or any single tool's) say-so.

## Per-file licensing
- `code/verify_a269021.py` — released to the public domain (CC0); reuse freely.
- `a269021_conjecture17.tex` / `.pdf` — © the author(s); all rights reserved.
- `assets/*` — Periapsis marks; not for reuse outside this paper.
