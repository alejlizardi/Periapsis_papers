# Closed forms for the hexagonal Hardin arrays via order polynomials

Preprint (unrefereed by a journal; confirmed by the conjecture authors M. Kauers and C. Koutschan).
Proves the empirical
closed forms for the side-2/3/4 hexagonal Hardin arrays — **OEIS A216938, A216939, A216940** —
uniformly, by identifying each count with the order polynomial of an explicit poset and applying
Stanley's a-priori degree bound plus an independent descent computation. (A216940 = Kauers–Koutschan
Conjecture 12.)

## Files
- `hexagonal_hardin_family.tex` — the manuscript (pdfLaTeX, two passes).
- `hexagonal_hardin_family.pdf` — compiled PDF (8 pp).
- `code/verify_family.py` — **standalone reproducibility verifier**. One command, Python 3 stdlib
  only, no external data files. Rebuilds each poset P_s (s=1,2,3,4) from first principles, counts
  order-preserving maps three independent ways, computes the order polynomial via the descent
  identity, checks the a-priori degree = |P_s| = 3s²−3s+1 and the leading coefficient e(P_s)/|P_s|!,
  and matches the OEIS empirical closed form for s=2,3,4. Exit 0 iff all pass (~1 s).
- `requirements.txt` — records that no third-party Python packages are required.
- `assets/` — Periapsis logo (title lockup + page-corner mark).

## Reproduce
```
python code/verify_family.py        # expect "FAMILY VERIFICATION: ALL PASS"
```
Compile the paper:
```
pdflatex hexagonal_hardin_family.tex && pdflatex hexagonal_hardin_family.tex
```

## Honest status
- **Computer-assisted; unrefereed by a journal.** Confirmed by the two authors of Conjecture 12,
  M. Kauers and C. Koutschan, who have read the manuscript and agree the argument is correct; not
  yet submitted to or refereed by a journal.
- The two order-polynomial theorems (degree = |P| with leading coeff e(P)/|P|!; the descent
  identity) are **cited from Stanley EC1 §3.15, not re-derived** here. Their numerical consequences
  are confirmed on the s=2,3,4 controls.
- The one load-bearing modeling step (poset P_s faithfully models the array) is defended by an
  acyclicity/identification argument plus an independent from-scratch reconstruction of P_2; it is
  the joint a referee should check.
- **Novelty:** the order-polynomial / order-polytope method is standard (Stanley) and has been
  applied to other OEIS sequences (arXiv:2412.18744). We are not aware of this hexagonal family
  having been treated, but make **no priority claim** — pointers to prior appearances welcome.
- **Author:** Alejandro Lizardi (independent researcher; recently B.S., Georgia Institute of
  Technology). Correspondence: alejlizardi05@gmail.com.

## Per-file licensing
- `code/verify_family.py` — released to the public domain (CC0); reuse freely.
- `hexagonal_hardin_family.tex` / `.pdf` — © the author(s); all rights reserved.
- `assets/*` — Periapsis marks; not for reuse outside this paper.
