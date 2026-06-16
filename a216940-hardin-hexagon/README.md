# A proof of the Kauers–Koutschan conjecture for the hexagonal Hardin array A216940

A short, computer-assisted proof that OEIS sequence
[**A216940**](https://oeis.org/A216940) — *"Number of side-4 hexagonal 0..n arrays
with values nondecreasing E, SW and SE"* — equals the closed form conjectured by
Kauers and Koutschan ([arXiv:2303.02793](https://arxiv.org/abs/2303.02793),
*J. Integer Seq.* **26** (2023), Conjecture 12). This upgrades the sequence's status
in that paper from **D** (D-finiteness proved, closed form not) to **P** (closed form proved).

> ⚠️ **Status: unrefereed first draft, AI-assisted. Read before relying on this.**
> The argument is computer-assisted and was developed and internally audited with AI
> assistance; no human mathematician has yet refereed it, and it has not been circulated
> to the original authors. The computational content (the 38-point match, reproduced by
> `code/verify.py`) is solid and re-runnable. Two links rest on strong evidence but are
> **not yet human-verified**: (i) that the poset `P` faithfully models the hexagonal
> array — the geometry needs a combinatorialist's eye; and (ii) that the two cited Stanley
> theorems (the degree bound and the descent identity) are applied with exactly the right
> conventions — cited, with consequences checked on small cases, not re-derived.
> A **prior-art / novelty check has not been done**: the method is standard P-partition
> theory, so it is possible this is already known or considered folklore. Treat this as a
> claim under verification, not a settled, peer-reviewed result. See *Caveats* below.

## The idea in one paragraph

The admissible fillings of the side-4 hexagon are exactly the order-preserving maps of an
explicit 37-element poset `P`, so the count `a(n)` is the **order polynomial** `Ω_P(n+1)`.
A classical theorem of Stanley gives, *a priori*, that `a` is a polynomial of degree
**exactly 37** (= the number of cells). One then computes 38 values `a(0..37)` directly from
`P` — independently of the conjectured formula — and matches them against Conjecture 12; two
degree-≤37 polynomials agreeing at 38 points are identical. The subtlety the proof is built
to avoid: the published OEIS b-file has only **37** terms, one short of even interpolating a
degree-37 polynomial, so a term-match against the data is *not* a proof. The certificate rests
on the rigorous degree bound plus the genuinely-new 38th value `a(0)=1`.

## Reproduce the proof

The computational content is checked end-to-end by a single standalone script that needs
only Python and `sympy`, and uses **no external data files**:

```bash
pip install -r requirements.txt
python code/verify.py
```

It prints `ALL CHECKS PASSED` iff:

1. the poset built **from geometry alone** reproduces the independent sibling sequences
   A216938 (side-2) and A216939 (side-3) on all published terms — a control that a wrong
   geometry would fail;
2. `deg a₄ = 37` and leading coefficient `= e(P)/37!` (the a-priori bound);
3. all 38 poset-derived values `a(0..37)` equal Conjecture 12 (and `a(1..37)` reproduce the
   b-file);
4. over-determination: `a(38..45)` also match.

Runtime: a few seconds.

## Contents

| Path | What it is |
|------|------------|
| `a216940.tex` | The paper (LaTeX source; compile with pdfLaTeX, e.g. on Overleaf). |
| `a216940.pdf` | The rendered paper. |
| `assets/mark.png` | The Periapsis mark used in the page footer. |
| `code/verify.py` | Standalone, from-scratch verification of the whole result. |
| `code/certificate.json` | The descent distribution `{A_d}` of `P` and the 38 values `a(0..37)`. |
| `requirements.txt` | `sympy`. |

## Building the PDF

The paper compiles with **pdfLaTeX** and only standard packages (`amsmath`, `amsthm`, `tikz`,
`algorithm`/`algpseudocode`, `fancyhdr`, `hyperref`). On [Overleaf](https://overleaf.com):
create a project, upload `a216940.tex` and the `assets/` folder, and compile. Locally:

```bash
pdflatex a216940.tex && pdflatex a216940.tex
```

(Two passes resolve cross-references.)

## Caveats (read before citing)

- **Computer-assisted.** The mathematical spine is a rigorous degree bound (a cited theorem)
  plus an independent 38-point computation. The two Stanley facts it cites — that an order
  polynomial has degree `|P|` with leading coefficient `e(P)/|P|!`, and the descent identity
  `Ω_P(m) = Σ_d A_d·C(m+|P|−1−d, |P|)` — are standard textbook results (Stanley, *Enumerative
  Combinatorics* Vol. 1, §3.15); we cite rather than re-derive them, and confirm their
  numerical consequences on the side-2/3 controls.
- **AI-assisted; not yet human-refereed.** The proof was developed and checked with AI
  assistance; all results were verified in exact arithmetic by the independent computations in
  `code/` (see the paper's Acknowledgements and Section 6). It is published here for
  transparency and to invite scrutiny prior to submission and contact with the original authors.

## License

- **Paper text and figures** (`a216940.tex`, `a216940.pdf`): CC BY 4.0.
- **Code** (`code/`): MIT.
- The Periapsis mark in `assets/` is a trademark of Periapsis and is not covered by the above
  licenses; it may not be reused to imply endorsement.

## Citation

> A. Lizardi. *A proof of the Kauers–Koutschan conjecture for the hexagonal Hardin array
> A216940.* Working paper, Periapsis, 2026.
> https://github.com/PeriapsisResearch/Periapsis_papers
