# arXiv submission checklist — hexagonal Hardin family

Goal: post the preprint to arXiv `math.CO`, get a permanent `arXiv:XXXX.XXXXX` URL, then use
that URL in the OEIS edits (see OEIS_UPDATES.md).

---

## 0. Account + endorsement (first-time-only hurdle)

- [ ] Create an arXiv account at https://arxiv.org/user/register using your **institutional or
      primary email**. If you still have a working `@gatech.edu` address, register with that —
      arXiv auto-endorses submitters from recognized academic domains, which can skip step 0
      entirely. Otherwise register with alejlizardi05@gmail.com and request endorsement (below).
- [ ] arXiv requires **endorsement** for a first submission to a category (here `math.CO`).
      Endorsement = an established `math.CO` author vouches that the submission is appropriate.
- [ ] Natural endorser: **Manuel Kauers** — he read the paper, confirmed it, and encouraged you to
      submit. Christoph Koutschan is an equally natural alternative. Either has published in
      math.CO and can endorse.
- [ ] When you start a submission, arXiv shows your personal **endorsement code** for the category.
      Send that code to the endorser with the request email below. They enter it at
      https://arxiv.org/auth/endorse and you're cleared (one-time, per category).

## 1. Metadata (paste into the arXiv submission form)

- **Title:** Closed forms for the hexagonal Hardin arrays via order polynomials
- **Authors:** Alejandro Lizardi
- **Primary category:** math.CO (Combinatorics)
- **Cross-list (optional):** math.NT (Number Theory) — minor; can omit
- **MSC class:** 05A15 (primary); 06A07, 52B20 (secondary)
- **Comments field:** 8 pages. Reproducibility script included as an ancillary file.
      Settles Conjecture 12 of Kauers and Koutschan (arXiv:2303.02793).
- **Abstract:** copy the abstract verbatim from the paper (the `\begin{abstract}` block in
      hexagonal_hardin_family.tex). arXiv wants plain text — replace the few TeX bits:
        - `\OEIS{A216938}` etc. -> `A216938`
        - `$a_s(n)$` -> `a_s(n)`, `$\{0,1,\dots,n\}$` -> `{0,1,...,n}`, etc.
        - `$\Om_{P_s}(n+1)$` -> `Omega_{P_s}(n+1)`
        - `$\lvert P_s\rvert = 3s^2-3s+1$` -> `|P_s| = 3s^2-3s+1`
        - `---` -> `--`  (em-dashes)

## 2. Files to upload

arXiv compiles from source. Upload a single .zip (or .tar.gz) containing:

- [ ] `hexagonal_hardin_family.tex`           (the manuscript — bibliography is inline
      `thebibliography`, so NO .bib/.bbl file is needed)
- [ ] `assets/lockup.png`                      (title-block logo; referenced)
- [ ] `assets/glyph_square.png`                (page-corner mark; referenced)
- [ ] `code/verify_family.py`  as an **ancillary file** (see step 4)

Do NOT upload: the compiled `hexagonal_hardin_family.pdf` (arXiv rebuilds it), `assets/mark.png`
(unused), `README.md`, `OEIS_UPDATES.md`, this file, `requirements.txt`, or any `.aux/.log/.out`.

A build-it script is provided: run `bash make_arxiv_zip.sh` in this directory to produce
`arxiv_upload.zip` with exactly the right contents.

## 3. Compile-cleanliness pre-checks (do before uploading)

- [ ] Compiles with pdfLaTeX in two passes locally with exit 0 (already verified after the last
      edits). arXiv uses a recent TeXLive; the packages used (amsmath, amssymb, amsthm, mathtools,
      microtype, booktabs, array, graphicx, xcolor, enumitem, eso-pic, listings, hyperref, lmodern)
      are all standard and present on arXiv.
- [ ] No `\input` of external files, no absolute paths — confirmed (assets referenced relatively).
- [ ] hyperref is loaded last among the conflicting packages — confirmed (it's the last \usepackage).

## 4. Ancillary file (the verifier)

arXiv supports "ancillary files" shown alongside the abstract. Put `verify_family.py` in an
`anc/` folder inside the upload (arXiv convention) OR mark it ancillary in the upload UI. This
lets reviewers/readers run the one-command reproducibility check. The build script handles this.

## 5. After it's live

- [ ] Note the assigned identifier `arXiv:XXXX.XXXXX` and abstract URL.
- [ ] Fill those into OEIS_UPDATES.md, then submit the three OEIS edits.
- [ ] (Optional) update the paper's own repro section / README to cite the arXiv ID.
- [ ] (Optional) submit to Journal of Integer Sequences, framed around the uniform family.

---

## Endorsement-request email (to Kauers, optional cc Koutschan)

> Subject: arXiv endorsement for math.CO — hexagonal Hardin / Conjecture 12 paper
>
> Dear Prof. Kauers,
>
> Thank you again for reading the manuscript and for the encouragement to submit it. As this
> will be my first arXiv submission to math.CO, the system requires an endorsement from an
> established author in that category, and I was hoping you might be willing.
>
> My arXiv endorsement code is [CODE], which I can enter at the start of the submission; you can
> endorse at https://arxiv.org/auth/endorse using that code. The paper is the one I shared,
> "Closed forms for the hexagonal Hardin arrays via order polynomials," which settles your
> Conjecture 12 together with the side-2 and side-3 cases.
>
> No obligation at all if you'd rather not — I can also seek endorsement elsewhere. Thank you
> for considering it.
>
> Best regards,
> Alejandro Lizardi

(arXiv generates [CODE] when you begin the submission; paste it in before sending.)
