# OEIS update text (ready to paste)

These are the edits to make on the three OEIS sequence pages once the arXiv preprint is
posted. Replace **`arXiv:XXXX.XXXXX`** below with the assigned arXiv identifier, and
**`https://arxiv.org/abs/XXXX.XXXXX`** with its abstract URL.

How to apply each one:
1. Log in at https://oeis.org (register a free account first if needed).
2. Open the sequence (e.g. https://oeis.org/A216940), click **edit** at the top.
3. Add the LINK and COMMENT lines below into the appropriate boxes. Do **not** delete the
   existing empirical-formula line — add the proof comment next to it (OEIS keeps the formula
   and records that it is now proved).
4. In the edit summary box, write something like: *"Empirical/conjectured closed form now proved;
   added reference and link to the proof."*
5. Submit. Edits enter the moderation queue and are reviewed before going live.

The author name renders as a link if written as `_Alejandro Lizardi_` (OEIS author-link syntax).
Date the comments with the day you submit.

---

## A216940 (side 4 — Kauers–Koutschan Conjecture 12)

**LINK** (add):
```
Alejandro Lizardi, <a href="https://arxiv.org/abs/XXXX.XXXXX">Closed forms for the hexagonal Hardin arrays via order polynomials</a>, arXiv:XXXX.XXXXX, 2026.
```

**COMMENT** (add):
```
The conjectured closed form of Kauers and Koutschan (their Conjecture 12) is proved in the Lizardi 2026 link: a(n) is the order polynomial Omega_{P}(n+1) of the side-4 hexagon poset P (37 cells), so by Stanley's theorem (EC1, Sec. 3.15.2) it is a polynomial of degree exactly 37; matching 38 poset-derived values forces equality. - _Alejandro Lizardi_, <date>
```

---

## A216938 (side 2)

**LINK** (add):
```
Alejandro Lizardi, <a href="https://arxiv.org/abs/XXXX.XXXXX">Closed forms for the hexagonal Hardin arrays via order polynomials</a>, arXiv:XXXX.XXXXX, 2026.
```

**COMMENT** (add):
```
The empirical formula is proved in the Lizardi 2026 link: a(n) is the order polynomial Omega_{P}(n+1) of the side-2 hexagon poset P (7 cells), so by Stanley's theorem (EC1, Sec. 3.15.2) it is a polynomial of degree exactly 7, equal to the closed form below. - _Alejandro Lizardi_, <date>
```

(Optionally, the existing "Empirical:" prefix on the formula line may be removed now that it is
proved; leave the formula text itself unchanged.)

---

## A216939 (side 3)

**LINK** (add):
```
Alejandro Lizardi, <a href="https://arxiv.org/abs/XXXX.XXXXX">Closed forms for the hexagonal Hardin arrays via order polynomials</a>, arXiv:XXXX.XXXXX, 2026.
```

**COMMENT** (add):
```
The empirical formula is proved in the Lizardi 2026 link: a(n) is the order polynomial Omega_{P}(n+1) of the side-3 hexagon poset P (19 cells), so by Stanley's theorem (EC1, Sec. 3.15.2) it is a polynomial of degree exactly 19, equal to the closed form below. - _Alejandro Lizardi_, <date>
```

(Optionally, the existing "Empirical:" prefix on the formula line may be removed now that it is
proved; leave the formula text itself unchanged.)

---

## Notes / conventions checked against the live pages (2026-06-18)

- All three sequences have **offset 1**. The paper's a_s(n) (n >= 0) matches the OEIS a(n) on the
  shared variable n; the OEIS empirical formulas evaluate correctly at n=1 (e.g. A216938 gives
  a(1)=10), so no reindexing note is needed. The "Omega_{P}(n+1)" in the comment refers to the
  order-polynomial argument shift, not to the OEIS offset.
- A216940 already links the Kauers–Koutschan paper (arXiv:2303.02793); the new comment ties the
  proof to their Conjecture 12 explicitly.
- Cell counts 3s^2-3s+1 = 7, 19, 37 for s = 2, 3, 4 equal the polynomial degrees.
- Stanley reference: Enumerative Combinatorics, Vol. 1, 2nd ed., Sec. 3.15.2 (degree = |P|,
  leading coeff e(P)/|P|!) and Thm. 3.15.8 (descent identity).
- Do NOT delete existing content; OEIS moderators prefer additive edits.
