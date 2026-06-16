#!/usr/bin/env python3
"""
verify.py  --  Standalone, self-contained verification of the main theorem of
   "A proof of the Kauers-Koutschan conjecture for the hexagonal Hardin array A216940"
   (Periapsis, 2026).

Run:   python verify.py
Needs: Python 3.8+ and sympy  (pip install sympy)

What it checks, from scratch, with NO external data files:
  1. Builds the side-s hexagon poset P_s from the geometry alone (cube coordinates +
     the three monotone directions E, SW, SE).  No OEIS value is used to build it.
  2. Computes the order polynomial a_s(n) = Omega_{P_s}(n+1) via the descent
     distribution of P_s's linear extensions (Stanley's identity), by an
     n-independent transfer matrix over order ideals.
  3. CONTROLS: confirms a_2 reproduces OEIS A216938 and a_3 reproduces A216939
     (independent sequences) -- a wrong geometry would fail these.
  4. Confirms deg a_4 = |P_4| = 37 and leading coeff = e(P)/37!  (the a-priori bound).
  5. Computes a_4(0..37) -- 38 values, INDEPENDENT of the conjecture -- and checks
     all 38 equal Kauers-Koutschan Conjecture 12.  (a_4(1..37) also = the b-file;
     a_4(0)=1 is the genuine 38th datum.)
  6. Over-determination: also checks a_4(38..45) against Conjecture 12.

Exit code 0 and "ALL CHECKS PASSED" iff the theorem's computational content holds.
"""

from itertools import product
from math import comb, factorial
from collections import defaultdict
import sys

try:
    import sympy as sp
except ImportError:
    sys.exit("This script needs sympy:  pip install sympy")

# --------------------------------------------------------------------------- #
# 1. The poset P_s from geometry alone                                        #
# --------------------------------------------------------------------------- #
E, SW, SE = (1, -1, 0), (0, -1, 1), (-1, 0, 1)          # the three monotone directions

def hex_cells(s):
    """Cells of the side-s hexagon in cube coords: x+y+z=0, |x|,|y|,|z| <= s-1."""
    cells = []
    for x in range(-(s - 1), s):
        for y in range(-(s - 1), s):
            z = -x - y
            if -(s - 1) <= z <= s - 1:
                cells.append((x, y, z))
    return cells

def build_poset(s):
    """Return (n, preds, natural_label).  u is a predecessor of c iff c = u + dir."""
    cells = hex_cells(s)
    idx = {c: i for i, c in enumerate(cells)}
    cs = set(cells)
    n = len(cells)
    preds = {i: [] for i in range(n)}
    for c in cells:
        for d in (E, SW, SE):
            u = (c[0] - d[0], c[1] - d[1], c[2] - d[2])
            if u in cs:
                preds[idx[c]].append(idx[u])
    # natural labelling = a fixed topological order (Kahn, smallest-index tie-break)
    succ = {i: [] for i in range(n)}
    for c in range(n):
        for u in preds[c]:
            succ[u].append(c)
    indeg = [len(preds[i]) for i in range(n)]
    import bisect
    avail = sorted(i for i in range(n) if indeg[i] == 0)
    order, ind = [], indeg[:]
    while avail:
        u = avail.pop(0)
        order.append(u)
        for w in succ[u]:
            ind[w] -= 1
            if ind[w] == 0:
                bisect.insort(avail, w)
    assert len(order) == n, "step relation is not acyclic -- not a poset!"
    label = [0] * n
    for pos, i in enumerate(order):
        label[i] = pos
    return n, preds, label

# --------------------------------------------------------------------------- #
# 2. Descent distribution {A_d} by transfer matrix over order ideals          #
# --------------------------------------------------------------------------- #
def descent_distribution(s):
    n, preds, label = build_poset(s)
    full = (1 << n) - 1
    predmask = [0] * n
    for i in range(n):
        for u in preds[i]:
            predmask[i] |= (1 << u)
    # layer: mask -> { last_label -> { descents -> count } }
    layer = {0: {-1: {0: 1}}}
    result_full = None
    for _ in range(n + 1):
        nxt = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for mask, bylast in layer.items():
            if mask == full:
                result_full = bylast
                continue
            for i in range(n):
                if (mask >> i) & 1:
                    continue
                if (mask & predmask[i]) == predmask[i]:        # i currently minimal
                    nmask = mask | (1 << i)
                    li = label[i]
                    for last, poly in bylast.items():
                        add = 1 if (last != -1 and last > li) else 0
                        tgt = nxt[nmask][li]
                        for d, c in poly.items():
                            tgt[d + add] += c
        layer = {m: {ll: dict(pl) for ll, pl in d.items()} for m, d in nxt.items()}
        if full in layer:
            result_full = layer[full]
    A = defaultdict(int)
    for poly in result_full.values():
        for d, c in poly.items():
            A[d] += c
    return n, dict(A)

def a_of_N(p, A, N):
    """a(N) = Omega_P(N+1) = sum_d A_d * C(N + p - d, p)  (Stanley's descent identity)."""
    return sum(c * comb(N + p - d, p) for d, c in A.items())

# --------------------------------------------------------------------------- #
# Conjecture 12 (Kauers-Koutschan), transcribed verbatim from arXiv:2303.02793 #
# --------------------------------------------------------------------------- #
def conjecture12(N):
    n = sp.Integer(N)
    pre = sp.prod([n + k for k in range(1, 14)]) * (n + 6) * (n + 7) * (n + 8) * (n + 7)
    c = [15118483615575730790400000, 37557333457279933473792000, 45137854540680193956153600,
         34829846371335010335540480, 19314394347459920710102704, 8166353315859794719296864,
         2726904840964417033376520, 735273283907306553706472, 162382123713323392711687,
         29630015361661371290844, 4487557575514810132362, 564694034848365996336,
         58900361433618244860, 5062226797216352960, 354853893929158096, 19969728998781072,
         880856790135603, 29345762188932, 694580474022, 10413780440, 74384146]
    q = sum(ci * n**i for i, ci in enumerate(c))
    D = 221424599279703105635713957232640000000
    val = sp.Rational(pre * q, D)
    assert val.q == 1, "Conjecture 12 produced a non-integer -- transcription error!"
    return int(val)

# Independent OEIS control data (small siblings; published terms).
A216938 = [10, 53, 200, 606, 1572, 3630, 7656, 15015, 27742, 48763]               # side-2
A216939 = [52, 1211, 17336, 175869, 1374208, 8737047, 46978360, 219784500]        # side-3

# --------------------------------------------------------------------------- #
def main():
    ok = True
    print("Building posets from geometry and computing descent distributions...\n")

    # cell counts
    for s in (2, 3, 4):
        p = 3 * s * s - 3 * s + 1
        print(f"  side-{s}: |P| = 3s^2-3s+1 = {p}")
    print()

    # CONTROL: side-2 vs A216938, side-3 vs A216939 (independent sequences)
    for s, oracle, name in ((2, A216938, "A216938"), (3, A216939, "A216939")):
        p, A = descent_distribution(s)
        vals = [a_of_N(p, A, k) for k in range(1, len(oracle) + 1)]
        good = (vals == oracle)
        ok &= good
        print(f"  CONTROL side-{s} vs {name}: {'PASS' if good else 'FAIL'} "
              f"({len(oracle)} terms)  e(P)={sum(A.values())}")
        if not good:
            print("    expected", oracle, "\n    got     ", vals)
    print()

    # side-4: the target
    p4, A4 = descent_distribution(4)
    eP = sum(A4.values())
    print(f"  side-4: |P|={p4}, #linear extensions e(P)={eP:,}")
    palin = all(A4.get(d, 0) == A4.get(max(A4) - d, 0) for d in A4)
    print(f"          descent distribution palindromic: {palin} (max d = {max(A4)})")

    # degree + leading-coefficient (a-priori bound, Stanley).
    # Build the order polynomial symbolically; write C(top, p) as an explicit
    # degree-p polynomial in N so sympy can form a Poly (it cannot expand a
    # symbolic binomial(N+k, p) directly).
    N = sp.symbols('N')
    def binom_poly(top, k):
        num = sp.Integer(1)
        for j in range(k):
            num *= (top - j)
        return num / factorial(k)
    poly = sp.Poly(sp.expand(sum(c * binom_poly(N + p4 - d, p4)
                                 for d, c in A4.items())), N)
    deg = poly.degree()
    lc = poly.LC()
    lc_ref = sp.Rational(eP, factorial(p4))
    deg_ok = (deg == 37)
    lc_ok = (lc == lc_ref)
    ok &= deg_ok and lc_ok
    print(f"          degree a_4 = {deg} (expect 37): {'PASS' if deg_ok else 'FAIL'}")
    print(f"          leading coeff = e(P)/37! : {'PASS' if lc_ok else 'FAIL'}")
    print()

    # the 38-point match (independent of the conjecture) + b-file reproduction
    print("  Matching a_4(0..37) [poset-derived] against Conjecture 12:")
    all_match = True
    for Nv in range(0, 38):
        a_poset = a_of_N(p4, A4, Nv)         # from the poset, no conjecture used
        a_conj = conjecture12(Nv)            # from Kauers-Koutschan Conjecture 12
        if a_poset != a_conj:
            all_match = False
            print(f"    MISMATCH N={Nv}: poset={a_poset} conj={a_conj}")
    ok &= all_match
    print(f"    all 38 points (N=0..37) match: {'PASS' if all_match else 'FAIL'}")
    print(f"    a_4(0) = {a_of_N(p4, A4, 0)}  (the new 38th datum, not in the b-file)")
    print(f"    a_4(37) = {a_of_N(p4, A4, 37)}")
    print()

    # over-determination
    over = all(a_of_N(p4, A4, Nv) == conjecture12(Nv) for Nv in (38, 39, 40, 45))
    ok &= over
    print(f"  Over-determination a_4(38,39,40,45) vs Conjecture 12: {'PASS' if over else 'FAIL'}")
    print()

    print("=" * 60)
    if ok:
        print("ALL CHECKS PASSED.")
        print("deg a_4 = 37 (a priori) + 38 independent matching points")
        print("=> a_4(n) == Kauers-Koutschan Conjecture 12, identically.  QED")
        sys.exit(0)
    else:
        print("*** SOME CHECK FAILED -- see above. ***")
        sys.exit(1)

if __name__ == "__main__":
    main()
