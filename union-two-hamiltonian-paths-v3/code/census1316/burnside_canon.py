#!/usr/bin/env python3
"""Independent coverage certificate for the WP-K counting censuses.

The censuses report `canon` = the number of canonical representatives they met
while scanning ALL of [0, n!).  That equals the number of orbits of S_n under
the 8-element symmetry group G = {sigma -> w^a sigma^eps w^b}, w(i) = n-1-i.
Burnside gives that orbit count in closed form, with NO enumeration:

  |orbits| = (1/8) [ n! + Fix(w s w) + 2*Fix(s = s^-1) + 2*Fix(s^2 = w) ]

  (the two elements s -> w s and s -> s w have no fixed points for n >= 2)

  Fix(w s w)  : centrally symmetric perms       = 2^floor(n/2) * floor(n/2)!
  Fix(s=s^-1) : involutions (order <= 2)        = I(n)
  Fix(s^2=w)  : square roots of w               = see sqrt_w()

If a census skipped or double-counted any rank range, its `canon` would differ
from this value.  Agreement to the digit certifies exact-once coverage.
"""
import math

def involutions(n):
    a, b = 1, 1          # I(0), I(1)
    for k in range(2, n + 1):
        a, b = b, b + (k - 1) * a
    return b

def centrally_symmetric(n):
    h = n // 2
    return 2**h * math.factorial(h)

def sqrt_w(n):
    """#{sigma : sigma^2 = w}. w = product of floor(n/2) 2-cycles (+1 fixed pt if n odd).
    Squaring maps 4-cycles -> two 2-cycles; 2-cycles -> two fixed points (forbidden
    unless they supply w's single odd fixed point, which a 1-cycle already does).
    So sigma = (one 1-cycle if n odd) + all 4-cycles, needing floor(n/2) even."""
    k = n // 2                       # number of 2-cycles in w
    if k % 2:                        # cannot pair them into 4-cycles
        return 0
    # pair the k 2-cycles: k!/(2^(k/2) (k/2)!) pairings, 2 interleavings each
    return math.factorial(k) // math.factorial(k // 2)

def orbits(n):
    return (math.factorial(n) + centrally_symmetric(n)
            + 2 * involutions(n) + 2 * sqrt_w(n)) // 8

# canon values reported by the censuses (WP-K RESULT records / audit reproductions)
REPORTED = {
    10: 456454,
    11: 4999004,
    12: 59916028,
    13: 778525516,
    14: 10897964660,
    15: 163461964024,
    16: 2615361578344,
}

if __name__ == "__main__":
    # sanity: brute-force the orbit count for small n against the formula
    from itertools import permutations
    def brute(n):
        w = lambda p: tuple(n - 1 - x for x in p)
        rev = lambda p: tuple(reversed(p))
        inv = lambda p: tuple(sorted(range(n), key=lambda i: p[i]))
        seen, orb = set(), 0
        for p in permutations(range(n)):
            if p in seen:
                continue
            orb += 1
            imgs = {p}
            for _ in range(3):
                imgs |= {w(q) for q in imgs} | {rev(q) for q in imgs} | {inv(q) for q in imgs}
            seen |= imgs
        return orb
    print("brute-force validation of the Burnside formula:")
    for n in range(2, 9):
        b, f = brute(n), orbits(n)
        print(f"  n={n:2d}  brute={b:8d}  burnside={f:8d}  {'OK' if b == f else 'MISMATCH'}")
        assert b == f, n

    print("\ncensus `canon` vs exact orbit count:")
    allok = True
    for n, rep in sorted(REPORTED.items()):
        exact = orbits(n)
        ok = (exact == rep)
        allok &= ok
        print(f"  n={n}  reported={rep:>16,}  exact={exact:>16,}  {'MATCH' if ok else 'MISMATCH'}")
    print("\nCOVERAGE CERTIFICATE:", "ALL MATCH" if allok else "FAILED")
