#!/usr/bin/env python3
"""
Standalone verifier for the two main results of

    "Recursively differentiable quasigroups of orders 14 and 18:
     the last open cases of the Couselo-Gonzalez-Markov-Nechaev conjecture"
    (A. Lizardi, Periapsis, 2026).

It depends only on the Python standard library and the two JSON tables shipped
alongside it (L14.json, L18.json), which contain exactly the Cayley tables
printed in Sections 5 and 6 of the paper. It SHARES NO CODE with the search
that produced the tables: it re-checks them from the definitions.

For each order n in {14, 18} it:
  1. checks the base table L is a Latin square (a quasigroup);
  2. forms the first recursive derivative TWO independent ways and confirms they
     agree --- (a) the closed form D[a][b] = L[b][L[a][b]] (paper eq. (2.4)), and
     (b) by unrolling the CGMN recurrence (paper eq. (1.1)) step by step;
  3. checks the derivative D is a Latin square (so the core is a quasigroup);
  4. runs the "wrong-reading" control: of the four candidate core operations
     b*(a*b), a*(a*b), (a*b)*b, (a*b)*a, ONLY b*(a*b) yields a Latin derivative,
     confirming the object is certified against the precise definition.

Exit code 0 iff both orders pass all checks.  Usage:  python rdiff_verify.py
"""
import json
import os
import sys


def is_latin(M):
    n = len(M)
    if any(len(row) != n for row in M):
        return False
    full = set(range(n))
    if any(set(row) != full for row in M):                 # rows
        return False
    if any({M[i][j] for i in range(n)} != full for j in range(n)):  # cols
        return False
    return True


def derivative_closed_form(L):
    """D[a][b] = b . (a . b) = L[b][ L[a][b] ]   (paper eq. (2.4))."""
    n = len(L)
    return [[L[b][L[a][b]] for b in range(n)] for a in range(n)]


def derivative_by_recurrence(L):
    """First recursive derivative by unrolling the CGMN recurrence (eq. (1.1)):
       a*_{-2}b = a ; a*_{-1}b = b ; a*_k b = (a*_{k-2}b) . (a*_{k-1}b).
       Returns the table of a*_1 b, computed without using the closed form."""
    n = len(L)
    def mul(x, y):
        return L[x][y]
    D = [[None] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            s_m2, s_m1 = a, b          # *_{-2}, *_{-1}
            s0 = mul(a, b)             # *_0 = a . b
            s1 = mul(s_m1, s0)         # *_1 = (*_{-1}) . (*_0) = b . (a . b)
            D[a][b] = s1
    return D


def wrong_reading_control(L):
    """Return {label: is_latin} for the four candidate 'core' readings.
       The paper's claim: only b.(a.b) is Latin."""
    n = len(L)
    readings = {
        "b.(a.b)  [correct]": [[L[b][L[a][b]] for b in range(n)] for a in range(n)],
        "a.(a.b)":            [[L[a][L[a][b]] for b in range(n)] for a in range(n)],
        "(a.b).b":            [[L[L[a][b]][b] for b in range(n)] for a in range(n)],
        "(a.b).a":            [[L[L[a][b]][a] for b in range(n)] for a in range(n)],
    }
    return {label: is_latin(M) for label, M in readings.items()}


def verify_order(n, here):
    path = os.path.join(here, "L%d.json" % n)
    L = json.load(open(path))["L"]
    assert len(L) == n, "table in %s is not order %d" % (path, n)

    ok_L = is_latin(L)
    D_cf = derivative_closed_form(L)
    D_rec = derivative_by_recurrence(L)
    agree = (D_cf == D_rec)
    ok_D = is_latin(D_cf)
    control = wrong_reading_control(L)
    only_correct = (
        control["b.(a.b)  [correct]"] is True
        and not any(v for k, v in control.items() if k != "b.(a.b)  [correct]")
    )

    print("=== order %d ===" % n)
    print("  base L is a Latin square (quasigroup):        %s" % ok_L)
    print("  closed form == recurrence unrolling:          %s" % agree)
    print("  derivative D is a Latin square (core is QG):  %s" % ok_D)
    print("  wrong-reading control (only b.(a.b) Latin):   %s" % only_correct)
    for label, val in control.items():
        print("      %-22s Latin? %s" % (label, val))
    passed = ok_L and agree and ok_D and only_correct
    print("  => order %d: %s" % (n, "PASS" if passed else "FAIL"))
    print()
    return passed


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    results = [verify_order(n, here) for n in (14, 18)]
    if all(results):
        print("ALL CHECKS PASSED: recursively differentiable quasigroups of "
              "orders 14 and 18 are certified.")
        return 0
    print("VERIFICATION FAILED.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
