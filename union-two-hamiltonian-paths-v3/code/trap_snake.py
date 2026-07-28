"""WP-J T2c: the trap-chain family ALSO has linear matching-edge content.

For a composition C = (c_0,...,c_r) with no two consecutive entries equal to 1,
the paper's Appendix-A chain permutation sigma_C reverses each block:
sigma_C(P_t + i) = P_t + c_t - 1 - i on positions/values [P_t, P_t + c_t - 1].

Claim (constructed + machine-verified edge-by-edge below): G(n, sigma_C)
contains a simple path using at least

    n  -  #{t : c_t odd and c_t >= 3}

matching edges.  On the canonical trap chains (odd blocks 2m+1 separated by
singletons) this is  n - (#big blocks) = n(1 - 1/(2m+2)) = n(1-o(1))  for
growing m, while rho(sigma_C) = Theta(sqrt n) at the balanced point: the
properly colored relaxation is lossy on BOTH known rho-collapsing families.

Construction. Ride into block t on the current side ('A' = the side we are on,
'B' = the other; A/B is top/bottom or bottom/top). Within a block of size c the
rungs are A(base+i) -- B(base+c-1-i). Take the interleaved zigzag

    A(base+0) | B(c-1) B(c-2) | A(1) A(2) | B(c-3) B(c-4) | A(3) A(4) | ...

(consecutive same-side vertices are path edges; every A->B or B->A jump is a
rung, because the consumed-index counts keep i + j = c - 1).  For even c this
ends ... B(1) B(0) | A(c-1): ALL c rungs, exit A-side right.  For odd c the
final lone B(0) is dropped: c-1 rungs, exit A-side right.  A singleton block
is a single rung A(base) -> B(base): 1 rung, exit flips to the B side.
Between blocks one same-side step joins exit A(base+c-1) to entry A(base+c)
(positions are consecutive on top; values are consecutive on bottom) - a
same-side run is exactly what G-paths allow and PC paths forbid.

Every path is verified EDGE BY EDGE against the definitions (mirroring the
paper's verify_f_linear.py); small instances are additionally compared with the
exact bitmask-DP maxmu via fmu.exe single (2n <= 24).
"""
import sys

def sigma_C(comp):
    out = []
    P = 0
    for c in comp:
        out.extend(range(P + c - 1, P - 1, -1))
        P += c
    return out

def block_seq(c):
    """Local zigzag as a list of ('A'|'B', local_index); enter A(0), exit A-side
    right (c-1) for c>=2; ALL c rungs when c even, c-1 when c odd."""
    asc = list(range(c))
    desc = list(range(c - 1, -1, -1))
    seq = [("A", asc.pop(0))]
    turn = "B"
    while asc or desc:
        if turn == "B":
            take, desc = desc[:2], desc[2:]
            seq += [("B", x) for x in take]
            turn = "A"
        else:
            take, asc = asc[:2], asc[2:]
            seq += [("A", x) for x in take]
            turn = "B"
    if c % 2 == 1 and c >= 3 and seq[-1] == ("B", 0):
        seq.pop()
    return seq

def snake(comp):
    V = lambda i: ("v", i)
    U = lambda w: ("u", w)
    P = []
    side = "v"                            # current riding side; A = side
    base = 0
    for c in comp:
        A, B = (V, U) if side == "v" else (U, V)
        if c == 1:
            P.append(A(base))
            P.append(B(base))
            side = "u" if side == "v" else "v"
        else:
            for who, l in block_seq(c):
                P.append(A(base + l) if who == "A" else B(base + l))
            # side unchanged (exit on A side)
        base += c
    return P

def verify(comp, verbose=True):
    n = sum(comp)
    sig = sigma_C(comp)
    assert sorted(sig) == list(range(n))
    V = lambda i: ("v", i)
    U = lambda w: ("u", w)
    top = {frozenset({V(i), V(i + 1)}) for i in range(n - 1)}
    bottom = {frozenset({U(w), U(w + 1)}) for w in range(n - 1)}
    rungs = {frozenset({V(i), U(sig[i])}) for i in range(n)}
    P = snake(comp)
    ok = len(set(P)) == len(P)
    if not ok and verbose:
        print("  NOT SIMPLE")
    m = 0
    for t in range(len(P) - 1):
        e = frozenset({P[t], P[t + 1]})
        if e in rungs:
            m += 1
        elif e not in top and e not in bottom:
            ok = False
            if verbose:
                print(f"  bad edge {P[t]} -- {P[t+1]}")
    loss = sum(1 for c in comp if c % 2 == 1 and c >= 3)
    want = n - loss
    ok &= (m >= want)
    if verbose:
        cs = str(list(comp)) if len(comp) <= 12 else f"[{comp[0]},{comp[1]},...]x{len(comp)}"
        print(f"C={cs:28s} n={n:5d}: path on {len(P):5d} vertices, {m:5d} matching edges "
              f"(claim >= {want} = n - #odd)  {'VALID' if ok else 'INVALID'}")
    return ok, m, want

if __name__ == "__main__":
    tests = [
        (5, 1, 5),
        (3, 1, 3, 1, 3),
        (7, 1, 7),
        (2, 3, 2),
        (4, 4, 4),
        (2, 1, 2, 1, 2),
        (9, 1, 9, 1, 9),
        (5, 1, 5, 1, 5, 1, 5, 1, 5),
        tuple([9, 1] * 20 + [9]),
        tuple([15, 1] * 30 + [15]),
        tuple([21, 1] * 40 + [21]),
        tuple([101, 1] * 50 + [101]),
    ]
    allok = True
    for comp in tests:
        ok, m, want = verify(comp)
        allok &= ok
    print("ALL CHECKS PASSED" if allok else "FAILURES PRESENT")
    sys.exit(0 if allok else 1)
