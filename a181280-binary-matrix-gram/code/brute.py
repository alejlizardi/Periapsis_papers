#!/usr/bin/env python3
"""
A181280 brute-force counter — ground truth for the modeling link.

a(n) = #{ M in Z2^{4 x n} :
            rows of M (as bit strings) lexicographically strictly INCREASING,
            rows of (M M^T mod 2) in Z2^{4x4} lexicographically strictly DECREASING }.

We enumerate all 16^n column choices for small n (n <= 6 feasible: 16^6 ~ 1.6e7),
build M, test both conditions exactly. This is the OBJECT, modeled as literally as the
definition reads. Compare to OEIS data:
  n: 1  2  3  4    5     6      7       8         ...
  a: 0  0  0  58   1629  28924  507052  8211776   ...
(offset: OEIS data starts at n=1.)
"""
import itertools

OEIS = {1:0, 2:0, 3:0, 4:58, 5:1629, 6:28924, 7:507052}

def row_bitstring(M, i, n):
    # row i as an integer with column 1 = most significant bit (lexicographic on bit strings)
    v = 0
    for k in range(n):
        v = (v << 1) | M[i][k]
    return v

def gram_row_bitstring(G, i):
    # G is 4x4; row i as 4-bit int, column 0 = MSB
    v = 0
    for j in range(4):
        v = (v << 1) | G[i][j]
    return v

def count(n):
    cols = list(itertools.product((0,1), repeat=4))  # 16 columns, each a 4-vector
    total = 0
    for choice in itertools.product(cols, repeat=n):
        # M[i][k]
        M = [[choice[k][i] for k in range(n)] for i in range(4)]
        # rows strictly increasing
        r = [row_bitstring(M, i, n) for i in range(4)]
        if not (r[0] < r[1] < r[2] < r[3]):
            continue
        # Gram mod 2
        G = [[0]*4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                s = 0
                for k in range(n):
                    s ^= (M[i][k] & M[j][k])
                G[i][j] = s
        g = [gram_row_bitstring(G, i) for i in range(4)]
        if not (g[0] > g[1] > g[2] > g[3]):
            continue
        total += 1
    return total

if __name__ == "__main__":
    for n in range(1, 7):
        c = count(n)
        exp = OEIS.get(n, "?")
        ok = "OK" if exp == c else ("??" if exp == "?" else "MISMATCH")
        print(f"n={n}: brute={c}  OEIS={exp}  [{ok}]")
