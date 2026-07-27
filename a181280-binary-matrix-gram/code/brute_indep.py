#!/usr/bin/env python3
"""
INDEPENDENT brute force for A181280 -- written from scratch, DIFFERENT code path
than the prover's brute.py.

Prover's brute.py iterates over 16^n COLUMN choices and builds M column-wise.
Here we iterate over the 4 ROWS directly as integers r0<r1<r2<r3 in [0, 2^n),
which automatically enforces the strict-increasing condition (column 1 = MSB
<=> the integer value, since we read the bit string MSB-first as the integer).

Gram entry G_{ij} = parity of (# columns where both rows have a 1) = popcount(r_i & r_j) mod 2.
Then read each Gram row as a 4-bit int (col 0 = MSB) and require strictly decreasing.

Because r0<r1<r2<r3 we only enumerate the C(2^n,4) increasing 4-tuples, much fewer
than 16^n, so this reaches n=6,7 comfortably.
"""
from itertools import combinations

OEIS = {4:58, 5:1629, 6:28924, 7:507052}

def gram_rows_strictly_decreasing(rows):
    # rows = (r0,r1,r2,r3) integers; G symmetric 4x4 parity matrix
    G = [[0]*4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            G[i][j] = bin(rows[i] & rows[j]).count("1") & 1
    grow = []
    for i in range(4):
        v = 0
        for j in range(4):
            v = (v << 1) | G[i][j]
        grow.append(v)
    return grow[0] > grow[1] > grow[2] > grow[3]

def count(n):
    N = 1 << n
    total = 0
    # strictly increasing 4-tuples of row-integers
    for rows in combinations(range(N), 4):
        if gram_rows_strictly_decreasing(rows):
            total += 1
    return total

if __name__ == "__main__":
    for n in range(4, 8):
        c = count(n)
        exp = OEIS.get(n, "?")
        ok = "OK" if exp == c else ("?" if exp == "?" else "MISMATCH")
        print(f"n={n}: brute_indep={c}  OEIS={exp}  [{ok}]")
