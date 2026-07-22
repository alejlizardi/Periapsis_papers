"""Verify the explicit PC path of Proposition 6.3 (rho(sigma_{a,b}) >= 2a+2b-7):
stage 1: (0,b-2) -V- (0,b-3) -P- ... alternating down to (0,2), then V to (0,b-1),
         P-cross to (1,0);
stage 2: for j = 1..a-2: (j,0) -V- (j,b-1) -P-cross- (j+1,0);
stage 3: (a-1,0) -V- (a-1,b-1) -P- (a-1,b-2) -V- (a-1,b-3) -P- ... down to (a-1,2).
Checks: vertices distinct, every consecutive pair is an edge of the claimed color,
colors alternate, and the count equals 2a+2b-7.
"""
import sys

def sigma_ab(a, b):
    out = []
    for j in range(a):
        blk = list(range(j * b, (j + 1) * b))
        blk[1], blk[b - 1] = j * b + b - 1, j * b + 1
        out.extend(blk)
    return out

def edge_colors(sig):
    n = len(sig)
    inv = [0] * n
    for i, v in enumerate(sig):
        inv[v] = i
    E = set()
    for i in range(n - 1):
        E.add((min(i, i + 1), max(i, i + 1), 0))
    for v in range(n - 1):
        x, y = inv[v], inv[v + 1]
        E.add((min(x, y), max(x, y), 1))
    return E

def build_path(a, b):
    P, C = [], []           # vertices; colors of edges
    g = lambda j, l: j * b + l
    # stage 1
    for l in range(b - 2, 1, -1):
        P.append(g(0, l))
    C.extend([1, 0] * ((b - 4) // 2))           # V,P,...,P  (b-4 edges)
    P.append(g(0, b - 1)); C.append(1)          # V {2, b-1}
    # stage 2 (express)
    for j in range(1, a - 1):
        P.append(g(j, 0));     C.append(0)      # P-cross
        P.append(g(j, b - 1)); C.append(1)      # V {0, b-1}
    # stage 3
    P.append(g(a - 1, 0));     C.append(0)      # P-cross
    P.append(g(a - 1, b - 1)); C.append(1)      # V {0, b-1}
    P.append(g(a - 1, b - 2)); C.append(0)      # P {b-2, b-1}
    col = 1
    for l in range(b - 3, 1, -1):
        P.append(g(a - 1, l)); C.append(col); col = 1 - col
    return P, C

if __name__ == "__main__":
    ok = True
    for (a, b) in [(2, 6), (3, 6), (4, 8), (5, 10), (10, 8), (7, 12), (20, 20)]:
        sig = sigma_ab(a, b)
        E = edge_colors(sig)
        P, C = build_path(a, b)
        good = len(set(P)) == len(P) == len(C) + 1
        for t in range(len(C)):
            u, w = P[t], P[t + 1]
            if (min(u, w), max(u, w), C[t]) not in E:
                good = False; print(f"  bad edge {u}-{w} color {C[t]}")
            if t and C[t] == C[t - 1]:
                good = False; print(f"  color repeat at {t}")
        want = 2 * a + 2 * b - 7
        good &= (len(P) == want)
        ok &= good
        print(f"a={a} b={b}: path size {len(P)} (claim {want}) "
              f"{'VALID PC PATH' if good else 'INVALID'}")
    print("ALL CHECKS PASSED" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)
