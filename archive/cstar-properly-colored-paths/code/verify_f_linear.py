"""The PC bound does NOT transfer to the matching-edge parameter: on the very
instances sigma_{a,b} where rho = Theta(sqrt n), the graph G(2n, sigma) (top path
v_0..v_{n-1}, bottom path u_0..u_{n-1}, matching v_i -- u_{sigma(i)}) contains a
simple path using at least a(b-3) = n - 3a matching edges. Construction: snake
through the identity ladder of each block (positions 2..b-2, where sigma(p) = p),
and travel between blocks along three consecutive top- or bottom-path edges
(values/positions b-1, 0, 1 of the gap, whose rungs the snake never uses).
This script builds the path explicitly and verifies it edge by edge.
"""
import sys

def sigma_ab(a, b):
    out = []
    for j in range(a):
        blk = list(range(j * b, (j + 1) * b))
        blk[1], blk[b - 1] = j * b + b - 1, j * b + 1
        out.extend(blk)
    return out

def check(a, b):
    n = a * b
    sig = sigma_ab(a, b)
    V = lambda i: ("v", i)          # top vertex at position i
    U = lambda w: ("u", w)          # bottom vertex at value w
    top    = {frozenset({V(i), V(i + 1)}) for i in range(n - 1)}
    bottom = {frozenset({U(w), U(w + 1)}) for w in range(n - 1)}
    rungs  = {frozenset({V(i), U(sig[i])}) for i in range(n)}

    P = []
    for j in range(a):
        base = j * b
        if j % 2 == 0:                                  # enter on top side
            for l in range(2, b - 1):
                if l % 2 == 0:
                    P += [V(base + l), U(base + l)]     # rung down
                else:
                    P += [U(base + l), V(base + l)]     # rung up
        else:                                           # mirrored: enter on bottom
            for l in range(2, b - 1):
                if l % 2 == 0:
                    P += [U(base + l), V(base + l)]
                else:
                    P += [V(base + l), U(base + l)]
        if j < a - 1:                                   # travel across the gap
            side = P[-1][0]
            mk = U if side == "u" else V
            for w in (base + b - 1, base + b, base + b + 1):
                P.append(mk(w))

    # validate: simple path, every edge exists; count matching edges
    ok = len(set(P)) == len(P)
    m = 0
    for t in range(len(P) - 1):
        e = frozenset({P[t], P[t + 1]})
        if e in rungs:
            m += 1
        elif e not in top and e not in bottom:
            ok = False
            print(f"  bad edge {P[t]} -- {P[t+1]}")
    want = a * (b - 3)
    ok &= (m >= want)
    print(f"a={a:3d} b={b:3d} n={n:4d}: path on {len(P)} vertices, "
          f"{m} matching edges (claim >= {want} = n - 3a)  "
          f"{'VALID' if ok else 'INVALID'}")
    return ok

if __name__ == "__main__":
    ok = all(check(a, b) for (a, b) in
             [(3, 6), (4, 8), (10, 8), (12, 12), (20, 20), (30, 20)])
    print("ALL CHECKS PASSED" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)
