"""Core solvers for PC paths in H_sigma (union of two Hamiltonian paths).

Vertices 0..n-1. P-edges {i,i+1}. V-edges {inv[v], inv[v+1]} where inv = sigma^{-1}.
A PC path alternates edge colors. Doubled pairs (both P and V) are parallel edges:
a traversal commits to one color.
Colors: 0 = P, 1 = V.
"""
import sys
sys.setrecursionlimit(100000)

def build(sigma):
    n = len(sigma)
    inv = [0]*n
    for i, v in enumerate(sigma):
        inv[v] = i
    # nbrs[color][vertex] = list of neighbors via that color
    P = [[] for _ in range(n)]
    V = [[] for _ in range(n)]
    for i in range(n-1):
        P[i].append(i+1); P[i+1].append(i)
    for v in range(n-1):
        a, b = inv[v], inv[v+1]
        V[a].append(b); V[b].append(a)
    return P, V, inv

def rho(sigma):
    """Longest PC path (number of vertices), exhaustive DFS."""
    n = len(sigma)
    P, V, inv = build(sigma)
    nb = (P, V)
    best = 1
    def dfs(u, last_color, mask, cnt):
        nonlocal best
        if cnt > best:
            best = cnt
        for c in (0, 1):
            if c == last_color:
                continue
            for w in nb[c][u]:
                if not (mask >> w) & 1:
                    dfs(w, c, mask | (1 << w), cnt + 1)
    for s in range(n):
        # first edge may be either color: encode "no last color" as -1
        dfs(s, -1, 1 << s, 1)
    return best

def constrained_longest(sigma, x, first_color, y=None, last_color=None):
    """Max vertices over PC paths starting at x whose first edge has color
    first_color; if y is not None, path must end at y with final edge color
    last_color. Returns 0 if no such path (a bare vertex doesn't count unless
    x==y handled by caller)."""
    n = len(sigma)
    P, V, inv = build(sigma)
    nb = (P, V)
    best = 0
    def dfs(u, last_c, mask, cnt):
        nonlocal best
        if y is None:
            if cnt > best:
                best = cnt
        else:
            if u == y and last_c == last_color and cnt > best:
                best = cnt
        nxt = 1 - last_c
        for w in nb[nxt][u]:
            if not (mask >> w) & 1:
                dfs(w, nxt, mask | (1 << w), cnt + 1)
        # also allow same color as last? no: PC condition forbids
    # start: must move with first_color
    for w in nb[first_color][x]:
        dfs(w, first_color, (1 << x) | (1 << w), 2)
    return best

def terminals(sigma):
    """Returns list of (role, vertex, cross_color)."""
    n = len(sigma)
    inv = [0]*n
    for i, v in enumerate(sigma):
        inv[v] = i
    return [("p0", 0, 0), ("p1", n-1, 0), ("v0", inv[0], 1), ("v1", inv[n-1], 1)]

def transit_M(sigma, detail=False):
    """M(sigma): max vertices collectible by an interior visit:
    enter at terminal (x, cross color c) -> first internal edge color != c;
    exit at terminal (y, cross color c') -> last internal edge color != c'.
    Includes 1-vertex visits when x==y and c != c'."""
    ts = terminals(sigma)
    best = 0
    arg = None
    for i in range(len(ts)):
        for j in range(i+1, len(ts)):
            r1, x, c1 = ts[i]
            r2, y, c2 = ts[j]
            if x == y:
                if c1 != c2 and 1 > best:
                    best = 1; arg = (r1, r2, 1)
                continue
            v = constrained_longest(sigma, x, 1 - c1, y, 1 - c2)
            if v > best:
                best = v; arg = (r1, r2, v)
    return (best, arg) if detail else best

def shallow_S(sigma):
    """S(sigma): max vertices over PC paths starting at a terminal with the
    forced first color (free end)."""
    ts = terminals(sigma)
    best = 1
    for r, x, c in ts:
        v = constrained_longest(sigma, x, 1 - c)
        if v > best:
            best = v
    return best

# ---- paper's permutations ----
def sigma_m(m):
    """Block-reversal family at n = 4m+3: blocks (2, 3,1,3,1,...,3, 2) reversed."""
    sizes = [2]
    for j in range(m):
        sizes.append(3)
        if j < m-1:
            sizes.append(1)
    sizes.append(2)
    out = []
    pos = 0
    for s in sizes:
        out.extend(range(pos+s-1, pos-1, -1))
        pos += s
    return tuple(out)

if __name__ == "__main__":
    import itertools, time
    # validate sigma_m
    for m in (1, 2, 3, 4):
        s = sigma_m(m)
        t0 = time.time()
        r = rho(s)
        print(f"sigma_{m} n={len(s)} = {s}: rho = {r} (expect {2*m+3})  [{time.time()-t0:.1f}s]")
    # validate sigma* at n=18
    sstar = (1,0,4,3,2,5,8,7,6,11,10,9,12,15,14,13,17,16)
    t0 = time.time()
    print("sigma* rho =", rho(sstar), "(expect 10)", f"[{time.time()-t0:.1f}s]")
