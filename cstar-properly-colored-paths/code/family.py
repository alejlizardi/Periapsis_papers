"""Test the candidate low-transit family T_b = (0, b-1, 2, 3, ..., b-2, 1)."""
import time
from pc import build, terminals

def constrained_longest_dedup(sigma, x, first_color, y=None, last_color=None):
    """Max vertices over PC paths from x (first edge color = first_color),
    optionally ending at y with last edge color = last_color.
    Iterative DFS with seen-state dedup on (mask, u, last_color)."""
    n = len(sigma)
    P, V, inv = build(sigma)
    nb = (P, V)
    best = 0
    seen = set()
    stack = []
    for w in nb[first_color][x]:
        st = ((1 << x) | (1 << w), w, first_color)
        if st not in seen:
            seen.add(st)
            stack.append(st)
    while stack:
        mask, u, lc = stack.pop()
        cnt = bin(mask).count("1")
        if y is None:
            if cnt > best: best = cnt
        else:
            if u == y and lc == last_color and cnt > best: best = cnt
        nxt = 1 - lc
        for w in nb[nxt][u]:
            if not (mask >> w) & 1:
                st = (mask | (1 << w), w, nxt)
                if st not in seen:
                    seen.add(st)
                    stack.append(st)
    return best

def transit_M_dedup(sigma, detail=False):
    ts = terminals(sigma)
    best = 0; arg = None
    for i in range(4):
        for j in range(i+1, 4):
            r1, xx, c1 = ts[i]; r2, yy, c2 = ts[j]
            if xx == yy:
                if c1 != c2 and best < 1:
                    best = 1; arg = (r1, r2, 1)
                continue
            v = constrained_longest_dedup(sigma, xx, 1-c1, yy, 1-c2)
            if v > best:
                best = v; arg = (r1, r2, v)
    return (best, arg) if detail else best

def shallow_S_dedup(sigma):
    ts = terminals(sigma)
    best = 1
    for r, xx, c in ts:
        v = constrained_longest_dedup(sigma, xx, 1-c)
        if v > best: best = v
    return best

def T_family(b):
    """identity with values at positions 1 and b-2... no: swap values 1 <-> b-1
    placed at positions 1 and b-1: (0, b-1, 2, 3, ..., b-2, 1)"""
    out = list(range(b))
    out[1], out[b-1] = b-1, 1
    return tuple(out)

if __name__ == "__main__":
    for b in (8, 10, 12, 16, 20, 28, 40):
        s = T_family(b)
        t0 = time.time()
        m, arg = transit_M_dedup(s, detail=True)
        S = shallow_S_dedup(s)
        print(f"b={b}: M(T_b) = {m} via {arg}, S = {S}  [{time.time()-t0:.2f}s]")
