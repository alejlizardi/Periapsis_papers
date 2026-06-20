import time
from pc import build
from family import T_family
from inflate import inflate

def rho_dedup(sigma, cap=None):
    """Exact longest PC path via iterative DFS with state dedup.
    If cap given, abort early returning cap+ if exceeded (lemma violation)."""
    n = len(sigma)
    P, V, inv = build(sigma)
    nb = (P, V)
    best = 1
    seen = set()
    stack = []
    for s in range(n):
        for c in (0, 1):
            for w in nb[c][s]:
                st = ((1 << s) | (1 << w), w, c)
                if st not in seen:
                    seen.add(st); stack.append(st)
    while stack:
        mask, u, lc = mask_u_lc = stack.pop()
        cnt = bin(mask).count("1")
        if cnt > best:
            best = cnt
            if cap is not None and best > cap:
                return best, len(seen), True
        nxt = 1 - lc
        for w in nb[nxt][u]:
            if not (mask >> w) & 1:
                st = (mask | (1 << w), w, nxt)
                if st not in seen:
                    seen.add(st); stack.append(st)
    return best, len(seen), False

if __name__ == "__main__":
    for (a, b) in [(3,6),(4,6),(5,6),(4,8),(6,8),(5,10),(8,8),(6,10),(10,8)]:
        tau = tuple(range(a))
        s = inflate(tau, [T_family(b)]*a)
        n = len(s)
        bound = 6*a + 2*b
        t0 = time.time()
        r, states, violated = rho_dedup(s, cap=bound)
        flag = "VIOLATION!!" if violated else ("bound vacuous" if bound >= n else "ok")
        print(f"a={a} b={b} n={n}: rho={r}  bound={bound}  rho/n={r/n:.3f}  "
              f"[{states} states, {time.time()-t0:.1f}s] {flag}", flush=True)
