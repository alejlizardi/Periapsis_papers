"""Quick local-optimality probe: hill-climb from the inflated construction,
plus short greedy descents from random starts."""
import random, time
from endtoend import rho_dedup
from family import T_family
from inflate import inflate

def hill(perm, budget=400, rnd=None):
    rnd = rnd or random.Random(0)
    cur = list(perm)
    cur_r, _, _ = rho_dedup(tuple(cur))
    for _ in range(budget):
        i = rnd.randrange(len(cur)-1)
        cur[i], cur[i+1] = cur[i+1], cur[i]
        r, _, _ = rho_dedup(tuple(cur))
        if r <= cur_r:
            cur_r = r
        else:
            cur[i], cur[i+1] = cur[i+1], cur[i]
    return cur_r, tuple(cur)

if __name__ == "__main__":
    rnd = random.Random(7)
    for (a, b) in [(4, 6), (4, 8)]:
        n = a*b
        infl = inflate(tuple(range(a)), [T_family(b)]*a)
        ir, _, _ = rho_dedup(infl)
        hr, _ = hill(infl, budget=300, rnd=rnd)
        print(f"n={n}: inflated rho={ir}; hill-climb from it reaches {hr}", flush=True)
        best = None
        t0 = time.time()
        for s in range(4):
            p = list(range(n)); rnd.shuffle(p)
            r, _ = hill(p, budget=250, rnd=rnd)
            best = r if best is None else min(best, r)
            if time.time() - t0 > 240: break
        print(f"n={n}: best from {s+1} random greedy descents: {best}", flush=True)
