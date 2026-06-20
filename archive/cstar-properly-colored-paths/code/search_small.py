import itertools, time
from pc import rho, transit_M, shallow_S

# validate rho_min(7)=5, rho_min(8)=6
for n in (7, 8):
    t0 = time.time()
    rm = min(rho(p) for p in itertools.permutations(range(n)))
    print(f"rho_min({n}) = {rm}  [{time.time()-t0:.1f}s]")

# exhaustive transit-number minimum over S_b
for b in range(4, 10):
    t0 = time.time()
    bestM = 10**9
    args = []
    for p in itertools.permutations(range(b)):
        m = transit_M(p)
        if m < bestM:
            bestM = m; args = [p]
        elif m == bestM and len(args) < 6:
            args.append(p)
    print(f"b={b}: min M = {bestM}  ratio {bestM/b:.3f}  examples: {args[:4]}  [{time.time()-t0:.1f}s]")
