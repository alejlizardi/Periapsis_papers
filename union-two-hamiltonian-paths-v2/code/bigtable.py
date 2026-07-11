import time
from family import T_family, transit_M_dedup
from inflate import inflate
from endtoend import rho_dedup

# M robustness at larger even b
for b in (50, 60, 80, 100):
    print(f"M(T_{b}) = {transit_M_dedup(T_family(b))}", flush=True)

print()
rows = []
for (a, b) in [(10,8),(12,8),(12,12),(16,12),(20,12),(16,16),(24,16),(20,20),(30,16),(24,20),(30,20),(36,20)]:
    tau = tuple(range(a))
    s = inflate(tau, [T_family(b)]*a)
    n = len(s)
    t0 = time.time()
    r, states, _ = rho_dedup(s)
    pred = 2*a + 2*b - 4
    print(f"a={a:3d} b={b:3d} n={n:4d}: rho={r:4d}  pred(2a+2b-4)={pred:4d}  "
          f"bound6a2b={6*a+2*b:4d}  rho/n={r/n:.3f}  rho/sqrt(n)={r/n**0.5:.2f}  "
          f"[{states} states {time.time()-t0:.1f}s]", flush=True)
