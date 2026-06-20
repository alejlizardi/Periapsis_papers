import random, time
from pc import transit_M, rho
from family import transit_M_dedup, constrained_longest_dedup, T_family
from inflate import inflate, M2

# 1) cross-validate dedup vs original on random perms
random.seed(1)
ok = True
for trial in range(300):
    b = random.randint(4, 8)
    p = list(range(b)); random.shuffle(p); p = tuple(p)
    if transit_M(p) != transit_M_dedup(p):
        print("MISMATCH", p); ok = False
print("dedup cross-validation:", "OK" if ok else "FAILED")

# 2) odd b should break the family (expect large M)
for b in (7, 9, 11, 13):
    print(f"odd b={b}: M(T_b) = {transit_M_dedup(T_family(b))} (expect ~b)")

# 3) M2 for even b
for b in (6, 8, 10, 12):
    s = T_family(b)
    print(f"even b={b}: M={transit_M_dedup(s)}, M2={M2(s)} (expect 3, 4)")
