"""Cross-check the trap-chain family from PROOF_upper_bound.md with this
session's independent solver."""
from endtoend import rho_dedup

def sigma_C(comp):
    out = []
    pos = 0
    for c in comp:
        out.extend(range(pos + c - 1, pos - 1, -1))
        pos += c
    return tuple(out)

for L in (5, 7):
    for k in (2, 3, 4):
        comp = [2] + [L, 1] * k + [3, 2]
        s = sigma_C(comp)
        n = len(s)
        r, states, _ = rho_dedup(s)
        pred = 2 * k + 2 * (L - 1)
        print(f"L={L} k={k} n={n}: rho={r} predicted={pred} {'OK' if r == pred else 'MISMATCH'}",
              flush=True)

# also cross-check b(sigma) refutation point at (a,b)=(20,12) once more
from family import T_family
from inflate import inflate
s = inflate(tuple(range(20)), [T_family(12)] * 20)
r, _, _ = rho_dedup(s)
print(f"sigma(20,12): n={len(s)} rho={r} (b(sigma)=80; rho<b: {r < 80})")
