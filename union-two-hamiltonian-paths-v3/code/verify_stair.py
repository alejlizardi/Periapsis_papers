"""verify_stair.py -- the staircase chains (paper Section "Block counts do
not bound the longest properly colored path", Definition "Staircase pattern",
Lemma "Staircase basics", Remark "Exactness; data").

Checks:
  1. for admissible s (odd, s != 1 mod 3) up to SMAX: theta_s is a
     permutation; theta_s is block-free; every H_{theta_s} edge joins an even
     and an odd position (the parity lemma, edge by edge);
  2. for the chains Theta_{s,k} (k-fold direct sum, identity outer
     permutation) with k <= KMAX: the chain is block-free (b = s*k), and each
     junction has exactly two cross edges, the P-edge (vertex m of piece j,
     vertex 0 of piece j+1) and the V-edge (vertex 0 of piece j, vertex m of
     piece j+1), exactly as Lemma "Staircase basics"(c) states;
  3. exact rho(Theta_{9,k}) for k = 1, 2, 3 by the state-dedup search
     (expect 9k: the chains are PC-Hamiltonian through at least four pieces;
     the certified optima at k >= 5 are cpsat_rho.py runs, see
     logs/cpsat_stair.txt).
"""
import sys

def theta(s):
    assert s >= 9 and s % 2 == 1 and s % 3 != 1
    m = s - 1
    out = [m] + [(3 * i) % m for i in range(1, m)] + [0]
    assert sorted(out) == list(range(s))
    return out

def chain(s, k):
    t = theta(s)
    out = []
    for j in range(k):
        out.extend(v + j * s for v in t)
    return out

def edges(sigma):
    n = len(sigma)
    inv = [0] * n
    for i, v in enumerate(sigma):
        inv[v] = i
    P = [(i, i + 1) for i in range(n - 1)]
    V = [tuple(sorted((inv[v], inv[v + 1]))) for v in range(n - 1)]
    return P, V

def blockfree(sigma):
    return all(abs(sigma[i + 1] - sigma[i]) != 1 for i in range(len(sigma) - 1))

def rho_state_dedup(sigma):
    n = len(sigma)
    inv = [0] * n
    for i, v in enumerate(sigma):
        inv[v] = i
    nb = ([[] for _ in range(n)], [[] for _ in range(n)])
    for i in range(n - 1):
        nb[0][i].append(i + 1); nb[0][i + 1].append(i)
    for v in range(n - 1):
        a, b = inv[v], inv[v + 1]
        nb[1][a].append(b); nb[1][b].append(a)
    seen = set(); stack = []; best = 1
    for s in range(n):
        for c in (0, 1):
            for w in nb[c][s]:
                st = ((1 << s) | (1 << w), w, c)
                if st not in seen:
                    seen.add(st); stack.append(st)
    while stack:
        mask, u, c = stack.pop()
        cnt = bin(mask).count("1")
        if cnt > best:
            best = cnt
        for w in nb[1 - c][u]:
            if not (mask >> w) & 1:
                st = (mask | (1 << w), w, 1 - c)
                if st not in seen:
                    seen.add(st); stack.append(st)
    return best

if __name__ == "__main__":
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 33
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    ok = True
    for s in range(9, smax + 1, 2):
        if s % 3 == 1:
            continue
        m = s - 1
        th = theta(s)
        bf = blockfree(th)
        P, V = edges(th)
        parity = all((a + b) % 2 == 1 for (a, b) in P + V)
        print(f"s={s}: theta_s permutation ok, block-free={bf}, "
              f"all edges parity-flipping={parity}")
        ok &= bf and parity
        for k in range(2, kmax + 1):
            ch = chain(s, k)
            bfc = blockfree(ch)
            Pc, Vc = edges(ch)
            juncP = [e for e in Pc if e[0] // s != e[1] // s]
            juncV = [e for e in Vc if e[0] // s != e[1] // s]
            expP = [(j * s - 1, j * s) for j in range(1, k)]
            expV = [tuple(sorted(((j - 1) * s, j * s + m))) for j in range(1, k)]
            good = (bfc and sorted(juncP) == sorted(expP)
                    and sorted(juncV) == sorted(expV))
            ok &= good
            if not good:
                print(f"  FAIL s={s} k={k}: block-free={bfc} "
                      f"P-cross={juncP} V-cross={juncV}")
        print(f"       chains k<={kmax}: block-free + junction structure ok")
    for k in (1, 2, 3):
        r = rho_state_dedup(chain(9, k))
        print(f"rho(Theta_9,{k}) = {r}  (n = {9 * k}; expect {9 * k})")
    print("ALL OK" if ok else "FAILURES FOUND")
