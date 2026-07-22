"""verify_n14.py -- the n = 14 threshold witness (paper Theorem "thm:n14").

Checks, for sigma-dagger = (2,0,3,1,7,5,9,4,8,6,12,10,13,11) in S_14:
  1. block-freeness: no consecutive positions carry consecutive values, so
     b(sigma) = 14 (every +-1-block a singleton);
  2. rho(sigma-dagger) = 13, by an exact longest-PC-path search over the state
     space (visited set, current vertex, last color), every reachable state
     visited once (same algorithm as endtoend.py; exact by construction).
The exhaustive censuses behind parts (b)-(c) of the theorem are separate
programs: cpp/bfcensus.cpp (all block-free sigma in S_13 are PC-Hamiltonian;
log logs/blockfree13.txt) and n14_doublebad.py with word_bruteforce.py (the
exact count of 16 failures at n = 14).
"""
SIGMA = (2, 0, 3, 1, 7, 5, 9, 4, 8, 6, 12, 10, 13, 11)

def rho_state_dedup(sigma):
    n = len(sigma)
    inv = [0] * n
    for i, v in enumerate(sigma):
        inv[v] = i
    nb = ([[] for _ in range(n)], [[] for _ in range(n)])  # 0 = P, 1 = V
    for i in range(n - 1):
        nb[0][i].append(i + 1); nb[0][i + 1].append(i)
    for v in range(n - 1):
        a, b = inv[v], inv[v + 1]
        nb[1][a].append(b); nb[1][b].append(a)
    seen = set()
    stack = []
    best = 1
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
    return best, len(seen)

if __name__ == "__main__":
    n = len(SIGMA)
    diffs = [SIGMA[i + 1] - SIGMA[i] for i in range(n - 1)]
    blockfree = all(abs(d) != 1 for d in diffs)
    print(f"sigma-dagger = {','.join(map(str, SIGMA))}")
    print(f"consecutive value differences: {diffs}")
    print(f"block-free (no +-1 difference): {blockfree}  =>  b(sigma) = {n}")
    r, states = rho_state_dedup(SIGMA)
    print(f"rho(sigma-dagger) = {r}  (expect 13; {states} reachable states)")
    print(f"rho >= b fails at n = 14: {r < n}")
