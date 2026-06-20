"""Inflation sigma = tau[pi_1..pi_a] and refined transit invariants."""
from pc import rho, transit_M, shallow_S, terminals, build, constrained_longest

def inflate(tau, pis):
    """tau in S_a (0-indexed tuple), pis = list of a permutations (tuples).
    All same size b assumed not required; blocks sized len(pis[i])."""
    a = len(tau)
    assert len(pis) == a
    sizes = [len(p) for p in pis]
    # value offset of block i = sum of sizes of blocks j with tau(j) < tau(i)
    order = sorted(range(a), key=lambda i: tau[i])
    voff = [0]*a
    acc = 0
    for i in order:
        voff[i] = acc
        acc += sizes[i]
    out = []
    for i in range(a):
        for r in range(sizes[i]):
            out.append(voff[i] + pis[i][r])
    return tuple(out)

def M2(sigma):
    """Max total vertices over TWO vertex-disjoint legitimate transits
    (or one, i.e. M2 >= M). Brute force: for each pair of role-pairs using
    4 distinct cross-edge slots, find max sum of disjoint constrained paths."""
    from itertools import permutations as perms
    ts = terminals(sigma)
    n = len(sigma)
    P, V, inv = build(sigma)
    nb = (P, V)

    def all_paths(x, first_color, y, last_color, forbidden):
        """Yield visited-masks of PC paths x->y meeting color constraints,
        avoiding 'forbidden' mask; record (mask, count). Collect maximal info."""
        res = []
        if (forbidden >> x) & 1 or (forbidden >> y) & 1:
            return res
        def dfs(u, last_c, mask, cnt):
            if u == y and last_c == last_color:
                res.append((mask, cnt))
            nxt = 1 - last_c
            for w in nb[nxt][u]:
                if not (mask >> w) & 1 and not (forbidden >> w) & 1:
                    dfs(w, nxt, mask | (1 << w), cnt + 1)
        for w in nb[first_color][x]:
            if not (forbidden >> w) & 1:
                dfs(w, first_color, (1 << x) | (1 << w), 2)
        return res

    roles = ts  # (name, vertex, crosscolor)
    best = 0
    # single transit
    for i in range(4):
        for j in range(i+1, 4):
            r1, x, c1 = roles[i]; r2, y, c2 = roles[j]
            if x == y:
                if c1 != c2: best = max(best, 1)
                continue
            v = constrained_longest(sigma, x, 1-c1, y, 1-c2)
            best = max(best, v)
    # two disjoint transits: partition 4 roles into two pairs (3 pairings)
    pairings = [((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2))]
    for (i,j),(k,l) in pairings:
        r1, x1, c1 = roles[i]; r2, y1, c2 = roles[j]
        r3, x2, c3 = roles[k]; r4, y2, c4 = roles[l]
        if len({x1,y1,x2,y2}) < 4:
            continue  # coinciding terminals can't host 2 disjoint visits cleanly; skip (conservative handled by 2M bound elsewhere)
        ps1 = all_paths(x1, 1-c1, y1, 1-c2, 0)
        if not ps1:
            continue
        for mask1, cnt1 in ps1:
            ps2 = all_paths(x2, 1-c3, y2, 1-c4, mask1)
            for mask2, cnt2 in ps2:
                if cnt1 + cnt2 > best:
                    best = cnt1 + cnt2
    return best

if __name__ == "__main__":
    # sanity: inflate identity
    tau = (0,1,2)
    pi = (0,1)
    s = inflate(tau, [pi]*3)
    print(s)  # expect (0,1,2,3,4,5)
    tau = (1,0)
    s = inflate(tau, [(1,0),(0,1)])
    print(s)  # block0 -> values {2,3} reversed: (3,2); block1 -> {0,1}: (0,1) => (3,2,0,1)
