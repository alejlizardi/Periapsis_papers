"""Minimal independent verifier for the transit number M(pi) and the
classification of admissible visits of H_{T_b} (paper Definitions 3.3-3.4,
Theorem 5.4). Written directly from the definitions, sharing no code with the
main solvers. Pure brute force: enumerates ALL properly colored simple paths
between terminals. Feasible for b <= 14.

Usage:  python verify_M.py
"""
import sys

def T(b):
    out = list(range(b))
    out[1], out[b - 1] = b - 1, 1
    return tuple(out)

def adjacency(sigma):
    """List of (neighbor, color) per vertex; color 0 = P, 1 = V.
    Doubled pairs appear as two entries of different colors."""
    b = len(sigma)
    inv = [0] * b
    for i, v in enumerate(sigma):
        inv[v] = i
    A = [[] for _ in range(b)]
    for i in range(b - 1):                      # P-edges
        A[i].append((i + 1, 0)); A[i + 1].append((i, 0))
    for v in range(b - 1):                      # V-edges
        x, y = inv[v], inv[v + 1]
        A[x].append((y, 1)); A[y].append((x, 1))
    return A

def roles(sigma):
    """(name, vertex, port color) for p0, p1, v0, v1."""
    b = len(sigma)
    inv = [0] * b
    for i, v in enumerate(sigma):
        inv[v] = i
    return [("p0", 0, 0), ("p1", b - 1, 0), ("v0", inv[0], 1), ("v1", inv[b - 1], 1)]

def admissible_visits(sigma):
    """Return every admissible visit as (vertices, rolepair, firstcolor, lastcolor).
    Type (ii): single vertex carrying two roles of different port colors.
    Type (i): PC path between terminals, first edge != pc(start role),
              last edge != pc(end role)."""
    A = adjacency(sigma)
    R = roles(sigma)
    found = []
    for i in range(4):                          # type (ii)
        for j in range(i + 1, 4):
            r1, x, c1 = R[i]; r2, y, c2 = R[j]
            if x == y and c1 != c2:
                found.append(((x,), (r1, r2), None, None))
    term_vs = {v for (_, v, _) in R}
    def dfs(path, used, firstc, lastc):
        u = path[-1]
        if len(path) >= 2 and u in term_vs:
            for (rn, rv, rc) in R:              # end role options
                if rv != u or lastc == rc:
                    continue
                for (sn, sv, sc) in R:          # start role options
                    if sv == path[0] and firstc != sc and sn != rn:
                        found.append((tuple(path), (sn, rn), firstc, lastc))
        for (w, c) in A[u]:
            if w not in used and c != lastc:
                path.append(w); used.add(w)
                dfs(path, used, firstc, c)
                path.pop(); used.remove(w)
    for v0 in term_vs:
        for (w, c) in A[v0]:
            dfs([v0, w], {v0, w}, c, c)
    return found

def M(sigma):
    vis = admissible_visits(sigma)
    return max((len(v[0]) for v in vis), default=0)

def canon(vs):
    return min(tuple(vs), tuple(reversed(vs)))

if __name__ == "__main__":
    ok = True
    for b in (6, 8, 10, 12, 14):
        sig = T(b)
        m = M(sig)
        kinds = sorted({(canon(v), frozenset(rp)) for (v, rp, _, _) in admissible_visits(sig)})
        expected = sorted({
            (canon((0,)), frozenset({"p0", "v0"})),
            (canon((0, b - 1)), frozenset({"p0", "p1"})),
            (canon((0, 1)), frozenset({"v0", "v1"})),
            (canon((b - 1, 0, 1)), frozenset({"p1", "v1"})),
            (canon((b - 1, 2, 1)), frozenset({"p1", "v1"})),
        })
        match = (kinds == expected)
        ok &= (m == 3) and match
        print(f"b={b}: M(T_b)={m} (expect 3); classification "
              f"{'EXACTLY the 5 listed visits' if match else 'MISMATCH: ' + str(kinds)}")
    for b in (7, 9, 11):
        m = M(T(b))
        ok &= (m == b)
        print(f"b={b} (odd): M(T_b)={m} (expect {b})")
    print("ALL CHECKS PASSED" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)
