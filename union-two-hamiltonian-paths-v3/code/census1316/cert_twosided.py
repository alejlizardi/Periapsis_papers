"""WP-K two-sided spot certification from a GPU run's reservoir samples.

Usage: py -3 cert_twosided.py <run_outdir> [per_side=100]

Reads checkpoint.json reservoirs (res_def / res_nondef = deterministic bottom-k
samples of canonical ranks the census called deficient / non-deficient).

Deficient side  : pc.py exact rho (state-dedup DFS) AND CP-SAT arc-MTZ (different
                  formalism) must BOTH return rho <= n-1, CP-SAT status OPTIMAL,
                  and agree with each other.
Non-deficient   : an explicit properly-colored Hamiltonian path is produced by a
                  dedicated DFS and re-verified by an independent structural
                  checker (adjacency rebuilt from scratch, alternation + coverage
                  asserted). An explicit verified path IS the certificate.
"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "R1-cpu-count"))
import pc
import cpsat_rho

def factorial(n):
    f = 1
    for i in range(2, n + 1): f *= i
    return f

def unrank(rank, n):
    avail = list(range(n)); s = [0] * n; f = factorial(n - 1)
    for i in range(n):
        q = rank // f; rank %= f; s[i] = avail.pop(q)
        if i < n - 1: f //= (n - 1 - i)
    return s

def adjacency(sigma):
    n = len(sigma)
    pos = [0] * n
    for i, v in enumerate(sigma): pos[v] = i
    adj = [[set() for _ in range(n)] for _ in range(2)]
    for i in range(n - 1):
        adj[0][i].add(i + 1); adj[0][i + 1].add(i)
    for v in range(n - 1):
        a, b = pos[v], pos[v + 1]
        adj[1][a].add(b); adj[1][b].add(a)
    return adj

def find_ham_path(sigma):
    """DFS with reach prune; returns (vertices, colors) covering all n vertices, or None."""
    n = len(sigma)
    adj = adjacency(sigma)
    adjm = [[0] * n for _ in range(2)]
    for c in range(2):
        for v in range(n):
            m = 0
            for w in adj[c][v]: m |= 1 << w
            adjm[c][v] = m
    full = (1 << n) - 1
    def reach(v, c, vis):
        seen = [0, 0]; frontier = [0, 0]; acc = 0
        first = adjm[c][v] & ~vis
        frontier[1 - c] = first; seen[1 - c] = first; acc |= first
        change = True
        while change:
            change = False
            for col in range(2):
                f = frontier[col]; frontier[col] = 0
                while f:
                    w = (f & -f).bit_length() - 1; f &= f - 1
                    nx = adjm[col][w] & ~vis & ~seen[1 - col]
                    if nx:
                        seen[1 - col] |= nx; frontier[1 - col] |= nx; acc |= nx; change = True
        return bin(acc).count("1")
    sys.setrecursionlimit(10000)
    def dfs(v, c, vis, path, cols):
        if len(path) == n: return True
        cand = adjm[c][v] & ~vis
        if not cand: return False
        if len(path) + reach(v, c, vis) < n: return False
        while cand:
            w = (cand & -cand).bit_length() - 1; cand &= cand - 1
            path.append(w); cols.append(c)
            if dfs(w, 1 - c, vis | (1 << w), path, cols): return True
            path.pop(); cols.pop()
        return False
    for sv in range(n):
        for sc in range(2):
            path = [sv]; cols = []
            if dfs(sv, sc, 1 << sv, path, cols):
                return path, cols
    return None

def check_path(sigma, path, cols):
    """Independent structural check: n distinct vertices, consecutive edges exist
    in the claimed color, colors strictly alternate."""
    n = len(sigma)
    if len(path) != n or len(set(path)) != n or len(cols) != n - 1: return False
    adj = adjacency(sigma)
    for i in range(n - 1):
        if path[i + 1] not in adj[cols[i]][path[i]]: return False
        if i > 0 and cols[i] == cols[i - 1]: return False
    return True

def main():
    outdir = sys.argv[1]
    per_side = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    st = json.load(open(os.path.join(outdir, "checkpoint.json")))
    n = st["n"]
    assert st["cursor"] >= st["total"], "run not complete"
    def_ranks = [r for _, r in st["res_def"][:per_side]]
    nd_ranks = [r for _, r in st["res_nondef"][:per_side]]
    print(f"[cert] n={n} deficient sample={len(def_ranks)} non-deficient sample={len(nd_ranks)}")
    ok = True
    for i, rank in enumerate(def_ranks):
        s = unrank(rank, n)
        t0 = time.time(); r_pc = pc.rho(tuple(s)); t_pc = time.time() - t0
        t0 = time.time(); status, r_cp = cpsat_rho.rho_cpsat(list(s), timeout=600.0); t_cp = time.time() - t0
        good = (r_pc <= n - 1) and (status == "OPTIMAL") and (r_cp == r_pc)
        ok &= good
        print(f"DEF {i} rank={rank} perm={','.join(map(str, s))} pc={r_pc}({t_pc:.1f}s) "
              f"cpsat={r_cp}[{status}]({t_cp:.1f}s) {'PASS' if good else 'FAIL'}", flush=True)
    for i, rank in enumerate(nd_ranks):
        s = unrank(rank, n)
        t0 = time.time(); res = find_ham_path(s); t_f = time.time() - t0
        good = res is not None and check_path(s, res[0], res[1])
        ok &= good
        pth = "-".join(map(str, res[0])) if res else "NONE"
        print(f"NONDEF {i} rank={rank} perm={','.join(map(str, s))} path={pth} "
              f"({t_f:.1f}s) {'PASS' if good else 'FAIL'}", flush=True)
    print("CERT SUMMARY:", "ALL PASS" if ok else "SOME FAIL",
          f"(n={n}, {len(def_ranks)}+{len(nd_ranks)} sigma)")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
