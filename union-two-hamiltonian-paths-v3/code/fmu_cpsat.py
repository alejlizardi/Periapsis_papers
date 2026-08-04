# fmu_cpsat.py — INDEPENDENT-FORMALISM exact max-mu for G(n,sigma): arc-based
# MILP with MTZ (modeled on the paper's cpsat_rho.py), maximizing the number of
# MATCHING edges used by a simple path (no alternation constraint — same-side
# runs are legal).  Cross-check for fmu.exe's bitmask DP and the scale probe
# beyond 2n = 24.
# Usage: py fmu_cpsat.py "sigma_comma" [timeout_s]
#        py fmu_cpsat.py --family  (runs the WP-J probe battery)
import sys
from ortools.sat.python import cp_model

def maxmu_cpsat(sigma, timeout=1200.0, workers=16):
    n = len(sigma)
    N = 2 * n                       # 0..n-1 top v_i ; n..2n-1 bottom u_j
    edges = []                      # (a, b, ismatch)
    for i in range(n - 1):
        edges.append((i, i + 1, 0))
        edges.append((n + i, n + i + 1, 0))
    for i in range(n):
        edges.append((i, n + sigma[i], 1))
    arcs = []
    for (a, b, mt) in edges:
        arcs.append((a, b, mt))
        arcs.append((b, a, mt))
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(arcs))]
    y = [m.NewBoolVar(f"y{v}") for v in range(N)]
    s = [m.NewBoolVar(f"s{v}") for v in range(N)]
    t = [m.NewBoolVar(f"t{v}") for v in range(N)]
    u = [m.NewIntVar(0, N - 1, f"u{v}") for v in range(N)]
    inA = [[] for _ in range(N)]
    outA = [[] for _ in range(N)]
    for i, (a, b, mt) in enumerate(arcs):
        outA[a].append(i)
        inA[b].append(i)
    m.Add(sum(s) == 1)
    m.Add(sum(t) == 1)
    for v in range(N):
        m.Add(sum(x[i] for i in inA[v]) == y[v] - s[v])
        m.Add(sum(x[i] for i in outA[v]) == y[v] - t[v])
        m.AddImplication(s[v], y[v])
        m.AddImplication(t[v], y[v])
    for i, (a, b, mt) in enumerate(arcs):
        m.Add(u[b] >= u[a] + 1 - N * (1 - x[i]))
    m.Maximize(sum(x[i] for i, (a, b, mt) in enumerate(arcs) if mt))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    solver.parameters.num_search_workers = workers
    r = solver.Solve(m)
    status = solver.StatusName(r)
    val = int(solver.ObjectiveValue()) if r in (cp_model.OPTIMAL, cp_model.FEASIBLE) else -1
    bound = int(solver.BestObjectiveBound())
    return status, val, bound

def sigma_ab(a, b):
    out = []
    for j in range(a):
        blk = list(range(j * b, (j + 1) * b))
        blk[1], blk[b - 1] = j * b + b - 1, j * b + 1
        out.extend(blk)
    return out

def sigma_C(comp):
    out = []
    P = 0
    for c in comp:
        out.extend(range(P + c - 1, P - 1, -1))
        P += c
    return out

if __name__ == "__main__":
    if sys.argv[1] == "--family":
        probes = [
            ("sigma_{2,6} n=12 (DP=12)", sigma_ab(2, 6)),
            ("sigma_{2,8} n=16 (DP=16)", sigma_ab(2, 8)),
            ("chain(3,1,3,1,3) n=11 (DP=10)", sigma_C([3, 1, 3, 1, 3])),
            ("chain(5,1,5) n=11 (DP=11)", sigma_C([5, 1, 5])),
            ("sigma_{3,6} n=18", sigma_ab(3, 6)),
            ("chain(5,1,5,1,5) n=17", sigma_C([5, 1, 5, 1, 5])),
            ("chain(9,1,9) n=19", sigma_C([9, 1, 9])),
            ("sigma_{4,6} n=24", sigma_ab(4, 6)),
            ("chain(7,1,7,1,7) n=23", sigma_C([7, 1, 7, 1, 7])),
            ("sigma_{3,8} n=24", sigma_ab(3, 8)),
            ("sigma_{6,6} n=36", sigma_ab(6, 6)),
        ]
        for name, sig in probes:
            status, val, bound = maxmu_cpsat(sig, timeout=float(sys.argv[2]) if len(sys.argv) > 2 else 1200.0)
            n = len(sig)
            print(f"{name:34s} maxmu {status}: {val} (bound {bound})  n={n}  n-1={n-1}", flush=True)
    else:
        sig = [int(z) for z in sys.argv[1].split(",")]
        timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 1200.0
        print(maxmu_cpsat(sig, timeout))
