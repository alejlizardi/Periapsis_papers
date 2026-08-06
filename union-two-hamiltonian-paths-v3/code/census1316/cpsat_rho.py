# cpsat_rho.py — independent CP-SAT verification of rho(sigma) = longest
# properly-colored (alternating) simple path in H_sigma, arc-based MILP with MTZ.
# Formalism: NOT the bitmask DP, NOT the IC/TX transfer — cross-formalism gate
# for the R4 killing-family chain points (sigma_n^{+k}).
# Usage: py -3 cpsat_rho.py "8,3,6,1,4,7,2,5,0" [timeout_s]
import sys
from ortools.sat.python import cp_model

def rho_cpsat(sigma, timeout=3600.0):
    n = len(sigma)
    inv = [0] * n
    for i, v in enumerate(sigma):
        inv[v] = i
    edges = []  # (u, v, color): P=0 position path, V=1 value path
    for i in range(n - 1):
        edges.append((i, i + 1, 0))
    for v in range(n - 1):
        edges.append((inv[v], inv[v + 1], 1))
    arcs = []
    for (u, v, c) in edges:
        arcs.append((u, v, c))
        arcs.append((v, u, c))
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(arcs))]
    y = [m.NewBoolVar(f"y{v}") for v in range(n)]
    s = [m.NewBoolVar(f"s{v}") for v in range(n)]
    t = [m.NewBoolVar(f"t{v}") for v in range(n)]
    u = [m.NewIntVar(0, n - 1, f"u{v}") for v in range(n)]
    inA = [[] for _ in range(n)]
    outA = [[] for _ in range(n)]
    for i, (a, b, c) in enumerate(arcs):
        outA[a].append(i)
        inA[b].append(i)
    m.Add(sum(s) == 1)
    m.Add(sum(t) == 1)
    for v in range(n):
        m.Add(sum(x[i] for i in inA[v]) == y[v] - s[v])
        m.Add(sum(x[i] for i in outA[v]) == y[v] - t[v])
        m.AddImplication(s[v], y[v])
        m.AddImplication(t[v], y[v])
        # alternation: in-arc and out-arc at v cannot share color
        for i in inA[v]:
            for j in outA[v]:
                if arcs[i][2] == arcs[j][2]:
                    m.Add(x[i] + x[j] <= 1)
    # MTZ subtour elimination on used arcs
    for i, (a, b, c) in enumerate(arcs):
        m.Add(u[b] >= u[a] + 1 - n * (1 - x[i]))
    m.Maximize(sum(y))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    solver.parameters.num_search_workers = 8
    st = solver.Solve(m)
    status = solver.StatusName(st)
    return status, int(solver.ObjectiveValue()) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else -1

if __name__ == "__main__":
    sigma = [int(z) for z in sys.argv[1].split(",")]
    timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0
    status, val = rho_cpsat(sigma, timeout)
    print(f"n={len(sigma)} status={status} rho={val}")
