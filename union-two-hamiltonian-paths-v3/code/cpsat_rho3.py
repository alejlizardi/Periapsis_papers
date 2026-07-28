# cpsat_rho3.py — INDEPENDENT-FORMALISM verification of rho3(id,s,t): longest
# properly-colored simple path in the union of THREE Hamiltonian paths (word
# convention: colors are the edge sets E(id), E(s), E(t); consecutive edges of
# the path must have distinct colors — arc-based MILP with MTZ, modeled on the
# paper's cpsat_rho.py, extended to 3 colors).  NOT the DFS bitmask engine.
# Usage: py cpsat_rho3.py defects_file  (reads DEF3 lines)  OR
#        py cpsat_rho3.py "s_comma" "t_comma"
import sys, re
from ortools.sat.python import cp_model

def rho3_cpsat(s, t, timeout=600.0):
    n = len(s)
    edges = []  # (u, v, color): 0 = id path, 1 = word s, 2 = word t
    for i in range(n - 1):
        edges.append((i, i + 1, 0))
    for i in range(n - 1):
        edges.append((s[i], s[i + 1], 1))
    for i in range(n - 1):
        edges.append((t[i], t[i + 1], 2))
    arcs = []
    for (u, v, c) in edges:
        arcs.append((u, v, c))
        arcs.append((v, u, c))
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(arcs))]
    y = [m.NewBoolVar(f"y{v}") for v in range(n)]
    st_ = [m.NewBoolVar(f"s{v}") for v in range(n)]
    tt = [m.NewBoolVar(f"t{v}") for v in range(n)]
    u = [m.NewIntVar(0, n - 1, f"u{v}") for v in range(n)]
    inA = [[] for _ in range(n)]
    outA = [[] for _ in range(n)]
    for i, (a, b, c) in enumerate(arcs):
        outA[a].append(i)
        inA[b].append(i)
    m.Add(sum(st_) == 1)
    m.Add(sum(tt) == 1)
    for v in range(n):
        m.Add(sum(x[i] for i in inA[v]) == y[v] - st_[v])
        m.Add(sum(x[i] for i in outA[v]) == y[v] - tt[v])
        m.AddImplication(st_[v], y[v])
        m.AddImplication(tt[v], y[v])
        for i in inA[v]:
            for j in outA[v]:
                if arcs[i][2] == arcs[j][2]:
                    m.Add(x[i] + x[j] <= 1)
    for i, (a, b, c) in enumerate(arcs):
        m.Add(u[b] >= u[a] + 1 - n * (1 - x[i]))
    m.Maximize(sum(y))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    solver.parameters.num_search_workers = 8
    r = solver.Solve(m)
    status = solver.StatusName(r)
    val = int(solver.ObjectiveValue()) if r in (cp_model.OPTIMAL, cp_model.FEASIBLE) else -1
    return status, val

if __name__ == "__main__":
    if len(sys.argv) == 3 and "," in sys.argv[1]:
        s = [int(z) for z in sys.argv[1].split(",")]
        t = [int(z) for z in sys.argv[2].split(",")]
        print(rho3_cpsat(s, t))
        sys.exit(0)
    ok = True
    with open(sys.argv[1]) as f:
        for line in f:
            m_ = re.match(r"DEF3 s=\[([\d,]+)\] t=\[([\d,]+)\] rho3=(\d+)", line)
            if not m_:
                continue
            s = [int(z) for z in m_.group(1).split(",")]
            t = [int(z) for z in m_.group(2).split(",")]
            claimed = int(m_.group(3))
            status, val = rho3_cpsat(s, t)
            tag = "PASS" if (status == "OPTIMAL" and val == claimed) else "FAIL"
            if tag == "FAIL":
                ok = False
            print(f"{tag} s={s} t={t} claimed={claimed} cpsat={status}:{val}", flush=True)
    print("ALL PASS" if ok else "*** MISMATCH ***")
