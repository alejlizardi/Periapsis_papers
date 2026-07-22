"""verify_floor.py -- mechanical validation of the transit-floor lower bound
(paper Theorem "Floor lower bound", Section "The transit floor").

CLAIM under test: for every pi in S_b, b>=3, the union of two canonical
near-perfect matchings (one inside the P-path, one inside the V-path), chosen
by the case tree below (the proof's designated construction), contains a
simple properly colored terminal-to-terminal path of size >= 3 whose boundary
edge colors are the complements of the endpoint roles' port colors, i.e. a
type-(i) admissible visit  ==>  M(pi) >= 3.

Case tree (the PROOF's designated construction -- any failure falsifies it):
  b ODD:
    D: M_P = positions {(1,2),(3,4),...,(b-2,b-1)}      exposes position 0 = p0
       c  = v1 if v1 != 0 else v0   (c != 0 always; value-index of c is even)
       M_V = value pairs exposing exactly value-index(c)
       => unique path 0 ~> c, first edge V, last edge P, odd size >= 3.
  b EVEN:
    coincidences := |{0,b-1} ∩ {v0,v1}| (as vertices)
    A: M_V = values {(0,1),(2,3),...,(b-2,b-1)} (perfect),
       M_P = positions {(1,2),...,(b-3,b-2)} (exposes {0, b-1})
       => unique path 0 ~> b-1, V boundary both ends, roles (p0,p1), even size;
          size>=4 iff NOT ({pi(0),pi(b-1)} == {2k,2k+1} for some k).
    C: M_P = positions {(1,2),...,(b-3,b-2)} (exposes {0,b-1}),
       M_V = values {(1,2),(3,4),...,(b-3,b-2)} (exposes {v0,v1})
       => components include admissible path(s); wins whenever
          NOT (spine pairing with both spines of size 2).
    designated: 2 coincidences -> A;  1 coincidence -> C;
                0 coincidences -> A if {pi(0),pi(b-1)} not an even-aligned
                consecutive pair {2k,2k+1}, else C.

Validation is FROM def:M directly (independent of the construction logic):
proper coloring, simplicity, endpoint roles exist, boundary colors are the
complements of the endpoint roles' port colors, size >= 3.
"""
import sys, itertools, random
from collections import Counter

P, V = 0, 1

def inv_of(pi):
    inv = [0]*len(pi)
    for i, v in enumerate(pi):
        inv[v] = i
    return inv

def p_edges(b):
    return {frozenset((i, i+1)) for i in range(b-1)}

def v_edges(pi):
    inv = inv_of(pi)
    return {frozenset((inv[v], inv[v+1])) for v in range(len(pi)-1)}

# ---------- canonical matchings ----------
def match_path_expose(order, exposed):
    """Matching inside the path visiting `order` (list of vertices) exposing
    exactly the vertices in `exposed` (subset of order). Greedy left-to-right:
    pair consecutive non-exposed vertices. Returns set of frozenset edges, or
    None if the exposure pattern is infeasible (a gap of odd length)."""
    m = set()
    run = []
    for x in order:
        if x in exposed:
            if len(run) % 2 == 1:
                return None
            for i in range(0, len(run), 2):
                m.add(frozenset((run[i], run[i+1])))
            run = []
        else:
            run.append(x)
    if len(run) % 2 == 1:
        return None
    for i in range(0, len(run), 2):
        m.add(frozenset((run[i], run[i+1])))
    return m

def walk_path(mA, mB, start, colA, colB):
    """From `start` (covered by exactly one of the matchings), alternate
    matchings; return (vertices, colors). Matchings are dicts v->partner."""
    path = [start]
    colors = []
    cur = start
    use_A = cur in mA
    assert (cur in mA) != (cur in mB), "start must be covered by exactly one"
    while True:
        m, col = (mA, colA) if use_A else (mB, colB)
        if cur not in m:
            break
        nxt = m[cur]
        path.append(nxt)
        colors.append(col)
        cur = nxt
        use_A = not use_A
    return path, colors

def as_dict(m):
    d = {}
    for e in m:
        a, b = tuple(e)
        d[a] = b; d[b] = a
    return d

# ---------- def:M validation (independent) ----------
def roles_list(pi):
    b = len(pi); inv = inv_of(pi)
    return [("p0", 0, P), ("p1", b-1, P), ("v0", inv[0], V), ("v1", inv[b-1], V)]

def validate_admissible(pi, path, colors):
    """True iff (path, colors) is a type-(i) admissible visit of size >= 3
    per def:M, checked from scratch."""
    b = len(pi)
    if len(path) < 3:
        return False, "size<3"
    if len(set(path)) != len(path):
        return False, "not simple"
    pe, ve = p_edges(b), v_edges(pi)
    for k in range(len(path)-1):
        e = frozenset((path[k], path[k+1]))
        if colors[k] == P and e not in pe:
            return False, f"edge {k} not a P-edge"
        if colors[k] == V and e not in ve:
            return False, f"edge {k} not a V-edge"
        if k > 0 and colors[k] == colors[k-1]:
            return False, "not properly colored"
    R = roles_list(pi)
    ok_start = [ (nm,pc) for (nm, vtx, pc) in R if vtx == path[0] and colors[0] != pc ]
    ok_end   = [ (nm,pc) for (nm, vtx, pc) in R if vtx == path[-1] and colors[-1] != pc ]
    for sn, _ in ok_start:
        for en, _ in ok_end:
            if sn != en:
                return True, f"{sn}->{en}"
    return False, "no valid role pair"

# ---------- the construction, exactly as the proof designates ----------
def construct(pi):
    """Returns (tag, path, colors) per the proof's designated case tree."""
    b = len(pi); inv = inv_of(pi)
    v0, v1 = inv[0], inv[b-1]
    val_order = [inv[v] for v in range(b)]          # V-path vertex order
    pos_order = list(range(b))                      # P-path vertex order

    if b % 2 == 1:
        c = v1 if v1 != 0 else v0
        mP = match_path_expose(pos_order, {0})
        mV = match_path_expose(val_order, {c})
        assert mP is not None and mV is not None, "odd-b exposure infeasible?!"
        dP, dV = as_dict(mP), as_dict(mV)
        # start at 0: covered by V only
        path, colors = walk_path(dV, dP, 0, V, P)
        return ("D", path, colors)

    # b even
    coinc = {0, b-1} & {v0, v1}
    evenpair = {pi[0], pi[b-1]}
    lo = min(evenpair)
    a_fails = (max(evenpair) - lo == 1) and (lo % 2 == 0)

    def build_A():
        mV = match_path_expose(val_order, set())            # perfect
        mP = match_path_expose(pos_order, {0, b-1})
        dP, dV = as_dict(mP), as_dict(mV)
        path, colors = walk_path(dV, dP, 0, V, P)
        return ("A", path, colors)

    def build_C():
        mP = match_path_expose(pos_order, {0, b-1})
        mV = match_path_expose(val_order, {v0, v1})
        dP, dV = as_dict(mP), as_dict(mV)
        # candidate path starts: every exposed-by-exactly-one vertex
        best = None
        for s in (0, b-1, v0, v1):
            inP, inV = s in dP, s in dV
            if inP == inV:
                continue                                    # deg 0 or 2
            if inV:
                path, colors = walk_path(dV, dP, s, V, P)
            else:
                path, colors = walk_path(dP, dV, s, P, V)
            ok, tag = validate_admissible(pi, path, colors)
            if ok:
                return ("C", path, colors)
            best = (path, colors)
        return ("C-none", *(best if best else ([], [])))

    if len(coinc) == 2:
        return build_A()
    if len(coinc) == 1:
        return build_C()
    return build_C() if a_fails else build_A()

# ---------- exhaustive / sampled verification ----------
def run_exhaustive(b):
    fails = []
    tags = Counter()
    for pi in itertools.permutations(range(b)):
        tag, path, colors = construct(pi)
        ok, info = validate_admissible(pi, path, colors)
        tags[tag] += 1
        if not ok:
            fails.append((pi, tag, path, colors, info))
            if len(fails) <= 5:
                print(f"  FAIL b={b} pi={pi} tag={tag} path={path} colors={colors}: {info}")
    return fails, tags

def run_sampled(b, n, seed=20260703):
    rng = random.Random(seed)
    fails = []
    tags = Counter()
    base = list(range(b))
    for _ in range(n):
        pi = base[:]
        rng.shuffle(pi)
        tag, path, colors = construct(tuple(pi))
        ok, info = validate_admissible(tuple(pi), path, colors)
        tags[tag] += 1
        if not ok:
            fails.append((tuple(pi), tag, info))
    return fails, tags

if __name__ == "__main__":
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    total_fail = 0
    for b in range(3, bmax+1):
        fails, tags = run_exhaustive(b)
        total_fail += len(fails)
        print(f"b={b}: {sum(tags.values())} perms, tags={dict(tags)}, FAILURES={len(fails)}")
    for b, n in ((11, 200000), (12, 200000), (13, 100000), (20, 50000), (31, 20000), (50, 5000)):
        fails, tags = run_sampled(b, n)
        total_fail += len(fails)
        print(f"b={b} sampled n={n}: tags={dict(tags)}, FAILURES={len(fails)}")
    print("TOTAL FAILURES:", total_fail)
