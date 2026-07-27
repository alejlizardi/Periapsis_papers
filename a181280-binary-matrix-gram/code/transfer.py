#!/usr/bin/env python3
"""
A181280 EXACT transfer matrix (n-independent automaton), the deep-pass engine.

State while reading columns left to right, column 1 = most significant:
  - ord:  the comparison status of the 4 rows as PREFIX bit strings so far.
          Represented as a tuple `rel` of the 6 pairs (i<j) in
          {'<','>','='}: rel[(i,j)] = '<' if prefix(row_i) < prefix(row_j) decided,
          '>' if decided the other way, '=' if still tied (equal prefixes).
          ('<'/'>' are FROZEN once set; '=' may later become '<' or '>'.)
  - g:    the 4x4 Gram parity matrix accumulated so far (G_{ij} = sum_k m_ik m_jk mod 2),
          stored as a 10-tuple (4 diag + 6 offdiag, i<=j).

Transition on reading a column c = (c0,c1,c2,c3) in {0,1}^4:
  - update Gram: G_{ij} ^= c_i & c_j  for all i<=j.
  - update ord: for each tied pair (i,j) with rel='=', if c_i != c_j then this column
    DECIDES the pair: rel becomes '<' if c_i < c_j (i.e. c_i=0,c_j=1), else '>'.
    Already-decided pairs are unchanged.
  - PRUNE: if any pair (i<j) is decided '>' (row_i > row_j) it can never give the
    required strict-increasing order row1<row2<row3<row4 -> dead, drop it.
    (We require rel[(i,j)] in {'<','='} for all i<j throughout; a '>' is fatal.)

Acceptance after n columns:
  - rows strictly increasing  <=>  rel[(i,j)] == '<' for ALL i<j (no ties left, all correct).
  - Gram rows strictly decreasing: build G (symmetric) from g, read each row i as 4-bit
    int (col 0 = MSB), require g0 > g1 > g2 > g3.

count(n) = number of length-n column sequences ending in an accepting state.

We build the reachable state graph once, then do a length-n DP (matrix power in spirit).
Verified against brute.py for n<=6, then extended.
"""
import itertools

COLS = list(itertools.product((0,1), repeat=4))  # 16 columns
PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
DIAGOFF = [(i,j) for i in range(4) for j in range(i,4)]  # 10 entries, i<=j

def gram_index(i,j):
    if i>j: i,j=j,i
    return DIAGOFF.index((i,j))

def initial_state():
    rel = {p:'=' for p in PAIRS}        # all rows tied at empty prefix
    g = tuple([0]*10)                   # all parities 0
    return (tuple(sorted(rel.items())), g)

def step(state, c):
    rel = dict(state[0]); g = list(state[1])
    # gram update
    for (i,j) in DIAGOFF:
        if c[i] & c[j]:
            g[gram_index(i,j)] ^= 1
    # ord update + prune
    for (i,j) in PAIRS:
        if rel[(i,j)] == '=':
            if c[i] != c[j]:
                rel[(i,j)] = '<' if c[i] < c[j] else '>'
        if rel[(i,j)] == '>':
            return None  # dead: row_i > row_j can't yield strictly increasing rows
    return (tuple(sorted(rel.items())), tuple(g))

def accepting(state):
    rel = dict(state[0]); g = state[1]
    # all pairs strictly '<'
    if any(rel[p] != '<' for p in PAIRS):
        return False
    # build symmetric G, read rows
    G = [[0]*4 for _ in range(4)]
    for (i,j) in DIAGOFF:
        G[i][j] = g[gram_index(i,j)]
        G[j][i] = g[gram_index(i,j)]
    grows = []
    for i in range(4):
        v = 0
        for j in range(4):
            v = (v<<1)|G[i][j]
        grows.append(v)
    return grows[0] > grows[1] > grows[2] > grows[3]

# Build reachable states and transition counts
def build():
    from collections import deque
    start = initial_state()
    states = {start: 0}
    order = [start]
    q = deque([start])
    trans = []  # list of (src_idx, dst_idx) per column, aggregated as counts
    edges = {}  # (src_idx, dst_idx) -> multiplicity
    while q:
        s = q.popleft(); si = states[s]
        for c in COLS:
            t = step(s, c)
            if t is None: continue
            if t not in states:
                states[t] = len(order); order.append(t); q.append(t)
            ti = states[t]
            edges[(si,ti)] = edges.get((si,ti),0)+1
    return order, states, edges

def count_via_dp(order, states, edges, n):
    # vector v[state] = number of length-k column sequences reaching state (from start, after k cols)
    N = len(order)
    # adjacency multiplicity matrix as dict
    from collections import defaultdict
    out = defaultdict(list)
    for (si,ti),m in edges.items():
        out[si].append((ti,m))
    start = 0
    v = [0]*N; v[start] = 1
    # after reading k columns
    res = {}
    for k in range(0, n+1):
        if k>0:
            nv = [0]*N
            for si in range(N):
                if v[si]:
                    for (ti,m) in out[si]:
                        nv[ti] += v[si]*m
            v = nv
        # count accepting at length k
        tot = 0
        for s,idx in states.items():
            if v[idx] and accepting(s):
                tot += v[idx]
        res[k] = tot
    return res

if __name__ == "__main__":
    order, states, edges = build()
    print(f"reachable states: {len(order)}, edges(distinct src-dst): {len(edges)}")
    res = count_via_dp(order, states, edges, 12)
    OEIS = {1:0,2:0,3:0,4:58,5:1629,6:28924,7:507052,8:8211776,9:133693904,
            10:2140571200,11:34361115072,12:549587348992}
    for n in range(1,13):
        e = OEIS.get(n,"?")
        ok = "OK" if e==res[n] else "MISMATCH"
        print(f"n={n}: tm={res[n]}  OEIS={e}  [{ok}]")
