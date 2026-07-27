#!/usr/bin/env python3
"""
INDEPENDENT CHECK (not part of the proof). From the EXACT transfer matrix, extract terms
a(0..K) (exact integers) and verify Conjecture 20's closed form against them: a(n) == C20(n)
for all n = 4..K, in exact rational arithmetic. This corroborates the result over a finite
range; the rigorous minimal-order argument (det H_11, shifted Hankel solve, residual
propagation) lives in proof_matrix.py. For reference, the closed form's characteristic
polynomial has roots with multiplicities {16:1, 8:2, -8:1, 4:3, -4:2, 2:2} (degree 11).

Exact rational arithmetic via fractions.Fraction. No floats anywhere.
"""
import itertools
from fractions import Fraction
from collections import defaultdict, deque

COLS = list(itertools.product((0,1), repeat=4))
PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
DIAGOFF = [(i,j) for i in range(4) for j in range(i,4)]
def gram_index(i,j):
    if i>j: i,j=j,i
    return DIAGOFF.index((i,j))
def initial_state():
    return (tuple(sorted({p:'=' for p in PAIRS}.items())), tuple([0]*10))
def step(state, c):
    rel = dict(state[0]); g = list(state[1])
    for (i,j) in DIAGOFF:
        if c[i] & c[j]: g[gram_index(i,j)] ^= 1
    for (i,j) in PAIRS:
        if rel[(i,j)] == '=' and c[i]!=c[j]:
            rel[(i,j)] = '<' if c[i]<c[j] else '>'
        if rel[(i,j)] == '>': return None
    return (tuple(sorted(rel.items())), tuple(g))
def accepting(state):
    rel = dict(state[0]); g = state[1]
    if any(rel[p]!='<' for p in PAIRS): return False
    G=[[0]*4 for _ in range(4)]
    for (i,j) in DIAGOFF:
        G[i][j]=g[gram_index(i,j)]; G[j][i]=g[gram_index(i,j)]
    gr=[]
    for i in range(4):
        v=0
        for j in range(4): v=(v<<1)|G[i][j]
        gr.append(v)
    return gr[0]>gr[1]>gr[2]>gr[3]

def build():
    start=initial_state(); states={start:0}; order=[start]; q=deque([start]); edges={}
    while q:
        s=q.popleft(); si=states[s]
        for c in COLS:
            t=step(s,c)
            if t is None: continue
            if t not in states: states[t]=len(order); order.append(t); q.append(t)
            ti=states[t]; edges[(si,ti)]=edges.get((si,ti),0)+1
    return order,states,edges

def terms(order,states,edges,K):
    N=len(order); out=defaultdict(list)
    for (si,ti),m in edges.items(): out[si].append((ti,m))
    acc=[idx for s,idx in states.items() if accepting(s)]
    v=[0]*N; v[0]=1; res=[]
    for k in range(0,K+1):
        if k>0:
            nv=[0]*N
            for si in range(N):
                if v[si]:
                    for ti,m in out[si]: nv[ti]+=v[si]*m
            v=nv
        res.append(sum(v[i] for i in acc))
    return res  # res[k] = a(k), k=0..K

# ---------- Conjecture 20 closed form, exact rational ----------
def conj20(n):
    n=Fraction(n)
    sgn = Fraction(-1)**int(n) if int(n)%2 else Fraction(1)  # (-1)^n
    def p2(e):  # 2^e, e can be negative integer -> exact Fraction
        e=int(e)
        return Fraction(2)**e if e>=0 else Fraction(1, 2**(-e))
    a = (Fraction(1,3)*p2(2*n-11)*(6*n*n-219*n+820)
         - Fraction(1,9)*p2(n-5)*(3*n+32)
         - Fraction(113,3)*sgn*p2(3*n-14)
         + p2(4*n-9)
         - Fraction(1,3)*sgn*p2(2*n-11)*(13*n-164)
         + Fraction(1,9)*p2(3*n-14)*(288*n-3473))
    return a

if __name__=="__main__":
    order,states,edges=build()
    K=60
    a=terms(order,states,edges,K)
    print("a(0..12) =", a[:13])

    # (A) Conjecture 20 vs exact terms, n>=4
    print("\n(A) Conjecture 20 vs transfer-matrix terms (exact), n=4..%d:"%K)
    allok=True
    for n in range(4,K+1):
        c=conj20(n)
        ok = (c==a[n])
        allok &= ok
        if not ok or n<=6 or n>=K-1:
            print(f"  n={n}: a={a[n]}  conj={c}  {'OK' if ok else 'MISMATCH'}")
    print("  ALL n=4..%d match:"%K, allok)

    # also report conj20 at n=1,2,3 (formula is only claimed n>=4)
    print("\n  formula at n=1,2,3 (NOT claimed; a=0 there):")
    for n in (1,2,3):
        print(f"   n={n}: conj={conj20(n)}  a={a[n]}")
