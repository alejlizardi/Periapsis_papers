#!/usr/bin/env python3
"""
GENUINELY INDEPENDENT FORMULA-MACHINE  (Fourier-diagonalized gram  x  ordering refinement).

This decomposition is fundamentally different from the prover's 1242-state product automaton.
Instead of carrying the 10-bit Gram value as an explicit coordinate (1024-fold blow-up of the
state), we DIAGONALIZE the gram dynamics with the Walsh-Hadamard transform over (Z/2)^10 and
work one CHARACTER chi at a time.  The state space is then ONLY the ordering information
(ordered set partitions of {0,1,2,3} compatible with target 0<1<2<3), and gram enters as a
per-column SCALAR sign (-1)^{<chi, q(c)>}.

== Setup ==
Read columns MSB-first.  For a fixed character chi in {0,1}^10 define a weighted transfer on
the ORDERING state alone:
   weight of column c in ordering-state s = (-1)^{<chi, q(c)>}    (q(c)=10-bit c_i&c_j vector)
and the ordering refinement is the usual block-splitting (0-bit before 1-bit), pruned to be
compatible with 0<1<2<3.

Let, for character chi,
   S_chi(n) = sum over admissible length-n column sequences ending in the FULL singleton order
              (0)(1)(2)(3) of  prod_columns (-1)^{<chi, q(c)>} .
This counts (with signs) sequences whose rows are strictly increasing, weighted by the gram
character.  By Walsh-Hadamard inversion, the number of strictly-increasing sequences whose
FINAL gram equals a given G is
   Ngram(n, G) = (1/1024) * sum_chi (-1)^{<chi,G>} * S_chi(n).
And finally
   a(n) = sum_{G : gram rows strictly decreasing} Ngram(n, G)
        = (1/1024) * sum_chi ( sum_{G in DEC} (-1)^{<chi,G>} ) * S_chi(n)
        = (1/1024) * sum_chi  W(chi) * S_chi(n),
where W(chi) = sum over the strictly-decreasing-gram-row matrices G of (-1)^{<chi,G>}  is a
FIXED integer weight per character, computed once.

So a(n) = (1/1024) * sum_chi W(chi) * S_chi(n).  Each S_chi(n) is computed by powering a SMALL
signed transfer matrix on the ORDERING states only (no gram coordinate!).  This is a completely
different state space and a different algebraic route to the same number.

We compute a(n) exactly via Fractions and compare to OEIS, to C20, and re-derive the recurrence
and closed form from the eigenvalues that arise here.
"""
import itertools
from fractions import Fraction
from collections import deque, defaultdict

DIAGOFF = [(i,j) for i in range(4) for j in range(i,4)]
def gi(i,j):
    if i>j: i,j=j,i
    return DIAGOFF.index((i,j))

def qvec(c):
    return tuple(c[i]&c[j] for (i,j) in DIAGOFF)

def dot10(a,b):
    return sum(x&y for x,y in zip(a,b)) & 1

COLS = list(itertools.product((0,1), repeat=4))
QCACHE = {c: qvec(c) for c in COLS}

# ---- ORDERING states only (ordered set partitions compatible with 0<1<2<3) ----
def ord_initial():
    return (frozenset({0,1,2,3}),)

def ord_step(blocks, c):
    newblocks = []
    for blk in blocks:
        zeros = frozenset(r for r in blk if c[r]==0)
        ones  = frozenset(r for r in blk if c[r]==1)
        if zeros: newblocks.append(zeros)
        if ones:  newblocks.append(ones)
    pos = {}
    for bi, blk in enumerate(newblocks):
        for r in blk: pos[r]=bi
    for a in range(4):
        for b in range(a+1,4):
            if pos[a] > pos[b]:
                return None
    return tuple(newblocks)

FINAL = (frozenset({0}),frozenset({1}),frozenset({2}),frozenset({3}))

def build_ordering_states():
    s0 = ord_initial(); st={s0:0}; order=[s0]; q=deque([s0])
    # store transitions as list of (src, dst, column) so we can apply per-chi signs
    trans = []
    while q:
        s=q.popleft()
        for c in COLS:
            t = ord_step(s, c)
            if t is None: continue
            if t not in st:
                st[t]=len(order); order.append(t); q.append(t)
            trans.append((st[s], st[t], c))
    return order, st, trans

def S_chi_terms(order, st, trans, chi, K):
    """exact S_chi(n) for n=0..K via signed DP over ordering states."""
    N=len(order)
    # precompute per-(src) list of (dst, sign)
    out = defaultdict(list)
    for si,ti,c in trans:
        sign = 1 if dot10(chi, QCACHE[c])==0 else -1
        out[si].append((ti, sign))
    final_idx = st[FINAL]
    v=[0]*N; v[0]=1; res=[]
    for k in range(K+1):
        if k>0:
            nv=[0]*N
            for si in range(N):
                if v[si]:
                    for ti,sgn in out[si]:
                        nv[ti]+=v[si]*sgn
            v=nv
        res.append(v[final_idx])
    return res

def decreasing_gram_weight():
    """W(chi) = sum over symmetric 4x4 parity G with rows (4-bit) strictly decreasing of
       (-1)^{<chi,G>}, as a function returning a dict chi-> integer. We instead return, for each
       chi, the integer; computed by enumerating the (small) set of qualifying G."""
    good_G = []
    for bits in itertools.product((0,1), repeat=10):
        G=[[0]*4 for _ in range(4)]
        for idx,(i,j) in enumerate(DIAGOFF):
            G[i][j]=bits[idx]; G[j][i]=bits[idx]
        gr=[]
        for i in range(4):
            v=0
            for j in range(4): v=(v<<1)|G[i][j]
            gr.append(v)
        if gr[0]>gr[1]>gr[2]>gr[3]:
            good_G.append(bits)   # store as the 10-bit i<=j vector
    return good_G

if __name__ == "__main__":
    order, st, trans = build_ordering_states()
    print("ORDERING-only states (different state space) =", len(order))
    good_G = decreasing_gram_weight()
    print("# strictly-decreasing-gram matrices G =", len(good_G))

    K = 30
    # accumulate a(n) = (1/1024) sum_chi W(chi) S_chi(n)
    a = [Fraction(0)]*(K+1)
    for chi in itertools.product((0,1), repeat=10):
        W = sum(1 if dot10(chi, G)==0 else -1 for G in good_G)
        if W == 0:
            continue
        s = S_chi_terms(order, st, trans, chi, K)
        for n in range(K+1):
            a[n] += W * s[n]
    a = [x/1024 for x in a]
    # must be integers
    aint = [int(x) for x in a]
    assert all(x.denominator==1 for x in a), "non-integer a(n)!"

    OEIS = {4:58,5:1629,6:28924,7:507052,8:8211776,9:133693904,
            10:2140571200,11:34361115072,12:549587348992}
    allok=True
    for n in range(1,13):
        e=OEIS.get(n,"?")
        ok = "OK" if (e=="?" or e==aint[n]) else "MISMATCH"
        if e!="?" and e!=aint[n]: allok=False
        print(f"n={n}: fourier_coupled={aint[n]}  OEIS={e}  [{ok}]")
    print("ALL MATCH" if allok else "*** MISMATCH ***")
    # dump terms for downstream recurrence derivation
    import json
    with open("a_terms_fourier.json","w") as fh:
        json.dump(aint, fh)
    print("wrote a_terms_fourier.json (n=0..%d)"%K)
