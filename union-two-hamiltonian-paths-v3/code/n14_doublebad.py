# n14_doublebad_hunt.py — exhaustive double-bad hunt at n=14 (m=7).
#
# The unique bad word class at m=7 (word_bruteforce.py) is
#   BAD = 0 1 0 1 2 3 4 3 4 2 5 6 5 6.
# gamma-BF at n=14 (even) holds IFF no block-free sigma has BOTH W(sigma) and
# W(sigma^{-1}) bad (Lemma W + uniqueness of the bad class). W(sigma) bad means
# W(sigma) is a relabeling of BAD; enumerate ALL such sigma:
#   - relabeling pi in S_7: word symbol s -> position pair {2 pi(s), 2 pi(s)+1}
#   - orientation bit per symbol: which vertex its first occurrence uses
# sigma^{-1}(v) = vertex assigned to word position v. Filter block-free
# (|seq[v+1]-seq[v]| != 1 for all v — inverse-invariant). For survivors,
# canonicalize W(sigma^{-1}) and compare against canonical BAD.
from itertools import permutations

BAD = [0, 1, 0, 1, 2, 3, 4, 3, 4, 2, 5, 6, 5, 6]
M = 7
N = 14

def canon(w):
    lab = {}
    out = []
    for s in w:
        if s not in lab:
            lab[s] = len(lab)
        out.append(lab[s])
    return tuple(out)

CANON_BAD = canon(BAD)
# occurrence index (0 = first, 1 = second) per word position
occ = []
seen = {}
for s in BAD:
    occ.append(seen.get(s, 0))
    seen[s] = seen.get(s, 0) + 1

total_real = blockfree_real = doublebad = 0
examples = []
for pi in permutations(range(M)):
    for bits in range(1 << M):
        seq = []
        for t in range(N):
            s = BAD[t]
            first = (bits >> s) & 1
            v = 2 * pi[s] + (first ^ occ[t])
            seq.append(v)
        total_real += 1
        if any(abs(seq[v + 1] - seq[v]) == 1 for v in range(N - 1)):
            continue
        blockfree_real += 1
        # sigma^{-1}(v) = seq[v]; sigma(p): invert
        sigma = [0] * N
        for v, p in enumerate(seq):
            sigma[p] = v
        winv = [sigma[t] // 2 for t in range(N)]  # W(sigma^{-1})
        if canon(winv) == CANON_BAD:
            doublebad += 1
            if len(examples) < 5:
                examples.append(sigma)

print(f"realizations={total_real} blockfree={blockfree_real} double_bad={doublebad}")
for e in examples:
    print("DOUBLE-BAD sigma:", ",".join(map(str, e)))
