# word_bruteforce.py — test the Word Conjecture (WC) behind gamma-BF (n even).
#
# WC: every double-occurrence word W of length 2m (m symbols, each exactly twice),
# with no immediate repeat (w[t] != w[t+1]), admits a set S of m-1 steps
# e_t = {w[t], w[t+1]} (t = 0..2m-2), pairwise non-adjacent (no t, t+1 both in S),
# whose edge set is a spanning tree of the m symbols (equivalently a Hamiltonian
# path on symbols; degree <= 2 is automatic since each symbol has 2 occurrences).
#
# Reduction (this run's Lemma W): for n = 2m, sigma with no value-adjacent position
# pair (2i, 2i+1), sigma has a PC-Ham path with P-majority IFF W(sigma) is good,
# where W(sigma)[t] = floor(sigma^{-1}(t) / 2).
#
# usage: py -3 word_bruteforce.py <m> [maxbad]
import sys
from itertools import combinations

def words(m):
    # all double-occurrence words up to symbol relabeling: fix first occurrence order
    # (canonical: symbols appear in first-occurrence order 0,1,2,...)
    def rec(seq, remaining, used):
        if not remaining and used == m:
            yield tuple(seq)
            return
        # place next: any symbol with one copy left (remaining) or a fresh symbol
        last = seq[-1] if seq else -1
        for s in set(remaining):
            if s != last:
                r2 = list(remaining); r2.remove(s)
                yield from rec(seq + [s], r2, used)
        if used < m and last != used:
            yield from rec(seq + [used], remaining + [used], used + 1)
    yield from rec([], [], 0)

def good(w, m):
    steps = [(w[t], w[t + 1]) for t in range(len(w) - 1)]
    T = len(steps)
    for S in combinations(range(T), m - 1):
        if any(S[i + 1] - S[i] == 1 for i in range(len(S) - 1)):
            continue
        # spanning tree check on m symbols with edges steps[t], t in S
        parent = list(range(m))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        ok = True
        for t in S:
            a, b = find(steps[t][0]), find(steps[t][1])
            if a == b:
                ok = False
                break
            parent[a] = b
        if ok:
            return True  # m-1 edges, acyclic => spanning tree
    return False

def main():
    m = int(sys.argv[1])
    maxbad = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    total = 0
    bad = []
    for w in words(m):
        total += 1
        if not good(w, m):
            bad.append(w)
            if len(bad) >= maxbad:
                break
    print(f"m={m} words_tested={total} bad={len(bad)}")
    for w in bad:
        print("BAD:", "".join(str(x) for x in w))

if __name__ == "__main__":
    main()
