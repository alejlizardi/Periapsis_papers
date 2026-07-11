"""verify_spectrum.py -- exhaustive transit-number spectrum for small b
(paper Remark "The spectrum of transit numbers" and the base cases b = 4, 5 of
the transit-floor upper bound).

For every b in 2..BMAX (default 8), computes M(pi) for EVERY pi in S_b using
the reference implementation pc.transit_M (an exact constrained search written
directly from the definition of an admissible visit), and prints:
  mu(b) = min M(pi), the attained spectrum {M(pi) : pi in S_b}, and the
  count of minimizers.
Expected (paper): mu(2) = 2; mu(b) = 3 for b >= 3; spectrum {3,...,b} for even
b and {3} u {5,...,b} for odd b (the value 4 unattained at odd b <= BMAX).
"""
import sys, itertools
from collections import Counter
import pc

def spectrum(b):
    cnt = Counter()
    for p in itertools.permutations(range(b)):
        cnt[pc.transit_M(p)] += 1
    return cnt

if __name__ == "__main__":
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for b in range(2, bmax + 1):
        cnt = spectrum(b)
        mu = min(cnt)
        vals = sorted(cnt)
        missing = [v for v in range(mu, b + 1) if v not in cnt]
        print(f"b={b}: mu(b)={mu}  spectrum={vals}  "
              f"unattained_in_[mu,b]={missing}  minimizers={cnt[mu]}")
