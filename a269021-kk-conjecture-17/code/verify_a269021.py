#!/usr/bin/env python3
"""
verify_a269021.py -- STANDALONE reproducibility verifier for the paper
  "Permutations of [2n] with an increasing subsequence of length n:
   a proof of the Kauers-Koutschan recurrence for OEIS A269021"

ONE COMMAND, NO EXTERNAL DATA FILES. Depends only on the Python 3 standard library.

What it checks (all arithmetic exact, Fractions/ints; no floats anywhere):

 (1) BRUTE FORCE. a(n) = #{s in S_{2n} : LIS(s) >= n} by enumerating ALL (2n)!
     permutations for n = 1,2,3 (LIS via patience sorting).
 (2) HOOK LENGTHS. a(n) = sum of (f^lambda)^2 over partitions lambda of 2n with
     lambda_1 >= n, f^lambda by the hook length formula, for n <= 10.
 (3) THE REPRESENTATION (Proposition 3 of the paper):
         a(n) = (2n)!^2 * sum_{rho=n}^{2n} (-1)^(rho-n) C(2rho-2,rho-n) / ((rho!)^2 (2n-rho)!)
     matches (1), (2), and the OEIS-listed terms a(0..15) (used only as anchors).
 (4) THE RECURRENCE ON DATA. The conjectured operator L (coefficients embedded
     verbatim from the OEIS entry) annihilates the representation values for n = 0..52.
 (5) THE WZ CERTIFICATE, identity (D) ["diamond"]: with the embedded certificate
     Rt(n,rho) = Xhat(n,rho) / (4(n+3)(n+4)^2(2n+7)^2 (rho-n-1)(rho-n-2)(rho-n-3)),
         sum_{i=0}^4 a_i(n) sigma_i(n,rho) = Rt(n,rho+1) rt(n,rho) - Rt(n,rho)
     as an identity of rational functions in QQ(n,rho). Verified DETERMINISTICALLY:
     the identity is cleared to a polynomial whose degree bound (BN, BR) is computed
     structurally in-script, and evaluated on a (BN+1) x (BR+1) integer grid off all
     denominator zeros -- a nonzero polynomial of bidegree <= (BN,BR) cannot vanish there.
 (6) THE JUNCTION IDENTITY (C) ["club"]:
         sum_{c=0}^{3} sum_i a_i(n) t(n+i,n+c)/t(n+4,n+4) = Rt(n,n+4)
     in QQ(n), verified the same deterministic way (plus a factorial-value guard).
 (7) TELESCOPING, EXHAUSTIVELY: for n = 5..40 and EVERY integer rho, the telescoped
     equation holds and sum_rho (L t)(n,rho) = 0; combined with (4) this instantiates
     the theorem's conclusion on a long range.
 (8) RATIO GUARDS: the closed forms for sigma_i and rt are checked against direct
     factorial ratios of t at many lattice points (so the hand algebra is not trusted).
 (9) PAPER-FORM EQUIVALENCE: the recurrence stated in arXiv:2303.02793 (Conjecture 17,
     in atilde_n = a_n/(2n)!^2) equals the OEIS form times 4(n+3)(n+4)^2(2n+7)^2.

Exit code 0 iff every check passes.  Runtime ~10 seconds.
"""
from fractions import Fraction
from math import comb, factorial
from bisect import bisect_left
from itertools import permutations
import sys

F = Fraction
PASS = []

def check(name, ok):
    PASS.append(bool(ok))
    print(("PASS  " if ok else "FAIL  ") + name)
    if not ok:
        sys.exit(1)

# ============================================================================
# The conjectured operator L: coefficients c_0..c_4 of a(n),...,a(n+4), verbatim
# from the OEIS A269021 entry ("Conjectured recurrence", Kauers & Koutschan,
# Mar 01 2023), with ^ replaced by **.  They are USED, never assumed: the point
# of the paper is that these exact polynomials admit a WZ certificate.
# ============================================================================
L_SRC = {
 0: "-64*(1 + n)**2*(2 + n)**2*(3 + n)*(1 + 2*n)**2*(3 + 2*n)**2*(5 + 2*n)**2*(549760 + 3266000*n + 7264534*n**2 + 8663374*n**3 + 6333869*n**4 + 3012795*n**5 + 952323*n**6 + 198469*n**7 + 26156*n**8 + 1968*n**9 + 64*n**10)",
 1: "16*(2 + n)**2*(3 + n)*(3 + 2*n)**2*(5 + 2*n)**2*(-5543040 - 487964*n + 78563984*n**2 + 229526554*n**3 + 325846005*n**4 + 284054698*n**5 + 165789363*n**6 + 67385886*n**7 + 19359535*n**8 + 3917758*n**9 + 545913*n**10 + 49788*n**11 + 2672*n**12 + 64*n**13)",
 2: "-4*(3 + n)*(5 + 2*n)**2*(-250467360 - 946158512*n - 562569136*n**2 + 3271711596*n**3 + 9313604242*n**4 + 12741784568*n**5 + 11118771121*n**6 + 6753094929*n**7 + 2966908118*n**8 + 959068672*n**9 + 228837227*n**10 + 39928763*n**11 + 4969164*n**12 + 419248*n**13 + 21568*n**14 + 512*n**15)",
 3: "2*(-1300242000 - 6360099840*n - 10812498240*n**2 - 778719870*n**3 + 27672902302*n**4 + 55047935941*n**5 + 60004107039*n**6 + 43766004538*n**7 + 22820793074*n**8 + 8747566435*n**9 + 2488583381*n**10 + 523718876*n**11 + 80260596*n**12 + 8677944*n**13 + 624800*n**14 + 26752*n**15 + 512*n**16)",
 4: "-3*(4 + n)*(8 + 3*n)*(10 + 3*n)*(-15900 - 61278*n - 68928*n**2 + 48699*n**3 + 204716*n**4 + 233810*n**5 + 143536*n**6 + 52389*n**7 + 11324*n**8 + 1328*n**9 + 64*n**10)",
}
L_CODE = {k: compile(v, "<L%d>" % k, "eval") for k, v in L_SRC.items()}
def cL(k, nv):
    return eval(L_CODE[k], {"n": F(nv)})
DEG_C = 21   # max total degree of the c_i in n (c_0: 2+2+1+2+2+2+10 = 21)

# ============================================================================
# The certificate numerator Xhat(n,rho): deg_rho = 11; coefficient of rho^k is a
# polynomial in n, listed ascending. Discovered by Maxima's zeilberger package;
# its correctness is NOT assumed -- it is (re)proved by checks (5)-(7) below.
# ============================================================================
XHAT = {
 0: [0], 1: [0],
 2: [-6015705600, -71744794432, -340218249872, -823169213800, -961021226020, 124265850242, 2438332391222, 4641328365191, 5309568469593, 4316280811201, 2643766419371, 1253332303153, 465624508027, 136031593893, 31134296511, 5520308936, 742809664, 73276736, 4994048, 209920, 4096],
 3: [76008398720, 608225620912, 2228174690440, 4896077714668, 7111706265622, 7048045048890, 4618545883056, 1588597165714, -343809919992, -869853751120, -643733673170, -308244143254, -107309318390, -28100494314, -5565832046, -823951480, -88497984, -6518016, -294400, -6144],
 4: [-40211798400, -176056470320, -189599110100, 437057345050, 1772677634896, 2995261489113, 3224697042433, 2468505910419, 1410599812940, 617445537593, 209962545312, 55838661747, 11624213792, 1884882502, 235063559, 21974184, 1463360, 62336, 1280],
 5: [-24646054240, -192080318780, -668277924786, -1380026212944, -1892614071868, -1822227675986, -1263993467979, -632786469622, -221171415365, -47180065738, -1696243173, 2827670532, 1152539831, 253497006, 35743416, 3233440, 171904, 4096],
 6: [12284419200, 66866806616, 156734550672, 202791566126, 144452346824, 27277926244, -56783535899, -72523162806, -48668124028, -22222515484, -7365118014, -1802389971, -323638993, -41520387, -3603444, -189392, -4544],
 7: [-152434080, 6693806184, 37615945898, 91771112674, 133495099583, 130461693542, 91005728125, 46880034025, 18167548913, 5335822035, 1184513340, 195994162, 23468539, 1923572, 96592, 2240],
 8: [-809164800, -6554503176, -20590657624, -36045877506, -40532796736, -31444382328, -17501663547, -7135652705, -2148067956, -475949216, -76462053, -8635557, -647660, -28848, -576],
 9: [183024640, 1256090516, 3403906826, 5143106268, 4959254786, 3260767977, 1511965679, 501187042, 118501244, 19613637, 2180637, 150284, 5424, 64],
 10: [-16492800, -104027360, -253862020, -339811094, -285313184, -160056409, -61710435, -16429623, -2967839, -346756, -23568, -704],
 11: [549760, 3266000, 7264534, 8663374, 6333869, 3012795, 952323, 198469, 26156, 1968, 64],
}
DEG_XHAT_N = max(len(v) - 1 for v in XHAT.values())   # = 20
DEG_XHAT_R = 11

def Xhat(nv, rv):
    nf, rf = F(nv), F(rv)
    tot = F(0)
    for k, coeffs in XHAT.items():
        s = F(0)
        for c in reversed(coeffs):
            s = s * nf + c
        tot += s * rf ** k
    return tot

def dpoly(nv):
    return F(4) * (nv + 3) * (nv + 4) ** 2 * (2 * nv + 7) ** 2

def Rt(nv, rv):
    nf, rf = F(nv), F(rv)
    return Xhat(nf, rf) / (dpoly(nf) * (rf - nf - 1) * (rf - nf - 2) * (rf - nf - 3))

# ============================================================================
# The summand t(n,rho) and the ratio closed forms
# ============================================================================
def comb_gen(a, b):
    """C(a,b) for any integer a, integer b (0 if b<0)."""
    if b < 0:
        return 0
    num = 1
    for i in range(b):
        num *= (a - i)
    return F(num, factorial(b))

def t_term(nv, rv):
    """t(n,rho); vanishes off n <= rho <= 2n."""
    if rv < nv or rv > 2 * nv:
        return F(0)
    return F((-1) ** (rv - nv)) * comb_gen(2 * rv - 2, rv - nv) * \
        F(factorial(2 * nv) ** 2, factorial(rv) ** 2 * factorial(2 * nv - rv))

def a_rep(nv):
    """a(n) from the representation (Proposition 3)."""
    v = sum(t_term(nv, rv) for rv in range(nv, 2 * nv + 1)) if nv > 0 else F(1)
    assert v.denominator == 1
    return int(v)

def sigma(i, nv, rv):
    """t(n+i,rho)/t(n+4,rho) as the closed rational form (paper, Section 3)."""
    nf, rf = F(nv), F(rv)
    e = F(1)
    for k in range(i, 4):
        e *= -(rf + nf + k - 1) * (2 * nf + 2 * k + 2 - rf) * (2 * nf + 2 * k + 1 - rf)
        e /= (2 * nf + 2 * k + 1) ** 2 * (2 * nf + 2 * k + 2) ** 2 * (rf - nf - k)
    return e

def rt_ratio(nv, rv):
    """t(n+4,rho+1)/t(n+4,rho) as the closed rational form."""
    nf, rf = F(nv), F(rv)
    return -(2 * rf) * (2 * rf - 1) * (2 * nf + 8 - rf) / ((rf - nf - 3) * (rf + nf + 3) * (rf + 1) ** 2)

# ============================================================================
# (1) brute force
# ============================================================================
def lis_len(perm):
    tails = []
    for x in perm:
        j = bisect_left(tails, x)
        if j == len(tails):
            tails.append(x)
        else:
            tails[j] = x
    return len(tails)

def a_brute(nv):
    return sum(1 for p in permutations(range(2 * nv)) if lis_len(p) >= nv)

for nv in (1, 2, 3):
    check("brute force over all (2*%d)! permutations equals representation" % nv,
          a_brute(nv) == a_rep(nv))

# ============================================================================
# (2) hook lengths
# ============================================================================
def partitions_min_first(m, minfirst):
    """partitions of m with largest part >= minfirst."""
    def gen(rem, maxpart, acc):
        if rem == 0:
            yield tuple(acc)
            return
        for p in range(min(maxpart, rem), 0, -1):
            acc.append(p)
            yield from gen(rem - p, p, acc)
            acc.pop()
    for first in range(minfirst, m + 1):
        yield from gen(m - first, first, [first])

def f_lambda(lam):
    m = sum(lam)
    conj = [0] * lam[0]
    for r in lam:
        for c in range(r):
            conj[c] += 1
    den = 1
    for i, r in enumerate(lam):
        for c in range(r):
            den *= (r - c) + (conj[c] - i) - 1
    v = F(factorial(m), den)
    assert v.denominator == 1
    return int(v)

ok = all(sum(f_lambda(l) ** 2 for l in partitions_min_first(2 * nv, nv)) == a_rep(nv)
         for nv in range(1, 11))
check("hook-length sum over {lambda |- 2n, lambda_1 >= n} equals representation, n=1..10", ok)

# ============================================================================
# (3) OEIS anchors (data listed on the entry; anchors only, not proof input)
# ============================================================================
OEIS_TERMS = [1, 2, 23, 588, 24553, 1438112, 108469917, 9996042284, 1086997811325,
              136102249609224, 19269396089593156, 3042212958893941456,
              529708789768374664407, 100813134967124531098768,
              20816198414187782633783462, 4634136282168760818748363080]
check("representation reproduces the 16 OEIS-listed terms a(0..15)",
      all(a_rep(i) == OEIS_TERMS[i] for i in range(16)))

# ============================================================================
# (4) the recurrence annihilates the representation values
# ============================================================================
AVAL = {nv: a_rep(nv) for nv in range(0, 57)}
ok = all(sum(cL(k, nv) * AVAL[nv + k] for k in range(5)) == 0 for nv in range(0, 53))
check("L (OEIS conjectured operator) annihilates representation values, n=0..52", ok)

# ============================================================================
# (8) ratio guards (before using the closed forms in the identities)
# ============================================================================
ok = True
for nv in (6, 9, 14):
    for rv in range(nv + 4, 2 * nv + 4):        # interior of supp t(n+4,.)
        base = t_term(nv + 4, rv)
        if base == 0:
            continue
        for i in range(5):
            ok &= (sigma(i, nv, rv) == t_term(nv + i, rv) / base)
        ok &= (rt_ratio(nv, rv) == t_term(nv + 4, rv + 1) / base)
check("closed ratio forms sigma_i, rt match direct factorial ratios of t", ok)

# ============================================================================
# generic deterministic rational-identity check by degree-bounded grids
# ============================================================================
def bound_from_terms(terms, var):
    """terms: list of (evaluator, num_deg(n,r), den_deg(n,r)). Bound for the
    cleared numerator N = sum_t num_t * prod_{u != t} den_u in variable var."""
    idx = 0 if var == "n" else 1
    Dall = sum(t[2][idx] for t in terms)
    return max(t[1][idx] + Dall - t[2][idx] for t in terms)

# ---- (5) identity (D): sum_i a_i sigma_i - Rt(rho+1) rt + Rt = 0 -----------
# structural degrees (num,den) per addend, derived from the factored forms:
#   a_i sigma_i : num (21 + 3(4-i), 3(4-i)),  den (5(4-i), (4-i))
#   Rt(rho+1)*rt: num (20+1, 11+3),           den (8+2, 3+4)
#   Rt          : num (20, 11),               den (8, 3)
terms_D = []
for i in range(5):
    terms_D.append((
        (lambda i=i: lambda nv, rv: cL(i, nv) * sigma(i, nv, rv))(),
        (21 + 3 * (4 - i), 3 * (4 - i)),
        (5 * (4 - i), (4 - i)),
    ))
terms_D.append((lambda nv, rv: -Rt(nv, rv + 1) * rt_ratio(nv, rv), (21, 14), (10, 7)))
terms_D.append((lambda nv, rv: Rt(nv, rv), (20, 11), (8, 3)))
BN = bound_from_terms(terms_D, "n")
BR = bound_from_terms(terms_D, "r")
N0, R0 = 500, 5000       # off every denominator zero line (rho = n+k, 2n+k, -n-k)
ok = True
for nv in range(N0, N0 + BN + 1):
    for rv in range(R0, R0 + BR + 1):
        s = F(0)
        for ev, _, _ in terms_D:
            s += ev(nv, rv)
        if s != 0:
            ok = False
            break
    if not ok:
        break
check("(D) WZ identity on a %dx%d grid (degree bounds BN=%d, BR=%d) -> holds in QQ(n,rho)"
      % (BN + 1, BR + 1, BN, BR), ok)

# ---- (6) identity (C): junction --------------------------------------------
def ratio_t_junc(i, c, nv):
    """t(n+i, n+c)/t(n+4, n+4) as the closed rational form (fixed offsets)."""
    if c < i:
        return F(0)
    nf = F(nv)
    e = F((-1) ** (c - i))
    for s in range(c - i):
        e *= (2 * nf + 2 * c - 2 - s)
    e /= factorial(c - i)
    d = F(1)
    for s in range(2 * i + 1, 9):
        d *= (2 * nf + s)
    e /= d ** 2
    p = F(1)
    for s in range(c + 1, 5):
        p *= (nf + s)
    e *= p ** 2
    q = F(1)
    for s in range(2 * i - c + 1, 5):
        q *= (nf + s)
    return e * q

ok = True
for nv in (5, 7, 12, 23, 40):
    base = t_term(nv + 4, nv + 4)
    for c in range(4):
        for i in range(5):
            ok &= (ratio_t_junc(i, c, nv) == t_term(nv + i, nv + c) / base)
check("junction ratio closed forms match direct factorial values", ok)

# per-addend structural degrees: a_i (21,0)/(0,-); ratio_t num deg <= 3+8+4=15, den deg 16;
# Rt(n,n+4): num 31 (Xhat with rho=n+4), den 8. 21 addends: use the generic bound helper
# in one variable by treating deg_r = 0 throughout.
terms_C = []
for c in range(4):
    for i in range(5):
        terms_C.append((
            (lambda i=i, c=c: lambda nv, rv: cL(i, nv) * ratio_t_junc(i, c, nv))(),
            (21 + 15, 0), (16, 0),   # num deg 15 and den deg 16 are conservative upper bounds
        ))
terms_C.append((lambda nv, rv: -Rt(nv, nv + 4), (31, 0), (8, 0)))  # 31/8 conservative (actual <=22/5); one-sided, soundness unaffected
BJ = bound_from_terms(terms_C, "n")
ok = True
for nv in range(300, 300 + BJ + 1):
    s = F(0)
    for ev, _, _ in terms_C:
        s += ev(nv, 0)
    if s != 0:
        ok = False
        break
check("(C) junction identity on %d points (degree bound BJ=%d) -> holds in QQ(n)"
      % (BJ + 1, BJ), ok)

# ============================================================================
# (7) exhaustive telescoping
# ============================================================================
def T_L(nv, rv):
    """(L t)(n,rho) = sum_k c_k(n) t(n+k,rho) -- the combination the certificate
    telescopes (the shipped Xhat is oriented for +L; identities (D),(C) above)."""
    return sum(cL(k, nv) * t_term(nv + k, rv) for k in range(5))

ok = True
for nv in range(5, 41):
    def Gstar(rv):
        if rv <= nv or rv >= 2 * nv + 9:
            return F(0)
        if rv >= nv + 4:
            return Rt(nv, rv) * t_term(nv + 4, rv)
        return sum(T_L(nv, q) for q in range(nv, rv))
    tot = F(0)
    for rv in range(nv - 2, 2 * nv + 12):
        v = T_L(nv, rv)
        if v != Gstar(rv + 1) - Gstar(rv):
            ok = False
        tot += v
    if tot != 0:
        ok = False
check("telescoping holds at EVERY integer rho and sums to zero, n=5..40", ok)

# ============================================================================
# (9) paper-form equivalence (arXiv:2303.02793 Conjecture 17 <-> OEIS form)
# ============================================================================
P_SRC = {
 0: "(-64*n**10 - 1968*n**9 - 26156*n**8 - 198469*n**7 - 952323*n**6 - 3012795*n**5 - 6333869*n**4 - 8663374*n**3 - 7264534*n**2 - 3266000*n - 549760)",
 1: "(64*n**13 + 2672*n**12 + 49788*n**11 + 545913*n**10 + 3917758*n**9 + 19359535*n**8 + 67385886*n**7 + 165789363*n**6 + 284054698*n**5 + 325846005*n**4 + 229526554*n**3 + 78563984*n**2 - 487964*n - 5543040)",
 2: "(-512*n**15 - 21568*n**14 - 419248*n**13 - 4969164*n**12 - 39928763*n**11 - 228837227*n**10 - 959068672*n**9 - 2966908118*n**8 - 6753094929*n**7 - 11118771121*n**6 - 12741784568*n**5 - 9313604242*n**4 - 3271711596*n**3 + 562569136*n**2 + 946158512*n + 250467360)",
 3: "2*(n+3)*(512*n**16 + 26752*n**15 + 624800*n**14 + 8677944*n**13 + 80260596*n**12 + 523718876*n**11 + 2488583381*n**10 + 8747566435*n**9 + 22820793074*n**8 + 43766004538*n**7 + 60004107039*n**6 + 55047935941*n**5 + 27672902302*n**4 - 778719870*n**3 - 10812498240*n**2 - 6360099840*n - 1300242000)",
 4: "-12*(n+4)**3*(n+3)*(2*n+7)**2*(3*n+8)*(3*n+10)*(64*n**10 + 1328*n**9 + 11324*n**8 + 52389*n**7 + 143536*n**6 + 233810*n**5 + 204716*n**4 + 48699*n**3 - 68928*n**2 - 61278*n - 15900)",
}
P_CODE = {k: compile(v, "<P%d>" % k, "eval") for k, v in P_SRC.items()}
def cP(k, nv):
    return eval(P_CODE[k], {"n": F(nv)})

# claim: cP_i(n) * [(2n+2i+1)...(2n+8)]^2 == dpoly(n) * cL_i(n) as polynomials
# (degrees <= 26 + 16 = 42 on the left, 5 + 21 = 26 on the right; check 50 points)
ok = True
for nv in range(100, 151):
    for i in range(5):
        fac = F(1)
        for s in range(2 * i + 1, 9):
            fac *= (2 * F(nv) + s)
        if cP(i, nv) * fac ** 2 != dpoly(nv) * cL(i, nv):
            ok = False
check("paper (atilde) form == OEIS form times 4(n+3)(n+4)^2(2n+7)^2 (51 points, degs <= 42)", ok)

print()
print("A269021 VERIFICATION: ALL %d CHECKS PASS" % len(PASS))
print("=> the representation is correct on all independent ground truths, the WZ")
print("   certificate identities (D) and (C) hold in QQ(n,rho), telescoping is exact,")
print("   and hence the Kauers-Koutschan recurrence (Conjecture 17) annihilates a(n).")
sys.exit(0)
