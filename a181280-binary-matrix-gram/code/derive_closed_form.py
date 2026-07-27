#!/usr/bin/env python3
"""
INDEPENDENT recovery of the recurrence + closed form for A181280, from the terms produced
by machine_fourier_coupled.py (the genuinely-different Fourier-coupled machine).  We run our
OWN exact Berlekamp-Massey (over Q) to recover, from this FINITE SAMPLE of terms, the shortest
linear recurrence consistent with it; factor it; do partial fractions to get the closed form;
and compare symbolically to Kauers-Koutschan C20 and to the char poly
p(x)=(x-16)(x-8)^2(x+8)(x-4)^3(x+4)^2(x-2)^2. NOTE: Berlekamp-Massey returns the minimal order
of the OBSERVED SAMPLE; that it equals the true minimal order of the infinite tail is proved
separately in proof_matrix.py (det H_11 + a-priori bound + residual propagation). This script
is an independent cross-check, not the minimality proof.
"""
import json
from fractions import Fraction
import sympy as sp

a = json.load(open("a_terms_fourier.json"))   # n=0..30, from the independent Fourier machine

# ---- our own exact Berlekamp-Massey over Q on the tail n>=4 (C-finite region) ----
def berlekamp_massey_Q(seq):
    seq = [Fraction(x) for x in seq]
    C=[Fraction(1)]; B=[Fraction(1)]; L=0; m=1; b=Fraction(1)
    for nidx in range(len(seq)):
        # discrepancy
        d = seq[nidx]
        for i in range(1, L+1):
            d += C[i]*seq[nidx-i]
        if d == 0:
            m += 1
        elif 2*L <= nidx:
            T=C[:]
            coef = d/b
            while len(C) < len(B)+m: C.append(Fraction(0))
            for i in range(len(B)):
                C[i+m] -= coef*B[i]
            L = nidx+1-L; B=T; b=d; m=1
        else:
            coef = d/b
            while len(C) < len(B)+m: C.append(Fraction(0))
            for i in range(len(B)):
                C[i+m] -= coef*B[i]
            m += 1
    return C, L  # recurrence: sum_{i=0}^{L} C[i] seq[n-i] = 0

tail = a[4:]
C, L = berlekamp_massey_Q(tail)
print("Berlekamp-Massey order of the finite sample (tail n>=4):", L)
# characteristic polynomial: x^L + C[1] x^{L-1} + ... + C[L]
x = sp.symbols('x')
charpoly = sum(sp.Rational(C[i].numerator, C[i].denominator)*x**(L-i) for i in range(L+1))
charpoly = sp.expand(charpoly)
print("char poly (monic):", charpoly)
print("factored:", sp.factor(charpoly))

# prover's claimed p(x)
p_claim = (x-16)*(x-8)**2*(x+8)*(x-4)**3*(x+4)**2*(x-2)**2
print("prover p(x) factored:", sp.factor(p_claim))
print("our charpoly == prover p(x)/(leading) ?",
      sp.simplify(sp.expand(p_claim) - charpoly*sp.LC(sp.Poly(p_claim,x))) == 0
      or sp.simplify(sp.factor(p_claim) - sp.factor(charpoly))==0)

# ---- independent closed form via solving the C-finite ansatz with the roots ----
# roots with multiplicities from OUR factored char poly
charP = sp.Poly(charpoly, x)
roots = charP.all_roots()
from collections import Counter
rootmult = {}
for r in sp.roots(charpoly, x).items():
    rootmult[r[0]] = r[1]
print("roots (with mult):", rootmult)

# Build ansatz a(n) = sum_root sum_{j<mult} c_{r,j} n^j r^n, solve with first L tail values
n = sp.symbols('n', integer=True)
coeffs = []
terms = []
for r, mult in rootmult.items():
    for j in range(mult):
        cc = sp.symbols(f'c_{sp.srepr(r)}_{j}'.replace("(","").replace(")","").replace("-","m").replace("Integer","I").replace(",","_").replace(" ",""))
        coeffs.append(cc)
        terms.append(cc * n**j * r**n)
ansatz = sum(terms)
eqs = []
for k in range(L):
    nval = 4+k
    eqs.append(sp.Eq(ansatz.subs(n, nval), tail[k]))
sol = sp.solve(eqs, coeffs, dict=True)[0]
closed = sp.simplify(ansatz.subs(sol))
print("independent closed form a(n) (n>=4):")
sp.pprint(closed)

# ---- compare to Kauers-Koutschan C20 symbolically ----
C20 = (sp.Rational(1,3)*2**(2*n-11)*(6*n**2-219*n+820)
       - sp.Rational(1,9)*2**(n-5)*(3*n+32)
       - sp.Rational(113,3)*(-1)**n*2**(3*n-14)
       + 2**(4*n-9)
       - sp.Rational(1,3)*(-1)**n*2**(2*n-11)*(13*n-164)
       + sp.Rational(1,9)*2**(3*n-14)*(288*n-3473))

diff = sp.simplify(closed - C20)
print("closed_form - C20 simplified to:", diff)
# robust check at many integer points
mism = [k for k in range(4,40) if sp.nsimplify(closed.subs(n,k)) != sp.nsimplify(C20.subs(n,k))]
print("integer points 4..39 where closed != C20:", mism)

# also confirm our closed form reproduces the independent terms exactly
mism2 = [k for k in range(4, len(a)) if int(closed.subs(n,k)) != a[k]]
print("points where closed != independent terms a[]:", mism2)
