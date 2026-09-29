r"""companion_k2.py (2026-09-29): the companion mechanism in the K^2 Hankel arena, exact content test.

Half-integer poles b = 1, 3, ..., 2K-1 (N = 0), as in padic_test.py.  Catalan machine: pole values
PV_cat(b) = A + G*B with A = -(chi(b) b/2) beta_alt(b) - 1/(4b), B = chi(b) b/2, moments 4^e (2^{2e+2}-1)|B_{2e+2}| (Genocchi).
Plain machine: PV_z(b) = -(b/4) beta(b) - 1/(8b) - 1/8 (+ (b/16)*Y, Y = pi^2/2), moments 4^e |B_{2e+2}|/2 (Bernoulli).
Companion machine = Catalan + c * plain, both in the pole values and in the moments (one kernel), with the plain unknown Y
set to 0: only the content of det(A + X B) is tested here (positivity and the second unknown are separate issues).
At a pole b > p the 1/p^2 coefficient of the mixed pole value is -(b/4)(2 chi(b) chi(p) + c):
  c = -2 cancels it exactly for the poles b = p (mod 4)  (half the poles of every prime).
Ledger prediction (top range K < p < 2K, n = #poles above p, every mirror pair straddles the classes):
  saved pair 2 (Vandermonde) - 1 (Hermite: both members p-integral) + 0 + 1 (von Staudt) = 2, unsaved 2 + 1 + 1 = 4,
  so e_p = 3n (+1 for b = p) again: an exact wash against Catalan's 3n + 1.
usage: python companion_k2.py K c"""
import sys, math
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, primerange
import flint
sys.set_int_max_str_digits(0)
K = int(sys.argv[1]); c = F(sys.argv[2]) if len(sys.argv) > 2 else F(-2); h = K
def bern(n): B = bernoulli(n); return F(int(B.p), int(B.q))
def polymul(a, b):
    r = [0]*(len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b): r[i + j] += ai*bj
    return r
labels = [2*j + 1 for j in range(K)]
mu_cat = lambda e: F(4)**e*(2**(2*e + 2) - 1)*abs(bern(2*e + 2))
mu_z = lambda e: F(4)**e*abs(bern(2*e + 2))/2
beta = {}; betaalt = {}; acc = F(0); acc2 = F(0); sgn = 1
for b in range(1, 2*K + 2, 2):
    beta[b] = acc; betaalt[b] = acc2; acc += F(1, b*b); acc2 += F(sgn, b*b); sgn = -sgn
def pv(b):
    s = (-1)**((b - 1)//2)
    A_cat = F(-s*b, 2)*betaalt[b] - F(1, 4*b); B_cat = F(s*b, 2)
    A_z = -F(b, 4)*beta[b] - F(1, 8*b) - F(1, 8)
    return (A_cat + c*A_z, B_cat)
D = [1]
for x in labels: D = polymul(D, [x*x, 1])
degD = len(D) - 1
Dp = {a: math.prod(b*b - a*a for b in labels if b != a) for a in labels}
PVs = {a: pv(a) for a in labels}
MU = [mu_cat(e) + c*mu_z(e) for e in range(2*h + 2)]
q = []; r = [1] + [0]*(degD - 1)
pw = {a: 1 for a in labels}; A = []; B = []
for s in range(2*h - 1):
    a_val = F(0); b_val = F(0)
    for e, cc in enumerate(q):
        if cc: a_val += MU[e]*cc
    for a in labels:
        cf = F(pw[a], Dp[a]); a_val += cf*PVs[a][0]; b_val += cf*PVs[a][1]
    A.append(a_val); B.append(b_val)
    cs = r[degD - 1]; q = [cs] + q
    r = [(r[i - 1] if i > 0 else 0) - cs*D[i] for i in range(degD)]
    for a in labels: pw[a] *= -a*a
Af = flint.fmpq_mat([[flint.fmpq(A[i + j].numerator, A[i + j].denominator) for j in range(h)] for i in range(h)])
Bf = flint.fmpq_mat([[flint.fmpq(B[i + j].numerator, B[i + j].denominator) for j in range(h)] for i in range(h)])
vals = []
for X in range(h + 1):
    d = (Af + Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p), int(d.q)))
xs = list(range(h + 1)); coef = vals[:]
for k in range(1, h + 1):
    for i in range(h, k - 1, -1): coef[i] = (coef[i] - coef[i - 1])/(xs[i] - xs[i - k])
poly = [coef[h]]
for k in range(h - 1, -1, -1):
    poly = polymul(poly, [F(-xs[k]), F(1)]); poly = [F(x) for x in poly]; poly[0] += coef[k]
g = 0; l = 1
for x in poly: g = gcd(g, x.numerator); l = l*x.denominator//gcd(l, x.denominator)
den = l//gcd(l, g)
rows = {1: [], 3: []}; tot = {1: [0, 0], 3: [0, 0]}
for p in primerange(3, 2*K):
    e = 0
    while den % p == 0: den //= p; e += 1
    if p > K:
        n = sum(1 for b in labels if b > p)
        rows[p % 4].append("%d:%d[%+d]" % (p, e, e - (3*n + 1)))
        tot[p % 4][0] += e; tot[p % 4][1] += 3*n + 1
print("companion K^2 machine, K = %d, c = %s: content exponents e_p for K < p < 2K, [e_p - (3n+1)] (Catalan's law)" % (K, c))
for cls in (1, 3):
    print("  p = %d mod 4: %s" % (cls, " ".join(rows[cls])))
    print("     sum e_p = %d  vs Catalan's sum 3n+1 = %d" % (tot[cls][0], tot[cls][1]))
