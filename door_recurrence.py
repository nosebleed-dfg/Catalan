r"""The non-congruence door, part 4 (2026-09-29): the exact recurrence of the index-9 (1,2,6) family, its limits to high
precision, and the targeted tests.

Exact u_n = [t^n] f and v_{k,n} = [t^n] f g_k (k = 1, 2, 3) from door_family.build(Fraction) (n <= M0), then the minimal
recurrence of u guessed over Q (order r, polynomial degree d; linear algebra on the first terms, checked on all of them).
v_k obeys the same recurrence away from its inhomogeneity; it is continued from exact terms.  Limits Lambda_k = lim v_k/u
to hundreds of digits; tests against small declared bases.
Usage: python door_recurrence.py [M0] [NMAX]"""
import sys, math
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

M0 = int(sys.argv[1]) if len(sys.argv) > 1 else 44
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 700
src = open(r"C:\Users\PC\Desktop\catalan\door_family.py", encoding="utf-8").read().split('print("(i) limits')[0]
sys.argv = ["door_family.py", str(M0)]
ns = {}; exec(src, ns)
P = ns['P']; build = ns['build']; expand_in = ns['expand_in']; mul = ns['mul']; inv = ns['inv']
t, f, Dt, tp, c = build(Fr)
u = expand_in(f, tp)
t0 = Fr(27, 256); b = Fr(27, 64)
p = [(-1/t0 - 1/b)*x for x in t]; p[0] += 1
tt = mul(t, t); p = [p[i] + tt[i]/(t0*b) for i in range(P)]
base = mul(mul(Dt, Dt), inv(mul(p, f)))
V = []
for k in (1, 2, 3):
    W = base if k == 1 else mul(base, tp[k - 1])
    g = [Fr(0)] + [W[n]/Fr(n*n) for n in range(1, P)]
    V.append(expand_in(mul(f, g), tp))

def guess(seq, r, d, start=0):
    """find c_{j,e} with sum_j sum_e c_{j,e} n^e seq[n+j] = 0 for n >= start"""
    ncol = (r + 1)*(d + 1)
    rows = []
    for n in range(start, len(seq) - r):
        rows.append([Fr(n)**e*seq[n + j] for j in range(r + 1) for e in range(d + 1)])
    M = sp.Matrix(rows)
    ns_ = M.nullspace()
    return ns_
for r, d in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4)):
    ker = guess(u, r, d)
    print("order %d degree %d: kernel dimension %d (from %d equations)" % (r, d, len(ker), len(u) - r))
    if len(ker) == 1:
        vec = ker[0]; vec = vec/[x for x in vec if x != 0][-1]
        cs = [[Fr(int(sp.fraction(vec[j*(d + 1) + e])[0]), int(sp.fraction(vec[j*(d + 1) + e])[1])) for e in range(d + 1)] for j in range(r + 1)]
        n_ = sp.symbols('n')
        for j in range(r + 1):
            print("   c_%d(n) = %s" % (j, sp.factor(sum(sp.Rational(cs[j][e].numerator, cs[j][e].denominator)*n_**e for e in range(d + 1)))))
        break
R, D = r, d
def coef(j, n): return sum(cs[j][e]*n**e for e in range(D + 1))
# check v_k against the same recurrence: residuals by n
for k in range(3):
    res = [sum(coef(j, n)*V[k][n + j] for j in range(R + 1)) for n in range(0, len(u) - R)]
    nz = [n for n, x in enumerate(res) if x != 0]
    print("   v_%d: recurrence residual nonzero at n = %s" % (k + 1, nz))
# continue u and v_k to NMAX in mpmath
mp.mp.dps = 400
def cont(seq, n_start):
    s = [mp.mpf(x.numerator)/x.denominator for x in seq[:n_start + R]]
    for n in range(n_start, NMAX):
        # coef(R, n) s[n+R] = - sum_{j<R} coef(j,n) s[n+j]
        acc = sum(mp.mpf(coef(j, n).numerator)/coef(j, n).denominator*s[n + j] for j in range(R))
        cR = coef(R, n); s.append(-acc/(mp.mpf(cR.numerator)/cR.denominator))
    return s
start = len(u) - R - 1
U = cont(u, start)
LAM = []
for k in range(3):
    Vk = cont(V[k], start)
    lam = Vk[NMAX - 1]/U[NMAX - 1]; lam2 = Vk[NMAX - 51]/U[NMAX - 51]
    LAM.append(lam)
    print("Lambda_%d = %s  (change over 50 steps: %s)" % (k + 1, mp.nstr(lam, 60), mp.nstr(abs(lam - lam2), 3)))
dig = int(-mp.log10(abs(LAM[0] - cont(V[0], start)[NMAX - 51]/U[NMAX - 51]) + mp.mpf(10)**-390))
mp.mp.dps = min(dig - 10, 300)
print("usable digits: %d" % mp.mp.dps)
print("relation Lambda_1, Lambda_2, Lambda_3, 1:", mp.pslq(LAM + [mp.mpf(1)], maxcoeff=10**15, maxsteps=10**7))
G = mp.catalan; pi = mp.pi; Om = mp.gamma(mp.mpf(1)/4)**4/pi**2
L3 = mp.nsum(lambda j: 1/(3*j + 1)**2 - 1/(3*j + 2)**2, [0, mp.inf]); L8 = mp.nsum(lambda j: 1/(8*j+1)**2 + 1/(8*j+3)**2 - 1/(8*j+5)**2 - 1/(8*j+7)**2, [0, mp.inf])
tests = [("1, G", [1, G]), ("1, pi^2", [1, pi**2]), ("1, G, pi^2", [1, G, pi**2]), ("1, Om", [1, Om]), ("1, 1/Om", [1, 1/Om]),
         ("1, pi, pi^2", [1, pi, pi**2]), ("1, G, pi^2, pi log2, log2^2", [1, G, pi**2, pi*mp.log(2), mp.log(2)**2]),
         ("1, L(2,chi-3), pi^2", [1, L3, pi**2]), ("1, L(2,chi-8), pi^2", [1, L8, pi**2]), ("1, G, pi^2, Om, 1/Om", [1, G, pi**2, Om, 1/Om]),
         ("1, zeta3, pi^2, G", [1, mp.zeta(3), pi**2, G]), ("1, log2, log3, pi*sqrt3", [1, mp.log(2), mp.log(3), pi*mp.sqrt(3)])]
for k in (0, 1):
    for name, vals in tests:
        rel = mp.pslq([LAM[k]] + [mp.mpf(v) for v in vals], maxcoeff=10**20, maxsteps=10**7)
        print("   Lambda_%d over {%s}: %s" % (k + 1, name, rel))
print("the span of the limits: does a combination of Lambda_1, Lambda_2 reach G?")
for name, vals in [("1, G", [1, G]), ("1, pi^2", [1, pi**2]), ("1, G, pi^2", [1, G, pi**2]), ("1, G, pi^2, Om", [1, G, pi**2, Om]),
                   ("1, G, pi^2, Om, 1/Om", [1, G, pi**2, Om, 1/Om])]:
    rel = mp.pslq([LAM[0], LAM[1]] + [mp.mpf(v) for v in vals], maxcoeff=10**20, maxsteps=10**7)
    print("   {Lambda_1, Lambda_2, %s}: %s" % (name, rel))
# the integer-coefficient form of the recurrence: w_n = (27/4)^n u_n
print("w_n = (27/4)^n u_n, first terms:", [str(Fr(27, 4)**n*u[n]) for n in range(8)])
