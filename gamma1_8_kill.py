"""Real-G push, step 6 (2026-09-29): cusp-killing on Gamma_1(8).
Frame T = (1 - theta3(q^2)/theta3(q))/2, f = theta3(q)^2 = sum u_n T^n; the ODE recurrence (5 terms) has singular points
T = 1/x, x in {4+2sqrt2 (dominant), 2, 4-2sqrt2, 1}, and the pole cusp at T = oo.
For inhomogeneities h (integral Taylor coefficients) solve  [recurrence] y = h  (y_0 = 0), Apery limit Lam(h) = lim y_n/u_n,
linear forms F_n(h) = y_n - Lam(h) u_n, amplitude at the cusp T = 1/2:  A(h) = lim n F_n / 2^n  (log singularity).
Identify Lam(h) in 1, G, L(2,chi_-8), pi^2, L(g,2), L(g,1) (g = the CM form eta(t)^2 eta(2t) eta(4t) eta(8t)^2).
Then look for rational combinations with Lam in Q + QG and A = 0 (the cusp 1/2 killed).
Usage: python gamma1_8_kill.py [N]"""
import sys, io, contextlib, math
from fractions import Fraction as Fr
import mpmath as mp
import numpy as np
sys.argv = [sys.argv[0], '8'] + sys.argv[1:]
import level12_modular as L
from guess_rec import guess
N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
mp.mp.dps = 260
P, M = 140, 120

def theta3(m):
    s = [0]*P; n = 0
    while m*n*n < P:
        s[m*n*n] += 1 if n == 0 else 2; n += 1
    return s

th1, th2 = theta3(1), theta3(2)
s = L.mul(th2, L.inv(th1, P), P)
T = [((1 if i == 0 else 0) - s[i])//2 for i in range(P)]
f = L.mul(th1, th1, P)
u = L.expand_in(f, T, M)
with contextlib.redirect_stdout(io.StringIO()):
    R, d, cs = guess(u, R=8, D=2, verbose=False)
ev = lambda c, n: sum(ci*n**e for e, ci in enumerate(c))
print('recurrence (forward, c_k(n) a_{n+k}):', cs)

def solve(h, upto):
    """sum_k c_k(n-R) y_{n-R+k} = h_n for n >= 1, y_0 = 0 (h_0 must be 0), zero extension."""
    y = [Fr(0)]
    for n in range(1, upto):
        m = n - R
        lead = ev(cs[R], m)
        s_ = sum(ev(cs[k], m)*(y[m + k] if m + k >= 0 else 0) for k in range(R))
        y.append((Fr(h[n]) - s_)/lead)
    return y

# u from the homogeneous recurrence with u_0 = 1
Uh = L.run(cs, [u[0]], N)
assert all(Uh[i] == u[i] for i in range(M))

def series_h(kind, j, upto):
    """Taylor coefficients of T^j, T^j/(1-2T), T^j/(1-T), T^j/(1-8T+8T^2)."""
    base = [0]*upto
    if kind == 'poly':
        base[j] = 1; return base
    den = {'p2': [1, -2], 'p1': [1, -1], 'pq': [1, -8, 8]}[kind]
    inv = [0]*upto; inv[0] = 1
    for n in range(1, upto):
        inv[n] = -sum(den[k]*inv[n - k] for k in range(1, min(len(den), n + 1)))
    for n in range(j, upto): base[n] = inv[n - j]
    return base

s2 = mp.sqrt(2); pi = mp.pi; G = +mp.catalan
L8 = mp.nsum(lambda k: 1/(8*k + 1)**2 + 1/(8*k + 3)**2 - 1/(8*k + 5)**2 - 1/(8*k + 7)**2, [0, mp.inf])
PC = 900
g = [0] + L.eta_series({1: 2, 2: 1, 4: 1, 8: 2}, PC)[:PC - 1]
def lval(a, k, Nl, sv, eps=1):
    A = mp.sqrt(Nl)/(2*mp.pi); sv = mp.mpf(sv); tot = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        x = 2*mp.pi*n/mp.sqrt(Nl)
        tot += a[n]*(mp.gammainc(sv, x)/mp.gamma(sv)/mp.mpf(n)**sv + eps*A**(k - 2*sv)*mp.gammainc(k - sv, x)/mp.gamma(sv)*mp.mpf(n)**(sv - k))
    return tot
Lg1, Lg2 = lval(g, 3, 8, 1), lval(g, 3, 8, 2)
BASIS = [('1', mp.mpf(1)), ('G', G), ('L8', L8), ('pi^2', pi**2), ('Lg2', Lg2), ('Lg1', Lg1), ('pi^2*s2', pi**2*s2), ('s2*G', s2*G)]

q = lambda x: mp.mpf(x.numerator)/x.denominator
Uq = [q(x) for x in Uh]
rows = []
for kind in ('poly', 'p2', 'p1'):          # (T^j/(1-8T+8T^2) has a pole at the dominant cusp: no Apery limit)
    for j in range(1, 5):
        h = series_h(kind, j, N)
        y = solve(h, N)
        Lam = q(y[-1])/Uq[-1]; e = abs(Lam - q(y[-31])/Uq[-31]); dig = int(-mp.log10(e)) if e > 0 else 250
        # amplitude at T = 1/2: n F_n / 2^n with Richardson in 1/n
        ns = list(range(150, 170))
        vals = [n*(q(y[n]) - Lam*Uq[n])/mp.mpf(2)**n for n in ns]
        # fit a + b/n + c/n^2 + d/n^3 by least squares
        Amat = mp.matrix([[1, 1/mp.mpf(n), 1/mp.mpf(n)**2, 1/mp.mpf(n)**3, 1/mp.mpf(n)**4] for n in ns])
        coef = mp.lu_solve(Amat.T*Amat, Amat.T*mp.matrix(vals))
        Aamp = coef[0]
        rel = None
        if dig >= 60:
            with mp.workdps(min(dig - 20, 200)):
                rel = mp.pslq([Lam] + [b for _, b in BASIS], maxcoeff=10**6, maxsteps=10**7)
        idtxt = ('%d*Lam = %s' % (rel[0], ' + '.join('%d*%s' % (-c, nm) for c, (nm, _) in zip(rel[1:], BASIS) if c))) if rel else 'no relation'
        den_growth = math.log(y[-1].denominator)/(N - 1)
        print('%-4s j=%d  Lam=%s (%d dig)  tau=%.3f  A(1/2)=%s   %s' % (kind, j, mp.nstr(Lam, 20), dig, den_growth, mp.nstr(Aamp, 15), idtxt))
        sys.stdout.flush()
        rows.append((kind, j, Lam, Aamp, rel))
