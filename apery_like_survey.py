"""Real-G push, step 2 (2026-09-29): which second-order Apery-like recurrences approximate something involving G?
For (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}   (Zagier's family; u_0 = 1, u_1 = b),
second solution v_0 = 0, v_1 = 1:  identify lim v_n/u_n by PSLQ against 1, G, pi^2, pi*ln2, ln2^2, L(2,chi_-3)*sqrt3... ,
and report: integrality of u_n (to n = 60), the characteristic roots (alpha, beta), ln(den v_n)/n, and the CDT-type net
   net = tau - ln r,  tau = ln(den v_n)/n (measured),  r = 1/|beta| (radius of the linear-form generating function),
compared with ln 16 = 2.773 (the rule read off CDT Remark 11.1.17; to be confirmed)."""
import sys, math, cmath
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 60
L3 = mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])   # L(2, chi_-3)
basis = [('1', mp.mpf(1)), ('G', +mp.catalan), ('pi^2', mp.pi**2), ('pi*ln2', mp.pi*mp.log(2)), ('ln2^2', mp.log(2)**2),
         ('L(2,chi-3)', L3), ('pi*sqrt3', mp.pi*mp.sqrt(3)), ('pi', mp.pi), ('ln2', mp.log(2)), ('ln3', mp.log(3))]
# Zagier's sporadic cases A-F and the four hypergeometric ones (a, b, c)
cases = {'A (7,2,-8)': (7, 2, -8), 'B (9,3,27)': (9, 3, 27), 'C (10,3,9)': (10, 3, 9), 'D (11,3,-1)': (11, 3, -1),
         'E (12,4,32)': (12, 4, 32), 'F (17,6,72)': (17, 6, 72),
         'h1 (0,0,-16)': (0, 0, -16), 'h2 (32,12,256)': (32, 12, 256), 'h3 (-16? 16,4,64)': (16, 4, 64), 'h4 (27,9,729)': (27, 9, 729)}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
for name, (a, b, c) in cases.items():
    u = [F(1), F(b)]; v = [F(0), F(1)]
    for n in range(1, N):
        A = a*n*n + a*n + b; C = c*n*n; D = (n + 1)**2
        u.append((A*u[n] - C*u[n - 1])/D); v.append((A*v[n] - C*v[n - 1])/D)
    integral = all(x.denominator == 1 for x in u[:61])
    # characteristic roots of x^2 - a x + c
    disc = a*a - 4*c
    r1 = (a + cmath.sqrt(disc))/2; r2 = (a - cmath.sqrt(disc))/2
    al, be = (r1, r2) if abs(r1) >= abs(r2) else (r2, r1)
    try:
        Lv = mp.mpf(v[N].numerator)/v[N].denominator/(mp.mpf(u[N].numerator)/u[N].denominator)
        Lv2 = mp.mpf(v[N - 30].numerator)/v[N - 30].denominator/(mp.mpf(u[N - 30].numerator)/u[N - 30].denominator)
        conv = mp.nstr(abs(Lv - Lv2), 3)
        rel = mp.pslq([Lv] + [x for _, x in basis], maxcoeff=10**5, maxsteps=10**6)
        relstr = ' + '.join('%d*%s' % (k, nm) for k, (nm, _) in zip(rel[1:], basis) if k) + '  (coeff of L: %d)' % rel[0] if rel else 'none found'
    except Exception as ex:
        Lv, conv, relstr = None, '-', 'error %s' % ex
    tau = math.log(v[N].denominator)/N if v[N].denominator > 1 else 0.0
    rr = 1/abs(be) if abs(be) > 1e-12 else float('inf')
    net = tau - math.log(rr) if rr != float('inf') else float('-inf')
    print('%-18s u integral: %-5s roots |%.3f|, |%.3f|  tau = %.3f  r = %.4f  net = %+.3f  (ln16 margin %+.3f)' % (
        name, integral, abs(al), abs(be), tau, rr, net, math.log(16) - net))
    print('    lim v/u = %s  (conv %s);  PSLQ: %s' % (mp.nstr(Lv, 25) if Lv is not None else '-', conv, relstr))
    sys.stdout.flush()
