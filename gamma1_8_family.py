"""Real-G push, step 5 (2026-09-29): escape the q -> -q symmetry of Gamma_0(4k) on Gamma_1(8).
Hauptmodul T = (1 - theta3(q^2)/theta3(q))/2 = q - 3q^2 + 6q^3 - ...  (integral; cusp values x = 1/T in {4+2sqrt2, 2, 4-2sqrt2, 1}
and the pole). Weight-1 forms theta3(q)^2 (chi_-4), theta3(q)theta3(q^2) (chi_-8), theta3(q^2)^2 (chi_-4).
For each form and frame T_a = T/(1 + aT): recurrence, Apery second solution v (L v = T), limit, tau, measured decay, CDT gate,
and PSLQ against 1, G, L(2,chi_-8), pi^2, pi^2*sqrt2, the CM period L(g,2) (g = eta(t)^2 eta(2t) eta(4t) eta(8t)^2, weight 3,
level 8, chi_-8), and sqrt2-multiples.
Usage: python gamma1_8_family.py [NN]"""
import sys, io, contextlib, math
import mpmath as mp
import numpy as np
sys.argv = [sys.argv[0], '8'] + sys.argv[1:]
import level12_modular as L
from guess_rec import guess
NN = int(sys.argv[2]) if len(sys.argv) > 2 else 700
mp.mp.dps = 170
P, M = 140, 120

def theta3(m):
    s = [0]*P; n = 0
    while m*n*n < P:
        s[m*n*n] += 1 if n == 0 else 2; n += 1
    return s

def lval(a, k, N, s, eps=1, terms=None):
    """L(f, s) by the approximate functional equation (self-dual newform, sign eps)."""
    A = mp.sqrt(N)/(2*mp.pi); s = mp.mpf(s); tot = mp.mpf(0)
    for n in range(1, terms or len(a)):
        if a[n] == 0: continue
        x = 2*mp.pi*n/mp.sqrt(N)
        tot += a[n]*(mp.gammainc(s, x)/mp.gamma(s)/mp.mpf(n)**s + eps*A**(k - 2*s)*mp.gammainc(k - s, x)/mp.gamma(s)*mp.mpf(n)**(s - k))
    return tot

th1, th2 = theta3(1), theta3(2)
s = L.mul(th2, L.inv(th1, P), P)
T = [((1 if i == 0 else 0) - s[i])//2 for i in range(P)]
forms = [('th3(q)^2', L.mul(th1, th1, P)), ('th3(q)th3(q^2)', L.mul(th1, th2, P)), ('th3(q^2)^2', L.mul(th2, th2, P))]

# constants
s2 = mp.sqrt(2); pi = mp.pi; G = +mp.catalan
L8 = mp.nsum(lambda k: 1/(8*k + 1)**2 + 1/(8*k + 3)**2 - 1/(8*k + 5)**2 - 1/(8*k + 7)**2, [0, mp.inf])   # L(2, chi_-8)
PC = 600
g = [0] + L.eta_series({1: 2, 2: 1, 4: 1, 8: 2}, PC)[:PC - 1]
Lg1, Lg2 = lval(g, 3, 8, 1), lval(g, 3, 8, 2)
print('L(2,chi_-8) = %s   L(g,1) = %s   L(g,2) = %s   (check L(g,1.5) = %s)' % (mp.nstr(L8, 20), mp.nstr(Lg1, 20), mp.nstr(Lg2, 20), mp.nstr(lval(g, 3, 8, 1.5), 12)))
BASIS = [('1', 1), ('G', G), ('L8', L8), ('pi^2', pi**2), ('pi^2*s2', pi**2*s2), ('Lg2', Lg2), ('Lg1', Lg1), ('pi', pi),
         ('pi*s2', pi*s2), ('s2', s2), ('s2*G', s2*G), ('s2*L8', s2*L8)]

for a in (0, -2, -1, 1):
    den = [1] + [a*x for x in T[1:]]
    Ta = L.mul(T, L.inv(den, P), P)
    for nm, f in forms:
        u = L.expand_in(f, Ta, M)
        with contextlib.redirect_stdout(io.StringIO()):
            sol = guess(u, R=8, D=2, verbose=False)
        if sol is None:
            print('a=%+d %-15s no recurrence' % (a, nm)); continue
        R, d, cs = sol
        lc = [c[2] if len(c) > 2 else 0 for c in cs]
        rts = sorted(np.roots(lc[::-1]), key=lambda z: -abs(z))
        U = L.run(cs, [u[0]], NN)
        if U is None or any(U[i] != u[i] for i in range(M)):
            print('a=%+d %-15s zero-extension fails' % (a, nm)); continue
        V = L.run(cs, [0, 1], NN)
        q = lambda x: mp.mpf(x.numerator)/x.denominator
        Lv = q(V[-1])/q(U[-1]); e = abs(Lv - q(V[-41])/q(U[-41]))
        dig = int(-mp.log10(e)) if e > 0 else 165
        # F_n = v_n - L u_n is only reliable while (x_2/x_1)^n stays above the accuracy of L
        ratio = abs(rts[0])/abs(rts[1]) if abs(rts[1]) > 0 else 10.0
        nmax = max(40, min(int(0.6*NN), int((dig - 20)/math.log10(ratio))))
        Fs = [abs(q(V[n]) - Lv*q(U[n])) for n in range(nmax)]
        rate = float((mp.log(Fs[nmax - 1]) - mp.log(Fs[nmax//2]))/(nmax - 1 - nmax//2))
        tau = math.log(V[-1].denominator)/(NN - 1)
        S = math.log(16) - rate - tau
        rel = None
        for k in range(2, len(BASIS) + 1):
            with mp.workdps(min(dig - 15, 150)):
                rr = mp.pslq([Lv] + [mp.mpf(b) for _, b in BASIS[:k]], maxcoeff=10**5, maxsteps=10**7)
            if rr and rr[0]:
                rel = '%d*L = %s' % (rr[0], ' + '.join('%d*%s' % (-c, n2) for c, (n2, _) in zip(rr[1:], BASIS[:k]) if c)); break
        print('a=%+d %-15s R=%d |x| = %s  decay e^%.4f  tau=%.3f  S=%+.3f  L=%s (%d dig)  %s' % (
            a, nm, R, ' '.join('%.3f' % abs(z) for z in rts), rate, tau, S, mp.nstr(Lv, 15), dig, rel or 'no relation'))
        sys.stdout.flush()
