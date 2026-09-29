"""Real-G push, step 6c (2026-09-29): the amplitude law at the cusp T = 1/2 on Gamma_1(8) (reads gamma1_8_amp.json).
Claim (verified by PSLQ at 38 digits for h = T, T^2, T/(1-T), T^2/(1-T), T^3/(1-T)):
    if  Lam(h) = c_Lg L(g,2) + c_G G + c_pi2 pi^2   then   A(h) = -6 c_Lg L(g,1) - (pi/2) c_G   (and nothing from pi^2),
where A(h) = lim n F_n / 2^n is the log amplitude of F = y(h) - Lam(h) f at T = 1/2 and g = eta(t)^2 eta(2t) eta(4t) eta(8t)^2.
So a family whose limit contains G is never analytic at T = 1/2 (pi and L(g,1) are Q-independent), i.e. r <= 1/2 on Gamma_1(8)."""
import sys, json
sys.argv = [sys.argv[0], '8']
import level12_modular as L
import mpmath as mp
mp.mp.dps = 60
PC = 1200
g = [0] + L.eta_series({1: 2, 2: 1, 4: 1, 8: 2}, PC)[:PC - 1]
def lval(a, k, Nl, sv, eps=1):
    A = mp.sqrt(Nl)/(2*mp.pi); sv = mp.mpf(sv); tot = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        x = 2*mp.pi*n/mp.sqrt(Nl)
        tot += a[n]*(mp.gammainc(sv, x)/mp.gamma(sv)/mp.mpf(n)**sv + eps*A**(k - 2*sv)*mp.gammainc(k - sv, x)/mp.gamma(sv)*mp.mpf(n)**(sv - k))
    return tot
Lg1, Lg2 = lval(g, 3, 8, 1), lval(g, 3, 8, 2)
print('L(g,1) =', mp.nstr(Lg1, 45)); print('L(g,2) =', mp.nstr(Lg2, 45))
res = json.load(open(r'C:\Users\PC\Desktop\catalan\gamma1_8_amp.json'))
pi = mp.pi; G = +mp.catalan
mp.mp.dps = 38
for k in ('poly1', 'poly2', 'poly3', 'p11', 'p12', 'p13'):
    Lam = mp.mpf(res[k][0]); A = mp.mpf(res[k][1])
    ra = mp.pslq([A, Lg1, pi], maxcoeff=10**6, maxsteps=10**6)
    rl = mp.pslq([Lam, Lg2, G, pi**2, mp.mpf(1)], maxcoeff=10**6, maxsteps=10**7)
    line = '%-6s  [A, L(g,1), pi]: %s   [Lam, L(g,2), G, pi^2, 1]: %s' % (k, ra, rl)
    if ra and rl and ra[0] and rl[0]:
        cLg1 = -mp.mpf(ra[1])/ra[0]; cpi = -mp.mpf(ra[2])/ra[0]
        cLg = -mp.mpf(rl[1])/rl[0]; cG = -mp.mpf(rl[2])/rl[0]
        ok = abs(cLg1 + 6*cLg) < 1e-30 and abs(cpi + cG/2) < 1e-30
        line += '   law holds: %s' % ok
    else:
        line += '   (a further constant enters; not covered)'
    print(line)
