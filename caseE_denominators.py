"""Real-G push, step 1 (2026-09-29): the true denominators of Zagier's sporadic case E, a holonomic family of
rational approximations to Catalan's constant G (Calegari-Dimitrov-Tang, arXiv 2408.15403, Remark 11.1.17).

Recurrence (Zagier 2009, case E):  (n+1)^2 u_{n+1} = (12n^2 + 12n + 4) u_n - 32 n^2 u_{n-1}.
  u: u_0 = 1, u_1 = 4 (integers: u_n = sum_k C(n,k) C(2k,k) C(2n-2k,n-k));  v: v_0 = 0, v_1 = 1 (rational).
CDT: ODE singularities {0, 1/8; 1/4, inf}; denominator type [1..n]^2 (growth e^2 per step); their method would need
growth below 16*(1/4) = 4, i.e. tau := lim log(den)/n < ln 4 = 1.386.  This script measures the TRUE denominators.
Usage: python caseE_denominators.py N"""
import sys, math
from fractions import Fraction as F
import mpmath as mp
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
u = [F(1), F(4)]; v = [F(0), F(1)]
for n in range(1, N):
    a = 12*n*n + 12*n + 4; b = 32*n*n; d = (n + 1)**2
    u.append((a*u[n] - b*u[n - 1])/d); v.append((a*v[n] - b*v[n - 1])/d)
assert all(x.denominator == 1 for x in u), 'u_n should be integers'
# 1. the limit of v_n/u_n, identified
mp.mp.dps = 80
L = mp.mpf(v[N].numerator)/v[N].denominator/(mp.mpf(u[N].numerator)/u[N].denominator)
Lm = mp.mpf(v[N - 20].numerator)/v[N - 20].denominator/(mp.mpf(u[N - 20].numerator)/u[N - 20].denominator)
print('v_N/u_N = %s   (change over the last 20 steps: %s)' % (mp.nstr(L, 40), mp.nstr(abs(L - Lm), 3)))
rel = mp.pslq([L, 1, mp.catalan, mp.pi**2, mp.pi*mp.log(2), mp.log(2)**2], maxcoeff=10**6, maxsteps=10**6)
print('PSLQ relation [L, 1, G, pi^2, pi*ln2, ln2^2]:', rel)
# 2. denominators: effective growth and per-prime exponents
def primes_upto(n):
    s = bytearray([1])*(n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def vp(x, p):
    c = 0
    while x % p == 0: x //= p; c += 1
    return c
print('\n  n    ln(den v_n)/n   2*ln(lcm(1..n))/n   [ratio]')
for n in [50, 100, 150, 200, 250, 300, 400, 500, 600, 800]:
    if n > N: break
    den = v[n].denominator
    lcm = 1
    for k in range(1, n + 1): lcm = lcm*k//math.gcd(lcm, k)
    print('%4d   %8.4f          %8.4f          %.4f' % (n, math.log(den)/n, 2*math.log(lcm)/n, math.log(den)/(2*math.log(lcm))))
n = N
den = v[n].denominator
ps = primes_upto(n)
tot = {'1mod4': [0.0, 0.0], '3mod4': [0.0, 0.0], '2': [0.0, 0.0]}
short = []
for p in ps:
    e = vp(den, p); full = 2*int(math.floor(math.log(n)/math.log(p) + 1e-12))
    key = '2' if p == 2 else ('1mod4' if p % 4 == 1 else '3mod4')
    tot[key][0] += e*math.log(p); tot[key][1] += full*math.log(p)
    if e < full: short.append((p, e, full))
print('\nn = %d: log-weight of den(v_n) vs the [1..n]^2 bound, by class:' % n)
for k in tot: print('  %-6s  actual %9.2f   bound %9.2f   (%.1f%%)' % (k, tot[k][0], tot[k][1], 100*tot[k][0]/max(tot[k][1], 1e-9)))
print('primes with exponent below 2*floor(log_p n): %d of %d; first few: %s' % (len(short), len(ps), short[:25]))
# where do the short primes sit? distribution of n/p for p > sqrt(n)
big = [(p, e, f) for (p, e, f) in short if p*p > n]
print('short primes above sqrt(n): %d; n/p values (rounded): %s' % (len(big), sorted(set(round(n/p, 2) for p, e, f in big))[:40]))
