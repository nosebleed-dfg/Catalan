r"""square_stuff.py (2026-10-02): the race of primes mod 4, played with squares.

1. Why 3 mod 4 leads: every odd prime square is 1 mod 4.  D(x) = pi(x;4,3) - pi(x;4,1) against (1/2) pi(sqrt x), and the
   fair race in which a prime power p^k counts 1/k (Riemann's weighting).
2. G as the race written in squares.  Over the odd primes,
       prod p^2/(p^2 - 1) = pi^2/8,   prod p^2/(p^2 - chi_-4(p)) = G,   prod p^2/(p^2 + 1) = pi^2/12.
   The bias of the partial products of G.
3. Pythagorean triples.  With T(n) the number of triples (a > b > 0, a^2 + b^2 = n^2):
       G = pi^2/12 + sum over all triples of 1/c^2,      6G/pi^2 = 1/2 + sum over primitive triples of 1/c^2.
4. Sums of squares are this program's modular forms: r_2(n) = 4 (d_1(n) - d_3(n)),
   r_6(n) = 16 sum_{d|n} chi(n/d) d^2 - 4 sum_{d|n} chi(d) d^2, and sum' 1/(a^2+b^2)^2 = 4 zeta(2) G.
5. Four slots: with x = d0 - d2, y = d1 - d3 and n* the string with d1 and d3 exchanged,
       n n* = x^2 + y^2 (mod 101),    n rev(n) = -10 (x^2 + y^2) (mod 101).
Usage: python square_stuff.py"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 25
G = float(mp.catalan); PI2 = math.pi**2
LIM = 10**8
sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(LIM**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
pr = np.nonzero(sieve)[0][1:]                    # odd primes
chi = np.where(pr % 4 == 1, 1, -1).astype(np.int64)
S = np.cumsum(chi)                               # sum of chi_-4(p), p <= x;  D(x) = -S
del sieve

def D(x):
    i = np.searchsorted(pr, x, side='right'); return -int(S[i - 1]) if i else 0
def PI(x): return int(np.searchsorted(pr, x, side='right')) + (1 if x >= 2 else 0)

print("1. the race and the prime squares")
print("    k   D(10^k) = pi(4,3) - pi(4,1)   (1/2) pi(sqrt x)   fair race (p^k counted 1/k)")
for k in range(2, 9):
    x = 10**k
    fair = D(x) - 0.5*(PI(int(x**0.5 + 1e-9)) - 1) + D(int(round(x**(1/3) - 0.5 + 1e-9)))/3 - 0.25*(PI(int(x**0.25 + 1e-9)) - 1)
    print("   %2d   %8d                     %8.1f           %8.1f" % (k, D(x), 0.5*(PI(int(x**0.5 + 1e-9)) - 1), fair))
lead1 = int(np.count_nonzero(S > 0)); tie = int(np.count_nonzero(S == 0))
first = int(pr[np.argmax(S > 0)])
print("   over the %d odd primes below 10^8: 1 mod 4 is ahead at %d of them (first at p = %d), tied at %d" % (len(pr), lead1, first, tie))

print("\n2. three products over the odd primes")
lp = pr.astype(np.float64)
t_minus = -np.log1p(-1/lp**2); t_plus = -np.log1p(1/lp**2); t_chi = -np.log1p(-chi/lp**2)
print("        x      all minus (-> pi^2/8)   by class (-> G)   all plus (-> pi^2/12)")
for k in (1, 2, 3, 4, 6, 8):
    i = np.searchsorted(pr, 10**k, side='right')
    print("   10^%d   %.10f            %.10f      %.10f" % (k, math.exp(t_minus[:i].sum()), math.exp(t_chi[:i].sum()), math.exp(t_plus[:i].sum())))
print("   limits %.10f            %.10f      %.10f" % (PI2/8, G, PI2/12))
tail = np.cumsum(t_chi[::-1])[::-1]              # tail[i] = sum over primes >= pr[i]
T = np.empty(len(pr)); T[:-1] = tail[1:]; T[-1] = 0.0          # T[i] = log G - log G_x for x = pr[i]
cut = np.searchsorted(pr, 3*10**6)
norm = T[:cut]*lp[:cut]**1.5*np.log(lp[:cut])
w = np.diff(np.log(lp[:cut + 1]))                # logarithmic weights
print("   log G - log(partial product up to x), times x^1.5 log x, for x <= 3*10^6: mean %.3f (log-weighted %.3f), std %.3f"
      % (norm.mean(), float((norm*w).sum()/w.sum()), norm.std()))
print("   partial product above G: %.1f%% of the primes x <= 3*10^6 (log-weighted %.1f%%)"
      % (100*np.count_nonzero(T[:cut] < 0)/cut, 100*float(w[T[:cut] < 0].sum()/w.sum())))

print("\n3. Pythagorean triples")
NMAX = 2*10**7
tot_prim = 0.0; cnt = 0
m = 2
while m*m + 1 <= NMAX:
    n0 = 1 if m % 2 == 0 else 2
    ns = np.arange(n0, m, 2, dtype=np.int64)
    c = m*m + ns*ns
    ns = ns[c <= NMAX]; c = c[c <= NMAX]
    if len(ns):
        g = np.gcd(m, ns); c = c[g == 1]
        tot_prim += float((1.0/c.astype(np.float64)**2).sum()); cnt += len(c)
    m += 1
tail_p = 1/(2*math.pi*NMAX)
print("   primitive triples with c <= %d: %d (N/(2 pi) = %.0f);  sum 1/c^2 = %.10f, plus tail 1/(2 pi N) = %.10f"
      % (NMAX, cnt, NMAX/(2*math.pi), tot_prim, tot_prim + tail_p))
print("   6G/pi^2 - 1/2 = %.10f" % (6*G/PI2 - 0.5))
print("   all triples: zeta(2) * that = %.10f;   G - pi^2/12 = %.10f" % ((tot_prim + tail_p)*PI2/6, G - PI2/12))
print("   first hypotenuses: 1/25 + 1/100 + 1/169 + 1/225 + 1/289 + 1/400 + 2/625 = %.6f" % (1/25 + 1/100 + 1/169 + 1/225 + 1/289 + 1/400 + 2/625))

print("\n4. sums of squares")
N2 = 2000
r2 = [0]*(N2 + 1)
R = int(N2**0.5) + 1
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        if a*a + b*b <= N2: r2[a*a + b*b] += 1
c4 = lambda d: (0, 1, 0, -1)[d % 4]
ok2 = all(r2[n] == 4*sum(c4(d) for d in range(1, n + 1) if n % d == 0) for n in range(1, N2 + 1))
print("   r_2(n) = 4 (d_1 - d_3) for n <= %d: %s" % (N2, ok2))
N6 = 80
th = [0]*(N6 + 1)
for a in range(-9, 10):
    if a*a <= N6: th[a*a] += 1
def pmul(p, q):
    r = [0]*(N6 + 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if i + j > N6: break
                r[i + j] += x*y
    return r
th2 = pmul(th, th); th6 = pmul(pmul(th2, th2), th2)
ok6 = all(th6[n] == 16*sum(c4(n//d)*d*d for d in range(1, n + 1) if n % d == 0) - 4*sum(c4(d)*d*d for d in range(1, n + 1) if n % d == 0)
          for n in range(1, N6 + 1))
print("   r_6(n) = 16 sum chi(n/d) d^2 - 4 sum chi(d) d^2 for n <= %d: %s" % (N6, ok6))
NL = 3*10**6
sig = np.zeros(NL + 1, dtype=np.int64)
for d in range(1, NL + 1, 2):
    sig[d::d] += (1 if d % 4 == 1 else -1)
nn = np.arange(1, NL + 1, dtype=np.float64)
lat = 4*float((sig[1:]/nn**2).sum()) + math.pi/NL
print("   sum over the square lattice of 1/(a^2+b^2)^2 = %.8f (to n = %d, tail pi/N added);   4 zeta(2) G = %.8f" % (lat, NL, 4*PI2/6*G))

s1 = float(mp.zeta(2, mp.mpf(1)/4)/16); s3 = float(mp.zeta(2, mp.mpf(3)/4)/16)
print("   odd squares by class: sum_{n = 1 mod 4} 1/n^2 = %.10f = pi^2/16 + G/2 = %.10f;   sum_{n = 3 mod 4} 1/n^2 = %.10f = pi^2/16 - G/2 = %.10f"
      % (s1, PI2/16 + G/2, s3, PI2/16 - G/2))

print("\n4b. primes of the form n^2 + 1 (Landau's problem): the same race")
import sympy as sp
C = float(np.exp(np.log1p(-chi/(lp - 1)).sum()))
print("   C = prod over odd primes (1 - chi_-4(p)/(p - 1)) = %.6f (primes below 10^8; known value 1.372813...)" % C)
for NN in (10**3, 10**4, 10**5):
    cnt = sum(1 for n in range(1, NN + 1) if sp.isprime(n*n + 1))
    pred = 0.5*C*float(mp.li(NN))
    print("   n <= %d: n^2 + 1 prime for %d values;  (C/2) li(N) = %.0f;  without the race factor: %.0f" % (NN, cnt, pred, 0.5*float(mp.li(NN))))

print("\n5. four slots: multiplying a string by its conjugate or its reversal")
digs = lambda n: (n//1000, n//100 % 10, n//10 % 10, n % 10)            # d3 d2 d1 d0
okc = okr = ok11 = ok9 = True
for n in range(10000):
    d3, d2, d1, d0 = digs(n)
    x, y = d0 - d2, d1 - d3
    conj = 1000*d1 + 100*d2 + 10*d3 + d0
    rv = 1000*d0 + 100*d1 + 10*d2 + d3
    okc &= (n*conj - (x*x + y*y)) % 101 == 0
    okr &= (n*rv + 10*(x*x + y*y)) % 101 == 0
    ok11 &= (n*rv + (d0 - d1 + d2 - d3)**2) % 11 == 0
    ok9 &= (n*rv - (d0 + d1 + d2 + d3)**2) % 9 == 0
print("   n n* = x^2 + y^2 (mod 101): %s;   n rev(n) = -10 (x^2 + y^2) (mod 101): %s" % (okc, okr))
print("   n rev(n) = -(alternating sum)^2 (mod 11): %s;   n rev(n) = (digit sum)^2 (mod 9): %s" % (ok11, ok9))
print("   example 1234: x, y = 2, 2;  1234 * 3214 mod 101 = %d;  1234 * 4321 mod 101 = %d = -80 mod 101" % (1234*3214 % 101, 1234*4321 % 101))
