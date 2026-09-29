"""Checks for the explainer paper (CATALAN_2ADIC_PAPER.md), 2026-09-29.
 (W1) finite differences of g(y) = (2y+1)^-2:  Delta^k g(y) = (-2)^k k! * (sum_{i<=k} 1/(2y+2i+1)) / prod_{i<=k} (2y+2i+1),
      and v2(sum of reciprocals of k+1 consecutive odd numbers) >= v2(k+1); hence v2(Delta^k g(y)) >= k + v2((k+1)!).
 (W2) the machine's formulas from the integral (mpmath): moments mu_k = int_0^inf (4t^2)^k w(t) dt = 4^k (2^(2k+2)-1)|B_(2k+2)|,
      pole values int_0^inf w(t)/(4t^2+(2j+1)^2) dt = F_j(G) = (-1)^j (2j+1)/2 (G - beta_j) - 1/(4(2j+1)),  w(t) = pi t^2 cosh(pi t)/sinh^2(pi t).
 (W3) CDH(1,1,0) recurrence in the v-scale: b_n = 4(2n+1)(n+1), lam_n = 16 n^3 (n+1) are the recurrence of the monic OPs of mu_k.
 (W4) C = 1/2 + sum (-1)^k (k+1) T_(2k+1) (tangent numbers) vs the Mahler-series C; binary digits of 2C.
 (W5) e2(K) - K^2/6 - (K/3) log2 K, divided by K, on the stored values (the O(K) claim)."""
import sys, math, json
from fractions import Fraction as F
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan'); sys.path.insert(0, r'C:\Users\PC\Desktop\catalan\proof_checks')
from lemlib import v2, s2, Cconst, monic_ops
from sympy import bernoulli
import mpmath as mp
# W1
bad_id = bad_sum = bad_bound = 0
for y in range(-12, 13):
    g = lambda t: F(1, (2*t + 1)**2)
    for k in range(0, 40):
        d = sum((-1)**(k - i)*math.comb(k, i)*g(y + i) for i in range(k + 1))
        odds = [2*y + 2*i + 1 for i in range(k + 1)]
        ssum = sum(F(1, o) for o in odds)
        rhs = F((-2)**k*math.factorial(k))*ssum/math.prod(odds)
        if d != rhs: bad_id += 1
        if ssum != 0 and v2(ssum) < v2(k + 1): bad_sum += 1
        if d != 0 and v2(d) < k + (k + 1 - s2(k + 1)): bad_bound += 1   # d = 0 (windows symmetric about -1/2) passes trivially
print('W1: identity violations %d, reciprocal-sum bound violations %d, v2(Delta^k g) >= k + v2((k+1)!) violations %d  (y in [-12,12], k < 40)'
      % (bad_id, bad_sum, bad_bound))
# W2
mp.mp.dps = 40
w = lambda t: mp.pi*t**2*mp.cosh(mp.pi*t)/mp.sinh(mp.pi*t)**2
G = +mp.catalan
worst_mu = worst_pole = mp.mpf(0)
for k in range(6):
    I = mp.quad(lambda t: (4*t*t)**k*w(t), [0, 1, 5, 20, 60])
    B_ = bernoulli(2*k + 2); exact = mp.mpf(4)**k*(2**(2*k + 2) - 1)*abs(mp.mpf(int(B_.p))/int(B_.q))
    worst_mu = max(worst_mu, abs(I - exact)/abs(exact))
beta = [F(0)]
for k in range(10): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
for j in range(6):
    I = mp.quad(lambda t: w(t)/(4*t*t + (2*j + 1)**2), [0, 1, 5, 20, 60])
    Fj = (-1)**j*mp.mpf(2*j + 1)/2*(G - mp.mpf(beta[j].numerator)/beta[j].denominator) - mp.mpf(1)/(4*(2*j + 1))
    worst_pole = max(worst_pole, abs(I - Fj))
print('W2: max relative error, moments k<6: %s;  max abs error, pole values j<6: %s' % (mp.nstr(worst_mu, 3), mp.nstr(worst_pole, 3)))
# W3
mu = []
for k in range(64):
    B_ = bernoulli(2*k + 2); mu.append(F(4)**k*(2**(2*k + 2) - 1)*abs(F(int(B_.p), int(B_.q))))
P, h = monic_ops(mu, 30)
bad_rec = 0
for n in range(28):
    lhs = [F(0)] + P[n]
    rhs = P[n + 1][:] + [F(0)]*(len(lhs) - len(P[n + 1]))
    for i, a in enumerate(P[n]): rhs[i] += 4*(2*n + 1)*(n + 1)*a
    if n >= 1:
        for i, a in enumerate(P[n - 1]): rhs[i] += 16*n**3*(n + 1)*a
    if lhs != rhs: bad_rec += 1
bad_h = sum(1 for n in range(30) if h[n] != F(16)**n*math.factorial(n)**3*math.factorial(n + 1)/2)
print('W3: CDH recurrence mismatches %d of 28; norm formula h_n = 16^n (n!)^3 (n+1)!/2 mismatches %d of 30' % (bad_rec, bad_h))
# W4
T = [1]   # tangent numbers T_1, T_3, ... via the Seidel/boustrophedon-free recurrence from tan' = 1 + tan^2
N = 140
tan = [F(0)]*(2*N + 2); tan[1] = F(1)                       # power series of tan x
for n in range(1, 2*N):                                      # tan' = 1 + tan^2  =>  (n+1) a_{n+1} = [n==0] + sum a_i a_{n-i}
    s = sum(tan[i]*tan[n - i] for i in range(n + 1))
    tan[n + 1] = s/(n + 1)
Tn = [int(tan[2*k + 1]*math.factorial(2*k + 1)) for k in range(N)]
print('W4: tangent numbers', Tn[:6])
Ctan = F(1, 2) + sum((-1)**k*(k + 1)*Tn[k] for k in range(N))
Cm = Cconst(220)
print('    v2(C_tangent - C_mahler) = %s  (precision ~ 2*%d)' % (v2(Ctan - Cm), N))
twoC = 2*Cm
mod = 2**40
r = (twoC.numerator*pow(twoC.denominator, -1, mod)) % mod
print('    2C mod 2^40 in binary (2-adic digits, least significant on the right): ...%s' % bin(r)[2:].zfill(40))
# W5
data = {int(k): v for k, v in json.load(open(r'C:\Users\PC\Desktop\catalan\hankel_2adic_e2.json'))['e2']['catalan'].items()}
for K in (32, 48, 64, 96, 128, 160, 192, 256, 320, 384):
    if K in data:
        print('W5: K=%3d  e2 = %6d   (e2 - K^2/6 - (K/3)log2 K)/K = %.3f' % (K, data[K], (data[K] - K*K/6 - K/3*math.log2(K))/K))
