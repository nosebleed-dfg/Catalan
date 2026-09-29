r"""The arithmetic side, exactly (2026-09-29): how much of lcm(1..n)^2 a G-family can save, prime by prime.

For a Beukers family the cumulative denominators of F = f (g - L) in any integral coordinate (T = q + O(q^2) with integral
inverse, f integral with f(0) = 1) are exactly those of the Eichler integral g = sum_m c(m) m^{-2} q^m.  So they are set by
the coefficients c(m) of W alone.  A prime p > sqrt(n) first enters at the least index k p with v_p(c(kp)) < 2.
For Eisenstein series that are multiplicative at p, c(kp) = sum_i coef_i A_i(p) A_i(k), and A_i(p) is a unit mod p:
   E_G  = E^{chi_-4,1}:  a(p)  = p^2 + chi_-4(p) = chi_-4(p)  (mod p)
   E_z  = E^{1,chi_-4}:  b(p)  = 1 + chi_-4(p) p^2 = 1
   E_-8 = E^{1,chi_-8}:  b8(p) = 1 + chi_-8(p) p^2 = 1
Weight-3 forms on Gamma_1(8) compatible with case E's harmless group <T, L_8> (vanishing at oo, no oo at the cusp 0):
   W_E = E_G(t) - 8 E_G(2t) = -E_G(t + 1/2)   (the only G-part, by the layer theorem)
   E_z(t), E_z(2t), E_-8                   (constant terms -1/4, -1/4, -3/2 must cancel)
   (E^{chi_-8,1} has an oo at the cusp 0, and the CM cusp form brings L(g,2): both are excluded.)
L-values at s = 2: W_E -> G/2, E_z(dt) -> pi^2/(12 d^2), E_-8 -> pi^2/6.  Pure G (no pi^2) leaves span{W_E, Z} with
   Z = E_z(t) + 8 E_z(2t) - (3/2) E_-8.
Families measured: W_E (case E); W*+ = W_E + 2Z (pure G, split primes saved at the first digit); W*- = W_E - 2Z (inert saved);
W_E -+ (E_z(t) - E_z(2t)) (limit G/2 -+ pi^2/16); W_B = W_E - E_z(t)/5 + 5 E_z(2t) - (4/5) E_-8 (limit G/2 - 11 pi^2/240,
split saved at p and 2p).
Part 1: cumulative denominators of g up to M (exact, by prime).  Part 2: the Gamma_1(8) family itself in the Hauptmodul
T = (1 - theta3(q^2)/theta3(q))/2 with f = theta3(q)^2: u_n, v_n exactly, v_n/u_n -> L, and the primes in denom(v_n).
Usage: python arithmetic_saving.py [M] [NMAX]"""
import sys, math
from fractions import Fraction
import numpy as np
import mpmath as mp
import sympy as sp

M = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 100

# ---------- coefficients by sieve ----------
n = np.arange(M + 1, dtype=np.int64)
chi4 = np.zeros(M + 1, dtype=np.int64); chi4[n % 4 == 1] = 1; chi4[n % 4 == 3] = -1
chi8 = np.zeros(M + 1, dtype=np.int64); chi8[(n % 8 == 1) | (n % 8 == 3)] = 1; chi8[(n % 8 == 5) | (n % 8 == 7)] = -1
a = np.zeros(M + 1, dtype=np.int64); b = np.zeros(M + 1, dtype=np.int64); b8 = np.zeros(M + 1, dtype=np.int64)
for d in range(1, M + 1):
    K = M//d; d2 = d*d
    a[d::d] += d2*chi4[1:K + 1]          # sum_{d | n} chi(n/d) d^2
    b[d::d] += d2*chi4[d]                # sum_{d | n} chi(d) d^2
    b8[d::d] += d2*chi8[d]
half = lambda arr: np.concatenate(([0], [arr[k//2] if k % 2 == 0 else 0 for k in range(1, M + 1)])).astype(np.int64)
bh = half(b)
sgn = np.where(n % 2 == 0, -1, 1).astype(np.int64)        # -(-1)^n
WE = sgn*a
Wz = np.where(n % 2 == 1, b, 0).astype(np.int64)         # E_z(t) - E_z(2t) = odd part of E_z
twoZ = 2*b + 16*bh - 3*b8
fams = {"W_E (case E)": WE, "W*+ = W_E + 2Z": WE + twoZ, "W*- = W_E - 2Z": WE - twoZ,
        "W_E - Wz (G/2 - pi^2/16)": WE - Wz, "W_E + Wz (G/2 + pi^2/16)": WE + Wz,
        "5 W_B (G/2 - 11pi^2/240)": 5*WE - b + 25*bh - 4*b8}

print("0. Checks of the forms (exact rationals)")
cterm = {"E_z": Fraction(-1, 4), "E_-8": Fraction(-3, 2)}
Zc = cterm["E_z"] + 8*cterm["E_z"] - Fraction(3, 2)*cterm["E_-8"]
Zpi = Fraction(1, 12) + Fraction(8, 48) - Fraction(3, 2)*Fraction(1, 6)
Bc = -Fraction(1, 5)*cterm["E_z"] + 5*cterm["E_z"] - Fraction(4, 5)*cterm["E_-8"]
Bpi = -Fraction(1, 5)*Fraction(1, 12) + 5*Fraction(1, 48) - Fraction(4, 5)*Fraction(1, 6)
print("   Z: constant term %s, pi^2-coefficient of L(Z,2) %s;   W_B: constant term %s, pi^2-coefficient %s"
      % (Zc, Zpi, Bc, Bpi))
for p in (5, 13, 17, 29, 3, 7, 11, 19):
    c = WE[p] + twoZ[p]
    print("   p = %2d (chi_-4 = %+d, chi_-8 = %+d):  W*+(p) = %d = %s p^2,  W*+(2p) mod p = %d"
          % (p, chi4[p], chi8[p], c, Fraction(int(c), p*p), (WE[2*p] + twoZ[2*p]) % p))

print("\n1. Cumulative denominators of g = sum c(m) m^-2 q^m for m <= %d" % M)
primes = list(sp.primerange(2, M + 1))
def vp(x, p):
    x = int(x)
    if x == 0: return 99
    k = 0
    while x % p == 0: x //= p; k += 1
    return k
full = sum(2*int(math.log(M, p) + 1e-9)*math.log(p) for p in primes)      # log lcm(1..M)^2
for name, c in fams.items():
    logD = 0.0; first = {}
    for p in primes:
        e = 0; k0 = None
        for m in range(p, M + 1, p):
            ex = 2*vp(m, p) - vp(c[m], p)
            if ex > e: e = ex
            if ex > 0 and k0 is None: k0 = m//p
        logD += e*math.log(p); first[p] = k0
    # the digit rule for large primes, by class
    big = [p for p in primes if p > 50 and p < M//4]
    rule = {}
    for cls in (1, 3):
        ks = [first[p] for p in big if p % 4 == cls]
        vals = sorted(set(ks), key=lambda x: (x is None, x))
        rule[cls] = {k: ks.count(k) for k in vals}
    print("   %-26s log D / M = %.4f   log D / log lcm^2 = %.4f   first digit k_p (50 < p < M/4): p=1 mod 4 %s, p=3 mod 4 %s"
          % (name, logD/M, logD/full, rule[1], rule[3]))

print("\n2. The Gamma_1(8) family: T = (1 - th3(q^2)/th3(q))/2, f = th3(q)^2, n <= %d" % NMAX)
P = NMAX + 1
def mul(x, y):
    r = [0]*P
    for i, xi in enumerate(x):
        if xi:
            for j in range(P - i): r[i + j] += xi*y[j]
    return r
def inv(x):
    r = [0]*P; r[0] = Fraction(1, 1)/x[0]
    for k in range(1, P): r[k] = -sum(x[i]*r[k - i] for i in range(1, k + 1))/x[0]
    return r
def th3(mm):
    s = [0]*P; k = 0
    while mm*k*k < P: s[mm*k*k] += 1 if k == 0 else 2; k += 1
    return s
t1, t2 = th3(1), th3(2)
ratio = mul(t2, [int(x) for x in inv(t1)])
T = [(1 if i == 0 else 0) - ratio[i] for i in range(P)]; T = [x//2 for x in T]
# re-expansion in T: powers T(q)^k (small integers), then peel off X = sum_n X_n T(q)^n from the q-series (O(P^2) per series)
Tpows = [[1] + [0]*(P - 1)]
for k in range(1, P): Tpows.append(mul(Tpows[-1], T))
def expand_in(X):
    res = list(X); out = [0]*P
    for k in range(P):
        ck = res[k]; out[k] = ck
        if ck:
            row = Tpows[k]
            for j in range(k, P): res[j] -= ck*row[j]
    return out
print("   T = %s ..." % (T[1:7],))
fq = mul(t1, t1)                              # theta3(q)^2 as a q-series
f = expand_in(fq)
Lc = 1
for k in range(1, P): Lc = Lc*k//math.gcd(Lc, k)
Lc2 = Lc*Lc
mp.mp.dps = 80
G = mp.catalan
checkpoints = [k for k in (60, 100, 150, 200, 250, 300, 400) if k <= NMAX]
buckets = [(1, 1.5), (1.5, 2), (2, 3), (3, 4), (4, 6), (6, 10)]
for name in ("W_E (case E)", "W*+ = W_E + 2Z", "W*- = W_E - 2Z"):
    c = fams[name]
    gq = [0] + [int(c[m])*(Lc2//(m*m)) for m in range(1, P)]     # Lc2 * Eichler integral, as a q-series
    fg = expand_in(mul(fq, gq))
    v = [Fraction(x, Lc2) for x in fg]
    NN = NMAX
    lim = mp.mpf(v[NN].numerator)/v[NN].denominator/f[NN]
    lim2 = mp.mpf(v[NN - 10].numerator)/v[NN - 10].denominator/f[NN - 10]
    print("   %-16s v_n/u_n at n = %d: %s ;  minus G/2: %s (n-10: %s)" % (name, NN, mp.nstr(lim, 25), mp.nstr(lim - G/2, 5), mp.nstr(lim2 - G/2, 5)))
    for k in checkpoints:
        den = v[k].denominator
        lk = 1
        for j in range(1, k + 1): lk = lk*j//math.gcd(lk, j)
        top = [(p, vp(den, p)) for p in sp.primerange(k//2 + 1, k + 1)]
        print("      n = %3d: log denom(v_n)/n = %.4f  (lcm^2: %.4f);  exponents of primes in (n/2, n]: split %s | inert %s"
              % (k, math.log(den)/k, 2*math.log(lk)/k, [e for p, e in top if p % 4 == 1], [e for p, e in top if p % 4 == 3]))
    # exponent profile at n = NMAX, by x = n/p and class, averaged over n in [NMAX-20, NMAX]
    prof = {}
    for k in range(NMAX - 20, NMAX + 1):
        den = v[k].denominator
        for p in sp.primerange(k//10 + 1, k + 1):
            x = k/p
            for lo, hi in buckets:
                if lo <= x < hi:
                    key = (lo, hi, 'split' if p % 4 == 1 else 'inert' if p % 4 == 3 else '2')
                    s0, c0 = prof.get(key, (0, 0)); prof[key] = (s0 + vp(den, p), c0 + 1)
    line = []
    for lo, hi in buckets:
        for cls in ('split', 'inert'):
            s0, c0 = prof.get((lo, hi, cls), (0, 0))
            if c0: line.append("x in [%g,%g) %s %.2f" % (lo, hi, cls, s0/c0))
    print("      mean exponent by x = n/p (n in [%d,%d]): %s" % (NMAX - 20, NMAX, "; ".join(line)))
