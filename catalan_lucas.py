"""catalan_lucas.py — Lucas / digit-law structure of Zudilin's Apéry-like forms for Catalan's constant (math/0201024; zudilin_ledger.py).
u_n G - v_n -> 0;  U_n = 2^{4n} u_n (integers);  V_n = 2^{4n} v_n, den(V_n) | D_{2n-1}^2 (odd part).
Tests, odd primes p:
  (1) U_n integral; Lucas U_{ap+b} = U_a U_b (mod p) with the base-p digits of n.
  (2) exponent e_p(n) = -v_p(V_n) against 2*floor(log_p(2n-1)) (D_{2n-1}^2), split by whether p <= n or n < p <= 2n-1.
  (3) leading digit X = p^e V_n mod p against candidate digit laws: (a) V_{n_L} prod_{i<L} U_{n_i} with the digits of n (L = digits of n - 1).
usage: python catalan_lucas.py N
"""
import sys
from fractions import Fraction as F
from collections import Counter
from sympy import primerange

def seqs(N):
    p_ = lambda n: 20*n*n - 8*n + 1
    q_ = lambda n: 3520*n**6 + 5632*n**5 + 2064*n**4 - 384*n**3 - 156*n**2 + 16*n + 7
    u = [F(1), F(7, 4)]; v = [F(0), F(13, 8)]
    for n in range(1, N):
        A = (2*n + 1)**2*(2*n + 2)**2*p_(n); B = q_(n); C = (2*n - 1)**2*(2*n)**2*p_(n + 1)
        u.append((B*u[n] + C*u[n - 1])/A); v.append((B*v[n] + C*v[n - 1])/A)
    U = [x*2**(4*n) for n, x in enumerate(u)]; V = [x*2**(4*n) for n, x in enumerate(v)]
    return U, V

def vp(x, p):
    x = F(x)
    if x == 0: return 10**9
    a, b, v = x.numerator, x.denominator, 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v

def modp(x, p):
    x = F(x); return x.numerator*pow(x.denominator, -1, p) % p

def digits(n, p):
    d = []
    while n: d.append(n % p); n //= p
    return d

def flog(x, p):
    L = 0; q = p
    while q <= x: L += 1; q *= p
    return L

def run(N):
    U, V = seqs(N)
    nonint = [n for n in range(N + 1) if U[n].denominator != 1]
    print('(1) U_n = 2^{4n} u_n integral for n <= %d: %s %s' % (N, not nonint, nonint[:8]))
    tot = bad = 0; badl = []
    for p in primerange(3, N + 1):
        for n in range(p, N + 1):
            rhs = 1
            for di in digits(n, p): rhs = rhs*int(U[di]) % p
            tot += 1
            if int(U[n]) % p != rhs:
                bad += 1
                if len(badl) < 10: badl.append((n, p))
    print('    Lucas U_n = prod U_{n_i} (mod p), digits of n: %d checks, %d failures %s' % (tot, bad, badl))
    # (2) exponents
    st = Counter(); ex = []
    for n in range(2, N + 1):
        for p in primerange(3, 2*n):
            e = -vp(V[n], p); pred = 2*flog(2*n - 1, p)
            key = ('p<=n' if p <= n else 'n<p<2n', e - pred)
            st[key] += 1
    print('(2) e_p(V_n) - 2 floor(log_p(2n-1)):', dict(sorted(st.items())))
    # (3) leading digit, candidate (a) with digits of n, only primes p <= n
    tot = agree = zero_ok = zero_bad = 0; worse = 0; mism = []
    for p in primerange(3, N + 1):
        for n in range(p, N + 1):
            d = digits(n, p); L = len(d) - 1
            e = -vp(V[n], p)
            if vp(V[d[L]], p) < 0: continue
            rhs = modp(V[d[L]], p)
            for di in d[:L]: rhs = rhs*int(U[di]) % p
            tot += 1
            lvl = 2*flog(2*n - 1, p)
            if e > lvl: worse += 1; continue
            if e == lvl:
                x = modp(V[n]*F(p)**lvl, p)
                if x == rhs: agree += 1
                else:
                    if len(mism) < 10: mism.append((n, p, d, x, rhs))
            else:
                if rhs == 0: zero_ok += 1
                else: zero_bad += 1
    print('(3a) digit law with the digits of n at level 2 floor(log_p(2n-1)): %d cases ; exponent above level %d ; at level: %d agree, %d disagree %s ; below level: %d with vanishing factor, %d without' % (
        tot, worse, agree, len(mism), mism[:6], zero_ok, zero_bad))

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
