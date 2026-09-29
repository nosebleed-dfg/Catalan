"""apery_digit_law.py — does the digit law of Zudilin's F~_n hold for Apéry's own approximations?
  zeta(3): a_n = sum_k C(n,k)^2 C(n+k,k)^2,  b_n = sum_k C(n,k)^2 C(n+k,k)^2 [ H3(n) + sum_{m<=k} (-1)^(m-1)/(2 m^3 C(n,m) C(n+m,m)) ],  D_n^3 b_n integral
  zeta(2): a_n = sum_k C(n,k)^2 C(n+k,k),    b_n = sum_k C(n,k)^2 C(n+k,k) [ 2 sum_{m<=n} (-1)^(m-1)/m^2 + sum_{m<=k} (-1)^(n+m-1)/(m^2 C(n,m) C(n+m,m)) ],  D_n^2 b_n integral
Digit law (w = 3 resp. 2): for every prime p <= n with n = sum n_i p^i (L = #digits - 1):
  v_p(b_n) >= -wL,  and  p^(wL) b_n = b_{n_L} * prod_{i<L} a_{n_i}  (mod p)   [so v_p(b_n) > -wL  iff  p divides the right side]
Also: Lucas a_n = prod a_{n_i} mod p (Gessel for zeta(3); known for zeta(2)).
usage: python apery_digit_law.py N
"""
import sys
from fractions import Fraction
from math import comb
from sympy import primerange

def seqs(N):
    out = {}
    for w in (3, 2):
        A, B = [], []
        for n in range(N + 1):
            a = 0; b = Fraction(0)
            if w == 3:
                H = sum(Fraction(1, m**3) for m in range(1, n + 1))
                for k in range(n + 1):
                    t = comb(n, k)**2*comb(n + k, k)**2
                    c = H + sum(Fraction((-1)**(m - 1), 2*m**3*comb(n, m)*comb(n + m, m)) for m in range(1, k + 1))
                    a += t; b += t*c
            else:
                H = 2*sum(Fraction((-1)**(m - 1), m**2) for m in range(1, n + 1))
                for k in range(n + 1):
                    t = comb(n, k)**2*comb(n + k, k)
                    c = H + sum(Fraction((-1)**(n + m - 1), m**2*comb(n, m)*comb(n + m, m)) for m in range(1, k + 1))
                    a += t; b += t*c
            A.append(a); B.append(b)
        out[w] = (A, B)
    return out

def vp(x, p):
    if x == 0: return 10**9
    x = Fraction(x); num, den = x.numerator, x.denominator; v = 0
    while num % p == 0: num //= p; v += 1
    while den % p == 0: den //= p; v -= 1
    return v

def modp(x, p):
    x = Fraction(x); return x.numerator*pow(x.denominator, -1, p) % p

def digits(n, p):
    d = []
    while n: d.append(n % p); n //= p
    return d

def test(N):
    S = seqs(N)
    for w, (A, B) in S.items():
        luc = bad_luc = 0; tot = worse = eq_bad = deep_ok = deep_bad = 0; ex = []
        for p in primerange(2, N + 1):
            for n in range(p, N + 1):
                d = digits(n, p); L = len(d) - 1
                rhs_a = 1
                for di in d: rhs_a = rhs_a*A[di] % p
                luc += 1
                if A[n] % p != rhs_a: bad_luc += 1
                if vp(B[d[L]], p) < 0: continue
                rhs = modp(B[d[L]], p)
                for di in d[:L]: rhs = rhs*A[di] % p
                v = vp(B[n], p); tot += 1
                if v < -w*L: worse += 1; ex.append(('worse', n, p, v, L)); continue
                if v == -w*L:
                    if modp(B[n]*Fraction(p)**(w*L), p) != rhs: eq_bad += 1; ex.append(('mismatch', n, p))
                else:
                    if rhs == 0: deep_ok += 1
                    else: deep_bad += 1; ex.append(('deep', n, p, v, L))
        print('zeta(%d), n <= %d: Lucas for a_n %d checks, %d failures ; digit law for b_n: %d cases, worse than -%dL: %d, mismatches at -%dL: %d, deeper: %d predicted / %d not   %s' % (
            w, N, luc, bad_luc, tot, w, worse, w, eq_bad, deep_ok, deep_bad, ex[:8]))

if __name__ == '__main__':
    test(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
