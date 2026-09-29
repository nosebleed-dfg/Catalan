"""catalan_doubled.py — the doubled-lattice Lucas law for Zudilin's Apéry-like forms for Catalan's constant (math/0201024).
Recurrence (x in (1/2)Z):  (2x+1)^2 (2x+2)^2 P(x) y(x+1) = Q(x) y(x) + (2x-1)^2 (2x)^2 P(x+1) y(x-1),
  P(x) = 20x^2 - 8x + 1,  Q(x) = 3520x^6 + 5632x^5 + 2064x^4 - 384x^3 - 156x^2 + 16x + 7.
Integer coset: u(0) = 1, u(1) = 7/4, v(0) = 0, v(1) = 13/8  (u(n) G - v(n) -> 0).
Half coset: the coefficient (2x-1)^2 vanishes at x = 1/2, so y(1/2) alone starts it:  u(1/2) = 1  (and any second solution is proportional).
Doubled sequences on N >= 0:  U(N) = 2^{2N} u(N/2),  V(N) = 2^{2N} v(N/2)  (V on odd N: c * U(N), c to be found).
Tests for odd primes p, N = sum N_i p^i:
  (1) Lucas  U(N) = prod U(N_i)  (mod p)  for every N (even and odd).
  (2) digit law  p^{2L} V(N) = V(N_L) prod_{i<L} U(N_i)  (mod p),  L = #digits(N) - 1,  N = 2n even (the forms);
      top digit N_L even: parameter-free; top digit odd: needs V(N_L) = c U(N_L); c reconstructed across primes.
usage: python catalan_doubled.py NMAX
"""
import sys, math
from fractions import Fraction as F
from collections import Counter, defaultdict
from sympy import primerange, factorint
from sympy.ntheory.modular import crt

def P(x): return 20*x*x - 8*x + 1
def Q(x): return 3520*x**6 + 5632*x**5 + 2064*x**4 - 384*x**3 - 156*x**2 + 16*x + 7

def solve(y0, y1, x0, count):
    """y at x0, x0+1, ..., via the recurrence; y1 = None means the degenerate start at x0 = 1/2."""
    ys = [F(y0)]
    if y1 is None:
        x = F(x0)
        ys.append(Q(x)*ys[0]/((2*x + 1)**2*(2*x + 2)**2*P(x)))
    else:
        ys.append(F(y1))
    for k in range(1, count - 1):
        x = F(x0) + k
        ys.append((Q(x)*ys[k] + (2*x - 1)**2*(2*x)**2*P(x + 1)*ys[k - 1])/((2*x + 1)**2*(2*x + 2)**2*P(x)))
    return ys

def build(NMAX):
    half = NMAX//2 + 2
    ui = solve(1, F(7, 4), 0, half); vi = solve(0, F(13, 8), 0, half)
    uh = solve(1, None, F(1, 2), half)
    U = []; V = []
    for N in range(NMAX + 1):
        if N % 2 == 0: U.append(ui[N//2]*2**(2*N)); V.append(vi[N//2]*2**(2*N))
        else: U.append(uh[N//2]*2**(2*N)); V.append(None)
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

def ratrec(x, M):
    B = math.isqrt(M//2); r0, r1 = M, x % M; s0, s1 = 0, 1
    while r1 > B:
        q = r0//r1; r0, r1 = r1, r0 - q*r1; s0, s1 = s1, s0 - q*s1
    if s1 == 0 or abs(s1) > B: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return r1, s1

def run(NMAX):
    U, V = build(NMAX)
    nonint = [N for N in range(NMAX + 1) if U[N].denominator % 2 == 1 and U[N].denominator != 1]
    print('U(N) = 2^{2N} u(N/2): odd-denominator entries for N <= %d: %s ; U(0..9) = %s' % (NMAX, nonint[:6], [str(U[N]) for N in range(10)]))
    # (1) Lucas on every N
    tot = Counter(); bad = Counter(); ex = []
    for p in primerange(3, NMAX + 1):
        for N in range(p, NMAX + 1):
            if any(vp(U[d], p) < 0 for d in digits(N, p)) or vp(U[N], p) < 0: continue
            rhs = 1
            for d in digits(N, p): rhs = rhs*modp(U[d], p) % p
            key = 'N even' if N % 2 == 0 else 'N odd'
            tot[key] += 1
            if modp(U[N], p) != rhs:
                bad[key] += 1
                if len(ex) < 8: ex.append((N, p, digits(N, p)))
    print('(1) Lucas U(N) = prod U(N_i) mod p (odd p <= N): %s checks, %s failures %s' % (dict(tot), dict(bad), ex))
    # (2) digit law for V on N = 2n
    tot_even = agree = mism = deeper_ok = deeper_bad = worse = 0; odd_top = defaultdict(list)
    for p in primerange(3, NMAX + 1):
        for N in range(2*((p + 1)//2), NMAX + 1, 2):
            d = digits(N, p); L = len(d) - 1
            if L < 1: continue
            if any(vp(U[x], p) < 0 for x in d[:L]): continue
            low = 1
            for x in d[:L]: low = low*modp(U[x], p) % p
            e = -vp(V[N], p)
            if d[L] % 2 == 0:
                if vp(V[d[L]], p) < 0: continue
                rhs = modp(V[d[L]], p)*low % p
                tot_even += 1
                if e > 2*L: worse += 1
                elif e == 2*L:
                    if modp(V[N]*F(p)**(2*L), p) == rhs: agree += 1
                    else: mism += 1
                else:
                    if rhs == 0: deeper_ok += 1
                    else: deeper_bad += 1
            else:
                # V(N_L) = c U(N_L):  c = p^{2L} V(N) / (U(N_L) * low)  mod p
                den = modp(U[d[L]], p)*low % p
                if e <= 2*L and den:
                    c = modp(V[N]*F(p)**(2*L), p)*pow(den, -1, p) % p
                    odd_top[p].append(c)
    print('(2) digit law, even top digit: %d cases ; exponent above 2L: %d ; at 2L: %d agree, %d mismatch ; below 2L: %d predicted by a vanishing factor, %d not' % (
        tot_even, worse, agree, mism, deeper_ok, deeper_bad))
    consistent = {p: len(set(cs)) == 1 for p, cs in odd_top.items() if cs}
    print('    odd top digit: c = p^{2L}V(N)/(U(N_L) prod U(N_i)) mod p constant within each p: %d/%d primes' % (sum(consistent.values()), len(consistent)))
    ps = [p for p in sorted(odd_top) if odd_top[p] and consistent.get(p)]
    if len(ps) > 6:
        train, test = ps[:-4], ps[-4:]
        M = 1
        for p in train: M *= p
        x, _ = crt(train, [odd_top[p][0] for p in train])
        rr = ratrec(int(x), M)
        if rr:
            a, b = rr
            ok = all((a*pow(b, -1, p) - odd_top[p][0]) % p == 0 for p in test)
            print('    c reconstructed across primes: %d/%d (num %s, den %s) ; predicts held-out primes %s: %s' % (a, b, factorint(abs(a)) if a else 0, factorint(b), test, ok))
        else:
            # try a character twist chi_{-4}(p)
            x2, _ = crt(train, [odd_top[p][0]*(1 if p % 4 == 1 else -1) % p for p in train])
            rr2 = ratrec(int(x2), M)
            print('    c: no p-independent rational ; with chi_{-4}(p) twist: %s' % (rr2,))

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
