r"""zu2003_lemma_check.py — exhaustive check of the combinatorial ingredients of the single-digit proof (results log, "Zudilin's 2003 question").
For F~_n (h0 = 3n+2, all h_j = n+1), window pole j = k - h1 in [p, n], top coefficient g_j = (-1)^n (n-2j) C(n+j,n)^3 C(2n-j,n)^3 C(n,j)^6.
For every n <= N and every prime p <= n with p^2 > 3n+2 (n = a p + b):
  (1) {j in [p,n] : v_p(g_j) = 0}  ==  {c p + j' : 1 <= c <= a, 0 <= j' <= b, 2b-p < j' < p-b, 2j' != b}      (Kummer)
  (2) g_{j*} == -g_j (mod p) for j* = c p + (b - j')                                                          (Lucas + centre factor)
  (3) the fixed point j' = b/2 (b even), when carry-free, has v_p(n - 2j) >= 1
  (4) the p-divisible distances seen from a leading pole depend only on (a, c):
      numerator (m - k)/p in {-(c+1),...,-(a+c)} (mult 3) and {a-c+1,...,2a-c} (mult 3); poles (k'-k)/p in [-c, a-c] \ {0} (mult 6).
usage: python zu2003_lemma_check.py N
"""
import sys
from sympy import primerange

def vfact(N, p):
    v = 0; q = p
    while q <= N: v += N//q; q *= p
    return v

def vbin(N, K, p):
    return vfact(N, p) - vfact(K, p) - vfact(N - K, p)

def vint(x, p):
    x = abs(x); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v

def binom_mod_p(N, K, p):
    """Lucas."""
    r = 1
    while N or K:
        n0, k0 = N % p, K % p
        if k0 > n0: return 0
        num = 1; den = 1
        for i in range(k0): num = num*(n0 - i) % p; den = den*(i + 1) % p
        r = r*num*pow(den, -1, p) % p
        N //= p; K //= p
    return r

def g_mod_p(n, j, p):
    s = -1 if n % 2 else 1
    return s*(n - 2*j)*binom_mod_p(n + j, n, p)**3*binom_mod_p(2*n - j, n, p)**3*binom_mod_p(n, j, p)**6 % p

def vg(n, j, p):
    return vint(n - 2*j, p) + 3*vbin(n + j, n, p) + 3*vbin(2*n - j, n, p) + 6*vbin(n, j, p)

def check(N):
    cases = 0; bad = []; empty_arcs = 0; stats = {'pairs': 0, 'fixed': 0}
    for n in range(2, N + 1):
        for p in primerange(3, n + 1):
            if p*p <= 3*n + 2: continue
            cases += 1
            a, b = divmod(n, p)
            pred = set()
            for c in range(1, a + 1):
                for jp in range(0, b + 1):
                    if 2*b - p < jp < p - b and 2*jp != b: pred.add(c*p + jp)
            act = {j for j in range(p, n + 1) if vg(n, j, p) == 0}
            if pred != act: bad.append(('set', n, p, sorted(pred ^ act)[:6]))
            if not act: empty_arcs += 1
            for j in act:
                c, jp = divmod(j, p)
                js = c*p + (b - jp)
                if js not in act or (g_mod_p(n, js, p) + g_mod_p(n, j, p)) % p != 0: bad.append(('mirror', n, p, j, js))
                stats['pairs'] += 1
                # (4) distance multisets
                k = n + 1 + j
                lo = sorted((m - k)//p for m in range(1, n + 1) if (m - k) % p == 0)
                hi = sorted((m - k)//p for m in range(2*n + 2, 3*n + 2) if (m - k) % p == 0)
                po = sorted((jj - j)//p for jj in range(0, n + 1) if jj != j and (jj - j) % p == 0)
                if lo != list(range(-(a + c), -c)) or hi != list(range(a - c + 1, 2*a - c + 1)) or po != [t for t in range(-c, a - c + 1) if t != 0]:
                    bad.append(('dist', n, p, j, lo, hi, po))
            if b % 2 == 0:
                for c in range(1, a + 1):
                    j = c*p + b//2
                    if b + b//2 < p and 2*b - b//2 < p:          # carry-free fixed point
                        stats['fixed'] += 1
                        if vint(n - 2*j, p) < 1: bad.append(('fixed', n, p, j))
    print('n <= %d: %d (n, p) cases with p^2 > 3n+2, p <= n ; mirror pairs checked %d ; carry-free fixed points %d ; cases with empty leading set %d' % (
        N, cases, stats['pairs'], stats['fixed'], empty_arcs))
    print('failures: %d %s' % (len(bad), bad[:10]))

if __name__ == '__main__':
    check(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
