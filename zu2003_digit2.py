"""zu2003_digit2.py — the second p-adic digit of Zudilin's F~_n at single-digit primes (p^2 > 3n+2).
The local mirror kills the p^-8 level of A_0; the exponent is 7 unless the p^-7 digit X(n,p) = (p^7 A_0) mod p also vanishes.
  compute N1 N2 OUT.json      exact A_0 for n in [N1, N2]; for each single-digit prime p: v_p(A_0), X(n,p) (None if v_p(A_0) > -7)
  analyse OUT.json [...]      (i) rank-1 test of X(ap+b, p) = F_p(a) G_p(b) over F_p; (ii) is X(ap+b, p) = R(a,b) mod p for one rational R(a,b)
                              independent of p (CRT + rational reconstruction over the primes with that (a,b))?
"""
import sys, json, math
from multiprocessing import Pool
from sympy import primerange
from sympy.ntheory.modular import crt
import zu2003 as Z3

def one(n):
    Z = Z3.build('B', n)
    A0 = Z.A0
    out = []
    for p in primerange(3, n + 1):
        if p*p <= 3*n + 2: continue
        v = Z3.vq(A0, p)
        vz = min(Z3.vq(a, p) for a in Z.A.values())
        X = None
        if v == -7:
            num, den = int(A0.p), int(A0.q)
            den_u = den // p**7
            X = (num*pow(den_u, -1, p)) % p
        elif v > -7:
            X = 0
        a, b = divmod(n, p)
        out.append(dict(n=n, p=p, a=a, b=b, v=v, vzeta=vz, X=X))
    return out

def compute(N1, N2, path):
    res = []
    with Pool(14) as pool:
        for R in pool.imap_unordered(one, range(N1, N2 + 1)):
            res += R
            print('n=%d done (%d rows)' % (R[0]['n'] if R else -1, len(R)), flush=True)
    json.dump(res, open(path, 'w'))

def ratrec(x, M):
    """rational reconstruction of x mod M: r/s with |r|, s < sqrt(M/2)."""
    import math as _m
    B = _m.isqrt(M//2)
    r0, r1 = M, x % M; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q*r1
        s0, s1 = s1, s0 - q*s1
    if s1 == 0 or abs(s1) > B: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return (r1, s1)

def arc_nonempty(p, b):
    arc = [jp for jp in range(0, b + 1) if 2*b - p < jp < p - b and 2*jp != b]
    fixed = (b % 2 == 0 and b + b//2 < p and 2*b - b//2 < p)
    return bool(arc or fixed)

def analyse(paths):
    rows = []
    for pth in paths: rows += json.load(open(pth))
    rows = [r for r in rows if arc_nonempty(r['p'], r['b'])]
    print('single-digit cases with the mirror prediction 7: %d ; v_p(A0) = -7: %d ; deeper (X = 0): %d ; v < -7 (would contradict): %d' % (
        len(rows), sum(1 for r in rows if r['v'] == -7), sum(1 for r in rows if r['v'] > -7), sum(1 for r in rows if r['v'] < -7)))
    # (i) rank-1 test per p
    byp = {}
    for r in rows: byp.setdefault(r['p'], {})[(r['a'], r['b'])] = r['X']
    bad = 0; tested = 0
    for p, M in byp.items():
        keys = list(M)
        for (a1, b1) in keys:
            for (a2, b2) in keys:
                if (a1, b2) in M and (a2, b1) in M:
                    tested += 1
                    if (M[(a1, b1)]*M[(a2, b2)] - M[(a1, b2)]*M[(a2, b1)]) % p: bad += 1
    print('rank-1 test X(a1,b1)X(a2,b2) = X(a1,b2)X(a2,b1) mod p: %d quadruples, %d failures' % (tested, bad))
    # (ii) fixed rational R(a,b)?
    byab = {}
    for r in rows: byab.setdefault((r['a'], r['b']), []).append((r['p'], r['X']))
    print('fixed-rational test: for each (a,b) with >= 8 primes, reconstruct from all but the last 3 primes, check the last 3')
    for (a, b) in sorted(byab):
        L = sorted(byab[(a, b)])
        if len(L) < 8: continue
        train, test = L[:-3], L[-3:]
        M = 1
        for p, _ in train: M *= p
        x, _ = crt([p for p, _ in train], [X for _, X in train])
        rr = ratrec(int(x), M)
        if rr is None:
            print('  (a,b)=(%d,%d): no small rational from %d primes' % (a, b, len(train))); continue
        r_, s_ = rr
        ok = all((r_*pow(s_, -1, p) - X) % p == 0 for p, X in test if s_ % p)
        f = __import__('sympy').factorint
        print('  (a,b)=(%d,%d): R = %d/%d  [num %s, den %s]  predicts held-out primes %s: %s' % (
            a, b, r_, s_, f(abs(r_)) if r_ else 0, f(s_), [p for p, _ in test], ok))

if __name__ == '__main__':
    if sys.argv[1] == 'compute': compute(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    else: analyse(sys.argv[2:])
