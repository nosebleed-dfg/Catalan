"""apery_p2.py — does the Apéry-like digit law hold modulo p^2?
  p^(wL) b_n  =?  b_{n_L} prod_{i<L} a_{n_i}   (mod p^2)        and, for comparison, Lucas for a_n mod p^2:  a_n =? prod a_{n_i} (mod p^2)
Families: Apéry zeta(3) (w = 3), Apéry zeta(2) (w = 2), Zudilin's F~_n (w = 7, a_n = u~_n/30, b_n = v~_n = A0).
Reports, per family: the distribution of v_p(LHS - RHS) (>= 1 always by the mod-p law), split by L = 1 / L >= 2 and by the last digit n_0 = 0 / != 0.
usage: python apery_p2.py N_apery N_zudilin
"""
import sys, json
from fractions import Fraction
from collections import Counter, defaultdict
from sympy import primerange
import apery_digit_law as ADL

def vp(x, p):
    return ADL.vp(x, p)

def digits(n, p):
    return ADL.digits(n, p)

def families(Na, Nz):
    S = ADL.seqs(Na)
    fam = [('Apery zeta(3)', 3, S[3][0], S[3][1]), ('Apery zeta(2)', 2, S[2][0], S[2][1])]
    rows = json.load(open(r"C:\Users\PC\Desktop\catalan\apery_lucas_180.json"))[:Nz + 1]
    def fr(s):
        if '/' in s: a, b = s.split('/'); return Fraction(int(a), int(b))
        return Fraction(int(s))
    U = [fr(r['A7'])/30 for r in rows]; B = [fr(r['A0']) for r in rows]
    fam.append(('Zudilin F~ (z7,z5)', 7, U, B))
    return fam

def run(Na, Nz):
    for name, w, A, B in families(Na, Nz):
        N = len(A) - 1
        stat_b = defaultdict(Counter); stat_a = defaultdict(Counter)
        for p in primerange(5, N + 1):
            for n in range(p, N + 1):
                d = digits(n, p); L = len(d) - 1
                key = ('L=1' if L == 1 else 'L>=2', 'n0=0' if d[0] == 0 else 'n0!=0')
                rhs = Fraction(B[d[L]])
                for di in d[:L]: rhs *= A[di]
                lhs = Fraction(B[n])*Fraction(p)**(w*L)
                stat_b[key][min(vp(lhs - rhs, p), 6)] += 1
                ra = Fraction(1)
                for di in d: ra *= A[di]
                stat_a[key][min(vp(Fraction(A[n]) - ra, p), 6)] += 1
        print('%s (p >= 5, n <= %d): v_p(LHS - RHS), 6 means >= 6' % (name, N))
        for key in sorted(stat_b):
            cb = stat_b[key]; ca = stat_a[key]
            tb = sum(cb.values())
            print('   %-12s digit law b_n: %s  -> mod p^2 holds in %d/%d | Lucas a_n: %s -> mod p^2 holds in %d/%d' % (
                key, dict(sorted(cb.items())), sum(v for k, v in cb.items() if k >= 2), tb,
                dict(sorted(ca.items())), sum(v for k, v in ca.items() if k >= 2), sum(ca.values())))

if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]))
