"""apery_lucas.py — Apéry-like (Lucas) structure of Zudilin's F~_n = u~_n z(7) + w~_n z(5) - v~_n  (zu2003 form B, r = 3, eta = (3; 1^9)).
  compute N OUT.json   exact u~_n = A7, w~_n = A5, v~_n = A0 for n = 0..N (n = 0: F~_0 = 30 z(7)), stored as strings
  test OUT.json        (1) U(n) = u~_n/30 integral?  Lucas U(ap+b) = U(a)U(b) mod p (all digits)?
                       (2) digit law: p^(7L) A0(n) = A0(n_L) prod_{i<L} U(n_i) (mod p), n = sum n_i p^i, L = #digits - 1
                       (3) U(b) = 0 mod p for 2p/3 <= b < p ?  (Phi~ from the Lucas factor)
"""
import sys, json
from multiprocessing import Pool
from flint import fmpq, fmpz
from sympy import primerange
import zu2003 as Z3

def one(n):
    if n == 0: return dict(n=0, A7='30', A5='0', A0='0')
    Z = Z3.build('B', n)
    return dict(n=n, A7=str(Z.A.get(7, fmpq(0))), A5=str(Z.A.get(5, fmpq(0))), A0=str(Z.A0))

def compute(N, path):
    res = {}
    with Pool(14) as pool:
        for R in pool.imap_unordered(one, range(0, N + 1)):
            res[R['n']] = R
    json.dump([res[n] for n in range(N + 1)], open(path, 'w'))
    print('stored n = 0..%d' % N)

def q(s):
    if '/' in s:
        a, b = s.split('/'); return fmpq(int(a), int(b))
    return fmpq(int(s))

def modp(x, p):
    """x (fmpq, p-integral) mod p"""
    num, den = int(x.p), int(x.q)
    return num*pow(den, -1, p) % p

def digits(n, p):
    d = []
    while n: d.append(n % p); n //= p
    return d

def test(path):
    rows = json.load(open(path))
    N = len(rows) - 1
    A7 = [q(r['A7']) for r in rows]; A0 = [q(r['A0']) for r in rows]; A5 = [q(r['A5']) for r in rows]
    U = [a/30 for a in A7]
    nonint = [n for n in range(N + 1) if U[n].q != 1]
    print('(1) U(n) = u~_n/30 integral for n <= %d: %s  (non-integral at %s)' % (N, not nonint, nonint[:10]))
    # Lucas for U
    bad = 0; tot = 0; badlist = []
    for p in primerange(2, N + 1):
        for n in range(p, N + 1):
            if U[n].q % p == 0: continue
            d = digits(n, p)
            rhs = 1
            for di in d: rhs = rhs*modp(U[di], p) % p
            tot += 1
            if modp(U[n], p) != rhs:
                bad += 1
                if len(badlist) < 12: badlist.append((n, p))
    print('    Lucas U(n) = prod U(n_i) mod p: %d checks (all p <= n), %d failures %s' % (tot, bad, badlist))
    # (2) digit law for A0
    tot = 0; bad = []; deeper_ok = 0; deeper_bad = []; worse = []
    for p in primerange(2, N + 1):
        for n in range(p, N + 1):
            d = digits(n, p); L = len(d) - 1
            v = Z3.vq(A0[n], p)
            rhs = modp(A0[d[L]], p) if A0[d[L]].q % p else None
            if rhs is None: continue
            for di in d[:L]:
                if U[di].q % p: rhs = rhs*modp(U[di], p) % p
                else: rhs = None; break
            if rhs is None: continue
            tot += 1
            if v < -7*L: worse.append((n, p, v, L)); continue
            if v == -7*L:
                lhs = modp(A0[n]*fmpq(p)**(7*L), p)
                if lhs != rhs: bad.append((n, p, lhs, rhs))
            else:
                if rhs == 0: deeper_ok += 1
                else: deeper_bad.append((n, p, v, L))
    print('(2) digit law p^(7L) A0(n) = A0(n_L) prod U(n_i) mod p: %d cases (all p <= n) ; exponent worse than 7L: %d %s' % (tot, len(worse), worse[:8]))
    print('    at exponent exactly 7L: %d mismatches %s' % (len(bad), bad[:8]))
    print('    deeper than 7L: %d predicted by a vanishing factor, %d not predicted %s' % (deeper_ok, len(deeper_bad), deeper_bad[:8]))
    # (3) U(b) mod p for b in [2p/3, p)
    viol = []
    for p in primerange(3, N + 1):
        for b in range(0, p):
            if 3*b >= 2*p and U[b].q % p:
                if modp(U[b], p) != 0: viol.append((p, b))
    print('(3) U(b) = 0 mod p for 2p/3 <= b < p: violations %d %s' % (len(viol), viol[:10]))
    zeros = {}
    for p in primerange(3, 60):
        zeros[p] = [b for b in range(0, p) if U[b].q % p == 0 and modp(U[b], p) == 0]
    print('    zeros of U(b) mod p, b < p (p < 60):')
    for p, z in zeros.items(): print('      p=%d: %s   (2p/3 = %.1f)' % (p, z, 2*p/3))

if __name__ == '__main__':
    if sys.argv[1] == 'compute': compute(int(sys.argv[2]), sys.argv[3])
    else: test(sys.argv[2])
