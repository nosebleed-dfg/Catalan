"""catalan_twostate_v.py — the rational part V_n = 2^{4n} v_n of Zudilin's Catalan forms under the two-state carry law.
Conjecture:  p^{2L} V(N) = Y_{c_{L-1}}(N_L) * prod_{i<L} X_{c_{i-1}}(N_i)   (mod p),   N = 2n, L = #digits(N) - 1,
  X_0, X_1 as in catalan_twostate.py; top digit: Y_0(even M) = V(M) (the form's own rational part at M/2);
  Y_1(odd M) (top digit reached with an incoming carry) inferred per prime from the data and tested for consistency,
  then reconstructed across primes.  Deeper cases (exponent < 2L) must coincide with a vanishing product.
usage: python catalan_twostate_v.py NMAX
"""
import sys
from collections import Counter, defaultdict
from sympy import primerange, factorint
from sympy.ntheory.modular import crt
import catalan_twostate as CT
import catalan_doubled as CD

def run(NMAX):
    U, V = CT.seq(NMAX)
    mp_, vpq = CT.mp_, CT.vpq
    Y1 = defaultdict(dict); incons = 0
    stats = Counter(); fails = []
    for p in primerange(5, NMAX):
        if 2*p*p + 1 > NMAX: break
        if mp_(U[1], p) == 0: continue
        X0, X1 = CT.families(U, p)
        # infer Y1(M) from two-digit N whose top digit is carried
        for N in range(p + 1, min(p*p, NMAX + 1), 2):
            Nd, cin = CT.carries_digits(N//2, p)
            if len(Nd) != 2 or not cin[1]: continue
            low = (X1 if cin[0] else X0)[Nd[0]]
            if low == 0 or -vpq(V[N], p) > 2: continue
            val = mp_(V[N]*p**2, p)*pow(low, -1, p) % p
            M = Nd[1]
            if M in Y1[p] and Y1[p][M] != val: incons += 1
            Y1[p][M] = val
        # test all even N
        for N in range(p + 1, NMAX + 1, 2):
            Nd, cin = CT.carries_digits(N//2, p)
            L = len(Nd) - 1
            low = 1
            for d, c in zip(Nd[:L], cin[:L]): low = low*(X1 if c else X0)[d] % p
            if cin[L]:
                if Nd[L] not in Y1[p]: continue
                top = Y1[p][Nd[L]]
            else:
                if vpq(V[Nd[L]], p) < 0: continue
                top = mp_(V[Nd[L]], p)
            chi = 1 if p % 4 == 1 else -1
            rhs = top*low*(chi**(L + cin[L])) % p
            e = -vpq(V[N], p)
            key = ('top carried' if cin[L] else 'top plain', L)
            if e > 2*L: stats[key + ('worse',)] += 1; fails.append((N, p, 'worse')); continue
            if e == 2*L:
                ok = mp_(V[N]*p**(2*L), p) == rhs
                stats[key + ('at 2L', ok)] += 1
                if not ok and len(fails) < 10: fails.append((N, p, Nd, cin))
            else:
                stats[key + ('deeper', rhs == 0)] += 1
    print('Y1 inference inconsistencies (same p, M, different values): %d' % incons)
    for k in sorted(stats, key=str): print('  ', k, stats[k])
    print('  first failures:', fails[:10])
    # reconstruct Y1(M) across primes for small odd M
    for M in (1, 3, 5, 7, 9):
        pts = [(p, Y1[p][M]) for p in sorted(Y1) if M in Y1[p]]
        if len(pts) < 8: continue
        train, test = pts[:-3], pts[-3:]
        Mod = 1
        for p, _ in train: Mod *= p
        x, _ = crt([p for p, _ in train], [y for _, y in train])
        rr = CD.ratrec(int(x), Mod)
        if rr:
            a, b = rr
            ok = all((a*pow(b, -1, p) - y) % p == 0 for p, y in test)
            print('  Y1(%d) = %d/%d  (num %s, den %s)  held-out ok: %s' % (M, a, b, factorint(abs(a)) if a else 0, factorint(b), ok))
        else:
            print('  Y1(%d): no small rational from %d primes' % (M, len(train)))

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 8000)
