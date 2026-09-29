r"""twoadic_law.py — the exact 2-adic structure of the twisted Catalan Hankel machine (2026-09-28, real date).

Machine: twisted_hankel.py 'twist' (weight pi t sinh(pi t)/cosh(pi t)^2, poles u = -a^2, a = 0..K-1), at X = 0 where the 2-adic content sits
(e2(K) = -v2(det A)); X = 0 is the 2-adic reading of the same kernel sum (the full-period partial sums beta_{2^N} -> 0 2-adically).

Ingredient laws (exact, verified):
  pole data      F(a) = -4(-1)^a beta_a :  v2 of Mahler coefficients Delta^n F(0) = n + 1                         (n = 1..63)
                 divided differences over consecutive square nodes z = -a^2:
                 v2([z_r, ..., z_{r+m}] F~) = 1 + s2(m)  for every start r and every m >= 1                     (r <= 16, m <= 23)
                 v2(F(r)) = 2 (r odd), 3 + 2 v2(r) (r even)
  moment data    (2e+1)|E_2e| : v2 of Mahler coefficients = n                                                  (n = 0..40)
                 omega(r, m) = phi0(prod_{l=r}^{r+m-1} (u + l^2)):  v2 = -(2m+1) for every start r
  nodes          natural order 0, 1, 2, ... is a 2-ordering of the squares; 2-sequence = v2((2k)!/2) (Bhargava's factorial of squares)
Structure (exact, verified K = 4..56, 60, 64 by the 2-adic Smith normal form; consistent with every e2 up to K = 384):
  the K elementary divisors of A over Z_2 are
    moment species   eps_n = -(2K-3) + 12 n - 4 s2(n),  n = 0 .. N-1,   N = N(K) = #{n : 12n - 4 s2(n) < 2K - 3}   (all odd, < 0)
    ones             K - 2N divisors equal to 1
    partners         N divisors >= 2 (2, 3 mostly; 4, 5, 6 near the transitions N -> N+1)
  so  e2(K) = sum_{n<N} [(2K-3) - 12n + 4 s2(n)] - (K - 2N) - (sum of partners).
  N ~ K/6 gives K^2/6 (the constant), 4 sum_{n<N} s2(n) ~ (K/3) log2(K/6) (the carry term, Delange's digit sum);
  partners = 3 on K = 3*2^j gives e2 = (3/2)4^j + (j-3)2^j exactly.  Open: the partner law (average 2.5 at 5*2^j, 3 at 3*2^j, -> 4 at 2^j).
usage:
  python twoadic_law.py dd N                 divided-difference table of the pole data over square nodes
  python twoadic_law.py ed K1,K2,...         2-adic elementary divisors of A and the check against the law
  python twoadic_law.py formula K1,K2,...    N(K), the moment/ones part, and (if e2 is in hankel_2adic_e2.json) the implied partner sum
"""
import sys, io, contextlib, math, json, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def s2(n): return bin(n).count('1')
def v2(x):
    x = F(x); a, b = x.numerator, x.denominator
    if a == 0: return None
    return ((abs(a) & -abs(a)).bit_length() - 1) - ((b & -b).bit_length() - 1)
def v2i(n):
    n = abs(n); return (n & -n).bit_length() - 1 if n else None

def moment_species(K):
    eps = []; n = 0
    while -(2*K - 3) + 12*n - 4*s2(n) < 0:
        eps.append(-(2*K - 3) + 12*n - 4*s2(n)); n += 1
    return eps

def base(K):
    eps = moment_species(K); N = len(eps)
    return N, -sum(eps) - (K - 2*N)

def elementary_divisors(K, machine='twist'):
    import twisted_hankel as TH, flint
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        A, B, h, labels, head, tail = TH.entries(machine, K)
    L = 1
    for x in A[:2*K - 1]: L = L*x.denominator//math.gcd(L, x.denominator)
    s = v2i(L)
    M = flint.fmpz_mat([[int(A[i + j]*L) for j in range(K)] for i in range(K)])
    S = M.snf()
    return sorted(v2i(int(S[i, i])) - s for i in range(K))

def dd_table(N):
    beta = [F(0)]
    for k in range(N + 2): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
    Fv = [-4*(-1)**a*beta[a] for a in range(N + 1)]; z = [-a*a for a in range(N + 1)]
    dd = {(r, 0): Fv[r] for r in range(N + 1)}
    for m in range(1, N + 1):
        for r in range(N + 1 - m): dd[(r, m)] = (dd[(r + 1, m - 1)] - dd[(r, m - 1)])/(z[r + m] - z[r])
    bad = sum(1 for (r, m), x in dd.items() if m >= 1 and v2(x) != 1 + s2(m))
    return dd, bad

def unified_snf(machine, K):
    """Unified v-frame law (2026-09-28): the 2-adic Smith form of A^v_K is the sorted union of
       pole band  kappa(K,m) + 4m - 2K + c,  m < K - M*   (kappa = 1 + s2(K-1-m) - s2(m); c = 0 integer poles, c = -1 half-integer poles)
       moment pivots h_q, q < M*  (twisted CDH(1/2,1/2,1/2): 8q - 4 s2(q) - 1;  Catalan CDH(1,1,0): 4q + 3 v2(q!) + v2((q+1)!) - 1)
    with M* the minimiser of the total; at a tie (two minimisers) the top divisor is raised (content: by exactly 1).
    Returns (list of minimisers, predicted multiset for the first minimiser, tropical minimum)."""
    vf = lambda n: n - s2(n)
    c = 0 if machine == 'twist' else -1
    pole = lambda m: 1 + s2(K - 1 - m) - s2(m) + 4*m - 2*K + c
    mom = (lambda q: 8*q - 4*s2(q) - 1) if machine == 'twist' else (lambda q: 4*q + 3*vf(q) + vf(q + 1) - 1)
    pc = [0]*(K + 1)
    for m in range(K): pc[m + 1] = pc[m] + pole(m)
    tots = []; acc = 0
    for M in range(K + 1):
        tots.append(pc[K - M] + acc)
        if M < K: acc += mom(M)
    mn = min(tots); Ms = [M for M, t in enumerate(tots) if t == mn]
    return Ms, sorted([pole(m) for m in range(K - Ms[0])] + [mom(q) for q in range(Ms[0])]), mn

def catalan_content_law(K):
    """Catalan machine (v-scale): e2 = -min_M [sum_{k<K-M} pi(K,k) + sum_{q<M} h_q] - [tie], pi = kappa + 4k - 2K - 1.
    Exact for all 134 known K (3..128 content, 160..384 det A)."""
    Ms, snf, mn = unified_snf('catalan', K)
    return -mn - (1 if len(Ms) > 1 else 0), Ms[0], len(Ms) > 1

def content_law(K):
    """THE 2-adic content law (2026-09-28): e2(K) = -min_M T(K, M) - [minimum attained twice],
    T(K, M) = sum_{n<M} eps_n + sum_{m<K-M} kappa_m,  eps_n = -(2K-3) + 12n - 4 s2(n)  (moment species, 'up 3 back 1'),
    kappa_m = 1 + s2(K-1-m) - s2(m)  (pole band, Kummer carries of m + (K-1-m)).
    Exact for all 126 content values K = 3..128 and every det-A value K = 160..384 (twisted machine)."""
    eps = lambda n: -(2*K - 3) + 12*n - 4*s2(n)
    kap = lambda m: 1 + s2(K - 1 - m) - s2(m)
    vals = []; acc_e = 0
    kum = [0]*(K + 1)
    for m in range(K): kum[m + 1] = kum[m] + kap(m)
    for M in range(0, K + 1):
        vals.append(acc_e + kum[K - M])
        if M < K: acc_e += eps(M)
    mn = min(vals); ties = vals.count(mn)
    return -mn - (1 if ties > 1 else 0), vals.index(mn), ties > 1

def cdh_check(nmax=30):
    """Both G machines' weights are continuous dual Hahn in u = t^2: twisted (1/2,1/2,1/2): b_n = 2n^2+2n+3/4, lambda_n = n^4, h_n = (n!)^4/2;
    Catalan (1,1,0): b_n = (2n+1)(n+1), lambda_n = n^3(n+1), h_n = (n!)^3 (n+1)!/2.  Checks orthogonality and norms against the exact moments."""
    from sympy import euler, bernoulli
    mu_t = [F((2*e + 1)*abs(int(euler(2*e))), 2**(2*e + 1)) for e in range(2*nmax + 2)]
    mu_c = [(2**(2*e + 2) - 1)*abs(F(int(bernoulli(2*e + 2).p), int(bernoulli(2*e + 2).q))) for e in range(2*nmax + 2)]
    fams = (('twisted CDH(1/2,1/2,1/2)', mu_t, lambda n: F(2*n*n + 2*n) + F(3, 4), lambda n: F(n)**4, lambda n: F(math.factorial(n))**4/2),
            ('Catalan CDH(1,1,0)', mu_c, lambda n: F((2*n + 1)*(n + 1)), lambda n: F(n**3*(n + 1)), lambda n: F(math.factorial(n))**3*math.factorial(n + 1)/2))
    for name, mu, b, lam, h in fams:
        P = [[F(1)], [-b(0), F(1)]]
        for k in range(1, nmax):
            nxt = [F(0)] + P[k]
            for i, c in enumerate(P[k]): nxt[i] -= b(k)*c
            for i, c in enumerate(P[k - 1]): nxt[i] -= lam(k)*c
            P.append(nxt)
        phi = lambda p: sum(c*mu[e] for e, c in enumerate(p))
        def mul(p, q):
            r = [F(0)]*(len(p) + len(q) - 1)
            for i, a in enumerate(p):
                for j, c in enumerate(q): r[i + j] += a*c
            return r
        bad = sum(1 for m in range(nmax) for n in range(m) if phi(mul(P[m], P[n])) != 0) + sum(1 for m in range(nmax) if phi(mul(P[m], P[m])) != h(m))
        print('%s (u = t^2): orthogonality + norms for n < %d: %d violations' % (name, nmax, bad))

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'cdh':
        cdh_check(int(sys.argv[2]) if len(sys.argv) > 2 else 30); sys.exit()
    if cmd == 'claw':
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hankel_2adic_e2.json')
        data = {int(k): v for k, v in json.load(open(path))['e2']['catalan'].items()}
        miss = [(K, data[K], catalan_content_law(K)[0]) for K in sorted(data) if K >= 3 and catalan_content_law(K)[0] != data[K]]
        print('Catalan content law vs hankel_2adic_e2.json: %d K checked, misses %s' % (len([K for K in data if K >= 3]), miss))
        sys.exit()
    if cmd == 'merge':
        # usage: python twoadic_law.py merge twist|catalan K1,K2,...   (exact 2-adic Smith form vs the unified union)
        import twisted_hankel as TH, flint
        mach = sys.argv[2]
        for K in [int(x) for x in sys.argv[3].split(',')]:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                A, B, h, labels, head, tail = TH.entries(mach, K)
            rows = [[A[i + j]*(F(4)**(i + j - K) if mach == 'twist' else 1) for j in range(K)] for i in range(K)]
            L = 1
            for r in rows:
                for x in r: L = L*x.denominator//math.gcd(L, x.denominator)
            s = v2i(L)
            S = flint.fmpz_mat([[int(x*L) for x in r] for r in rows]).snf()
            snf = sorted(v2i(int(S[i, i])) - s for i in range(K))
            Ms, pred, mn = unified_snf(mach, K)
            tie = len(Ms) > 1
            ok = snf == pred or (tie and snf[:-1] == pred[:-1] and snf[-1] > pred[-1])
            print('%s K=%d: M*=%s %s  SNF == union%s' % (mach, K, Ms, 'OK' if ok else 'FAIL', ' (tie: top divisor +%d)' % (snf[-1] - pred[-1]) if tie else ''))
        sys.exit()
    if cmd == 'law':
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hankel_2adic_e2.json')
        data = {int(k): v for k, v in json.load(open(path))['e2']['twist'].items()} if os.path.exists(path) else {}
        Ks = sorted(k for k in data if k >= 3) if sys.argv[2] == 'all' else [int(x) for x in sys.argv[2].split(',')]
        miss = []
        for K in Ks:
            e, Mstar, tie = content_law(K)
            if K in data and data[K] != e: miss.append((K, data[K], e))
            if sys.argv[2] != 'all': print('K=%d: e2 = %d (M* = %d%s)%s' % (K, e, Mstar, ', tie' if tie else '', ('  data %d' % data[K]) if K in data else ''))
        if sys.argv[2] == 'all': print('content law vs hankel_2adic_e2.json: %d K checked, misses %s' % (len(Ks), miss))
        sys.exit()
    if cmd == 'dd':
        N = int(sys.argv[2]); dd, bad = dd_table(N)
        print('divided differences over square nodes, N = %d: violations of v2 = 1 + s2(m): %d of %d' % (N, bad, sum(1 for k in dd if k[1] >= 1)))
        for r in range(min(N, 12)):
            print('r=%2d:' % r, ' '.join('%3s' % v2(dd[(r, m)]) for m in range(0, min(24, N + 1 - r))))
    elif cmd == 'ed':
        for K in [int(x) for x in sys.argv[2].split(',')]:
            ed = elementary_divisors(K)
            neg = [e for e in ed if e < 0]; ones = ed.count(1); partners = [e for e in ed if e >= 2]
            ok = neg == moment_species(K) and ones == K - 2*len(neg) and len(partners) == len(neg)
            print('K=%d: e2 = %d | N = %d, moment species %s, ones %d (K-2N = %d), partners %s | law %s' % (
                K, -sum(ed), len(neg), 'match' if neg == moment_species(K) else neg, ones, K - 2*len(neg), partners, 'OK' if ok else 'FAIL'), flush=True)
    elif cmd == 'formula':
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hankel_2adic_e2.json')
        e2 = {int(k): v for k, v in json.load(open(path))['e2']['twist'].items()} if os.path.exists(path) else {}
        for K in [int(x) for x in sys.argv[2].split(',')]:
            N, b = base(K)
            if K in e2: print('K=%d: N=%d, moment+ones part %d, e2 %d, implied partner sum %d (%.2f per partner)' % (K, N, b, e2[K], b - e2[K], (b - e2[K])/N))
            else: print('K=%d: N=%d, moment+ones part %d; e2 = %d - (partner sum, N values in [2, 6])' % (K, N, b, b))
