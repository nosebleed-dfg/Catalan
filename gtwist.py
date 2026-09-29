r"""gtwist.py — Catalan's character in the SUMMATION, poles on the INTEGER lattice.
F = sum_{t >= t0} chi(t + sigma) R(t),  chi = chi_{-4},  R = P/Q with poles at t = -k (k in [K0, K1], order S), P = (2t + c)^cen prod (t + r) sum_j x_j (t(t+c))^j.
Partial fractions R = sum c_{i,k} (t+k)^{-i}.  With d = (sigma - k) mod 4:
  sum_{t>=t0} chi(t+sigma)(t+k)^{-i} = L_d(i) - sum_{u=1}^{k+t0-1} chi(u + sigma - k) u^{-i},
  L_0 = beta(i), L_2 = -beta(i), L_1 = -2^{-i} eta(i), L_3 = +2^{-i} eta(i),   eta(i) = (1 - 2^{1-i}) zeta(i), eta(1) = log 2.
Pure G (linear in x): beta-coefficient Bt_i = 0 (i != 2) and eta-coefficient Et_i = 0 (all i);  then F = Bt_2 G - V.
Everything exact; numeric check by direct summation; denominators per prime.
usage: python gtwist.py DESIGN n1 n2      (designs in DESIGNS)
"""
import sys, math
from flint import fmpq, fmpz, fmpq_poly, fmpz_mat
from mpmath import mp, mpf, catalan, log as mlog, fabs, nsum, inf
from sympy import primerange
import gcon

def chi(u):
    u %= 4
    return 1 if u == 1 else (-1 if u == 3 else 0)

def design_eq(n, S=2, e=1, sigma=None, gap=1.0, span=1.0):
    """all-equal palindromic weight: poles k in [G+1, G+1+N] (N = span n, G = gap n) order S, integer roots t = -1..-G and the mirror
    -(c-1)..-(c-G), centre c = 2G + N + 2 ; free symmetric degree e ; sum from t0 = -G (the roots absorb t in [-G, -1])."""
    N = max(1, int(round(span*n))); G = int(round(gap*n)); c = 2*G + N + 2
    poles = {fmpq(k): S for k in range(G + 1, G + N + 2)}
    roots = [fmpq(i) for i in range(1, G + 1)] + [fmpq(c - i) for i in range(1, G + 1)]
    if sigma is None: sigma = 0
    return dict(poles=poles, roots=roots, c=c, cen=1, e=e, sigma=sigma, t0=-G)

def design_ap(n, nu=1.0, mu=1.0, S=2, e=1, sigma=0):
    """Apery layout, twisted: poles t = -k, k in [1, N] (order S); zeros INSIDE the summation range at t = 0..M-1 and the mirror zeros
    t = -c - i (c = N + 1, so the pole set is symmetric under k -> c - k); centre factor; free symmetric degree e; sum from t0 = 0."""
    N = max(1, int(round(nu*n))); M = int(round(mu*n)); c = N + 1
    poles = {fmpq(k): S for k in range(1, N + 1)}
    roots = [fmpq(-i) for i in range(M)] + [fmpq(c + i) for i in range(M)]
    return dict(poles=poles, roots=roots, c=c, cen=1, e=e, sigma=sigma, t0=0)

DESIGNS = {}
for S in (2, 3):
    for e in (1, 2, 3):
        for sg in (0, 1, 2, 3):
            DESIGNS['eq_S%d_e%d_s%d' % (S, e, sg)] = (lambda S, e, sg: (lambda n: design_eq(n, S, e, sg)))(S, e, sg)

def design_gen(n, nu=1.0, mu=1.0, S=2, sigma=0, extra=0):
    """NO symmetry: poles t = -k, k in [1, N] (order S); zeros inside the summation range t = 0..M-1; free general part
    sum_j x_j (t - M)^j of the minimal degree E = (#pure-G conditions) + extra; sum from t0 = 0."""
    N = max(1, int(round(nu*n))); M = int(round(mu*n))
    poles = {fmpq(k): S for k in range(1, N + 1)}
    roots = [fmpq(-i) for i in range(M)]
    E = (2*S - 1) + extra
    return dict(poles=poles, roots=roots, c=0, cen=0, e=E, sigma=sigma, t0=0, general=True, M=M)

def build(spec):
    poles, roots, c, cen, e, sigma, t0 = spec['poles'], spec['roots'], spec['c'], spec['cen'], spec['e'], spec['sigma'], spec['t0']
    t = fmpq_poly([0, 1]); base = fmpq_poly([1])
    if cen: base *= (2*t + c)
    for r in roots: base *= (t + r)
    basis = []; b = base
    if spec.get('general'):
        step = t - spec['M']
    else:
        step = t*(t + c)
    for j in range(e + 1): basis.append(b); b = b*step
    degQ = sum(poles.values())
    assert basis[-1].degree() <= degQ - 1, 'degree condition fails (deg P %d, deg Q %d)' % (basis[-1].degree(), degQ)
    maxS = max(poles.values())
    other = {s: gcon.laurent_other(poles, s, S) for s, S in poles.items()}
    ps = {}
    rows = []
    for bj in basis:
        Bt = [fmpq(0)]*(maxS + 2); Et = [fmpq(0)]*(maxS + 2); V = fmpq(0)
        for s, S in poles.items():
            k = int(s)
            tay = gcon.taylor(bj, -s, S); oth = other[s]
            ser = [sum((tay[a]*oth[m - a] for a in range(m + 1)), fmpq(0)) for m in range(S)]
            d = (sigma - k) % 4
            for i in range(1, S + 1):
                cik = ser[S - i]
                if cik == 0: continue
                if d == 0: Bt[i] += cik
                elif d == 2: Bt[i] -= cik
                elif d == 1: Et[i] -= cik/fmpq(2)**i
                else: Et[i] += cik/fmpq(2)**i
                top = k + t0 - 1
                key = (i, k, d, top)
                if key not in ps:
                    acc = fmpq(0)
                    for uu in range(1, top + 1):
                        cv = chi(uu + sigma - k)
                        if cv: acc += fmpq(cv, uu**i)
                    ps[key] = acc
                V += cik*ps[key]
        rows.append((Bt, Et, V))
    return basis, rows, maxS

def solve(spec):
    basis, rows, maxS = build(spec)
    nb = len(basis); conds = []
    for i in range(1, maxS + 1):
        if i != 2: conds.append([rows[j][0][i] for j in range(nb)])
        conds.append([rows[j][1][i] for j in range(nb)])
    Mz = []
    for r in conds:
        L = 1
        for x in r: L = L*int(x.q)//math.gcd(L, int(x.q))
        row = [int(x*L) for x in r]
        if any(row): Mz.append(row)
    if not Mz:
        ker = [[1 if i == j else 0 for i in range(nb)] for j in range(nb)]
    else:
        K, nul = fmpz_mat(Mz).nullspace()
        ker = []
        for cc in range(nul):
            v = [int(K[i, cc]) for i in range(nb)]; g = 0
            for a in v: g = math.gcd(g, a)
            ker.append([a//g for a in v])
    return basis, rows, ker

def measure(spec, n, check=True):
    basis, rows, ker = solve(spec)
    if len(ker) != 1: return dict(n=n, kdim=len(ker)), None
    x = ker[0]
    U = sum((x[j]*rows[j][0][2] for j in range(len(x))), fmpq(0)); V = sum((x[j]*rows[j][2] for j in range(len(x))), fmpq(0))
    if U == 0: return dict(n=n, kdim=1, degenerate=True), None
    den = int(U.q)*int(V.q)//math.gcd(int(U.q), int(V.q))
    g = math.gcd(int(U.p)*(den//int(U.q)), int(V.p)*(den//int(V.q)))
    digits = len(str(abs(int(U.p)))) + len(str(int(U.q))) + len(str(abs(int(V.p)))) + len(str(int(V.q)))
    mp.dps = digits + 80
    F = mpf(int(U.p))/int(U.q)*catalan - mpf(int(V.p))/int(V.q)
    out = dict(n=n, kdim=1, logF=float(mlog(fabs(F)))/n, height=(math.log(den) - (math.log(g) if g > 1 else 0))/n)
    out['net'] = out['logF'] + out['height']
    if check and n <= 5:
        P = sum((basis[j]*x[j] for j in range(len(x))), fmpq_poly([0]))
        cf = [mpf(int(q.p))/int(q.q) for q in P.coeffs()]
        def term(tt):
            tt = int(tt); cv = chi(tt + spec['sigma'])
            if cv == 0: return mpf(0)
            num = sum(a*mpf(tt)**kk for kk, a in enumerate(cf)); den_ = mpf(1)
            for s, S in spec['poles'].items(): den_ *= (mpf(tt) + int(s))**S
            return cv*num/den_
        tot = mpf(0)
        for tt in range(spec['t0'], 4000): tot += term(tt)
        # tail by pairing the 4-periodic signs (alternating in effect): extrapolate with nsum over blocks of 4
        tail = nsum(lambda b: term(4000 + 4*int(b)) + term(4001 + 4*int(b)) + term(4002 + 4*int(b)) + term(4003 + 4*int(b)), [0, inf])
        out['check'] = float(fabs(tot + tail - F)/fabs(F))
    return out, (U, V)

def run(name, n1, n2):
    print('design %s' % name)
    for n in range(n1, n2 + 1):
        spec = DESIGNS[name](n)
        try:
            out, UV = measure(spec, n)
        except AssertionError as ex:
            print('  n=%2d: %s' % (n, ex)); continue
        if UV is None:
            print('  n=%2d: %s' % (n, 'degenerate (G-coefficient 0)' if out.get('degenerate') else 'kernel dimension %d' % out['kdim'])); continue
        U, V = UV
        ex = gcon.prime_exps(U, V, 4*n + 8)
        rng = {lab: sorted(set(ex[p] for p in ex if lo < p <= hi)) for lab, lo, hi in (('(n/2,n]', n//2, n), ('(n,2n]', n, 2*n), ('(2n,4n]', 2*n, 4*n + 8))}
        print('  n=%2d: logF/n %+.4f  height %.4f  NET %+.4f  2-adic %d  odd %s%s' % (n, out['logF'], out['height'], out['net'], ex.get(2, 0), rng,
              (' chk %.0e' % out['check']) if 'check' in out else ''))

if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
