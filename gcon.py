r"""gcon.py — construct linear forms in 1 and Catalan's constant G from scratch, by linear algebra on the quarter lattice.

Space.  R(t) = P(t) / Q(t),   Q(t) = prod_{s in POLES} (t + s)^{mult(s)},   every pole s in 1/4 + Z (class A) or 3/4 + Z (class B), s > 0.
        P(t) = (2t + c)^{cen} * prod_{r in ROOTS} (t + r) * sum_j x_j (t (t + c))^j        (ROOTS on (1/4)Z; the unknowns x_j)
        With an integer centre c the reflection t -> -c - t swaps the classes (k + 1/4 <-> c - 1 - k + 3/4).
Sum.    F = sum_{t >= 0} R(t) = sum_i [A_i zeta(i,1/4) + B_i zeta(i,3/4)] - V,   A_i, B_i = class sums of the Laurent coefficients a_{i,s},
        V = sum_{i,s} a_{i,s} sum_{0 <= l < floor(s)} (l + {s})^{-i}.
        zeta(i,1/4) + zeta(i,3/4) = (4^i - 2^i) zeta(i),  zeta(i,1/4) - zeta(i,3/4) = 4^i beta(i),  so a PURE G form needs, linearly in x:
        A_i + B_i = 0 for every i  and  A_i - B_i = 0 for every i != 2;   then  F = 8 (A_2 - B_2) G - V.
Solve.  conditions (+ optional extra vanishing rows) -> exact integer kernel (fmpz_mat.nullspace on the denominator-cleared matrix);
        kernel dimension 1 = a determined weight.  Everything downstream is exact (fmpq); F = U G - V evaluated with mpmath.
Report. per n: kernel dimension, log|F|/n, primitive height log(den/content)/n, NET, the exponent of each prime in the denominator;
        the sequences U(n) (G-coefficient) and V(n) are returned for recurrence search (guess_rec.py).
Designs are functions n -> spec dict (see DESIGNS); run:  python gcon.py DESIGN n1 n2 [rec]
"""
import sys, math, json
from flint import fmpq, fmpz, fmpq_poly, fmpz_mat
from mpmath import mp, mpf, catalan, log as mlog, fabs

Q4 = fmpq(1, 4); Q34 = fmpq(3, 4); HALF = fmpq(1, 2)

# ------------------------------------------------------------------ designs (all scaled by n)
def d_qz(n, alpha=2, halfroots=False):
    """quarter analogue of the all-equal shape: poles k in [n, 2n] in both classes (order alpha), centre c = 3n+1,
    roots at t = -i, i in [1, n] u [2n+1, 3n] (integers; halfroots=True puts them at -i-1/2 instead: i in [0,n-1] u [2n+1,3n])."""
    c = 3*n + 1
    poles = {}
    for k in range(n, 2*n + 1):
        poles[fmpq(k) + Q4] = alpha; poles[fmpq(k) + Q34] = alpha
    if halfroots:
        roots = [fmpq(i) + HALF for i in list(range(0, n)) + list(range(2*n + 1, 3*n + 1))]
    else:
        roots = [fmpq(i) for i in list(range(1, n + 1)) + list(range(2*n + 1, 3*n + 1))]
    return dict(poles=poles, roots=roots, c=c, cen=1, e=0 if alpha <= 2 else 1)

def d_low(n, alpha=2, span=1.0, gap=0.0):
    """poles k in [G, G+N] (N = span*n, G = gap*n) both classes, integer roots filling [1, G] and the mirror, centre c = 2G+N+1."""
    N = int(round(span*n)); G = int(round(gap*n)); c = 2*G + N + 1
    poles = {}
    for k in range(G, G + N + 1):
        poles[fmpq(k) + Q4] = alpha; poles[fmpq(k) + Q34] = alpha
    roots = [fmpq(i) for i in list(range(1, G + 1)) + list(range(G + N + 1, 2*G + N + 1))]
    return dict(poles=poles, roots=roots, c=c, cen=1, e=0 if alpha <= 2 else 1)

def d_apery(n, nu=1.0, mu=1.0, alpha=2, e=None):
    """Apery-style layout: poles k in [0, N] in both classes (order alpha), zeros INSIDE the summation range at t = 0..M-1,
    mirror zeros at t = -c - i (c = N + 1 swaps the classes), centre factor; free symmetric part of degree e (default: 0, or 1 if alpha = 3)."""
    N = max(1, int(round(nu*n))); M = int(round(mu*n)); c = N + 1
    poles = {}
    for k in range(0, N + 1):
        poles[fmpq(k) + Q4] = alpha; poles[fmpq(k) + Q34] = alpha
    roots = [fmpq(-i) for i in range(M)] + [fmpq(c + i) for i in range(M)]
    if e is None: e = 0 if alpha <= 2 else 1
    return dict(poles=poles, roots=roots, c=c, cen=1, e=e)

def d_half(n, nu=1.0, mu=0.5, alpha=2, e=None):
    """HALF-INTEGER centre c = N + 1/2: t -> -c - t maps each quarter class to itself (1/4: k -> N-k ; 3/4: k -> N-1-k) and integers to
    half-integers.  Poles: class 1/4 k in [0,N], class 3/4 k in [0,N-1] (order alpha); zeros at t = 0..M-1 (inside the summation range)
    mirrored to t = -(N + 1/2 + i); centre factor (2t + c) sits on the self-mirror pole."""
    N = max(1, int(round(nu*n))); M = int(round(mu*n)); c = fmpq(2*N + 1, 2)
    poles = {}
    for k in range(0, N + 1): poles[fmpq(k) + Q4] = alpha
    for k in range(0, N): poles[fmpq(k) + Q34] = alpha
    roots = [fmpq(-i) for i in range(M)] + [c + i for i in range(M)]
    if e is None: e = 0 if alpha <= 2 else 1
    return dict(poles=poles, roots=roots, c=c, cen=1, e=e)

DESIGNS = {
    'half_1_05_2': lambda n: d_half(n, 1.0, 0.5, 2), 'half_1_1_2': lambda n: d_half(n, 1.0, 1.0, 2), 'half_1_025_2': lambda n: d_half(n, 1.0, 0.25, 2),
    'half_1_05_3': lambda n: d_half(n, 1.0, 0.5, 3), 'half_1_1_3': lambda n: d_half(n, 1.0, 1.0, 3),
    'ap_1_1_2': lambda n: d_apery(n, 1.0, 1.0, 2), 'ap_1_05_2': lambda n: d_apery(n, 1.0, 0.5, 2), 'ap_1_15_2': lambda n: d_apery(n, 1.0, 1.5, 2),
    'ap_1_1_3': lambda n: d_apery(n, 1.0, 1.0, 3), 'ap_1_2_3': lambda n: d_apery(n, 1.0, 2.0, 3), 'ap_1_15_3': lambda n: d_apery(n, 1.0, 1.5, 3),
    'ap_1_1_1': lambda n: d_apery(n, 1.0, 1.0, 1), 'ap_1_05_1': lambda n: d_apery(n, 1.0, 0.5, 1),
    'qz2': lambda n: d_qz(n, 2), 'qz3': lambda n: d_qz(n, 3), 'qz2h': lambda n: d_qz(n, 2, True), 'qz3h': lambda n: d_qz(n, 3, True),
    'low2': lambda n: d_low(n, 2, 1.0, 0.0), 'low3': lambda n: d_low(n, 3, 1.0, 0.0),
    'gap2': lambda n: d_low(n, 2, 1.0, 0.5), 'gap3': lambda n: d_low(n, 3, 1.0, 0.5),
}

# ------------------------------------------------------------------ partial fractions of one basis numerator
def laurent_other(poles, s, S):
    """series of prod_{s' != s} (t + s')^{-m'} at t = -s + e, to order S-1 (list of fmpq)."""
    # log-series: -sum m' log(1 + e/(s'-s)) ;  power sums p_q = sum m' (s'-s)^{-q}
    const = fmpq(1); ps = [fmpq(0)]*S
    for s2, m2 in poles.items():
        if s2 == s: continue
        d = s2 - s
        const /= d**m2
        inv = 1/d; pw = fmpq(1)
        for q in range(1, S):
            pw *= inv; ps[q] += m2*pw
    # exp of  sum_q (-1)^q p_q e^q / q  * (-1)  ... log(1+e/d)^{-m} = -m sum_q (-1)^{q+1} (e/d)^q / q
    lg = [fmpq(0)]*S
    for q in range(1, S): lg[q] = -((-1)**(q + 1))*ps[q]/q
    ex = [fmpq(0)]*S; ex[0] = fmpq(1)
    for k in range(1, S):                      # ex' = lg' ex  ->  k ex_k = sum_{q=1}^{k} q lg_q ex_{k-q}
        acc = fmpq(0)
        for q in range(1, k + 1): acc += q*lg[q]*ex[k - q]
        ex[k] = acc/k
    return [const*x for x in ex]

def taylor(poly, x0, S):
    """[poly(x0), poly'(x0), poly''(x0)/2, ...] up to order S-1."""
    out = []; p = poly
    for m in range(S):
        out.append(p(x0)/math.factorial(m)); p = p.derivative()
    return out

def build(spec):
    poles, roots, c, cen, e = spec['poles'], spec['roots'], spec['c'], spec['cen'], spec['e']
    t = fmpq_poly([0, 1])
    base = fmpq_poly([1])
    if cen: base *= (2*t + c)
    for r in roots: base *= (t + r)
    u = t*(t + c)
    basis = []; b = base
    for j in range(e + 1):
        basis.append(b); b = b*u
    degQ = sum(poles.values())
    assert basis[-1].degree() <= degQ - 2, 'degree condition fails: deg P = %d, deg Q = %d' % (basis[-1].degree(), degQ)
    maxS = max(poles.values())
    other = {s: laurent_other(poles, s, S) for s, S in poles.items()}
    # per basis element: class sums A_i, B_i and V
    rows = []
    for bj in basis:
        A = [fmpq(0)]*(maxS + 1); B = [fmpq(0)]*(maxS + 1); V = fmpq(0)
        for s, S in poles.items():
            tay = taylor(bj, -s, S); oth = other[s]
            ser = [sum((tay[a]*oth[m - a] for a in range(m + 1)), fmpq(0)) for m in range(S)]
            k = int(s.p) // int(s.q); frac = s - k
            for i in range(1, S + 1):
                a_is = ser[S - i]
                if a_is == 0: continue
                if frac == Q4: A[i] += a_is
                else: B[i] += a_is
                # partial sum  sum_{l<k} (l + frac)^{-i}
                ps = spec.setdefault('_ps', {})
                key = (i, k, frac)
                if key not in ps:
                    acc = fmpq(0)
                    for l in range(k): acc += 1/(fmpq(l) + frac)**i
                    ps[key] = acc
                V += a_is*ps[key]
        rows.append((A, B, V))
    return basis, rows, maxS

def top_coeffs(spec, basis, which):
    """top-order Laurent coefficients a_{S,s} of each basis numerator at the poles in `which` (list of s): matrix [pole][basis]."""
    poles = spec['poles']
    out = []
    for s in which:
        S = poles[s]
        oth = laurent_other(poles, s, 1)[0]            # value of the other factors at -s
        out.append([bj(-s)*oth for bj in basis])
    return out

def kernel_modular(M, ncol, nprimes=None):
    """integer kernel of an integer matrix M (list of rows) by nmod_mat nullspaces mod 62-bit primes + CRT + rational reconstruction,
    normalised to a primitive integer vector (kernel assumed one-dimensional); verified exactly."""
    from flint import nmod_mat
    from sympy import nextprime
    from sympy.ntheory.modular import crt
    if not M:
        if ncol == 1: return [1], 0
        raise ValueError('no conditions and %d unknowns: kernel is not one-dimensional' % ncol)
    P = (1 << 62); mods, vecs = [], []
    tries = 0
    idx = None; sol = None
    while sol is None:
        P = nextprime(P)
        A = nmod_mat([[int(x) % P for x in r] for r in M], P)
        K, nul = A.nullspace()
        tries += 1
        if nul != 1:
            if tries > 6: raise ValueError('kernel dimension %d mod several primes (not one-dimensional)' % nul)
            continue
        v = [int(K[i, 0]) for i in range(ncol)]
        if idx is None: idx = max(i for i in range(ncol) if v[i])
        if v[idx] == 0: continue
        inv = pow(v[idx], -1, P); v = [x*inv % P for x in v]
        mods.append(P); vecs.append(v)
        Mod = 1
        for m in mods: Mod *= m
        cand = []
        ok = True
        for j in range(ncol):
            x, _ = crt(mods, [w[j] for w in vecs])
            rr = _ratrec(int(x), Mod)
            if rr is None: ok = False; break
            cand.append(rr)
        if not ok: continue
        L = 1
        for a, b in cand: L = L*b//math.gcd(L, b)
        ints = [a*(L//b) for a, b in cand]
        g = 0
        for x in ints: g = math.gcd(g, x)
        ints = [x//g for x in ints]
        if all(sum(r[j]*ints[j] for j in range(ncol)) == 0 for r in M): sol = ints
    return sol, len(mods)

def _ratrec(x, M):
    B = math.isqrt(M//2); r0, r1 = M, x % M; s0, s1 = 0, 1
    while r1 > B:
        q = r0//r1; r0, r1 = r1, r0 - q*r1; s0, s1 = s1, s0 - q*s1
    if s1 == 0 or abs(s1) > B: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return r1, s1

def solve_kill(spec, kill):
    """pure-G conditions + vanishing of the top-order coefficient at the poles in `kill`; modular kernel."""
    basis, rows, maxS = build(spec)
    nb = len(basis)
    conds = []
    for i in range(1, maxS + 1):
        conds.append([rows[j][0][i] + rows[j][1][i] for j in range(nb)])
        if i != 2: conds.append([rows[j][0][i] - rows[j][1][i] for j in range(nb)])
    conds += top_coeffs(spec, basis, kill)
    M = []
    for r in conds:
        L = 1
        for x in r: L = L*int(x.q)//math.gcd(L, int(x.q))
        row = [int(x*L) for x in r]
        if any(row): M.append(row)
    x, used = kernel_modular(M, nb)
    return basis, rows, [x], used

def measure_kill(spec, n, kill):
    basis, rows, ker, used = solve_kill(spec, kill)
    x = ker[0]
    A2 = sum((x[j]*rows[j][0][2] for j in range(len(x))), fmpq(0)); B2 = sum((x[j]*rows[j][1][2] for j in range(len(x))), fmpq(0))
    U = 8*(A2 - B2); V = sum((x[j]*rows[j][2] for j in range(len(x))), fmpq(0))
    if U == 0: return dict(n=n, degenerate=True), None
    den = int(U.q)*int(V.q)//math.gcd(int(U.q), int(V.q))
    g = math.gcd(int(U.p)*(den//int(U.q)), int(V.p)*(den//int(V.q)))
    digits = len(str(abs(int(U.p)))) + len(str(int(U.q))) + len(str(abs(int(V.p)))) + len(str(int(V.q)))
    mp.dps = digits + 60
    F = mpf(int(U.p))/int(U.q)*catalan - mpf(int(V.p))/int(V.q)
    lf = float(mlog(fabs(F)))/n
    h = (math.log(den) - (math.log(g) if g > 1 else 0))/n
    return dict(n=n, logF=lf, height=h, net=h + lf, primes_used=used), (U, V)

def solve(spec):
    basis, rows, maxS = build(spec)
    nb = len(basis)
    conds = []
    for i in range(1, maxS + 1):
        conds.append([rows[j][0][i] + rows[j][1][i] for j in range(nb)])            # zeta(i) coefficient
        if i != 2: conds.append([rows[j][0][i] - rows[j][1][i] for j in range(nb)])  # beta(i) coefficient
    # clear denominators row by row, exact integer kernel
    M = []
    for r in conds:
        L = 1
        for x in r: L = L*int(x.q)//math.gcd(L, int(x.q))
        M.append([int(x*L) for x in r])
    if all(all(v == 0 for v in r) for r in M):
        ker = [[1 if j == jj else 0 for j in range(nb)] for jj in range(nb)]
    else:
        K, nullity = fmpz_mat(M).nullspace()
        ker = [[int(K[i, j]) for i in range(nb)] for j in range(nullity)]
    return basis, rows, ker

def measure(spec, n):
    basis, rows, ker = solve(spec)
    out = dict(n=n, kdim=len(ker))
    if len(ker) != 1: return out, None
    x = ker[0]
    A2 = sum((x[j]*rows[j][0][2] for j in range(len(x))), fmpq(0)); B2 = sum((x[j]*rows[j][1][2] for j in range(len(x))), fmpq(0))
    U = 8*(A2 - B2); V = sum((x[j]*rows[j][2] for j in range(len(x))), fmpq(0))
    if U == 0:
        out['degenerate'] = 'G-coefficient vanishes: the sum is the rational number %s' % ('0' if V == 0 else 'V')
        return out, None
    den = int(U.q)*int(V.q)//math.gcd(int(U.q), int(V.q))
    g = math.gcd(int(U.p)*(den//int(U.q)), int(V.p)*(den//int(V.q)))
    digits = len(str(abs(int(U.p)))) + len(str(int(U.q))) + len(str(abs(int(V.p)))) + len(str(int(V.q)))
    mp.dps = digits + 60
    F = mpf(int(U.p))/int(U.q)*catalan - mpf(int(V.p))/int(V.q)
    lf = float(mlog(fabs(F)))/n if F != 0 else float('-inf')
    h = (math.log(den) - (math.log(g) if g > 1 else 0))/n
    out.update(logF=lf, height=h, net=h + lf, U=U, V=V, den=den, content=g)
    if n <= 6:
        # numeric check: sum_{t>=0} P(t)/Q(t) directly
        from mpmath import nsum, inf
        P = sum((basis[j]*x[j] for j in range(len(x))), fmpq_poly([0]))
        coeffs = [mpf(int(cq.p))/int(cq.q) for cq in P.coeffs()]
        def Rn(tt):
            num = sum(cf*tt**k for k, cf in enumerate(coeffs))
            den_ = mpf(1)
            for s, m in spec['poles'].items(): den_ *= (tt + mpf(int(s.p))/int(s.q))**m
            return num/den_
        t_start = 0
        while t_start < 10*n + 10 and Rn(mpf(t_start)) == 0: t_start += 1
        Fnum = nsum(Rn, [t_start, inf])
        out['check'] = float(fabs(Fnum - F)/fabs(F))
    return out, (U, V)

def prime_exps(U, V, P):
    den = int(U.q)*int(V.q)//math.gcd(int(U.q), int(V.q))
    g = math.gcd(int(U.p)*(den//int(U.q)), int(V.p)*(den//int(V.q)))
    res = {}
    from sympy import primerange
    for p in primerange(2, P):
        e = 0; d = den
        while d % p == 0: d //= p; e += 1
        cc = 0; gg = g
        while gg and gg % p == 0: gg //= p; cc += 1
        res[p] = e - cc
    return res

def run(name, n1, n2, rec=False):
    seqU, seqV = [], []
    print('design %s' % name)
    for n in range(n1, n2 + 1):
        spec = DESIGNS[name](n)
        out, UV = measure(spec, n)
        if UV is None:
            print('  n=%2d: kernel dimension %d (not a determined weight)' % (n, out['kdim'])); continue
        U, V = UV
        ex = prime_exps(U, V, 8*n + 8)
        top = {lab: sorted(set(ex[p] for p in ex if lo < p <= hi)) for lab, lo, hi in (('(n,2n]', n, 2*n), ('(2n,4n]', 2*n, 4*n), ('(4n,8n]', 4*n, 8*n + 8))}
        print('  n=%2d: kernel 1 ; log|F|/n = %+.4f ; height = %.4f ; NET = %+.4f ; 2-adic %d ; odd exponents by range %s%s' % (
            n, out['logF'], out['height'], out['net'], ex.get(2, 0), top, (' ; numeric check rel.err %.1e' % out['check']) if 'check' in out else ''))
        seqU.append(U); seqV.append(V)
    if rec and seqU:
        import guess_rec as GR
        L = 1
        for u in seqU: L = L*int(u.q)//math.gcd(L, int(u.q))
        print('  U(n) denominators lcm = %d' % L)
    return seqU, seqV

if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), 'rec' in sys.argv)
