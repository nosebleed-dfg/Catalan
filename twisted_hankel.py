r"""twisted_hankel.py — the K^2 Hankel machine with Catalan's character in the WEIGHT and the poles on the INTEGER lattice
(2026-09-28, real date).

Weight (t > 0, u = t^2):   w(t) = pi t sinh(pi t)/cosh(pi t)^2 = (2t^2/pi) sum_{n>=0} (-1)^n (2n+1)/(t^2 + (n+1/2)^2)^2,
i.e. double poles of the kernel at t = i(n+1/2) with signs (-1)^n = chi_{-4}(2n+1): the character lives on the kernel's lattice.
  moments      mu(e) = (2e+1)|E_{2e}|/2^{2e+1}          (Euler numbers: von Staudt-free, only the prime 2 in the denominator)
  pole values  phi(1/(u+a^2)) = sum_{n>=0} (-1)^n/(n+a+1/2)^2 = 4(-1)^a (G - beta_a),   beta_a = sum_{k<a} (-1)^k/(2k+1)^2,
               for EVERY integer a >= 0 (a = 0 gives 4G); no rational part.
Both identities verified to 40 digits by quadrature (mu(0..5), a = 0..6).
Decay e^{-pi t}, pole spacing 1: rho = pi, so the closeness is U(pi) = -ln 9 by the rho-law, the same as the Catalan machine; the change is
all on the arithmetic side (integer poles instead of half-integer poles, Euler instead of Genocchi moments).

Machines (same engine):
  twist    : the weight above, poles a = a0, a0+step, ... (K of them).
  catalan  : control, identical to profiles.py 'catalan' (poles at b = 2j+1 in v = 4t^2, alternating kernel on the integers).
  zeta2int, zeta2half : the two zeta(2) machines of profiles.py ('zeta2int', 'zeta2'), as 2-adic controls (target pi^2/6, pi^2/2).
Determinant polynomial: det(A + XB) = det(B) * charpoly(-B^{-1}A)(X) (one solve + one characteristic polynomial; identical to the
h+1 exact determinants + Newton interpolation of profiles.py, checked at K = 20, 40, 60; 15x faster at K = 60). detpoly_interp keeps
the old route.
Entries  phi(u^s Num(u)/D(u)),  D = prod_{tail}(u + a^2),  Num = prod_{head}(u + a^2)^copies * Z(u)^r,
Z(u) = prod_{n<M} (4u + (2n+1)^2) (zeros on the kernel's lattice, i.e. on the character's support; twist only).
Hankel pencil det(A + XB) of size h = K - N, primitive part P, value at X = G (ball arithmetic), full content factorisation.
usage: python twisted_hankel.py MACHINE K [a0=0] [N=0] [copies=0] [M=0] [r=1] [step=1]
writes twisted_<machine>_K.._a.._N.._c.._M.._r.._s...json next to this file.
"""
import sys, math, time, json, os; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, euler, primerange, factorint
import flint

def polymul(a, b):
    r = [0]*(len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b): r[i + j] += ai*bj
    return r

def ev(p, v):
    s = 0
    for c in reversed(p): s = s*v + c
    return s

def polydivmod_monic(num, den):
    num = list(num); q = [0]*max(0, len(num) - len(den) + 1)
    for k in range(len(num) - len(den), -1, -1):
        c = num[k + len(den) - 1]; q[k] = c
        if c:
            for i, d in enumerate(den): num[k + i] -= c*d
    r = num[:len(den) - 1]; r += [0]*(len(den) - 1 - len(r))
    return q, r

def setup(machine, K, a0=0, N=0, copies=0, M=0, r=1, step=1):
    """returns (poles-in-u list of 'square' values, head, tail, mu(e), polevalue(label)->(A,B), numerator poly, label->square)."""
    if machine == 'twist':
        labels = [a0 + step*j for j in range(K)]
        sq = lambda a: a*a
        def mu(e): return F((2*e + 1)*abs(int(euler(2*e))), 2**(2*e + 1))
        beta = {}; acc = F(0)
        for k in range(max(labels) + 1):
            beta[k] = acc; acc += F((-1)**k, (2*k + 1)**2)
        def polevalue(a): return (F(-4*(-1)**a)*beta[a], F(4*(-1)**a))
        extra = [1]
        if M:
            Z = [1]
            for n in range(M): Z = polymul(Z, [(2*n + 1)**2, 4])
            for _ in range(r): extra = polymul(extra, Z)
    elif machine == 'catalan':
        labels = [2*j + 1 for j in range(K)]
        sq = lambda b: b*b
        def mu(e): return F(4)**e*(2**(2*e + 2) - 1)*abs(F(int(bernoulli(2*e + 2).p), int(bernoulli(2*e + 2).q)))
        betaalt = {}; acc = F(0); sgn = 1
        for b in range(1, max(labels) + 2, 2): betaalt[b] = acc; acc += F(sgn, b*b); sgn = -sgn
        def polevalue(b):
            sign = (-1)**((b - 1)//2)
            return (F(-sign*b, 2)*betaalt[b] - F(1, 4*b), F(sign*b, 2))
        extra = [1]
        assert M == 0
    elif machine == 'zeta2int':
        labels = [j + 1 for j in range(K)]
        sq = lambda j: j*j
        def mu(e): return abs(F(int(bernoulli(2*e + 2).p), int(bernoulli(2*e + 2).q)))/2
        H2 = {0: F(0)}
        for j in range(1, K + 1): H2[j] = H2[j - 1] + F(1, j*j)
        def polevalue(j): return (-F(j, 2)*H2[j - 1] - F(1, 2) - F(1, 4*j), F(j, 2))
        extra = [1]
        assert M == 0
    elif machine == 'zeta2half':
        labels = [2*j + 1 for j in range(K)]
        sq = lambda b: b*b
        def mu(e): return F(4)**e*abs(F(int(bernoulli(2*e + 2).p), int(bernoulli(2*e + 2).q)))/2
        beta = {}; acc = F(0)
        for b in range(1, max(labels) + 2, 2): beta[b] = acc; acc += F(1, b*b)
        def polevalue(b): return (-F(b, 4)*beta[b] - F(1, 8*b) - F(1, 8), F(b, 16))
        extra = [1]
        assert M == 0
    else:
        raise SystemExit('machine must be twist | catalan | zeta2int | zeta2half')
    head = labels[:N]; tail = labels[N:]
    numer = [1]
    for x in head: numer = polymul(numer, [sq(x), 1])
    Nm = [1]
    for _ in range(copies): Nm = polymul(Nm, numer)
    Nm = polymul(Nm, extra)
    return labels, head, tail, mu, polevalue, Nm, sq

def entries(machine, K, a0=0, N=0, copies=0, M=0, r=1, step=1, verbose=True):
    t0 = time.time()
    labels, head, tail, mu, polevalue, Nm, sq = setup(machine, K, a0, N, copies, M, r, step)
    h = K - N
    D = [1]
    for x in tail: D = polymul(D, [sq(x), 1])
    degD = len(D) - 1
    Dp = {a: math.prod(sq(b) - sq(a) for b in tail if b != a) for a in tail}
    Nm_at = {a: ev(Nm, -sq(a)) for a in tail}
    pv = {a: polevalue(a) for a in tail}
    q, rem = polydivmod_monic(Nm, D)
    maxe = len(q) + 2*h
    MU = [mu(e) for e in range(maxe + 1)]
    pw = {a: 1 for a in tail}
    A = []; B = []
    for s in range(2*h - 1):
        a_val = F(0); b_val = F(0)
        for e, c in enumerate(q):
            if c: a_val += MU[e]*c
        for a in tail:
            c = F(Nm_at[a]*pw[a], Dp[a]); a_val += c*pv[a][0]; b_val += c*pv[a][1]
        A.append(a_val); B.append(b_val)
        cs = rem[degD - 1] if degD > 0 else 0
        q = [cs] + q
        rem = [(rem[i - 1] if i > 0 else 0) - cs*D[i] for i in range(degD)]
        for a in tail: pw[a] *= -sq(a)
    if verbose: print('%s K=%d a0=%d N=%d copies=%d M=%d r=%d step=%d h=%d: entries %.0fs' % (machine, K, a0, N, copies, M, r, step, h, time.time() - t0), flush=True)
    return A, B, h, labels, head, tail

def detpoly(A, B, h, verbose=True):
    """det(A + XB) = det(B) * charpoly(-B^{-1}A)(X), coefficients low -> high as Fractions."""
    t0 = time.time()
    Af = flint.fmpq_mat([[flint.fmpq(A[i + j].numerator, A[i + j].denominator) for j in range(h)] for i in range(h)])
    Bf = flint.fmpq_mat([[flint.fmpq(B[i + j].numerator, B[i + j].denominator) for j in range(h)] for i in range(h)])
    dB = Bf.det()
    cp = (-Bf.solve(Af)).charpoly()
    dBF = F(int(dB.p), int(dB.q))
    poly = [F(int(cp[i].p), int(cp[i].q))*dBF for i in range(h + 1)]
    if verbose: print('  solve + charpoly %.0fs' % (time.time() - t0), flush=True)
    return poly

def content_exponents(machine, K, a0=0, N=0, copies=0, M=0, r=1, step=1):
    """exponent of every prime in the content of det(A + XB): returns (den dict, num dict) with the content = num/den."""
    buf = __import__('io').StringIO()
    with __import__('contextlib').redirect_stdout(buf):
        A, B, h, labels, head, tail = entries(machine, K, a0, N, copies, M, r, step)
        poly = detpoly(A, B, h)
    g = 0; l = 1
    for c in poly:
        g = gcd(g, c.numerator); l = l*c.denominator//gcd(l, c.denominator)
    den = l//gcd(l, g); num = g//gcd(l, g)
    top = 4*max(labels) + 200
    return factor_small(den, top), factor_small(num, top)

def detpoly_interp(A, B, h, verbose=True):
    t0 = time.time()
    Af = flint.fmpq_mat([[flint.fmpq(A[i + j].numerator, A[i + j].denominator) for j in range(h)] for i in range(h)])
    Bf = flint.fmpq_mat([[flint.fmpq(B[i + j].numerator, B[i + j].denominator) for j in range(h)] for i in range(h)])
    vals = []
    for X in range(h + 1):
        d = (Af + Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p), int(d.q)))
    if verbose: print('  %d determinants %.0fs' % (h + 1, time.time() - t0), flush=True)
    xs = list(range(h + 1)); coef = vals[:]
    for k in range(1, h + 1):
        for i in range(h, k - 1, -1): coef[i] = (coef[i] - coef[i - 1])/(xs[i] - xs[i - k])
    poly = [coef[h]]
    for k in range(h - 1, -1, -1):
        poly = polymul(poly, [F(-xs[k]), F(1)]); poly = [F(c) for c in poly]; poly[0] += coef[k]
    return poly

def factor_small(n, bound):
    f = {}
    for p in primerange(2, bound):
        if n % p == 0:
            e = 0
            while n % p == 0: n //= p; e += 1
            f[p] = e
    if n != 1: f.update({int(p): int(e) for p, e in factorint(n).items()})
    return f

def measure(machine, K, a0=0, N=0, copies=0, M=0, r=1, step=1, write=True):
    t0 = time.time()
    A, B, h, labels, head, tail = entries(machine, K, a0, N, copies, M, r, step)
    poly = detpoly(A, B, h)
    g = 0; l = 1
    for c in poly:
        g = gcd(g, c.numerator); l = l*c.denominator//gcd(l, c.denominator)
    den = l//gcd(l, g); num = g//gcd(l, g)
    top = 4*max(labels) + 200
    fd = factor_small(den, top); fn = factor_small(num, top)
    P = [c/F(num, den) for c in poly]
    assert all(c.denominator == 1 for c in P)
    H = max(abs(c.numerator) for c in P)
    logH = math.log(H)
    logden = sum(e*math.log(p) for p, e in fd.items()); lognum = sum(e*math.log(p) for p, e in fn.items())
    def lg(fr): return math.log(abs(fr.numerator)) - math.log(fr.denominator) if fr != 0 else float('-inf')
    logmax = max(lg(c) for c in poly)
    flint.ctx.prec = 2*H.bit_length() + 8*K*K + 4096
    if machine == 'zeta2int': xi = flint.arb.pi()**2/6
    elif machine == 'zeta2half': xi = flint.arb.pi()**2/2
    else:
        try: xi = flint.arb.const_catalan()
        except Exception:
            from mpmath import mp, catalan
            mp.dps = int(flint.ctx.prec*0.302) + 10; xi = flint.arb(str(+catalan))
    val = flint.arb(0)
    for c in reversed(P): val = val*xi + flint.arb(c.numerator)
    lv = float(abs(val).log())
    try: pos = bool(val > 0)
    except Exception: pos = None
    closeness = (lv - logH)/K**2; height = logH/K**2; net = lv/K**2
    two = fd.get(2, 0) - fn.get(2, 0)
    rng = lambda lo, hi: sum(e*math.log(p) for p, e in fd.items() if lo < p <= hi and p > 2)/K**2
    L = max(labels)
    split = {'two': two*math.log(2)/K**2, 'odd<=L/2': rng(2, L/2), '(L/2,2L/3]': rng(L/2, 2*L/3), '(2L/3,L]': rng(2*L/3, L), '(L,2L]': rng(L, 2*L), '>2L': rng(2*L, 10**9)}
    print('%s K=%d a0=%d N=%d copies=%d M=%d r=%d step=%d h=%d | closeness/K^2=%.4f height/K^2=%.4f net/K^2=%.4f | log P(G)=%.1f P>0:%s | %.0fs' % (
        machine, K, a0, N, copies, M, r, step, h, closeness, height, net, lv, pos, time.time() - t0))
    print('  height = intrinsic %.4f + den %.4f - num %.4f ; 2-adic exponent %d (%.4f per K^2)' % (logmax/K**2, logden/K**2, lognum/K**2, two, split['two']))
    print('  den split per K^2 (primes by position against the largest pole label L=%d): %s' % (L, ' '.join('%s %.4f' % (k, v) for k, v in split.items())))
    print('  e_p (den) :', ' '.join('%d:%d' % (p, e) for p, e in sorted(fd.items())))
    if fn: print('  e_p (num) :', ' '.join('%d:%d' % (p, e) for p, e in sorted(fn.items())))
    out = dict(machine=machine, K=K, a0=a0, N=N, copies=copies, M=M, r=r, step=step, h=h, closeness=closeness, height=height, net=net, logP=lv,
               positive=pos, log_intrinsic=logmax, log_den=logden, log_num=lognum, split=split,
               den={str(p): e for p, e in sorted(fd.items())}, num={str(p): e for p, e in sorted(fn.items())}, time=time.time() - t0)
    if write:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twisted_%s_K%d_a%d_N%d_c%d_M%d_r%d_s%d.json' % (machine, K, a0, N, copies, M, r, step))
        json.dump(out, open(path, 'w'), indent=1)
        print('  written', path)
    return out

def e2_detA(machine, K, **kw):
    """-v2(det A), A = the Hankel matrix at X = 0 (the constant coefficient of det(A + XB)).
    This is a lower bound for the 2-adic content exponent e2 and equals it for most K: for twist it differs (content = this + 1) only at
    K = 11, 61, 127 among K = 2..128 (the content's minimum sits at a higher power of X there).  Not for zeta2half, whose minimum is in the
    middle of the polynomial."""
    buf = __import__('io').StringIO()
    with __import__('contextlib').redirect_stdout(buf):
        A, B, h, labels, head, tail = entries(machine, K, **kw)
    Af = flint.fmpq_mat([[flint.fmpq(A[i + j].numerator, A[i + j].denominator) for j in range(h)] for i in range(h)])
    d = Af.det()
    if d == 0: raise ValueError('det A = 0 (e.g. twist K = 1)')
    v = lambda n: (abs(n) & -abs(n)).bit_length() - 1
    return v(int(d.q)) - v(int(d.p))

if __name__ == '__main__':
    if sys.argv[1] == 'e2':
        # usage: python twisted_hankel.py e2 MACHINE K1,K2,...
        for K in [int(x) for x in sys.argv[3].split(',')]:
            t0 = time.time(); e = e2_detA(sys.argv[2], K)
            print('%s K=%d e2 = -v2(det A) = %d (%.5f K^2) %.0fs' % (sys.argv[2], K, e, e/K**2, time.time() - t0), flush=True)
        sys.exit()
    machine = sys.argv[1]; K = int(sys.argv[2])
    ints = [int(x) for x in sys.argv[3:]]
    names = ['a0', 'N', 'copies', 'M', 'r', 'step']
    kw = dict(zip(names, ints))
    measure(machine, K, **kw)
