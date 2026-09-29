"""cubic.py — the cubic (χ₋₃) analogue of Zudilin's β-forms: linear forms in 1, L(2,χ₋₃), L(4,χ₋₃), … from Hurwitz
differences of Lai–Zhou/Sprang-type sums with roots on the thirds lattices.

R(t) = (2t + m2 n) · ∏_{i=1}^{m1 n} (t − i)(t + m2 n + i) · (t − m1 n + 1/3)_L · (t − m1 n + 2/3)_L / ∏_{j=1}^{q} (t + δ_j n)_{(m2 − 2δ_j) n + 1},
L = (2 m1 + m2) n, q odd, n even.  Poles at t = −k, k ∈ ∪_j [δ_j n, (m2 − δ_j) n], multiplicity = #covering blocks.
Even fold R(−t − m2 n) = R(t)  ⇒  a_{i, m2 n − k} = (−1)^i a_{i,k}  ⇒  ρ_i = Σ_k a_{i,k} = 0 for odd i.
S_θ = Σ_{t ≥ −m} R(t + θ), m = m2 n/2 (the sum may start anywhere in [−(m1+m2) n, 1] since the 2/3-block kills R(t + 1/3)
there and vice versa);  S_{1/3} − S_{2/3} = a0 + Σ_{i even} ρ_i · 3^i L(i, χ₋₃),  L(i,χ₋₃) = 3^{−i}(ζ(i,1/3) − ζ(i,2/3)).
Partial sums: (ℓ + 1/3)^{−i} − (ℓ + 2/3)^{−i} for ℓ ∈ [k − m, −1] (k < m) or −[0, k − m − 1] (k > m): denominators 3ℓ ± 1
up to 3(m − k): the base-3 lattice extent (1.5× the β case).
usage: python cubic.py m1 m2 d1,d2,...,dq n [nocheck]
"""
import sys, math, time
sys.set_int_max_str_digits(0)
from fractions import Fraction
from flint import fmpq
from mpmath import mp, mpf, zeta as mzeta, log as mlog
from sympy import primerange, factorint

def cubic_forms(m1, m2, deltas, n, want_coeffs=False):
    q = len(deltas); assert q % 2 == 1 and n % 2 == 0
    assert (q - 2)*m2 - 6*m1 - 2*sum(deltas) >= 0, 'degree condition (q-2) m2 - 6 m1 - 2 sum(delta) >= 0 violated: the sums diverge'
    assert all(0 <= d < m2/2 for d in deltas)
    L = (2*m1 + m2)*n; m = m2*n//2
    blocks = [(d*n, (m2 - 2*d)*n + 1) for d in deltas]          # (start h, length)
    K = sorted(set(k for h, Lj in blocks for k in range(h, h + Lj)))
    third = fmpq(1, 3); twothird = fmpq(2, 3)
    A = {}; mult = {}
    for k in K:
        # numerator factor values at t = -k
        alphas = []
        cen = m2*n - 2*k                 # (2t + m2 n) = 2 (t + m2 n/2): monic factor with value cen/2, constant 2
        z = 0
        if cen == 0: z += 1
        else: alphas.append(fmpq(cen, 2))
        for i in range(1, m1*n + 1):
            alphas.append(fmpq(-k - i)); alphas.append(fmpq(m2*n + i - k))
        for th in (third, twothird):
            for i in range(L): alphas.append(th + (i - k - m1*n))
        betas = []; order = 0
        for h, Lj in blocks:
            for i in range(Lj):
                d = h + i - k
                if d == 0: order += 1
                else: betas.append(fmpq(d))
        sk = order - z
        if sk <= 0: continue
        C = fmpq(2)           # the constant 2 of the factor 2(t + m2 n/2) (= 2u at the centre)
        for a in alphas: C *= a
        for b in betas: C /= b
        # log-series
        Lc = [None]
        for r in range(1, sk):
            S = fmpq(0)
            for a in alphas: S += (1/a)**r
            for b in betas: S -= (1/b)**r
            Lc.append((-1)**(r + 1)*S/r)
        bser = [C]
        for r in range(1, sk):
            acc = fmpq(0)
            for j in range(1, r + 1): acc += j*Lc[j]*bser[r - j]
            bser.append(acc/r)
        for i in range(1, sk + 1): A[(i, k)] = bser[sk - i]
        mult[k] = sk
    smax = max(mult.values())
    rho = {i: fmpq(0) for i in range(1, smax + 1)}
    for (i, k), c in A.items(): rho[i] += c
    # a0: partial sums of (l+1/3)^-i - (l+2/3)^-i
    def ps_diff(i, lo, hi):   # sum_{l=lo}^{hi}
        s = fmpq(0)
        for l in range(lo, hi + 1): s += (1/(l + third))**i - (1/(l + twothird))**i
        return s
    a0 = fmpq(0)
    for (i, k), c in A.items():
        if k < m: a0 += c*ps_diff(i, k - m, -1)
        elif k > m: a0 -= c*ps_diff(i, 0, k - m - 1)
    info = dict(L=L, m=m, K=K, mult=mult, blocks=blocks, smax=smax, q=q)
    if want_coeffs: info['A'] = A
    return rho, a0, info

def R_exact(m1, m2, deltas, n, info, t):
    L = info['L']; blocks = info['blocks']
    v = fmpq(2)*t + m2*n
    for i in range(1, m1*n + 1): v *= (t - i)*(t + m2*n + i)
    for th in (fmpq(1, 3), fmpq(2, 3)):
        for i in range(L): v *= (t - m1*n + th + i)
    for h, Lj in blocks:
        for i in range(Lj): v /= (t + h + i)
    return v
def identity_check(m1, m2, deltas, n, info, A):
    import random
    random.seed(3)
    for _ in range(3):
        t = fmpq(random.randint(1, 10**5), random.randint(1, 10**5)) + random.randint(0, 5)
        lhs = R_exact(m1, m2, deltas, n, info, t)
        rhs = fmpq(0)
        for (i, k), c in A.items(): rhs += c/(t + k)**i
        if lhs != rhs: return False
    return True
def numeric_diff(m1, m2, deltas, n, info, dps):
    """S_{1/3} - S_{2/3}: mpf products for the terms, mpmath nsum (Richardson/Shanks) for the power-law tail."""
    mp.dps = dps
    L = info['L']; blocks = info['blocks']
    def R(t):
        v = mpf(2*t + m2*n)
        for i in range(1, m1*n + 1): v *= (t - i)*(t + m2*n + i)
        for th in (mpf(1)/3, mpf(2)/3):
            for i in range(L): v *= (t - m1*n + th + i)
        for h, Lj in blocks:
            for i in range(Lj): v /= (t + h + i)
        return v
    T0 = -(m1 + m2)*n
    f = lambda t: R(t + mpf(1)/3) - R(t + mpf(2)/3)
    head = mpf(0); mx = mpf(0); T1 = 4*(m1 + m2)*n + 20
    for T in range(T0, T1):
        term = f(mpf(T)); head += term; mx = max(mx, abs(term))
    tail = mp.nsum(lambda t: f(t), [T1, mp.inf])
    return head + tail, mx

def run(m1, m2, deltas, n, check=True):
    t0 = time.time()
    rho, a0, info = cubic_forms(m1, m2, deltas, n, want_coeffs=True)
    q = info['q']; smax = info['smax']
    odd = [i for i in rho if i % 2 == 1 and rho[i] != 0]
    print('cubic forms: m1=%d m2=%d deltas=%s n=%d : L=%d, poles k in [%d,%d] (%d), max multiplicity %d  [%.1fs]' % (m1, m2, deltas, n, info['L'], min(info['K']), max(info['K']), len(info['K']), smax, time.time() - t0))
    print('  partial fractions reproduce R(t) exactly: %s ; odd rho_i vanish (even fold): %s ; nonzero even rho: %s' % (identity_check(m1, m2, deltas, n, info, info['A']), 'yes' if not odd else 'NO %s' % odd, [i for i in rho if i % 2 == 0 and rho[i] != 0]))
    sys.stdout.flush()
    if check:
        dps = 80 + 3*n*(2*m1 + m2)
        diff, mx = numeric_diff(m1, m2, deltas, n, info, dps)
        form = mpf(int(a0.p))/int(a0.q)
        for i in range(2, smax + 1, 2):
            Li = (mzeta(i, mpf(1)/3) - mzeta(i, mpf(2)/3))/mpf(3)**i
            form += mpf(int(rho[i].p))/int(rho[i].q)*mpf(3)**i*Li
        print('  numeric S_{1/3}-S_{2/3} = %s ; form = %s ; rel diff %s ; largest term %s ; log|S|/n = %.4f' % (mp.nstr(diff, 15), mp.nstr(form, 15), mp.nstr(abs(diff - form)/abs(diff), 3) if diff != 0 else 'n/a', mp.nstr(mx, 4), float(mlog(abs(diff)))/n if diff != 0 else float('nan')))
    # primitive integer form
    coeffs = [a0] + [rho[i] for i in range(2, smax + 1, 2)]
    D = 1
    for c in coeffs: D = D*int(c.q)//math.gcd(D, int(c.q))
    ints = [int(c.p)*(D//int(c.q)) for c in coeffs]
    g = 0
    for c in ints: g = math.gcd(g, c)
    def facsmall(N):
        f = {}
        for p in primerange(2, 6*(m2*n) + 100):
            e = 0
            while N % p == 0: N //= p; e += 1
            if e: f[p] = e
        if N != 1: f.update({int(p): int(e) for p, e in factorint(N).items()})
        return f
    fac = facsmall(D); facg = facsmall(g)
    logD = sum(e*math.log(p) for p, e in fac.items()); logg = sum(e*math.log(p) for p, e in facg.items())
    print('  height: log D/n = %.4f (D = lcm of denominators), content g: %s (log g/n = %.4f), primitive height log(D/g)/n = %.4f ; den(a0)/n = %.4f' % (logD/n, facg if g > 1 else 1, logg/n, (logD - logg)/n, math.log(int(a0.q))/n))
    Mp = (m2 - 2*min(deltas))*n
    byrange = {'p<=Mp/2': 0.0, 'Mp/2<p<=Mp': 0.0, 'Mp<p<=1.5Mp': 0.0, 'p>1.5Mp': 0.0}
    for p, e in fac.items():
        key = 'p<=Mp/2' if p <= Mp/2 else ('Mp/2<p<=Mp' if p <= Mp else ('Mp<p<=1.5Mp' if p <= 1.5*Mp else 'p>1.5Mp'))
        byrange[key] += e*math.log(p)
    print('  height by prime range (M\' = (m2-2 delta_min) n = %d): ' % Mp + '  '.join('%s: %.4f' % (k, v/n) for k, v in byrange.items()))
    print('  largest primes: ' + ' '.join('%d^%d' % (p, e) for p, e in sorted(fac.items())[-12:]))
    return rho, a0, info

def decay(m1, m2, deltas):
    """-lim log|S|/n for the unnormalised R (Stirling saddle): max over kappa>0 of F(kappa),
    F = sum over Gamma ratios of (a ln a - a) in units of n, with R(t) evaluated at t = (m1+kappa) n."""
    import mpmath
    mpmath.mp.dps = 30
    ph = lambda a: a*mpmath.log(a) - a if a > 0 else mpmath.mpf(0)
    def F(kap):
        v = (ph(m1 + kap) - ph(kap)) + (ph(2*m1 + m2 + kap) - ph(m1 + m2 + kap)) + 2*(ph(2*m1 + m2 + kap) - ph(kap))
        for d in deltas: v -= ph(m1 + m2 - d + kap) - ph(m1 + d + kap)
        return v
    # maximise on (0, inf): F'(kap) = sum log-ratios; unimodal in practice — golden section on log scale then refine
    best = None
    for k in [mpmath.mpf(10)**(e/20) for e in range(-120, 80)]:
        v = F(k)
        if best is None or v > best[0]: best = (v, k)
    k = best[1]
    # refine by ternary search around k
    lo, hi = k/3, k*3
    for _ in range(200):
        a = lo + (hi - lo)/3; b = hi - (hi - lo)/3
        if F(a) < F(b): lo = a
        else: hi = b
    k = (lo + hi)/2
    return float(-F(k)), float(k)

def table(m1, m2, deltas, ns):
    d, k0 = decay(m1, m2, deltas)
    L = (2*m1 + m2)
    print('m1=%d m2=%d deltas=%s : Stirling decay of the unnormalised R: %.4f per n (saddle kappa = %.4f); 3-adic normalisation 3^{2L} contributes %.4f per n' % (m1, m2, deltas, d, k0, 2*L*math.log(3)))
    for n in ns:
        rho, a0, info = cubic_forms(m1, m2, deltas, n, want_coeffs=False)
        smax = info['smax']
        coeffs = [a0] + [rho[i] for i in range(2, smax + 1, 2)]
        D = 1
        for c in coeffs: D = D*int(c.q)//math.gcd(D, int(c.q))
        ints = [int(c.p)*(D//int(c.q)) for c in coeffs]
        g = 0
        for c in ints: g = math.gcd(g, c)
        lam = Fraction(D, g)
        # per-prime exponents of lambda
        num, den = lam.numerator, lam.denominator
        def facsmall(N):
            f = {}
            for p in primerange(2, 6*(m2*n) + 100):
                e = 0
                while N % p == 0: N //= p; e += 1
                if e: f[p] = e
            if N != 1: f.update({int(p): int(e) for p, e in factorint(N).items()})
            return f
        fn, fd = facsmall(num), facsmall(den)
        ex = {p: fn.get(p, 0) - fd.get(p, 0) for p in set(fn) | set(fd)}
        loglam = sum(e*math.log(p) for p, e in ex.items())
        v3 = ex.get(3, 0); v2 = ex.get(2, 0)
        dps = 60 + 3*n*(2*m1 + m2)
        S, mx = numeric_diff(m1, m2, deltas, n, info, dps)
        net = (loglam + float(mlog(abs(S))))/n if S != 0 else float('nan')
        Mp = (m2 - 2*min(deltas))*n
        big = {p: e for p, e in ex.items() if p > 3}
        byr = [sum(e*math.log(p) for p, e in big.items() if lo < p <= hi) for lo, hi in ((0, Mp/2), (Mp/2, Mp), (Mp, 1.5*Mp), (1.5*Mp, 10**9))]
        print('  n=%2d: log(lambda)/n = %.4f [v2=%d, v3=%d: %.4f per n; primes>3: %.4f = (%.3f | %.3f | %.3f | %.3f) by p<=M\'/2, <=M\', <=1.5M\', >1.5M\', M\'=%d] ; log|S|/n = %.4f (naive Stirling %.4f) ; NET log|lambda S|/n = %+.4f' % (
            n, loglam/n, v2, v3, (v2*math.log(2) + v3*math.log(3))/n, sum(byr)/n, byr[0]/n, byr[1]/n, byr[2]/n, byr[3]/n, Mp, float(mlog(abs(S)))/n, -d, net))
        sys.stdout.flush()

if __name__ == '__main__':
    if sys.argv[1] == 'table':
        m1 = int(sys.argv[2]); m2 = int(sys.argv[3]); deltas = [int(x) for x in sys.argv[4].split(',')]; ns = [int(x) for x in sys.argv[5].split(',')]
        table(m1, m2, deltas, ns); sys.exit()
    m1 = int(sys.argv[1]); m2 = int(sys.argv[2]); deltas = [int(x) for x in sys.argv[3].split(',')]; n = int(sys.argv[4])
    run(m1, m2, deltas, n, check=not (len(sys.argv) > 5 and sys.argv[5] == 'nocheck'))
