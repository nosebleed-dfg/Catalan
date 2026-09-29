"""thirds.py — Zudilin-type ("type A") construction with poles on the thirds lattice (the cubic analogue of multibeta.py).

R(t) = (2t+h0) (t+1)_{h0-1} / prod_j [(t + eta_j n + 1/3)_{K_j} (t + eta_j n + 2/3)_{K_j}],   h0 = eta0 n + 1,  K_j = (eta0 - 2 eta_j) n + 1,  n even.
Palindrome: R(-t-h0) = -R(t).  Untwisted sum  r = sum_{nu>=0} R(nu)  (the full-line sum vanishes, so r is a genuine half-line object):
   r = a0 + sum_{i even} 3^i A_i L(i,chi_-3) + sum_{i odd>=3} (3^i - 1) A_i zeta(i),     A_i = sum over the 1/3-coset poles of a_{i,k}.
(The alternating sum is a pi-identity; see the results log, "The cubic construction, type A".)
The derivative forms sum_nu R^{(2d)}(nu) shift i -> i+2d and kill the lowest zeta: d=1 gives a form in 1, L(4), zeta(5), L(6), zeta(7), ...

usage:  python thirds.py forms ETA n            exact partial fractions, a0, primitive height, per-prime ledger by range, numeric r, NET
        python thirds.py peak ETA [n...]        Stirling closeness (peak of log R(xn)/n) and boundary value; finite-n check
        python thirds.py ledger ETA n           per-prime table of v_p(primitive coefficients) for the cubic tax analysis
"""
import sys, time, math
from fractions import Fraction
from flint import fmpq, fmpz
from mpmath import mp, mpf, zeta as mzeta, nsum, inf, log, fabs, loggamma, exp
import sympy

DOUBLE = False   # double-root variant: numerator (t+1)_{h0-1} (t+1/2)_{h0}, no (2t+h0) factor (sums over Z and Z+1/2 share the coefficients)
COSET = 3        # pole cosets 1/COSET and (COSET-1)/COSET: 3 -> L(i,chi_-3) (with zeta(odd)); 4 -> beta(i) (with zeta(odd))
CUBE = 1         # numerator multiplicity: (t+1)_{h0-1}^CUBE (CUBE=3 allows the second-derivative form sum R''(nu), which kills zeta(3) and drops L(2))

def q2mp(v): return mpf(int(v.p))/mpf(int(v.q))

def poles(eta, n):
    eta0, es = eta[0], eta[1:]
    P = {}
    for ej in es:
        K = (eta0 - 2*ej)*n + 1
        for th in (Fraction(1, COSET), Fraction(COSET - 1, COSET)):
            for l in range(K):
                pos = -(ej*n + th + l)
                P[pos] = P.get(pos, 0) + 1
    return P

def partial_fractions(eta, n):
    """a[(pos, i)] with R(t) = sum a/(t - pos)^i (exact fmpq), via log-series expansion at each pole."""
    eta0, es = eta[0], eta[1:]; h0 = eta0*n + 1
    P = poles(eta, n)
    num = [Fraction(m) for m in range(1, h0)]*CUBE
    if DOUBLE: num += [Fraction(2*m + 1, 2) for m in range(0, h0)]
    den = []
    for ej in es:
        K = (eta0 - 2*ej)*n + 1
        for th in (Fraction(1, COSET), Fraction(COSET - 1, COSET)):
            for l in range(K): den.append(ej*n + th + l)
    a = {}
    for pos, s in P.items():
        C = fmpq(1) if DOUBLE else fmpq(2)
        series = [fmpq(0)]*s; series[0] = fmpq(1)
        def mul_series(ser, alpha, sign):
            out = [fmpq(0)]*s
            inva = 1/alpha
            if sign > 0:
                for j in range(s):
                    out[j] += ser[j]
                    if j + 1 < s: out[j + 1] += ser[j]*inva
            else:
                pw = [fmpq(1)]
                for r in range(1, s): pw.append(pw[-1]*(-inva))
                for j in range(s):
                    for r in range(s - j): out[j + r] += ser[j]*pw[r]
            return out
        p0 = fmpq(pos.numerator, pos.denominator)
        if not DOUBLE:
            alph = p0 + fmpq(h0, 2); C *= alph; series = mul_series(series, alph, +1)
        for m in num:
            alph = p0 + fmpq(m.numerator, m.denominator); C *= alph; series = mul_series(series, alph, +1)
        for d in den:
            alph = p0 + fmpq(d.numerator, d.denominator)
            if alph == 0: continue
            C /= alph; series = mul_series(series, alph, -1)
        for i in range(1, s + 1): a[(pos, i)] = C*series[s - i]
    return a, P

def R_exact(eta, n, t):
    eta0, es = eta[0], eta[1:]; h0 = eta0*n + 1
    tt = fmpq(t.numerator, t.denominator)
    v = fmpq(1) if DOUBLE else 2*tt + h0
    for m in range(1, h0): v *= (tt + m)**CUBE
    if DOUBLE:
        for m in range(0, h0): v *= (tt + fmpq(2*m + 1, 2))
    for ej in es:
        K = (eta0 - 2*ej)*n + 1
        for th in (fmpq(1, COSET), fmpq(COSET - 1, COSET)):
            for l in range(K): v /= (tt + ej*n + th + l)
    return v

def check_identity(eta, n, a, trials=2):
    import random
    for _ in range(trials):
        t = Fraction(random.randint(1, 10**6), random.randint(1, 10**3)) + Fraction(1, 7)
        lhs = R_exact(eta, n, t); rhs = fmpq(0); tt = fmpq(t.numerator, t.denominator)
        for (pos, i), c in a.items(): rhs += c/(tt - fmpq(pos.numerator, pos.denominator))**i
        if lhs != rhs: return False
    return True

def logR_num(eta, n, t):
    """log R(t) for real t > 0 via loggamma (mpf)."""
    eta0, es = eta[0], eta[1:]; h0 = eta0*n + 1
    t = mpf(t)
    v = (loggamma(t + h0 + mpf(1)/2) - loggamma(t + mpf(1)/2)) if DOUBLE else log(2*t + h0)
    v += CUBE*(loggamma(t + h0) - loggamma(t + 1))
    for ej in es:
        K = (eta0 - 2*ej)*n + 1
        for th in (mpf(1)/COSET, mpf(COSET - 1)/COSET):
            v -= loggamma(t + ej*n + th + K) - loggamma(t + ej*n + th)
    return v

def numeric_r(eta, n, deriv=0, shift=0):
    """sum_{nu>=0} R^{(deriv)}(nu + shift) numerically (deriv even); direct summation plus nsum tail."""
    mp.dps = 50
    sh = mpf(shift.numerator)/shift.denominator if isinstance(shift, Fraction) else mpf(shift)
    if deriv == 0:
        f = lambda nu: exp(logR_num(eta, n, nu + sh))
    else:
        from mpmath import diff
        f = lambda nu: diff(lambda t: exp(logR_num(eta, n, t)), nu + sh, deriv)
    head = sum(f(nu) for nu in range(0, 3000))
    tail = nsum(lambda k: f(k), [3000, inf])
    return head + tail

def cmd_combo(eta, n):
    """Two forms from the same R: r0 = sum R(nu), rh = sum R(nu+1/2); eliminate zeta(3) (q=3) or report both; primitive heights and NETs."""
    t0 = time.time()
    a, P = partial_fractions(eta, n)
    h0 = eta[0]*n + 1
    a0_0, out0 = linear_form(eta, n, a, 0, Fraction(0))
    a0_h, outh = linear_form(eta, n, a, 0, Fraction(1, 2))
    mp.dps = 50
    r0 = numeric_r(eta, n, 0, Fraction(0)); rh = numeric_r(eta, n, 0, Fraction(1, 2))
    def evalform(a0, out):
        v = q2mp(a0)
        for (kind, j), c in out.items():
            if c == 0: continue
            v += q2mp(c)*(L_chi3(j) if kind == 'L' else mzeta(j))
        return v
    print('combo eta=%s n=%d (h0=%d): r0 = %s (form check |diff| %s) ; r_half = %s (|diff| %s) ; ratio r_half/r0 = %s  [%.1fs]' % (
        eta, n, h0, mp.nstr(r0, 12), mp.nstr(fabs(r0 - evalform(a0_0, out0)), 3), mp.nstr(rh, 12), mp.nstr(fabs(rh - evalform(a0_h, outh)), 3), mp.nstr(rh/r0, 10), time.time() - t0))
    keys = sorted(set(out0) | set(outh), key=lambda k: k[1])
    for k in keys:
        c0 = out0.get(k, fmpq(0)); ch = outh.get(k, fmpq(0))
        if c0 == 0 and ch == 0: continue
        print('   %s(%d): weight ratio r_half/r0 coefficient = %s' % (k[0], k[1], str(ch/c0) if c0 != 0 else 'inf'))
    # eliminate zeta(3) when present
    kz = ('zeta', 3)
    if kz in out0 and out0[kz] != 0:
        w0 = outh[kz]; wh = out0[kz]     # w0*r0 - wh*rh kills zeta(3)
        g = fmpq(1)  # scale to make weights coprime integers
        from math import gcd
        num0, numh = int((w0).p*(wh).q), int((wh).p*(w0).q)
        gg = gcd(num0, numh); num0 //= gg; numh //= gg
        comb_a0 = num0*a0_0 - numh*a0_h
        comb = {k: num0*out0.get(k, fmpq(0)) - numh*outh.get(k, fmpq(0)) for k in keys}
        comb = {k: v for k, v in comb.items() if v != 0}
        rc = num0*r0 - numh*rh
        vc = q2mp(comb_a0) + sum(q2mp(c)*(L_chi3(j) if kind == 'L' else mzeta(j)) for (kind, j), c in comb.items())
        print('   combination %d*r0 - %d*r_half: remaining %s ; value %s (form check |diff| %s) ; log|value|/n = %.6f' % (
            num0, numh, ', '.join('%s(%d)' % k for k in sorted(comb, key=lambda k: k[1])), mp.nstr(rc, 12), mp.nstr(fabs(rc - vc), 3), log(fabs(rc))/n))
        for name, (aa, oo, rr) in (('r0', (a0_0, out0, r0)), ('r_half', (a0_h, outh, rh)), ('combination', (comb_a0, comb, rc))):
            coeffs = [aa] + [v for v in oo.values() if v != 0]
            D, gc = primitive(coeffs)
            hD = math.log(int(D)); hg = math.log(int(gc)) if gc > 0 else 0.0
            M = (eta[0] - 2*min(eta[1:]))*n + 1
            ranges = [('p<=M', 1, M), ('M<p<=2M', M, 2*M), ('2M<p<=3M', 2*M, 3*M), ('3M<p<=6M', 3*M, 6*M), ('p>6M', 6*M, 10**12)]
            sD, fD = factor_split(D, ranges)
            print('   %-12s primitive height log(D/g)/n = %.4f ; NET = %+.4f ; D by range: %s' % (name, (hD - hg)/n, (hD - hg)/n + float(log(fabs(rr)))/n, '  '.join('%s: %.3f' % (lab, sD[lab]/n) for lab, _, _ in ranges)))
            if name == 'combination': print('   combination D factorization:', ' '.join('%d^%d' % (p, e) for p, e in sorted(fD.items())))

def L_chi3(i): return (mzeta(i, mpf(1)/COSET) - mzeta(i, mpf(COSET - 1)/COSET))/mpf(COSET)**i   # L(i,chi_-3) for COSET=3, beta(i) for COSET=4

def hurwitz_split(j, c):
    """zeta(j, c) for c in {1/3, 2/3, 1/6, 5/6} as (coefficient of zeta(j), coefficient of L(j,chi_-3)) — exact fmpq.
    zeta(j,1/3) = ((3^j-1) z + 3^j L)/2,  zeta(j,2/3) = ((3^j-1) z - 3^j L)/2,
    zeta(j,1/6) = 6^j((1-2^-j)(1-3^-j) z + (1+2^-j) L)/2,  zeta(j,5/6) = 6^j((1-2^-j)(1-3^-j) z - (1+2^-j) L)/2."""
    t3 = fmpz(3)**j; t2 = fmpz(2)**j; t6 = fmpz(6)**j
    t4 = fmpz(4)**j
    if c == Fraction(1, 4): return (fmpq(t2*(t2 - 1), 2), fmpq(t4, 2))
    if c == Fraction(3, 4): return (fmpq(t2*(t2 - 1), 2), -fmpq(t4, 2))
    if c == Fraction(1, 3): return (fmpq(t3 - 1, 2), fmpq(t3, 2))
    if c == Fraction(2, 3): return (fmpq(t3 - 1, 2), -fmpq(t3, 2))
    if c == Fraction(1, 6): return (fmpq(t6, 2)*(1 - fmpq(1, t2))*(1 - fmpq(1, t3)), fmpq(t6, 2)*(1 + fmpq(1, t2)))
    if c == Fraction(5, 6): return (fmpq(t6, 2)*(1 - fmpq(1, t2))*(1 - fmpq(1, t3)), -fmpq(t6, 2)*(1 + fmpq(1, t2)))
    raise ValueError(c)

def linear_form(eta, n, a, deriv=0, shift=Fraction(0)):
    """Exact coefficients of r = sum_{nu>=0} R^{(deriv)}(nu + shift) as a0 + sum_j c_j X_j, X_j = L(j,chi_-3) or zeta(j).
    shift in {0, 1/2}. R^{(2d)} has coefficient (i)_{2d} a_{i,k} at order i+2d."""
    coef = {}
    a0 = fmpq(0)
    for (pos, i), c in a.items():
        j = i + deriv
        if deriv:
            m = fmpq(1)
            for r in range(i, i + deriv): m *= r
            c = c*m
        theta = Fraction(1, COSET) if (pos % 1) == Fraction(COSET - 1, COSET) else Fraction(COSET - 1, COSET)
        k = int(-pos - theta); assert Fraction(k) + theta == -pos
        # sum_{nu>=0} (nu + shift + k + theta)^{-j} = zeta(j, cc) - sum_{l<k'} (l + cc)^{-j},  cc = frac(shift+theta) in (0,1), k' = k + floor(shift+theta)
        c0 = shift + theta; kk = k + (c0.numerator // c0.denominator); cc = c0 - (c0.numerator // c0.denominator)
        coef.setdefault(j, {})
        coef[j][cc] = coef[j].get(cc, fmpq(0)) + c
        ccq = fmpq(cc.numerator, cc.denominator)
        ps = fmpq(0)
        for l in range(kk): ps += 1/(l + ccq)**j
        a0 -= c*ps
    out = {}
    for j, d in coef.items():
        cz = fmpq(0); cL = fmpq(0)
        for cc, A in d.items():
            z, L = hurwitz_split(j, cc)
            cz += A*z; cL += A*L
        if j == 1: assert cz == 0, 'divergent'
        else: out[('zeta', j)] = cz
        out[('L', j)] = cL
    return a0, out

def primitive(coeffs):
    """coeffs: list of fmpq -> (D = lcm of denominators, g = gcd of numerators of D*coeffs)"""
    D = fmpz(1)
    for c in coeffs: D = D*c.q//D.gcd(c.q) if hasattr(D, 'gcd') else fmpz(int(sympy.ilcm(int(D), int(c.q))))
    g = fmpz(0)
    for c in coeffs:
        v = (c*D).p
        g = fmpz(int(sympy.igcd(int(g), int(v))))
    return D, g

def factor_split(N, ranges):
    """log-contributions of the factorization of N by prime ranges; ranges = list of (label, lo, hi]"""
    f = sympy.factorint(int(N))
    out = {lab: 0.0 for lab, lo, hi in ranges}
    for p, e in f.items():
        for lab, lo, hi in ranges:
            if lo < p <= hi: out[lab] += e*math.log(p); break
    return out, f

def cmd_forms(eta, n, deriv=0):
    t0 = time.time()
    a, P = partial_fractions(eta, n)
    q = len(eta) - 1; h0 = eta[0]*n + 1; Kmin = (eta[0] - 2*max(eta[1:]))*n + 1
    print('thirds type A: eta=%s n=%d deriv=%d : h0=%d, %d poles, max multiplicity %d ; identity %s  [%.1fs]' % (
        eta, n, deriv, h0, len(P), max(P.values()), check_identity(eta, n, a), time.time() - t0))
    a0, out = linear_form(eta, n, a, deriv)
    keys = sorted(out.keys(), key=lambda k: k[1])
    nz = [(k, v) for k, v in out.items() if v != 0]
    print('  nonzero coefficients:', ', '.join('%s(%d)' % k for k, v in sorted(nz, key=lambda kv: kv[0][1])))
    # numeric check
    mp.dps = 50
    r = numeric_r(eta, n, deriv)
    val = q2mp(a0)
    for (kind, j), c in out.items():
        if c == 0: continue
        val += q2mp(c)*(L_chi3(j) if kind == 'L' else mzeta(j))
    print('  numeric r = %s ; form = %s ; |diff| = %s ; log|r|/n = %.6f' % (mp.nstr(r, 15), mp.nstr(val, 15), mp.nstr(fabs(r - val), 3), log(fabs(r))/n))
    coeffs = [a0] + [v for k, v in nz]
    D, g = primitive(coeffs)
    M = (eta[0] - 2*min(eta[1:]))*n + 1     # longest block length K_max
    ranges = [('p<=M', 1, M), ('M<p<=2M', M, 2*M), ('2M<p<=3M', 2*M, 3*M), ('p>3M', 3*M, 10**12)]
    sD, fD = factor_split(D, ranges); sg, fg = factor_split(g, ranges)
    hD = math.log(int(D)); hg = math.log(int(g)) if g > 0 else 0.0
    # normalisation N_n = prod_j (K_j!)^2 3^{h0-1-2 sum K_j}/(h0-1)!  (top coefficients ~ integral up to the cross-coset defects)
    lN = CUBE*(-float(loggamma(h0)) + (h0 - 1)*math.log(3))
    if DOUBLE: lN += -float(loggamma(h0 + 1)) + h0*math.log(3)
    for ej in eta[1:]:
        K = (eta[0] - 2*ej)*n + 1
        lN += 2*float(loggamma(K + 1)) - 2*K*math.log(3)
    print('  height: log D/n = %.4f, content log g/n = %.4f, primitive log(D/g)/n = %.4f ; normalised: height %.4f, closeness %.4f ; NET log|r D/g|/n = %+.4f   (M = K_max = %d)' % (
        hD/n, hg/n, (hD - hg)/n, (hD - hg - lN)/n, (float(log(fabs(r))) + lN)/n, (hD - hg)/n + float(log(fabs(r)))/n, M))
    print('  D by prime range (per n): ' + '  '.join('%s: %.4f' % (lab, sD[lab]/n) for lab, _, _ in ranges) + ' ; content by range: ' + '  '.join('%s: %.4f' % (lab, sg[lab]/n) for lab, _, _ in ranges))
    print('  D factorization:', ' '.join('%d^%d' % (p, e) for p, e in sorted(fD.items())))
    return a, out, a0, r, D, g

def cmd_peak(eta, ns):
    """closeness: lim (1/n) log max_nu R(nu) with the normalisation N_n = prod_j (K_j!)^2 * 3^{h0-1-2 sum K_j} / (h0-1)! (makes top coefficients ~integral);
    reports the normalised peak, the boundary value, and the finite-n numeric log|r|/n with the same normalisation."""
    eta0, es = eta[0], eta[1:]
    mp.dps = 30
    def logN(n):
        h0 = eta0*n + 1
        v = CUBE*(-loggamma(h0) + (h0 - 1)*log(3))
        if DOUBLE: v += -loggamma(h0 + 1) + h0*log(3)
        for ej in es:
            K = (eta0 - 2*ej)*n + 1
            v += 2*loggamma(K + 1) - 2*K*log(3)
        return v
    # asymptotic (Stirling) normalised log R(xn)/n
    def phi(x):
        x = mpf(x)
        v = CUBE*((x + eta0)*log(x + eta0) - (x*log(x) if x > 0 else 0) - eta0)   # (t+1)_{h0-1}^CUBE ~ (Gamma(t+h0)/Gamma(t+1))^CUBE, minus the n ln n part handled by normalisation
        v += (CUBE - 1)*(eta0*log(3) - (eta0*log(eta0) - eta0))
        if DOUBLE: v += (x + eta0)*log(x + eta0) - (x*log(x) if x > 0 else 0) - eta0 + eta0*log(3) - (eta0*log(eta0) - eta0)
        for ej in es:
            K = eta0 - 2*ej
            for _ in range(2):
                v -= (x + ej + K)*log(x + ej + K) - (x + ej)*log(x + ej) - K
                v += K*log(K) - K       # K_j!  (Stirling constant part)
                v -= K*log(3)
        v += eta0*log(3) - (eta0*log(eta0) - eta0)   # 3^{h0-1}/(h0-1)!
        return v
    xs = [mpf(k)/2000 for k in range(0, 40000)]
    vals = [(phi(x), x) for x in xs]
    best = max(vals)
    # refine
    lo, hi = max(best[1] - mpf(1)/2000, 0), best[1] + mpf(1)/2000
    for _ in range(100):
        m1 = lo + (hi - lo)/3; m2 = hi - (hi - lo)/3
        if phi(m1) < phi(m2): lo = m1
        else: hi = m2
    x0 = (lo + hi)/2
    print('peak closeness (normalised) eta=%s : max_x phi(x) = %.6f at x0 = %.5f ; boundary phi(0) = %.6f' % (eta, phi(x0), x0, phi(0)))
    for n in ns:
        r = numeric_r(eta, n)
        print('   n=%4d: log|r|/n (normalised) = %.6f ; log N_n/n = %.4f ; raw log|r|/n = %.6f' % (n, (log(fabs(r)) + logN(n))/n, logN(n)/n, log(fabs(r))/n))
    return float(phi(x0)), float(x0)

if __name__ == '__main__':
    if 'double' in sys.argv:
        DOUBLE = True; sys.argv.remove('double')
    if 'cube' in sys.argv:
        CUBE = 3; sys.argv.remove('cube')
    if 'quarter' in sys.argv:
        COSET = 4; sys.argv.remove('quarter')
    cmd = sys.argv[1]
    eta = [int(x) for x in sys.argv[2].split(',')]
    if cmd == 'forms':
        n = int(sys.argv[3]); deriv = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        cmd_forms(eta, n, deriv)
    elif cmd == 'peak':
        ns = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else []
        cmd_peak(eta, ns)
    elif cmd == 'combo':
        for n in [int(x) for x in sys.argv[3].split(',')]: cmd_combo(eta, n)
    elif cmd == 'ledger':
        # per-prime exponents of the primitive denominator D/g of the untwisted form, for a list of n; prints p, n, x=n/p, p mod 3, exponent
        import json
        rows = []
        for n in [int(x) for x in sys.argv[3].split(',')]:
            a, P = partial_fractions(eta, n)
            a0, out = linear_form(eta, n, a, 0)
            coeffs = [a0] + [v for v in out.values() if v != 0]
            D, g = primitive(coeffs)
            f = sympy.factorint(int(D)); fg = sympy.factorint(int(g)) if g > 1 else {}
            for p in sorted(set(f) | set(fg)):
                rows.append(dict(n=n, p=int(p), e=int(f.get(p, 0)) - int(fg.get(p, 0)), x=n/p, pmod3=p % 3))
            print('n=%d done: %s' % (n, ' '.join('%d^%d' % (p, e) for p, e in sorted(f.items()))))
            sys.stdout.flush()
        json.dump(rows, open(sys.argv[4], 'w'))
