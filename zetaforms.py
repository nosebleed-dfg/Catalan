"""zetaforms.py — Zudilin's general construction of linear forms in odd zeta values (Zu04, J. Théor. Nombres Bordeaux 16, §8), exact.

Parameters: odd r, odd q >= r+4, directions eta = (eta0; eta1 <= ... <= eta_q < eta0/2) with sum eta_j <= eta0 (q-r)/2, n >= 1.
  h0 = eta0 n + 2,  h_j = eta_j n + 1.
  R(t) = (h0+2t) * prod_{j<=r} [(t+1)_{h_j-1} (t+1+h0-h_j)_{h_j-1} / (h_j-1)!^2] * prod_{j>r} [(h0-2h_j)! / (t+h_j)_{1+h0-2h_j}]
  F = (1/(r-1)!) sum_{t >= 1-h_1} R^{(r-1)}(t)  =  sum_{m odd, r+2 <= m <= q-2} A_m zeta(m) - A_0        (Lemma 19)
  with A_m = C(m-1, r-1) sum_k B_{m-r+1,k},  A_0 = sum_{k,i} B_{i,k} C(i+r-2, r-1) sum_{l=1}^{k-h_1} l^{-(i+r-1)},
  B_{i,k} the Laurent coefficient of (t+k)^{-i} at the pole -k  (k in [h_{r+1}, h0-h_{r+1}]).
Zudilin's arithmetic (Lemma 19): D_{m_1}^r D_{m_2}...D_{m_{q-r}} Phi^{-1} F in Z zeta(q-2)+...+Z zeta(r+2)+Z, m_0 = max(h_r-1, h0-2h_{r+1}),
  m_j = max(m_0, h0-h_1-h_{r+j}), Phi = prod_{sqrt(h0) < p <= m_{q-r}} p^{nu_p}, nu_p = min_k nu_{k,p}.
Asymptotics (Lemma 20, r=3 stated; used for all r): tau0 = root of (tau-eta0)^r prod(tau-eta_j) - tau^r prod(tau-eta0+eta_j) with Im>0, max Re;
  C0 = -Re f0(tau0);  C2 = r m_1 + m_2 + ... + m_{q-r} - kappa,  kappa = int_0^1 phi dpsi - int_0^{1/m_{q-r}} phi dx/x^2,  phi(x) = min_y phi0(x,y).
usage:  python zetaforms.py forms r ETA n [nocheck]       exact form, true denominators per prime vs Zudilin's bound, NET
        python zetaforms.py asym r ETA                    C0, kappa, C2, net = C2 - C0 (per n)
ETA = eta0,eta1,...,eta_q   e.g.  python zetaforms.py asym 3 91,27,27,27,29,30,31,32,33,34,35,36,37,38
"""
import sys, math, time
from fractions import Fraction
from flint import fmpq, fmpz
import sympy
from mpmath import mp, mpf, mpc, zeta as mzeta, loggamma, exp, log, diff, nsum, inf, fabs, psi as mpsi, polyroots, binomial as mbinom, gamma as mgamma
import numpy as np

class ZetaForm:
    def __init__(self, r, eta, n):
        self.r, self.n = r, n
        self.eta0, self.es = eta[0], list(eta[1:])
        self.q = len(self.es)
        assert r % 2 == 1 and self.q % 2 == 1 and self.q >= r + 4
        assert all(self.es[i] <= self.es[i + 1] for i in range(self.q - 1)) and 2*self.es[-1] < self.eta0
        self.h0 = self.eta0*n + 2
        self.h = [ej*n + 1 for ej in self.es]
        assert sum(self.h) <= self.h0*(self.q - r)//2, 'degree condition fails'
        # numerator factor multiplicities mu[m] for (t+m), m in 1..h0-1
        mu = {}
        for j in range(r):
            hj = self.h[j]
            for m in range(1, hj): mu[m] = mu.get(m, 0) + 1
            for m in range(self.h0 - hj + 1, self.h0): mu[m] = mu.get(m, 0) + 1
        self.mu = mu
        # pole multiplicities s[k] for (t+k), k in [h_{r+1}, h0 - h_{r+1}]
        s = {}
        for j in range(r, self.q):
            hj = self.h[j]
            for k in range(hj, self.h0 - hj + 1): s[k] = s.get(k, 0) + 1
        self.s = s
        # constant
        C = fmpq(1)
        for j in range(r): C /= fmpz(math.factorial(self.h[j] - 1))**2
        for j in range(r, self.q): C *= fmpz(math.factorial(self.h0 - 2*self.h[j]))
        self.C = C
        self.centre = Fraction(self.h0, 2)   # numerator factor 2 (t + h0/2)

    def partial_fractions(self):
        """B[(k, i)] = coefficient of (t+k)^{-i}; the centre factor cancels one order at k = h0/2 (h0 even)."""
        B = {}
        h0 = self.h0
        cen_int = (h0 % 2 == 0)
        for k, sk in self.s.items():
            s_eff = sk - (1 if (cen_int and k == h0//2) else 0) - self.mu.get(k, 0)
            if s_eff <= 0: continue
            S = s_eff
            series = [fmpq(0)]*S; series[0] = fmpq(1)
            Cc = self.C*2
            def mul_pos(ser, alpha, power):
                for _ in range(power):
                    out = [fmpq(0)]*S
                    inva = 1/alpha
                    for j in range(S):
                        out[j] += ser[j]
                        if j + 1 < S: out[j + 1] += ser[j]*inva
                    ser = out
                return ser
            def mul_neg(ser, alpha, power):
                inva = 1/alpha
                pw = [fmpq(1)]
                for t in range(1, S): pw.append(pw[-1]*(-inva))
                for _ in range(power):
                    out = [fmpq(0)]*S
                    for j in range(S):
                        sj = ser[j]
                        if sj == 0: continue
                        for t in range(S - j): out[j + t] += sj*pw[t]
                    ser = out
                return ser
            # centre factor (t + h0/2) -> alpha = h0/2 - k
            alpha = fmpq(h0, 2) - k
            if alpha != 0: Cc *= alpha; series = mul_pos(series, alpha, 1)
            for m, mm in self.mu.items():
                alpha = fmpq(m - k)
                if alpha == 0: continue
                Cc *= alpha**mm; series = mul_pos(series, alpha, mm)
            for k2, sk2 in self.s.items():
                if k2 == k: continue
                alpha = fmpq(k2 - k)
                Cc /= alpha**sk2; series = mul_neg(series, alpha, sk2)
            for i in range(1, S + 1): B[(k, i)] = Cc*series[S - i]
        self.B = B
        return B

    def R_exact(self, t):
        """R(t) exactly at a rational t (Fraction)."""
        tt = fmpq(t.numerator, t.denominator)
        v = self.C*(2*tt + self.h0)
        for m, mm in self.mu.items(): v *= (tt + m)**mm
        for k, sk in self.s.items(): v /= (tt + k)**sk
        return v
    def identity_ok(self, trials=2):
        import random
        for _ in range(trials):
            t = Fraction(random.randint(1, 10**5), random.randint(1, 10**3)) + Fraction(1, 7)
            lhs = self.R_exact(t); tt = fmpq(t.numerator, t.denominator); rhs = fmpq(0)
            for (k, i), b in self.B.items(): rhs += b/(tt + k)**i
            if lhs != rhs: return False
        return True

    def linear_form(self):
        """A_m (m = i + r - 1) and A_0 with Zudilin's partial sums sum_{l=1}^{k-h_1}."""
        r = self.r; h1 = self.h[0]
        A = {}; A0 = fmpq(0)
        for (k, i), b in self.B.items():
            m = i + r - 1
            w = b*fmpz(math.comb(i + r - 2, r - 1))
            A[m] = A.get(m, fmpq(0)) + w
            ps = fmpq(0)
            for l in range(1, k - h1 + 1): ps += fmpq(1, fmpz(l)**m)
            A0 += w*ps
        self.A = {m: v for m, v in A.items() if v != 0}
        self.A0 = A0
        return self.A, A0

    def R_num(self, t):
        """R(t) numerically via loggamma for real t (mp precision set by caller)."""
        h0 = self.h0; r = self.r
        t = mpf(t)
        v = log(h0 + 2*t)
        for j in range(r):
            hj = self.h[j]
            v += loggamma(t + hj) - loggamma(t + 1) - 2*loggamma(hj)
            v += loggamma(t + h0) - loggamma(t + 1 + h0 - hj)
        for j in range(r, self.q):
            hj = self.h[j]
            v += loggamma(h0 - 2*hj + 1) + loggamma(t + hj) - loggamma(t + 1 + h0 - hj)
        return exp(v)

    def needed_dps(self):
        """digits needed to evaluate a0 + sum A_m zeta(m) through the cancellation: size of the coefficients plus the size of F."""
        big = max([(int(a.p).bit_length() - int(a.q).bit_length())*math.log10(2) for a in self.A.values()] + [0.0])
        return int(60 + 2.2*max(big, 0.0)) + 40

    def numeric(self, dps=None, T=400):
        mp.dps = dps or self.needed_dps()
        r = self.r
        f = lambda t: diff(self.R_num, t, r - 1)/math.factorial(r - 1)
        head = sum(f(mpf(t)) for t in range(0, T))
        tail = nsum(lambda t: f(t), [T, inf])
        return head + tail

    def form_value(self, dps=None):
        mp.dps = dps or self.needed_dps()
        v = -mpf(int(self.A0.p))/int(self.A0.q)
        for m, a in self.A.items(): v += mpf(int(a.p))/int(a.q)*mzeta(m)
        return v

    # ---- Zudilin's arithmetic
    def m_params(self):
        r = self.r; h0 = self.h0; h = self.h
        m0 = max(h[r - 1] - 1, h0 - 2*h[r])
        ms = [max(m0, h0 - h[0] - h[r + j - 1]) for j in range(1, self.q - r + 1)]
        return m0, ms
    def nu_kp(self, k, p):
        r = self.r; h0 = self.h0; h = self.h
        v = 0
        for j in range(r):
            v += (k - 1)//p + (h0 - k - 1)//p - (k - h[j])//p - (h0 - h[j] - k)//p - 2*((h[j] - 1)//p)
        for j in range(r, self.q):
            v += (h0 - 2*h[j])//p - (k - h[j])//p - (h0 - h[j] - k)//p
        return v
    def zudilin_bound_exponent(self, p):
        """exponent of p in Zudilin's common denominator D_{m_1}^r D_{m_2}...D_{m_{q-r}} / Phi for the whole form."""
        m0, ms = self.m_params()
        def vD(M): return int(math.floor(math.log(M)/math.log(p))) if M >= p else 0
        e = self.r*vD(ms[0]) + sum(vD(M) for M in ms[1:])
        if p*p > self.h0 and p <= ms[-1]:
            nu = min(self.nu_kp(k, p) for k in range(self.h[self.r], self.h0 - self.h[self.r] + 1))
            e -= max(nu, 0)
        return e

def vp_int(x, p):
    x = abs(int(x)); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v
def vp_fact(N, p):
    v = 0; q = p
    while q <= N: v += N//q; q *= p
    return v

def refined_exponent(Z, p):
    """zeta-analogue of Theorem 1 (finite n): predicted exponent of p in the primitive denominator of the form.
    For each pole k (effective order S) and order i (1..S): v_p(B_{i,k}) >= nu_hat(k) - min(S-i, D_k) [or -(S-i) if some other pole k' has p | k'-k],
    nu_hat(k) = exact valuation of the top coefficient; partial sums cost m*floor(log_p(k-h_1)), m = i+r-1.  Returns (pred, details)."""
    r = Z.r; h0 = Z.h0; h1 = Z.h[0]
    vC = sum(vp_fact(h0 - 2*Z.h[j], p) for j in range(r, Z.q)) - 2*sum(vp_fact(Z.h[j] - 1, p) for j in range(r))
    worst = 10**9
    cen_int = (h0 % 2 == 0)
    for k, sk in Z.s.items():
        S = sk - (1 if (cen_int and k == h0//2) else 0)
        if S <= 0: continue
        nu = vC + vp_int(h0 - 2*k, p) if (h0 - 2*k) != 0 else vC   # centre factor 2(t+h0/2) contributes v_p(h0-2k) (p odd)
        if p == 2: nu = vC + (vp_int(h0 - 2*k, 2) - 1 if h0 - 2*k != 0 else 0)
        Dk = 0
        for m, mm in Z.mu.items():
            v = vp_int(m - k, p)
            if v: nu += mm*v; Dk += mm*v
        if (h0 - 2*k) != 0 and vp_int(h0 - 2*k, p) > (1 if p == 2 else 0): Dk += 1
        den_div = False
        for k2, sk2 in Z.s.items():
            if k2 == k: continue
            v = vp_int(k2 - k, p)
            if v: nu -= sk2*v; den_div = True
        lg = int(math.floor(math.log(k - h1)/math.log(p))) if k - h1 >= p else 0
        for i in range(1, S + 1):
            loss = (S - i) if den_div else min(S - i, Dk)
            m = i + r - 1
            val = nu - loss - m*lg
            worst = min(worst, val)
    return max(0, -worst)

def primitive(coeffs):
    D = fmpz(1)
    for c in coeffs: D = fmpz(int(sympy.ilcm(int(D), int(c.q))))
    g = fmpz(0)
    for c in coeffs: g = fmpz(int(sympy.igcd(int(g), int((c*D).p))))
    return D, g

def cmd_forms(r, eta, n, check=True):
    t0 = time.time()
    Z = ZetaForm(r, eta, n)
    Z.partial_fractions()
    A, A0 = Z.linear_form()
    ms = Z.m_params()[1]
    allA = {}
    for (k, i), b in Z.B.items():
        m = i + r - 1; allA[m] = allA.get(m, fmpq(0)) + b*fmpz(math.comb(i + r - 2, r - 1))
    print('Zudilin zeta forms: r=%d q=%d eta=%s n=%d : h0=%d, poles k in [%d,%d] (%d), max order %d ; zeta slots %s  [%.1fs]' % (
        r, Z.q, eta, n, Z.h0, min(Z.s), max(Z.s), len(Z.s), max(Z.s.values()), sorted(A), time.time() - t0))
    print('  partial fractions reproduce R(t) at random t: %s ; vanishing of A_m for m even and m = r: %s' % (
        Z.identity_ok(), all(allA.get(m, fmpq(0)) == 0 for m in allA if (m % 2 == 0 or m == r))))
    if check:
        fv = Z.form_value(); nv = Z.numeric()
        print('  numeric F = %s ; form = %s ; |diff| = %s ; log|F|/n = %.6f' % (mp.nstr(nv, 15), mp.nstr(fv, 15), mp.nstr(fabs(nv - fv), 3), float(log(fabs(fv)))/n))
    else:
        fv = Z.form_value(); print('  form = %s ; log|F|/n = %.6f' % (mp.nstr(fv, 15), float(log(fabs(fv)))/n))
    coeffs = [A0] + list(A.values())
    D, g = primitive(coeffs)
    hD = math.log(int(D)); hg = math.log(int(g)) if g > 1 else 0.0
    print('  height: log D/n = %.4f, content log g/n = %.4f, primitive log(D/g)/n = %.4f ; NET log|F D/g|/n = %+.4f ; m-parameters %s' % (
        hD/n, hg/n, (hD - hg)/n, (hD - hg)/n + float(log(fabs(fv)))/n, ms))
    fD = sympy.factorint(int(D)); fg = sympy.factorint(int(g)) if g > 1 else {}
    # per-prime comparison with Zudilin's bound
    print('  p : true exponent in D/g | Zudilin bound | refined prediction   (primes up to m_1 = %d; refined is a lower bound on v_p, i.e. upper bound on the exponent)' % ms[0])
    zb_total = 0.0; true_total = 0.0; rf_total = 0.0; viol = []; single = []
    rows = []
    for p in sympy.primerange(2, ms[0] + 1):
        te = int(fD.get(p, 0)) - int(fg.get(p, 0)); zb = Z.zudilin_bound_exponent(p); rf = refined_exponent(Z, p)
        rows.append((p, te, zb, rf)); zb_total += zb*math.log(p); true_total += te*math.log(p); rf_total += rf*math.log(p)
        if te > rf: viol.append(p)
        if p*p > Z.h0: single.append((p, te, zb, rf))
    print('  ' + '  '.join('%d:%d|%d|%d' % row for row in rows))
    print('  per n: Zudilin bound %.4f ; refined %.4f ; true %.4f   (all primes <= m_1)' % (zb_total/n, rf_total/n, true_total/n))
    sz = sum(r_[2]*math.log(r_[0]) for r_ in single); sr = sum(r_[3]*math.log(r_[0]) for r_ in single); st = sum(r_[1]*math.log(r_[0]) for r_ in single)
    print('  single-digit primes (p^2 > h0): Zudilin %.4f ; refined %.4f ; true %.4f per n ; refined exact at %d/%d primes ; VIOLATIONS (true > refined): %s' % (
        sz/n, sr/n, st/n, sum(1 for r_ in single if r_[1] == r_[3]), len(single), viol))
    extra = {p: e for p, e in fD.items() if p > ms[0]}
    if extra: print('  primes above m_1 in D:', extra)
    return Z

def asym(r, eta, verbose=True):
    eta0, es = eta[0], list(eta[1:]); q = len(es)
    mp.dps = 30
    # saddle point: roots of P(tau) = (tau-eta0)^r prod(tau-eta_j) - tau^r prod(tau-eta0+eta_j)
    import sympy as sp
    tau = sp.symbols('tau')
    P = sp.expand((tau - eta0)**r*sp.prod([tau - e for e in es]) - tau**r*sp.prod([tau - eta0 + e for e in es]))
    coeffs = [complex(c) for c in sp.Poly(P, tau).all_coeffs()]
    roots = np.roots(coeffs)
    cand = [z for z in roots if z.imag > 1e-9]
    tau0 = max(cand, key=lambda z: z.real)
    tau0 = mpc(tau0.real, tau0.imag)
    def f0(t):
        v = r*eta0*log(eta0 - t)
        for e in es: v += e*log(t - e) - (eta0 - e)*log(t - eta0 + e)
        for e in es[:r]: v -= 2*e*log(e)
        for e in es[r:]: v += (eta0 - 2*e)*log(eta0 - 2*e)
        return v
    C0 = -float(f0(tau0).real)
    # m parameters (scaled)
    m0 = max(es[r - 1], eta0 - 2*es[r])
    ms = [max(m0, eta0 - es[0] - es[r + j - 1]) for j in range(1, q - r + 1)]
    # phi(x) = min_y phi0(x,y), evaluated on the Farey grid at midpoints
    def phi0(x, y):
        v = np.zeros_like(y)
        for e in es[:r]:
            v += np.floor(y) + np.floor(eta0*x - y) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y) - 2*math.floor(e*x)
        for e in es[r:]:
            v += math.floor((eta0 - 2*e)*x) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y)
        return v
    def phi(x):
        b = [0.0, (eta0*x) % 1] + [(e*x) % 1 for e in es] + [((eta0 - e)*x) % 1 for e in es]
        b = np.unique(np.array(b)); b2 = np.append(b, b[0] + 1); mids = (b2[:-1] + b2[1:])/2
        return int(phi0(x, mids).min())
    Dg = 2*eta0
    pts = sorted(set(Fraction(a, d) for d in range(1, Dg + 1) for a in range(0, d + 1)))
    kap = mpf(0); steps = []
    lo = Fraction(1, ms[-1])
    for a, b in zip(pts[:-1], pts[1:]):
        v = phi(float((a + b)/2))
        if v:
            fa, fb = mpf(a.numerator)/a.denominator, mpf(b.numerator)/b.denominator
            kap += v*(mpsi(0, fb) - mpsi(0, fa))
            if a < lo:
                bb = min(fb, mpf(lo.numerator)/lo.denominator)
                kap -= v*(1/fa - 1/bb)
            steps.append((a, b, v))
    C2 = r*ms[0] + sum(ms[1:]) - float(kap)
    if verbose:
        comp = []
        for a, b, v in steps:
            if comp and comp[-1][2] == v and comp[-1][1] == a: comp[-1] = (comp[-1][0], b, v)
            else: comp.append((a, b, v))
        print('zeta forms asym: r=%d q=%d eta=%s : slots zeta(%s)' % (r, q, eta, ','.join(str(m) for m in range(r + 2, q - 1, 2))))
        print('  tau0 = %s ; C0 (closeness) = %.8f' % (mp.nstr(tau0, 12), C0))
        print('  m = %s ; kappa = %.8f ; C2 (height) = %d - kappa = %.8f ; net C2 - C0 = %+.6f per n' % (ms, float(kap), r*ms[0] + sum(ms[1:]), C2, C2 - C0))
        print('  phi steps: ' + ' '.join('[%s,%s)=%d' % (a, b, v) for a, b, v in comp))
    return C0, float(kap), C2

def refined_rebate(r, eta, verbose=False):
    """Asymptotic rebate of the refined (pole-coupled) exponent below Zudilin's bound, per n, over the single-digit range x = n/p in (0, 1]:
    E_Z(x) = r 1[x >= 1/m_1] + sum_{j>=2} 1[x >= 1/m_j] - phi(x),
    E_ref(x) = max(0, max_y { S(y)+r-1 - phi0(x,y) if the pole is in the window (y - eta_1 x >= 1) else (S(y)-1 if den_div else min(S(y)-1, D(x,y))) - phi0(x,y) }),
    rebate = int_0^1 (E_Z - E_ref) dx/x^2 (E_Z >= E_ref pointwise up to the centre term)."""
    eta0, es = eta[0], list(eta[1:]); q = len(es)
    m0 = max(es[r - 1], eta0 - 2*es[r]); ms = [max(m0, eta0 - es[0] - es[r + j - 1]) for j in range(1, q - r + 1)]
    def pieces(x):
        # breakpoints in y within the pole range [eta_{r+1} x, (eta0 - eta_{r+1}) x]
        lo_y, hi_y = es[r]*x, (eta0 - es[r])*x
        b = set([lo_y, hi_y, es[0]*x + 1.0])
        for e in es:
            for base in (e*x, (eta0 - e)*x):
                t0 = math.floor(lo_y - base) - 1
                while base + t0 <= hi_y + 1:
                    if lo_y <= base + t0 <= hi_y: b.add(base + t0)
                    t0 += 1
        for base in (0.0, eta0*x):   # floor(y), floor(eta0 x - y)
            t0 = math.floor(lo_y - base) - 1
            while base + t0 <= hi_y + 1:
                if lo_y <= base + t0 <= hi_y: b.add(base + t0)
                t0 += 1
        b = np.array(sorted(b))
        return (b[:-1] + b[1:])/2
    def phi0(x, y):
        v = np.zeros_like(y)
        for e in es[:r]:
            v += np.floor(y) + np.floor(eta0*x - y) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y) - 2*math.floor(e*x)
        for e in es[r:]:
            v += math.floor((eta0 - 2*e)*x) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y)
        return v
    def Eref(x):
        y = pieces(x)
        if len(y) == 0: return 0.0
        S = np.zeros_like(y); den = np.zeros_like(y, dtype=bool)
        for e in es[r:]:
            cov = (y >= e*x) & (y <= (eta0 - e)*x)
            S += cov
            den |= cov & ((np.floor((eta0 - e)*x - y) + np.floor(y - e*x)) >= 1)
        D = np.zeros_like(y)
        for e in es[:r]:
            D += np.floor(y) - np.floor(y - e*x) + np.floor(eta0*x - y) - np.floor((eta0 - e)*x - y)
        # exact top valuation: Zudilin's phi0 overstates by one per denominator block NOT covering the pole
        ph = phi0(x, y) - ((q - r) - S)
        win = (y - es[0]*x) >= 1.0
        cost = np.where(win, S + r - 1 - ph, np.where(den, S - 1 - ph, np.minimum(S - 1, D) - ph))
        cost = np.where(S >= 1, cost, -1e9)
        return float(cost.max())          # may be negative: content of the whole form
    def EZ(x, his=False):
        """Zudilin's bound exponent: D-count minus phi = min_y phi0; his convention applies phi only for p <= m_{q-r} n (x >= 1/m_{q-r})."""
        y = pieces(x)
        ph = phi0(x, y); S = np.zeros_like(y)
        for e in es[r:]: S += (y >= e*x) & (y <= (eta0 - e)*x)
        phmin = float(ph[S >= 1].min()) if (S >= 1).any() else 0.0
        dcount = r*(1 if x >= 1/ms[0] else 0) + sum(1 for M in ms[1:] if x >= 1/M)
        if his and x < 1/ms[-1]: return float(dcount)
        return dcount - phmin
    Dg = 3*eta0
    pts = sorted(set(Fraction(a, d) for d in range(1, Dg + 1) for a in range(0, d + 1)))
    reb = 0.0; hZ = 0.0; hR = 0.0; hH = 0.0; steps = []
    for a, b in zip(pts[:-1], pts[1:]):
        if a == 0: continue
        x = float((a + b)/2)
        eh, ez, er = EZ(x, his=True), EZ(x), Eref(x)
        w = float(1/a - 1/b)
        reb += (eh - er)*w; hZ += ez*w; hR += er*w; hH += eh*w
        if eh != er: steps.append((a, b, eh, er))
    if verbose:
        comp = []
        for a, b, ez, er in steps:
            if comp and comp[-1][2:] == (ez, er) and comp[-1][1] == a: comp[-1] = (comp[-1][0], b, ez, er)
            else: comp.append((a, b, ez, er))
        print('  exponent functions on x in (0,1], integrated dx/x^2: Zudilin (his convention) %.4f ; Zudilin with phi for all p %.4f ; refined %.4f per n ; rebate (his -> refined) %.4f per n' % (hH, hZ, hR, reb))
        print('  steps where his and refined differ (x-interval: his -> refined): ' + ' '.join('[%s,%s): %g->%g' % (a, b, ez, er) for a, b, ez, er in comp[:40]) + (' ...' if len(comp) > 40 else ''))
    return reb, hH, hR

def net_of(r, eta, refined=False):
    try:
        C0, kap, C2 = asym(r, eta, verbose=False)
        net = C2 - C0
        if refined: net -= refined_rebate(r, eta)[0]
        return net
    except Exception:
        return None

def valid(r, eta):
    eta0, es = eta[0], eta[1:]; q = len(es)
    return all(1 <= es[i] <= es[i + 1] for i in range(q - 1)) and 2*es[-1] < eta0 and sum(es) <= eta0*(q - r)//2

def search(r, q, eta0s, refined=False, seed=None, max_sweeps=12):
    """coordinate descent over eta for fixed r, q at each eta0: single +-1 moves and paired moves; objective net = C2 - C0 (Zudilin's bound) or the refined net.
    seed: an explicit starting eta (its eta0 is replaced by each eta0 in eta0s, other entries scaled)."""
    best_all = None
    for eta0 in eta0s:
        if seed is not None:
            es = sorted(max(1, int(round(e*eta0/seed[0]))) for e in seed[1:])
        else:
            # seed: numerator directions ~0.30 eta0, denominator directions spread 0.32..0.42 eta0
            es = [max(1, int(round(0.30*eta0)))]*r + [max(1, int(round((0.32 + 0.10*i/max(1, q - r - 1))*eta0))) for i in range(q - r)]
            es = sorted(es)
        while sum(es) > eta0*(q - r)//2: es[es.index(max(es))] -= 1; es = sorted(es)
        while 2*es[-1] >= eta0: es[-1] -= 1; es = sorted(es)
        cur = [eta0] + es; cv = net_of(r, cur, refined)
        if cv is None:
            print('eta0=%d: seed failed' % eta0); continue
        improved = True; sweeps = 0
        while improved and sweeps < max_sweeps:
            improved = False; sweeps += 1
            for j in range(1, q + 1):
                for dl in (1, -1):
                    cand = cur[:]; cand[j] += dl; cand = [cand[0]] + sorted(cand[1:])
                    if not valid(r, cand): continue
                    v = net_of(r, cand, refined)
                    if v is not None and v < cv - 1e-9: cur, cv, improved = cand, v, True
            for j in range(1, q + 1):
                for l in range(1, q + 1):
                    if j == l: continue
                    cand = cur[:]; cand[j] += 1; cand[l] -= 1; cand = [cand[0]] + sorted(cand[1:])
                    if not valid(r, cand): continue
                    v = net_of(r, cand, refined)
                    if v is not None and v < cv - 1e-9: cur, cv, improved = cand, v, True
        C0, kap, C2 = asym(r, cur, verbose=False)
        reb = refined_rebate(r, cur)[0]
        print('r=%d q=%d eta0=%3d: best %snet %+.4f per n (%+.4f per unit eta0) at eta=%s ; C0 %.3f, C2 %.3f (kappa %.3f), Zudilin net %+.4f, refined rebate %.4f, refined net %+.4f' % (
            r, q, eta0, 'refined ' if refined else '', cv, cv/eta0, cur, C0, C2, kap, C2 - C0, reb, C2 - C0 - reb))
        sys.stdout.flush()
        if best_all is None or cv < best_all[0]: best_all = (cv, cur)
    print('best: net %+.4f at eta=%s' % best_all)

if __name__ == '__main__':
    cmd = sys.argv[1]; r = int(sys.argv[2])
    if cmd == 'search':
        q = int(sys.argv[3]); eta0s = [int(x) for x in sys.argv[4].split(',')]
        seed = None; sweeps = 12
        for a in sys.argv[5:]:
            if a.startswith('seed='): seed = [int(x) for x in a[5:].split(',')]
            if a.startswith('sweeps='): sweeps = int(a[7:])
        search(r, q, eta0s, refined=('ref' in sys.argv), seed=seed, max_sweeps=sweeps); sys.exit()
    eta = [int(x) for x in sys.argv[3].split(',')]
    if cmd == 'forms':
        n = int(sys.argv[4]); cmd_forms(r, eta, n, check=('nocheck' not in sys.argv))
    elif cmd == 'asym':
        C0, kap, C2 = asym(r, eta)
        reb, hZ, hR = refined_rebate(r, eta, verbose=True)
        print('  REFINED: height %.6f ; net %+.6f per n   (Zudilin net %+.6f)' % (C2 - reb, C2 - reb - C0, C2 - C0))
