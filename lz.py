"""lz.py — Lai–Zhou's two-value construction (arXiv 2103.00904) exact, with the two eliminated zeta values fixed.

Parameters (s odd >= 5; m1, m2 >= 1; 0 <= delta_j < m2/2, sum delta_j < ((s-2)m2 - 8m1)/2; n even):
  L = (2m1+m2) n,  R_n(t) = N (2t+m2 n)(t-m1 n) prod_{theta in {1,1/2,1/3,2/3}} (t-m1 n+theta)_L / [ (t)_{m2 n+1} prod_j (t+delta_j n)_{(m2-2delta_j)n+1} ],
  N = 2^{2L} 3^{3L} prod_j ((m2-2delta_j)n)! / n!^{8m1+3m2}.
  S_theta = sum_{t>=1} R_n(t+theta) = rho_{0,theta} + sum_{i odd>=3} rho_i zeta(i,theta),  rho_i independent of theta;
  S^_b = sum_{k=1}^b S_{k/b} = rho^_{0,b} + sum_i rho_i b^i zeta(i);  S~ = sum_b w_b S^_b with w orthogonal to (1,2^i,3^i) for i = i1, i2.
  With (i1,i2) = (3,5): w = (45,-9,1), S~ in Q + sum_{i odd >= 7} Q zeta(i): net < 0 proves "one of zeta(7),...,zeta(s) is irrational".
Asymptotics (their Lemmas 3.3, 4.1): C1 = int_0^1 nu0 dpsi - int_0^{1/M'} nu0 dx/x^2 (the Phi_n rate, a gain), C2 = sum_j max(m2-2 delta_min, m2-delta_j) + log g(x0);
  net = C2 - C1 (needs < 0).  Their s=35: m1=209, m2=243, delta = 4^5, 5..10, 12..52 step 2, 56,60,64,68: x0 = 2.89493833, C1 = 16779.9312, C2 = 16779.2826.
usage:  python lz.py asym s m1 m2 d1,...,d_{s+1}
        python lz.py forms s m1 m2 d1,...,d_{s+1} n [nocheck]
"""
import sys, math, time
from fractions import Fraction
from flint import fmpq, fmpz
import sympy
import numpy as np
from mpmath import mp, mpf, zeta as mzeta, loggamma, exp, log, nsum, inf, fabs, psi as mpsi

def fr2q(f): return fmpq(f.numerator, f.denominator)
def q2mp(v): return mpf(int(v.p))/mpf(int(v.q))

THETAS = [Fraction(1), Fraction(1, 2), Fraction(1, 3), Fraction(2, 3)]

class LZ:
    def __init__(self, s, m1, m2, deltas, n):
        assert s % 2 == 1 and n % 2 == 0 and len(deltas) == s + 1
        assert all(0 <= d and 2*d < m2 for d in deltas) and 2*sum(deltas) < (s - 2)*m2 - 8*m1
        self.s, self.m1, self.m2, self.n = s, m1, m2, n
        self.deltas = sorted(deltas)
        self.L = (2*m1 + m2)*n
        # numerator factors (t + c) with multiplicity 1: integer c in [-m1 n, -1] U [m2 n + 1, (m1+m2) n]; rational c = theta + l - m1 n for theta in {1/2,1/3,2/3}, l < L
        num = [Fraction(c) for c in range(-m1*n, 0)] + [Fraction(c) for c in range(m2*n + 1, (m1 + m2)*n + 1)]
        for th in THETAS[1:]:
            num += [th + l - m1*n for l in range(self.L)]
        self.num = num
        self.centre = Fraction(m2*n, 2)      # factor 2 (t + m2 n/2)
        # poles: k in [delta_j n, (m2 - delta_j) n] per block
        poles = {}
        for d in self.deltas:
            for k in range(d*n, (m2 - d)*n + 1): poles[k] = poles.get(k, 0) + 1
        self.poles = poles
        # constant N
        C = fmpq(fmpz(2)**(2*self.L)*fmpz(3)**(3*self.L))
        for d in self.deltas: C *= fmpz(math.factorial((m2 - 2*d)*n))
        C /= fmpz(math.factorial(n))**(8*m1 + 3*m2)
        self.C = C

    def partial_fractions(self):
        a = {}
        for k, sk in self.poles.items():
            S = sk - (1 if Fraction(k) == self.centre else 0)
            zeros_here = sum(1 for c in self.num if c == k)   # numerator zero at the pole (does not happen: integer zeros avoid [0, m2 n])
            S -= zeros_here
            if S <= 0: continue
            series = [fmpq(0)]*S; series[0] = fmpq(1)
            Cc = self.C*2
            def mul_pos(ser, alpha):
                out = [fmpq(0)]*S; inva = 1/alpha
                for j in range(S):
                    out[j] += ser[j]
                    if j + 1 < S: out[j + 1] += ser[j]*inva
                return out
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
            alpha = fr2q(self.centre) - k
            if alpha != 0: Cc *= alpha; series = mul_pos(series, alpha)
            for c in self.num:
                alpha = fr2q(c) - k
                if alpha == 0: continue
                Cc *= alpha; series = mul_pos(series, alpha)
            for k2, sk2 in self.poles.items():
                if k2 == k: continue
                alpha = fmpq(k2 - k)
                Cc /= alpha**sk2; series = mul_neg(series, alpha, sk2)
            for i in range(1, S + 1): a[(k, i)] = Cc*series[S - i]
        self.a = a
        return a

    def R_exact(self, t):
        tt = fr2q(t)
        v = self.C*2*(tt + fr2q(self.centre))
        for c in self.num: v *= (tt + fr2q(c))
        for k, sk in self.poles.items(): v /= (tt + k)**sk
        return v
    def identity_ok(self, trials=2):
        import random
        for _ in range(trials):
            t = Fraction(random.randint(1, 10**5), random.randint(1, 10**3)) + Fraction(1, 7)
            lhs = self.R_exact(t); tt = fr2q(t); rhs = fmpq(0)
            for (k, i), c in self.a.items(): rhs += c/(tt + k)**i
            if lhs != rhs: return False
        return True

    def forms(self):
        """rho_i (theta-independent) and rho_{0,theta} for the four thetas; S^_b and the (3,5)-combination."""
        rho = {}
        for (k, i), c in self.a.items(): rho[i] = rho.get(i, fmpq(0)) + c
        rho0 = {}
        for th in THETAS:
            thq = fr2q(th); tot = fmpq(0)
            # rho_{0,theta} = - sum_{k,i} a_{i,k} sum_{l=0}^{k} (l + theta)^{-i}
            for (k, i), c in self.a.items():
                ps = fmpq(0)
                for l in range(0, k + 1): ps += 1/(l + thq)**i
                tot -= c*ps
            rho0[th] = tot
        # S^_b constant parts
        hat0 = {1: rho0[Fraction(1)], 2: rho0[Fraction(1, 2)] + rho0[Fraction(1)], 3: rho0[Fraction(1, 3)] + rho0[Fraction(2, 3)] + rho0[Fraction(1)]}
        self.rho, self.rho0, self.hat0 = rho, rho0, hat0
        return rho, rho0, hat0

    def combination(self, i1=3, i2=5):
        w = (2**i1*3**i2 - 3**i1*2**i2, 3**i1 - 3**i2, 2**i2 - 2**i1)
        g = math.gcd(math.gcd(abs(w[0]), abs(w[1])), abs(w[2])); w = tuple(x//g for x in w)
        assert w[0] + 2*w[1] + 3*w[2] != 0
        const = sum(fmpz(w[b - 1])*self.hat0[b] for b in (1, 2, 3))
        coef = {}
        for i, r in self.rho.items():
            cw = sum(w[b - 1]*b**i for b in (1, 2, 3))
            if cw != 0 and r != 0: coef[i] = r*fmpz(cw)
        self.w = w
        return const, coef

    def R_num(self, t):
        """R(t) numerically for real t > 0 via loggamma (mp precision set by caller)."""
        m1, m2, n, L = self.m1, self.m2, self.n, self.L
        t = mpf(t)
        v = log(2*t + m2*n) + log(t - m1*n) if t > m1*n else None
        # use signed product: compute log|R| and sign separately
        sign = 1
        def lg(x):
            nonlocal sign
            if x < 0: sign = -sign
            return log(fabs(x))
        v = lg(2*t + m2*n) + lg(t - m1*n)
        for th in THETAS:
            thf = mpf(th.numerator)/th.denominator
            v += loggamma(t - m1*n + thf + L) - loggamma(t - m1*n + thf)   # (t - m1 n + theta)_L, all factors positive for t > m1 n - theta... handle sign below
        # the theta = 1 block and (t)_{m2 n + 1}: for t >= 1 non-integer shifts, t - m1 n + 1 + l may be negative for small t: use loggamma of positive shifted args instead:
        return v, sign

    def numeric(self, dps=None):
        """direct numerical S~ = sum_b w_b sum_{k<=b} sum_{t>=1} R(t + k/b), summing R exactly at rational points (fmpq) converted to mpf."""
        mp.dps = dps or self.needed_dps()
        tot = mpf(0)
        T = 60 + 3*self.n
        for b in (1, 2, 3):
            for kk in range(1, b + 1):
                th = Fraction(kk, b)
                sub = mpf(0)
                for t in range(1, T):
                    sub += q2mp(self.R_exact(t + th))
                # tail: geometric-like decay ~ t^{deg}; use nsum on floats of the exact function
                tail = nsum(lambda tt: q2mp(self.R_exact(Fraction(int(tt)) + th)), [T, inf])
                tot += self.w[b - 1]*(sub + tail)
        return tot

    def needed_dps(self):
        big = max([(int(c.p).bit_length() - int(c.q).bit_length())*math.log10(2) for c in list(self.coef.values()) + [self.const]] + [0.0])
        return int(60 + 2.2*max(big, 0.0)) + 40

    def form_value(self, dps=None):
        mp.dps = dps or self.needed_dps()
        v = q2mp(self.const)
        for i, c in self.coef.items(): v += q2mp(c)*mzeta(i)
        return v

def vp_int(x, p):
    x = abs(int(x)); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v
def vp_fact(N, p):
    v = 0; q = p
    while q <= N: v += N//q; q *= p
    return v
def count_ap(lo, hi, resid, p):
    """number of integers c in [lo, hi] with c = resid (mod p)"""
    if hi < lo: return 0
    first = lo + ((resid - lo) % p)
    return 0 if first > hi else (hi - first)//p + 1

def refined_exponent_lz(Z, p):
    """Predicted exponent of p in the primitive denominator of S~ (single-digit primes p > 3, p^2 > 3L): exact top valuation nu_hat(k),
    saturation min(S-i, D_k) (S-i if another pole is = k mod p), windows for the theta partial sums (k >= l0(theta)),
    and Lai-Zhou's cross-pole cap for theta != 1: v_p(rho_{0,theta}) >= nu_LZ - (s+1) 1[p <= M]."""
    s, m1, m2, n, L = Z.s, Z.m1, Z.m2, Z.n, Z.L
    M = (m2 - 2*Z.deltas[0])*n
    vN = sum(vp_fact((m2 - 2*d)*n, p) for d in Z.deltas) - (8*m1 + 3*m2)*vp_fact(n, p)
    inv2 = pow(2, -1, p); inv3 = pow(3, -1, p)
    # numerator content at pole k: integer part c in [-m1 n, -1] U [m2 n + 1, (m1+m2) n] with c = k (mod p); theta blocks: l with theta + l - m1 n - k = 0 (mod p), l in [0, L)
    def Dk(k):
        d = count_ap(-m1*n, -1, k % p, p) + count_ap(m2*n + 1, (m1 + m2)*n, k % p, p)
        for th in THETAS[1:]:
            a, b = th.numerator, th.denominator
            inv = inv2 if b == 2 else inv3
            # a/b + l - m1 n - k = 0 mod p  <=>  l = m1 n + k - a*inv (mod p)
            d += count_ap(0, L - 1, (m1*n + k - a*inv) % p, p)
        return d
    poles = Z.poles
    worst_coef = 10**9; worst_c1 = 10**9; worst_cth = {th: 10**9 for th in THETAS[1:]}
    nuLZ = 10**9
    for k, sk in poles.items():
        S = sk - (1 if Fraction(k) == Z.centre else 0)
        if S <= 0: continue
        D = Dk(k)
        cen = vp_int(m2*n - 2*k, p) if (m2*n - 2*k) != 0 else 0
        den = 0; den_div = False
        for k2, sk2 in poles.items():
            if k2 == k: continue
            v = vp_int(k2 - k, p)
            if v: den += sk2*v; den_div = True
        nu_hat = vN + D + cen - den
        # Lai-Zhou nu (their convention: full loss s+1-i; overstates by one per non-covering block): nu_LZ(k) = nu_hat + (s+1 - sk) (- centre term)
        nuLZ = min(nuLZ, nu_hat - cen + (s + 1 - sk))
        for i in range(1, S + 1):
            loss = (S - i) if den_div else min(S - i, D + (1 if cen else 0))
            Lik = nu_hat - loss
            worst_coef = min(worst_coef, Lik)
            worst_c1 = min(worst_c1, Lik - i*(1 if k >= p - 1 else 0))
            for th in THETAS[1:]:
                a, b = th.numerator, th.denominator
                l0 = ((-a)*(inv2 if b == 2 else inv3)) % p
                worst_cth[th] = min(worst_cth[th], Lik - i*(1 if k >= l0 else 0))
    # LZ cap for theta != 1
    cap = nuLZ - (s + 1)*(1 if p <= M else 0)
    worst = min([worst_coef, worst_c1] + [max(worst_cth[th], cap) for th in THETAS[1:]])
    return -worst

def lz_bound_exponent(Z, p):
    """Lai-Zhou's exponent: sum_j 1[p <= max(M, (m2-delta_j)n)+1] - nu0(n/p) 1[sqrt(3L) < p <= M]."""
    s, m2, n = Z.s, Z.m2, Z.n
    M = (m2 - 2*Z.deltas[0])*n
    e = sum(1 for d in Z.deltas if p <= max(M, (m2 - d)*n) + 1)
    if p*p > 3*Z.L and p <= M:
        nu = 10**9
        m1 = Z.m1
        for k in Z.poles:
            v = sum(math.floor((m2 - 2*d)*n/p) - math.floor((k - d*n)/p) - math.floor(((m2 - d)*n - k)/p) for d in Z.deltas)
            v += math.floor((2*m1*n + 2*k)/p) - math.floor((m1*n + k)/p) + math.floor((2*(m1 + m2)*n - 2*k)/p) - math.floor(((m1 + m2)*n - k)/p)
            v += math.floor((3*m1*n + 3*k)/p) + math.floor((3*(m1 + m2)*n - 3*k)/p) - math.floor(k/p) - math.floor((m2*n - k)/p) - (8*m1 + 3*m2)*math.floor(n/p)
            nu = min(nu, v)
        e -= max(nu, 0)
    return e

def lz_exponents(s, m1, m2, deltas, x):
    """Asymptotic exponent functions at x = n/p (single-digit range x <= 1, p > 3):
    E_LZ(x) = sum_j 1[x >= 1/max(M', m2-delta_j)] - nu0(x) 1[x >= 1/M'],
    E_ref(x) = max over coefficient types of the pole-coupled costs (numerator content credited, saturation, windows, LZ cross-pole cap for theta != 1).
    Returns (E_LZ, E_ref, nu0)."""
    deltas = sorted(deltas); dmin = deltas[0]; Mp = m2 - 2*dmin; q = len(deltas)
    lo_y, hi_y = dmin*x, (m2 - dmin)*x
    # y breakpoints inside the pole range
    b = set([lo_y, hi_y, 1.0, 0.5, 1.0/3, 2.0/3, dmin*x + 1.0, (m2 - dmin)*x - 1.0])
    def add_shifts(base, mult=1):
        # points where mult*(y - base) or mult*(base - y) is an integer: y = base + t/mult
        t0 = math.floor((lo_y - base)*mult) - 1
        while base + t0/mult <= hi_y + 1e-12:
            yv = base + t0/mult
            if lo_y - 1e-12 <= yv <= hi_y + 1e-12: b.add(yv)
            t0 += 1
    add_shifts(0.0); add_shifts(m2*x)
    for mult in (1, 2, 3): add_shifts(-m1*x, mult); add_shifts((m1 + m2)*x, mult)
    for d in deltas: add_shifts(d*x); add_shifts((m2 - d)*x)
    b = np.array(sorted(v for v in b if lo_y - 1e-12 <= v <= hi_y + 1e-12))
    if len(b) < 2: return 0.0, 0.0, 0.0
    y = (b[:-1] + b[1:])/2
    A = m1*x + y; B = (m1 + m2)*x - y
    D = (np.floor(2*A) - np.floor(A) + np.floor(2*B) - np.floor(B)) + (np.floor(3*A) + np.floor(3*B)) - np.floor(y) - np.floor(m2*x - y) - (8*m1 + 3*m2)*math.floor(x)
    S = np.zeros_like(y); nuLZ = D.copy()
    for d in deltas:
        cov = (y >= d*x) & (y <= (m2 - d)*x)
        S += cov
        nuLZ += math.floor((m2 - 2*d)*x) - np.floor(y - d*x) - np.floor((m2 - d)*x - y)
    nuhat = nuLZ - (q - S)
    den_div = ((y - dmin*x) >= 1.0) | (((m2 - dmin)*x - y) >= 1.0)
    valid = S >= 1
    nu0 = float(nuLZ[valid].min()) if valid.any() else 0.0
    lossmax = np.where(den_div, S - 1, np.minimum(S - 1, D))
    cost_coef = lossmax - nuhat
    def cost_window(thr):
        win = y >= thr
        return np.where(win, S - nuhat, lossmax - nuhat)
    c1 = cost_window(1.0); c2 = cost_window(0.5); c3 = cost_window(1.0/3)
    big = -1e9
    m_coef = float(np.where(valid, cost_coef, big).max())
    m1c = float(np.where(valid, c1, big).max())
    cap = (q if x >= 1/Mp else 0) - nu0     # LZ: v_p(rho_{0,theta}) >= nu0 - (s+1) 1[p <= M]  ->  exponent <= (s+1)1[...] - nu0
    m2c = min(float(np.where(valid, c2, big).max()), cap)
    m3c = min(float(np.where(valid, c3, big).max()), cap)
    E_ref = max(m_coef, m1c, m2c, m3c)
    E_LZ = sum(1 for d in deltas if x >= 1/max(Mp, m2 - d)) - (nu0 if x >= 1/Mp else 0)
    return E_LZ, E_ref, nu0

def lz_refined(s, m1, m2, deltas, du=0.02, umax=None, verbose=True):
    """integrate E_LZ and E_ref over x in (0,1] via u = 1/x (uniform grid): heights per n on the single-digit range, and the rebate."""
    deltas = sorted(deltas); dmin = deltas[0]
    if umax is None: umax = (m2 - dmin) + 2.0
    us = np.arange(1.0 + du/2, umax, du)
    hL = 0.0; hR = 0.0; steps = []
    for u in us:
        x = 1.0/u
        eL, eR, _ = lz_exponents(s, m1, m2, deltas, x)
        hL += eL*du; hR += eR*du
        steps.append((u, eL, eR))
    reb = hL - hR
    if verbose:
        print('  single-digit range (p >= n): LZ exponent integral %.3f, refined %.3f per n ; rebate %.3f per n' % (hL, hR, reb))
        # summarise where they differ, by u-intervals
        diff = [(u, eL, eR) for u, eL, eR in steps if abs(eL - eR) > 1e-9]
        if diff:
            comp = []
            for u, eL, eR in diff:
                if comp and comp[-1][2:] == (eL, eR) and abs(comp[-1][1] - (u - du)) < 1e-9: comp[-1] = (comp[-1][0], u, eL, eR)
                else: comp.append((u, u, eL, eR))
            print('  differences (p/n interval: LZ -> refined): ' + ' '.join('[%.2f,%.2f] %g->%g' % (a, bb, eL, eR) for a, bb, eL, eR in comp[:30]) + (' ...' if len(comp) > 30 else ''))
    return reb, hL, hR

def primitive(coeffs):
    D = fmpz(1)
    for c in coeffs: D = fmpz(int(sympy.ilcm(int(D), int(c.q))))
    g = fmpz(0)
    for c in coeffs: g = fmpz(int(sympy.igcd(int(g), int((c*D).p))))
    return D, g

def cmd_forms(s, m1, m2, deltas, n, check=True):
    t0 = time.time()
    Z = LZ(s, m1, m2, deltas, n)
    Z.partial_fractions()
    ok = Z.identity_ok()
    rho, rho0, hat0 = Z.forms()
    const, coef = Z.combination(3, 5)
    Z.const, Z.coef = const, coef
    even_zero = all(rho.get(i, fmpq(0)) == 0 for i in range(2, s + 2, 2))
    print('LZ forms: s=%d m1=%d m2=%d deltas=%s n=%d : L=%d, poles k in [%d,%d] (%d), max order %d ; identity %s ; even rho vanish %s ; w=%s ; remaining zeta slots %s  [%.1fs]' % (
        s, m1, m2, Z.deltas, n, Z.L, min(Z.poles), max(Z.poles), len(Z.poles), max(Z.poles.values()), ok, even_zero, Z.w, sorted(coef), time.time() - t0))
    fv = Z.form_value()
    if check:
        nv = Z.numeric()
        print('  numeric S~ = %s ; form = %s ; |diff| = %s ; log|S~|/n = %.6f' % (mp.nstr(nv, 15), mp.nstr(fv, 15), mp.nstr(fabs(nv - fv), 3), float(log(fabs(fv)))/n))
    else:
        print('  form = %s ; log|S~|/n = %.6f' % (mp.nstr(fv, 15), float(log(fabs(fv)))/n))
    coeffs = [const] + list(coef.values())
    D, g = primitive(coeffs)
    hD = math.log(int(D)); hg = math.log(int(g)) if g > 1 else 0.0
    print('  height: log D/n = %.4f, content log g/n = %.4f, primitive log(D/g)/n = %.4f ; NET log|S~ D/g|/n = %+.4f' % (hD/n, hg/n, (hD - hg)/n, (hD - hg)/n + float(log(fabs(fv)))/n))
    fD = sympy.factorint(int(D)); fg = sympy.factorint(int(g)) if g > 1 else {}
    print('  D/g factorization: ' + ' '.join('%d^%d' % (p, int(fD.get(p, 0)) - int(fg.get(p, 0))) for p in sorted(set(fD) | set(fg)) if int(fD.get(p, 0)) - int(fg.get(p, 0)) != 0))
    M1 = (m2 - Z.deltas[0])*n + 1
    rows = []
    for p in sympy.primerange(5, M1 + 1):
        if p*p <= 3*Z.L: continue
        te = int(fD.get(p, 0)) - int(fg.get(p, 0)); rf = refined_exponent_lz(Z, p); lb = lz_bound_exponent(Z, p)
        rows.append((p, te, lb, rf))
    if rows:
        print('  single-digit primes (p^2 > 3L, p <= M_1): p:true|LZ|refined  ' + '  '.join('%d:%d|%d|%d' % r_ for r_ in rows))
        print('  refined exact at %d/%d ; violations (true > refined): %s ; per n: LZ %.3f refined %.3f true %.3f' % (
            sum(1 for r_ in rows if r_[1] == r_[3]), len(rows), [r_[0] for r_ in rows if r_[1] > r_[3]],
            sum(r_[2]*math.log(r_[0]) for r_ in rows)/n, sum(r_[3]*math.log(r_[0]) for r_ in rows)/n, sum(r_[1]*math.log(r_[0]) for r_ in rows)/n))
    return Z

def asym(s, m1, m2, deltas, verbose=True):
    deltas = sorted(deltas); dmin = deltas[0]; Mp = m2 - 2*dmin
    mp.dps = 30
    # x0: unique positive root of f(X) = 1
    def f(X): return ((2*m1 + m2 + X)/X)**4*(m1 + X)/(m1 + m2 + X)*math.prod([(m1 + d + X)/(m1 + m2 - d + X) for d in deltas])
    lo, hi = mpf('1e-9'), mpf('1e9')
    for _ in range(300):
        mid = (lo*hi)**0.5
        if f(mid) > 1: lo = mid
        else: hi = mid
    x0 = (lo*hi)**0.5
    lg = (2*m1 + m2)*log(108) + sum((m2 - 2*d)*log(m2 - 2*d) for d in deltas) + 4*(2*m1 + m2)*log(2*m1 + m2 + x0) + m1*log(m1 + x0) - (m1 + m2)*log(m1 + m2 + x0)
    for d in deltas: lg += (m1 + d)*log(m1 + d + x0) - (m1 + m2 - d)*log(m1 + m2 - d + x0)
    Dpart = sum(max(m2 - 2*dmin, m2 - d) for d in deltas)
    C2 = Dpart + float(lg)
    # C1 = int_0^1 nu0 dpsi - int_0^{1/M'} nu0 dx/x^2 with nu0(x) = min_y nu(x,y)
    def nu(x, y):
        v = np.floor(2*m1*x + 2*y) - np.floor(m1*x + y) + np.floor(2*(m1 + m2)*x - 2*y) - np.floor((m1 + m2)*x - y)
        v = v + np.floor(3*m1*x + 3*y) + np.floor(3*(m1 + m2)*x - 3*y) - np.floor(y) - np.floor(m2*x - y) - (8*m1 + 3*m2)*math.floor(x)
        for d in deltas: v = v + math.floor((m2 - 2*d)*x) - np.floor(y - d*x) - np.floor((m2 - d)*x - y)
        return v
    def nu0(x):
        b = [0.0, (m2*x) % 1] + [(-m1*x + k/2) % 1 for k in range(2)] + [((m1 + m2)*x - k/2) % 1 for k in range(2)]
        b += [(-m1*x + k/3) % 1 for k in range(3)] + [((m1 + m2)*x - k/3) % 1 for k in range(3)]
        for d in deltas: b += [(d*x) % 1, ((m2 - d)*x) % 1]
        b = np.unique(np.array(b)); b2 = np.append(b, b[0] + 1); mids = (b2[:-1] + b2[1:])/2
        return int(nu(x, mids).min())
    Dg = 3*max(3*(m1 + m2), m2) + 3
    Dg = min(Dg, 2*(m1 + m2) + 3*m2)   # breakpoints in x: multiples of 1/(k*m) with the coefficients above; cap for speed
    pts = sorted(set(Fraction(a, d) for d in range(1, Dg + 1) for a in range(0, d + 1)))
    C1 = mpf(0)
    lo_x = Fraction(1, Mp)
    for a, b in zip(pts[:-1], pts[1:]):
        v = nu0(float((a + b)/2))
        if v:
            fa, fb = mpf(a.numerator)/a.denominator, mpf(b.numerator)/b.denominator
            C1 += v*(mpsi(0, fb) - mpsi(0, fa))
            if a < lo_x:
                bb = min(fb, mpf(lo_x.numerator)/lo_x.denominator)
                if bb > fa and a > 0: C1 -= v*(1/fa - 1/bb)
    net = C2 - float(C1)
    if verbose:
        print('LZ asym: s=%d m1=%d m2=%d deltas=%s : x0 = %s ; log g(x0) = %.6f ; D-part %d ; C2 = %.6f ; C1 (Phi rate) = %.6f ; net C2 - C1 = %+.6f per n' % (
            s, m1, m2, deltas, mp.nstr(x0, 10), float(lg), Dpart, C2, float(C1), net))
    return net, float(C1), C2, float(x0)

def nu0_func(m1, m2, deltas):
    deltas = sorted(deltas)
    def nu(x, y):
        v = np.floor(2*m1*x + 2*y) - np.floor(m1*x + y) + np.floor(2*(m1 + m2)*x - 2*y) - np.floor((m1 + m2)*x - y)
        v = v + np.floor(3*m1*x + 3*y) + np.floor(3*(m1 + m2)*x - 3*y) - np.floor(y) - np.floor(m2*x - y) - (8*m1 + 3*m2)*math.floor(x)
        for d in deltas: v = v + math.floor((m2 - 2*d)*x) - np.floor(y - d*x) - np.floor((m2 - d)*x - y)
        return v
    dx_list = [d*1.0 for d in deltas] + [(m2 - d)*1.0 for d in deltas]
    def nu0(x):
        b = [0.0, (m2*x) % 1] + [(-m1*x + k/2) % 1 for k in range(2)] + [((m1 + m2)*x - k/2) % 1 for k in range(2)]
        b += [(-m1*x + k/3) % 1 for k in range(3)] + [((m1 + m2)*x - k/3) % 1 for k in range(3)]
        b += [(c*x) % 1 for c in dx_list]
        b = np.unique(np.array(b)); b2 = np.append(b, b[0] + 1); mids = (b2[:-1] + b2[1:])/2
        return float(nu(x, mids).min())
    return nu0

FAST = {'du': 0.002, 'nx': 20000}

def fast_C1(m1, m2, deltas, du=None, nx=None):
    du = du or FAST['du']; nx = nx or FAST['nx']
    """C1 = int_{1/M'}^{inf} nu0(x) dx/x^2 = int_{u=1}^{M'} nu0(1/u) du + int_0^1 nu0(x) psi'(1+x) dx  (nu0 1-periodic in x)."""
    deltas = sorted(deltas); Mp = m2 - 2*deltas[0]
    nu0 = nu0_func(m1, m2, deltas)
    us = np.arange(1.0 + du/2, Mp, du)
    part1 = sum(nu0(1.0/u) for u in us)*du
    from mpmath import psi as _psi
    xs = (np.arange(nx) + 0.5)/nx
    w = np.array([float(_psi(1, 1.0 + x)) for x in xs[::200]])   # psi'(1+x) sampled coarsely and interpolated
    wfull = np.interp(xs, xs[::200], w)
    part2 = sum(nu0(float(x))*float(ww) for x, ww in zip(xs, wfull))/nx
    return part1 + part2

def refined_net(s, m1, m2, deltas, verbose=False, exact_C1=False):
    """(s+1) M' + log g(x0) - C1 - (residual top-range cost if the content does not fully cover the theta=1 tails)."""
    deltas = sorted(deltas); dmin = deltas[0]; Mp = m2 - 2*dmin
    mp.dps = 30
    def f(X): return ((2*m1 + m2 + X)/X)**4*(m1 + X)/(m1 + m2 + X)*math.prod([(m1 + d + X)/(m1 + m2 - d + X) for d in deltas])
    lo, hi = mpf('1e-9'), mpf('1e9')
    for _ in range(300):
        mid = (lo*hi)**0.5
        if f(mid) > 1: lo = mid
        else: hi = mid
    x0 = (lo*hi)**0.5
    lg = (2*m1 + m2)*log(108) + sum((m2 - 2*d)*log(m2 - 2*d) for d in deltas) + 4*(2*m1 + m2)*log(2*m1 + m2 + x0) + m1*log(m1 + x0) - (m1 + m2)*log(m1 + m2 + x0)
    for d in deltas: lg += (m1 + d)*log(m1 + d + x0) - (m1 + m2 - d)*log(m1 + m2 - d + x0)
    lg = float(lg)
    C1 = fast_C1(m1, m2, deltas)
    Dpart_LZ = sum(max(Mp, m2 - d) for d in deltas)
    # residual top-range cost: integrate the refined exponent over u in (M', m2 - dmin]
    resid = 0.0; du = 0.01
    for u in np.arange(Mp + du/2, m2 - dmin, du):
        eL, eR, _ = lz_exponents(s, m1, m2, deltas, 1.0/u)
        resid += max(eR, 0.0)*du
    net_LZ = Dpart_LZ + lg - C1
    net_ref = (s + 1)*Mp + resid + lg - C1
    if verbose:
        print('  s=%d m1=%d m2=%d dmin=%d M\'=%d : x0 %.6f, log g %.4f, C1(fast) %.4f ; LZ D-part %d -> refined (s+1)M\' + resid = %d + %.3f ; LZ net %+.4f ; refined net %+.4f' % (
            s, m1, m2, dmin, Mp, float(x0), lg, C1, Dpart_LZ, (s + 1)*Mp, resid, net_LZ, net_ref))
    return net_ref, net_LZ, C1, lg, resid

def valid_params(s, m1, m2, deltas):
    return m1 >= 1 and m2 >= 1 and all(0 <= d and 2*d < m2 for d in deltas) and 2*sum(deltas) < (s - 2)*m2 - 8*m1 and len(deltas) == s + 1

def search_lz(s, m1, m2, deltas, sweeps=3, moves_m=True):
    """coordinate descent on the refined net: single +-1 moves on each delta_j (kept sorted), +-1 on m1, m2, and paired delta moves."""
    cur = (m1, m2, sorted(deltas))
    def val(p):
        if not valid_params(s, *p): return None
        try: return refined_net(s, *p)[0]
        except Exception: return None
    cv = val(cur)
    print('seed s=%d: refined net %+.4f at m1=%d m2=%d deltas=%s' % (s, cv, cur[0], cur[1], cur[2])); sys.stdout.flush()
    for sw in range(sweeps):
        improved = False
        cand_list = []
        for j in range(s + 1):
            for dl in (1, -1):
                ds = cur[2][:]; ds[j] += dl; cand_list.append((cur[0], cur[1], sorted(ds)))
        if moves_m:
            for dm in (1, -1, 2, -2):
                cand_list.append((cur[0] + dm, cur[1], cur[2][:]))
                cand_list.append((cur[0], cur[1] + dm, cur[2][:]))
        for c in cand_list:
            v = val(c)
            if v is not None and v < cv - 1e-6:
                cur, cv, improved = c, v, True
                print('  sweep %d: %+.4f at m1=%d m2=%d deltas=%s' % (sw, cv, cur[0], cur[1], cur[2])); sys.stdout.flush()
        if not improved: break
    print('best s=%d: refined net %+.4f at m1=%d m2=%d deltas=%s' % (s, cv, cur[0], cur[1], cur[2]))
    return cv, cur

def cmd_check_asym(s, m1, m2, deltas, n):
    """compare the asymptotic exponent functions with the finite-n predictor / LZ bound at n for single-digit primes."""
    Z = LZ(s, m1, m2, deltas, n)
    M1 = (m2 - Z.deltas[0])*n + 1
    print('n=%d: p | u=p/n | LZ finite | E_LZ(x) | refined finite | E_ref(x)' % n)
    agree = 0; tot = 0
    for p in sympy.primerange(5, M1 + 1):
        if p*p <= 3*Z.L: continue
        x = n/p
        lb = lz_bound_exponent(Z, p); rf = refined_exponent_lz(Z, p)
        eL, eR, _ = lz_exponents(s, m1, m2, deltas, x)
        tot += 1; agree += (rf == round(eR))
        print('  %5d | %7.3f | %3d | %4g | %3d | %4g %s' % (p, p/n, lb, eL, rf, eR, '' if rf == round(eR) else '<--'))
    print('finite refined = asymptotic refined at %d/%d primes' % (agree, tot))

if __name__ == '__main__':
    cmd = sys.argv[1]; s = int(sys.argv[2])
    if cmd == 'scan':
        FAST['du'] = 0.01; FAST['nx'] = 4000
        a1 = [int(x) for x in sys.argv[3].split(',')]; a2 = [int(x) for x in sys.argv[4].split(',')]
        deltas = [int(x) for x in sys.argv[5].split(',')]
        best = None
        for mm1 in range(a1[0], a1[1] + 1, a1[2]):
            for mm2 in range(a2[0], a2[1] + 1, a2[2]):
                if not valid_params(s, mm1, mm2, deltas): continue
                try: v = refined_net(s, mm1, mm2, deltas)
                except Exception as e: continue
                print('  s=%d m1=%d m2=%d : refined net %+.3f (LZ %+.3f, C1 %.2f, log g %.2f, resid %.2f)' % (s, mm1, mm2, v[0], v[1], v[2], v[3], v[4])); sys.stdout.flush()
                if best is None or v[0] < best[0]: best = (v[0], mm1, mm2)
        print('best: %+.3f at m1=%d m2=%d' % best); sys.exit()
    m1 = int(sys.argv[3]); m2 = int(sys.argv[4]); deltas = [int(x) for x in sys.argv[5].split(',')]
    if cmd == 'asym':
        net, C1, C2, x0 = asym(s, m1, m2, deltas)
        reb, hL, hR = lz_refined(s, m1, m2, deltas)
        print('  REFINED net = %+.4f per n  (LZ net %+.4f, rebate %.4f)' % (net - reb, net, reb))
    elif cmd == 'forms': cmd_forms(s, m1, m2, deltas, int(sys.argv[6]), check=('nocheck' not in sys.argv))
    elif cmd == 'checkasym': cmd_check_asym(s, m1, m2, deltas, int(sys.argv[6]))
    elif cmd == 'refined':
        reb, hL, hR = lz_refined(s, m1, m2, deltas)
    elif cmd == 'refnet':
        refined_net(s, m1, m2, deltas, verbose=True)
    elif cmd == 'search':
        sw = 3
        for a in sys.argv[6:]:
            if a.startswith('sweeps='): sw = int(a[7:])
            if a == 'coarse': FAST['du'] = 0.01; FAST['nx'] = 4000
        search_lz(s, m1, m2, deltas, sweeps=sw)
    elif cmd == 'scan':
        # scan s m1lo,m1hi,step m2lo,m2hi,step deltas  (deltas fixed; coarse grid); prints refined net per (m1, m2)
        FAST['du'] = 0.01; FAST['nx'] = 4000
        a1 = [int(x) for x in sys.argv[3].split(',')]; a2 = [int(x) for x in sys.argv[4].split(',')]
        deltas = [int(x) for x in sys.argv[5].split(',')]
        best = None
        for mm1 in range(a1[0], a1[1] + 1, a1[2]):
            for mm2 in range(a2[0], a2[1] + 1, a2[2]):
                if not valid_params(s, mm1, mm2, deltas): continue
                try: v = refined_net(s, mm1, mm2, deltas)
                except Exception as e: continue
                print('  s=%d m1=%d m2=%d : refined net %+.3f (LZ %+.3f, C1 %.2f, log g %.2f, resid %.2f)' % (s, mm1, mm2, v[0], v[1], v[2], v[3], v[4])); sys.stdout.flush()
                if best is None or v[0] < best[0]: best = (v[0], mm1, mm2)
        print('best: %+.3f at m1=%d m2=%d' % best)
    elif cmd == 'drop':
        # drop the largest-delta blocks two at a time from the given vector, report the refined net at each s
        ds = sorted(deltas); cur_s = s
        while cur_s >= 7:
            refined_net(cur_s, m1, m2, ds, verbose=True); sys.stdout.flush()
            ds = ds[:-2]; cur_s -= 2
