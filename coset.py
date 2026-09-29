"""coset.py — Zudilin-type forms with poles on chosen cosets of (1/b)Z and a periodic rational twist on the summation lattice.

R(t) = (2t+h0) (t+1)_{h0-1} / prod_{j=1..q} prod_{a in A} (t + eta_j n + a/b)_{K_j},   h0 = eta0 n + 1,  K_j = (eta0 - 2 eta_j) n + 1,  n even.
Form:  r = sum_{nu>=0} tw(nu mod T) R(nu),  tw a rational vector of length T.
   sum_{nu>=0} tw(nu) (nu + k + a/b)^{-i} = sum_{r<T} tw(r) T^{-i} [ zeta(i, u/(bT)) - sum_{m<k''} (m + u/(bT))^{-i} ],   (r+k) b + a = k'' bT + u.
So the form is  a0 + sum_i sum_u C_{i,u} zeta(i, u/(bT)),  u in (Z/bT) with u = a mod b  — the "conductor" is bT.
   Catalan (Zudilin 2018):  b=2, A={1}, T=2, tw=(1,-1)      -> conductor 4  (beta values)
   thirds type A:           b=3, A={1,2}, T=1                -> conductor 3  (L(i,chi_-3) and zeta(odd))
   9 = 9 x 1:               b=9, A={1,8}                     -> conductor 9 by placement (one class {+-1} of (Z/9)^*)
   9 = 3 x 3:               b=3, A={1,2}, T=3, tw=[nu = h0 mod 3]  -> conductor 9 by factorization (thirds poles x period-3 sublattice)
Interesting-sum criterion: with R(-t-h0) = eps R(t) (eps = -(-1)^{|A| q}) and tw(-nu-h0) = sigma tw(nu), the half-line sum is a genuine
form iff eps sigma = -1 (else it is half of a residue sum = pi-identity).

usage: python coset.py forms b A eta n [T tw]        e.g. python coset.py forms 9 1,8 3,1,1,1 4
                                                     python coset.py forms 3 1,2 3,1,1,1 8 3 1,0,0
       python coset.py peak b A eta                  normalised Stirling closeness
       python coset.py ledger b A eta n1,n2,... OUT.json [T tw]   per-prime exponents of D/g with p mod bT
"""
import sys, time, math, json
from fractions import Fraction
from flint import fmpq, fmpz
from mpmath import mp, mpf, zeta as mzeta, nsum, inf, log, fabs, loggamma, exp
import sympy

def q2mp(v): return mpf(int(v.p))/mpf(int(v.q))
def fr2q(f): return fmpq(f.numerator, f.denominator)

class Coset:
    def __init__(self, b, A, eta, n, T=1, tw=None):
        self.b, self.A, self.eta, self.n, self.T = b, list(A), list(eta), n, T
        self.tw = [Fraction(1)] if tw is None else [Fraction(x) for x in tw]
        assert len(self.tw) == T
        self.eta0, self.es = eta[0], eta[1:]
        self.h0 = self.eta0*n + 1
        self.K = [(self.eta0 - 2*ej)*n + 1 for ej in self.es]
        assert sorted(self.A) == sorted((b - a) % b for a in self.A), 'coset set must be closed under a -> b-a'
        assert all(0 < a < b for a in self.A)
        deg = self.h0 - len(self.A)*sum(self.K)
        assert deg <= -2, 'degree %d: sum diverges' % deg
        self.deg = deg
    # ---- poles and partial fractions
    def poles(self):
        P = {}
        for ej, K in zip(self.es, self.K):
            for a in self.A:
                for l in range(K):
                    pos = -(ej*self.n + l + Fraction(a, self.b))
                    P[pos] = P.get(pos, 0) + 1
        return P
    def partial_fractions(self):
        h0 = self.h0
        P = self.poles()
        num = [Fraction(m) for m in range(1, h0)]
        den = []
        for ej, K in zip(self.es, self.K):
            for a in self.A:
                for l in range(K): den.append(ej*self.n + l + Fraction(a, self.b))
        coef = {}
        for pos, s0 in P.items():
            p0 = fr2q(pos)
            # numerator zeros at the pole reduce its order (Zudilin's centre pole -h0/2 on the half lattice)
            zeros = (1 if p0 + fmpq(h0, 2) == 0 else 0) + sum(1 for m in num if p0 + fr2q(m) == 0)
            s = s0 - zeros
            if s <= 0: continue
            C = fmpq(2)
            series = [fmpq(0)]*s; series[0] = fmpq(1)
            def mul(ser, alpha, sign):
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
            alph = p0 + fmpq(h0, 2)
            if alph != 0: C *= alph; series = mul(series, alph, +1)
            for m in num:
                alph = p0 + fr2q(m)
                if alph != 0: C *= alph; series = mul(series, alph, +1)
            for d in den:
                alph = p0 + fr2q(d)
                if alph == 0: continue
                C /= alph; series = mul(series, alph, -1)
            for i in range(1, s + 1): coef[(pos, i)] = C*series[s - i]
        self.a = coef; self.P = P
        return coef
    def R_exact(self, t):
        tt = fr2q(t); v = 2*tt + self.h0
        for m in range(1, self.h0): v *= (tt + m)
        for ej, K in zip(self.es, self.K):
            for a in self.A:
                for l in range(K): v /= (tt + ej*self.n + l + fmpq(a, self.b))
        return v
    def identity_ok(self, trials=2):
        import random
        for _ in range(trials):
            t = Fraction(random.randint(1, 10**6), random.randint(1, 10**3)) + Fraction(1, 7)
            lhs = self.R_exact(t); rhs = fmpq(0); tt = fr2q(t)
            for (pos, i), c in self.a.items(): rhs += c/(tt - fr2q(pos))**i
            if lhs != rhs: return False
        return True
    def palindrome_sign(self):
        """eps with R(-t-h0) = eps R(t), checked exactly at a random point."""
        t = Fraction(3, 7)
        v1 = self.R_exact(-t - self.h0); v0 = self.R_exact(t)
        if v1 == v0: return +1
        if v1 == -v0: return -1
        return 0
    # ---- the linear form
    def linear_form(self):
        b, T = self.b, self.T
        N = b*T
        a0 = fmpq(0); C = {}
        for (pos, i), c in self.a.items():
            a = int((-pos) % 1 * b); k = int(-pos - Fraction(a, b)); assert Fraction(k) + Fraction(a, b) == -pos
            for r in range(T):
                w = self.tw[r]
                if w == 0: continue
                tot = (r + k)*b + a
                kk, u = divmod(tot, N)
                cc = c*fr2q(w)/fmpz(T)**i
                C[(i, u)] = C.get((i, u), fmpq(0)) + cc
                ps = fmpq(0); uq = fmpq(u, N)
                for m in range(kk): ps += 1/(m + uq)**i
                a0 -= cc*ps
        C = {key: v for key, v in C.items() if v != 0}
        self.a0, self.C = a0, C
        return a0, C
    def numeric(self, dps=50):
        mp.dps = dps
        h0, n, b = self.h0, self.n, self.b
        def logR(t):
            t = mpf(t)
            v = log(2*t + h0) + loggamma(t + h0) - loggamma(t + 1)
            for ej, K in zip(self.es, self.K):
                for a in self.A:
                    v -= loggamma(t + ej*n + mpf(a)/b + K) - loggamma(t + ej*n + mpf(a)/b)
            return v
        tot = mpf(0)
        for r in range(self.T):
            w = self.tw[r]
            if w == 0: continue
            wq = mpf(w.numerator)/w.denominator
            f = lambda m: exp(logR(self.T*m + r))
            head = sum(f(m) for m in range(0, 2000))
            tail = nsum(f, [2000, inf])
            tot += wq*(head + tail)
        self.r = tot
        return tot
    def form_value(self):
        """a0 + sum C_{i,u} zeta(i,u/N); at i = 1 the symbols are the regularised values -psi(u/N) (sum_u C_{1,u} = 0 kills the pole)."""
        from mpmath import psi
        N = self.b*self.T
        v = q2mp(self.a0)
        s1 = sum((c for (i, u), c in self.C.items() if i == 1), fmpq(0))
        assert s1 == 0, 'level-1 coefficients do not sum to zero: divergent'
        for (i, u), c in self.C.items():
            v += q2mp(c)*(mzeta(i, mpf(u)/N) if i > 1 else -psi(0, mpf(u)/N))
        return v
    def height(self):
        coeffs = [self.a0] + list(self.C.values())
        D = fmpz(1)
        for c in coeffs: D = fmpz(int(sympy.ilcm(int(D), int(c.q))))
        g = fmpz(0)
        for c in coeffs: g = fmpz(int(sympy.igcd(int(g), int((c*D).p))))
        self.D, self.g = D, g
        return D, g
    def logN(self):
        """normalisation making the top coefficients ~integral: b^{(h0-1) - |A| sum K} prod K_j!^{|A|} / (h0-1)!"""
        b = self.b; nA = len(self.A)
        v = -float(loggamma(self.h0)) + (self.h0 - 1)*math.log(b)
        for K in self.K: v += nA*float(loggamma(K + 1)) - nA*K*math.log(b)
        return v

def factor_split(N, ranges):
    f = sympy.factorint(int(N)) if N > 1 else {}
    out = {lab: 0.0 for lab, lo, hi in ranges}
    for p, e in f.items():
        for lab, lo, hi in ranges:
            if lo < p <= hi: out[lab] += e*math.log(p); break
    return out, f

def describe_values(C, b, T):
    """group the Hurwitz symbols by order, show palindromic pairs u <-> N-u and the coefficient ratio."""
    N = b*T
    by_i = {}
    for (i, u), c in C.items(): by_i.setdefault(i, {})[u] = c
    lines = []
    for i in sorted(by_i):
        d = by_i[i]; seen = set(); parts = []
        if i == 1: parts.append('LEVEL 1 (digamma -psi(u/N): logs of cyclotomic units, coefficients sum to 0)')
        for u in sorted(d):
            if u in seen: continue
            v = (N - u) % N
            if v in d and v != u:
                seen |= {u, v}
                ratio = d[v]/d[u]
                parts.append('zeta(%d,%d/%d)%s' % (i, u, N, (' & zeta(%d,%d/%d) ratio %s' % (i, v, N, str(ratio)))))
            else:
                seen.add(u); parts.append('zeta(%d,%d/%d) [unpaired]' % (i, u, N))
        lines.append('    i=%d: ' % i + ' ; '.join(parts))
    return '\n'.join(lines)

def cmd_forms(b, A, eta, n, T=1, tw=None, verbose=True):
    t0 = time.time()
    X = Coset(b, A, eta, n, T, tw)
    X.partial_fractions()
    ok = X.identity_ok(); eps = X.palindrome_sign()
    a0, C = X.linear_form()
    r = X.numeric(); fv = X.form_value()
    D, g = X.height()
    hD = math.log(int(D)); hg = math.log(int(g)) if g > 1 else 0.0
    lN = X.logN()
    M = max(X.K)
    ranges = [('p<=M', 1, M), ('M<p<=2M', M, 2*M), ('2M<p<=3M', 2*M, 3*M), ('3M<p<=6M', 3*M, 6*M), ('p>6M', 6*M, 10**12)]
    sD, fD = factor_split(D, ranges); sg, fg = factor_split(g, ranges)
    net = (hD - hg)/n + float(log(fabs(r)))/n
    if verbose:
        print('coset forms: b=%d A=%s eta=%s n=%d T=%d tw=%s : h0=%d, deg %d, %d poles, max mult %d ; identity %s ; palindrome eps=%+d ; conductor %d  [%.1fs]' % (
            b, A, eta, n, T, [str(x) for x in X.tw], X.h0, X.deg, len(X.P), max(X.P.values()), ok, eps, b*T, time.time() - t0))
        print(describe_values(C, b, T))
        print('  numeric r = %s ; form = %s ; |diff| = %s ; log|r|/n = %.6f (normalised closeness %.4f)' % (mp.nstr(r, 15), mp.nstr(fv, 15), mp.nstr(fabs(r - fv), 3), float(log(fabs(r)))/n, (float(log(fabs(r))) + lN)/n))
        print('  height: log D/n = %.4f, content log g/n = %.4f, primitive log(D/g)/n = %.4f (normalised %.4f) ; NET = %+.4f   (M = K_max = %d)' % (hD/n, hg/n, (hD - hg)/n, (hD - hg - lN)/n, net, M))
        print('  D by range (per n): ' + '  '.join('%s: %.3f' % (lab, sD[lab]/n) for lab, _, _ in ranges) + ' ; content: ' + '  '.join('%s: %.3f' % (lab, sg[lab]/n) for lab, _, _ in ranges if sg[lab]))
        print('  D factorization:', ' '.join('%d^%d' % (p, e) for p, e in sorted(fD.items())))
    return X, net

def cmd_peak(b, A, eta):
    """normalised Stirling closeness: max_x phi(x), phi = lim (log R(xn) + log N_n)/n."""
    mp.dps = 30
    eta0, es = eta[0], eta[1:]; nA = len(A)
    def phi(x):
        x = mpf(x)
        v = (x + eta0)*log(x + eta0) - (x*log(x) if x > 0 else 0) - eta0 + eta0*log(b) - (eta0*log(eta0) - eta0)
        for ej in es:
            kap = eta0 - 2*ej
            v -= nA*((x + ej + kap)*log(x + ej + kap) - (x + ej)*log(x + ej) - kap)
            v += nA*(kap*log(kap) - kap) - nA*kap*log(b)
        return v
    xs = [mpf(k)/2000 for k in range(0, 40000)]
    best = max((phi(x), x) for x in xs)
    lo, hi = max(best[1] - mpf(1)/2000, 0), best[1] + mpf(1)/2000
    for _ in range(100):
        m1 = lo + (hi - lo)/3; m2 = hi - (hi - lo)/3
        if phi(m1) < phi(m2): lo = m1
        else: hi = m2
    x0 = (lo + hi)/2
    print('peak: b=%d A=%s eta=%s : normalised closeness max_x phi = %.6f at x0 = %.5f ; boundary phi(0) = %.6f' % (b, A, eta, phi(x0), x0, phi(0)))
    return float(phi(x0))

def cmd_ledger(b, A, eta, ns, out, T=1, tw=None):
    rows = []
    for n in ns:
        X = Coset(b, A, eta, n, T, tw); X.partial_fractions(); X.linear_form(); D, g = X.height()
        f = sympy.factorint(int(D)); fg = sympy.factorint(int(g)) if g > 1 else {}
        for p in sorted(set(f) | set(fg)):
            rows.append(dict(n=n, p=int(p), e=int(f.get(p, 0)) - int(fg.get(p, 0)), x=n/p, pmod=p % (b*T)))
        print('n=%d: %s' % (n, ' '.join('%d^%d' % (p, e) for p, e in sorted(f.items())))); sys.stdout.flush()
    json.dump(rows, open(out, 'w'))

if __name__ == '__main__':
    cmd = sys.argv[1]; b = int(sys.argv[2]); A = [int(x) for x in sys.argv[3].split(',')]; eta = [int(x) for x in sys.argv[4].split(',')]
    if cmd == 'forms':
        n = int(sys.argv[5]); T = int(sys.argv[6]) if len(sys.argv) > 6 else 1
        tw = [Fraction(x) for x in sys.argv[7].split(',')] if len(sys.argv) > 7 else None
        cmd_forms(b, A, eta, n, T, tw)
    elif cmd == 'peak':
        cmd_peak(b, A, eta)
    elif cmd == 'ledger':
        ns = [int(x) for x in sys.argv[5].split(',')]; out = sys.argv[6]
        T = int(sys.argv[7]) if len(sys.argv) > 7 else 1
        tw = [Fraction(x) for x in sys.argv[8].split(',')] if len(sys.argv) > 8 else None
        cmd_ledger(b, A, eta, ns, out, T, tw)
