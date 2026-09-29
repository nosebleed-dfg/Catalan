"""multibeta.py — the multi-beta route (W. Zudilin, "Arithmetic of Catalan's constant and its relatives", arXiv 1804.09922)
under the ledger: exact linear forms in 1, beta(2), beta(4), ..., beta(s-1), their true denominators prime by prime
against the proven bound, and the asymptotic closeness / height / net per n.

Construction (his (rat2)/(eq:5)): eta = (eta0; eta1..etas) integers, 0 < eta_j < eta0/2, sum eta_j <= (s-1) eta0/2, s odd,
  h0 = eta0 n + 1 (h0-1 even), h_j = eta_j n + 1/2,
  R_n(t) = gamma_n (2t+h0) (t+1)_{h0-1} / prod_j (t+h_j)_{1+h0-2h_j},   gamma_n = 4^{h0-1} prod_{j>=2}(h0-2h_j)! / (h_1-1/2)!^2,
  r_n = sum_{nu>=0} (-1)^nu R_n(nu) = sum_{i even} a_i beta(i) + a_0,
  a_{i,k} = partial-fraction coefficients at t = -(k+1/2), N <= k <= h0-N-1;  a_i = 2^i sum_k (-1)^k a_{i,k};
  a_0 from the partial sums of (-1)^l/(l+1/2)^i (his Lemma 5);  odd a_i vanish (palindrome a_{i,k} = (-1)^i a_{i,h0-1-k}).
Proven bound (Lemma 4/5): Phi_n^{-1} d_M^{s-i} a_i in Z, N = min(h_j-1/2), M = max(h0-2N-1, h_1-1/2),
  Phi_n = prod_{sqrt(2h0) < p <= M} p^{phi_0(n/p)}, phi_0(x) = min_{0<=y<1} phi(x,y), phi = the carry count of Lemma 4.
Asymptotics: closeness = log|r_n|/n -> -decay (Lemma 2/3), height = log den/n -> s M/n - kappa, kappa = int phi_0(x) dx/x^2.
The simple construction of his Section 2 is eta = (3; 1,...,1) (his r_n^{simple} = -r_n here for even n).

usage:
  python multibeta.py forms ETA n [nocheck]   exact forms, verification, per-prime ledger (ETA = 31,10,10,10,10,10,11,11,11,11,12,12,12,12)
  python multibeta.py asym ETA                closed-form closeness / height / net per n, kappa ledger by digit range
  python multibeta.py search s eta0 [eta0max] integer search over eta for odd s (coordinate descent + window enumeration)
"""
import sys, math, time
sys.set_int_max_str_digits(0)
from fractions import Fraction
import flint
from flint import fmpq, fmpz
from mpmath import mp, mpf, zeta as mzeta, log as mlog, psi as mpsi, findroot
from sympy import primerange, factorint

def parse_eta(s): return [int(x) for x in s.split(',')]

# ---------------------------------------------------------------- carry count phi and phi_0 (Lemma 4)
def phi_fn(eta):
    eta0, es = eta[0], eta[1:]; e1 = es[0]
    def phi(x, y):
        fl = math.floor
        v = fl(2*(eta0*x - y)) + fl(2*y) - fl(eta0*x - y) - fl(y) - 2*fl(e1*x) - fl((eta0 - 2*e1)*x)
        for ej in es:
            v += fl((eta0 - 2*ej)*x) - fl(y - ej*x) - fl((eta0 - ej)*x - y)
        return v
    return phi
def _ypoints(eta, x, half):
    """breakpoints of phi(x, .) in y (mod 1) plus the midpoints between consecutive ones: terms like floor(2(eta0 x - y))
    drop just AFTER a breakpoint, so the minimum over y lives in the open intervals."""
    eta0, es = eta[0], eta[1:]
    bps = {0*x, half, (eta0*x) % 1, (eta0*x - half) % 1}
    for ej in es: bps.add((ej*x) % 1); bps.add(((eta0 - ej)*x) % 1)
    bps = sorted(bps); pts = list(bps)
    for u, v in zip(bps, bps[1:] + [bps[0] + 1]): pts.append((u + v)/2)
    return pts
def phi0_fn(eta):
    """phi_0(x) = min over y in [0,1) of phi(x,y) (floats)."""
    phi = phi_fn(eta)
    def phi0(x):
        return min(phi(x, y) for y in _ypoints(eta, x, 0.5))
    return phi0
def phi0_exact(eta, x):
    """phi_0 at a rational x (Fraction), exact floors."""
    eta0, es = eta[0], eta[1:]; e1 = es[0]; fl = math.floor
    x = Fraction(x)
    best = None
    for y in _ypoints(eta, x, Fraction(1, 2)):
        v = fl(2*(eta0*x - y)) + fl(2*y) - fl(eta0*x - y) - fl(y) - 2*fl(e1*x) - fl((eta0 - 2*e1)*x)
        for ej in es: v += fl((eta0 - 2*ej)*x) - fl(y - ej*x) - fl((eta0 - ej)*x - y)
        best = v if best is None else min(best, v)
    return best
def phi_exact(eta, x, y):
    eta0, es = eta[0], eta[1:]; e1 = es[0]; fl = math.floor
    v = fl(2*(eta0*x - y)) + fl(2*y) - fl(eta0*x - y) - fl(y) - 2*fl(e1*x) - fl((eta0 - 2*e1)*x)
    for ej in es: v += fl((eta0 - 2*ej)*x) - fl(y - ej*x) - fl((eta0 - ej)*x - y)
    return v
import numpy as np
def _phi_np(eta, x, y):
    """phi(x, y) for float x and numpy array y (safe away from breakpoints)."""
    eta0, es = eta[0], eta[1:]; e1 = es[0]
    v = np.floor(2*(eta0*x - y)) + np.floor(2*y) - np.floor(eta0*x - y) - np.floor(y) - 2*math.floor(e1*x) - math.floor((eta0 - 2*e1)*x)
    for ej in es: v += math.floor((eta0 - 2*ej)*x) - np.floor(y - ej*x) - np.floor((eta0 - ej)*x - y)
    return v
def _ymids(eta, x):
    """midpoints between consecutive breakpoints (mod 1) of phi(x, .) — the minimum over y lives on open intervals."""
    eta0, es = eta[0], eta[1:]
    b = [0.0, 0.5, (eta0*x) % 1, (eta0*x - 0.5) % 1] + [(ej*x) % 1 for ej in es] + [((eta0 - ej)*x) % 1 for ej in es]
    b = np.unique(np.array(b)); b2 = np.append(b, b[0] + 1)
    return (b2[:-1] + b2[1:])/2
def phi0_float(eta, x):
    return int(_phi_np(eta, x, _ymids(eta, x)).min())
_FAREY = {}
def farey(D):
    if D not in _FAREY:
        pts = set()
        for d in range(1, D + 1):
            for m in range(0, d + 1): pts.add(Fraction(m, d))
        _FAREY[D] = sorted(pts)
    return _FAREY[D]
def phi0_grid(eta, xs):
    """phi_0 at an array of x values (all away from breakpoints), fully vectorised."""
    eta0, es = eta[0], eta[1:]; e1 = es[0]
    xs = np.asarray(xs, dtype=float)
    cols = [np.zeros_like(xs), np.full_like(xs, 0.5), (eta0*xs) % 1, (eta0*xs - 0.5) % 1]
    for ej in es: cols.append((ej*xs) % 1); cols.append(((eta0 - ej)*xs) % 1)
    B = np.sort(np.stack(cols, axis=1), axis=1)                      # (N, m)
    B2 = np.concatenate([B, B[:, :1] + 1], axis=1)
    Y = (B2[:, :-1] + B2[:, 1:])/2                                    # midpoints, (N, m)
    X = xs[:, None]
    v = np.floor(2*(eta0*X - Y)) + np.floor(2*Y) - np.floor(eta0*X - Y) - np.floor(Y) - 2*np.floor(e1*X) - np.floor((eta0 - 2*e1)*X)
    for ej in es: v += np.floor((eta0 - 2*ej)*X) - np.floor(Y - ej*X) - np.floor((eta0 - ej)*X - Y)
    return v.min(axis=1).astype(int)
def phi0_intervals(eta):
    """phi_0 as a step function on [0,1): list of (a, b, value), a,b Fractions (breakpoints m/d, d <= 2 eta0), values at midpoints."""
    pts = farey(2*eta[0])
    mids = np.array([float((a + b)/2) for a, b in zip(pts[:-1], pts[1:])])
    vals = phi0_grid(eta, mids)
    out = []
    for (a, b), v in zip(zip(pts[:-1], pts[1:]), vals):
        v = int(v)
        if out and out[-1][2] == v and out[-1][1] == a: out[-1] = (out[-1][0], b, v)
        else: out.append((a, b, v))
    return out
def kappa_from_intervals(iv, mu):
    """kappa = int_{1/mu}^inf phi_0(x) dx/x^2 (phi_0 1-periodic), mu = M/n. Returns kappa and a ledger by digit range."""
    mp.dps = 30
    kap = mpf(0); ledger = []
    lo = Fraction(1, int(mu)) if float(mu) == int(mu) else Fraction(mu).limit_denominator(10**6)
    lo = 1/Fraction(mu)
    for a, b, v in iv:
        if v == 0: continue
        # m >= 1 periods
        c1 = v*(mpsi(0, mpf(1) + mpf(b.numerator)/b.denominator) - mpsi(0, mpf(1) + mpf(a.numerator)/a.denominator))
        # m = 0 part on [max(a,1/mu), b)
        c0 = mpf(0)
        if b > lo:
            aa = max(a, lo); c0 = v*(mpf(1)/mpf(aa.numerator)*aa.denominator - mpf(1)/mpf(b.numerator)*b.denominator)
        kap += c0 + c1; ledger.append((a, b, v, float(c0), float(c1)))
    return kap, ledger

# ---------------------------------------------------------------- refined bound (this work)
# At a pole k of multiplicity s_k with top coefficient C_k (= a_{s_k,k}), the lower Laurent coefficients are Taylor
# coefficients of C_k prod(1+u/alpha)/prod(1+u/beta); a derivative can cost a factor p only through a factor divisible by p,
# and if no DENOMINATOR factor at that pole is divisible by p the p-singular part is a polynomial of degree D_k = #{p-divisible
# numerator factors}, so v_p(a_{i,k}) >= v_p(C_k) - min(s_k - i, D_k).  Partial sums of order i carry p^{-i} only at poles
# with |k - centre| >= (p+1)/2 (the "PS window").  With v_p(C_k) = phi(n/p,k/p) - (s - s_k) (Legendre, p > sqrt(2h0)):
#   v_p(a_0) >= -s + min( min_{k in PS} phi(x,y_k),  min_k [phi(x,y_k) + max(1, s_k - D_k)] ),
#   v_p(a_i) >= -s + min_k [phi(x,y_k) + max(i, s_k - D_k)]   (i >= 2).
def kref_point(eta, x):
    """refined rebate exponent at digit x = n/p (float, away from breakpoints): s - (a_0 exponent bound).
    Equals phi_0(x) once the PS window covers a period and every pole has a p-divisible denominator factor (x >= 2/mu')."""
    x = float(x); eta0, es = eta[0], eta[1:]; emin = min(es)
    ylo, yhi = emin*x, (eta0 - emin)*x
    cands = [ylo, yhi, eta0*x/2 - 0.5, eta0*x/2 + 0.5]
    def addm(base, step):   # base + m*step for all integers m with value in [ylo, yhi]
        a, b = (ylo - base)/step, (yhi - base)/step
        for m in range(math.ceil(min(a, b)), math.floor(max(a, b)) + 1): cands.append(base + m*step)
    addm(eta0*x, -0.5); addm(0.0, 0.5)
    for ej in es: addm(ej*x, 1.0); addm((eta0 - ej)*x, -1.0)
    cs = np.unique(np.array([c for c in cands if ylo <= c <= yhi]))
    ys = (cs[:-1] + cs[1:])/2 if len(cs) > 1 else cs
    ph = _phi_np(eta, x, ys)
    best_ps = None; best_all = None
    for y, p_ in zip(ys, ph):
        sy = sum(1 for ej in es if ej*x <= y <= (eta0 - ej)*x)
        D = None
        for ej in es:
            lo, hi = ej*x - y, (eta0 - ej)*x - y
            m = math.ceil(lo)
            if m == 0: m = 1
            if m <= hi: D = 10**9; break
        if D is None: D = math.floor(eta0*x - y + 0.5) + math.floor(y + 0.5)
        v_all = p_ + max(1, sy - D)
        best_all = v_all if best_all is None else min(best_all, v_all)
        if y <= eta0*x/2 - 0.5 or y >= eta0*x/2 + 0.5:
            best_ps = p_ if best_ps is None else min(best_ps, p_)
    v = best_all if best_ps is None else min(best_ps, best_all)
    return int(v)
def kappa_refined(eta, mu, verbose=False):
    """kappa_ref = int_{1/mu}^inf kref(x) dx/x^2: kref on the Farey grid over [1/mu, X], X = 2/mu' (beyond it kref = phi_0),
    then phi_0 on [X, ceil(X)) and the digamma tail from ceil(X)."""
    eta0, es = eta[0], eta[1:]; emin = min(es); mup = eta0 - 2*emin
    X = Fraction(2, mup); lo = Fraction(1, mu)
    if X <= lo: X = lo
    XI = math.ceil(X)
    base = farey(2*eta0)
    pts = sorted(set([lo, X] + [q + m for m in range(0, XI) for q in base if lo <= q + m <= X]))
    mp.dps = 30
    kap = mpf(0); ledger = []
    for a, b in zip(pts[:-1], pts[1:]):
        v = kref_point(eta, float((a + b)/2))
        if v: kap += v*(mpf(1)/mpf(a.numerator)*a.denominator - mpf(1)/mpf(b.numerator)*b.denominator)
        ledger.append((a, b, v))
    iv = phi0_intervals(eta)
    fr = lambda q: mpf(q.numerator)/q.denominator
    for a, b, v in iv:
        if not v: continue
        # [X, XI): the period containing X
        if XI > X:
            aa, bb = max(a + XI - 1, X), b + XI - 1
            if bb > aa: kap += v*(1/fr(aa) - 1/fr(bb))
        # [XI, inf): digamma tail
        kap += v*(mpsi(0, mpf(XI) + fr(b)) - mpsi(0, mpf(XI) + fr(a)))
    if verbose:
        comp = []
        for a, b, v in ledger:
            if comp and comp[-1][2] == v and comp[-1][1] == a: comp[-1] = (comp[-1][0], b, v)
            else: comp.append((a, b, v))
        print('  refined rebate steps on [1/mu, 2/mu\') = [%s, %s): ' % (lo, X) + ' '.join('[%s,%s)=%d' % (a, b, v) for a, b, v in comp if v))
    return kap

# ---------------------------------------------------------------- asymptotics (Lemma 2/3)
def decay(eta):
    """-lim log|r_n|/n (positive) with the saddle x_0 of Lemma 3; returns (decay, x0, xs)."""
    eta0, es = eta[0], eta[1:]; s = len(es); e1 = es[0]
    mp.dps = 30
    def P(x):
        return x*math.prod([(eta0 - ej) - ej*x for ej in es]) - math.prod([ej - (eta0 - ej)*x for ej in es])
    # unique zero in (0,1) (Lemma 3 hypothesis): count sign changes on a fine grid, then bisect
    grid = [mpf(i)/20000 for i in range(1, 20000)]
    signs = [1 if P(g) > 0 else -1 for g in grid]
    changes = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    if changes != 1: print('  WARNING: Lemma 3 polynomial has %d sign changes in (0,1); decay not reliable' % changes)
    lo, hi = mpf('1e-12'), mpf(1) - mpf('1e-12')
    assert P(lo) < 0 < P(hi), 'no sign change for Lemma 3 polynomial'
    for _ in range(200):
        mid = (lo + hi)/2
        if P(mid) < 0: lo = mid
        else: hi = mid
    x0 = (lo + hi)/2
    xs = [(ej - (eta0 - ej)*x0)/((eta0 - ej) - ej*x0) for ej in es]
    logmax = sum(ej*mlog(xj) + (eta0 - 2*ej)*mlog(1 - xj) for ej, xj in zip(es, xs)) - eta0*mlog(1 + math.prod(xs))
    pref = eta0*mlog(4*eta0) - 2*e1*mlog(e1) - (eta0 - 2*e1)*mlog(eta0 - 2*e1)
    return -(pref + logmax), x0, xs
def asym(eta, verbose=True, refined=True):
    eta0, es = eta[0], eta[1:]; s = len(es); e1 = es[0]
    emin = min(es); mu = max(eta0 - 2*emin, e1)
    d, x0, xs = decay(eta)
    iv = phi0_intervals(eta)
    kap, ledger = kappa_from_intervals(iv, mu)
    naive = s*mu; height = naive - kap; net = height - d
    kr = kappa_refined(eta, mu, verbose=verbose) if refined else kap
    height_r = naive - kr; net_r = height_r - d
    if verbose:
        print('eta = (%d; %s)  s = %d  M/n = %d' % (eta0, ','.join(map(str, es)), s, mu))
        print('  closeness = %.10f   height = %d (naive s M/n) - %.10f (carry rebate kappa, Zudilin) = %.10f   net = %+.6f per n  (%+.5f per unit eta0)'
              % (-d, naive, kap, height, net, net/eta0))
        if refined:
            print('  refined:  height = %d - %.10f (kappa_ref) = %.10f   net = %+.6f per n  (%+.5f per unit eta0)   extra rebate %.6f' % (naive, kr, height_r, net_r, net_r/eta0, kr - kap))
        print('  saddle x0 = %s ; x_j = %s' % (mp.nstr(x0, 10), ' '.join(mp.nstr(v, 6) for v in xs)))
        by_val = {}
        for a, b, v, c0, c1 in ledger: by_val[v] = by_val.get(v, 0) + c0 + c1
        print('  kappa by phi_0 value: ' + '  '.join('%d: %.4f' % (v, by_val[v]) for v in sorted(by_val)))
        print('  phi_0 steps: ' + ' '.join('[%s,%s)=%d' % (a, b, v) for a, b, v in iv if v))
    return dict(decay=float(d), kappa=float(kap), height=float(height), net=float(net), mu=mu, intervals=iv,
                kappa_ref=float(kr), height_ref=float(height_r), net_ref=float(net_r))

# ---------------------------------------------------------------- exact forms
def prod_range(a, b):
    """product of the integers a..b (inclusive) skipping zeros; returns (product, zeros)."""
    z = 0; P = 1
    if a <= 0 <= b: z = 1
    if b >= 1: P *= math.prod(range(max(a, 1), b + 1))
    if a <= -1: P *= math.prod(range(a, min(b, -1) + 1))
    return P, z
def forms(eta, n, want_coeffs=False):
    eta0, es = eta[0], eta[1:]; s = len(es)
    h0 = eta0*n + 1; assert (h0 - 1) % 2 == 0, 'need eta0 n even'
    hs = [ej*n for ej in es]                    # h_j - 1/2
    L = [h0 - 2*h for h in hs]                  # pole count of factor j
    N = min(hs); M = max(h0 - 2*N - 1, hs[0])
    gamma = fmpq(4)**(h0 - 1)*math.prod(math.factorial(h0 - 2*h - 1) for h in hs[1:])/fmpq(math.factorial(hs[0]))**2
    # prefix sums of reciprocal powers: integers 1..(h0) and odd numbers 1..(2h0)
    smax = s - 1
    IP = [None] + [[fmpq(0)]*(h0 + 2) for _ in range(smax)]
    OP = [None] + [[fmpq(0)]*(2*h0 + 2) for _ in range(smax)]
    for m in range(1, smax + 1):
        acc = fmpq(0)
        for j in range(1, h0 + 2): acc += fmpq(1, j)**m; IP[m][j] = acc
        acc = fmpq(0)
        for o in range(1, 2*h0 + 2):
            if o % 2: acc += fmpq(1, o)**m
            OP[m][o] = acc
    def int_sum(m, a, b):        # sum_{d=a}^{b}, d != 0, of d^{-m}
        tot = fmpq(0)
        if b >= 1: lo = max(a, 1); tot += IP[m][b] - IP[m][lo - 1]
        if a <= -1: hi = min(b, -1); tot += (-1)**m*(IP[m][-a] - IP[m][-hi - 1])
        return tot
    def odd_sum(m, a, b):        # sum over odd o in [a,b] of o^{-m}
        tot = fmpq(0)
        if b >= 1: lo = max(a, 1); tot += OP[m][b] - OP[m][lo - 1]
        if a <= -1: hi = min(b, -1); tot += (-1)**m*(OP[m][-a] - OP[m][-hi - 1])
        return tot
    center = (h0 - 1)//2
    A = {}   # (i,k) -> a_{i,k}
    mult = {}
    for k in range(N, h0 - N):
        sk = sum(1 for h in hs if h <= k <= h0 - h - 1)
        if sk == 0: continue
        z = 1 if k == center else 0
        ordr = sk - z                        # effective pole order
        # constant: gamma * cf * prod_o (o/2) / prod_d d
        cf = 2 if z else (h0 - 2*k - 1)
        num_odd = math.prod(range(1 - 2*k, 2*h0 - 2 - 2*k, 2))
        den_int = 1
        for h, Lj in zip(hs, L):
            P, zz = prod_range(h - k, h + Lj - 1 - k); den_int *= P
        C = gamma*cf*fmpq(num_odd, 2**(h0 - 1))/den_int
        if ordr == 0:
            continue
        # log-series coefficients
        Lc = [None]
        for m in range(1, ordr):
            S = fmpq(2)**m*odd_sum(m, 1 - 2*k, 2*h0 - 3 - 2*k)
            if not z: S += fmpq(2, h0 - 2*k - 1)**m
            for h, Lj in zip(hs, L): S -= int_sum(m, h - k, h + Lj - 1 - k)
            Lc.append((-1)**(m + 1)*S/m)
        b = [C]
        for m in range(1, ordr):
            acc = fmpq(0)
            for r in range(1, m + 1): acc += r*Lc[r]*b[m - r]
            b.append(acc/m)
        for i in range(1, ordr + 1):
            A[(i, k)] = b[ordr - i]
        mult[k] = (sk, z)
    # a_i and a_0
    a = {i: fmpq(0) for i in range(1, s + 1)}
    for (i, k), c in A.items(): a[i] += c*(-1)**k
    for i in a: a[i] *= fmpq(2)**i
    PS = {i: [fmpq(0)] for i in range(1, s + 1)}; NS = {i: [fmpq(0)] for i in range(1, s + 1)}
    for i in range(1, s + 1):
        for l in range(0, center + 2): PS[i].append(PS[i][-1] + fmpq((-1)**l*2**i, (2*l + 1)**i))
        for m in range(1, center + 2): NS[i].append(NS[i][-1] + fmpq((-1)**(m + i)*2**i, (2*m - 1)**i))
    a0 = fmpq(0)
    for (i, k), c in A.items():
        if k < center: a0 += (-1)**k*c*NS[i][center - k]
        elif k > center: a0 -= (-1)**k*c*PS[i][k - center]
    info = dict(h0=h0, hs=hs, L=L, N=N, M=M, s=s, gamma=gamma, mult=mult, center=center)
    if want_coeffs: info['A'] = A
    return a, a0, info

def R_exact(eta, n, info, t):
    """R_n(t) exactly at a rational t (fmpq)."""
    h0, hs, L, gamma = info['h0'], info['hs'], info['L'], info['gamma']
    v = gamma*(2*t + h0)
    for i in range(1, h0): v *= (t + i)
    for h, Lj in zip(hs, L):
        for i in range(Lj): v /= (t + fmpq(2*h + 1, 2) + i)
    return v
def identity_check(eta, n, info, A, trials=3):
    """exact check R(t) = sum a_{i,k}/(t+k+1/2)^i at random rational t."""
    import random
    random.seed(7)
    for _ in range(trials):
        t = fmpq(random.randint(1, 10**6), random.randint(1, 10**6)) + random.randint(0, 3)
        lhs = R_exact(eta, n, info, t)
        rhs = fmpq(0)
        for (i, k), c in A.items(): rhs += c/(t + fmpq(2*k + 1, 2))**i
        if lhs != rhs: return False
    return True
def numeric_r(eta, n, info, dps, maxterms=200000):
    """r_n = sum_{nu>=0} (-1)^nu R_n(nu): direct summation by the term ratio when the terms decay fast enough,
    else mpmath nsum (Richardson/Shanks on the alternating series). Returns (value, largest term, method)."""
    mp.dps = dps
    h0, hs, L, gamma = info['h0'], info['hs'], info['L'], info['gamma']
    R0 = gamma*h0*math.factorial(h0 - 1)
    for h, Lj in zip(hs, L):
        R0 /= fmpq(math.prod(range(2*h + 1, 2*h + 2*Lj, 2)), 2**Lj)
    R0m = mpf(int(R0.p))/mpf(int(R0.q))
    powerdecay = sum(L) - h0 - 1
    def ratio(nu):
        r = mpf(2*nu + 2 + h0)/(2*nu + h0)*mpf(nu + h0)/(nu + 1)
        for h in hs: r *= mpf(2*nu + 2*h + 1)/(2*nu + 2 + 2*h0 - 2*h - 1)
        return r
    term = R0m; tot = term; mx = abs(term); nu = 0
    while nu < maxterms:
        term = -term*ratio(nu); nu += 1; tot += term; mx = max(mx, abs(term))
        if abs(term) < mx*mpf(10)**(-dps - 5) and nu > h0: return tot, mx, 'direct(%d terms)' % nu
    # slow polynomial decay: accelerate. terms via the ratio recurrence cached up to nu, then nsum on an mpf function
    cache = [R0m]
    def Rm(nu):
        nu = int(nu)
        while len(cache) <= nu: cache.append(-cache[-1]*ratio(len(cache) - 1))
        return cache[nu]
    val = mp.nsum(lambda nu: Rm(nu), [0, mp.inf])
    return val, mx, 'nsum (power decay nu^-%d)' % powerdecay
def beta_val(i):
    return (mzeta(i, mpf(1)/4) - mzeta(i, mpf(3)/4))/mpf(4)**i

def factor_den(q, bound):
    f = {}
    for p in primerange(2, bound + 1):
        if q % p == 0:
            e = 0
            while q % p == 0: q //= p; e += 1
            f[p] = e
    if q != 1:
        for p, e in factorint(q).items(): f[int(p)] = f.get(int(p), 0) + int(e)
    return f
def vp(p, n):
    e = 0
    while n % p == 0: n //= p; e += 1
    return e
def run_forms(eta, n, check=True):
    t0 = time.time()
    a, a0, info = forms(eta, n, want_coeffs=True)
    s, h0, M, N = info['s'], info['h0'], info['M'], info['N']
    print('eta = (%d; %s)  n = %d : h0 = %d, poles k = %d..%d (%d), M = %d, s = %d   [%.1fs]' % (eta[0], ','.join(map(str, eta[1:])), n, h0, N, h0 - N - 1, h0 - 2*N, M, s, time.time() - t0))
    odd = [i for i in range(1, s + 1, 2) if a[i] != 0]
    print('  odd a_i vanish: %s ; partial fractions reproduce R(t) exactly at random t: %s' % ('yes' if not odd else 'NO %s' % odd, identity_check(eta, n, info, info['A'])))
    sys.stdout.flush()
    if check:
        dps = 60 + int(0.5*n*eta[0]*s)  # generous
        rn, mx, meth = numeric_r(eta, n, info, dps)
        form = mpf(int(a0.p))/mpf(int(a0.q))
        for i in range(2, s + 1, 2): form += mpf(int(a[i].p))/mpf(int(a[i].q))*beta_val(i)
        print('  numeric r_n = %s [%s]; form = %s ; |diff|/|r_n| = %s ; largest term %s' % (mp.nstr(rn, 15), meth, mp.nstr(form, 15), mp.nstr(abs(rn - form)/abs(rn), 3), mp.nstr(mx, 5)))
        print('  closeness log|r_n|/n = %.6f' % (float(mlog(abs(rn)))/n))
        sys.stdout.flush()
    # denominators
    dens = {0: int(a0.q)}
    for i in range(2, s + 1, 2): dens[i] = int(a[i].q)
    fac = {i: factor_den(q, 2*M + 100) for i, q in dens.items()}
    primes = sorted(set().union(*[set(f) for f in fac.values()]))
    common = {p: max(f.get(p, 0) for f in fac.values()) for p in primes}
    logden = sum(e*math.log(p) for p, e in common.items())
    print('  height log den/n = %.6f   (a_0: %.6f, a_2: %.6f)' % (logden/n, sum(e*math.log(p) for p, e in fac[0].items())/n, sum(e*math.log(p) for p, e in fac[2].items())/n))
    # bound
    dM = {p: int(math.log(M)/math.log(p) + 1e-9) for p in primes}
    lo = math.sqrt(2*h0)
    print('  per prime: p | n/p | phi_0 | bound a_0 (s vp(dM) - phi_0) | actual a_0 | actual common | max_i excess | bound a_2 | actual a_2')
    tot_b = 0.0; tot_a = 0.0; gain = {}
    for p in primes:
        x = Fraction(n, p); ph = phi0_exact(eta, x % 1) if (lo < p <= M) else 0
        b0 = s*dM[p] - ph; b2 = (s - 2)*dM[p] - ph
        e0 = fac[0].get(p, 0); ec = common[p]; e2 = fac[2].get(p, 0)
        tot_b += max(b0, 0)*math.log(p); tot_a += ec*math.log(p)
        flag = '' if ec <= b0 else '  <-- EXCEEDS BOUND'
        print('   %4d | %6s | %2d | %3d | %3d | %3d | %+3d | %3d | %3d%s' % (p, x, ph, b0, e0, ec, b0 - ec, b2, e2, flag))
        key = 'p<=sqrt(2h0)' if p <= lo else ('sqrt<p<=M/2' if p <= M/2 else ('M/2<p<=M' if p <= M else 'p>M'))
        gain[key] = gain.get(key, 0.0) + (b0 - ec)*math.log(p)
    print('  bound height (a_0) per n = %.6f ; actual common per n = %.6f ; rebate below the bound per n = %.6f' % (tot_b/n, tot_a/n, (tot_b - tot_a)/n))
    print('  rebate by prime range per n: ' + '  '.join('%s: %.4f' % (k, v/n) for k, v in gain.items()))
    # refined predictor (single-digit primes): per pole k: s_k, v_p(C_k) exact from the top coefficient, D_k(p), PS window
    A = info['A']; hs, L, center = info['hs'], info['L'], info['center']
    print('  refined bound test (sqrt(2h0) < p <= M): p | phi_0 | Zudilin bound | refined bound a_0 | actual a_0 | refined a_2 | actual a_2 | (PS term, non-PS term) | kref(x) asymptotic')
    tot_ref = 0.0
    for p in primes:
        if not (lo < p <= M): continue
        x = Fraction(n, p); ph = phi0_exact(eta, x % 1)
        ps_term = None; np_term = {i: None for i in range(1, s + 1)}
        for k, (sk, z) in info['mult'].items():
            ordr = sk - z
            C = A[(ordr, k)]; vC = vp(p, int(C.p)) - vp(p, int(C.q))
            # p-divisible analytic factors at this pole
            odd_lo, odd_hi = 1 - 2*k, 2*h0 - 3 - 2*k
            Dnum = sum(1 for o in range(odd_lo, odd_hi + 1, 2) if o % p == 0)
            if not z and (h0 - 2*k - 1) % p == 0: Dnum += 1
            Dden = 0
            for h, Lj in zip(hs, L):
                for dd in range(h - k, h + Lj - k):
                    if dd != 0 and dd % p == 0: Dden += 1
            D = 10**9 if Dden else Dnum
            psact = (k <= (h0 - 2 - p)//2) or (k >= (h0 + p)/2)
            if psact:
                v = ordr - vC
                ps_term = v if ps_term is None else max(ps_term, v)
            for i in range(1, ordr + 1):
                v = min(ordr - i, D) - vC
                np_term[i] = v if np_term[i] is None else max(np_term[i], v)
        ref0 = max([v for v in (ps_term, np_term[1]) if v is not None])
        ref2 = np_term[2] if np_term[2] is not None else -99
        e0 = fac[0].get(p, 0); e2 = fac[2].get(p, 0); zb = s*dM[p] - ph
        tot_ref += max(ref0, 0)*math.log(p)
        flag = '' if e0 <= ref0 else '  <-- ACTUAL EXCEEDS REFINED BOUND'
        kl, kr_ = s - kref_point(eta, float(x)*(1 - 1e-7)), s - kref_point(eta, float(x)*(1 + 1e-7))
        print('   %4d | %2d | %3d | %3d | %3d | %3d | %3d | (%s, %s) | %s%s' % (p, ph, zb, ref0, e0, ref2, e2, ps_term, np_term[1], '%d' % kl if kl == kr_ else '%d|%d' % (kl, kr_), flag))
    print('  refined height (single-digit primes only) per n = %.6f ; Zudilin bound on the same primes = %.6f ; actual on the same primes = %.6f' % (
        tot_ref/n, sum(max(s*dM[p] - (phi0_exact(eta, Fraction(n, p) % 1)), 0)*math.log(p) for p in primes if lo < p <= M)/n,
        sum(common[p]*math.log(p) for p in primes if lo < p <= M)/n))
    return a, a0, info, fac, common

# ---------------------------------------------------------------- search
def valid(eta):
    eta0, es = eta[0], eta[1:]; s = len(es)
    return all(0 < e < eta0/2 for e in es) and 2*sum(es) <= (s - 1)*eta0
def search(s, eta0, eta0max=None, key='net_ref'):
    """coordinate descent over integer eta (eta_1 = the smallest, sorted) minimising key ('net' Zudilin bound, 'net_ref' refined)."""
    eta0max = eta0max or eta0
    results = []; cache = {}
    def val(e):
        t = tuple(e)
        if t not in cache: cache[t] = asym(e, verbose=False, refined=(key == 'net_ref'))
        return cache[t][key]
    for e0 in range(eta0, eta0max + 1):
        base = max(1, round(e0/3))
        best = None
        starts = [[base]*s, [max(1, base - 1)]*s, [min(base + 1, (e0 - 1)//2)]*s]
        for st in starts:
            cur = [e0] + sorted(st)
            if not valid(cur): continue
            cv = val(cur)
            improved = True
            while improved:
                improved = False
                moves = []
                for j in range(1, s + 1):
                    for dlt in (1, -1):
                        cand = cur[:]; cand[j] += dlt
                        cand = [cand[0]] + sorted(cand[1:])
                        if valid(cand): moves.append(cand)
                # also paired moves (one up, one down) to walk along the sum constraint
                for j in range(1, s + 1):
                    for l in range(1, s + 1):
                        if j == l: continue
                        cand = cur[:]; cand[j] += 1; cand[l] -= 1
                        cand = [cand[0]] + sorted(cand[1:])
                        if valid(cand): moves.append(cand)
                for cand in moves:
                    v = val(cand)
                    if v < cv - 1e-9: cur, cv, improved = cand, v, True
            if best is None or cv < best[0]: best = (cv, cur)
        r = cache[tuple(best[1])]
        results.append((best[0], best[1], r))
        print('eta0 = %2d : best %s %+.4f per n (%+.5f per unit eta0)  eta = %s  [closeness %.4f; Zudilin height %.4f = %d - %.4f, net %+.4f; refined height %.4f, net %+.4f]' % (
            e0, key, best[0], best[0]/e0, best[1], -r['decay'], r['height'], s*r['mu'], r['kappa'], r['net'], r['height_ref'], r['net_ref']))
        sys.stdout.flush()
    results.sort(key=lambda t: t[0]/t[1][0])
    print('\nbest per unit eta0: %s %+.5f, eta = %s' % (key, results[0][0]/results[0][1][0], results[0][1]))
    return results

def mech2(eta, n, p):
    """anatomy of the p-singular part of a_0: which (i,k) terms attain the minimal p-order and how their leading residues sum."""
    a, a0, info = forms(eta, n, want_coeffs=True)
    A = info['A']; s = info['s']; h0 = info['h0']; center = info['center']
    PS = {i: [fmpq(0)] for i in range(1, s + 1)}
    for i in range(1, s + 1):
        for l in range(0, center + 2): PS[i].append(PS[i][-1] + fmpq((-1)**l*2**i, (2*l + 1)**i))
    def vpq(q):  # p-adic order of an fmpq
        if q == 0: return 10**9
        return vp(p, int(q.p)) - vp(p, int(q.q))
    e0 = -vpq(a0)
    print('eta = (%d; %s) n = %d p = %d : a_0 has p-order %d (exponent %d in the denominator)' % (eta[0], ','.join(map(str, eta[1:])), n, p, -e0, e0))
    # a_0 = 2 sum_i (-1)^{i+1} sum_{k<center} (-1)^k a_{i,k} PS_i(center-k)
    terms = []
    for (i, k), c in A.items():
        if k >= center: continue
        t = 2*(-1)**(i + 1)*(-1)**k*c*PS[i][center - k]
        terms.append((vpq(t), i, k, t))
    terms.sort()
    vmin = terms[0][0]
    lead = [t for t in terms if t[0] == vmin]
    print('  minimal term order %d (would give exponent %d); %d terms attain it:' % (vmin, -vmin, len(lead)))
    tot = fmpq(0); byi = {}
    for v, i, k, t in lead:
        r = (t*fmpq(p)**(-vmin))
        res = (int(r.p)*pow(int(r.q), -1, p)) % p
        tot += t; byi.setdefault(i, []).append((k, res))
    for i in sorted(byi):
        ks = byi[i]
        print('   i = %2d : k = %s' % (i, ' '.join('%d' % k for k, _ in ks)))
        print('            residues mod p of the leading coefficients: %s   sum mod p = %d' % (' '.join('%d' % r for _, r in ks), sum(r for _, r in ks) % p))
    print('  order of the sum of the leading terms: %d ; order of the full a_0: %d ; next term orders: %s' % (vpq(tot), vpq(a0), sorted(set(v for v, *_ in terms))[:4]))
    # also the a_{i,k} alone over the PS window for the top i
    for i in sorted(byi):
        S = fmpq(0)
        for k, _ in byi[i]: S += (-1)**k*A[(i, k)]
        print('  sum_{k in leading set} (-1)^k a_{%d,k}: p-order %d' % (i, vpq(S)))

def refined_pred(eta, n, info, A, p):
    """finite-n refined bound for the a_0 exponent at prime p (single-digit primes), from exact v_p(C_k), D_k(p), s_k, window."""
    s, h0, hs, L = info['s'], info['h0'], info['hs'], info['L']
    ps_term = None; np_term = None
    def odd_mult(a, b):   # odd multiples of p in [a, b] (a, b integers)
        # q odd with a/p <= q <= b/p:  #{odd q <= B} = floor((floor(B)+1)/2)
        hi = (b // p + 1)//2; lo = ((a - 1)//p + 1)//2
        return hi - lo
    def has_nonzero_mult(a, b):   # a nonzero multiple of p in [a, b]
        if b < a: return False
        cnt = b//p - (a - 1)//p
        if a <= 0 <= b: cnt -= 1
        return cnt > 0
    for k, (sk, z) in info['mult'].items():
        ordr = sk - z
        C = A[(ordr, k)]; vC = vp(p, int(C.p)) - vp(p, int(C.q))
        Dnum = odd_mult(1 - 2*k, 2*h0 - 3 - 2*k)
        if not z and (h0 - 2*k - 1) % p == 0: Dnum += 1
        Dden = any(has_nonzero_mult(h - k, h + Lj - 1 - k) for h, Lj in zip(hs, L))
        D = 10**9 if Dden else Dnum
        if (k <= (h0 - 2 - p)//2) or (k >= (h0 + p)/2):
            v = ordr - vC; ps_term = v if ps_term is None else max(ps_term, v)
        v = min(ordr - 1, D) - vC; np_term = v if np_term is None else max(np_term, v)
    return max(v for v in (ps_term, np_term) if v is not None)
def event_type(eta, n, info, A, a0, p):
    """classify the extra cancellation at p: 'none' (a_0 exponent = order of the leading terms), 'per-order', 'cross', 'other'."""
    s, center = info['s'], info['center']
    PS = {i: [fmpq(0)] for i in range(1, s + 1)}
    for i in range(1, s + 1):
        for l in range(0, center + 2): PS[i].append(PS[i][-1] + fmpq((-1)**l*2**i, (2*l + 1)**i))
    def vpq(q): return 10**9 if q == 0 else vp(p, int(q.p)) - vp(p, int(q.q))
    terms = []
    for (i, k), c in A.items():
        if k >= center: continue
        t = 2*(-1)**(i + 1)*(-1)**k*c*PS[i][center - k]
        if t != 0: terms.append((vpq(t), i, k, t))
    vmin = min(v for v, *_ in terms)
    if -vpq(a0) == -vmin: return 'none', vmin, 0
    lead = [t for t in terms if t[0] == vmin]
    byi = {}
    for v, i, k, t in lead: byi.setdefault(i, fmpq(0)); byi[i] += t
    per_order = all(vpq(S) > vmin for S in byi.values())
    return ('per-order' if per_order else 'cross'), vmin, len(lead)
def sweep(eta, n1, n2, step, out=None):
    eta0 = eta[0]; rows = []
    for n in range(n1, n2 + 1, step):
        if (eta0*n) % 2: continue
        a, a0, info = forms(eta, n, want_coeffs=True)
        A = info['A']; s, h0, M = info['s'], info['h0'], info['M']
        q = int(a0.q); f0 = factor_den(q, 2*M + 100)
        for p in primerange(3, M + 1):
            if p*p <= 2*h0: continue
            pred = refined_pred(eta, n, info, A, p); act = f0.get(p, 0)
            typ, vmin, nlead = ('none', None, 0)
            if act < pred:
                typ, vmin, nlead = event_type(eta, n, info, A, a0, p)
            rows.append(dict(n=n, p=p, x=n/p, fl=n//p, n0=n % p, pmod3=p % 3, pmod4=p % 4, pred=pred, act=act, gap=pred - act, type=typ, nlead=nlead))
        ev = [r for r in rows if r['n'] == n and r['gap'] > 0]
        print('n=%d: %d primes in range, %d events: %s' % (n, sum(1 for r in rows if r['n'] == n), len(ev), ' '.join('%d(x=%.2f,%s,gap %d)' % (r['p'], r['x'], r['type'][:3], r['gap']) for r in ev)))
        sys.stdout.flush()
    if out:
        import json
        with open(out, 'w') as fh: json.dump(rows, fh)
    # summary
    tot = len(rows); ev = [r for r in rows if r['gap'] > 0]
    print('\nTOTAL: %d (n,p) pairs, %d events (%.1f%%); gap>1: %d; exceed (act>pred): %d' % (tot, len(ev), 100*len(ev)/max(tot, 1), sum(1 for r in ev if r['gap'] > 1), sum(1 for r in rows if r['gap'] < 0)))
    def rate(key, vals):
        parts = []
        for v in vals:
            sub = [r for r in rows if key(r) == v]
            if sub: parts.append('%s: %d/%d = %.2f' % (v, sum(1 for r in sub if r['gap'] > 0), len(sub), sum(1 for r in sub if r['gap'] > 0)/len(sub)))
        return '  '.join(parts)
    print('by floor(n/p): ' + rate(lambda r: r['fl'], sorted(set(r['fl'] for r in rows))))
    print('by p mod 3:   ' + rate(lambda r: r['pmod3'], [1, 2]))
    print('by p mod 4:   ' + rate(lambda r: r['pmod4'], [1, 3]))
    print('by n0 = n mod p side: ' + rate(lambda r: 'n0<p/2' if r['n0'] < r['p']/2 else 'n0>p/2', ['n0<p/2', 'n0>p/2']))
    print('by type: ' + ', '.join('%s: %d' % (t, sum(1 for r in ev if r['type'] == t)) for t in ('per-order', 'cross', 'other')))
    print('by frac(n/p) tenths: ' + rate(lambda r: int(10*(r['x'] - r['fl'])), list(range(10))))
    return rows

def _vp_int(p, n):
    e = 0
    while n and n % p == 0: n //= p; e += 1
    return e
def digit_exponent(eta, n, p):
    """Theorems 1 + 5 as an algorithm: the exponent of p in den(a_0) from base-p digits alone (odd p, p^2 > 2h0, p not | n).
    Returns dict(exp=predicted exponent, e1=Theorem-1 bound, phimin, residue, lead=#leading poles, window_dominated, perorder_pairs)."""
    eta0, es = eta[0], eta[1:]; s = len(es); h0 = eta0*n + 1; m = (h0 - 1)//2
    hs = [e*n for e in es]; L = [h0 - 2*h for h in hs]; N = min(hs); W = (h0 - 2 - p)//2
    vgam = sum(((eta0 - 2*e)*n)//p for e in es[1:]) - 2*((es[0]*n)//p)
    inv = lambda a: pow(a % p, -1, p)
    def vC(k):            # Lemma B exact count, and s_k, D_k
        D = m - k
        Nk = (2*eta0*n - 2*k + p)//(2*p) + (2*k + p)//(2*p)
        v = vgam + (1 if (k != m and D % p == 0) else 0) + Nk
        sk = 0; dden = False
        for h, Lj in zip(hs, L):
            cover = h <= k <= h + Lj - 1; sk += cover
            lo, hi = h - k, h + Lj - 1 - k
            cnt = hi//p - (lo - 1)//p          # multiples of p in [lo, hi]
            if lo <= 0 <= hi: cnt -= 1          # the zero (pole) is not a factor
            v -= cnt
            if cnt > 0: dden = True
        sk -= (k == m)
        Dk = 10**9 if dden else Nk + (1 if (k != m and D % p == 0) else 0)
        return v, sk, Dk
    # window poles: Phi and leading set
    Phi = {}
    for k in range(N, W + 1):
        v, sk, Dk = vC(k); Phi[k] = (v + s - sk, v, sk)
    if not Phi: return None
    phimin = min(t[0] for t in Phi.values()); Lset = [k for k, t in Phi.items() if t[0] == phimin]
    ps_term = s - phimin
    # non-PS term over all poles
    np_term = -10**9
    for k in range(N, h0 - N):
        v, sk, Dk = vC(k)
        np_term = max(np_term, min(sk - 1, Dk) - v)
    e1 = max(ps_term, np_term)
    window_dominated = ps_term > np_term
    # unit(C_k) mod p for all window poles by the boundary recurrence, starting from a direct product at k = N
    def unit_of(x):        # x nonzero integer or Fraction-like (num, den): returns (v_p, unit mod p)
        if isinstance(x, tuple):
            a, b = x; va, ua = unit_of(a); vb, ub = unit_of(b); return va - vb, (ua*inv(ub)) % p
        v = _vp_int(p, abs(x)); u = (x // p**v) % p; return v, u
    k = N; v0 = 0; u0 = 1
    # gamma: 4^{h0-1} prod_{j>=2} (L_j-1)! / (h_1*)!^2  -> unit via p-adic factorial units
    def ufact(M):          # (v_p(M!), unit(M!) mod p) by Anton/Wilson: M! = p^{floor(M/p)} floor(M/p)! prod_{v<=M, p not|v} v, and prod_{v<=M,p not|v} v = (-1)^{floor(M/p)} (M mod p)! mod p
        v = 0; u = 1
        while M > 0:
            q, r = divmod(M, p)
            u = (u * pow(-1, q, p) * math.factorial(r)) % p
            v += q; M = q
        return v, u
    v0 = 0; u0 = pow(4, h0 - 1, p)
    for h, Lj in zip(hs[1:], L[1:]):
        vv, uu = ufact(Lj - 1); v0 += vv; u0 = (u0*uu) % p
    vv, uu = ufact(hs[0]); v0 -= 2*vv; u0 = (u0*inv(uu*uu % p)) % p
    # centre factor and numerator at k = N: prod_{i=1}^{h0-1} (2i-2k-1)/2 ; denominator families prod (h + i - k), i in [0, L-1], nonzero
    D = m - k
    vv, uu = unit_of(2*D); v0 += vv; u0 = (u0*uu) % p
    # numerator: odd o from 1-2k to 2h0-3-2k: use products of odd numbers via double factorials of units: do it directly but with p-free skipping (O(h0), once)
    for o in range(1 - 2*k, 2*h0 - 2 - 2*k, 2):
        vv, uu = unit_of(o); v0 += vv; u0 = (u0*uu) % p
    u0 = (u0*pow(inv(2), h0 - 1, p)) % p
    for h, Lj in zip(hs, L):
        for i in range(Lj):
            d = h + i - k
            if d: vv, uu = unit_of(d); v0 -= vv; u0 = (u0*inv(uu)) % p
    units = {N: (v0, u0)}
    vk, uk = v0, u0
    for k in range(N, W):
        # move k -> k+1: centre 2D -> 2(D-1); numerator gains alpha_0(k) = -(2k+1)/2, loses alpha_{h0-1}(k) = (2h0-3-2k)/2;
        # each denominator family gains beta_{j,-1}(k) = h-1-k (if nonzero), loses beta_{j,L-1}(k) = h+L-1-k (if nonzero)
        D = m - k
        if D - 1 != 0 and D != 0:
            vv, uu = unit_of(2*(D - 1)); vk += vv; uk = (uk*uu) % p
            vv, uu = unit_of(2*D); vk -= vv; uk = (uk*inv(uu)) % p
        vv, uu = unit_of(-(2*k + 1)); vk += vv; uk = (uk*uu) % p
        vv, uu = unit_of(2*h0 - 3 - 2*k); vk -= vv; uk = (uk*inv(uu)) % p
        for h, Lj in zip(hs, L):
            g = h - 1 - k
            if g: vv, uu = unit_of(g); vk -= vv; uk = (uk*inv(uu)) % p
            l = h + Lj - 1 - k
            if l: vv, uu = unit_of(l); vk += vv; uk = (uk*uu) % p
        units[k + 1] = (vk, uk)
    # sanity: recurrence valuation must match Lemma B's count
    for k in Lset:
        assert units[k][0] == Phi[k][1], ('valuation mismatch at k', k, units[k][0], Phi[k][1])
    # residue
    R = 0; cells = {}
    for k in Lset:
        sk = Phi[k][2]; D = m - k
        # digits
        qlo = -((2*k - 1)//p)            # smallest q with q p >= 1-2k  -> ceil((1-2k)/p)
        qlo = (1 - 2*k + p - 1)//p if (1 - 2*k) > 0 else -((2*k - 1)//p)
        # simpler: enumerate q from ceil((1-2k)/p) to floor((2h0-3-2k)/p), odd only
        qa = math.ceil((1 - 2*k)/p); qb = math.floor((2*h0 - 3 - 2*k)/p)
        A = [q for q in range(qa, qb + 1) if q % 2]
        if D % p == 0: A.append(2*(D//p))          # centre factor digit in the 2*alpha/p convention (even q allowed)
        B = []
        for h, Lj in zip(hs, L):
            blo = math.ceil((h - k)/p); bhi = math.floor((h + Lj - 1 - k)/p)
            B += [b for b in range(blo, bhi + 1) if b != 0]
        key = (tuple(A), tuple(sorted(B)), sk)
        if key not in cells:
            # Taylor coefficients mod p of F(w) = prod_A (1 + 2w/q) / prod_B (1 + w/b), degree < sk
            num = [1]
            for q in A:
                c = (2*inv(q)) % p; new = [0]*(len(num) + 1)
                for i, x in enumerate(num): new[i] = (new[i] + x) % p; new[i+1] = (new[i+1] + x*c) % p
                num = new[:sk]
            den = [1]
            for b in B:
                c = inv(b); new = [0]*(len(den) + 1)
                for i, x in enumerate(den): new[i] = (new[i] + x) % p; new[i+1] = (new[i+1] + x*c) % p
                den = new[:sk]
            num += [0]*(sk - len(num)); den += [0]*(sk - len(den))
            hcoef = [0]*sk
            for r in range(sk):
                acc = num[r]
                for j in range(1, r + 1): acc -= den[j]*hcoef[r - j]
                hcoef[r] = acc % p
            cells[key] = hcoef
        hcoef = cells[key]
        J = (2*D - 1 + p)//(2*p)
        T = 0
        for q in range(1, 2*J, 2):
            eps = -1 if ((q*p + 1)//2) % 2 else 1
            w = (-q*inv(2)) % p
            poly = 0
            for r in range(sk - 1, -1, -1): poly = (poly*w + hcoef[r]) % p
            T = (T + eps*pow((-2*inv(q)) % p, sk, p)*poly) % p
        R = (R + (-1)**k*units[k][1]*T) % p
    exp = ps_term - (1 if R == 0 else 0)
    return dict(exp=max(exp, np_term) if not window_dominated else exp, e1=e1, phimin=phimin, residue=R, lead=len(Lset), window_dominated=window_dominated, ps_term=ps_term, np_term=np_term)
def validate_digit(json_files):
    import json
    for name, eta in json_files:
        rows = json.load(open(name)); ok = 0; tot = 0; bad = []
        for r in rows:
            if r['n'] % r['p'] == 0: continue
            d = digit_exponent(eta, r['n'], r['p'])
            if d is None: continue
            tot += 1
            pred = d['exp']; act = r['act']
            good = (act == pred) or (d['residue'] == 0 and act < d['e1'] and act <= pred)
            if good: ok += 1
            else: bad.append((r['n'], r['p'], act, pred, d['e1'], d['residue'], d['window_dominated']))
        print('%s: digit-only exponent agrees with the exact exponent in %d/%d pairs' % (name, ok, tot))
        if bad: print('   disagreements (n,p,actual,predicted,e1,residue,window_dom):', bad[:10])
        sys.stdout.flush()

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'validate_digit':
        validate_digit([('sweep_s5_simple.json', [3,1,1,1,1,1]), ('sweep_s7_simple.json', [3,1,1,1,1,1,1,1]), ('sweep_s9_simple.json', [3,1,1,1,1,1,1,1,1,1]),
                        ('sweep_s11_winner_b.json', [33,9,10,10,11,11,12,12,12,13,13,14]), ('sweep_s11_eta36.json', [36,10,10,11,11,12,12,13,13,14,14,15]), ('sweep_s13_zudilin.json', [31,10,10,10,10,10,11,11,11,11,12,12,12,12])]); sys.exit()
    if cmd == 'digitexp':
        eta = parse_eta(sys.argv[2]); n = int(sys.argv[3]); p = int(sys.argv[4]); print(digit_exponent(eta, n, p)); sys.exit()
    if cmd == 'mech2':
        mech2(parse_eta(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])); sys.exit()
    if cmd == 'sweep':
        sweep(parse_eta(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6] if len(sys.argv) > 6 else None); sys.exit()
    if cmd == 'forms':
        eta = parse_eta(sys.argv[2]); n = int(sys.argv[3]); chk = not (len(sys.argv) > 4 and sys.argv[4] == 'nocheck')
        run_forms(eta, n, chk)
    elif cmd == 'asym':
        asym(parse_eta(sys.argv[2]))
    elif cmd == 'search':
        search(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else None, sys.argv[5] if len(sys.argv) > 5 else 'net_ref')
