"""zfast.py — fast evaluator for Zudilin's (2004, §8) odd-zeta forms with general odd r, for searching directions eta.

Objective (his Proposition 5 generalised to odd r, arithmetic from Lemma 19 which is stated for all odd r):
  net = C2 - C0,   C2 = r m_1 + m_2 + ... + m_{q-r} - kappa,   C0 = -Re f0(tau0),
  kappa = int_0^1 phi(x) dpsi(x) - int_0^{1/m_{q-r}} phi(x) dx/x^2  (= int_{1/m}^{inf} phi(x) dx/x^2 by periodicity),
computed on a uniform u = 1/x grid (as lz.fast_C1) instead of the exact Farey grid; the saddle polynomial is built with numpy.
Refined rebate (zetaforms.refined_rebate) is added optionally.  Validated against zetaforms.asym at Zudilin's r=3, q=13 point.
"""
import sys, math, time
import numpy as np
from mpmath import mp, mpf, mpc, log, psi as mpsi

FAST = {'du': 0.01, 'nx': 4000}

def m_params(r, eta):
    eta0, es = eta[0], list(eta[1:]); q = len(es)
    m0 = max(es[r - 1], eta0 - 2*es[r])
    return [max(m0, eta0 - es[0] - es[r + j - 1]) for j in range(1, q - r + 1)]

def _polymul(a, b):
    out = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b): out[i + j] += x*y
    return out

def saddle_poly(r, eta):
    """exact integer coefficients (highest first) of P(tau) = (tau-eta0)^r prod(tau-eta_j) - tau^r prod(tau-eta0+eta_j)."""
    eta0, es = eta[0], list(eta[1:])
    P1 = [1]; P2 = [1]
    for _ in range(r): P1 = _polymul(P1, [1, -eta0]); P2 = _polymul(P2, [1, 0])
    for e in es: P1 = _polymul(P1, [1, -e]); P2 = _polymul(P2, [1, -(eta0 - e)])
    P = [a - b for a, b in zip(P1, P2)]
    while P and P[0] == 0: P = P[1:]
    return P

def f0(r, eta, t):
    eta0, es = eta[0], list(eta[1:])
    v = r*eta0*log(eta0 - t)
    for e in es: v += e*log(t - e) - (eta0 - e)*log(t - eta0 + e)
    for e in es[:r]: v -= 2*e*log(e)
    for e in es[r:]: v += (eta0 - 2*e)*log(eta0 - 2*e)
    return v

_ROOT_CACHE = {}
def all_roots(r, eta):
    """all roots of the saddle polynomial, high precision (mpmath polyroots on exact integer coefficients)."""
    key = (r, tuple(eta))
    if key in _ROOT_CACHE: return _ROOT_CACHE[key]
    P = saddle_poly(r, eta)
    deg = len(P) - 1
    mp.dps = max(50, 2*deg)
    roots = None
    for extra in (100, 300, 800):
        try:
            roots = mp.polyroots([mpf(c) for c in P], maxsteps=400 + 20*deg, extraprec=extra)
            break
        except Exception:
            continue
    if roots is None:
        _ROOT_CACHE[key] = None; return None
    roots = [complex(z) for z in roots]
    _ROOT_CACHE[key] = roots
    return roots

def newton_roots(r, eta, dps=40):
    """roots with Im > 0 found by Newton iteration (exact integer coefficients, mp precision) from a fan of starting points
    around eta0 (where the r-fold cluster sits), along the real axis beyond eta0, and on the symmetry line Re = eta0/2."""
    import cmath
    P = saddle_poly(r, eta)
    mp.dps = dps
    eta0 = eta[0]
    coeffs = [mpf(c) for c in P]
    def pv(t):
        p = mpc(0); d = mpc(0)
        for c in coeffs: d = d*t + p; p = p*t + c
        return p, d
    starts = [eta0 + rho*cmath.exp(1j*th) for rho in (0.5, 2.0, 6.0) for th in (0.5, 1.2, 1.9, 2.6)]
    starts += [eta0*0.75 + 3j, eta0*0.75 + 15j, eta0*0.9 + 8j, eta0 + 15.0 + 8j]
    found = []
    tol = mpf(10)**(-12)
    for s in starts:
        t = mpc(s.real, s.imag); ok = False
        for _ in range(40):
            p, d = pv(t)
            if d == 0: break
            step = p/d; t = t - step
            if abs(step) < tol: ok = True; break
        if ok and t.imag > 1e-9 and all(abs(t - f) > 1e-6 for f in found): found.append(t)
    return [complex(f) for f in found]

def cluster_roots(r, eta, dps=40):
    """The (r-1)/2 roots with Im > 0 of the r-fold cluster near eta0 (the roots with Re > eta0/2), found by Newton iteration from
    starting points on circles of the cluster radius rho = eta0 (prod eta_j / prod (eta0 - eta_j))^{1/r} around eta0.
    If the count is not (r-1)/2 the exact all-roots solver is used instead.  (A plain Newton fan missed the max-Re root at
    (130; 35,36,...,56): it found 125.16+2.44i but not 129.30+7.25i, overstating C0 by 30.7.)"""
    import cmath
    eta0, es = eta[0], list(eta[1:])
    lr = r*math.log(eta0) + sum(math.log(e) - math.log(eta0 - e) for e in es)
    rho = math.exp(lr/r)
    P = saddle_poly(r, eta)
    mp.dps = dps
    coeffs = [mpf(c) for c in P]
    def pv(t):
        p = mpc(0); d = mpc(0)
        for c in coeffs: d = d*t + p; p = p*t + c
        return p, d
    starts = [eta0 + rho*f*cmath.exp(1j*(math.pi*k/r + 0.13)) for f in (0.5, 1.0, 1.7, 2.8) for k in range(2*r)]
    found = []
    tol = mpf(10)**(-12)
    for s in starts:
        t = mpc(s.real, s.imag); ok = False
        for _ in range(60):
            p, d = pv(t)
            if d == 0: break
            step = p/d; t = t - step
            if abs(step) < tol: ok = True; break
        if ok and t.imag > 1e-9 and t.real > eta0/2 + 1e-6 and all(abs(t - f) > 1e-6 for f in found): found.append(t)
    found = [complex(f) for f in found]
    expected = (r - 1)//2
    if len(found) != expected:
        roots = all_roots(r, eta)
        if not roots: return None
        found = [z for z in roots if z.imag > 1e-9 and z.real > eta0/2 + 1e-6]
    return found

def saddle(r, eta, full=False):
    """tau0: root with Im > 0 and maximal real part (Zudilin's rule); C0 = -Re f0(tau0).  full=True uses all roots (slow, verification)."""
    if full:
        roots = all_roots(r, eta)
        cand = [z for z in roots if z.imag > 1e-9] if roots else []
    else:
        cand = cluster_roots(r, eta)
    if not cand: return None, None
    tau0 = max(cand, key=lambda z: z.real)
    mp.dps = 30
    v = f0(r, eta, mpc(tau0.real, tau0.imag))
    return tau0, -float(v.real)

def phi_func(r, eta):
    eta0, es = eta[0], list(eta[1:])
    num = np.array(es[:r], dtype=float); den = np.array(es[r:], dtype=float)
    def phi0(x, y):
        # y array
        v = np.zeros_like(y)
        for e in num:
            v += np.floor(y) + np.floor(eta0*x - y) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y) - 2*math.floor(e*x)
        for e in den:
            v += math.floor((eta0 - 2*e)*x) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y)
        return v
    allc = np.concatenate([[0.0, eta0], num, eta0 - num, den, eta0 - den])
    def phi(x):
        b = np.unique(np.mod(allc*x, 1.0)); b2 = np.append(b, b[0] + 1.0); mids = (b2[:-1] + b2[1:])/2
        return float(phi0(x, mids).min())
    return phi

_PSI_CACHE = {}
def _psi_weights(nx):
    if nx not in _PSI_CACHE:
        xs = (np.arange(nx) + 0.5)/nx
        w = np.array([float(mpsi(1, 1.0 + x)) for x in xs[::100]])
        _PSI_CACHE[nx] = (xs, np.interp(xs, xs[::100], w))
    return _PSI_CACHE[nx]

def phi_vec(r, eta, xs):
    """phi(x) = min_y phi0(x, y) for an array of x, vectorised: breakpoints of y are frac(c x) for the 2q+2 shifts c; evaluate at midpoints."""
    eta0, es = eta[0], list(eta[1:])
    num = np.array(es[:r], dtype=float); den = np.array(es[r:], dtype=float)
    allc = np.concatenate([[0.0, eta0], num, eta0 - num, den, eta0 - den])
    xs = np.asarray(xs, dtype=float)
    B = np.sort(np.mod(np.outer(xs, allc), 1.0), axis=1)              # (nx, nc)
    B2 = np.concatenate([B, B[:, :1] + 1.0], axis=1)
    Y = (B2[:, :-1] + B2[:, 1:])/2                                       # midpoints (nx, nc)
    X = xs[:, None]
    v = np.zeros_like(Y)
    fl = np.floor
    for e in num:
        v += fl(Y) + fl(eta0*X - Y) - fl(Y - e*X) - fl((eta0 - e)*X - Y) - 2*fl(e*X)
    for e in den:
        v += fl((eta0 - 2*e)*X) - fl(Y - e*X) - fl((eta0 - e)*X - Y)
    return v.min(axis=1)

def kappa(r, eta, du=None, nx=None):
    du = du or FAST['du']; nx = nx or FAST['nx']
    ms = m_params(r, eta); M = ms[-1]
    us = np.arange(1.0 + du/2, M, du)
    part1 = phi_vec(r, eta, 1.0/us).sum()*du
    xs, w = _psi_weights(nx)
    part2 = float((phi_vec(r, eta, xs)*w).sum())/nx
    return part1 + part2

def fast_rebate(r, eta, du=None):
    """u-grid version of zetaforms.refined_rebate: int (E_Z(his) - E_ref) dx/x^2 over x = 1/u, u in [1, m_1 + 2] (same exponent functions)."""
    du = du or FAST.get('du_reb', 0.02)
    eta0, es = eta[0], list(eta[1:]); q = len(es)
    m0 = max(es[r - 1], eta0 - 2*es[r]); ms = [max(m0, eta0 - es[0] - es[r + j - 1]) for j in range(1, q - r + 1)]
    num = np.array(es[:r], dtype=float); den = np.array(es[r:], dtype=float)
    def pieces(x):
        lo_y, hi_y = es[r]*x, (eta0 - es[r])*x
        bases = np.concatenate([[0.0, eta0*x, es[0]*x + 1.0], num*x, (eta0 - num)*x, den*x, (eta0 - den)*x])
        pts = [np.array([lo_y, hi_y])]
        for base in bases:
            k0 = math.ceil(lo_y - base); k1 = math.floor(hi_y - base)
            if k1 >= k0: pts.append(base + np.arange(k0, k1 + 1))
        b = np.unique(np.concatenate(pts))
        return (b[:-1] + b[1:])/2
    def phi0(x, y):
        v = np.zeros_like(y)
        for e in num: v += np.floor(y) + np.floor(eta0*x - y) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y) - 2*math.floor(e*x)
        for e in den: v += math.floor((eta0 - 2*e)*x) - np.floor(y - e*x) - np.floor((eta0 - e)*x - y)
        return v
    def both(x):
        y = pieces(x)
        if len(y) == 0: return 0.0, 0.0
        S = np.zeros_like(y); dv = np.zeros_like(y, dtype=bool)
        for e in den:
            cov = (y >= e*x) & (y <= (eta0 - e)*x); S += cov
            dv |= cov & ((np.floor((eta0 - e)*x - y) + np.floor(y - e*x)) >= 1)
        D = np.zeros_like(y)
        for e in num: D += np.floor(y) - np.floor(y - e*x) + np.floor(eta0*x - y) - np.floor((eta0 - e)*x - y)
        ph0 = phi0(x, y)
        ph = ph0 - ((q - r) - S)
        win = (y - es[0]*x) >= 1.0
        cost = np.where(win, S + r - 1 - ph, np.where(dv, S - 1 - ph, np.minimum(S - 1, D) - ph))
        cost = np.where(S >= 1, cost, -1e9)
        er = float(cost.max())
        inside = S >= 1
        phmin = float(ph0[inside].min()) if inside.any() else 0.0
        dcount = r*(1 if x >= 1/ms[0] else 0) + sum(1 for M in ms[1:] if x >= 1/M)
        eh = float(dcount) if x < 1/ms[-1] else dcount - phmin
        return eh, er
    us = np.arange(1.0 + du/2, ms[0] + 2.0, du)
    reb = 0.0
    for u in us:
        eh, er = both(1.0/u)
        reb += (eh - er)*du
    return reb

def valid(r, eta):
    eta0, es = eta[0], eta[1:]; q = len(es)
    return all(1 <= es[i] <= es[i + 1] for i in range(q - 1)) and 2*es[-1] < eta0 and sum(es) <= eta0*(q - r)//2

def net(r, eta, refined=False, parts=False, full=False):
    if not valid(r, eta): return None
    tau0, C0 = saddle(r, eta, full=full)
    if C0 is None: return None
    ms = m_params(r, eta)
    kap = kappa(r, eta)
    Dpart = r*ms[0] + sum(ms[1:])
    C2 = Dpart - kap
    v = C2 - C0
    reb = 0.0
    if refined == 'exact':
        import zetaforms as zf
        reb = zf.refined_rebate(r, eta)[0]
        v -= reb
    elif refined:
        reb = fast_rebate(r, eta)
        v -= reb
    if parts: return v, C0, Dpart, kap, reb, tau0
    return v

def search(r, eta, sweeps=6, refined=False, moves='single,shift,pair', verbose=True):
    """coordinate descent on net: single +-1 moves on each eta_j (j>=1), +-1/2 on eta0, block shifts of all denominator or all numerator directions; optional pair moves."""
    cur = [eta[0]] + sorted(eta[1:]); cv = net(r, cur, refined)
    q = len(cur) - 1
    if verbose: print('seed r=%d q=%d: net %+.4f at eta=%s' % (r, q, cv, cur)); sys.stdout.flush()
    for sw in range(sweeps):
        improved = False
        cands = []
        if 'single' in moves:
            for j in range(1, q + 1):
                for dl in (1, -1):
                    c = cur[:]; c[j] += dl; cands.append([c[0]] + sorted(c[1:]))
            for d0 in (1, -1, 2, -2):
                c = cur[:]; c[0] += d0; cands.append(c)
        if 'shift' in moves:
            for dl in (1, -1):
                c = cur[:]; c[1:r + 1] = [e + dl for e in c[1:r + 1]]; cands.append([c[0]] + sorted(c[1:]))
                c = cur[:]; c[r + 1:] = [e + dl for e in c[r + 1:]]; cands.append([c[0]] + sorted(c[1:]))
        if 'pair' in moves:
            for j in range(1, q + 1):
                for l in range(j + 1, q + 1):
                    c = cur[:]; c[j] += 1; c[l] -= 1; cands.append([c[0]] + sorted(c[1:]))
                    c = cur[:]; c[j] -= 1; c[l] += 1; cands.append([c[0]] + sorted(c[1:]))
        seen = set()
        for c in cands:
            key = tuple(c)
            if key in seen or not valid(r, c): continue
            seen.add(key)
            v = net(r, c, refined)
            if v is not None and v < cv - 1e-6:
                cur, cv, improved = c, v, True
                if verbose: print('  sweep %d: %+.4f at eta=%s' % (sw, cv, cur)); sys.stdout.flush()
        if not improved: break
    if verbose:
        v, C0, Dpart, kap, reb, tau0 = net(r, cur, True, parts=True)
        print('best r=%d q=%d: refined net %+.4f (Zudilin net %+.4f, rebate %.4f) at eta=%s ; C0 %.3f, D-part %d, kappa %.3f, tau0 %s' % (r, q, v, v + reb, reb, cur, C0, Dpart, kap, tau0))
    return cv, cur

if __name__ == '__main__':
    cmd = sys.argv[1]; r = int(sys.argv[2]); eta = [int(x) for x in sys.argv[3].split(',')]
    if cmd == 'net':
        refmode = 'exact' if 'refexact' in sys.argv else ('ref' in sys.argv)
        t0 = time.time(); v, C0, Dpart, kap, reb, tau0 = net(r, eta, refined=refmode, parts=True, full=('full' in sys.argv))
        print('r=%d q=%d eta=%s : C0 %.4f ; D-part %d ; kappa %.4f ; C2 %.4f ; net %+.4f ; rebate %.4f ; tau0 %s  [%.1fs]' % (r, len(eta) - 1, eta, C0, Dpart, kap, Dpart - kap, v, reb, tau0, time.time() - t0))
    elif cmd == 'search':
        sw = 6; moves = 'single,shift'
        for a in sys.argv[4:]:
            if a.startswith('sweeps='): sw = int(a[7:])
            if a.startswith('moves='): moves = a[6:]
        search(r, eta, sweeps=sw, refined=('ref' in sys.argv), moves=moves)
