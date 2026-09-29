"""cubic_asym.py — asymptotic ledger for the cubic (χ₋₃) forms of cubic.py, Lai–Zhou style (arXiv 2103.00904, §3–4, half block removed).

Normalised R̂(t) = 3^{3L} ∏_j ((m2−2δ_j)n)! / n!^{6m1+2m2} · (2t+m2 n) ∏_{i=1}^{m1 n}(t−i)(t+m2 n+i) (t−m1 n+1/3)_L (t−m1 n+2/3)_L / ∏_j (t+δ_j n)_{(m2−2δ_j)n+1},
L = (2m1+m2) n, q = #δ odd.  Forms Ŝ_{1/3} − Ŝ_{2/3} ∈ Q + Σ_{i even ≤ q−1} Q·L(i, χ₋₃).
Arithmetic (Lai–Zhou Lemma 3.1/3.2 with the half-block terms removed): Φ_n^{-1} D_M^{q} (coefficients) ∈ Z, M = (m2 − 2 δ_min) n,
  ν(x,y) = Σ_j (⌊(m2−2δ_j)x⌋ − ⌊y−δ_j x⌋ − ⌊(m2−δ_j)x − y⌋) + ⌊3m1 x + 3y⌋ + ⌊3(m1+m2)x − 3y⌋ − ⌊y⌋ − ⌊m2 x − y⌋ − (6m1+2m2)⌊x⌋,
  ν0(x) = min_y ν(x,y),  log Φ_n / n → κ = Σ_pieces ν0·(ψ(b) − ψ(a)) − ∫_0^{1/M'} ν0 dx/x²  (M' = m2 − 2δ_min).
Analysis (their Lemma 4.2 with three numerator blocks): lim |Ŝ|^{1/n} = g(x0), f(x0) = 1,
  f(X) = ((2m1+m2+X)/X)^3 (m1+X)/(m1+m2+X) ∏_j (m1+δ_j+X)/(m1+m2−δ_j+X),
  g(X) = 27^{2m1+m2} ∏_j (m2−2δ_j)^{m2−2δ_j} (2m1+m2+X)^{3(2m1+m2)} (m1+X)^{m1}/(m1+m2+X)^{m1+m2} ∏_j (m1+δ_j+X)^{m1+δ_j}/(m1+m2−δ_j+X)^{m1+m2−δ_j}.
net per n = q M' − κ + log g(x0)  (needs < 0).
usage: python cubic_asym.py m1 m2 d1,...,dq        |   python cubic_asym.py search q m2min m2max
"""
import sys, math
from fractions import Fraction
import numpy as np
from mpmath import mp, mpf, psi as mpsi, log as mlog, findroot

def nu_np(m1, m2, deltas, x, y):
    v = np.floor(3*m1*x + 3*y) + np.floor(3*(m1 + m2)*x - 3*y) - np.floor(y) - np.floor(m2*x - y) - (6*m1 + 2*m2)*math.floor(x)
    for d in deltas: v += math.floor((m2 - 2*d)*x) - np.floor(y - d*x) - np.floor((m2 - d)*x - y)
    return v
def nu0(m1, m2, deltas, x):
    b = [0.0, (m2*x) % 1] + [(-3*m1*x + k/3) % 1 for k in range(3)] + [(3*(m1 + m2)*x - k/3) % 1 for k in range(3)]
    for d in deltas: b += [(d*x) % 1, ((m2 - d)*x) % 1]
    b = np.unique(np.array(b)); b2 = np.append(b, b[0] + 1); mids = (b2[:-1] + b2[1:])/2
    return int(nu_np(m1, m2, deltas, x, mids).min())
_F = {}
def farey(D):
    if D not in _F:
        pts = set()
        for d in range(1, D + 1):
            for m in range(0, d + 1): pts.add(Fraction(m, d))
        _F[D] = sorted(pts)
    return _F[D]
def kappa(m1, m2, deltas, verbose=False):
    Mp = m2 - 2*min(deltas)
    D = 3*max(3*(m1 + m2), m2) + 3
    pts = farey(D)
    mp.dps = 25
    kap = mpf(0); pieces = []
    for a, b in zip(pts[:-1], pts[1:]):
        v = nu0(m1, m2, deltas, float((a + b)/2))
        if v:
            fa, fb = mpf(a.numerator)/a.denominator, mpf(b.numerator)/b.denominator
            kap += v*(mpsi(0, fb) - mpsi(0, fa))            # sum over m >= 0 of the period pieces (includes [0,1))
            lo = Fraction(1, Mp)
            if a < lo:                                        # remove x < 1/M' (primes > M')
                aa = fa; bb = min(fb, mpf(lo.numerator)/lo.denominator)
                if bb > aa: kap -= v*(1/aa - 1/bb) if a > 0 else v*(mpf(10)**30)   # a = 0 piece: nu0 must be 0 there
        pieces.append((a, b, v))
    if verbose:
        comp = []
        for a, b, v in pieces:
            if comp and comp[-1][2] == v and comp[-1][1] == a: comp[-1] = (comp[-1][0], b, v)
            else: comp.append((a, b, v))
        print('  nu0 steps: ' + ' '.join('[%s,%s)=%d' % (a, b, v) for a, b, v in comp if v))
    return float(kap)
def decay(m1, m2, deltas):
    mp.dps = 30
    def f(X): return ((2*m1 + m2 + X)/X)**3*(m1 + X)/(m1 + m2 + X)*math.prod([(m1 + d + X)/(m1 + m2 - d + X) for d in deltas])
    lo, hi = mpf('1e-9'), mpf('1e9')
    for _ in range(300):
        mid = (lo*hi)**0.5
        if f(mid) > 1: lo = mid
        else: hi = mid
    x0 = (lo*hi)**0.5
    lg = (2*m1 + m2)*mlog(27) + sum((m2 - 2*d)*mlog(m2 - 2*d) for d in deltas) + 3*(2*m1 + m2)*mlog(2*m1 + m2 + x0) + m1*mlog(m1 + x0) - (m1 + m2)*mlog(m1 + m2 + x0)
    for d in deltas: lg += (m1 + d)*mlog(m1 + d + x0) - (m1 + m2 - d)*mlog(m1 + m2 - d + x0)
    return float(lg), float(x0)
def ledger(m1, m2, deltas, verbose=True):
    q = len(deltas); Mp = m2 - 2*min(deltas)
    assert q % 2 == 1 and (q - 2)*m2 - 6*m1 - 2*sum(deltas) >= 0 and all(0 <= d < m2/2 for d in deltas)
    lg, x0 = decay(m1, m2, deltas)
    kap = kappa(m1, m2, deltas, verbose=verbose)
    naive = q*Mp; height = naive - kap; net = height + lg
    if verbose:
        print('cubic ledger m1=%d m2=%d deltas=%s (q=%d, L-values L(2..%d, chi_-3)): closeness log g(x0) = %.4f (x0 = %.5f) ; height = %d (q M\') - %.4f (kappa) = %.4f ; NET = %+.4f per n' % (m1, m2, deltas, q, q - 1, lg, x0, naive, kap, height, net))
    return dict(closeness=lg, kappa=kap, height=height, net=net, x0=x0)

def search(q, m2min, m2max):
    best_all = None
    for m2 in range(m2min, m2max + 1):
        best = None
        for m1 in range(1, max(2, (q - 2)*m2//6) + 1):
            room = (q - 2)*m2 - 6*m1
            if room < 0: continue
            # coordinate descent over deltas starting from equal deltas
            d0 = max(0, min((m2 - 1)//2, room//(2*q) ))
            cur = [d0]*q
            def val(ds):
                if sum(ds)*2 > room or any(d < 0 or 2*d >= m2 for d in ds): return None
                return ledger(m1, m2, ds, verbose=False)['net']
            cv = val(cur)
            if cv is None: continue
            improved = True
            while improved:
                improved = False
                for j in range(q):
                    for dl in (1, -1):
                        cand = cur[:]; cand[j] += dl; cand = sorted(cand)
                        v = val(cand)
                        if v is not None and v < cv - 1e-9: cur, cv, improved = cand, v, True
                for j in range(q):
                    for l in range(q):
                        if j == l: continue
                        cand = cur[:]; cand[j] += 1; cand[l] -= 1; cand = sorted(cand)
                        v = val(cand)
                        if v is not None and v < cv - 1e-9: cur, cv, improved = cand, v, True
            if best is None or cv < best[0]: best = (cv, m1, cur)
        if best:
            r = ledger(best[1], m2, best[2], verbose=False)
            print('q=%d m2=%2d: best net %+.4f per n (m1=%d, deltas=%s; closeness %.3f, height %.3f = %d - %.3f) ; per unit m2: %+.4f' % (q, m2, best[0], best[1], best[2], r['closeness'], r['height'], q*(m2 - 2*min(best[2])), r['kappa'], best[0]/m2))
            sys.stdout.flush()
            if best_all is None or best[0]/m2 < best_all[0]: best_all = (best[0]/m2, m2, best[1], best[2])
    print('best per unit m2: %+.4f at m2=%d m1=%d deltas=%s' % best_all)

if __name__ == '__main__':
    if sys.argv[1] == 'search':
        search(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])); sys.exit()
    m1 = int(sys.argv[1]); m2 = int(sys.argv[2]); deltas = [int(x) for x in sys.argv[3].split(',')]
    ledger(m1, m2, deltas)
