"""siegel.py — Siegel-disk / continued-fraction portrait of a constant.

For theta in (0,1): continued fraction [0; a1, a2, ...], convergents p_n/q_n, Gauss-map iterates theta_n, beta_n = theta_0...theta_n,
  Brjuno function  B(theta) = sum_{n>=0} beta_{n-1} log(1/theta_n)  (finite  <=>  the quadratic polynomial P_theta(z) = e^{2 pi i theta} z + z^2
  has a Siegel disk at 0 with rotation number theta; Yoccoz: |log r(theta) + B(theta)| <= C, r = conformal radius of the disk),
  Khinchin (geometric mean of a_n -> 2.6854...), Levy (q_n^{1/n} -> e^{pi^2/(12 ln 2)} = 3.2758...), Gauss-Kuzmin frequencies,
  Dedekind sums s(p_n, q_n) of the convergents (exact, via reciprocity) against the alternating sum of partial quotients (Barkan-Hickerson).
Siegel disk: linearisation h(lambda z) = P(h(z)), h(z) = z + sum h_n z^n, h_n = c sum_{k} g_k g_{n-k}/(lambda^n - lambda) with rescaling c;
  conformal radius r = 1/limsup |h_n|^{1/n}.
usage:  python siegel.py cf NAME [terms]          python siegel.py disk NAME [N]          python siegel.py dedekind NAME [n]
NAME in: G, 2G/pi, G/pi, ln(2+sqrt3), R/3 (Ramanujan number (8G - pi ln(2+sqrt3))/3), e^{2G/pi}, golden, sqrt2, e, pi, zeta3, L2chi3
"""
import sys, math
from fractions import Fraction
from mpmath import mp, mpf, mpc, pi, catalan, log, sqrt, exp, zeta, e as E, floor, nstr, expj

def constant(name, dps):
    mp.dps = dps
    G = catalan
    table = {
        'G': G, '2G/pi': 2*G/pi, 'G/pi': G/pi, 'ln(2+sqrt3)': log(2 + sqrt(3)), 'R/3': (8*G - pi*log(2 + sqrt(3)))/3,
        'e^{2G/pi}': exp(2*G/pi), 'golden': (sqrt(5) - 1)/2, 'sqrt2': sqrt(2), 'e': E, 'pi': pi, 'zeta3': zeta(3),
        'L2chi3': (zeta(2, mpf(1)/3) - zeta(2, mpf(2)/3))/9,
    }
    x = table[name]
    return x - floor(x)

def cf(theta, terms):
    """partial quotients a_1..a_terms of theta in (0,1) (theta_0 = theta), plus the Gauss iterates as floats."""
    a = []; th = theta; iters = []
    for _ in range(terms):
        iters.append(float(th))
        inv = 1/th; ai = int(floor(inv)); a.append(ai); th = inv - ai
        if th == 0: break
    return a, iters

def convergents(a):
    p0, q0, p1, q1 = 0, 1, 1, a[0]
    out = [(p1, q1)]
    for ai in a[1:]:
        p0, q0, p1, q1 = p1, q1, ai*p1 + p0, ai*q1 + q0
        out.append((p1, q1))
    return out

def brjuno_partial(theta, terms):
    """partial sums of B(theta) = sum beta_{n-1} log(1/theta_n) with theta_0 = theta (mp precision)."""
    B = mpf(0); beta = mpf(1); th = theta; sums = []
    for n in range(terms):
        B += beta*log(1/th); sums.append(float(B))
        beta *= th; inv = 1/th; th = inv - floor(inv)
        if th == 0: break
    return sums

def dedekind(h, k):
    """Dedekind sum s(h,k) exactly via reciprocity + Euclid (h, k > 0 coprime)."""
    s = Fraction(0); sign = 1
    h %= k
    while h:
        s += sign*(Fraction(-1, 4) + Fraction(h*h + k*k + 1, 12*h*k))
        h, k = k % h, h
        sign = -sign
    return s

def cmd_cf(name, terms):
    theta = constant(name, int(terms*0.6) + 200)
    a, iters = cf(theta, terms)
    conv = convergents(a)
    n = len(a)
    khin = math.exp(sum(math.log(x) for x in a)/n)
    levy = math.exp(math.log(conv[-1][1])/n)
    big = sorted(((ai, i + 1) for i, ai in enumerate(a)), reverse=True)[:6]
    freq = {k: sum(1 for x in a if x == k)/n for k in (1, 2, 3, 4)}
    gk = {k: math.log2((k + 1)**2/(k*(k + 2))) for k in (1, 2, 3, 4)}
    Bs = brjuno_partial(theta, min(terms, 400))
    print('%s = %s...  [%d partial quotients]' % (name, nstr(theta, 25), n))
    print('  a_1..a_30: %s' % a[:30])
    print('  Khinchin geometric mean %.4f (limit 2.6854) ; Levy q_n^(1/n) %.4f (limit e^{pi^2/(12 ln 2)} = 3.2758)' % (khin, levy))
    print('  Gauss-Kuzmin: ' + '  '.join('P(a=%d) %.3f (%.3f)' % (k, freq[k], gk[k]) for k in (1, 2, 3, 4)))
    print('  largest partial quotients (value, position): %s' % big)
    print('  Brjuno sum B(theta): %.6f after 5 terms, %.6f after 20, %.6f after 100, %.6f after %d (tail beyond is < %.1e)' % (Bs[4], Bs[19], Bs[99], Bs[-1], len(Bs), abs(Bs[-1] - Bs[-2])))
    return a, conv

def cmd_disk(name, N):
    """linearisation coefficients of P_theta(z) = lambda z + z^2, lambda = e^{2 pi i theta}; conformal radius estimate."""
    import numpy as np
    theta = constant(name, 60)
    Bs = brjuno_partial(theta, 200)
    lam = complex(expj(2*pi*theta))
    th = float(theta)
    # first pass with c = 1 up to n where overflow threatens, to estimate r; then rescaled pass
    def run(c, N):
        g = np.zeros(N + 1, dtype=complex); g[1] = 1.0
        logs = np.zeros(N + 1)
        for n in range(2, N + 1):
            div = lam**n - lam
            s = np.dot(g[1:n], g[n - 1:0:-1])
            g[n] = c*s/div
        return g
    # estimate: crude r from small N
    g = run(1.0, 300)
    mags = np.abs(g[2:301]); ns = np.arange(2, 301)
    r_est = float(np.exp(-np.max(np.log(mags)/ns)))
    c = r_est
    g = run(c, N)
    mags = np.abs(g[1:N + 1]); ns = np.arange(1, N + 1)
    # |h_n| = |g_n| / c^{n-1}  ->  log|h_n|/n = (log|g_n| - (n-1) log c)/n
    loghn = (np.log(np.maximum(mags, 1e-300)) - (ns - 1)*math.log(c))/ns
    tail = loghn[N//2:]
    r = math.exp(-float(np.max(tail)))
    # dips: positions of small divisors (n-1 = q_k) show as spikes in |g_n|
    peaks = sorted(((float(loghn[i]), int(ns[i])) for i in range(len(ns))), reverse=True)[:8]
    print('Siegel disk of P_theta, theta = %s = %.15f: lambda = e^{2 pi i theta}' % (name, th))
    print('  conformal radius estimate r = 1/limsup|h_n|^(1/n) over n in [%d, %d]: r ~ %.5f, log r ~ %.4f' % (N//2, N, r, math.log(r)))
    print('  Brjuno B(theta) = %.4f (200 terms); Yoccoz: log r + B(theta) should be bounded: %.4f' % (Bs[-1], math.log(r) + Bs[-1]))
    print('  largest log|h_n|/n (n): %s' % ', '.join('%.3f (%d)' % pk for pk in peaks))
    return r, Bs[-1]

def cmd_dedekind(name, nmax):
    theta = constant(name, int(nmax*0.6) + 200)
    a, iters = cf(theta, nmax)
    conv = convergents(a)
    print('Dedekind sums of the convergents p_n/q_n of %s versus the alternating sum of partial quotients (Barkan-Hickerson):' % name)
    print('   n   a_n   q_n(digits)   12 s(p_n,q_n)          sum_{i<=n} (-1)^{i+1} a_i    difference')
    alt = 0
    for i, (ai, (p, q)) in enumerate(zip(a, conv)):
        n = i + 1
        alt += (-1)**(n + 1)*ai
        s12 = 12*dedekind(p, q)
        if n <= 12 or n % 25 == 0 or n == len(a):
            print('  %3d  %4d   %4d   %-22s %-8d    %s' % (n, ai, len(str(q)), nstr(mpf(s12.numerator)/s12.denominator, 12), alt, nstr(mpf(s12.numerator)/s12.denominator - alt, 8)))

if __name__ == '__main__':
    cmd = sys.argv[1]; name = sys.argv[2]
    if cmd == 'cf': cmd_cf(name, int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    elif cmd == 'disk': cmd_disk(name, int(sys.argv[3]) if len(sys.argv) > 3 else 3000)
    elif cmd == 'dedekind': cmd_dedekind(name, int(sys.argv[3]) if len(sys.argv) > 3 else 200)
