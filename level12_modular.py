"""Real-G push, step 4 (2026-09-29; the user: "Might need to attack 3G/6 or 6G/12 more modular arithmetic").
Level-N modular families (Beukers' construction), the modular version of the level-12 search.
  t : Hauptmodul of Gamma_0(N), an eta quotient with a simple zero at the cusp oo (t = q + O(q^2), integral);
      every integral Hauptmodul with t = q + O(q^2) is a frame t_a = t/(1 + a t), a an integer.
  f : weight-1 Eisenstein form with character chi_-4 (theta_3(q^m)^2) or chi_-3 (a(q^m), the cubic theta), or a small
      integer combination; f = sum u_n t_a^n with u_n integers.
For each (f, a):
  1. the recurrence of u_n (guess_rec) and its characteristic roots x_i: the singular points of the ODE are 1/x_i;
  2. the Apery second solution v (L v = t: v_0 = 0, v_1 = 1, the same recurrence for n >= 2) and L = lim v_n/u_n;
  3. PSLQ of L against 1, G, L(2,chi_-3), pi^2, pi^2*sqrt3 (the weight-3 Eisenstein values at level 12);
  4. tau = ln den(v_n)/n (measured), r = the radius of sum (v_n - L u_n) t^n (measured from the decay), and the CDT gate
     S = ln(16 r) - tau (necessary; realistically S > 0.8 is needed).
Usage: python level12_modular.py N [amax] [Nn]"""
import sys, math, io, contextlib
from fractions import Fraction as F
import mpmath as mp
import numpy as np
from sympy import Matrix, divisors, Rational
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from guess_rec import guess

NLEV = int(sys.argv[1]) if len(sys.argv) > 1 else 12
AMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 16
NN = int(sys.argv[3]) if len(sys.argv) > 3 else 500        # length of u, v used for the limit
M = 110                                                     # terms used to guess the recurrence
mp.mp.dps = 160
DIVS = divisors(NLEV)

def cusp_order(r, c):
    """Ligozat: order of prod eta(delta tau)^r_delta at a cusp with denominator c (c | N), in the local uniformizer."""
    s = sum(Rational(math.gcd(c, d)**2*r[d], d) for d in DIVS)
    return Rational(NLEV, 24)*s/(math.gcd(c, NLEV//c)*c)

def hauptmoduln():
    """eta quotients with order 1 at oo, order -1 at one cusp c != oo, 0 elsewhere (solve the Ligozat system)."""
    A = Matrix([[Rational(NLEV, 24)*Rational(math.gcd(c, d)**2, d)/(math.gcd(c, NLEV//c)*c) for d in DIVS] for c in DIVS])
    out = []
    for pole in DIVS[:-1]:
        e = Matrix([1 if c == NLEV else (-1 if c == pole else 0) for c in DIVS])
        sol = A.LUsolve(e)
        if all(x.q == 1 for x in sol):
            r = {d: int(x) for d, x in zip(DIVS, sol)}
            if sum(r.values()) == 0: out.append((pole, r))
    return out

def eta_series(r, P):
    """prod_d prod_n (1 - q^{dn})^{r_d} as integer coefficients (the q^{sum d r_d/24} prefactor is handled by the caller)."""
    s = [0]*P; s[0] = 1
    for d, e in r.items():
        for n in range(1, (P - 1)//d + 1):
            m = d*n
            for _ in range(abs(e)):
                if e > 0:
                    for i in range(P - 1, m - 1, -1): s[i] -= s[i - m]
                else:
                    for i in range(m, P): s[i] += s[i - m]
    return s

def mul(a, b, P):
    out = [0]*P
    for i, x in enumerate(a[:P]):
        if x:
            for j in range(0, P - i):
                out[i + j] += x*b[j]
    return out

def inv(a, P):
    """1/a for a[0] = +-1."""
    out = [0]*P; out[0] = a[0]
    for n in range(1, P):
        out[n] = -a[0]*sum(a[k]*out[n - k] for k in range(1, n + 1))
    return out

def chi(D, n):
    if D == -4: return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
    if D == -3: return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)

def eis1(D, m, P):
    """theta_3(q^m)^2 (D = -4) or a(q^m) (D = -3): 1 + c * sum_n (sum_{d|n} chi_D(d)) q^{mn}."""
    c = 4 if D == -4 else 6
    s = [0]*P; s[0] = 1
    for n in range(1, (P - 1)//m + 1):
        s[m*n] = c*sum(chi(D, d) for d in range(1, n + 1) if n % d == 0)
    return s

def expand_in(f, t, P):
    """u_n with f = sum u_n t^n, t = q + O(q^2) integral."""
    u = []; rest = f[:P]; T = [1] + [0]*(P - 1)
    for n in range(P):
        c = rest[n]                       # T = t^n = q^n + ...
        u.append(c)
        if c:
            rest = [x - c*y for x, y in zip(rest, T)]
        T = mul(T, t, P)
    return u

def run(cs, init, upto):
    """forward: sum_k c_k(n) a_{n+k} = 0 in the backward reading a_m for m >= len(init), zero-extended at negative
    indices (valid for the ODE recurrence)."""
    R = len(cs) - 1
    ev = lambda c, n: sum(ci*n**e for e, ci in enumerate(c))
    a = [F(x) for x in init]
    for m in range(len(init), upto):
        n = m - R
        lead = ev(cs[R], n)
        if lead == 0: return None
        s = sum(ev(cs[k], n)*(a[n + k] if n + k >= 0 else 0) for k in range(R))
        a.append(-s/lead)
    return a

L3 = mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])
BASIS = [('1', mp.mpf(1)), ('G', +mp.catalan), ('L3', L3), ('pi^2', mp.pi**2), ('pi^2*sqrt3', mp.pi**2*mp.sqrt(3))]

def analyse(f, t, label):
    u = expand_in(f, t, M)
    # the second-order ODE gives quadratic coefficients (degree 2 in n); a shorter recurrence of higher degree belongs
    # to a higher-order operator with spurious extra solutions, so degree 2 is tried first
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sol = guess(u, R=8, D=2, verbose=False)
        if sol is None:
            sol = guess(u, R=8, D=3, verbose=False)
    if sol is None:
        return None
    R, d, cs = sol
    deg = max(len(c) - 1 for c in cs)
    lc = [(c[deg] if len(c) > deg else 0) for c in cs]      # sum_k lc_k x^k: growth rates x of a_n
    roots = np.roots(lc[::-1]) if any(lc[1:]) else np.array([])
    roots = sorted(roots, key=lambda z: -abs(z))
    U = run(cs, [u[0]], NN)                                  # validates zero-extension: must reproduce u
    if U is None or any(U[i] != u[i] for i in range(M)):
        return ('bad-ext', R, d, roots)
    V = run(cs, [0, 1], NN)
    if V is None:
        return ('bad-v', R, d, roots)
    info = {'R': R, 'd': d, 'roots': roots}
    if len(roots) < 2 or abs(abs(roots[0]) - abs(roots[1])) < 1e-9*abs(roots[0]):
        info['L'] = None
        return info
    mpq = lambda x: mp.mpf(x.numerator)/x.denominator
    Lv = mpq(V[-1])/mpq(U[-1]); Lv2 = mpq(V[-31])/mpq(U[-31])
    err = abs(Lv - Lv2)
    info['digits'] = int(-mp.log10(err)) if err > 0 else 150
    info['L'] = Lv
    # decay of the linear forms F_n = v_n - L u_n: use n up to where L is still accurate
    # F_n is only reliable while (x_2/x_1)^n stays above the accuracy of L
    ratio = abs(roots[0])/abs(roots[1])
    nmax = max(40, min(NN - 1, int(0.6*NN), int((info['digits'] - 20)/math.log10(ratio))))
    Fs = [abs(mpq(V[n]) - Lv*mpq(U[n])) for n in range(nmax)]
    lo, hi = nmax//2, nmax - 1
    if Fs[hi] > 0 and Fs[lo] > 0:
        info['rate'] = float((mp.log(Fs[hi]) - mp.log(Fs[lo]))/(hi - lo))   # ln |x| of the linear forms
    info['tau'] = math.log(V[-1].denominator)/(NN - 1)
    if info['digits'] >= 60:
        with mp.workdps(55):
            rel = mp.pslq([Lv] + [b for _, b in BASIS], maxcoeff=10**5, maxsteps=10**6)
        if rel and rel[0]:
            resid = abs(sum(k*b for k, b in zip(rel, [Lv] + [b for _, b in BASIS])))
            if resid < mp.mpf(10)**(-(info['digits'] - 5)):
                info['rel'] = rel
    return info

def fmt_rel(rel):
    return '%d*L = %s' % (rel[0], ' + '.join('%d*%s' % (-k, nm) for k, (nm, _) in zip(rel[1:], BASIS) if k))

if __name__ == '__main__':
    P = M + 2
    print('Gamma_0(%d): Hauptmoduln (eta quotients, zero at oo)' % NLEV)
    hs = hauptmoduln()
    for pole, r in hs:
        print('  pole at cusp 1/%d: r = %s' % (pole, r))
    if not hs:
        print('  none'); sys.exit()
    pole, r = hs[0]
    A0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0       # centre of the frame scan
    base = eta_series(r, P)
    t0 = [0] + base[:P - 1]                                  # q^1 * prod
    forms = []
    for D in (-4, -3):
        cond = 4 if D == -4 else 3
        if NLEV % cond: continue
        ms = [m for m in DIVS if (NLEV//cond) % m == 0]
        for m in ms:
            forms.append(('%s(q^%d)' % ('th3^2' if D == -4 else 'a', m), D, eis1(D, m, P)))
    print('weight-1 forms: %s' % ', '.join(nm for nm, _, _ in forms))
    rows = []
    for a in range(A0 - AMAX, A0 + AMAX + 1):
        # t_a = t0/(1 + a t0)
        den = [1] + [a*x for x in t0[1:]]
        ta = mul(t0, inv(den, P), P)
        for nm, D, f in forms:
            info = analyse(f, ta, nm)
            if info is None or isinstance(info, tuple):
                print('a=%+3d %-10s %s' % (a, nm, 'no recurrence' if info is None else info[0])); continue
            rts = ' '.join('%.4g' % abs(z) for z in info['roots'])
            if info.get('L') is None:
                print('a=%+3d %-10s R=%d d=%d |x| = %s  (no dominant root)' % (a, nm, info['R'], info['d'], rts)); continue
            x2 = math.exp(info.get('rate', float('nan')))
            r_meas = 1/x2
            S = math.log(16*r_meas) - info['tau']
            relS = fmt_rel(info['rel']) if 'rel' in info else '(no small relation; %d digits)' % info['digits']
            print('a=%+3d %-10s R=%d d=%d |x| = %s  |F_n|^(1/n) -> %.4f  tau=%.3f  S=%+.3f  %s' % (
                a, nm, info['R'], info['d'], rts, x2, info['tau'], S, relS))
            sys.stdout.flush()
            rows.append((S, a, nm, info))
    print('\n=== families whose limit involves G, best gate first ===')
    for S, a, nm, info in sorted(rows, key=lambda z: -z[0]):
        if 'rel' in info and info['rel'][2] != 0:
            print('S=%+.3f  a=%+3d  %-10s tau=%.3f  %s' % (S, a, nm, info['tau'], fmt_rel(info['rel'])))
