"""gdesign.py — design and measure our own weights for Catalan's constant.
F = sum_{t >= t0} (-1)^t R(t),  R(t) = C * (2t + c) * prod_{r in ROOTS} (t + r) / prod_{s in POLES} (t + s)^{mult(s)},
with every pole s a half-integer (s = k + 1/2).  Partial fractions R = sum A_{i,s} (t+s)^{-i} (exact, fmpq series), then
  sum_{t>=t0} (-1)^t (t+k+1/2)^{-i} = 2^i (-1)^k [beta(i) - S_i(t0 + k)],  S_i(L) = sum_{0<=u<L} (-1)^u (2u+1)^{-i}  (L < 0: minus the sum over L<=u<0),
  F = sum_i U_i beta(i) - V,   U_i = 2^i sum (-1)^k A_{i,s},   V = sum 2^i (-1)^k A_{i,s} S_i(t0+k).
Reports: U_i for every i (odd ones must vanish for a pure G form; beta(4), ... must vanish too), the numeric check, log|F|/n,
the true denominators of the primitive form (U_2, V) per prime, and the top-coefficient valuations.
Designs (name, n):
  zudilin       a = 3n+1, h = n+1/2 (three blocks), h4 = n+1   (Zudilin's Apery-like forms; control, = vwp_forms)
  shift         the zeta-template (3;1^4) shifted by 1/2: poles k in [n,2n] order 3, roots at m+1/2 for m in [0,n-1] u [2n+1,3n], c = 3n+1, C = n!
  shift_t0=K    the same with the sum starting at t0 = K (default 0); K may be 'n' (t0 = -n)
usage: python gdesign.py DESIGN n1 n2 [t0=...]
"""
import sys, math
from flint import fmpq, fmpz
from mpmath import mp, mpf, catalan, nsum, inf, log as mlog, fabs
from sympy import primerange, factorint

def design(name, n, t0spec=None):
    half = fmpq(1, 2)
    if name == 'zudilin':
        a = 3*n + 1; h = fmpq(2*n + 1, 2); h4 = n + 1
        roots = [fmpq(i) for i in range(1, h4)] + [fmpq(i) for i in range(a - h4 + 1, a)]
        poles = {}
        for _ in range(3):
            for i in range(0, int(a - 2*h) + 1):
                s = h + i; poles[s] = poles.get(s, 0) + 1
        C = fmpq(math.factorial(int(a - 2*h))**3, math.factorial(h4 - 1)**2)
        c = fmpq(a); t0 = 0
    elif name == 'shift' or name.startswith('shiftL'):
        # half-integer roots [0, L-1] u [L+N+1, 2L+N], poles [L, L+N] order 3, centre 2t + 2L + N + 1; L = lam*n (name 'shiftL0.5'), N = n
        lam = float(name[6:]) if name.startswith('shiftL') else 1.0
        L = int(round(lam*n)); N = n
        roots = [fmpq(m) + half for m in list(range(0, L)) + list(range(L + N + 1, 2*L + N + 1))]
        poles = {fmpq(k) + half: 3 for k in range(L, L + N + 1)}
        C = fmpq(math.factorial(N)); c = fmpq(2*L + N + 1); t0 = 0
    else:
        raise ValueError(name)
    if t0spec is not None:
        t0 = -n if t0spec == 'n' else int(t0spec)
    return dict(roots=roots, poles=poles, C=C, c=c, t0=t0)

def partial_fractions(D):
    roots, poles, C, c = D['roots'], D['poles'], D['C'], D['c']
    A = {}
    for s, S0 in poles.items():
        # Laurent at t = -s: R = (t+s)^{-S} g(t), g(-s + e) expanded to order S-1; the centre factor 2(t + c/2) cancels one order at s = c/2
        centre_here = (c/2 - s == 0)
        S = S0 - 1 if centre_here else S0
        if S <= 0: continue
        ser = [fmpq(0)]*S; ser[0] = fmpq(1)
        const = C*2 if centre_here else C
        def mulpos(ser, alpha):                  # multiply by (alpha + e) = alpha (1 + e/alpha)
            out = [fmpq(0)]*S
            for j in range(S):
                out[j] += ser[j]
                if j + 1 < S: out[j + 1] += ser[j]/alpha
            return out
        def mulneg(ser, alpha, mult):            # multiply by (alpha + e)^{-mult}
            inv = 1/alpha
            pw = [fmpq(1)]
            for _ in range(1, S): pw.append(pw[-1]*(-inv))
            for _ in range(mult):
                out = [fmpq(0)]*S
                for j in range(S):
                    if ser[j] == 0: continue
                    for tt in range(S - j): out[j + tt] += ser[j]*pw[tt]
                ser = out
            return ser
        # centre factor (2t + c) at t = -s + e: (c - 2s) + 2e = 2 * ((c/2 - s) + e)
        if not centre_here:
            alpha = c/2 - s
            const *= 2*alpha; ser = mulpos(ser, alpha)
        for r in roots:
            alpha = r - s
            if alpha == 0: raise ValueError('root at a pole')
            const *= alpha; ser = mulpos(ser, alpha)
        for s2, S2 in poles.items():
            if s2 == s: continue
            alpha = s2 - s
            const /= alpha**S2; ser = mulneg(ser, alpha, S2)
        for i in range(1, S + 1): A[(i, s)] = const*ser[S - i]
    return A

def Ssum(i, L):
    if L >= 0:
        v = fmpq(0)
        for u in range(L): v += fmpq((-1)**u, (2*u + 1)**i)
        return v
    v = fmpq(0)
    for u in range(L, 0): v -= fmpq((-1)**u, (2*u + 1)**i)
    return v

def form(D):
    A = partial_fractions(D)
    t0 = D['t0']
    U = {}; V = fmpq(0); cache = {}
    for (i, s), a in A.items():
        k = int(s - fmpq(1, 2))
        w = a*2**i*(-1)**k
        U[i] = U.get(i, fmpq(0)) + w
        key = (i, t0 + k)
        if key not in cache: cache[key] = Ssum(i, t0 + k)
        V += w*cache[key]
    return A, U, V

def qm(x): return mpf(int(x.p))/mpf(int(x.q))

def R_num(D, t):
    v = qm(D['C'])*(2*t + qm(D['c']))
    for r in D['roots']: v *= (t + qm(r))
    for s, S in D['poles'].items(): v /= (t + qm(s))**S
    return v

def vp(x, p):
    if x == 0: return 10**9
    a, b, v = int(x.p), int(x.q), 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v

def run(name, n1, n2, t0spec=None):
    print('design %s, t0 = %s' % (name, t0spec))
    for n in range(n1, n2 + 1):
        D = design(name, n, t0spec)
        A, U, V = form(D)
        odd_ok = all(U.get(i, 0) == 0 for i in U if i % 2 == 1)
        high = {i: U[i] for i in U if i >= 4 and U[i] != 0}
        U2 = U.get(2, fmpq(0))
        mp.dps = 60 + 4*n
        Fform = qm(U2)*catalan - qm(V)
        t0 = D['t0']
        Fnum = nsum(lambda t: (-1)**int(t)*R_num(D, t), [t0, inf], method='alternating') if n <= 8 else None
        # primitive denominators of (U2, V)
        den = int(U2.q)*int(V.q)//math.gcd(int(U2.q), int(V.q))
        num_g = math.gcd(int(U2.p)*(den//int(U2.q)), int(V.p)*(den//int(V.q)))
        hden = math.log(den)/n if den > 1 else 0.0; hnum = math.log(num_g)/n if num_g > 1 else 0.0
        print('  n=%2d: odd beta coefficients vanish: %s ; beta(>=4) terms: %s ; log|F|/n = %+.4f%s ; log den/n = %.4f, log content/n = %.4f, primitive height/n = %.4f ; NET = %+.4f' % (
            n, odd_ok, {i: 'nonzero' for i in high}, float(mlog(fabs(Fform)))/n,
            (' (numeric check diff %.1e)' % float(fabs(Fnum - Fform))) if Fnum is not None else '',
            hden, hnum, hden - hnum, hden - hnum + float(mlog(fabs(Fform)))/n))
    return D, A, U, V

if __name__ == '__main__':
    name = sys.argv[1]; n1, n2 = int(sys.argv[2]), int(sys.argv[3])
    t0spec = None
    for a in sys.argv[4:]:
        if a.startswith('t0='): t0spec = a[3:]
    run(name, n1, n2, t0spec)
