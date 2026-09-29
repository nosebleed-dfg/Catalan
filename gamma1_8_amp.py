"""Real-G push, step 6b (2026-09-29): high-precision cusp amplitudes on Gamma_1(8) (companion of gamma1_8_kill.py).
For inhomogeneities h in {T^j, T^j/(1-2T), T^j/(1-T)} (j = 1..3): Lam(h) = lim y_n/u_n and the amplitude of the log singularity
of F(h) = y(h) - Lam(h) f at the cusp T = 1/2:  n F_n / 2^n -> A(h)  (polynomial extrapolation in 1/n).
Then: the Q-rank of the amplitudes, and the combinations with Lam in Q + QG and A = 0.
Usage: python gamma1_8_amp.py [N]"""
import sys, io, contextlib, math, json
from fractions import Fraction as Fr
import mpmath as mp
sys.argv = [sys.argv[0], '8'] + sys.argv[1:]
import level12_modular as L
from guess_rec import guess
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1400
mp.mp.dps = 900
P, M = 140, 120

def theta3(m):
    s = [0]*P; n = 0
    while m*n*n < P:
        s[m*n*n] += 1 if n == 0 else 2; n += 1
    return s
th1, th2 = theta3(1), theta3(2)
s = L.mul(th2, L.inv(th1, P), P)
T = [((1 if i == 0 else 0) - s[i])//2 for i in range(P)]
f = L.mul(th1, th1, P)
u = L.expand_in(f, T, M)
with contextlib.redirect_stdout(io.StringIO()):
    R, d, cs = guess(u, R=8, D=2, verbose=False)
ev = lambda c, n: sum(ci*n**e for e, ci in enumerate(c))

def solve(h, upto):
    y = [Fr(0)]
    for n in range(1, upto):
        m = n - R
        s_ = sum(ev(cs[k], m)*(y[m + k] if m + k >= 0 else 0) for k in range(R))
        y.append((Fr(h[n]) - s_)/ev(cs[R], m))
    return y

def series_h(kind, j, upto):
    base = [0]*upto
    if kind == 'poly':
        base[j] = 1; return base
    ratio = {'p2': 2, 'p1': 1}[kind]
    for n in range(j, upto): base[n] = ratio**(n - j)
    return base

Uh = L.run(cs, [u[0]], N)
q = lambda x: mp.mpf(x.numerator)/x.denominator
Uq = [q(x) for x in Uh]
res = {}
for kind in ('poly', 'p2', 'p1'):
    for j in (1, 2, 3):
        y = solve(series_h(kind, j, N), N)
        Lam = q(y[-1])/Uq[-1]
        ns = list(range(N//2 - 60, N//2, 2))
        vals = [n*(q(y[n]) - Lam*Uq[n])/mp.mpf(2)**n for n in ns]
        # polynomial extrapolation in x = 1/n (Neville at x = 0)
        xs = [1/mp.mpf(n) for n in ns]
        pts = list(zip(xs, vals))
        def neville(pts):
            p = [v for _, v in pts]; x = [a for a, _ in pts]
            for k in range(1, len(p)):
                for i in range(len(p) - k):
                    p[i] = ((0 - x[i + k])*p[i] + (x[i] - 0)*p[i + 1])/(x[i] - x[i + k])
            return p[0]
        A1 = neville(pts); A2 = neville(pts[:-4])
        dig = int(-mp.log10(abs(A1 - A2))) if A1 != A2 else 300
        res['%s%d' % (kind, j)] = (Lam, A1, dig)
        print('%s j=%d  Lam = %s   A(1/2) = %s  (~%d digits)' % (kind, j, mp.nstr(Lam, 25), mp.nstr(A1, 40), dig))
        sys.stdout.flush()
# h with a pole at T = 1/2 (kind p2) give a pole there, not a log: their "amplitude" does not converge and is skipped
keys = [k for k in res if res[k][2] >= 30]
mind = min(res[k][2] for k in keys)
print('\nQ-relations among the log amplitudes (%s; PSLQ at %d digits):' % (', '.join(keys), mind - 5))
with mp.workdps(mind - 5):
    vec = [res[k][1] for k in keys]
    for i in range(1, len(keys)):
        rel = mp.pslq([vec[0], vec[i]], maxcoeff=10**6, maxsteps=10**6)
        print('  A(%s) : A(%s)  relation %s' % (keys[0], keys[i], rel))
print('(the amplitude law in terms of L(g,1) and pi: gamma1_8_amplaw.py)')
json.dump({k: [str(v[0]), str(v[1]), v[2]] for k, v in res.items()}, open(r'C:\Users\PC\Desktop\catalan\gamma1_8_amp.json', 'w'), indent=1)
