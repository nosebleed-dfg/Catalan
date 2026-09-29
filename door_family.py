r"""The non-congruence door, part 3 (2026-09-29): the index-9 group with harmless pair (oo: width 1, 0: width 6), built exactly.

Group: the index-3 subgroup of Gamma_0(2) cut out by the non-Galois cubic cover  y = t / (1 - 256 t/27)^3,
y = (eta(2tau)/eta(tau))^24 = q prod (1+q^n)^24 (Hauptmodul of Gamma_0(2), j = (1+256y)^3/y).
Ramification of the cover: y = 0 <- {t = 0 (simple, width 1), t = oo (double, width 2)}; y = oo <- {t0 = 27/256 (triple,
width 6: the cusp 0)}; y = -1/64 (the elliptic point of Gamma_0(2)) <- {b = 27/64 (simple: elliptic of order 2),
beta = -27/512 (double: regular)}.  Monodromy of order 54, non-congruence (door_groups.py).
Weight 1: f = E4^{1/4} (A(0)/A(t))^{1/4}, A(t) = 256 t + (1 - 256t/27)^3 (the points over j = 0); f(0) = 1, f nonzero at
the cusp 0 and at b, order 3/4 at the width-2 cusp (multiplier -i there).
Family: W_k = (Dt)^2 t^(k-1) / (p(t) f), p = (1 - t/t0)(1 - t/b), g_k = D^{-2} W_k, y_k = f g_k;
u_n = [t^n] f, v_{k,n} = [t^n] y_k, Lambda_k = lim v_{k,n}/u_n (rate (t0/b)^n = 4^{-n}).  Harmless pair (oo, 0): log R* >= 2.68.
Part (i): the limits, tested against a declared basis.  Part (ii): exact 3-adic and 2-adic denominators of t(q) and of v_n.
Usage: python door_family.py [N]"""
import sys
import mpmath as mp
from fractions import Fraction as Fr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
mp.mp.dps = 260
P = N + 1
def mul(a, b):
    r = [0]*P
    for i, x in enumerate(a):
        if x:
            for j in range(P - i): r[i + j] += x*b[j]
    return r
def powser(a, e):          # a[0] = 1, a^e for rational e via log/exp recursion: (a^e)' a = e a' a^e
    r = [0]*P; r[0] = type(a[0])(1) if not isinstance(a[0], int) else 1
    r[0] = a[0]**0
    for n in range(1, P):
        s = 0
        for k in range(1, n + 1):
            s += (e*k - (n - k))*a[k]*r[n - k]
        r[n] = s/(n*a[0])
    return r
def inv(a):
    r = [0]*P; r[0] = 1/a[0]
    for n in range(1, P): r[n] = -sum(a[k]*r[n - k] for k in range(1, n + 1))/a[0]
    return r

def build(field):
    """field = Fr (exact) or mp.mpf.  Returns q-series of t, f, Dt and the t-power table."""
    one = field(1)
    # y = q prod (1+q^n)^24
    pr = [one] + [field(0)]*(P - 1)
    for n in range(1, P):
        fac = [field(0)]*P; fac[0] = one
        # (1 + q^n)^24 via binomial
        from math import comb
        for k in range(1, 25):
            if n*k < P: fac[n*k] = field(comb(24, k))
        pr = mul(pr, fac)
    y = [field(0)] + pr[:P - 1]
    # t = y (1 - c t)^3, c = 256/27, fixed point
    c = field(256)/field(27)
    t = y[:]
    for it in range(P):
        u = [(-c)*x for x in t]; u[0] += one
        u3 = mul(mul(u, u), u)
        tn = mul(y, u3)
        if all(tn[i] == t[i] for i in range(P)) if field is Fr else max(abs(tn[i] - t[i]) for i in range(P)) < mp.mpf(10)**(-mp.mp.dps + 20):
            t = tn; break
        t = tn
    # E4 = 1 + 240 sum sigma3(n) q^n
    E4 = [field(0)]*P; E4[0] = one
    for n in range(1, P):
        E4[n] = field(240*sum(d**3 for d in range(1, n + 1) if n % d == 0))
    E4q = powser(E4, Fr(1, 4) if field is Fr else mp.mpf(1)/4)
    # A(t) = 256 t + (1 - c t)^3, as a q-series
    u = [(-c)*x for x in t]; u[0] += one
    A = mul(mul(u, u), u)
    A = [A[i] + field(256)*t[i] for i in range(P)]
    Am = powser(A, Fr(-1, 4) if field is Fr else -mp.mpf(1)/4)
    f = mul(E4q, Am)
    Dt = [field(n)*t[n] for n in range(P)]
    tp = [[one] + [field(0)]*(P - 1)]
    for k in range(1, P): tp.append(mul(tp[-1], t))
    return t, f, Dt, tp, c

def expand_in(X, tp):
    res = list(X); out = [0]*P
    for k in range(P):
        ck = res[k]; out[k] = ck
        if ck:
            row = tp[k]
            for j in range(k, P): res[j] -= ck*row[j]
    return out

print("(i) limits (mpmath, %d digits, N = %d)" % (mp.mp.dps, N))
t, f, Dt, tp, c = build(mp.mpf)
t0 = mp.mpf(27)/256; b = mp.mpf(27)/64
p = [(-1/t0 - 1/b)*x for x in t]; p[0] += 1
tt = mul(t, t); p = [p[i] + tt[i]/(t0*b) for i in range(P)]
fu = expand_in(f, tp)
Dt2 = mul(Dt, Dt); base = mul(Dt2, inv(mul(p, f)))
LAM = []
for k in (1, 2, 3):
    W = base if k == 1 else mul(base, tp[k - 1])
    assert abs(W[0]) < mp.mpf(10)**-200
    g = [mp.mpf(0)] + [W[n]/mp.mpf(n)**2 for n in range(1, P)]
    yk = expand_in(mul(f, g), tp)
    lam = yk[N]/fu[N]; lam2 = yk[N - 10]/fu[N - 10]
    LAM.append(lam)
    print("   k = %d: Lambda = %s   (change over 10 steps %s)" % (k, mp.nstr(lam, 40), mp.nstr(abs(lam - lam2), 3)))
G = mp.catalan; pi = mp.pi
digits = int(-mp.log10(abs(LAM[0] - (expand_in(mul(f, [mp.mpf(0)] + [base[n]/mp.mpf(n)**2 for n in range(1, P)]), tp)[N - 10]/fu[N - 10]))))
print("   working precision of the limits: about %d digits" % digits)
mp.mp.dps = max(30, digits - 8)
rel = mp.pslq(LAM + [mp.mpf(1)], maxcoeff=10**12, maxsteps=10**6)
print("   relation among Lambda_1, Lambda_2, Lambda_3, 1: %s" % (rel,))
Om = mp.gamma(mp.mpf(1)/4)**4/pi**2
tests = [("1, G", [1, G]), ("1, pi^2", [1, pi**2]), ("1, G, pi^2", [1, G, pi**2]), ("1, Om", [1, Om]), ("1, pi, pi^2", [1, pi, pi**2]),
         ("1, G, pi^2, pi*log2", [1, G, pi**2, pi*mp.log(2)]), ("1, pi*log2, log2^2", [1, pi*mp.log(2), mp.log(2)**2]),
         ("1, L(2,chi_-3)", [1, mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])]),
         ("1, L(2,chi_-8)", [1, mp.nsum(lambda k: 1/(8*k+1)**2 + 1/(8*k+3)**2 - 1/(8*k+5)**2 - 1/(8*k+7)**2, [0, mp.inf])]),
         ("1, G, pi^2, Om, pi*sqrt3", [1, G, pi**2, Om, pi*mp.sqrt(3)])]
for k, lam in enumerate(LAM[:2], 1):
    for name, vals in tests:
        rel = mp.pslq([lam] + [mp.mpf(v) for v in vals], maxcoeff=10**10, maxsteps=10**6)
        print("   Lambda_%d over {%s}: %s" % (k, name, rel))
