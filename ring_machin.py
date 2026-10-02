r"""ring_machin.py (2026-10-02): which same-way ring sums add up to a rational multiple of G?  ("Machin formulas of weight two".)

W(n) = sum_k sin(k theta_n)/k^2 = Cl2(theta_n) = D(z_n): the same-way walk around ring n with corner k weighted by 1/k^2
(D = Bloch-Wigner dilogarithm, z_n = e^(i theta_n) the step of ring n; same_vs_alternating.py part 5).
For n = m^2 + 1 the step is z = (m + i)/(m - i) in Q(i) and 1 - z = -2i/(m - i), so z and 1 - z only involve the primes
of 2 and of m^2 + 1.  If xi = sum n_m [z_m] has sum n_m z_m ^ (1 - z_m) = 0 in Lambda^2(Q(i)^x) (x) Q, then xi lies in the
Bloch group of Q(i), which has rank 1 and is spanned by [i] (Borel, Bloch, Suslin), with D(i) = G.  So D(xi) = r G, r rational.
1. S = {2, 5}: the rings m = 2, 3, 7 (n = 5, 10, 50).  One relation: 2 W(5) + 6 W(10) + W(50) = 9 G.
2. Proof of it from four five-term relations D(x) + D(y) + D((1-x)/(1-xy)) + D(1-xy) + D((1-y)/(1-xy)) = 0,
   reduced by D(x) = D(1/(1-x)) = D(1-1/x) = -D(1/x) = -D(1-x) = -D(x/(x-1)) = -D(conj x); exact linear algebra.
3. Larger prime sets: Stormer rings (m^2 + 1 smooth), the lattice of relations, and the multiples of G they give.
4. All points i^j u/conj(u) of the unit circle of Q(i) for S = {2, 5} (the 3-4-5 angles): Lewin's identity appears.
5. Machin's pair 5, 239 (S = {2, 13}): no relation between the two rings alone; with a quarter-turn twist there is one.
6. W(m^2 + 1) = sum_k 2^k Im((1 + m i)^k)/((m^2 + 1)^k k^2) + theta_m log(sqrt(m^2 + 1)/2): the first relation in base 5.
Status: relation 1 is proved (part 2).  For the others the multiple of G is rational by the theorem and is read off
from 48 digits.  Ring lists: m <= 10^8 searched (for primes below 100 the count agrees with Luca's 156 values).
Usage: python ring_machin.py   (about 75 seconds)"""
import math
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
G = mp.catalan
EPS = mp.mpf(10)**(-48)

# ---------- Gaussian primes and the wedge x ^ (1 - x) ----------
def two_squares(p):
    y = 1
    while True:
        x = math.isqrt(p - y*y)
        if x*x + y*y == p: return max(x, y), min(x, y)
        y += 1
def gfac(a, b):
    """exponents of the Gaussian primes in a + bi, units dropped: (2,0) = 1+i, (p,1) = x+iy (x > y > 0), (p,-1) its conjugate."""
    out = {}
    for p, e in sp.factorint(a*a + b*b).items():
        if p == 2: out[(2, 0)] = e
        elif p % 4 == 3: out[(p, 0)] = e//2
        else:
            x, y = two_squares(p); k = 0; c, d = a, b
            while True:
                re, im = c*x + d*y, d*x - c*y
                if re % p or im % p: break
                c, d = re//p, im//p; k += 1
            if k: out[(p, 1)] = k
            if e - k: out[(p, -1)] = e - k
    return out
def vsub(u, v):
    out = dict(u)
    for k, e in v.items(): out[k] = out.get(k, 0) - e
    return {k: e for k, e in out.items() if e}
def wedge(v, w):
    keys = sorted(set(v) | set(w)); out = {}
    for i, p in enumerate(keys):
        for q in keys[i + 1:]:
            c = v.get(p, 0)*w.get(q, 0) - v.get(q, 0)*w.get(p, 0)
            if c: out[(p, q)] = c
    return out
def beta(al, be):
    """x = al/be (Gaussian integers): the vector x ^ (1 - x), and the primes of 1 - x."""
    fb = gfac(*be); one_minus = vsub(gfac(be[0] - al[0], be[1] - al[1]), fb)
    return wedge(vsub(gfac(*al), fb), one_minus), one_minus
def int_kernel(rows):
    """integer relations among the rows (dicts): a reduced basis of the full lattice, and the rank of the rows."""
    cols = sorted(set().union(*[set(r) for r in rows]))
    M = [[r.get(c, 0) for c in cols] for r in rows]; R = len(rows)
    try:
        from flint import fmpz_mat
        rank = fmpz_mat(M).rank()
    except Exception:
        rank = sp.Matrix(M).rank()
    try:
        from flint import fmpz_mat
        A = fmpz_mat([[int(i == j) for j in range(R)] + [10**40*x for x in M[i]] for i in range(R)]).lll()
        ker = [[int(A[i, j]) for j in range(R)] for i in range(R) if all(int(A[i, j]) == 0 for j in range(R, R + len(cols)))]
    except Exception:
        ker = []
        for v in sp.Matrix(M).T.nullspace():
            den = math.lcm(*[sp.Rational(x).q for x in v], 1); w = [int(x*den) for x in v]
            g = math.gcd(*w, 0); ker.append([x//g for x in w])
    assert len(ker) == R - rank
    return ker, rank
def W(m): return mp.clsin(2, 2*mp.atan(mp.mpf(1)/m))
def as_rational(val):
    fr = Fr(mp.nstr(val, 45)).limit_denominator(10**6)
    assert abs(val - mp.mpf(fr.numerator)/fr.denominator) < EPS, "not rational to 48 digits"
    return fr
def show(ms, vec):
    return " ".join("%+d W(%d)" % (c, m*m + 1) for m, c in zip(ms, vec) if c)
def ring_row(m): return beta((m, 1), (m, -1))[0]

print("1. S = {2, 5}: rings m = 2, 3, 7")
for m in (2, 3, 7):
    print("   m = %d, ring %d: step (m+i)/(m-i), m + i = %s;  z ^ (1-z) = %s" % (m, m*m + 1, gfac(m, 1), ring_row(m)))
ker, rank = int_kernel([ring_row(m) for m in (2, 3, 7)])
for v in ker:
    val = sum(c*W(m) for m, c in zip((2, 3, 7), v))/G
    if val < 0: v = [-c for c in v]; val = -val
    print("   relation: %s = %s G   (ratio to G: %s)" % (show((2, 3, 7), v), as_rational(val), mp.nstr(val, 50)))
th5, th10, th50 = [2*mp.atan(mp.mpf(1)/m) for m in (2, 3, 7)]
print("   angles: theta_5 + theta_10 - pi/2 = %s;  theta_5 - theta_10 - theta_50 = %s" % (mp.nstr(th5 + th10 - mp.pi/2, 3), mp.nstr(th5 - th10 - th50, 3)))

print("\n2. proof from five-term relations")
def qmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def qinv(x):
    n = x[0]*x[0] + x[1]*x[1]; return (x[0]/n, -x[1]/n)
def qsub(x, y): return (x[0] - y[0], x[1] - y[1])
def q(a, b): return (Fr(a), Fr(b))
ONE = q(1, 0)
def orbit(x):
    inv = qinv(x); om = qsub(ONE, x)
    six = [(x, 1), (qinv(om), 1), (qsub(ONE, inv), 1), (inv, -1), (om, -1), (qmul(x, qinv(qsub(x, ONE))), -1)]
    return six + [((y[0], -y[1]), -s) for y, s in six]
def canon(x):
    """(rep, s) with D(x) = s D(rep); real x has D = 0."""
    if x[1] == 0: return None, 0
    return min(orbit(x), key=lambda t: t[0])
def five(x, y):
    d = qinv(qsub(ONE, qmul(x, y)))
    return [x, y, qmul(qsub(ONE, x), d), qsub(ONE, qmul(x, y)), qmul(qsub(ONE, y), d)]
def relvec(terms):
    v = {}
    for t, c in terms:
        rep, s = canon(t)
        if s: v[rep] = v.get(rep, 0) + c*s
    return {k: c for k, c in v.items() if c}
def Dnum(x):
    z = mp.mpc(mp.mpf(x[0].numerator)/x[0].denominator, mp.mpf(x[1].numerator)/x[1].denominator)
    if abs(z) > 1 + EPS: return -Dnum(qinv(x))
    if abs(abs(z) - 1) < EPS: return mp.clsin(2, mp.arg(z))
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z)*mp.log(abs(z))
I_, z5, z10, z50, u = q(0, 1), (Fr(3, 5), Fr(4, 5)), (Fr(4, 5), Fr(3, 5)), (Fr(24, 25), Fr(7, 25)), q(2, 1)
names = {}
for lab, x in (("G", I_), ("A", z5), ("B", z10), ("C", z50), ("d", u), ("e", q(-2, -1)), ("f", q(0, 2))):
    rep, s = canon(x); names[rep] = (lab, s)
def pretty(v):
    out = []
    for rep, c in v.items():
        lab, s = names.get(rep, ("D(%s%+si)" % rep, 1)); out.append("%+d %s" % (c*s, lab))
    return " ".join(sorted(out, key=lambda t: t.split()[1])) + " = 0"
print("   names: G = D(i), A = W(5), B = W(10), C = W(50), d = D(2+i), e = D(-2-i), f = D(2i)")
pairs = [("R1", I_, qmul(q(0, -1), z5)), ("R2", q(1, 1), q(1, 1)), ("R3", I_, q(-2, -1)), ("R4", qmul(q(0, -1), z5), z5)]
rels = []
for lab, x, y in pairs:
    terms = five(x, y); num = sum(Dnum(t) for t in terms if t[1] != 0)
    v = relvec([(t, 1) for t in terms]); rels.append(v)
    print("   %s: (x, y) = (%s%+si, %s%+si): five-term sum = %s;  reduced: %s" % (lab, x[0], x[1], y[0], y[1], mp.nstr(num, 3), pretty(v)))
target = relvec([(z5, 2), (z10, 6), (z50, 1), (I_, -9)])
cols = sorted(set().union(*[set(v) for v in rels + [target]]))
M = sp.Matrix([[v.get(c, 0) for c in cols] for v in rels]); Mt = sp.Matrix([[v.get(c, 0) for c in cols] for v in rels + [target]])
sol = sp.Matrix([[v.get(c, 0) for c in cols] for v in rels]).T.solve_least_squares(sp.Matrix([target.get(c, 0) for c in cols])) if M.rank() == Mt.rank() else None
print("   target 2A + 6B + C - 9G: %s" % pretty(target))
print("   rank of the four relations = %d, with the target added = %d;  target = %s . (R1, R2, R3, R4)" % (M.rank(), Mt.rank(), None if sol is None else list(sol)))

print("\n3. larger prime sets: rings with m^2 + 1 smooth (m <= 10^8 searched)")
ALL = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101]
BOUND, CH = 10**8, 10**7
smooth = []
for lo in range(1, BOUND + 1, CH):
    mm = np.arange(lo, min(lo + CH, BOUND + 1), dtype=np.int64); vv = mm*mm + 1
    for p in [2] + ALL:
        idx = np.nonzero(vv % p == 0)[0]
        while idx.size:
            vv[idx] //= p; idx = idx[vv[idx] % p == 0]
    smooth += [int(x) for x in mm[vv == 1]]
del mm, vv
top_prime = {m: max(sp.factorint(m*m + 1)) for m in smooth}
def signed(v, fr):
    if fr < 0 or (fr == 0 and next(c for c in v if c) < 0): return [-c for c in v], -fr
    return v, fr
print("   (z ^ (1-z) of a ring lies in a space of dimension s^2 + s for s primes p = 1 mod 4: (x_p - y_p)^(1+i), x_p^y_p,")
print("    and for each pair x_p^y_q + x_q^y_p, x_p^x_q - y_p^y_q.  Luca: 156 values of m have m^2 + 1 with all primes < 100.)")
prev = 0
for top in ALL:
    s = ALL.index(top) + 1
    ms = [m for m in smooth if m > 1 and top_prime[m] <= top]
    ker, rank = int_kernel([ring_row(m) for m in ms])
    out = [signed(v, as_rational(sum(c*W(m) for m, c in zip(ms, v))/G)) for v in ker]
    rs = sorted(fr for v, fr in out)
    g = math.gcd(*[int(fr) for fr in rs], 0) if all(fr.denominator == 1 for fr in rs) else None
    print("   primes 2, 5..%d: %d rings, largest m = %d;  rank %d of at most %d;  %d relations;  gcd of the multiples of G: %s"
          % (top, len(ms), ms[-1], rank, s*s + s, len(ker), g))
    if top <= 41: print("      m = %s" % ms)
    if top <= 17:
        for v, fr in sorted(out, key=lambda t: t[1]):
            print("      %s = %s G" % (show(ms, v), fr))
    elif len(ker) > prev:
        new = [(v, fr) for v, fr in out if any(c and top_prime[m] == top for m, c in zip(ms, v))]
        v, fr = min(new, key=lambda t: max(abs(c) for c in t[0]))
        print("      %d new; the tamest new reduced relation has %d terms, largest coefficient %d, and gives %s G"
              % (len(ker) - prev, sum(1 for c in v if c), max(abs(c) for c in v), fr))
    prev = len(ker)
    if top == 17:                                           # a readable basis of the three relations
        basis = [{2: 2, 3: 6, 7: 1},
                 {3: 3, 4: 1, 5: 1, 7: -2, 8: -1, 13: -2, 18: -1, 38: -1, 47: -2, 57: 1, 268: 1},
                 {2: 36, 3: 60, 5: 24, 8: 12, 57: 4, 239: -1}]
        V = sp.Matrix([[b.get(m, 0) for m in ms] for b in basis]); K = sp.Matrix(ker)
        X = V*K.T*(K*K.T).inv()
        assert X*K == V and all(x.is_integer for x in X) and abs(X.det()) == 1
        print("      a basis of the lattice (checked: unimodular change from the reduced basis):")
        for b in basis:
            val = sum(c*W(m) for m, c in b.items())/G
            print("      (%s) %s = %s G   (ratio %s)" % ("I" * (basis.index(b) + 1), show(sorted(b), [b[m] for m in sorted(b)]), as_rational(val), mp.nstr(val, 45)))

print("\n4. all unit-circle points for S = {2, 5} (the 3-4-5 angles)")
def circle_points(p, E):
    x, y = two_squares(p); allowed = {(2, 0), (p, 1), (p, -1)}; pts = []
    a, b = 1, 0
    for e in range(1, E + 1):
        a, b = a*x - b*y, a*y + b*x                         # u = (x + iy)^e
        al = (a, b)
        for j in range(4):
            w, om = beta(al, (a, -b))
            if set(om) <= allowed:
                th = mp.arg(mp.mpc(al[0], al[1])/mp.mpc(a, -b))
                pts.append((w, th))
            al = (-al[1], al[0])                            # multiply by i
    return pts
def angle_name(th):
    """D of the point as (sign, name)."""
    t = abs(th); s = -1 if th < 0 else 1
    for m in (2, 3, 5, 7, 239):
        if abs(t - 2*mp.atan(mp.mpf(1)/m)) < EPS: return s, "W(%d)" % (m*m + 1)
    return s, "Cl2(%s deg)" % mp.nstr(t*180/mp.pi, 7)
for p, E in ((5, 40), (13, 40)):
    pts = circle_points(p, E); nm = [angle_name(th) for w, th in pts]
    if p == 13: print("\n5. Machin's pair: unit-circle points for S = {2, 13}")
    print("   u = %d + %di, powers up to %d: %d points with 1 - x smooth: %s" % (two_squares(p) + (E, len(pts), [n for s, n in nm])))
    ker, rank = int_kernel([w for w, th in pts])
    for v in ker:
        v, fr = signed(v, as_rational(sum(c*mp.clsin(2, th) for c, (w, th) in zip(v, pts))/G))
        print("      " + " ".join("%+d %s" % (c*s, n) for c, (s, n) in zip(v, nm) if c) + " = %s G" % fr)
    if p == 13:
        ker2, rank2 = int_kernel([ring_row(5), ring_row(239)])
        print("   the two rings 26 and 57122 alone: rank %d, %d relations" % (rank2, len(ker2)))
Ti2 = lambda x: mp.im(mp.polylog(2, mp.mpc(0, x)))
print("\n   Lewin's inverse-tangent identity 6 Ti2(1) - 4 Ti2(1/2) - 2 Ti2(1/3) - Ti2(3/4) - pi log 2 = %s"
      % mp.nstr(6*Ti2(1) - 4*Ti2(mp.mpf(1)/2) - 2*Ti2(mp.mpf(1)/3) - Ti2(mp.mpf(3)/4) - mp.pi*mp.log(2), 3))
print("   (in Clausen form it is the relation 3 W(5) + 2 W(10) + Cl2(pi - theta_5) = 6 G of part 4)")

print("\n6. the first relation as a series in powers of 1/5")
def im_powers(m, K):
    """Im((1 + m i)^k) for k = 1..K, exact."""
    a, b, out = 1, 0, []
    for _ in range(K):
        a, b = a - m*b, m*a + b; out.append(b)
    return out
def series(m, K):
    return sum(mp.mpf(2**k*b)/(mp.mpf(m*m + 1)**k*k*k) for k, b in enumerate(im_powers(m, K), 1))
for m, K in ((7, 300), (239, 60)):
    th = 2*mp.atan(mp.mpf(1)/m)
    print("   W(%d) - [sum_k 2^k Im((1+%di)^k)/(%d^k k^2) + theta log(sqrt(%d)/2)] = %s"
          % (m*m + 1, m, m*m + 1, m*m + 1, mp.nstr(W(m) - series(m, K) - th*mp.log(mp.sqrt(m*m + 1)/2), 3)))
K = 1400; I2, I3, I7 = im_powers(2, K), im_powers(3, K), im_powers(7, K)
S5 = sum((mp.mpf(2**(k + 1)*I2[k - 1] + 6*I3[k - 1])/mp.mpf(5)**k + mp.mpf(I7[k - 1])/mp.mpf(25)**k)/(k*k) for k in range(1, K + 1))
print("   9 G - pi log 5 + (5 pi/4) log 2 - sum_k [(2^(k+1) Im(1+2i)^k + 6 Im(1+3i)^k)/5^k + Im(1+7i)^k/25^k]/k^2 = %s"
      % mp.nstr(9*G - mp.pi*mp.log(5) + 5*mp.pi/4*mp.log(2) - S5, 3))
