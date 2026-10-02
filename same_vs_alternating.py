r"""same_vs_alternating.py (2026-10-02): the two ways of zooming in on the start of a ring walk, compared exactly.

N = pi/arcsin(1/sqrt n) = lines per turn of ring n (ring_walk.py).
Alternating sides: regular continued fraction N = a0 + 1/(a1 + 1/(a2 + ...)).
Same side:         minus continued fraction   N = b0 - 1/(b1 - 1/(b2 - ...)),  all b_j >= 2.
1. Conversion rule: b = (a0 + 1), then (a1 - 1) twos, then (a2 + 2), then (a3 - 1) twos, then (a4 + 2), ...
   So the same side needs 1 + a1 + a3 + ... + a_(k-1) stages to cover the alternating stages 0..k.
2. "Same divided by alternating" as stage counts: there is no limit constant.  The ratio grows like (1/2) log2(k).
3. The exact link is a difference.  On the convergent [a0; a1, ..., a_2m] = [[b0, ..., bs]]:
       sum (b_j - 3) = a0 - a1 + a2 - ... + a_2m - 3.
4. That number is a Dedekind sum.  For P/Q = [[b0, ..., bs]] (P lines, Q turns) and Q Q* = 1 mod P:
       12 s(Q, P) = sum (b_j - 3) + (Q + Q*)/P.
   In alternating digits (Barkan, Hickerson), for h/k = [0; a1, ..., ar] with previous denominator q_(r-1):
       12 s(h, k) = sum (-1)^(i+1) a_i + (h + (-1)^(r+1) q_(r-1))/k - 3 [r odd].
5. Weighted walks: weight corner k of the walk by 1/k^2 and add the heights (on the unit circle).
   Same way: sum sin(k theta)/k^2 = Cl2(theta); the heights are odd in k.
   Turning around at every crossing makes the heights even in k, and the sum is pi^2 times an algebraic number
   (psi1(x) + psi1(1 - x) = pi^2/sin^2(pi x)).
   Square ring (n = 2): same way = G, turning around = L(2, chi_8) = pi^2/(8 sqrt 2).
   Hexagon ring (n = 4): same way = Cl2(pi/3) = (3 sqrt 3/4) L(2, chi_-3), turning around = pi^2 (9 + sqrt 3)/108.
Usage: python same_vs_alternating.py"""
import math
from fractions import Fraction as Fr
import mpmath as mp
import numpy as np

mp.mp.dps = 220
def Nring(n): return mp.pi/mp.asin(1/mp.sqrt(n))
def regular(x, k):
    a = []
    for _ in range(k):
        f = mp.floor(x); a.append(int(f)); x = 1/(x - f)
    return a
def minus(x, k):
    b = []
    for _ in range(k):
        c = mp.ceil(x); b.append(int(c)); x = 1/(c - x)
    return b
def minus_q(fr):
    b = []
    while True:
        c = -((-fr.numerator)//fr.denominator); b.append(c); r = c - fr
        if r == 0: return b
        fr = 1/r
def convergents(a):
    p0, q0, p1, q1 = 1, 0, a[0], 1; out = [(p1, q1)]
    for x in a[1:]:
        p0, q0, p1, q1 = p1, q1, x*p1 + p0, x*q1 + q0; out.append((p1, q1))
    return out

print("1. conversion rule")
ok = True
for n in [3] + list(range(5, 60)):
    N = Nring(n); a = regular(N, 12)
    pred = [a[0] + 1]
    for i in range(1, 11, 2):
        pred += [2]*(a[i] - 1) + [a[i + 1] + 2]
    b = minus(N, len(pred))
    ok &= (b == pred)
a3 = regular(Nring(3), 9)
print("   same-side digits = (a0 + 1), (a1 - 1) twos, (a2 + 2), (a3 - 1) twos, ...: holds for rings 3, 5..59: %s" % ok)
print("   ring 3: alternating %s  ->  same side %s" % (a3, minus(Nring(3), 14)))

print("\n2. stage counts, rings 3..1500 (ring 4 excluded): same-side stages needed to cover k alternating stages")
depths = (6, 12, 24, 48, 96)
rat = {k: [] for k in depths}; cen = {k: [] for k in depths}
for n in range(3, 1501):
    if n == 4: continue
    a = regular(Nring(n), 98)
    for k in depths:
        odd = sum(a[i] for i in range(1, k, 2))             # a1 + a3 + ... + a_(k-1): k/2 digits
        rat[k].append((1 + odd)/k); cen[k].append(odd/(k/2) - math.log2(k/2))
print("   alternating stages k:                    " + "  ".join("%7d" % k for k in depths))
print("   median of same/alternating:              " + "  ".join("%7.2f" % float(np.median(rat[k])) for k in depths))
print("   (1/2) log2(k):                           " + "  ".join("%7.2f" % (0.5*math.log2(k)) for k in depths))
print("   mean of same/alternating:                " + "  ".join("%7.2f" % float(np.mean(rat[k])) for k in depths))
print("   median of (digit sum)/(k/2) - log2(k/2): " + "  ".join("%7.2f" % float(np.median(cen[k])) for k in depths))
print("   (the mean is carried by rare giant digits; the digit sum has no finite average)")

print("\n3. the difference, with 3 as the zero of the same-side digits")
ok3 = True; ex = None
for n in [3] + list(range(5, 200)):
    a = regular(Nring(n), 9)                                # a0..a8, ending at an even index
    P, Q = convergents(a)[-1]
    b = minus_q(Fr(P, Q))
    lhs = sum(x - 3 for x in b); rhs = sum((-1)**i*x for i, x in enumerate(a)) - 3
    ok3 &= (lhs == rhs)
    if n == 3: ex = (a, b, lhs, rhs)
print("   sum(b_j - 3) = a0 - a1 + a2 - ... + a8 - 3 on the convergent [a0; a1, ..., a8], rings 3, 5..199: %s" % ok3)
print("   ring 3: alternating %s, same side %s: %d = %d" % ex)

print("\n4. the Dedekind sum")
def dedekind12(h, k):
    """12 s(h, k), exact."""
    tot = 0
    for i in range(1, k):
        tot += (2*i - k)*(2*((h*i) % k) - k)
    return Fr(12*tot, 4*k*k)
okA = True; okB = True; nA = nB = 0; shown = []
for n in (3, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15):
    a = regular(Nring(n), 9)
    for r in (3, 4, 5, 6):                                  # alternating digits: h/k = [0; a1, ..., ar]
        if a[r] == 1: continue
        cv = convergents([0] + a[1:r + 1]); (h, k), (_, qprev) = cv[-1], cv[-2]
        if k > 40000: continue
        rhs = sum((-1)**(i + 1)*a[i] for i in range(1, r + 1)) + Fr(h + (-1)**(r + 1)*qprev, k) - (3 if r % 2 else 0)
        okA &= (dedekind12(h, k) == rhs); nA += 1
    for r in (2, 4, 6):                                     # same-side digits: P/Q = [a0; a1, ..., ar] = [[b0, ..., bs]]
        P, Q = convergents(a[:r + 1])[-1]
        if P > 40000 or Q == 1: continue
        b = minus_q(Fr(P, Q)); lhs = dedekind12(Q, P)
        okB &= (lhs == sum(x - 3 for x in b) + Fr(Q + pow(Q, -1, P), P)); nB += 1
        if len(shown) < 3: shown.append((n, P, Q, b, str(lhs)))
print("   alternating digits: 12 s(h, k) = sum (-1)^(i+1) a_i + (h +- q_(r-1))/k - 3[r odd]: %s (%d cases)" % (okA, nA))
print("   same-side digits:   12 s(Q, P) = sum (b_j - 3) + (Q + Q*)/P: %s (%d cases)" % (okB, nB))
print("   examples (ring, lines P, turns Q, same-side digits, 12 s(Q, P)): %s" % shown)

print("\n5. weighted walks: corner k weighted by 1/k^2, heights added (unit circle)")
mp.mp.dps = 50
def psum(h, M, s=2):
    """sum over k >= 1 of h(k)/k^s for h periodic mod M."""
    return sum(h(r)*mp.zeta(s, mp.mpf(r)/M) for r in range(1, M + 1) if h(r) != 0)/mp.mpf(M)**s
def tri(k, m):
    t = k % (2*m)
    return t if t <= m else 2*m - t
G = mp.catalan
L3 = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
for n, name in ((2, "square"), (4, "hexagon"), (3, "ring 3"), (5, "ring 5")):
    th = 2*mp.asin(1/mp.sqrt(n)); N = 2*mp.pi/th
    m = int(mp.nint(N)) if abs(N - mp.nint(N)) < mp.mpf(10)**(-40) else int(mp.ceil(N))
    same = mp.clsin(2, th)
    turn = psum(lambda k: mp.sin(th*tri(k, m)) if tri(k, m) % m or (n not in (2, 4)) else 0, 2*m)
    closed = mp.pi**2/(4*m*m)*(sum(mp.sin(r*th)/mp.sin(mp.pi*r/(2*m))**2 for r in range(1, m)) + mp.sin(m*th)/2)
    print("   %-8s (%d lines per lap): same way = %s;  turning around = %s = %s pi^2  (sine formula agrees: %s)"
          % (name, m, mp.nstr(same, 20), mp.nstr(turn, 20), mp.nstr(turn/mp.pi**2, 15), abs(turn - closed) < mp.mpf(10)**(-40)))
    print("            same/turning = %s" % mp.nstr(same/turn, 20))
sq_turn = psum(lambda k: mp.sin(mp.pi/2*tri(k, 4)) if tri(k, 4) % 2 else 0, 8)
hx_turn = psum(lambda k: mp.sin(mp.pi/3*tri(k, 6)) if tri(k, 6) % 3 else 0, 12)
print("   square:  same way - G = %s;  turning around - pi^2/(8 sqrt 2) = %s" % (mp.nstr(mp.clsin(2, mp.pi/2) - G, 3), mp.nstr(sq_turn - mp.pi**2/(8*mp.sqrt(2)), 3)))
print("   hexagon: same way - (3 sqrt 3/4) L(2, chi_-3) = %s;  turning around - pi^2 (9 + sqrt 3)/108 = %s"
      % (mp.nstr(mp.clsin(2, mp.pi/3) - 3*mp.sqrt(3)/4*L3, 3), mp.nstr(hx_turn - mp.pi**2*(9 + mp.sqrt(3))/108, 3)))
print("   square: same/turning = 8 sqrt(2) G/pi^2 = %s" % mp.nstr(8*mp.sqrt(2)*G/mp.pi**2, 25))

print("\n   the two sign patterns of the square at other weights k^(-s)  (same way: + - + -;  turning around: + - - +)")
chi4 = {1: 1, 3: -1}; chi8 = {1: 1, 3: -1, 5: -1, 7: 1}
def Lval(chi, f, s):
    if s == 1: return -sum(c*mp.digamma(mp.mpf(a)/f) for a, c in chi.items())/f
    return sum(c*mp.zeta(s, mp.mpf(a)/f) for a, c in chi.items())/mp.mpf(f)**s
known4 = {-1: "0", 0: "1/2", 1: "pi/4", 2: "G", 3: "pi^3/32"}
known8 = {-1: "-1", 0: "0", 1: "log(1 + sqrt 2)/sqrt 2", 2: "pi^2/(8 sqrt 2)", 3: "no closed form known"}
for s in (-1, 0, 1, 2, 3):
    print("   s = %2d: same way %s (%s);  turning around %s (%s)"
          % (s, mp.nstr(mp.chop(Lval(chi4, 4, s)), 15), known4[s], mp.nstr(mp.chop(Lval(chi8, 8, s)), 15), known8[s]))
