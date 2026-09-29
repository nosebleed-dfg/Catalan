r"""Why G's width doubles (2026-09-29): the harmless group is never a lattice, and at the threshold the sign of oo sits at the prime 2.

Exact checks behind the results-log section "The harmless group is never a lattice".
1. Two-parabolic groups <T, L_w>, L_w = (1 0; w 1).  T*L_w^{-1} = (1-w 1; -w 1) has trace 2 - w:
   elliptic for w <= 3, -(parabolic) at w = 4, hyperbolic for w >= 5.  With the Fricke involution the group is
   Hecke's H(sqrt w), a lattice exactly when sqrt w <= 2.  At w = 4, conjugating Sanov's free group
   <(1 2; 0 1), (1 0; 2 1)> by diag(2^-1/2, 2^1/2) gives <T, L_4> = Gamma_0(4)/{+-1}.
2. Where the valence formula puts the weight-1 zero: at the fixed point of the product element.
   Gamma_0(3), chi_-3: the elliptic point e, automorphy factor j = -3e + 1, a primitive cube root of 1 (exact in Q(sqrt-3)).
   Gamma_0(4), chi_-4: sigma^-1 (T L_4^-1) sigma = -T^-1 with sigma = (1 0; 2 1), so f|sigma picks up (-1)^k = chi_-4(-1).
3. Ligozat orders of theta = eta(2t)^5 / (eta(t)^2 eta(4t)^2) on Gamma_0(4), and numerical checks of both forced zeros.
Usage: python lattice_principle.py"""
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp

def mul(A, B):
    a, b, c, d = A; e, f, g, h = B
    return (a*e + b*g, a*f + b*h, c*e + d*g, c*f + d*h)

def inv(A):
    a, b, c, d = A
    assert a*d - b*c == 1
    return (d, -b, -c, a)

T = (1, 1, 0, 1)
def L(w): return (1, 0, w, 1)

print("1. The product of the two harmless parabolics, T * L_w^-1")
for w in range(1, 11):
    M = mul(T, inv(L(w)))
    tr = M[0] + M[3]
    kind = "elliptic" if abs(tr) < 2 else ("-(parabolic)" if abs(tr) == 2 else "hyperbolic")
    cov = "lattice" if w <= 4 else "infinite covolume"
    print("   w = %2d   %-18s trace %3d   %-13s H(sqrt %d): %s" % (w, str(M), tr, kind, w, cov))
# Sanov at w = 4: diag(a, 1/a) with a^2 = 1/2 sends (1 b; 0 1) -> (1 a^2 b; 0 1) and (1 0; c 1) -> (1 0; c/a^2 1)
a2 = Fr(1, 2)
print("   Sanov conjugated (a^2 = 1/2): (1 2; 0 1) -> (1 %s; 0 1),  (1 0; 2 1) -> (1 0; %s 1)" % (a2*2, Fr(2)/a2))

# exact arithmetic in Q(sqrt -3): pairs (x, y) = x + y*s, s^2 = -3
def qmul(u, v): return (u[0]*v[0] - 3*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def qadd(u, v): return (u[0] + v[0], u[1] + v[1])
def qsc(k, u): return (k*u[0], k*u[1])
def qinv(u):
    n = u[0]**2 + 3*u[1]**2
    return (u[0]/n, -u[1]/n)

print("\n2a. Gamma_0(3), weight 1, chi_-3")
g3 = mul(T, inv(L(3)))
e = (Fr(1, 2), Fr(1, 6))                                  # e = 1/2 + sqrt(-3)/6 = (3 + i sqrt 3)/6
num = qadd(qsc(g3[0], e), (Fr(g3[1]), Fr(0)))
den = qadd(qsc(g3[2], e), (Fr(g3[3]), Fr(0)))
print("   product element", g3, "  fixes e = 1/2 + sqrt(-3)/6:", qmul(num, qinv(den)) == e)
j = den                                                   # j(g, e) = c e + d
j3 = qmul(qmul(j, j), j)
print("   automorphy factor j = -3e + 1 = %s %s*sqrt(-3);  j^3 = 1: %s,  j != 1: %s" % (j[0], j[1], j3 == (1, 0), j != (1, 0)))
print("   chi_-3(d) = chi_-3(1) = 1, so f(e) = j^k f(e): weight 1 forces f(e) = 0; weight 3 does not (j^3 = 1)")
print("   valence: weight 1, index [PSL2(Z) : Gamma_0(3)] = 4, total zero order", Fr(1*4, 12))

print("\n2b. Gamma_0(4), weight 1, chi_-4")
g4 = mul(T, inv(L(4)))
sig = (1, 0, 2, 1)
loc = mul(mul(inv(sig), g4), sig)
print("   product element", g4, "; sigma^-1 * it * sigma =", loc, "= -T^-1:", loc == (-1, 1, 0, -1))
print("   chi_-4(d) = chi_-4(1) = 1, so h = f|sigma satisfies h|(-T^-1) = h, i.e. h(tau+1) = (-1)^k h(tau)")
print("   same element as -(product) = %s with d = -1: the sign is chi_-4(-1) = -1 = (-1)^k" % (tuple(-x for x in g4),))
print("   odd k: only exponents in 1/2 + Z at the cusp 1/2.  valence: index 6, total zero order", Fr(1*6, 12))

def ligozat(N, r, dd):
    """order at the cusp c/dd (dd | N) of prod eta(delta tau)^r_delta on Gamma_0(N), in its own local parameter"""
    s = sum(Fr(gcd(dd, de)**2 * rd, de) for de, rd in r.items())
    return Fr(N, 24) * s / (gcd(dd, N // dd) * dd)

rth = {1: -2, 2: 5, 4: -2}
print("\n3. theta = eta(2t)^5/(eta(t)^2 eta(4t)^2) on Gamma_0(4): orders at cusps 0 (d=1), 1/2 (d=2), oo (d=4):",
      [str(ligozat(4, rth, dd)) for dd in (1, 2, 4)], " -> theta^2 has order 1/2 at 1/2 and nowhere else")

mp.mp.dps = 50
q = mp.exp(2j*mp.pi*mp.mpc(mp.mpf(1)/2, mp.sqrt(3)/6))
chi3 = lambda n: (0, 1, -1)[n % 3]
a_e = 1 + 6*mp.fsum(chi3(n)*q**n/(1 - q**n) for n in range(1, 400))
print("   a(q) = sum q^(m^2+mn+n^2) at e: |a(e)| = %s   (|q| = %s)" % (mp.nstr(abs(a_e), 3), mp.nstr(abs(q), 6)))
for y in (mp.mpf('0.05'), mp.mpf('0.02'), mp.mpf('0.01')):
    th = mp.fsum((-1)**n * mp.exp(-2*mp.pi*n*n*y) for n in range(-200, 201))    # theta(1/2 + i y)
    ratio = th * mp.sqrt(2*y) / (2*mp.exp(-mp.pi/(8*y)))                          # |q'|^(1/4) = exp(-pi/(8y)), |q'| = exp(-pi/(2y))
    print("   theta(1/2 + i y) / (2 (2y)^-1/2 |q'|^(1/4)) at y = %s: %s" % (mp.nstr(y, 3), mp.nstr(ratio, 12)))

print("\n4. Idele-class signs of chi_-4 on -1: chi_oo(-1) = sgn(-1) = -1, chi_p(-1) = 1 (p odd, unramified),")
print("   so the product formula forces chi_2(-1) = -1: the half turn at the 2-adic cusp 1/2 carries the sign of oo.")

print("\n5. Smallest allowed width w for Gamma_0(N)/Gamma_1(N) families (w a multiple of N, cond | N, w >= 5 by part 1)")
print("   log R*(w) from uniformize_radius.py; R* decreases in w, so it bounds the two-parabolic families of that conductor")
table = [("zeta(2), Apery (Gamma_1(5))", 5, "5.13", "attained"), ("L(2,chi_-3), case C (CDT)", 3, "2.68", "attained"),
         ("G = L(2,chi_-4), case E", 4, "1.232", "attained"), ("L(2,chi_-7)", 7, "1.73", "bound only"),
         ("L(2,chi_-8)", 8, "1.232", "bound only")]
for name, cond, rad, status in table:
    w = cond
    while w < 5: w += cond
    print("   %-30s conductor %d  ->  w >= %d  (w/cond = %s)   log R* <= %-6s (%s)   vs tau = 1.95" % (name, cond, w, Fr(w, cond), rad, status))
