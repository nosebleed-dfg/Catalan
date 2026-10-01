r"""ramanujan_dodecahedron.py (2026-10-01, real date): how Ramanujan summation sits on the dodecahedron, checked exactly.

The dodecahedron is the modular curve X(5) = H/Gamma(5) with its cusps filled: PSL2(F_5) = A_5 (60 rotations); the 12 cusps
(width 5) are the face centres, the 20 points over j = 0 the vertices, the 30 points over j = 1728 the edge midpoints.
1. Euler's formula through zeta(-1):  chi(X(N)) = |PSL2(Z/N)| (1/N - 1/6),  1/6 = -2 zeta(-1).
2. The fractional powers of q are Ramanujan sums over residue classes mod 5:
     sum_{n = a mod 5} n := 5 zeta(-1, a/5) = -(5/2) B_2(a/5);   G: (1/2)(S_1 + S_4) = -1/60,  H: (1/2)(S_2 + S_3) = 11/60,
     and their difference 1/5 = -(1/2) L(-1, chi_5) is the exponent of the Rogers-Ramanujan continued fraction r.
3. With those exponents (q^-1/60 G, q^11/60 H) transforms linearly under tau -> -1/tau, and r = q^1/5 H/G is Klein's
   icosahedral coordinate:  r(-1/tau) = (1 - phi r)/(phi + r)  (the half turn about an edge),  r(tau + 1) = e(1/5) r.
   Ramanujan's r(i) = sqrt((5 + sqrt5)/2) - phi is the fixed point of that half turn (an edge midpoint).
4. The face equation:  1/r^5 - 11 - r^5 = (eta(tau)/eta(5 tau))^6   (exact q-series).
5. t = r^5 is the Hauptmodul of Apery's zeta(2) family on Gamma_1(5): sum u_n t^n = 1 + sum c(n) q^n/(1 - q^n),
   c = (3, 1, -1, -3) for n = 1, 2, 3, 4 mod 5, u_n = sum_k C(n,k)^2 C(n+k,k).  Its singular points t^2 + 11t - 1 = 0 are
   phi^-5 and -phi^5: the two rings of five faces.  Apery's inequality is phi^5 > e^2.
6. The three Platonic faces against the lattice principle: tr(T L_w^-1) = 2 - w is elliptic, parabolic, hyperbolic at 3, 4, 5.
7. Volumes: the regular ideal octahedron has volume 4G, the regular ideal tetrahedron (3 sqrt3/4) L(2, chi_-3).
Usage: python ramanujan_dodecahedron.py"""
from fractions import Fraction as Fr
from math import comb, gcd
import mpmath as mp

mp.mp.dps = 40
phi = (1 + mp.sqrt(5))/2

print("1. Euler's formula with Ramanujan's 1/6 = -2 zeta(-1)")
def psl2_order(N):
    o = N**3
    for p in range(2, N + 1):
        if N % p == 0 and all(p % d for d in range(2, p)): o = o*(p*p - 1)//(p*p)
    return o//2 if N > 2 else o
for N in range(2, 9):
    mu = psl2_order(N); chi = mu*(Fr(1, N) - Fr(1, 6))
    print("   N = %d: |PSL2(Z/N)| = %4d, cusps %3d, chi = mu (1/N - 1/6) = %4s, genus %s"
          % (N, mu, mu//N, chi, (2 - chi)/2))

print("\n2. Ramanujan sums over residue classes mod 5 (exact)")
B2 = lambda x: x*x - x + Fr(1, 6)
S = {a: -Fr(5, 2)*B2(Fr(a, 5)) for a in range(1, 5)}; S[0] = 5*Fr(-1, 12)
print("   S_a = sum_{n = a mod 5} n: %s" % {a: str(S[a]) for a in range(5)})
print("   total = %s (zeta(-1));  G: (S_1 + S_4)/2 = %s;  H: (S_2 + S_3)/2 = %s;  difference = %s;  -L(-1, chi_5)/2 = %s"
      % (sum(S.values()), (S[1] + S[4])/2, (S[2] + S[3])/2, (S[2] + S[3] - S[1] - S[4])/2, -(S[1] + S[4] - S[2] - S[3])/2))

print("\n3. The icosahedral coordinate")
def GH(tau):
    q = mp.exp(2j*mp.pi*tau)
    G = mp.mpc(1); H = mp.mpc(1)
    for n in range(1, 400):
        G /= (1 - q**(5*n - 4))*(1 - q**(5*n - 1)); H /= (1 - q**(5*n - 3))*(1 - q**(5*n - 2))
    return mp.matrix([mp.exp(2j*mp.pi*tau*mp.mpf(-1)/60)*G, mp.exp(2j*mp.pi*tau*mp.mpf(11)/60)*H])
t1, t2, t3 = mp.mpc(0.1, 1.1), mp.mpc(-0.2, 0.8), mp.mpc(0.37, 0.6)
V = mp.matrix(2, 2); W = mp.matrix(2, 2)
for j, tt in enumerate((t1, t2)):
    a = GH(tt); b = GH(-1/tt)
    V[0, j], V[1, j] = a[0], a[1]; W[0, j], W[1, j] = b[0], b[1]
M = W*V**-1
err = mp.norm(GH(-1/t3) - M*GH(t3))
print("   S-matrix of (q^-1/60 G, q^11/60 H), times sqrt5/2: [[%s, %s], [%s, %s]]   (sin 2pi/5 = %s, sin pi/5 = %s)"
      % (mp.nstr(M[0, 0]*mp.sqrt(5)/2, 12), mp.nstr(M[0, 1]*mp.sqrt(5)/2, 12), mp.nstr(M[1, 0]*mp.sqrt(5)/2, 12),
         mp.nstr(M[1, 1]*mp.sqrt(5)/2, 12), mp.nstr(mp.sin(2*mp.pi/5), 12), mp.nstr(mp.sin(mp.pi/5), 12)))
print("   linearity tested at a third point: error %s" % mp.nstr(err, 3))
r = lambda tau: GH(tau)[1]/GH(tau)[0]
print("   r(-1/tau) - (1 - phi r)/(phi + r) at tau = %s: %s" % (t3, mp.nstr(abs(r(-1/t3) - (1 - phi*r(t3))/(phi + r(t3))), 3)))
print("   r(tau + 1)/r(tau) - e(1/5): %s" % mp.nstr(abs(r(t3 + 1)/r(t3) - mp.exp(2j*mp.pi/5)), 3))
print("   r(i) = %s;  sqrt((5 + sqrt5)/2) - phi = %s" % (mp.nstr(r(mp.mpc(0, 1)).real, 25), mp.nstr(mp.sqrt((5 + mp.sqrt(5))/2) - phi, 25)))

print("\n4-5. Exact q-series")
P = 60
def prod_series(expo):
    s = [0]*P; s[0] = 1
    for n in range(1, P):
        e = expo(n)
        for _ in range(abs(e)):
            if e > 0:
                for k in range(P - 1, n - 1, -1): s[k] -= s[k - n]
            else:
                for k in range(n, P): s[k] += s[k - n]
    return s
def mul(a, b):
    r_ = [0]*P
    for i, x in enumerate(a):
        if x:
            for j in range(P - i): r_[i + j] += x*b[j]
    return r_
def inv(a):
    r_ = [0]*P; r_[0] = 1
    for n in range(1, P): r_[n] = -sum(a[k]*r_[n - k] for k in range(1, n + 1))
    return r_
leg = lambda n: (0, 1, -1, -1, 1)[n % 5]
Pser = prod_series(lambda n: 5*leg(n))                       # t = r^5 = q * Pser
lhs = inv(Pser)                                                # q/t
lhs = [lhs[k] - (11 if k == 1 else 0) - (Pser[k - 2] if k >= 2 else 0) for k in range(P)]   # q (1/t - 11 - t)
e1 = prod_series(lambda n: 6); e5 = [0]*P                     # prod (1 - q^n)^6
for k in range(P):
    if 5*k < P: e5[5*k] = e1[k]                               # prod (1 - q^{5n})^6
rhs = mul(e1, inv(e5))                                         # q (eta(tau)/eta(5tau))^6
print("   1/r^5 - 11 - r^5 = (eta(tau)/eta(5tau))^6: %s  (%d coefficients)" % (lhs == rhs, P))
t = [0] + Pser[:P - 1]
c = lambda n: (0, 3, 1, -1, -3)[n % 5]
f = [0]*P; f[0] = 1
for d in range(1, P):
    for m in range(d, P, d): f[m] += c(d)
u = [sum(comb(n, k)**2*comb(n + k, k) for k in range(n + 1)) for n in range(P)]
acc = [0]*P; tp = [1] + [0]*(P - 1)
for n in range(P):
    acc = [x + u[n]*y for x, y in zip(acc, tp)]; tp = mul(tp, t)
print("   Apery numbers %s ...: sum u_n t^n = 1 + sum c(n) q^n/(1-q^n) with t = r^5: %s" % (u[:6], acc == f))
print("   recurrence (n+1)^2 u_{n+1} = (11n^2+11n+3) u_n + n^2 u_{n-1}: %s"
      % all((n + 1)**2*u[n + 1] == (11*n*n + 11*n + 3)*u[n] + n*n*u[n - 1] for n in range(1, P - 1)))
print("   roots of t^2 + 11t - 1: %s, %s;  phi^-5 = %s, -phi^5 = %s;  phi^5 = %s > e^2 = %s"
      % (mp.nstr((-11 + 5*mp.sqrt(5))/2, 12), mp.nstr((-11 - 5*mp.sqrt(5))/2, 12), mp.nstr(phi**-5, 12), mp.nstr(-phi**5, 12),
         mp.nstr(phi**5, 8), mp.nstr(mp.e**2, 8)))

print("\n6. The Platonic faces and the two-parabolic group <T, L_w>")
for w, face in ((3, "triangle (tetrahedron)"), (4, "square (cube)"), (5, "pentagon (dodecahedron)"), (6, "hexagon (flat)"), (8, "octagon")):
    tr = 2 - w
    kind = "elliptic" if abs(tr) < 2 else "parabolic" if abs(tr) == 2 else "hyperbolic"
    print("   w = %d  %-24s tr(T L_w^-1) = %2d  %s" % (w, face, tr, kind))
print("   w = 5: eigenvalues of T L_5^-1 are -phi^2, -phi^-2: %s" % mp.nstr(abs(phi**2 + phi**-2 - 3), 3))

print("\n7. Volumes of the regular ideal solids")
G = mp.catalan; L3 = mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])
print("   octahedron: 8 Lambda(pi/4) = 4 Cl_2(pi/2) = %s;  4G = %s" % (mp.nstr(4*mp.clsin(2, mp.pi/2), 20), mp.nstr(4*G, 20)))
print("   tetrahedron: 3 Lambda(pi/3) = %s;  (3 sqrt3/4) L(2, chi_-3) = %s"
      % (mp.nstr(mp.mpf(3)/2*mp.clsin(2, 2*mp.pi/3), 20), mp.nstr(3*mp.sqrt(3)/4*L3, 20)))

print("\n8. Arithmetic only: 1 - 1/24 - 1/60 = %s" % (1 - Fr(1, 24) - Fr(1, 60)))
