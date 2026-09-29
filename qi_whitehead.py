r"""The Q(i) side, made exact (2026-09-29): the Whitehead link W and the figure-eight knot 4_1 against the lattice principle.

A. vol/(2 pi) = L'(-1, chi_D) for M = 4_1 (D = -3) and W (D = -4): Humbert's covolume |D|^{3/2} zeta_K(2)/(4 pi^2) equals
   -2 pi zeta_K'(-1) = -2 pi zeta(-1) L'(-1, chi_D), and both groups have index 12 = -1/zeta(-1) in PSL_2(O_K).
B. Riley: the 2-bridge groups b(5,3) = 4_1 and b(8,3) = W are two-parabolic groups <T, (1 0; u 1)>, the same family as
   our harmless groups <T, (1 0; w 1)> (real w).  Solve for u exactly (sympy; Riley's word needs q odd).
C. Inside W's group: a pair of parabolics with parabolic product, i.e. a thrice-punctured sphere, i.e. the w = 4 point
   <T, L_4> = Gamma_0(4) up to conjugacy over Q(i).
D. Kashaev invariant of W (Murakami-Murakami-Okamoto-Takata-Yokota, Exp. Math. 11 (2002), formula for J_N(W)):
   J_N = sum_{k} A_k^2 / (q)_k^4,  A_k = sum_{i=k}^{N-1} (qbar)_i^2 / (qbar)_{i-k},  q = e^{2 pi i/N}
   (the i- and j-sums of their triple sum factor for fixed k; their printed prefactor q^{-(N-1)N/2} = (-1)^{N-1} is
   dropped, which reproduces their Table 1).  Checks: Table 1; growth 2G/pi; constant and first corrections of
   J_N / (N^{3/2} exp(N (4G + i pi^2/4)/(2 pi))); norms over Q(zeta_N) are rational integers.
Usage: python qi_whitehead.py"""
import itertools
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
G = mp.catalan
print("A. Volumes as derivatives at the trivial zero")
for name, k, chi, vol in (("4_1", 3, [0, 1, -1], 2*mp.im(mp.polylog(2, mp.exp(1j*mp.pi/3)))),
                          ("Whitehead", 4, [0, 1, 0, -1], 4*G)):
    Lp = mp.dirichlet(-1, chi, 1)
    L2 = mp.dirichlet(2, chi)
    zK2 = mp.zeta(2)*L2
    covol = mp.mpf(k)**1.5*zK2/(4*mp.pi**2)
    print("   %-9s vol = %s   vol/(2 pi) = %s   L'(-1,chi) = %s   k^{3/2} L(2,chi)/(4 pi) = %s"
          % (name, mp.nstr(vol, 15), mp.nstr(vol/(2*mp.pi), 15), mp.nstr(Lp, 15), mp.nstr(k**1.5*L2/(4*mp.pi), 15)))
    print("             covol PSL2(O_K) = %s,  12*covol = %s,  -2 pi zeta(-1) L'(-1,chi) = %s"
          % (mp.nstr(covol, 12), mp.nstr(12*covol, 12), mp.nstr(-2*mp.pi*mp.zeta(-1)*Lp, 12)))

print("\nB. Riley: 2-bridge groups as <x = T, y = (1 0; u 1)>")
u = sp.symbols('u')
X = sp.Matrix([[1, 1], [0, 1]])
Y = sp.Matrix([[1, 0], [u, 1]])
def word(p, q):
    eps = [(-1)**((k*q)//p) for k in range(1, p)]
    # knots (p odd): w = x^e1 y^e2 x^e3 ...; links (p even): w = y^e1 x^e2 y^e3 ...
    first, second = (X, Y) if p % 2 else (Y, X)
    M = sp.eye(2)
    for k, e in enumerate(eps):
        g = first if k % 2 == 0 else second
        M = M * (g if e == 1 else g.inv())
    return sp.simplify(M)
riley = {}
for name, p, q in (("4_1 = b(5,3)", 5, 3), ("Whitehead = b(8,3)", 8, 3)):     # Riley's word needs q odd
    Wm = word(p, q)
    rel = (Wm*X - Y*Wm) if p % 2 else (Wm*X - X*Wm)
    g = sp.gcd_list([sp.expand(e) for e in rel if sp.expand(e) != 0])
    poly = sp.factor(g)
    roots = [r for r in sp.solve(g, u) if r != 0]
    riley[name] = roots
    print("   %-20s relation polynomial %s ; nonzero roots %s" % (name, poly, roots))

print("\nC. A thrice-punctured sphere (the w = 4 point, Gamma_0(4)) inside the Whitehead group")
uW = [r for r in riley["Whitehead = b(8,3)"] if sp.im(r) > 0][0]
x, y = X, Y.subs(u, uW)
gens = {"x": x, "X": x.inv(), "y": y, "Y": y.inv()}
found = None
for L in range(1, 4):
    for wd in itertools.product("xXyY", repeat=L):
        gm = sp.eye(2)
        for c in wd: gm = gm*gens[c]
        P2 = sp.simplify(gm*y*gm.inv())
        if P2 == y or P2 == y.inv(): continue
        for s in (1, -1):
            prod = sp.simplify(y*(P2 if s == 1 else P2.inv()))
            if sp.simplify(prod.trace()) == -2:
                found = ("".join(wd), s, P2); break
        if found: break
    if found: break
wd, s, P2 = found
print("   u_W = %s;  P1 = y, P2 = g y g^-1 with g = %s;  tr(P1 * P2^%d) = -2" % (uW, wd, s))
# normalize: y fixes 0, P2 fixes g(0); send 0 -> oo and g(0) -> 0
g0 = sp.simplify((gm[0, 1])/(gm[1, 1])) if gm[1, 1] != 0 else sp.oo
M = sp.Matrix([[1, -g0], [1, 0]]) if g0 != sp.oo else sp.Matrix([[0, -1], [1, 0]])
M = M / sp.sqrt(M.det())
A1 = sp.simplify(M*y*M.inv()); A2 = sp.simplify(M*(P2 if s == 1 else P2.inv())*M.inv())
beta, gamma = A1[0, 1], A2[1, 0]
print("   conjugated: P1 -> (1 %s; 0 1), P2^%d -> (%s %s; %s %s);  invariant beta*gamma = %s"
      % (sp.simplify(A1[0, 1]), s, *[sp.simplify(e) for e in A2], sp.simplify(beta*gamma)))

print("\nD. Kashaev invariant of the Whitehead link")
def J(N, a=1, dps=None):
    with mp.workdps(dps or (40 + N//3)):
        q = mp.exp(2j*mp.pi*a/N); qb = 1/q
        P = [mp.mpc(1)]; Pb = [mp.mpc(1)]
        for j in range(1, N):
            P.append(P[-1]*(1 - q**j)); Pb.append(Pb[-1]*(1 - qb**j))
        tot = mp.mpc(0)
        for k in range(N):
            A = mp.fsum(Pb[i]**2/Pb[i - k] for i in range(k, N))
            tot += A*A/P[k]**4
        # MMOTY print a prefactor q^{-(N-1)N/2} = (-1)^{N-1}; their Table 1 matches the sum without it
        return +tot
table = {40: "3.892920359101811097809525583+2.457483997330866045812504703j",
         50: "3.848161466402914225154530180+2.461039474018016569869745301j"}
for N, ref in table.items():
    val = 2*mp.pi*mp.log(J(N + 1)/J(N))
    print("   N = %d: 2 pi log(J_{N+1}/J_N) = %s   (MMOTY Table 1: %s)" % (N, mp.nstr(val, 20), ref))
V = 4*G + 1j*mp.pi**2/4
Ns = list(range(160, 401, 20))
Rs = []
for N in Ns:
    JN = J(N)
    R = JN/(mp.mpf(N)**1.5*mp.exp(N*V/(2*mp.pi)))
    Rs.append(R)
    print("   N = %3d  log|J_N|/N = %s  (2G/pi = %s)   R_N = J_N/(N^1.5 e^{N V/2pi}) = %s"
          % (N, mp.nstr(mp.log(abs(JN))/N, 10), mp.nstr(2*G/mp.pi, 10), mp.nstr(R, 15)))
# Richardson in 1/N: fit R_N = C (1 + b1/N + b2/N^2 + ...) with len(Ns) terms
n = len(Ns)
Amat = mp.matrix([[mp.mpf(1)/mp.mpf(N)**j for j in range(n)] for N in Ns])
coef = mp.lu_solve(Amat, mp.matrix(Rs))
C = coef[0]
print("   extrapolated C = %s ;  |C| = %s ;  C^2 = %s ;  C^4 = %s" % (mp.nstr(C, 12), mp.nstr(abs(C), 12), mp.nstr(C**2, 12), mp.nstr(C**4, 12)))
h1 = coef[1]/C/(2j*mp.pi); h2 = coef[2]/C/(2j*mp.pi)**2; h3 = coef[3]/C/(2j*mp.pi)**3
print("   in h = 2 pi i/N:  R = C (1 + a1 h + a2 h^2 + a3 h^3 + ...)")
for lab, val in (("a1", h1), ("a2", h2), ("a3", h3)):
    print("   %s = %s   x 96 = %s   x 96^2 = %s" % (lab, mp.nstr(val, 16), mp.nstr(96*val, 14), mp.nstr(96**2*val, 12)))

print("\nE. Integrality at roots of unity: norms of J_N(zeta_N) over Q")
from math import gcd
for N in (3, 4, 5, 6, 7, 8, 9, 10, 12):
    prod = mp.mpc(1)
    for a in range(1, N):
        if gcd(a, N) == 1: prod *= J(N, a, dps=60)
    nr = int(mp.nint(prod.real))
    print("   N = %2d  Norm = %s  (imag %s)  = %s" % (N, nr, mp.nstr(prod.imag, 3), sp.factorint(nr)))
