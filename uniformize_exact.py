"""Uniformizing at infinity (2026-09-29; the user chose: "Uniformize at infinity" -- the pi in the cost is the 2 pi i of the
cusp coordinate; do the bound in the modular variable).
Part 1 (this file, section A): on Gamma_1(8), decompose the pure-G weight-3 form W = T g - g/2 (the class of h = T^2 - T/2)
exactly in the Eisenstein basis, and evaluate its additively twisted values L(W, e(x.), 2) at cusp representatives x by the
exact finite formulas (Ramanujan-regularized: Li_0(e(y)) = -1/2 + (i/2)cot(pi y), Hurwitz zeta(0, a) = 1/2 - a):
  E^{psi,1}(d t):  d^-2 sum_{r mod P} psi(r) Li_0(e(r d x)) zeta(2, r/P)/P^2
  E^{1,psi}(d t):  d^-2 sum_{r=1..P'} psi(r) Li_2(e(r d x)) (1/2 - r/P')
A representative x is harmless iff W vanishes there (true for W = T g - g/2 at every cusp but 1/2) and L(W, e(x.), 2) = Lambda.
Part 2 (section B): the conformal radius of the surface on which F lives, R* = exp(2 pi^2 sum 1/c^2) over the double cosets of
Gamma' = <translation, parabolics at harmless representatives> (Green's function = 2 pi E_{Gamma'}(tau, 1))."""
import sys, math, itertools
from fractions import Fraction as Fr
import mpmath as mp
sys.argv = [sys.argv[0], '8'] + sys.argv[1:]
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
import level12_modular as L
from sympy import Matrix, Rational
mp.mp.dps = 40

def chi_m4(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi_m8(n): return {1: 1, 3: 1, 5: -1, 7: -1}.get(n % 8, 0)
def one(n): return 1
CHARS = {'chi-4': (chi_m4, 4), 'chi-8': (chi_m8, 8)}

def section_A():
    P = 80
    def theta3(m):
        s = [0]*P; k = 0
        while m*k*k < P:
            s[m*k*k] += 1 if k == 0 else 2; k += 1
        return s
    th1, th2 = theta3(1), theta3(2)
    s_ = L.mul(th2, L.inv(th1, P), P)
    T = [((1 if i == 0 else 0) - s_[i])//2 for i in range(P)]
    g = [0] + L.eta_series({1: 2, 2: 1, 4: 1, 8: 2}, P)[:P - 1]
    W = [Fr(a) - Fr(b, 2) for a, b in zip(L.mul(T, g, P), g)]
    # Eisenstein basis of M_3(Gamma_1(8)) (odd characters mod 8) and the cusp form g
    def E_psi1(psi, d):
        s = [Fr(0)]*P
        for n in range(1, (P - 1)//d + 1):
            s[d*n] = Fr(sum(psi(n//j)*j*j for j in range(1, n + 1) if n % j == 0))
        return s
    def E_1psi(psi, f, d):
        B3 = sum(psi(a)*(Fr(a, f)**3 - Fr(3, 2)*Fr(a, f)**2 + Fr(a, 2*f)) for a in range(1, f + 1))*f*f   # B_{3,psi}
        s = [Fr(0)]*P
        s[0] = -B3/3/2                                  # L(-2, psi)/2 = -B_{3,psi}/6
        for n in range(1, (P - 1)//d + 1):
            s[d*n] = Fr(sum(psi(j)*j*j for j in range(1, n + 1) if n % j == 0))
        return s
    basis = [('E^{chi-4,1}(t)', E_psi1(chi_m4, 1), ('psi1', 'chi-4', 1)), ('E^{chi-4,1}(2t)', E_psi1(chi_m4, 2), ('psi1', 'chi-4', 2)),
             ('E^{1,chi-4}(t)', E_1psi(chi_m4, 4, 1), ('1psi', 'chi-4', 1)), ('E^{1,chi-4}(2t)', E_1psi(chi_m4, 4, 2), ('1psi', 'chi-4', 2)),
             ('E^{chi-8,1}(t)', E_psi1(chi_m8, 1), ('psi1', 'chi-8', 1)), ('E^{1,chi-8}(t)', E_1psi(chi_m8, 8, 1), ('1psi', 'chi-8', 1)),
             ('g', [Fr(x) for x in g], None)]
    R = P - 5
    M = Matrix([[Rational(b[1][n].numerator, b[1][n].denominator) for b in basis] for n in range(R)])
    rhs = Matrix([Rational(W[n].numerator, W[n].denominator) for n in range(R)])
    sol, params = M.gauss_jordan_solve(rhs)
    coef = [Fr(int(x.p), int(x.q)) for x in sol]
    print('W = T g - g/2 =', ' + '.join('(%s) %s' % (c, b[0]) for c, b in zip(coef, basis) if c), '   free parameters:', params.shape)

    def li0(y):                     # Sum_{d>=1} e(d y), Ramanujan/Abel; y a Fraction
        y = y - (y.numerator // y.denominator)
        if y == 0: return mp.mpf(-0.5)
        return mp.mpf(-0.5) + 0.5j*mp.cot(mp.pi*mp.mpf(y.numerator)/y.denominator)
    def twisted(spec, x):
        kind, cname, d = spec
        psi, fcond = CHARS[cname]
        y = Fr(d)*x
        den = y.denominator
        if kind == 'psi1':
            Pp = fcond*den//math.gcd(fcond, den)
            tot = mp.mpc(0)
            for r in range(1, Pp + 1):
                if psi(r) == 0: continue
                tot += psi(r)*li0(Fr(r)*y)*mp.zeta(2, mp.mpf(r)/Pp)/Pp**2
            return tot/d**2
        else:
            Pp = fcond*den//math.gcd(fcond, den)
            tot = mp.mpc(0)
            for r in range(1, Pp + 1):
                if psi(r) == 0: continue
                z = Fr(r)*y
                tot += psi(r)*mp.polylog(2, mp.expjpi(2*mp.mpf(z.numerator)/z.denominator))*(mp.mpf(1)/2 - mp.mpf(r)/Pp)
            return tot/d**2
    G = +mp.catalan
    xs = [Fr(0), Fr(1, 4), Fr(3, 4), Fr(1, 3), Fr(2, 3), Fr(3, 8), Fr(5, 8), Fr(1, 5), Fr(1, 7), Fr(3, 7), Fr(1, 9), Fr(1, 11)]
    vals = {}
    for x in xs:
        v = sum(mp.mpf(c.numerator)/c.denominator*twisted(b[2], x) for c, b in zip(coef, basis) if c and b[2] is not None)
        vals[x] = v
    Lam = vals[Fr(0)]
    print('Lambda = L(W,2) =', mp.nstr(Lam, 20), '   -G/4 =', mp.nstr(-G/4, 20))
    for x in xs[1:]:
        dv = vals[x] - Lam
        print('  x = %-5s  L(W,e(x.),2) - Lambda = %s   %s' % (x, mp.nstr(dv, 15), 'HARMLESS' if abs(dv) < mp.mpf(10)**-30 else 'pays'))
    return coef

if __name__ == '__main__':
    section_A()
