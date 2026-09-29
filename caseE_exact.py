"""Case E exactly: W = (Dt)^2 / (t p(t) f) in the Eisenstein basis of M_3(Gamma_0(8), chi_-4), by exact rational linear algebra.
Prediction from the Eichler-integral derivation:  Apery limit = L(W,2);  log amplitude at the cusp 1/4 = -alpha*pi/8
(alpha = coefficient of E_G(tau)).  Result: W = E_G(tau) - 8 E_G(2 tau) = -E_G(tau + 1/2), limit G/2, amplitude -pi/8
(the check n F_n / 4^n -> -pi/8 is recorded in the results log, section "Case E exactly")."""
import sys
from fractions import Fraction as Fr
sys.argv = ['x', '8']
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
import level12_modular as L
P = 60
chi = lambda n: L.chi(-4, n)
# t: Hauptmodul with pole at the cusp 1/2 (case E frame), f = theta_3^2
r = {1: 4, 2: -10, 4: 2, 8: 4}
t = [0] + L.eta_series(r, P)[:P - 1]
f = L.eis1(-4, 1, P)
Dt = [n*x for n, x in enumerate(t)]
p = L.mul([1, -12, 32] + [0]*(P - 3), [1] + [0]*(P - 1), P)       # p(t) as a polynomial in t: 1 - 12t + 32t^2
pt = [0]*P
tp = [1] + [0]*(P - 1)
for c in (1, -12, 32):
    pt = [a + c*b for a, b in zip(pt, tp)]
    tp = L.mul(tp, t, P)
num = L.mul(Dt, Dt, P)                                               # q^2 + ...
den = L.mul(L.mul(t, pt, P), f, P)                                   # q + ...
# W = num/den: divide by q in both
num1, den1 = num[1:] + [0], den[1:] + [0]
W = L.mul(num1, [Fr(x) for x in L.inv(den1, P)], P) if den1[0] in (1, -1) else None
print('W =', [str(x) for x in W[:12]])
EG = [0] + [sum(chi(n//d)*d*d for d in range(1, n + 1) if n % d == 0) for n in range(1, P)]
EZ = [Fr(-1, 4)] + [sum(chi(d)*d*d for d in range(1, n + 1) if n % d == 0) for n in range(1, P)]
def resc(s, k): return [s[i//k] if i % k == 0 else 0 for i in range(P)]
basis = {'E_G(t)': EG, 'E_G(2t)': resc(EG, 2), 'E_Z(t)': EZ, 'E_Z(2t)': resc(EZ, 2)}
from sympy import Matrix, Rational
names = list(basis)
R = P - 3                                                            # the shifted quotient is exact below P - 1
M = Matrix([[Rational(str(basis[nm][n])) for nm in names] for n in range(R)])
rhs = Matrix([Rational(str(W[n])) for n in range(R)])
sol, params = M.gauss_jordan_solve(rhs)
print('W =', ' + '.join('(%s) %s' % (sol[i], nm) for i, nm in enumerate(names)), '  free params:', params)
alpha, beta, gamma, delta = [Fr(int(x.p), int(x.q)) for x in sol]
print('predicted limit L(W,2) = -(alpha + beta/4) G/2 + (3 gamma/4)(pi^2/12)  ->  G-coefficient %s, pi^2-coefficient %s' % (-(alpha + beta/4)/2, Fr(3, 4)*gamma/12 if gamma == -delta else 'n/a'))
print('predicted log amplitude at t = 1/4:  lim n F_n / 4^n = -alpha pi / 8 = %s * pi' % (-alpha/8))
