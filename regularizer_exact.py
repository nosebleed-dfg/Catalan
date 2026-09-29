"""The "infinity minus 1" regularizer, exactly (2026-09-29; the user chose: use the blow-up itself, 1/(1-e(y)) is the
infinity, subtracting 1 regularizes it to -1/2 which carries G; build the approximations from the divergent part,
subtract the exact blow-up, bound what is left for every n by exact linear algebra).
All statements are exact (rational linear algebra on recurrences and q-expansions); no PSLQ, no fitted decay.
  1. Case E: the cokernel of the recurrence operator on t*Q[t] is 1-dimensional, so every inhomogeneity h gives
     y_h = kappa_h * y_t + P_h(t) - P_h(0) f  (P_h a polynomial): one limit, one blow-up, fixed ratio.
  2. Case E: Casoratian C_n = u_n v_{n+1} - u_{n+1} v_n = 32^n/(n+1)^2 exactly, hence
     G/2 - v_n/u_n = sum_{k>=n} 32^k / ((k+1)^2 u_k u_{k+1})   for every n.
  3. Gamma_1(8): the cokernel is 3-dimensional (T, T^2, T^3), with the exact reduction of T^4, T^5, T^6; and the
     weight-3 form of the class T is exactly the CM cusp form g = eta(t)^2 eta(2t) eta(4t) eta(8t)^2."""
import sys, io, contextlib
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], '8']
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
import level12_modular as L
from guess_rec import guess
from sympy import Matrix, Rational

def ev(c, n): return sum(ci*n**e for e, ci in enumerate(c))

def cokernel(cs, hdeg):
    """Recurrence sum_k c_k(n-R) y_{n-R+k} = h_n (backward form).  Acting on polynomials P = sum p_j T^j gives
    (LP)_n = sum_k c_k(n-R) p_{n-R+k}.  For each monomial h = T^m (1 <= m <= hdeg) solve LP + sum_{i in basis} k_i T^i = T^m
    with basis = the smallest set of monomials that spans the cokernel."""
    R = len(cs) - 1
    def L_of_monomial(j, top):
        out = [Fr(0)]*(top + 1)
        for k in range(R + 1):
            n = j + R - k                      # p_{n-R+k} = p_j  <=>  n = j + R - k
            if 0 <= n <= top:
                out[n] += ev(cs[k], n - R)
        return out
    top = hdeg + R + 2
    images = [L_of_monomial(j, top) for j in range(0, hdeg + 1)]
    # find the cokernel basis greedily: monomials T^m (m >= 1) not in span(images + previous basis)
    basis, rows = [], [im[1:] for im in images]   # restrict to T^1..T^top (h(0) = 0; image constant terms are 0 anyway?)
    def rank(vecs): return Matrix([[Rational(x.numerator, x.denominator) for x in v] for v in vecs]).rank() if vecs else 0
    cur = rows[:]
    r0 = rank(cur)
    for m in range(1, hdeg + 1):
        e = [Fr(0)]*top; e[m - 1] = Fr(1)
        if rank(cur + [e]) > r0:
            basis.append(m); cur.append(e); r0 += 1
    return basis, images

# ---- 1, 2: case E ----------------------------------------------------------------------------------------------
# forward form of (n+1)^2 a_{n+1} = (12n^2+12n+4) a_n - 32 n^2 a_{n-1}:  c_0(n) a_n + c_1(n) a_{n+1} + c_2(n) a_{n+2} = 0
csE = [[32, 64, 32], [-28, -36, -12], [4, 4, 1]]      # 32(n+1)^2, -(12(n+1)^2+12(n+1)+4), (n+2)^2
basisE, imE = cokernel(csE, 6)
print('case E: cokernel basis (monomials t^m):', basisE)
# explicit reduction of t^2: L(P) + kappa t = t^2 with P constant
# L(1) as a polynomial:
print('   L(1) =', [str(x) for x in imE[0][:4]], ' L(t) =', [str(x) for x in imE[1][:5]])
# Casoratian
N = 40
u = [Fr(1), Fr(4)]; v = [Fr(0), Fr(1)]
for n in range(1, N):
    u.append(((12*n*n + 12*n + 4)*u[n] - 32*n*n*u[n - 1])/((n + 1)**2))
    v.append(((12*n*n + 12*n + 4)*v[n] - 32*n*n*v[n - 1])/((n + 1)**2))
ok = all(u[n]*v[n + 1] - u[n + 1]*v[n] == Fr(32**n, (n + 1)**2) for n in range(N))
print('case E Casoratian u_n v_{n+1} - u_{n+1} v_n = 32^n/(n+1)^2 for n < %d: %s' % (N, ok))

# ---- 3: Gamma_1(8) ---------------------------------------------------------------------------------------------
P = 120
def theta3(m):
    s = [0]*P; k = 0
    while m*k*k < P:
        s[m*k*k] += 1 if k == 0 else 2; k += 1
    return s
th1, th2 = theta3(1), theta3(2)
s_ = L.mul(th2, L.inv(th1, P), P)
T = [((1 if i == 0 else 0) - s_[i])//2 for i in range(P)]
f = L.mul(th1, th1, P)
u8 = L.expand_in(f, T, 100)
with contextlib.redirect_stdout(io.StringIO()):
    R8, d8, cs8 = guess(u8, R=8, D=2, verbose=False)
basis8, im8 = cokernel(cs8, 8)
print('Gamma_1(8): recurrence order %d; cokernel basis (monomials T^m): %s' % (R8, basis8))
# reductions of T^m for m not in the basis: solve  sum_j p_j L(T^j) + sum_{i in basis} k_i T^i = T^m
top = 8 + R8 + 2
for m in range(1, 9):
    if m in basis8: continue
    cols = [im[1:top + 1] for im in im8[:m]] + [[Fr(1) if n == i - 1 else Fr(0) for n in range(top)] for i in basis8]
    A = Matrix([[Rational(c[n].numerator, c[n].denominator) for c in cols] for n in range(top)])
    b = Matrix([1 if n == m - 1 else 0 for n in range(top)])
    sol, params = A.gauss_jordan_solve(b)
    ks = sol[m:]
    print('   T^%d = L(P) + %s   (P of degree < %d)' % (m, ' + '.join('(%s) T^%d' % (k, i) for k, i in zip(ks, basis8)), m))
# the weight-3 form of the class T: W0 = (DT)^2 / (T p(T) f), p(T) = sum_k lc(c_{R-k}) T^k from the n^2 coefficients
pc = [cs8[R8 - k][2] if len(cs8[R8 - k]) > 2 else 0 for k in range(R8 + 1)]
print('   p(T) coefficients (T^0..T^%d):' % R8, pc)
pT = [0]*P; Tk = [1] + [0]*(P - 1)
for c in pc:
    pT = [a + c*b for a, b in zip(pT, Tk)]; Tk = L.mul(Tk, T, P)
DT = [n*x for n, x in enumerate(T)]
num = L.mul(DT, DT, P)[2:] + [0, 0]                   # (DT)^2 / q^2
den = L.mul(L.mul(T, pT, P), f, P)[1:] + [0]           # T p f / q
W0 = [0] + L.mul(num, L.inv(den, P), P)[:P - 1]        # times q
g = [0] + L.eta_series({1: 2, 2: 1, 4: 1, 8: 2}, P)[:P - 1]
print('   W0 =', W0[:12])
print('   g  =', g[:12])
print('   W0 == g on %d coefficients: %s' % (P - 5, W0[:P - 5] == g[:P - 5]))
