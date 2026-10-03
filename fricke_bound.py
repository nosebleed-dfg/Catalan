r"""fricke_bound.py (2026-10-02): the holonomy bound on the Fricke quotient X'' = H/<T, w_8> (filled cusp).

fricke_lever.py found a family F = f (g - Lambda), f = theta^2 + sqrt2 theta(2tau)^2, W = W_E + (E2 - E2(2tau))/sqrt2,
Lambda = G/2 + pi^2/(16 sqrt 2), invariant under the Hecke group Gamma'' = <T, w_8> (free product Z * Z/2), so F is a
holomorphic function on X''.  log R'' = 2 pi^2 (S_even + S_odd) = 4.2171 is the conformal radius of X'' at the cusp in the
coordinate s = t (1 - 8t)/(1 - 4t) (t the Hauptmodul of Gamma_0(8), s that of Gamma_0(8)+).
Calegari-Dimitrov-Tang (Theorem 2.5.1): m functions of rates sigma_i, meromorphic on the disc after a map phi, satisfy
    m <= I(phi)/(log|phi'(0)| - tau(b)),   I(phi) = double integral of log|phi(z) - phi(w)| over the torus,
    tau(b) = (1/m^2) sum (2i - 1) sigma_i;  the rearrangement integral int 2t (log|phi(e^{2 pi i t})|)* dt majorises I.
Here phi(z) = s(tau(r z)) with tau(z) the inverse of the uniformizing map z: X'' -> D, computed as
    z(tau) = exp(2 pi i (Z(tau) + C)),   Z(tau) = (tau - i) - (pi/8) Phi_tau(0),
    Phi(u) = K(u) + sum_{a != 0} (8 (u + a)^2)^-1 Phi(-1/(8(u + a))),   K_tau(u) = cot(pi(tau + u)) - cot(pi(i + u)),
from Z(tau) = sum over cosets (gamma tau - gamma i) = (tau - i) - pi sum over double cosets c^-2 [cot(pi(tau + x)) - cot(pi(i + x))],
x = d/c, and the double cosets of the free product are the words w_8 T^a1 w_8 ... with the ratio u = d/c evolving by
u -> -1/(8(u + a)) and c -> sqrt8 c (u + a).  The functional equation is solved by Chebyshev collocation with exact
Hurwitz-zeta sums (radius_transfer.py).  The same with K = 1 gives 8 (S_even + S_odd), a check.
Then for several radii r: the level curve |z| = r is traced by Newton continuation, the s-values on it give the
rearrangement and Bost-Charles integrals, and the bound is evaluated for the function systems
    {1, F} (tau(b) = 3/2), {1, sq, F, F sq} (3/2), {1, sq, F, F', F'', and their products with sq} (15/8), sq = sqrt(s - s_e2).
Usage: python fricke_bound.py"""
import sys
import mpmath as mp

mp.mp.dps = 55                                   # the monomial Vandermonde needs about N/3 digits of headroom
W = 8
rho = (1 - mp.sqrt(1 - mp.mpf(4)/W))/2            # the ratio u = d/c stays in [-rho, rho]
N = 100
nodes = [rho*mp.cos(mp.pi*(2*i + 1)/(2*N)) for i in range(N)]
V = mp.matrix(N, N)
for i in range(N):
    for j in range(N): V[i, j] = (nodes[i]/rho)**j
Vinv = V**-1
H = mp.matrix(N, N)
for i in range(N):
    u = nodes[i]
    for j in range(N):
        n = 2 + j
        H[i, j] = (-1)**j*(mp.zeta(n, 1 + u) + (-1)**j*mp.zeta(n, 1 - u))/(W*(W*rho)**j)
A = Vinv*(mp.eye(N) - H*Vinv)**-1                 # Phi coefficients = A * K-values
a0 = [A[0, i] for i in range(N)]                  # Phi(0) = sum a0[i] K(u_i)
Phi1 = sum(a0)
S_total = Phi1/W
logR2 = 2*mp.pi**2*S_total
print("collocation: 8 (S_even + S_odd) = %s;  log R'' = %s" % (mp.nstr(Phi1, 20), mp.nstr(logR2, 20)))
cot = lambda w: 1/mp.tan(w)
I_ = mp.mpc(0, 1)
base = [cot(mp.pi*(I_ + u)) for u in nodes]
def Z(tau):
    return (tau - I_) - (mp.pi/W)*sum(a0[i]*(cot(mp.pi*(tau + nodes[i])) - base[i]) for i in range(N))
def dZ(tau):
    return 1 + (mp.pi**2/W)*sum(a0[i]/mp.sin(mp.pi*(tau + nodes[i]))**2 for i in range(N))
def E1(tau):
    """Eisenstein series E(tau, 1) of Gamma'' by the same operator."""
    x, y = mp.re(tau), mp.im(tau)
    return y + (mp.pi/W)*sum(a0[i]*mp.sinh(2*mp.pi*y)/(mp.cosh(2*mp.pi*y) - mp.cos(2*mp.pi*(x + nodes[i]))) for i in range(N))
c0 = Z(20*I_) - 20*I_
C = -c0 + I_*logR2/(2*mp.pi)
print("c0 = Z(tau) - tau at the cusp = %s;  Im c0 + E(i, 1) - log R''/(2 pi) = %s" % (mp.nstr(c0, 12), mp.nstr(mp.im(c0) + E1(I_) - logR2/(2*mp.pi), 3)))
def reduce(tau):
    """a representative of the Gamma''-orbit with |x| <= 1/2 and |tau| >= 1/sqrt8 (largest imaginary part)."""
    for _ in range(200):
        tau = tau - mp.nint(mp.re(tau))
        if abs(tau) >= 1/mp.sqrt(8) - mp.mpf(10)**(-25): return tau
        tau = -1/(W*tau)
    return tau
def z(tau): return mp.exp(2*mp.pi*I_*(Z(reduce(tau)) + C))
def dz(tau): return z(tau)*2*mp.pi*I_*dZ(tau)        # derivative at tau itself: only used at reduced points
t0 = mp.mpc(0.13, 0.41)
print("invariance of z: |z(tau+1)/z(tau) - 1| = %s, |z(-1/(8 tau))/z(tau) - 1| = %s, |z(tau/(8tau+1))/z(tau) - 1| = %s"
      % (mp.nstr(abs(z(t0 + 1)/z(t0) - 1), 3), mp.nstr(abs(z(-1/(8*t0))/z(t0) - 1), 3), mp.nstr(abs(z(t0/(8*t0 + 1))/z(t0) - 1), 3)))
print("|z(tau)| = exp(-2 pi E(tau,1)): %s;  z(tau) ~ q/R'' at the cusp: %s"
      % (mp.nstr(abs(z(t0)) - mp.exp(-2*mp.pi*E1(t0)), 3), mp.nstr(z(6*I_)/mp.exp(2*mp.pi*I_*6*I_)*mp.exp(logR2) - 1, 3)))
print("z on the imaginary axis is real: Im z(0.7 i) = %s;  z(i/sqrt8) = %s" % (mp.nstr(mp.im(z(0.7*I_)), 3), mp.nstr(z(I_/mp.sqrt(8)), 12)))

# ---------- the Hauptmoduln: t = (theta^2 - theta(2tau)^2)/(4 theta^2) is case E's coordinate (f = 1 + 4t + 20t^2 + ...),
# with t(0) = 1/8; u = theta(2tau)^2/theta^2 = 1 - 4t has u(w8 tau) = 1/(2u); s = -(1 - u)(1 - 2u)/(4u) = t(1 - 8t)/(1 - 4t).
def theta(tau):
    q = mp.exp(2*mp.pi*I_*tau)
    return 1 + 2*sum(q**(n*n) for n in range(1, 80))
def u_haupt(tau):
    tau = reduce(tau)
    return theta(2*tau)**2/theta(tau)**2
def t_haupt(tau): return (1 - u_haupt(tau))/4
def s_haupt(tau):
    u = u_haupt(tau); return -(1 - u)*(1 - 2*u)/(4*u)
tt = (1 - theta(2*t0)**2/theta(t0)**2)/4; tw = (1 - theta(-2/(8*t0))**2/theta(-1/(8*t0))**2)/4
print("Hauptmodul checks: t(w8 tau) - mu(t) = %s (mu(t) = (8t-1)/(32t-8));  s(w8 tau) - s(tau) = %s;  t = q + ...: %s"
      % (mp.nstr(tw - (8*tt - 1)/(32*tt - 8), 3), mp.nstr(s_haupt(-1/(8*t0)) - s_haupt(t0), 3), mp.nstr(t_haupt(8*I_)/mp.exp(2*mp.pi*I_*8*I_), 8)))
print("t at the cusp 0 (tau = 0.001 i): %s;  s there: %s" % (mp.nstr(t_haupt(0.001*I_), 8), mp.nstr(s_haupt(0.001*I_), 8)))
e1 = I_/mp.sqrt(8); e2 = mp.mpf(1)/3 + I_*mp.sqrt(2)/12
g2 = lambda tau: (8*tau - 3)/(24*tau - 8)
print("e2 = 1/3 + i sqrt2/12 is fixed by (8 -3; 24 -8): %s" % mp.nstr(abs(g2(e2) - e2), 3))
s_e1, s_e2 = s_haupt(e1), s_haupt(e2)
print("s(e1) = %s  ((3 - 2 sqrt2)/4 = %s);  s(e2) = %s  ((3 + 2 sqrt2)/4 = %s)"
      % (mp.nstr(s_e1, 15), mp.nstr((3 - 2*mp.sqrt(2))/4, 15), mp.nstr(s_e2, 15), mp.nstr((3 + 2*mp.sqrt(2))/4, 15)))
print("univalent slit plane C \\ [s(e2), oo): log(4 s(e2)) = %s" % mp.nstr(mp.log(4*mp.re(s_e2)), 10))
r1 = mp.re(z(e1))
print("the elliptic point e1 sits at |z| = %s" % mp.nstr(r1, 10))
sys.stdout.flush()

# ---------- level curves and the integrals ----------
def on_axis(r):
    """tau on the sigma-fixed curve (imaginary axis, then the arc) with z(tau) = r > 0."""
    if r < r1:
        lo, hi = mp.mpf(1)/mp.sqrt(8), mp.mpf(12)
        for _ in range(80):
            mid = (lo + hi)/2
            if mp.re(z(mid*I_)) > r: lo = mid
            else: hi = mid
        return mid*I_
    lo, hi = mp.mpf(0.02), mp.pi/2                  # arc tau = e^{i psi}/sqrt8, psi from pi/2 (e1) down to 0 (funnel)
    for _ in range(80):
        mid = (lo + hi)/2
        if mp.re(z(mp.exp(I_*mid)/mp.sqrt(8))) > r: lo = mid
        else: hi = mid
    return mp.exp(I_*mid)/mp.sqrt(8)
def level_curve(r, M):
    taus = []; tau = reduce(on_axis(r))
    for k in range(M + 1):
        target = r*mp.exp(I_*mp.pi*k/M)
        for it in range(60):
            step = (z(tau) - target)/dz(tau)
            tau = reduce(tau - step)
            if abs(step) < mp.mpf(10)**(-22): break
        if abs(z(tau) - target) > mp.mpf(10)**(-15):
            raise RuntimeError("Newton failed at r = %s, k = %d" % (r, k))
        taus.append(tau)
    return taus
def integrals(r, M=160):
    taus = level_curve(r, M)
    svals = [s_haupt(tau) for tau in taus]
    full = svals[:-1] + [mp.conj(v) for v in reversed(svals[1:])]   # the lower half by the reflection symmetry
    n = len(full)
    logs = sorted(mp.log(abs(v)) for v in full)
    rear = sum((2*(k + mp.mpf(1)/2)/n)*logs[k] for k in range(n))/n
    bc = sum(mp.log(abs(full[k] - full[l])) for k in range(n) for l in range(n) if k != l)/(n*(n - 1))
    smax = max(abs(v) for v in full); smin = min(abs(v) for v in full)
    ymin = min(mp.im(tau) for tau in taus)
    return rear, bc, smax, smin, ymin
print("\nradius r | log|phi'(0)| = log R'' + log r | rearrangement I_r | Bost-Charles I_bc | max|s|, min|s| on the circle | lowest Im tau")
rows = []
for r in (0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75):
    r = mp.mpf(r)
    try:
        rear, bc, smax, smin, ymin = integrals(r)
    except RuntimeError as ex:
        print("   r = %s: %s" % (mp.nstr(r, 3), ex)); continue
    L = logR2 + mp.log(r)
    rows.append((r, L, rear, bc))
    print("   %4s | %8s | %9s | %9s | %9s %9s | %6s" % (mp.nstr(r, 3), mp.nstr(L, 6), mp.nstr(rear, 6), mp.nstr(bc, 6), mp.nstr(smax, 4), mp.nstr(smin, 4), mp.nstr(ymin, 4)))
    sys.stdout.flush()
print("\nbounds m <= I/(L - tau):  systems {1, F} and {1, sq, F, F sq} have tau = 3/2 (contradiction if < 2, resp. < 4);")
print("   the 8-function system has tau = 15/8 (contradiction if < 8).  Using the rearrangement integral (an upper bound for I):")
for r, L, rear, bc in rows:
    b1 = rear/(L - mp.mpf(3)/2) if L > mp.mpf(3)/2 else mp.inf
    b2 = rear/(L - mp.mpf(15)/8) if L > mp.mpf(15)/8 else mp.inf
    print("   r = %4s: L = %6s;  I/(L - 3/2) = %7s (need < 2 or < 4);  I/(L - 15/8) = %7s (need < 8)" % (mp.nstr(r, 3), mp.nstr(L, 5), mp.nstr(b1, 5), mp.nstr(b2, 5)))
