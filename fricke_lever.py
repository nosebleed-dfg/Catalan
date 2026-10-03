r"""fricke_lever.py (2026-10-02): the Fricke lever for G.

The harmless group of case E is Gamma' = <T, L_8>, with conformal radius log R' = 2 pi^2 S(8) = 1.2323 at the cusp
(radius_transfer.py).  The Fricke involution w_8 = (0 -1; 8 0) normalises Gamma' (w_8 T w_8^-1 = L_8^-1) and swaps the two
harmless cusps oo and 0.  If a family's F = f (g - Lambda) is invariant under w_8 as well, F lives on the bigger quotient
X'' = H/<T, w_8> (the Hecke group H(sqrt 8)), whose conformal radius is
    log R'' = 2 pi^2 (S_even + S_odd),   S_even = S(8),   S_odd = (1/8) sum 1/d^2 over <T>\Gamma'/<L_8>
(the odd double cosets are gamma w_8, with lower-left entry 8 d(gamma) before normalising the determinant).
Part A computes S_odd by the transfer operator of radius_transfer.py with the constant term moved:
    h(x) = sum_{m != 0} v^2 Phi(v),  v = 1/(x + w m),   Phi(u) = sum_{k != 0} (k + u)^-2 (1 + h(1/(k + u))),
    S_odd = (1 + h(0))/w.
Part B builds a Fricke-eigen family: f_eps = theta^2 -+ sqrt2 theta(2 tau)^2 (eigenvalues +-i under |_1 w_8), and W in the
4-dimensional space E^{chi,1}(tau), E^{chi,1}(2 tau), E^{1,chi}(tau), E^{1,chi}(2 tau) (chi = chi_-4) that is a w_8-eigenform
with zero constant term.  With Lambda = L(W, 2) the function F = f (g - Lambda), g = sum c(m) m^-2 q^m, should satisfy
F(tau/(8 tau + 1)) = F(tau) (cusp 0 harmless) and F(-1/(8 tau)) = eps_f eps_W F(tau).  Case E is run as a control.
Usage: python fricke_lever.py"""
import sys
import mpmath as mp

mp.mp.dps = 60

# ---------- Part A: the radii ----------
def chebyshev_setup(w, N):
    w = mp.mpf(w)
    rho = (1 - mp.sqrt(1 - 4/w))/2; z = 1/(1 - rho)
    s = [mp.cos(mp.pi*(2*i + 1)/(2*N)) for i in range(N)]
    V = mp.matrix(N, N)
    for i in range(N):
        for j in range(N): V[i, j] = s[i]**j
    Vinv = V**-1
    HA = mp.matrix(N, N); HB = mp.matrix(N, N)
    for i in range(N):
        x = z*s[i]; u = rho*s[i]
        for j in range(N):
            n = 2 + j; sg = (-1)**j
            HA[i, j] = (mp.zeta(n, 1 + x/w) + sg*mp.zeta(n, 1 - x/w))/(w**n*rho**j)
            HB[i, j] = (mp.zeta(n, 1 + u) + sg*mp.zeta(n, 1 - u))/z**j
    return w, Vinv, HA, HB
def S_even(w, N):
    w, Vinv, HA, HB = chebyshev_setup(w, N)
    hA0 = mp.matrix([HA[i, 0] for i in range(N)])
    MA = Vinv*HA*Vinv; cA = Vinv*hA0
    psi = mp.lu_solve(mp.eye(N) - HB*MA, HB*cA)
    return (MA*psi + cA)[0]
def S_odd(w, N):
    w, Vinv, HA, HB = chebyshev_setup(w, N)
    hB0 = mp.matrix([HB[i, 0] for i in range(N)])
    MA = Vinv*HA*Vinv
    phi = mp.lu_solve(mp.eye(N) - HB*MA, hB0)       # Phi values at the nodes
    h0 = (MA*phi)[0]
    return (1 + h0)/w

def enumerate_odd(w, X):
    """direct check: (1/w) sum 1/d^2 over words B^m1 A^k1 ... B^mj A^kj with |entries| <= X (plus the identity)."""
    tot = mp.mpf(1); stack = [(0, 1)]                 # (c, d) after an even number of steps, before a B-step
    while stack:
        c, d = stack.pop()
        m = 1
        while True:
            hit = False
            for mm in (m, -m):
                c1 = c + d*w*mm
                if abs(c1) > X: continue
                hit = True
                k = 1
                while True:
                    hit2 = False
                    for kk in (k, -k):
                        d1 = d + c1*kk
                        if abs(d1) > X: continue
                        hit2 = True; tot += mp.mpf(1)/d1**2; stack.append((c1, d1))
                    if not hit2: break
                    k += 1
            if not hit: break
            m += 1
    return tot/w

print("A. radii")
for w in (5, 6, 8):
    se = [S_even(w, N) for N in (28, 40)]; so = [S_odd(w, N) for N in (28, 40)]
    print("   w = %d: S_even = %s (N-change %s), S_odd = %s (N-change %s)" % (w, mp.nstr(se[1], 20), mp.nstr(abs(se[0] - se[1]), 2), mp.nstr(so[1], 20), mp.nstr(abs(so[0] - so[1]), 2)))
    print("        log R' = 2 pi^2 S_even = %s;  log R'' = 2 pi^2 (S_even + S_odd) = %s;  gain = %s"
          % (mp.nstr(2*mp.pi**2*se[1], 15), mp.nstr(2*mp.pi**2*(se[1] + so[1]), 15), mp.nstr(2*mp.pi**2*so[1], 15)))
    if w == 8:
        for X in (200, 2000):
            print("        direct enumeration of S_odd with entries <= %d: %s" % (X, mp.nstr(enumerate_odd(w, X), 12)))
sys.stdout.flush()

# ---------- Part B: the Fricke-eigen family ----------
mp.mp.dps = 50
NQ = 700
def chi4(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]
c_E1 = [0] + [sum(chi4(n//d)*d*d for d in divisors(n)) for n in range(1, NQ + 1)]      # E^{chi,1}: no constant term
c_E2 = [mp.mpf(-1)/4] + [sum(chi4(d)*d*d for d in divisors(n)) for n in range(1, NQ + 1)]  # E^{1,chi}: constant -1/4
def qser(coefs, tau, step=1):
    q = mp.exp(2j*mp.pi*tau*step)
    return sum(c*q**n for n, c in enumerate(coefs))
def theta(tau):
    q = mp.exp(2j*mp.pi*tau)
    return 1 + 2*sum(q**(n*n) for n in range(1, 60))
W8 = lambda tau: -1/(8*tau)
L8 = lambda tau: tau/(8*tau + 1)
basis3 = [lambda t: qser(c_E1, t), lambda t: qser(c_E1, t, 2), lambda t: qser(c_E2, t), lambda t: qser(c_E2, t, 2)]
basis1 = [lambda t: theta(t)**2, lambda t: theta(2*t)**2]
def fricke_matrix(basis, k, pts):
    """M with (E_i |_k w_8)(tau) = sum_j M[i, j] E_j(tau), fitted on the points; returns M and the fit residual."""
    n = len(basis)
    A = mp.matrix(len(pts), n); Bm = mp.matrix(len(pts), n)
    for r, tau in enumerate(pts):
        for j in range(n):
            A[r, j] = basis[j](tau)
            Bm[r, j] = mp.mpf(8)**(mp.mpf(k)/2)*(8*tau)**(-k)*basis[j](W8(tau))
    M = mp.matrix(n, n); res = mp.mpf(0)
    for i in range(n):
        rhs = mp.matrix([Bm[r, i] for r in range(len(pts))])
        sol = mp.lu_solve(A.T*A, A.T*rhs)                  # least squares
        for j in range(n): M[i, j] = sol[j]
        res = max(res, mp.norm(A*sol - rhs))
    return M, res
pts = [mp.mpc(x, y) for x, y in ((0.05, 0.36), (-0.17, 0.34), (0.23, 0.37), (0.41, 0.35), (-0.33, 0.38), (0.11, 0.33))]
M3, r3 = fricke_matrix(basis3, 3, pts)
M1, r1 = fricke_matrix(basis1, 1, pts)
print("\nB. the Fricke action")
print("   weight 3 (basis E1(tau), E1(2tau), E2(tau), E2(2tau)): fit residual %s" % mp.nstr(r3, 3))
for i in range(4): print("      " + "  ".join(mp.nstr(M3[i, j], 12) for j in range(4)))
print("   weight 1 (basis theta^2, theta(2tau)^2): fit residual %s" % mp.nstr(r1, 3))
for i in range(2): print("      " + "  ".join(mp.nstr(M1[i, j], 12) for j in range(2)))
print("   M3^2 = -I: %s;  M1^2 = -I: %s" % (mp.nstr(mp.norm(M3*M3 + mp.eye(4)), 3), mp.nstr(mp.norm(M1*M1 + mp.eye(2)), 3)))
G = mp.catalan
def eichler(coefs):
    return [mp.mpf(0)] + [coefs[m]/mp.mpf(m)**2 for m in range(1, len(coefs))]
def family(fvec, wvec, Lam):
    # W = a1 E1(tau) + a2 E1(2tau) + a3 E2(tau) + a4 E2(2tau) as a q-series
    cW = [mp.mpf(0)]*(NQ + 1)
    for n in range(NQ + 1):
        v = wvec[0]*c_E1[n] + wvec[2]*c_E2[n]
        if n % 2 == 0: v += wvec[1]*c_E1[n//2] + wvec[3]*c_E2[n//2]
        cW[n] = v
    cg = eichler(cW)
    def F(tau):
        f = fvec[0]*theta(tau)**2 + fvec[1]*theta(2*tau)**2
        return f*(qser(cg, tau) - Lam)
    return F, cW
test = [mp.mpc(0.07, 0.31), mp.mpc(-0.21, 0.29), mp.mpc(0.37, 0.33)]
def report(name, F):
    dL = max(abs(F(L8(t)) - F(t)) for t in test); dW = [F(W8(t))/F(t) for t in test]
    print("   %s: |F(L8 tau) - F(tau)| <= %s;  F(w8 tau)/F(tau) = %s" % (name, mp.nstr(dL, 3), [mp.nstr(r, 12) for r in dW]))
FE, _ = family([1, 0], [1, -8, 0, 0], G/2)
report("case E (f = theta^2, W = E1 - 8 E1(2tau), Lambda = G/2)", FE)
# eigenvectors
ev3, EV3 = mp.eig(M3.T)        # right eigenvectors of M3^T: coefficient vectors a with (sum a_i E_i)|w8 = eps sum a_i E_i
ev1, EV1 = mp.eig(M1.T)
print("   eigenvalues of w_8 on weight 3: %s" % [mp.nstr(e, 8) for e in ev3])
print("   eigenvalues of w_8 on weight 1: %s" % [mp.nstr(e, 8) for e in ev1])
fam = {}
for sgn in (1, -1):
    eps = mp.mpc(0, sgn)
    idx = [i for i in range(4) if abs(ev3[i] - eps) < 1e-20]
    vecs = [mp.matrix([EV3[j, i] for j in range(4)]) for i in idx]
    # combine the two eigenvectors to kill the constant term a3 + a4
    v1, v2 = vecs
    s1, s2 = v1[2] + v1[3], v2[2] + v2[3]
    W = v1*s2 - v2*s1
    W = W/W[0]
    a = [W[j] for j in range(4)]
    Lam = -(a[0] + a[1]/4)*G/2 + a[2]*mp.pi**2/16
    fi = [i for i in range(2) if abs(ev1[i] - mp.conj(eps)) < 1e-20][0]
    fv = [EV1[j, fi] for j in range(2)]; fv = [x/fv[0] for x in fv]
    print("   eps_W = %s: W = E1 + (%s) E1(2tau) + (%s) (E2 - E2(2tau));  f = theta^2 + (%s) theta(2tau)^2 (eps_f = %s)"
          % (mp.nstr(eps, 3), mp.nstr(a[1], 12), mp.nstr(a[2], 12), mp.nstr(fv[1], 12), mp.nstr(ev1[fi], 3)))
    print("      a2/a3 = %s,  a3^2 = %s,  fv1^2 = %s;  Lambda = L(W, 2) = %s" % (mp.nstr(a[1]/a[2], 12), mp.nstr(a[2]**2, 12), mp.nstr(fv[1]**2, 12), mp.nstr(Lam, 15)))
    F, cW = family(fv, a, Lam)
    report("eigen family eps_W = %s" % mp.nstr(eps, 2), F)
    fam[sgn] = (fv, a, Lam, cW)
