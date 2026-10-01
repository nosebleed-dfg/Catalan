r"""radius_transfer.py (2026-10-01, real date): S(w) = sum of 1/c^2 over the double cosets of <T, L_w>, to many digits,
by a transfer operator instead of word enumeration (uniformize_radius.py converges like X^(2 delta - 2)).

Words B^{m1} A^{k1} ... B^{mj} (A = T, B = L_w) have bottom rows built by (c, d) -> (c + d w m, d) and (c, d) -> (c, d + c k).
With G(c, d) = d^-2 g(c/d) the total weight of all continuations of a prefix,
    g(x)   = sum_{m != 0} v^2 (1 + Psi(v)),      v = 1/(x + w m),
    Psi(u) = sum_{k != 0} (k + u)^-2 g(1/(k + u)),
and S(w) = g(0).  Psi lives on [-rho, rho], g on [-z, z], with rho = (1 - sqrt(1 - 4/w))/2 and z = 1/(1 - rho); both are
analytic well beyond these intervals, so polynomial collocation at Chebyshev nodes converges geometrically.  For a
polynomial the sums over m and k are exact Hurwitz zeta values:
    sum_{m != 0} (x + w m)^-n = w^-n [zeta(n, 1 + x/w) + (-1)^n zeta(n, 1 - x/w)].
So the fixed point is one linear solve.  log R* = 2 pi^2 S(w).
Usage: python radius_transfer.py [N] [w1 w2 ...]"""
import sys
import mpmath as mp

mp.mp.dps = 70

def S_of_w(w, N):
    w = mp.mpf(w)
    rho = (1 - mp.sqrt(1 - 4/w))/2; z = 1/(1 - rho)
    s = [mp.cos(mp.pi*(2*i + 1)/(2*N)) for i in range(N)]            # Chebyshev nodes in the scaled variable
    V = mp.matrix(N, N)
    for i in range(N):
        for j in range(N): V[i, j] = s[i]**j
    Vinv = V**-1
    HA = mp.matrix(N, N); HB = mp.matrix(N, N)
    for i in range(N):
        x = z*s[i]; u = rho*s[i]
        for j in range(N):
            n = 2 + j; sg = (-1)**j
            HA[i, j] = (mp.zeta(n, 1 + x/w) + sg*mp.zeta(n, 1 - x/w))/(w**n*rho**j)      # psi coefficient j (scaled) -> g(x_i)
            HB[i, j] = (mp.zeta(n, 1 + u) + sg*mp.zeta(n, 1 - u))/z**j                   # g coefficient j (scaled) -> Psi(u_i)
    hA0 = mp.matrix([HA[i, 0] for i in range(N)])                                         # the "1 +" term
    MA = Vinv*HA*Vinv                      # Psi values -> g coefficients (linear part)
    cA = Vinv*hA0
    A = HB*MA; b = HB*cA
    I = mp.eye(N)
    psi = mp.lu_solve(I - A, b)
    gc = MA*psi + cA
    return gc[0], rho

if __name__ == '__main__':
    args = [a for a in sys.argv[1:]]
    Ns = [int(args[0])] if args else [28, 36, 44]
    ws = [int(a) for a in args[1:]] or [5, 6, 7, 8, 9]
    res = {}
    for w in ws:
        vals = [S_of_w(w, N)[0] for N in Ns]
        res[w] = vals[-1]
        print("w = %d: S = %s   (change from N = %s: %s)   log R* = %s"
              % (w, mp.nstr(vals[-1], 25), Ns[:-1], [mp.nstr(abs(v - vals[-1]), 3) for v in vals[:-1]], mp.nstr(2*mp.pi**2*vals[-1], 20)))
        sys.stdout.flush()
    if all(w in res for w in (5, 6, 8, 9)):
        tot = sum(res[w] for w in (5, 6, 8, 9))
        print("sum over w = 5, 6, 8, 9: S = %s,  2 pi^2 S = %s;   pi^2 = %s,  10 = 10" % (mp.nstr(tot, 20), mp.nstr(2*mp.pi**2*tot, 20), mp.nstr(mp.pi**2, 12)))
