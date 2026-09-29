"""Square-root-free check of section 9, step 4 (rewritten 2026-09-29).
Node set X = {0..N}, binomial weight, monic Krawtchouk Q_n, nu(n) = v2 ||Q_n||^2.
Finite monic Jacobi matrix Jm (column convention): Jm[c+1,c] = 1, Jm[c,c] = N/2, Jm[c-1,c] = lam_c  (x Q_c = Q_{c+1} + (N/2) Q_c + lam_c Q_{c-1}, Q_{N+1} = 0 on X).
Balanced valuation of a matrix: beta(M) = min_{a,c} v2(M[a][c]) + (nu(a) - nu(c))/2.
Checks: (1) matrix of multiplication by x on X == Jm; (2) beta(C(Jm,k)) >= -k - v2(k!) (N odd) / -k/2 - v2(k!) (N even), all k <= N;
        (3) the multiplication-by-S matrix == sum_k s_k C(Jm,k) (finite Mahler); (4) its balanced (a,c)-valuations on the criterion pairs >= 1."""
import sys, math
from lemlib import *
C = Cconst()
beta_ = betas(200)
G = lambda t: (-1)**abs(t)*beta_[abs(t)]
G1 = lambda t: G(t) - C*(-1)**abs(t)
vf = lambda n: n - s2(n)
for K in [int(x) for x in sys.argv[1].split(',')]:
    for mach in ['catalan', 'twist']:
        N = 2*K - 1 if mach == 'catalan' else 2*K - 2
        ys = [F(2*x - N, 2) for x in range(N + 1)]; wts = [math.comb(N, x) for x in range(N + 1)]
        lam = [F(0)] + [F(n*(N + 1 - n), 4) for n in range(1, N + 1)]
        Q = [[F(1)]*(N + 1), ys[:]]
        for n in range(1, N): Q.append([ys[i]*Q[n][i] - lam[n]*Q[n - 1][i] for i in range(N + 1)])
        ip = lambda f, g: sum(wts[x]*f[x]*g[x] for x in range(N + 1))
        nrm = [ip(Q[n], Q[n]) for n in range(N + 1)]; nu = [v2(t) for t in nrm]
        def mult_matrix(fvals):
            return [[ip(Q[a], [fvals[x]*Q[cc][x] for x in range(N + 1)])/nrm[a] for cc in range(N + 1)] for a in range(N + 1)]
        xs = [F(x) for x in range(N + 1)]
        Jm = [[F(0)]*(N + 1) for _ in range(N + 1)]
        for cc in range(N + 1):
            if cc + 1 <= N: Jm[cc + 1][cc] = F(1)
            Jm[cc][cc] = F(N, 2)
            if cc >= 1: Jm[cc - 1][cc] = lam[cc]
        ok1 = mult_matrix(xs) == Jm
        def matmul(A_, B_): return [[sum(A_[i][k]*B_[k][j] for k in range(N + 1)) for j in range(N + 1)] for i in range(N + 1)]
        def bal(M): return min((v2(M[a][cc]) + F(nu[a] - nu[cc], 2)) for a in range(N + 1) for cc in range(N + 1) if M[a][cc] != 0)
        Ck = [[F(int(i == j)) for j in range(N + 1)] for i in range(N + 1)]; prod_ = [r[:] for r in Ck]
        bad2 = []; Cks = [Ck]
        for k in range(1, N + 1):
            prod_ = matmul(prod_, [[Jm[i][j] - (k - 1)*(1 if i == j else 0) for j in range(N + 1)] for i in range(N + 1)])
            Ck = [[t/math.factorial(k) for t in r] for r in prod_]; Cks.append(Ck)
            bound = -k - vf(k) if N % 2 else F(-k, 2) - vf(k)
            if bal(Ck) < bound: bad2.append(k)
        if mach == 'catalan':
            Sx = lambda x: (G1(x - K + 1) - G1(x - K))/(2*C)
            pairs = [(a, cc) for a in range(1, N + 1, 2) for cc in range(0, N + 1, 2)]
        else:
            Sx = lambda x: G1(x - K + 1)/C
            kk = K//2
            pairs = [(a, N - bb) for a in range(0, 2*kk - 1, 2) for bb in range(0, 2*kk - 1, 2)]
        Sv = [Sx(x) for x in range(N + 1)]
        s = [sum((-1)**(k - i)*math.comb(k, i)*Sv[i] for i in range(k + 1)) for k in range(N + 1)]
        MS = [[sum(s[k]*Cks[k][a][cc] for k in range(N + 1)) for cc in range(N + 1)] for a in range(N + 1)]
        ok3 = MS == mult_matrix(Sv)
        marg = min(v2(MS[a][cc]) + F(nu[a] - nu[cc], 2) for (a, cc) in pairs if MS[a][cc] != 0)
        print('%s K=%d N=%d: (1) mult-by-x == Jm %s  (2) balanced bound violations %s  (3) M^S == sum s_k C(Jm,k) %s  (4) min balanced margin on criterion pairs = %s'
              % (mach, K, N, ok1, bad2, ok3, marg))
        sys.stdout.flush()
