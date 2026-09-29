"""Verify each link of the Lemma A proof (Krawtchouk mirror + smoothness).
Binomial on x = 0..N (N' = N), monic Krawtchouk Q_n in y = x - N/2 (recurrence y Q_n = Q_{n+1} + lam_n Q_{n-1}, lam_n = n(N+1-n)/4).
 (i)   mirror: (-1)^x Q_b(y) = (-1)^N 2^(N-2b) b!/(N-b)! Q_{N-b}(y) at every node
 (ii)  nu(n) = v2 ||Q_n||^2 = N - 2n + v2(n!) + v2(N!) - v2((N-n)!)
 (iii) Jacobi entries: v2(lam_n) >= -2 (N odd), >= -1 (N even)
 (iv)  Mahler: S(x) = DeltaG1(x-K) (Catalan, N = 2K-1) has v2(Delta^k S(0)) >= k + v2((k+2)!);
       S(x) = G1(x-(K-1))/C (twisted, N = 2K-2) has v2(Delta^k S(0)) >= k + v2((k+1)!) for k >= 1, S(0)... (value at y=0 is -1)
 (v)   orthonormal multiplication matrix: min over a odd, c even (Catalan) of v2<Q_a,S Q_c> - (nu(a)+nu(c))/2 >= 1;
       twisted: a even <= 2k-2 (k = floor(K/2)), c = N - b with b even <= 2k-2: margin >= 1."""
import sys, math
from lemlib import *
C = Cconst()
beta = betas(400)
G = lambda t: (-1)**abs(t)*beta[abs(t)]            # even extension on Z
G1 = lambda t: G(t) - C*(-1)**abs(t)                # smooth part
def kraw(N):
    ys = [F(2*x - N, 2) for x in range(N + 1)]
    lam = [F(0)] + [F(n*(N + 1 - n), 4) for n in range(1, N + 1)]
    Q = [[F(1)]*(N + 1), ys[:]]
    for n in range(1, N):
        Q.append([ys[i]*Q[n][i] - lam[n]*Q[n - 1][i] for i in range(N + 1)])
    return ys, lam, Q
def fd(h, k, x0): return sum((-1)**(k - i)*math.comb(k, i)*h(x0 + i) for i in range(k + 1))
for K in [int(x) for x in sys.argv[1].split(',')]:
    for mach in ['catalan', 'twist']:
        N = 2*K - 1 if mach == 'catalan' else 2*K - 2
        ys, lam, Q = kraw(N)
        b = [math.comb(N, x) for x in range(N + 1)]
        ip = lambda f, g: sum(b[x]*f[x]*g[x] for x in range(N + 1))
        # (i)
        bad_i = sum(1 for bb in range(N + 1) for x in range(N + 1)
                    if (-1)**x*Q[bb][x] != (-1)**N*F(2)**(N - 2*bb)*F(math.factorial(bb), math.factorial(N - bb))*Q[N - bb][x])
        # (ii)
        nu = [v2(ip(Q[n], Q[n])) for n in range(N + 1)]
        vf = lambda n: n - s2(n)
        bad_ii = sum(1 for n in range(N + 1) if nu[n] != N - 2*n + vf(n) + vf(N) - vf(N - n))
        # (iii)
        mlam = min(v2(lam[n]) for n in range(1, N + 1))
        # (iv), (v)
        if mach == 'catalan':
            Sx = lambda x: G1(x - K + 1) - G1(x - K)          # DeltaG1(x-K)
            need = lambda k: k + vf(k + 2)
            pairs = [(a, c) for a in range(1, N + 1, 2) for c in range(0, N + 1, 2)]
        else:
            Sx = lambda x: G1(x - (K - 1))/C
            need = lambda k: k + vf(k + 1)
            kk = K//2
            pairs = [(a, N - bb) for a in range(0, 2*kk - 1, 2) for bb in range(0, 2*kk - 1, 2)]
        worst_iv = min(v2(fd(Sx, k, 0)) - need(k) for k in range(1, N + 1) if fd(Sx, k, 0) != 0)
        Svals = [Sx(x) for x in range(N + 1)]
        marg = min(v2(ip(Q[a], [Svals[x]*Q[c][x] for x in range(N + 1)])) - F(nu[a] + nu[c], 2)
                   for (a, c) in pairs if ip(Q[a], [Svals[x]*Q[c][x] for x in range(N + 1)]) != 0)
        print('%s K=%d N=%d: (i) mirror violations %d  (ii) norm-formula violations %d  (iii) min v2(lam) = %d  (iv) min v2(D^k S(0)) - bound = %d  (v) min margin = %s'
              % (mach, K, N, bad_i, bad_ii, mlam, worst_iv, marg))
        sys.stdout.flush()
