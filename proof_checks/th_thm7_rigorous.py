"""Square-root-free check of Theorem 7 (twoadic_law_proof.md section 12, rewritten 2026-09-29).
Monic CDH Jacobi matrix J (column convention): J[n+1,n] = 1, J[n,n] = b_n, J[n-1,n] = lam_n, so v^k P_q = sum_p (J^k)[p,q] P_p.
Rescaled: Jhat = S^-1 (J - v* I) S, S = diag(rho^-n), rho = 4, v* = 0 (Catalan) / rho = 2, v* = 1 (twisted).
Checks, all in Q_2 (exact rationals):
 (1) Jhat entries are integers divisible by rho (on a large index window);
 (2) G[p][q] := phi0^(2)(P_p P_q pi) == h_p rho^(q-p) Pi[p][q], Pi = sum_k c_k Jhat^k (computed independently);
 (3) Gaussian elimination on Pi_[K]: L0 == I, U0 == c0 I (mod rho); L := H L0 H^-1 integral with v2(L[p][q]) >= 3(p-q) + (2 Catalan | 1 twisted);
 (4) G == L diag(Lam) L^T and v2(Lam_p) == v2(h_p)."""
import sys, math
from lemlib import *
from sympy import bernoulli, euler
mach = sys.argv[1]; K = int(sys.argv[2]); KM = int(sys.argv[3]) if len(sys.argv) > 3 else 8*K + 80
if mach == 'catalan':
    rho, vstar = 4, 0
    b = lambda n: 4*(2*n + 1)*(n + 1); lam = lambda n: 16*n**3*(n + 1)
    mu = []
    for k in range(KM + 2*K + 2):
        B_ = bernoulli(2*k + 2); mu.append(F(4)**k*(2**(2*k + 2) - 1)*abs(F(int(B_.p), int(B_.q))))
    nodes = [-(2*j + 1)**2 for j in range(K)]
else:
    rho, vstar = 2, 1
    b = lambda n: 8*n*n + 8*n + 3; lam = lambda n: 16*n**4
    mu = [F((2*i + 1)*abs(int(euler(2*i))), 2) for i in range(KM + 2*K + 2)]
    nodes = [-4*a*a for a in range(K)]
h = [F(1, 2)]
for n in range(1, K + 1): h.append(h[-1]*lam(n))
W = K + KM + 2                                   # index window for the infinite matrices
Jh = [[F(0)]*W for _ in range(W)]
for n in range(W):
    if n + 1 < W: Jh[n + 1][n] = F(rho)                                    # rho^(n+1) * 1 * rho^-n
    Jh[n][n] = F(b(n) - vstar)
    if n >= 1: Jh[n - 1][n] = F(lam(n), rho)                               # rho^(n-1) * lam_n * rho^-n
ok1 = all(x.denominator == 1 and x.numerator % rho == 0 for r in Jh for x in r)
# pi = 1/N_K expanded at v*: N_K(v* + e) = prod (e + v* - w)
Np = [F(1)]
for w in nodes:
    a0 = vstar - w; new = [F(0)]*(len(Np) + 1)
    for i, cc in enumerate(Np): new[i] += a0*cc; new[i + 1] += cc
    Np = new
c = [F(1)/Np[0]]
for k in range(1, KM):
    c.append(-sum(Np[i]*c[k - i] for i in range(1, min(k, K) + 1))/Np[0])
# Pi = sum_k c_k Jhat^k, leading K x K block (vectors e_q propagated)
Pi = [[F(0)]*K for _ in range(K)]
for q in range(K):
    vec = [F(0)]*W; vec[q] = F(1)
    for k in range(KM):
        for p in range(K): Pi[p][q] += c[k]*vec[p]
        vec = [sum(Jh[i][j]*vec[j] for j in (i - 1, i, i + 1) if 0 <= j < W) for i in range(W)]
# centred moments m_k = phi0((v - v*)^k)
m = [sum(math.comb(k, i)*(-vstar)**(k - i)*mu[i] for i in range(k + 1)) for k in range(KM + 2*K)]
# CDH monic OPs, then G directly from the definition phi0^(2)(P_p P_q pi)
Pc, hc = monic_ops(mu[:2*K + 1], K)
def polymul(p, q):
    r = [F(0)]*(len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): r[i + j] += x*y
    return r
def recentre(p):        # coefficients in (v - v*)
    out = [F(0)]*len(p)
    for i, a in enumerate(p):
        for j in range(i + 1): out[j] += a*math.comb(i, j)*F(vstar)**(i - j)
    return out
G = [[F(0)]*K for _ in range(K)]
for p in range(K):
    for q in range(K):
        f = recentre(polymul(Pc[p], Pc[q]))
        G[p][q] = sum(f[i]*c[k]*m[i + k] for i in range(len(f)) for k in range(KM - len(f)))
ok2 = min(v2(G[p][q] - h[p]*F(rho)**(q - p)*Pi[p][q]) or 10**6 for p in range(K) for q in range(K))
# (3) LU of Pi
Mx = [r[:] for r in Pi]; L0 = [[F(int(i == j)) for j in range(K)] for i in range(K)]
for t in range(K):
    for i in range(t + 1, K):
        l = Mx[i][t]/Mx[t][t]; L0[i][t] = l
        for j in range(t, K): Mx[i][j] -= l*Mx[t][j]
U0 = Mx
c0 = c[0]
ok3a = all(v2(L0[i][j]) is None or v2(L0[i][j]) >= (2 if rho == 4 else 1) for i in range(K) for j in range(i))
ok3b = all(v2(U0[i][i] - c0) >= (2 if rho == 4 else 1) for i in range(K)) and all(v2(U0[i][j]) is None or v2(U0[i][j]) >= (2 if rho == 4 else 1) for i in range(K) for j in range(i + 1, K))
Hd = [h[p]*F(rho)**(-p) for p in range(K)]
L = [[Hd[i]*L0[i][j]/Hd[j] for j in range(K)] for i in range(K)]
extra = 2 if rho == 4 else 1
ok3c = all(v2(L[i][j]) is None or v2(L[i][j]) >= 3*(i - j) + extra for i in range(K) for j in range(i))
Lam = [h[p]*U0[p][p] for p in range(K)]
LLt = [[sum(L[i][k]*Lam[k]*L[j][k] for k in range(K)) for j in range(K)] for i in range(K)]
ok4 = min(v2(LLt[i][j] - G[i][j]) or 10**6 for i in range(K) for j in range(K))
ok4b = all(v2(Lam[p]) == v2(h[p]) for p in range(K))
print('%s K=%d KM=%d: (1) Jhat integral, divisible by %d: %s  (2) min v2(G - h rho^(q-p) Pi) = %s' % (mach, K, KM, rho, ok1, ok2))
print('   (3) L0 == I, U0 == c0 I mod %d: %s %s;  L integral with v2(L_pq) >= 3(p-q)+%d: %s' % (rho, ok3a, ok3b, extra, ok3c))
print('   (4) min v2(L Lam L^T - G) = %s (precision of truncation);  v2(Lam_p) == v2(h_p): %s' % (ok4, ok4b))
