"""Catalan pole measure = Krawtchouk base omega0 + perturbation. In the omega0 monic-OP basis (v-scale, 2-integral):
Gram = diag(h0) + E,  E_mn = sum_j (omega_j - omega0_j) p_m(w_j) p_n(w_j).
Report margin(m,n) = v2(E_mn) - (b_m + b_n)/2 with b = band; (star) holds iff all margins > 0.
Also check: p_n in Z_2[v] (2-integrality of the v-scale Krawtchouk OPs)."""
import sys, math
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from fractions import Fraction as F
from th_lemmaA2 import v2, s2
def Cconst(nmax):
    C = F(0)
    for n in range(nmax + 1):
        d = sum((-1)**(n - i)*math.comb(n, i)*F(1, (2*i + 1)**2) for i in range(n + 1))
        C += F((-1)**n, 2**(n + 1))*d
    return C
C = Cconst(int(sys.argv[2]) if len(sys.argv) > 2 else 160)
beta = [F(0)]
for k in range(400): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
def polys(nodes, w, nmax):
    """monic OPs as coefficient lists (low->high) plus values at nodes"""
    K = len(nodes)
    P = [[F(1)]]; vals = [[F(1)]*K]; h = []
    prevc = [F(0)]; prevv = [F(0)]*K
    for n in range(nmax):
        hn = sum(w[i]*vals[-1][i]**2 for i in range(K)); h.append(hn)
        if n == nmax - 1: break
        bn = sum(w[i]*nodes[i]*vals[-1][i]**2 for i in range(K))/hn
        lam = hn/h[-2] if n > 0 else F(0)
        cur = P[-1]
        newc = [F(0)] + cur
        for i in range(len(cur)): newc[i] -= bn*cur[i]
        for i in range(len(prevc)): newc[i] -= lam*prevc[i]
        newv = [(nodes[i] - bn)*vals[-1][i] - lam*prevv[i] for i in range(K)]
        prevc, prevv = cur, vals[-1]
        P.append(newc); vals.append(newv)
    return P, vals, h
for K in [int(x) for x in sys.argv[1].split(',')]:
    nodes = [F(-(2*j + 1)**2) for j in range(K)]
    DG = lambda j: (-1)**(j + 1)*beta[j + 1] - (-1)**j*beta[j]
    den = lambda j: 4**K*math.factorial(K - 1 - j)*math.factorial(K + j)
    w = [F((2*j + 1)**2)*(-1)**j*DG(j)/den(j) for j in range(K)]
    w0 = [F((2*j + 1)**2)*(-2*C)/den(j) for j in range(K)]
    P, vals, h0 = polys(nodes, w0, K)
    integral = all(v2(c) is None or v2(c) >= 0 for p in P for c in p)
    b = [1 + s2(K - 1 - n) - s2(n) + 4*n - 2*K - 1 for n in range(K)]
    okh = all(v2(h0[n]) == b[n] for n in range(K))
    dw = [w[j] - w0[j] for j in range(K)]
    marg = {}
    for m in range(K):
        for n in range(m, K):
            e = sum(dw[j]*vals[m][j]*vals[n][j] for j in range(K))
            marg[(m, n)] = (v2(e) if e != 0 else 10**6) - F(b[m] + b[n], 2)
    worst = min(marg.values()); arg = min(marg, key=marg.get)
    diag = [float(marg[(n, n)]) for n in range(K)]
    print('K=%d: OPs 2-integral: %s; v2(h0)=band: %s; min margin %s at %s; diag margins %s' % (K, integral, okh, worst, arg, diag))
    sys.stdout.flush()
