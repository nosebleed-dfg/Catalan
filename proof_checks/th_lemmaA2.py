"""Lemma A reformulated: the discrete pole measure omega_K = sum_a F(a)/D'(z_a) delta_{z_a}.
Compute its monic-OP norms h_n (n < K) exactly (Stieltjes on K nodes) and compare v2(h_n) with
  twisted (u-scale, nodes -a^2):        kappa(K,n) = 1 + s2(K-1-n) - s2(n)
  Catalan (v-scale, nodes -(2j+1)^2):   printed, to find its law.
Also the pure Krawtchouk parts: twisted weight C(2K-2,K-1-a)*eps_a, Catalan weight C(2K-1,K-1-j)*(2j+1)^2."""
from fractions import Fraction as F
import math, sys
s2 = lambda n: bin(n).count('1')
def v2(x):
    x = F(x); a, b = x.numerator, x.denominator
    if a == 0: return None
    return ((abs(a) & -abs(a)).bit_length() - 1) - ((b & -b).bit_length() - 1)
def norms(nodes, w, nmax):
    K = len(nodes)
    p_prev = [F(0)]*K; p = [F(1)]*K; h = []; lam = F(0)
    for n in range(nmax):
        hn = sum(w[i]*p[i]*p[i] for i in range(K))
        if hn == 0: break
        h.append(hn)
        bn = sum(w[i]*nodes[i]*p[i]*p[i] for i in range(K))/hn
        lam = hn/h[-2] if n > 0 else F(0)
        p_new = [(nodes[i] - bn)*p[i] - lam*p_prev[i] for i in range(K)]
        p_prev, p = p, p_new
    return h
def twist_measure(K):
    beta = [F(0)]
    for k in range(K): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
    nodes = [F(-a*a) for a in range(K)]
    w = []
    for a in range(K):
        Fa = F(-4*(-1)**a)*beta[a]
        Dp = math.prod(b*b - a*a for b in range(K) if b != a)   # D'(z_a) with D(u) = prod(u + b^2) -> D'(-a^2) = prod_{b != a}(b^2 - a^2)
        w.append(Fa/Dp)
    return nodes, w
def cat_measure(K):
    beta = [F(0)]
    for k in range(K): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
    nodes = [F(-(2*j + 1)**2) for j in range(K)]
    w = []
    for j in range(K):
        b = 2*j + 1; sign = (-1)**j
        Fb = F(-sign*b, 2)*beta[j] - F(1, 4*b)
        Dp = math.prod((2*i + 1)**2 - b*b for i in range(K) if i != j)
        w.append(Fb/Dp)
    return nodes, w
if __name__ == '__main__':
    Ks = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [8, 12, 16, 20, 24, 32]
    for K in Ks:
        nodes, w = twist_measure(K)
        h = norms(nodes, w, K)
        vs = [v2(x) for x in h]
        kap = [1 + s2(K - 1 - n) - s2(n) for n in range(K)]
        bad = [n for n in range(len(vs)) if vs[n] != kap[n]]
        # pure Krawtchouk part
        wk = [F(math.comb(2*K - 2, K - 1 - a)*(1 if a == 0 else 2)) for a in range(K)]
        hk = norms(nodes, wk, K)
        vk = [v2(x) for x in hk]
        print('twist K=%d: v2(h_n) %s' % (K, vs))
        print('          kappa   %s   first mismatch n = %s (of %d)' % (kap, bad[0] if bad else None, K))
        print('          Kraw    %s' % vk)
        nodes, w = cat_measure(K)
        h = norms(nodes, w, K)
        vs = [v2(x) for x in h]
        wk = [F(math.comb(2*K - 1, K - 1 - j)*(2*j + 1)**2) for j in range(K)]
        hk = norms(nodes, wk, K)
        vk = [v2(x) for x in hk]
        print('cat   K=%d: v2(h_n) %s' % (K, vs))
        print('          Kraw    %s   diff %s' % (vk, [vs[i] - vk[i] for i in range(min(len(vs), len(vk)))]))
        sys.stdout.flush()
