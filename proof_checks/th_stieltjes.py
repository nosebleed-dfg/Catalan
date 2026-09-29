"""Test: Catalan pole values at X = C equal the 2-adic Stieltjes transform of the CDH moments at the nodes.
F_j(C) = (2j+1) DeltaG1(j) / 4     vs     S(w_j) = -sum_k mu_k w_j^(-k-1),  w_j = -(2j+1)^2,  mu_k = 4^k (2^(2k+2)-1)|B_(2k+2)| (v-scale).
Also the twisted analogue in u-scale for the record (moments do not decay there)."""
import sys
from lemlib import *
from sympy import bernoulli
C = Cconst(200)
beta = betas(80)
G = lambda t: (-1)**abs(t)*beta[abs(t)]
G1 = lambda t: G(t) - C*(-1)**abs(t)
KM = int(sys.argv[1]) if len(sys.argv) > 1 else 90
mu = []
for k in range(KM):
    B = bernoulli(2*k + 2); mu.append(F(4)**k*(2**(2*k + 2) - 1)*abs(F(int(B.p), int(B.q))))
print('v2(mu_k), k<8:', [v2(m) for m in mu[:8]])
for j in range(12):
    w = F(-(2*j + 1)**2)
    S = -sum(mu[k]*w**(-k - 1) for k in range(KM))
    Fj = F(2*j + 1, 4)*(G1(j + 1) - G1(j))
    print('j=%2d  v2(F_j(C)) = %s   v2(S(w_j)) = %s   v2(difference) = %s' % (j, v2(Fj), v2(S), v2(S - Fj)))
