"""Twisted analogue: v-scale pole values at X = C are -G1(a); test
   -G1(a) == sum_k (-1)^k m_k / (1 - w_a)^(k+1),   w_a = -4a^2,   m_k = phi0((v-1)^k) = sum_i C(k,i)(-1)^(k-i) mu_v(i),  mu_v(i) = (2i+1)|E_2i|/2.
Then end-to-end: A(C) = A(0) + C B (v-scale), SNF(A(C)) vs twisted h_q = 8q - 4 s2(q) - 1, transition mod 2, and the coefficient law for all M."""
import sys, math, io, contextlib
from lemlib import *
from sympy import euler
import twisted_hankel as TH
C = Cconst(200)
beta = betas(80)
G = lambda t: (-1)**abs(t)*beta[abs(t)]
G1 = lambda t: G(t) - C*(-1)**abs(t)
KM = int(sys.argv[1]) if len(sys.argv) > 1 else 200
muv = [F((2*i + 1)*abs(int(euler(2*i))), 2) for i in range(KM)]
m = [sum(math.comb(k, i)*(-1)**(k - i)*muv[i] for i in range(k + 1)) for k in range(KM)]
print('v2 of centred moments m_k, k<10:', [v2(x) for x in m[:10]])
for a in range(10):
    w = F(-4*a*a)
    S = sum((-1)**k*m[k]/(1 - w)**(k + 1) for k in range(KM))
    print('a=%d  v2(-G1(a)) = %s  v2(S) = %s  v2(diff) = %s' % (a, v2(-G1(a)), v2(S), v2(S + G1(a))))
