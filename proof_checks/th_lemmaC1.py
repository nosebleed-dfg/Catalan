"""Coefficient law (Q): q(Y) = det(A(C) + Y B) = p(Y + C), p(X) = det(A + X B) (native scale of twisted_hankel).
t_k := v2(q_k) - k  should satisfy  t_k >= T(K, K-k)  with equality at the tropical minimiser(s), and all non-minimising t_k >= m0 + 1.
T is the native-scale tropical function: twisted u-scale T^u(K,M) (content_law), Catalan v-scale (unified_snf)."""
import sys, io, contextlib, math
from lemlib import *
import twisted_hankel as TH
import twoadic_law as TL
C = Cconst()
def Tfun(mach, K):
    if mach == 'twist':
        eps = lambda n: -(2*K - 3) + 12*n - 4*s2(n)
        kap = lambda m: 1 + s2(K - 1 - m) - s2(m)
        return [sum(eps(n) for n in range(M)) + sum(kap(m) for m in range(K - M)) for M in range(K + 1)]
    c = -1
    pole = lambda m: 1 + s2(K - 1 - m) - s2(m) + 4*m - 2*K + c
    return [sum(pole(m) for m in range(K - M)) + sum(hq(mach, q) for q in range(M)) for M in range(K + 1)]
for mach in sys.argv[1].split(','):
    for K in [int(x) for x in sys.argv[2].split(',')]:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            A, B, h, labels, head, tail = TH.entries(mach, K)
            p = TH.detpoly(A, B, h)
        # q_k = sum_j p_j C(j,k) C^(j-k)
        q = [sum(p[j]*math.comb(j, k)*C**(j - k) for j in range(k, len(p))) for k in range(len(p))]
        t = [(v2(q[k]) - k) if q[k] != 0 else None for k in range(len(q))]
        T = Tfun(mach, K)
        m0 = min(T); Ms = [M for M in range(K + 1) if T[M] == m0]
        ks = [K - M for M in Ms]
        viol = [(k, t[k], T[K - k]) for k in range(K + 1) if t[k] is not None and t[k] < T[K - k]]
        eqmin = all(t[k] == m0 for k in ks)
        others = min(t[k] for k in range(K + 1) if k not in ks and t[k] is not None)
        cont = -min(v2(x) for x in p if x != 0)
        print('%s K=%d: minimisers M=%s  t at minimisers = %s (m0 = %d) %s;  min other t - m0 = %d;  t >= T violations %s;  content e2 = %d, law = %d'
              % (mach, K, Ms, [t[k] for k in ks], m0, 'OK' if eqmin else 'FAIL', others - m0, viol, cont, -m0 - (1 if len(Ms) > 1 else 0)))
        sys.stdout.flush()
