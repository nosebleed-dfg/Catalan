"""Twisted structure theorem, end-to-end (v-scale, v = 4u, nodes w_a = -4a^2).
 A(C)_s = phi^(2)( v^s / N_K(v) ), phi^(2) = expansion around v = 1 with centred moments m_k = phi0((v-1)^k)   -- no C needed.
 (1) A(C) == A^v(0) + C B^v      (2) C == phi^(2)(1/v) = sum (-1)^k m_k
 (3) SNF(A(C)) == {8q - 4 s2(q) - 1}
 (4) Krawtchouk (B) OPs == v^m mod 2;  CDH OPs == (v+1)^q mod 2;  T == Pascal mod 2
 (5) v2(q_{K-M}) == H_M + sum_{m<K-M} (b_m + 1)  for all M, q(Y) = det(A(C) + Y B)."""
import sys, math, io, contextlib
from lemlib import *
from sympy import euler
import twisted_hankel as TH
K = int(sys.argv[1]); KM = int(sys.argv[2]) if len(sys.argv) > 2 else 10*K + 60  # relative precision only needs KM > max h_q ~ 8K
muv = [F((2*i + 1)*abs(int(euler(2*i))), 2) for i in range(max(KM, 2*K + 2))]
m = [sum(math.comb(k, i)*(-1)**(k - i)*muv[i] for i in range(k + 1)) for k in range(KM)]
# Taylor expansion of 1/N_K(v) at v = 1:  N_K(1 + e) as a polynomial in e
Npoly = [F(1)]
for a in range(K):
    c0 = 1 + 4*a*a                       # v + 4a^2 = (1 + 4a^2) + e
    new = [F(0)]*(len(Npoly) + 1)
    for i, cc in enumerate(Npoly): new[i] += c0*cc; new[i + 1] += cc
    Npoly = new
pi = [F(1)/Npoly[0]]
for k in range(1, KM):
    s = sum(Npoly[i]*pi[k - i] for i in range(1, min(k, K) + 1))
    pi.append(-s/Npoly[0])
def expand_vs(s):
    # v^s = (1+e)^s -> coefficients in e, times pi(e), truncated
    vs = [F(math.comb(s, i)) for i in range(s + 1)]
    out = [F(0)]*KM
    for i, a in enumerate(vs):
        for k in range(KM - i): out[i + k] += a*pi[k]
    return out
AC = [sum(cf*m[k] for k, cf in enumerate(expand_vs(s))) for s in range(2*K - 1)]
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    A0, B0, h, labels, head, tail = TH.entries('twist', K)
sc = lambda s: F(4)**(s - K)
A0v = [A0[s]*sc(s) for s in range(2*K - 1)]; B0v = [B0[s]*sc(s) for s in range(2*K - 1)]
C200 = Cconst(200)
cons = min(v2(AC[s] - A0v[s] - C200*B0v[s]) for s in range(2*K - 1))
Cst = sum((-1)**k*m[k] for k in range(KM))
print('K=%d KM=%d: (1) min v2(A(C) - A(0) - C B) = %d   (2) v2(C_moments - C_mahler) = %d' % (K, KM, cons, v2(Cst - C200)))
H = hankel(AC, K)
hqs = [hq('twist', q) for q in range(K)]
print('   (3) SNF(A(C)) == h_q: %s' % (snf_v2(H) == sorted(hqs)))
Pk, hB = monic_ops(B0v, K)
Pc, hc = monic_ops(muv[:2*K + 1], K)
mod = lambda x, mm: (x.numerator*pow(x.denominator, -1, mm)) % mm
kr_ok = all(all(mod(Pk[mm][i], 2) == (1 if i == mm else 0) for i in range(mm + 1)) for mm in range(K))
cdh_ok = all(all(mod(Pc[q][i], 2) == math.comb(q, i) % 2 for i in range(q + 1)) for q in range(K))
T = [[F(0)]*K for _ in range(K)]
for mm in range(K):
    rem = Pk[mm][:] + [F(0)]*(K - len(Pk[mm]))
    for q in range(mm, -1, -1):
        coef = rem[q]; T[mm][q] = coef
        for i in range(q + 1): rem[i] -= coef*Pc[q][i]
T_int = all(v2(x) is None or v2(x) >= 0 for r in T for x in r)
pascal = all(mod(T[mm][q], 2) == math.comb(mm, q) % 2 for mm in range(K) for q in range(K))
print('   (4) CDH pivots %s; Krawtchouk == v^m mod 2 %s; CDH == (v+1)^q mod 2 %s; T integral %s; T == Pascal mod 2 %s'
      % (all(v2(hc[q]) == hqs[q] for q in range(K)), kr_ok, cdh_ok, T_int, pascal))
Af = fm(H); Bf = fm(hankel(B0v, K))
dB = Bf.det(); cp = (-Bf.solve(Af)).charpoly()
q = [F(int((cp[i]*dB).p), int((cp[i]*dB).q)) for i in range(K + 1)]
b = band('twist', K)
bad = []
for M in range(K + 1):
    pred = sum(hqs[:M]) + sum(b[mm] + 1 for mm in range(K - M))
    got = v2(q[K - M]) if q[K - M] != 0 else None
    if got != pred: bad.append((M, got, pred))
print('   (5) coefficient law for all M: violations %s' % bad)
