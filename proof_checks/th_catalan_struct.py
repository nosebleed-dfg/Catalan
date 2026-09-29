"""Catalan structure theorem, end-to-end check (v-scale).
 A(C)_s = phi0( v^s / N_K(v) ) computed 2-adically as sum_k c_k mu_{s+k}  (c_k = Taylor coefficients of 1/N_K at 0) -- no C needed.
 (1) consistency: A(C) == A(0) + C*B  (C from the Mahler series, precision ~2^190)
 (2) C = 1/2 + 2 sum_k (-1)^k mu_k  agrees with the Mahler-series C
 (3) SNF(A(C)) == {h_q : q < K}; leading pivots v2 == h_q
 (4) Krawtchouk polys p_m == (v+1)^m mod 2; CDH polys P_q == v^q mod 4; T (p in P) == Pascal mod 2
 (5) q(Y) = det(A(C) + Y B):  v2(q_{K-M}) == H_M + sum_{m<K-M} (b_m + 1)   for EVERY M = 0..K."""
import sys, math, io, contextlib
from lemlib import *
from sympy import bernoulli
import twisted_hankel as TH
K = int(sys.argv[1]); KM = int(sys.argv[2]) if len(sys.argv) > 2 else 6*K + 60   # relative precision only needs KM > max h_q
mu = []
for k in range(KM + 2*K):
    B_ = bernoulli(2*k + 2); mu.append(F(4)**k*(2**(2*k + 2) - 1)*abs(F(int(B_.p), int(B_.q))))
# Taylor coefficients of 1/N_K(v), N_K(v) = prod (v + (2j+1)^2)
Npoly = [F(1)]
for j in range(K):
    a = (2*j + 1)**2; new = [F(0)]*(len(Npoly) + 1)
    for i, c in enumerate(Npoly): new[i] += a*c; new[i + 1] += c
    Npoly = new
c = [F(1)/Npoly[0]]
for k in range(1, KM):
    s = sum(Npoly[i]*c[k - i] for i in range(1, min(k, K) + 1))
    c.append(-s/Npoly[0])
AC = [sum(c[k]*mu[s + k] for k in range(KM)) for s in range(2*K - 1)]
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    A0, B0, h, labels, head, tail = TH.entries('catalan', K)
C200 = Cconst(200)
cons = min(v2(AC[s] - A0[s] - C200*B0[s]) for s in range(2*K - 1))
Cst = F(1, 2) + 2*sum((-1)**k*mu[k] for k in range(KM))
print('K=%d KM=%d: (1) min v2(A(C) - A(0) - C B) = %d   (2) v2(C_stieltjes - C_mahler) = %d' % (K, KM, cons, v2(Cst - C200)))
H = hankel(AC, K)
hqs = [hq('catalan', q) for q in range(K)]
print('   (3) SNF(A(C)) == h_q: %s' % (snf_v2(H) == sorted(hqs)))
Pk, hB = monic_ops(B0[:2*K - 1], K)          # Krawtchouk (B functional)
# CDH monic OPs from the true moments mu (v-scale)
Pc, hc = monic_ops(mu[:2*K + 1], K)
pc_ok = all(v2(hc[q]) == hqs[q] for q in range(K))
mod = lambda x, m: (x.numerator*pow(x.denominator, -1, m)) % m
def binom_coeffs(m): return [math.comb(m, i) for i in range(m + 1)]
kr_ok = all(all(mod(Pk[m][i], 2) == math.comb(m, i) % 2 for i in range(m + 1)) for m in range(K))
cdh_ok = all(all(mod(Pc[q][i], 4) == (1 if i == q else 0) for i in range(q + 1)) for q in range(K))
# T: p_m = sum_q T_mq P_q  (solve triangular)
T = [[F(0)]*K for _ in range(K)]
for m in range(K):
    rem = Pk[m][:] + [F(0)]*(K - len(Pk[m]))
    for q in range(m, -1, -1):
        coef = rem[q]; T[m][q] = coef
        for i in range(q + 1): rem[i] -= coef*Pc[q][i]
T_int = all(v2(x) is None or v2(x) >= 0 for r in T for x in r)
pascal = all(mod(T[m][q], 2) == math.comb(m, q) % 2 for m in range(K) for q in range(K))
print('   (4) CDH pivots h_q %s; Krawtchouk == (v+1)^m mod 2 %s; CDH == v^q mod 4 %s; T integral %s, T == Pascal mod 2 %s' % (pc_ok, kr_ok, cdh_ok, T_int, pascal))
# (5)
Af = fm(H); Bf = fm(hankel(B0[:2*K - 1], K))
dB = Bf.det(); cp = (-Bf.solve(Af)).charpoly()
q = [F(int((cp[i]*dB).p), int((cp[i]*dB).q)) for i in range(K + 1)]   # det(A(C) + Y B) = sum q_i Y^i
b = band('catalan', K)
bad = []
for M in range(K + 1):
    pred = sum(hqs[:M]) + sum(b[m] + 1 for m in range(K - M))
    got = v2(q[K - M]) if q[K - M] != 0 else None
    if got != pred: bad.append((M, got, pred))
print('   (5) v2(q_{K-M}) == H_M + sum(b+1) for all M = 0..K: violations %s' % bad)
