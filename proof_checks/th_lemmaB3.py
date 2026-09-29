"""Top block of P' in reversed order i = K-1-m: Z_ij = Psi(p_{K-1-i} p_{K-1-j}).
Compare its leading pivots (order i = 0,1,2,..) with the CDH moment pivots h_q, and find the range of i where they agree.
Also test the Pell-type identity v p_{K-1}^2 + c^2 = N_K Q (Catalan)."""
import sys
from lemlib import *
C = Cconst()
for mach in sys.argv[1].split(','):
    for K in [int(x) for x in sys.argv[2].split(',')]:
        parts = machine_parts(mach, K, C)
        P, hB = monic_ops(parts['B'], K)
        Pext = parts['Pext']
        # reversed top block
        idx = list(range(K - 1, -1, -1))
        Z = [[sum(a*b*Pext[i + j] for i, a in enumerate(P[m]) for j, b in enumerate(P[n])) for n in idx] for m in idx]
        piv = []; prev = F(1)
        for k in range(1, K + 1):
            d = fm([r[:k] for r in Z[:k]]).det(); d = F(int(d.p), int(d.q))
            if d == 0 or prev == 0: piv.append(None); prev = d; continue
            piv.append(v2(d/prev)); prev = d
        h = [hq(mach, q) for q in range(K)]
        agree = 0
        while agree < K and piv[agree] == h[agree]: agree += 1
        print('%s K=%d  pivots(Z) %s\n            h_q      %s   agree for i < %d  (K/3 = %.1f)' % (mach, K, piv, h, agree, K/3))
        sys.stdout.flush()
