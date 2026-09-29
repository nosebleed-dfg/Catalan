"""Check Theorem 4's explicit integer recurrences for the v-scale Krawtchouk polynomials of B:
 twisted: b = -2(4n(K-1-n) + K - 1), lam = 4n(2n-1)(K-n)(2K-1-2n);  Catalan: b = -[(2n+1)(2K-1-2n) + 4(n+1)(K-1-n)], lam = 4n(K-n)(2n+1)(2K-1-2n);
 and the norm valuations v2(beta_n) = kappa(K,n) + 4n - 2K + 1 + c."""
import sys
from lemlib import *
for K in [int(x) for x in sys.argv[1].split(',')]:
    for mach in ['twist', 'catalan']:
        parts = machine_parts(mach, K)
        P, hB = monic_ops(parts['B'], K)
        if mach == 'twist':
            bf = lambda n: -2*(4*n*(K - 1 - n) + K - 1); lf = lambda n: 4*n*(2*n - 1)*(K - n)*(2*K - 1 - 2*n); c = 0
        else:
            bf = lambda n: -((2*n + 1)*(2*K - 1 - 2*n) + 4*(n + 1)*(K - 1 - n)); lf = lambda n: 4*n*(K - n)*(2*n + 1)*(2*K - 1 - 2*n); c = -1
        bad = 0
        for n in range(K - 1):
            lhs = [F(0)] + P[n]                                     # v p_n
            rhs = P[n + 1][:] + [F(0)]*(len(lhs) - len(P[n + 1]))
            for i, a in enumerate(P[n]): rhs[i] += bf(n)*a
            if n >= 1:
                for i, a in enumerate(P[n - 1]): rhs[i] += lf(n)*a
            if lhs != rhs: bad += 1
        nb = sum(1 for n in range(K) if v2(hB[n]) != 1 + s2(K - 1 - n) - s2(n) + 4*n - 2*K + 1 + c)
        print('%s K=%d: recurrence mismatches %d of %d; norm-valuation mismatches %d of %d' % (mach, K, bad, K - 1, nb, K))
