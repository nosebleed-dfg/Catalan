"""Lemma B in the Krawtchouk basis.  Model G = -C*B + P_ext written in the monic OPs p_n of B (v-scale, 2-integral):
G = D + P',  D = diag(-C h^B_n)  (v2 = band b_n),  P'_mn = <p_m, p_n>_{P_ext}  (0 for m+n < K (Catalan), < K-1 (twisted)).
At the tropical split s = K - M*:  T = [0, s), Bo = [s, K).
 (1) SNF(G_TT) == band[0:s]?   (2) min v2 of G_BT G_TT^{-1} (>= 0: block elimination unimodular)
 (3) SNF(S), S = G_BB - G_BT G_TT^{-1} G_TB,  vs moment values h_q, q < M*    (4) SNF(P'_BB)   (5) SNF(G) vs union."""
import sys
from lemlib import *
import twoadic_law as TL
C = Cconst()
for mach in sys.argv[1].split(','):
    for K in [int(x) for x in sys.argv[2].split(',')]:
        parts = machine_parts(mach, K, C)
        P, hB = monic_ops(parts['B'], K)
        integral = all(v2(c) is None or v2(c) >= 0 for p in P for c in p)
        bd = band(mach, K)
        Dv = [-C*x for x in hB]
        okD = [v2(x) for x in Dv] == bd
        Pp = gram(P, parts['Pext'])
        G = [[(Dv[m] if m == n else F(0)) + Pp[m][n] for n in range(K)] for m in range(K)]
        Ms, pred, mn = TL.unified_snf(mach, K)
        Mstar = Ms[0]; s = K - Mstar
        GTT = [r[:s] for r in G[:s]]; GTB = [r[s:] for r in G[:s]]; GBT = [r[:s] for r in G[s:]]; GBB = [r[s:] for r in G[s:]]
        iTT = fm(GTT).inv()
        X = tofr(fm(GBT)*iTT)
        minX = min((v2(x) for r in X for x in r if x != 0), default=None)
        S = tofr(fm(GBB) - fm(GBT)*iTT*fm(GTB))
        PBB = [r[s:] for r in Pp[s:]]
        sG = snf_v2(G)
        print('%s K=%d M*=%s s=%d  OPs integral %s  v2(D)=band %s' % (mach, K, Ms, s, integral, okD))
        print('   (1) SNF(G_TT) %s  band[:s] %s  %s' % (snf_v2(GTT), sorted(bd[:s]), 'OK' if snf_v2(GTT) == sorted(bd[:s]) else 'DIFF'))
        print('   (2) min v2(G_BT G_TT^-1) = %s' % minX)
        print('   (3) SNF(S) %s   h_q %s' % (snf_v2(S), [hq(mach, q) for q in range(Mstar)]))
        print('   (4) SNF(P\'_BB) %s   D_B %s' % (snf_v2(PBB), bd[s:]))
        print('   (5) SNF(G) %s  union %s  %s' % (sG, pred, 'OK' if sG == pred else 'DIFF'))
        sys.stdout.flush()
