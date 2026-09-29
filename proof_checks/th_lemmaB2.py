"""Structure of P' (moment corner in the Krawtchouk OP basis) and of its principal minors.
Prints the v2-matrix of P', then for every M: v2(det P'_[K-M,K)) vs H_M = sum_{q<M} h_q,
and the minimum over ALL M-subsets U (small K) of v2(det P'_U) - sum_{m in U} b_m, with its argmin."""
import sys, itertools
from lemlib import *
C = Cconst()
mach = sys.argv[1]; K = int(sys.argv[2]); allsub = len(sys.argv) > 3
parts = machine_parts(mach, K, C)
P, hB = monic_ops(parts['B'], K)
Pp = gram(P, parts['Pext'])
bd = band(mach, K)
print('%s K=%d  v2 of P\'_mn (. = 0):' % (mach, K))
for m in range(K):
    print('  %2d | ' % m + ' '.join(('%4d' % v2(x)) if x != 0 else '   .' for x in Pp[m]))
H = [0]
for q in range(K): H.append(H[-1] + hq(mach, q))
print('top blocks: M, v2(det P\'_[K-M,K)), H_M')
for M in range(1, K + 1):
    sub = [r[K - M:] for r in Pp[K - M:]]
    d = fm(sub).det()
    dv = v2(F(int(d.p), int(d.q))) if d != 0 else None
    print('  M=%2d  %s  %d  %s' % (M, dv, H[M], 'EQ' if dv == H[M] else ''))
if allsub:
    for M in range(1, min(K, 6) + 1):
        best = None
        for U in itertools.combinations(range(K), M):
            sub = [[Pp[i][j] for j in U] for i in U]
            d = fm(sub).det()
            if d == 0: continue
            val = v2(F(int(d.p), int(d.q))) - sum(bd[m] for m in U)
            if best is None or val < best[0]: best = (val, [U])
            elif val == best[0]: best[1].append(U)
        print('M=%d: min over U of v2(det P\'_U) - sum b_U = %s at %s   (top: %d)' % (M, best[0], best[1][:4], H[M] - sum(bd[K - M:])))
