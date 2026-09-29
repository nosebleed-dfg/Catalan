"""Twisted model with the 2-adic mu(-1) = 4C:  A_model = -C*B + P + P_C,  P_C[s] = 4C h_{s-K+1}(z) (z_a = -a^2),
i.e. the constant part of the smooth pole values -4 G1(a) = 4C - 4(G1(a) - G1(0)).
Compare SNF(A_model) with SNF(A) and the union law (v-scale)."""
import sys, io, contextlib, math
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from fractions import Fraction as F
import twisted_hankel as TH
import twoadic_law as TL
from th_lemmaA2 import twist_measure, v2, s2
from lemlib import Cconst, snf_v2
C = Cconst(160)
def hcomp(zs, tmax):
    """complete homogeneous h_t(zs), t = 0..tmax"""
    h = [F(1)] + [F(0)]*tmax
    for z in zs:
        for t in range(1, tmax + 1): h[t] += z*h[t - 1]
    return h
for K in [int(x) for x in sys.argv[1].split(',')]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        A, B, h, labels, head, tail = TH.entries('twist', K)
    nodes, w = twist_measure(K)
    pole = [sum(w[a]*nodes[a]**s for a in range(K)) for s in range(2*K - 1)]
    mom = [A[s] - pole[s] for s in range(2*K - 1)]
    wB = [F(4*(-1)**a)/math.prod(b*b - a*a for b in range(K) if b != a) for a in range(K)]
    Bseq = [sum(wB[a]*nodes[a]**s for a in range(K)) for s in range(2*K - 1)]
    assert all(Bseq[s] == B[s] for s in range(2*K - 1)), 'B is not the Krawtchouk Hankel?'
    hz = hcomp([F(-a*a) for a in range(K)], K)
    PC = [4*C*hz[s - K + 1] if s >= K - 1 else F(0) for s in range(2*K - 1)]
    # check PC equals the Hankel of constant pole values 4C
    wc = [4*C/math.prod(b*b - a*a for b in range(K) if b != a) for a in range(K)]
    assert all(sum(wc[a]*nodes[a]**s for a in range(K)) == PC[s] for s in range(2*K - 1))
    sc = lambda s: F(4)**(s - K)
    H = lambda seq: [[seq[i + j]*sc(i + j) for j in range(K)] for i in range(K)]
    model = [-C*Bseq[s] + mom[s] + PC[s] for s in range(2*K - 1)]
    sA, sM = snf_v2(H(A)), snf_v2(H(model))
    Ms, pred, mn = TL.unified_snf('twist', K)
    print('twist K=%d M*=%s\n  SNF(A)     %s\n  SNF(model) %s  %s\n  union      %s' % (K, Ms, sA, sM, 'SAME' if sA == sM else 'DIFF', pred))
    sys.stdout.flush()
