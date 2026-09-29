"""Krawtchouk base model: replace the pole weights by their 'C-part' and keep the true moment corner.
twisted (u-scale, nodes -a^2):   omega_a = -4 eps_a beta_a / ((K-1-a)!(K-1+a)!),   beta_a -> C
Catalan (v-scale, nodes -(2j+1)^2): omega_j = (2j+1)^2 (-1)^j DeltaG(j) / (4^K (K-1-j)!(K+j)!),  (-1)^j DeltaG(j) -> -2C
C = 2-adic constant with beta_a = C + (-1)^a G1(a), G1 2-adically smooth:  C = sum_n (-1)^n Delta^n g(0) / 2^(n+1), g(a) = (2a+1)^-2.
Checks: formula for omega vs direct; smoothness of G1; SNF of (Kraw + corner) vs SNF(A) vs union law; SNF(Kraw) vs band."""
import sys, io, contextlib, math
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from fractions import Fraction as F
import flint
import twisted_hankel as TH
import twoadic_law as TL
from th_lemmaA2 import twist_measure, cat_measure, v2, s2
from lemlib import snf_v2, band
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 160
def Cconst(nmax):
    g = [F(1, (2*i + 1)**2) for i in range(nmax + 1)]
    C = F(0)
    for n in range(nmax + 1):
        d = sum((-1)**(n - i)*math.comb(n, i)*g[i] for i in range(n + 1))
        C += F((-1)**n, 2**(n + 1))*d
    return C
C = Cconst(NC)
# precision estimate: v2 of the last terms
print('C ~ 2-adic, v2(C) = %d; truncated at n = %d' % (v2(C), NC))
beta = [F(0)]
for k in range(200): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
G1 = [(-1)**a*(beta[a] - C) for a in range(60)]
def fd(h, n, x): return sum((-1)**(n - i)*math.comb(n, i)*h[x + i] for i in range(n + 1))
worst = min(v2(fd(G1, n, x)) - (n - 1 + (math.factorial(n + 1).bit_length() and (n + 1 - s2(n + 1)))) for n in range(0, 30) for x in range(0, 25))
print('G1 smoothness: min over n<30,x<25 of v2(Delta^n G1(x)) - (n-1+v2((n+1)!)) = %d  (>= 0 expected)' % worst)
for K in [int(x) for x in sys.argv[1].split(',')]:
    for mach in ['twist', 'catalan']:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            A, B, h, labels, head, tail = TH.entries(mach, K)
        nodes, w = (twist_measure if mach == 'twist' else cat_measure)(K)
        if mach == 'twist':
            eps = lambda a: 1 if a == 0 else 2
            wf = [F(-4*eps(a))*beta[a]/(math.factorial(K - 1 - a)*math.factorial(K - 1 + a)) for a in range(K)]
            wk = [F(-4*eps(a))*C/(math.factorial(K - 1 - a)*math.factorial(K - 1 + a)) for a in range(K)]
        else:
            DG = lambda j: (-1)**(j + 1)*beta[j + 1] - (-1)**j*beta[j]
            wf = [F((2*j + 1)**2)*(-1)**j*DG(j)/(4**K*math.factorial(K - 1 - j)*math.factorial(K + j)) for j in range(K)]
            wk = [F((2*j + 1)**2)*(-2*C)/(4**K*math.factorial(K - 1 - j)*math.factorial(K + j)) for j in range(K)]
        assert wf == w, 'closed-form weights disagree with direct F/D\''
        pole = [sum(w[a]*nodes[a]**s for a in range(K)) for s in range(2*K - 1)]
        kraw = [sum(wk[a]*nodes[a]**s for a in range(K)) for s in range(2*K - 1)]
        mom = [A[s] - pole[s] for s in range(2*K - 1)]
        sc = (lambda s: F(4)**(s - K)) if mach == 'twist' else (lambda s: F(1))
        H = lambda seq: [[seq[i + j]*sc(i + j) for j in range(K)] for i in range(K)]
        sA = snf_v2(H(A)); sK = snf_v2(H(kraw)); sKP = snf_v2(H([kraw[s] + mom[s] for s in range(2*K - 1)]))
        Ms, pred, mn = TL.unified_snf(mach, K)
        print('%s K=%d M*=%s\n  SNF(A)        %s\n  SNF(Kraw+mom) %s  %s\n  union law     %s\n  SNF(Kraw)==band: %s' % (
            mach, K, Ms, sA, sKP, 'SAME' if sKP == sA else 'DIFF', pred, sK == sorted(band(mach, K))))
        sys.stdout.flush()
