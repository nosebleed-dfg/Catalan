"""Shared helpers for the Lemma A/B/C work (no side effects on import)."""
import sys, io, contextlib, math
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from fractions import Fraction as F
import flint
s2 = lambda n: bin(n).count('1')
def v2(x):
    x = F(x); a, b = x.numerator, x.denominator
    if a == 0: return None
    return ((abs(a) & -abs(a)).bit_length() - 1) - ((b & -b).bit_length() - 1)
def v2i(n):
    n = abs(n); return (n & -n).bit_length() - 1 if n else None
_C = {}
def Cconst(nmax=160):
    """2-adic constant with beta_a = C + (-1)^a G1(a), G1 smooth; exact rational truncation (error O(2^~nmax))."""
    if nmax in _C: return _C[nmax]
    C = F(0)
    for n in range(nmax + 1):
        d = sum((-1)**(n - i)*math.comb(n, i)*F(1, (2*i + 1)**2) for i in range(n + 1))
        C += F((-1)**n, 2**(n + 1))*d
    _C[nmax] = C
    return C
def betas(n):
    beta = [F(0)]
    for k in range(n): beta.append(beta[-1] + F((-1)**k, (2*k + 1)**2))
    return beta
def twist_measure(K):
    beta = betas(K)
    nodes = [F(-a*a) for a in range(K)]
    w = [F(-4*(-1)**a)*beta[a]/math.prod(b*b - a*a for b in range(K) if b != a) for a in range(K)]
    return nodes, w
def cat_measure(K):
    beta = betas(K)
    nodes = [F(-(2*j + 1)**2) for j in range(K)]
    w = []
    for j in range(K):
        b = 2*j + 1
        Fb = F(-(-1)**j*b, 2)*beta[j] - F(1, 4*b)
        w.append(Fb/math.prod((2*i + 1)**2 - b*b for i in range(K) if i != j))
    return nodes, w
def snf_v2(rows):
    K = len(rows); L = 1
    for r in rows:
        for x in r: L = L*x.denominator//math.gcd(L, x.denominator)
    s = v2i(L)
    S = flint.fmpz_mat([[int(x*L) for x in r] for r in rows]).snf()
    return sorted((v2i(int(S[i, i])) - s) if int(S[i, i]) else 10**9 for i in range(K))
def hcomp(zs, tmax):
    h = [F(1)] + [F(0)]*tmax
    for z in zs:
        for t in range(1, tmax + 1): h[t] += z*h[t - 1]
    return h
def fm(rows): return flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in r] for r in rows])
def tofr(M): return [[F(int(M[i, j].p), int(M[i, j].q)) for j in range(M.ncols())] for i in range(M.nrows())]
def machine_parts(mach, K, C=None):
    """v-scale Hankel sequences (length 2K-1): A (X=0), B (X-coefficient = Krawtchouk), Pext (moment corner, twisted incl. mu(-1)=4C),
    pole (true X=0 pole part), mom (true moment corner)."""
    import twisted_hankel as TH
    if C is None: C = Cconst()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        A, B, h, labels, head, tail = TH.entries(mach, K)
    nodes, w = (twist_measure if mach == 'twist' else cat_measure)(K)
    pole = [sum(w[a]*nodes[a]**s for a in range(K)) for s in range(2*K - 1)]
    mom = [A[s] - pole[s] for s in range(2*K - 1)]
    if mach == 'twist':
        hz = hcomp(nodes, K)
        ext = [mom[s] + (4*C*hz[s - K + 1] if s >= K - 1 else 0) for s in range(2*K - 1)]
        sc = lambda s: F(4)**(s - K)
    else:
        ext = mom[:]
        sc = lambda s: F(1)
    V = lambda seq: [seq[s]*sc(s) for s in range(2*K - 1)]
    return dict(A=V(A), B=V(B[:2*K - 1]), Pext=V(ext), pole=V(pole), mom=V(mom))
def hankel(seq, K): return [[seq[i + j] for j in range(K)] for i in range(K)]
def monic_ops(mom, K):
    """monic OPs (coefficient lists low->high) and norms for the functional v^s -> mom[s] (Chebyshev/Stieltjes, exact)."""
    def L(poly, pw=0): return sum(c*mom[i + pw] for i, c in enumerate(poly))
    def mul(p, q):
        r = [F(0)]*(len(p) + len(q) - 1)
        for i, a in enumerate(p):
            for j, b in enumerate(q): r[i + j] += a*b
        return r
    P = [[F(1)]]; h = []; prev = [F(0)]
    for n in range(K):
        pp = mul(P[-1], P[-1]); hn = L(pp); h.append(hn)
        if n == K - 1: break
        bn = L(pp, 1)/hn; lam = hn/h[-2] if n else F(0)
        new = [F(0)] + P[-1]
        for i in range(len(P[-1])): new[i] -= bn*P[-1][i]
        for i in range(len(prev)): new[i] -= lam*prev[i]
        prev = P[-1]; P.append(new)
    return P, h
def gram(P, seq):
    K = len(P)
    return [[sum(a*b*seq[i + j] for i, a in enumerate(P[m]) for j, b in enumerate(P[n])) for n in range(K)] for m in range(K)]
def band(mach, K):
    c = 0 if mach == 'twist' else -1
    return [1 + s2(K - 1 - m) - s2(m) + 4*m - 2*K + c for m in range(K)]
def hq(mach, q):
    return (8*q - 4*s2(q) - 1) if mach == 'twist' else (4*q + 3*(q - s2(q)) + (q + 1 - s2(q + 1)) - 1)
