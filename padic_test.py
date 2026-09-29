"""padic_test.py — mechanism test for the p-adic exponent laws (HANDOFF open item 2).

Builds the half-integer-pole Hankel machine (poles v=-b^2, b odd, N=0) with a chosen pair
(pole-value type, moment type) — positivity is irrelevant for the content, so hybrids are allowed —
and prints the content exponents e_p for p > K/2 against the two laws
    plain law      4n+2 = 2(2K-p)          (von Staudt 1/p in the moment B_{p-1} survives)
    alternating law 3n+1 = (3(2K-p)-1)/2   (Genocchi moments: no von Staudt term)
n = #{poles b > p}.  Prediction: the law is decided by the MOMENTS alone for p > K.
usage: python padic_test.py K PV MOM [POLES]     PV, MOM in {zeta2, catalan}; POLES in {half (default), int}
  POLES=int: integer poles u=-j^2 (j=1..K) with the zeta2int pole values (j/2)(zeta(2)-H2_{j-1}) - 1/2 - 1/(4j)
             (PV is ignored); moments |B_{2e+2}|/2 (zeta2) or (2^{2e+2}-1)|B_{2e+2}| (catalan/Genocchi).
"""
import sys, math; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, primerange
import flint
K=int(sys.argv[1]); PV=sys.argv[2]; MOM=sys.argv[3]; POLES=sys.argv[4] if len(sys.argv)>4 else 'half'; h=K
def bern(n): B=bernoulli(n); return F(int(B.p),int(B.q))
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b): r[i+j]+=ai*bj
    return r
def ev(p,v):
    s=0
    for c in reversed(p): s=s*v+c
    return s
if POLES=='int':
    labels=list(range(1,K+1))
    if MOM=='zeta2': mu=lambda e: abs(bern(2*e+2))/2
    else:            mu=lambda e: (2**(2*e+2)-1)*abs(bern(2*e+2))
    H2={0:F(0)}
    for j in range(1,K+1): H2[j]=H2[j-1]+F(1,j*j)
    pv=lambda j: (-F(j,2)*H2[j-1]-F(1,2)-F(1,4*j), F(j,2))
else:
    labels=[2*j+1 for j in range(K)]
    if MOM=='zeta2':   mu=lambda e: F(4)**e*abs(bern(2*e+2))/2
    else:              mu=lambda e: F(4)**e*(2**(2*e+2)-1)*abs(bern(2*e+2))
    if PV=='zeta2':
        beta={}; acc=F(0)
        for b in range(1,2*K+2,2): beta[b]=acc; acc+=F(1,b*b)
        pv=lambda b: (-F(b,4)*beta[b]-F(1,8*b)-F(1,8), F(b,16))
    else:
        betaalt={}; acc=F(0); sgn=1
        for b in range(1,2*K+2,2): betaalt[b]=acc; acc+=F(sgn,b*b); sgn=-sgn
        def pv(b):
            sign=(-1)**((b-1)//2); return (F(-sign*b,2)*betaalt[b]-F(1,4*b), F(sign*b,2))
D=[1]
for x in labels: D=polymul(D,[x*x,1])
degD=len(D)-1
Dp={a:math.prod(b*b-a*a for b in labels if b!=a) for a in labels}
PVs={a:pv(a) for a in labels}
q,r=[],[0]*degD; r[:1]=[1]          # v^0 / D : quotient 0, remainder 1
r=[1]+[0]*(degD-1)
MU=[mu(e) for e in range(2*h+2)]
pw={a:1 for a in labels}; A=[]; B=[]
for s in range(2*h-1):
    a_val=F(0); b_val=F(0)
    for e,c in enumerate(q):
        if c: a_val+=MU[e]*c
    for a in labels:
        c=F(pw[a],Dp[a]); a_val+=c*PVs[a][0]; b_val+=c*PVs[a][1]
    A.append(a_val); B.append(b_val)
    cs=r[degD-1]; q=[cs]+q
    r=[(r[i-1] if i>0 else 0)-cs*D[i] for i in range(degD)]
    for a in labels: pw[a]*=-a*a
Af=flint.fmpq_mat([[flint.fmpq(A[i+j].numerator,A[i+j].denominator) for j in range(h)] for i in range(h)])
Bf=flint.fmpq_mat([[flint.fmpq(B[i+j].numerator,B[i+j].denominator) for j in range(h)] for i in range(h)])
vals=[]
for X in range(h+1):
    d=(Af+Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p),int(d.q)))
xs=list(range(h+1)); coef=vals[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1): coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
poly=[coef[h]]
for k in range(h-1,-1,-1):
    poly=polymul(poly,[F(-xs[k]),F(1)]); poly=[F(c) for c in poly]; poly[0]+=coef[k]
g=0; l=1
for c in poly: g=gcd(g,c.numerator); l=l*c.denominator//gcd(l,c.denominator)
den=l//gcd(l,g)
out=[]
for p in primerange(3,2*K):
    e=0
    while den%p==0: den//=p; e+=1
    if p>K/2:
        if POLES=='int':
            vV=sum(1 for a in labels for b in labels if a<b and ((a+b)%p==0 or (b-a)%p==0))   # v_p(det V), single p per pair here
            out.append((p,e,e-2*vV,0))
        else:
            n=sum(1 for b in labels if b>p)
            out.append((p,e,e-(4*n+2),e-(3*n+1)))
if POLES=='int':
    print('POLES=int MOM=%s K=%d   p:e_p [resid vs 2 v_p(det V)]'%(MOM,K))
    print(' '.join('%d:%d[%+d]'%(t[0],t[1],t[2]) for t in out))
    sel=[t for t in out if t[0]>K]
    print('  primes K<p<2K: mean e_p - 2v_p(det V) = %.2f'%(sum(t[2] for t in sel)/max(1,len(sel))))
else:
    print('PV=%s MOM=%s K=%d   p:e_p [resid vs 4n+2 | vs 3n+1]'%(PV,MOM,K))
    print(' '.join('%d:%d[%+d|%+d]'%t for t in out))
    print('  primes K<p<2K:  mean resid vs 4n+2 = %.2f   vs 3n+1 = %.2f'%(sum(t[2] for t in out if t[0]>K)/max(1,sum(1 for t in out if t[0]>K)), sum(t[3] for t in out if t[0]>K)/max(1,sum(1 for t in out if t[0]>K))))
