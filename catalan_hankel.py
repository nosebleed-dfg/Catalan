import sys, time, math
sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from sympy import bernoulli
import flint
from mpmath import mp, mpf, catalan, log as mlog
t0=time.time()
n=int(sys.argv[1]); m=int(sys.argv[2]); K=40*n; N=int(sys.argv[3])*n; h=K-N
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b): r[i+j]+=ai*bj
    return r
odd=lambda j: 2*j+1
EN=[1]
for j in range(N): EN=polymul(EN,[odd(j)**2,1])
ENm=[1]
for _ in range(m): ENm=polymul(ENm,EN)
tail=list(range(N,K))                      # poles at -(2j+1)^2, j=N..K-1
Et=[1]
for j in tail: Et=polymul(Et,[odd(j)**2,1])
def polydivmod(num,den):
    num=[F(c) for c in num]; den=[F(c) for c in den]
    q=[F(0)]*(max(0,len(num)-len(den)+1)); r=num[:]
    for k in range(len(num)-len(den),-1,-1):
        c=r[k+len(den)-1]/den[-1]; q[k]=c
        for i,d in enumerate(den): r[k+i]-=c*d
    return q, r
def evalpoly(p,v):
    s=0
    for c in reversed(p): s=s*v+c
    return s
Etp={}
for a in tail:
    prod=1
    for b in tail:
        if b!=a: prod*=(odd(b)**2-odd(a)**2)
    Etp[a]=prod
beta={0:F(0)}
for j in range(1,K+1): beta[j]=beta[j-1]+F((-1)**(j-1),odd(j-1)**2)   # beta_j = sum_{k<j} (-1)^k/(2k+1)^2
def phi_mono(e):
    B=bernoulli(2*e+2); return F(4)**e*(2**(2*e+2)-1)*abs(F(int(B.p),int(B.q)))
tp={}
for s in range(2*h-1):
    q,_=polydivmod([0]*s+ENm,Et)
    a_val=sum(c*phi_mono(e) for e,c in enumerate(q) if c); b_val=F(0)
    for j in tail:
        b=odd(j); c=F(evalpoly(ENm,-b*b)*((-b*b)**s), Etp[j])
        sign=(-1)**j
        a_val+=c*(F(-sign*b,2)*beta[j]-F(1,4*b))     # phi(1/(u+b^2)) = (b/2)(-1)^j (X - beta_j) - 1/(4b)
        b_val+=c*F(sign*b,2)
    tp[s]=(a_val,b_val)
Af=flint.fmpq_mat([[flint.fmpq(tp[i+j][0].numerator,tp[i+j][0].denominator) for j in range(h)] for i in range(h)])
Bf=flint.fmpq_mat([[flint.fmpq(tp[i+j][1].numerator,tp[i+j][1].denominator) for j in range(h)] for i in range(h)])
vals=[]
for X in range(h+1):
    d=(Af+Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p),int(d.q)))
xs=list(range(h+1)); coef=vals[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1): coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
poly=[coef[h]]
for k in range(h-1,-1,-1):
    poly=polymul(poly,[F(-xs[k]),F(1)]); poly=[F(c) for c in poly]; poly[0]+=coef[k]
from math import gcd
g=0; l=1
for c in poly:
    g=gcd(g,c.numerator); l=l*c.denominator//gcd(l,c.denominator)
P=[c/F(g,l) for c in poly]
mp.dps=max(3000,8*max(len(str(abs(c.numerator))) for c in P))
Gv=+catalan
val=mpf(0)
for c in reversed(P): val=val*Gv+mpf(c.numerator)/mpf(c.denominator)
vD=mpf(0)
for c in reversed(poly): vD=vD*Gv+mpf(c.numerator)/mpf(c.denominator)
print('  log Delta_K(G)/K^2 = %.4f   -log content/K^2 = %.4f'%(float(mlog(abs(vD)))/K**2, -float(mlog(abs(F(g,l))))/K**2))
print('Catalan Hankel test: n=%d m=%d K=%d h=%d  deg=%d  P_K(G)>0: %s  log P_K(G) = %.1f  per n^2 = %.1f  max coeff digits %d  time %.0fs'%(n,m,K,h,len(P)-1,val>0,float(mlog(abs(val))),float(mlog(abs(val)))/n**2,max(len(str(abs(c.numerator))) for c in P),time.time()-t0))
