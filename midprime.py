import sys, math, time; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from sympy import bernoulli, factorint
import flint
K=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); h=K-N
t0=time.time()
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
tail=list(range(N,K)); Et=[1]
for j in tail: Et=polymul(Et,[odd(j)**2,1])
def polydiv(num,den):
    num=[F(c) for c in num]; den=[F(c) for c in den]
    q=[F(0)]*(max(0,len(num)-len(den)+1)); r=num[:]
    for k in range(len(num)-len(den),-1,-1):
        c=r[k+len(den)-1]/den[-1]; q[k]=c
        for i,d in enumerate(den): r[k+i]-=c*d
    return q
def ev(p,v):
    s=0
    for c in reversed(p): s=s*v+c
    return s
Etp={a:math.prod(odd(b)**2-odd(a)**2 for b in tail if b!=a) for a in tail}
beta={0:F(0)}
for j in range(1,K+1): beta[j]=beta[j-1]+F((-1)**(j-1),odd(j-1)**2)
def phi_mono(e):
    B=bernoulli(2*e+2); return F(4)**e*(2**(2*e+2)-1)*abs(F(int(B.p),int(B.q)))
ENm_at={j:ev(ENm,-odd(j)**2) for j in tail}
tp={}
for s in range(2*h-1):
    q=polydiv([0]*s+ENm,Et)
    a_val=sum(c*phi_mono(e) for e,c in enumerate(q) if c); b_val=F(0)
    for j in tail:
        b=odd(j); c=F(ENm_at[j]*((-b*b)**s), Etp[j]); sign=(-1)**j
        a_val+=c*(F(-sign*b,2)*beta[j]-F(1,4*b)); b_val+=c*F(sign*b,2)
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
den=l//gcd(l,g); num=g//gcd(l,g)
fd=factorint(den)
mid={p:e for p,e in fd.items() if K/2<p<2*K}
small={p:e for p,e in fd.items() if p<=K/2}
logden=sum(e*math.log(p) for p,e in fd.items())
P=[c/F(num,den) for c in poly]
H=max(abs(c.numerator) for c in P)
print('K=%d N=%d m=%d h=%d | log(den of content)/K^2=%.3f  mid-prime part/K^2=%.3f  small-prime part/K^2=%.3f | log H(P)/K^2=%.3f | time %.0fs'%(K,N,m,h,logden/K**2,sum(e*math.log(p) for p,e in mid.items())/K**2,sum(e*math.log(p) for p,e in small.items())/K**2,math.log(H)/K**2,time.time()-t0))
print('   mid primes (p:e):',dict(sorted(mid.items())))
print('   small primes (p:e):',dict(sorted(small.items())))
