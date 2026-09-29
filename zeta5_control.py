import sys, time, math; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from sympy import bernoulli, factorint
import flint
from mpmath import mp, mpf, zeta, log as mlog
t0=time.time(); n=int(sys.argv[1]); K=40*n; N=3*n; h=37*n; m=6
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b): r[i+j]+=ai*bj
    return r
def mu_mono(e):
    B=bernoulli(2*e+2); B=F(int(B.p),int(B.q))
    return (-1)**e*B*F((2*e+3)*(2*e+4)*(2*e+5),24)
DN=[1]
for j in range(1,N+1): DN=polymul(DN,[j*j,1])
DNm=[1]
for _ in range(m-1): DNm=polymul(DNm,DN)      # D_N^(m-1): one copy cancels against D_K
tail=list(range(N+1,K+1)); Dt=[1]
for j in tail: Dt=polymul(Dt,[j*j,1])
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
Dtp={a:math.prod(b*b-a*a for b in tail if b!=a) for a in tail}
H5={0:F(0)}
for j in range(1,K+1): H5[j]=H5[j-1]+F(1,j**5)
DNm_at={j:ev(DNm,-j*j) for j in tail}
tp={}
for s in range(2*h-1):
    q=polydiv([0]*s+DNm,Dt)
    a_val=sum(c*mu_mono(e) for e,c in enumerate(q) if c); b_val=F(0)
    for j in tail:
        c=F(DNm_at[j]*((-j*j)**s), Dtp[j])
        a_val+=c*(-(j**4)*H5[j]+F(1,2*j)-F(1,4)); b_val+=c*(j**4)
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
lead_pred=F((-1)**(h*(h-1)//2))
for j in tail: lead_pred*= j**4 * F(ev(DN,-j*j))**5
print('leading coefficient matches Lemma 2.3 analogue:', poly[-1]==lead_pred)
from math import gcd
g=0; l=1
for c in poly:
    g=gcd(g,c.numerator); l=l*c.denominator//gcd(l,c.denominator)
den=l//gcd(l,g); num=g//gcd(l,g)
P=[c/F(num,den) for c in poly]
mp.dps=max(3000,4*max(len(str(abs(c.numerator))) for c in P))
z5=zeta(5); val=mpf(0)
for c in reversed(P): val=val*z5+mpf(c.numerator)/mpf(c.denominator)
H=max(abs(c.numerator) for c in P)
fd=factorint(den); fn=factorint(num)
mid={p:e for p,e in fd.items() if K/2<p<2*K}
print('ZETA(5) CONTROL K=%d: log P_K(zeta5) = %.1f (mo271: -265 at K=40, -833 at K=80)  | closeness/K^2=%.3f height/K^2=%.3f net/K^2=%.3f | content den: mid-prime part/K^2=%.3f, all/K^2=%.3f | time %.0fs'%(K,float(mlog(abs(val))),(float(mlog(abs(val)))-math.log(H))/K**2,math.log(H)/K**2,float(mlog(abs(val)))/K**2,sum(e*math.log(p) for p,e in mid.items())/K**2,sum(e*math.log(p) for p,e in fd.items())/K**2,time.time()-t0))
print('   den primes:',dict(sorted(fd.items())))
print('   num primes:',dict(sorted(fn.items())))
