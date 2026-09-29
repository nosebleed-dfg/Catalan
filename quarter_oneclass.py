"""One-unknown quarter-integer functional.

Plain (full-rate) kernel w(t) = pi t^2/(2 sinh^2(pi t)), poles at a single class a = j + OFF/4 (OFF = 1 or 3),
j = 0..K-1: spacing 1 at full rate, like the zeta(2) half-integer control, but the pole values are affine in
    xi = zeta(2, OFF/4) = pi^2 + 8G (OFF=1)  or  pi^2 - 8G (OFF=3)
alone (Hurwitz recurrence zeta(2,a) - zeta(2,a+1) = 1/a^2).  Working variable v = 16 t^2, v-poles at -b^2, b = 4j+OFF.
Pole value: phi(1/(v+b^2)) = (b/128) zeta(2,b/4) - 1/32 - 1/(16b),  zeta(2,b/4) = xi - 16 sum_{0<k<b, k=b mod 4} 1/k^2.
Moments: phi(v^e) = 16^e |B_{2e+2}|/2.   Entries A + xi B.   Output per K^2 as in the other scripts.

usage: python quarter_oneclass.py K N m OFF
"""
import sys, math, time; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, factorint
import flint
from mpmath import mp, mpf, catalan, log as mlog, quad, sinh, inf

K=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); OFF=int(sys.argv[4]); h=K-N
t0=time.time()
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b): r[i+j]+=ai*bj
    return r
bb=lambda j: 4*j+OFF
EN=[1]
for j in range(N): EN=polymul(EN,[bb(j)**2,1])
ENm=[1]
for _ in range(m): ENm=polymul(ENm,EN)
tail=list(range(N,K)); Et=[1]
for j in tail: Et=polymul(Et,[bb(j)**2,1])
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
Etp={a:math.prod(bb(b)**2-bb(a)**2 for b in tail if b!=a) for a in tail}
def phi_mono(e):
    B=bernoulli(2*e+2); return F(16)**e*abs(F(int(B.p),int(B.q)))/2
ENm_at={j:ev(ENm,-bb(j)**2) for j in tail}
part={}; acc=F(0)
for b in range(OFF,4*K+1,4): part[b]=acc; acc+=F(1,b*b)
tp={}
for s in range(2*h-1):
    q=polydiv([0]*s+ENm,Et)
    a_val=sum(c*phi_mono(e) for e,c in enumerate(q) if c); b_val=F(0)
    for j in tail:
        b=bb(j); c=F(ENm_at[j]*((-b*b)**s), Etp[j])
        a_val+=c*(-F(b,8)*part[b]-F(1,32)-F(1,16*b)); b_val+=c*F(b,128)
    tp[s]=(a_val,b_val)

# entry check against the direct integral (positive integrand; exact side needs the cancellation digits)
def w(t): return mp.pi*t*t/(2*sinh(mp.pi*t)**2)
worst=0
for s in (0,2*h-2):
    A,B=tp[s]
    mp.dps=30
    direct=quad(lambda t: w(t)*ev(ENm,16*t*t)*(16*t*t)**s/ev(Et,16*t*t),[0,0.5,1,2,3,4,6,8,10,12,15,18,22,27,33,40,50,inf])
    mp.dps=100+3*K+max(len(str(x.numerator))+len(str(x.denominator)) for x in (A,B))
    xi=mp.pi**2+(8 if OFF==1 else -8)*catalan
    exact=mpf(A.numerator)/A.denominator+xi*mpf(B.numerator)/B.denominator
    worst=max(worst,abs((direct-exact)/exact))
mp.dps=30
print('entry check vs direct integral (s=0,2h-2): worst relative error %s'%mp.nstr(worst,3), flush=True)

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
g=0; l=1
for c in poly:
    g=gcd(g,c.numerator); l=l*c.denominator//gcd(l,c.denominator)
den=l//gcd(l,g); num=g//gcd(l,g)
fd=factorint(den); fn=factorint(num)
mid={p:e for p,e in fd.items() if K/2<p<2*K}
small={p:e for p,e in fd.items() if p<=K/2}
big={p:e for p,e in fd.items() if p>=2*K}
logden=sum(e*math.log(p) for p,e in fd.items())
P=[c/F(num,den) for c in poly]
H=max(abs(c.numerator) for c in P)
mp.dps=max(3000,3*len(str(H)))
xi=mp.pi**2+(8 if OFF==1 else -8)*catalan
val=mpf(0)
for c in reversed(P): val=val*xi+mpf(c.numerator)/mpf(c.denominator)
lv=float(mlog(abs(val))); lH=math.log(H)
print('[quarter-integer one class a=j+%d/4, xi = pi^2 %s 8G] K=%d N=%d m=%d h=%d | closeness/K^2 = %.3f  height/K^2 = %.3f  net/K^2 = %.3f | P(xi)>0: %s | time %.0fs'%(OFF,'+' if OFF==1 else '-',K,N,m,h,(lv-lH)/K**2,lH/K**2,lv/K**2,val>0,time.time()-t0))
print('   controls (K=40): zeta(5) -3.457/3.291/-0.166 ; zeta(2) half-poles m=0 -2.905/2.775/-0.130, m=1 -3.035/2.891/-0.145 ; Catalan m=1 -2.15/2.53/+0.38')
print('   content: log(den)/K^2=%.3f  small(p<=K/2)/K^2=%.3f  mid(K/2<p<2K)/K^2=%.3f  big(p>=2K)/K^2=%.3f'%(logden/K**2,sum(e*math.log(p) for p,e in small.items())/K**2,sum(e*math.log(p) for p,e in mid.items())/K**2,sum(e*math.log(p) for p,e in big.items())/K**2))
print('   den primes:',dict(sorted(fd.items())))
print('   num primes:',dict(sorted(fn.items())))
