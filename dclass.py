"""dclass.py — one-class Hankel machine on the shifted lattice a = j + OFF/d (j = 0..K-1), full-rate kernel.

Weight w(t) = pi t^2/(2 sinh^2 pi t) (moments |B_{2e+2}|/2 in t), Binet: int w/(t^2+a^2) = (a/2) zeta(2,a) - 1/2 - 1/(4a).
Variable v = d^2 t^2, integer poles -b^2 with b = d*j + OFF;  a = b/d.
  moments  phi(v^e) = d^{2e} |B_{2e+2}|/2,
  pole value phi(1/(v+b^2)) = (b/(2d^3)) zeta(2,b/d) - 1/(2d^2) - 1/(4db),   zeta(2,b/d) = xi - d^2 sum_{0<k<b, k=b mod d} 1/k^2,
  xi := zeta(2, OFF/d)  (d=2: 3 zeta(2) = pi^2/2;  d=4, OFF=1: pi^2 + 8G;  d=3, OFF=1: 2pi^2/3 + (9/2) L(2,chi_-3) = psi'(1/3)).
Entries A + xi B.  Prints closeness/height/net per K^2, the content by residue class of p mod d (ledger), and the rho-law
prediction U(2pi) = -(9/2) ln 2 for every d (spacing 1 at full rate).
usage: python dclass.py d OFF K [N] [m]
"""
import sys, math, time; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, factorint
import flint
from mpmath import mp, mpf, zeta, log as mlog

d=int(sys.argv[1]); OFF=int(sys.argv[2]); K=int(sys.argv[3]); N=int(sys.argv[4]) if len(sys.argv)>4 else 0; m=int(sys.argv[5]) if len(sys.argv)>5 else 0
h=K-N; t0=time.time()
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
def bern(n): B=bernoulli(n); return F(int(B.p),int(B.q))
bb=lambda j: d*j+OFF
labels=[bb(j) for j in range(K)]; head=labels[:N]; tail=labels[N:]
numer=[1]
for x in head: numer=polymul(numer,[x*x,1])
Nm=[1]
for _ in range(m): Nm=polymul(Nm,numer)
D=[1]
for x in tail: D=polymul(D,[x*x,1])
degD=len(D)-1
Dp={a:math.prod(b*b-a*a for b in tail if b!=a) for a in tail}
Nm_at={a:ev(Nm,-a*a) for a in tail}
part={}; acc=F(0)
for b in range(OFF,d*K+1,d): part[b]=acc; acc+=F(1,b*b)
def pv(b): return (-F(b,2*d)*part[b]-F(1,2*d*d)-F(1,4*d*b), F(b,2*d**3))
PV={a:pv(a) for a in tail}
def mu(e): return F(d)**(2*e)*abs(bern(2*e+2))/2
def polydivmod_int(num,den):
    num=list(num); q=[0]*max(0,len(num)-len(den)+1)
    for k in range(len(num)-len(den),-1,-1):
        c=num[k+len(den)-1]; q[k]=c
        for i,dd in enumerate(den): num[k+i]-=c*dd
    r=num[:len(den)-1]; r+=[0]*(len(den)-1-len(r)); return q,r
q,r=polydivmod_int(Nm,D); MU=[mu(e) for e in range(len(q)+2*h+1)]
pw={a:1 for a in tail}; A=[]; B=[]
for s in range(2*h-1):
    a_val=F(0); b_val=F(0)
    for e,c in enumerate(q):
        if c: a_val+=MU[e]*c
    for a in tail:
        c=F(Nm_at[a]*pw[a],Dp[a]); a_val+=c*PV[a][0]; b_val+=c*PV[a][1]
    A.append(a_val); B.append(b_val)
    cs=r[degD-1]; q=[cs]+q; r=[(r[i-1] if i>0 else 0)-cs*D[i] for i in range(degD)]
    for a in tail: pw[a]*=-a*a
Af=flint.fmpq_mat([[flint.fmpq(A[i+j].numerator,A[i+j].denominator) for j in range(h)] for i in range(h)])
Bf=flint.fmpq_mat([[flint.fmpq(B[i+j].numerator,B[i+j].denominator) for j in range(h)] for i in range(h)])
vals=[]
for X in range(h+1):
    dd=(Af+Bf*flint.fmpq(X)).det(); vals.append(F(int(dd.p),int(dd.q)))
xs=list(range(h+1)); coef=vals[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1): coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
poly=[coef[h]]
for k in range(h-1,-1,-1):
    poly=polymul(poly,[F(-xs[k]),F(1)]); poly=[F(c) for c in poly]; poly[0]+=coef[k]
g=0; l=1
for c in poly: g=gcd(g,c.numerator); l=l*c.denominator//gcd(l,c.denominator)
den=l//gcd(l,g); num=g//gcd(l,g)
def factor_small(n,bound):
    f={}
    for p in __import__('sympy').primerange(2,bound):
        if n%p==0:
            e=0
            while n%p==0: n//=p; e+=1
            f[p]=e
    if n!=1: f.update({int(p):int(e) for p,e in factorint(n).items()})
    return f
fd=factor_small(den,d*K+200); fn=factor_small(num,d*K+200)
P=[c/F(num,den) for c in poly]; H=max(abs(c.numerator) for c in P)
logH=math.log(H); logden=sum(e*math.log(p) for p,e in fd.items()); lognum=sum(e*math.log(p) for p,e in fn.items())
flint.ctx.prec=2*H.bit_length()+8*K*K+4096
mp.dps=int(flint.ctx.prec*0.302)+20
xi_mp=zeta(2,mpf(OFF)/d); xi=flint.arb(str(xi_mp))
val=flint.arb(0)
for c in reversed(P): val=val*xi+flint.arb(c.numerator)
lv=float(abs(val).log())
closeness=(lv-logH)/K**2; height=logH/K**2; net=lv/K**2
print('[one class a = j + %d/%d, xi = zeta(2,%d/%d) = %s]  K=%d N=%d m=%d h=%d | closeness/K^2 = %.4f  height/K^2 = %.4f  net/K^2 = %.4f | P(xi)>0: %s | time %.0fs'%(OFF,d,OFF,d,mp.nstr(xi_mp,12),K,N,m,h,closeness,height,net,bool(val>0) if val!=0 else None,time.time()-t0))
print('   rho-law prediction U(2pi) = -(9/2) ln 2 = -3.1192 (N/K -> 0); controls K=40 N=3: zeta2 half -2.905/2.775/-0.130, quarter one-class -2.912/3.835/+0.924')
print('   height = intrinsic %.4f + den %.4f - num %.4f ; den by p mod %d:'%((logH-logden+lognum)/K**2,logden/K**2,lognum/K**2,d),
      ' '.join('%d:%.4f'%(r,sum(e*math.log(p) for p,e in fd.items() if p%d==r)/K**2) for r in range(d)),
      '| primes <= K/2: %.4f, (K/2,dK): %.4f'%(sum(e*math.log(p) for p,e in fd.items() if p<=K/2)/K**2,sum(e*math.log(p) for p,e in fd.items() if p>K/2)/K**2))
print('   den primes (p:e), p>K/2:',' '.join('%d:%d'%(p,e) for p,e in sorted(fd.items()) if p>K/2))
print('   num primes:',dict(sorted(fn.items())))
