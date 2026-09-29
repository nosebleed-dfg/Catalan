from fractions import Fraction as F
from math import gcd
import math, cmath

def saw(x):  # ((x)) for rational x
    if x.denominator == 1: return F(0)
    return x - (x.numerator // x.denominator) - F(1,2)

def ded(h,k):  # carry form: sum ((a/k))((ha/k))
    return sum(saw(F(a,k))*saw(F(h*a,k)) for a in range(1,k))

def ded_recip(h,k):  # Euclid via reciprocity
    h %= k
    if h == 0: return F(0)
    if h == 1: return F((k-1)*(k-2),12*k)
    return (F(h,k)+F(k,h)+F(1,h*k))/12 - F(1,4) - ded_recip(k,h)

def ded_cot(h,k):
    return sum(1/math.tan(math.pi*a/k)/math.tan(math.pi*h*a/k) for a in range(1,k))/(4*k)

def ordmod(b,p):
    x, n = b % p, 1
    while x != 1: x = x*b % p; n += 1
    return n
from sympy import primerange

def cycles(p, b=10):
    seen=set(); cyc=[]
    for a in range(1,p):
        if a in seen: continue
        c=[]; r=a
        while r not in seen:
            seen.add(r); c.append(r); r=r*b%p
        cyc.append(c)
    return cyc

def digitsum(c,p,b=10):  # digits of long division along remainder cycle
    return sum((b*r)//p for r in c)

rows=[]
for p in primerange(3,400):
    if p in (5,): continue
    n=ordmod(10,p)
    orb=sum(ded_recip(pow(10,j,p),p) for j in range(n))
    cyc=cycles(p)
    rhs=sum(F(digitsum(c,p)*2-9*n,18)**2 for c in cyc)
    rows.append((p,n,orb,orb==rhs))
print('identity holds for all:', all(r[3] for r in rows))
print('even period -> 0:', all(r[2]==0 for r in rows if r[1]%2==0))
print('odd period -> >0:', all(r[2]>0 for r in rows if r[1]%2==1))
for p,n,orb,_ in rows:
    if n%2==1: print(p,n,orb)
from sympy import primerange, primitive_root
import cmath
def check(p,b=10):
    g=primitive_root(p); ind={}; x=1
    for e in range(p-1): ind[x]=e; x=x*g%p
    n=ordmod(b,p); m=(p-1)//n
    tot=0; terms=[]
    for t in range(p-1):           # chi(g)=exp(2pi i t/(p-1))
        chi=lambda a: cmath.exp(2j*cmath.pi*t*ind[a%p]/(p-1))
        if abs(chi(p-1)+1)>1e-9: continue      # odd only
        if abs(chi(b)-1)>1e-9: continue        # trivial on <b>
        B=sum(a*chi(a) for a in range(1,p))/p
        terms.append(abs(B)**2)
    orb=sum(ded_recip(pow(b,j,p),p) for j in range(n))
    return p,n,m,float(orb), sum(terms)/m, terms
bad=0
for p in primerange(7,700):
    if p==5: continue
    r=check(p)
    if abs(r[3]-r[4])>1e-6: bad+=1; print('BAD',r[:5])
print('spectral identity mismatches:',bad)
for p in [37,41,79,239,271]:
    r=check(p); print(r[0],'n',r[1],'index',r[2],'orb',r[3],'|B|^2 terms',[round(v,4) for v in r[5]])
from sympy import cyclotomic_poly, factorint, symbols
x=symbols('x')
def orb(k,b=10):
    n=ordmod(b,k); return n, sum(ded_recip(pow(b,j,k),k) for j in range(n))
for n in range(1,12):
    k=int(cyclotomic_poly(n,x).subs(x,10))
    if k<3: continue
    o=orb(k); print(n,k,factorint(k),'period',o[0],'orbit sum',o[1])
from sympy import primerange, primitive_root, Matrix
import cypari2; pari=cypari2.Pari()
def detN(p,b=10):
    n=ordmod(b,p); m=(p-1)//n
    if n%2==0: return None
    g=primitive_root(p)
    reps=[pow(g,i,p) for i in range(m//2)]     # reps of (G/H)/{±1}
    Hset={pow(b,j,p) for j in range(n)}
    def c(a):  # centered remainder sum over the cycle through a
        return sum(saw(F(a*h%p,p)) for h in Hset)
    N=Matrix(m//2,m//2,lambda i,j: c(reps[i]*pow(reps[j],-1,p)%p))
    return n,m,N.det()
def hno(p,d):
    if d==1: return 1
    return int(pari(f'my(v=polsubcyclo({p},{d})); if(type(v)=="t_VEC",v=v[1]); bnfinit(v).no'))
def hminus(p,m):
    return hno(p,m)//hno(p,m//2)
for p in primerange(7,1200):
    r=detN(p)
    if r is None: continue
    n,m,d=r
    if m>16: continue
    hm=hminus(p,m)
    print(p,'period',n,'index',m,'2|det N|=',abs(2*d),'h-=',hm,'OK' if abs(2*d)==hm else 'XX')
import gmpy2, math, time
from mpmath import mp, catalan, mpf, floor, log, pi, sqrt
t=time.time()
mp.dps = 12500
G = +catalan
K = 12000                         # base-9 digits
X = int(floor(G*mpf(9)**K))
s = gmpy2.digits(X,9)             # '0' prefix omitted: G<1
s = s.rjust(K,'0')
print('computed',len(s),'base-9 digits in',round(time.time()-t,1),'s')
print('G in base 9 = 0.'+s[:60]+'...')
from collections import Counter
c=Counter(s); print('digit counts:',[c[str(d)] for d in range(9)],' expected each ~',K//9)
def maxrun(ch):
    best=cur=0; pos=-1
    for i,x in enumerate(s):
        cur=cur+1 if x==ch else 0
        if cur>best: best,pos=cur,i-cur+1
    return best,pos
for ch in '08': print('longest run of',ch,':',maxrun(ch))
print('random expectation for longest run ~ log_9(K*8/9) =',round(math.log(K*8/9,9),2))
open('G_base9.txt','w').write(s)
import math
from mpmath import mp, catalan, mpf, floor, log, pi, sqrt, binomial
mp.dps=3100; K=3000
G=+catalan; s=open('G_base9.txt').read()[:K]
N9=mpf(9)
def prefix(S):          # base-9 digit prefix agreement via monotone binary search
    lo,hi=0,K
    while lo<hi:
        k=(lo+hi+1)//2
        if floor(S*N9**k)==floor(G*N9**k): lo=k
        else: hi=k-1
    return lo
lead=pi/8*log(2+sqrt(3)); tail=mpf(0); term=mpf(1)
events=[]
for n in range(0,4700):
    tail+= term/(2*n+1)**2          # term = 1/C(2n,n)
    term*= mpf(n+1)/(2*(2*n+1))     # 1/C(2n+2,n+1) from 1/C(2n,n)
    if n%3: continue
    S=lead+mpf(3)/8*tail
    err=G-S
    E=int(floor(-log(abs(err),9)))
    T=prefix(S)
    if E-T>=2: events.append((n,E,T,E-T,s[T:E+1]))
print('partial sums checked: every 3rd n up to 4700 (~',int(4700*math.log(4,9)),'base-9 digits)')
print('carry-mistake events (value right to E digits, digit string right to only T):',len(events))
for e in events[:15]: print('n=%d value-digits=%d string-digits=%d gap=%d  G digits there: %s'%e)
import sys; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
s=open('G_base9.txt').read(); K=len(s)
X=int(s,9); lo=F(X,9**K); hi=F(X+1,9**K)
def cf(x,maxn=100000):
    out=[]
    for _ in range(maxn):
        a=x.numerator//x.denominator; out.append(a); x-=a
        if x==0: break
        x=1/x
    return out
a,b=cf(lo),cf(hi); k=0
while k<min(len(a),len(b)) and a[k]==b[k]: k+=1
p0,q0,p1,q1=0,1,1,0
for t in a[:k-1]: p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
print('guaranteed continued-fraction terms of G:',k-1,' largest term:',max(a[:k-1]))
print('first terms:',a[:20])
print('denominator of last certified convergent has',len(str(q1)),'decimal digits')
from mpmath import mp, mpf, exp, pi, binomial, clsin, polylog
mp.dps=30
def s9(n):
    t=0
    while n: t+=n%9; n//=9
    return t
S9=[s9(n) for n in range(6000)]
def F_all(z, K=1300):
    A=[sum((mpf(d)**k if k else 1)*z**d for d in range(9)) for k in range(K+1)]
    zp={}; 
    F={}
    for s in range(8, K+10):
        N=6000 if s<12 else 400
        F[s]=sum(z**S9[n]/mpf(n)**s for n in range(1,N))
    for s in range(7,1,-1):
        rhs=sum(z**d/mpf(d)**s for d in range(1,9))
        rhs+=sum(binomial(-s,k)*mpf(9)**(-s-k)*A[k]*F[s+k] for k in range(1,K+1))
        F[s]=rhs/(1-mpf(9)**(-s)*A[0])
    return F
if __name__=='__main__':
    z=exp(1j*pi/2); F=F_all(z)
    print('check vs Li_s(i):',[mp.nstr(abs(F[s]-polylog(s,z)),3) for s in (2,3,5,9)])
    deg=pi/180
    for name,th in [('90 (Catalan)',90*deg),('45',45*deg),('10/9',mpf(10)/9*deg),('20/9 (2 steps)',mpf(20)/9*deg),('40',40*deg),('60',60*deg),('120',120*deg)]:
        z=exp(1j*th); F=F_all(z); tw=F[2].imag; cl=clsin(2,th)
        print('%-15s twisted=%s  Cl2=%s  carry cost=%s'%(name,mp.nstr(tw,16),mp.nstr(cl,16),mp.nstr(tw-cl,8)))
import numpy as np, math
N=10**8
s=np.ones(N//2,dtype=bool); s[0]=False          # s[i] <-> 2i+1
for i in range(1,int(N**0.5)//2+1):
    if s[i]: p=2*i+1; s[p*p//2::p]=False
P=2*np.nonzero(s)[0]+1; P=P[P>3]                 # drop 2,3 (not coprime to 36)
r=P%36
classes=[a for a in range(36) if math.gcd(a,36)==1]
c={a:sum(1 for x in range(36) if x*x%36==a) for a in classes}
xs=np.unique(np.logspace(4,8,500).astype(np.int64))
idx=np.searchsorted(P,xs,side='right')
pi_x=idx+2
E={a:[] for a in classes}
cum={a:np.concatenate([[0],np.cumsum(r==a)]) for a in classes}
for x,k,px in zip(xs,idx,pi_x):
    for a in classes:
        E[a].append(math.log(x)/math.sqrt(x)*(12*cum[a][k]-px))
print('class a mod 36 | base-9 last digit | Pythagorean? | count to 1e8 | square roots c | predicted mean E (1-c) | observed mean E (log-avg)')
for a in sorted(classes,key=lambda a:(a%9,a%4)):
    print('%2d | %d | %s | %8d | %d | %+d | %+.2f'%(a,a%9,'yes' if a%4==1 else 'no ',cum[a][-1],c[a],1-c[a],np.mean(E[a])))
print()
for d in (1,2,4,5,7,8):
    a1=[a for a in classes if a%9==d and a%4==1][0]; a3=[a for a in classes if a%9==d and a%4==3][0]
    lead=np.mean(np.array(E[a3])>np.array(E[a1]))
    print('last base-9 digit %d: non-Pythagorean (3 mod 4) ahead of Pythagorean at %.1f%% of sampled x'%(d,100*lead))
tot1=sum(cum[a][-1] for a in classes if a%4==1); tot3=sum(cum[a][-1] for a in classes if a%4==3)
print('overall to 1e8: 1 mod 4:',tot1,' 3 mod 4:',tot3)
import numpy as np, math, time, sys
X=int(sys.argv[1]); t0=time.time()
isp=np.ones(X+1,dtype=bool); isp[:2]=False
for p in range(2,int(X**0.5)+1):
    if isp[p]: isp[p*p::p]=False
primes=np.nonzero(isp)[0]
par=np.zeros(X+1,dtype=np.int8); sqf=np.ones(X+1,dtype=bool)
for p in primes:
    p=int(p)
    if p==2 or p%4==1:
        q=p
        while q<=X: par[q::q]^=1; q*=p
        if p*p<=X: sqf[p*p::p*p]=False
    else:
        q=p*p
        while q<=X: par[q::q]^=1; q*=p*p
        if p**4<=X: sqf[p**4::p**4]=False
R=int(X**0.5); ns=[]
for a in range(0,R+1):
    bmax=int((X-a*a)**0.5); b=np.arange(0,bmax+1,dtype=np.int32); n=(a*a+b*b).astype(np.int32)
    mult=np.full(n.shape,4,dtype=np.int8)      # symmetry: (±a,±b)
    if a==0: mult//=2
    mult[b==0]//=2 if a>0 else 1
    if a==0: mult[b==0]=0
    ns.append((n[n>0],mult[n>0]))
n=np.concatenate([x[0] for x in ns]); m=np.concatenate([x[1] for x in ns]).astype(np.int64)
w=np.bincount(n,weights=m*(1-2*par[n].astype(np.int64)),minlength=X+1)
cnt=np.bincount(n,weights=m,minlength=X+1); sq=np.bincount(n,weights=m*sqf[n],minlength=X+1)
L=np.cumsum(w).astype(np.int64); C=np.cumsum(cnt); S=np.cumsum(sq)
print('X=%d  points=%d  time %.0fs'%(X,C[-1],time.time()-t0))
for x in [10**k for k in range(3,int(math.log10(X))+1)]:
    print('x=%9d  L(x)=%8d  L/sqrt(x)=%.3f  Gaussian-squarefree density=%.5f'%(x,L[x],L[x]/math.sqrt(x),S[x]/C[x]))
pos=np.nonzero(L[2:]>0)[0]
print('first x>=2 with L(x)>0:', (pos[0]+2) if len(pos) else 'none', '  max L(x>=2)=',L[2:].max())
r=L[1000:]/np.sqrt(np.arange(1000,X+1)); print('L/sqrt(x) on [1e3,X]: min %.3f max %.3f'%(r.min(),r.max()))
np.save('L_gauss_%d.npy'%X,L[::max(1,X//100000)])
import numpy as np, math, cmath, time
from sympy import primerange
t0=time.time()
def primary(p):
    for b in range(0,int(2*math.sqrt(p))+2,3):
        d=4*p-3*b*b
        if d<0: break
        r=int(math.isqrt(d))
        if r*r==d:
            for a in ((b+r)//2,(b-r)//2):
                if (b+r)%2==0 and a*a-a*b+b*b==p and a%3==2: return a,b
    return None
def powmod_vec(base,e,p):
    res=np.ones_like(base); b=base.copy()
    while e:
        if e&1: res=(res*b)%p
        b=(b*b)%p; e>>=1
    return res
rows=[]
for p in primerange(7,120000):
    if p%3!=1: continue
    ab=primary(p)
    if ab is None: continue
    a,b=ab
    r=(-a*pow(b,-1,p))%p            # omega mod pi
    n=np.arange(1,p,dtype=np.int64)
    c=powmod_vec(n,(p-1)//3,p)     # in {1,r,r^2}
    k=np.where(c==1,0,np.where(c==r,1,2))
    w=np.exp(2j*np.pi*k/3)*np.exp(2j*np.pi*n/p)
    g=w.sum()/math.sqrt(p)
    rows.append((p,p%9,g.real,g.imag,cmath.phase(g)))
rows=np.array(rows)
print('primes p=1 mod 3 up to 120000:',len(rows),' time %.0fs'%(time.time()-t0))
print('check |g|=1: max dev',np.abs(np.hypot(rows[:,2],rows[:,3])-1).max())
for cls in (1,4,7):
    sel=rows[rows[:,1]==cls]
    print('p = %d mod 9: n=%4d  mean Re(g/sqrt p)=%+.4f  mean cos(3 theta)=%+.4f'%(cls,len(sel),sel[:,2].mean(),np.cos(3*sel[:,4]).mean()))
print('all:            mean Re=%+.4f   x^(-1/6) scale = %.3f'%(rows[:,2].mean(),120000**(-1/6)))
# cumulative bias vs x
for x in (10000,30000,60000,120000):
    sel=rows[rows[:,0]<=x]; print('x=%6d  sum Re(g/sqrt p) = %+.2f   n=%d'%(x,sel[:,2].sum(),len(sel)))
np.save('cubic_gauss.npy',rows)
from mpmath import mp, mpf, mpc, zeta, binomial, zetazero, diff, exp
mp.dps=20
def digit_series(b, s0, K=600):
    """D(s)=sum s_b(n) n^{-s} at s0 (complex), via the base-b self-similar recursion."""
    A=[sum((mpf(d)**k if k else 1) for d in range(b)) for k in range(K+2)]
    def sb(n):
        t=0
        while n: t+=n%b; n//=b
        return t
    memo={}
    def D(k):   # returns D(s0+k)
        if k in memo: return memo[k]
        s=s0+k
        if s.real>=12:
            N=400; v=sum(sb(n)*mpf(n)**(-s) for n in range(1,N))
        else:
            rhs=sum(mpf(d)**(1-s) for d in range(1,b))
            rhs+=mpf(b)**(-s)*A[1]*zeta(s)
            for j in range(1,K+1):
                c=binomial(-s,j)*mpf(b)**(-s-j)
                rhs+=c*(A[j]*D(k+j)+A[j+1]*zeta(s+j))
            v=rhs/(1-b*mpf(b)**(-s))
        memo[k]=v; return v
    return D(0)
# sanity: base 9 at s=3 vs direct
s=mpc(3,0)
direct=sum(sum(int(c) for c in __import__('numpy').base_repr(n,9))*mpf(n)**(-3) for n in range(1,200000))
print('check D_9(3): recursion',digit_series(9,s),' direct(2e5 terms)',direct)
for b in (2,9,10):
    print('base',b)
    for j in range(1,4):
        rho=zetazero(j); Dv=digit_series(b,rho)
        zp=diff(zeta,rho)
        v=-1j*Dv/zp          # d rho / d theta at theta=0
        print('  zero %d at %s : d(rho)/d(theta) = %s   (Re = push off the line, Im = slide along it)'%(j, mp.nstr(rho,8), mp.nstr(v,6)))
import sys, math, time, numpy as np
from mpmath import mp, mpf, mpc, zeta, binomial, zetazero, diff
mp.dps=30
def moments(b, s0, thr=4, K=700, N=60000):
    A=[sum((mpf(d)**k if k else 1) for d in range(b)) for k in range(K+3)]
    n=np.arange(1,N,dtype=float); m=np.arange(1,N,dtype=np.int64); sb=np.zeros(N-1,dtype=np.int64)
    while m.any(): sb+=m%b; m//=b
    sb=sb.astype(float); logn=np.log(n)
    s0c=complex(s0)
    zt={}; memo={}
    def Z(k):
        if k not in zt: zt[k]=zeta(s0+k)
        return zt[k]
    def D(k):
        if k in memo: return memo[k]
        s=s0+k
        if s.real>=thr:
            w=np.exp(-(s0c+k)*logn); v=(mpc(np.sum(sb*w)), mpc(np.sum(sb*sb*w)))
        else:
            r1=sum(mpf(d)**(1-s) for d in range(1,b)); r2=sum(mpf(d)**(2-s) for d in range(1,b))
            bs=mpf(b)**(-s); den=1-b*bs
            S1=r1+bs*A[1]*Z(k); S2=r2+bs*A[2]*Z(k)
            for j in range(1,K+1):
                c=binomial(-s,j)*mpf(b)**(-s-j); d1,d2=D(k+j); zj=Z(k+j)
                S1+=c*(A[j]*d1+A[j+1]*zj); S2+=c*(A[j]*d2+2*A[j+1]*d1+A[j+2]*zj)
            D1=S1/den; D2=(S2+bs*2*A[1]*D1)/den
            v=(D1,D2)
        memo[k]=v; return v
    return D(0)
if __name__=="__main__":
    if sys.argv[1]=='check':
        for s in (mpc(3,0),mpc(2.5,3),mpc(1.5,14)):
            print(s, moments(9,s)[0], ' thr=6:',moments(9,s,thr=6)[0])
        r=zetazero(3); print('zero3 D_9 thr4',moments(9,r)[0],' thr6',moments(9,r,thr=6)[0])
        r=zetazero(4); print('zero4 D_5 thr4',moments(5,r)[0],' thr6',moments(5,r,thr=6)[0])
        sys.exit()
    bases=[int(x) for x in sys.argv[1].split(',')]; nz=int(sys.argv[2])
    for b in bases:
        for j in range(1,nz+1):
            t=time.time(); rho=zetazero(j)
            z1=diff(zeta,rho); z2=diff(zeta,rho,2)
            D1,D2=moments(b,rho); h=mpf('1e-2')
            D1p=(moments(b,rho+h)[0]-moments(b,rho-h)[0])/(2*h)
            v=-1j*D1/z1; rpp=-(-D2+2*(1j*D1p)*v+z2*v*v)/z1
            print('%2d | %d | %+.4f %+.4fi | %+.4f %+.4fi | %.0fs'%(b,j,float(v.real),float(v.imag),float(rpp.real),float(rpp.imag),time.time()-t), flush=True)
from mpmath import mp, mpf, mpc, zeta, findroot, nstr
mp.dps=20
def G(s, J=200):   # sum_{j>=1} (zeta(js)-1): Dirichlet series of the all-base carry count per step, divided by zeta(s)
    return sum(zeta(j*s)-1 for j in range(1,J+1))
# real zeros
print('real zeros of G:')
for a,b in ((0.55,0.99),(0.34,0.49),(0.26,0.33),(0.21,0.249)):
    try: print('  in (%.2f,%.2f):'%(a,b), nstr(findroot(lambda s: G(s), (mpf(a),mpf(b)), solver='bisect'),10))
    except Exception as e: print('  ',a,b,'fail',e)
# complex zeros in 0.2<Re<1.6, 0<Im<45 : seed grid + Newton, dedupe
def hunt(f, res, ims, tol=1e-8):
    found=[]
    for re in res:
        for im in ims:
            try: z=findroot(f, mpc(re,im), tol=1e-14, maxsteps=60)
            except Exception: continue
            if not (0.05<z.real<2.5 and 0<z.imag<50): continue
            if abs(f(z))>1e-9: continue
            if all(abs(z-w)>1e-6 for w in found): found.append(z)
    return sorted(found,key=lambda z:z.imag)
import numpy as np
res=list(np.arange(0.3,1.7,0.2)); ims=list(np.arange(1,45,1.5))
zG=hunt(G,res,ims)
z1=hunt(lambda s: zeta(s)-1,res,ims)
print('\ncomplex zeros of G (carry zeros), 0<Im<45:')
for z in zG: print('  ',nstr(z,8))
print('\n1-points of zeta (zeta(s)=1), 0<Im<45:')
for z in z1: print('  ',nstr(z,8))
import sys, numpy as np
from mpmath import mp, mpf, mpc, zeta, findroot, nstr
mp.dps=15
def G(s, J=70): return sum(zeta(j*s)-1 for j in range(1,J+1))
def hunt(f, res, ims):
    found=[]
    for re in res:
        for im in ims:
            try: z=findroot(f, mpc(re,im), tol=1e-12, maxsteps=40)
            except Exception: continue
            if not (0.3<z.real<3 and 0.5<z.imag<90): continue
            if abs(f(z))>1e-8: continue
            if all(abs(z-w)>1e-5 for w in found): found.append(z)
    return sorted(found,key=lambda z:z.imag)
which=sys.argv[1]
res=[0.5,0.8,1.1,1.4]; ims=list(np.arange(42,80,2.0)) if len(sys.argv)>2 else list(np.arange(2,42,2.0))
if which=='G':
    for z in hunt(G,res,ims): print('G-zero ',nstr(z,8), flush=True)
else:
    for z in hunt(lambda s: zeta(s)-1,res,ims): print('1-point',nstr(z,8), flush=True)
import numpy as np, math, sys
from mpmath import clsin, mp
mp.dps=15
def carry_cost(b, k, K=None):
    th=2*math.pi*k/12; z=np.exp(1j*th)
    if K is None: K=int(25*math.log(10)/math.log(b/(b-1)))+50
    d=np.arange(b,dtype=float); zd=z**d
    A=lambda j: np.sum(d**j*zd) if j else np.sum(zd)
    N=8000; n=np.arange(1,N,dtype=np.int64); sb=np.zeros(N-1,dtype=np.int64); m=n.copy()
    while m.any(): sb+=m%b; m//=b
    zs=z**sb; nf=n.astype(float)
    F={}
    for s in range(K+12,1,-1):
        if s>=12: F[s]=np.sum(zs*nf**(-float(s)))
        else:
            rhs=np.sum(zd[1:]*d[1:]**(-float(s)))
            for j in range(1,K+1):
                if s+j>K+12: break
                c=(-1)**j*math.comb(s+j-1,j)*b**(-float(s))
                Aj=np.sum((d/b)**j*zd)
                rhs+=c*Aj*F[s+j]
            F[s]=rhs/(1-b**(-float(s))*A(0))
    tw=F[2].imag; cl=float(clsin(2,th))
    return tw-cl
bases=[int(x) for x in sys.argv[1].split(',')]
ks=[1,2,3,4,5]   # 30,60,90,120,150 degrees (k=0,6 trivial; k>6 mirror)
print('base | gcd(b,12) dead | gcd(b-1,12) invisible | carry cost at 30,60,90,120,150 deg')
for b in bases:
    row=[]
    for k in ks:
        try: row.append(carry_cost(b,k))
        except Exception as e: row.append(float('nan'))
    print('%2d | %2d | %2d | '%(b,math.gcd(b,12),math.gcd(b-1,12))+' '.join('%+.5f'%v for v in row), flush=True)
