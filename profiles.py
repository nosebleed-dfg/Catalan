"""profiles.py — HANDOFF open item 1: exponent profiles e_p/K vs p/K and closeness/height/net vs K,
one pipeline for the three machines.

usage: python profiles.py MACHINE K [N] [m]
  zeta5   : Fauzan weight (moments (-1)^e B_{2e+2}(2e+3)(2e+4)(2e+5)/24), integer poles u=-j^2 (j=1..K),
            numerator D_N^(m-1) (one copy cancels against D_K), pole value j^4(zeta5 - H5_j) + 1/(2j) - 1/4.
            defaults N = round(3K/40), m = 6 (Fauzan's shape; K=40n reproduces zeta5_control.py exactly).
  zeta2   : full-rate weight pi t^2/(2 sinh^2 pi t), v=4t^2, v-poles -(2j+1)^2 (j=0..K-1), numerator E_N^m,
            moments 4^e|B_{2e+2}|/2, pole value (b/16) zeta(2,b/2) - 1/(8b) - 1/8, X = pi^2/2.  defaults N=3, m=0.
  catalan : alternating weight pi t^2 cosh/sinh^2, same poles, moments 4^e(2^{2e+2}-1)|B_{2e+2}|,
            pole value (b/2)(-1)^j (G - beta_j) - 1/(4b), X = G.  defaults N=3, m=0.
Entries are built by the shift recurrence  u^{s+1}N = (u q_s + c_s) D + (u r_s - c_s D)  (O(hK) instead of O(h K^2)).
Writes profile_<machine>_K<K>_N<N>_m<m>.json next to this file.
"""
import sys, math, time, json, os; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd
from sympy import bernoulli, primerange, factorint
import flint

machine=sys.argv[1]; K=int(sys.argv[2])
if machine=='zeta5':
    N=int(sys.argv[3]) if len(sys.argv)>3 else round(3*K/40); m=int(sys.argv[4]) if len(sys.argv)>4 else 6
else:
    N=int(sys.argv[3]) if len(sys.argv)>3 else 3; m=int(sys.argv[4]) if len(sys.argv)>4 else 0
STEP=int(sys.argv[5]) if len(sys.argv)>5 else 2      # zeta2 / catalan only: poles b = STEP*j+1 (STEP even), a = b/2, spacing STEP/2 in a
SUBSET=eval(sys.argv[6]) if len(sys.argv)>6 else None  # zeta2 / catalan only: explicit list of odd pole numerators b (overrides K and STEP)
TAG=sys.argv[7] if len(sys.argv)>7 else 'sub'
if SUBSET is not None: K=len(SUBSET)
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
def bern(n):
    B=bernoulli(n); return F(int(B.p),int(B.q))
sq=lambda x: x*x

if machine=='zeta5':
    labels=list(range(1,K+1))
    def mu(e): return (-1)**e*bern(2*e+2)*F((2*e+3)*(2*e+4)*(2*e+5),24)
    H5={0:F(0)}
    for j in range(1,K+1): H5[j]=H5[j-1]+F(1,j**5)
    def polevalue(j): return (-(j**4)*H5[j]+F(1,2*j)-F(1,4), F(j**4))
    copies=m-1
elif machine=='zeta2':
    labels=sorted(SUBSET) if SUBSET is not None else [STEP*j+1 for j in range(K)]
    assert all(b%2==1 for b in labels)
    def mu(e): return F(4)**e*abs(bern(2*e+2))/2
    beta={}; acc=F(0)
    for b in range(1,max(labels)+2,2): beta[b]=acc; acc+=F(1,b*b)
    def polevalue(b): return (-F(b,4)*beta[b]-F(1,8*b)-F(1,8), F(b,16))
    copies=m
elif machine=='catalan':
    labels=sorted(SUBSET) if SUBSET is not None else [STEP*j+1 for j in range(K)]
    assert all(b%2==1 for b in labels)
    def mu(e): return F(4)**e*(2**(2*e+2)-1)*abs(bern(2*e+2))
    betaalt={}; acc=F(0); sgn=1
    for b in range(1,max(labels)+2,2): betaalt[b]=acc; acc+=F(sgn,b*b); sgn=-sgn
    def polevalue(b):
        sign=(-1)**((b-1)//2)
        return (F(-sign*b,2)*betaalt[b]-F(1,4*b), F(sign*b,2))
    copies=m
elif machine=='zeta2int':
    # same full-rate weight pi t^2/(2 sinh^2 pi t) as 'zeta2', but INTEGER poles u=-j^2 (u=t^2), j=1..K:
    # phi(1/(u+j^2)) = (j/2) zeta(2,j) - 1/2 - 1/(4j) = (j/2)(zeta(2) - H2_{j-1}) - 1/2 - 1/(4j);  moments |B_{2e+2}|/2;  X = pi^2/6
    labels=list(range(1,K+1))
    def mu(e): return abs(bern(2*e+2))/2
    H2={0:F(0)}
    for j in range(1,K+1): H2[j]=H2[j-1]+F(1,j*j)
    def polevalue(j): return (-F(j,2)*H2[j-1]-F(1,2)-F(1,4*j), F(j,2))
    copies=m
else:
    raise SystemExit('machine must be zeta5 | zeta2 | catalan | zeta2int')

head=labels[:N]; tail=labels[N:]
numer=[1]
for x in head: numer=polymul(numer,[sq(x),1])
Nm=[1]
for _ in range(copies): Nm=polymul(Nm,numer)
D=[1]
for x in tail: D=polymul(D,[sq(x),1])
degD=len(D)-1
Dp={a:math.prod(sq(b)-sq(a) for b in tail if b!=a) for a in tail}
Nm_at={a:ev(Nm,-sq(a)) for a in tail}
pv={a:polevalue(a) for a in tail}

def polydivmod_int(num,den):
    num=list(num); q=[0]*max(0,len(num)-len(den)+1)
    for k in range(len(num)-len(den),-1,-1):
        c=num[k+len(den)-1]; q[k]=c
        for i,d in enumerate(den): num[k+i]-=c*d
    r=num[:len(den)-1]; r+=[0]*(len(den)-1-len(r))
    return q,r
q,r=polydivmod_int(Nm,D)
maxe=len(q)+2*h
MU=[mu(e) for e in range(maxe+1)]
pw={a:1 for a in tail}
A=[]; B=[]
for s in range(2*h-1):
    a_val=F(0); b_val=F(0)
    for e,c in enumerate(q):
        if c: a_val+=MU[e]*c
    for a in tail:
        c=F(Nm_at[a]*pw[a],Dp[a]); a_val+=c*pv[a][0]; b_val+=c*pv[a][1]
    A.append(a_val); B.append(b_val)
    cs=r[degD-1]
    q=[cs]+q
    r=[(r[i-1] if i>0 else 0)-cs*D[i] for i in range(degD)]
    for a in tail: pw[a]*=-sq(a)
print('%s K=%d N=%d m=%d h=%d: entries %.0fs'%(machine,K,N,m,h,time.time()-t0), flush=True)

Af=flint.fmpq_mat([[flint.fmpq(A[i+j].numerator,A[i+j].denominator) for j in range(h)] for i in range(h)])
Bf=flint.fmpq_mat([[flint.fmpq(B[i+j].numerator,B[i+j].denominator) for j in range(h)] for i in range(h)])
vals=[]
for X in range(h+1):
    d=(Af+Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p),int(d.q)))
print('  %d determinants %.0fs'%(h+1,time.time()-t0), flush=True)
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
def factor_small(n,bound):
    f={}
    for p in primerange(2,bound):
        if n%p==0:
            e=0
            while n%p==0: n//=p; e+=1
            f[p]=e
    if n!=1: f.update({int(p):int(e) for p,e in factorint(n).items()})
    return f
fd=factor_small(den,4*K+200); fn=factor_small(num,4*K+200)
P=[c/F(num,den) for c in poly]
assert all(c.denominator==1 for c in P)
H=max(abs(c.numerator) for c in P)
logH=math.log(H); logden=sum(e*math.log(p) for p,e in fd.items()); lognum=sum(e*math.log(p) for p,e in fn.items())
def lg(fr): return math.log(abs(fr.numerator))-math.log(fr.denominator) if fr!=0 else float('-inf')
logmax=max(lg(c) for c in poly)

# value at the target, ball arithmetic
flint.ctx.prec=2*H.bit_length()+8*K*K+4096
if machine=='zeta5': xi=flint.arb(5).zeta()
elif machine=='zeta2': xi=flint.arb.pi()**2/2
elif machine=='zeta2int': xi=flint.arb.pi()**2/6
else:
    try: xi=flint.arb.const_catalan()
    except Exception:
        from mpmath import mp, catalan
        mp.dps=int(flint.ctx.prec*0.302)+10; xi=flint.arb(str(+catalan))
val=flint.arb(0)
for c in reversed(P): val=val*xi+flint.arb(c.numerator)
lv=float(abs(val).log())
try: pos=bool(val>0)
except Exception: pos=None
relrad=float(val.rad()/abs(val.mid())) if val.mid()!=0 else float('nan')
closeness=(lv-logH)/K**2; height=logH/K**2; net=lv/K**2
mid={p:e for p,e in fd.items() if K/2<p<2*K}; small={p:e for p,e in fd.items() if p<=K/2}; big={p:e for p,e in fd.items() if p>=2*K}
print('%s K=%d N=%d m=%d h=%d | closeness/K^2=%.4f height/K^2=%.4f net/K^2=%.4f | log P(xi)=%.1f  P>0:%s  ball rel.rad %.1e | time %.0fs'%(machine,K,N,m,h,closeness,height,net,lv,pos,relrad,time.time()-t0))
print('  height = intrinsic %.4f + den %.4f - num %.4f   (per K^2; intrinsic = log max|coeff| of det polynomial)'%(logmax/K**2,logden/K**2,lognum/K**2))
print('  den: small(p<=K/2) %.4f  mid(K/2<p<2K) %.4f  big(p>=2K) %.4f ; num primes %s'%(sum(e*math.log(p) for p,e in small.items())/K**2,sum(e*math.log(p) for p,e in mid.items())/K**2,sum(e*math.log(p) for p,e in big.items())/K**2,dict(sorted(fn.items()))))
print('  e_p/K at p/K (den):', ' '.join('%d:%.2f@%.3f'%(p,e/K,p/K) for p,e in sorted(fd.items())))
out={'machine':machine,'K':K,'N':N,'m':m,'step':STEP,'h':h,'closeness':closeness,'height':height,'net':net,'logP':lv,'positive':pos,
     'log_intrinsic':logmax,'log_den':logden,'log_num':lognum,'den':{str(p):e for p,e in sorted(fd.items())},'num':{str(p):e for p,e in sorted(fn.items())},'time':time.time()-t0}
path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'profile_%s_K%d_N%d_m%d%s%s.json'%(machine,K,N,m,'' if STEP==2 else '_s%d'%STEP,'' if SUBSET is None else '_'+TAG))
json.dump(out,open(path,'w'),indent=1)
print('  written',path)
