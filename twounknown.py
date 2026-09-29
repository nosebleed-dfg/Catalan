"""Two-unknown full-field functional (HANDOFF open item 1).

Plain (full-rate) kernel, weight w(t) = pi t^2 / (2 sinh^2(pi t))  (moments |B_{2e+2}|/2 in t),
poles at the quarter-integers a = b/4, b odd, i.e. (j+1/4)^2 and (j+3/4)^2.
Working variable v = 16 t^2, so the poles sit at v = -b^2 (b odd), same pole set as zeta2_halfpoles.py.

Pole value (Binet, exact for every a>0):  int_0^inf w(t)/(t^2+a^2) dt = (a/2) zeta(2,a) - 1/2 - 1/(4a)
  => phi(1/(v+b^2)) = (b/128) zeta(2,b/4) - 1/32 - 1/(16b),
     zeta(2,b/4) = pi^2 + 8 chi_{-4}(b) G - 16 sum_{0<k<b, k=b mod 4} 1/k^2.
  => every entry is A + pi^2 B + G C with A,B,C rational.
Moments: phi(v^e) = 16^e |B_{2e+2}| / 2.

det(A + X B + Y C) is a polynomial of total degree h in (X,Y); recovered exactly from the
triangular grid {(a,b): a+b<=h} by 2-D forward differences.  Reported per K^2 at (X,Y)=(pi^2,G).

usage: python twounknown.py K N m
"""
import sys, math, time; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from math import gcd, factorial
from sympy import bernoulli, factorint
import flint
from mpmath import mp, mpf, catalan, log as mlog, quad, sinh, inf

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
def phi_mono(e):
    B=bernoulli(2*e+2); return F(16)**e*abs(F(int(B.p),int(B.q)))/2
ENm_at={j:ev(ENm,-odd(j)**2) for j in tail}
part={}; s1=F(0); s3=F(0)
for b in range(1,2*K+2,2):
    if b%4==1: part[b]=s1; s1+=F(1,b*b)
    else: part[b]=s3; s3+=F(1,b*b)
chi=lambda b: 1 if b%4==1 else -1
tp={}
for s in range(2*h-1):
    q=polydiv([0]*s+ENm,Et)
    a_val=sum(c*phi_mono(e) for e,c in enumerate(q) if c); b_val=F(0); c_val=F(0)
    for j in tail:
        b=odd(j); c=F(ENm_at[j]*((-b*b)**s), Etp[j])
        a_val+=c*(-F(b,8)*part[b]-F(1,32)-F(1,16*b))
        b_val+=c*F(b,128)
        c_val+=c*F(chi(b)*b,16)
    tp[s]=(a_val,b_val,c_val)
print('entries built %.0fs'%(time.time()-t0), flush=True)

# --- check a few entries against the direct integral of w(t) R(16 t^2) ---
def w(t): return mp.pi*t*t/(2*sinh(mp.pi*t)**2)
worst=0
for s in (0,1,h,2*h-2):
    A,B,C=tp[s]
    mp.dps=30      # the integrand is positive: no cancellation, 30 digits are 30 digits
    direct=quad(lambda t: w(t)*ev(ENm,16*t*t)*(16*t*t)**s/ev(Et,16*t*t),[0,0.5,1,2,3,4,6,8,10,12,15,18,22,27,33,40,50,inf])
    # A, B, C cancel by hundreds of digits at K=40+: evaluate the exact side well beyond that
    mp.dps=100+3*K+max(len(str(x.numerator))+len(str(x.denominator)) for x in (A,B,C))
    exact=mpf(A.numerator)/A.denominator+mp.pi**2*mpf(B.numerator)/B.denominator+catalan*mpf(C.numerator)/C.denominator
    rel=abs((direct-exact)/exact); worst=max(worst,rel)
mp.dps=30
print('entry check vs direct integral (s=0,1,h,2h-2): worst relative error %s'%mp.nstr(worst,3), flush=True)

# --- determinants on the triangular grid ---
Af=flint.fmpq_mat([[flint.fmpq(tp[i+j][0].numerator,tp[i+j][0].denominator) for j in range(h)] for i in range(h)])
Bf=flint.fmpq_mat([[flint.fmpq(tp[i+j][1].numerator,tp[i+j][1].denominator) for j in range(h)] for i in range(h)])
Cf=flint.fmpq_mat([[flint.fmpq(tp[i+j][2].numerator,tp[i+j][2].denominator) for j in range(h)] for i in range(h)])
vals={}
for a in range(h+1):
    for b in range(h+1-a):
        d=(Af+Bf*flint.fmpq(a)+Cf*flint.fmpq(b)).det(); vals[(a,b)]=F(int(d.p),int(d.q))
print('%d determinants done %.0fs'%(len(vals),time.time()-t0), flush=True)

# --- 2-D forward differences -> falling-factorial coefficients -> monomials ---
dX={}
for b in range(h+1):
    cur=[vals[(a,b)] for a in range(h+1-b)]; dX[(0,b)]=cur[0]
    for i in range(1,len(cur)+0):
        cur=[cur[k+1]-cur[k] for k in range(len(cur)-1)]; dX[(i,b)]=cur[0]
c={}
for i in range(h+1):
    cur=[dX[(i,b)] for b in range(h+1-i)]; c[(i,0)]=cur[0]/factorial(i)
    for j in range(1,len(cur)+0):
        cur=[cur[k+1]-cur[k] for k in range(len(cur)-1)]; c[(i,j)]=cur[0]/(factorial(i)*factorial(j))
def newton_to_monomial(coefs):     # nodes 0,1,2,...; coefs = divided differences
    n=len(coefs)-1; poly=[coefs[n]]
    for k in range(n-1,-1,-1):
        poly=polymul(poly,[F(-k),F(1)]); poly=[F(x) for x in poly]; poly[0]+=coefs[k]
    return poly
q={}
for j in range(h+1):
    mono=newton_to_monomial([c[(i,j)] for i in range(h+1-j)])
    for ip,v in enumerate(mono): q[(ip,j)]=v
p={}
for ip in range(h+1):
    mono=newton_to_monomial([q[(ip,j)] for j in range(h+1-ip)])
    for jp,v in enumerate(mono): p[(ip,jp)]=v
print('interpolation done %.0fs'%(time.time()-t0), flush=True)

# --- primitive part, height, factorization of the content ---
g=0; l=1
for cc in p.values():
    g=gcd(g,cc.numerator); l=l*cc.denominator//gcd(l,cc.denominator)
den=l//gcd(l,g); num=g//gcd(l,g)
P={k:cc/F(num,den) for k,cc in p.items()}
assert all(cc.denominator==1 for cc in P.values())
H=max(abs(cc.numerator) for cc in P.values())
fd=factorint(den); fn=factorint(num)
mid={pp:e for pp,e in fd.items() if K/2<pp<2*K}; small={pp:e for pp,e in fd.items() if pp<=K/2}
logden=sum(e*math.log(pp) for pp,e in fd.items())

# --- structure ---
nz={k:v for k,v in P.items() if v!=0}
degX=max(i for i,j in nz); degY=max(j for i,j in nz); degT=max(i+j for i,j in nz)
oddY=max([abs(v.numerator) for (i,j),v in nz.items() if j%2==1],default=0)
evenY=max([abs(v.numerator) for (i,j),v in nz.items() if j%2==0],default=0)
onlyX=max([abs(v.numerator) for (i,j),v in nz.items() if j==0],default=0)
onlyY=max([abs(v.numerator) for (i,j),v in nz.items() if i==0],default=0)

# --- evaluation at (pi^2, G) and at the mirror points ---
mp.dps=max(3000,3*len(str(H)))
X=mp.pi**2; Y=+catalan
def evalP(X,Y):
    tot=mpf(0)
    for ip in range(h,-1,-1):
        inner=mpf(0)
        for jp in range(h-ip,-1,-1): inner=inner*Y+mpf(P[(ip,jp)].numerator)
        tot=tot*X+inner
    return tot
val=evalP(X,Y); val_m=evalP(X,-Y); val_x=evalP(X,mpf(0)); val_y=evalP(mpf(0),Y)
lv=float(mlog(abs(val))); lH=math.log(H)
print('K=%d N=%d m=%d h=%d | closeness/K^2 = %.3f  height/K^2 = %.3f  net/K^2 = %.3f | P(pi^2,G)>0: %s | time %.0fs'%(K,N,m,h,(lv-lH)/K**2,lH/K**2,lv/K**2,val>0,time.time()-t0))
print('   controls (K=40): zeta(5) -3.457/3.291/-0.166 ; zeta(2) half-poles m=0 -2.905/2.775/-0.130 ; Catalan m=1 -2.15/2.53/+0.38')
print('   content: log(den)/K^2=%.3f  mid-prime part/K^2=%.3f  small-prime part/K^2=%.3f'%(logden/K**2,sum(e*math.log(pp) for pp,e in mid.items())/K**2,sum(e*math.log(pp) for pp,e in small.items())/K**2))
print('   den primes:',dict(sorted(fd.items())))
print('   num primes:',dict(sorted(fn.items())))
print('   structure: deg_X(pi^2)=%d deg_Y(G)=%d total=%d  nonzero monomials %d of %d'%(degX,degY,degT,len(nz),len(P)))
print('   log10 max|coeff|: all %.1f | odd in G %.1f | even in G %.1f | pure pi^2 (G^0) %.1f | pure G (pi^0) %.1f'%(math.log10(H),math.log10(oddY) if oddY else float('-inf'),math.log10(evenY),math.log10(onlyX) if onlyX else float('-inf'),math.log10(onlyY) if onlyY else float('-inf')))
print('   log|P|/K^2 at (pi^2,G) %.3f | (pi^2,-G) %.3f | (pi^2,0) %.3f | (0,G) %.3f'%(lv/K**2,float(mlog(abs(val_m)))/K**2,float(mlog(abs(val_x)))/K**2,float(mlog(abs(val_y)))/K**2))

# --- independent check: numeric determinant of A + pi^2 B + G C (arb) vs content * P(pi^2,G) ---
try:
    flint.ctx.prec=max(20000,30*K*K)
    pi2=flint.arb.pi()**2; Gc=flint.arb(str(+catalan))
    def toarb(fr): return flint.arb(fr.numerator)/flint.arb(fr.denominator)
    M=flint.arb_mat([[toarb(tp[i+j][0])+pi2*toarb(tp[i+j][1])+Gc*toarb(tp[i+j][2]) for j in range(h)] for i in range(h)])
    d=M.det(); ld=float(abs(d).log())
    print('   arb check: log|det(A+pi^2 B+G C)| = %.6f  vs  log|content*P(pi^2,G)| = %.6f'%(ld,float(mlog(abs(F(num,den))))+lv))
except Exception as e:
    print('   arb check failed:',repr(e))

# --- dump the primitive polynomial ---
import os
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'P_twounknown_K%d_N%d_m%d.txt'%(K,N,m))
with open(out,'w') as f:
    f.write('# coefficient of (pi^2)^i G^j ; content = %d/%d\n'%(num,den))
    for (i,j),v in sorted(P.items()): f.write('%d %d %d\n'%(i,j,v.numerator))
print('   polynomial written to',out)
