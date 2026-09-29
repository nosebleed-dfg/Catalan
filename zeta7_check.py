import sys; sys.set_int_max_str_digits(0)
import sys, time, math
from fractions import Fraction as F
from sympy import bernoulli, Poly, symbols, Rational
import flint
from mpmath import mp, mpf, zeta, log as mlog
t0=time.time()
n=int(sys.argv[1]); K=40*n; N=3*n; h=37*n
x=symbols('x')
# Bernoulli moments mu(t^e)
def mu_mono(e):
    B=F(int(bernoulli(2*e+2).p),int(bernoulli(2*e+2).q))
    return (-1)**e*B*F((2*e+3)*(2*e+4)*(2*e+5)*(2*e+6)*(2*e+7),720)
# D_N^7 as integer polynomial coefficients (low->high) in t
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b): r[i+j]+=ai*bj
    return r
DN=[1]
for j in range(1,N+1): DN=polymul(DN,[j*j,1])
DN7=[1]
for _ in range(7): DN7=polymul(DN7,DN)
tail=list(range(N+1,K+1))            # poles at -j^2
Dtail=[1]
for j in tail: Dtail=polymul(Dtail,[j*j,1])
def polydivmod(num,den):             # exact rational division, coefficients Fractions
    num=[F(c) for c in num]; den=[F(c) for c in den]
    q=[F(0)]*(max(0,len(num)-len(den)+1))
    r=num[:]
    for k in range(len(num)-len(den),-1,-1):
        c=r[k+len(den)-1]/den[-1]; q[k]=c
        for i,d in enumerate(den): r[k+i]-=c*d
    return q, r[:len(den)-1]
def evalpoly(p,v):
    s=0
    for c in reversed(p): s=s*v+c
    return s
# residues: c_{j'} = DN7(-j'^2) (-j'^2)^{i+j} / Dtail'(-j'^2)
Dtp={}
for a in tail:
    prod=1
    for b in tail:
        if b!=a: prod*=(b*b-a*a)
    Dtp[a]=prod
H7={0:F(0)}
for j in range(1,K+1): H7[j]=H7[j-1]+F(1,j**7)
mom={}
def mu_poly(q):
    return sum(c*mu_mono(e) for e,c in enumerate(q) if c)
A=[[F(0)]*h for _ in range(h)]; Bm=[[F(0)]*h for _ in range(h)]
tpow={}
for s in range(2*h-1):
    num=DN7+[0]*0
    num=[0]*s+DN7           # t^s * DN^7
    q,_=polydivmod(num,Dtail)
    polypart=mu_poly(q)
    a_val=polypart; b_val=F(0)
    for jp in tail:
        c=F(evalpoly(DN7,-jp*jp)*((-jp*jp)**s), Dtp[jp])
        a_val+=c*(-(jp**6)*H7[jp]+F(1,2*jp)-F(1,6))
        b_val+=c*(jp**6)
    tpow[s]=(a_val,b_val)
for i in range(h):
    for j in range(h):
        A[i][j],Bm[i][j]=tpow[i+j]
print('entries built %.0fs'%(time.time()-t0), flush=True)
Af=flint.fmpq_mat([[flint.fmpq(v.numerator,v.denominator) for v in row] for row in A])
Bf=flint.fmpq_mat([[flint.fmpq(v.numerator,v.denominator) for v in row] for row in Bm])
# evaluate det(A+XB) at X=0..h and interpolate
vals=[]
for X in range(h+1):
    d=(Af+Bf*flint.fmpq(X)).det(); vals.append(F(int(d.p),int(d.q)))
print('dets done %.0fs'%(time.time()-t0), flush=True)
# Newton interpolation
xs=list(range(h+1)); coef=vals[:]
for k in range(1,h+1):
    for i in range(h,k-1,-1):
        coef[i]=(coef[i]-coef[i-1])/(xs[i]-xs[i-k])
# convert Newton form to monomial
poly=[F(0)]
for k in range(h,-1,-1):
    poly=polymul(poly,[F(-xs[k]),F(1)]) if k<h else [F(0)]
    poly=[F(c) for c in poly]
    poly[0]+=coef[k]
# poly now = sum coef[k] prod_{i<k}(X-x_i)? redo properly:
poly=[coef[h]]
for k in range(h-1,-1,-1):
    poly=polymul(poly,[F(-xs[k]),F(1)]); poly=[F(c) for c in poly]; poly[0]+=coef[k]
assert len(poly)==h+1
lead=poly[-1]
lead_pred=F((-1)**(h*(h-1)//2))
for j in tail: lead_pred*= j**6 * F(evalpoly(DN,-j*j))**7
print('leading coefficient matches Lemma 2.3:', lead==lead_pred)
# S_K
def fact(m): return math.factorial(m)
SK=F(fact(K)**(2*h)*4**(h-1), fact(N)**(16*h))
for i in range(1,h): SK/= fact(2*i)**2
Fpoly=[c*SK for c in poly]
from math import gcd
num_g=0; den_l=1
for c in Fpoly:
    num_g=gcd(num_g,c.numerator); den_l=den_l*c.denominator//gcd(den_l,c.denominator)
content=F(num_g,den_l)
P=[c/content for c in Fpoly]
mp.dps=20000
z7=zeta(7)
def ev(p):
    s=mpf(0)
    for c in reversed(p): s=s*z7+mpf(c.numerator)/mpf(c.denominator)
    return s
FK=ev(Fpoly); PK=ev(P)
print('K=%d h=%d: log F_K(zeta7)/K^2 = %.4f   -log content(F_K)/K^2 = %.4f   log P_K(zeta7) = %.1f   /n^2 = %.1f   P_K>0: %s'%(K,h,float(mlog(abs(FK)))/K**2, -float(mlog(abs(content)))/K**2, float(mlog(abs(PK))), float(mlog(abs(PK)))/n**2, PK>0))
print('max coefficient digits of P_K:', max(len(str(abs(c.numerator))) for c in P), ' time %.0fs'%(time.time()-t0))
