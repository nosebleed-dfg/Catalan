"""ledger.py — extreme bookkeeping of the p-adic content of P_K, in base-p (hybrid-base) terms.

For a machine (zeta2 | catalan | zeta2int), K (N=0, m=0) and every prime 3 <= p < 2K:
  * poles written in base p (leading digit d1 = floor(b/p), low digit d0 = b mod p);
  * complement pairs  b + b' = 0 mod p^k  (a carry with zero result digit; Midy / p-complement pairs) and
    coincidence pairs b' - b = 0 mod p^k  (equal digits), both with depth k;
  * v_p(det V) = sum over pairs of [v_p(b'+b) + v_p(b'-b)]   (checked against the exact product prod (b'^2-b^2));
  * per pole the diagonal valuation  delta_b = v_p(D'(-b^2)) + min(v_p A_b, v_p B_b)  and the pole cost sum max(0,-delta_b);
  * von Staudt moments: e with p | denom(mu_e) (Lagrange entries use e <= K-2), naive rank K-1-e each;
  * Lagrange (no-cancellation) prediction  e_naive = 2 v_p(det V) + pole cost + sum of von Staudt ranks,
    the measured e_p (from profile_<machine>_K<K>_N0_m0.json), and the residual = measured - naive
    (negative = Hermite / Midy cancellation), next to the empirical digit rules.
usage: python ledger.py MACHINE K [p_detail]     (p_detail: print the per-pole lines for that prime)
"""
import sys, os, json, math; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from sympy import bernoulli, primerange
machine=sys.argv[1]; K=int(sys.argv[2]); pdet=int(sys.argv[3]) if len(sys.argv)>3 else None
def bern(n): B=bernoulli(n); return F(int(B.p),int(B.q))
def vp(n,p):
    if n==0: return 10**9
    n=abs(n); v=0
    while n%p==0: n//=p; v+=1
    return v
def vpF(x,p): return vp(x.numerator,p)-vp(x.denominator,p)
if machine=='zeta2':
    labels=[2*j+1 for j in range(K)]
    mu=lambda e: F(4)**e*abs(bern(2*e+2))/2
    beta={}; acc=F(0)
    for b in range(1,2*K+2,2): beta[b]=acc; acc+=F(1,b*b)
    pv=lambda b: (-F(b,4)*beta[b]-F(1,8*b)-F(1,8), F(b,16))
elif machine=='catalan':
    labels=[2*j+1 for j in range(K)]
    mu=lambda e: F(4)**e*(2**(2*e+2)-1)*abs(bern(2*e+2))
    betaalt={}; acc=F(0); sgn=1
    for b in range(1,2*K+2,2): betaalt[b]=acc; acc+=F(sgn,b*b); sgn=-sgn
    def pv(b):
        sign=(-1)**((b-1)//2); return (F(-sign*b,2)*betaalt[b]-F(1,4*b), F(sign*b,2))
elif machine=='zeta2int':
    labels=list(range(1,K+1))
    mu=lambda e: abs(bern(2*e+2))/2
    H2={0:F(0)}
    for j in range(1,K+1): H2[j]=H2[j-1]+F(1,j*j)
    pv=lambda j: (-F(j,2)*H2[j-1]-F(1,2)-F(1,4*j), F(j,2))
else: raise SystemExit('machine: zeta2 | catalan | zeta2int')
here=os.path.dirname(os.path.abspath(__file__))
path=os.path.join(here,'profile_%s_K%d_N0_m0.json'%(machine,K))
meas=json.load(open(path))['den'] if os.path.exists(path) else None
if meas is None: print('no measured exponents (run profiles.py %s %d 0 0 first); ledger prints predictions only'%(machine,K))
PVs={b:pv(b) for b in labels}
MU=[mu(e) for e in range(K-1)]
prodV=1
for i,b in enumerate(labels):
    for b2 in labels[i+1:]: prodV*=(b2*b2-b*b)
top=max(labels)
print('%s K=%d N=0: poles %d, top pole %d.  Columns: u=p/K | d1-digit histogram of poles | pairs compl/coinc (with depth) | 2v(detV) [check] | pole cost | von Staudt e:rank | naive | measured | residual | digit rule'%(machine,K,len(labels),top))
rows=[]
for p in primerange(2,2*top+1):
    vV=0; compl=0; coinc=0; vD={b:0 for b in labels}; depth={}
    for i,b in enumerate(labels):
        for b2 in labels[i+1:]:
            s=vp(b2+b,p); d=vp(b2-b,p)
            vV+=s+d; compl+=(s>0); coinc+=(d>0); vD[b]+=s+d; vD[b2]+=s+d
            depth[s+d]=depth.get(s+d,0)+1
    if p==2:
        chk=vp(prodV,2)
        delta={b:vD[b]+min(vpF(PVs[b][0],2),vpF(PVs[b][1],2)) for b in labels}
        cost=sum(max(0,-delta[b]) for b in labels)
        vs=[(e,vpF(MU[e],2)) for e in range(K-1) if vpF(MU[e],2)<0]
        m=meas.get('2',0) if meas is not None else None
        print('p=  2 (base-2 ledger) | pairs by v_2(b\'^2-b^2) = k: %s | 2v(detV)=%d [%s] | pole cost %d | moments with v_2<0: %d (min %d) | measured e_2 = %s | e_2/K^2 = %s | 2v/K^2 = %.3f'%(
            ' '.join('%d:%d'%(k,c) for k,c in sorted(depth.items())),2*vV,'ok' if chk==vV else 'MISMATCH',cost,len(vs),min([v for e,v in vs],default=0),m,('%.4f'%(m/K**2)) if m is not None else '-',2*vV/K**2))
        if pdet==2:
            print('   per-pole 2-adic ledger (b | b in base 2 | v_2 D\' | v_2 A_b, v_2 B_b | delta)')
            for b in labels: print('   %4d | %10s | %3d | %3d,%3d | %3d'%(b,bin(b)[2:],vD[b],vpF(PVs[b][0],2),vpF(PVs[b][1],2),delta[b]))
            print('   moments v_2(mu_e), e=0..K-2: '+' '.join('%d:%d'%(e,vpF(MU[e],2)) for e in range(K-1)))
        continue
    chk=vp(prodV,p)
    delta={b:vD[b]+min(vpF(PVs[b][0],p),vpF(PVs[b][1],p)) for b in labels}
    cost=sum(max(0,-delta[b]) for b in labels)
    vs=[(e,K-1-e) for e in range(K-1) if vpF(MU[e],p)<0]
    vsrank=sum(r for e,r in vs)
    naive=2*vV+cost+vsrank
    m=meas.get(str(p),0) if meas is not None else None
    hist={}
    for b in labels: hist[b//p]=hist.get(b//p,0)+1
    if machine=='zeta2': rule=2*sum(1 for n in range(2*K) if (n//p)%2==1)-(1 if p<K else 0)
    elif machine=='catalan':
        n=sum(1 for b in labels if b>p); rule=min(3*n+1,2*K-2)
    else: rule=2*vV
    if p>K/2 or pdet==p:
        print('p=%3d u=%.3f | d1:%-18s | pairs %3d/%-3d | 2v=%4d [%s] | cost %3d | vS %-16s | naive %4d | meas %4s | resid %5s | rule %4d'%(
            p,p/K,' '.join('%d:%d'%(k,v) for k,v in sorted(hist.items())),compl,coinc,2*vV,'ok' if chk==vV else 'MISMATCH',cost,
            ' '.join('%d:%d'%(e,r) for e,r in vs) if vs else '-',naive,m if m is not None else '-',(m-naive) if m is not None else '-',rule))
    if pdet==p:
        print('   per-pole ledger for p=%d (b | base-p digits d1,d0 | complement partners (depth) | coincidence partners (depth) | v_p D\' | v_p A_b, v_p B_b | delta | cost)'%p)
        for b in labels:
            cp=[(b2,vp(b+b2,p)) for b2 in labels if b2!=b and (b+b2)%p==0]
            cc=[(b2,vp(abs(b2-b),p)) for b2 in labels if b2!=b and (b2-b)%p==0]
            print('   %4d | %2d,%3d | %-26s | %-26s | %2d | %3d,%3d | %3d | %d'%(b,b//p,b%p,' '.join('%d(%d)'%t for t in cp),' '.join('%d(%d)'%t for t in cc),vD[b],vpF(PVs[b][0],p),vpF(PVs[b][1],p),delta[b],max(0,-delta[b])))
        print('   moments with p in the denominator (e: v_p mu_e): '+' '.join('%d:%d'%(e,vpF(MU[e],p)) for e in range(K-1) if vpF(MU[e],p)<0))
