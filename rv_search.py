"""rv_search.py — Rhin–Viola-type search over the very-well-poised family for Catalan's constant.

Parameters (Zudilin math/0210423, (15)-(16)): h0, h4 in Z, h1, h2, h3 in Z+1/2, all positive, 1+h0-hj-hl > 0.
Directions: h0 = n (normalise alpha = 1), hj = eta_j n.  Derived 10-set c (per n):
  c00 = 1 - eta3 - eta4,  c11 = 1 - eta1 - eta2,  c21 = eta3,  c31 = eta4,  c12 = eta2,  c13 = eta1,
  c22 = 1 - eta1 - eta3,  c23 = 1 - eta2 - eta3,  c32 = 1 - eta1 - eta4,  c33 = 1 - eta2 - eta4.
Linear form  H(c) = iint x^{c21}(1-x)^{c22} y^{c31}(1-y)^{c33} / (1-xy)^{c11+1} dx dy  in Q G + Q.
  decay rate  = max_{(x,y) in (0,1)^2} [c21 ln x + c22 ln(1-x) + c31 ln y + c33 ln(1-y) - c11 ln(1-xy)]   (Laplace),
  denominator rate (Zudilin's expectation (24), exact on his direction) = 2 M ln 2 + m1 + m2,  M = c22 + c31,  m1 >= m2 the two largest of 2c.
  miss = denominator rate + decay rate  (irrationality needs miss < 0).   Zudilin's direction eta = (1,1,1,1)/3: miss = 4 ln2 + 4 - 2.406 = 4.37.
usage: python rv_search.py
"""
import math, itertools
import numpy as np
xs=np.linspace(0.002,0.998,300); X,Y=np.meshgrid(xs,xs,indexing='ij')
LX,L1X,LY,L1Y,LXY=np.log(X),np.log(1-X),np.log(Y),np.log(1-Y),np.log(1-X*Y)
def cset(e):
    e1,e2,e3,e4=e
    return {'00':1-e3-e4,'11':1-e1-e2,'21':e3,'31':e4,'12':e2,'13':e1,'22':1-e1-e3,'23':1-e2-e3,'32':1-e1-e4,'33':1-e2-e4}
def decay(c):
    f=c['21']*LX+c['22']*L1X+c['31']*LY+c['33']*L1Y-c['11']*LXY
    i=np.unravel_index(np.argmax(f),f.shape); x0,y0=X[i],Y[i]
    # local refinement (Newton-free: coordinate golden search)
    def F(x,y): return c['21']*math.log(x)+c['22']*math.log(1-x)+c['31']*math.log(y)+c['33']*math.log(1-y)-c['11']*math.log(1-x*y)
    h=0.01
    for _ in range(40):
        best=(F(x0,y0),x0,y0)
        for dx,dy in ((h,0),(-h,0),(0,h),(0,-h),(h,h),(-h,-h),(h,-h),(-h,h)):
            x,y=x0+dx,y0+dy
            if 0<x<1 and 0<y<1:
                v=F(x,y)
                if v>best[0]: best=(v,x,y)
        if best[1]==x0 and best[2]==y0: h*=0.5
        else: _,x0,y0=best
    return best[0],x0,y0
def denom(c):
    vals=sorted((2*v for v in c.values()),reverse=True); M=c['22']+c['31']
    return 2*M*math.log(2)+vals[0]+vals[1], M, vals[0], vals[1]
def miss(e):
    c=cset(e)
    if min(c.values())<0: return None
    d,x0,y0=decay(c); dn,M,m1,m2=denom(c)
    return dn+d,d,dn,M,m1,m2,x0,y0
# Zudilin's direction
e0=(1/3,1/3,1/3,1/3); r=miss(e0)
print('Zudilin direction eta=(1,1,1,1)/3 (alpha=1): decay %.4f  den %.4f (M=%.3f, m1=%.3f, m2=%.3f)  miss %.4f  saddle (%.4f,%.4f)   [x3 for his normalisation: decay %.3f, miss %.3f]'%(r[1],r[2],r[3],r[4],r[5],r[0],r[6],r[7],3*r[1],3*r[0]))
# coarse grid over eta in (0,1/2]^4 with the positivity constraints, then refine the best
grid=np.linspace(0.02,0.5,13); best=[]
for e in itertools.product(grid,repeat=4):
    r=miss(e)
    if r is None: continue
    best.append((r[0]/max(1e-9,sum(e)),r,e))   # normalise by the total parameter size so directions are comparable per 'n'
best.sort(key=lambda t:t[0])
print('\nbest 10 directions by miss per unit of sum(eta) (coarse grid):')
for s,r,e in best[:10]:
    print('  eta=%s  miss/sum=%.4f   (raw miss %.4f = den %.4f + decay %.4f; M=%.3f m1=%.3f m2=%.3f; saddle %.3f,%.3f)'%(tuple(round(v,3) for v in e),s,r[0],r[2],r[1],r[3],r[4],r[5],r[6],r[7]))
# refine around the best
e=np.array(best[0][2]); step=0.03
for it in range(60):
    improved=False
    for k in range(4):
        for sgn in (1,-1):
            f=e.copy(); f[k]+=sgn*step
            if f.min()<=0.005: continue
            r=miss(tuple(f))
            if r is None: continue
            if r[0]/sum(f)<best[0][0]-1e-12:
                best[0]=(r[0]/sum(f),r,tuple(f)); e=f; improved=True
    if not improved: step*=0.5
    if step<1e-4: break
s,r,e=best[0]
print('\nrefined optimum: eta=%s  miss/sum(eta)=%.4f  raw: den %.4f (2M ln2 = %.4f, m1+m2 = %.4f) + decay %.4f = %.4f'%(tuple(round(v,4) for v in e),s,r[2],2*r[3]*math.log(2),r[4]+r[5],r[1],r[0]))
print('per-unit ledger (divide by sum(eta)=%.4f): 2-adic %.4f, primes %.4f, decay %.4f'%(sum(e),2*r[3]*math.log(2)/sum(e),(r[4]+r[5])/sum(e),r[1]/sum(e)))
