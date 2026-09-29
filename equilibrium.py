"""equilibrium.py — the analytic side U(rho) from potential theory (Brown, Prop. 3.5 = Heine's formula).

det M_K = (1/K!) int prod_{i<j} (v_i - v_j)^2 prod dmu(v_i),  dmu = w(t)/D(4t^2) dt,  v = 4t^2.
With x = t/K and K poles a_j = s(j+1/2) (spacing s in a-units), w ~ t^2 e^{-rate t}:
  log det M_K = -K^2 min_nu E[nu] + o(K^2),
  E[nu] = iint [-log|x-y| - log(x+y)] dnu dnu + int [rate*x + (1/s) int_0^s log(x^2+y^2) dy] dnu,   nu prob. measure on [0,inf).
Scaling x -> x/s shows E depends on rho = rate*s only.  Symmetrising nu to an even measure on [-R,R] turns this into
the classical weighted log-energy 2[ I(nu~) + int Q dnu~ ],  Q(x) = (rho/2)|x| + (1/2) g(|x|),  g(x) = log(1+x^2) - 2 + 2x arctan(1/x)
( = int_0^1 log(x^2+y^2) dy ).  Chebyshev inversion of the log kernel on [-R,R] gives, with Q(Ru) = sum_k q_k T_k(u):
  density  psi(Ru) R = [1/pi + sum_{k>=2} c_k T_k(u)] / sqrt(1-u^2),  c_k = -k q_k/(2 pi);   soft edges  <=>  sum_{k>=2} k q_k = 2  (fixes R);
  E_min = 2 q_0 - 2 log(R/2) - (1/4) sum_{k>=2} k q_k^2.
The |x| part is done analytically (|u| = 2/pi - (4/pi) sum_m (-1)^m T_{2m}(u)/(4m^2-1); its sum k q_k = rho R/pi, its sum k q_k^2 = (rho R/pi)^2),
the smooth part g by a cosine transform.  Closeness_inf = log det M(xi)/K^2 - intrinsic/K^2 -> -E_min  (intrinsic -> 0 for the half-integer machines).
usage: python equilibrium.py [rho1 rho2 ...]   (multiples of pi may be written like 2pi)
"""
import sys, math
import numpy as np
def Gab(x,a,b):
    """int_a^b log(x^2+y^2) dy, vectorised in x>=0 (x=0 handled)."""
    x=np.abs(np.asarray(x,dtype=float))
    def F(y):
        out=np.where(y>0, y*np.log(x**2+y**2+1e-300)-2*y, 0.0)
        with np.errstate(divide='ignore',invalid='ignore'):
            at=np.where(x>0, 2*x*np.arctan(np.where(x>0,y/np.where(x>0,x,1),0)), 0.0)
        return out+at
    return F(b)-F(a)
ALPHA=0.0; MEFF=0      # head fraction N/K and numerator power: density ∝ E_head^MEFF / D_tail, h = K-N particles
def g(x):
    """pole part of the external field in h-units: tail poles on (t1,t2) = (a/(1-a), 1/(1-a)), head poles on (0,t1) with power -MEFF"""
    t1=ALPHA/(1-ALPHA); t2=1.0/(1-ALPHA)
    return Gab(x,t1,t2)-MEFF*Gab(x,0.0,t1)
M=1<<14; theta=np.pi*(np.arange(M)+0.5)/M; u=np.cos(theta)
KMAX=4000; ks=np.arange(2,KMAX+1,2)
def cheb_g(R):
    """Chebyshev coefficients (k=0 and even k>=2) of (1/2) g(R u)."""
    f=0.5*g(R*u)
    q0=f.mean()
    qk=np.array([2*np.mean(f*np.cos(k*theta)) for k in ks])
    return q0,qk
def solve(rho, verbose=True):
    def edge(R):
        q0,qk=cheb_g(R); return rho*R/np.pi+np.sum(ks*qk)-2.0
    lo,hi=1e-3,50.0
    flo=edge(lo)
    for _ in range(200):
        mid=0.5*(lo+hi); fm=edge(mid)
        if (fm>0)==(flo>0): lo,flo=mid,fm
        else: hi=mid
        if hi-lo<1e-12: break
    R=0.5*(lo+hi)
    q0g,qkg=cheb_g(R)
    m=ks//2; qkx=-(2*rho*R/np.pi)*((-1.0)**m)/(4*m**2-1)     # |x| part: (rho/2) R |u|
    q0=rho*R/np.pi+q0g
    S=(rho*R/np.pi)**2+2*np.sum(ks*qkx*qkg)+np.sum(ks*qkg**2)   # sum k q_k^2 (|x| part exact)
    E=2*q0-2*np.log(R/2)-0.25*S
    # density check on (0,1): psi R sqrt(1-u^2) = 1/pi + sum c_k T_k(u), c_k = -k q_k/(2pi)
    ug=np.cos(np.linspace(0.02,np.pi/2-0.02,400)); ck=-ks*(qkx+qkg)/(2*np.pi)
    dens=1/np.pi+np.array([np.sum(ck*np.cos(ks*t)) for t in np.arccos(ug)])
    if verbose: print('rho=%.6f (%.4f pi): R=%.6f  E_min=%.6f  ->  U = -E_min = %.6f   (density min on (0,1): %.4f, edge residual %.1e)'%(rho,rho/np.pi,R,E,-E,dens.min(),edge(R)))
    return R,E
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='fit':
        # closed-form hunt: E(n), n = rho/pi, on a grid; R formula check; least squares E ~ sum_j c_j ln(n+j) + c
        ns=[0.5,1,1.5,2,2.5,3,4,5,6,8,10,12,16,24,32]
        rows=[]
        for n in ns:
            R,E=solve(n*math.pi,verbose=False); rows.append((n,R,E))
            print('n=%6.2f  R=%.8f  R*n(n+2)/(2(n+1))=%.8f  E=%.8f  E-(1+2ln n)=%.6f'%(n,R,R*n*(n+2)/(2*(n+1)),E,E-(1+2*math.log(n))))
        A=np.array([[math.log(n),math.log(n+1),math.log(n+2),math.log(n+3),1.0] for n,_,_ in rows]); b=np.array([E for _,_,E in rows])
        c,res,rk,sv=np.linalg.lstsq(A,b,rcond=None)
        print('least squares E ~ c0 ln n + c1 ln(n+1) + c2 ln(n+2) + c3 ln(n+3) + c4 :',' '.join('%.6f'%v for v in c),' max|resid| = %.2e'%np.max(np.abs(A@c-b)))
        sys.exit()
    if len(sys.argv)>1 and sys.argv[1]=='shape':
        # fixed-ratio numerator shapes: python equilibrium.py shape rho alpha meff
        rho=float(sys.argv[2][:-2])*math.pi if sys.argv[2].endswith('pi') else float(sys.argv[2])
        ALPHA=float(sys.argv[3]); MEFF=int(sys.argv[4])
        globals()['ALPHA']=ALPHA; globals()['MEFF']=MEFF
        R,E=solve(rho)
        print('      prediction for log det M(xi)/K^2 = -(1-alpha)^2 E = %.4f   (compare closeness + intrinsic, both per K^2)'%(-(1-ALPHA)**2*E)); sys.exit()
    args=sys.argv[1:] or ['1pi','2pi','3pi','4pi','6pi','8pi']
    print('pure |x| check (no poles): E = 3 + 2 log(rho/pi):', end=' ')
    print(', '.join('%.4f'%(3+2*math.log(r)) for r in (1,2)))
    meas={'1pi':'-2.1973 (Catalan N=0 fit; -ln 9 = -2.1972)','2pi':'-3.1192 (zeta2 half N=0 fit; -(9/2) ln 2 = -3.1192)','4pi':'about -4.0 (zeta2 spacing 2: -3.869 at K=40, sweep running)'}
    for a in args:
        rho=float(a[:-2])*math.pi if a.endswith('pi') else float(a)
        R,E=solve(rho)
        if a in meas: print('      measured closeness limit: %s'%meas[a])
