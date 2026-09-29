"""gapgas.py — the log-gas of the Hankel machine in z = v/(4K^2) = (t/K)^2 with a possible gap at the origin.

Kernel -log|z-z'| (the Vandermonde in v), h = K-N particles, external field
  W(z) = rho sqrt(z) + int_{t1}^{t2} log(z+y^2) dy - m int_0^{t1} log(z+y^2) dy,   t1 = a/(1-a), t2 = 1/(1-a), a = N/K
(tail poles carry -1, head poles carry +m: density ∝ E_head^m / D_tail).  Equilibrium measure sigma on [za, zb]:
  Frostman  2 U^sigma(z) + W(z) = l on the support.  Chebyshev inversion with z = c + d u, W(c+du) = sum w_k T_k(u):
  sigma = sum c_k T_k(u) du / sqrt(1-u^2),  c_0 = 1/pi,  c_k = -k w_k/(2 pi);
  soft edge at zb:  1/pi + sum_{k>=1} c_k = 0;  soft edge at za (gap): 1/pi + sum (-1)^k c_k = 0;  hard edge at za = 0: no condition.
  E = -log(d/2) + w_0 - (1/8) sum_{k>=1} k w_k^2.
closeness_inf = -(1-a)^2 E  (log det M(xi)/K^2 up to the normalisation of the numerator, which cancels in the closeness).
usage: python gapgas.py rho alpha m      (rho like 2pi)
"""
import sys, math
import numpy as np
M=1<<16; theta=np.pi*(np.arange(M)+0.5)/M; u=np.cos(theta); KMAX=12000; ks=np.arange(1,KMAX+1)
TW=np.exp(-1j*np.pi*ks/(2*M))     # DCT-II via a 2M-point FFT: w_k = (2/M) Re( e^{-i pi k/(2M)} FFT([f,0])[k] )
def Gab(z,a,b):
    z=np.asarray(z,dtype=float); s=np.sqrt(np.maximum(z,1e-300))
    def F(y):
        return (y*np.log(z+y*y+1e-300)-2*y if y>0 else 0.0)+(2*s*np.arctan(y/s) if y>0 else 0.0)
    return F(b)-F(a)
def make_W(rho,alpha,m):
    t1=alpha/(1-alpha); t2=1.0/(1-alpha)
    def W(z): return rho*np.sqrt(np.maximum(z,0))+Gab(z,t1,t2)-m*Gab(z,0.0,t1)
    return W
def cheb(W,za,zb):
    c=0.5*(za+zb); d=0.5*(zb-za); f=W(c+d*u)
    F=np.fft.fft(np.concatenate([f,np.zeros(M)]))
    w0=f.mean(); wk=(2.0/M)*np.real(TW*F[1:KMAX+1])
    return w0,wk,d
def energy(W,za,zb):
    w0,wk,d=cheb(W,za,zb); E=-math.log(d/2)+w0-np.sum(ks*wk**2)/8
    ck=-ks*wk/(2*np.pi)
    ug=np.linspace(-0.995,0.995,600); ang=np.arccos(ug)
    dens=1/np.pi+np.cos(np.outer(ang,ks[:3000]))@ck[:3000]
    return E,dens.min(),w0,wk
def edge_b(W,za,zb):   # soft-edge residual at zb
    w0,wk,d=cheb(W,za,zb); return np.sum(ks*wk)-2.0
def edge_a(W,za,zb):   # soft-edge residual at za
    w0,wk,d=cheb(W,za,zb); return np.sum(((-1.0)**ks)*ks*wk)+2.0*0  # sum (-1)^k k w_k must equal -2*... see below
def solve(rho,alpha,m,verbose=True):
    W=make_W(rho,alpha,m)
    # hard edge at 0 first: find zb with sum k w_k = 2
    lo,hi=1e-4,20.0; flo=edge_b(W,0.0,lo)
    for _ in range(100):
        mid=0.5*(lo+hi); fm=edge_b(W,0.0,mid)
        if (fm>0)==(flo>0): lo,flo=mid,fm
        else: hi=mid
    zb=0.5*(lo+hi); E,dmin,w0,wk=energy(W,0.0,zb)
    beta=rho-m*np.pi          # W(z) = W(0) + beta sqrt(z) + ... : beta>0 attractive cusp (hard edge), beta<0 repulsive (gap)
    kind='hard edge at 0 (beta=%.3f)'%beta
    if beta<-1e-9 or dmin<-1e-3:
        # gap: soft-soft, unknowns (za,zb): conditions  sum_k k w_k = 2 (edge zb)  and  sum_k (-1)^k k w_k = 2 (edge za)
        def F(p):
            za,zb=p; w0,wk,d=cheb(W,za,zb); return np.array([np.sum(ks*wk)-2.0, np.sum(((-1.0)**ks)*ks*wk)-2.0])
        t1=alpha/(1-alpha); p=np.array([0.5*t1*t1 if t1>0 else 0.01,zb])
        for it in range(60):
            f=F(p);
            if np.max(np.abs(f))<1e-10: break
            J=np.zeros((2,2)); hstep=1e-6
            for j in range(2):
                dp=p.copy(); dp[j]+=hstep; J[:,j]=(F(dp)-f)/hstep
            step=np.linalg.solve(J,-f)
            # damping to keep 0<za<zb
            lam=1.0
            while True:
                q=p+lam*step
                if q[0]>1e-9 and q[1]>q[0]+1e-6: break
                lam*=0.5
                if lam<1e-6: break
            p=p+lam*step
        za,zb=p; E,dmin,w0,wk=energy(W,za,zb); kind='gap (beta=%.3f): support [%.6f, %.6f], Frostman residuals %.1e %.1e'%(beta,za,zb,*F(p))
    else: za=0.0
    if verbose:
        print('rho=%.5f alpha=%.3f m=%d : %s  E=%.6f  density min %.4f  -> closeness_inf = -(1-alpha)^2 E = %.4f'%(rho,alpha,m,kind,E,dmin,-(1-alpha)**2*E))
    return E,za,zb
if __name__=='__main__':
    if len(sys.argv)<4:
        print('checks against the closed form ((n+3)/2)ln(n+2) - ((n-1)/2)ln n:')
        for n in (1,2,4):
            E,za,zb=solve(n*math.pi,0.0,0,verbose=False)
            print('  n=%d: E=%.6f  closed form %.6f  zb=%.6f (=R^2=%.6f)'%(n,E,(n+3)/2*math.log(n+2)-(n-1)/2*math.log(n),zb,(2*(n+1)/(n*(n+2)))**2))
        sys.exit()
    rho=float(sys.argv[1][:-2])*math.pi if sys.argv[1].endswith('pi') else float(sys.argv[1])
    solve(rho,float(sys.argv[2]),int(sys.argv[3]))
