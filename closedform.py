"""closedform.py — E(n) for rho = n pi from the exact Chebyshev coefficients of the pole field, high precision, then PSLQ.

With tau(y) = (sqrt(y^2+R^2) - y)/R:  log(y^2 + R^2 u^2) = 2 log(R/(2 tau)) - 4 sum_m (-1)^m tau^{2m} T_{2m}(u)/(2m),
so the Chebyshev coefficients of Q(Ru) = (rho/2) R|u| + (1/2) int_0^1 log(R^2u^2+y^2) dy are
  q_0 = rho R/pi + int_0^1 log(R/(2 tau)) dy,   q_{2m} = -(-1)^m [ (2 rho R/pi)/(4m^2-1) + J_m/m ],  J_m = int_0^1 tau^{2m} dy.
Edge condition sum_{k>=2} k q_k = 2  <=>  rho R/pi + R - sqrt(1+R^2) = 1  <=>  R = 2(n+1)/(n(n+2)).
E = 2 q_0 - 2 log(R/2) - (1/4) sum_k k q_k^2,  sum_k k q_k^2 = (rho R/pi)^2 + (8 rho R/pi) I1 + 2 I2,
  I1 = int_0^1 (1/2)[1 - artanh(tau)(1-tau^2)/tau] dy,   I2 = sum_m J_m^2/m = -int_0^1 int_0^1 log(1 - tau(y)^2 tau(y')^2) dy dy'.
usage: python closedform.py [n ...]
"""
import sys
from mpmath import mp, mpf, sqrt, atanh, log, quad, pi, pslq, nstr
mp.dps=30
def E_of_n(n):
    n=mpf(n); rho=n*pi; R=2*(n+1)/(n*(n+2))
    tau=lambda y: (sqrt(y*y+R*R)-y)/R
    q0=rho*R/pi+quad(lambda y: log(R/(2*tau(y))),[0,1])
    I1=quad(lambda y: (1-atanh(tau(y))*(1-tau(y)**2)/tau(y))/2,[0,1])
    I2=-quad(lambda y: quad(lambda z: log(1-tau(y)**2*tau(z)**2),[0,1]),[0,1])
    S=(rho*R/pi)**2+(8*rho*R/pi)*I1+2*I2
    return 2*q0-2*log(R/2)-S/4, R
ns=[int(a) for a in sys.argv[1:]] or [1,2,3,4,5,6,8]
for n in ns:
    E,R=E_of_n(n)
    print('n=%d  R=%s  E=%s'%(n,nstr(R,12),nstr(E,20)))
    if n==1: print('   ln 9        =',nstr(log(9),20))
    if n==2: print('   (9/2) ln 2  =',nstr(mpf(9)/2*log(2),20))
    primes=[p for p in (2,3,5,7,11,13) if p<=n+2]
    basis=[E,mpf(1)]+[log(p) for p in primes]
    rel=pslq(basis,maxcoeff=2000,maxsteps=10**6)
    print('   PSLQ over {E, 1, ln p for p<=n+2}:',rel,' i.e. E = %s'%('none' if rel is None else ' '.join('%+g ln%d'%(-mpf(c)/rel[0],p) for c,p in zip(rel[2:],primes))+' %+g'%(-mpf(rel[1])/rel[0])))
