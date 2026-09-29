"""zudilin_ledger.py — carry ledger for Zudilin's Apéry-like linear forms in Catalan's constant (math/0201024).

Recurrence (2)-(4): (2n+1)^2 (2n+2)^2 p(n) u_{n+1} - q(n) u_n - (2n-1)^2 (2n)^2 p(n+1) u_{n-1} = 0,
  p(n) = 20n^2 - 8n + 1,  q(n) = 3520n^6 + 5632n^5 + 2064n^4 - 384n^3 - 156n^2 + 16n + 7,
  u_0 = 1, u_1 = 7/4,  v_0 = 0, v_1 = 13/8,   v_n/u_n -> G,   |u_n G - v_n|^{1/n} -> ((sqrt5-1)/2)^5 = e^{-2.406}.
Experimental inclusions (Zudilin): 2^{4n} u_n in Z,  2^{4n} D_{2n-1}^2 v_n in Z.
For each n in the list: exact factorization of den(u_n), den(v_n); the 2-adic exponent vs 4n; the exponent of every odd
prime p <= 2n in den(v_n), tabulated by range (p <= n, n < p <= 2n) and by p mod 4; the totals log den / n against the
decay 2.406; and the decay itself from G at high precision.
usage: python zudilin_ledger.py [n1 n2 ...]
"""
import sys, math; sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from sympy import primerange
from mpmath import mp, mpf, catalan, log as mlog
ns=[int(a) for a in sys.argv[1:]] or [50,100,200,300,400]
N=max(ns)
p=lambda n: 20*n*n-8*n+1
q=lambda n: 3520*n**6+5632*n**5+2064*n**4-384*n**3-156*n**2+16*n+7
u=[F(1),F(7,4)]; v=[F(0),F(13,8)]
for n in range(1,N):
    A=(2*n+1)**2*(2*n+2)**2*p(n); B=q(n); C=(2*n-1)**2*(2*n)**2*p(n+1)
    u.append((B*u[n]+C*u[n-1])/A); v.append((B*v[n]+C*v[n-1])/A)
mp.dps=int(1.2*N)+60; G=+catalan
def vp(x,pp):
    c=0
    while x%pp==0: x//=pp; c+=1
    return c
def fac(d,bound):
    f={}
    for pp in primerange(2,bound):
        if d%pp==0:
            e=vp(d,pp); f[pp]=e; d//=pp**e
    return f,d
print(' n | log|u_n G - v_n|/n | den(u_n): e_2 (4n) rest | den(v_n): e_2 (4n), log(odd part)/n, of which p<=n and n<p<=2n [D_n^2 = 2 per n, D_{2n}^2 ~ 4 per n] | exponents of primes in (n,2n]: count by exponent, split p=1 mod 4 / p=3 mod 4')
for n in ns:
    un,vn=u[n],v[n]
    form=float(mlog(abs(un*G-vn)))/n
    fu,ru=fac(un.denominator,2*n+3); fv,rv=fac(vn.denominator,2*n+3)
    e2u=fu.get(2,0); e2v=fv.get(2,0)
    odd=sum(e*math.log(pp) for pp,e in fv.items() if pp>2)/n
    low=sum(e*math.log(pp) for pp,e in fv.items() if 2<pp<=n)/n
    high=sum(e*math.log(pp) for pp,e in fv.items() if n<pp<=2*n)/n
    hist={}
    for pp,e in fv.items():
        if n<pp<=2*n: hist[(e,pp%4)]=hist.get((e,pp%4),0)+1
    nprimes_hi=sum(1 for pp in primerange(n+1,2*n+1))
    print('%3d | %+.4f | %d (%d) rest=%s | %d (%d), %.4f, %.4f + %.4f | %d primes in (n,2n]: %s'%(n,form,e2u,4*n,ru,e2v,4*n,odd,low,high,nprimes_hi,' '.join('e=%d,p%%4=%d:%d'%(k[0],k[1],c) for k,c in sorted(hist.items()))))
    if rv!=1 or ru!=1: print('    leftover cofactors (primes > 2n+2):',ru,rv)
# per-prime detail for the largest n: exponent of each prime p <= 2n in den(v_n), as e_p vs the Kummer-type guesses
n=N; fv,rv=fac(v[n].denominator,2*n+3)
print('\nden(v_%d): e_p for odd p (p:e), D_{2n-1}^2 would give e=2 for all p<=2n-1, e=2*floor(log_p(2n-1)) with powers:'%n)
print(' '.join('%d:%d'%(pp,fv.get(pp,0)) for pp in primerange(3,2*n+1)))
print('2-adic: e_2 = %d = 4n %+d'%(fv.get(2,0),fv.get(2,0)-4*n))
