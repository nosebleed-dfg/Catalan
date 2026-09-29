"""vwp_forms.py — exact linear forms u G - v from the very-well-poised 6F5(-1) family with demi-integral parameters
(h0 = a, h4 integers; h1, h2, h3 half-integers; 1 + a - hj - hl > 0), by partial fractions (Zudilin's Lemma 1, generalised).

Summand as a rational function of t (the (-1)^t stripped):
  R(t) = [(2t+a)/a] * [1/(a-1)!] * [G(1+a-h4)/G(h4)] * prod_j [G(1+a-hj)/G(hj)]
         * prod_{i=1}^{h4-1}(t+i) * prod_{i=a-h4+1}^{a-1}(t+i) / prod_{j=1..3} prod_{i=0}^{a-2hj}(t+hj+i)
Poles at t = -(k+1/2), multiplicity <= 3.  sum_{t>=0} (-1)^t/(t+k+1/2)^m = 2^m (-1)^k [beta(m) - S_{m,k}],  S_{m,k} = sum_{l<k} (-1)^l/(2l+1)^m.
So  F = sum_{m,k} A_{m,k} 2^m (-1)^k (beta(m) - S_{m,k}) = U1 beta(1) + U2 G + U3 beta(3) - V;  the well-poised symmetry gives U1 = U3 = 0.
usage: python vwp_forms.py a h1 h2 h3 h4      (h's as decimals like 10.5)   or   python vwp_forms.py zudilin n
"""
import sys; sys.set_int_max_str_digits(0)
from sympy import Rational, symbols, gamma, factorial, apart, Poly, together, fraction, nsimplify, S, Add
from fractions import Fraction as F
from mpmath import mp, mpf, catalan, log as mlog
t=symbols('t')
def forms(a,h1,h2,h3,h4):
    a=int(a); h4=int(h4); hs=[Rational(h) for h in (h1,h2,h3)]
    for hj in hs+[Rational(h4)]:
        assert hj>0
    for hj in hs+[Rational(h4)]:
        for hl in hs+[Rational(h4)]:
            assert 1+a-hj-hl>0, 'constraint 1+a-hj-hl>0 violated'
    # Zudilin's normalisation: R^ = (2t+a) * P1 * P2 * Q1 Q2 Q3 with integer-valued P's and Q_j = (a-2hj)!/prod(t+hj+i);
    # for his parameters this is R_n(t+n) exactly, so (U2, V) = (-1)^n (U_n', V_n) = (-1)^n 8 (u_n, v_n).
    num=(2*t+a)
    for i in range(1,h4): num*=(t+i)
    for i in range(a-h4+1,a): num*=(t+i)
    num=num/factorial(h4-1)**2
    den=S(1)
    for hj in hs:
        for i in range(0,int(a-2*hj)+1): den*=(t+hj+i)
        num*=factorial(int(a-2*hj))
    R=num/den
    pf=apart(R,t,full=False)
    # collect coefficients A_{m,k}: terms c/(t+k+1/2)^m
    A={}
    for term in Add.make_args(pf):
        n_,d_=fraction(together(term))
        dp=Poly(d_,t)
        if dp.degree()==0:
            assert nsimplify(term)==0 or Poly(n_,t).degree()<0, 'polynomial part should vanish: %s'%term
            continue
        # d_ = c*(t + r)^m
        roots=dp.all_roots()
        r=roots[0]; m=len(roots)
        assert all(rt==r for rt in roots), 'unexpected denominator %s'%d_
        lead=dp.LC()
        coeff=Rational(n_)/lead if Poly(n_,t).degree()==0 else None
        assert coeff is not None, 'numerator not constant: %s'%term
        k=-r-Rational(1,2)
        assert k.is_integer and k>=0, 'pole not at a negative half-integer: %s'%r
        A[(m,int(k))]=A.get((m,int(k)),0)+coeff
    U={1:0,2:0,3:0}; V=0
    kmax=max(k for m,k in A)
    PS={m:[Rational(0)] for m in (1,2,3)}
    for m in (1,2,3):
        for l in range(kmax+1): PS[m].append(PS[m][-1]+Rational((-1)**l,(2*l+1)**m))
    for (m,k),c in A.items():
        U[m]+=c*2**m*(-1)**k
        V+=c*2**m*(-1)**k*PS[m][k]
    return U[2],V,U[1],U[3]
def fmt(x): return '%s'%x
if __name__=='__main__':
    if sys.argv[1]=='zudilin':
        n=int(sys.argv[2]); a=3*n+1; h=Rational(2*n+1,2); h4=n+1
        U2,V,U1,U3=forms(a,h,h,h,h4)
        print('Zudilin n=%d: U1=%s U3=%s  U2=%s  V=%s   V/U2=%s  (recurrence gives u_n=v_n scaled: u1=7/4, v1=13/8 -> ratio 13/14)'%(n,U1,U3,U2,V,Rational(V)/Rational(U2)))
        mp.dps=40; print('   U2 G - V = %s'%mlog(abs(mpf(U2.p)/U2.q*catalan-mpf(V.p)/V.q)))
        sys.exit()
    a=int(sys.argv[1]); h1,h2,h3=[Rational(x) for x in sys.argv[2:5]]; h4=int(float(sys.argv[5]))
    U2,V,U1,U3=forms(a,h1,h2,h3,h4)
    from sympy import factorint, lcm
    import math
    # the 10-set c and Zudilin's expected denominator 2^{2M} D_{m1} D_{m2}
    a1=1+a-h1-h2; a2=h3; a3=h4; b2=1+a-h1; b3=1+a-h2
    c={'00':b2+b3-a1-a2-a3-1,'11':a1-1,'21':a2-1,'31':a3-1,'12':b2-a1-1,'13':b3-a1-1,'22':b2-a2-1,'23':b3-a2-1,'32':b2-a3-1,'33':b3-a3-1}
    M=c['22']+c['31']; twoc=sorted((2*v for v in c.values()),reverse=True); m1,m2=int(twoc[0]),int(twoc[1])
    D=lambda N: lcm(list(range(1,N+1))) if N>=1 else 1
    expected=2**(2*int(M))*D(m1)*D(m2)
    dU=Rational(U2).q; dV=Rational(V).q
    print('a=%d h=(%s,%s,%s) h4=%d :  U1=%s  U3=%s  | c = %s'%(a,h1,h2,h3,h4,U1,U3,{k:str(v) for k,v in c.items()}))
    print('  den(U2) = %s   [2^{2M}=2^%d]'%(factorint(dU),2*int(M)))
    print('  den(V)  = %s'%factorint(dV))
    fe=factorint(expected); fv=factorint(dV)
    excess={p:e-fe.get(p,0) for p,e in fv.items() if e>fe.get(p,0)}; deficit={p:fe[p]-fv.get(p,0) for p in fe if fe[p]>fv.get(p,0)}
    print('  expected 2^{2M} D_%d D_%d = %s'%(m1,m2,fe))
    print('  den(V) vs expected: EXCESS primes (den has more) %s ; DEFICIT (carry cancellations) %s ; log den(V)/a = %.4f vs expected rate %.4f'%(excess,deficit,math.log(dV)/a,math.log(expected)/a))
    mp.dps=80; form=mpf(Rational(U2).p)/Rational(U2).q*catalan-mpf(Rational(V).p)/Rational(V).q
    print('  log|U2 G - V| = %.4f ;  log den(V) = %.4f ;  integer form log|den*(U2 G - V)| = %.4f   (per a: %.4f, %.4f, %.4f)'%(mlog(abs(form)),math.log(dV),mlog(abs(form))+math.log(dV),mlog(abs(form))/a,math.log(dV)/a,(mlog(abs(form))+math.log(dV))/a))
