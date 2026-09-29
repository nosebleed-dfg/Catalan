"""profiles_analyze.py — read profile_*.json (from profiles.py) and report
  (1) closeness/height/net vs K per machine and fixed (N, m), with 1/K and 1/K^2 extrapolations,
  (2) the decomposition height = intrinsic + den(small/mid/big) - num,
  (3) exponent profiles: exact e_p for p > K/2 against the alternating-ramp carry model
        M_c(p) = c * sum_{k>=1} (-1)^(k-1) max(0, 2K - k p)
      (c = 2 for zeta2, 3/2 for catalan; residuals printed), and the binned collapse e_p/K vs u = p/K over K.
usage: python profiles_analyze.py [machine ...] [--exact]
"""
import sys, os, json, math, glob
here=os.path.dirname(os.path.abspath(__file__))
args=[a for a in sys.argv[1:] if not a.startswith('--')]; EXACT='--exact' in sys.argv
runs={}
for f in glob.glob(os.path.join(here,'profile_*.json')):
    d=json.load(open(f)); runs.setdefault(d['machine'],[]).append(d)
machines=args or sorted(runs)
def extrap(Ks,vals):
    out=[]
    if len(Ks)>=2:
        (K1,v1),(K2,v2)=(Ks[-2],vals[-2]),(Ks[-1],vals[-1])
        out.append(('a+b/K (last 2)',(v2*K2-v1*K1)/(K2-K1)))
    if len(Ks)>=3:
        import numpy as np
        A=np.array([[1,1/K,1/K**2] for K in Ks[-3:]]); b=np.array(vals[-3:])
        out.append(('a+b/K+c/K^2 (last 3)',float(np.linalg.solve(A,b)[0])))
    if len(Ks)>=3:
        import numpy as np
        A=np.array([[1,math.log(K)/K,1/K] for K in Ks[-3:]]); b=np.array(vals[-3:])
        out.append(('a+b logK/K+c/K (last 3)',float(np.linalg.solve(A,b)[0])))
    return out
def model(mch,K,p):
    """exact-exponent model for the N=0 machines (found empirically at K=40,60,80, p>K/2):
       zeta2  : 2*sum_k (-1)^(k-1) (2K-kp)_+  - [p<K]     (poles with odd leading base-p digit count 4 each)
       catalan: min( 3*#{b>p} + 1 , 2(K-1) ) = min(1.5(2K-p)-0.5, 2K-2)   (no mirror cancellation; trivial cap)"""
    if mch=='zeta2':
        s=0; k=1
        while k*p<2*K: s+=(-1)**(k-1)*(2*K-k*p); k+=1
        return 2*s-(1 if p<K else 0)
    if mch=='catalan':
        return min(1.5*(2*K-p)-0.5, 2*K-2)
    if mch=='zeta2int':      # conjectured analogue: 4 per pole with odd leading base-p digit, poles j<=K
        s=0; k=1
        while k*p<K: s+=(-1)**(k-1)*(K-k*p); k+=1
        return 4*s
    return float('nan')
ASYMPT={'zeta2':4*math.log(2),'catalan':8/3,'zeta2int':2*math.log(2)}   # int phi(u) du for the models
for mch in machines:
    rs=runs.get(mch,[])
    groups={}
    for d in rs: groups.setdefault((d['N'],d['m'],d.get('step',2)),[]).append(d)
    for (N,m,step),seq in sorted(groups.items()):
        seq.sort(key=lambda d:d['K'])
        print('\n=== %s  N=%d  m=%d  spacing=%s ==='%(mch,N,m,'%g (rho=%s)'%(step/2,{'zeta2':'%gpi'%step,'catalan':'%gpi'%(step/2)}.get(mch,'?')) if mch in ('zeta2','catalan') else '-'))
        print(' K    h   closeness   height      net   | intrinsic   den    small    mid     big     num  | time')
        for d in seq:
            K=d['K']; den=d['den']
            sm=sum(e*math.log(int(p)) for p,e in den.items() if int(p)<=K/2)/K**2
            md=sum(e*math.log(int(p)) for p,e in den.items() if K/2<int(p)<2*K)/K**2
            bg=sum(e*math.log(int(p)) for p,e in den.items() if int(p)>=2*K)/K**2
            print('%4d %4d   %8.4f  %8.4f  %8.4f  |  %7.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f | %5.0fs'%(K,d['h'],d['closeness'],d['height'],d['net'],d['log_intrinsic']/K**2,d['log_den']/K**2,sm,md,bg,d['log_num']/K**2,d['time']))
        Ks=[d['K'] for d in seq]
        if len(Ks)>=2:
            for key in ('closeness','height','net'):
                for lab,a in extrap(Ks,[d[key] for d in seq]): print('   %-9s -> %-26s %8.4f'%(key,lab,a))
        if len(seq)>=2:
            print('  binned e_p/K vs u=p/K (den), mean over primes in bin:')
            print('   u-bin     '+' '.join('K=%-4d'%d['K'] for d in seq))
            for i in range(15):
                lo,hi=0.5+0.1*i,0.6+0.1*i; row='  [%.1f,%.1f) '%(lo,hi)
                for d in seq:
                    K=d['K']; xs=[e/K for p,e in d['den'].items() if lo<=int(p)/K<hi]
                    row+=' %6s'%('%.3f'%(sum(xs)/len(xs)) if xs else '  -  ')
                print(row)
        for d in seq:
            K=d['K']; pts=[(int(p),e) for p,e in d['den'].items() if int(p)>K/2]
            res=[e-model(mch,K,p) for p,e in pts]
            if pts and mch in ('zeta2','catalan'):
                print('   K=%d: e_p - model(p) for p>K/2: mean %.2f  max|.| %.2f   (mid+big sum e_p log p / K^2 = %.4f)'%(K,sum(res)/len(res),max(abs(r) for r in res),sum(e*math.log(p) for p,e in pts)/K**2))
            if EXACT and pts:
                print('      p:e_p:resid  '+' '.join('%d:%d:%+.1f'%(p,e,r) for (p,e),r in zip(pts,res)))
            if N==0 and mch in ('zeta2','catalan'):
                # whole-content comparison: actual vs model over all primes < 2K (model extended to small p as is)
                import sympy
                allp=[int(p) for p in d['den']]; act=d['log_den']
                mod=sum(max(0,model(mch,K,p))*math.log(p) for p in sympy.primerange(2,2*K))
                devs=sorted(((e-model(mch,K,int(p)),int(p),e) for p,e in d['den'].items()),key=lambda t:-abs(t[0])*math.log(t[1]))
                print('   K=%d (N=0): den actual %.4f  model %.4f  model-asymptote %.4f  | largest deviations (e_p - model, weighted): %s'%(K,act/K**2,mod/K**2,ASYMPT[mch],' '.join('p=%d:%+.0f(e=%d)'%(p,r,e) for r,p,e in devs[:6])))
                print('      small primes actual e_p: '+' '.join('%d:%d'%(int(p),e) for p,e in sorted(d['den'].items(),key=lambda t:int(t[0])) if int(p)<=K/2))
