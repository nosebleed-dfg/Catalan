"""Real-G push, step 3 (2026-09-29; the user: "Might need to attack 3G/6 or 6G/12 more modular arithmetic").
Level-6/12 search: Apery-like binomial sums  u_n(z) = sum_k T(n,k) z^k  for templates T built from 2- and 3-structure,
twists z in {2^a 3^b, sign}. For each family:
  1. U_n = q^n u_n (integers, z = p/q); guess the recurrence (guess_rec.py, order <= 2 here; order 3 reported only);
  2. second solution V (V_0 = 0, V_1 = 1), Apery limit L = lim V_n/U_n;
  3. PSLQ of L against the level-12 basis 1, pi^2, G, L(2,chi_-3), pi*ln2, pi*ln3, ln2^2, ln2*ln3, ln3^2, pi*sqrt3, pi, ln2, ln3;
  4. if G occurs: true tau = ln(den V_n)/n, r = 1/|recessive root|, gate S = ln(16 r) - tau (CDT; need S > 0, realistically > 0.8).
Usage: python level12_search.py [N]"""
import sys, math, cmath, io, contextlib
from fractions import Fraction as F
import mpmath as mp
sys.path.insert(0, r'C:\Users\PC\Desktop\catalan')
from guess_rec import guess
C = math.comb
N = int(sys.argv[1]) if len(sys.argv) > 1 else 360
M = 70                                    # terms used to guess the recurrence
mp.mp.dps = 90
L3 = mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])
BASIS = [('1', mp.mpf(1)), ('pi^2', mp.pi**2), ('G', +mp.catalan), ('L3', L3), ('pi*ln2', mp.pi*mp.log(2)),
         ('pi*ln3', mp.pi*mp.log(3)), ('ln2^2', mp.log(2)**2), ('ln2*ln3', mp.log(2)*mp.log(3)), ('ln3^2', mp.log(3)**2),
         ('pi*sqrt3', mp.pi*mp.sqrt(3)), ('pi', mp.pi), ('ln2', mp.log(2)), ('ln3', mp.log(3))]
TEMPLATES = {
    'C(n,k)^2 C(2k,k)':            lambda n, k: C(n, k)**2*C(2*k, k),
    'C(n,k)^2 C(n+k,k)':           lambda n, k: C(n, k)**2*C(n + k, k),
    'C(n,k) C(2k,k) C(2n-2k,n-k)': lambda n, k: C(n, k)*C(2*k, k)*C(2*n - 2*k, n - k),
    'C(n,k)^3':                    lambda n, k: C(n, k)**3,
    'C(n,k) C(2k,k)^2':            lambda n, k: C(n, k)*C(2*k, k)**2 if 2*k <= 2*n else 0,
    'C(2k,k)^2 C(2n-2k,n-k)^2':    lambda n, k: C(2*k, k)**2*C(2*n - 2*k, n - k)**2,
    'C(n,k) C(2k,k) C(3k,k)':      lambda n, k: C(n, k)*C(2*k, k)*C(3*k, k),
    'C(n,k)^2 C(3k,k)':            lambda n, k: C(n, k)**2*C(3*k, k),
    'C(n,k) C(n+k,k) C(2k,k)':     lambda n, k: C(n, k)*C(n + k, k)*C(2*k, k),
    'C(n,2k) C(2k,k)^2':           lambda n, k: C(n, 2*k)*C(2*k, k)**2,
    'C(n,3k) C(2k,k) C(3k,k)':     lambda n, k: C(n, 3*k)*C(2*k, k)*C(3*k, k),
}
ZS = sorted({s*F(2)**a*F(3)**b for s in (1, -1) for a in range(-4, 5) for b in range(-2, 3) if abs(a) + abs(b) <= 4}, key=lambda x: (abs(x), x))
def seq_U(T, z, upto):
    p, q = z.numerator, z.denominator
    out = []
    for n in range(upto):
        s = sum(T(n, k)*p**k*q**(n - k) for k in range(0, n + 1))   # q^n * sum T z^k
        out.append(s)
    return out
def run_rec(cs, init, upto):
    """solve sum_k c_k(n) a_{n+k} = 0 forward (order r = len(cs)-1) from init (length r)."""
    r = len(cs) - 1
    a = [F(x) for x in init]
    ev = lambda c, n: sum(ci*n**e for e, ci in enumerate(c))
    for n in range(0, upto - r):
        lead = ev(cs[r], n)
        if lead == 0: return None
        a.append(-sum(ev(cs[k], n)*a[n + k] for k in range(r))/lead)
    return a
hits = []
for tname, T in TEMPLATES.items():
    for z in ZS:
        U = seq_U(T, z, M)
        if all(x == 0 for x in U[1:]): continue
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            sol = guess(U, R=2, D=8, verbose=False)
        if sol is None: continue
        r, d, cs = sol
        if r != 2: continue
        # characteristic roots from the leading coefficients
        lc = [c[d] if len(c) > d else 0 for c in cs]
        deg = max(len(c) - 1 for c in cs)
        lc = [(c[deg] if len(c) > deg else 0) for c in cs]
        if lc[2] == 0: continue
        disc = lc[1]**2 - 4*lc[2]*lc[0]
        rt = [(-lc[1] + cmath.sqrt(disc))/(2*lc[2]), (-lc[1] - cmath.sqrt(disc))/(2*lc[2])]
        rt.sort(key=abs, reverse=True)
        if abs(abs(rt[0]) - abs(rt[1])) < 1e-9: continue            # equal moduli: no Apery limit
        Uf = run_rec(cs, [U[0], U[1]], N)
        V = run_rec(cs, [0, 1], N)
        if Uf is None or V is None or Uf[-1] == 0: continue
        Lv = mp.mpf(V[-1].numerator)/V[-1].denominator/(mp.mpf(Uf[-1].numerator)/Uf[-1].denominator)
        Lv2 = mp.mpf(V[-40].numerator)/V[-40].denominator/(mp.mpf(Uf[-40].numerator)/Uf[-40].denominator)
        err = abs(Lv - Lv2)
        digits = int(-mp.log10(err)) if err > 0 else 90
        if digits < 70: continue                                     # need >= 70 correct digits
        # PSLQ at 60 digits, small height only; then verify the relation to the full accuracy with a safety margin
        with mp.workdps(60):
            rel = mp.pslq([Lv] + [b for _, b in BASIS], maxcoeff=2000, maxsteps=10**6)
        if not rel or rel[0] == 0: continue
        resid = abs(sum(k*b for k, b in zip(rel, [Lv] + [b for _, b in BASIS])))
        if resid > mp.mpf(10)**(-(digits - 5)): continue              # spurious at higher precision
        terms = {nm: k for k, (nm, _) in zip(rel[1:], BASIS) if k}
        tau = math.log(V[-1].denominator)/(N - 1) if V[-1].denominator > 1 else 0.0
        rr = 1/abs(rt[1]) if abs(rt[1]) > 1e-12 else float('inf')
        S = math.log(16*rr) - tau if rr != float('inf') else float('inf')
        desc = ' + '.join('%d*%s' % (k, nm) for nm, k in terms.items())
        tag = 'G' if 'G' in terms else '-'
        if 'G' in terms: hits.append((S, tname, str(z), desc, rel[0], tau, rr))
        print('%-30s z=%-6s roots |%.3f| |%.3f|  tau=%.3f r=%.4f S=%+.3f  %s  %d*L = %s' % (
            tname, str(z), abs(rt[0]), abs(rt[1]), tau, rr, S, tag, rel[0], desc))
        sys.stdout.flush()
print('\n=== families whose Apery limit contains G, sorted by the gate S ===')
for h in sorted(hits, reverse=True):
    print('S=%+.3f  %-30s z=%-6s  tau=%.3f r=%.4f   %d*L = %s' % (h[0], h[1], h[2], h[5], h[6], h[4], h[3]))
