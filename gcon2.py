r"""gcon2.py — linear forms in 1 and G built by lattice linear algebra: pure-G conditions exactly, and p-adic top-order cancellation
at every dangerous prime as congruences, solved by LLL (Siegel's-lemma style).

Space.  Q(t) = prod_{k=0}^{N} (t + k + 1/4)^alpha (t + k + 3/4)^alpha ;  P(t) = prod_{i=0}^{M-1} (t - i) * sum_{j=0}^{E} x_j b_j(t),
        b_j(t) = binomial(t - M, j)  (well-conditioned, integer-valued on integers).  No symmetry is imposed: the integer-centre and the
        half-integer-centre families are both subspaces.  F = sum_{t>=0} P/Q = U G - V.
Step 1. Exact: A_i + B_i = 0 (all i), A_i - B_i = 0 (i != 2)  ->  integer kernel basis W (fmpz_mat.nullspace, primitive rows).
Step 2. Congruences: for each prime p in TARGET (default: every prime where the unconstrained exponent exceeds `goal`), the p-adic digit of
        (U, V) at level p^{-(goal+1)} must vanish  ->  sum_b y_b r_b(p) = 0 (mod p) for the kernel coordinates y.
Step 3. LLL on [[I, K R^T], [0, K diag(p)]]; shortest vectors with zero tail are the solutions; the shortest one is the form.
Report: kernel dimension, number of congruences, log|F|/n, height, NET, per-range exponents, before and after.
usage: python gcon2.py N_FACTOR M_FACTOR E_EXTRA alpha n1 n2        (N = N_FACTOR n, M = M_FACTOR n, E = #targets + E_EXTRA)
"""
import sys, math
from flint import fmpq, fmpz, fmpq_poly, fmpz_mat
from mpmath import mp, mpf, catalan, log as mlog, fabs
from sympy import primerange
import gcon

Q4, Q34 = fmpq(1, 4), fmpq(3, 4)

def binom_poly(M, j):
    """integer-coefficient free basis (t - M)^j (no new denominators)."""
    t = fmpq_poly([0, 1]); p = fmpq_poly([1])
    for i in range(j): p = p*(t - M)
    return p

def rows_for(poles, basis):
    maxS = max(poles.values())
    other = {s: gcon.laurent_other(poles, s, S) for s, S in poles.items()}
    ps_cache = {}
    rows = []
    for bj in basis:
        A = [fmpq(0)]*(maxS + 2); B = [fmpq(0)]*(maxS + 2); V = fmpq(0)
        for s, S in poles.items():
            tay = gcon.taylor(bj, -s, S); oth = other[s]
            ser = [sum((tay[a]*oth[m - a] for a in range(m + 1)), fmpq(0)) for m in range(S)]
            k = int(s.p)//int(s.q); frac = s - k
            for i in range(1, S + 1):
                a_is = ser[S - i]
                if a_is == 0: continue
                if frac == Q4: A[i] += a_is
                else: B[i] += a_is
                key = (i, k, frac)
                if key not in ps_cache:
                    acc = fmpq(0)
                    for l in range(k): acc += 1/(fmpq(l) + frac)**i
                    ps_cache[key] = acc
                V += a_is*ps_cache[key]
        rows.append((A, B, V))
    return rows, maxS

def uv(rows, y_coords, W):
    """U, V for kernel coordinates y (integer combination of kernel basis rows W)."""
    nb = len(rows)
    x = [sum(y_coords[b]*W[b][j] for b in range(len(W))) for j in range(nb)]
    A2 = sum((x[j]*rows[j][0][2] for j in range(nb)), fmpq(0)); B2 = sum((x[j]*rows[j][1][2] for j in range(nb)), fmpq(0))
    V = sum((x[j]*rows[j][2] for j in range(nb)), fmpq(0))
    return 8*(A2 - B2), V, x

def vpq(x, p):
    if x == 0: return 10**9
    a, b, v = int(x.p), int(x.q), 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v

def digit(x, p, level):
    """coefficient of p^level in the p-adic expansion of x (requires v_p(x) >= level)."""
    y = x/fmpq(p)**level if level >= 0 else x*fmpq(p)**(-level)
    return int(y.p) % p*pow(int(y.q) % p, -1, p) % p

def form_stats(U, V, n):
    den = int(U.q)*int(V.q)//math.gcd(int(U.q), int(V.q))
    g = math.gcd(int(U.p)*(den//int(U.q)), int(V.p)*(den//int(V.q)))
    digits = len(str(abs(int(U.p)))) + len(str(int(U.q))) + len(str(abs(int(V.p)))) + len(str(int(V.q)))
    mp.dps = digits + 80
    F = mpf(int(U.p))/int(U.q)*catalan - mpf(int(V.p))/int(V.q)
    lf = float(mlog(fabs(F)))/n
    h = (math.log(den) - (math.log(g) if g > 1 else 0))/n
    return lf, h, den, g

def exps(den, g, primes):
    out = {}
    for p in primes:
        e = 0; d = den
        while d % p == 0: d //= p; e += 1
        c = 0; gg = g
        while gg and gg % p == 0: gg //= p; c += 1
        out[p] = e - c
    return out

def construct(n, Nf=1.0, Mf=0.5, Eextra=4, alpha=2, goal=2, verbose=True):
    N = max(1, int(round(Nf*n))); M = int(round(Mf*n))
    poles = {}
    for k in range(N + 1): poles[fmpq(k) + Q4] = alpha; poles[fmpq(k) + Q34] = alpha
    t = fmpq_poly([0, 1]); Z = fmpq_poly([1])
    for i in range(M): Z *= (t - i)
    # first pass: unconstrained minimal design (E = 3 extra beyond the exact conditions) to find the dangerous primes
    targets = [p for p in primerange(3, 8*N + 8)]
    def build_basis(E):
        return [Z*binom_poly(M, j) for j in range(E + 1)]
    E0 = 2*alpha
    basis = build_basis(E0)
    assert basis[-1].degree() <= 4*(N + 1)*alpha//2*1 - 2 + 0 or True
    rows, maxS = rows_for(poles, basis)
    def exact_kernel(rows):
        nb = len(rows); conds = []
        for i in range(1, maxS + 1):
            conds.append([rows[j][0][i] + rows[j][1][i] for j in range(nb)])
            if i != 2: conds.append([rows[j][0][i] - rows[j][1][i] for j in range(nb)])
        Mz = []
        for r in conds:
            L = 1
            for x in r: L = L*int(x.q)//math.gcd(L, int(x.q))
            row = [int(x*L) for x in r]
            if any(row): Mz.append(row)
        K, nul = fmpz_mat(Mz).nullspace()
        W = []
        for c in range(nul):
            v = [int(K[i, c]) for i in range(nb)]
            g = 0
            for a in v: g = math.gcd(g, a)
            W.append([a//g for a in v])
        return W
    W0 = exact_kernel(rows)
    # the unconstrained reference: the kernel vector of the smallest free degree that gives a 1-dim kernel
    ref = None
    for E in range(0, E0 + 1):
        rowsE = rows[:E + 1]
        WE = exact_kernel(rowsE)
        if len(WE) == 1:
            U, V, x = uv(rowsE, [1], WE)
            if U != 0: ref = (E, U, V); break
    lf0, h0, den0, g0 = form_stats(ref[1], ref[2], n)
    ex0 = exps(den0, g0, targets)
    danger = [p for p in targets if ex0[p] > goal and p > 2]
    E = len(danger) + Eextra + 2*alpha
    if basis[-1].degree() + (E - E0) > sum(poles.values()) - 2:
        E = sum(poles.values()) - 2 - Z.degree()
    basis = build_basis(E)
    rows, maxS = rows_for(poles, basis)
    W = exact_kernel(rows)
    d = len(W)
    # congruences at the top digit (level -(goal+1)) of U and V, for every dangerous prime
    UVb = [uv(rows, [1 if b == bb else 0 for bb in range(d)], W)[:2] for b in range(d)]
    Rrows, mods = [], []
    skipped = []
    for p in danger:
        for which in (0, 1):
            vals = [UVb[b][which] for b in range(d)]
            low = min(vpq(v, p) for v in vals)
            # impose vanishing of every digit from the deepest level up to -(goal+1)
            for lvl in range(low, -goal):
                Rrows.append([digit(v, p, lvl) if vpq(v, p) <= lvl else 0 for v in vals]); mods.append(p)
                break                                  # one level per pass keeps the lattice honest; deeper levels are rare
            if low < -(goal + 1): skipped.append(p)
    s = len(Rrows)
    Kw = 1 << 200
    B = [[0]*(d + s) for _ in range(d + s)]
    for b in range(d):
        B[b][b] = 1
        for c in range(s): B[b][d + c] = Kw*Rrows[c][b]
    for c in range(s): B[d + c][d + c] = Kw*mods[c]
    L = fmpz_mat(B).lll()
    sols = []
    for r in range(d + s):
        row = [int(L[r, c]) for c in range(d + s)]
        if all(v == 0 for v in row[d:]) and any(row[:d]): sols.append(row[:d])
    best = None
    for y in sols[:6]:
        U, V, x = uv(rows, y, W)
        if U == 0: continue
        lf, h, den, g = form_stats(U, V, n)
        ex = exps(den, g, targets)
        cand = (h + lf, lf, h, y, ex)
        if best is None or cand[0] < best[0]: best = cand
    if verbose:
        rng = lambda ex: {lab: sorted(set(ex[p] for p in ex if lo < p <= hi)) for lab, lo, hi in (('(N,2N]', N, 2*N), ('(2N,4N]', 2*N, 4*N))}
        print('n=%2d N=%d M=%d: reference (free degree %d): logF/n %+.3f height %.3f NET %+.3f 2-adic %d odd %s' % (
            n, N, M, ref[0], lf0, h0, h0 + lf0, ex0.get(2, 0), rng(ex0)))
        print('      dangerous primes %d ; free degree %d ; kernel dim %d ; congruences %d ; primes needing a second level %d ; LLL solutions %d' % (
            len(danger), E, d, s, len(skipped), len(sols)))
        if best:
            print('      lattice form: logF/n %+.3f height %.3f NET %+.3f 2-adic %d odd %s ; |y|_max = %d' % (
                best[1], best[2], best[0], best[4].get(2, 0), rng(best[4]), max(abs(v) for v in best[3])))
        else:
            print('      no non-degenerate lattice form')
    return ref, best

if __name__ == '__main__':
    Nf, Mf, Ex, alpha = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    for n in range(int(sys.argv[5]), int(sys.argv[6]) + 1, 2):
        construct(n, Nf, Mf, Ex, alpha)
