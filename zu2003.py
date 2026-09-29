"""zu2003.py — Zudilin's 2003 question (JTNB 15 (2003) 593-626, §7): "What is a trick that makes arithmetic as it is?"

Forms (both are zetaforms.ZetaForm instances after the shift t -> tau + n + 1, h0 = 3n+2, all h_j = n+1):
  B:  F~_n = (1/2) sum_{t>=1} d^2/dt^2 [ (2t+n) ((t-1)...(t-n)(t+n+1)...(t+2n))^3 / (t(t+1)...(t+n))^6 ]  = u~ z(7) + w~ z(5) - v~
      = r = 3, eta = (3; 1^9).  Zu04's exact degree condition fails at the boundary (sum h_j = 9n+9 > 9n+6), but deg R = -5, which is
      all the construction needs; hence LooseForm below.
  A:  F_{k,n} = n!^{k-1} sum_{t>=1} (2t+n)(t-1)...(t-n)(t+n+1)...(t+2n) / (t(t+1)...(t+n))^{k+1},  k odd  = r = 1, eta = (3; 1^{k+2}).
Zudilin: proved  D_n^8 Phi~^-3 F~ in Z-span (B),  D_n^{k+1} Phi~^-1 F_k (A);  observed to n = 1000:  D_n^7 Phi~^-2 (B),  D_n^k (A);
  Phi~_n = prod_{p < n, {n/p} in [2/3, 1)} p.   Here e_p = exponent of p in the common denominator of all coefficients (not the primitive part).
usage:  python zu2003.py control
        python zu2003.py table A|B n1 n2 [k=5] [OUT.json]     per-prime true exponent vs Zudilin proved / observed vs refined predictor
        python zu2003.py anatomy A|B n p [k=5]               leading terms of A_0 at p: per pole, per order, mirror partners, cancellation
"""
import sys, math, json, time
from multiprocessing import Pool
from flint import fmpq, fmpz
import sympy
import zetaforms as zf


class LooseForm(zf.ZetaForm):
    """ZetaForm without Zu04's exact degree condition; requires deg R <= -2 (convergence, no polynomial part, A_r = 0)."""
    def __init__(self, r, eta, n):
        self.r, self.n = r, n
        self.eta0, self.es = eta[0], list(eta[1:])
        self.q = len(self.es)
        assert r % 2 == 1 and self.q % 2 == 1
        self.h0 = self.eta0*n + 2
        self.h = [ej*n + 1 for ej in self.es]
        degR = 1 + 2*sum(self.h[j] - 1 for j in range(r)) - sum(1 + self.h0 - 2*self.h[j] for j in range(r, self.q))
        assert degR <= -2, 'deg R = %d' % degR
        mu = {}
        for j in range(r):
            hj = self.h[j]
            for m in range(1, hj): mu[m] = mu.get(m, 0) + 1
            for m in range(self.h0 - hj + 1, self.h0): mu[m] = mu.get(m, 0) + 1
        self.mu = mu
        s = {}
        for j in range(r, self.q):
            hj = self.h[j]
            for k in range(hj, self.h0 - hj + 1): s[k] = s.get(k, 0) + 1
        self.s = s
        C = fmpq(1)
        for j in range(r): C /= fmpz(math.factorial(self.h[j] - 1))**2
        for j in range(r, self.q): C *= fmpz(math.factorial(self.h0 - 2*self.h[j]))
        self.C = C


def make(kind, n, k=5):
    if kind == 'B': return LooseForm(3, [3] + [1]*9, n)
    if kind == 'A': return LooseForm(1, [3] + [1]*(k + 2), n)
    raise ValueError(kind)


def vq(x, p):
    """p-adic valuation of an fmpq (large for 0)."""
    if x == 0: return 10**9
    return zf.vp_int(int(x.p), p) - zf.vp_int(int(x.q), p)


def phi_tilde(n, p):
    return p < n and (n % p)*3 >= 2*p          # {n/p} in [2/3, 1)


def vD(n, p):
    e = 0; q = p
    while q <= n: e += 1; q *= p
    return e


def bounds(kind, n, p, k=5):
    """(proved, observed) exponents of p per Zudilin."""
    d = vD(n, p); f = 1 if phi_tilde(n, p) else 0
    if kind == 'B': return 8*d - 3*f, 7*d - 2*f
    return (k + 1)*d - f, k*d


def build(kind, n, k=5):
    Z = make(kind, n, k)
    Z.partial_fractions()
    Z.linear_form()
    return Z


def true_exponents(Z):
    coeffs = [Z.A0] + list(Z.A.values())
    D = fmpz(1)
    for c in coeffs: D = fmpz(int(sympy.ilcm(int(D), int(c.q))))
    return sympy.factorint(int(D))


def row(args):
    kind, n, k = args
    t0 = time.time()
    Z = build(kind, n, k)
    fD = true_exponents(Z)
    out = []
    for p in sympy.primerange(2, n + 1):
        te = int(fD.get(p, 0)); pr, ob = bounds(kind, n, p, k); rf = zf.refined_exponent(Z, p)
        e0 = max(0, -vq(Z.A0, p)); em = max([0] + [-vq(a, p) for a in Z.A.values()])
        out.append(dict(p=p, x=n/p, true=te, proved=pr, observed=ob, refined=rf, a0=e0, zeta=em, phi=phi_tilde(n, p)))
    extra = {int(q): int(e) for q, e in fD.items() if q > n}
    return dict(kind=kind, n=n, k=k, rows=out, extra=extra, slots=sorted(Z.A), secs=time.time() - t0)


def cmd_table(kind, n1, n2, k=5, out=None):
    ns = list(range(n1, n2 + 1))
    res = []
    with Pool(min(14, len(ns))) as pool:
        for R in pool.imap(row, [(kind, n, k) for n in ns]):
            res.append(R)
            rows = R['rows']
            big = [r_ for r_ in rows if r_['p']**2 > 3*R['n'] + 2]
            ob_ok = all(r_['true'] <= r_['observed'] for r_ in rows)
            ob_eq = sum(1 for r_ in rows if r_['true'] == r_['observed'])
            rf_ok = all(r_['true'] <= r_['refined'] for r_ in rows)
            rf_eq_big = sum(1 for r_ in big if r_['true'] == r_['refined'])
            below_rf = [(r_['p'], r_['true'], r_['refined']) for r_ in big if r_['true'] < r_['refined']]
            print('%s n=%3d slots %s | true <= observed: %s (equal at %d/%d primes) | true <= refined: %s | single-digit primes: refined exact %d/%d, true below refined at %s | extra %s [%.0fs]' % (
                kind, R['n'], R['slots'], ob_ok, ob_eq, len(rows), rf_ok, rf_eq_big, len(big), below_rf, R['extra'], R['secs']))
            print('      ' + ' '.join('%d:%d/%d/%d/%d%s' % (r_['p'], r_['true'], r_['observed'], r_['proved'], r_['refined'], '*' if r_['phi'] else '') for r_ in rows))
            sys.stdout.flush()
    if out:
        json.dump(res, open(out, 'w'), indent=0)
    return res


def cmd_anatomy(kind, n, p, k=5):
    Z = build(kind, n, k)
    r = Z.r; h1 = Z.h[0]; h0 = Z.h0
    ps_cache = {}
    def ps(J, m):
        key = (J, m)
        if key not in ps_cache:
            v = fmpq(0)
            for l in range(1, J + 1): v += fmpq(1, fmpz(l)**m)
            ps_cache[key] = v
        return ps_cache[key]
    terms = []
    for (kk, i), b in Z.B.items():
        m = i + r - 1
        w = b*fmpz(math.comb(i + r - 2, r - 1))
        t = w*ps(kk - h1, m)
        if t != 0: terms.append((kk, i, m, t, vq(t, p), vq(b, p)))
    vmin = min(t[4] for t in terms)
    vA0 = vq(Z.A0, p)
    print('%s n=%d p=%d (n/p = %.3f, {n/p} = %.3f, Phi~: %s): h0 = %d, h1 = %d, window k - h1 >= p ; v_p(A0) = %d ; min term valuation %d ; proved/observed %s ; refined %d' % (
        kind, n, p, n/p, (n % p)/p, phi_tilde(n, p), h0, h1, vA0, vmin, bounds(kind, n, p, k), zf.refined_exponent(Z, p)))
    lead = [t for t in terms if t[4] == vmin]
    pk = fmpz(p)
    def unit_mod_p(x, v):
        y = x/(fmpq(pk)**v) if v >= 0 else x*(fmpq(pk)**(-v))
        num = int(y.p) % p; den = int(y.q) % p
        return (num*pow(den, -1, p)) % p
    print('  leading terms (valuation %d): %d' % (vmin, len(lead)))
    by_i = {}
    for (kk, i, m, t, v, vb) in sorted(lead):
        u = unit_mod_p(t, vmin)
        by_i.setdefault(i, []).append((kk, u))
        mirror = [kk2 for kk2 in range(h0 - kk - 3*p, h0 - kk + 3*p + 1) if (kk2 + kk - h0) % p == 0 and kk2 != kk and any(L[0] == kk2 and L[1] == i for L in lead)]
        print('    pole k=%d (j = k-h1 = %d = %d*p + %d) order i=%d (m=%d): v_p(B)=%d, unit mod p = %d ; leading partners k* = h0 - k mod p: %s' % (
            kk, kk - h1, (kk - h1)//p, (kk - h1) % p, i, m, vb, u, mirror))
    tot = 0
    for i, L in sorted(by_i.items()):
        s = sum(u for _, u in L) % p; tot = (tot + s) % p
        print('  order i=%d: %d leading poles, sum of units mod p = %d %s' % (i, len(L), s, '(cancels)' if s == 0 else ''))
    print('  all orders: sum of units mod p = %d %s ; v_p(A0) - vmin = %d' % (tot, '(cancels)' if tot == 0 else '', vA0 - vmin))


def cmd_control():
    # A, k = 5: Zudilin's initial values u1 = 18, w1 = 66, v1 = 98; u2 = 938, w2 = 6125/2, v2 = 74463/16 (F = u z(5) + w z(3) - v)
    for n, (u, w, v) in ((1, (fmpq(18), fmpq(66), fmpq(98))), (2, (fmpq(938), fmpq(6125, 2), fmpq(74463, 16)))):
        Z = build('A', n, 5)
        print('A k=5 n=%d: A5 = %s, A3 = %s, A0 = %s ; Zudilin u, w, v = %s, %s, %s ; match: %s ; identity %s' % (
            n, Z.A.get(5), Z.A.get(3), Z.A0, u, w, v, (Z.A.get(5) == u and Z.A.get(3) == w and Z.A0 == v), Z.identity_ok()))
    for n in (1, 2, 3, 4):
        Z = build('B', n)
        fv = Z.form_value(); nv = Z.numeric()
        print('B n=%d: slots %s, identity %s, numeric %s vs form %s, |diff| %s' % (
            n, sorted(Z.A), Z.identity_ok(), zf.mp.nstr(nv, 12), zf.mp.nstr(fv, 12), zf.mp.nstr(abs(nv - fv), 3)))


if __name__ == '__main__':
    cmd = sys.argv[1]
    kv = {a.split('=')[0]: int(a.split('=')[1]) for a in sys.argv[2:] if '=' in a}
    pos = [a for a in sys.argv[2:] if '=' not in a]
    k = kv.get('k', 5)
    if cmd == 'control': cmd_control()
    elif cmd == 'table': cmd_table(pos[0], int(pos[1]), int(pos[2]), k, pos[3] if len(pos) > 3 else None)
    elif cmd == 'anatomy': cmd_anatomy(pos[0], int(pos[1]), int(pos[2]), k)
