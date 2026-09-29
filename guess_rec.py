"""guess_rec.py — find the linear recurrence with polynomial coefficients satisfied by an integer (or rational) sequence:
    sum_{k=0}^{r} c_k(n) a_{n+k} = 0,   deg c_k <= d.
Search order r = 1..R and degree d = 0..D: the first (r, d) with a one-dimensional kernel (mod a 62-bit prime) is solved exactly
(kernel mod several primes, normalised, CRT + rational reconstruction, cleared to primitive integers) and verified on every term.
Module use:  guess(seq, R=6, D=40)  ->  (r, d, [c_0(n), ..., c_r(n)] as integer coefficient lists low->high) or None
"""
import sys, math, random
from sympy import nextprime, factorint, Poly, symbols, factor
from sympy.ntheory.modular import crt

def _rows(seq, r, d, P):
    rows = []
    for n in range(0, len(seq) - r):
        row = []
        for k in range(r + 1):
            a = seq[n + k] % P
            pw = 1
            for e in range(d + 1):
                row.append(a*pw % P); pw = pw*n % P
        rows.append(row)
    return rows

def _kernel(rows, ncol, P):
    M = [r[:] for r in rows]
    piv_cols = []; rk = 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], -1, P)
        M[rk] = [x*inv % P for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(x - f*y) % P for x, y in zip(M[i], M[rk])]
        piv_cols.append(c); rk += 1
        if rk == len(M): break
    free = [c for c in range(ncol) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0]*ncol; v[fc] = 1
        for i, pc in enumerate(piv_cols): v[pc] = (-M[i][fc]) % P
        basis.append(v)
    return basis

def _ratrec(x, M):
    B = math.isqrt(M//2); r0, r1 = M, x % M; s0, s1 = 0, 1
    while r1 > B:
        q = r0//r1; r0, r1 = r1, r0 - q*r1; s0, s1 = s1, s0 - q*s1
    if s1 == 0 or abs(s1) > B: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return r1, s1

def guess(seq, R=6, D=40, verbose=True):
    seq = [int(x) for x in seq]
    P0 = (1 << 61) - 1
    for r in range(1, R + 1):
        for d in range(0, D + 1):
            ncol = (r + 1)*(d + 1)
            if len(seq) - r < ncol + 8: break
            ker = _kernel(_rows(seq, r, d, P0), ncol, P0)
            if len(ker) == 0: continue
            if len(ker) > 1:
                if verbose: print('  order %d degree %d: kernel dimension %d (not minimal in degree; try lower order)' % (r, d, len(ker)))
                break
            # exact solution: kernel mod many primes, normalised at the last nonzero coordinate of the mod-P0 kernel
            v0 = ker[0]; norm_idx = max(i for i, x in enumerate(v0) if x)
            vals, mods = [], []
            P = 1 << 62
            modulus = 1; sol = None
            for _ in range(200):
                P = nextprime(P + random.randint(1, 1 << 20))
                kk = _kernel(_rows(seq, r, d, P), ncol, P)
                if len(kk) != 1 or kk[0][norm_idx] == 0: continue
                inv = pow(kk[0][norm_idx], -1, P)
                vals.append([x*inv % P for x in kk[0]]); mods.append(P); modulus *= P
                if len(mods) % 4 == 0:
                    cand = []
                    ok = True
                    for j in range(ncol):
                        x, _ = crt(mods, [v[j] for v in vals])
                        rr = _ratrec(int(x), modulus)
                        if rr is None: ok = False; break
                        cand.append(rr)
                    if not ok: continue
                    L = 1
                    for a, b in cand: L = L*b//math.gcd(L, b)
                    ints = [a*(L//b) for a, b in cand]
                    g = 0
                    for x in ints: g = math.gcd(g, x)
                    ints = [x//g for x in ints]
                    # exact verification on every term
                    good = True
                    for n in range(0, len(seq) - r):
                        s = 0
                        for k in range(r + 1):
                            ck = sum(ints[k*(d + 1) + e]*n**e for e in range(d + 1))
                            s += ck*seq[n + k]
                        if s: good = False; break
                    if good:
                        sol = [ints[k*(d + 1):(k + 1)*(d + 1)] for k in range(r + 1)]
                        break
            if sol is None:
                if verbose: print('  order %d degree %d: kernel found mod p but exact lift failed' % (r, d))
                return None
            if verbose: print('  found: order %d, degree %d, verified on all %d terms' % (r, d, len(seq)))
            return r, d, sol
    return None

def show(sol, var='n'):
    n = symbols(var)
    r, d, cs = sol
    out = []
    for k, c in enumerate(cs):
        poly = sum(ci*n**e for e, ci in enumerate(c))
        out.append((k, factor(poly)))
    return out
