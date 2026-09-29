"""catalan_twostate.py — the two-state (carry) Lucas law for Zudilin's Apéry-like numbers for Catalan's constant.
U_n = 2^{4n} u_n (u_n G - v_n -> 0, math/0201024).  Doubled lattice N = 2n with base-p digits N_i; carries of the doubling n -> 2n:
c_{-1} = 0, c_i = floor((2 b_i + c_{i-1}) / p) with b_i the digits of n, so N_i = 2 b_i + c_{i-1} - p c_i.
Law:  U_n = prod_i X_{c_{i-1}}(N_i)  (mod p),  X_0 = cal U (recurrence family: U_{M/2} at even M, half-coset solution at odd M),
                                              X_1 = cal F (carried family: h at even M, f at odd M), cal F(0..12) = 1,4,32,304,3136,...
cal F(M) mod p is read from the data: f(M) = U(Mp+1)/U(1), h(M) = U(p^2+Mp+1)/(f(1) g(1)) (two-digit and three-digit identities).
usage: python catalan_twostate.py NMAX     (tests every even N <= NMAX and every prime 5 <= p with p^2 + p^2 <= NMAX)
"""
import sys
from flint import fmpq
from collections import Counter
from sympy import primerange

def P(x): return 20*x*x - 8*x + 1
def Q(x): return 3520*x**6 + 5632*x**5 + 2064*x**4 - 384*x**3 - 156*x**2 + 16*x + 7

def seq(NMAX):
    half = NMAX//2 + 2
    def run(y0, y1, x0):
        ys = [fmpq(y0)]
        if y1 is None:
            x = fmpq(x0); ys.append(Q(x)*ys[0]/((2*x + 1)**2*(2*x + 2)**2*P(x)))
        else: ys.append(fmpq(y1))
        for k in range(1, half - 1):
            x = fmpq(x0) + k
            ys.append((Q(x)*ys[k] + (2*x - 1)**2*(2*x)**2*P(x + 1)*ys[k - 1])/((2*x + 1)**2*(2*x + 2)**2*P(x)))
        return ys
    ui = run(1, fmpq(7, 4), 0); vi = run(0, fmpq(13, 8), 0); uh = run(1, None, fmpq(1, 2))
    U = [(ui[N//2] if N % 2 == 0 else uh[N//2])*fmpq(2)**(2*N) for N in range(NMAX + 1)]
    V = [(vi[N//2]*fmpq(2)**(2*N) if N % 2 == 0 else None) for N in range(NMAX + 1)]
    return U, V

def mp_(x, p):
    return int(x.p) % p*pow(int(x.q) % p, -1, p) % p

def vpq(x, p):
    if x == 0: return 10**9
    a, b, v = int(x.p), int(x.q), 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v

def carries_digits(n, p):
    """digits N_i of 2n and incoming carries c_{i-1} (from doubling n digitwise)."""
    b = []
    m = n
    while m: b.append(m % p); m //= p
    Nd, cin = [], []
    c = 0
    for bi in b:
        t = 2*bi + c
        Nd.append(t % p); cin.append(c); c = t//p
    if c: Nd.append(c); cin.append(c)     # final carry-out becomes a new top digit (value 1), entered with carry 1
    return Nd, cin

def families(U, p):
    X0 = {M: mp_(U[M], p) for M in range(0, p)}
    u1 = mp_(U[1], p)
    f = {M: mp_(U[M*p + 1], p)*pow(u1, -1, p) % p for M in range(1, p, 2)}
    g1 = X0[1]
    h = {M: mp_(U[p*p + M*p + 1], p)*pow(f[1]*g1 % p, -1, p) % p for M in range(0, p, 2)}
    X1 = dict(h); X1.update(f)
    return X0, X1

def run(NMAX):
    U, V = seq(NMAX)
    tot = Counter(); bad = Counter(); ex = []
    for p in primerange(5, NMAX):
        if 2*p*p + 1 > NMAX: break
        if mp_(U[1], p) == 0: continue
        X0, X1 = families(U, p)
        for N in range(p + 1, NMAX + 1, 2):
            n = N//2
            Nd, cin = carries_digits(n, p)
            assert sum(d*p**i for i, d in enumerate(Nd)) == N
            rhs = 1
            for d, c in zip(Nd, cin): rhs = rhs*(X1 if c else X0)[d] % p
            key = (len(Nd), sum(cin))
            tot[key] += 1
            if mp_(U[N], p) != rhs:
                bad[key] += 1
                if len(ex) < 8: ex.append((N, p, Nd, cin))
    print('two-state Lucas law U_n = prod X_{c_(i-1)}(N_i) mod p, p >= 5, even N <= %d:' % NMAX)
    for key in sorted(tot): print('   digits %d, carried digits %d: %d checks, %d failures' % (key[0], key[1], tot[key], bad[key]))
    print('   total: %d checks, %d failures %s' % (sum(tot.values()), sum(bad.values()), ex))
    return U, V

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 6000)
