"""carried_family.py — reconstruct the carried family F(M) of the Catalan two-state law as exact integers, M = 0..MMAX.
F(odd M) = f(M) = U(Mp+1)/U(1) mod p   (two-digit identity),   F(even M) = h(M) = U(p^2 + Mp + 1)/(f(1) U(1)) mod p   (three-digit identity),
over all primes p > M with the needed N <= NMAX; CRT over the primes (integers: symmetric residue once the modulus exceeds 2|F(M)|
with margin; the reconstruction is confirmed when adding primes no longer changes it).  Writes carried_family.json.
usage: python carried_family.py NMAX MMAX
"""
import sys, json
from flint import fmpq
from sympy import primerange
from sympy.ntheory.modular import crt
import catalan_twostate as CT

def main(NMAX, MMAX):
    U, V = CT.seq(NMAX)
    mp_ = CT.mp_
    out = {}
    for M in range(0, MMAX + 1):
        res = []
        for p in primerange(max(M + 1, 5), NMAX):
            u1 = mp_(U[1], p)
            if u1 == 0: continue
            if M % 2 == 1:
                N = M*p + 1
                if N > NMAX: break
                res.append((p, mp_(U[N], p)*pow(u1, -1, p) % p))
            else:
                N = p*p + M*p + 1
                if N > NMAX: break
                f1 = mp_(U[p + 1], p)*pow(u1, -1, p) % p
                den = f1*u1 % p
                if den == 0: continue
                res.append((p, mp_(U[N], p)*pow(den, -1, p) % p))
        # CRT with symmetric representative; stable if the last 3 primes do not change it
        vals = []
        for k in range(len(res) - 3, len(res) + 1):
            if k < 2: continue
            mods = [p for p, _ in res[:k]]; rs = [r for _, r in res[:k]]
            x, Mod = crt(mods, rs); x, Mod = int(x), int(Mod)
            if x > Mod//2: x -= Mod
            vals.append((x, Mod))
        stable = len(vals) >= 2 and len(set(v for v, _ in vals)) == 1 and vals[-1][1] > 4*abs(vals[-1][0]) + 1
        out[M] = dict(value=str(vals[-1][0]) if vals else None, primes=len(res), stable=bool(stable))
        print('F(%d) = %s  [%d primes, stable %s]' % (M, vals[-1][0] if vals else None, len(res), stable), flush=True)
    json.dump(out, open(r"C:\Users\PC\Desktop\catalan\carried_family.json", 'w'), indent=1)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
