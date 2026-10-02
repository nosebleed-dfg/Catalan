r"""digit_games.py (2026-10-02): Wieferich primes, palindromes, reversal, cross multiplication and "four slots", tested exactly.

1. Wieferich primes (2^(p-1) = 1 mod p^2): p - 1 as a repdigit, digit reversals, Eisenstein's formula
       (2^(p-1) - 1)/p = (1/2) sum_{k=1}^{p-1} (-1)^(k-1)/k   (mod p),
   and Wieferich primes in other bases (p < 2*10^6).
2. Four slots: the cusp widths of the six index-12 families; products and cross products.
3. Zagier's case E (limit G/2): u_n integers, v_n rationals, (n+1)^2 x_{n+1} = (12n^2+12n+4) x_n - 32 n^2 x_{n-1}.
   a. Cross multiplication: u_n v_{n+1} - u_{n+1} v_n = 32^n/(n+1)^2 exactly.
   b. Lucas: u_n = prod u_{digit} (mod p) over the base-p digits of n, so u_n mod p is the same for every arrangement of
      the digits (reversal, palindromes, all permutations).
   c. Digit law for the rational part: p^(2L) v_n = v_{top digit} * prod u_{lower digits} (mod p), L = number of lower digits.
      Over the arrangements of the same digits only the leading digit matters, and two arrangements differ by a cross
      product u_a v_b - u_b v_a.
4. Wieferich-type primes for case E: the p-adic depth of p^2 v_p - v_1 and of u_p - u_1, against E_{p-3} mod p.
Usage: python digit_games.py"""
from fractions import Fraction as Fr
from itertools import permutations
from math import gcd
import sympy as sp

def digits(n, b):
    d = []
    while n: d.append(n % b); n //= b
    return d or [0]                       # least significant first
def undigits(d, b):
    return sum(x*b**i for i, x in enumerate(d))
def rev(n, b): return undigits(digits(n, b)[::-1], b)
def vp(x, p):
    x = Fr(x)
    if x == 0: return 99
    a, b, v = x.numerator, x.denominator, 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v
def modp(x, p):
    x = Fr(x); return x.numerator*pow(x.denominator, -1, p) % p

print("1. Wieferich primes")
for p in (1093, 3511):
    assert pow(2, p - 1, p*p) == 1
    alt = sum(Fr((-1)**(k - 1), k) for k in range(1, p))
    print("   p = %d: p - 1 = %s (base 2) = %s (base 8) = %s (base 16);  order of 2 mod p = %d"
          % (p, "".join(map(str, digits(p - 1, 2)[::-1])), "".join(map(str, digits(p - 1, 8)[::-1])),
             "".join("0123456789abcdef"[x] for x in digits(p - 1, 16)[::-1]), sp.n_order(2, p)))
    print("      p - 1 = (2^12 - 1) * %s;  alternating harmonic sum 1 - 1/2 + ... - 1/(p-1) = 0 mod p: %s"
          % (Fr(p - 1, 4095), modp(alt, p) == 0))
    for b in (2, 10):
        r = rev(p, b)
        print("      reversed in base %2d: %d  %s" % (b, r, "prime" if sp.isprime(r) else "= " + str(sp.factorint(r))))
ok = all(modp(Fr(pow(2, p - 1) - 1, p), p) == modp(Fr(1, 2)*sum(Fr((-1)**(k - 1), k) for k in range(1, p)), p) for p in sp.primerange(3, 200))
print("   Eisenstein's formula for all odd p < 200: %s" % ok)
print("   Wieferich primes in base b, p < 2*10^6, and whether p - 1 is a palindrome in base b, b^2, b^3, b^4:")
primes = list(sp.primerange(3, 2*10**6))
for b in range(2, 13):
    found = [p for p in primes if b % p and pow(b, p - 1, p*p) == 1]
    notes = []
    for p in found:
        pal = [k for k in (1, 2, 3, 4) if digits(p - 1, b**k) == digits(p - 1, b**k)[::-1]]
        r = rev(p, b)
        notes.append("%d (p-1 palindrome in b^k for k = %s; reversal %d %s)" % (p, pal or "none", r, "prime" if sp.isprime(r) else "composite"))
    print("      b = %2d: %s" % (b, "; ".join(notes) if notes else "none"))

print("\n2. Four slots: cusp widths of the six families")
for name, w in (("Apery zeta(2)", (5, 5, 1, 1)), ("CDT chi_-3", (6, 3, 2, 1)), ("G, 2tau frame", (4, 4, 2, 2)), ("G, Gamma_0(8)", (8, 2, 1, 1)),
                ("level 9, 3tau frame", (3, 3, 3, 3)), ("level 9, Gamma_0(9)", (9, 1, 1, 1))):
    prod = w[0]*w[1]*w[2]*w[3]; root = sp.integer_nthroot(prod, 2)
    crosses = sorted(set((w[i]*w[j], w[k]*w[l]) for i, j, k, l in permutations(range(4)) if w[i]*w[j] == w[k]*w[l] and i < j and k < l and i < k))
    print("   %-20s %s  product %3d = %d^2 (%s)   equal cross products: %s" % (name, w, prod, root[0], root[1], [c[0] for c in crosses] or "none"))

print("\n3. Case E")
N = 2500
u = [Fr(1), Fr(4)]; v = [Fr(0), Fr(1)]
for n in range(1, N):
    c1 = 12*n*n + 12*n + 4; c0 = 32*n*n; c2 = (n + 1)**2
    u.append((c1*u[n] - c0*u[n - 1])/c2); v.append((c1*v[n] - c0*v[n - 1])/c2)
assert all(x.denominator == 1 for x in u)
print("   a. u_n v_{n+1} - u_{n+1} v_n = 32^n/(n+1)^2 for all n < %d: %s" % (N, all(u[n]*v[n + 1] - u[n + 1]*v[n] == Fr(32**n, (n + 1)**2) for n in range(N))))
chi4 = lambda p: 1 if p % 4 == 1 else -1
print("   b. u_{p-1} = chi_-4(p) (mod p) for all odd p < %d: %s" % (N, all((int(u[p - 1]) - chi4(p)) % p == 0 for p in sp.primerange(3, N))))
for p in (3, 5, 7, 11, 13):
    top = min(N, p**4)
    luc = sum(1 for n in range(top) if int(u[n]) % p != (sp.prod([int(u[d]) for d in digits(n, p)]) % p))
    plain = 0; signed = 0; tested = 0
    for n in range(p, top):
        d = digits(n, p); L = len(d) - 1
        lhs = v[n]*p**(2*L); tested += 1
        if vp(lhs, p) < 0: plain += 1; signed += 1; continue
        rhs = v[d[-1]]*sp.prod([int(u[x]) for x in d[:-1]])
        if modp(lhs, p) != modp(rhs, p): plain += 1
        if modp(lhs, p) != modp(rhs*chi4(p)**L, p): signed += 1
    print("   p = %2d, n < %d: Lucas failures for u: %d;  digit law for v without the sign: %d failures of %d;  with chi_-4(p)^L: %d"
          % (p, top, luc, plain, tested, signed))
# four slots: all arrangements of four digits
p = 7; ds = (1, 2, 4, 5)
vals = {}
for perm in set(permutations(ds)):
    n = undigits(list(perm), p)                      # perm[3] is the leading digit
    vals.setdefault(perm[3], set()).add((int(u[n]) % p, modp(v[n]*p**6, p)))
print("   c. p = 7, digits %s in the four slots, all 24 arrangements: (u_n mod p, p^6 v_n mod p) by leading digit: %s"
      % (ds, {k: sorted(x) for k, x in sorted(vals.items())}))
bad = 0; cnt = 0
for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
    for a in range(1, p):
        for b in range(a + 1, p):
            n1, n2 = a*p + b, b*p + a
            if n2 >= N: continue
            cross = u[a]*v[b] - u[b]*v[a]
            series = u[a]*u[b]*sum(Fr(32**k, (k + 1)**2)/(u[k]*u[k + 1]) for k in range(a, b))
            assert cross == series
            cnt += 1
            if modp(p*p*(v[n2] - v[n1]), p) != modp(chi4(p)*cross, p): bad += 1
print("      reversal: p^2 (v_{bp+a} - v_{ap+b}) = chi_-4(p) (u_a v_b - u_b v_a)  (mod p), and")
print("      u_a v_b - u_b v_a = u_a u_b sum_{k=a}^{b-1} 32^k/((k+1)^2 u_k u_{k+1}) exactly: %d failures of %d" % (bad, cnt))

PMAX = 700
print("\n4. Wieferich-type primes for case E (5 <= p < %d)" % PMAX)
from collections import Counter
rows = []
for p in sp.primerange(5, PMAX):
    e = int(sp.euler(p - 3)) % p
    du = vp(u[p] - u[1], p); dv = vp(p*p*v[p] - chi4(p)*v[1], p)
    ru = modp((u[p] - u[1])/p**2, p); rv = modp((p*p*v[p] - chi4(p)*v[1])/p**2, p) if dv >= 2 else None
    rows.append((p, e, du, dv, ru, rv))
print("   depth of u_p - u_1:                 %s" % dict(Counter(r[2] for r in rows)))
print("   depth of p^2 v_p - chi_-4(p) v_1:   %s" % dict(Counter(r[3] for r in rows)))
print("   primes with E_{p-3} = 0 (mod p):    %s" % [r[0] for r in rows if r[1] == 0])
print("   primes with u_p = u_1 (mod p^3):    %s" % [r[0] for r in rows if r[2] >= 3])
print("   primes with p^2 v_p = chi v_1 (mod p^3): %s" % [r[0] for r in rows if r[3] >= 3])
def centred(x, p): return x if x <= p//2 else x - p
cu = Counter(); cv = Counter()
for p, e, du, dv, ru, rv in rows:
    if e:
        cu[centred(ru*pow(e, -1, p) % p, p)] += 1
        if rv is not None: cv[(chi4(p), centred(rv*pow(e, -1, p) % p, p))] += 1
print("   (u_p - u_1)/p^2 divided by E_{p-3} (mod p), most common values: %s" % cu.most_common(3))
print("   (p^2 v_p - chi v_1)/p^2 divided by E_{p-3} (mod p), by class, most common: %s" % cv.most_common(4))
ok_u = all((u[p] - 4 - 8*chi4(p)*p*p*e) % p**3 == 0 for p, e, *_ in rows)
ok_v = all(vp(p*p*v[p] - chi4(p) - p*p*e, p) >= 3 for p, e, *_ in rows)
print("   u_p = 4 + 8 chi_-4(p) p^2 E_{p-3} (mod p^3) for all 5 <= p < %d: %s" % (PMAX, ok_u))
print("   p^2 v_p = chi_-4(p) + p^2 E_{p-3} (mod p^3) for all 5 <= p < %d: %s" % (PMAX, ok_v))
# Lehmer-type quarter sum, to reach large primes: sum_{0<k<p/4} k^-2 against E_{p-3}
qs = Counter()
for p, e, *_ in rows:
    s = sum(pow(k, -2, p) for k in range(1, p//4 + 1)) % p
    if e: qs[(chi4(p), centred(s*pow(e, -1, p) % p, p))] += 1
print("   sum_{0<k<p/4} k^-2 divided by E_{p-3} (mod p), by class: %s" % qs.most_common(4))
for p in (2946901, 2946907):
    s = 0
    for k in range(1, p//4 + 1): s += pow(k, -2, p)
    print("   p = %d: sum_{0<k<p/4} k^-2 = %d (mod p)  %s" % (p, s % p, "-> p divides E_{p-3}" if s % p == 0 else ""))
for p in (149, 241):
    r10 = rev(p, 10); r2 = rev(p, 2)
    print("   p = %d: reversed in base 10: %d (%s); in base 2: %d (%s)" % (p, r10, "prime" if sp.isprime(r10) else sp.factorint(r10), r2, "prime" if sp.isprime(r2) else sp.factorint(r2)))
