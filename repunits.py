r"""repunits.py (2026-10-02): numbers written with ones only, R_k = 11...1 (k ones) = (10^k - 1)/9, as a ruler.

1. Running totals.  R_k is the number of digit strings with fewer than k slots (the empty string included).  So in the list
   of all strings by length, string number n has k slots exactly when R_k <= n < R_(k+1), and it is the string n - R_k.
   The same statement is zero-free ("bijective") base ten, with digits 1..10.
2. Legendre.  Read the digits of n with repunit weights instead of powers: sum d_j R_j = (n - digit sum)/9, and in a prime
   base p this is the exponent of p in n!.  lcm(1..n) uses the other ruler: its exponent is (number of digits) - 1.
3. The sieve written in ones.  1/9 + 1/99 + 1/999 + ... = sum d(n)/10^n: the n-th decimal is the number of divisors of n.
   The alternating odd version gives the two-squares counts r_2(n)/4, and theta_3(1/10) = 1.2002000020000002...
4. Squares and powers: R_k^2 is the palindrome 123...k...321 for k <= 9; 11^n is Pascal's row for n <= 4; R_k = 3 mod 4
   for k >= 2, so it is never a square nor a sum of two squares.
5. Infinitely many ones: R_k(b) = -1/(b - 1) mod b^k.
6. Primes: p divides R_k iff the period of 1/p divides k (p not 2, 3, 5).  Factors of R_k up to k = 18.
7. The no-repeat totals: sum_{j<=k} 10!/(10-j)!, ending at floor(e 10!).
Usage: python repunits.py"""
from fractions import Fraction as Fr
from itertools import product
from math import comb, factorial, gcd, e
import sympy as sp

R = lambda k, b=10: (b**k - 1)//(b - 1)
def digits(n, b):
    d = []
    while n: d.append(n % b); n //= b
    return d                                     # least significant first

print("1. running totals")
strings = [""]
for L in range(1, 5): strings += ["".join(t) for t in product("0123456789", repeat=L)]
ok = True
for n, s in enumerate(strings):
    k = max(j for j in range(0, 6) if R(j) <= n)
    ok &= (len(s) == k) and (s == (str(n - R(k)).zfill(k) if k else ""))
print("   all %d strings with at most 4 slots, listed by length: string number n has k slots iff R_k <= n < R_(k+1), and it is n - R_k: %s" % (len(strings), ok))
def bij(n, b=10):
    d = []
    while n > 0:
        r = n % b or b; d.append(r); n = (n - r)//b
    return d[::-1]
okb = all(len(bij(n)) == max(j for j in range(0, 7) if R(j) <= n) for n in range(1, 200000))
print("   zero-free base ten (digits 1..10): n has k digits iff R_k <= n < R_(k+1), for n < 200000: %s" % okb)
print("   2345 - 1111 = %d;  1000 in zero-free digits: %s;  10^18 - R_18 = %d" % (2345 - 1111, bij(1000), 10**18 - R(18)))

print("\n2. Legendre: digits weighed by repunits")
def vfact(n, p):
    v = 0; q = p
    while q <= n: v += n//q; q *= p
    return v
okL = all(vfact(n, p) == sum(d*R(j, p) for j, d in enumerate(digits(n, p))) for p in (2, 3, 5, 7, 11, 13) for n in range(1, 5001))
okT = all(sum(n//10**j for j in range(1, 8)) == sum(d*R(j) for j, d in enumerate(digits(n, 10))) == (n - sum(digits(n, 10)))//9 for n in range(1, 50001))
print("   v_p(n!) = sum d_j R_j(p) for p <= 13, n <= 5000: %s;   base ten: sum floor(n/10^j) = sum d_j R_j = (n - digit sum)/9: %s" % (okL, okT))
print("   5^4 = 625: v_5(625!) = %d = 1111 in base 5;  2345 = %s in base 5, so 2345! ends in %d zeros" % (vfact(625, 5), digits(2345, 5)[::-1], vfact(2345, 5)))
n = 2345; rep = []
for k in range(4, 0, -1):
    rep.append(n//R(k)); n %= R(k)
print("   2345 = %s in repunits 1111, 111, 11, 1;  and sum floor(21110/10^j) = %d" % (rep, sum(21110//10**j for j in range(1, 6))))
from math import lcm
okl = all(int(sp.multiplicity(p, lcm(*range(1, n + 1)))) == len(digits(n, p)) - 1 for p in (2, 3, 5, 7) for n in range(1, 400))
print("   lcm(1..n): exponent of p = (number of base-p digits of n) - 1, for p <= 7, n < 400: %s" % okl)

print("\n3. the sieve written in ones")
K = 70; PL = 60
S = sum(Fr(1, 10**k - 1) for k in range(1, K + 1))
dig = str(S.numerator*10**PL//S.denominator).zfill(PL)
dn = [sp.divisor_count(n) for n in range(1, PL + 1)]
print("   1/9 + 1/99 + 1/999 + ... = 0.%s..." % dig)
print("   number of divisors d(n):    %s" % "".join(str(x) if x < 10 else "X" for x in dn))
first_carry = next(n for n in range(1, 200) if sp.divisor_count(n) >= 10)
print("   digits agree with d(n) up to n = %d; first d(n) >= 10 at n = %d (carry);  places with digit 2 below %d: %s"
      % (next(i for i in range(PL) if int(dig[i]) != dn[i]), first_carry, first_carry, [i + 1 for i in range(first_carry - 2) if dig[i] == "2"]))
A = sum(Fr((-1)**j, 10**(2*j + 1) - 1) for j in range(0, 40))
adig = str(A.numerator*10**PL//A.denominator).zfill(PL)
c4 = lambda d: (0, 1, 0, -1)[d % 4]
r2q = [sum(c4(d) for d in range(1, n + 1) if n % d == 0) for n in range(1, PL + 1)]
print("   1/9 - 1/999 + 1/99999 - ... = 0.%s..." % adig)
print("   r_2(n)/4:                      %s   (agree: %s)" % ("".join(map(str, r2q)), adig == "".join(map(str, r2q))))
th = 1 + 2*sum(Fr(1, 10**(m*m)) for m in range(1, 9))
print("   theta_3(1/10) = %s...;  theta_3^2 = 1 + 4 * (the line above): %s"
      % (str(th.numerator*10**40//th.denominator)[0] + "." + str(th.numerator*10**40//th.denominator)[1:], abs(th*th - 1 - 4*A) < Fr(1, 10**55)))

print("\n4. squares and powers of ones")
for k in (2, 3, 4, 9, 10):
    s = str(R(k)**2); print("   R_%d^2 = %s  %s" % (k, s, "palindrome" if s == s[::-1] else "not a palindrome (carry)"))
for n_ in (2, 3, 4, 5):
    print("   11^%d = %d;  Pascal row %s" % (n_, 11**n_, [comb(n_, j) for j in range(n_ + 1)]))
print("   R_k mod 4 for k = 2..18: %s" % sorted(set(R(k) % 4 for k in range(2, 19))))

print("\n5. infinitely many ones in base b: R_k(b) = -1/(b-1) mod b^k")
for b in (2, 3, 7, 10, 13):
    print("   b = %2d: (b - 1) R_k + 1 = 0 mod b^k for k <= 30: %s   so ...111 = -1/%d" % (b, all(((b - 1)*R(k, b) + 1) % b**k == 0 for k in range(1, 31)), b - 1))

print("\n6. primes")
seen = set()
for k in range(2, 19):
    f = sp.factorint(R(k)); new = sorted(p for p in f if p not in seen); seen |= set(f)
    n3 = sum(e_ for p, e_ in f.items() if p % 4 == 3)
    print("   R_%-2d = %-45s new primes: %s   (prime factors 3 mod 4, with multiplicity: %d)" % (k, " * ".join("%d%s" % (p, "^%d" % e_ if e_ > 1 else "") for p, e_ in sorted(f.items())), new, n3))
okp = all(next(k for k in range(1, p) if R(k) % p == 0) == sp.n_order(10, p) for p in sp.primerange(7, 1000))
print("   the first repunit divisible by p has as many ones as the period of 1/p, for 7 <= p < 1000: %s" % okp)
print("   487 (Wieferich in base ten): first repunit divisible by 487 is R_%d, and 487^2 divides it: %s" % (sp.n_order(10, 487), R(486) % 487**2 == 0))
print("   repunit primes R_k, k <= 60: k = %s" % [k for k in range(2, 61) if sp.isprime(R(k))])

print("\n7. totals without repeated digits")
tot = [sum(factorial(10)//factorial(10 - j) for j in range(0, k + 1)) for k in range(0, 11)]
print("   strings with at most k different-digit slots: %s" % tot)
print("   the last one is floor(e * 10!) = %d;  with at most 4 slots: %d of %d = %.4f" % (int(e*factorial(10)), tot[4], R(5), tot[4]/R(5)))
