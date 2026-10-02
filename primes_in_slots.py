r"""primes_in_slots.py (2026-10-02): primes among digit strings.  pi(100) = 25, pi(1000) = 168, pi(10^4) = 1229.

1. The counts pi(10^k) and three readings: k pi/10^k (the "half measure"), 10^k/pi (which grows by log 10 per digit),
   and the same in other bases, where k pi(b^k)/b^k -> 1/log b.
2. Four slots: primes among the 10^4 strings, among the 5040 strings with four different digits, and by number of
   different digits; the bias from multiples of 3.
3. The exact structure of four slots: 10^4 - 1 = 9 * 11 * 101.
   mod 9: all arrangements of the same digits agree.   mod 11: reversal flips the sign, so abba is a multiple of 11.
   mod 101: 10^2 = -1, so 10 is a square root of -1 and the string d3 d2 d1 d0 is (d0 - d2) + 10 (d1 - d3).
4. Arrangements of the same four digits: how many are prime, set by set.
5. Sisters (a prime whose reversal is another prime).
6. Primes 1 mod 4 against 3 mod 4.
Usage: python primes_in_slots.py"""
import math
from itertools import combinations, permutations
from collections import Counter
import numpy as np

LIM = 10**8
sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(LIM**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
cum = None
def pi(x): return int(np.count_nonzero(sieve[:x + 1]))

print("1. pi(10^k)")
print("    k      pi(10^k)   k*pi/10^k   10^k/pi   step    k*log(10) - 1")
prev = None
for k in range(1, 9):
    c = pi(10**k); r = 10**k/c
    print("   %2d  %12d   %.5f    %8.4f  %s   %8.4f" % (k, c, k*c/10**k, r, ("%.4f" % (r - prev)) if prev else "      ", k*math.log(10) - 1))
    prev = r
print("   limit of k*pi/10^k: 1/log(10) = %.5f;   it would be 1/2 in base e^2 = %.4f;   3 * 168 = %d" % (1/math.log(10), math.e**2, 3*168))
print("   other bases, k*pi(b^k)/b^k at the largest b^k <= 10^8, against 1/log b:")
for b in range(2, 13):
    k = int(math.log(LIM)/math.log(b) + 1e-9)
    print("      b = %2d  k = %2d  %.4f   1/log b = %.4f" % (b, k, k*pi(b**k)/b**k, 1/math.log(b)))

print("\n2. four slots (strings 0000 .. 9999)")
P = sieve[:10000]
digs = lambda n: (n//1000, n//100 % 10, n//10 % 10, n % 10)
by = Counter(); tot = Counter()
for n in range(10000):
    d = len(set(digs(n))); tot[d] += 1
    if P[n]: by[d] += 1
print("   primes in all strings: %d of 10000 (%.4f)" % (sum(by.values()), sum(by.values())/10000))
for d in (4, 3, 2, 1):
    print("   %d different digits: %4d strings, %4d primes, density %.4f" % (d, tot[d], by[d], by[d]/tot[d]))
m3_all = sum(1 for n in range(10000) if sum(digs(n)) % 3 == 0)/10000
m3_dis = sum(1 for n in range(10000) if len(set(digs(n))) == 4 and sum(digs(n)) % 3 == 0)/tot[4]
print("   multiples of 3: %.4f of all strings, %.4f of the strings with four different digits" % (m3_all, m3_dis))
print("   0.504 * 1229 = %.1f;  with the multiples-of-3 correction: %.1f;  actual: %d" % (0.504*1229, 0.504*1229*(1 - m3_dis)/(1 - m3_all), by[4]))

groups = {d: [n for n in range(10000) if len(set(digs(n))) == d] for d in (4, 3, 2, 1)}
print("   multiples of 11: %s" % ", ".join("%d different: %d of %d" % (d, sum(1 for n in v if n % 11 == 0), len(v)) for d, v in groups.items()))
for label, ps in (("3", (3,)), ("3 and 11", (3, 11))):
    f4 = sum(1 for n in groups[4] if all(n % p for p in ps))/len(groups[4]); fa = sum(1 for n in range(10000) if all(n % p for p in ps))/10000
    print("   predicted primes with four different digits from the rule(s) for %s: %.1f" % (label, 1229*0.504*f4/fa))

print("\n3. structure: 10^4 - 1 = 9 * 11 * 101")
print("   palindromes abba divisible by 11: %s;   multiples of 101 are exactly the strings abab: %s"
      % (all(n % 11 == 0 for n in range(10000) if digs(n) == digs(n)[::-1]),
         all((n % 101 == 0) == (digs(n)[0] == digs(n)[2] and digs(n)[1] == digs(n)[3]) for n in range(10000))))
print("   10^2 = %d (mod 101);  n = (d0 - d2) + 10 (d1 - d3) (mod 101) for every string: %s"
      % (100 % 101 - 101, all((n - ((digs(n)[3] - digs(n)[1]) + 10*(digs(n)[2] - digs(n)[0]))) % 101 == 0 for n in range(10000))))
rev = lambda n: int("".join(map(str, digs(n)[::-1])))
print("   reversal: n + rev(n) = 0 (mod 11) and n = rev(n) (mod 9) for every string: %s"
      % all((n + rev(n)) % 11 == 0 and (n - rev(n)) % 9 == 0 for n in range(10000)))
pal_primes = [n for n in range(10000) if P[n] and digs(n) == digs(n)[::-1]]
print("   four-slot palindromes that are prime: %s" % pal_primes)

print("\n4. arrangements of the same four different digits (210 digit sets, 24 arrangements each)")
dist = Counter(); best = []; kill11 = 0
for s in combinations(range(10), 4):
    arr = [1000*a + 100*b + 10*c + d for a, b, c, d in permutations(s)]
    k = sum(1 for n in arr if P[n]); dist[k] += 1; best.append((k, s))
    cls = Counter(n % 11 for n in arr)
    if cls.get(0): kill11 += 1
    assert len(set(n % 9 for n in arr)) == 1 and all(v % 4 == 0 for v in cls.values())
print("   number of prime arrangements -> number of digit sets: %s" % dict(sorted(dist.items())))
print("   sets with digit sum divisible by 3 (no primes possible except 3 itself): %d" % sum(1 for s in combinations(range(10), 4) if sum(s) % 3 == 0))
print("   sets where two pairs of digits have equal sums mod 11 (8 arrangements are multiples of 11): %d" % kill11)
best.sort(reverse=True)
print("   most prime arrangements: %s" % [(k, "".join(map(str, s))) for k, s in best[:5]])
k, s = best[0]
print("   the %d primes from the digits %s: %s" % (k, s, sorted(1000*a + 100*b + 10*c + d for a, b, c, d in permutations(s) if P[1000*a + 100*b + 10*c + d])))

print("\n5. sisters: primes whose four-slot reversal is a different prime")
pairs = sorted((n, rev(n)) for n in range(10000) if P[n] and P[rev(n)] and n < rev(n))
print("   %d pairs (%d primes);  e.g. %s;  1153 <-> 3511: %s" % (len(pairs), 2*len(pairs), pairs[:6] + pairs[-3:], (1153, 3511) in pairs))
print("   expected if reversal were independent of primality, among strings ending and starting in 1, 3, 7, 9: about %.0f pairs"
      % (0.5*1600*(sum(1 for n in range(10000) if P[n] and digs(n)[0] in (1, 3, 7, 9))/1600)**2))

print("\n6. primes 1 mod 4 against 3 mod 4 up to 10^k")
idx = np.nonzero(sieve)[0]
for k in range(2, 9):
    q = idx[idx <= 10**k]
    a = int(np.count_nonzero(q % 4 == 1)); c = int(np.count_nonzero(q % 4 == 3))
    print("   k = %d: %9d  %9d   difference %5d   sqrt(10^k)/log(10^k) = %.0f" % (k, a, c, c - a, 10**(k/2)/(k*math.log(10))))
