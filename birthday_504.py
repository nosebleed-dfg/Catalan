r"""birthday_504.py (2026-10-02): what 5040/10000 = 0.504 is, exactly.

10^4 strings of four decimal digits; 10*9*8*7 = 5040 with no repeated digit.  The ratio is the chance that four digits are
all different: the birthday problem with 4 people and 10 days.
1. The counts, 5040 = 7!, and 10! = 6! 7!.
2. The carry: b(b-1)(b-2)(b-3)/b^4 = 1 - 6/b + 11/b^2 - 6/b^3 (Stirling numbers of the first kind).  At b = 10 the 11 does
   not fit in one digit: 1 - 5/10 = 1/2 and (11 - 10)/100 - 6/1000 = 4/1000.
3. Other bases and slot counts; the base at which four slots cross 1/2.
4. Robin's criterion: sigma(n) < e^gamma n log log n for all n > 5040 is equivalent to the Riemann hypothesis; the
   exceptions up to 10^6 are listed.
5. 504 = -2/zeta(-5), and 5040 = 2 lcm(1..10).
Usage: python birthday_504.py"""
from fractions import Fraction as Fr
from math import factorial, lcm
import mpmath as mp
import sympy as sp

print("1. counts: with repeats %d, without %d = 7! (%s);  10! = 6! 7!: %s;  ratio %s = %s"
      % (10**4, 10*9*8*7, 10*9*8*7 == factorial(7), factorial(10) == factorial(6)*factorial(7), Fr(5040, 10000), float(Fr(5040, 10000))))

x = sp.symbols('x')
print("\n2. x(x-1)(x-2)(x-3) = %s" % sp.expand(x*(x - 1)*(x - 2)*(x - 3)))
b = 10
print("   at b = 10: 1 - 6/10 + 11/100 - 6/1000 = %s;   (1 - 5/10) + ((11 - 10)/100 - 6/1000) = %s + %s"
      % (1 - Fr(6, b) + Fr(11, b*b) - Fr(6, b**3), 1 - Fr(5, b), Fr(11 - b, b*b) - Fr(6, b**3)))
print("   in base b the same split is (1 - 5/b) + ((11 - b) b - 6)/b^3; the first part is 1/2 only for b = 10")

print("\n3. chance that k symbols out of b are all different")
def P(b, k):
    r = Fr(1)
    for j in range(k): r *= Fr(b - j, b)
    return r
print("        " + "  ".join("b=%-5d" % b for b in range(7, 14)))
for k in (3, 4, 5):
    print("   k=%d  " % k + "  ".join("%.4f " % float(P(b, k)) for b in range(7, 14)))
mp.mp.dps = 20
bstar = mp.findroot(lambda t: (1 - 1/t)*(1 - 2/t)*(1 - 3/t) - mp.mpf(1)/2, 9.9)
print("   four slots cross 1/2 at b = %s;  23 people in 365 days: %.4f;  e^(-6/10) = %.4f"
      % (mp.nstr(bstar, 10), float(P(365, 23)), float(mp.e**(-mp.mpf(6)/10))))

print("\n4. Robin's inequality sigma(n) < e^gamma n log log n")
NMAX = 10**6
sig = [0]*(NMAX + 1)
for d in range(1, NMAX + 1):
    for m in range(d, NMAX + 1, d): sig[m] += d
eg = float(mp.e**mp.euler)
import math
exc = [n for n in range(3, NMAX + 1) if sig[n] >= eg*n*math.log(math.log(n))]
print("   exceptions for 3 <= n <= 10^6: %s" % exc)
print("   sigma(5040)/5040 = %.5f  against  e^gamma log log 5040 = %.5f" % (sig[5040]/5040, eg*math.log(math.log(5040))))

print("\n5. 504 = -2/zeta(-5): %s;   5040 = 2 lcm(1..10): %s;   504 = 7*8*9: %s"
      % (mp.nstr(-2/mp.zeta(-5), 12), 5040 == 2*lcm(*range(1, 11)), 504 == 7*8*9))
print("   sum of the four radius sums S(5) + S(6) + S(8) + S(9) = 0.5048685776 (radius_transfer.py): not 0.504")
