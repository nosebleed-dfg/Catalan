r"""The non-congruence door, part 2 (2026-09-29): the Belyi map of the index-9 group with cusp widths (1, 2, 6), exactly.

Found by factoring through Gamma_0(2): the permutation data (door_groups.py) are imprimitive with blocks of size 3, and the
group is an index-3 subgroup of Gamma_0(2) (widths 1, 2; one elliptic point of order 2).  With y = (eta(2tau)/eta(tau))^24,
j = (1 + 256y)^3 / y  and  j - 1728 = (1 - 512y)^2 (1 + 64y) / y.  The cubic cover y = t / (1 - t/t6)^3 has a simple zero at
t = 0 (width 1), a double zero at t = oo (width 2) and a triple pole at t6 (width 6); requiring y + 1/64 to have a double root
forces t6 = 27/256 (double root beta = -t6/2, simple root b = 4 t6 = 27/64, the elliptic point).  This script checks the
ramification of j = R(t) over 0, 1728 and oo exactly (sympy).
Usage: python door_belyi.py"""
import sympy as sp

t = sp.symbols('t')
t6 = sp.Rational(27, 256)
y = t/(1 - t/t6)**3
R = sp.together((1 + 256*y)**3/y)
num, den = sp.fraction(sp.factor(R))
print("j = R(t):  numerator %s" % sp.factor(num))
print("           denominator %s" % sp.factor(den))
num2, den2 = sp.fraction(sp.factor(sp.together(R - 1728)))
print("j - 1728:  numerator %s" % sp.factor(num2))
print("degree of R: %d" % max(sp.degree(num, t), sp.degree(den, t)))
print("y + 1/64 = %s" % sp.factor(sp.together(y + sp.Rational(1, 64))))
