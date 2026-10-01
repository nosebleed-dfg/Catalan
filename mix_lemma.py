r"""mix_lemma.py (2026-10-01, real date): mixing families in the Calegari-Dimitrov-Tang holonomy bound.

Theorem 2.5.1 of arXiv 2408.15403: for m functions (Q(x)-linearly independent, holonomic) with denominator rates
sigma_1 <= ... <= sigma_m, all meromorphic on the disc after composing with phi,
    m <= I(phi) / (L - tau),   L = log|phi'(0)|,   tau = (1/m^2) sum (2i - 1) sigma_i,   I(phi) >= L.
So different families may be mixed, and tau is a weighted average of their rates.
Lemma (proved in the results log; tested here): with Phi(S) = |S| tau(S) = (1/|S|) sum (2i - 1) sigma_i,
    Phi(S + {sigma}) - Phi(S) >= sigma      for every multiset S and every sigma >= 0.
Consequence: if S0 is a set of functions that exist unconditionally and S1 a set that exists only under a hypothesis, the
bound can be violated by S0 + S1 only if L > the average rate of S1.
Usage: python mix_lemma.py"""
from fractions import Fraction as Fr
import random

def Phi(types):
    s = sorted(types); m = len(s)
    return Fr(sum((2*i + 1)*x for i, x in enumerate(s)), m) if m else Fr(0)
def tau(types): return Phi(types)/len(types)

print("1. The formula against the paper's own examples")
for name, types, expect in (("log: (0, 1)", [0, 1], Fr(3, 4)), ("Thm 2.7.2: (0, 1, 1)", [0, 1, 1], Fr(8, 9)),
                            ("Thm 2.8.4: (0, 1, 3/2, 3/2, 3/2)", [0, 1, Fr(3, 2), Fr(3, 2), Fr(3, 2)], Fr(69, 50)),
                            ("Thm A: (0, 2, 2, 4 x 11)", [0, 2, 2] + [4]*11, Fr(191, 49))):
    print("   %-36s tau = %-7s (paper: %s)  %s" % (name, tau([Fr(x) for x in types]), expect, tau([Fr(x) for x in types]) == expect))

L = 1.2322548520321622                       # log R* for G (radius_transfer.py); no admissible phi is larger
print("\n2. Mixes for G (L <= %.10f)" % L)
h = Fr(3, 2)
for name, types in (("2 pure-G (3/2) + 1 + sqrt(1-4t)", [0, 0, h, h]), ("1, sqrt, one pure-G", [0, 0, h]),
                    ("1, sqrt, 4 log-type (rate 1), one pure-G", [0, 0, 1, 1, 1, 1, h]),
                    ("1, sqrt, 4 log-type, no G (unconditional)", [0, 0, 1, 1, 1, 1]),
                    ("1, sqrt, two pi^2-mixed (4/3)", [0, 0, Fr(4, 3), Fr(4, 3)])):
    t_ = float(tau([Fr(x) for x in types])); m = len(types)
    q = L/(L - t_) if L > t_ else float('inf')
    print("   %-44s m = %d  tau = %.4f  L - tau = %+.4f  bound if I = L: m <= %.2f" % (name, m, t_, L - t_, q))

print("\n3. Lemma: Phi(S + {sigma}) - Phi(S) >= sigma")
random.seed(3); worst = None
for _ in range(200000):
    S = [Fr(random.randint(0, 12), random.choice((1, 2, 3, 4))) for _ in range(random.randint(0, 9))]
    s = Fr(random.randint(0, 12), random.choice((1, 2, 3, 4)))
    d = Phi(S + [s]) - Phi(S) - s
    if worst is None or d < worst: worst = d
print("   200000 random multisets: minimum of Phi(S + {sigma}) - Phi(S) - sigma = %s" % worst)

print("\n4. What a contradiction would need: L > average rate of the G-dependent functions")
for name, rate in (("case E", 2.0), ("best pure-G families (W*+-)", 1.5), ("pi^2-mixed (W_B)", 4/3)):
    print("   %-28s rate %.4f   L - rate = %+.4f" % (name, rate, L - rate))
