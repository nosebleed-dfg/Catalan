r"""rings.py (2026-10-02): concentric rings of equal area, starting from the unit disc.

The disc of radius r has area pi r^2, so equal areas pi give radii sqrt(1), sqrt(2), sqrt(3), ...: ring n lies between
sqrt(n-1) and sqrt(n).  Its width is sqrt(n) - sqrt(n-1) = 1/(sqrt(n) + sqrt(n-1)).
1. Radii and widths; 2k+1 rings between radius k and k+1; the sum of the approximate widths 1/(2 sqrt n).
2. Lattice points: x^2 + y^2 is an integer, so every lattice point lies exactly on a ring boundary, and boundary n carries
   r_2(n) = 4 (d_1(n) - d_3(n)) of them.  The average is pi, the area of a ring (Gauss circle problem for the error).
3. Which boundaries are occupied: the Landau-Ramanujan constant K, and K^4 = pi^2/(32 G) * prod_{p = 3 mod 4} (1 - p^-4)^-1.
4. Weights: sum r_2(n)/n^2 = 4 zeta(2) G over all rings; over the integer radii, sum r_2(k^2)/k^2 = 8 G.
5. Squaring the plane, z -> z^2, sends boundary n to the circle of radius n and lattice points to Pythagorean points.
6. Prime-numbered rings: 8 points if p = 1 mod 4, none if p = 3 mod 4.
Usage: python rings.py"""
import math
import numpy as np
import mpmath as mp

G = float(mp.catalan); PI = math.pi

print("1. radii and widths")
print("   ring n:   " + "  ".join("%6d" % n for n in range(1, 10)))
print("   radius:   " + "  ".join("%6.4f" % math.sqrt(n) for n in range(1, 10)))
print("   width:    " + "  ".join("%6.4f" % (math.sqrt(n) - math.sqrt(n - 1)) for n in range(1, 10)))
print("   width * (inner + outer radius) = 1 for n <= 10^5: %s"
      % all(abs((math.sqrt(n) - math.sqrt(n - 1))*(math.sqrt(n) + math.sqrt(n - 1)) - 1) < 1e-9 for n in range(1, 10**5)))
print("   rings between radius k and k+1: %s (the odd numbers)" % [(k + 1)**2 - k**2 for k in range(0, 8)])
for N in (10**2, 10**4, 10**6):
    s = float(np.sum(0.5/np.sqrt(np.arange(1, N + 1, dtype=np.float64))))
    print("   sum of 1/(2 sqrt n) for n <= %d, minus sqrt N: %.6f" % (N, s - math.sqrt(N)))
print("   zeta(1/2)/2 = %.6f" % float(mp.zeta(0.5)/2))

print("\n2. lattice points sit on the boundaries")
NMAX = 10**7
q = np.zeros(NMAX + 1, dtype=np.int32)                  # q[n] = d_1(n) - d_3(n) = r_2(n)/4
for d in range(1, NMAX + 1, 2):
    q[d::d] += (1 if d % 4 == 1 else -1)
r2 = 4*q.astype(np.int64)
print("   points on boundary n = 1..25: %s" % r2[1:26].tolist())
cum = np.cumsum(r2)
for k in range(2, 8):
    N = 10**k
    print("   n <= 10^%d: mean points per ring %.6f (pi = %.6f);  points in the disc minus its area: %8.1f;  N^(1/4) = %.1f"
          % (k, cum[N]/N, PI, 1 + cum[N] - PI*N, N**0.25))

print("\n3. occupied boundaries (n a sum of two squares)")
LIM = 10**8
occ = np.zeros(LIM + 1, dtype=bool)
A = int(LIM**0.5)
bs = np.arange(0, A + 1, dtype=np.int64)
for a in range(0, A + 1):
    v = a*a + bs[a:]*bs[a:]
    occ[v[v <= LIM]] = True
occ[0] = False
sieve = np.ones(10**7 + 1, dtype=bool); sieve[:2] = False
for i in range(2, 3163):
    if sieve[i]: sieve[i*i::i] = False
pr = np.nonzero(sieve)[0]
p3 = pr[pr % 4 == 3].astype(np.float64)
K = math.exp(-0.5*np.log1p(-1/p3**2).sum())/math.sqrt(2)
p4 = math.exp(-np.log1p(-1/p3**4).sum())
print("   Landau-Ramanujan constant K = (1/sqrt 2) prod_{p = 3 mod 4} (1 - p^-2)^(-1/2) = %.8f" % K)
print("   K^4 = %.8f;   pi^2/(32 G) * prod_{p = 3 mod 4} (1 - p^-4)^-1 = %.8f;   (1/sqrt 2)(pi^2/(8G))^(1/4) = %.6f"
      % (K**4, PI**2/(32*G)*p4, (PI**2/(8*G))**0.25/math.sqrt(2)))
c = np.cumsum(occ)
for k in range(2, 9):
    x = 10**k
    print("   n <= 10^%d: %9d occupied (share %.4f);   K/sqrt(log x) = %.4f" % (k, c[x], c[x]/x, K/math.sqrt(math.log(x))))
del occ, c

print("\n4. weights")
n = np.arange(1, NMAX + 1, dtype=np.float64)
print("   sum r_2(n)/n^2 over all rings n <= 10^7 (+ tail pi/N) = %.8f;   4 zeta(2) G = %.8f" % (float((r2[1:]/n**2).sum()) + PI/NMAX, 4*PI**2/6*G))
# r_2(k^2)/4 = prod over p = 1 mod 4 of (2 a_p + 1): multiplicative sieve up to KM
KM = 10**6
g = np.ones(KM + 1, dtype=np.int64)
for p in pr[(pr % 4 == 1) & (pr <= KM)]:
    pk = int(p); a = 1
    while pk <= KM:
        g[pk::pk] = g[pk::pk]//(2*a - 1)*(2*a + 1)
        pk *= int(p); a += 1
kk = np.arange(1, KM + 1, dtype=np.float64)
part = np.cumsum(4*g[1:]/kk**2)
for k in (2, 4, 6):
    print("   sum over integer radii k <= 10^%d of r_2(k^2)/k^2 = %.6f" % (k, part[10**k - 1]))
print("   8 G = %.6f;   points on the circles of radius 1..13: %s" % (8*G, (4*g[1:14]).tolist()))

NW = 10**6
wsum = np.zeros(NW + 1, dtype=np.float64)
for d in range(1, NW + 1, 2):
    wsum[d::d] += (1.0 if d % 4 == 1 else -1.0)/d
print("   mean over rings n <= 10^6 of sum_{d|n} chi_-4(d)/d (each divisor's vote weighted 1/d) = %.6f;   G = %.6f" % (wsum[1:].mean(), G))

print("\n5. squaring the plane")
ok = True; cnt = 0
for a in range(-12, 13):
    for b in range(-12, 13):
        z = complex(a, b); w = z*z
        ok &= abs(abs(w) - (a*a + b*b)) < 1e-9 and (int(round(w.real))**2 + int(round(w.imag))**2 == (a*a + b*b)**2)
        cnt += 1
print("   (a + bi)^2 = (a^2 - b^2) + 2ab i lies on the circle of radius a^2 + b^2, for %d lattice points: %s" % (cnt, ok))
print("   e.g. (2 + i)^2 = 3 + 4i on radius 5;  (3 + 2i)^2 = 5 + 12i on radius 13;  (4 + i)^2 = 15 + 8i on radius 17")

print("\n6. prime-numbered rings")
print("   p:        " + " ".join("%3d" % p for p in pr[:25]))
print("   points:   " + " ".join("%3d" % r2[p] for p in pr[:25]))
print("   p mod 4:  " + " ".join("%3d" % (p % 4) for p in pr[:25]))
