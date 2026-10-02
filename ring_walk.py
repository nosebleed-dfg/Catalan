r"""ring_walk.py (2026-10-02): walking around an equal-area ring with straight lines.

Ring n lies between the circles of radius sqrt(n-1) and sqrt(n).  A line that starts on the outer circle and may not enter
the inner disc goes furthest when it is tangent to the inner circle.  Facts checked here:
1. Every such chord has length exactly 2 (half-chord^2 = n - (n-1) = 1), and it turns by theta_n with cos theta_n = 1 - 2/n.
   The number of lines for one full turn is N_n = 2 pi/theta_n = pi/arcsin(1/sqrt n) = pi sqrt n - pi/(6 sqrt n) - ...
2. The first rings: lines needed, and the overshoot past the starting point.
3. Exact closure needs e^(i theta_n) = ((n-2) + 2 i sqrt(n-1))/n to be a root of unity.  It lies in Q(sqrt(1-n)), whose only
   roots of unity are +-1, +-i and the sixth roots, so cos theta_n in {0, +-1/2, +-1}: n = 1, 2, 4 only (and n = 4/3).
   For n = m^2 + 1 the step is (m + i)/(m - i), a Pythagorean rotation.
4. Whole-number radius r (ring r^2): N = pi/arcsin(1/r).  At r = 10^k its integer part reads off the digits of pi.
   Near closures happen at the continued-fraction denominators of pi: 7, 106, 113, 33102.
5. The tail (overshoot as a fraction of one line) over many rings: its distribution and mean.
6. Dictionary: the step has trace 2 cos theta_n = 2 - 4/n, the trace of T L_w^-1 with w = 4/n (lattice_principle.py).
Usage: python ring_walk.py"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 40
print("1-2. the first rings")
print("    n   chord   cos(theta)   theta (deg)   lines for one turn   lines used   overshoot (deg)   overshoot (fraction of a line)")
for n in range(2, 17):
    rin, rout = math.sqrt(n - 1), math.sqrt(n)
    th = 2*math.acos(rin/rout)
    chord = 2*math.sqrt(rout**2 - rin**2)
    N = 2*math.pi/th
    m = round(N) if abs(N - round(N)) < 1e-9 else math.ceil(N)
    over = m*th - 2*math.pi
    print("   %2d   %.4f   %8.5f     %8.4f      %9.5f           %3d          %8.4f          %.4f"
          % (n, chord, math.cos(th), math.degrees(th), N, m, math.degrees(over), over/th))
print("   cos(theta_n) = 1 - 2/n for n <= 10^4: %s" % all(abs(math.cos(2*math.asin(1/math.sqrt(n))) - (1 - 2/n)) < 1e-12 for n in range(2, 10001)))
print("   expansion check at n = 10^4: N = %.8f;  pi sqrt n - pi/(6 sqrt n) = %.8f" % (math.pi/math.asin(0.01), math.pi*100 - math.pi/600))
print("   three lines: outer radius twice the inner one (cos theta = -1/2, n = 4/3); around the unit disc that is the circle of radius 2")

print("\n3. Pythagorean rings n = m^2 + 1: the step is (m + i)/(m - i)")
for m_ in (1, 2, 3, 4):
    n = m_*m_ + 1; z = complex(m_, 1)/complex(m_, -1)
    print("   n = %2d: step = (%d + i)/(%d - i) = %s;  cos, sin = %d/%d, %d/%d" % (n, m_, m_, z, m_*m_ - 1, n, 2*m_, n))

print("\n4. whole-number radius r: N = pi/arcsin(1/r)")
for k in range(1, 9):
    r = mp.mpf(10)**k
    N = mp.pi/mp.asin(1/r)
    print("   r = 10^%d: N = %s   (pi * r = %s)" % (k, mp.nstr(N, 14 + k), mp.nstr(mp.pi*r, 14 + k)))
best = []
rec = 1.0
for r in range(2, 200001):
    N = math.pi/math.asin(1/r)
    d = math.ceil(N) - N
    if d < 1e-9: continue                         # exact closure (r = 2, the hexagon)
    if d < rec:
        rec = d; best.append((r, math.ceil(N), d))
print("   record small overshoots (radius, lines, fraction of a line): %s" % [(r, m, float("%.3g" % d)) for r, m, d in best])
for r, p in ((7, 22), (106, 333), (113, 355), (33102, 103993)):
    N = mp.pi/mp.asin(mp.mpf(1)/r)
    print("   r = %d: N = %s, so %d lines overshoot by %s of a line (arc %s)" % (r, mp.nstr(N, 12), p, mp.nstr(p - N, 4), mp.nstr((p - N)*2*r*mp.asin(mp.mpf(1)/r), 4)))

print("\n5. the tail over many rings")
NMAX = 10**7
n = np.arange(5, NMAX + 1, dtype=np.float64)
N = np.pi/np.arcsin(1/np.sqrt(n))
tail = np.ceil(N) - N
arc = tail*2*np.sqrt(n)*np.arcsin(1/np.sqrt(n))
print("   rings 5..10^7: mean tail = %.5f of a line;  mean overshoot along the circle = %.5f units" % (tail.mean(), arc.mean()))
hist = np.histogram(tail, bins=10, range=(0, 1))[0]/len(tail)
print("   share of rings by tail in tenths: %s" % np.round(hist, 4).tolist())
for lo, hi in ((5, 10**2), (10**2, 10**4), (10**4, 10**6)):
    sel = tail[(n >= lo) & (n < hi)]
    print("   rings %d..%d: mean tail %.4f" % (lo, hi, sel.mean()))

print("\n6. dictionary with the two-parabolic groups: trace 2 - 4/n = 2 - w")
for nn, name in ((4, "hexagon, 6 lines"), (2, "square, 4 lines"), (4/3, "triangle, 3 lines"), (1, "diameter, 2 lines"), (4/5, "none"), (2/3, "none"), (1/2, "none")):
    w = 4/nn; tr = 2 - w
    kind = "rotation of finite order" if tr in (1, 0, -1) else ("half turn (threshold)" if tr == -2 else "no rotation (|trace| > 2)")
    print("   ring n = %-6s  w = 4/n = %-4s  trace %5.2f   %-28s %s" % (("%.4g" % nn), ("%.4g" % w), tr, kind, name))
