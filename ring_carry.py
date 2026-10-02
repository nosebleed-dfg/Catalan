r"""ring_carry.py (2026-10-02): crossing the starting point more than once, in the two ways.

N = N_n = pi/arcsin(1/sqrt n) is the number of lines for one turn of ring n (ring_walk.py); t = ceil(N) - N is the tail.
Way 1, keep going the same direction.  After j turns the walk has used ceil(jN) lines and its tail is ceil(jN) - jN, which
   is j t reduced mod 1: the tails add, and each time they pass one whole line the lap is one line shorter (a carry).
Way 2, turn around.
   With the same lines the walk retraces itself and lands exactly on the start: the tail cancels.
   With the leftover gap as the new step it is Euclid's algorithm.  Rounding down each time and alternating sides gives
   the regular continued fraction N = a0 + 1/(a1 + 1/(a2 + ...)); always overshooting (same side every time) gives the
   "minus" continued fraction N = b0 - 1/(b1 - 1/(b2 - ...)).
Compared here: the misses |q N - p| stage by stage, the number of stages needed, and the digit statistics over many rings
against the universal laws (Gauss-Kuzmin, Khinchin, Levy: (1/k) log q_k -> pi^2/(12 log 2)).
Usage: python ring_carry.py"""
import math
import mpmath as mp

mp.mp.dps = 160
def Nring(n): return mp.pi/mp.asin(1/mp.sqrt(n))

print("1. same direction: laps and carries")
for n in (3, 5, 7):
    N = Nring(n); t = mp.ceil(N) - N
    laps = [int(mp.ceil(j*N)) - int(mp.ceil((j - 1)*N)) for j in range(1, 21)]
    tails = [float(mp.ceil(j*N) - j*N) for j in range(1, 7)]
    short = sum(1 for j in range(1, 20001) if int(mp.ceil(j*N)) - int(mp.ceil((j - 1)*N)) == int(mp.floor(N)))
    print("   ring %d: N = %s, tail t = %s" % (n, mp.nstr(N, 8), mp.nstr(t, 6)))
    print("      lines per lap, first 20 laps: %s" % laps)
    print("      tail after 1..6 crossings: %s   (t, 2t, 3t, ... mod 1)" % ["%.4f" % x for x in tails])
    print("      share of short laps over 20000 laps: %.5f   (t = %.5f)" % (short/20000, float(t)))

def regular_cf(x, k):
    a = [];
    for _ in range(k):
        f = mp.floor(x); a.append(int(f)); r = x - f
        if r < mp.mpf(10)**(-mp.mp.dps + 20): break
        x = 1/r
    return a
def minus_cf(x, k):
    b = []
    for _ in range(k):
        c = mp.ceil(x); b.append(int(c)); r = c - x
        if r < mp.mpf(10)**(-mp.mp.dps + 20): break
        x = 1/r
    return b
def conv_regular(a):
    p0, q0, p1, q1 = 1, 0, a[0], 1; out = [(p1, q1)]
    for x in a[1:]:
        p0, q0, p1, q1 = p1, q1, x*p1 + p0, x*q1 + q0; out.append((p1, q1))
    return out
def conv_minus(b):
    p0, q0, p1, q1 = 1, 0, b[0], 1; out = [(p1, q1)]
    for x in b[1:]:
        p0, q0, p1, q1 = p1, q1, x*p1 - p0, x*q1 - q0; out.append((p1, q1))
    return out

print("\n2. turning around: Euclid on the leftover")
for n in (3, 5, 7):
    N = Nring(n)
    a = regular_cf(N, 9); b = minus_cf(N, 14)
    print("   ring %d" % n)
    print("      alternating sides: digits %s" % a)
    print("         (lines, turns, miss in lines): %s" % [(p, q, float("%.2g" % float(q*N - p))) for p, q in conv_regular(a)[:7]])
    print("      same side every time: digits %s" % b)
    print("         (lines, turns, miss in lines): %s" % [(p, q, float("%.2g" % float(p - q*N))) for p, q in conv_minus(b)[:7]])

print("\n3. stages needed to bring the miss below 10^-6 of a line, rings 3..400 (ring 4 excluded)")
st_r = []; st_m = []
for n in range(3, 401):
    if n == 4: continue
    N = Nring(n)
    cr = conv_regular(regular_cf(N, 40)); k = next(i for i, (p, q) in enumerate(cr) if abs(q*N - p) < 1e-6); st_r.append(k + 1)
    x = N; p0, q0, p1, q1 = 1, 0, int(mp.ceil(N)), 1; k = 1
    r = mp.ceil(x) - x
    while p1 - q1*N >= 1e-6 and k < 200000:
        x = 1/r; c = int(mp.ceil(x)); r = c - x
        p0, q0, p1, q1 = p1, q1, c*p1 - p0, c*q1 - q0; k += 1
    st_m.append(k)
st_r.sort(); st_m.sort()
print("   alternating sides: median %d, mean %.1f, max %d" % (st_r[len(st_r)//2], sum(st_r)/len(st_r), st_r[-1]))
print("   same side:         median %d, mean %.1f, max %d" % (st_m[len(st_m)//2], sum(st_m)/len(st_m), st_m[-1]))

print("\n4. digit statistics of the alternating method over rings 3..2000 (first 60 digits each)")
cnt = {}; tot = 0; logsum = mp.mpf(0); levy = []
for n in range(3, 2001):
    if n == 4: continue
    N = Nring(n)
    a = regular_cf(N, 61)
    for d in a[1:]:
        cnt[d] = cnt.get(d, 0) + 1; tot += 1; logsum += mp.log(d)
    q = conv_regular(a)[-1][1]
    levy.append(float(mp.log(q))/ (len(a) - 1))
gk = lambda k: math.log2(1 + 1/(k*(k + 2)))
print("   share of digit 1, 2, 3, 4: %s   (Gauss-Kuzmin: %s)"
      % (["%.4f" % (cnt.get(k, 0)/tot) for k in (1, 2, 3, 4)], ["%.4f" % gk(k) for k in (1, 2, 3, 4)]))
print("   geometric mean of the digits: %.4f   (Khinchin 2.6854)" % float(mp.exp(logsum/tot)))
print("   mean of (1/k) log q_k: %.4f   (Levy: pi^2/(12 log 2) = %.4f);  misses shrink by e^(2 * that) = %.2f per stage"
      % (sum(levy)/len(levy), math.pi**2/(12*math.log(2)), math.exp(2*math.pi**2/(12*math.log(2)))))
