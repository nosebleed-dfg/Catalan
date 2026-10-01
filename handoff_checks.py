r"""handoff_checks.py (2026-10-01, real date): independent checks of the chat-side handoff of 2026-09-30
("Catalan session ... widths and 2-adic"), run against this repo's own definitions.

A. Width criterion.  -I lies in <T, L_w> mod 2^a iff 4 does not divide w.  Proof (handoff): mod 2^a, L_w generates the same
   cyclic group as L at the 2-part of w; w odd: T, L_1 generate everything; w = 2 odd: (T^-1 L_2)^2 = -I exactly;
   4 | w: every element is (1 *; 0 1) mod 4.  Here: the exact identity, and brute force for w <= 32, a = 2..6.
   Same shape at 3: -I in <T, L_w> mod 3^a iff 3 does not divide w (w <= 27, a = 1..3).
B. The weight-1 form of the index-9 (1,2,6) family:  E4^(1/4) (A(0)/A(t))^(1/4) = theta4(2 tau)^2 (1 - 256t/27)^(-3/4),
   theta4(2 tau) = sum (-1)^n q^(n^2).  Exact q-series comparison.
C. The integer sequence a_n = 27^n u_n = 256^n [s^n] f  (s = 256t/27):
   (n+1)^2 a_{n+1} = 4(80n^2 + 72n + 21) a_n - 1024 (4n-1)^2 a_{n-1},  a = 1, 84, 12228, ...
   Proved in the results log: a_n is an integer, and a_n = 64^n prod_{i<=n}(4i-1)/n!  (mod 27).  Checked here:
   exact divisibility, v_2(a_n) = 2 s_2(n), v_3(a_n) = number of borrows in (...2020)_3 - n  (Kummer for C(-3/4, n)).
   So v_3(u_n) = -3n + O(log n): the 3-adic cost is 3 ln 3 = 3.296 per n (not 3.2).
D. The handoff's L = lim v_n/w_n (v_0 = 0, v_1 = 1) and Lambda_1 = (1701 L - 81)/1024 against door_recurrence_out.txt.
E. C = 1/2 + sum (-1)^k (k+1) T_{2k+1} (this repo's 2-adic Catalan constant): v_2(E_{2^N - 2} - 2C) = N, N = 3..9,
   i.e. C is the 2-adic limit of beta(-2k) = E_{2k}/2 as 2k -> -2: the Kubota-Leopoldt zeta_2(2) (Beukers 2008, section 1).
Usage: python handoff_checks.py [NMAX]"""
import sys
from collections import deque
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000

# ---------------- A ----------------
print("A. Width criterion")
def mul2(A, B, m=None):
    r = (A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3], A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3])
    return tuple(x % m for x in r) if m else r
Tinv = (1, -1, 0, 1); L2 = (1, 0, 2, 1)
M = mul2(Tinv, L2)
print("   (T^-1 L_2)^2 = %s  (= -I exactly: %s)" % (mul2(M, M), mul2(M, M) == (-1, 0, 0, -1)))
def has_minus_I(w, m):
    I = (1, 0, 0, 1); target = (m - 1, 0, 0, m - 1)
    gens = [(1, 1, 0, 1), (1, 0, w % m, 1)]
    seen = {I}; dq = deque([I])
    while dq:
        x = dq.popleft()
        for g in gens:
            y = mul2(x, g, m)
            if y not in seen:
                if y == target: return True
                seen.add(y); dq.append(y)
    return target in seen
bad = [(w, a) for a in range(2, 7) for w in range(1, 33) if has_minus_I(w, 2**a) != (w % 4 != 0)]
print("   mod 2^a, a = 2..6, w = 1..32: -I in <T, L_w> iff 4 does not divide w; exceptions: %s" % bad)
bad3 = [(w, a) for a in range(1, 4) for w in range(1, 28) if has_minus_I(w, 3**a) != (w % 3 != 0)]
print("   mod 3^a, a = 1..3, w = 1..27: -I in <T, L_w> iff 3 does not divide w; exceptions: %s" % bad3)

# ---------------- B ----------------
P0 = 36
src = open(r"C:\Users\PC\Desktop\catalan\door_family.py", encoding="utf-8").read().split('print("(i) limits')[0]
argv_keep = sys.argv; sys.argv = ["door_family.py", str(P0)]
ns = {}; exec(src, ns); sys.argv = argv_keep
P = ns['P']; build = ns['build']; expand_in = ns['expand_in']; mul = ns['mul']; powser = ns['powser']
t, f, Dt, tp, c = build(Fr)
th = [Fr(0)]*P; k = 0
while k*k < P:
    th[k*k] += (-1)**k*(1 if k == 0 else 2); k += 1
u1 = [(-c)*x for x in t]; u1[0] += 1
rhs = mul(mul(th, th), powser(u1, Fr(-3, 4)))
print("\nB. f = E4^(1/4) A(t)^(-1/4) equals theta4(2tau)^2 (1 - 256t/27)^(-3/4): %s  (%d q-coefficients, exact)"
      % (all(f[i] == rhs[i] for i in range(P)), P))
u = expand_in(f, tp)

# ---------------- C ----------------
print("\nC. a_n = 27^n u_n, n <= %d" % NMAX)
def s2(n): return bin(n).count("1")
def vp(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
def borrows(n):
    b = 0; cnt = 0; i = 0
    while n > 0 or b:
        ai = 0 if i % 2 == 0 else 2
        if ai - n % 3 - b < 0: b = 1; cnt += 1
        else: b = 0
        n //= 3; i += 1
    return cnt
a = [1, 84]
ok_div = True
for n in range(1, NMAX):
    num = 4*(80*n*n + 72*n + 21)*a[n] - 1024*(4*n - 1)**2*a[n - 1]
    if num % (n + 1)**2: ok_div = False
    a.append(num//(n + 1)**2)
print("   first terms: %s" % a[:7])
print("   integrality of the recurrence for all n: %s" % ok_div)
print("   a_n = 27^n u_n for n < %d (exact expansion of f in t): %s" % (P, all(Fr(27)**n*u[n] == a[n] for n in range(P))))
print("   v_2(a_n) = 2 s_2(n): %s" % all(vp(a[n], 2) == 2*s2(n) for n in range(1, NMAX + 1)))
cc = 1; ok27 = True
v3mis = []; low = 0; lowok = 0; maxv3 = 0
for n in range(1, NMAX + 1):
    cc = cc*64*(4*n - 1)
    assert cc % n == 0
    cc //= n
    if (a[n] - cc) % 27: ok27 = False
    v3 = vp(a[n], 3); bo = borrows(n); maxv3 = max(maxv3, v3)
    if bo <= 2:
        low += 1; lowok += (v3 == bo)
    if v3 != bo: v3mis.append((n, v3, bo))
print("   a_n = 64^n prod(4i-1)/n!  (mod 27): %s" % ok27)
print("   v_3(a_n) = borrows(n) when borrows <= 2 (proved): %d of %d" % (lowok, low))
print("   v_3(a_n) = borrows(n) for every n <= %d: %s   (max v_3 = %d; first exceptions %s)" % (NMAX, not v3mis, maxv3, v3mis[:6]))
import math
print("   so v_3(u_n) = -3n + v_3(a_n): cost 3 ln 3 = %.4f per n" % (3*math.log(3)))

# ---------------- D ----------------
mp.mp.dps = 90
w0, w1 = mp.mpf(1), mp.mpf(21); v0, v1 = mp.mpf(0), mp.mpf(1)
for n in range(1, 220):
    w0, w1 = w1, ((80*n*n + 72*n + 21)*w1 - 64*(4*n - 1)**2*w0)/(n + 1)**2
    v0, v1 = v1, ((80*n*n + 72*n + 21)*v1 - 64*(4*n - 1)**2*v0)/(n + 1)**2
L = v1/w1
lam_log = mp.mpf("0.0190568597202821465163513994623608870150249754754331853445047")
print("\nD. L = lim v_n/w_n = %s" % mp.nstr(L, 45))
print("   handoff:          0.05909125476400289125969655088151531352")
print("   (1701 L - 81)/1024 - Lambda_1 (door_recurrence_out.txt, 60 digits) = %s" % mp.nstr((1701*L - 81)/1024 - lam_log, 5))

# ---------------- E ----------------
MOD = 2**200
S = 0
for k in range(0, 104):
    B = sp.bernoulli(2*k + 2)
    T = Fr(2**(2*k + 2)*(2**(2*k + 2) - 1)*abs(int(B.p)), int(B.q)*(2*k + 2))
    assert T.denominator == 1
    S += (-1)**k*(k + 1)*int(T)
twoC = (1 + 2*S) % MOD
res = []
for N in range(3, 10):
    E = int(sp.euler(2**N - 2))
    d = (E - twoC) % MOD
    res.append((N, vp(d, 2) if d else 200))
print("\nE. v_2(E_{2^N - 2} - 2C) for N = 3..9 (signed Euler numbers, C from the tangent-number series mod 2^200): %s" % res)
