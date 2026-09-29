r"""Uniformizing at infinity, part B (2026-09-29): the best conformal radius of the surface on which the linear forms live.
For Gamma' = <A, P_1, ..., P_r> (A = translation by 1, P_i parabolics at the harmless cusp representatives), the Green's function
of X' = H/Gamma' with the harmless cusps filled, pole at the filled cusp oo, is 2 pi E_{Gamma'}(tau, 1) (converges: infinite
covolume).  Its constant term gives the conformal radius in the coordinate q = e^{2 pi i tau}:
    log R* = 2 pi^2 * S,   S = sum over double cosets Gamma'_oo \ Gamma' / Gamma'_oo, gamma != 1, of 1/c(gamma)^2.
The CDT/Apery gate in these terms: log R* > tau (the denominators grow like e^{tau n}).
S is summed over reduced words P-syllable (A-syllable P-syllable)*, the last P-syllable in closed form:
    sum_{m != 0} 1/(c + d w m)^2 = (pi^2/sin^2(pi theta) - 1/theta^2)/(d w)^2,  theta = c/(d w)   (or 2 zeta(2)/(d w)^2 if c = 0).
Usage: python uniformize_radius.py w X   (Gamma' = <A, (1 0; w 1)>, words cut at |c| > X)"""
import sys, math
import mpmath as mp

def S_two_parabolics(w, X, kmax=None):
    """S for <A, B_w>, B_w = (1 0; w 1).  Enumerate prefixes B^{m1} A^{k1} ... B^{m_j} A^{k_j} with |c| <= X."""
    pi2 = math.pi**2
    def last_sum(c, d):
        # sum over m != 0 of 1/(c + d w m)^2 ; here (c, d) = bottom row of the prefix ending in A^k
        if d == 0: return 0.0
        dw = d*w
        if c % dw == 0:
            # theta integer: the m = -c/(dw) term would be c + dwm = 0 -> impossible in a free group except identity; skip it
            th = c//dw
            return (2*pi2/6)/dw**2 - (1.0/(dw*th)**2 if th != 0 else 0.0) if th != 0 else (2*pi2/6)/dw**2
        th = c/dw
        return (pi2/math.sin(math.pi*th)**2 - 1.0/th**2)/dw**2
    total = 2*(pi2/6)/w**2                    # single syllables B^m
    count = 0
    # stack of (c, d) bottom rows of prefixes ending in a B-syllable
    stack = []
    for m in range(1, X//w + 1):
        for s in (1, -1):
            stack.append((w*m*s, 1))           # B^m = (1 0; wm 1): bottom row (wm, 1)
    while stack:
        c, d = stack.pop()
        # append A^k (k != 0): (c, d) -> (c, d + c k)
        kr = (X + abs(d))//max(abs(c), 1) + 2
        for k in range(-kr, kr + 1):
            if k == 0: continue
            d2 = d + c*k
            # then B^{m}: bottom row (c + d2 w m, d2); sum over the final m in closed form
            total += last_sum(c, d2)
            count += 1
            # continue: B^m with |c + d2 w m| <= X
            if d2 == 0: continue
            mlo = math.ceil((-X - c)/(d2*w)) if d2*w > 0 else math.ceil((X - c)/(d2*w))
            mhi = math.floor((X - c)/(d2*w)) if d2*w > 0 else math.floor((-X - c)/(d2*w))
            for m in range(mlo, mhi + 1):
                if m == 0: continue
                c3 = c + d2*w*m
                if c3 == 0: continue
                stack.append((c3, d2))
    return total, count

if __name__ == '__main__':
    w = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for X in [int(a) for a in sys.argv[2:]] or [200, 800, 3200]:
        S, cnt = S_two_parabolics(w, X)
        print('w=%d X=%d  S=%.8f  (%d prefixes)   log R* = 2 pi^2 S = %.6f   R* = %.4f' % (w, X, S, cnt, 2*math.pi**2*S, math.exp(2*math.pi**2*S)))
        sys.stdout.flush()
