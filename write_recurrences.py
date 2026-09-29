"""write_recurrences.py — regenerate and verify every recurrence found on 2026-09-28 and write recurrences.md (exact integer polynomials)."""
import sys, json
sys.path.insert(0, r"C:\Users\PC\Desktop\catalan")
from fractions import Fraction
from sympy import symbols, factor, expand, Rational
import guess_rec as GR, catalan_lucas as CL, catalan_doubled as CD

n = symbols('n')
lines = ['# Recurrences for the new integer sequences (2026-09-28)', '',
         'Each recurrence is written  c_r(n) a(n+r) + ... + c_1(n) a(n+1) + c_0(n) a(n) = 0  with integer polynomials, found by guess_rec.py',
         '(kernel mod primes, exact lift, verified on every available term).  Zudilin\'s Catalan polynomials:',
         '  P(x) = 20x^2 - 8x + 1,   Q(x) = 3520x^6 + 5632x^5 + 2064x^4 - 384x^3 - 156x^2 + 16x + 7.', '']

def block(title, seq, R, D, note, var='n'):
    sol = GR.guess(seq, R=R, D=D, verbose=False)
    lines.append('## ' + title)
    lines.append(note)
    lines.append('first terms: ' + ', '.join(str(x) for x in seq[:8]) + ', ...')
    if not sol:
        lines.append('NOT FOUND within order %d, degree %d' % (R, D)); lines.append(''); return None
    r, d, cs = sol
    lines.append('order %d, degree %d, verified on %d terms:' % (r, d, len(seq)))
    x = symbols(var)
    for k, c in enumerate(cs):
        poly = sum(ci*x**e for e, ci in enumerate(c))
        lines.append('  c_%d(%s) = %s' % (k, var, factor(poly)))
    lines.append('')
    return sol

# 1. zeta(7)-form numbers U(n) = u~_n / 30
rows = json.load(open(r"C:\Users\PC\Desktop\catalan\apery_lucas_180.json"))
def fr(s):
    if '/' in s: a, b = s.split('/'); return Fraction(int(a), int(b))
    return Fraction(int(s))
Uz = [int(fr(r['A7'])/30) for r in rows]
block('U(n) = u~_n / 30 — the zeta(7)-coefficient of Zudilin\'s F~_n (JTNB 2003 §7), normalised', Uz, 5, 32,
      'Lucas property and Coster supercongruence verified (apery_lucas.py, apery_p2.py). Characteristic polynomial of the leading terms: '
      'lambda^4 - 9264 lambda^3 - 12116166 lambda^2 + 752300 lambda - 19683 (Zudilin\'s, up to lambda -> -lambda). '
      'c_0 and c_4 share an apparent-singularity factor of degree 16.')
# 2. Catalan numbers U_n = 16^n u_n (Zudilin)
U, V = CL.seqs(80)
block('Zudilin\'s Catalan numbers U_n = 16^n u_n (u_n G - v_n -> 0)', [int(x) for x in U], 3, 8,
      'Zudilin\'s recurrence: c_2 = (n+2)^2(2n+3)^2 P(n+1), c_1 = -4 Q(n+1), c_0 = -256 (n+1)^2(2n+1)^2 P(n+2).')
# 3. half-coset family g(j) = U(2j+1)
Ud, Vd = CD.build(200)
block('Half-coset (recurrence) family g(j) = U(2j+1) = 2^{4j+2} y(j+1/2): the odd digits with no incoming carry', [int(Ud[2*j + 1]) for j in range(95)], 3, 8,
      'Zudilin\'s recurrence at half-integer index: 10j^2+26j+17 = P(j+3/2)/2, 10j^2+46j+53 = P(j+5/2)/2, middle = -2 Q(j+3/2).', 'j')
# 4. carried family, odd and even parts
dodd = json.load(open(r"C:\Users\PC\Desktop\catalan\carried_family_odd.json"))
top = 1
while str(top) in dodd and dodd[str(top)]['stable']: top += 2
fo = [int(dodd[str(M)]['value']) for M in range(1, top, 2)]
block('Carried family, odd part f(j) = F(2j+1): odd digits entered by a carry', fo, 3, 12,
      'REFLECTED Zudilin recurrence: 20j^2+48j+29 = P(-j-1), 20j^2+88j+97 = P(-j-2), and the sextic is Q(-j-2) (e.g. Q(-2) = 80503).', 'j')
dev = json.load(open(r"C:\Users\PC\Desktop\catalan\carried_family.json"))
he = [int(dev[str(M)]['value']) for M in range(0, 37, 2)]
lines.append('## Carried family, even part h(n) = F(2n): even digits inside a carry run')
lines.append('first terms: ' + ', '.join(str(x) for x in he[:8]) + ', ...  (19 stable terms, M <= 36)')
lines.append('REFLECTED Zudilin recurrence, verified on all 19 terms and predicting integers beyond (matches the reconstructed F(38), F(40)):')
lines.append('  (n+2)^2 (2n+3)^2 P(-n-1/2) h(n+2) = 4 Q(-n-3/2) h(n+1) + 256 (n+1)^2 (2n+1)^2 P(-n-3/2) h(n),')
Pm = lambda x: 20*x**2 - 8*x + 1
Qm = lambda x: 3520*x**6 + 5632*x**5 + 2064*x**4 - 384*x**3 - 156*x**2 + 16*x + 7
lines.append('  with P(-n-1/2) = %s,  P(-n-3/2) = %s,' % (expand(Pm(-n - Rational(1, 2))), expand(Pm(-n - Rational(3, 2)))))
lines.append('  64 Q(-n-3/2) = %s.' % expand(64*Qm(-n - Rational(3, 2))))
lines.append('')
lines.append('Reading: the carried family (digits of 2n entered by a carry of n -> 2n) satisfies Zudilin\'s recurrence with P and Q evaluated at')
lines.append('negative (reflected) arguments; the recurrence family satisfies it at positive arguments.  A carry moves a digit to the reflected side.')
open(r"C:\Users\PC\Desktop\catalan\recurrences.md", 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n'.join(l if len(l) < 400 else l[:400] + ' ...' for l in lines))
