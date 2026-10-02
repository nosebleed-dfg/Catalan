# Carry correlations, decimal periods, and class numbers — run 1 (2026-09-22)

## Setup
Dedekind sum as a carry correlation:  s(h,k) = Σ_{a=1}^{k-1} ((a/k))((ha/k)),  ((x)) = frac(x) − 1/2.
With h = 10^j and k = p, 10^j·a mod p is the remainder j steps into the long division of a/p.
So s(10^j, p) = correlation between a long-division remainder and the remainder j carries later.

Sanity: carry form = cotangent form = reciprocity form, exactly, all h,k < 60.

## Result 1 — orbit identity (verified exactly, all primes < 700)
Let n = decimal period of 1/p, m = (p−1)/n. For each long-division cycle C of remainders, let D_C = digit sum of that repetend.

  Σ_{j=0}^{n−1} s(10^j, p)  =  Σ_C ((D_C − 9n/2)/9)²  =  (1/m) Σ_{χ odd, χ(10)=1} |B_{1,χ}|²

- Middle: squared deviation of repetend digit sums from the fair value 4.5 per digit.
- Right: B_{1,χ} = (1/p)Σ a·χ(a), the "power spectrum" (|L(1,χ)|² up to π²/p).
- Even period ⇒ sum = 0 (the rotation hits −1: Midy's 9-complement halves). Odd period ⇒ strictly > 0.
- Same vanishing holds on composite base-10 cyclotomic moduli Φ_n(10), n even (11, 101, 91, 10001, 9091).

## Result 2 — period (p−1)/2, p ≡ 3 mod 4
Orbit sum = h(−p)²/2, so the digit sum of the period of 1/p is (9/2)(n ± h(−p)).
Checked: 31, 43, 67, 71, 83, 107, 163, 191, 199, 227, 283, 307, 311, 347, 359.
Example 1/31: period 15, digit sum 54 = 67.5 − 13.5, h(−31) = 3.

## Result 3 — class number bound (derived, via analytic class number formula, Q=1, w=2)
AM–GM on Result 1: for odd period, orbit sum ≥ 2·(h⁻/2)^{4/m}, where h⁻ is the relative class number of the index-m subfield of Q(ζ_p). Equality at m = 2 and m = 4.

## Novelty check
Core connection (digits of 1/p ↔ Dedekind sums ↔ class numbers) is K. Girstmair's program: Digit variance and Dedekind sums (1997); variance for period (p−1)/2 (Ramanujan J., 2026, arXiv 2601.08416); mean values and variances for period (p−1)/2^m (arXiv 2606.29930, June 2026).
Possibly beyond those abstracts (needs full-text check): odd periods with index not a power of 2 (e.g. 37 idx 12, 79 idx 6, 239 idx 34, 271 idx 54); the sum-of-squares/Parseval form; the AM–GM bound; composite cyclotomic moduli Φ_n(10).

## Orbit sums on Φ_n(10), odd n
n=1 (9): 14/27 · n=3 (111): 55/6 · n=5 (11111): 2255/2 · n=7 (1111111): 226325/2 · n=9 (1001001): 101750 · n=11: 2263374465/2

---
# Run 2 — from a sum to a product (determinant)

Cycle deviation for the long-division cycle through a:  c(a) = Σ_{r in cycle} ((r/p)) = (D_cycle − 9n/2)/9.
For odd period n, index m = (p−1)/n, pick representatives a_1..a_{m/2} of the cycles modulo ±1 and set N_ij = c(a_i / a_j).

  h⁻(K) = 2·|det N|,   K = subfield of Q(ζ_p) fixed by ⟨10⟩ (degree m)

Eigenvalues of N (up to 2) are the individual B_{1,χ}; the determinant is their product = relative class number.
Verified against PARI class numbers (independent computation) for all 39 odd-period primes 31…719 with index ≤ 16,
incl. index 6 (79: 5, 547: 36, 643: 9) and index 12 (37: 1, 613: 2425).

Novelty: this is Maillet's determinant (Carlitz–Olson 1955) as generalized by Girstmair to any imaginary abelian field;
the digit/cycle reading is his 1994 "digits of 1/p and class number factors" line (extended by Hirabayashi 2005). Known.

Structural wall: the same machinery on the Catalan side needs B_{2,χ} with χ odd, which is identically 0.
The digit/carry determinants cannot see L(2, χ₋₄) = G at all.

---
# Run 3 — Catalan via walk + inverse + carry (all verified numerically)

1. Reflection (walk + inverse = constant + carry):  Li₂(z) + Li₂(1−z) = π²/6 − log z·log(1−z).
   Bloch–Wigner D(z) = Im Li₂(z) + arg(1−z)·log|z| (carry folded in): D(1/z) = D(1−z) = −D(z), D(i) = G.
   Volumes: 4G = Whitehead link complement; (3√3/2)·L(2,χ₋₃) = figure-eight knot complement.
2. Pure inverse:  G·(6/π²) = 1/2 + Σ over primitive Pythagorean triples 1/c².   (agrees to ~4e-9)
3. Three-step form:  G = (π²/12) · ∏_{p≡1 mod 4} (p²+1)/(p²−1).   (primes to 10^7, agrees to ~5e-9)
   π²/12 = 1 − 1/4 + 1/9 − … ; (product − 1)/2 = the triple sum.
4. Better use of π (Ramanujan):  G = (π/8)·ln(2+√3) + (3/8)·Σ_{n≥0} 1/((2n+1)²·C(2n,n)).   (30 digits)
   π/4 = L(1,χ₋₄), ln(2+√3)/√3 = L(1,χ₁₂), and χ₋₄·χ₁₂ = χ₋₃: the three quadratic subfields of Q(ζ₁₂).
   Carry = product of level-1 values from your basis (π, ln(2−√3)); tail = Apéry-shaped, shrinks ~4× per term.
   Open for Lemma E: tail alone converges too slowly vs. its denominators to prove irrationality; needs an Apéry-style partner sequence.

---
# Run 4 — Pythagorean sort in base 9, carry mistakes in G

- Base 9 is carry-proof for the Pythagorean test: 9 ≡ 1 mod 8, so n ≡ (base-9 digit sum of n) mod 8 and every carry changes the digit sum by −8.
  Whether n is ≡ 1 or 3 mod 4 (hypotenuse prime or not, the χ₋₄ sign in G's series) is read off the digit sum; carries never touch it. Base 10 has no such property.
- G in base 9 (12,000 digits): 0.821657603356012411640328675581512826616725250646082732236313…
  Digit counts 1289–1368 each (expected ~1333). Longest run of 0s: 4 (at digit 10831). Longest run of 8s: 3.
- Carry mistakes: Ramanujan partial sums checked to ~2960 base-9 digits. 68 events where the value was right to E digits but the digit string only to E−2 or E−3.
  Every one sits exactly on a 0 (or 00) in G just past the frontier: partial sums come from below, so they read …88 where G reads …00.
  Depth = length of that 0-run + 1. Fully predicted by G's own digits; they never compound.
- Certified continued fraction: 11,085 terms (largest 9075). If G = a/b, then b has more than 5,724 decimal digits.

---
# Run 5 — the 10011001-style transform for G

Sign pattern of G's series (0, +1, 0, −1, repeating) as a polynomial: x − x³.  Generating function: Σ χ₋₄(n)xⁿ = x/Φ₄(x).
All three characters of Q(ζ₁₂):
- χ₋₄ ↔ x/Φ₄(x)                 at x=1/10: 10/101   = 0.0990 0990…   (Φ₄(10) = 101)
- χ₋₃ ↔ x/Φ₃(x)                 at x=1/10: 10/111   = 0.090 090…     (Φ₃(10) = 111)
- χ₁₂ ↔ x(1−x⁴)(1−x⁶)/(1−x¹²) = x(1−x²)/Φ₁₂(x)   at x=1/10: 990/9901 (Φ₁₂(10) = 9901, prime, period 12)
Ladder: 10/9 (all +, ζ), 10/11 (alternating), 10/101 (quarter-turn, Catalan's character).
Level 0 → level 2 (xⁿ → 1/n²):  G = −∫₀¹ ln x / Φ₄(x) dx = ½∫₀^∞ t/cosh t dt   (25 digits)

---
# Run 6 — the carry, z-wise (carry-twisted Catalan in base 9)

Weight each n by z^(base-9 digit sum of n) instead of zⁿ.  A base-9 carry drops the digit sum by 8, so each carry costs a factor z^(−8).
Carry-twisted Clausen:  C_z = Im Σ z^{s₉(n)}/n²,  z = e^{iθ}.  Carry cost = C_z − Cl₂(θ).
Computed exactly via the base-9 self-similar recursion (n = 9m + d), checked against Li_s(i) to 1e-31.

| θ | carry-twisted | plain Cl₂(θ) | carry cost |
|---|---|---|---|
| 90° (Catalan) | 0.915965594177219 | 0.915965594177219 | 0 (to 1e-31) |
| 45° | 0.981872151050203 | same | 0 |
| 10/9° | 0.065583662210746 | 0.095854872433802 | −0.0302712 |
| 20/9° | 0.130528528232848 | 0.164826573728388 | −0.0342980 |
| 40° (z⁹=1) | 0.953764145633952 | 0.953741016638851 | +0.0000231 |
| 60° | 1.038657636077585 | 1.014941606409654 | +0.0237160 |

- Carry cost is exactly 0 iff z⁸ = 1 (multiples of 45°): carries are free. Catalan sits in that set; in the literature this is the "improper" case of strongly 9-multiplicative sequences (Alkauskas 2004; Allouche–Mendès France–Peyrière 2000 framework).
- 81-fold split of G at 10/9° = roots-of-unity filter keeping only n ≡ 0 mod 81: the carry positions.
- Observed, unexplained: at 9th roots of unity (40°, 80°, 160°) the self-term of the recursion vanishes and the carry cost is small (2.3e-5, −1.2e-3, −3.5e-4).

---
# Run 7 — primes sorted by base-9 last digit × Pythagorean class (mod 36), up to 10^8

In base 9: p mod 9 = last digit; p mod 4 = digit sum mod 4 (carry-proof). Together: p mod 36.
E(x) = (log x/√x)(12·π(x;36,a) − π(x)); theory predicts mean 1 − c(a), c(a) = # square roots of a mod 36.

| a | last digit | Pythagorean | count | predicted | observed (log-avg 1e4–1e8) |
|---|---|---|---|---|---|
| 1 | 1 | yes | 479929 | −3 | −3.32 |
| 19 | 1 | no | 480030 | +1 | +0.37 |
| 29 | 2 | yes | 480067 | +1 | +0.68 |
| 11 | 2 | no | 480273 | +1 | +1.86 |
| 13 | 4 | yes | 480058 | −3 | −3.47 |
| 31 | 4 | no | 480302 | +1 | +1.31 |
| 5 | 5 | yes | 480073 | +1 | +1.44 |
| 23 | 5 | no | 480205 | +1 | +1.21 |
| 25 | 7 | yes | 480034 | −3 | −3.97 |
| 7 | 7 | no | 480164 | +1 | +1.75 |
| 17 | 8 | yes | 480343 | +1 | +1.28 |
| 35 | 8 | no | 479975 | +1 | +0.27 |

Non-Pythagorean ahead of Pythagorean (share of sampled x): digit 1: 80.6%, 4: 87.8%, 7: 92.2% | digit 2: 55.2%, 5: 50.0%, 8: 40.8%.
The whole Chebyshev bias against Pythagorean primes sits on last digits 1, 4, 7 (the classes p ≡ 1 mod 12 = primes splitting completely in Q(ζ₁₂)).
On digits 2, 5, 8 the race is fair. Consistent with Rubinstein–Sarnak square-class theory; the base-9 reading is the new presentation.

---
# Run 8 — inverting out: Möbius and Liouville over the Gaussian integers Z[i]

Inverse of ζ(s)L(s,χ₋₄) = ζ_{Q(i)}(s) is Σ μ_{Z[i]}(a)/N(a)^s over Gaussian ideals.
- Gaussian-squarefree density = 1/ζ_{Q(i)}(2) = 6/(π²G) = 0.663701.  Counted by Möbius inclusion–exclusion over ideals up to norm 10^6: 0.66375. ✓
- Gaussian Liouville λ_K(a) = (−1)^{Ω_K(a)}, with Ω_K(a) = e₂ + Σ_{p≡1(4)} e_p + Σ_{p≡3(4)} e_p/2 (e_p = exponent of p in N(a)). Depends only on the norm.
  Dirichlet series ζ_K(2s)/ζ_K(s); pole at s=1/2 gives L_K(x) = Σ_{N(a)≤x} λ_K(a) (over ideals) ~ (π/4)/ζ_K(1/2)·√x = −0.8055√x,
  with ζ_K(1/2) = ζ(1/2)·L(1/2,χ₋₄) = (−1.46035)(0.66769) = −0.97507.
- Numerics to norm 10^6 (elements; ideals = /4): L/√x at 10^3..10^6: −2.15, −3.20, −2.95, −1.11; range on [10^3,10^6]: −7.14 to +0.84.
- Gaussian Pólya problem fails early: L_K(x) > 0 first near x ≈ 2210 (vs 906,150,257 for Z), max +28 at x = 466,370, last positive (≤10^6) at 485,810.
- Closest literature found: Humphries–Shekatkar–Wong 2019 (Liouville refined to primes in arithmetic progressions). The Z[i]-lattice version with r₂ weighting: not found.

---
# Run 9 — forward (i^{z²}, quadratic) and the new backward (primes mod 9, cubic)

Forward: quadratic Gauss sum Σ e^{2πi n²/p} = √p (p≡1 mod 4) or i√p (p≡3 mod 4). Clean; i pops out exactly for non-Pythagorean primes.
Backward: cubic Gauss sums g(π)/√p for primes p ≡ 1 mod 3 (π primary in Z[ω]), 5634 primes to 120,000.
- |g/√p| = 1 to 1e-15. Angles: mean cos(3θ) = 0.000 (tripled angle equidistributed, Heath-Brown–Patterson).
- Bias: mean Re(g/√p) = +0.140 overall; x^{-1/6} scale = 0.142. Cumulative Σ Re(g/√p): 123, 301, 495, 786 at x = 1e4, 3e4, 6e4, 1.2e5
  (Patterson's bias, proved under GRH by Dunn–Radziwiłł, Annals 2024; first-order prediction ~1470 at 1.2e5, same order, slow convergence).
- Sorted by p mod 9: mean Re = +0.136 (1), +0.144 (4), +0.139 (7). The cubic bias is uniform across p mod 9.
  Contrast: the quadratic (Chebyshev) bias split cleanly by p mod 36 (Run 7). The backward step is class-blind mod 9.

---
# Run 10 — the anchor at the negative rungs: von Staudt–Clausen

ζ(1−k) = −B_k/k, and B_k = (integer) − Σ 1/p over primes p with (p−1) | k, i.e. primes whose rotation order divides k.
k=2: primes 2,3 → −1/12 = −(1 − 1/2 − 1/3)/2.  k=4: 2,3,5.  k=6: 2,3,7.  k=12: 2,3,5,7,13; numerator 691.  k=14: ζ(−13) = −1/12 again.
Integer parts: 1,1,1,1,1,1,2,−6,56,−528,6193,−86579 (k=2..24).
Denominators are always squarefree: 9 never appears; μ(denominator) ≠ 0. The inverse is counted by rotation order, never by size, and never overflows.
Numerators carry the irregular primes (691, 3617, 43867, …): Kummer — p divides a numerator ⟺ p divides the class number of Q(ζ_p).
Squaring side: ζ(2s)/ζ(s) vanishes at every negative odd s (trivial zeros of ζ(2s)); the μ side gives Eisenstein integers 1/ζ(−1) = −12 → −24, 240, −504.

---
# Run 11 — −1/12 on the i side: the rung is a zero, the slope is a volume (verified 20 digits)

L(−1,χ₋₄) = 0, so ζ_{Q(i)}(−1) = (−1/12)·0 = 0. The information moves to the slope:
  L'(−1,χ₋₄) = 2G/π;   ζ_{Q(i)}'(−1) = (−1/12)(2G/π) = −G/(6π) = −0.0485935
  L'(−1,χ₋₃) = (3√3/4π)·L(2,χ₋₃);   ζ_{Q(ω)}'(−1) = −(√3/16π)·L(2,χ₋₃) = −0.0269222
Uniform: ζ_K'(−1) = −vol(H³/PSL₂(O_K))/(2π) = −vol(link complement)/(24π), with vol(Whitehead) = 4G, vol(figure-eight) = (3√3/2)L(2,χ₋₃); both link groups are index-12 subgroups.
Compare Z: ζ'(−1) = 1/12 − ln A (Glaisher A = 1.28243). Gaussian analogue of Glaisher's constant: e^{G/(6π)} = 1.04979.
Reframe: G irrational ⟺ the Gaussian lattice zeta's slope at the −1/12 rung is not a rational multiple of 1/π.

---
# Run 12 — what counting in base b does to the zeros (first-order response)

Carry-twisted zeta F_b(θ,s) = Σ e^{iθ·s_b(n)} n^{-s} (s_b = base-b digit sum). At θ=0 it is ζ(s). Each zero ρ moves with velocity
  dρ/dθ = −i·D_b(ρ)/ζ'(ρ),   D_b(s) = Σ s_b(n) n^{-s}  (continued to Re s = 1/2 via the base-b self-similar recursion; checked at s=3 and 2.5+3i).
Re(dρ/dθ) = push off the critical line; Im = slide along it.

| zero | base 2 | base 9 | base 10 |
|---|---|---|---|
| ρ₁ = 1/2+14.13i | −0.222 + 0.904i | +0.833 + 1.345i | +0.773 + 2.178i |
| ρ₂ = 1/2+21.02i | +0.088 + 0.629i | +0.874 + 2.323i | +0.877 + 1.409i |
| ρ₃ = 1/2+25.01i | +0.068 + 0.905i | +0.282 + 1.616i | +0.433 + 2.905i |

Bases 9 and 10 push all three zeros to the right (off the line); base 2 pushes the first left, the next two right, and much more weakly.
Not found in the literature; the twisted series itself is Alkauskas-type, the zero-response is the new quantity.

---
# Run 12 (revised, stable continuation) — zero response to base-b carry twist, zeros 1–4
Velocity v = dρ/dθ (Re = push off the line, Im = slide along it) and second-order Re ρ''(0) (net drift under ±θ alternation).
Slide is the base clock: Im v ≈ (b−1)/(2 ln b) (mean digit sum). Push is the carry fluctuation: it grows with b.

| base | v(ρ₁) | v(ρ₂) | v(ρ₃) | v(ρ₄) | Re ρ''(0) at ρ₁..ρ₄ | mean Re v | (b−1)/(2 ln b) |
|---|---|---|---|---|---|---|---|
| 2 | -0.2218 +0.9041i | +0.0876 +0.6294i | +0.0683 +0.9048i | -0.3999 +0.6951i | -0.0199, -0.1156, -0.1015, -0.3039 | -0.12 | 0.72 |
| 3 | +0.4455 +1.3623i | -0.3438 +0.8633i | -0.0333 +0.8417i | +0.2557 +0.7773i | -0.4163, -0.2137, -0.3075, -0.1790 | +0.08 | 0.91 |
| 4 | +0.0180 +2.1868i | +0.5701 +1.1676i | +0.2916 +1.2414i | -0.9215 +0.8899i | +0.4953, -1.0907, -1.0562, -1.2749 | -0.01 | 1.08 |
| 5 | +0.3661 +1.2143i | +0.8245 +2.3096i | +0.6736 +2.2232i | +0.0347 +0.9821i | +0.1055, +0.0283, -0.7770, -0.6619 | +0.47 | 1.24 |
| 6 | +0.4126 +1.8359i | +0.3417 +1.8887i | -0.0315 +2.9827i | +1.2739 +1.4501i | +1.4782, +1.6783, +1.6380, -3.5880 | +0.50 | 1.40 |
| 7 | +0.5655 +1.8290i | +0.8767 +1.8932i | +0.0701 +1.5427i | +1.4597 +2.8487i | +0.3807, -0.6239, +0.2832, -0.3625 | +0.74 | 1.54 |
| 8 | +0.5612 +1.1581i | +0.6518 +1.6470i | +0.3800 +2.7660i | +0.7351 +2.6176i | +0.9995, +2.8745, +1.5086, -1.4800 | +0.58 | 1.68 |
| 9 | +0.8326 +1.3449i | +0.8741 +2.3234i | +0.2822 +1.6157i | +1.3973 +1.9245i | +3.6722, +0.4158, +1.0542, -0.7916 | +0.85 | 1.82 |
| 10 | +0.7726 +2.1785i | +0.8767 +1.4090i | +0.3856 +2.8664i | +1.0251 +2.8567i | +2.0592, +0.8413, +1.6617, -4.8340 | +0.77 | 1.95 |
| 11 | +0.8923 +1.8382i | +1.0801 +1.9938i | +0.6768 +2.0259i | +1.4958 +2.0233i | +0.9900, +6.0973, +0.6158, +0.7031 | +1.04 | 2.09 |
| 12 | +0.9715 +1.3466i | +1.0620 +2.3764i | +0.5985 +1.8379i | +1.3645 +2.2901i | +1.6331, +1.0631, +3.3198, -4.1311 | +1.00 | 2.21 |
| 16 | +1.3927 +2.2330i | +1.3909 +2.4244i | +1.3213 +2.4563i | n/a | +5.0633, +2.9911, +14.6980 | +1.37 | 2.71 | (ρ₄ unstable, excluded)

Findings: (1) bases 2–4 push zeros off the line weakly and in mixed directions; from base 5 up every zero is pushed right, and the push grows with the base (mean Re v: 2: −0.12, 5: 0.47, 10: 0.77, 12: 1.00, 16: 1.37).
(2) The slide along the line matches the mean-digit-sum clock (b−1)/(2 ln b) within the Delange fluctuation, so the off-line push is purely the carry fluctuation.
(3) Second order: under move-then-correct (±θ), the even part is not zero. Base 2 (and mostly 3, 4) drifts LEFT at second order; bases 6 and up drift RIGHT; the crossover is near base 5. The twist is closest to pure sine (drift-free) around base 5.
Earlier Run 12 numbers for ρ₃ were ~10% off (shallow continuation); this table supersedes it.

---
# Run 13 — carries over ALL bases: the divisor sum, ζ(1/2), and a new set of zeros

Per step: carries when going n−1 → n, summed over every base b ≥ 2, = Σ_b v_b(n) = #{(b,j): b^j | n} = Σ_{j≥1}(τ_j(n) − 1), τ_j = number of j-th-power divisors.
Dirichlet series: ζ(s)·G(s), G(s) = Σ_{j≥1}(ζ(js) − 1) (the perfect-power series; Goldbach–Euler 1737: Σ_{j≥2}(ζ(j)−1) = 1).
Cumulative: C(n) = Σ_{b≥2}Σ_{j≥1}⌊n/b^j⌋ = D(n) − Σ_{j≥2}Σ_{b≥2}{n/b^j}  (exact), so
  C(n) = n ln n + (2γ − 1)n + ζ(1/2)√n + ζ(1/3)n^{1/3} + … + Δ(n)   [checked: n=10^5, D−C = 530 vs −ζ(1/2)√n − ζ(1/3)n^{1/3} − … ≈ 545]
  The square-position carries contribute ζ(1/2)√n = −1.4604√n (Pólya's constant); Δ(n) = Dirichlet divisor error term (spectral).
Carry zeros (zeros of G):
  real: 0.68046 (in (1/2,1)), 0.40113, 0.28552, 0.22186 (one in each (1/(j+1),1/j)).
  complex, 0.35<Re<2.3, Im<40: exactly five (winding number 5): 0.928+5.928i, 0.508+13.625i, 1.974+23.094i, 1.134+31.614i, 1.398+38.969i.
  40<Im<80 (Newton hunt, completeness not certified): 0.670+42.171i, 1.456+52.063i, 0.502+59.989i, 0.661+66.109i, 1.590+68.554i.
  Not on any line; real parts scatter from 0.50 to 1.97. Two sit at Re ≈ 0.50 (Im 13.62, 59.99). They are NOT the 1-points of ζ (1.408+23.33i, 0.411+31.72i, 1.368+38.42i).
  No literature found on the zeros of G.

Third layer (Run 13b): D(n) − C(n) = Σ_{j≥2}Σ_{b≥2} {n/b^j} splits exactly into
  continuous tails Σ_{b^j>n} n/b^j  +  anchor debts (1/2 per perfect-power pair b^j ≤ n)  +  sawtooth Σ_{b^j≤n} ((n/b^j)).
  n = 10^5: 352.0 + 203.5 + (−26.4) = 529.1 vs 530. The sawtooth ((x)) = −Σ sin(2πkx)/(πk) is the Dedekind kernel from Run 1: a sum of sines with one period per perfect power.
  Per carry position: ζ(1/j) = −1/(j−1) − 1/2 + r_j, with r_j = 0.040, 0.027, 0.020, 0.016 (j=2..5): pole (tail) + anchor (ζ(0) = −1/2) + Bernoulli remainder.

Pressing the 0.04 (Run 13c): r_j = −s∫₁^∞ ((x)) x^{−s−1} dx at s = 1/j — the Mellin transform of the sawtooth. Verified: 0.03965, 0.02664, 0.02005 (j=2,3,4).
  Bernoulli expansion: r_j = (1/12)(1/j) − (1/720)(1/j)(1/j+1)(1/j+2) + … ; leading term 1/24 at j=2, so the 0.04 is the −1/12 constant divided by the carry position.
  Catalan analogue: no pole, so no tail: β(1/2) = 1/2 + (1/2)∫₀^∞ (−1)^{⌊x⌋}(2x+1)^{−3/2}dx = 1/2 + 0.16769 — a SQUARE wave (pure over/under sign) where ζ has a sawtooth (ramp).
  Gaussian Pólya constant ζ(1/2)·β(1/2) = (−1 − 1/2 + 0.0396)(1/2 + 0.1677) = −0.97507.

---
# Run 14 — the 12-fold: carry cost on the 12th-root lattice (level 2), bases 2–13 and 14, 17, 21, 25

Per base: invisible rotations z^{b−1}=1 (cost exactly 0), dead rotations (digit polynomial (z^b−1)/(z−1)=0; self-similar term killed, cost small). On the 12-lattice these are Φ_d for d | gcd(b−1,12) and d | gcd(b,12).

| base | dead gcd(b,12) | invisible gcd(b−1,12) | 30° | 60° | 90° | 120° | 150° |
|---|---|---|---|---|---|---|---|
| 2 | 2 | 1 | +0.07488 | +0.31353 | +0.36264 | +0.29512 | +0.16249 |
| 3 | 3 | 2 | +0.14354 | +0.23140 | +0.12978 | +0.00504 | −0.04613 |
| 4 | 4 | 3 | +0.15697 | +0.12269 | +0.00197 | 0 | +0.04211 |
| 5 | 1 | 4 | +0.14372 | +0.03967 | 0 | +0.03524 | −0.00981 |
| 6 | 6 | 1 | +0.11817 | +0.00025 | +0.02984 | −0.00127 | +0.00745 |
| 7 | 1 | 6 | +0.08855 | 0 | +0.01947 | 0 | +0.00449 |
| 8 | 4 | 1 | +0.05992 | +0.01574 | −0.00132 | +0.01262 | −0.00354 |
| 9 | 3 | 4 | +0.03546 | +0.02372 | 0 | −0.00096 | +0.00695 |
| 10 | 2 | 3 | +0.01701 | +0.01785 | +0.01016 | 0 | −0.00365 |
| 11 | 1 | 2 | +0.00532 | +0.00634 | +0.00723 | +0.00639 | +0.00374 |
| 12 | 12 | 1 | +0.00006 | −0.00086 | −0.00088 | −0.00065 | −0.00034 |
| 13 | 1 | 12 | 0 | 0 | 0 | 0 | 0 |
| 14 (≡2) | 2 | 1 | +0.00327 | +0.00519 | +0.00505 | +0.00384 | +0.00205 |
| 17 (≡5) | 1 | 4 | +0.01402 | +0.00233 | 0 | +0.00256 | −0.00087 |
| 21 (≡9) | 3 | 4 | +0.00595 | +0.00419 | 0 | −0.00026 | +0.00117 |
| 25 (≡1) | 1 | 12 | 0 | 0 | 0 | 0 | 0 |

Fold confirmed: the zero / small / live pattern repeats with period 12 (17 matches 5, 21 matches 9, 25 matches 13, 14 matches 2). Magnitudes do not repeat: live costs shrink roughly like 1/b².
The spectrum closes at 12/13: base 12 has every lattice rotation dead (all costs < 0.001), base 13 has every rotation invisible (all exactly 0). Min = base 2 (nothing dead or invisible beyond the trivial half-turn).
Item 3 (base-9 dip at 9th roots) is resolved: those are zeros of the digit polynomial. Item 2's base-5 crossover: base 5 is the first base where the quarter-turn is invisible (4 | b−1).

Run 14b — closed form of the carry cost (leading order): cost_b(z) ≈ Im{ A(z)·[Li₂(z) − Li₂(z^b)] } / b², A(z) = (z^b−1)/(z−1) the digit polynomial.
  Reproduces the zero pattern exactly (A=0 dead; z^b=z invisible) and the 1/b² scale; magnitudes right to ~×1.5–2, signs sometimes off at 30°/150°.
  The normalizer (the "other half"): the full tower Σ_k C(−s,k) b^{−s−k} (z d/dz)^k A(z) · F(z, s+k) — every Euler-derivative of the digit polynomial, each coupled to a HIGHER level s+k.
  Truncating at 1–2 derivatives oscillates around the truth (alternating, ratio (b−1)/b). It is the Taylor tower of the Hurwitz shift (m + d/b)^{−s}: the exact object is the b-fold Hurwitz decomposition ζ(s) = b^{−s} Σ_{d<b} ζ(s, d/b), twisted by z^{digit sum}.
  So the pull-back does not live in the z-plane: it lives along the level axis s → s+k, i.e. in the digit-offset coordinate a = d/b of ζ(s, a).

---
# Check of "ζ(7) is irrational" (Prince Anand PDF), exact arithmetic

Context: ζ(5) irrationality preprint (A. Fauzan, 17 Sep 2026) has a complete Lean/Mathlib formalization (github.com/mo271/Zeta5, irrational_five without sorry, with deviations from the paper). The ζ(7) PDF copies that architecture (D_N^6 → D_N^8, pole functional shifted), defers the p-adic core to "[5] with modifications", and claims constants A_M ≈ −13.75 (normalizer), U = +1.76 (real side), net −12 per K².
Exact recomputation of the paper's own Δ_K(X), S_K, content, P_K = F_K/content (leading coefficient matches its Lemma 2.3; P_K(ζ(7)) > 0 as claimed):

| K | h | log F_K(ζ7)/K² | −log content(F_K)/K² | log P_K(ζ7) | per n² |
|---|---|---|---|---|---|
| 40 | 37 | −0.6554 | +1.4519 | +1274.4 | +1274 |
| 80 | 74 | −0.7073 | +1.5785 | +5575.8 | +1394 |

For the criterion one needs log P_K(ζ7) → −∞. It grows (net ≈ +0.87 per K² and worsening). Compare the verified ζ(5) data (mo271): −1.276 + 1.110 → net −0.17 at K=40, tending to −0.017.
The paper's real-side sign (+1.76) and normalizer sign (−13.75) are both contradicted by exact data (−0.7, +1.5): the determinant is not more divisible; it has denominators. The ζ(7) argument fails at the arithmetic normalization (Section 5). Script: zeta7_check.py.

---
# Extracting from ζ(5): a Catalan version of the Fauzan machine (built and tested)

Catalan analogue of the moment functional (all three identities verified numerically to 20 digits):
  weight w(t) = π t² cosh(πt)/sinh²(πt) > 0   (kernel 1/sinh(πt): same poles as 1/(e^{2πt}−1), alternating residues — the square-wave version)
  ∫₀^∞ w(t)/(t²+a²) dt = a·ζ_alt(2,a) − 1/(2a),  ζ_alt(2,a) = Σ(−1)^n/(n+a)²;  at a = j+1/2: = 4(−1)^j (G − β_j) ⇒ pole values affine in G
  ∫₀^∞ t^{2e} w(t) dt = (2^{2e+2} − 1)|B_{2e+2}|  (rational)
Hankel determinant Δ_K(X) with poles at −(2j+1)², E_N^m numerator, K = 40n, N = 3n (Fauzan's shape), primitive polynomial P_K:
  K=40, m=6:  log P_K(G) = +1591 (log Δ/K² = +0.77);  K=80: log P_K(G) = +6718 (+1.30/K²). Larger m or N: worse (up to +40076).
  Compare verified ζ(5): −265, −833. So the direct transplant fails; same signature as the bad ζ(7) paper.
Mechanism: alternation halves the decay rate of the kernel (e^{−πt} vs e^{−2πt}) with the same pole spacing, so the external field driving the real-side decay is half as strong. That is the Catalan wall inside this architecture. Doubling pole density brings π² in as a second unknown.
Script: catalan_hankel.py.

Catalan machine, refined diagnosis (K=40): closeness log(P(G)/H(P))/K² = −2.34 (ζ(5): ≈ −2.2) — the analytic side is FINE, better than ζ(5).
Failure is height: log H(P)/K² = 3.33 (ζ(5): ≈ 2.0), rising to 4.5 at K=80. Content of Δ: numerator 2^917·3^404·5^146·7^71·11^38·13^15 (small primes: divisibility, good);
denominator primes 17..79 with exponents up to 55 (845 digits ≈ 1.2/K²) — from the residue factors ∏(b'−b)(b'+b) over odd b,b' ≤ 79. That ≈ the whole gap.
So the "half field" story was wrong; the wall for Catalan in this architecture is p-adic, at the mid-range primes K/2 < p < 2K.

Mid-prime scan and parameter scan (Catalan Hankel machine, half-integer poles). All numbers per K²; net = log P_K(G)/K² (needs < 0).
| config | K | closeness | height | mid-prime denominators | net |
|---|---|---|---|---|---|
| N=3n, m=6 (Fauzan shape) | 40/60/80 | −2.34 | 3.33/3.32/3.40 | 1.15/1.29/1.26 | +1.00 |
| b ≡ 1 mod 4 poles only, m=6 | 40 | −3.53 | 5.97 | 1.83 | +2.44 |
| m=3 / m=2 | 40 | −2.25/−2.20 | 2.82/2.64 | 1.18/1.19 | +0.57/+0.44 |
| N=3, m=1 | 40 | −2.15 | 2.53 | 1.25 | **+0.38** |
| N=4, m=1 | 60 | −2.16 | 2.60 | 1.42 | +0.45 |
| m=0 (pure Hankel) | 40/60/80 | −2.19 | 2.56/2.64/2.65 | 1.36/1.53/1.58 | +0.38/+0.45/+0.46 |
| N=10, m=1 / N=15, m=1 | 40 | −1.86/−1.57 | 2.35/2.12 | 0.96/0.78 | +0.48/+0.55 |
Verdict: the gap is stable at ≈ +0.4 per K² at the best settings (ζ(5) closes at −0.017). Mid-prime denominators scale like K² (a constant in the balance, not a K log K effect).
Basis tricks (reflection/palindrome) cannot help: the primitive polynomial is basis-invariant; only the pole set and functional matter. Restricting poles to one class mod 4 makes it worse.
Untested and requiring new Hermite-type formulas: poles at Gaussian norms (2D lattice), and poles on ζ₁₂ norms (12-fold).

Analytic press on the lattice/remainder route (2026-09-24):
- Norm-factorized lattice functional at integer poles: pole values are Σ_d χ₋₄(d) ψ'(1+n/d)/d², i.e. Hurwitz values at rational shifts n/d ⇒ every L(2,χ mod d) enters. Proliferation of unknowns; route closed.
- Prime-indexed poles (b = odd primes): closeness improves to −3.1 but height 4.9, net +1.77 (K=40), +2.10 (K=60); denominators from the larger pole values dominate. Sparse pole sets are worse.
- Invariant across all pole sets/weights tried: net ≥ +0.38 per K²; mid-prime cost is set by pole spacing (consecutive odd squares) and is ≈1.2.

---
# ζ(5) control (same pipeline) and the corrected diagnosis — factorization focus

Control reproduces the verified numbers exactly: log P_K(ζ5) = −265.1 (K=40), −833.3 (K=80) [mo271: −265, −833]. Pipeline validated.

| per K² | ζ(5), K=40 | ζ(5), K=80 | Catalan m=6, K=40 | Catalan m=1, K=40 |
|---|---|---|---|---|
| closeness log(P(ξ)/H) | −3.457 | −3.469 | −2.34 | −2.15 |
| height log H | 3.291 | 3.338 | 3.33 | 2.53 |
| net | −0.166 | −0.130 | +1.00 | +0.38 |
| content denominators, total | 1.18 | 1.29 | 1.22 | 1.95 |
| of which mid primes (K/2,2K) | 0.45 | 0.70 | 1.15 | 1.25 |
ζ(5) denominator primes at K=40: 5^31, 11^83, 13^93, 17^123, 19^113, 23^110, 29^60, 31^46, 37^4 (all ≤ K); numerator 2^324·3^87.
Catalan (m=6) denominator primes: 17..79, i.e. K/2 < p < 2K; numerator 2^917·3^404·5^146·7^71·11^38·13^15.

Correction: heights are the SAME (3.29 vs 3.33); the whole Catalan deficit is closeness: −3.46 vs −2.34 (1.1 per K²). The arithmetic side is a wash (same total denominator cost, shifted from small/medium primes in ζ(5) to mid primes for Catalan). So the wall IS the analytic side after all: the alternating kernel 1/sinh(πt) decays at half the rate of 1/(e^{2πt}−1); the machine's closeness scales with the field. Earlier "height" diagnosis withdrawn.
Implication: fixes must raise closeness (a full-rate kernel with G-affine pole values), not lower height. Any full-rate kernel found so far brings π².

---
# Oscillating factors (2026-09-24)
- The lowest nontrivial zero of L(s,χ₋₄) is at 1/2 + 6.0209489i (ζ's first is 14.1347). Catalan's character owns the fastest arithmetic oscillation: one cycle per factor e^{2π/6.02} ≈ 2.8 in x.
- Correction to Run 8: over ideals, the Gaussian Liouville sum first turns positive at x = 914 (exact; 2,210 came from a step-10 sample). The explicit formula (main term + 8 β-zeros + 3 ζ-zeros) tracks it: at x=2210 predicts −5.7 vs actual +2; at 6000, −25.5 vs −33.
- New testbed: Mertens over Z[i]. M_K(x) = Σ_{n≤x} (μ ∗ μχ₋₄)(n) computed to 10^7: M_K/√x = −0.03, −0.04, −0.14, +0.51, +0.47 at 10^3..10^7; no |M_K| > √x beyond x = 17. Ordinary Mertens max 0.57. No early counterexample.

---
# Deep dive on the ζ(5) proof (2026-09-24)

Status of the ζ(5) proof: two Lean projects. mo271/Zeta5 (Firsching) claims complete without sorry. danromik/zeta5-irrationality (Romik, Claude Code, with a 22-page referee audit): everything on the paper's route is machine-checked EXCEPT the logarithmic-energy inequality (6.14) — the real-side bound — left as a sorry; two classical axioms (Hermite's formula, PNT). The audit found one gap in the paper (b_a ≤ 6αx+3 on p.10, asserted without proof) and supplied the proof. So: arithmetic half verified; analytic half (6.14) is the soft spot, numerically checked at K=40,80,120 only.
Hidden mechanisms: (i) the distribution formula is Raabe's multiplication theorem m⁴τ(P)=Σ τ(P(a+mz)) — the functional is a p-adic measure (Kubota–Leopoldt in disguise; the Catalan analogue would be the p-adic L-function of χ₋₄); (ii) local bases reduce mod p to Hermite-interpolation bases at the mirror pairs ±a — unimodularity; (iii) the rational remainders −1/4 + 1/(2j) are load-bearing in the divided-difference congruence (the anchor terms matter p-adically); (iv) the pullback x⁵R(−x²) is antisymmetric under x ↦ −1−x (the mirror); (v) the real side is a transfinite-diameter/capacity statement (F. Brown 2026 gives the general criterion).

DECISIVE CONTROL: ζ(2) functional on the SAME half-integer poles with the full-rate kernel 1/(e^{2πt}−1) (moments 4^e|B_{2e+2}|/2; pole values (b/16)ζ(2,b/2) − 1/(8b) − 1/8; entries checked against the direct integral). Per K²:
| config | closeness | height | net |
|---|---|---|---|
| ζ(2), m=0, K=40 | −2.905 | 2.775 | **−0.130** |
| ζ(2), m=1, K=40 | −3.035 | 2.891 | **−0.145** |
| ζ(2), m=1, K=60 | −3.054 | 2.951 | **−0.103** |
| Catalan, m=1, K=40 | −2.15 | 2.53 | +0.38 |
Same poles, same numerators, same denominators structure: plain kernel closes, alternating kernel fails. The alternating twist alone costs ≈0.5–0.8 per K² of closeness. The obstruction is isolated to the sign carry in the kernel.

---
# Tail-defined (pole-only) G² functional — tested (2026-09-24)
Poles at −n², pole value (n/2)·(X − Σ_{m<n} a(m)/m²) + e_n, pole-only (numerator degree < K, h ≤ K/2). Per K²:
| functional | closeness | consecutive gcd |
|---|---|---|
| ζ(2) tails with exact elementary term (honest functional restricted) | −0.35 to −0.40 | trivial |
| ζ(2) tails WITHOUT the elementary term | −0.08 | trivial |
| G² tails (a = (r₂/4)²), with or without elementary term | −0.002 | trivial |
Closeness is produced by the integral representation, not by the affine algebraic structure. Remove the analytic form and the cancellation vanishes. Coprimality was never the bottleneck. Door 2 (tails + coprimality) is closed.

---
# Quarter-integer poles with the full-rate kernel — open item 1 measured and closed (2026-09-25)

Pipeline re-validated on this machine first (Python 3.12.10, python-flint 0.9.0): ζ(5) control log P_K(ζ5) = −265.1 at K=40 (−3.457/3.291/−0.166), Catalan m=1 log P_K(G) = +608.7 = +0.380/K², midprime height 2.525, ζ(2) half-poles m=0 net −0.130 — all to the digit, 2 s each.

Weight: the full-rate ζ(2)-half-pole weight is w(t) = πt²/(2 sinh²(πt)) = −t² d/dt[1/(e^{2πt}−1)] (moments ∫t^{2e}w = |B_{2e+2}|/2). Binet's second formula gives, for EVERY a>0,
  ∫₀^∞ w(t)/(t²+a²) dt = (a/2)ζ(2,a) − 1/2 − 1/(4a)        (checked to 40 digits at a = 1/4, 3/4, 5/4, 7/4, 77/4, 1/2, 3),
which at a = b/2 is exactly the (b/16)ζ(2,b/2) − 1/(8b) − 1/8 used in zeta2_halfpoles.py. At quarter-integers ζ(2,1/4) = π² + 8G, ζ(2,3/4) = π² − 8G, and ζ(2,a) − ζ(2,a+1) = 1/a², so
  both classes a=(2j+1)/4  ⇒ entries A + π²B + GC  (two unknowns; working variable v=16t², v-poles −(2j+1)², moments 16^e|B_{2e+2}|/2, pole value (b/128)ζ(2,b/4) − 1/32 − 1/(16b));
  one class a=j+1/4 (or j+3/4) ⇒ entries A + ξB with ξ = π² ± 8G (one unknown), spacing 1 at full rate — the ζ(2) machine with 1/2 replaced by 1/4.
Scripts: twounknown.py (exact bivariate det(A+XB+YC) from the triangular grid {a+b≤h} by 2-D forward differences; independent arb determinant at (π²,G) agrees to 1e-6 in log; entries vs direct integrals 1e-8…1e-19), quarter_oneclass.py (K N m OFF).

| functional (N=3) | ρ = rate×spacing | K | closeness | height | net |
|---|---|---|---|---|---|
| two-unknown (π²,G), both classes, m=0 | π | 40 / 60 | −2.072 / −2.119 | 2.530 / 2.637 | +0.458 / +0.519 |
| two-unknown, m=1 | π | 40 | −2.140 | 2.648 | +0.508 |
| Catalan, alternating kernel, m=0 (control, earlier) | π | 40 / 60 | −2.19 | 2.56 / 2.64 | +0.38 / +0.45 |
| one class a=j+1/4, ξ=π²+8G, m=0 | 2π | 40 / 60 | −2.912 / −2.986 | 3.835 / 3.985 | +0.924 / +0.999 |
| one class a=j+3/4, ξ=π²−8G, m=0 | 2π | 40 | −2.912 | 3.819 | +0.907 |
| one class, m=1 (1/4 ; 3/4) | 2π | 40 | −3.044 ; −3.040 | 4.050 ; 4.029 | +1.006 ; +0.989 |
| ζ(2) half-integer, m=0 (control) | 2π | 40 / 60 | −2.905 / −2.983 | 2.775 / 2.878 | −0.130 / −0.104 |
| ζ(2) half-integer, m=1 / m=6 | 2π | 40 | −3.035 / −3.530 | 2.891 / 3.946 | −0.145 / +0.415 |
| ζ(2) one class b≡1 mod 4 (spacing 2), m=0 | 4π | 40 | −3.869 | 4.836 | +0.967 |
| ζ(5) control, m=6 | 2π | 40 | −3.457 | 3.291 | −0.166 |

Two-unknown functional: fails exactly like the alternating Catalan functional — closeness −2.07 (Catalan −2.19), height 2.53 (Catalan 2.53–2.56), net +0.46 → +0.52 from K=40 to 60. The sign carry is not free after all: in v = 16t² the pole set is the same −(2j+1)² as the half-integer machine, and the weight now decays at half the rate per pole spacing (e^{−(π/2)√v} instead of e^{−π√v}) — the same half field the alternating kernel has, reached by doubling the pole density instead of by alternating signs. Structure of P(X,Y): all (h+1)(h+2)/2 monomials nonzero (741 at K=40, 1711 at K=60); the odd-in-G part is as large as the even part (log10 1758.0 vs 1758.2) — nothing lives in (π², G²); |P(π²,−G)|, |P(π²,0)|, |P(0,G)| are all at height size (2.52–2.53 per K²) — the cancellation is specific to +G. Polynomials saved: P_twounknown_K40_N3_m0.txt, _m1.txt, K60_N3_m0.txt.

One-class functional: closeness equals the ζ(2) machine's to three decimals at both K (−2.912 vs −2.905; −2.986 vs −2.983), as ρ=2π predicts — but the height is +1.06 (K=40) / +1.11 (K=60) higher, and the excess is exactly the content-denominator excess: primes in (2K,4K) appear (0.88/K², absent for half-integer poles) and the mid-prime part rises by 0.2. Cause: the pole set {j+1/4} has denominator 4, so its pairwise sums (i+j+1/2) and its Hurwitz partial sums 1/k² (k ≡ 1 mod 4) involve odd integers up to 4K−3, versus 2K−1 for {j+1/2}. Net +0.92 → +1.00, growing with K. Both quarter-integer routes are closed.

Law (13 runs, four different functionals): closeness depends only on ρ = kernel decay rate × pole spacing and on the numerator power m —
  ρ=π:  −2.07…−2.19 (m=0), −2.14/−2.15 (m=1), −2.34 (m=6)      [two-unknown, Catalan]
  ρ=2π: −2.905/−2.912 (m=0), −3.035/−3.044 (m=1), −3.46…−3.53 (m=6)   [ζ(2), one-class quarter, ζ(5), Catalan on one class mod 4]
  ρ=4π: −3.87 (m=0)   [ζ(2) on b≡1 mod 4]
Each doubling of ρ buys ≈ 0.85–0.95 per K². Height depends on the arithmetic of the pole set (v-poles −(2j+1)²: 2.5–2.9; −(4j+1)²: 3.8–4.8) and grows with m. Closing needs closeness + height < 0: the half-integer ζ(2) set with ρ=2π does it (−2.9 + 2.8); every G-affine construction so far has either ρ=π (closeness −2.1 against height 2.5) or ρ=2π on the denominator-4 set (closeness −2.9 against height 3.8).

Why every G-affine construction sits at ρ=π (argument, not a proof): take a kernel with poles on iZ and residue pattern r(n) of period p — p=1 plain (rate 2π), p=2 alternating (rate π), p=4 the χ₋₄ pattern sech(2πt)-type (rate π/2 on iZ). Rational moments need Σ r(n) n^{−2e−2} ∈ Q·π^{2e+2}, which holds for the ζ pattern and the η pattern only (the χ₋₄ pattern gives β(2e+2) — G itself sits in the moments). For those two patterns the pole values are Hurwitz values at (a+j)/p, and G appears iff (a+j)/p ≡ ±1/4 mod 1: using both signs of G the pole spacing is p/2, so ρ = (2π/p)(p/2) = π independent of p; using one sign class the spacing is p (ρ = 2π) but the integer parametrization of the pole set reaches 4K instead of 2K — the height penalty measured above (and already seen in the Catalan b≡1 mod 4 run: closeness −3.53, height 5.97). ζ(2) has no such constraint: spacing 1 at p=1, ρ=2π, integer pairwise sums. That is the whole difference between the machine that closes and the ones that do not. Item 3 of the handoff is thereby sharpened: a G-affine weight with ρ=2π on a lattice-summed pole set has to come from outside the periodic-residue kernel class, or the machine has to change.

---
# Exponent profiles and the K → ∞ limit — open item "asymptotics" (2026-09-25)

Pipeline: profiles.py (one script for the three machines; entries by the shift recurrence u^{s+1}N = (uq_s + c_s)D + (ur_s − c_sD), exact flint determinants, ball-arithmetic evaluation; reproduces zeta5_control, zeta2_halfpoles, catalan_hankel/midprime to the digit; writes profile_*.json), profiles_analyze.py (tables, extrapolations, exact exponent models). Lenses: hybrid base (e_p as base-p digit/carry counts), palindrome (mirror pairs ±a mod p).

## Sweep (per K²; N=3, m=0 unless stated; height = intrinsic + den − num, intrinsic = log max|coeff| of the det polynomial before content removal)
| machine | K | closeness | height | net | intrinsic | den (small / mid) |
|---|---|---|---|---|---|---|
| ζ(2) half-integer | 40 / 60 / 80 / 100 | −2.905 / −2.983 / −3.019 / −3.041 | 2.775 / 2.878 / 2.939 / 2.970 | −0.130 / −0.104 / −0.081 / −0.071 | 0.059 → 0.033 | 2.72 → 2.94 (1.24 / 1.48 → 1.28 / 1.66) |
| ζ(2), N=0 | 40 / 60 / 80 | −3.096 / −3.104 / −3.108 | 2.950 / 2.993 / 3.022 | −0.145 / −0.111 / −0.086 | 0.059 → 0.038 | 2.89 → 2.98 |
| Catalan | 40 / 60 / 80 / 100 | −2.076 / −2.122 / −2.143 / −2.155 | 2.411 / 2.537 / 2.576 / 2.615 | +0.335 / +0.416 / +0.433 / +0.460 | 0.081 → 0.042 | 2.33 → 2.57 (1.03 / 1.30 → 1.06 / 1.51) |
| Catalan, N=0 | 40 / 60 / 80 | −2.185 / −2.189 / −2.192 | 2.564 / 2.641 / 2.650 | +0.379 / +0.452 / +0.458 | 0.082 → 0.050 | 2.48 → 2.60 |
| ζ(5), N=3K/40, m=6 | 40 / 80 | −3.457 / −3.469 | 3.291 / 3.338 | −0.166 / −0.130 | 2.31 → 2.67 | 1.18 → 1.29; num 0.20 → 0.62 |
(K=120 and 160 running for all three at the time of writing; appended below when done.)
N-scan at K=40, m=0: Catalan N=0/3/5/8/10/15/20 → net +0.379/+0.335/+0.328/+0.340/+0.362/+0.458/+0.542 (shallow optimum N/K≈0.1, worth 0.05); ζ(2) N=0/3/5/10 → −0.145/−0.130/−0.108/+0.015 (dropping poles only hurts ζ(2)). The closeness law (ρ, m) holds at fixed N/K; all earlier comparisons had N/K = 0.075.

## Exact exponent laws (N=0; every prime p > K/2; zero residual at K = 40, 60, 80)
  ζ(2):    e_p = 2·#{0 ≤ n < 2K : ⌊n/p⌋ odd} − [p < K] = 2 Σ_{k≥1} (−1)^{k−1} (2K − kp)₊ − [p < K]
  Catalan: e_p = min( 3·#{poles b > p} + 1 , 2(K−1) ) = min( (3(2K−p) − 1)/2 , 2K − 2 )
(N=3 shifts these by a constant −3 resp. −6 on the affected ranges; same limits.)
Reading. ζ(2): each pole costs p⁴ when its leading base-p digit ⌊b/p⌋ is odd and nothing when it is even — the poles in an even block are the mirror images (2kp − b) of poles in the preceding odd block and cancel p-adically: the palindrome ±a mod p, i.e. the "Hermite interpolation at mirror pairs" mechanism of the ζ(5) proof, seen directly in the content. Catalan: the alternating kernel destroys the mirror cancellation — every pole above p costs p³ (not p⁴), nothing cancels, and the exponent saturates at the trivial bound 2 per pole (each row carries at most p² in its denominators; one p² leaves with the content). Scaling limits e_p ≈ K·φ(p/K):
  φ_ζ2(u) = 2 Σ_k (−1)^{k−1} (2 − ku)₊   (dips at u = 1, 2/3, …; equals 2u on (2/3,1), 2(2−u) on (1,2)),  ∫₀² φ_ζ2 = 4 ln 2 = 2.7726
  φ_cat(u) = min( 3(2−u)/2 , 2 )          (plateau 2 up to u = 2/3, one ramp to 0 at u = 2),           ∫₀² φ_cat = 8/3 = 2.6667
So Σ_{p} e_p log p / K² → ∫φ by the prime number theorem, from below (θ(x)/x < 1): the mid-prime denominators keep rising with K in both machines — this is why the finite-K nets drift upward. Catalan's arithmetic is slightly cheaper than ζ(2)'s (8/3 < 4 ln 2); the deficit is not arithmetic.

## The 2-adic K² term (the base-2 digit of the poles has a price)
The ramp models describe p ≥ 3. The prime 2 is a separate K²-scale term: e₂ = 1004, 2121, 3646 (ζ(2), K=40/60/80) → e₂/K² = 0.628, 0.589, 0.570 = 0.51 + 4.6/K, i.e. 0.51·ln2 ≈ 0.36 per K² of height forever; Catalan e₂ = 393, 782, 1337 → e₂/K² → ≈ 0.17, i.e. ≈ 0.12 per K². (Half-integer poles b/2: every residue denominator ∏(b'−b)(b'+b) is divisible by 2³ per pair, 3K² in total, of which 0.5K² survives the determinant.) p = 3 leaves ≈ 0.05 (ζ(2)) / 0.02 (Catalan); p ≥ 5 are O(K log K) and vanish per K². Integer poles (ζ(5)) have no such term — there 2 and 3 sit in the numerator content (0.20 → 0.62 per K²).

## Asymptotic estimates (exact profile + 2-adic term for the height; closeness from the sequences)
| machine | closeness_∞ | height_∞ = ∫φ + δ₂ + δ₃ | net_∞ |
|---|---|---|---|
| ζ(2) half-integer, m=0 | −3.12 (N=0: −3.096, −3.104, −3.108) | 2.773 + 0.36 + 0.05 ≈ 3.18 | ≈ +0.05 ± 0.1 |
| Catalan, m=0 | −2.20 (N=0: −2.185, −2.189, −2.192) | 2.667 + 0.12 + 0.02 ≈ 2.81 | ≈ +0.6 |
| ζ(5), m=6 (Fauzan) | −3.47 | ≈ 3.4 (paper: A+U = −0.017) | ≈ −0.02 |
Direct extrapolations of net(K) agree in sign: ζ(2) −0.03 / −0.06 / −0.11 (N=3; 1/K, 1/K², logK/K) and −0.01 / +0.02 / +0.07 (N=0); Catalan +0.57 / +0.69 / +0.85 (N=3), +0.48 / +0.36 / +0.19 (N=0).
Conclusions. (1) The Catalan gap is structural and grows with K: +0.335 (K=40) → +0.460 (K=100) → ≈ +0.6. (2) The "decisive control" needs a correction: the ζ(2) half-integer machine closes at every K computed (−0.13 → −0.07) but its margin erodes like the PNT deficit and the asymptotic net is ≈ 0 ± 0.1 — marginal, not a clean close; only the ζ(5) machine (integer poles, m=6, numerator content) has a genuine negative limit. (3) Asymptotically Catalan has the LOWEST height of the three (2.8 vs 3.2 vs 3.4); the whole deficit is closeness: −2.2 (ρ=π) against −3.1 / −3.5 (ρ=2π). Restoring the mirror cancellation in the alternating kernel would save at most ∫(φ_cat − φ_ζ2-like) ≈ 0.1–0.3 per K² — not a route.

## Addendum (2026-09-25, later): K=120, N=0 at K=100, the integer-pole control, the per-pole cost law
Sweep additions (per K²). ζ(5) K=120 (N=9, m=6): −3.4729 / 3.3476 / −0.1253; K=160 (N=12, m=6, 149 determinants, 3.8 h): −3.4749 / 3.3480 / −0.1269 — the sequence −0.166, −0.130, −0.125, −0.127 has converged at ≈ −0.126, NOT at the paper's −0.017: the exact primitive part beats the m_{K,M}·S_K normalizer by ≈ 0.11 per K² (num content 0.20 → 0.62 → 0.91 → 1.11; intrinsic 2.31 → 2.67 → 2.90 → 3.07; den 1.18 → 1.29 → 1.35 → 1.39; height 3.29 → 3.34 → 3.348 → 3.348; closeness → −3.476). The per-pole law e_p = 7(K−p) − (6N−1) holds at K=160 for 0.67 < p/K < 0.87 (p = 107…137 exact; deviations grow below u ≈ 0.65 as at K=80). ζ(2) half-poles K=120 (N=3): −3.0546 / 2.9921 / −0.0625 (nets −0.130, −0.104, −0.081, −0.071, −0.063 at K=40…120). Catalan K=120 (N=3): −2.1624 / 2.6326 / +0.4702 (nets +0.335, +0.416, +0.433, +0.460, +0.470). N=0 at K=100: ζ(2) −3.1102 / 3.0375 / −0.0727; Catalan −2.1927 / 2.6745 / +0.4817.

Integer-pole control (machine zeta2int in profiles.py: the SAME full-rate weight πt²/(2sinh²πt) on integer poles u=−j², j=1..K; pole value (j/2)(ζ(2) − H2_{j−1}) − 1/2 − 1/(4j); moments |B_{2e+2}|/2; ρ=2π):
| K (N=0) | closeness | height | net | den small / mid |
|---|---|---|---|---|
| 40 | −3.127 | 2.274 | −0.852 | 1.11 / 1.09 |
| 60 | −3.125 | 2.289 | −0.836 | 0.99 / 1.25 |
| 80 | −3.124 | 2.287 | −0.836 | — |
Closeness identical to the half-integer machine (−3.096, −3.104, −3.108): the ρ-law does not see the arithmetic of the pole set. Height 0.68–0.70 lower: the half-integer pole set costs ≈ 0.7 per K² at equal closeness (its Hurwitz partial sums 1/k² reach k = 2K instead of K, and it pays the 2-adic K² term). Relative to this best level-2 machine the Catalan deficit at K=40, N=0 is 1.23 = 0.94 closeness (ρ=π vs 2π) + 0.29 height (half-integer poles, partly offset by the alternating kernel's cheaper pairs).

Per-pole p-adic cost law (top range of p, exact): e_p = c·#{poles above p} + O(1) with
  c = s + 2 for plain kernels, s = level: ζ(5) integer poles c = 7 (e_p = 7(K−p) − 2 at m=1; 7(K−p) − (6N−1) at m=6, K=40 and 80); ζ(2) half-integer poles c = 4; ζ(2) integer poles above K: pairs only, e_p = 2K − p − 2;
  c = s + 1 for the alternating kernel: Catalan c = 3.
The level enters the exponent linearly — each rung of the ladder costs one more power of p per pole (hyperoperation lens: level ↔ p-adic multiplicity); the pair term is 2 for plain and 1 for alternating kernels; the mirror cancellation (poles with even leading base-p digit drop out) exists only for plain kernels.

What the numerator power m does (K=40, N=3; height = intrinsic + den − num): Catalan m=6: 3.10 + 1.22 − 0.99 = 3.33 (m=0: 0.08 + 2.33 − 0 = 2.41); ζ(2) m=6: 3.04 + 1.83 − 0.92 = 3.95; ζ(5) m=1 (no numerator): 0.29 + 2.64 − 0 = 2.93, closeness −2.94, net −0.017 (m=6: 2.31 + 1.18 − 0.20 = 3.29, closeness −3.46, net −0.166). E_N^m removes the small-prime denominators (small part 0.06 at m=6) and creates ≈ 1 per K² of numerator content, but inflates the coefficient size by ≈ 3: a good trade at integer poles (ζ(5): +0.15 net), a bad one at half-integer poles.

Convergence. Σ_p e_p log p / K² approaches ∫φ as θ(x)/x → 1, i.e. from below with an O(K^{−1/2}) deficit (10% at K=100: model sum 2.45 vs 2.77 for ζ(2)); the multi-digit small-prime terms (Legendre-type, p² ≤ 2K) decay from above like K^{−1/2}; the 2-adic term is a genuine constant (e₂/K² → 0.51 for ζ(2), ≈ 0.18 for Catalan). So 1/K extrapolations of height and net are unreliable, and the profile-based limits are the ones to quote: ζ(2) half-poles den → 4 ln 2 + 0.36 ≈ 3.13, closeness → −3.112, net → 0.00 ± 0.05 (marginal); Catalan den → 8/3 + 0.12 ≈ 2.79, closeness → −2.197, net → ≈ +0.6; ζ(5) → −0.12; integer-pole ζ(2) → ≈ −0.7. K=160 (158 determinants, 3.8–3.9 h each, added 2026-09-26): ζ(2) half-poles N=3: −3.0716 / 3.0235 / −0.0481 (nets −0.130, −0.104, −0.081, −0.071, −0.063, −0.048 at K = 40…160 — still rising toward 0; height 2.775 → 3.024 toward the predicted 3.13; the K=160 exponent profile shows the full φ_ζ2 shape: 2u on (2/3,1) with the dip to 1.33 at p=107, 2(2−u) on (1,2), e₂/K² = 0.519). Catalan N=3: −2.1717 / 2.6784 / +0.5067 (nets +0.335, +0.416, +0.433, +0.460, +0.470, +0.507; height 2.411 → 2.678 toward the predicted 2.79; closeness → −2.20). zeta2int K=80: −3.1235 / 2.2873 / −0.8362 (flat). Direct 1/K extrapolations from the last two points now give ζ(2) −0.005 and Catalan +0.62, in line with the profile-based limits 0.00 ± 0.05 and ≈ +0.6. All sweep data are in the profile_*.json files (37 runs).

---
# Theory of the p-adic laws — open item 2 (2026-09-26)

Frame. Change basis from the monomials v^i to the Lagrange basis L_b(v) = D(v)/(v+b²) at the poles (P_K is basis-invariant up to content; the change of basis has det ±∏_{b<b'}(b'²−b²) = ±det V, V the Vandermonde in the nodes −b²). Then det M = det G / det(V)² with the Lagrange Gram matrix
  G_{bb'} = μ(D/((v+b²)(v+b'²)))   (b ≠ b': a polynomial of degree K−2 — moments only, no pole values),
  G_{bb}  = μ(S_b) + D'(−b²)·w_b,   w_b = A_b + X·B_b   (the pole value enters only on the diagonal, multiplied by the residue product D'(−b²) = ∏_{b'≠b}(b'²−b²)),
so  e_p = −v_p(content det M) = 2·v_p(det V) − v_p(det G).

Top range K < p < 2K, half-integer poles, N=0; n = #{poles b > p} = (2K−p−1)/2.
 (i) Vandermonde: p | b'²−b² only through b+b' = 2p (differences are < p); each pole above p has exactly one mirror partner 2p−b below p ⇒ v_p(det V) = n ⇒ 2n.
 (ii) Pole values: for b > p the partial sum (Σ_{odd k<b} 1/k² or its alternating version) contains 1/p² and D'(−b²) contains p once (the partner) ⇒ v_p(D'(−b²)w_b) = −1; for b = p, w_p has 1/(8b) (resp. 1/(4b)) and D'(−p²) has no p ⇒ −1. The diagonal carries n+1 entries of valuation −1.
 (iii) Moments — von Staudt–Clausen, the anchor at the negative rungs (Run 10): p | denom(B_{2e+2}) iff (p−1) | (2e+2). The Lagrange entries use degrees ≤ K−2, so for p ∈ (K,2K) exactly one moment, e = (p−3)/2, has v_p = −1. Its coefficient matrix N_{bb'} = [v^{(p−3)/2}] D/((v+b²)(v+b'²)) is, in the monomial basis, the shifted Hankel matrix (a_{i+j−K−m}) of the expansion 1/D = Σ a_n v^{−K−n}: anti-triangular with 1's, rank K−1−m = K−(p−1)/2 = n+1.
 Plain kernel (ζ(2)): the 1/p-part of G is c·N + diag δ, generic rank 2(n+1) ⇒ v_p(det G) = −2(n+1) ⇒ e_p = 2n + 2n + 2 = 4n + 2 = 2(2K−p).  (exact at K = 40, 60, 80, 160)
 Alternating kernel (Catalan): the moments 4^e(2^{2e+2}−1)|B_{2e+2}| are Genocchi numbers up to powers of 2 — (2^{p−1}−1)B_{p−1} is p-integral by Fermat (checked e ≤ 119) — so c = 0 and the rank is n+1 ⇒ e_p = 2n + n + 1 = 3n + 1 = (3(2K−p)−1)/2.
 Hybrid confirmation (padic_test.py, K = 40 and 60; positivity is irrelevant for the content): ζ(2) pole values with Genocchi moments give the Catalan exponents EXACTLY for every p ≥ 2K/3, plateau 2K−2 included; Catalan pole values with von Staudt moments give 4n+2 exactly for p > K and a flat 2K−1 below K. Above K the moments alone decide the law; the sign pattern of the pole values is irrelevant.

Middle range 2K/3 < p < K, alternating kernel: the Lagrange accounting stays exact with the extra pairs b+b' = 4p and b'−b = 2p: v_p(det V) = (p−1)/2 + 2(K−p); a pole b > p costs 1 iff it has exactly one partner in range (b < 4p−2K), 0 if two (v_p(D') = 2 absorbs the 1/p²), plus 1 for b = p when 3p > 2K−1. Total (3(2K−p)−1)/2 again — checked by hand at K=40: p=31: 7 + 66 = 73, p=37: 16 + 48 = 64. Below 2K/3 the count overshoots (100 vs 78 at p=23) and the trivial cap 2(K−1) (one p² per row, one p² out with the content) takes over.
Middle range, plain kernel: e_p = 2p−1 on (2K/3, K), far below the Lagrange count — det G is divisible by p⁵ at K=40, p=31 — and this cancellation needs BOTH plain pole values AND von Staudt moments (each hybrid loses it). This is the Hermite-basis mechanism of the ζ(5) proof (rows of mirror poles ±a mod p become congruent after the divided-difference change of basis); its net effect is the digit rule "4 per pole with odd leading base-p digit". Not derived here.

ζ(5), integer poles: the moment factor (2e+3)(2e+4)(2e+5)/24 = (n+1)(n+2)(n+3)/24 at n = 2e+2 = k(p−1) is ≡ (1−k)(2−k)(3−k) mod p, so the von Staudt p is cancelled for k = 1, 2, 3 (checked: Fauzan's moments for e ≤ 119 contain no prime > 61, the first p with 4(p−1) ≤ 240). For p > K/2 the ζ(5) moments are therefore p-integral — Fauzan's weight is Genocchi-like exactly where it matters — and the per-pole cost 7 = (5 − 2 partners) + 2·2 pairs has no von Staudt term.

Integer poles, top range K < p < 2K (padic_test.py POLES=int, N=0, K = 40 and 60; all mirror pairs (a, p−a) lie below p and have p-integral pole values):
  von Staudt moments (= zeta2int):  e_p = 2·v_p(det V)  exactly   (the earlier remark "2K−p−2" came from an N=3 run; with N=0 it is 2K−p+1 = 2·#pairs);
  Genocchi moments:                  e_p = v_p(det V)    exactly   (one p per pair: the Lagrange rows of a mirror pair with p-integral pole values are congruent mod p, so det G gains one p per pair).
Universal statement for K < p < 2K, both geometries: a mirror pair costs 2 (Vandermonde) − 1 (Hermite congruence, only if both members have p-integral pole values) + 1 (if one member carries the level's 1/p^s) + 1 (von Staudt moment B_{p−1}, plain kernels only):
  half-integers, Genocchi (Catalan): 2 + 1 = 3 per pair (all pairs straddle p) → 3n (+1 for b=p);   half-integers, von Staudt (ζ(2)): 4 per pair → 4n (+2);
  integers, Genocchi: 2 − 1 = 1 per pair;   integers, von Staudt (zeta2int): 2 per pair.
So the von Staudt moment costs exactly one p per mirror pair everywhere, which is what the (2e+3)(2e+4)(2e+5) factor of Fauzan's weight and the Genocchi structure of the alternating kernel avoid.
Per-pole law restated: c = (level − #partners) + 2·(pairs per pole) + [von Staudt]: half-integers 2−1+2+1 = 4 (ζ(2)) and 2−1+2+0 = 3 (Catalan); integers at level 5: 5−2+4+0 = 7.
Not derived: the digit rule below K (Hermite bases) and the 2-adic K² term (e₂ ≈ 0.51K² for ζ(2), 0.18K² for Catalan): for odd b, b' every (b'−b)(b'+b) is divisible by 8, so det V² alone carries ≈ 4K² powers of 2 and the mod-2 Hermite cancellation (all poles congruent) removes most but not all.

# Consequence for open item 1 — floor and wall
The Catalan machine sits at BOTH limits of this machine class. Arithmetic floor: per pole above p it pays 2 (mirror pairs — unavoidable for any K odd numbers below 2K) + 1 (the 1/p² of a level-2 Hurwitz partial sum — unavoidable for G-affine values) and nothing for the moments (Genocchi); 3 is the minimum and the alternating kernel achieves it, which is why its ∫φ = 8/3 is below ζ(2)'s 4 ln 2. Analytic wall: ρ = π is forced for two-class G-affine values from any lattice kernel with rational moments (rational moments ⇔ ζ or η residue pattern ⇔ G reachable only through the base-2/base-4 refinement of the shift, which halves rate × spacing). Closeness −2.2 against height ≥ 8/3 + 2-adic ≈ 2.8: the class is exhausted for G by ≈ 0.6 per K². Exits are outside the class: a weight with non-periodic residues and rational moments (none known), or a different criterion (Apéry-style partner sequence for the Ramanujan series, Run 3 item 4).

---
# The ledger (2026-09-26): every p-adic contribution named, in base p, summed against the measurement
Tool: ledger.py MACHINE K [p]. For each prime it lists the poles' base-p digits (leading digit d₁ = ⌊b/p⌋), the complement pairs b+b' ≡ 0 mod p^k (a carry with zero result digit — the Midy / p-complement pairs of Run 1) and the coincidence pairs b' ≡ b mod p^k, v_p(det V) = Σ_pairs [v_p(b'+b) + v_p(b'−b)] (checked against the exact product), the diagonal valuation δ_b = v_p(D'(−b²)) + min(v_p A_b, v_p B_b) of every pole, the von Staudt moments with their naive ranks, the Lagrange prediction naive = 2v_p(det V) + pole cost + Σ ranks, the measured e_p, and residual = measured − naive (the Hermite / Midy cancellation).

## Ledger, K = 40, N = 0 (columns: 2v(det V) | pole cost | von Staudt ranks | naive | measured | residual)
| p | u | ζ(2) half | Catalan | integer poles (zeta2int) |
|---|---|---|---|---|
| 23 | 0.575 | 100 + 0 + 54 = 154 → 67 (−87) | 100 + 0 = 100 → 78 (−22, cap 2K−2) | 102 + 1 + 54 = 157 → 79 (−78) |
| 29 | 0.725 | 72 + 4 + 38 = 114 → 57 (−57) | 72 + 4 = 76 → 76 (0) | 72 + 1 + 38 = 111 → 73 (−38 = −ranks) |
| 31 | 0.775 | 66 + 7 + 35 = 108 → 61 (−47) | 66 + 7 = 73 → 73 (0) | 66 + 1 + 35 = 102 → 67 (−35 = −ranks) |
| 37 | 0.925 | 48 + 16 + 26 = 90 → 73 (−17) | 48 + 16 = 64 → 64 (0) | 48 + 1 + 26 = 75 → 49 (−26 = −ranks) |
| 41 | 1.025 | 38 + 20 + 20 = 78 → 78 (0) | 38 + 20 = 58 → 58 (0) | 40 + 0 + 20 = 60 → 40 (−20 = −rank) |
| 43…79 | >1 | residual 0 throughout | residual 0 throughout | residual = −rank throughout (e_p = 2v_p(det V)) |
Reading: the Catalan ledger is exact for u ≥ 0.725 (no cancellation at all); the integer-pole ledger is exact once the von Staudt term is deleted (the pairing absorbs it completely, e_p = 2v_p(det V) + [b=p]); the plain half-integer ledger is exact for p > K and shows the Hermite cancellation below K, which removes the von Staudt ranks, the pole costs and a little more (p=31: 35 + 7 + 5).
Per-pole ledger at p = 31 (both half-integer machines, identical rows): poles with d₁ = 0 have a complement partner (2p−b, d₁ = 1) and, for b ≤ 17, a coincidence partner (b+2p, d₁ = 2): δ = 1 or 2, cost 0. Poles with d₁ = 1: b = 31 (δ = −1, from 1/(8b)), b = 33…43 (one complement partner below, δ = 1 − 2 = −1, cost 1 each), b = 45…61 (partners 17…1 below and 79…63 above, δ = 0). Poles with d₁ = 2 (63…79): complement below and coincidence below, δ = 0. Total pole cost 7 = the seven poles whose base-31 digits are (1, d₀ ≤ 12). The von Staudt moments for p=31 are e = 14 (2e+2 = 30) and e = 29 (2e+2 = 60).

## The rules on the doubled lattice
Half-integer poles a = b/2 are the odd points of the lattice n = 2a ∈ [0, 2K); the base-p leading digit d₁(n) = ⌊n/p⌋ is the block index. Exactly (p > K/2, all K tested):
  ζ(2):    e_p = 2 · #{n < 2K : d₁(n) odd} − [p < K]          (every lattice point in an odd block, odd or even n, is charged 2)
  Catalan: e_p = (3 · #{n < 2K : d₁(n) ≥ 1} − 1)/2, capped at 2K−2   (every lattice point outside block 0 is charged 3/2)
For p > K both counts are 2K − p and the charges 2 vs 3/2 are the von Staudt (+1/2 per point) difference. Below K the plain kernel's even blocks are cancelled against the odd ones by the complement map n ↦ 2kp − n (block 2k−1 ↔ block 2k−... the Midy halving), the alternating kernel's are not.

## Base-2 ledger (all poles share the low digit 1)
Pairs by depth k = v₂(b'²−b²): n_k = K²/2^{k−1} (K=40: 400, 200, 100, 48, 24, 8 for k = 3…8), so 2v₂(det V) = 4K²(1 + O(1/K)) (3.78, 3.84, 3.88, 3.90 per K² at K = 40…100). Pole cost 0 (the 2-adic denominators of the pole values are absorbed by D'). What survives the mod-2 Hermite cancellation, per pair: ζ(2) e₂/#pairs = 1.287, 1.198, 1.154, 1.124 → 1 (e₂ ≈ K(K−1)/2 + K log₂K); Catalan 0.504, 0.442, 0.423, 0.404 → 1/3 (1/K fit 0.328); integer poles e₂ ≈ 2.7·K log K → 0 per K². Hence the 2-adic constants are exact: ½ ln 2 = 0.3466 (ζ(2) half), ⅙ ln 2 = 0.1155 (Catalan), 0 (integer poles). The alternating sign (−1)^j is the base-2 digit d₁ of b (b ≡ 1 or 3 mod 4), which is why the alternating kernel keeps only a third of a power of 2 per pair.

## Asymptotic bookkeeping, closed form where the ledger allows it
| machine | arithmetic side A (exact) | analytic side U (fitted) | A + U |
|---|---|---|---|
| ζ(2) half-integer, m=0 | 4 ln 2 + ½ ln 2 = (9/2) ln 2 = 3.1192 | −3.1192 ± 0.002 (N=0: −3.0955, −3.1038, −3.1078, −3.1102; 1/K² fit −3.1192) | 0.000 ± 0.003 |
| Catalan, m=0 | 8/3 + ⅙ ln 2 = 2.7822 | −2.1973 ± 0.002 (N=0 fit; −ln 9 = −2.1972) | +0.585 |
| integer-pole ζ(2), m=0 | ≈ 2.3 (ledger: 2v_p(det V)-driven, no 2-adic term) | −3.12 | ≈ −0.8 |
| ζ(5), Fauzan | 3.348 (intrinsic 3.07 + den 1.39 − num 1.11) | −3.476 | −0.126 |
The half-integer ζ(2) machine is critical to four decimals: A = U = (9/2) ln 2 within the fit error. [Both identities are now proved by the potential theory below: U(2π) = −(9/2) ln 2 and U(π) = −ln 9 exactly.] Sub-leading term: log P_K(ζ(2)) = −232.4, −400.7, −549.6, −727.0 at K = 40, 60, 80, 100 (N=0), i.e. −1.575, −1.631, −1.568, −1.579 per K ln K — consistent with ≈ −1.58·K ln K but the θ(x) − x fluctuations of the mid-prime content (±20 in log P at these K) defeat any closed-form identification (π/2 = 1.5708 is within the scatter, not confirmed); a 3-parameter fit is unstable. Left as an observation.

---
# The analytic side in closed form — U(ρ) from potential theory (2026-09-26)
Brown (arXiv 2604.20741, Prop. 3.5) is Heine's identity: det Q_N = (1/N!)∫(det V)² ω^{⊠N}. For our machines the Hankel matrix at X = ξ is the moment matrix of the positive measure dμ = w(t)/D(4t²) dt in v = 4t², so det M_K(ξ) is a β=2 log-gas partition function. With x = t/K, K poles a_j = s(j+½) (spacing s), w ~ t²e^{−rate·t}:
  log det M_K(ξ) = −K² · min_ν E[ν] + o(K²),   E[ν] = ∬ [−log|x−y| − log(x+y)] dν dν + ∫ [rate·x + (1/s)∫₀^s log(x²+y²) dy] dν,   ν a probability measure on [0,∞).
The substitution x → x/s leaves E invariant up to nothing: E depends on ρ = rate × s only — the ρ-law is a theorem. Symmetrising ν to [−R,R] gives the classical weighted log-energy 2[I(ν̃) + ∫Q dν̃], Q(x) = (ρ/2)|x| + ½∫₀¹ log(x²+y²) dy, solved by Chebyshev inversion of the log kernel (equilibrium.py). With τ(y) = (√(y²+R²) − y)/R the Chebyshev coefficients of the pole field are exact, q_{2m} ∝ (−1)^m ∫₀¹ τ^{2m} dy, and the soft-edge condition Σ_k k q_k = 2 becomes  (n+1)R − 1 = √(1+R²),  n = ρ/π:
  R(n) = 2(n+1)/(n(n+2))   (4/3, 3/2, 8/5, 5/3, 7/4, 9/5 for n = 1, 2, 3, 4, 6, 8 — the numerically observed values).
The energy E = 2q₀ − 2log(R/2) − ¼Σ_k k q_k² reduces to three explicit integrals (closedform.py, 30 digits); PSLQ identifies E(1) = 2 ln 3, E(2) = (9/2) ln 2, E(3) = 3 ln 5 − ln 3, E(4) = (7/2) ln 3 + ½ ln 2, E(5) = 4 ln 7 − 2 ln 5, E(6) = 11 ln 2 − (5/2) ln 3, i.e.
  **E(n) = ((n+3)/2) ln(n+2) − ((n−1)/2) ln n,     U(ρ) = closeness_∞ = −E(ρ/π) = ((ρ/π−1)/2) ln(ρ/π) − ((ρ/π+3)/2) ln(ρ/π+2).**
Checks: U(π) = −ln 9 (Catalan), U(2π) = −(9/2) ln 2 (ζ(2) half-poles, integer poles, two-unknown, Catalan at spacing 2), U(4π) = −½ln 2 − (7/2) ln 3 = −4.1917, and U → −1 − 2 ln(ρ/π) for ρ → ∞ (the pure-|x| gas 3 + 2 ln n minus the pole field's value −2 at the origin). Each doubling of ρ buys ≈ 0.9, as the law said.
Sweep confirmation (N=0, m=0, K = 40/60/80; the finite-K offset is +0.01 to +0.03 as for ρ = 2π): ζ(2) spacing 2 (ρ=4π): −4.149, −4.164, −4.171 → −4.1917; spacing 3 (6π): −4.827, −4.844, −4.853 → −4.8781; spacing 4 (8π): −5.330, −5.349, −5.359 → −5.3862; Catalan spacing 2 (2π): −3.093, −3.102, −3.107 → −3.1192; spacing 3 (3π): −3.697, −3.709, −3.714 → −3.7297; spacing 4 (4π): −4.155, −4.168, −4.174 → −4.1917. The alternating kernel at spacing 2 lands on the plain kernel's constant, as the scaling demands. K=100 points (sweep complete): ζ(2) spacing 2/3/4: −4.1751 / −4.8582 / −5.3644 against −4.1917 / −4.8781 / −5.3862; Catalan spacing 2/3/4: −3.1092 / −3.7173 / −4.1779 against −3.1192 / −3.7297 / −4.1917 — every offset is +0.010 to +0.022 and shrinking like 1/K.
Dropped poles (N ∝ K, no numerator): h = K−N particles with the tail poles on (α/(1−α), 1/(1−α)), α = N/K, so closeness = −(1−α)²·E'. At α = 0.075: predicted −2.923 (ζ(2)) and −2.083 (Catalan) against measured −2.905 and −2.076 at K=40 (N=3). With a numerator power m the head field −m∫ log(x²+y²)dy repels the gas from the origin, the support opens a gap around 0, and the single-interval solver is invalid (its m ≥ 1 numbers are wrong); the gap problem is the content of Fauzan's inequality (6.14) and is treated next.
Bookkeeping consequence: for every m=0 machine the whole asymptotic ledger is now closed-form: net_∞ = A(lattice) + U(ρ) — ζ(2) half-poles (9/2)ln 2 − (9/2)ln 2 = 0 exactly; Catalan 8/3 + ⅙ ln 2 − ln 9 = +0.585; integer-pole ζ(2) ≈ 2.3 − (9/2) ln 2 ≈ −0.8. For G at ρ = π the analytic side is fixed at −ln 9 = −2.197 by the theorem, so no G-affine construction in the class can beat it, and the arithmetic floor 8/3 + ⅙ ln 2 = 2.782 exceeds it by 0.585.

## Numerator shapes (m ≥ 1, N ∝ K): the raw determinant is predicted, the coefficient size is not
gapgas.py solves the same gas in z = (t/K)² with kernel −log|z−z'| and field W(z) = ρ√z + ∫_{t1}^{t2} log(z+y²)dy − m∫₀^{t1} log(z+y²)dy, t1 = α/(1−α), t2 = 1/(1−α); near 0, W = W(0) + (ρ − mπ)√z + …, so the cusp is attractive (hard edge, density ∝ z^{−1/2}) for ρ > mπ and repulsive for ρ < mπ, where the support opens a gap [za, zb] with two soft edges. The gap is tiny (za ≈ 8·10⁻⁵ for Fauzan's shape) and changes E by < 10⁻⁴; the m ≥ 1 energies are E'(2π, 0.075, 5) = 4.6568 (ζ(5)), E'(π, 0.075, 6) = 3.4098 (Catalan m=6), E'(2π, 0.075, 6) = 4.8679, E'(2π, 0.075, 1) = 3.6883, E'(π, 0.075, 1) = 2.6355.
What they predict is the RAW determinant, up to the trivial normalisation of the numerator E_head(v)^m ≈ (4K²)^{mN} at v ≈ 4K²x²:
  log det M(ξ)/K² = h·m·N·log(4K²)/K² − (1−α)²E'      (integer poles: log(K²) and m−1 copies).
| machine | K | prediction | measured closeness + intrinsic |
|---|---|---|---|
| ζ(2) half, N=3, m=6 | 40 | 3.648 − 4.165 = −0.517 | −3.5305 + 3.0378 = −0.493 |
| ζ(5) Fauzan (N=3K/40, m=6) | 40 / 80 / 120 / 160 | −1.43 / −0.94 / −0.66 / −0.46 | −1.148 / −0.802 / −0.574 / −0.403 |
(the ζ(5) differences 0.28, 0.14, 0.09, 0.06 shrink like 1/K). So Heine + potential theory account for det M(ξ) in every machine tested, numerator or not, and the K² log K drift of the ζ(5) raw determinant is just the normalisation.
The closeness, however, subtracts log max|coeff| of det(A + XB). For m = 0 that is o(K²) (measured intrinsic 0.059 → 0.024 ∝ K log K/K²): the polynomial's size is set by its leading coefficient ∏_b B_b = e^{O(K log K)}, so closeness = −E exactly. For m ≥ 1 with N ∝ K the max coefficient is a K² quantity below the normalisation (ζ(2) m=6: 3.04 = 3.65 − 0.61; leading coefficient alone would be 3.65 − 1.32) and the closeness is −(1−α)²E' + C_m with C_m ≈ 0.6 (ζ(2), Catalan m=6), 0.51 (ζ(5): −3.98 + 0.51 = −3.47 ✓). C_m is a two-species gas: the coefficient of X^k is the mixed partition function of h−k continuous particles on v > 0 and k discrete particles at the poles v = −b² (Heine for the signed measure μ + (X−ξ)ν_B). Not solved here. For the criterion itself only net = log det M(ξ) − log content matters, so the coefficient size is a bookkeeping quantity, not an obstacle; but a closed-form net for Fauzan-type shapes needs the numerator's content (the ledger with numerator) as well.
Dropped poles without numerator (m = 0, N ∝ K) need no correction: closeness = −(1−α)²E'(ρ, α, 0); at α = 0.075: −2.923 (ζ(2)) and −2.083 (Catalan) vs −2.905 and −2.076 at K=40.

---
# Working backwards from the miss (2026-09-26): the budget, the restricted search, the target
Budget at ρ = π (fixed by the theorem): A < ln 9 = 2.1972. Catalan spends A = 8/3 + ⅙ ln 2 = 2.7822, itemised per K²:
| line | cost | origin |
|---|---|---|
| plateau u < 2/3 | 1.333 | trivial cap: one p² per row (every pole above p carries 1/p²) |
| top-range mirror pairs | 0.889 | Vandermonde, 2 per pole above p (the complement b ↦ 2p−b) |
| top-range partial sums | 0.444 | the 1/p² of the Hurwitz tail Σ_{k<b} 1/k² above p |
| 2-adic tax | 0.116 | all poles share the low base-2 digit; ⅓ of a power of 2 survives per pair |
Miss = 0.585. Sub-budgets that would close it: partial sums + 2-adic = 0.560 (0.025 short); adding the ζ(2)-type digit cancellation on the plateau (worth 1.333 − 1.217 = 0.117) closes it; pairs alone (0.889) would also close it but the complement pairs are unavoidable for any K odd numerators below 2K.

Restricted search inside the class — the only free knob left is the pole subset S (same kernel, m=0). Model: for a subset of density 1/ℓ on extent ℓ (in units of K), U = −E(ℓ) (potential theory sees only the density) and, from the per-pole rule "2 without a mirror partner, 3 with one" plus the cap 2, A_random(ℓ) = 4ℓ(1+ℓ)/(2ℓ+1), A_mirror-closed(ℓ) = 8ℓ/3 (both = 8/3 at ℓ=1). d(net)/dℓ at ℓ=1 is +2.22 − 1.22 = +1.0: thinning always loses. Exact machine (profiles.py with an explicit subset, K=40):
| design | closeness | height | net |
|---|---|---|---|
| full half-integer lattice | −2.185 | 2.564 | +0.379 |
| random half density, extent 4K (partners w.p. ½) | −3.053 | 4.551 | +1.498 |
| one class mod 4, extent 4K (mirror-closed) | −3.093 | 4.923 | +1.831 |
| 1 mod 4 below 80 then 3 mod 4 to 160 (no straddling pairs) | −3.108 | 4.799 | +1.690 |
| dense core (odd b<60) + every other odd b to 140, K=50 | −2.318 | 3.055 | +0.737 (full K=50: +0.398) |
The partner rule is visible (random 4.55 < closed 4.92) and every thinning loses more height than it gains closeness. Together with the N-scan (optimum N/K ≈ 0.1, worth 0.05) and the m-results, the within-class minimum is the plain Catalan machine at ≈ +0.585 (finite K: +0.33 to +0.51).

Target specification (what an outside construction must have), read off the ledger:
1. ρ = π is fine IF the arithmetic drops below ln 9; the closeness need not improve.
2. Delete the partial-sum line: G-affine pole values whose "tails" are p-integral above p — i.e. denominators governed by carries (Kummer/Legendre type, as in binomial sums Σ 1/((2n+1)²C(2n,n)) where v_p(C(2n,n)) = number of carries of n+n in base p) instead of by lcm(1..2K)² (every prime squared, as in every Hurwitz tail Σ_{k<b} 1/k²). This is the hybrid-base lens as a design rule: cost = carries, not digits.
3. Delete the 2-adic tax: poles not all in one base-2 class (integer poles pay nothing), which for G means leaving the half-integer lattice — only possible outside lattice kernels.
4. Keep Genocchi-type (von Staudt-free) moments, and the plateau's digit cancellation is a bonus worth 0.117.
The Ramanujan/Apéry-type representation of G is the natural carrier of 2 (its tails are binomial), but it brings π ln(2+√3) as a second unknown (results Run 3 item 4) — the "partner sequence" problem. That is the point from which to work backwards next.
Structural addendum: carry-cheap tails cannot live inside a Stieltjes (Heine-type) machine. Hurwitz tails come from a kernel with residues of size 1 at t = i(n+½) (periodic pattern) — that is what makes it decay like e^{−πt} and its moments Bernoulli-rational. Binomial tails need residues ∝ 1/C(2n,n) ~ 4^{−n}: the kernel then decays only like 1/t² (ρ = 0, no K² closeness) and its moments Σ n^k/C(2n,n) lie in Q + Q·π√3 (level-1 χ₋₃ values), not Q. So the K² arena is closed for G and the carries move to the K¹ arena (linear forms / recurrences).

---
# The K¹ arena (2026-09-26): Zudilin's Apéry-like forms for G under the same ledger
Source: W. Zudilin, "An Apéry-like difference equation for Catalan's constant" (math/0201024; Electron. J. Combin. 10 (2003)). Forms F_n = Σ_t (−1)^t R_n(t), R_n(t) = n!(2t+n+1)·t(t−1)…(t−n+1)(t+n+1)…(t+2n)/((t+½)(t+3/2)…(t+n+½))³ (very-well-poised ₆F₅ at −1) = u_nG − v_n (times 8), with u_nG − v_n = ((−1)^n/4)∫∫ x^{n−½}(1−x)^n y^n(1−y)^{n−½}/(1−xy)^{n+1} dx dy. Recurrence (2n+1)²(2n+2)²p(n)u_{n+1} − q(n)u_n − (2n−1)²(2n)²p(n+1)u_{n−1} = 0, p = 20n²−8n+1, q = 3520n⁶+5632n⁵+2064n⁴−384n³−156n²+16n+7, u₀=1, u₁=7/4, v₀=0, v₁=13/8.
Analytic side: characteristic polynomial λ² − 11λ − 1, the same as Apéry's ζ(2) equation, so |u_nG − v_n|^{1/n} → ((√5−1)/2)⁵ = e^{−2.406} (measured −2.58, −2.51, −2.46, −2.45, −2.44 at n = 25, 50, 100, 150, 200). In the K¹ arena the half-integer shift costs nothing analytically — the Beukers-type integral has the same saddle as for ζ(2) — where in the K² arena it halved ρ.
Arithmetic side, exact ledger (zudilin_ledger.py, n ≤ 300): den(u_n) = 2^{4n−6}; den(v_n) = 2^{4n−5}·D_{2n−1}² exactly — every odd prime p ≤ 2n−1 with exponent 2⌊log_p(2n−1)⌋, both classes mod 4 alike (n=300: 3:10, 5:6, 7:6, 11…23: 4, 29…599: 2), no cancellation anywhere. Per n: 2-adic 4 ln 2 = 2.773 + primes ≤ n: 2.00 (= D_n², what Apéry's ζ(2) pays and can afford) + primes in (n, 2n]: 2.00 (the lattice doubling) = 6.77 against 2.406. Miss 4.37 per n.
Reading: the two taxes of the Hankel ledger reappear at full strength. The half-integer shift (t+k+½) doubles the prime range (+2) and puts 2^{4n} in every denominator (+2.77). The Hankel machine pays these more gently (lattice ≈ +0.4–0.7, 2-adic ⅙ ln 2) but loses 0.92 on the analytic side (ρ = π); the linear-form machine keeps the analytic side and pays 4.4 in arithmetic. Zudilin's "no chance to prove that Catalan's constant is irrational" is this ledger; the Rivoal–Zudilin escape (many β(2k) at once) trades one unknown for a dimension argument.
Backwards target in the K¹ arena: forms with denominators ≤ e^{2.4n} need integer-lattice denominators (D_n², 2.0 per n) and no 2-power — the same specification as in K²: leave the half-integer lattice while staying G-affine.

---
# The untested handoff items, analysed and run (2026-09-27)
The original handoff listed "poles at Gaussian norms (2D lattice)" and "poles on ζ₁₂ norms (12-fold)" as untested. Analytical factorization first:
- 2D (Z[i]) machine: structurally impossible in Fauzan form. The K-dimensional rational family of the 1D machine comes from the Hurwitz recurrence ζ(2,a) − ζ(2,a+1) = 1/a²; the Epstein series E(z,2) = Σ_ω|ω+z|⁻⁴ is periodic under lattice shifts (E(1/2,2) = 4π²G, E((1+i)/2,2) = 2π²G, E(0,2) = (2π²/3)G are all Q·π²G, but shifting by ω gives the same value), so there is no partial-sum structure; partial lattice sums (half-planes, quadrants) have hyperbolic-function remainders; and the theta kernel Σr₂(N)e^{−2πNt} has moments ∝ ζ(2e+2)L(2e+2,χ₋₄) ∝ β(2e+2), not rational. Not run, by analysis.
- One-class machines on a = j + OFF/d (dclass.py; validated: d=2 reproduces the ζ(2) half-pole run, d=4 the quarter one-class run). Target ξ_d = ζ(2, OFF/d): d=3, OFF=1 gives ψ'(1/3) = 2π²/3 + (9/2)L(2,χ₋₃) (irrationality open); d=6 gives ζ(2,1/6) = 2π² + (45/2)L(2,χ₋₃). ρ = 2π for every d (spacing 1), so closeness → −(9/2) ln 2 for all of them; measured at K=40, N=0: −3.0955 (d=2), −3.0934 / −3.1035 (d=3, OFF=1/2), −2.9117 (d=4, N=3).
| d, OFF | K, N | closeness | height | net | den by class of p mod d |
|---|---|---|---|---|---|
| 2, 1 (ζ(2)) | 40, 0 | −3.0955 | 2.950 | −0.145 | — |
| 3, 1 (ψ'(1/3)) | 40, 0 | −3.0934 | 3.787 | +0.694 | p≡1: 1.70, p≡2: 1.95, 3: 0.08 |
| 3, 2 | 40, 0 | −3.1035 | 3.818 | +0.715 | p≡1: 1.59, p≡2: 2.12 |
| 3, 1 | 40, 3 | −2.9067 | 3.484 | +0.577 | p≡1: 1.59, p≡2: 1.76 |
| 4, 1 (π²+8G) | 40, 3 | −2.9117 | 3.835 | +0.924 | p≡1: 2.08, p≡3: 1.51, 2: 0.19 |
Parity of the numerators is a new ledger line. On thirds the numerators b = 3j+1 alternate in parity, so b + b' can be an odd prime: p ≡ 2 mod 3 gets mirror pairs b+b' = p up to 6K (costing 1 each: both members have p-integral pole values, Hermite congruence), p ≡ 1 mod 3 gets pairs b+b' = 2p up to 3K plus the partial-sum 1/p² (3 per pole above p, +1): e.g. p=197 (≡2): 7 pairs → e = 7 exactly; p=97 (≡1): 7 pairs + 7 poles above + 1 → 22 exactly. So a mixed-parity lattice of extent L admits primes up to 2L, an all-odd lattice (d even) only up to L: the effective lattice of the thirds machine is 6K, and its height (3.79) lands next to the quarter machine's (3.84). No one-class machine beats d = 2, and none closes; the lattice-extent law is monotone: d=2: 2.95, d=3: 3.79, d=4: 3.84 at K=40.
| 6, 1 (ζ(2,1/6)) | 40, 0 | −3.1061 | 5.225 | +2.119 | p≡1 mod 6: 3.35, p≡5: 1.30, 2: 0.44, 3: 0.08 |
Sixths (b = 6j+1, all odd, extent 6K): the class p ≡ 5 mod 6 is not free either — b + b' ≡ 2 mod 6 admits b + b' = 4p for p ≡ 2 mod 3 (pairs up to 3K; exponents 50, 34, 32, 31, 21, 12, 10, 6, 4, 2 at p = 41…113, none above 113), while p ≡ 1 mod 6 pays pairs 2p up to 6K plus partial sums (p = 127…229 with exponents 55…4). Rule: for a one-class lattice b ≡ OFF mod d, prime p acquires mirror pairs b + b' = 2kp whenever 2kp ≡ 2·OFF mod d for some k ≥ 1 with 2kp ≤ 2·max b — every class pays at some multiple; only the multiple (hence the range) differs. Height 5.23: the extent law is monotone through d = 6, and the one-class route is closed for all d.
Zudilin's follow-up (math/0210423, "A few remarks on linear forms involving Catalan's constant") proves 2^{4n+o(n)}u_n ∈ Z and 2^{4n+o(n)}D_{2n−1}²v_n ∈ Z (o(n) = O(log₂ n)), matching the exact ledger above. Its partial-fraction proof (formula 8) puts the poles at t + l − n/2 − 3/4 and t + l − n/2 − 1/4 — the two quarter-integer classes — and builds v_n from the partial sums Σ 1/(4μ ± ε)²: the base-4 Hurwitz decomposition of G, i.e. the same two-class quarter-integer lattice that the Hankel ledger found (sign carry = base-2 refinement of the half-integer shift). A second recursion (13), a Rhin–Viola-type 120-element permutation group, and a stated expectation 2^{2M+o(M)}D_{m₁}D_{m₂}H(c) ∈ ZG + Z are all "beyond reach". Section 4 states the Conjecture: a second-order recursion x_{n+1} + a(n)x_n + b(n)x_{n−1} = 0 with a(n), b(n) → a₀, b₀ ∈ Q whose two Perron solutions are rational with geometrically bounded denominators must have rational characteristic roots; it would imply G irrational (if G ∈ Q, ũ_n and ũ_nG − ṽ_n are such solutions with roots ((1±√5)/2)⁵).

## Two arenas, one lattice (synthesis)
| | K¹ arena (Apéry-like forms, per n) | K² arena (Hankel machine, per K²) |
|---|---|---|
| analytic side | e^{−2.406}, unchanged from ζ(2) | U(π) = −ln 9 vs U(2π) = −(9/2) ln 2: −0.92 lost |
| primes in the ζ(2) range | D_n²: 2.00 (affordable) | integer-pole A ≈ 2.3 |
| lattice doubling (primes to 2n / 2K) | +2.00 | ≈ +0.4 to +0.7 |
| 2-adic tax | 4 ln 2 = 2.77 | ⅙ ln 2 = 0.12 |
| miss | 4.37 | 0.585 |
Rhin–Viola-type search over the whole very-well-poised family (rv_search.py; parameters h₀ = n, h_j = η_j n, decay from the Laplace saddle of ∫∫ x^{c₂₁}(1−x)^{c₂₂}y^{c₃₁}(1−y)^{c₃₃}/(1−xy)^{c₁₁+1}, denominators from Zudilin's expectation 2^{2M}D_{m₁}D_{m₂}, exact on his direction): on the α=1 normalisation Zudilin's direction has denominators 2.258 against decay −0.802 (miss 1.456); the best of 13⁴ grid directions plus local refinement, η ≈ (0.42, 0.42, 0.42, 0.20), has 2.176 against −0.662 (miss 1.514; per unit parameter 2-adic 0.34 + primes 1.15 − decay 0.45). Denominators exceed the decay by a factor ≈ 3 everywhere on the grid: no re-parametrisation of this family closes the gap, and Zudilin's "beyond reach" is quantified as a factor of three.
Exact check of the denominators off Zudilin's direction (vwp_forms.py: general partial fractions of the very-well-poised summand in his normalisation; validated: U₁ = U₃ = 0 exactly by the well-poised palindrome, and (U₂, V) = ±8(u_n, v_n) on his direction, e.g. V/U₂ = 13/14, 10699/11682). Best search direction η = (0.425, 0.425, 0.425, 0.2):
| a | den(V) | vs Zudilin's expectation 2^{2M}D_{m₁}D_{m₂} | log den/a | log|U₂G−V|/a | integer form per a |
|---|---|---|---|---|---|
| 20 | 2⁶·3³·5²·7²·11²·13² | ⊂ 2^{18}D₁₅² (deficit 2^{12}·3) | 1.225 | −0.680 | +0.545 |
| 60 | 2^{36}·3⁵·5⁵·11²·13²·17²·19²·23·29²…43² | excess 5¹; deficit 2^{14}, 3, 7⁴, 23, 47² | 1.648 | −0.624 | +1.024 |
| 100 | 2^{68}·3⁴·5⁴·7²·13²·17³·19³·23²…73² | excess 17, 19; deficit 2^{12}, 3⁴, 7², 11², 79², 83² | 1.839 | −0.607 | +1.232 |
| Zudilin's direction, a = 76 (n = 25) | 2^{92}·D₄₉² exactly | deficit 2^{18} only | 2.050 | −0.821 | +1.229 |
So the expectation (24) is neither an upper nor a lower bound off the symmetric direction: den(U₂) itself acquires 5, 17, 19 (the derivatives of the degree-(h₄−1) integer-valued polynomials bring D_{h₄−1}, which the symmetric direction cancels), while carries remove 7, 11, 79, 83 and lower 3 and 23 (mostly primes ≡ 3 mod 4, but 19, 23, 31, 43, … ≡ 3 stay full — a carry rule, not a character rule; not derived). The integer linear form grows at +1.23 per a on both directions: the K¹ family is closed at a miss ≈ 1.2 per unit parameter (≈ 3.7 per Zudilin n), three to four times the K² machine's 0.585.
The determinant absorbs the lattice taxes an order of magnitude better than the linear forms do; its price is the halved ρ. The common root is that G lives on the base-4 refinement of the shift lattice (quarter-integers, two classes), which every construction must pay for either in the prime range or in the decay rate. The K² arena is where the miss is smallest, and there both sides are now closed-form; the exits are (a) a non-lattice positive weight with rational moments and G-affine values (none known; carry-cheap tails are incompatible with exponential decay), (b) the multi-β dimension route (Rivoal–Zudilin), (c) Zudilin's recursion conjecture. Consistently, log P_K(ζ(2)) = −232, −401, −550, −727 at K = 40, 60, 80, 100 is not quadratic but ≈ −1.58·K ln K (ratios 1.58, 1.63, 1.57, 1.58) — enough for the criterion (it beats −K log q), but with no K² margin at all. Whether U = −(9/2) ln 2 exactly, and whether Catalan's U = −ln 9 (two four-decimal coincidences), is a question for the potential-theory side (F. Brown's criterion) — flagged, not claimed.

---
# The multi-β route (2026-09-28): Zudilin's forms under the ledger, a refined arithmetic lemma, and five β values
Source: W. Zudilin, "Arithmetic of Catalan's constant and its relatives" (arXiv 1804.09922, 2018/19): at least one of β(2), β(4), …, β(12) is irrational (six values; Rivoal–Zudilin 2003 had seven). Script: multibeta.py (commands forms / asym / search / mech2). The user asked for this route alone.

## The construction and its ledger (per n)
η = (η₀; η₁…η_s), s odd, 0 < η_j < η₀/2, Ση_j ≤ (s−1)η₀/2; h₀ = η₀n+1, h_j = η_jn+½; R_n(t) = γ_n(2t+h₀)(t+1)_{h₀−1}/∏_j(t+h_j)_{1+h₀−2h_j}; r_n = Σ_{ν≥0}(−1)^νR_n(ν) = Σ_{i even} a_iβ(i) + a₀. Poles at t = −(k+½), N ≤ k ≤ h₀−N−1, multiplicity s_k = #{j : η_jn ≤ k ≤ h₀−1−η_jn} (a tower: the edge poles have few floors, the centre has all s). Palindrome R(−t−h₀) = R(t) ⇒ a_{i,k} = (−1)^i a_{i,h₀−1−k} ⇒ odd β's vanish. a_i = 2^iΣ_k(−1)^k a_{i,k}; a₀ from partial sums Σ_ℓ(−1)^ℓ/(ℓ+½)^i whose denominators are odd numbers < h₀−2N−1 (his shift trick; RZ03 had twice the range).
- closeness: −decay, closed form (Lemma 3: one root x₀ of a degree-(s+1) polynomial; uniqueness now checked numerically on a 20000-grid for every η used).
- height (his Lemma 4/5): naive s·M/n, M = max(h₀−2N−1, η₁n), minus the carry rebate κ = ∫_{1/μ}^∞ φ₀(x)dx/x², φ₀(x) = min_y φ(x,y), φ = the base-p carry count of the digit pair (x,y) = (n/p, k/p) (I verified φ(x,y) = v_p(a_{s_k,k}) + (s − s_k) exactly: the top coefficient's Legendre count plus one unit per floor missing from the tower at that pole).
Reproduced to all printed digits: his η = (31; 10⁵,11⁴,12⁴), s = 13: closeness −100.7396631711, κ = 42.7664565011, height 100.2335434989, net −0.506 per n; his Section-2 construction η = (3;1¹⁷): −16.1123070755, κ = 0.9411124762. (A first implementation sampled φ only at breakpoints and got κ too large by 1.1; the minimum over the pole digit lives on the open intervals — terms like ⌊2(η₀x−y)⌋ drop just after a breakpoint.)

## Exact forms and the true denominators (multibeta.py forms)
Exact partial fractions at every pole to full multiplicity (log-series of the analytic part, fmpq), checked three ways: R(t) reproduced exactly at random rational t; odd a_i = 0; numeric r_n (direct alternating sum, 250–850 digits) equals Σa_iβ(i)+a₀ to all digits. Per prime p ≤ M: exponent of p in den(a₀), den(a₂), … against his bound s·v_p(d_M) − φ₀(n/p).
| η, n | closeness log|r_n|/n | height log den/n | rebate below his bound per n, by range: p ≤ √(2h₀) / middle / (M/2, M] |
|---|---|---|---|
| (31;10⁵,11⁴,12⁴), n=2 | −109.52 | 83.19 | 26.63 / – / 1.28 |
| same, n=4 | −106.24 | 90.21 | 18.88 / 0.00 / 6.55 |
| same, n=6 | −104.84 | 87.33 | 23.88 / 0.00 / 2.06 |
| same, n=8 | −104.05 | 91.74 | 15.67 / 0.39 / 2.21 |
| (3;1⁵), n=60 | +0.28 | 3.28 | 0.58 / 0.11 / 0.06 |
| (3;1⁵), n=100 | – | 3.63 | 0.44 / 0.11 / 0.00 |
| (3;1⁷), n=80 | – | 5.39 | 0.65 / 0.13 / 0.14 |
| (5;2⁵), n=20 | +1.74 | 4.78 | 0.89 / 0.00 / 1.38 |
Reading: in the middle range his bound is exactly tight (no hidden cancellation — the K¹ Catalan forms' "one power less" does not recur for s ≥ 5); the small-prime rebate is large at these n but is O(√n log n), asymptotically nothing; the rebate near p ≈ M is structural (below).

## The refined lemma (this work; proof sketch, then exact verification)
At a pole k of multiplicity s_k the Laurent coefficients are a_{i,k} = [u^{s_k−i}] C_k ∏(1+u/α)/∏(1+u/β) (α, β the other factors at the pole; C_k = a_{s_k,k}). A power of u can cost a factor p only through a factor divisible by p, so v_p(a_{i,k}) ≥ v_p(C_k) − (s_k − i), and if NO denominator factor β at that pole is divisible by p the p-singular part is a polynomial of degree D_k = #{p-divisible numerator factors} and v_p(a_{i,k}) ≥ v_p(C_k) − min(s_k − i, D_k) (derivative saturation). A partial sum of order i carries p^{−i} only at poles with |k − centre| ≥ (p+1)/2 (the PS window: the odd numbers in the partial sum are < 2|k−centre|). With v_p(C_k) = φ(x,y_k) − (s − s_k) (Legendre, p > √(2h₀)):
  v_p(a₀) ≥ −s + min( min_{k∈PS(p)} φ(x,y_k), min_k [φ(x,y_k) + max(1, s_k − D_k)] ),   v_p(a_i) ≥ −s + min_k[φ(x,y_k) + max(i, s_k − D_k)].
Zudilin's bound is the special case "every pole has all s floors, every derivative costs p, every pole feeds the partial sums". Asymptotically (x = n/p): PS window y ∈ [η_min x, η₀x/2 − ½] ∪ mirror, covering a full period once x ≥ 2/μ′ (μ′ = η₀ − 2η_min); D(x,y) = ∞ iff some denominator interval [η_jx − y, (η₀−η_j)x − y] contains a nonzero integer, else D = ⌊η₀x − y + ½⌋ + ⌊y + ½⌋ (the odd multiples of p among the numerator factors). So the refinement changes only x ∈ [1/μ, 2/μ′) — the primes in (M/2, M] — and there it removes several powers at once, because the partial sums at those primes are fed only by edge poles with few floors, while the centre poles, which set φ₀, have no p-divisible denominator factor and their derivatives cannot buy p.
Exact verification (multibeta.py forms prints the finite-n predictor next to the actual exponent; every single-digit prime √(2h₀) < p ≤ M): Zudilin's η, n = 6: all 10 primes exact; n = 8: 14 of 15 exact, p = 23 one unit better than predicted; (5;2⁵) n = 20: all exact; (33; 9,10,10,11,11,12,12,12,13,13,14) n = 2,4,6,8 and (36; 10,10,11,11,12,12,13,13,14,14,15) n = 3,5: never exceeded, exact at n = 3,4,5. Example, Zudilin's η, n = 8, p = 83: his bound p⁷, actual p³ — the only poles feeding 1/83 are k = 80,81,82 of multiplicity 5 with v_p(C_k) = 2, giving 5 − 2 = 3; the multiplicity-13 poles have D = 3 (numerator factors −83, 83, 249 only), so their a_{i,k} are bounded by p^{−3}·p^{−0}.
Asymptotics: for Zudilin's own η the refined rebate is κ_ref − κ = 11/3 exactly (gain 4 on [1/11, 3/31), 3 on [3/31, 1/10)): net −4.17 per n instead of −0.51. Heights per n: his η: 143 − 46.43 = 96.57 against 100.74 decay.

## Five β values (under the refined lemma)
Integer search over η (coordinate descent with single and paired moves from three starts, per η₀; objective net_ref = s·M/n − κ_ref − decay):
| s | β's | best per unit η₀ found | η | closeness | refined height | net per n |
|---|---|---|---|---|---|---|
| 13 | β(2)…β(12) | −0.408 (η₀=35) | (35; 11,11,12,12,12,13,13,13,14,14,14,15,15) | −101.29 | 87.02 | −14.27 |
| 11 | β(2)…β(10) | −0.142 (η₀=40) | (40; 11,12,12,13,13,14,14,15,15,16,17) | −94.71 | 89.04 | −5.66 |
| 11 | | −0.080 (η₀=36) | (36; 10,10,11,11,12,12,13,13,14,14,15) | −86.63 | 83.75 | −2.88 |
| 11 | | −0.021 (η₀=33) | (33; 9,10,10,11,11,12,12,12,13,13,14) | −76.81 | 76.13 | −0.69 |
| 9 | β(2)…β(8) | +0.273 (η₀=36), trend 0.79 → 0.27 from η₀ = 12 to 36 | (36; 9,10,11,11,12,13,13,14,15) | −56.65 | 66.48 | +9.83 |
| 7 | β(2)…β(6) | +0.670 (η₀=29), flat | (29; 7,8,9,10,10,11,12) | −20.28 | 39.7 | +19.4 |
The s = 11 directions close: under the refined lemma, at least one of β(2), β(4), β(6), β(8), β(10) is irrational — six → five. Under Zudilin's bound the same η have net +17 to +26 (they spread the tower: η_j from 9 to 14, so the edge poles have one floor over a long stretch and the top primes (14n, 15n] cost nothing at all: at n = 8 the primes 97…109 have exponent 1 where his bound says 7–9). The exact runs confirm the bound at every prime for these η (above) and the numeric r_n at n = 2, 4 (−84.4, −81.6 per n) is heading to the closed-form −76.8. Searches at larger η₀ (s = 11 to 56, s = 9 to 66) are running; the s = 9 trend is flattening and three β's (s = 7) are out of reach in this family.

## Mechanism II — one more power, sometimes (multibeta.py mech2)
Where the actual exponent is one below the refined prediction, the anatomy is always the same: the leading terms (i,k) of a₀ at p^{−e} have residues whose sum vanishes mod p. In the simple family (3;1⁷), n = 80, at p = 29, 31, 37, 43 the alternating sum Σ_{k∈window}(−1)^k a_{i,k} gains a power of p for EVERY order i separately; at p = 41, 53 only the cross-order total cancels; at p = 47, 59 nothing cancels. For x = n/p ∈ [2,3) the per-i case is explained: a_{i,k} mod p depends only on k mod p on the window, the window contains the full digit block [0, n₀] (n₀ = n mod p), and the global palindrome k → h₀−1−k descends to the local mirror k′ → n₀ − k′ under which (−1)^{k′}a_{i,k′} is odd — so the window sum vanishes and the fixed point k′ = n₀/2 is divisible by p (observed: k′ = 11 missing from the leading set at p = 29). Frequency: 6 of 14 primes in range for (3;1⁷) n = 80, 2 of 16 for (3;1⁵) n = 100, 1 of 15 for Zudilin's η at n = 8 (p = 23, cross-order). Not yet a rule for general η; it can only lower heights further.

## Lenses
- Hyperoperations (level axis): the tower of multiplicities s_k is the whole story of the refinement — the partial sums live on the lower floors near the edges, the carry minimum φ₀ on the top floor at the centre, and Zudilin's bound charges every floor everywhere.
- Hybrid base: φ is a base-p carry count of the digit pair (n/p, k/p); D_k is the count of p-divisible factors, a digit count; the refined κ is again an integral of a step function of the digit x with rational breakpoints (denominators ≤ 2η₀), giving rational rebates (11/3, 26, 20, 86/3).
- Palindrome: R(−t−h₀) = R(t) kills the odd β's and, mod p, its local shadow k′ → n₀−k′ is mechanism II.

## Addendum (later on 2026-09-28): larger η₀, larger n
- s = 11, η₀ = 41…56: every η₀ closes; best per unit η₀: −0.159 at η₀ = 47, (47; 14,15,15,16,16,17,17,18,18,19,20), closeness −106.05, refined net −7.47 per n; η₀ = 54: (54; 16,17,17,18,19,19,20,21,21,22,23) net −8.56 (−0.158 per unit). The per-unit margin saturates near −0.15; the five-β conclusion does not depend on any single direction.
- s = 9, η₀ = 37…66: best +0.222 per unit η₀ (η₀ = 51, (51; 14,15,16,16,17,18,19,20,21)); the per-unit miss sits at 0.22–0.29 for every η₀ from 37 to 66, i.e. ≈ 11–16 per n. Four β values are not reached by this family under the refined lemma; the gap is far beyond the observed mechanism-II rebates.
- Exact runs of the winner (33; 9,10,10,11,11,12,12,12,13,13,14) at n = 12 and 16: numeric r_n equals the exact form to 1147 and 1510 digits; closeness −78.90, −78.48 per n (closed form −76.81); the refined bound holds at every single-digit prime and is attained exactly at all of them for n = 16 (n = 12: p = 37, 43 one unit lower, mechanism II). Heights 72.8, 73.6 per n, still 17 and 13 per n of small-prime (p ≤ √(2h₀)) share above the asymptotic 76.1 − … structure; the middle range is tight, (M/2, M] rebates 17.7 and 14.3 per n.
- Mechanism II at the winner (n = 8, p = 41) is of the cross-order kind (no per-i vanishing).

## The proof (2026-09-28, later): refined_lemma_proof.md
Written out in full in the folder: Lemma A (derivative saturation, term-by-term valuations), Lemma B (the exact identity v_p(a_{s_k,k}) + (s − s_k) = φ(n/p, k/p) + [k ≡ m mod p], by Legendre counts — the odd-number count needs three parity observations), Lemma C (window), Theorem 1 (per-prime bound for a₀ and the a_i, with the centre pole and the k ≡ m poles handled), Theorem 2 (denominator sequence D_n and (log D_n)/n → ∫ max(0, s − κ̂) dx/x²; Zudilin's Lemma 4 is cited only for p ≤ √(2h₀), an O(√n) contribution), Theorem 3 (five β values, with the window table: gains 9, 8, 6, 3 on [1/15,1/14), [1/14,1/13), [1/13,1/12), [1/12,1/11), total 26; brute-force check of both κ's). Two points the numerics had glossed over are settled there: the height is ∫ max(0, s − κ̂), equal to sμ − κ_ref only because κ̂ ≤ s (true for all η used; the maximum κ̂ = s = 11 sits on [1/15, 1/14)), and the minima over y must be over closed sets — evaluation on open pieces is exact unless p | n (double breakpoints), a negligible set of primes.

## Caveats (honest ledger)
1. The refined lemma is now proved in refined_lemma_proof.md (self-contained except for Zudilin's Lemma 2/3 for the decay and the first part of his Lemma 4 for primes p ≤ √(2h₀)); it has not been refereed by anyone else. 2. RZ03 §7 (Math. Ann., not on arXiv) was not read; Zudilin 2018 bounds the partial sums crudely by d^i_{h₀−2N−2} for all poles, so the refinement is not in that paper. 3. The search is heuristic (descent); better η may exist. 4. Everything is the PNT limit; the finite-n heights are still dominated by small primes at n ≤ 8.

---
# Mechanism II resolved to a digit rule in the simple family; the odd-floor residue; literature and transfer checks (2026-09-28, continued)
Sweeps (multibeta.py sweep; JSON in the folder): (3;1⁵) n = 20…140 (818 (n,p) pairs, 192 events), (3;1⁷) n = 20…120 (594 pairs, 139 events), s = 11 winner n = 2…14 (147 pairs, 8 events). Refined bound never exceeded. Content of the winner's integer forms (gcd of D·a_i) is 1 at n = 4, 8, 12, 16.
## The rule (simple family (3;1^s), primes √(2h₀) < p ≤ n, x = n/p, n₁ = ⌊x⌋, n₀ = n mod p; n even so n₀ ≡ n₁ mod 2)
| n₁ | rate | type | detail |
|---|---|---|---|
| 2 | 142/191 (74%) | all per-order | 100% on every tenth of {x} except {x} ∈ [0.3, 0.5) (2/37) and the case p | n (x integer, 0%) |
| 4 | 13/21 (62%) | per-order | same shape: .1–.2 and .5–.7 all events, .3–.4 none |
| 1 | 35/537 (7%) | mostly cross-order | clustered at one prime over runs of consecutive n: p = 29 for n = 30…42, p = 41 for n = 62…80 (s = 7), p = 71 for n = 108, 110 |
| 3 | 2/69 (3%) | – | |
p mod 3 and p mod 4 are irrelevant (19–28% each); gap 2 or 3 occurs 3 times per family, always at n₁ = 2. Same rule for s = 5 and s = 7.
## Mechanism, verified pole by pole
- Even n₁ = the half fold. n = 80, p = 31 (n₀ = 18): the 18 leading poles k = 80…88, 90…98 are paired by the local mirror k ↦ h₀−1−k (mod p) — (80,98), (81,97), …, (88,90) — and every pair sums to 0 mod p at every order i = 1…5; the fixed point k = 89 ≡ m (mod 31) is the one pole absent from the leading set (its top coefficient is divisible by p). So the per-order event is the global palindrome reduced mod p, acting on the digit block, for all orders at once.
- Odd n₁: no mirror. n = 70, p = 41 (leading k = 70…84, residues 29…40, 0…2, m ≡ 23); n = 36, p = 29 (k = 36…39); n = 110, p = 71 (k = 110…129): the leading set is one cell, its image under the mirror lies outside the window, per-order sums are nonzero, and only the total over orders vanishes. These are the events beyond the fold.
- Exact reformulation for the odd case: a₀ (lower half) = Σ_{k<m} Σ_{ν=−m}^{−k−1} (−1)^ν block_k(ν), block_k(t) = Σ_i a_{i,k}(t+k+½)^{−i}, ν integer. A block is p-singular at ν iff 2ν + 2k + 1 ∈ {±p, ±3p, …}, i.e. iff ν is the midpoint between the pole k and a pole an odd multiple of p away; and since R(ν) = 0 at every such ν, the p-singular parts of the blocks of one residue class mod p on the two sides of ν cancel exactly. The cross-order total of a leading pole k is its block at the midpoint ν_k = −k − (p+1)/2, which equals minus the singular parts of the blocks at k ± p, k ± 2p, … at the same point. Mechanism II in the odd case is therefore a statement about poles one prime apart interacting at their midpoint — the midpoint–prime interaction in the literal sense — and the rule for it (which primes, which n) is not yet found; the winner η's events (n ≤ 14: x = 0.11–0.45, 4 per-order, 4 cross) are too few, sweeps at larger n are running.
## Asymptotic weight of the rule (simple family)
The even-floor events give one power at the primes with ⌊n/p⌋ ∈ {2, 4, 6, …} outside the dead zones: Σ_{j even}(1/j − 1/(j+1)) = 1 − ln 2 = 0.307 per n before dead zones, ≈ 0.25 after — small against heights of 4–16 per n, and it is the x > 1 regime; the Zudilin-type η live at x < 1 where the rule is not yet known.
## Checks
- Fischler 2019 (arXiv 1904.02402): Padé approximation, Shidlovsky lemma, Siegel's criterion — a different method; no overlap with the refined Lemma 4. RZ03 (Math. Ann. 326) still unread: no PDF on Zudilin's publication page (entry 26), HAL search blocked from here.
- Zudilin's SIGMA 2018 ζ-construction (arXiv 1801.09895) read: R_n(t) = 2^{6n} n!^{s−5} ∏_{j=0}^{6n}(t−n+j/2)/∏_{j=0}^n(t+j)^{s+1}, poles at integers of one multiplicity s (no tower, no η, no Φ_n), odd fold R(−t−n) = −R(t), the ζ(3) term removed by the twist-by-half combination 7r_n − r̂_n. Only the window and saturation would transfer, and Zudilin states that the Zu04 arithmetic already reduces s far below 25 without reaching 9; so applying the refinement to the elementary ζ-forms cannot touch "one of ζ(5), ζ(7), ζ(9), ζ(11)". The target worth the machinery would be the refinement on top of the general Zu04 construction (towers, Φ, integer poles, odd fold) — a separate project, not started.
- s = 9 closed at +0.22 per unit η₀ (η₀ = 51), flat to η₀ = 66.

## The per-order rule, stated and tested (2026-09-28, late)
Anatomy at the pole level (winner η, n = 18, p = 37; n = 20, p = 61; n = 32, p = 97): in a per-order event the leading poles form blocks of consecutive k, each block centred on a pole c ≡ m (mod p) which is itself absent (Lemma B's +1), and within each block k ↔ 2c − k cancel exactly at every order; in a mirror-closed non-event (n = 16, p = 37: residues {12,13,14} ∪ {33,34,35}, m ≡ 5) the blocks are mirror images of each other but centred elsewhere, and their class sums do not vanish; at n = 16, p = 59 the only cancelling residue pair is (27, 29), the one straddling m ≡ 28.
**Rule (conjecture, computable from φ alone).** Per-order mechanism II at (n, p) ⟺ the set where φ(n/p, ·) attains its minimum over the lower window arc y ∈ [η_min x, η₀x/2 − ½], reduced mod 1, is a single arc containing the fixed point y₀ = η₀x/2 (mod 1) in its interior. (φ(x,·) is symmetric under y ↦ η₀x − y, so such an arc is automatically symmetric; the excluded centre is the pole k ≡ m whose top coefficient carries one more p.)
Scores against the exact sweeps (per-order events only; pairs with p | n excluded):
| data | predicted & event | predicted & none | event & not predicted | precision | recall |
|---|---|---|---|---|---|
| (3;1⁵) n = 20…140 | 155 | 5 | 13 | 0.97 | 0.92 |
| (3;1⁷) n = 20…120 | 104 | 4 | 12 | 0.96 | 0.90 |
| (5;2⁵) n = 20…200 | 143 | 3 | 12 | 0.98 | 0.92 |
| (33;9,…,14) n = 2…40 | 36 | 0 | 0 | 1.00 | 1.00 |
| (36;10,…,15) n = 2…24 | 8 | 0 | 1 | 1.00 | 0.89 |
| (31;10⁵,11⁴,12⁴) n = 2…20 | 11 | 1 | 0 | 0.92 | 1.00 |
The false positives sit at x = 2.316–2.326 and 2.194–2.197 (cell edges, finite-n window offset 1/(2p)); the false negatives are the sporadic odd-floor per-order events (x ∈ (1,2) ∪ (3,4)), a separate structure. For the simple family the rule gives exactly: even ⌊x⌋ with {x} ∉ [1/3, 1/2) — the sweep's dead zone is [1/3, 1/2), not [0.3, 0.5).
Asymptotic worth of the rule, ∫[rule holds] dx/x²: winner (33;…) 0.717 per n (x-ranges (0.308,0.333), (0.445,0.455), (0.485,0.5), (0.667,0.697), (0.95,1.0), …), so its margin would be 0.686 + 0.717 = 1.40 per n; (36;…) 0.658; (40;…) 0.559; Zudilin's η 0.918 (his −0.506 → −5.09 with both refinements); (3;1^s) 0.2545 = Σ_{j even}(1/j − 1/(j+⅓) + 1/(j+½) − 1/(j+1)). It does not reach s = 9 (needs ≈ 11 per n).
What a proof needs: for a block centred at c = m − jp, the pair (c−d, c+d) is the image of the exact palindromic pair (m−d, m+d) under a shift by jp; the residues of all non-divisible factors are unchanged by the shift, and the unit parts of the p-divisible factors shift by 2j (numerator) or j (denominator) — the sign of the pairing is the Lucas/Wilson bookkeeping of those shifted digits. Not done.
The cross-order events (odd ⌊x⌋ in the simple family: p = 29 for n = 30…42, p = 41 for n = 62…80, p = 71 for n = 108, 110; and x ≈ 0.08–0.2 for the winner) have no mirror structure at all (every mirror class sum nonzero) and remain the events beyond the fold.

## Cross-order events: two kinds, and the primes that carry them (2026-09-28, last)
Per-pole test (block_k evaluated at its own midpoint ν_k = −k − (p+1)/2, i.e. halfway to the pole k + p): in the simple family's odd-floor cross events every leading pole vanishes on its own — (3;1⁷) n = 70, p = 41: 15/15 poles; n = 64, p = 41: 12/12; (3;1⁵) n = 36, p = 29: 4/4; n = 110, p = 71: 20/20. Since R(ν_k) = 0, block_k(ν_k) ≡ −block_{k+p}(ν_k) modulo p-integral terms: the statement is symmetric in the pair (k, k+p) and is a local supercongruence at the midpoint of two poles one prime apart. In the winner's small-x cross events (n = 24, p = 311; n = 8, p = 41) no pole vanishes on its own; only the running total over the cell returns to 0 at its last pole. So mechanism II has three shapes: per-order (mirror blocks around k ≡ m; rule found), per-pole midpoint (local), and cell-sum.
Which primes carry the per-pole kind (all cross events in the sweeps, x ∈ (1, 2) unless stated):
| family | prime (p mod 3, mod 4, mod 9) | n-run | x-range |
|---|---|---|---|
| (3;1⁵) | 29 (2, 1, 2) | 30…42 | 1.03–1.45 |
| (3;1⁵) | 71 (2, 3, 8) | 108…140 (sweep end) | 1.52–1.97 |
| (3;1⁷) | 41 (2, 1, 5) | 62…80 | 1.51–1.95 |
| (3;1⁷) | 53 (2, 1, 8) | 80…104 | 1.51–1.96 |
| (5;2⁵) | 71 (2, 3, 8) | 72…88 | 1.01–1.24 |
| (5;2⁵) | 89 (2, 1, 8) | 156…176 (and n = 64, x = 0.72) | 1.75–1.98 |
| (5;2⁵) | 83, 227 (2, 3, 2) | single n at x ≈ 0.67–0.70 | |
| (5;2⁵) | 43 (1, 3, 7) | 112, 116 | 2.60–2.70 |
Every prime carrying a run is ≡ 2 (mod 3), but most primes ≡ 2 (mod 3) in the same x-range carry nothing (for (3;1⁵): 17, 23, 41, 47, 53, 59, 83, 89, 101, 107, 113 all 0 events over 7–29 pairs each), and the carrying primes change with s (29, 71 for s = 5; 41, 53 for s = 7). So "p ≡ 2 mod 3" is at most necessary; the selecting condition (a supercongruence for the pair (k, k+p) at its midpoint, depending on s and on the digit n₀) is open. This is the "beyond the fold" structure: not a reflection, a local identity between two poles a prime apart.

## Rivoal–Zudilin 2003 read (user-supplied PDF, 2026-09-28)
Construction (their §2): q even ≥ 4, η = (η₀; η₁, …, η_{q−1}; η_q) with an extra INTEGER parameter h_q = η_q n + 1; R_n(t) = γ_n (2t+h₀)(t+1)_{h_q−1}(t+1+h₀−h_q)_{h_q−1}/∏_{j<q}(t+h_j)_{1+h₀−2h_j} — two numerator blocks with a gap, so the sum can only be started at t = 1 − h_q and the partial sums reach odd numbers up to 2(η₀−η_j−η_q)n (this is the d_{2n} of Zudilin's 2018 Remark 1; his full block (t+1)_{h₀−1} starts the sum at the midpoint and halves them). Forms in 1, β(2), …, β(q−2); odd β's killed by the ±-reciprocity of P_{jn}. Theorem 1 (Nesterenko): dim ≥ (1+o(1)) log a/(2 + log 2). Theorem 2: q = 16, η = (101; 28,29,…,42; 45) — a steep consecutive tower — φ = e^{−594.58616762}, denominators e^{593.98582857}: margin 0.60 per n, seven values β(2)…β(14).
Arithmetic: Lemma 6 gives d_{m₁n}⋯d_{m_{q−1−j}n}·16^{η_q n}·P_{jn} ∈ Z[z] with m₀ = max(η₀−2η₁, η_q), m_j = max(m₀, 2(η₀−η_j−η_q), 2(η_q−η_j)); via (24) d^{q−1−j}_{m₀n} 16^{η_q n} c_{kj} ∈ Z (uniform derivative cost) and (25) ∏_{i≤j} d_{max(2h₀−2h_i−2h_q, 2h_q−2h_i)} × (partial sum of order j) ∈ Z — the partial-sum lcm range is tied to the pole class j (order-j coefficients exist only for k in the j-th factor's range), a partial tower-awareness at the level of lcm ranges. §7 (Lemmas 10–13): per-factor p-adic estimates ord_p(H(t)(t+k))^{(j)} ≥ −j + ⌊(b−a−1)/p⌋ − ⌊(k−a)/p⌋ − ⌊(b−1−k)/p⌋ and ord_p G^{(j)}(−k+½) ≥ −j + ⌊⌊(a−k)/p⌋⌋ − ⌊⌊(b−k)/p⌋⌋ − ⌊(a−b)/p⌋ (⌊⌊x⌋⌋ = ⌊2x⌋ − ⌊x⌋), combined by Leibniz into ord_p D^λ(R·(t+k−½)^{q−1}) ≥ −λ + ϖ₀(n/p, (k−1)/p) ≥ −λ + ϖ(n/p), ϖ = min over all y; Π_n = ∏ p^{ϖ(n/p)} and Lemma 13 gives log Π_n/n → ∫ϖ dψ − ∫_0^{1/m} ϖ dx/x². So: every derivative costs one p (no saturation), every pole may feed the partial sums (no window), and the minimum is over all y. **The refined lemma is not in RZ03 either**; its per-factor Lemmas 10–11 are the ingredients Lemma B here re-derives as an exact identity.
Two remarks of theirs bear on our ledger. (i) §9: for q = 4 (Catalan alone) the group ⟨b, S₃⟩ of order 24 (a subgroup of Rhin–Viola's 120) acts on the parameters, but "the denominators of the linear forms given by Lemma 6 are too large to use the arithmetic arguments in [RV], [Zu2]"; our exact K¹ ledger shows why no group can help: on the symmetric direction the exact denominators 2^{4n−5}D_{2n−1}² are the truth, not a bound. (ii) They verified numerically to n = 1000 that 16^n u_n ∈ Z and d_{2n}²16^n v_n ∈ Z — matching the exact ledger (den u_n = 2^{4n−6}, den v_n = 2^{4n−5}D_{2n−1}²).
Not run: the refined lemma on the RZ03 family itself (different zero geometry: window measured from t = 1 − h_q, not from the midpoint). Zudilin 2018 argues the midpoint construction dominates it on the a₀ axis; the extra integer parameter η_q is the one freedom it has that the 2018 family lacks.

## Mechanism II, all three shapes, from one formula (2026-09-28, last)
On a digit cell (all p-divisible factor multisets constant; digits a = α/p, b = β/p), Lemma A gives a_{i,k} = unit(C_k)·p^{v_p(C_k)−(s_k−i)}·h_{s_k−i} + higher, with h_r the Taylor coefficients of the digit profile F(w) = ∏_{α∈A^p}(1 + 2w/q_α)/∏_{β∈B^p}(1 + w/b_β) (q = 2α/p odd). Hence the p-singular part of a₀ contributed by pole k is, at leading order,
  (−1)^k · unit(C_k) · p^{v_p(C_k) − s_k} · T_k,   T_k = Σ_{q odd ≤ 2J−1} (−1)^{(qp+1)/2} (−2/q)^{s_k} Σ_{r<s_k} h_r (−q/2)^r,
J = number of odd multiples of p below 2(m−k). T_k depends only on the cell's digits, s_k and J. Three ways for Σ_{k leading} (−1)^k unit(C_k) T_k to vanish mod p:
1. **Per-order (Theorem 4, proved in refined_lemma_proof.md §8):** the leading set is a union of cells symmetric about fixed-point poles c ≡ m; within such a cell unit(C_{2c−k}) ≡ −unit(C_k) (the centre factor 2(m−k) gives −1, every other family gives (−1)^{#p-unit members} and the total parity is even) while T is constant, so the sum cancels pairwise — for every order i.
2. **Per-pole:** p | numerator(T_cell). For the simple family at x ∈ (1,2) (J = 1, cell at the window edge) the configurations and T are: s=5: A={−1,1,3}, B={1}⁵: T = 145/24 = 5·29/24; A={−3,−1,1,3,5}: T = 781/120 = 11·71/120. s=7: A={−3,−1,1,3,5}, B={1}⁷: T = 15211/720 = 7·41·53/720; A={−1,1,3,5}: 63/4 (no large prime). s=9: A={−1,1,3}: 8723/128 = 11·13·61/128. s=11: A={−3,−1,1,3,5}: 66781/256 = 11·13·467/256; A={−1,1,3,5}: 187187/960 = 7·11²·13·17/960. The carrier primes observed in the sweeps — 29, 71 (s=5); 41, 53 (s=7); 61 (s=9); 13, 17 (s=11) — are exactly the large prime factors of these numerators, and primes sharing a configuration but not dividing (43, 31, 73, 47) never fire. So the "beyond the fold" events are not a fold at all: they are the primes dividing a fixed hypergeometric-type rational attached to the digit configuration; sparse, s-dependent, no residue-class rule (the earlier p ≡ 2 mod 3 reading was a small-sample accident: s = 9 has 61, s = 11 has 13 and 43). Note F(−1/2) = 0 exactly whenever the digit q = 1 is present (the numerator factor α = p/2), so T is a truncation tail.
3. **Cell-sum:** the weighted sum over several cells with different T vanishes (the winner's small-x events). Defined by the same formula; an accident of the weights.
Full recomputation of the leading residue of a₀ from these ingredients (cross_theory.py) over every (n, p) pair of all seven sweeps: 5663/5663 agree (100%) — simple s = 5, 7, 9, 11 (818, 594, 1673, 1070 pairs), winner n ≤ 40 (883), (36;…) (421), Zudilin's η (204). Per-pole share of the cross events: 24/24, 23/23, 14/24, 7/11, 1/3, 1/6, 2/3; the rest cell-sum. Written up as Theorem 5 (with Lemmas A′, C′) in refined_lemma_proof.md §9: the exponent of every prime in den(a₀) is now computable from base-p digits alone. Mechanism II is closed as a description; what remains open about it is only which x-ranges the cell-sum accidents favour, and they carry no density.

## The digit-only exponent (2026-09-28, last): Theorems 1 + 5 as an algorithm
`multibeta.py digitexp ETA n p` computes the exponent of p in den(a₀) from base-p digits only (Lemma B counts for Φ_p, leading set, unit(C_k) by a boundary recurrence along the poles from an Anton–Wilson start, digit profiles mod p, the residue of Theorem 5). Validation against the exact exponents of all sweep pairs with p ∤ n: 4538/4538. Theorem 6 (proof doc §10) turns Theorem 4 into the rigorous rebate ∫G dx/x² = 0.7172 per n for the five-β direction: margin 1.4035 per n. The calculator is being used to scan the winner at n = 60…200 over all primes (digit_scan_winner.json) for the total mechanism-II rebate and for any density in the cell-sum accidents.
Scan result (digit_scan_winner.json; winner η, all primes p² > 2h₀, p ≤ M, p ∤ n): n = 60: 136 primes, rebate below Theorem 1 = 0.686 per n = 0.413 per-order (5 primes) + 0.273 cell-sum (3 primes at x = 0.335, 0.258, 0.193); n = 100: 217 primes, 0.473, all per-order (9); n = 150: 309 primes, 0.534, all per-order (14); n = 200: 400 primes, 0.581, all per-order (19), no other drops. The per-order rebate rises toward its density 0.717 (the G-intervals are narrow, so the prime counts grow slowly) and nothing else survives: the total mechanism-II rebate for the five-β direction is the Theorem 6 density, and the rigorous margin is 1.4035 per n. Primes not window-dominated (the top end x < 2/μ′): 8, 17, 17, 23 — there Theorem 1 stands and no drop was seen. Mechanism II is closed.

---
# Priority correction (2026-09-28, last): Lai–Zhou 2021 already have five β values
Li Lai and Li Zhou, "At least two of ζ(5), ζ(7), …, ζ(35) are irrational" (arXiv 2103.00904, 2021), Theorem 6.1: at least one of β(2), β(4), β(6), β(8), β(10) is irrational — Zudilin's construction with s = 11, η = (94; 32,32,32,32,33,34,35,36,37,38,39), and a different normalization γ̃_n = 4^{h₀−1}∏_{j=1}^s(h₀−2h_j)!/n!^{η₀} (the numerator (t+1)_{h₀−1} split into η₀ blocks of length n instead of Zudilin's three blocks), which changes the provable carry function to φ̃ (−η₀⌊x⌋ in place of −2⌊η₁x⌋ − ⌊(η₀−2η₁)x⌋) and M̃ = h₀ − 2N − 1; their margin ≈ 0.21 per n. So the statement of Theorem 3 in refined_lemma_proof.md is theirs; the handoff and memory are corrected.
What the present work adds, checked on their direction: (i) the normalization is a constant factor of the forms, so the primitive integer form and its true denominators do not depend on it — only provable bounds do; the refined bound is covariant under normalization and near-exact, so it measures the truth: for their η, Zudilin's bound (his γ) gives net +16.45 per n, the refined bound gives −6.21 per n (30× their margin); the exact denominators at n = 2, 4 obey the refined bound at every single-digit prime (exact at n = 2). (ii) Since the refined bound is essentially the truth (Theorem 5), the failure of s = 9 in the search (+0.22 per unit η₀, ≈ 11 per n, flat over η₀ = 12…66) is a statement about the true denominators: four even β values cannot be reached by this family, for any normalization. (iii) The simpler direction (33; 9,…,14) with margin 1.40 per n (Theorem 6) is a second proof of their theorem. Lai–Zhou's Remark 4.3 asks whether the arithmetic of the coefficients is "even better" than their bound via hypergeometric transformations (the Krattenthaler–Rivoal denominator conjecture); Theorems 1, 4, 5 answer this for the β-forms: the true exponents are s − φ_min minus the residue drop, with Theorem 4's density and the two accidental shapes — nothing else.

# The cubic construction (2026-09-28/29): thirds instead of halves, both types, measured and closed
Motivation: the user's 3-fold stance (fold over three terms, cubic base, the Pythagorean/rational picture, not the complex plane). Concretely: replace the half-integer lattice of Zudilin's β forms by the thirds lattice; the character is then χ₋₃ and the values L(i, χ₋₃) = 3^{−i}(ζ(i,1/3) − ζ(i,2/3)).

## Literature first (checked 2026-09-28)
- Calegari–Dimitrov–Tang, "The linear independence of 1, ζ(2), and L(2,χ₋₃)" (arXiv 2408.15403, Aug 2024): L(2,χ₋₃) is irrational, and 1, ζ(2), L(2,χ₋₃) are Q-linearly independent — by an arithmetic holonomy bound applied to Zagier's construction, not by linear forms. So any "one of L(2,χ₋₃), …, L(2a,χ₋₃)" statement is subsumed; a new theorem would have to exclude L(2) or give two values.
- Fischler 2019 (arXiv 1904.02402): asymptotic counts for L-values of any Dirichlet character, no small cases.
- Lai–Zhou 2021 (arXiv 2103.00904): their numerator roots at 1/3, 2/3 (and 1/2) produce ζ-values (S_{1/3} + S_{2/3}-type combinations with weights 1, 2^i, 3^i for the Vandermonde elimination), not L(i,χ₋₃); their Lemma 3.1 carry function ν(x,y) and Lemma 4.2 Stirling formula are the templates used below.

## Type B — roots at thirds, poles at integers, S_{1/3} − S_{2/3} (cubic.py, cubic_asym.py): dead
R(t) = (2t+m₂n)∏_{i≤m₁n}(t−i)(t+m₂n+i)·(t−m₁n+1/3)_L(t−m₁n+2/3)_L/∏_{j=1}^q(t+δ_jn)_{(m₂−2δ_j)n+1}, L = (2m₁+m₂)n, q odd, degree condition (q−2)m₂ − 6m₁ − 2Σδ_j ≥ 0. Verified exactly (partial fractions reproduce R; odd ρ_i vanish; S_{1/3} − S_{2/3} = a₀ + Σ_{i even}ρ_i3^iL(i,χ₋₃) to 100+ digits). NET = log|S·D/g|/n is normalisation-invariant; the Lai–Zhou normalisation 3^{3L}∏K_j!/n!^{6m₁+2m₂} makes both sides exponential.
- (1,4,[0]⁵) (q = 5, L(2), L(4)): asymptotic closeness +32.435 (Stirling peak, x₀ = 0.464), height q·M′ − κ = 20 − 21.984 = −1.984 (the normalised form is nearly integral: the size is the whole problem), NET +30.45; exact NET +25.02 (n = 12), +25.69 (n = 16), converging upward. (1,3,[0]⁷): +24.38. Search q = 7: best +28.34 (m₂ = 4) rising to +53 at m₂ = 11; q = 9, 13 not finished (irrelevant).
- The cancellation between S_{1/3} and S_{2/3} is real but weak: c_n = (1/n)log|S_{1/3} − S_{2/3}|/S_{1/3} = −1.058, −1.354, −1.240, −1.303, −1.301, −1.328, −1.335 at n = 4, 8, 12, 16, 24, 32, 48 → ≈ −1.34 per n, with an oscillating sign (a Fourier mode at frequency 1 of a Gaussian peak: −2π²/|φ̂″(x₀)|, |φ̂″| ≈ 14.8). Nowhere near the +30.
- Why it is bloated: the thirds root blocks must span the whole pole range (length (2m₁+m₂)n each), so the summand has an interior peak of size e^{+32n}; everything else is small corrections. Verdict: the type-B cubic difference cannot reach net < 0 for any q tried.

## Type A — poles on the two third-cosets, sum over the integers (thirds.py): the right cubic analogue of Zudilin's β forms
R(t) = (2t+h₀)(t+1)_{h₀−1}/∏_{j=1}^q[(t+η_jn+1/3)_{K_j}(t+η_jn+2/3)_{K_j}], h₀ = η₀n+1, K_j = (η₀−2η_j)n+1, n even. Palindrome R(−t−h₀) = −R(t), pairing the 1/3-pole −(k+1/3) with the 2/3-pole −(k′+2/3), k+k′ = h₀−1, so B_i = −(−1)^iA_i (A_i, B_i = sums of a_{i,k} over the two cosets).
**Which sums are interesting (full-line criterion).** For a twist c(ν) with c(−ν−h₀) = σc(ν) and R(−t−h₀) = εR(t), the numerator zeros at −1, …, −(h₀−1) give Σ_{ν∈Z}c(ν)R(ν) = (1 + εσ)·Σ_{ν≥0}c(ν)R(ν). If εσ = +1 the half-line sum is half the full-line sum, hence a sum of residues of R·(kernel), i.e. an explicit identity in powers of π (and √3) with no rational part — trivial. If εσ = −1 the full-line sum vanishes and the half-line sum is a genuine linear form with a rational part a₀ from partial sums. Here ε = −1 (n even), so the untwisted sum r = Σ_{ν≥0}R(ν) is the interesting one and the alternating sum is a π-identity; for h₀ even (n odd, η₀ odd) ε = +1 and nothing on the integer lattice is interesting. Twists must be functions of ν mod 2 only: a twist that depends on the 2-adic valuation of 3ν+1 keeps the 1/3-coset sums in span(ζ, L(χ₋₃)) but not the 2/3-coset sums (3ν+2 has the opposite parity), so period-4 or coset-valuation twists bring in characters mod 12 (β(i), L(i,χ₁₂)). Consequence: exactly one clean interesting form per rational function. Checked numerically for (3;1,1,1), n = 2, 4: Σ(−1)^νR(ν) = −½·Σ_poles Res[R·π/sin(πt)] to 10⁻⁴⁷ (the π-identity), and the untwisted full-line sum Σ_{ν∈Z}R(ν) = −Σ_poles Res[R·π cot(πt)] = 0 to 10⁻⁴⁷.
**Structure (verified exactly at n = 2…8, several η, form = numeric to 10⁻⁴⁷):** r = a₀ + Σ_{i even}3^iA_i·L(i,χ₋₃) + Σ_{i odd ≥ 3}(3^i−1)A_i·ζ(i). The thirds lattice mixes the odd zeta values in: the antisymmetric combination of the two cosets is χ₋₃, the symmetric one is ζ; the palindrome supplies one relation per parity, so both survive. On the half lattice there is one coset, the twist (−1)^ν carries the character, and the symmetric part (β(odd) ∈ π^{odd}Q) is what the palindrome kills. This is the precise sense in which "the third residue is ζ": a form in 1, L(2,χ₋₃), ζ(3), L(4,χ₋₃), ζ(5), … . Trivially true as an "or" statement (ζ(3)); to be new the ζ(3) coefficient must be removed.
**Size.** The untwisted positive sum is dominated by the peak of R(νn) at x₀ (Stirling); no twist cancellation is lost, because Zudilin's β forms have none either: for his η = (31;10⁵,11⁴,12⁴), log R_n(0)/n → −100.7454 (n = 4000) and the interior peak (x ≈ 0.010) is −100.7350, against the true closeness −100.7397 — the alternating sum is the size of its own terms; the smallness is in the terms, the twist is only for the arithmetic (β). Measured for the thirds forms: normalised closeness converges to the Stirling peak (e.g. (3;1,1,1): −11.16 at n = 20, −11.36 at n = 60, peak −11.5625 at x₀ = 0.0561, boundary −11.6136).
**Ledger (exact primitive forms; NET = log|r·D/g|/n, invariant; normalisation N_n = ∏_j(K_j!)²·3^{h₀−1−2ΣK_j}/(h₀−1)! for the split):**
- (3;1,1,1) — form in 1, L(2,χ₋₃), ζ(3): NET −2.55, −2.46, −1.29, −1.57, −0.16, −0.48, +0.11, +0.16, +0.13, +0.20, +0.31 at n = 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 30; normalised height 8.66 → 11.55 (still creeping up, ~log n), closeness → −11.56. Asymptotic NET ≈ +0.3 to +0.5 per n. Half-lattice control, same η, Zudilin's β form (3;1,1,1): closeness +3.007, height 2.059, NET +5.066 per n. The two-coset construction is 4.7 per n better than the half construction at the same η — the second coset doubles the denominator degree (closeness) while the pole order (hence D^q) stays q — but it is not enough.
- (5;2⁵) — 1, L(2), ζ(3), L(4), ζ(5): peak −24.494 (x₀ = 0.105); NET −6.65, −1.99, −0.82, +0.32, +2.23, +2.00, +2.62 at n = 2, 4, 6, 8, 10, 12, 16 (control +10.09, refined +8.09). (5;1,1,1): peak −27.78, NET +9.85 (n = 2) → +13.17 (n = 12) (control +10.34). (5;2,2,2): +2.47 (n = 8), peak −11.88. (5;1,2,2): +2.13 (n = 8), peak −17.56. Zudilin's own η on thirds, n = 2: normalised height 611.15, closeness −574.16 (peak −582.87), NET +36.99. (7;3,3,3), (9;4,4,4): degree condition fails (6η₀ − 12η ≥ η₀ needed).
- Quarter lattice (poles at ±1/4, form in 1, G = β(2), ζ(3)), (3;1,1,1): NET +0.27, +1.59, +2.50, +3.34, +2.78 at n = 4, 8, 12, 16, 20 — worse than thirds (the 2-adic side of the quarter lattice costs more than the 3-adic side of thirds).
- Removing ζ(3). (i) Second derivative with a cubed numerator (t+1)³_{h₀−1}, sum of R″(ν): the orders shift by 2, the residue slot i = 1 (automatically zero) lands on ζ(3), L(2) is absent — verified exactly: (5;1,1,1) gives a form in 1, L(4,χ₋₃), ζ(5) — but the cube triples the numerator carries: NET +28.7 (n = 2), +30.1 (n = 4), peak −15.37 against heights 42–45. (A plain R″ of the simple-zero numerator loses the zero-cancellation of the partial sums: primes up to 3h₀ enter with exponent q — (3;1,1,1), n = 4: NET +13.1.) (ii) A second lattice: Σ_ν R(ν+1/2) has the same A_i with weights −(2^i+1) [L] and (2^i−1) [ζ], so 7r₀ − r_{1/2} is a form in 1 and L(2,χ₋₃) alone (coefficient 108A₂), size ≈ 6r₀ — but r_{1/2} has no zero-cancellation on its lattice (partial sums over 6l+5 up to 6h₀, exponent q): height 28–38 per n, NET +21 → +29. (iii) Double roots (numerator (t+1)_{h₀−1}(t+1/2)_{h₀}, centre factor dropped so ε = −1 and both untwisted sums are interesting and clean): (3;1,1,1) then has degree −5 (power-law terms, NET +9 → +13); (4;1,1,1): r₀ +13 → +17, combination +20 → +23. All three routes to a pure L(2,χ₋₃) form (or to L(4), ζ(5), …) are dead by 20–30 per n.

## The cubic arithmetic: where the tax sits and how χ₋₃(p) enters (thirds_ledger_3111.json, n = 2…40)
Exponent e of p in D/g for (3;1,1,1) as a function of x = n/p (p ≥ 11; class = p mod 3): 0 for x < 1/3 — except e = 3 exactly at n = (p−1)/3 for p ≡ 1 mod 3 (observed at p = 7, 13, 19, 31, 37, 43, 61: the cross-coset difference (3d+1)/3 with d = n is p/3, so each of the three blocks of the other coset contributes p⁻¹ at the pole l = 0 and nothing compensates); 1 on (1/3, 2/3); 3 on (2/3, 1); 6 on (1, 4/3); 7 on (4/3, 2) (class 1: 6 at x = 1.69, 1.85); 9 on (2, 3) (class 2: 8 at x = 2.00–2.24); 12 on (3, 10/3). In the normalised scale ẽ = e − v_p(N_n): 2, 5 | 3, 5, 6 | 3, 4, 5 | 3, … on the thirds of the unit intervals (1/3,1) | (1,2) | (2,3) | (3,…) — bounded, so the normalised height converges (∫ẽ dx/x² ≈ 10 from these steps, plus p = 2 (ẽ₂ = 9 at n = 8) and higher levels; measured 11.55 at n = 30).
Level-1 formula for the top coefficient at the 1/3-pole −(k+1/3), k = n+l, with x = n/p, y = l/p and ρ₁ = (3⁻¹ mod p)/p, which is (p+1)/(3p) ≈ 1/3 for p ≡ 2 mod 3 and (2p+1)/(3p) ≈ 2/3 for p ≡ 1 mod 3: v_p(a_{3,k}) = [⌊2x−y−ρ₁⌋ + ⌊x+y+ρ₁⌋ + 1] − 3[⌊x−y+ρ₁⌋ + ⌊y−ρ₁⌋ + 1] − 3[⌊y⌋ + ⌊x−y⌋] (numerator AP window of length 3n minus three reciprocal AP windows of the other coset minus the same-coset factorials). At x ∈ (1/3, 1) its minimum over y is −1 (class 2: poles with y ∈ (1/3, x)): each cross window holds one multiple of p that the numerator window does not cover — the AP defect. The lower orders add up to 2 more on (2/3, 1) (reciprocals of the cross factors 3d+1 = p) and nothing on (1/3, 2/3); the partial sums add nothing for p > 3n (complete cancellation, the analogue of the RZ/Lai–Zhou contradiction argument). So the character enters the carries exactly through 3⁻¹ mod p, i.e. through χ₋₃(p): the thirds lattice reads a prime by its residue mod 3. Full proof of ẽ(x; χ₋₃(p)) not written (the cubic Lemma B); the data above is the specification.

## What the 3-fold does, in one paragraph
Two cosets give the character without a twist and double the denominator degree for one pole order; the palindrome exchanges them, so the symmetric part (ζ at odd orders) survives next to the antisymmetric part (χ₋₃ at even orders) — all three residues mod 3 are present in every form: the integers carry the sum, the two other cosets carry ζ and L. The price is paid in the arithmetic: the reciprocal blocks of the other coset are arithmetic progressions with step 3, their p-adic content is aligned by 3⁻¹ mod p, and the misalignment (one multiple per window per level) is a denominator of the order D_{3M}^q that the half lattice never pays (all its poles are one coset). At the smallest parameters the trade is nearly a wash (net ≈ +0.3 for thirds against +5.07 for halves at (3;1,1,1)); at every larger parameter the tax wins. Closed as a route to a theorem; kept as the cleanest instance so far of "the base reads the prime" in the carries.

## Files
cubic.py, cubic_asym.py (type B, exact and asymptotic), thirds.py (type A: forms / peak / combo / ledger; flags double, cube, quarter), thirds_ledger_3111.json (per-prime exponents to n = 40); tmp: cancel_test.py, zud_peak.py, ledger_analyze.py.

# The 9 (2026-09-29): every placement of 3² in the linear-form machine, measured
The user's correction: "3 cubed and ζ(3) are being taken too literally; the simple integer 9 = 3² makes more sense in placement — placement and factorization unknown." What follows is the systematic answer: where a 9 can sit in a Zudilin-type form, what each placement does to the values and to the carries, and the numbers.

## Where 9 already was
- K² arena: the analytic side U(π) = −ln 9 is E(1) with E(n) = ((n+3)/2)ln(n+2) − ((n−1)/2)ln n at n = ρ/π = 1: 9 = (1+2)^{(1+3)/2}. The square is the exponent (n+3)/2 = 2, the 3 is n+2. Understood, not a base-9 fact.
- Runs 4, 6, 7: base 9 is carry-proof for χ₋₄ (9 ≡ 1 mod 8: n ≡ digit sum mod 8), so (−1)^ν = (−1)^{s₉(ν)} exactly (Run 6: zero carry cost at z = −1); primes sorted by base-9 last digit × Pythagorean class = p mod 36 (Run 7).
- Thirds forms: the L(2,χ₋₃) coefficient is 9A₂ = 3²A₂ (Σ_ν (ν+k+a/3)^{−2} = 3²Σ(3ν+3k+a)^{−2}); ζ(3) enters with 26A₃ = (3³−1)A₃. So 9 is the weight of the second moment (conductor²), and 27 only appears as 3³ − 1 in the symmetric (ζ) part — that is what "too literal" meant: the cube in the ζ(3) coefficient and the cubed numerator (the "cube" flag) are devices, not structure.

## The general machine: coset.py (poles on cosets a/b of (1/b)Z, sum over ν ≥ 0 with a rational twist of period T)
Σ_{ν≥0} tw(ν)(ν + k + a/b)^{−i} = Σ_{r<T} tw(r)T^{−i}[ζ(i, u/(bT)) − Σ_{m<k″}(m + u/(bT))^{−i}] with (r+k)b + a = k″bT + u: the form lives at conductor bT. Catalan (Zudilin 2018) is b = 2, A = {1}, T = 2, tw = (1,−1): conductor 4 = 2·2. Thirds type A is b = 3, A = {1,2}, T = 1: conductor 3. So 9 has two factorizations here: 9·1 (poles at ±1/9) and 3·3 (thirds poles × sublattice-3 sum) — plus the Catalan-side placement 4·9 = 36 (the alternating Catalan form split by ν mod 9). Controls reproduced exactly: b = 2 gives Zudilin's form (r = 176.3156…/4096 at (3;1,1,1), n = 2, the γ_n normalization), b = 3 gives thirds.py's numbers to all digits. The centre factor (2t+h₀) reduces the order of the centre pole −h₀/2 on the half lattice (the pole coincides with the numerator zero — the "midpoint requires a square"); on odd lattices there is no centre pole, but the factor is still what makes ε = −1 (without it the untwisted sum is a π-identity and the alternating sum has only π-values: checked).
Interesting-sum criterion with twists: R(−t−h₀) = εR(t), tw(−ν−h₀) = σ tw(ν) ⇒ Σ_{ν∈Z} = (1+εσ)Σ_{ν≥0}; a period-T twist must be (anti)symmetric under ν ↦ −ν−h₀ mod T. For T = 3 exactly one sublattice class is self-symmetric, c* ≡ h₀ (mod 3) (for η₀ = 3, c* = 1 for every n); S_{c*} and the full sum r are the two genuine forms, S_0 ≡ S_2 modulo a π-identity.
Level 1 appears in every sublattice form: the simple-pole coefficients C_{1,u} no longer cancel symbol by symbol (only Σ_u C_{1,u} = 0), so ζ(1,u/N) is replaced by its regularised value −ψ(u/N) (form check then passes to 10⁻⁴⁸): by Gauss's digamma theorem and the symmetric pairing u ↔ N−u (cot parts cancel, γ and ln 2N cancel because Σ C = 0) the level-1 part is Σ_k w_k ln sin(πk/N) with w_k in the real cyclotomic field — for N = 9 the cyclic cubic field Q(cos 2π/9) of conductor 9, discriminant 81 = 9², units sin(2π/9)/sin(π/9) etc. This is the one place where 9 = 3² is genuinely the object: the 3·3 factorization drops a hyperoperation level (Σ 1/n² → logs), and the level-1 constants are the units of the conductor-9 cubic field.

## Measurements at (3;1,1,1), NET = log|r·D/g|/n (normalization-free), n = 2, 4, 6, 8
| placement | conductor | values | NET |
|---|---|---|---|
| thirds, untwisted (control) | 3 | 1, L(2,χ₋₃), ζ(3) | −2.55, −2.46, −1.29, −1.57 |
| poles ±1/9 (9·1) | 9 | 1, X₂ = Σ_{n≡±1(9)}±n⁻², Y₃ = Σ_{n≡±1(9)}n⁻³ | +6.92, +7.30, +10.86, +8.94 |
| thirds × sublattice c* = 1 (3·3) | 9 | 1, logs of Q(ζ₉)⁺-units, three antisymmetric ninths at i=2, three symmetric at i=3 | −1.32, −0.53, +0.45, +1.05 |
| thirds × sublattice 0 | 9 | same symbols, no palindromic pairing | −0.05, −0.08, +1.25, +1.32 |
| half lattice, alternating (Zudilin, control) | 4 | 1, β(2) | +1.89 (n=2) … +3.22 (n=8) |
| half lattice × base-9 class d* (4·9) | 36 | 1, logs of Q(ζ₃₆)⁺-units, β(2), L(2,χ₋₃), L(2,χ₋₄ρ) (ρ cubic mod 9), ζ(3)-type at i=3 | +4.62, +3.78, +3.53, +3.47 |
Normalised Stirling peaks: ±1/9 −14.858 = thirds −11.562 − 3 ln 3 (pure normalization: same decay), all six ninths −57.9 (degree 6q; not run exactly — the height scales the same way), (5;2⁵) at ±1/9 −29.99. Poles at ±1/9 lose ≈ 12 per n against thirds: the Hurwitz tails Σ_{m<k}(m+1/9)^{−i} = 9^iΣ(9m+1)^{−i} reach 9k ≈ 9M (D by range at n = 8: 3M<p≤6M 4.1, p>6M 3.2 per n; primes 43, 61, 79 to the 3rd power), and the entering exponents are 5 where thirds has 1–3.

## The base reads the prime: the entry law (exact, 43/43 primes)
For two cosets ±a/b, the cross-coset factors at a pole are (b·d + (b−2a))/b, d = k′−k ∈ [−n, n] (K = n+1 blocks). A prime p first enters the denominator when b·d + (b−2a) = ±jp has a solution with |d| ≤ n, j the least positive integer with jp ≡ ±(b−2a) (mod b): x_entry = j/b, boundary exactly at p = (bn ± (b−2a))/j.
- b = 3, a = 1 (thirds): b−2a = 1, j = 1 for every p: entry at x = 1/3 for all classes; sign + for p ≡ 1 (mod 3) at n = (p−1)/3 (the "sporadic e = 3" of the earlier ledger, now explained: the partner pole sits exactly at the window edge, exponent q = 3 uncompensated), sign − for p ≡ 2 at n = (p+1)/3 with exponent 1. So thirds read χ₋₃(p).
- b = 9, a = 1: b−2a = 7: j = 1 for p ≡ 2, 7 (mod 9), j = 2 for p ≡ 1, 8, j = 4 for p ≡ 4, 5 — the coset of p in (Z/9)*/{±1} ≅ Z/3, i.e. the cubic character mod 9. Verified against coset_ledger_9_18_3111.json: predicted first n = observed first n for all 43 primes 11 ≤ p ≤ 367 (entry_law.py). Entering exponent 3 for j = 1 (the partner pole is at distance exactly p/9 and the pair's leading parts cancel two powers) and 5 = 2q−1 for j = 2, 4 (full derivative saturation). Beyond entry the exponent function e(x; p mod 9) differs class by class (ledger9_analyze.py): class {1,8}: 5 on (2/9, 0.65), 4, 3, 6 at x ≈ 1, 7 on (1.26, 1.9); class {4,5}: 5 from 4/9, then 10–14 beyond x = 1.4 for digit 4; class {2,7}: 3 on (1/9, 0.43), 1 on (0.46, 0.64). So the ±1/9 ledger reads the base-9 last digit of p (Run 7's variable), through 9⁻¹ mod p — and it reads it expensively.
- b = 2 (Catalan): one coset, no cross factors; the ledger is read by x = n/p and the base-p digits only (Zudilin's φ), never by p mod 4.

## The placement that touches Catalan: 4·9 = 36
Because (−1)^ν = (−1)^{s₉(ν)}, the alternating Catalan form is the sum of nine class sums S_d = Σ_{ν≡d(9)}(−1)^νR(ν); the self-symmetric class is d* ≡ −5h₀ (mod 9) (d* = 1, 7, 4 cycling with n for η₀ = 3). Each S_d is a form at conductor 36 = Run 7's modulus: at even i it contains β(i) (u = 9, 27), L(i,χ₋₃) (u = 3, 15, 21, 33 → conductor 12), and the conductor-36 odd characters χ₋₄ρ, χ₋₄ρ̄ (ρ the cubic character mod 9); at odd i the even characters (ζ-like unknowns); at level 1 the logs of Q(ζ₃₆)⁺-units. Ledgers (cat_full_ledger vs cat9_ledger, n = 10…40 at (3;1,1,1), n = 12…32 at (5;2⁵)): every top-range prime enters the class sum with exponent s − 2 instead of s (3 → 1 at s = 3: p = 17, 19, 23, 29, 31, 37; 5 → 3 at s = 5: p = 17, 19, 23, 29, 31), the odd-prime height drops from 1.38 to 0.99 per n (s = 3, n = 40) and from 2.16 to 1.84 (s = 5, n = 32), the 2-adic part rises by 0.09–0.17 per n (2^{6n} vs 2^{≈5.8n}), 3-adic unchanged: net height gain 0.3–0.6 per n (s = 3), 0.06–0.6 (s = 5). The mechanism is the K² budget's "partial-sum line": the class sum's Hurwitz tails are progressions 36m + u, so a prime in (M/2, 2M] appears in a tail only if jp ≡ u (mod 36) has a solution with j ≤ 2M/p — the base decides which primes pay. The size does not change for exponentially decaying η (no twist cancellation: for Zudilin's own η, log R_n(0)/n = −100.745 and the interior peak −100.735 bracket the true closeness −100.740), so for such η the class sum would be a genuine 0.3–0.6 per n better on the ledger — at the price of five new constants per even order and level-1 logs. At the power-law parameters measured here the class sum is larger than the total (weaker alternation), so NET is worse (+3.47 vs +3.22 at n = 8). Not a Catalan tool; it is the exact sense in which base 9 is Catalan's base: it makes the tails sparse and the tail tax a function of p mod 36.

## Verdict on the hint
9 = 3² has three placements in this machine and each was measured. As a pole lattice (9·1) it is the conductor-9 analogue of Catalan's 4 = 2² and reads primes by their base-9 last digit (cubic character mod 9), but costs 12 per n. As lattice × sublattice (3·3) it is the conductor of the cyclic cubic field Q(ζ₉)⁺ and drops the form to level 1 (logs of its units) — the only place where the square is the object itself — at +1 to +2.6 per n. As Catalan's carry-proof base (4·9 = 36) it thins the Hurwitz tails and cuts the top-range tail tax by exactly two powers of every prime, at the price of conductor-36 constants. None yields a theorem; the second and third are the ones to remember: the base of the lattice is the base in which the carries read the prime, and the square of the base is where the level drops.

## Files
coset.py (general machine: forms / peak / ledger; controls b = 2 T = 2 and b = 3 reproduce multibeta.py and thirds.py), coset_ledger_9_18_3111.json (±1/9 ledger to n = 40), cat9_ledger_3111.json, cat_full_ledger_3111.json, cat9_ledger_522222.json, cat_full_ledger_522222.json (Catalan full vs base-9 class sums); tmp: cat9.py (class-sum driver), entry_law.py, ledger9_analyze.py, compare_ledgers.py.

## Overlay of the two dimensions: tenths, sixths, twelfths, and the odd/even law (2026-09-29, later)
The user asked whether the key is "dimension 2 and 3 overlaid by base 9 and base 10, cross factorization, with 5i in the middle", and added "3² / μ(9) = 0". Measured with coset.py at (3;1,1,1), n = 2, 4, 6, 8:
| lattice | constants at i = 2 | entry law reads | NET |
|---|---|---|---|
| ±1/10 (base 10) | 100·Re[(1 − i/4)L(2,ψ₅)], ψ₅ the quartic character mod 5, values {±1, ±i} | (p/5): last decimal digit 1, 9 enter at x = 1/5 with exponent 5; digits 3, 7 at 2/5 with exponent 1 (23/24 primes exact; p = 37 enters late) | +4.26, +7.62, +8.73, +9.44 |
| ±1/6 | 45·L(2,χ₋₃) (thirds' constant with the Euler factor 1 + 2⁻²) | every p at x = 1/3 with exponent 5 (j = 2 for all) | +2.99, +5.29, +7.15, +8.33 |
| ±1/12 | ζ(2,1/12) − ζ(2,11/12) = 80·G + 90·L(2,χ₋₃) (one constant: (12²/2)[(1+3⁻²)β(2) + (1+2⁻²)L(2,χ₋₃)]); at i = 3: 864[(7/8)(26/27)ζ(3) + L(3,χ₁₂)] | χ₁₂(p) = χ₋₄(p)χ₋₃(p): p ≡ ±1 (12) enter with exponent 5, p ≡ ±5 with exponent 3, both at x = 1/6 (all 31 primes to 149) | +9.37, +11.24, +14.31, +14.74 |
| {1,5,7,11}/12 (four cosets) | β(2) and L(2,χ₋₃) separately (pairs (1,11), (5,7)); level-1 logs of the units of Q(√3) (ln(2+√3) — Run 3's "partner" constant) | every p with exponent 5 | +4.43, +8.48, +9.51, +11.22 |
(thirds control −2.55 → −1.57; half-lattice Catalan control +1.89 → +3.22.)
The law behind all of it: the values of a two-coset form are the ODD characters of the conductor (β from χ₋₄, L(χ₋₃) from χ₋₃, the quartic ψ₅ from 10, the sextic from 9), while the arithmetic — the entry point and entering exponent of a prime — is read by the EVEN characters of the conductor through the least j with jp ≡ ±(b−2a) (mod b): the cubic character for 9, (p/5) for 10, χ₁₂ = χ₋₄χ₋₃ for 12, nothing for 3, 4, 6 (which have no even character; thirds read χ₋₃ only through the sign of the boundary n = (p∓1)/3). The "cross" is exact for the overlay: the twelfths form's constant is the sum of dimension 2 and dimension 3 with weights 80 = 8·10 and 90 = 9·10, and its ledger reads their product χ₁₂. The overlay is the most expensive lattice tried: each extra factor in the conductor multiplies the tail range (tails to bM) and raises the entering exponents to 2q − 1 = 5; the smallest lattice that carries a character is always the cheapest (b = 2 for χ₋₄, b = 3 for χ₋₃).
μ(9) = 0, read on the ledger: split log(D/g) over odd primes p > 3 into the squarefree kernel Σ_{e_p ≥ 1} ln p (what Möbius sees) and the square part Σ (e_p − 1) ln p. Catalan full form, n = 40: kernel 0.696, square part 0.685 per n; its base-9 class sum: kernel 0.696 (identical), square part 0.298. The class sum changes nothing that Möbius sees; its entire gain is one square (p²) removed from every top-range prime — which is exactly the "3 → 1, 5 → 3" law above. The thirds ledger at n = 40 is kernel 2.42, squares 7.79 (raw, normalization included). Verdict on the overlay: measured and dead as a form; alive as the statement "odd characters are the constants, even characters are the carries, and a base-b class sum removes one square".
The wrap, in the user's terms (2026-09-29, later). The user's picture — "3 − 2 = 1 and 10 − 9 = 1; base 9 is relative to 2D, base 10 to 3D; a square in the middle, wrapped together is i; four parts squared in the middle with the midpoint ½ of 10, giving 5i" — is the conductor-12 lattice, i.e. Q(ζ₁₂) = Q(i, ω) (Run 5's three characters), made exact: 9 = 2³ + 1 and 10 = 3² + 1, so each base reads the other prime's power by digit sums (base 9 reads mod 8 ⊃ conductor 4 of χ₋₄; base 10 reads mod 9 = conductor of the cubic side); 9 = 2³ + 1³ factors through the Eisenstein norm a² − ab + b² = 3·3, 10 = 3² + 1² through the Gaussian norm = (3+i)(3−i), 5 = 2² + 1² = (2+i)(2−i), and (2+i)(1+2i) = 5i. In the twelfths form the single even-order constant is, verified to 30 digits,
  ζ(2,1/12) − ζ(2,11/12) = 10·(8·G + 9·L(2,χ₋₃)) = 10·(2³·G + 3²·L(2,χ₋₃)),
the 8 = 9 − 1 on the 2D constant, the 9 = 10 − 1 on the 3D constant, the 10 in front; the two weights come from the Euler factors at the inert primes, (1 + 3⁻²) = 10/9 = N(3+i)/3² and (1 + 2⁻²) = 5/4 = N(2+i)/2² (3 is inert in Z[i], 2 is inert in Z[ω]; 1 + p⁻² = (p² + 1)/p² is "the square plus one, wrapped by i"), times 12²/2 = 72 = 8·9. Order 4: ζ(4,1/12) − ζ(4,11/12) = 128·82·β(4) + 648·17·L(4,χ₋₃), 82 = 3⁴ + 1, 17 = 2⁴ + 1 (verified). Order 3: ζ(3,1/12) + ζ(3,11/12) = 864[(7/8)(26/27)ζ(3) + L(3,χ₁₂)]. The "four parts" are (Z/12)* = {1, 5, 7, 11} ≅ (Z/2)², the middle pair 5, 7 = 6 ∓ 1 around the midpoint 6 = 12/2. Mod 36 is not part of this picture: it was the modulus of "base 9 reading the Catalan form" (conductor 4 × 9), which is Run 7's modulus. Cost of the wrap: the twelfths ledger reads χ₁₂ = χ₋₄χ₋₃ with entering exponents 5 (p ≡ ±1 mod 12) and 3 (p ≡ ±5), NET +9 to +15 per n — wrapping two characters into one pole lattice multiplies the conductor and the tail range; Zudilin's way of carrying a character for free is the twist on the smallest lattice, and a twist can carry only one odd character at a time (χ₁₂ itself is even, so its twist gives π-values at even order).

## The ladder from the pole (2026-09-29, later): "compare ∞ to 10, ∞−1 to 9, ∞−2 to 8"
Exact version: put the pole ζ(1) at the top ("10") and step down: rung m is ζ(1−m). Then (von Staudt–Clausen with the m-part; checked m = 1…12, tmp/rungs.py)
  den ζ(1−m) = lcm{N ≥ 1 : φ(N) | m} = ∏_{(p−1)|m} p^{1+v_p(m)},
the conductor of the compositum of every cyclotomic field whose degree divides m. So each rung down wraps exactly the cyclotomic fields of that degree; odd rungs wrap nothing (only N = 1, 2) and vanish.
| digit 10−m | m | ζ(1−m) | conductor | what enters |
|---|---|---|---|---|
| 10 | 0 | pole | — | the carry |
| 9 | 1 | −1/2 | 2 | the half |
| 8 | 2 | −1/12 | 12 = 4·3 | Q(i) and Q(ω): the 2D–3D wrap of the twelfths identity 10(8G + 9L(2,χ₋₃)) |
| 6 | 4 | 1/120 | 8·3·5 | Q(ζ₅) (quartic, values ±i: the "5i"), Q(ζ₈) |
| 4 | 6 | −1/252 | 4·9·7 | Q(ζ₉) (the cubic/sextic side): 4·9 = 36 first appears here |
| 2 | 8 | 1/240 | 16·3·5 | Q(ζ₁₆) |
| 0 | 10 | −1/132 | 4·3·11 | Q(ζ₁₁) |
Reading: the forms pay D_M = lcm(1…M) — every N up to M — while the rungs pay only the N whose cyclotomic degree divides m; the base-b class sums (one square removed per top-range prime) are a step from the first kind of denominator toward the second. Dictionary, not a construction.

# Ramanujan's weight cannot break π out (2026-09-29, later): the locking theorem in data
The user's aim: use the ladder/wrap to separate G from π ln(2+√3) in Ramanujan's G = (π/8)ln(2+√3) + (3/8)Σ_{n≥0} 1/((2n+1)²C(2n,n)) — the "partner sequence problem" of Run 3 / the K² target spec.
## Where the term comes from (polygon table, tmp/polygons.py)
S₂(x) := Σ_t (2x)^{2t}/((2t+1)²C(2t,t)) = (1/x)[φ ln tan(φ/2) + Cl₂(φ) + Cl₂(π−φ)], φ = arcsin x. The boundary term is the angle times the log of tan(φ/2):
| x = sin φ | φ | (2x)²/4 per term | Clausen part | boundary φ ln tan(φ/2) | tan(φ/2) |
|---|---|---|---|---|---|
| 1 | π/2 | 1 (no decay) | 2G | 0 | 1 (Q(i) has no unit) |
| √3/2 | π/3 | 3/4 | (5√3/4)L(2,χ₋₃) | −(π/6) ln 3 | 1/√3 (not a unit) |
| 1/√2 | π/4 | 1/2 | conductor 8 | (π/4)ln(√2−1) | unit of Q(√2) |
| 1/2 | π/6 | 1/4 | 4G/3 (the χ₋₃ part cancels in the odd-n sum) | (π/6)ln(2−√3) | unit of Q(√3) = Q(ζ₁₂)⁺ |
So G alone appears only at the edge x = 1 (t^{−3/2} convergence, no closeness); inside the disk the first G-carrying point is the 12-gon angle, and it carries the regulator of Q(√3): π ln(2+√3) = 4L(1,χ₋₄)·√3L(1,χ₁₂) = 9L(1,χ₋₃)L(1,χ₁₂) (verified 60 digits), the product of the level-1 values of the two characters whose product is χ₋₄. In Run 11's language: π/4 = L(1,χ₋₄) sits over the rational L(0,χ₋₄) = 1/2, and ln(2+√3) = −L′(0,χ₁₂) is the slope at the zero of the even character. Ramanujan's identity is level 1 × level 1 = level 2 + binomial sum.
## Locking (tmp/ram_basis3.py, ram_basis4.py; PSLQ at 100 digits, residuals 10⁻⁸⁸…10⁻¹⁰¹)
For the Ramanujan weight w_t = 1/C(2t,t) and half-integer poles, T(m,k) := Σ_t w_t/(2t+2m+1)^k:
- k = 1: T(m,1) ∈ Q + Q·π/√3 for m = 0…4 (3T(0,1) = 2π/√3; 3T(1,1) = −24 + 14π/√3; 9T(2,1) = −400 + 222π/√3; 75T(3,1) = −16072 + 8870π/√3; 3675T(4,1) = −3602528 + 1986530π/√3).
- k = 2: T(m,2) ∈ Q + Q·π/√3 + Q·R with R := 8G − π ln(2+√3), for m = 0…6: 3T(0,2) = R; 3T(1,2) = −36 + 6π/√3 + 8R; 27T(2,2) = −1844 + 342π/√3 + 384R; 375T(3,2) = −125684 + 24270π/√3 + 25600R; 128625T(4,2) = −199336092 + 39295410π/√3 + 40140800R; 6251175T(5,2) = …; 2773437975T(6,2) = ….
- k = 3: T(m,3) ∈ Q + Qπ/√3 + QR + Q·R₃ with R₃ := T(0,3) = 1.02002080065…, for m = 1…4 (3T(1,3) = −48 + 6π/√3 + 4R + 24R₃; 27T(2,3) = −2488 + 330π/√3 + 224R + 1152R₃; …).
So every linear form Σ_t Q(t)/C(2t,t) with poles at half-integers of order ≤ k lives in the tower {1, π/√3, R, R₃, …}: exactly one new constant per level, as in Zudilin's β tower, and the level-2 constant is Ramanujan's number R = 8G − π ln(2+√3), never G. G and π ln(2+√3) are locked in the ratio 8 : −1 by the weight itself (the shifts reduce to derivatives of S₂ at x = 1/2, which are level-1 and rational). Consequences: (i) no Ramanujan-weighted construction can separate them — an irrationality result from this weight would be about R, i.e. "one of G, π ln(2+√3) is irrational" at best (π ln(2+√3) alone is, as far as found, not known irrational, so that would be new but says nothing about G); (ii) G-pure binomial series (Lupaş 2000; Lima, arXiv 1207.3139, with (4n)!-type weights) carry the conductor 4 inside the weight (quarter-integer factorials), not in the evaluation point — that is the only way G appears alone with geometric decay, and it is the very-well-poised world of Zudilin's u_nG − v_n with its 2^{4n}D_{2n}² ledger (K¹ arena, miss 4.37 per n). Integer poles with this weight give the 3-gon constants (π², √3L(2,χ₋₃)) plus level-1 logs (ln 2, ln 3 — not identified, spurious PSLQ at 100 digits; not needed).
## Literature found in passing (2026-09-29)
- **Zhi-Wei Sun, "Catalan's constant is irrational", arXiv 2609.04176 (3 Sep 2026, math.GM, 20 pp).** A full-proof claim in this program's K² arena: weighted tails u_m = T_m/(2m+1), T_m = Σ_r(−1)^r/(2m+2r+1)² (a Stieltjes moment sequence), pole factors Π_i = ∏_{h≤B}(2(h+i)+1)², binomial finite differences of order a+2B (i.e. an external field (1−z)^{2B}), an S×S residual determinant with S = B/20, Cauchy–Binet into a Pascal alternant (Vandermonde) × Cauchy determinant, a p-adic ledger in the three ranges Q ≤ S, S < p < B, p > B, and a claimed net −δ₀B² with δ₀ = 0.00966: raw archimedean 39/200 = 4ρ − 2ρ² (ρ = 1/20) against p-adic gains −c_odd + Λ_mid + 83/2400 = −0.00628 + 0.17636 + 0.03458 = 0.20466 — a 5% cancellation. Acknowledgments: the proof was developed in conversation with AI, the numerical data were produced by AI, and "the whole proof has passed the verification of Chatgpt 5.6 Solar". Category math.GM (moderator reclassification is likely). Cites Nesterenko 2016 "On Catalan's constant", Rivoal 2006, RZ03, Zudilin 2019. Not checked here; it is exactly the kind of object this program's closed-form K² theory (Heine + potential theory for the size, per-prime exponent laws for the height) is built to audit, and the audit is the natural next task (see handoff).
- Eskandari–Murty–Nemoto, "Mixed motives and linear forms in the Catalan constant" (arXiv 2510.20648, 2025/26): a 2-dimensional mixed motive with period G, a supply of linear forms in 1 and G with explicit coefficients; Eskandari, "On two families of period integrals related to Catalan's constant" (arXiv 2609.26308, 22 Sep 2026): the boundary family is Rivoal–Zudilin–Nesterenko's hypergeometric family.
- Krattenthaler–Rivoal, "On a linear form for Catalan's constant" (arXiv 0810.1927): Andrews' multidimensional Watson transformation proves Rivoal's conjecture on the arithmetic of the coefficients of a very-well-poised linear form in 1 and G (denominators, not irrationality).
- Lima, "A rapidly converging Ramanujan-type series for Catalan's constant" (arXiv 1207.3139): G-pure central binomial series faster than Ramanujan's and Lupaş's, from a hypergeometric identity; the author suggests an Apéry-like proof might follow.
- Sun's paper, per the user (2026-09-29): already read; one or two undefined quantities (the parameter D of §3 is indeed never defined in the text) and one major omission; not to be treated as established.

# Ramanujan summation and Catalan's ladder (2026-09-29, later): what "breaking π out" means exactly
The user meant Ramanujan summation (the −1/12 mechanism), not Ramanujan's series and not the K² arena. Made exact (tmp/ramsum.py):
- Catalan's ladder is the Euler numbers: β(1−m) for m = 1…12 is 1/2, 0, −1/2, 0, 5/2, 0, −61/2, 0, 1385/2, 0, −50521/2, 0 — β(−2k) = E_{2k}/2, the odd negative arguments are trivial zeros. The only denominator is 2: unlike ζ's rungs (lcm{N : φ(N) | m}, the von Staudt ladder), the χ₋₄ ladder is von-Staudt-free — the "Genocchi-type moments" of Run 10 — and the prime 2 is the whole arithmetic of the negative side, as the 2-adic tax is on the positive side.
- G sits at the wrong parity for the functional equation: β(2) is paired with the trivial zero β(−1) = 0, so π is not "rationalized" as for ζ(2) = −2π²ζ(−1); instead G/π is the slope: β′(−1) = 2G/π (verified 30 digits; Run 11). Ramanujan/Abel summation realises it as the regularised divergent sum −Σ_{k≥0}(−1)^k(2k+1)ln(2k+1) = 2G/π, and the Euler transform (Hasse–Sondow globally convergent series for β(s), differentiated at s = −1) gives the convergent base-2 expansion
    2G/π = −Σ_{n≥0} 2^{−n−1} ln r_n,   r_n = ∏_{k≤n}(2k+1)^{(−1)^k C(n,k)(2k+1)} ∈ Q  (r₁ = 1/27, r₂ = 3125/729, r₃ = 5¹⁵/(3⁹7⁷), …),
  verified: partial sums reach 2G/π to 7·10⁻²¹ at n = 60; the digits ln r_n decay like 1/n, the weights like 2⁻ⁿ. So π is broken out in this sense: G/π is a level-1 object, a dyadic expansion whose digits are logarithms of rationals built from the odd numbers ≤ 2n+1 (the χ₋₄ lattice) with binomial exponents (the carries). Probably known in the Sondow–Hadjicostas circle (their ln(4/π) is the same construction at s = 0: γ(−1) = ln(π/4)); not found stated for 2G/π.
- What it does not give: an irrationality route. The truncations are linear forms in ~N/ln N logarithms of odd primes with denominators 2^{N+1} and error 2^{−N}/N; Baker–Wüstholz lower bounds are exponentially weak in the number of logarithms, and there is no approximation sequence with a single unknown. Ramanujan summation also cannot extend the Zudilin-type forms to rational functions of degree ≥ −1: the polynomial part sums to Euler numbers times huge coefficients and the smallness is lost. Verdict: a true and pretty dictionary — the descent to the negative rungs converts level 2 to level 1 at the price of unboundedly many logarithms — not a construction.
## Adding infinity back (2026-09-29, later; tmp/infinity.py) — the divergent side made exact
The user: "infinity can be added to both sides … maybe i can be added … play with −1/12 and +1/6". Three verified identities.
1. Boole (the infinity itself). With f(k) = (2k+1)ln(2k+1) and A_N = Σ_{k≤N}(−1)^k f(k): A_N = −2G/π + (−1)^N D_N, D_N = f(N)/2 − (E₁(0)/2)f′(N) − (E₃(0)/12)f‴(N) − …, Euler-polynomial constants E_j(0) = 2(1−2^{j+1})B_{j+1}/(j+1): E₁(0) = −3B₂ = −1/2, E₃(0) = 1/4, E₅(0) = −1/2, E₇(0) = 17/8. Explicitly D_N = (N+1)ln(2N+1) + 1/2 + 1/(6(2N+1)²) + O(N⁻⁴); the residual A_N − (−1)^N D_N + 2G/π is 1.0·10⁻⁸ with two Boole terms and 1.6·10⁻¹⁵ with three at N = 2000. So the "infinity" is elementary: (N+1)ln(2N+1) + 1/2 with Bernoulli-number corrections; B₂ = 1/6 is the coefficient that converts the half-term f(N)/2 into (N+1)ln(2N+1) + 1/2. It carries no arithmetic; all content is in the constant.
2. Adding i. Li_s(i) = −2^{−s}η(s) + iβ(s), so the s-derivative at the rung −1 is d/ds Li_s(i)|_{s=−1} = (7/6)ln 2 + 1/2 − 6 ln A + i·(2G/π) (A = Glaisher; verified to 3·10⁻⁴¹): the real part is ζ's slope ζ′(−1) = 1/12 − ln A at the same rung (this is where −1/12 enters, together with η(−1) = (1−2²)ζ(−1) = 1/4 in the (ln 2)/2 term), the imaginary part is Catalan's. In words: Σ^{Abel}_{m≥1} i^m m ln m = (6 ln A − 1/2 − (7/6)ln 2) − i·2G/π — Glaisher and Catalan are the real and imaginary parts of one regularised Gaussian sum, Run 11's pairing made literal. The Hasse–Sondow series at z = i is an expansion in the Gaussian base (1−i)/2 (norm 1/2) with digits ln ρ_n, ρ_n = ∏_{k≤n}(k+1)^{(−1)^kC(n,k)(k+1)} ∈ Q; it converges (6.6·10⁻¹⁰ after 200 terms) but slowly.
3. Multiplicative (Glaisher–Kinkelin analogue). H_N := ∏_{k≤N}(2k+1)^{(−1)^k(2k+1)} (the alternating odd hyperfactorial) satisfies, for even N, H_N = (2N+1)^{N+1}·e^{1/2 − 2G/π}·(1 + O(N⁻²)); i.e. √e·(2N+1)^{N+1}/H_N → e^{2G/π}, verified (1.08667415… = e^{2G/π−1/2} at N = 2000). e^{2G/π} is to the odd alternating hyperfactorial what Glaisher's A is to ∏k^k. As approximations these are Stirling-type: the height of H_N is Σ(2k+1)ln(2k+1) ≈ N² ln N against an error O(N⁻²) — algebraic convergence, no closeness per unit height — exactly as Stirling's formula does not prove √(2π) irrational. Playground recorded; nothing here is an irrationality route.

# The refined lemma on Zudilin's odd-zeta forms (2026-09-29/30; zetaforms.py)
Landscape (2026-09): Fauzan's "ζ(5) is irrational" (Zenodo, 17 Sep 2026; Lean formalization; vouched for on Calegari's blog; its Hankel/weight machine is the one this program's K² arena measured at net −0.126 per K²). So Zudilin's "one of ζ(5), ζ(7), ζ(9), ζ(11)" is vacuous and the live tower starts at ζ(7); the current record there is Lai–Zhou's "two of ζ(5)…ζ(35)" plus Fauzan, i.e. one of ζ(7)…ζ(35). A ζ(7) Zenodo paper exists (the user has checked one: architecture copied from Fauzan, problems).
## The construction (Zu04 §8) and the machine
R(t) = (h₀+2t)·∏_{j≤r}[(t+1)_{h_j−1}(t+1+h₀−h_j)_{h_j−1}/(h_j−1)!²]·∏_{j>r}[(h₀−2h_j)!/(t+h_j)_{1+h₀−2h_j}], h₀ = η₀n+2, h_j = η_jn+1, r odd, q ≥ r+4 odd, Σh_j ≤ h₀(q−r)/2; F = (1/(r−1)!)Σ_{t≥1−h₁}R^{(r−1)}(t) = Σ_{m odd, r+2≤m≤q−2}A_mζ(m) − A₀ (the (r−1)-th derivative kills ζ(3)…ζ(r); the numerator to the r-th power is the "cube" of thirds.py rediscovered). Integer poles k ∈ [h_{r+1}, h₀−h_{r+1}], order = number of covering denominator blocks, centre pole −h₀/2 reduced by the centre factor. zetaforms.py: `forms r ETA n` exact partial fractions (identity-checked), A_m, A₀ with Zudilin's partial sums Σ_{l≤k−h₁}, numeric check by high-precision differentiation (needs ~2×digits of the coefficients: at his parameters the coefficients are 10⁹⁷ and F is 10⁻⁹⁹), true denominators per prime vs his bound vs the refined predictor; `asym r ETA` his C₀ (saddle point, Lemma 20), κ, C₂ and the refined height; `search r q eta0s [ref] [seed=…]`.
Validation: his parameters (91; 27³, 29…38), r = 3, q = 13: C₀ = 227.58020 (his 227.58019641), C₂ = 226.24944266 (his digits), net −1.3308; exact forms at n = 1, 2 reproduce R and the parity vanishing; small case (8;1,1,1,2,2,2,2): numeric = exact to 10⁻¹⁰⁹.
## Theorem 1 transfers (finite-n predictor, refined_exponent)
Exact top valuation ν̂(k) = v_p(constant) + v_p(h₀−2k) + Σμ_m v_p(m−k) − Σ_{k'≠k}s_{k'}v_p(k'−k) (for p > √(2h₀) all factors are single-digit, so Zudilin's floor formula is exact up to the centre term — with one twist: his ν_{k,p} exceeds ν̂ by exactly one per denominator block NOT covering the pole, and his lemma is still right because he charges the full loss q−r at every pole; the two errors cancel); derivative saturation v_p(B_{i,k}) ≥ ν̂ − min(S−i, D_k) (D_k = p-divisible numerator factors with multiplicity; S−i if another pole is ≡ k mod p); partial-sum window: Σ_{l≤k−h₁}l^{−m} carries p^m iff k − h₁ ≥ p. Against the truth at his parameters: n = 1 exact at 7/7 single-digit primes; n = 2 exact at 11/13 (37, 67: truth one lower — the residue-drop shape); n = 3 exact at 21/21; n = 4 exact at 23/26 (41, 53, 59 one lower); no violations above √h₀ (below it the predictor is not claimed). Single-digit heights per n (Zudilin | refined | true): n=1 132.2 | 128.8 | 128.8; n=2 151.2 | 139.7 | 135.8; n=3 166.1 | 160.5 | 160.5; n=4 177.4 | 171.1 | 168.1 — the rebate fluctuates around the asymptotic 4.67.
## Asymptotics (refined_rebate; validated: integrating his exponent function over all x reproduces C₂ = 226.2494 exactly; the asymptotic refined exponent matches the finite predictor at 35/36 primes at n = 20)
| family | slots | best η (Zudilin's bound) | Zudilin net | refined rebate | refined net |
|---|---|---|---|---|---|
| r=3, q=13 | ζ(5,7,9,11) | his (91;27³,29…38) | −1.331 | 4.667 | −5.997 |
| r=3, q=13 | ζ(5,7,9,11) | (91;27,27,28,29…38) (search) | −3.083 | — | — |
| r=3, q=11 | ζ(5,7,9) | (50;15³,16,17,17,18,19,20,20,21) | +20.398 | 1.667 | +18.731 |
| r=3, q=11 | ζ(5,7,9) | (70;21³,22…29) | +24.881 | 1.333 | +23.547 |
| r=3, q=11 | ζ(5,7,9) | (90;27,27,28,29,…,37) | +29.823 | — | — |
| r=5, q=13 | ζ(7,9,11) | (60;17⁵,18,18,19,19,20,20,20,21) | +152.5 (C₀ = 46.2, C₂ = 198.8) | — | — |
Where the rebate sits for his direction: x = n/p ∈ [1/33, 1/31): 5 → 4; [3/91, 1/30): 4 → 2; [1/30, 1/29): 4 → 3; [1/27, 1/26): 8 → 7 — the same shape as in the β case (primes just below the largest block length), worth 4.67 per n against his margin 1.33 (a 4.5× margin, as the β case gave 30×).
## Verdict for the tower
The refined lemma transfers to the odd-zeta forms exactly as to β, but it moves nothing that matters now: ζ(5,7,9) is +18.7 per n even refined and rises with η₀ (the family needs ζ(11), as Zudilin found); towers starting at ζ(7) through the fourth derivative (r = 5) collapse the closeness (0.8 per unit η₀ against 2.5 for r = 3) and sit near +150 per n; and the {5,7,9,11} form, now at −6.0 per n, proves nothing new after Fauzan. The one direction left in this family is Lai–Zhou's two-value construction (rational roots at 1/2, 1/3, 2/3 and Vandermonde elimination), where exact denominators could lower their a = 17 and give "one of ζ(7), …, ζ(2a+1)" with a smaller a; not started. (The r = 5 closeness formula is Zudilin's Lemma 20 applied with r = 5, which he states for r = 3; an exact check on (16; 2⁵, 3⁴) — pure ζ(7) forms, saddle C₀ = −2.29, i.e. the forms grow — was running when this was written; the structural reason, the numerator to the fifth power dominating the summand, does not depend on it.)

# Lai–Zhou's two-value construction with exact denominators (2026-09-30; lz.py) — the live direction
Target: with Fauzan's ζ(5), Lai–Zhou's "two of ζ(5)…ζ(35)" is "one of ζ(7)…ζ(35)"; their elimination uses three forms sharing the coefficients ρ_i, Ŝ_b = Σ_{k≤b}S_{n,k/b} = ρ̂_{0,b} + Σ_i ρ_i b^i ζ(i), b = 1, 2, 3, combined with integer weights w ⊥ (1, 2^{i₁}, 3^{i₁}), (1, 2^{i₂}, 3^{i₂}). Fixing the eliminated pair at (3, 5) (both known irrational) gives w = (45, −9, 1) and S̃ = 45Ŝ₁ − 9Ŝ₂ + Ŝ₃ ∈ Q + Σ_{i odd ≥ 7} Q ζ(i): every s with net < 0 is the theorem "one of ζ(7), …, ζ(s)". Their s = 35: m₁ = 209, m₂ = 243, δ = (4⁵, 5, 6, …, 10, 12, 14, …, 52, 56, 60, 64, 68), C₁ = 16779.9312 (the Φ_n content rate — a gain), C₂ = 16779.2826 (D-part 8486 + log g(x₀) 8293.28), margin 0.65 per n.
## Machine (lz.py) and validations
`forms s m1 m2 deltas n`: exact partial fractions of R_n (four root blocks θ ∈ {1, 1/2, 1/3, 2/3}, the (t)_{m₂n+1} cancellation, integer poles from the δ-blocks, centre factor), ρ_i and ρ_{0,θ} with LZ's partial sums Σ_{ℓ≤k}(ℓ+θ)^{−i}, the (3,5)-combination, numeric check by direct summation at rational points, true denominators vs LZ's bound vs the refined predictor. Toys (s = 7, 9): identity ✓, even ρ vanish ✓, numeric = exact to 10⁻¹³⁵, refined predictor exact at 8/10 single-digit primes (other two: truth lower), including negative exponents (content). `asym`: x₀ = 2.894938334 (theirs 2.89493833), C₂ = 16779.282641 (their digits), C₁ = 16779.912 (theirs 16779.931; my Farey cap), net −0.629.
## The refined lemma here, and a 26-per-n rebate at their own parameters
Exact top valuation ν̂ = ν_LZ − (number of non-covering blocks) (same compensation as in Zudilin: they charge the full loss s+1−i at every pole); numerator content D(x,y) = v_p(G(−k)) = (⌊2A⌋−⌊A⌋+⌊2B⌋−⌊B⌋) + ⌊3A⌋ + ⌊3B⌋ − ⌊y⌋ − ⌊m₂x−y⌋ with A = m₁x+y, B = (m₁+m₂)x−y (the three shifted root blocks and the θ = 1 remainder); saturation; windows y ≥ 1 (θ = 1), 1/2, 1/3 (the θ ≠ 1 pair, symmetric in p mod 3); LZ's cross-pole cap for θ ≠ 1 (their Lemma 3.2: the numerator roots at θ-shifted points make the θ ≠ 1 tails cancel between poles — no such roots survive for θ = 1 in [−m₂n, 0], which is why only ρ_{0,1} carries the D_{M_j} beyond M). Asymptotic exponent functions validated against the finite predictor at n = 2 on their parameters (69/74 single-digit primes). Result at s = 35: rebate exactly 26.000 per n, all on p/n ∈ (235, 239) (LZ exponents 8, 7, 6, 5 → 0): for p > M = (m₂−2δ_min)n every coefficient is p-integral with content D ≈ 9 at the window poles k ≥ p−1, whose order S ≤ 8, so the θ = 1 partial-sum tax (their ∏_j D_{M_j} excess Σ_j (m₂−δ_j−M′)₊ = 26) is fully covered. This follows from their own Lemma 3.1 (its per-pole bound holds for all p > √(3L)); they apply Φ_n only up to p ≤ M. Closed form: the refined D-part is (s+1)·M′ instead of Σ_j max(M′, m₂−δ_j), i.e. refined net = (s+1)(m₂−2δ_min) + log g(x₀) − C₁ (+ a residual if D < S somewhere in the top range; zero here). Their −0.63 becomes −26.63 per n.
## Lowering s (2026-09-30/10-01): measured and closed for this family
Objective: refined net = (s+1)(m₂−2δ_min) + resid + log g(x₀) − C₁ (fast_C1 on a u = p/n grid, du = 0.01, nx = 4000; coarse values differ from the Farey ones by ≈ 1.6 at s = 35: −28.2 vs −26.65).
- Budget at their s = 35 point (lz_budget.py): D-part 8460 + log g 8293.28 (normalisation 38747.60 + saddle −30454.31) − C₁ 16780.37 (small primes p < n: 1028.91; single-digit range n < p < M′n: 15751.46, of which (8m₁+3m₂)·log M′ = 13108 is the numerator content ∝ 1/u and ≈ 2643 the block carries at the worst pole) = −27.09. The three big terms are each ~10⁴ and cancel to O(10); ν₀(1/u) averages 1678, 846, 430, 222, 118, 65, 37, 24 on u ∈ [1,2), [2,4), …, [128,235).
- Scale law: the objective is homogeneous of degree ≈ 1 in (m₁, m₂, δ): (209,243,δ) → (150,190,δ) at s = 33 gives +258.3 → +184.2 (ratio 1.40 ≈ 209/150). So the sign is a property of the shape; shrinking the scale only shrinks |net|. This is why the s = 33 scan's best sits at the corner (150,190): the shape is positive.
- s = 33 (34 blocks), their δ minus (64, 68): (m₁, m₂) scan m₁ ∈ [150,270], m₂ ∈ [190,310] step 10: minimum +183.3 at (150,190), rising in all directions inside the box; uniform δ-shifts +k: k = 4 gives +180.4 (resid 5), k ≥ 8 worse (the content window closes as δ_min grows). Coordinate descent (±1 on every δ_j, ±1/±2 on m₁, m₂; 4 sweeps, coarse) from three seeds: their shape → +176.76 at (150,190), δ = (4⁵,5,…,10,12,…,52,54,58); a shape with eight blocks at 4 and pairs at 5, 6, 7 → +260.9 at (150,182); a uniformly spread shape (4,6,8,9,11,…,60) → +380.9 at (150,194). The landscape in shape space is shallow (±10 for moves that keep the block count) against a gap of +177.
- Marginal cost of one block at their s = 35 point (remove block j, keep m's; lz_marginal.py): every block costs +140 to +167 per n (dD −235, dlog g +228…+312, dC₁ −70…−148); the cheapest to lose are δ = 64, 68, 42 (+140…+142), the dearest δ = 12…20 (+162…+167). After re-tuning (m₁, m₂) the pair costs ≈ +210, i.e. ≈ 105 per block. The refined rebate is 26 = Σ_j (2δ_min − δ_j)₊, a quarter of a block; the uniform shift to δ_min = 6 buys 3.5 more and nothing beyond.
- Where the s = 35 landscape sits: the coarse (m₁, m₂) scan at s = 35 (their δ) has its valley bottom near their point (−28.2 coarse at (209,243); −22.6 at (190,230); +2.8 at (240,260)) with the rebate flat at 26 wherever resid = 0.
Verdict: with the refined lemma the Lai–Zhou family gives one of ζ(7), …, ζ(35) with margin 26.6 per n instead of 0.65, and nothing below 35: s = 33 is ≈ 180 per n away, seven times the total arithmetic gain available (the refined predictor is exact at the single-digit primes, so there is no hidden denominator gain; the min-over-poles loss in Φ_n — about 4 carries per unit u, ≈ 900 per n — is the only large slack and it is structural). The two-value elimination is what makes this family expensive: its three shifted root blocks are 4·(2m₁+m₂)·log(2m₁+m₂+x₀) ≈ 17000 per n of form size, paid for by the numerator content ≈ 13100 and the small primes ≈ 1000.

# Zudilin's r = 5 tower (2026-10-01; zfast.py, zetaforms.py): one of ζ(7), …, ζ(23) — conditional on his saddle lemma for r = 5
## The literature, re-read
Zudilin already proved "one of ζ(7), ζ(9), …, ζ(35)" in 2001 (Math. Notes 70, Theorem 2: a = 32, b = 5, c = 5, d = 12 in his single-block family R_n(t) = t·[(t±(cn+1))⋯(t±(cn+dn))]^b/[t(t±1)⋯(t±cn)]^a·(2cn)!^a/(dn)!^{2b}, I_n = (1/(b−1)!)Σ_{t>cn}R^{(b−1)}(t); and ζ(9)…ζ(51) with a = 46, b = 7). So the s = 35 record for the ζ(7)-tower is 25 years old, and Lai–Zhou + Fauzan only re-derive it. In Zu04 §8, after Theorem 3 (r = 3, q = 13: ζ(5)…ζ(11)), he writes that "our previous results [Zu4] on … ζ(7), …, ζ(35) [and] ζ(9), …, ζ(51) can be improved. We are not able to demonstrate the general case of Lemma 20, although this lemma (after removing the hypothesis Re τ₀ < η₀) remains true for odd r > 3 and for any suitable choice of directions (cf. [Zu3], Section 2)." Lemma 19 (the arithmetic: D_{m₁}^r D_{m₂}⋯D_{m_{q−r}}·Φ⁻¹·F ∈ Zζ(q−2) + … + Zζ(r+2) + Z) is stated and proved for every odd r. The saddle-point lemma for a derivative kernel of general odd order b is proved in his single-block family (arXiv math/0104249, §4, Lemmas 2–6: F = −(1/2πi)∫π^b cot_b(πt)R(t)dt, sin^b·cot_b = V_b(cos) so only odd λ with |λ| ≤ b−2 contribute, saddles f′(τ) = λπi, and under a geometric condition (19) on the real root the dominant one is the root of the polynomial in Im τ > 0 with maximal real part — the same rule as Lemma 20). Nobody seems to have carried this out for the §8 family with r = 5 (Lai–Zhou 2021 cite only Zu01/Zu04 for one-value results).
## What the earlier verdict got wrong
The line "towers starting at ζ(7) through r = 5 collapse the closeness … ≈ +150 per n" was measured at q ≤ 15 (ζ(7)…ζ(13)) — the wrong regime: closeness in this family grows faster than linearly with the number of blocks (C₀ = 46, 82, 185, 239, 297, 344 at q = 13, 15, 19, 21, 23, 25 for η₀ ≈ 60) while the D-part grows linearly, so the family only becomes viable at q ≈ 25. In addition, this session's first r = 5 ladder used numpy roots for the degree-41 saddle polynomial and produced garbage (roots violating the τ ↔ η₀ − τ̄ symmetry); the fix is exact integer coefficients with mpmath polyroots (verification) and a Newton fan from the r-fold cluster near η₀ (search); the two agree on every case checked, and reproduce his r = 3, q = 13 point (C₀ 227.5802, C₂ 226.261, net −1.32; refined −5.99).
## Ladder (Zudilin's objective C₂ − C₀ with Lemma 19's D-part r·m₁ + Σm_j and κ = Φ rate; refined rebate from zetaforms.refined_rebate; coordinate descent from default shapes, numerator 0.30η₀, denominators spread 0.32–0.45η₀)
| q | slots | best η | C₀ | D-part | κ | Zudilin net | refined net |
|---|---|---|---|---|---|---|---|
| 37 | ζ(7…35) | (101; 30⁵, 32,32,33,33,34,34,35³,36,36,37³,38,38,39,39,40³,41,41,42³,43,43,44,44,45,45) | 978.42 | 1346 | 537.4 | −169.8 | −272.8 (rebate 103) |
| 37 | | search at η₀ = 100…104, still improving | | | | −232 | |
| 25 | ζ(7…23) | (62; 18⁵, 19³, 20³, 21³, 22,22,23,23,24³,25,25,26,26) | 344.34 | 583 | 249.27 | −10.61 | −15.61 |
| 25 | ζ(7…23) | (83; 24⁵, 25,25,26,26,27,27,28,28,29,29,30,31,31,32,32,33,33,34,34,35) | 460.39 | 798 | 351.8 | −14.21 | −26.88 |
| 23 | ζ(7…21) | (65; 19⁵, 20,20,21³,22,22,23,23,24,24,25,25,26,26,27,27,28) | 296.70 | 556 | 242.3 | +16.96 | +12.29 |
| 23 | ζ(7…21) | (82; 24⁵, 25,25,26,26,27,28,28,29,29,30,31,31,32,32,33,34,34,35) | 375.68 | 710 | 313.7 | +20.58 | +13.91 |
| 21 | ζ(7…19) | (65; 19⁵, 20,20,21,21,22,22,23,23,24,25,25,26,26,27,27,28) | 238.73 | 506 | 224.2 | +43.10 | +40.44 |
| 19 | ζ(7…17) | (62; 18⁵, 19,19,20,20,21,21,22,22,23,24,24,25,25,26) | 185.26 | 438 | 188.3 | +64.43 | +64.43 |
| 15, 13, 11, 9 | | | 82, 46, 6, −35 | | | +110, +146, +173, +231 | |
The nets scale with η₀ (−10.6/62 = −0.171, −14.2/83 = −0.171 per unit at q = 25; +0.26, +0.25 per unit at q = 23), so the sign is a shape property and q = 25 is the boundary: one of ζ(7), ζ(9), …, ζ(23).
## Exact checks of the q = 25 direction (62; 18⁵, …) (zetaforms.py forms)
n = 1 (h₀ = 64, 25 poles, order ≤ 20): partial fractions reproduce R ✓, A_m = 0 for even m and m = 5 ✓, slots exactly ζ(7), …, ζ(23); log|F|/n = −384.24; true primitive height 316.60 per n (Lemma 19 bound 410.12, refined 312.58); NET log|F·D|/n = −67.64. n = 2 (h₀ = 126): log|F|/n = −369.11 → toward −344.34; true height 319.55 (bound 382.76, refined 319.06); NET −49.56; refined predictor exact at 10/10 single-digit primes; the only violations are p = 2, 3 (below √h₀, outside the lemma). n = 3, 4 running. Saddle: τ₀ = 61.9837 + 2.4790i, Re τ₀ < η₀ (Lemma 20's hypothesis holds as stated), Im f₀(τ₀)/π = −182.1155 ∉ Z; the other root near η₀ (60.31 + 0.91i) has Re f₀ = −355.38 (the λ = 1 saddle; the dominant one is λ = 3 as in his Lemma 6), and the roots on Re τ = η₀/2 have Im f₀ ∈ πZ exactly and are the even-λ saddles the odd kernel does not see. Same picture for the q = 37 direction: n = 1, 2, 3 give log|F|/n = −1019.06, −1008.15 (→ −978.42 predicted), true height 703.2, 692.2 per n, NET −315.9 per n at n = 2 and 3; refined predictor exact at 14/15 and 23/23 single-digit primes.
## q = 23 under the refined lemma (later on 2026-10-01; zfast `fast_rebate`, a u-grid version of zetaforms.refined_rebate, agrees with the Farey version to 0.03)
Zudilin's own bound never goes negative at q = 23 (+5.0 is the best seen, at (102; 31⁴, 32, 34,34,35,35,…,43)), but the refined rebate is strongly shape-dependent (4.7 … 60.8 per n) and the refined objective does: (100; 30⁴,31, 32,32,33,34,34,35,36,36,37,38,39,39,40,41,41,42,43,44): C₀ 452.89, D-part 805, κ 340.55, Zudilin net +11.55, rebate 35.36 → refined −23.8; refined-objective descent → (98; 29,30,30,31,31, 32,33,33,34,34,35,36,36,37,38,38,39,40,40,41,42,43,44): Zudilin +29.67, rebate 60.82, refined −31.15; at η₀ = 128: (128; 36,36,37,38,39, 40,41,42,42,43,…,56) refined −44.4 and still descending. Exact checks: (100; 30⁴,31,…) n = 1, 2, 3: identity ✓, vanishing ✓, slots ζ(7)…ζ(21); log|F|/n = −468.72, −451.18, −444.15 (predicted −452.89; here the finite values sit above the asymptote — the prefactor's sign differs from the q = 25 case); true height 396.4, 431.7, 426.2 per n against Zudilin's bound 549.8, 568.8, 543.0 and the refined 395.7, 431.7, 423.4; refined predictor exact at 7/8, 15/15, 20/23 single-digit primes (misses only at p ≤ 11 < √h₀). (98; 29,30,30,31,31,…) n = 2, 3: single-digit heights 229.51 (= refined, 8/8) and 286.54 (refined 288.23, 14/15); log|F|/n = −423.0 at n = 3 (predicted −426.63). So, with the refined lemma in place of Lemma 19's D-part: **one of ζ(7), …, ζ(21)** (margin ≈ 31 per n at η₀ = 98). q = 21 stays positive: +20.5 refined at η₀ = 107 (descent still running), +19.4 at 81; the Zudilin-bound values are +42 to +65.
## Validation of the r = 5 saddle rule at larger n
(30; 7⁵, 9,9,10,10,11,11,12,12), q = 13 (slots ζ(7…11)); the exact degree condition Σh_j ≤ h₀(q−r)/2 fails for n ≤ 4 here (asymptotic condition holds), so n = 5, 6, 7, 8, 10, 12: log|F_n|/n = −30.18, −29.08, −28.54, −28.34, −27.47, −27.08; least squares on log|F_n| = −C₀n − A log n + B gives C₀ = 24.77, A = 1.36 (residuals ±1.4, the cosine of the conjugate pair), against the predicted 24.1268 from the max-Re root 28.932 + 2.937i; the other near-η₀ root (28.030 + 0.822i, Re f₀ = −35.36) is excluded. Same rule, same outcome, as in the q = 9 toy.
## CORRECTION (2026-10-01, later): the fast saddle picked the wrong root on the q = 23 shapes
The search evaluator's Newton fan (zfast.newton_roots) found only one of the two cluster roots near η₀ with Im > 0 — the one with the smaller real part (Im ≈ 2–3) — and missed the max-Re root (Im ≈ 7) on every q = 23 shape and on (124; 37⁴, 38, 40, …) at q = 25. The max-Re root has the larger Re f₀ (it is the dominant, λ = r − 2 saddle), so the fan overstated C₀ by 20–31 per n and every q = 23 "negative" above is an artifact; the polyroots listings (zf_roots.py) had both roots all along, only the selection was wrong. Fixed: zfast.cluster_roots starts Newton on circles of the cluster radius ρ = η₀·(∏η_j/∏(η₀−η_j))^{1/r} around η₀, requires exactly (r−1)/2 roots with Im > 0, Re > η₀/2, and falls back to polyroots otherwise; re-verified against polyroots on all twelve reported shapes (tmp/zf_recheck.py, all agree). Corrected values (Zudilin net | refined net per n):
| shape | before (wrong) | corrected |
|---|---|---|
| q=25 (62; 18⁵, …) | −10.61 \| −15.61 | unchanged |
| q=25 (83; 24⁵, …) | −14.21 \| −26.87 | unchanged |
| q=25 (124; 37⁴, 38, 40, 40, 41, 41, 42, 42, 43, 44, …, 55) | −32.33 \| −77.99 | −2.96 \| −48.66 |
| q=23 (100; 30⁴, 31, 32, 32, 33, 34, 34, 35, …, 44) | +11.55 \| −23.78 | +39.33 \| +3.97 |
| q=23 (98; 29, 30, 30, 31, 31, 32, …, 44) | +29.67 \| −31.15 | +59.84 \| −0.98 |
| q=23 (102; 31⁴, 32, 34, 34, 35, 35, …, 43) | +5.05 \| −23.95 | +33.13 \| +4.13 |
| q=23 (128; 37⁴, 38, 40, 40, 41, 42, 43, 43, 44, …, 55) | −20.08 \| −35.08 | +10.28 \| −4.70 |
| q=23 (130; 35, 36, 37, …, 44, 44, 45, …, 56) | −5.61 \| −47.97 | +25.09 \| −17.27 |
| q=21 (107; 29, 29, 30, 31, 31, 32, …, 47) | +40.84 \| +20.52 | +70.52 \| +50.20 |
| q=31 (107; …), q=37 (104; …) | −152.4 \| −187.1 ; −238.9 \| −307.9 | unchanged |
The exact forms decide between the two roots and confirm the max-Re rule: for (100; 30⁴, 31, …) the exact log|F_n|/n = −468.7, −451.2, −444.2 (n = 1, 2, 3) has already crossed the sub-dominant root's level −452.9 and is heading for the max-Re level −425.1; for (98; 29, …) the values −439.7, −423.0, −415.5 crossed −426.6 at n = 2 (max-Re level −396.5); for the q = 13 shape the fit gave 24.77 against the max-Re 24.13 and the other root's 35.36. So the forms are as large as Zudilin's rule says, not as small as the fan claimed. Consequences: "one of ζ(7), …, ζ(23)" stands (q = 25, margins 10.6 / 14.2 / 3.0 per n under his bound at η₀ = 62 / 83 / 124; refined 15.6 / 26.9 / 48.7). "one of ζ(7), …, ζ(21)": Zudilin's bound is ≥ +10 on every q = 23 shape seen; the refined net is −17.3 at (130; 35, …, 56) and −4.7 at (128; …) — negative, but these shapes were found with the wrong objective, so the corrected searches decide the true q = 23 margin (next paragraph); nothing about q = 21.
## Corrected searches (2026-10-01, last; cluster saddle verified against polyroots at every reported point)
- q = 25, Zudilin's objective from (124; 37⁴, 38, 40, …): → **(131; 37⁴, 38, 40, 40, 41, 41, 42, 43, 44, 44, 45, 46, 47, 48, 49, 49, 50, 51, 52, 53, 54, 55)**: C₀ 737.761, D-part 1247, κ 549.594, C₂ 697.406, **Zudilin net −40.35**, exact Farey rebate 14.667, refined −55.02. Caveat: τ₀ = 131.028 + 4.952i has Re τ₀ > η₀ by 0.03 (his Lemma 20 hypothesis as stated fails; he says it can be removed); the shapes (62; 18⁵, …) (−10.6) and (83; 24⁵, …) (−14.2) satisfy Re τ₀ < η₀. Exact forms at n = 1, 2: identity ✓, vanishing ✓, slots ζ(7)…ζ(23); log|F|/n = −787.96, −767.74 (→ −737.76); true heights 677.06, 674.49 per n (his bound 848.2, 807.8); single-digit range: refined 450.82 / 538.04 vs true 448.26 / 536.06, exact at 10/11 and 21/22 primes, violations only at p ≤ 11 < √h₀.
- q = 25 (62; 18⁵, …) exact to n = 6: log|F_n|/n = −384.24, −369.11, −363.02, −359.30, −356.70, **−354.98** — it has crossed the sub-dominant root's level −355.38 and heads for the max-Re level −344.34. The discriminating check for the direction itself: the forms are governed by Zudilin's root. Height at n = 6: true 321.16 per n (his bound 376.24, refined 321.65); single-digit refined exact at 25/27 primes, the two misses at p = 13 < √h₀ (h₀ = 374).
- q = 23, refined objective: (134; 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56): C₀ 642.58, D-part 1223, κ 566.70, Zudilin +13.72, rebate 38.66, **refined −24.94**; (148; 41, 41, 41, 42, 43, 45, 45, 46, …, 62): Zudilin +8.37, refined −9.97; (102; 29, 30, 30, 31, 31, …, 43): +22.60 / −7.40; (104; 30, 30, 30, 31, 31, …): +12.51 / −0.81. His bound never goes below +8 at q = 23; the refined margin grows with the scale (−17 at 130, −25 at 134, −29.46 at (170; 46,47,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,65,66,67,68,70,71) after 5 sweeps, when the system killed the run for low memory — the memory was held by stale exact-form jobs from the earlier session, which were killed at the same time; restarted from that shape it converged to **(170; 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 66, 67, 68, 70, 71): C₀ 816.12, D-part 1531, κ 703.60, Zudilin +11.29, rebate 44.00, refined −32.71** — about −0.19 per unit η₀, the same scale law as at q = 25; his own bound stays at +11). The (134; …) shape confirmed with polyroots (C₀ 642.5757, exact Farey rebate 38.667 vs fast 38.66) and exact forms: n = 1 (h₀ = 136, 55 poles, order ≤ 18): identity ✓, vanishing ✓, slots ζ(7)…ζ(21), log|F|/n = −689.45, true height 598.59 per n (his bound 805.23; refined 598.02), single-digit refined exact at 12/12 primes; n = 2 (h₀ = 270, 109 poles): log|F|/n = −670.51 (→ −642.58), true 611.18 (bound 762.03), single-digit refined 486.13 = true, exact at 24/24 primes; violations only at p ∈ {2, 5, 11} < √h₀.
- q = 21: (104; 30⁴, 31, 32, 32, 33, 34, 35, 36, 36, 37, …, 44): Zudilin +53.97, refined +48.31. Positive at every seed and scale tried.
Summary of what stands after the correction: **one of ζ(7), …, ζ(23)** on Lemma 19 + the r = 5 saddle lemma (margin 40 per n at η₀ = 131, 14 at η₀ = 83 with Re τ₀ < η₀); **one of ζ(7), …, ζ(21)** only if the refined lemma is added (margin ≈ 25 per n); ζ(7)…ζ(19) out of reach for this family (+48 refined).
## Status
Arithmetic: Zudilin's Lemma 19 (all odd r) plus the refined lemma (finite predictor exact at every single-digit prime above √h₀ in every run: q = 25 n = 1…4, q = 37 n = 2, 3, q = 23 n = 1…3; proof written for the β-forms in refined_lemma_proof.md, the odd-zeta transfer is verified but not yet written as a proof). Analysis: lim sup log|F_n|/n = −C₀ with C₀ from the max-Re root is Zudilin's Lemma 20 read with r = 5 — asserted by him for odd r > 3, proved by him for r = 3 in this family and for every odd derivative order b in the single-block family (math/0104249 §4: kernel π^b cot_b, sin^b·cot_b = V_b(cos), contributing saddles λ odd with |λ| ≤ b−2, dominant λ = b−2 = the max-Re root under a geometric condition on the saddle curve); not written for r = 5 in this family; our exact values approach the predicted constants (q = 13 fit 24.77 vs 24.13; q = 25: −384, −369, −363, −359 → −344; q = 37: −1019, −1008 → −978; q = 9 toy → +2.29). Everything else is standard (A_m ∈ Q, non-vanishing from the lim sup). Statements, by what they rest on (numbers as corrected above and in "Corrected searches"):
- **one of ζ(7), ζ(9), …, ζ(23) is irrational** — Lemma 19 as printed + the r = 5 saddle lemma (margin 40.4 per n at (131; 37⁴, 38, 40, …, 55); 14.2 at (83; 24⁵, …) where Re τ₀ < η₀ holds as he states it);
- **one of ζ(7), …, ζ(21)** — only with the refined lemma (margin 24.9 per n at (134; 35, 36, …, 56)); his own bound stays ≥ +8 at q = 23.
The first improves Zudilin 2001's ζ(7)…ζ(35) and Lai–Zhou + Fauzan. Before any claim: (i) write the r = 5 asymptotic lemma (or verify the geometric condition of math/0104249 Lemma 4 for the chosen η and cite); (ii) write the refined lemma for the integer-pole family; (iii) confirm the final η with `zfast.py net 5 ETA full refexact` and exact forms at n = 4; (iv) literature check that nothing below 35 has been published for the ζ(7)-tower since 2001.

# Siegel disk (2026-09-29, later): G as a rotation number (siegel.py)
The Siegel disk of P_θ(z) = e^{2πiθ}z + z² exists iff θ is a Brjuno number, B(θ) = Σ β_{n−1}log(1/θ_n) < ∞ (θ_n the Gauss-map iterates, β_n = θ₀…θ_n, so the n-th term is ≈ log(a_{n+1})/q_n: the carry at depth n weighted by its scale); Yoccoz: log r(θ) + B(θ) = O(1) with r the conformal radius — the disk is exactly as large as the arithmetic of the carries allows, the closeness = height law of this world, with no margin. B satisfies the self-similar equation B(x) = −log x + xB(1/x), the same shape as the ladder's functional equations.
## Portraits (2000 partial quotients each; siegel.py cf)
| θ (mod 1) | first quotients | Khinchin (2.6854) | Lévy (3.2758) | largest a_n (position) | B(θ) |
|---|---|---|---|---|---|
| G | 1, 10, 1, 8, 1, 88, 4, 1, 1, 7, 22 | 2.685 | 3.261 | 6328 (984) | 2.4883 |
| 2G/π | 1, 1, 2, 1, 1, 32, 4, 20 | 2.578 | 3.172 | 24749 (594) | 1.5234 |
| G/π | 3, 2, 3, 16, 8, 10, 3 | 2.857 | 3.446 | 12374 (586) | 1.7388 |
| ln(2+√3) | 3, 6, 2, 4, 1, 2, 3, 3, 3 | 2.632 | 3.228 | 1716 (77) | 1.8202 |
| R/3 (Ramanujan number) | 15, 1, 3, 7, 1, 1, 3, 2, 4 | 2.695 | 3.285 | 1192 (1551) | 2.8639 |
| e^{2G/π} | 1, 3, 1, 3, 1, 38, 1, 3, 1, 9 | 2.659 | 3.262 | 6669 (907) | 1.7191 |
| L(2,χ₋₃) | 1, 3, 1, 1, 2, 1, 17, 1, 10 | 2.570 | 3.164 | 4081 (1353) | 1.5936 |
| ζ(3) | 4, 1, 18, 1, 1, 1, 4, 1, 9, 9 | 2.713 | 3.302 | 26550 (1852) | 2.1862 |
| e (control, patterned) | 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8 | 3.563 | 4.429 | 7982 (885) | 1.4420 |
| π (control) | 7, 15, 1, 292, 1, 1, 1, 2 | 2.624 | 3.213 | 20776 (431) | 2.3975 |
Gauss–Kuzmin frequencies of 1, 2, 3, 4 match the law (0.415, 0.170, 0.093, 0.059) within noise for every constant except e. G is a typical irrational to 2000 quotients (the user's certified 11,085 terms say the same); the extreme values (6328 in 2000 digits) are ordinary heavy-tail statistics. B(G) = 2.488 is large mostly because G is close to 11/12: the third convergent, |G − 11/12| = 7.0·10⁻⁴, so θ₁ = 0.0917 and the term β₀log(1/θ₁) = 2.19 carries most of B; in the ladder's language 1 − G ≈ 1/12 (0.0840 vs 0.0833) is the nearest small-denominator neighbour of Catalan's constant.
## The disks (siegel.py disk; linearisation coefficients to n = 4000, r = 1/limsup|h_n|^{1/n} over n ∈ [2000, 4000], upper estimates where a large quotient lies beyond 4000)
| θ | r | log r | B | log r + B (Yoccoz, bounded) |
|---|---|---|---|---|
| golden mean | 0.3263 | −1.120 | 1.2598 | 0.14 |
| 2G/π | 0.2945 | −1.223 | 1.5234 | 0.30 |
| ln(2+√3) | 0.3106 | −1.169 | 1.8202 | 0.65 |
| G | 0.2469 | −1.399 | 2.4883 | 1.09 (r is an upper estimate: the small divisor at n − 1 = q₆ = 10579 is beyond the range) |
| R/3 | 0.2271 | −1.482 | 2.8639 | 1.38 |
## Dedekind carry (siegel.py dedekind): the second cocycle of the same digits
Exact Dedekind sums of G's convergents against the alternating sum of the partial quotients (Barkan–Hickerson): 12 s(p_n, q_n) − Σ_{i≤n}(−1)^{i+1}a_i stays in (−2.1, 0.9) for n ≤ 200 (q₂₀₀ has 106 digits). The Dedekind sums of the convergents swing with the carries: −102 at n = 6 (the 88), +466 at n = 125, +87 at n = 200. So the two natural cocycles of the Euclidean carries are the Brjuno sum (log-weighted by scale: the Siegel disk) and the Dedekind sum (alternating-signed: the modular side, with its 1/12 normalisation). Nothing here bears on the irrationality of G — a Siegel-disk argument presupposes an infinite continued fraction — but the portrait is recorded: G behaves like a generic Brjuno number whose disk is shrunk by its proximity to 11/12.
√3-test (user's question, tmp/root3_test.py): the 31 quantities of the portrait (B, r, log r, Yoccoz differences, Khinchin, Lévy, 1−G, 2G/π) against ~600 forms q·√3·9^k (q = a/b ≤ 12, k = −2…1) and 2√3−1, (2√3−1)/2, √3−1, 2−√3: best matches have relative error 2·10⁻³–4·10⁻² except two chance hits at the expected density (B(e) ≈ (5/6)√3 to 9·10⁻⁴, U(golden) ≈ (8/11)√3/9 to 5·10⁻⁴). The only exact closed form is B(golden) = φ²ln φ = 1.259829 (all θ_n equal; ≈ ∛2 by coincidence, 3·10⁻⁶ apart); B(√3−1) = 1.431123, B(√2−1) = 1.504599 exactly by the same periodic-CF formula. Where √3·9⁻¹ is exact in this program: L(1,χ₋₃) = π√3/9, hence π/√3 = 3L(1,χ₋₃) (the level-1 constant of the Ramanujan tower), ln(2+√3) = √3·L(1,χ₁₂), and π ln(2+√3) = 9·L(1,χ₋₃)·L(1,χ₁₂); also 2G/π = 0.5831 sits 1% from 1/√3 = 0.5774 — a near miss, not an identity (β′(−1) = 2G/π is transcendental-looking and 1/√3 is algebraic).

# Rerun of the r = 5 tower from the machine (2026-09-28, real date; the day labels above ran ahead of the calendar and were all written by 2026-09-26)
New session. zfast.py `net r ETA` with full polyroots and the exact Farey rebate (scratchpad driver zcheck.py, six shapes in parallel, 22 to 262 s each) reproduces every recorded value to the printed digits:
| shape | slots | C₀ | D-part | κ | Zudilin net | rebate | refined net | τ₀ |
|---|---|---|---|---|---|---|---|---|
| r = 3 control (91; 27³, 29, …, 38) | ζ(5…11) | 227.5802 | 403 | 176.7390 | −1.319 (his −1.331) | 4.667 | −5.986 | 87.479 + 3.328i |
| (62; 18⁵, …) | ζ(7…23) | 344.3436 | 583 | 249.2694 | −10.613 | 5.000 | −15.613 | 61.984 + 2.479i |
| (83; 24⁵, …) | ζ(7…23) | 460.3881 | 798 | 351.8246 | −14.213 | 12.667 | −26.879 | 82.977 + 3.319i |
| (131; 37⁴, 38, 40, …, 55) | ζ(7…23) | 737.7608 | 1247 | 549.5940 | −40.355 | 14.667 | −55.021 | 131.028 + 4.952i |
| (134; 35, …, 45, 45, …, 56) | ζ(7…21) | 642.5757 | 1223 | 566.7012 | +13.723 | 38.667 | −24.944 | 133.723 + 6.168i |
| (170; 46, …, 60, 62, 63, 64, 66, 67, 68, 70, 71) | ζ(7…21) | 816.1156 | 1531 | 703.5974 | +11.287 | 44.000 | −32.713 | 169.658 + 7.815i |
Not recorded before: both q = 23 shapes have Re τ₀ < η₀ (by 0.28 and 0.34), so the ζ(7)…ζ(21) statement meets Lemma 20's hypothesis as Zudilin printed it; only (131; …) needs the hypothesis removed (Re τ₀ − η₀ = +0.028). The standing statements are unchanged: one of ζ(7), …, ζ(23) on Lemma 19 + the r = 5 saddle lemma (margin 14.2 per n at η₀ = 83, 40.4 at 131); one of ζ(7), …, ζ(21) with the refined lemma added (margin 24.9 at 134, 32.7 at 170). The published record for the ζ(7)-tower is Zudilin 2001, ζ(7), …, ζ(35); Lai–Zhou 2021 + Fauzan re-derive 35, and their construction lowered to s = 33 is closed (≥ +177 per n). The "33" of the final report was the ζ(21) margin, 32.7 per n, not a ζ index.

# Literature check, and a parallel proof of ζ(7)…ζ(21) (2026-09-28, real date)
## Literature (subagent sweep; key items read in the source)
- Refereed/arXiv record for the ζ(7)-tower is still 35: Zudilin, Math. Notes 70 (2001) 426–431, Thm 2 (ζ(7)…ζ(35); ζ(9)…ζ(51) is Thm 3); Lai–Zhou, Publ. Math. Debrecen 101 (2022) 353–372 (two of ζ(5)…ζ(35)); Rivoal–Zudilin, Sém. Lothar. Combin. 81 (2020) B81b (two of ζ(5)…ζ(69)). Later work is asymptotic or p-adic: Lai 2501.05321, Lai 2407.14236, Fischler 2109.10136, Lai–Sprang 2306.10393, Lai–Lupu–Sprang 2505.23088, Lai–Sprang–Zudilin 2505.05005. Lemma 20 for r > 3: not in the refereed/arXiv literature.
- Denominator refinement: not found in refereed form. Nearest items: Zudilin JTNB 2004 §9 conjecture, one fewer D factor, still open per Rivoal 2006. Zudilin JTNB 15 (2003) §7: he proves D_n⁸Φ̃_n^{−3} for a second-derivative form in ζ(7), ζ(5), observes D_n⁷Φ̃_n^{−2} up to n = 1000, and asks "What is a trick that makes arithmetic as it is?" Krattenthaler–Rivoal Mem. AMS 2007 cover equal parameters only.
- News: Fauzan, "ζ(5) is irrational", Zenodo 10.5281/zenodo.22826419 (17 Sep 2026). Lean: mo271/zeta5 (no sorry) and danromik/zeta5-irrationality (PNT as an axiom). Calegari, Persiflage blog, 24 Sep 2026: verified in Lean; the method gives ζ(k) for k ≤ 5, not directly beyond. Anand, "Zeta 7 is Irrational", Zenodo 10.5281/zenodo.22920911 (23 Sep 2026): the v1 record admits an "optimality mismatch error" in Thm 1.3.
## gmDevi/zeta-7-21-lean
Created 2026-09-28 10:15 UTC. AI-produced by Claude models under a maintainer's direction, unrefereed; the work is dated 2026-09-24…26.
- Lean 4 + Mathlib theorem ZetaWindow.zeta_7_to_21_not_all_rational: at least one of ζ(7), …, ζ(21) is irrational. Axioms: propext, Classical.choice, Quot.sound.
- Same construction: Zudilin §8, r = 5, q = 23, η = (160; 47⁴, 48, 50, 50, 51, 52, …, 66), Lemma 19 unrefined.
- Lemma 20 is replaced by two results:
  - Theorem U, a contour upper bound |F̃_n| ≤ K n¹⁶ e^{−C₀′n} with C₀′ ≥ 750.7117 (sharp line); the Lean proof uses 748.1.
  - Theorem N: for n even with ℓ = 63n − 1 prime, v_ℓ(c_s) ≥ 0 for every ζ coefficient and v_ℓ(c₀) = −5 exactly. This is the Lai–Sprang real auxiliary-prime criterion; Dirichlet gives infinitely many such ℓ.
- The −5 comes from the two far poles k = 110n, 110n + 1. There k − h₁ = ℓ, ℓ + 1, and the two η = 50 blocks give pole order 2. This is the partial-sum window at its edge: only the constant term sees ℓ.
- Their docs/proof.md also refines Lemma 19 (Theorem H′1: far poles k − h₁ ≥ p; near and flat poles by residue mod p). It gives C₂′ ≤ 729.0513, a gain of 18, and a certified margin C₀′ − C₂′ ≥ 21.66 per n. It is not used in the Lean proof.
## Our machine on their η (zfast.py net 5 … full refexact, 202 s)
C₀ 750.7117 (their sharp C₀′), D-part 1341, κ 593.9388 (u-grid; C₂ 747.0612 against their 747.0513), Zudilin net −3.65 per n, exact Farey rebate 18.0000 (their H′1 gain), refined net −21.65 (their certified 21.66), τ₀ = 159.356 + 8.394i (Re τ₀ < η₀).
1. CORRECTION: "Zudilin's bound ≥ +8 at every q = 23 shape" was a search failure. Their shape is −3.65 under his plain arithmetic, 12 per n better than our best, (148; …) at +8.37. Its numerators are 47⁴, 48 ≈ 0.294η₀; its denominators 50, 50, 51, …, 66 = 0.3125…0.4125η₀, with the first one doubled. So ζ(7)…ζ(21) never needed the refined lemma, and our q = 21 verdict (+48 refined) came from the same search and is not reliable.
2. Two independent refinements of Lemma 19 give the same rebate at this η, 18.000 per n: ours (window, saturation, exact top valuation) and their H′1 (far/near/flat).
3. Their Theorem N turns the partial-sum window into a nonvanishing certificate. One prime just below m₁n sees only the constant term (m₁ = 63, so ℓ ≡ −1 mod 63, the top base-63 digit). Together with an upper bound this removes the saddle lemma, and the same device should work for any shape in this family that has such a top-range prime.
Open, ranked: (a) ζ(7)…ζ(19) (q = 21), needing a better search seeded with the doubled-first-denominator pattern, plus the refined lemma and the one-prime certificate; (b) Zudilin's 2003 question, whether the refined lemma predicts his observed D_n⁷Φ̃_n^{−2}; (c) the ζ(9)-tower with r = 7, where nothing is in the literature.

# Zudilin's 2003 question, "What is a trick that makes arithmetic as it is?" (2026-09-28; zu2003.py)
## The question (JTNB 15 (2003) 593–626, §7; www.math.ru.nl/~zudilin/PS/z234.pdf)
Zudilin looks at three forms.
- F̃_n = ½ Σ_{t≥1} (d/dt)²[(2t+n)((t−1)⋯(t−n)(t+n+1)⋯(t+2n))³/(t(t+1)⋯(t+n))⁶] = ũζ(7) + w̃ζ(5) − ṽ. He proves the denominator divides D_n⁸Φ̃_n^{−3}; his recursion up to n = 1000 shows D_n⁷Φ̃_n^{−2}. The question is asked right after this form.
- The Ball-type forms F_{k,n}, k = 5, 7: proved D_n^{k+1}Φ̃_n^{−1}, observed D_n^k.
- Φ̃_n = ∏_{p<n, {n/p}∈[2/3,1)} p, with density ψ(1) − ψ(2/3) − ½ = 0.2410.
## The machine
Under t = τ + n + 1 these are zetaforms' Zu04 §8 forms with h₀ = 3n + 2 and every h_j = n + 1:
- F̃_n is r = 3, η = (3; 1⁹). It sits on the boundary of Zu04's exact degree condition (Σh_j = 9n+9 against 9n+6), but deg R = −5, which is all the construction needs. zu2003.LooseForm drops that assertion.
- F_{k,n} is r = 1, η = (3; 1^{k+2}).

Controls:
- F_{5,n} reproduces his printed u, w, v at n = 1, 2 exactly (18, 66, 98; 938, 6125/2, 74463/16).
- F̃_n equals a direct numeric sum of his expression to 28 digits at n = 1, 2, 3.
- Slots: ζ(5), ζ(7) for F̃; ζ(3), ζ(5) (k = 5) and ζ(3), ζ(5), ζ(7) (k = 7) for F.
## Tables
Command: zu2003.py table, n = 2…40, e_p = exponent of p in the common denominator. Output files: zu2003_B_2_40.json, zu2003_A5_2_40.json, zu2003_A7_2_40.json.
- Zudilin's observed exponent (7v_p(D_n) − 2[p∈Φ̃], resp. k·v_p(D_n)) holds at every prime for every n, small primes included.
- At the large primes (p² > 3n + 2) outside Φ̃, the refined predictor, which is exact term by term, says 8 (resp. k + 1) and the truth is 7 (resp. k).
- At Φ̃ primes: truth = observed = proved = refined = 5 (form F̃).
- So the saving is not a valuation effect on single terms. It is a cancellation between terms: Mechanism II.
## The mechanism, per prime p with p² > 3n + 2 (n = ap + b, 0 ≤ b < p, a ≥ 1; window pole j = k − h₁ = cp + j′, 1 ≤ c ≤ a)
1. Top coefficient: B_{6,k} = (−1)ⁿ(n − 2j)·C(n+j,n)³·C(2n−j,n)³·C(n,j)⁶, an integer.
2. Carries (Kummer) cut out the leading arc. The pole is a p-unit at top order iff:
   - no borrow in n − j, so j′ ≤ b (C(n,j));
   - no carry in n + j, so b + j′ < p (C(n+j,n));
   - no carry in n + (n − j), so 2b − j′ < p (C(2n−j,n));
   - 2j′ ≠ b (centre factor).
   The leading set in each cell is L_c = {j′ ∈ (2b − p, p − b) ∩ [0, b], 2j′ ≠ b}. It is symmetric under j′ ↦ b − j′, and empty iff b ≥ 2p/3, which is exactly Φ̃. The two cutting conditions are mirror images of each other, and the three binomials put the threshold at a third.
3. The palindrome mod p. The local mirror j′ ↦ b − j′ in cell c is k ↦ h₀ + (2c − a)p − k, i.e. the global well-poised symmetry k ↦ h₀ − k shifted by a multiple of p. By Lucas, C(n,j) ≡ C(a,c)C(b,j′) is unchanged. The pair C(n+j,n), C(2n−j,n) ≡ C(a+c,a)C(b+j′,b), C(2a−c,a)C(2b−j′,b) is swapped. n − 2j ≡ b − 2j′ changes sign. So B_{6,k*} ≡ −B_{6,k} (mod p), and at the fixed point 2j′ = b the centre factor carries the p.
4. The Taylor profile is constant on the cell. The p-divisible distances from a leading pole are:
   - lower numerator block: −tp for t = c+1, …, a+c (multiplicity 3);
   - upper numerator block: +tp for t = a−c+1, …, 2a−c (multiplicity 3);
   - poles: tp for t ∈ [−c, a−c], t ≠ 0 (multiplicity 6).
   These depend only on (a, c). So B_{i,k}·p^{6−i} ≡ B_{6,k}·[u^{6−i}]G_{a,c}(u) (mod p) with one universal G_{a,c}. The partial-sum weight is p^m·PS_m(j) ≡ H_m(c).
5. Hence, per order i and per cell c, the level-p^{−8} coefficient of A₀ is a constant times Σ_{j′∈L_c} B_{6,k}, which is ≡ 0: the pairs cancel and the fixed point is not leading. The ζ-coefficients have exponent ≤ 3. So v_p(den F̃_n) ≤ 7, and ≤ 5 when L_c = ∅. Every constant above is p-integral because p > 3a (from p² > 3n + 2).

For p > √(3n+2) this is Zudilin's observed exponent, derived. The same argument with r = 1 gives k for F_{k,n}. anatomy output (per order, per pole, with partners): F̃ at n = 40, p = 29 (one cell, 12 poles, pairs k ↔ 151 − k); p = 17 (two cells); p = 13 (three cells); n = 37, p = 31 (fixed point excluded); F_{5}, F_{7} at n = 40, p = 29, 17. All per-order sums ≡ 0. Control at a Φ̃ prime (F_7, n = 38, p = 13): the leading sums do not vanish and the exponent stays at the proved 7.
## Reading
- Zudilin's "trick" is Mechanism II, the same fact as Theorem 4 of refined_lemma_proof.md for the β-forms, now in the ζ-forms. The well-poised palindrome, read modulo p inside each base-p digit cell of the partial-sum window, pairs the uncarried poles with opposite signs at every derivative order.
- The half-fold (the mirror) does the cancelling. Three carry conditions from three binomials decide where it acts, and they put the threshold at 2/3, which is Φ̃'s ψ(2/3).
- Refined lemma, gmDevi's H′1 and Zudilin's Φ are all term-by-term valuation bounds and cannot see this; it lives in the sum.
## Open
- Small primes p ≤ √(3n+2). The data satisfy the observed bound there too, often strictly, but multi-digit cells are not yet derived. They contribute O(√n) to log D, so they do not affect asymptotics.
- Whether Krattenthaler–Rivoal (Mem. AMS 2007) already cover F_{k,n} as an equal-parameter series (check). F̃_n, the derivative form the question is attached to, is not covered there.
- The general Zu04 family with distinct η_j: Mechanism II fires when the argmin arc is mirror-symmetric (per-order rule of the β-work). How much it is worth for the r = 5 tower (e.g. q = 21) is not yet measured.
## Checks (2026-09-28, later)
- Single-digit lemma, checked exhaustively (zu2003_lemma_check.py 600): 30,533 cases (n ≤ 600, primes p ≤ n with p² > 3n + 2), 1,750,916 mirror pairs, 0 failures, for all four ingredients:
  - the carry-free window poles are exactly {cp + j′ : 2b − p < j′ < p − b, j′ ≤ b, 2j′ ≠ b} (7,785 cases have an empty leading set);
  - g_{j*} ≡ −g_j mod p;
  - the carry-free fixed point has p | n − 2j (27,987 fixed points);
  - the p-divisible distances seen from a leading pole depend only on (a, c).
- Exact tables extended to n = 80 (zu2003_B_41_80.json): Zudilin's observed exponent holds at every prime for every n ≤ 80.
- Small primes (p² ≤ 3n + 2), scratchpad digit_box.py, n ≤ 80, 358 cases: the observed bound 7·v_p(D_n) − 2[p ∈ Φ̃] holds in all 358.
  - Clean examples where the leading poles are exactly the product of per-digit carry-free arcs, and each lower digit's flip d_i ↦ n_i − d_i cancels one level:
    - n = 40 = 1111₃: 8 leading poles 1d₂d₁d₀ with d_i ∈ {0, 1}; 3⁻²³ → 3⁻²⁰, depth 3 = L;
    - n = 36 = 121₅: 6 poles; 5⁻¹⁶ → 5⁻¹⁴, depth 2 = L;
    - n = 30 = 110₅: the last-digit fixed point supplies the centre-factor p.
  - This "digit box" holds in only 76 of the 358 cases (the leading set lies inside the box in 151). In the rest, carries at some digit already raise the leading level, and the depth is below L in 211 cases, with the bound still met.
  - So at small primes the bound comes from carries and digit-wise palindromes together. Not derived. Asymptotically irrelevant (O(√n)).

# The stand-out primes: a Lucas law for the second digit, and a priority correction (2026-09-28, later; zu2003_digit2.py)
## CORRECTION (literature)
Krattenthaler–Rivoal, Mem. AMS 186 (2007) no. 875, arXiv math/0311114, already prove both of Zudilin's §7 observations for ALL primes:
- Théorème 1 covers F_{k,n} = 2·S_{n,k+1,1,0,1}(1): D_n^k F_{k,n} ∈ Z-span.
- Théorème 5 (r = 1, A = 6, B = 3, C = 2) covers F̃_n: D_n⁷Φ̃_n^{−2}. Their text says this case "répond à la question sur la série F̃n à la fin du paragraphe 7 dans [55] [= Zudilin 2003], à un facteur 2 près".

Their proof is global:
- Andrews' identity turns the series into multiple binomial sums.
- Proposition 7, "divisible par k", removes the top power.
- Lemme 8, "Φ_n divise C(n+j,n)C(2n−j,n)", via carries for {n/p} ≥ 2/3; Lemme 12 gives Φ̃.

Hence the lines above that say "F̃_n … is not covered there" and "small primes not derived" are superseded: the statement is theirs, for every prime. What is ours is the per-prime reading for p² > 3n + 2: the local mirror of the palindrome inside each base-p digit cell. Our carry arcs are their Lemme 8 in cell form. The local-mirror formulation was not found in the literature (subagent sweep: citers of Zudilin 2003 and the Memoir; Fischler's Bourbaki Remark 2.13 only says "pour z = 1 on a des compensations particulières"). Beyond the symmetric family the denominators conjecture is open: Zudilin's asymmetric conjecture, KR §17.1; Lai–Zhou Remark 4.4 "tremendously difficult"; Marcovecchio–Zudilin 1905.12579 has only the most symmetric ζ(4) case.
## The second digit (the user's question: "is Wieferich type the stand-out case?")
The mirror kills the p⁻⁸ level. The exponent is 7 unless the next digit X(n,p) = (p⁷A₀) mod p vanishes too. Data: exact forms n = 5…160, every single-digit prime; zu2003_digit2_5_160.json.
- 1872 cases have the mirror prediction 7: 1777 are exactly 7, 95 go deeper, 0 are worse.
- Two kinds of deep drop:
  - whole rows in a = ⌊n/p⌋: (p, a) = (7, 1), (11, 2), (37, 2);
  - whole columns in b = n mod p: e.g. p = 19 at b = 4, 8 for every a; p = 23, b = 4; p = 37, b = 4, 11; p = 41, b = 20, 22, 23; p = 11, b = 1.
  - Where a row and a column cross, the exponent drops twice, to 5: n = 23 at p = 11; n = 78 at p = 37.
- Rank-one law: X(ap+b, p) ≡ F_p(a)·G_p(b) mod p holds in 118,846 of 118,846 quadruple checks.
- The factors are p-independent rationals: CRT plus rational reconstruction over the primes, with held-out primes predicted.
- **Lucas law: p⁷·A₀(ap + b) ≡ A₀(a) · ũ_b / 30 (mod p)**, where A₀(a) is the rational part of F̃_a and ũ_b = A₇ is the ζ(7)-coefficient of F̃_b:
  - A₀(1) = 3633 = 3·7·173;
  - A₀(2) = 375268245/2⁷ = 3·5·11·37·61469/2⁷;
  - ũ₁…ũ₄ = 1320, 1260360, 2669460000, 8604317181000 = 30·(44, 42012, 88982000, 286810572700);
  - ũ₄ = 2³·3·5³·19·23·37·177383.
  - These match the reconstructed R(a, b) exactly.
- So a prime cancels one digit deeper exactly when it divides the form at one of n's two base-p digits: the rational part at the top digit, or the ζ(7)-coefficient at the bottom digit. The digits of n select smaller copies of the same form.
- The shape is Wieferich-like (a congruence one digit deeper than forced), but no Fermat or Wilson quotient enters. By Lucas mod p², C(ap+b, cp+d) ≡ C(a,c)C(b,d)(1 + p[aH_b − cH_d − (a−c)H_{b−d}]); the quotients cancel in binomial ratios. 11 being base-3 Wieferich is a coincidence: 11 | A₀(2) and 11 | ũ₁.
- Out-of-sample test: 173 | A₀(1) predicts exponent 6 at p = 173 for n = 173…176 (a = 1, b = 0…3). Confirmed at all four, with every other prime at those n at 7 (zu2003_digit2_173_176.json). The next row prime is 61469 (a = 2), out of computational reach.
- Lenses. Base: the base-p digits of n pick the sub-forms. Palindrome: the mirror kills the first digit. Carries: arcs and Φ̃. Hyperoperation: a level of the p-adic tower is read by the form one level down. Literature on this Lucas law for rational parts of linear forms was not checked. Lucas congruences for Apéry-like numbers (Gessel 1982; Malik–Straub, mod p²) are the nearest known relatives.

# The Apéry-like digit law (2026-09-28, later; apery_lucas.py, apery_digit_law.py). The user: "Back to apery like I expected. That's exactly the direction."
## Zudilin's F̃_n = ũ_n ζ(7) + w̃_n ζ(5) − ṽ_n (exact, n ≤ 180, every prime p ≤ n, p = 2 and 3 included; apery_lucas_180.json)
1. U(n) = ũ_n/30 is an integer for every n ≤ 180: 1, 44, 42012, 88982000, 286810572700, 1177434247946544, … (not in the OEIS).
   - Closed form (derived from the partial fractions and checked at n = 1): U(n) = (−1)ⁿ Σ_j C(n,j)⁶C(n+j,n)³C(2n−j,n)³·[1 + (3/2)(n−2j)(H(2n−j) − H(n+j) + 3H(j) − 3H(n−j))].
   - Lucas property U(n) ≡ ∏ U(n_i) (mod p) over the base-p digits: 4155 checks, 0 failures.
2. Digit law for the rational part: p^{7L}·A₀(n) ≡ A₀(n_L)·∏_{i<L} U(n_i) (mod p), with n = Σ n_i pⁱ and L = v_p(D_n).
   - 4155 cases: exponent never below −7L, 0 mismatches.
   - All 1620 deeper cases are predicted by a vanishing factor.
   - This gives Zudilin's D_n⁷ at EVERY prime (7 per base-p digit), small primes included; the p² > 3n + 2 local-mirror analysis was the L = 1 case.
3. Zeros of the Lucas factor:
   - U(b) ≡ 0 mod p² for every 2p/3 ≤ b < p (1057 cases v_p = 2, 10 cases v_p = 3). This is a supercongruence, and it is Φ̃'s double drop.
   - The boundary b = (2p−1)/3 (p ≡ 2 mod 3, where the carry-free arc has no integer) also has v_p = 2.
   - Accidental simple zeros: (11,1), (19,4), (19,8), (23,4), (37,4), (37,11), (41,20), (41,22), (41,23), (53,32), (71,32), (73,10), (79,39), … These are the "stand-out" columns.
## Apéry's own approximations (apery_digit_law.py, n ≤ 150, every prime p ≤ n)
- ζ(3): a_n = Σ C(n,k)²C(n+k,k)², b_n = Σ C(n,k)²C(n+k,k)²[H₃(n) + Σ_{m≤k}(−1)^{m−1}/(2m³C(n,m)C(n+m,m))].
  - Lucas for a_n (Gessel 1982): 3009/3009.
  - Digit law p^{3L}b_n ≡ b_{n_L}·∏_{i<L} a_{n_i} (mod p): 3009 cases, 0 below −3L, 0 mismatches, 596 deeper, all predicted.
- ζ(2): a_n = Σ C(n,k)²C(n+k,k): the same law with exponent 2L. 3009 cases, 0 failures, 414 deeper, all predicted.
## Reading
The rational part of an Apéry-like approximation has p-adic leading term p^{−wL} · (rational part at the top base-p digit) · (Apéry numbers at the lower digits). Equivalently, the p-adic top of the approximant b_n/a_n is p^{−wL} times the approximant b_{n_L}/a_{n_L} at the top digit. The approximations are self-similar under n ↦ ⌊n/p^L⌋, one level of the p-adic tower at a time.

Denominators: w per base-p digit, which is the D_n^w of the classical proofs. A prime removes powers exactly when it divides a sub-form at one of n's digits. The carry region b ≥ 2p/3 is a forced zero (a supercongruence mod p²); the rest are accidents, like the carriers of the β-work.

Lenses:
- hybrid base: the base-p digits of n pick the sub-forms;
- palindrome: the local mirror kills the first digit at large p;
- carries: the zeros;
- hyperoperation: each p-power level is read by the form one level down.

Literature for the digit law not yet checked. Nearest: Gessel 1982 (Lucas for a_n); Beukers 1985 / Coster (congruences and supercongruences for a_n); Krattenthaler–Rivoal 2007 (denominators). Open: a proof (Lucas + the p-adic top of the partial sums H_w(n) ≈ p^{−wL}H_w(n_L)); the Catalan forms (Zudilin's Apéry-like forms for G, whose K¹ miss had a 2-adic part of 2.77 per n); the r = 5 tower.
## Modulo p² and beyond (the user: "Does it satisfy p²?"; apery_p2.py, scratchpad coster2.py; p ≥ 5)
- General digits: no. The digit law p^{wL}b_n ≡ b_{n_L}∏a_{n_i} holds mod p but lifts to p² only sometimes, for L = 1 and last digit ≠ 0:

| family | law lifts to p² | Lucas for a_n lifts to p² |
|---|---|---|
| Apéry ζ(3) | 291/2354 | 225/2354 |
| Apéry ζ(2) | 271/2354 | 264/2354 |
| Zudilin F̃ | 1073/3317 | 1094/3317 |

  This is the same behaviour as the Apéry numbers' own Lucas congruence (cf. Malik–Straub 2016, "Lucas congruences for the Apéry numbers modulo p²").
  - Mod p² Lucas on binomials, C(ap+b, cp+d) ≡ C(a,c)C(b,d)(1 + p[aH_b − cH_d − (a−c)H_{b−d}]), gives for Apéry's a_n (T(n,k) = C(n,k)²C(n+k,k)²):
    a_{ap+b} ≡ a_a a_b + 2p[a·a_a·E(b) + F(a)·G(b)], with E(b) = Σ_d T(b,d)(H_{b+d} − H_{b−d}) = ½ a′(b), F(a) = Σ_c c·T(a,c), G(b) = Σ_d T(b,d)(H_{b+d} + H_{b−d} − 2H_d).
  - So the second digit is carried by derivative ("second solution") sums.
- Lower digits zero, n = m·p^r (m < p): yes, and to p^{3r}. The rational parts satisfy the twin of Coster's supercongruence:
  p^{wr}·b_{mp^r} ≡ p^{w(r−1)}·b_{mp^{r−1}} (mod p^{3r}).
  - Held in 109/109 cases (ζ(3), n ≤ 150), 109/109 (ζ(2)) and 131/131 (Zudilin F̃, n ≤ 180). Excess over 3r is 0 in most cases (sharp).
  - At r = 1: p^w b_{mp} ≡ b_m (mod p³).
  - The a_n side, a_{mp^r} ≡ a_{mp^{r−1}} (mod p^{3r}), is Coster's theorem for Apéry's numbers. It also holds for Zudilin's U(n) = ũ_n/30: 131/131, a new sequence with the Coster supercongruence.
## The rational-part supercongruence, stated in general (scratchpad coster4.py; the user wrote it as v_p(p^{wr}b_{mp^r} − p^{w(r−1)}b_{mp^{r−1}}) ≥ 3r)
General form, every n divisible by p, any cofactor m (with L(n) = ⌊log_p n⌋, the number of base-p digits minus one):
  v_p( p^{w·L(n)}·b_n − p^{w·L(n/p)}·b_{n/p} ) ≥ 3·v_p(n).
For n = m·p^r with m < p this is the user's formula (L(n) = r). For m ≥ p the scaling must use L, not r. The earlier "failures" with p^{wr} scaling at m ≥ p were artifacts of that.

Results (n ≤ 150 for Apéry, n ≤ 180 for Zudilin):

| family | p ≥ 5 | p = 3 | p = 2 |
|---|---|---|---|
| Apéry ζ(3) | holds, 144 cases | holds | holds |
| Apéry ζ(2) | holds, 144 cases | holds only mod p^{3r−1} | fails by 1–2 |
| Zudilin F̃ | holds, 179 cases | holds, excess ≥ 1 | fails by 1 only at n = 2^r; holds with excess ≥ 1 for m ≥ 2 |

The a_n side (Coster) holds at p ≥ 5 for all m in every family.

Source of the 3, as for Coster: Jacobsthal–Kazandzidis, C(ap, bp) ≡ C(a, b) (mod p^{3 + v_p(ab(a−b))}) for p ≥ 5, i.e. Wolstenholme's theorem. The exponent does not depend on the weight w (2, 3, 7): the cube belongs to the binomials, not to the ζ-value. The modulus is (p^r)³, the cube of the base-p scale. The small-prime failures follow the known weakening of Jacobsthal at p = 2, 3. Literature for the rational-part version not yet checked; a proof route is Coster's (formal groups / Dwork congruences) applied to the second solution.

# Designing weights: the Catalan check (2026-09-28, later; catalan_lucas.py, scratchpad cat_probe.py, cat_companion.py). The user: "We can create our own weights now right"
Design rules from the Apéry-like laws, for a well-poised weight whose leading numbers are Lucas:
- the denominator exponent is w per base-p digit (w = weight of the top ζ-value; the palindrome removes the naive extra power);
- rebates come at the digits where the weight's Apéry numbers vanish mod p, forced on whole regions by carries (Φ̃: b ≥ 2p/3), plus accidents;
- the mod-p³ supercongruence along n = mp^r is free (Wolstenholme).
The decay is not given by these laws; it comes from the saddle machinery.

Check on the β-lattice: Zudilin's Apéry-like forms for G (math/0201024, u_nG − v_n → 0, recurrence in zudilin_ledger.py), U_n = 2^{4n}u_n = 1, 28, 2596, 311536, 41759524, … (integers, n ≤ 150).
- Lucas U_{ap+b} ≡ U_aU_b (mod p) holds EXACTLY when the lower digit b < p/2: 1550 two-digit cases, 0 failures, for both p mod 4.
- It fails for b > p/2: 817 of 910 cases, with ratio neither ±1.
- This is the lattice doubling: for b > p/2, 2n = (2a+1)p + (2b−p), and the digits of 2n turn odd.
- Companion (reflection) law: U_{(p+1)/2+j} ≡ W_j (mod p) for every prime, with a fixed integer sequence W = 16, 1024, 111616, 14286848, 1982611456, 288848084992, 43474477907968, 6697633473101824 (j = 0…7). Reconstructed by CRT over the primes to 277–283 and confirmed on held-out primes up to 317. W_j = 16·2^{4j}·(rationals with 2-power denominators) from the 2^{4n} normalisation.
- So the G-forms' Lucas law runs on the digits of 2n, interleaving U (even digits) and a companion W (odd digits).
- Denominators: e_p(V_n) = 2⌊log_p(2n−1)⌋ in 4892 of 5044 (n, p) cases (the "no cancellation" recorded in the K¹ arena). U_b ≡ 0 mod p only 21 times in 1562 (b < p < 120): these forms have no forced-zero region, hence no Φ̃-type rebate.
- This is the Catalan reading of the K¹ miss's "lattice doubling 2.00 per n": the arithmetic lives on 2n's digits.

Open (design): the top-digit companion and the rational-part law on the doubled lattice. Then choose G-weights (half-integer poles, alternating sums) whose U and W have a carry-forced zero region, i.e. a Φ̃-type rebate for G, and run them through the saddle for the net.

# Step 1 done: the two-state carry Lucas law for Zudilin's Catalan forms (2026-09-28, later; catalan_doubled.py, catalan_twostate.py, catalan_twostate_v.py)
Setting: u_nG − v_n → 0 (Zudilin math/0201024), U_n = 2^{4n}u_n, V_n = 2^{4n}v_n. N = 2n with base-p digits N_i. The carries of doubling are c_{−1} = 0, c_i = ⌊(2b_i + c_{i−1})/p⌋ (b_i the digits of n), so N_i = 2b_i + c_{i−1} − p·c_i.
- Two families on the doubled lattice:
  - X₀ = 𝒰, the recurrence family: U_{M/2} at even M. At odd M, the half-coset solution y(x), x ∈ ½ + Z, of Zudilin's recurrence, started at y(½) = 1. The coefficient (2x−1)² vanishes at x = ½, so this solution is unique. Values 4, 256, 27904, 3571712, 495652864, …
  - X₁ = ℱ, the carried family: 1, 4, 32, 304, 3136, 34064, 382976, 4412608, 51782656, 616387984, 7420936192, 90176553152, 1104287432704 (M = 0…12). Its odd entries were read from two-digit N and its even entries from three-digit N, then reconstructed across primes as p-independent integers.
  - ℱ is not a solution of the recurrence, and none of these sequences is in the OEIS. ℱ/𝒰 → about 1.27 (4/π?) slowly.
- LAW for U: U_n ≡ ∏_i X_{c_{i−1}}(N_i) (mod p). Each digit of 2n takes the recurrence family if no carry enters it and the carried family if one does.
  - 113,656 checks (every even N ≤ 12000, every prime 5 ≤ p ≤ 73, 2–6 digits, 0–5 carried digits), 0 failures. The 4–6-digit cases are out of sample.
- LAW for V: p^{2L}V(2n) ≡ χ₋₄(p)^{L + c_{L−1}} · Y_{c_{L−1}}(N_L) · ∏_{i<L} X_{c_{i−1}}(N_i) (mod p), with L = #digits − 1.
  - Top digit: Y₀(even M) = V_{M/2}, the form's own rational part; Y₁ is the carried rational family at odd M, read per prime (0 inconsistencies). Its p-independence is not yet established: reconstruction from 16 primes fails, probably for size.
  - Catalan's character enters exactly as the top level of Σ(−1)^k/(2k+1)² predicts: 2k + 1 = p^L(2t+1) gives (−1)^k = χ₋₄(p)^L(−1)^t.
  - About 114,000 checks, 0 failures.
  - The exponent is never above 2L, and every deeper case coincides with a vanishing product.
- Reading. The "lattice doubling" of the K¹ miss is the doubled digit law, with the arithmetic on the digits of 2n. The carries of n ↦ 2n switch each digit between two Apéry-like families. The palindrome-odd constant G shows its character χ₋₄ only in the rational part and only through the top level. Zudilin's forms have essentially no zeros of X₀, X₁ mod p (U_b ≡ 0 in 21 of 1562 cases), hence 2 per digit of 2n with no rebate.
- Design target (step 2): a G-weight whose X₀ and X₁ (and Y) vanish on carry-forced regions of digits. Then the rebate is readable from carries, exactly as Φ̃ was for Zudilin's ζ(7) form.

# The recurrences of the new sequences (2026-09-28, later; guess_rec.py, write_recurrences.py -> recurrences.md). The user: "can you write the polynomials for the new integer sequences"
All recurrences were found by kernel mod primes plus exact lift, and verified on every available term. Full polynomials are in recurrences.md.
- U(n) = ũ_n/30 (Zudilin's ζ(7) form, integral, Lucas, Coster): order 4, degree 30, verified on 181 terms.
  - c₀ = −27(n+1)⁵(n+2)³(3n+2)³(3n+4)³·S(n) and c₄ = (n+3)⁵(n+4)⁹·S̃(n), where S, S̃ are degree-16 apparent-singularity polynomials with the same leading coefficient 467578191732948.
  - The leading terms give λ⁴ − 9264λ³ − 12116166λ² + 752300λ − 19683: Zudilin's characteristic polynomial (λ → −λ). This is the "cumbersome" recurrence he did not print.
- Catalan (P(x) = 20x² − 8x + 1, Q(x) = 3520x⁶ + 5632x⁵ + 2064x⁴ − 384x³ − 156x² + 16x + 7):
  - U_n = 16ⁿu_n: (n+2)²(2n+3)²P(n+1)U_{n+2} = 4Q(n+1)U_{n+1} + 256(n+1)²(2n+1)²P(n+2)U_n (Zudilin).
  - Half-coset family g(j) = U(2j+1) (odd digits, no carry in): (j+2)²(2j+5)²·½P(j+3/2)·g(j+2) = 2Q(j+3/2)g(j+1) + 256(j+1)²(2j+3)²·½P(j+5/2)·g(j). Verified on 95 terms.
  - Carried family, odd part f(j) = ℱ(2j+1): (j+2)²(2j+5)²P(−j−1)f(j+2) = 4Q(−j−2)f(j+1) + 256(j+1)²(2j+3)²P(−j−2)f(j). Verified on 48 terms (ℱ odd reconstructed to M = 95).
  - Carried family, even part h(n) = ℱ(2n): (n+2)²(2n+3)²P(−n−½)h(n+2) = 4Q(−n−3/2)h(n+1) + 256(n+1)²(2n+1)²P(−n−3/2)h(n). Verified on 19 terms; predicts integers, matching the reconstructed ℱ(38), ℱ(40).
- Reading: the carried family is Zudilin's recurrence with P, Q evaluated at NEGATIVE (reflected) arguments; the recurrence family has them at positive arguments. A carry of n ↦ 2n moves a digit to the reflected side of the recurrence. This fits the reflection congruence U_{(p+1)/2+j} ≡ g(0)g(j) found earlier, and it is the carries lens and the palindrome lens in one statement.

# Step 2, first screen: G-weights by carries (2026-09-28, later; multibeta.py search 3 3 40 net_ref)
Zudilin's 2018 very-well-poised construction at s = 3 gives forms in 1 and β(2) = G. Its carry count φ(x, y) is the doubled-lattice one (⌊2·⌋ − ⌊·⌋ terms), and net_ref includes the refined lemma, i.e. exact carry accounting per prime.
- Exhaustive integer search, η₀ = 3…40: best +1.316 per unit η₀ at η = (28; 7, 9, 11). Positive everywhere and linear in the scale, so the sign is a shape property.
- The rebates are large (η₀ = 28: naive 42 → Zudilin 12.69 → refined 8.69 per n), but the forms do not decay: closeness is positive (+28.1 per n at η₀ = 28). With three denominator blocks the numerator (t+1)_{h₀−1} wins.
- So for G alone this family is out. The family that decays, Zudilin's ₆F₅(−1) Apéry-like forms (decay 2.406), was closed in the K¹ arena at about +1.2 per unit parameter with exact off-direction checks.
- Why Zudilin's G-weight has no forced zeros (the two-state law's finding), made structural: its top coefficients contain half-integer binomials C(m − ½, r). For odd p:
  v_p(C(m − ½, r)) = v_p(C(2m, 2r)) − v_p(C(m, r)) + v_p(C(2r, r)).
  That is, the carries of doubling minus the carries of the half-size binomial. The negative term cancels the forced carries, so no digit region is forced to zero. In the ζ(7) all-equal form every carry entered with a plus sign (C(n,j)⁶C(n+j,n)³C(2n−j,n)³), hence the forced region b ≥ 2p/3.
- Design principle for our own G-weights: the carry count of the top coefficients must be a positive combination on a whole digit region, e.g. numerator blocks on the half-integer lattice (odd products, +doubling carries) against denominator blocks, while keeping (i) forms in 1 and G only (well-poised parity) and (ii) decay. That is the next construction to build and test.

# Building our own G-weights (2026-09-28, later; gdesign.py). The user: "Absolutely"
Machine: gdesign.py computes exact F = Σ_{t≥t₀}(−1)^tR(t) for any R with half-integer poles (partial fractions by fmpq series, centre-cancellation handled). It returns the β(i) coefficients, a numeric check, log|F|/n and the true primitive denominators per prime. Control: Zudilin's ₆F₅(−1) forms (a = 3n+1, h = n+½, h₄ = n+1) give decay −2.53…−2.55 (→ −2.406), height 5.4 at n = 14, NET +2.9 (→ +4.37).
Design "shift": Zudilin's ζ(7) all-equal template translated by ½. R = n!(2t+2L+N+1)·∏(t+m+½)/∏_{k∈[L,L+N]}(t+k+½)³, with half-integer roots m ∈ [0, L−1] ∪ [L+N+1, 2L+N].
- Top coefficients ∝ (N−2j)·C(N,j)³·C(L+j,j)·C(L+N−j,N−j): every carry with a plus sign, for any L (verified structure).
- The forms are PURE G: odd β-coefficients vanish, no β(≥4), numeric check to 65 digits.

| L | decay log\|F\|/n (n = 16) | height (n = 16) | NET (n = 16) |
|---|---|---|---|
| n (n = 14) | −2.59 | 9.03 | +6.44 |
| 0 (no roots) | −3.81 (growing) | 8.61 | +4.80 |
| n/2 | −4.46 | 9.53 | +5.07 |
| 1.5n | stops decaying | — | — |

Per-prime anatomy:

| form | exponent on (n, 2n] | exponent on (2n, 4n] | 2-adic |
|---|---|---|---|
| Zudilin, n = 16 | 2 | 0 | 60 |
| shift L = 0, n = 16 | 3 | 0 | 28 |
| shift L = n, n = 12 | 2 | 3 | — |

STRUCTURAL FINDING (why Catalan is harder than ζ in these families):
- Zudilin's exponent is 2 although his poles have order 3. The order-3 partial sums cancel in mirror pairs: the palindrome, the p-adic shadow of β(3)'s vanishing coefficient. This needs the sum to start at the mirror position. That in turn needs numerator roots at the INTEGERS t = −1…−n, which let the partial sums be counted from −n (length 2n).
- Integer roots at half-integer poles produce the half-integer binomials C(k−½, r), whose carry count has a minus sign: no forced zeros.
- Half-integer roots give all-plus carries (forced zeros) but cannot be absorbed into the start. Then the top order does not cancel (exponent 3 instead of 2) and the partial sums run to 4n.
- On the half-integer lattice the palindrome and the positive carries pull against each other. On the integer lattice (Zudilin's ζ(7) form) they did not: exponent 7 instead of 8 AND the forced region b ≥ 2p/3.
- A reflection cannot mix the lattices while keeping poles on half-integers (a half-integer centre maps half-integers to integers).
Next candidates: the quarter lattice (poles in both classes 1/4, 3/4 as in Zudilin's Apéry-like forms, where a reflection exchanges the classes); combinations of two sums; derivative forms with simple poles (pure β(2) from order-1 poles).

# Constructing G-forms from scratch by linear algebra (2026-09-28, later; gcon.py). The user: "We don't have to use zudilin form exactly. Use better modular linear algebra and construct the entire thing."
Constructor (gcon.py, no templates):
- R = P/Q on the quarter lattice. Poles at t = −(k+¼), −(k+¾) with any orders. P = (2t+c)·∏(prescribed zeros on (¼)Z)·Σ x_j (t(t+c))^j, with x unknown. An integer centre c makes t ↦ −c−t swap the classes.
- F = Σ_{t≥0} R(t) = Σ_i [A_i ζ(i,¼) + B_i ζ(i,¾)] − V. Pure G ⟺ A_i + B_i = 0 for all i and A_i − B_i = 0 for i ≠ 2, which is linear in x. Then F = 8(A₂ − B₂)G − V.
- Extra linear conditions: vanishing of top-order Laurent coefficients at chosen poles.
- Kernels by nmod_mat mod 62-bit primes, CRT, rational reconstruction and exact verification (kernel_modular). Numeric checks agree to 10⁻⁸⁴ … 10⁻¹⁵⁷.
Designs:
- "qz2" (Zudilin-like root placement on the quarter lattice) and "low2": NET +6.6 and +6.7 at n = 28 (Zudilin's forms +3.48 at n = 28). Decay and height both grow like log n: an unbalanced normalisation, which the scale-invariant NET absorbs.
- Apéry-style layout (poles k ∈ [0, N] in both classes, zeros inside the summation range at t = 0..M−1, mirror zeros at −c−i):

| N | M | order | NET at n = 20 |
|---|---|---|---|
| n | n/2 | 2 | +4.02 |
| n | n | 2 | +6.09 (= qz2 by translation) |
| n | 1.5n | 2 | +6.67 |
| n | n | 3 (ζ(3) killed linearly) | +9.8 |
| n | 1.5n | 3 | +10.4 |
| n | 2n | 3 | no decay |

- Mechanics-driven kills (Apéry N = n, M = n/2, order 2, free degree K; top-order coefficient killed at the K farthest class-¼ poles; modular kernels). n = 20: K = 0, 1, 2, 3, 5, 10 gives NET +4.02, +4.27, +4.23, +3.80, +4.11, +4.70. The height drops (13.8 → 10.1) and the 2-adic exponent goes from +30 to −8 (the quarter lattice can shed the whole 2-adic tax), but the decay weakens as much (−9.75 → −5.37). The odd exponents stay 2–3 on (n, 2n].
Finding: in every quarter-lattice layout with an integer centre, the odd-prime exponent on (n, 2n] stays 3, never the weight 2. The class-swapping mirror pairs a ¼-pole with a ¾-pole whose partial sums run over different progressions (4l+1 against 4l+3), so the top order cannot cancel termwise. Every lever (zeros, kills) moves decay and height together, and NET stays about +4.
Next principled design: a HALF-INTEGER centre. Then t ↦ −c−t maps each quarter class to itself (the top-order pairs share a progression, so they can cancel) and maps integer zeros to half-integer zeros: a mixed-lattice numerator with the symmetric start on the integer side.
Half-integer centre (d_half, c = N + ½; t ↦ −c−t keeps each quarter class, integer zeros ↔ half-integer zeros):
- Order 2, free degree 0: no pure-G weight (kernel dimension 0).
- Order 2, free degree 1: the unique solution has G-coefficient ZERO. The sum is a rational number (primitive value ±1, NET "0.000"), so the within-class palindrome kills β(2) itself.
- Order 2, free degree 2: kernel dimension 2, one direction of which is that degenerate rational form.
- Order 3: pure-G forms exist but are worse. NET +6.6 (n = 12, M = n/2) and +8.9 (n = 16, M = n), with odd exponents 3–5 on (n, 2n].
- gcon.measure now flags U = 0 as degenerate.
NO-GO in these symmetric single-sum layouts on the quarter lattice:
- integer centre (the mirror swaps the classes): pure G, but the top order never cancels (exponent 3 on (n, 2n]), NET about +4;
- half-integer centre (the mirror keeps each class): the mirror that could cancel the top order kills G itself at order 2, and at order 3 it gives worse forms.
The palindrome that cancels the top order is the same palindrome that decides which L-values survive, and on the quarter lattice the parity that keeps G is the one that does not cancel. On the integer lattice (ζ) the two parities coincide, which is why the ζ(7) form got both.

# Two weights by lattice linear algebra (2026-09-28, later; gcon2.py). The user: "Exactly what I'm thinking"
Construction:
- Q = both quarter classes k ∈ [0, N], order 2; P = ∏_{i<M}(t−i) · Σ x_j (t−M)^j. No symmetry is imposed, so the integer-centre and half-centre families are both subspaces.
- Step 1: pure-G conditions solved exactly (integer kernel W).
- Step 2: for every prime p where the unconstrained exponent exceeds 2, the p-adic digit of (U, V) at the lowest level must vanish. That is a congruence Σ y_b r_b(p) ≡ 0 (mod p) on the kernel coordinates.
- Step 3: LLL on [[I, K·Rᵀ], [0, K·diag(p)]], taking the shortest solution. (Siegel's-lemma style: the index is ∏p ≈ e^N, so the short vectors cost only e^{N/dim}.)
Results (N = n, M = n/2):

| n | reference NET | lattice NET, minimal freedom | lattice NET, extra freedom 2 |
|---|---|---|---|
| 8 | +3.77 | +5.03 | +6.37 |
| 12 | +4.35 | +6.13 | +6.83 |
| 16 | +5.12 | +7.35 | +8.18 |
| 20 | +4.89 | +7.07 | +7.64 |

- The congruences work: the targeted exponents drop (e.g. n = 12, exponents on (N, 2N] go from 2, 3 to 2), and the height falls 0.9–1.3 per n.
- But the decay falls about 3 per n (n = 20: −8.43 → −5.34). Most dangerous primes need a second p-adic level, so each costs more than one parameter.
- First run (binomial basis, only 1–3 congruences imposed): NET +6.4 against +4.3 at n = 12.
THE BALANCE LAW (the mechanics of freedom):
- Every free parameter spent on arithmetic (a kill or a congruence) weakens the decay by about as much as, or more than, it saves in height, roughly log n per parameter on each side.
  - Kill experiment, n = 20, K = 10: +8.8 decay loss per parameter against +7.4 height gain.
  - Lattice run: worse than par.
- Arithmetic bought with generic linear-algebra freedom does not pay. It has to come FREE, from structure: a palindrome that cancels a top order, or carries that force zeros. That is why hypergeometric weights work at all.
- Combined with the quarter-lattice no-go (the free palindrome and the parity that keeps G are incompatible in symmetric single sums), Catalan needs a new source of free arithmetic.

# Catalan's character in the summation, poles on the integer lattice (2026-09-28, later; gtwist.py). The user: "continue" (their "x-³" was taken as χ₋₃: not used here; χ₋₄ is used)
Construction: F = Σ_{t≥t₀} χ₋₄(t+σ)·R(t) with R = P/Q, poles at integers t = −k.
- With d = (σ − k) mod 4, the tail Σ_{t≥t₀} χ₋₄(t+σ)(t+k)^{−i} is L_d(i) − (partial sum over u < k + t₀), where L₀ = β(i), L₂ = −β(i), L₁ = −2^{−i}η(i), L₃ = +2^{−i}η(i), and η(i) = (1 − 2^{1−i})ζ(i), η(1) = log 2.
- Pure G: the β(i≠2) and η(all i) coefficients vanish. This is linear in the free numerator.
- The mirror t ↦ −c−t is compatible with the twist only if c ≡ 2σ (mod 4).
Results:
1. Zudilin-like root placement (roots −1…−G and the mirror, sum from −G): height only 1.6–2.6 per n and no prime above n (no lattice doubling), but almost no decay (log|F|/n −0.2…−0.6). NET +1.3…+2.2 at n ≤ 12, rising. The period-4 twist breaks the whole-line cancellation of the balanced weight.
2. Symmetric Apéry layouts:
   - even n (incompatible phase): kernel 0 or degree failure;
   - odd n with the compatible phase σ = (N+1)/2: order 2 has kernel 0, or U = 0 once there is freedom (the palindrome kills G again);
   - order 3, e = 2: pure-G forms, but NET +3.8 (M = n/2) and +4.1 (M = n) at n = 21.
3. BEST: no symmetry, order 2, poles k ∈ [1, n], zeros t = 0…n/2−1, free part of minimal degree 3 (exactly the three pure-G conditions: β(1), η(1), η(2)).

| n | 24 | 32 | 40 | 48 | 56 |
|---|---|---|---|---|---|
| NET, σ = 1 | +2.14 | +2.47 | +2.84 | +2.75 | +2.83 |
| NET, σ = 0 | +2.63 | +2.82 | +2.82 | +2.89 | +2.89 |

   - At n = 56 (σ = 1): decay −2.58, height 5.41, 2-adic 77 (= 0.95 per n against Zudilin's 2.77), and no prime above n.
   - Zudilin's forms at the same n: +3.91 (n = 40), +4.00 (n = 56), heading to +4.37.
   - Verified exact by direct summation: 23 digits (n = 8) and 37 digits (n = 12).
   - The remaining excess sits on the middle primes (√n, n/2], with exponents 4–9. Without symmetry the pole-distance denominators have no binomial cancellation.
Reading: moving χ₋₄ from the pole lattice into the summation removes the lattice doubling and most of the 2-adic tax: about 1.1 per n better than Zudilin at equal n (asymptotically about +2.9 against 4.37). The palindrome that would give binomial structure and top-order cancellation still kills G, even in the twisted setting. The parity no-go is robust across lattices and twists: for G, the free palindrome and the surviving L-value are always on opposite parities.
## Twisted construction refined by the character's support (2026-09-28, later; scratchpad gtwist_odd.py, gtwist_oddscan.py, gtwist_beta.py, gtwist_apery2.py)
- Apéry's ζ(2) weight ((t−1)⋯(t−n+1))²/((t)(t+1)⋯(t+n))², twisted by χ₋₄ (3 free parameters for pure G):
  - near-perfect arithmetic: height 0.8–1.1 per n, exponent 1 on every middle and large prime, 2-adic 4–13;
  - but the forms GROW (log|F|/n +4…+6.5); the twist destroys Apéry's cancellation. NET +6.8 (n = 32).
- ZEROS ONLY WHERE χ₋₄(t+σ) ≠ 0 (odd class over [0, span·n); general free part of minimal degree 3; poles at all integers [1, N]):
  - N = n, span = n: NET about +1.26 (σ = 1) and +1.43 (σ = 0) at n = 48.
  - Scan at n = 32/40/48: N = 0.75n, span 0.5n gives +0.82/+0.86/+0.92 (σ = 0) and +0.86/+0.96/+0.87 (σ = 1); per unit pole block about 1.15, against 1.27 at N = n and 1.36 at N = 1.25n.
  - The middle primes still carry 7–10 (no binomial structure without symmetry).
- POLES ONLY WHERE THEY GIVE β (k ≡ σ mod 2, so no η terms; pure G needs only β(1) = 0 and one free parameter):
  - poles k ∈ [1, n] of one parity, zeros on the support below n/2, σ = 1: NET +0.21, +0.47, +0.62, +0.89, +0.92 (n = 16…48);
  - the middle-prime exponents drop to 2–5, but per pole it is no better than the odd-zero design (0.92 per n/2 poles).
- Status: the best constructed G-forms have NET about +0.9 per n at n = 48, still rising slowly with n (the asymptote is not known). Zudilin's Apéry-like forms: +4.0 at n = 56, heading to +4.37.
- Mechanics learned:
  - (a) χ₋₄ in the summation instead of in the pole lattice: no lattice doubling, small 2-adic tax.
  - (b) zeros where the character vanishes are wasted, so put them on its support (same degree, twice the vanishing range).
  - (c) poles whose phase gives η are wasted, so put them on the β phase (fewer conditions, less freedom spent).
  - (d) the binomial structure (Apéry) gives height about 1 per n, but its decay needs the untwisted summation. Arithmetic structure and analytic cancellation sit on different sides of the twist.
- Open: the asymptote of the best design (saddle point for the decay, prime ranges for the height); a structure that is binomial within the character's support.
## (1) Asymptotics and (2) binomial structure on the support (2026-09-28, later; gtwist_asym.jsonl, scratchpad gtwist_hybrid.py)
(1) Exact computation to n = 160 (NET per n):

| n | 48 | 64 | 80 | 96 | 112 | 128 | 144 | 160 |
|---|---|---|---|---|---|---|---|---|
| odd-zero N = 0.75n, span 0.5n, σ = 0 | 0.924 | 1.076 | 1.132 | 1.167 | 1.210 | 1.207 | 1.266 | 1.252 |
| same, σ = 1 | 0.867 | 1.018 | 1.134 | 1.143 | 1.227 | 1.240 | 1.248 | 1.302 |
| β-poles K = n, span 0.5n, σ = 1 | 0.920 | 1.008 | 1.037 | 1.177 | 1.185 | 1.120 | 1.230 | 1.214 |

- Decay and height both drift like log n: odd σ = 0 decay −2.75 → −4.15, height 3.67 → 5.40.
- Height split at n = 160 (odd, σ = 0): 2-adic 0.96, small primes ≤ √n 1.74, middle 2.01, (n/2, n] 0.69, above n 0. The small and middle parts are factorial-type denominators (no binomial cancellation), matched by factorial-type decay.
- Extrapolation: A − C/n gives A ≈ 1.38; A − C/log n gives A ≈ 2.0. So the asymptotic NET is about +1.4…+2.0 per n: positive, and better than Zudilin's +4.37.
(2) In s = (t − t₁)/2 the twisted β-pole design is an alternating sum with half-integer poles (the β-phase) and integer zeros (the character's support): Zudilin's lattice geometry. So on-support zeros give smallness with mixed-sign (half-integer) binomials, and off-support zeros (the poles' own lattice) give plus-sign binomials but kill no term.
- Hybrids at n = 64 (on-support zeros below n/2 plus off-support zeros), NET against the base +1.01:
  - inside the summation range: +1.47 (0.25n), +1.40 (0.5n), +1.37 (0.75n);
  - beyond the poles: +1.67, +1.96, +2.30.
  - Middle exponents barely move (3–4 against 3–5).
- So the balance law again: off-support zeros are degree spent without smallness.
- SUPPORT NO-GO: on the character's support, zeros buy smallness and pay mixed carries; off it, they buy plus carries and no smallness.
Summary of the Catalan construction round:
- Three mechanical ideas took the net from +4.37 (Zudilin) to about +1.4…+2.0: the character in the summation, zeros on its support, poles on the β phase.
- Three walls remain: the parity no-go (palindrome against G), the support no-go (smallness against plus carries), and the balance law (freedom costs decay at par).
- Next idea: carry the three mechanical ideas back to the K² determinant arena, where Catalan's deficit was only +0.585 per K² and part of the arithmetic floor 8/3 + ⅙ ln 2 was the half-integer-lattice and 2-adic tax. There the integer-pole ζ(2) machine had floor ≈ 2.3 and net −0.8.

# The twisted Hankel machine: Catalan's character in the weight, poles on the integers (2026-09-28, later; twisted_hankel.py). The user: "Ok go"
Construction (t > 0, u = t²):
- Weight: w(t) = πt·sinh(πt)/cosh²(πt) = (2t²/π) Σ_{n≥0} (−1)^n (2n+1)/(t² + (n+½)²)².
  - The kernel has double poles at t = i(n+½) with signs (−1)^n = χ₋₄(2n+1). The character sits on the kernel's lattice ℤ+½; the poles sit on ℤ.
  - Moments: μ(e) = (2e+1)|E_{2e}|/2^{2e+1}. These are Euler numbers, von Staudt-free like Catalan's Genocchi moments; the only prime in their denominators is 2.
  - Pole values: φ(1/(u+a²)) = Σ_{n≥0} (−1)^n/(n+a+½)² = 4(−1)^a (G − β_a), with β_a = Σ_{k<a} (−1)^k/(2k+1)².
    - This holds for every integer a ≥ 0 (a = 0 gives 4G), and there is no rational part.
    - A pole plus a kernel point is always a half-integer, so no integer pole lands on the η (ζ(2)) phase: "poles on the β phase" is automatic.
  - Checked by quadrature to 40 digits: μ(0…5) and a = 0…6.
- Decay e^{−πt} and pole spacing 1 give ρ = π. By the ρ-law the closeness is −ln 9, the same as Catalan's.
- Engine: the profiles.py pipeline (shift-recurrence entries, exact determinants at h+1 points, Newton interpolation, full content factorisation, ball value at G).
  - The Catalan control inside twisted_hankel.py reproduces profiles.py at K = 40 exactly: closeness −2.1847, height 2.5639, net +0.3792, the same exponents.
  - Options: pole start a0, dropped head N with Fauzan copies, and support zeros Z(u)^r = ∏_{n<M}(4u + (2n+1)²)^r.
Results (N = 0, m = 0; per K²):

| K | net, twisted a0 = 0 | net, twisted a0 = 1 | net, Catalan | e₂/K², twisted | e₂/K², Catalan |
|---|---|---|---|---|---|
| 12 | +0.1047 | +0.1006 | +0.0996 | 0.146 | 0.361 |
| 20 | +0.2200 | +0.2010 | +0.2111 | 0.180 | 0.313 |
| 28 | +0.3169 | +0.3148 | +0.3117 | 0.168 | 0.264 |
| 36 | +0.3528 | +0.3539 | +0.3536 | 0.174 | 0.249 |
| 40 | +0.3782 | +0.3816 | +0.3792 | 0.180 | 0.246 |
| 48 | +0.4063 | +0.3974 | +0.3983 | 0.174 | 0.230 |
| 56 | +0.4319 | +0.4348 | +0.4322 | 0.173 | 0.221 |
| 60 | +0.4524 | +0.4530 | +0.4516 | 0.1725 | 0.217 |
| 80 | +0.4595 | +0.4587 | +0.4582 | 0.1755 | 0.209 |
| 100 | +0.4794 | +0.4838 | +0.4817 | 0.1732 | 0.200 |
| 120 | +0.4910 | +0.4914 | +0.4908 | 0.1719 | 0.194 |

- Further runs: K = 16, 24, 32, 44, 52 (all three machines) and K = 64 (+0.4555) and 96 (+0.4711) for the twisted machine, a0 = 0. The three nets agree to within ±0.01 at every K, and to 0.0006 at K = 120.
- Closeness at K = 40…120 (all tend to −ln 9 = −2.1972):
  - twisted a0 = 0: −2.1314, −2.1535, −2.1645, −2.1711, −2.1755 (from above);
  - twisted a0 = 1: −2.2142, −2.2087, −2.2058, −2.2041, −2.2030 (from below);
  - Catalan: −2.1847 … −2.1936.
- Height at K = 120: 2.6664 (twisted a0 = 0), 2.6944 (a0 = 1), 2.6844 (Catalan). P > 0 in every run.
Per-prime ledger at K = 60 (twisted a0 = 0 against Catalan):
- Top range K < p < 2K: e_p = 3·#{a ≥ (p+1)/2} exactly, for example p = 61: 87 and p = 113: 9. That is Catalan − 1 at every prime, because there is no fixed-point pole b = p (whose rational part 1/(4b) costs Catalan the +1).
- Middle range 2K/3 < p < K: identical (p = 41…59).
- Plateau: a flat 2K = 120, against Catalan's 2K − 2 with small cancellations (115–118).
- Primes 3, 5, 7: 417, 240, 226 against 388, 226, 217.
- Totals per K²: odd primes 2.453 against 2.429; 2-adic 0.120 against 0.151; sum 2.573 against 2.579.
The 2-adic exponent:
- It is identical for a0 = 0 and a0 = 1 at every K.
- It is a genuine K² term: e₂/K² stays at 0.172–0.176 for K = 36…120. Integer-pole ζ(2), by contrast, falls 0.246 → 0.125 (K = 40 → 100), which is K log K.
- It oscillates with the binary digits of K:
  - K = 2^j: 0.156, 0.168, 0.172 (K = 16, 32, 64);
  - K = 5·2^j: 0.180, 0.180, 0.1755 (K = 20, 40, 80);
  - K = 3·2^j: e₂ = 400 at K = 48 and 1600 at K = 96, exactly (5K/12)², i.e. 25/144.
- Constant ≈ 0.1736 powers of 2 per K², i.e. 0.120 in log units, against Catalan's fitted ⅙ (0.116). Catalan's finite-K values (0.246 → 0.194 at K = 40 → 120) carry an extra K log K that the twisted machine does not have.
Numerators at K = 40 (twisted, a0 = 0):
- Support zeros (Z vanishes at the kernel's lattice points u = −(n+½)², the character's support):
  - r = 1, M = 2, 5, 10, 20, 30, 40: net +0.387, +0.444, +0.618, +1.238, +2.102, +3.005;
  - r = 2: +0.414, +0.555, +0.979, +2.331, +3.982, +5.684.
  - The denominators collapse (den 2.46 → 0.28 per K², the 2-adic exponent 288 → −40), but the intrinsic coefficient size grows faster (0.05 → 7.3). This is the balance law in the K² arena.
- Fauzan head copies:
  - N = 3 with one copy: +0.368 against +0.378, the known N/K ≈ 0.1 effect;
  - N = 3, 5, 8 with 1–6 copies: otherwise +0.39 to +2.16.
Reading:
1. The twisted machine is Catalan's machine in the dual geometry.
   - G needs a half-offset between the kernel's lattice and the pole lattice.
   - Catalan puts it on the poles (kernel ℤ with signs (−1)^n, poles ℤ+½); the twisted machine puts it on the kernel (kernel ℤ+½, poles ℤ).
   - A χ₋₄-kernel on ℤ with period 4 and poles on the β phase is the twisted machine rescaled by 2. So up to scale these are the only two lattice Heine machines for G.
2. Odd primes cannot tell the two geometries apart: ½ is a p-adic unit, so ℤ and ℤ+½ are the same lattice p-adically.
   - The mirror pairs t + t' ≡ 0 (mod p) (from u = t²) and the partial-sum threshold t = p/2 (G's odd denominators) fall in the same places in both.
   - The threshold p/2 is the centre of every mirror pair (t, p − t), so every pair straddles it and costs 2 + 1 = 3. This is the arithmetic content of G's "lattice doubling", and it does not depend on the lattice.
   - Integer-pole ζ(2) has its threshold at p: both members of a pair lie below it and each pair costs 2 − 1 + 1 (von Staudt) = 2.
3. The prime 2 sees the offset, but its cost is conserved when the offset changes hands.
   - Catalan pays in the Vandermonde of its odd numerators; the twisted machine pays in the moments (Euler numbers over 2^{2e+1}).
   - The K² rate is ≈ ⅙ either way (0.1736 against ⅙).
4. The asymptotic net is 8/3 + ⅙ ln 2 − ln 9 ≈ +0.585 in both geometries; the measured nets agree to ±0.01 from K = 12 to 120.
   - The twist closes nothing, but it completes the classification: the K² arena is exhausted for G in both lattice geometries, not just within the half-integer class.
5. Why the three K¹ mechanisms of gtwist do not transfer:
   - (a) The character in the weight: in K¹ it removed the lattice doubling because a linear form has no t ↦ −t mirror. In K² the u = t² structure forces the mirror pairs, and G's odd denominators put the partial-sum threshold at their centre.
   - (b) Zeros on the support buy denominators at a loss.
   - (c) The β phase is automatic.
   - The 2-adic line is conserved, not removed.

# The 2-adic law of the lattice Hankel machines (2026-09-28, later; twisted_hankel.py e2, hankel_2adic_e2.json). The user: "12 times 2 plus 1 / 12² … 6 times 4 plus 1 / (6 times 2)² … this seems very very carry like getting pushed to the right of the decimal. Ramanujan summation again is what I'm seeing."
Prompt: the twisted machine had e₂ = 400 = (5K/12)² at K = 48 and 1600 at K = 96, i.e. 25/144 = (2·12 + 1)/12² = 1/6 + 1/144.
Tools:
- det(A + XB) = det(B)·charpoly(−B⁻¹A)(X). It gives the identical polynomial, 15× faster at K = 60, and is now the default in twisted_hankel.py. The controls zeta2int and zeta2half were added.
- Where the 2-adic content sits:
  - For the twisted, Catalan and integer-pole ζ(2) machines, the minimum 2-adic valuation over the coefficients of det(A + XB) is at X⁰. At K = 8, 12, 16, 24 the valuations rise smoothly up to v₂(det B) = 2K at the top.
  - For half-integer ζ(2) the minimum is in the middle of the polynomial.
  - So e₂ = −v₂(det A), where A is the Hankel matrix at X = 0: the regularized moments plus the reflected finite sums φ₀(1/(u+a²)) = −Σ_{n=−a}^{−1} (−1)ⁿ/(n+a+½)². CLI: `twisted_hankel.py e2 MACHINE K1,K2,…`.
  - [CORRECTED later on 2026-09-28] This was first written as "checked at every K = 2…128", but only 13 values of K had been checked.
    - The full comparison (section "Proof route, first stage") gives the content equal to −v₂(det A) at every K ≤ 128 except K = 11, 61 and 127, where the content is exactly one power of 2 larger (its minimum sits at a higher power of X).
    - In general e₂(content) ≥ −v₂(det A).
    - The e₂ values for K ≥ 160 in hankel_2adic_e2.json are −v₂(det A).
Data (exact; hankel_2adic_e2.json): every K = 1…128 for twisted, Catalan, zeta2int and zeta2half; K = 160, 192, 224, 256 for twisted and Catalan.
Results:
1. The K² constant is ⅙ for both geometries, not 25/144.
   - Least squares on K = 32…256 (basis K², K log₂K, K, 1): twisted 0.1679, Catalan 0.1676.
   - With K and 1 only (K = 24…112): 0.16623 and 0.16659.
   - Controls: half-integer ζ(2) ≈ ½ (0.48–0.50); integer-pole ζ(2) ≈ 0 (K log K).
2. Exact law on K = 3·2^j:
   - twisted: e₂ = (3/2)·4^j + (j − 3)·2^j = K²/6 + (K/3)·log₂(K/24). Exact at j = 4, 5, 6 (400, 1600, 6336); j ≤ 3 is off by 1–2.
   - Catalan: e₂ = (3/2)·4^j + (j + 5)·2^j + 1. Exact at j = 3…6 (161, 529, 1857, 6849).
   - Prediction for K = 384: 25088 (twisted) and 26113 (Catalan), running.
   - Hence e₂/K² = 1/6 + (j − 3)/(9·2^j). The correction is 1/144 at j = 4 and again at j = 5 (1/16 = 2/32), 1/192 at j = 6, and 1/288 at j = 7 (predicted). The user's 1/144 is this carry term: it coincides at two consecutive doublings and then decays like (log K)/K.
3. The first differences d(K) = e₂(K) − e₂(K−1) are constant on runs and jump at binary carries of K.
   - Twisted: d = 14 for K = 41…48 and 30 for K = 84…96, then 34 at K = 97.
   - Catalan: 17 and 33 on the same runs.
   - So the lower-order term is a carry count: K log K on the 3·2^j family, with coefficient ⅓ per binary level.
4. Catalan − twisted = (8/3)K + O(1) at every K = 8…256: 125 values, offset in [−3.7, +2.7], exactly (8/3)K + 1 at K = 3·2^j. The twist saves 8/3 powers of 2 per pole, a linear saving; the K² constant and the carry term are shared.
5. The 2^j family converges more slowly: (e₂ − K²/6)/K = −0.17, 0.04, 0.33, 0.58, 0.90 at K = 16…256 (≈ 0.21–0.32 per doubling).
   - No exact law for general K was found among small combinations of digit functions. Searched by exact integer nullspace, up to 4 extra functions: Σ_{n<K}s₂(n), s₂(K), Σ⌊K/2^k⌋², Σ⌊K/2^k⌋⌊K/2^{k+1}⌋, residues mod 3.
   - So the general-K law needs more than digit sums of K.
Reading (established against conjectural):
- Established: ⅙ is the constant for both geometries, and the finite-K excess is a carry term (flat increments between binary carries; exactly (K/3)log₂(K/24) on 3·2^j).
- The Ramanujan-summation link is structural in two places.
  - (a) The moments are the regularized sums over the kernel's lattice. For the twisted machine μ(e) = (2e+1)|E_{2e}|/2^{2e+1} is the regularized value of Σ(−1)ⁿ(n+½)^{2e}, i.e. β(−2e) up to normalisation. The whole 2-adic tax sits at X = 0, in those moments and the reflected finite sums.
  - (b) Digit-sum (carry) terms are governed by Delange's theorem: Σ_{n<N} s₂(n) = (N/2)log₂N + N·F(log₂N), with F continuous and 1-periodic.
    - The mean of F is log₂√(2π) − 3/4 − 1/(2 ln 2) = −0.145599; checked numerically as −0.145596 over N ∈ [2^16, 2^17).
    - √(2π) = e^{−ζ′(0)} is Ramanujan's regularized ∞!, and the Fourier coefficients of F are ∝ ζ(2πik/ln 2).
- Conjectural reading, not derived:
  - ⅙ = ½·⅓, with ⅓ of a power of 2 surviving per pair of poles (the ledger's Catalan figure, now the twisted one too).
  - ⅓ = 0.0101…₂ = Σ_{k≥1} 4^{−k} to the right of the binary point, while 2-adically −⅓ = 1 + 4 + 16 + ⋯ = …0101₂, the Ramanujan value of that divergent series: the same "01" carry pattern mirrored across the point.
  - Survivors per pair: 1 = 0.111…₂ (half-integer ζ(2)); ⅓ = 0.0101…₂ (both G machines: every other binary level, the alternating sign being the second binary digit); 0 (integer ζ(2)).
- Structural (algebra, not yet a theorem about minima):
  - In the Lagrange basis, det M = det(Δ″ + C′)/det(V)², with C′ = WᵀTW (W the Vandermonde rows, T an anti-triangular Hankel matrix of the moments). So e₂ = 2v₂(det V) − v₂(det(Δ″₀ + C′)).
  - Cauchy–Binet expands the second term over which poles are discrete and which rows are continuous: the 2-adic version of the two-species gas of open item 2.
  - 2v₂(det V) = 2(K−1)² − 2Σ_{b<K} s₂(b) exactly (integer poles 0…K−1), so the Vandermonde brings Delange's sum in directly.
## K = 360, the user's test, and the pre-registered predictions (2026-09-28, later). The user: "K=360 would be natural?"
Why 360 separates the readings:
- 360 = 15·24 = 45·8 = 101101000₂.
- The constant reading gives 25K²/144 = 22500. The carry reading gives K²/6 + (K/3)log₂(K/24) = 21600 + 120·log₂15 ≈ 22069.
- The gap between them is K²/144 − carry term ≈ 431 at K = 360, against only 51 at K = 180.
Predictions were written before the runs:
- K = 384: 25088 (twisted) and 26113 (Catalan), from the 3·2^j law.
- K = 320, twisted: 17469–17477, from Catalan − twisted = (8/3)K + O(1).
- K = 360, Catalan: twisted + 960 ± 4.
- K = 360, twisted: ≈ 22085–22105 from the 45·2^j family.
Exact results (det A route; 1.5–2.3 h each):

| K | twisted e₂ | Catalan e₂ | check |
|---|---|---|---|
| 180 | 5592 | 6069 | difference 477 = (8/3)K − 3 |
| 320 | 17472 | 18326 | twisted inside the predicted 17469–17477 |
| 360 | 22080 = 360²/6 + (4/3)·360 | 23041 = twisted + (8/3)K + 1 | carry reading 22069 (+11); 25/144 reading 22500 (−420) |
| 384 | 25088 | 26113 | both exactly the 3·2^j law (j = 7) |

- e₂/K² at 360 is 1/6 + 1/270. The constant is ⅙, and the 1/144 was carry.
- The 3·2^j law is now exact at j = 4…7 (twisted) and j = 3…7 (Catalan).
Gap to the law K²/6 + (K/3)log₂(K/24), in units of K/3, by family:
- 3·2^j: exactly 0 from K = 48 to 384.
- 5·2^j: 1.063, 0.863, 0.376, 0.213, 0.063 (K = 20…320).
- 45·2^j: 0.460, 0.426, 0.293, 0.093 (K = 45…360).
- 2^j: 0.085, −0.290, −0.415, −0.673, −0.728 (K = 16…256).
So the law is exact on 3·2^j and is the attractor of the 5·2^j and 45·2^j families. Only pure powers of two drift below it: the pole set {0, …, 2^j − 1} is a complete residue system mod 2^j, a candidate for extra cancellation, not checked.
Catalan − twisted = (8/3)K + 1 whenever 3 | K (K = 3·2^j and 360), and (8/3)K + O(1) otherwise.
## Pieces for an exact 2-adic law: coherence (2026-09-28, later; scratchpad th_pieces.py, th_coherence.py). The user: "Do you have the pieces for exact?"
Exact reduction:
- Rows u^i and columns 1/(u − z_b) give det M(X) = ±det(diag(F_b(X)) + J), with J = V⁻ᵀLVᵀ and L the strictly lower-triangular Toeplitz matrix of the moments. Equivalently J_{ab} = φ((ℓ_a(u) − δ_{ab})/(u − z_b)), which involves moments only.
- So det(A + XB) ∝ charpoly(diag(β_a) − ¼·diag((−1)^a)·J)(X).
- In the scale v = 4t² the moment part is ½·ŴᵀUŴ, where Ŵ is the Vandermonde matrix at −(2a)² and U is the anti-triangular Hankel matrix of the odd integers ν_r = Σ_k e_k((2a)²)(2j+1)E_{2j}, with j = K−2−k−r.
Coarse data do not determine e₂ (twisted, K = 48, true value 400):
- Euler numbers → 1: 126.
- Euler numbers perturbed at the 2¹ / 2³ digit: 1034 / 982.
- (2e+1)|E_{2e}| → random odd, with the right powers of 2: 1057.
- β_a → 2^{v₂(β_a)}: 712. β_a perturbed at digit v+1: 703. β_a → random odd × 2^{v₂}: 1153.
- β_a → 0: det A = 0 (the moment part has rank K−1).
- Same pattern at K = 8…40.
- Integer-pole ζ(2) with random odd moments of the right 2-adic size: 178, 371, 643, 1355 at K = 16, 24, 32, 48, against the true 111, 204, 287, 500 (a K log K term).
- Incoherent data give ≈ ½K² to 0.6K². The true data cancel two thirds of that for G (⅙) and all of it at the K² scale for integer-pole ζ(2) (0). Half-integer ζ(2) keeps ½.
The coherence is 2-adic continuity:
- Partial sums: v₂(β_{a+2^N} − β_a) ≥ N + 2 (min over a ≤ 64, N = 1…9), and v₂(β_{2^N}) = 2N + 1 exactly (N = 1…9).
- Euler numbers: v₂(E_{2e+2^N} − E_{2e}) = N (min over e ≤ 40, N = 1…7), a Kummer-type congruence.
- So the moments (the regularized sums, β(−2e)) and the reflected partial sums behave as integrals of one 2-adic measure. That is the χ₋₄ (Kubota–Leopoldt/Euler) measure, and the handoff already records "the functional is a p-adic measure" at odd primes.
- The handoff's open item 3 note ("needs both plain pole values and von Staudt moments — each hybrid loses it") is the same phenomenon at odd p.
State:
- In hand: the location (det A), the exact Vandermonde (Delange sum), the integral model (U, Ŵ), and the coherence (one measure).
- Missing: the p-adic Heine step. Write det A as the 2-adic Gram determinant of that measure and evaluate its valuation by an equilibrium on the support (p-orderings / Bhargava's generalized factorials of the square nodes).
- This is the 2-adic twin of the real-side computation that gave U(ρ) = −E(n). If it works it applies at every prime, so the whole ledger would become one log-gas per place.

# The exact 2-adic structure of the twisted machine (2026-09-28, later; twoadic_law.py; scratchpad th_padic1–9.py). The user: "Sounds good"
Why X = 0 is the 2-adic reading:
- At X = 0 the pole values are −4(−1)^aβ_a. These are the 2-adic limits of the same kernel sums Σ_{n≥0}(−1)ⁿ/(n+a+½)², because the full-period partial sums vanish 2-adically: v₂(β_{2^N}) = 2N + 1.
- The moments are the Ramanujan values. For p = 2 the Teichmüller character is χ₋₄, so E_{2e}/2 = L(−2e, χ₋₄) = ζ₂(−2e), a value of the Kubota–Leopoldt 2-adic zeta function (the Euler factor is 1 because χ₋₄(2) = 0).
- So det A = P(0) is the 2-adic twin of the Gram determinant P(G).
Ingredient laws (exact):
- Pole data F(a) = −4(−1)^aβ_a:
  - v₂(ΔⁿF(0)) = n + 1 for n = 1…63.
  - Divided differences over consecutive square nodes z = −a²: v₂([z_r, …, z_{r+m}]F̃) = 1 + s₂(m) for every start r and every m ≥ 1. 0 violations in 820 entries (N = 40).
  - v₂(F(r)) = 2 for r odd and 3 + 2v₂(r) for r even.
- Moments: the Mahler coefficients of (2e+1)|E_{2e}| have v₂ = n (n ≤ 40). The moment functional on Newton products over consecutive squares, ω(r, m) = φ₀(∏_{l=r}^{r+m−1}(u + l²)), has v₂ = −(2m+1) for every r.
- Nodes: the natural order is a 2-ordering of the squares (greedy check), and the 2-sequence is v₂((2k)!/2), Bhargava's factorial of the squares.
Frames that fail, and why:
- Lagrange/configuration expansion: at K = 12 the smallest term has v₂ = −67, attained by 8 subsets, while v₂(det A) = −21. Even the 1×1 minor Σw_a has v₂ = 5 against weights of v₂ = −23 (K = 16). This is the divided-difference cancellation from the smoothness of F.
- Mixed Newton basis (top nodes × bottom nodes): the entries are exactly the divided differences, the pole values on the anti-diagonal, and ω below it, with det G = det A. But the tropical determinant is attained by (K/2)! permutations, so the leading digits cancel.
- The invariant that resolves the ties is the 2-adic Smith normal form.
THE STRUCTURE. The elementary divisors of A over ℤ₂ are exactly:
- moment species: ε_n = −(2K−3) + 12n − 4s₂(n) = −(2K−3) + 8n + 4v₂(n!), for n = 0…N−1, with N(K) = #{n : 12n − 4s₂(n) < 2K − 3}. All are odd and negative.
- ones: K − 2N(K) divisors equal to 1.
- partners: N(K) divisors ≥ 2 (mostly 2 and 3, with 4, 5, 6 near the transitions N → N + 1). From K = 40 to 47 each new pole turns one partner from 2 into 3.
Verified exactly by the Smith normal form at K = 4…56, 60 and 64: 0 failures.
- The formula e₂(K) = Σ_{n<N}[(2K−3) − 12n + 4s₂(n)] − (K − 2N) − Σpartners is consistent with all 133 exact e₂ up to K = 384.
- The implied partner average is 2.00–3.91: exactly 3.00 at K = 48, 96, 192, 384; 3.67, 3.83, 3.91 at K = 64, 128, 256; 2.43–2.49 at 160, 320; 2.62, 2.74 at 180, 360.
Consequences (derived, not fitted):
- ⅙: N ≈ 2K/12 = K/6 moment terms, each ≈ −(2K − 12n), sum to K²/6.
- ⅓·K·log₂K: 4Σ_{n<N} s₂(n) with N ≈ K/6 is Delange's digit sum, ≈ 2N·log₂N = (K/3)·log₂(K/6), with Delange's log-periodic fluctuation (mean log₂√(2π) − 3/4 − 1/(2 ln 2), where √(2π) = e^{−ζ′(0)}).
- The 3·2^j law: partners all equal to 3 there give exactly (3/2)4^j + (j−3)2^j. Predicted: K = 768 gives 99584.
- The powers-of-two anomaly: partners tend to 4 there, so e₂ sits below the smooth law.
The user's reading, made literal:
- Each moment term steps by 12 and gives back 4 per binary digit of its index (the carries). There are N ≈ 2K/12 of them, which is where ⅙ comes from.
- Their digit sum carries the ζ-regularized fluctuation.
Open:
- The partner law. Its average is a log-periodic function of K: 2.5 at 5·2^j, 2.6–2.75 at 45·2^j, 3 at 3·2^j, → 4 at 2^j.
- A proof of the moment-species law (the p-adic Heine statement to prove).
- The same Smith normal form for the Catalan machine (expected: shifted moment species, since Catalan − twisted = (8/3)K + O(1)).
- Odd primes: the same method is the candidate for open item 3.
## Both G machines are continuous dual Hahn weights (2026-09-28, later; scratchpad th_go1.py, th_go2.py). The user: "Up 3 back 1. That's Pythagorean. Anyways. Go"
"Up 3, back 1": ε_n = −(2K−3) + 4(3n − s₂(n)). Each moment step goes up 3, and each binary digit of n takes 1 back (×3 then ÷2, the Pythagorean move, in base 2).
1. The moment species is intrinsic to the moments.
   - The pure moment Hankel [μ(i+j)]_{i,j<n} (twisted, u-scale) has 2-adic elementary divisors exactly −(4n−3) + 12k − 4s₂(k), k < n (checked n = 4, 8, 12, 16, 24).
   - The machine's moment species is its negative part with effective n = K/2.
2. The twisted weight is the continuous dual Hahn weight with a = b = c = ½:
   - π·t·sinh(πt)/cosh²(πt) = |Γ(½+it)|⁶ / (4π|Γ(2it)|²).
   - Monic orthogonal polynomials in u = t²: b_n = 2n² + 2n + ¾, λ_n = n⁴, norms h_n = (n!)⁴/2 (orthogonality and norms exact for n < 30).
   - In v = 4t² the recurrence is integral (8n² + 8n + 3, 16n⁴), so the v-scale moment Hankel has elementary divisors equal to its pivots, 8k − 4s₂(k) − 1 = 4v₂(k!) + 4k − 1 (checked n = 8, 16, 24).
   - The u-scale list is the v-scale list with the scaling applied in reverse order: ED^u_k = ED^v_k − 4(n−1−k).
   - So "up 3" = 4k from (k!)⁴ via Legendre, + 4k from the 16 in λ^v (the half-integer kernel lattice), + 4k from the reversed u↔v scaling; "back 1" = −4s₂(k) from Legendre.
3. Catalan's weight is continuous dual Hahn too, with (a, b, c) = (1, 1, 0):
   - By Γ(2it) duplication, π·t²·cosh(πt)/sinh²(πt) ∝ |Γ(1+it)²Γ(it)/Γ(2it)|².
   - The recurrence would be b_n = (2n+1)(n+1), λ_n = n³(n+1). Checks: b₀ = 1 and λ₁ = 2 from the moments.
   - The predicted pivots 4n + 3v₂(n!) + v₂((n+1)!) − 1 = −1, 4, 11, 17, 27, 32, … equal the computed Catalan moment-Hankel elementary divisors (v-scale).
   - So the two G machines are CDH(½,½,½) with integer poles and CDH(1,1,0) with half-integer poles. Moment integrality and pole integrality always sit a factor 4 apart: the half-offset.
4. The pole values of the orthogonal polynomials are explicit:
   - p_n(−b²) = (−1)ⁿ(n!)²·₃F₂(−n, ½−b, ½+b; 1, 1; 1), with v₂ = −2n for every b (n < 20, b < 30, 0 violations).
   - Hence P = L·U with L_{nk} = (−1)^{n+k}(n!)²C(n,k)/(k!)² and U_{kb} = ∏_{j<k}((j+½)² − b²) = 4^{−k}∏_{j<k}((2j+1)² − 4b²), all odd numerators.
   - This is the Newton basis at the kernel's lattice evaluated at the pole lattice: the two lattices meet in one matrix.
5. Uvarov form (exact, K = 4…10, X = 0 and 1): det M(X) = ±det[p_i(−b²)F_b(X) + R_i(−b²)]_{i,b<K} / det V, with R_i(z) = φ((p_i(u) − p_i(z))/(u − z)) the associated polynomials. The machine is a rational modification of a classical hypergeometric weight, which is the route to a proof.
6. The Catalan machine's own Smith form has a different, paired pattern (values in pairs 2 apart: −76, −74, −70, −68, … at K = 40), to be analysed with CDH(1,1,0).
Open:
- The partner law.
- A proof via the ₃F₂ structure (L·U with the cross-lattice products).
- Catalan's Smith form.
- The ζ(2) machines, which are Wilson (1, 1, ½, 0) by the same duplication.
## Proof route, first stage (2026-09-28, later; scratchpad th_proof1.py, th_proof2.py, th_detA_all.py). The user: "Ok"
Correction first:
- The content exponent e₂ equals −v₂(det A) at every K ≤ 128 except K = 11, 61 and 127, where the content is exactly 1 larger.
- The elementary-divisor law of the previous sections is a law for det A = P(0).
1. Exact factorisation (K = 4…16), in the v-scale:
   - det M^v(0) = ±det(F̂ + L̂). Use divided differences over the pole nodes w_b = −4b² on every row of the Uvarov matrix, then the Leibniz rule.
   - F̂_{rm} = [w_r, …, w_m]F̃ is upper triangular: the pole data, with diagonal F_m.
   - L̂ = A⁻¹R̂ is strictly lower triangular: the moments. Here A_{ir} = [w_0..w_r]p_i is unit lower triangular.
   - 2-adic pattern: v₂(F̂_{r,r+d}) = −(2d + 1 − s₂(d)) (Toeplitz); the diagonal is 0 for odd m and 1 + 2v₂(m) for even m; every entry of L̂ has v₂ = −1.
2. One infinite matrix.
   - The lower entries are L̂_{rm} = φ₀(∏_{j=m+1}^{r−1}(v + 4j²)) and the upper entries are φ₀(1/∏_{j=r}^{m}(v + 4j²)). So with N_k = ∏_{j<k}(v + 4j²):
     T_{rm} = φ₀(N_r(v)/N_{m+1}(v)),
     and no entry depends on K.
   - The K-pole machine is the K×K leading block: −v₂(det T_[K]) − 2K = −v₂(det A_K) for every K = 2…32.
   - So −v₂(det A_K) is a running sum of the pivot valuations of T: d(K) = −v₂(τ_{K−1}) − 2.
3. The flat runs are explained exactly.
   - With N fixed, each moment term gains 2 per new pole, so d = 2N − 1 − ΔΠ. For example, N = 8 and ΔΠ = +1 give 14 on K = 41…48; N = 16 gives 30 on 84…96; N = 12 with ΔΠ = −2 gives 25 on 65…68.
   - ΔΠ(K) for K ≤ 80 is in blocks: −1 on 17–20, +1 on 21–23, +1 on 41–47, −2 on 65–68, and so on, with jumps where N increases.
4. The pivots as Newton–Padé residuals.
   - The biorthogonal systems of T are: P_k, the degree-k polynomial with φ₀(P_k/(v − w_j)) = 0 for j < k; and the partial fractions at the nodes.
   - So −R[P_k]/P_k is the multipoint Padé (rational) interpolant of the Stieltjes function S(w) = φ₀(1/(v − w)) at the square nodes w_0…w_{k−1}.
   - The pivot is τ_k = ρ_k/∏_{j<k}(w_k − w_j), where ρ_k = P_k(w_k)F̃(w_k) + R[P_k](w_k) is the residual at the next node. Hence v₂(τ_k) = v₂(ρ_k) − 4k + 1 + s₂(k).
State:
- The whole 2-adic law of det A is reduced to one lemma: the 2-adic valuation of the Newton–Padé (Thiele-type) residuals ρ_k of the CDH(½,½,½) Stieltjes function at the squares.
- The moment species, the ones and the partners are its consequences, and the partner law is ΔΠ = 2N − 1 − d.
- Not yet proved.
## The residuals computed (2026-09-28, later; scratchpad th_lemma1.py, lemma1.json). The user: "Ok go"
- P_k is solved exactly for k = 2…48 from φ₀(P_k/(v − w_j)) = 0 for j < k.
  - Equivalently, P_k is the k-th orthogonal polynomial of the varying functional φ₀/N_k: the denominator of the two-point (k nodes + k moments at ∞) Padé approximant of S(w) = φ₀(1/(v − w)).
  - At k = 1 there is no solution, because F₀ = 0 at X = 0 (the a = 0 pole is silent), so the first step is a 2×2 block.
- The pivot identity is exact for every k = 2…48: d(k+1) = −v₂(τ_k) − 2, with τ_k = ρ_k/∏_{j<k}(w_k − w_j). This includes the K = 11/12 steps of det A.
- v₂(ρ_k) = 3, 7, 10, 13, 17, 20, 24, 27, 35, 33, 38, 42, … , 166, 169 (k = 2…48).
  - Increments are mostly 3 or 4 (mean ≈ 11/3), with jumps at carries.
  - v₂(ρ_k) = v₂(P_k(w_k)) (≈ 1.6k) + v₂(error at the next node) (≈ 2k).
  - The Newton coefficients of P_k have v₂ ≈ 1.6(k − i): P_k is not 2-adically close to N_k.
- No closed law for v₂(ρ_k) yet.
- The partners change one at a time, in runs:
  - K = 16…20: {3⁴} → {2,3³} → {2²,3²} → {2³,3} → {2⁴} (one 3 → 2 per step);
  - K = 21…23: one 2 → 3 per step;
  - K = 41…47: one 2 → 3 per step;
  - K = 65…68: ΔΠ = −2 per step, consistent with one 5 → 3 per step.
  - New partners (4, 5, 6) appear when N increases.
  - So each partner behaves like a two-level switch attached to one pole mode, toggled as K passes thresholds set by the binary structure. Labelling partners by pole is the next step (Schur complements / sub-block Smith forms).

# THE PARTNER LAW, PINNED DOWN: the complete 2-adic content law (2026-09-28, later; twoadic_law.py law; scratchpad th_frames.py, th_partner_search.py, th_halves.py, th_halfschur.py, th_midband.py). The user: "Go for the partner law. pin down the rule. try a bit more"
The law. With s₂ the binary digit sum:
- e₂(K) = −min_{0≤M≤K} T(K, M) − [the minimum is attained twice], where
  - T(K, M) = Σ_{n<M} ε_n + Σ_{m<K−M} κ_m;
  - ε_n = −(2K−3) + 12n − 4s₂(n), the moment species ("up 3, back 1");
  - κ_m = 1 + s₂(K−1−m) − s₂(m), the pole band (Kummer: s₂(m) + s₂(K−1−m) − s₂(K−1) is the number of carries in m + (K−1−m)).
- Exact for all 126 content values K = 3…128, and for every det-A value K = 160, 180, 192, 224, 256, 320, 360, 384 (134 values, 0 misses).
- It predicts 99584 at K = 768 and 396288 at K = 1536, which is the 3·2^j law.
How it was found:
1. Frames. The Smith form of the Uvarov matrix is one stable, K-independent sequence: −1, 0, 3, 7, 10, … = the residual valuations v₂(ρ_k). The scale-v and T frames agree exactly. The u-frame "partners" are a regrouping, so the partner law is the law of the per-pole increments.
2. The "set bits of N" rule for ΔΠ is exact only for N = 2, 4, 6, 8, 12, 16. The best block/phase triangle-wave rule reaches 56 of 104 steps, so it is not the structure.
3. First-half Kummer law: the natural-order pivots have v₂ = 1 + s₂(K−1−k) − s₂(k) for every k < (K−1)/2 (0 violations in 400).
4. Pivots read from the corner are exactly the moment species. The leading half-block has v₂(det) = K exactly (the first-half law sums to K). The Schur complement of the trailing half has elementary divisors {ε_n} ∪ {p_n − 2} ∪ zeros.
5. Eliminate the pole half [0, L) and the moment corner [K−N, K) together.
   - They decouple exactly: v₂(det A_EE) = K + Σε at every K tested.
   - The middle band's pivots follow the same Kummer law (exact at K = 12, 16, 20, 28, 32, 36, 40, 44; for example K = 40, m = 20…31 gives 2, 0, 0, −2, 3, 1, 1, −1, 1, −1, −1, −3).
   - So the pole band obeys κ_m for all m < K − N once the moment modes are out.
6. Summing gives F(K) = N(2K−3) − 6N(N−1) + 5D(N) − (K−N) − D(K) + D(K−N), with D(x) = Σ_{j<x} s₂(j) Delange's digit sum. This is exact at 117 of 135 K.
   - All exceptions sit at the moment/pole boundary: births (newborn depth ε = −1) and the steps just before births (incoming ε = +1 or +3).
   - The 2-adic min rule settles them: the smaller valuation takes the boundary slot. At births of n = 2^j the Kummer value κ = 2 − s₂(K−N) wins, which is the j − 2 seen at K = 48, 96, 192, 384.
   - At the 9 ties (ε = κ: K = 11, 24, 34, 61, 70, 82, 117, 127, 224) the units cancel. The content gains exactly one power of 2; det A gains two at K = 11, 61, 127, which is exactly where content ≠ det A.
Reading:
- The 2-adic tax is a tropical two-species equilibrium. Moment modes cost ε_n each and pole slots cost κ_m each, and the split M* = N(K) or N ± 1 minimises the total. The boundary mode (depth ±1) goes to whichever species has the smaller valuation.
- The partners are no separate species: they are the Kummer carries of the pole band as the Smith form regroups them, with Σ(p_n − 2) = Σ_{m<K−N} κ_m − K = D(K) − D(N) − D(K−N) − N in the generic case (10 at K = 32, 0 at K = 40, 4 at K = 16).
- Every constant is now accounted for:
  - ⅙ comes from N ≈ K/6 moment modes;
  - the (K/3)·log₂K carry term comes from 4Σ_{n<N} s₂(n);
  - the pole band adds −(K−N) − D(K) + D(N) + D(K−N): Delange sums again, the ζ-regularized fluctuation;
  - ties cost 1.
- Not yet a proof: the moment-species law is intrinsic (continuous dual Hahn Hankel), the pole-band law is the first-half Kummer law, and the decoupling and tie rule are verified, not derived.

# Unified law and Catalan; what is proved (2026-09-28, later; twoadic_law_proof.md; twoadic_law.py merge/claw; scratchpad th_catalan_law.py, th_catalan_split.py, th_merge.py). The user: "Prove all and check catalan"
Catalan check: the same tropical form, with the roles of the two species swapped.
- Catalan pole band (v-scale): π(K, k) = −2K + 4k + s₂(K−1−k) − s₂(k) = κ(K, k) + 4k − 2K − 1.
  - It holds for all pole slots, well past K/2: at K = 40, k = 0…31 gives −76, −74, …, −5, 1, 3, 7, 9, 18, 20, 24, 26, 32, 34, 38, 40.
  - The half-block determinant equals the formula sum at K = 8…48.
- Catalan moment pivots: the bare continuous dual Hahn (1,1,0) values h_q = 4q + 3v₂(q!) + v₂((q+1)!) − 1 = −1, 4, 11, 17, 27, 32, 39, 46, … These are K-independent.
- The Smith form is the union: at K = 40, 32 pole values + 8 moment pivots = the full Smith form, exactly.
- Content law: e₂^C(K) = −min_M[Σ_{k<K−M} π(K,k) + Σ_{q<M} h_q] − [tie]. Exact at all 134 values of K (3…128 content, 160…384 det A); `twoadic_law.py claw`.
- Catalan's content is not the twisted content plus an exact linear term: 3(e₂^C − twisted law) − 8K ranges over −11…8, rising by 1 per step in runs.
UNIFIED LAW (v-scale, both machines): the 2-adic Smith form of A^v_K is the sorted union of
- the pole band κ(K,m) + 4m − 2K + c for m < K − M* (c = 0 for integer poles, c = −1 for half-integer poles), and
- the machine's own continuous dual Hahn moment pivots h_q for q < M* (twisted 8q − 4s₂(q) − 1; Catalan as above),
with M* = the tropical minimiser.
- Exact multiset equality at K = 8, 11, 12, 16, 20, 24, 28, 32, 34, 36, 40, 44, 48 for both machines.
- At the ties (twisted K = 11, 24, 34; Catalan K = 12) everything matches except the top divisor, which is +1 (+2 at twisted K = 11).
- For the twisted machine the u-frame tropical function equals the v-frame one plus 2K for every M (proved). The u-frame "moment species / ones / partners" were a scale artifact: the v-frame union is the invariant description.
Proved today (twoadic_law_proof.md):
- Theorem 1: both weights are continuous dual Hahn, (½,½,½) and (1,1,0); recurrences and norms follow.
- Theorem 2: the moment-Hankel Smith form in the v-scale is the pivots (integral recurrence, hence unimodular).
- Lemma 3: Mahler smoothness v₂(ΔⁿF(r)) = n + 1, from G(a) + G(a+1) = −(2a+1)⁻² and Stirling expansion.
- The algebraic reductions: content ≥ −v₂(det A), the Uvarov form, the K-independent infinite matrix T, the pivots as Newton–Padé residuals, u↔v invariance.
- The closed form [z₀..z_m]F̃ = −(8/(2m)!)·Σ_j C(2m, m−j)β_j. With the identities Σ_j (−1)^j C(2m, m−j)C(j,n) = (−1)ⁿC(2m−n−1, m−n) and Σ_{p≤n} C(n+p, p)2^{−p} = 2ⁿ, the r = 0 divided-difference law becomes the congruence Σ_j C(2m, m−j)β_j ≡ 4^{m−1} (mod 2^{2m−1}). Verified, not proved: the k ≥ 1 terms cancel among themselves.
Not proved:
- Lemma A: the pole band, via the divided-difference law D_{r,m} = 2C(2m,m)·unit (0 violations for r, m ≤ 30).
- Lemma B: the merge (decoupling verified).
- Lemma C: the tie (+1 in the content).
- The Main Theorem is verified at every tested K with no failures.

# THE THREE LEMMAS PROVED: the 2-adic law is a pencil, Krawtchouk against continuous dual Hahn (2026-09-28, later still; twoadic_law_proof.md sections 8–12; proof_checks/). The user: "Prove the 3 lemmas."
Result: the Main Theorem (the unified Smith-form law, the content law with the tie rule, the det law) is proved for both machines and every K.
- The continuity step in Theorem 6(b), twisted only, was sketched at first; it is now written out in full, so no gaps remain. Unrefereed.

Five facts do all the work.

**1. The X-coefficient is Krawtchouk (Theorem 4).** B_K is the Hankel matrix of a symmetric binomial measure on the poles.
- Twisted: C(2K−2, K−1+y) on the integers y, pushed to u = −y².
- Catalan: C(2K−1, ·)·(2y)² on the half-integers.
- Its v-scale orthogonal polynomials are the even and odd Krawtchouk polynomials, with 2-integral coefficients.
- Norms ‖Q_n‖² = 2^N·n!·N!/((N−n)!·4ⁿ). Their valuations give the Kummer band exactly:
  - v₂ = κ(K,m) + 4m − 2K + 1 + c, with κ = 1 + s₂(K−1−m) − s₂(m);
  - the carries of m + (K−1−m) are the carries of the binomial norms.
- The Catalan c = −1 is the odd Krawtchouk family (a Christoffel factor (2y)²).

**2. The 2-adic Catalan constant (Theorem 5).** G(y) = (−1)^yβ_{|y|} splits as G₁ + C(−1)^y.
- G₁ = −Σ(−1)ⁿΔⁿg/2^{n+1} is 2-adically smooth.
- C = −G₁(0), with v₂(C) = −1.
- β_{2^N} → 0 is G₁(2^N) → −C.
- A(0) = −C·B + A(C) exactly.

**3. At X = C the pole values are the moments' own 2-adic Stieltjes transform (Theorem 6).**
- Catalan: F_j(C) = −Σ_k μ_k w_j^{−k−1}, with w_j = −(2j+1)² and v₂(μ_k) = 2k − 1 (the moments shrink, the nodes are units).
  - Proof by Euler–Boole: G₁ equals its derivative form −½Σ E_n(0)g^{(n)}/n!, and E_{2k+1}(0) ~ B_{2k+2} ~ the Catalan moments.
- Twisted: F_a(C)/4 = Σ_k (−1)^k m_k (1+4a²)^{−k−1}, with centred moments m_k = φ₀((v−1)^k) and v₂(m_k) = k − 1.
  - Here the moments live near v = 1 and the nodes near 0.
  - Proof: the Euler moments are Boole-regularised residue sums over the double poles of sech² at v = −(2m+1)². Applied to 1/(u+a²), this gives 2(h(a) + h(−a)) = −4G₁(a).
- Checked to 2^{−179} and 2^{−197}.
- Formulas for the 2-adic Catalan constant:
  - C = ½ + 2Σ(−1)^k μ_k = ½ + Σ(−1)^k (k+1)T_{2k+1} (tangent numbers);
  - C = Σ(−1)^k m_k: the twisted (−1)-st moment. This is the "μ(−1) = 4C" found earlier.
- The real identity F_j(G) = ∫w/(v − w_j) has an exact 2-adic twin, with G replaced by C.

**4. The moment side (Theorem 7).**
- A(C) = Hankel(φ₀·(1/N_K)), where 1/N_K is a unit power series.
- In the monic CDH basis, rescaled by the rational diagonal diag(4^{−n}) (Catalan) or diag(2^{−n}) (twisted), multiplication by v (resp. v − 1) is an integral matrix ≡ 0 mod 4 (resp. mod 2). (This line originally said "orthonormal CDH basis"; see the 2026-09-29 rigour pass below.)
- So A(C) = LΛLᵀ with L ≡ I and v₂(Λ_q) = h_q: the moment pivots, all K of them.

**5. The Pascal matrix (Theorem 8).**
- Catalan: Krawtchouk ≡ (v+1)^m mod 2 and CDH ≡ v^q mod 4.
- Twisted: Krawtchouk ≡ v^m and CDH ≡ (v+1)^q mod 2.
- Either way the transition is ≡ [C(m,q)] mod 2, and every corner minor is det[C(a+i, j)] = 1.

**Then Cauchy–Binet.**
- det(RΛRᵀ + Y·diag(β̂)) = Σ (det R[U^c, Q])²·ΠΛ_Q·Π(Yβ̂)_U.
- h and b are strictly increasing, so each Y-coefficient has a unique minimal term, with a unit coefficient:
  - v₂(q_{K−M}) = Σ_{q<M} h_q + Σ_{m<K−M} (b_m + 1) for EVERY M.
  - This was verified end to end: Catalan K = 6, 8, 12, 16, 20, 24 and twisted K = 6, 8, 11, 12, 16, ties included.
  - The apparent violations at large M seen earlier were artefacts of truncating C.
- Content (p(X) = q(X − C) with v₂(C) = −1): e₂ = −min_M T(K, M) − [tie].
  - T is strictly convex, so a tie is always two consecutive M, of opposite parity.
  - That forces the X¹ coefficient to m₀ + 1 exactly. This is the tie rule (Lemma C).
  - det A can gain +1 or +2 (the next 2-adic digit), which is where K = 11, 61, 127 come from.
- Smith form: the K smallest of {band} ∪ {moment pivots}. This is Lemma B: the "2-adic orthogonality" of the two species is the Pascal matrix mod 2.
- Lemma A is proved twice:
  - as the Krawtchouk band of B;
  - as the band of the X = 0 pole Hankel itself (Catalan: all m < K; twisted: m < ⌊K/2⌋). This second proof uses the Krawtchouk mirror (−1)^x·K_n = K_{N−n} plus the Mahler smoothness of G₁, with margin exactly v₂((k+1)(k+2)) ≥ 1.

Reading:
- The two species really are two measures, a binomial on the lattice of poles and the CDH weight at the 2-adic Catalan constant.
- The 2-adic tax is "take the K smallest of two increasing lists". That is the tropical minimum, and M* ≈ K/6 is where the lists cross.
- The constant ⅙, the (K/3)·log₂K carry term and the Delange fluctuations are now theorems about these two lists.
- The reason X = 0 is the 2-adic reading: A(0) = A(C) − C·B, and 2-adically C is a half-unit (v₂ = −1). The shift by the 2-adic Catalan constant is what the content sees.

# Rigour pass on the pencil proof: no square roots (2026-09-29; twoadic_law_proof.md sections 0, 8, 9, 12; proof_checks/). The user: "Gram matrix is D^{1/2}·π(J)·D^{1/2} … the entries of L are W_qq′·√(h_q/h_q′) … the individual √h_q and ratios √(h_q/h_q′) need not exist in Q₂ … Can you re write with more rigour."
The objection was correct.
- The old Theorem 7 factored the CDH Gram matrix through the orthonormal Jacobi matrix, with entries √λ_n, and asserted that L = D^{1/2}WD^{−1/2} is rational, integral and ≡ I.
- The square roots live in a ramified extension of ℚ₂, and integrality of L was not shown. Section 9 had the same device twice: "D^{1/2}(I + E″)D^{1/2}" and "the orthonormal Krawtchouk basis".

Rewritten, all over ℚ₂ with integral monic bases:

**Lemma 12.1.**
- Use the monic Jacobi matrix 𝒥 (subdiagonal 1, diagonal b_n, superdiagonal λ_n); then φ₀(f·P_pP_q) = h_p·f(𝒥)_pq.
- The rational diagonal S = diag(4^{−n}) (Catalan) or diag(2^{−n}) (twisted) balances it: S^{−1}(𝒥 − v*)S has entries
  - Catalan: 4, 4(2n+1)(n+1), 4n³(n+1);
  - twisted: 2, 2(2n+1)², 8n⁴.
  - So it is integral and ≡ 0 mod 4 (resp. 2).
- The moment decay v₂(μ_k) ≥ 2k−1 and v₂(m_k) ≥ k−1 falls out of the same matrix, which makes the 2-adic extension φ₀^{(2)} well defined with no appeal to von Staudt.

**Theorem 7.**
- G_pq = h_p·ρ^{q−p}·Π̃_pq, where Π̃ = Σ c_k𝒥̂^k ≡ c₀I (mod ρ).
- Gaussian elimination in ℤ₂ gives Π̃ = 𝐋𝐔 with 𝐋 ≡ I and 𝐔 ≡ c₀I.
- Conjugating back by the rational diagonals: L = H𝐋H^{−1}, with L_pq = (h_p/h_q)·ρ^{q−p}·𝐋_pq.
- The growth of h (v₂(λ_n) ≥ 5, resp. ≥ 4) beats the ρ^{q−p} loss, so v₂(L_pq) ≥ 3(p−q) + 2 (Catalan) or 3(p−q) + 1 (twisted). L is integral and ≡ I mod 16.
- Uniqueness of LDU turns L·U′ into LΛLᵀ with Λ_p = h_p·(unit).
- Hence R = TL ≡ T ≡ Pascal (mod 2) is now proved, not asserted.

**Section 9.**
- The criterion is Lemma 9.1, proved by determinant expansion: every off-diagonal factor beats the average valuation (δ_i + δ_σ(i))/2, and integrality turns "strictly greater" into "at least +1".
- The smoothness step uses balanced valuations β(M) = min v₂(M_ac) + ½(ν(a) − ν(c)). These are rational numbers attached to ℚ₂-entries; they telescope under products (β(MM′) ≥ β(M) + β(M′)), and they give exactly the old bounds.

**Also tightened.**
- Theorem 4 now states the explicit integer recurrences of the v-scale Krawtchouk polynomials:
  - twisted: b = −2(4n(K−1−n) + K − 1), λ = 4n(2n−1)(K−n)(2K−1−2n);
  - Catalan: b = −[(2n+1)(2K−1−2n) + 4(n+1)(K−1−n)], λ = 4n(K−n)(2n+1)(2K−1−2n).
- The twisted norm line had a slip. The u-scale valuation is ν(2n) + v₂(4/N!) = κ + 1; it had been written ν(2n) − 1 + v₂(4/N!). The stated band was right.
- Theorem 6(a) now spells out the convergence of the derivative form (Gauss norm on |y| ≤ r < 2), the functional equation (absolute convergence of the double series), and uniqueness (Mahler excess → ∞).
- Theorem 8 proves the Pascal reduction and det[C(a+i, j)] = 1 (row differences plus Pascal's rule).

**Checks.**
- th_thm7_rigorous.py (Catalan and twisted, K = 8, 12): all the rescaled integrality claims, G = h·ρ^{q−p}·Π̃, the mod-ρ elimination, the v₂(L) bounds, G = LΛLᵀ and v₂(Λ) = v₂(h).

# The readable paper; the 2-adic C is Calegari's L₂(2, χ₋₄); Lemma 3′ (2026-09-29)
The user asked: "a more readable md file everyone can read with proofs … for a reddit specific crowd … more of a proof paper … a single proof of Catalan … if anything isn't worked or shown, show it. If more work needs done let's do it." Then: "The obvious place to push is an irrationality proof." Then: "Go ahead and work on real G after write up."

**The paper.** catalan-2adic-paper.md, "The Price of the Prime 2 in Catalan's Constant".
- It proves one theorem for the Catalan machine: e₂(K) = −(sum of the K smallest of the two lists) − [tie].
  - band b_m = 4m − 2K + s₂(K−1−m) − s₂(m);
  - moments h_q = 8q − 3s₂(q) − s₂(q+1).
- It also gives the fingerprint (Smith form) and the growth K²/6 + (K/3)·log₂K + O(K).
- It is self-contained except for four classical facts, all cited: the continuous dual Hahn recurrence (and checked exactly for n < 30), von Staudt–Clausen, Cauchy–Binet, and Delange.
- It is explicitly NOT an irrationality proof. §12 gives the honest ledger (gain ln 9, odd-prime cost 8/3, prime-2 cost (1/6)·ln 2; net +0.585).
- A 10-line Python version of the theorem matches all 134 stored values.

**A simpler route to C.**
- C is defined directly by C := ½ + 2Σ(−1)^k μ_k = ½ + 1 − 4 + 48 − 1088 + ⋯ (tangent numbers).
- ρ₁(y) := −½Σ a_n(2y+1)^{−n−2}, with a_n = E_n(0)(−2)ⁿ(n+1).
- Boole's formula gives ρ₁(y) + ρ₁(y+1) = −g(y).
- Then β_y = C + (−1)^y ρ₁(y), and the 2-adic Stieltjes identity F_j(C) = Σ(−1)^kμ_k(2j+1)^{−2k−2} is a direct computation.
- This removes the Mahler-series uniqueness argument from the main line. In the technical file it stays in Theorem 6(a).

**Gap found and fixed (Lemma 3′).**
- In the proof file, §4 Lemma 3 asserted v₂(Δ^k g) ≥ k + v₂((k+1)!) from a termwise argument that only gives k + v₂(k!).
- New proof: Δ^k g(y) = (−2)^k·k!·[Σ_{i≤k} 1/(2y+2i+1)]/Π_{i≤k}(2y+2i+1), because g = −f′/2 with f = 1/(2y+1).
- The sum of reciprocals of k + 1 consecutive odd numbers has v₂ ≥ v₂(k+1): blocks of 2^t consecutive odd numbers are complete odd residue systems mod 2^{t+1}, and inversion permutes residues.
- Checked for y ∈ [−12, 12], k < 40. Lemma 3 itself only needed ≥ k + 1. The sharp bound feeds §9 (Lemma A) and Theorem 5's smoothness.

**Literature.**
- The 2-adic Catalan constant is known and proven irrational: Calegari, IMRN 2005 no. 20, 1235–1249 (arXiv math/0408214), Theorem 4.2, L₂(2, χ) irrational. Our C agrees with his approximation 783269/13060350 to 2-adic order 2³⁵ (he guarantees 2³⁴).
- Beukers (Acta Math. Sinica 2008, arXiv math/0603277) proves such p-adic irrationalities by Stieltjes continued fractions, the same objects as our CDH J-fraction.
- So the 2-adic analogue of "is G irrational?" is settled; the real G is open.
- Sun's claim (arXiv 2609.04176) has a critique: D. Wachs, arXiv 2609.22339 (15 Sept 2026), "the proof of the main theorem is incomplete". Its main quantitative gap is an unaccounted prime-2 contribution (≈ (2·ln 2)B²).

**Real-G push, step 1: the frontier, and Zagier's case E measured (caseE_denominators.py).**
- **The frontier.** Calegari–Dimitrov–Tang (arXiv 2408.15403) state it in Remark 11.1.17 (p. 165), for the two known holonomic families of rational approximations to G:
  - Zudilin/Rivoal/Nesterenko: the integrand peak is φ⁻⁵, the same geometry as Apéry's ζ(2), but the half-integral exponents give denominators 16ⁿ[1..2n]².
  - Zagier's case E: ODE singularities {0, 1/8; 1/4, ∞}, denominators [1..n]².
  - "Since e² > 16·(1/4) and e⁴ > ((1+√5)/2)⁵ (by a wide margin!), this definitely precludes an approach to the irrationality of the Catalan constant by our method using either of these particular families … unless some completely new idea is discovered."
- **In this program's net language** (net = log(denominator growth) − log(radius of the linear forms' generating function)), both comparisons read "net < ln 16 = 2.773":
  - Zudilin net = ln 16 + 4 − 5 ln φ = +4.367 (the +4.37 recorded here), which fails by 1.594;
  - case E net = 2 + ln 4 = 3.386, which fails by 0.614;
  - their successful L(2, χ₋₃) case (Zagier case C, r = 1) has net 2, a margin of 0.773.
  - To be confirmed against their holonomy theorem; a background reading of §2/§11–15 is running.
- **Case E measured.** (n+1)²u_{n+1} = (12n² + 12n + 4)u_n − 32n²u_{n−1}, with u₀ = 1, u₁ = 4 (integers) and v₀ = 0, v₁ = 1.
  - v_n/u_n → G/2 exactly: PSLQ gives [−2, 0, 1] against (1, G).
  - The true denominators ARE the full [1..n]². ln(den v_n)/n = 1.94 at n = 600, 98.2% of 2·ln lcm(1..n).
  - By class at n = 600: primes ≡ 1 mod 4 carry 98.4% of their bound; primes ≡ 3 mod 4 carry 100%.
  - The prime 2 never divides den(v_n): v_n is 2-integral. That is irrelevant to the growth.
  - Only sporadic primes fall short (61 and 137 at n = 600, with exponent 1 instead of 2).
  - No hidden savings, so the case-E gap of 0.614 per step is genuine.
- **Consequence for this program's best G-forms.** The gtwist designs (net +1.3 at n = 160, asymptote +1.4…2.0) are below ln 16. But their decay and height both drift like log n (factorial-type denominators), and their free parts are solved per n. So they are not G-function or holonomic families, and the CDT criterion does not apply as it stands.
- **The target** is therefore a HOLONOMIC family for G with geometric denominators and net < ln 16.
- **Survey of Zagier's second-order recurrences** (apery_like_survey.py, n = 400).
  - Setup: (n+1)²u_{n+1} = (an² + an + b)u_n − cn²u_{n−1}, with lim v_n/u_n identified by PSLQ, and net = τ(measured) − ln(1/|β|).
  - A: π²/24, net 1.95 (ln16 − net = +0.82).
  - C: L(2,χ₋₃)/2, net 1.95 (+0.83), CDT's successful case.
  - D: π²/30 (Apéry), net −0.45 (classical Apéry works).
  - E: G/2, net 3.33 (−0.56).
  - B (complex roots of equal modulus: v/u does not converge) and F (roots 9, 8): no relation with 1, G, π², π ln 2, ln²2, L(2,χ₋₃), … at height 10⁵.
  - So E is the only G-family among them, and it is the only failing case with an L-value limit. The rule "net < ln 16" separates C from E exactly as CDT's remark says.

**Checks.** proof_checks/th_paper_checks.py:
- W1: Lemma 3′;
- W2: Appendix A integrals to 40 digits;
- W3: the continuous dual Hahn recurrence and norms, exact for n < 30;
- W4: the tangent-number C agrees with the Mahler C to 2^{−215};
- W5: the growth table, (e₂ − K²/6 − (K/3)log₂K)/K ∈ [0.89, 1.22].
- th_lemmaA_balanced.py (K = 4, 5, 8, 11): multiplication by x = the monic Jacobi matrix, the balanced bounds, M^S = Σ s_k·C(𝒥, k), margins ≥ 1.
- th_kraw_recurrence.py (K = 3…24): the recurrences and norm valuations are exact.
- No statement changed; only the proofs.

# Real G by modular arithmetic: levels 6, 12, 16 and Γ₁(8) (2026-09-29; level12_search.py, level12_modular.py, gamma1_8_family.py, gamma1_8_kill.py, gamma1_8_amp.py, gamma1_8_amplaw.py). The user: "Might need to attack 3G/6 or 6G/12 more modular arithmetic"
**Reading.** Case E's limit is G/2 = 3G/6 = 6G/12. So attack it at level 6, where CDT proved L(2,χ₋₃) irrational, or at level 12 = lcm(3, 4), where χ₋₃ and χ₋₄ both live.
- Mod 12: S₁+S₅+S₇+S₁₁ = (2/3)ζ(2), S₁+S₅−S₇−S₁₁ = (10/9)G, S₁−S₅+S₇−S₁₁ = (5/4)L(2,χ₋₃), S₁−S₅−S₇+S₁₁ = L(2,χ₁₂) = π²/(6√3), with S_r = Σ_{n≡r (12)} n⁻².
- Level 6 alone cannot carry G: χ₋₄ has conductor 4, so every modular G-family lives on a level divisible by 4.

**The gate, sharpened.** S = ln(16r) − τ (the λ-function bound) is necessary.
- It is meaningful only in a coordinate whose point at infinity is a cusp.
- A Möbius frame t/(1 + at) that makes ∞ a regular point only adds an omitted point (maps must avoid t = −1/a), so a better score there is an artifact. At level 8 such frames scored +0.14 for (1 + 6G)/16 until this was seen.
- A symmetric pair ±s of singular points lowers the bound from 16|s| to 4|s| (square root of λ(q²)).

**Step 3: binomial sums at level 6/12 (level12_search.py, fixed).**
- Setup: 11 templates × twists z = ±2^a3^b (|a| + |b| ≤ 4). Relations were accepted only if verified to ≥ 70 digits with height ≤ 2000.
- The earlier "hits" came from a float bug in the twists plus PSLQ overfitting.
- Only G-families: case E (S = −0.56) and its rescaling Σ C(n,2k)C(2k,k)²16^{−k} (S = −1.95).

**Step 4: modular families (level12_modular.py).** Beukers' construction:
- Hauptmodul t (an eta quotient with a zero at ∞) and a weight-1 Eisenstein form f: θ₃(q^m)² for χ₋₄, a(q^m) for χ₋₃. Then f = Σ u_n tⁿ.
- v from the inhomogeneous recurrence (L v = t), then the limit, τ and the measured r.
- Cusp values x = 1/t, in frames with the pole at a cusp:
  - Γ₀(6): {9, 1} (case C, r = 1);
  - Γ₀(8): {8, 4} (case E, r = 1/4);
  - Γ₀(12): {6, 4, 3, 2} (pole at 1/2), {4, 2, −2, 1} (pole at 1/3), {±3, ±1} (pole at 1/4);
  - Γ₀(16): {4, 2 ± 2i, 2}.
- **Why 4 | N hurts.** τ ↦ τ + ½ normalizes Γ₀(N) exactly when 4 | N. It acts as q ↦ −q, and on x as x ↦ −c − x (c = twice the q² coefficient of t).
  - The cusp values come in mirror pairs.
  - The pole's mirror x = −c is the dominant one, and each other pair has max ≥ |c|/2. So the next singularity is within a factor 2 of the dominant (8 : 4, 6 : 4, 4 : 2√2).
  - Otherwise the frame has ± pairs and the bound drops to 4r.
- Level 12, best frame: r = 1/4, S = −0.55, the same as case E.
- **New periods at level 12.** M₃(Γ₀(12), χ₋₄) has dimension 6: 4 Eisenstein series plus 2 cusp forms.
  - The cusp forms are a non-CM conjugate pair: a₅ = −2 (T₅: (λ−26)⁴(λ+2)²) and a₇ = ±4√−3 (T₇ factor λ² + 48).
  - No level-12 limit matched any basis of G, L(2,χ₋₃), π², π²√3, π, π√3, logs or 1/π powers (140+ digits). These cusp periods are the likely extra constants.
- **Level 16.** One family (θ₃(q)², pole at 1/2) converges to L(η(4τ)⁶, 2) = Γ(1/4)⁴/(64π) (checked; L(η(4τ)⁶, 1) = Γ(1/4)⁴/(32π²)).
  - Its τ is 0.98, not 2: the CM form lives only on primes ≡ 1 mod 4, so its Eichler integral carries half the denominators.
  - Its naive S is +0.81. This is a CM period (transcendental by Chudnovsky), not G. It shows what an arithmetic discount looks like; G's Eisenstein series has none (w_p = p² + χ₋₄(p) ≠ 0 at every p).

**Step 5: breaking the symmetry on Γ₁(8) (gamma1_8_family.py).** τ ↦ τ + ½ does not normalize Γ₁(8).
- Hauptmodul T = (1 − θ₃(q²)/θ₃(q))/2 = q − 3q² + 6q³ − 11q⁴ + …, integral.
- Cusp values x ∈ {4+2√2, 2, 4−2√2, 1}: the dominant-to-next ratio is 3.41 instead of 2.
- **New G-family** (frame T/(1 − 2T), pole at the x = 2 cusp; f = θ₃(q)²): the limit is exactly (1 + 2G)/8, τ = 1.949, r = 1/2, S = +0.141.
  - In that frame the coordinate is T/(1 − 2T) = q − q² − 2q³ + 3q⁴ + ….
  - u = 1, 4, 8, 24, 88, 336, 1352, 5584, 23576, … (integers) and v = 0, 1, 3, 76/9, 280/9, 26792/225, ….
  - Recurrence order 5 (6 terms), degree 2, leading coefficient (n+5)². Characteristic polynomial (x+2)²(x+1)(x² − 4x − 4): the dominant root is 2+2√2, and the next singularity is the old pole cusp at t = −½.
  - This is the best G score of any holonomic family in the program (case E −0.55, Zudilin −1.59).
- Frame 0 converges to L(g, 2), where g = η(τ)²η(2τ)η(4τ)η(8τ)² (weight 3, level 8, χ₋₈, CM by Q(√−2)). τ = 1.03 (CM lacunarity again), S = +1.05; not G.
- Frames with the pole at the x = 1 cusp score ≈ +0.8 naively, but their next singularities are the pair ±1, so the true bound is ≤ ln 4 − τ < 0.

**Step 6: the cusp T = ½ cannot be killed for G (gamma1_8_kill.py, gamma1_8_amp.py, gamma1_8_amplaw.py).** Inhomogeneities h ∈ {T^j, T^j/(1−2T), T^j/(1−T)} in the frame-0 recurrence.
- Limits:
  - Λ(T) = L(g,2);
  - Λ(T²) = (2L(g,2) − G)/4;
  - Λ(T/(1−2T)) = π²/16;
  - Λ(T/(1−T)) = (96L(g,2) − π²)/72;
  - Λ(T²/(1−T)) = (24L(g,2) − π²)/72;
  - Λ(T³/(1−T)) = (18G − π² − 12L(g,2))/72.
- Hence h = T² − T/2 is a pure G-family (limit −G/4).
- **Amplitude law at T = ½.** Let A = lim nF_n/2ⁿ, the log amplitude of F = y(h) − Λ(h)f there (40 digits; PSLQ-exact on 5 of 6 cases).
  - If Λ = c_Lg·L(g,2) + c_G·G + c_π²·π², then A = −6c_Lg·L(g,1) − (π/2)c_G.
  - The ζ(2) part is analytic at T = ½.
- **Consequence.** A family whose limit contains G is singular at T = ½: π and L(g,1) are Q-independent, since L(g,1)/π is a CM period ratio. So r ≤ ½ and S ≤ +0.14 for every such G-family on this curve with this f.
- The G part pays exactly π/2 = 2L(1,χ₋₄) there. That is the s = 1 critical value of the same Eisenstein series whose s = 2 value is G: a modular-symbol statement.
- T³ brings a further constant (not identified), so it is not covered.

**Verdict on 3G/6 and 6G/12.**
- Level 6 cannot hold G.
- Level 12 has case E's geometry (r = 1/4) and adds cusp-form periods.
- The obstruction is the prime 2 again: χ₋₄ forces 4 | N, and 4 | N brings the half-translation that pairs the cusps.
- Breaking the pairing (Γ₁(8)) lifts the best G score from −0.55 to +0.14, and the amplitude law caps it there. CDT's working margin is ≈ 0.77.

| family | limit | τ | r | naive S |
|---|---|---|---|---|
| case E, Γ₀(8) | G/2 | 1.95 | 1/4 | −0.55 |
| Γ₀(12), best frame | unidentified (cusp periods) | 1.95 | 1/4 | −0.55 |
| Γ₁(8), frame T/(1−2T) | (1 + 2G)/8 | 1.95 | 1/2 | +0.14 (capped by the amplitude law) |
| case C, Γ₀(6) (CDT) | L(2,χ₋₃)/2 | 1.95 | 1 | +0.77 |
| Γ₀(16), CM | Γ(1/4)⁴/(64π) | 0.98 | ≈ 0.35 | +0.81 (not G) |
| Γ₁(8), frame 0, CM | L(η(τ)²η(2τ)η(4τ)η(8τ)², 2) | 1.03 | 1/2 | +1.05 (not G) |

**Open, in order.**
1. Other weight-1 forms on Γ₁(8) (θ₃(q)θ₃(q²), θ₃(q²)², combinations), and Γ₁(12), Γ₁(16): does some curve put G's π/2 cusp far away?
2. Prove the amplitude law as a modular-symbol identity: the Eichler integral of the χ₋₄ Eisenstein series between the dominant cusp and T = ½.
3. The CM discount (τ ≈ 1) is the only arithmetic saving seen. G's Eisenstein series has none, so any real gain for G must come from the geometry (r).

# Case E exactly: G/2 is G times the Ramanujan value of ±1, and the obstruction is the quarter turn (2026-09-29; caseE_exact.py). The user: "It's not pure 6 though. It's more like ramanujan summation like before and -1 plus1 like everything else. We can't book keep our way to a better bound than cdt. We have to be analytically correct. Exact thru modular linear algebra like everything else."
**Correction of reading.** "3G/6 or 6G/12" meant G/2 = G × ½ with ½ = 3·⅙ = 6·(1/12), a Ramanujan-summation value. It did not mean the modular level 6 or 12. The previous section's "level 6 cannot hold G" answered a question that was not asked. Its gate scores are bookkeeping and are not a route past CDT.

**Setting (exact).** Take a weight-1 form f and a Hauptmodul t (t(∞) = 0), with the ODE p(t)θ²y + ⋯ = 0 and p(0) = 1.
- Apéry's second solution v (v₀ = 0, v₁ = 1) solves L y = t.
- So y = f·g with D²g = W := (Dt)²/(t·p(t)·f), a weight-3 form, and g = Σ w_n n⁻² qⁿ (D = q d/dq).

**Lemma (derived).** Let the cusp be a/c = γ∞, with γ = (a b; c d) ∈ SL₂(ℤ), and let W vanish there. Put τ = γτ′. Then
    g(γτ′) = L(W, e(a/c·), 2) − 2πi·L(W, e(a/c·), 1)/(c(cτ′ + d)) + g_γ(τ′)/(cτ′ + d),
where L(W, e(x·), s) = Σ w_n e(nx) n⁻ˢ (continued) and g_γ is the Eichler integral of W|γ. So
    F = f·(g − L) = (f|γ)(τ′)·[c·(L(W, e(a/c·), 2) − L)·τ′ + analytic].
- The dominant cusp fixes L: at the cusp 0, L = L(W, 2).
- Every other cusp carries the log amplitude c·[L(W, e(a/c·), 2) − L]; in terms of t it is this divided by 2πi.
- The s = 1 values (the −1/12 rung) cancel out of every amplitude.

**Ramanujan summation evaluates it.** For E_G(τ) = Σ_n (Σ_{d|n} χ₋₄(n/d)d²) qⁿ (L-function L(s,χ₋₄)ζ(s−2)):
    L(E_G, e(x·), 2) = Σ_m χ₋₄(m) m⁻² · Σ_{d≥1} e(mdx),   Σ_{d≥1} e(dy) = e(y)/(1 − e(y)) = −½ + (i/2)cot πy,
with ζ(0) = −½ when y ∈ ℤ.
- The −½ is the same for every twist, so it carries −G/2 at every cusp and cancels between cusps.
- What is left is (i/2)Σ_m χ₋₄(m)cot(πmx)/m².
  - A half turn (x ∈ ½ + ℤ, twist ±1) gives zero, because −1 + 1 − 1 + ⋯ = 1 + 1 + 1 + ⋯ = −½.
  - A quarter turn (twist ±i) gives cot(πm/4) = χ₋₄(m), hence iπ²/16.
- E_G(dτ) sees the cusp a/c through the twist d·a/c.

**Case E, exact.**
- Data: t = η(τ)⁴η(2τ)⁻¹⁰η(4τ)²η(8τ)⁴ (pole at the cusp ½, t(0) = ⅛ dominant, t(¼) = ¼), f = θ₃², p = (1 − 4t)(1 − 8t).
- W = E_G(τ) − 8E_G(2τ) = −E_G(τ + ½), by exact rational linear algebra on q-expansions (coefficients 1, −4, 8, −16, 26, −32, 48, −64, 73, …). Case E's weight-3 form is G's Eisenstein series with q ↦ −q.
- Limit: L(W, 2) = −Σ_m χ₋₄(m)m⁻² · [Σ_d (−1)^{md}] = −Σ_m χ₋₄(m)m⁻² · (−½) = G/2. This is exactly case E's limit: G times the Ramanujan value of the ±1 sum.
- Obstruction at t = ¼ (cusp ¼, c = 4): the E_G(τ) part sees a quarter turn and gives iπ²/16. The −8E_G(2τ) part sees 2·¼ = ½, a half turn, and gives 0.
  - So F = (π/8)·log(1 − 4t) + analytic, i.e. nF_n/4ⁿ → −π/8.
  - Check of the formula (not a scan): −0.38883, −0.39075, −0.39140, −0.39158 at n = 100, 200, 300, 350, against −π/8 = −0.39270, approaching like 0.39/n.
- In words: the quarter turn i adds (i/2)χ₋₄(m) to the ±1 sum, and that turns G's series into the odd part of ζ(2), π²/8.

**Consequence.** A G-component E_G(dτ) is invisible at a cusp exactly when it sees a half turn there (2da ≡ c mod 2c).
- So a G-family can be analytic at a cusp only if every G-component meets that cusp at a half turn, and no other component spoils it.
- This is an exact linear-algebra problem over these Ramanujan values: no PSLQ and no measured gates.
- The Γ₁(8) amplitude law (the G part pays π/2 at T = ½) should follow from the same Lemma. Not yet derived.

# ∞ − 1 as the regularizer: the cusp dictionary, the cokernel, the every-n bound (2026-09-29; regularizer_exact.py). The user: "We want it to blow up in a way. We prove a bound out to infinity minus 1. We just have to do it through linear algebra." Asked to choose, they took: use the blow-up itself; 1/(1 − e(y)) is the infinity, and subtracting 1 regularizes it to −½, which carries G.
**The regularizer at every twist.** Σ_{k≥0} e(ky) = 1/(1 − e(y)) = ½ + (i/2)cot πy for y ∉ ℤ, so Σ_{k≥1} e(ky) = (½ − 1) + (i/2)cot πy.
- The real part of the ∞ is ½ at every twist, and the −1 makes it −½ = ζ(0). All of G sits there.
- The blow-up (i/2)cot πy is purely imaginary and never carries G.

**Cusp dictionary for G's Eisenstein series.** Take E_G(dτ) at the cusp a/c, and put y = d·a/c mod 1, c″ = the denominator of y. From L(E_G(d·), e(a/c·), s) = d⁻ˢ Σ_m χ₋₄(m) m⁻ˢ Li_{s−2}(e(mdy)):
- The ∞ (the pole at s = 3, i.e. the constant term of W at the cusp) is d⁻³χ₋₄(c″)c″⁻³L(3,χ₋₄). It is nonzero exactly when c″ is odd and zero for every even c″. Where it is nonzero, F has a log² blow-up.
- The finite part at s = 2 is d⁻²·(−G/2) at every cusp: the same −½ everywhere.
- The log blow-up is d⁻²(i/2)Σ_m χ₋₄(m)cot(πmy)/m²:
  - zero at a half turn (c″ = 2);
  - iπ²/16 at a quarter turn;
  - π² times an algebraic number otherwise.
- So a G-component is completely invisible at a cusp exactly when it meets it at a half turn.

**Case E and Γ₁(8), read through it.**
- Γ₀(8): cusps ∞ (t = 0), 0 (t = ⅛, dominant), ¼ (t = ¼), ½ (the pole).
  - W = E_G(τ) − 8E_G(2τ). Its ∞'s cancel at the dominant cusp: 1 − 8·(1/8) = 0.
  - At ½, E_G(2τ) is untwisted, so W ≠ 0 there, and ½ is exactly where the frame puts t = ∞.
  - At ¼, E_G(τ) meets a quarter turn: the (π/8)log(1 − 4t) blow-up. E_G(2τ) meets a half turn there and is invisible.
- Γ₁(8), from the θ-asymptotics of T = (1 − θ₃(q²)/θ₃(q))/2:
  - T(0) = (2 − √2)/4 (dominant, untwisted);
  - T(½) = ∞ (the pole, a half turn);
  - T(¼) = ½, the blocking cusp, a quarter turn.

**Exact linear algebra (regularizer_exact.py).**
- **Case E cokernel.** The cokernel of the recurrence operator on t·ℚ[t] is one-dimensional (L(1) = −4t + 32t², L(t) = t − 28t² + 128t³).
  - So every inhomogeneity h gives y_h = κ_h y_t + P_h − P_h(0)f: one limit and one blow-up, in a fixed ratio.
  - Example: h = t² has κ = ⅛, P = 1/32, limit G/16 − 1/32.
- **Case E Casoratian.** u_n v_{n+1} − u_{n+1} v_n = 32ⁿ/(n+1)² exactly (checked for n < 40). Hence, for every n,
    G/2 − v_n/u_n = Σ_{k≥n} 32^k/((k+1)² u_k u_{k+1}),   G = 2Σ_{k≥0} 32^k/((k+1)² u_k u_{k+1}),
  with u_k = Σ_j C(k,j)C(2j,j)C(2k−2j,k−j). This is the bound out to infinity, exact term by term.
- **Γ₁(8).**
  - The recurrence has order 4, with leading polynomial p(T) = (1 − T)(1 − 2T)(1 − 8T + 8T²).
  - The cokernel is 3-dimensional (T, T², T³), with exact reductions T⁴ ≡ T/16 − 5T²/8 + 3T³/2, T⁵ ≡ 5T/48 − 67T²/72 + 121T³/72, and so on.
  - The weight-3 form of the class T, W₀ = (DT)²/(T·p·f), equals the CM cusp form g = η(τ)²η(2τ)η(4τ)η(8τ)² exactly (115 coefficients). Hence Λ(T) = L(g,2).

**What "blow up only at infinity" requires.** Suppose a pure-G weight-3 form W had, at every finite cusp, no ∞ and zero blow-up.
- Then F = f(g − L) blows up only at t = ∞, so F is entire and |F_n| falls below every geometric rate.
- F is not a polynomial, since its Eisenstein projection is nonzero.
- So q·d_n²·F_n ∈ ℤ ∖ {0} would force G irrational.
- The whole question is therefore the linear algebra of the ∞'s and the twists. For pure-G Eisenstein combinations:
  - **Γ₀(8).** αE_G(τ) + βE_G(2τ); the dominant ∞ forces α : β = 1 : −8, and that form meets ¼ at a quarter turn. It fails by one condition.
  - **Γ₀(12)** (the "6G/12" check), with d ∈ {1, 3}. There are three conditions:
    - the ∞ at ⅓ forces α₃ = α₁ (χ₋₄(3) = −1 against 3⁻³);
    - the quarter turns at ¼ cancel only if α₃ = 9α₁ (cot(3πm/4) = −χ₋₄(m));
    - the dominant ∞ forces α₁ + α₃/27 = 0.
    - Any two of them kill everything, and at most one of 0, ⅓, ¼ can be the pole. It fails.
  - **Γ₀(16).**
    - Quarter turns at both ¼ and ¾ force α₁ = 0 (only one of them can be the pole).
    - The ∞'s of E_G(4τ) at ¼ and ¾ force α₄ = 0.
    - E_G(2τ) has an ∞ at ½ and a quarter turn at ⅛, so α₂ = 0.
    - It fails.
- In each case the ∞-cancellation at the dominant cusp (weights d⁻³) collides with the cancellations needed at the other cusps.
- Next: the full classification, allowing E_ζ (1, χ₋₄) components and cusp forms that cancel among themselves, on the Γ₀(4k) and Γ₁(N) curves. It either finds a family that blows up only at infinity, or it explains exactly why G always keeps one ∞ or one quarter turn.

# The widened linear algebra: the layer theorem, the one survivor, and why it is "infinity minus 1" (2026-09-29). The user: "Ok"
**A missing ingredient: the weight-1 form f.** F = f·(g − Λ) also inherits f's behaviour at each finite cusp.
- At a regular cusp where f ≠ 0, the dictionary above is complete.
- Where f vanishes (a square-root zero at an irregular cusp), write F(γτ′) = (f|γ)·[c(L₂ − Λ)τ′ − 2πiL₁/c + …]. The s = 1 term −2πiL₁/c, which only shifted the analytic part at a regular cusp, now multiplies (t − t_c)^{1/2}.
- So at such a cusp the s = 1 value (the −1/12 rung) is the blow-up.
- Realizability: every weight-3 E that vanishes at ∞ and at the finite cusps comes from a polynomial inhomogeneity h = t·E/W₀, provided f has no zeros at finite cusps.

**Decoupling (exact, given the standard independence of π² from the L(2, odd ψ) Clausen values).**
- At the representative x = a/c, E_G's ∞ is real and E_ζ's is imaginary: μ_{χ₋₄}(r/4) = (i/2)χ₋₄(r), and the residue is 2^{−3(v−2)}(i/2)χ₋₄(r)L(3,χ₋₄) for E_ζ(2^f τ) at a cusp r/2^v with v ≥ f + 2.
- E_G's blow-ups are i·π²·(algebraic), while E_ζ's values are π²·(rational) plus i·(Clausen values).
- So the G-type part must satisfy the conditions by itself.
- On Γ₁(N), the other odd quadratic characters decouple by square class: L(3,χ₋₈) = 3√2π³/128 and L(3,χ₋₃) ∈ √3·ℚπ³, against L(3,χ₋₄) = π³/32; the cot sums follow the same pattern. With rational coefficients, χ₋₄ stands alone again, on more cusps.

**Layer theorem (Γ₀(N), N = 2^a M, M odd; G-components E_G(2^f d′τ), 0 ≤ f ≤ a − 2, d′ | M).** Sort the cusps by layer v = v₂(c), with c = 2^v c′.
- **∞'s.** In layer v the components with f ≥ v carry an ∞. The condition matrix T_{c′,d′} = d′⁻³χ₋₄(c′_{d′})c′_{d′}⁻³ (c′_{d′} = c′/gcd(c′,d′)) factors over the primes of M. After scaling column j by p^{3j}, each prime block has 1 on and above the diagonal and ρ^{i−j} below it (ρ = χ₋₄(p)p⁻³), so its determinant is (1 − ρ)^e ≠ 0.
  - Hence B^{(v)} := (Σ_{f≥v} α_{f,d′} 8^{−f})_{d′} = 0 in every layer v ≤ a − 2 that holds no pole.
  - The pole can remove a row only where a cusp stands alone, i.e. layers v₀ ∈ {0, 1}.
  - That forces W = Σ_{d′} β*_{d′}(8^{v₀}E_G(2^{v₀}d′τ) − 8^{v₀−1}E_G(2^{v₀−1}d′τ)), with β* = T⁻¹e_{pole}. For v₀ = 1, d′ = 1 this is exactly case E's form.
- **Blow-up.** In both cases (for v₀ = 0 when N ≠ 4) the cusp ¼ is finite and meets the lower component at a quarter turn. Its condition is Σ_{d′} β*_{d′}χ₋₄(d′)d′⁻² = (R·T⁻¹)_{pole} = 0.
  - That fails for every pole choice: M = 1 gives 1 ≠ 0; M = p gives 1 − p⁻² or 1 ∓ p; general M by multiplicativity.
- **The one survivor** is N = 4, pole at cusp 0, W = E_G. Its only finite cusp is ½, a half turn: no ∞ and no s = 2 blow-up.
  - But ½ is irregular, and the only weight-1 form, θ₃², vanishes there: θ₃(τ + ½) = θ₄(τ), and Σ_{n∈ℤ}(−1)ⁿ = 1 + 2·(−½) = 0.
  - So ½ pays through the s = 1 value: L(1,χ₋₄)·Σ_d d(−1)^d = (π/4)(−¼), the alternating cousin of −1/12 (η(−1) = (1 − 2²)ζ(−1) = ¼).
  - The same ±1 sum that frees G at the half turn (Σ(−1)^d = −½) makes θ₃ vanish there.
- **Check of the survivor (numbers only to test the derived statement).** h = t + 16t² gives W = (1 + 16t)W₀ = E_G exactly (80 coefficients). y_n/u_n → −G/2, as predicted (L(E_G, e(½·), 2) = −G/2), but only logarithmically: (y_n/u_n + G/2)·log n = 1.22, 1.39, 1.49, 1.57, 1.60 at n = 10², 10³, 10⁴, 10⁵, 4·10⁵. That is the square-root branch at t = −1/16.

**Conclusion: "infinity minus 1" is exactly the shape.**
- Every pure-G Beukers family (Γ₀(N); Γ₁(N) with quadratic characters) blows up at at least one finite cusp. Blowing up only at infinity is impossible, and one paying cusp class is the best possible. Case E already achieves it (the cusp ¼).
- The bound then runs out to that one cusp, and its distance decides everything:
  - The x = 1/t values of the cusps are algebraic integers (8, 4 / 6, 4, 3, 2 / 4 ± 2√2, 2, 1 / 2 ± 2i on the curves computed).
  - Galois-conjugate cusps share the twist type (a quarter turn stays a quarter turn), so a paying cusp's whole orbit pays.
  - The norm is ≥ 1, so the nearest paying cusp has |t| ≤ 1 in the integral frame.
- Apéry's direct criterion needs |t| > e^τ ≈ e² ≈ 7.4. The deficit is at least e² per step, and it is structural.
- CDT's λ-uniformization is worth at most a factor 16. So the whole question is R ≥ e²/16 ≈ 0.46 for the single paying orbit: case E has R = ¼ (fails), Γ₁(8) has R = ½ (the naive margin +0.1).
- Not yet covered: Atkin–Lehner quotients (their elliptic fixed points bring square-root branches) and non-quadratic character orbits on Γ₁(N).

# Uniformized at infinity: the exact radius is 2π²·Σ1/c², and G sits at width 8 (2026-09-29; uniformize_exact.py, uniformize_radius.py). The user: "Well pi is baked in at infinity. For the proof. Pi is extracted from modular arithmetic just like ramanujan summation." Asked how to use it, they chose: uniformize at infinity.
**The surface the forms live on.**
- F = f·(g − Λ) is holomorphic on H.
- By the Lemma it is invariant under the parabolic at every harmless representative x (no ∞ there, and L(W, e(x·), 2) = Λ), and under the translation T.
- So F lives on X′ = H/Γ′ with the harmless cusps filled in, where Γ′ = ⟨T, harmless parabolics⟩.
- X′ is simply connected, because Γ′ is generated by the parabolics being filled.

**The exact radius.**
- The Green's function of X′ with its pole at the filled cusp ∞ is 2π·E_{Γ′}(τ, 1) = 2πΣ_{Γ′_∞\Γ′} Im(γτ). It converges at s = 1 because Γ′ has infinite covolume.
- Its constant term at ∞ gives the conformal radius in the coordinate t = q + O(q²):
    log R* = 2π² · Σ 1/c(γ)²   (over the double cosets Γ′_∞\Γ′/Γ′_∞, γ ≠ 1),
  because the constant term of Σ y/|cτ + d|² over d mod c is π/c².
- π is "baked in at infinity": q = e^{2πiτ}, and each class contributes π/c².
- The sum is pure modular arithmetic. For ⟨T, (1 0; w 1)⟩ its leading part is Σ_{m≠0} 1/(wm)² = 2ζ(2)/w² = π²/(3w²), which is Ramanujan's B₂ = ⅙ once more.
- The gate is log R* > τ.

**The Γ₁(8) family is case E.**
- The pure-G form of the class h = T² − T/2 decomposes exactly as W = −½E_G(τ) + 4E_G(2τ), which is −½ times case E's form. f = θ₃² in both.
- So it is the same function on H, only expanded in T instead of t.
- The radius is intrinsic (dt/dq = dT/dq = 1 at the cusp). The "+0.14 on Γ₁(8)" was an artifact of the coordinate.

**Harmless representatives are an orbit.**
- For case E's F, x is harmless exactly when the cotangent sums balance: K(x) = 2K(2x), with K(y) = Σ_m χ₋₄(m)cot(πmy)/m².
- Checked exactly: 0, 1/7, 1/9, 1/15, 2/15, 1/17, 2/17, 1/23, 3/23, … are harmless; ¼, ⅓, ⅜, 1/5, 3/7, 1/11 pay.
- Every harmless odd-denominator x < 1 with c < 40 lies in the orbit of 0 under ⟨T, (1 0; 8 1)⟩, and every orbit point is harmless. So the forms live exactly on H/⟨T, (1 0; 8 1)⟩.
- Case C's dominant singularity t = 1/9 is the width-6 cusp 0 of Γ₀(6) (Fricke value ∏δ^{−r_δ/2} = 1/9). So its group is ⟨T, (1 0; 6 1)⟩.

**The numbers** (word sums cut at |c| ≤ X = 25600, geometric tail extrapolated; the partial sums are lower bounds).

| width invariant w | example | S = Σ1/c² | log R* | λ-bound log(16·r_pay) | vs τ |
|---|---|---|---|---|---|
| 5 | Apéry ζ(2) (Γ₁(5), cusp 0) | ≈ 0.260 (≥ 0.2448) | ≈ 5.13 | log 177 = 5.18 | passes directly |
| 6 | case C, L(2,χ₋₃) (CDT) | 0.1357 (≥ 0.1338) | 2.68 (≥ 2.64) | log 16 = 2.77 | +0.74 |
| 7 | — | 0.0877 (≥ 0.0871) | 1.73 | — | −0.22 |
| 8 | case E and Γ₁(8), G | 0.0624 (≥ 0.0622) | 1.232 | log 4 = 1.386 | −0.72 |

In every case R* sits just under the λ-bound. The threshold lies between widths 6 and 7.

**Why G is at 8: the prime 2 one last time.**
- χ₋₄ forces 4 | N.
- The only level-4 configuration (Γ₀(4), W = E_G, pole at 0) has the invariant-4 pair (∞, ½). There ⟨T, p_½⟩ = Γ₀(4) has finite covolume, so the forms would blow up only at infinity.
- But θ₃² vanishes at ½, so F has a square root there and is invariant only under p_½². That doubles the width: invariant 8, the same as case E.
- The same ±1 sum that frees G at the half turn (Σ(−1)^d = −½) makes θ₃ vanish there (1 + 2·(−½) = 0), and that turns width 4 into width 8.
- Result: G's forms reach log R* ≈ 1.23 against denominators e^{1.95n}. L(2,χ₋₃) sits at width 6 (2.68) and ζ(2) at 5.

**Open.**
- A G-form whose harmless group has a second fat parabolic. S must rise from 0.0624 to about 0.099: one more invariant-8 parabolic adds 0.051 at first order.
- Γ₀(N) larger than 8 makes the widths larger. The layer-theorem candidates on N ≤ 48 should be scanned exactly for their harmless groups.

# The half measure and the missing factor 2 (2026-09-29). The user: "G will be found at a half measure. 1/2G. What truly breaks this is a prime factor over infinity. This gives us the offset needed. Anyways. Continue … Maybe something like 113/120 or something"
**The constant, sharpened.**
- S(8) = 0.062046, 0.062226, 0.062320 at X = 6400, 25600, 102400. The tail ratio is 0.529 per ×4, so S(8) = 0.062427.
- Hence log R* = 1.2323 and R* = 3.429. This is close to π²/8 = 1.2337 (S = 1/16), but not equal: the gap 0.00007 in S is outside the tail estimate.

**The deficit is one factor of the prime 2.**
- τ − log R* = 1.948 − 1.232 = 0.716, i.e. R* is 1/2.05 of e^τ. G is found at a half measure, literally.
- log 2 = 0.693 accounts for all but 0.023, which is 1.2% of τ.
- Consequence: if the radius were doubled (log R* → 1.925), a denominator saving of 113/120 (τ → 1.834) would pass with margin 0.09. A saving of 1.2% would already suffice.

**Where the factor 2 is lost (exact).**
- The two-parabolic groups are Hecke groups in disguise. ⟨T, (1 0; w 1)⟩ is normalized by the Fricke involution (0, −1/√w; √w, 0), and the two together are conjugate (τ ↦ √w·τ) to the Hecke group H(√w).
  - Apéry ζ(2) ↔ H(√5); case C, L(2,χ₋₃) ↔ H(√6); G ↔ H(√8).
  - The boundary case H(2) is θ₃'s theta group, of finite covolume (radius ∞).
- The only level-4 G-configuration is Γ₀(4), W = E_G, pole at 0. Its harmless pair (∞, ½) has invariant 4, which would be H(2).
- By the valence formula, the only weight-1 form with character χ₋₄ on Γ₀(4) (θ₃², total zero order 6/12 = ½) must put a zero of order ½ at the one irregular cusp, ½. The forms therefore carry a square root there and are invariant only under p_½², which doubles the invariant 4 → 8: H(2) → H(√8).
- On Γ₀(8) the cusp ½ is regular, but its width is 2, so the invariant is 8 again. Coordinate changes (Γ₁(8), or the double cover σ = (√(1 + 16t) − 1)/8, where F becomes analytic at t = −1/16) leave R* unchanged, since R* is intrinsic.
- Larger levels only raise the invariants: every parabolic in Γ₀(N) has lower-left entry ≡ 0 (mod N), so each harmless pair has invariant ≥ N.
- So among congruence Beukers G-forms, case E's H(√8) is optimal, and the missing factor 2 is exactly the half-order zero of θ₃² at the half-turn cusp.

**What would break it.** A G-form invariant under p_½ itself: a weight-1 χ₋₄ form with no half-order zero at a half-turn cusp of width 1, or the half-translation τ ↦ τ + ½ in the harmless group, which halves the invariant 8 → 4.
- Checked: symmetrizing F(τ) + F(τ + ½) brings back a log² singularity at t = ⅛. W = E_G(2τ) (the only τ+½-invariant G-form) is the Γ₀(4) survivor in the variable 2τ, with the same square root. Both are circular.
- So the factor 2 is not available inside congruence Beukers families. It needs a new ingredient: a non-congruence structure, an adelic contribution (a place where the forms converge beyond radius 1), or a genuine arithmetic saving of 37% (1.2% once the factor 2 is found).

## Rational at infinity, irrational at every prime? (2026-09-29; the user: "Could it be rational at infinity. Irrational for every prime?")
**Logically, yes.**
- The real G = L(2,χ₋₄) and the p-adic Catalan constants (values of Kubota–Leopoldt p-adic L-functions at s = 2) are different numbers, not one number seen in different completions.
- Both are built from the same Ramanujan values β(−2k) = E_{2k}/2. At infinity they are continued through the functional equation, which brings in π: G = (π/2)·β′(−1). At p they are interpolated p-adically, with Euler factors and no π.
- No theorem links the rationality of the one to the other.
- Status:
  - The 2-adic value is proven irrational (Calegari 2005, IMRN no. 20; our C). Beukers 2008 treats such p-adic values by Stieltjes continued fractions.
  - A few small primes and characters are proven; general p is open. The p-adic methods also stop at a gate: overconvergence radius against denominators, which fails for large p.
  - The real G is open.

**Consequence for the adelic "prime factor".** The product formula couples the places only through one and the same rational.
- If G = 2r/q were rational, the integers N_n = q·d_n²·(v_n − r·u_n) would be small at ∞.
- At p = 2 their size is set by |r − L₂|₂, where L₂ is the 2-adic limit of v_n/u_n. That is a fixed nonzero distance if L₂ is irrational, so there is no extra 2-adic decay.
- "Rational at infinity, irrational at every prime" is therefore exactly the scenario that a p-adic factor cannot exclude.
- The missing factor 2 has to come from structure the forms carry at every place independently of the value: their true denominators (the arithmetic side, where 113/120 would live), or a larger archimedean domain.

# Scholze and the Langlands groups: G at every place in one object (2026-09-29; the user: "Scholze and langlands groups")
**What G is, at every place.**
- G = D(i), the Bloch–Wigner dilogarithm at i. That is the Borel regulator of the class [i] generating K₃(ℚ(i)) ⊗ ℚ, which has rank 1.
- Geometry:
  - G = ¼·vol(Whitehead link complement). The regular ideal octahedron is 4 tetrahedra of shape i.
  - G = 3·covol(PSL₂(ℤ[i])), by Humbert: |d|^{3/2}ζ_{ℚ(i)}(2)/(4π²) = G/3.
  - The Whitehead link group has index 12 in the Picard group: 12·G/3 = 4G.
- p-adically: Coleman's p-adic dilogarithm D_p(i) is the p-adic regulator of the same class (Besser–de Jeu).
  - Garoufalidis–Scholze–Wheeler–Zagier (arXiv 2412.04241) Lemma 3.1: D_p(ζ) ∈ p²ℤ_p[ζ] for roots of unity.
  - Their Prop. 3.2: p⁻²D_p(ζ^p) ≡ (ζ − 1)^{−p}·Σ_{k<p} ζ^k/k² (mod p), Kontsevich's finite polylogarithm. So the p-adic Catalan values mod p are finite, π-free sums over residues.

**The one object (GSWZ, "The Habiro ring of a number field", v2 Aug 2025).**
- Modules over the Habiro ring of K, graded by K₃(K). Their elements are power series at every root of unity that are integral and glue p-adically after a Frobenius twist, after dividing at each prime by a series depending only on the Bloch element.
- The main theorems involve the Borel, p-adic and étale regulators together.
- The perturbative Chern–Simons series of hyperbolic 3-manifolds are elements. For the Whitehead link this packages 4G at infinity and D_p(i) at each p in one integral object.
- It is built, like our cusp dictionary, from data at roots of unity.

**The Langlands side.**
- The Picard/Bianchi group PSL₂(ℤ[i]) (GL₂ over ℚ(i)). Scholze (Annals 2015) attaches Galois representations to torsion classes in the cohomology of such groups.
- Bergeron–Venkatesh (arXiv 1004.1083) conjecture log|H₁(Γ_n)_tors| ~ vol(Γ_n)/(6π) for trivial coefficients. For congruence towers over ℚ(i) that would make the logarithms of the torsion primes add up to index·G/(18π).
  - Proved: exponential growth for strongly acyclic coefficient systems (Bergeron–Venkatesh), and for Bianchi congruence towers with symmetric-power coefficients (Pfaff, arXiv 1302.3079).
  - CORRECTION (2026-09-29, later): this line first called the trivial-coefficient law proved. It is a conjecture.
- The exact finite form is Cheeger–Müller: Reidemeister torsion (primes) = analytic torsion (the archimedean place). This is the most literal "prime factor over infinity" in existing mathematics, with G on the archimedean side.

**What none of this does.**
- None of it transfers (ir)rationality between places. Borel's theorem (the regulators vanish only on torsion; GSWZ §1.7) is nonvanishing, not irrationality.
- "Rational at ∞, irrational at every p" stays in the territory of the period conjecture.

**Why it matters here.** The missing factor 2 must be structural at every place. A Habiro-module element for [i] is exactly such a structure: integral at every root of unity and glued across primes by Frobenius. An adelic holonomy argument on that object, rather than on one power series at one cusp, is where the factor could live.

**Next.**
- Write the Whitehead-link / [i] element in Nahm form (shapes among i, 1 + i, (1 + i)/2; the Neumann–Zagier data in GSWZ's (A, z) form).
- Check its q → 1 growth, 4G/(2π), and its integrality at roots of unity.
- Compare it with our case-E forms.

## Langlands dual groups and automorphic forms (the user, same day: "Langland dual groups and automotphic forms")
**The adelic reading of the user's hints.**
- An automorphic representation is π = ⊗_v π_v over all places, and L(s,π) = L(s,π_∞)·∏_p L(s,π_p).
- "π baked in at infinity" is the archimedean factor, a Γ-factor carrying π.
- The "prime factors" are the Euler factors, i.e. the Satake parameters in the dual group at each p.

**G's automorphic homes.**
- GL₁: χ₋₄ (dual group ℂ^×).
- GL₂ over ℚ, by automorphic induction from ℚ(i): θ₃² = E₁(1, χ₋₄), with L-function ζ_{ℚ(i)}(s) = ζ(s)L(s,χ₋₄).
- The weight-3 E_G is the Borel Eisenstein series of (χ₋₄, |·|²) (dual: the diagonal torus of GL₂(ℂ)).
- GL₂ over ℚ(i): the Picard group, of covolume G/3.

**Our cusp dictionary is Langlands' constant-term theory.**
- The ∞ at a cusp is the Langlands–Shahidi constant term, carrying L(3,χ₋₄) = π³/32. That is a critical value, rational times π³: "rational at infinity".
- The finite part carrying G is L(2,χ₋₄), a non-critical value.
- The quarter-turn blow-up is π²·(rational), from critical data again.

**Why G is hard, in this language.**
- s = 2 is non-critical for the odd character χ₋₄ (the Γ-factor has the wrong parity). The functional equation pairs L(2,χ₋₄) with the trivial zero at s = −1, so G is the derivative there: β′(−1) = 2G/π.
- That is Beilinson's regulator (the K₃ class [i]). Borel gives nonvanishing only.
- Criticality is intrinsic to the χ₋₄-piece. So no functoriality (Rankin–Selberg, adjoint, base change) can make G a critical value, which would be π^k × algebraic.
- Every automorphic L-function containing χ₋₄ sees G at a non-critical point.

**What this means for the proof.**
- The critical data (the ∞'s, the π² blow-ups) are exactly the parts we can control rationally.
- G is the one non-critical piece, located at the half measure between the critical points s = 1 and 3.
- Any factor-2 gain must come from structure attached to the regulator itself: the K₃ class [i], which is what the Habiro-module element packages at every place.

# The harmless group is never a lattice: why G's width is 8, and the sign of ∞ at the prime 2 (2026-09-29; lattice_principle.py, lattice_principle_out.txt). The user: "Langland dual groups and automotphic forms."
**The automorphic statement.** Let Γ be the group of the forms (f and W modular on it). Take χ(d) = 1 on Γ, e.g. Γ₁(N).
- For γ ∈ Γ, F∘γ = F + f·p_γ, with p_γ a polynomial of degree ≤ 1. The reason: y∘γ − y solves the homogeneous equation, whose solutions on H are f and τf.
- p is a 1-cocycle with values in V₁. Its class is the Eichler–Shimura class of W; Λ only adds a coboundary.
- The harmless group Γ′ is exactly the stabilizer {γ : p_γ = 0}, which is a subgroup. On the parabolic at x, p vanishes iff x is harmless. So the cusp dictionary is this cocycle evaluated on parabolics.
- **Lattice principle: Γ′ always has infinite index in Γ (it is thin).**
  - Proof 1 (Eichler–Shimura). If p vanished on a finite-index Γ″, the class of W would be a coboundary on Γ″. But W ↦ [p_W] is injective on M₃(Γ″).
  - Proof 2 (Liouville). F would be Γ″-invariant, and invariance under each cusp's parabolic kills every τ′-growing term in the Lemma, so F is bounded at every cusp. Then F is a bounded holomorphic function on a compact curve, hence constant. But F ≡ c is impossible:
    - c = 0 gives g ≡ Λ, but W ≠ 0.
    - c ≠ 0 gives W = c·D²(1/f), and by Bol's identity that has a pole wherever f vanishes. A weight-1 form always vanishes somewhere: its total order is [PSL₂(ℤ) : Γ̄]/12 > 0.
- This is "∞ − 1" at the level of groups. On a genus-0 curve without elliptic points, Γ/±1 is free on the parabolics of all cusps but one. So if every finite cusp were harmless, the product relation would force p = 0 at the pole as well.
  - The pole (the ∞) therefore needs at least one paying finite cusp (the −1).
  - The earlier TARGET, a W with no ∞ and no blow-up at any finite cusp, was impossible for this reason, for every L-value and not only for G.

**Two-parabolic consequence (exact; part 1 of the script).** Up to a real translation, a harmless pair (∞, x) of invariant w generates ⟨T, L_w⟩ with L_w = (1 0; w 1).
- T·L_w⁻¹ = (1−w 1; −w 1) has trace 2 − w. It is elliptic for w ≤ 3, −(parabolic) at w = 4, and hyperbolic for w ≥ 5.
- With the Fricke involution this is Hecke's H(√w), a lattice exactly when w ≤ 4.
  - At w = 4: Sanov's free group ⟨(1 2; 0 1), (1 0; 2 1)⟩, conjugated by diag(2^{−1/2}, 2^{1/2}), gives ⟨T, L₄⟩ = Γ₀(4)/±1.
- So w ≥ 5 always.
- On Γ₀(N) and Γ₁(N) every invariant is a multiple of N (it is lcm(N, c²) at a/c), and the conductor divides N. So w_min is the least multiple of the conductor that is ≥ 5.

| L-value | conductor | w_min | w_min/cond | log R*(w_min) | vs τ ≈ 1.95 |
|---|---|---|---|---|---|
| ζ(2), Apéry (Γ₁(5)) | 5 | 5 | 1 | 5.13 | passes (attained) |
| L(2,χ₋₃), case C (CDT) | 3 | 6 | 2 | 2.68 | passes (attained) |
| G = L(2,χ₋₄), case E | 4 | 8 | 2 | 1.232 | fails (attained) |
| L(2,χ₋₇) | 7 | 7 | 1 | ≤ 1.73 | fails (bound only) |
| L(2,χ₋₈) | 8 | 8 | 1 | ≤ 1.232 | fails (bound only) |

- **The factor 2 is 8/4.** Conductor 4 sits exactly on Hecke's threshold λ = 2 (trace 2, the theta group). So the first allowed width is 2·4.
- χ₋₃ is doubled too (3 → 6), but 6 still passes.
- With two parabolics, χ₋₃ is the only odd quadratic character whose L(2,χ) reaches the gate.

**What breaks at the threshold: the valence formula.** The lattice principle says something must break at w ≤ 4; this is what does.
- **Γ₀(3), weight 1, χ₋₃.** The total zero order is 4/12 = ⅓.
  - The product (−2 1; −3 1) fixes e = ½ + √−3/6. Its automorphy factor there is j = −3e + 1 = −½ − ½√−3, a primitive cube root of 1 (exact in ℚ(√−3)).
  - So f(e) = j·f(e) forces the zero. In weight 3 it would not be forced, because j³ = 1.
  - Check: |a(e)| = 1.05·10⁻⁵¹.
- **Γ₀(4), weight 1, χ₋₄.** The total is ½.
  - σ⁻¹(T·L₄⁻¹)σ = −T⁻¹ exactly, with σ = (1 0; 2 1). So f|σ picks up (−1)^k and has exponents in ½ + ℤ.
  - Ligozat: θ has orders 0, ¼, 0 at 0, ½, ∞.
  - Check: θ(½ + iy)·(2y)^{1/2}/(2|q′|^{1/4}) = 1.0 to 12 digits at y = 0.05, 0.02, 0.01.
- For χ₋₃ the zero comes from a rotation of order 3 at a CM point; for χ₋₄ it comes from a sign at a cusp.

**The sign, adelically.**
- At ½ the sign is (−1)^k = χ₋₄(−1). The same element read as (3 −1; 4 −1) has d = −1.
- As an idele-class character, χ₋₄ has χ_∞ = sgn and is unramified at odd p. The product formula on −1 then forces χ₂(−1) = χ_∞(−1) = −1.
- So the forced half-order zero at the 2-adic cusp ½ is the sign of the real place, carried by the prime 2. That is the exact form of "a prime factor over infinity".
- The same sign makes s = 2 non-critical: odd χ has the Γ-factor Γ_ℝ(s + 1), so G sits between the critical points 1 and 3.
- It also makes the weight-1 parameter 1 ⊕ χ₋₄ odd, so no twist puts it in SL₂(ℂ): squares of characters are even.
- So the lost factor 2 and G's non-criticality are one sign, seen at 2 and at ∞.

**Consequences for the search.**
- Congruence Beukers G-families (Γ₀(N), Γ₁(N)) cannot reach w ∈ {5, 6, 7}. With the layer theorem (no second harmless class), case E's H(√8) is final there.
- Leaving congruence costs denominators: Calegari–Dimitrov–Tang (arXiv 2109.09040, JAMS 2025) prove that bounded denominators force congruence.
- Over ℚ(i) the sign is gone: there is no real place, and χ₋₄ becomes trivial by base change. There G is a covolume (covol PSL₂(ℤ[i]) = G/3).
  - But GL₂ over ℚ(i) has no holomorphic forms (Bianchi forms live on H³), so the Beukers shape does not transfer as it is.
- What is left for the factor 2:
  - an arithmetic saving of 37%;
  - a construction outside the Beukers shape (the Bianchi/Habiro side);
  - or a proof that none exists.

# The ℚ(i) side, pushed exactly: the Whitehead link in the same two-parabolic family, and what it can carry (2026-09-29; qi_whitehead.py, qi_whitehead_out.txt). The user: "I've worked this angle a lot. But let's push it anyways."
**A. Volumes are derivatives at the trivial zero (exact).**
- Humbert gives covol PSL₂(O_K) = |D|^{3/2}ζ_K(2)/(4π²).
- For odd χ_D the functional equation gives L′(−1,χ_D) = |D|^{3/2}L(2,χ_D)/(4π). So covol = −2π·ζ(−1)·L′(−1,χ_D).
- The figure-eight group and the Whitehead link group both have index 12 = −1/ζ(−1) in their Bianchi groups. Hence vol/(2π) = L′(−1,χ_D):
  - 4₁: 0.323065947219451 = L′(−1,χ₋₃);
  - W: 0.583121808061638 = β′(−1) = 2G/π (15 digits).
- So Ramanujan's −1/12 is the index, and the Kashaev growth rate is the derivative of β at its trivial zero s = −1.

**B. Riley: the same two-parabolic family.**
- b(5,3) = 4₁ is ⟨T, (1 0; u 1)⟩ with u² − u + 1 = 0.
- b(8,3) = W is the same group shape with u² + 2u + 2 = 0, i.e. u = −1 ± i.
- (Riley's word needs q odd; b(5,2) gives nothing in this convention.)
- Our harmless groups are the real points u = w of this family: 5 (ζ(2)), 6 (χ₋₃), 8 (G). The constants' 3D homes are complex points of it.

**C. The forbidden lattice sits inside the Whitehead group (exact).**
- Take P₁ = y and P₂ = g·y·g⁻¹ with g = xyx. Then tr(P₁P₂) = −2.
- Conjugating over ℚ(i) (send 0 ↦ ∞ and g(0) ↦ 0, then scale by λ = (1+i)/2, λ² = i/2) gives (1 −2i; 0 1) and (1 0; −2i 1), with invariant −4.
- That is ⟨T, L₄⁻¹⟩ = Γ₀(4). So the thrice-punctured sphere in the Whitehead link complement is exactly the w = 4 lattice that the Eichler–Shimura principle forbids for a holomorphic F.
- Case E's harmless group ⟨T, L₈⟩ is ⟨P₁, P₂²⟩ after that conjugation, so it is a thin subgroup of the Whitehead group.

**D. q → 1: the Kashaev invariant.**
- The formula is Murakami–Murakami–Okamoto–Takata–Yokota's (Exp. Math. 11, 2002); it becomes O(N²) because the i- and j-sums factor for fixed k.
- Their Table 1 is reproduced to all printed digits at N = 40 and 50 (the printed prefactor (−1)^{N−1} dropped).
- A 13-point Richardson fit over N = 160…400 gives J_N = C·N^{3/2}·e^{N(4G + iπ²/4)/(2π)}·(1 + a₁h + a₂h² + a₃h³ + ⋯), h = 2πi/N, with:
  - C² = (−1 + i)/8 (C⁴ = −i/32 to 10⁻²²);
  - a₁ = (39 − 14i)/96;
  - a₂ = (965 − 1452i)/(2·96²);
  - a₃ = (−26505 − 325922i)/(30·96³) (10 digits).
- The denominators are 96^k·(1, 2, 30): the pattern of Garoufalidis–Zagier's 4₁ series with 72√−3 replaced by 96. Only the primes 2, 3, 5 occur.

**E. Integrality at roots of unity.**
- The norms of J_N(ζ_N) over ℚ are rational integers for N = 3…12. Examples: 981 = 3²·109; 3200 = 2⁷·5²; 196075625 = 5⁴·313721; 5461721005027301 = 7⁶·71·653857219.
- Their p-parts give N | J_N(ζ_N) every time: p^{p−1} ∥ Norm for N = p = 3, 5, 7; 3¹² for N = 9; extra powers of 2.

**Verdict: what this angle can carry.**
- Everything is integral and exact, and G enters only as β′(−1) = 2G/π: as a growth rate e^{N·β′(−1)}, a volume, or the exponent of the completion e^{V/h}.
- The natural series are resurgent (factorially growing, as for 4₁). They are not G-functions, so no arithmetic holonomy bound applies to them.
- They carry G/π, never a linear form in 1 and G.
- The one exact contact with our construction is geometric. The w = 4 lattice Γ₀(4), forbidden in 2D, is the totally geodesic thrice-punctured sphere of G's 3D home. G's volume is its thickening, and its holomorphic shadow is exactly the lattice we cannot use.
- This matches the user's own experience of the angle.
- Possible continuations (none aimed at a proof):
  - the series at the roots of unity m = 2, 4, 8, and GSWZ's Frobenius gluing with D_p(i) (the Whitehead link is not among their examples);
  - the Stokes constants.

# The arithmetic side, exactly: a pure-G family with ¾ of the denominators, and why ¾ is the floor (2026-09-29; arithmetic_saving.py, arithmetic_saving_out.txt). The user: "2" (the arithmetic lever, the 37% saving)
**The reduction (exact).**
- Take any integral coordinate T = q + O(q²) with integral inverse, and f integral with f(0) = 1. Then the cumulative denominators of F = f(g − Λ) are exactly those of the Eichler integral g = Σc(m)m⁻²qᵐ.
- So W alone sets the arithmetic. The frame and f do not matter.
- A prime p > √n first enters at the least multiple kp with v_p(c(kp)) < 2.
- For Eisenstein series multiplicative at p, c(kp) ≡ Σᵢ coefᵢ·Aᵢ(p)·Aᵢ(k) mod p, and each Aᵢ(p) is a unit:
  - E_G: p² + χ₋₄(p) ≡ χ₋₄(p);
  - E^{1,χ₋₄}: 1 + χ₋₄(p)p² ≡ 1;
  - E^{1,χ₋₈}: ≡ 1.
- So the p-digit at the k-th multiple is a class function R_k(χ₋₄(p)).

**The admissible space (Γ₁(8), harmless group ⟨T, L₈⟩).**
- Weight 3 on Γ₁(8) is 7-dimensional: E_G(τ), E_G(2τ), E_ζ(τ), E_ζ(2τ), E^{χ₋₈,1}, E^{1,χ₋₈}, and the CM cusp form g.
- Two are excluded:
  - E^{χ₋₈,1}: it has an ∞ at the cusp 0 that nothing can cancel (the √2 square class);
  - g: it brings in L(g,2).
- The G-part must be W_E = E_G(τ) − 8E_G(2τ) (layer theorem).
- The constant terms −¼, −¼, −3/2 must cancel.
- L-values at s = 2: W_E → G/2; E_ζ(dτ) → π²/(12d²); E^{1,χ₋₈} → π²/6.
- Pure G leaves span{W_E, Z}, with Z = E_ζ(τ) + 8E_ζ(2τ) − (3/2)E^{1,χ₋₈}. Z has constant term 0 and π²-part 0 (exact).

**The digits.** For W = αW_E + γZ: R₁(ε) = αε − γ/2 and R₂(ε) = −4αε + (15/2)γ.
- At most one class is saved at digit 1, since R₁(+1) − R₁(−1) = 2α ≠ 0.
- γ = 2α saves the split class. Then R₂(+1) = 11α ≠ 0, so digit 2 is always paid. Its −4 is case E's half turn (τ + ½).
- So τ = 1 (the unsaved class, from digit 1) + ½ (the saved class, from digit 2) = 3/2. This is the floor for pure G.
- With π² allowed there is one more parameter, which clears digit 2 for one class: τ = 1 + ⅓ = 4/3.
  - W_B = W_E − E_ζ(τ)/5 + 5E_ζ(2τ) − (4/5)E^{1,χ₋₈}, limit G/2 − 11π²/240.
  - Digit 3 is paid (8/5).

**The new family, built and checked.**
- W*₊ = W_E + 2Z saves the split class; W*₋ = W_E − 2Z saves the inert class.
- Both live on Γ₁(8), with T = (1 − θ₃(q²)/θ₃(q))/2 and f = θ₃(q)².
- Exact values: W*₊(p) = 3p²(1 − χ₋₈(p)) for split p (that is, 0 or 6p²), and W*₊(2p) ≡ 11 mod p.
- v_n/u_n → G/2 for both, to 80 digits at n = 300 (the working precision). This is the same limit as case E.
- Denominators of v_n (exact, n ≤ 300):
  - every split (resp. inert) prime in (n/2, n] has exponent 0, and everything else is ≈ 2;
  - log denom(v_n)/n = 1.45 at n = 300, against 1.995 for lcm² and 1.92 for case E.
- Cumulative, to 200000 (exact per prime): log D / log lcm² is
  - 0.9999 for case E;
  - 0.7499 for W*₊ and 0.7503 for W*₋;
  - 0.6662 for W_B.
- The first-digit rule holds for all 5118 primes in (50, 50000), with no exceptions.
- At n ≈ 300 the per-n bucket averages show no further digit savings (≈ 2 outside the saved class).

**Gate.**

| family | limit | τ | log R* | log R* − τ |
|---|---|---|---|---|
| case E | G/2 | 2 | 1.232 | −0.77 |
| W*₊, W*₋ (new) | G/2 | 3/2 | 1.232 | −0.27 |
| W_B | G/2 − 11π²/240 | 4/3 | 1.232 | −0.10 |

- The needed 37% (τ ≤ 1.232) is out of reach for Beukers G-families. The floor is a 25% saving for pure G and 33% with π².
- Higher levels break the harmless group: every invariant is ≥ N, and log R* ≈ 65/N² loses more than any extra digit saves.
- In carry language: the saving is a first-digit cancellation on one class of primes mod 4.
  - G's digit is χ₋₄(p) and the companions' digit is 1, so only one class cancels.
  - The second digit carries case E's half-turn weight −4, which the admissible companions cannot match for pure G.

**What this closes and what it opens.**
- Beukers G-families are now pinned at both ends: analytically by H(√8) (log R* = 1.232), and arithmetically by τ ≥ 3/2. The pure-G deficit is 0.27, down from 0.77.
- The mechanism is general. A companion Eisenstein series whose p-digit is 1 cancels G's digit χ₋₄(p) on one class, and a third character (χ₋₈) removes the companion's π².
- It can be tried wherever G's linear forms have Eisenstein-like rational parts, for example the hypergeometric forms.

# Pushing the companion mechanism hard: the exact modular floor, and why it does not transfer (2026-09-29; companion_k2.py). The user: "Push this hard."
**1. The modular floor.** For every congruence Beukers family whose limit is pure G, τ − log R* ≥ 3/2 − 1.2323 = 0.268. This assumes the standard independence of π², G and the Clausen values already used in the layer theorem.
- The group side:
  - the G-part is case E's pattern (layer theorem);
  - by decoupling, the harmless set lies inside the G-part's harmless set;
  - T must be in the group.
  - So the harmless group lies in ⟨T, L₈⟩ (or in its level-raised copies ⟨T, L_{8d′}⟩, which are worse).
- With L₈ in the harmless group, log R* = 1.2323 and every component must be L₈-invariant, i.e. of level 8. That is the 7-dimensional weight-3 space of Γ₁(8), whose optimum is τ = 3/2 (previous section).
- Without L₈, the single-syllable term 2ζ(2)/64 = 0.051 of S = 0.062 is lost, so log R* < 1. Yet τ ≥ 1 still holds, because one class always pays at digit 1: G's digit χ₋₄(p) is odd under the class swap, and every pure-G companion's digit is class-independent.
- So the minimum is 0.268, attained by W*±. With π² allowed it is 4/3 − 1.232 = 0.101 (W_B).

**2. Why the mechanism is free only here.** It needs a companion with three properties:
- (a) first digit 1 at every prime;
- (b) a constant in ℚ·π²;
- (c) the same homogeneous part and singular set, so the decay is unchanged.

In GL₂/ℚ these are exactly the Eisenstein series E^{1,ψ}. The divisor d = 1 gives the p-coefficient 1 + ψ(p)p² ≡ 1, while the L-value ζ(2)L(0,ψ) is rational × π². Only Hecke eigenforms have this "trivial first digit with a character-flavoured constant".

**3. The transfers, checked.**
- **K¹, Zudilin's G-recurrence.** It has order 2 and no singular index (A(n) = (2n+1)²(2n+2)²p(n) ≠ 0). So every sequence obeying it from some index on is a ℚ-combination of u and v, and every inhomogeneous companion has its limit in ℚ + ℚG. A π² companion does not exist.
- **K², Catalan's Hankel machine (exact, companion_k2.py).** The kernel is Catalan + c·plain, in both the pole values and the moments; the plain unknown is set to 0, since only the content is tested.
  - At a pole b > p the 1/p² coefficient is −(b/4)(2χ₋₄(b)χ₋₄(p) + c). So c = −2 cancels it exactly for the poles b ≡ p (mod 4), half the poles of every prime.
  - Content exponents for K < p < 2K at K = 30, as sums over the primes of each class (Catalan's law 3n + 1 gives 72 for p ≡ 1 mod 4 and 88 for p ≡ 3 mod 4):
    - c = 0 (control): 72 / 88, Catalan's law exactly;
    - c = −2: 75 / 88 (every p ≡ 1 mod 4 costs one more);
    - c = +2: 72 / 92;
    - c = −1, 1, −4, ½: 98 / 120 (the plain machine's cost).
  - No weight gains. The ledger rule predicted this: per mirror pair,
    - saved: 2 − 1 (Hermite, both members p-integral) + 1 (von Staudt from the Bernoulli moments) = 2;
    - unsaved: 2 + 1 + 1 = 4;
    - average 3, Catalan's own cost.
  - The digit-1 companion is forced to be the plain residue pattern, whose Bernoulli moments bring the von Staudt prime back. It also puts π² into the pole values with pole-dependent weights, i.e. a second unknown.
  - Companions with a nontrivial character are von Staudt-free, but their digit is ψ(p) and their constant is √d·π² (even ψ) or a new non-critical value (odd ψ).
- Outside the modular world the balance law holds exactly: the saving and its cost are the same prime, counted twice.

**4. The one door left on this route.** A holomorphic family carrying G with a harmless parabolic of invariant w ∈ {5, 6, 7}.
- log R* would be 5.13, 2.68, 1.73. The first two pass even at τ = 2; w = 7 passes with the companion saving (1.73 > 3/2).
- The conductor 4 forbids this for congruence forms.
- A non-congruence family pays in denominators (Calegari–Dimitrov–Tang). Whether that payment can stay below the analytic gain is the only quantitative question left on the Beukers route.

# The non-congruence door, gone through at its smallest instance (2026-09-29; door_groups.py, door_belyi.py, door_family.py, door_recurrence.py). The user: "Go for it"
**A. Bounded denominators cannot do it (exact).**
- In SL₂(ℤ/2^a), −I lies in ⟨T, L_w⟩ exactly when 4 ∤ w. Checked for w ≤ 16 and a = 2…6.
  - UPDATE (2026-10-01, real date): now proved for all w and all a ≥ 2. See the section "The chat-side handoff of 2026-09-30" below.
- G needs the χ₋₄ component at 2, on which −I acts by −1 (odd weight).
- A multiplier that is trivial on T and L_w must therefore kill that component whenever 4 ∤ w.
- A family with bounded denominators is congruence (Calegari–Dimitrov–Tang), so it carries G only when 4 | w, i.e. w ≥ 8 by the lattice principle. This is the sign of ∞ at the prime 2 once more.
- So every family at widths 5, 6, 7 has unbounded denominators.

**B. The smallest candidates** (genus 0, three cusps, width 1 at ∞, width w at 0 = S∞; all non-congruence by the order of the monodromy group):

| index | elliptic points | widths | monodromy group |
|---|---|---|---|
| 9 | one of order 2 | (1, 1, 7) ×2 | order 504 |
| 9 | one of order 2 | (1, 2, 6) | order 54 |
| 9 | one of order 2 | (1, 3, 5) | A₉ |
| 10 | one of order 3 | (1, 2, 7) | S₁₀ |
| 10 | one of order 3 | (1, 4, 5) | S₁₀ |

**C. The (1, 2, 6) group, built exactly.**
- It is the index-3 subgroup of Γ₀(2) cut out by the non-Galois cubic cover y = t/(1 − 256t/27)³, with y = (η(2τ)/η(τ))^24 and j = (1 + 256y)³/y.
- The Belyi map is over ℚ, and its passport 3³ / 2⁴1 / (1, 6, 2) is checked by factorisation:
  - cusps at t = 0 (width 1), 27/256 (width 6, the cusp 0) and ∞ (width 2);
  - the elliptic point at 27/64.
- Weight 1: f = E₄^{1/4}(A(0)/A(t))^{1/4}, with A = 256t + (1 − 256t/27)³.
  - f is nonzero at the harmless cusps.
  - At the width-2 cusp it has order ¾, i.e. multiplier −i there: the quarter turn again.
- The family u_n = [tⁿ]f obeys an exact Apéry-like recurrence (checked on 43 terms). With w_n = (27/4)ⁿu_n:
  (n+1)² w_{n+1} = (80n² + 72n + 21) w_n − 64 (4n−1)² w_{n−1},   w = 1, 21, 3057/4, 135469/4, …
  - The characteristic roots are 256/27 and 64/27, i.e. exactly the width-6 cusp and the elliptic point.
  - The quarter shift (4n − 1)² is the −i multiplier.
- The harmless pair (∞: 1, 0: 6) gives log R* ≥ 2.68.

**D. What it carries.**
- **Limits.** One new constant, Λ₁ = 0.019056859720282146516351399462360887015024975475…, computed to 390 digits (n = 700).
  - The inhomogeneities t² and t³ give 440208Λ₁ − 1404928Λ₂ = 6561 and 729Λ₁ − 13095Λ₂ + 30976Λ₃ = 0.
  - At 300 digits, no combination of Λ₁, Λ₂, 1 lies in the ℚ-span of G, π², Γ(¼)⁴/π² or its inverse. Singly, Λ₁ and Λ₂ also avoid L(2,χ₋₃), L(2,χ₋₈), ζ(3), π log 2, log² 2, log 3 and π√3.
  - **G does not appear.**
- **Denominators.**
  - t(q), u_n and v_n carry 3-adic denominators ≈ 3^{2.93n}, from the cubic cover: about 3.2 per n. They are 2-adically integral.
  - So τ ≈ 5.2 against log R* ≥ 2.68: a failure by about 2.5 even before G.
  - CORRECTION (2026-10-01, real date; found by the user's chat-side session): the rate is exactly 3 per n, i.e. 3·ln 3 = 3.296, so τ ≈ 5.3 and the failure is about 2.6. The figure 2.93 was a fit at n = 40. The exact law is v₃(u_n) = −3n + (a borrow count); see the section "The chat-side handoff of 2026-09-30" below.

**Reading.**
- This family is non-congruence only at 3. At 2 it is bounded, so part A's 2-adic argument applies (−I ∈ ⟨T, L₆⟩ mod 2^a), and G is indeed absent.
- A G-family at width 5–7 would have to be non-congruence at 2 itself. There is no mechanism that makes such a family's periods Tate, and the example's single constant is new.
- The door is closed at its smallest instance on both counts: no G, and the arithmetic cost exceeds the analytic gain.
- The remaining unbuilt candidate, (1, 3, 5) with monodromy A₉ and log R* ≥ 5.13, is predicted to behave the same way. It would need unbounded 2-adic denominators and G among its periods, and neither has a mechanism.

# The chat-side handoff of 2026-09-30, checked and merged (2026-10-01, real date; handoff_checks.py, handoff_checks_out.txt). The user: "A small handoff doc on some side work."
**Source.** The user's chat-side session of 2026-09-30 recomputed parts of this repo with separate code and wrote a handoff. In it, "verified" means recomputed from the written definitions, and "inference" means argued, not proved. This section records it, re-checks what could be re-checked here, and answers its open items where this log already has the answer. (Earlier day labels in this log ran ahead of the calendar, up to "10-01". This section is dated by the real clock.)

**1. The 2-adic proof (twoadic_law_proof.md): an independent check.**
- The chat side rebuilt the Catalan machine in the v-scale: moments μ_k = (k+1)T_{2k+1}/2, nodes −(2j+1)², pole values F_j(X) = −(−1)^j(2j+1)(β_j + β_{j+1})/4 + X·(−1)^j(2j+1)/2.
- Verified there for every K = 3…14:
  - the content law e₂(K), ties included (ties at K = 6 and 12);
  - the Smith form of A(0): the K smallest of band ∪ moment pivots, top divisor raised at ties;
  - the Smith form of A(C) = {v₂(h_q)}, which tests Theorems 6 and 7 end to end;
  - the coefficient law v₂(q_{K−M}) for every M;
  - the tangent-number and Boole/Mahler formulas for C agree to 2^200.
- Read line by line there, no gap found: Lemma 12.1, Theorem 7, Theorem 8, the Cauchy–Binet consequences.
- Not checked there: the twisted machine numerically, Theorem 6's analytic steps, section 9 in detail.
- Reported there: A(X) = A(C) + (X − C)B gives the machine polynomial a 2-adic near-root at C of depth 204, 693, 2159, 9039 bits at K = 8, 14, 24, 48, tending to about 4K².
  - The K = 8 value agrees with this repo's formulas: Σ_{q<8} h_q = 175 (the paper's 2¹⁷⁵) plus e₂(8) = 29.
- Their reading matches ours: the law proves the 2-part of a height that was already measured, and the real gap is closeness. (Their "+0.38 per K²" is the K = 40 measurement; the asymptotic value here is +0.585.)

**2. The name of C, and whose continued fraction the twisted machine is (credit).**
- In the Kubota–Leopoldt convention, C = ζ₂(2).
  - At p = 2 the Teichmüller character is χ₋₄, so L₂(1−n, 1) = L(1−n, χ₋₄) = β(1−n) for odd n. Hence ζ₂(2) is the 2-adic limit of β(−2k) = E_{2k}/2 as 2k → −2.
  - Kubota–Leopoldt's L₂(s, χ₋₄) vanishes identically (odd character). Beukers says exactly this about Calegari's "2-adic Catalan constant" (arXiv math/0603277, §1, read in the source).
  - So "Calegari's L₂(2, χ₋₄)" in this log and ζ₂(2) are the same number under two names.
- Checked here with this repo's C (the tangent-number series): v₂(E_{2^N−2} − 2C) = N for N = 3…9.
- Chat side: C matches Beukers' ζ₂(2) = H₂(2,1,4) + H₂(2,3,4) to 318 bits.
- **Credit.** The twisted machine's J-fraction (λ_n = n⁴, b_n = 2n² + 2n + ¾) is Beukers' continued fraction for Θ(x) at x = ½.
  - His recurrence is U_{n+1} = (2n² + 2n + 1 − x + x²)U_n − n⁴U_{n−1} (§6, read in the source).
  - He states that it converges 2-adically to Θ₂(½) = −8ζ₂(2), and that its convergents are Calegari's approximations.
  - So Theorem 1 and Theorem 6 for the twisted machine re-derive a known framework.
  - The pencil A(X) = A(C) + (X − C)B and the Smith-form and content laws are not in Beukers' paper. That is the only paper either side checked.
- Also in Beukers §6 (his remark, not re-derived here): the same Padé approximants at x = −n + ½ are Rivoal's Catalan approximations. So the twisted K² machine and the Rivoal–Zudilin K¹ forms are two specialisations of one Padé table, that of Θ(x).

**3. The width criterion is now a theorem.**
- **Theorem.** For every w ≥ 1 and a ≥ 2, −I lies in ⟨T, L_w⟩ ⊂ SL₂(ℤ/2^a) exactly when 4 ∤ w.
- Proof (the chat side's). Write w = 2^k·m with m odd. Mod 2^a, L_w = (L_{2^k})^m and m is invertible, so L_w and L_{2^k} generate the same cyclic group.
  - k = 0: T and L₁ generate SL₂(ℤ), which maps onto SL₂(ℤ/2^a).
  - k = 1: (T⁻¹L₂)² = −I exactly, already in SL₂(ℤ).
  - k ≥ 2: every generated matrix is (1 ∗; 0 1) mod 4, and −I is not.
- Brute force agrees for w ≤ 32 and a = 2…6 (both sides).
- **Lemma (product closure; added here).** The closure of ⟨T, L_w⟩ in SL₂(ℤ̂) is A × B, where A is its closure in SL₂(ℤ₂) and B its closure in ∏_{p odd} SL₂(ℤ_p).
  - By Goursat the closure is a fibre product of A and B over a common quotient Q.
  - The kernel of A → SL₂(𝔽₂) is pro-2, and the image is S₃ (w odd) or of order 2 (w even). So every nontrivial finite quotient of A has a quotient of order 2.
  - B has no quotient of order 2. In ∏_{p odd} SL₂(ℤ_p) the square roots T^{1/2} = (1 ½; 0 1) and L_w^{1/2} are limits of powers of T and L_w, so they lie in B. A continuous map B → ℤ/2 therefore kills T and L_w, hence all of B.
  - So Q = 1.
- **Corollary (the congruence case).** Let W be a congruence form of odd weight with nebentypus χ₋₄, fixed by T and L_w. If 4 ∤ w, then W = 0.
  - The image of ⟨T, L_w⟩ mod the level N = 2^a·M contains (−I mod 2^a, I mod M).
  - That element comes from some γ ∈ Γ₀(N) with d ≡ −1 mod 2^a and d ≡ 1 mod M, and it acts on W by χ₋₄(d) = −1.
- The premise "G needs nebentypus χ₋₄" is the independence assumption of the layer theorem: L(E^{χ₁,χ₂}, 2) = L(2,χ₁)·L(0,χ₂), which is a rational multiple of G only for χ₁ = χ₋₄ and χ₂ principal. The chat side lists this premise as not checked. It is an assumption here too, not a theorem.
- **The same at 3.** −I lies in ⟨T, L_w⟩ ⊂ SL₂(ℤ/3^a) exactly when 3 ∤ w.
  - If 3 ∤ w, the group is everything. If 3 | w, it is unipotent mod 3.
  - Checked for w ≤ 27, a ≤ 3. The product lemma holds with cube roots in place of square roots.
  - So χ₋₃ survives exactly at widths divisible by 3. The least such width ≥ 5 is 6, which is Γ₀(6), CDT's case. For G it is 8.
  - For the harmless pair (∞, 0), i.e. for groups containing T and L_w, this proves the "least multiple of the conductor that is ≥ 5" rule on every congruence group, for conductors 3 and 4. The earlier table covered only Γ₀(N) and Γ₁(N). Other harmless pairs (∞, a/c) are not covered by this argument.
- The chat side's stronger inference, that non-congruence families at widths 5–7 lose G as well, stays an inference. A non-congruence form is not in the span of the congruence Eisenstein series, so the corollary does not reach it. The status is unchanged: no mechanism, and one worked example.

**4. The index-9 (1, 2, 6) family: simplification, an integer sequence, two carry laws, and a correction.**
- **Simplification (chat side; checked here exactly on 37 q-coefficients).** f = θ₄(2τ)²·(1 − s)^{−3/4}, with s = 256t/27 and θ₄(2τ) = Σ(−1)ⁿq^{n²} = η(τ)²/η(2τ).
  - Reason: E₄ = (1 + 256y)·θ₄(2τ)⁸ on Γ₀(2), and 1 + 256y = A(t)/(1 − s)³.
  - θ₄(2τ)² is the χ₋₄ theta-square with its whole zero at the cusp 0 (it is θ₃(τ + ½)²). It has order ¾ at the width-6 cusp, and the twist (1 − s)^{−3/4} cancels exactly that. So f is a nonzero constant there.
  - 64t + (1 − s)³ = −(s + ½)²(s − 4) (chat side). So the cover's last branch point is the elliptic point of Γ₀(2), and the group is a genuine index-3 subgroup with widths (1, 2, 6).
- **An integer sequence.** a_n := 27ⁿu_n = 256ⁿ[sⁿ]f = 1, 84, 12228, 2167504, 423881700, …, with
  (n+1)²a_{n+1} = 4(80n² + 72n + 21)a_n − 1024(4n − 1)²a_{n−1}.
  - **Proved: a_n ∈ ℤ.**
    - θ₄(2τ)² = Σ_k b_k y^k with b_k ∈ ℤ, because q ∈ y + y²ℤ[[y]].
    - y = (27/256)·s·(1 − s)^{−3}, so y^k = (27/256)^k Σ_{m≥k} C(m + 2k − 1, m − k) s^m.
    - (1 − s)^{−3/4} = Σ_j c_j s^j with c_j = ∏_{i≤j}(4i − 1)/(4i). Here 8^j c_j ∈ ℤ: its 2-adic valuation is s₂(j), and c_j = ±C(−3/4, j) is p-integral for odd p.
    - So a_n = 256ⁿc_n + Σ_{k≥1} Σ_{m=k}^{n} b_k·27^k·256^{n−k}·C(m + 2k − 1, m − k)·c_{n−m}, and every term is an integer.
  - **Proved: a_n ≡ 256ⁿc_n = 64ⁿ∏_{i≤n}(4i − 1)/n! (mod 27).** The terms with k ≥ 1 carry 27^k.
  - **3-adic carry law.** By Kummer, v₃(C(−3/4, n)) is the number of borrows in the base-3 subtraction (…2020)₃ − n, since −3/4 = (…2020)₃.
    - Proved: v₃(a_n) equals that borrow count whenever the count is ≤ 2 (862 of the n ≤ 3000).
    - Checked: equality for every n ≤ 3000 (maximum 7).
  - **2-adic carry law (checked, n ≤ 3000).** v₂(a_n) = 2·s₂(n), the Kummer count for C(2n, n)².
- **CORRECTION of the door section.** v₃(u_n) = −3n + v₃(a_n), so the 3-adic rate is exactly 3 per n: 3·ln 3 = 3.296, not "about 3.2", and τ ≈ 5.3.
  - The chat side found this from n < 58. The congruence above proves it.
  - Equivalently, the family is integral in the variable t/27, where the conformal radius is 27 times smaller.
- **The chat side's L.** L := lim v_n/w_n with v₀ = 0, v₁ = 1 is 0.059091254764002891259696550881531352…, and Λ₁ = (1701·L − 81)/1024. Confirmed here to 62 digits.
  - Their relation search at 50 digits (coefficients up to 10⁶) against 1, G, π², Γ(¼)⁴/π² is negative, like ours at 300 digits.
- Not checked on the chat side: log R* ≥ 2.68. It comes from uniformize_radius.py at w = 6, and it uses that F is analytic at the width-6 cusp: f is nonzero there (above) and W vanishes there.

**5. Width 8 against width 6: the handoff's open item 1 is answered here.**
- Their comparison (re-verified on their side to 40 and 37 digits):

| family | group, cusp widths | singular points | ratio (convergence rate) |
|---|---|---|---|
| χ₋₃, Zagier (10, 9, 3), limit L(2,χ₋₃)/2 | Γ₀(6): 6, 3, 2, 1 | 1/9 and 1 | 9 |
| χ₋₄, Zagier (12, 32, 4), limit G/2 | Γ₀(8): 8, 2, 1, 1 | 1/8 and 1/4 | 2 |
| index-9 family | (1, 2, 6) | 1/64 and 1/16 | 4 |

- Their question: does the Calegari–Dimitrov–Tang bound close on Γ₀(8) with ratio 2 instead of 9? It needs log R* for the (1, 8) pair with the two extra cusps, "not computed".
- It is computed in this log ("Uniformized at infinity", "The half measure"): log R* = 2π²Σ1/c² over ⟨T, L₈⟩ = 1.2323 (S = 0.062427).
  - The two extra cusps are not in the harmless group: ¼ pays (the quarter turn) and ½ is the pole.
- **Answer: no.** The necessary gate log R* > τ fails.
  - Zagier's (12, 32, 4) has τ = 2 (its true denominators are the full lcm², checked to n = 600): 1.2323 − 2 = −0.77.
  - The best pure-G family on the same curve (W*±, τ = 3/2) gives −0.27, and that is the floor for every congruence pure-G family.
  - CDT's case at width 6: 2.68 − 2 = +0.68.
- Their open item 2, the two-unknown (π², G) functional: unchanged. Its K¹ counterpart here is W_B (limit G/2 − 11π²/240, τ = 4/3), which misses the gate by 0.10.
- Their open item 3, a full-rate kernel whose pole values involve G only: unchanged. The companion test in this log (plain-kernel companion, exact at K = 30) washes exactly.

**6. Closed on the chat side (recorded so that they are not repeated).**
- Evens, odds and primes with Ramanujan summation.
  - The parity sums regularise to Euler numbers, and the log-weighted prime sum is β′/β.
  - The surviving identity, 3·log 3 − 5·log 5 + 7·log 7 − … = 2G/π, is the functional equation of β (this log's β′(−1) = 2G/π). It gives no approximations.
  - Bare prime sums stop at the natural boundary. Krapivsky–Luck (arXiv 2606.24536) get past it only by a chosen regularisation.
- Splitting G into two rational parts plus one irrational part: the rational parts carry nothing. The only three-term version with content is 1, π², G, the existing two-unknown construction.

**7. Standing constraint (the user, through the handoff).** No p-adic Catalan targets as a direction. Items 1 and 2 above are recorded because the Hankel machine's own denominators led there, not as a proposal.

# Ramanujan summation on the dodecahedron (2026-10-01, real date; ramanujan_dodecahedron.py, ramanujan_dodecahedron_out.txt). The user: "How does ramanujan summation apply to the dodecahedron?"
**The solid.** The dodecahedron is the modular curve X(5) = H/Γ(5) with its cusps filled in (Klein).
- PSL₂(𝔽₅) ≅ A₅, the 60 rotations.
- The 12 cusps, each of width 5, are the face centres. The 20 points over j = 0 are the vertices. The 30 points over j = 1728 are the edge midpoints.
- T is the rotation by 72° about a face, S the half turn about an edge, ST the rotation by 120° about a vertex.
- The cusp of the modular group is a face with infinitely many sides, because T has infinite order. Read mod 5 it is a pentagon: 1 + 1 + 1 + ⋯ closes after 5.

**1. Euler's formula is Ramanujan's value.**
- χ(X(N)) = |PSL₂(ℤ/N)|·(1/N − 1/6), with 1/6 = −2ζ(−1). The orbifold Euler characteristic of the modular group is 2ζ(−1) = −1/6, and each cusp gives back 1.
- N = 5: 60·(1/5 − 1/6) = 2. N = 3 and 4 also give 2 (tetrahedron, octahedron). N = 6 gives 0, and N = 7 gives −4 (Klein's quartic).
- So the Platonic solids exist exactly while 1/N > 1/6. The modular group itself is the case N = ∞: 1/3 − 1/2 + 0 = −1/6.

**2. The coordinates are Ramanujan sums.**
- Regularised sums over residue classes mod 5: Σ_{n ≡ a (5)} n = 5ζ(−1, a/5) = −(5/2)B₂(a/5) (exact).
  - 1 + 6 + 11 + ⋯ = 4 + 9 + 14 + ⋯ = −1/60.
  - 2 + 7 + 12 + ⋯ = 3 + 8 + 13 + ⋯ = 11/60.
  - 5 + 10 + 15 + ⋯ = −5/12, and the total is −1/12.
- These are the exponents of the Rogers–Ramanujan functions: q^{−1/60}G(q) and q^{11/60}H(q), with G = ∏1/((1 − q^{5n+1})(1 − q^{5n+4})) and H = ∏1/((1 − q^{5n+2})(1 − q^{5n+3})).
- Their difference, 1/5 = 12/60 = −½L(−1, χ₅), is the exponent of Ramanujan's continued fraction r(τ) = q^{1/5}H/G = q^{1/5}∏(1 − qⁿ)^{(n/5)}.
- Checked to 40 digits: with these exponents the pair transforms under τ ↦ −1/τ by (2/√5)·(sin 2π/5, sin π/5; sin π/5, −sin 2π/5).
  - Hence r(−1/τ) = (1 − φr)/(φ + r) and r(τ + 1) = e(1/5)·r. So r is Klein's icosahedral coordinate.
  - Ramanujan's r(i) = √((5 + √5)/2) − φ is the fixed point of that half turn: an edge midpoint.
- The face equation 1/r⁵ − 11 − r⁵ = (η(τ)/η(5τ))⁶ holds exactly (60 coefficients).

**3. In this program the dodecahedron is ζ(2)'s solid.**
- t = r⁵ is the Hauptmodul of Apéry's ζ(2) family on Γ₁(5), the quotient of the dodecahedron by the 72° rotation.
  - Σu_ntⁿ = 1 + Σc(n)qⁿ/(1 − qⁿ), with c = (3, 1, −1, −3) for n ≡ 1, 2, 3, 4 mod 5 and u_n = Σ_k C(n,k)²C(n+k,k) (exact, 60 coefficients).
- Its singular points, the roots of t² + 11t − 1, are φ⁻⁵ and −φ⁵. They are the two rings of five faces around the top and bottom faces (|r| = φ⁻¹ and φ).
- Apéry's inequality is φ⁵ = 11.09 > e² = 7.39. In the uniformized form: width 5 and log R* = 5.13 against τ = 2.

**4. The three Platonic faces are the three regimes of the lattice principle.** tr(T·L_w⁻¹) = 2 − w.

| face | w | T·L_w⁻¹ | ⟨T, L_w⟩ | constant living at that conductor | outcome |
|---|---|---|---|---|---|
| triangle (tetrahedron) | 3 | elliptic | lattice Γ₀(3) | L(2,χ₋₃) | doubled to 6: passes (CDT) |
| square (cube, octahedron) | 4 | parabolic | lattice Γ₀(4) | G | doubled to 8: fails |
| pentagon (dodecahedron) | 5 | hyperbolic, eigenvalues −φ^{±2} | thin | ζ(2) on Γ₁(5) | usable as it is: Apéry |

- The pentagon is the only Platonic face that is not a lattice, so it is the only one usable without doubling.
- G's solid is the cube and octahedron, in two senses:
  - χ₋₄ has conductor 4, and X(4) is the octahedron;
  - the regular ideal octahedron has volume 4G (it is the Whitehead link complement).
- χ₋₃'s solid is the tetrahedron in the same two senses: X(3), and the regular ideal tetrahedron has volume (3√3/4)·L(2,χ₋₃).
- So the dodecahedron does not carry G: its conductor is 5, and G needs 4. It shows what G lacks, a face of width ≥ 5 at its own conductor.
- The doubled triangle is the hexagon, N = 6, the flat case where 1/N equals Ramanujan's 1/6; CDT's proof passes narrowly there. The gate's own threshold (2π²S(w) = 2) lies between w = 6 and 7. The two thresholds are close, and I have no proof that they are the same thing.

**Arithmetic remark, with no derivation.** The user's 113/120 equals 1 − 1/24 − 1/60: η's exponent and the dodecahedron's exponent.

# Squaring the faces: the six families form one flip tree, and G is one flip past CDT (2026-10-01, real date; flip_families.py, flip_families_out.txt). The user: "Can we not square each face and times by 12 to normalize. The carry will be the equalateral triangle part left over after normalization. Should be 12 times pi plus that bit left over"
**Reading.** A face is a cusp, and its size is the cusp width: the number of ideal-triangle corners that meet there. Squaring a face means bringing its width to 4. The carry is the corners removed, and they have to go to other faces.

**The normalisation by 12 is exact.**
- Each of our families lives on a sphere with 4 cusps and no elliptic points, of index 12 in PSL₂(ℤ). It is tiled by 4 ideal triangles, each of area π, with 12 corners in all.
- So the widths add up to 12: 5 + 5 + 1 + 1 (Apéry), 6 + 3 + 2 + 1 (CDT), 8 + 2 + 1 + 1 (G).
- Summing each face's fan of triangles gives 12π for every family. Nothing is left over in the total: a carry only moves corners from one face to another.
- In angle measure the user's split is exact for one face: a pentagon is a square plus a triangle, 3π = 2π + π. In Euclidean area it is not (regular pentagon − unit square = 0.7205, equilateral triangle = 0.4330).

**All such surfaces (exact enumeration).** There are exactly six up to conjugacy, and all six are congruence: Beauville's six families.

| widths | group | constant |
|---|---|---|
| (5, 5, 1, 1) | Γ₁(5) | ζ(2), Apéry |
| (6, 3, 2, 1) | Γ₀(6) | L(2,χ₋₃), CDT |
| (8, 2, 1, 1) | Γ₀(8) | G, case E |
| (4, 4, 2, 2) | Γ₀(4) ∩ Γ⁰(2) | G again: case E in the variable 2τ |
| (9, 1, 1, 1) | Γ₀(9) | level 9 (Zagier's case B) |
| (3, 3, 3, 3) | Γ(3) | level 9 in the variable 3τ; the tetrahedron |

**The carry is a flip.** A flip replaces one edge by the other diagonal of its quadrilateral. The two ends of the old edge lose a corner each, and the two ends of the new edge gain one. (A loop loses two at its one vertex.)
- The flip graph of the six is a tree: (3,3,3,3) – (4,4,2,2) – (6,3,2,1) – (5,5,1,1), with the branch (6,3,2,1) – (8,2,1,1) – (9,1,1,1).
- Up to the isogenies τ ↦ 2τ and τ ↦ 3τ it is a path: Apéry – CDT – G – level 9.
- Squaring both pentagons, (5,5,1,1) → (4,4,2,2), takes two flips and passes through CDT's hexagon. The two carried triangles land on the two unit cusps and turn them into digons.
- So the squared dodecahedral family is exactly G's own family.

**What the carry costs.** The gate number is the invariant of an adjacent pair of cusps, which is the product of their widths.
- Smallest products: 1·5 = 5 (Apéry); 1·6 = 2·3 = 6 (CDT); 1·8 = 2·4 = 8 (G, both frames); 1·9 = 3·3 = 9 (level 9).
- log R*: 5.13, 2.68, 1.23, 0.93 (the last from uniformize_radius.py at w = 9, S ≈ 0.0471), against τ = 2, or τ = 3/2 for the best G-family.
- Each flip away from Apéry's family raises the invariant: 5, 6, 8, 9. The gate passes at 5 and 6 and fails from 8 on.
- The carry is the whole loss: it turns the pair (1, 5) into (2, 4).
- The 2-adic form of the same statement is the width criterion: L₅ ≡ L₁ mod 4, and the remainder after removing squares is what brings −I into the group.

**Scope.** A flip is combinatorial. It does not map the forms of one family to another, so it does not carry a proof across. It shows where G sits: one flip past the family that CDT proved, and two past Apéry's.

# The radii to 25 digits by a transfer operator; their sum is not a relation (2026-10-01, real date; radius_transfer.py, radius_transfer_out.txt). The user: "Those log R add up to 9.97. is that relevant at all?"
**Answer: no.**
- The four numbers belong to four different surfaces. Each is compared with its own family's τ, and no criterion adds them.
- With exact values the sum is 9.9657, not 10 and not π² = 9.8696. (The rounded table values, one of them an extrapolation, added to 9.97.)

**The exact values: a new method.** S(w) = g(0), where g and Ψ solve the linear functional equations
    g(x) = Σ_{m≠0} v²·(1 + Ψ(v)), v = 1/(x + wm);   Ψ(u) = Σ_{k≠0} (k + u)⁻²·g(1/(k + u)).
- Derivation. Let G(c, d) be the total of 1/c² over all continuations of a word whose bottom row is (c, d) and whose next syllable is a power of L_w. The two moves are (c, d) ↦ (c + dwm, d) and (c, d) ↦ (c, d + ck). G is homogeneous of degree −2, so G(c, d) = d⁻²g(c/d), and the recursion for G is the pair of equations above.
- Ψ lives on [−ρ, ρ], with ρ = (1 − √(1 − 4/w))/2, and g on [−1/(1 − ρ), 1/(1 − ρ)]. Both are analytic well beyond these intervals.
- For a polynomial the sums over m and k are exact Hurwitz zeta values: Σ_{m≠0}(x + wm)⁻ⁿ = w⁻ⁿ[ζ(n, 1 + x/w) + (−1)ⁿζ(n, 1 − x/w)]. So the fixed point is one linear solve at Chebyshev nodes.
- Between 36 and 44 nodes the values agree to 26 digits (w = 5) and to 42 digits (w = 9).

| w | S(w) | log R* = 2π²S(w) |
|---|---|---|
| 5 | 0.25965825266580560724 | 5.12544846657921541 |
| 6 | 0.13565555560902939819 | 2.67773333734219755 |
| 7 | 0.08765062438052809051 | 1.73015397628858056 |
| 8 | 0.06242676007845623084 | 1.23225485203216218 |
| 9 | 0.04712800926140793363 | 0.93026961524194359 |

- These confirm the word-enumeration values (0.260, 0.1357, 0.0877, 0.062427, 0.0471), which were extrapolated tails.
- The modular floor is 3/2 − 1.2322548520 = 0.2677451480.
- Sum over w = 5, 6, 8, 9: S = 0.5048685776, 2π²S = 9.9657062712.

**What is meaningful instead: the threshold width.** The operator works for real w, so the gate log R*(w) = τ can be solved for w.

| τ | which families | passes below w = |
|---|---|---|
| 2 | full lcm² denominators (Apéry, CDT, case E) | 6.6348 |
| 3/2 | the best pure-G family (W*±) | 7.3949 |
| 4/3 | the π²-mixed family (W_B) | 7.7485 |

- G's invariant is at best 8 in a congruence family. It misses the pure-G threshold by 0.61 in width and the π²-mixed threshold by 0.25.
- The integers below the thresholds are 5 and 6 (τ = 2), and also 7 (τ ≤ 3/2). No congruence G-family has invariant 7: by the width criterion 4 must divide it.

# Mixing families in the Calegari–Dimitrov–Tang bound: allowed, but it cannot beat the gate (2026-10-01, real date; mix_lemma.py, mix_lemma_out.txt). The user: "Can you use 2 pure G families and 1 other family and another 1 other family. Mix and match over 4 parts"
**Reading.** Use several functions at once, with different denominator rates, instead of one family. This is how the CDT bound is built, so the test is their theorem itself.

**The bound (arXiv 2408.15403, Theorem 2.5.1, read in the source).**
- Take m holonomic functions in ℚ[[x]], ℚ(x)-linearly independent, with denominator rates σ₁ ≤ … ≤ σ_m, all meromorphic on the disc after composing with a map φ, φ(0) = 0.
- Then m ≤ I(φ)/(L − τ), where:
  - L = log|φ′(0)|;
  - τ = (1/m²)·Σ(2i − 1)σ_i, a weighted average of the rates (it needs L > τ);
  - I(φ) = ∬ log|φ(z) − φ(w)| over the torus, which is ≥ L, with equality only for univalent φ.
- mix_lemma.py reproduces the paper's own values of τ: 3/4, 8/9, 69/50 (Theorem 2.8.4) and 191/49 (the 14 functions of Theorem A).
- For G, L ≤ log R* = 1.2322548520 for every admissible φ: the functions live on X′ = H/⟨T, L₈⟩, and φ lifts to it.

**The four-part mix.** Two pure-G functions of rate 3/2 (F₊ and F₋, from W*±), and two free functions with integer coefficients (1 and √(1 − 4t), which is holomorphic on X′ because t never takes the value 1/4 there).
- τ = (5 + 7)·(3/2)/16 = 9/8 = 1.125, which is below L. So the averaged denominators do drop under the radius.
- But the bound then allows I/(L − τ) ≥ 1.2323/0.1073 = 11.5 functions, and the mix has 4. No contradiction.
- F₊ − F₋ = 4·f·D⁻²Z has limit exactly 0, so it is a function that exists whether or not G is rational. The two pure-G families are one G-dependent function plus one free function of rate 2.

**Lemma (no rescue by mixing).** Put Φ(S) = |S|·τ(S) = (1/|S|)·Σ(2i − 1)σ_i for a multiset S of rates. Then Φ(S ∪ {σ}) − Φ(S) ≥ σ for every S and every σ ≥ 0.
- Proof. Let S have m elements and insert σ at position k. The claim is equivalent to 2mC − B − A ≥ m(m − 2k + 2)σ, with A = Σ_{i<k}(2i − 1)σ_i, B = Σ_{i≥k}(2i − 1)σ_i, C = Σ_{i≥k}σ_i.
  - 2mC − B = Σ_{i≥k}(2m − 2i + 1)σ_i ≥ σ·(m − k + 1)², since σ_i ≥ σ there and the coefficients are the first m − k + 1 odd numbers.
  - A ≤ σ·(k − 1)².
  - (m − k + 1)² − (k − 1)² = m(m − 2k + 2).
- Tested on 200000 random multisets: the minimum of Φ(S ∪ {σ}) − Φ(S) − σ is 0.
- **Consequence.** Split a set of functions as S₀ ⊔ S₁, where S₀ exists unconditionally and S₁ only under a hypothesis (here: G rational).
  - The theorem applied to S₀ gives I ≥ |S₀|(L − τ(S₀)). (If L ≤ τ(S₀) this is trivial.)
  - A contradiction needs |S|(L − τ(S)) > I. Subtracting, |S₁|·L > Φ(S) − Φ(S₀) ≥ Σ_{F∈S₁} σ_F.
  - So L must exceed the average rate of the hypothesis-dependent functions. If S₀ is empty, I ≥ L gives the same with |S₁| − 1 in place of |S₁|.
- In words: every free function added to dilute the average also adds to the count the bound has to beat, and the two effects cancel exactly at best.

**Verdict for G.**

| G-dependent functions | rate | L − rate |
|---|---|---|
| case E | 2 | −0.7677 |
| best pure-G families (W*±) | 3/2 | −0.2677 |
| π²-mixed (W_B; hypothesis: 1, π², G dependent) | 4/3 | −0.1011 |

- No mixture can produce a contradiction from Theorem 2.5.1, whatever the other families are. The margins are the same as for a single family.
- This is the precise form of CDT's Remark 11.1.17 ("e² > 16·(1/4) … definitely precludes").
- So mixing is not a third lever. The two levers stay the same: the radius, and the G-dependent functions' own denominators.
- Scope: proved for the bound of Theorem 2.5.1. The paper's refined bounds (Theorems 6.0.2, 7.1.6, 8.0.1) change the denominator term and the numerator; they are not covered by this lemma. CDT's own proof has L ≈ 5.08 against a top rate 4 in its variable, on the right side of the same inequality.

**An illustration of why the count matters.** The six free functions 1, √(1 − 4t), log(1 − 4t), log((1 + √(1 − 4t))/2) and the two logs times √(1 − 4t) have rates (0, 0, 1, 1, 1, 1).
- With one pure-G function, m = 7 and τ = 1.051, and the best-case bound L/(L − τ) is 6.80 < 7. That looks like a contradiction.
- It is not one: the six free functions alone give L/(L − τ) = 3.59 < 6 in the same best case. So I(φ) = L is impossible for any φ that carries them, and in fact I ≥ 6(L − 8/9) = 2.06.

**Noted in passing (from the source).**
- CDT's Lemma 2.11.7: the Bost–Charles integral of their bivalent map 8(z + z³)/(1 + z)⁴ is log 8 + 4G/π. Catalan's constant is a term in their own bound.
- Their Theorem 2.8.4 is about exactly the rate 3/2 (type [1..n][1..n/2]) with singular points {0, δ, 1, ∞}. Our W*± families have that shape in x = 4t, with δ = ½, but the change of variable costs 4ⁿ in the denominators. That is the same gap in another form.

# Digit games made exact: reversal, four slots, cross multiplication, and G's Wieferich primes (2026-10-02; digit_games.py, digit_games_out.txt). The user: "I was thinking about weiferich primes and palindromes. And foundational number theory. And cross multiplication to extract squares. And writing primes reversed to get its sister primes. But also using the same numbers to fill 4 slots in every combo."
All statements are for Zagier's case E (limit G/2): u_n integers, v_n rationals, (n+1)²x_{n+1} = (12n² + 12n + 4)x_n − 32n²x_{n−1}, u = 1, 4, 20, 112, …, v = 0, 1, 7, ….

**1. The same digits in every slot: the residues do not see the order.**
- Lucas: u_n ≡ ∏ u_d over the base-p digits d of n (mod p). No failures for p = 3, 5, 7, 11, 13 and n < min(p⁴, 2500).
- So u_n mod p is unchanged by every rearrangement of n's base-p digits. Reversal is one of them, and palindromes are its fixed points.
- **Digit law for the rational part.** With L lower digits, p^{2L}·v_n ≡ χ₋₄(p)^L · v_{leading digit} · ∏ u_{lower digits} (mod p). No failures (same range).
  - Without the sign it fails for p ≡ 3 mod 4 (1110 of 2394 at p = 7).
  - The sign is u_{p−1} ≡ χ₋₄(p) (mod p), which holds for every odd p < 2500.
  - Derivation: v_n = u_n·Σ_{k<n} 32^k/((k+1)²u_ku_{k+1}). The terms with p^L | k + 1 dominate, and Lucas gives u_{mp^L − 1} ≡ u_{m−1}·u_{p−1}^L. So the character whose L-value is the limit enters the rational parts through u_{p−1}.
- So over all arrangements of the same digits, v depends only on which digit leads.
  - Example: p = 7, digits 1, 2, 4, 5 in four slots, 24 arrangements. u_n ≡ 2 every time, and p⁶v_n takes exactly four values: 3, 0, 2, 6 for leading digit 1, 2, 4, 5.

**2. Cross multiplication extracts squares.**
- Neighbours: u_n·v_{n+1} − u_{n+1}·v_n = 32ⁿ/(n+1)², exactly, for all n < 2500.
- Reversal: p²·(v_{bp+a} − v_{ap+b}) ≡ χ₋₄(p)·(u_a v_b − u_b v_a) (mod p), and u_a v_b − u_b v_a = u_a u_b·Σ_{k=a}^{b−1} 32^k/((k+1)²u_ku_{k+1}) exactly. No failures in 4755 cases (p ≤ 47).
  - Reversing a two-digit index is a cross multiplication, and the cross product is a sum over the squares between the two digits.
- The four cusp widths of the six families: the products are 25, 36, 64, 16, 81, 9, all squares. (They are the squared orders of the torsion groups of the elliptic surfaces.)
  - In the four frames where the widths form a proportion a : b = c : d, the cross product ad = bc is 5, 6, 8, 9: the gate invariant of the flip section.

**3. Wieferich primes (verified).**
- 1093 − 1 = 444 in base 16, and 3511 − 1 = 6666 in base 8: both are repdigits, hence palindromes. They are (2¹² − 1)·4/15 and (2¹² − 1)·6/7.
- Reversed in base 2, both give primes: 1093 ↔ 1297 and 3511 ↔ 3803. In base 10: 3511 ↔ 1153 (prime), and 1093 ↔ 3901 = 47·83.
- Eisenstein: (2^{p−1} − 1)/p ≡ ½·(1 − 1/2 + 1/3 − ⋯ − 1/(p−1)) (mod p), checked for p < 200. So a Wieferich prime is one where the alternating harmonic sum vanishes mod p; it does for 1093 and 3511.
- In other bases (p < 2·10⁶, b ≤ 12) the palindrome pattern does not persist: none in bases 5 and 6. It is a feature of the two base-2 primes, not a law.

**4. G's Wieferich primes.** For every prime 5 ≤ p < 700:
    u_p ≡ 4 + 8·χ₋₄(p)·p²·E_{p−3} (mod p³),   p²·v_p ≡ χ₋₄(p) + p²·E_{p−3} (mod p³),
with E_n the Euler numbers.
- Hence v_p/u_p ≡ χ₋₄(p)/(4p²) − E_{p−3}/4 (mod p).
- The Euler number E_{p−3} plays the part of the Fermat quotient. The primes with one more power are exactly those dividing E_{p−3}: 149 and 241 below 700.
- The next is 2946901, confirmed through the quarter sum: Σ_{0<k<p/4} k⁻² ≡ 4χ₋₄(p)E_{p−3} (mod p) (checked for p < 700; Lehmer's congruence).
- E_{p−3}/2 = β(3 − p), and 3 − p ≡ 2 mod (p − 1). By Kummer's congruence this is G read mod p. (Recorded as bookkeeping. No p-adic target is proposed.)
- Sisters: 149 ↔ 941 (prime) in base 10. 241 ↔ 142 is not prime.
- Literature not checked. The u-congruence is probably known (supercongruences for Apéry-like numbers in the work of Z.-W. Sun and Z.-H. Sun). The v-congruence and the digit law with the character sign may be new.

**5. What this does and does not do.**
- It explains how a prime's class enters the rational parts (through u_{p−1}), and it gives the exact second digit.
- It does not change the denominator growth. The digit symmetries concern residues mod p, not exponents, and the Wieferich-type primes are far too rare (two below 700) to save anything.

# 5040/10000 = 0.504: the birthday number of base ten, and its carry (2026-10-02; birthday_504.py, birthday_504_out.txt). The user: ".504. with that .004 being that carry I see pop up all the time. I got this from 10000 4 number conbos with repeaters. And I got 5040 with non repeaters. Take 5040/10000 gives us something universal I feel."
**What it is.** 10⁴ strings of four decimal digits, and 10·9·8·7 = 5040 with no repeated digit. So 0.504 = 63/125 is the chance that four decimal digits are all different: the birthday problem with 4 people and 10 days.
- 5040 = 7! because 10·9·8 = 6!. That is the identity 10! = 6!·7!.

**The carry, exactly.** b(b−1)(b−2)(b−3)/b⁴ = 1 − 6/b + 11/b² − 6/b³, with the Stirling numbers 6, 11, 6.
- At b = 10 this reads 1 − 0.6 + 0.11 − 0.006. The 11 does not fit in one decimal digit, so it carries.
- 1 − 6/10 + 10/100 = ½, and what is left is 1/100 − 6/1000 = 4/1000.
- So "½ plus a carry of .004" is literally correct in base ten. In base b the same split is (1 − 5/b) + ((11 − b)b − 6)/b³, and the first part is ½ only for b = 10.

**Universal or not.**
- The number is not universal. Four slots give 0.3499, 0.4102, 0.4609, 0.5040, 0.5409, 0.5729 in bases 7 to 12. Three and five slots in base ten give 0.72 and 0.3024.
- The law is universal: the chance is ∏_{j<k}(1 − j/b) ≈ e^{−k(k−1)/(2b)}, and it crosses ½ near k ≈ 1.18·√b.
- For four slots the crossing is at b = 9.9003. Base ten is the first base in which four different digits are more likely than not, which is why the value is ½ plus a little. Four is to ten what 23 is to 365 (0.4927).

**Where 5040 itself is a boundary.** Robin's criterion: the Riemann hypothesis is equivalent to σ(n) < e^γ·n·log log n for every n > 5040.
- The exceptions for 3 ≤ n ≤ 10⁶ are 3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24, 30, 36, 48, 60, 72, 84, 120, 180, 240, 360, 720, 840, 2520, 5040 (recomputed here).
- σ(5040)/5040 = 3.8381 against 3.8169.

**Also true, with no derivation linking it to the digit count.** 504 = 7·8·9 = −2/ζ(−5), the coefficient of E₆, and 5040 = 2·lcm(1..10).

**In this program.** 0.504 does not occur in the G numbers. The sum of the four radius sums is 0.50487, a different number. The other small excesses in this log (0.023 over log 2, 1/144 over 1/6) have their own causes; there is no single carry constant.

# Primes among digit strings: the count is the prime number theorem, the structure is 9·11·101 (2026-10-02; primes_in_slots.py, primes_in_slots_out.txt). The user: "Well look at primes. 1/4 of 100 are primes. And 168 in 1000. 1229 is first 9999 digits. Something is happening here. Push this a bit hard"
**1. How many.**

| digits k | π(10^k) | k·π/10^k | 10^k/π | step |
|---|---|---|---|---|
| 2 | 25 | 0.5000 | 4.000 | |
| 3 | 168 | 0.5040 | 5.952 | 1.95 |
| 4 | 1229 | 0.4916 | 8.137 | 2.18 |
| 5 | 9592 | 0.4796 | 10.425 | 2.29 |
| 6 | 78498 | 0.4710 | 12.739 | 2.31 |
| 7 | 664579 | 0.4652 | 15.047 | 2.31 |
| 8 | 5761455 | 0.4609 | 17.357 | 2.31 |

- The user's pattern is "one in 2k": ¼, about 1/6, about 1/8. It is exact at k = 2, off by 1.3 at k = 3 (166.7 against 168), by 21 at k = 4, and it drifts after that.
- The law is the prime number theorem: 10^k/π(10^k) ≈ k·log 10 − 1, so each extra digit adds log 10 = 2.3026 to the ratio. The measured steps are 1.95, 2.18, 2.29, 2.31, 2.31, 2.31.
- In the "half" form, k·π(10^k)/10^k tends to 1/log 10 = 0.4343, not ½. It would tend to ½ in base e² = 7.389. In bases 7 and 8 the value near 10⁸ is 0.547 and 0.514.
- 3·168 = 504 exactly, the same digits as 5040/10000. No derivation: 168 is a count of primes and 504 = 7·8·9.

**2. Which ones: primes against repeated digits (strings 0000 to 9999).**

| different digits | strings | primes | density |
|---|---|---|---|
| 4 | 5040 | 593 | 0.1177 |
| 3 | 4320 | 587 | 0.1359 |
| 2 | 630 | 49 | 0.0778 |
| 1 | 10 | 0 | 0 |

- Primality and repetition are not independent: 0.504·1229 = 619, and the count is 593.
- Two digit rules account for it.
  - Multiples of 3 (digit sum): 34.29% of the no-repeat strings, against 33.34% of all strings.
  - Multiples of 11 (alternating digit sum): exactly 1/9 of the no-repeat strings (560 of 5040), 1/27 of the strings with one repeated digit (160 of 4320), and 2/7 of those with two different digits (180 of 630).
  - With both rules the prediction is 597 (610.6 with the rule for 3 alone); the count is 593.
- Why 1/27: if the two equal digits sit at an odd distance, the alternating sum is the difference of the other two digits, which is never 0 mod 11.

**3. The exact structure of four slots.** 10⁴ − 1 = 9·11·101, and 10 is 1, −1 and a square root of −1 modulo the three factors (10² ≡ −1 mod 101).
- So a four-slot number mod 9999 is the Fourier transform of its digit string:
  - mod 9: the digit sum;
  - mod 11: the alternating sum;
  - mod 101: (d₀ − d₂) + i·(d₁ − d₃), with i = 10.
- Checked on all 10⁴ strings:
  - all arrangements of the same digits agree mod 9;
  - n + rev(n) ≡ 0 mod 11 and n ≡ rev(n) mod 9;
  - every palindrome abba is a multiple of 11, so there is no four-slot palindromic prime;
  - the multiples of 101 are exactly the strings abab;
  - rotating the digits multiplies by i mod 101.
- These are the untwisted sum, the half turn and the quarter turn of the cusp dictionary, here for the base ten.

**4. The same four different digits in every arrangement (210 digit sets, 24 arrangements each).**
- Number of prime arrangements → number of sets: 0 → 83, 1 → 8, 2 → 19, 3 → 20, 4 → 23, 5 → 13, 6 → 13, 7 → 14, 8 → 6, 9 → 6, 10 → 3, 11 → 2.
- 72 sets have digit sum divisible by 3, so no arrangement is prime. 70 sets lose 8 arrangements to 11 (two pairs of digits with equal sums mod 11).
- The maximum is 11 primes out of 24, for the digits 1, 2, 7, 9 and for 1, 2, 3, 7. From 1, 2, 7, 9: 1279, 1297, 2179, 2719, 2791, 2917, 2971, 7129, 7219, 9127, 9721.

**5. Sisters.** 102 pairs of four-digit primes are each other's reversal (204 primes), 1153 ↔ 3511 among them.
- Chance alone would give about 70 pairs. Reversal keeps the residue mod 9 and only flips the sign mod 11, which raises the expectation to about 115.
- The four-digit circular primes are the rotations of 1193 and 3779.

**6. The race mod 4.** Primes ≡ 3 mod 4 lead primes ≡ 1 mod 4 at every power of ten: by 2, 7, 10, 25, 147, 218, 446 for k = 2…8. The lead is of the size √x/log x (2, 5, 11, 27, 72, 196, 543). This is Chebyshev's bias, a statement about χ₋₄.

**What is happening.**
- The counts are the prime number theorem. The "half" in ¼, 1/6, 1/8 is 1/log 10 = 0.4343 seen at small size.
- Which strings are prime is governed by the three characters of four slots (9, 11, 101). Repeats, reversal and palindromes act on them in a simple exact way.

# Square stuff: the race mod 4, G between π²/12 and π²/8, and the Pythagorean identity (2026-10-02; square_stuff.py, square_stuff_out.txt). The user: "On that number 4. Try more stuff. Play with square stuffs"
**0. The race is a question about squares.** χ₋₄(p) = +1 exactly when −1 is a square mod p, which is exactly when p = a² + b². So "1 mod 4 against 3 mod 4" is "is −1 a square?".

**1. Why 3 mod 4 leads: the prime squares.** Every odd prime square is 1 mod 4.

| k | D(10^k) = π(x;4,3) − π(x;4,1) | ½·π(√x) | fair race (p^j counted 1/j) |
|---|---|---|---|
| 2 | 2 | 1.5 | 0.6 |
| 3 | 7 | 5.0 | 1.8 |
| 4 | 10 | 12.0 | −2.4 |
| 5 | 25 | 32.0 | −8.2 |
| 6 | 147 | 83.5 | 61.7 |
| 7 | 218 | 222.5 | −6.9 |
| 8 | 446 | 614.0 | −173.7 |

- The lead of 3 mod 4 is the half-count of prime squares. With prime powers counted as Riemann does, the race is level and oscillates around 0.
- Among the 5,761,454 odd primes below 10⁸, 1 mod 4 is ahead at only 1940 of them (first at p = 26861) and tied at 284.

**2. G is the race written in squares.** Over the odd primes:
    ∏ p²/(p² − 1) = π²/8,   ∏ p²/(p² − χ₋₄(p)) = G,   ∏ p²/(p² + 1) = π²/12.
- So G is the product in which each prime picks its sign by its class: p²/(p² − 1) if −1 is a square mod p, and p²/(p² + 1) if not. If every prime were 1 mod 4 the value would be π²/8 = 1.2337; if every prime were 3 mod 4 it would be π²/12 = 0.8225.
- G = (π²/8)·∏_{p≡3}(p² − 1)/(p² + 1) = (π²/12)·∏_{p≡1}(p² + 1)/(p² − 1). G is below 1 because 3 comes before 5.
- Partial products at x = 10, 10², 10³, 10⁴, 10⁶, 10⁸: 0.91875, 0.9160246, 0.9159650, 0.9159657, 0.9159655941, 0.9159655942.
- **The bias reaches G.** Let T(x) = log G − log(partial product up to x). For primes x ≤ 3·10⁶:
  - T(x)·x^{3/2}·log x has mean −0.307 and standard deviation 0.387;
  - the partial product is above G at 79.3% of those primes (73.8% with logarithmic weights).
  - Prediction of the mean: Σ_{p≤t}χ₋₄(p) ≈ −√t/log t from the prime squares, and T(x) = −S(x)/x² + 2∫_x^∞ S(t)t⁻³dt, which gives −⅓.
- Odd squares by class: Σ_{n≡1 (4)} 1/n² = π²/16 + G/2 and Σ_{n≡3 (4)} 1/n² = π²/16 − G/2. So G/2 is how far each class of odd squares sits from its fair share π²/16.

**3. The Pythagorean identity.**
- Let g(n) = r₂(n²)/4 = Σ_{d|n²}χ₋₄(d). It is multiplicative, with g(p^a) = 2a + 1 for p ≡ 1 (4) and 1 otherwise. So
    Σ g(n)·n⁻ˢ = ζ(s)²·L(s,χ₋₄)/((1 + 2⁻ˢ)·ζ(2s)).
- At s = 2 this is 2G. Since g(n) = 1 + 2·T(n), with T(n) the number of triples a > b > 0, a² + b² = n²:
    **G = π²/12 + Σ over all Pythagorean triples of 1/c²,   6G/π² = ½ + Σ over primitive triples of 1/c².**
- Second derivation: Σ over coprime pairs of 1/(m² + n²)² = 4ζ(2)G/ζ(4) = 60G/π². The both-odd pairs are a quarter of the opposite-parity ones, so the opposite-parity pairs give 48G/π². Remove the 4 axis points and divide by 8.
- Check: 3,183,100 primitive triples with c ≤ 2·10⁷ (N/2π = 3,183,099). Their sum plus the tail 1/(2πN) is 0.0568403091, and 6G/π² − ½ = 0.0568403091. For all triples: 0.0934985608 = G − π²/12.
- Reading in the user's terms: G/ζ(2) is a half plus a carry, and the carry is the sum of 1/c² over the primitive Pythagorean triples.
- Probably classical (the hypotenuse Dirichlet series goes back to Lehmer). Not checked.

**4. Sums of squares are this program's forms (checked by counting).**
- Two squares: r₂(n) = 4·(d₁(n) − d₃(n)), the race among the divisors of n. This is f = θ₃², the weight-1 form of case E.
- Six squares: r₆(n) = 16·Σ_{d|n}χ₋₄(n/d)d² − 4·Σ_{d|n}χ₋₄(d)d², i.e. 16·E_G − 4·E_ζ. G's own Eisenstein series counts sums of six squares.
- Σ over the square lattice of 1/(a² + b²)² = 4ζ(2)G = 6.02681204 (checked to 8 digits).

**5. Primes of the form n² + 1 (Landau's open problem) are the same race.** n² + 1 is never divisible by a prime that is 3 mod 4.
- The conjectured count is (C/2)·li(N), with C = ∏_{p odd}(1 − χ₋₄(p)/(p − 1)) = 1.3728.
- n ≤ 10³, 10⁴, 10⁵: 112, 841, 6656 primes; predicted 122, 855, 6610; without the factor C it would be 89, 623, 4815.

**6. Four slots: multiplying sisters extracts squares.** With x = d₀ − d₂, y = d₁ − d₃, and n* the string with d₁ and d₃ exchanged (true for all 10⁴ strings):
- n·n* ≡ x² + y² (mod 101);
- n·rev(n) ≡ −10·(x² + y²) (mod 101), ≡ −(alternating sum)² (mod 11), ≡ (digit sum)² (mod 9).
- Example: 1234·4321 ≡ 21 ≡ −80 (mod 101), with x = y = 2.

**What it gives.** These are identities and they locate G: between π²/12 and π²/8, a half of ζ(2) plus the Pythagorean sum. They do not give new rational approximations to G.

# Numbers written with ones only: the ruler behind combos, factorials and the sieve (2026-10-02; repunits.py, repunits_out.txt). The user: "I feel like more needs to be played with when it comes to putting solely 1s as digits. Like 111. … So like say you want 2345. It's bigger than 1111. So you need another digit. So 11111. So you have to compare 2345 and 11111 and 1111 in a way. Something foundational going on when just using 1s"
Notation: R_k = 11…1 (k ones) = (10^k − 1)/9, and R_k(b) = (b^k − 1)/(b − 1) in base b.

**1. Ones are the running totals of all combos.**
- R_k = 1 + 10 + ⋯ + 10^{k−1} is the number of digit strings with fewer than k slots (the empty one included).
- List all strings by length. String number n has k slots exactly when R_k ≤ n < R_{k+1}, and it is the string n − R_k. Checked on all 11111 strings with at most 4 slots.
- The user's example: 1111 ≤ 2345 < 11111, so 2345 is a 4-slot string, and 2345 − 1111 = 1234.
- The same statement is zero-free base ten (digits 1 to 10): n has k digits there exactly when R_k ≤ n < R_{k+1} (checked for n < 200000). For instance 1000 = (9, 9, 10), three digits, because 1000 < 1111.
- A power of ten always sits 8/9 of the way between two repunits: 10^18 − R₁₈ = 888888888888888889.

**2. Ones are the weights of factorials (Legendre).**
- Read the digits of n with the weights 1, 11, 111, … in place of 1, 10, 100, …. The result is Σ d_j·R_j = (n − digit sum)/9 = Σ_{j≥1}⌊n/10^j⌋.
- In a prime base p this is the exponent of p in n!: v_p(n!) = Σ d_j·R_j(p). Checked for p ≤ 13, n ≤ 5000.
  - v₅(625!) = 156, which is 1111 in base 5: the factorial of a power of the base has a repunit exponent.
  - 2345 = (3, 3, 3, 4, 0) in base 5, so 2345! ends in 3·156 + 3·31 + 3·6 + 4·1 = 583 zeros.
- The other ruler is lcm(1..n): its exponent of p is (number of base-p digits of n) − 1.
- So the two denominators of this program use the two rulers: lcm(1..n) counts digits against 1, 10, 100, …, and n! weighs digits against 1, 11, 111, ….
- In repunits, 2345 = 2·1111 + 1·111 + 1·11 + 1·1, and those digits form 21110, the number with Σ⌊21110/10^j⌋ = 2345.

**3. The sieve written in ones.**
- 1/9 + 1/99 + 1/999 + ⋯ = 0.1223242434262445262644283446282644492448282664…. The n-th decimal is d(n), the number of divisors of n, because 1/(10^k − 1) puts a 1 at every multiple of k.
  - So the primes are exactly the places holding a 2: 2, 3, 5, 7, 11, …, 43.
  - It holds up to n = 46. d(48) = 10 is the first carry, and it spoils the digit at 47.
- With signs on the odd lengths, 1/9 − 1/999 + 1/99999 − ⋯ = 0.110120011200200121020000320020010201…. The n-th decimal is r₂(n)/4, the count of ways to write n as a sum of two squares (60 digits checked).
- θ₃(1/10) = 1.2002000020000002000000002…, with its 2s at the squares, and θ₃(1/10)² = 1 + 4·(the previous number). This is the weight-1 form f of case E at q = 1/10, readable digit by digit.

**4. Squares and powers of ones.**
- R_k² is the palindrome 123…k…321 for k ≤ 9 (R₉² = 12345678987654321). At ten ones a carry breaks it: R₁₀² = 1234567900987654321.
- 11ⁿ is the n-th row of Pascal's triangle for n ≤ 4 (11⁴ = 14641). 11⁵ = 161051 against the row 1, 5, 10, 10, 5, 1: carries again.
- R_k ≡ 3 mod 4 for every k ≥ 2. So a number made of ones is never a square and never a sum of two squares, and it always has an odd number of prime factors that are 3 mod 4 (checked to k = 18).

**5. Infinitely many ones.** R_k(b) ≡ −1/(b − 1) mod b^k, so in the b-adic numbers …111 = −1/(b − 1): Euler's value of 1 + b + b² + ⋯.
- Base ten: −1/9. Base 3: −½. Base 7: −1/6. Base 13: −1/12.
- These equal ζ(0), −B₂ and ζ(−1) as numbers. The series are different, so this is a match of values, not a derivation.

**6. Ones and primes.** For p ≠ 2, 3, 5 the first repunit divisible by p has as many ones as the period of 1/p (checked for p < 1000).

| ones | new primes | ones | new primes |
|---|---|---|---|
| 2 | 11 | 11 | 21649, 513239 |
| 3 | 3, 37 | 12 | 9901 |
| 4 | 101 | 13 | 53, 79, 265371653 |
| 5 | 41, 271 | 14 | 909091 |
| 6 | 7, 13 | 15 | 31, 2906161 |
| 7 | 239, 4649 | 16 | 17, 5882353 |
| 8 | 73, 137 | 17 | 2071723, 5363222357 |
| 9 | 333667 | 18 | 19, 52579 |
| 10 | 9091 | | |

- 487 is a Wieferich prime in base ten: the first repunit it divides is R₄₈₆, and 487² divides it too.
- R_k is prime for k = 2, 19, 23 (k ≤ 60). Whether there are infinitely many repunit primes is open.

**7. The no-repeat totals.** Strings with at most k slots and no repeated digit: 1, 11, 101, 821, 5861, 36101, …, ending at 9864101 = ⌊e·10!⌋.
- With repeats the totals are the partial sums of a geometric series (repunits). Without repeats they are the partial sums of the series for e, times 10!.
- With at most 4 slots: 5861 of 11111, a share of 0.5275.

**What is foundational here.** The number with k ones is the count of everything shorter than k slots. That one fact gives the three rulers above: combos (totals), factorials (Legendre weights) and primes (periods and divisor counts).

