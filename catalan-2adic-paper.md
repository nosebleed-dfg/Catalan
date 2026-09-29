# The Price of the Prime 2 in Catalan's Constant

**An exact formula for the powers of 2 in the Hankel-determinant approximations to Catalan's constant, with a complete proof written for non-specialists**

nos3bl33d · [github.com/nosebleed-dfg/Catalan](https://github.com/nosebleed-dfg/Catalan) · 29 September 2026

*Worked out with AI assistance (Claude, by Anthropic). Every step is checked by a script in [`proof_checks/`](proof_checks). This has not been peer-reviewed; corrections are welcome as GitHub issues.*

> **Read this first.** This is **not** a proof that Catalan's constant is irrational. Nobody knows whether it is; that is a famous open problem. This paper answers a smaller, exact question that every attempt of a certain kind runs into: *how many factors of 2 pile up in the denominators?* Section 12 explains what this does and does not mean for the irrationality problem.

---

## Abstract

Catalan's constant G = 1 − 1/3² + 1/5² − 1/7² + ⋯ is believed to be irrational, but no proof is known.

A standard way to attack such questions is Apéry's method. You build rational approximations and weigh how fast they converge against how fast their denominators grow. We study one natural family of such approximations. For each size K there is a K×K determinant D_K(X) = det(A_K + X·B_K), a polynomial in X with rational coefficients. At X = G it equals a positive Gram determinant.

**Main result.** We prove an exact formula for e₂(K), the largest power of 2 occurring in the denominators of D_K. It is minus the sum of the K smallest numbers in two explicit lists of integers, with one extra factor of 2 when there is a tie.
- One list is built from binary digit sums (the "band").
- The other is built from factorials (the "moments").

**Method.** The proof shifts X by the *2-adic Catalan constant*

  C = ½ + 1 − 4 + 48 − 1088 + 39680 − ⋯

(a series of tangent numbers that diverges as ordinary numbers but converges 2-adically). After the shift the matrix splits into a "coin-flip" (Krawtchouk) matrix and a classical moment matrix. Pascal's triangle mod 2 glues the two together. This constant C is the 2-adic L-value that Calegari proved irrational in 2005.

**Consequence.** e₂(K) = K²/6 + (K/3)·log₂K + O(K). This accounts exactly for the ⅙·ln 2 per K² that the prime 2 costs this approach.

---

## The result at a glance

Write **s₂(n)** for the number of 1s in the binary expansion of n; for example s₂(6) = s₂(110₂) = 2. For each size K, form two lists:

| list | formula | the list for K = 8 |
|---|---|---|
| **band** | b_m = 4m − 2K + s₂(K−1−m) − s₂(m), for m = 0, 1, …, K−1 | −13, −11, −7, −5, 1, 3, 7, 9 |
| **moments** | h_q = 8q − 3·s₂(q) − s₂(q+1), for q = 0, 1, 2, … | −1, 4, 11, 17, 27, 32, 39, 46, … |

> **Main Theorem.** The largest power of 2 that divides a denominator of a coefficient of D_K(X) is 2^e₂(K), where
>
> **e₂(K) = −(the sum of the K smallest numbers among b₀, …, b_{K−1}, h₀, …, h_{K−1})**,
>
> minus 1 more if the K-th and (K+1)-th smallest of these 2K numbers are equal.

**Example, K = 8.**
- The 8 smallest numbers are −13, −11, −7, −5, −1, 1, 3, 4: six from the band and two from the moments.
- Their sum is −29. The 9th smallest number is 7, which differs from the 8th (4), so there is no tie.
- Hence e₂(8) = 29. Some coefficient of D₈ has exactly 2²⁹ in its denominator, and none has more.

**Example, K = 12 (a tie).**
- Band: −21, −19, −15, −13, −6, −4, 0, 2, 9, 11, 15, 17. Moments: −1, 4, 11, 17, …
- The 12 smallest sum to −53.
- The 12th and 13th smallest are both 11, one from each list. So e₂(12) = 53 − 1 = 52.

| K | e₂(K), computed from the determinant | e₂(K), from the formula | K²/6 |
|---:|---:|---:|---:|
| 3 | 8 | 8 | 1.5 |
| 4 | 11 | 11 | 2.7 |
| 5 | 16 | 16 | 4.2 |
| 6 | 18 | 18 | 6.0 |
| 7 | 24 | 24 | 8.2 |
| 8 | 29 | 29 | 10.7 |
| 9 | 36 | 36 | 13.5 |
| 10 | 43 | 43 | 16.7 |
| 11 | 48 | 48 | 20.2 |
| 12 | 52 | 52 | 24.0 |
| 16 | 82 | 82 | 42.7 |
| 32 | 258 | 258 | 170.7 |
| 64 | 874 | 874 | 682.7 |
| 128 | 3146 | 3146 | 2730.7 |
| 256 | 11835 | 11835 | 10922.7 |
| 384 | 26113 | 26113 | 24576.0 |

The middle column was computed directly from D_K(X) in exact rational arithmetic, for 134 values of K: every K from 3 to 128, and K = 160, 180, 192, 224, 256, 320, 360, 384. The formula agrees every time. Here is the whole theorem as a program:

```python
def s2(n):                      # number of 1s in the binary expansion of n
    return bin(n).count("1")

def e2(K):
    band    = [4*m - 2*K + s2(K - 1 - m) - s2(m) for m in range(K)]
    moments = [8*q - 3*s2(q) - s2(q + 1)          for q in range(K)]
    both = sorted(band + moments)
    tie = both[K - 1] == both[K]          # K-th and (K+1)-th smallest equal?
    return -sum(both[:K]) - (1 if tie else 0)

print([e2(K) for K in range(3, 13)])      # [8, 11, 16, 18, 24, 29, 36, 43, 48, 52]
```

**Why K²/6.**
- The band climbs by about 4 per step, from about −2K up to about +2K. The moment list climbs by about 8 per step, starting near 0.
- Taking the K smallest overall means using M moments and K − M band values, where the next band value (about 2K − 4M) meets the next moment value (about 8M).
- That happens at M ≈ K/6, and the total works out to about −K²/6.
- The binary digit sums add the (K/3)·log₂K term.

**How to read this paper.**
- In 5 minutes: this page.
- In 30 minutes: add §1–3.
- For the complete proof: §4–10.
- For what it means for irrationality: §12.

---

## Contents

1. Background
2. The results
3. The idea of the proof, in one page
4. The 2-adic Catalan constant
5. Piece one: B_K is a coin-flip matrix
6. Piece two: shifted by C, the matrix is a moment matrix
7. The glue: Pascal's triangle mod 2
8. The cheapest term wins: the coefficient law
9. Back to X: the content and the tie rule
10. The fingerprint of A_K
11. Growth: where K²/6 and (K/3)·log₂K come from
12. What this means for the irrationality of G

Appendices:
- A. Where the machine's numbers come from
- B. The classical orthogonal polynomials
- C. Check it yourself
- D. Glossary

References.

---

## 1. Background

### 1.1 Catalan's constant

  G = Σ_{k≥0} (−1)^k/(2k+1)² = 1 − 1/9 + 1/25 − 1/49 + ⋯ = 0.9159655941772190150546…

It turns up in many places:
- ∫₀¹ arctan(x)/x dx = G.
- The number of domino tilings of an m×n rectangle grows like e^{(G/π)·mn} (Kasteleyn; Temperley–Fisher, 1961).
- The hyperbolic volume of the Whitehead link complement is 4G.

It is the simplest constant of its kind whose irrationality is unknown. Compare:
- ζ(2) = π²/6 has long been known to be irrational;
- ζ(3) was proved irrational by Apéry in 1978;
- the "sister" constant L(2, χ₋₃) = 1 − 1/2² + 1/4² − 1/5² + ⋯ (period 3) was proved irrational only in 2024 (Calegari, Dimitrov and Tang).

### 1.2 How irrationality proofs work (Apéry's recipe)

Suppose G = p/q. Then for any integers a and b, the number a·G − b = (ap − bq)/q is either 0 or at least 1/q in size. So if you can find integers a_n, b_n with a_nG − b_n ≠ 0 and a_nG − b_n → 0, then G cannot be rational.

In practice the natural approximations have rational coefficients, and you must multiply through by their denominators. Everything comes down to a contest:
- **gain:** how fast the approximation error shrinks;
- **cost:** how fast the denominators grow.

The cost splits prime by prime: for each prime p, how many factors of p appear. Apéry won this contest for ζ(3). For G, every construction tried so far has lost it. This paper computes the prime-2 part of the cost exactly, for one natural family.

### 1.3 The machine

**The weight and its averages.** Let ω(t) = π·t²·cosh(πt)/sinh²(πt) for t > 0. For a function f define the average

  φ(f) = ∫₀^∞ f(4t²)·ω(t) dt.

We write v for the variable 4t². Two facts, both proved in Appendix A:

- **Moments.** μ_k := φ(v^k) = 4^k·(2^{2k+2} − 1)·|B_{2k+2}|, where B_n are the Bernoulli numbers. So μ₀ = 1/2, μ₁ = 2, μ₂ = 24, μ₃ = 544, … In terms of the tangent numbers 1, 2, 16, 272, 7936, … (tan x = x + 2x³/3! + 16x⁵/5! + ⋯), μ_k = (k+1)·T_{2k+1}/2.
- **Pole values.** For j ≥ 0, φ(1/(v + (2j+1)²)) = F_j(G), where
  - F_j(X) := (−1)^j·(2j+1)/2·(X − β_j) − 1/(4(2j+1));
  - β_j := Σ_{k<j} (−1)^k/(2k+1)² are the partial sums of G's series.
  - For example F₀(X) = X/2 − 1/4 and F₁(X) = −(3/2)(X − 1) − 1/12.

**The matrix.** Fix K. Put w_j := −(2j+1)² and N_K(v) := (v − w₀)(v − w₁)⋯(v − w_{K−1}) = (v + 1²)(v + 3²)⋯(v + (2K−1)²). For 0 ≤ i, k < K, the (i, k) entry of the K×K matrix M_K(X) is obtained as follows:
1. Split v^{i+k}/N_K(v) by partial fractions into a polynomial P(v) plus Σ_j r_j/(v − w_j).
2. Replace φ(P) by its value, computed from the moments.
3. Replace each φ(1/(v − w_j)) by the linear function F_j(X).

The result is M_K(X) = A_K + X·B_K with rational matrices A_K and B_K. At X = G every replacement is the true integral, so M_K(G) = [∫₀^∞ (4t²)^{i+k}·ω(t)/N_K(4t²) dt]. This is the Gram matrix of a positive weight, so its determinant is positive.

**The object of study.** D_K(X) := det(A_K + X·B_K), a polynomial of degree K in X with rational coefficients d₀, …, d_K. Its **2-adic content exponent** is

  e₂(K) := max over j of (the power of 2 in the reduced denominator of d_j).

### 1.4 2-adic numbers in five minutes

**Valuation.** For a nonzero rational x, v₂(x) is the number of factors of 2 in its numerator minus the number in its denominator. For example v₂(12) = 2, v₂(3/8) = −3, v₂(5) = 0. So e₂(K) = −min_j v₂(d_j). "x is 2-adically small" means that v₂(x) is large.

**The two rules.**
- v₂(xy) = v₂(x) + v₂(y).
- v₂(x + y) ≥ min(v₂(x), v₂(y)), **with equality when v₂(x) ≠ v₂(y)**.

The second rule is the heart of every argument below. A sum is exactly as divisible by 2 as its least divisible term, *unless* two terms tie for least divisible, in which case they may cancel.

**The 2-adic numbers ℚ₂.** Complete the rationals with respect to this notion of size, the way the reals complete them with respect to the usual one. In ℚ₂:
- a series Σa_n converges exactly when v₂(a_n) → ∞;
- for example, 1 + 2 + 4 + 8 + ⋯ = −1, because (1 − 2)(1 + 2 + ⋯ + 2^{n−1}) = 1 − 2ⁿ → 1;
- the numbers with v₂ ≥ 0 are the 2-adic integers ℤ₂, which include every fraction with an odd denominator.

**Two classical facts about binary digits.**
- Legendre: v₂(n!) = n − s₂(n).
- Kummer: v₂ of the binomial coefficient C(a+b, a) equals the number of carries when a and b are added in binary.

**The fingerprint (Smith form).** For an invertible matrix M over ℚ₂ there are integers d₁ ≤ d₂ ≤ ⋯ ≤ d_K such that, for each k, the smallest v₂ among all k×k minors of M is d₁ + ⋯ + d_k. Call (d₁, …, d_K) the *fingerprint* of M.
- It does not change when M is multiplied on either side by a matrix with 2-adic-integer entries and odd determinant.
- v₂(det M) = d₁ + ⋯ + d_K.

**Cauchy–Binet.** For a K×n matrix X and an n×K matrix Y, det(XY) = Σ_S det X[:, S]·det Y[S, :], summing over K-element subsets S of {1, …, n}.

---

## 2. The results

**Definition (the two lists).** For K ≥ 1 put

  b_m := 4m − 2K + s₂(K−1−m) − s₂(m)   (0 ≤ m < K),
  h_q := 8q − 3·s₂(q) − s₂(q+1)     (q ≥ 0).

Equivalently, by Legendre's formula, h_q = 4q + 3·v₂(q!) + v₂((q+1)!) − 1.

**Lemma 2.1.** Both lists are strictly increasing: b_{m+1} − b_m ≥ 2 and h_{q+1} − h_q ≥ 4.

*Proof.* Adding 1 to n changes s₂ by 1 − v₂(n+1), which is at most 1. Hence
- b_{m+1} − b_m = 4 − [s₂(K−1−m) − s₂(K−2−m)] − [s₂(m+1) − s₂(m)] ≥ 4 − 1 − 1 = 2;
- h_{q+1} − h_q = 8 − 3(1 − v₂(q+1)) − (1 − v₂(q+2)) = 4 + 3v₂(q+1) + v₂(q+2) ≥ 4. ∎

Let L_K be the K smallest entries of the combined list b₀, …, b_{K−1}, h₀, …, h_{K−1} (2K numbers). Let τ_K = 1 if the K-th and (K+1)-th smallest entries are equal, and τ_K = 0 otherwise.

**Theorem 1 (the power of 2 in D_K).** e₂(K) = −(sum of L_K) − τ_K.

**Theorem 2 (the fingerprint).** The 2-adic fingerprint of A_K = M_K(0) is L_K, in increasing order, except that when τ_K = 1 its largest entry is strictly larger.

**Theorem 3 (growth).** e₂(K) = K²/6 + (K/3)·log₂K + O(K).

For example, the fingerprint of A₈ is (−13, −11, −7, −5, −1, 1, 3, 4), which is literally L₈. For K = 12 it is (−21, −19, −15, −13, −6, −4, −1, 0, 2, 4, 9, 12): the top entry is 12 rather than 11 because of the tie.

---

## 3. The idea of the proof, in one page

1. **Shift X by a 2-adic Catalan constant.** There is a 2-adic number C with v₂(C) = −1 (§4) such that A_K + X·B_K = (A_K + C·B_K) + (X − C)·B_K, and the two matrices A_K + C·B_K and B_K are each *simple*.
2. **B_K is a coin-flip matrix (§5).** It is the moment matrix of the binomial distribution: the number of heads in 2K − 1 fair coin tosses, recentred and squared. Its natural polynomials are the Krawtchouk polynomials. In their basis B_K is diagonal, and the powers of 2 on the diagonal are exactly the band list (plus one). The binary digit sums in b_m are Kummer's carries in the coin-flip normalisations.
3. **A_K + C·B_K is a moment matrix (§6).** At X = C the pole values become the 2-adic expansion of the same averages, so A_K + C·B_K is the moment matrix of the weight ω, multiplied by the harmless factor 1/N_K. In the basis of ω's own orthogonal polynomials it is (up to a triangular change of basis with odd diagonal) diagonal, with powers of 2 equal to the moment list.
4. **The glue (§7).** The two bases differ by a triangular matrix which, mod 2, is Pascal's triangle, i.e. Sierpinski's triangle. Every "corner" determinant of Pascal's triangle equals 1.
5. **Count (§8).** Expand det((A_K + C·B_K) + Y·B_K) with Cauchy–Binet. For each power of Y, the terms correspond to choosing some directions from each list. The term that picks the smallest numbers is the only cheapest term, and its coefficient is a Pascal corner, which is odd, so nothing can cancel it.
6. **Shift back (§9).** Undo X = Y + C. Because C has exactly one factor 1/2, the bookkeeping turns "sum of the K smallest numbers" into e₂(K). A tie produces exactly one extra factor of 2, by a parity argument.

---

## 4. The 2-adic Catalan constant

**Euler numbers at 0.** The values E_n(0) of the Euler polynomials are defined by

  2/(e^t + 1) = Σ_{n≥0} E_n(0)·tⁿ/n!.

Since 2/(e^t + 1) = 1 − tanh(t/2), and tanh x = Σ_{n≥1} 2^{2n}(2^{2n} − 1)B_{2n}·x^{2n−1}/(2n)!, we get:
- E₀(0) = 1;
- E_n(0) = 0 for even n ≥ 2;
- E_{2k+1}(0) = −(2^{2k+2} − 1)·B_{2k+2}/(k+1).

**The coefficients a_n.** Put a_n := E_n(0)·(−2)ⁿ·(n+1). Then a₀ = 1, a_n = 0 for even n ≥ 2, and, using B_{2k+2} = (−1)^k|B_{2k+2}|,

  a_{2k+1} = 2^{2k+2}·(2^{2k+2} − 1)·B_{2k+2} = 4·(−1)^k·μ_k.

By von Staudt–Clausen, each Bernoulli number B_{2n} has exactly one factor 2 in its denominator. So v₂(μ_k) = 2k − 1 and v₂(a_{2k+1}) = 2k + 1.

**Definition 4.1.** For a 2-adic number y with v₂(y) ≥ 0 put

  ρ₁(y) := −½ Σ_{n≥0} a_n·(2y+1)^{−n−2} = −½(2y+1)^{−2} − 2·Σ_{k≥0} (−1)^k·μ_k·(2y+1)^{−2k−3},

  **C := −ρ₁(0) = ½ + 2·Σ_{k≥0} (−1)^k μ_k = ½ + 1 − 4 + 48 − 1088 + 39680 − ⋯.**

Both series converge 2-adically: 2y + 1 is odd, and v₂(μ_k) = 2k − 1 → ∞. Every term of C after ½ is an integer, so C = ½ + (an odd integer) + (a multiple of 4), and therefore **v₂(C) = −1**. The 2-adic digits of 2C end in …0100010111011011.

**Lemma 4.2 (Boole's formula).** ρ₁(y) + ρ₁(y+1) = −g(y), where g(y) := (2y+1)^{−2}.

*Proof.* The derivatives of g are g^{(n)}(y)/n! = (−2)ⁿ(n+1)·(2y+1)^{−n−2}, so ρ₁ = −½ Σ_n E_n(0)·g^{(n)}/n!. By Taylor's theorem ρ₁(y+1) = Σ_{m≥0} ρ₁^{(m)}(y)/m!, and hence

  ρ₁(y) + ρ₁(y+1) = −½ Σ_N g^{(N)}(y)·[E_N(0)/N! + Σ_{n+m=N} E_n(0)/(n!·m!)].

The bracket is the coefficient of t^N in (1 + e^t)·Σ_n E_n(0)tⁿ/n! = 2, so it equals 2 for N = 0 and 0 otherwise. ∎

*Why the formal steps are legitimate.* For |y′ − y| < 2 in the 2-adic sense, 2y′ + 1 is a 2-adic unit. So every term of ρ₁ is a convergent power series on that disk, and so is ρ₁ itself, because the terms shrink uniformly. The point y + 1 lies in the disk, so the Taylor expansion converges there. The double series converges absolutely: |E_n(0)|₂ ≤ 2(n+1) and |g^{(N)}(y)/N!|₂ ≤ 2^{−N}. So the rearrangement is allowed.

**Corollary 4.3 (the splitting of the partial sums).** For every integer y ≥ 0,

  β_y = C + (−1)^y·ρ₁(y).

*Proof.* Put ρ(y) := (−1)^y β_y. Since β_{y+1} = β_y + (−1)^y g(y), ρ satisfies the same relation ρ(y) + ρ(y+1) = −g(y). The difference ρ − ρ₁ therefore changes sign at each step, so ρ(y) − ρ₁(y) = (−1)^y·(ρ(0) − ρ₁(0)) = (−1)^y·(0 + C). ∎

*In plain words.* In the real numbers the partial sums β_y converge to G. In the 2-adic numbers they do not converge, but they split into a constant C plus an alternating part that varies smoothly. C plays the role of G in the 2-adic world. Indeed β_{2^N} → 0 2-adically, because ρ₁(2^N) → ρ₁(0) = −C.

**Corollary 4.4 (the 2-adic Stieltjes identity).** For every j ≥ 0,

  F_j(C) = Σ_{k≥0} (−1)^k·μ_k/(2j+1)^{2k+2}.

*Proof.*
- By Corollary 4.3, C − β_j = −(−1)^j ρ₁(j), so F_j(C) = −(2j+1)/4·(2ρ₁(j) + g(j)).
- In 2ρ₁(j) + g(j), the n = 0 term of ρ₁ cancels g(j), leaving 2ρ₁(j) + g(j) = −Σ_{n≥1} a_n(2j+1)^{−n−2} = −4·Σ_k (−1)^kμ_k(2j+1)^{−2k−3}.
- Hence F_j(C) = (2j+1)·Σ_k (−1)^kμ_k(2j+1)^{−2k−3} = Σ_k (−1)^kμ_k/(2j+1)^{2k+2}. ∎

*In plain words.* The real number F_j(G) is the integral φ(1/(v + (2j+1)²)). Expanding 1/(v + (2j+1)²) = Σ_k (−v)^k/(2j+1)^{2k+2} and integrating term by term gives the series Σ(−1)^kμ_k/(2j+1)^{2k+2}. In the real numbers that series diverges badly: it is only an asymptotic series, because μ_k grows like a factorial. **In the 2-adic numbers the same series converges, and its value is the same formula F_j with G replaced by C.** The pole values at X = C are "the 2-adic integral".

**Credit.** C is a known constant. It agrees with the 2-adic value L₂(2, χ₋₄) studied by Calegari (2005) to 2-adic order 2⁻³⁵: his Apéry-type approximation 783269/13060350 is guaranteed correct to order 2⁻³⁴, and differs from our C by a number divisible by 2³⁵. **Calegari proved that it is irrational**, and Beukers (2008) gave a proof by Stieltjes continued fractions, the classical cousins of the orthogonal polynomials in §6. What is new here is how C enters the Hankel machine, and the consequences below.

---

## 5. Piece one: B_K is a coin-flip matrix

### 5.1 B_K is the moment matrix of a discrete weight

Only the pole terms of the partial-fraction split contain X, and the X-coefficient of F_j is γ_j := (−1)^j(2j+1)/2. Therefore

  (B_K)_{ik} = Σ_{j<K} θ_j·w_j^{i+k},  with θ_j := γ_j/N_K′(w_j).

Now N_K′(w_j) = Π_{i≠j}((2i+1)² − (2j+1)²) = Π_{i≠j} 4(i−j)(i+j+1). Using Π_{i≠j}(i − j) = (−1)^j·j!·(K−1−j)! and Π_{i≠j}(i+j+1) = (K+j)!/(j!·(2j+1)), this equals

  N_K′(w_j) = 4^{K−1}·(−1)^j·(K−1−j)!·(K+j)!/(2j+1),

and hence

  θ_j = (2j+1)²/(2·4^{K−1}·(K−1−j)!·(K+j)!) = C(2K−1, K−1−j)·(2j+1)²/(2·4^{K−1}·(2K−1)!) > 0.

**In plain words.** Toss 2K − 1 fair coins, subtract (2K−1)/2 from the number of heads to get y ∈ {±½, ±3/2, …}, and record v = −4y². The weights θ_j are those probabilities (times y² and a constant). B_K is the table of averages of v^{i+k} under that distribution. Precisely, with N := 2K − 1 and x ∈ {0, …, N}, y = x − N/2:

  Σ_j θ_j·f(w_j) = (1/(4^{K−1}·N!))·Σ_{x=0}^{N} C(N, x)·y²·f(−4y²),

because each pair ±y lands on the same node.

### 5.2 The coin-flip (Krawtchouk) polynomials

**The generating function.** Define polynomials K_n(x) by

  Σ_{n=0}^{N} C(N, n)·K_n(x)·tⁿ = (1 − t)^x·(1 + t)^{N−x}.

**Orthogonality.** Multiplying two such generating functions and summing against C(N, x):

  Σ_x C(N,x)·[(1−t)(1−s)]^x·[(1+t)(1+s)]^{N−x} = [(1−t)(1−s) + (1+t)(1+s)]^N = 2^N·(1 + ts)^N.

Comparing coefficients of t^m s^n gives Σ_x C(N,x)·K_m(x)·K_n(x) = 0 for m ≠ n, and = 2^N/C(N,n) for m = n.

**The monic version.** K_n has degree n and leading coefficient (−2)ⁿ/(n!·C(N,n)), read off from the x-degree-n part of (1−t)^x(1+t)^{N−x}. Let Q_n be the monic multiple of K_n, written in y = x − N/2. Then

  ‖Q_n‖² := Σ_x C(N,x)·Q_n(y)² = 2^N·n!·N!/((N−n)!·4ⁿ).

**The recurrence.** The ratios λ_n := ‖Q_n‖²/‖Q_{n−1}‖² = n(N+1−n)/4 are the coefficients of the three-term recurrence y·Q_n = Q_{n+1} + λ_n·Q_{n−1}. There is no Q_n term, because the coin-flip distribution is symmetric under y → −y.

### 5.3 From y to v: integral polynomials

Our weight is y²·C(N,x), and it only sees y². By symmetry the odd polynomials have the form Q_{2n+1}(y) = y·R_n(y²), and the R_n are exactly the monic orthogonal polynomials of the weight y²·C(N,x) in the variable y². Rescale to v = −4y²:

  p_n(v) := (−4)ⁿ·R_n(−v/4).

These are monic in v and orthogonal for B_K. From s·R_n = R_{n+1} + (λ_{2n+1} + λ_{2n+2})R_n + λ_{2n}λ_{2n+1}R_{n−1} (with s = y²) one gets

  v·p_n = p_{n+1} + α′_n·p_n + κ′_n·p_{n−1}, with
  α′_n = −[(2n+1)(2K−1−2n) + 4(n+1)(K−1−n)],  κ′_n = 4n(K−n)(2n+1)(2K−1−2n).

Two consequences:
- **Integrality.** The coefficients are integers, so every p_n has integer coefficients. The change of basis from 1, v, v², … to p₀, p₁, p₂, … is unitriangular with integer entries, so it preserves fingerprints. In this basis B_K = diag(θ̂₀, …, θ̂_{K−1}) with θ̂_n := Σ_j θ_j·p_n(w_j)².
- **Mod 2.** α′_n is odd and κ′_n is divisible by 4, so p_{n+1} ≡ (v + 1)·p_n (mod 2). Hence p_n ≡ (v + 1)ⁿ (mod 2).

### 5.4 The norms: the band list appears

From §5.1–5.3, θ̂_n = 16ⁿ·‖Q_{2n+1}‖²/(4^{K−1}·(2K−1)!). With Legendre's formula and s₂(2m+1) = 1 + s₂(m):

  v₂(θ̂_n) = 4n + [N + v₂((2n+1)!) + v₂(N!) − v₂((N−2n−1)!) − 2(2n+1)] − [2K − 2 + v₂(N!)]
     = 4n − 2K + 1 + s₂(K−1−n) − s₂(n) = **b_n + 1**.

*In plain words.* The band list is the list of powers of 2 in the normalisations of the coin-flip polynomials. The digit sums s₂(K−1−n) − s₂(n) are the binary carries that Kummer's theorem counts in binomial coefficients.

---

## 6. Piece two: shifted by C, the matrix is a moment matrix

### 6.1 The weight's own orthogonal polynomials (a classical fact)

Let P₀, P₁, … be the monic orthogonal polynomials of the averages φ, in the variable v. They satisfy

  v·P_n = P_{n+1} + α_n·P_n + κ_n·P_{n−1},  α_n = 4(2n+1)(n+1),  κ_n = 16n³(n+1),

with norms H_n := φ(P_n²) = 16ⁿ·(n!)³·(n+1)!/2.

These are the *continuous dual Hahn polynomials* with parameters (1, 1, 0), rescaled (Wilson 1980; Koekoek–Lesky–Swarttouw 2010, §9.3). Appendix B shows the specialisation. The recurrence was also checked against the exact moments μ_k for all n < 30.

By Legendre's formula,

  v₂(H_q) = 4q + 3v₂(q!) + v₂((q+1)!) − 1 = **h_q**.

This is the moment list. Also P_n ≡ vⁿ (mod 4), because α_n ≡ 0 (mod 4) and κ_n ≡ 0 (mod 16).

### 6.2 Multiplication by v, rescaled

Let 𝒥 be the infinite matrix with 𝒥_{n+1,n} = 1, 𝒥_{n,n} = α_n and 𝒥_{n−1,n} = κ_n, so that v·P_q = Σ_p 𝒥_{pq}P_p. By induction, f(v)·P_q = Σ_p f(𝒥)_{pq}·P_p for every polynomial f, and by orthogonality

  φ(f·P_p·P_q) = H_p·f(𝒥)_{pq}.

The entries of 𝒥 are not 2-adically small, because the 1s below the diagonal are not. **Rescale by a rational diagonal matrix**, 𝒥̂ := S^{−1}𝒥S with S = diag(1, 4^{−1}, 4^{−2}, …). The entries of 𝒥̂ are

  𝒥̂_{n+1,n} = 4,  𝒥̂_{n,n} = 4(2n+1)(n+1),  𝒥̂_{n−1,n} = κ_n/4 = 4n³(n+1),

all divisible by 4. Each entry of 𝒥̂^k is a finite sum of products of k entries, so it is divisible by 4^k. Undoing the rescaling gives (𝒥^k)_{pq} = 4^{q−p}·(𝒥̂^k)_{pq}. For example μ_k = H₀·(𝒥^k)₀₀ = ½·(𝒥̂^k)₀₀, which gives v₂(μ_k) ≥ 2k − 1 once more.

### 6.3 A_K + C·B_K is the moment matrix of φ times 1/N_K

For a power series f = Σ f_k v^k with 2-adically bounded coefficients, put φ⁽²⁾(f) := Σ_k f_k·μ_k. This converges because μ_k → 0 2-adically, and it agrees with φ on polynomials.

Corollary 4.4 says exactly that F_j(C) = φ⁽²⁾(1/(v − w_j)), since 1/(v − w_j) = Σ_k (−1)^k v^k/(2j+1)^{2k+2}. By the partial-fraction definition of the matrix,

  (A_K + C·B_K)_{ik} = φ⁽²⁾(v^{i+k}·π(v)),  π(v) := 1/N_K(v) = Σ_k c_k v^k.

Here c_k ∈ ℤ₂, and c₀ = 1/(1²·3²⋯(2K−1)²) is odd.

### 6.4 The moment side is clean

**Proposition 6.1.** Let G_{pq} := φ⁽²⁾(P_p·P_q·π) for p, q < K. This is the matrix A_K + C·B_K written in the basis P₀, …, P_{K−1}. Then G = LΛLᵀ, where:
- L is lower triangular with 1s on the diagonal and entries in ℤ₂, with v₂(L_{pq}) ≥ 3(p − q) + 2 for p > q (in particular L ≡ I mod 32);
- Λ = diag(H₀u₀, …, H_{K−1}u_{K−1}) with every u_q odd. So v₂(Λ_q) = h_q.

*Proof.*
1. **Expand in the rescaled matrix.** Rearranging the absolutely convergent double series and using §6.2:

   G_{pq} = Σ_k c_k·φ(v^k P_p P_q) = H_p·Σ_k c_k(𝒥^k)_{pq} = H_p·4^{q−p}·Π̃_{pq}, with Π̃ := Σ_k c_k·𝒥̂^k.

   Every entry of Π̃ is a 2-adic integer, and Π̃ ≡ c₀·I (mod 4).
2. **Eliminate in ℤ₂.** Perform Gaussian elimination on the K×K block Π̃_{[K]}. Each pivot is ≡ c₀ (odd), and each multiplier is ≡ 0 (mod 4), so the remaining block stays ≡ c₀·I (mod 4) throughout. Result: Π̃_{[K]} = 𝐋𝐔, with 𝐋 lower unitriangular ≡ I and 𝐔 upper triangular ≡ c₀·I (mod 4), both with 2-adic-integer entries.
3. **Undo the rescaling.** With H := diag(H_p·4^{−p}) and Q := diag(4^q), G = H·Π̃_{[K]}·Q = (H𝐋H^{−1})·(H𝐔Q) =: L·U′.
   - The entries of L are L_{pq} = (H_p/H_q)·4^{q−p}·𝐋_{pq}.
   - v₂(H_p/H_q) = Σ_{q<n≤p} v₂(κ_n) ≥ 5(p − q), because v₂(16n³(n+1)) ≥ 5.
   - Hence v₂(L_{pq}) ≥ 5(p−q) − 2(p−q) + 2 = 3(p−q) + 2.
   - U′ is upper triangular with diagonal entries H_p·𝐔_{pp}.
4. **Symmetrise.** G is symmetric and its leading minors Π_{p<k} H_p𝐔_{pp} are nonzero. The factorisation "(unit lower)·(diagonal)·(unit upper)" of such a matrix is unique, so U′ = ΛLᵀ with Λ_p = H_p𝐔_{pp}. ∎

(No square roots are needed anywhere: the balancing is done by the rational diagonal matrices S, H and Q.)

---

## 7. The glue: Pascal's triangle mod 2

Write each Krawtchouk polynomial in terms of the moment polynomials: p_m = Σ_{q≤m} T_{mq}·P_q. Both families are monic with integer coefficients, so **T is lower unitriangular with 2-adic-integer entries**. Mod 2, p_m ≡ (v+1)^m (§5.3) and P_q ≡ v^q (§6.1), so

  p_m ≡ (v + 1)^m = Σ_q C(m, q)·v^q ≡ Σ_q C(m, q)·P_q (mod 2),

and therefore **T ≡ Pascal's triangle (mod 2)**:

```
1
1 1
1 0 1
1 1 1 1
1 0 0 0 1
1 1 0 0 1 1
1 0 1 0 1 0 1
1 1 1 1 1 1 1 1          (C(m, q) mod 2: Sierpinski's triangle)
```

**Lemma 7.1 (corner lemma).** For all a ≥ 0 and M ≥ 1, det[C(a+i, j)]_{0≤i,j<M} = 1.

*Proof.*
- For i = M−1, M−2, …, 1 in turn, subtract row i−1 from row i. By Pascal's rule C(a+i, j) − C(a+i−1, j) = C(a+i−1, j−1), so row i becomes (0, C(a+i−1, 0), C(a+i−1, 1), …).
- Column 0 is now (1, 0, …, 0). Expanding along it leaves det[C(a+i′, j′)]_{0≤i′,j′<M−1}, the same determinant one size smaller.
- By induction the value is 1. ∎

*Example (a = 5, M = 3).* det [[1, 5, 10], [1, 6, 15], [1, 7, 21]] = 1.

So every "corner" block of T (rows a, …, a+M−1, columns 0, …, M−1) has **odd** determinant.

Set R := T·L. By Proposition 6.1, R is lower unitriangular with 2-adic-integer entries, and R ≡ T (mod 2). In the Krawtchouk basis we now have

  A_K + C·B_K = R·Λ·Rᵀ,  B_K = diag(θ̂₀, …, θ̂_{K−1}),  v₂(Λ_q) = h_q,  v₂(θ̂_m) = b_m + 1.

---

## 8. The cheapest term wins: the coefficient law

**Proposition 8.1.** Let q(Y) := det(A_K + C·B_K + Y·B_K) = Σ_k q_k Y^k. For every M with 0 ≤ M ≤ K,

  v₂(q_{K−M}) = Σ_{q<M} h_q + Σ_{m<K−M} (b_m + 1).

*Proof.*
- **Rewrite as a product.** A_K + C·B_K + Y·B_K = [R | I]·diag(Λ₀, …, Λ_{K−1}, Y·θ̂₀, …, Y·θ̂_{K−1})·[R | I]ᵀ, where [R | I] is the K×2K matrix made of R and the identity.
- **Apply Cauchy–Binet.** The determinant becomes a sum over ways of choosing K of the 2K columns: a set I of columns of R (moment directions) and a set J of columns of the identity (band directions), with |I| + |J| = K. The term for (I, J) is

   (det R[rows not in J, columns I])² · Π_{q∈I} Λ_q · Π_{m∈J} Y·θ̂_m.

- **Bound every term.** The terms with |J| = K − M make up q_{K−M}. Each has v₂ ≥ Σ_{q∈I} h_q + Σ_{m∈J} (b_m + 1), because R has 2-adic-integer entries, v₂(Λ_q) = h_q and v₂(θ̂_m) = b_m + 1.
- **Find the unique cheapest term.** Both lists are strictly increasing (Lemma 2.1). So the smallest possible value of that bound, over |I| = M and |J| = K − M, is attained only by I = {0, …, M−1} and J = {0, …, K−M−1}.
- **Its coefficient is odd.** For that one term the minor is det R[{K−M, …, K−1}, {0, …, M−1}]. Mod 2 this is a corner of Pascal's triangle, which is 1 by Lemma 7.1. So this term has exactly the bound as its valuation, and every other term is strictly more divisible by 2.
- **Conclude.** A sum is exactly as divisible as its unique least divisible term. ∎

*In plain words.* Picking the M cheapest moment directions and the K − M cheapest band directions gives one term. The corner lemma guarantees it is not accidentally even, so nothing can cancel it.

---

## 9. Back to X: the content and the tie rule

**The tropical total.** Put Φ_K(M) := Σ_{m<K−M} b_m + Σ_{q<M} h_q. Since both lists are increasing, the K smallest numbers of the combined list are a beginning segment of each list. Hence

  sum of L_K = min over M of Φ_K(M) =: m₀.

Moreover Φ_K(M+1) − Φ_K(M) = h_M − b_{K−1−M}, which increases in M by at least 4 + 2 = 6 at each step. So Φ_K is strictly convex, and it attains its minimum either at one M or at two consecutive values M, M + 1. The second case happens exactly when h_M = b_{K−1−M}, i.e. when τ_K = 1.

**The shift.** D_K(X) = q(X − C). Put t_k := v₂(q_k) − k. By Proposition 8.1, t_{K−M} = Φ_K(M) for every M. Write C = ε/2 with ε odd (v₂(C) = −1). Expanding q(X − C) = Σ_k q_k(X − C)^k, the coefficient of X^j is

  d_j = 2^j·(−ε)^{−j}·Σ_{k≥j} Q_k·(−ε)^k·C(k, j),  where Q_k := q_k/2^k and v₂(Q_k) = t_k.

**Proof of Theorem 1.**
- **Every coefficient.** Every term in d_j has v₂ ≥ j + m₀, so v₂(d_j) ≥ j + m₀.
- **No tie (τ_K = 0).** Exactly one k has t_k = m₀, and all others have t_k ≥ m₀ + 1. So v₂(d₀) = m₀. Hence min_j v₂(d_j) = m₀ and e₂(K) = −m₀.
- **Tie (τ_K = 1).** The minimum is attained at k₁ and k₂ = k₁ + 1, and all other t_k ≥ m₀ + 6.
  - In d₀ the two main terms are each 2^{m₀}·(odd), so their sum is 2^{m₀}·(even) and v₂(d₀) ≥ m₀ + 1.
  - In d₁ the main terms are multiplied by k₁ and k₂. Exactly one of these is odd, so the pair contributes 2^{m₀}·(odd) and v₂(d₁) = 1 + m₀.
  - Every d_j with j ≥ 2 has v₂ ≥ m₀ + 2.
  - So min_j v₂(d_j) = m₀ + 1, and e₂(K) = −m₀ − 1. ∎

*Remark.* The determinant det A_K = d₀ itself gains **at least** one factor of 2 at a tie. It can gain more, which depends on the next binary digit of the two tied units. That is why det A_K and e₂(K) can differ by more than 1 in the companion "twisted" machine described in the repository.

---

## 10. The fingerprint of A_K (Theorem 2)

A_K = (A_K + C·B_K) − C·B_K = R·Λ·Rᵀ + diag(δ₀, …, δ_{K−1}), with δ_m := −C·θ̂_m and v₂(δ_m) = b_m. Write A_K = [R | I]·diag(Λ, δ)·[R | I]ᵀ and apply Cauchy–Binet to any k×k minor with rows I′ and columns J′. The terms are indexed by a moment set I and a band set J ⊆ I′ ∩ J′ with |I| + |J| = k:

  det A_K[I′, J′] = Σ ± det R[I′∖J, I]·det R[J′∖J, I]·Π_{q∈I} Λ_q·Π_{m∈J} δ_m.

**Lower bound.** Every term has v₂ ≥ Σ_{q∈I} h_q + Σ_{m∈J} b_m, which is at least the sum of the k smallest numbers of the combined list. So the smallest v₂ over all k×k minors is at least that sum.

**It is attained.**
- If the k-th and (k+1)-th smallest numbers differ, take the leading minor I′ = J′ = {0, …, k−1}. Its unique cheapest term has coefficient det R[{k−M, …, k−1}, {0, …, M−1}]², which is odd by the corner lemma.
- If they are equal, h_M = b_{k−1−M}, and two terms tie. Take I′ = J′ = {0, …, k−M−2} ∪ {k−M, …, k}. This excludes the band direction k − M − 1, and the surviving cheapest term has the odd coefficient det R[{k−M, …, k}, {0, …, M}]². This needs k < K, since row k must exist.

So d₁ + ⋯ + d_k equals the sum of the k smallest numbers for every k < K, and the fingerprint agrees with L_K except possibly in its last entry. For k = K there is only one minor (the determinant), whose valuation is m₀ or at least m₀ + 1 at a tie (§9). ∎

---

## 11. Growth: where K²/6 and (K/3)·log₂K come from (Theorem 3)

Let S(n) := Σ_{i<n} s₂(i) be the total number of binary 1s below n. Summing the definitions,

  Φ_K(M) = −2(K−M)(M+1) + 4M(M−1) + S(K) − 4S(M) − S(M+1) − S(K−M).

**Delange's theorem.** S(n) = (n/2)·log₂n + n·F(log₂n) with F bounded (Delange 1975; the elementary estimate S(n) = (n/2)·log₂n + O(n) suffices here). Write M = xK. Then, uniformly for 0 ≤ x ≤ 1,

  Φ_K(xK) = K²(6x² − 2x) − 2x·K·log₂K + O(K),

because the digit-sum terms contribute (K/2)·log₂K·(1 − 5x − (1 − x)) + O(K).

**The minimum.**
- Over real x, the minimum of the main part is at x = 1/6 + O(log K/K). Its value is −(K + log₂K)²/6 = −K²/6 − (K/3)·log₂K + O(log²K), which gives the lower bound.
- Taking M = round(K/6) gives the matching upper bound.
- Hence min_M Φ_K(M) = −K²/6 − (K/3)·log₂K + O(K), and Theorem 3 follows from Theorem 1.

| K | e₂(K) | (e₂ − K²/6 − (K/3)log₂K)/K |
|---:|---:|---:|
| 32 | 258 | 1.06 |
| 64 | 874 | 0.99 |
| 128 | 3146 | 0.91 |
| 192 | 6849 | 1.14 |
| 256 | 11835 | 0.90 |
| 384 | 26113 | 1.14 |

The ratio in the last column wobbles but stays bounded. That is Delange's log-periodic fluctuation of binary digit sums.

---

## 12. What this means for the irrationality of G

**The contest for this family.** Here is the ledger, per K², in natural-log units. The analytic and odd-prime entries are the repository's own asymptotic analysis (see `HANDOFF_catalan_hankel-1.md`), not part of the proof above.

| | per K² |
|---|---:|
| **gain:** how fast the approximations converge | ln 9 ≈ 2.197 |
| **cost at odd primes** | 8/3 ≈ 2.667 |
| **cost at the prime 2** (Theorems 1 and 3 of this paper) | (1/6)·ln 2 ≈ 0.116 |
| **net** (cost − gain) | ≈ +0.585 |

A positive net means the denominators win. **This family of determinants cannot prove G irrational.** It is not the prime 2's fault: even the odd primes alone (2.667) outweigh the gain (2.197). What this paper contributes is an exact, proved account of the prime-2 share. That is the kind of accounting any claimed proof has to get right. A recent claimed proof of the irrationality of G (Sun, arXiv 2609.04176) was found incomplete, and the published critique (Wachs, arXiv 2609.22339) locates most of the gap in an unaccounted prime-2 contribution.

**The 2-adic side is already known.** At X = C, Proposition 8.1 with M = K shows that D_K(C) is divisible by exactly 2^{h₀ + h₁ + ⋯ + h_{K−1}}, about 2^{4K²}. For K = 8 that is exactly 2¹⁷⁵. That is enormous 2-adic smallness. The 2-adic constant C itself is known to be irrational (Calegari 2005; Beukers 2008). So the same machine lives in two worlds:
- at the real point X = G its value is small in the usual sense;
- at the 2-adic point X = C its value is small 2-adically.

The real-world question is the open one.

**What a proof of the irrationality of G would need.** A construction whose arithmetic cost, summed over all primes, is less than its analytic gain. The closest known success is the 2024 proof that the sister constant L(2, χ₋₃) is irrational, by Calegari, Dimitrov and Tang, using "arithmetic holonomy bounds". G = L(2, χ₋₄) remains open. The tools in this paper give exact prime-by-prime bookkeeping, and they extend to other primes (Pascal's triangle mod p is again unitriangular with corner determinants 1, by Lucas's theorem). They are the instruments with which such a construction would have to be measured.

---

## Appendix A. Where the machine's numbers come from

**A.1 Moments.** Since d/dt(−1/sinh(πt)) = π·cosh(πt)/sinh²(πt), integration by parts gives

  μ_k = ∫₀^∞ (4t²)^k·t²·π cosh(πt)/sinh²(πt) dt = 4^k·(2k+2)·∫₀^∞ t^{2k+1}/sinh(πt) dt.

The boundary terms vanish. Expanding 1/sinh u = 2Σ_{n≥0} e^{−(2n+1)u} gives

  ∫₀^∞ u^{s−1}/sinh(u) du = 2(1 − 2^{−s})·Γ(s)·ζ(s).

With ζ(2n) = |B_{2n}|(2π)^{2n}/(2·(2n)!) this yields μ_k = 4^k·(2^{2k+2} − 1)·|B_{2k+2}|.

**A.2 Pole values.** Let I(s) := ∫₀^∞ ω(t)/(4t² + s²) dt for s > 0.

1. **Integrate by parts** as in A.1: I(s) = 2s²·∫₀^∞ t/((4t²+s²)²·sinh(πt)) dt = −s·J′(s), where J(s) := ∫₀^∞ t/((4t²+s²)·sinh(πt)) dt.
2. **Expand the sinh.** The Mittag-Leffler expansion πt/sinh(πt) = 1 + 2Σ_{n≥1} (−1)ⁿ·t²/(t² + n²) has partial sums bounded by 1, so it may be integrated term by term. The standard integrals are ∫₀^∞ dt/(4t²+s²) = π/(4s) and ∫₀^∞ t²/((t²+n²)(4t²+s²)) dt = π/(4(s+2n)). Together they give

   J(s) = ½·Σ_{n≥0} (−1)ⁿ/(s + 2n) − 1/(4s),  so  I(s) = (s/2)·Σ_{n≥0} (−1)ⁿ/(s + 2n)² − 1/(4s).

3. **Evaluate at s = 2j + 1.** The sum is Σ_{m≥j} (−1)^{m−j}/(2m+1)² = (−1)^j·(G − β_j), so I(2j+1) = F_j(G).

(Both formulas were also checked numerically to 40 digits, for k < 6 and j < 6.)

**A.3 Why M_K(G) is a Gram matrix.** At X = G every pole value is the true integral, so the (i, k) entry is ∫₀^∞ (4t²)^{i+k}·ω(t)/N_K(4t²) dt. That is the Gram matrix of the positive weight ω(t)/N_K(4t²) on (0, ∞), so it is positive definite and D_K(G) > 0.

---

## Appendix B. The classical orthogonal polynomials

**The weight.**
- Using |Γ(1+it)|² = πt/sinh(πt), |Γ(it)|² = π/(t·sinh(πt)) and |Γ(2it)|² = π/(2t·sinh(2πt)),

  |Γ(1+it)² Γ(it)|²/|Γ(2it)|² = 4π²t²·cosh(πt)/sinh²(πt) = 4π·ω(t).

- So ω is the continuous dual Hahn weight with parameters (a, b, c) = (1, 1, 0), divided by 4π.
- The case c = 0 is the limit c → 0⁺. There the moments converge by dominated convergence, and the orthogonality relations, which are polynomial identities in the moments and the recurrence coefficients, pass to the limit.

**The recurrence.**
- In the variable u = t², the standard three-term relation (Koekoek–Lesky–Swarttouw, §9.3) for the monic polynomials p_n = (−1)ⁿS_n reads

  u·p_n = p_{n+1} + (A_n + C_n − a²)·p_n + A_{n−1}·C_n·p_{n−1},

  with A_n = (n+a+b)(n+a+c) and C_n = n(n+b+c−1).
- For (1, 1, 0) this gives (2n+1)(n+1) and n³(n+1). The rescaling v = 4u multiplies these by 4 and by 16.
- The norms follow from H_n = μ₀·κ₁⋯κ_n with μ₀ = ½.
- Independent check: the recurrence and the norms agree exactly with the moments μ_k for all n < 30 (proof_checks/th_paper_checks.py).

---

## Appendix C. Check it yourself

Everything is in the repository and runs on Python 3.12 with `python-flint`, `sympy`, `mpmath` and `numpy`.

```
python twoadic_law.py claw                      # the formula vs the 134 computed values of e2(K)
python proof_checks/th_catalan_struct.py 12     # the whole proof chain for K = 12
python proof_checks/th_thm7_rigorous.py catalan 8
python proof_checks/th_paper_checks.py          # Appendices A-B, the constant C, Theorem 3's table
```

The end-to-end script `th_catalan_struct.py` checks, for a given K:
- that A_K + C·B_K is the moment matrix of §6.3;
- the fingerprints of both pieces;
- that T ≡ Pascal (mod 2);
- Proposition 8.1 for **every** M.

It passes for K = 6, 8, 12, 16, 20, 24. The companion file `twoadic_law_proof.md` is the full technical version, including the twisted machine.

---

## Appendix D. Glossary

- **2-adic valuation v₂(x):** the number of factors of 2 in x (numerator minus denominator).
- **2-adic numbers ℚ₂:** numbers built so that divisibility by high powers of 2 means "small".
- **Content exponent e₂(K):** the largest power of 2 in any denominator of the coefficients of D_K.
- **Fingerprint (Smith form):** the list of powers of 2 that a matrix "really" has, invariant under odd-determinant integer changes of basis.
- **Hankel matrix:** a matrix whose (i, k) entry depends only on i + k, such as a table of moments.
- **Orthogonal polynomials:** polynomials P₀, P₁, … with φ(P_m·P_n) = 0 for m ≠ n. In their basis a moment matrix becomes diagonal.
- **Krawtchouk polynomials:** the orthogonal polynomials of the binomial (coin-flip) distribution.
- **Continuous dual Hahn polynomials:** a classical family from the Askey scheme of hypergeometric orthogonal polynomials.
- **Cauchy–Binet:** the determinant of a product, expanded over choices of columns.
- **Tropical minimum:** the minimum over M of a sum of two list-prefixes. In valuation arguments it plays the role that a determinant plays in ordinary algebra.

---

## References

- R. Apéry, *Irrationalité de ζ(2) et ζ(3)*, Astérisque 61 (1979), 11–13.
- F. Beukers, *Irrationality of some p-adic L-values*, Acta Math. Sinica (Engl. Ser.) 24 (2008), 663–686; arXiv:math/0603277.
- J. M. Borwein, N. J. Calkin, D. Manna, *Euler–Boole summation revisited*, Amer. Math. Monthly 116 (2009), 387–412.
- F. Calegari, *Irrationality of certain p-adic periods for small p*, Int. Math. Res. Not. 2005, no. 20, 1235–1249; arXiv:math/0408214.
- F. Calegari, V. Dimitrov, Y. Tang, *The linear independence of 1, ζ(2), and L(2, χ₋₃)*, arXiv:2408.15403 (2024).
- H. Delange, *Sur la fonction sommatoire de la fonction « somme des chiffres »*, L'Enseignement Math. 21 (1975), 31–47.
- P. W. Kasteleyn, *The statistics of dimers on a lattice*, Physica 27 (1961), 1209–1225; H. N. V. Temperley, M. E. Fisher, *Dimer problem in statistical mechanics — an exact result*, Phil. Mag. 6 (1961), 1061–1063.
- R. Koekoek, P. A. Lesky, R. F. Swarttouw, *Hypergeometric Orthogonal Polynomials and Their q-Analogues*, Springer, 2010 (§9.3 continuous dual Hahn, §9.11 Krawtchouk).
- E. E. Kummer, *Über die Ergänzungssätze zu den allgemeinen Reciprocitätsgesetzen*, J. reine angew. Math. 44 (1852), 93–146.
- Z.-W. Sun, *Catalan's constant is irrational*, arXiv:2609.04176 (2026), and D. Wachs, *A note on a recent claimed proof of the irrationality of Catalan's constant*, arXiv:2609.22339 (2026).
- J. A. Wilson, *Some hypergeometric orthogonal polynomials*, SIAM J. Math. Anal. 11 (1980), 690–701.
